#!/usr/bin/env python3
"""Stage 2 -- classify junctions, fingerprint the delivery, build the locus table.

Takes the identity map from `construct_identity.py` plus one or more junction
tables and produces:

  * a per-junction classification: host-vector / vector-vector /
    vector-unmapped, and unique / shared / vector-end on the vector side
  * artefact flags, so cross-mapping and library noise are separated from real
    integration junctions BEFORE anything is counted
  * a **delivery fingerprint**: whether host-vector junctions cluster at a
    fixed pair of vector coordinates (integrase- or transposase-mediated, one
    copy per site) or are scattered (end-joining, concatemer-capable). This
    decides which copy-number arithmetic is valid downstream, and getting it
    wrong is the single largest error available in this analysis.
  * a locus table with flanks paired where both were recovered, and the
    target-site footprint of each pair checked against the modality
  * cross-sample support ratios with confidence intervals, which act as a
    vector-of-origin discriminator when constructs share sequence

Usage
-----
    python classify_junctions.py --identity identity.json --mapping pAlpha \\
        --modality integrase \\
        --table parental=junctions_p.tsv --table adapted=junctions_a.tsv \\
        --norm adapted=0.87 \\
        --out-junctions junctions.tsv --out-loci loci.tsv --md stage2.md

`--vector-ends` lists the defined ends of the delivered molecule in that
construct's own numbering: the restriction site for linearized plasmid, the
LTR termini for a provirus, the ITR termini for a transposon. Its meaning
depends on `--modality` -- for end-joining delivery a junction there is
ambiguous and is flagged, while for integrase or transposase delivery it is
the canonical, expected junction and its ABSENCE is what deserves suspicion.
"""
from __future__ import annotations

import argparse
import math
import os
import sys
from collections import defaultdict
from typing import Dict, List, Optional, Sequence, Tuple

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from via_common import (Junction, NameMatcher, fmt_warnings,  # noqa: E402
                        load_json, natural_key, parse_clip_map,
                        parse_column_spec, read_junction_table)


# --------------------------------------------------------------------------
# modalities
# --------------------------------------------------------------------------
#
# `footprint` is the expected value of
#       delta = host_pos(right flank) - host_pos(left flank) - 1
# A target-site duplication of d bp makes the two flanks overlap in reference
# coordinates, so delta = -d. A target-site deletion of m bp gives delta = +m.
# A clean join gives 0.

MODALITIES: Dict[str, Dict] = {
    "nhej": {
        "label": "end-joining (linearized plasmid, naked DNA, minicircle)",
        "fixed_ends": False,
        "concatemers": True,
        "footprint": None,            # variable: small deletions are the norm
        "footprint_range": (-3, 40),
        "copies_per_site": "variable (1 to tens; tandem arrays are the rule)",
    },
    "integrase": {
        "label": "integrase-mediated (lentivirus, gammaretrovirus)",
        "fixed_ends": True,
        "concatemers": False,
        "footprint": -5,              # HIV-1; MLV is -4
        "footprint_range": (-6, -3),
        "copies_per_site": "1",
    },
    "transposon-pb": {
        "label": "PiggyBac transposase",
        "fixed_ends": True,
        "concatemers": False,
        "footprint": -4,              # TTAA duplicated
        "footprint_range": (-4, -4),
        "copies_per_site": "1",
    },
    "transposon-sb": {
        "label": "Sleeping Beauty transposase",
        "fixed_ends": True,
        "concatemers": False,
        "footprint": -2,              # TA duplicated
        "footprint_range": (-3, -1),
        "copies_per_site": "1",
    },
    "hdr": {
        "label": "homology-directed targeted integration (knock-in)",
        "fixed_ends": True,
        "concatemers": True,          # concatemers at the target locus happen
        "footprint": 0,
        "footprint_range": (-2, 2),
        "copies_per_site": "1 per targeted allele, but check for concatemers",
    },
    "auto": {
        "label": "undeclared -- inferred from the data",
        "fixed_ends": None,
        "concatemers": None,
        "footprint": None,
        "footprint_range": None,
        "copies_per_site": "unknown until the modality is settled",
    },
}


# --------------------------------------------------------------------------
# classification
# --------------------------------------------------------------------------

def block_at(pos: int, mapdata: Dict) -> Optional[Dict]:
    for b in mapdata["blocks"]:
        if b["start"] <= pos <= b["end"]:
            return b
    return None


def classify(js: Sequence[Junction], ident: Dict, mapping: str,
             ends: Sequence[int], end_tol: int) -> None:
    md = ident["constructs"][mapping]
    for j in js:
        if not j.vector:
            continue
        b = block_at(j.vector_pos, md)
        if any(abs(j.vector_pos - e) <= end_tol for e in ends):
            j.region = "vector-end"
        elif b is None:
            j.region = "off-reference"
        elif b["unique"]:
            j.region = "unique"
        else:
            j.region = "shared"
        j.membership = " + ".join(b["membership"]) if b else ""
        j.feature = b["features"] if b else "off-reference"


def flag_artefacts(js: Sequence[Junction], ident: Dict, mapping: str,
                   modality: str, ends: Sequence[int], end_tol: int,
                   boundary_tol: int, long_filler: int,
                   inverted_span: Tuple[int, int]) -> None:
    md = ident["constructs"][mapping]
    ref_len = md["length"]
    circular = md.get("circular", False)
    mod = MODALITIES[modality]
    boundaries = md["boundaries"]

    def near_boundary(pos: int) -> Optional[Dict]:
        for bd in boundaries:
            if min(abs(pos - bd["pos"]), abs(pos - bd["next_pos"])) <= boundary_tol:
                return bd
        return None

    for j in js:
        if j.kind == "vector-unmapped":
            bd = near_boundary(j.vector_pos)
            if bd is not None and j.filler >= long_filler:
                who = ",".join(bd["clipping_constructs"]) or "?"
                j.artefact = "boundary-crossmap(@%d;%s)" % (bd["pos"], who)
            elif j.filler >= long_filler:
                j.artefact = "long-filler"

        elif j.kind == "vector-vector":
            same_clip = bool(j.clip_a) and j.clip_a == j.clip_b
            span = abs(j.pos_a - j.pos_b)
            at_both_ends = (
                len(ends) >= 2
                and any(abs(j.pos_a - e) <= end_tol for e in ends)
                and any(abs(j.pos_b - e) <= end_tol for e in ends)
                and j.pos_a != j.pos_b)
            at_ref_edges = (not circular
                            and min(j.pos_a, j.pos_b) <= 1 + end_tol
                            and max(j.pos_a, j.pos_b) >= ref_len - end_tol)
            # A self-junction (span ~0, same clip side) is a palindromic
            # head-to-head or tail-to-tail concatemer join and is REAL. The
            # artefact class has a short but non-zero span.
            if same_clip and inverted_span[0] <= span <= inverted_span[1] \
                    and j.filler < 0:
                j.artefact = "short-range-inverted"
            elif at_both_ends or at_ref_edges:
                j.artefact = ("vector-end-religation" if mod["concatemers"] is not False
                              else "vector-end-religation(unexpected-for-modality)")

        elif j.kind == "host-vector":
            if j.region == "vector-end":
                # A host junction at a declared vector end is the canonical
                # integration junction, not an artefact -- and for a provirus
                # or transposon mapped to a linear reference those ends sit at
                # position 1 and the last base, exactly where the
                # reference-edge rule would otherwise fire. Checking the
                # declared ends first is what keeps integrase-mediated
                # delivery from having every real junction thrown away.
                continue
            if not circular and (j.vector_pos <= 1 or j.vector_pos >= ref_len):
                j.artefact = "reference-edge"
            elif mod["fixed_ends"] and ends and j.region != "vector-end":
                # For integrase/transposase delivery the vector breakpoint is
                # not free. A host junction elsewhere in the vector is either a
                # rearranged provirus, backbone carryover, or an artefact -- in
                # every case it is not a canonical integration junction.
                j.artefact = "off-end-host-junction"


# --------------------------------------------------------------------------
# delivery fingerprint
# --------------------------------------------------------------------------

def fingerprint(js: Sequence[Junction], ends: Sequence[int], end_tol: int,
                cluster_tol: int = 10) -> Dict:
    """Do host-vector junctions cluster at a fixed pair of vector coordinates?

    Integrase and transposase cut the vector at defined termini, so every
    integration junction sits at one of two vector coordinates. End-joining
    breaks the molecule wherever it happens to break, so vector breakpoints
    scatter. The difference is visible in the junction table alone, before any
    assumption about what was transfected -- which matters, because the
    delivery recorded in a line's paperwork is not always what is in the cells.
    """
    hv = [j for j in js if j.kind == "host-vector"]
    if not hv:
        return {"verdict": "no host-vector junctions", "n": 0}

    buckets: Dict[int, List[Junction]] = defaultdict(list)
    for j in sorted(hv, key=lambda x: x.vector_pos):
        for key in buckets:
            if abs(j.vector_pos - key) <= cluster_tol:
                buckets[key].append(j)
                break
        else:
            buckets[j.vector_pos].append(j)

    ranked = sorted(buckets.items(), key=lambda kv: -len(kv[1]))
    top2 = ranked[:2]
    in_top2 = sum(len(v) for _, v in top2)
    frac = in_top2 / len(hv)
    reads_top2 = sum(j.support for _, v in top2 for j in v)
    reads_all = sum(j.support for j in hv) or 1

    if len(hv) < 4:
        verdict = "too few host-vector junctions (%d) to fingerprint" % len(hv)
        suggests = None
    elif frac >= 0.8 and len(buckets) <= max(3, len(hv) // 3):
        verdict = ("FIXED vector ends: %.0f%% of host-vector junctions sit at "
                   "%d vector coordinate(s)" % (100 * frac, len(top2)))
        suggests = "integrase / transposase (one copy per site)"
    elif frac <= 0.5:
        verdict = ("SCATTERED vector breakpoints: %d distinct vector positions "
                   "across %d junctions" % (len(buckets), len(hv)))
        suggests = "end-joining (concatemer arrays possible)"
    else:
        verdict = ("MIXED: %.0f%% at the two commonest vector coordinates, "
                   "%d distinct positions overall" % (100 * frac, len(buckets)))
        suggests = ("either delivery plus rearrangement, or two delivery "
                    "routes in one line")

    return {
        "verdict": verdict,
        "suggests": suggests,
        "n": len(hv),
        "distinct_vector_positions": len(buckets),
        "fraction_at_top2": round(frac, 3),
        "read_fraction_at_top2": round(reads_top2 / reads_all, 3),
        "top_positions": [{"vector_pos": k, "junctions": len(v),
                           "reads": sum(x.support for x in v),
                           "at_declared_end": any(abs(k - e) <= end_tol for e in ends)}
                          for k, v in ranked[:5]],
        "declared_ends": list(ends),
    }


# --------------------------------------------------------------------------
# locus assembly
# --------------------------------------------------------------------------

def build_loci(js: Sequence[Junction], pair_window: int,
               clip_usable: bool) -> List[Dict]:
    """Group host junctions into loci and pair flanks.

    A complete insertion gives an `R` junction (reference runs up to it, then
    insert: the LEFT flank) and an `L` junction (insert, then reference
    resumes: the RIGHT flank). Only an R+L pair within `pair_window` is a
    matched pair; two junctions of the same clip side are two left flanks or
    two right flanks and must not be paired. Getting this wrong inflates the
    site count.

    Without a usable clip side nothing can be paired, so every junction
    becomes its own single-flank locus. That is a real loss of a consistency
    check, not a silent degradation -- the caller reports it.
    """
    gen = [j for j in js if j.kind == "host-vector" and not j.artefact]
    by_contig: Dict[str, List[Junction]] = defaultdict(list)
    for j in gen:
        by_contig[j.host_contig].append(j)

    loci: List[Dict] = []
    for contig in sorted(by_contig, key=natural_key):
        items = sorted(by_contig[contig], key=lambda j: j.host_pos)
        used: set = set()
        for i, j in enumerate(items):
            if i in used:
                continue
            partner_i = None
            if clip_usable:
                for n in range(i + 1, len(items)):
                    if n in used:
                        continue
                    if items[n].host_pos - j.host_pos > pair_window:
                        break
                    if {j.host_clip, items[n].host_clip} == {"R", "L"}:
                        partner_i = n
                        break
            if partner_i is not None:
                other = items[partner_i]
                used.update({i, partner_i})
                left = j if j.host_clip == "R" else other
                right = other if left is j else j
                loci.append({
                    "contig": contig,
                    "start": left.host_pos, "end": right.host_pos,
                    "flanks": 2, "rep": left, "left": left, "right": right,
                    "footprint": right.host_pos - left.host_pos - 1,
                    "support": left.support + right.support,
                })
            else:
                used.add(i)
                loci.append({
                    "contig": contig, "start": j.host_pos, "end": j.host_pos,
                    "flanks": 1, "rep": j,
                    "left": j if j.host_clip == "R" else None,
                    "right": j if j.host_clip == "L" else None,
                    "footprint": None, "support": j.support,
                })
    return loci


def check_footprints(loci: Sequence[Dict], modality: str) -> Dict:
    """Compare observed target-site footprints against the modality's expectation.

    A duplicated target site is a positive, independent confirmation that a
    paired flank really is one insertion: it is the enzyme's signature, and a
    library chimera has no reason to carry it.
    """
    mod = MODALITIES[modality]
    obs = [L["footprint"] for L in loci if L["footprint"] is not None]
    out: Dict = {"n_paired": len(obs), "observed": sorted(obs),
                 "expected": mod["footprint"],
                 "expected_range": mod["footprint_range"]}
    if not obs or mod["footprint_range"] is None:
        out["verdict"] = ("no paired flanks to check" if not obs else
                          "modality undeclared: no footprint expectation to test")
        return out
    lo, hi = mod["footprint_range"]
    inside = [v for v in obs if lo <= v <= hi]
    out["n_consistent"] = len(inside)
    out["fraction_consistent"] = round(len(inside) / len(obs), 3)
    if len(inside) == len(obs):
        out["verdict"] = ("all %d paired flanks carry the expected footprint "
                          "(%d..%d) -- strong independent support" % (len(obs), lo, hi))
    elif inside:
        out["verdict"] = ("%d of %d paired flanks carry the expected footprint "
                          "(%d..%d); the others need explaining individually"
                          % (len(inside), len(obs), lo, hi))
    else:
        out["verdict"] = ("NO paired flank carries the expected footprint "
                          "(%d..%d). Either the modality is wrong, the clip "
                          "convention is inverted, or these pairs are not pairs."
                          % (lo, hi))
    return out


# --------------------------------------------------------------------------
# cross-sample ratios
# --------------------------------------------------------------------------

def ratio_ci(a: int, b: int, norm_factor: float = 1.0
             ) -> Tuple[float, float, float]:
    """Support ratio b/a with a 95% interval from the Poisson approximation."""
    if a <= 0 or b <= 0:
        return (float("nan"),) * 3
    r = (b * norm_factor) / a
    se = math.sqrt(1.0 / a + 1.0 / b)
    return r, r * math.exp(-1.96 * se), r * math.exp(1.96 * se)


def cross_sample(js: Sequence[Junction], ref_sample: str,
                 norms: Dict[str, float]) -> Dict[Tuple, Dict]:
    idx: Dict[Tuple, Dict[str, int]] = defaultdict(dict)
    for j in js:
        if j.kind != "host-vector":
            continue
        idx[(j.host_contig, j.host_pos, j.vector_pos)][j.sample] = j.support
    out: Dict[Tuple, Dict] = {}
    for k, d in idx.items():
        a = d.get(ref_sample, 0)
        per = {}
        for s, b in d.items():
            if s == ref_sample:
                continue
            r, lo, hi = ratio_ci(a, b, norms.get(s, 1.0))
            per[s] = {"support": b, "ratio": r, "lo": lo, "hi": hi}
        out[k] = {"ref_support": a, "others": per}
    return out


# --------------------------------------------------------------------------
# reporting
# --------------------------------------------------------------------------

def render_md(mapping: str, modality: str, samples: Sequence[str],
              infos: Dict[str, Dict], fp: Dict, foot: Dict,
              loci: Sequence[Dict], counts: Dict, warnings: Sequence[str]) -> str:
    mod = MODALITIES[modality]
    L = ["# Stage 2 -- junction classification (%s)" % mapping, "",
         "**Declared modality:** %s" % mod["label"],
         "", "**Expected copies per site:** %s" % mod["copies_per_site"], ""]

    L += ["## How the input was read", "",
          "| Sample | rows | junctions | clip side | columns |",
          "|---|---|---|---|---|"]
    for s in samples:
        i = infos[s]
        clip = i["clip"]
        L.append("| %s | %d | %d | %s | %s |"
                 % (s, i["rows_read"], i["junctions"],
                    "usable (%s)" % clip.reason if clip.ok
                    else "**UNUSABLE** -- %s" % clip.reason,
                    ", ".join("%s=%s" % (k, v) for k, v in i["columns"].items())))
    L.append("")
    if warnings:
        L += ["### Warnings", ""] + ["- %s" % w for w in warnings] + [""]

    L += ["## Delivery fingerprint", "",
          "> %s" % fp["verdict"], ""]
    if fp.get("suggests"):
        L += ["Consistent with: **%s**." % fp["suggests"], ""]
        if mod["fixed_ends"] is not None:
            agrees = (mod["fixed_ends"] and "FIXED" in fp["verdict"]) or \
                     (not mod["fixed_ends"] and "SCATTERED" in fp["verdict"])
            L += ["Declared modality %s the fingerprint.%s"
                  % ("**agrees with**" if agrees else "**DISAGREES with**",
                     "" if agrees else " Settle this before trusting any copy "
                     "number or site count: the two modalities use different "
                     "arithmetic."), ""]
    if fp.get("top_positions"):
        L += ["| Vector position | junctions | reads | at a declared end? |",
              "|---|---|---|---|"]
        for t in fp["top_positions"]:
            L.append("| %d | %d | %d | %s |"
                     % (t["vector_pos"], t["junctions"], t["reads"],
                        "yes" if t["at_declared_end"] else "no"))
        L.append("")

    L += ["## Junction census (%s)" % samples[0], "",
          "| Class | junctions | flagged artefact |", "|---|---|---|"]
    for kind in ("host-vector", "vector-vector", "vector-unmapped", "host-host"):
        L.append("| %s | %d | %d |" % (kind, counts["kind"].get(kind, 0),
                                       counts["artefact_by_kind"].get(kind, 0)))
    L.append("")
    L += ["Host-vector junctions surviving artefact flags, by vector region:", "",
          "| Region | junctions |", "|---|---|"]
    for reg in ("unique", "shared", "vector-end", "off-reference"):
        L.append("| %s | %d |" % (reg, counts["region"].get(reg, 0)))
    L.append("")
    if counts["artefact_kinds"]:
        L += ["Artefact flags raised:", "", "| Flag | junctions | reads |",
              "|---|---|---|"]
        for k, (n, reads) in sorted(counts["artefact_kinds"].items()):
            L.append("| `%s` | %d | %d |" % (k, n, reads))
        L += ["", "Flags are starting points. Confirm each by the test in "
              "`references/artefacts.md` before accepting it.", ""]

    L += ["## Target-site footprint", "", "> %s" % foot["verdict"], ""]
    if foot.get("observed"):
        L += ["Observed footprints: %s (negative = target-site duplication, "
              "positive = bases deleted at the target)."
              % ", ".join(str(v) for v in foot["observed"]), ""]

    paired = sum(1 for x in loci if x["flanks"] == 2)
    L += ["## Loci", "",
          "%d locus/loci: %d with both flanks, %d single-flank."
          % (len(loci), paired, len(loci) - paired), "",
          "A one-sided site list is normal -- a flank landing in repetitive "
          "sequence cannot be placed uniquely. It is not evidence that the "
          "site is false.", ""]
    return "\n".join(L)


# --------------------------------------------------------------------------

def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--identity", required=True)
    ap.add_argument("--mapping", required=True,
                    help="construct name the reads were mapped against")
    ap.add_argument("--modality", default="auto", choices=sorted(MODALITIES),
                    help="delivery route. Decides which arithmetic is valid "
                         "downstream and what footprint to expect. 'auto' "
                         "infers it from the data and refuses to assume.")
    ap.add_argument("--table", action="append", required=True,
                    dest="table", metavar="SAMPLE=PATH",
                    help="junction table; repeatable, any number of samples")
    ap.add_argument("--ref-sample", default="",
                    help="sample that cross-sample ratios are taken against "
                         "(default: the first --table)")
    ap.add_argument("--norm", action="append", default=[], metavar="SAMPLE=FACTOR",
                    help="scale factor for a sample, the ratio of deduplicated "
                         "read counts to the reference sample; repeatable")
    ap.add_argument("--columns", default="",
                    help="explicit column mapping when auto-detection fails, "
                         "e.g. 'ref1=chrom,pos1=3,support=n_reads' "
                         "(header names or 1-based indices)")
    ap.add_argument("--clip-map", default="",
                    help="declare what the clip column's values mean, e.g. "
                         "'rightmost=R,leftmost=L'. R = reference continues up "
                         "to this base then novel sequence begins; L = novel "
                         "sequence then reference resumes. Confirm against a "
                         "junction of known orientation before trusting it.")
    ap.add_argument("--construct-alias", action="append", default=[],
                    metavar="CONSTRUCT=TABLENAME",
                    help="map a construct name in the identity map onto the "
                         "name used in the junction table; repeatable")
    ap.add_argument("--unmapped-token", default="",
                    help="extra token(s) your caller uses for an unmapped "
                         "partner, comma separated")
    ap.add_argument("--vector-ends", "--ends", dest="ends", default="",
                    help="comma-separated coordinates of the delivered "
                         "molecule's defined ends, in this construct's numbering")
    ap.add_argument("--end-tol", type=int, default=3)
    ap.add_argument("--min-support", type=int, default=1,
                    help="drop junctions below this support. Default 1: the "
                         "lowest-support junctions are often the most "
                         "informative, so filtering is left to you, explicitly.")
    ap.add_argument("--pair-window", type=int, default=500,
                    help="maximum host distance between two flanks of one "
                         "insertion (default 500). Raise it only if you expect "
                         "a large target-site deletion; a wide window pairs "
                         "unrelated junctions and deflates the site count.")
    ap.add_argument("--boundary-tol", type=int, default=3)
    ap.add_argument("--long-filler", type=int, default=0,
                    help="filler at or above this is a cross-mapping "
                         "signature. Default 0 = derive it from --read-length "
                         "as read_length - 2*25.")
    ap.add_argument("--read-length", type=int, default=150)
    ap.add_argument("--inverted-min-span", type=int, default=15)
    ap.add_argument("--inverted-max-span", type=int, default=200)
    ap.add_argument("--no-prefix-match", action="store_true",
                    help="require construct names to match exactly or via an "
                         "explicit --construct-alias")
    ap.add_argument("--out-junctions")
    ap.add_argument("--out-loci")
    ap.add_argument("--md")
    a = ap.parse_args(argv)

    try:
        ident = load_json(a.identity)
    except (OSError, ValueError) as exc:
        print("error: cannot read identity map: %s" % exc, file=sys.stderr)
        return 2
    if a.mapping not in ident["constructs"]:
        print("error: %r is not in the identity map (%s)"
              % (a.mapping, ", ".join(ident["constructs"])), file=sys.stderr)
        return 2

    ends = []
    for x in a.ends.split(","):
        x = x.strip()
        if not x:
            continue
        if not x.lstrip("-").isdigit():
            print("error: --vector-ends takes integers, got %r" % x, file=sys.stderr)
            return 2
        ends.append(int(x))

    long_filler = a.long_filler or max(40, a.read_length - 50)

    aliases: Dict[str, str] = {}
    for item in a.construct_alias:
        k, _, v = item.partition("=")
        if not v:
            print("error: --construct-alias needs CONSTRUCT=TABLENAME, got %r"
                  % item, file=sys.stderr)
            return 2
        aliases[k.strip()] = v.strip()

    norms: Dict[str, float] = {}
    for item in a.norm:
        k, _, v = item.partition("=")
        try:
            norms[k.strip()] = float(v)
        except ValueError:
            print("error: --norm needs SAMPLE=FACTOR, got %r" % item, file=sys.stderr)
            return 2

    try:
        colspec = parse_column_spec(a.columns) if a.columns else {}
        clipmap = parse_clip_map(a.clip_map) if a.clip_map else None
    except ValueError as exc:
        print("error: %s" % exc, file=sys.stderr)
        return 2

    matcher = NameMatcher(list(ident["constructs"]), aliases,
                          allow_prefix=not a.no_prefix_match)

    all_js: List[Junction] = []
    samples: List[str] = []
    infos: Dict[str, Dict] = {}
    warnings: List[str] = []
    clip_usable = True

    for spec in a.table:
        sample, _, path = spec.partition("=")
        if not path:
            print("error: --table needs SAMPLE=PATH, got %r" % spec, file=sys.stderr)
            return 2
        sample = sample.strip()
        if sample in infos:
            print("error: duplicate sample name %r" % sample, file=sys.stderr)
            return 2
        try:
            js, info = read_junction_table(
                path, sample, a.mapping, matcher,
                min_support=a.min_support, columns=colspec,
                unmapped_tokens=[t for t in a.unmapped_token.split(",") if t],
                clip_map=clipmap)
        except (OSError, ValueError) as exc:
            print("error: %s" % exc, file=sys.stderr)
            return 2

        if not info["vector_side_found"]:
            print("error: no reference name in %s matched a construct in the "
                  "identity map (%s). Every junction would be misclassified as "
                  "host-to-host, which is silent and catastrophic, so this is "
                  "fatal. Fix it with --construct-alias CONSTRUCT=TABLENAME, "
                  "using the name as it appears in the table."
                  % (path, ", ".join(ident["constructs"])), file=sys.stderr)
            return 2

        if not info["clip"].ok:
            clip_usable = False
        classify(js, ident, a.mapping, ends, a.end_tol)
        flag_artefacts(js, ident, a.mapping, a.modality, ends, a.end_tol,
                       a.boundary_tol, long_filler,
                       (a.inverted_min_span, a.inverted_max_span))
        all_js += js
        samples.append(sample)
        infos[sample] = info
        warnings += ["[%s] %s" % (sample, w) for w in info["warnings"]]

        print("  %s: %d junctions from %d rows (%s%s)"
              % (sample, info["junctions"], info["rows_read"],
                 "header" if info["header"] else "positional columns",
                 ", %d unparsable" % info["rows_unparsable"]
                 if info["rows_unparsable"] else ""), file=sys.stderr)

    seen = {w for w in warnings}
    warnings = sorted(seen)
    if warnings:
        print(fmt_warnings(warnings), file=sys.stderr)

    ref_sample = a.ref_sample or samples[0]
    if ref_sample not in infos:
        print("error: --ref-sample %r is not one of %s"
              % (ref_sample, ", ".join(samples)), file=sys.stderr)
        return 2

    primary = [j for j in all_js if j.sample == ref_sample]
    fp = fingerprint(primary, ends, a.end_tol)
    loci = build_loci(primary, a.pair_window, clip_usable)
    foot = check_footprints(loci, a.modality)
    xs = cross_sample(all_js, ref_sample, norms) if len(samples) > 1 else {}

    counts = {"kind": defaultdict(int), "artefact_by_kind": defaultdict(int),
              "region": defaultdict(int), "artefact_kinds": {}}
    for j in primary:
        counts["kind"][j.kind] += 1
        if j.artefact:
            counts["artefact_by_kind"][j.kind] += 1
            key = j.artefact.split("(")[0]
            n, reads = counts["artefact_kinds"].get(key, (0, 0))
            counts["artefact_kinds"][key] = (n + 1, reads + j.support)
        elif j.kind == "host-vector":
            counts["region"][j.region] += 1

    if a.out_junctions:
        cols = ["sample", "mapping", "kind", "host_contig", "host_pos",
                "host_clip", "vector", "vector_pos", "vector_clip", "filler",
                "support", "region", "membership", "feature", "artefact", "flag"]
        with open(a.out_junctions, "w", encoding="utf-8") as fh:
            fh.write("\t".join(cols) + "\n")
            for j in sorted(all_js, key=lambda x: (x.sample, -x.support)):
                d = j.as_dict()
                fh.write("\t".join(str(d[c]) for c in cols) + "\n")

    if a.out_loci:
        others = [s for s in samples if s != ref_sample]
        with open(a.out_loci, "w", encoding="utf-8") as fh:
            hdr = ["contig", "start", "end", "flanks", "footprint",
                   "vector_pos", "region", "membership", "feature",
                   "evidence_tier", "support_" + ref_sample]
            for s in others:
                hdr += ["support_" + s, "ratio_" + s, "lo_" + s, "hi_" + s]
            fh.write("\t".join(hdr) + "\n")
            for x in sorted(loci, key=lambda L: (natural_key(L["contig"]),
                                                 L["start"])):
                j = x["rep"]
                # The tier is set by how many constructs carry the vector base
                # the junction sits on, not by the region label. For
                # integrase-mediated delivery every junction sits at a vector
                # end, so a region-based tier would mark everything "inferred"
                # even where the end sequence is unique to one vector.
                carriers = [c for c in j.membership.split(" + ") if c]
                tier = ("certain" if len(carriers) == 1
                        else "inferred" if carriers else "unknown")
                row: List[object] = [
                    x["contig"], x["start"], x["end"], x["flanks"],
                    x["footprint"] if x["footprint"] is not None else "",
                    j.vector_pos, j.region, j.membership, j.feature, tier,
                    x["support"]]
                k = (j.host_contig, j.host_pos, j.vector_pos)
                rec = xs.get(k)
                for s in others:
                    o = (rec or {}).get("others", {}).get(s)
                    if o and o["support"]:
                        row += [o["support"], round(o["ratio"], 2),
                                round(o["lo"], 2), round(o["hi"], 2)]
                    else:
                        row += ["", "", "", ""]
                fh.write("\t".join(str(v) for v in row) + "\n")

    md = render_md(a.mapping, a.modality, samples, infos, fp, foot, loci,
                   counts, warnings)
    if a.md:
        with open(a.md, "w", encoding="utf-8") as fh:
            fh.write(md)

    print("[%s] modality=%s | %s" % (a.mapping, a.modality, fp["verdict"]),
          file=sys.stderr)
    if fp.get("suggests") and MODALITIES[a.modality]["fixed_ends"] is not None:
        agrees = ((MODALITIES[a.modality]["fixed_ends"] and "FIXED" in fp["verdict"])
                  or (not MODALITIES[a.modality]["fixed_ends"]
                      and "SCATTERED" in fp["verdict"]))
        if not agrees:
            print("  WARNING: the data do not look like the declared modality "
                  "(%s). Settle this before trusting any site count."
                  % a.modality, file=sys.stderr)
    for kind in ("host-vector", "vector-vector", "vector-unmapped", "host-host"):
        print("    %-18s %4d  (%d flagged artefact)"
              % (kind, counts["kind"].get(kind, 0),
                 counts["artefact_by_kind"].get(kind, 0)), file=sys.stderr)
    for reg in ("unique", "shared", "vector-end"):
        print("    host junctions in %-14s %d" % (reg, counts["region"].get(reg, 0)),
              file=sys.stderr)
    paired = sum(1 for x in loci if x["flanks"] == 2)
    print("    loci: %d  (%d with both flanks, %d single-flank)%s"
          % (len(loci), paired, len(loci) - paired,
             "  [clip side unusable: NOTHING could be paired]"
             if not clip_usable else ""), file=sys.stderr)
    print("    footprint: %s" % foot["verdict"], file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
