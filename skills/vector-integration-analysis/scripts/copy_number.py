#!/usr/bin/env python3
"""Stage 3/4 -- copy number, architecture and site accounting.

Answers the two questions a vendor insertion-site report never answers:

  1. How much of each construct is present, per haploid genome?
  2. How many integration sites does that amount of DNA require, and does that
     number agree with the number of loci actually detected?

The agreement between those two is what turns a list of loci into a defensible
site count. If the detected loci cannot hold the measured DNA, either sites are
missing or there is an extrachromosomal pool.

Two routes to copy number, in strict order of preference
--------------------------------------------------------
**Measured (preferred).** With per-base depth over each vector reference
(`samtools depth -a`), every membership class in the identity map gives one
equation -- depth over bases carried by constructs S equals the summed depth of
those constructs -- and the whole set is solved at once by weighted least
squares. This works for any number of constructs, needs no bespoke estimator,
and is a measurement rather than an inference. One `samtools depth` command
produces it; always ask the provider for it.

**Inferred (fallback).** Without depth, mean coverage over a shared reference
is a blend of every construct carrying that sequence, so it is an **upper
bound** on the mapped construct, not its copy number. Where exactly two
constructs are in play, cross-mapping boundary junctions give their ratio and
the blend can be split. Every such estimate is labelled, and the caveats are
printed at the top of the report rather than buried in a JSON field -- an
undecomposed number that looks like a measured one is how this analysis goes
wrong quietly.

Usage
-----
    python copy_number.py config.json --md accounting.md [--json out.json]

See `assets/copy_number_config.example.json` for the schema and
`references/quantification.md` for the derivations and their assumptions.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from typing import Dict, List, Optional, Sequence, Tuple

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from via_common import (load_json, read_depth_file, solve_wls)  # noqa: E402

# Keep in step with classify_junctions.MODALITIES
CONCATEMER_CAPABLE = {"nhej": True, "hdr": True, "integrase": False,
                      "transposon-pb": False, "transposon-sb": False,
                      "auto": None}
FIXED_COPIES_PER_SITE = {"integrase": 1, "transposon-pb": 1, "transposon-sb": 1}


# --------------------------------------------------------------------------

def read_classified(path: str) -> List[Dict]:
    """Read the junction TSV written by classify_junctions.py."""
    rows: List[Dict] = []
    with open(path, encoding="utf-8") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) != len(hdr):
                continue
            d: Dict = dict(zip(hdr, f))
            for k in ("host_pos", "vector_pos", "filler", "support"):
                try:
                    d[k] = int(d.get(k, 0))
                except (TypeError, ValueError):
                    d[k] = 0
            d.setdefault("artefact", "")
            d.setdefault("kind", "")
            rows.append(d)
    if not rows:
        raise ValueError("%s: no junction rows" % path)
    return rows


def genome_coverage(dedup_pairs: float, read_length: int,
                    genome_size: float) -> float:
    return dedup_pairs * 2 * read_length / genome_size


def conversion_factor(C: float, read_length: int, anchor: int) -> float:
    """Reads expected to support one junction present at one copy per haploid."""
    return max(0.0, C * (read_length - 2 * anchor) / read_length)


# --------------------------------------------------------------------------
# measured route: solve all constructs from per-base depth
# --------------------------------------------------------------------------

def class_depths(ident: Dict, mapping: str, depth: Dict[int, float],
                 edge_trim: int, min_class_bp: int
                 ) -> List[Tuple[List[str], float, int]]:
    """Mean depth over each membership class of one vector reference.

    `edge_trim` bases are dropped from each end of every block. Within roughly
    one read length of a block edge a read cannot lie wholly inside the block,
    so cross-mapping tapers off and the observed depth there is neither the
    shared nor the unique value. Trimming is what makes the measured and
    inferred routes agree; untrimmed shared blocks are the reason naive
    decomposition comes out several-fold high.
    """
    out: List[Tuple[List[str], float, int]] = []
    for cl in ident["constructs"][mapping]["membership_classes"]:
        vals: List[float] = []
        for a, b in cl["blocks"]:
            lo, hi = a + edge_trim, b - edge_trim
            if hi - lo + 1 < 20:
                continue
            vals += [depth[p] for p in range(lo, hi + 1) if p in depth]
        if len(vals) >= min_class_bp:
            out.append((cl["membership"], sum(vals) / len(vals), len(vals)))
    return out


def solve_measured(cfg: Dict, ident: Dict, sample: str, edge_trim: int,
                   min_class_bp: int) -> Optional[Dict]:
    """Solve every construct's depth at once from per-base depth files."""
    order: List[str] = list(cfg["mappings"])
    files = (cfg["samples"][sample].get("depth_files") or {})
    rows: List[List[float]] = []
    y: List[float] = []
    w: List[float] = []
    used: List[Dict] = []

    base = cfg["_base"]
    for mapping, path in files.items():
        if mapping not in ident["constructs"]:
            continue
        p = path if os.path.isabs(path) else os.path.join(base, path)
        try:
            depth = read_depth_file(p, cfg["mappings"].get(mapping, {}).get("contig"))
        except (OSError, ValueError) as exc:
            print("  warning: %s" % exc, file=sys.stderr)
            continue
        for membership, mean, nbp in class_depths(ident, mapping, depth,
                                                  edge_trim, min_class_bp):
            row = [1.0 if n in membership else 0.0 for n in order]
            if not any(row):
                continue
            rows.append(row)
            y.append(mean)
            w.append(float(nbp))
            used.append({"mapping": mapping, "membership": membership,
                         "mean_depth": round(mean, 2), "bases": nbp})

    if len(rows) < len(order):
        return None
    x = solve_wls(rows, y, w)
    if x is None:
        return None

    resid = []
    for row, yi, u in zip(rows, y, used):
        pred = sum(r * v for r, v in zip(row, x))
        resid.append({**u, "predicted_depth": round(pred, 2),
                      "residual": round(yi - pred, 2)})
    denom = sum(wi * yi * yi for wi, yi in zip(w, y)) or 1.0
    rss = sum(wi * (yi - sum(r * v for r, v in zip(row, x))) ** 2
              for row, yi, wi in zip(rows, y, w))
    return {"method": "measured (per-base depth, weighted least squares)",
            "depths": {n: round(v, 2) for n, v in zip(order, x)},
            "equations": len(rows), "unknowns": len(order),
            "weighted_relative_rss": round(rss / denom, 4),
            "classes_used": resid}


# --------------------------------------------------------------------------
# inferred route: boundary-junction estimator, or an honest upper bound
# --------------------------------------------------------------------------

def boundary_signal(rows: Sequence[Dict], sample: str) -> Optional[float]:
    """Mean support of cross-mapping boundary junctions in one sample.

    These are made by reads belonging to a construct that was NOT in the
    mapping reference, so their number is proportional to the copy number of
    that other construct -- which is what makes them useful.
    """
    vals = [r["support"] for r in rows
            if r["sample"] == sample
            and str(r.get("artefact", "")).startswith("boundary-crossmap")]
    return sum(vals) / len(vals) if vals else None


def shared_fraction(ident: Dict, name: str, other: str) -> float:
    c = ident["constructs"][name]
    bp = sum(cl["bp"] for cl in c["membership_classes"]
             if other in cl["membership"])
    return bp / c["length"]


def solve_inferred(cfg: Dict, ident: Dict, sample: str,
                   jrows: Dict[str, List[Dict]]) -> Dict:
    order = list(cfg["mappings"])
    caveats: List[str] = []

    if len(order) == 1:
        n = order[0]
        mean = cfg["mappings"][n].get("mean_coverage", {}).get(sample)
        if mean is None:
            return {"method": "none", "depths": {},
                    "caveats": ["no depth file and no mean_coverage for sample "
                                "%r: copy number cannot be estimated" % sample]}
        return {"method": "inferred (mean coverage, single construct)",
                "depths": {n: float(mean)},
                "caveats": ["Single construct, so mean coverage is its depth "
                            "directly. This still assumes the vector contains "
                            "no sequence shared with the host genome; check the "
                            "identity map against the host if the vector "
                            "carries genomic elements."]}

    if len(order) == 2:
        A, B = order
        sig_A = boundary_signal(jrows.get(A, []), sample)
        sig_B = boundary_signal(jrows.get(B, []), sample)
        meanA = cfg["mappings"][A].get("mean_coverage", {}).get(sample)
        if sig_A and sig_B and meanA:
            r = sig_B / sig_A                     # copies(A) / copies(B)
            fA = shared_fraction(ident, A, B)
            PA = meanA / (1 + fA / r)
            return {"method": "inferred (boundary-junction estimator, 2 constructs)",
                    "depths": {A: round(PA, 2), B: round(PA / r, 2)},
                    "ratio": round(r, 3),
                    "boundary_signal": {A: round(sig_A, 1), B: round(sig_B, 1)},
                    "caveats": [
                        "Copy numbers come from the boundary-junction "
                        "estimator, not from per-base depth. The RATIO (%.2f) "
                        "is the robust quantity; the absolute values are a "
                        "lower bound, because the effective cross-mappable "
                        "fraction is smaller than the nominal shared fraction. "
                        "Build conclusions on the ratio, and request "
                        "`samtools depth -a` to convert this from inferred to "
                        "measured." % r]}
        caveats.append(
            "The boundary-junction estimator needs cross-mapping boundary "
            "junctions in BOTH mappings and a mean coverage figure; one or more "
            "was missing for sample %r." % sample)

    depths: Dict[str, float] = {}
    for n in order:
        mean = cfg["mappings"][n].get("mean_coverage", {}).get(sample)
        if mean is not None:
            depths[n] = float(mean)
    caveats.append(
        "NOT DECOMPOSED. With %d constructs sharing sequence, mean coverage "
        "over one reference is the summed depth of every construct carrying "
        "that sequence. Each figure below is therefore an UPPER BOUND on its "
        "construct, inflated roughly in proportion to the shared fraction, and "
        "every number derived from it -- array length, implied sites, flank "
        "detection -- inherits that inflation. Supply per-base depth "
        "(`samtools depth -a`) for each vector reference to replace this with a "
        "solved measurement." % len(order))
    return {"method": "NOT DECOMPOSED (upper bound from mean coverage)",
            "depths": depths, "caveats": caveats}


# --------------------------------------------------------------------------
# architecture and site accounting
# --------------------------------------------------------------------------

def architecture(cfg: Dict, jrows: Dict[str, List[Dict]], mapping: str,
                 sample: str, copies: float, modality: str,
                 include_end_religation: bool) -> Dict:
    """Copies per site, and the site count that implies."""
    rows = jrows.get(mapping, [])
    out_notes: List[str] = []
    present = {r["sample"] for r in rows}
    if rows and sample not in present:
        # Silent zeros here would set the read counts to nothing and make
        # copies-per-site and the whole reconciliation meaningless while still
        # printing numbers, so this is called out rather than absorbed.
        out_notes.append(
            "Sample %r does not appear in the junction table for %s, which "
            "holds sample(s) %s. Flank and internal read counts are therefore "
            "ZERO and every quantity derived from them below is meaningless. "
            "The sample keys in this config must match the names passed to "
            "classify_junctions.py as --table NAME=PATH."
            % (sample, mapping, ", ".join(sorted(present)) or "none"))
    flank = sum(r["support"] for r in rows
                if r["sample"] == sample and r["kind"] == "host-vector"
                and not r["artefact"])

    def internal_ok(r: Dict) -> bool:
        if r["sample"] != sample or r["kind"] != "vector-vector":
            return False
        if not r["artefact"]:
            return True
        return (include_end_religation
                and str(r["artefact"]).startswith("vector-end-religation"))

    internal = sum(r["support"] for r in rows if internal_ok(r))
    out: Dict = {"flank_reads": flank, "internal_reads": internal,
                 "end_religation_included": include_end_religation,
                 "modality": modality, "notes": out_notes}

    capable = CONCATEMER_CAPABLE.get(modality)
    if capable is False:
        n = float(FIXED_COPIES_PER_SITE.get(modality, 1))
        out["copies_per_site"] = n
        out["copies_per_site_source"] = (
            "fixed at %g by the delivery mechanism (%s). The internal:flank "
            "ratio estimator was NOT applied: it assumes a tandem array, and "
            "this mechanism inserts one copy per site. Applying it here would "
            "divide the copy number by a spurious array length and understate "
            "the site count several-fold." % (n, modality))
        if internal:
            out["notes"].append(
                "%d reads support vector-vector junctions even though this "
                "mechanism inserts single copies. That is worth explaining: "
                "candidates are a rearranged or internally deleted insert, "
                "two insertions adjacent in the genome, carried-over producer "
                "plasmid, or unintegrated vector." % internal)
    elif capable is True and flank:
        out["copies_per_site"] = round(1 + 2 * (internal / flank), 2)
        out["copies_per_site_source"] = (
            "internal:flank read ratio. An array of n copies at one site has "
            "n-1 internal junctions and 2 host flanks, so "
            "n = 1 + 2*(internal/flank). This is an average over all sites and "
            "cannot tell one 20-copy array from twenty single copies.")
    else:
        out["copies_per_site"] = None
        out["copies_per_site_source"] = (
            "cannot be estimated: " + ("no host-vector flank reads" if not flank
                                       else "modality undeclared, so it is "
                                       "unknown whether arrays are possible"))

    n = out["copies_per_site"]
    out["implied_sites"] = round(copies / n, 2) if n else None
    return out


def flank_detection(cfg: Dict, C: float, flank_reads: int,
                    implied_sites: Optional[float],
                    anchors: Sequence[int]) -> List[Dict]:
    sweep = []
    for anchor in anchors:
        R = conversion_factor(C, cfg["read_length"], anchor)
        detected = flank_reads / R if R else None
        sweep.append({
            "anchor": anchor,
            "reads_per_copy": round(R, 2),
            "flank_copies_detected": round(detected, 2) if detected else None,
            "flank_detection_rate": (round(detected / (2 * implied_sites), 3)
                                     if detected and implied_sites else None),
        })
    return sweep


def reconcile(implied: Optional[float], detected_loci: Optional[int],
              tol: float = 0.25) -> Dict:
    if implied is None or detected_loci is None:
        return {"verdict": "not reconciled: implied site count or detected "
                           "locus count unavailable"}
    if detected_loci == 0:
        return {"verdict": "no loci detected, so nothing to reconcile"}
    ratio = implied / detected_loci
    if abs(ratio - 1) <= tol:
        v = ("AGREES within %d%%: the detected loci account for the DNA "
             "present, so there is no large undetected reservoir."
             % int(100 * tol))
    elif ratio > 1:
        v = ("SHORTFALL: the measured DNA needs %.1f sites but only %d were "
             "detected (%.1fx). Either sites are missing -- flanks lost to "
             "repeats are the usual reason -- or part of the vector is "
             "extrachromosomal. Test with a plasmid-safe nuclease or by "
             "passaging without selection." % (implied, detected_loci, ratio))
    else:
        v = ("EXCESS: %d loci were detected but the measured DNA only supports "
             "%.1f (%.1fx). Loci are probably being double-counted -- check "
             "flank pairing first, then whether one site is being reported "
             "under two nearby coordinates." % (detected_loci, implied, ratio))
    return {"implied_sites": implied, "detected_loci": detected_loci,
            "ratio": round(ratio, 2), "verdict": v}


def calibrate(cfg: Dict, jrows: Dict[str, List[Dict]], C: float,
              sample: str, anchors: Sequence[int]) -> Optional[Dict]:
    cal = cfg.get("calibration")
    if not cal:
        return None
    contig, _, pos = cal["locus"].partition(":")
    if not pos.strip().isdigit():
        return {"error": "calibration.locus must be CONTIG:POSITION, got %r"
                         % cal["locus"]}
    rows = jrows.get(cal["mapping"], [])
    hits = [r["support"] for r in rows
            if r["sample"] == sample and r["kind"] == "host-vector"
            and r["host_contig"] == contig
            and abs(r["host_pos"] - int(pos)) <= cal.get("window", 50)]
    if not hits:
        return {"locus": cal["locus"], "error": "no junction found at this locus "
                                                "in sample %r" % sample}
    mean = sum(hits) / len(hits)
    est = []
    for anchor in anchors:
        R = conversion_factor(C, cfg["read_length"], anchor)
        est.append({"anchor": anchor, "reads_per_copy": round(R, 2),
                    "implied_copies": round(mean / R, 3) if R else None})
    exp = cal["expected_copies_per_haploid"]
    mid = est[len(est) // 2]["implied_copies"]
    verdict = "no estimate"
    if mid:
        f = mid / exp if exp else float("inf")
        verdict = (
            "within 30%%: the conversion factor, the anchor assumption and the "
            "ploidy estimate are all consistent" if abs(f - 1) <= 0.3 else
            "off by %.1fx. One of three things is wrong: the anchor length, the "
            "ploidy, or the assumed structure of this locus. Ploidy is the "
            "usual culprit and the calibration can be inverted to estimate it "
            "(implied ploidy ~ %.1f)." % (f, cfg.get("ploidy", 2) * f))
    return {"locus": cal["locus"], "flanks_found": len(hits),
            "mean_support": round(mean, 1),
            "expected_copies_per_haploid": exp,
            "estimates": est, "verdict": verdict}


# --------------------------------------------------------------------------

def render(cfg: Dict, result: Dict) -> str:
    primary = result["primary_mapping"]
    L = ["# Copy number and site accounting", "",
         "Primary construct: **%s**. Declared modality: **%s**."
         % (primary, result["modality"]), ""]

    caveats = result["caveats"]
    if caveats:
        L += ["## Read this before using any number below", ""]
        L += ["> %s" % c for c in caveats]
        L.append("")
    else:
        L += ["All copy numbers below are solved from per-base depth, not "
              "inferred.", ""]

    L += ["## Construct abundance", "",
          "Units: copies per haploid genome equivalent. Multiply by the ploidy "
          "(%s) for copies per cell, and treat ploidy as measured rather than "
          "assumed -- many immortalised lines are aneuploid and drift."
          % cfg.get("ploidy", "?"), ""]
    order = list(cfg["mappings"])
    L += ["| Sample | Genome cov | " + " | ".join(order) + " | Route |",
          "|---" * (len(order) + 3) + "|"]
    for s, d in result["samples"].items():
        cells = [s, "%.1fx" % d["genome_coverage"]]
        for n in order:
            c = d["copies_per_haploid"].get(n)
            cells.append("%.2f" % c if c is not None else "-")
        cells.append(d["solve"]["method"])
        L.append("| " + " | ".join(cells) + " |")
    L.append("")

    for s, d in result["samples"].items():
        sol = d["solve"]
        if sol.get("classes_used"):
            L += ["### %s -- the solve" % s, "",
                  "%d equations, %d unknowns, weighted relative residual %s."
                  % (sol["equations"], sol["unknowns"],
                     sol["weighted_relative_rss"]),
                  "",
                  "| Reference | Bases carried by | bp used | observed | predicted | residual |",
                  "|---|---|---|---|---|---|"]
            for c in sol["classes_used"]:
                L.append("| %s | %s | %s | %.2f | %.2f | %+.2f |"
                         % (c["mapping"], " + ".join(c["membership"]),
                            format(c["bases"], ","), c["mean_depth"],
                            c["predicted_depth"], c["residual"]))
            L += ["", "Large residuals mean the block map does not explain the "
                  "depth: a repeat shared with the host genome, an unlisted "
                  "construct, or non-uniform library coverage.", ""]

    L += ["## Architecture and site accounting (%s)" % primary, "",
          "| Sample | copies/haploid | flank reads | internal reads | "
          "copies/site | implied sites | loci detected |",
          "|---|---|---|---|---|---|---|"]
    for s, d in result["samples"].items():
        ar = d["architecture"]
        L.append("| %s | %s | %d | %d | %s | %s | %s |"
                 % (s,
                    "%.2f" % d["copies_per_haploid"][primary]
                    if d["copies_per_haploid"].get(primary) is not None else "-",
                    ar["flank_reads"], ar["internal_reads"],
                    ar["copies_per_site"] if ar["copies_per_site"] else "-",
                    ar["implied_sites"] if ar["implied_sites"] else "-",
                    d.get("detected_loci", "-")))
    L.append("")
    first = next(iter(result["samples"].values()), None)
    if first:
        L += ["**Copies per site:** %s"
              % first["architecture"]["copies_per_site_source"], ""]
        for n in first["architecture"]["notes"]:
            L += ["> %s" % n, ""]

    L += ["## Reconciliation -- this is the result", ""]
    for s, d in result["samples"].items():
        L.append("- **%s**: %s" % (s, d["reconciliation"]["verdict"]))
    L += ["", "A flank detection rate well below 100% is normal and needs no "
          "exotic explanation: roughly half of a mammalian genome is "
          "repetitive, so a flank landing in a LINE, SINE, satellite or "
          "segmental duplication cannot be placed uniquely; an array edge may "
          "join more vector rather than returning to the chromosome; and a "
          "junction whose filler exceeds what a read can span is never called.",
          ""]
    L += ["### Sensitivity to the aligner anchor", "",
          "The minimum anchor is a property of the junction caller and is "
          "rarely documented, so it is swept rather than assumed.", "",
          "| Sample | anchor | reads/copy | flank copies detected | detection rate |",
          "|---|---|---|---|---|"]
    for s, d in result["samples"].items():
        for row in d["flank_detection"]:
            L.append("| %s | %d bp | %s | %s | %s |"
                     % (s, row["anchor"], row["reads_per_copy"],
                        row["flank_copies_detected"] if row["flank_copies_detected"] else "-",
                        row["flank_detection_rate"] if row["flank_detection_rate"] else "-"))
    L.append("")

    cals = {s: d["calibration"] for s, d in result["samples"].items()
            if d.get("calibration")}
    if cals:
        L += ["## Calibration against a locus of known structure", ""]
        for s, c in cals.items():
            if c.get("error"):
                L.append("- **%s**: %s" % (s, c["error"]))
                continue
            L += ["- **%s**, locus `%s`: %d flank(s), mean support %s, "
                  "expected %s copies/haploid -> %s"
                  % (s, c["locus"], c["flanks_found"], c["mean_support"],
                     c["expected_copies_per_haploid"], c["verdict"])]
        L.append("")
    else:
        L += ["## Calibration", "",
              "No calibration locus supplied. If any locus has independently "
              "known structure -- both flanks recovered in sequence unique to "
              "one construct, consistent with a single insertion on one "
              "homologue -- put it in the config. It validates the conversion "
              "factor, the anchor assumption and the ploidy estimate at once, "
              "and costs nothing.", ""]
    return "\n".join(L)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("config")
    ap.add_argument("--md")
    ap.add_argument("--json")
    ap.add_argument("--edge-trim", type=int, default=0,
                    help="bases trimmed from each end of every block before "
                         "averaging depth. Default 0 = use read_length.")
    ap.add_argument("--min-class-bp", type=int, default=200,
                    help="membership classes with fewer usable bases than this "
                         "are left out of the solve (default 200)")
    ap.add_argument("--anchors", default="15,20,25",
                    help="aligner anchor lengths to sweep")
    ap.add_argument("--include-end-religation", action="store_true",
                    help="count junctions joining the two vector ends as "
                         "internal array junctions. OFF by default: against a "
                         "linear reference they are also produced by uncut "
                         "input molecule and by any co-resident construct whose "
                         "equivalent site is intact, so they are ambiguous.")
    a = ap.parse_args(argv)

    try:
        cfg = load_json(a.config)
    except (OSError, ValueError) as exc:
        print("error: cannot read config: %s" % exc, file=sys.stderr)
        return 2
    cfg = {k: v for k, v in cfg.items() if not k.startswith("_")}
    base = os.path.dirname(os.path.abspath(a.config))
    cfg["_base"] = base

    for key in ("identity", "genome_size", "read_length", "samples", "mappings"):
        if key not in cfg:
            print("error: config is missing %r" % key, file=sys.stderr)
            return 2

    try:
        ident = load_json(cfg["identity"] if os.path.isabs(cfg["identity"])
                          else os.path.join(base, cfg["identity"]))
    except (OSError, ValueError) as exc:
        print("error: cannot read identity map: %s" % exc, file=sys.stderr)
        return 2

    unknown = [n for n in cfg["mappings"] if n not in ident["constructs"]]
    if unknown:
        print("error: mapping(s) %s are not in the identity map (%s). The "
              "identity map must cover every construct in the line, including "
              "ones the provider did not map against."
              % (unknown, ", ".join(ident["constructs"])), file=sys.stderr)
        return 2

    jrows: Dict[str, List[Dict]] = {}
    for name, m in cfg["mappings"].items():
        p = m.get("junctions")
        if not p:
            continue
        try:
            jrows[name] = read_classified(
                p if os.path.isabs(p) else os.path.join(base, p))
        except (OSError, ValueError) as exc:
            print("error: %s" % exc, file=sys.stderr)
            return 2

    primary = cfg.get("primary_mapping", list(cfg["mappings"])[0])
    if primary not in cfg["mappings"]:
        print("error: primary_mapping %r is not one of %s"
              % (primary, ", ".join(cfg["mappings"])), file=sys.stderr)
        return 2
    modality = cfg.get("modality", "auto")
    if modality not in CONCATEMER_CAPABLE:
        print("error: modality %r unknown; expected one of %s"
              % (modality, ", ".join(sorted(CONCATEMER_CAPABLE))), file=sys.stderr)
        return 2

    edge_trim = a.edge_trim or int(cfg["read_length"])
    anchors = [int(x) for x in a.anchors.split(",") if x.strip()]

    result: Dict = {"primary_mapping": primary, "modality": modality,
                    "samples": {}, "caveats": []}
    caveats: List[str] = []
    if modality == "auto":
        caveats.append(
            "The delivery modality is undeclared. Whether one site carries one "
            "copy or a tandem array is a property of the mechanism, not of the "
            "data, and it changes the site count several-fold. Settle it from "
            "the Stage 2 fingerprint and set `modality` in the config.")

    for sample, s in cfg["samples"].items():
        if "dedup_pairs" not in s:
            print("error: sample %r has no dedup_pairs" % sample, file=sys.stderr)
            return 2
        C = genome_coverage(float(s["dedup_pairs"]), int(cfg["read_length"]),
                            float(cfg["genome_size"]))
        sol = solve_measured(cfg, ident, sample, edge_trim, a.min_class_bp)
        if sol is None:
            sol = solve_inferred(cfg, ident, sample, jrows)
        for c in sol.get("caveats", []):
            if c not in caveats:
                caveats.append("[%s] %s" % (sample, c) if len(cfg["samples"]) > 1 else c)

        copies = {n: (v / C if C else None) for n, v in sol["depths"].items()}
        cn = copies.get(primary)
        arch = architecture(cfg, jrows, primary, sample, cn if cn else 0.0,
                            modality, a.include_end_religation)
        for n in arch["notes"]:
            if "does not appear in the junction table" in n and n not in caveats:
                caveats.append(n)
        detected = s.get("detected_loci")
        if detected is None:
            detected = cfg["mappings"][primary].get("detected_loci")
        result["samples"][sample] = {
            "genome_coverage": round(C, 2),
            "solve": sol,
            "depths": sol["depths"],
            "copies_per_haploid": {n: (round(v, 3) if v is not None else None)
                                   for n, v in copies.items()},
            "copies_per_cell": {n: (round(v * float(cfg.get("ploidy", 2)), 3)
                                    if v is not None else None)
                                for n, v in copies.items()},
            "architecture": arch,
            "detected_loci": detected,
            "reconciliation": reconcile(arch["implied_sites"], detected),
            "flank_detection": flank_detection(cfg, C, arch["flank_reads"],
                                               arch["implied_sites"], anchors),
            "calibration": calibrate(cfg, jrows, C, sample, anchors),
        }

    result["caveats"] = caveats
    md = render(cfg, result)
    if a.md:
        with open(a.md, "w", encoding="utf-8") as fh:
            fh.write(md)
    if a.json:
        with open(a.json, "w", encoding="utf-8") as fh:
            json.dump(result, fh, indent=2)
    if not a.md and not a.json:
        print(md)

    for c in caveats:
        print("  CAVEAT: %s" % c, file=sys.stderr)
    for s, d in result["samples"].items():
        parts = " | ".join("%s %s" % (n, v) for n, v in
                           d["copies_per_haploid"].items())
        print("%s: genome %.1fx | %s | route: %s"
              % (s, d["genome_coverage"], parts, d["solve"]["method"]),
              file=sys.stderr)
        ar = d["architecture"]
        print("    copies/site %s -> %s sites implied | %s"
              % (ar["copies_per_site"], ar["implied_sites"],
                 d["reconciliation"]["verdict"].split(":")[0]), file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
