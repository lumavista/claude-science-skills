#!/usr/bin/env python3
"""Final stage -- reconcile the per-mapping loci into ONE called-site table.

This is the authoritative site producer. Report its output; do not re-assemble
loci, re-apply thresholds or re-argue attribution by hand. Every step is
deterministic and every threshold is a parameter, so the same inputs always
give the same sites.

Inputs
------
* one `loci_<mapping>.tsv` per mapping run, from classify_junctions.py
* the copy_number.py JSON, for each construct's copies per haploid
  (or explicit `--copies CONSTRUCT=N`)
* the identity map, `--vector-ends CONSTRUCT=a,b` per construct, the ploidy,
  and the retention parameters

What it does
------------
1. COLLAPSE across mappings. One physical host-vector junction is reported by
   every mapping run that saw it, at the same host coordinate. Entries from
   DIFFERENT mappings whose host intervals lie within `--window` are linked,
   and each linked group is ONE unified locus -- never two sites. Support is
   not summed across mappings, since that is the same reads counted twice; per
   sample it is the largest value any single mapping reports.

2. RETAIN. A unified locus is called only if its support is >= `--min-reads`
   AND it is supported in >= `--min-samples` independently prepared samples.
   Everything else goes to the sub-threshold table as the single-read/noise
   tail. It is never in the called set.

3. ATTRIBUTE, by strict priority -- this is the judgement being removed from
   the caller:

     a. unique       the vector side lies, in some mapping, in sequence only
                     one construct carries -> that construct, `certain`.
     b. vector-end   the junction sits within `--end-tol` of exactly one
                     construct's declared end -> that construct, `inferred`.
     c. copy-ceiling no unique and no end signal. Each construct X admits at
                     most ceil(copies_X * ploidy / --min-copies-per-site)
                     sites, and that cap is HARD. Among the candidates whose
                     cap is not yet full, the locus goes to the highest-copy
                     one only when it leads the next by a factor of
                     `--cn-margin`; otherwise nobody gets it.
     d. unresolved   everything left: listed ONCE, with its candidates named,
                     attributed to no construct.

   No support ratio and no cross-sample scaling enters this, so neither can
   override the ceiling. A construct measured at N copies cannot be given more
   copy-ceiling sites than its ceiling allows, whatever the support looks like.

4. COUNT, per construct: a floor of certain + vector-end loci, a bounded range
   up to the copy-number ceiling, and the copy-ceiling and unresolved counts
   stated separately. One range per construct -- not per sample, not per
   mapping.

Usage
-----
    python reconcile_sites.py --identity identity.json \\
        --loci pAlpha=loci_pAlpha.tsv --loci pBeta=loci_pBeta.tsv \\
        --copy-number accounting.json --ploidy 2 \\
        --vector-ends pAlpha=1,8412 --vector-ends pBeta=1,7650 \\
        --out-called called_sites.tsv \\
        --out-subthreshold subthreshold_sites.tsv \\
        --summary site_summary.json --md site_summary.md

With `--run run.json` the identity map, the loci tables, the vector ends, the
copy-number JSON and the ploidy are taken from the run manifest that the
earlier stages recorded them in, and the four outputs are named for you:

    python reconcile_sites.py --run run.json --strict

Rules, priority order and the ceiling argument: `references/attribution.md`.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
from typing import Dict, FrozenSet, List, Optional, Sequence, Tuple

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from via_common import (MODALITIES, load_json, natural_key,  # noqa: E402
                        parse_int_list, parse_kv, read_loci)
from via_run import (FAIL, PASS, WARN, Checks, InputError,  # noqa: E402
                     RunError, add_check_arguments, load_run,
                     reconcile_inputs)


# --------------------------------------------------------------------------
# input
# --------------------------------------------------------------------------

def parse_ends(items: Sequence[str]) -> Dict[str, List[int]]:
    """`--vector-ends CONSTRUCT=a,b`, one per construct."""
    return {name: parse_int_list(spec, "--vector-ends", name)
            for name, spec in parse_kv(items, "--vector-ends").items()}


def load_copies(path: Optional[str], explicit: Dict[str, str],
                cn_sample: Optional[str]
                ) -> Tuple[Dict[str, Optional[float]], Optional[float], str]:
    """Per-construct copies per haploid, and the ploidy if the JSON implies it.

    With several samples and no `--cn-sample`, each construct's MAXIMUM across
    samples is used: the called set is a union over samples, so the ceiling has
    to admit every site any one sample could hold. `--copies` overrides the
    JSON for the constructs it names.
    """
    copies: Dict[str, Optional[float]] = {}
    ploidy: Optional[float] = None
    how = "no copy number supplied"
    if path:
        d = load_json(path)
        samples = d.get("samples") or {}
        if not samples:
            raise ValueError("%s: no `samples` block; is this a copy_number.py "
                             "--json output?" % path)
        if cn_sample:
            if cn_sample not in samples:
                raise ValueError("--cn-sample %r is not in %s (%s)"
                                 % (cn_sample, path, ", ".join(samples)))
            use = [cn_sample]
            how = "copies per haploid from sample %s" % cn_sample
        else:
            use = list(samples)
            how = ("copies per haploid from sample %s" % use[0]
                   if len(use) == 1 else
                   "copies per haploid: per-construct maximum over samples "
                   "%s (the called set is a union over samples, so the "
                   "ceiling must admit every site any one of them holds)"
                   % ", ".join(use))
        for s in use:
            cph = samples[s].get("copies_per_haploid") or {}
            cpc = samples[s].get("copies_per_cell") or {}
            for n, v in cph.items():
                if v is None:
                    copies.setdefault(n, None)
                    continue
                prev = copies.get(n)
                copies[n] = float(v) if prev is None else max(prev, float(v))
                if ploidy is None and cpc.get(n) and float(v) > 0:
                    ploidy = round(float(cpc[n]) / float(v), 6)
    for n, v in explicit.items():
        try:
            copies[n] = float(v)
        except ValueError:
            raise ValueError("--copies needs CONSTRUCT=NUMBER, got %r=%r"
                             % (n, v))
    if explicit:
        how += ("; overridden from --copies for %s"
                % ", ".join(sorted(explicit)))
    return copies, ploidy, how


# --------------------------------------------------------------------------
# 1. collapse across mappings
# --------------------------------------------------------------------------

def _gap(a: Dict, b: Dict) -> int:
    """Distance between two host intervals on one contig; 0 if they overlap."""
    return max(0, max(a["start"], b["start"]) - min(a["end"], b["end"]))


class _UnionFind:
    """Union-find over entry indices: the smallest thing that expresses
    "these reports are the same physical junction"."""

    def __init__(self, n: int) -> None:
        self.parent = list(range(n))

    def find(self, i: int) -> int:
        while self.parent[i] != i:
            self.parent[i] = self.parent[self.parent[i]]
            i = self.parent[i]
        return i

    def union(self, i: int, j: int) -> None:
        a, b = self.find(i), self.find(j)
        if a != b:
            self.parent[max(a, b)] = min(a, b)

    def groups(self, items: Sequence) -> List[List]:
        out: Dict[int, List] = {}
        for i, x in enumerate(items):
            out.setdefault(self.find(i), []).append(x)
        return list(out.values())


def _link_across_mappings(entries: Sequence[Dict], window: int) -> _UnionFind:
    """Link entries of DIFFERENT mappings whose host intervals are within
    `window`. Entries of the same mapping are never linked to each other:
    classify_junctions.py has already decided whether they are one locus or
    two, and re-merging them here would deflate the site count."""
    g = _UnionFind(len(entries))
    by_contig: Dict[str, List[int]] = {}
    for i, e in enumerate(entries):
        by_contig.setdefault(e["contig"], []).append(i)
    for idx in by_contig.values():
        idx.sort(key=lambda i: (entries[i]["start"], entries[i]["mapping"]))
        for a in range(len(idx)):
            for b in range(a + 1, len(idx)):
                ea, eb = entries[idx[a]], entries[idx[b]]
                if eb["start"] - ea["end"] > window:
                    break
                if ea["mapping"] != eb["mapping"] and _gap(ea, eb) <= window:
                    g.union(idx[a], idx[b])
    return g


def _unify(members: List[Dict], samples: Sequence[str]) -> Dict:
    """One unified locus from the entries that report it.

    Support is combined per sample as the MAXIMUM over entries, never the sum:
    the mappings report the same reads, so summing inflates every shared site
    in proportion to the number of mappings.
    """
    members.sort(key=lambda e: (e["mapping"], e["start"]))
    support = {s: max(e["support"].get(s, 0) for e in members) for s in samples}
    maps = sorted({e["mapping"] for e in members})
    foot = [e["footprint"] for e in members if e["footprint"] is not None]
    return {
        "contig": members[0]["contig"],
        "start": min(e["start"] for e in members),
        "end": max(e["end"] for e in members),
        "flanks": max(e["flanks"] for e in members),
        "footprint": foot[0] if foot else None,
        "entries": members,
        "mappings": maps,
        "same_mapping_merge": len(members) > len(maps),
        "support": support,
        "total_support": sum(support.values()),
        "n_samples": sum(1 for v in support.values() if v > 0),
    }


def collapse(entries: Sequence[Dict], window: int, samples: Sequence[str]
             ) -> List[Dict]:
    """Collapse the per-mapping entries into unified loci: one per junction.

    The same physical junction appears once per mapping run that detected it,
    at the same host coordinate but a different vector coordinate. Linking
    those into one locus is what stops a two-construct line reporting every
    shared-backbone site twice, once under each construct.

    A group can still hold two entries of one mapping by being linked through
    a third mapping's entry; that is recorded as `same_mapping_merge` so it is
    visible rather than silent.
    """
    groups = _link_across_mappings(entries, window)
    loci = [_unify(members, samples)
            for members in groups.groups(list(entries))]
    loci.sort(key=lambda locus: (natural_key(locus["contig"]), locus["start"]))
    return loci


# --------------------------------------------------------------------------
# 2. retain
# --------------------------------------------------------------------------

def retain(loci: Sequence[Dict], min_reads: int, min_samples: int
           ) -> Tuple[List[Dict], List[Dict]]:
    """Split the unified loci into the called set and the sub-threshold tail.

    Two independent conditions, because they fail for different reasons: read
    support separates signal from the single-read noise tail, and the sample
    count asks whether the event survived an independent library preparation.
    A locus seen once, in one sample, at one read is the shape of a library
    chimera; requiring both is what keeps that tail out of a site count.
    """
    called: List[Dict] = []
    sub: List[Dict] = []
    for locus in loci:
        why: List[str] = []
        if locus["total_support"] < min_reads:
            why.append("support %d < --min-reads %d"
                       % (locus["total_support"], min_reads))
        if locus["n_samples"] < min_samples:
            why.append("in %d sample(s) < --min-samples %d"
                       % (locus["n_samples"], min_samples))
        if why:
            locus["subthreshold_reason"] = "; ".join(why)
            sub.append(locus)
        else:
            called.append(locus)
    return called, sub


# --------------------------------------------------------------------------
# 3. attribute
# --------------------------------------------------------------------------

def ceilings(copies: Dict[str, Optional[float]], constructs: Sequence[str],
             ploidy: float, min_copies_per_site: float
             ) -> Dict[str, Optional[int]]:
    """How many sites each construct's measured copy number can source.

        cap_X = ceil(copies_per_haploid_X * ploidy / min_copies_per_site)

    `min_copies_per_site` is the fewest copies a real site is assumed to
    carry: 1 by default, which is the most permissive honest value, and the
    only one valid for integrase or transposase delivery. Raising it for a
    concatemer line (where each site holds an array) tightens every cap.

    A construct with no reported copy number gets `None`: unknown, which is
    not the same as unlimited. The copy-ceiling step refuses to attribute to
    it at all, because the whole point of the step is that the copy number
    is what sources the site.
    """
    caps: Dict[str, Optional[int]] = {}
    for n in constructs:
        v = copies.get(n)
        caps[n] = (None if v is None else
                   int(math.ceil(max(0.0, v) * ploidy / min_copies_per_site)))
    return caps


def candidates(locus: Dict, constructs: Sequence[str]) -> FrozenSet[str]:
    """Constructs that carry the vector sequence this junction sits in.

    The union over the locus's per-mapping entries, so a construct named by
    any mapping stays a candidate: over-inclusive on purpose, since the cost
    of dropping a true candidate is a confidently wrong attribution, while
    the cost of keeping a false one is an honest `unresolved`.
    """
    out: set = set()
    for e in locus["entries"]:
        out |= {m for m in e["membership"] if m in constructs}
    return frozenset(out)


def unique_attribution(locus: Dict,
                       constructs: Sequence[str]) -> Optional[str]:
    """(a) CERTAIN: some mapping puts the vector side in one construct's own
    sequence. Singleton membership in the identity map is sequence evidence,
    and nothing downstream may overrule it."""
    for e in locus["entries"]:
        mem = [m for m in e["membership"] if m in constructs]
        if len(mem) == 1 and len(e["membership"]) == 1:
            return mem[0]
    return None


def end_attribution(locus: Dict, ends: Dict[str, List[int]], end_tol: int
                    ) -> Tuple[Optional[str], List[str], List[str]]:
    """(b) VECTOR-END: at exactly one construct's declared end.

    Each entry's `vector_pos` is in its own mapping's coordinates, so a
    construct's end can only be tested where that construct was mapped
    against. Returns the attribution, the constructs whose end matched, and
    the constructs that could be tested at all -- an untested construct is
    reported, not assumed away.
    """
    matched: List[str] = []
    tested: List[str] = []
    for e in locus["entries"]:
        m = e["mapping"]
        if m not in ends or not ends[m] or e["vector_pos"] is None:
            continue
        tested.append(m)
        if any(abs(e["vector_pos"] - x) <= end_tol for x in ends[m]):
            matched.append(m)
    matched = sorted(set(matched))
    tested = sorted(set(tested))
    return (matched[0] if len(matched) == 1 else None), matched, tested


def ceiling_attribution(cands: Sequence[str], copies: Dict[str, Optional[float]],
                        caps: Dict[str, Optional[int]], used: Dict[str, int],
                        margin: float
                        ) -> Tuple[Optional[str], str]:
    """(c) COPY-CEILING, the hard cap. Returns (construct or None, reason).

    Admissible = a candidate with a known copy number whose cap is not yet
    full. The locus goes to the highest-copy admissible candidate only when
    it leads the next by at least `margin` times; a construct whose cap is
    full is not admissible at all, so no amount of read support can put a
    site on it. That is the whole point: a construct measured at N copies
    cannot be given more sites than its ceiling allows.
    """
    known = [c for c in cands if caps.get(c) is not None]
    unknown = [c for c in cands if caps.get(c) is None]
    if not known:
        return None, ("no candidate has a reported copy number (%s), so the "
                      "ceiling cannot source this locus"
                      % (", ".join(sorted(cands)) or "none"))
    full = [c for c in known if used.get(c, 0) >= (caps[c] or 0)]
    open_ = [c for c in known if c not in full]
    if not open_:
        return None, ("every candidate's copy-number ceiling is already full "
                      "(%s), so no construct can source this locus"
                      % ", ".join("%s %d/%d" % (c, used.get(c, 0), caps[c] or 0)
                                  for c in sorted(full)))
    ranked = sorted(open_, key=lambda c: (-(copies.get(c) or 0.0), c))
    top = ranked[0]
    top_cn = copies.get(top) or 0.0
    note_full = ("; %s excluded, ceiling full"
                 % ", ".join(sorted(full)) if full else "")
    note_unknown = ("; %s excluded, no copy number"
                    % ", ".join(sorted(unknown)) if unknown else "")
    if len(ranked) == 1:
        return top, ("only %s is admissible under the ceilings%s%s"
                     % (top, note_full, note_unknown))
    nxt = ranked[1]
    nxt_cn = copies.get(nxt) or 0.0
    if top_cn >= margin * max(nxt_cn, 1e-12):
        return top, ("%s at %.2f copies leads %s at %.2f by >= %.2gx "
                     "(--cn-margin)%s%s"
                     % (top, top_cn, nxt, nxt_cn, margin, note_full,
                        note_unknown))
    return None, ("%s at %.2f and %s at %.2f copies are within %.2gx "
                  "(--cn-margin), so neither is preferred%s%s"
                  % (top, top_cn, nxt, nxt_cn, margin, note_full, note_unknown))


def _init_tiers(called: Sequence[Dict], constructs: Sequence[str]) -> None:
    for locus in called:
        locus["candidates"] = sorted(candidates(locus, constructs))
        locus["construct"] = ""
        locus["confidence"] = "unresolved"
        locus["evidence"] = ""
        locus["attribution_note"] = ""


def _apply_unique(called: Sequence[Dict], constructs: Sequence[str],
                  used: Dict[str, int]) -> None:
    """(a) Sequence evidence. Not capped: nothing downstream overrules it."""
    for locus in called:
        n = unique_attribution(locus, constructs)
        if n:
            locus.update(construct=n, confidence="certain", evidence="unique",
                         attribution_note="vector side in sequence unique "
                                          "to %s" % n)
            used[n] = used.get(n, 0) + 1


def _apply_vector_end(called: Sequence[Dict], ends: Dict[str, List[int]],
                      end_tol: int, used: Dict[str, int]) -> None:
    """(b) Structural evidence: at exactly one construct's declared end.

    Also not capped -- it is what lets a junction inside SHARED sequence still
    name its construct, when the constructs were linearized at different sites.
    """
    for locus in called:
        if locus["construct"]:
            continue
        n, matched, tested = end_attribution(locus, ends, end_tol)
        if n:
            note = "within %d bp of %s's declared end" % (end_tol, n)
            untested = [c for c in locus["candidates"] if c not in tested]
            if untested:
                note += ("; %s could not be tested (no mapping or no "
                         "--vector-ends for it)" % ", ".join(untested))
            locus.update(construct=n, confidence="inferred",
                         evidence="vector-end", attribution_note=note)
            used[n] = used.get(n, 0) + 1
        elif len(matched) > 1:
            locus["attribution_note"] = (
                "at the declared end of more than one construct (%s)"
                % ", ".join(matched))


def _apply_ceiling(called: Sequence[Dict], copies: Dict[str, Optional[float]],
                   caps: Dict[str, Optional[int]], used: Dict[str, int],
                   margin: float) -> None:
    """(c) The hard cap, then (d) unresolved.

    Loci are taken best-supported first, so cap consumption is deterministic.
    Support decides only the ORDER in which interchangeable loci consume
    capacity; it never decides which construct a locus goes to.
    """
    rest = [locus for locus in called if not locus["construct"]]
    rest.sort(key=lambda locus: (-locus["total_support"],
                                 natural_key(locus["contig"]),
                                 locus["start"]))
    for locus in rest:
        n, why = ceiling_attribution(locus["candidates"], copies, caps,
                                     used, margin)
        if n:
            locus.update(construct=n, confidence="inferred",
                         evidence="copy-ceiling", attribution_note=why)
            used[n] = used.get(n, 0) + 1
        else:
            locus.update(construct="", confidence="unresolved",
                         evidence="", attribution_note=why)


def _ceiling_conflicts(constructs: Sequence[str], caps: Dict[str, Optional[int]],
                       used: Dict[str, int]) -> List[str]:
    """Sequence evidence that overfills a ceiling is a contradiction to report,
    not a locus to suppress."""
    out: List[str] = []
    for n in constructs:
        cap = caps.get(n)
        if cap is not None and used.get(n, 0) > cap:
            out.append(
                "%s carries %d attributed site(s) but its copy number only "
                "ceils at %d. Sequence-level evidence (unique sequence, a "
                "declared vector end) is not overruled by the ceiling, so the "
                "sites stand and the CONTRADICTION is reported: either the "
                "copy number is underestimated, the ploidy is wrong, or those "
                "sites carry fewer than --min-copies-per-site copies each."
                % (n, used[n], cap))
    return out


def attribute(called: List[Dict], constructs: Sequence[str],
              copies: Dict[str, Optional[float]], caps: Dict[str, Optional[int]],
              ends: Dict[str, List[int]], end_tol: int, margin: float
              ) -> Tuple[Dict[str, int], List[str]]:
    """Run the priority a -> b -> c -> d over every called locus, in place.

    Sequence evidence (a) and the declared vector ends (b) outrank the
    ceiling, so they go first and are not capped. The ceiling then strictly
    governs (c), the only step where a locus could otherwise be placed by
    argument -- which is why no support ratio is computed anywhere here.
    """
    used: Dict[str, int] = {n: 0 for n in constructs}
    _init_tiers(called, constructs)
    _apply_unique(called, constructs, used)
    _apply_vector_end(called, ends, end_tol, used)
    _apply_ceiling(called, copies, caps, used, margin)
    return used, _ceiling_conflicts(constructs, caps, used)


# --------------------------------------------------------------------------
# 4. site counts
# --------------------------------------------------------------------------

def site_counts(called: Sequence[Dict], constructs: Sequence[str],
                copies: Dict[str, Optional[float]],
                caps: Dict[str, Optional[int]]) -> Dict[str, Dict]:
    """One bounded site range per construct -- not per sample, not per mapping.

    The floor counts only loci attributed on sequence evidence or a declared
    vector end, the two tiers that do not depend on the copy number. The upper
    bound is the copy-number ceiling, which is what the measured DNA can
    source. Copy-ceiling attributions sit between the two: they are inside the
    range but are not part of the floor, because the argument for them is
    arithmetic, not sequence.
    """
    out: Dict[str, Dict] = {}
    for n in constructs:
        certain = [locus for locus in called
                   if locus["construct"] == n
                   and locus["evidence"] == "unique"]
        vend = [locus for locus in called
                if locus["construct"] == n
                and locus["evidence"] == "vector-end"]
        ceil_ = [locus for locus in called
                 if locus["construct"] == n
                 and locus["evidence"] == "copy-ceiling"]
        unres = [locus for locus in called
                 if not locus["construct"] and n in locus["candidates"]]
        floor = len(certain) + len(vend)
        cap = caps.get(n)
        hi = None if cap is None else max(cap, floor + len(ceil_))
        out[n] = {
            "copies_per_haploid": copies.get(n),
            "ceiling": cap,
            "certain": len(certain),
            "vector_end": len(vend),
            "copy_ceiling": len(ceil_),
            "attributed": floor + len(ceil_),
            "unresolved_candidate_for": len(unres),
            "site_floor": floor,
            "site_high": hi,
            "site_range": ("%d-%s" % (floor, hi if hi is not None else "?")
                           if hi is None or hi != floor else str(floor)),
        }
    return out


# --------------------------------------------------------------------------
# output
# --------------------------------------------------------------------------

def _vector_positions(locus: Dict) -> str:
    return ";".join("%s:%s" % (e["mapping"], e["vector_pos"])
                    for e in locus["entries"])


def _regions(locus: Dict) -> str:
    return ";".join("%s:%s" % (e["mapping"], e["region"] or "-")
                    for e in locus["entries"])


CALLED_COLS = ["contig", "start", "end", "flanks", "footprint", "construct",
               "confidence", "evidence", "candidates", "mappings",
               "vector_positions", "regions", "total_support", "n_samples"]


def write_table(path: str, loci: Sequence[Dict], samples: Sequence[str],
                extra: Sequence[str] = (), tail: Sequence[str] = ()) -> None:
    """`tail` columns (the host annotation) go after `note`, taken from the
    locus dict, so a table without them is unchanged."""
    cols = list(CALLED_COLS) + ["support_" + s for s in samples] \
        + list(extra) + ["note"] + list(tail)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\t".join(cols) + "\n")
        for locus in loci:
            row: Dict[str, object] = {
                "contig": locus["contig"], "start": locus["start"],
                "end": locus["end"],
                "flanks": locus["flanks"],
                "footprint": ("" if locus["footprint"] is None
                              else locus["footprint"]),
                "construct": locus.get("construct", ""),
                "confidence": locus.get("confidence", ""),
                "evidence": locus.get("evidence", ""),
                "candidates": ",".join(locus.get("candidates", [])),
                "mappings": ",".join(locus["mappings"]),
                "vector_positions": _vector_positions(locus),
                "regions": _regions(locus),
                "total_support": locus["total_support"],
                "n_samples": locus["n_samples"],
                "note": locus.get("attribution_note", ""),
                "reason": locus.get("subthreshold_reason", ""),
            }
            for s in samples:
                row["support_" + s] = locus["support"].get(s, 0)
            for c in tail:
                row[c] = locus.get(c, "")
            if locus.get("same_mapping_merge"):
                row["note"] = (str(row["note"]) + " | two entries of one "
                               "mapping were linked through another mapping's "
                               "entry: check this group").strip(" |")
            fh.write("\t".join(str(row.get(c, "")) for c in cols) + "\n")


def _md_header(result: Dict) -> List[str]:
    p = result["parameters"]
    return ["# Called integration sites", "",
            "Produced by `reconcile_sites.py`. Locus collapse, retention and "
            "attribution were done here, deterministically. **Report this "
            "table; do not re-derive it.**", "",
            "Mappings reconciled: %s. Samples: %s."
            % (", ".join(result["mappings"]), ", ".join(result["samples"])), "",
            "Parameters: window %d bp, --min-reads %d, --min-samples %d, "
            "--end-tol %d bp, ploidy %s, --min-copies-per-site %s, "
            "--cn-margin %s. %s."
            % (p["window"], p["min_reads"], p["min_samples"], p["end_tol"],
               p["ploidy"], p["min_copies_per_site"], p["cn_margin"],
               result["copy_number_source"]), ""]


def _md_per_construct(result: Dict) -> List[str]:
    lines = ["## Sites per construct", "",
         "| Construct | copies/haploid | ceiling | certain | vector-end | "
         "copy-ceiling | site count | unresolved candidate for |",
         "|---|---|---|---|---|---|---|---|"]
    for n, d in result["per_construct"].items():
        lines.append("| %s | %s | %s | %d | %d | %d | **%s** | %d |"
                 % (n,
                    "%.2f" % d["copies_per_haploid"]
                    if d["copies_per_haploid"] is not None else "not reported",
                    d["ceiling"] if d["ceiling"] is not None else "unknown",
                    d["certain"], d["vector_end"], d["copy_ceiling"],
                    d["site_range"], d["unresolved_candidate_for"]))
    return lines + ["",
                "The floor is certain + vector-end loci: the tiers that do "
                "not depend on the copy number. The upper bound is the "
                "copy-number ceiling, ceil(copies x ploidy / "
                "min_copies_per_site), which is what the measured DNA can "
                "source. A construct can never hold more copy-ceiling sites "
                "than its ceiling allows, whatever the read support ratios "
                "suggest.", ""]


def _md_accounting(result: Dict) -> List[str]:
    rows = result.get("site_accounting") or []
    if not rows:
        return []
    lines = ["## Do the called sites account for the DNA?", "",
         "| Construct | copies/haploid | sites the DNA implies | attributed | "
         "unaccounted for | |", "|---|---|---|---|---|---|"]
    for r in rows:
        lines.append("| %s | %s | %s | %d | %s | %s |"
                 % (r["construct"],
                    "%.2f" % r["copies_per_haploid"]
                    if r["copies_per_haploid"] is not None else "not reported",
                    r["implied_sites"] if r["implied_sites"] is not None
                    else "unknown",
                    r["attributed"],
                    r["unaccounted"] if r["unaccounted"] is not None else "-",
                    r["verdict"]))
    return lines + ["", "The sites the DNA implies is the copy-number ceiling, "
                "ceil(copies x ploidy / min_copies_per_site). A shortfall "
                "means missing sites, an extrachromosomal pool, or more "
                "copies per site than the divisor assumes.", ""]


def _md_collapse(result: Dict) -> List[str]:
    p = result["parameters"]
    return ["## Collapse and retention", "",
            "| Step | loci |", "|---|---|",
            "| entries read from all mapping tables | %d |" % result["n_entries"],
            "| unified loci after cross-mapping collapse | %d |"
            % result["n_unified"],
            "| called (total support >= %d, in >= %d samples) | %d |"
            % (p["min_reads"], p["min_samples"], result["n_called"]),
            "| sub-threshold (single-read / noise tail) | %d |" % result["n_sub"],
            "",
            "%d duplicate report(s) of one junction were collapsed across "
            "mappings. A site seen in two mapping runs is one event; listing "
            "it under both constructs would double the site count."
            % (result["n_entries"] - result["n_unified"]), ""]


def _md_unresolved(result: Dict) -> List[str]:
    unres = result["unresolved"]
    if not unres:
        return ["## Unresolved loci", "",
                "None: every called locus was attributed.", ""]
    lines = ["## Unresolved loci", "",
         "Listed once each, attributed to no construct. These are the "
         "validation priorities (`references/validation.md`).", "",
         "| Locus | candidates | support | samples | why |",
         "|---|---|---|---|---|"]
    for u in unres:
        lines.append("| %s:%d-%d | %s | %d | %d | %s |"
                 % (u["contig"], u["start"], u["end"],
                    ", ".join(u["candidates"]) or "none",
                    u["total_support"], u["n_samples"], u["note"]))
    return lines + [""]


def render_md(result: Dict) -> str:
    lines = _md_header(result)
    if result["caveats"]:
        lines += ["## Read this before using any count below", ""]
        lines += ["> %s" % c for c in result["caveats"]] + [""]
    lines += _md_per_construct(result)
    lines += _md_accounting(result)
    lines += _md_collapse(result)
    lines += _md_unresolved(result)
    return "\n".join(lines)


# --------------------------------------------------------------------------

def read_all_loci(tables: Dict[str, str]
                  ) -> Tuple[List[Dict], List[str], Dict[str, List[str]]]:
    """Every mapping's loci table, as one flat list of entries."""
    entries: List[Dict] = []
    samples: List[str] = []
    per_table: Dict[str, List[str]] = {}
    for mapping, path in tables.items():
        rows, ss = read_loci(path, mapping)
        entries += rows
        per_table[mapping] = ss
        for s in ss:
            if s not in samples:
                samples.append(s)
        print("  %s: %d loci, sample(s) %s"
              % (mapping, len(rows), ", ".join(ss)), file=sys.stderr)
    if not entries:
        raise ValueError("no loci in any table; nothing to reconcile")
    return entries, samples, per_table


def preflight(args: argparse.Namespace, tables: Dict[str, str],
              samples: Sequence[str], per_table: Dict[str, List[str]],
              copies: Dict[str, Optional[float]], assumed_ploidy: bool
              ) -> List[Tuple[str, str, str]]:
    """Everything the inputs cannot deliver, said before any count is shown.

    Each of these silently weakens a site count if it is not stated: an
    assumed ploidy scales every ceiling, a missing copy number removes a
    construct from the ceiling tier entirely, and an impossible
    `--min-samples` returns an empty site list that looks like a result.

    Returns (check id, status, message). Every message is a caveat in the
    report; the two that leave no usable site result are failing checks.
    """
    out: List[Tuple[str, str, str]] = []
    if assumed_ploidy:
        out.append(("ploidy", WARN, 
            "PLOIDY ASSUMED to be 2. It multiplies every copy-number ceiling, "
            "so it scales the upper bound of every site count directly. Many "
            "immortalised lines are aneuploid and drift on passage: measure "
            "it and pass --ploidy."))
    odd = {m: ss for m, ss in per_table.items() if list(ss) != list(samples)}
    if odd:
        out.append(("samples_consistent", WARN,
            "The mapping tables do not carry the same samples (%s). Support "
            "is combined over the union, so a locus can look sample-specific "
            "only because one mapping run never saw that sample."
            % "; ".join("%s: %s" % (m, ", ".join(ss) or "none")
                        for m, ss in odd.items())))
    if args.min_samples > len(samples):
        out.append(("min_samples", FAIL,
            "--min-samples is %d but only %d sample(s) were supplied, so NO "
            "locus can be called. Lower it deliberately (and say so), or add "
            "the other samples' loci tables."
            % (args.min_samples, len(samples))))
    ends = parse_ends(args.vector_ends)
    missing_ends = [m for m in tables if m not in ends or not ends[m]]
    if missing_ends:
        out.append(("vector_ends", WARN,
            "No --vector-ends for %s, so the vector-end tier could not be "
            "tested there. Loci that sit at those constructs' declared ends "
            "fall through to the copy-number ceiling or to unresolved."
            % ", ".join(sorted(missing_ends))))
    no_cn = [n for n in tables if copies.get(n) is None]
    if no_cn:
        out.append(("copy_number", FAIL,
            "No copy number for %s, so it has no ceiling and the "
            "copy-ceiling tier will not attribute anything to it. Supply "
            "copy_number.py --json output, or --copies %s=N."
            % (", ".join(sorted(no_cn)), sorted(no_cn)[0])))
    return out


#: preflight checks, in a fixed order, with what a pass means
PREFLIGHT_CHECKS = (
    ("min_samples", "--min-samples can be met by the samples supplied"),
    ("copy_number", "every mapped construct has a copy number, so a ceiling"),
    ("ploidy", "ploidy was supplied, not assumed"),
    ("vector_ends", "declared ends supplied for every mapped construct"),
    ("samples_consistent", "every mapping table carries the same samples"),
)


def ceiling_divisor_finding(args: argparse.Namespace,
                            measured: Optional[float],
                            modality: Optional[str]) -> Optional[str]:
    """Is the ceiling nominal where the run's own measurement could bind it?

    `--min-copies-per-site` is never taken from the measurement silently: it
    changes every site range, so it stays a decision. But leaving it at 1 on
    a line whose copy-number stage measured arrays is worth saying out loud.
    """
    capable = (MODALITIES.get(modality or "auto") or {}).get("concatemers")
    if (measured is None or measured <= 1.0 or capable is False
            or args.min_copies_per_site != 1.0):
        return None
    return ("--min-copies-per-site was left at 1, but the copy-number stage "
            "measured %.4g copies per site on this line. With a divisor of 1 "
            "every ceiling is nominal: the upper bound of each site range is "
            "the copy number itself. To make the ceiling binding, rerun with "
            "--min-copies-per-site %.4g (references/attribution.md, section "
            "3c); to keep 1 deliberately, say so in the report."
            % (measured, measured))


#: `references/quantification.md` section 7: agreement within ~25% means the
#: called sites account for the DNA, so a shortfall needs fewer than this
#: fraction of the implied sites attributed...
SHORTFALL_FRACTION = 0.75
#: ...and at least one whole site's worth of DNA unaccounted for. Measured
#: against the implied count before it is rounded up to the ceiling, so a
#: copy number just over an integer does not manufacture a missing site.
SHORTFALL_SITES = 1.0


def site_accounting(per_construct: Dict[str, Dict], unresolved: int,
                    sites_per_copy: float) -> Tuple[List[Dict], Optional[str]]:
    """Do the called sites account for the DNA each construct was measured at?

    The ceiling IS the implied site count -- ceil(copies x ploidy /
    min_copies_per_site) -- so the comparison the report needs is the called
    sites against the ceiling, per construct. Done here because it is the
    step most easily skipped by hand, and because doing it for the primary
    construct alone misses a shortfall on any other one.

    `sites_per_copy` is ploidy / min_copies_per_site, the factor the ceiling
    was built with. The shortfall test uses copies x sites_per_copy unrounded:
    the ceiling rounds up, and a rounded-up site is not missing DNA.
    """
    rows: List[Dict] = []
    short: List[str] = []
    for n, d in per_construct.items():
        cap, got = d["ceiling"], d["attributed"]
        implied = (None if cap is None else
                   max(0.0, d["copies_per_haploid"]) * sites_per_copy)
        row = {"construct": n, "copies_per_haploid": d["copies_per_haploid"],
               "implied_sites": cap, "attributed": got,
               "unaccounted": None if cap is None else cap - got,
               "verdict": ""}
        if cap is None:
            row["verdict"] = ("no copy number, so the DNA cannot be "
                              "accounted for")
        elif got > cap:
            row["verdict"] = ("more sites than the copy number implies; see "
                              "the ceiling contradiction above")
        elif (implied - got >= SHORTFALL_SITES
              and got < SHORTFALL_FRACTION * implied):
            row["verdict"] = ("SHORTFALL: %d of the %d sites the measured DNA "
                              "implies are not accounted for"
                              % (cap - got, cap))
            short.append("%s: %d of %d" % (n, cap - got, cap))
        else:
            row["verdict"] = "the called sites account for the measured DNA"
        rows.append(row)
    if not short:
        return rows, None
    return rows, (
        "SITE SHORTFALL (%s). The measured copy number implies more sites "
        "than were attributed%s. Either sites are missing -- flanks lost to "
        "repeats are the usual reason -- or part of the vector is "
        "extrachromosomal, or each site carries more copies than "
        "--min-copies-per-site assumes, which is the first thing to check on "
        "a concatemer line. `references/quantification.md` section 7 and "
        "`references/validation.md` section 4 give the tests."
        % ("; ".join(short),
           ", and %d called locus/loci are unresolved" % unresolved
           if unresolved else ""))


def stage_checks(checks: Checks, found: Sequence[Tuple[str, str, str]],
                 attr_caveats: Sequence[str], divisor: Optional[str],
                 n_called: int, shortfall: Optional[str] = None) -> Checks:
    by_id = {cid: (status, msg) for cid, status, msg in found}
    for cid, ok in PREFLIGHT_CHECKS:
        status, msg = by_id.get(cid, (PASS, ok))
        checks.add(cid, status, msg)
    if attr_caveats:
        checks.add("ceiling_contradiction", FAIL, " ".join(attr_caveats))
    else:
        checks.add("ceiling_contradiction", PASS,
                   "no construct holds more sites than its ceiling admits")
    checks.add("ceiling_divisor", WARN if divisor else PASS,
               divisor or "ceiling divisor needs no attention")
    checks.add("sites_called", PASS if n_called else WARN,
               "%d locus/loci called" % n_called)
    checks.add("site_accounting", WARN if shortfall else PASS,
               shortfall or "the called sites account for every construct's "
                            "measured DNA")
    return checks


def build_result(args: argparse.Namespace, tables: Dict[str, str],
                 samples: Sequence[str], cn_how: str, ploidy: float,
                 ploidy_source: str, entries: Sequence[Dict],
                 unified: Sequence[Dict], called: Sequence[Dict],
                 sub: Sequence[Dict], caps: Dict[str, Optional[int]],
                 used: Dict[str, int], per_construct: Dict[str, Dict],
                 caveats: Sequence[str]) -> Dict:
    """The whole result, including every parameter it was produced with.

    The parameters are part of the output because a site count is only
    reproducible alongside the thresholds that produced it.
    """
    return {
        "mappings": sorted(tables),
        "samples": list(samples),
        "copy_number_source": cn_how,
        "parameters": {"window": args.window, "min_reads": args.min_reads,
                       "min_samples": args.min_samples,
                       "end_tol": args.end_tol,
                       "ploidy": ploidy,
                       "min_copies_per_site": args.min_copies_per_site,
                       "cn_margin": args.cn_margin,
                       "ploidy_source": ploidy_source},
        "n_entries": len(entries), "n_unified": len(unified),
        "n_called": len(called), "n_sub": len(sub),
        "ceilings": caps,
        "attributed": used,
        "per_construct": per_construct,
        "unresolved": [{"contig": locus["contig"], "start": locus["start"],
                        "end": locus["end"], "candidates": locus["candidates"],
                        "total_support": locus["total_support"],
                        "n_samples": locus["n_samples"],
                        "note": locus["attribution_note"]}
                       for locus in called if not locus["construct"]],
        "caveats": list(caveats),
    }


def write_outputs(args: argparse.Namespace, result: Dict,
                  called: Sequence[Dict],
                  sub: Sequence[Dict], samples: Sequence[str]) -> None:
    if args.out_called:
        tail: Sequence[str] = ()
        if args.genes:
            from host_reference import ANNOT_COLS, annotate_rows
            called = annotate_rows(called, args.genes)
            tail = ANNOT_COLS
            regions: Dict[str, int] = {}
            for locus in called:
                reg = locus["host_region"]
                regions[reg] = regions.get(reg, 0) + 1
            result["host_annotation"] = {"genes": args.genes,
                                         "regions": regions}
        write_table(args.out_called, called, samples, tail=tail)
    if args.out_subthreshold:
        write_table(args.out_subthreshold, sub, samples, extra=["reason"])
    md = render_md(result)
    if args.md:
        with open(args.md, "w", encoding="utf-8") as fh:
            fh.write(md)
    if args.summary:
        with open(args.summary, "w", encoding="utf-8") as fh:
            json.dump(result, fh, indent=2)
    if not args.md and not args.summary:
        print(md)


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--identity", default="",
                    help="identity map (required unless --run supplies it)")
    ap.add_argument("--loci", action="append", default=[],
                    metavar="MAPPING=PATH",
                    help="loci table from classify_junctions.py --out-loci, "
                         "one per mapping run; repeatable (required unless "
                         "--run supplies them)")
    ap.add_argument("--copy-number", default="",
                    help="copy_number.py --json output, for each construct's "
                         "copies per haploid")
    ap.add_argument("--cn-sample", default="",
                    help="take copy numbers from this sample only (default: "
                         "the per-construct maximum over all samples)")
    ap.add_argument("--copies", action="append", default=[],
                    metavar="CONSTRUCT=N",
                    help="copies per haploid for a construct, overriding the "
                         "copy-number JSON; repeatable")
    ap.add_argument("--vector-ends", action="append", default=[],
                    metavar="CONSTRUCT=a,b",
                    help="the delivered molecule's declared ends in that "
                         "construct's own numbering; repeatable, one per "
                         "construct")
    ap.add_argument("--ploidy", type=float, default=None,
                    help="haploid equivalents per cell. Default: taken from "
                         "the copy-number JSON if it implies one, else 2 with "
                         "a caveat. It scales every ceiling, so measure it.")
    ap.add_argument("--min-copies-per-site", type=float, default=1.0,
                    help="fewest copies a real site is assumed to carry, used "
                         "as the ceiling's divisor (default 1.0, the most "
                         "permissive honest value and the only one valid for "
                         "integrase or transposase delivery)")
    ap.add_argument("--cn-margin", type=float, default=2.0,
                    help="how far the highest-copy admissible candidate must "
                         "lead the next before a shared locus is attributed "
                         "to it (default 2.0 = twice the copy number). Below "
                         "this the locus is unresolved.")
    ap.add_argument("--window", type=int, default=500,
                    help="host distance within which entries from DIFFERENT "
                         "mappings are the same junction (default 500, the "
                         "classify_junctions --pair-window)")
    ap.add_argument("--min-reads", type=int, default=5,
                    help="a called locus needs at least this much total "
                         "support (default 5)")
    ap.add_argument("--min-samples", type=int, default=2,
                    help="a called locus must appear in at least this many "
                         "independently prepared samples (default 2)")
    ap.add_argument("--end-tol", type=int, default=3)
    ap.add_argument("--out-called")
    ap.add_argument("--out-subthreshold")
    ap.add_argument("--summary")
    ap.add_argument("--md")
    ap.add_argument("--genes", default="",
                    help="host gene table from host_reference.py: adds the "
                         "host_region / host_gene / nearest_gene columns to "
                         "the called table (default: the host stage's, if "
                         "it recorded one)")
    add_check_arguments(ap)
    return ap


def merge_manifest(args: argparse.Namespace) -> Tuple[Optional[object], Dict]:
    """Fill unset arguments from the run manifest and name this stage's
    outputs in it. Explicit flags win. Returns the run and what the earlier
    stages recorded that is not an argument (measured copies per site, the
    modality, and where the ploidy came from)."""
    given: Dict = {"copies_per_site": None, "modality": None,
                   "ploidy_source": "--ploidy" if args.ploidy is not None
                                    else ""}
    try:
        run = load_run(args.run)
    except RunError as exc:
        raise InputError(str(exc))
    if run is None:
        return None, given
    try:
        recorded = reconcile_inputs(run)
        args.identity = args.identity or recorded["identity"] or ""
        if not args.genes and run.stage("host").get("genes"):
            args.genes = run.resolve(run.stage("host")["genes"])
        args.loci = args.loci or ["%s=%s" % kv
                                  for kv in recorded["loci"].items()]
        args.vector_ends = args.vector_ends or [
            "%s=%s" % (m, ",".join(str(e) for e in ends))
            for m, ends in recorded["vector_ends"].items()]
        # --copies overrides the JSON only for the constructs it names, so
        # the manifest's copy-number JSON is still read underneath it.
        args.copy_number = args.copy_number or recorded["copy_number"] or ""
        if args.ploidy is None and recorded["ploidy"] is not None:
            args.ploidy = float(recorded["ploidy"])
            given["ploidy_source"] = "run manifest"
        given["copies_per_site"] = recorded["copies_per_site"]
        given["modality"] = recorded["modality"]
        args.out_called = run.output(args.out_called, "called_sites.tsv")
        args.out_subthreshold = run.output(args.out_subthreshold,
                                           "subthreshold_sites.tsv")
        args.summary = run.output(args.summary, "site_summary.json")
        args.md = run.output(args.md, "site_summary.md")
    except RunError as exc:
        raise InputError(str(exc))
    return run, given


def read_inputs(args: argparse.Namespace) -> Dict:
    """Every input this stage reads from disk, validated against each other.

    Raises `InputError` for anything that makes the stage unanswerable: a
    missing identity map or loci table, a construct not in the map, an
    unreadable copy-number JSON, a ploidy that is not positive.
    """
    if not args.identity or not args.loci:
        raise InputError("need --identity and at least one --loci "
                         "MAPPING=PATH, or --run with a manifest that Stages "
                         "1 and 2 have recorded them in")
    if args.genes and not os.path.exists(args.genes):
        raise InputError("host gene table %s not found: rerun "
                         "host_reference.py or drop --genes" % args.genes)
    try:
        ident = load_json(args.identity)
        tables = parse_kv(args.loci, "--loci")
        ends = parse_ends(args.vector_ends)
        explicit = parse_kv(args.copies, "--copies")
    except (OSError, ValueError) as exc:
        raise InputError(str(exc))

    known = list(ident.get("construct_order") or ident["constructs"])
    for name in list(tables) + list(ends) + list(explicit):
        if name not in known:
            raise InputError(
                "%r is not in the identity map (%s). Every construct in the "
                "line must be in the map, including ones the provider did "
                "not map against." % (name, ", ".join(known)))
    try:
        entries, samples, per_table = read_all_loci(tables)
        copies, implied_ploidy, cn_how = load_copies(
            args.copy_number or None, explicit, args.cn_sample or None)
    except (OSError, ValueError) as exc:
        raise InputError(str(exc))
    return {"ident": ident, "known": known, "tables": tables, "ends": ends,
            "entries": entries, "samples": samples, "per_table": per_table,
            "copies": copies, "implied_ploidy": implied_ploidy,
            "cn_how": cn_how}


def resolve_ploidy(args: argparse.Namespace, implied: Optional[float],
                   source: str) -> Tuple[float, str, bool]:
    """The ploidy to use, where it came from, and whether it was assumed."""
    ploidy = args.ploidy if args.ploidy is not None else implied
    if ploidy is not None and ploidy <= 0:
        raise InputError("ploidy must be greater than 0, got %g (from %s)"
                         % (ploidy, source or "the copy-number JSON"))
    if ploidy is None:
        return 2.0, "assumed 2", True
    return ploidy, source or "copy-number JSON", False


def compute_sites(args: argparse.Namespace, inputs: Dict, ploidy: float,
                  given: Dict) -> Dict:
    """Collapse, retain, attribute, count -- the stage's whole computation.

    Reads nothing and writes nothing, so the site result is testable on its
    inputs alone.
    """
    known, tables = inputs["known"], inputs["tables"]
    copies, samples = inputs["copies"], inputs["samples"]
    found = preflight(args, tables, samples, inputs["per_table"], copies,
                      given["ploidy_assumed"])
    caveats = [msg for _, _, msg in found]

    unified = collapse(inputs["entries"], args.window, samples)
    called, sub = retain(unified, args.min_reads, args.min_samples)
    caps = ceilings(copies, known, ploidy, args.min_copies_per_site)
    used, attr_caveats = attribute(called, known, copies, caps,
                                   inputs["ends"], args.end_tol,
                                   args.cn_margin)
    caveats += attr_caveats
    divisor = ceiling_divisor_finding(args, given["copies_per_site"],
                                      given["modality"])
    if divisor:
        caveats.append(divisor)
    constructs = [n for n in known
                  if n in tables or copies.get(n) is not None]
    per_construct = site_counts(called, constructs, copies, caps)
    n_unresolved = sum(1 for locus in called if not locus["construct"])
    accounting, shortfall = site_accounting(
        per_construct, n_unresolved, ploidy / args.min_copies_per_site)
    if shortfall:
        caveats.append(shortfall)
    return {"found": found, "unified": unified, "called": called, "sub": sub,
            "caps": caps, "used": used, "per_construct": per_construct,
            "attr_caveats": attr_caveats, "divisor": divisor,
            "shortfall": shortfall, "accounting": accounting,
            "caveats": caveats}


def report_stderr(result: Dict, per_construct: Dict[str, Dict],
                  n_entries: int, n_unified: int, n_called: int,
                  n_sub: int) -> int:
    """The stage's one-screen summary; returns the unresolved count."""
    print("collapsed %d entries -> %d unified loci; %d called, %d "
          "sub-threshold" % (n_entries, n_unified, n_called, n_sub),
          file=sys.stderr)
    for n, d in per_construct.items():
        print("    %-12s sites %-8s (certain %d, vector-end %d, copy-ceiling "
              "%d; ceiling %s)"
              % (n, d["site_range"], d["certain"], d["vector_end"],
                 d["copy_ceiling"],
                 d["ceiling"] if d["ceiling"] is not None else "unknown"),
              file=sys.stderr)
    n_unresolved = len(result["unresolved"])
    print("    unresolved   %d (listed once each, attributed to no construct)"
          % n_unresolved, file=sys.stderr)
    return n_unresolved


def run_stage(args: argparse.Namespace) -> int:
    if args.min_copies_per_site <= 0:
        raise InputError("--min-copies-per-site must be greater than 0, got "
                         "%g" % args.min_copies_per_site)
    run, given = merge_manifest(args)
    inputs = read_inputs(args)
    ploidy, ploidy_source, assumed = resolve_ploidy(
        args, inputs["implied_ploidy"], given["ploidy_source"])
    given["ploidy_assumed"] = assumed

    computed = compute_sites(args, inputs, ploidy, given)
    checks = stage_checks(Checks("reconcile_sites", args.accept),
                          computed["found"], computed["attr_caveats"],
                          computed["divisor"], len(computed["called"]),
                          computed["shortfall"])

    result = build_result(args, inputs["tables"], inputs["samples"],
                          inputs["cn_how"], ploidy, ploidy_source,
                          inputs["entries"], computed["unified"],
                          computed["called"], computed["sub"],
                          computed["caps"], computed["used"],
                          computed["per_construct"], computed["caveats"])
    result["site_accounting"] = computed["accounting"]
    result["checks"] = checks.as_list()
    write_outputs(args, result, computed["called"], computed["sub"],
                  inputs["samples"])

    n_unresolved = report_stderr(result, computed["per_construct"],
                                 len(inputs["entries"]),
                                 len(computed["unified"]),
                                 len(computed["called"]),
                                 len(computed["sub"]))
    if run is not None:
        run.record("reconcile", {
            "called": run.relative(args.out_called),
            "subthreshold": run.relative(args.out_subthreshold),
            "summary": run.relative(args.summary),
            "md": run.relative(args.md),
            "sites": {n: d["site_range"]
                      for n, d in computed["per_construct"].items()},
            "unresolved": n_unresolved, "annotated": bool(args.genes),
            "checks": checks.as_list()})
    return checks.finish(args.strict)


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    try:
        return run_stage(args)
    except InputError as exc:
        print("error: %s" % exc, file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
