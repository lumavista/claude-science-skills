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
from via_common import load_json, natural_key  # noqa: E402


# --------------------------------------------------------------------------
# input
# --------------------------------------------------------------------------

def _int(v: Optional[str]) -> Optional[int]:
    try:
        return int(float(v))                          # type: ignore[arg-type]
    except (TypeError, ValueError):
        return None


def read_loci(path: str, mapping: str) -> Tuple[List[Dict], List[str]]:
    """Read one loci table written by classify_junctions.py --out-loci.

    Returns its loci as per-mapping *entries*, plus the sample names taken
    from the `support_<sample>` columns.
    """
    with open(path, encoding="utf-8") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        need = {"contig", "start", "end", "vector_pos", "membership"}
        missing = need - set(hdr)
        if missing:
            raise ValueError(
                "%s: not a classify_junctions loci table (no %s column). Pass "
                "the file written by --out-loci."
                % (path, ", ".join(sorted(missing))))
        samples = [h[len("support_"):] for h in hdr if h.startswith("support_")]
        if not samples:
            raise ValueError(
                "%s: no support_<sample> column, so retention by sample count "
                "is impossible. Re-run classify_junctions.py with --out-loci."
                % path)
        entries: List[Dict] = []
        for n, line in enumerate(fh, start=2):
            f = line.rstrip("\n").split("\t")
            if len(f) != len(hdr):
                continue
            d = dict(zip(hdr, f))
            start, end = _int(d["start"]), _int(d["end"])
            if start is None or end is None:
                continue
            entries.append({
                "mapping": mapping, "row": n, "contig": d["contig"],
                "start": min(start, end), "end": max(start, end),
                "flanks": _int(d.get("flanks")) or 1,
                "footprint": _int(d.get("footprint")),
                "vector_pos": _int(d.get("vector_pos")),
                "region": d.get("region", ""),
                "membership": frozenset(
                    m.strip() for m in d.get("membership", "").split("+")
                    if m.strip()),
                "support": {s: (_int(d.get("support_" + s)) or 0)
                            for s in samples},
            })
    return entries, samples


def parse_kv(items: Sequence[str], flag: str) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for item in items:
        k, sep, v = item.partition("=")
        if not sep or not k.strip():
            raise ValueError("%s needs NAME=VALUE, got %r" % (flag, item))
        out[k.strip()] = v.strip()
    return out


def parse_ends(items: Sequence[str]) -> Dict[str, List[int]]:
    """`--vector-ends CONSTRUCT=a,b`, one per construct."""
    out: Dict[str, List[int]] = {}
    for name, spec in parse_kv(items, "--vector-ends").items():
        vals: List[int] = []
        for x in spec.split(","):
            x = x.strip()
            if not x:
                continue
            if not x.lstrip("-").isdigit():
                raise ValueError("--vector-ends takes integers, got %r for %s"
                                 % (x, name))
            vals.append(int(x))
        out[name] = vals
    return out


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


class _Groups:
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


def _link_across_mappings(entries: Sequence[Dict], window: int) -> _Groups:
    """Link entries of DIFFERENT mappings whose host intervals are within
    `window`. Entries of the same mapping are never linked to each other:
    classify_junctions.py has already decided whether they are one locus or
    two, and re-merging them here would deflate the site count."""
    g = _Groups(len(entries))
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
    loci.sort(key=lambda L: (natural_key(L["contig"]), L["start"]))
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
    for L in loci:
        why: List[str] = []
        if L["total_support"] < min_reads:
            why.append("support %d < --min-reads %d"
                       % (L["total_support"], min_reads))
        if L["n_samples"] < min_samples:
            why.append("in %d sample(s) < --min-samples %d"
                       % (L["n_samples"], min_samples))
        if why:
            L["subthreshold_reason"] = "; ".join(why)
            sub.append(L)
        else:
            called.append(L)
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


def unique_attribution(locus: Dict, constructs: Sequence[str]) -> Optional[str]:
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
    for L in called:
        L["candidates"] = sorted(candidates(L, constructs))
        L["construct"] = ""
        L["confidence"] = "unresolved"
        L["evidence"] = ""
        L["attribution_note"] = ""


def _apply_unique(called: Sequence[Dict], constructs: Sequence[str],
                  used: Dict[str, int]) -> None:
    """(a) Sequence evidence. Not capped: nothing downstream overrules it."""
    for L in called:
        n = unique_attribution(L, constructs)
        if n:
            L.update(construct=n, confidence="certain", evidence="unique",
                     attribution_note="vector side in sequence unique to %s" % n)
            used[n] = used.get(n, 0) + 1


def _apply_vector_end(called: Sequence[Dict], ends: Dict[str, List[int]],
                      end_tol: int, used: Dict[str, int]) -> None:
    """(b) Structural evidence: at exactly one construct's declared end.

    Also not capped -- it is what lets a junction inside SHARED sequence still
    name its construct, when the constructs were linearized at different sites.
    """
    for L in called:
        if L["construct"]:
            continue
        n, matched, tested = end_attribution(L, ends, end_tol)
        if n:
            note = "within %d bp of %s's declared end" % (end_tol, n)
            untested = [c for c in L["candidates"] if c not in tested]
            if untested:
                note += ("; %s could not be tested (no mapping or no "
                         "--vector-ends for it)" % ", ".join(untested))
            L.update(construct=n, confidence="inferred", evidence="vector-end",
                     attribution_note=note)
            used[n] = used.get(n, 0) + 1
        elif len(matched) > 1:
            L["attribution_note"] = ("at the declared end of more than one "
                                     "construct (%s)" % ", ".join(matched))


def _apply_ceiling(called: Sequence[Dict], copies: Dict[str, Optional[float]],
                   caps: Dict[str, Optional[int]], used: Dict[str, int],
                   margin: float) -> None:
    """(c) The hard cap, then (d) unresolved.

    Loci are taken best-supported first, so cap consumption is deterministic.
    Support decides only the ORDER in which interchangeable loci consume
    capacity; it never decides which construct a locus goes to.
    """
    rest = [L for L in called if not L["construct"]]
    rest.sort(key=lambda L: (-L["total_support"], natural_key(L["contig"]),
                             L["start"]))
    for L in rest:
        n, why = ceiling_attribution(L["candidates"], copies, caps, used, margin)
        if n:
            L.update(construct=n, confidence="inferred", evidence="copy-ceiling",
                     attribution_note=why)
            used[n] = used.get(n, 0) + 1
        else:
            L.update(construct="", confidence="unresolved", evidence="",
                     attribution_note=why)


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
        certain = [L for L in called
                   if L["construct"] == n and L["evidence"] == "unique"]
        vend = [L for L in called
                if L["construct"] == n and L["evidence"] == "vector-end"]
        ceil_ = [L for L in called
                 if L["construct"] == n and L["evidence"] == "copy-ceiling"]
        unres = [L for L in called
                 if not L["construct"] and n in L["candidates"]]
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

def _vector_positions(L: Dict) -> str:
    return ";".join("%s:%s" % (e["mapping"], e["vector_pos"])
                    for e in L["entries"])


def _regions(L: Dict) -> str:
    return ";".join("%s:%s" % (e["mapping"], e["region"] or "-")
                    for e in L["entries"])


CALLED_COLS = ["contig", "start", "end", "flanks", "footprint", "construct",
               "confidence", "evidence", "candidates", "mappings",
               "vector_positions", "regions", "total_support", "n_samples"]


def write_table(path: str, loci: Sequence[Dict], samples: Sequence[str],
                extra: Sequence[str] = ()) -> None:
    cols = list(CALLED_COLS) + ["support_" + s for s in samples] \
        + list(extra) + ["note"]
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\t".join(cols) + "\n")
        for L in loci:
            row: Dict[str, object] = {
                "contig": L["contig"], "start": L["start"], "end": L["end"],
                "flanks": L["flanks"],
                "footprint": "" if L["footprint"] is None else L["footprint"],
                "construct": L.get("construct", ""),
                "confidence": L.get("confidence", ""),
                "evidence": L.get("evidence", ""),
                "candidates": ",".join(L.get("candidates", [])),
                "mappings": ",".join(L["mappings"]),
                "vector_positions": _vector_positions(L),
                "regions": _regions(L),
                "total_support": L["total_support"],
                "n_samples": L["n_samples"],
                "note": L.get("attribution_note", ""),
                "reason": L.get("subthreshold_reason", ""),
            }
            for s in samples:
                row["support_" + s] = L["support"].get(s, 0)
            if L.get("same_mapping_merge"):
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
    L = ["## Sites per construct", "",
         "| Construct | copies/haploid | ceiling | certain | vector-end | "
         "copy-ceiling | site count | unresolved candidate for |",
         "|---|---|---|---|---|---|---|---|"]
    for n, d in result["per_construct"].items():
        L.append("| %s | %s | %s | %d | %d | %d | **%s** | %d |"
                 % (n,
                    "%.2f" % d["copies_per_haploid"]
                    if d["copies_per_haploid"] is not None else "not reported",
                    d["ceiling"] if d["ceiling"] is not None else "unknown",
                    d["certain"], d["vector_end"], d["copy_ceiling"],
                    d["site_range"], d["unresolved_candidate_for"]))
    return L + ["",
                "The floor is certain + vector-end loci: the tiers that do "
                "not depend on the copy number. The upper bound is the "
                "copy-number ceiling, ceil(copies x ploidy / "
                "min_copies_per_site), which is what the measured DNA can "
                "source. A construct can never hold more copy-ceiling sites "
                "than its ceiling allows, whatever the read support ratios "
                "suggest.", ""]


def _md_collapse(result: Dict) -> List[str]:
    p = result["parameters"]
    return ["## Collapse and retention", "",
            "| Step | loci |", "|---|---|",
            "| entries read from all mapping tables | %d |" % result["n_entries"],
            "| unified loci after cross-mapping collapse | %d |"
            % result["n_unified"],
            "| called (support >= %d in >= %d samples) | %d |"
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
    L = ["## Unresolved loci", "",
         "Listed once each, attributed to no construct. These are the "
         "validation priorities (`references/validation.md`).", "",
         "| Locus | candidates | support | samples | why |",
         "|---|---|---|---|---|"]
    for u in unres:
        L.append("| %s:%d-%d | %s | %d | %d | %s |"
                 % (u["contig"], u["start"], u["end"],
                    ", ".join(u["candidates"]) or "none",
                    u["total_support"], u["n_samples"], u["note"]))
    return L + [""]


def render_md(result: Dict) -> str:
    L = _md_header(result)
    if result["caveats"]:
        L += ["## Read this before using any count below", ""]
        L += ["> %s" % c for c in result["caveats"]] + [""]
    L += _md_per_construct(result)
    L += _md_collapse(result)
    L += _md_unresolved(result)
    return "\n".join(L)


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


def preflight_caveats(a: argparse.Namespace, tables: Dict[str, str],
                      samples: Sequence[str], per_table: Dict[str, List[str]],
                      copies: Dict[str, Optional[float]], assumed_ploidy: bool
                      ) -> List[str]:
    """Everything the inputs cannot deliver, said before any count is shown.

    Each of these silently weakens a site count if it is not stated: an
    assumed ploidy scales every ceiling, a missing copy number removes a
    construct from the ceiling tier entirely, and an impossible
    `--min-samples` returns an empty site list that looks like a result.
    """
    out: List[str] = []
    if assumed_ploidy:
        out.append(
            "PLOIDY ASSUMED to be 2. It multiplies every copy-number ceiling, "
            "so it scales the upper bound of every site count directly. Many "
            "immortalised lines are aneuploid and drift on passage: measure "
            "it and pass --ploidy.")
    odd = {m: ss for m, ss in per_table.items() if list(ss) != list(samples)}
    if odd:
        out.append(
            "The mapping tables do not carry the same samples (%s). Support "
            "is combined over the union, so a locus can look sample-specific "
            "only because one mapping run never saw that sample."
            % "; ".join("%s: %s" % (m, ", ".join(ss) or "none")
                        for m, ss in odd.items()))
    if a.min_samples > len(samples):
        out.append(
            "--min-samples is %d but only %d sample(s) were supplied, so NO "
            "locus can be called. Lower it deliberately (and say so), or add "
            "the other samples' loci tables."
            % (a.min_samples, len(samples)))
    ends = parse_ends(a.vector_ends)
    missing_ends = [m for m in tables if m not in ends or not ends[m]]
    if missing_ends:
        out.append(
            "No --vector-ends for %s, so the vector-end tier could not be "
            "tested there. Loci that sit at those constructs' declared ends "
            "fall through to the copy-number ceiling or to unresolved."
            % ", ".join(sorted(missing_ends)))
    no_cn = [n for n in tables if copies.get(n) is None]
    if no_cn:
        out.append(
            "No copy number for %s, so it has no ceiling and the "
            "copy-ceiling tier will not attribute anything to it. Supply "
            "copy_number.py --json output, or --copies %s=N."
            % (", ".join(sorted(no_cn)), sorted(no_cn)[0]))
    return out


def build_result(a: argparse.Namespace, tables: Dict[str, str],
                 samples: Sequence[str], cn_how: str, ploidy: float,
                 implied_ploidy: Optional[float], entries: Sequence[Dict],
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
        "parameters": {"window": a.window, "min_reads": a.min_reads,
                       "min_samples": a.min_samples, "end_tol": a.end_tol,
                       "ploidy": ploidy,
                       "min_copies_per_site": a.min_copies_per_site,
                       "cn_margin": a.cn_margin,
                       "ploidy_source": ("--ploidy" if a.ploidy is not None
                                         else "copy-number JSON"
                                         if implied_ploidy else "assumed 2")},
        "n_entries": len(entries), "n_unified": len(unified),
        "n_called": len(called), "n_sub": len(sub),
        "ceilings": caps,
        "attributed": used,
        "per_construct": per_construct,
        "unresolved": [{"contig": L["contig"], "start": L["start"],
                        "end": L["end"], "candidates": L["candidates"],
                        "total_support": L["total_support"],
                        "n_samples": L["n_samples"],
                        "note": L["attribution_note"]}
                       for L in called if not L["construct"]],
        "caveats": list(caveats),
    }


def write_outputs(a: argparse.Namespace, result: Dict, called: Sequence[Dict],
                  sub: Sequence[Dict], samples: Sequence[str]) -> None:
    if a.out_called:
        write_table(a.out_called, called, samples)
    if a.out_subthreshold:
        write_table(a.out_subthreshold, sub, samples, extra=["reason"])
    md = render_md(result)
    if a.md:
        with open(a.md, "w", encoding="utf-8") as fh:
            fh.write(md)
    if a.summary:
        with open(a.summary, "w", encoding="utf-8") as fh:
            json.dump(result, fh, indent=2)
    if not a.md and not a.summary:
        print(md)


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--identity", required=True)
    ap.add_argument("--loci", action="append", required=True,
                    metavar="MAPPING=PATH",
                    help="loci table from classify_junctions.py --out-loci, "
                         "one per mapping run; repeatable")
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
    return ap


def main(argv=None) -> int:
    a = build_parser().parse_args(argv)

    try:
        ident = load_json(a.identity)
        tables = parse_kv(a.loci, "--loci")
        ends = parse_ends(a.vector_ends)
        explicit = parse_kv(a.copies, "--copies")
    except (OSError, ValueError) as exc:
        print("error: %s" % exc, file=sys.stderr)
        return 2

    known = list(ident.get("construct_order") or ident["constructs"])
    for name in list(tables) + list(ends) + list(explicit):
        if name not in known:
            print("error: %r is not in the identity map (%s). Every construct "
                  "in the line must be in the map, including ones the "
                  "provider did not map against."
                  % (name, ", ".join(known)), file=sys.stderr)
            return 2

    try:
        entries, samples, per_table = read_all_loci(tables)
    except (OSError, ValueError) as exc:
        print("error: %s" % exc, file=sys.stderr)
        return 2

    try:
        copies, implied_ploidy, cn_how = load_copies(
            a.copy_number or None, explicit, a.cn_sample or None)
    except (OSError, ValueError) as exc:
        print("error: %s" % exc, file=sys.stderr)
        return 2

    ploidy = a.ploidy if a.ploidy is not None else implied_ploidy
    caveats = preflight_caveats(a, tables, samples, per_table, copies,
                                ploidy is None)
    if ploidy is None:
        ploidy = 2.0

    unified = collapse(entries, a.window, samples)
    called, sub = retain(unified, a.min_reads, a.min_samples)
    caps = ceilings(copies, known, ploidy, a.min_copies_per_site)
    used, attr_caveats = attribute(called, known, copies, caps, ends,
                                   a.end_tol, a.cn_margin)
    caveats += attr_caveats

    constructs = [n for n in known if n in tables or copies.get(n) is not None]
    per_construct = site_counts(called, constructs, copies, caps)

    result = build_result(a, tables, samples, cn_how, ploidy, implied_ploidy,
                          entries, unified, called, sub, caps, used,
                          per_construct, caveats)
    write_outputs(a, result, called, sub, samples)

    for c in caveats:
        print("  CAVEAT: %s" % c, file=sys.stderr)
    print("collapsed %d entries -> %d unified loci; %d called, %d "
          "sub-threshold" % (len(entries), len(unified), len(called), len(sub)),
          file=sys.stderr)
    for n, d in per_construct.items():
        print("    %-12s sites %-8s (certain %d, vector-end %d, copy-ceiling "
              "%d; ceiling %s)"
              % (n, d["site_range"], d["certain"], d["vector_end"],
                 d["copy_ceiling"],
                 d["ceiling"] if d["ceiling"] is not None else "unknown"),
              file=sys.stderr)
    nu = len(result["unresolved"])
    print("    unresolved   %d (listed once each, attributed to no construct)"
          % nu, file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
