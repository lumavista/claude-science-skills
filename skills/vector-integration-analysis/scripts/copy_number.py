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
(`samtools depth -a`), each construct's copy number is the depth over the
bases only it carries -- its singleton membership class in the identity map --
read directly, with no solve. Every membership class also gives one equation
(depth over bases carried by constructs S = summed depth of those constructs);
that stacked system is solved UNCONSTRAINED by weighted least squares as a
cross-check, and the shared classes are checked for closure against the
unique-block estimates. A negative solution or a shared-block deficit means
the additive model is incomplete (typically backbone loss in concatemers): it
is reported, never clamped to zero. One `samtools depth` command produces the
input; always ask the provider for it, and for the genome-median depth to use
as the single-genome baseline (`baseline_depth`).

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
    python copy_number.py --run run.json --strict

With `--run` no config file is written by hand: the run manifest holds the
facts, and the identity map, each mapping's classified junctions and the
detected-locus count are taken from what Stages 1 and 2 recorded in it.

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
from via_common import (MODALITIES, load_json, read_depth_file,  # noqa: E402
                        solve_wls)
from via_run import (FAIL, PASS, WARN, Checks, InputError,  # noqa: E402
                     RunError, add_check_arguments, copy_number_config,
                     load_run)

#: Copies per site for the routes that fix it, from the one modality table --
#: so the stages cannot disagree about whether a route can form arrays.
FIXED_COPIES_PER_SITE = {name: 1 for name, m in MODALITIES.items()
                         if m["concatemers"] is False}


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


def choose_baseline(sample_cfg: Dict, coverage: float
                    ) -> Tuple[float, bool, str]:
    """Single-genome baseline for copy number, and a line saying which.

    The genome-median per-base depth is filtered exactly as the vector depth
    is (same BAM, same samtools filters); the read-count estimate is not, and
    typically runs high, biasing every copy number low. That estimate is still
    what the junction conversion factor uses, since junction support is a read
    count.
    """
    bd = sample_cfg.get("baseline_depth")
    if bd is not None and float(bd) > 0:
        return float(bd), True, (
            "copy-number baseline: genome-median per-base depth %.2fx "
            "(baseline_depth; read-count estimate would be %.2fx). The "
            "genome median is preferred because it is filtered identically "
            "to the vector depth." % (float(bd), coverage))
    return coverage, False, (
        "copy-number baseline: read-count estimate %.2fx "
        "(dedup_pairs x 2 x read_length / genome_size). A genome-median "
        "per-base depth (`baseline_depth`) is preferred when available "
        "because it is filtered identically to the vector depth; the "
        "read-count figure typically runs high and biases copy numbers "
        "low." % coverage)


def conversion_factor(coverage: float, read_length: int, anchor: int) -> float:
    """Reads expected to support one junction present at one copy per haploid."""
    return max(0.0, coverage * (read_length - 2 * anchor) / read_length)


# --------------------------------------------------------------------------
# measured route: solve all constructs from per-base depth
# --------------------------------------------------------------------------

TRIM_POLICY = "trim_boundaries"
#: Deprecated `trim_mode` values. `block-edge` is what trim_boundaries already
#: does; `contig-end` is kept as a legacy policy so existing configs still run.
LEGACY_TRIM_MODES = ("block-edge", "contig-end")
NEGATIVE_COPIES_TOL = -0.1      # unconstrained solve below this is "material"
MIN_UNIQUE_BP = 500             # retained unique bp below this is flagged
MAX_UNIQUE_SPREAD = 0.25        # per-block relative range above this is flagged


def _trim_flags(blocks: Sequence[Dict], circular: bool
                ) -> List[Tuple[bool, bool]]:
    """For each block, whether its start and its end should be trimmed.

    An end is trimmed when it abuts a block of a DIFFERENT membership class
    (cross-mapping tapers there) or a contig terminus (coverage falls off
    there). It is not trimmed when it abuts a block of the same membership,
    because depth runs straight through such a join. On a circular reference
    there are no termini: the first and last blocks are neighbours.
    """
    n = len(blocks)
    out: List[Tuple[bool, bool]] = []
    for i, b in enumerate(blocks):
        if i > 0:
            prev = blocks[i - 1]
        else:
            prev = blocks[-1] if circular and n > 1 else None
        if i < n - 1:
            nxt = blocks[i + 1]
        else:
            nxt = blocks[0] if circular and n > 1 else None
        trim_lo = prev is None or prev["membership"] != b["membership"]
        trim_hi = nxt is None or nxt["membership"] != b["membership"]
        out.append((trim_lo, trim_hi))
    return out


def block_depths(ident: Dict, mapping: str, depth: Dict[int, float],
                 trim: int, policy: str = TRIM_POLICY) -> List[Dict]:
    """Mean depth over every block of one vector reference, after trimming.

    Returns one dict per block that keeps at least one base: membership,
    block span, the retained span, mean depth and retained bp. NO block is
    dropped for being small -- a block is lost only when trimming leaves it
    with zero bases. Small blocks then carry little weight because every
    downstream average is bp-weighted, which is the honest way to let them
    count for little rather than for nothing.

    `policy`:

    * ``trim_boundaries`` (default): trim `trim` bases inward only from a block
      end that abuts a different membership class or a contig terminus. Within
      about one read length of a membership boundary a read cannot lie wholly
      inside the block, so cross-mapping tapers and the depth is neither the
      shared nor the unique value; at a linear terminus coverage falls off.
    * ``contig-end`` (DEPRECATED legacy policy): trim only within `trim` of
      the linear reference termini, leaving the boundary taper in.
    """
    c = ident["constructs"][mapping]
    length = int(c["length"])
    blocks = c["blocks"]
    flags = _trim_flags(blocks, bool(c.get("circular", False)))
    out: List[Dict] = []
    for b, (tlo, thi) in zip(blocks, flags):
        a, z = int(b["start"]), int(b["end"])
        if policy == "contig-end":
            lo, hi = max(a, 1 + trim), min(z, length - trim)
        else:
            lo = a + trim if tlo else a
            hi = z - trim if thi else z
        if hi < lo:
            continue
        # A position absent from the depth file is a position at ZERO depth,
        # not a position to exclude from the average. `samtools depth` without
        # -a omits zero-depth bases, so averaging only the present positions
        # silently deletes the evidence of a deletion: 700 uncovered bp of a
        # 3000 bp class reported 2.00 copies/haploid where the truth was 1.53,
        # with no caveat. Sum over the full span and count the span.
        span = hi - lo + 1
        observed = 0
        total = 0.0
        for p in range(lo, hi + 1):
            v = depth.get(p)
            if v is None:
                continue
            observed += 1
            total += v
        if observed == 0:
            continue
        out.append({"membership": list(b["membership"]), "start": a, "end": z,
                    "lo": lo, "hi": hi, "mean": total / span,
                    "bases": span, "observed_bases": observed})
    return out


def bp_weighted(items: Sequence, mean, bases) -> float:
    """bp-weighted mean of `mean(x)` over `items`, weighted by `bases(x)`.

    Every depth average in this module is bp-weighted, so that a block or
    class contributes in proportion to the sequence it actually covers. That
    is what lets small blocks be kept rather than dropped: they count for
    little instead of for nothing.
    """
    bp = sum(bases(x) for x in items)
    return (sum(mean(x) * bases(x) for x in items) / bp) if bp else 0.0


def _collect_classes(cfg: Dict, ident: Dict, sample: str, order: List[str],
                     trim: int, policy: str
                     ) -> Tuple[List[Dict], Dict[str, List[Dict]]]:
    """Read each depth file; return every membership class and unique block.

    Each class is {mapping, membership, row, mean, bases}; every class with a
    non-zero retained bp is kept. Also returns, per construct, its individual
    unique blocks, so the estimate's spread across blocks can be reported.
    """
    files = (cfg["samples"][sample].get("depth_files") or {})
    classes: List[Dict] = []
    ublocks: Dict[str, List[Dict]] = {}
    for mapping, path in files.items():
        if mapping not in ident["constructs"]:
            continue
        p = path if os.path.isabs(path) else os.path.join(cfg["_base"], path)
        try:
            depth = read_depth_file(p, cfg["mappings"].get(mapping, {}).get("contig"))
        except (OSError, ValueError) as exc:
            print("  warning: %s" % exc, file=sys.stderr)
            continue
        grouped: Dict[Tuple[str, ...], List[Dict]] = {}
        for b in block_depths(ident, mapping, depth, trim, policy):
            mem = b["membership"]
            if len(mem) == 1 and mem[0] in order:
                ublocks.setdefault(mem[0], []).append(dict(b, mapping=mapping))
            grouped.setdefault(tuple(mem), []).append(b)
        for mem_t, blocks in grouped.items():
            membership = list(mem_t)
            row = [1.0 if n in membership else 0.0 for n in order]
            nbp = sum(b["bases"] for b in blocks)
            if not nbp or not any(row):
                continue
            classes.append({
                "mapping": mapping, "membership": membership, "row": row,
                "mean": bp_weighted(blocks, _mean, _bases), "bases": int(nbp)})
    return classes, ublocks


def _mean(x: Dict) -> float:
    return x["mean"]


def _bases(x: Dict) -> float:
    return float(x["bases"])


def _unique_depths(classes: List[Dict], order: List[str]) -> Dict[str, float]:
    """bp-weighted mean depth of each construct's singleton class."""
    out: Dict[str, float] = {}
    for n in order:
        hits = [c for c in classes if c["membership"] == [n]]
        if hits:
            out[n] = bp_weighted(hits, _mean, _bases)
    return out


def unique_block_spread(ublocks: Dict[str, List[Dict]], baseline: float
                        ) -> Dict[str, Dict]:
    """Per construct: each unique block's own copy-number estimate, and spread.

    The reported copy number is the bp-weighted mean over these blocks. When
    they disagree -- a partially deleted cassette, a block too short to escape
    the boundary taper, host-homologous sequence -- the mean hides it, so the
    individual estimates are always reported beside it.
    """
    out: Dict[str, Dict] = {}
    base = baseline or 1.0
    for n, bl in ublocks.items():
        rows = [{"mapping": b["mapping"], "start": b["start"], "end": b["end"],
                 "retained": [b["lo"], b["hi"]], "bp": b["bases"],
                 "depth": round(b["mean"], 2),
                 "copies": round(b["mean"] / base, 3)}
                for b in bl]
        bp = sum(r["bp"] for r in rows)
        mean = bp_weighted(bl, _mean, _bases) / base
        cps = [b["mean"] / base for b in bl]
        if len(cps) > 1 and mean > 0:
            var = bp_weighted(bl, lambda b: (b["mean"] / base - mean) ** 2,
                              _bases)
            rel_range = (max(cps) - min(cps)) / mean
            cv = var ** 0.5 / mean
        else:
            rel_range = cv = None
        out[n] = {"blocks": rows, "n_blocks": len(rows), "retained_bp": bp,
                  "copies_bp_weighted": round(mean, 3),
                  "min": round(min(cps), 3) if cps else None,
                  "max": round(max(cps), 3) if cps else None,
                  "relative_range": (round(rel_range, 3)
                                     if rel_range is not None else None),
                  "cv_bp_weighted": round(cv, 3) if cv is not None else None}
    return out


def _choose_depths(order: List[str], unique: Dict[str, float],
                   wls: Optional[Dict[str, float]], baseline: float,
                   caveats: List[str]
                   ) -> Tuple[Dict[str, Optional[float]], Dict[str, str]]:
    """Unique-block depth where there is one; otherwise the flagged solve."""
    depths: Dict[str, Optional[float]] = {}
    source: Dict[str, str] = {}
    for n in order:
        if n in unique:
            depths[n], source[n] = unique[n], "unique block"
            continue
        v = wls[n]                                    # type: ignore[index]
        if baseline and v / baseline < NEGATIVE_COPIES_TOL:
            depths[n] = None
            source[n] = "none: no usable unique block, and the solve is negative"
            caveats.append(
                "%s has no usable unique block, and the stacked solve gives it "
                "%.2f copies -- negative, so no copy number is reported for it "
                "rather than a clamped zero. See the closure check below."
                % (n, v / baseline))
        else:
            depths[n] = v
            source[n] = "WLS solve (FLAGGED: no usable unique block)"
            caveats.append(
                "%s has no usable unique block (none in the depth supplied, or "
                "none left after trimming), so its copy number comes from "
                "the stacked weighted least-squares solve "
                "over shared classes. That inherits any non-additivity of the "
                "shared blocks (e.g. backbone loss); treat it as provisional."
                % n)
    return depths, source


def _thin_unique_caveat(n: str, bp: int, trim: int, policy: str) -> str:
    return ("%s's copy number rests on only %d bp of unique sequence after "
            "trimming (%s, %d bp). The estimate is low-confidence: report it "
            "as a flagged answer, and compare the per-block estimates below "
            "(references/quantification.md, section 2)."
            % (n, bp, policy, trim))


def _spread_caveat(n: str, sp: Dict, limit: float) -> str:
    return ("%s's unique blocks disagree: %d blocks give %s copies per "
            "haploid (relative range %.0f%%, above the %.0f%% flag). The "
            "reported %.2f is their bp-weighted mean. Disagreeing blocks mean "
            "part of the construct is at a different copy number (partial "
            "deletion, truncated arrays) or a block's depth is biased "
            "(taper, host homology); report the range, not just the mean."
            % (n, sp["n_blocks"],
               ", ".join("%.2f" % b["copies"] for b in sp["blocks"]),
               100 * sp["relative_range"], 100 * limit,
               sp["copies_bp_weighted"]))


def _unique_block_caveats(order: Sequence[str], unique: Dict[str, float],
                          spread: Dict[str, Dict], unique_bp: Dict[str, int],
                          trim: int, policy: str, min_unique_bp: int,
                          max_unique_spread: float) -> List[str]:
    """Flag a construct resting on little unique sequence, or on blocks that
    disagree. Either makes the estimate low-confidence, and both are reported
    rather than left for the reader to notice."""
    out: List[str] = []
    for n in order:
        if n not in unique:
            continue
        if unique_bp[n] < min_unique_bp:
            out.append(_thin_unique_caveat(n, unique_bp[n], trim, policy))
        rr = spread[n]["relative_range"]
        if rr is not None and rr > max_unique_spread:
            out.append(_spread_caveat(n, spread[n], max_unique_spread))
    return out


def _closure(classes: List[Dict], depths: Dict[str, Optional[float]],
             order: List[str]) -> Optional[Dict]:
    """Shared-class depth against the sum of the members' unique depths."""
    rows: List[Dict] = []
    for c in classes:
        mem = [m for m in c["membership"] if m in order]
        if len(c["membership"]) < 2 or any(depths.get(m) is None for m in mem):
            continue
        pred = sum(depths[m] for m in mem)            # type: ignore[misc]
        if pred <= 0:
            continue
        rows.append({"mapping": c["mapping"], "membership": c["membership"],
                     "bases": c["bases"], "observed": round(c["mean"], 2),
                     "predicted": round(pred, 2),
                     "deficit_pct": round(100 * (c["mean"] - pred) / pred, 1)})
    if not rows:
        return None
    bp = sum(r["bases"] for r in rows)
    ob = sum(r["observed"] * r["bases"] for r in rows)
    pr = sum(r["predicted"] * r["bases"] for r in rows)
    return {"classes": rows, "observed": round(ob / bp, 2),
            "predicted": round(pr / bp, 2),
            "deficit_pct": round(100 * (ob - pr) / pr, 1)}


def _wls_residuals(classes: List[Dict], x: List[float]
                   ) -> Tuple[List[Dict], float]:
    """Per-class residuals of the stacked solve, and weighted relative RSS."""
    resid: List[Dict] = []
    rss = denom = 0.0
    for c in classes:
        pred = sum(r * v for r, v in zip(c["row"], x))
        resid.append({"mapping": c["mapping"], "membership": c["membership"],
                      "mean_depth": round(c["mean"], 2), "bases": c["bases"],
                      "predicted_depth": round(pred, 2),
                      "residual": round(c["mean"] - pred, 2)})
        rss += c["bases"] * (c["mean"] - pred) ** 2
        denom += c["bases"] * c["mean"] ** 2
    return resid, round(rss / (denom or 1.0), 4)


def _negative_caveat(wls: Dict[str, float], baseline: float,
                     closure: Optional[Dict]) -> Optional[str]:
    """A materially negative solve means the additive model is incomplete."""
    negative = {n: v / baseline for n, v in wls.items()
                if baseline and v / baseline < NEGATIVE_COPIES_TOL}
    if not negative:
        return None
    return (
        "ADDITIVE MODEL INCOMPLETE. The unconstrained stacked solve "
        "gives %s copies per haploid -- below zero, which is not "
        "physical. Shared blocks carry %s depth than the sum of their "
        "constructs' unique blocks predicts (closure %s%%), most likely "
        "partial loss of shared backbone in concatemers. Copy numbers "
        "below are each construct's own unique-block depth, which does "
        "not depend on the shared blocks; the solve is shown only as a "
        "diagnostic and has NOT been clamped to zero."
        % (", ".join("%s %.2f" % kv for kv in negative.items()),
           "less" if closure and closure["deficit_pct"] < 0 else "different",
           "%+.1f" % closure["deficit_pct"] if closure else "n/a"))


def solve_measured(cfg: Dict, ident: Dict, sample: str, trim: int,
                   baseline: float, policy: str = TRIM_POLICY,
                   min_unique_bp: int = MIN_UNIQUE_BP,
                   max_unique_spread: float = MAX_UNIQUE_SPREAD
                   ) -> Optional[Dict]:
    """Per-construct depth from per-base depth files.

    The ESTIMATOR is each construct's own singleton membership class -- the
    bases only that construct carries -- read directly, with no solve. The
    all-class weighted least-squares solve is computed UNCONSTRAINED and kept
    as a consistency check: it is used for a construct only when that
    construct has no usable unique block, and that is flagged.

    Shared classes are then checked for closure: observed depth against the
    sum of their members' unique-block depths. A deficit there, and a
    negative unconstrained solution, both mean the additive model is
    incomplete -- most often partial loss of shared backbone in concatemers.
    They are reported, never clamped away.

    Each construct's individual unique-block estimates and their spread are
    reported too, with a caveat when the retained unique bp is low or the
    blocks disagree.
    """
    order: List[str] = list(cfg["mappings"])
    classes, ublocks = _collect_classes(cfg, ident, sample, order, trim, policy)
    if not classes:
        return None
    unique = _unique_depths(classes, order)
    x = (solve_wls([c["row"] for c in classes], [c["mean"] for c in classes],
                   [float(c["bases"]) for c in classes], non_negative=False)
         if len(classes) >= len(order) else None)
    wls = dict(zip(order, x)) if x is not None else None
    if wls is None and any(n not in unique for n in order):
        return None          # some construct has neither: fall back, as before

    caveats: List[str] = []
    # Positions missing from the depth file are now averaged in as zero, which
    # is correct, but a file with many missing positions still means the caller
    # ran `samtools depth` without -a. Say so: the user cannot otherwise tell a
    # genuinely uncovered region from an omitted one, and that distinction is
    # the difference between a deletion and a reporting artefact.
    span_bp = sum(int(c.get("bases") or 0) for c in classes)
    seen_bp = sum(int(c.get("observed_bases", c.get("bases") or 0))
                  for c in classes)
    if span_bp and seen_bp < 0.99 * span_bp:
        caveats.append(
            "DEPTH FILE INCOMPLETE: %d of %d positions in the measured classes "
            "(%.1f%%) are absent from the depth file and were averaged in as "
            "zero depth, which is what absence means. If the file was produced "
            "without `samtools depth -a`, zero-coverage bases were omitted and "
            "these means are correct only if every omitted base is genuinely "
            "uncovered. Re-run with -a to remove the ambiguity; until then a "
            "low mean here cannot be separated from a reporting artefact."
            % (span_bp - seen_bp, span_bp, 100.0 * (span_bp - seen_bp) / span_bp))
    depths, source = _choose_depths(order, unique, wls, baseline, caveats)
    spread = unique_block_spread(ublocks, baseline)
    unique_bp = {n: spread[n]["retained_bp"] if n in spread else 0 for n in order}
    caveats += _unique_block_caveats(order, unique, spread, unique_bp, trim,
                                     policy, min_unique_bp, max_unique_spread)
    closure = _closure(classes, depths, order)
    resid, rss_rel = _wls_residuals(classes, x) if x is not None else ([], None)
    neg = _negative_caveat(wls, baseline, closure) if wls is not None else None
    if neg:
        caveats.insert(0, neg)

    return {"method": "measured (unique-block depth; weighted least squares "
                      "as cross-check)",
            "depths": {n: (round(v, 2) if v is not None else None)
                       for n, v in depths.items()},
            "source": source,
            "unique_bp": unique_bp,
            "unique_blocks": {n: spread[n] for n in order if n in spread},
            "trim_policy": policy, "trim_boundaries": trim,
            # legacy keys, kept so existing consumers of the JSON still work
            "trim_mode": policy, "edge_trim": trim,
            "wls_depths": ({n: round(v, 2) for n, v in wls.items()}
                           if wls is not None else None),
            "wls_note": (None if wls is not None else
                         "stacked solve unavailable (singular or too few "
                         "classes); unique blocks need no solve"),
            "equations": len(classes), "unknowns": len(order),
            "weighted_relative_rss": rss_rel,
            "classes_used": resid,
            "closure": closure,
            "caveats": caveats}


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
                 sample: str, copies: Optional[float], modality: str,
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

    capable = MODALITIES[modality]["concatemers"]
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
        if not internal:
            # 1 + 2*(0/flank) is 1 whatever the line is. On a concatemer-
            # capable modality that is the absence of evidence for arrays, not
            # evidence of single copies, and the two look identical in the
            # number alone.
            out["notes"].append(
                "Copies per site came out at 1 because no vector-vector "
                "junction survived flagging, not because the sites were shown "
                "to carry one copy each: 1 + 2*(0/flank) is 1 by "
                "construction. On %s delivery arrays are possible, so treat "
                "this as a LOWER bound on copies per site and the implied "
                "site count as an upper bound. Internal junctions can be "
                "absent because the array junctions are religation junctions "
                "excluded by default (--include-end-religation), because "
                "fillers exceed what a read spans, or because the line really "
                "is single-copy per site; long reads separate these."
                % modality)
    else:
        out["copies_per_site"] = None
        out["copies_per_site_source"] = (
            "cannot be estimated: " + ("no host-vector flank reads" if not flank
                                       else "modality undeclared, so it is "
                                       "unknown whether arrays are possible"))

    n = out["copies_per_site"]
    # `copies` is per HAPLOID genome; a site count is per cell. Stage 4's
    # ceiling uses ceil(copies * ploidy / divisor) (reconcile_sites.caps) and
    # SKILL.md states that form, so omitting ploidy here made implied_sites a
    # factor of `ploidy` too small: a correct, self-consistent diploid line
    # (4 single-copy loci = 2.0 copies/haploid) was reported as
    # "EXCESS: 4 loci detected but the measured DNA only supports 2.0 (0.5x)",
    # i.e. a correct result flagged as a 2x double-counting artefact, with
    # remediation advice pointing at flank pairing.
    ploidy = float(cfg.get("ploidy", 2))
    out["implied_sites"] = (round(copies * ploidy / n, 2)
                            if n and copies is not None else None)
    return out


def flank_detection(cfg: Dict, coverage: float, flank_reads: int,
                    implied_sites: Optional[float],
                    anchors: Sequence[int]) -> List[Dict]:
    sweep = []
    for anchor in anchors:
        R = conversion_factor(coverage, cfg["read_length"], anchor)
        detected = flank_reads / R if R else None
        rate = (detected / (2 * implied_sites)
                if detected and implied_sites else None)
        entry = {
            "anchor": anchor,
            "reads_per_copy": round(R, 2),
            "flank_copies_detected": round(detected, 2) if detected else None,
            "flank_detection_rate": round(rate, 3) if rate is not None else None,
        }
        # A detection rate above 1 means more flanks were detected than the
        # structure admits (2 per site). That is a unit or assumption
        # inconsistency, not a sensitivity figure -- say so rather than
        # printing an impossible fraction.
        if rate is not None and rate > 1.05:
            entry["note"] = (
                "detection rate %.2f exceeds 1: more flank copies detected "
                "than 2 per implied site allows. Suspect the anchor "
                "assumption, the conversion factor, or an implied-site count "
                "in the wrong unit." % rate
            )
        sweep.append(entry)
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


def calibrate(cfg: Dict, jrows: Dict[str, List[Dict]], coverage: float,
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
        R = conversion_factor(coverage, cfg["read_length"], anchor)
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

def _render_closure(cl: Optional[Dict]) -> List[str]:
    if not cl:
        return []
    L: List[str] = []
    L += ["**Shared-block closure: %+.1f%%** (observed %.2f vs %.2f "
          "predicted from the members' unique blocks, bp-weighted)."
          % (cl["deficit_pct"], cl["observed"], cl["predicted"]), "",
          "| Reference | Bases carried by | bp used | observed | "
          "predicted (sum of unique) | deficit |",
          "|---|---|---|---|---|---|"]
    for r in cl["classes"]:
        L.append("| %s | %s | %s | %.2f | %.2f | %+.1f%% |"
                 % (r["mapping"], " + ".join(r["membership"]),
                    format(r["bases"], ","), r["observed"],
                    r["predicted"], r["deficit_pct"]))
    L += ["", "A shared-block deficit means the constructs' depths do "
          "not add over shared sequence. In concatemer lines the usual "
          "cause is partial loss of shared backbone; a surplus points to "
          "host-homologous sequence or an unlisted construct.", ""]
    return L


def _render_unique_blocks(ub: Dict[str, Dict], order: Sequence[str]) -> List[str]:
    """Each construct's individual unique-block estimates, beside their mean.

    The mean alone hides a construct whose blocks disagree, which is a real
    result (a partially deleted cassette, truncated arrays) and not noise.
    """
    if not ub:
        return []
    L = ["Per-unique-block estimates (the reported copy number is their "
         "bp-weighted mean; a wide spread is flagged above):", "",
         "| Construct | block | retained bp | copies | spread |",
         "|---|---|---|---|---|"]
    for n in order:
        sp = ub.get(n)
        if not sp:
            continue
        spread = ("%.2f-%.2f (range %.0f%%, CV %.0f%%)"
                  % (sp["min"], sp["max"], 100 * sp["relative_range"],
                     100 * sp["cv_bp_weighted"])
                  if sp["relative_range"] is not None else "single block")
        for i, b in enumerate(sp["blocks"]):
            L.append("| %s | %s:%d-%d | %s | %.2f | %s |"
                     % (n, b["mapping"], b["start"], b["end"],
                        format(b["bp"], ","), b["copies"],
                        spread if i == 0 else ""))
        L.append("| %s | **total** | **%s** | **%.2f** | |"
                 % (n, format(sp["retained_bp"], ","),
                    sp["copies_bp_weighted"]))
    return L + [""]


def render_measured(s: str, d: Dict, order: List[str]) -> List[str]:
    """Per-sample report section for the measured route: estimates, the
    WLS cross-check, and shared-block closure. Empty for other routes."""
    sol = d["solve"]
    if "source" not in sol:
        return []
    lines: List[str] = []
    baseline_depth = d["baseline_depth"] or 1.0
    lines += ["### %s -- per-construct estimate" % s, "",
          "Each construct's copy number is the depth over the bases only it "
          "carries, read directly (trim %d bp, policy `%s`; every block with "
          "bases left after trimming is kept, bp-weighted). The stacked "
          "weighted least-squares solve over every class is shown beside "
          "it, unconstrained, as a cross-check only."
          % (sol["trim_boundaries"], sol["trim_policy"]), "",
          "| Construct | reported copies | source | unique bp | "
          "WLS cross-check copies |",
          "|---|---|---|---|---|"]
    wls = sol.get("wls_depths") or {}
    for n in order:
        c = d["copies_per_haploid"].get(n)
        lines.append("| %s | %s | %s | %s | %s |"
                 % (n, "%.2f" % c if c is not None else "not reported",
                    sol["source"].get(n, "-"),
                    format(sol["unique_bp"].get(n, 0), ","),
                    "%.2f" % (wls[n] / baseline_depth) if n in wls else "-"))
    lines.append("")
    lines += _render_unique_blocks(sol.get("unique_blocks") or {}, order)
    if sol.get("wls_note"):
        lines += ["Stacked solve: %s." % sol["wls_note"], ""]
    lines += _render_closure(sol.get("closure"))
    if sol.get("classes_used"):
        lines += ["Stacked solve: %d equations, %d unknowns, weighted relative "
              "residual %s." % (sol["equations"], sol["unknowns"],
                                sol["weighted_relative_rss"]),
              "",
              "| Reference | Bases carried by | bp used | observed | WLS predicted | residual |",
              "|---|---|---|---|---|---|"]
        for c in sol["classes_used"]:
            lines.append("| %s | %s | %s | %.2f | %.2f | %+.2f |"
                     % (c["mapping"], " + ".join(c["membership"]),
                        format(c["bases"], ","), c["mean_depth"],
                        c["predicted_depth"], c["residual"]))
        lines += ["", "Large residuals mean the block map does not explain the "
              "depth: a repeat shared with the host genome, an unlisted "
              "construct, non-uniform library coverage, or shared-backbone "
              "loss.", ""]
    return lines


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
        L += ["All copy numbers below are measured from per-base depth, not "
              "inferred.", ""]

    L += ["## Construct abundance", "",
          "Units: copies per haploid genome equivalent. Multiply by the ploidy "
          "(%s) for copies per cell, and treat ploidy as measured rather than "
          "assumed -- many immortalised lines are aneuploid and drift."
          % cfg.get("ploidy", "?"), ""]
    order = list(cfg["mappings"])
    L += ["| Sample | Genome cov (reads) | CN baseline | " + " | ".join(order)
          + " | Route |",
          "|---" * (len(order) + 4) + "|"]
    for s, d in result["samples"].items():
        cells = [s, "%.1fx" % d["genome_coverage"],
                 "%.1fx (%s)" % (d["baseline_depth"], d["baseline_source"])]
        for n in order:
            c = d["copies_per_haploid"].get(n)
            cells.append("%.2f" % c if c is not None else "-")
        cells.append(d["solve"]["method"])
        L.append("| " + " | ".join(cells) + " |")
    L.append("")
    for s, d in result["samples"].items():
        L.append("- **%s** %s" % (s, d["baseline_note"]))
    L.append("")

    for s, d in result["samples"].items():
        L += render_measured(s, d, order)

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


def stage_checks(checks: Checks, result: Dict,
                 skipped_samples: Sequence[str] = ()) -> Checks:
    """The conditions under which these copy numbers must not be reported."""
    primary = result["primary_mapping"]
    samples = result["samples"]

    if result["modality"] == "auto":
        checks.add("modality_declared", FAIL,
                   "the delivery modality is undeclared, so copies per site -- "
                   "and with it the implied site count -- cannot be chosen. "
                   "Settle it from the Stage 2 fingerprint and declare it.")
    else:
        checks.add("modality_declared", PASS, "modality %s" % result["modality"])

    routes = {s: d["solve"]["method"] for s, d in samples.items()}
    per_sample = {s: (PASS if m.startswith("measured") else
                      WARN if m.startswith("inferred") else FAIL)
                  for s, m in routes.items()}
    listing = "; ".join("%s: %s" % kv for kv in routes.items())
    measured = [s for s, v in per_sample.items() if v == PASS]
    # A sample with no depth does not invalidate the samples that have it, so
    # the check fails only when NO sample was measured. Where some were, it
    # warns and names which, because a reader will otherwise read one route
    # label as covering the whole table.
    if not samples:
        checks.add("copy_number_route", WARN, "no samples")
    elif not measured:
        checks.add("copy_number_route", FAIL,
                   "copy number was measured in no sample (%s). The figures "
                   "are upper bounds or absent, and everything derived from "
                   "them inherits that. Request `samtools depth -a` per vector "
                   "reference (references/quantification.md)." % listing)
    elif len(measured) < len(samples):
        checks.add("copy_number_route", WARN,
                   "copy number is measured from per-base depth in %s and not "
                   "in %s (%s). Quote the measured samples as measured and the "
                   "others as bounds; a sample with no depth is still usable "
                   "for junction support in Stage 4."
                   % (", ".join(measured),
                      ", ".join(s for s in samples if s not in measured),
                      listing))
    else:
        checks.add("copy_number_route", PASS, listing)

    mismatched = [s for s, d in samples.items()
                  if any("does not appear in the junction table" in n
                         for n in d["architecture"]["notes"])]
    if mismatched:
        checks.add("sample_names", FAIL,
                   "sample(s) %s are not in the junction table for %s, so "
                   "flank and internal read counts are zero and copies per "
                   "site, implied sites and the reconciliation are "
                   "meaningless. Sample names must match the names given to "
                   "classify_junctions.py." % (", ".join(mismatched), primary))
    else:
        checks.add("sample_names", PASS,
                   "every sample is present in the junction table")

    withheld = [s for s, d in samples.items()
                if d["copies_per_haploid"].get(primary) is None]
    incomplete = any("ADDITIVE MODEL INCOMPLETE" in c
                     for c in result["caveats"])
    notes = []
    if withheld:
        notes.append("no copy number is reported for %s in %s"
                     % (primary, ", ".join(withheld)))
    if incomplete:
        notes.append("the additive model is incomplete (shared-block deficit "
                     "or negative solve); copy numbers are the unique-block "
                     "estimates")
    if skipped_samples:
        notes.append("sample(s) %s declare no dedup_pairs and were left out "
                     "of the copy-number stage" % ", ".join(skipped_samples))
    other = len(result["caveats"])
    checks.add("estimate_caveats", WARN if (notes or other) else PASS,
               ("; ".join(notes) + ". " if notes else "")
               + "%d caveat(s) at the top of the report; quote them with the "
                 "numbers." % other if (notes or other) else "no caveats")

    verdicts = {s: d["reconciliation"]["verdict"].split(":")[0]
                for s, d in samples.items()}
    agree = all(v.startswith("AGREES") for v in verdicts.values())
    checks.add("reconciliation", PASS if agree and verdicts else WARN,
               "; ".join("%s: %s" % kv for kv in verdicts.items())
               or "no samples")
    return checks


def copies_per_site_summary(result: Dict) -> Optional[float]:
    """The largest measured copies per site over the samples, for Stage 4."""
    vals = [d["architecture"]["copies_per_site"]
            for d in result["samples"].values()
            if d["architecture"]["copies_per_site"]]
    return max(vals) if vals else None


def trim_from_config(cfg: Dict, default: int = 150
                     ) -> Tuple[int, List[str]]:
    """The boundary trim a copy-number config asks for, and any deprecation
    notes it earned. Used by `resolve_trim` here and by make_figures.py, so a
    coverage track is trimmed exactly as the copy number was."""
    notes: List[str] = []
    if cfg.get("trim_boundaries") is not None:
        return int(cfg["trim_boundaries"]), notes
    if cfg.get("edge_trim") is not None:
        notes.append("config key `edge_trim` is deprecated; rename it "
                     "`trim_boundaries`.")
        return int(cfg["edge_trim"]), notes
    return int(cfg.get("read_length") or default), notes


def resolve_trim(args: argparse.Namespace, cfg: Dict
                 ) -> Tuple[int, str, List[str]]:
    """Trim length and policy from CLI and config, mapping deprecated keys.

    New: `--trim-boundaries` / `trim_boundaries`. Old, still accepted with a
    deprecation warning: `--edge-trim` / `edge_trim`, `--trim-mode` /
    `trim_mode`, `--min-class-bp` / `min_class_bp`.
    """
    notes: List[str] = []
    trim = args.trim_boundaries
    if trim is None and args.edge_trim:
        trim = args.edge_trim
        notes.append("--edge-trim is deprecated; use --trim-boundaries "
                     "(same meaning).")
    if not trim:
        trim, cfg_notes = trim_from_config(cfg)
        notes += cfg_notes

    mode = args.trim_mode
    if mode is not None:
        notes.append("--trim-mode is deprecated; trimming now follows the "
                     "trim_boundaries policy.")
    elif "trim_mode" in cfg:
        mode = cfg["trim_mode"]
        notes.append("config key `trim_mode` is deprecated; delete it -- "
                     "trimming now follows the trim_boundaries policy.")
    if mode is not None and mode not in LEGACY_TRIM_MODES:
        raise ValueError("trim_mode %r unknown; expected one of %s"
                         % (mode, ", ".join(LEGACY_TRIM_MODES)))
    policy = TRIM_POLICY
    if mode == "contig-end":
        policy = "contig-end"
        notes.append("trim_mode contig-end runs as a LEGACY policy: it leaves "
                     "the cross-mapping taper at membership boundaries in. "
                     "Prefer the default trim_boundaries policy, which no "
                     "longer drops small blocks.")
    if args.min_class_bp is not None or "min_class_bp" in cfg:
        notes.append("min_class_bp is deprecated and ignored: every class with "
                     "retained bases is used, bp-weighted.")
    return trim, policy, notes


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("config", nargs="?", default="",
                    help="config JSON (omit it when --run supplies the run "
                         "manifest)")
    ap.add_argument("--md")
    ap.add_argument("--json")
    ap.add_argument("--trim-boundaries", type=int, default=None,
                    help="bases trimmed inward from a block end that abuts a "
                         "different membership class or a contig terminus, "
                         "before averaging depth. Ends abutting a same-"
                         "membership block are not trimmed. Default: "
                         "`trim_boundaries` in the config, else read_length.")
    ap.add_argument("--edge-trim", type=int, default=None,
                    help="DEPRECATED alias of --trim-boundaries")
    ap.add_argument("--trim-mode", choices=LEGACY_TRIM_MODES, default=None,
                    help="DEPRECATED. block-edge is the trim_boundaries "
                         "policy; contig-end is kept as a legacy policy that "
                         "trims only the reference termini.")
    ap.add_argument("--min-unique-bp", type=int, default=MIN_UNIQUE_BP,
                    help="flag any construct whose retained unique bp after "
                         "trimming is below this (default %d)" % MIN_UNIQUE_BP)
    ap.add_argument("--max-unique-spread", type=float, default=MAX_UNIQUE_SPREAD,
                    help="flag any construct whose per-unique-block copy "
                         "numbers span more than this fraction of their "
                         "bp-weighted mean (default %.2f)" % MAX_UNIQUE_SPREAD)
    ap.add_argument("--min-class-bp", type=int, default=None,
                    help="DEPRECATED and ignored: no block or class is dropped "
                         "on size any more; every class with retained bases "
                         "is used, bp-weighted.")
    ap.add_argument("--anchors", default="15,20,25",
                    help="aligner anchor lengths to sweep")
    ap.add_argument("--include-end-religation", action="store_true",
                    help="count junctions joining the two vector ends as "
                         "internal array junctions. OFF by default: against a "
                         "linear reference they are also produced by uncut "
                         "input molecule and by any co-resident construct whose "
                         "equivalent site is intact, so they are ambiguous.")
    add_check_arguments(ap)
    return ap


def read_config(args: argparse.Namespace
                ) -> Tuple[Optional[object], Dict, List[str]]:
    """The config, from `--config` or the run manifest, with its base
    directory recorded in `_base` and this stage's outputs named."""
    try:
        run = load_run(args.run)
        skipped_samples: List[str] = []
        if args.config:
            cfg = load_json(args.config)
            base = os.path.dirname(os.path.abspath(args.config))
        elif run is not None:
            cfg, skipped_samples = copy_number_config(run)
            base = run.base
        else:
            raise InputError("give a config file, or --run RUN.json")
    except RunError as exc:
        raise InputError(str(exc))
    except (OSError, ValueError) as exc:
        raise InputError("cannot read config: %s" % exc)
    if run is not None:
        args.md = run.output(args.md, "accounting.md")
        args.json = run.output(args.json, "accounting.json")
    cfg = {k: v for k, v in cfg.items() if not k.startswith("_")}
    cfg["_base"] = base
    for key in ("identity", "genome_size", "read_length", "samples",
                "mappings"):
        if key not in cfg:
            raise InputError("config is missing %r%s" % (
                key, " (declare it in the manifest, or run host_reference.py "
                     "--run first to take it from the host assembly)"
                if key == "genome_size" else ""))
    return run, cfg, skipped_samples


def read_stage_inputs(cfg: Dict) -> Tuple[Dict, Dict[str, List[Dict]], str,
                                          str]:
    """The identity map, each mapping's classified junctions, the primary
    mapping and the modality -- validated against each other."""
    base = cfg["_base"]
    try:
        ident = load_json(cfg["identity"] if os.path.isabs(cfg["identity"])
                          else os.path.join(base, cfg["identity"]))
    except (OSError, ValueError) as exc:
        raise InputError("cannot read identity map: %s" % exc)

    unknown = [n for n in cfg["mappings"] if n not in ident["constructs"]]
    if unknown:
        raise InputError(
            "mapping(s) %s are not in the identity map (%s). The identity map "
            "must cover every construct in the line, including ones the "
            "provider did not map against."
            % (unknown, ", ".join(ident["constructs"])))

    jrows: Dict[str, List[Dict]] = {}
    for name, m in cfg["mappings"].items():
        path = m.get("junctions")
        if not path:
            continue
        try:
            jrows[name] = read_classified(
                path if os.path.isabs(path) else os.path.join(base, path))
        except (OSError, ValueError) as exc:
            raise InputError(str(exc))

    primary = cfg.get("primary_mapping", list(cfg["mappings"])[0])
    if primary not in cfg["mappings"]:
        raise InputError("primary_mapping %r is not one of %s"
                         % (primary, ", ".join(cfg["mappings"])))
    modality = cfg.get("modality", "auto")
    if modality not in MODALITIES:
        raise InputError("modality %r unknown; expected one of %s"
                         % (modality, ", ".join(sorted(MODALITIES))))
    return ident, jrows, primary, modality


def measure_sample(args: argparse.Namespace, cfg: Dict, ident: Dict,
                   jrows: Dict[str, List[Dict]], sample: str, primary: str,
                   modality: str, trim: int, policy: str,
                   anchors: Sequence[int]) -> Tuple[Dict, List[str]]:
    """One sample's copy numbers, architecture and reconciliation, plus the
    caveats it earned. Reads nothing and writes nothing."""
    sample_cfg = cfg["samples"][sample]
    if "dedup_pairs" not in sample_cfg:
        raise InputError("sample %r has no dedup_pairs" % sample)
    coverage = genome_coverage(float(sample_cfg["dedup_pairs"]),
                               int(cfg["read_length"]),
                               float(cfg["genome_size"]))
    baseline, use_median, baseline_note = choose_baseline(sample_cfg,
                                                          coverage)
    print("%s: %s" % (sample, baseline_note), file=sys.stderr)
    sol = solve_measured(cfg, ident, sample, trim, baseline, policy,
                         args.min_unique_bp, args.max_unique_spread)
    if sol is None:
        sol = solve_inferred(cfg, ident, sample, jrows)
    caveats = ["[%s] %s" % (sample, c) if len(cfg["samples"]) > 1 else c
               for c in sol.get("caveats", [])]

    copies = {n: (v / baseline if baseline and v is not None else None)
              for n, v in sol["depths"].items()}
    cn = copies.get(primary)
    if cn is None and primary in sol["depths"]:
        # measured route withheld the number: do not let it become a 0
        arch_cn: Optional[float] = None
    else:
        arch_cn = cn if cn else 0.0
    arch = architecture(cfg, jrows, primary, sample, arch_cn, modality,
                        args.include_end_religation)
    caveats += [n for n in arch["notes"]
                if "does not appear in the junction table" in n]
    detected = sample_cfg.get("detected_loci")
    if detected is None:
        detected = cfg["mappings"][primary].get("detected_loci")
    return {
        "genome_coverage": round(coverage, 2),
        "baseline_depth": round(baseline, 2),
        "baseline_source": ("genome-median per-base depth" if use_median
                            else "read count"),
        "baseline_note": baseline_note,
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
        "flank_detection": flank_detection(cfg, coverage, arch["flank_reads"],
                                           arch["implied_sites"], anchors),
        "calibration": calibrate(cfg, jrows, coverage, sample, anchors),
    }, caveats


def write_outputs(args: argparse.Namespace, cfg: Dict, result: Dict) -> None:
    md = render(cfg, result)
    if args.md:
        with open(args.md, "w", encoding="utf-8") as fh:
            fh.write(md)
    if args.json:
        with open(args.json, "w", encoding="utf-8") as fh:
            json.dump(result, fh, indent=2)
    if not args.md and not args.json:
        print(md)


def report_stderr(result: Dict, caveats: Sequence[str]) -> None:
    for c in caveats:
        print("  CAVEAT: %s" % c, file=sys.stderr)
    for sample, d in result["samples"].items():
        parts = " | ".join("%s %s" % (n, v) for n, v in
                           d["copies_per_haploid"].items())
        print("%s: genome %.1fx | %s | route: %s"
              % (sample, d["genome_coverage"], parts, d["solve"]["method"]),
              file=sys.stderr)
        arch = d["architecture"]
        print("    copies/site %s -> %s sites implied | %s"
              % (arch["copies_per_site"], arch["implied_sites"],
                 d["reconciliation"]["verdict"].split(":")[0]),
              file=sys.stderr)


def run_stage(args: argparse.Namespace) -> int:
    run, cfg, skipped_samples = read_config(args)
    ident, jrows, primary, modality = read_stage_inputs(cfg)

    try:
        trim, policy, deprecations = resolve_trim(args, cfg)
    except ValueError as exc:
        raise InputError(str(exc))
    for d in deprecations:
        print("  DEPRECATED: %s" % d, file=sys.stderr)
    anchors = [int(x) for x in args.anchors.split(",") if x.strip()]

    result: Dict = {"primary_mapping": primary, "modality": modality,
                    "samples": {}, "caveats": [],
                    "trim": {"policy": policy, "trim_boundaries": trim},
                    "deprecations": deprecations}
    caveats: List[str] = []
    if modality == "auto":
        caveats.append(
            "The delivery modality is undeclared. Whether one site carries "
            "one copy or a tandem array is a property of the mechanism, not "
            "of the data, and it changes the site count several-fold. Settle "
            "it from the Stage 2 fingerprint and set `modality` in the "
            "config.")

    for sample in cfg["samples"]:
        measured, sample_caveats = measure_sample(
            args, cfg, ident, jrows, sample, primary, modality, trim, policy,
            anchors)
        result["samples"][sample] = measured
        caveats += [c for c in sample_caveats if c not in caveats]
        # A read-count baseline biases every copy number LOW (the note says
        # "typically runs high and biases copy numbers low"). That note used to
        # land only on stderr and in a bullet below the section 6 table, so
        # result["caveats"] stayed empty, the markdown claimed "All copy
        # numbers below are measured from per-base depth, not inferred", and
        # report.html's caveat table said no caveat was raised -- directly
        # above the biased number. It is a caveat; record it as one.
        if measured.get("baseline_source") == "read count":
            note = measured.get("baseline_note")
            if note and note not in caveats:
                caveats.append("Sample %s -- %s" % (sample, note))

    result["caveats"] = caveats
    checks = stage_checks(Checks("copy_number", args.accept), result,
                          skipped_samples)
    result["checks"] = checks.as_list()
    write_outputs(args, cfg, result)
    report_stderr(result, caveats)
    if run is not None:
        run.record("copy_number", {
            "json": run.relative(args.json), "md": run.relative(args.md),
            "primary_mapping": primary, "modality": modality,
            "copies_per_haploid": {s: d["copies_per_haploid"]
                                   for s, d in result["samples"].items()},
            "copies_per_site": copies_per_site_summary(result),
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
