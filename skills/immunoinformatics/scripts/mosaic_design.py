#!/usr/bin/env python3
"""Mosaic and cocktail antigen design by k-mer coverage maximisation.

Implements the EpiGraph / mosaic-vaccine coverage objective (Fischer et al.,
Nat Med 2007; Theiler & Korber, Sci Rep 2016). The original LANL web tools
could not be confirmed live as of 2026-10, so the algorithm is implemented
here directly — which also removes the licensing unknowns of depending on a
hosted tool.

SKILL.md states that string-of-beads designs create junctional neo-epitopes at
every join and must be checked. This script is what does that check, and in
`mosaic` mode it constructs sequences that cannot contain them.

Modes
-----
  cocktail  Greedily select N natural panel sequences whose combined k-mer
            content covers the most of the panel. Every component is a real
            strain sequence, so foldability and processing are not in question.

  mosaic    Construct N synthetic sequences by beam search over the panel's
            k-mer de Bruijn graph. Every (k-1)-overlap transition used is one
            observed in a real sequence, so no k-mer absent from the panel can
            be created. That is the junction guarantee: any peptide the
            construct presents is one some real strain presents too.

  evaluate  Score an externally supplied design set for coverage and scan it
            for novel k-mers. Use this on any string-of-beads or fusion
            construct, however it was designed.

Coverage is frequency-weighted by default — a k-mer in 200 strains counts more
than one in 2. `--unweighted` counts distinct k-mers instead, favouring
rare-variant breadth over common-variant depth. Neither is correct in the
abstract; pick deliberately and say which you used.

Writes a sidecar <out>.meta.json for reproducibility.

Usage
-----
  python mosaic_design.py cocktail --fasta panel.fasta -n 4 -k 9 --out work/cocktail
  python mosaic_design.py mosaic   --fasta panel.fasta -n 2 -k 9 --out work/mosaic
  python mosaic_design.py evaluate --fasta panel.fasta --designs c.fasta -k 9 \
      --out work/eval

Input should be UNALIGNED sequences; gaps are stripped if present.
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")
GAPS = set("-.~")
AMBIG = set("XBZJUO*")


def log(msg: str) -> None:
    print(msg, file=sys.stderr)


def die(msg: str, code: int = 1) -> None:
    print(f"error: {msg}", file=sys.stderr)
    raise SystemExit(code)


def write_meta(out: Path, payload: dict) -> None:
    out.with_suffix(out.suffix + ".meta.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n"
    )


def read_fasta(path: Path) -> list[tuple[str, str]]:
    records, name, chunks = [], None, []
    with path.open() as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(chunks)))
                name, chunks = line[1:].split()[0], []
            else:
                chunks.append(line)
    if name is not None:
        records.append((name, "".join(chunks)))
    out = [(n, "".join(c for c in s.upper() if c not in GAPS))
           for n, s in records]
    out = [(n, s) for n, s in out if s]
    if not out:
        die(f"no sequences in {path}")
    return out


def write_fasta(path: Path, named: list[tuple[str, str]]) -> None:
    with path.open("w") as fh:
        for name, seq in named:
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")


def kmers(seq: str, k: int):
    """Yield k-mers, skipping any containing an ambiguity code.

    Ambiguous k-mers are skipped rather than imputed: a k-mer containing X is
    not a real epitope, and counting it inflates numerator and denominator in
    ways that become impossible to reason about downstream.
    """
    for i in range(len(seq) - k + 1):
        km = seq[i:i + k]
        if set(km) & AMBIG:
            continue
        yield km


def panel_weights(seqs, k, unweighted=False) -> dict:
    counts = Counter()
    for _, s in seqs:
        for km in set(kmers(s, k)):
            counts[km] += 1
    if unweighted:
        return {km: 1.0 for km in counts}
    return {km: float(c) for km, c in counts.items()}


def coverage(selected: set, weights: dict) -> float:
    tot = sum(weights.values())
    if tot == 0:
        return 0.0
    return sum(w for km, w in weights.items() if km in selected) / tot


def per_strain_coverage(selected: set, seqs, k):
    """Fraction of each strain's own k-mers the design set covers.

    This is the number that matters and it is routinely far worse than the
    aggregate. A design can reach 90% of weighted panel k-mers while leaving
    one divergent clade at 40% — which means no coverage there at all. Always
    report the minimum and the distribution, never the aggregate alone.
    """
    out = []
    for sid, s in seqs:
        ks = set(kmers(s, k))
        out.append((sid, len(ks & selected) / len(ks) if ks else float("nan")))
    return out


# ---------------------------------------------------------------- cocktail


def greedy_cocktail(seqs, k, n, weights):
    """Greedy max-coverage over natural sequences.

    Max-coverage is NP-hard; greedy carries the standard (1-1/e) ~ 63%
    approximation guarantee and in practice lands close to optimal on real
    panels. Not worth an exact solver at these sizes.
    """
    seq_kmers = {sid: set(kmers(s, k)) for sid, s in seqs}
    chosen, covered = [], set()
    for _ in range(n):
        best, best_gain = None, -1.0
        for sid, ks in seq_kmers.items():
            if sid in chosen:
                continue
            gain = sum(weights.get(km, 0.0) for km in ks - covered)
            if gain > best_gain:
                best, best_gain = sid, gain
        if best is None or best_gain <= 0:
            log(f"no further coverage gain after {len(chosen)} sequences")
            break
        chosen.append(best)
        covered |= seq_kmers[best]
        log(f"+ {best}  cumulative coverage {coverage(covered, weights):.4f}")
    return chosen, covered


# ------------------------------------------------------------------ mosaic


def build_graph(seqs, k):
    nxt, starts = defaultdict(set), Counter()
    for _, s in seqs:
        ks = list(kmers(s, k))
        if ks:
            starts[ks[0][:k - 1]] += 1
        for km in ks:
            nxt[km[:k - 1]].add(km)
    return nxt, starts


def beam_mosaic(seqs, k, weights, covered, beam_width, max_len, min_len):
    """Beam search for one mosaic maximising newly covered k-mer weight.

    Expansion only follows edges that exist in the panel, so the constructed
    sequence can never contain an unobserved k-mer. Beam search rather than
    exact optimisation because the graph is large and the objective submodular
    — a good path found quickly beats an optimal path found never. Widen
    --beam if coverage plateaus early.
    """
    nxt, starts = build_graph(seqs, k)
    if not starts:
        return "", set()

    seeds = [s for s, _ in starts.most_common(max(1, beam_width // 4))]
    beam = [(s, frozenset(), 0.0) for s in seeds]
    best = ("", frozenset(), -1.0)

    while beam:
        new_beam = []
        for seq, used, score in beam:
            tail = seq[-(k - 1):] if len(seq) >= k - 1 else seq
            cands = nxt.get(tail, ())
            if not cands or len(seq) >= max_len:
                if len(seq) >= min_len and score > best[2]:
                    best = (seq, used, score)
                continue
            for km in cands:
                nseq = seq + km[-1]
                if km in used or km in covered:
                    gain, nused = 0.0, used
                else:
                    gain, nused = weights.get(km, 0.0), used | {km}
                new_beam.append((nseq, nused, score + gain))
        if not new_beam:
            break
        new_beam.sort(key=lambda t: -t[2])
        seen, pruned = set(), []
        for item in new_beam:
            key = (item[0][-(k - 1):], len(item[1]))
            if key in seen:
                continue
            seen.add(key)
            pruned.append(item)
            if len(pruned) >= beam_width:
                break
        beam = pruned
        for item in beam:
            if len(item[0]) >= min_len and item[2] > best[2]:
                best = item

    return best[0], set(best[1])


def design_mosaics(seqs, k, n, weights, beam_width, max_len, min_len):
    covered, out = set(), []
    for i in range(n):
        seq, new = beam_mosaic(seqs, k, weights, covered, beam_width,
                               max_len, min_len)
        if not seq or not new:
            log(f"mosaic {i + 1} added no new coverage; stopping")
            break
        covered |= set(kmers(seq, k))
        out.append(seq)
        log(f"mosaic {i + 1}: {len(seq)} aa, cumulative coverage "
            f"{coverage(covered, weights):.4f}")
    return out, covered


# -------------------------------------------------------- junction checking


def novel_kmers(designs, seqs, k):
    """k-mers in the designs that occur nowhere in the panel.

    Each one is a potential junctional neo-epitope: a peptide the construct
    can present that no real strain presents. These divert the response
    without adding protection, and in a self-similar case could be worse than
    neutral — so they are flagged, never silently tolerated.
    """
    panel = set()
    for _, s in seqs:
        panel |= set(kmers(s, k))
    out = []
    for name, d in designs:
        for i in range(len(d) - k + 1):
            km = d[i:i + k]
            if (set(km) & AMBIG) or km in panel:
                continue
            out.append((name, i + 1, km))
    return out


def report(selected, seqs, k, weights, label, out_path: Path):
    agg = coverage(selected, weights)
    per = [(sid, c) for sid, c in per_strain_coverage(selected, seqs, k)
           if c == c]
    vals = sorted(c for _, c in per)
    worst = sorted(per, key=lambda t: t[1])[:10]

    with out_path.open("w") as fh:
        fh.write("strain\tcoverage\n")
        for sid, c in sorted(per, key=lambda t: t[1]):
            fh.write(f"{sid}\t{c:.4f}\n")

    log(f"--- {label} ---")
    log(f"aggregate weighted {k}-mer coverage: {100 * agg:.1f}%")
    below = 0
    if vals:
        n = len(vals)
        below = sum(1 for v in vals if v < 0.5)
        log(f"per-strain coverage: median {100 * vals[n // 2]:.1f}%, "
            f"min {100 * vals[0]:.1f}%, "
            f"10th pct {100 * vals[max(0, n // 10)]:.1f}%")
        log(f"strains below 50% coverage: {below} / {n}")
        if below:
            log("  ^ these strains are effectively uncovered. The aggregate "
                "figure hides them — decide deliberately whether to add a "
                "component, rebalance the panel, or accept the gap.")
    for sid, c in worst:
        log(f"  {100 * c:6.1f}%  {sid}")
    return {"aggregate_coverage": round(agg, 4),
            "median_per_strain": round(vals[len(vals) // 2], 4) if vals else None,
            "min_per_strain": round(vals[0], 4) if vals else None,
            "strains_below_50pct": below,
            "n_strains": len(vals)}


def main() -> None:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("mode", choices=["cocktail", "mosaic", "evaluate"])
    p.add_argument("--fasta", required=True, help="strain panel (unaligned)")
    p.add_argument("--designs", help="FASTA of designs (evaluate mode)")
    p.add_argument("-k", type=int, default=9,
                   help="k-mer length: 9 for class I, 15 for class II")
    p.add_argument("-n", type=int, default=4, help="number of components")
    p.add_argument("--unweighted", action="store_true")
    p.add_argument("--beam", type=int, default=200)
    p.add_argument("--max-len", type=int, default=0, dest="max_len")
    p.add_argument("--min-len", type=int, default=0, dest="min_len")
    p.add_argument("--out", required=True,
                   help="output prefix; writes <prefix>.fasta, "
                        "<prefix>.coverage.tsv and a .meta.json sidecar")
    args = p.parse_args()

    if args.k < 3:
        die("-k must be at least 3")

    seqs = read_fasta(Path(args.fasta))
    mean_len = sum(len(s) for _, s in seqs) / len(seqs)
    max_len = args.max_len or int(1.2 * mean_len)
    min_len = args.min_len or int(0.8 * mean_len)
    weights = panel_weights(seqs, args.k, args.unweighted)
    log(f"panel: {len(seqs)} sequences, mean length {mean_len:.0f}, "
        f"{len(weights)} distinct {args.k}-mers, "
        f"{'unweighted' if args.unweighted else 'frequency-weighted'}")

    prefix = Path(args.out)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    cov_path = Path(str(prefix) + ".coverage.tsv")
    fasta_path = Path(str(prefix) + ".fasta")

    meta = {"panel": args.fasta, "n_panel_sequences": len(seqs),
            "k": args.k, "mode": args.mode,
            "weighting": "unweighted" if args.unweighted else "frequency",
            "algorithm": "greedy max-coverage" if args.mode == "cocktail"
            else "de Bruijn beam search" if args.mode == "mosaic" else "scoring only"}

    if args.mode == "cocktail":
        chosen, covered = greedy_cocktail(seqs, args.k, args.n, weights)
        by_id = dict(seqs)
        named = [(sid, by_id[sid]) for sid in chosen]
        write_fasta(fasta_path, named)
        meta.update(report(covered, seqs, args.k, weights,
                           f"cocktail of {len(chosen)} natural sequences",
                           cov_path))
        meta["components"] = chosen
        log(f"wrote {fasta_path}")

    elif args.mode == "mosaic":
        mosaics, covered = design_mosaics(seqs, args.k, args.n, weights,
                                          args.beam, max_len, min_len)
        named = [(f"mosaic_{i + 1}", m) for i, m in enumerate(mosaics)]
        write_fasta(fasta_path, named)
        meta.update(report(covered, seqs, args.k, weights,
                           f"mosaic set of {len(mosaics)} designs", cov_path))
        nov = novel_kmers(named, seqs, args.k)
        meta["novel_kmers"] = len(nov)
        if nov:
            log(f"ERROR: {len(nov)} novel {args.k}-mers found. This is "
                f"impossible by construction — do not use these designs, and "
                f"report this as a bug.")
            for nm, pos, km in nov[:10]:
                log(f"  {nm}:{pos} {km}")
        else:
            log(f"junction check passed: every {args.k}-mer in every design "
                f"occurs in the panel")
        meta["beam_width"] = args.beam
        log(f"wrote {fasta_path}")

    else:  # evaluate
        if not args.designs:
            die("evaluate mode requires --designs")
        designs = read_fasta(Path(args.designs))
        covered = set()
        for _, d in designs:
            covered |= set(kmers(d, args.k))
        meta.update(report(covered, seqs, args.k, weights,
                           f"supplied design set ({len(designs)} sequences)",
                           cov_path))
        nov = novel_kmers(designs, seqs, args.k)
        meta["novel_kmers"] = len(nov)
        meta["designs"] = args.designs
        log(f"junction check: {len(nov)} {args.k}-mers absent from the panel")
        if nov:
            log("  Each is a potential junctional neo-epitope — a peptide the "
                "construct presents that no real strain does. Review them "
                "before proceeding; they divert the response without adding "
                "protection.")
            nov_path = Path(str(prefix) + ".novel_kmers.tsv")
            with nov_path.open("w") as fh:
                fh.write("design\tposition\tkmer\n")
                for nm, pos, km in nov:
                    fh.write(f"{nm}\t{pos}\t{km}\n")
            log(f"  wrote {nov_path}")

    write_meta(cov_path, meta)


if __name__ == "__main__":
    main()
