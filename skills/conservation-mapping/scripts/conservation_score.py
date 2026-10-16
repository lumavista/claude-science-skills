#!/usr/bin/env python3
"""Per-position conservation across a strain panel, in reference numbering.

Computes Shannon entropy with Henikoff position-based sequence weighting, maps
every alignment column onto a reference sequence's residue numbering, and
optionally writes the score into a PDB B-factor column.

Why weighting is the default: strain panels reflect surveillance effort, not
pathogen diversity. A panel with 200 near-identical clade-A sequences and 8
clade-C sequences reports clade-A-specific residues as "conserved" unless
sequences are downweighted by redundancy. Henikoff weighting corrects this, and
the effective sample size printed at the top is the diagnostic that tells you
whether the panel supports any conservation claim at all.

PDB numbering is resolved by aligning the structure's own SEQRES-derived
sequence to the reference, not by assuming FASTA offsets equal residue IDs.
They routinely differ — construct tags, disordered gaps, author numbering
schemes — and that mismatch is a silent, invisible error in every downstream
structure-based step.

Writes a sidecar <out>.meta.json recording the panel, weighting and versions.

Usage
-----
  python conservation_score.py --aln panel.afa --ref STRAIN_01 \
      --out conservation.tsv

  python conservation_score.py --aln panel.afa --ref STRAIN_01 \
      --out conservation.tsv \
      --pdb antigen.pdb --pdb-chain A --pdb-out conservation.pdb

Input must be an aligned FASTA. Produce one with the phylogenetics skill
(scripts/phylogenetic_analysis.py) rather than aligning here — alignment
belongs there and this script will refuse unaligned input.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from collections import Counter
from pathlib import Path

try:
    import numpy as np
except ImportError:  # pragma: no cover
    print("error: numpy required", file=sys.stderr)
    raise SystemExit(1)

VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")
GAPS = set("-.~")

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "MSE": "M", "SEC": "U", "PYL": "O",
}

# Grouped classes for the class-entropy column. A position that varies only
# within a class is far less likely to matter structurally than one that swaps
# charge or size, and the two columns together distinguish those cases.
AA_CLASS = {
    **{c: "hydrophobic" for c in "AVLIMC"},
    **{c: "aromatic" for c in "FWY"},
    **{c: "polar" for c in "STNQ"},
    **{c: "positive" for c in "KRH"},
    **{c: "negative" for c in "DE"},
    **{c: "special" for c in "GP"},
}


def log(msg: str) -> None:
    print(msg, file=sys.stderr)


def die(msg: str, code: int = 1) -> None:
    print(f"error: {msg}", file=sys.stderr)
    raise SystemExit(code)


def write_meta(out: Path, payload: dict) -> None:
    out.with_suffix(out.suffix + ".meta.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n"
    )


# --------------------------------------------------------------------------- IO


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
    return records


def read_alignment(path: Path) -> list[tuple[str, str]]:
    recs = read_fasta(path)
    if not recs:
        die(f"no sequences found in {path}")
    lengths = {len(s) for _, s in recs}
    if len(lengths) != 1:
        die(
            f"sequences are not aligned ({len(lengths)} distinct lengths, "
            f"e.g. {sorted(lengths)[:4]}). Align first with the phylogenetics "
            f"skill: scripts/phylogenetic_analysis.py <in.fasta> --aligner famsa"
        )
    return recs


def pdb_sequence(path: Path, chain: str | None) -> list[tuple[int, str, str]]:
    """Return [(resseq, icode, one_letter)] in file order, CA atoms only."""
    out, seen = [], set()
    with path.open() as fh:
        for line in fh:
            if not line.startswith("ATOM"):
                continue
            if line[12:16].strip() != "CA":
                continue
            ch = line[21]
            if chain and ch != chain:
                continue
            resn = line[17:20].strip().upper()
            one = THREE_TO_ONE.get(resn)
            if one is None:
                continue
            try:
                resseq = int(line[22:26])
            except ValueError:
                continue
            icode = line[26].strip()
            key = (ch, resseq, icode)
            if key in seen:
                continue
            seen.add(key)
            out.append((resseq, icode, one))
    return out


# --------------------------------------------------------------------- scoring


def henikoff_weights(seqs: list[str]) -> np.ndarray:
    """Position-based sequence weights (Henikoff & Henikoff, JMB 1994).

    Each sequence's weight is the mean over columns of 1/(r*s), where r is the
    number of distinct residue types in the column and s the count of this
    sequence's residue there. Rare residues in diverse columns contribute most,
    which is the redundancy correction a skewed strain panel needs.
    """
    n, L = len(seqs), len(seqs[0])
    w = np.zeros(n)
    for j in range(L):
        col = [s[j] for s in seqs]
        counts = Counter(c for c in col if c not in GAPS)
        r = len(counts)
        if r <= 1:
            continue
        for i, c in enumerate(col):
            if c not in GAPS:
                w[i] += 1.0 / (r * counts[c])
    if w.sum() == 0:
        return np.ones(n) / n
    return w / w.sum()


def column_stats(col, weights, gap_cutoff):
    """(entropy, class_entropy, gap_frac, consensus, consensus_frac, n_variants).

    Entropy is normalized by log2(20) to [0,1] so it is comparable across
    datasets. Gapped positions are excluded from the residue distribution but
    reported separately: a column that is 60% gap is an indel region, not a
    conserved one, and the gap fraction is what tells you that.
    """
    wsum = weights.sum()
    gap_w = sum(w for c, w in zip(col, weights) if c in GAPS)
    gap_frac = gap_w / wsum if wsum else 1.0

    dist, cls_dist = {}, {}
    for c, w in zip(col, weights):
        if c in GAPS:
            continue
        c = c.upper()
        if c not in VALID_AA:
            continue  # X/B/Z excluded, never imputed
        dist[c] = dist.get(c, 0.0) + w
        k = AA_CLASS.get(c, "other")
        cls_dist[k] = cls_dist.get(k, 0.0) + w

    tot = sum(dist.values())
    if tot <= 0 or gap_frac > gap_cutoff:
        return (float("nan"), float("nan"), gap_frac, "-", 0.0, 0)

    ent = -sum((p / tot) * math.log2(p / tot) for p in dist.values() if p > 0)
    ent = max(0.0, ent / math.log2(20))

    ctot = sum(cls_dist.values())
    cent = -sum((p / ctot) * math.log2(p / ctot)
                for p in cls_dist.values() if p > 0)
    cent = max(0.0, cent / math.log2(6))

    cons = max(dist, key=dist.get)
    return (ent, cent, gap_frac, cons, dist[cons] / tot, len(dist))


def nglyc_sequons(seq: str) -> set[int]:
    """N-X-S/T sequon positions (1-based N), X != P.

    Sequence-level only, and deliberately duplicated from the glycoengineering
    skill's find_n_glycosylation_sequons only to annotate this table inline.
    For any actual glycan work — occupancy, shielding, site addition or removal
    — use the glycoengineering skill, which is the owner of that logic.
    """
    hits = set()
    s = seq.upper()
    for i in range(len(s) - 2):
        if s[i] == "N" and s[i + 1] != "P" and s[i + 2] in "ST":
            hits.add(i + 1)
    return hits


# ------------------------------------------------------------- PDB numbering


def map_reference_to_pdb(ref_ungapped: str, pdb_res: list[tuple[int, str, str]]):
    """Align the reference sequence to the structure's own residues.

    Returns {ref_pos (1-based) -> (resseq, icode)} plus a short report. This
    exists because FASTA offsets are not PDB residue IDs: construct tags,
    disordered gaps and author numbering all shift them, and assuming they
    match is a silent error that corrupts every structure-based step
    downstream without ever raising.
    """
    try:
        from Bio import Align
    except ImportError:
        die("biopython required for --pdb (pip install biopython)")

    pdb_seq = "".join(one for _, _, one in pdb_res)
    if not pdb_seq:
        die("no CA atoms with standard residues found in the PDB "
            "(wrong --pdb-chain?)")

    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.open_gap_score = -11
    aligner.extend_gap_score = -1
    aligner.substitution_matrix = None
    aligner.match_score = 2
    aligner.mismatch_score = -1

    aln = aligner.align(ref_ungapped, pdb_seq)[0]
    mapping, n_match, n_mismatch = {}, 0, 0
    for (rs, re_), (ps, pe) in zip(aln.aligned[0], aln.aligned[1]):
        for off in range(re_ - rs):
            ref_i = rs + off          # 0-based into ref_ungapped
            pdb_i = ps + off          # 0-based into pdb_seq
            resseq, icode, one = pdb_res[pdb_i]
            mapping[ref_i + 1] = (resseq, icode)
            if ref_ungapped[ref_i] == one:
                n_match += 1
            else:
                n_mismatch += 1

    ident = n_match / len(pdb_seq) if pdb_seq else 0.0
    report = {
        "pdb_residues": len(pdb_seq),
        "reference_length": len(ref_ungapped),
        "aligned_positions": len(mapping),
        "identical": n_match,
        "mismatched": n_mismatch,
        "identity_over_pdb": round(ident, 4),
    }
    return mapping, report


def write_pdb(pdb_in, pdb_out, chain, score_by_refpos, ref_to_pdb):
    """Write 100*(1-entropy) into B-factors: high B = conserved.

    Inverted on purpose, so that the usual blue-to-red B-factor coloring puts
    conserved surface in the hot color — which is the thing you want to look
    at. Unscored residues get -1 so they are visibly distinct from genuinely
    variable ones rather than rendering as a misleading zero.

    Note the glycoengineering skill's GlycoSHIELD workflow also encodes
    shielding in B-factors. Keep the two PDBs in separate files and say which
    is which; they are not comparable scales.
    """
    by_res = {}
    for refpos, val in score_by_refpos.items():
        key = ref_to_pdb.get(refpos)
        if key is not None:
            by_res[key] = val

    n_set = n_miss = 0
    with open(pdb_in) as fi, open(pdb_out, "w") as fo:
        for line in fi:
            if line.startswith(("ATOM", "HETATM")):
                ch = line[21]
                if chain and ch != chain:
                    fo.write(line)
                    continue
                try:
                    resseq = int(line[22:26])
                except ValueError:
                    fo.write(line)
                    continue
                key = (resseq, line[26].strip())
                if key in by_res:
                    val, n_set = by_res[key], n_set + 1
                else:
                    val, n_miss = -1.0, n_miss + 1
                fo.write(line[:60] + f"{val:6.2f}" + line[66:])
            else:
                fo.write(line)
    log(f"pdb: wrote {pdb_out} ({n_set} atoms scored, {n_miss} unscored)")


# ------------------------------------------------------------------------ main


def main() -> None:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--aln", required=True, help="aligned FASTA (strain panel)")
    p.add_argument("--ref", help="reference sequence ID (default: first)")
    p.add_argument("--out", required=True)
    p.add_argument("--no-weight", action="store_true",
                   help="disable Henikoff weighting (not recommended for "
                        "unevenly sampled panels)")
    p.add_argument("--gap-cutoff", type=float, default=0.5,
                   help="columns above this weighted gap fraction are reported "
                        "but not scored (default 0.5)")
    p.add_argument("--pdb", help="PDB to receive conservation in B-factors")
    p.add_argument("--pdb-chain", default=None, dest="pdb_chain")
    p.add_argument("--pdb-out", default=None, dest="pdb_out")
    args = p.parse_args()

    recs = read_alignment(Path(args.aln))
    ids = [i for i, _ in recs]
    seqs = [s for _, s in recs]

    if args.ref:
        if args.ref not in ids:
            die(f"reference {args.ref!r} not in alignment. First IDs: {ids[:5]}")
        ridx = ids.index(args.ref)
    else:
        ridx = 0
        log(f"no --ref given, using {ids[0]!r}")

    ref_aln = seqs[ridx]
    ref_ungapped = "".join(c for c in ref_aln if c not in GAPS)
    sequons = nglyc_sequons(ref_ungapped)

    colmap, n = {}, 0
    for j, c in enumerate(ref_aln):
        if c not in GAPS:
            n += 1
            colmap[j] = n

    weights = (np.ones(len(seqs)) / len(seqs) if args.no_weight
               else henikoff_weights(seqs))
    eff_n = float(1.0 / np.sum(weights ** 2))

    log(f"panel: {len(seqs)} sequences, {len(ref_aln)} columns, "
        f"effective N = {eff_n:.1f}")
    if eff_n < 10:
        log("WARNING: effective sample size below 10. Conservation estimates "
            "are noisy and clade-biased; report this figure alongside any "
            "conservation claim rather than the raw sequence count.")

    rows = []
    for j in range(len(ref_aln)):
        col = [s[j] for s in seqs]
        ent, cent, gapf, cons, consf, nvar = column_stats(
            col, weights, args.gap_cutoff)
        refnum = colmap.get(j)
        refaa = ref_aln[j] if ref_aln[j] not in GAPS else "-"
        rows.append({
            "aln_col": j + 1,
            "ref_pos": refnum if refnum else "",
            "ref_aa": refaa,
            "entropy": "" if ent != ent else f"{ent:.4f}",
            "class_entropy": "" if cent != cent else f"{cent:.4f}",
            "gap_frac": f"{gapf:.3f}",
            "consensus_aa": cons,
            "consensus_frac": f"{consf:.3f}",
            "n_variants": nvar,
            "ref_differs_from_consensus": int(
                refaa not in GAPS and cons != "-" and refaa.upper() != cons),
            "nglyc_sequon_N": int(refnum in sequons) if refnum else 0,
        })

    out = Path(args.out)

    # Resolve PDB numbering BEFORE writing the table, so the reconciled
    # residue IDs are persisted as a column. Without this the mapping lives
    # only inside the B-factors and every downstream join has to rediscover
    # it — which is the exact trap glycoshield.md warns about.
    ref_to_pdb, pdb_report = {}, None
    if args.pdb:
        pdb_res = pdb_sequence(Path(args.pdb), args.pdb_chain)
        ref_to_pdb, pdb_report = map_reference_to_pdb(ref_ungapped, pdb_res)
        log(f"pdb numbering: {pdb_report['aligned_positions']} positions "
            f"mapped, identity over structure "
            f"{pdb_report['identity_over_pdb']:.3f}")
        if pdb_report["identity_over_pdb"] < 0.8:
            log("WARNING: structure matches the reference poorly. Check that "
                "the PDB is the same protein and the right chain before "
                "trusting the mapping.")
        for r in rows:
            key = ref_to_pdb.get(int(r["ref_pos"])) if r["ref_pos"] else None
            r["pdb_resnum"] = f"{key[0]}{key[1]}" if key else ""

    cols = list(rows[0].keys())
    with out.open("w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")

    scored = [float(r["entropy"]) for r in rows if r["entropy"]]
    stats = {}
    if scored:
        arr = np.array(scored)
        stats = {
            "median_entropy": round(float(np.median(arr)), 4),
            "near_invariant_columns": int((arr < 0.1).sum()),
            "highly_variable_columns": int((arr > 0.5).sum()),
            "scored_columns": len(scored),
        }
        log(f"entropy: median {stats['median_entropy']:.3f}, "
            f"{stats['near_invariant_columns']} near-invariant, "
            f"{stats['highly_variable_columns']} highly variable")
    n_diff = sum(r["ref_differs_from_consensus"] for r in rows)
    log(f"reference differs from consensus at {n_diff} positions -> {out}")

    if args.pdb:
        pdb_out = args.pdb_out or str(out.with_suffix(".pdb"))
        score_by_refpos = {
            int(r["ref_pos"]): 100.0 * (1.0 - float(r["entropy"]))
            for r in rows if r["ref_pos"] and r["entropy"]
        }
        write_pdb(args.pdb, pdb_out, args.pdb_chain, score_by_refpos, ref_to_pdb)

    write_meta(out, {
        "alignment": args.aln,
        "reference": ids[ridx],
        "n_sequences": len(seqs),
        "effective_n": round(eff_n, 3),
        "weighting": "none" if args.no_weight else "henikoff",
        "gap_cutoff": args.gap_cutoff,
        "entropy_normalisation": "log2(20), range 0-1",
        "pdb_mapping": pdb_report,
        **stats,
    })


if __name__ == "__main__":
    main()
