#!/usr/bin/env python3
"""Local sequence-level developability liabilities for protein constructs.

These are the failures that kill constructs at the bench more often than design
errors do, and all of them are visible from sequence alone — no GPU, no cloud
job, no structure. Running this before submitting anything to a folding or
design service is close to free and routinely saves a round.

Checks performed, per sequence:

  free_cysteines      Unpaired Cys count. An odd count is the strong signal;
                      free thiols drive disulfide scrambling, dimerisation and
                      heterogeneity. Resurfacing creates these by accident.
  n_glyc_sequons      N-X-S/T, X != P. Flagged as ADDED when a parent is given
                      and the sequon is absent there — surface redesign creates
                      unintended glycosylation sites routinely, and an added
                      sequon on or beside a target epitope can mask it.
  deamidation         NG, NS, NT motifs. NG is the fast one.
  isomerisation       DG, DS motifs.
  oxidation           Met and Trp counts (exposure matters, but sequence is the
                      screen).
  hydrophobic_patches Windows of high local hydrophobicity, an aggregation
                      proxy.
  gravy               Grand average of hydropathy (Kyte-Doolittle).
  pI                  Isoelectric point. Constructs with pI near the
                      purification or formulation pH aggregate; report it
                      against the buffer rather than in the abstract.
  rk_ratio            Charge composition, context for pI.

Every flag is a liability to weigh, not a verdict. A construct with three
deamidation motifs is not disqualified — it is a construct whose stability
study should look for deamidation. Report them, do not silently filter.

Writes a TSV plus a sidecar <out>.meta.json.

Usage
-----
  python developability.py --fasta constructs.fasta --out work/dev.tsv
  python developability.py --fasta constructs.fasta --parent parent.fasta \
      --out work/dev.tsv
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")

KD = {  # Kyte & Doolittle hydropathy
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5,
    "Q": -3.5, "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5,
    "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6,
    "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}

# EMBOSS pKa set — stated explicitly because pI values differ by a few tenths
# between references and a bare number invites false precision.
PKA_SIDE = {"D": 3.65, "E": 4.25, "C": 8.18, "Y": 10.07, "H": 6.00,
            "K": 10.53, "R": 12.48}
PKA_NTERM, PKA_CTERM = 8.6, 3.6


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
    out = [(n, "".join(c for c in s.upper() if c in VALID_AA))
           for n, s in records]
    out = [(n, s) for n, s in out if s]
    if not out:
        die(f"no usable sequences in {path}")
    return out


def net_charge(seq: str, ph: float) -> float:
    q = 1.0 / (1.0 + 10 ** (ph - PKA_NTERM))
    q -= 1.0 / (1.0 + 10 ** (PKA_CTERM - ph))
    for aa in "KRH":
        n = seq.count(aa)
        if n:
            q += n / (1.0 + 10 ** (ph - PKA_SIDE[aa]))
    for aa in "DECY":
        n = seq.count(aa)
        if n:
            q -= n / (1.0 + 10 ** (PKA_SIDE[aa] - ph))
    return q


def isoelectric_point(seq: str) -> float:
    lo, hi = 0.0, 14.0
    for _ in range(100):
        mid = (lo + hi) / 2
        if net_charge(seq, mid) > 0:
            lo = mid
        else:
            hi = mid
    return round((lo + hi) / 2, 2)


def motifs(seq: str, pairs: tuple[str, ...]) -> list[int]:
    return [i + 1 for i in range(len(seq) - 1) if seq[i:i + 2] in pairs]


def sequons(seq: str) -> list[int]:
    """N-X-S/T, X != P. 1-based position of the N.

    Duplicated from glycoengineering only to flag sequons created by design.
    For any glycan engineering — adding, removing, occupancy, shielding — use
    the glycoengineering skill, which owns that logic.
    """
    return [i + 1 for i in range(len(seq) - 2)
            if seq[i] == "N" and seq[i + 1] != "P" and seq[i + 2] in "ST"]


def hydrophobic_patches(seq: str, window: int = 7, cutoff: float = 2.0):
    """Windows whose mean hydropathy exceeds a cutoff, merged into runs."""
    hits = []
    for i in range(len(seq) - window + 1):
        w = seq[i:i + window]
        if sum(KD[c] for c in w) / window >= cutoff:
            hits.append(i + 1)
    runs, start, prev = [], None, None
    for h in hits:
        if start is None:
            start = prev = h
        elif h == prev + 1:
            prev = h
        else:
            runs.append((start, prev + window - 1))
            start = prev = h
    if start is not None:
        runs.append((start, prev + window - 1))
    return runs


def main() -> None:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--fasta", required=True, help="construct sequences")
    p.add_argument("--parent",
                   help="parent/reference FASTA; enables ADDED-sequon and "
                        "added-cysteine detection")
    p.add_argument("--ph", type=float, default=7.4,
                   help="pH to report net charge at (default 7.4)")
    p.add_argument("--patch-window", type=int, default=7, dest="patch_window")
    p.add_argument("--patch-cutoff", type=float, default=2.0,
                   dest="patch_cutoff")
    p.add_argument("--out", required=True)
    args = p.parse_args()

    seqs = read_fasta(Path(args.fasta))
    parent_seq = None
    parent_sequon_pos: set[int] = set()
    parent_sequon_count = parent_cys = None
    if args.parent:
        par = read_fasta(Path(args.parent))
        if len(par) > 1:
            log(f"parent file has {len(par)} sequences; using the first "
                f"({par[0][0]}) as the comparison reference")
        parent_seq = par[0][1]
        parent_sequon_pos = set(sequons(parent_seq))
        parent_sequon_count = len(parent_sequon_pos)
        parent_cys = parent_seq.count("C")
        log(f"parent {par[0][0]}: {len(parent_seq)} aa, "
            f"{parent_sequon_count} sequons, {parent_cys} cysteines")

    rows = []
    for name, s in seqs:
        sq = sequons(s)
        # Attribute added sequons by POSITION when the construct is the same
        # length as the parent (the substitution-only case, which covers
        # resurfacing and point stabilisation). Comparing motif strings
        # instead would miss a new NGT in a protein that already has an NGT
        # somewhere else — a silent false negative on exactly the check that
        # matters. When lengths differ, an indel has shifted numbering, so
        # report the count delta and say positions need an alignment.
        added, added_note = [], ""
        if args.parent:
            if len(s) == len(parent_seq):
                added = [i for i in sq if i not in parent_sequon_pos]
            else:
                added_note = "length-differs"
        cys = s.count("C")
        patches = hydrophobic_patches(s, args.patch_window, args.patch_cutoff)
        deam = motifs(s, ("NG", "NS", "NT"))
        iso = motifs(s, ("DG", "DS"))
        rows.append({
            "construct": name,
            "length": len(s),
            "cysteines": cys,
            "cys_odd": int(cys % 2 == 1),
            "cys_added_vs_parent": (cys - parent_cys) if parent_cys is not None else "",
            "n_glyc_sequons": len(sq),
            "sequon_positions": ",".join(map(str, sq)),
            "sequon_count_delta": (len(sq) - parent_sequon_count)
            if args.parent else "",
            "sequons_added": (added_note or len(added)) if args.parent else "",
            "added_sequon_positions": (added_note or ",".join(map(str, added)))
            if args.parent else "",
            "deamidation_motifs": len(deam),
            "deamidation_positions": ",".join(map(str, deam)),
            "isomerisation_motifs": len(iso),
            "met_trp": s.count("M") + s.count("W"),
            "hydrophobic_patches": len(patches),
            "patch_spans": ";".join(f"{a}-{b}" for a, b in patches),
            "gravy": round(sum(KD[c] for c in s) / len(s), 3),
            "pI": isoelectric_point(s),
            f"net_charge_pH{args.ph}": round(net_charge(s, args.ph), 2),
        })

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    cols = list(rows[0].keys())
    with out.open("w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")

    n_odd = sum(r["cys_odd"] for r in rows)
    n_added = sum(1 for r in rows
                  if isinstance(r["sequons_added"], int) and r["sequons_added"])
    n_shift = sum(1 for r in rows if r["sequons_added"] == "length-differs")
    n_patch = sum(1 for r in rows if r["hydrophobic_patches"] > 0)
    log(f"{len(rows)} constructs -> {out}")
    if n_odd:
        log(f"  {n_odd} with an ODD cysteine count — unpaired thiol, expect "
            f"scrambling and heterogeneity unless deliberate")
    if n_added:
        log(f"  {n_added} with N-glycosylation sequons ADDED relative to the "
            f"parent — check none sits on or beside a target epitope")
    if n_shift:
        log(f"  {n_shift} differ in length from the parent — sequon positions "
            f"not attributable without an alignment; see sequon_count_delta")
    if n_patch:
        log(f"  {n_patch} with hydrophobic patches (aggregation proxy)")
    near_pi = [r["construct"] for r in rows if abs(r["pI"] - args.ph) < 1.0]
    if near_pi:
        log(f"  {len(near_pi)} with pI within 1 unit of pH {args.ph} — "
            f"solubility minimum, check against your actual buffer")

    write_meta(out, {
        "fasta": args.fasta,
        "parent": args.parent,
        "n_constructs": len(rows),
        "ph": args.ph,
        "hydropathy_scale": "Kyte-Doolittle",
        "pka_set": "EMBOSS",
        "patch_window": args.patch_window,
        "patch_cutoff": args.patch_cutoff,
        "note": "sequence-level liabilities only; no structural or "
                "experimental verdict is implied",
    })


if __name__ == "__main__":
    main()
