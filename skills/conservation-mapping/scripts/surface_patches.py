#!/usr/bin/env python3
"""Find spatially contiguous conserved surface patches on a structure.

Conservation alone is not a target. Buried residues are conserved because the
fold requires it, and they are useless as antibody epitopes. What matters for
immunogen design is conserved residues that are *exposed* and *adjacent* —
a patch large enough for an antibody footprint, not a scattered list of
conserved surface positions.

This takes the conservation table from conservation_score.py plus the
B-factor-mapped PDB it wrote, computes relative solvent accessibility, and
clusters low-entropy exposed residues by CA distance into patches.

Patch area is reported because that is what decides whether a patch is a
plausible epitope: a typical Fab buries roughly 600-900 A^2. A residue count
alone does not tell you that.

Requires freesasa (conda-forge) and biopython. If DSSP is preferred and
available, pass --rsa-tsv with your own two-column table instead.

Writes a sidecar <out>.meta.json.

Usage
-----
  python surface_patches.py --pdb conservation.pdb --chain A \
      --conservation conservation.tsv --out patches.tsv

Note on glycans: for a viral envelope glycoprotein an exposed conserved patch
may still be antibody-inaccessible under the glycan shield. Run the
glycoengineering skill's GlycoSHIELD workflow and read its shielding B-factors
alongside these patches. This script sees geometry, not glycans.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

# Tien et al. (2013) theoretical maximum ASA, the usual choice for RSA.
MAX_ASA = {
    "A": 129.0, "R": 274.0, "N": 195.0, "D": 193.0, "C": 167.0,
    "Q": 225.0, "E": 223.0, "G": 104.0, "H": 224.0, "I": 197.0,
    "L": 201.0, "K": 236.0, "M": 224.0, "F": 240.0, "P": 159.0,
    "S": 155.0, "T": 172.0, "W": 285.0, "Y": 263.0, "V": 174.0,
}

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "MSE": "M",
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


def read_ca(path: Path, chain: str | None):
    """[(resseq, icode, one_letter, (x,y,z))] for CA atoms."""
    out, seen = [], set()
    with path.open() as fh:
        for line in fh:
            if not line.startswith("ATOM") or line[12:16].strip() != "CA":
                continue
            ch = line[21]
            if chain and ch != chain:
                continue
            one = THREE_TO_ONE.get(line[17:20].strip().upper())
            if one is None:
                continue
            try:
                resseq = int(line[22:26])
                xyz = (float(line[30:38]), float(line[38:46]),
                       float(line[46:54]))
            except ValueError:
                continue
            key = (resseq, line[26].strip())
            if key in seen:
                continue
            seen.add(key)
            out.append((resseq, line[26].strip(), one, xyz))
    return out


def compute_rsa(pdb: Path, chain: str | None) -> dict:
    """{(resseq, icode): rsa} via freesasa."""
    try:
        import freesasa
    except ImportError:
        die("freesasa required (conda install -c conda-forge freesasa), "
            "or supply --rsa-tsv")
    structure = freesasa.Structure(str(pdb))
    result = freesasa.calc(structure)
    areas = {}
    for i in range(structure.nAtoms()):
        ch = structure.chainLabel(i)
        if chain and ch != chain:
            continue
        resn = structure.residueName(i).strip().upper()
        one = THREE_TO_ONE.get(resn)
        if one is None:
            continue
        raw = structure.residueNumber(i).strip()
        icode = raw[-1] if raw and raw[-1].isalpha() else ""
        try:
            resseq = int(raw[:-1] if icode else raw)
        except ValueError:
            continue
        key = (resseq, icode)
        areas.setdefault(key, [0.0, one])
        areas[key][0] += result.atomArea(i)
    return {k: min(1.0, v[0] / MAX_ASA.get(v[1], 200.0))
            for k, v in areas.items()}


def read_rsa_tsv(path: Path) -> dict:
    out = {}
    with path.open() as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            try:
                out[(int(parts[0]), "")] = float(parts[1])
            except ValueError:
                continue
    return out


def read_conservation(path: Path) -> tuple[dict, str]:
    """{(resseq, icode): entropy} from conservation_score.py output.

    Joins on the `pdb_resnum` column when present — that column holds the
    structure's own residue IDs as reconciled by alignment, which is the only
    reliable key. Falls back to `ref_pos` only when the table was written
    without --pdb, and says which key it used so a silent mis-join is
    impossible.
    """
    out = {}
    with path.open() as fh:
        header = fh.readline().rstrip("\n").split("\t")
        if "entropy" not in header:
            die("conservation table has no entropy column — was it written "
                "by conservation_score.py?")
        ie = header.index("entropy")
        if "pdb_resnum" in header:
            ik, key_name = header.index("pdb_resnum"), "pdb_resnum"
        elif "ref_pos" in header:
            ik, key_name = header.index("ref_pos"), "ref_pos"
        else:
            die("conservation table has neither pdb_resnum nor ref_pos")
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= max(ik, ie) or not parts[ik] or not parts[ie]:
                continue
            raw = parts[ik].strip()
            icode = raw[-1] if raw and raw[-1].isalpha() else ""
            try:
                resseq = int(raw[:-1] if icode else raw)
            except ValueError:
                continue
            out[(resseq, icode)] = float(parts[ie])
    return out, key_name


def cluster(residues, coords, linkage):
    """Single-linkage clustering of residues by CA distance."""
    import math
    parent = {r: r for r in residues}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    rl = list(residues)
    for i in range(len(rl)):
        for j in range(i + 1, len(rl)):
            xi, xj = coords[rl[i]], coords[rl[j]]
            d = math.dist(xi, xj)
            if d <= linkage:
                union(rl[i], rl[j])

    groups = {}
    for r in rl:
        groups.setdefault(find(r), []).append(r)
    return list(groups.values())


def main() -> None:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--pdb", required=True,
                   help="structure, ideally the B-factor-mapped PDB from "
                        "conservation_score.py")
    p.add_argument("--chain", default=None)
    p.add_argument("--conservation", required=True,
                   help="conservation.tsv from conservation_score.py")
    p.add_argument("--rsa-tsv", dest="rsa_tsv",
                   help="optional two-column resseq<TAB>rsa table (e.g. from "
                        "DSSP) in place of freesasa")
    p.add_argument("--entropy-max", type=float, default=0.2, dest="entropy_max",
                   help="conserved if entropy <= this (default 0.2)")
    p.add_argument("--rsa-min", type=float, default=0.25, dest="rsa_min",
                   help="exposed if RSA >= this (default 0.25)")
    p.add_argument("--linkage", type=float, default=9.0,
                   help="CA-CA single-linkage distance in angstroms "
                        "(default 9.0)")
    p.add_argument("--min-residues", type=int, default=4, dest="min_residues")
    p.add_argument("--out", required=True)
    args = p.parse_args()

    pdb = Path(args.pdb)
    ca = read_ca(pdb, args.chain)
    if not ca:
        die("no CA atoms found (wrong --chain?)")

    rsa = (read_rsa_tsv(Path(args.rsa_tsv)) if args.rsa_tsv
           else compute_rsa(pdb, args.chain))
    cons, key_name = read_conservation(Path(args.conservation))
    log(f"joining conservation on '{key_name}'")
    if key_name == "ref_pos":
        log("NOTE: the conservation table carries no pdb_resnum column, so "
            "this join assumes reference numbering equals PDB numbering. "
            "Re-run conservation_score.py with --pdb to reconcile them "
            "properly if the match rate below looks wrong.")

    coords, keep = {}, []
    matched = 0
    for resseq, icode, one, xyz in ca:
        key = (resseq, icode)
        r = rsa.get(key)
        e = cons.get(key)
        if e is not None:
            matched += 1
        if r is None or e is None:
            continue
        coords[key] = xyz
        if e <= args.entropy_max and r >= args.rsa_min:
            keep.append(key)

    log(f"residues: {len(ca)} in structure, {matched} matched to conservation "
        f"table, {len(keep)} conserved+exposed")
    if matched == 0:
        die("no residues matched the conservation table. Run "
            "conservation_score.py with --pdb and use the PDB it wrote, so "
            "numbering is already reconciled.")

    clusters = [c for c in cluster(keep, coords, args.linkage)
                if len(c) >= args.min_residues]
    clusters.sort(key=len, reverse=True)

    rows = []
    for i, cl in enumerate(sorted(clusters, key=len, reverse=True), 1):
        cl_sorted = sorted(cl)
        area = sum(rsa.get(k, 0.0) * MAX_ASA.get(
            next((o for rs, ic, o, _ in ca if (rs, ic) == k), "A"), 200.0)
            for k in cl_sorted)
        ents = [cons[k] for k in cl_sorted if k in cons]
        rows.append([
            f"patch_{i}", len(cl_sorted), round(area, 1),
            round(sum(ents) / len(ents), 4) if ents else "",
            ",".join(str(k[0]) + k[1] for k in cl_sorted),
        ])

    out = Path(args.out)
    with out.open("w") as fh:
        fh.write("patch_id\tn_residues\tapprox_area_A2\tmean_entropy\t"
                 "residues\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")

    log(f"patches: {len(rows)} found -> {out}")
    big = [r for r in rows if isinstance(r[2], float) and r[2] >= 600]
    if big:
        log(f"{len(big)} patch(es) at or above 600 A^2, the rough lower bound "
            f"for a Fab footprint")
    else:
        log("no patch reaches 600 A^2 — the conserved surface here may be too "
            "fragmented for a single antibody footprint. Consider relaxing "
            "--entropy-max or widening --linkage, and say which you changed.")

    write_meta(out, {
        "pdb": args.pdb,
        "chain": args.chain,
        "conservation_table": args.conservation,
        "rsa_source": "dssp-tsv" if args.rsa_tsv else "freesasa",
        "max_asa_reference": "Tien et al. 2013 theoretical",
        "entropy_max": args.entropy_max,
        "rsa_min": args.rsa_min,
        "linkage_angstrom": args.linkage,
        "min_residues": args.min_residues,
        "n_patches": len(rows),
        "residues_matched": matched,
    })


if __name__ == "__main__":
    main()
