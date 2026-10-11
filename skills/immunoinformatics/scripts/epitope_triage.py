#!/usr/bin/env python3
"""Epitope triage helpers for protective vaccine design.

Subcommands:
  prepare       Clean and validate an antigen FASTA.
  mhc1          MHC class I binding prediction (MHCflurry, or NetMHCpan if present).
  mhc2          MHC class II binding prediction (NetMHCIIpan; reports if absent).
  iedb          Retrieve experimentally validated epitopes from IEDB.
  conservation  Score candidates against a set of circulating variant sequences.
  selfcheck     Flag candidates with close identity to the human proteome.
  coverage      HLA population coverage via the IEDB tool.

Every subcommand writes a TSV and prints a one-line summary to stderr.
Tool versions are recorded in a sidecar <out>.meta.json for reproducibility.
"""

from __future__ import annotations

import argparse
import json
import re
import shutil
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")


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
                chunks.append(line.upper())
    if name is not None:
        records.append((name, "".join(chunks)))
    if not records:
        die(f"no sequences found in {path}")
    return records


def write_tsv(path: Path, header: list[str], rows: list[list]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as fh:
        fh.write("\t".join(header) + "\n")
        for row in rows:
            fh.write("\t".join("" if v is None else str(v) for v in row) + "\n")


def read_tsv(path: Path) -> tuple[list[str], list[dict]]:
    lines = path.read_text().splitlines()
    if not lines:
        die(f"{path} is empty")
    header = lines[0].split("\t")
    rows = [dict(zip(header, ln.split("\t"))) for ln in lines[1:] if ln.strip()]
    return header, rows


def peptide_column(rows: list[dict]) -> str:
    for key in ("peptide", "sequence", "Peptide"):
        if rows and key in rows[0]:
            return key
    die("input has no 'peptide' column; run mhc1/mhc2 first")
    return ""  # unreachable


# ----------------------------------------------------------------------- prepare


def cmd_prepare(args: argparse.Namespace) -> None:
    records = read_fasta(Path(args.fasta))
    kept, problems = [], []
    for name, seq in records:
        bad = sorted(set(seq) - VALID_AA)
        if bad:
            cleaned = "".join(c for c in seq if c in VALID_AA)
            problems.append((name, f"removed non-standard residues: {','.join(bad)}"))
            seq = cleaned
        if len(seq) < 9:
            problems.append((name, f"too short after cleaning ({len(seq)} aa), dropped"))
            continue
        if seq.startswith("M") and len(seq) > 40:
            hydro = sum(c in "AILMFWVC" for c in seq[:25]) / 25
            if hydro > 0.5:
                problems.append(
                    (name, "N-terminus looks like a signal peptide; confirm the "
                           "mature-protein boundary before predicting")
                )
        kept.append((name, seq))

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w") as fh:
        for name, seq in kept:
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i : i + 60] + "\n")

    for name, note in problems:
        log(f"  {name}: {note}")
    log(f"prepare: {len(kept)}/{len(records)} sequences written to {out}")
    write_meta(out, {"input": args.fasta, "kept": len(kept), "notes": problems})


# -------------------------------------------------------------------------- mhc1


def enumerate_peptides(records, lengths):
    seen = {}
    for name, seq in records:
        for k in lengths:
            for i in range(len(seq) - k + 1):
                pep = seq[i : i + k]
                seen.setdefault(pep, (name, i + 1))
    return seen


def cmd_mhc1(args: argparse.Namespace) -> None:
    records = read_fasta(Path(args.fasta))
    lengths = [int(x) for x in args.lengths.split(",")]
    alleles = [a.strip() for a in args.alleles.split(",") if a.strip()]
    peptides = enumerate_peptides(records, lengths)
    log(f"mhc1: {len(peptides)} unique peptides x {len(alleles)} alleles")

    try:
        from mhcflurry import Class1PresentationPredictor
    except ImportError:
        die("mhcflurry not installed. pip install 'mhcflurry>=2.1' && "
            "mhcflurry-downloads fetch models_class1_presentation")

    predictor = Class1PresentationPredictor.load()
    supported = set(predictor.supported_alleles)
    usable = [a for a in alleles if a in supported]
    for a in sorted(set(alleles) - supported):
        log(f"  warning: {a} not supported by this MHCflurry release, skipped")
    if not usable:
        die("no requested alleles are supported")

    # mhcflurry >=2.2 treats a LIST of alleles as one individual's genotype
    # (max 6, collapsed to the best allele per peptide). For epitope triage we
    # want a score per allele, so pass a dict with one sample per allele.
    allele_samples = {a: [a] for a in usable}

    pep_list = sorted(peptides)
    df = predictor.predict(peptides=pep_list, alleles=allele_samples, verbose=0)

    def num(rec, key):
        try:
            return float(rec.get(key))
        except (TypeError, ValueError):
            return float("nan")

    rows = []
    for rec in df.to_dict("records"):
        pep = rec["peptide"]
        if pep not in peptides:
            continue
        source, pos = peptides[pep]
        allele = rec.get("best_allele") or rec.get("sample_name")
        score = num(rec, "presentation_score")
        rows.append([
            pep, len(pep), source, pos, allele,
            round(score, 4),
            round(num(rec, "presentation_percentile"), 3),
            round(num(rec, "affinity"), 1),
        ])
    rows.sort(key=lambda r: (r[5] != r[5], -r[5]))  # NaNs last
    if args.top:
        rows = rows[: args.top]

    out = Path(args.out)
    write_tsv(out, [
        "peptide", "length", "source", "position", "allele",
        "presentation_score", "presentation_percentile", "affinity_nM",
    ], rows)
    log(f"mhc1: {len(rows)} rows -> {out}  (rank on percentile, not affinity)")
    import mhcflurry
    write_meta(out, {"predictor": "mhcflurry", "version": mhcflurry.__version__,
                     "alleles": usable, "lengths": lengths})


# -------------------------------------------------------------------------- mhc2


def cmd_mhc2(args: argparse.Namespace) -> None:
    exe = shutil.which("netMHCIIpan")
    if not exe:
        die("netMHCIIpan not found on PATH.\n"
            "  Class II prediction has no permissively licensed substitute.\n"
            "  Register and download from https://services.healthtech.dtu.dk/\n"
            "  (free for academic use; commercial use needs a license).")

    alleles = [a.strip() for a in args.alleles.split(",") if a.strip()]
    fasta = Path(args.fasta).resolve()
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)

    raw = out.with_suffix(".netmhciipan.txt")
    cmd = [exe, "-f", str(fasta), "-a", ",".join(alleles),
           "-length", str(args.length), "-BA"]
    log(f"mhc2: running {' '.join(cmd[:1])} on {len(alleles)} alleles")
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        die(f"netMHCIIpan failed:\n{proc.stderr[:2000]}")
    raw.write_text(proc.stdout)

    rows = []
    for line in proc.stdout.splitlines():
        parts = line.split()
        if len(parts) < 8 or not parts[0].isdigit():
            continue
        try:
            rows.append([parts[2], parts[1], parts[0], float(parts[-2])])
        except (ValueError, IndexError):
            continue
    rows.sort(key=lambda r: r[3])
    write_tsv(out, ["peptide", "allele", "position", "rank_percent"], rows)
    log(f"mhc2: {len(rows)} rows -> {out}  (class II ranks are less reliable "
        f"than class I; raw output kept at {raw})")
    write_meta(out, {"predictor": "netMHCIIpan", "alleles": alleles,
                     "length": args.length})


# -------------------------------------------------------------------------- iedb


def cmd_iedb(args: argparse.Namespace) -> None:
    import requests

    base = "https://query-api.iedb.org/epitope_search"
    params = {
        "order": "structure_iri",
        "limit": str(args.limit),
        "select": ("structure_descriptions,source_organism_names,"
                   "parent_source_antigen_names,structure_type,"
                   "e_related_object_source_molecule_iri"),
    }
    if args.organism.isdigit():
        params["parent_source_antigen_source_org_ids"] = f"cs.{{{args.organism}}}"
    else:
        params["source_organism_names"] = f"cs.{{{args.organism}}}"

    log(f"iedb: querying for {args.organism}")
    try:
        resp = requests.get(base, params=params, timeout=90)
        resp.raise_for_status()
        data = resp.json()
    except Exception as exc:  # noqa: BLE001
        die(f"IEDB query failed: {exc}\n"
            "  The IEDB API surface changes; check https://query-api.iedb.org/docs\n"
            "  You can also export a search from https://www.iedb.org/ manually.")

    rows = []
    for rec in data:
        desc = rec.get("structure_descriptions")
        if isinstance(desc, list):
            desc = desc[0] if desc else None
        if not desc or not re.fullmatch(r"[A-Z]+", str(desc)):
            continue
        orgs = rec.get("source_organism_names") or []
        ags = rec.get("parent_source_antigen_names") or []
        rows.append([
            desc, len(desc), rec.get("structure_type"),
            "; ".join(orgs) if isinstance(orgs, list) else orgs,
            "; ".join(ags) if isinstance(ags, list) else ags,
        ])

    out = Path(args.out)
    write_tsv(out, ["peptide", "length", "structure_type",
                    "source_organism", "source_antigen"], rows)
    log(f"iedb: {len(rows)} validated linear epitopes -> {out}  "
        f"(these outrank any prediction)")
    write_meta(out, {"source": "IEDB query API", "query": args.organism,
                     "retrieved": len(rows)})


# ------------------------------------------------------------------ conservation


def cmd_conservation(args: argparse.Namespace) -> None:
    _, rows = read_tsv(Path(args.candidates))
    pep_col = peptide_column(rows)
    variants = read_fasta(Path(args.variants))
    total = len(variants)
    log(f"conservation: {len(rows)} candidates vs {total} variant sequences")

    concat = [seq for _, seq in variants]
    out_rows = []
    for row in rows:
        pep = row[pep_col]
        hits = sum(1 for seq in concat if pep in seq)
        out_rows.append([
            pep, row.get("source"), row.get("position"),
            hits, total, round(hits / total, 4) if total else None,
        ])
    out_rows.sort(key=lambda r: (-(r[5] or 0), r[0]))

    out = Path(args.out)
    write_tsv(out, ["peptide", "source", "position",
                    "n_intact", "n_sequences", "conserved_fraction"], out_rows)
    fully = sum(1 for r in out_rows if r[5] == 1.0)
    log(f"conservation: {fully}/{len(out_rows)} intact in all sampled sequences "
        f"-> {out}  (check sampling depth before trusting a high fraction)")
    write_meta(out, {"variant_file": args.variants, "n_sequences": total})


# --------------------------------------------------------------------- selfcheck


def best_identity(pep: str, proteome: list[str], min_identity: float) -> tuple[float, str | None]:
    """Exact-substring first, then ungapped sliding-window identity."""
    k = len(pep)
    for seq in proteome:
        if pep in seq:
            return 1.0, seq[:0] or pep
    best, best_match = 0.0, None
    threshold = int(min_identity * k)
    for seq in proteome:
        if len(seq) < k:
            continue
        for i in range(len(seq) - k + 1):
            window = seq[i : i + k]
            matches = 0
            for a, b in zip(pep, window):
                if a == b:
                    matches += 1
            if matches > best:
                best, best_match = matches, window
                if matches == k:
                    return 1.0, window
        if best < threshold and best / k < 0.4:
            continue
    return best / k, best_match


def cmd_selfcheck(args: argparse.Namespace) -> None:
    _, rows = read_tsv(Path(args.candidates))
    pep_col = peptide_column(rows)
    human = [seq for _, seq in read_fasta(Path(args.human_proteome))]
    log(f"selfcheck: {len(rows)} candidates vs {len(human)} human proteins "
        f"(this is a safety filter, not an optimisation)")

    # Index exact k-mers for speed; fall back to windowed identity for near-matches.
    out_rows = []
    for n, row in enumerate(rows, 1):
        pep = row[pep_col]
        identity, match = best_identity(pep, human, args.threshold)
        flag = "EXACT_SELF_MATCH" if identity == 1.0 else (
            "NEAR_SELF_MATCH" if identity >= args.threshold else "ok"
        )
        out_rows.append(list(row.values()) + [round(identity, 3), match, flag])
        if n % 50 == 0:
            log(f"  {n}/{len(rows)}")

    header = list(rows[0].keys()) + ["human_max_identity", "human_best_match", "self_flag"]
    out = Path(args.out)
    write_tsv(out, header, out_rows)
    flagged = sum(1 for r in out_rows if r[-1] != "ok")
    log(f"selfcheck: {flagged} candidates flagged -> {out}")
    if flagged:
        log("  Review every flagged candidate with someone qualified to assess "
            "autoimmunity risk. Do not silently rank them.")
    write_meta(out, {"human_proteome": args.human_proteome,
                     "threshold": args.threshold, "flagged": flagged})


# ---------------------------------------------------------------------- coverage


def cmd_coverage(args: argparse.Namespace) -> None:
    _, rows = read_tsv(Path(args.candidates))
    pep_col = peptide_column(rows)
    allele_col = next((c for c in ("best_allele", "allele") if c in rows[0]), None)
    if not allele_col:
        die("input has no allele column; run mhc1 or mhc2 first")

    pairs = defaultdict(set)
    for row in rows:
        if row.get("self_flag", "ok") != "ok":
            continue
        pairs[row[pep_col]].add(row[allele_col])

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    tool_input = out.with_suffix(".population_input.txt")
    with tool_input.open("w") as fh:
        for pep, alleles in sorted(pairs.items()):
            fh.write(f"{pep}\t{','.join(sorted(alleles))}\n")

    write_tsv(out, ["peptide", "alleles", "n_alleles"],
              [[p, ",".join(sorted(a)), len(a)] for p, a in sorted(pairs.items())])

    log(f"coverage: {len(pairs)} self-clean peptide/allele sets -> {out}")
    log(f"  Submit {tool_input} to the IEDB population coverage tool:")
    log("    https://nextgen-tools.iedb.org/pipeline?tool=population_coverage")
    log(f"  Request population '{args.population}' AND per-region breakdowns — a "
        "high world figure can hide a low one in a specific population.")
    write_meta(out, {"population_requested": args.population,
                     "n_peptides": len(pairs),
                     "tool": "IEDB population coverage (manual submission)"})


# ------------------------------------------------------------------------- main


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("prepare"); s.set_defaults(func=cmd_prepare)
    s.add_argument("--fasta", required=True); s.add_argument("--out", required=True)

    s = sub.add_parser("mhc1"); s.set_defaults(func=cmd_mhc1)
    s.add_argument("--fasta", required=True); s.add_argument("--alleles", required=True)
    s.add_argument("--lengths", default="8,9,10,11")
    s.add_argument("--top", type=int, default=0)
    s.add_argument("--out", required=True)

    s = sub.add_parser("mhc2"); s.set_defaults(func=cmd_mhc2)
    s.add_argument("--fasta", required=True); s.add_argument("--alleles", required=True)
    s.add_argument("--length", type=int, default=15); s.add_argument("--out", required=True)

    s = sub.add_parser("iedb"); s.set_defaults(func=cmd_iedb)
    s.add_argument("--organism", required=True)
    s.add_argument("--limit", type=int, default=5000); s.add_argument("--out", required=True)

    s = sub.add_parser("conservation"); s.set_defaults(func=cmd_conservation)
    s.add_argument("--candidates", required=True); s.add_argument("--variants", required=True)
    s.add_argument("--out", required=True)

    s = sub.add_parser("selfcheck"); s.set_defaults(func=cmd_selfcheck)
    s.add_argument("--candidates", required=True)
    s.add_argument("--human-proteome", required=True, dest="human_proteome")
    s.add_argument("--threshold", type=float, default=0.8)
    s.add_argument("--out", required=True)

    s = sub.add_parser("coverage"); s.set_defaults(func=cmd_coverage)
    s.add_argument("--candidates", required=True)
    s.add_argument("--population", default="World"); s.add_argument("--out", required=True)

    args = p.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
