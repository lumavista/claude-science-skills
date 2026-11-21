#!/usr/bin/env python3
"""Epitope triage helpers for protective vaccine design.

Subcommands:
  prepare       Clean and validate an antigen FASTA.
  mhc1          MHC class I binding prediction (MHCflurry, or NetMHCpan if present).
  mhc2          MHC class II binding prediction (NetMHCIIpan; reports if absent).
  iedb          Retrieve experimentally validated epitopes from IEDB.
  conservation  Count candidates intact across variant sequences (fallback —
                see note below).
  selfcheck     Flag candidates with close identity to the human proteome.
  coverage      HLA population coverage via the IEDB tool.

Every subcommand writes a TSV and prints a one-line summary to stderr.
Tool versions are recorded in a sidecar <out>.meta.json for reproducibility.

Two capabilities live outside this script:

  conservation  The `conservation` subcommand here counts exact peptide
                occurrences in unaligned sequences. It does not weight for
                redundant sampling and cannot see per-position or structural
                context, so a clade-skewed panel inflates its fractions. When
                an alignment is available, prefer the `conservation-mapping`
                skill (scripts/conservation_score.py), which applies Henikoff
                weighting, reports effective sample size, and maps scores onto
                structure. Keep this subcommand for the alignment-free case.

  mosaic        Mosaic and cocktail construct design, and junctional
                neo-epitope checking, are in mosaic_design.py beside this
                file.
"""

from __future__ import annotations

import argparse
import json
import os
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
    # Every subcommand prepends '#' banner lines to its own output, so the
    # header is not necessarily line 0. Dropping comments is what makes the
    # documented mhc1 -> conservation -> selfcheck chain runnable.
    lines = [ln for ln in path.read_text().splitlines() if not ln.startswith("#")]
    if not lines:
        die(f"{path} is empty")
    header = lines[0].split("\t")
    rows = [dict(zip(header, ln.split("\t"))) for ln in lines[1:] if ln.strip()]
    return header, rows


def _optional_float(value: str) -> float | None:
    """'' disables a numeric cutoff; anything else must parse as a float."""
    if value is None or value.strip() == "":
        return None
    return float(value)


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

    # setup_models.sh downloads into <skill>/models/mhcflurry, but its `export`
    # died with that shell. Point MHCflurry at the bundled cache here so the
    # skill really does use it automatically, as the docs claim.
    if "MHCFLURRY_MODELS_DIR" not in os.environ:
        bundled = Path(__file__).resolve().parent.parent / "models" / "mhcflurry"
        if bundled.is_dir() and any(bundled.iterdir()):
            os.environ["MHCFLURRY_MODELS_DIR"] = str(bundled)
            log(f"mhc1: using skill-bundled MHCflurry models at {bundled}")

    try:
        from mhcflurry import Class1PresentationPredictor
    except ImportError:
        die("mhcflurry not installed. pip install 'mhcflurry>=2.3,<2.4' && "
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
        # Column order must match the header written below: percentile first,
        # because that is the ranking metric the skill mandates.
        rows.append([
            pep, len(pep), source, pos, allele,
            round(num(rec, "presentation_percentile"), 3),
            round(score, 4),
            round(num(rec, "affinity"), 1),
        ])
    # Ascending: a lower presentation percentile is a better candidate. NaNs last.
    rows.sort(key=lambda r: (r[5] != r[5], r[5]))
    if args.top:
        rows = rows[: args.top]

    out = Path(args.out)
    write_tsv(out, [
        "peptide", "length", "source", "position", "allele",
        "presentation_percentile", "presentation_score", "affinity_nM",
    ], rows)

    # Add header comments to clarify ranking metric
    lines = out.read_text().splitlines()
    with out.open("w") as fh:
        fh.write("# EPITOPE RANKING: Sort by 'presentation_percentile' (≤2%), NOT by 'affinity_nM'\n")
        fh.write("# presentation_percentile is the ranking metric (lower is better).\n")
        fh.write("# affinity_nM and presentation_score are shown for reference only.\n")
        for line in lines:
            fh.write(line + "\n")

    log(f"mhc1: {len(rows)} rows -> {out}  (rank on percentile, not affinity)")
    import mhcflurry
    write_meta(out, {
        "predictor": "mhcflurry",
        "version": mhcflurry.__version__,
        "alleles": usable,
        "lengths": lengths,
        "ranking_metric": "presentation_percentile (≤2% = strong binder)",
        "warning": "Rank ordering is NOT calibrated. Rank 1 vs rank 20 are not reliably distinguishable. Carry ≥15 candidates into screening.",
        "note": "affinity_nM and presentation_score are reference only; do not rank on them."
    })


# -------------------------------------------------------------------------- mhc2


def _dtu_columns(fields: list[str], tool: str) -> dict[str, int | None]:
    """Column indices for a DTU netMHC*-family header line, by NAME.

    The eluted-ligand percentile is mandatory: with -BA the trailing columns
    are "Score_BA Affinity(nM) %Rank_BA", so a positional parse silently ranks
    on affinity in nM, which this skill forbids. Refuse rather than guess.
    """
    def col(*names):
        for n in names:
            if n in fields:
                return fields.index(n)
        return None

    idx = {
        "pos": col("Pos"),
        "pep": col("Peptide"),
        "mhc": col("MHC", "Allele"),
        "el": col("%Rank_EL", "%Rank"),
        "ba": col("%Rank_BA"),
        "aff": col("Affinity(nM)", "Affinity"),
    }
    if idx["el"] is None:
        die("%s output has no %%Rank_EL/%%Rank column; refusing to rank on "
            "affinity. Columns seen: %s" % (tool, " ".join(fields)))
    if idx["pep"] is None or idx["mhc"] is None or idx["pos"] is None:
        die("%s output is missing Pos/Peptide/MHC columns: %s"
            % (tool, " ".join(fields)))
    return idx


def parse_netmhciipan(stdout: str, raw_path: Path | None = None) -> list[list]:
    """Parse netMHCIIpan tabular output by column NAME, never by offset.

    With -BA the trailing columns are "Score_BA Affinity(nM) %Rank_BA", so a
    positional parse (parts[-2]) picks up affinity in nM and ranks class II on
    exactly the metric this skill forbids. Returns rows of
    [peptide, allele, position, %Rank_EL, %Rank_BA, affinity_nM].
    """
    fields = None
    idx: dict[str, int | None] = {}
    rows: list[list] = []
    for line in stdout.splitlines():
        parts = line.split()
        if not parts:
            continue
        if parts[0] == "Pos" and "Peptide" in parts:
            fields = parts
            idx = _dtu_columns(parts, "netMHCIIpan")
            continue
        if fields is None or not parts[0].isdigit():
            continue

        def val(key):
            i = idx.get(key)
            if i is None or i >= len(parts):
                return None
            try:
                return float(parts[i])
            except ValueError:
                return None

        rank_el = val("el")
        if rank_el is None:
            continue
        rows.append([parts[idx["pep"]], parts[idx["mhc"]], parts[idx["pos"]],
                     rank_el, val("ba"), val("aff")])
    if fields is None:
        die("could not find netMHCIIpan's column header in its output"
            + (f"; raw output kept at {raw_path}" if raw_path else ""))
    return rows


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

    rows = parse_netmhciipan(proc.stdout, raw)
    rows.sort(key=lambda r: r[3])
    write_tsv(out, ["peptide", "allele", "position", "rank_percent",
                    "rank_percent_ba", "affinity_nM"], rows)

    # Add header comments to clarify ranking metric and limitations
    lines = out.read_text().splitlines()
    with out.open("w") as fh:
        fh.write("# CLASS II EPITOPE RANKING: Sort by 'rank_percent' = %Rank_EL "
                 "(eluted-ligand percentile, lower is better)\n")
        fh.write("# rank_percent_ba and affinity_nM are reference columns. NEVER rank on affinity_nM.\n")
        fh.write("# ⚠️  CLASS II PREDICTIONS ARE LESS ACCURATE than class I (open groove, shifting registers)\n")
        fh.write("# Do not present class II rankings with equal confidence to class I predictions.\n")
        for line in lines:
            fh.write(line + "\n")

    log(f"mhc2: {len(rows)} rows -> {out}  (class II ranks are less reliable "
        f"than class I; raw output kept at {raw})")
    write_meta(out, {
        "predictor": "netMHCIIpan",
        "alleles": alleles,
        "length": args.length,
        "ranking_metric": "rank_percent = %Rank_EL eluted-ligand percentile (lower is better)",
        "warning": "Class II prediction is substantially less accurate than class I. Do not present with equal confidence.",
        "note": "Negative predictions are weak; absence from this list is not evidence a peptide is not an epitope."
    })


# -------------------------------------------------------------------------- iedb


def cmd_iedb(args: argparse.Namespace) -> None:
    import requests

    base = "https://query-api.iedb.org/epitope_search"
    params = {
        "order": "structure_iri",
        "limit": str(args.limit),
        # Every name here must be a real column: the API rejects the whole
        # query with HTTP 400 on one bad field, and this subcommand did exactly
        # that for every query while it asked for the nonexistent
        # "e_related_object_source_molecule_iri". Verify any addition against
        # `curl 'https://query-api.iedb.org/epitope_search?limit=1&select=*'`.
        "select": ("structure_descriptions,source_organism_names,"
                   "parent_source_antigen_names,structure_type,"
                   "qualitative_measures"),
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

    # An IEDB record is evidence of an assay, not evidence of a POSITIVE assay.
    # qualitative_measures carries the outcome, and negative records were
    # previously written into a table banner-labelled "these outrank any
    # computational prediction" -- the inverse of what that record shows.
    rows = []
    dropped_negative = 0
    for rec in data:
        desc = rec.get("structure_descriptions")
        if isinstance(desc, list):
            desc = desc[0] if desc else None
        if not desc or not re.fullmatch(r"[A-Z]+", str(desc)):
            continue
        measures = rec.get("qualitative_measures") or []
        if isinstance(measures, str):
            measures = [measures]
        positive = [m for m in measures if str(m).lower().startswith("positive")]
        if not args.include_negative and not positive:
            dropped_negative += 1
            continue
        orgs = rec.get("source_organism_names") or []
        ags = rec.get("parent_source_antigen_names") or []
        rows.append([
            desc, len(desc), rec.get("structure_type"),
            "; ".join(orgs) if isinstance(orgs, list) else orgs,
            "; ".join(ags) if isinstance(ags, list) else ags,
            "; ".join(str(m) for m in measures) or "unreported",
        ])

    if dropped_negative:
        log(f"iedb: dropped {dropped_negative} records with no positive assay "
            f"outcome (use --include-negative to keep them, labelled)")
    if len(data) >= args.limit:
        log(f"iedb: WARNING: the response hit --limit ({args.limit}) and this "
            f"endpoint is not paginated here, so the list is TRUNCATED and not "
            f"a complete picture of the evidence. Raise --limit or narrow the "
            f"query before concluding an epitope is unreported.")

    out = Path(args.out)
    write_tsv(out, ["peptide", "length", "structure_type",
                    "source_organism", "source_antigen",
                    "qualitative_measures"], rows)

    # Add header comment to emphasize importance of IEDB validation
    lines = out.read_text().splitlines()
    with out.open("w") as fh:
        fh.write("# EXPERIMENTALLY TESTED EPITOPES from IEDB (positive assays only\n")
        fh.write("#   unless --include-negative; see the qualitative_measures column)\n")
        fh.write("# ✓ A positive assay outranks any computational prediction\n")
        fh.write("# A validated epitope in a well-studied antigen with no IEDB support is worth skepticism\n")
        for line in lines:
            fh.write(line + "\n")

    log(f"iedb: {len(rows)} linear epitopes with a positive assay -> {out}  "
        f"(a positive assay outranks any prediction; read the "
        f"qualitative_measures column, several records carry both outcomes)")
    write_meta(out, {"source": "IEDB query API", "query": args.organism,
                     "retrieved": len(rows),
                     "records_returned": len(data),
                     "limit": args.limit,
                     "truncated": len(data) >= args.limit,
                     "assay_outcome_filter": ("positive only"
                                              if not args.include_negative
                                              else "all outcomes kept"),
                     "dropped_non_positive": dropped_negative,
                     "note": "A POSITIVE assay outranks predictions. Check this list BEFORE ranking. "
                             "Absence here is not evidence against an epitope, especially if truncated=true."})


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

    # Add header comment about sampling bias
    lines = out.read_text().splitlines()
    with out.open("w") as fh:
        fh.write("# CONSERVATION (alignment-free fallback, no Henikoff weighting)\n")
        fh.write("# ⚠️  This is a raw count without correcting for redundant sampling\n")
        fh.write("# A clade-dominated panel will report clade-specific peptides as highly conserved\n")
        fh.write("# Check WHY a region is conserved: functional constraint or shallow sampling?\n")
        for line in lines:
            fh.write(line + "\n")

    log(f"conservation: {fully}/{len(out_rows)} intact in all sampled sequences "
        f"-> {out}  (check sampling depth before trusting a high fraction)")
    write_meta(out, {
        "variant_file": args.variants,
        "n_sequences": total,
        "method": "alignment-free count (no Henikoff weighting)",
        "warning": "Does not correct for redundant sampling. Clade-skewed panel inflates fractions.",
        "note": "For publication, prefer conservation-mapping skill with alignment + Henikoff weighting + effective sample size"
    })


# --------------------------------------------------------------------- selfcheck


def build_proteome_scan(proteome: list[str]):
    """One flat uint8 array of the whole proteome, built once for all peptides.

    The per-peptide pure-Python sliding window was ~1e8 character comparisons
    per candidate (11.4 M residues x 9 offsets), i.e. 30-60 s each and hours for
    a routine 200-candidate run -- which is how a mandatory safety filter ends
    up skipped. Returns None when numpy is unavailable, and the caller falls
    back to the slow path.

    Sequences are joined with a NUL sentinel so a window cannot straddle two
    proteins and report a chimeric self-match.
    """
    try:
        import numpy as np
    except ImportError:
        return None
    joined = "\x00".join(proteome)
    return np.frombuffer(joined.encode("latin-1", "replace"), dtype=np.uint8)


def best_identity_scan(pep: str, arr) -> tuple[float, str | None]:
    """Ungapped best-identity against the flat proteome array, vectorised."""
    import numpy as np

    k = len(pep)
    n = int(arr.size) - k + 1
    if n <= 0:
        return 0.0, None
    matches = np.zeros(n, dtype=np.int16)
    spans_boundary = np.zeros(n, dtype=bool)
    for j, ch in enumerate(pep):
        col = arr[j : j + n]
        matches += (col == ord(ch))
        spans_boundary |= (col == 0)
    matches[spans_boundary] = -1
    idx = int(matches.argmax())
    best = int(matches[idx])
    if best <= 0:
        return 0.0, None
    window = arr[idx : idx + k].tobytes().decode("latin-1")
    return best / k, window


def best_identity(pep: str, proteome: list[str]) -> tuple[float, str | None]:
    """Exact-substring first, then ungapped sliding-window identity.

    Pure-Python fallback used only when numpy is absent. Ungapped: a self-match
    separated by a single indel is not detected by either path.
    """
    k = len(pep)
    for seq in proteome:
        if pep in seq:
            return 1.0, pep
    best, best_match = 0.0, None
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
    return best / k, best_match


def cmd_selfcheck(args: argparse.Namespace) -> None:
    _, rows = read_tsv(Path(args.candidates))
    pep_col = peptide_column(rows)
    human = [seq for _, seq in read_fasta(Path(args.human_proteome))]
    log(f"selfcheck: {len(rows)} candidates vs {len(human)} human proteins "
        f"(this is a safety filter, not an optimisation)")

    # Build the flat proteome array once, not per candidate.
    scan = build_proteome_scan(human)
    if scan is None:
        log("numpy unavailable: using the pure-Python window (minutes per "
            "candidate). Install numpy to make this screen practical.")
    out_rows = []
    for n, row in enumerate(rows, 1):
        pep = row[pep_col]
        identity, match = (best_identity_scan(pep, scan) if scan is not None
                           else best_identity(pep, human))
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

    # Add header comment about mandatory review
    lines = out.read_text().splitlines()
    with out.open("w") as fh:
        fh.write("# SAFETY FILTER: Self-similarity screening (human proteome cross-reactivity check)\n")
        if flagged:
            fh.write("# ⚠️  FLAGGED CANDIDATES REQUIRE EXPERT REVIEW\n")
            fh.write("# Exact or near-exact self-matches carry theoretical autoimmunity risk\n")
            fh.write("# Do NOT silently rank flagged candidates. Review with someone qualified.\n")
        else:
            fh.write("# ✓ No exact or near-exact matches to human proteome\n")
        fh.write("# Also screen separately against known allergens (AllergenOnline, FAO/WHO criterion)\n")
        for line in lines:
            fh.write(line + "\n")

    log(f"selfcheck: {flagged} candidates flagged -> {out}")
    if flagged:
        log("  Review every flagged candidate with someone qualified to assess "
            "autoimmunity risk. Do not silently rank them.")
    write_meta(out, {
        "human_proteome": args.human_proteome,
        "threshold": args.threshold,
        "flagged": flagged,
        "screening_status": "Safety filter (mandatory for anything entering a human)",
        "warning": "Flag and review any candidate with close identity to human protein. Escalate, do not silently rank.",
        "note": "Also screen separately against allergens (AllergenOnline or FAO/WHO ≥35% identity over 80 aa)"
    })


# ---------------------------------------------------------------------- coverage


def _binder_pairs(rows: list[dict], pep_col: str, allele_col: str,
                  rank_col: str | None, max_rank: float | None):
    """{peptide: {alleles}} over predicted binders, plus the drop counts.

    Coverage must be computed over BINDERS. mhc1 applies no cutoff of its own
    (--top defaults to 0 = keep all), so passing its output through whole put
    every enumerated peptide into the coverage set and returned a trivially
    ~100% figure -- on the step that most determines deployability.
    """
    pairs: dict[str, set] = defaultdict(set)
    kept = dropped_self = dropped_rank = 0
    for row in rows:
        if row.get("self_flag", "ok") != "ok":
            dropped_self += 1
            continue
        if max_rank is not None:
            try:
                if float(row[rank_col]) > max_rank:
                    dropped_rank += 1
                    continue
            except (TypeError, ValueError):
                dropped_rank += 1
                continue
        kept += 1
        pairs[row[pep_col]].add(row[allele_col])
    return pairs, kept, dropped_rank, dropped_self


def cmd_coverage(args: argparse.Namespace) -> None:
    _, rows = read_tsv(Path(args.candidates))
    pep_col = peptide_column(rows)
    allele_col = next((c for c in ("best_allele", "allele") if c in rows[0]), None)
    if not allele_col:
        die("input has no allele column; run mhc1 or mhc2 first")

    # Coverage must be computed over PREDICTED BINDERS. mhc1 applies no cutoff
    # (--top defaults to 0 = keep everything), so feeding its output in whole
    # put every enumerated peptide into the coverage set and produced a
    # trivially ~100% figure -- on the step that most determines deployability.
    rank_col = next((c for c in ("presentation_percentile", "rank_percent")
                     if c in rows[0]), None)
    if rank_col is None and args.max_rank is not None:
        die("input has no percentile-rank column (presentation_percentile or "
            "rank_percent); run mhc1/mhc2 first, or pass --max-rank '' to "
            "accept every row and say so in the write-up.")

    pairs, kept, dropped_rank, dropped_self = _binder_pairs(
        rows, pep_col, allele_col, rank_col, args.max_rank)

    if args.max_rank is not None:
        log(f"coverage: kept {kept} rows at {rank_col} <= {args.max_rank}; "
            f"dropped {dropped_rank} above the cutoff and {dropped_self} "
            f"self-flagged")
    else:
        log(f"coverage: NO binder cutoff applied ({kept} rows). Coverage over "
            f"non-binders is not a deployability figure -- state this.")
    if not pairs:
        die("no peptide/allele pairs survived the cutoff; raise --max-rank or "
            "check that mhc1/mhc2 produced binders")

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    tool_input = out.with_suffix(".population_input.txt")
    with tool_input.open("w") as fh:
        for pep, alleles in sorted(pairs.items()):
            fh.write(f"{pep}\t{','.join(sorted(alleles))}\n")

    write_tsv(out, ["peptide", "alleles", "n_alleles"],
              [[p, ",".join(sorted(a)), len(a)] for p, a in sorted(pairs.items())])

    # Add header comment about per-region reporting
    lines = out.read_text().splitlines()
    with out.open("w") as fh:
        fh.write("# POPULATION COVERAGE INPUT: Ready for IEDB analysis\n")
        fh.write("# ⚠️  Request BOTH global and per-region breakdowns\n")
        fh.write("# A construct at 95% world coverage can sit far lower in a specific region\n")
        fh.write("# Report minimum and 10th percentile, not just the aggregate\n")
        for line in lines:
            fh.write(line + "\n")

    log(f"coverage: {len(pairs)} self-clean peptide/allele sets -> {out}")
    log(f"  Submit {tool_input} to the IEDB population coverage tool:")
    log("    https://nextgen-tools.iedb.org/pipeline?tool=population_coverage")
    log(f"  Request population '{args.population}' AND per-region breakdowns — a "
        "high world figure can hide a low one in a specific population.")
    write_meta(out, {
        "population_requested": args.population,
        "n_peptides": len(pairs),
        "binder_cutoff": (f"{rank_col} <= {args.max_rank}"
                          if args.max_rank is not None else "NONE APPLIED"),
        "rows_kept": kept,
        "rows_dropped_above_cutoff": dropped_rank,
        "rows_dropped_self_flagged": dropped_self,
        "tool": "IEDB population coverage (manual submission)",
        "reporting_requirement": "Report per-region coverage, not just global aggregate",
        "warning": "A high global figure can mask low coverage in a specific region. Report minimum and 10th percentile.",
        "note": "The aggregate is the misleading number. Name uncovered strains and state whether gaps were deliberate."
    })


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
    s.add_argument("--include-negative", dest="include_negative",
                   action="store_true",
                   help="keep records whose assays were negative or unreported. "
                        "They are labelled in qualitative_measures, but a "
                        "negative record is evidence AGAINST an epitope.")

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
    s.add_argument("--max-rank", dest="max_rank", type=_optional_float,
                   default=2.0,
                   help="keep only rows at percentile rank <= this (default "
                        "2.0, the conventional class I binder cutoff). Pass an "
                        "empty string to disable, which makes the coverage "
                        "figure meaningless as deployability evidence.")

    args = p.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
