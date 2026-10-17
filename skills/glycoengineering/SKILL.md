---
name: glycoengineering
description: Analyzes and engineers protein glycosylation by scanning canonical N-glycosylation sequons, describing S/T-rich regions, checking curated glycan evidence, and preparing NetNGlyc, NetOGlyc and GlycoSHIELD workflows. Use for glycoprotein engineering, antibody Fc glycosylation, glycan shielding, and site-specific glycoproteomics interpretation.
license: Unknown
compatibility: Local sequence helpers require Python 3.10+. Public GlyTouCan lookup requires requests and network access. DTU predictors use a browser; optional GlycoSHIELD needs a separate source installation, glycan libraries, and GROMACS for SASA analysis.
metadata:
  version: "1.4"
  skill-author: Kuan-lin Huang
  last-reviewed: "2026-10-01"
---

# Glycoengineering

## When to use

For canonical N-glycosylation candidate analysis, O-GalNAc candidate triage,
antibody glycoform comparison, or structural glycan shielding. **Keep four kinds
of result separate** — sequence motif, prediction, experimentally supported
occupancy, glycan structure/composition — since none establishes the others.

Scope is eukaryotic secretory-pathway N-glycosylation and mucin-type O-GalNAc;
O-GlcNAc, O-mannose, O-fucose and other O-linked modifications need their own
evidence and predictors, and S/T enrichment does not identify the type.

## Workflow

1. Record accession, isoform/version, sequence, expression host, signal
   peptide/transmembrane topology and construct boundaries, plus a mapping from
   submitted coordinates to mature protein, PDB chain/residue IDs and antibody
   numbering.
2. Scan canonical N-X-[S/T], X not Pro, **retaining overlaps**. A missing
   canonical motif does not rule out unusual N-glycosylation; a present motif
   does not establish occupancy. Prioritize secreted/luminal/extracellular
   regions by topology and curated evidence — a sequence-only hit in a
   cytoplasmic region is not evidence of secretory-pathway glycosylation.
3. For O-GalNAc, treat S/T density as descriptive, then inspect a predictor plus
   cell and transferase context. **Do not exclude SP/TP motifs** — adjacent Pro
   can be favorable depending on the GalNAc-transferase
   ([substrate study](https://doi.org/10.1093/glycob/cwu089)).
4. Propose explicitly numbered mutations, list every changed residue, and rescan
   the product for lost or created overlapping motifs. Evaluate structural and
   expression effects independently of glycosylation.
5. Validate occupancy and glycoforms experimentally. For biosimilar comparisons,
   match sample handling, analytical coverage and quantification before
   comparing glycan percentages.

## Local sequence analysis

Use [scripts/glycoengineering_tools.py](scripts/glycoengineering_tools.py)
(standard library only) with this skill's `scripts/` on `PYTHONPATH`. It takes
raw canonical 20-amino-acid sequences and rejects FASTA headers, gaps, unknown
residues and empty input, so parse FASTA first and do not strip unknown residues
— that shifts the 1-based coordinates it reports.

```python
from glycoengineering_tools import (normalize_sequence,
    find_n_glycosylation_sequons, eliminate_glycosite, add_glycosite,
    find_st_rich_sites)

sequence = normalize_sequence("NNST")
assert [s["position"] for s in find_n_glycosylation_sequons(sequence)] == [1, 2]
assert eliminate_glycosite(sequence, 1, "Q") == "QNST"  # sequon at 2 remains
assert add_glycosite("APA", 1, "T", allow_proline_substitution=True) == "NAT"
```

`eliminate_glycosite` needs a complete canonical sequon and a one-residue
replacement other than N. `add_glycosite` may alter up to three residues (here
A1N, P2A, A3T — the P→A must be explicit) and retains an existing S/T at +2.
**Neither predicts a mutant's fold, occupancy, function or tolerability**, and
an N-to-Q substitution can change protein behavior with no glycan effect.

`find_st_rich_sites(seq, window=7, min_st_fraction=0.4)` reports a fraction of
residues, not an O-glycosylation probability. The window must be a positive odd
integer, terminal windows shorten and use their actual length as denominator,
and isolated real O-GalNAc sites can be missed. In batch work report sequons per
100 residues, keeping identifiers and lengths.

## Prediction services

Use the official forms (input limits are on them) — there is no stable
submission API or job-polling endpoint, and saving a CGI URL is not a submitted
job. Preserve output, input and service version together. The documentation and
sample output were reviewed; no new jobs were submitted.

- [NetNGlyc 1.0](https://services.healthtech.dtu.dk/services/NetNGlyc-1.0/),
  human N-glycosylation potential with jury agreement. Its default 0.5 threshold
  is **not a calibrated occupancy probability**, SignalP runs but extracellular
  topology is not checked, and it may score N-P-S/T, which should be excluded
  absent independent evidence.
- [NetOGlyc 4.0](https://services.healthtech.dtu.dk/services/NetOGlyc-4.0/),
  mammalian mucin-type O-GalNAc, GFF2 scores positive at **greater than** 0.5.
  Prefer the full sequence including signal peptide, and give an isolated site
  15 residues of flanking context each side. A positive supports regional
  likelihood, not occupancy or glycan type beyond O-GalNAc.

## Structural shielding with GlycoSHIELD

[GlycoSHIELD](https://gitlab.mpcdf.mpg.de/dioscuri-biophysics/glycoshield-md/)
grafts pre-simulated glycan conformers onto structures and filters steric
clashes; it does not predict which sequons are occupied and runs no fresh MD per
query ([Tsai et al., 2024](https://doi.org/10.1016/j.cell.2024.01.034)). Pinned
source, installation, direct Python API, input mapping and SASA analysis:
[references/glycoshield.md](references/glycoshield.md).

**The reviewed upstream CLI silently ignores several parsed options, including
`--mode`, `--threshold` and `--dryrun`; use the documented direct API
workaround.** Outputs are per-site PDB/XTC ensembles plus shielding in a PDB
B-factor column, which are not crystallographic temperature factors.

## Database evidence and glycan notation

[references/glycan_databases.md](references/glycan_databases.md) has the live
GlyTouCan SPARQL lookup, GlyConnect access limits, curated resources and
notation rules. Query exact accessions and preserve dataset dates, species,
tissue or cell context and supporting publications. **Missing records or service
errors cannot establish that a protein is unglycosylated.**

A composition such as `Hex:5 HexNAc:4 dHex:1` resolves neither linkages nor
branch positions and does not distinguish GlcNAc from GalNAc, so record an
actual sequence (WURCS/GlycoCT/IUPAC, uncertainty retained) plus accession where
available. Core fucose attaches to the innermost GlcNAc of an N-glycan, not to
core mannose.

## Antibody engineering decisions

Number Fc mutations in an explicit scheme (commonly EU), then map them to the
actual construct — EU N297 is not residue 297 of an isolated Fc FASTA. Fc and
any Fab glycans must be distinguished analytically.

| Objective | Strategy | Interpretation |
| --- | --- | --- |
| ADCC / FcγRIIIa engagement | Reduce Fc core fucosylation, preserving the glycan | No universal fold gain — effect depends on antibody, receptor, assay ([evidence](https://pubmed.ncbi.nlm.nih.gov/22023369/)) |
| Remove canonical Fc N-glycan | N297Q/A/D, or disrupt +2 with T299A | Sequon disruption only; verify structure and effector function |
| Alter serum persistence | Characterize glycoform-specific clearance | High-mannose clearance can increase; sialylation is no universal half-life recipe ([PK study](https://pubmed.ncbi.nlm.nih.gov/21421994/)) |
| Anti-inflammatory Fc activity | Compare defined sialylated glycoforms | Linkage, preparation and model matter; no universal IVIG mechanism follows ([study](https://pubmed.ncbi.nlm.nih.gov/18420934/)) |
| Reduce non-human epitopes | Measure α-Gal and Neu5Gc, select production conditions | Sequence editing alone does not control host glycan processing |
| Epitope accessibility | Add or remove a mapped surface sequon | Confirm occupancy, folding, binding, antigenicity; shielding is a geometric hypothesis |

Fc variants such as S298A/E333A/K334A or F243L combinations can alter both
receptor interaction and host-dependent glycan processing — **do not label F243L
alone a deterministic defucosylation switch**
([2026 study](https://pubmed.ncbi.nlm.nih.gov/41873859/)).

## Interpretation and verification boundary

Glycopeptide searches support peptide identity, composition and sometimes site
localization without resolving a full glycan structure, so review localization
fragments, competing assignments, search space and error control. Composition
does not determine topology, and O-glycosite ambiguity can remain even with an
identified glycopeptide
([Byonic](https://support.proteinmetrics.com/hc/en-us/articles/18139992247060-Byonic-O-Linked-Glycopeptide-Analysis)).
**Relative signal among detected glycoforms is not absolute site occupancy** —
that needs the appropriate modified and unmodified denominator plus analytical
response considerations.

Verification boundary: synthetic tests cover local sequence behavior and public
accession lookup was executed with mocked error handling, but conformer
grafting, GROMACS SASA, DTU jobs, commercial MS software and biological
performance were not.
