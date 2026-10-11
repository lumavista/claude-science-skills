# vector-integration-analysis

A Claude skill for characterising integrated vector sequence in stable cell
lines, from WGS or targeted insertion-site sequencing.

Answers: how many integration sites are there, how many copies sit at each,
which construct does each junction belong to, and how much of a provider's
report to believe.

## Scope

Built for the cases that break naive analysis:

- **A line carrying two or more constructs that share vector backbone**, where
  reads cannot be assigned to a construct without an explicit sequence-identity
  analysis.
- **Delivery routes with different arithmetic.** Copies per site is a property
  of the mechanism, not of the data. Plasmid transfection makes tandem arrays;
  lentivirus and transposase insert one copy per site. Applying one modality's
  estimator to the other rescales the site count several-fold, and the skill
  routes each to the arithmetic that is valid for it.
- **Copy number from a shared reference**, which is a blend of every construct
  carrying that sequence rather than the copy number of the one you mapped
  against.

**In scope:** clonal and oligoclonal stable lines; `nhej` (linearized plasmid,
naked DNA), `integrase` (lentivirus, gammaretrovirus), `transposon-pb`,
`transposon-sb`, `hdr` knock-in; multi-construct disambiguation; copy number
and VCN; array architecture; backbone carryover; provider-report QC.

**Out of scope:** polyclonal clonal-tracking studies (LAM-PCR with UMIs, clonal
abundance, CIS and hotspot statistics), off-target editing analysis, vector
design, AAV episomal configuration, expression analysis. See
`references/scope.md`.

## What it does

| Stage | Script | Output |
|---|---|---|
| 0 | — | deliverable QC (`references/report-qc.md`) |
| 1 | `construct_identity.py` | per-base construct membership map, attributable fraction, predicted cross-mapping boundaries naming which construct clips where, coordinate correspondence |
| 2 | `classify_junctions.py` | classified junctions, artefact flags, delivery fingerprint, target-site footprint check, locus table with flanks paired and evidence tiered, cross-sample ratios with CIs |
| 3–4 | `copy_number.py` | per-construct copy number solved from per-base depth, copies per site, implied site count, reconciliation, flank detection, calibration |
| 5 | — | evidence tiering and reporting (`assets/report_template.md`) |

## Quick start

```bash
# 1. which parts of each construct are attributable, and where will
#    cross-mapping artefacts appear?
python scripts/construct_identity.py pAlpha.gb pBeta.gb pGamma.gb \
    --min-block 150 --out identity.json --md identity_map.md

# 2. classify junctions and fingerprint the delivery, per mapped reference
python scripts/classify_junctions.py --identity identity.json \
    --mapping pAlpha --modality integrase --vector-ends 1,8412 \
    --table parental=junctions_parental.tsv \
    --table treated=junctions_treated.tsv --norm treated=0.87 \
    --out-junctions junctions_pAlpha.tsv --out-loci loci_pAlpha.tsv \
    --md stage2_pAlpha.md
#    ... repeat for each construct mapped against separately

# 3. copy number, architecture, site accounting
cp assets/copy_number_config.example.json config.json   # fill in every value
python scripts/copy_number.py config.json --md accounting.md --json accounting.json
```

## Copy number: measured, not estimated

Given per-base depth over each vector reference, every membership class in the
identity map yields one equation — depth over bases carried by constructs *S*
equals the summed depth of those constructs — and the whole system is solved at
once by weighted least squares. This works for any number of constructs.

```bash
samtools depth -a -r <vector_contig> SAMPLE.bam > depth_<construct>.txt
```

Ask for these files. Without them the skill falls back to the cross-mapping
boundary estimator (two constructs only) or to an explicit upper bound, and
**prints the caveat at the top of the report** — an undecomposed blend that
looks like a measured number is how this analysis goes wrong quietly.

The test suite verifies the solve against known ground truth: three constructs
sharing 71% of their sequence, recovered to within 0.02 copies per haploid
genome.

## Input formats

Any delimited junction table — tab, comma, semicolon or pipe, headed or
headerless. Columns resolve by header name against a list of common spellings,
positionally against the canonical order

```
ref1  pos1  clip1  filler  ref2  pos2  clip2  support  flag
```

or explicitly via `--columns 'ref1=chromosome,support=n_reads'`.

Clip conventions are resolved from the values present, declared with
`--clip-map`, or — where the column turns out to hold a strand rather than a
clip side — reported as unusable, with the analysis degrading to single-flank
mode rather than inventing flank pairings.

Which side is the host genome is decided by matching reference names against
the constructs in the identity map, so RefSeq accessions, Ensembl names,
scaffolds, alt contigs and non-human assemblies all work without
configuration.

References are GenBank (preferred — feature annotation makes loci
interpretable) or FASTA, one construct per file.

Full details: `references/input-formats.md`.

## Requirements

Python 3.8+. No third-party packages.

## Tests

```bash
tests/run_tests.sh
```

58 assertions over synthetic fixtures. Every sequence used in testing is
pseudo-random and generated at run time — **no real vector, plasmid or genome
sequence is stored in this repository.**

The suite asserts on content, not just exit status, because a regression here
is a wrong biological answer rather than a cosmetic bug. It covers: one, two
and three constructs; one and two samples; five clip vocabularies including a
strand column and an absent column; host contigs whose names resemble construct
names; each artefact class, and the self-junction that must *not* be flagged;
integrase delivery, where LTR-terminus junctions must not be discarded as
reference-edge artefacts; a declared modality contradicted by the fingerprint;
the copy-number solve against ground truth; and every path where a caveat must
reach the rendered report rather than only the JSON.

## Reference material

| File | Read when |
|---|---|
| `references/modalities.md` | settling the delivery route, or the fingerprint disagrees with the paperwork |
| `references/input-formats.md` | a table fails to parse, columns resolve wrongly, or the clip side is rejected |
| `references/quantification.md` | deriving or defending any number |
| `references/artefacts.md` | classifying junctions; confirming or overriding a flag |
| `references/report-qc.md` | Stage 0, and whenever a report disagrees with the tables |
| `references/validation.md` | designing junction PCR, ddPCR or long-read follow-up |
| `references/scope.md` | the data turns out to be polyclonal clonal-tracking |
