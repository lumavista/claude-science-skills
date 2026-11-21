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
| 3 | `copy_number.py` | per-construct copy number solved from per-base depth, per-unique-block estimates and their spread, copies per site, implied site count, flank detection, calibration |
| 4 | `reconcile_sites.py` | **the called-site table**: cross-mapping collapse, retention, construct attribution by strict priority, bounded site count per construct, sub-threshold tail |
| 5 | — | reporting what Stage 4 produced (`assets/report_template.md`) |

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

# 4. ONE called-site table: collapse across mappings, retain, attribute
python scripts/reconcile_sites.py --identity identity.json \
    --loci pAlpha=loci_pAlpha.tsv --loci pBeta=loci_pBeta.tsv \
    --copy-number accounting.json --ploidy 2 \
    --vector-ends pAlpha=1,8412 --vector-ends pBeta=1,7650 \
    --out-called called_sites.tsv --out-subthreshold subthreshold_sites.tsv \
    --summary site_summary.json --md site_summary.md
```

Stage 4 is the site result — report it rather than re-deriving loci by hand.

## The site result is produced, not argued

`reconcile_sites.py` owns locus collapse, retention and attribution, because
leaving them to the caller has produced wrong answers in both available
directions: shared-backbone loci over-assigned to a low-copy construct by
support ratio, and the two mapping runs never collapsed, so every shared site
was listed twice, once under each construct.

- **Collapse.** A junction reported by two mapping runs is one event. Support
  is combined as the per-sample maximum, never summed — the runs report the
  same reads.
- **Retain.** Support ≥ `--min-reads` (5) in ≥ `--min-samples` (2)
  independently prepared samples. The rest is reported as the single-read tail,
  not quietly dropped and not counted.
- **Attribute, by strict priority.** Unique sequence → `certain`. Else exactly
  one construct's declared vector end → `inferred`. Else the copy-number
  ceiling, `ceil(copies × ploidy ÷ min_copies_per_site)`, as a **hard cap**.
  Else `unresolved`: listed once, candidates named, attributed to nobody.
- **The ceiling beats support.** No support ratio or cross-sample scaling
  enters attribution, so none can override it. A construct measured at N copies
  cannot be given more ceiling-attributed sites than N permits, and with two
  constructs of similar copy number the shared loci come out `unresolved` —
  which is the correct answer, not a gap.
- **Count.** Per construct: a floor of certain + vector-end loci, a bounded
  range up to the ceiling, unresolved stated separately.

Every threshold is a documented parameter; nothing is tuned to a dataset. See
`references/attribution.md`.

## Copy number: measured, not estimated

Given per-base depth over each vector reference, each construct's copy number
is read directly from the depth over its own unique blocks, the bases no other
construct carries. Every membership class also yields one equation: depth over
bases carried by constructs *S* equals the summed depth of those constructs.
That stacked system is solved by weighted least squares as a cross-check. A
negative solution, or a deficit in shared-block depth, is reported as evidence
that the additive model is incomplete (typically backbone loss in concatemers),
never clamped to zero. Supply the genome-median depth as `baseline_depth` so the
baseline is filtered the same way as the vector depth.

```bash
samtools depth -a -r <vector_contig> SAMPLE.bam > depth_<construct>.txt
```

Ask for these files. Without them the skill falls back to the cross-mapping
boundary estimator (two constructs only) or to an explicit upper bound, and
**prints the caveat at the top of the report** — an undecomposed blend that
looks like a measured number is how this analysis goes wrong quietly.

The test suite checks the estimates against known ground truth. Three
constructs sharing 71% of their sequence are recovered to within 0.02 copies
per haploid genome. In a two-construct concatemer case (1.5 and 28 copies,
shared backbone 15% under-represented), the low-copy construct is recovered
at 1.5, not clamped to 0, and the deficit is reported at −15%.

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

144 assertions over synthetic fixtures. Every sequence used in testing is
pseudo-random and generated at run time — **no real vector, plasmid or genome
sequence is stored in this repository.**

The suite asserts on content, not just exit status, because a regression here
is a wrong biological answer rather than a cosmetic bug. It covers: one, two
and three constructs; one and two samples; five clip vocabularies including a
strand column and an absent column; host contigs whose names resemble construct
names; each artefact class, and the self-junction that must *not* be flagged;
integrase delivery, where LTR-terminus junctions must not be discarded as
reference-edge artefacts; a declared modality contradicted by the fingerprint;
the copy-number solve against ground truth; a construct whose unique blocks
disagree, where both must contribute and the spread must be reported; and every
path where a caveat must reach the rendered report rather than only the JSON.

For the reconciliation it asserts the properties that matter rather than the
numbers: no locus appears under two constructs; no construct exceeds its
copy-number ceiling; cross-mapping duplicates collapse to one site; the
single-read tail stays out of the called set; a locus at one construct's
declared vector end goes to that construct; read support inflated a
hundred-fold does not move a locus past a full ceiling; and with two
**equal-copy** constructs the shared loci come out unresolved — proof that the
procedure declines to attribute when the data cannot support it.

## Reference material

| File | Read when |
|---|---|
| `references/modalities.md` | settling the delivery route, or the fingerprint disagrees with the paperwork |
| `references/input-formats.md` | a table fails to parse, columns resolve wrongly, or the clip side is rejected |
| `references/quantification.md` | deriving or defending any number |
| `references/attribution.md` | the site priority order, the copy-number ceiling, why a locus is unresolved |
| `references/artefacts.md` | classifying junctions; confirming or overriding a flag |
| `references/report-qc.md` | Stage 0, and whenever a report disagrees with the tables |
| `references/validation.md` | designing junction PCR, ddPCR or long-read follow-up |
| `references/scope.md` | the data turns out to be polyclonal clonal-tracking |
