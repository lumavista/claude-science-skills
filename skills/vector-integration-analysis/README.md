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
| 0 | — | deliverable QC (`references/report-qc.md`), and the run manifest (`assets/run.example.json`) |
| 1 | `construct_identity.py` | per-base construct membership map, attributable fraction, predicted cross-mapping boundaries naming which construct clips where, coordinate correspondence |
| 2 | `classify_junctions.py` | classified junctions, artefact flags, delivery fingerprint, target-site footprint check, locus table with flanks paired and evidence tiered, cross-sample ratios with CIs |
| 2b | `host_reference.py` | the host organism and assembly (declared, hinted by cell line, or inferred from the host contigs; it stops to have you ask the user when the evidence does not decide), the NCBI assembly report and RefSeq gene annotation, the host `.fai`, and the haploid `genome_size` |
| 3 | `copy_number.py` | per-construct copy number solved from per-base depth, per-unique-block estimates and their spread, copies per site, implied site count, flank detection, calibration |
| 4 | `reconcile_sites.py` | **the called-site table**: cross-mapping collapse, retention, construct attribution by strict priority, bounded site count per construct, sub-threshold tail |
| 5 | `make_figures.py` | the report's standard SVG figures (coverage tracks, genome overview, copy number, site-count bound) |
| 6 | `build_report.py` | **the report**: `assets/report_template.md` filled with the figures, tables and values from the stage outputs, as one self-contained HTML file (figures inline, prints to PDF); a `.docx` too when an SVG rasteriser is installed |

## Quick start

```bash
# 0. state the run's facts once
cp assets/run.example.json run.json      # fill in every value

# 1. which parts of each construct are attributable, and where will
#    cross-mapping artefacts appear?
python scripts/construct_identity.py pAlpha.gb pBeta.gb pGamma.gb --run run.json --strict

# 2. classify junctions and fingerprint the delivery, per mapped reference
python scripts/classify_junctions.py --run run.json --mapping pAlpha --strict
python scripts/classify_junctions.py --run run.json --mapping pBeta  --strict

# 2b. which host and build; fetch its .fai and gene annotation
python scripts/host_reference.py --run run.json --strict

# 3. copy number, architecture, site accounting
python scripts/copy_number.py --run run.json --strict

# 4. ONE called-site table: collapse across mappings, retain, attribute
python scripts/reconcile_sites.py --run run.json --strict

# 5. the report's figures (SVG)
python scripts/make_figures.py --run run.json

# 6. the report: write the prose into a copy of the template, leaving every
#    [[FIG:]], [[TABLE:]] and [[VAL:]] placeholder in place; the script fills
#    them and writes report.html (and report.docx if a rasteriser is present)
cp assets/report_template.md report.md
python scripts/build_report.py --run run.json --report report.md --strict
```

With a host annotation from Stage 2b, `called_sites.tsv` gains, after its
existing columns, each site's host context: `host_region` (exonic, intronic or
intergenic), `host_gene`, and for an intergenic site the `nearest_gene` with its
distance and side.

Stage 4 is the site result — report it rather than re-deriving loci by hand.
The figures come from `make_figures.py`, never drawn by hand: a coverage track
per construct per sample, a genome-wide integration-site overview, copy number
with the shared-block closure, and the site-count bound per construct. Each is
driven only by the standard outputs, so it works for any genome and any number
of constructs. A figure whose input is missing is skipped with a note, not a
failure; an input that exists but cannot be parsed is exit status 2.

The report comes from `build_report.py`, never assembled by hand. It fills the
template's figure, table and value placeholders from the stage outputs and
embeds `make_figures.py`'s SVGs unchanged in one self-contained HTML file: no
external file, no script, no third-party package, so it builds in a locked-down
environment and prints to PDF from any browser. The `.docx` is an optional
extra built only by rasterising those same SVGs (cairosvg or `rsvg-convert`).
Without a rasteriser it is skipped with the reason, and no substitute chart is
ever drawn. An image in the prose that the skill did not draw is removed and
fails the `foreign_figure` check.

Every flag each stage accepts still works; the manifest only removes the need
to pass them. `--help` on any script lists them. The exceptions are
`copy_number.py`'s deprecated trim flags, which are still accepted but print a
`DEPRECATED:` line when used: `--edge-trim` is an alias of `--trim-boundaries`,
`--trim-mode` is superseded by the trim_boundaries policy (`contig-end` runs as
a legacy policy), and `--min-class-bp` is **ignored**, since no class is
dropped on size any more.

## One file states the run; the exit status says whether to continue

`run.json`, the run manifest, holds the facts a run rests on — genome size,
read length, ploidy, modality, each construct's vector ends, and per sample the
junction tables, read counts and depth files. Each stage reads what it needs
from it, writes its outputs beside it under fixed names, and records them in it
for the next stage. Nothing is retyped between stages: not the identity map,
not the classified-junction tables, not the loci tables, not the copy numbers,
not the detected-locus count.

Each stage also emits **checks** — the conditions under which its output must
not be used as it stands — in three places at once: a banner at the top of its
report, a `checks` list in its JSON (`id`, `status`, `message`), and stderr.
With `--strict` a failing check makes the stage exit 3, so a run stops at the
stage that found the problem rather than carrying it into the site count.
`--accept ID` acknowledges one failing check after investigation; it stays in
the output marked `accepted`.

Checks exist for a declared modality the junctions contradict, a fixed-end
modality with no vector ends, an uninterpretable clip side, a footprint that
matches no pair, an undeclared modality, a copy number that was not measured, a
sample name absent from the junction table, an unmeetable `--min-samples`, a
construct with no copy-number ceiling, a ceiling contradicted by sequence
evidence, and a ceiling left nominal on a concatemer line. `SKILL.md` lists what
each one means and what resolves it.

## The site result is produced, not argued

`reconcile_sites.py` owns locus collapse, retention and attribution, because
leaving them to the caller has produced wrong answers in both available
directions: shared-backbone loci over-assigned to a low-copy construct by
support ratio, and the two mapping runs never collapsed, so every shared site
was listed twice, once under each construct.

- **Collapse.** A junction reported by two mapping runs is one event. Support
  is combined as the per-sample maximum, never summed — the runs report the
  same reads.
  Stage 2's locus table is already a union over samples, with each sample's
  support counted over both flanks, so a site present in only one sample
  reaches this stage.
- **Retain.** Total support ≥ `--min-reads` (5), summed over samples, and
  presence in ≥ `--min-samples` (2) independently prepared samples. The rest is
  reported as the single-read tail, not quietly dropped and not counted.
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

291 assertions over synthetic fixtures, including the host-stage suite
(`tests/test_host.sh`), which `run_tests.sh` runs as its last step and which
can also be run on its own. Every sequence used in testing is
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

For the figures it asserts that every SVG `make_figures.py` writes is
well-formed with an `<svg>` root, for the one-sample three-construct case and
the two-construct, two-trim-policy case. The genome overview must draw every
contig of the fixture `.fai`, with one tick per called locus and the
confidence tiers intact. Each coverage track must shade every identity-map
block and show the baseline and the trimmed bases. A missing input must skip
its figure with a reason rather than fail the run. `--fill` must leave no
figure slot unfilled in the report body.

For the report it asserts that `build_report.py` runs under `python -I -S` (no
site-packages) and imports nothing outside the standard library, and that its
HTML embeds every figure inline exactly as `make_figures.py` wrote it (the
coverage tracks and the genome overview included). The HTML must fetch nothing
external, carry the called-site table value for value, and hold no unfilled
`[[...]]` placeholder. A rebuild must be byte-identical. Without a rasteriser
the docx is skipped with its reason; with one, the docx holds one image per
skill figure and no other. An image pasted into the prose fails `--strict`
and cannot be accepted away. A dropped figure slot appends the figure rather
than losing it, and a missing `figures.json` stops the build.

For the run manifest and the checks it asserts that the manifest route
reproduces the flag-and-config route byte for byte on the loci, called-site and
classified-junction tables and value for value on the copy numbers and the site
summary; that every stage records its outputs while leaving the declared facts
untouched; that an explicit flag still overrides the manifest; that each
blocking condition exits 3 under `--strict`, exits 0 without it, reaches the
report as a banner and the JSON as a status, and unblocks under `--accept`; and
that a sample with depth is not spoiled by a sample without it; that the
shortfall between a construct's measured DNA and its called sites is reported
per construct rather than for the primary mapping alone; and that the ceiling's
divisor is the copies per site the copy-number JSON measured, falls back to 1
when there is no estimate, and yields to an explicit `--min-copies-per-site`.

`tests/test_host.sh` covers the host stage offline, against a synthetic
catalogue cache. It checks that a RefSeq accession is inferred, that
`chr`-style names fitting several assemblies stop the run instead of being
guessed, and that a cell-line hint settles it. It checks that a contig missing
from the host, or a position past a contig's end, fails, and that decoys only
warn. It checks that the assembly's `genome_size` reproduces the declared
one's copy numbers, and that the annotation only appends columns, leaving the
called table otherwise byte-identical. Finally it checks that exonic,
intronic and intergenic sites are labelled correctly, that `make_figures.py
--run` draws the genome overview from the host `.fai`, that the overview
carries gene labels and the species and assembly when the host stage supplied
them and is the plain `.fai` figure when it did not, that crowded labels never
overlap, and that a declared genome size contradicting the assembly fails.

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

`SKILL.md` is the workflow: the stages in order, what each one writes, what its
checks mean, and the domain rules that govern the analysis.
