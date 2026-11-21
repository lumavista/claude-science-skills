---
name: vector-integration-analysis
description: Characterises integrated vector sequence in stable cell lines from WGS or insertion-site sequencing — site count, per-construct copy number or VCN, copies per site, construct attribution when vectors share backbone, junction artefacts. Use for integration or insertion sites, junction or breakpoint tables, soft-clipped reads, clonal line characterisation, concatemers, or QC of a provider's report.
license: MIT
metadata:
  skill-author: LumaVista.BIO
---

# Vector integration analysis for stable cell lines

Given a stable line and a set of junction calls: how many integration sites are
there really, how many copies sit at each, which construct does each belong to,
and how much of the provider's report to believe.

## Scope

**Use for** clonal or oligoclonal lines characterised by WGS or targeted
insertion-site sequencing, copy-number/VCN estimation from depth, and QC of a
provider's report.

**Not for** polyclonal clonal-tracking (LAM-PCR/LM-PCR with UMIs, gene-therapy
insertional-mutagenesis surveillance), which needs clonal abundance, Gini
coefficients and Common Integration Site statistics — a different framework not
reimplemented here, and applying this one produces confident nonsense. Read
`references/scope.md` before redirecting.

**This skill never aligns reads.** It consumes junction calls and `samtools
depth` output, so the work is table arithmetic and IO — CPU-bound, with no GPU
path (and Parabricks, were alignment in scope, needs ≥16 GB VRAM and ships
Linux-only via Docker).

## Settle the delivery modality first

Copies per site is a property of the delivery mechanism, not of the data, and it
changes the site count several-fold.

| Modality | Flag | Copies per site | Tandem array estimator |
|---|---|---|---|
| Linearized plasmid, naked DNA, minicircle | `nhej` | variable, 1 to tens | **valid** |
| Lentivirus, gammaretrovirus | `integrase` | 1 | **invalid** |
| PiggyBac transposase | `transposon-pb` | 1 | **invalid** |
| Sleeping Beauty transposase | `transposon-sb` | 1 | **invalid** |
| HDR knock-in | `hdr` | 1 per targeted allele, check for concatemers | check |

The array estimator on integrase/transposase delivery understates the site count
several-fold; single-copy-per-site logic on a transfected plasmid line
overstates it by the same factor.

**Do not take the modality from the paperwork alone.** Declare what the
paperwork says; Stage 2 fingerprints the junctions and tests the declaration
against them, because integrase and transposase cut at defined vector termini,
so every integration junction sits at one of two vector coordinates, while
end-joining breakpoints scatter. A disagreement fails the
`modality_fingerprint` check. Per-mechanism signatures:
`references/modalities.md`.

## How a run works

Five scripts run in a fixed order, and the scripts carry the analysis: every
number, flag, tier and site in the report is a script output. A sixth,
`make_figures.py`, draws the report's figures from what the five wrote.

```
construct_identity  ->  classify_junctions (per mapping)  ->  host_reference
                                                                   |
                                                              copy_number
                                                                   |
                                   reconcile_sites  <--------------+
                                   (the authoritative site result)
                                          |
                                     make_figures  ->  report, filled from
                                     (the figures)     assets/report_template.md
```

The report is the last stage's output, not a free composition: it is the
template filled with the script outputs and the figures `make_figures.py`
drew. No figure is drawn by hand, and no section structure is invented.

**One file states the run.** `run.json`, the run manifest, holds the facts you
supply once at Stage 0. Every stage takes `--run run.json`, reads what it needs
from it, writes its outputs beside it under fixed names, and records them in it
for the next stage, so no path and no number is carried between stages by hand.
An explicit flag overrides the manifest for that one invocation; `--help` on
each script lists them.

**The exit status says whether to continue.** Run every stage with `--strict`:

| Exit | Meaning | Do |
|---|---|---|
| 0 | every check passed, or only warned — stderr shows each as `CHECK WARNING [id]: …` | read the warnings, then go to the next stage |
| 3 | a check failed; stderr names it `CHECK FAILED [id]: …` | stop, fix the cause it states, rerun this stage |
| 2 | an input could not be read, or was contradictory | fix the input named in the `error:` line |

A check is a condition under which the stage's output must not be used as it
stands. Each stage's checks are also in its JSON (`checks`: `id`, `status`,
`message`) and in a banner at the top of its report. `--accept ID` acknowledges
one failing check so it stops blocking: use it only after investigating that
check, and quote every accepted check in the report. What each one means is
under "Checks" below.

Script paths here are relative to this skill's folder; `run.json` and the data
paths inside it are yours. `via_common.py` and `via_run.py` beside the scripts
are shared modules, not stages: nothing calls them directly.

Formats, column mapping, clip conventions, naming: `references/input-formats.md`
— read it if a table fails to parse or the clip side is uninterpretable.

### Stage 0 — inventory the deliverable, then write `run.json`

Four items about the deliverable gate the rest. Each either gets an answer or
goes into the report as an open item: the **genome build** (often unstated, and
builds differ by megabases); **which constructs were in the mapping reference**
(fewer than the line carries leaves every shared-sequence junction ambiguous
until Stage 1); **whether junction support is post-deduplication** (if not, copy
numbers inflate several-fold); and **per-base depth**, which separates a
measured copy number from a bounded guess — ask for `samtools depth -a` if it
was not delivered, together with the **genome-median per-base depth** from the
same BAM. The genome median (`baseline_depth`) is the better single-genome
baseline: it is filtered the same way as the vector depth, while the read-count
estimate `dedup_pairs × 2 × read_length ÷ genome_size` is not and typically runs
high.

The genome build comes with the organism. If the paperwork names the host
(species, cell line or assembly), put it in the manifest's `host` block. If it
does not, leave the block out: Stage 2b infers the host from the junction
tables and stops to have you ask the user when it cannot. Do not fill in
a human default.

Then copy `assets/run.example.json` to `run.json` in the directory you want the
outputs in — every path inside it, and every file the stages write, resolves
against the manifest's own directory, so a manifest beside the data can name
the data by bare filename. The `_`-prefixed comment keys are ignored; delete
them or leave them. Replace every value: the host (or delete the block), read length, ploidy, modality, each construct's
`vector_ends`, and per sample the junction tables, `dedup_pairs`,
`baseline_depth` and `depth_files`.

**Done when** `run.json` parses, holds no example value, and each of the four
items has an answer or is listed as open.

The full deliverable checklist, and what to request from a provider:
`references/report-qc.md`.

### Stage 1 — construct identity map

```bash
python scripts/construct_identity.py pAlpha.gb pBeta.gb pGamma.gb --run run.json
```

Writes `identity.json` and `identity_map.md`. Three flags are judgements:
`--min-block`, which is your read length (default 150); `--circular`, for
references numbered as the uncut circular plasmid — where the delivered
molecule was cut is a `vector_ends` matter, not this one; and `--k`, the seed
length, which you lower for diverged constructs and raise to ignore short
incidental matches such as restriction and primer sites.

Yields per construct a block map labelled by **which constructs carry that
sequence** (not merely whether it is shared), the attributable fraction, and the
**predicted boundary positions** where cross-mapping artefacts will appear.
`--correspond pAlpha:4809` maps a coordinate onto every other construct — how
you recognise one event described twice against two references.

Run this **before looking at any junction**, since only construct-unique
positions can be attributed on sequence alone, and supply a reference for every
construct in the line including ones the provider did not map against — their
sequence explains the artefacts their absence created.

Add `--correspond pAlpha:4809` to map a coordinate onto every other construct —
how you recognise one event described twice against two references.

**Done when** the exit status is 0 and stderr lists every construct in the line
with its attributable percentage.

### Stage 2 — classify junctions, fingerprint delivery, build the locus table

```bash
python scripts/classify_junctions.py --run run.json --mapping pAlpha --strict
```

Run it once for **each construct that was mapped against separately** —
comparing mappings is what resolves ambiguous junctions. Writes
`classified_<mapping>.tsv`, `loci_<mapping>.tsv`, `stage2_<mapping>.md` and
`stage2_<mapping>.json`.

These `loci_<mapping>.tsv` tables are **intermediates**: Stage 4 collapses them
into the one called-site table. A site in shared sequence appears in every
mapping's table, so never count them directly.

The manifest's `vector_ends` are the defined ends of the delivered molecule in
that construct's numbering (restriction site, LTR or ITR termini), and **their
meaning depends on the modality**: under end-joining a junction there is
ambiguous and gets flagged, while under integrase or transposase it is canonical
and a host junction *elsewhere* in the vector deserves the suspicion.

The stage stops when the junctions contradict the declared modality, when a
fixed-end modality has no vector ends, when the clip side cannot be interpreted,
or when no paired flank carries the modality's target-site footprint. Each of
those changes the site count, so resolve it here rather than downstream.

**Review the artefact flags rather than trusting them** — the boundary and
filler-length rules are reliable, the short-range inverted rule is a heuristic;
`references/artefacts.md` gives the confirming test for each.

Flank pairing is strict — only an `R` + `L` pair within the window is a match.
Two junctions of the same clip side are two left flanks, not a bracketed site,
and conflating them inflates the site count.

**Done when** every mapping has exited 0 and its `stage2_<mapping>.md` reports
how each column and the clip side were read.

### Stage 2b — settle the host genome and fetch its reference

```bash
python scripts/host_reference.py --run run.json --strict
```

Writes `host_reference/host.fai`, `host_reference/host_genes.tsv.gz` and
`host_reference/host.json`, and records the species, assembly, accession,
`genome_size` and how the host was decided under `stages.host`.

The host is decided in this order. **Declared** means the manifest's `host`
block or `--species`, `--assembly`, `--cell-line` names it. **Inferred** means
the host contigs in the Stage 2 tables, and the largest junction coordinate on
each, are matched against the catalogue under every NCBI naming convention
(RefSeq, GenBank, UCSC, bare). A cell-line or species hint narrows the
candidates. A provider `.fai` in `host.fai` decides between assemblies that
share contig names, because lengths are compared too. A declared assembly is
still checked: every host contig must exist in it and every position must fall
inside its contig.

The catalogue holds the current, most-used assembly of each common host:
GRCh38 (HEK293, HeLa, PER.C6 and other human lines), GRCm39 (NS0, Sp2/0),
CriGri-PICRH-1.0 (CHO lines), mRatBN7.2, Chlorocebus_sabeus_1.1 (Vero, COS)
and ROS_Cfam_1.0 (MDCK). `--list` prints it. For any other NCBI assembly,
declare its accession (`"assembly": "GCF_000001405.25"`). For one NCBI does not
hold, declare local `fai` and `annotation` (GTF) files.

**When `host_species` fails, ask the user.** It fails when several assemblies
fit (`chr1`/`chr12` alone fits human, mouse, rat, dog and monkey) or when none
fits (another organism, or another build). The message names the candidates.
Ask which organism the line is and which assembly the reads were mapped
against, record the answer in `host`, and rerun. Never pick a candidate
yourself, and never `--accept` this check.

The stage downloads the NCBI assembly report and RefSeq GTF into
`host_reference/`, once. `--offline` uses only that cache. The `genome_size` it
records is the haploid primary-assembly length, and Stage 3 uses it whenever
the manifest declares none. A declared `genome_size` more than 10% away from it
fails `host_genome_size`.

**Done when** the exit status is 0 and `stages.host` names a species and an
assembly, together with how they were decided (`declared`, `hint` or
`inferred`). Quote that line in the report's Host field.

### Stage 3 — copy number and architecture

```bash
python scripts/copy_number.py --run run.json --strict
```

Writes `accounting.md` and `accounting.json`. `genome_size` comes from the
manifest, or from Stage 2b when the manifest declares none.

The abundance table names two different numbers and divides copy number by only
the second: **genome cov (reads)**, the read-count estimate, which the junction
conversion factor needs because junction support is itself a read count; and
**CN baseline**, the genome-median per-base depth where you supplied one. They
differ by a few percent because they are filtered differently
(`quantification.md` section 1).

Mean coverage over a shared reference is an upper bound, not a copy number — it
sums the depth of every construct carrying that sequence.

**Prefer the measured route.** With per-base depth, each construct's copy number
is the depth over the bases only it carries (its unique blocks), read directly.
Every membership class also gives one equation (depth over bases carried by
constructs *S* = summed depth of those constructs). That system is solved by
weighted least squares, unconstrained, **as a cross-check only**. A negative
solution, or shared blocks carrying less depth than their members' unique
blocks predict, means the additive model is incomplete. In concatemers that is
usually backbone loss. It is reported as a caveat, never clamped to zero.
Each construct's **individual unique-block estimates and their spread** are
reported beside the mean, and a construct resting on little sequence, or on
blocks that disagree, is a flagged answer rather than a silent one.
Without depth the script falls back to the cross-mapping
boundary estimator (only with exactly two constructs) or an explicit upper
bound; **every fallback caveat is printed at the top of the report**, because an
undecomposed number that looks measured is how this analysis goes wrong quietly.

**Calibrate** when any locus has independently known structure (both flanks,
construct-unique sequence, single copy): put it in the manifest's `calibration`
block. Agreement validates the conversion factor, the anchor assumption and the
ploidy estimate at once, and disagreement beyond ~2× means one of the three is
wrong, usually ploidy.

**Done when** the exit status is 0 and every caveat under "Read this before
using any number below" in `accounting.md` is in your notes for the report.

Derivations: `references/quantification.md`.

### Stage 4 — reconcile the loci into one called-site table

```bash
python scripts/reconcile_sites.py --run run.json --strict
```

Writes `called_sites.tsv`, `subthreshold_sites.tsv`, `site_summary.json` and
`site_summary.md`.

**This script is the site result. Report its output; do not re-derive it.**
Locus collapse, retention and construct attribution are done here, by the
script, not by the caller — deterministically, with every threshold a
parameter. Do not assemble a site list from the per-mapping `loci_<mapping>.tsv`
tables: those are intermediates, and every site in shared sequence appears in
all of them.

It does four things, in this order:

1. **Collapses** across mappings. One junction reported by two mapping runs is
   one event, never two. Support is combined as the per-sample maximum, not the
   sum — the runs report the same reads.
2. **Retains** only loci with support ≥ `--min-reads` (5) in ≥ `--min-samples`
   (2) independently prepared samples. The rest go to a sub-threshold table as
   the single-read/noise tail, reported as a tail.
3. **Attributes** by strict priority: **unique** sequence → `certain`;
   else a **declared vector end** of exactly one construct → `inferred`; else
   the **copy-number ceiling**, `ceil(copies × ploidy ÷ min_copies_per_site)`,
   as a **hard cap**, which places a locus only when the leading candidate's
   copy number leads the next by at least `--cn-margin` (2.0); else
   **unresolved**, listed once with its candidates named and assigned to
   nobody. No support ratio or cross-sample scaling
   enters attribution, so none can override the ceiling: a construct measured
   at N copies cannot be given more ceiling-attributed sites than N permits.
4. **Counts**: per construct a floor of certain + vector-end loci, a bounded
   range up to the ceiling, and the unresolved count stated separately.

Three parameters are decisions, so the script asks rather than choosing:

- **One sample only** fails the `min_samples` check, because nothing can be
  confirmed in a second library. Rerun with `--min-samples 1` and say in the
  report that no locus had that confirmation.
- **A concatemer line** raises `ceiling_divisor`: with the default
  `--min-copies-per-site 1` every ceiling is nominal, since the upper bound
  becomes the copy number itself. Rerun with the value the warning names, from
  Stage 3's measured copies per site, to make the ceiling binding — or keep 1
  deliberately and say so.
- **`--cn-margin`** (2.0) sets how far the leading construct must lead on copy
  number before a shared locus goes to it rather than to `unresolved`, so on a
  line whose constructs differ by less than that it decides the split between
  attributed and unresolved. Raise it to be more conservative; check the
  per-construct spread in `accounting.md` before lowering it
  (`attribution.md` section 5).

**Whether the called sites account for the DNA is the last question, and the
script answers it** — per construct, under "Do the called sites account for the
DNA?" in `site_summary.md`, and as the `site_accounting` check. The sites the
DNA implies is the ceiling, so the comparison is the ceiling against what was
attributed:

```
ceil(copies per haploid × ploidy ÷ copies per site)   vs   sites attributed
```

The script calls a shortfall when at least one whole site's worth of DNA is
unaccounted for and fewer than 75% of the implied sites were attributed, both
measured before the ceiling rounds up.

A shortfall means missing sites — flanks lost to repeats are the usual reason —
or an extrachromosomal pool, or more copies per site than the divisor assumes,
which is the first thing to check on a concatemer line. An excess over the
ceiling is the `ceiling_contradiction` check instead. Stage 3's own
reconciliation line uses one mapping's Stage 2 locus count, before collapse, so
treat it as a first look and this table as the result.

**Done when** the exit status is 0, every locus in `called_sites.tsv` carries a
`confidence`, and every construct's line in the accounting table is in your
notes — including any shortfall, with what you think absorbs it.

Rules, parameters and the argument for each: `references/attribution.md`.

### Stage 5 — report what the reconciliation produced

Report `called_sites.tsv` and the site summary as they stand. The evidence
tiers are already in them: **certain** (vector side in construct-unique
sequence), **inferred** (vector-end or copy-ceiling, with the reason recorded
per locus), **unresolved** (candidates named, attributed to nobody). Do not
re-tier, re-argue or redistribute the unresolved loci — they are the validation
priorities, and one junction PCR settles each. Quote the bounded site range and
the unresolved count together; never a single undifferentiated number.

**Draw the figures with the script, then fill the template.** Never draw,
re-plot or describe a figure by hand: every figure in the report comes from
`make_figures.py`, so it is the same for any model and any dataset.

```bash
cp assets/report_template.md report.md     # fill the text from the outputs
python scripts/make_figures.py --run run.json --fill report.md
```

With `--run`, every input comes from the stages' records: the identity map, the
depth files, `accounting.json`, the called sites, the site summary, and the
host `.fai` from Stage 2b. An explicit flag (`--identity`, `--config`,
`--copy-number`, `--called`, `--site-summary`, `--fai`) overrides its input.

It writes standard-library SVG (no plotting package), each figure driven only by
the standard outputs, for any genome and any number of constructs:

| Figure | Placeholder | From |
|---|---|---|
| per-construct coverage track, one per contig per sample: depth, membership blocks shaded unique vs shared and labelled by carrier, trimmed bases, baseline, block means in copies | `[[FIG:coverage_<construct>]]` (all samples) or `[[FIG:coverage_<sample>_<construct>]]` | depth files (via `--config`, `--depth SAMPLE:CONSTRUCT=PATH` or `--depth-dir`), identity.json, the copy-number JSON's baseline and trim |
| genome overview: every chromosome as a length bar, each called locus a tick coloured by construct and styled by confidence | `[[FIG:genome_overview]]` | the host `.fai` (Stage 2b, or `--fai`), called_sites.tsv |
| copy number: bars per construct per sample, one per copy-number JSON (pass two `--copy-number` runs, e.g. two trim policies, to set them side by side), with per-block spread and the shared-block closure deficit | `[[FIG:copy_number]]` | copy_number.py `--json` |
| architecture: site floor, attributed, ceiling and unresolved per construct, and the copies per site that implies | `[[FIG:architecture]]` | site_summary.json |

`--fill` replaces each placeholder with an image link and caption. When a
figure's input is missing, the script skips that figure, records why in
`figures/figures.json` and puts a visible "not generated" note in its slot:
supply the input and rerun rather than substituting a figure of your own.

Structure in `assets/report_template.md`, which ends with the validation
priorities (ddPCR copy number, junction PCR, long reads); bench follow-up in
`references/validation.md`.

**Done when** every section of the template is filled or deleted, every
`[[FIG:...]]` placeholder is filled by `make_figures.py --fill`, the site
range and the unresolved count appear together, and every caveat and every
accepted check from the stage outputs is in the report.

## Checks

What each check means when it fails, and what resolves it. The message on
stderr carries the specifics; this is the standing meaning.

**Stage 1 — `construct_identity.py`**

| Check | Failing means | Resolve by |
|---|---|---|
| `attributable` | a construct carries under 5% unique sequence, so its junctions cannot be attributed on sequence alone and the copy-number solve is weakly determined | checking that every construct has its own reference and that `--k` suits how diverged they are; accepting it, and saying so, where that is simply what the line is |
| `constructs` | *(warns only)* one construct in the map | adding the others' references, including any the provider did not map against |

**Stage 2 — `classify_junctions.py`**

| Check | Failing means | Resolve by |
|---|---|---|
| `vector_ends` | a fixed-end modality was declared with no vector ends, so every genuine junction at a terminus is flagged `reference-edge` and dropped from the locus table | putting the termini in the manifest's `vector_ends`; the message gives the construct's length from the identity map |
| `modality_fingerprint` | the junctions do not look like the declared modality, and the two kinds of delivery use different copies-per-site arithmetic | settling which is right (`modalities.md`): a mislabelled sample, two delivery routes in one line, or a rearranged provirus |
| `clip_side` | the clip column could not be interpreted, so no flank could be paired and every locus is single-flank | declaring the convention with `--clip-map` (`input-formats.md`), or accepting the check when the caller truly reports no clip side |
| `footprint` | no paired flank carries the modality's target-site footprint | checking whether the modality is wrong, the clip convention inverted, or the pairs not pairs |

**Stage 2b — `host_reference.py`**

| Check | Failing means | Resolve by |
|---|---|---|
| `host_species` | the host could not be decided: several catalogue assemblies fit the host contigs equally well, none fits, or there is no evidence and no declaration | asking the user which organism and assembly, and recording it in the manifest's `host` block. Never accept this check |
| `host_contigs` | the junction tables do not fit the decided assembly: a host contig it lacks, a position past a contig's end, or `.fai` lengths that differ. The reads were mapped to another build. It warns only for decoy, viral or HLA contigs | confirming the build with the provider and declaring it in `host.assembly` |
| `host_genome_size` | the declared `genome_size` is more than 10% from the assembly's primary length, so the read-count baseline is wrong by as much | correcting it, or deleting it so Stage 3 uses the assembly's |
| `host_annotation` | *(warns only)* no gene annotation, so sites carry no gene context | `host.annotation`, or rerunning online |
| `host_catalogue` | *(warns only)* some catalogue assemblies could not be fetched and were left out of the inference | rerunning online, or declaring the host |

**Stage 3 — `copy_number.py`**

| Check | Failing means | Resolve by |
|---|---|---|
| `modality_declared` | the modality is `auto`, so copies per site cannot be chosen and the implied site count is unavailable | declaring it from the Stage 2 fingerprint |
| `copy_number_route` | **no** sample was measured from per-base depth: the figures are upper bounds, and everything derived from them inherits that | requesting `samtools depth -a` per vector reference. Where some samples have depth and others do not it warns instead, naming both: the measured samples stand, and a sample with no depth is still usable for junction support in Stage 4 |
| `sample_names` | a sample in the manifest is absent from the junction table, so its flank and internal read counts are zero while numbers are still printed | using the same sample names in Stage 2 and the manifest |
| `estimate_caveats` | *(warns only)* caveats stand against the estimate — a withheld number, an incomplete additive model, a sample left out | quoting them with the numbers |
| `reconciliation` | *(warns only)* implied sites and detected loci disagree | the Stage 4 comparison, which supersedes this one |

**Stage 4 — `reconcile_sites.py`**

| Check | Failing means | Resolve by |
|---|---|---|
| `min_samples` | `--min-samples` cannot be met, so no locus can be called and the site list would be empty | lowering it deliberately and saying so, or adding the other samples' loci tables |
| `copy_number` | a mapped construct has no copy number, so it has no ceiling and the copy-ceiling tier attributes nothing to it | supplying the Stage 3 JSON, or `--copies CONSTRUCT=N` |
| `ceiling_contradiction` | sequence-level evidence gives a construct more sites than its ceiling admits; the sites stand and the disagreement is real | checking the copy number, the ploidy, and whether those sites carry fewer copies each than assumed |
| `ploidy` | *(warns only)* ploidy was assumed to be 2, and it scales every ceiling | measuring it and passing `--ploidy` |
| `vector_ends` | *(warns only)* a construct has no declared ends, so the vector-end tier could not be tested there | putting them in the manifest |
| `samples_consistent` | *(warns only)* the mapping tables do not carry the same samples, so a locus can look sample-specific | noting which run saw which sample |
| `ceiling_divisor` | *(warns only)* the ceiling is nominal where Stage 3 measured arrays | rerunning with the `--min-copies-per-site` the message names |
| `site_accounting` | *(warns only)* a construct's called sites do not account for the copy number measured for it | the shortfall's three candidates above: missing flanks, an extrachromosomal pool, or more copies per site than `--min-copies-per-site` assumes |
| `sites_called` | *(warns only)* no locus was called | the retention parameters, or the input |

## Rules and failure modes

- **Shared backbone makes junctions unattributable** — constructs from one
  vector family can be over half identical, and an aligner given only one files
  the other's reads under it. Stage 1 exists to recover from that.
- **Cross-mapping manufactures junctions** — reads crossing shared into
  construct-unique sequence get soft-clipped and reported as junctions to
  unknown DNA. Not rearrangements, but quantitatively useful, so keep them.
- **Junction count ≠ site count ≠ copy number.** An array of *n* copies gives 2
  host flanks and *n*−1 internal junctions; reconcile the three.
- **One junction, two mapping runs, one site.** A site in shared sequence is
  reported by every mapping run that saw it. Collapsing those is
  `reconcile_sites.py`'s first step; counting the per-mapping tables instead
  doubles the site count and lists the same locus under two constructs.
- **The copy-number ceiling beats every support argument.** A construct
  measured at N copies per haploid cannot hold more ceiling-attributed sites
  than `ceil(N × ploidy ÷ min_copies_per_site)`. Support ratios and
  cross-sample scaling may break ties among constructs that are all admissible
  under the ceiling; they may never fill a construct past it.
- **Most sites show only one flank** — a flank in repetitive sequence cannot be
  placed uniquely, so assuming full recovery demands twice the real site count.
- **Clip convention.** `R` = reference continues to that base, then novel
  sequence (an insert's **left** flank); `L` = the **right** flank. Callers
  differ and a strand is not a clip side, so confirm against a junction of known
  orientation and declare it with `--clip-map` (`input-formats.md`).
- **Footprint sign.** `footprint = host_pos(right) − host_pos(left) − 1`;
  negative is a target-site duplication, positive a deletion, zero a clean join.
- **A junction joining the two ends of a linearized reference is ambiguous** —
  concatemer, uncut input molecule, or a co-resident construct with that site
  intact. Excluded from array-length estimates by default; say which you did.
- **Long fillers mean cross-mapping, not biology**, and a **self-junction with a
  short non-zero span** is the artefact while a zero-span one is a real
  palindromic concatemer junction (`artefacts.md`).
- **More junctions in a deeper sample is not more rearrangement** — depth pushes
  sub-threshold events over the cutoff, so compare host–vector junction counts.
- **Read support scales with copy number, not clone abundance**; a clonal line
  has one clone, so relative-abundance and Gini metrics do not apply.
- **Never state a site count without its uncertainty**, or without saying which
  loci are certain and which inferred. `reconcile_sites.py` reports the range
  and the unresolved count; pass both through.
- **A stage that failed a check has not produced a usable result.** Exit status
  3 means the output is there but must not be reported as it stands; an
  accepted check is a judgement you have made and owe the reader.

## Reference material

- `modalities.md` — delivery route; fingerprint vs paperwork conflict
- `input-formats.md` — parse failures, column or clip-side resolution
- `quantification.md` — deriving or defending any number
- `attribution.md` — the site priority order, the copy-number ceiling
- `artefacts.md` — classifying junctions, confirming a flag
- `report-qc.md` — Stage 0; report disagrees with the tables
- `validation.md` — junction PCR, ddPCR, long-read follow-up
- `scope.md` — the data turns out to be polyclonal clonal-tracking

`assets/run.example.json` is the run manifest template, with every field
documented in place. Scripts need Python 3.8+ and no third-party packages.
`tests/run_tests.sh` exercises every path against synthetic fixtures, and runs
`tests/test_host.sh` (the host stage, offline) as its last step.
