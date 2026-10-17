---
name: vector-integration-analysis
description: Characterises integrated vector sequence in stable cell lines from WGS or insertion-site sequencing — site count, per-construct copy number or VCN, copies per site, construct attribution when vectors share backbone, junction artefacts. Use for integration or insertion sites, junction or breakpoint tables, soft-clipped reads, clonal line characterisation, concatemers, or QC of a provider's report.
license: Apache-2.0
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

**Do not take the modality from the paperwork alone.** Stage 2 fingerprints it:
integrase and transposase cut at defined vector termini, so every integration
junction sits at one of two vector coordinates, while end-joining breakpoints
scatter. If fingerprint and declared modality disagree, stop and resolve it.
Per-mechanism signatures: `references/modalities.md`.

## Workflow

Formats, column mapping, clip conventions, naming: `references/input-formats.md`
— read it if a table fails to parse or the clip side is uninterpretable.

### Stage 0 — inventory and QC

Work `references/report-qc.md` as a checklist. Four items gate the rest: the
**genome build** (often unstated, and builds differ by megabases); **which
constructs were in the mapping reference** (fewer than the line carries leaves
every shared-sequence junction ambiguous until Stage 1); **whether junction
support is post-deduplication** (if not, copy numbers inflate several-fold); and
**per-base depth**, which separates a measured copy number from a bounded guess
— ask for `samtools depth -a` if it was not delivered. Then compute
`genome_coverage = dedup_pairs × 2 × read_length ÷ genome_size`; every
downstream quantity is relative to it.

### Stage 1 — construct identity map

```bash
python scripts/construct_identity.py pAlpha.gb pBeta.gb pGamma.gb \
    --min-block 150 --out identity.json --md identity_map.md
```

Yields per construct a block map labelled by **which constructs carry that
sequence** (not merely whether it is shared), the attributable fraction, and the
**predicted boundary positions** where cross-mapping artefacts will appear.
`--correspond pAlpha:4809` maps a coordinate onto every other construct — how
you recognise one event described twice against two references.

Run this **before looking at any junction**, since only construct-unique
positions can be attributed on sequence alone, and supply a reference for every
construct in the line including ones the provider did not map against — their
sequence explains the artefacts their absence created.

### Stage 2 — classify junctions, fingerprint delivery, build the locus table

```bash
python scripts/classify_junctions.py --identity identity.json \
    --mapping pAlpha --modality integrase --vector-ends 1,8412 \
    --read-length 150 --norm treated=0.87 \
    --table parental=junctions_parental.tsv \
    --table treated=junctions_treated.tsv \
    --out-junctions junctions_pAlpha.tsv --out-loci loci_pAlpha.tsv \
    --md stage2_pAlpha.md
```

Repeat per separately mapped construct — comparing mappings is what resolves
ambiguous junctions.

`--vector-ends` are the defined ends of the delivered molecule in that
construct's numbering (restriction site, LTR or ITR termini), and **their
meaning depends on the modality**: under end-joining a junction there is
ambiguous and gets flagged, while under integrase or transposase it is canonical
and a host junction *elsewhere* in the vector deserves the suspicion.

**Review the artefact flags rather than trusting them** — the boundary and
filler-length rules are reliable, the short-range inverted rule is a heuristic;
`references/artefacts.md` gives the confirming test for each.

Check before counting: the fingerprint against the declared modality; the
**target-site footprint** of each paired flank, independent evidence that a pair
is one insertion since a library chimera has no reason to carry one; and
**whether the clip side was interpretable at all**, without which every site
looks single-flanked and the count loses its consistency check.

Flank pairing is strict — only an `R` + `L` pair within the window is a match.
Two junctions of the same clip side are two left flanks, not a bracketed site,
and conflating them inflates the site count.

### Stage 3 — copy number and architecture

Copy `assets/copy_number_config.example.json`, fill it in, then:

```bash
python scripts/copy_number.py config.json --md accounting.md --json accounting.json
```

Mean coverage over a shared reference is an upper bound, not a copy number — it
sums the depth of every construct carrying that sequence.

**Prefer the measured route.** With per-base depth, each membership class in the
identity map gives one equation (depth over bases carried by constructs *S* =
summed depth of those constructs), solved by weighted least squares for any
number of constructs. Without depth the script falls back to the cross-mapping
boundary estimator (only with exactly two constructs) or an explicit upper
bound; **every fallback caveat is printed at the top of the report**, because an
undecomposed number that looks measured is how this analysis goes wrong quietly.
Derivations: `references/quantification.md`.

### Stage 4 — reconcile, then count

```
copies per haploid ÷ copies per site ≈ number of loci detected
```

Agreement within ~25% means the detected loci account for the DNA present; a
large shortfall means missing sites or an extrachromosomal pool, a large excess
means double-counting — check flank pairing first. **Calibrate** against any
locus of independently known structure (both flanks, construct-unique sequence,
single copy): agreement validates the conversion factor, the anchor assumption
and the ploidy estimate at once, and disagreement beyond ~2× means one of the
three is wrong, usually ploidy.

### Stage 5 — tier the evidence, then report

Never present one undifferentiated site list. Split into **certain** (vector
side in construct-unique sequence) and **inferred** (vector side shared,
assigned by argument — state the argument, its margin, and what would overturn
it). Record which strands support each inferred locus: position at a defined
vector end, a footprint matching the modality, cross-sample scaling consistent
with that construct's cohort, or none. Single-strand loci are the validation
priorities. Structure in `assets/report_template.md`; bench follow-up in
`references/validation.md`.

## Rules and failure modes

- **Shared backbone makes junctions unattributable** — constructs from one
  vector family can be over half identical, and an aligner given only one files
  the other's reads under it. Stage 1 exists to recover from that.
- **Cross-mapping manufactures junctions** — reads crossing shared into
  construct-unique sequence get soft-clipped and reported as junctions to
  unknown DNA. Not rearrangements, but quantitatively useful, so keep them.
- **Junction count ≠ site count ≠ copy number.** An array of *n* copies gives 2
  host flanks and *n*−1 internal junctions; reconcile the three.
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
  loci are certain and which inferred.

## Reference material

- `modalities.md` — delivery route; fingerprint vs paperwork conflict
- `input-formats.md` — parse failures, column or clip-side resolution
- `quantification.md` — deriving or defending any number
- `artefacts.md` — classifying junctions, confirming a flag
- `report-qc.md` — Stage 0; report disagrees with the tables
- `validation.md` — junction PCR, ddPCR, long-read follow-up
- `scope.md` — the data turns out to be polyclonal clonal-tracking

Scripts need Python 3.8+ and no third-party packages. `tests/run_tests.sh`
exercises every path against synthetic fixtures.
