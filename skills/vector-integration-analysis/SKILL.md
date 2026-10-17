---
name: vector-integration-analysis
description: Characterise integrated vector sequence in stable cell lines from WGS or targeted insertion-site sequencing — how many integration sites exist, which construct each junction belongs to when a line carries more than one vector, copy number or VCN per construct, how many copies sit at each site, and which reported junctions are artefacts rather than rearrangements. Covers plasmid transfection, lentiviral and retroviral transduction, PiggyBac and Sleeping Beauty transposition, and HDR knock-in, routing each to the arithmetic that is valid for it. Use this whenever a user mentions integration or insertion sites, junction or breakpoint tables, soft-clipped reads, stable or clonal line characterisation, transgene or vector copy number, VCN, proviral integrity, concatemers or tandem arrays, backbone carryover, or shares an insertion-site report from a sequencing provider — and also when they ask you to interpret, sanity-check or second-guess such a report, even if they never use the word "integration". Especially important for lines built by sequential delivery of constructs that share backbone sequence, where reads cannot be assigned to a construct without the analysis in this skill.
license: Apache-2.0
---

# Vector integration analysis for stable cell lines

A framework for establishing where a construct went in the genome, how much of
it is there, how it is arranged, and how much of a provider's report to
believe.

**The question this skill exists to answer:** given a stable line and a set of
junction calls, how many integration sites are there really, how many copies
sit at each, and which construct does each one belong to?

## When this applies, and when it does not

**Use it for** clonal or oligoclonal stable lines characterised by WGS or
targeted insertion-site sequencing; lines carrying two or more constructs that
share backbone sequence; copy-number and VCN estimation from sequencing depth;
interpretation or QC of a provider's insertion-site report.

**Do not use it for** polyclonal clonal-tracking studies — LAM-PCR or LM-PCR
with UMIs, gene-therapy insertional-mutagenesis surveillance. Those need clonal
abundance, Gini coefficients and Common Integration Site statistics across many
samples: a genuinely different framework that this skill does not reimplement.
Read `references/scope.md` before redirecting; the two are easy to confuse and
the wrong one produces confident nonsense.

**This skill never aligns reads.** It consumes a provider's junction calls and
`samtools depth` output, so the work is table arithmetic and IO — CPU-bound, and
GPU aligners simply do not apply. Even if alignment were in scope, NVIDIA
Parabricks requires ≥16 GB VRAM and ships Linux-only via Docker, which fails both
gates on an 8 GiB Windows card.

## Settle the delivery modality first

**This is the first decision, and the most consequential one in the skill.**
How many copies sit at one site is a property of the delivery mechanism, not of
the data, and it changes the site count several-fold.

| Modality | Flag | Copies per site | Tandem array estimator |
|---|---|---|---|
| Linearized plasmid, naked DNA, minicircle | `nhej` | variable, 1 to tens | **valid** |
| Lentivirus, gammaretrovirus | `integrase` | 1 | **invalid** |
| PiggyBac transposase | `transposon-pb` | 1 | **invalid** |
| Sleeping Beauty transposase | `transposon-sb` | 1 | **invalid** |
| HDR knock-in | `hdr` | 1 per targeted allele, but check for concatemers at the target | check |

Applying the tandem-array estimator to integrase- or transposase-mediated
delivery divides the copy number by a spurious array length and understates the
site count several-fold. Applying single-copy-per-site logic to a transfected
plasmid line overstates it by the same factor.

**Do not take the modality from the paperwork alone.** Stage 2 fingerprints it
from the junction table: integrase and transposase cut the vector at defined
termini, so every integration junction sits at one of two vector coordinates,
while end-joining breaks the molecule wherever it happens to break and vector
breakpoints scatter. When the fingerprint disagrees with the declared modality,
stop and resolve it. `references/modalities.md` covers each mechanism's
junction signature, expected target-site footprint, and characteristic failure.

## The failure modes this skill is built around

Work through these in order.

1. **Shared backbone makes junctions unattributable.** Two constructs from the
   same vector family can be more than half identical. An aligner given only
   one of them files the other's reads under it. Much of the analysis is
   recovering from this.
2. **Cross-mapping manufactures junctions.** Reads crossing from shared into
   construct-unique sequence get soft-clipped and reported as junctions to
   unknown DNA. They are not rearrangements — but they are quantitatively
   useful, so do not simply discard them.
3. **Junction count is not site count, and site count is not copy number.** A
   tandem array of *n* copies gives 2 host flanks and *n*−1 internal junctions.
   These three numbers must be reconciled, not conflated.
4. **Most sites show only one flank.** A one-sided site list is normal: a flank
   landing in repetitive sequence cannot be placed uniquely. An analysis
   assuming full recovery will demand roughly twice as many sites as exist.
5. **Mean coverage over a shared reference is not that construct's copy
   number.** It is the summed depth of every construct carrying that sequence,
   so it is an upper bound. Per-base depth turns it into a measurement.
6. **Provider reports show a curated subset.** Interactive tables that look
   complete routinely display a fraction of the detected junctions, and the
   omitted ones are often the most informative. Work from the full tables.

## Workflow

Input formats, column mapping, clip conventions and naming are covered in
`references/input-formats.md`. Read it if any table fails to parse, or if the
scripts report that they could not interpret the clip side.

### Stage 0 — inventory and QC

Establish what you have before analysing it. Use `references/report-qc.md` as
the checklist. Non-negotiable items:

- **Which genome build?** Often absent from provider reports. Without it no
  locus can be annotated, and builds differ by megabases.
- **Which constructs were in the mapping reference?** If fewer than the line
  carries, every shared-sequence junction is ambiguous until Stage 1.
- **Is junction support post-deduplication?** If not, copy numbers are inflated
  several-fold.
- **Are the delivered tables complete?** Compare the stated detection count
  against the rows in the table.
- **Was per-base depth delivered?** If not, ask. One `samtools depth -a`
  command is the difference between a measured copy number and a bounded guess.
- **Compute genome coverage** — `dedup_pairs × 2 × read_length ÷ genome_size`.
  Every downstream quantity is expressed relative to it.

### Stage 1 — construct identity map

```bash
python scripts/construct_identity.py pAlpha.gb pBeta.gb pGamma.gb \
    --min-block 150 --out identity.json --md identity_map.md
```

Produces, for each construct, a block map labelled by **which constructs carry
that sequence** — not merely whether it is shared — plus the attributable
fraction and the **predicted boundary positions** where cross-mapping artefacts
will appear, naming which construct's reads clip at each one.

`--correspond pAlpha:4809` maps a coordinate onto every other construct, which
is how you recognise that a junction reported against two references is one
event described twice rather than two events.

Run this **before looking at any junction.** Only positions inside a
construct-unique block can be attributed on sequence alone.

Supply a reference for every construct in the line, including ones the provider
did not map against — their sequence is what explains the artefacts their
absence created. Use real plasmid IDs as names; one- and two-letter names are
rejected because they cannot be matched safely against contig names.

### Stage 2 — classify junctions, fingerprint the delivery, build the locus table

```bash
python scripts/classify_junctions.py --identity identity.json \
    --mapping pAlpha --modality integrase \
    --vector-ends 1,8412 --read-length 150 \
    --table parental=junctions_parental.tsv \
    --table treated=junctions_treated.tsv \
    --norm treated=0.87 \
    --out-junctions junctions_pAlpha.tsv --out-loci loci_pAlpha.tsv \
    --md stage2_pAlpha.md
```

Repeat for each construct that was mapped against separately — the comparison
between mappings is what resolves ambiguous junctions.

`--vector-ends` are the defined ends of the delivered molecule in that
construct's numbering: the restriction site for a linearized plasmid, the LTR
termini for a provirus, the ITR termini for a transposon. **Its meaning depends
on the modality.** For end-joining delivery a junction there is ambiguous and
gets flagged; for integrase or transposase delivery it is the canonical
junction, and a host junction *elsewhere* in the vector is what deserves
suspicion.

The script flags several artefact classes automatically. **Review the flags
rather than trusting them** — the boundary and filler-length rules are
reliable; the short-range inverted rule is a heuristic.
`references/artefacts.md` gives the confirming test for each.

Three outputs deserve attention before anything is counted:

- **The delivery fingerprint**, against the declared modality.
- **The target-site footprint** of each paired flank. A duplicated target site
  is the enzyme's signature — 4 bp for PiggyBac and gammaretrovirus, 5 bp for
  HIV-1, 2 bp for Sleeping Beauty — and a library chimera has no reason to
  carry one. Footprints consistent with the modality are positive, independent
  evidence that a pair really is one insertion.
- **Whether the clip side was interpretable at all.** Without it nothing can be
  paired, every site looks single-flanked, and the site count loses its main
  consistency check. The script says so loudly rather than degrading quietly.

Flank pairing is strict: only an `R` + `L` pair within the window is a matched
pair. Two junctions of the same clip side are two left flanks, not a bracketed
site. Getting this wrong inflates the site count.

### Stage 3 — copy number and architecture

Copy `assets/copy_number_config.example.json`, fill it in, then:

```bash
python scripts/copy_number.py config.json --md accounting.md --json accounting.json
```

**Prefer the measured route.** Given per-base depth over each vector reference,
every membership class in the identity map yields one equation — depth over
bases carried by constructs *S* equals the summed depth of those constructs —
and the whole set is solved at once by weighted least squares. This works for
any number of constructs and is a measurement, not an inference.

Without depth the script falls back: the cross-mapping boundary estimator where
exactly two constructs are in play, otherwise an explicit upper bound. **Every
fallback caveat is printed at the top of the report**, because an undecomposed
number that looks like a measured one is how this analysis goes wrong quietly.

Derivations, assumptions and failure conditions: `references/quantification.md`.

### Stage 4 — reconcile, and only then count

The site count is credible only when the numbers agree:

```
copies per haploid ÷ copies per site ≈ number of loci detected
```

Agreement within ~25% means the detected loci account for the DNA present. A
large shortfall means missing sites or an extrachromosomal pool. A large excess
means loci are being double-counted — check flank pairing first.

**Calibrate.** If any locus has independently known structure (both flanks, in
sequence unique to one construct, single copy), run the conversion factor
against it. Agreement validates the conversion factor, the anchor assumption
and the ploidy estimate at once; disagreement by more than ~2× means one of the
three is wrong, and ploidy is the usual culprit.

### Stage 5 — tier the evidence, then report

Never present one undifferentiated site list. Split it:

- **Certain** — vector side in sequence unique to one construct.
- **Inferred** — vector side in shared sequence, assigned by argument. State
  the argument, its margin, and what would overturn it.

For each inferred locus record which supporting strands apply: position at a
defined vector end, a target-site footprint matching the modality, cross-sample
scaling consistent with that construct's cohort, or none. Loci resting on a
single strand are the validation priorities.

`assets/report_template.md` carries the structure. `references/validation.md`
turns each open question into a bench experiment.

## Rules that are easy to get wrong

- **Clip convention.** `R` = reference continues up to that base, then novel
  sequence (the **left** flank of an insert). `L` = novel sequence, then
  reference resumes (the **right** flank). Callers differ, and some emit a
  strand instead — which cannot be converted to a clip side, because it says
  which way the read aligned, not which side of the breakpoint the reference
  continues on. Confirm the convention against a junction of known orientation
  and declare it with `--clip-map`.
- **Target-site footprint sign.** `footprint = host_pos(right flank) −
  host_pos(left flank) − 1`. Negative means a target-site duplication of that
  many bases; positive means bases were deleted at the target; zero is a clean
  join.
- **A junction joining the two ends of a linearized reference is ambiguous.**
  It is produced by a genuine head-to-tail concatemer, by uncut input molecule,
  *and* by any co-resident construct whose equivalent site is intact. Excluded
  from array-length estimates by default; say which you did.
- **Long fillers mean cross-mapping, not biology.** Real end-joining junctions
  carry a few bases of microhomology or filler. Fillers near read length minus
  the aligner anchor are unmapped read remainder.
- **More junctions in a deeper sample is not more rearrangement.** Doubling
  depth pushes sub-threshold events over the cutoff. Compare host–vector
  junction counts, not total junction counts.
- **Read support scales with copy number, not clone abundance.** In a clonal
  line there is one clone; relative-abundance and Gini metrics do not apply.
- **A self-junction is real.** A vector position joined to itself with the same
  clip side is a palindromic head-to-head or tail-to-tail concatemer junction.
  The artefact class is the one with a short but non-zero span.
- **Never state a site count without its uncertainty**, and without saying
  which loci are certain and which inferred.

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

Scripts need Python 3.8+ and no third-party packages. `tests/run_tests.sh`
exercises every path against synthetic fixtures.
