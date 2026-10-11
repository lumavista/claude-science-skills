# Turning open questions into bench experiments

Short-read sequencing reaches a hard limit in this analysis. These are the
experiments that close each remaining gap, in rough order of cost.

Tie each experiment to the specific claim it tests, and state what result would
overturn the current interpretation. **A validation plan that cannot fail is
not a validation plan.**

## 1. Junction PCR + Sanger — confirms a specific site

**Design.** One primer in flanking host DNA, one in the vector, oriented so the
product spans the breakpoint. Sanger-sequence the amplicon.

**Primer placement is the whole trick.** Put the vector primer in
**construct-unique sequence** — take the coordinates straight from
`identity.json` — never in shared backbone. A primer in shared backbone
amplifies from either construct and proves nothing about which one is at that
locus. For a line where the constructs share most of their length, this single
decision determines whether the experiment answers the question.

**Orientation.** Work it out from the clip sides before ordering primers:

- host `R` → host runs up to the breakpoint, then the insert begins
- host `L` → insert, then host resumes
- the vector-side clip tells you which way the vector reads from the junction,
  and therefore which direction to walk toward unique sequence

**Prioritise:** loci assigned by inference with no supporting evidence strand;
then the highest-support certain locus as a positive control; then any locus
with both flanks, where the two junctions can be confirmed independently and
the target-site footprint checked directly in the Sanger trace.

**Use a known-good locus as a method control.** A site with both flanks in
unique sequence should amplify cleanly. If it does not, the assay is at fault,
not the site list.

## 2. ddPCR copy number — tests the whole attribution argument

**Design.** One amplicon inside each construct's **unique** region. Never in
shared backbone: a shared amplicon measures the sum and tells you nothing about
the ratio, which is the quantity the attribution argument rests on.

Reference assay: a stable single-copy autosomal locus. In an aneuploid line,
validate the reference against karyotype or SNP array, or the copy number
inherits the ploidy error.

**What it settles.** The absolute copy number of each construct, the ratio
between them, and therefore whether a low-copy construct could plausibly
account for the ambiguous loci. Usually the single highest-value experiment,
and it is cheap.

**For proviral and transposon lines** this is the VCN assay, and because copies
per site is 1 it is also a direct site-count estimate — an independent check on
the sequencing-derived number rather than a different quantity.

## 3. Backbone carryover — a second ddPCR amplicon, almost free

**Design.** One amplicon **outside** the LTRs, ITRs or homology arms — in the
plasmid backbone, origin, or bacterial selection marker — run alongside the
cassette amplicon on the same sample.

**What it settles.** Whether the copy number you measured is all intended
insert. A backbone:cassette ratio near zero is a clean line. A ratio near one
means whole plasmid integrated, the delivery was not what the paperwork says,
and the regulatory and stability picture changes.

For HDR knock-in lines this is the check most often skipped and most often
informative. Run it before anything else on a knock-in line.

## 4. Episomal versus integrated

**Option A.** Plasmid-safe nuclease or exonuclease treatment of the gDNA prep
before ddPCR. Circular extrachromosomal DNA is degraded; integrated copies are
not. Include an undigested control and a spiked circular plasmid control to
confirm the digestion worked.

**Option B.** Passage 15–20 generations without selection and requantify. A
drop indicates an unstable extrachromosomal pool. Slower, needs no enzymes, and
answers the stability question directly, which is usually what was really being
asked.

Run this whenever the site accounting cannot absorb the measured copy number,
or when concatemer junctions scale differently from integration junctions
between samples. For a freshly transduced lentiviral line, run it before
reporting any VCN at all — unintegrated vector persists for days to weeks.

## 5. Functional readout — the question behind the question

If the coding sequence is interrupted at internal positions in the integrated
material, short reads cannot tell you what fraction of copies carries an intact
expressible cassette. Induction plus a western against an epitope tag, or a
reporter assay, answers it directly.

Where expression is what the line is for, this outranks every sequencing
experiment above. A line with twelve copies of which two are intact behaves
like a two-copy line.

## 6. Clonality

A site list means something different in a clonal line than in a pool. If
clonality was assumed rather than established, confirm it before relying on
copies-per-site arithmetic:

- junction support ratios between sites, which should be consistent with
  integer copy numbers in a true clone;
- single-cell cloning and genotyping a panel of subclones at two or three
  junctions by PCR;
- SNP or STR profiling, which also confirms line identity.

Sites present in a subset of subclones mean the population is mixed, and the
copies-per-site estimate was averaging across genotypes.

## 7. Long-read sequencing — the definitive experiment

Nanopore or PacBio on high-molecular-weight DNA. Reads spanning tens of
kilobases cross whole insertion sites in one piece and resolve, simultaneously:

- both flanks of each site, including those lost to repeats in short reads
- array architecture, orientation and true copy number **per site**, which no
  short-read estimator can give
- integrated versus extrachromosomal
- cassette integrity across the whole insert
- which construct, without any inference — a single read carries unique
  sequence from both the construct and the flanking host genome

Use HMW extraction, not a column prep. Where the vector is present at many
copies per genome it is abundant enough that even modest coverage is
informative.

Commission this when architecture matters, when the short-read accounting fails
to reconcile, or when the line is going into anything that will be reviewed
externally.

## Choosing, in practice

| Open question | Cheapest experiment that settles it |
|---|---|
| Is this specific locus real? | junction PCR + Sanger |
| Which construct is at this locus? | junction PCR with a unique-region primer |
| How many copies of each construct? | ddPCR on unique amplicons |
| Is the copy number all intended insert? | backbone ddPCR amplicon |
| Is some of it episomal? | plasmid-safe nuclease, or passage without selection |
| How many copies at *this* site? | long reads |
| Is the cassette intact and expressed? | western or reporter assay |
| Is the line actually clonal? | subclone and genotype |
