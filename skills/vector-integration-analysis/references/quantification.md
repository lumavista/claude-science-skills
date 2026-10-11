# Quantification — derivations, assumptions and failure conditions

Every quantity `copy_number.py` produces is derived here. Read this before
defending a number to a reviewer, and before changing a default.

Worked examples use round illustrative numbers, not figures from any dataset.

## Contents

1. Genome coverage and units
2. Copy number, route 1: solved from per-base depth
3. Copy number, route 2: the boundary-junction estimator
4. Copy number, route 3: an honest upper bound
5. Conversion factor (reads per junction per copy)
6. Copies per site
7. Site count, reconciliation and flank detection
8. Calibration
9. Cross-sample scaling as a vector-of-origin discriminator
10. Assumptions, and what breaks each one

---

## 1. Genome coverage and units

```
C = dedup_pairs × 2 × read_length ÷ genome_size
```

Use **deduplicated** pairs. PCR duplicates inflate junction support without
adding independent evidence, and the junction caller has normally removed them,
so the denominator must match.

`genome_size` is the haploid genome size of the host — set it for the organism
the line is, not for human by default.

Everything downstream is in **copies per haploid genome equivalent**.
Converting to copies per cell needs the ploidy: a pseudo-triploid line carries
about 3 haploid equivalents, so multiply by 3. Many immortalised lines are
aneuploid and drift on passage, so treat ploidy as a measured quantity, not a
constant, and validate it against karyotype or SNP array if the number matters.

**State which unit you are using.** Conflating copies per haploid genome with
copies per cell is the commonest error in copy-number reporting, and it is a
silent factor of 2 to 4.

---

## 2. Copy number, route 1: solved from per-base depth

**This is the preferred route. Ask for the depth files.**

The identity map labels every base of every vector reference with the set of
constructs carrying it — its *membership class*. Over the reference for
construct *i*, the expected depth at a base whose membership class is *S* is

```
depth(S) = Σ_{j ∈ S} P_j
```

where `P_j` is the true depth of construct *j*. Each membership class in each
reference therefore gives one linear equation, and with *n* constructs you have
far more equations than unknowns. `copy_number.py` solves the stacked system by
weighted least squares, weighting each class by the number of bases it
contributes, and clamps negative solutions to zero.

Copies per haploid genome: divide each `P_j` by `C`.

### Why block edges must be trimmed

Within roughly one read length of a membership-class boundary, a read cannot
lie wholly inside the block, so cross-mapping tapers off and the observed depth
is neither the shared nor the unique value. Averaging across untrimmed blocks
is the reason naive decomposition comes out several-fold high. The script trims
`read_length` bases from each end of every block by default (`--edge-trim`) and
drops classes with too few usable bases (`--min-class-bp`).

This also explains why nominal and effective shared fractions differ: shared
blocks shorter than a read contribute nothing to cross-mapping at all.

### Reading the residuals

The report prints observed against predicted depth per class. Large residuals
mean the block map does not explain the depth. Candidates, in order of
likelihood:

- a vector element also present in the host genome — a human promoter, a poly(A)
  signal, a genomic fragment used as a homology arm — so host reads pile onto
  the vector reference;
- an unlisted construct or helper plasmid in the line;
- non-uniform library coverage, typically GC-driven;
- a partially integrated construct, where the cassette and the backbone
  genuinely carry different copy numbers. This one is a result, not an error.

### Getting the depth

```bash
samtools depth -a -r <vector_contig> SAMPLE.bam > depth_<construct>.txt
```

Tiny output, one command, and it converts the copy number from bounded to
measured. Request it for every vector reference that was mapped against.

---

## 3. Copy number, route 2: the boundary-junction estimator

Used when there is no depth file and **exactly two** constructs are in play.

Reads from construct B that cross out of shared sequence into B-unique sequence
get clipped at a membership boundary of A. Their number is proportional to the
copy number of **B**, not A. So with both mappings available:

```
boundary_signal(mapping A) ∝ copies(B)
boundary_signal(mapping B) ∝ copies(A)

r = boundary_signal(B) / boundary_signal(A) = P_A / P_B
```

Over the construct-A reference, unique blocks receive `P_A` and shared blocks
receive `P_A + P_B`. With `f` the fraction of A shared with B:

```
mean_A = (1 − f)·P_A + f·(P_A + P_B) = P_A · (1 + f/r)

P_A = mean_A / (1 + f/r)
P_B = P_A / r
```

**Report the ratio `r` as the robust quantity and the absolute values as a
lower bound**, because the effective cross-mappable fraction is smaller than
the nominal shared fraction (section 2). Build downstream conclusions on the
ratio wherever possible.

This estimator requires both mappings to exist. If only one was done, ask the
provider to run the other — it is cheap, and it is the difference between a
bounded and an unbounded answer.

---

## 4. Copy number, route 3: an honest upper bound

With three or more constructs and no depth, there is no estimator. Mean
coverage over one reference is the summed depth of every construct carrying
that sequence, so each figure is an **upper bound**, inflated roughly in
proportion to the shared fraction.

`copy_number.py` reports these as upper bounds, prints the caveat at the top of
the report, and labels the route in the abundance table. Everything derived
from them — copies per site, implied sites, flank detection rate — inherits the
inflation.

A report that presents an undecomposed blend as a copy number is the failure
mode this route exists to prevent. If you reach route 3, the correct next action
is to request per-base depth, not to publish the number.

---

## 5. Conversion factor

A read supports a junction only if it spans the breakpoint with enough sequence
on each side to align. With read length `L` and minimum aligner anchor `m`:

```
R = C × (L − 2m) / L          reads per junction per copy per haploid genome
```

Illustrative: `C` = 30×, `L` = 150, `m` = 20 gives `R` = 22 reads per copy.

`m` is a property of the caller and is rarely documented. Treat it as uncertain;
the script sweeps 15–25 bp by default (`--anchors`). The site count is
insensitive to `m`; the flank detection rate is not.

---

## 6. Copies per site

**Mechanism-dependent — see `references/modalities.md` first.**

For integrase- and transposase-mediated delivery, copies per site is **1** and
the estimator below must not be used.

For end-joining delivery, a tandem array of `n` copies at one site has `n − 1`
internal vector–vector junctions and 2 host flanks. Summing read support across
all sites:

```
internal_reads / flank_reads = (n − 1) / 2

n = 1 + 2 × (internal_reads / flank_reads)
```

**Exclude vector-end religation junctions** from `internal_reads` by default.
Against a linear reference they are also produced by uncut input molecule and
by any co-resident construct whose equivalent site is intact — three sources,
indistinguishable by short reads. Including them inflates `n` and deflates the
site count. The script excludes them unless `--include-end-religation` is
passed; report which you used.

Exclude artefact-flagged junctions from both numerator and denominator.

This estimate is an **average over all sites**. It cannot tell one 20-copy
array from twenty single copies. Only long reads, or a per-site assay, can.

---

## 7. Site count, reconciliation and flank detection

```
sites                    = copies_per_haploid ÷ copies_per_site
expected_flank_junctions = 2 × sites
detected_flank_copies    = flank_reads ÷ R
flank_detection_rate     = detected ÷ expected
```

A detection rate well below 100% is normal and needs no exotic explanation:

- roughly half of a mammalian genome is repetitive, so a flank landing in a
  LINE, SINE, satellite or segmental duplication cannot be placed uniquely;
- an array edge may join more vector rather than returning to the chromosome,
  which is a vector–vector junction and invisible as a flank;
- a junction carrying a filler longer than a read can span is never called.

**The reconciliation is the result.** Compare `sites` against the number of
loci actually detected. Agreement within ~25% means the detected loci account
for the DNA present and there is no large undetected reservoir.

- **Shortfall** (more DNA than the loci can hold): sites are missing, or part
  of the vector is extrachromosomal. Distinguish with a plasmid-safe nuclease
  or by passaging without selection — see `references/validation.md`.
- **Excess** (more loci than the DNA supports): loci are being double-counted.
  Check flank pairing first, then whether one site is reported under two nearby
  coordinates in different samples.

If detection had been near 100%, the arithmetic would demand roughly twice as
many sites as were found, and that gap would need explaining. Knowing the
expected detection rate is what stops a normal result being treated as a
discrepancy.

---

## 8. Calibration

Find a locus with independently known structure — ideally both flanks
recovered, in sequence unique to one construct, consistent with a single
insertion on one homologue. In a line of ploidy `p` that is `1/p` copies per
haploid equivalent.

Run its mean flank support through `R` and compare. Agreement within ~30%
validates the conversion factor, the anchor assumption and the ploidy estimate
simultaneously. Disagreement by more than ~2× means one of the three is wrong;
ploidy is the usual culprit, and the calibration can be inverted to estimate it.

This costs nothing and catches a whole class of systematic error. Supply a
calibration locus whenever one exists.

---

## 9. Cross-sample scaling as a vector-of-origin discriminator

When the same line is sequenced in two or more states and genome coverage is
comparable, junction support ratios separate into cohorts. If one construct's
copy number changed and another's did not, their junctions scale with different
ratios.

Calibrate the cohorts on junctions of known origin — those in unique sequence —
then apply to ambiguous ones. Ratios with a 95% interval from the Poisson
approximation:

```
ratio = support₂ × norm / support₁
CI    = ratio × exp(± 1.96 × √(1/support₁ + 1/support₂))
```

`norm` is the ratio of deduplicated read counts between samples, passed as
`--norm SAMPLE=FACTOR`.

**This is suggestive per junction, not decisive.** At support in the tens the
intervals are wide and overlap. Use it as one strand of evidence, recorded per
locus, never alone.

---

## 10. Assumptions, and what breaks each one

| Assumption | Breaks when | Symptom |
|---|---|---|
| Uniform mappability across the vector | vector contains sequence shared with the host genome | depth spikes; large positive residuals; junctions inside the shared element |
| Reads deduplicated before junction calling | provider reports raw support | copy numbers inflated several-fold |
| The identity map covers every construct present | a helper plasmid, third vector or bacterial backbone was carried over | unexplained junctions; depth residuals that no block class explains |
| Copies per site is a mechanism constant | the line received two delivery routes | fingerprint reads MIXED; neither arithmetic fits |
| All vector copies are integrated | episomal or unintegrated vector persists | reconciliation shortfall; copy number falls on passage without selection |
| Ploidy known | aneuploid or drifting line | calibration off by the ploidy error |
| Genome coverage comparable between samples | different library depth | cross-sample ratios meaningless without `norm` |
| Arrays are tandem | inverted or interspersed architecture | the internal:flank ratio still estimates copies per site, but `n` is not a tandem repeat number |
| One flank per site is missing at random | systematic — e.g. the vector end always lands in repeats | detection rate misestimated; site count biased |

When an assumption is untestable from the data, say so in the report rather
than letting it sit silently inside a number.
