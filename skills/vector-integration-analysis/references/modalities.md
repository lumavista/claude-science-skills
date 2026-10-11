# Delivery modalities

How the vector got into the genome determines what its junctions look like, how
many copies sit at one site, and which arithmetic is valid. Settle this before
counting anything.

All numbers in this file are mechanism constants from the published literature,
not measurements from any particular line.

## Why it has to come first

Two quantities are linked by the mechanism, not by the data:

```
integration sites = total copies per haploid genome ÷ copies per site
```

`copies per site` is 1 for integrase- and transposase-mediated delivery and
variable for end-joining delivery. Guessing it wrong rescales the site count by
the whole array length — typically several-fold, sometimes tenfold.

## Fingerprinting it from the data

`classify_junctions.py` reports this automatically. The logic:

**Integrase and transposase cut the vector at defined termini.** Every
integration junction therefore sits at one of exactly two vector coordinates —
the two ends of the LTR or ITR pair. Across twenty sites you see twenty
different host coordinates and two vector coordinates.

**End-joining breaks the molecule wherever it breaks.** Even when a plasmid is
linearized at a single site, the ends are resected and the array internal
junctions form at varied positions, so vector breakpoints scatter.

A fingerprint saying FIXED with a declared modality of `nhej`, or SCATTERED
with a declared modality of `integrase`, is a contradiction that must be
resolved before any number is reported. The commonest innocent explanations are
a mislabelled sample, a line that received both a transfection and a
transduction, and a provirus that has rearranged.

## Target-site footprints

Each enzyme leaves a signature at the host target site. With both flanks
recovered:

```
footprint = host_pos(right flank) − host_pos(left flank) − 1
```

Negative means the two flanks overlap in reference coordinates: a target-site
duplication of that many bases. Positive means host bases were deleted.

| Mechanism | Target site | Expected footprint |
|---|---|---|
| HIV-1 / lentiviral integrase | 5 bp duplicated | −5 |
| Gammaretroviral (MLV) integrase | 4 bp duplicated | −4 |
| PiggyBac transposase | TTAA duplicated | −4 |
| Sleeping Beauty transposase | TA duplicated | −2 |
| Tol2 transposase | 8 bp duplicated | −8 |
| End-joining (plasmid, naked DNA) | variable deletion | 0 to tens positive; small negative where microhomology was used |
| HDR knock-in | designed | ~0 |

**This is one of the few genuinely independent confirmations available.** A
library chimera — vector accidentally ligated to a genomic fragment during prep
— has no reason to carry a target-site duplication of exactly the enzyme's
footprint. A paired flank whose footprint matches the mechanism is strong
evidence it is one real insertion; a set of pairs where none matches means
either the modality is wrong, the clip convention is inverted, or the pairs are
not pairs.

For PiggyBac and Sleeping Beauty you can go further: the duplicated bases
should literally be TTAA or TA. Pull the host sequence at the junction and
check.

---

## `nhej` — linearized plasmid, naked DNA, minicircle

**What happens.** Free DNA ends are joined to each other and to a chromosomal
break by non-homologous end joining. Multiple molecules commonly concatemerise
before or during integration, so one site carries a tandem array of anything
from one to tens of copies, often with inverted members and internal deletions.

**Junction signature.** Vector breakpoints scattered; a few bases of
microhomology or untemplated filler; host target deletions of 0 to tens of
bases. Head-to-tail, head-to-head and tail-to-tail internal junctions all
occur, and a palindromic head-to-head junction appears as a vector position
joined to itself.

**Arithmetic.** The tandem-array estimator applies:

```
copies per site = 1 + 2 × (internal reads / host-flank reads)
```

**Characteristic failures.**
- Uncut or religated input plasmid persisting in the prep, indistinguishable by
  short reads from a genuine head-to-tail concatemer junction.
- Extrachromosomal plasmid surviving selection, which inflates copy number
  without contributing a site. The reconciliation shortfall is the symptom.
- Bacterial backbone — origin, selection marker — integrating alongside the
  cassette. Check whether backbone-only regions carry the same copy number as
  the cassette; a mismatch means partial integration.

---

## `integrase` — lentivirus and gammaretrovirus

**What happens.** Integrase excises the provirus at the LTR termini and inserts
**exactly one copy** per site, with a fixed target-site duplication. There is
no array. Two proviruses can land close together, but each is one copy at one
site.

**Junction signature.** Host–vector junctions at the two LTR termini and
nowhere else; a 4–5 bp target-site duplication; typically 2 bp lost from each
LTR end relative to the plasmid sequence, because integrase trims the terminal
dinucleotide. Build the reference from the proviral form if you can; a transfer
plasmid reference puts the expected junction a couple of bases inside the LTR
and makes every real junction look slightly off-end.

**Arithmetic.** Copies per site is 1, so copy number per haploid genome **is**
the vector copy number, and:

```
integration sites ≈ VCN per haploid genome
```

The tandem-array estimator is invalid and `copy_number.py` refuses to apply it.

**Characteristic failures.**
- **Transfer-plasmid backbone carryover.** Sequence outside the LTRs should not
  be present at all. Depth over the backbone that tracks depth over the
  cassette means plasmid DNA came along, and the copy number is not all
  proviral. This is the single most useful QC check for a lentiviral line and
  the identity map gives it for free: compare depth inside the LTRs against
  depth over the backbone.
- **Unintegrated vector.** Episomal linear and circular forms persist for
  days to weeks after transduction and are counted by any depth- or
  qPCR-based assay. Quantify a line only after enough passages, or test
  directly.
- **Internal deletion or rearrangement.** A provirus that lost part of its
  cassette still produces canonical LTR junctions, so the site count is right
  and the functional conclusion is wrong. Compare depth across the cassette
  against depth at the LTRs.
- **Vector-vector junctions at all.** With one copy per site there should be
  none. Their presence means a rearranged provirus, two adjacent insertions,
  carried-over producer plasmid, or unintegrated vector — worth explaining
  rather than filtering.
- **Self-inactivating vectors.** A ΔU3 SIN LTR makes the two LTRs
  non-identical in the provirus but identical in the plasmid; which reference
  you used changes which junction is attributable.

---

## `transposon-pb` and `transposon-sb` — PiggyBac, Sleeping Beauty

**What happens.** Transposase excises the cassette at the ITRs and inserts one
copy per site into a short target motif — TTAA for PiggyBac, TA for Sleeping
Beauty. Transposition is usually multi-copy across the genome rather than
multi-copy at one site, so a line commonly carries many single-copy sites.

**Junction signature.** Host–vector junctions at the two ITR termini; target
site duplicated exactly (TTAA or TA); no array.

**Arithmetic.** As for `integrase`: copies per site is 1, so copy number
approximates site count.

**Characteristic failures.**
- **Helper plasmid integration.** The transposase-expressing plasmid is
  delivered as DNA and can integrate by end joining — scattered vector
  breakpoints, no TTAA footprint. A line with persistent transposase expression
  keeps remobilising. Include the helper plasmid in the identity map and look
  for it explicitly; its absence from the mapping reference is a common reason
  for unexplained junctions.
- **Donor backbone integration.** Sequence outside the ITRs should not be
  present. Same check as lentiviral backbone carryover.
- **Remobilisation between samples.** Sites that appear or disappear between
  timepoints are the signature. Cross-sample ratios catch it.
- **Local hopping.** Transposons re-insert near the donor site at elevated
  frequency, so several sites clustered within a few hundred kb is expected
  rather than suspicious.

---

## `hdr` — targeted knock-in

**What happens.** The cassette is inserted at a designed locus by
homology-directed repair. One copy per targeted allele — but the same
transfection also produces random end-joining integrations elsewhere, and
concatemers at the target locus are common.

**Junction signature.** At the target, junctions exactly at the designed
homology-arm boundaries, footprint ~0. Elsewhere, an end-joining signature.

**Arithmetic.** Expect copy number ≈ number of targeted alleles. Anything above
that is either a concatemer at the target or a random integration, and the two
must be separated.

**What to check, in order.**
1. **On-target junctions present on both sides**, at the designed coordinates.
   A single-sided junction can mean one homology arm integrated and the other
   did not.
2. **Backbone at the target locus.** Plasmid backbone sequence should be absent
   from a clean HDR event. Backbone present at the target means the donor
   integrated whole rather than by homology. This is the commonest silent
   failure in knock-in line building and it is usually not reported.
3. **Copy number against allele count.** More copies than targeted alleles
   means a concatemer at the target or extra random integrations. Both need
   locating, because they behave differently on passage and in regulatory
   review.
4. **Random integrations elsewhere.** Run the whole workflow with
   `--modality nhej` on the non-target loci; they followed end-joining rules.
5. **Zygosity.** Both alleles targeted, one allele targeted, or one targeted
   and one indel-disrupted are three different lines. Short-read depth at the
   target plus the junction evidence usually distinguishes them; long reads
   settle it.

---

## When you genuinely do not know

Run Stage 2 with `--modality auto`. The fingerprint is reported, no footprint
expectation is tested, and `copy_number.py` refuses to pick a copies-per-site
value and says so in the report. An honest "undetermined" is worth more than a
site count resting on an assumed mechanism.
