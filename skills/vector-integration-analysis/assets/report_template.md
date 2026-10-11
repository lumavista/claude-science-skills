# Vector integration analysis — <line name>

**Project:** <project ID, pipeline and version>
**Samples:** <names and what distinguishes them>
**Host:** <cell line background, organism, genome build, measured ploidy>
**Delivery modality:** <declared, and whether the junction fingerprint agrees>
**Constructs, in order of delivery:**

1. <name, backbone, size, delivery route, defined ends, cargo, selection>
2. <...>
3. <any helper or donor plasmid, even if not expected to persist>

---

## Summary

| | Construct A | Construct B |
|---|---|---|
| Integration sites | | |
| Both flanks recovered | | |
| Copies per haploid genome, <sample1> | | |
| Copies per haploid genome, <sample1> | | |
| Copies per site | | |
| Architecture | | |
| Backbone carryover | | |
| Cassette integrity | | |
| Effect of <the variable distinguishing the samples> | | |

**Headline result.** <two or three sentences>

**How the site count is derived.** <N> junctions, minus <N> cross-mapping
through shared sequence and <N> other flagged artefacts, leaves <N> for this
construct; <N> of those are the paired flanks of one site; giving **<N> distinct
loci**. State any alternative reading and the count it would give.

**How copy number was obtained.** Solved from per-base depth / inferred from
boundary junctions / reported as an upper bound — say which, and if it is not
the first, say what was requested to make it so.

**Main caveat.** <N> loci are proven by sequence; <N> are assigned by argument
(section 5) and should be confirmed before being relied on.

## 1. What was delivered, and how

The modality, the evidence for it from the junction fingerprint, and what
follows arithmetically. If the fingerprint disagreed with the paperwork, this
section is where that is resolved — not a footnote.

## 2. The central problem: which construct did a read come from?

Mechanism, what it did to the delivered data, how it was resolved. Include the
shared percentage, the construct-unique coordinate blocks, and the predicted
artefact boundaries. One worked example of an artefactual junction is worth
more than a paragraph of explanation.

Omit this section only if the line carries a single construct.

## 3. <Construct B> — <result>

## 4. <Construct A> — <N> sites

Locus table: `#`, locus, vector position, target-site footprint, evidence tier
(certain / inferred), support per sample.

### Are these junctions real?

The five arguments from `references/artefacts.md §8`, applied to this dataset
with its numbers. Say explicitly which loci carry a target-site footprint
matching the modality, since that is the strongest single argument where it
applies.

## 5. What is certain, what is inferred

**Certain (<N> loci).** Vector side in sequence unique to one construct.

**Inferred (<N> loci).** State the argument as numbered steps, give its margin,
and say what would overturn it. Then a per-locus table of which evidence
strands apply — position at a defined vector end, footprint matching the
modality, cross-sample scaling — and name the loci resting on a single strand.
Those are the validation priorities.

## 6. Is <N> the whole picture?

The reconciliation: total copies, copies per site, implied sites, loci
detected, and whether they agree. Calibration against a locus of known
structure. Flank detection rate and why it is below 100%. If there is a
shortfall, say what would absorb it and how to test that.

## 7. Integrity of what integrated

Backbone or helper-plasmid carryover. Cassette completeness across the insert.
Whether the copy number is all intended sequence. For a knock-in line, zygosity
and whether the donor integrated by homology or whole.

## 8. Differences between samples

Separate real copy-number or site change from sensitivity change due to
sequencing depth (`artefacts.md §7`).

## 9. What this means for the line

Not just what the data say — what follows for the line's use. Expression
expectations, genetic stability on passage, regulatory or publication exposure,
whether it should be recloned.

## 10. Corrections to the delivered report

Each correction stated as: what the report says, what the data show, why. Keep
it factual and specific; the purpose is a correct record, not a critique of the
provider.

## 11. Recommended next steps

Numbered, in priority order, each tied to the claim it tests and the result
that would overturn the current interpretation.
