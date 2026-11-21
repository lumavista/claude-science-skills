<!--
HOW THIS REPORT IS BUILT. Fill this template; do not invent a structure.

* Every number comes from a script output: identity.json / identity_map.md
  (Stage 1), the classify_junctions.py tables (Stage 2), the copy_number.py
  JSON and Markdown (Stage 3), and called_sites.tsv + site_summary.json from
  reconcile_sites.py (Stage 4). Quote them; do not re-derive them.
* Every figure comes from scripts/make_figures.py. Do not draw, describe-as-
  drawn or substitute any figure by hand. Each [[FIG:<id>]] placeholder below
  is replaced by the script:

      python scripts/make_figures.py --run run.json --fill report.md

  Placeholder ids: [[FIG:coverage_<construct>]] (every sample's track for that
  construct; [[FIG:coverage_<sample>_<construct>]] for one), [[FIG:genome_overview]],
  [[FIG:copy_number]], [[FIG:architecture]]. Replace <construct> with the
  construct's name exactly as in identity.json. A figure that could not be drawn
  is replaced by a visible "not generated" note naming the missing input: leave
  that note in, or supply the input and rerun.
* Keep every section. One per construct for the per-construct results
  (sections 3-4: repeat for a third construct, delete section 3 for one). A
  section that does not apply says so in one line rather than disappearing.
-->

# Vector integration analysis — <line name>

**Project:** <project ID, pipeline and version>
**Samples:** <names and what distinguishes them>
**Host:** <cell line background; organism and assembly (accession) from stages.host, and how they were decided: declared / hint / inferred; measured ploidy>
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
| Copies per haploid genome, <sample2> | | |
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

**How copy number was obtained.** Measured from unique-block depth / inferred
from boundary junctions / reported as an upper bound. Say which, and if it is
not the first, say what was requested to make it so. Name the baseline
(genome-median depth or read count) and the shared-block closure. If the stacked
solve went negative, say so, because it means the additive model was incomplete.

**Main caveat.** <N> loci are proven by sequence; <N> are assigned by argument
(section 5) and should be confirmed before being relied on.

**Checks and caveats carried forward.** Every caveat from the stage reports,
and every check accepted with `--accept`, with what was investigated and why it
does not invalidate the result.

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

Same structure as section 4, scaled to what this construct carries: its sites
(locus table from called_sites.tsv), its copy number, and what its coverage
track shows.

[[FIG:coverage_<construct B>]]

Read the track before quoting the copy number: the unique blocks should sit
flat at one level, the shared blocks at the sum of their carriers, and the
trimmed edges should be the only places the depth tapers. A step inside a
block, or a unique block at a different level from its sibling, is a deletion,
a rearrangement or a mis-specified reference, and belongs in section 7.

## 4. <Construct A> — <N> sites

[[FIG:coverage_<construct A>]]

State the unique-block copy number shown on the track, the number of unique
bases it rests on, and whether the per-block estimates agree.

Locus table: `#`, locus, vector position, target-site footprint, evidence tier
(certain / inferred), support per sample, and the host context from
called_sites.tsv: `host_region` (exonic / intronic / intergenic), `host_gene`, or
`nearest_gene` with its distance and side for an intergenic site. If there was
no host annotation, leave the host-context columns out and say why.

### Are these junctions real?

The five arguments from `references/artefacts.md §8`, applied to this dataset
with its numbers. Say explicitly which loci carry a target-site footprint
matching the modality, since that is the strongest single argument where it
applies.

## 5. What is certain, what is inferred

<!-- Take this straight from reconcile_sites.py: called_sites.tsv and the site
     summary already carry the tier, the evidence and the reason per locus.
     Quote the bounded range and the unresolved count together. Do not re-tier
     or redistribute the unresolved loci. -->

[[FIG:genome_overview]]

Where the loci sit, by attributed construct and confidence. Name any contig
carrying more than one locus, and any locus whose position (centromeric,
subtelomeric, near a gene of known consequence) matters for section 9.

**Certain (<N> loci).** Vector side in sequence unique to one construct.

**Inferred (<N> loci).** State the argument as numbered steps, give its margin,
and say what would overturn it. Then a per-locus table of which evidence
strands apply — position at a defined vector end, footprint matching the
modality, cross-sample scaling — and name the loci resting on a single strand.
Those are the validation priorities.

## 6. Is <N> the whole picture?

### Copy number

[[FIG:copy_number]]

Copies per haploid genome per construct and sample, from copy_number.py. Say
which route produced each number (measured from unique-block depth, boundary
estimator, or upper bound), the baseline, and the per-block spread. Quote the
shared-block closure shown under the figure: near zero means the additive model
holds; well below zero means shared sequence is under-represented (typically
backbone loss in concatemers), and the stacked solve is a diagnostic only. If
two trim policies are shown side by side, say whether they agree.

### Copies per site and the site-count bound

[[FIG:architecture]]

Per construct: the floor (certain + vector-end loci), the attributed count, the
copy-number ceiling and the unresolved candidates, from reconcile_sites.py, and
the copies per site that range implies. Quote the bounded range and the
unresolved count together.

### Reconciliation

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
that would overturn the current interpretation. Analytical steps (data to
request, reruns, re-mapping) go here; bench experiments go in section 12.

## 12. Validation priorities

Orthogonal wet-lab checks of the claims above, from `references/validation.md`.
A table, one row per check, in priority order:

| # | Check | Claim it tests | Design (amplicon / primers / platform) | Result that would overturn it |
|---|---|---|---|---|

Always consider, and say why any is omitted:

1. **ddPCR copy number** — one amplicon in each construct's *unique* region
   (coordinates from identity.json), against a validated single-copy reference.
   Tests the measured copy numbers, their ratio, and so every copy-ceiling
   attribution. Add a backbone amplicon for carryover.
2. **Junction PCR + Sanger** — one per `inferred` and every `unresolved` locus
   from called_sites.tsv, vector primer in construct-unique sequence, plus the
   highest-support `certain` locus as a positive control. Each one settles a
   single locus's attribution and confirms its footprint.
3. **Long-read sequencing** (ONT/PacBio, or targeted Cas9 enrichment) — resolves
   array length, copies per site, internal rearrangement and which construct
   sits at each shared-backbone locus at once. The definitive experiment when
   the site range or the closure deficit leaves the architecture open.

Name the loci resting on a single evidence strand (section 5) first: they are
the cheapest claims to overturn.
