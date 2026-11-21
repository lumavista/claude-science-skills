<!--
HOW THIS REPORT IS BUILT. Fill this template; do not invent a structure.

* scripts/build_report.py builds the report from this template, and is the
  only thing that does:

      cp assets/report_template.md report.md      # write the prose into it
      python scripts/make_figures.py --run run.json
      python scripts/build_report.py --run run.json --report report.md --strict

  It writes report.html (the deliverable: self-contained, figures inline),
  report_filled.md, and report.docx only when an SVG rasteriser is present.
* You write ONLY the interpretive prose: replace the <...> text, keep the
  [[...]] placeholders exactly where they are. The script fills them:
    [[FIG:<id>]]    the skill's figures from make_figures.py, embedded as drawn
    [[TABLE:<id>]]  tables from the stage outputs (called sites, copy numbers,
                    site bounds, accounting, checks, caveats, ...)
    [[VAL:<id>]]    numbers and facts (site range, unresolved count, closure
                    deficit, host, modality, ...)
  Never retype a number a placeholder carries, never draw or paste in a
  figure (an image the skill did not draw is removed and fails the
  foreign_figure check), and never move a table into prose.
  `build_report.py --list-placeholders` lists every id.
* <construct A>, <construct B>, ... inside a placeholder are resolved by the
  script: A is the primary mapping, the rest follow identity.json's order. You
  may also write a construct's name, e.g. [[FIG:coverage_pAlpha]]. A slot whose
  input is missing gets a visible "not generated" note naming what is missing:
  leave it, or supply the input and rerun.
* Keep every section. One per construct for the per-construct results
  (sections 3-4: repeat for a third construct with <construct C>, delete
  section 3 for one). A section that does not apply says so in one line rather
  than disappearing.
-->

# Vector integration analysis — <line name>

**Project:** <project ID, pipeline and version>
**Samples:** [[VAL:samples]]. <what distinguishes them>
**Host:** [[VAL:host]]; ploidy [[VAL:ploidy]]. <cell line background, if more is known>
**Delivery modality:** [[VAL:modality]]. <whether the paperwork agrees>
**Constructs, in order of delivery:**

[[TABLE:constructs]]

1. <name, backbone, size, delivery route, defined ends, cargo, selection>
2. <...>
3. <any helper or donor plasmid, even if not expected to persist>

---

## Summary

[[TABLE:summary]]

| | [[VAL:name:<construct A>]] | [[VAL:name:<construct B>]] |
|---|---|---|
| Architecture | <...> | <...> |
| Backbone carryover | <...> | <...> |
| Cassette integrity | <...> | <...> |
| Effect of <the variable distinguishing the samples> | <...> | <...> |

**Site result.** [[VAL:site_range]].

**Headline result.** <two or three sentences>

**How the site count is derived.** [[VAL:site_derivation]]. <State any
alternative reading and the count it would give.>

**How copy number was obtained.** [[VAL:copy_number_route]]. Shared-block
closure: [[VAL:deficit]]. <If the route is not "measured", say what was
requested to make it so. If the stacked solve went negative, say what that
means for this line.>

**Main caveat.** Loci proven by sequence: [[VAL:n_certain]]; assigned by
argument (section 5): [[VAL:n_inferred]]; unresolved: [[VAL:n_unresolved]].
Confirm the last two before relying on them.

**Checks and caveats carried forward.**

[[TABLE:checks]]

[[TABLE:caveats]]

<For every check accepted with `--accept`: what was investigated and why it
does not invalidate the result.>

## 1. What was delivered, and how

[[TABLE:fingerprint]]

<The modality, the evidence for it from the junction fingerprint, and what
follows arithmetically. If the fingerprint disagreed with the paperwork, this
section is where that is resolved — not a footnote.>

## 2. The central problem: which construct did a read come from?

[[TABLE:identity]]

<Mechanism, what it did to the delivered data, how it was resolved. One worked
example of an artefactual junction is worth more than a paragraph of
explanation.>

Omit this section only if the line carries a single construct.

## 3. [[VAL:name:<construct B>]] — <result>

<Same structure as section 4, scaled to what this construct carries.>

[[FIG:coverage_<construct B>]]

**Copy number.** [[VAL:copy_number:<construct B>]].

[[TABLE:called_sites:<construct B>]]

<Read the track before quoting the copy number: the unique blocks should sit
flat at one level, the shared blocks at the sum of their carriers, and the
trimmed edges should be the only places the depth tapers. A step inside a
block, or a unique block at a different level from its sibling, is a deletion,
a rearrangement or a mis-specified reference, and belongs in section 7.>

## 4. [[VAL:name:<construct A>]] — [[VAL:sites:<construct A>]] sites

[[FIG:coverage_<construct A>]]

**Copy number.** [[VAL:copy_number:<construct A>]].

<Whether the track bears that number out, and whether the per-block estimates
agree.>

[[TABLE:called_sites:<construct A>]]

<The host-context columns come from called_sites.tsv when the host stage
annotated it; if they are absent, say why.>

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

Contigs with more than one called locus: [[VAL:multi_locus_contigs]].
<Any locus whose position (centromeric, subtelomeric, near a gene of known
consequence) matters for section 9.>

[[TABLE:called_sites]]

**Certain: [[VAL:n_certain]].** Vector side in sequence unique to one
construct.

**Inferred: [[VAL:n_inferred]].** <State the argument as numbered steps,
give its margin, and say what would overturn it. Name the loci resting on a
single evidence strand — position at a defined vector end, footprint matching
the modality, cross-sample scaling. Those are the validation priorities.>

**Unresolved: [[VAL:n_unresolved]].** Attributed to no construct; the
candidates are named.

[[TABLE:unresolved]]

## 6. Is the site range the whole picture?

### Copy number

[[FIG:copy_number]]

[[TABLE:copy_number]]

Shared-block closure: [[VAL:deficit]]. <Near zero means the additive model
holds; well below zero means shared sequence is under-represented (typically
backbone loss in concatemers), and the stacked solve is a diagnostic only. If
two trim policies are shown side by side, say whether they agree.>

### Copies per site and the site-count bound

[[FIG:architecture]]

[[TABLE:site_bounds]]

The ceiling divides by [[VAL:ceiling_divisor]]. Site result:
[[VAL:site_range]]. <What copies per site that range implies, and what the
upper bound rests on.>

### Reconciliation

[[TABLE:site_accounting]]

Flank recovery: [[VAL:flank_detection]]. <Calibration against a locus of known
structure. Why flank detection is below 100%. If there is a shortfall, what
would absorb it and how to test that.>

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

Junction-PCR priorities from called_sites.tsv:

[[TABLE:validation_loci]]
