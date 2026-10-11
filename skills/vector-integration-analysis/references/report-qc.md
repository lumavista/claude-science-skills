# Insertion-site deliverable — QC checklist

A provider's report is a summary built from pipeline output. The pipeline output
is the data. Work from the tables and treat the report as a claim to be
checked.

Run this at Stage 0, before analysis. None of the items below imply bad work by
a provider — they are the normal gaps between a standard pipeline and a
project-specific question, and most are fixed by asking.

## Blocking items — resolve before analysing

- [ ] **Genome build.** Often unstated. Without it no locus can be annotated to
      a gene, and builds differ by megabases.
- [ ] **Constructs in the mapping reference.** List them. If the line carries
      more constructs than the reference held, every shared-sequence junction
      is unattributable until the identity map is built, and ideally the run
      should be repeated against a combined reference.
- [ ] **Reference form for proviral or transposon delivery.** Was the reference
      the integrated form or the delivery plasmid? This shifts every expected
      junction (see `input-formats.md`).
- [ ] **Deduplication.** Confirm junction support is post-deduplication. If
      not, copy numbers are inflated several-fold.
- [ ] **Support threshold.** What cutoff separates the reported table from the
      full table? Get both.
- [ ] **Clip-side convention.** What do the clip values mean, and has it been
      confirmed against a junction of known orientation?

## Completeness

- [ ] **Does the report show all junctions?** Compare the stated detection
      count against the rows displayed. Interactive tables with a search box
      and a "show N entries" control routinely display a fraction of what was
      detected; the body text often says so, the widget does not.
- [ ] **How were the displayed junctions selected?** Often unstated, and often
      not by read support. If the criteria are unstated, the selection carries
      no meaning and the full table is the only usable input.
- [ ] **Lowest-support junctions.** Frequently the most informative: a
      single-copy insertion in a multi-copy background has the *fewest* reads.
      Check they are in the full table even when absent from the report.
- [ ] **Zero-byte or empty deliverables.** Files that would have provided
      evidence independent of split reads — read-pair support tables, for
      instance — are occasionally shipped empty. Check file sizes; ask for a
      rerun rather than concluding the evidence was absent.
- [ ] **Per-base depth.** Usually not delivered unless requested. Ask.

## Internal consistency

- [ ] **Counts in the summary against rows in the tables.** A mismatch usually
      means two different filters rather than an error — work out which, and
      which one you want.
- [ ] **Coverage figures against any low-coverage region file.** An empty
      low-coverage file contradicting visible gaps in a coverage plot usually
      means the plot is a rendering artefact, not that coverage is patchy.
- [ ] **Mean against median vector coverage.** A large gap means non-uniform
      coverage. In a multi-construct line this is the signature of
      cross-mapping and is quantitatively useful, not a defect.
- [ ] **Metrics reading zero or constant across every sample.** On-target
      percentages, enrichment factors and similar that are identical or zero
      everywhere are pipeline configuration problems; note them so nobody
      builds a conclusion on them.
- [ ] **Junction counts against sequencing depth.** If the deeper sample has
      proportionally more junctions, you are looking at sensitivity
      (`artefacts.md §7`).

## Clerical — cheap to check, expensive to miss

- [ ] Filename matches the construct inside the file.
- [ ] Table and figure captions name the right construct. A caption carrying a
      construct ID that is not yours is a template or copy-paste error; confirm
      which table the content belongs to rather than assuming.
- [ ] Table and figure numbering is consistent — repeated numbers, missing
      numbers, or boilerplate text referring to a different study all suggest a
      template was filled in by hand.
- [ ] Sample labels match your own sample sheet, in both directions.

## What to request

When asking for a rerun or supplementary files, ask for exactly this:

1. **A single combined reference** — host genome plus **every** construct in
   the line, including helper and donor plasmids — with alignments in shared
   regions flagged ambiguous rather than left to the aligner's mapping-quality
   tie-breaking, which punches false coverage holes through precisely the
   regions of interest.
2. **Per-base depth over each vector reference** (`samtools depth -a`). Tiny
   output, and it converts the copy-number estimate from bounded to measured.
3. **The genome build**, named explicitly.
4. **The full junction and breakpoint tables**, not the report subset, with the
   support threshold stated.
5. **Junction consensus or clipped sequence** per call, which is what makes the
   artefact confirming tests in `artefacts.md` possible.
6. **Both mappings** where two constructs share sequence and only one was
   mapped against — it is cheap and it bounds the answer.
