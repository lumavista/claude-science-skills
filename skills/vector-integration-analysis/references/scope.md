# Scope — when this is the wrong framework

This skill targets **clonal or oligoclonal engineered cell lines** characterised
by WGS or targeted insertion-site sequencing. A superficially similar body of
work — integration-site analysis for gene therapy — needs a different
framework, and applying either one to the other produces confident nonsense.

## Telling them apart

| | This skill | Clonal tracking |
|---|---|---|
| Material | clonal or oligoclonal line, a few to a few dozen sites | polyclonal population, 10³–10⁶ sites |
| Assay | WGS, or targeted junction capture | LAM-PCR / LM-PCR / linear-amplification protocols, usually with UMIs |
| Read support means | copy number × local depth | abundance of the clone carrying that site |
| Central question | how many sites, which construct, how many copies, is the cassette intact | which clones are expanding, is any site over-represented |
| Risk being assessed | cassette integrity, line stability, gene disruption | insertional mutagenesis, clonal dominance |

**The quickest discriminator:** are there UMIs, and is the sample a population
or a clone? UMIs plus a population means clonal tracking. No UMIs plus a clone
means this skill.

A **clonal** lentiviral line — a dCas9, shRNA or reporter line taken to a
single clone — belongs in **this** skill, with `--modality integrase`. It is
only the polyclonal, abundance-focused case that does not. The delivery route
being lentiviral is not what decides it; the question being asked is.

## Methods that do not transfer into this skill

**Relative abundance and the Gini coefficient.** These measure clonal dominance
across a population. In a clonal line there is one clone, so a Gini coefficient
is meaningless even though it computes. Without UMIs, read counts track copy
number and local depth, not abundance, so the metric is measuring the wrong
thing twice over.

**Common Integration Site and hotspot detection.** Sliding-window Poisson tests
with multiple-testing correction are built for thousands of sites across many
samples, to ask whether integration is non-random with respect to the genome.
With a few dozen sites in one clone the test is hopelessly underpowered, and
the question is wrong anyway — three sites megabases apart on one chromosome is
what chance looks like, not a hotspot. If you want to know whether integration
favours open chromatin, that needs a matched random-site background model and
far more sites.

**Network or pathway analysis of neighbouring genes.** Building an interaction
graph over genes near a couple of dozen loci yields hub genes that reflect the
database's degree distribution, not the integrations. It produces an attractive
figure with no statistical content. If gene context matters, annotate the loci
individually and say what each one is — that is honest and more useful. The
questions worth asking per locus are concrete: is the insertion inside a gene,
in which orientation, and is that gene essential or dosage-sensitive in this
cell type.

## Methods that do not transfer the other way

Do not take this skill's copy-number solve, copies-per-site estimator or site
accounting into a clonal-tracking study. They assume a single genotype shared
by all cells. In a population, a junction's read count reflects how many cells
carry it, and the arithmetic collapses.

## Also out of scope

- **Off-target editing analysis** for nuclease-edited lines — indels at
  predicted off-target sites, translocations between cut sites. Different
  question, different tools; this skill covers where the *vector* went, not
  where the nuclease cut.
- **Designing** vectors, cassettes or targeting constructs.
- **AAV** episomal persistence and genome configuration, which is dominated by
  concatemeric episomes rather than integration and needs its own treatment.
  The integration fraction of an AAV-transduced line can be analysed here, but
  the dominant species usually is not integrated at all, so start by
  establishing what fraction is.
- **Expression analysis.** Copy number is not expression; a silenced multi-copy
  array and a single active copy are indistinguishable by sequencing, and the
  difference is usually what matters. See `validation.md §5`.

## If the data turns out to be clonal tracking

Say so plainly and point at the established tools maintained by the gene-therapy
community rather than improvising. The statistics there are well developed and
this skill does not reimplement them.
