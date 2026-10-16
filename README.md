# claude-science-skills

A curated set of Claude skills for virology, immunology, antigen & epitope
work, and vaccine design. Most are a trimmed subset of
[K-Dense-AI/scientific-agent-skills](https://github.com/K-Dense-AI/scientific-agent-skills)
(MIT), selected for these workflows; two skills are original to this repo and
two are locally modified.

The full upstream repo exceeds the 100 MB skill-import limit, mostly from
documentation screenshots and a demo GIF. This subset keeps only the skills
relevant to the workflows above, so it imports cleanly.

## Provenance

- **Original to this repo:** `immunoinformatics`, `conservation-mapping`, `construct-assembly`, `vector-integration-analysis`
- **From upstream, locally modified:** `tamarind`, `phylogenetics`
- **From upstream, unchanged:** the remaining 15 skills

Upstream is MIT-licensed. Each skill records its own license in the `license`
field of its `SKILL.md`.

## Skills

21 skills. Several reference each other by name — `immunoinformatics` hands off
to `pathogen-variant-surveillance` for circulating variants, to
`conservation-mapping` for per-position conservation, to `glycoengineering` for
glycan shielding, to `construct-assembly` for linkers and developability, and to
`scientific-writing` for the write-up — so removing one can degrade another's
workflow.

| Domain | Skill | What it does |
|---|---|---|
| Epitopes & antigen design | `immunoinformatics` | T-/B-cell epitope prediction and triage: MHC class I/II binding, conservation across variants, self-similarity screening, HLA population coverage |
| | `conservation-mapping` | Henikoff-weighted per-position conservation across a strain panel, reconciled to PDB numbering, written into B-factors, clustered into conserved surface patches |
| | `construct-assembly` | Linkers, tags, trimerisation and display domains, junction scanning, ThermoMPNN stability, and local developability screening for orderable constructs |
| | `glycoengineering` | Analyzes and engineers protein glycosylation — N-glycosylation sequons, glycan shielding, curated glycan databases |
| | `esm` | ESM3 protein generation, ESMC embeddings, and ESMFold2 folding via the Biohub ESM SDK |
| | `tamarind` | De novo, binder, and antibody/nanobody design and developability on the Tamarind Bio platform |
| | `adaptyv` | Designs and submits protein characterization experiments via the Adaptyv Bio Foundry API |
| Sequences & evolution | `biopython` | Sequence manipulation, FASTA/GenBank/PDB parsing, and programmatic NCBI/PubMed access |
| | `phylogenetics` | Multiple sequence alignment and phylogenetic tree inference |
| | `pathogen-variant-surveillance` | Genomic surveillance via GenSpectrum LAPIS — lineage nomenclature, sequence proportions, mutation profiles |
| Structure & dynamics | `molecular-dynamics` | Runs and analyzes MD simulations with OpenMM and MDAnalysis |
| Transcriptomics | `scanpy` | Single-cell RNA-seq QC, normalization, dimensionality reduction, clustering, and pseudobulk prep |
| | `anndata` | Annotated matrices (`.h5ad`, Zarr) for the scverse ecosystem |
| | `pydeseq2` | Bulk RNA-seq differential expression with PyDESeq2 |
| | `pathway-enrichment` | Pathway and gene-set enrichment analysis on gene lists or ranked data |
| Vector / integration | `vector-integration-analysis` | Characterizes integrated vector sequence in stable cell lines from WGS or targeted insertion-site sequencing |
| Literature & writing | `paper-lookup` | Multi-API scholarly search (18 sources) with provenance — papers, preprints, citations, OA full text |
| | `scientific-writing` | Drafts and audits manuscripts with evidence provenance, reporting-guideline coverage, and authorship tracking |
| | `scientific-visualization` | Truthful, accessible, publication-ready figures with Matplotlib, Seaborn, or Plotly |
| Study design | `experimental-design` | Chooses designs, randomizes, and blocks before data collection |
| | `statistical-power` | Sample-size and power calculations for study planning |

## Original skills

### `conservation-mapping`

Written from scratch for this repo. Scores per-position conservation with
Henikoff sequence weighting and reports effective sample size, so a
clade-skewed panel cannot masquerade as a conserved one. Reconciles alignment
positions to PDB residue numbering by aligning the structure's own CA-derived
sequence — the hazard `glycoengineering/references/glycoshield.md` flags
without fixing — and persists the result as a `pdb_resnum` column so every
downstream structural join has a reliable key. Clusters conserved exposed
residues into surface patches sized against a Fab footprint.

Alignment is delegated to `phylogenetics`; this skill does not align.

### `construct-assembly`

Written from scratch for this repo. Covers the parts of a construct that are
not the antigen: linker family and length, trimerisation and display domains,
tag placement, signal peptides, and junction scanning. Adds ThermoMPNN point-
mutant stability and a dependency-free local developability screen (free
cysteines, sequons created by design, deamidation and isomerisation motifs,
hydrophobic patches, GRAVY, pI).

Inverse folding, de novo design and refolding are delegated to `tamarind`
rather than reimplemented, so they cannot drift out of sync with its live
schema discovery.

### `immunoinformatics`

Written from scratch for this repo. Covers MHC class I/II binding prediction,
B-cell epitope prediction, IEDB cross-referencing, conservation across
circulating variants, self-similarity screening against the human proteome, and
HLA population coverage. Its `scripts/epitope_triage.py` encodes the step
ordering and writes a `.meta.json` sidecar per output so each result carries its
tool versions.

`scripts/mosaic_design.py` adds mosaic and cocktail construct design by k-mer
coverage maximisation, plus the junctional neo-epitope check the workflow has
always called for. Mosaic mode walks a de Bruijn graph built from the strain
panel, so no k-mer absent from the panel can be constructed. The LANL Mosaic
Vaccine Designer and EpiGraph could not be confirmed live as of 2026-10, so the
algorithm is implemented directly.

`references/tool-licensing.md` records current versions and licences for every
predictor, folding model and design tool the collection touches, including the
code-versus-weights splits that a single licence field would misrepresent — and
corrects NetMHCpan to 4.2e. `references/compute-routes.md` covers what runs on
CPU locally, what goes to Tamarind, and the free hosted GPU routes.

### `vector-integration-analysis`

Written from scratch for this repo. Characterizes integrated vector sequence in
stable cell lines — how many integration sites exist, where they are, and copy
number — from whole-genome or targeted insertion-site sequencing.

## Local modifications

### `tamarind`

Scope narrowed to de novo design, binder design, and antibody design and
developability, with trigger keywords reduced and a `## Scope` section added. As
shipped upstream it competed with the Claude Science featured skills for plain
structure prediction, docking, and embeddings (AlphaFold2, Boltz, Chai-1,
ESMFold2, OpenFold3, DiffDock, ESM-2). It now defers to those when they apply
and submits a cloud job only under a stated condition — explicit request, no
local GPU, or batch size.

### `phylogenetics`

Made cross-platform. Upstream assumes MAFFT and IQ-TREE executables, which are
not reliably available on every host. The skill now prefers wheel-based aligners
(`pyfamsa` for MSA, `pyhmmer` for profile alignment, `parasail`/Biopython for
pairwise work) and runs end to end with FastTree where IQ-TREE is unavailable,
stating what is lost when it falls back (ModelFinder, ultrafast bootstrap,
SH-aLRT) rather than letting a write-up imply a full ML analysis. Several
upstream defects found by running the code — a temp-file rename failure, a
silently skipped tree-collapse check, and an amino-acid alphabet missing J, U,
and O — were fixed, with regression tests added.

## Not included

Deliberately excluded because they overlap Claude Science built-ins, which do
the same job better. Skills that reference these by name resolve to the
built-ins.

- `scvi-tools` — Claude Science ships a featured scvi-tools skill (scVI, scANVI,
  totalVI, MultiVI).
- `literature-review` — Claude Science ships a featured literature-review skill.
  `paper-lookup` is kept because it is raw multi-API search with provenance, not
  review methodology, so it does not collide.
- `database-lookup` — overlaps the Claude Science featured connectors.

## Specialists

`specialists.yaml` documents the Claude Science specialists built on these
skills — each one's description, instructions, selected skills, and connectors.

| Specialist | Covers | Skills |
|---|---|---|
| **Claude Science** | Default profile — the full skill and connector set | all |
| **Immuno-Designer** | Epitope prediction, conservation, population coverage, candidate triage | `immunoinformatics`, `conservation-mapping`, `construct-assembly`, `glycoengineering`, `tamarind`, `biopython`, `phylogenetics`, `pathogen-variant-surveillance` |
| **Literature-Search** | Provenance-tracked multi-API search, evidence-traceable write-ups, study sizing | `paper-lookup`, `scientific-writing`, `scientific-visualization`, `experimental-design`, `statistical-power` |
| **Pathogen-Surveillance** | Circulating lineages, growth, mutations, phylogenies | `pathogen-variant-surveillance`, `phylogenetics`, `biopython`, `scientific-visualization` |
| **Immune-Profiler** | Systems-vaccinology transcriptomics: composition, DE, enrichment | `scanpy`, `anndata`, `pydeseq2`, `pathway-enrichment`, `statistical-power`, `scientific-visualization` |
| **Antigen-Structure** | Structure retrieval, epitope accessibility, glycan shielding, embeddings, light MD, cloud design | `tamarind`, `esm`, `glycoengineering`, `conservation-mapping`, `construct-assembly`, `molecular-dynamics`, `biopython`, `adaptyv`, `scientific-visualization` |

Claude Science has no specialist import/export, so `specialists.yaml` is
documentation, not something the app can load. Re-enter a specialist by hand
under **Settings > Specialists**. Re-importing this repo clears each
specialist's skill selections, so re-select them afterward.

## Usage

Import this repo into Claude Science as a skill source:

1. Clone or download the repo.
2. In Claude Science, add the repo (or its `skills/` folder) as a skill
   directory.
3. The skills become available for selection in the default profile and in
   per-specialist skill lists.

Each skill is self-contained under `skills/<name>/`, with a `SKILL.md` and its
own `references/`, `scripts/`, and `assets/` as needed. The skills follow the
standard Claude skill layout and also work anywhere else that loads
`SKILL.md`-based skills.

## License

Upstream skills are MIT ([K-Dense-AI/scientific-agent-skills](https://github.com/K-Dense-AI/scientific-agent-skills)).
Each skill states its own license in its `SKILL.md`. The two original skills
(`immunoinformatics`, `conservation-mapping`, `construct-assembly`,
`vector-integration-analysis`) carry their licenses in the
same field.
