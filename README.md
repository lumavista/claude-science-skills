# claude-science-skills

21 Claude skills for virology, immunology, antigen & epitope work, and vaccine
design. Most are a trimmed subset of
[K-Dense-AI/scientific-agent-skills](https://github.com/K-Dense-AI/scientific-agent-skills)
(MIT); the full upstream repo exceeds the 100 MB skill-import limit, so this
subset keeps only what these workflows need and imports cleanly.

- **Original to this repo:** `immunoinformatics`, `conservation-mapping`, `construct-assembly`, `vector-integration-analysis`
- **Locally modified:** `tamarind`, `phylogenetics`
- **Unchanged from upstream:** the remaining 15

Each skill states its own license in its `SKILL.md`.

## Usage

Clone the repo and add it (or its `skills/` folder) as a skill directory in
Claude Science; the skills then appear in the default profile and in
per-specialist skill lists. Each skill is self-contained under `skills/<name>/`
with a `SKILL.md` plus its own `references/`, `scripts/`, and `assets/`, and
works anywhere else that loads `SKILL.md`-based skills.

## Skills

Several reference each other by name — `immunoinformatics` hands off to
`pathogen-variant-surveillance`, `conservation-mapping`, `glycoengineering`,
`construct-assembly`, and `scientific-writing` — so removing one can degrade
another's workflow.

| Domain | Skill | What it does |
|---|---|---|
| Epitopes & antigen design | `immunoinformatics` | T-/B-cell epitope prediction and triage: MHC class I/II binding, conservation across variants, self-similarity screening, HLA population coverage |
| | `conservation-mapping` | Henikoff-weighted per-position conservation across a strain panel, reconciled to PDB numbering, written into B-factors, clustered into conserved surface patches |
| | `construct-assembly` | Linkers, tags, trimerisation and display domains, junction scanning, ThermoMPNN stability, local developability screening |
| | `glycoengineering` | N-glycosylation sequons, glycan shielding, curated glycan databases |
| | `esm` | ESM3 generation, ESMC embeddings, ESMFold2 folding via the Biohub ESM SDK |
| | `tamarind` | De novo, binder, and antibody/nanobody design and developability on Tamarind Bio |
| | `adaptyv` | Designs and submits protein characterization experiments via the Adaptyv Bio Foundry API |
| Sequences & evolution | `biopython` | Sequence manipulation, FASTA/GenBank/PDB parsing, programmatic NCBI/PubMed access |
| | `phylogenetics` | Multiple sequence alignment and phylogenetic tree inference |
| | `pathogen-variant-surveillance` | Genomic surveillance via GenSpectrum LAPIS — lineages, sequence proportions, mutation profiles |
| Structure & dynamics | `molecular-dynamics` | Runs and analyzes MD simulations with OpenMM and MDAnalysis |
| Transcriptomics | `scanpy` | scRNA-seq QC, normalization, dimensionality reduction, clustering, pseudobulk prep |
| | `anndata` | Annotated matrices (`.h5ad`, Zarr) for the scverse ecosystem |
| | `pydeseq2` | Bulk RNA-seq differential expression |
| | `pathway-enrichment` | Pathway and gene-set enrichment on gene lists or ranked data |
| Vector / integration | `vector-integration-analysis` | Characterizes integrated vector sequence in stable cell lines from WGS or targeted insertion-site sequencing |
| Literature & writing | `paper-lookup` | Multi-API scholarly search (18 sources) with provenance |
| | `scientific-writing` | Drafts and audits manuscripts with evidence provenance, reporting-guideline coverage, authorship tracking |
| | `scientific-visualization` | Truthful, accessible, publication-ready figures (Matplotlib, Seaborn, Plotly) |
| Study design | `experimental-design` | Chooses designs, randomizes, and blocks before data collection |
| | `statistical-power` | Sample-size and power calculations |

### Not included

Excluded because Claude Science built-ins do the same job better; skills that
reference them by name resolve to the built-ins: `scvi-tools`,
`literature-review`, `database-lookup`. `paper-lookup` is kept because it is raw
multi-API search with provenance, not review methodology.

## Compute

Every skill was audited for whether a GPU genuinely helps it. Where it does, the
skill selects one automatically and falls back to CPU; where it does not, the
skill says so and why, so the question does not get re-asked. **CPU always
produces a correct result** — the GPU is an accelerator, never a requirement,
and no skill fails because a GPU is absent.

`SCI_DEVICE` overrides autodetection everywhere: `cpu` forces a reproducible CPU
run, `cuda` demands a GPU. The value is honored verbatim rather than validated,
so a forced device fails loudly instead of quietly running somewhere else.
`molecular-dynamics` also accepts `SCI_PLATFORM` and `SCI_PRECISION`.

Six skills have a GPU route: `esm` (torch), `molecular-dynamics` (OpenMM CUDA),
`scanpy` (rapids-singlecell), `construct-assembly` (ThermoMPNN, for large mutant
scans), `tamarind` (cloud GPU), and `immunoinformatics` (routing only). Each
carries its own copy of the device selector — the skills import independently,
so the duplication is deliberate; `skills/esm/scripts/_device.py` is the
canonical version and `grep -rl SCI_DEVICE skills/` lists the copies.

GPU availability is a packaging question as much as a hardware one: several
accelerated stacks publish no wheel for every platform, so a GPU route can be
unavailable on an otherwise capable machine. The skills detect this and fall
back rather than failing. `skills/immunoinformatics/references/compute-routes.md`
is the canonical routing document — which steps run locally on CPU, which use a
local GPU, and which exceed it and go to a hosted or cloud route.

**GPU results are not bitwise identical to CPU results.** Mainly in `scanpy`:
cuML's UMAP and cuGraph's Leiden are separate implementations that a shared seed
does not reconcile, and cuVS neighbors is exact where scanpy's default is
approximate. Cluster IDs can differ between backends at the same resolution, so
annotations keyed to cluster IDs do not survive a backend switch — pin
`SCI_DEVICE` and record which backend produced a result.

The other fourteen skills are CPU-bound for reasons that are not
interchangeable — network-bound, IO-bound, dependent on a non-GPU external
executable, already fast, or too small to pay off — and each states its own.
Three have a real speedup available that is not a GPU, left opt-in rather than
changed silently: `pydeseq2` defaults to `--n-cpus 1` and `pathway-enrichment`
to `--threads 4`, and `statistical-power` benefits from vectorisation in the
closed-form cases and multiprocessing across replicates.

## Original skills

**`immunoinformatics`** — MHC class I/II binding, B-cell epitope prediction,
IEDB cross-referencing, conservation across circulating variants,
self-similarity screening against the human proteome, HLA population coverage.
`scripts/epitope_triage.py` encodes the step ordering and writes a `.meta.json`
sidecar per output, so every result carries its tool versions.
`scripts/mosaic_design.py` adds mosaic and cocktail construct design by k-mer
coverage maximisation plus a junctional neo-epitope check; mosaic mode walks a
de Bruijn graph built from the strain panel, so no k-mer absent from the panel
can be constructed. `references/tool-licensing.md` records versions and licences
for every predictor, folding model and design tool the collection touches,
including the code-versus-weights splits a single licence field would
misrepresent.

**`conservation-mapping`** — Per-position conservation with Henikoff sequence
weighting, reporting effective sample size so a clade-skewed panel cannot
masquerade as a conserved one. Reconciles alignment positions to PDB residue
numbering via the structure's own CA-derived sequence — the hazard
`glycoengineering/references/glycoshield.md` flags without fixing — and persists
a `pdb_resnum` column so downstream structural joins have a reliable key.
Clusters conserved exposed residues into patches sized against a Fab footprint.
Alignment is delegated to `phylogenetics`.

**`construct-assembly`** — The parts of a construct that are not the antigen:
linker family and length, trimerisation and display domains, tag placement,
signal peptides, junction scanning. Adds ThermoMPNN point-mutant stability and a
dependency-free developability screen (free cysteines, designed-in sequons,
deamidation and isomerisation motifs, hydrophobic patches, GRAVY, pI). Inverse
folding, de novo design and refolding are delegated to `tamarind` so they cannot
drift out of sync with its live schema discovery.

**`vector-integration-analysis`** — Characterizes integrated vector sequence in
stable cell lines — how many integration sites, where, and copy number — from
whole-genome or targeted insertion-site sequencing.

## Local modifications

**`tamarind`** — Scope narrowed to de novo, binder, and antibody design and
developability, with fewer trigger keywords and an added `## Scope` section. As
shipped upstream it competed with the Claude Science featured skills for plain
structure prediction, docking and embeddings (AlphaFold2, Boltz, Chai-1,
ESMFold2, OpenFold3, DiffDock, ESM-2); it now defers to those and submits a
cloud job only under a stated condition — explicit request, insufficient local
VRAM, or batch size. Also records Tamarind's account-free key provisioning,
since the missing credential is what otherwise makes the cloud tier unreachable.

**`phylogenetics`** — Made cross-platform. Upstream assumes MAFFT and IQ-TREE
executables; the skill now prefers wheel-based aligners (`pyfamsa` for MSA,
`pyhmmer` for profile alignment, `parasail`/Biopython for pairwise work) and
runs end to end with FastTree where IQ-TREE is unavailable, stating what is lost
(ModelFinder, ultrafast bootstrap, SH-aLRT) rather than letting a write-up imply
a full ML analysis. Several upstream defects found by running the code — a
temp-file rename failure, a silently skipped tree-collapse check, and an amino
acid alphabet missing J, U and O — were fixed, with regression tests added.

## Specialists

`specialists.yaml` documents the Claude Science specialists built on these
skills — descriptions, instructions, skill selections, connectors.

| Specialist | Covers | Skills |
|---|---|---|
| **Claude Science** | Default profile — full skill and connector set | all |
| **Immuno-Designer** | Epitope prediction, conservation, population coverage, candidate triage | `immunoinformatics`, `conservation-mapping`, `construct-assembly`, `glycoengineering`, `tamarind`, `biopython`, `phylogenetics`, `pathogen-variant-surveillance` |
| **Literature-Search** | Provenance-tracked search, evidence-traceable write-ups, study sizing | `paper-lookup`, `scientific-writing`, `scientific-visualization`, `experimental-design`, `statistical-power` |
| **Pathogen-Surveillance** | Circulating lineages, growth, mutations, phylogenies | `pathogen-variant-surveillance`, `phylogenetics`, `biopython`, `scientific-visualization` |
| **Immune-Profiler** | Systems-vaccinology transcriptomics: composition, DE, enrichment | `scanpy`, `anndata`, `pydeseq2`, `pathway-enrichment`, `statistical-power`, `scientific-visualization` |
| **Antigen-Structure** | Structure retrieval, epitope accessibility, glycan shielding, embeddings, light MD, cloud design | `tamarind`, `esm`, `glycoengineering`, `conservation-mapping`, `construct-assembly`, `molecular-dynamics`, `biopython`, `adaptyv`, `scientific-visualization` |

Claude Science has no specialist import/export, so `specialists.yaml` is
documentation, not something the app can load. Re-enter a specialist by hand
under **Settings > Specialists**. Re-importing this repo clears each
specialist's skill selections, so re-select them afterward.

## Conventions

The skills were trimmed in 2026-10 to carry what is required and little else:
SKILL.md bodies went from 318k to 234k characters (−26.6%). Three layers load at
different times, and the weight that matters is in the first two:

| Layer | Loads | Budget |
|---|---|---|
| `description`, `compatibility` | every session, for skill selection | ≤400 / ≤300 chars |
| SKILL.md body | when the skill activates | orient and point, do not explain |
| `references/`, `assets/` | on demand only | detail belongs here |

1. **SKILL.md orients; references hold detail.** A passage restating upstream
   library documentation belongs in a reference or behind a link — it is also
   the main maintenance liability, since it goes stale silently.
2. **Caveats are one sentence.** Correctness, determinism and "what was not
   tested" notes stay, compressed to the claim. They are the point of this
   collection; do not delete them to save space.
3. **Every reference file is cited by its SKILL.md.** There are no orphans; a
   reference nobody links to is dead weight.
4. **No tutorial prose.** The reader already knows Python and the science.

Two skills deliberately sit above the trend: `molecular-dynamics`, whose four
simulation stages are distinct functions sharing one simulation object rather
than the copy-paste they resemble, and `scientific-writing`, which stopped at
−21% because the remainder is normative safety rules, where cutting further
means deleting rules rather than compressing them.

## License

Upstream skills are MIT
([K-Dense-AI/scientific-agent-skills](https://github.com/K-Dense-AI/scientific-agent-skills)).
Each skill, the original ones included, states its own license in its
`SKILL.md`.
