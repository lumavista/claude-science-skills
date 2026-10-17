# claude-science-skills

A curated set of Claude skills for virology, immunology, antigen & epitope
work, and vaccine design. Most are a trimmed subset of
[K-Dense-AI/scientific-agent-skills](https://github.com/K-Dense-AI/scientific-agent-skills)
(MIT), selected for these workflows; four skills are original to this repo and
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

## Compute: GPU first, CPU fallback

Every skill was audited for whether a GPU genuinely helps it. Where one does,
the skill now selects a GPU automatically and falls back to CPU; where one does
not, the skill says so and gives the reason, so the question does not get
re-asked. **CPU always produces a correct result** — the GPU is an accelerator,
never a requirement, and no skill fails because a GPU is absent.

Reference host: Windows 11, **an 8 GiB consumer GPU**.

### How to control it

`SCI_DEVICE` overrides autodetection everywhere: `SCI_DEVICE=cpu` forces a
reproducible CPU run, `SCI_DEVICE=cuda` demands the GPU. The value is honored
verbatim rather than validated, so a forced device fails loudly instead of
silently running somewhere you did not ask for. `molecular-dynamics` also
accepts `SCI_PLATFORM` and `SCI_PRECISION`.

Each skill carries its own copy of the small selector — they are imported
independently, so the duplication is deliberate. `skills/esm/scripts/_device.py`
is the canonical version; `grep -rl SCI_DEVICE skills/` lists the copies.

### Where the GPU is used

| Skill | GPU route | Notes |
|---|---|---|
| `esm` | ESM3 / ESMC / ESMFold2 via torch | `scripts/_device.py`; bf16 on Ampere+ |
| `molecular-dynamics` | OpenMM CUDA platform | `scripts/select_platform.py`; mixed precision |
| `scanpy` | rapids-singlecell: PCA, neighbors, UMAP, Leiden | **WSL2 only** — see below |
| `construct-assembly` | ThermoMPNN when scanning many mutants | CPU fine for a handful |
| `immunoinformatics` | routing only | `references/compute-routes.md` |
| `tamarind` | cloud GPU | for jobs exceeding 8 GiB |

### Three constraints worth knowing before you debug

**The default PyPI torch wheel on Windows is CPU-only.** Everything works and
is merely slow, which is the hardest failure to notice. `torch.cuda.is_available()`
returning `False` on this machine is almost always the wheel, not the hardware.
Install from the `download.pytorch.org/whl/cu###` index.

**On native Windows, wheel availability gates more than VRAM does.** This is the
constraint that surprises people, because the card looks adequate and the
failure is an install error rather than an out-of-memory:

| Package | Windows wheel? | What it gates |
|---|---|---|
| `OpenMM-CUDA-12` / `-13` | **yes** | `molecular-dynamics` — works natively |
| `cupy-cuda12x` | yes | CuPy arrays |
| `rapids-singlecell` | **no** (sdist only) | `scanpy`'s GPU path |
| `jax-cuda12/13-plugin` | **no** (manylinux only) | AlphaFold2/ColabFold — no CUDA JAX at all |
| `flash-attn`, `deepspeed` | **no** (sdist only) | ESMFold2 fused kernels, scGPT, evo2 |
| `triton` | **no** upstream | same (a community `triton-windows` fork exists) |

So `scanpy`'s GPU path never activates on native Windows, and local folding is
blocked for a platform reason rather than a memory one. One further Windows
specific: RAPIDS managed memory cannot oversubscribe VRAM under Windows or WSL,
so an oversized run OOMs rather than spilling to host RAM — hence scanpy's
free-VRAM guard, which diverts a too-large matrix to the CPU instead.

**WSL2 is not a fix when these skills run under Claude Science.** It is the
usual advice, and it is wrong in this context: the app executes skill code in a
native Windows sandbox, with ACL-granted `C:\` roots, a private
desktop and a loopback filtering proxy for egress. WSL is not among the granted
roots, and shelling into the WSL VM would leave that sandbox's filesystem and
egress policy behind — which is what the sandbox exists to prevent. So a CUDA
stack installed inside WSL is reachable from a terminal, but not from a skill
running in the app. Treat WSL2 as a route for manual work outside Claude
Science, not as a way to unblock a skill inside it. For in-app work the real
alternatives are the hosted clients, the free notebooks, the structure
connector, or a cloud GPU route.

**8 GiB is the binding constraint, and it binds differently per workload.** For
language models it is the weights: ESMC-300M/600M and `esm3-sm-open-v1` fit in
bf16, ESMC-6B (~12 GiB of weights alone) does not and must use the hosted
client. For folding it is chain length, through the activations rather than the
checkpoint — fold a protomer before a trimer. For MD, VRAM is rarely the limit
at all. `fits_in_vram()` applies the sizing arithmetic; it is a planning
heuristic, not an allocation guarantee.

### GPU results are not bitwise identical to CPU results

This matters most in `scanpy`: cuML's UMAP and cuGraph's Leiden are separate
implementations that a shared seed does not reconcile, and cuVS neighbors is
exact where scanpy's default is approximate. Cluster IDs can therefore differ
between backends at the same resolution, so **annotations keyed to cluster IDs
do not survive a backend switch** — pin `SCI_DEVICE` and record which backend
produced a result.

### Where a GPU deliberately is not used

Fourteen skills are CPU-bound. The eight where a reader might plausibly expect
GPU acceleration now say why in their own `SKILL.md`; the rest are listed here
only, because "it is an HTTP client" needs no note in the skill itself. The
reasons are not interchangeable:

- **Network-bound:** `paper-lookup`, `adaptyv`, `pathogen-variant-surveillance`
  (aggregation happens server-side at GenSpectrum)
- **IO / format-bound:** `anndata` (HDF5/Zarr codecs), `biopython`
- **External non-GPU executable:** `phylogenetics` — neither FastTree nor
  IQ-TREE has a CUDA path. BEAGLE's CUDA kernels serve Bayesian MCMC in
  BEAST/MrBayes, a different inference method, so using them would mean
  replacing the method rather than accelerating it.
- **Already fast, or no GPU backend exists:** `pathway-enrichment` (GSEApy's
  permutation core is a multi-threaded Rust extension), `pydeseq2` (no GPU
  DESeq2 equivalent exists in Python), `scientific-visualization` (Matplotlib
  rasterises through Agg)
- **Too small to pay off:** `conservation-mapping`, `glycoengineering`,
  `experimental-design`, `scientific-writing`
- **Unbatchable by design:** `statistical-power` — the simulation harness
  deliberately calls your exact statsmodels model per replicate
- **Not applicable:** `vector-integration-analysis` never aligns reads; it
  consumes provider junction calls and `samtools depth` output

Three of these have a real speedup available that is *not* a GPU, left as
opt-in rather than changed silently: `pydeseq2` defaults `--n-cpus 1` and
`pathway-enrichment` defaults `--threads 4` — raising both toward the core
count is the entire available win — while `statistical-power` benefits from
numpy vectorisation in the closed-form cases and multiprocessing across
replicates.

### Verification status

The CPU fallback paths were executed and behave correctly: with no GPU
libraries installed, the selectors return CPU, the helpers return `None`
instead of raising, `SCI_DEVICE` is honored, and `SCI_DEVICE=cuda` fails loudly
with an install hint. All 65 Python files compile and all 21 frontmatter blocks
parse. **The CUDA paths themselves have not been executed** — no GPU library
was installed at the time of writing — and no benchmark numbers or speedups are
quoted anywhere in this repo. The VRAM tables are arithmetic, not measurements.
Confirm on first real GPU run.

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
Each skill states its own license in its `SKILL.md`. The four original skills
(`immunoinformatics`, `conservation-mapping`, `construct-assembly`,
`vector-integration-analysis`) carry their licenses in the
same field.
