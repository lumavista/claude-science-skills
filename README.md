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

## Compute: GPU first, CPU fallback

Every skill was audited for whether a GPU genuinely helps. Where it does, the
skill selects a GPU automatically and falls back to CPU; where it does not, the
skill says so and why. **CPU always produces a correct result** — the GPU is an
accelerator, never a requirement.

`SCI_DEVICE` overrides autodetection everywhere (`cpu` forces a reproducible CPU
run, `cuda` demands the GPU). It is honored verbatim rather than validated, so a
forced device fails loudly. `molecular-dynamics` also accepts `SCI_PLATFORM` and
`SCI_PRECISION`. Each skill carries its own copy of the selector — they import
independently, so the duplication is deliberate;
`skills/esm/scripts/_device.py` is canonical and `grep -rl SCI_DEVICE skills/`
lists the copies.

Reference host: Windows 11, **an 8 GiB consumer GPU**.

| Skill | GPU route | Notes |
|---|---|---|
| `esm` | ESM3 / ESMC / ESMFold2 via torch | `scripts/_device.py`; bf16 on Ampere+ |
| `molecular-dynamics` | OpenMM CUDA platform | `scripts/select_platform.py`; mixed precision |
| `scanpy` | rapids-singlecell: PCA, neighbors, UMAP, Leiden | **WSL2 only** — see below |
| `construct-assembly` | ThermoMPNN when scanning many mutants | CPU fine for a handful |
| `immunoinformatics` | routing only | `references/compute-routes.md` |
| `tamarind` | cloud GPU | for jobs exceeding the ~6 GiB usable budget |

### Four constraints worth knowing before you debug

**The default PyPI torch wheel on Windows is CPU-only.** Everything works and is
merely slow — the hardest failure to notice. `torch.cuda.is_available()` being
`False` here is almost always the wheel, not the hardware; install from the
`download.pytorch.org/whl/cu###` index.

**On native Windows, wheel availability gates more than VRAM does.** The card
looks adequate and the failure is an install error, not an OOM:

| Package | Windows wheel? | What it gates |
|---|---|---|
| `OpenMM-CUDA-12` / `-13` | **yes** | `molecular-dynamics` — works natively |
| `cupy-cuda12x` | yes | CuPy arrays |
| `rapids-singlecell` | **no** (sdist only) | `scanpy`'s GPU path |
| `jax-cuda12/13-plugin` | **no** (manylinux only) | AlphaFold2/ColabFold — no CUDA JAX at all |
| `flash-attn`, `deepspeed` | **no** (sdist only) | ESMFold2 fused kernels, scGPT, evo2 |
| `triton` | **no** upstream | same (a community `triton-windows` fork exists) |

So `scanpy`'s GPU path never activates on native Windows, and local folding is
blocked for a platform reason rather than a memory one. Also: RAPIDS managed
memory cannot oversubscribe VRAM under Windows or WSL, so an oversized run OOMs
rather than spilling to host RAM — hence scanpy's free-VRAM guard, which diverts
a too-large matrix to the CPU.

**WSL2 is not a fix when these skills run under Claude Science.** The app
executes skill code in a native Windows sandbox (ACL-granted `C:\`
roots, private desktop, loopback egress proxy). WSL is not a granted root, and
shelling into the WSL VM would leave that sandbox — which is what it exists to
prevent. Treat WSL2 as a route for manual work; in-app, use the hosted clients,
free notebooks, structure connector, or a cloud GPU route.

**8 GiB binds differently per workload.** For language models it is the weights:
ESMC-300M/600M and `esm3-sm-open-v1` fit in bf16; ESMC-6B (~12 GiB of weights
alone) does not and must use the hosted client. For folding it is chain length,
through activations rather than the checkpoint — fold a protomer before a
trimer. For MD, VRAM is rarely the limit. `fits_in_vram()` applies the sizing
arithmetic; it is a planning heuristic, not an allocation guarantee.

### GPU results are not bitwise identical to CPU results

Mainly in `scanpy`: cuML's UMAP and cuGraph's Leiden are separate
implementations a shared seed does not reconcile, and cuVS neighbors is exact
where scanpy's default is approximate. Cluster IDs can differ between backends at
the same resolution, so **annotations keyed to cluster IDs do not survive a
backend switch** — pin `SCI_DEVICE` and record which backend produced a result.

### Where a GPU deliberately is not used

Fourteen skills are CPU-bound; the reasons are not interchangeable. The eight
where a reader might plausibly expect acceleration say why in their own
`SKILL.md`.

- **Network-bound:** `paper-lookup`, `adaptyv`, `pathogen-variant-surveillance`
  (aggregation happens server-side at GenSpectrum)
- **IO / format-bound:** `anndata` (HDF5/Zarr codecs), `biopython`
- **External non-GPU executable:** `phylogenetics` — neither FastTree nor
  IQ-TREE has a CUDA path. BEAGLE's CUDA kernels serve Bayesian MCMC in
  BEAST/MrBayes, a different inference method, so using them would replace the
  method rather than accelerate it.
- **Already fast, or no GPU backend exists:** `pathway-enrichment` (GSEApy's
  permutation core is a multi-threaded Rust extension), `pydeseq2` (no GPU
  DESeq2 in Python), `scientific-visualization` (Matplotlib rasterises via Agg)
- **Too small to pay off:** `conservation-mapping`, `glycoengineering`,
  `experimental-design`, `scientific-writing`
- **Unbatchable by design:** `statistical-power` — the simulation harness
  deliberately calls your exact statsmodels model per replicate
- **Not applicable:** `vector-integration-analysis` never aligns reads; it
  consumes provider junction calls and `samtools depth` output

Three have a real non-GPU speedup left opt-in rather than changed silently:
`pydeseq2` defaults `--n-cpus 1` and `pathway-enrichment` defaults `--threads 4`
(raising both toward the core count is the entire win), and `statistical-power`
benefits from numpy vectorisation in closed-form cases and multiprocessing
across replicates.

### Verification status

CPU fallback paths were executed and behave correctly: with no GPU libraries
installed, selectors return CPU, helpers return `None` instead of raising,
`SCI_DEVICE` is honored, and `SCI_DEVICE=cuda` fails loudly with an install
hint. All 65 Python files compile and all 21 frontmatter blocks parse.

The CUDA paths were measured on the reference host (8 GiB consumer GPU, 8 GiB, Windows,
2026-10-17). These numbers are measured; anything not listed is still arithmetic.

| Measured | Result |
|---|---|
| MD, ~23k-atom solvated lysozyme, NPT, `Precision=mixed` | **a much higher rate CUDA vs 13.6 CPU — 23.6×** |
| fp32 matmul (4096³) | a higher throughput vs 0.55 CPU — 14.9× |
| tf32 / fp16 / bf16 matmul | bf16 several times fp32 |
| Usable VRAM | ~6 GiB working set (7074 of 8 GiB free; desktop holds about a GiB) |
| OpenMM platforms | 3 registered raw, 4 with the DLL fix, 0 plugin failures |
| torch CUDA | a CUDA build, , `is_available()` true inside the sandbox |
| `torch.compile` | unavailable — Triton has no Windows wheel |

A ns/day figure means nothing without its system, hence the atom count and
ensemble; that run held 300 K and 1.03 g/mL, so it is physically sound and not
merely fast. Overnight it is ~157 ns in 12 h against 6.7 ns on the CPU. It also
needs a quiet machine — an otherwise identical pass taken while a 2 GB download
ran reported a measurably lower rate with a drifting speed column. **bf16 is 4.2× fp32 here**,
which is why `autocast_dtype()` prefers it on Ampere and newer: size a model
against a higher throughput, not 8.3.

Two measured findings changed how this repo is written, both failures that never
raise an exception:

1. **OpenMM's CUDA platform needs the DLL fix** or it never registers and OpenMM
   quietly serves OpenCL instead, so `select_platform` must be imported before
   `openmm`.
2. **`torch.cuda.OutOfMemoryError` never fires on Windows.** an oversized allocation was
   allocated while `mem_get_info` reported 0 free, because WDDM pages GPU
   allocations into host RAM rather than failing. What breaks is throughput: the
   matmul held full speed near budget, then far slower past it — a 19–50×
   silent cliff. So `fits_in_vram()` is the only guard available, and
   `try/except OutOfMemoryError` is not a usable fallback here. RAPIDS managed
   memory fails the opposite way, with a real OOM; the two are not
   interchangeable.

Still unexecuted: the `esm` checkpoints themselves, the scanpy GPU backend
(unavailable on native Windows), and the per-checkpoint VRAM table.

## Original skills

**`immunoinformatics`** — MHC class I/II binding, B-cell epitope prediction,
IEDB cross-referencing, conservation across circulating variants, self-similarity
screening against the human proteome, HLA population coverage.
`scripts/epitope_triage.py` encodes the step ordering and writes a `.meta.json`
sidecar per output so each result carries its tool versions.
`scripts/mosaic_design.py` adds mosaic/cocktail construct design by k-mer
coverage maximisation plus a junctional neo-epitope check; mosaic mode walks a de
Bruijn graph built from the strain panel, so no k-mer absent from the panel can
be constructed (the LANL Mosaic Vaccine Designer and EpiGraph could not be
confirmed live as of 2026-10, so the algorithm is implemented directly).
`references/tool-licensing.md` records versions and licences for every predictor,
folding model and design tool in the collection, including code-versus-weights
splits, and corrects NetMHCpan to 4.2e. `references/compute-routes.md` is the
collection's canonical compute-routing doc.

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
dependency-free local developability screen (free cysteines, designed-in
sequons, deamidation and isomerisation motifs, hydrophobic patches, GRAVY, pI).
Inverse folding, de novo design and refolding are delegated to `tamarind` so they
cannot drift out of sync with its live schema discovery.

**`vector-integration-analysis`** — Characterizes integrated vector sequence in
stable cell lines (how many integration sites, where, and copy number) from
whole-genome or targeted insertion-site sequencing.

## Local modifications

**`tamarind`** — Scope narrowed to de novo, binder, and antibody design and
developability, with fewer trigger keywords and an added `## Scope` section. As
shipped upstream it competed with Claude Science featured skills for plain
structure prediction, docking and embeddings (AlphaFold2, Boltz, Chai-1,
ESMFold2, OpenFold3, DiffDock, ESM-2); it now defers to those and submits a cloud
job only under a stated condition — explicit request, insufficient local VRAM, or
batch size. Also records Tamarind's account-free key provisioning (a single POST
to `/api/agent/provision`), since the missing credential is what otherwise makes
the cloud tier unreachable.

**`phylogenetics`** — Made cross-platform. Upstream assumes MAFFT and IQ-TREE
executables; the skill now prefers wheel-based aligners (`pyfamsa` for MSA,
`pyhmmer` for profile alignment, `parasail`/Biopython for pairwise) and runs end
to end with FastTree where IQ-TREE is unavailable, stating what is lost
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
under **Settings > Specialists**. Re-importing this repo clears each specialist's
skill selections, so re-select them afterward.

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
   library docs belongs in a reference or behind a link — it is also the main
   maintenance liability, since it goes stale silently.
2. **Caveats are one sentence.** Correctness, determinism and "what was not
   tested" notes stay, compressed to the claim. They are the point of this
   collection; do not delete them to save space.
3. **Every reference file is cited by its SKILL.md.** There are no orphans; a
   reference nobody links to is dead weight.
4. **No tutorial prose.** The reader already knows Python and the science.

Two skills deliberately sit above the trend: `molecular-dynamics` (its four
simulation stages are distinct functions sharing one simulation object, not the
copy-paste they resemble) and `scientific-writing` (stopped at −21% because the
remainder is normative safety rules, where cutting means deleting rules).

## License

Upstream skills are MIT
([K-Dense-AI/scientific-agent-skills](https://github.com/K-Dense-AI/scientific-agent-skills)).
Each skill, original ones included, states its own license in its `SKILL.md`.
