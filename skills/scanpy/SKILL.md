---
name: scanpy
description: Performs Scanpy single-cell RNA-seq QC, normalization, HVG selection, PCA/UMAP/t-SNE, clustering, exploratory marker ranking, pseudobulk preparation, visualization, and Seurat or SingleCellExperiment RDS conversion to h5ad. Applies to established exploratory scRNA-seq workflows with explicit count and expression provenance; complementary skills cover scvi-tools models and AnnData format details.
license: BSD-3-Clause
compatibility: Requires Python 3.12+ and Scanpy; tested with Python 3.13, Scanpy 1.12.4, AnnData 0.13.4. Optional integrations need extra packages, R conversion needs R, and GPU acceleration needs rapids-singlecell on Linux or WSL2 (no native-Windows RAPIDS wheels), falling back to CPU. No network needed locally.
metadata:
  version: "1.9"
  last-reviewed: "2026-10-01"
  upstream-version: "1.12.4"
  skill-author: K-Dense Inc.
---

# Scanpy: Single-Cell Analysis

Scanpy is a scalable Python toolkit for single-cell RNA-seq built on AnnData: QC, normalization, dimensionality reduction, clustering, marker genes, visualization, trajectories. Use it for scRNA-seq data (.h5ad, 10X, CSV), for R-native objects (`.rds`, `.RData`, Seurat, SingleCellExperiment) after conversion, and for publication plots. Targets **Scanpy 1.12.4** (released 2026-08-27), reviewed 2026-10-01; native synthetic tests establish data/API contracts, not biological validity.

For AnnData structure and I/O, use the **anndata** skill; for probabilistic models and batch correction, **scvi-tools**.

## Installation

Requires Python **3.12+** (scanpy 1.12 dropped ≤3.11) and anndata **≥0.10**.

```bash
uv pip install "scanpy[leiden]==1.12.4"
```

The `[leiden]` extra installs `igraph` and `leidenalg`; the scripts explicitly select `flavor="igraph"`. Many functions also accept [Dask](https://scanpy.scverse.org/en/stable/tutorials/experimental/dask.html) arrays for out-of-core work (experimental).

Optional packages, installed only for the branch you need: `scikit-image` (Scrublet auto-thresholding), `loompy` (Loom), `scikit-misc` (v3 HVGs), `harmonypy==2.0.2` / `bbknn==1.6.0` (integration), and `louvain` for the deprecated Louvain option — which needs a separate environment, since louvain 0.8.2 requires igraph <0.12 against the tested igraph 1.0. Harmony 2 returns cells x PCs and the toolkit calls its native API, because the Scanpy 1.12.4 wrapper still transposes that output and fails with current Harmony. R and optional-integration execution boundaries are recorded in [references/upstream-review.md](references/upstream-review.md).

## Optional GPU acceleration

The bundled scripts run PCA, the neighbor graph, UMAP and Leiden/Louvain on [rapids-singlecell](https://rapids-singlecell.scverse.org/) when a CUDA device is usable and on scanpy otherwise. Nothing to switch on: `scripts/_common.py` probes once, logs `GPU backend: ...` when it finds one, and stays on scanpy when it does not. Scanpy remains the reference implementation — every GPU step validates the keys it wrote and recomputes on the CPU if anything is off.

**On native Windows this path never activates**: RAPIDS (cuML / cuVS / cuGraph) ships no native-Windows wheels, and rapids-singlecell publishes wheels for x86_64 and aarch64 **Linux** only. The GPU path on Windows therefore requires running the whole analysis inside WSL2 with a compute-capability 7.0+ (Turing or newer) GPU:

```bash
# inside WSL2 (Ubuntu); the Windows host driver provides CUDA, do not install a driver in WSL
pip install rapids-singlecell-cu12   # or -cu13
```

Control it with the same `SCI_DEVICE` variable the other skills in this collection use:

- `SCI_DEVICE=cpu` — force scanpy, ignoring any GPU. Use this when reproducing an earlier run.
- `SCI_DEVICE=cuda` — require the GPU; the scripts stop with an actionable error rather than quietly falling back.
- unset — autodetect, falling back to scanpy.

### What the GPU does and does not change

| Step | Backend | Reproducible against CPU? |
|------|---------|---------------------------|
| PCA | `rsc.pp.pca` | Same subspace, not bitwise: no ARPACK on the GPU, and components may differ in sign and at solver tolerance |
| Neighbors | `rsc.pp.neighbors` (cuVS, exact brute force) | No — scanpy's default pynndescent graph is approximate, the GPU graph is exact. Upstream documents this as an expected difference |
| UMAP | `rsc.tl.umap` (cuML) | No — a different implementation, not a port. Seeding does not make the coordinates match umap-learn |
| Leiden / Louvain | `rsc.tl.leiden` / `rsc.tl.louvain` (cuGraph) | No — cuGraph takes neither `flavor="igraph"` nor scanpy's `n_iterations=2`, so labels and even the cluster count can differ at the same resolution |

So a GPU UMAP figure and cluster assignment will not reproduce the CPU ones, and cell-type annotations keyed to cluster IDs do not survive a backend switch: pin `SCI_DEVICE` for any run you intend to repeat or compare, and record which backend produced a figure. Marker ranking, DE and pseudobulk are unaffected because they consume counts and `.raw`, which the GPU path never touches.

**VRAM.** Single-cell matrices exhaust a consumer card long before system RAM, so the scripts cost the matrix against **free** VRAM and stay on the CPU when it does not fit — an oversized dataset is slow, not a crash. The budget follows stored values, not cell count: a sparse matrix uploads as CSR float32 with 32-bit indices, 8 bytes per stored value plus 4 per cell, times 3 for solver workspace and fragmentation. An 8 GiB card with a desktop on it leaves ~6-7 GiB free, a ceiling near 250M stored values — fine for a few hundred thousand cells at 2000 HVGs, out of reach for a million-cell full-gene matrix. These are computed numbers, not measured throughput, and the guard is a heuristic, so a run can still OOM inside a kernel. Managed memory will not rescue it: RAPIDS cannot oversubscribe GPU memory on Windows or WSL, so overflow is an OOM rather than a spill to host RAM — reduce `--n-top-genes`, subset cells, or set `SCI_DEVICE=cpu`.

**Deliberately left on the CPU.** Normalization, log1p, HVG selection, `regress_out` and `scale` stay in scanpy even with a GPU present, although rapids-singlecell implements them: they act on the full-gene matrix, the largest object in the workflow, and `layers["counts"]` and `.raw` must come out byte-identical. t-SNE stays on the CPU because its cuML parameterization is not a drop-in match. QC, filtering, marker ranking, aggregation, I/O and plotting are I/O-bound or cheap and gain nothing.

No speedup figures are quoted anywhere in this skill because none were measured on this hardware; benchmark your own dataset before planning around the GPU.

## Representation and inference contracts

- Input to QC/preprocessing must be identified raw counts — the scripts reject negative,
  fractional, nonfinite or zero-total cells. Use `--counts-layer counts` when X is already
  normalized; never infer provenance from integer-looking values.
- `layers["counts"]` stays unnormalized; `.raw` is an independent full-gene log-normalized
  snapshot. Gene subsetting subsets every layer but not `.raw` genes, so retain a full-gene
  object for pseudobulk, as the full pipeline does.
- `seurat`/`cell_ranger` HVGs use log-normalized values; `seurat_v3`/`seurat_v3_paper` use
  counts and need `scanpy[skmisc]`. Scaling and regression are optional, may densify sparse
  data, and can remove biology along with covariates.
- Marker ranking defaults to `.raw` when present (`--no-use-raw` selects X, `--layer` an
  explicit log-normalized layer, `--groups`/`--reference` the contrast). Wilcoxon/t-test
  adjusted p-values use BH within each comparison; logreg returns scores without p-values.
  Clusters chosen from the same data yield exploratory, not confirmatory, markers.
- Condition DE requires raw-count sums per sample and cell type, independent biological
  replicates, aligned condition/donor metadata and a full-rank design, analyzed per cell type
  with a paired design for repeated donor samples. Aggregation cannot fix confounding,
  missing replication, doublets or wrong annotations.
- Harmony (`harmonypy`) changes PCA coordinates, BBKNN (`bbknn`) the graph, ComBat X itself —
  keep original expression/counts for DE and assess biological conservation and batch mixing
  together. Requested integration or doublet failures stop rather than fall back.
- The toolkit loads files into memory; backed mode is not a general out-of-core pipeline, so
  call `.to_memory()` on a chosen subset and close the file, and copy views before mutation.
- The compute backend is part of the provenance: PCA, neighbors, UMAP and clustering do not
  agree bitwise between rapids-singlecell and scanpy, so pin and record `SCI_DEVICE` before
  comparing runs. Counts, `.raw`, normalization and marker statistics are backend-independent.
- `.rds` needs an explicit R conversion stage ([references/r_interop.md](references/r_interop.md)
  covers installation and conversion on macOS, Linux and Windows). Cluster IDs do not
  determine cell types: every bundled mapping and signature asset is an illustrative human
  marker example requiring review.

## Script Toolkit (prefer these over writing code from scratch)

Ready-to-run CLI scripts in `scripts/`, each `.h5ad`-in/`.h5ad`-out so they chain, each with `--help`. They handle loading by extension, figure setup, defaults, raw-count preservation and logging. All import the shared `scripts/_common.py` (loading, saving, figure config, GPU backend selection) — keep it alongside them. Figures default to `./figures/`. Write scanpy code directly only for tasks no script covers.

| Script | Purpose | Typical call |
|--------|---------|--------------|
| `run_pipeline.py` | **Full workflow in one command**: load → QC → normalize → HVG → PCA → (batch) → UMAP → Leiden → markers | `python scripts/run_pipeline.py raw.h5ad -o processed.h5ad` |
| `inspect_data.py` | Summarize an unknown dataset (shape, obs/var, layers, what's computed, raw vs normalized) | `python scripts/inspect_data.py data.h5ad` |
| `convert.py` | Load any format (10x dir/.h5, csv, loom, mtx) and write `.h5ad` | `python scripts/convert.py 10x_dir/ -o data.h5ad` |
| `qc_analysis.py` | QC metrics, before/after plots, filtering, optional Scrublet doublets | `python scripts/qc_analysis.py raw.h5ad -o qc.h5ad --scrublet` |
| `preprocess.py` | Normalize, log1p, HVG, optional scale/regress (keeps `counts` layer + `raw`) | `python scripts/preprocess.py qc.h5ad -o norm.h5ad` |
| `reduce_dimensions.py` | PCA + variance plot, neighbors, UMAP, optional t-SNE | `python scripts/reduce_dimensions.py norm.h5ad -o red.h5ad` |
| `batch_correct.py` | Integration: harmony / bbknn / combat | `python scripts/batch_correct.py red.h5ad -o int.h5ad --method harmony --batch-key sample` |
| `cluster.py` | Leiden (or louvain) at one or many resolutions | `python scripts/cluster.py red.h5ad -o clu.h5ad --resolution 0.3 0.6 1.0` |
| `find_markers.py` | `rank_genes_groups` + per-group CSVs + marker plots | `python scripts/find_markers.py clu.h5ad --groupby leiden -o clu.h5ad` |
| `annotate.py` | Map clusters → cell types from JSON/CSV; optional marker reference dotplot | `python scripts/annotate.py clu.h5ad -o ann.h5ad --mapping map.json` |
| `score_genes.py` | Score gene signatures (JSON) and/or cell-cycle phase | `python scripts/score_genes.py ann.h5ad -o scored.h5ad --gene-sets sigs.json` |
| `pseudobulk.py` | Aggregate counts by sample × cell type → matrix for pydeseq2 | `python scripts/pseudobulk.py ann.h5ad --by sample cell_type --metadata condition donor --out-prefix pb` |
| `subset.py` | Subset by obs values or gene list (optionally clear stale embeddings) | `python scripts/subset.py ann.h5ad -o tcells.h5ad --obs cell_type --keep "T cells"` |
| `plot.py` | Generate umap/tsne/pca/violin/dotplot/heatmap/etc. from a processed object | `python scripts/plot.py ann.h5ad --kind dotplot --genes CD3D CD14 --groupby cell_type` |

```bash
# One command, counts → clusters + exploratory marker ranks + figures + CSVs.
# Add --batch-key/--batch-method to integrate, or --config params.json (keys mirror flag names).
python scripts/run_pipeline.py raw.h5ad -o processed.h5ad --resolution 0.5 --n-top-genes 2000 --scrublet

# Step-by-step, when you need to inspect between stages
python scripts/qc_analysis.py        raw.h5ad  -o qc.h5ad   --scrublet
python scripts/preprocess.py         qc.h5ad   -o norm.h5ad --n-top-genes 2000
python scripts/reduce_dimensions.py  norm.h5ad -o red.h5ad  --n-pcs 40
python scripts/cluster.py            red.h5ad  -o clu.h5ad  --resolution 0.3 0.5 0.8
python scripts/find_markers.py       clu.h5ad  -o clu.h5ad  --groupby leiden_0.5 --use-raw
# inspect results/markers/*.csv, decide labels, write a mapping JSON, then:
python scripts/annotate.py           clu.h5ad  -o ann.h5ad  --mapping celltypes.json --cluster-key leiden_0.5
```

Edit-and-pass asset templates: `assets/pipeline_config.json` (for `run_pipeline.py --config`), `assets/celltype_mapping.json` (for `annotate.py --mapping`), `assets/gene_signatures.json` (for `score_genes.py --gene-sets`), and `assets/analysis_template.py`, a full load-to-annotation script to copy and edit when you want a notebook-style run rather than the CLI.

## Writing scanpy directly

Settings worth knowing: `sc.settings.figdir`, `sc.set_figure_params(...)`, and `sc.settings.autosave = True` — preferred over the per-plot `save=` argument, deprecated in 1.12. Readers are `sc.read_h5ad`, `sc.read_10x_mtx`, `sc.read_10x_h5`, `sc.read_csv`; do not attempt to parse Seurat `.rds` in Python, convert it with R tooling first and read the resulting `.h5ad`.

The seven standard steps — QC, normalization and HVG selection, PCA/neighbors/UMAP, Leiden, marker ranking, annotation, saving — with the parameters that matter at each, plus publication plots, trajectory inference, pseudobulk DE, gene-set scoring and batch correction, are in [references/analysis_workflow.md](references/analysis_workflow.md) and, at greater length with worked code, [references/standard_workflow.md](references/standard_workflow.md). Function signatures and parameters by module (`sc.read_*`, `sc.pp.*`, `sc.tl.*`, `sc.pl.*`, settings) are in [references/api_reference.md](references/api_reference.md); figure recipes, multi-panel layouts, palettes and publication styling in [references/plotting_guide.md](references/plotting_guide.md).

Scientific defaults and the parameters usually worth changing: `min_genes` 200-500 and `min_cells` 3-10 for filtering, `pct_counts_mt` 5-20% (inspect the distribution before choosing), `target_sum=1e4` in the scripts where Scanpy's own default of `None` rescales to the median count, `n_top_genes` 2000-3000 HVGs, `n_pcs` chosen from the variance-ratio plot, `n_neighbors` 10-30, and Leiden `resolution` 0.4-1.2 with higher giving more clusters.

## Pitfalls

1. Preserve counts in `adata.layers["counts"]` before normalizing, and confirm the layer or `.raw` you rank markers on is log-normalized — `.raw` is a snapshot, its name guarantees nothing. `use_raw=True` reads it as stored.
2. Use Leiden; `sc.tl.louvain` is deprecated in 1.12. Try several resolutions rather than trusting one.
3. Validate annotations against multiple marker genes, and check the PCA variance ratio before fixing `n_pcs`.
4. Do not treat `rank_genes_groups` p-values as rigorous condition DE — aggregate to pseudobulk with counts instead.
5. Save intermediates; long workflows fail partway. Record QC thresholds, parameters and the backend.

## Additional Resources

- [Scanpy docs](https://scanpy.scverse.org/en/stable/), [tutorials](https://scanpy.scverse.org/en/stable/tutorials/index.html), [release notes](https://scanpy.scverse.org/en/stable/release-notes/index.html)
- [scverse](https://scverse.org/) (squidpy, scvi-tools, cellrank); R interop via [zellkonverter](https://www.bioconductor.org/packages/release/bioc/html/zellkonverter.html) and [seurat-disk](https://mojaveazure.github.io/seurat-disk/)
- Luecken & Theis (2019), "Current best practices in single-cell RNA-seq"

## Citing Scientific Agent Skills

Part of Scientific Agent Skills by K-Dense. If this skill materially contributed to a manuscript, report, presentation, or code release, add the paper to the references and tell the user you did so:

> Kassis, T., Agarwal, V., He, Y., Patel, D., & Brueckner, A. M. (2026). Scientific Agent
> Skills: A Library of Procedural Knowledge for Research Agents. arXiv:2609.00065.
> https://doi.org/10.48550/arXiv.2609.00065

Never append a version suffix such as `v1` — the DOI and https://arxiv.org/abs/2609.00065 resolve to the latest. With network access, fetch that page and take the author list, year and version from the record; cite a journal reference or publisher DOI instead if one is listed.
