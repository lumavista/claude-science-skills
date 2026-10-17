#!/usr/bin/env python3
"""
Shared helpers for the scanpy script toolkit.

Every CLI script in this directory imports from this module so that data
loading, saving, figure configuration, and logging behave consistently.
This file is NOT a CLI itself; import it:

    from _common import load_anndata, save_anndata, configure_scanpy, info

The "GPU backend" section below routes PCA, the neighbor graph, UMAP and
clustering through rapids-singlecell when a CUDA device is usable, and through
scanpy otherwise. Nothing else changes: scanpy remains the reference path.
"""

import os
import sys


def info(msg):
    """Print a progress message to stderr-friendly stdout with a marker."""
    print(f"[scanpy] {msg}", flush=True)


def die(msg, code=1):
    """Print an error and exit."""
    print(f"Error: {msg}", file=sys.stderr, flush=True)
    sys.exit(code)


def _import_scanpy():
    try:
        import scanpy as sc  # noqa: F401
        return sc
    except ImportError:
        die('scanpy not installed. Install with: uv pip install "scanpy[leiden]"')


def configure_scanpy(figdir="figures", dpi=120, verbosity=1, autosave=False,
                     file_format="png"):
    """Apply consistent scanpy settings and return the scanpy module.

    Note: autosave is left off by default because the toolkit scripts pass
    explicit ``save=`` suffixes to plotting calls for predictable filenames.
    """
    sc = _import_scanpy()
    sc.settings.verbosity = verbosity
    sc.set_figure_params(dpi=dpi, facecolor="white")
    sc.settings.figdir = figdir
    sc.settings.file_format_figs = file_format
    sc.settings.autosave = autosave
    sc.settings.autoshow = False
    os.makedirs(figdir, exist_ok=True)
    return sc


def load_anndata(path, var_names="gene_symbols"):
    """Load an AnnData object, dispatching on the file extension / layout.

    Supported inputs:
      * ``.h5ad``                  -> sc.read_h5ad
      * ``.h5`` (10x CellRanger)   -> sc.read_10x_h5
      * ``.csv`` / ``.tsv`` / ``.txt`` -> sc.read_csv / read_text
      * ``.loom``                  -> sc.read_loom
      * ``.mtx``                   -> sc.read (matrix market)
      * a directory                -> sc.read_10x_mtx (10x mtx folder)
    """
    sc = _import_scanpy()
    if not os.path.exists(path):
        die(f"input not found: {path}")

    path = os.fspath(path)
    lower = path.lower()
    if os.path.isdir(path):
        info(f"Reading 10x mtx directory: {path}")
        return sc.read_10x_mtx(path, var_names=var_names)
    if lower.endswith(".h5ad"):
        return sc.read_h5ad(path)
    if lower.endswith(".h5"):
        info("Reading 10x HDF5 (.h5)")
        return sc.read_10x_h5(path)
    if lower.endswith(".loom"):
        return sc.read_loom(path)
    if lower.endswith(".csv"):
        return sc.read_csv(path)
    if lower.endswith((".tsv", ".txt")):
        return sc.read_text(path)
    if lower.endswith(".mtx") or lower.endswith(".mtx.gz"):
        return sc.read_mtx(path)
    die(f"unrecognized input format: {path}")


def save_anndata(adata, path):
    """Write an AnnData object to .h5ad, creating parent dirs as needed."""
    parent = os.path.dirname(os.path.abspath(path))
    os.makedirs(parent, exist_ok=True)
    adata.write_h5ad(path)
    info(f"Wrote {path}  ({adata.n_obs} cells x {adata.n_vars} genes)")


def add_io_args(parser, default_output=None):
    """Attach the standard input/output/figdir arguments to an argparse parser."""
    parser.add_argument("input", help="Input file (.h5ad, .h5, .csv, .loom, or 10x mtx dir)")
    parser.add_argument("-o", "--output", default=default_output,
                        help="Output .h5ad path" +
                             (f" (default: {default_output})" if default_output else ""))
    parser.add_argument("--figdir", default="figures",
                        help="Directory for saved figures (default: figures)")
    return parser


def _named_keys(mapping):
    """Named keys of an AnnData mapping, in order.

    anndata >= 0.13 reports an unnamed `None` key on `.layers` standing for X
    itself. Joining that into a string raises TypeError, so filter it out.
    """
    return [key for key in mapping.keys() if isinstance(key, str)]


def summarize(adata):
    """Return a short human-readable summary string of an AnnData object."""
    lines = [f"{adata.n_obs} cells x {adata.n_vars} genes"]
    if len(adata.obs.columns):
        lines.append("obs: " + ", ".join(adata.obs.columns[:20]))
    obsm = _named_keys(adata.obsm)
    if obsm:
        lines.append("obsm: " + ", ".join(obsm))
    layers = _named_keys(adata.layers)
    if layers:
        lines.append("layers: " + ", ".join(layers))
    return "\n".join(lines)


def validate_counts(matrix, label="counts", require_nonzero_cells=True):
    """Validate finite, nonnegative integer counts without densifying sparse input.

    Values alone cannot establish provenance: the caller must identify the actual
    unnormalized assay. Floating dtypes are accepted when every value is integral.
    """
    import numpy as np
    from scipy import sparse

    if matrix is None or len(matrix.shape) != 2 or min(matrix.shape) == 0:
        die(f"{label}: expected a nonempty cells x genes count matrix")
    values = matrix.data if sparse.issparse(matrix) else np.asarray(matrix)
    if (not np.isfinite(values).all() or (values < 0).any()
            or not np.equal(values, np.floor(values)).all()):
        die(f"{label}: expected finite, nonnegative integer counts; verify assay provenance")
    if require_nonzero_cells and (np.asarray(matrix.sum(axis=1)).ravel() <= 0).any():
        die(f"{label}: contains zero-total cells; filter them before normalization")


def prepare_counts(adata, layer=None):
    """Select explicitly identified counts and preserve an independent copy."""
    import numpy as np
    from scipy import sparse

    if layer is not None and layer not in adata.layers:
        die(f"counts layer '{layer}' is missing")
    matrix = adata.X if layer is None else adata.layers[layer]
    validate_counts(matrix, "X" if layer is None else f"layers[{layer!r}]")
    if layer is None and "log1p" in adata.uns:
        die("X has log1p provenance; select the original --counts-layer explicitly")
    if layer is None and "counts" in adata.layers:
        old = adata.layers["counts"]
        delta = old - matrix
        differs = delta.nnz != 0 if sparse.issparse(delta) else np.any(delta != 0)
        if differs:
            die("X differs from the preserved counts layer; select --counts-layer counts")
    adata.X = matrix.copy()
    adata.layers["counts"] = matrix.copy()
    adata.uns.pop("log1p", None)


def normalize_hvg(sc, adata, *, target_sum=1e4, n_top_genes=2000,
                  flavor="seurat", batch_key=None):
    """Keep all genes/counts, select HVGs on the flavor's required representation."""
    if target_sum <= 0 or n_top_genes < 1:
        die("target_sum and n_top_genes must be positive")
    if batch_key is not None:
        ensure_categories(adata, batch_key)
    count_flavor = flavor in {"seurat_v3", "seurat_v3_paper"}
    if count_flavor:
        sc.pp.highly_variable_genes(adata, layer="counts", n_top_genes=n_top_genes,
                                    flavor=flavor, batch_key=batch_key)
    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)
    if not count_flavor:
        sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes,
                                    flavor=flavor, batch_key=batch_key)
    if not adata.var["highly_variable"].any():
        die("HVG selection returned no genes; inspect expression/filtering parameters")
    adata.raw = adata.copy()


def ensure_categories(adata, key):
    """Use observed, nonmissing categories; numeric labels become lexical strings."""
    if key not in adata.obs or adata.obs[key].isna().any():
        die(f"obs column '{key}' must exist and contain no missing labels")
    adata.obs[key] = adata.obs[key].astype(str).astype("category")


# --------------------------------------------------------------------------
# GPU backend
#
# Adapted from the canonical device snippet in skills/esm/scripts/_device.py,
# which selects a torch device; here the equivalent choice is a whole module.
# Same contract: scanpy on the CPU always produces a correct result, the GPU is
# only an accelerator, nothing raises when the GPU stack is absent, and an
# explicit SCI_DEVICE is honored rather than silently downgraded.
#
# Reality check before relying on this: RAPIDS (cuML/cuVS/cuGraph), which
# rapids-singlecell is built on, ships no native-Windows wheels. NVIDIA's
# install matrix lists Linux with glibc>=2.28 and "Windows 11 using a WSL2
# specific install (WSL1 not supported)"; rapids-singlecell ships prebuilt
# wheels for x86_64/aarch64 Linux only. On native Windows the import below
# fails and every caller runs on the CPU, by design.
# --------------------------------------------------------------------------

GPU_INSTALL_HINT = (
    "rapids-singlecell needs an NVIDIA GPU of compute capability 7.0+ (Turing or "
    "newer) and Linux, or Windows 11 via WSL2 -- RAPIDS has no native-Windows "
    "wheels. Install with: pip install rapids-singlecell-cu12"
)

_UNPROBED = object()
_backend_cache = _UNPROBED


def _probe_backend():
    """Import rapids-singlecell, but only once a CUDA device actually answers.

    CuPy is probed first because it is the cheap question: importing
    rapids-singlecell pulls in the whole cuML/cuVS/cuGraph stack. Every failure
    mode -- package absent (native Windows, CPU-only box), driver/runtime
    mismatch, no visible device -- returns None instead of propagating.
    """
    try:
        import cupy

        if cupy.cuda.runtime.getDeviceCount() < 1:
            return None
    except Exception:
        return None
    try:
        import rapids_singlecell as rsc
    except Exception:
        return None
    return rsc


def select_backend(prefer=None):
    """Return the rapids-singlecell module to compute with, or None for scanpy.

    rapids-singlecell reimplements the expensive scanpy steps -- PCA, the
    neighbor graph, UMAP, Leiden/Louvain, HVG selection, regress_out -- on CuPy
    and NVIDIA CUDA-X Data Science. This toolkit uses it for PCA, neighbors,
    UMAP and clustering only; see the module docstring of each helper.

    Priority: explicit argument > SCI_DEVICE env var > autodetection. A request
    that does not name CUDA ("cpu", "mps", ...) forces the scanpy path and
    returns None, which is the supported way to switch the GPU off. A request
    that *does* name CUDA is honored rather than quietly downgraded: if the
    backend is unusable this calls die() with an actionable message, matching
    _device.py's refusal to run somewhere the caller did not ask for.
    Autodetection itself never raises -- a missing package, a broken driver or
    a CUDA build mismatch all yield None and the caller proceeds on the CPU.

    The probe is cached because importing cuML costs seconds.
    """
    global _backend_cache
    requested = (prefer or os.environ.get("SCI_DEVICE") or "").strip().lower()
    if requested and not requested.startswith(("cuda", "gpu")):
        return None
    if _backend_cache is _UNPROBED:
        _backend_cache = _probe_backend()
        if _backend_cache is not None:
            info(f"GPU backend: {backend_report(_backend_cache)}")
    if requested and _backend_cache is None:
        die(f"SCI_DEVICE={requested} requested a GPU backend that is not usable. "
            + GPU_INSTALL_HINT)
    return _backend_cache


def vram_gb(free=False):
    """Device memory in GiB: total by default, free when ``free=True``.

    Total is the stable number to plan against; free is the number that decides
    whether the next allocation survives, so the size guards below use it.
    Returns None when there is no CUDA device to ask.
    """
    try:
        import cupy

        available, total = cupy.cuda.Device().mem_info
        return (available if free else total) / (1024**3)
    except Exception:
        return None


def matrix_gpu_gb(matrix):
    """GiB a cells x genes matrix occupies on the device as float32.

    Sparse input is costed as the CSR form `_to_device` uploads: float32 data
    plus 32-bit indices plus one indptr entry per cell. The 32-bit indices also
    cap a single device matrix at 2**31-1 stored values.
    """
    from scipy import sparse

    if sparse.issparse(matrix):
        return (int(matrix.nnz) * 8 + (matrix.shape[0] + 1) * 4) / (1024**3)
    return float(matrix.shape[0]) * float(matrix.shape[1]) * 4 / (1024**3)


def fits_in_vram(matrix, *, overhead=3.0, reserve=0.0):
    """Rough check that a matrix plus workspace fits. None when memory is unknown.

    `overhead` covers the solver's scratch buffers and allocator fragmentation;
    `reserve` adds GiB for a result the caller already knows the size of, such
    as a kNN graph. This is a planning heuristic, not an allocation guarantee.
    It exists because an 8 GiB consumer card runs out on single-cell matrices
    long before a datacenter card does, and spilling does not rescue it: RAPIDS
    managed memory cannot oversubscribe GPU memory on Windows or WSL, so the
    failure there is an OOM rather than a slowdown.
    """
    available = vram_gb(free=True)
    if available is None:
        return None
    return matrix_gpu_gb(matrix) * overhead + reserve <= available


def backend_report(backend=None):
    """One-line human-readable summary for logs and provenance sidecars."""
    backend = backend if backend is not None else _backend_cache
    if backend is None or backend is _UNPROBED:
        return f"scanpy on CPU ({os.cpu_count()} logical cores)"
    name, total = "CUDA device", vram_gb()
    try:
        import cupy

        properties = cupy.cuda.runtime.getDeviceProperties(cupy.cuda.Device().id)
        name = properties["name"].decode() if isinstance(properties["name"], bytes) \
            else str(properties["name"])
    except Exception:
        pass
    version = getattr(backend, "__version__", "unknown")
    suffix = f", {total:.1f} GiB" if total else ""
    return f"rapids-singlecell {version} on {name}{suffix}"


def _to_host(value):
    """Bring one CuPy array or sparse matrix back to NumPy/SciPy; pass anything else through."""
    if type(value).__module__.startswith("cupy") and hasattr(value, "get"):
        return value.get()
    if isinstance(value, dict):
        return {key: _to_host(item) for key, item in value.items()}
    return value


def _results_to_host(*mappings):
    """Rewrite the named AnnData mappings in place so the object stays writable.

    GPU results land as CuPy arrays, which no h5ad writer accepts and no CPU
    plotting call understands. Every GPU helper calls this before returning, so
    callers downstream cannot tell which backend produced the result.
    """
    for mapping in mappings:
        for key in _named_keys(mapping):
            mapping[key] = _to_host(mapping[key])


def _to_device(matrix):
    """Copy a cells x genes or cells x components matrix to the device as float32."""
    import cupy
    from cupyx.scipy import sparse as cusparse
    from scipy import sparse

    if sparse.issparse(matrix):
        return cusparse.csr_matrix(matrix.tocsr().astype("float32"))
    return cupy.asarray(matrix, dtype="float32")


def _call_seeded(fn, *args, seed=0, **kwargs):
    """Call a rapids-singlecell function with a fixed seed under either spelling.

    rapids-singlecell renamed `random_state` to `rng`. The old name is still
    accepted on current releases and is the only one older releases know, so it
    is tried first and `rng` is the fallback.
    """
    try:
        return fn(*args, random_state=seed, **kwargs)
    except TypeError:
        return fn(*args, rng=seed, **kwargs)


def _gpu_failed(step, exc):
    """Report a GPU step that did not complete and hand the caller back to scanpy."""
    info(f"GPU {step} failed ({type(exc).__name__}: {exc}); using scanpy on the CPU")
    return False


def _gpu_pca(rsc, adata, n_comps):
    """Run PCA on the device. True only when it produced every key the toolkit reads.

    X is swapped for a device copy and the original object is put back
    afterwards, so the caller's expression matrix keeps its own dtype and
    sparsity and only the PCA outputs come from the GPU. Components are the same
    subspace as ARPACK's up to sign and solver tolerance, not bitwise equal.
    """
    if fits_in_vram(adata.X) is False:
        info(f"GPU: {matrix_gpu_gb(adata.X):.2f} GiB matrix does not fit in "
             f"{vram_gb(free=True):.2f} GiB free VRAM with workspace; PCA on the CPU")
        return False
    host_x = adata.X
    try:
        adata.X = _to_device(host_x)
        _call_seeded(rsc.pp.pca, adata, n_comps=n_comps)
    except Exception as exc:
        return _gpu_failed("PCA", exc)
    finally:
        adata.X = host_x
    _results_to_host(adata.obsm, adata.varm)
    adata.uns["pca"] = _to_host(adata.uns.get("pca", {}))
    if ("X_pca" not in adata.obsm or "PCs" not in adata.varm
            or "variance_ratio" not in adata.uns["pca"] or "params" not in adata.uns["pca"]):
        info("GPU PCA did not populate the expected keys; recomputing on the CPU")
        return False
    return True


def _gpu_neighbors(rsc, adata, *, use_rep, n_pcs, n_neighbors):
    """Build the kNN graph on the device. True when connectivities were written.

    cuVS computes an exact brute-force graph where scanpy's default pynndescent
    is approximate (upstream notes this as an expected difference), so the graph
    -- and therefore UMAP and cluster labels -- can differ between backends.
    """
    rep = adata.obsm[use_rep]
    # cuVS batches the distance computation, so the peak is the representation
    # plus the kNN result: float32 distances and int64 indices for every edge.
    graph_gb = adata.n_obs * n_neighbors * 12 / (1024**3)
    if fits_in_vram(rep, reserve=graph_gb) is False:
        info("GPU: representation plus kNN graph does not fit in free VRAM; "
             "neighbors on the CPU")
        return False
    try:
        adata.obsm[use_rep] = _to_device(rep)
        _call_seeded(rsc.pp.neighbors, adata, use_rep=use_rep, n_pcs=n_pcs,
                     n_neighbors=n_neighbors)
    except Exception as exc:
        return _gpu_failed("neighbors", exc)
    finally:
        adata.obsm[use_rep] = rep
    _results_to_host(adata.obsm, adata.obsp)
    if "connectivities" not in adata.obsp or "neighbors" not in adata.uns:
        info("GPU neighbors did not populate the expected keys; recomputing on the CPU")
        return False
    return True


def compute_pca(sc, adata, n_comps=50):
    """ARPACK requires fewer components than both cells and selected features."""
    n_features = int(adata.var["highly_variable"].sum()) if "highly_variable" in adata.var else adata.n_vars
    actual = min(n_comps, adata.n_obs - 1, n_features - 1)
    if actual < 1:
        die("PCA requires at least two cells and two selected genes")
    rsc = select_backend()
    if rsc is None or not _gpu_pca(rsc, adata, actual):
        sc.pp.pca(adata, n_comps=actual, svd_solver="arpack", random_state=0)
    return actual


def build_neighbors(sc, adata, *, use_rep="X_pca", n_pcs=40, n_neighbors=15):
    if use_rep not in adata.obsm:
        die(f"embedding '{use_rep}' is missing")
    if adata.n_obs < 3 or n_pcs < 1 or n_neighbors < 2:
        die("neighbors requires at least three cells, n_pcs >= 1 and n_neighbors >= 2")
    actual = min(n_pcs, adata.obsm[use_rep].shape[1])
    k = min(n_neighbors, adata.n_obs - 1)
    rsc = select_backend()
    if rsc is None or not _gpu_neighbors(rsc, adata, use_rep=use_rep, n_pcs=actual,
                                         n_neighbors=k):
        sc.pp.neighbors(adata, use_rep=use_rep, n_pcs=actual,
                        n_neighbors=k, random_state=0)


def compute_umap(sc, adata):
    """UMAP on the existing neighbor graph, on the device when one is usable.

    cuML's UMAP is a different implementation, not a port: coordinates are not
    reproducible against scanpy's umap-learn even with the seed fixed, and the
    two embeddings are comparable in layout rather than equal. Clusters and
    marker statistics do not depend on it; only the picture does.
    """
    if "neighbors" not in adata.uns:
        die("UMAP requires a neighborhood graph; run neighbors first")
    rsc = select_backend()
    if rsc is not None:
        try:
            _call_seeded(rsc.tl.umap, adata)
        except Exception as exc:
            _gpu_failed("UMAP", exc)
        else:
            _results_to_host(adata.obsm)
            if "X_umap" in adata.obsm:
                return
            info("GPU UMAP did not populate X_umap; recomputing on the CPU")
    sc.tl.umap(adata)


def cluster_cells(sc, adata, *, algorithm="leiden", resolution=0.5, key=None):
    """Leiden/Louvain on the existing graph, on the device when one is usable.

    The CPU path is scanpy's `flavor="igraph"` Leiden with `n_iterations=2`. The
    GPU path is cuGraph's, which takes neither argument and counts iterations
    differently, so cluster labels and even the cluster count can differ between
    backends at the same resolution. Fix the backend with SCI_DEVICE before
    comparing two runs, and re-derive annotations after switching.
    """
    if algorithm not in ("leiden", "louvain"):
        die(f"unknown clustering algorithm: {algorithm}")
    key = key or algorithm
    rsc = select_backend()
    if rsc is not None:
        gpu_fn = rsc.tl.leiden if algorithm == "leiden" else rsc.tl.louvain
        try:
            _call_seeded(gpu_fn, adata, resolution=resolution, key_added=key)
        except Exception as exc:
            _gpu_failed(algorithm, exc)
        else:
            if key in adata.obs:
                return key
            info(f"GPU {algorithm} did not populate obs['{key}']; recomputing on the CPU")
    kwargs = {"resolution": resolution, "key_added": key}
    if algorithm == "leiden":
        # flavor='igraph' is the scanpy 1.12 default-recommended Leiden backend.
        kwargs.update(flavor="igraph", n_iterations=2, directed=False)
    (sc.tl.leiden if algorithm == "leiden" else sc.tl.louvain)(adata, **kwargs)
    return key


def clear_graph(adata, clear_pca=False):
    """Remove standard derived results after changing expression or cell selection."""
    for key in ("neighbors", "umap", "paga", "rank_genes_groups", "diffmap_evals", "iroot"):
        adata.uns.pop(key, None)
    for key in list(adata.uns):
        if key.startswith("dendrogram_"):
            adata.uns.pop(key)
    for key in ("connectivities", "distances"):
        adata.obsp.pop(key, None)
    for key in ("X_umap", "X_tsne", "X_diffmap", "X_draw_graph_fa"):
        adata.obsm.pop(key, None)
    if clear_pca:
        adata.uns.pop("pca", None)
        adata.varm.pop("PCs", None)
        for key in ("X_pca", "X_pca_harmony"):
            adata.obsm.pop(key, None)


def integrate_harmony(adata, key):
    """Use harmonypy 2's cells x PCs output; Scanpy 1.12.4 transposes it incorrectly."""
    import numpy as np
    import harmonypy
    from importlib.metadata import version

    if int(version("harmonypy").split(".", 1)[0]) != 2:
        die("this Harmony path targets harmonypy 2.x; install harmonypy==2.0.2")
    pcs = np.asarray(adata.obsm["X_pca"], dtype=np.float64)
    result = harmonypy.run_harmony(pcs, adata.obs, key, random_state=0, ncores=1)
    corrected = np.asarray(result.Z_corr)
    if corrected.shape != pcs.shape or not np.isfinite(corrected).all():
        die("Harmony returned invalid cells x PCs coordinates")
    adata.obsm["X_pca_harmony"] = corrected
