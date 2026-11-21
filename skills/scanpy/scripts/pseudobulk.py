#!/usr/bin/env python3
"""Sum raw counts by biological sample and cell type for downstream DE.

Exports genes x pseudobulk samples and aligned metadata. Transpose counts for
PyDESeq2's samples x genes API. Analyze one cell type at a time with independent
biological replicates; cells and sample x cell-type rows are not new donors.

Example:
    python pseudobulk.py annotated.h5ad --by sample cell_type --metadata condition donor --out-prefix results/pb
"""

import argparse
import os

from _common import configure_scanpy, die, info, load_anndata, validate_counts


def aggregate_counts(adata, by, layer="counts", metadata=(), min_cells=1):
    """Collision-free groups, exact layer selection, constant group metadata."""
    import numpy as np
    import pandas as pd
    import scanpy as sc

    if min_cells < 1 or not by or len(set(by)) != len(by):
        die("--by must contain unique grouping columns and --min-cells must be positive")
    for key in list(by) + list(metadata):
        if key not in adata.obs or adata.obs[key].isna().any():
            die(f"obs column '{key}' must exist and contain no missing values")
    if layer is not None and layer not in adata.layers:
        die(f"raw-count layer '{layer}' is missing; no fallback to X is permitted")
    matrix = adata.X if layer is None else adata.layers[layer]
    validate_counts(matrix, "pseudobulk counts", require_nonzero_cells=False)
    if not adata.var_names.is_unique:
        die("gene identifiers must be unique before aggregation")
    # Scanpy joins multiple keys with underscores; ambiguous labels can collide.
    # Aggregate by opaque IDs and retain the original grouping columns separately.
    codes, _ = pd.factorize(pd.MultiIndex.from_frame(adata.obs[list(by)]), sort=False)
    ids = pd.Index([f"pb{i:06d}" for i in range(int(codes.max()) + 1)], name="sample_id")
    cols = list(dict.fromkeys(list(by) + list(metadata)))
    # One grouped pass. The previous per-group `obs.iloc[flatnonzero(...)]` loop
    # was O(n_groups x n_cells): 100 donors x 20 cell types over 500k cells is
    # ~1e9 element scans plus 2000 DataFrame copies, minutes where the
    # aggregation itself is seconds.
    grouped = adata.obs.groupby(pd.Categorical.from_codes(codes, ids), observed=True)
    if metadata:
        varying = grouped[list(metadata)].nunique(dropna=False).max(axis=1) != 1
        if varying.any():
            die("metadata is not constant within pseudobulk group(s): "
                + ", ".join(map(str, varying.index[varying][:5])))
    meta = grouped[cols].first()
    meta["n_cells"] = grouped.size()
    meta.index = pd.Index(meta.index.astype(str), name="sample_id")
    meta = meta.reindex(ids)
    # No copy of expression and no mutation of the source object.
    import anndata as ad
    work = ad.AnnData(matrix, obs=pd.DataFrame({"group": pd.Categorical(ids[codes])},
                                              index=adata.obs_names), var=adata.var.copy())
    pb = sc.get.aggregate(work, by="group", func="sum")
    mat = pb.layers["sum"]  # aggregate returns layers, not an expression matrix in X
    arr = mat.toarray() if hasattr(mat, "toarray") else np.asarray(mat)
    ordered = pb.obs["group"].astype(str).to_numpy()
    meta = meta.loc[ordered]
    keep = meta["n_cells"].to_numpy() >= min_cells
    if not keep.any():
        die("no pseudobulk groups meet --min-cells")
    if not keep.all():
        info(f"dropped {int((~keep).sum())} of {len(keep)} profiles below "
             f"--min-cells {min_cells}")
    counts = pd.DataFrame(arr[keep].T, index=adata.var_names,
                          columns=meta.index[keep])
    counts.index.name = "gene_id"
    return counts, meta.iloc[np.flatnonzero(keep)]


def report_replication(meta, by, metadata, allow_underpowered=False):
    """Count replicates per arm per stratum and refuse designs that cannot support DE.

    Convention: --by <replicate column> [<stratum columns>...]. Each pseudobulk
    profile is one replicate within one stratum, so counting profiles per level
    of a design factor counts donors per arm.
    """
    import pandas as pd

    strata = [c for c in list(by)[1:] if c in meta]
    worst = None
    for factor in metadata:
        levels = meta[factor].astype(str)
        n_levels = levels.nunique()
        # Skip identifiers (one level per profile) and constants.
        if n_levels < 2 or n_levels > max(2, len(meta) // 2):
            continue
        keys = strata + [factor]
        sizes = meta.assign(**{factor: levels}).groupby(keys, observed=True).size()
        for label, n in sizes.items():
            if worst is None or n < worst[1]:
                worst = (f"{factor}={label}" if not strata else str(label), int(n))
        info(f"replicates per arm for '{factor}':\n"
             + pd.DataFrame({"n_replicates": sizes}).to_string())
    if worst is None:
        info("no design factor found in --metadata; replication not assessed")
        return
    label, n = worst
    if n < 2:
        msg = (f"smallest arm has {n} replicate ({label}); DE across replicates is "
               "not defined below 2 per arm")
        if not allow_underpowered:
            die(msg + ". Re-run with --allow-underpowered to emit the table anyway.")
        info(f"WARNING: {msg}")
    elif n < 4:
        info(f"WARNING: smallest arm has {n} replicates ({label}). At 4 per arm and "
             "donor SD 0.35, recovery of true log2FC 1.0 was ~34% and log2FC 0.5 "
             "~0%; a gene not called is uninformative, not unchanged. Report the "
             "replicate count with every conclusion.")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("input", help="Input .h5ad with full-gene raw counts")
    p.add_argument("--by", nargs="+", required=True, help="Sample and cell-type obs columns")
    source = p.add_mutually_exclusive_group()
    source.add_argument("--layer", default="counts", help="Raw counts layer (default counts)")
    source.add_argument("--use-x", action="store_true", help="Explicitly declare X as raw counts")
    p.add_argument("--func", choices=["sum"], default="sum", help="Count-based DE requires sums")
    p.add_argument("--metadata", nargs="+", default=[], help="Group-constant donor/condition/design columns")
    p.add_argument("--min-cells", type=int, default=10,
                   help="Minimum contributing cells per profile (default 10)")
    p.add_argument("--allow-underpowered", action="store_true",
                   help="Emit the table even when an arm has fewer than 2 replicates")
    p.add_argument("--out-prefix", default="pseudobulk")
    args = p.parse_args()
    configure_scanpy()
    adata = load_anndata(args.input)
    if "highly_variable" in adata.var and bool(adata.var["highly_variable"].all()):
        info("WARNING: every gene is flagged highly_variable; this input looks "
             "HVG-subset. Pseudobulk from the full-gene checkpoint instead.")
    counts, meta = aggregate_counts(adata, args.by, None if args.use_x else args.layer,
                                    args.metadata, args.min_cells)
    report_replication(meta, args.by, args.metadata, args.allow_underpowered)
    os.makedirs(os.path.dirname(os.path.abspath(args.out_prefix)), exist_ok=True)
    counts.to_csv(f"{args.out_prefix}_counts.csv")
    meta.to_csv(f"{args.out_prefix}_samples.csv")
    info(f"Wrote {counts.shape[0]} genes x {counts.shape[1]} pseudobulk profiles")
    info("Transpose counts for PyDESeq2; preserve IDs as strings, align metadata, and fit each cell type with biological replicates.")


if __name__ == "__main__":
    main()
