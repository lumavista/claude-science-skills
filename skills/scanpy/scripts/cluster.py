#!/usr/bin/env python3
"""
Leiden clustering on a precomputed neighborhood graph.

Runs Leiden at one or more resolutions and writes a UMAP colored by each
clustering. Requires that ``sc.pp.neighbors`` has already been run
(use reduce_dimensions.py first).

Clustering runs on rapids-singlecell (cuGraph) when a CUDA device is usable and
on scanpy otherwise; set ``SCI_DEVICE=cpu`` to force scanpy. Labels and cluster
counts are not identical between the two backends at the same resolution.

Examples:
    python cluster.py reduced.h5ad -o clustered.h5ad --resolution 0.5
    python cluster.py reduced.h5ad -o clustered.h5ad --resolution 0.3 0.5 0.8 1.0
    python cluster.py reduced.h5ad -o clustered.h5ad --algorithm louvain
"""

import argparse

from _common import add_io_args, cluster_cells, configure_scanpy, die, info, load_anndata, save_anndata


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    add_io_args(p, default_output="clustered.h5ad")
    p.add_argument("--resolution", type=float, nargs="+", default=[0.5],
                   help="One or more resolutions (default 0.5). Higher = more clusters")
    p.add_argument("--algorithm", default="leiden", choices=["leiden", "louvain"],
                   help="Clustering algorithm (default leiden)")
    p.add_argument("--key", default=None,
                   help="obs key for the result (single resolution only; "
                        "default '<algorithm>'). Multiple resolutions use '<algorithm>_<res>'")
    p.add_argument("--no-plots", action="store_true", help="Skip UMAP plots")
    args = p.parse_args()

    sc = configure_scanpy(figdir=args.figdir)
    adata = load_anndata(args.input)
    if "neighbors" not in adata.uns:
        die("no neighborhood graph found. Run reduce_dimensions.py first.")

    keys = []
    for res in args.resolution:
        if len(args.resolution) == 1:
            key = args.key or args.algorithm
        else:
            key = f"{args.algorithm}_{res}"
        cluster_cells(sc, adata, algorithm=args.algorithm, resolution=res, key=key)
        n = adata.obs[key].nunique()
        info(f"{key}: {n} clusters at resolution {res}")
        keys.append(key)

    if not args.no_plots:
        sc.pl.umap(adata, color=keys, legend_loc="on data",
                   show=False, save="_clusters.png")

    save_anndata(adata, args.output)


if __name__ == "__main__":
    main()
