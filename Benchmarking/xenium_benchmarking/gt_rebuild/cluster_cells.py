"""Level-1 Leiden clustering of Xenium RNA (all cells), the first step of the GT rebuild."""

import argparse
import json
import re
from pathlib import Path
from typing import Tuple

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.sparse import csc_matrix

from . import gt_config, manifest

CONTROL_PATTERN = re.compile(r"^(NegControl|BLANK|antisense|UnassignedCodeword|Intergenic)", re.IGNORECASE)


def load_xenium_rna(h5_path: Path, cells_parquet: Path) -> Tuple[ad.AnnData, dict]:
    """Load raw RNA counts, drop control probes and zero-count cells; verify barcodes match cells.parquet."""
    recipe = gt_config.RECIPE
    with h5py.File(h5_path, "r") as f:
        m = f["matrix"]
        names = m["features"]["name"][:].astype(str)
        ftypes = m["features"]["feature_type"][:].astype(str)
        barcodes = m["barcodes"][:].astype(str)
        X = csc_matrix((m["data"][:], m["indices"][:], m["indptr"][:]), shape=tuple(m["shape"][:]))

    cell_ids = pd.read_parquet(cells_parquet, columns=["cell_id"])["cell_id"]
    cell_ids = np.array([c.decode() if isinstance(c, bytes) else str(c) for c in cell_ids])
    if len(cell_ids) != len(barcodes) or not np.array_equal(cell_ids, barcodes):
        raise ValueError("h5 barcodes differ from cells.parquet cell_id in content or order")

    keep = np.where(ftypes == recipe["feature_type"])[0]
    X = X[keep, :].T.tocsr().astype(np.float32)
    genes = names[keep]
    is_control = np.array([bool(CONTROL_PATTERN.match(g)) for g in genes])
    n_control = int(is_control.sum()) if recipe["remove_control_probes"] else 0
    if recipe["remove_control_probes"]:
        X = X[:, np.where(~is_control)[0]]
        genes = genes[~is_control]

    totals = np.asarray(X.sum(axis=1)).ravel()
    zero = totals <= 0
    dropped = barcodes[zero].tolist() if recipe["drop_zero_count_cells"] else []
    if recipe["drop_zero_count_cells"]:
        X = X[np.where(~zero)[0], :]
        barcodes = barcodes[~zero]

    adata = ad.AnnData(X)
    adata.obs_names = barcodes
    adata.var_names = genes
    qc = {
        "n_cells_total": int(len(cell_ids)),
        "n_zero_count_dropped": len(dropped),
        "n_genes": int(adata.n_vars),
        "n_control_removed": n_control,
        "dropped_ids": dropped,
    }
    return adata, qc


def cluster_level1(
    adata: ad.AnnData,
    resolution: float,
    seed: int,
    n_pcs: int = gt_config.RECIPE["level1"]["n_pcs"],
    n_neighbors: int = gt_config.RECIPE["level1"]["n_neighbors"],
) -> pd.Series:
    """Normalise, PCA, kNN and Leiden on a copy; return string cluster labels indexed by cell id."""
    recipe = gt_config.RECIPE
    lv = recipe["level1"]
    work = ad.AnnData(adata.X.copy())
    work.obs_names = adata.obs_names.copy()
    work.var_names = adata.var_names.copy()
    sc.pp.normalize_total(work, target_sum=recipe["normalize_total_target_sum"])
    if recipe["log1p"]:
        sc.pp.log1p(work)
    sc.pp.pca(work, n_comps=n_pcs, random_state=seed, svd_solver=recipe["pca_solver"])
    sc.pp.neighbors(work, n_neighbors=n_neighbors, n_pcs=n_pcs, random_state=seed, metric=recipe["neighbors_metric"])
    sc.tl.leiden(
        work,
        resolution=resolution,
        key_added="L1",
        flavor=lv["leiden_flavor"],
        n_iterations=lv["leiden_n_iterations"],
        directed=lv["leiden_directed"],
        random_state=seed,
    )
    return work.obs["L1"].astype(str).rename("cluster")


def main(argv=None) -> None:
    lv = gt_config.RECIPE["level1"]
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--resolution", type=float, default=lv["leiden_resolution"])
    p.add_argument("--seed", type=int, default=gt_config.RECIPE["seed"])
    p.add_argument("--out-dir", type=Path, required=True)
    p.add_argument("--xenium-dir", type=Path, default=Path(gt_config.XENIUM_DIR))
    p.add_argument("--manifest", type=Path, default=None, help="default: <out-dir>/../manifest.json if it exists")
    p.add_argument("--replace-section", action="store_true", help="--resume rerun: overwrite the 'cluster' section")
    args = p.parse_args(argv)

    adata, qc = load_xenium_rna(args.xenium_dir / "cell_feature_matrix.h5", args.xenium_dir / "cells.parquet")
    labels = cluster_level1(adata, args.resolution, args.seed)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    labels.rename_axis("cell_id").to_csv(args.out_dir / "level1_clusters.csv")
    (args.out_dir / "level1_qc.json").write_text(json.dumps(qc, indent=2))

    mpath = args.manifest or args.out_dir.parent / "manifest.json"
    if mpath.exists():
        existing = json.loads(mpath.read_text())
        if "versions" not in existing:
            manifest.append_section(mpath, "versions", manifest.collect_versions())
        manifest.put_section(
            mpath,
            "cluster",
            {
                "resolution": args.resolution,
                "seed": args.seed,
                "n_clusters": int(labels.nunique()),
                "n_cells": int(len(labels)),
                "clusters_sha256": manifest.sha256_file(args.out_dir / "level1_clusters.csv"),
            },
            args.replace_section,
        )
    else:
        print(f"manifest not found at {mpath}; wrote outputs only")


if __name__ == "__main__":
    main()
