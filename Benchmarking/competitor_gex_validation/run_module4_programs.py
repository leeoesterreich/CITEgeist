#!/usr/bin/env python
"""
Task 8: Standalone NMF program discovery for competitor GEX validation.

Runs sklearn NMF on each cell-type GEX layer, extracts top genes,
checks for MDK, and computes spatial Moran's I on program activities.
Does NOT import the CITEgeist model package (avoids cuOPT dependency).
"""

import argparse
import json
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from esda.moran import Moran
from libpysal import weights
from sklearn.decomposition import NMF

sys.path.insert(0, str(Path(__file__).resolve().parent))
from constants import CANDIDATE_GENES, MORANS_K_NEIGH, NMF_K_PROGRAMS, NMF_LAMBDA_SPARSITY, OUTPUT_ROOT

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
)
logger = logging.getLogger(__name__)


def run_nmf_for_celltype(X, gene_names, coords, k=5, alpha=0.01):
    """Run NMF on a cell-type expression matrix and return program results."""
    X_dense = X.toarray() if hasattr(X, "toarray") else np.array(X, dtype=np.float64)
    X_dense = np.clip(X_dense, 0, None)

    # Filter to spots with nonzero expression
    spot_mask = X_dense.sum(axis=1) > 0
    if spot_mask.sum() < 10:
        return None
    X_filt = X_dense[spot_mask]
    coords_filt = coords[spot_mask]

    model = NMF(n_components=k, alpha_W=alpha, init="nndsvda", max_iter=500, random_state=42)
    W = model.fit_transform(X_filt)  # (spots x k) — activities
    H = model.components_  # (k x genes) — gene loadings

    # Build spatial weights for Moran's I
    w = weights.KNN.from_array(coords_filt, k=min(MORANS_K_NEIGH, len(coords_filt) - 1))
    w.transform = "R"

    total_var = np.var(X_filt)
    programs = []
    for prog_idx in range(k):
        loadings = H[prog_idx]
        top_indices = np.argsort(loadings)[::-1][:50]
        top_50 = [gene_names[i] for i in top_indices]
        top_20 = top_50[:20]

        activity = W[:, prog_idx]
        var_explained = np.var(np.outer(activity, loadings)) / total_var if total_var > 0 else 0

        # Spatial Moran's I on program activity
        moran_I, moran_p = 0.0, 1.0
        if np.std(activity) > 1e-10:
            m = Moran(activity, w, permutations=999)
            moran_I, moran_p = m.I, m.p_sim

        programs.append(
            {
                "program_idx": prog_idx,
                "top_50_genes": top_50,
                "top_20_genes": top_20,
                "variance_explained": float(var_explained),
                "spatial_moran_I": float(moran_I),
                "spatial_moran_p": float(moran_p),
                "mdk_in_top20": "MDK" in top_20,
                "mdk_in_top50": "MDK" in top_50,
                "candidate_genes_in_top50": [g for g in top_50 if g in set(CANDIDATE_GENES)],
                "n_spots": int(spot_mask.sum()),
            }
        )

    return programs


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--method", required=True, choices=["cell2location", "tangram"])
    parser.add_argument("--sample", required=True)
    args = parser.parse_args()

    method, sample = args.method, args.sample
    sample_dir = OUTPUT_ROOT / method / sample

    gex_candidates = [
        sample_dir / f"{sample}_gex_layers.h5ad",
        sample_dir / "gex_layers.h5ad",
        sample_dir / "c2l_gex_layers.h5ad",
    ]
    prop_candidates = [
        sample_dir / f"{sample}_proportions.csv",
        sample_dir / "proportions.csv",
    ]
    layers_path = next((p for p in gex_candidates if p.exists()), None)
    prop_path = next((p for p in prop_candidates if p.exists()), None)
    out_dir = OUTPUT_ROOT / method / "programs"
    out_dir.mkdir(parents=True, exist_ok=True)

    if layers_path is None:
        logger.error("GEX layers not found in %s", sample_dir)
        sys.exit(1)
    if prop_path is None:
        logger.error("Proportions not found in %s", sample_dir)
        sys.exit(1)

    logger.info("=== Module 4 Programs: method=%s sample=%s ===", method, sample)

    adata = sc.read_h5ad(layers_path)
    adata.var_names_make_unique()
    prop_df = pd.read_csv(prop_path, index_col=0)
    coords = adata.obsm["spatial"]
    gene_names = list(adata.var_names)

    logger.info("Loaded: %d spots x %d genes, %d layers", adata.n_obs, adata.n_vars, len(adata.layers))

    # Align spots
    common = adata.obs_names.intersection(prop_df.index)
    if len(common) < adata.n_obs:
        adata = adata[common].copy()
        prop_df = prop_df.loc[common]
        coords = adata.obsm["spatial"]

    # Strip C2L prefixes from proportion columns for matching
    C2L_PROP_PREFIX = "q05cell_abundance_w_sf_means_per_cluster_mu_fg_"
    C2L_LAYER_PREFIX = "means_per_cluster_mu_fg_"

    def strip_prefix(name):
        for pfx in [C2L_PROP_PREFIX, C2L_LAYER_PREFIX]:
            if name.startswith(pfx):
                return name[len(pfx) :]
        return name

    prop_clean = {strip_prefix(c): c for c in prop_df.columns}

    # Run NMF per cell-type layer
    summary = {}
    for layer_name in adata.layers:
        ct = strip_prefix(layer_name.replace("_genes_pass1", ""))
        # Only run for types with sufficient proportion
        prop_col = prop_clean.get(ct)
        mean_prop = prop_df[prop_col].mean() if prop_col else 0
        if mean_prop < 0.1:
            logger.info("  Skipping %s (mean proportion %.3f < 0.1)", ct, mean_prop)
            continue

        X = adata.layers[layer_name]
        logger.info("  Running NMF for %s (mean prop=%.3f) ...", ct, mean_prop)
        programs = run_nmf_for_celltype(X, gene_names, coords, k=NMF_K_PROGRAMS, alpha=NMF_LAMBDA_SPARSITY)

        if programs is None:
            logger.warning("  %s: too few nonzero spots, skipped", ct)
            continue

        summary[ct] = programs
        for p in programs:
            if p["mdk_in_top20"] or p["mdk_in_top50"]:
                logger.info(
                    "    Prog %d: MDK in top20=%s top50=%s, Moran I=%.3f (p=%.4f)",
                    p["program_idx"],
                    p["mdk_in_top20"],
                    p["mdk_in_top50"],
                    p["spatial_moran_I"],
                    p["spatial_moran_p"],
                )

    out_path = out_dir / f"{sample}_programs.json"
    with open(out_path, "w") as f:
        json.dump({"method": method, "sample": sample, "programs": summary}, f, indent=2)
    logger.info("Saved to %s (%d cell types)", out_path, len(summary))


if __name__ == "__main__":
    main()
