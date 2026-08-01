#!/usr/bin/env python
"""
Persist bivariate Moran's I (MDK vs secretory score) for all 12 samples.

Computes bivariate Moran's I using esda.Moran_BV with queen contiguity
from Visium spatial coordinates. Saves per-sample results as CSV.

Output: mdk_analysis/if_analysis/bivariate_morans_results.csv
"""

import os
import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq
from scipy.spatial import cKDTree
import warnings

warnings.filterwarnings("ignore")

PROJECT_DIR = "/path/to/CITEgeist_analysis"
DATA_ROOT = "/path/to/CITEgeist_public_data/processed_files"
OUTPUT_DIR = os.path.join(PROJECT_DIR, "mdk_analysis/outputs/if_analysis")
MORPH_ASSIGN_DIR = os.path.join(PROJECT_DIR, "output/morphology_assignment_v3")

SAMPLES = [
    "HCC22-088-P1-S1",
    "HCC22-088-P1-S2",
    "HCC22-088-P2-S1",
    "HCC22-088-P2-S2",
    "HCC22-088-P3-S1_A",
    "HCC22-088-P3-S2",
    "HCC22-088-P4-S1",
    "HCC22-088-P4-S2_1i_rep",
    "HCC22-088-P5-S1",
    "HCC22-088-P5-S2_F_rep",
    "HCC22-088-P6-S1",
    "HCC22-088-P6-S2_D",
]

SECRETORY_GENES = ["HSP90B1", "HSPA5", "CALR", "CANX", "PDIA4", "PDIA6", "SEC23A", "SEC61B", "ATF6", "MAN1A1", "XBP1"]

N_PERM = 999


def compute_bivariate_moran(x, y, coords, n_perm=N_PERM, seed=42):
    """Compute bivariate Moran's I with permutation test."""
    try:
        from esda.moran import Moran_BV
        from libpysal.weights import KNN
    except ImportError:
        # Fallback: manual computation
        return _manual_bivariate_moran(x, y, coords, n_perm, seed)

    w = KNN.from_array(coords, k=6)
    w.transform = "r"
    np.random.seed(seed)
    bv = Moran_BV(x, y, w, permutations=n_perm)
    return float(bv.I), float(bv.p_sim)


def _manual_bivariate_moran(x, y, coords, n_perm=N_PERM, seed=42):
    """Manual bivariate Moran's I if esda not available."""
    tree = cKDTree(coords)
    # k=6 nearest neighbors
    dists, indices = tree.query(coords, k=7)  # includes self
    indices = indices[:, 1:]  # remove self

    n = len(x)
    x_std = (x - x.mean()) / x.std()
    y_std = (y - y.mean()) / y.std()

    # Bivariate Moran's I = sum_i(x_i * sum_j(w_ij * y_j)) / n
    lag_y = np.array([y_std[indices[i]].mean() for i in range(n)])
    I_obs = np.sum(x_std * lag_y) / n

    # Permutation test
    rng = np.random.RandomState(seed)
    count = 0
    for _ in range(n_perm):
        perm_y = rng.permutation(y_std)
        lag_perm = np.array([perm_y[indices[i]].mean() for i in range(n)])
        I_perm = np.sum(x_std * lag_perm) / n
        if abs(I_perm) >= abs(I_obs):
            count += 1

    return float(I_obs), float(count / n_perm)


def main():
    results = []

    for sample in SAMPLES:
        print(f"\n=== {sample} ===")
        sc_path = os.path.join(MORPH_ASSIGN_DIR, sample, f"{sample}_single_cell.h5ad")
        if not os.path.isfile(sc_path):
            print(f"  SKIP: single_cell.h5ad not found at {sc_path}")
            results.append(
                {
                    "sample": sample,
                    "morans_I": np.nan,
                    "p_value": np.nan,
                    "n_cells": 0,
                    "mdk_nonzero": 0,
                    "significant": False,
                    "note": "single_cell.h5ad not found",
                }
            )
            continue

        import scanpy as sc_

        adata = sc_.read_h5ad(sc_path)
        # Cell centroid coordinates (x, y) stored in obsm["spatial"]
        coords = adata.obsm["spatial"].copy()
        # Drop cells with NaN/inf coords
        valid_mask = np.isfinite(coords).all(axis=1)
        if not valid_mask.all():
            print(f"  Dropping {(~valid_mask).sum()} cells with NaN/inf coords")
            adata = adata[valid_mask]
            coords = coords[valid_mask]
        print(f"  Loaded: {adata.shape[0]} cells, {adata.shape[1]} genes")

        # Get MDK expression
        gene_mask = adata.var_names.isin(["MDK"])
        if gene_mask.sum() == 0:
            print("  SKIP: MDK not in gene list")
            results.append(
                {
                    "sample": sample,
                    "morans_I": np.nan,
                    "p_value": np.nan,
                    "n_cells": adata.shape[0],
                    "mdk_nonzero": 0,
                    "significant": False,
                    "note": "MDK not in single_cell var",
                }
            )
            continue

        mdk_x = adata[:, "MDK"].X
        if hasattr(mdk_x, "toarray"):
            mdk = mdk_x.toarray().ravel().astype(float)
        else:
            mdk = np.asarray(mdk_x).ravel().astype(float)
        mdk_nonzero = int(np.sum(mdk > 0))

        if mdk_nonzero < 10:
            print(f"  SKIP: only {mdk_nonzero} nonzero MDK spots")
            results.append(
                {
                    "sample": sample,
                    "morans_I": 0.0,
                    "p_value": 1.0,
                    "n_cells": adata.shape[0],
                    "mdk_nonzero": mdk_nonzero,
                    "significant": False,
                    "note": f"too few nonzero ({mdk_nonzero})",
                }
            )
            continue

        # Compute secretory score
        present = [g for g in SECRETORY_GENES if g in adata.var_names]
        print(f"  Secretory genes present: {len(present)}/{len(SECRETORY_GENES)}")

        sec_x = adata[:, present].X
        if hasattr(sec_x, "toarray"):
            sec_mat = sec_x.toarray().astype(float)
        else:
            sec_mat = np.asarray(sec_x).astype(float)
        # Z-score each gene, then mean
        from scipy.stats import zscore

        sec_z = zscore(sec_mat, axis=0, nan_policy="omit")
        sec_z = np.nan_to_num(sec_z, 0)
        sec_score = sec_z.mean(axis=1)

        try:
            I_val, p_val = compute_bivariate_moran(mdk, sec_score, coords)
        except Exception as e:
            print(f"  ERROR computing Moran's I: {e}")
            results.append(
                {
                    "sample": sample,
                    "morans_I": np.nan,
                    "p_value": np.nan,
                    "n_cells": adata.shape[0],
                    "mdk_nonzero": mdk_nonzero,
                    "significant": False,
                    "note": f"computation error: {e}",
                }
            )
            continue
        sig = p_val <= 0.001  # min p with 999 perms is 1/999≈0.001; use <= not <

        print(f"  Moran's I = {I_val:.4f}, p = {p_val:.4f}, significant = {sig}")

        results.append(
            {
                "sample": sample,
                "morans_I": round(I_val, 4),
                "p_value": round(p_val, 4),
                "n_cells": adata.shape[0],
                "mdk_nonzero": mdk_nonzero,
                "significant": sig,
                "n_secretory_genes": len(present),
                "n_permutations": N_PERM,
                "note": "",
            }
        )

    df = pd.DataFrame(results)
    out_path = os.path.join(OUTPUT_DIR, "bivariate_morans_results.csv")
    df.to_csv(out_path, index=False)
    print(f"\n\nSaved to {out_path}")
    print(df.to_string(index=False))

    # Summary stats
    sig_mask = df["significant"] == True
    print(f"\n{sig_mask.sum()}/{len(df)} samples significant (p < 0.001)")
    print(f"Mean I (significant): {df.loc[sig_mask, 'morans_I'].mean():.3f}")


if __name__ == "__main__":
    main()
