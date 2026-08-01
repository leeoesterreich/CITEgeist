#!/usr/bin/env python
"""
spatial_morans.py — Unified bivariate Moran's I: MDK vs secretory genes.

NOTE ON `cancer_cell_morans.csv`: despite its name, the composite written by
run_cancer_cell_level() is computed over ALL cell types, not cancer cells only
(mean Moran's I 0.290). The cancer-only value reported in the manuscript
(mean 0.254) is produced by persist_bivariate_morans.py and written to
if_analysis/bivariate_morans_results.csv. Use that file for the cancer-only
figure; this script's trajectory_deltas.csv is what feeds the biopsy->surgery
panel.

Three resolutions:
  1. Spot-level per-gene (MDK vs each of 11 secretory genes, + 11 random controls)
  2. Cancer-cell-level composite (MDK vs mean-z secretory score)
  3. Trajectory deltas (biopsy→surgery changes in Moran's I)

Reads sample_manifest.csv for inclusion/exclusion, sample_pairs.csv for pairing.

Usage:
    python -u spatial_morans.py --output-dir mdk_analysis/outputs/spatial
"""

import argparse
import json
import logging
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import zscore
from statsmodels.stats.multitest import multipletests

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

# Root of your local copy of the CITEgeist analysis outputs (see README).
PROJECT_DIR = Path("/path/to/CITEgeist_analysis")
DATA_ROOT = Path("/path/to/CITEgeist_public_data/processed_files")
MORPH_ASSIGN_DIR = PROJECT_DIR / "output" / "morphology_assignment_v3"

SECRETORY_GENES = [
    "HSP90B1",
    "HSPA5",
    "CALR",
    "CANX",
    "PDIA4",
    "PDIA6",
    "SEC23A",
    "SEC61B",
    "ATF6",
    "MAN1A1",
    "XBP1",
]

SHORT_TO_FULL = {
    "P1-S1": "HCC22-088-P1-S1",
    "P1-S2": "HCC22-088-P1-S2",
    "P2-S1": "HCC22-088-P2-S1",
    "P2-S2": "HCC22-088-P2-S2",
    "P3-S1": "HCC22-088-P3-S1_A",
    "P3-S2": "HCC22-088-P3-S2",
    "P4-S1": "HCC22-088-P4-S1",
    "P4-S2": "HCC22-088-P4-S2_1i_rep",
    "P5-S1": "HCC22-088-P5-S1",
    "P5-S2": "HCC22-088-P5-S2_F_rep",
    "P6-S1": "HCC22-088-P6-S1",
    "P6-S2": "HCC22-088-P6-S2_D",
}

N_PERM = 999


def load_manifest(manifest_path: str) -> list:
    """Load sample manifest and return list of included full sample IDs."""
    df = pd.read_csv(manifest_path)
    short_ids = df.loc[df["included"] == True, "sample_id"].tolist()
    return [SHORT_TO_FULL[s] for s in short_ids if s in SHORT_TO_FULL]


def load_pairs(pairs_path: str) -> pd.DataFrame:
    """Load sample pairing CSV."""
    return pd.read_csv(pairs_path)


def compute_bivariate_moran(x, y, coords, n_perm=N_PERM, seed=42):
    """Compute bivariate Moran's I with permutation test.

    Uses esda.Moran_BV if available, otherwise falls back to manual KNN computation.
    """
    try:
        from esda.moran import Moran_BV
        from libpysal.weights import KNN

        k = min(6, len(x) - 1)
        w = KNN.from_array(coords, k=k)
        w.transform = "r"
        np.random.seed(seed)
        bv = Moran_BV(x, y, w, permutations=n_perm)
        return float(bv.I), float(bv.p_sim)
    except ImportError:
        return _manual_bivariate_moran(x, y, coords, n_perm, seed)


def _manual_bivariate_moran(x, y, coords, n_perm=N_PERM, seed=42):
    """Manual bivariate Moran's I fallback when esda is unavailable."""
    tree = cKDTree(coords)
    k = min(6, len(x) - 1)
    _, indices = tree.query(coords, k=k + 1)
    indices = indices[:, 1:]

    n = len(x)
    x_std = (x - x.mean()) / x.std()
    y_std = (y - y.mean()) / y.std()

    lag_y = np.array([y_std[indices[i]].mean() for i in range(n)])
    I_obs = np.sum(x_std * lag_y) / n

    rng = np.random.RandomState(seed)
    count = 0
    for _ in range(n_perm):
        perm_y = rng.permutation(y_std)
        lag_perm = np.array([perm_y[indices[i]].mean() for i in range(n)])
        I_perm = np.sum(x_std * lag_perm) / n
        if abs(I_perm) >= abs(I_obs):
            count += 1

    return float(I_obs), float(count / n_perm)


def _extract_gene_vector(adata, gene):
    """Extract a gene expression vector from AnnData, handling sparse matrices."""
    x = adata[:, gene].X
    if hasattr(x, "toarray"):
        return x.toarray().ravel().astype(float)
    return np.asarray(x).ravel().astype(float)


def run_spot_level(samples, output_dir):
    """Run per-gene bivariate Moran's I at spot resolution (raw Visium).

    For each sample: MDK vs each secretory gene + 11 random controls.
    FDR correction across all secretory tests.
    """
    import squidpy as sq

    rows = []
    for sample in samples:
        logger.info(f"[spot] {sample}")
        visium_dir = str(DATA_ROOT / sample / "outs")
        if not os.path.isdir(visium_dir):
            logger.warning(f"  SKIP: {visium_dir} not found")
            continue

        adata = sq.read.visium(visium_dir, counts_file="filtered_feature_bc_matrix.h5", load_images=True)
        coords = adata.obsm["spatial"].copy()

        valid = np.isfinite(coords).all(axis=1)
        if not valid.all():
            adata = adata[valid]
            coords = coords[valid]

        if "feature_types" in adata.var.columns:
            adata = adata[:, adata.var["feature_types"] == "Gene Expression"]

        dup = adata.var_names.duplicated(keep="first")
        if dup.any():
            adata = adata[:, ~dup]

        if "MDK" not in adata.var_names:
            logger.warning(f"  SKIP: no MDK in {sample}")
            continue

        mdk = _extract_gene_vector(adata, "MDK")
        mdk_nz = int(np.sum(mdk > 0))
        if mdk_nz < 10:
            logger.warning(f"  SKIP: only {mdk_nz} nonzero MDK in {sample}")
            continue

        logger.info(f"  {adata.shape[0]} spots, MDK nonzero: {mdk_nz}")

        for gene in SECRETORY_GENES:
            if gene not in adata.var_names:
                rows.append(
                    {
                        "sample": sample,
                        "gene": gene,
                        "gene_type": "secretory",
                        "morans_I": np.nan,
                        "p_value": np.nan,
                        "mdk_nonzero": mdk_nz,
                        "gene_nonzero": 0,
                        "n_spots": adata.shape[0],
                        "note": "gene not found",
                    }
                )
                continue

            gene_expr = _extract_gene_vector(adata, gene)
            gene_nz = int(np.sum(gene_expr > 0))

            try:
                I_val, p_val = compute_bivariate_moran(mdk, gene_expr, coords)
            except Exception:
                I_val, p_val = np.nan, np.nan

            rows.append(
                {
                    "sample": sample,
                    "gene": gene,
                    "gene_type": "secretory",
                    "morans_I": round(I_val, 4) if not np.isnan(I_val) else np.nan,
                    "p_value": round(p_val, 4) if not np.isnan(p_val) else np.nan,
                    "mdk_nonzero": mdk_nz,
                    "gene_nonzero": gene_nz,
                    "n_spots": adata.shape[0],
                    "note": "",
                }
            )
            logger.info(f"  {gene}: I={I_val:.4f} p={p_val:.4f}")

        rng = np.random.RandomState(123)
        non_sec = [g for g in adata.var_names if g not in SECRETORY_GENES and g != "MDK"]
        ctrl_genes = rng.choice(non_sec, size=min(11, len(non_sec)), replace=False)
        for gene in ctrl_genes:
            gene_expr = _extract_gene_vector(adata, gene)
            gene_nz = int(np.sum(gene_expr > 0))
            try:
                I_val, p_val = compute_bivariate_moran(mdk, gene_expr, coords)
            except Exception:
                I_val, p_val = np.nan, np.nan
            rows.append(
                {
                    "sample": sample,
                    "gene": gene,
                    "gene_type": "control",
                    "morans_I": round(I_val, 4) if not np.isnan(I_val) else np.nan,
                    "p_value": round(p_val, 4) if not np.isnan(p_val) else np.nan,
                    "mdk_nonzero": mdk_nz,
                    "gene_nonzero": gene_nz,
                    "n_spots": adata.shape[0],
                    "note": "",
                }
            )

    df = pd.DataFrame(rows)
    if len(df) == 0:
        logger.warning("No spot-level results produced")
        return df

    sec_mask = (df["gene_type"] == "secretory") & df["p_value"].notna()
    if sec_mask.sum() > 0:
        _, fdr_q, _, _ = multipletests(df.loc[sec_mask, "p_value"].values, method="fdr_bh")
        df.loc[sec_mask, "fdr_q"] = np.round(fdr_q, 4)

    ctrl_mask = (df["gene_type"] == "control") & df["p_value"].notna()
    if ctrl_mask.sum() > 0:
        _, fdr_q_ctrl, _, _ = multipletests(df.loc[ctrl_mask, "p_value"].values, method="fdr_bh")
        df.loc[ctrl_mask, "fdr_q"] = np.round(fdr_q_ctrl, 4)

    out_path = os.path.join(output_dir, "spot_level_pergene_morans.csv")
    df.to_csv(out_path, index=False)
    logger.info(f"Spot-level per-gene results saved to {out_path}")

    sec_df = df[df["gene_type"] == "secretory"]
    summary_rows = []
    for sample in sec_df["sample"].unique():
        sdf = sec_df[sec_df["sample"] == sample]
        mean_I = sdf["morans_I"].mean()
        n_sig = int((sdf["fdr_q"] < 0.05).sum()) if "fdr_q" in sdf.columns else 0
        summary_rows.append(
            {
                "sample": sample,
                "morans_I": round(mean_I, 4),
                "n_sig_genes": n_sig,
                "n_tested": len(sdf),
                "significant": n_sig > len(sdf) * 0.5,
            }
        )
    summary_df = pd.DataFrame(summary_rows)
    summary_path = os.path.join(output_dir, "spot_level_morans.csv")
    summary_df.to_csv(summary_path, index=False)
    logger.info(f"Spot-level summary saved to {summary_path}")

    return df


def run_cancer_cell_level(samples, output_dir):
    """Run composite bivariate Moran's I at cancer-cell resolution.

    Uses single_cell.h5ad from morphology assignment. MDK vs mean-z secretory score.
    """
    import scanpy as sc

    results = []
    for sample in samples:
        logger.info(f"[cell] {sample}")
        sc_path = str(MORPH_ASSIGN_DIR / sample / f"{sample}_single_cell.h5ad")
        if not os.path.isfile(sc_path):
            logger.warning(f"  SKIP: {sc_path} not found")
            results.append(
                {
                    "sample": sample,
                    "morans_I": np.nan,
                    "p_value": np.nan,
                    "n_cells": 0,
                    "mdk_nonzero": 0,
                    "significant": False,
                    "n_secretory_genes": 0,
                    "note": "single_cell.h5ad not found",
                }
            )
            continue

        adata = sc.read_h5ad(sc_path)
        coords = adata.obsm["spatial"].copy()

        valid_mask = np.isfinite(coords).all(axis=1)
        if not valid_mask.all():
            logger.info(f"  Dropping {(~valid_mask).sum()} cells with NaN/inf coords")
            adata = adata[valid_mask]
            coords = coords[valid_mask]

        logger.info(f"  Loaded: {adata.shape[0]} cells, {adata.shape[1]} genes")

        if "MDK" not in adata.var_names:
            logger.warning(f"  SKIP: MDK not in gene list")
            results.append(
                {
                    "sample": sample,
                    "morans_I": np.nan,
                    "p_value": np.nan,
                    "n_cells": adata.shape[0],
                    "mdk_nonzero": 0,
                    "significant": False,
                    "n_secretory_genes": 0,
                    "note": "MDK not in single_cell var",
                }
            )
            continue

        mdk = _extract_gene_vector(adata, "MDK")
        mdk_nonzero = int(np.sum(mdk > 0))

        if mdk_nonzero < 10:
            logger.info(f"  SKIP: only {mdk_nonzero} nonzero MDK")
            results.append(
                {
                    "sample": sample,
                    "morans_I": 0.0,
                    "p_value": 1.0,
                    "n_cells": adata.shape[0],
                    "mdk_nonzero": mdk_nonzero,
                    "significant": False,
                    "n_secretory_genes": 0,
                    "note": f"too few nonzero ({mdk_nonzero})",
                }
            )
            continue

        present = [g for g in SECRETORY_GENES if g in adata.var_names]
        logger.info(f"  Secretory genes present: {len(present)}/{len(SECRETORY_GENES)}")

        sec_x = adata[:, present].X
        if hasattr(sec_x, "toarray"):
            sec_mat = sec_x.toarray().astype(float)
        else:
            sec_mat = np.asarray(sec_x).astype(float)
        sec_z = zscore(sec_mat, axis=0, nan_policy="omit")
        sec_z = np.nan_to_num(sec_z, 0)
        sec_score = sec_z.mean(axis=1)

        try:
            I_val, p_val = compute_bivariate_moran(mdk, sec_score, coords)
        except Exception as e:
            logger.error(f"  ERROR computing Moran's I: {e}")
            results.append(
                {
                    "sample": sample,
                    "morans_I": np.nan,
                    "p_value": np.nan,
                    "n_cells": adata.shape[0],
                    "mdk_nonzero": mdk_nonzero,
                    "significant": False,
                    "n_secretory_genes": len(present),
                    "note": f"computation error: {e}",
                }
            )
            continue

        sig = p_val <= 0.001
        logger.info(f"  Moran's I = {I_val:.4f}, p = {p_val:.4f}, significant = {sig}")

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
    out_path = os.path.join(output_dir, "cancer_cell_morans.csv")
    df.to_csv(out_path, index=False)
    logger.info(f"Cancer-cell results saved to {out_path}")

    sig_mask = df["significant"] == True
    logger.info(f"{sig_mask.sum()}/{len(df)} samples significant (p <= 0.001)")
    if sig_mask.any():
        logger.info(f"Mean I (significant): {df.loc[sig_mask, 'morans_I'].mean():.3f}")

    return df


def compute_trajectory_deltas(cell_df, pairs_path, output_dir):
    """Compute biopsy→surgery deltas in Moran's I using sample_pairs.csv."""
    pairs = load_pairs(pairs_path)
    deltas = []
    for _, pair in pairs.iterrows():
        biopsy_full = SHORT_TO_FULL.get(pair["biopsy_sample"], pair["biopsy_sample"])
        surgery_full = SHORT_TO_FULL.get(pair["surgery_sample"], pair["surgery_sample"])
        biopsy_row = cell_df[cell_df["sample"] == biopsy_full]
        surgery_row = cell_df[cell_df["sample"] == surgery_full]
        if len(biopsy_row) == 0 or len(surgery_row) == 0:
            continue
        b_I = biopsy_row.iloc[0]["morans_I"]
        s_I = surgery_row.iloc[0]["morans_I"]
        if pd.isna(b_I) or pd.isna(s_I):
            continue
        deltas.append(
            {
                "patient": pair["patient_id"],
                "biopsy_sample": pair["biopsy_sample"],
                "surgery_sample": pair["surgery_sample"],
                "biopsy_I": round(b_I, 4),
                "surgery_I": round(s_I, 4),
                "delta_I": round(s_I - b_I, 4),
            }
        )

    delta_df = pd.DataFrame(deltas)
    if len(delta_df) > 0:
        out_path = os.path.join(output_dir, "trajectory_deltas.csv")
        delta_df.to_csv(out_path, index=False)
        logger.info(f"Trajectory deltas saved to {out_path}")
        logger.info(f"Mean delta: {delta_df['delta_I'].mean():.4f}")
    else:
        logger.warning("No trajectory deltas computed (missing paired samples)")

    return delta_df


def main():
    parser = argparse.ArgumentParser(description="Unified spatial Moran's I analysis")
    parser.add_argument(
        "--output-dir",
        type=str,
        default=str(PROJECT_DIR / "mdk_analysis" / "outputs" / "spatial"),
    )
    parser.add_argument(
        "--manifest",
        type=str,
        default=str(PROJECT_DIR / "mdk_analysis" / "data" / "sample_manifest.csv"),
    )
    parser.add_argument(
        "--pairs",
        type=str,
        default=str(PROJECT_DIR / "mdk_analysis" / "data" / "sample_pairs.csv"),
    )
    parser.add_argument(
        "--skip-spot",
        action="store_true",
        help="Skip spot-level per-gene analysis (slow)",
    )
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    samples = load_manifest(args.manifest)
    logger.info(f"Loaded {len(samples)} samples from manifest")

    if not args.skip_spot:
        logger.info("=== Spot-level per-gene Moran's I ===")
        run_spot_level(samples, args.output_dir)

    logger.info("=== Cancer-cell-level composite Moran's I ===")
    cell_df = run_cancer_cell_level(samples, args.output_dir)

    logger.info("=== Trajectory deltas ===")
    compute_trajectory_deltas(cell_df, args.pairs, args.output_dir)

    logger.info("=== Spatial Moran's I analysis complete ===")


if __name__ == "__main__":
    main()
