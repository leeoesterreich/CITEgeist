#!/usr/bin/env python
"""
Per-gene bivariate Moran's I: MDK vs each of 11 secretory genes across 12 samples.

132 tests total, FDR-corrected. Produces:
1. CSV with all per-gene per-sample Moran's I values + p-values + FDR q-values
2. Summary: how many genes significant per sample, how many samples per gene
3. Negative control: 11 random non-secretory genes for comparison

Output: mdk_analysis/if_analysis/pergene_bivariate_morans.csv
        mdk_analysis/if_analysis/pergene_morans_summary.txt
"""

import os
import numpy as np
import pandas as pd
import squidpy as sq
from esda.moran import Moran_BV
from libpysal.weights import KNN
from statsmodels.stats.multitest import multipletests
import warnings

warnings.filterwarnings("ignore")

PROJECT_DIR = "/path/to/CITEgeist_analysis"
DATA_ROOT = "/path/to/CITEgeist_public_data/processed_files"
OUTPUT_DIR = os.path.join(PROJECT_DIR, "mdk_analysis/outputs/if_analysis")

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


def load_sample(sample):
    """Load Visium data, filter to GEX, deduplicate, return adata + coords."""
    visium_dir = os.path.join(DATA_ROOT, sample, "outs")
    if not os.path.isdir(visium_dir):
        return None, None

    adata = sq.read.visium(visium_dir, counts_file="filtered_feature_bc_matrix.h5", load_images=True)
    coords = adata.obsm["spatial"].copy()

    # Drop NaN/inf coords
    valid = np.isfinite(coords).all(axis=1)
    if not valid.all():
        adata = adata[valid]
        coords = coords[valid]

    # GEX only
    if "feature_types" in adata.var.columns:
        adata = adata[:, adata.var["feature_types"] == "Gene Expression"]

    # Deduplicate
    dup = adata.var_names.duplicated(keep="first")
    if dup.any():
        adata = adata[:, ~dup]

    return adata, coords


def compute_moran_bv(x, y, coords, n_perm=N_PERM):
    """Compute bivariate Moran's I using esda with KNN(k=6)."""
    k = min(6, len(x) - 1)
    w = KNN.from_array(coords, k=k)
    w.transform = "r"
    bv = Moran_BV(x, y, w, permutations=n_perm)
    return float(bv.I), float(bv.p_sim)


def main():
    rows = []

    for sample in SAMPLES:
        print(f"\n=== {sample} ===")
        adata, coords = load_sample(sample)
        if adata is None:
            print("  SKIP: not found")
            continue

        # MDK expression
        if "MDK" not in adata.var_names:
            print("  SKIP: no MDK")
            continue
        mdk = np.asarray(adata[:, "MDK"].X.todense()).ravel().astype(float)
        mdk_nz = int(np.sum(mdk > 0))
        print(f"  {adata.shape[0]} spots, MDK nonzero: {mdk_nz}")

        if mdk_nz < 10:
            print("  SKIP: too few MDK spots")
            continue

        # Build spatial weights once per sample
        np.random.seed(42)
        w = KNN.from_array(coords, k=min(6, len(mdk) - 1))
        w.transform = "r"

        # Test each secretory gene
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

            gene_expr = np.asarray(adata[:, gene].X.todense()).ravel().astype(float)
            gene_nz = int(np.sum(gene_expr > 0))

            try:
                bv = Moran_BV(mdk, gene_expr, w, permutations=N_PERM)
                I_val, p_val = float(bv.I), float(bv.p_sim)
            except Exception as e:
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
            print(f"  {gene}: I={I_val:.4f} p={p_val:.4f} (nz={gene_nz})")

        # Negative controls: 11 random non-secretory, non-MDK genes
        rng = np.random.RandomState(123)
        non_sec = [g for g in adata.var_names if g not in SECRETORY_GENES and g != "MDK"]
        ctrl_genes = rng.choice(non_sec, size=min(11, len(non_sec)), replace=False)
        for gene in ctrl_genes:
            gene_expr = np.asarray(adata[:, gene].X.todense()).ravel().astype(float)
            gene_nz = int(np.sum(gene_expr > 0))
            try:
                bv = Moran_BV(mdk, gene_expr, w, permutations=N_PERM)
                I_val, p_val = float(bv.I), float(bv.p_sim)
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

    # Build dataframe
    df = pd.DataFrame(rows)

    # FDR correction across ALL secretory tests
    sec_mask = (df["gene_type"] == "secretory") & df["p_value"].notna()
    _, fdr_q, _, _ = multipletests(df.loc[sec_mask, "p_value"].values, method="fdr_bh")
    df.loc[sec_mask, "fdr_q"] = np.round(fdr_q, 4)

    # FDR for controls separately
    ctrl_mask = (df["gene_type"] == "control") & df["p_value"].notna()
    if ctrl_mask.sum() > 0:
        _, fdr_q_ctrl, _, _ = multipletests(df.loc[ctrl_mask, "p_value"].values, method="fdr_bh")
        df.loc[ctrl_mask, "fdr_q"] = np.round(fdr_q_ctrl, 4)

    # Save full results
    out_path = os.path.join(OUTPUT_DIR, "pergene_bivariate_morans.csv")
    df.to_csv(out_path, index=False)
    print(f"\n\nSaved to {out_path}")

    # Summary
    sec_df = df[df["gene_type"] == "secretory"].copy()
    print("\n" + "=" * 80)
    print("SECRETORY GENES: Per-gene bivariate Moran's I with MDK")
    print("=" * 80)

    # Per-gene summary
    print("\nPer-gene (across samples):")
    for gene in SECRETORY_GENES:
        gdf = sec_df[sec_df["gene"] == gene]
        sig = (gdf["fdr_q"] < 0.05).sum()
        mean_I = gdf["morans_I"].mean()
        print(f"  {gene:>8}: mean I={mean_I:.3f}, significant in {sig}/{len(gdf)} samples (FDR<0.05)")

    # Per-sample summary
    print("\nPer-sample (across genes):")
    for sample in SAMPLES:
        sdf = sec_df[sec_df["sample"] == sample]
        if len(sdf) == 0:
            continue
        sig = (sdf["fdr_q"] < 0.05).sum()
        mean_I = sdf["morans_I"].mean()
        print(f"  {sample}: mean I={mean_I:.3f}, {sig}/{len(sdf)} genes significant (FDR<0.05)")

    # Overall
    n_sig = (sec_df["fdr_q"] < 0.05).sum()
    n_total = sec_df["fdr_q"].notna().sum()
    print(f"\nOverall: {n_sig}/{n_total} secretory gene-sample pairs significant (FDR<0.05)")
    print(f"Mean I (secretory): {sec_df['morans_I'].mean():.3f}")

    # Control comparison
    ctrl_df = df[df["gene_type"] == "control"]
    n_sig_ctrl = (ctrl_df["fdr_q"] < 0.05).sum()
    n_total_ctrl = ctrl_df["fdr_q"].notna().sum()
    print(f"\nControl: {n_sig_ctrl}/{n_total_ctrl} control gene-sample pairs significant (FDR<0.05)")
    print(f"Mean I (control): {ctrl_df['morans_I'].mean():.3f}")

    # Save summary
    summary_path = os.path.join(OUTPUT_DIR, "pergene_morans_summary.txt")
    with open(summary_path, "w") as f:
        f.write("Per-gene bivariate Moran's I: MDK vs secretory genes\n")
        f.write(f"Method: esda.Moran_BV, KNN(k=6), {N_PERM} permutations, all spots\n")
        f.write(f"FDR correction: Benjamini-Hochberg across {n_total} secretory tests\n\n")
        f.write(f"Secretory: {n_sig}/{n_total} pairs significant (FDR<0.05)\n")
        f.write(f"Control: {n_sig_ctrl}/{n_total_ctrl} pairs significant (FDR<0.05)\n")
        f.write(f"Mean I secretory: {sec_df['morans_I'].mean():.3f}\n")
        f.write(f"Mean I control: {ctrl_df['morans_I'].mean():.3f}\n")
    print(f"\nSummary saved to {summary_path}")


if __name__ == "__main__":
    main()
