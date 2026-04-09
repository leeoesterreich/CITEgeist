#!/usr/bin/env python
"""Compare MDK bivariate Moran's I under two cancer-cell grouping strategies.

Strategy A (current): Cancer_Basal only — matches the existing manuscript numbers
Strategy B (proposed): Cancer_Basal + Cancer_Luminal merged as "Epithelial"

For each sample, computes the composite bivariate Moran's I (MDK vs secretory
score) under both strategies and reports the delta.  If the mean absolute delta
across samples exceeds RERUN_THRESHOLD (default 0.02), the script exits with a
non-zero code and a recommendation to do a full pipeline re-run with a single
Epithelial QP type.

Also uses func_KRT5_*_gate (from SACE protein run) to verify that KRT5+ cells
within the merged Epithelial population reproduce the Cancer_Basal signal.

Usage:
    python analyze_mdk_epithelial_comparison.py
    python analyze_mdk_epithelial_comparison.py --rerun-threshold 0.01
"""

import argparse
import logging
import sys
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import scanpy as sc
from esda import Moran_BV
from libpysal.weights import KNN

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)

# =============================================================================
# Configuration
# =============================================================================
REPO = Path(__file__).resolve().parents[2]

V3_DIR = REPO / "output" / "mdk_v3_reanalysis_expanded"
MORPH_DIR = REPO / "output" / "morphology_assignment_v3"
SACE_PROTEIN_DIR = REPO / "output" / "sace_protein_12patient"
OUTPUT_DIR = REPO / "output" / "mdk_epithelial_comparison"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

SAMPLES = [
    "HCC22-088-P1-S1", "HCC22-088-P1-S2",
    "HCC22-088-P2-S1", "HCC22-088-P2-S2",
    "HCC22-088-P3-S1_A", "HCC22-088-P3-S2",
    "HCC22-088-P4-S1", "HCC22-088-P4-S2_1i_rep",
    "HCC22-088-P5-S1", "HCC22-088-P5-S2_F_rep",
    "HCC22-088-P6-S1", "HCC22-088-P6-S2_D",
]

PATIENT_PAIRS = {
    "P1": ("HCC22-088-P1-S1", "HCC22-088-P1-S2"),
    "P2": ("HCC22-088-P2-S1", "HCC22-088-P2-S2"),
    "P3": ("HCC22-088-P3-S1_A", "HCC22-088-P3-S2"),
    "P4": ("HCC22-088-P4-S1", "HCC22-088-P4-S2_1i_rep"),
    "P5": ("HCC22-088-P5-S1", "HCC22-088-P5-S2_F_rep"),
    "P6": ("HCC22-088-P6-S1", "HCC22-088-P6-S2_D"),
}

SECRETORY_GENES = [
    "HSP90B1", "HSPA5", "CALR", "CANX", "PDIA4",
    "PDIA6", "SEC23A", "SEC61B", "ATF6", "MAN1A1", "XBP1",
]

KNN_K = 6
N_PERM = 999
RERUN_THRESHOLD = 0.02  # Mean |delta| that triggers full-rerun recommendation


# =============================================================================
# Helpers
# =============================================================================

def _load_spatial_coords(sample_name: str) -> pd.DataFrame:
    """Load spot spatial coords from tissue_positions.csv without loading images."""
    pos_path = (
        Path("/ix1/alee/LO_LAB/General/Lab_Data/"
             "20250210_CITEGeistPublicData_GEO_Alex/processed_files")
        / sample_name / "outs" / "spatial" / "tissue_positions.csv"
    )
    df = pd.read_csv(pos_path, index_col=0)
    return df[["pxl_col_in_fullres", "pxl_row_in_fullres"]].rename(
        columns={"pxl_col_in_fullres": "x", "pxl_row_in_fullres": "y"}
    )


def _aggregate_to_spots(cell_adata, type_mask, coords_df):
    """Aggregate per-cell GEX to spot level for a given cell type mask.

    Returns (spot_X, spot_barcodes, spot_coords, gene_names) or None if < 20 spots.
    """
    subset = cell_adata[type_mask].copy()
    if subset.n_obs == 0:
        return None

    X = subset.X.toarray() if hasattr(subset.X, "toarray") else np.asarray(subset.X)
    spot_ids = subset.obs["spot_barcode"].values
    unique_spots = np.unique(spot_ids)

    spot_sums = np.zeros((len(unique_spots), X.shape[1]))
    for i, spot in enumerate(unique_spots):
        m = spot_ids == spot
        spot_sums[i] = X[m].sum(axis=0)

    # Align spatial coords
    valid = [s for s in unique_spots if s in coords_df.index]
    if len(valid) < 20:
        return None

    valid = np.array(valid)
    keep = np.array([np.where(unique_spots == s)[0][0] for s in valid])
    spot_sums = spot_sums[keep]
    coords = coords_df.loc[valid, ["x", "y"]].values.astype(float)

    return spot_sums, valid, coords, list(cell_adata.var_names)


def _bivariate_morans_I(x_vals, y_vals, coords):
    """Compute bivariate Moran's I with permutation p-value."""
    if np.std(x_vals) == 0 or np.std(y_vals) == 0:
        return np.nan, np.nan
    w = KNN.from_array(coords, k=KNN_K)
    w.transform = "r"
    bv = Moran_BV(x_vals, y_vals, w, permutations=N_PERM)
    return float(bv.I), float(bv.p_sim)


def _secretory_score(X, gene_names):
    """Mean log1p expression across secretory genes present in matrix."""
    import scanpy as sc
    import anndata
    adata = anndata.AnnData(X=X.copy(), var=pd.DataFrame(index=gene_names))
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    present = [g for g in SECRETORY_GENES if g in gene_names]
    if not present:
        return None
    idx = [gene_names.index(g) for g in present]
    Xn = adata.X.toarray() if hasattr(adata.X, "toarray") else np.asarray(adata.X)
    return Xn[:, idx].mean(axis=1)


def _mdk_vector(X, gene_names):
    """Log1p-normalized MDK expression vector."""
    if "MDK" not in gene_names:
        return None
    import anndata, scanpy as sc
    adata = anndata.AnnData(X=X.copy(), var=pd.DataFrame(index=gene_names))
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    Xn = adata.X.toarray() if hasattr(adata.X, "toarray") else np.asarray(adata.X)
    mdk_idx = gene_names.index("MDK")
    return Xn[:, mdk_idx]


# =============================================================================
# Per-sample analysis
# =============================================================================

def analyze_sample(sample_name: str) -> dict | None:
    """Compute composite bivariate Moran's I under both grouping strategies.

    Returns dict with keys: sample, I_basal_only, p_basal_only,
    I_epithelial, p_epithelial, delta, n_spots_basal, n_spots_epithelial,
    krt5_gate_available, I_krt5_pos (Moran's I for KRT5+ cells only).
    """
    sc_path = V3_DIR / sample_name / f"{sample_name}_single_cell.h5ad"
    if not sc_path.exists():
        logger.warning("Missing v3 h5ad for %s", sample_name)
        return None

    cell_adata = sc.read_h5ad(sc_path)

    # Load KRT5 gate columns from the updated morphology h5ad (written by SACE protein run)
    morph_path = MORPH_DIR / sample_name / f"{sample_name}_single_cell.h5ad"
    krt5_gate = None
    if morph_path.exists():
        morph_adata = sc.read_h5ad(morph_path)
        krt5_cols = [c for c in morph_adata.obs.columns if "KRT5" in c and c.endswith("_gate")]
        if krt5_cols:
            # Merge KRT5 gate into cell_adata via cell_id alignment
            krt5_series = morph_adata.obs[krt5_cols].max(axis=1)  # any KRT5 gate positive
            krt5_gate = krt5_series.reindex(cell_adata.obs_names)
            logger.info("  KRT5 gates available: %s", krt5_cols)
        else:
            logger.warning("  No KRT5 gate cols in %s (SACE protein may not have run yet)", sample_name)

    # Load spatial coords
    try:
        coords_df = _load_spatial_coords(sample_name)
    except Exception as e:
        logger.error("  Could not load spatial coords for %s: %s", sample_name, e)
        return None

    # -----------------------------------------------------------------------
    # Strategy A: Cancer_Basal only
    # -----------------------------------------------------------------------
    mask_basal = cell_adata.obs["cell_type"] == "Cancer_Basal"
    result_a = _aggregate_to_spots(cell_adata, mask_basal, coords_df)
    if result_a is None:
        logger.warning("  Too few Cancer_Basal spots — skipping %s", sample_name)
        return None
    X_a, spots_a, coords_a, genes = result_a
    mdk_a = _mdk_vector(X_a, genes)
    sec_a = _secretory_score(X_a, genes)
    if mdk_a is None or sec_a is None:
        logger.warning("  MDK or secretory genes missing (strategy A) — skipping %s", sample_name)
        return None
    I_a, p_a = _bivariate_morans_I(mdk_a, sec_a, coords_a)

    # -----------------------------------------------------------------------
    # Strategy B: Cancer_Basal + Cancer_Luminal merged (Epithelial)
    # -----------------------------------------------------------------------
    mask_epi = cell_adata.obs["cell_type"].isin(["Cancer_Basal", "Cancer_Luminal"])
    result_b = _aggregate_to_spots(cell_adata, mask_epi, coords_df)
    if result_b is None:
        logger.warning("  Too few Epithelial spots — skipping %s", sample_name)
        return None
    X_b, spots_b, coords_b, _ = result_b
    mdk_b = _mdk_vector(X_b, genes)
    sec_b = _secretory_score(X_b, genes)
    if mdk_b is None or sec_b is None:
        logger.warning("  MDK or secretory genes missing (strategy B) — skipping %s", sample_name)
        return None
    I_b, p_b = _bivariate_morans_I(mdk_b, sec_b, coords_b)

    # -----------------------------------------------------------------------
    # Strategy C (optional): KRT5+ cells within merged Epithelial
    # -----------------------------------------------------------------------
    I_c, p_c = np.nan, np.nan
    if krt5_gate is not None:
        krt5_pos = krt5_gate.fillna(0).astype(bool)
        mask_krt5 = mask_epi & krt5_pos
        result_c = _aggregate_to_spots(cell_adata, mask_krt5, coords_df)
        if result_c is not None:
            X_c, spots_c, coords_c, _ = result_c
            mdk_c = _mdk_vector(X_c, genes)
            sec_c = _secretory_score(X_c, genes)
            if mdk_c is not None and sec_c is not None:
                I_c, p_c = _bivariate_morans_I(mdk_c, sec_c, coords_c)
                logger.info(
                    "  KRT5+ Epithelial: %d spots, I=%.4f, p=%.3f",
                    len(spots_c), I_c, p_c,
                )

    delta = float(I_b - I_a) if not (np.isnan(I_a) or np.isnan(I_b)) else np.nan
    logger.info(
        "  A (Basal-only): n=%d spots, I=%.4f, p=%.3f",
        len(spots_a), I_a, p_a,
    )
    logger.info(
        "  B (Epithelial): n=%d spots, I=%.4f, p=%.3f  | delta=%.4f",
        len(spots_b), I_b, p_b, delta,
    )

    return {
        "sample": sample_name,
        "I_basal_only": I_a,
        "p_basal_only": p_a,
        "n_spots_basal": len(spots_a),
        "I_epithelial": I_b,
        "p_epithelial": p_b,
        "n_spots_epithelial": len(spots_b),
        "delta": delta,                      # positive = epithelial > basal-only
        "krt5_gate_available": krt5_gate is not None and len([c for c in (morph_adata.obs.columns if morph_path.exists() else []) if "KRT5" in c]) > 0,
        "I_krt5_pos": I_c,
        "p_krt5_pos": p_c,
    }


# =============================================================================
# Summary + trajectory
# =============================================================================

def trajectory_delta(rows: list[dict], strategy: str) -> pd.DataFrame:
    """Compute biopsy→surgery trajectory delta for a given strategy."""
    I_col = "I_basal_only" if strategy == "A" else "I_epithelial"
    df = pd.DataFrame(rows).set_index("sample")
    records = []
    for patient, (s1, s2) in PATIENT_PAIRS.items():
        if s1 in df.index and s2 in df.index:
            I1 = df.loc[s1, I_col]
            I2 = df.loc[s2, I_col]
            records.append({
                "patient": patient,
                "biopsy": s1,
                "surgery": s2,
                "I_biopsy": I1,
                "I_surgery": I2,
                "traj_delta": I2 - I1,
                "strategy": strategy,
            })
    return pd.DataFrame(records)


# =============================================================================
# Main
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="MDK Moran's I: Cancer_Basal-only vs merged Epithelial comparison"
    )
    parser.add_argument(
        "--rerun-threshold",
        type=float,
        default=RERUN_THRESHOLD,
        help="Mean |delta| above which a full re-run is recommended (default: 0.02)",
    )
    args = parser.parse_args()

    rows = []
    for sample in SAMPLES:
        logger.info("=== %s ===", sample)
        result = analyze_sample(sample)
        if result is not None:
            rows.append(result)

    if not rows:
        logger.error("No samples processed — check input paths.")
        sys.exit(1)

    df = pd.DataFrame(rows)
    out_csv = OUTPUT_DIR / "mdk_epithelial_comparison.csv"
    df.to_csv(out_csv, index=False)
    logger.info("Results written to %s", out_csv)

    # --- Summary table ---
    print("\n" + "=" * 80)
    print("MDK Bivariate Moran's I: Cancer_Basal-only (A) vs Merged Epithelial (B)")
    print("=" * 80)
    print(df[["sample", "n_spots_basal", "I_basal_only", "p_basal_only",
              "n_spots_epithelial", "I_epithelial", "p_epithelial", "delta",
              "I_krt5_pos"]].to_string(index=False, float_format="{:.4f}".format))

    # Aggregate stats
    valid = df.dropna(subset=["delta"])
    mean_delta = valid["delta"].mean()
    mean_I_a = valid["I_basal_only"].mean()
    mean_I_b = valid["I_epithelial"].mean()
    n_sig_a = (valid["p_basal_only"] < 0.05).sum()
    n_sig_b = (valid["p_epithelial"] < 0.05).sum()

    print(f"\nSummary across {len(valid)} samples:")
    print(f"  Strategy A (Basal-only):  mean I = {mean_I_a:.4f}, N significant (p<0.05) = {n_sig_a}/{len(valid)}")
    print(f"  Strategy B (Epithelial):  mean I = {mean_I_b:.4f}, N significant (p<0.05) = {n_sig_b}/{len(valid)}")
    print(f"  Mean delta (B - A):       {mean_delta:+.4f}  (positive = merged is stronger)")

    krt5_rows = valid.dropna(subset=["I_krt5_pos"])
    if not krt5_rows.empty:
        mean_I_c = krt5_rows["I_krt5_pos"].mean()
        print(f"  Strategy C (KRT5+ only): mean I = {mean_I_c:.4f}, N samples = {len(krt5_rows)}/{len(valid)}")

    # --- Trajectory comparison ---
    traj_a = trajectory_delta(rows, "A")
    traj_b = trajectory_delta(rows, "B")
    if not traj_a.empty:
        traj = pd.merge(
            traj_a[["patient", "I_biopsy", "I_surgery", "traj_delta"]].rename(
                columns={"I_biopsy": "I_biopsy_A", "I_surgery": "I_surgery_A", "traj_delta": "delta_A"}),
            traj_b[["patient", "I_biopsy", "I_surgery", "traj_delta"]].rename(
                columns={"I_biopsy": "I_biopsy_B", "I_surgery": "I_surgery_B", "traj_delta": "delta_B"}),
            on="patient",
        )
        print("\nTrajectory (biopsy → surgery):")
        print(traj[["patient", "delta_A", "delta_B"]].to_string(index=False, float_format="{:+.4f}".format))
        traj.to_csv(OUTPUT_DIR / "mdk_trajectory_comparison.csv", index=False)

    # --- Re-run recommendation ---
    print("\n" + "-" * 80)
    abs_delta = valid["delta"].abs().mean()
    if abs_delta > args.rerun_threshold:
        print(
            f"RECOMMENDATION: Full re-run warranted.\n"
            f"  Mean |delta| = {abs_delta:.4f} exceeds threshold {args.rerun_threshold:.4f}.\n"
            f"  Re-run QP with single Epithelial type (EPCAM only, KRT5 freed as functional marker)."
        )
        sys.exit(2)  # Non-zero exit signals caller that re-run is needed
    else:
        print(
            f"RECOMMENDATION: Post-hoc merge is acceptable.\n"
            f"  Mean |delta| = {abs_delta:.4f} is within threshold {args.rerun_threshold:.4f}.\n"
            f"  Cancer_Basal and merged Epithelial give equivalent MDK signal.\n"
            f"  Rename types to Epithelial in downstream analysis without re-running QP."
        )


if __name__ == "__main__":
    main()
