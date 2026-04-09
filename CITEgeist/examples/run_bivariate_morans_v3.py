#!/usr/bin/env python
"""
Recompute bivariate Moran's I on v3 single-cell data.

Two analyses:
  1. Per-gene bivariate Moran's I (MDK vs each of 11 secretory genes) — matches manuscript
  2. Composite bivariate Moran's I (MDK vs secretory score) — for trajectory analysis

Operates on cancer-cell spot-level aggregates from v3 SACE-refilter h5ad files.

Usage:
  python run_bivariate_morans_v3.py           # run all 12 samples
  python run_bivariate_morans_v3.py --sample 0 # run single sample (array jobs)
"""

import argparse
import logging
import sys
import warnings
from pathlib import Path

warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq
from esda import Moran_BV
from libpysal.weights import KNN

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# =============================================================================
# Configuration
# =============================================================================
V3_EXPANDED_DIR = Path("/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/output/mdk_v3_reanalysis_expanded")
OUTPUT_DIR = Path("/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/output/bivariate_morans_v3")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

PATIENT_DATA_ROOT = Path(
    "/ix1/alee/LO_LAB/General/Lab_Data/"
    "20250210_CITEGeistPublicData_GEO_Alex/processed_files"
)

SAMPLES = [
    "HCC22-088-P1-S1", "HCC22-088-P1-S2",
    "HCC22-088-P2-S1", "HCC22-088-P2-S2",
    "HCC22-088-P3-S1_A", "HCC22-088-P3-S2",
    "HCC22-088-P4-S1", "HCC22-088-P4-S2_1i_rep",
    "HCC22-088-P5-S1", "HCC22-088-P5-S2_F_rep",
    "HCC22-088-P6-S1", "HCC22-088-P6-S2_D",
]

# Patient pairing for trajectory analysis
PATIENT_PAIRS = {
    "P1": ("HCC22-088-P1-S1", "HCC22-088-P1-S2"),
    "P2": ("HCC22-088-P2-S1", "HCC22-088-P2-S2"),
    "P3": ("HCC22-088-P3-S1_A", "HCC22-088-P3-S2"),
    "P4": ("HCC22-088-P4-S1", "HCC22-088-P4-S2_1i_rep"),
    "P5": ("HCC22-088-P5-S1", "HCC22-088-P5-S2_F_rep"),
    "P6": ("HCC22-088-P6-S1", "HCC22-088-P6-S2_D"),
}

# 11 secretory genes from manuscript (T1-2 in evidence ledger)
SECRETORY_GENES_MANUSCRIPT = [
    "HSP90B1", "HSPA5", "CALR", "CANX", "PDIA4",
    "PDIA6", "SEC23A", "SEC61B", "ATF6", "MAN1A1", "XBP1",
]

# Control genes (housekeeping, for comparison — same count as secretory)
CONTROL_GENES = [
    "ACTB", "GAPDH", "RPL13", "RPS18", "TUBB",
    "EEF1A1", "UBC", "ALDOA", "ENO1", "TPI1", "LDHA",
]

# Both QP subtypes pooled as a single Epithelial population for Moran's I
CANCER_TYPES = ['Cancer_Basal', 'Cancer_Luminal']
N_PERM = 999
KNN_K = 6  # match original methodology


def load_spot_spatial(sample_name):
    """Load raw Visium to get spot spatial coordinates."""
    patient_dir = PATIENT_DATA_ROOT / sample_name / "outs"
    adata = sq.read.visium(
        str(patient_dir),
        counts_file='filtered_feature_bc_matrix.h5',
        load_images=True,
        gex_only=True,
    )
    return adata


def aggregate_cancer_to_spots(cell_adata, spot_coords_adata):
    """Aggregate per-cell cancer GEX back to spot level.

    Returns spot-level AnnData with only spots that have cancer cells,
    using spatial coords from the Visium data.
    """
    cancer_mask = cell_adata.obs['cell_type'].isin(CANCER_TYPES)
    cancer_cells = cell_adata[cancer_mask].copy()

    if cancer_cells.n_obs == 0:
        return None

    # Get raw counts
    X = cancer_cells.X.toarray() if hasattr(cancer_cells.X, 'toarray') else np.asarray(cancer_cells.X)

    # Sum per spot
    spot_ids = cancer_cells.obs['spot_barcode'].values
    unique_spots = np.unique(spot_ids)

    spot_sums = np.zeros((len(unique_spots), X.shape[1]))
    spot_counts = np.zeros(len(unique_spots), dtype=int)
    for i, spot in enumerate(unique_spots):
        mask = spot_ids == spot
        spot_sums[i] = X[mask].sum(axis=0)
        spot_counts[i] = mask.sum()

    # Build spot-level AnnData
    obs = pd.DataFrame({
        'spot_barcode': unique_spots,
        'n_cancer_cells': spot_counts,
    }, index=unique_spots)

    import anndata
    spot_adata = anndata.AnnData(
        X=spot_sums,
        obs=obs,
        var=cancer_cells.var.copy(),
    )

    # Attach spatial coords from Visium
    common = spot_adata.obs_names.intersection(spot_coords_adata.obs_names)
    spot_adata = spot_adata[common].copy()
    spot_adata.obsm['spatial'] = spot_coords_adata[common].obsm['spatial'].copy()

    return spot_adata


def compute_bivariate_morans(spot_adata, gene_x, gene_y, n_perm=N_PERM, knn_k=KNN_K):
    """Compute bivariate Moran's I between two genes at spot level."""
    if gene_x not in spot_adata.var_names or gene_y not in spot_adata.var_names:
        return {'I': np.nan, 'p_value': np.nan, 'present': False}

    X = spot_adata.X if not hasattr(spot_adata.X, 'toarray') else spot_adata.X.toarray()
    x_idx = list(spot_adata.var_names).index(gene_x)
    y_idx = list(spot_adata.var_names).index(gene_y)

    x_vals = X[:, x_idx].flatten().astype(float)
    y_vals = X[:, y_idx].flatten().astype(float)

    # Skip if either has zero variance
    if np.std(x_vals) == 0 or np.std(y_vals) == 0:
        return {'I': 0.0, 'p_value': 1.0, 'present': True, 'note': 'zero_variance'}

    coords = spot_adata.obsm['spatial']
    w = KNN.from_array(coords, k=knn_k)
    w.transform = 'r'

    bv = Moran_BV(x_vals, y_vals, w, permutations=n_perm)
    return {
        'I': float(bv.I),
        'p_value': float(bv.p_sim),
        'present': True,
    }


def analyze_sample(sample_name):
    """Run bivariate Moran's I for one sample."""
    logger.info(f"=== {sample_name} ===")

    # Load v3 single-cell h5ad
    sc_path = V3_EXPANDED_DIR / sample_name / f"{sample_name}_single_cell.h5ad"
    if not sc_path.exists():
        logger.error(f"  Missing: {sc_path}")
        return None
    cell_adata = sc.read_h5ad(sc_path)
    logger.info(f"  Loaded {cell_adata.n_obs} cells x {cell_adata.n_vars} genes")

    # Load spot spatial coords
    spot_adata_raw = load_spot_spatial(sample_name)

    # Aggregate cancer cells to spots
    spot_adata = aggregate_cancer_to_spots(cell_adata, spot_adata_raw)
    if spot_adata is None or spot_adata.n_obs < 20:
        logger.warning(f"  Too few cancer spots ({spot_adata.n_obs if spot_adata else 0}) — skipping")
        return None

    n_spots = spot_adata.n_obs
    mdk_present = "MDK" in spot_adata.var_names
    logger.info(f"  Cancer spots: {n_spots}, MDK present: {mdk_present}")

    if not mdk_present:
        logger.warning(f"  MDK not in spot-aggregated genes — skipping")
        return None

    # Normalize (log1p for Moran's I, matching standard practice)
    sc.pp.normalize_total(spot_adata, target_sum=1e4)
    sc.pp.log1p(spot_adata)

    # --- Per-gene bivariate Moran's I (MDK vs each secretory gene) ---
    per_gene_results = []
    for gene in SECRETORY_GENES_MANUSCRIPT:
        result = compute_bivariate_morans(spot_adata, "MDK", gene)
        per_gene_results.append({
            'sample': sample_name,
            'gene': gene,
            'gene_type': 'secretory',
            'morans_I': result['I'],
            'p_value': result['p_value'],
            'present': result['present'],
        })
        if result['present']:
            logger.info(f"    MDK vs {gene}: I={result['I']:.4f}, p={result['p_value']:.3f}")

    # --- Control genes ---
    for gene in CONTROL_GENES:
        result = compute_bivariate_morans(spot_adata, "MDK", gene)
        per_gene_results.append({
            'sample': sample_name,
            'gene': gene,
            'gene_type': 'control',
            'morans_I': result['I'],
            'p_value': result['p_value'],
            'present': result['present'],
        })

    # --- Composite bivariate Moran's I (MDK vs secretory score) ---
    X = spot_adata.X if not hasattr(spot_adata.X, 'toarray') else spot_adata.X.toarray()
    gene_names = list(spot_adata.var_names)

    secretory_present = [g for g in SECRETORY_GENES_MANUSCRIPT if g in gene_names]
    if len(secretory_present) >= 3:
        sec_idx = [gene_names.index(g) for g in secretory_present]
        secretory_score = X[:, sec_idx].mean(axis=1).flatten().astype(float)

        mdk_idx = gene_names.index("MDK")
        mdk_vals = X[:, mdk_idx].flatten().astype(float)

        if np.std(secretory_score) > 0 and np.std(mdk_vals) > 0:
            coords = spot_adata.obsm['spatial']
            w = KNN.from_array(coords, k=KNN_K)
            w.transform = 'r'
            bv = Moran_BV(mdk_vals, secretory_score, w, permutations=N_PERM)
            composite_I = float(bv.I)
            composite_p = float(bv.p_sim)
        else:
            composite_I = 0.0
            composite_p = 1.0
    else:
        composite_I = np.nan
        composite_p = np.nan
        secretory_present = []

    mdk_nonzero_raw = int((spot_adata[:, "MDK"].X.toarray().flatten() > 0).sum()) if "MDK" in spot_adata.var_names else 0

    composite_result = {
        'sample': sample_name,
        'composite_I': composite_I,
        'composite_p': composite_p,
        'n_spots': n_spots,
        'n_cancer_cells': int(spot_adata.obs['n_cancer_cells'].sum()),
        'mdk_nonzero_spots': mdk_nonzero_raw,
        'n_secretory_present': len(secretory_present),
        'secretory_present': ','.join(secretory_present),
    }

    logger.info(f"  Composite: I={composite_I:.4f}, p={composite_p:.4f}")
    logger.info(f"  Secretory genes present: {len(secretory_present)}/{len(SECRETORY_GENES_MANUSCRIPT)}")

    return per_gene_results, composite_result


def run_all_and_summarize(sample_indices=None):
    """Run all samples and produce summary CSVs."""
    if sample_indices is None:
        sample_indices = range(len(SAMPLES))

    all_per_gene = []
    all_composite = []

    for idx in sample_indices:
        sample = SAMPLES[idx]
        result = analyze_sample(sample)
        if result is None:
            continue
        per_gene, composite = result
        all_per_gene.extend(per_gene)
        all_composite.append(composite)

    # Save per-gene results
    per_gene_df = pd.DataFrame(all_per_gene)
    per_gene_df.to_csv(OUTPUT_DIR / "per_gene_bivariate_morans_v3.csv", index=False)
    logger.info(f"Saved per-gene results: {len(per_gene_df)} rows")

    # Save composite results
    composite_df = pd.DataFrame(all_composite)
    composite_df.to_csv(OUTPUT_DIR / "composite_bivariate_morans_v3.csv", index=False)
    logger.info(f"Saved composite results: {len(composite_df)} rows")

    # --- Summary statistics ---
    if len(per_gene_df) == 0:
        logger.warning("No results to summarize")
        return

    # FDR correction on per-gene results (Benjamini-Hochberg)
    from statsmodels.stats.multitest import multipletests
    valid = per_gene_df['present'] & per_gene_df['p_value'].notna()
    if valid.sum() > 0:
        _, fdr, _, _ = multipletests(per_gene_df.loc[valid, 'p_value'], method='fdr_bh')
        per_gene_df.loc[valid, 'fdr'] = fdr
    else:
        per_gene_df['fdr'] = np.nan

    per_gene_df.to_csv(OUTPUT_DIR / "per_gene_bivariate_morans_v3.csv", index=False)

    # Stats for manuscript
    sec_mask = (per_gene_df['gene_type'] == 'secretory') & per_gene_df['present']
    ctrl_mask = (per_gene_df['gene_type'] == 'control') & per_gene_df['present']

    sec_df = per_gene_df[sec_mask]
    ctrl_df = per_gene_df[ctrl_mask]

    n_sec_tests = len(sec_df)
    n_sec_sig = (sec_df['fdr'] < 0.05).sum() if 'fdr' in sec_df.columns else 0
    mean_sec_I = sec_df['morans_I'].mean()

    n_ctrl_tests = len(ctrl_df)
    n_ctrl_sig = (ctrl_df['fdr'] < 0.05).sum() if 'fdr' in ctrl_df.columns else 0
    mean_ctrl_I = ctrl_df['morans_I'].mean()

    # Per-gene: how many samples each secretory gene is significant
    gene_sample_sig = sec_df[sec_df['fdr'] < 0.05].groupby('gene').size()

    # Composite trajectory deltas
    trajectory_rows = []
    for patient, (biopsy, surgery) in PATIENT_PAIRS.items():
        biopsy_row = composite_df[composite_df['sample'] == biopsy]
        surgery_row = composite_df[composite_df['sample'] == surgery]
        if len(biopsy_row) > 0 and len(surgery_row) > 0:
            delta = surgery_row.iloc[0]['composite_I'] - biopsy_row.iloc[0]['composite_I']
            trajectory_rows.append({
                'patient': patient,
                'biopsy_sample': biopsy,
                'surgery_sample': surgery,
                'biopsy_I': float(biopsy_row.iloc[0]['composite_I']),
                'surgery_I': float(surgery_row.iloc[0]['composite_I']),
                'delta_I': float(delta),
            })

    trajectory_df = pd.DataFrame(trajectory_rows)
    trajectory_df.to_csv(OUTPUT_DIR / "trajectory_deltas_v3.csv", index=False)

    # Write human-readable summary
    summary_lines = [
        "=" * 80,
        "BIVARIATE MORAN'S I — v3 SINGLE-CELL REANALYSIS",
        "Cancer cell spot-level aggregates from SACE-refilter h5ad files",
        f"KNN k={KNN_K}, {N_PERM} permutations, log1p normalized",
        "=" * 80,
        "",
        "PER-GENE BIVARIATE MORAN'S I (MDK vs each secretory gene)",
        "-" * 60,
        f"  Total tests (secretory): {n_sec_tests} ({len(SECRETORY_GENES_MANUSCRIPT)} genes x {len(set(sec_df['sample']))} samples)",
        f"  Significant (FDR < 0.05): {n_sec_sig}/{n_sec_tests} ({100*n_sec_sig/max(n_sec_tests,1):.0f}%)",
        f"  Mean Moran's I (secretory): {mean_sec_I:.3f}",
        "",
        f"  Total tests (control): {n_ctrl_tests}",
        f"  Significant (FDR < 0.05): {n_ctrl_sig}/{n_ctrl_tests} ({100*n_ctrl_sig/max(n_ctrl_tests,1):.0f}%)",
        f"  Mean Moran's I (control): {mean_ctrl_I:.3f}",
        "",
        "  Per-gene significance (N samples with FDR < 0.05):",
    ]
    for gene in SECRETORY_GENES_MANUSCRIPT:
        n = gene_sample_sig.get(gene, 0)
        n_total = len(set(sec_df[sec_df['gene'] == gene]['sample']))
        summary_lines.append(f"    {gene:12s}: {n}/{n_total} samples")

    summary_lines.extend([
        "",
        "COMPOSITE BIVARIATE MORAN'S I (MDK vs mean secretory score)",
        "-" * 60,
    ])
    for _, row in composite_df.iterrows():
        sig = "***" if row['composite_p'] <= 0.001 else ("**" if row['composite_p'] <= 0.01 else ("*" if row['composite_p'] <= 0.05 else ""))
        summary_lines.append(
            f"  {row['sample']:30s}  I={row['composite_I']:.4f}  p={row['composite_p']:.3f} {sig}"
            f"  ({row['n_spots']} spots, {row['n_cancer_cells']} cells, {row['mdk_nonzero_spots']} MDK+ spots)"
        )
    if len(composite_df) > 0:
        mean_comp = composite_df['composite_I'].mean()
        n_sig_comp = (composite_df['composite_p'] <= 0.001).sum()
        summary_lines.extend([
            f"  Mean composite I: {mean_comp:.4f}",
            f"  Significant (p<=0.001): {n_sig_comp}/{len(composite_df)}",
        ])

    summary_lines.extend([
        "",
        "PATIENT TRAJECTORIES (biopsy → surgery)",
        "-" * 60,
    ])
    for _, row in trajectory_df.iterrows():
        direction = "↑" if row['delta_I'] > 0 else "↓"
        summary_lines.append(
            f"  {row['patient']}: {row['biopsy_I']:.4f} → {row['surgery_I']:.4f}  "
            f"(delta = {row['delta_I']:+.4f} {direction})"
        )

    summary_lines.extend([
        "",
        "MANUSCRIPT COMPARISON",
        "-" * 60,
        f"  Current manuscript (v7.md): mean I = 0.224, {123 if n_sec_tests >= 132 else '?'}/132 sig",
        f"  v3 recomputation:           mean I = {mean_sec_I:.3f}, {n_sec_sig}/{n_sec_tests} sig",
        f"  Old persisted CSV (spot-level): mean I = 0.367 (composite, different method)",
        "",
    ])

    summary_text = "\n".join(summary_lines)
    with open(OUTPUT_DIR / "bivariate_morans_v3_summary.txt", 'w') as f:
        f.write(summary_text)
    logger.info(f"\n{summary_text}")

    # Manuscript numbers comparison
    logger.info("=" * 60)
    logger.info("NUMBERS FOR MANUSCRIPT UPDATE:")
    logger.info(f"  Per-gene: mean I = {mean_sec_I:.3f}")
    logger.info(f"  Per-gene: {n_sec_sig}/{n_sec_tests} significant (FDR < 0.05)")
    logger.info(f"  Control: mean I = {mean_ctrl_I:.3f}")
    logger.info(f"  Control: {n_ctrl_sig}/{n_ctrl_tests} significant (FDR < 0.05)")
    for _, row in trajectory_df.iterrows():
        logger.info(f"  {row['patient']}: delta = {row['delta_I']:+.4f}")
    logger.info("=" * 60)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--sample', type=int, default=None,
                        help='Single sample index (0-11)')
    args = parser.parse_args()

    if args.sample is not None:
        run_all_and_summarize(sample_indices=[args.sample])
    else:
        run_all_and_summarize()
