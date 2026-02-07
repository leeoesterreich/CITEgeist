#!/usr/bin/env python
"""
Sweep lambda_gex_reg values on region 0 to find optimal L2 regularization for GEX.
"""
import argparse
import json
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import pearsonr

REPO_ROOT = Path(__file__).parent.parent.parent.parent.parent
sys.path.insert(0, str(REPO_ROOT))

BENCHMARK_ROOT = Path(__file__).parent.parent.parent
sys.path.insert(0, str(BENCHMARK_ROOT))
from benchmark_constants import (
    ACHIEVABLE_7_CELL_PROFILE_DICT,
    ACHIEVABLE_7_MARKER_SIGNATURES,
    GT_TO_ACHIEVABLE_7_MAPPING,
)

from CITEgeist.model.citegeist_model import CitegeistModel

logging.basicConfig(level=logging.WARNING, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

# Paths
PSEUDOVISIUM_DIR = Path("/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_pseudovisium/data_protein_gt")
GT_GEX_DIR = Path("/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_pseudovisium/data_protein_gt")

CELL_TYPE_ORDER = ["B_cells", "CD4+_T_cells", "CD8+_T_cells", "Macrophages", "Endothelial", "Epithelial", "Fibroblasts"]

# GT file naming uses "pos" instead of "+" for cell type names
GT_FILENAME_MAP = {
    "B_cells": "B_cells_GT.csv",
    "CD4+_T_cells": "CD4pos_T_cells_GT.csv",
    "CD8+_T_cells": "CD8pos_T_cells_GT.csv",
    "Macrophages": "Macrophages_GT.csv",
    "Endothelial": "Endothelial_GT.csv",
    "Epithelial": "Epithelial_GT.csv",
    "Fibroblasts": "Fibroblasts_GT.csv",
}


def normalize_per_spot_cpm(df):
    """Normalize each spot (row) to CPM."""
    row_sums = df.sum(axis=1)
    row_sums = row_sums.replace(0, 1)
    return df.div(row_sums, axis=0) * 1e6


def evaluate_gex(layers_dir, region_id):
    """Evaluate GEX layers against GT for a single region."""
    sample_name = f"Xenium_region_{region_id}"

    per_ct_r = {}
    all_pred = []
    all_gt = []

    for ct in CELL_TYPE_ORDER:
        # Load prediction
        pred_path = layers_dir / f"{ct}_layer_pass1.csv"
        if not pred_path.exists():
            continue
        pred_df = pd.read_csv(pred_path, index_col=0)

        # Load GT (genes × spots format, need to transpose to spots × genes)
        gt_path = GT_GEX_DIR / "ground_truth_gex" / f"Xenium_region_{region_id}" / GT_FILENAME_MAP[ct]
        if not gt_path.exists():
            continue
        gt_df = pd.read_csv(gt_path, index_col=0).T

        # Align
        common_genes = sorted(set(pred_df.columns) & set(gt_df.columns))
        common_spots = sorted(set(pred_df.index) & set(gt_df.index))
        if len(common_genes) == 0 or len(common_spots) == 0:
            continue

        pred_aligned = pred_df.loc[common_spots, common_genes]
        gt_aligned = gt_df.loc[common_spots, common_genes]

        # CPM + log1p
        pred_log = np.log1p(normalize_per_spot_cpm(pred_aligned))
        gt_log = np.log1p(normalize_per_spot_cpm(gt_aligned))

        pred_flat = pred_log.values.flatten()
        gt_flat = gt_log.values.flatten()

        mask = np.isfinite(pred_flat) & np.isfinite(gt_flat)
        if mask.sum() > 10:
            r, _ = pearsonr(pred_flat[mask], gt_flat[mask])
            per_ct_r[ct] = r
            all_pred.append(pred_flat[mask])
            all_gt.append(gt_flat[mask])

    if all_pred:
        all_pred_cat = np.concatenate(all_pred)
        all_gt_cat = np.concatenate(all_gt)
        overall_r, _ = pearsonr(all_pred_cat, all_gt_cat)
    else:
        overall_r = float('nan')

    return overall_r, per_ct_r


def run_single_lambda(region_id, lambda_val, output_base):
    """Run GEX deconvolution with a specific lambda and evaluate."""
    sample_name = f"Xenium_region_{region_id}"
    output_dir = output_base / f"lambda_{lambda_val}"
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load data (separate GEX and CITE files)
    gex_adata = sc.read_h5ad(PSEUDOVISIUM_DIR / "h5ad_objects" / f"{sample_name}_GEX.h5ad")
    protein_adata = sc.read_h5ad(PSEUDOVISIUM_DIR / "h5ad_objects" / f"{sample_name}_CITE.h5ad")

    model = CitegeistModel(
        sample_name=sample_name,
        output_folder=str(output_dir),
        simulation=True,
        gene_expression_adata=gex_adata,
        antibody_capture_adata=protein_adata,
    )
    model.filter_gex(min_counts=25)
    model.preprocess_gex(target_sum=10000)
    model.preprocess_antibody()

    # Load profiles
    model.load_cell_profile_dict(ACHIEVABLE_7_CELL_PROFILE_DICT)

    # Load existing proportions (don't re-run)
    prop_path = Path("/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/CITEgeist/output_protein_gt") / f"{sample_name}_cell_prop_finetuned_results.csv"
    if prop_path.exists():
        prop_df = pd.read_csv(prop_path, index_col=0)
        # Align to model's spots
        common_spots = [s for s in model.gene_expression_adata.obs_names if s in prop_df.index]
        prop_aligned = prop_df.loc[common_spots]
        model.results["cell_prop"] = prop_aligned
    else:
        logger.error(f"No proportions found at {prop_path}")
        return None

    # Run GEX with specific lambda
    model.run_cell_expression_pass1(
        radius=4.0,
        alpha=0.5,
        checkpoint_interval=100,
        output_dir=str(output_dir / "checkpoints"),
        rerun=True,
        continuous_relaxation=True,
        lambda_gex_reg=lambda_val,
    )

    # Evaluate
    layers_dir = output_dir / f"{sample_name}_pass1" / "layers" / "pass1"
    overall_r, per_ct_r = evaluate_gex(layers_dir, region_id)

    return overall_r, per_ct_r


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--region-id", type=int, default=0)
    parser.add_argument("--lambdas", type=str, default="0.0,0.01,0.05,0.1,0.5,1.0")
    args = parser.parse_args()

    lambda_values = [float(x) for x in args.lambdas.split(",")]
    output_base = Path("/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/CITEgeist/output_protein_gt/lambda_sweep")
    output_base.mkdir(parents=True, exist_ok=True)

    results = {}
    for lam in lambda_values:
        logger.info(f"\n{'='*60}")
        logger.info(f"Testing lambda_gex_reg = {lam}")
        logger.info(f"{'='*60}")

        overall_r, per_ct_r = run_single_lambda(args.region_id, lam, output_base)
        results[str(lam)] = {
            "overall_r": overall_r,
            "per_ct_r": per_ct_r,
        }

        logger.info(f"lambda={lam}: overall r={overall_r:.4f}")
        for ct, r in per_ct_r.items():
            logger.info(f"  {ct}: r={r:.4f}")

    # Summary
    print("\n" + "=" * 70)
    print("LAMBDA SWEEP RESULTS")
    print("=" * 70)
    print(f"{'Lambda':<12} {'Overall r':<12} " + " ".join(f"{ct[:8]:<10}" for ct in CELL_TYPE_ORDER))
    print("-" * 70)
    for lam in lambda_values:
        r = results[str(lam)]
        ct_vals = " ".join(f"{r['per_ct_r'].get(ct, float('nan')):<10.4f}" for ct in CELL_TYPE_ORDER)
        print(f"{lam:<12.4f} {r['overall_r']:<12.4f} {ct_vals}")

    # Save
    with open(output_base / "sweep_results.json", "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved to {output_base / 'sweep_results.json'}")


if __name__ == "__main__":
    main()
