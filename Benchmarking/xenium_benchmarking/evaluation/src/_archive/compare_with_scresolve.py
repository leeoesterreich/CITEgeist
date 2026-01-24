#!/usr/bin/env python
"""
Compare all methods including scResolve on region 0.

Since scResolve only completed for region 0, this script does a
targeted comparison for that single region.
"""

import json
import logging
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from scipy.spatial.distance import jensenshannon

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

BASE_DIR = Path("/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking")
GT_DIR = Path("/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_pseudovisium/data_granular_gt")

# Method configurations for region 0
METHODS = {
    "CITEgeist": {
        "pred_path": BASE_DIR / "CITEgeist/output_granular/Xenium_region_0/Xenium_region_0_deconv_predictions.csv",
    },
    "Cell2Location": {
        "pred_path": BASE_DIR / "Cell2Location/output_granular/Xenium_region_0/Xenium_region_0_cell_prop_finetuned_results.csv",
    },
    "Tangram": {
        "pred_path": BASE_DIR / "Tangram/output_granular/Xenium_region_0/Xenium_region_0_cell_prop_finetuned_results.csv",
        "normalize": True,  # Tangram outputs need normalization
    },
    "RCTD": {
        "pred_path": BASE_DIR / "RCTD/output_granular/Xenium_region_0/Xenium_region_0_cell_prop_finetuned_results.csv",
    },
    "Seurat": {
        "pred_path": BASE_DIR / "Seurat/output_granular/Xenium_region_0/Xenium_region_0_cell_prop_finetuned_results.csv",
        "fix_colnames": True,  # Seurat uses R-style dot notation
    },
    "scResolve": {
        "pred_path": BASE_DIR / "scResolve/output_morphology/Xenium_region_0/Xenium_region_0_deconv_predictions.csv",
    },
}

# R-style name mapping for Seurat
SEURAT_NAME_MAPPING = {
    "B.cells": "B cells",
    "CD8..T.cells": "CD8+ T cells",
    "CD8+.T.cells": "CD8+ T cells",
    "Macrophages": "Macrophages",
    "Mixed.Immune": "Mixed Immune",
    "Myofibroblasts": "Myofibroblasts",
    "Stromal": "Stromal",
    "Epithelial": "Epithelial",
    "Endothelial": "Endothelial",
    "Proliferating.T": "Proliferating T",
    "Vascular.Stromal": "Vascular Stromal",
}


def normalize_celltype_names(df):
    """Normalize Seurat R-style column names."""
    new_columns = {}
    for col in df.columns:
        if col in SEURAT_NAME_MAPPING:
            new_columns[col] = SEURAT_NAME_MAPPING[col]
        else:
            normalized = col.replace(".", " ").replace("  ", "+ ")
            normalized = normalized.replace("CD8  T cells", "CD8+ T cells")
            new_columns[col] = normalized
    return df.rename(columns=new_columns)


def normalize_proportions(df):
    """Normalize rows to sum to 1."""
    row_sums = df.sum(axis=1).replace(0, 1)
    return df.div(row_sums, axis=0)


def calculate_metrics(gt_df, pred_df, cell_types):
    """Calculate evaluation metrics."""
    metrics = {}

    common_types = [c for c in cell_types if c in pred_df.columns and c in gt_df.columns]
    if not common_types:
        return {"error": "No common cell types"}

    gt_matrix = gt_df[common_types].values
    pred_matrix = pred_df[common_types].values

    # Per-spot JSD
    jsd_values = []
    for i in range(len(gt_df)):
        gt_row = gt_matrix[i] + 1e-10
        pred_row = pred_matrix[i] + 1e-10
        gt_row = gt_row / gt_row.sum()
        pred_row = pred_row / pred_row.sum()
        jsd = jensenshannon(gt_row, pred_row)
        if not np.isnan(jsd):
            jsd_values.append(jsd)

    metrics["jsd"] = np.mean(jsd_values) if jsd_values else np.nan

    # Overall metrics (flattened)
    gt_flat = gt_matrix.flatten()
    pred_flat = pred_matrix.flatten()

    if np.std(gt_flat) > 0 and np.std(pred_flat) > 0:
        r, _ = stats.pearsonr(gt_flat, pred_flat)
        metrics["pearson_r"] = r
    else:
        metrics["pearson_r"] = np.nan

    metrics["rmse"] = np.sqrt(np.mean((gt_flat - pred_flat) ** 2))
    metrics["mae"] = np.mean(np.abs(gt_flat - pred_flat))

    return metrics


def main():
    logger.info("=" * 60)
    logger.info("Comparing All Methods Including scResolve - Region 0")
    logger.info("=" * 60)

    # Load ground truth
    gt_path = GT_DIR / "ground_truth" / "Xenium_region_0_prop.csv"
    gt_df = pd.read_csv(gt_path, index_col=0)

    # Get cell types from GT
    metadata_cols = ["n_cells", "spot_x", "spot_y", "x", "y", "Unknown"]
    cell_types = [c for c in gt_df.columns if c not in metadata_cols]
    logger.info(f"Cell types: {cell_types}")

    results = {}

    for method, config in METHODS.items():
        logger.info(f"\nEvaluating {method}...")

        pred_path = config["pred_path"]
        if not pred_path.exists():
            logger.warning(f"  Predictions not found: {pred_path}")
            continue

        try:
            pred_df = pd.read_csv(pred_path, index_col=0)

            # Apply method-specific preprocessing
            if config.get("fix_colnames"):
                pred_df = normalize_celltype_names(pred_df)

            if config.get("normalize"):
                pred_df = normalize_proportions(pred_df)

            # Align spots
            common_spots = pred_df.index.intersection(gt_df.index)
            if len(common_spots) == 0:
                logger.warning(f"  No common spots found")
                continue

            pred_aligned = pred_df.loc[common_spots]
            gt_aligned = gt_df.loc[common_spots]

            metrics = calculate_metrics(gt_aligned, pred_aligned, cell_types)
            metrics["n_spots"] = len(common_spots)
            metrics["n_cell_types"] = len([c for c in cell_types if c in pred_df.columns])

            results[method] = metrics

            logger.info(f"  Pearson r: {metrics.get('pearson_r', 'N/A'):.4f}")
            logger.info(f"  JSD: {metrics.get('jsd', 'N/A'):.4f}")
            logger.info(f"  RMSE: {metrics.get('rmse', 'N/A'):.4f}")
            logger.info(f"  MAE: {metrics.get('mae', 'N/A'):.4f}")

        except Exception as e:
            logger.error(f"  Error: {e}")

    # Create comparison table
    logger.info("\n" + "=" * 60)
    logger.info("COMPARISON RESULTS (Region 0)")
    logger.info("=" * 60)

    comparison = []
    for method, metrics in results.items():
        comparison.append({
            "Method": method,
            "Pearson r": f"{metrics.get('pearson_r', np.nan):.4f}",
            "JSD": f"{metrics.get('jsd', np.nan):.4f}",
            "RMSE": f"{metrics.get('rmse', np.nan):.4f}",
            "MAE": f"{metrics.get('mae', np.nan):.4f}",
        })

    df = pd.DataFrame(comparison)

    # Sort by Pearson r descending
    df['pearson_sort'] = df['Pearson r'].astype(float)
    df = df.sort_values('pearson_sort', ascending=False).drop('pearson_sort', axis=1)

    print("\n")
    print(df.to_string(index=False))

    # Save results
    output_dir = BASE_DIR / "evaluation/results_with_scresolve"
    output_dir.mkdir(parents=True, exist_ok=True)

    df.to_csv(output_dir / "region_0_comparison.csv", index=False)

    with open(output_dir / "region_0_metrics.json", "w") as f:
        # Convert numpy types
        def convert(obj):
            if isinstance(obj, (np.floating, np.float64)):
                return float(obj)
            elif isinstance(obj, (np.integer, np.int64)):
                return int(obj)
            return obj

        results_converted = {k: {kk: convert(vv) for kk, vv in v.items()} for k, v in results.items()}
        json.dump(results_converted, f, indent=2)

    logger.info(f"\nResults saved to {output_dir}")


if __name__ == "__main__":
    main()
