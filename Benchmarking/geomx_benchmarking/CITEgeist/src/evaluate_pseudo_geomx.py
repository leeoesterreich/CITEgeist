"""
evaluate_pseudo_geomx.py

Evaluation script for pseudo-GeoMx CITEgeist predictions.
Compares predicted proportions against ground truth AOI proportions.

Input
-----
  --pred-dir   : Directory containing pipeline outputs (e.g., predicted_proportions.csv)
  --gt-dir     : Directory containing ground truth (e.g., aoi_gt_proportions.csv)
  --output     : Path to write evaluation results (JSON)

Output
------
  Evaluation metrics (JSON):
    - per_type_metrics: {type_name: {r, rmse}}
    - overall_metrics: {r, rmse}
    - jsd_metrics: {mean, median}
    - mean_r: mean Pearson r across all types
"""

import argparse
import json
import logging
import sys
import warnings
from pathlib import Path
from typing import Dict, Tuple

import numpy as np
import pandas as pd
from scipy.spatial.distance import jensenshannon
from scipy.stats import pearsonr

# ── logging ──────────────────────────────────────────────────────────────────
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(name)s %(levelname)s %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)


def load_predictions(pred_dir: Path) -> pd.DataFrame:
    """Load predicted proportions from pipeline output directory.

    Parameters
    ----------
    pred_dir : Path
        Directory containing predicted_proportions.csv

    Returns
    -------
    pd.DataFrame
        Shape (n_aois, n_types), index = AOI names, columns = type names
    """
    pred_path = pred_dir / "predicted_proportions.csv"
    if not pred_path.exists():
        raise FileNotFoundError(f"Predictions not found: {pred_path}")

    df = pd.read_csv(pred_path, index_col=0)
    logger.info(f"Loaded predictions: {df.shape[0]} AOIs × {df.shape[1]} types")
    return df


def load_ground_truth(gt_dir: Path) -> pd.DataFrame:
    """Load ground truth proportions from evaluation directory.

    Parameters
    ----------
    gt_dir : Path
        Directory containing aoi_gt_proportions.csv

    Returns
    -------
    pd.DataFrame
        Shape (n_aois, n_types), index = AOI names, columns = type names
    """
    gt_path = gt_dir / "aoi_gt_proportions.csv"
    if not gt_path.exists():
        raise FileNotFoundError(f"Ground truth not found: {gt_path}")

    df = pd.read_csv(gt_path, index_col=0)
    logger.info(f"Loaded ground truth: {df.shape[0]} AOIs × {df.shape[1]} types")
    return df


def align_proportions(pred: pd.DataFrame, gt: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Align predictions and ground truth to common AOIs and types.

    If column names don't match (e.g., M1/M2 discovers Type_0 vs GT "Epithelial"),
    returns empty intersection and triggers fallback to flattened correlation.

    Parameters
    ----------
    pred : pd.DataFrame
        Predicted proportions (n_aois × n_types_pred)
    gt : pd.DataFrame
        Ground truth proportions (n_aois × n_types_gt)

    Returns
    -------
    Tuple[pd.DataFrame, pd.DataFrame]
        Aligned (pred_aligned, gt_aligned) for common AOIs and types.
        If no common types, returns full datasets for flattened fallback.
    """
    # Align to common AOIs
    common_aois = pred.index.intersection(gt.index)
    if len(common_aois) == 0:
        raise ValueError("No common AOIs between predictions and ground truth")

    pred_aoi = pred.loc[common_aois]
    gt_aoi = gt.loc[common_aois]

    # Align to common types
    common_types = pred_aoi.columns.intersection(gt_aoi.columns)

    if len(common_types) == 0:
        logger.warning(
            f"No common type names found. "
            f"Predicted types: {list(pred_aoi.columns)}, "
            f"GT types: {list(gt_aoi.columns)}. "
            f"Will fall back to flattened overall correlation."
        )
        return pred_aoi, gt_aoi

    pred_aligned = pred_aoi[common_types]
    gt_aligned = gt_aoi[common_types]

    logger.info(f"Aligned to {len(common_aois)} AOIs and {len(common_types)} types")
    return pred_aligned, gt_aligned


def evaluate_proportions(pred_path: Path, gt_path: Path) -> Dict:
    """Compare predicted and ground truth proportions.

    Parameters
    ----------
    pred_path : Path
        Directory containing predicted_proportions.csv
    gt_path : Path
        Directory containing aoi_gt_proportions.csv

    Returns
    -------
    Dict
        Evaluation results:
        - per_type_metrics: {type_name: {r, rmse}}
        - overall_metrics: {r, rmse}
        - jsd_metrics: {mean, median}
        - mean_r: mean Pearson r across types
    """
    # Load data
    pred = load_predictions(pred_path)
    gt = load_ground_truth(gt_path)

    # Align
    pred_aligned, gt_aligned = align_proportions(pred, gt)

    results = {}

    # ────────────────────────────────────────────────────────────────────────
    # Per-type Pearson r and RMSE
    # ────────────────────────────────────────────────────────────────────────
    per_type_metrics = {}

    types_match = set(pred_aligned.columns) == set(gt_aligned.columns)
    common_cols = list(set(pred_aligned.columns) & set(gt_aligned.columns))

    if common_cols:
        for col in common_cols:
            pred_col = pred_aligned[col].values
            gt_col = gt_aligned[col].values

            r, _ = pearsonr(pred_col, gt_col)
            rmse = np.sqrt(np.mean((pred_col - gt_col) ** 2))

            per_type_metrics[col] = {"r": float(r), "rmse": float(rmse)}

        results["per_type_metrics"] = per_type_metrics
        logger.info(f"Computed per-type metrics for {len(per_type_metrics)} types")
    else:
        logger.warning("No matching type names — per-type metrics unavailable. " "Using flattened correlation only.")
        results["type_name_mismatch"] = {
            "predicted": list(pred_aligned.columns),
            "ground_truth": list(gt_aligned.columns),
        }

    # ────────────────────────────────────────────────────────────────────────
    # Overall flattened Pearson r and RMSE
    # ────────────────────────────────────────────────────────────────────────
    pred_flat = pred_aligned.values.flatten()
    gt_flat = gt_aligned.values.flatten()

    overall_r, _ = pearsonr(pred_flat, gt_flat)
    overall_rmse = np.sqrt(np.mean((pred_flat - gt_flat) ** 2))

    results["overall_metrics"] = {
        "r": float(overall_r),
        "rmse": float(overall_rmse),
    }
    logger.info(f"Overall: r={overall_r:.4f}, rmse={overall_rmse:.4f}")

    # ────────────────────────────────────────────────────────────────────────
    # Per-AOI Jensen-Shannon divergence
    # ────────────────────────────────────────────────────────────────────────
    jsd_values = []
    for i in range(len(pred_aligned)):
        pred_row = pred_aligned.iloc[i].values
        gt_row = gt_aligned.iloc[i].values

        # Ensure valid probability distributions (sum to ~1)
        pred_row = np.clip(pred_row, 0, np.inf)
        gt_row = np.clip(gt_row, 0, np.inf)

        pred_sum = pred_row.sum()
        gt_sum = gt_row.sum()

        if pred_sum > 0:
            pred_row = pred_row / pred_sum
        else:
            pred_row = np.ones_like(pred_row) / len(pred_row)

        if gt_sum > 0:
            gt_row = gt_row / gt_sum
        else:
            gt_row = np.ones_like(gt_row) / len(gt_row)

        jsd = jensenshannon(pred_row, gt_row)
        jsd_values.append(jsd)

    jsd_mean = np.mean(jsd_values)
    jsd_median = np.median(jsd_values)

    results["jsd_metrics"] = {
        "mean": float(jsd_mean),
        "median": float(jsd_median),
    }
    logger.info(f"JSD: mean={jsd_mean:.4f}, median={jsd_median:.4f}")

    # ────────────────────────────────────────────────────────────────────────
    # Mean per-type Pearson r
    # ────────────────────────────────────────────────────────────────────────
    if len(per_type_metrics) > 0:
        mean_r = np.mean([m["r"] for m in per_type_metrics.values()])
        results["mean_r"] = float(mean_r)
        logger.info(f"Mean per-type r: {mean_r:.4f}")
    else:
        logger.warning("Mean per-type r skipped: no common types found")

    return results


def main():
    """CLI entry point."""
    parser = argparse.ArgumentParser(description="Evaluate pseudo-GeoMx CITEgeist predictions")
    parser.add_argument(
        "--pred-dir",
        type=Path,
        required=True,
        help="Directory with pipeline outputs (contains predicted_proportions.csv)",
    )
    parser.add_argument(
        "--gt-dir",
        type=Path,
        required=True,
        help="Directory with ground truth (contains aoi_gt_proportions.csv)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        required=True,
        help="Path to write evaluation results (JSON)",
    )

    args = parser.parse_args()

    # Validate inputs
    if not args.pred_dir.is_dir():
        logger.error(f"Predictions directory not found: {args.pred_dir}")
        sys.exit(1)

    if not args.gt_dir.is_dir():
        logger.error(f"Ground truth directory not found: {args.gt_dir}")
        sys.exit(1)

    # Evaluate
    try:
        results = evaluate_proportions(args.pred_dir, args.gt_dir)
    except Exception as e:
        logger.error(f"Evaluation failed: {e}", exc_info=True)
        sys.exit(1)

    # Write output
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w") as f:
        json.dump(results, f, indent=2)

    logger.info(f"Results written to {args.output}")


if __name__ == "__main__":
    main()
