"""
Evaluate CITEgeist benchmark results against ground truth.

This module calculates metrics comparing CITEgeist predictions to
ground truth cell type proportions.

Reference for RNA-based ground truth:
    Zhao et al. (2025). "Benchmarking cell type annotation methods for 10x
    Xenium spatial transcriptomics data." BMC Bioinformatics, 26(1), 25.
    https://doi.org/10.1186/s12859-025-06044-0
"""

import argparse
import json
import logging
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd
from scipy import stats
from scipy.spatial.distance import jensenshannon

logger = logging.getLogger(__name__)


def calculate_metrics(
    gt_df: pd.DataFrame,
    pred_df: pd.DataFrame,
    cell_types: List[str],
) -> Dict[str, float]:
    """
    Calculate evaluation metrics for a single region.

    Args:
        gt_df: Ground truth proportions (spots x cell_types)
        pred_df: Predicted proportions (spots x cell_types)
        cell_types: List of cell type column names

    Returns:
        Dict with metrics per cell type and overall
    """
    metrics = {}

    # Per-cell-type metrics
    for ct in cell_types:
        if ct in pred_df.columns and ct in gt_df.columns:
            gt_vals = gt_df[ct].values
            pred_vals = pred_df[ct].values

            # Pearson correlation
            r, p = stats.pearsonr(gt_vals, pred_vals)
            metrics[f"{ct}_pearson_r"] = r
            metrics[f"{ct}_pearson_p"] = p

            # RMSE
            rmse = np.sqrt(np.mean((gt_vals - pred_vals) ** 2))
            metrics[f"{ct}_rmse"] = rmse

            # MAE
            mae = np.mean(np.abs(gt_vals - pred_vals))
            metrics[f"{ct}_mae"] = mae

    # Overall metrics (across all cell types)
    gt_matrix = gt_df[cell_types].values
    pred_matrix = pred_df[[c for c in cell_types if c in pred_df.columns]].values

    # Ensure same shape
    common_types = [c for c in cell_types if c in pred_df.columns]
    gt_matrix = gt_df[common_types].values
    pred_matrix = pred_df[common_types].values

    # Per-spot JSD (Jensen-Shannon Divergence)
    jsd_values = []
    for i in range(len(gt_df)):
        gt_row = gt_matrix[i] + 1e-10  # Add epsilon to avoid log(0)
        pred_row = pred_matrix[i] + 1e-10
        # Normalize to sum to 1
        gt_row = gt_row / gt_row.sum()
        pred_row = pred_row / pred_row.sum()
        jsd = jensenshannon(gt_row, pred_row)
        if not np.isnan(jsd):
            jsd_values.append(jsd)

    metrics["mean_jsd"] = np.mean(jsd_values) if jsd_values else np.nan
    metrics["median_jsd"] = np.median(jsd_values) if jsd_values else np.nan

    # Overall Pearson (flattened)
    gt_flat = gt_matrix.flatten()
    pred_flat = pred_matrix.flatten()
    r_overall, _ = stats.pearsonr(gt_flat, pred_flat)
    metrics["overall_pearson_r"] = r_overall

    # Overall RMSE
    metrics["overall_rmse"] = np.sqrt(np.mean((gt_flat - pred_flat) ** 2))

    # Overall MAE
    metrics["overall_mae"] = np.mean(np.abs(gt_flat - pred_flat))

    return metrics


def evaluate_region(
    region_id: int,
    gt_dir: str,
    pred_dir: str,
    prefix: str = "Xenium",
) -> Dict[str, any]:
    """
    Evaluate a single region.

    Args:
        region_id: Region identifier
        gt_dir: Directory containing ground_truth/ folder
        pred_dir: Directory containing prediction CSVs
        prefix: Filename prefix

    Returns:
        Dict with region metrics
    """
    gt_dir = Path(gt_dir)
    pred_dir = Path(pred_dir)

    # Load ground truth
    gt_path = gt_dir / "ground_truth" / f"{prefix}_region_{region_id}_prop.csv"
    if not gt_path.exists():
        raise FileNotFoundError(f"Ground truth not found: {gt_path}")
    gt_df = pd.read_csv(gt_path, index_col=0)

    # Load predictions (try finetuned first, then global)
    pred_path = pred_dir / f"{prefix}_region_{region_id}_cell_prop_finetuned_results.csv"
    if not pred_path.exists():
        pred_path = pred_dir / f"{prefix}_region_{region_id}_cell_prop_global_results.csv"
    if not pred_path.exists():
        raise FileNotFoundError(f"Predictions not found in {pred_dir}")
    pred_df = pd.read_csv(pred_path, index_col=0)

    # Align spots
    common_spots = pred_df.index.intersection(gt_df.index)
    if len(common_spots) == 0:
        raise ValueError(f"No common spots between GT and predictions for region {region_id}")

    pred_aligned = pred_df.loc[common_spots]
    gt_aligned = gt_df.loc[common_spots]

    # Get cell type columns (exclude metadata)
    metadata_cols = ["n_cells", "spot_x", "spot_y", "x", "y"]
    cell_types = [c for c in gt_aligned.columns if c not in metadata_cols]

    logger.info(f"Region {region_id}: {len(common_spots)} spots, {len(cell_types)} cell types")

    # Calculate metrics
    metrics = calculate_metrics(gt_aligned, pred_aligned, cell_types)
    metrics["region_id"] = region_id
    metrics["n_spots"] = len(common_spots)
    metrics["n_cell_types"] = len(cell_types)
    metrics["cell_types"] = cell_types

    return metrics


def evaluate_all_regions(
    gt_dir: str,
    pred_dir: str,
    n_regions: int = 5,
    output_path: Optional[str] = None,
    prefix: str = "Xenium",
) -> Dict[str, any]:
    """
    Evaluate all regions and compute summary statistics.

    Args:
        gt_dir: Directory containing ground_truth/ folder
        pred_dir: Directory containing prediction CSVs
        n_regions: Number of regions
        output_path: Optional path to save results JSON
        prefix: Filename prefix

    Returns:
        Dict with all metrics and summary
    """
    all_metrics = []

    for region_id in range(n_regions):
        try:
            metrics = evaluate_region(region_id, gt_dir, pred_dir, prefix)
            all_metrics.append(metrics)
            logger.info(f"Region {region_id}: Pearson r = {metrics['overall_pearson_r']:.4f}")
        except Exception as e:
            logger.error(f"Failed to evaluate region {region_id}: {e}")

    if not all_metrics:
        raise ValueError("No regions evaluated successfully")

    # Get cell types from first region
    cell_types = all_metrics[0]["cell_types"]

    # Compute summary statistics
    summary = {
        "n_regions": len(all_metrics),
        "total_spots": sum(m["n_spots"] for m in all_metrics),
        "cell_types": cell_types,
    }

    # Per-cell-type summary
    for ct in cell_types:
        key = f"{ct}_pearson_r"
        values = [m[key] for m in all_metrics if key in m]
        if values:
            summary[f"{ct}_mean_r"] = np.mean(values)
            summary[f"{ct}_std_r"] = np.std(values)

    # Overall summary
    overall_r_values = [m["overall_pearson_r"] for m in all_metrics]
    summary["overall_mean_pearson_r"] = np.mean(overall_r_values)
    summary["overall_std_pearson_r"] = np.std(overall_r_values)

    jsd_values = [m["mean_jsd"] for m in all_metrics if not np.isnan(m["mean_jsd"])]
    if jsd_values:
        summary["overall_mean_jsd"] = np.mean(jsd_values)

    rmse_values = [m["overall_rmse"] for m in all_metrics]
    summary["overall_mean_rmse"] = np.mean(rmse_values)

    mae_values = [m["overall_mae"] for m in all_metrics]
    summary["overall_mean_mae"] = np.mean(mae_values)

    results = {
        "summary": summary,
        "regions": all_metrics,
    }

    # Save if output path provided
    if output_path:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        # Convert numpy types for JSON serialization
        def convert_numpy(obj):
            if isinstance(obj, np.floating):
                return float(obj)
            elif isinstance(obj, np.integer):
                return int(obj)
            elif isinstance(obj, np.ndarray):
                return obj.tolist()
            elif isinstance(obj, dict):
                return {k: convert_numpy(v) for k, v in obj.items()}
            elif isinstance(obj, list):
                return [convert_numpy(v) for v in obj]
            return obj

        with open(output_path, "w") as f:
            json.dump(convert_numpy(results), f, indent=2)
        logger.info(f"Results saved to {output_path}")

    return results


def print_summary(results: Dict[str, any]) -> None:
    """Print a formatted summary of results."""
    summary = results["summary"]

    print("\n" + "=" * 60)
    print("BENCHMARK EVALUATION RESULTS")
    print("=" * 60)

    print(f"\nRegions evaluated: {summary['n_regions']}")
    print(f"Total spots: {summary['total_spots']}")
    print(f"Cell types: {', '.join(summary['cell_types'])}")

    print("\n--- Per-Cell-Type Correlations ---")
    for ct in summary["cell_types"]:
        mean_key = f"{ct}_mean_r"
        std_key = f"{ct}_std_r"
        if mean_key in summary:
            mean_r = summary[mean_key]
            std_r = summary.get(std_key, 0)
            print(f"  {ct:20s}: r = {mean_r:.4f} ± {std_r:.4f}")

    print("\n--- Overall Metrics ---")
    print(f"  Mean Pearson r: {summary['overall_mean_pearson_r']:.4f} ± {summary['overall_std_pearson_r']:.4f}")
    if "overall_mean_jsd" in summary:
        print(f"  Mean JSD:       {summary['overall_mean_jsd']:.4f}")
    print(f"  Mean RMSE:      {summary['overall_mean_rmse']:.4f}")
    print(f"  Mean MAE:       {summary['overall_mean_mae']:.4f}")

    print("\n--- Per-Region Results ---")
    for region in results["regions"]:
        print(f"  Region {region['region_id']}: r = {region['overall_pearson_r']:.4f}, "
              f"JSD = {region['mean_jsd']:.4f}, spots = {region['n_spots']}")

    print("=" * 60)


def main():
    parser = argparse.ArgumentParser(
        description="Evaluate CITEgeist benchmark results"
    )
    parser.add_argument(
        "--gt-dir",
        type=str,
        required=True,
        help="Directory containing ground_truth/ folder",
    )
    parser.add_argument(
        "--pred-dir",
        type=str,
        required=True,
        help="Directory containing prediction CSVs",
    )
    parser.add_argument(
        "--n-regions",
        type=int,
        default=5,
        help="Number of regions to evaluate",
    )
    parser.add_argument(
        "--output",
        type=str,
        default=None,
        help="Output JSON path for results",
    )
    parser.add_argument(
        "--prefix",
        type=str,
        default="Xenium",
        help="Filename prefix",
    )

    args = parser.parse_args()

    # Setup logging
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    )

    # Run evaluation
    results = evaluate_all_regions(
        gt_dir=args.gt_dir,
        pred_dir=args.pred_dir,
        n_regions=args.n_regions,
        output_path=args.output,
        prefix=args.prefix,
    )

    # Print summary
    print_summary(results)


if __name__ == "__main__":
    main()
