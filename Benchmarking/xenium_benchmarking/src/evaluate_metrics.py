"""
Evaluate benchmark metrics for CITEgeist predictions.

This module computes JSD, RMSE, MAE, and Pearson correlation
between CITEgeist predictions and ground truth proportions.
"""

import os
import sys
import argparse
import logging
import json
from pathlib import Path
from typing import Dict, Any, List, Optional

import numpy as np
import pandas as pd
from scipy.spatial.distance import jensenshannon
from scipy.stats import pearsonr
from sklearn.metrics import mean_absolute_error, mean_squared_error

# Add simulation_benchmarking/src to path for reusing existing code
BENCH_SRC = Path(__file__).parent.parent.parent / "simulation_benchmarking" / "src"
sys.path.insert(0, str(BENCH_SRC))

logger = logging.getLogger(__name__)


def benchmark_proportions(
    ground_truth: np.ndarray,
    predictions: np.ndarray,
    column_names: List[str],
) -> Dict[str, Any]:
    """
    Calculate benchmark metrics for cell type proportions.

    Args:
        ground_truth: (n_spots, n_cell_types) ground truth proportions
        predictions: (n_spots, n_cell_types) predicted proportions
        column_names: List of cell type names

    Returns:
        Dict with JSD, RMSE, MAE, and correlation metrics
    """
    n_spots, n_types = ground_truth.shape

    # Jensen-Shannon Divergence per spot
    jsd_values = np.zeros(n_spots)
    for i in range(n_spots):
        gt = ground_truth[i, :]
        pred = predictions[i, :]

        # Ensure non-negative values
        gt = np.maximum(gt, 0)
        pred = np.maximum(pred, 0)

        # Check if both distributions are valid (sum > 0)
        gt_sum = np.sum(gt)
        pred_sum = np.sum(pred)

        if gt_sum > 0 and pred_sum > 0:
            # Normalize to probability distributions
            gt_norm = gt / gt_sum
            pred_norm = pred / pred_sum
            jsd_values[i] = jensenshannon(gt_norm, pred_norm, base=2)
        elif gt_sum == 0 and pred_sum == 0:
            jsd_values[i] = 0.0  # Both empty = perfect match
        else:
            jsd_values[i] = 1.0  # Max divergence if one is empty

    # RMSE and MAE per cell type
    rmse_per_type = {}
    mae_per_type = {}
    total_mse = 0

    for i, ct in enumerate(column_names):
        gt = ground_truth[:, i]
        pred = predictions[:, i]

        mse = np.sum((gt - pred) ** 2)
        total_mse += mse
        rmse_per_type[ct] = np.sqrt(mse / n_spots)
        mae_per_type[ct] = mean_absolute_error(gt, pred)

    # Global metrics
    global_rmse = np.sqrt(total_mse / (n_spots * n_types))
    global_mae = np.mean(list(mae_per_type.values()))

    # Pearson correlation
    flat_gt = ground_truth.flatten()
    flat_pred = predictions.flatten()
    correlation, p_value = pearsonr(flat_gt, flat_pred)

    return {
        "JSD_median": np.median(jsd_values),
        "JSD_mean": np.mean(jsd_values),
        "JSD_q25": np.percentile(jsd_values, 25),
        "JSD_q75": np.percentile(jsd_values, 75),
        "RMSE_per_type": rmse_per_type,
        "RMSE_global": global_rmse,
        "MAE_per_type": mae_per_type,
        "MAE_global": global_mae,
        "Pearson_r": correlation,
        "Pearson_p": p_value,
        "n_spots": n_spots,
        "n_cell_types": n_types,
    }


def evaluate_region(
    ground_truth_path: str,
    predictions_path: str,
    output_path: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Evaluate predictions against ground truth for a single region.

    Args:
        ground_truth_path: Path to ground truth CSV
        predictions_path: Path to predictions CSV
        output_path: Optional path to save metrics

    Returns:
        Dict with benchmark metrics
    """
    # Load data
    gt_df = pd.read_csv(ground_truth_path, index_col=0)
    pred_df = pd.read_csv(predictions_path, index_col=0)

    logger.info(f"Ground truth: {gt_df.shape}")
    logger.info(f"Predictions: {pred_df.shape}")

    # Remove metadata columns
    metadata_cols = ["n_cells", "spot_x", "spot_y"]
    gt_df = gt_df.drop(columns=[c for c in metadata_cols if c in gt_df.columns])
    pred_df = pred_df.drop(columns=[c for c in metadata_cols if c in pred_df.columns])

    # Sort by spot names and cell types
    gt_df = gt_df.sort_index().sort_index(axis=1)
    pred_df = pred_df.sort_index().sort_index(axis=1)

    # Find common spots and cell types
    common_spots = gt_df.index.intersection(pred_df.index)
    common_types = gt_df.columns.intersection(pred_df.columns)

    logger.info(f"Common spots: {len(common_spots)}")
    logger.info(f"Common cell types: {list(common_types)}")

    if len(common_spots) == 0:
        raise ValueError("No common spots between ground truth and predictions!")

    if len(common_types) == 0:
        raise ValueError("No common cell types between ground truth and predictions!")

    # Subset to common spots and types
    gt_df = gt_df.loc[common_spots, common_types]
    pred_df = pred_df.loc[common_spots, common_types]

    # Convert to numpy
    ground_truth = gt_df.values
    predictions = pred_df.values
    column_names = list(gt_df.columns)

    # Calculate metrics
    metrics = benchmark_proportions(ground_truth, predictions, column_names)

    # Add file info
    metrics["ground_truth_path"] = str(ground_truth_path)
    metrics["predictions_path"] = str(predictions_path)

    # Save if path provided
    if output_path:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        # Save summary metrics as JSON
        summary = {k: v for k, v in metrics.items()
                   if not isinstance(v, dict)}
        with open(output_path.with_suffix(".json"), "w") as f:
            json.dump(summary, f, indent=2)

        # Save per-cell-type metrics as CSV
        per_type_df = pd.DataFrame({
            "cell_type": column_names,
            "RMSE": [metrics["RMSE_per_type"][ct] for ct in column_names],
            "MAE": [metrics["MAE_per_type"][ct] for ct in column_names],
        })
        per_type_df.to_csv(output_path.with_suffix(".csv"), index=False)

        logger.info(f"Saved metrics to {output_path}")

    return metrics


def evaluate_all_regions(
    input_dir: str,
    output_dir: str,
    predictions_dir: str,
    n_regions: int = 5,
    prefix: str = "Xenium",
) -> Dict[str, Any]:
    """
    Evaluate all regions and aggregate metrics.

    Args:
        input_dir: Directory containing ground_truth/
        output_dir: Directory to save metrics
        predictions_dir: Directory containing predictions
        n_regions: Number of regions
        prefix: Filename prefix

    Returns:
        Dict with aggregated metrics
    """
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)
    predictions_dir = Path(predictions_dir)

    output_dir.mkdir(parents=True, exist_ok=True)

    all_metrics = []

    for region_id in range(n_regions):
        logger.info(f"\n{'='*60}")
        logger.info(f"Evaluating region {region_id}")
        logger.info(f"{'='*60}")

        gt_path = input_dir / "ground_truth" / f"{prefix}_region_{region_id}_prop.csv"
        pred_path = predictions_dir / f"{prefix}_region_{region_id}" / f"{prefix}_region_{region_id}_deconv_predictions.csv"

        # Alternative path if predictions are in a different structure
        if not pred_path.exists():
            pred_path = predictions_dir / f"{prefix}_region_{region_id}" / "cell_proportions.csv"

        if not pred_path.exists():
            logger.warning(f"Predictions not found for region {region_id}")
            continue

        try:
            metrics = evaluate_region(
                ground_truth_path=str(gt_path),
                predictions_path=str(pred_path),
                output_path=str(output_dir / f"metrics_region_{region_id}"),
            )
            metrics["region_id"] = region_id
            all_metrics.append(metrics)
        except Exception as e:
            logger.error(f"Error evaluating region {region_id}: {e}")
            continue

    if not all_metrics:
        raise ValueError("No regions were successfully evaluated!")

    # Aggregate metrics
    aggregated = aggregate_metrics(all_metrics)

    # Save aggregated metrics
    with open(output_dir / "aggregated_metrics.json", "w") as f:
        json.dump(aggregated, f, indent=2)

    # Save summary table
    summary_df = pd.DataFrame([
        {
            "region_id": m["region_id"],
            "n_spots": m["n_spots"],
            "JSD_median": m["JSD_median"],
            "RMSE_global": m["RMSE_global"],
            "MAE_global": m["MAE_global"],
            "Pearson_r": m["Pearson_r"],
        }
        for m in all_metrics
    ])
    summary_df.to_csv(output_dir / "summary_table.csv", index=False)

    logger.info(f"\n{'='*60}")
    logger.info("Evaluation Summary")
    logger.info(f"{'='*60}")
    logger.info(f"Regions evaluated: {len(all_metrics)}")
    logger.info(f"Mean JSD: {aggregated['JSD_median_mean']:.4f} ± {aggregated['JSD_median_std']:.4f}")
    logger.info(f"Mean RMSE: {aggregated['RMSE_global_mean']:.4f} ± {aggregated['RMSE_global_std']:.4f}")
    logger.info(f"Mean MAE: {aggregated['MAE_global_mean']:.4f} ± {aggregated['MAE_global_std']:.4f}")
    logger.info(f"Mean Pearson r: {aggregated['Pearson_r_mean']:.4f} ± {aggregated['Pearson_r_std']:.4f}")

    return aggregated


def aggregate_metrics(metrics_list: List[Dict]) -> Dict[str, Any]:
    """
    Aggregate metrics across regions.

    Args:
        metrics_list: List of metrics dicts from each region

    Returns:
        Dict with mean and std of each metric
    """
    # Scalar metrics to aggregate
    scalar_keys = ["JSD_median", "JSD_mean", "RMSE_global", "MAE_global", "Pearson_r"]

    aggregated = {}
    for key in scalar_keys:
        values = [m[key] for m in metrics_list if key in m]
        if values:
            aggregated[f"{key}_mean"] = np.mean(values)
            aggregated[f"{key}_std"] = np.std(values)
            aggregated[f"{key}_min"] = np.min(values)
            aggregated[f"{key}_max"] = np.max(values)

    aggregated["n_regions"] = len(metrics_list)
    aggregated["total_spots"] = sum(m["n_spots"] for m in metrics_list)

    return aggregated


def benchmark_gex(
    ground_truth_dir: str,
    predictions_parquet: str,
    normalize: str = 'range',
) -> Dict[str, Any]:
    """
    Benchmark gene expression deconvolution.

    Compares CITEgeist deconvolved GEX against ground truth per-cell-type
    expression.

    Args:
        ground_truth_dir: Directory with {CellType}_GT.csv files
        predictions_parquet: Path to CITEgeist gene_expression_pass1.parquet
        normalize: Normalization method ('range' or 'mean')

    Returns:
        Dict with GEX benchmark metrics per cell type and overall
    """
    from pathlib import Path

    ground_truth_dir = Path(ground_truth_dir)
    predictions_parquet = Path(predictions_parquet)

    if not predictions_parquet.exists():
        logger.warning(f"Predictions file not found: {predictions_parquet}")
        return {"error": "predictions_not_found"}

    # Load CITEgeist predictions
    # Format: index = "spot_X:::CellType", columns = genes
    pred_df = pd.read_parquet(predictions_parquet)

    # Parse the index to get spot and cell type
    spot_celltype = pred_df.index.str.split(":::", expand=True)
    pred_df["spot"] = spot_celltype.get_level_values(0)
    pred_df["cell_type"] = spot_celltype.get_level_values(1)

    metrics_per_cell_type = {}

    # Get list of ground truth files
    gt_files = list(ground_truth_dir.glob("*_GT.csv"))
    if not gt_files:
        logger.warning(f"No ground truth files found in {ground_truth_dir}")
        return {"error": "no_ground_truth_files"}

    # Cell types to exclude from aggregate metrics (garbage bin)
    exclude_from_aggregate = ["Unknown", "Unassigned"]

    for gt_file in gt_files:
        # Extract cell type from filename
        ct_safe = gt_file.stem.replace("_GT", "")
        # Convert back to original format
        ct_original = ct_safe.replace("_", " ").replace("pos", "+")

        # Also try variations
        ct_variants = [
            ct_original,
            ct_safe,
            ct_original.replace("CD4+ T", "CD4+ T"),  # Keep as-is
            ct_original.replace("CD8+ T", "CD8+ T"),  # Keep as-is
        ]

        # Load ground truth (genes x spots)
        gt_df = pd.read_csv(gt_file, index_col=0)

        # Find matching cell type in predictions
        matching_ct = None
        for ct in ct_variants:
            if ct in pred_df["cell_type"].values:
                matching_ct = ct
                break

        if matching_ct is None:
            logger.debug(f"Cell type {ct_original} not found in predictions")
            continue

        # Get predictions for this cell type
        ct_pred = pred_df[pred_df["cell_type"] == matching_ct].copy()
        ct_pred = ct_pred.set_index("spot")
        ct_pred = ct_pred.drop(columns=["cell_type"])

        # Transpose to (genes x spots) to match ground truth
        ct_pred_t = ct_pred.T

        # Find common genes and spots
        common_genes = gt_df.index.intersection(ct_pred_t.index)
        common_spots = gt_df.columns.intersection(ct_pred_t.columns)

        if len(common_genes) == 0 or len(common_spots) == 0:
            logger.debug(f"No common genes/spots for {ct_original}")
            continue

        # Subset to common
        gt_aligned = gt_df.loc[common_genes, common_spots]
        pred_aligned = ct_pred_t.loc[common_genes, common_spots]

        # Log-transform (add 1 to avoid log(0))
        gt_log = np.log1p(gt_aligned.values)
        pred_log = np.log1p(pred_aligned.values)

        # Calculate metrics
        mse = mean_squared_error(gt_log, pred_log)
        rmse = np.sqrt(mse)
        mae = mean_absolute_error(gt_log, pred_log)

        # Normalized RMSE
        if normalize == 'range':
            range_gt = gt_log.max() - gt_log.min()
            nrmse = rmse / range_gt if range_gt != 0 else np.nan
        elif normalize == 'mean':
            mean_gt = gt_log.mean()
            nrmse = rmse / mean_gt if mean_gt != 0 else np.nan
        else:
            nrmse = np.nan

        # Pearson correlation (flatten to 1D)
        if gt_log.std() > 0 and pred_log.std() > 0:
            corr, _ = pearsonr(gt_log.flatten(), pred_log.flatten())
        else:
            corr = np.nan

        metrics_per_cell_type[ct_original] = {
            'RMSE': rmse,
            'NRMSE': nrmse,
            'MAE': mae,
            'Pearson_r': corr,
            'n_genes': len(common_genes),
            'n_spots': len(common_spots),
        }

        logger.info(
            f"  {ct_original}: RMSE={rmse:.4f}, NRMSE={nrmse:.4f}, "
            f"MAE={mae:.4f}, r={corr:.4f}"
        )

    if not metrics_per_cell_type:
        return {"error": "no_matching_cell_types"}

    # Aggregate metrics - EXCLUDE Unknown/Unassigned from summary (garbage bin)
    known_types = {k: v for k, v in metrics_per_cell_type.items()
                   if k not in exclude_from_aggregate}

    all_rmse = [m['RMSE'] for m in known_types.values()]
    all_nrmse = [m['NRMSE'] for m in known_types.values() if not np.isnan(m['NRMSE'])]
    all_mae = [m['MAE'] for m in known_types.values()]
    all_corr = [m['Pearson_r'] for m in known_types.values() if not np.isnan(m['Pearson_r'])]

    return {
        "per_cell_type": metrics_per_cell_type,  # Include all for reference
        "RMSE_mean": np.mean(all_rmse) if all_rmse else np.nan,
        "RMSE_median": np.median(all_rmse) if all_rmse else np.nan,
        "NRMSE_mean": np.nanmean(all_nrmse) if all_nrmse else np.nan,
        "NRMSE_median": np.nanmedian(all_nrmse) if all_nrmse else np.nan,
        "MAE_mean": np.mean(all_mae) if all_mae else np.nan,
        "MAE_median": np.median(all_mae) if all_mae else np.nan,
        "Pearson_r_mean": np.mean(all_corr) if all_corr else np.nan,
        "Pearson_r_median": np.median(all_corr) if all_corr else np.nan,
        "n_cell_types": len(known_types),  # Count of known types only
        "n_cell_types_total": len(metrics_per_cell_type),  # Total including Unknown
    }


def print_metrics_summary(metrics: Dict[str, Any]) -> None:
    """Pretty-print metrics summary."""
    print("\n" + "=" * 60)
    print("BENCHMARK METRICS SUMMARY")
    print("=" * 60)

    print(f"\nSpots: {metrics['n_spots']}")
    print(f"Cell types: {metrics['n_cell_types']}")

    print(f"\nJensen-Shannon Divergence:")
    print(f"  Median: {metrics['JSD_median']:.4f}")
    print(f"  Mean:   {metrics['JSD_mean']:.4f}")
    print(f"  IQR:    [{metrics['JSD_q25']:.4f}, {metrics['JSD_q75']:.4f}]")

    print(f"\nGlobal Metrics:")
    print(f"  RMSE:      {metrics['RMSE_global']:.4f}")
    print(f"  MAE:       {metrics['MAE_global']:.4f}")
    print(f"  Pearson r: {metrics['Pearson_r']:.4f}")

    print(f"\nPer-Cell-Type RMSE:")
    for ct, rmse in sorted(metrics["RMSE_per_type"].items()):
        print(f"  {ct}: {rmse:.4f}")

    print("=" * 60)


def main():
    parser = argparse.ArgumentParser(
        description="Evaluate CITEgeist predictions against ground truth"
    )
    parser.add_argument(
        "--ground-truth",
        type=str,
        required=True,
        help="Path to ground truth CSV",
    )
    parser.add_argument(
        "--predictions",
        type=str,
        required=True,
        help="Path to predictions CSV",
    )
    parser.add_argument(
        "--output",
        type=str,
        default=None,
        help="Path to save metrics (optional)",
    )

    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO)

    metrics = evaluate_region(
        ground_truth_path=args.ground_truth,
        predictions_path=args.predictions,
        output_path=args.output,
    )

    print_metrics_summary(metrics)


if __name__ == "__main__":
    main()
