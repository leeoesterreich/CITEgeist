"""
Evaluate CITEgeist gene expression deconvolution results against ground truth.

This module calculates metrics comparing CITEgeist predictions to
ground truth per-cell-type gene expression.

Metrics (matching simulation benchmarking pattern):
- RMSE: Root Mean Square Error on log1p-transformed data
- NRMSE: Range-normalized RMSE
- MAE: Mean Absolute Error on log1p-transformed data

Reference for RNA-based ground truth:
    Zhao et al. (2025). "Benchmarking cell type annotation methods for 10x
    Xenium spatial transcriptomics data." BMC Bioinformatics, 26(1), 25.
    https://doi.org/10.1186/s12859-025-06044-0
"""

import argparse
import json
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from sklearn.metrics import mean_squared_error, mean_absolute_error

logger = logging.getLogger(__name__)


def load_gex_csv(filepath: Path, transpose: bool = False) -> pd.DataFrame:
    """
    Load a gene expression CSV file.

    Args:
        filepath: Path to CSV file
        transpose: If True, transpose the data (spots x genes -> genes x spots)

    Returns:
        DataFrame with genes as rows and spots as columns
    """
    df = pd.read_csv(filepath, index_col=0)
    if transpose:
        df = df.T
    return df


def find_matching_celltypes(
    gt_files: List[Path],
    pred_files: List[Path],
) -> Tuple[Dict[str, Tuple[Path, Path]], List[str], List[str]]:
    """
    Match ground truth and prediction files by cell type name.

    GT files: {CellType}_GT.csv
    Pred files: {CellType}_{markers}_layer_pass1.csv

    Returns:
        matched: Dict[celltype -> (gt_path, pred_path)]
        spurious: List of cell types with predictions but no GT
        missed: List of cell types with GT but no predictions
    """
    # Extract cell type names from GT files
    gt_celltypes = {}
    for f in gt_files:
        name = f.stem.replace("_GT", "")
        gt_celltypes[name] = f

    # Extract cell type names from prediction files
    pred_celltypes = {}
    for f in pred_files:
        # Format: {CellType}_{markers}_layer_pass1.csv
        # Try to match against known GT cell types
        stem = f.stem.replace("_layer_pass1", "")

        # Try exact match first
        matched_ct = None
        for gt_ct in gt_celltypes.keys():
            # Check if the prediction file starts with the GT cell type
            # Handle underscore variations (e.g., "B_cells" vs "B cells")
            gt_normalized = gt_ct.replace(" ", "_")
            if stem.startswith(gt_normalized):
                matched_ct = gt_ct
                break

        if matched_ct:
            pred_celltypes[matched_ct] = f
        else:
            # Try more lenient matching
            for gt_ct in gt_celltypes.keys():
                if gt_ct.replace(" ", "_") in stem or stem.startswith(gt_ct.split()[0]):
                    pred_celltypes[gt_ct] = f
                    break
            else:
                # No match found, add as spurious
                pred_celltypes[stem] = f

    # Find matches, spurious, and missed
    matched = {}
    spurious = []
    missed = []

    all_celltypes = set(gt_celltypes.keys()) | set(pred_celltypes.keys())

    for ct in all_celltypes:
        if ct in gt_celltypes and ct in pred_celltypes:
            matched[ct] = (gt_celltypes[ct], pred_celltypes[ct])
        elif ct in pred_celltypes:
            spurious.append(ct)
        else:
            missed.append(ct)

    return matched, spurious, missed


def calculate_pair_metrics(
    gt_df: pd.DataFrame,
    pred_df: pd.DataFrame,
    cell_type: str,
) -> Dict[str, float]:
    """
    Calculate metrics for a single cell type.

    Uses log1p transformation and range normalization (matching simulation benchmarking).
    """
    # Find common genes and spots
    common_genes = gt_df.index.intersection(pred_df.index)
    common_spots = gt_df.columns.intersection(pred_df.columns)

    if len(common_genes) == 0 or len(common_spots) == 0:
        logger.warning(f"No overlap for {cell_type}: genes={len(common_genes)}, spots={len(common_spots)}")
        return {
            "rmse": np.nan,
            "nrmse": np.nan,
            "mae": np.nan,
            "n_genes": 0,
            "n_spots": 0,
        }

    gt_subset = gt_df.loc[common_genes, common_spots]
    pred_subset = pred_df.loc[common_genes, common_spots]

    # Log1p transformation (matching simulation benchmarking)
    gt_log = np.log1p(gt_subset.values)
    pred_log = np.log1p(pred_subset.values)

    # Flatten for metric calculation
    gt_flat = gt_log.flatten()
    pred_flat = pred_log.flatten()

    # RMSE
    mse = mean_squared_error(gt_flat, pred_flat)
    rmse = np.sqrt(mse)

    # MAE
    mae = mean_absolute_error(gt_flat, pred_flat)

    # NRMSE (range normalization)
    range_gt = gt_log.max() - gt_log.min()
    if range_gt > 0:
        nrmse = rmse / range_gt
    else:
        nrmse = np.nan

    return {
        "rmse": rmse,
        "nrmse": nrmse,
        "mae": mae,
        "n_genes": len(common_genes),
        "n_spots": len(common_spots),
    }


def evaluate_region_gex(
    region_id: int,
    gt_dir: str,
    pred_dir: str,
    prefix: str = "Xenium",
) -> Dict[str, any]:
    """
    Evaluate gene expression deconvolution for a single region.

    Args:
        region_id: Region identifier
        gt_dir: Directory containing ground_truth_gex/Xenium_region_X/
        pred_dir: Directory containing CITEgeist output with pass1/layers/
        prefix: Filename prefix

    Returns:
        Dict with per-cell-type metrics and summary
    """
    gt_dir = Path(gt_dir)
    pred_dir = Path(pred_dir)

    # Find GT files
    gt_region_dir = gt_dir / "ground_truth_gex" / f"{prefix}_region_{region_id}"
    if not gt_region_dir.exists():
        raise FileNotFoundError(f"GT directory not found: {gt_region_dir}")

    gt_files = list(gt_region_dir.glob("*_GT.csv"))
    logger.info(f"Found {len(gt_files)} GT files in {gt_region_dir}")

    # Find prediction files - layers are in layers/pass1/ subdirectory
    pred_layers_dir = pred_dir / f"{prefix}_region_{region_id}_pass1" / "layers" / "pass1"
    if not pred_layers_dir.exists():
        raise FileNotFoundError(f"Prediction layers not found: {pred_layers_dir}")

    pred_files = list(pred_layers_dir.glob("*_layer_pass1.csv"))
    logger.info(f"Found {len(pred_files)} prediction files in {pred_layers_dir}")

    # Match cell types
    matched, spurious, missed = find_matching_celltypes(gt_files, pred_files)
    logger.info(f"Matched: {len(matched)}, Spurious: {len(spurious)}, Missed: {len(missed)}")

    metrics = {
        "region_id": region_id,
        "n_matched": len(matched),
        "n_spurious": len(spurious),
        "n_missed": len(missed),
        "spurious_celltypes": spurious,
        "missed_celltypes": missed,
        "per_celltype": {},
    }

    # Calculate metrics for matched cell types
    # Note: CITEgeist exports layers as spots×genes, so we transpose predictions to genes×spots
    for ct, (gt_path, pred_path) in matched.items():
        logger.info(f"  Evaluating {ct}...")
        gt_df = load_gex_csv(gt_path)
        pred_df = load_gex_csv(pred_path, transpose=True)  # Transpose: spots×genes -> genes×spots

        ct_metrics = calculate_pair_metrics(gt_df, pred_df, ct)
        metrics["per_celltype"][ct] = ct_metrics

        logger.info(f"    RMSE={ct_metrics['rmse']:.4f}, NRMSE={ct_metrics['nrmse']:.4f}")

    # Calculate metrics for spurious predictions (compare to zeros)
    for ct in spurious:
        logger.info(f"  [SPURIOUS] {ct}")
        # Find the prediction file for this spurious cell type
        for f in pred_files:
            if ct in f.stem:
                pred_df = load_gex_csv(f, transpose=True)
                # Create zero GT of same shape
                zero_gt = pd.DataFrame(0.0, index=pred_df.index, columns=pred_df.columns)
                ct_metrics = calculate_pair_metrics(zero_gt, pred_df, f"[SPURIOUS] {ct}")
                metrics["per_celltype"][f"[SPURIOUS] {ct}"] = ct_metrics
                break

    # Calculate metrics for missed cell types (compare zeros to GT)
    for ct in missed:
        logger.info(f"  [MISSED] {ct}")
        gt_path = gt_dir / "ground_truth_gex" / f"{prefix}_region_{region_id}" / f"{ct.replace(' ', '_')}_GT.csv"
        if gt_path.exists():
            gt_df = load_gex_csv(gt_path)
            # Create zero predictions of same shape
            zero_pred = pd.DataFrame(0.0, index=gt_df.index, columns=gt_df.columns)
            ct_metrics = calculate_pair_metrics(gt_df, zero_pred, f"[MISSED] {ct}")
            metrics["per_celltype"][f"[MISSED] {ct}"] = ct_metrics

    # Calculate holistic metrics (mean across all cell types including spurious/missed)
    all_rmse = [m["rmse"] for m in metrics["per_celltype"].values() if not np.isnan(m["rmse"])]
    all_nrmse = [m["nrmse"] for m in metrics["per_celltype"].values() if not np.isnan(m["nrmse"])]
    all_mae = [m["mae"] for m in metrics["per_celltype"].values() if not np.isnan(m["mae"])]

    metrics["holistic_rmse"] = np.mean(all_rmse) if all_rmse else np.nan
    metrics["holistic_nrmse"] = np.mean(all_nrmse) if all_nrmse else np.nan
    metrics["holistic_mae"] = np.mean(all_mae) if all_mae else np.nan

    # Calculate matched-only metrics (excluding spurious/missed)
    matched_rmse = [metrics["per_celltype"][ct]["rmse"] for ct in matched.keys()
                    if not np.isnan(metrics["per_celltype"][ct]["rmse"])]
    matched_nrmse = [metrics["per_celltype"][ct]["nrmse"] for ct in matched.keys()
                     if not np.isnan(metrics["per_celltype"][ct]["nrmse"])]
    matched_mae = [metrics["per_celltype"][ct]["mae"] for ct in matched.keys()
                   if not np.isnan(metrics["per_celltype"][ct]["mae"])]

    metrics["matched_mean_rmse"] = np.mean(matched_rmse) if matched_rmse else np.nan
    metrics["matched_mean_nrmse"] = np.mean(matched_nrmse) if matched_nrmse else np.nan
    metrics["matched_mean_mae"] = np.mean(matched_mae) if matched_mae else np.nan

    return metrics


def evaluate_all_regions_gex(
    gt_dir: str,
    pred_dir: str,
    n_regions: int = 5,
    output_path: Optional[str] = None,
    prefix: str = "Xenium",
) -> Dict[str, any]:
    """
    Evaluate gene expression deconvolution for all regions.

    Args:
        gt_dir: Directory containing ground_truth_gex/
        pred_dir: Directory containing CITEgeist outputs
        n_regions: Number of regions
        output_path: Optional path to save results JSON
        prefix: Filename prefix

    Returns:
        Dict with all metrics and summary
    """
    all_metrics = []

    for region_id in range(n_regions):
        try:
            metrics = evaluate_region_gex(region_id, gt_dir, pred_dir, prefix)
            all_metrics.append(metrics)
            logger.info(f"Region {region_id}: Holistic RMSE = {metrics['holistic_rmse']:.4f}")
        except Exception as e:
            logger.error(f"Failed to evaluate region {region_id}: {e}")

    if not all_metrics:
        raise ValueError("No regions evaluated successfully")

    # Compute summary statistics
    summary = {
        "n_regions": len(all_metrics),
        "total_matched": sum(m["n_matched"] for m in all_metrics),
        "total_spurious": sum(m["n_spurious"] for m in all_metrics),
        "total_missed": sum(m["n_missed"] for m in all_metrics),
    }

    # Holistic metrics across all regions
    holistic_rmse_values = [m["holistic_rmse"] for m in all_metrics if not np.isnan(m["holistic_rmse"])]
    holistic_nrmse_values = [m["holistic_nrmse"] for m in all_metrics if not np.isnan(m["holistic_nrmse"])]
    holistic_mae_values = [m["holistic_mae"] for m in all_metrics if not np.isnan(m["holistic_mae"])]

    summary["overall_mean_holistic_rmse"] = np.mean(holistic_rmse_values) if holistic_rmse_values else np.nan
    summary["overall_std_holistic_rmse"] = np.std(holistic_rmse_values) if holistic_rmse_values else np.nan
    summary["overall_mean_holistic_nrmse"] = np.mean(holistic_nrmse_values) if holistic_nrmse_values else np.nan
    summary["overall_mean_holistic_mae"] = np.mean(holistic_mae_values) if holistic_mae_values else np.nan

    # Matched-only metrics
    matched_rmse_values = [m["matched_mean_rmse"] for m in all_metrics if not np.isnan(m["matched_mean_rmse"])]
    matched_nrmse_values = [m["matched_mean_nrmse"] for m in all_metrics if not np.isnan(m["matched_mean_nrmse"])]
    matched_mae_values = [m["matched_mean_mae"] for m in all_metrics if not np.isnan(m["matched_mean_mae"])]

    summary["overall_mean_matched_rmse"] = np.mean(matched_rmse_values) if matched_rmse_values else np.nan
    summary["overall_mean_matched_nrmse"] = np.mean(matched_nrmse_values) if matched_nrmse_values else np.nan
    summary["overall_mean_matched_mae"] = np.mean(matched_mae_values) if matched_mae_values else np.nan

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
    """Print a formatted summary of GEX evaluation results."""
    summary = results["summary"]

    print("\n" + "=" * 60)
    print("GENE EXPRESSION DECONVOLUTION EVALUATION")
    print("=" * 60)

    print(f"\nRegions evaluated: {summary['n_regions']}")
    print(f"Total matched cell types: {summary['total_matched']}")
    print(f"Total spurious predictions: {summary['total_spurious']}")
    print(f"Total missed cell types: {summary['total_missed']}")

    print("\n--- Holistic Metrics (all cell types including spurious/missed) ---")
    print(f"  Mean RMSE:  {summary['overall_mean_holistic_rmse']:.4f} ± {summary['overall_std_holistic_rmse']:.4f}")
    print(f"  Mean NRMSE: {summary['overall_mean_holistic_nrmse']:.4f}")
    print(f"  Mean MAE:   {summary['overall_mean_holistic_mae']:.4f}")

    print("\n--- Matched-Only Metrics (excluding spurious/missed) ---")
    print(f"  Mean RMSE:  {summary['overall_mean_matched_rmse']:.4f}")
    print(f"  Mean NRMSE: {summary['overall_mean_matched_nrmse']:.4f}")
    print(f"  Mean MAE:   {summary['overall_mean_matched_mae']:.4f}")

    print("\n--- Per-Region Results ---")
    for region in results["regions"]:
        print(f"  Region {region['region_id']}: "
              f"RMSE = {region['holistic_rmse']:.4f}, "
              f"matched = {region['n_matched']}, "
              f"spurious = {region['n_spurious']}, "
              f"missed = {region['n_missed']}")

    print("=" * 60)


def main():
    parser = argparse.ArgumentParser(
        description="Evaluate CITEgeist gene expression deconvolution results"
    )
    parser.add_argument(
        "--gt-dir",
        type=str,
        required=True,
        help="Directory containing ground_truth_gex/ folder",
    )
    parser.add_argument(
        "--pred-dir",
        type=str,
        required=True,
        help="Directory containing CITEgeist output with pass1/layers/",
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
    results = evaluate_all_regions_gex(
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
