"""
Unified evaluation for all deconvolution methods on Xenium benchmark.

This script evaluates and compares:
- CITEgeist
- Cell2Location
- Tangram
- RCTD
- Seurat

For proportions: JSD, RMSE, MAE, Pearson correlation
For GEX layers: RMSE, NRMSE, MAE (Cell2Location, Tangram, CITEgeist only)

Usage:
    python evaluate_all_methods.py --gt-dir /path/to/ground_truth --output-dir /path/to/output
"""

import argparse
import json
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.spatial.distance import jensenshannon

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)

# Method configurations
METHODS = {
    "CITEgeist": {
        "output_dir": "CITEgeist/output_rna_gt",
        "pred_pattern": "{sample_name}_cell_prop_finetuned_results.csv",
        "has_gex": True,
        "gex_dir": "{sample_name}_pass1/layers/pass1",
        "gex_pattern": "underscore_pass1",  # B_cells_layer_pass1.csv
    },
    "Cell2Location": {
        "output_dir": "Cell2Location/output_rna_gt",
        "pred_pattern": "{sample_name}/{sample_name}_cell_prop_finetuned_results.csv",
        "has_gex": True,
        "gex_dir": "{sample_name}/layers",
        "gex_pattern": "space_layer",  # B cells_layer.csv
    },
    "Tangram": {
        "output_dir": "Tangram/output_rna_gt",
        "pred_pattern": "{sample_name}/{sample_name}_cell_prop_finetuned_results.csv",
        "has_gex": True,
        "gex_dir": "{sample_name}/layers",
        "gex_pattern": "space_layer",  # B cells_layer.csv
    },
    "RCTD": {
        "output_dir": "RCTD/output_rna_gt",
        "pred_pattern": "{sample_name}/{sample_name}_cell_prop_finetuned_results.csv",
        "has_gex": False,
    },
    "Seurat": {
        "output_dir": "Seurat/output_rna_gt",
        "pred_pattern": "{sample_name}/{sample_name}_cell_prop_finetuned_results.csv",
        "has_gex": False,
    },
}

# Target cell types
TARGET_CELL_TYPES = [
    "B cells",
    "T cells",
    "Macrophages",
    "Fibroblasts",
    "Epithelial",
    "Endothelial",
]


def calculate_proportion_metrics(
    gt_df: pd.DataFrame,
    pred_df: pd.DataFrame,
    cell_types: List[str],
) -> Dict[str, float]:
    """
    Calculate metrics for cell type proportions.

    Args:
        gt_df: Ground truth proportions
        pred_df: Predicted proportions
        cell_types: List of cell type columns

    Returns:
        Dict with metrics
    """
    metrics = {}

    # Align columns
    common_types = [c for c in cell_types if c in pred_df.columns and c in gt_df.columns]
    if not common_types:
        logger.warning("No common cell types found")
        return metrics

    gt_matrix = gt_df[common_types].values
    pred_matrix = pred_df[common_types].values

    # Per-cell-type metrics
    for i, ct in enumerate(common_types):
        gt_vals = gt_matrix[:, i]
        pred_vals = pred_matrix[:, i]

        # Pearson correlation
        if np.std(gt_vals) > 0 and np.std(pred_vals) > 0:
            r, _ = stats.pearsonr(gt_vals, pred_vals)
            metrics[f"{ct}_pearson_r"] = r
        else:
            metrics[f"{ct}_pearson_r"] = np.nan

        # RMSE
        metrics[f"{ct}_rmse"] = np.sqrt(np.mean((gt_vals - pred_vals) ** 2))

        # MAE
        metrics[f"{ct}_mae"] = np.mean(np.abs(gt_vals - pred_vals))

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

    metrics["mean_jsd"] = np.mean(jsd_values) if jsd_values else np.nan
    metrics["median_jsd"] = np.median(jsd_values) if jsd_values else np.nan

    # Overall metrics
    gt_flat = gt_matrix.flatten()
    pred_flat = pred_matrix.flatten()

    if np.std(gt_flat) > 0 and np.std(pred_flat) > 0:
        r_overall, _ = stats.pearsonr(gt_flat, pred_flat)
        metrics["overall_pearson_r"] = r_overall
    else:
        metrics["overall_pearson_r"] = np.nan

    metrics["overall_rmse"] = np.sqrt(np.mean((gt_flat - pred_flat) ** 2))
    metrics["overall_mae"] = np.mean(np.abs(gt_flat - pred_flat))

    return metrics


def calculate_gex_metrics(
    gt_layers_dir: Path,
    pred_layers_dir: Path,
    cell_types: List[str],
    gex_pattern: str = "space_layer",
) -> Dict[str, float]:
    """
    Calculate metrics for per-cell-type gene expression.

    Args:
        gt_layers_dir: Directory with ground truth GEX layers
        pred_layers_dir: Directory with predicted GEX layers
        cell_types: List of cell types
        gex_pattern: Naming pattern for prediction files

    Returns:
        Dict with GEX metrics
    """
    metrics = {}

    for ct in cell_types:
        # Ground truth uses underscores: B_cells_GT.csv, T_cells_GT.csv
        ct_underscore = ct.replace(" ", "_")

        # Try different GT naming patterns
        gt_path = None
        for pattern in [f"{ct_underscore}_GT.csv", f"{ct}_GT.csv", f"{ct_underscore}.csv"]:
            p = gt_layers_dir / pattern
            if p.exists():
                gt_path = p
                break

        # Prediction file patterns depend on method
        pred_path = None
        if gex_pattern == "underscore_pass1":
            # CITEgeist: B_cells_layer_pass1.csv
            patterns = [f"{ct_underscore}_layer_pass1.csv", f"{ct_underscore}_layer.csv"]
        else:
            # Cell2Location/Tangram: B cells_layer.csv (with space)
            patterns = [f"{ct}_layer.csv", f"{ct_underscore}_layer.csv"]

        for pattern in patterns:
            p = pred_layers_dir / pattern
            if p.exists():
                pred_path = p
                break

        if gt_path is None or pred_path is None:
            logger.debug(f"GEX layer not found for {ct}: gt={gt_path}, pred={pred_path}")
            continue

        try:
            gt_df = pd.read_csv(gt_path, index_col=0)
            pred_df = pd.read_csv(pred_path, index_col=0)

            # Ground truth has genes as rows, spots as columns
            # Predictions may have spots as rows, genes as columns
            # Check orientation and transpose if needed
            if gt_df.index[0].startswith('spot_') or gt_df.index[0].startswith('Spot'):
                # GT is already spots x genes
                pass
            else:
                # GT is genes x spots, transpose to spots x genes
                gt_df = gt_df.T

            if pred_df.index[0].startswith('spot_') or pred_df.index[0].startswith('Spot'):
                # Pred is already spots x genes
                pass
            else:
                # Pred is genes x spots, transpose to spots x genes
                pred_df = pred_df.T

            # Normalize column names to uppercase for case-insensitive matching
            gt_df.columns = gt_df.columns.str.upper()
            pred_df.columns = pred_df.columns.str.upper()

            # Align spots and genes
            common_spots = gt_df.index.intersection(pred_df.index)
            common_genes = gt_df.columns.intersection(pred_df.columns)

            if len(common_spots) == 0 or len(common_genes) == 0:
                logger.debug(f"No common spots/genes for {ct}: spots={len(common_spots)}, genes={len(common_genes)}")
                continue

            gt_aligned = gt_df.loc[common_spots, common_genes].values
            pred_aligned = pred_df.loc[common_spots, common_genes].values

            # Clip negative values to 0 (some methods like Tangram can produce negatives)
            gt_aligned = np.clip(gt_aligned, 0, None)
            pred_aligned = np.clip(pred_aligned, 0, None)

            # Log transform for comparison
            gt_log = np.log1p(gt_aligned)
            pred_log = np.log1p(pred_aligned)

            # RMSE
            rmse = np.sqrt(np.mean((gt_log - pred_log) ** 2))
            metrics[f"{ct}_gex_rmse"] = rmse

            # NRMSE (normalized by range)
            data_range = gt_log.max() - gt_log.min()
            if data_range > 0:
                metrics[f"{ct}_gex_nrmse"] = rmse / data_range
            else:
                metrics[f"{ct}_gex_nrmse"] = np.nan

            # MAE
            metrics[f"{ct}_gex_mae"] = np.mean(np.abs(gt_log - pred_log))

        except Exception as e:
            logger.warning(f"Failed to evaluate GEX for {ct}: {e}")

    # Average GEX metrics
    rmse_vals = [v for k, v in metrics.items() if k.endswith("_gex_rmse")]
    if rmse_vals:
        metrics["mean_gex_rmse"] = np.mean(rmse_vals)
        metrics["mean_gex_mae"] = np.mean([v for k, v in metrics.items() if k.endswith("_gex_mae")])

    return metrics


def find_prediction_file(
    base_dir: Path,
    method_config: Dict,
    sample_name: str,
) -> Optional[Path]:
    """Find the prediction file for a method."""
    pred_pattern = method_config["pred_pattern"].format(sample_name=sample_name)
    pred_path = base_dir / method_config["output_dir"] / pred_pattern

    if pred_path.exists():
        return pred_path

    # Try alternative pattern
    if "alt_pattern" in method_config:
        alt_pattern = method_config["alt_pattern"].format(sample_name=sample_name)
        alt_path = base_dir / method_config["output_dir"] / alt_pattern
        if alt_path.exists():
            return alt_path

    return None


def evaluate_method_region(
    method: str,
    region_id: int,
    base_dir: Path,
    gt_dir: Path,
    prefix: str = "Xenium",
) -> Optional[Dict]:
    """
    Evaluate a single method on a single region.

    Returns:
        Dict with metrics, or None if evaluation failed
    """
    config = METHODS[method]
    sample_name = f"{prefix}_region_{region_id}"

    # Find prediction file
    pred_path = find_prediction_file(base_dir, config, sample_name)
    if pred_path is None:
        logger.warning(f"{method} predictions not found for region {region_id}")
        return None

    # Load ground truth
    gt_path = gt_dir / "ground_truth" / f"{sample_name}_prop.csv"
    if not gt_path.exists():
        logger.warning(f"Ground truth not found: {gt_path}")
        return None

    try:
        pred_df = pd.read_csv(pred_path, index_col=0)
        gt_df = pd.read_csv(gt_path, index_col=0)

        # Align spots
        common_spots = pred_df.index.intersection(gt_df.index)
        if len(common_spots) == 0:
            logger.warning(f"No common spots for {method} region {region_id}")
            return None

        pred_aligned = pred_df.loc[common_spots]
        gt_aligned = gt_df.loc[common_spots]

        # Get cell types
        metadata_cols = ["n_cells", "spot_x", "spot_y", "x", "y", "Unknown"]
        cell_types = [c for c in gt_aligned.columns if c not in metadata_cols]

        # Calculate proportion metrics
        metrics = calculate_proportion_metrics(gt_aligned, pred_aligned, cell_types)
        metrics["region_id"] = region_id
        metrics["method"] = method
        metrics["n_spots"] = len(common_spots)
        metrics["n_cell_types"] = len(cell_types)

        # Calculate GEX metrics if available
        if config.get("has_gex", False):
            gt_gex_dir = gt_dir / "ground_truth_gex" / sample_name
            if "gex_dir" in config:
                gex_dir_pattern = config["gex_dir"].format(sample_name=sample_name)
                pred_gex_dir = base_dir / config["output_dir"] / gex_dir_pattern
            else:
                pred_gex_dir = None

            if gt_gex_dir.exists() and pred_gex_dir and pred_gex_dir.exists():
                gex_pattern = config.get("gex_pattern", "space_layer")
                gex_metrics = calculate_gex_metrics(gt_gex_dir, pred_gex_dir, cell_types, gex_pattern)
                metrics.update(gex_metrics)

        return metrics

    except Exception as e:
        logger.error(f"Failed to evaluate {method} region {region_id}: {e}")
        return None


def evaluate_all(
    base_dir: Path,
    gt_dir: Path,
    output_dir: Path,
    n_regions: int = 5,
    prefix: str = "Xenium",
) -> Dict:
    """
    Evaluate all methods on all regions.

    Returns:
        Complete evaluation results
    """
    results = {
        "methods": {},
        "comparison": {},
    }

    all_metrics = []

    for method in METHODS:
        logger.info(f"Evaluating {method}...")
        method_results = []

        for region_id in range(n_regions):
            metrics = evaluate_method_region(
                method, region_id, base_dir, gt_dir, prefix
            )
            if metrics:
                method_results.append(metrics)
                all_metrics.append(metrics)

        if method_results:
            # Compute method summary
            summary = {
                "n_regions": len(method_results),
                "total_spots": sum(m["n_spots"] for m in method_results),
            }

            # Average overall metrics
            for key in ["overall_pearson_r", "overall_rmse", "overall_mae", "mean_jsd"]:
                values = [m[key] for m in method_results if key in m and not np.isnan(m[key])]
                if values:
                    summary[f"mean_{key}"] = np.mean(values)
                    summary[f"std_{key}"] = np.std(values)

            # GEX metrics
            for key in ["mean_gex_rmse", "mean_gex_mae"]:
                values = [m[key] for m in method_results if key in m and not np.isnan(m[key])]
                if values:
                    summary[f"mean_{key}"] = np.mean(values)

            results["methods"][method] = {
                "summary": summary,
                "regions": method_results,
            }

            logger.info(f"  {method}: r={summary.get('mean_overall_pearson_r', 'N/A'):.4f}")
        else:
            logger.warning(f"  {method}: No results")

    # Create comparison table
    if all_metrics:
        comparison_df = pd.DataFrame([
            {
                "method": m["method"],
                "region": m["region_id"],
                "pearson_r": m.get("overall_pearson_r", np.nan),
                "jsd": m.get("mean_jsd", np.nan),
                "rmse": m.get("overall_rmse", np.nan),
                "mae": m.get("overall_mae", np.nan),
                "gex_rmse": m.get("mean_gex_rmse", np.nan),
            }
            for m in all_metrics
        ])

        # Save comparison table
        comparison_df.to_csv(output_dir / "comparison_table.csv", index=False)

        # Summary by method
        method_summary = comparison_df.groupby("method").agg({
            "pearson_r": ["mean", "std"],
            "jsd": ["mean", "std"],
            "rmse": ["mean", "std"],
            "mae": ["mean", "std"],
            "gex_rmse": ["mean", "std"],
        }).round(4)

        method_summary.to_csv(output_dir / "method_summary.csv")
        # Flatten multi-level columns for JSON serialization
        method_summary_flat = method_summary.copy()
        method_summary_flat.columns = ['_'.join(col).strip() for col in method_summary_flat.columns.values]
        results["comparison"]["method_summary"] = method_summary_flat.to_dict()

    return results


def create_comparison_figures(
    results: Dict,
    output_dir: Path,
):
    """Create comparison figures for all methods."""
    if "methods" not in results or not results["methods"]:
        logger.warning("No methods to compare")
        return

    methods = list(results["methods"].keys())
    metrics_to_plot = ["mean_overall_pearson_r", "mean_mean_jsd", "mean_overall_rmse"]
    labels = ["Pearson r", "JSD", "RMSE"]

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    for ax, metric, label in zip(axes, metrics_to_plot, labels):
        values = []
        errors = []
        for method in methods:
            summary = results["methods"][method]["summary"]
            values.append(summary.get(metric, 0))
            std_key = metric.replace("mean_", "std_")
            errors.append(summary.get(std_key, 0))

        x = np.arange(len(methods))
        ax.bar(x, values, yerr=errors, capsize=5, color=plt.cm.tab10.colors[:len(methods)])
        ax.set_xticks(x)
        ax.set_xticklabels(methods, rotation=45, ha="right")
        ax.set_ylabel(label)
        ax.set_title(label)

    plt.tight_layout()
    plt.savefig(output_dir / "method_comparison.png", dpi=150, bbox_inches="tight")
    plt.close()

    logger.info(f"Saved comparison figure to {output_dir / 'method_comparison.png'}")


def main():
    parser = argparse.ArgumentParser(
        description="Evaluate all deconvolution methods on Xenium benchmark"
    )
    parser.add_argument(
        "--base-dir",
        type=str,
        default=".",
        help="Base directory containing method output folders",
    )
    parser.add_argument(
        "--gt-dir",
        type=str,
        required=True,
        help="Directory containing ground_truth/ folder",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="evaluation/results",
        help="Output directory for results",
    )
    parser.add_argument(
        "--n-regions",
        type=int,
        default=5,
        help="Number of regions",
    )
    parser.add_argument(
        "--prefix",
        type=str,
        default="Xenium",
        help="Filename prefix",
    )

    args = parser.parse_args()

    base_dir = Path(args.base_dir)
    gt_dir = Path(args.gt_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("=" * 60)
    logger.info("Unified Xenium Benchmark Evaluation")
    logger.info("=" * 60)
    logger.info(f"Base directory: {base_dir}")
    logger.info(f"Ground truth: {gt_dir}")
    logger.info(f"Output: {output_dir}")

    # Run evaluation
    results = evaluate_all(
        base_dir=base_dir,
        gt_dir=gt_dir,
        output_dir=output_dir,
        n_regions=args.n_regions,
        prefix=args.prefix,
    )

    # Save results
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

    with open(output_dir / "full_results.json", "w") as f:
        json.dump(convert_numpy(results), f, indent=2)

    # Create figures
    create_comparison_figures(results, output_dir)

    # Print summary
    logger.info("\n" + "=" * 60)
    logger.info("EVALUATION SUMMARY")
    logger.info("=" * 60)

    for method, data in results["methods"].items():
        summary = data["summary"]
        r = summary.get("mean_overall_pearson_r", "N/A")
        jsd = summary.get("mean_mean_jsd", "N/A")
        rmse = summary.get("mean_overall_rmse", "N/A")
        gex = summary.get("mean_mean_gex_rmse", "N/A")

        r_str = f"{r:.4f}" if isinstance(r, float) else r
        jsd_str = f"{jsd:.4f}" if isinstance(jsd, float) else jsd
        rmse_str = f"{rmse:.4f}" if isinstance(rmse, float) else rmse
        gex_str = f"{gex:.4f}" if isinstance(gex, float) else gex

        logger.info(f"{method:15s}: r={r_str}, JSD={jsd_str}, RMSE={rmse_str}, GEX_RMSE={gex_str}")

    logger.info("=" * 60)
    logger.info(f"Results saved to {output_dir}")


if __name__ == "__main__":
    main()
