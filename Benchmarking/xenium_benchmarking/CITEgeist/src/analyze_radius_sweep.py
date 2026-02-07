#!/usr/bin/env python
"""
Analyze radius sweep results and generate comparison metrics.

Aggregates results across all radii and regions, computes evaluation metrics
against ground truth, and generates visualization.

Usage:
    python analyze_radius_sweep.py
    python analyze_radius_sweep.py --output-dir /path/to/sweep/results
"""

import argparse
import json
import logging
import sys
from pathlib import Path
from typing import Any, Dict, List, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial.distance import jensenshannon
from scipy.stats import pearsonr

# Add paths
REPO_ROOT = Path(__file__).parent.parent.parent.parent.parent
sys.path.insert(0, str(REPO_ROOT))

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)

# Radius configuration
RADII = [50, 105, 205, 305]
RING_NAMES = {50: "0 rings", 105: "1 ring", 205: "2 rings", 305: "3 rings"}
N_REGIONS = 5


def load_ground_truth(region_id: int, gt_dir: Path) -> pd.DataFrame:
    """Load ground truth proportions for a region."""
    gt_path = gt_dir / f"Xenium_region_{region_id}_ground_truth.csv"
    if not gt_path.exists():
        # Try alternative path
        gt_path = gt_dir / "ground_truth" / f"Xenium_region_{region_id}_ground_truth.csv"

    if not gt_path.exists():
        logger.warning(f"Ground truth not found: {gt_path}")
        return None

    return pd.read_csv(gt_path, index_col=0)


def load_predictions(region_id: int, radius: int, sweep_dir: Path) -> pd.DataFrame:
    """Load predictions for a specific region and radius."""
    pred_path = (
        sweep_dir
        / f"radius_{radius}"
        / f"Xenium_region_{region_id}"
        / f"Xenium_region_{region_id}_deconv_predictions.csv"
    )

    if not pred_path.exists():
        logger.warning(f"Predictions not found: {pred_path}")
        return None

    return pd.read_csv(pred_path, index_col=0)


def compute_metrics(
    pred: pd.DataFrame, gt: pd.DataFrame
) -> Dict[str, float]:
    """
    Compute evaluation metrics between predictions and ground truth.

    Returns:
        Dict with JSD, Pearson, RMSE metrics
    """
    # Align columns and indices
    common_cols = list(set(pred.columns) & set(gt.columns))
    common_idx = list(set(pred.index) & set(gt.index))

    if not common_cols or not common_idx:
        return {"jsd": np.nan, "pearson": np.nan, "rmse": np.nan}

    pred_aligned = pred.loc[common_idx, common_cols].values
    gt_aligned = gt.loc[common_idx, common_cols].values

    # Normalize rows to sum to 1 (for JSD)
    pred_norm = pred_aligned / (pred_aligned.sum(axis=1, keepdims=True) + 1e-10)
    gt_norm = gt_aligned / (gt_aligned.sum(axis=1, keepdims=True) + 1e-10)

    # Per-spot JSD, then average
    jsd_values = []
    for i in range(len(common_idx)):
        jsd = jensenshannon(pred_norm[i], gt_norm[i])
        if not np.isnan(jsd):
            jsd_values.append(jsd)
    jsd_mean = np.mean(jsd_values) if jsd_values else np.nan

    # Flatten for Pearson and RMSE
    pred_flat = pred_aligned.flatten()
    gt_flat = gt_aligned.flatten()

    # Pearson correlation
    if len(pred_flat) > 1:
        pearson, _ = pearsonr(pred_flat, gt_flat)
    else:
        pearson = np.nan

    # RMSE
    rmse = np.sqrt(np.mean((pred_flat - gt_flat) ** 2))

    return {
        "jsd": jsd_mean,
        "pearson": pearson,
        "rmse": rmse,
    }


def analyze_sweep(
    sweep_dir: Path,
    gt_dir: Path,
) -> pd.DataFrame:
    """
    Analyze all sweep results and compute metrics.

    Returns:
        DataFrame with metrics for each region/radius combination
    """
    results = []

    for region_id in range(N_REGIONS):
        gt = load_ground_truth(region_id, gt_dir)
        if gt is None:
            logger.warning(f"Skipping region {region_id}: no ground truth")
            continue

        for radius in RADII:
            pred = load_predictions(region_id, radius, sweep_dir)
            if pred is None:
                logger.warning(f"Skipping region {region_id}, radius {radius}: no predictions")
                continue

            metrics = compute_metrics(pred, gt)

            results.append({
                "region": region_id,
                "radius": radius,
                "rings": RING_NAMES[radius],
                "jsd": metrics["jsd"],
                "pearson": metrics["pearson"],
                "rmse": metrics["rmse"],
            })

    return pd.DataFrame(results)


def plot_sweep_results(df: pd.DataFrame, output_path: Path):
    """Generate visualization of sweep results."""
    fig, axes = plt.subplots(1, 3, figsize=(12, 4))

    metrics = [("jsd", "JSD (lower is better)"),
               ("pearson", "Pearson r (higher is better)"),
               ("rmse", "RMSE (lower is better)")]

    for ax, (metric, title) in zip(axes, metrics):
        # Group by radius and compute mean/std
        grouped = df.groupby("radius")[metric].agg(["mean", "std"])

        x = range(len(RADII))
        ax.errorbar(
            x,
            grouped["mean"],
            yerr=grouped["std"],
            marker="o",
            capsize=5,
            linewidth=2,
            markersize=8,
        )

        ax.set_xticks(x)
        ax.set_xticklabels([RING_NAMES[r] for r in RADII])
        ax.set_xlabel("Neighborhood Size")
        ax.set_ylabel(metric.upper())
        ax.set_title(title)
        ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()

    logger.info(f"Plot saved to {output_path}")


def find_optimal_radius(df: pd.DataFrame) -> Dict[str, Any]:
    """Identify optimal radius for each metric."""
    grouped = df.groupby("radius").agg({
        "jsd": "mean",
        "pearson": "mean",
        "rmse": "mean",
    })

    optimal = {
        "best_jsd": {
            "radius": int(grouped["jsd"].idxmin()),
            "value": grouped["jsd"].min(),
        },
        "best_pearson": {
            "radius": int(grouped["pearson"].idxmax()),
            "value": grouped["pearson"].max(),
        },
        "best_rmse": {
            "radius": int(grouped["rmse"].idxmin()),
            "value": grouped["rmse"].min(),
        },
    }

    # Check if there's a clear winner
    best_radii = [
        optimal["best_jsd"]["radius"],
        optimal["best_pearson"]["radius"],
        optimal["best_rmse"]["radius"],
    ]

    if len(set(best_radii)) == 1:
        optimal["recommendation"] = f"radius={best_radii[0]} ({RING_NAMES[best_radii[0]]}) is optimal across all metrics"
    else:
        optimal["recommendation"] = "Trade-offs exist between metrics; see detailed results"

    return optimal


def main():
    parser = argparse.ArgumentParser(description="Analyze radius sweep results")

    parser.add_argument(
        "--sweep-dir",
        type=str,
        default=str(REPO_ROOT / "Benchmarking/xenium_benchmarking/CITEgeist/output/radius_sweep"),
        help="Directory containing sweep results",
    )
    parser.add_argument(
        "--gt-dir",
        type=str,
        default=str(REPO_ROOT / "Benchmarking/xenium_pseudovisium/data_protein_gt"),
        help="Directory containing ground truth",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default=None,
        help="Output directory for analysis results (default: sweep-dir)",
    )

    args = parser.parse_args()

    sweep_dir = Path(args.sweep_dir)
    gt_dir = Path(args.gt_dir)
    output_dir = Path(args.output_dir) if args.output_dir else sweep_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("Analyzing radius sweep results...")
    logger.info(f"Sweep dir: {sweep_dir}")
    logger.info(f"GT dir: {gt_dir}")

    # Analyze results
    df = analyze_sweep(sweep_dir, gt_dir)

    if df.empty:
        logger.error("No results found to analyze")
        sys.exit(1)

    # Save detailed results
    csv_path = output_dir / "radius_sweep_metrics.csv"
    df.to_csv(csv_path, index=False)
    logger.info(f"Detailed metrics saved to {csv_path}")

    # Generate plot
    plot_path = output_dir / "radius_sweep_comparison.png"
    plot_sweep_results(df, plot_path)

    # Find optimal
    optimal = find_optimal_radius(df)

    # Save summary
    summary = {
        "metrics_by_radius": df.groupby("radius").agg({
            "jsd": ["mean", "std"],
            "pearson": ["mean", "std"],
            "rmse": ["mean", "std"],
        }).to_dict(),
        "optimal": optimal,
    }

    summary_path = output_dir / "radius_sweep_summary.json"
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2, default=str)
    logger.info(f"Summary saved to {summary_path}")

    # Print results
    print("\n" + "=" * 60)
    print("RADIUS SWEEP RESULTS")
    print("=" * 60)

    print("\nMean metrics by radius:")
    print(df.groupby("radius")[["jsd", "pearson", "rmse"]].mean().round(4).to_string())

    print("\nOptimal radii:")
    print(f"  Best JSD:     radius={optimal['best_jsd']['radius']} ({RING_NAMES[optimal['best_jsd']['radius']]}), value={optimal['best_jsd']['value']:.4f}")
    print(f"  Best Pearson: radius={optimal['best_pearson']['radius']} ({RING_NAMES[optimal['best_pearson']['radius']]}), value={optimal['best_pearson']['value']:.4f}")
    print(f"  Best RMSE:    radius={optimal['best_rmse']['radius']} ({RING_NAMES[optimal['best_rmse']['radius']]}), value={optimal['best_rmse']['value']:.4f}")
    print(f"\nRecommendation: {optimal['recommendation']}")


if __name__ == "__main__":
    main()
