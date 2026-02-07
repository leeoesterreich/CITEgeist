#!/usr/bin/env python
"""
Evaluate VIM-Epithelial experiment results and compare with baseline.

Compares:
  - Baseline achievable-7 (output/manual or output_achievable_7)
  - VIM-Epithelial experiment (output/vim_epithelial)

Outputs side-by-side metrics to identify improvements in Fibroblast and
Epithelial cell type correlations.

Usage:
    python evaluate_vim_epithelial.py
    python evaluate_vim_epithelial.py --baseline-dir /path/to/baseline --experiment-dir /path/to/experiment
"""

import argparse
import json
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).parent.parent.parent.parent.parent
sys.path.insert(0, str(REPO_ROOT / "Benchmarking/xenium_benchmarking/evaluation"))

from src.evaluate_benchmark import (
    evaluate_all_regions,
)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)


def compare_experiments(baseline_results: dict, experiment_results: dict) -> str:
    """Generate a comparison report between baseline and experiment."""

    lines = []
    lines.append("=" * 80)
    lines.append("VIM-EPITHELIAL EXPERIMENT COMPARISON")
    lines.append("=" * 80)
    lines.append("")
    lines.append("Change: Epithelial ['PanCK'] → ['PanCK', 'Vimentin']")
    lines.append("Hypothesis: Adding VIM to Epithelial reduces fibroblast proportion over-attribution")
    lines.append("")

    # Summary comparison
    bs = baseline_results["summary"]
    es = experiment_results["summary"]

    lines.append("-" * 80)
    lines.append(f"{'Metric':<35} {'Baseline':>15} {'VIM-Epi':>15} {'Delta':>12}")
    lines.append("-" * 80)

    # Overall metrics
    for metric_key, label in [
        ("overall_mean_pearson_r", "Overall Pearson r"),
        ("overall_mean_rmse", "Overall RMSE"),
        ("overall_mean_mae", "Overall MAE"),
        ("overall_mean_jsd", "Overall JSD"),
    ]:
        bv = bs.get(metric_key, np.nan)
        ev = es.get(metric_key, np.nan)
        delta = ev - bv if not (np.isnan(bv) or np.isnan(ev)) else np.nan
        sign = "+" if delta > 0 else ""
        lines.append(f"{label:<35} {bv:>15.4f} {ev:>15.4f} {sign}{delta:>11.4f}")

    lines.append("")
    lines.append("-" * 80)
    lines.append("Per-Cell-Type Pearson r (mean across regions)")
    lines.append("-" * 80)

    cell_types = bs.get("cell_types", [])
    for ct in cell_types:
        bv = bs.get(f"{ct}_mean_r", np.nan)
        ev = es.get(f"{ct}_mean_r", np.nan)
        delta = ev - bv if not (np.isnan(bv) or np.isnan(ev)) else np.nan
        sign = "+" if delta > 0 else ""
        marker = " ***" if ct in ["Epithelial", "Fibroblasts"] else ""
        lines.append(f"  {ct:<33} {bv:>15.4f} {ev:>15.4f} {sign}{delta:>11.4f}{marker}")

    lines.append("")
    lines.append("*** = cell types directly affected by profile change")
    lines.append("")

    # Per-region comparison for Fibroblasts and Epithelial
    lines.append("-" * 80)
    lines.append("Per-Region Detail: Fibroblasts & Epithelial Pearson r")
    lines.append("-" * 80)
    lines.append(f"{'Region':<10} {'Fib (base)':>12} {'Fib (exp)':>12} {'Epi (base)':>12} {'Epi (exp)':>12}")
    lines.append("-" * 80)

    for br, er in zip(baseline_results["regions"], experiment_results["regions"]):
        rid = br["region_id"]
        bf = br.get("Fibroblasts_pearson_r", np.nan)
        ef = er.get("Fibroblasts_pearson_r", np.nan)
        be = br.get("Epithelial_pearson_r", np.nan)
        ee = er.get("Epithelial_pearson_r", np.nan)
        lines.append(f"  {rid:<8} {bf:>12.4f} {ef:>12.4f} {be:>12.4f} {ee:>12.4f}")

    lines.append("")
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description="Evaluate VIM-Epithelial experiment")
    parser.add_argument(
        "--baseline-dir",
        type=str,
        default=str(
            REPO_ROOT / "Benchmarking/xenium_benchmarking/CITEgeist/output/manual"
        ),
    )
    parser.add_argument(
        "--experiment-dir",
        type=str,
        default=str(
            REPO_ROOT / "Benchmarking/xenium_benchmarking/CITEgeist/output/vim_epithelial"
        ),
    )
    parser.add_argument(
        "--gt-dir",
        type=str,
        default=str(REPO_ROOT / "Benchmarking/xenium_pseudovisium/data_granular_gt"),
    )
    parser.add_argument("--n-regions", type=int, default=5)
    parser.add_argument(
        "--output-dir",
        type=str,
        default=str(
            REPO_ROOT / "Benchmarking/xenium_benchmarking/evaluation/results_vim_epithelial"
        ),
    )
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Check if baseline results exist, try alternative path
    baseline_dir = Path(args.baseline_dir)
    if not baseline_dir.exists():
        alt = REPO_ROOT / "Benchmarking/xenium_benchmarking/CITEgeist/output_achievable_7"
        if alt.exists():
            baseline_dir = alt
            logger.info(f"Using alternative baseline dir: {baseline_dir}")
        else:
            logger.error(f"Baseline dir not found: {args.baseline_dir} or {alt}")
            sys.exit(1)

    experiment_dir = Path(args.experiment_dir)
    if not experiment_dir.exists():
        logger.error(f"Experiment dir not found: {experiment_dir}")
        sys.exit(1)

    # Evaluate baseline
    logger.info("Evaluating BASELINE (achievable-7)...")
    baseline_results = evaluate_all_regions(
        gt_dir=args.gt_dir,
        pred_dir=str(baseline_dir),
        n_regions=args.n_regions,
        output_path=str(output_dir / "baseline_results.json"),
    )

    # Evaluate experiment
    logger.info("Evaluating VIM-EPITHELIAL experiment...")
    experiment_results = evaluate_all_regions(
        gt_dir=args.gt_dir,
        pred_dir=str(experiment_dir),
        n_regions=args.n_regions,
        output_path=str(output_dir / "vim_epithelial_results.json"),
    )

    # Compare
    report = compare_experiments(baseline_results, experiment_results)
    print(report)

    # Save report
    report_path = output_dir / "comparison_report.txt"
    with open(report_path, "w") as f:
        f.write(report)
    logger.info(f"Comparison report saved to {report_path}")

    # Save combined results
    combined = {
        "baseline": baseline_results,
        "vim_epithelial": experiment_results,
    }
    combined_path = output_dir / "combined_results.json"

    def convert_numpy(obj):
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, dict):
            return {k: convert_numpy(v) for k, v in obj.items()}
        if isinstance(obj, list):
            return [convert_numpy(v) for v in obj]
        return obj

    with open(combined_path, "w") as f:
        json.dump(convert_numpy(combined), f, indent=2)
    logger.info(f"Combined results saved to {combined_path}")


if __name__ == "__main__":
    main()
