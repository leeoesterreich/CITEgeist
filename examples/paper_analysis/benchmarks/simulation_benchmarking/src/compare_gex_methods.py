#!/usr/bin/env python
"""Unified GEX method comparison for simulation benchmarks.

Evaluates Cell2Location, Tangram, and CITEgeist-SACE against
ground truth GEX layers across conditions and replicates. Produces per-condition
CSV summaries and a printed summary table grouped by method.

Usage:
    python compare_gex_methods.py [--conditions high_seg mixed] [--methods ...] \
        [--reps 0 1 2 3 4] [--output_dir ...]

Output:
    benchmarks/simulation_benchmarking/evaluation/gex_comparison_{condition}.csv
"""

import argparse
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BENCHMARK = Path("/path/to/CITEgeist_analysis" "/benchmarks/simulation_benchmarking")
GT_BASE = Path("/path/to/CITEgeist_analysis/replicates")

# Ground truth layers: GT_BASE / {condition} / ST_GEX_sim / sample_{rep} / layers
GT_TEMPLATE = str(GT_BASE / "{condition}" / "ST_GEX_sim" / "sample_{rep}" / "layers")

METHOD_PATHS = {
    "Cell2Location": str(BENCHMARK / "Cell2Location" / "{condition}" / "cell2location_map_{rep}" / "layers"),
    "Tangram": str(BENCHMARK / "Tangram" / "{condition}" / "tangram_map_{rep}" / "layers"),
    "CITEgeist-SACE": str(BENCHMARK / "CITEgeist" / "results_sace" / "{condition}" / "Wu_rep_{rep}" / "layers"),
}

ALL_METHODS = list(METHOD_PATHS.keys())
ALL_CONDITIONS = ["high_seg", "mixed"]
ALL_REPS = list(range(5))

# ---------------------------------------------------------------------------
# Import calculate_rmse from sibling module
# ---------------------------------------------------------------------------
_SRC_DIR = Path(__file__).parent
sys.path.insert(0, str(_SRC_DIR))
from benchmarking_gex import calculate_rmse  # noqa: E402


def evaluate_one(method: str, condition: str, rep: int) -> dict | None:
    """Run calculate_rmse for one method/condition/rep combination.

    Returns a flat dict of scalar metrics, or None if the layers directory
    does not exist / has no evaluable files.
    """
    gt_path = GT_TEMPLATE.format(condition=condition, rep=rep)
    layers_path = METHOD_PATHS[method].format(condition=condition, rep=rep)

    if not os.path.isdir(gt_path):
        print(f"  [SKIP] GT not found: {gt_path}")
        return None

    if not os.path.isdir(layers_path):
        print(f"  [SKIP] layers not found: {layers_path}")
        return None

    try:
        result = calculate_rmse(gt_path, layers_path)
    except Exception as exc:
        print(f"  [ERROR] {method}/{condition}/rep{rep}: {exc}")
        return None

    if result is None:
        return None

    row: dict = {
        "Method": method,
        "Condition": condition,
        "Replicate": rep,
        "Mean_Pearson_r": result["mean_pearson_r"],
        "Mean_Spot_Layer_r": result["mean_spot_layer_r"],
        "Avg_RMSE": result["average_rmse"],
        "Avg_NRMSE": result["average_nrmse"],
        "Avg_MAE": result["average_mae"],
    }

    # Per-cell-type breakdown columns
    for cell_type, metrics in result["metrics_per_cell_type"].items():
        safe = cell_type.replace(" ", "_")
        row[f"{safe}_pearson_r"] = metrics["pearson_r"]
        row[f"{safe}_spot_layer_r"] = metrics["spot_layer_r"]
        row[f"{safe}_RMSE"] = metrics["RMSE"]
        row[f"{safe}_NRMSE"] = metrics["NRMSE"]
        row[f"{safe}_MAE"] = metrics["MAE"]

    return row


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser(description="Compare GEX deconvolution methods on simulated data.")
    parser.add_argument(
        "--conditions",
        nargs="+",
        default=ALL_CONDITIONS,
        choices=ALL_CONDITIONS,
        help="Conditions to evaluate (default: all)",
    )
    parser.add_argument(
        "--methods",
        nargs="+",
        default=ALL_METHODS,
        choices=ALL_METHODS,
        help="Methods to evaluate (default: all)",
    )
    parser.add_argument(
        "--reps",
        nargs="+",
        type=int,
        default=ALL_REPS,
        help="Replicate indices to evaluate (default: 0-4)",
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        default=str(BENCHMARK / "evaluation"),
        help="Directory for output CSVs",
    )
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    all_rows: list[dict] = []

    for condition in args.conditions:
        print(f"\n{'='*60}")
        print(f"Condition: {condition}")
        print(f"{'='*60}")

        cond_rows: list[dict] = []

        for method in args.methods:
            for rep in args.reps:
                print(f"\n  Method={method}  rep={rep}")
                row = evaluate_one(method, condition, rep)
                if row is not None:
                    cond_rows.append(row)
                    all_rows.append(row)

        if not cond_rows:
            print(f"  No results for condition={condition}; skipping CSV.")
            continue

        cond_df = pd.DataFrame(cond_rows)
        out_csv = os.path.join(args.output_dir, f"gex_comparison_{condition}.csv")
        cond_df.to_csv(out_csv, index=False)
        print(f"\nSaved: {out_csv}")

    if not all_rows:
        print("\nNo results collected — check paths and run SACE benchmarks first.")
        sys.exit(1)

    # Summary table grouped by Method (across all conditions + reps)
    summary_df = pd.DataFrame(all_rows)
    summary_cols = ["Method", "Condition", "Mean_Pearson_r", "Mean_Spot_Layer_r", "Avg_RMSE", "Avg_NRMSE", "Avg_MAE"]
    agg = (
        summary_df[summary_cols]
        .groupby(["Method", "Condition"])
        .agg(
            Mean_Pearson_r=("Mean_Pearson_r", "mean"),
            Mean_Spot_Layer_r=("Mean_Spot_Layer_r", "mean"),
            Avg_RMSE=("Avg_RMSE", "mean"),
            Avg_NRMSE=("Avg_NRMSE", "mean"),
            Avg_MAE=("Avg_MAE", "mean"),
            N_reps=("Mean_Pearson_r", "count"),
        )
        .reset_index()
        .sort_values(["Condition", "Mean_Spot_Layer_r"], ascending=[True, False])
    )

    print("\n" + "=" * 80)
    print("SUMMARY: Mean across replicates (sorted by Mean_Spot_Layer_r per condition)")
    print("=" * 80)
    pd.set_option("display.float_format", "{:.4f}".format)
    pd.set_option("display.max_columns", 20)
    pd.set_option("display.width", 120)
    print(agg.to_string(index=False))

    # Save overall summary
    overall_csv = os.path.join(args.output_dir, "gex_comparison_summary.csv")
    agg.to_csv(overall_csv, index=False)
    print(f"\nOverall summary saved: {overall_csv}")


if __name__ == "__main__":
    main()
