#!/usr/bin/env python
"""Unified proportion method comparison for simulation benchmarks.

Evaluates CITEgeist-SACE, Cell2Location, Tangram, RCTD, Seurat, and CARD
against ground truth cell type proportions across conditions and replicates.
Produces per-condition CSV summaries and a printed summary table grouped by
method.

Usage:
    python compare_prop_methods.py [--conditions high_seg mixed] [--methods ...] \
        [--reps 0 1 2 3 4] [--output_dir ...]

Output:
    benchmarks/simulation_benchmarking/evaluation/prop_comparison_{condition}.csv
    benchmarks/simulation_benchmarking/evaluation/prop_comparison_summary.csv
"""

import argparse
import os
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BENCHMARK = Path("/path/to/CITEgeist_analysis" "/benchmarks/simulation_benchmarking")
GT_BASE = Path("/path/to/CITEgeist_analysis" "/benchmarks/simulation_benchmarking/replicates")

# Ground truth proportions: GT_BASE / {condition} / ST_sim / Wu_ST_{rep}_prop.csv
GT_TEMPLATE = str(GT_BASE / "{condition}" / "ST_sim" / "Wu_ST_{rep}_prop.csv")

METHOD_PATHS = {
    "CITEgeist-SACE": str(
        BENCHMARK
        / "CITEgeist"
        / "results_sace"
        / "{condition}"
        / "Wu_rep_{rep}"
        / "citegeist"
        / "Wu_rep_{rep}_cell_prop_finetuned_results.csv"
    ),
    "Cell2Location": str(
        BENCHMARK / "Cell2Location" / "{condition}" / "cell2location_map_{rep}" / "cell2loc_deconv_predictions.csv"
    ),
    "Tangram": str(BENCHMARK / "Tangram" / "{condition}" / "tangram_map_{rep}" / "cell_type_proportions.csv"),
    "RCTD": str(BENCHMARK / "RCTD" / "{condition}" / "Wu_ST_{rep}_RCTD_deconv_predictions.csv"),
    "Seurat": str(BENCHMARK / "Seurat" / "{condition}" / "output" / "Wu_rep_{rep}_Seurat_deconv_predictions.csv"),
    "CARD": str(BENCHMARK / "CARD" / "{condition}" / "Wu_ST_{rep}_CARD_deconv_predictions.csv"),
}

ALL_METHODS = list(METHOD_PATHS.keys())
ALL_CONDITIONS = ["high_seg", "mixed"]
ALL_REPS = list(range(5))

# ---------------------------------------------------------------------------
# Import benchmarking_spot_deconv from sibling module
# ---------------------------------------------------------------------------
_SRC_DIR = Path(__file__).parent
sys.path.insert(0, str(_SRC_DIR))
from benchmarking_spot_deconv import main as spot_deconv_main  # noqa: E402

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# Columns that are NOT cell types and should be dropped before comparison.
_NON_CELLTYPE_COLS = {"spot_x", "spot_y", "recon_error"}


def _align_and_call(gt_file: str, pred_file: str) -> dict:
    """Load GT and prediction CSVs, align to shared cell type columns, then
    call benchmarking_spot_deconv.main() via temp files.

    benchmarking_spot_deconv.main() requires that both CSVs have identical
    sorted column names after dropping spot_x/spot_y. Predictions may carry
    extra columns (e.g. recon_error) that must be stripped first.

    Returns the result dict from benchmark_performance().
    """
    gt_df = pd.read_csv(gt_file, index_col=0)
    pred_df = pd.read_csv(pred_file, index_col=0)

    # Drop known non-cell-type columns
    gt_cols = set(gt_df.columns) - _NON_CELLTYPE_COLS
    pred_cols = set(pred_df.columns) - _NON_CELLTYPE_COLS

    # Use intersection of cell type columns
    shared = sorted(gt_cols & pred_cols)
    if not shared:
        raise ValueError(
            f"No shared cell type columns between GT ({sorted(gt_cols)}) " f"and prediction ({sorted(pred_cols)})"
        )

    gt_only = gt_cols - pred_cols
    pred_only = pred_cols - gt_cols
    if gt_only:
        print(f"    [WARN] GT-only columns (dropped): {sorted(gt_only)}")
    if pred_only:
        print(f"    [WARN] Pred-only columns (dropped): {sorted(pred_only)}")

    gt_df = gt_df[shared]
    pred_df = pred_df[shared]

    # Align row indices (spots)
    common_idx = gt_df.index.intersection(pred_df.index)
    if len(common_idx) == 0:
        raise ValueError("No shared spot indices between GT and prediction.")
    if len(common_idx) < len(gt_df):
        print(f"    [WARN] Only {len(common_idx)}/{len(gt_df)} GT spots matched in prediction")

    gt_df = gt_df.loc[common_idx]
    pred_df = pred_df.loc[common_idx]

    # Write to temp files and call main()
    tmp_dir = tempfile.mkdtemp(prefix="prop_bench_")
    try:
        gt_tmp = os.path.join(tmp_dir, "gt.csv")
        pred_tmp = os.path.join(tmp_dir, "pred.csv")
        gt_df.to_csv(gt_tmp)
        pred_df.to_csv(pred_tmp)
        result = spot_deconv_main(gt_tmp, pred_tmp)
    finally:
        # Clean up temp files
        for f in [gt_tmp, pred_tmp]:
            if os.path.exists(f):
                os.remove(f)
        os.rmdir(tmp_dir)

    return result


def evaluate_one(method: str, condition: str, rep: int) -> dict | None:
    """Run proportion benchmark for one method/condition/rep combination.

    Returns a flat dict of scalar metrics, or None if files are missing.
    """
    gt_path = GT_TEMPLATE.format(condition=condition, rep=rep)
    pred_path = METHOD_PATHS[method].format(condition=condition, rep=rep)

    if not os.path.isfile(gt_path):
        print(f"  [SKIP] GT not found: {gt_path}")
        return None

    if not os.path.isfile(pred_path):
        print(f"  [SKIP] Prediction not found: {pred_path}")
        return None

    try:
        result = _align_and_call(gt_path, pred_path)
    except Exception as exc:
        print(f"  [ERROR] {method}/{condition}/rep{rep}: {exc}")
        return None

    row: dict = {
        "Method": method,
        "Condition": condition,
        "Replicate": rep,
        "Pearson_r": result["corr"],
        "JSD_median": result["JSD"],
        "RMSE": result["Sum_RMSE"],
        "MAE": result["Sum_MAE"],
    }

    # Per-cell-type breakdown columns
    for cell_type, rmse_val in result["RMSE"].items():
        safe = cell_type.replace(" ", "_")
        row[f"{safe}_RMSE"] = rmse_val

    for cell_type, mae_val in result["MAE"].items():
        safe = cell_type.replace(" ", "_")
        row[f"{safe}_MAE"] = mae_val

    return row


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser(
        description="Compare cell type proportion deconvolution methods on simulated data."
    )
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
        out_csv = os.path.join(args.output_dir, f"prop_comparison_{condition}.csv")
        cond_df.to_csv(out_csv, index=False)
        print(f"\nSaved: {out_csv}")

    if not all_rows:
        print("\nNo results collected -- check paths and ensure benchmark outputs exist.")
        sys.exit(1)

    # Summary table grouped by Method (across all conditions + reps)
    summary_df = pd.DataFrame(all_rows)
    summary_cols = ["Method", "Condition", "Pearson_r", "JSD_median", "RMSE", "MAE"]
    agg = (
        summary_df[summary_cols]
        .groupby(["Method", "Condition"])
        .agg(
            Pearson_r=("Pearson_r", "mean"),
            JSD_median=("JSD_median", "mean"),
            RMSE=("RMSE", "mean"),
            MAE=("MAE", "mean"),
            N_reps=("Pearson_r", "count"),
        )
        .reset_index()
        .sort_values(["Condition", "Pearson_r"], ascending=[True, False])
    )

    print("\n" + "=" * 80)
    print("SUMMARY: Mean across replicates (sorted by Pearson_r per condition)")
    print("=" * 80)
    pd.set_option("display.float_format", "{:.4f}".format)
    pd.set_option("display.max_columns", 20)
    pd.set_option("display.width", 120)
    print(agg.to_string(index=False))

    # Save overall summary
    overall_csv = os.path.join(args.output_dir, "prop_comparison_summary.csv")
    agg.to_csv(overall_csv, index=False)
    print(f"\nOverall summary saved: {overall_csv}")


if __name__ == "__main__":
    main()
