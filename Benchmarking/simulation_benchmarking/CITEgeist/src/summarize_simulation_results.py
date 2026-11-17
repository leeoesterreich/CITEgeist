#!/usr/bin/env python
"""Aggregate simulation benchmark results into summary table."""

import argparse
import json
from pathlib import Path
import pandas as pd


def collect_results(results_dir: Path) -> pd.DataFrame:
    rows = []
    for condition in ["high_seg", "mixed"]:
        for mode in ["dapi", "h_and_e", "gt_ablation"]:
            for rep in range(5):
                run_dir = results_dir / condition / mode / f"Wu_rep_{rep}"
                row = {"condition": condition, "mode": mode, "replicate": rep}
                # Try unified benchmark_results.json first, then individual files
                unified = run_dir / "benchmark_results.json"
                if unified.exists():
                    with open(unified) as f:
                        data = json.load(f)
                    # Extract proportion metrics
                    pm = data.get("proportion_metrics", {})
                    for k, v in pm.items():
                        if isinstance(v, (int, float)):
                            row[f"prop_{k}"] = v
                    # Extract SC accuracy
                    if "hungarian_sc_accuracy" in data:
                        row["sc_accuracy"] = data["hungarian_sc_accuracy"]
                    # Extract segmentation metrics
                    seg = data.get("segmentation", {})
                    for k, v in seg.items():
                        if isinstance(v, (int, float)):
                            row[f"seg_{k}"] = v
                    # Extract GEX metrics
                    gex = data.get("gex_metrics", {})
                    if isinstance(gex, dict) and "error" not in gex:
                        for k, v in gex.items():
                            if isinstance(v, (int, float)):
                                row[f"gex_{k}"] = v
                else:
                    for metric_file in [
                        "proportion_metrics.json",
                        "singlecell_metrics.json",
                        "segmentation_metrics.json",
                        "gex_metrics.json",
                    ]:
                        path = run_dir / metric_file
                        if path.exists():
                            with open(path) as f:
                                fdata = json.load(f)
                            for k, v in fdata.items():
                                if isinstance(v, (int, float)):
                                    row[f"{metric_file.replace('_metrics.json', '')}_{k}"] = v
                rows.append(row)
    return pd.DataFrame(rows)


def summarize(df: pd.DataFrame) -> pd.DataFrame:
    """Mean +/- std across replicates, grouped by condition x mode."""
    numeric_cols = df.select_dtypes(include="number").columns.drop("replicate")
    grouped = df.groupby(["condition", "mode"])[numeric_cols]
    mean = grouped.mean()
    std = grouped.std()
    summary = mean.round(4).astype(str) + " +/- " + std.round(4).astype(str)
    return summary


def main():
    parser = argparse.ArgumentParser(description="Summarize simulation benchmark results")
    parser.add_argument("--results-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args()

    df = collect_results(args.results_dir)
    summary = summarize(df)
    print(summary.to_string())
    if args.output:
        summary.to_csv(args.output)
        print(f"\nSaved to {args.output}")


if __name__ == "__main__":
    main()
