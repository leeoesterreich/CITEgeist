#!/usr/bin/env python3
"""
Aggregate cellularity ablation results into a summary table.

Reads result_region{0-4}_{baseline,stardist_dapi,oracle_gt}.json files
and produces a comparison table + ablation_summary.csv.

Usage:
    python evaluate_cellularity_ablation.py --results-dir output/cellularity_ablation
"""
import argparse
import json
import sys
from pathlib import Path

import pandas as pd

CONDITIONS = ["baseline", "stardist_dapi", "oracle_gt", "oracle_abundance_prior", "oracle_rowsum_bounds"]
REGIONS = range(5)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--results-dir", type=str, required=True)
    args = parser.parse_args()

    results_dir = Path(args.results_dir)
    rows = []

    for region in REGIONS:
        for cond in CONDITIONS:
            f = results_dir / f"result_region{region}_{cond}.json"
            if not f.exists():
                print(f"MISSING: {f.name}", file=sys.stderr)
                continue
            with open(f) as fh:
                d = json.load(fh)
            row = {
                "region": region,
                "condition": cond,
                "flat_r": d["flat_r"],
                "flat_rmse": d["flat_rmse"],
                "n_spots": d["n_spots"],
            }
            for t, m in d.get("per_type", {}).items():
                row[f"r_{t}"] = m["r"]
                row[f"rmse_{t}"] = m["rmse"]
            rows.append(row)

    df = pd.DataFrame(rows)
    if df.empty:
        print("No results found!", file=sys.stderr)
        return

    out_csv = results_dir / "ablation_summary.csv"
    df.to_csv(out_csv, index=False)
    print(f"Saved: {out_csv}")

    print("\n=== SUMMARY (mean across 5 regions) ===\n")
    summary = df.groupby("condition")[["flat_r", "flat_rmse"]].mean()
    summary = summary.reindex(CONDITIONS)
    print(summary.to_string())

    if "baseline" in summary.index:
        baseline_r = summary.loc["baseline", "flat_r"]
        print(f"\n=== DELTA from baseline (flat_r={baseline_r:.4f}) ===\n")
        for cond in CONDITIONS:
            if cond == "baseline" or cond not in summary.index:
                continue
            delta = summary.loc[cond, "flat_r"] - baseline_r
            print(f"  {cond}: delta_r = {delta:+.4f}")

    r_cols = [c for c in df.columns if c.startswith("r_")]
    if r_cols:
        print("\n=== PER-TYPE r (mean across regions) ===\n")
        type_summary = df.groupby("condition")[r_cols].mean().reindex(CONDITIONS)
        print(type_summary.to_string())


if __name__ == "__main__":
    main()
