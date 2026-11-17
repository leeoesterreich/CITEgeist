#!/usr/bin/env python
"""Aggregate SACE ablation results across 8 variants × 2 conditions × 5 reps.

Reads per-rep sace_metrics.json from results_sace_{variant}/{condition}/Wu_rep_{id}/layers/
and produces sace_ablation_summary.csv in the evaluation directory.

Usage:
    python aggregate_sace_ablation.py [--output_dir ...]
"""

import sys as _sys
from pathlib import Path as _Path

_sys.path.insert(0, str(_Path(__file__).resolve().parents[2]))
import _paths as _cg  # noqa: E402

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

BENCHMARK = Path(str(_cg.REPO_ROOT / "Benchmarking" / "simulation_benchmarking"))

VARIANTS = [
    "marker_i1_cont", "marker_i1_lr", "marker_i5_cont", "marker_i5_lr",
    "prop_i1_cont", "prop_i1_lr", "prop_i5_cont", "prop_i5_lr",
]
CONDITIONS = ["high_seg", "mixed"]
REPS = list(range(5))


def load_metrics(variant: str, condition: str, rep: int) -> "dict | None":
    metrics_path = (
        BENCHMARK / "CITEgeist" / f"results_sace_{variant}"
        / condition / f"Wu_rep_{rep}" / "layers" / "sace_metrics.json"
    )
    if not metrics_path.exists():
        print(f"MISSING: {metrics_path}", file=sys.stderr)
        return None
    with open(metrics_path) as f:
        return json.load(f)


def main():
    parser = argparse.ArgumentParser(description="Aggregate SACE ablation results")
    parser.add_argument(
        "--output_dir",
        type=Path,
        default=BENCHMARK / "evaluation",
        help="Output directory for summary CSV",
    )
    args = parser.parse_args()

    rows = []
    for variant in VARIANTS:
        parts = variant.split("_")
        init_method = parts[0]
        max_iter = int(parts[1][1:])
        integer_alloc = parts[2] == "lr"

        for condition in CONDITIONS:
            rep_metrics = []
            for rep in REPS:
                m = load_metrics(variant, condition, rep)
                if m is not None:
                    rep_metrics.append(m)

            if not rep_metrics:
                print(f"No metrics for {variant}/{condition}", file=sys.stderr)
                continue

            rows.append({
                "variant": variant,
                "init_method": init_method,
                "max_iter": max_iter,
                "integer_alloc": integer_alloc,
                "condition": condition,
                "mean_pearson_r": pd.Series([m["mean_pearson_r"] for m in rep_metrics]).mean(),
                "mean_spot_layer_r": pd.Series([m["mean_spot_layer_r"] for m in rep_metrics]).mean(),
                "avg_rmse": pd.Series([m["average_rmse"] for m in rep_metrics]).mean(),
                "avg_nrmse": pd.Series([m.get("average_nrmse", float("nan")) for m in rep_metrics]).mean(),
                "avg_mae": pd.Series([m.get("average_mae", float("nan")) for m in rep_metrics]).mean(),
                "n_reps": len(rep_metrics),
            })

    df = pd.DataFrame(rows)
    out_path = args.output_dir / "sace_ablation_summary.csv"
    args.output_dir.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, index=False)
    print(f"\nSaved: {out_path} ({len(df)} rows)")

    print("\n=== SACE Ablation Summary ===\n")
    for condition in CONDITIONS:
        sub = df[df["condition"] == condition].sort_values("mean_pearson_r", ascending=False)
        print(f"--- {condition} ---")
        print(sub[["variant", "mean_pearson_r", "mean_spot_layer_r", "avg_rmse", "n_reps"]].to_string(index=False))
        print()


if __name__ == "__main__":
    main()
