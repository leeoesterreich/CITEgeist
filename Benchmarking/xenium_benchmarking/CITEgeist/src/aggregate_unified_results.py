#!/usr/bin/env python
"""Aggregate unified GEX benchmark results across regions/replicates.

Reads per-region JSON files and computes mean ± std for each method.
Works for both Xenium (6 methods, percell_r + spot_layer_r) and
simulation (5 methods, spot_layer_r only) benchmarks.

Usage:
    # Xenium
    python aggregate_unified_results.py \
        --input_dir results_unified_gex/ \
        --prefix xenium_unified_gex_region \
        --n 5

    # Simulation
    python aggregate_unified_results.py \
        --input_dir results_unified_gex/ \
        --prefix simulation_unified_gex_rep \
        --n 5
"""

import argparse
import json
from pathlib import Path

import numpy as np


def aggregate(input_dir, prefix, n):
    """Load per-region/rep JSONs and aggregate."""
    input_dir = Path(input_dir)
    records = []
    missing = []

    for i in range(n):
        path = input_dir / f"{prefix}_{i}.json"
        if path.exists():
            with open(path) as f:
                records.append(json.load(f))
        else:
            missing.append(i)

    if not records:
        print(f"ERROR: No result files found in {input_dir} with prefix {prefix}")
        return

    if missing:
        print(f"WARNING: Missing indices: {missing}")

    # Discover methods (exclude metadata keys)
    meta_keys = {
        "region",
        "replicate",
        "n_common_genes",
        "sace_params",
        "sace_diagnostics",
        "sace_runtime_seconds",
    }
    methods = [k for k in records[0] if k not in meta_keys]

    # Aggregate scalar metrics per method
    agg = {"n_records": len(records), "missing": missing, "methods": {}}

    for method in methods:
        method_agg = {}
        # Collect all scalar metrics
        sample_rec = records[0][method]
        for metric_key in sample_rec:
            if isinstance(sample_rec[metric_key], (int, float)):
                values = []
                for rec in records:
                    if method in rec and metric_key in rec[method]:
                        v = rec[method][metric_key]
                        if isinstance(v, (int, float)) and np.isfinite(v):
                            values.append(float(v))
                if values:
                    method_agg[metric_key] = {
                        "mean": round(float(np.mean(values)), 4),
                        "std": round(float(np.std(values)), 4),
                        "n": len(values),
                        "values": [round(v, 4) for v in values],
                    }
            elif isinstance(sample_rec[metric_key], dict):
                # Per-type metrics
                sub_agg = {}
                for type_name in sample_rec[metric_key]:
                    values = []
                    for rec in records:
                        if method in rec and metric_key in rec[method] and type_name in rec[method][metric_key]:
                            v = rec[method][metric_key][type_name]
                            if isinstance(v, (int, float)) and np.isfinite(v):
                                values.append(float(v))
                    if values:
                        sub_agg[type_name] = {
                            "mean": round(float(np.mean(values)), 4),
                            "std": round(float(np.std(values)), 4),
                            "n": len(values),
                        }
                if sub_agg:
                    method_agg[metric_key] = sub_agg

        agg["methods"][method] = method_agg

    # Save aggregate JSON
    out_name = prefix.rsplit("_", 1)[0] + "_aggregate.json"
    out_path = input_dir / out_name
    with open(out_path, "w") as f:
        json.dump(agg, f, indent=2)
    print(f"Saved: {out_path}")

    # Print summary table
    print()
    print("=" * 70)
    # Detect if this is Xenium (has percell) or simulation (spot_layer only)
    has_percell = any("percell_mean_gene_r" in agg["methods"].get(m, {}) for m in methods)

    if has_percell:
        print(f"{'Method':<18} {'percell_r':>14} {'spot_layer_r':>14}")
        print("-" * 48)
        for m in methods:
            ma = agg["methods"].get(m, {})
            pc = ma.get("percell_mean_gene_r", {})
            sl = ma.get("spot_layer_overall_r", {})
            pc_str = f"{pc.get('mean', 0):.4f}±{pc.get('std', 0):.4f}" if pc else "N/A"
            sl_str = f"{sl.get('mean', 0):.4f}±{sl.get('std', 0):.4f}" if sl else "N/A"
            print(f"{m:<18} {pc_str:>14} {sl_str:>14}")
    else:
        print(f"{'Method':<18} {'spot_layer_r':>14}")
        print("-" * 34)
        for m in methods:
            ma = agg["methods"].get(m, {})
            sl = ma.get("spot_layer_overall_r", {})
            sl_str = f"{sl.get('mean', 0):.4f}±{sl.get('std', 0):.4f}" if sl else "N/A"
            print(f"{m:<18} {sl_str:>14}")

    print("=" * 70)


def main():
    parser = argparse.ArgumentParser(description="Aggregate unified GEX benchmark results")
    parser.add_argument("--input_dir", required=True, help="Directory with per-region/rep JSONs")
    parser.add_argument(
        "--prefix",
        required=True,
        help="File prefix (e.g., xenium_unified_gex_region or simulation_unified_gex_rep)",
    )
    parser.add_argument("--n", type=int, default=5, help="Number of regions/replicates")
    args = parser.parse_args()
    aggregate(args.input_dir, args.prefix, args.n)


if __name__ == "__main__":
    main()
