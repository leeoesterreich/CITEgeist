#!/usr/bin/env python
"""Aggregate Module 3.5 benchmark results across all 5 Xenium regions.

Reads per-region outputs from output_module3_5_benchmark_array/ and produces:
  summary_table.csv     -- one row per (region, cell_type, marker) with all metrics
  benchmark_summary.csv -- one row per region with pass/fail and pair counts

If singlecell_validation.json exists in a region dir, sc_* columns are joined.

Usage:
    python aggregate_module3_5_results.py
    python aggregate_module3_5_results.py --benchmark-dir /path/to/dir
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent.parent.parent.parent.parent
DEFAULT_BENCHMARK_DIR = (
    REPO_ROOT / "Benchmarking" / "xenium_benchmarking" / "CITEgeist" / "output_module3_5_benchmark_array"
)
N_REGIONS = 5


def load_region(benchmark_dir: Path, region_id: int) -> tuple[dict | None, pd.DataFrame | None]:
    """Load summary JSON and pair metrics CSV for one region."""
    region_dir = benchmark_dir / f"Xenium_region_{region_id}"
    summary_path = region_dir / "module3_5_summary.json"
    metrics_path = region_dir / "module3_5_pair_metrics.csv"
    if not summary_path.exists() or not metrics_path.exists():
        print(f"[WARN] Region {region_id}: missing files in {region_dir}", file=sys.stderr)
        return None, None
    summary = json.loads(summary_path.read_text())
    metrics = pd.read_csv(metrics_path)
    metrics["region_id"] = region_id
    return summary, metrics


def load_singlecell_validation(benchmark_dir: Path, region_id: int) -> pd.DataFrame | None:
    """Load single-cell validation JSON for one region, if present."""
    sc_path = benchmark_dir / f"Xenium_region_{region_id}" / "singlecell_validation.json"
    if not sc_path.exists():
        return None
    data = json.loads(sc_path.read_text())
    pairs = data.get("pairs", [])
    if not pairs:
        return None
    df = pd.DataFrame(pairs)
    df["region_id"] = region_id
    # sc_* columns already prefixed by the validation script — no rename needed
    return df


def build_tables(benchmark_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return (summary_table, benchmark_summary) DataFrames."""
    all_metrics: list[pd.DataFrame] = []
    region_rows: list[dict] = []

    for region_id in range(N_REGIONS):
        summary, metrics = load_region(benchmark_dir, region_id)
        if summary is None:
            continue
        region_rows.append(
            {
                "region_id": region_id,
                "benchmark_passed": summary.get("benchmark_passed", False),
                "n_pairs": summary.get("n_pairs", 0),
                "n_evaluable_pairs": summary.get("n_evaluable_pairs", 0),
                "n_validated_pairs": len(summary.get("validated_pairs", [])),
                "n_headline_failures": len(summary.get("headline_failures", [])),
            }
        )
        sc_df = load_singlecell_validation(benchmark_dir, region_id)
        if sc_df is not None:
            join_on = [
                c
                for c in ("pair_id", "cell_type", "marker", "region_id")
                if c in metrics.columns and c in sc_df.columns
            ]
            metrics = metrics.merge(sc_df, on=join_on, how="left")
        all_metrics.append(metrics)

    summary_table = pd.concat(all_metrics, ignore_index=True) if all_metrics else pd.DataFrame()
    benchmark_summary = pd.DataFrame(region_rows)
    return summary_table, benchmark_summary


def print_table(benchmark_summary: pd.DataFrame) -> None:
    print("\n=== Module 3.5 Benchmark Summary ===")
    print(f"{'Region':<8} {'Passed':<8} {'Pairs':<8} {'Eval':<8} {'Valid':<8} {'HeadFail'}")
    print("-" * 52)
    for _, row in benchmark_summary.iterrows():
        status = "PASS" if row["benchmark_passed"] else "FAIL"
        print(
            f"{int(row['region_id']):<8} {status:<8} {int(row['n_pairs']):<8} "
            f"{int(row['n_evaluable_pairs']):<8} {int(row['n_validated_pairs']):<8} "
            f"{int(row['n_headline_failures'])}"
        )
    overall = benchmark_summary["benchmark_passed"].all()
    n_pass = int(benchmark_summary["benchmark_passed"].sum())
    n_total = len(benchmark_summary)
    print("-" * 52)
    print(f"Overall: {'PASS' if overall else 'FAIL'} ({n_pass}/{n_total} regions passed)\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--benchmark-dir", type=Path, default=DEFAULT_BENCHMARK_DIR)
    args = parser.parse_args()

    if not args.benchmark_dir.exists():
        print(f"ERROR: {args.benchmark_dir} does not exist", file=sys.stderr)
        sys.exit(1)

    summary_table, benchmark_summary = build_tables(args.benchmark_dir)
    out_summary = args.benchmark_dir / "summary_table.csv"
    out_bench = args.benchmark_dir / "benchmark_summary.csv"
    summary_table.to_csv(out_summary, index=False)
    benchmark_summary.to_csv(out_bench, index=False)
    print(f"Wrote {out_summary}")
    print(f"Wrote {out_bench}")
    print_table(benchmark_summary)


if __name__ == "__main__":
    main()
