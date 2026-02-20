#!/usr/bin/env python3
"""
Consolidate method results into canonical_results.json.

Reads individual method result files from xenium_benchmarking evaluation
and produces a single canonical JSON with standardized schema.

Usage:
    python consolidate_results.py
    python consolidate_results.py --input-dir /path/to/results --output /path/to/output.json
"""

import argparse
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np

BENCHMARKING_ROOT = Path(__file__).parent.parent


def get_git_commit() -> str:
    """Get current git commit hash (first 12 characters)."""
    try:
        result = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            capture_output=True,
            text=True,
            cwd=BENCHMARKING_ROOT,
        )
        return result.stdout.strip()[:12]
    except Exception:
        return "unknown"


def load_method_result(file_path: Path) -> Dict[str, Any]:
    """Load a single method result JSON file."""
    with open(file_path) as f:
        return json.load(f)


def extract_per_celltype_metrics(summary: Dict) -> Dict[str, Dict[str, float]]:
    """
    Extract per-cell-type Pearson r and std from summary.

    Looks for keys like "B cells_mean_r" and "B cells_std_r".
    """
    per_celltype = {}
    cell_types = summary.get("cell_types", [])

    for ct in cell_types:
        mean_r_key = f"{ct}_mean_r"
        std_r_key = f"{ct}_std_r"

        if mean_r_key in summary:
            mean_val = summary[mean_r_key]
            std_val = summary.get(std_r_key, 0.0)

            # Handle NaN values
            if mean_val is not None and not (isinstance(mean_val, float) and np.isnan(mean_val)):
                per_celltype[ct] = {
                    "pearson_r": mean_val,
                    "pearson_r_std": std_val if std_val and not np.isnan(std_val) else 0.0,
                }

    return per_celltype


def convert_to_canonical_format(
    method_name: str,
    result: Dict,
    source_file: str,
    settings: Optional[Dict] = None,
) -> Dict[str, Any]:
    """
    Convert a method result to canonical schema format.

    Args:
        method_name: Name of the method (e.g., "CITEgeist_Hybrid")
        result: Raw result dict loaded from JSON
        source_file: Relative path to source file
        settings: Method-specific settings to document

    Returns:
        Canonical format dict with standardized schema
    """
    summary = result.get("summary", {})

    canonical = {
        "version": get_git_commit(),
        "run_date": datetime.now(timezone.utc).strftime("%Y-%m-%d"),
        "settings": settings or {},
        "proportion_metrics": {
            "per_celltype": extract_per_celltype_metrics(summary),
            "overall": {
                "pearson_r": summary.get("overall_mean_pearson_r"),
                "pearson_r_std": summary.get("overall_std_pearson_r"),
                "rmse": summary.get("overall_mean_rmse"),
                "mae": summary.get("overall_mean_mae"),
                "jsd": summary.get("overall_mean_jsd"),
            },
        },
        "gex_metrics": None,  # Will be populated if GEX results exist
        "source_file": source_file,
    }

    return canonical


def load_gex_results(results_dir: Path) -> Dict[str, Dict]:
    """
    Load GEX comparison results from full_comparison_gex.json.

    The GEX file contains per-cell-type, per-region results.
    We compute the overall mean Pearson r and RMSE.

    Args:
        results_dir: Directory containing full_comparison_gex.json

    Returns:
        Dict mapping method name to GEX metrics
    """
    gex_file = results_dir / "full_comparison_gex.json"
    if not gex_file.exists():
        return {}

    with open(gex_file) as f:
        data = json.load(f)

    gex_metrics = {}

    for method, results_list in data.get("results", {}).items():
        if not results_list:
            continue

        # Compute aggregate metrics across all cell types and regions
        pearson_rs = [r.get("pearson_r", 0) for r in results_list if r.get("pearson_r") is not None]
        rmses = [r.get("rmse", 0) for r in results_list if r.get("rmse") is not None]
        raw_rmses = [r.get("raw_rmse", 0) for r in results_list if r.get("raw_rmse") is not None]
        n_genes = results_list[0].get("n_genes") if results_list else None

        gex_metrics[method] = {
            "pearson_r": float(np.mean(pearson_rs)) if pearson_rs else None,
            "rmse_cpm": float(np.mean(rmses)) if rmses else None,
            "raw_rmse": float(np.mean(raw_rmses)) if raw_rmses else None,
            "n_genes": n_genes,
        }

    return gex_metrics


def main():
    parser = argparse.ArgumentParser(
        description="Consolidate benchmark results into canonical JSON format"
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=BENCHMARKING_ROOT / "xenium_benchmarking" / "evaluation" / "results" / "method_comparison",
        help="Directory containing method result JSONs",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=BENCHMARKING_ROOT / "canonical_results.json",
        help="Output canonical results file",
    )
    parser.add_argument(
        "--metadata",
        type=Path,
        default=BENCHMARKING_ROOT / "xenium_benchmarking" / "ground_truth" / "metadata.json",
        help="Ground truth metadata file",
    )
    args = parser.parse_args()

    print(f"Consolidating results from: {args.input_dir}")
    print(f"Output: {args.output}")

    # Load ground truth metadata
    if args.metadata.exists():
        with open(args.metadata) as f:
            gt_metadata = json.load(f)
    else:
        print(f"Warning: Metadata file not found at {args.metadata}, using defaults")
        gt_metadata = {
            "benchmark_dataset": "xenium_pseudovisium_protein_gated",
            "n_spots_total": 7054,
            "n_regions": 5,
            "cell_types": [
                "B cells", "CD4+ T cells", "CD8+ T cells",
                "Macrophages", "Endothelial", "Epithelial", "Fibroblasts"
            ],
        }

    # Build canonical structure
    canonical = {
        "metadata": {
            "generated_at": datetime.now(timezone.utc).isoformat(),
            "benchmark_dataset": gt_metadata.get("benchmark_dataset", "xenium_pseudovisium_protein_gated"),
            "n_spots": gt_metadata.get("n_spots_total", 7054),
            "n_regions": gt_metadata.get("n_regions", 5),
            "cell_types": gt_metadata.get("cell_types", []),
            "ground_truth_type": gt_metadata.get("ground_truth_type", "protein_gated"),
        },
        "methods": {},
    }

    # Method settings documentation
    method_settings = {
        "CITEgeist_Hybrid": {"preprocessing": "CLR", "discretization": "largest_remainder"},
        "CITEgeist_Continuous": {"preprocessing": "CLR", "optimization": "QP"},
        "CITEgeist_Discrete": {"preprocessing": "per_marker", "optimization": "IQP"},
        "Cell2Location": {"n_cells_per_location": 30},
        "Tangram": {"mode": "cells"},
        "RCTD": {"mode": "full"},
        "Seurat": {"method": "label_transfer"},
        "CARD": {"exclude_vim": True},
    }

    # Load GEX results
    gex_metrics = load_gex_results(args.input_dir)
    if gex_metrics:
        print(f"Loaded GEX metrics for: {list(gex_metrics.keys())}")

    # Find and process method result files
    result_files = sorted(args.input_dir.glob("*_results.json"))

    methods_processed = []
    methods_skipped = []

    for result_file in result_files:
        # Skip aggregate/comparison files
        if any(skip in result_file.name for skip in ["full_", "comparison", "marker_genes"]):
            methods_skipped.append(result_file.name)
            continue

        # Skip generic "CITEgeist" in favor of specific variants (Hybrid, Continuous, Discrete)
        if result_file.name == "CITEgeist_results.json":
            methods_skipped.append(f"{result_file.name} (use specific variant)")
            continue

        method_name = result_file.stem.replace("_results", "")
        print(f"  Processing: {method_name}")

        try:
            result = load_method_result(result_file)
            source_file = str(result_file.relative_to(BENCHMARKING_ROOT))
            settings = method_settings.get(method_name, {})

            method_canonical = convert_to_canonical_format(
                method_name, result, source_file, settings
            )

            # Add GEX metrics if available
            if method_name in gex_metrics:
                method_canonical["gex_metrics"] = gex_metrics[method_name]

            canonical["methods"][method_name] = method_canonical
            methods_processed.append(method_name)

        except Exception as e:
            print(f"    ERROR processing {method_name}: {e}")

    # Write output
    with open(args.output, "w") as f:
        json.dump(canonical, f, indent=2)

    print(f"\nWrote canonical results to: {args.output}")
    print(f"Methods included: {methods_processed}")
    if methods_skipped:
        print(f"Files skipped: {methods_skipped}")

    # Print summary table
    print("\n" + "=" * 70)
    print("CANONICAL BENCHMARK RESULTS")
    print("=" * 70)
    print(f"{'Method':<24} {'Prop r':>10} {'GEX r':>10} {'RMSE':>10}")
    print("-" * 54)

    # Sort by proportion Pearson r (descending)
    sorted_methods = sorted(
        canonical["methods"].items(),
        key=lambda x: x[1]["proportion_metrics"]["overall"].get("pearson_r") or 0,
        reverse=True,
    )

    for method_name, data in sorted_methods:
        prop_r = data["proportion_metrics"]["overall"].get("pearson_r")
        gex_r = data["gex_metrics"]["pearson_r"] if data["gex_metrics"] else None
        rmse = data["proportion_metrics"]["overall"].get("rmse")

        prop_str = f"{prop_r:.4f}" if prop_r is not None else "N/A"
        gex_str = f"{gex_r:.4f}" if gex_r is not None else "N/A"
        rmse_str = f"{rmse:.4f}" if rmse is not None else "N/A"

        print(f"{method_name:<24} {prop_str:>10} {gex_str:>10} {rmse_str:>10}")


if __name__ == "__main__":
    main()
