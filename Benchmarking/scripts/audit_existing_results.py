#!/usr/bin/env python3
"""
Audit existing benchmark results across the Benchmarking directory.

Crawls for *results*.json and benchmark_results.json files, extracts key
metadata, and generates audit_report.csv showing what exists.
"""

import csv
import hashlib
import json
import os
import re
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional

# Directories to search
BENCHMARKING_ROOT = Path(__file__).parent.parent
SEARCH_DIRS = [
    BENCHMARKING_ROOT / "xenium_benchmarking",
    BENCHMARKING_ROOT / "simulation_benchmarking",
]

# Output
OUTPUT_CSV = BENCHMARKING_ROOT / "audit_report.csv"


def extract_method_from_path(file_path: Path) -> str:
    """Infer method name from file path."""
    path_str = str(file_path)
    filename = file_path.stem

    # Check for method names in path or filename
    methods = ["CITEgeist", "Cell2Location", "RCTD", "Tangram", "Seurat", "CARD", "scResolve"]

    # First check filename for explicit method names
    for method in methods:
        if method.lower() in filename.lower():
            # Check for variants in filename
            if "hybrid" in filename.lower():
                return f"{method}_Hybrid"
            elif "discrete" in filename.lower():
                return f"{method}_Discrete"
            elif "continuous" in filename.lower():
                return f"{method}_Continuous"
            return method

    # Then check path for method directories
    for method in methods:
        if f"/{method.lower()}/" in path_str.lower() or f"/{method}/" in path_str:
            # Check for variants in path
            if "hybrid" in path_str.lower():
                return f"{method}_Hybrid"
            elif "discrete" in path_str.lower():
                return f"{method}_Discrete"
            elif "continuous" in path_str.lower():
                return f"{method}_Continuous"
            return method

    # If filename has _results suffix, extract method name
    if "_results" in filename:
        method_candidate = filename.replace("_results", "")
        # Clean up common suffixes
        method_candidate = re.sub(r'_region_\d+', '', method_candidate)
        if method_candidate and method_candidate not in ["full", "combined", "benchmark"]:
            return method_candidate

    return "Unknown"


def extract_dataset_type(file_path: Path) -> str:
    """Determine if this is xenium or simulation benchmark."""
    path_str = str(file_path)
    if "xenium" in path_str.lower():
        return "xenium"
    elif "simulation" in path_str.lower():
        return "simulation"
    return "unknown"


def extract_ground_truth_type(file_path: Path, content: Any) -> str:
    """Determine ground truth type (protein_gated vs rna_based)."""
    path_str = str(file_path)

    # Check path patterns
    if "protein_gt" in path_str.lower() or "protein_gated" in path_str.lower():
        return "protein_gated"
    elif "rna_gt" in path_str.lower() or "rna_based" in path_str.lower():
        return "rna_based"
    elif "granular" in path_str.lower():
        return "rna_granular"

    # Check for method_comparison directory (uses protein-gated GT)
    if "method_comparison" in path_str:
        return "protein_gated"

    # Handle list content
    if isinstance(content, list):
        return "unknown"

    # Check content for hints
    if isinstance(content, dict) and "summary" in content:
        summary = content["summary"]
        if isinstance(summary, dict):
            cell_types = summary.get("cell_types", [])
            # Protein-gated uses specific naming with + signs
            if any("CD4+ T cells" in ct or "CD8+ T cells" in ct for ct in cell_types):
                return "protein_gated"
            # Simulation uses different naming
            if any("B-cells" in ct or "CAFs" in ct for ct in cell_types):
                return "simulation_gt"

    # Check for simulation-specific fields
    if isinstance(content, dict) and ("replicate_id" in content or "condition" in content):
        return "simulation_gt"

    return "unknown"


def extract_metrics(content: Any) -> Dict[str, Optional[float]]:
    """Extract key metrics from results JSON."""
    metrics = {
        "prop_pearson_overall": None,
        "gex_pearson": None,
        "n_celltypes": None,
    }

    # Handle list content
    if isinstance(content, list):
        # Try to extract from first element if it's a dict with metrics
        if len(content) > 0 and isinstance(content[0], dict):
            if "pearson_r" in content[0]:
                pearson_vals = [r.get("pearson_r") for r in content if isinstance(r, dict) and "pearson_r" in r]
                if pearson_vals:
                    metrics["gex_pearson"] = sum(pearson_vals) / len(pearson_vals)
        return metrics

    # Try summary format (from evaluate_benchmark.py output)
    if "summary" in content:
        summary = content["summary"]
        metrics["prop_pearson_overall"] = summary.get("overall_mean_pearson_r")
        if "cell_types" in summary:
            metrics["n_celltypes"] = len(summary.get("cell_types", []))

    # Try direct format (simulation benchmarks)
    if "proportion_correlation" in content:
        metrics["prop_pearson_overall"] = content["proportion_correlation"]
    if "overall_mean_pearson_r" in content:
        metrics["prop_pearson_overall"] = content["overall_mean_pearson_r"]

    # Try results nested format (full_comparison files)
    if "results" in content and isinstance(content["results"], dict):
        for method_name, method_data in content["results"].items():
            if isinstance(method_data, dict) and "summary" in method_data:
                metrics["prop_pearson_overall"] = method_data["summary"].get("overall_mean_pearson_r")
                break
            # GEX comparison format - list of per-region results
            if isinstance(method_data, list) and len(method_data) > 0:
                # Calculate average pearson_r across regions
                pearson_vals = [r.get("pearson_r") for r in method_data if "pearson_r" in r]
                if pearson_vals:
                    metrics["gex_pearson"] = sum(pearson_vals) / len(pearson_vals)
                break

    # GEX metrics - various formats
    if "gex_metrics" in content:
        metrics["gex_pearson"] = content["gex_metrics"].get("pearson_r")
    elif "gex_pearson_r" in content:
        metrics["gex_pearson"] = content["gex_pearson_r"]
    elif "gex_average_rmse" in content:
        # Simulation format - has GEX but no single pearson
        metrics["gex_pearson"] = None  # Could compute from per-celltype if needed

    # Count cell types from various sources
    if metrics["n_celltypes"] is None:
        if "gex_per_celltype" in content:
            metrics["n_celltypes"] = len(content["gex_per_celltype"])
        elif "regions" in content and len(content["regions"]) > 0:
            first_region = content["regions"][0]
            if "cell_types" in first_region:
                metrics["n_celltypes"] = len(first_region["cell_types"])

    return metrics


def compute_settings_hash(content: Any) -> str:
    """Compute a short hash of the settings for deduplication."""
    # Handle case where content is a list (some files are arrays)
    if isinstance(content, list):
        settings_str = json.dumps(content[:3] if len(content) > 3 else content, sort_keys=True)
        return hashlib.md5(settings_str.encode()).hexdigest()[:8]

    # Extract settings or use key identifying fields
    settings_dict = content.get("settings", {})

    # If no settings, use some identifying fields
    if not settings_dict:
        identifying = {
            "n_regions": content.get("summary", {}).get("n_regions") if isinstance(content.get("summary"), dict) else None,
            "total_spots": content.get("summary", {}).get("total_spots") if isinstance(content.get("summary"), dict) else None,
            "cell_types": content.get("summary", {}).get("cell_types") if isinstance(content.get("summary"), dict) else None,
            "condition": content.get("condition"),
            "mode": content.get("mode"),
        }
        settings_dict = {k: v for k, v in identifying.items() if v is not None}

    settings_str = json.dumps(settings_dict, sort_keys=True)
    return hashlib.md5(settings_str.encode()).hexdigest()[:8]


def find_result_files() -> List[Path]:
    """Find all result JSON files in benchmark directories."""
    patterns = ["*results*.json", "benchmark_results.json", "full_comparison*.json"]
    files = []

    for search_dir in SEARCH_DIRS:
        if not search_dir.exists():
            print(f"Warning: Directory not found: {search_dir}")
            continue
        for pattern in patterns:
            found = list(search_dir.rglob(pattern))
            files.extend(found)

    # Deduplicate and sort
    return sorted(set(files))


def audit_file(file_path: Path) -> Optional[Dict[str, Any]]:
    """Extract audit information from a single results file."""
    try:
        with open(file_path) as f:
            content = json.load(f)
    except (json.JSONDecodeError, IOError) as e:
        return {
            "file_path": str(file_path.relative_to(BENCHMARKING_ROOT)),
            "method": "PARSE_ERROR",
            "dataset": extract_dataset_type(file_path),
            "ground_truth": "unknown",
            "n_celltypes": None,
            "has_proportion": "error",
            "has_gex": "error",
            "prop_pearson_overall": None,
            "gex_pearson": None,
            "run_date": None,
            "settings_hash": None,
            "error": str(e),
        }

    metrics = extract_metrics(content)

    try:
        mtime = datetime.fromtimestamp(file_path.stat().st_mtime).strftime("%Y-%m-%d")
    except:
        mtime = None

    return {
        "file_path": str(file_path.relative_to(BENCHMARKING_ROOT)),
        "method": extract_method_from_path(file_path),
        "dataset": extract_dataset_type(file_path),
        "ground_truth": extract_ground_truth_type(file_path, content),
        "n_celltypes": metrics["n_celltypes"],
        "has_proportion": "yes" if metrics["prop_pearson_overall"] is not None else "no",
        "has_gex": "yes" if metrics["gex_pearson"] is not None else "no",
        "prop_pearson_overall": metrics["prop_pearson_overall"],
        "gex_pearson": metrics["gex_pearson"],
        "run_date": mtime,
        "settings_hash": compute_settings_hash(content),
        "error": None,
    }


def main():
    print(f"Auditing benchmark results in {BENCHMARKING_ROOT}")
    print("=" * 70)

    files = find_result_files()
    print(f"Found {len(files)} result files")

    audit_results = []
    for f in files:
        result = audit_file(f)
        if result:
            audit_results.append(result)

    # Write CSV
    if audit_results:
        fieldnames = list(audit_results[0].keys())
        with open(OUTPUT_CSV, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(audit_results)
        print(f"\nWrote {len(audit_results)} rows to {OUTPUT_CSV}")

    # Summary statistics
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)

    xenium_protein = [r for r in audit_results if r["dataset"] == "xenium" and r["ground_truth"] == "protein_gated"]
    xenium_rna = [r for r in audit_results if r["dataset"] == "xenium" and r["ground_truth"] not in ("protein_gated", "unknown")]
    xenium_unknown = [r for r in audit_results if r["dataset"] == "xenium" and r["ground_truth"] == "unknown"]
    simulation = [r for r in audit_results if r["dataset"] == "simulation"]
    errors = [r for r in audit_results if r["error"]]

    print(f"Xenium (protein-gated):   {len(xenium_protein):3d} files")
    print(f"Xenium (RNA/other):       {len(xenium_rna):3d} files")
    print(f"Xenium (unknown GT):      {len(xenium_unknown):3d} files")
    print(f"Simulation:               {len(simulation):3d} files")
    print(f"Parse errors:             {len(errors):3d} files")

    # Show methods found
    methods_found = set(r["method"] for r in audit_results if r["method"] != "PARSE_ERROR")
    print(f"\nMethods detected: {sorted(methods_found)}")

    # Show xenium protein-gated files (what we care about)
    if xenium_protein:
        print(f"\n--- Xenium Protein-Gated Results (canonical benchmark) ---")
        for r in sorted(xenium_protein, key=lambda x: (x["method"], x["file_path"])):
            prop_r = f"{r['prop_pearson_overall']:.4f}" if r['prop_pearson_overall'] else "N/A"
            gex_r = f"{r['gex_pearson']:.4f}" if r['gex_pearson'] else "N/A"
            print(f"  {r['method']:<24} prop_r={prop_r:<8} gex_r={gex_r:<8} {r['file_path']}")

    # Show parse errors if any
    if errors:
        print(f"\n--- Parse Errors ---")
        for r in errors:
            print(f"  {r['file_path']}: {r['error'][:60]}")


if __name__ == "__main__":
    main()
