# Benchmark Cleanup Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Consolidate scattered benchmark results into a single `canonical_results.json` file.

**Architecture:** Audit existing results first, then build consolidation infrastructure, then migrate. Existing `compare_all_methods.py` and `evaluate_benchmark.py` provide the core evaluation logic - we're adding orchestration on top.

**Tech Stack:** Python 3.10, pandas, numpy, json

---

## Task 1: Create Directory Structure

**Files:**
- Create: `Benchmarking/scripts/` (directory)
- Create: `Benchmarking/xenium_benchmarking/ground_truth/` (directory)

**Step 1: Create scripts directory**

```bash
mkdir -p Benchmarking/scripts
```

**Step 2: Create ground_truth symlink**

The ground truth already exists at `Benchmarking/xenium_pseudovisium/data_protein_gt/ground_truth/`. Create a symlink for cleaner access:

```bash
mkdir -p Benchmarking/xenium_benchmarking/ground_truth
ln -sf ../../xenium_pseudovisium/data_protein_gt/ground_truth Benchmarking/xenium_benchmarking/ground_truth/proportions
ln -sf ../../xenium_pseudovisium/data_protein_gt/ground_truth_gex Benchmarking/xenium_benchmarking/ground_truth/gex
```

**Step 3: Verify structure**

```bash
ls -la Benchmarking/scripts/
ls -la Benchmarking/xenium_benchmarking/ground_truth/
```

Expected: Both directories exist, symlinks point to correct locations.

**Step 4: Commit**

```bash
git add Benchmarking/scripts/.gitkeep Benchmarking/xenium_benchmarking/ground_truth/
git commit -m "chore: create benchmark scripts directory and ground_truth symlinks"
```

---

## Task 2: Create Ground Truth Metadata File

**Files:**
- Create: `Benchmarking/xenium_benchmarking/ground_truth/metadata.json`

**Step 1: Write metadata file**

```python
# Save as Benchmarking/xenium_benchmarking/ground_truth/metadata.json
{
    "benchmark_dataset": "xenium_pseudovisium_protein_gated",
    "description": "Xenium all-genes pseudo-Visium regions with protein-gated cell type ground truth",
    "n_regions": 5,
    "n_spots_total": 7054,
    "n_spots_per_region": [1407, 1417, 1352, 1464, 1414],
    "cell_types": [
        "B cells",
        "CD4+ T cells",
        "CD8+ T cells",
        "Macrophages",
        "Endothelial",
        "Epithelial",
        "Fibroblasts"
    ],
    "ground_truth_source": "Single-cell protein gating (hierarchical)",
    "proportion_files": "proportions/Xenium_region_{0-4}_prop.csv",
    "gex_files": "gex/Xenium_region_{0-4}/"
}
```

**Step 2: Verify file is valid JSON**

```bash
python -c "import json; json.load(open('Benchmarking/xenium_benchmarking/ground_truth/metadata.json'))"
```

Expected: No error

**Step 3: Commit**

```bash
git add Benchmarking/xenium_benchmarking/ground_truth/metadata.json
git commit -m "chore: add ground truth metadata file"
```

---

## Task 3: Build Audit Script

**Files:**
- Create: `Benchmarking/scripts/audit_existing_results.py`

**Step 1: Write the audit script**

```python
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

    # Check for method names in path
    methods = ["CITEgeist", "Cell2Location", "RCTD", "Tangram", "Seurat", "CARD", "scResolve"]
    for method in methods:
        if method.lower() in path_str.lower():
            # Check for variants
            if "hybrid" in path_str.lower():
                return f"{method}_Hybrid"
            elif "discrete" in path_str.lower():
                return f"{method}_Discrete"
            elif "continuous" in path_str.lower():
                return f"{method}_Continuous"
            return method

    # Try to get from filename
    filename = file_path.stem
    if "_results" in filename:
        return filename.replace("_results", "")

    return "Unknown"


def extract_dataset_type(file_path: Path) -> str:
    """Determine if this is xenium or simulation benchmark."""
    path_str = str(file_path)
    if "xenium" in path_str.lower():
        return "xenium"
    elif "simulation" in path_str.lower():
        return "simulation"
    return "unknown"


def extract_ground_truth_type(file_path: Path, content: Dict) -> str:
    """Determine ground truth type (protein_gated vs rna_based)."""
    path_str = str(file_path)

    if "protein_gt" in path_str.lower() or "protein_gated" in path_str.lower():
        return "protein_gated"
    elif "rna_gt" in path_str.lower() or "rna_based" in path_str.lower():
        return "rna_based"
    elif "granular" in path_str.lower():
        return "rna_granular"

    # Check content for hints
    if "summary" in content:
        cell_types = content["summary"].get("cell_types", [])
        # Protein-gated uses specific naming
        if "CD4+ T cells" in cell_types or "CD8+ T cells" in cell_types:
            return "protein_gated"

    return "unknown"


def extract_metrics(content: Dict) -> Dict[str, Optional[float]]:
    """Extract key metrics from results JSON."""
    metrics = {
        "prop_pearson_overall": None,
        "gex_pearson": None,
        "n_celltypes": None,
    }

    # Try summary format (from evaluate_benchmark.py output)
    if "summary" in content:
        summary = content["summary"]
        metrics["prop_pearson_overall"] = summary.get("overall_mean_pearson_r")
        metrics["n_celltypes"] = len(summary.get("cell_types", []))

    # Try direct format
    if "overall_mean_pearson_r" in content:
        metrics["prop_pearson_overall"] = content["overall_mean_pearson_r"]

    # Try results nested format
    if "results" in content and isinstance(content["results"], dict):
        for method_name, method_data in content["results"].items():
            if "summary" in method_data:
                metrics["prop_pearson_overall"] = method_data["summary"].get("overall_mean_pearson_r")
                break

    # GEX metrics
    if "gex_metrics" in content:
        metrics["gex_pearson"] = content["gex_metrics"].get("pearson_r")
    elif "gex_pearson_r" in content:
        metrics["gex_pearson"] = content["gex_pearson_r"]

    return metrics


def compute_settings_hash(content: Dict) -> str:
    """Compute a short hash of the settings for deduplication."""
    settings_str = json.dumps(content.get("settings", {}), sort_keys=True)
    return hashlib.md5(settings_str.encode()).hexdigest()[:8]


def find_result_files() -> List[Path]:
    """Find all result JSON files in benchmark directories."""
    patterns = ["*results*.json", "benchmark_results.json", "full_comparison*.json"]
    files = []

    for search_dir in SEARCH_DIRS:
        if not search_dir.exists():
            continue
        for pattern in patterns:
            files.extend(search_dir.rglob(pattern))

    # Deduplicate
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
    mtime = datetime.fromtimestamp(file_path.stat().st_mtime).strftime("%Y-%m-%d")

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
    print("=" * 60)

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
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)

    xenium_protein = [r for r in audit_results if r["dataset"] == "xenium" and r["ground_truth"] == "protein_gated"]
    xenium_rna = [r for r in audit_results if r["dataset"] == "xenium" and r["ground_truth"] != "protein_gated"]
    simulation = [r for r in audit_results if r["dataset"] == "simulation"]

    print(f"Xenium (protein-gated): {len(xenium_protein)} files")
    print(f"Xenium (RNA/other):     {len(xenium_rna)} files")
    print(f"Simulation:             {len(simulation)} files")
    print(f"Parse errors:           {len([r for r in audit_results if r['error']])}")

    # Show xenium protein-gated files (what we care about)
    if xenium_protein:
        print(f"\n--- Xenium Protein-Gated Results (canonical benchmark) ---")
        for r in sorted(xenium_protein, key=lambda x: x["file_path"]):
            prop_r = f"{r['prop_pearson_overall']:.4f}" if r['prop_pearson_overall'] else "N/A"
            print(f"  {r['method']:<24} prop_r={prop_r:<8} {r['file_path']}")


if __name__ == "__main__":
    main()
```

**Step 2: Make executable and test**

```bash
chmod +x Benchmarking/scripts/audit_existing_results.py
cd Benchmarking && python scripts/audit_existing_results.py
```

Expected: Prints summary and creates `audit_report.csv`

**Step 3: Review audit_report.csv**

```bash
head -20 Benchmarking/audit_report.csv
```

**Step 4: Commit**

```bash
git add Benchmarking/scripts/audit_existing_results.py
git commit -m "feat: add audit script for existing benchmark results"
```

---

## Task 4: Run Audit and Analyze Results

**Step 1: Run the audit**

```bash
cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist
python Benchmarking/scripts/audit_existing_results.py
```

**Step 2: Analyze xenium protein-gated results**

```bash
grep "protein_gated" Benchmarking/audit_report.csv | grep "xenium"
```

These are the files that matter for the canonical benchmark.

**Step 3: Document findings**

Note which files should go into canonical_results.json vs which should be archived.

---

## Task 5: Build Consolidation Script

**Files:**
- Create: `Benchmarking/scripts/consolidate_results.py`

**Step 1: Write the consolidation script**

```python
#!/usr/bin/env python3
"""
Consolidate method results into canonical_results.json.

Reads individual method result files from xenium_benchmarking evaluation
and produces a single canonical JSON with standardized schema.
"""

import argparse
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List

BENCHMARKING_ROOT = Path(__file__).parent.parent


def get_git_commit() -> str:
    """Get current git commit hash."""
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
    """Load a single method result file."""
    with open(file_path) as f:
        return json.load(f)


def extract_per_celltype_metrics(summary: Dict) -> Dict[str, Dict[str, float]]:
    """Extract per-cell-type metrics from summary."""
    per_celltype = {}
    cell_types = summary.get("cell_types", [])

    for ct in cell_types:
        mean_r_key = f"{ct}_mean_r"
        std_r_key = f"{ct}_std_r"

        if mean_r_key in summary:
            per_celltype[ct] = {
                "pearson_r": summary[mean_r_key],
                "pearson_r_std": summary.get(std_r_key, 0.0),
            }

    return per_celltype


def convert_to_canonical_format(
    method_name: str,
    result: Dict,
    source_file: str,
    settings: Dict = None,
) -> Dict[str, Any]:
    """Convert a method result to canonical schema."""
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
    """Load GEX comparison results if available."""
    gex_file = results_dir / "full_comparison_gex.json"
    if not gex_file.exists():
        return {}

    with open(gex_file) as f:
        data = json.load(f)

    gex_metrics = {}
    for method, metrics in data.get("methods", {}).items():
        gex_metrics[method] = {
            "pearson_r": metrics.get("cpm_pearson_r"),
            "rmse_cpm": metrics.get("cpm_rmse"),
            "n_genes": metrics.get("n_genes"),
        }

    return gex_metrics


def main():
    parser = argparse.ArgumentParser(description="Consolidate benchmark results")
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
        gt_metadata = {
            "benchmark_dataset": "xenium_pseudovisium_protein_gated",
            "n_spots": 7054,
            "n_regions": 5,
            "cell_types": ["B cells", "CD4+ T cells", "CD8+ T cells", "Macrophages", "Endothelial", "Epithelial", "Fibroblasts"],
        }

    # Build canonical structure
    canonical = {
        "metadata": {
            "generated_at": datetime.now(timezone.utc).isoformat(),
            "benchmark_dataset": gt_metadata.get("benchmark_dataset", "xenium_pseudovisium_protein_gated"),
            "n_spots": gt_metadata.get("n_spots_total", 7054),
            "n_regions": gt_metadata.get("n_regions", 5),
            "cell_types": gt_metadata.get("cell_types", []),
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

    # Find and process method result files
    result_files = sorted(args.input_dir.glob("*_results.json"))

    for result_file in result_files:
        # Skip aggregate files
        if "full_" in result_file.name or "comparison" in result_file.name:
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

        except Exception as e:
            print(f"    ERROR: {e}")

    # Write output
    with open(args.output, "w") as f:
        json.dump(canonical, f, indent=2)

    print(f"\nWrote canonical results to: {args.output}")
    print(f"Methods included: {list(canonical['methods'].keys())}")

    # Print summary table
    print("\n" + "=" * 70)
    print("CANONICAL BENCHMARK RESULTS")
    print("=" * 70)
    print(f"{'Method':<24} {'Prop r':>10} {'GEX r':>10} {'RMSE':>10}")
    print("-" * 54)

    sorted_methods = sorted(
        canonical["methods"].items(),
        key=lambda x: x[1]["proportion_metrics"]["overall"].get("pearson_r") or 0,
        reverse=True,
    )

    for method_name, data in sorted_methods:
        prop_r = data["proportion_metrics"]["overall"].get("pearson_r")
        gex_r = data["gex_metrics"]["pearson_r"] if data["gex_metrics"] else None
        rmse = data["proportion_metrics"]["overall"].get("rmse")

        prop_str = f"{prop_r:.4f}" if prop_r else "N/A"
        gex_str = f"{gex_r:.4f}" if gex_r else "N/A"
        rmse_str = f"{rmse:.4f}" if rmse else "N/A"

        print(f"{method_name:<24} {prop_str:>10} {gex_str:>10} {rmse_str:>10}")


if __name__ == "__main__":
    main()
```

**Step 2: Test the consolidation script**

```bash
python Benchmarking/scripts/consolidate_results.py
```

Expected: Creates `Benchmarking/canonical_results.json` and prints summary table.

**Step 3: Verify canonical_results.json**

```bash
python -c "import json; d=json.load(open('Benchmarking/canonical_results.json')); print(json.dumps(d['metadata'], indent=2))"
```

**Step 4: Commit**

```bash
git add Benchmarking/scripts/consolidate_results.py
git commit -m "feat: add consolidation script for canonical results"
```

---

## Task 6: Generate Canonical Results

**Step 1: Run consolidation**

```bash
python Benchmarking/scripts/consolidate_results.py
```

**Step 2: Review output**

```bash
cat Benchmarking/canonical_results.json | python -m json.tool | head -60
```

**Step 3: Verify all methods present**

Expected methods: CITEgeist_Hybrid, Cell2Location, RCTD, Tangram, Seurat, CARD

**Step 4: Commit canonical results**

```bash
git add Benchmarking/canonical_results.json
git commit -m "feat: generate canonical benchmark results (protein-gated GT)"
```

---

## Task 7: Create Archive Directory and Document

**Files:**
- Create: `Benchmarking/_archive/.gitkeep`
- Create: `Benchmarking/_archive/README.md`

**Step 1: Create archive structure**

```bash
mkdir -p Benchmarking/_archive/2026-02-20_pre_cleanup
```

**Step 2: Write archive README**

```markdown
# Archived Benchmark Results

This directory contains benchmark results from before the cleanup on 2026-02-20.

## Why Archived

These results were generated during exploratory development with inconsistent:
- Ground truth datasets (RNA-based vs protein-gated)
- Evaluation metrics
- Method configurations

## Canonical Results

The canonical benchmark results are now in:
- `Benchmarking/canonical_results.json`

## Contents

- `2026-02-20_pre_cleanup/` - Snapshot before cleanup (not committed, for local reference only)

## Regenerating Results

To regenerate canonical results:
```bash
python Benchmarking/scripts/consolidate_results.py
```
```

**Step 3: Add to gitignore**

Add to `.gitignore`:
```
Benchmarking/_archive/2026-02-20_pre_cleanup/
```

**Step 4: Commit**

```bash
git add Benchmarking/_archive/README.md Benchmarking/_archive/.gitkeep
git commit -m "chore: create archive directory for old benchmark results"
```

---

## Task 8: Update MEMORY.md

**Files:**
- Modify: `/ihome/alee/alc376/.claude/projects/-ix1-alee-LO-LAB-Personal-Alexander-Chang-alc376-CITEgeist/memory/MEMORY.md`

**Step 1: Add canonical results section**

Add after "## Data Locations":

```markdown
## Canonical Benchmark Results (SINGLE SOURCE OF TRUTH)

**ALWAYS** use `Benchmarking/canonical_results.json` for benchmark metrics.

**DO NOT** cite metrics from:
- Individual `*_results.json` files
- `audit_report.csv`
- `full_comparison*.json` files
- Scattered `benchmark_results.json` in experiment directories

**To get metrics:**
```python
import json
with open("Benchmarking/canonical_results.json") as f:
    results = json.load(f)

# CITEgeist proportion Pearson r
prop_r = results["methods"]["CITEgeist_Hybrid"]["proportion_metrics"]["overall"]["pearson_r"]
```

**To regenerate:**
```bash
python Benchmarking/scripts/consolidate_results.py
```
```

**Step 2: Commit**

```bash
git add /ihome/alee/alc376/.claude/projects/-ix1-alee-LO-LAB-Personal-Alexander-Chang-alc376-CITEgeist/memory/MEMORY.md
git commit -m "docs: update MEMORY.md with canonical results location"
```

---

## Task 9: Final Verification

**Step 1: Run full pipeline**

```bash
# Audit
python Benchmarking/scripts/audit_existing_results.py

# Consolidate
python Benchmarking/scripts/consolidate_results.py

# Verify
python -c "
import json
with open('Benchmarking/canonical_results.json') as f:
    d = json.load(f)
print('Methods:', list(d['methods'].keys()))
print('CITEgeist_Hybrid prop_r:', d['methods']['CITEgeist_Hybrid']['proportion_metrics']['overall']['pearson_r'])
"
```

**Step 2: Verify idempotency**

```bash
# Run twice, compare
python Benchmarking/scripts/consolidate_results.py
cp Benchmarking/canonical_results.json /tmp/first.json
python Benchmarking/scripts/consolidate_results.py
diff Benchmarking/canonical_results.json /tmp/first.json
```

Expected: No diff (except generated_at timestamp)

**Step 3: Final commit**

```bash
git status
# If any uncommitted changes:
git add -A
git commit -m "chore: benchmark cleanup complete"
```

---

## Success Criteria

- [ ] `Benchmarking/canonical_results.json` exists with all methods
- [ ] `Benchmarking/audit_report.csv` shows inventory of all existing results
- [ ] `Benchmarking/scripts/` contains audit and consolidation scripts
- [ ] `Benchmarking/xenium_benchmarking/ground_truth/metadata.json` exists
- [ ] MEMORY.md updated to point to canonical results
- [ ] Running consolidation twice produces identical output (idempotent)
