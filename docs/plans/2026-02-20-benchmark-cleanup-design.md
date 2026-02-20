# Benchmark Cleanup Design

**Date:** 2026-02-20
**Status:** Approved
**Goal:** Consolidate scattered benchmark results into a single source of truth

## Problem Statement

The benchmarking directory is a mess with 100+ scattered result files across 12+ experiment directories. Different agents report different Pearson r values (0.4 to 0.7) because there's no canonical benchmark. This blocks meaningful model improvement work.

## Requirements

### Primary Benchmark
- **Dataset:** Xenium pseudo-Visium
- **Ground truth:** Protein-gated (not RNA-based)
- **Cell types:** 7 "achievable" panel
- **Metrics:** Per-cell-type proportion Pearson r + GEX Pearson r (all genes)
- **Methods:** CITEgeist (hybrid) vs Cell2Location, RCTD, Tangram, Seurat, etc.

### Infrastructure
- Single `canonical_results.json` as source of truth
- Per-method evaluation scripts with standardized interface
- Consolidation script that merges results into canonical JSON

### Cleanup Approach
- Audit existing experiments first (Approach 3: Hybrid Audit)
- Archive valuable experiments, delete garbage
- Re-run only what's necessary

---

## Design

### Section 1: Canonical JSON Schema

```json
{
  "metadata": {
    "generated_at": "2026-02-20T12:00:00Z",
    "benchmark_dataset": "xenium_pseudovisium_protein_gated",
    "n_spots": 7054,
    "n_regions": 5,
    "cell_types": ["B_cells", "CD4_T", "CD8_T", "Endothelial", "Epithelial", "Fibroblasts", "Macrophages"]
  },
  "methods": {
    "CITEgeist_Hybrid": {
      "version": "commit_hash_or_tag",
      "run_date": "2026-02-18",
      "settings": {"preprocessing": "CLR", "discretization": "largest_remainder"},
      "proportion_metrics": {
        "per_celltype": {
          "B_cells": {"pearson_r": 0.68, "rmse": 0.10},
          "CD4_T": {"pearson_r": 0.39, "rmse": 0.12}
        },
        "overall": {"pearson_r": 0.726, "rmse": 0.106}
      },
      "gex_metrics": {
        "pearson_r": 0.38,
        "rmse_cpm": 4.43,
        "n_genes": 500
      },
      "source_file": "path/to/raw/results.json"
    }
  }
}
```

**Design decisions:**
- `per_celltype` metrics are first-class (primary metric of interest)
- `source_file` provides provenance back to raw results
- `settings` captures method-specific config for reproducibility
- `version` tracks code version

### Section 2: Audit Script

**Location:** `Benchmarking/scripts/audit_existing_results.py`

**Input:** Walks `Benchmarking/` recursively for `*results*.json` and `benchmark_results.json`

**Output:** `Benchmarking/audit_report.csv` with columns:
- `file_path`
- `method`
- `dataset` (xenium vs simulation)
- `ground_truth` (protein_gated vs rna_based)
- `n_celltypes`
- `has_proportion` (yes/no)
- `has_gex` (yes/no)
- `prop_pearson_overall`
- `gex_pearson`
- `run_date`
- `settings_hash`

**Logic:**
1. Find all JSON files matching result patterns
2. Extract method name, dataset type, ground truth type, key metrics
3. Flag files that can't be parsed or have missing fields
4. Generate summary statistics

### Section 3: Per-Method Evaluation Script Interface

**Location:** `Benchmarking/xenium_benchmarking/{method}/src/evaluate_{method}.py`

**Interface:**
```bash
python evaluate_citegeist.py \
  --predictions /path/to/predictions.csv \
  --ground-truth /path/to/protein_gated_gt.csv \
  --output /path/to/results.json \
  --method-name "CITEgeist_Hybrid" \
  --settings '{"preprocessing": "CLR"}'
```

**Output format:** Method-level subset of canonical schema

**Ground truth files (standardized):**
- `Benchmarking/xenium_benchmarking/ground_truth/protein_gated_proportions.csv`
- `Benchmarking/xenium_benchmarking/ground_truth/protein_gated_gex.csv`

**Constraint:** All methods MUST use the same ground truth files.

### Section 4: Consolidation Script

**Location:** `Benchmarking/scripts/consolidate_results.py`

**Usage:**
```bash
python consolidate_results.py \
  --input-dir Benchmarking/xenium_benchmarking/evaluation/results/method_comparison/ \
  --output Benchmarking/canonical_results.json \
  --ground-truth-meta Benchmarking/xenium_benchmarking/ground_truth/metadata.json
```

**Logic:**
1. Read ground truth metadata
2. Scan input directory for `*_results.json`
3. Validate each file matches schema
4. Merge into canonical structure
5. Compute summary statistics and rank methods
6. Write `canonical_results.json`
7. Print human-readable summary

**Property:** Idempotent (same inputs = same output)

### Section 5: Directory Structure

```
Benchmarking/
├── canonical_results.json          # THE source of truth
├── audit_report.csv                # Output of audit script
│
├── scripts/
│   ├── audit_existing_results.py
│   └── consolidate_results.py
│
├── xenium_benchmarking/
│   ├── ground_truth/
│   │   ├── metadata.json
│   │   ├── protein_gated_proportions.csv
│   │   └── protein_gated_gex.csv
│   │
│   ├── CITEgeist/
│   │   ├── src/
│   │   │   └── evaluate_citegeist.py
│   │   ├── predictions/
│   │   └── results/
│   │
│   ├── Cell2Location/
│   │   └── ... (same structure)
│   │
│   └── ... (RCTD, Tangram, Seurat, etc.)
│
├── simulation_benchmarking/        # Keep, lower priority
│
└── _archive/                       # Old experiments
    └── 2026-02-20_pre_cleanup/
```

---

## Success Criteria

1. `canonical_results.json` exists and contains all methods
2. Any agent asked "what's CITEgeist's performance?" reads from canonical JSON
3. Running consolidation script twice produces identical output
4. Old experiment directories are archived, not cluttering active paths
5. Each method's evaluation can be re-run independently

## Non-Goals

- Improving model performance (that's after cleanup)
- Simulation benchmarking (lower priority, keep but don't clean yet)
- Automated CI/CD for benchmarks (future work)
