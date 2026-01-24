# Xenium Benchmark Cleanup & Hierarchical Validation Design

**Date:** 2026-01-23
**Status:** Design approved, ready for implementation
**Branch:** `hierarchical_approach`

## Problem Statement

The Xenium benchmarking codebase has accumulated complexity from iterative development:
- 10 CITEgeist output directories with different configurations
- 7 evaluation results directories
- 18 SLURM scripts with overlapping functionality
- 16 Python source files with duplication

The achievable-7 benchmark shows CITEgeist (manual profiles) beats Cell2Location (r=0.41 vs r=0.38), but the codebase needs cleanup before:
1. Making achievable-7 the default for fair method comparison
2. Validating hierarchical autodiscovery on Xenium regions
3. Preparing artifacts for paper publication

## Goals

1. **Clean canonical structure** - Two modes: `manual` (achievable-7) and `hierarchical`
2. **Archive, don't delete** - Preserve old work in `_archive/` subdirectories
3. **Unified entry points** - Single `run_benchmark.py` and `evaluate_benchmark.py`
4. **Hierarchical validation** - Test autodiscovery on all 14 Xenium regions
5. **Paper-ready artifacts** - Figures and reports for publication

## Design

### Canonical Directory Structure

```
Benchmarking/xenium_benchmarking/
├── CITEgeist/
│   ├── src/
│   │   ├── run_benchmark.py          # Unified entry point (--mode manual|hierarchical)
│   │   ├── analyze_concordance.py    # Keep - useful analysis
│   │   └── _archive/                  # Old scripts
│   │       ├── run_benchmark_autodiscovery.py
│   │       ├── run_full_pipeline.py
│   │       ├── compare_benchmark_modes.py
│   │       └── evaluate_profile_discovery.py
│   ├── slurm/
│   │   ├── run_benchmark.sh          # Main SLURM script (--mode flag)
│   │   ├── run_hierarchical_benchmark.sh  # Keep for now (hierarchical-specific)
│   │   └── _archive/                  # Old SLURM scripts
│   │       ├── xenium_benchmark.sh
│   │       ├── xenium_benchmark_granular.sh
│   │       ├── xenium_benchmark_achievable.sh
│   │       ├── xenium_benchmark_rna_gt.sh
│   │       ├── run_all_benchmarks.sh
│   │       ├── run_autodiscovery_benchmark.sh
│   │       ├── compare_benchmarks.sh
│   │       └── evaluate_profile_discovery.sh
│   └── output/
│       ├── manual/                    # Achievable-7 results (renamed from output_achievable_7)
│       ├── hierarchical/              # Autodiscovery results (new)
│       └── _archive/                  # Old output directories
│           ├── autodiscovery/
│           ├── granular/
│           ├── achievable/
│           ├── profile_discovery/
│           └── rna_gt/
├── evaluation/
│   ├── src/
│   │   ├── evaluate_benchmark.py     # Unified evaluation
│   │   ├── evaluate_metrics.py       # Keep - metric implementations
│   │   └── _archive/
│   │       ├── evaluate_all_methods.py
│   │       ├── compare_with_scresolve.py
│   │       └── compare_achievable_7_with_scresolve.py
│   └── results/
│       ├── method_comparison/         # CITEgeist vs Cell2Location vs others
│       ├── autodiscovery_validation/  # Hierarchical profile quality analysis
│       └── _archive/                  # Old results directories
│           ├── results_granular/
│           ├── results_achievable_7/
│           └── results_with_scresolve/
└── {Cell2Location,RCTD,Tangram,Seurat}/  # Other methods unchanged
```

### Unified run_benchmark.py Interface

```python
"""
Unified CITEgeist benchmark runner for Xenium pseudo-Visium data.

Usage:
    python run_benchmark.py --region-id 0 --mode manual       # Achievable-7 profiles
    python run_benchmark.py --region-id 0 --mode hierarchical # Auto-discover with hierarchy
"""

# Canonical achievable-7 profiles (fair benchmark baseline)
ACHIEVABLE_7_CELL_PROFILE_DICT = {
    "B cells": {"Major": ["CD20"], "Minor": ["CD45RA"]},
    "CD4+ T cells": {"Major": ["CD3E", "CD4"], "Minor": ["CD45RO"]},
    "CD8+ T cells": {"Major": ["CD3E", "CD8A"], "Minor": ["GranzymeB"]},
    "Macrophages": {"Major": ["CD68", "CD163"], "Minor": ["CD16"]},
    "Endothelial": {"Major": ["CD31"], "Minor": []},
    "Epithelial": {"Major": ["PanCK", "E-Cadherin"], "Minor": ["Beta-catenin"]},
    "Fibroblasts": {"Major": ["alphaSMA", "Vimentin"], "Minor": []},
}

# Key arguments
parser.add_argument("--mode", choices=["manual", "hierarchical"], default="manual")
parser.add_argument("--region-id", type=int, required=True)
parser.add_argument("--output-dir", type=str, default=None)  # Auto-set based on mode
parser.add_argument("--run-gex", action="store_true")

# Hierarchical-specific options
parser.add_argument("--fdr-alpha", type=float, default=0.05)
parser.add_argument("--top-k", type=int, default=10)
parser.add_argument("--improvement-threshold", type=float, default=0.05)
```

### Ground Truth and Evaluation

**GT Collapse Mapping (10 → 7):**

```python
GT_COLLAPSE_MAPPING = {
    "B cells": "B cells",
    "Mixed Immune": "CD4+ T cells",
    "CD8+ T cells": "CD8+ T cells",
    "Proliferating T": "CD8+ T cells",
    "Macrophages": "Macrophages",
    "Endothelial": "Endothelial",
    "Vascular Stromal": "Endothelial",
    "Epithelial": "Epithelial",
    "Myofibroblasts": "Fibroblasts",
    "Stromal": "Fibroblasts",
}
```

**Autodiscovery Compression (N → 7):**

For hierarchical mode, autodiscovered profiles get compressed to achievable-7 for fair benchmarking:

```python
ACHIEVABLE_7_MARKER_SIGNATURES = {
    "B cells": {"CD20", "CD45RA"},
    "CD4+ T cells": {"CD3E", "CD4", "CD45RO"},
    "CD8+ T cells": {"CD3E", "CD8A", "GranzymeB"},
    "Macrophages": {"CD68", "CD163", "CD16"},
    "Endothelial": {"CD31"},
    "Epithelial": {"PanCK", "E-Cadherin", "Beta-catenin"},
    "Fibroblasts": {"alphaSMA", "Vimentin"},
}

def compress_to_achievable_7(discovered_profiles: Dict) -> Dict:
    """
    Map N discovered profiles → 7 achievable types via Jaccard similarity.

    Multiple discovered profiles may map to same achievable type.
    Proportions get summed when compressed.

    Example:
      - "Profile_1" (CD3, CD4, CD45RO) → "CD4+ T cells"
      - "Profile_2" (CD3, CD8A) → "CD8+ T cells"
      - "Profile_3" (CD3, CD8A, GranzymeB) → "CD8+ T cells"  # Also maps here
    """
    compressed = {}
    for profile_name, markers in discovered_profiles.items():
        best_match = max(
            ACHIEVABLE_7_MARKER_SIGNATURES.keys(),
            key=lambda ct: jaccard(set(markers), ACHIEVABLE_7_MARKER_SIGNATURES[ct])
        )
        compressed[profile_name] = best_match
    return compressed
```

**Key insight:** Autodiscovery finding more granular profiles than 7 is a STRENGTH of CITEgeist, not a failure. The compression mapping enables fair benchmarking while the raw profiles showcase discovery capability.

### Two-Tier Evaluation

| Tier | Purpose | Comparison |
|------|---------|------------|
| **Fair benchmark** | Method comparison | Compress autodiscovered → achievable-7, compare to Cell2Location etc. |
| **Discovery showcase** | Highlight CITEgeist strength | Show the richer profiles autodiscovery finds |

### Hierarchical Validation Artifacts

```
results/autodiscovery_validation/
├── discovered_profiles_raw.json         # Full granularity (N profiles per region)
├── discovered_profiles_compressed.json  # Mapped to achievable-7
├── granularity_analysis.csv             # How many profiles per achievable type?
├── profile_consistency_matrix.csv       # Cross-region profile similarity
├── marker_hierarchy_summary.csv         # Which markers are shared vs specific
├── novel_populations.md                 # Profiles that don't map cleanly (interesting!)
└── figures/
    ├── profile_hierarchy_tree.png       # Full discovered hierarchy
    ├── compression_sankey.png           # N profiles → 7 types flow diagram
    ├── marker_sharing_heatmap.png       # Shared marker patterns
    └── region_consistency.png           # Cross-region agreement
```

### Validation Questions

| Question | What it shows |
|----------|---------------|
| How many profiles does autodiscovery find? | CITEgeist discovers richer structure |
| Which achievable types get subdivided? | E.g., T cells split into 3+ subtypes |
| Are there novel populations? | Profiles that don't map to any achievable type |
| Does hierarchy detect T cell subtypes? | CD3 shared, CD4/CD8 distinguish subtypes |
| Does hierarchy detect stromal subtypes? | Vimentin shared, αSMA distinguishes |
| Are profiles consistent across regions? | >70% profile overlap across regions |
| Does compression still beat other methods? | Fair benchmark with autodiscovery |

## Implementation Plan

### Phase 1: Cleanup

**Step 1: Create archive directories**
```bash
mkdir -p CITEgeist/src/_archive
mkdir -p CITEgeist/slurm/_archive
mkdir -p CITEgeist/output/_archive
mkdir -p evaluation/src/_archive
mkdir -p evaluation/results/_archive
```

**Step 2: Archive obsolete files**

Move to `CITEgeist/src/_archive/`:
- run_benchmark_autodiscovery.py
- run_full_pipeline.py
- compare_benchmark_modes.py
- evaluate_profile_discovery.py

Move to `CITEgeist/slurm/_archive/`:
- xenium_benchmark.sh
- xenium_benchmark_granular.sh
- xenium_benchmark_achievable.sh
- xenium_benchmark_rna_gt.sh
- run_all_benchmarks.sh
- run_autodiscovery_benchmark.sh
- compare_benchmarks.sh
- evaluate_profile_discovery.sh

Move to `CITEgeist/output/_archive/`:
- autodiscovery/ (rename from output_autodiscovery)
- granular/ (rename from output_granular)
- achievable/ (rename from output_achievable)
- profile_discovery/ (rename from output_profile_discovery)
- rna_gt/ (rename from output_rna_gt)

Move to `evaluation/src/_archive/`:
- evaluate_all_methods.py
- compare_with_scresolve.py
- compare_achievable_7_with_scresolve.py

Move to `evaluation/results/_archive/`:
- results_granular/
- results_achievable_7/ (copy, keep as reference)
- results_with_scresolve/

**Step 3: Consolidate run_benchmark.py**
- Merge run_benchmark_hierarchical.py logic into run_benchmark.py
- Add `--mode manual|hierarchical` flag
- Move ACHIEVABLE_7_CELL_PROFILE_DICT to top as canonical reference
- Default output dirs: `output/manual/` and `output/hierarchical/`

**Step 4: Rename output directories**
- output_achievable_7/ → output/manual/

**Step 5: Update SLURM script**
- Rename xenium_benchmark_achievable_7.sh → run_benchmark.sh
- Add --mode parameter support

**Step 6: Update evaluation**
- Create evaluate_benchmark.py with compression mapping
- Default to achievable-7 GT for all methods

### Phase 2: Validation (after cleanup)

1. Run hierarchical discovery on all 14 Xenium regions
2. Collect discovered profiles per region
3. Analyze consistency across regions
4. Generate validation figures and reports
5. Compare compressed results against other methods

## Files Changed

| Action | File | Notes |
|--------|------|-------|
| Archive | ~15 Python/bash files | Move to _archive/ |
| Archive | ~6 output directories | Move to output/_archive/ |
| Consolidate | run_benchmark.py | Merge hierarchical logic |
| Rename | output_achievable_7/ | → output/manual/ |
| Rename | xenium_benchmark_achievable_7.sh | → run_benchmark.sh |
| Create | evaluate_benchmark.py | Unified evaluation with compression |
| Create | output/hierarchical/ | New autodiscovery output |
| Create | results/method_comparison/ | New results structure |
| Create | results/autodiscovery_validation/ | Validation artifacts |

## Success Criteria

1. **Clean structure** - Single entry point for each task (run, evaluate)
2. **Backward compatible** - Old results preserved in archives
3. **Fair benchmark** - Manual achievable-7 beats Cell2Location (r > 0.38)
4. **Hierarchical works** - Discovers biologically sensible profiles
5. **Compression works** - Hierarchical compressed to 7 types performs within 0.10 JSD of manual

---

**Last Updated:** 2026-01-23
**Author:** Claude + Alex
