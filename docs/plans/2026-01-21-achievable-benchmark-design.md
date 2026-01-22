# Achievable Ground Truth Benchmarking & Module 3 Debug

**Date:** 2026-01-21
**Status:** Design approved, ready for implementation
**Branch:** `hierarchical_approach`

## Problem Statement

CITEgeist is being evaluated against 10 granular cell types, but some types are **indistinguishable** with the 27-antibody panel:

| Current Performance | CITEgeist | Cell2Location |
|---------------------|-----------|---------------|
| Pearson r (10 types) | 0.15 | 0.27 |
| JSD | 0.63 | 0.51 |

**Root Cause:** We're penalizing CITEgeist for failing to distinguish cell types that are fundamentally indistinguishable with available markers.

**Indistinguishable types (with 27-antibody panel):**
- Stromal + Vascular Stromal: Both VIM+, missing DCN/LUM/PDGFRB to distinguish
- Myofibroblasts + Stromal: Both VIM+, alphaSMA varies but overlaps
- Mixed Immune: CD3E+HLA-DR interface, ambiguous
- Proliferating T: CD3E+PCNA, overlaps with CD8+ T cells

## Solution: 7 Achievable Cell Types

After analyzing Xenium single-cell protein expression:

| Cluster | Cell Type | alphaSMA mean | Vimentin mean |
|---------|-----------|---------------|---------------|
| 5 | Myofibroblasts | 108 (high) | 374 (high) |
| 6 | Stromal | 20 (low) | 198 (high) |

Both express high Vimentin, so Vimentin alone cannot distinguish Stromal. Merging to "Fibroblasts" is the cleanest solution.

### Final 7 Achievable Cell Types

| Cell Type | Major Markers | Minor Markers | GT Types Absorbed |
|-----------|---------------|---------------|-------------------|
| B cells | CD20 | CD45RA | B cells |
| CD4+ T cells | CD3E, CD4 | CD45RO | Mixed Immune |
| CD8+ T cells | CD3E, CD8A | GranzymeB | CD8+ T cells, Proliferating T |
| Macrophages | CD68, CD163 | CD16 | Macrophages |
| Endothelial | CD31 | | Endothelial, Vascular Stromal |
| Epithelial | PanCK, E-Cadherin | Beta-catenin | Epithelial |
| Fibroblasts | alphaSMA, Vimentin | | Myofibroblasts, Stromal |

### GT Collapse Mapping (10 → 7)

```python
GT_TO_ACHIEVABLE_7_MAPPING = {
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

ACHIEVABLE_7_CELL_TYPES = [
    "B cells",
    "CD4+ T cells",
    "CD8+ T cells",
    "Macrophages",
    "Endothelial",
    "Epithelial",
    "Fibroblasts",
]

ACHIEVABLE_7_CELL_PROFILE_DICT = {
    "B cells": {
        "Major": ["CD20"],
        "Minor": ["CD45RA"],
    },
    "CD4+ T cells": {
        "Major": ["CD3E", "CD4"],
        "Minor": ["CD45RO"],
    },
    "CD8+ T cells": {
        "Major": ["CD3E", "CD8A"],
        "Minor": ["GranzymeB"],
    },
    "Macrophages": {
        "Major": ["CD68", "CD163"],
        "Minor": ["CD16"],
    },
    "Endothelial": {
        "Major": ["CD31"],
        "Minor": [],
    },
    "Epithelial": {
        "Major": ["PanCK", "E-Cadherin"],
        "Minor": ["Beta-catenin"],
    },
    "Fibroblasts": {
        "Major": ["alphaSMA", "Vimentin"],
        "Minor": [],
    },
}
```

## Implementation Plan

### Phase 1: Establish Fair Benchmark Baseline

**Step 1a: Run CITEgeist with 7-type profiles**

Create new benchmark run with `ACHIEVABLE_7_CELL_PROFILE_DICT`:
- Output: `output_achievable_7/`
- Uses 7 cell type profiles as input to Module 3

**Step 1b: Evaluate all methods against 7-type GT**

| Method | Input Profiles | Evaluate Against |
|--------|----------------|------------------|
| CITEgeist (7-type manual) | ACHIEVABLE_7_CELL_PROFILE_DICT | 7-type collapsed GT |
| CITEgeist (autodiscovery) | Auto-discovered | 7-type collapsed GT |
| Cell2Location | Reference signatures | 7-type collapsed GT |
| RCTD | Reference signatures | 7-type collapsed GT |
| Tangram | Reference signatures | 7-type collapsed GT |
| Seurat | Reference signatures | 7-type collapsed GT |

**Step 1c: Analyze gaps**

| Comparison | What it tells us |
|------------|------------------|
| CITEgeist 7-type vs Cell2Location | Is Module 3 competitive on fair benchmark? |
| CITEgeist autodiscovery vs 7-type manual | How much does Module 2 cost? |

### Phase 2: Autodiscovery Profile Mapping

For autodiscovery evaluation, compare two mapping approaches:

**Approach A: Existing mapping → collapse to 7**
```
Autodiscovery Profile → 10-type celltype → 7-type achievable
```

**Approach B: Direct re-mapping to 7 types**
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

def map_profile_to_achievable_7(profile_markers: List[str]) -> str:
    """Match profile to best achievable type by Jaccard similarity."""
    best_match = "Unknown"
    best_score = 0
    for celltype, signature in ACHIEVABLE_7_MARKER_SIGNATURES.items():
        overlap = len(set(profile_markers) & signature)
        union = len(set(profile_markers) | signature)
        jaccard = overlap / union if union > 0 else 0
        if jaccard > best_score:
            best_score = jaccard
            best_match = celltype
    return best_match
```

Compare both approaches to understand mapping impact.

### Phase 3: Module 3 Enhancement (if needed)

Based on Phase 1-2 findings:
- If Module 3 competitive with fair benchmark → Done
- If still underperforming → Investigate beta optimization, spatial regularization

## File Changes Required

### 1. `evaluation/src/evaluate_all_methods.py`

Add 7-type configurations:

```python
GT_TO_ACHIEVABLE_7_MAPPING = {
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

ACHIEVABLE_7_CELL_TYPES = [
    "B cells", "CD4+ T cells", "CD8+ T cells", "Macrophages",
    "Endothelial", "Epithelial", "Fibroblasts",
]

def collapse_gt_to_achievable_7(gt_df: pd.DataFrame) -> pd.DataFrame:
    """Collapse 10 GT cell types to 7 achievable types."""
    # Similar to existing collapse_gt_to_achievable but with 7-type mapping
    ...

# Add method configs
"CITEgeist_achievable_7": {
    "output_dir": "CITEgeist/output_achievable_7",
    "uses_achievable_7_gt": True,
},
"CITEgeist_autodiscovery_achievable_7": {
    "output_dir": "CITEgeist/output_autodiscovery",
    "uses_achievable_7_gt": True,
},
"Cell2Location_achievable_7": {
    "output_dir": "Cell2Location/output_granular",
    "uses_achievable_7_gt": True,
},
# ... etc for RCTD, Tangram, Seurat
```

### 2. `CITEgeist/slurm/xenium_benchmark_achievable_7.sh` (new)

SLURM script to run CITEgeist with 7-type profiles.

### 3. `evaluation/slurm/evaluate_achievable_7.sh` (new)

SLURM script to evaluate all methods against 7-type GT.

## Success Criteria

1. **Match manual profile quality** - Autodiscovery performs within 0.05 JSD of manual profiles
2. **Beat or match Cell2Location** - CITEgeist r >= Cell2Location r on 7-type benchmark
3. **Reduce variance** - Std of Pearson r across regions < 0.15
4. **Biologically interpretable** - Discovered profiles map cleanly to 7 achievable types

## Output Structure

```
evaluation/
├── results_granular/           # Current (10-type GT)
├── results_achievable/         # 8-type GT (existing, if any)
├── results_achievable_7/       # NEW (7-type GT)
│   ├── full_results.json
│   ├── method_summary.csv
│   ├── comparison_table.csv
│   └── mapping_comparison.csv  # Approach A vs B for autodiscovery
```

## Timeline Estimate

- Phase 1a (run 7-type benchmark): 2-4 hours (SLURM)
- Phase 1b (evaluate all methods): 1-2 hours (SLURM)
- Phase 1c (analysis): 1 hour
- Phase 2 (autodiscovery mapping): 2 hours
- Phase 3 (if needed): TBD based on findings

## Previous Attempts (Context)

| Approach | Outcome |
|----------|---------|
| EM Joint Optimization | Failed - JSD +59%, r -60% |
| NMF Profile Discovery | Replaced with hierarchical |
| Over-fragmentation fix | Improved modularity 0.19→0.68 |
| 8-type achievable | Defined but not fully benchmarked |

This design addresses the fair evaluation problem first before attempting algorithmic changes.
