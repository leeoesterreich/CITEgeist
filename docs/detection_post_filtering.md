# Detection-Based Post-Filtering for Hybrid Model

## Overview

Detection-based post-filtering improves the hybrid CITEgeist model by **+8% relative** (0.687 → 0.743 Pearson r) by correcting spurious cell type allocations after discretization.

## The Problem

The continuous CITEgeist model sometimes allocates small proportions to cell types that aren't actually present in a spot. After discretization, these become 1-2 cells that shouldn't exist.

## The Solution

Apply confidence-weighted filtering **after** discretization:

1. Run GMM detection to identify which cell types are present per spot
2. For non-detected cell types, reduce their proportion by a reliability factor (r²)
3. Renormalize proportions to sum to 1

```python
# For each cell type
if not detected[spot, cell_type]:
    proportion[spot, cell_type] *= (1 - reliability[cell_type])

# Renormalize
proportions = proportions / proportions.sum(axis=1)
```

## Critical: Post-Processing vs Pre-Processing

| Approach | When Applied | Effect | Mean r |
|----------|--------------|--------|--------|
| **Pre-filtering** | Before discretization | **Hurts** | 0.63 (-0.05) |
| **Post-filtering** | After discretization | **Helps** | 0.74 (+0.06) |

**Why pre-filtering fails:**
- Filtering modifies the continuous proportions
- Discretization then re-allocates cells to fill the total
- The filtered types get cells back during discretization
- Net effect: model's learned allocations are broken

**Why post-filtering works:**
- Discretization preserves the model's learned allocations
- Filtering corrects obvious errors (cells in non-detected types)
- No re-allocation occurs after filtering

## Reliability Values

Per-cell-type reliability (r²) derived from detection-estimation model validation:

| Cell Type | r² | Interpretation |
|-----------|-----|----------------|
| B cells | 0.85 | High confidence - strongly filter non-detected |
| Endothelial | 0.66 | Good confidence |
| CD8+ T cells | 0.63 | Good confidence |
| Epithelial | 0.60 | Moderate confidence |
| Macrophages | 0.50 | Moderate confidence |
| Fibroblasts | 0.45 | Lower confidence |
| CD4+ T cells | 0.35 | Low confidence - weak filtering |

Higher r² means more aggressive filtering for non-detected spots.

## Results

Validated on 5 Xenium pseudo-Visium regions (Job 1602642):

| Region | Baseline | Post-Filtered | Improvement |
|--------|----------|---------------|-------------|
| 0 | 0.6766 | 0.7371 | +0.0605 |
| 1 | 0.6156 | 0.6916 | +0.0760 |
| 2 | 0.6673 | 0.7280 | +0.0607 |
| 3 | 0.7373 | 0.7841 | +0.0468 |
| 4 | 0.7379 | 0.7733 | +0.0354 |
| **Mean** | **0.6869** | **0.7428** | **+0.056** |

**+8% relative improvement** over baseline hybrid.

## Usage

Post-filtering is **enabled by default** in `benchmark_hybrid_cellpose.py`.

```bash
# Default: with post-filtering (recommended)
python benchmark_hybrid_cellpose.py --region 0 --output-dir ./output

# Disable post-filtering (not recommended)
python benchmark_hybrid_cellpose.py --region 0 --output-dir ./output --no-detection-filter
```

## Implementation

The filtering is implemented in Step 5.5 of `benchmark_hybrid_cellpose.py`:

```python
# Step 5.5: Detection-based post-filtering
if use_detection_filter:
    # Run GMM detection
    detected = detect_cell_types(ab_matrix, marker_groups, adaptive_threshold=True)

    # Apply confidence-weighted filtering
    for col in proportions_df.columns:
        r2 = DETECTION_RELIABILITY.get(col, 0.5)
        not_detected = ~detected_df[col].values
        proportions_df.loc[not_detected, col] *= (1 - r2)

    # Renormalize
    proportions_df = proportions_df.div(proportions_df.sum(axis=1), axis=0)
```

## Dependencies

- `CITEgeist.model.detection.detect_cell_types` - GMM-based cell type detection
- Requires antibody data with markers for each cell type

## References

- Detection algorithm: `docs/detection_estimation_algorithm.md`
- GMM implementation: `CITEgeist/model/detection.py`

---

*Added: 2026-02-27*
