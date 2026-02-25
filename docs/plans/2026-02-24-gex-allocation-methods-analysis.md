# GEX Allocation Methods Analysis

**Date**: 2026-02-24
**Status**: Analysis Complete, Next Steps Identified

## Problem Statement

GEX deconvolution (Module 3 Pass 2) allocates total gene expression at each spot to cell types based on their proportions. The question: how should we distribute gene counts when multiple cell types are present?

## Methods Tested

### 1. L2 Regularization (Baseline)

**Implementation**: Gurobi QP with L2 penalty on allocation variables
```
minimize: ||allocation - target||² + λ * ||allocation||²
```

**Result**: Creates **uniform spreading** across all cell types regardless of enrichment.

| Metric | Value |
|--------|-------|
| per_gene_pearson | 0.273 |
| marker_pearson | 0.636 |
| variance_ratio | **3.02** (BAD - 3x spurious variance) |

**Problem**: L2 penalty is convex and symmetric - it penalizes deviation from zero equally for all cell types, causing genes to spread ~14.3% to each of 7 cell types regardless of biological enrichment.

### 2. Softmax Allocation (Expression-Weighted)

**Implementation**: Direct proportional allocation without Gurobi
```python
# Step 1: Compute enrichment per gene per cell type
enrichment = weighted_props / background_props  # expression-weighted

# Step 2: Softmax allocation
log_enrichment = np.log(enrichment + 1e-10)
weights = softmax(log_enrichment / temperature)
allocation = total_counts * weights
```

**Key parameters**:
- `temperature = 0.3` (lower = more concentrated)
- Expression-weighted enrichment (no percentile threshold)

| Metric | Value | vs Baseline |
|--------|-------|-------------|
| per_gene_pearson | 0.288 | **+5.5%** |
| marker_pearson | 0.644 | +1.2% |
| variance_ratio | **1.53** | **-49%** (much closer to 1.0) |

**Improvement**: Cuts spurious variance in half by concentrating allocation to enriched cell types.

## Key Insight: The 2.3x Gap

**Marker genes**: r = 0.64
**General genes**: r = 0.27
**Gap**: 2.3x

Both methods show this same gap. Why?

### Root Cause Analysis

The enrichment calculation uses **cell type proportions**:
```python
# What we're doing:
enrichment[gene, celltype] = f(proportions_where_gene_is_expressed)

# What this means:
# - High expression in spots with high B-cell proportion → allocate to B cells
# - Works great for CD19 (B-cell marker) because CD19 expression = B-cell presence
# - Fails for weak genes that don't have cell-type-specific expression patterns
```

**The fundamental problem**: A gene with weak/uniform expression across cell types will have low enrichment variance. Even perfect allocation won't create spatial pattern that doesn't exist in the original data.

## What Would Close the Gap

The original pipeline vision:
1. **Protein → Proportions** ✓ (working well, r=0.73)
2. **Proportions → Marker Genes** ✓ (markers have r=0.64)
3. **Marker Genes → Gene Modules** ✗ (NOT IMPLEMENTED)

### Missing Step: Marker-Gene-Guided Allocation

For weak genes, instead of using proportion-based enrichment, use **spatial correlation with marker genes**:

```python
# CURRENT (proportion-based):
enrichment[gene, B_cells] = correlation(gene_expression, B_cell_proportion)

# PROPOSED (marker-guided):
enrichment[gene, B_cells] = correlation(gene_expression, CD19_expression)
```

**Why this would work**:
- CD19 has strong spatial pattern (proven: r=0.70)
- BANK1 might not correlate with B-cell proportions directly
- But BANK1 might spatially co-express with CD19
- Use CD19's spatial pattern to guide BANK1's allocation

## Benchmark Results Summary

### Per-Cell-Type Marker Gene Correlation (Softmax T=0.3)

| Cell Type | Markers | r |
|-----------|---------|---|
| B_cells | CD19, MS4A1, CD79A | 0.74 |
| Endothelial | PECAM1, VWF | 0.74 |
| CD8+ T cells | CD3E, CD8A, GZMB | 0.71 |
| Macrophages | CD68, CD163 | 0.65 |
| Fibroblasts | PDGFRA | 0.58 |
| CD4+ T cells | CD3E, CD4, IL7R | 0.54 |
| Epithelial | EPCAM, KRT8 | 0.53 |

### Top Individual Markers

| Marker | Cell Type | r |
|--------|-----------|---|
| MS4A1 | B_cells | 0.87 |
| VWF | Endothelial | 0.83 |
| CD8A | CD8+ T cells | 0.81 |
| IL7R | CD4+ T cells | 0.78 |
| CD19 | B_cells | 0.77 |
| CD68 | Macrophages | 0.71 |

## Files Created

- `Benchmarking/xenium_benchmarking/evaluation/src/evaluate_gex_spatial.py`
  - Spatial-aware GEX evaluation
  - per_gene_pearson, per_spot_pearson, variance_ratio metrics
  - Marker gene analysis per cell type

## Next Steps

1. **Design doc needed**: Marker-gene-guided allocation mode
   - How to identify marker genes per cell type (top N by proportion correlation)
   - How to compute spatial correlation efficiently
   - How to weight multiple markers per cell type

2. **Implementation**: Add `allocation_mode="marker_guided"` to `deconvolute_spot_with_neighbors_with_prior()`

3. **Expected improvement**: Weak genes that co-express with markers should gain ~0.1-0.2 correlation boost

## Code Reference

Softmax allocation (not yet merged, reverted):
```python
# In deconvolute_spot_with_neighbors_with_prior()
if allocation_mode == "softmax":
    from scipy.special import softmax as scipy_softmax

    result = np.zeros((T, M))
    center_counts = deconvolution_expression_data[spot_idx, :]

    for k in range(M):
        total_counts = center_counts[k]
        if total_counts <= 0:
            continue

        enrichment = gene_specific_enrichment[k]

        if enrichment.sum() > 0 and softmax_temperature > 0:
            log_enrichment = np.log(enrichment + 1e-10)
            weights = scipy_softmax(log_enrichment / softmax_temperature)
            result[:, k] = total_counts * weights
        else:
            result[:, k] = total_counts / T

    return result
```

Expression-weighted enrichment (replaces percentile-based):
```python
def compute_expression_aware_enrichment(expression_data, cell_type_props, gene_idx):
    gene_expr = expression_data[:, gene_idx]
    total_expr = gene_expr.sum()

    if total_expr < 1e-10:
        return np.ones(cell_type_props.shape[1]) / cell_type_props.shape[1]

    # Weight each spot's contribution by its expression level
    weights = gene_expr / total_expr

    # Normalize by global cell type frequency
    normalized_props = cell_type_props / (celltype_frequencies + 1e-10)

    # Weighted mean of cell type proportions
    weighted_props = np.sum(normalized_props * weights[:, np.newaxis], axis=0)
    background_props = np.mean(normalized_props, axis=0)

    enrichment = weighted_props / (background_props + 1e-10)
    return enrichment / (np.sum(enrichment) + 1e-10)
```
