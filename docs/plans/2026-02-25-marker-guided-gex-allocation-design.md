# Marker-Guided GEX Allocation Design

**Date**: 2026-02-25
**Status**: Approved
**Author**: Brainstorming session

## Problem Statement

GEX deconvolution (Module 3 Pass 2) suffers from uniform spreading due to:
1. L2 regularization pulling allocations toward uniform distribution
2. 80/20 smoothing diluting enrichment signal
3. Weak proportion-based enrichment for non-marker genes

**Current performance:**
- Marker genes: r = 0.64
- General genes: r = 0.27
- Variance ratio: 3.0 (3x spurious variance)

## Design Decisions

| Decision | Choice | Rationale |
|----------|--------|-----------|
| Marker source | RNA anchor genes (gex_modules.py) | 50-100 genes per type, proven r>0.5 |
| Aggregation | Weighted mean by anchor strength | Trusted anchors (high r) count more |
| Spatial scale | Local neighborhood | Captures tissue heterogeneity |
| Fallback | Always use marker enrichment | Simplicity; let optimization handle weak signals |
| Integration | Adaptive blend (enrichment variance) | Trust proportions when confident, use anchors when needed |
| Scope | Remove L2 + 80/20 smoothing | Clean break from problematic components |

## Architecture Overview

**What changes:**
1. `compute_expression_aware_enrichment()` → new `compute_adaptive_marker_enrichment()`
2. Remove L2 penalty term from Gurobi objective
3. Remove 80/20 uniform smoothing

**What stays the same:**
- Gurobi QP structure with non-negativity and sum constraints
- Neighborhood-based optimization scope
- Per-gene allocation loop
- Output format (cell type × gene matrix)

**Data flow:**
```
Module 3 Pass 1 (proportions)
    ↓
discover_anchor_genes() → anchor_genes dict {cell_type: [gene_list]}
    ↓
Module 3 Pass 2 (GEX) with adaptive marker enrichment
```

## Algorithm

### Step 1: Compute proportion-based enrichment (existing, minus smoothing)

```python
def compute_proportion_enrichment(gene_expr, cell_type_props, neighborhood):
    """Same as current, but WITHOUT 80/20 smoothing."""
    weights = gene_expr / (gene_expr.sum() + 1e-10)
    weighted_props = np.sum(cell_type_props * weights[:, np.newaxis], axis=0)
    background_props = np.mean(cell_type_props, axis=0)
    enrichment = weighted_props / (background_props + 1e-10)
    return enrichment / (enrichment.sum() + 1e-10)  # normalize
```

### Step 2: Compute marker-guided enrichment (new)

```python
def compute_marker_enrichment(gene_expr, anchor_expr, anchor_weights, neighborhood):
    """
    gene_expr: (N_neighbors,) expression of target gene
    anchor_expr: (N_neighbors, T) mean anchor expression per cell type
    anchor_weights: (T,) weight per cell type = mean(anchor correlations with proportions)
    """
    enrichment = np.zeros(T)
    for t in range(T):
        r = pearsonr(gene_expr, anchor_expr[:, t])[0]
        enrichment[t] = max(0, r) * anchor_weights[t]  # weighted by anchor strength

    if enrichment.sum() < 1e-10:
        return np.ones(T) / T  # fallback to uniform
    return enrichment / enrichment.sum()
```

### Step 3: Adaptive blend (per-gene)

```python
def compute_adaptive_enrichment(gene_expr, cell_type_props, anchor_expr, anchor_weights):
    prop_enrich = compute_proportion_enrichment(gene_expr, cell_type_props)
    marker_enrich = compute_marker_enrichment(gene_expr, anchor_expr, anchor_weights)

    # Enrichment variance: high = peaked (trust proportions), low = flat (use anchors)
    # Variance calculated per gene across cell types
    variance = np.var(prop_enrich)
    max_variance = 0.25  # theoretical max for normalized distribution
    anchor_weight = 1 - min(1, variance / max_variance)

    return (1 - anchor_weight) * prop_enrich + anchor_weight * marker_enrich
```

**Key detail**: Variance is calculated per gene. Each gene gets its own `anchor_weight`:
- **CD19**: `prop_enrich = [0.6, 0.05, 0.05, ...]` → high variance → low anchor_weight
- **COL3A1**: `prop_enrich = [0.15, 0.14, 0.14, ...]` → low variance → high anchor_weight

## Gurobi Objective Modification

**Current objective:**
```python
obj = gp.quicksum(
    enrichment[k, j] * center_props[j] * X[j, k]
    for j in range(T) for k in range(M)
) - lambda_gex_reg * gp.quicksum(
    X[j, k] * X[j, k] for j in range(T) for k in range(M)
)
```

**New objective (remove L2 entirely):**
```python
obj = gp.quicksum(
    enrichment[k, j] * center_props[j] * X[j, k]
    for j in range(T) for k in range(M)
)
```

**Constraints remain unchanged:**
- Non-negativity: `X[j, k] >= 0`
- Sum constraint: `sum_j(X[j, k]) == total_counts[k]`

**Note**: We keep multiplying by `center_props[j]` to respect spot composition. This is distinct from the enrichment calculation and ensures we don't allocate to absent cell types. May revisit if issues arise.

## Integration Points

### Files to modify:

1. **`gurobi_impl.py`**
   - Add `compute_marker_enrichment()` function
   - Add `compute_adaptive_enrichment()` function
   - Modify `deconvolute_spot_with_neighbors_with_prior()`:
     - Accept `anchor_expr` and `anchor_weights` parameters
     - Replace enrichment calculation with adaptive version
     - Remove L2 penalty from objective
     - Remove 80/20 smoothing

2. **`citegeist_model.py`**
   - Modify `run_cell_expression_pass1()`:
     - Call `discover_anchor_genes()` before GEX allocation (or accept pre-computed)
     - Compute per-cell-type anchor expression matrix
     - Pass to deconvolution function

3. **`gex_modules.py`** (minimal changes)
   - Ensure `discover_anchor_genes()` returns anchor weights (correlation with proportions)
   - May already do this; need to verify

### New parameters for `run_cell_expression_pass1()`:

```python
def run_cell_expression_pass1(
    ...,
    use_marker_guidance: bool = True,  # enable/disable for A/B testing
    anchor_genes: dict = None,  # optional pre-computed anchors
)
```

## Testing Strategy

**Evaluation metrics** (using existing `evaluate_gex_spatial.py`):
- `per_gene_pearson`: target >0.35 (up from 0.27)
- `variance_ratio`: target <2.0 (down from 3.0)
- `marker_pearson`: maintain >0.60

**Test plan:**

1. **Unit test**: Verify `compute_adaptive_enrichment()` produces expected outputs
   - High-variance input → low anchor weight
   - Low-variance input → high anchor weight
   - Correlation with synthetic anchors works correctly

2. **Regression test**: Run on 1 Xenium region (region 2, not the outlier region 1)
   - Compare against baseline hybrid_cellpose
   - Verify marker genes don't degrade

3. **Full benchmark**: Run on all 5 regions via SLURM
   - Use existing evaluation harness
   - Compare all metrics against baseline

**A/B toggle**: `use_marker_guidance=True/False` allows direct comparison without code changes.

**Success criteria:**
- per_gene_pearson improvement of at least +0.05 (0.27 → 0.32+)
- variance_ratio reduction (3.0 → <2.5)
- marker_pearson stable (no more than -0.02 degradation)

## Future Considerations

1. **Proportion multiplication in objective**: Currently kept; may test removing if double-weighting suspected
2. **Explicit upper bounds**: If extreme allocations become an issue without L2, can add `X[j,k] <= 2 * total * prop[j]`
3. **Region 1 outlier**: Has variance_ratio=8.37; may need separate investigation

## Related Documents

- `docs/plans/2026-02-24-gex-allocation-methods-analysis.md` - Previous experiments
- `docs/plans/2026-02-25-gex-l2-alternatives-exploration.md` - Root cause analysis
