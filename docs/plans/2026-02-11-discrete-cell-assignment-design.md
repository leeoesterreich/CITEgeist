# Discrete Cell Assignment Design for Module 3

**Date:** 2026-02-11
**Branch:** `feat/discrete-cell-assignment`
**Parent Branch:** `feat/cellpose-nuclei-prior`
**Status:** Design

## Overview

Replace continuous proportion estimation with discrete cell identity assignment using nuclei counts from Cellpose segmentation. Instead of estimating "40% macrophage, 40% fibroblast, 20% cancer," we assign each detected nucleus to a specific cell type, then derive proportions from those counts.

### Motivation

1. **Biological realism**: Cells are discrete entities, not continuous proportions
2. **Improved accuracy**: Discrete assignments constrain the solution space
3. **Better GEX deconvolution**: Fixed cell identities in Phase 2 eliminate joint optimization complexity
4. **Natural cellularity handling**: Spots with more nuclei contribute more signal (no artificial normalization)

## Architecture

```
┌─────────────────────────────────────────────────────────────────┐
│                        PHASE 1                                  │
│              Discrete Cell Assignment (ILP + EM)                │
├─────────────────────────────────────────────────────────────────┤
│  Inputs:                                                        │
│    - Protein signal X[i,m] (winsorized, NOT normalized)         │
│    - Nuclei counts N_i from Cellpose                            │
│    - Cell type profiles profile[t,m]                            │
│                                                                 │
│  Variables:                                                     │
│    - c[i,t] ∈ Z≥0  (integer count of type-t cells in spot i)   │
│    - β[m] ∈ R+     (per-marker scaling factor)                  │
│                                                                 │
│  Constraints:                                                   │
│    - Σₜ c[i,t] = N_i  (counts sum to nuclei count)             │
│                                                                 │
│  Objective:                                                     │
│    - min Σᵢ Σₘ (X[i,m] - Σₜ c[i,t] × profile[t,m] × β[m])²    │
│                                                                 │
│  Output: Locked integer assignments c[i,t]                      │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼ LOCK ASSIGNMENTS
┌─────────────────────────────────────────────────────────────────┐
│                        PHASE 2                                  │
│              GEX Deconvolution (Fixed Counts)                   │
├─────────────────────────────────────────────────────────────────┤
│  Inputs:                                                        │
│    - GEX data G[i,g]                                            │
│    - Fixed cell counts c[i,t] from Phase 1                      │
│                                                                 │
│  Variables:                                                     │
│    - E[i,t,g] (expression of gene g in type t at spot i)       │
│                                                                 │
│  Regularization:                                                │
│    - Global prior: E_global[t,g] from all spots                │
│    - Local prior: E_local[i,t,g] from neighborhood             │
│                                                                 │
│  Output: Deconvolved expression layers per cell type            │
└─────────────────────────────────────────────────────────────────┘
```

## Phase 1: Discrete Cell Assignment

### Mathematical Formulation

**Integer Quadratic Program (IQP):**

```
minimize    Σᵢ Σₘ (X[i,m] - Σₜ c[i,t] × profile[t,m] × β[m])²
subject to  Σₜ c[i,t] = N_i          ∀ spots i with N_i > 0
            c[i,t] ∈ Z≥0             ∀ i, t
            β[m] > 0                  ∀ markers m
```

**Key insight:** We don't need per-nucleus variables. Since nuclei within a spot are interchangeable, we only need to know HOW MANY of each type, not WHICH SPECIFIC nucleus is which.

For a spot with 5 nuclei and 6 cell types, we're partitioning 5 into 6 bins:
- Number of combinations: C(5+6-1, 6-1) = C(10,5) = 252
- Much more tractable than 6⁵ = 7,776 individual assignments

### EM Algorithm for β Optimization

The per-marker scaling factors β[m] are optimized via EM:

**E-step (given β):**
- Solve IQP to get optimal integer counts c[i,t]
- Gurobi handles integer quadratic programs

**M-step (given c):**
- Update β[m] via least squares regression:
```
β[m] = Σᵢ X[i,m] × pred[i,m] / Σᵢ pred[i,m]²
where pred[i,m] = Σₜ c[i,t] × profile[t,m]
```

**Convergence:** Iterate until β changes < ε or max iterations reached.

### Preprocessing Changes

**Remove per-spot normalization:**
- Current: CLR transform normalizes each spot to unit scale
- New: Raw counts (or log-transformed) preserve cellularity information
- A spot with 10 macrophages SHOULD have ~2x the CD68 signal of a spot with 5 macrophages

**Keep winsorization:**
- Clip extreme values (e.g., 1st/99th percentile) to handle outliers
- Prevents single aberrant measurements from dominating the fit

### Edge Cases

**Spots with N_i = 0 (no detected nuclei):**
- Exclude from Phase 1 optimization entirely
- These spots get no cell assignments
- Phase 2 can skip them or handle separately
- Rationale: If we can't see nuclei, we shouldn't pretend we know what's there

**Spots with very high nuclei counts:**
- IQP complexity scales with partition count, not exponentially
- C(N+T-1, T-1) is polynomial in N
- For N=50, T=10: ~10 billion combinations (may need heuristics)
- Consider: cap N_i at reasonable maximum (e.g., 30) or use relaxation + rounding

### No Spatial Smoothing in Phase 1

**Decision:** Spatial smoothing is NOT applied during discrete assignment.

**Rationale:**
1. Discrete counts can't smoothly interpolate between neighbors
2. Mixing continuous penalty with integer variables complicates the IQP
3. Spatial coherence should emerge naturally from biology
4. Keeps Phase 1 clean: pure signal-driven assignment

Spatial structure is preserved in Phase 2 via local neighborhood priors.

## Phase 2: GEX Deconvolution

### Fixed-Count Formulation

With locked cell counts c[i,t], Phase 2 becomes a constrained regression:

```
G[i,g] = Σₜ c[i,t] × E[i,t,g] + noise
```

We solve for E[i,t,g] (expression of gene g in cell type t at spot i).

### Global and Local Priors

**Global prior E_global[t,g]:**
- Average expression of gene g in cell type t across all spots
- Computed from current E estimates weighted by cell counts
- Prevents biologically implausible assignments (e.g., collagen → T cells)

**Local prior E_local[i,t,g]:**
- Average expression in the spatial neighborhood of spot i
- Captures local microenvironment effects
- Weighted by neighbor cell counts and distance

**Regularization:**
```
Loss = reconstruction_error + λ_global × ||E - E_global||² + λ_local × ||E - E_local||²
```

### Advantages Over Current Approach

1. **No joint optimization:** Proportions are fixed, only optimize expression
2. **Cleaner attribution:** Each cell contributes discretely to the GEX signal
3. **Better scaling:** More cells = more signal, naturally weighted
4. **Downstream compatibility:** Module 4 program discovery still works on E[i,t,g] layers

## Implementation Plan

### New/Modified Files

| File | Change |
|------|--------|
| `CITEgeist/model/discrete_assignment.py` | NEW: IQP solver, EM loop |
| `CITEgeist/model/gurobi_impl.py` | ADD: `solve_discrete_cell_counts()` |
| `CITEgeist/model/citegeist_model.py` | ADD: `run_discrete_cell_assignment()` |
| `CITEgeist/model/citegeist_model.py` | MODIFY: `run_cell_expression_pass1()` to accept fixed counts |
| `CITEgeist/model/preprocessing.py` | NEW: Non-normalized protein preprocessing |
| `tests/test_discrete_assignment.py` | NEW: Unit tests |

### API Changes

**New method on CitegeistModel:**
```python
def run_discrete_cell_assignment(
    self,
    nuclei_counts: Optional[pd.Series] = None,  # Or use precomputed from Cellpose
    max_em_iterations: int = 20,
    beta_convergence_tol: float = 1e-3,
    max_nuclei_cap: int = 30,  # Cap for computational tractability
) -> pd.DataFrame:
    """
    Phase 1: Assign discrete cell identities using IQP.

    Returns:
        DataFrame with columns for each cell type, values are integer counts.
        Index is spot names.
    """
```

**Modified Phase 2:**
```python
def run_cell_expression_pass1(
    self,
    cell_counts: Optional[pd.DataFrame] = None,  # NEW: Use fixed counts if provided
    # ... existing parameters ...
) -> AnnData:
    """
    If cell_counts provided, uses fixed discrete counts instead of proportions.
    """
```

### Computational Considerations

**IQP Complexity:**
- Gurobi handles IQP well for small-medium problems
- Per-spot: O(C(N+T-1, T-1)) combinations to consider
- Total: O(n_spots × combinations)
- For 1000 spots, 10 types, ~20 nuclei avg: tractable

**Potential optimizations:**
1. Warm-start IQP from previous EM iteration
2. Parallelize across spots (independent subproblems)
3. Use continuous relaxation + rounding for very high N_i
4. Cache profile × β products

## Testing Strategy

### Unit Tests

1. **IQP solver correctness:**
   - Known assignments → verify reconstruction
   - Edge cases: N=0, N=1, all same type

2. **EM convergence:**
   - Synthetic data with known β → verify recovery
   - Check monotonic loss decrease

3. **Preprocessing:**
   - Winsorization preserves scale relationships
   - No per-spot normalization applied

### Integration Tests

1. **End-to-end on simulated data:**
   - Known ground truth cell counts
   - Compare discrete assignments to truth

2. **Xenium benchmark:**
   - Compare discrete vs continuous on existing benchmark
   - Metrics: JSD, correlation with ground truth

## Success Criteria

1. **Proportion accuracy:** JSD improvement over continuous baseline on Xenium benchmark
2. **GEX deconvolution:** Correlation improvement for marker genes
3. **Computational:** Phase 1 completes in < 2x time of continuous approach
4. **Biological plausibility:** No impossible assignments (verified manually on real data)

## Open Questions

1. **High nuclei count spots:** At what N_i does IQP become impractical? Need benchmarking.
2. **Nuclei count uncertainty:** Cellpose isn't perfect. Should we model count uncertainty?
3. **Mixed-type nuclei:** Some cells may express markers of multiple types. Handle how?
4. **Phase 2 spots with 0 counts:** Impute from neighbors? Skip? Flag for review?

## References

- Parent feature: `feat/cellpose-nuclei-prior` (Cellpose segmentation infrastructure)
- Current Module 3: `CITEgeist/model/gurobi_impl.py`
- Gurobi IQP documentation: https://www.gurobi.com/documentation/
