# CITEgeist Single-Cell Output Mode Design

**Date**: 2026-02-17
**Status**: Draft
**Authors**: Alexander Chang, Claude

## Overview

Transform CITEgeist output from spot-level layers to true single-cell AnnData where each observation represents one nucleus. This requires first fixing the discrete IQP performance regression, then implementing the single-cell assembly pipeline.

## Problem Statement

1. **Current output is confusing**: `expand_prop_gex_adata()` creates pseudo-single-cells by duplicating spots — it's a workaround, not true single-cell
2. **Layer-based output is hard to use**: `{celltype}_genes_pass1` layers don't integrate with standard scanpy workflows
3. **Discrete IQP underperforms continuous**: Proportion correlation dropped from ~0.95 (continuous) to ~0.78 (discrete) — must fix before building on top

## Goals

1. Fix discrete IQP to recover ~0.95 proportion correlation
2. Generate true single-cell AnnData (each nucleus = one observation)
3. Use nucleus-scaled GEX optimization with GEX-informed refinement
4. Simplify API to one method call
5. Deprecate legacy layer-based output

## Non-Goals (This Phase)

- Morphology-informed cell type assignment (future work)
- Position-informed cell type assignment (future work)
- Removing legacy continuous mode (kept for backwards compatibility)

## Success Criteria

- Proportion correlation ≥ 0.90 on high_seg, mixed, and xenium benchmarks
- GEX metrics improved or maintained
- Clean single-cell AnnData output compatible with scanpy/squidpy

---

## Phase 0: Fix Discrete IQP (Prerequisite)

### Problem Analysis

| Dataset | Continuous Mode | Discrete Mode | Gap |
|---------|-----------------|---------------|-----|
| high_seg | 0.948 | 0.78 | -18% |
| mixed | 0.869 | 0.44 | -49% |

**Root causes identified:**
1. No sparsity regularization in IQP (continuous has it)
2. 2.6% of spots get zero cells assigned (proportion = 0)
3. Integer rounding loses precision for low-nuclei spots

### Proposed Fix

Add sparsity penalty to IQP objective:

```
Current:  min Σ (X[i,m] - α[m] - Σₜ c[i,t] × profile[t,m] × β[m])²

Proposed: min Σ (X[i,m] - α[m] - Σₜ c[i,t] × profile[t,m] × β[m])²
              + λ_sparse × Σₜ 1{c[i,t] > 0}  (count non-zero types)
```

The sparsity term encourages fewer cell types per spot, matching the biological reality that most spots have 2-4 dominant cell types.

### Implementation

1. Modify `solve_discrete_cell_counts()` in `gurobi_impl.py`
2. Add `lambda_sparse` parameter (default: 0.1)
3. Use indicator constraints in Gurobi: `y[t] = 1 if c[t] > 0`
4. Add penalty term: `+ λ_sparse × Σ y[t]`

### Validation

- Re-run benchmarks on high_seg, mixed (5 replicates each)
- Target: proportion correlation ≥ 0.90
- If not achieved, investigate additional fixes (preprocessing, beta estimation)

---

## Phase 1: Initial Cell Type Count Estimation (IQP)

**No changes** — use existing IQP/EM algorithm (after Phase 0 fix)

- Input: Antibody data, cell profiles, nuclei counts per spot
- Output: Integer cell counts per type per spot
- Method: `solve_discrete_cell_counts()` with sparsity penalty

---

## Phase 2: Nucleus-Scaled GEX Optimization

**Mathematical Formulation:**

```
For each spot i with N_i nuclei:
    minimize    ||Y_i - Σ_t (n_{i,t} × g_t)||² + λ × ||g_t||²
    subject to  g_t ≥ 0  (non-negative expression)

Where:
    Y_i     = observed gene expression (1 × M genes)
    n_{i,t} = cell count for type t from Phase 1
    g_t     = per-cell expression vector for type t (decision variable)
```

**Key change from current**: Uses integer counts `n_t` instead of proportions `p_t`, so expression magnitude scales with cellularity.

**Output**: Per-cell-type expression vectors + reconstruction error per spot

---

## Phase 3: GEX-Informed Count Refinement

**Mathematical Formulation:**

```
For each spot i:
    minimize    α × ||Y_i - Σ_t (n_t × g_t^*)||² + (1-α) × ||n - n^{IQP}||²
    subject to  Σ_t n_t = N_i        (nuclei constraint from Cellpose)
                n_t ∈ Z≥0           (non-negative integers)
                |n_t - n^{IQP}_t| ≤ δ  (optional: bounded deviation)

Where:
    g_t^*      = expression vectors from Phase 2 (fixed)
    n^{IQP}    = original counts from Phase 1
    α          = GEX trust weight (tunable, default 0.5)
    δ          = max deviation from antibody counts (optional)
```

**Purpose**: Allow GEX signal to correct antibody-based estimates when there's strong evidence.

**Convergence**: Iterate Phases 2-3 until Δcounts < threshold or max_iterations reached (typically 2-5 iterations).

---

## Phase 4: Single-Cell Assembly

**Input:**
- Final cell counts from Phase 3
- Final per-cell expressions from Phase 2
- Cellpose nuclei centroids

**Method:**
1. For each spot, get actual nucleus centroid positions from Cellpose
2. Randomly assign cell types to nuclei (respecting counts from Phase 3)
3. Assign expression vector based on cell type (equal split within type)

**Output**: Single-cell AnnData with structure:

```python
adata.X                    # Gene expression (n_cells × n_genes)
adata.obs['barcode']       # Parent spot barcode
adata.obs['cell_type']     # Assigned cell type
adata.obs['parent_spot']   # Parent spot ID
adata.obs['nuclei_in_spot'] # Total nuclei in parent spot
adata.obs['cell_index_in_spot'] # 0-indexed position within spot
adata.obsm['spatial']      # Actual nucleus centroid from Cellpose
```

**Expression distribution**: Equal split within type. If spot has 3 Cancer cells, each gets (deconvolved Cancer GEX) / 3.

---

## API Design

### New Method

```python
def run_single_cell_deconvolution(
    self,
    nuclei_counts: pd.Series,
    max_refinement_iterations: int = 5,
    gex_trust_weight: float = 0.5,
    lambda_sparse: float = 0.1,
    output_formats: List[str] = ["h5ad", "parquet", "csv"],
) -> ad.AnnData:
    """
    Run complete single-cell deconvolution pipeline.

    Combines: IQP cell assignment + nucleus-scaled GEX + refinement + assembly

    Returns:
        Single-cell AnnData with one observation per nucleus
    """
```

### Deprecations

- `expand_prop_gex_adata()` — raise DeprecationWarning
- `run_discrete_cell_assignment()` — raise DeprecationWarning
- `run_cell_expression_pass1()` — raise DeprecationWarning

All warnings point to `run_single_cell_deconvolution()`.

---

## Output Files

| File | Format | Purpose |
|------|--------|---------|
| `{sample}_single_cell.h5ad` | AnnData | Primary, scanpy-compatible |
| `{sample}_single_cell.parquet` | Parquet | Interoperability (pandas, R) |
| `{sample}_single_cell_metadata.csv` | CSV | Easy inspection of cell-level metadata |

---

## Legacy Mode

Keep continuous proportions mode for backwards compatibility:
- Existing papers used it
- Datasets without histology images can fall back to it
- No removal planned, just soft-deprecated in docs

---

## Benchmarking Plan

### Datasets
1. **high_seg** — 5 replicates, high segmentation quality
2. **mixed** — 5 replicates, mixed cell types per spot
3. **xenium** — Xenium pseudo-Visium regions

### Metrics
- Proportion correlation (Pearson r)
- Proportion RMSE
- GEX RMSE, NRMSE, MAE
- JSD (Jensen-Shannon divergence)

### Process
1. Implement Phase 0 (sparsity fix)
2. Benchmark on all datasets
3. If proportion correlation < 0.90, investigate further
4. Once Phase 0 passes, implement Phases 1-4
5. Final benchmark with full pipeline

---

## Implementation Timeline

| Phase | Description | Estimated Effort |
|-------|-------------|------------------|
| 0 | Fix discrete IQP (sparsity) | 1 week |
| 0.5 | Benchmark & validate | 2-3 days |
| 1-4 | Single-cell pipeline | 2-3 weeks |
| Testing | Integration tests, edge cases | 1 week |
| **Total** | | **5-6 weeks** |

---

## Open Questions

1. Should we add position-informed cell type assignment in a future phase?
2. What morphology features would be most informative (size, eccentricity)?
3. Should the refinement loop have an early stopping criterion based on GEX improvement?

---

## Appendix: Architecture Diagram

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                    run_single_cell_deconvolution()                          │
└─────────────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│  PHASE 1: Cell Type Count Estimation (IQP + Sparsity)                       │
│  Input:  Antibody data, cell profiles, nuclei counts per spot               │
│  Output: Integer cell counts per type per spot                              │
└─────────────────────────────────────────────────────────────────────────────┘
                                    │
                          ┌─────────┴─────────┐
                          ▼                   │
┌─────────────────────────────────────────────────────────────────────────────┐
│  PHASE 2: Nucleus-Scaled GEX Optimization                                   │
│  Input:  GEX data, current cell counts                                      │
│  Output: Per-cell-type expression vectors + reconstruction error            │
└─────────────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│  PHASE 3: GEX-Informed Count Refinement                                     │
│  Input:  GEX error, current counts, nuclei constraint                       │
│  Output: Refined cell counts                                                │
│  Convergence: Δcounts < threshold OR max_iterations                         │
└─────────────────────────────────────────────────────────────────────────────┘
                                    │
                          ┌─────────┘
                          │ (iterate until convergence)
                          ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│  PHASE 4: Single-Cell Assembly                                              │
│  Input:  Final counts, final expressions, Cellpose centroids                │
│  Output: Single-cell AnnData (h5ad + parquet + csv)                         │
└─────────────────────────────────────────────────────────────────────────────┘
```
