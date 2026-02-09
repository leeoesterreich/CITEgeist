# Fibroblast & Endothelial Benchmark Fix Design

**Date:** 2026-02-09
**Problem:** Poor Xenium benchmark performance for Fibroblasts (R=0.42) and Endothelial (R=0.58) due to single-marker profiles being underestimated
**Goal:** Implement two independent approaches (A': asymmetric loss, C: adaptive Laplacian) plus combined A'+C

---

## Problem Analysis

| Cell Type | Current Markers | Pearson R | Issue |
|-----------|-----------------|-----------|-------|
| Fibroblasts | alphaSMA only | 0.42 | Underestimated where should be high |
| Endothelial | CD31 only | 0.58 | Underestimated where should be high |
| T cells, Macrophages | 2-3 markers | 0.60-0.80 | Good performance |

**Root cause:** Single-marker profiles lose the "competition" for proportion budget against multi-marker profiles that have stronger, more consistent signal. The optimization allocates proportions to cell types with redundant marker signal first, underestimating single-marker profiles.

**Key constraint:** Beta normalization must be preserved. The per-marker beta/alpha EM loop normalizes for absolute signal differences (e.g., CD68 at 20,000 vs alphaSMA at 500). Any fix must work within the beta-normalized residual framework, not raw signal values.

---

## Approach A': Asymmetric Marker-Count Loss

### Concept

Penalize underestimation more heavily for cell types with fewer markers. The asymmetric boost scales continuously with marker count - no special-case code for single markers.

### Mathematical Formulation

```python
# Reference: max marker count across all cell types
max_markers = max(len(profile["Major"]) for profile in cell_profile_dict.values())

# Per-cell-type boost factor (scales inversely with marker count)
n_markers_j = len(cell_profile_dict[j]["Major"])
underestimation_boost = (max_markers / n_markers_j) ** lambda_coverage

# In the loss function (works with beta-normalized residuals)
residual = S[i,m] - alpha[m] - beta[m] * Y[i,j]

if residual > 0:  # Underestimating - marker signal not explained by proportion
    loss += base_weight * underestimation_boost * residual^2
else:  # Overestimating - normal loss
    loss += base_weight * residual^2
```

### Example Values

With `lambda_coverage = 1.0` (linear inverse scaling):

| Cell Type | Markers | Boost Factor | Effect |
|-----------|---------|--------------|--------|
| Macrophages | 3 | 1.0x | No change |
| CD4+ T cells | 3 | 1.0x | No change |
| CD8+ T cells | 3 | 1.0x | No change |
| B cells | 2 | 1.5x | Mild boost |
| Fibroblasts | 1 | 3.0x | Strong boost |
| Endothelial | 1 | 3.0x | Strong boost |
| Epithelial | 1 | 3.0x | Strong boost |

### Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `lambda_coverage` | 1.0 | Exponent for marker-count scaling. 0=symmetric, 1=linear, 2=aggressive |

### Implementation Location

Modify `optimize_cell_proportions_per_marker()` in `gurobi_impl.py` (lines 645-703, the E-step loss construction).

---

## Approach C: Adaptive Spatial Regularization

### Concept

Single-marker profiles are noisier, so lean harder on spatial coherence. Fibroblasts (stroma) and endothelial cells (vessels) form spatially coherent domains, making stronger spatial smoothing a reasonable prior.

### Mathematical Formulation

```python
# Reference: max marker count across all cell types
max_markers = max(len(profile["Major"]) for profile in cell_profile_dict.values())

# Per-cell-type lambda (scales inversely with marker count)
n_markers_j = len(cell_profile_dict[j]["Major"])
lambda_laplacian_j = lambda_base * (max_markers / n_markers_j) ** laplacian_scale

# In the Laplacian term (per cell type)
laplacian_term = sum_j lambda_laplacian_j * sum_{(i,k) in edges} L[i,k] * Y[i,j] * Y[k,j]
```

### Example Values

With `lambda_base = 0.1`, `laplacian_scale = 1.0`:

| Cell Type | Markers | Effective Lambda |
|-----------|---------|------------------|
| Macrophages | 3 | 0.10 |
| CD4+ T cells | 3 | 0.10 |
| CD8+ T cells | 3 | 0.10 |
| B cells | 2 | 0.15 |
| Fibroblasts | 1 | 0.30 |
| Endothelial | 1 | 0.30 |
| Epithelial | 1 | 0.30 |

### Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `adaptive_laplacian` | False | Enable per-celltype lambda scaling |
| `laplacian_scale` | 1.0 | Exponent for marker-count scaling |

### Implementation Location

Modify `optimize_cell_proportions_per_marker()` in `gurobi_impl.py` (lines 442-455, Laplacian term construction).

---

## Approach A'+C: Combined

Both fixes use the same `(max_markers / n_markers_j)` scaling formula, so they're naturally compatible. No additional code beyond A' and C - just enable both flags.

Combined effect for single-marker profiles:
- 3x boost on underestimation penalty (A')
- 3x stronger spatial smoothing (C)

---

## Git Worktree Strategy

Three isolated worktrees branching from `dev`:

```
CITEgeist/                              # Main worktree (dev) - untouched
../CITEgeist-asymmetric-loss/           # Approach A'
../CITEgeist-adaptive-laplacian/        # Approach C
../CITEgeist-combined-fix/              # A' + C merged
```

### Branch Names

| Worktree | Branch | Approach |
|----------|--------|----------|
| CITEgeist-asymmetric-loss | `feature/asymmetric-marker-loss` | A' only |
| CITEgeist-adaptive-laplacian | `feature/adaptive-laplacian` | C only |
| CITEgeist-combined-fix | `feature/marker-count-normalization` | A' + C |

---

## Implementation Tasks

### Worktree 1: Approach A' (Asymmetric Loss)

1. [ ] Add `lambda_coverage` parameter to `optimize_cell_proportions_per_marker()`
2. [ ] Compute `max_markers` and per-celltype `underestimation_boost`
3. [ ] Modify loss construction to apply asymmetric weighting when residual > 0
4. [ ] Pass parameter through `run_cell_proportion_model()`
5. [ ] Add unit test for asymmetric loss behavior
6. [ ] Update `benchmark_constants.py` or benchmark script to enable feature
7. [ ] Run Xenium benchmark, record Fibroblast/Endothelial R

### Worktree 2: Approach C (Adaptive Laplacian)

1. [ ] Add `adaptive_laplacian` and `laplacian_scale` parameters
2. [ ] Compute per-celltype `lambda_laplacian_j` based on marker count
3. [ ] Modify Laplacian term to use per-celltype lambda vector
4. [ ] Pass parameters through `run_cell_proportion_model()`
5. [ ] Add unit test for adaptive lambda calculation
6. [ ] Run Xenium benchmark, record Fibroblast/Endothelial R

### Worktree 3: Approach A'+C (Combined)

1. [ ] Create branch from dev
2. [ ] Cherry-pick or merge A' changes
3. [ ] Cherry-pick or merge C changes
4. [ ] Resolve any conflicts (should be minimal - different code sections)
5. [ ] Run Xenium benchmark with both features enabled
6. [ ] Compare all three approaches

---

## Evaluation Plan

Run Xenium benchmark for each approach on all 5 regions:

```bash
# In each worktree
cd Benchmarking/xenium_benchmarking/CITEgeist
sbatch src/run_benchmark.sh  # with appropriate flags

# Then evaluate
cd ../evaluation
python evaluate_pipeline_stages.py
```

### Success Metrics

| Cell Type | Current R | Target R | Stretch Goal |
|-----------|-----------|----------|--------------|
| Fibroblasts | 0.42 | 0.55 | 0.65 |
| Endothelial | 0.58 | 0.65 | 0.70 |

### Comparison Matrix

| Approach | Fibroblast R | Endothelial R | Overall R | Notes |
|----------|--------------|---------------|-----------|-------|
| Baseline | 0.42 | 0.58 | 0.60 | Current |
| A' only | ? | ? | ? | |
| C only | ? | ? | ? | |
| A' + C | ? | ? | ? | |

---

## Risks and Mitigations

| Risk | Mitigation |
|------|------------|
| Asymmetric loss destabilizes EM convergence | Start with low `lambda_coverage` (0.5), increase gradually |
| Adaptive lambda over-smooths rare populations | Cap effective lambda at 3x base (laplacian_scale ≤ 1.0) |
| Combined approach over-corrects → overestimation | Test independently first, tune parameters |
| Gurobi solver issues with asymmetric loss | Loss is still convex (squared residuals), should be fine |

---

## Parameter Tuning Strategy

If initial results are underwhelming:

1. **A' tuning:** Try `lambda_coverage` in [0.5, 1.0, 1.5, 2.0]
2. **C tuning:** Try `laplacian_scale` in [0.5, 1.0, 1.5]
3. **Combined tuning:** Grid search over both parameters

---

## Files to Modify

| File | A' Changes | C Changes |
|------|------------|-----------|
| `gurobi_impl.py` | Asymmetric loss in E-step | Per-celltype Laplacian term |
| `citegeist_model.py` | Add `lambda_coverage` param | Add `adaptive_laplacian`, `laplacian_scale` params |
| Benchmark scripts | Enable new features | Enable new features |

---

**Created:** 2026-02-09
**Status:** Ready for implementation
