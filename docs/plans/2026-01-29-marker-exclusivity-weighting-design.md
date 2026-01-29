# Design: Reconstruction Exclusivity Weighting for Marker Quality

**Date**: 2026-01-29
**Status**: Approved
**Branch**: `hierarchical_approach`

## Problem

Per-marker beta optimization (Bug 1 + Bug 2 fixes) improved 6/7 cell types but regressed Fibroblasts (r=0.360 → 0.156). Root cause: VIM (Vimentin) is ubiquitously expressed across many cell types, making it a poor discriminator for Fibroblasts. The loss normalization fix (Bug 1) removed an accidental 2x weight advantage that previously compensated for VIM's non-specificity.

The optimizer treats all markers equally in the loss function. A marker like CD68 (specific to Macrophages) gets the same weight as VIM (expressed everywhere). This inflates Fibroblast predictions at spots where VIM signal is present but Fibroblasts are not.

## Approach: Reconstruction Exclusivity

After the global EM pass converges, compute how exclusively each marker correlates with its assigned cell type versus other cell types. Use this score to weight the finetuning pass loss function.

### Exclusivity Score

For each marker m assigned to cell type j:

- **Owner correlation**: `r_owner = corr(S[:,m], Y[:,j])` — how well the marker tracks its assigned cell type's proportion
- **Best non-owner correlation**: `r_best_other = max over k != owners(m) of corr(S[:,m], Y[:,k])` — how well the marker tracks the best alternative cell type

```
exclusivity[m] = max(r_owner, 0) / (max(r_owner, 0) + max(r_best_other, 0) + epsilon)
```

**Examples (expected):**
- CD68 (Macrophages): r_owner=0.8, r_best_other=0.1 → exclusivity ≈ 0.89
- VIM (Fibroblasts): r_owner=0.3, r_best_other=0.25 → exclusivity ≈ 0.55
- EPCAM (Epithelial): r_owner=0.7, r_best_other=0.15 → exclusivity ≈ 0.82

### Pipeline Placement

```
Global EM (unweighted) → Y_global, beta
    → compute_marker_exclusivity(S, Y_global, assignment_matrix)
        → Finetuning (exclusivity-weighted) → Y_final
```

The global pass runs unweighted to produce an honest first estimate. Exclusivity is computed from those estimates, then applied only during finetuning. This avoids circular reasoning (exclusivity depending on proportions that depend on exclusivity).

### Loss Weight Modification

In `deconvolute_local_cell_proportions_per_marker()`, the error term weight changes from:

```python
weight = 1.0 / (n_owners * markers_per_celltype[j])
```

to:

```python
weight = exclusivity[m] / (n_owners * markers_per_celltype[j])
```

### Edge Cases and Guardrails

1. **Floor at 0.3**: Even the worst marker retains 30% weight. Prevents collapsing cell types with only one weak marker (the Fibroblast scenario — we want to downweight VIM, not remove it).

2. **Shared markers**: For markers like CD3E (CD4+ T and CD8+ T), `r_owner` is the correlation with the combined owner signal (sum of owner proportions). `r_best_other` excludes all owner cell types.

3. **Near-zero variance markers**: Epsilon in the denominator prevents division by zero. Floor ensures these markers don't vanish.

4. **Unknown cell type**: Unaffected directly (no markers). Indirectly affected via proportion redistribution.

5. **Backward compatibility**: `marker_exclusivity=None` default preserves current behavior (uniform weight).

## Implementation Scope

### New Function

**`compute_marker_exclusivity()`** in `gurobi_impl.py`
- Input: `marker_level_data (N,M)`, `Y_values (N,T)`, `marker_owners`, `assignment_matrix (M,T)`
- Output: `exclusivity (M,)` array, values in [0.3, 1.0]
- ~30 lines, pure numpy

### Modified Functions

**`finetune_cell_proportions_per_marker()`** in `gurobi_impl.py`
- New parameter: `marker_exclusivity: Optional[np.ndarray] = None`
- Passes through to `deconvolute_local_cell_proportions_per_marker()`

**`deconvolute_local_cell_proportions_per_marker()`** in `gurobi_impl.py`
- New parameter: `marker_exclusivity: Optional[np.ndarray] = None`
- Multiplies loss weight by `exclusivity[m]` when constructing error terms
- Falls back to 1.0 when `None`

**`run_cell_proportion_model()`** in `citegeist_model.py`
- After global EM returns, calls `compute_marker_exclusivity()`
- Logs exclusivity scores
- Passes exclusivity to `finetune_cell_proportions_per_marker()`

### Test

- New test verifying exclusivity computation for known marker patterns
- Specific markers → high scores, ubiquitous markers → low scores
- Backward compatibility: `None` exclusivity produces identical results to current code

### Files Changed

| File | Change |
|------|--------|
| `CITEgeist/model/gurobi_impl.py` | New function + 2 modified functions |
| `CITEgeist/model/citegeist_model.py` | Orchestration in `run_cell_proportion_model()` |
| `tests/test_marker_exclusivity.py` | New test file |

## Success Criteria

- Fibroblast correlation improves from 0.156 toward pre-fix baseline (0.360)
- Other 6 cell types maintain their post-fix improvements
- Overall mean r improves from 0.382 toward or above pre-fix 0.412
