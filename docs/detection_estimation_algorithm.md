# Detection-Estimation Algorithm

This document describes the two-stage detection + estimation model for cell type proportion estimation from spatial protein markers.

## Overview

The algorithm replaces continuous proportion optimization with discrete cell assignment:

1. **Stage 1: Detection** - Identify which cell types are present in each spot
2. **Stage 2: Estimation** - Assign integer cell counts given the detection mask

This enforces true zeros for non-detected cell types, improving accuracy over continuous relaxation.

## Stage 1: Adaptive GMM Detection

**File**: `CITEgeist/model/detection.py`

For each cell type k with markers M_k:

```
1. Extract marker data X[:, M_k]
2. Log-transform: X_log = log1p(X)
3. Fit 2-component GMM (background vs signal)
4. Identify signal cluster (higher mean across markers)
5. Compute adaptive threshold based on cluster weights:
   - w < 0.3 (rare type)   → threshold = 0.3
   - w < 0.5 (uncommon)    → threshold = 0.4
   - w < 0.7 (moderate)    → threshold = 0.5 (base)
   - w ≥ 0.7 (common type) → threshold = 0.6
6. detected[i,k] = (posterior > threshold)
```

### Key Design Decisions

**Log-transform before GMM**: Raw antibody data has heavy right tails. Log-transform compresses dynamic range, helping GMM find the true signal/background boundary rather than splitting within the signal population.

**Adaptive threshold**: Rare cell types (small signal cluster) get lower thresholds for sensitivity. Common cell types (large signal cluster) get higher thresholds to avoid over-detection (GMM may be splitting moderate vs high expression rather than signal vs background).

## Stage 2: Masked IQP Estimation

**File**: `CITEgeist/model/masked_iqp.py`

### Model

For spot i with markers m:
```
X[i,m] = alpha[m] + sum_k counts[i,k] * profile[k,m] * beta[m] + noise
```

Where:
- `alpha[m]` = baseline (background signal per marker)
- `beta[m]` = signal-per-cell per marker
- `counts[i,k]` = integer cell count for type k in spot i
- `profile[k,m]` = 1 if marker m defines cell type k, else 0

### Constraints

1. `counts[i,k] = 0` if `detected[i,k] = False`
2. `counts[i,k] in {0, 1, ..., N_i}` if `detected[i,k] = True`
3. `sum_k counts[i,k] = N_i` (nuclei sum constraint)
4. `alpha[m] >= 0`
5. `beta[m] >= beta_floor[m]` (floor constraint - see below)

### Key Innovation: Direct Beta Estimation

The chicken-and-egg problem: OLS learns beta from IQP counts, but IQP needs good beta to allocate correctly. If beta starts wrong, IQP makes bad allocations, and OLS learns bad beta, causing collapse.

**Solution**: Estimate beta directly from single-type spots.

For spots where ONLY cell type k is detected, all nuclei belong to that type:
```
beta_m = (X_m - alpha_m) / nuclei_count
```

This provides unbiased estimates without needing to solve the allocation problem first.

### Block Coordinate Descent Algorithm

```
1. DIRECT BETA INITIALIZATION
   For each marker m used by cell type k:
     - Find spots where ONLY type k is detected
     - alpha_m = median(X_m in non-detected spots)
     - beta_m = median((X_m - alpha_m) / nuclei_count in single-type spots)

   SAVE beta_floor = beta.copy()  # This is the FLOOR

2. FOR each BCD iteration:

   a. FIX BETA → SOLVE IQP FOR COUNTS
      For each spot i with detected types T_i:
        - If |T_i| = 0: skip
        - If |T_i| = 1: counts[i, T_i[0]] = nuclei_count[i]
        - If |T_i| > 1: solve small IQP (only ~7 integer variables)

      Update alpha from residuals: alpha = max(0, median(X - predicted))

   b. FIX COUNTS → SOLVE OLS FOR ALPHA, BETA
      For each marker m:
        - Weighted least squares: X[:,m] = alpha_m + beta_m * effective_counts
        - beta_ols = OLS estimate

      FLOOR CONSTRAINT: beta = max(beta_ols, beta_floor)

   c. CHECK CONVERGENCE
      If max|beta_new - beta_old| / max|beta_old| < 0.01: STOP
```

### The Floor Constraint

**Key robustness mechanism**: Beta can only go UP from the direct estimate, never collapse below.

Why this works:
- Direct estimate is unbiased (from single-type spots)
- OLS may push beta down due to allocation errors (mixing gives lower per-cell signal)
- Floor prevents runaway collapse without hyperparameter tuning
- No ridge regularization needed - the floor is the regularization

### Per-Spot IQP Solver

For spots with >1 detected types, we solve a small IQP:

```
minimize sum_m weights[m] * (X[i,m] - alpha[m] - sum_k c[k] * signal_coeff[k,m])^2
subject to:
  c[k] in {0, 1, ..., N_i} for k in detected_types
  sum_k c[k] = N_i
```

For very small problems (2 types, ≤20 nuclei), we enumerate all valid splits.

## Key Files

| File | Purpose |
|------|---------|
| `detection.py` | Stage 1: Adaptive GMM detection |
| `masked_iqp.py` | Stage 2: Masked IQP with direct beta |
| `detection_estimation.py` | Orchestrator combining both stages |

## Critical Bug Fixes Applied

| Bug | Symptom | Root Cause | Fix |
|-----|---------|------------|-----|
| Beta collapse | CD68 beta → 0.001 | OLS learned from bad allocations | Floor constraint: beta ≥ direct_estimate |
| Sigma explosion | weights → 1M | Near-zero MAD from good fits | Use uniform weights (skip sigma EM) |
| Beta reset | Lost learning per EM iter | Re-initialized beta each iteration | Pass beta_init, now simplified with floor |

## Results

### Xenium Region 0 (1407 spots, 7 cell types)

| Cell Type | Pearson r | Detection Rate |
|-----------|-----------|----------------|
| B cells | 0.95 | 9.7% |
| CD8+ T cells | 0.83 | 79.6% |
| Macrophages | 0.77 | 83.4% |
| Endothelial | 0.66 | 71.3% |
| Epithelial | 0.53 | 54.8% |
| Fibroblasts | 0.49 | 70.9% |
| CD4+ T cells | 0.44 | 58.0% |
| **Overall** | **0.68** | - |

### Before vs After Floor Constraint Fix

| Metric | Before (broken) | After (fixed) |
|--------|-----------------|---------------|
| Overall r | 0.24 | 0.68 |
| Macrophage r | 0.28 | 0.77 |
| CD68 beta | 0.001 | 19.3 |

## Comparison to Other Approaches

| Method | Prop r | Notes |
|--------|--------|-------|
| CITEgeist Hybrid | 0.726 | Continuous → discretize |
| **Detection-Estimation** | **0.68** | Direct IQP |
| CITEgeist Continuous | 0.725 | No integer constraint |
| Discrete IQP + CLR | 0.63 | No detection stage |
| Discrete IQP only | 0.43 | No detection, bad preprocessing |

## When to Use

**Detection-Estimation is preferred when:**
- You need explicit cell counts (not just proportions)
- Cell types have varying prevalences (adaptive threshold helps)
- You want to avoid continuous-to-discrete artifacts

**CITEgeist Hybrid is preferred when:**
- You need maximum accuracy (0.726 vs 0.68)
- Continuous proportions are acceptable
- You can discretize post-hoc if needed

## Implementation Notes

### Computational Cost

Per-spot IQP is much faster than global IQP:
- Global IQP: O(n_spots × n_types) integer variables
- Per-spot IQP: O(n_types) integer variables per solve

For 1407 spots × 7 types, per-spot takes ~2 seconds total vs ~300 seconds for global.

### Memory Usage

Negligible - all operations are spot-wise or marker-wise, no large matrices.

### Dependencies

- `gurobipy` for IQP solver
- `sklearn.mixture.GaussianMixture` for GMM detection

---

*Last updated: 2026-02-27*
