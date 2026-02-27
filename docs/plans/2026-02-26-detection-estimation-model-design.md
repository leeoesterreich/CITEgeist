# Design: Two-Stage Detection + Estimation Model for Cell Type Deconvolution

**Date**: 2026-02-26
**Status**: Approved
**Author**: Claude + Alex

## Problem Statement

Current proportion estimation over-allocates rare cell types due to:
1. L2 regularization pushing toward uniform distributions
2. No distinction between background noise and true signal
3. Continuous model can't represent true zeros

**Evidence (Region 1 B cells)**:
- Ground truth: 0.3% of spots have B cells
- Continuous model: 59.6% of spots have B cells (>1%)
- Discrete IQP: 42.9% of spots have B cells
- Overestimate: 143-200x

**Root cause**: Weak marker signal (background noise) gets interpreted as "low but nonzero" proportion instead of "absent."

## Solution Overview

Two-stage model:
1. **Detection**: Binary classification per cell type per spot (present/absent)
2. **Estimation**: Global IQP for cell counts only among detected types, with learned noise model

## Stage 1: Multivariate GMM Detection

### Rationale

Cell types are defined by multiple markers. A multivariate GMM in the joint marker space:
- Captures covariance (e.g., CD4+ T cells have BOTH CD3 AND CD4 elevated)
- Distinguishes true signal from single-marker noise
- Handles markers with different scales

### Algorithm

```python
def detect_cell_types(X, marker_groups, threshold=0.5):
    """
    Binary detection per cell type using multivariate GMM.

    Args:
        X: (n_spots, n_markers) antibody signal matrix
        marker_groups: dict mapping cell_type_name -> list of marker indices
        threshold: posterior probability threshold for detection

    Returns:
        detected: (n_spots, n_types) boolean mask
    """
    n_spots = X.shape[0]
    n_types = len(marker_groups)
    detected = np.zeros((n_spots, n_types), dtype=bool)

    for k, (cell_type, marker_indices) in enumerate(marker_groups.items()):
        # Extract markers for this cell type
        marker_data = X[:, marker_indices]  # (n_spots, n_markers_k)

        # Fit 2-component GMM (background vs signal)
        gmm = GaussianMixture(
            n_components=2,
            covariance_type='full',  # capture marker correlations
            random_state=42,
        )
        gmm.fit(marker_data)

        # Identify signal cluster (higher mean across markers)
        cluster_means = gmm.means_.sum(axis=1)
        signal_cluster = np.argmax(cluster_means)

        # Get posterior probability of signal
        posteriors = gmm.predict_proba(marker_data)[:, signal_cluster]

        # Binary detection
        detected[:, k] = posteriors > threshold

    return detected
```

### Output

- `detected[i, k] = 1` if cell type k is present in spot i
- `detected[i, k] = 0` if cell type k is absent (HARD ZERO enforced in Stage 2)

## Stage 2: Global IQP with Learned Noise Variance

### Observation Model

```
X[i,m] = α[m] + Σ_k (c[i,k] × profile[k,m] × β[m]) + ε[i,m]

Where:
  X[i,m]      = observed antibody signal for spot i, marker m
  α[m]        = baseline/background for marker m (global, learned)
  c[i,k]      = cell count of type k in spot i (integer, ≥0)
  profile[k,m]= 1 if marker m defines type k, 0 otherwise (fixed)
  β[m]        = signal-per-cell for marker m (global, learned)
  ε[i,m]      ~ N(0, σ²[m]) heteroskedastic noise
```

### Optimization Problem

```
minimize  Σ_i Σ_m (1/σ²[m]) × (X[i,m] - α[m] - Σ_k c[i,k]·profile[k,m]·β[m])²

subject to:
  c[i,k] = 0              if detected[i,k] = 0     # DETECTION MASK
  c[i,k] ∈ {0,1,...,N_i}  if detected[i,k] = 1     # integer counts
  Σ_k c[i,k] = N_i                                  # sum to nuclei count
  α[m] ≥ 0                                          # non-negative baseline
  β[m] ≥ ε_min                                      # positive signal-per-cell
```

### EM Algorithm for Joint Estimation

```python
def solve_detection_estimation(X, nuclei_counts, profile, marker_groups, max_iter=10):
    """
    Two-stage detection + estimation with learned noise variance.

    Args:
        X: (n_spots, n_markers) antibody signal
        nuclei_counts: (n_spots,) from Cellpose segmentation
        profile: (n_types, n_markers) binary assignment matrix
        marker_groups: dict mapping cell_type -> marker indices
        max_iter: maximum EM iterations

    Returns:
        detected: (n_spots, n_types) binary presence mask
        counts: (n_spots, n_types) integer cell counts
        alpha: (n_markers,) learned baselines
        beta: (n_markers,) learned signal-per-cell
        sigma_sq: (n_markers,) learned noise variances
    """
    n_spots, n_markers = X.shape
    n_types = profile.shape[0]

    # Stage 1: Detection (one-time)
    detected = detect_cell_types(X, marker_groups)

    # Stage 2: Estimation with EM
    sigma_sq = np.ones(n_markers)  # initialize uniform

    for iteration in range(max_iter):
        # E-step: Solve IQP with current weights
        weights = 1.0 / sigma_sq
        counts, alpha, beta = solve_masked_iqp(
            X, nuclei_counts, profile, detected, weights
        )

        # M-step: Update sigma_sq from residuals
        predicted = alpha + (counts @ profile) * beta
        residuals = X - predicted

        # Robust variance estimation (MAD-based)
        sigma_sq_new = (1.4826 * np.median(np.abs(residuals), axis=0)) ** 2
        sigma_sq_new = np.maximum(sigma_sq_new, 1e-6)  # floor for stability

        # Check convergence
        if np.allclose(sigma_sq, sigma_sq_new, rtol=0.01):
            break
        sigma_sq = sigma_sq_new

    return detected, counts, alpha, beta, sigma_sq
```

### Gurobi IQP Formulation

```python
def solve_masked_iqp(X, nuclei_counts, profile, detected, weights):
    """
    Solve IQP for counts, alpha, beta with detection mask.
    """
    n_spots, n_markers = X.shape
    n_types = profile.shape[0]

    model = gp.Model("detection_estimation")

    # Variables
    c = {}  # c[i,k] = cell count
    for i in range(n_spots):
        for k in range(n_types):
            if detected[i, k]:
                c[i, k] = model.addVar(
                    vtype=GRB.INTEGER,
                    lb=0,
                    ub=int(nuclei_counts[i]),
                    name=f"c_{i}_{k}"
                )
            else:
                c[i, k] = 0  # HARD ZERO - not a variable

    alpha = model.addVars(n_markers, lb=0, name="alpha")
    beta = model.addVars(n_markers, lb=1e-3, name="beta")

    # Objective: weighted sum of squared residuals
    obj = gp.QuadExpr()
    for i in range(n_spots):
        for m in range(n_markers):
            # predicted = alpha[m] + sum_k c[i,k] * profile[k,m] * beta[m]
            pred = alpha[m]
            for k in range(n_types):
                if detected[i, k] and profile[k, m] > 0:
                    pred = pred + c[i, k] * profile[k, m] * beta[m]

            residual = X[i, m] - pred
            obj += weights[m] * residual * residual

    model.setObjective(obj, GRB.MINIMIZE)

    # Constraints: cell counts sum to nuclei count
    for i in range(n_spots):
        n_i = int(nuclei_counts[i])
        if n_i > 0:
            detected_types = [k for k in range(n_types) if detected[i, k]]
            if detected_types:
                model.addConstr(
                    gp.quicksum(c[i, k] for k in detected_types) == n_i,
                    name=f"sum_{i}"
                )

    model.optimize()

    # Extract solution
    counts = np.zeros((n_spots, n_types), dtype=int)
    for i in range(n_spots):
        for k in range(n_types):
            if detected[i, k]:
                counts[i, k] = int(round(c[i, k].X))

    alpha_vals = np.array([alpha[m].X for m in range(n_markers)])
    beta_vals = np.array([beta[m].X for m in range(n_markers)])

    return counts, alpha_vals, beta_vals
```

## Key Improvements Over Current Model

| Aspect | Current Model | New Model |
|--------|---------------|-----------|
| Rare cell types | 5-6% floor (L2 regularization) | True zeros via detection mask |
| Background noise | Not modeled | α[m] baseline absorbs noise |
| Marker differences | Per-marker β only | Per-marker α, β, AND σ² |
| Sparsity | Regularization (soft) | Detection mask (hard) |
| Multi-marker types | Independent markers | Multivariate GMM captures covariance |

## Expected Benefits

1. **B cells in Region 1**: 42.9% → <5% (closer to 0.3% GT)
2. **Interpretable**: α[m] shows background level, β[m] shows signal-per-cell
3. **No hyperparameter tuning**: σ²[m] learned from data
4. **Downstream GEX**: True zeros mean no spurious gene expression allocation

## Implementation Plan

1. **New file**: `CITEgeist/model/detection_estimation.py`
2. **Functions**:
   - `detect_cell_types()` - Stage 1 GMM detection
   - `solve_masked_iqp()` - Stage 2 Gurobi optimization
   - `solve_detection_estimation()` - Combined pipeline
3. **Integration**: Add `run_detection_estimation()` to `CitegeistModel`
4. **Testing**: Benchmark on Xenium pseudo-Visium regions

## Open Questions

1. **Detection threshold**: 0.5 default, but may need tuning per dataset
2. **Spatial smoothing**: Should detection use spatial context? (neighbors)
3. **Edge case**: What if NO types detected in a spot with nuclei?

## Validation Metrics

- Proportion Pearson r (overall and per-cell-type)
- Sparsity: % spots with true zeros per cell type
- False positive rate: % spots with predicted presence but GT absence
- GEX correlation after using new proportions
