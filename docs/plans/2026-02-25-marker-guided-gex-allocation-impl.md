# Marker-Guided GEX Allocation Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Replace L2 regularization and weak proportion-based enrichment with adaptive marker-guided allocation to improve GEX deconvolution from r=0.27 to r=0.32+.

**Architecture:** Modify enrichment calculation to blend proportion-based and anchor-correlation-based signals adaptively. Remove L2 penalty and 80/20 smoothing. Keep Gurobi QP structure.

**Tech Stack:** Python 3.10, NumPy, SciPy (pearsonr), Gurobi

---

## Task 1: Extend discover_anchor_genes to Return Correlation Weights

**Files:**
- Modify: `CITEgeist/model/gex_modules.py:19-137`
- Test: `CITEgeist/tests/test_gex_modules.py` (create)

**Step 1: Write the failing test**

```python
# CITEgeist/tests/test_gex_modules.py
import numpy as np
import pytest
from CITEgeist.model.gex_modules import discover_anchor_genes


def test_discover_anchor_genes_returns_weights():
    """Test that discover_anchor_genes returns anchor weights (correlations)."""
    np.random.seed(42)
    N, M, T = 100, 50, 3

    # Create synthetic data where gene 0 correlates with type 0, etc.
    proportions = np.random.dirichlet(np.ones(T), N)
    gene_expr = np.zeros((N, M))

    # First 3 genes are markers for each cell type
    for t in range(T):
        gene_expr[:, t] = proportions[:, t] * 100 + np.random.randn(N) * 5

    # Rest are noise
    gene_expr[:, T:] = np.random.rand(N, M - T) * 10

    anchors, thresholds, weights = discover_anchor_genes(
        gene_expr, proportions, min_anchors=1, max_anchors=5
    )

    # Should return 3 values now
    assert len(weights) == T, "Should have weights for each cell type"

    # Weights should be dicts mapping gene_idx -> correlation
    for t in range(T):
        assert isinstance(weights[t], dict), f"weights[{t}] should be dict"
        # Marker gene t should be in anchors for type t with high weight
        if t in anchors[t]:
            assert weights[t][t] > 0.3, f"Marker gene {t} should have r > 0.3"
```

**Step 2: Run test to verify it fails**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest CITEgeist/tests/test_gex_modules.py::test_discover_anchor_genes_returns_weights -v`

Expected: FAIL with "cannot unpack" or "not enough values"

**Step 3: Modify discover_anchor_genes to return weights**

In `CITEgeist/model/gex_modules.py`, change the function signature and return statement:

```python
def discover_anchor_genes(
    gene_expression: np.ndarray,
    cell_proportions: np.ndarray,
    min_anchors: int = 5,
    max_anchors: int = 10,
    initial_min_correlation: float = 0.3,
    min_expressing_spots: float = 0.1,
) -> Tuple[Dict[int, List[int]], Dict[int, float], Dict[int, Dict[int, float]]]:
    """
    ...existing docstring...

    Returns:
        anchors: Dict mapping cell type index to list of gene indices
        thresholds_used: Dict mapping cell type index to correlation threshold
        anchor_weights: Dict mapping cell type index to dict of {gene_idx: correlation}
    """
    # ... existing code until line 122 ...

    anchors = {}
    thresholds_used = {}
    anchor_weights = {}  # NEW

    for t in range(T):
        # ... existing loop code ...

        # After selecting anchors, store their correlations
        anchors[t] = final_candidates[:max_anchors]
        thresholds_used[t] = selected_threshold

        # NEW: Store correlation for each anchor
        anchor_weights[t] = {}
        for g, r in all_correlations:
            if g in anchors[t]:
                anchor_weights[t][g] = r

    return anchors, thresholds_used, anchor_weights
```

**Step 4: Run test to verify it passes**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest CITEgeist/tests/test_gex_modules.py::test_discover_anchor_genes_returns_weights -v`

Expected: PASS

**Step 5: Commit**

```bash
git add CITEgeist/model/gex_modules.py CITEgeist/tests/test_gex_modules.py
git commit -m "feat(gex_modules): return anchor weights from discover_anchor_genes"
```

---

## Task 2: Add compute_proportion_enrichment Without Smoothing

**Files:**
- Modify: `CITEgeist/model/gurobi_impl.py:2410-2430`
- Test: `CITEgeist/tests/test_gex_modules.py`

**Step 1: Write the failing test**

```python
# Add to CITEgeist/tests/test_gex_modules.py
def test_compute_proportion_enrichment_no_smoothing():
    """Proportion enrichment should NOT have 80/20 smoothing."""
    from CITEgeist.model.gurobi_impl import compute_proportion_enrichment

    np.random.seed(42)
    N, T = 20, 4

    # Create data where gene is only expressed in spots with high type-0 proportion
    cell_props = np.random.dirichlet(np.ones(T), N)
    gene_expr = np.zeros(N)
    high_type0_spots = cell_props[:, 0] > 0.5
    gene_expr[high_type0_spots] = 100

    enrichment = compute_proportion_enrichment(gene_expr, cell_props)

    # Without smoothing, type 0 should dominate
    assert enrichment[0] > 0.5, "Type 0 should have majority enrichment"

    # Should sum to 1
    assert abs(enrichment.sum() - 1.0) < 1e-6, "Enrichment should sum to 1"

    # Should NOT be smoothed toward uniform (1/T = 0.25)
    # Old code would have: 0.8 * enrichment + 0.2 * 0.25 = some smoothed value
    # Check that non-enriched types can go below 0.05 (smoothing floor was ~0.05)
    assert enrichment[1:].min() < 0.1, "Non-enriched types should be low (no smoothing)"
```

**Step 2: Run test to verify it fails**

Run: `python -m pytest CITEgeist/tests/test_gex_modules.py::test_compute_proportion_enrichment_no_smoothing -v`

Expected: FAIL (function doesn't exist yet as module-level function)

**Step 3: Add compute_proportion_enrichment function**

Add to `CITEgeist/model/gurobi_impl.py` around line 2400:

```python
def compute_proportion_enrichment(
    gene_expr: np.ndarray,
    cell_type_props: np.ndarray,
    celltype_frequencies: np.ndarray = None,
) -> np.ndarray:
    """
    Compute proportion-based enrichment WITHOUT smoothing.

    Args:
        gene_expr: (N,) expression values for one gene across spots
        cell_type_props: (N, T) cell type proportions per spot
        celltype_frequencies: (T,) global cell type frequencies for normalization

    Returns:
        enrichment: (T,) enrichment scores summing to 1
    """
    T = cell_type_props.shape[1]

    # Handle zero expression
    total_expr = gene_expr.sum()
    if total_expr < 1e-10:
        return np.ones(T) / T

    # Expression-weighted cell type proportions
    weights = gene_expr / total_expr

    # Normalize by global frequency if provided
    if celltype_frequencies is not None:
        normalized_props = cell_type_props / (celltype_frequencies + 1e-10)
    else:
        normalized_props = cell_type_props

    weighted_props = np.sum(normalized_props * weights[:, np.newaxis], axis=0)
    background_props = np.mean(normalized_props, axis=0)

    epsilon = 1e-10
    enrichment = weighted_props / (background_props + epsilon)

    # Normalize to sum to 1 (NO 80/20 smoothing)
    return enrichment / (enrichment.sum() + epsilon)
```

**Step 4: Run test to verify it passes**

Run: `python -m pytest CITEgeist/tests/test_gex_modules.py::test_compute_proportion_enrichment_no_smoothing -v`

Expected: PASS

**Step 5: Commit**

```bash
git add CITEgeist/model/gurobi_impl.py CITEgeist/tests/test_gex_modules.py
git commit -m "feat(gurobi_impl): add compute_proportion_enrichment without smoothing"
```

---

## Task 3: Add compute_marker_enrichment Function

**Files:**
- Modify: `CITEgeist/model/gurobi_impl.py`
- Test: `CITEgeist/tests/test_gex_modules.py`

**Step 1: Write the failing test**

```python
def test_compute_marker_enrichment():
    """Marker enrichment uses correlation with anchor genes."""
    from CITEgeist.model.gurobi_impl import compute_marker_enrichment

    np.random.seed(42)
    N, T = 20, 3

    # Gene correlates with type-0 anchors
    anchor_expr = np.random.rand(N, T) * 10  # mean anchor expr per type
    gene_expr = anchor_expr[:, 0] + np.random.randn(N) * 0.1  # correlates with type 0

    anchor_weights = np.array([0.7, 0.5, 0.3])  # type 0 has strongest anchor

    enrichment = compute_marker_enrichment(gene_expr, anchor_expr, anchor_weights)

    # Type 0 should have highest enrichment (correlates + highest anchor weight)
    assert enrichment[0] > enrichment[1], "Type 0 should have highest enrichment"
    assert enrichment[0] > enrichment[2], "Type 0 should have highest enrichment"
    assert abs(enrichment.sum() - 1.0) < 1e-6, "Should sum to 1"
```

**Step 2: Run test to verify it fails**

Run: `python -m pytest CITEgeist/tests/test_gex_modules.py::test_compute_marker_enrichment -v`

Expected: FAIL (function doesn't exist)

**Step 3: Implement compute_marker_enrichment**

Add to `CITEgeist/model/gurobi_impl.py`:

```python
def compute_marker_enrichment(
    gene_expr: np.ndarray,
    anchor_expr: np.ndarray,
    anchor_weights: np.ndarray,
) -> np.ndarray:
    """
    Compute marker-guided enrichment via correlation with anchor genes.

    Args:
        gene_expr: (N,) expression of target gene across neighborhood
        anchor_expr: (N, T) mean anchor expression per cell type in neighborhood
        anchor_weights: (T,) weight per cell type (mean anchor correlation with proportions)

    Returns:
        enrichment: (T,) enrichment scores summing to 1
    """
    from scipy.stats import pearsonr

    T = anchor_expr.shape[1]
    enrichment = np.zeros(T)

    # Skip if gene has no variance
    if np.std(gene_expr) < 1e-10:
        return np.ones(T) / T

    for t in range(T):
        anchor_t = anchor_expr[:, t]

        # Skip if anchor has no variance
        if np.std(anchor_t) < 1e-10:
            continue

        r, _ = pearsonr(gene_expr, anchor_t)

        if not np.isnan(r):
            # Only positive correlations, weighted by anchor strength
            enrichment[t] = max(0, r) * anchor_weights[t]

    # Normalize
    if enrichment.sum() < 1e-10:
        return np.ones(T) / T

    return enrichment / enrichment.sum()
```

**Step 4: Run test to verify it passes**

Run: `python -m pytest CITEgeist/tests/test_gex_modules.py::test_compute_marker_enrichment -v`

Expected: PASS

**Step 5: Commit**

```bash
git add CITEgeist/model/gurobi_impl.py CITEgeist/tests/test_gex_modules.py
git commit -m "feat(gurobi_impl): add compute_marker_enrichment function"
```

---

## Task 4: Add compute_adaptive_enrichment Function

**Files:**
- Modify: `CITEgeist/model/gurobi_impl.py`
- Test: `CITEgeist/tests/test_gex_modules.py`

**Step 1: Write the failing test**

```python
def test_compute_adaptive_enrichment():
    """Adaptive enrichment blends proportion and marker based on variance."""
    from CITEgeist.model.gurobi_impl import compute_adaptive_enrichment

    np.random.seed(42)
    N, T = 20, 4

    # Case 1: High variance proportion enrichment -> trust proportions
    prop_enrich_peaked = np.array([0.7, 0.1, 0.1, 0.1])
    marker_enrich = np.array([0.25, 0.25, 0.25, 0.25])

    result = compute_adaptive_enrichment(prop_enrich_peaked, marker_enrich)

    # Should be close to prop_enrich (high variance = low anchor weight)
    assert result[0] > 0.5, "High-variance case should trust proportions"

    # Case 2: Low variance proportion enrichment -> use markers
    prop_enrich_flat = np.array([0.26, 0.25, 0.25, 0.24])
    marker_enrich_peaked = np.array([0.6, 0.2, 0.1, 0.1])

    result2 = compute_adaptive_enrichment(prop_enrich_flat, marker_enrich_peaked)

    # Should lean toward marker_enrich (low variance = high anchor weight)
    assert result2[0] > 0.4, "Low-variance case should use markers"
```

**Step 2: Run test to verify it fails**

Run: `python -m pytest CITEgeist/tests/test_gex_modules.py::test_compute_adaptive_enrichment -v`

Expected: FAIL (function doesn't exist)

**Step 3: Implement compute_adaptive_enrichment**

Add to `CITEgeist/model/gurobi_impl.py`:

```python
def compute_adaptive_enrichment(
    prop_enrichment: np.ndarray,
    marker_enrichment: np.ndarray,
    max_variance: float = 0.25,
) -> np.ndarray:
    """
    Adaptively blend proportion and marker enrichment based on proportion variance.

    High proportion variance = trust proportions (peaked distribution)
    Low proportion variance = use marker guidance (flat distribution needs help)

    Args:
        prop_enrichment: (T,) proportion-based enrichment (normalized)
        marker_enrichment: (T,) marker-guided enrichment (normalized)
        max_variance: Maximum expected variance for normalized distribution

    Returns:
        blended: (T,) adaptive blend summing to 1
    """
    # Variance of proportion enrichment across cell types
    variance = np.var(prop_enrichment)

    # anchor_weight: 0 when variance is high (peaked), 1 when variance is low (flat)
    anchor_weight = 1 - min(1.0, variance / max_variance)

    # Blend
    blended = (1 - anchor_weight) * prop_enrichment + anchor_weight * marker_enrichment

    # Normalize (should already sum to ~1, but ensure)
    return blended / (blended.sum() + 1e-10)
```

**Step 4: Run test to verify it passes**

Run: `python -m pytest CITEgeist/tests/test_gex_modules.py::test_compute_adaptive_enrichment -v`

Expected: PASS

**Step 5: Commit**

```bash
git add CITEgeist/model/gurobi_impl.py CITEgeist/tests/test_gex_modules.py
git commit -m "feat(gurobi_impl): add compute_adaptive_enrichment with variance-based blending"
```

---

## Task 5: Precompute Anchor Expression Matrix

**Files:**
- Modify: `CITEgeist/model/gurobi_impl.py` (inside `deconvolute_spot_with_neighbors_with_prior`)

**Step 1: Write the failing test**

```python
def test_precompute_anchor_expression():
    """Test that anchor expression matrix is computed correctly."""
    from CITEgeist.model.gurobi_impl import precompute_anchor_expression

    np.random.seed(42)
    N, M, T = 50, 100, 3

    gene_expr = np.random.rand(N, M) * 100

    # Anchors: type 0 -> genes [0,1], type 1 -> genes [2,3], type 2 -> genes [4,5]
    anchor_genes = {0: [0, 1], 1: [2, 3], 2: [4, 5]}
    anchor_weights = {
        0: {0: 0.8, 1: 0.6},  # gene 0 has r=0.8, gene 1 has r=0.6
        1: {2: 0.7, 3: 0.5},
        2: {4: 0.9, 5: 0.4},
    }

    anchor_expr, type_weights = precompute_anchor_expression(
        gene_expr, anchor_genes, anchor_weights
    )

    # Shape: (N, T) - mean anchor expression per type
    assert anchor_expr.shape == (N, T)

    # Type weights: mean correlation per type
    assert len(type_weights) == T
    assert type_weights[0] == pytest.approx(0.7, abs=0.01)  # mean(0.8, 0.6)
    assert type_weights[2] == pytest.approx(0.65, abs=0.01)  # mean(0.9, 0.4)
```

**Step 2: Run test to verify it fails**

Run: `python -m pytest CITEgeist/tests/test_gex_modules.py::test_precompute_anchor_expression -v`

Expected: FAIL (function doesn't exist)

**Step 3: Implement precompute_anchor_expression**

Add to `CITEgeist/model/gurobi_impl.py`:

```python
def precompute_anchor_expression(
    gene_expression: np.ndarray,
    anchor_genes: Dict[int, List[int]],
    anchor_weights: Dict[int, Dict[int, float]],
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Precompute weighted mean anchor expression per cell type.

    Args:
        gene_expression: (N, M) expression matrix
        anchor_genes: {cell_type: [gene_indices]}
        anchor_weights: {cell_type: {gene_idx: correlation}}

    Returns:
        anchor_expr: (N, T) weighted mean anchor expression per cell type
        type_weights: (T,) mean anchor weight per cell type
    """
    N = gene_expression.shape[0]
    T = len(anchor_genes)

    anchor_expr = np.zeros((N, T))
    type_weights = np.zeros(T)

    for t in range(T):
        if t not in anchor_genes or not anchor_genes[t]:
            # No anchors for this type - use uniform
            anchor_expr[:, t] = 1.0
            type_weights[t] = 0.1  # low weight
            continue

        genes = anchor_genes[t]
        weights_t = anchor_weights.get(t, {})

        # Weighted mean of anchor expressions
        total_weight = 0.0
        for g in genes:
            w = weights_t.get(g, 0.3)  # default weight if missing
            anchor_expr[:, t] += w * gene_expression[:, g]
            total_weight += w

        if total_weight > 0:
            anchor_expr[:, t] /= total_weight
            type_weights[t] = total_weight / len(genes)  # mean weight
        else:
            anchor_expr[:, t] = 1.0
            type_weights[t] = 0.1

    return anchor_expr, type_weights
```

**Step 4: Run test to verify it passes**

Run: `python -m pytest CITEgeist/tests/test_gex_modules.py::test_precompute_anchor_expression -v`

Expected: PASS

**Step 5: Commit**

```bash
git add CITEgeist/model/gurobi_impl.py CITEgeist/tests/test_gex_modules.py
git commit -m "feat(gurobi_impl): add precompute_anchor_expression for marker guidance"
```

---

## Task 6: Integrate Adaptive Enrichment into Deconvolution

**Files:**
- Modify: `CITEgeist/model/gurobi_impl.py:2432-2455` (enrichment loop)

**Step 1: Understand current code**

Current enrichment loop (lines 2432-2454):
```python
# Compute expression-aware enrichment for each gene
gene_specific_enrichment = np.zeros((M, T))

for k in range(M):
    local_enrich = compute_expression_aware_enrichment(...)
    global_enrich = compute_expression_aware_enrichment(...)
    gene_specific_enrichment[k] = local_enrichment_weight * local_enrich + global_enrichment_weight * global_enrich

# Apply module-aware enrichment if anchors provided
if anchor_genes is not None and module_weight > 0:
    gene_specific_enrichment = compute_module_aware_enrichment(...)
```

**Step 2: Replace with adaptive enrichment**

Replace the enrichment calculation block with:

```python
# Compute expression-aware enrichment for each gene using adaptive marker guidance
gene_specific_enrichment = np.zeros((M, T))

# Precompute anchor expression for neighborhood if anchors provided
if anchor_genes is not None and anchor_weights is not None:
    neighborhood_anchor_expr, type_weights = precompute_anchor_expression(
        neighborhood_expression_data, anchor_genes, anchor_weights
    )
else:
    neighborhood_anchor_expr = None
    type_weights = None

for k in range(M):
    # Proportion-based enrichment (no smoothing)
    prop_enrich = compute_proportion_enrichment(
        gene_expr=neighborhood_expression_data[:, k],
        cell_type_props=neighborhood_cell_type_numbers,
        celltype_frequencies=celltype_frequencies,
    )

    # If anchors available, compute adaptive blend
    if neighborhood_anchor_expr is not None:
        marker_enrich = compute_marker_enrichment(
            gene_expr=neighborhood_expression_data[:, k],
            anchor_expr=neighborhood_anchor_expr,
            anchor_weights=type_weights,
        )
        gene_specific_enrichment[k] = compute_adaptive_enrichment(
            prop_enrich, marker_enrich
        )
    else:
        gene_specific_enrichment[k] = prop_enrich
```

**Step 3: Update function signature**

Add `anchor_weights` parameter to `deconvolute_spot_with_neighbors_with_prior()`:

```python
def deconvolute_spot_with_neighbors_with_prior(
    ...,
    anchor_genes: Dict[int, List[int]] = None,
    anchor_weights: Dict[int, Dict[int, float]] = None,  # NEW
    ...
)
```

**Step 4: Manual verification**

Run existing test to ensure no regression:
```bash
python -m pytest CITEgeist/tests/test_citegeist_simulated.py -v -k "gex" --timeout=300
```

**Step 5: Commit**

```bash
git add CITEgeist/model/gurobi_impl.py
git commit -m "feat(gurobi_impl): integrate adaptive marker-guided enrichment"
```

---

## Task 7: Remove L2 Penalty

**Files:**
- Modify: `CITEgeist/model/gurobi_impl.py:2549-2550`

**Step 1: Locate L2 penalty**

Current code (lines 2549-2550):
```python
if lambda_gex_reg > 0:
    obj_terms.append(-lambda_gex_reg * X[j, k] * X[j, k])
```

**Step 2: Remove L2 penalty block**

Delete or comment out lines 2549-2550:
```python
# L2 regularization removed - adaptive enrichment provides sufficient guidance
# if lambda_gex_reg > 0:
#     obj_terms.append(-lambda_gex_reg * X[j, k] * X[j, k])
```

**Step 3: Update docstring**

Update the function docstring to note that L2 is deprecated:
```python
"""
...
lambda_gex_reg (float): DEPRECATED - L2 regularization weight. No longer used
    when use_marker_guidance=True. Kept for backward compatibility.
...
"""
```

**Step 4: Commit**

```bash
git add CITEgeist/model/gurobi_impl.py
git commit -m "feat(gurobi_impl): remove L2 penalty when using marker guidance"
```

---

## Task 8: Update citegeist_model.py to Pass Anchor Weights

**Files:**
- Modify: `CITEgeist/model/citegeist_model.py:1477-1488`

**Step 1: Update anchor discovery call**

Change from:
```python
anchor_genes, anchor_thresholds = discover_anchor_genes(...)
```

To:
```python
anchor_genes, anchor_thresholds, anchor_weights = discover_anchor_genes(...)
```

**Step 2: Pass anchor_weights to optimize_gene_expression**

Add to the `optimize_gene_expression()` call around line 1564:
```python
spotwise_profiles = optimize_gene_expression(
    ...
    anchor_genes=anchor_genes,
    anchor_weights=anchor_weights,  # NEW
    ...
)
```

**Step 3: Update optimize_gene_expression signature**

In `gurobi_impl.py`, add `anchor_weights` parameter to `optimize_gene_expression()` and pass it through to `deconvolute_spot_with_neighbors_with_prior()`.

**Step 4: Commit**

```bash
git add CITEgeist/model/citegeist_model.py CITEgeist/model/gurobi_impl.py
git commit -m "feat: wire anchor_weights through to deconvolution"
```

---

## Task 9: Add use_marker_guidance Parameter

**Files:**
- Modify: `CITEgeist/model/citegeist_model.py:1353-1374`

**Step 1: Add parameter**

Add to `run_cell_expression_pass1()` signature:
```python
def run_cell_expression_pass1(
    ...,
    use_marker_guidance: bool = True,  # NEW - enable adaptive marker enrichment
    ...
)
```

**Step 2: Use parameter to control behavior**

```python
# Discover anchor genes if using marker guidance
anchor_genes = None
anchor_weights = None
if use_marker_guidance and not (self.resolution == "cell"):
    from .gex_modules import discover_anchor_genes
    ...
    anchor_genes, anchor_thresholds, anchor_weights = discover_anchor_genes(...)
```

**Step 3: Update docstring**

```python
"""
...
use_marker_guidance (bool): If True, use adaptive marker-guided enrichment
    instead of proportion-only enrichment. Recommended for improved weak gene
    allocation. Default True.
...
"""
```

**Step 4: Commit**

```bash
git add CITEgeist/model/citegeist_model.py
git commit -m "feat: add use_marker_guidance parameter for A/B testing"
```

---

## Task 10: Integration Test on Xenium Region

**Files:**
- Test script: `Benchmarking/xenium_benchmarking/CITEgeist/slurm/sbatch_marker_guidance_test.sh`

**Step 1: Create test SLURM script**

```bash
#!/bin/bash
#SBATCH --job-name=marker_guidance_test
#SBATCH --output=logs/marker_guidance_%j.out
#SBATCH --error=logs/marker_guidance_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

module load gurobi/12.0.3

eval "$(conda shell.bash hook)"
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist

# Run on region 2 (not the outlier region 1)
python Benchmarking/xenium_benchmarking/CITEgeist/src/run_multimodal_xenium.py \
    --region 2 \
    --output_dir Benchmarking/xenium_benchmarking/CITEgeist/output/marker_guidance_test \
    --use_marker_guidance
```

**Step 2: Submit and monitor**

```bash
sbatch Benchmarking/xenium_benchmarking/CITEgeist/slurm/sbatch_marker_guidance_test.sh
```

**Step 3: Evaluate results**

```bash
python Benchmarking/xenium_benchmarking/evaluation/src/evaluate_gex_spatial.py \
    --method CITEgeist \
    --experiment marker_guidance_test \
    --region 2
```

**Step 4: Compare metrics**

Expected improvements:
- per_gene_pearson: 0.27 → 0.32+ (+0.05)
- variance_ratio: 3.0 → <2.5 (-0.5)
- marker_pearson: ≥0.60 (stable)

**Step 5: Commit test script**

```bash
git add Benchmarking/xenium_benchmarking/CITEgeist/slurm/sbatch_marker_guidance_test.sh
git commit -m "test: add integration test for marker-guided allocation"
```

---

## Summary

| Task | Description | Files |
|------|-------------|-------|
| 1 | Extend discover_anchor_genes | gex_modules.py |
| 2 | Add compute_proportion_enrichment | gurobi_impl.py |
| 3 | Add compute_marker_enrichment | gurobi_impl.py |
| 4 | Add compute_adaptive_enrichment | gurobi_impl.py |
| 5 | Add precompute_anchor_expression | gurobi_impl.py |
| 6 | Integrate into deconvolution loop | gurobi_impl.py |
| 7 | Remove L2 penalty | gurobi_impl.py |
| 8 | Wire anchor_weights through | citegeist_model.py, gurobi_impl.py |
| 9 | Add use_marker_guidance param | citegeist_model.py |
| 10 | Integration test | SLURM script |

**Total estimated time**: 2-3 hours

**Success criteria**:
- All unit tests pass
- per_gene_pearson improves ≥0.05
- variance_ratio decreases
- marker_pearson stable
