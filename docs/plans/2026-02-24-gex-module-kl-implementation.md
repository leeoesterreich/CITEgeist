# GEX Module-Aware Enrichment + KL Regularization Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Improve GEX deconvolution by replacing L2 regularization with softmax KL-divergence and adding module-aware enrichment based on anchor gene spatial correlation.

**Architecture:** Three-stage pipeline: (1) discover anchor genes per cell type via correlation with proportions, (2) compute per-spot module-aware enrichment by correlating each gene with anchors in neighborhood, (3) allocate gene counts using KL-divergence from softmax target instead of L2.

**Tech Stack:** Python 3.10, NumPy, SciPy (pearsonr), Gurobi (QP), pytest

**Design Doc:** `docs/plans/2026-02-24-gex-module-kl-design.md`

---

## Task 1: Implement Anchor Gene Discovery

**Files:**
- Create: `CITEgeist/model/gex_modules.py`
- Test: `tests/test_gex_modules.py`

### Step 1.1: Write failing test for anchor discovery

```python
# tests/test_gex_modules.py
"""Tests for GEX module-aware enrichment and KL regularization."""

import numpy as np
import pytest
from scipy.stats import pearsonr


class TestAnchorGeneDiscovery:
    """Tests for discover_anchor_genes function."""

    def test_discovers_correlated_genes(self):
        """Anchor discovery finds genes correlated with cell type proportions."""
        np.random.seed(42)
        n_spots, n_genes, n_types = 100, 50, 3

        # Create cell type proportions
        proportions = np.random.dirichlet([1, 1, 1], size=n_spots)

        # Create gene expression where first 5 genes correlate with type 0
        gene_expression = np.random.rand(n_spots, n_genes) * 10
        for g in range(5):
            gene_expression[:, g] = proportions[:, 0] * 100 + np.random.randn(n_spots) * 5

        from CITEgeist.model.gex_modules import discover_anchor_genes

        anchors, thresholds = discover_anchor_genes(
            gene_expression=gene_expression,
            cell_proportions=proportions,
            min_anchors=3,
            max_anchors=10,
        )

        # Type 0 should have anchors including genes 0-4
        assert 0 in anchors
        assert len(anchors[0]) >= 3
        # At least some of genes 0-4 should be anchors for type 0
        anchor_set = set(anchors[0])
        expected_anchors = set(range(5))
        assert len(anchor_set & expected_anchors) >= 2

    def test_threshold_adapts_for_weak_celltypes(self):
        """Threshold drops until min anchors reached for weak cell types."""
        np.random.seed(42)
        n_spots, n_genes, n_types = 100, 50, 3

        proportions = np.random.dirichlet([1, 1, 1], size=n_spots)

        # Type 0: strong signal (5 genes with r>0.3)
        # Type 1: weak signal (only 3 genes with r>0.1)
        # Type 2: medium signal
        gene_expression = np.random.rand(n_spots, n_genes) * 10

        # Type 0: 5 strongly correlated genes
        for g in range(5):
            gene_expression[:, g] = proportions[:, 0] * 100 + np.random.randn(n_spots) * 5

        # Type 1: 5 weakly correlated genes (r ~ 0.15-0.25)
        for g in range(5, 10):
            noise = np.random.randn(n_spots) * 30
            gene_expression[:, g] = proportions[:, 1] * 50 + noise

        from CITEgeist.model.gex_modules import discover_anchor_genes

        anchors, thresholds = discover_anchor_genes(
            gene_expression=gene_expression,
            cell_proportions=proportions,
            min_anchors=5,
            max_anchors=10,
        )

        # Type 0 should have higher threshold than type 1
        assert thresholds[0] >= thresholds[1]
        # Both should have at least min_anchors
        assert len(anchors[0]) >= 5
        assert len(anchors[1]) >= 5

    def test_respects_max_anchors(self):
        """Never returns more than max_anchors per cell type."""
        np.random.seed(42)
        n_spots, n_genes, n_types = 100, 50, 2

        proportions = np.random.dirichlet([1, 1], size=n_spots)
        gene_expression = np.random.rand(n_spots, n_genes) * 10

        # Make all genes correlate with type 0
        for g in range(n_genes):
            gene_expression[:, g] = proportions[:, 0] * 100 + np.random.randn(n_spots) * 10

        from CITEgeist.model.gex_modules import discover_anchor_genes

        anchors, _ = discover_anchor_genes(
            gene_expression=gene_expression,
            cell_proportions=proportions,
            min_anchors=5,
            max_anchors=10,
        )

        assert len(anchors[0]) <= 10
```

### Step 1.2: Run test to verify it fails

```bash
cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist
module load gurobi/12.0.3
eval "$(conda shell.bash hook)" && conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env
pytest tests/test_gex_modules.py::TestAnchorGeneDiscovery -v
```

Expected: FAIL with `ModuleNotFoundError: No module named 'CITEgeist.model.gex_modules'`

### Step 1.3: Implement discover_anchor_genes

```python
# CITEgeist/model/gex_modules.py
"""
GEX module-aware enrichment and KL regularization for gene expression deconvolution.

This module implements:
1. Anchor gene discovery - find genes correlated with cell type proportions
2. Module-aware enrichment - boost gene enrichment based on anchor correlation
3. Softmax KL allocation - replace L2 regularization with KL-divergence
"""

import logging
from typing import Dict, List, Tuple

import numpy as np
from scipy.stats import pearsonr

logger = logging.getLogger(__name__)


def discover_anchor_genes(
    gene_expression: np.ndarray,
    cell_proportions: np.ndarray,
    min_anchors: int = 5,
    max_anchors: int = 10,
    initial_min_correlation: float = 0.3,
    min_expressing_spots: float = 0.1,
) -> Tuple[Dict[int, List[int]], Dict[int, float]]:
    """
    Discover anchor genes for each cell type based on correlation with proportions.

    Args:
        gene_expression: (N_spots x M_genes) expression matrix
        cell_proportions: (N_spots x T_celltypes) proportion matrix
        min_anchors: Minimum anchors per cell type (threshold adapts to reach this)
        max_anchors: Maximum anchors per cell type
        initial_min_correlation: Starting correlation threshold
        min_expressing_spots: Minimum fraction of spots where gene must be expressed

    Returns:
        anchors: {cell_type_idx: [gene_indices]}
        thresholds_used: {cell_type_idx: correlation_threshold}
    """
    N, M = gene_expression.shape
    T = cell_proportions.shape[1]

    # Threshold sequence: 0.30 -> 0.25 -> 0.20 -> 0.15 -> 0.10
    threshold_sequence = [0.30, 0.25, 0.20, 0.15, 0.10]

    anchors = {}
    thresholds_used = {}

    for t in range(T):
        prop_vector = cell_proportions[:, t]

        # Compute correlations for all genes
        all_correlations = []
        for g in range(M):
            gene_vector = gene_expression[:, g]

            # Filter: gene must be expressed in enough spots
            expressing_fraction = (gene_vector > 0).mean()
            if expressing_fraction < min_expressing_spots:
                continue

            # Filter: need variance for correlation
            if np.std(gene_vector) < 1e-10 or np.std(prop_vector) < 1e-10:
                continue

            r, p = pearsonr(gene_vector, prop_vector)
            if p < 0.05 and not np.isnan(r):
                all_correlations.append((g, r))

        # Sort by correlation descending
        all_correlations.sort(key=lambda x: -x[1])

        # Step through thresholds until we have min_anchors
        selected_threshold = threshold_sequence[-1]  # default to floor
        for threshold in threshold_sequence:
            candidates = [g for g, r in all_correlations if r >= threshold]
            if len(candidates) >= min_anchors:
                selected_threshold = threshold
                break

        # Select anchors at chosen threshold (capped at max)
        final_candidates = [g for g, r in all_correlations if r >= selected_threshold]
        anchors[t] = final_candidates[:max_anchors]
        thresholds_used[t] = selected_threshold

        if len(anchors[t]) < min_anchors:
            logger.warning(
                f"Cell type {t}: only {len(anchors[t])} anchors found "
                f"(min={min_anchors}, threshold={selected_threshold:.2f})"
            )

    return anchors, thresholds_used
```

### Step 1.4: Run test to verify it passes

```bash
pytest tests/test_gex_modules.py::TestAnchorGeneDiscovery -v
```

Expected: All 3 tests PASS

### Step 1.5: Commit

```bash
git add CITEgeist/model/gex_modules.py tests/test_gex_modules.py
git commit -m "feat: add anchor gene discovery for GEX module enrichment"
```

---

## Task 2: Implement Module-Aware Enrichment

**Files:**
- Modify: `CITEgeist/model/gex_modules.py`
- Test: `tests/test_gex_modules.py`

### Step 2.1: Write failing test for module-aware enrichment

```python
# Add to tests/test_gex_modules.py

class TestModuleAwareEnrichment:
    """Tests for compute_module_aware_enrichment function."""

    def test_boosts_enrichment_for_correlated_genes(self):
        """Genes correlated with anchors get boosted enrichment."""
        np.random.seed(42)
        n_neighbors, n_genes, n_types = 20, 10, 3

        # Neighborhood expression
        neighborhood_expr = np.random.rand(n_neighbors, n_genes) * 10

        # Make gene 5 correlate with genes 0,1 (anchors for type 0)
        neighborhood_expr[:, 0] = np.random.rand(n_neighbors) * 10
        neighborhood_expr[:, 1] = neighborhood_expr[:, 0] + np.random.randn(n_neighbors) * 0.5
        neighborhood_expr[:, 5] = neighborhood_expr[:, 0] * 0.8 + np.random.randn(n_neighbors) * 1

        # Base enrichment: gene 5 has uniform enrichment
        base_enrichment = np.ones((n_genes, n_types)) / n_types

        # Anchors: type 0 has genes 0,1
        anchor_genes = {0: [0, 1], 1: [2, 3], 2: [4]}

        from CITEgeist.model.gex_modules import compute_module_aware_enrichment

        adjusted = compute_module_aware_enrichment(
            spot_idx=0,
            neighborhood_expression=neighborhood_expr,
            base_enrichment=base_enrichment,
            anchor_genes=anchor_genes,
            module_weight=0.5,
        )

        # Gene 5 should have higher enrichment for type 0 now
        assert adjusted[5, 0] > base_enrichment[5, 0]
        # Should still sum to 1
        assert np.allclose(adjusted[5, :].sum(), 1.0)

    def test_skips_small_neighborhoods(self):
        """Returns base enrichment if neighborhood too small."""
        n_neighbors, n_genes, n_types = 5, 10, 3  # < min_neighbors_for_corr

        neighborhood_expr = np.random.rand(n_neighbors, n_genes) * 10
        base_enrichment = np.ones((n_genes, n_types)) / n_types
        anchor_genes = {0: [0, 1], 1: [2, 3], 2: [4]}

        from CITEgeist.model.gex_modules import compute_module_aware_enrichment

        adjusted = compute_module_aware_enrichment(
            spot_idx=0,
            neighborhood_expression=neighborhood_expr,
            base_enrichment=base_enrichment,
            anchor_genes=anchor_genes,
            module_weight=0.5,
            min_neighbors_for_corr=10,
        )

        # Should return base enrichment unchanged
        assert np.allclose(adjusted, base_enrichment)

    def test_only_positive_correlations_boost(self):
        """Negative correlations don't boost enrichment."""
        np.random.seed(42)
        n_neighbors, n_genes, n_types = 20, 10, 3

        neighborhood_expr = np.random.rand(n_neighbors, n_genes) * 10

        # Make gene 5 NEGATIVELY correlate with type 0 anchors
        neighborhood_expr[:, 0] = np.random.rand(n_neighbors) * 10
        neighborhood_expr[:, 5] = 10 - neighborhood_expr[:, 0] + np.random.randn(n_neighbors) * 0.5

        base_enrichment = np.ones((n_genes, n_types)) / n_types
        anchor_genes = {0: [0], 1: [2], 2: [4]}

        from CITEgeist.model.gex_modules import compute_module_aware_enrichment

        adjusted = compute_module_aware_enrichment(
            spot_idx=0,
            neighborhood_expression=neighborhood_expr,
            base_enrichment=base_enrichment,
            anchor_genes=anchor_genes,
            module_weight=0.5,
        )

        # Gene 5 should NOT have boosted enrichment for type 0
        # (negative correlation gets clipped to 0)
        assert adjusted[5, 0] <= base_enrichment[5, 0] + 0.01
```

### Step 2.2: Run test to verify it fails

```bash
pytest tests/test_gex_modules.py::TestModuleAwareEnrichment -v
```

Expected: FAIL with `ImportError: cannot import name 'compute_module_aware_enrichment'`

### Step 2.3: Implement compute_module_aware_enrichment

```python
# Add to CITEgeist/model/gex_modules.py

def compute_module_aware_enrichment(
    spot_idx: int,
    neighborhood_expression: np.ndarray,
    base_enrichment: np.ndarray,
    anchor_genes: Dict[int, List[int]],
    module_weight: float = 0.5,
    min_neighbors_for_corr: int = 10,
) -> np.ndarray:
    """
    Compute module-aware enrichment by correlating genes with anchors in neighborhood.

    Args:
        spot_idx: Index of center spot (for logging)
        neighborhood_expression: (N_neighbors x M_genes) expression in neighborhood
        base_enrichment: (M_genes x T_celltypes) base enrichment scores
        anchor_genes: {cell_type_idx: [gene_indices]} anchor genes per type
        module_weight: Weight for module signal vs base enrichment (0-1)
        min_neighbors_for_corr: Minimum neighbors for reliable correlation

    Returns:
        adjusted_enrichment: (M_genes x T_celltypes) module-adjusted enrichment
    """
    M, T = base_enrichment.shape
    adjusted = base_enrichment.copy()

    # Skip if neighborhood too small for reliable correlation
    if neighborhood_expression.shape[0] < min_neighbors_for_corr:
        return adjusted

    for g in range(M):
        gene_expr = neighborhood_expression[:, g]

        # Skip if gene has no variance in neighborhood
        if np.std(gene_expr) < 1e-6:
            continue

        # Compute correlation with each cell type's anchor genes
        module_scores = np.zeros(T)
        for t in range(T):
            if t not in anchor_genes or not anchor_genes[t]:
                continue

            # Average correlation with this cell type's anchors
            anchor_corrs = []
            for anchor_idx in anchor_genes[t]:
                if anchor_idx >= neighborhood_expression.shape[1]:
                    continue
                anchor_expr = neighborhood_expression[:, anchor_idx]
                if np.std(anchor_expr) > 1e-6:
                    r, _ = pearsonr(gene_expr, anchor_expr)
                    if not np.isnan(r):
                        anchor_corrs.append(max(0, r))  # only positive correlations

            if anchor_corrs:
                module_scores[t] = np.mean(anchor_corrs)

        # Normalize module scores to sum to 1 (if any signal)
        if module_scores.sum() > 0:
            module_scores = module_scores / module_scores.sum()
        else:
            module_scores = np.ones(T) / T  # fallback to uniform

        # Combine base enrichment with module signal
        adjusted[g, :] = (
            (1 - module_weight) * base_enrichment[g, :] +
            module_weight * module_scores
        )

        # Re-normalize to sum to 1
        row_sum = adjusted[g, :].sum()
        if row_sum > 0:
            adjusted[g, :] = adjusted[g, :] / row_sum

    return adjusted
```

### Step 2.4: Run test to verify it passes

```bash
pytest tests/test_gex_modules.py::TestModuleAwareEnrichment -v
```

Expected: All 3 tests PASS

### Step 2.5: Commit

```bash
git add CITEgeist/model/gex_modules.py tests/test_gex_modules.py
git commit -m "feat: add module-aware enrichment via anchor correlation"
```

---

## Task 3: Implement Softmax KL Objective Builder

**Files:**
- Modify: `CITEgeist/model/gex_modules.py`
- Test: `tests/test_gex_modules.py`

### Step 3.1: Write failing test for softmax target computation

```python
# Add to tests/test_gex_modules.py

class TestSoftmaxKLAllocation:
    """Tests for softmax KL-divergence allocation."""

    def test_compute_softmax_target(self):
        """Softmax target concentrates on high-enrichment cell types."""
        enrichment = np.array([0.4, 0.3, 0.1, 0.1, 0.1])

        from CITEgeist.model.gex_modules import compute_softmax_target

        # Sharp temperature
        target_sharp = compute_softmax_target(enrichment, temperature=0.3)

        # Soft temperature
        target_soft = compute_softmax_target(enrichment, temperature=1.0)

        # Both should sum to 1
        assert np.allclose(target_sharp.sum(), 1.0)
        assert np.allclose(target_soft.sum(), 1.0)

        # Sharp should concentrate more on top cell type
        assert target_sharp[0] > target_soft[0]

        # Top cell type should have highest probability
        assert target_sharp[0] > target_sharp[1] > target_sharp[2]

    def test_compute_kl_penalty_terms(self):
        """KL penalty penalizes deviation from target."""
        from CITEgeist.model.gex_modules import compute_kl_penalty_coefficients

        target = np.array([0.5, 0.3, 0.2])
        total_counts = 100
        lambda_kl = 0.1

        coeffs = compute_kl_penalty_coefficients(
            target_distribution=target,
            total_counts=total_counts,
            lambda_kl=lambda_kl,
        )

        # Should return coefficients for quadratic penalty
        assert 'target_counts' in coeffs
        assert 'penalty_weight' in coeffs

        # Target counts should match target * total
        assert np.allclose(coeffs['target_counts'], target * total_counts)

    def test_kl_replaces_l2_in_objective(self):
        """Verify KL penalty structure differs from L2."""
        from CITEgeist.model.gex_modules import compute_kl_penalty_coefficients

        # Uniform target (like L2 would produce)
        uniform_target = np.array([0.33, 0.33, 0.34])

        # Concentrated target (KL should produce)
        concentrated_target = np.array([0.6, 0.3, 0.1])

        coeffs_uniform = compute_kl_penalty_coefficients(uniform_target, 100, 0.1)
        coeffs_concentrated = compute_kl_penalty_coefficients(concentrated_target, 100, 0.1)

        # Target counts should differ
        assert not np.allclose(
            coeffs_uniform['target_counts'],
            coeffs_concentrated['target_counts']
        )
```

### Step 3.2: Run test to verify it fails

```bash
pytest tests/test_gex_modules.py::TestSoftmaxKLAllocation -v
```

Expected: FAIL with `ImportError: cannot import name 'compute_softmax_target'`

### Step 3.3: Implement softmax target and KL penalty functions

```python
# Add to CITEgeist/model/gex_modules.py

def compute_softmax_target(
    enrichment: np.ndarray,
    temperature: float = 0.3,
) -> np.ndarray:
    """
    Compute softmax target distribution from enrichment scores.

    Args:
        enrichment: (T,) enrichment scores per cell type
        temperature: Softmax temperature (lower = sharper)

    Returns:
        target: (T,) probability distribution summing to 1
    """
    logits = enrichment / temperature
    logits = logits - logits.max()  # numerical stability
    exp_logits = np.exp(logits)
    target = exp_logits / exp_logits.sum()

    # Clip to avoid zeros (for KL stability)
    target = np.clip(target, 1e-6, 1.0)
    target = target / target.sum()

    return target


def compute_kl_penalty_coefficients(
    target_distribution: np.ndarray,
    total_counts: int,
    lambda_kl: float = 0.1,
) -> Dict[str, np.ndarray]:
    """
    Compute coefficients for KL-divergence penalty in Gurobi objective.

    The KL penalty is approximated as quadratic: λ × (X[j] - target[j])²
    This pulls allocations toward the target distribution.

    Args:
        target_distribution: (T,) target probabilities
        total_counts: Total counts to allocate for this gene
        lambda_kl: Penalty weight

    Returns:
        Dict with 'target_counts' and 'penalty_weight' for Gurobi objective
    """
    target_counts = target_distribution * total_counts

    # Normalize penalty by total_counts to make lambda_kl scale-invariant
    penalty_weight = lambda_kl / (total_counts + 1)

    return {
        'target_counts': target_counts,
        'penalty_weight': penalty_weight,
    }
```

### Step 3.4: Run test to verify it passes

```bash
pytest tests/test_gex_modules.py::TestSoftmaxKLAllocation -v
```

Expected: All 3 tests PASS

### Step 3.5: Commit

```bash
git add CITEgeist/model/gex_modules.py tests/test_gex_modules.py
git commit -m "feat: add softmax target and KL penalty computation"
```

---

## Task 4: Integrate into gurobi_impl.py

**Files:**
- Modify: `CITEgeist/model/gurobi_impl.py:2252-2451` (deconvolute_spot_with_neighbors_with_prior)
- Test: `tests/test_gex_modules.py`

### Step 4.1: Write integration test

```python
# Add to tests/test_gex_modules.py

class TestGurobiIntegration:
    """Tests for integration with Gurobi optimization."""

    @pytest.mark.requires_gurobi
    def test_kl_objective_produces_valid_allocation(self):
        """KL objective produces valid count allocations."""
        import gurobipy as gp
        from gurobipy import GRB
        from CITEgeist.model.gex_modules import (
            compute_softmax_target,
            compute_kl_penalty_coefficients,
        )

        n_types = 3
        total_counts = 100
        enrichment = np.array([0.5, 0.3, 0.2])
        props = np.array([0.4, 0.4, 0.2])

        target = compute_softmax_target(enrichment, temperature=0.3)
        kl_coeffs = compute_kl_penalty_coefficients(target, total_counts, lambda_kl=0.1)

        # Build simple Gurobi model
        model = gp.Model("test_kl")
        model.setParam("OutputFlag", 0)

        X = {}
        for j in range(n_types):
            X[j] = model.addVar(vtype=GRB.CONTINUOUS, lb=0, ub=total_counts, name=f"X_{j}")

        # Conservation constraint
        model.addConstr(gp.quicksum(X[j] for j in range(n_types)) == total_counts)

        # KL-style objective
        obj_terms = []
        for j in range(n_types):
            # Enrichment term
            obj_terms.append(enrichment[j] * props[j] * X[j])
            # KL penalty (quadratic deviation from target)
            target_j = kl_coeffs['target_counts'][j]
            penalty = kl_coeffs['penalty_weight'] * (X[j] - target_j) * (X[j] - target_j)
            obj_terms.append(-penalty)

        model.setObjective(gp.quicksum(obj_terms), GRB.MAXIMIZE)
        model.optimize()

        assert model.status == GRB.OPTIMAL

        # Extract solution
        allocation = np.array([X[j].X for j in range(n_types)])

        # Should sum to total_counts
        assert np.isclose(allocation.sum(), total_counts)

        # Should be closer to target than uniform
        uniform = np.ones(n_types) * total_counts / n_types
        dist_to_target = np.sum((allocation - kl_coeffs['target_counts'])**2)
        dist_to_uniform = np.sum((allocation - uniform)**2)
        assert dist_to_target < dist_to_uniform
```

### Step 4.2: Run test to verify it passes (this one should pass with existing code)

```bash
pytest tests/test_gex_modules.py::TestGurobiIntegration -v
```

Expected: PASS (uses existing Gurobi, just tests our helper functions)

### Step 4.3: Commit

```bash
git add tests/test_gex_modules.py
git commit -m "test: add Gurobi integration test for KL objective"
```

---

## Task 5: Modify deconvolute_spot_with_neighbors_with_prior

**Files:**
- Modify: `CITEgeist/model/gurobi_impl.py:2252-2451`

### Step 5.1: Add imports and new parameters

Add to top of `gurobi_impl.py` imports section:

```python
from .gex_modules import (
    discover_anchor_genes,
    compute_module_aware_enrichment,
    compute_softmax_target,
    compute_kl_penalty_coefficients,
)
```

### Step 5.2: Modify function signature

Change `deconvolute_spot_with_neighbors_with_prior` signature at line ~2252:

```python
def deconvolute_spot_with_neighbors_with_prior(
    spot_idx: int,
    adata: sc.AnnData,
    cell_type_numbers_array: np.ndarray,
    radius: float,
    global_prior: Optional[np.ndarray] = None,
    lambda_prior_weight: float = 0.0,
    local_enrichment_weight: float = 0.5,
    global_enrichment_weight: float = 0.5,
    continuous_relaxation: bool = True,
    lambda_gex_reg: float = 0.01,
    # NEW parameters
    anchor_genes: Optional[Dict[int, List[int]]] = None,
    module_weight: float = 0.5,
    use_kl_regularization: bool = False,
    kl_temperature: float = 0.3,
    lambda_kl: float = 0.1,
) -> Optional[np.ndarray]:
```

### Step 5.3: Replace L2 objective with KL (in objective building section ~line 2400)

Replace the objective building loop (lines ~2400-2425):

```python
        # Compute module-aware enrichment if anchors provided
        if anchor_genes is not None and module_weight > 0:
            gene_specific_enrichment = compute_module_aware_enrichment(
                spot_idx=spot_idx,
                neighborhood_expression=neighborhood_expression_data,
                base_enrichment=gene_specific_enrichment,
                anchor_genes=anchor_genes,
                module_weight=module_weight,
            )

        # Build objective
        obj_terms = []
        for k in range(M):
            total_counts = int(center_counts[k])
            if total_counts > 0:
                enrichment_for_gene = gene_specific_enrichment[k, :]

                if use_kl_regularization:
                    # NEW: Softmax KL-divergence regularization
                    target = compute_softmax_target(enrichment_for_gene, temperature=kl_temperature)
                    kl_coeffs = compute_kl_penalty_coefficients(target, total_counts, lambda_kl)

                    for j in range(T):
                        # Enrichment term
                        base_term = enrichment_for_gene[j] * center_props[j] * X[j, k]
                        obj_terms.append(base_term)

                        # KL penalty (pulls toward target)
                        target_j = kl_coeffs['target_counts'][j]
                        penalty = kl_coeffs['penalty_weight'] * (X[j, k] - target_j) * (X[j, k] - target_j)
                        obj_terms.append(-penalty)

                        # Prior terms (unchanged)
                        if global_prior is not None and lambda_prior_weight > 0:
                            try:
                                prior_value = float(global_prior[j, k])
                                prior_penalty = lambda_prior_weight * (1 - prior_value) * X[j, k]
                                obj_terms.append(-prior_penalty)
                            except Exception as e:
                                logging.warning(f"Error accessing prior at [{j}, {k}]: {e}")
                                continue
                else:
                    # OLD: L2 regularization (backward compatible)
                    for j in range(T):
                        enrichment_weight = enrichment_for_gene[j]
                        cell_type_weight = center_props[j]

                        base_term = enrichment_weight * cell_type_weight * X[j, k]
                        obj_terms.append(base_term)

                        if lambda_gex_reg > 0:
                            obj_terms.append(-lambda_gex_reg * X[j, k] * X[j, k])

                        if global_prior is not None and lambda_prior_weight > 0:
                            try:
                                prior_value = float(global_prior[j, k])
                                prior_penalty = lambda_prior_weight * (1 - prior_value) * X[j, k]
                                obj_terms.append(-prior_penalty)
                            except Exception as e:
                                logging.warning(f"Error accessing prior at [{j}, {k}]: {e}")
                                continue
```

### Step 5.4: Commit

```bash
git add CITEgeist/model/gurobi_impl.py
git commit -m "feat: integrate module enrichment and KL regularization into GEX deconv"
```

---

## Task 6: Modify run_cell_expression_pass1 in citegeist_model.py

**Files:**
- Modify: `CITEgeist/model/citegeist_model.py:1352-1450`

### Step 6.1: Add new parameters to function signature

```python
def run_cell_expression_pass1(
    self,
    radius=None,
    alpha=0.5,
    global_enrichment_weight=0.5,
    local_enrichment_weight=0.5,
    max_workers=None,
    checkpoint_interval=100,
    output_dir="checkpoints",
    rerun=True,
    continuous_relaxation=True,
    lambda_gex_reg=0.01,
    cell_counts: Optional[pd.DataFrame] = None,
    use_discrete_mode: bool = False,
    # NEW parameters
    use_module_enrichment: bool = True,
    module_weight: float = 0.5,
    use_kl_regularization: bool = True,
    kl_temperature: float = 0.3,
    lambda_kl: float = 0.1,
    n_anchor_genes: Tuple[int, int] = (5, 10),
):
```

### Step 6.2: Add anchor discovery before optimization loop

Add after data preparation, before calling optimize_gene_expression:

```python
        # Discover anchor genes if using module enrichment
        anchor_genes = None
        if use_module_enrichment:
            from .gex_modules import discover_anchor_genes

            gene_expr_matrix = self.gene_expression_adata.X
            if scipy.sparse.issparse(gene_expr_matrix):
                gene_expr_matrix = gene_expr_matrix.toarray()

            anchor_genes, anchor_thresholds = discover_anchor_genes(
                gene_expression=gene_expr_matrix,
                cell_proportions=self.cell_prop_results.values,
                min_anchors=n_anchor_genes[0],
                max_anchors=n_anchor_genes[1],
            )

            total_anchors = sum(len(v) for v in anchor_genes.values())
            logging.info(f"Discovered {total_anchors} anchor genes across {len(anchor_genes)} cell types")
            for t, genes in anchor_genes.items():
                ct_name = self.cell_type_names[t] if hasattr(self, 'cell_type_names') else f"Type_{t}"
                logging.info(f"  {ct_name}: {len(genes)} anchors (r>={anchor_thresholds[t]:.2f})")
```

### Step 6.3: Pass new parameters to optimize_gene_expression

Update the call to `optimize_gene_expression`:

```python
        results = optimize_gene_expression(
            sample_name=self.sample_name,
            deconvolution_expression_data=gene_expression_data,
            cell_type_numbers_array=cell_type_proportions,
            filtered_adata=self.gene_expression_adata,
            radius=radius,
            global_enrichment_weight=global_enrichment_weight,
            local_enrichment_weight=local_enrichment_weight,
            global_prior=global_prior,
            lambda_prior_weight=lambda_prior_weight,
            max_workers=max_workers,
            checkpoint_interval=checkpoint_interval,
            output_dir=output_dir,
            rerun=rerun,
            continuous_relaxation=continuous_relaxation,
            lambda_gex_reg=lambda_gex_reg if not use_kl_regularization else 0.0,
            # NEW parameters
            anchor_genes=anchor_genes,
            module_weight=module_weight if use_module_enrichment else 0.0,
            use_kl_regularization=use_kl_regularization,
            kl_temperature=kl_temperature,
            lambda_kl=lambda_kl,
        )
```

### Step 6.4: Commit

```bash
git add CITEgeist/model/citegeist_model.py
git commit -m "feat: add module enrichment and KL params to run_cell_expression_pass1"
```

---

## Task 7: Update optimize_gene_expression signature

**Files:**
- Modify: `CITEgeist/model/gurobi_impl.py:2507-2520`

### Step 7.1: Add new parameters

```python
def optimize_gene_expression(
    sample_name: str,
    deconvolution_expression_data: np.ndarray,
    cell_type_numbers_array: np.ndarray,
    filtered_adata: sc.AnnData,
    radius: float = 2,
    global_enrichment_weight: float = 0.5,
    local_enrichment_weight: float = 0.5,
    global_prior: Optional[np.ndarray] = None,
    lambda_prior_weight: float = 0.0,
    max_workers: Optional[int] = None,
    checkpoint_interval: int = 100,
    output_dir: str = "checkpoints",
    rerun: bool = False,
    continuous_relaxation: bool = True,
    lambda_gex_reg: float = 0.01,
    # NEW parameters
    anchor_genes: Optional[Dict[int, List[int]]] = None,
    module_weight: float = 0.5,
    use_kl_regularization: bool = False,
    kl_temperature: float = 0.3,
    lambda_kl: float = 0.1,
) -> Dict[str, Any]:
```

### Step 7.2: Pass parameters to deconvolute_spot_with_neighbors_with_prior

Update the executor.submit call (~line 2601):

```python
                        future = executor.submit(
                            deconvolute_spot_with_neighbors_with_prior,
                            spot_idx,
                            filtered_adata,
                            cell_type_numbers_array,
                            radius,
                            global_prior,
                            lambda_prior_weight,
                            local_enrichment_weight,
                            global_enrichment_weight,
                            continuous_relaxation,
                            lambda_gex_reg,
                            # NEW parameters
                            anchor_genes,
                            module_weight,
                            use_kl_regularization,
                            kl_temperature,
                            lambda_kl,
                        )
```

### Step 7.3: Commit

```bash
git add CITEgeist/model/gurobi_impl.py
git commit -m "feat: thread module/KL params through optimize_gene_expression"
```

---

## Task 8: Add exports to __init__.py

**Files:**
- Modify: `CITEgeist/model/__init__.py`

### Step 8.1: Add exports

```python
from .gex_modules import (
    discover_anchor_genes,
    compute_module_aware_enrichment,
    compute_softmax_target,
    compute_kl_penalty_coefficients,
)
```

### Step 8.2: Commit

```bash
git add CITEgeist/model/__init__.py
git commit -m "feat: export gex_modules functions from model package"
```

---

## Task 9: Create benchmark SLURM script

**Files:**
- Create: `Benchmarking/xenium_benchmarking/CITEgeist/slurm/sbatch_gex_module_kl.sh`

### Step 9.1: Write SLURM script

```bash
#!/bin/bash
#SBATCH --job-name=gex_module_kl
#SBATCH --output=logs/gex_module_kl_%A_%a.out
#SBATCH --error=logs/gex_module_kl_%A_%a.err
#SBATCH --array=0-4
#SBATCH --time=06:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --partition=htc
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

module load gurobi/12.0.3

eval "$(conda shell.bash hook)"
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist

# Module-aware enrichment + KL regularization
python Benchmarking/xenium_benchmarking/CITEgeist/src/benchmark_hybrid_cellpose.py \
    --region ${SLURM_ARRAY_TASK_ID} \
    --output-dir Benchmarking/xenium_benchmarking/CITEgeist/output/gex_module_kl \
    --lambda-laplacian 0.1 \
    --spot-diameter-um 55.0 \
    --use-module-enrichment \
    --module-weight 0.5 \
    --use-kl-regularization \
    --kl-temperature 0.3 \
    --lambda-kl 0.1
```

### Step 9.2: Commit

```bash
git add Benchmarking/xenium_benchmarking/CITEgeist/slurm/sbatch_gex_module_kl.sh
git commit -m "feat: add benchmark script for module+KL GEX deconvolution"
```

---

## Task 10: Add CLI arguments to benchmark script

**Files:**
- Modify: `Benchmarking/xenium_benchmarking/CITEgeist/src/benchmark_hybrid_cellpose.py`

### Step 10.1: Add argparse arguments

Add after existing arguments (~line 444):

```python
    # Module enrichment and KL regularization
    parser.add_argument("--use-module-enrichment", action="store_true", default=False,
                        help="Use module-aware enrichment based on anchor genes")
    parser.add_argument("--module-weight", type=float, default=0.5,
                        help="Weight for module signal vs base enrichment (0-1)")
    parser.add_argument("--use-kl-regularization", action="store_true", default=False,
                        help="Use KL-divergence regularization instead of L2")
    parser.add_argument("--kl-temperature", type=float, default=0.3,
                        help="Softmax temperature for KL target")
    parser.add_argument("--lambda-kl", type=float, default=0.1,
                        help="KL penalty weight")
```

### Step 10.2: Pass to model.run_cell_expression_pass1

Update the call in run_hybrid_benchmark:

```python
    model.run_cell_expression_pass1(
        ...,
        use_module_enrichment=args.use_module_enrichment,
        module_weight=args.module_weight,
        use_kl_regularization=args.use_kl_regularization,
        kl_temperature=args.kl_temperature,
        lambda_kl=args.lambda_kl,
    )
```

### Step 10.3: Commit

```bash
git add Benchmarking/xenium_benchmarking/CITEgeist/src/benchmark_hybrid_cellpose.py
git commit -m "feat: add module/KL CLI arguments to benchmark script"
```

---

## Task 11: Run ablation benchmarks

**Files:**
- Run scripts, evaluate results

### Step 11.1: Submit baseline (control)

```bash
cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist
sbatch Benchmarking/xenium_benchmarking/CITEgeist/slurm/sbatch_hybrid_cellpose.sh
```

### Step 11.2: Submit KL-only variant

Create and submit:
```bash
# sbatch_gex_kl_only.sh with --use-kl-regularization but NOT --use-module-enrichment
```

### Step 11.3: Submit Module-only variant

Create and submit:
```bash
# sbatch_gex_module_only.sh with --use-module-enrichment but NOT --use-kl-regularization
```

### Step 11.4: Submit full variant

```bash
sbatch Benchmarking/xenium_benchmarking/CITEgeist/slurm/sbatch_gex_module_kl.sh
```

### Step 11.5: Evaluate all variants

```bash
for variant in hybrid_cellpose gex_kl_only gex_module_only gex_module_kl; do
    python Benchmarking/xenium_benchmarking/evaluation/src/evaluate_gex_spatial.py \
        --output-subdir $variant \
        --method-name $variant
done
```

### Step 11.6: Commit evaluation results

```bash
git add Benchmarking/xenium_benchmarking/CITEgeist/output/*/spatial_gex_summary.json
git commit -m "results: GEX module+KL ablation benchmark results"
```

---

## Summary

| Task | Description | Files |
|------|-------------|-------|
| 1 | Anchor gene discovery | gex_modules.py, test_gex_modules.py |
| 2 | Module-aware enrichment | gex_modules.py, test_gex_modules.py |
| 3 | Softmax KL functions | gex_modules.py, test_gex_modules.py |
| 4 | Gurobi integration test | test_gex_modules.py |
| 5 | Modify deconvolute_spot | gurobi_impl.py |
| 6 | Modify run_cell_expression_pass1 | citegeist_model.py |
| 7 | Update optimize_gene_expression | gurobi_impl.py |
| 8 | Add exports | __init__.py |
| 9 | Create benchmark SLURM | sbatch_gex_module_kl.sh |
| 10 | Add CLI arguments | benchmark_hybrid_cellpose.py |
| 11 | Run ablation benchmarks | Evaluation |

**Success criteria:**
- Per-gene spatial r > 0.30 (baseline: 0.273)
- Close marker-general gap
- Variance ratio closer to 1.0 (baseline: 3.02)
