# Module 3 Multimodal Refinement Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add RNA-based refinement to Module 3 that learns gene signatures from protein-confident spots and iteratively refines cell type proportions using both protein priors and RNA reconstruction.

**Architecture:** Pass 1.5 learns anchor genes (top-20 per cell type by Pearson correlation). Pass 2 EM iterates between E-step (estimate expression profiles with anchors locked) and M-step (refine proportions with protein prior). This closes the gap for cells with RNA signal but low protein expression.

**Tech Stack:** Python 3.10, NumPy, SciPy (pearsonr, optimize), Gurobi (QP for M-step), existing CITEgeist infrastructure.

---

## Task 1: Create Test File and First Failing Test

**Files:**
- Create: `CITEgeist/tests/test_multimodal_refinement.py`

**Step 1: Write the failing test for anchor gene selection**

```python
"""Tests for multimodal refinement (Pass 1.5 + Pass 2 EM)."""

import numpy as np
import pytest
from scipy.stats import pearsonr


class TestSelectAnchorGenes:
    """Test anchor gene selection from protein-confident spots."""

    def test_select_anchor_genes_returns_correct_structure(self):
        """Anchor selection returns dict of gene lists and weights."""
        # Import will fail until we create the module
        from CITEgeist.model.multimodal_refinement import select_anchor_genes

        # Synthetic data: 100 spots, 50 genes, 3 cell types
        np.random.seed(42)
        n_spots, n_genes, n_types = 100, 50, 3

        # Create Y_protein with clear cell type assignments
        Y_protein = np.zeros((n_spots, n_types))
        Y_protein[:33, 0] = 0.8  # First 33 spots are type 0
        Y_protein[33:66, 1] = 0.8  # Next 33 spots are type 1
        Y_protein[66:, 2] = 0.8  # Last 34 spots are type 2

        # Create GEX where first 5 genes correlate with each type
        GEX = np.random.randn(n_spots, n_genes) * 0.1
        GEX[:33, 0:5] += 2.0  # Genes 0-4 high in type 0 spots
        GEX[33:66, 5:10] += 2.0  # Genes 5-9 high in type 1 spots
        GEX[66:, 10:15] += 2.0  # Genes 10-14 high in type 2 spots

        gene_names = [f"Gene_{i}" for i in range(n_genes)]
        cell_type_names = ["TypeA", "TypeB", "TypeC"]

        anchors, weights = select_anchor_genes(
            GEX=GEX,
            Y_protein=Y_protein,
            gene_names=gene_names,
            cell_type_names=cell_type_names,
            n_anchors=5,
            min_correlation=0.3,
        )

        # Check structure
        assert isinstance(anchors, dict)
        assert isinstance(weights, dict)
        assert set(anchors.keys()) == set(cell_type_names)
        assert set(weights.keys()) == set(cell_type_names)

        # Check each type has up to n_anchors genes
        for ct in cell_type_names:
            assert len(anchors[ct]) <= 5
            assert len(weights[ct]) == len(anchors[ct])

        # Check correct genes selected (genes 0-4 for TypeA, etc.)
        assert "Gene_0" in anchors["TypeA"] or "Gene_1" in anchors["TypeA"]
        assert "Gene_5" in anchors["TypeB"] or "Gene_6" in anchors["TypeB"]
        assert "Gene_10" in anchors["TypeC"] or "Gene_11" in anchors["TypeC"]
```

**Step 2: Run test to verify it fails**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest CITEgeist/tests/test_multimodal_refinement.py::TestSelectAnchorGenes::test_select_anchor_genes_returns_correct_structure -v`

Expected: FAIL with `ModuleNotFoundError: No module named 'CITEgeist.model.multimodal_refinement'`

**Step 3: Commit test file**

```bash
git add CITEgeist/tests/test_multimodal_refinement.py
git commit -m "test: add failing test for anchor gene selection"
```

---

## Task 2: Implement select_anchor_genes Function

**Files:**
- Create: `CITEgeist/model/multimodal_refinement.py`

**Step 1: Create module with select_anchor_genes**

```python
"""
Multimodal refinement for CITEgeist Module 3.

This module implements Pass 1.5 (anchor gene learning) and Pass 2 EM
(RNA-based proportion refinement) to handle cells with RNA signal
but low protein expression.
"""

import logging
from typing import Dict, List, Optional, Tuple

import numpy as np
from scipy.stats import pearsonr

logger = logging.getLogger(__name__)


def select_anchor_genes(
    GEX: np.ndarray,
    Y_protein: np.ndarray,
    gene_names: List[str],
    cell_type_names: List[str],
    n_anchors: int = 20,
    min_correlation: float = 0.3,
) -> Tuple[Dict[str, List[str]], Dict[str, Dict[str, float]]]:
    """
    Select top-N anchor genes per cell type based on correlation and specificity.

    Pass 1.5: Learn gene signatures from protein-confident spots.

    Args:
        GEX: Gene expression matrix (N_spots × G_genes)
        Y_protein: Cell type proportions from Pass 1 (N_spots × T_types)
        gene_names: List of gene names (length G)
        cell_type_names: List of cell type names (length T)
        n_anchors: Number of anchor genes per cell type (default: 20)
        min_correlation: Minimum Pearson r to be considered (default: 0.3)

    Returns:
        anchors: Dict[cell_type] -> List[gene_names]
        weights: Dict[cell_type] -> Dict[gene_name -> correlation_weight]
    """
    n_spots, n_genes = GEX.shape
    n_types = Y_protein.shape[1]

    if n_genes != len(gene_names):
        raise ValueError(f"GEX has {n_genes} genes but {len(gene_names)} names provided")
    if n_types != len(cell_type_names):
        raise ValueError(f"Y_protein has {n_types} types but {len(cell_type_names)} names provided")

    # Compute correlation matrix (genes × cell types)
    correlations = np.zeros((n_genes, n_types))
    for t in range(n_types):
        y_t = Y_protein[:, t]
        # Skip if no variance in proportions
        if np.std(y_t) < 1e-10:
            continue
        for g in range(n_genes):
            gex_g = GEX[:, g]
            if np.std(gex_g) < 1e-10:
                continue
            r, _ = pearsonr(gex_g, y_t)
            if not np.isnan(r):
                correlations[g, t] = r

    logger.info(f"Computed correlations for {n_genes} genes × {n_types} cell types")

    # Compute specificity: r[g,t] - max(r[g, other_types])
    specificity = np.zeros((n_genes, n_types))
    for g in range(n_genes):
        for t in range(n_types):
            other_corrs = np.delete(correlations[g, :], t)
            other_max = np.max(other_corrs) if len(other_corrs) > 0 else 0
            specificity[g, t] = correlations[g, t] - other_max

    # Combined score: correlation × specificity (only positive specificity)
    score = correlations * np.clip(specificity, 0, None)

    # Select top-N per cell type
    anchors = {}
    weights = {}

    for t, ct_name in enumerate(cell_type_names):
        # Filter by minimum correlation
        valid_mask = correlations[:, t] >= min_correlation
        valid_indices = np.where(valid_mask)[0]

        if len(valid_indices) == 0:
            logger.warning(f"No genes pass min_correlation={min_correlation} for {ct_name}")
            anchors[ct_name] = []
            weights[ct_name] = {}
            continue

        # Rank by score among valid genes
        valid_scores = score[valid_indices, t]
        ranked_order = np.argsort(valid_scores)[::-1]  # Descending
        top_indices = valid_indices[ranked_order[:n_anchors]]

        # Build output
        anchors[ct_name] = [gene_names[i] for i in top_indices]
        weights[ct_name] = {
            gene_names[i]: float(correlations[i, t]) for i in top_indices
        }

        logger.info(
            f"{ct_name}: selected {len(anchors[ct_name])} anchors, "
            f"top gene={anchors[ct_name][0] if anchors[ct_name] else 'None'}"
        )

    return anchors, weights
```

**Step 2: Run test to verify it passes**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest CITEgeist/tests/test_multimodal_refinement.py::TestSelectAnchorGenes::test_select_anchor_genes_returns_correct_structure -v`

Expected: PASS

**Step 3: Commit**

```bash
git add CITEgeist/model/multimodal_refinement.py
git commit -m "feat: implement select_anchor_genes for Pass 1.5"
```

---

## Task 3: Add Test for E-Step (Expression Profile Estimation)

**Files:**
- Modify: `CITEgeist/tests/test_multimodal_refinement.py`

**Step 1: Add failing test for E-step**

```python
class TestComputeExpressionProfiles:
    """Test E-step: estimate expression profiles given proportions."""

    def test_compute_expression_profiles_anchor_locked(self):
        """Anchor genes should be assigned to their cell type only."""
        from CITEgeist.model.multimodal_refinement import compute_expression_profiles

        np.random.seed(42)
        n_spots, n_genes, n_types = 50, 20, 3

        # Y: proportions
        Y = np.random.dirichlet([1, 1, 1], size=n_spots)

        # GEX: random expression
        GEX = np.abs(np.random.randn(n_spots, n_genes))

        gene_names = [f"Gene_{i}" for i in range(n_genes)]
        cell_type_names = ["TypeA", "TypeB", "TypeC"]

        # Anchors: Gene_0,1 for TypeA, Gene_2,3 for TypeB, Gene_4,5 for TypeC
        anchors = {
            "TypeA": ["Gene_0", "Gene_1"],
            "TypeB": ["Gene_2", "Gene_3"],
            "TypeC": ["Gene_4", "Gene_5"],
        }
        weights = {
            "TypeA": {"Gene_0": 0.8, "Gene_1": 0.7},
            "TypeB": {"Gene_2": 0.75, "Gene_3": 0.65},
            "TypeC": {"Gene_4": 0.9, "Gene_5": 0.6},
        }

        E = compute_expression_profiles(
            GEX=GEX,
            Y=Y,
            gene_names=gene_names,
            cell_type_names=cell_type_names,
            anchors=anchors,
            weights=weights,
        )

        # E should be (n_types × n_genes)
        assert E.shape == (n_types, n_genes)

        # Anchor genes should have non-zero expression ONLY for their cell type
        # Gene_0 is anchor for TypeA (index 0)
        gene_0_idx = 0
        assert E[0, gene_0_idx] > 0  # TypeA has expression
        # For locked anchors, other types should be 0
        assert E[1, gene_0_idx] == 0  # TypeB should be 0
        assert E[2, gene_0_idx] == 0  # TypeC should be 0
```

**Step 2: Run test to verify it fails**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest CITEgeist/tests/test_multimodal_refinement.py::TestComputeExpressionProfiles -v`

Expected: FAIL with `ImportError` (function not defined)

**Step 3: Commit**

```bash
git add CITEgeist/tests/test_multimodal_refinement.py
git commit -m "test: add failing test for E-step expression profiles"
```

---

## Task 4: Implement compute_expression_profiles (E-Step)

**Files:**
- Modify: `CITEgeist/model/multimodal_refinement.py`

**Step 1: Add compute_expression_profiles function**

```python
def compute_expression_profiles(
    GEX: np.ndarray,
    Y: np.ndarray,
    gene_names: List[str],
    cell_type_names: List[str],
    anchors: Dict[str, List[str]],
    weights: Dict[str, Dict[str, float]],
) -> np.ndarray:
    """
    E-step: Estimate gene expression profiles per cell type.

    Anchor genes are LOCKED to their assigned cell type.
    Non-anchor genes are solved via least squares (can load on any type).

    Args:
        GEX: Gene expression matrix (N_spots × G_genes)
        Y: Current cell type proportions (N_spots × T_types)
        gene_names: List of gene names (length G)
        cell_type_names: List of cell type names (length T)
        anchors: Dict[cell_type] -> List[anchor_gene_names]
        weights: Dict[cell_type] -> Dict[gene_name -> weight]

    Returns:
        E: Expression profiles (T_types × G_genes)
    """
    n_spots, n_genes = GEX.shape
    n_types = len(cell_type_names)

    # Build gene name to index mapping
    gene_to_idx = {g: i for i, g in enumerate(gene_names)}
    type_to_idx = {t: i for i, t in enumerate(cell_type_names)}

    # Build set of anchor genes and their assignments
    anchor_assignments = {}  # gene_name -> cell_type_name
    for ct_name, gene_list in anchors.items():
        for g in gene_list:
            anchor_assignments[g] = ct_name

    # Initialize E
    E = np.zeros((n_types, n_genes))

    for g_idx, g_name in enumerate(gene_names):
        gex_g = GEX[:, g_idx]  # Expression of gene g across spots

        if g_name in anchor_assignments:
            # LOCKED: anchor gene assigned to one cell type only
            ct_name = anchor_assignments[g_name]
            t_idx = type_to_idx[ct_name]

            # Weighted mean of expression where cell type is present
            y_t = Y[:, t_idx]
            if np.sum(y_t) > 1e-10:
                E[t_idx, g_idx] = np.sum(gex_g * y_t) / np.sum(y_t)
            else:
                E[t_idx, g_idx] = 0.0

            # Other cell types get 0 for this anchor gene
            # (already initialized to 0)

        else:
            # FREE: non-anchor gene can load on any cell type
            # Solve least squares: GEX[:, g] ≈ Y @ E[:, g]
            # E[:, g] = (Y^T Y)^{-1} Y^T GEX[:, g]
            try:
                E[:, g_idx], _, _, _ = np.linalg.lstsq(Y, gex_g, rcond=None)
                # Clip negative values (expression should be non-negative)
                E[:, g_idx] = np.clip(E[:, g_idx], 0, None)
            except np.linalg.LinAlgError:
                # Fallback: uniform assignment
                E[:, g_idx] = np.mean(gex_g) / n_types

    return E
```

**Step 2: Run test to verify it passes**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest CITEgeist/tests/test_multimodal_refinement.py::TestComputeExpressionProfiles -v`

Expected: PASS

**Step 3: Commit**

```bash
git add CITEgeist/model/multimodal_refinement.py
git commit -m "feat: implement compute_expression_profiles (E-step)"
```

---

## Task 5: Add Test for M-Step (Proportion Refinement)

**Files:**
- Modify: `CITEgeist/tests/test_multimodal_refinement.py`

**Step 1: Add failing test for M-step**

```python
class TestRefineProportions:
    """Test M-step: refine proportions given expression profiles."""

    def test_refine_proportions_respects_constraints(self):
        """Refined proportions should be non-negative and sum to <= 1."""
        from CITEgeist.model.multimodal_refinement import refine_proportions

        np.random.seed(42)
        n_spots, n_genes, n_types = 30, 15, 3

        # Y_protein: initial proportions from Pass 1
        Y_protein = np.random.dirichlet([2, 2, 2], size=n_spots)

        # Y_current: current estimate
        Y_current = Y_protein.copy()

        # E: expression profiles
        E = np.abs(np.random.randn(n_types, n_genes))

        # GEX: observed expression (with some noise)
        GEX = Y_protein @ E + np.random.randn(n_spots, n_genes) * 0.1

        gene_names = [f"Gene_{i}" for i in range(n_genes)]
        cell_type_names = ["TypeA", "TypeB", "TypeC"]

        anchors = {"TypeA": ["Gene_0"], "TypeB": ["Gene_1"], "TypeC": ["Gene_2"]}
        weights = {"TypeA": {"Gene_0": 0.8}, "TypeB": {"Gene_1": 0.7}, "TypeC": {"Gene_2": 0.9}}

        Y_refined = refine_proportions(
            GEX=GEX,
            Y_current=Y_current,
            E=E,
            Y_protein=Y_protein,
            gene_names=gene_names,
            cell_type_names=cell_type_names,
            anchors=anchors,
            weights=weights,
            lambda_prior=1.0,
        )

        # Check shape
        assert Y_refined.shape == (n_spots, n_types)

        # Check non-negativity
        assert np.all(Y_refined >= -1e-10)

        # Check sum <= 1 (with small tolerance)
        assert np.all(Y_refined.sum(axis=1) <= 1.0 + 1e-6)

    def test_refine_proportions_stays_near_prior_with_high_lambda(self):
        """With high lambda_prior, Y_refined should stay close to Y_protein."""
        from CITEgeist.model.multimodal_refinement import refine_proportions

        np.random.seed(42)
        n_spots, n_genes, n_types = 20, 10, 3

        Y_protein = np.random.dirichlet([2, 2, 2], size=n_spots)
        Y_current = Y_protein.copy()
        E = np.abs(np.random.randn(n_types, n_genes))
        GEX = np.abs(np.random.randn(n_spots, n_genes))  # Random, doesn't match

        gene_names = [f"Gene_{i}" for i in range(n_genes)]
        cell_type_names = ["TypeA", "TypeB", "TypeC"]
        anchors = {"TypeA": [], "TypeB": [], "TypeC": []}
        weights = {"TypeA": {}, "TypeB": {}, "TypeC": {}}

        Y_refined = refine_proportions(
            GEX=GEX,
            Y_current=Y_current,
            E=E,
            Y_protein=Y_protein,
            gene_names=gene_names,
            cell_type_names=cell_type_names,
            anchors=anchors,
            weights=weights,
            lambda_prior=100.0,  # Very high prior
        )

        # Should stay very close to Y_protein
        assert np.allclose(Y_refined, Y_protein, atol=0.1)
```

**Step 2: Run test to verify it fails**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest CITEgeist/tests/test_multimodal_refinement.py::TestRefineProportions -v`

Expected: FAIL with `ImportError`

**Step 3: Commit**

```bash
git add CITEgeist/tests/test_multimodal_refinement.py
git commit -m "test: add failing tests for M-step proportion refinement"
```

---

## Task 6: Implement refine_proportions (M-Step)

**Files:**
- Modify: `CITEgeist/model/multimodal_refinement.py`

**Step 1: Add refine_proportions function**

```python
def refine_proportions(
    GEX: np.ndarray,
    Y_current: np.ndarray,
    E: np.ndarray,
    Y_protein: np.ndarray,
    gene_names: List[str],
    cell_type_names: List[str],
    anchors: Dict[str, List[str]],
    weights: Dict[str, Dict[str, float]],
    lambda_prior: float = 1.0,
) -> np.ndarray:
    """
    M-step: Refine proportions Y given expression profiles E.

    Objective per spot i:
        minimize: Σ_g w[g] * (GEX[i,g] - Y[i,:] @ E[:,g])²
                  + λ * ||Y[i,:] - Y_protein[i,:]||²

    Constraints:
        Y[i, :] >= 0
        sum(Y[i, :]) <= 1

    Args:
        GEX: Gene expression matrix (N_spots × G_genes)
        Y_current: Current proportions (N_spots × T_types)
        E: Expression profiles (T_types × G_genes)
        Y_protein: Original protein-based proportions (N_spots × T_types)
        gene_names: List of gene names
        cell_type_names: List of cell type names
        anchors: Dict of anchor genes per cell type
        weights: Dict of weights per anchor gene
        lambda_prior: Weight of protein prior (default: 1.0)

    Returns:
        Y_refined: Updated proportions (N_spots × T_types)
    """
    from scipy.optimize import minimize, Bounds

    n_spots, n_genes = GEX.shape
    n_types = len(cell_type_names)

    # Build gene weights vector
    gene_to_idx = {g: i for i, g in enumerate(gene_names)}
    gene_weights = np.ones(n_genes)
    for ct_name, gene_dict in weights.items():
        for g_name, w in gene_dict.items():
            if g_name in gene_to_idx:
                gene_weights[gene_to_idx[g_name]] = w

    Y_refined = np.zeros((n_spots, n_types))

    for i in range(n_spots):
        gex_i = GEX[i, :]
        y_protein_i = Y_protein[i, :]

        def objective(y):
            # Reconstruction error (weighted)
            reconstruction = E.T @ y  # (G,)
            recon_error = np.sum(gene_weights * (gex_i - reconstruction) ** 2)

            # Prior toward protein proportions
            prior_error = lambda_prior * np.sum((y - y_protein_i) ** 2)

            return recon_error + prior_error

        # Bounds: 0 <= y[t] <= 1
        bounds = Bounds(lb=np.zeros(n_types), ub=np.ones(n_types))

        # Constraint: sum(y) <= 1
        constraints = {"type": "ineq", "fun": lambda y: 1.0 - np.sum(y)}

        # Initial guess: current proportions
        y0 = Y_current[i, :].copy()

        result = minimize(
            objective,
            y0,
            method="SLSQP",
            bounds=bounds,
            constraints=constraints,
            options={"maxiter": 100, "ftol": 1e-8},
        )

        if result.success:
            Y_refined[i, :] = result.x
        else:
            # Fallback to current
            Y_refined[i, :] = Y_current[i, :]

    return Y_refined
```

**Step 2: Run test to verify it passes**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest CITEgeist/tests/test_multimodal_refinement.py::TestRefineProportions -v`

Expected: PASS

**Step 3: Commit**

```bash
git add CITEgeist/model/multimodal_refinement.py
git commit -m "feat: implement refine_proportions (M-step)"
```

---

## Task 7: Add Test for Full EM Loop

**Files:**
- Modify: `CITEgeist/tests/test_multimodal_refinement.py`

**Step 1: Add failing test for full EM**

```python
class TestMultimodalEMRefinement:
    """Test full Pass 1.5 + Pass 2 EM refinement."""

    def test_multimodal_em_improves_reconstruction(self):
        """EM refinement should improve GEX reconstruction error."""
        from CITEgeist.model.multimodal_refinement import multimodal_em_refinement

        np.random.seed(42)
        n_spots, n_genes, n_types = 50, 30, 3

        # True proportions
        Y_true = np.random.dirichlet([2, 2, 2], size=n_spots)

        # True expression profiles
        E_true = np.abs(np.random.randn(n_types, n_genes)) + 0.5

        # Observed GEX
        GEX = Y_true @ E_true + np.random.randn(n_spots, n_genes) * 0.1

        # Y_protein: noisy version of Y_true (simulating protein estimation error)
        Y_protein = Y_true + np.random.randn(n_spots, n_types) * 0.15
        Y_protein = np.clip(Y_protein, 0, 1)
        Y_protein = Y_protein / Y_protein.sum(axis=1, keepdims=True)

        gene_names = [f"Gene_{i}" for i in range(n_genes)]
        cell_type_names = ["TypeA", "TypeB", "TypeC"]

        Y_refined, E_final, anchors = multimodal_em_refinement(
            GEX=GEX,
            Y_protein=Y_protein,
            gene_names=gene_names,
            cell_type_names=cell_type_names,
            n_anchors=5,
            lambda_prior=1.0,
            max_iterations=10,
            tolerance=1e-4,
        )

        # Check shapes
        assert Y_refined.shape == Y_protein.shape
        assert E_final.shape == (n_types, n_genes)
        assert isinstance(anchors, dict)

        # Reconstruction error should be lower with Y_refined than Y_protein
        recon_protein = np.sum((GEX - Y_protein @ E_true) ** 2)
        recon_refined = np.sum((GEX - Y_refined @ E_final) ** 2)

        # Refined should have lower or equal reconstruction error
        assert recon_refined <= recon_protein * 1.1  # Allow 10% tolerance

    def test_multimodal_em_converges(self):
        """EM should converge within max_iterations."""
        from CITEgeist.model.multimodal_refinement import multimodal_em_refinement

        np.random.seed(123)
        n_spots, n_genes, n_types = 30, 20, 2

        Y_protein = np.random.dirichlet([3, 3], size=n_spots)
        GEX = np.abs(np.random.randn(n_spots, n_genes))

        gene_names = [f"Gene_{i}" for i in range(n_genes)]
        cell_type_names = ["TypeA", "TypeB"]

        # Should complete without error
        Y_refined, E_final, anchors = multimodal_em_refinement(
            GEX=GEX,
            Y_protein=Y_protein,
            gene_names=gene_names,
            cell_type_names=cell_type_names,
            n_anchors=3,
            lambda_prior=1.0,
            max_iterations=20,
            tolerance=1e-4,
        )

        assert Y_refined is not None
        assert not np.any(np.isnan(Y_refined))
```

**Step 2: Run test to verify it fails**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest CITEgeist/tests/test_multimodal_refinement.py::TestMultimodalEMRefinement -v`

Expected: FAIL with `ImportError`

**Step 3: Commit**

```bash
git add CITEgeist/tests/test_multimodal_refinement.py
git commit -m "test: add failing tests for full EM refinement"
```

---

## Task 8: Implement multimodal_em_refinement (Full EM Loop)

**Files:**
- Modify: `CITEgeist/model/multimodal_refinement.py`

**Step 1: Add multimodal_em_refinement function**

```python
def multimodal_em_refinement(
    GEX: np.ndarray,
    Y_protein: np.ndarray,
    gene_names: List[str],
    cell_type_names: List[str],
    n_anchors: int = 20,
    min_correlation: float = 0.3,
    lambda_prior: float = 1.0,
    max_iterations: int = 20,
    tolerance: float = 1e-4,
) -> Tuple[np.ndarray, np.ndarray, Dict[str, List[str]]]:
    """
    Full multimodal EM refinement: Pass 1.5 + Pass 2 EM.

    1. Pass 1.5: Learn anchor genes from protein-confident spots
    2. Pass 2 EM: Iterate E-step (expression profiles) and M-step (proportions)

    Args:
        GEX: Gene expression matrix (N_spots × G_genes)
        Y_protein: Cell type proportions from Pass 1 (N_spots × T_types)
        gene_names: List of gene names (length G)
        cell_type_names: List of cell type names (length T)
        n_anchors: Number of anchor genes per cell type (default: 20)
        min_correlation: Minimum Pearson r for anchor selection (default: 0.3)
        lambda_prior: Weight of protein prior in M-step (default: 1.0)
        max_iterations: Maximum EM iterations (default: 20)
        tolerance: Convergence tolerance (default: 1e-4)

    Returns:
        Y_refined: Refined proportions (N_spots × T_types)
        E_final: Final expression profiles (T_types × G_genes)
        anchors: Dict of anchor genes per cell type
    """
    n_spots, n_genes = GEX.shape
    n_types = len(cell_type_names)

    logger.info(f"Starting multimodal EM refinement: {n_spots} spots, {n_genes} genes, {n_types} types")

    # Pass 1.5: Learn anchor genes
    logger.info("Pass 1.5: Selecting anchor genes...")
    anchors, weights = select_anchor_genes(
        GEX=GEX,
        Y_protein=Y_protein,
        gene_names=gene_names,
        cell_type_names=cell_type_names,
        n_anchors=n_anchors,
        min_correlation=min_correlation,
    )

    total_anchors = sum(len(v) for v in anchors.values())
    logger.info(f"Selected {total_anchors} total anchor genes across {n_types} cell types")

    # Initialize
    Y = Y_protein.copy()
    E = None

    # EM loop
    for iteration in range(max_iterations):
        logger.info(f"EM iteration {iteration + 1}/{max_iterations}")

        # E-step: Estimate expression profiles
        E = compute_expression_profiles(
            GEX=GEX,
            Y=Y,
            gene_names=gene_names,
            cell_type_names=cell_type_names,
            anchors=anchors,
            weights=weights,
        )

        # M-step: Refine proportions
        Y_new = refine_proportions(
            GEX=GEX,
            Y_current=Y,
            E=E,
            Y_protein=Y_protein,
            gene_names=gene_names,
            cell_type_names=cell_type_names,
            anchors=anchors,
            weights=weights,
            lambda_prior=lambda_prior,
        )

        # Check convergence
        max_change = np.max(np.abs(Y_new - Y))
        logger.info(f"  Max proportion change: {max_change:.6f}")

        if max_change < tolerance:
            logger.info(f"Converged after {iteration + 1} iterations")
            Y = Y_new
            break

        Y = Y_new

    else:
        logger.warning(f"Did not converge within {max_iterations} iterations")

    # Final E-step to get consistent E
    E_final = compute_expression_profiles(
        GEX=GEX,
        Y=Y,
        gene_names=gene_names,
        cell_type_names=cell_type_names,
        anchors=anchors,
        weights=weights,
    )

    logger.info("Multimodal EM refinement complete")

    return Y, E_final, anchors
```

**Step 2: Run test to verify it passes**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest CITEgeist/tests/test_multimodal_refinement.py::TestMultimodalEMRefinement -v`

Expected: PASS

**Step 3: Commit**

```bash
git add CITEgeist/model/multimodal_refinement.py
git commit -m "feat: implement multimodal_em_refinement (full EM loop)"
```

---

## Task 9: Integrate with CitegeistModel

**Files:**
- Modify: `CITEgeist/model/citegeist_model.py`

**Step 1: Add import at top of file (after existing imports)**

Find the imports section (~line 20-45) and add:

```python
from .multimodal_refinement import multimodal_em_refinement
```

**Step 2: Add run_multimodal_refinement method**

Add this method to the CitegeistModel class (after `run_cell_proportion_model`):

```python
def run_multimodal_refinement(
    self,
    n_anchors: int = 20,
    min_correlation: float = 0.3,
    lambda_prior: float = 1.0,
    max_iterations: int = 20,
    tolerance: float = 1e-4,
) -> pd.DataFrame:
    """
    Run Pass 1.5 + Pass 2 EM multimodal refinement.

    Uses RNA expression to refine protein-based proportions, addressing
    cells with RNA signal but low protein expression.

    Prerequisites:
        - run_cell_proportion_model() must be called first (Pass 1)
        - GEX data must be preprocessed

    Args:
        n_anchors: Number of anchor genes per cell type (default: 20)
        min_correlation: Minimum Pearson r for anchor selection (default: 0.3)
        lambda_prior: Trust in protein prior vs RNA (default: 1.0)
        max_iterations: Maximum EM iterations (default: 20)
        tolerance: Convergence tolerance (default: 1e-4)

    Returns:
        DataFrame with refined cell type proportions
    """
    if not hasattr(self, "cell_prop_global_results") or self.cell_prop_global_results is None:
        raise ValueError("Must run run_cell_proportion_model() first (Pass 1)")

    if self.gene_expression_adata is None:
        raise ValueError("Gene expression data not loaded")

    self.logger.info("=" * 60)
    self.logger.info("Running Multimodal Refinement (Pass 1.5 + Pass 2 EM)")
    self.logger.info("=" * 60)

    # Get GEX matrix
    GEX = self.gene_expression_adata.X
    if hasattr(GEX, "toarray"):
        GEX = GEX.toarray()
    GEX = np.array(GEX, dtype=np.float64)

    gene_names = list(self.gene_expression_adata.var_names)
    spot_names = list(self.gene_expression_adata.obs_names)

    # Get Y_protein from Pass 1 results
    cell_type_names = [c for c in self.cell_prop_global_results.columns if c not in ["spot_id"]]
    Y_protein = self.cell_prop_global_results[cell_type_names].values

    self.logger.info(f"Input: {len(spot_names)} spots, {len(gene_names)} genes, {len(cell_type_names)} cell types")

    # Run multimodal EM
    Y_refined, E_final, anchors = multimodal_em_refinement(
        GEX=GEX,
        Y_protein=Y_protein,
        gene_names=gene_names,
        cell_type_names=cell_type_names,
        n_anchors=n_anchors,
        min_correlation=min_correlation,
        lambda_prior=lambda_prior,
        max_iterations=max_iterations,
        tolerance=tolerance,
    )

    # Store results
    self.cell_prop_refined_results = pd.DataFrame(
        Y_refined,
        index=spot_names,
        columns=cell_type_names,
    )
    self.expression_profiles = pd.DataFrame(
        E_final,
        index=cell_type_names,
        columns=gene_names,
    )
    self.anchor_genes = anchors

    # Log anchor gene summary
    self.logger.info("Anchor genes per cell type:")
    for ct, genes in anchors.items():
        self.logger.info(f"  {ct}: {genes[:5]}{'...' if len(genes) > 5 else ''}")

    # Save results
    output_path = self.output_folder / f"{self.sample_name}_cell_prop_refined_results.csv"
    self.cell_prop_refined_results.to_csv(output_path)
    self.logger.info(f"Saved refined proportions to {output_path}")

    anchors_path = self.output_folder / f"{self.sample_name}_anchor_genes.json"
    import json
    with open(anchors_path, "w") as f:
        json.dump(anchors, f, indent=2)
    self.logger.info(f"Saved anchor genes to {anchors_path}")

    return self.cell_prop_refined_results
```

**Step 3: Run all tests to verify nothing broke**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest CITEgeist/tests/test_multimodal_refinement.py -v`

Expected: All tests PASS

**Step 4: Commit**

```bash
git add CITEgeist/model/citegeist_model.py CITEgeist/model/multimodal_refinement.py
git commit -m "feat: integrate multimodal refinement with CitegeistModel"
```

---

## Task 10: Add Integration Test with Simulated Data

**Files:**
- Modify: `CITEgeist/tests/test_multimodal_refinement.py`

**Step 1: Add integration test**

```python
class TestIntegrationWithCitegeistModel:
    """Integration tests with CitegeistModel."""

    def test_run_multimodal_refinement_end_to_end(self, tmp_path):
        """Test full pipeline: Pass 1 -> Multimodal refinement."""
        import scanpy as sc
        from CITEgeist.model.citegeist_model import CitegeistModel

        np.random.seed(42)
        n_spots, n_genes, n_proteins, n_types = 50, 100, 10, 3

        # Create synthetic GEX AnnData
        gex_X = np.abs(np.random.randn(n_spots, n_genes))
        gex_adata = sc.AnnData(X=gex_X)
        gex_adata.obs_names = [f"spot_{i}" for i in range(n_spots)]
        gex_adata.var_names = [f"Gene_{i}" for i in range(n_genes)]
        gex_adata.obsm["spatial"] = np.random.rand(n_spots, 2) * 100

        # Create synthetic protein AnnData
        prot_X = np.abs(np.random.randn(n_spots, n_proteins))
        prot_adata = sc.AnnData(X=prot_X)
        prot_adata.obs_names = [f"spot_{i}" for i in range(n_spots)]
        prot_adata.var_names = [f"Protein_{i}" for i in range(n_proteins)]
        prot_adata.obsm["spatial"] = gex_adata.obsm["spatial"].copy()

        # Initialize model
        model = CitegeistModel(
            sample_name="test_multimodal",
            output_folder=str(tmp_path),
            simulation=True,
            gene_expression_adata=gex_adata,
            antibody_capture_adata=prot_adata,
        )

        # Mock Pass 1 results (normally from run_cell_proportion_model)
        cell_type_names = ["TypeA", "TypeB", "TypeC"]
        Y_protein = np.random.dirichlet([2, 2, 2], size=n_spots)
        model.cell_prop_global_results = pd.DataFrame(
            Y_protein,
            index=gex_adata.obs_names,
            columns=cell_type_names,
        )

        # Run multimodal refinement
        result = model.run_multimodal_refinement(
            n_anchors=5,
            lambda_prior=1.0,
            max_iterations=5,
        )

        # Verify outputs
        assert result is not None
        assert result.shape == (n_spots, n_types)
        assert hasattr(model, "cell_prop_refined_results")
        assert hasattr(model, "anchor_genes")
        assert hasattr(model, "expression_profiles")

        # Check files were saved
        assert (tmp_path / "test_multimodal_cell_prop_refined_results.csv").exists()
        assert (tmp_path / "test_multimodal_anchor_genes.json").exists()
```

**Step 2: Run integration test**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest CITEgeist/tests/test_multimodal_refinement.py::TestIntegrationWithCitegeistModel -v`

Expected: PASS

**Step 3: Commit**

```bash
git add CITEgeist/tests/test_multimodal_refinement.py
git commit -m "test: add integration test for multimodal refinement"
```

---

## Task 11: Run Full Test Suite and Final Commit

**Step 1: Run all multimodal refinement tests**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest CITEgeist/tests/test_multimodal_refinement.py -v`

Expected: All tests PASS

**Step 2: Run existing CITEgeist tests to ensure no regression**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest CITEgeist/tests/ -v --ignore=CITEgeist/tests/test_multimodal_refinement.py -x`

Expected: Existing tests still PASS

**Step 3: Final commit with feature complete message**

```bash
git add -A
git commit -m "feat(module3): complete multimodal refinement implementation

- Pass 1.5: Anchor gene selection via Pearson correlation + specificity
- Pass 2 EM: Iterative E-step (expression profiles) and M-step (proportions)
- Anchors locked during EM, weighted by correlation
- Integrated with CitegeistModel.run_multimodal_refinement()
- Full test coverage with unit and integration tests"
```

---

## Summary

| Task | Description | Files |
|------|-------------|-------|
| 1 | Create test file + first failing test | `tests/test_multimodal_refinement.py` |
| 2 | Implement `select_anchor_genes` | `model/multimodal_refinement.py` |
| 3 | Add E-step test | `tests/test_multimodal_refinement.py` |
| 4 | Implement `compute_expression_profiles` | `model/multimodal_refinement.py` |
| 5 | Add M-step test | `tests/test_multimodal_refinement.py` |
| 6 | Implement `refine_proportions` | `model/multimodal_refinement.py` |
| 7 | Add full EM test | `tests/test_multimodal_refinement.py` |
| 8 | Implement `multimodal_em_refinement` | `model/multimodal_refinement.py` |
| 9 | Integrate with CitegeistModel | `model/citegeist_model.py` |
| 10 | Add integration test | `tests/test_multimodal_refinement.py` |
| 11 | Final test suite + commit | - |

**Estimated time:** 2-3 hours for implementation
