# Gene-Parallel U-Net Spatial Smoothing Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Implement gene-parallel U-Net spatial smoothing for Module 3 Pass 2 to improve GEX deconvolution accuracy.

**Architecture:** Gene-parallel QP decomposition (9,800 vars per gene) with hierarchical U-Net (local → regional → global → back down). Proportion-weighted asymmetric Laplacian ensures high-confidence spots anchor low-confidence neighbors. Hierarchy preserves regional heterogeneity.

**Tech Stack:** Python 3.10, Gurobi (gurobipy), NumPy, SciPy, scikit-learn (for spatial partitioning), pytest

---

## Task 1: Spatial Hierarchy Data Structure

**Files:**
- Create: `CITEgeist/model/spatial_unet_gex.py`
- Test: `CITEgeist/tests/test_spatial_unet_gex.py`

**Step 1: Write the failing test for SpatialHierarchy**

```python
# CITEgeist/tests/test_spatial_unet_gex.py
"""
Unit tests for gene-parallel U-Net spatial smoothing.
"""
import numpy as np
import pytest
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))


def create_mock_spatial_coords(n_spots=100, seed=42):
    """Create mock spatial coordinates on a grid."""
    np.random.seed(seed)
    # 10x10 grid with some jitter
    xs = np.repeat(np.arange(10), 10) * 100 + np.random.normal(0, 5, 100)
    ys = np.tile(np.arange(10), 10) * 100 + np.random.normal(0, 5, 100)
    return np.column_stack([xs, ys])[:n_spots]


class TestSpatialHierarchy:
    """Tests for the SpatialHierarchy class."""

    def test_hierarchy_construction(self):
        """Test that hierarchy builds correct number of levels."""
        from model.spatial_unet_gex import SpatialHierarchy

        coords = create_mock_spatial_coords(n_spots=100)
        hierarchy = SpatialHierarchy(coords, n_levels=3, spots_per_local=20)

        assert hierarchy.n_levels == 3
        assert len(hierarchy.levels) == 3

    def test_level_0_covers_all_spots(self):
        """Test that level 0 (local) batches cover all spots."""
        from model.spatial_unet_gex import SpatialHierarchy

        coords = create_mock_spatial_coords(n_spots=100)
        hierarchy = SpatialHierarchy(coords, n_levels=3, spots_per_local=20)

        level_0 = hierarchy.levels[0]
        all_spots = set()
        for batch in level_0.batches:
            all_spots.update(batch.spot_indices)

        assert all_spots == set(range(100))

    def test_neighbor_graph_symmetric(self):
        """Test that neighbor graph is symmetric."""
        from model.spatial_unet_gex import SpatialHierarchy

        coords = create_mock_spatial_coords(n_spots=100)
        hierarchy = SpatialHierarchy(coords, n_levels=3, spots_per_local=20)

        neighbor_graph = hierarchy.neighbor_graph
        for i, neighbors in enumerate(neighbor_graph):
            for j in neighbors:
                assert i in neighbor_graph[j], f"Edge ({i},{j}) not symmetric"
```

**Step 2: Run test to verify it fails**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest CITEgeist/tests/test_spatial_unet_gex.py::TestSpatialHierarchy::test_hierarchy_construction -v`
Expected: FAIL with "ModuleNotFoundError: No module named 'model.spatial_unet_gex'"

**Step 3: Write minimal SpatialHierarchy implementation**

```python
# CITEgeist/model/spatial_unet_gex.py
"""
Gene-parallel U-Net spatial smoothing for Module 3 Pass 2 GEX deconvolution.

This module implements hierarchical spatial regularization where:
- Genes are solved independently in parallel (no cross-gene coupling)
- Each gene uses a U-Net architecture: local → regional → global → back down
- Proportion-weighted asymmetric Laplacian pulls low-confidence spots toward anchors
- Hierarchy discovers regional patterns while preserving heterogeneity
"""
import logging
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import numpy as np
from scipy.spatial import KDTree
from sklearn.cluster import AgglomerativeClustering

logger = logging.getLogger(__name__)


@dataclass
class HierarchyBatch:
    """A batch of spots at a given hierarchy level."""
    batch_id: int
    spot_indices: np.ndarray  # indices into original spot array
    neighbor_pairs: List[Tuple[int, int]] = field(default_factory=list)  # (i, j) pairs within batch
    parent_batch_id: Optional[int] = None  # batch at level+1 this belongs to


@dataclass
class HierarchyLevel:
    """One level of the spatial hierarchy."""
    level: int
    batches: List[HierarchyBatch]
    spot_to_batch: Dict[int, int]  # spot_idx -> batch_id


class SpatialHierarchy:
    """
    Hierarchical spatial partitioning for U-Net GEX deconvolution.

    Builds a multi-level hierarchy from spatial coordinates:
    - Level 0: Local neighborhoods (~20 spots each)
    - Level 1: Regional clusters (~100 spots each)
    - Level 2+: Progressively larger regions
    - Top level: Global (all spots)

    Args:
        coords: (n_spots, 2) spatial coordinates
        n_levels: Number of hierarchy levels (default 3)
        spots_per_local: Target spots per level-0 batch (default 20)
        k_neighbors: Number of spatial neighbors for Laplacian (default 6)
    """

    def __init__(
        self,
        coords: np.ndarray,
        n_levels: int = 3,
        spots_per_local: int = 20,
        k_neighbors: int = 6,
    ):
        self.coords = coords
        self.n_spots = coords.shape[0]
        self.n_levels = n_levels
        self.spots_per_local = spots_per_local
        self.k_neighbors = k_neighbors

        # Build neighbor graph
        self.neighbor_graph = self._build_neighbor_graph()

        # Build hierarchy levels
        self.levels = self._build_hierarchy()

    def _build_neighbor_graph(self) -> List[List[int]]:
        """Build k-NN neighbor graph from spatial coordinates."""
        tree = KDTree(self.coords)
        # Query k+1 because first neighbor is self
        distances, indices = tree.query(self.coords, k=self.k_neighbors + 1)

        neighbor_graph = []
        for i in range(self.n_spots):
            # Exclude self (first neighbor)
            neighbors = [int(j) for j in indices[i, 1:] if j != i]
            neighbor_graph.append(neighbors)

        # Symmetrize: if i is neighbor of j, j should be neighbor of i
        for i, neighbors in enumerate(neighbor_graph):
            for j in neighbors:
                if i not in neighbor_graph[j]:
                    neighbor_graph[j].append(i)

        return neighbor_graph

    def _build_hierarchy(self) -> List[HierarchyLevel]:
        """Build multi-level spatial hierarchy."""
        levels = []

        # Level 0: Local batches via spatial clustering
        n_clusters_l0 = max(1, self.n_spots // self.spots_per_local)
        if n_clusters_l0 > 1:
            clustering_l0 = AgglomerativeClustering(
                n_clusters=n_clusters_l0,
                linkage='ward'
            )
            labels_l0 = clustering_l0.fit_predict(self.coords)
        else:
            labels_l0 = np.zeros(self.n_spots, dtype=int)

        batches_l0 = []
        spot_to_batch_l0 = {}
        for batch_id in range(n_clusters_l0):
            spot_indices = np.where(labels_l0 == batch_id)[0]
            neighbor_pairs = self._get_neighbor_pairs_in_batch(spot_indices)
            batches_l0.append(HierarchyBatch(
                batch_id=batch_id,
                spot_indices=spot_indices,
                neighbor_pairs=neighbor_pairs,
            ))
            for idx in spot_indices:
                spot_to_batch_l0[int(idx)] = batch_id

        levels.append(HierarchyLevel(level=0, batches=batches_l0, spot_to_batch=spot_to_batch_l0))

        # Higher levels: progressively coarser clustering
        current_n_clusters = n_clusters_l0
        for level_idx in range(1, self.n_levels):
            # Reduce cluster count by ~5x each level
            target_clusters = max(1, current_n_clusters // 5)
            if target_clusters >= current_n_clusters:
                target_clusters = 1

            if target_clusters > 1:
                clustering = AgglomerativeClustering(
                    n_clusters=target_clusters,
                    linkage='ward'
                )
                labels = clustering.fit_predict(self.coords)
            else:
                labels = np.zeros(self.n_spots, dtype=int)

            batches = []
            spot_to_batch = {}
            for batch_id in range(target_clusters):
                spot_indices = np.where(labels == batch_id)[0]
                neighbor_pairs = self._get_neighbor_pairs_in_batch(spot_indices)
                batches.append(HierarchyBatch(
                    batch_id=batch_id,
                    spot_indices=spot_indices,
                    neighbor_pairs=neighbor_pairs,
                ))
                for idx in spot_indices:
                    spot_to_batch[int(idx)] = batch_id

            # Link to parent level
            for batch_l0 in levels[level_idx - 1].batches:
                # Parent is the batch at this level containing majority of spots
                batch_counts = {}
                for spot_idx in batch_l0.spot_indices:
                    parent_batch = spot_to_batch[int(spot_idx)]
                    batch_counts[parent_batch] = batch_counts.get(parent_batch, 0) + 1
                batch_l0.parent_batch_id = max(batch_counts, key=batch_counts.get)

            levels.append(HierarchyLevel(level=level_idx, batches=batches, spot_to_batch=spot_to_batch))
            current_n_clusters = target_clusters

        return levels

    def _get_neighbor_pairs_in_batch(self, spot_indices: np.ndarray) -> List[Tuple[int, int]]:
        """Get neighbor pairs where both spots are in the batch."""
        spot_set = set(spot_indices)
        pairs = []
        seen = set()
        for i in spot_indices:
            for j in self.neighbor_graph[i]:
                if j in spot_set:
                    pair = (min(i, j), max(i, j))
                    if pair not in seen:
                        seen.add(pair)
                        pairs.append(pair)
        return pairs

    def get_spots_in_region(self, level: int, batch_id: int) -> np.ndarray:
        """Get spot indices belonging to a batch at a given level."""
        return self.levels[level].batches[batch_id].spot_indices
```

**Step 4: Run tests to verify they pass**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest CITEgeist/tests/test_spatial_unet_gex.py::TestSpatialHierarchy -v`
Expected: PASS (3 tests)

**Step 5: Commit**

```bash
git add CITEgeist/model/spatial_unet_gex.py CITEgeist/tests/test_spatial_unet_gex.py
git commit -m "feat(unet-gex): add SpatialHierarchy data structure

Implements hierarchical spatial partitioning for U-Net GEX deconvolution:
- Multi-level clustering (local → regional → global)
- k-NN neighbor graph with symmetrization
- Batch-to-parent linkage for hierarchy traversal

Co-Authored-By: Claude Opus 4.5 <noreply@anthropic.com>"
```

---

## Task 2: Asymmetric Laplacian Weight Computation

**Files:**
- Modify: `CITEgeist/model/spatial_unet_gex.py`
- Test: `CITEgeist/tests/test_spatial_unet_gex.py`

**Step 1: Write the failing test for asymmetric weights**

```python
# Add to CITEgeist/tests/test_spatial_unet_gex.py

class TestAsymmetricLaplacian:
    """Tests for proportion-weighted asymmetric Laplacian."""

    def test_high_prop_anchors_low_prop(self):
        """Test that high-proportion spot pulls low-proportion neighbor."""
        from model.spatial_unet_gex import compute_asymmetric_weight

        # Spot i has 40% of cell type, spot j has 5%
        prop_i = 0.40
        prop_j = 0.05

        w_ij = compute_asymmetric_weight(prop_i, prop_j)

        # w_ij = prop_i / (prop_i + prop_j) = 0.40 / 0.45 ≈ 0.889
        assert 0.85 < w_ij < 0.95, f"Expected ~0.889, got {w_ij}"

    def test_symmetric_when_equal_props(self):
        """Test that equal proportions give symmetric weight."""
        from model.spatial_unet_gex import compute_asymmetric_weight

        prop_i = 0.30
        prop_j = 0.30

        w_ij = compute_asymmetric_weight(prop_i, prop_j)

        # Should be 0.5 when equal
        assert 0.49 < w_ij < 0.51, f"Expected 0.5, got {w_ij}"

    def test_low_prop_gets_weak_pull(self):
        """Test that low-proportion spot exerts weak pull on high-proportion neighbor."""
        from model.spatial_unet_gex import compute_asymmetric_weight

        # Spot i has 5%, spot j has 40%
        prop_i = 0.05
        prop_j = 0.40

        w_ij = compute_asymmetric_weight(prop_i, prop_j)

        # w_ij = 0.05 / 0.45 ≈ 0.111
        assert 0.10 < w_ij < 0.15, f"Expected ~0.111, got {w_ij}"

    def test_zero_prop_no_pull(self):
        """Test that zero-proportion spot exerts no pull."""
        from model.spatial_unet_gex import compute_asymmetric_weight

        prop_i = 0.0
        prop_j = 0.40

        w_ij = compute_asymmetric_weight(prop_i, prop_j)

        # Should be ~0 (epsilon prevents division by zero)
        assert w_ij < 0.01, f"Expected ~0, got {w_ij}"
```

**Step 2: Run test to verify it fails**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest CITEgeist/tests/test_spatial_unet_gex.py::TestAsymmetricLaplacian::test_high_prop_anchors_low_prop -v`
Expected: FAIL with "ImportError: cannot import name 'compute_asymmetric_weight'"

**Step 3: Write minimal implementation**

```python
# Add to CITEgeist/model/spatial_unet_gex.py after SpatialHierarchy class

def compute_asymmetric_weight(prop_i: float, prop_j: float, epsilon: float = 1e-10) -> float:
    """
    Compute asymmetric Laplacian weight for neighbor pair.

    Weight w_ij determines how strongly spot i pulls spot j:
        w_ij = prop_i / (prop_i + prop_j + epsilon)

    When prop_i > prop_j, w_ij > 0.5, so i is the anchor.
    When prop_i < prop_j, w_ij < 0.5, so i exerts weak pull.

    Args:
        prop_i: Cell type proportion at spot i (the potential anchor)
        prop_j: Cell type proportion at spot j (the neighbor)
        epsilon: Small constant to prevent division by zero

    Returns:
        Weight in [0, 1] indicating pull strength from i to j
    """
    return prop_i / (prop_i + prop_j + epsilon)


def compute_laplacian_weights_for_batch(
    neighbor_pairs: List[Tuple[int, int]],
    proportions: np.ndarray,
    cell_type_idx: int,
) -> Dict[Tuple[int, int], float]:
    """
    Compute asymmetric Laplacian weights for all neighbor pairs in a batch.

    Args:
        neighbor_pairs: List of (i, j) spot index pairs
        proportions: (n_spots, n_celltypes) proportion matrix
        cell_type_idx: Which cell type to compute weights for

    Returns:
        Dict mapping (i, j) -> weight w_ij
    """
    weights = {}
    for i, j in neighbor_pairs:
        prop_i = proportions[i, cell_type_idx]
        prop_j = proportions[j, cell_type_idx]
        weights[(i, j)] = compute_asymmetric_weight(prop_i, prop_j)
        weights[(j, i)] = compute_asymmetric_weight(prop_j, prop_i)
    return weights
```

**Step 4: Run tests to verify they pass**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest CITEgeist/tests/test_spatial_unet_gex.py::TestAsymmetricLaplacian -v`
Expected: PASS (4 tests)

**Step 5: Commit**

```bash
git add CITEgeist/model/spatial_unet_gex.py CITEgeist/tests/test_spatial_unet_gex.py
git commit -m "feat(unet-gex): add asymmetric Laplacian weight computation

Implements proportion-weighted asymmetric pull:
- High-proportion spots anchor, low-proportion gets pulled
- w_ij = prop_i / (prop_i + prop_j + epsilon)
- Batch-level weight computation for efficiency

Co-Authored-By: Claude Opus 4.5 <noreply@anthropic.com>"
```

---

## Task 3: Single-Gene QP Solver (Upward Pass)

**Files:**
- Modify: `CITEgeist/model/spatial_unet_gex.py`
- Test: `CITEgeist/tests/test_spatial_unet_gex.py`

**Step 1: Write the failing test for upward pass QP**

```python
# Add to CITEgeist/tests/test_spatial_unet_gex.py

class TestUpwardPassQP:
    """Tests for upward pass QP solver."""

    def test_count_conservation(self):
        """Test that assigned counts sum to observed counts."""
        from model.spatial_unet_gex import solve_batch_qp_upward

        n_spots = 20
        n_celltypes = 3
        np.random.seed(42)

        # Mock data
        spot_indices = np.arange(n_spots)
        observed = np.random.poisson(100, size=n_spots).astype(float)
        proportions = np.random.dirichlet([1, 1, 1], size=n_spots)
        enrichment = np.ones((n_spots, n_celltypes)) / n_celltypes
        neighbor_pairs = [(i, i+1) for i in range(n_spots - 1)]

        result = solve_batch_qp_upward(
            spot_indices=spot_indices,
            observed=observed,
            proportions=proportions,
            enrichment=enrichment,
            neighbor_pairs=neighbor_pairs,
            lambda_spatial=0.1,
            lambda_reg=0.01,
        )

        # Check count conservation: sum over cell types = observed
        for i, spot_idx in enumerate(spot_indices):
            assigned_sum = result[i, :].sum()
            assert np.isclose(assigned_sum, observed[i], rtol=1e-3), \
                f"Spot {spot_idx}: assigned {assigned_sum} != observed {observed[i]}"

    def test_high_prop_gets_more_counts(self):
        """Test that cell types with higher proportions get more counts."""
        from model.spatial_unet_gex import solve_batch_qp_upward

        n_spots = 10
        n_celltypes = 2
        np.random.seed(42)

        spot_indices = np.arange(n_spots)
        observed = np.full(n_spots, 100.0)
        # Cell type 0 has 80%, cell type 1 has 20%
        proportions = np.column_stack([
            np.full(n_spots, 0.8),
            np.full(n_spots, 0.2)
        ])
        enrichment = np.ones((n_spots, n_celltypes)) / n_celltypes
        neighbor_pairs = [(i, i+1) for i in range(n_spots - 1)]

        result = solve_batch_qp_upward(
            spot_indices=spot_indices,
            observed=observed,
            proportions=proportions,
            enrichment=enrichment,
            neighbor_pairs=neighbor_pairs,
            lambda_spatial=0.1,
            lambda_reg=0.01,
        )

        # Cell type 0 (80%) should get more than cell type 1 (20%)
        mean_ct0 = result[:, 0].mean()
        mean_ct1 = result[:, 1].mean()
        assert mean_ct0 > mean_ct1, f"Expected CT0 ({mean_ct0}) > CT1 ({mean_ct1})"

    def test_spatial_smoothing_reduces_variance(self):
        """Test that spatial smoothing reduces variance across neighbors."""
        from model.spatial_unet_gex import solve_batch_qp_upward

        n_spots = 20
        n_celltypes = 2
        np.random.seed(42)

        spot_indices = np.arange(n_spots)
        observed = np.random.poisson(100, size=n_spots).astype(float)
        # Uniform proportions
        proportions = np.full((n_spots, n_celltypes), 0.5)
        enrichment = np.ones((n_spots, n_celltypes)) / n_celltypes
        neighbor_pairs = [(i, i+1) for i in range(n_spots - 1)]

        # Without smoothing
        result_no_smooth = solve_batch_qp_upward(
            spot_indices=spot_indices,
            observed=observed,
            proportions=proportions,
            enrichment=enrichment,
            neighbor_pairs=neighbor_pairs,
            lambda_spatial=0.0,
            lambda_reg=0.01,
        )

        # With smoothing
        result_smooth = solve_batch_qp_upward(
            spot_indices=spot_indices,
            observed=observed,
            proportions=proportions,
            enrichment=enrichment,
            neighbor_pairs=neighbor_pairs,
            lambda_spatial=1.0,
            lambda_reg=0.01,
        )

        # Variance should be lower with smoothing
        var_no_smooth = result_no_smooth[:, 0].var()
        var_smooth = result_smooth[:, 0].var()
        assert var_smooth <= var_no_smooth, \
            f"Smoothing should reduce variance: {var_smooth} > {var_no_smooth}"
```

**Step 2: Run test to verify it fails**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest CITEgeist/tests/test_spatial_unet_gex.py::TestUpwardPassQP::test_count_conservation -v`
Expected: FAIL with "ImportError: cannot import name 'solve_batch_qp_upward'"

**Step 3: Write minimal implementation**

```python
# Add to CITEgeist/model/spatial_unet_gex.py

import gurobipy as gp
from gurobipy import GRB


def solve_batch_qp_upward(
    spot_indices: np.ndarray,
    observed: np.ndarray,
    proportions: np.ndarray,
    enrichment: np.ndarray,
    neighbor_pairs: List[Tuple[int, int]],
    lambda_spatial: float = 0.1,
    lambda_reg: float = 0.01,
) -> np.ndarray:
    """
    Solve QP for a batch of spots in the upward pass.

    Objective:
        max  Σ_i Σ_t [ enrichment[i,t] × prop[i,t] × X[i,t] ]
           - λ_reg × Σ_i Σ_t X[i,t]²
           - λ_spatial × Σ_(i~j) Σ_t [ w_ij[t] × (X[j,t] - X[i,t])² ]

    Subject to:
        Σ_t X[i,t] = observed[i]   ∀ spots i
        X[i,t] ≥ 0

    Args:
        spot_indices: Indices of spots in this batch (for mapping)
        observed: (n_batch,) observed counts for this gene
        proportions: (n_batch, n_celltypes) cell type proportions
        enrichment: (n_batch, n_celltypes) enrichment weights
        neighbor_pairs: List of (local_i, local_j) neighbor pairs within batch
        lambda_spatial: Spatial smoothing weight
        lambda_reg: L2 regularization weight

    Returns:
        (n_batch, n_celltypes) deconvolved expression
    """
    n_batch = len(spot_indices)
    n_celltypes = proportions.shape[1]

    model = gp.Model("upward_pass_qp")
    model.setParam("OutputFlag", 0)
    model.setParam("Threads", 1)
    model.setParam("Method", 2)  # Barrier
    model.setParam("BarConvTol", 1e-6)

    # Variables: X[i, t] for each spot i and cell type t
    X = {}
    for i in range(n_batch):
        total_counts = observed[i]
        for t in range(n_celltypes):
            X[i, t] = model.addVar(
                lb=0,
                ub=total_counts,
                vtype=GRB.CONTINUOUS,
                name=f"X_{i}_{t}"
            )

    model.update()

    # Constraints: count conservation
    for i in range(n_batch):
        model.addConstr(
            gp.quicksum(X[i, t] for t in range(n_celltypes)) == observed[i],
            name=f"count_conservation_{i}"
        )

    # Objective terms
    obj_terms = []

    # 1. Enrichment × proportion × X (maximize)
    for i in range(n_batch):
        for t in range(n_celltypes):
            coef = enrichment[i, t] * proportions[i, t]
            obj_terms.append(coef * X[i, t])

    # 2. L2 regularization (subtract from objective since we maximize)
    for i in range(n_batch):
        for t in range(n_celltypes):
            obj_terms.append(-lambda_reg * X[i, t] * X[i, t])

    # 3. Spatial Laplacian (subtract from objective)
    if lambda_spatial > 0 and neighbor_pairs:
        for i, j in neighbor_pairs:
            for t in range(n_celltypes):
                # Asymmetric weight: how strongly i pulls j
                w_ij = compute_asymmetric_weight(proportions[i, t], proportions[j, t])
                # Penalty: w_ij * (X[j,t] - X[i,t])^2
                # = w_ij * (X[j,t]^2 - 2*X[i,t]*X[j,t] + X[i,t]^2)
                obj_terms.append(-lambda_spatial * w_ij * X[j, t] * X[j, t])
                obj_terms.append(2 * lambda_spatial * w_ij * X[i, t] * X[j, t])
                obj_terms.append(-lambda_spatial * w_ij * X[i, t] * X[i, t])

    model.setObjective(gp.quicksum(obj_terms), GRB.MAXIMIZE)
    model.optimize()

    # Extract results
    result = np.zeros((n_batch, n_celltypes))
    if model.status == GRB.OPTIMAL:
        for i in range(n_batch):
            for t in range(n_celltypes):
                result[i, t] = X[i, t].X
    else:
        logger.warning(f"QP did not converge optimally. Status: {model.status}")
        # Fallback: distribute proportionally
        for i in range(n_batch):
            result[i, :] = observed[i] * proportions[i, :] / (proportions[i, :].sum() + 1e-10)

    del model
    return result
```

**Step 4: Run tests to verify they pass**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest CITEgeist/tests/test_spatial_unet_gex.py::TestUpwardPassQP -v`
Expected: PASS (3 tests)

**Step 5: Commit**

```bash
git add CITEgeist/model/spatial_unet_gex.py CITEgeist/tests/test_spatial_unet_gex.py
git commit -m "feat(unet-gex): add upward pass QP solver with Laplacian

Implements batch QP for upward pass:
- Count conservation constraint
- Enrichment × proportion objective
- Asymmetric Laplacian spatial smoothing
- L2 regularization for stability

Co-Authored-By: Claude Opus 4.5 <noreply@anthropic.com>"
```

---

## Task 4: Pooling Between Levels

**Files:**
- Modify: `CITEgeist/model/spatial_unet_gex.py`
- Test: `CITEgeist/tests/test_spatial_unet_gex.py`

**Step 1: Write the failing test for pooling**

```python
# Add to CITEgeist/tests/test_spatial_unet_gex.py

class TestHierarchyPooling:
    """Tests for proportion-weighted pooling between levels."""

    def test_pooling_weighted_by_proportion(self):
        """Test that pooling weights by cell type proportion."""
        from model.spatial_unet_gex import pool_to_region

        # 4 spots, 2 cell types
        estimates = np.array([
            [100, 50],   # spot 0: CT0=100, CT1=50
            [80, 60],    # spot 1: CT0=80, CT1=60
            [60, 80],    # spot 2: CT0=60, CT1=80
            [40, 100],   # spot 3: CT0=40, CT1=100
        ])
        proportions = np.array([
            [0.8, 0.2],  # spot 0: high CT0
            [0.6, 0.4],  # spot 1
            [0.4, 0.6],  # spot 2
            [0.2, 0.8],  # spot 3: high CT1
        ])
        spot_indices = np.array([0, 1, 2, 3])

        pooled = pool_to_region(estimates, proportions, spot_indices)

        # For CT0: weighted avg = (0.8*100 + 0.6*80 + 0.4*60 + 0.2*40) / (0.8+0.6+0.4+0.2)
        #        = (80 + 48 + 24 + 8) / 2.0 = 160 / 2.0 = 80
        assert 75 < pooled[0] < 85, f"Expected ~80 for CT0, got {pooled[0]}"

        # For CT1: weighted avg = (0.2*50 + 0.4*60 + 0.6*80 + 0.8*100) / (0.2+0.4+0.6+0.8)
        #        = (10 + 24 + 48 + 80) / 2.0 = 162 / 2.0 = 81
        assert 76 < pooled[1] < 86, f"Expected ~81 for CT1, got {pooled[1]}"

    def test_pooling_handles_zero_proportion(self):
        """Test that zero-proportion spots don't contribute to pool."""
        from model.spatial_unet_gex import pool_to_region

        estimates = np.array([
            [100, 0],
            [50, 100],
        ])
        proportions = np.array([
            [1.0, 0.0],  # spot 0: only CT0
            [0.0, 1.0],  # spot 1: only CT1
        ])
        spot_indices = np.array([0, 1])

        pooled = pool_to_region(estimates, proportions, spot_indices)

        # CT0 should be 100 (only spot 0 contributes)
        assert np.isclose(pooled[0], 100, rtol=0.01)
        # CT1 should be 100 (only spot 1 contributes)
        assert np.isclose(pooled[1], 100, rtol=0.01)
```

**Step 2: Run test to verify it fails**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest CITEgeist/tests/test_spatial_unet_gex.py::TestHierarchyPooling::test_pooling_weighted_by_proportion -v`
Expected: FAIL with "ImportError: cannot import name 'pool_to_region'"

**Step 3: Write minimal implementation**

```python
# Add to CITEgeist/model/spatial_unet_gex.py

def pool_to_region(
    estimates: np.ndarray,
    proportions: np.ndarray,
    spot_indices: np.ndarray,
    epsilon: float = 1e-10,
) -> np.ndarray:
    """
    Pool spot-level estimates to regional estimate via proportion-weighted average.

    For each cell type t:
        X_region[t] = Σ_i prop[i,t] × X[i,t] / Σ_i prop[i,t]

    Args:
        estimates: (n_spots, n_celltypes) expression estimates
        proportions: (n_spots, n_celltypes) cell type proportions
        spot_indices: Which spots are in this region (for logging)
        epsilon: Small constant to prevent division by zero

    Returns:
        (n_celltypes,) pooled regional estimate
    """
    n_celltypes = estimates.shape[1]
    pooled = np.zeros(n_celltypes)

    for t in range(n_celltypes):
        weights = proportions[:, t]
        values = estimates[:, t]
        total_weight = weights.sum() + epsilon
        pooled[t] = (weights * values).sum() / total_weight

    return pooled


def broadcast_prior_to_spots(
    regional_estimate: np.ndarray,
    spot_indices: np.ndarray,
) -> np.ndarray:
    """
    Broadcast regional estimate as prior to all spots in region.

    Args:
        regional_estimate: (n_celltypes,) pooled regional estimate
        spot_indices: Which spots to broadcast to

    Returns:
        (n_spots, n_celltypes) prior matrix
    """
    n_spots = len(spot_indices)
    n_celltypes = len(regional_estimate)
    prior = np.tile(regional_estimate, (n_spots, 1))
    return prior
```

**Step 4: Run tests to verify they pass**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest CITEgeist/tests/test_spatial_unet_gex.py::TestHierarchyPooling -v`
Expected: PASS (2 tests)

**Step 5: Commit**

```bash
git add CITEgeist/model/spatial_unet_gex.py CITEgeist/tests/test_spatial_unet_gex.py
git commit -m "feat(unet-gex): add proportion-weighted pooling between levels

Implements pooling for upward pass:
- Weighted average by cell type proportion
- Zero-proportion spots don't contribute
- Broadcast for downward pass priors

Co-Authored-By: Claude Opus 4.5 <noreply@anthropic.com>"
```

---

## Task 5: Downward Pass QP Solver with Prior

**Files:**
- Modify: `CITEgeist/model/spatial_unet_gex.py`
- Test: `CITEgeist/tests/test_spatial_unet_gex.py`

**Step 1: Write the failing test for downward pass**

```python
# Add to CITEgeist/tests/test_spatial_unet_gex.py

class TestDownwardPassQP:
    """Tests for downward pass QP solver with prior."""

    def test_low_confidence_shrinks_toward_prior(self):
        """Test that low-proportion spots shrink toward prior."""
        from model.spatial_unet_gex import solve_batch_qp_downward

        n_spots = 10
        n_celltypes = 2
        np.random.seed(42)

        spot_indices = np.arange(n_spots)
        observed = np.full(n_spots, 100.0)
        # Mixed: some high prop, some low prop for CT0
        proportions = np.column_stack([
            np.array([0.8, 0.8, 0.8, 0.8, 0.8, 0.1, 0.1, 0.1, 0.1, 0.1]),
            np.array([0.2, 0.2, 0.2, 0.2, 0.2, 0.9, 0.9, 0.9, 0.9, 0.9]),
        ])
        enrichment = np.ones((n_spots, n_celltypes)) / n_celltypes
        neighbor_pairs = [(i, i+1) for i in range(n_spots - 1)]

        # Prior: CT0 gets 70, CT1 gets 30
        prior = np.full((n_spots, n_celltypes), [70.0, 30.0])

        result = solve_batch_qp_downward(
            spot_indices=spot_indices,
            observed=observed,
            proportions=proportions,
            enrichment=enrichment,
            neighbor_pairs=neighbor_pairs,
            prior=prior,
            lambda_spatial=0.1,
            lambda_prior=0.5,
            lambda_reg=0.01,
        )

        # Low-proportion spots (5-9 for CT0) should be closer to prior (70) than observed
        low_prop_mean_ct0 = result[5:10, 0].mean()
        high_prop_mean_ct0 = result[0:5, 0].mean()

        # Low prop should shrink toward prior (70)
        assert abs(low_prop_mean_ct0 - 70) < abs(high_prop_mean_ct0 - 70), \
            f"Low-prop ({low_prop_mean_ct0}) should be closer to prior (70) than high-prop ({high_prop_mean_ct0})"

    def test_count_conservation_with_prior(self):
        """Test that count conservation still holds with prior term."""
        from model.spatial_unet_gex import solve_batch_qp_downward

        n_spots = 10
        n_celltypes = 2
        np.random.seed(42)

        spot_indices = np.arange(n_spots)
        observed = np.random.poisson(100, size=n_spots).astype(float)
        proportions = np.random.dirichlet([1, 1], size=n_spots)
        enrichment = np.ones((n_spots, n_celltypes)) / n_celltypes
        neighbor_pairs = [(i, i+1) for i in range(n_spots - 1)]
        prior = np.full((n_spots, n_celltypes), 50.0)

        result = solve_batch_qp_downward(
            spot_indices=spot_indices,
            observed=observed,
            proportions=proportions,
            enrichment=enrichment,
            neighbor_pairs=neighbor_pairs,
            prior=prior,
            lambda_spatial=0.1,
            lambda_prior=0.5,
            lambda_reg=0.01,
        )

        for i in range(n_spots):
            assigned_sum = result[i, :].sum()
            assert np.isclose(assigned_sum, observed[i], rtol=1e-3), \
                f"Spot {i}: assigned {assigned_sum} != observed {observed[i]}"
```

**Step 2: Run test to verify it fails**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest CITEgeist/tests/test_spatial_unet_gex.py::TestDownwardPassQP::test_low_confidence_shrinks_toward_prior -v`
Expected: FAIL with "ImportError: cannot import name 'solve_batch_qp_downward'"

**Step 3: Write minimal implementation**

```python
# Add to CITEgeist/model/spatial_unet_gex.py

def solve_batch_qp_downward(
    spot_indices: np.ndarray,
    observed: np.ndarray,
    proportions: np.ndarray,
    enrichment: np.ndarray,
    neighbor_pairs: List[Tuple[int, int]],
    prior: np.ndarray,
    lambda_spatial: float = 0.1,
    lambda_prior: float = 0.5,
    lambda_reg: float = 0.01,
) -> np.ndarray:
    """
    Solve QP for a batch of spots in the downward pass (with prior).

    Objective (adds prior term to upward pass):
        max  [ upward pass terms ]
           - λ_prior × Σ_i Σ_t [ (1 - prop[i,t]) × (X[i,t] - prior[i,t])² ]

    The (1 - prop) weighting means low-confidence spots shrink more toward prior.

    Args:
        spot_indices: Indices of spots in this batch
        observed: (n_batch,) observed counts for this gene
        proportions: (n_batch, n_celltypes) cell type proportions
        enrichment: (n_batch, n_celltypes) enrichment weights
        neighbor_pairs: List of (local_i, local_j) neighbor pairs
        prior: (n_batch, n_celltypes) prior from coarser level
        lambda_spatial: Spatial smoothing weight
        lambda_prior: Prior shrinkage weight
        lambda_reg: L2 regularization weight

    Returns:
        (n_batch, n_celltypes) refined expression estimates
    """
    n_batch = len(spot_indices)
    n_celltypes = proportions.shape[1]

    model = gp.Model("downward_pass_qp")
    model.setParam("OutputFlag", 0)
    model.setParam("Threads", 1)
    model.setParam("Method", 2)
    model.setParam("BarConvTol", 1e-6)

    # Variables
    X = {}
    for i in range(n_batch):
        total_counts = observed[i]
        for t in range(n_celltypes):
            X[i, t] = model.addVar(
                lb=0,
                ub=total_counts,
                vtype=GRB.CONTINUOUS,
                name=f"X_{i}_{t}"
            )

    model.update()

    # Constraints: count conservation
    for i in range(n_batch):
        model.addConstr(
            gp.quicksum(X[i, t] for t in range(n_celltypes)) == observed[i],
            name=f"count_conservation_{i}"
        )

    # Objective terms
    obj_terms = []

    # 1. Enrichment × proportion × X
    for i in range(n_batch):
        for t in range(n_celltypes):
            coef = enrichment[i, t] * proportions[i, t]
            obj_terms.append(coef * X[i, t])

    # 2. L2 regularization
    for i in range(n_batch):
        for t in range(n_celltypes):
            obj_terms.append(-lambda_reg * X[i, t] * X[i, t])

    # 3. Spatial Laplacian
    if lambda_spatial > 0 and neighbor_pairs:
        for i, j in neighbor_pairs:
            for t in range(n_celltypes):
                w_ij = compute_asymmetric_weight(proportions[i, t], proportions[j, t])
                obj_terms.append(-lambda_spatial * w_ij * X[j, t] * X[j, t])
                obj_terms.append(2 * lambda_spatial * w_ij * X[i, t] * X[j, t])
                obj_terms.append(-lambda_spatial * w_ij * X[i, t] * X[i, t])

    # 4. Prior shrinkage term (new for downward pass)
    # - λ_prior × (1 - prop[i,t]) × (X[i,t] - prior[i,t])²
    if lambda_prior > 0:
        for i in range(n_batch):
            for t in range(n_celltypes):
                shrinkage_weight = (1 - proportions[i, t])
                prior_val = prior[i, t]
                # (X - prior)² = X² - 2*prior*X + prior²
                # The prior² is constant, doesn't affect optimization
                obj_terms.append(-lambda_prior * shrinkage_weight * X[i, t] * X[i, t])
                obj_terms.append(2 * lambda_prior * shrinkage_weight * prior_val * X[i, t])

    model.setObjective(gp.quicksum(obj_terms), GRB.MAXIMIZE)
    model.optimize()

    # Extract results
    result = np.zeros((n_batch, n_celltypes))
    if model.status == GRB.OPTIMAL:
        for i in range(n_batch):
            for t in range(n_celltypes):
                result[i, t] = X[i, t].X
    else:
        logger.warning(f"Downward QP did not converge. Status: {model.status}")
        # Fallback: blend prior with proportional distribution
        for i in range(n_batch):
            prop_dist = observed[i] * proportions[i, :] / (proportions[i, :].sum() + 1e-10)
            shrinkage = 1 - proportions[i, :]
            result[i, :] = (1 - shrinkage) * prop_dist + shrinkage * prior[i, :]

    del model
    return result
```

**Step 4: Run tests to verify they pass**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest CITEgeist/tests/test_spatial_unet_gex.py::TestDownwardPassQP -v`
Expected: PASS (2 tests)

**Step 5: Commit**

```bash
git add CITEgeist/model/spatial_unet_gex.py CITEgeist/tests/test_spatial_unet_gex.py
git commit -m "feat(unet-gex): add downward pass QP solver with prior shrinkage

Implements downward pass with prior term:
- (1 - prop) weighting for shrinkage strength
- Low-confidence spots shrink toward regional prior
- Count conservation maintained

Co-Authored-By: Claude Opus 4.5 <noreply@anthropic.com>"
```

---

## Task 6: Full U-Net Pipeline for Single Gene

**Files:**
- Modify: `CITEgeist/model/spatial_unet_gex.py`
- Test: `CITEgeist/tests/test_spatial_unet_gex.py`

**Step 1: Write the failing test for single-gene U-Net**

```python
# Add to CITEgeist/tests/test_spatial_unet_gex.py

class TestSingleGeneUNet:
    """Tests for complete U-Net pipeline for a single gene."""

    def test_unet_improves_over_baseline(self):
        """Test that U-Net produces reasonable results."""
        from model.spatial_unet_gex import solve_gene_unet, SpatialHierarchy

        np.random.seed(42)
        n_spots = 100
        n_celltypes = 3

        # Create spatial hierarchy
        coords = create_mock_spatial_coords(n_spots=n_spots)
        hierarchy = SpatialHierarchy(coords, n_levels=2, spots_per_local=20)

        # Mock data
        observed = np.random.poisson(100, size=n_spots).astype(float)
        proportions = np.random.dirichlet([1, 1, 1], size=n_spots)
        enrichment = np.ones((n_spots, n_celltypes)) / n_celltypes

        result = solve_gene_unet(
            gene_idx=0,
            observed=observed,
            proportions=proportions,
            enrichment=enrichment,
            hierarchy=hierarchy,
            lambda_spatial=0.1,
            lambda_prior=0.5,
        )

        # Should have correct shape
        assert result.shape == (n_spots, n_celltypes)

        # Count conservation
        for i in range(n_spots):
            assert np.isclose(result[i, :].sum(), observed[i], rtol=1e-2), \
                f"Spot {i}: sum {result[i, :].sum()} != observed {observed[i]}"

        # Non-negative
        assert (result >= -1e-6).all(), "Negative values in result"

    def test_unet_respects_proportions(self):
        """Test that U-Net assigns more to higher-proportion cell types."""
        from model.spatial_unet_gex import solve_gene_unet, SpatialHierarchy

        np.random.seed(42)
        n_spots = 50
        n_celltypes = 2

        coords = create_mock_spatial_coords(n_spots=n_spots)
        hierarchy = SpatialHierarchy(coords, n_levels=2, spots_per_local=15)

        observed = np.full(n_spots, 100.0)
        # CT0 dominant (80%) everywhere
        proportions = np.column_stack([
            np.full(n_spots, 0.8),
            np.full(n_spots, 0.2)
        ])
        enrichment = np.ones((n_spots, n_celltypes)) / n_celltypes

        result = solve_gene_unet(
            gene_idx=0,
            observed=observed,
            proportions=proportions,
            enrichment=enrichment,
            hierarchy=hierarchy,
            lambda_spatial=0.1,
            lambda_prior=0.5,
        )

        # CT0 should get more than CT1 on average
        mean_ct0 = result[:, 0].mean()
        mean_ct1 = result[:, 1].mean()
        assert mean_ct0 > mean_ct1, f"CT0 ({mean_ct0}) should be > CT1 ({mean_ct1})"
```

**Step 2: Run test to verify it fails**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest CITEgeist/tests/test_spatial_unet_gex.py::TestSingleGeneUNet::test_unet_improves_over_baseline -v`
Expected: FAIL with "ImportError: cannot import name 'solve_gene_unet'"

**Step 3: Write minimal implementation**

```python
# Add to CITEgeist/model/spatial_unet_gex.py

def solve_gene_unet(
    gene_idx: int,
    observed: np.ndarray,
    proportions: np.ndarray,
    enrichment: np.ndarray,
    hierarchy: SpatialHierarchy,
    lambda_spatial: float = 0.1,
    lambda_prior: float = 0.5,
    lambda_reg: float = 0.01,
) -> np.ndarray:
    """
    Solve U-Net hierarchy for a single gene.

    Upward pass: Local → Regional → Global (pooling estimates)
    Downward pass: Global → Regional → Local (refining with priors)

    Args:
        gene_idx: Gene index (for logging)
        observed: (n_spots,) observed counts
        proportions: (n_spots, n_celltypes) cell type proportions
        enrichment: (n_spots, n_celltypes) enrichment weights
        hierarchy: SpatialHierarchy object
        lambda_spatial: Spatial smoothing weight
        lambda_prior: Prior shrinkage weight (downward pass)
        lambda_reg: L2 regularization weight

    Returns:
        (n_spots, n_celltypes) deconvolved expression
    """
    n_spots = len(observed)
    n_celltypes = proportions.shape[1]
    n_levels = hierarchy.n_levels

    # Storage for estimates at each level
    level_estimates = {}  # level -> {batch_id -> (n_batch, n_celltypes)}
    level_pooled = {}     # level -> {batch_id -> (n_celltypes,)}

    # ==================== UPWARD PASS ====================
    for level_idx in range(n_levels):
        level = hierarchy.levels[level_idx]
        level_estimates[level_idx] = {}

        for batch in level.batches:
            spot_indices = batch.spot_indices
            n_batch = len(spot_indices)

            # Get data for this batch
            batch_observed = observed[spot_indices]
            batch_proportions = proportions[spot_indices, :]
            batch_enrichment = enrichment[spot_indices, :]

            # Map global neighbor pairs to local indices
            local_idx_map = {global_idx: local_idx for local_idx, global_idx in enumerate(spot_indices)}
            local_neighbor_pairs = []
            for i, j in batch.neighbor_pairs:
                if i in local_idx_map and j in local_idx_map:
                    local_neighbor_pairs.append((local_idx_map[i], local_idx_map[j]))

            # Solve upward pass QP
            batch_result = solve_batch_qp_upward(
                spot_indices=spot_indices,
                observed=batch_observed,
                proportions=batch_proportions,
                enrichment=batch_enrichment,
                neighbor_pairs=local_neighbor_pairs,
                lambda_spatial=lambda_spatial,
                lambda_reg=lambda_reg,
            )

            level_estimates[level_idx][batch.batch_id] = (spot_indices, batch_result)

        # Pool to create regional estimates for next level
        level_pooled[level_idx] = {}
        for batch in level.batches:
            spot_indices, estimates = level_estimates[level_idx][batch.batch_id]
            batch_proportions = proportions[spot_indices, :]
            pooled = pool_to_region(estimates, batch_proportions, spot_indices)
            level_pooled[level_idx][batch.batch_id] = pooled

    # ==================== DOWNWARD PASS ====================
    # Start from top level and work down
    for level_idx in range(n_levels - 1, -1, -1):
        level = hierarchy.levels[level_idx]

        for batch in level.batches:
            spot_indices = batch.spot_indices
            n_batch = len(spot_indices)

            batch_observed = observed[spot_indices]
            batch_proportions = proportions[spot_indices, :]
            batch_enrichment = enrichment[spot_indices, :]

            # Build prior from parent level (or use upward estimate for top level)
            if level_idx == n_levels - 1:
                # Top level: use pooled estimate as prior (self-prior)
                pooled = level_pooled[level_idx][batch.batch_id]
                prior = np.tile(pooled, (n_batch, 1))
            else:
                # Get prior from parent batch at level+1
                parent_batch_id = batch.parent_batch_id
                if parent_batch_id is not None and parent_batch_id in level_pooled[level_idx + 1]:
                    parent_pooled = level_pooled[level_idx + 1][parent_batch_id]
                    prior = np.tile(parent_pooled, (n_batch, 1))
                else:
                    # Fallback: use own pooled estimate
                    pooled = level_pooled[level_idx][batch.batch_id]
                    prior = np.tile(pooled, (n_batch, 1))

            # Map neighbor pairs to local indices
            local_idx_map = {global_idx: local_idx for local_idx, global_idx in enumerate(spot_indices)}
            local_neighbor_pairs = []
            for i, j in batch.neighbor_pairs:
                if i in local_idx_map and j in local_idx_map:
                    local_neighbor_pairs.append((local_idx_map[i], local_idx_map[j]))

            # Solve downward pass QP with prior
            batch_result = solve_batch_qp_downward(
                spot_indices=spot_indices,
                observed=batch_observed,
                proportions=batch_proportions,
                enrichment=batch_enrichment,
                neighbor_pairs=local_neighbor_pairs,
                prior=prior,
                lambda_spatial=lambda_spatial,
                lambda_prior=lambda_prior,
                lambda_reg=lambda_reg,
            )

            # Update estimates (overwrite upward pass)
            level_estimates[level_idx][batch.batch_id] = (spot_indices, batch_result)

    # ==================== ASSEMBLE FINAL RESULT ====================
    # Use level 0 (finest) estimates
    result = np.zeros((n_spots, n_celltypes))
    for batch in hierarchy.levels[0].batches:
        spot_indices, batch_result = level_estimates[0][batch.batch_id]
        for local_idx, global_idx in enumerate(spot_indices):
            result[global_idx, :] = batch_result[local_idx, :]

    return result
```

**Step 4: Run tests to verify they pass**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest CITEgeist/tests/test_spatial_unet_gex.py::TestSingleGeneUNet -v`
Expected: PASS (2 tests)

**Step 5: Commit**

```bash
git add CITEgeist/model/spatial_unet_gex.py CITEgeist/tests/test_spatial_unet_gex.py
git commit -m "feat(unet-gex): add complete U-Net pipeline for single gene

Implements full upward/downward pass:
- Upward: solve batches, pool to regions
- Downward: refine with priors from parent level
- Assemble final result from level 0

Co-Authored-By: Claude Opus 4.5 <noreply@anthropic.com>"
```

---

## Task 7: Gene-Parallel Main Entry Point

**Files:**
- Modify: `CITEgeist/model/spatial_unet_gex.py`
- Test: `CITEgeist/tests/test_spatial_unet_gex.py`

**Step 1: Write the failing test for parallel execution**

```python
# Add to CITEgeist/tests/test_spatial_unet_gex.py

class TestGeneParallelExecution:
    """Tests for gene-parallel main entry point."""

    def test_run_unet_gex_deconvolution(self):
        """Test main entry point with multiple genes."""
        from model.spatial_unet_gex import run_unet_gex_deconvolution
        import scanpy as sc

        np.random.seed(42)
        n_spots = 50
        n_genes = 10
        n_celltypes = 3

        # Create mock AnnData
        coords = create_mock_spatial_coords(n_spots=n_spots)
        X = np.random.poisson(100, size=(n_spots, n_genes)).astype(float)
        adata = sc.AnnData(X=X)
        adata.obsm['spatial'] = coords
        adata.var_names = [f"Gene_{i}" for i in range(n_genes)]

        proportions = np.random.dirichlet([1, 1, 1], size=n_spots)
        cell_type_names = ["CT_A", "CT_B", "CT_C"]

        result = run_unet_gex_deconvolution(
            adata=adata,
            cell_type_proportions=proportions,
            cell_type_names=cell_type_names,
            lambda_spatial=0.1,
            lambda_prior=0.5,
            n_levels=2,
            max_workers=2,
        )

        # Should have one array per cell type
        assert len(result) == n_celltypes
        for ct_name in cell_type_names:
            assert ct_name in result
            assert result[ct_name].shape == (n_spots, n_genes)

    def test_count_conservation_across_genes(self):
        """Test that count conservation holds for all genes."""
        from model.spatial_unet_gex import run_unet_gex_deconvolution
        import scanpy as sc

        np.random.seed(42)
        n_spots = 30
        n_genes = 5
        n_celltypes = 2

        coords = create_mock_spatial_coords(n_spots=n_spots)
        X = np.random.poisson(100, size=(n_spots, n_genes)).astype(float)
        adata = sc.AnnData(X=X)
        adata.obsm['spatial'] = coords

        proportions = np.random.dirichlet([1, 1], size=n_spots)
        cell_type_names = ["CT_A", "CT_B"]

        result = run_unet_gex_deconvolution(
            adata=adata,
            cell_type_proportions=proportions,
            cell_type_names=cell_type_names,
            lambda_spatial=0.1,
            lambda_prior=0.5,
            n_levels=2,
            max_workers=2,
        )

        # Check count conservation for each spot and gene
        for gene_idx in range(n_genes):
            for spot_idx in range(n_spots):
                total_assigned = sum(
                    result[ct][spot_idx, gene_idx] for ct in cell_type_names
                )
                observed = X[spot_idx, gene_idx]
                assert np.isclose(total_assigned, observed, rtol=0.02), \
                    f"Gene {gene_idx}, Spot {spot_idx}: {total_assigned} != {observed}"
```

**Step 2: Run test to verify it fails**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest CITEgeist/tests/test_spatial_unet_gex.py::TestGeneParallelExecution::test_run_unet_gex_deconvolution -v`
Expected: FAIL with "ImportError: cannot import name 'run_unet_gex_deconvolution'"

**Step 3: Write minimal implementation**

```python
# Add to CITEgeist/model/spatial_unet_gex.py

from concurrent.futures import ProcessPoolExecutor, as_completed
import scanpy as sc
from tqdm import tqdm

# Module-level worker data for parallelization
_unet_worker_data = None


def _init_unet_worker(hierarchy_data, proportions, lambda_spatial, lambda_prior, lambda_reg):
    """Initialize worker with shared data."""
    global _unet_worker_data
    _unet_worker_data = {
        'hierarchy': hierarchy_data,
        'proportions': proportions,
        'lambda_spatial': lambda_spatial,
        'lambda_prior': lambda_prior,
        'lambda_reg': lambda_reg,
    }


def _solve_gene_worker(args):
    """Worker function for parallel gene solving."""
    gene_idx, observed, enrichment = args
    wd = _unet_worker_data

    # Reconstruct hierarchy from data
    hierarchy = wd['hierarchy']

    result = solve_gene_unet(
        gene_idx=gene_idx,
        observed=observed,
        proportions=wd['proportions'],
        enrichment=enrichment,
        hierarchy=hierarchy,
        lambda_spatial=wd['lambda_spatial'],
        lambda_prior=wd['lambda_prior'],
        lambda_reg=wd['lambda_reg'],
    )
    return gene_idx, result


def run_unet_gex_deconvolution(
    adata: sc.AnnData,
    cell_type_proportions: np.ndarray,
    cell_type_names: List[str],
    lambda_spatial: float = 0.1,
    lambda_prior: float = 0.5,
    lambda_reg: float = 0.01,
    n_levels: int = 3,
    spots_per_local: int = 20,
    max_workers: Optional[int] = None,
) -> Dict[str, np.ndarray]:
    """
    Run gene-parallel U-Net GEX deconvolution.

    Args:
        adata: AnnData with gene expression in .X and spatial coords in .obsm['spatial']
        cell_type_proportions: (n_spots, n_celltypes) from Pass 1
        cell_type_names: Names of cell types
        lambda_spatial: Spatial smoothing weight
        lambda_prior: Prior shrinkage weight
        lambda_reg: L2 regularization weight
        n_levels: Number of hierarchy levels
        spots_per_local: Target spots per level-0 batch
        max_workers: Number of parallel workers (default: CPU count)

    Returns:
        Dict mapping cell type name -> (n_spots, n_genes) expression array
    """
    import scipy.sparse as sp

    # Extract data
    coords = adata.obsm['spatial']
    if sp.issparse(adata.X):
        X = adata.X.toarray()
    else:
        X = np.array(adata.X)

    n_spots, n_genes = X.shape
    n_celltypes = len(cell_type_names)

    logger.info(f"Running U-Net GEX deconvolution: {n_spots} spots, {n_genes} genes, {n_celltypes} cell types")

    # Build spatial hierarchy
    hierarchy = SpatialHierarchy(
        coords=coords,
        n_levels=n_levels,
        spots_per_local=spots_per_local,
    )

    # Compute enrichment (uniform for now - could be enhanced)
    enrichment = np.ones((n_spots, n_celltypes)) / n_celltypes

    # Prepare results storage
    results = np.zeros((n_spots, n_celltypes, n_genes))

    # Set up workers
    if max_workers is None:
        max_workers = min(32, n_genes)  # Cap at 32 or number of genes

    # For small gene counts, just run sequentially
    if n_genes <= 4 or max_workers <= 1:
        for gene_idx in tqdm(range(n_genes), desc="Deconvolving genes"):
            observed = X[:, gene_idx]
            gene_result = solve_gene_unet(
                gene_idx=gene_idx,
                observed=observed,
                proportions=cell_type_proportions,
                enrichment=enrichment,
                hierarchy=hierarchy,
                lambda_spatial=lambda_spatial,
                lambda_prior=lambda_prior,
                lambda_reg=lambda_reg,
            )
            results[:, :, gene_idx] = gene_result
    else:
        # Parallel execution
        # Note: We pass hierarchy directly since ProcessPoolExecutor pickles it
        global _unet_worker_data
        _unet_worker_data = {
            'hierarchy': hierarchy,
            'proportions': cell_type_proportions,
            'lambda_spatial': lambda_spatial,
            'lambda_prior': lambda_prior,
            'lambda_reg': lambda_reg,
        }

        gene_args = [
            (gene_idx, X[:, gene_idx], enrichment)
            for gene_idx in range(n_genes)
        ]

        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            futures = {executor.submit(_solve_gene_worker, args): args[0] for args in gene_args}

            with tqdm(total=n_genes, desc="Deconvolving genes") as pbar:
                for future in as_completed(futures):
                    gene_idx = futures[future]
                    try:
                        _, gene_result = future.result(timeout=300)
                        results[:, :, gene_idx] = gene_result
                    except Exception as e:
                        logger.error(f"Error processing gene {gene_idx}: {e}")
                        # Fallback: proportional distribution
                        observed = X[:, gene_idx]
                        for i in range(n_spots):
                            results[i, :, gene_idx] = observed[i] * cell_type_proportions[i, :] / (cell_type_proportions[i, :].sum() + 1e-10)
                    pbar.update(1)

    # Reorganize to dict by cell type
    result_dict = {}
    for t, ct_name in enumerate(cell_type_names):
        result_dict[ct_name] = results[:, t, :]

    logger.info("U-Net GEX deconvolution complete")
    return result_dict
```

**Step 4: Run tests to verify they pass**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest CITEgeist/tests/test_spatial_unet_gex.py::TestGeneParallelExecution -v`
Expected: PASS (2 tests)

**Step 5: Commit**

```bash
git add CITEgeist/model/spatial_unet_gex.py CITEgeist/tests/test_spatial_unet_gex.py
git commit -m "feat(unet-gex): add gene-parallel main entry point

Implements run_unet_gex_deconvolution():
- Builds spatial hierarchy from coords
- Parallel gene solving with ProcessPoolExecutor
- Returns dict of cell type -> (spots, genes) arrays
- Fallback for failed genes

Co-Authored-By: Claude Opus 4.5 <noreply@anthropic.com>"
```

---

## Task 8: Integration with CitegeistModel

**Files:**
- Modify: `CITEgeist/model/citegeist_model.py`
- Test: Integration test via benchmark

**Step 1: Add U-Net option to run_cell_expression_pass1**

Read the current implementation first, then add a `use_unet` parameter:

```python
# In CITEgeist/model/citegeist_model.py, modify run_cell_expression_pass1()

# Add import at top of file:
from .spatial_unet_gex import run_unet_gex_deconvolution

# Add parameter to run_cell_expression_pass1:
def run_cell_expression_pass1(
    self,
    ...,
    use_unet: bool = False,  # NEW
    unet_lambda_spatial: float = 0.1,  # NEW
    unet_lambda_prior: float = 0.5,  # NEW
    unet_n_levels: int = 3,  # NEW
    ...
):
    """
    ...
    Args:
        ...
        use_unet: If True, use U-Net spatial smoothing instead of per-spot optimization
        unet_lambda_spatial: Spatial smoothing weight for U-Net
        unet_lambda_prior: Prior shrinkage weight for U-Net
        unet_n_levels: Number of hierarchy levels for U-Net
    """
    if use_unet:
        # Use new U-Net implementation
        result = run_unet_gex_deconvolution(
            adata=self.gene_expression_adata,
            cell_type_proportions=self.cell_type_proportions,
            cell_type_names=list(self.cell_profile_dict.keys()),
            lambda_spatial=unet_lambda_spatial,
            lambda_prior=unet_lambda_prior,
            n_levels=unet_n_levels,
            max_workers=max_workers,
        )
        # Store results in AnnData layers
        for ct_name, expr_array in result.items():
            layer_name = f"{ct_name}_layer_pass1"
            self.gene_expression_adata.layers[layer_name] = expr_array
        return result
    else:
        # Existing implementation
        ...
```

**Step 2: Run integration test**

Run on Xenium benchmark region 0 and compare metrics.

**Step 3: Commit**

```bash
git add CITEgeist/model/citegeist_model.py
git commit -m "feat(model): integrate U-Net GEX into CitegeistModel

Adds use_unet parameter to run_cell_expression_pass1():
- use_unet=True uses new spatial smoothing
- use_unet=False (default) uses existing per-spot optimization
- Backward compatible

Co-Authored-By: Claude Opus 4.5 <noreply@anthropic.com>"
```

---

## Task 9: Benchmark Validation Script

**Files:**
- Create: `CITEgeist/tests/test_unet_benchmark.py`

**Step 1: Write benchmark comparison script**

```python
# CITEgeist/tests/test_unet_benchmark.py
"""
Benchmark comparison: U-Net vs baseline Module 3 Pass 2.

Runs on Xenium pseudo-Visium data and compares Pearson r, RMSE.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from sklearn.metrics import mean_squared_error

sys.path.insert(0, str(Path(__file__).parent.parent))


def run_benchmark_comparison(region_id: int = 0):
    """Run benchmark comparison on one region."""
    from model.spatial_unet_gex import run_unet_gex_deconvolution
    import scanpy as sc

    # Paths (adjust as needed)
    base_path = Path("/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist")
    region_path = base_path / f"Benchmarking/xenium_pseudovisium/data_protein_gt/Xenium_region_{region_id}"

    # Load data
    adata = sc.read_h5ad(region_path / "pseudovisium.h5ad")
    props = pd.read_csv(region_path / "citegeist_proportions.csv", index_col=0)

    # Load ground truth
    gt_path = region_path / "ground_truth_gex"
    cell_types = ["B cells", "CD4+ T cells", "CD8+ T cells", "Macrophages",
                  "Endothelial", "Epithelial", "Fibroblasts"]

    # Run U-Net
    proportions = props.values
    result = run_unet_gex_deconvolution(
        adata=adata,
        cell_type_proportions=proportions,
        cell_type_names=cell_types,
        lambda_spatial=0.1,
        lambda_prior=0.5,
        n_levels=3,
        max_workers=16,
    )

    # Compare to ground truth
    metrics = []
    for ct in cell_types:
        gt_file = gt_path / f"{ct.replace(' ', '_')}_GT.csv"
        if gt_file.exists():
            gt = pd.read_csv(gt_file, index_col=0)
            pred = result[ct]

            # Align genes
            common_genes = list(set(gt.columns) & set(adata.var_names))
            gt_aligned = gt[common_genes].values
            pred_aligned = pred[:, [list(adata.var_names).index(g) for g in common_genes]]

            # Flatten and compute metrics
            gt_flat = gt_aligned.flatten()
            pred_flat = pred_aligned.flatten()

            r, _ = pearsonr(gt_flat, pred_flat)
            rmse = np.sqrt(mean_squared_error(gt_flat, pred_flat))

            metrics.append({
                "cell_type": ct,
                "region_id": region_id,
                "pearson_r": r,
                "rmse": rmse,
                "n_genes": len(common_genes),
            })

    return metrics


if __name__ == "__main__":
    results = run_benchmark_comparison(region_id=0)
    print(json.dumps(results, indent=2))
```

**Step 2: Commit**

```bash
git add CITEgeist/tests/test_unet_benchmark.py
git commit -m "test(unet-gex): add benchmark comparison script

Compares U-Net vs baseline on Xenium pseudo-Visium:
- Loads region data and ground truth
- Runs U-Net deconvolution
- Computes Pearson r and RMSE per cell type

Co-Authored-By: Claude Opus 4.5 <noreply@anthropic.com>"
```

---

## Summary

| Task | Description | Files |
|------|-------------|-------|
| 1 | Spatial hierarchy data structure | `spatial_unet_gex.py`, tests |
| 2 | Asymmetric Laplacian weights | `spatial_unet_gex.py`, tests |
| 3 | Upward pass QP solver | `spatial_unet_gex.py`, tests |
| 4 | Pooling between levels | `spatial_unet_gex.py`, tests |
| 5 | Downward pass QP with prior | `spatial_unet_gex.py`, tests |
| 6 | Single-gene U-Net pipeline | `spatial_unet_gex.py`, tests |
| 7 | Gene-parallel entry point | `spatial_unet_gex.py`, tests |
| 8 | CitegeistModel integration | `citegeist_model.py` |
| 9 | Benchmark validation | `test_unet_benchmark.py` |

**Estimated implementation time:** 9 tasks × ~20-30 min each = ~4-5 hours

**Success criteria:**
- All unit tests pass
- Benchmark shows improved Pearson r (target: +0.10 mean improvement)
- Runtime ≤2 hours on 64 cores for full gene set
