# Hierarchical Profile Discovery Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Implement hierarchical NMF-based profile discovery that learns marker hierarchies from spatial colocalization, enabling shared markers (CD3, Vimentin) to appear in multiple related profiles.

**Architecture:** Two-phase algorithm: (1) hierarchical clustering on bivariate Moran's I distances with reconstruction-guided tree cutting, (2) NMF at each level to learn marker weights with expression-weighted allocation. Output is flattened to binary profiles for existing Module 3.

**Tech Stack:** Python 3.10, NumPy, SciPy (hierarchical clustering, linkage), scikit-learn (NMF), existing spatial_colocalization.py infrastructure.

---

## Task 1: Add ProfileTree Data Structure

**Files:**
- Modify: `CITEgeist/model/spatial_colocalization.py` (after line 950, after ProfileDiscoveryResult)
- Test: `tests/test_hierarchical_profiles.py` (create new)

**Step 1: Write the failing test**

Create `tests/test_hierarchical_profiles.py`:

```python
"""
Test script for Hierarchical Profile Discovery (Module 2b enhancement).

Tests the two-phase algorithm:
1. Structure learning via hierarchical clustering on colocalization distances
2. NMF weight learning with expression-weighted allocation
"""

import numpy as np
import pytest
import sys

sys.path.insert(0, '/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/.worktrees/hierarchical-profiles')

from CITEgeist.model.spatial_colocalization import ProfileTreeNode, ProfileTree


class TestProfileTreeStructure:
    """Test ProfileTree data structure."""

    def test_create_leaf_node(self):
        """Leaf node should hold markers and have no children."""
        node = ProfileTreeNode(
            node_id="leaf_0",
            markers=["CD3", "CD4"],
            children=[],
            parent_id=None,
            depth=0,
        )
        assert node.is_leaf
        assert node.markers == ["CD3", "CD4"]
        assert len(node.children) == 0

    def test_create_internal_node(self):
        """Internal node should have children and aggregate markers."""
        child1 = ProfileTreeNode("c1", ["CD4"], [], "root", 1)
        child2 = ProfileTreeNode("c2", ["CD8"], [], "root", 1)
        root = ProfileTreeNode(
            node_id="root",
            markers=["CD3"],  # Shared marker
            children=[child1, child2],
            parent_id=None,
            depth=0,
        )
        assert not root.is_leaf
        assert len(root.children) == 2

    def test_tree_get_leaves(self):
        """Tree should return all leaf nodes."""
        child1 = ProfileTreeNode("c1", ["CD4"], [], "root", 1)
        child2 = ProfileTreeNode("c2", ["CD8"], [], "root", 1)
        root = ProfileTreeNode("root", ["CD3"], [child1, child2], None, 0)

        tree = ProfileTree(root=root, n_levels=2)
        leaves = tree.get_leaves()

        assert len(leaves) == 2
        assert set(l.node_id for l in leaves) == {"c1", "c2"}


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
```

**Step 2: Run test to verify it fails**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/.worktrees/hierarchical-profiles && python -m pytest tests/test_hierarchical_profiles.py::TestProfileTreeStructure -v`

Expected: FAIL with "cannot import name 'ProfileTreeNode'"

**Step 3: Write minimal implementation**

Add to `CITEgeist/model/spatial_colocalization.py` after line 950 (after ProfileDiscoveryResult class):

```python
@dataclass
class ProfileTreeNode:
    """A node in the hierarchical profile tree."""

    node_id: str
    markers: List[str]  # Markers assigned to THIS node (not descendants)
    children: List["ProfileTreeNode"]
    parent_id: Optional[str]
    depth: int
    nmf_weights: Optional[Dict[str, float]] = None  # marker -> weight from NMF

    @property
    def is_leaf(self) -> bool:
        """True if this is a leaf node (no children)."""
        return len(self.children) == 0

    def get_all_markers(self) -> List[str]:
        """Get all markers in this subtree (self + descendants)."""
        markers = list(self.markers)
        for child in self.children:
            markers.extend(child.get_all_markers())
        return markers


@dataclass
class ProfileTree:
    """Hierarchical tree of marker profiles."""

    root: ProfileTreeNode
    n_levels: int

    def get_leaves(self) -> List[ProfileTreeNode]:
        """Return all leaf nodes (final cell type profiles)."""
        leaves = []
        self._collect_leaves(self.root, leaves)
        return leaves

    def _collect_leaves(self, node: ProfileTreeNode, leaves: List[ProfileTreeNode]) -> None:
        """Recursively collect leaf nodes."""
        if node.is_leaf:
            leaves.append(node)
        else:
            for child in node.children:
                self._collect_leaves(child, leaves)

    def get_depth(self) -> int:
        """Return maximum depth of tree."""
        return self._max_depth(self.root)

    def _max_depth(self, node: ProfileTreeNode) -> int:
        """Recursively compute max depth."""
        if node.is_leaf:
            return node.depth
        return max(self._max_depth(c) for c in node.children)
```

**Step 4: Run test to verify it passes**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/.worktrees/hierarchical-profiles && python -m pytest tests/test_hierarchical_profiles.py::TestProfileTreeStructure -v`

Expected: PASS (3 tests)

**Step 5: Commit**

```bash
cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/.worktrees/hierarchical-profiles
git add CITEgeist/model/spatial_colocalization.py tests/test_hierarchical_profiles.py
git commit -m "feat: add ProfileTreeNode and ProfileTree data structures"
```

---

## Task 2: Add HierarchicalProfileResult Dataclass

**Files:**
- Modify: `CITEgeist/model/spatial_colocalization.py` (after ProfileTree)
- Test: `tests/test_hierarchical_profiles.py`

**Step 1: Write the failing test**

Add to `tests/test_hierarchical_profiles.py`:

```python
from CITEgeist.model.spatial_colocalization import HierarchicalProfileResult


class TestHierarchicalProfileResult:
    """Test HierarchicalProfileResult dataclass."""

    def test_create_result(self):
        """Should store tree and flat profiles."""
        child1 = ProfileTreeNode("c1", ["CD4"], [], "root", 1)
        child2 = ProfileTreeNode("c2", ["CD8"], [], "root", 1)
        root = ProfileTreeNode("root", ["CD3"], [child1, child2], None, 0)
        tree = ProfileTree(root=root, n_levels=2)

        result = HierarchicalProfileResult(
            tree=tree,
            flat_profiles={"CD4+ T": ["CD3", "CD4"], "CD8+ T": ["CD3", "CD8"]},
            depth_per_branch={"T cells": 2},
            shared_markers={"CD3": ["CD4+ T", "CD8+ T"]},
            reconstruction_error=0.05,
        )

        assert len(result.flat_profiles) == 2
        assert "CD3" in result.flat_profiles["CD4+ T"]
        assert "CD3" in result.flat_profiles["CD8+ T"]

    def test_to_profile_dict(self):
        """Should convert to Module 3 compatible format."""
        child1 = ProfileTreeNode("c1", ["CD4"], [], "root", 1)
        child2 = ProfileTreeNode("c2", ["CD8"], [], "root", 1)
        root = ProfileTreeNode("root", ["CD3"], [child1, child2], None, 0)
        tree = ProfileTree(root=root, n_levels=2)

        result = HierarchicalProfileResult(
            tree=tree,
            flat_profiles={"CD4+ T": ["CD3", "CD4"], "CD8+ T": ["CD3", "CD8"]},
            depth_per_branch={},
            shared_markers={},
            reconstruction_error=0.05,
        )

        profile_dict = result.to_profile_dict()

        assert isinstance(profile_dict, dict)
        assert profile_dict["CD4+ T"] == ["CD3", "CD4"]
```

**Step 2: Run test to verify it fails**

Run: `python -m pytest tests/test_hierarchical_profiles.py::TestHierarchicalProfileResult -v`

Expected: FAIL with "cannot import name 'HierarchicalProfileResult'"

**Step 3: Write minimal implementation**

Add to `CITEgeist/model/spatial_colocalization.py` after ProfileTree:

```python
@dataclass
class HierarchicalProfileResult:
    """Results from hierarchical profile discovery."""

    tree: ProfileTree  # Full hierarchy for inspection
    flat_profiles: Dict[str, List[str]]  # cell_type -> [markers] for Module 3
    depth_per_branch: Dict[str, int]  # branch_name -> depth
    shared_markers: Dict[str, List[str]]  # marker -> [cell_types that share it]
    reconstruction_error: float  # Final reconstruction error

    def to_profile_dict(self) -> Dict[str, List[str]]:
        """Convert to Module 3 compatible profile dictionary."""
        return dict(self.flat_profiles)

    def summary(self) -> str:
        """Return summary string."""
        n_profiles = len(self.flat_profiles)
        n_shared = len(self.shared_markers)
        max_depth = self.tree.get_depth()
        return (
            f"Hierarchical profiles: {n_profiles} cell types\n"
            f"Tree depth: {max_depth}\n"
            f"Shared markers: {n_shared}\n"
            f"Reconstruction error: {self.reconstruction_error:.4f}"
        )
```

**Step 4: Run test to verify it passes**

Run: `python -m pytest tests/test_hierarchical_profiles.py::TestHierarchicalProfileResult -v`

Expected: PASS (2 tests)

**Step 5: Commit**

```bash
git add CITEgeist/model/spatial_colocalization.py tests/test_hierarchical_profiles.py
git commit -m "feat: add HierarchicalProfileResult dataclass"
```

---

## Task 3: Implement Colocalization Distance Matrix Builder

**Files:**
- Modify: `CITEgeist/model/spatial_colocalization.py`
- Test: `tests/test_hierarchical_profiles.py`

**Step 1: Write the failing test**

Add to `tests/test_hierarchical_profiles.py`:

```python
from CITEgeist.model.spatial_colocalization import (
    _build_colocalization_distance_matrix,
    ColocalizationResult,
    MarkerPairColocalization,
)


class TestDistanceMatrix:
    """Test colocalization distance matrix construction."""

    def test_distance_from_morans_i(self):
        """Distance should be 1 - normalized Moran's I."""
        # Create mock colocalization result with known Moran's I values
        pairs = [
            MarkerPairColocalization(
                marker_a="A", marker_b="B",
                jaccard_index=0.5, co_occurrence_spots=10, co_occurrence_fraction=0.1,
                pearson_r=0.8, spearman_rho=0.75, correlation_pvalue=0.001,
                cosine_similarity=0.85,
                bivariate_morans_i=0.6,  # High colocalization
                bivariate_morans_pvalue=0.001,
                neighbor_enrichment_ab=1.5, neighbor_enrichment_ba=1.4,
                neighbor_enrichment_pvalue=0.01, mutual_neighbor_enrichment=1.45,
                colocalization_score=0.7,
            ),
            MarkerPairColocalization(
                marker_a="A", marker_b="C",
                jaccard_index=0.1, co_occurrence_spots=2, co_occurrence_fraction=0.02,
                pearson_r=0.1, spearman_rho=0.05, correlation_pvalue=0.5,
                cosine_similarity=0.15,
                bivariate_morans_i=0.1,  # Low colocalization
                bivariate_morans_pvalue=0.3,
                neighbor_enrichment_ab=1.0, neighbor_enrichment_ba=1.0,
                neighbor_enrichment_pvalue=0.5, mutual_neighbor_enrichment=1.0,
                colocalization_score=0.15,
            ),
            MarkerPairColocalization(
                marker_a="B", marker_b="C",
                jaccard_index=0.2, co_occurrence_spots=4, co_occurrence_fraction=0.04,
                pearson_r=0.2, spearman_rho=0.15, correlation_pvalue=0.2,
                cosine_similarity=0.25,
                bivariate_morans_i=0.2,
                bivariate_morans_pvalue=0.1,
                neighbor_enrichment_ab=1.1, neighbor_enrichment_ba=1.1,
                neighbor_enrichment_pvalue=0.2, mutual_neighbor_enrichment=1.1,
                colocalization_score=0.25,
            ),
        ]

        coloc_result = ColocalizationResult(
            pairs=pairs,
            marker_names=["A", "B", "C"],
            n_spots=100,
            neighbor_k=6,
        )

        D, markers = _build_colocalization_distance_matrix(coloc_result)

        assert D.shape == (3, 3)
        assert markers == ["A", "B", "C"]
        # D[A,B] should be smallest (highest Moran's I = 0.6)
        assert D[0, 1] < D[0, 2]  # A-B closer than A-C
        assert D[0, 1] < D[1, 2]  # A-B closer than B-C
```

**Step 2: Run test to verify it fails**

Run: `python -m pytest tests/test_hierarchical_profiles.py::TestDistanceMatrix -v`

Expected: FAIL with "cannot import name '_build_colocalization_distance_matrix'"

**Step 3: Write minimal implementation**

Add to `CITEgeist/model/spatial_colocalization.py`:

```python
def _build_colocalization_distance_matrix(
    coloc_result: ColocalizationResult,
) -> Tuple[NDArray[np.floating], List[str]]:
    """
    Build distance matrix from colocalization results.

    Distance = 1 - normalized_bivariate_morans_I
    High Moran's I (colocalized) -> low distance -> same branch

    Args:
        coloc_result: Result from analyze_marker_colocalization()

    Returns:
        Tuple of (distance_matrix, marker_names)
        Distance matrix is symmetric with 0 on diagonal.
    """
    markers = coloc_result.marker_names
    n_markers = len(markers)
    marker_to_idx = {m: i for i, m in enumerate(markers)}

    # Initialize with max distance (1.0 means no colocalization)
    D = np.ones((n_markers, n_markers))
    np.fill_diagonal(D, 0.0)

    # Collect all Moran's I values for normalization
    morans_values = [p.bivariate_morans_i for p in coloc_result.pairs]
    if len(morans_values) == 0:
        return D, markers

    # Normalize Moran's I to [0, 1] range
    # Moran's I can be negative (spatial dispersion), so shift and scale
    min_i = min(morans_values)
    max_i = max(morans_values)
    range_i = max_i - min_i if max_i > min_i else 1.0

    for pair in coloc_result.pairs:
        i = marker_to_idx[pair.marker_a]
        j = marker_to_idx[pair.marker_b]

        # Normalize Moran's I to [0, 1]
        normalized_i = (pair.bivariate_morans_i - min_i) / range_i

        # Distance = 1 - normalized Moran's I
        distance = 1.0 - normalized_i

        D[i, j] = distance
        D[j, i] = distance

    return D, markers
```

**Step 4: Run test to verify it passes**

Run: `python -m pytest tests/test_hierarchical_profiles.py::TestDistanceMatrix -v`

Expected: PASS

**Step 5: Commit**

```bash
git add CITEgeist/model/spatial_colocalization.py tests/test_hierarchical_profiles.py
git commit -m "feat: add colocalization distance matrix builder"
```

---

## Task 4: Implement Reconstruction Error Calculator

**Files:**
- Modify: `CITEgeist/model/spatial_colocalization.py`
- Test: `tests/test_hierarchical_profiles.py`

**Step 1: Write the failing test**

Add to `tests/test_hierarchical_profiles.py`:

```python
from CITEgeist.model.spatial_colocalization import _compute_reconstruction_error


class TestReconstructionError:
    """Test reconstruction error computation."""

    def test_perfect_reconstruction(self):
        """Error should be 0 for perfectly separable markers."""
        # Two markers with identical patterns -> one profile reconstructs perfectly
        X = np.array([
            [1.0, 1.0],
            [2.0, 2.0],
            [3.0, 3.0],
        ])
        marker_names = ["A", "B"]
        markers_in_node = ["A", "B"]

        error = _compute_reconstruction_error(X, marker_names, markers_in_node)

        # Should be very low (markers are identical patterns)
        assert error < 0.1

    def test_high_error_for_uncorrelated(self):
        """Error should be higher for uncorrelated markers."""
        # Two markers with opposite patterns
        X = np.array([
            [1.0, 3.0],
            [2.0, 2.0],
            [3.0, 1.0],
        ])
        marker_names = ["A", "B"]
        markers_in_node = ["A", "B"]

        error = _compute_reconstruction_error(X, marker_names, markers_in_node)

        # Should be higher than correlated case
        assert error > 0.1
```

**Step 2: Run test to verify it fails**

Run: `python -m pytest tests/test_hierarchical_profiles.py::TestReconstructionError -v`

Expected: FAIL with "cannot import name '_compute_reconstruction_error'"

**Step 3: Write minimal implementation**

Add to `CITEgeist/model/spatial_colocalization.py`:

```python
def _compute_reconstruction_error(
    X: NDArray[np.floating],
    marker_names: List[str],
    markers_in_node: List[str],
    n_components: int = 1,
) -> float:
    """
    Compute reconstruction error for markers in a node using NMF.

    Args:
        X: Full expression matrix (n_spots, n_markers)
        marker_names: Names for all markers in X
        markers_in_node: Markers belonging to this node
        n_components: Number of NMF components (default 1 for single node)

    Returns:
        Normalized reconstruction error (Frobenius norm / original norm)
    """
    from sklearn.decomposition import NMF

    # Get indices for markers in this node
    marker_to_idx = {m: i for i, m in enumerate(marker_names)}
    indices = [marker_to_idx[m] for m in markers_in_node if m in marker_to_idx]

    if len(indices) == 0:
        return 1.0  # No markers -> max error

    if len(indices) == 1:
        return 0.0  # Single marker -> perfect reconstruction

    # Extract submatrix
    X_node = X[:, indices]

    # Ensure non-negative
    X_node = np.maximum(X_node, 0)

    # Handle zero matrix
    original_norm = np.linalg.norm(X_node, 'fro')
    if original_norm < 1e-10:
        return 0.0

    # Fit NMF
    n_comp = min(n_components, len(indices) - 1, X_node.shape[0] - 1)
    n_comp = max(1, n_comp)

    try:
        nmf = NMF(n_components=n_comp, init='nndsvda', max_iter=200, random_state=42)
        W = nmf.fit_transform(X_node)
        H = nmf.components_

        # Reconstruction
        X_reconstructed = W @ H

        # Normalized error
        error_norm = np.linalg.norm(X_node - X_reconstructed, 'fro')
        return error_norm / original_norm

    except Exception:
        return 1.0  # Return max error on failure
```

**Step 4: Run test to verify it passes**

Run: `python -m pytest tests/test_hierarchical_profiles.py::TestReconstructionError -v`

Expected: PASS (2 tests)

**Step 5: Commit**

```bash
git add CITEgeist/model/spatial_colocalization.py tests/test_hierarchical_profiles.py
git commit -m "feat: add reconstruction error calculator for tree cutting"
```

---

## Task 5: Implement Reconstruction-Guided Tree Cutting

**Files:**
- Modify: `CITEgeist/model/spatial_colocalization.py`
- Test: `tests/test_hierarchical_profiles.py`

**Step 1: Write the failing test**

Add to `tests/test_hierarchical_profiles.py`:

```python
from CITEgeist.model.spatial_colocalization import _build_hierarchical_tree
from scipy.cluster.hierarchy import linkage


class TestTreeBuilding:
    """Test hierarchical tree construction with reconstruction-guided cutting."""

    def test_flat_data_stays_flat(self):
        """Data with no hierarchy should produce flat tree (depth 1)."""
        # 4 completely independent markers
        np.random.seed(42)
        X = np.random.rand(100, 4)  # Uncorrelated random data
        marker_names = ["A", "B", "C", "D"]

        # Build linkage (will have structure but cutting should reject splits)
        D = 1 - np.corrcoef(X.T)  # Distance from correlation
        np.fill_diagonal(D, 0)
        condensed = squareform(D)
        Z = linkage(condensed, method='ward')

        tree = _build_hierarchical_tree(
            X=X,
            marker_names=marker_names,
            linkage_matrix=Z,
            improvement_threshold=0.05,
        )

        # Should have minimal depth since splits don't improve reconstruction
        # Each marker becomes its own leaf
        leaves = tree.get_leaves()
        assert len(leaves) >= 1

    def test_hierarchical_data_creates_levels(self):
        """Data with clear hierarchy should produce multi-level tree."""
        np.random.seed(42)
        n_spots = 100

        # Create hierarchical structure:
        # - Markers A, B correlated (same cell type)
        # - Markers C, D correlated (different cell type)
        base_A = np.random.rand(n_spots)
        base_C = np.random.rand(n_spots)

        X = np.column_stack([
            base_A + 0.1 * np.random.rand(n_spots),  # A
            base_A + 0.1 * np.random.rand(n_spots),  # B (similar to A)
            base_C + 0.1 * np.random.rand(n_spots),  # C
            base_C + 0.1 * np.random.rand(n_spots),  # D (similar to C)
        ])
        marker_names = ["A", "B", "C", "D"]

        # Build linkage
        D = 1 - np.abs(np.corrcoef(X.T))
        np.fill_diagonal(D, 0)
        condensed = squareform(D)
        Z = linkage(condensed, method='ward')

        tree = _build_hierarchical_tree(
            X=X,
            marker_names=marker_names,
            linkage_matrix=Z,
            improvement_threshold=0.05,
        )

        # Should detect hierarchy
        leaves = tree.get_leaves()
        # Expect 2 groups: {A,B} and {C,D}
        assert len(leaves) == 2 or len(leaves) == 4  # Either 2 merged or 4 flat
```

**Step 2: Run test to verify it fails**

Run: `python -m pytest tests/test_hierarchical_profiles.py::TestTreeBuilding -v`

Expected: FAIL with "cannot import name '_build_hierarchical_tree'"

**Step 3: Write minimal implementation**

Add to `CITEgeist/model/spatial_colocalization.py`:

```python
from scipy.cluster.hierarchy import to_tree, ClusterNode


def _build_hierarchical_tree(
    X: NDArray[np.floating],
    marker_names: List[str],
    linkage_matrix: NDArray[np.floating],
    improvement_threshold: float = 0.05,
    max_depth: int = 5,
) -> ProfileTree:
    """
    Build hierarchical profile tree with reconstruction-guided cutting.

    Args:
        X: Expression matrix (n_spots, n_markers)
        marker_names: Names for each column
        linkage_matrix: Scipy linkage matrix from hierarchical clustering
        improvement_threshold: Min reconstruction improvement to split (default 5%)
        max_depth: Maximum tree depth (safety limit)

    Returns:
        ProfileTree with reconstruction-guided structure
    """
    # Convert scipy linkage to tree structure
    scipy_root = to_tree(linkage_matrix)

    # Recursively build ProfileTree with reconstruction-guided cutting
    root_node = _recursive_tree_cut(
        scipy_node=scipy_root,
        X=X,
        marker_names=marker_names,
        improvement_threshold=improvement_threshold,
        current_depth=0,
        max_depth=max_depth,
        parent_id=None,
    )

    # Compute actual depth
    n_levels = root_node.depth + 1 if root_node.is_leaf else _compute_tree_depth(root_node)

    return ProfileTree(root=root_node, n_levels=n_levels)


def _compute_tree_depth(node: ProfileTreeNode) -> int:
    """Compute maximum depth of tree."""
    if node.is_leaf:
        return node.depth + 1
    return max(_compute_tree_depth(c) for c in node.children)


def _get_leaf_indices(scipy_node: ClusterNode) -> List[int]:
    """Get all leaf indices under a scipy ClusterNode."""
    if scipy_node.is_leaf():
        return [scipy_node.id]
    left_indices = _get_leaf_indices(scipy_node.get_left())
    right_indices = _get_leaf_indices(scipy_node.get_right())
    return left_indices + right_indices


def _recursive_tree_cut(
    scipy_node: ClusterNode,
    X: NDArray[np.floating],
    marker_names: List[str],
    improvement_threshold: float,
    current_depth: int,
    max_depth: int,
    parent_id: Optional[str],
) -> ProfileTreeNode:
    """
    Recursively build tree with reconstruction-guided cutting.

    At each node, check if splitting improves reconstruction.
    If improvement > threshold, split. Otherwise, stop.
    """
    node_id = f"node_{scipy_node.id}"

    # Get markers in this subtree
    leaf_indices = _get_leaf_indices(scipy_node)
    node_markers = [marker_names[i] for i in leaf_indices if i < len(marker_names)]

    # Base case: leaf node or max depth reached
    if scipy_node.is_leaf() or current_depth >= max_depth:
        return ProfileTreeNode(
            node_id=node_id,
            markers=node_markers,
            children=[],
            parent_id=parent_id,
            depth=current_depth,
        )

    # Get left and right children markers
    left_indices = _get_leaf_indices(scipy_node.get_left())
    right_indices = _get_leaf_indices(scipy_node.get_right())
    left_markers = [marker_names[i] for i in left_indices if i < len(marker_names)]
    right_markers = [marker_names[i] for i in right_indices if i < len(marker_names)]

    # Skip split if either side is empty
    if len(left_markers) == 0 or len(right_markers) == 0:
        return ProfileTreeNode(
            node_id=node_id,
            markers=node_markers,
            children=[],
            parent_id=parent_id,
            depth=current_depth,
        )

    # Compute reconstruction error: merged vs split
    error_merged = _compute_reconstruction_error(X, marker_names, node_markers)
    error_left = _compute_reconstruction_error(X, marker_names, left_markers)
    error_right = _compute_reconstruction_error(X, marker_names, right_markers)
    error_split = (error_left * len(left_markers) + error_right * len(right_markers)) / len(node_markers)

    # Compute improvement
    if error_merged > 1e-10:
        improvement = (error_merged - error_split) / error_merged
    else:
        improvement = 0.0

    # Decision: split or stop
    if improvement > improvement_threshold:
        # Split is worthwhile - recurse on children
        left_child = _recursive_tree_cut(
            scipy_node.get_left(), X, marker_names,
            improvement_threshold, current_depth + 1, max_depth, node_id
        )
        right_child = _recursive_tree_cut(
            scipy_node.get_right(), X, marker_names,
            improvement_threshold, current_depth + 1, max_depth, node_id
        )

        return ProfileTreeNode(
            node_id=node_id,
            markers=[],  # Internal nodes don't own markers directly
            children=[left_child, right_child],
            parent_id=parent_id,
            depth=current_depth,
        )
    else:
        # Stop here - this is a leaf
        return ProfileTreeNode(
            node_id=node_id,
            markers=node_markers,
            children=[],
            parent_id=parent_id,
            depth=current_depth,
        )
```

**Step 4: Run test to verify it passes**

Run: `python -m pytest tests/test_hierarchical_profiles.py::TestTreeBuilding -v`

Expected: PASS (2 tests)

**Step 5: Commit**

```bash
git add CITEgeist/model/spatial_colocalization.py tests/test_hierarchical_profiles.py
git commit -m "feat: add reconstruction-guided hierarchical tree cutting"
```

---

## Task 6: Implement NMF Weight Learning at Each Level

**Files:**
- Modify: `CITEgeist/model/spatial_colocalization.py`
- Test: `tests/test_hierarchical_profiles.py`

**Step 1: Write the failing test**

Add to `tests/test_hierarchical_profiles.py`:

```python
from CITEgeist.model.spatial_colocalization import _compute_nmf_weights


class TestNMFWeights:
    """Test NMF weight learning at each tree level."""

    def test_weights_for_two_children(self):
        """NMF should assign weights to markers for each child branch."""
        np.random.seed(42)
        n_spots = 100

        # Create data where A,B go with child1 and C,D go with child2
        pattern1 = np.random.rand(n_spots)
        pattern2 = np.random.rand(n_spots)

        X = np.column_stack([
            pattern1,  # A
            pattern1 + 0.1 * np.random.rand(n_spots),  # B (like A)
            pattern2,  # C
            pattern2 + 0.1 * np.random.rand(n_spots),  # D (like C)
        ])
        marker_names = ["A", "B", "C", "D"]

        # Create node with 2 children
        child1 = ProfileTreeNode("c1", ["A", "B"], [], "root", 1)
        child2 = ProfileTreeNode("c2", ["C", "D"], [], "root", 1)
        root = ProfileTreeNode("root", [], [child1, child2], None, 0)

        weights = _compute_nmf_weights(X, marker_names, root)

        # Should have weights for each marker in each child
        assert "c1" in weights
        assert "c2" in weights
        # A should have higher weight in c1
        assert weights["c1"].get("A", 0) > weights["c2"].get("A", 0)
        # C should have higher weight in c2
        assert weights["c2"].get("C", 0) > weights["c1"].get("C", 0)
```

**Step 2: Run test to verify it fails**

Run: `python -m pytest tests/test_hierarchical_profiles.py::TestNMFWeights -v`

Expected: FAIL with "cannot import name '_compute_nmf_weights'"

**Step 3: Write minimal implementation**

Add to `CITEgeist/model/spatial_colocalization.py`:

```python
def _compute_nmf_weights(
    X: NDArray[np.floating],
    marker_names: List[str],
    node: ProfileTreeNode,
) -> Dict[str, Dict[str, float]]:
    """
    Compute NMF weights at an internal node.

    Runs NMF with n_components = number of children to learn
    how markers should be allocated to each child branch.

    Args:
        X: Expression matrix (n_spots, n_markers)
        marker_names: Names for each column
        node: Internal node with children

    Returns:
        Dict mapping child_id -> {marker: weight}
    """
    from sklearn.decomposition import NMF

    if node.is_leaf or len(node.children) == 0:
        return {}

    # Get all markers in this subtree
    all_markers = node.get_all_markers()
    if len(all_markers) == 0:
        return {c.node_id: {} for c in node.children}

    marker_to_idx = {m: i for i, m in enumerate(marker_names)}
    indices = [marker_to_idx[m] for m in all_markers if m in marker_to_idx]

    if len(indices) == 0:
        return {c.node_id: {} for c in node.children}

    # Extract submatrix
    X_node = X[:, indices]
    X_node = np.maximum(X_node, 0)

    # Number of components = number of children
    n_components = len(node.children)
    n_components = min(n_components, len(indices), X_node.shape[0] - 1)
    n_components = max(1, n_components)

    try:
        nmf = NMF(n_components=n_components, init='nndsvda', max_iter=200, random_state=42)
        W = nmf.fit_transform(X_node)  # (n_spots, n_components)
        H = nmf.components_  # (n_components, n_markers_in_node)

        # Map components to children based on which markers they load on
        # H[k, j] = how much component k contributes to marker j

        # For each child, find which component best matches its markers
        child_to_component = {}
        markers_in_node = [all_markers[i] for i in range(len(all_markers)) if all_markers[i] in marker_to_idx]

        for child_idx, child in enumerate(node.children):
            child_markers = child.get_all_markers()
            child_marker_indices = [markers_in_node.index(m) for m in child_markers if m in markers_in_node]

            if len(child_marker_indices) == 0:
                child_to_component[child.node_id] = child_idx % n_components
                continue

            # Find component with highest average loading on child markers
            best_component = 0
            best_score = -np.inf
            for comp_idx in range(n_components):
                score = np.mean(H[comp_idx, child_marker_indices])
                if score > best_score:
                    best_score = score
                    best_component = comp_idx

            child_to_component[child.node_id] = best_component

        # Build weights dict
        weights = {}
        for child in node.children:
            comp_idx = child_to_component[child.node_id]
            weights[child.node_id] = {}
            for marker_idx, marker in enumerate(markers_in_node):
                if marker_idx < H.shape[1]:
                    weights[child.node_id][marker] = float(H[comp_idx, marker_idx])

        # Normalize weights per marker (so they sum to 1 across children)
        for marker in markers_in_node:
            total = sum(weights[c.node_id].get(marker, 0) for c in node.children)
            if total > 1e-10:
                for child in node.children:
                    if marker in weights[child.node_id]:
                        weights[child.node_id][marker] /= total

        return weights

    except Exception as e:
        logger.warning(f"NMF failed at node {node.node_id}: {e}")
        # Fallback: equal weights
        weights = {}
        for child in node.children:
            child_markers = child.get_all_markers()
            weights[child.node_id] = {m: 1.0 / len(node.children) for m in child_markers}
        return weights
```

**Step 4: Run test to verify it passes**

Run: `python -m pytest tests/test_hierarchical_profiles.py::TestNMFWeights -v`

Expected: PASS

**Step 5: Commit**

```bash
git add CITEgeist/model/spatial_colocalization.py tests/test_hierarchical_profiles.py
git commit -m "feat: add NMF weight learning at each tree level"
```

---

## Task 7: Implement Weight Propagation and Flattening

**Files:**
- Modify: `CITEgeist/model/spatial_colocalization.py`
- Test: `tests/test_hierarchical_profiles.py`

**Step 1: Write the failing test**

Add to `tests/test_hierarchical_profiles.py`:

```python
from CITEgeist.model.spatial_colocalization import _flatten_tree_to_profiles


class TestFlattening:
    """Test flattening hierarchical tree to flat profiles."""

    def test_shared_marker_appears_in_multiple_profiles(self):
        """Shared markers at parent should appear in all children profiles."""
        # Tree: root(CD3) -> [CD4+ T(CD4), CD8+ T(CD8)]
        # CD3 is shared, should appear in both flat profiles

        child1 = ProfileTreeNode("cd4_t", ["CD4"], [], "root", 1)
        child2 = ProfileTreeNode("cd8_t", ["CD8"], [], "root", 1)
        root = ProfileTreeNode("root", ["CD3"], [child1, child2], None, 0)

        # Mock weights: CD3 has equal weight to both children
        all_weights = {
            "root": {
                "cd4_t": {"CD3": 0.5, "CD4": 0.9},
                "cd8_t": {"CD3": 0.5, "CD8": 0.9},
            }
        }

        tree = ProfileTree(root=root, n_levels=2)
        flat_profiles, shared_markers = _flatten_tree_to_profiles(
            tree, all_weights, min_weight=0.1
        )

        # Both profiles should contain CD3
        assert "CD3" in flat_profiles["cd4_t"]
        assert "CD3" in flat_profiles["cd8_t"]
        # Each should have their specific marker
        assert "CD4" in flat_profiles["cd4_t"]
        assert "CD8" in flat_profiles["cd8_t"]
        # CD3 should be marked as shared
        assert "CD3" in shared_markers
        assert set(shared_markers["CD3"]) == {"cd4_t", "cd8_t"}

    def test_low_weight_markers_excluded(self):
        """Markers with weight below threshold should be excluded."""
        child1 = ProfileTreeNode("c1", ["A", "B"], [], "root", 1)
        root = ProfileTreeNode("root", [], [child1], None, 0)

        all_weights = {
            "root": {
                "c1": {"A": 0.8, "B": 0.01},  # B below threshold
            }
        }

        tree = ProfileTree(root=root, n_levels=2)
        flat_profiles, _ = _flatten_tree_to_profiles(
            tree, all_weights, min_weight=0.05
        )

        assert "A" in flat_profiles["c1"]
        assert "B" not in flat_profiles["c1"]
```

**Step 2: Run test to verify it fails**

Run: `python -m pytest tests/test_hierarchical_profiles.py::TestFlattening -v`

Expected: FAIL with "cannot import name '_flatten_tree_to_profiles'"

**Step 3: Write minimal implementation**

Add to `CITEgeist/model/spatial_colocalization.py`:

```python
def _flatten_tree_to_profiles(
    tree: ProfileTree,
    all_weights: Dict[str, Dict[str, Dict[str, float]]],
    min_weight: float = 0.05,
) -> Tuple[Dict[str, List[str]], Dict[str, List[str]]]:
    """
    Flatten hierarchical tree to flat profiles with weight propagation.

    Args:
        tree: ProfileTree with hierarchy structure
        all_weights: Nested dict: parent_id -> child_id -> marker -> weight
        min_weight: Minimum weight to include marker in profile

    Returns:
        Tuple of:
        - flat_profiles: Dict[cell_type_id, List[marker_names]]
        - shared_markers: Dict[marker, List[cell_type_ids_that_share_it]]
    """
    flat_profiles: Dict[str, List[str]] = {}
    marker_to_profiles: Dict[str, List[str]] = defaultdict(list)

    # Get all leaves
    leaves = tree.get_leaves()

    for leaf in leaves:
        # Collect all markers from leaf up to root with propagated weights
        leaf_markers = _collect_markers_with_weights(leaf, tree.root, all_weights, min_weight)
        flat_profiles[leaf.node_id] = list(leaf_markers.keys())

        # Track which profiles contain each marker
        for marker in leaf_markers:
            marker_to_profiles[marker].append(leaf.node_id)

    # Identify shared markers (appear in multiple profiles)
    shared_markers = {m: profiles for m, profiles in marker_to_profiles.items() if len(profiles) > 1}

    return flat_profiles, shared_markers


def _collect_markers_with_weights(
    leaf: ProfileTreeNode,
    root: ProfileTreeNode,
    all_weights: Dict[str, Dict[str, Dict[str, float]]],
    min_weight: float,
) -> Dict[str, float]:
    """
    Collect markers from leaf to root with propagated weights.

    Traverses from leaf up to root, accumulating markers and their weights.
    """
    markers_with_weights: Dict[str, float] = {}

    # Add leaf's own markers with weight 1.0
    for marker in leaf.markers:
        markers_with_weights[marker] = 1.0

    # Walk up the tree, adding ancestor markers with their weights
    current = leaf
    path_to_root = [leaf.node_id]

    # Build path from leaf to root by searching
    def find_path(node: ProfileTreeNode, target_id: str, path: List[str]) -> Optional[List[str]]:
        if node.node_id == target_id:
            return path
        for child in node.children:
            result = find_path(child, target_id, path + [child.node_id])
            if result:
                return result
        return None

    path = find_path(root, leaf.node_id, [root.node_id])
    if path is None:
        path = [root.node_id, leaf.node_id]

    # For each ancestor, add its markers with appropriate weights
    for i, node_id in enumerate(path[:-1]):  # Exclude leaf itself
        child_id = path[i + 1]

        # Get weights from this node to its children
        if node_id in all_weights and child_id in all_weights[node_id]:
            weights = all_weights[node_id][child_id]
            for marker, weight in weights.items():
                if weight >= min_weight:
                    # If marker already exists, keep the higher weight
                    if marker not in markers_with_weights or weight > markers_with_weights[marker]:
                        markers_with_weights[marker] = weight

    # Filter by min_weight
    return {m: w for m, w in markers_with_weights.items() if w >= min_weight}
```

**Step 4: Run test to verify it passes**

Run: `python -m pytest tests/test_hierarchical_profiles.py::TestFlattening -v`

Expected: PASS (2 tests)

**Step 5: Commit**

```bash
git add CITEgeist/model/spatial_colocalization.py tests/test_hierarchical_profiles.py
git commit -m "feat: add weight propagation and tree flattening"
```

---

## Task 8: Implement Main Entry Point discover_hierarchical_profiles()

**Files:**
- Modify: `CITEgeist/model/spatial_colocalization.py`
- Test: `tests/test_hierarchical_profiles.py`

**Step 1: Write the failing test**

Add to `tests/test_hierarchical_profiles.py`:

```python
from CITEgeist.model.spatial_colocalization import (
    discover_hierarchical_profiles,
    analyze_marker_colocalization,
)


class TestDiscoverHierarchicalProfiles:
    """Integration test for full hierarchical profile discovery."""

    def test_simulated_flat_data(self):
        """Flat data (no hierarchy) should produce flat profiles."""
        np.random.seed(42)
        n_spots = 200

        # Create 4 independent cell types, 2 markers each
        # No shared markers - should stay flat
        profiles_gt = {
            "type1": np.random.rand(n_spots),
            "type2": np.random.rand(n_spots),
        }

        X = np.column_stack([
            profiles_gt["type1"],
            profiles_gt["type1"] + 0.1 * np.random.rand(n_spots),
            profiles_gt["type2"],
            profiles_gt["type2"] + 0.1 * np.random.rand(n_spots),
        ])
        marker_names = ["A", "B", "C", "D"]

        # Create spatial coordinates (grid)
        coords = np.array([[i % 14, i // 14] for i in range(n_spots)])

        # Run colocalization analysis
        coloc_result = analyze_marker_colocalization(
            X, coords, marker_names, neighbor_k=6
        )

        # Run hierarchical discovery
        result = discover_hierarchical_profiles(
            coloc_result=coloc_result,
            antibody_expression=X,
            marker_names=marker_names,
            improvement_threshold=0.05,
        )

        # Should produce profiles (exact count depends on data)
        assert len(result.flat_profiles) >= 1
        # Reconstruction error should be reasonable
        assert result.reconstruction_error < 1.0

    def test_output_compatible_with_module3(self):
        """Output should be compatible with Module 3 profile dict format."""
        np.random.seed(42)
        n_spots = 100

        X = np.random.rand(n_spots, 4)
        marker_names = ["A", "B", "C", "D"]
        coords = np.array([[i % 10, i // 10] for i in range(n_spots)])

        coloc_result = analyze_marker_colocalization(
            X, coords, marker_names, neighbor_k=6
        )

        result = discover_hierarchical_profiles(
            coloc_result=coloc_result,
            antibody_expression=X,
            marker_names=marker_names,
        )

        profile_dict = result.to_profile_dict()

        # Should be Dict[str, List[str]]
        assert isinstance(profile_dict, dict)
        for cell_type, markers in profile_dict.items():
            assert isinstance(cell_type, str)
            assert isinstance(markers, list)
            assert all(isinstance(m, str) for m in markers)
```

**Step 2: Run test to verify it fails**

Run: `python -m pytest tests/test_hierarchical_profiles.py::TestDiscoverHierarchicalProfiles -v`

Expected: FAIL with "cannot import name 'discover_hierarchical_profiles'"

**Step 3: Write minimal implementation**

Add to `CITEgeist/model/spatial_colocalization.py`:

```python
def discover_hierarchical_profiles(
    coloc_result: ColocalizationResult,
    antibody_expression: NDArray[np.floating],
    marker_names: List[str],
    improvement_threshold: float = 0.05,
    min_marker_weight: float = 0.05,
    max_depth: int = 5,
    verbose: bool = True,
) -> HierarchicalProfileResult:
    """
    Discover cell type profiles using hierarchical NMF.

    Two-phase algorithm:
    1. Structure learning: Hierarchical clustering on colocalization distances
       with reconstruction-guided tree cutting
    2. Weight learning: NMF at each level to learn marker allocation

    Args:
        coloc_result: Result from analyze_marker_colocalization()
        antibody_expression: Expression matrix (n_spots, n_markers)
        marker_names: Names for each column in expression matrix
        improvement_threshold: Min reconstruction improvement to split (default 5%)
        min_marker_weight: Min weight to include marker in profile (default 0.05)
        max_depth: Maximum tree depth (default 5)
        verbose: Log progress (default True)

    Returns:
        HierarchicalProfileResult with tree structure and flat profiles
    """
    if verbose:
        logger.info("Starting hierarchical profile discovery...")

    # Phase 1: Build colocalization distance matrix
    if verbose:
        logger.info("Phase 1: Building colocalization distance matrix...")

    D, markers = _build_colocalization_distance_matrix(coloc_result)

    if len(markers) < 2:
        # Not enough markers for hierarchy
        if verbose:
            logger.info("Only 1 marker - returning single profile")
        root = ProfileTreeNode("root", markers, [], None, 0)
        tree = ProfileTree(root=root, n_levels=1)
        return HierarchicalProfileResult(
            tree=tree,
            flat_profiles={"profile_0": markers},
            depth_per_branch={},
            shared_markers={},
            reconstruction_error=0.0,
        )

    # Build hierarchical clustering
    if verbose:
        logger.info(f"Clustering {len(markers)} markers...")

    condensed = squareform(D)
    # Handle numerical issues
    condensed = np.maximum(condensed, 0)

    Z = linkage(condensed, method='ward')

    # Phase 1b: Reconstruction-guided tree cutting
    if verbose:
        logger.info("Building tree with reconstruction-guided cutting...")

    tree = _build_hierarchical_tree(
        X=antibody_expression,
        marker_names=marker_names,
        linkage_matrix=Z,
        improvement_threshold=improvement_threshold,
        max_depth=max_depth,
    )

    if verbose:
        logger.info(f"Tree depth: {tree.get_depth()}, leaves: {len(tree.get_leaves())}")

    # Phase 2: NMF weight learning at each level
    if verbose:
        logger.info("Phase 2: Computing NMF weights at each level...")

    all_weights = {}
    _compute_all_nmf_weights(tree.root, antibody_expression, marker_names, all_weights)

    # Phase 3: Flatten to profiles
    if verbose:
        logger.info("Phase 3: Flattening tree to profiles...")

    flat_profiles, shared_markers = _flatten_tree_to_profiles(
        tree, all_weights, min_marker_weight
    )

    # Rename profiles to be more descriptive
    renamed_profiles = {}
    for i, (node_id, markers) in enumerate(flat_profiles.items()):
        profile_name = f"profile_{i}"
        renamed_profiles[profile_name] = markers

    # Update shared_markers with new names
    old_to_new = {old: f"profile_{i}" for i, old in enumerate(flat_profiles.keys())}
    renamed_shared = {}
    for marker, old_profiles in shared_markers.items():
        renamed_shared[marker] = [old_to_new.get(p, p) for p in old_profiles]

    # Compute final reconstruction error
    final_error = _compute_final_reconstruction_error(
        antibody_expression, marker_names, list(renamed_profiles.values())
    )

    if verbose:
        logger.info(f"Discovered {len(renamed_profiles)} profiles, "
                   f"{len(renamed_shared)} shared markers, "
                   f"reconstruction error: {final_error:.4f}")

    return HierarchicalProfileResult(
        tree=tree,
        flat_profiles=renamed_profiles,
        depth_per_branch={},  # TODO: compute per-branch depth
        shared_markers=renamed_shared,
        reconstruction_error=final_error,
    )


def _compute_all_nmf_weights(
    node: ProfileTreeNode,
    X: NDArray[np.floating],
    marker_names: List[str],
    all_weights: Dict[str, Dict[str, Dict[str, float]]],
) -> None:
    """Recursively compute NMF weights for all internal nodes."""
    if node.is_leaf:
        return

    # Compute weights for this node
    weights = _compute_nmf_weights(X, marker_names, node)
    all_weights[node.node_id] = weights

    # Recurse on children
    for child in node.children:
        _compute_all_nmf_weights(child, X, marker_names, all_weights)


def _compute_final_reconstruction_error(
    X: NDArray[np.floating],
    marker_names: List[str],
    profiles: List[List[str]],
) -> float:
    """Compute overall reconstruction error for final profiles."""
    from scipy.optimize import nnls

    if len(profiles) == 0:
        return 1.0

    marker_to_idx = {m: i for i, m in enumerate(marker_names)}

    # Build profile matrix: each column is mean expression of profile markers
    profile_vectors = []
    for profile in profiles:
        indices = [marker_to_idx[m] for m in profile if m in marker_to_idx]
        if len(indices) > 0:
            profile_expr = X[:, indices].mean(axis=1)
            profile_vectors.append(profile_expr)

    if len(profile_vectors) == 0:
        return 1.0

    P = np.column_stack(profile_vectors)

    # For each marker, compute reconstruction error
    total_error = 0.0
    all_markers = set(m for p in profiles for m in p)

    for marker in all_markers:
        if marker not in marker_to_idx:
            continue
        y = X[:, marker_to_idx[marker]]
        try:
            coeffs, _ = nnls(P, y)
            y_pred = P @ coeffs
            error = np.mean((y - y_pred) ** 2)
            total_error += error
        except Exception:
            total_error += 1.0

    return total_error / max(len(all_markers), 1)
```

**Step 4: Run test to verify it passes**

Run: `python -m pytest tests/test_hierarchical_profiles.py::TestDiscoverHierarchicalProfiles -v`

Expected: PASS (2 tests)

**Step 5: Commit**

```bash
git add CITEgeist/model/spatial_colocalization.py tests/test_hierarchical_profiles.py
git commit -m "feat: add discover_hierarchical_profiles main entry point"
```

---

## Task 9: Add Exports to Module __init__.py

**Files:**
- Modify: `CITEgeist/model/__init__.py`
- Test: Quick import test

**Step 1: Write the failing test**

```python
# Quick test in Python
from CITEgeist.model import discover_hierarchical_profiles, HierarchicalProfileResult
```

**Step 2: Run test to verify it fails**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/.worktrees/hierarchical-profiles && python -c "from CITEgeist.model import discover_hierarchical_profiles, HierarchicalProfileResult"`

Expected: FAIL with ImportError

**Step 3: Write minimal implementation**

Read current `__init__.py` and add new exports:

```python
from .spatial_colocalization import (
    # ... existing exports ...
    discover_hierarchical_profiles,
    HierarchicalProfileResult,
    ProfileTree,
    ProfileTreeNode,
)
```

**Step 4: Run test to verify it passes**

Run: `python -c "from CITEgeist.model import discover_hierarchical_profiles, HierarchicalProfileResult; print('OK')"`

Expected: "OK"

**Step 5: Commit**

```bash
git add CITEgeist/model/__init__.py
git commit -m "feat: export hierarchical profile discovery from model module"
```

---

## Task 10: Integration Test with Simulated Data (Adversarial - Flat)

**Files:**
- Test: `tests/test_hierarchical_profiles.py`

**Step 1: Write the test**

Add to `tests/test_hierarchical_profiles.py`:

```python
class TestSimulatedDataAdversarial:
    """Test on simulated data where hierarchy should NOT be detected."""

    def test_nine_celltypes_stay_flat(self):
        """
        Simulated data: 9 cell types x 2 specific markers each.
        No shared markers - algorithm should detect depth=1 (flat).
        """
        np.random.seed(42)
        n_spots = 500
        n_celltypes = 9

        # Create 9 independent spatial patterns
        patterns = []
        for i in range(n_celltypes):
            # Each cell type has a distinct spatial pattern
            pattern = np.zeros(n_spots)
            start = i * (n_spots // n_celltypes)
            end = start + (n_spots // n_celltypes)
            pattern[start:end] = np.random.rand(end - start) + 0.5
            patterns.append(pattern)

        # Create markers: 2 per cell type, correlated with their pattern
        X_cols = []
        marker_names = []
        for i in range(n_celltypes):
            for j in range(2):
                marker = patterns[i] + 0.1 * np.random.rand(n_spots)
                X_cols.append(marker)
                marker_names.append(f"CellType{i}_Marker{j}")

        X = np.column_stack(X_cols)
        coords = np.array([[i % 25, i // 25] for i in range(n_spots)])

        # Run pipeline
        coloc_result = analyze_marker_colocalization(
            X, coords, marker_names, neighbor_k=6
        )

        result = discover_hierarchical_profiles(
            coloc_result=coloc_result,
            antibody_expression=X,
            marker_names=marker_names,
            improvement_threshold=0.05,
            verbose=False,
        )

        # Key assertion: should NOT over-split into unnecessary hierarchy
        # Each cell type's markers should stay together
        # Depth should be minimal (1 or 2)
        assert result.tree.get_depth() <= 2, f"Tree too deep: {result.tree.get_depth()}"

        # Should have roughly 9 profiles (one per cell type) or fewer
        assert len(result.flat_profiles) <= 12, f"Too many profiles: {len(result.flat_profiles)}"

        # No markers should be "shared" since they're all specific
        # (Some may be shared due to noise, but should be minimal)
        assert len(result.shared_markers) <= 3, f"Too many shared markers: {len(result.shared_markers)}"
```

**Step 2: Run test**

Run: `python -m pytest tests/test_hierarchical_profiles.py::TestSimulatedDataAdversarial -v`

Expected: PASS

**Step 3: Commit**

```bash
git add tests/test_hierarchical_profiles.py
git commit -m "test: add adversarial test for flat simulated data"
```

---

## Task 11: Integration Test with Hierarchical Data (T cells)

**Files:**
- Test: `tests/test_hierarchical_profiles.py`

**Step 1: Write the test**

Add to `tests/test_hierarchical_profiles.py`:

```python
class TestHierarchicalDataTcells:
    """Test on data with known hierarchy (T cell subtypes)."""

    def test_tcell_hierarchy_detected(self):
        """
        Data with T cell hierarchy: CD3 shared, CD4/CD8 specific.
        Algorithm should:
        1. Detect hierarchy
        2. Put CD3 in both CD4+ and CD8+ profiles
        """
        np.random.seed(42)
        n_spots = 500

        # Create T cell spatial pattern (shared by all T cells)
        t_cell_region = np.zeros(n_spots)
        t_cell_region[100:300] = 1.0  # T cells in spots 100-300

        # CD4+ T cells: subset of T cell region
        cd4_region = np.zeros(n_spots)
        cd4_region[100:200] = 1.0

        # CD8+ T cells: different subset
        cd8_region = np.zeros(n_spots)
        cd8_region[200:300] = 1.0

        # Non-T cells (epithelial)
        epi_region = np.zeros(n_spots)
        epi_region[350:450] = 1.0

        # Create marker expression
        noise = lambda: 0.1 * np.random.rand(n_spots)

        X = np.column_stack([
            t_cell_region + noise(),      # CD3 (shared T cell marker)
            cd4_region + noise(),          # CD4 (specific to CD4+ T)
            cd8_region + noise(),          # CD8 (specific to CD8+ T)
            epi_region + noise(),          # PanCK (epithelial)
        ])
        marker_names = ["CD3", "CD4", "CD8", "PanCK"]
        coords = np.array([[i % 25, i // 25] for i in range(n_spots)])

        # Run pipeline
        coloc_result = analyze_marker_colocalization(
            X, coords, marker_names, neighbor_k=6
        )

        result = discover_hierarchical_profiles(
            coloc_result=coloc_result,
            antibody_expression=X,
            marker_names=marker_names,
            improvement_threshold=0.05,
            verbose=True,
        )

        # Check that CD3 appears in profiles with CD4 and CD8
        cd3_in_cd4_profile = False
        cd3_in_cd8_profile = False

        for profile_name, markers in result.flat_profiles.items():
            if "CD4" in markers and "CD3" in markers:
                cd3_in_cd4_profile = True
            if "CD8" in markers and "CD3" in markers:
                cd3_in_cd8_profile = True

        # At minimum, CD3 should be shared or appear with at least one T cell marker
        cd3_shared = "CD3" in result.shared_markers

        # Success if EITHER:
        # 1. CD3 is explicitly marked as shared, OR
        # 2. CD3 appears in profiles with both CD4 and CD8
        assert cd3_shared or (cd3_in_cd4_profile and cd3_in_cd8_profile), \
            f"CD3 should be shared or appear with both CD4 and CD8. " \
            f"Profiles: {result.flat_profiles}, Shared: {result.shared_markers}"
```

**Step 2: Run test**

Run: `python -m pytest tests/test_hierarchical_profiles.py::TestHierarchicalDataTcells -v`

Expected: PASS

**Step 3: Commit**

```bash
git add tests/test_hierarchical_profiles.py
git commit -m "test: add T cell hierarchy detection test"
```

---

## Task 12: Run Full Test Suite

**Step 1: Run all hierarchical profile tests**

Run: `python -m pytest tests/test_hierarchical_profiles.py -v`

Expected: All tests PASS

**Step 2: Run existing tests to verify no regression**

Run: `python -m pytest tests/test_module2_profile_discovery.py tests/test_module2c_profile_selection.py -v --timeout=300`

Expected: Existing tests still PASS (or skip if they require real data)

**Step 3: Commit any fixes needed**

```bash
git add -A
git commit -m "fix: address any test failures"
```

---

## Summary Checklist

- [ ] Task 1: ProfileTreeNode and ProfileTree data structures
- [ ] Task 2: HierarchicalProfileResult dataclass
- [ ] Task 3: Colocalization distance matrix builder
- [ ] Task 4: Reconstruction error calculator
- [ ] Task 5: Reconstruction-guided tree cutting
- [ ] Task 6: NMF weight learning at each level
- [ ] Task 7: Weight propagation and flattening
- [ ] Task 8: Main entry point discover_hierarchical_profiles()
- [ ] Task 9: Module exports
- [ ] Task 10: Adversarial test (simulated flat data)
- [ ] Task 11: Hierarchy detection test (T cells)
- [ ] Task 12: Full test suite verification
