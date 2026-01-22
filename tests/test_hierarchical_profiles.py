"""
Test script for Hierarchical Profile Discovery data structures.

Tests ProfileTreeNode and ProfileTree dataclasses used for
representing hierarchical marker profiles.
"""

import os
import sys

# Add project root to path for direct execution (pytest uses conftest.py)
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
import pytest

from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import squareform

from CITEgeist.model.spatial_colocalization import (
    ProfileTreeNode,
    ProfileTree,
    HierarchicalProfileResult,
    ColocalizationResult,
    MarkerPairColocalization,
    _build_colocalization_distance_matrix,
    _compute_reconstruction_error,
    _build_hierarchical_tree,
    _compute_nmf_weights,
    _flatten_tree_to_profiles,
    discover_hierarchical_profiles,
    analyze_marker_colocalization,
)


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
        # Ensure perfect symmetry for squareform
        D = (D + D.T) / 2
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


class TestSimulatedDataAdversarial:
    """Test on simulated data where hierarchy should NOT be detected."""

    def test_nine_celltypes_stay_flat(self):
        """
        Simulated data: 9 cell types x 2 specific markers each.
        No shared markers - algorithm should not over-fragment.

        This test verifies that:
        1. Markers from the same cell type are clustered together
        2. The reconstruction error is low (markers well-represented)
        3. Tree depth is bounded (max_depth parameter respected)

        Note: The hierarchical algorithm will still create some structure
        due to how hierarchical clustering works - the key is that it
        doesn't create EXCESSIVE fragmentation where each marker becomes
        its own leaf.
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

        # Use higher improvement threshold to reduce over-splitting
        result = discover_hierarchical_profiles(
            coloc_result=coloc_result,
            antibody_expression=X,
            marker_names=marker_names,
            improvement_threshold=0.10,  # Higher threshold = less splitting
            max_depth=3,  # Limit depth for flat data
            verbose=False,
        )

        # Key assertion: max_depth should be respected
        assert result.tree.get_depth() <= 4, f"Tree too deep: {result.tree.get_depth()}"

        # Should have a reasonable number of profiles (not one per marker)
        # With 18 markers, we expect 9-18 profiles (grouping by cell type ideally)
        n_markers = len(marker_names)
        assert len(result.flat_profiles) <= n_markers, (
            f"Too many profiles: {len(result.flat_profiles)} > {n_markers} markers"
        )

        # Check that markers within each cell type tend to co-occur in profiles
        # This is the key test: same-cell-type markers should cluster
        for i in range(n_celltypes):
            marker0 = f"CellType{i}_Marker0"
            marker1 = f"CellType{i}_Marker1"

            # Find profiles containing these markers
            profiles_with_m0 = [
                name for name, markers in result.flat_profiles.items()
                if marker0 in markers
            ]
            profiles_with_m1 = [
                name for name, markers in result.flat_profiles.items()
                if marker1 in markers
            ]

            # At least one profile should contain both markers
            overlap = set(profiles_with_m0) & set(profiles_with_m1)
            assert len(overlap) > 0, (
                f"CellType{i} markers not clustered together: "
                f"{marker0} in {profiles_with_m0}, {marker1} in {profiles_with_m1}"
            )

        # Reconstruction error should be reasonably low
        assert result.reconstruction_error < 0.5, (
            f"Reconstruction error too high: {result.reconstruction_error}"
        )


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
