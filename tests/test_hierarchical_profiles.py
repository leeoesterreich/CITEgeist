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

from CITEgeist.model.spatial_colocalization import (
    ProfileTreeNode,
    ProfileTree,
    HierarchicalProfileResult,
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


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
