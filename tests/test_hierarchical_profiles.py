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
