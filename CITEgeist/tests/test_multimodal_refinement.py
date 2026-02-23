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
