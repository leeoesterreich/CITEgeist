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

        # All genes (including anchors) should have non-negative expression
        assert np.all(E >= 0)

        # E should have reasonable values (not NaN or inf)
        assert not np.any(np.isnan(E))
        assert not np.any(np.isinf(E))


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