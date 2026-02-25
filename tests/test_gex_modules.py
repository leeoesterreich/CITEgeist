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

    def test_filters_low_expressing_genes(self):
        """Genes not expressed in enough spots should be filtered out."""
        np.random.seed(42)
        n_spots, n_genes, n_types = 100, 20, 2

        proportions = np.random.dirichlet([1, 1], size=n_spots)
        gene_expression = np.zeros((n_spots, n_genes))

        # Gene 0: expressed in only 5% of spots (below default 10% threshold)
        expressing_spots = np.random.choice(n_spots, size=5, replace=False)
        gene_expression[expressing_spots, 0] = proportions[expressing_spots, 0] * 100

        # Genes 1-5: expressed in 50% of spots, correlated with type 0
        for g in range(1, 6):
            expressing_spots = np.random.choice(n_spots, size=50, replace=False)
            gene_expression[expressing_spots, g] = proportions[expressing_spots, 0] * 100 + np.random.randn(50) * 5

        from CITEgeist.model.gex_modules import discover_anchor_genes

        anchors, _ = discover_anchor_genes(
            gene_expression=gene_expression,
            cell_proportions=proportions,
            min_anchors=3,
            max_anchors=10,
            min_expressing_spots=0.1,
        )

        # Gene 0 should NOT be an anchor (too sparse)
        if 0 in anchors:
            assert 0 not in anchors[0]

    def test_handles_zero_variance_genes(self):
        """Genes with zero variance should not cause errors."""
        np.random.seed(42)
        n_spots, n_genes, n_types = 100, 10, 2

        proportions = np.random.dirichlet([1, 1], size=n_spots)
        gene_expression = np.random.rand(n_spots, n_genes) * 10

        # Gene 0: constant expression (zero variance)
        gene_expression[:, 0] = 5.0

        # Genes 1-5: correlated with type 0
        for g in range(1, 6):
            gene_expression[:, g] = proportions[:, 0] * 100 + np.random.randn(n_spots) * 5

        from CITEgeist.model.gex_modules import discover_anchor_genes

        # Should not raise
        anchors, thresholds = discover_anchor_genes(
            gene_expression=gene_expression,
            cell_proportions=proportions,
            min_anchors=3,
            max_anchors=10,
        )

        # Gene 0 should NOT be an anchor
        if 0 in anchors:
            assert 0 not in anchors[0]

    def test_returns_empty_for_no_correlations(self):
        """Returns empty anchors when no genes correlate significantly."""
        np.random.seed(42)
        n_spots, n_genes, n_types = 100, 10, 2

        proportions = np.random.dirichlet([1, 1], size=n_spots)
        # Pure noise - no correlation with proportions
        gene_expression = np.random.rand(n_spots, n_genes) * 10

        from CITEgeist.model.gex_modules import discover_anchor_genes

        anchors, thresholds = discover_anchor_genes(
            gene_expression=gene_expression,
            cell_proportions=proportions,
            min_anchors=3,
            max_anchors=10,
        )

        # Should still return dict structure, may have few/no anchors
        assert isinstance(anchors, dict)
        assert isinstance(thresholds, dict)
        assert len(anchors) == n_types


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
