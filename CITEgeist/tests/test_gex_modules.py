# CITEgeist/tests/test_gex_modules.py
"""Tests for GEX module-aware enrichment functions."""

import numpy as np
import pytest
from CITEgeist.model.gex_modules import discover_anchor_genes


def test_discover_anchor_genes_returns_weights():
    """Test that discover_anchor_genes returns anchor weights (correlations)."""
    np.random.seed(42)
    N, M, T = 100, 50, 3

    # Create synthetic data where gene 0 correlates with type 0, etc.
    proportions = np.random.dirichlet(np.ones(T), N)
    gene_expr = np.zeros((N, M))

    # First 3 genes are markers for each cell type
    for t in range(T):
        gene_expr[:, t] = proportions[:, t] * 100 + np.random.randn(N) * 5

    # Rest are noise
    gene_expr[:, T:] = np.random.rand(N, M - T) * 10

    anchors, thresholds, weights = discover_anchor_genes(
        gene_expr, proportions, min_anchors=1, max_anchors=5
    )

    # Should return 3 values now
    assert len(weights) == T, "Should have weights for each cell type"

    # Weights should be dicts mapping gene_idx -> correlation
    for t in range(T):
        assert isinstance(weights[t], dict), f"weights[{t}] should be dict"
        # Marker gene t should be in anchors for type t with high weight
        if t in anchors[t]:
            assert weights[t][t] > 0.3, f"Marker gene {t} should have r > 0.3"


def test_discover_anchor_genes_basic():
    """Test basic functionality of discover_anchor_genes."""
    np.random.seed(123)
    N, M, T = 50, 20, 2

    # Create data with clear markers
    proportions = np.random.dirichlet(np.ones(T), N)
    gene_expr = np.zeros((N, M))

    # Gene 0 strongly correlates with type 0
    gene_expr[:, 0] = proportions[:, 0] * 50 + np.random.randn(N) * 2
    # Gene 1 strongly correlates with type 1
    gene_expr[:, 1] = proportions[:, 1] * 50 + np.random.randn(N) * 2
    # Rest are noise
    gene_expr[:, 2:] = np.random.rand(N, M - 2) * 5

    anchors, thresholds, weights = discover_anchor_genes(
        gene_expr, proportions, min_anchors=1, max_anchors=3
    )

    # Check structure
    assert len(anchors) == T
    assert len(thresholds) == T
    assert len(weights) == T

    # Gene 0 should be anchor for type 0
    assert 0 in anchors[0], "Gene 0 should be anchor for type 0"
    # Gene 1 should be anchor for type 1
    assert 1 in anchors[1], "Gene 1 should be anchor for type 1"

    # Weights should contain these genes with high correlation
    assert 0 in weights[0], "Gene 0 should have weight for type 0"
    assert 1 in weights[1], "Gene 1 should have weight for type 1"
    assert weights[0][0] > 0.5, "Gene 0 should have high correlation with type 0"
    assert weights[1][1] > 0.5, "Gene 1 should have high correlation with type 1"


def test_discover_anchor_genes_empty_cell_type():
    """Test handling of cell types with zero variance."""
    np.random.seed(456)
    N, M, T = 50, 10, 2

    # Type 1 has zero variance (absent cell type)
    proportions = np.zeros((N, T))
    proportions[:, 0] = np.random.rand(N)
    proportions[:, 1] = 0.0  # Zero variance
    proportions = proportions / proportions.sum(axis=1, keepdims=True)

    gene_expr = np.random.rand(N, M) * 10

    anchors, thresholds, weights = discover_anchor_genes(
        gene_expr, proportions, min_anchors=1, max_anchors=5
    )

    # Type 1 should have empty anchors and weights
    assert anchors[1] == [], "Type 1 should have no anchors"
    assert weights[1] == {}, "Type 1 should have no weights"


def test_compute_proportion_enrichment_no_smoothing():
    """Proportion enrichment should NOT have 80/20 smoothing."""
    from CITEgeist.model.gurobi_impl import compute_proportion_enrichment

    np.random.seed(42)
    N, T = 20, 4

    # Create data where gene is only expressed in spots with high type-0 proportion
    cell_props = np.random.dirichlet(np.ones(T), N)
    gene_expr = np.zeros(N)
    high_type0_spots = cell_props[:, 0] > 0.5
    gene_expr[high_type0_spots] = 100

    enrichment = compute_proportion_enrichment(gene_expr, cell_props)

    # Without smoothing, type 0 should dominate
    assert enrichment[0] > 0.5, "Type 0 should have majority enrichment"

    # Should sum to 1
    assert abs(enrichment.sum() - 1.0) < 1e-6, "Enrichment should sum to 1"

    # Should NOT be smoothed toward uniform (1/T = 0.25)
    # Check that non-enriched types can go below 0.05 (smoothing floor was ~0.05)
    assert enrichment[1:].min() < 0.1, "Non-enriched types should be low (no smoothing)"


def test_compute_marker_enrichment():
    """Marker enrichment uses correlation with anchor genes."""
    from CITEgeist.model.gurobi_impl import compute_marker_enrichment

    np.random.seed(42)
    N, T = 20, 3

    # Gene correlates with type-0 anchors
    anchor_expr = np.random.rand(N, T) * 10  # mean anchor expr per type
    gene_expr = anchor_expr[:, 0] + np.random.randn(N) * 0.1  # correlates with type 0

    anchor_weights = np.array([0.7, 0.5, 0.3])  # type 0 has strongest anchor

    enrichment = compute_marker_enrichment(gene_expr, anchor_expr, anchor_weights)

    # Type 0 should have highest enrichment (correlates + highest anchor weight)
    assert enrichment[0] > enrichment[1], "Type 0 should have highest enrichment"
    assert enrichment[0] > enrichment[2], "Type 0 should have highest enrichment"
    assert abs(enrichment.sum() - 1.0) < 1e-6, "Should sum to 1"


def test_compute_adaptive_enrichment():
    """Adaptive enrichment blends proportion and marker based on variance."""
    from CITEgeist.model.gurobi_impl import compute_adaptive_enrichment

    np.random.seed(42)
    N, T = 20, 4

    # Case 1: High variance proportion enrichment -> trust proportions
    prop_enrich_peaked = np.array([0.7, 0.1, 0.1, 0.1])
    marker_enrich = np.array([0.25, 0.25, 0.25, 0.25])

    result = compute_adaptive_enrichment(prop_enrich_peaked, marker_enrich)

    # Should be close to prop_enrich (high variance = low anchor weight)
    assert result[0] > 0.5, "High-variance case should trust proportions"

    # Case 2: Low variance proportion enrichment -> use markers
    prop_enrich_flat = np.array([0.26, 0.25, 0.25, 0.24])
    marker_enrich_peaked = np.array([0.6, 0.2, 0.1, 0.1])

    result2 = compute_adaptive_enrichment(prop_enrich_flat, marker_enrich_peaked)

    # Should lean toward marker_enrich (low variance = high anchor weight)
    assert result2[0] > 0.4, "Low-variance case should use markers"


def test_precompute_anchor_expression():
    """Test that anchor expression matrix is computed correctly."""
    from CITEgeist.model.gurobi_impl import precompute_anchor_expression

    np.random.seed(42)
    N, M, T = 50, 100, 3

    gene_expr = np.random.rand(N, M) * 100

    # Anchors: type 0 -> genes [0,1], type 1 -> genes [2,3], type 2 -> genes [4,5]
    anchor_genes = {0: [0, 1], 1: [2, 3], 2: [4, 5]}
    anchor_weights = {
        0: {0: 0.8, 1: 0.6},  # gene 0 has r=0.8, gene 1 has r=0.6
        1: {2: 0.7, 3: 0.5},
        2: {4: 0.9, 5: 0.4},
    }

    anchor_expr, type_weights = precompute_anchor_expression(
        gene_expr, anchor_genes, anchor_weights
    )

    # Shape: (N, T) - mean anchor expression per type
    assert anchor_expr.shape == (N, T)

    # Type weights: mean correlation per type
    assert len(type_weights) == T
    assert type_weights[0] == pytest.approx(0.7, abs=0.01)  # mean(0.8, 0.6)
    assert type_weights[2] == pytest.approx(0.65, abs=0.01)  # mean(0.9, 0.4)
