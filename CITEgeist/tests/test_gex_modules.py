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
