"""Tests for gex_qc.py."""

import numpy as np
import pandas as pd
import pytest

from CITEgeist.model.qc import QCResult
from CITEgeist.model.qc.gex_qc import (
    analyze_per_gene,
    compare_nmf_programs,
    compute_gex_correlations,
    run_gex_qc,
)


@pytest.mark.unit
class TestGEXCorrelations:
    """Tests for per-type GEX correlation computation."""

    def test_perfect_layers(self, mock_gt_gex_layers):
        """Identical predicted and GT layers should yield r~1."""
        result = compute_gex_correlations(mock_gt_gex_layers, mock_gt_gex_layers)
        for ct in mock_gt_gex_layers:
            assert result[ct]["pearson_r"] > 0.99

    def test_returns_expected_keys(self, mock_gt_gex_layers):
        rng = np.random.default_rng(42)
        pred = {ct: df + rng.random(df.shape) * 0.1
                for ct, df in mock_gt_gex_layers.items()}
        result = compute_gex_correlations(pred, mock_gt_gex_layers)
        for ct in mock_gt_gex_layers:
            for key in ["pearson_r", "spearman_rho", "pearson_spearman_gap"]:
                assert key in result[ct]


@pytest.mark.unit
class TestPerGeneAnalysis:
    """Tests for per-gene scaling analysis."""

    def test_returns_dataframe(self, mock_gt_gex_layers):
        ct = list(mock_gt_gex_layers.keys())[0]
        pred_layer = mock_gt_gex_layers[ct] * 1.5  # Scaled
        result = analyze_per_gene(pred_layer, mock_gt_gex_layers[ct])
        assert isinstance(result, pd.DataFrame)
        assert "spatial_r" in result.columns
        assert "log_mean_ratio" in result.columns
        assert "classification" in result.columns

    def test_scaled_data_classified_correctly(self, mock_gt_gex_layers):
        """Pure scaling should classify genes as magnitude_only."""
        ct = list(mock_gt_gex_layers.keys())[0]
        pred_layer = mock_gt_gex_layers[ct] * 2.0
        result = analyze_per_gene(pred_layer, mock_gt_gex_layers[ct])
        magnitude_only = (result["classification"] == "magnitude_only").sum()
        assert magnitude_only > len(result) * 0.5  # Most should be magnitude-only


@pytest.mark.unit
class TestNMFPropagation:
    """Tests for NMF program recovery diagnostic."""

    def test_returns_cosine_matrix(self, mock_gt_gex_layers):
        ct = list(mock_gt_gex_layers.keys())[0]
        result = compare_nmf_programs(
            mock_gt_gex_layers[ct], mock_gt_gex_layers[ct], k=3
        )
        assert result.shape == (3, 3)
        # Self-comparison should have high diagonal
        assert np.mean(np.diag(result)) > 0.5


@pytest.mark.unit
class TestRunGEXQC:
    """Tests for the main entry point."""

    def test_benchmark_mode(self, mock_gt_gex_layers):
        rng = np.random.default_rng(42)
        pred = {ct: df + rng.random(df.shape) * 0.5
                for ct, df in mock_gt_gex_layers.items()}
        result = run_gex_qc(pred_gex_layers=pred, gt_gex_layers=mock_gt_gex_layers)
        assert isinstance(result, QCResult)
        assert "correlations" in result.metrics
        assert "bar_pearson_spearman" in result.figures

    def test_self_consistency_mode(self, mock_gt_gex_layers):
        """Self-consistency mode with reference should return cosine similarities."""
        result = run_gex_qc(
            pred_gex_layers=mock_gt_gex_layers,
            gt_gex_layers=None,
            reference_profiles=mock_gt_gex_layers,
        )
        assert isinstance(result, QCResult)
        assert "reference_cosine" in result.metrics
