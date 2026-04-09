"""Tests for proportion_qc.py."""

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
from anndata import AnnData

from CITEgeist.model.qc import QCResult
from CITEgeist.model.qc.proportion_qc import (
    compute_discrete_round_trip,
    compute_error_confusion,
    compute_proportion_correlations,
    decompose_scaling_error,
    run_proportion_qc,
)


@pytest.mark.unit
class TestCorrelationMetrics:
    """Tests for per-type correlation computation."""

    def test_perfect_correlation(self):
        """Identical predictions should yield r=1, rho=1."""
        gt = pd.DataFrame({"A": [0.5, 0.3, 0.2], "B": [0.5, 0.7, 0.8]})
        result = compute_proportion_correlations(gt, gt)
        for ct in ["A", "B"]:
            assert abs(result[ct]["pearson_r"] - 1.0) < 1e-6
            assert abs(result[ct]["spearman_rho"] - 1.0) < 1e-6

    def test_returns_all_expected_metrics(self, mock_proportions, mock_gt_proportions):
        result = compute_proportion_correlations(mock_proportions, mock_gt_proportions)
        for ct in mock_proportions.columns:
            for key in ["pearson_r", "spearman_rho", "rmse", "mae"]:
                assert key in result[ct], f"Missing {key} for {ct}"
        assert "overall" in result


@pytest.mark.unit
class TestErrorDecomposition:
    """Tests for scaling vs misallocation decomposition."""

    def test_pure_scaling_has_low_residual_fraction(self):
        """If predictions are just scaled GT, residual fraction should be small."""
        gt = np.array([0.1, 0.2, 0.3, 0.4, 0.5])
        pred = gt * 2.0  # Pure scaling
        result = decompose_scaling_error(pred, gt)
        assert result["residual_fraction"] < 0.1

    def test_random_noise_has_high_residual_fraction(self):
        """Random predictions should have high residual fraction."""
        rng = np.random.default_rng(42)
        gt = rng.random(50)
        pred = rng.random(50)
        result = decompose_scaling_error(pred, gt)
        assert result["residual_fraction"] > 0.3


@pytest.mark.unit
class TestDiscreteRoundTrip:
    """Tests for discrete assignment round-trip."""

    def test_round_trip_from_assignments(self, mock_sc_adata, mock_proportions):
        result = compute_discrete_round_trip(mock_sc_adata, mock_proportions)
        assert "discrete_proportions" in result
        assert isinstance(result["discrete_proportions"], pd.DataFrame)

    def test_zero_nuclei_spots_skipped(self):
        """Spots with no cells should be excluded."""
        adata = AnnData(
            X=sp.csr_matrix(np.ones((3, 10))),
            obs=pd.DataFrame({
                "cell_type": ["A", "A", "B"],
                "spot_id": ["s1", "s1", "s2"],  # s3 has no cells
            }, index=["c1", "c2", "c3"]),
        )
        props = pd.DataFrame(
            {"A": [0.5, 0.5, 0.5], "B": [0.5, 0.5, 0.5]},
            index=["s1", "s2", "s3"],
        )
        result = compute_discrete_round_trip(adata, props)
        assert "s3" not in result["discrete_proportions"].index


@pytest.mark.unit
class TestConfusionAnalysis:
    """Tests for over/under-prediction confusion matrix."""

    def test_returns_square_matrix(self, mock_proportions, mock_gt_proportions):
        mat = compute_error_confusion(mock_proportions, mock_gt_proportions)
        assert mat.shape[0] == mat.shape[1]
        assert list(mat.index) == list(mat.columns)


@pytest.mark.unit
class TestRunProportionQC:
    """Tests for the main entry point."""

    def test_returns_qcresult(self, mock_sc_adata, mock_proportions, mock_gt_proportions):
        result = run_proportion_qc(mock_sc_adata, mock_proportions, mock_gt_proportions)
        assert isinstance(result, QCResult)

    def test_produces_expected_figures(self, mock_sc_adata, mock_proportions, mock_gt_proportions):
        result = run_proportion_qc(mock_sc_adata, mock_proportions, mock_gt_proportions)
        assert "scatter_grid" in result.figures
        assert "spatial_errors" in result.figures
        assert "round_trip" in result.figures

    def test_flags_large_spearman_pearson_gap(self):
        """Types with large Spearman-Pearson gap should be flagged."""
        rng = np.random.default_rng(42)
        n = 100
        # Create data where outliers cause Pearson-Spearman divergence
        gt_vals = rng.random(n) * 0.5
        pred_vals = gt_vals.copy()
        pred_vals[:3] = 5.0  # Outliers

        gt = pd.DataFrame({"TypeA": gt_vals, "TypeB": 1 - gt_vals},
                          index=[f"s{i}" for i in range(n)])
        pred = pd.DataFrame({"TypeA": pred_vals, "TypeB": 1 - pred_vals},
                            index=[f"s{i}" for i in range(n)])

        adata = AnnData(
            X=sp.csr_matrix(np.ones((n, 10))),
            obs=pd.DataFrame({
                "cell_type": ["TypeA"] * n,
                "spot_id": [f"s{i}" for i in range(n)],
            }, index=[f"c{i}" for i in range(n)]),
        )
        result = run_proportion_qc(adata, pred, gt)
        gap_flags = [f for f in result.flags if "Spearman-Pearson gap" in f]
        assert len(gap_flags) > 0
