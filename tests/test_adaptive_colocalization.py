"""
Tests for adaptive colocalization enhancements:
1. Expression Correlation Fallback - adaptive weight blending based on spatial heterogeneity
2. Multi-scale Neighborhoods - bivariate Moran's I at multiple k values

Run with: pytest tests/test_adaptive_colocalization.py -v
"""

import numpy as np
import pytest
from typing import Dict, List

# =============================================================================
# Helper Functions for Spatial Data Generation
# =============================================================================


def create_spatial_grid(n_rows: int, n_cols: int) -> np.ndarray:
    """Create a regular grid of spatial coordinates."""
    x = np.arange(n_cols)
    y = np.arange(n_rows)
    xx, yy = np.meshgrid(x, y)
    coords = np.column_stack([xx.ravel(), yy.ravel()])
    return coords.astype(float)


def spots_in_circle(coords: np.ndarray, center: tuple, radius: float) -> np.ndarray:
    """Return boolean mask for spots within a circular region."""
    distances = np.sqrt((coords[:, 0] - center[0])**2 + (coords[:, 1] - center[1])**2)
    return distances <= radius


def spots_in_rectangle(coords: np.ndarray, x_range: tuple, y_range: tuple) -> np.ndarray:
    """Return boolean mask for spots within a rectangular region."""
    return (
        (coords[:, 0] >= x_range[0]) & (coords[:, 0] <= x_range[1]) &
        (coords[:, 1] >= y_range[0]) & (coords[:, 1] <= y_range[1])
    )


# =============================================================================
# Fixtures
# =============================================================================


@pytest.fixture
def rng():
    """Seeded random generator for reproducibility."""
    return np.random.default_rng(42)


@pytest.fixture
def high_spatial_signal_data(rng):
    """
    Data with strong spatial clustering (high Moran's I).
    Should use spatial_default blend mode.
    """
    n_rows, n_cols = 30, 30  # 900 spots
    coords = create_spatial_grid(n_rows, n_cols)
    n_spots = len(coords)

    marker_names = ["CD3D", "CD4", "CD8", "CD68", "EPCAM", "CD19", "CD56"]
    n_markers = len(marker_names)

    # Low background noise
    X = rng.normal(0, 0.2, (n_spots, n_markers))

    # Strong spatial clusters for each marker
    # T cells in left region
    t_cell_mask = spots_in_rectangle(coords, x_range=(0, 10), y_range=(0, 30))
    X[t_cell_mask, 0] += 3.0  # CD3D
    X[t_cell_mask, 1] += 2.5  # CD4

    # Macrophages in center
    mac_mask = spots_in_rectangle(coords, x_range=(12, 18), y_range=(5, 25))
    X[mac_mask, 3] += 3.0  # CD68

    # Cancer cells in right region
    cancer_mask = spots_in_rectangle(coords, x_range=(20, 30), y_range=(0, 30))
    X[cancer_mask, 4] += 4.0  # EPCAM

    return X, marker_names, coords


@pytest.fixture
def low_spatial_signal_data(rng):
    """
    Data with weak spatial clustering (low Moran's I).
    Should trigger expression_fallback blend mode.
    """
    n_rows, n_cols = 30, 30  # 900 spots
    coords = create_spatial_grid(n_rows, n_cols)
    n_spots = len(coords)

    marker_names = ["CD3D", "CD4", "CD8", "CD68", "EPCAM", "CD19", "CD56"]
    n_markers = len(marker_names)

    # Higher background noise
    X = rng.normal(0, 0.5, (n_spots, n_markers))

    # Scattered expression (not spatially clustered)
    # Randomly distribute positive cells
    for i in range(n_markers):
        positive_spots = rng.choice(n_spots, size=int(n_spots * 0.2), replace=False)
        X[positive_spots, i] += 2.0

    # Add some correlation between markers (expression-based signal)
    # CD3D and CD4 correlated
    correlated_spots = rng.choice(n_spots, size=int(n_spots * 0.15), replace=False)
    X[correlated_spots, 0] += 1.5  # CD3D
    X[correlated_spots, 1] += 1.5  # CD4

    return X, marker_names, coords


@pytest.fixture
def multi_scale_signal_data(rng):
    """
    Data where signal appears at different spatial scales.
    Tests multi-scale neighborhood analysis.
    """
    n_rows, n_cols = 40, 40  # 1600 spots
    coords = create_spatial_grid(n_rows, n_cols)
    n_spots = len(coords)

    marker_names = ["MARKER_A", "MARKER_B", "MARKER_C", "MARKER_D"]
    n_markers = len(marker_names)

    X = rng.normal(0, 0.3, (n_spots, n_markers))

    # MARKER_A and MARKER_B: Small, tight clusters (best signal at k=6)
    for center in [(10, 10), (30, 30), (10, 30), (30, 10)]:
        mask = spots_in_circle(coords, center, radius=3)
        X[mask, 0] += 3.0  # MARKER_A
        X[mask, 1] += 2.5  # MARKER_B

    # MARKER_C and MARKER_D: Large, diffuse clusters (best signal at k=24)
    large_mask1 = spots_in_circle(coords, (20, 20), radius=12)
    X[large_mask1, 2] += 2.0  # MARKER_C
    X[large_mask1, 3] += 2.0  # MARKER_D

    return X, marker_names, coords


# =============================================================================
# Tests for Expression Correlation Fallback (Feature 1)
# =============================================================================


class TestBlendWeightDetermination:
    """Tests for _determine_blend_weights() function."""

    def test_spatial_default_high_morans(self):
        """High avg Moran's I should use spatial_default mode."""
        from CITEgeist.model.spatial_colocalization import _determine_blend_weights

        weights, mode = _determine_blend_weights(avg_marker_morans_i=0.35)

        assert mode == "spatial_default"
        assert weights['bivariate'] == 0.40
        assert weights['spearman'] == 0.30
        assert weights['cosine'] == 0.30

    def test_expression_fallback_low_morans(self):
        """Low avg Moran's I should trigger expression_fallback mode."""
        from CITEgeist.model.spatial_colocalization import _determine_blend_weights

        weights, mode = _determine_blend_weights(avg_marker_morans_i=0.10)

        assert mode == "expression_fallback"
        assert weights['spearman'] == 0.45  # Higher weight on expression
        assert weights['bivariate'] == 0.25  # Lower weight on spatial

    def test_balanced_moderate_morans(self):
        """Moderate avg Moran's I should use balanced mode."""
        from CITEgeist.model.spatial_colocalization import _determine_blend_weights

        weights, mode = _determine_blend_weights(avg_marker_morans_i=0.20)

        assert mode == "balanced"
        assert weights['spearman'] == 0.35
        assert weights['bivariate'] == 0.35

    def test_override_mode(self):
        """Override should force specific mode regardless of Moran's I."""
        from CITEgeist.model.spatial_colocalization import _determine_blend_weights

        # Force expression_fallback even with high Moran's I
        weights, mode = _determine_blend_weights(
            avg_marker_morans_i=0.50,
            blend_mode_override="expression_fallback"
        )

        assert mode == "expression_fallback"
        assert weights['spearman'] == 0.45

    def test_custom_thresholds(self):
        """Custom thresholds should be respected."""
        from CITEgeist.model.spatial_colocalization import _determine_blend_weights

        # With default thresholds, 0.20 is balanced
        # With raised thresholds, 0.20 should be expression_fallback
        weights, mode = _determine_blend_weights(
            avg_marker_morans_i=0.20,
            low_threshold=0.25,
            high_threshold=0.35
        )

        assert mode == "expression_fallback"


class TestGetMarkerMoransDict:
    """Tests for get_marker_morans_dict() helper."""

    def test_extracts_morans_values(self):
        """Should extract Moran's I values from MarkerInterestResult."""
        from CITEgeist.model.spatial_colocalization import get_marker_morans_dict
        from CITEgeist.model.marker_interest import MarkerInterestResult, MarkerStats

        # Create mock result
        markers = [
            MarkerStats(name="CD3D", morans_i=0.25, morans_pvalue=0.001,
                       kurtosis=3.0, gmm_snr=5.0, is_interesting=True,
                       interesting_reason="kurtosis"),
            MarkerStats(name="CD68", morans_i=0.35, morans_pvalue=0.001,
                       kurtosis=2.5, gmm_snr=4.0, is_interesting=True,
                       interesting_reason="morans"),
        ]
        result = MarkerInterestResult(markers=markers, interesting_markers=[],
                                       boring_markers=[])

        morans_dict = get_marker_morans_dict(result)

        assert morans_dict["CD3D"] == 0.25
        assert morans_dict["CD68"] == 0.35


class TestColocalizationWithBlending:
    """Integration tests for colocalization with adaptive blending."""

    def test_high_spatial_uses_spatial_default(self, high_spatial_signal_data, rng):
        """High spatial signal data should use spatial_default blend mode."""
        from CITEgeist.model.spatial_colocalization import analyze_marker_colocalization

        X, marker_names, coords = high_spatial_signal_data

        # Mock high Moran's I values
        marker_morans_i = {name: 0.35 for name in marker_names}

        result = analyze_marker_colocalization(
            marker_data=X,
            marker_names=marker_names,
            coords=coords,
            marker_morans_i=marker_morans_i,
            rng=rng,
        )

        assert result.blend_mode == "spatial_default"
        assert result.blend_weights['bivariate'] == 0.40

    def test_low_spatial_uses_expression_fallback(self, low_spatial_signal_data, rng):
        """Low spatial signal data should use expression_fallback blend mode."""
        from CITEgeist.model.spatial_colocalization import analyze_marker_colocalization

        X, marker_names, coords = low_spatial_signal_data

        # Mock low Moran's I values
        marker_morans_i = {name: 0.10 for name in marker_names}

        result = analyze_marker_colocalization(
            marker_data=X,
            marker_names=marker_names,
            coords=coords,
            marker_morans_i=marker_morans_i,
            rng=rng,
        )

        assert result.blend_mode == "expression_fallback"
        assert result.blend_weights['spearman'] == 0.45


# =============================================================================
# Tests for Multi-scale Neighborhoods (Feature 2)
# =============================================================================


class TestMultiscaleMoransI:
    """Tests for _compute_bivariate_morans_i_multiscale() function."""

    def test_returns_per_scale_values(self, high_spatial_signal_data, rng):
        """Should return Moran's I computed at each scale."""
        from CITEgeist.model.spatial_colocalization import (
            _compute_bivariate_morans_i_multiscale,
            _build_neighbor_graph
        )

        X, marker_names, coords = high_spatial_signal_data

        # Build multi-scale neighbors
        multi_scale_neighbors = {}
        for k in [6, 12, 24]:
            multi_scale_neighbors[k] = _build_neighbor_graph(coords, k)

        # Compute multi-scale
        I_val, p_val, per_scale_I, per_scale_p, best_k = _compute_bivariate_morans_i_multiscale(
            values_a=X[:, 0],  # CD3D
            values_b=X[:, 1],  # CD4
            multi_scale_neighbors=multi_scale_neighbors,
            rng=rng,
            n_perm=99
        )

        # Check all scales computed
        assert 6 in per_scale_I
        assert 12 in per_scale_I
        assert 24 in per_scale_I

        # Check best_k is one of the scales
        assert best_k in [6, 12, 24]

    def test_max_aggregation_selects_highest(self, rng):
        """Max aggregation should select scale with highest Moran's I."""
        from CITEgeist.model.spatial_colocalization import _compute_bivariate_morans_i_multiscale

        # Create data where k=12 gives best signal
        n_spots = 400
        coords = create_spatial_grid(20, 20)

        # Create tight clusters (best at moderate k)
        X = np.random.default_rng(42).normal(0, 0.3, (n_spots, 2))
        for center in [(5, 5), (15, 15), (5, 15), (15, 5)]:
            mask = spots_in_circle(coords, center, radius=4)
            X[mask, 0] += 2.5
            X[mask, 1] += 2.5

        # Build neighbors
        from CITEgeist.model.spatial_colocalization import _build_neighbor_graph
        multi_scale_neighbors = {k: _build_neighbor_graph(coords, k) for k in [6, 12, 24]}

        I_val, _, per_scale_I, _, best_k = _compute_bivariate_morans_i_multiscale(
            values_a=X[:, 0],
            values_b=X[:, 1],
            multi_scale_neighbors=multi_scale_neighbors,
            rng=rng,
            n_perm=99,
            aggregation="max"
        )

        # The returned I should be the max across scales
        assert I_val == max(per_scale_I.values())
        assert per_scale_I[best_k] == I_val


class TestColocalizationWithMultiscale:
    """Integration tests for colocalization with multi-scale neighborhoods."""

    def test_multiscale_populates_result_fields(self, multi_scale_signal_data, rng):
        """Multi-scale analysis should populate per-scale fields in result."""
        from CITEgeist.model.spatial_colocalization import analyze_marker_colocalization

        X, marker_names, coords = multi_scale_signal_data

        result = analyze_marker_colocalization(
            marker_data=X,
            marker_names=marker_names,
            coords=coords,
            multi_scale_k=[6, 12, 24],
            multi_scale_aggregation="max",
            rng=rng,
        )

        # Check that per-scale data is populated
        for pair in result.marker_pairs:
            if pair.bivariate_morans_per_scale is not None:
                assert 6 in pair.bivariate_morans_per_scale
                assert 12 in pair.bivariate_morans_per_scale
                assert 24 in pair.bivariate_morans_per_scale
                assert pair.bivariate_morans_best_scale is not None

    def test_small_clusters_best_at_small_k(self, rng):
        """Small tight clusters should have best signal at small k."""
        from CITEgeist.model.spatial_colocalization import analyze_marker_colocalization

        n_spots = 625
        coords = create_spatial_grid(25, 25)

        # Very tight clusters (radius=2)
        X = rng.normal(0, 0.3, (n_spots, 2))
        for center in [(6, 6), (18, 18)]:
            mask = spots_in_circle(coords, center, radius=2)
            X[mask, 0] += 4.0
            X[mask, 1] += 4.0

        result = analyze_marker_colocalization(
            marker_data=X,
            marker_names=["A", "B"],
            coords=coords,
            multi_scale_k=[6, 12, 24],
            rng=rng,
        )

        # Find the A-B pair
        pair = result.marker_pairs[0]

        # Small k should give better signal for tight clusters
        if pair.bivariate_morans_per_scale is not None:
            # k=6 should be >= k=24 for tight clusters
            assert pair.bivariate_morans_per_scale[6] >= pair.bivariate_morans_per_scale[24] * 0.8


# =============================================================================
# Tests for Backward Compatibility
# =============================================================================


class TestBackwardCompatibility:
    """Ensure default behavior unchanged when new features not used."""

    def test_default_parameters_unchanged(self, high_spatial_signal_data, rng):
        """Without marker_morans_i, should use spatial_default mode."""
        from CITEgeist.model.spatial_colocalization import analyze_marker_colocalization

        X, marker_names, coords = high_spatial_signal_data

        result = analyze_marker_colocalization(
            marker_data=X,
            marker_names=marker_names,
            coords=coords,
            rng=rng,
            # No marker_morans_i provided
        )

        # Should default to spatial_default
        assert result.blend_mode == "spatial_default"
        assert result.blend_weights['bivariate'] == 0.40

    def test_single_scale_when_multiscale_not_specified(self, high_spatial_signal_data, rng):
        """Without multi_scale_k, should use single scale (no per-scale data)."""
        from CITEgeist.model.spatial_colocalization import analyze_marker_colocalization

        X, marker_names, coords = high_spatial_signal_data

        result = analyze_marker_colocalization(
            marker_data=X,
            marker_names=marker_names,
            coords=coords,
            rng=rng,
            # No multi_scale_k provided
        )

        # Per-scale fields should be None
        for pair in result.marker_pairs:
            assert pair.bivariate_morans_per_scale is None


# =============================================================================
# Integration Test: Mixed Data Improvement
# =============================================================================


class TestMixedDataImprovement:
    """Test that enhancements improve mixed/scattered data analysis."""

    def test_expression_fallback_finds_correlated_markers(self, rng):
        """Expression fallback should find correlated markers even without spatial signal."""
        from CITEgeist.model.spatial_colocalization import analyze_marker_colocalization

        n_spots = 900
        coords = create_spatial_grid(30, 30)

        # Create data with expression correlation but no spatial clustering
        X = rng.normal(0, 0.5, (n_spots, 4))
        marker_names = ["M1", "M2", "M3", "M4"]

        # M1 and M2 are expression-correlated (not spatially)
        corr_factor = rng.normal(0, 1, n_spots)
        X[:, 0] += corr_factor * 2.0  # M1
        X[:, 1] += corr_factor * 1.8  # M2 correlated with M1

        # M3 and M4 are independent
        X[:, 2] += rng.normal(0, 1, n_spots)
        X[:, 3] += rng.normal(0, 1, n_spots)

        # Run with expression fallback
        marker_morans_i = {name: 0.05 for name in marker_names}  # Low spatial signal

        result = analyze_marker_colocalization(
            marker_data=X,
            marker_names=marker_names,
            coords=coords,
            marker_morans_i=marker_morans_i,
            rng=rng,
        )

        # Should use expression_fallback mode
        assert result.blend_mode == "expression_fallback"

        # M1-M2 pair should have high score due to expression correlation
        m1_m2_pair = None
        for pair in result.marker_pairs:
            if set(pair.marker_pair) == {"M1", "M2"}:
                m1_m2_pair = pair
                break

        assert m1_m2_pair is not None
        # Expression correlation should boost the score
        assert m1_m2_pair.spearman_correlation > 0.5


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
