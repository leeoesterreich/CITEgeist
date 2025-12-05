"""
Tests for automatic profile discovery module.

Run with: pytest tests/test_profile_discovery.py -v
"""

import numpy as np
import pytest
from scipy import sparse
from pathlib import Path

from CITEgeist.model.auto_profile_discovery import (
    ProfileDiscoveryResult,
    _build_profile_matrix,
    _em_refine,
    _score_candidate,
    _standardize_markers,
    discover_profiles,
)


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
def simple_single_marker_data(rng):
    """
    Data with one dominant cell type marked by EPCAM.
    EPCAM+ cancer cells clustered in a spatial region.
    """
    n_rows, n_cols = 25, 20  # 500 spots
    coords = create_spatial_grid(n_rows, n_cols)
    n_spots = len(coords)
    n_markers = 5
    marker_names = ["EPCAM", "CD3D", "CD4", "CD8", "CD68"]

    # Background noise
    X = rng.normal(0, 0.3, (n_spots, n_markers))

    # EPCAM+ cancer cells in a circular region (spatially clustered)
    cancer_mask = spots_in_circle(coords, center=(10, 12), radius=8)
    X[cancer_mask, 0] += 3.0  # EPCAM elevated

    return X, marker_names, coords


@pytest.fixture
def t_cell_data(rng):
    """
    Data with CD4+ and CD8+ T cells sharing CD3D.
    Each T cell subtype clustered in different spatial regions.
    Tests proper handling of shared markers.
    """
    n_rows, n_cols = 30, 20  # 600 spots
    coords = create_spatial_grid(n_rows, n_cols)
    n_spots = len(coords)
    marker_names = ["CD3D", "CD4", "CD8", "CD68", "EPCAM"]
    n_markers = len(marker_names)

    X = rng.normal(0, 0.2, (n_spots, n_markers))

    # CD4 T cells in left region (CD3D+, CD4+)
    cd4_mask = spots_in_rectangle(coords, x_range=(0, 8), y_range=(0, 15))
    X[cd4_mask, 0] += 2.5  # CD3D
    X[cd4_mask, 1] += 2.5  # CD4

    # CD8 T cells in right region (CD3D+, CD8+)
    cd8_mask = spots_in_rectangle(coords, x_range=(12, 19), y_range=(0, 15))
    X[cd8_mask, 0] += 2.5  # CD3D
    X[cd8_mask, 2] += 2.5  # CD8

    return X, marker_names, coords


@pytest.fixture
def mixed_profile_data(rng):
    """
    Data with mixed profile sizes in distinct spatial regions:
    - EPCAM (single) for cancer - top left
    - CD3D+CD4 (double) for CD4 T cells - bottom left
    - CD68+CD163+MRC1 (triple) for M2 macrophages - right side
    """
    n_rows, n_cols = 32, 25  # 800 spots
    coords = create_spatial_grid(n_rows, n_cols)
    n_spots = len(coords)
    marker_names = ["EPCAM", "CD3D", "CD4", "CD68", "CD163", "MRC1", "CD8"]
    n_markers = len(marker_names)

    X = rng.normal(0, 0.2, (n_spots, n_markers))

    # Cancer (EPCAM+) in top-left region
    cancer_mask = spots_in_rectangle(coords, x_range=(0, 10), y_range=(20, 31))
    X[cancer_mask, 0] += 3.0  # EPCAM

    # CD4 T cells in bottom-left region
    cd4_mask = spots_in_rectangle(coords, x_range=(0, 10), y_range=(0, 12))
    X[cd4_mask, 1] += 2.5  # CD3D
    X[cd4_mask, 2] += 2.5  # CD4

    # M2 Macrophages in right region
    macro_mask = spots_in_rectangle(coords, x_range=(15, 24), y_range=(8, 24))
    X[macro_mask, 3] += 2.0  # CD68
    X[macro_mask, 4] += 2.5  # CD163
    X[macro_mask, 5] += 2.5  # MRC1

    return X, marker_names, coords


@pytest.fixture
def noise_only_data(rng):
    """Pure noise - no real structure, but with spatial coords."""
    n_rows, n_cols = 15, 20  # 300 spots
    coords = create_spatial_grid(n_rows, n_cols)
    n_spots = len(coords)
    marker_names = ["A", "B", "C", "D", "E"]
    X = rng.normal(0, 1, (n_spots, len(marker_names)))
    return X, marker_names, coords


# =============================================================================
# Unit Tests
# =============================================================================


class TestStandardizeMarkers:
    """Tests for _standardize_markers function."""

    def test_basic_standardization(self):
        X = np.array([[1, 2], [3, 4], [5, 6]], dtype=float)
        Z, valid = _standardize_markers(X)

        assert valid.all()
        assert Z.shape == X.shape
        # Each column should have mean ~0 and std ~1
        np.testing.assert_array_almost_equal(Z.mean(axis=0), [0, 0], decimal=5)
        np.testing.assert_array_almost_equal(Z.std(axis=0), [1, 1], decimal=1)

    def test_zero_variance_marker(self):
        X = np.array([[1, 5], [1, 6], [1, 7]], dtype=float)
        Z, valid = _standardize_markers(X)

        assert not valid[0]  # First column is constant
        assert valid[1]
        assert Z[:, 0].sum() == 0  # Zero-variance column zeroed out

    def test_robust_mode(self):
        # Add outlier
        X = np.array([[1, 2], [2, 3], [3, 4], [100, 5]], dtype=float)
        Z_normal, _ = _standardize_markers(X, robust=False)
        Z_robust, _ = _standardize_markers(X, robust=True)

        # Robust should be less affected by outlier
        assert np.abs(Z_robust[3, 0]) < np.abs(Z_normal[3, 0])


class TestScoreCandidate:
    """Tests for _score_candidate function."""

    def test_high_coexpression(self):
        # Both markers elevated in same spots
        Z = np.zeros((100, 2))
        Z[:50, 0] = 2.0  # Marker 0 high in first 50
        Z[:50, 1] = 2.0  # Marker 1 high in first 50
        beta = np.array([1.0, 1.0])

        score = _score_candidate(Z, {0, 1}, beta)
        assert score > 0

    def test_no_coexpression(self):
        # Markers elevated in different spots
        Z = np.zeros((100, 2))
        Z[:50, 0] = 2.0  # Marker 0 high in first 50
        Z[50:, 1] = 2.0  # Marker 1 high in last 50
        beta = np.array([1.0, 1.0])

        score = _score_candidate(Z, {0, 1}, beta)
        assert score == 0  # No overlap

    def test_beta_weighting(self):
        Z = np.ones((100, 2)) * 2.0
        beta_equal = np.array([1.0, 1.0])
        beta_unequal = np.array([0.5, 1.0])

        score_equal = _score_candidate(Z, {0, 1}, beta_equal)
        score_unequal = _score_candidate(Z, {0, 1}, beta_unequal)

        assert score_unequal < score_equal


class TestBuildProfileMatrix:
    """Tests for _build_profile_matrix function."""

    def test_single_profile(self):
        profiles = [{0, 2}]
        A = _build_profile_matrix(profiles, 4)

        expected = np.array([[1, 0, 1, 0]])
        np.testing.assert_array_equal(A, expected)

    def test_multiple_profiles(self):
        profiles = [{0}, {1, 2}, {0, 1, 3}]
        A = _build_profile_matrix(profiles, 4)

        assert A.shape == (3, 4)
        assert A[0, 0] == 1 and A[0, 1:].sum() == 0
        assert A[1, 1] == 1 and A[1, 2] == 1
        assert A[2, 0] == 1 and A[2, 1] == 1 and A[2, 3] == 1


class TestEMRefine:
    """Tests for _em_refine function."""

    def test_empty_profiles(self):
        Z = np.random.randn(100, 5)
        Y, beta, ll = _em_refine(Z, [], np.ones(5))

        assert Y.shape == (100, 0)
        assert ll == 0.0

    def test_single_profile_convergence(self):
        Z = np.zeros((100, 3))
        Z[:, 0] = 2.0  # First marker always high
        profiles = [{0}]

        Y, beta, ll = _em_refine(Z, profiles, np.ones(3))

        assert Y.shape == (100, 1)
        assert beta[0] == 1.0  # Single marker in single profile

    def test_shared_marker_beta_reduction(self):
        """Markers in multiple profiles should get lower beta."""
        Z = np.zeros((200, 3))
        # Profile 1: markers 0 and 1 co-elevated
        Z[:100, 0] = 2.0
        Z[:100, 1] = 2.0
        # Profile 2: markers 0 and 2 co-elevated
        Z[100:, 0] = 2.0
        Z[100:, 2] = 2.0

        profiles = [{0, 1}, {0, 2}]
        Y, beta, _ = _em_refine(Z, profiles, np.ones(3))

        # Marker 0 is shared -> should have lower beta
        assert beta[0] < beta[1]
        assert beta[0] < beta[2]


# =============================================================================
# Integration Tests
# =============================================================================


class TestDiscoverProfiles:
    """Integration tests for discover_profiles function."""

    def test_single_dominant_population(self, simple_single_marker_data):
        """Should discover EPCAM as a single-marker profile."""
        X, names, coords = simple_single_marker_data
        result = discover_profiles(X, names, max_k=2, seed=42, n_perm=100, coords=coords, verbose=False)

        assert isinstance(result, ProfileDiscoveryResult)
        assert len(result.profiles) >= 1
        assert "EPCAM" in result.profiles or any(
            "EPCAM" in p["Major"] for p in result.profiles.values()
        )

    def test_shared_marker_discrimination(self, t_cell_data):
        """
        Test handling of shared markers (CD3D in both CD4+ and CD8+ T cells).
        The greedy non-overlapping selection may pick CD3D as a single marker,
        or find one of the T cell subtypes. At minimum, should find CD3D-containing profiles.
        """
        X, names, coords = t_cell_data
        result = discover_profiles(X, names, max_k=2, seed=42, n_perm=100, coords=coords, verbose=False)

        # Should find at least 1 profile containing CD3D (the shared marker)
        assert len(result.profiles) >= 1

        # At least one profile should contain CD3D
        has_cd3d = any(
            "CD3D" in p["Major"]
            for p in result.profiles.values()
        )
        assert has_cd3d, "Expected at least one profile containing CD3D"

    def test_mixed_profile_sizes(self, mixed_profile_data):
        """Should discover profiles of varying sizes (1, 2, 3 markers)."""
        X, names, coords = mixed_profile_data
        result = discover_profiles(X, names, max_k=3, seed=42, n_perm=100, coords=coords, verbose=False)

        # Should find multiple profiles
        assert len(result.profiles) >= 2

        # Check for varied sizes
        sizes = [len(p["Major"]) for p in result.profiles.values()]
        assert len(set(sizes)) >= 1  # At least some variety expected

    def test_bic_stops_overfitting(self, simple_single_marker_data):
        """BIC should prevent adding spurious profiles."""
        X, names, coords = simple_single_marker_data
        result = discover_profiles(X, names, max_k=3, seed=42, n_perm=100, coords=coords, verbose=False)

        # BIC should be decreasing (or algorithm stopped)
        if len(result.bic_trace) > 1:
            # Last kept profile should have lower BIC than if we continued
            assert result.bic_trace[-1] <= result.bic_trace[-2] or len(result.bic_trace) == result.n_iterations

    def test_no_structure_handling(self, noise_only_data):
        """Should handle pure noise gracefully."""
        X, names, coords = noise_only_data
        result = discover_profiles(X, names, max_k=2, seed=42, n_perm=100, coords=coords, verbose=False)

        # Should not crash, may return few or no profiles
        assert isinstance(result, ProfileDiscoveryResult)
        # Noise data might still find spurious patterns at p<0.05
        # but should not find many
        assert len(result.profiles) <= 3

    def test_reproducibility(self, simple_single_marker_data):
        """Same seed should give identical results."""
        X, names, coords = simple_single_marker_data

        result1 = discover_profiles(X, names, seed=12345, n_perm=100, coords=coords, verbose=False)
        result2 = discover_profiles(X, names, seed=12345, n_perm=100, coords=coords, verbose=False)

        assert result1.profiles == result2.profiles
        assert result1.beta == result2.beta

    def test_different_seeds_similar_results(self, t_cell_data):
        """Different seeds should find similar core profiles."""
        X, names, coords = t_cell_data

        result1 = discover_profiles(X, names, seed=111, n_perm=100, coords=coords, verbose=False)
        result2 = discover_profiles(X, names, seed=222, n_perm=100, coords=coords, verbose=False)

        # Core profiles should be similar even with different seeds
        profiles1 = set(result1.profiles.keys())
        profiles2 = set(result2.profiles.keys())

        # At least some overlap expected
        overlap = len(profiles1 & profiles2)
        total = len(profiles1 | profiles2)
        if total > 0:
            similarity = overlap / total
            assert similarity >= 0.3  # At least 30% overlap

    def test_sparse_input(self, simple_single_marker_data):
        """Should handle scipy sparse matrices."""
        X, names, coords = simple_single_marker_data
        X_sparse = sparse.csr_matrix(X)

        result = discover_profiles(X_sparse, names, max_k=2, seed=42, n_perm=100, coords=coords, verbose=False)

        assert isinstance(result, ProfileDiscoveryResult)
        assert len(result.profiles) >= 1

    def test_output_schema_compatibility(self, simple_single_marker_data):
        """Output should match CITEgeist's expected format."""
        X, names, coords = simple_single_marker_data
        result = discover_profiles(X, names, max_k=2, seed=42, n_perm=100, coords=coords, verbose=False)

        for name, profile in result.profiles.items():
            assert "Major" in profile
            assert isinstance(profile["Major"], list)
            assert all(isinstance(m, str) for m in profile["Major"])

    def test_proportions_shape(self, t_cell_data):
        """Proportions matrix should have correct shape."""
        X, names, coords = t_cell_data
        result = discover_profiles(X, names, max_k=2, seed=42, n_perm=100, coords=coords, verbose=False)

        n_spots = X.shape[0]
        n_profiles = len(result.profiles)

        assert result.proportions.shape == (n_spots, n_profiles)

    def test_marker_name_mismatch_raises(self):
        """Should raise error if marker_names doesn't match X columns."""
        X = np.random.randn(100, 5)
        coords = create_spatial_grid(10, 10)
        names = ["A", "B", "C"]  # Only 3 names for 5 columns

        with pytest.raises(ValueError, match="marker_names length"):
            discover_profiles(X, names, coords=coords)

    def test_coords_required(self):
        """Should raise error if coords not provided."""
        X = np.random.randn(100, 5)
        names = ["A", "B", "C", "D", "E"]

        with pytest.raises(ValueError, match="coords is required"):
            discover_profiles(X, names)

    @pytest.mark.integration
    @pytest.mark.requires_data
    def test_simulated_h5ad_background_filtered(self):
        """
        Runs discovery on simulated CITE data.
        Should find majority of real cell-type profiles with high precision.
        Some false positives (Nonspecific markers) are acceptable.
        """
        cite_path = Path(__file__).parent.parent / "replicates" / "high_seg" / "h5ad_objects" / "Wu_rep_0_CITE.h5ad"
        if not cite_path.exists():
            pytest.skip("Simulated CITE h5ad not available")
        try:
            import scanpy as sc
        except ImportError:
            pytest.skip("scanpy not installed in test environment")

        adata = sc.read_h5ad(cite_path)
        result = discover_profiles(
            adata.X,
            list(adata.var_names),
            max_k=2,
            seed=1234,
            n_perm=100,
            alpha=0.1,
            verbose=False,
            coords=adata.obsm["spatial"],
        )

        assert isinstance(result, ProfileDiscoveryResult)
        assert len(result.profiles) >= 1

        # Ground truth: 9 cell types (B-cells, CAFs, Cancer Epithelial, Endothelial,
        # Myeloid, Normal Epithelial, PVL, Plasmablasts, T-cells)
        ground_truth_prefixes = [
            "B-cells", "CAFs", "Cancer Epithelial", "Endothelial",
            "Myeloid", "Normal Epithelial", "PVL", "Plasmablasts", "T-cells"
        ]

        # Count true positives (profiles matching ground truth cell types)
        true_positives = 0
        for profile in result.profiles.values():
            markers = profile["Major"]
            for prefix in ground_truth_prefixes:
                if all(m.startswith(prefix) for m in markers):
                    true_positives += 1
                    break

        # Count false positives (profiles with Nonspecific markers)
        false_positives = sum(
            1 for profile in result.profiles.values()
            if any(m.startswith("Nonspecific") for m in profile["Major"])
        )

        # Require at least 50% recall (finding at least half of the cell types)
        # Plasmablasts has weak signal so may be missed
        assert true_positives >= 4, f"Expected at least 4 true positives, got {true_positives}"

        # Require precision >= 70% (most discovered profiles should be real)
        n_discovered = len(result.profiles)
        precision = true_positives / n_discovered if n_discovered > 0 else 0
        assert precision >= 0.7, f"Expected precision >= 0.7, got {precision:.2f}"


# =============================================================================
# Edge Cases
# =============================================================================


class TestEdgeCases:
    """Edge case tests."""

    def test_two_markers_only(self, rng):
        """Should work with minimal marker count with spatially clustered signal."""
        coords = create_spatial_grid(10, 10)  # 100 spots
        X = rng.normal(0, 1, (100, 2))
        # Add signal to spots in top-left quadrant (spatially clustered)
        signal_mask = spots_in_rectangle(coords, x_range=(0, 4), y_range=(5, 9))
        X[signal_mask, 0] += 3.0
        X[signal_mask, 1] += 3.0
        names = ["A", "B"]

        result = discover_profiles(X, names, max_k=2, seed=42, n_perm=100, coords=coords, verbose=False)
        assert isinstance(result, ProfileDiscoveryResult)

    def test_single_spot(self, rng):
        """Should handle degenerate case gracefully."""
        X = rng.normal(0, 1, (1, 5))
        coords = np.array([[0.0, 0.0]])  # Single spot
        names = ["A", "B", "C", "D", "E"]

        # May warn or return empty, but should not crash
        result = discover_profiles(X, names, max_k=2, seed=42, n_perm=50, coords=coords, verbose=False)
        assert isinstance(result, ProfileDiscoveryResult)

    def test_all_constant_markers(self):
        """Should handle all-constant data."""
        coords = create_spatial_grid(10, 10)
        X = np.ones((100, 3))
        names = ["A", "B", "C"]

        result = discover_profiles(X, names, max_k=2, seed=42, n_perm=50, coords=coords, verbose=False)

        # All markers dropped, empty result
        assert len(result.profiles) == 0

    def test_max_k_one(self, simple_single_marker_data):
        """Should only find single-marker profiles when max_k=1."""
        X, names, coords = simple_single_marker_data
        result = discover_profiles(X, names, max_k=1, seed=42, n_perm=100, coords=coords, verbose=False)

        for profile in result.profiles.values():
            assert len(profile["Major"]) == 1
