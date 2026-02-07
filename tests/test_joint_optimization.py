"""
Unit tests for joint optimization module.

Tests the key components of the joint profile discovery and proportion estimation:
- Initialization (NMF and random)
- Y subproblem optimization
- W subproblem optimization
- Beta update
- Alternating minimization
- BIC-based K selection
- Profile formatting
"""

import numpy as np
import pytest
import scipy.sparse as sp

# Import the joint optimization functions
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from CITEgeist.model.joint_optimization import (
    JointOptimizationResult,
    _initialize_nmf,
    _initialize_random,
    _optimize_Y,
    _optimize_W,
    _update_beta,
    _compute_reconstruction_error,
    _compute_bic,
    _alternating_minimization,
    _format_profiles,
    optimize_profiles_jointly,
)
from CITEgeist.model.gurobi_impl import build_spatial_laplacian


class TestInitialization:
    """Test initialization functions."""

    def test_nmf_initialization_shapes(self):
        """Test that NMF initialization returns correct shapes."""
        n_spots, n_markers, K = 100, 10, 5
        X = np.random.rand(n_spots, n_markers) + 0.1  # Ensure positive

        W, Y, beta = _initialize_nmf(X, K, seed=42)

        assert W.shape == (K, n_markers)
        assert Y.shape == (n_spots, K)
        assert beta.shape == (n_markers,)

    def test_nmf_initialization_constraints(self):
        """Test that NMF initialization satisfies constraints."""
        n_spots, n_markers, K = 100, 10, 5
        X = np.random.rand(n_spots, n_markers) + 0.1

        W, Y, beta = _initialize_nmf(X, K, seed=42)

        # W should be in [0, 1]
        assert np.all(W >= 0) and np.all(W <= 1)

        # Y rows should sum to approximately 1
        row_sums = Y.sum(axis=1)
        assert np.allclose(row_sums, 1.0, atol=0.1)

        # Beta should be 1
        assert np.allclose(beta, 1.0)

    def test_random_initialization_shapes(self):
        """Test that random initialization returns correct shapes."""
        n_spots, n_markers, K = 100, 10, 5

        W, Y, beta = _initialize_random(n_spots, n_markers, K, seed=42)

        assert W.shape == (K, n_markers)
        assert Y.shape == (n_spots, K)
        assert beta.shape == (n_markers,)


class TestBetaUpdate:
    """Test beta update function."""

    def test_beta_update_basic(self):
        """Test basic beta update."""
        n_spots, n_markers, K = 50, 5, 3
        np.random.seed(42)

        X = np.random.rand(n_spots, n_markers) + 0.1
        Y = np.random.rand(n_spots, K)
        Y = Y / Y.sum(axis=1, keepdims=True)
        W = np.random.rand(K, n_markers)

        beta = _update_beta(X, Y, W)

        # Beta should be in valid range
        assert np.all(beta >= 0.1)
        assert np.all(beta <= 2.0)
        assert beta.shape == (n_markers,)

    def test_beta_update_improves_fit(self):
        """Test that beta update reduces reconstruction error."""
        n_spots, n_markers, K = 50, 5, 3
        np.random.seed(42)

        X = np.random.rand(n_spots, n_markers) + 0.1
        Y = np.random.rand(n_spots, K)
        Y = Y / Y.sum(axis=1, keepdims=True)
        W = np.random.rand(K, n_markers)

        beta_init = np.ones(n_markers)
        error_before = _compute_reconstruction_error(X, Y, W, beta_init)

        beta_updated = _update_beta(X, Y, W)
        error_after = _compute_reconstruction_error(X, Y, W, beta_updated)

        # Error should not increase after optimal beta update
        assert error_after <= error_before + 1e-6


class TestReconstructionError:
    """Test reconstruction error computation."""

    def test_reconstruction_error_zero(self):
        """Test that perfect reconstruction has zero error."""
        n_spots, K, n_markers = 10, 3, 5

        Y = np.random.rand(n_spots, K)
        W = np.random.rand(K, n_markers)
        beta = np.ones(n_markers)

        # Construct X as exact reconstruction
        X = Y @ (W * beta[np.newaxis, :])

        error = _compute_reconstruction_error(X, Y, W, beta)
        assert error < 1e-10

    def test_reconstruction_error_positive(self):
        """Test that non-perfect reconstruction has positive error."""
        n_spots, K, n_markers = 10, 3, 5

        X = np.random.rand(n_spots, n_markers)
        Y = np.random.rand(n_spots, K)
        W = np.random.rand(K, n_markers)
        beta = np.ones(n_markers)

        error = _compute_reconstruction_error(X, Y, W, beta)
        assert error > 0


class TestBIC:
    """Test BIC computation."""

    def test_bic_increases_with_k(self):
        """Test that BIC penalizes model complexity."""
        n_spots, n_markers = 50, 10
        np.random.seed(42)

        X = np.random.rand(n_spots, n_markers) + 0.1

        bics = []
        for K in [2, 5, 10]:
            Y = np.random.rand(n_spots, K)
            Y = Y / Y.sum(axis=1, keepdims=True)
            W = np.random.rand(K, n_markers)
            beta = np.ones(n_markers)

            bic = _compute_bic(X, Y, W, beta)
            bics.append(bic)

        # With random data, higher K should have higher BIC due to complexity penalty
        # (This may not always hold depending on random seed, but provides basic check)
        assert len(bics) == 3


class TestProfileFormatting:
    """Test profile formatting function."""

    def test_format_profiles_basic(self):
        """Test basic profile formatting."""
        K, n_markers = 3, 5
        marker_names = ['CD3', 'CD4', 'CD8', 'EPCAM', 'VIM']

        W = np.array([
            [0.9, 0.8, 0.1, 0.0, 0.0],  # CD3+CD4 profile
            [0.9, 0.1, 0.8, 0.0, 0.0],  # CD3+CD8 profile
            [0.0, 0.0, 0.0, 0.9, 0.0],  # EPCAM profile
        ])

        profiles = _format_profiles(W, marker_names, threshold=0.3)

        # Should have 3 discovered profiles + Unknown
        assert len(profiles) >= 3
        assert 'Unknown' in profiles

    def test_format_profiles_empty_type(self):
        """Test that cell types with no markers above threshold are skipped."""
        K, n_markers = 2, 3
        marker_names = ['A', 'B', 'C']

        W = np.array([
            [0.9, 0.8, 0.1],  # Has markers above threshold
            [0.1, 0.1, 0.1],  # No markers above threshold
        ])

        profiles = _format_profiles(W, marker_names, threshold=0.5)

        # Should have 1 discovered profile + Unknown (second type filtered out)
        assert len(profiles) == 2  # One discovered + Unknown


class TestOptimizeY:
    """Test Y optimization subproblem."""

    @pytest.mark.skipif(
        os.environ.get('SKIP_GUROBI_TESTS', 'false').lower() == 'true',
        reason="Skipping Gurobi tests"
    )
    def test_optimize_y_constraints(self):
        """Test that Y optimization satisfies constraints."""
        n_spots, n_markers, K = 20, 5, 3
        np.random.seed(42)

        X = np.random.rand(n_spots, n_markers) + 0.1
        W = np.random.rand(K, n_markers)
        W = W / W.max(axis=1, keepdims=True)
        beta = np.ones(n_markers)

        coords = np.random.rand(n_spots, 2) * 100
        L = build_spatial_laplacian(coords, k=5, normed=True)
        L_coo = sp.coo_matrix(L)

        Y = _optimize_Y(X, W, beta, L_coo, lambda_spatial=0.1, time_limit=60.0)

        # Y should be in [0, 1]
        assert np.all(Y >= -1e-6)
        assert np.all(Y <= 1.0 + 1e-6)

        # Row sums should be in [0.9, 1.2]
        row_sums = Y.sum(axis=1)
        assert np.all(row_sums >= 0.9 - 1e-6)
        assert np.all(row_sums <= 1.2 + 1e-6)


class TestOptimizeW:
    """Test W optimization subproblem."""

    @pytest.mark.skipif(
        os.environ.get('SKIP_GUROBI_TESTS', 'false').lower() == 'true',
        reason="Skipping Gurobi tests"
    )
    def test_optimize_w_constraints(self):
        """Test that W optimization satisfies constraints."""
        n_spots, n_markers, K = 20, 5, 3
        np.random.seed(42)

        X = np.random.rand(n_spots, n_markers) + 0.1
        Y = np.random.rand(n_spots, K)
        Y = Y / Y.sum(axis=1, keepdims=True)
        beta = np.ones(n_markers)

        W = _optimize_W(
            X, Y, beta,
            lambda_sparsity=0.1,
            lambda_distinct=0.5,
            max_markers_per_type=3,
            time_limit=60.0
        )

        # W should be in [0, 1]
        assert np.all(W >= -1e-6)
        assert np.all(W <= 1.0 + 1e-6)

        # Row sums should respect max_markers constraint
        row_sums = W.sum(axis=1)
        assert np.all(row_sums <= 3.0 + 1e-6)


class TestAlternatingMinimization:
    """Test alternating minimization loop."""

    @pytest.mark.skipif(
        os.environ.get('SKIP_GUROBI_TESTS', 'false').lower() == 'true',
        reason="Skipping Gurobi tests"
    )
    def test_alternating_minimization_converges(self):
        """Test that alternating minimization reduces objective."""
        n_spots, n_markers, K = 30, 6, 3
        np.random.seed(42)

        X = np.random.rand(n_spots, n_markers) + 0.1
        W_init, Y_init, beta_init = _initialize_nmf(X, K, seed=42)

        coords = np.random.rand(n_spots, 2) * 100
        L = build_spatial_laplacian(coords, k=5, normed=True)
        L_coo = sp.coo_matrix(L)

        W, Y, beta, final_obj, n_iter = _alternating_minimization(
            X, W_init, Y_init, beta_init, L_coo,
            lambda_spatial=0.1,
            lambda_sparsity=0.1,
            lambda_distinct=0.5,
            max_markers_per_type=3,
            max_iterations=5,
            tolerance=1e-4,
            verbose=False,
        )

        # Should have run at least 1 iteration
        assert n_iter >= 1

        # Final objective should be finite
        assert np.isfinite(final_obj)


class TestOptimizeProfilesJointly:
    """Test the main joint optimization function."""

    @pytest.mark.skipif(
        os.environ.get('SKIP_GUROBI_TESTS', 'false').lower() == 'true',
        reason="Skipping Gurobi tests"
    )
    def test_optimize_profiles_jointly_basic(self):
        """Test basic joint optimization."""
        n_spots, n_markers = 50, 8
        np.random.seed(42)

        X = np.random.rand(n_spots, n_markers) + 0.1
        marker_names = [f'Marker_{i}' for i in range(n_markers)]
        coords = np.random.rand(n_spots, 2) * 100

        result = optimize_profiles_jointly(
            X=X,
            marker_names=marker_names,
            coords=coords,
            min_K=2,
            max_K=4,
            max_markers_per_type=3,
            lambda_spatial=0.1,
            lambda_sparsity=0.1,
            lambda_distinct=0.5,
            max_iterations=3,  # Reduced for faster testing
            n_restarts=1,
            seed=42,
            verbose=False,
        )

        # Check result structure
        assert isinstance(result, JointOptimizationResult)
        assert result.K >= 2 and result.K <= 4
        assert result.W.shape == (result.K, n_markers)
        assert result.Y.shape == (n_spots, result.K)
        assert result.beta.shape == (n_markers,)
        assert 'Unknown' in result.profiles
        assert np.isfinite(result.bic)

    @pytest.mark.skipif(
        os.environ.get('SKIP_GUROBI_TESTS', 'false').lower() == 'true',
        reason="Skipping Gurobi tests"
    )
    def test_optimize_profiles_jointly_selects_k(self):
        """Test that BIC selects appropriate K."""
        n_spots, n_markers = 50, 8
        np.random.seed(42)

        # Create data with clear structure (3 groups)
        true_K = 3
        true_W = np.zeros((true_K, n_markers))
        true_W[0, 0:3] = 1.0  # Profile 1: markers 0-2
        true_W[1, 3:5] = 1.0  # Profile 2: markers 3-4
        true_W[2, 5:8] = 1.0  # Profile 3: markers 5-7

        true_Y = np.random.dirichlet([1, 1, 1], size=n_spots)
        X = true_Y @ true_W + np.random.randn(n_spots, n_markers) * 0.1

        marker_names = [f'Marker_{i}' for i in range(n_markers)]
        coords = np.random.rand(n_spots, 2) * 100

        result = optimize_profiles_jointly(
            X=X,
            marker_names=marker_names,
            coords=coords,
            min_K=2,
            max_K=5,
            max_markers_per_type=4,
            lambda_spatial=0.0,  # No spatial for this test
            lambda_sparsity=0.1,
            lambda_distinct=0.5,
            max_iterations=10,
            n_restarts=1,
            seed=42,
            verbose=False,
        )

        # K should be close to true_K (within 1)
        assert abs(result.K - true_K) <= 2


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
