"""
Test negative marker functionality in Module 2d and Module 3.

This test validates that:
1. Module 2d correctly detects residual signal for shared markers
2. Module 2d creates hierarchical profiles with negative markers
3. Module 3 optimization correctly uses negative markers to separate overlapping populations
"""

import numpy as np
import pytest
import logging

# Set up logging
logging.basicConfig(level=logging.INFO)

# Import the new functions
import sys
sys.path.insert(0, "/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist")

from CITEgeist.model.spatial_colocalization import (
    detect_residual_signal,
    resolve_hierarchical_profiles,
    ResolvedProfile,
    HierarchicalProfileResult,
)
from CITEgeist.model.gurobi_impl import (
    optimize_cell_proportions_with_negatives,
    build_assignment_matrices_from_profiles,
)


def create_synthetic_emt_stromal_data(
    n_spots: int = 500,
    n_emt: int = 100,
    n_stromal: int = 100,
    noise_level: float = 0.1,
    seed: int = 42,
):
    """
    Create synthetic data mimicking EMT/Stromal separation problem.

    Creates three populations:
    - EMT cells: Vimentin HIGH, E-Cadherin HIGH
    - Stromal cells: Vimentin HIGH, E-Cadherin LOW
    - Other cells: Vimentin LOW, E-Cadherin LOW

    Returns:
        X: (n_spots, 2) marker expression matrix [Vimentin, E-Cadherin]
        ground_truth: (n_spots, 3) proportion matrix [EMT, Stromal, Other]
    """
    np.random.seed(seed)

    # Initialize
    n_other = n_spots - n_emt - n_stromal
    X = np.zeros((n_spots, 2))  # [Vimentin, E-Cadherin]
    ground_truth = np.zeros((n_spots, 3))  # [EMT, Stromal, Other]

    # EMT spots: high Vimentin, high E-Cadherin
    emt_indices = np.arange(0, n_emt)
    X[emt_indices, 0] = 0.8 + np.random.normal(0, noise_level, n_emt)  # Vimentin
    X[emt_indices, 1] = 0.8 + np.random.normal(0, noise_level, n_emt)  # E-Cadherin
    ground_truth[emt_indices, 0] = 1.0  # EMT

    # Stromal spots: high Vimentin, LOW E-Cadherin
    stromal_indices = np.arange(n_emt, n_emt + n_stromal)
    X[stromal_indices, 0] = 0.8 + np.random.normal(0, noise_level, n_stromal)  # Vimentin HIGH
    X[stromal_indices, 1] = 0.1 + np.random.normal(0, noise_level, n_stromal)  # E-Cadherin LOW
    ground_truth[stromal_indices, 1] = 1.0  # Stromal

    # Other spots: low both
    other_indices = np.arange(n_emt + n_stromal, n_spots)
    X[other_indices, 0] = 0.1 + np.random.normal(0, noise_level, n_other)
    X[other_indices, 1] = 0.1 + np.random.normal(0, noise_level, n_other)
    ground_truth[other_indices, 2] = 1.0  # Other

    # Clip to valid range
    X = np.clip(X, 0, 1)

    return X, ground_truth


class TestResidualSignalDetection:
    """Test Module 2d residual signal detection."""

    def test_detects_vimentin_residual(self):
        """Vimentin should have residual signal outside EMT profile."""
        X, _ = create_synthetic_emt_stromal_data(n_spots=500, n_emt=100, n_stromal=100)
        marker_names = ["Vimentin", "E-Cadherin"]

        # Simulated discovered profile from Module 2b: EMT = {Vimentin, E-Cadherin}
        profiles = [["Vimentin", "E-Cadherin"]]

        residual = detect_residual_signal(
            X, marker_names, profiles,
            residual_threshold=0.10,
            verbose=True,
        )

        # Vimentin should have residual signal (appears in Stromal spots without E-Cadherin)
        assert "Vimentin" in residual, "Vimentin should have residual signal"
        assert residual["Vimentin"]["profile_0"] > 0.10, "Residual fraction should be > 10%"

        # E-Cadherin should NOT have much residual (only appears with Vimentin in EMT)
        if "E-Cadherin" in residual:
            assert residual["E-Cadherin"]["profile_0"] < 0.30, \
                "E-Cadherin residual should be lower than Vimentin"


class TestHierarchicalProfileResolution:
    """Test Module 2d hierarchical profile resolution."""

    def test_creates_stromal_profile_with_negative(self):
        """Should create Stromal profile with E-Cadherin as negative marker."""
        X, _ = create_synthetic_emt_stromal_data(n_spots=500, n_emt=100, n_stromal=100)
        marker_names = ["Vimentin", "E-Cadherin"]
        profiles = [["Vimentin", "E-Cadherin"]]
        profile_names = ["EMT"]

        result = resolve_hierarchical_profiles(
            X, marker_names, profiles, profile_names,
            residual_threshold=0.10,
            verbose=True,
        )

        # Should create a new profile for Vimentin-only (Stromal)
        assert result.n_resolved_profiles > result.n_original_profiles, \
            "Should create additional profile(s)"

        # Check that a Vimentin-only profile exists with E-Cadherin as negative
        vimentin_only_found = False
        for name, profile in result.resolved_profiles.items():
            if "Vimentin" in profile.positive_markers and len(profile.positive_markers) == 1:
                vimentin_only_found = True
                assert "E-Cadherin" in profile.negative_markers, \
                    "Vimentin-only profile should have E-Cadherin as negative marker"
                break

        assert vimentin_only_found, "Should have created a Vimentin-only profile"

        print("\nResolved profiles:")
        print(result.summary())


class TestBuildAssignmentMatrices:
    """Test building assignment matrices from profile dict."""

    def test_positive_negative_assignment(self):
        """Test that positive/negative markers are correctly assigned."""
        marker_names = ["Vimentin", "E-Cadherin", "CD3E"]
        profile_dict = {
            "EMT": {"positive": ["Vimentin", "E-Cadherin"], "negative": []},
            "Stromal": {"positive": ["Vimentin"], "negative": ["E-Cadherin"]},
            "T cells": {"positive": ["CD3E"], "negative": []},
        }

        pos, neg, cell_types = build_assignment_matrices_from_profiles(
            marker_names, profile_dict
        )

        assert pos.shape == (3, 3), f"Expected (3, 3), got {pos.shape}"
        assert neg.shape == (3, 3), f"Expected (3, 3), got {neg.shape}"

        # Check EMT positive markers
        emt_idx = cell_types.index("EMT")
        vim_idx = marker_names.index("Vimentin")
        ecad_idx = marker_names.index("E-Cadherin")

        assert pos[vim_idx, emt_idx] == 1.0, "Vimentin should be positive for EMT"
        assert pos[ecad_idx, emt_idx] == 1.0, "E-Cadherin should be positive for EMT"

        # Check Stromal
        stromal_idx = cell_types.index("Stromal")
        assert pos[vim_idx, stromal_idx] == 1.0, "Vimentin should be positive for Stromal"
        assert neg[ecad_idx, stromal_idx] == 1.0, "E-Cadherin should be negative for Stromal"


class TestOptimizationWithNegatives:
    """Test Module 3 optimization with negative markers."""

    @pytest.mark.skipif(
        not _gurobi_available(),
        reason="Gurobi not available"
    )
    def test_separates_emt_from_stromal(self):
        """Optimization should correctly separate EMT from Stromal using negative markers."""
        X, ground_truth = create_synthetic_emt_stromal_data(
            n_spots=200, n_emt=50, n_stromal=50, noise_level=0.05
        )
        marker_names = ["Vimentin", "E-Cadherin"]

        # Profile dict with negative markers
        profile_dict = {
            "EMT": {"positive": ["Vimentin", "E-Cadherin"], "negative": []},
            "Stromal": {"positive": ["Vimentin"], "negative": ["E-Cadherin"]},
            "Other": {"positive": [], "negative": []},  # Catch-all
        }

        pos, neg, cell_types = build_assignment_matrices_from_profiles(
            marker_names, profile_dict
        )

        # Run optimization
        Y, beta, _ = optimize_cell_proportions_with_negatives(
            X, marker_names, pos, neg, cell_types,
            lambda_neg=1.0,
            lambda_reg=0.1,
            max_iterations=20,
            warn_only=True,
        )

        # Compute correlation with ground truth for EMT and Stromal
        emt_idx = cell_types.index("EMT")
        stromal_idx = cell_types.index("Stromal")

        emt_corr = np.corrcoef(Y[:, emt_idx], ground_truth[:, 0])[0, 1]
        stromal_corr = np.corrcoef(Y[:, stromal_idx], ground_truth[:, 1])[0, 1]

        print(f"\nEMT correlation with ground truth: {emt_corr:.3f}")
        print(f"Stromal correlation with ground truth: {stromal_corr:.3f}")

        # Both should have positive correlation with ground truth
        assert emt_corr > 0.5, f"EMT correlation too low: {emt_corr:.3f}"
        assert stromal_corr > 0.5, f"Stromal correlation too low: {stromal_corr:.3f}"

        # Check that EMT is not predicted where Stromal should be
        # In Stromal spots (indices 50-100), EMT should be low
        stromal_spots = np.arange(50, 100)
        mean_emt_in_stromal = Y[stromal_spots, emt_idx].mean()
        mean_stromal_in_stromal = Y[stromal_spots, stromal_idx].mean()

        print(f"Mean EMT proportion in Stromal spots: {mean_emt_in_stromal:.3f}")
        print(f"Mean Stromal proportion in Stromal spots: {mean_stromal_in_stromal:.3f}")

        assert mean_stromal_in_stromal > mean_emt_in_stromal, \
            "Stromal should be higher than EMT in Stromal spots"


def _gurobi_available():
    """Check if Gurobi is available."""
    try:
        import gurobipy as gp
        # Try to create a model to verify license
        model = gp.Model()
        model.dispose()
        return True
    except Exception:
        return False


if __name__ == "__main__":
    # Run tests manually for debugging
    print("=" * 60)
    print("Testing Residual Signal Detection")
    print("=" * 60)
    test = TestResidualSignalDetection()
    test.test_detects_vimentin_residual()
    print("PASSED\n")

    print("=" * 60)
    print("Testing Hierarchical Profile Resolution")
    print("=" * 60)
    test = TestHierarchicalProfileResolution()
    test.test_creates_stromal_profile_with_negative()
    print("PASSED\n")

    print("=" * 60)
    print("Testing Assignment Matrix Building")
    print("=" * 60)
    test = TestBuildAssignmentMatrices()
    test.test_positive_negative_assignment()
    print("PASSED\n")

    if _gurobi_available():
        print("=" * 60)
        print("Testing Optimization with Negative Markers")
        print("=" * 60)
        test = TestOptimizationWithNegatives()
        test.test_separates_emt_from_stromal()
        print("PASSED\n")
    else:
        print("Skipping optimization test (Gurobi not available)")

    print("\n" + "=" * 60)
    print("ALL TESTS PASSED")
    print("=" * 60)
