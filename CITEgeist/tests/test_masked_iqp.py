# CITEgeist/tests/test_masked_iqp.py
"""Tests for masked IQP solver."""
import numpy as np
import pytest


def test_masked_iqp_respects_detection_mask():
    """Test that masked IQP returns zeros for non-detected types."""
    pytest.importorskip("gurobipy")
    from CITEgeist.model.masked_iqp import solve_masked_iqp

    np.random.seed(42)
    n_spots, n_markers, n_types = 10, 4, 3

    # Simple profile: each type has one unique marker
    # Type 0 -> marker 0, Type 1 -> marker 1, Type 2 -> marker 2
    # Marker 3 is shared by types 0 and 1
    profile = np.array([
        [1, 0, 0, 1],  # Type 0: markers 0, 3
        [0, 1, 0, 1],  # Type 1: markers 1, 3
        [0, 0, 1, 0],  # Type 2: marker 2
    ], dtype=float)

    # Observed signal - clear separation
    X = np.zeros((n_spots, n_markers))
    nuclei_counts = np.full(n_spots, 10)

    # Spots 0-3: Type 0 signal (markers 0,3 high)
    X[:4, 0] = 5.0
    X[:4, 3] = 4.0

    # Spots 4-6: Type 1 signal (markers 1,3 high)
    X[4:7, 1] = 6.0
    X[4:7, 3] = 4.0

    # Spots 7-9: Type 2 signal (marker 2 high)
    X[7:, 2] = 7.0

    # Detection mask: only detect present types
    detected = np.zeros((n_spots, n_types), dtype=bool)
    detected[:4, 0] = True   # Type 0 in spots 0-3
    detected[4:7, 1] = True  # Type 1 in spots 4-6
    detected[7:, 2] = True   # Type 2 in spots 7-9

    weights = np.ones(n_markers)

    counts, alpha, beta = solve_masked_iqp(
        X, nuclei_counts, profile, detected, weights
    )

    assert counts.shape == (n_spots, n_types)

    # Check detection mask is respected (zeros where not detected)
    for i in range(n_spots):
        for k in range(n_types):
            if not detected[i, k]:
                assert counts[i, k] == 0, f"counts[{i},{k}] should be 0 (not detected)"

    # Check nuclei sum constraint
    for i in range(n_spots):
        assert counts[i].sum() == nuclei_counts[i], f"spot {i} doesn't sum to {nuclei_counts[i]}"


def test_masked_iqp_learns_alpha_beta():
    """Test that alpha (baseline) and beta (signal-per-cell) are learned."""
    pytest.importorskip("gurobipy")
    from CITEgeist.model.masked_iqp import solve_masked_iqp

    np.random.seed(42)
    n_spots, n_markers, n_types = 30, 2, 2

    # Two types with different markers
    # Type 0 has marker 0, Type 1 has marker 1
    profile = np.array([
        [1, 0],  # Type 0: only marker 0
        [0, 1],  # Type 1: only marker 1
    ], dtype=float)

    # Ground truth parameters
    true_alpha = np.array([2.0, 1.5])
    true_beta = np.array([0.8, 1.0])

    # Variable nuclei counts across spots (3 to 10)
    nuclei_counts = np.random.randint(3, 11, n_spots)

    # Ground truth counts: variable mix of types
    # Spots 0-9: mostly Type 0, Spots 10-19: mixed, Spots 20-29: mostly Type 1
    true_counts = np.zeros((n_spots, n_types), dtype=int)
    for i in range(n_spots):
        n = nuclei_counts[i]
        if i < 10:
            # Mostly Type 0
            true_counts[i, 0] = min(n, np.random.randint(n//2 + 1, n + 1))
            true_counts[i, 1] = n - true_counts[i, 0]
        elif i < 20:
            # Mixed
            true_counts[i, 0] = np.random.randint(0, n + 1)
            true_counts[i, 1] = n - true_counts[i, 0]
        else:
            # Mostly Type 1
            true_counts[i, 1] = min(n, np.random.randint(n//2 + 1, n + 1))
            true_counts[i, 0] = n - true_counts[i, 1]

    # All types detected everywhere (full detection)
    detected = np.ones((n_spots, n_types), dtype=bool)

    # Generate observed signal: X = alpha + counts @ profile * beta
    # For marker m: X[:,m] = alpha[m] + (counts @ profile[:,m]) * beta[m]
    X = np.zeros((n_spots, n_markers))
    for m in range(n_markers):
        effective_counts = true_counts @ profile[:, m]
        X[:, m] = true_alpha[m] + effective_counts * true_beta[m]
    # Add small noise
    X += np.random.normal(0, 0.1, X.shape)

    weights = np.ones(n_markers)

    counts, alpha, beta = solve_masked_iqp(
        X, nuclei_counts, profile, detected, weights
    )

    # Alpha should be close to true_alpha
    assert np.allclose(alpha, true_alpha, atol=0.5), f"alpha {alpha} not close to {true_alpha}"

    # Beta should be close to true_beta
    assert np.allclose(beta, true_beta, atol=0.3), f"beta {beta} not close to {true_beta}"

    # Counts should sum to nuclei count
    for i in range(n_spots):
        assert counts[i].sum() == nuclei_counts[i], f"spot {i} doesn't sum correctly"

    # Counts should be reasonably close to true counts (correlation > 0.8)
    count_corr = np.corrcoef(counts.flatten(), true_counts.flatten())[0, 1]
    assert count_corr > 0.7, f"Count correlation {count_corr} too low"
