# CITEgeist/tests/test_detection_estimation.py
"""Tests for combined detection + estimation pipeline."""
import numpy as np
import pytest


def test_solve_detection_estimation_full_pipeline():
    """Test full pipeline: detection -> estimation with learned variance."""
    pytest.importorskip("gurobipy")
    from CITEgeist.model.detection_estimation import solve_detection_estimation

    np.random.seed(42)
    n_spots = 50
    n_markers = 4
    n_types = 2

    # Profile: Type 0 has markers [0,1], Type 1 has markers [2,3]
    profile = np.array([
        [1, 1, 0, 0],
        [0, 0, 1, 1],
    ], dtype=float)

    marker_groups = {
        "TypeA": [0, 1],
        "TypeB": [2, 3],
    }

    # Ground truth: Type 0 in spots 0-24, Type 1 in spots 25-49
    # 10 nuclei per spot
    nuclei_counts = np.full(n_spots, 10)

    # Generate signal with baseline
    true_alpha = np.array([1.0, 1.0, 1.5, 1.5])
    true_beta = np.array([0.5, 0.4, 0.6, 0.5])

    X = np.zeros((n_spots, n_markers))

    # Background everywhere
    X[:, :] = true_alpha

    # Type 0 signal in first 25 spots
    X[:25, 0] += 10 * true_beta[0]  # 10 cells * beta
    X[:25, 1] += 10 * true_beta[1]

    # Type 1 signal in last 25 spots
    X[25:, 2] += 10 * true_beta[2]
    X[25:, 3] += 10 * true_beta[3]

    # Add noise
    X += np.random.normal(0, 0.2, X.shape)

    detected, counts, alpha, beta, sigma_sq = solve_detection_estimation(
        X, nuclei_counts, profile, marker_groups, max_iter=5
    )

    # Check shapes
    assert detected.shape == (n_spots, n_types)
    assert counts.shape == (n_spots, n_types)
    assert alpha.shape == (n_markers,)
    assert beta.shape == (n_markers,)
    assert sigma_sq.shape == (n_markers,)

    # Check detection pattern
    # Type 0 should be detected mostly in first 25 spots
    assert detected[:25, 0].sum() >= 20
    assert detected[25:, 0].sum() <= 10

    # Type 1 should be detected mostly in last 25 spots
    assert detected[25:, 1].sum() >= 20
    assert detected[:25, 1].sum() <= 10

    # Check that counts respect detection mask
    for i in range(n_spots):
        for k in range(n_types):
            if not detected[i, k]:
                assert counts[i, k] == 0

    # Check nuclei sum (where detected)
    for i in range(n_spots):
        if detected[i].any():
            assert counts[i].sum() == nuclei_counts[i]


def test_solve_detection_estimation_learns_reasonable_sigma():
    """Test that learned sigma_sq reflects actual noise level."""
    pytest.importorskip("gurobipy")
    from CITEgeist.model.detection_estimation import solve_detection_estimation

    np.random.seed(42)
    n_spots = 100
    n_markers = 2
    n_types = 1

    profile = np.array([[1, 1]], dtype=float)
    marker_groups = {"TypeA": [0, 1]}

    nuclei_counts = np.full(n_spots, 5)

    # Different noise levels per marker
    true_sigma = np.array([0.5, 2.0])

    X = np.zeros((n_spots, n_markers))
    X[:, 0] = 1.0 + 5 * 0.5 + np.random.normal(0, true_sigma[0], n_spots)
    X[:, 1] = 1.0 + 5 * 0.5 + np.random.normal(0, true_sigma[1], n_spots)

    detected, counts, alpha, beta, sigma_sq = solve_detection_estimation(
        X, nuclei_counts, profile, marker_groups, max_iter=10
    )

    # Learned sigma_sq should reflect that marker 1 is noisier
    assert sigma_sq[1] > sigma_sq[0], "Marker 1 should have higher variance"
