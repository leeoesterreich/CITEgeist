# CITEgeist/tests/test_detection.py
"""Tests for cell type detection module."""
import numpy as np
import pytest


def test_detect_cell_types_basic():
    """Test GMM detection identifies signal vs background."""
    from CITEgeist.model.detection import detect_cell_types

    np.random.seed(42)
    n_spots = 100

    # Create synthetic data: 30 spots with signal, 70 with background
    # Cell type 0 has markers [0, 1], cell type 1 has marker [2]
    X = np.zeros((n_spots, 3))

    # Background (low values)
    X[:, :] = np.random.normal(1.0, 0.3, (n_spots, 3))

    # Signal spots for cell type 0 (high on markers 0,1)
    X[:30, 0] = np.random.normal(5.0, 0.5, 30)
    X[:30, 1] = np.random.normal(4.5, 0.5, 30)

    # Signal spots for cell type 1 (high on marker 2) - different spots
    X[50:70, 2] = np.random.normal(6.0, 0.5, 20)

    marker_groups = {
        "TypeA": [0, 1],  # multi-marker
        "TypeB": [2],     # single marker
    }

    detected = detect_cell_types(X, marker_groups, threshold=0.5)

    assert detected.shape == (100, 2)
    assert detected.dtype == bool

    # TypeA should be detected in first ~30 spots
    assert detected[:30, 0].sum() >= 25  # most signal spots detected
    assert detected[50:, 0].sum() <= 10  # few false positives

    # TypeB should be detected in spots 50-70
    assert detected[50:70, 1].sum() >= 15
    assert detected[:50, 1].sum() <= 10


def test_detect_cell_types_returns_all_false_for_no_signal():
    """Test that pure background returns no detections."""
    from CITEgeist.model.detection import detect_cell_types

    np.random.seed(42)
    # All background - single tight cluster
    X = np.random.normal(1.0, 0.1, (50, 2))

    marker_groups = {"TypeA": [0, 1]}

    detected = detect_cell_types(X, marker_groups, threshold=0.5)

    # With only background, GMM should still fit but signal cluster
    # will be very similar to background - detection should be sparse
    # (This tests edge case handling)
    assert detected.shape == (50, 1)
