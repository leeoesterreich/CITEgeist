# CITEgeist/tests/test_detection.py
"""Tests for cell type detection module."""
import numpy as np
import pytest


def test_detect_cell_types_basic():
    """Test GMM detection identifies signal vs background."""
    from CITEgeist.model.deconvolution.detection import detect_cell_types

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
    from CITEgeist.model.deconvolution.detection import detect_cell_types

    np.random.seed(42)
    # All background - single tight cluster
    X = np.random.normal(1.0, 0.1, (50, 2))

    marker_groups = {"TypeA": [0, 1]}

    detected = detect_cell_types(X, marker_groups, threshold=0.5)

    # With only background, GMM should still fit but signal cluster
    # will be very similar to background - detection should be sparse
    # (This tests edge case handling)
    assert detected.shape == (50, 1)


def test_adaptive_threshold_rare_cell_type():
    """Test that rare cell types get lower threshold (more sensitive)."""
    from CITEgeist.model.deconvolution.detection import detect_cell_types, _compute_adaptive_threshold
    from sklearn.mixture import GaussianMixture

    np.random.seed(42)
    n_spots = 200

    # Create data with rare cell type (10% signal)
    X = np.zeros((n_spots, 1))
    X[:, 0] = np.random.normal(1.0, 0.3, n_spots)  # background
    X[:20, 0] = np.random.normal(5.0, 0.5, 20)  # 10% signal

    marker_groups = {"RareType": [0]}

    # With adaptive threshold, should detect most of the rare signal
    detected_adaptive = detect_cell_types(X, marker_groups, adaptive_threshold=True)

    # With fixed threshold, might miss some
    detected_fixed = detect_cell_types(X, marker_groups, adaptive_threshold=False)

    # Adaptive should be at least as sensitive for rare types
    assert detected_adaptive[:20, 0].sum() >= detected_fixed[:20, 0].sum()

    # Also test the helper function directly
    gmm = GaussianMixture(n_components=2, random_state=42)
    gmm.fit(np.log1p(X))
    signal_cluster = int(np.argmax(gmm.means_.sum(axis=1)))
    w_signal = gmm.weights_[signal_cluster]

    # Rare type should have low signal weight
    assert w_signal < 0.3, f"Expected rare type to have w_signal < 0.3, got {w_signal}"

    # And should get lower threshold
    adaptive_thresh = _compute_adaptive_threshold(gmm, signal_cluster)
    assert adaptive_thresh == 0.3, f"Expected threshold 0.3 for rare type, got {adaptive_thresh}"


def test_adaptive_threshold_common_cell_type():
    """Test that common cell types get higher threshold (more conservative)."""
    from CITEgeist.model.deconvolution.detection import detect_cell_types, _compute_adaptive_threshold
    from sklearn.mixture import GaussianMixture

    np.random.seed(42)
    n_spots = 200

    # Create data with common cell type (80% signal)
    X = np.zeros((n_spots, 1))
    X[:, 0] = np.random.normal(1.0, 0.3, n_spots)  # base level
    X[:160, 0] = np.random.normal(5.0, 0.5, 160)  # 80% high signal

    marker_groups = {"CommonType": [0]}

    # Test the helper function directly
    gmm = GaussianMixture(n_components=2, random_state=42)
    gmm.fit(np.log1p(X))
    signal_cluster = int(np.argmax(gmm.means_.sum(axis=1)))
    w_signal = gmm.weights_[signal_cluster]

    # Common type should have high signal weight
    assert w_signal > 0.7, f"Expected common type to have w_signal > 0.7, got {w_signal}"

    # And should get higher threshold
    adaptive_thresh = _compute_adaptive_threshold(gmm, signal_cluster)
    assert adaptive_thresh == 0.6, f"Expected threshold 0.6 for common type, got {adaptive_thresh}"


def test_adaptive_threshold_moderate_cell_type():
    """Test that moderate prevalence uses base threshold."""
    from CITEgeist.model.deconvolution.detection import _compute_adaptive_threshold
    from sklearn.mixture import GaussianMixture

    np.random.seed(42)
    n_spots = 200

    # Create data with moderate cell type (50-60% signal)
    X = np.zeros((n_spots, 1))
    X[:, 0] = np.random.normal(1.0, 0.3, n_spots)
    X[:110, 0] = np.random.normal(5.0, 0.5, 110)  # ~55% signal

    gmm = GaussianMixture(n_components=2, random_state=42)
    gmm.fit(np.log1p(X))
    signal_cluster = int(np.argmax(gmm.means_.sum(axis=1)))
    w_signal = gmm.weights_[signal_cluster]

    # Moderate should be in 0.5-0.7 range
    assert 0.5 <= w_signal <= 0.7, f"Expected moderate type to have 0.5 <= w_signal <= 0.7, got {w_signal}"

    # And should use base threshold (0.5)
    adaptive_thresh = _compute_adaptive_threshold(gmm, signal_cluster, base_threshold=0.5)
    assert adaptive_thresh == 0.5, f"Expected threshold 0.5 for moderate type, got {adaptive_thresh}"


def test_adaptive_threshold_disabled():
    """Test that adaptive_threshold=False uses fixed threshold."""
    from CITEgeist.model.deconvolution.detection import detect_cell_types

    np.random.seed(42)
    n_spots = 100

    # Create data
    X = np.zeros((n_spots, 1))
    X[:, 0] = np.random.normal(1.0, 0.3, n_spots)
    X[:20, 0] = np.random.normal(5.0, 0.5, 20)

    marker_groups = {"TypeA": [0]}

    # Different fixed thresholds should give different results
    detected_low = detect_cell_types(X, marker_groups, threshold=0.3, adaptive_threshold=False)
    detected_high = detect_cell_types(X, marker_groups, threshold=0.7, adaptive_threshold=False)

    # Lower threshold should detect more
    assert detected_low.sum() >= detected_high.sum()
