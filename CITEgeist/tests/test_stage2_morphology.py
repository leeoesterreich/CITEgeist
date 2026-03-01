"""Tests for Stage 2 morphology-based assignment."""
import numpy as np
import pytest


class TestMorphologyFeatures:
    """Test morphology feature extraction."""

    def test_extract_extended_features_shape(self):
        """Features should be (12,) array."""
        from CITEgeist.model.morphology_features import extract_extended_features

        # Create synthetic 2-channel patch (DAPI, boundary)
        patch = np.random.rand(2, 64, 64).astype(np.float32)
        features = extract_extended_features(patch)

        assert features.shape == (12,)
        assert features.dtype == np.float32

    def test_extract_extended_features_no_nan(self):
        """Features should not contain NaN for valid input."""
        from CITEgeist.model.morphology_features import extract_extended_features

        patch = np.random.rand(2, 64, 64).astype(np.float32) * 255
        features = extract_extended_features(patch)

        assert not np.any(np.isnan(features))

    def test_extract_extended_features_handles_zeros(self):
        """Should handle all-zero patches gracefully."""
        from CITEgeist.model.morphology_features import extract_extended_features

        patch = np.zeros((2, 64, 64), dtype=np.float32)
        features = extract_extended_features(patch)

        assert features.shape == (12,)
        # May have NaN for correlation, but should be replaced with 0
        assert not np.any(np.isnan(features))


class TestMorphologyClassifier:
    """Test GMM-based morphology classifier."""

    def test_fit_creates_gmm_per_class(self):
        """Fit should create one GMM per cell type."""
        from CITEgeist.model.morphology_classifier import MorphologyClassifier

        # Synthetic data: 100 samples, 12 features, 3 classes
        features = np.random.rand(100, 12).astype(np.float32)
        labels = np.array([0] * 40 + [1] * 35 + [2] * 25)
        cell_types = ["TypeA", "TypeB", "TypeC"]

        clf = MorphologyClassifier(cell_types=cell_types, n_components=2)
        clf.fit(features, labels)

        assert len(clf.gmms) == 3
        assert clf.scaler is not None

    def test_log_likelihood_shape(self):
        """log_likelihood should return (N, K) array."""
        from CITEgeist.model.morphology_classifier import MorphologyClassifier

        features = np.random.rand(100, 12).astype(np.float32)
        labels = np.array([0] * 40 + [1] * 35 + [2] * 25)
        cell_types = ["TypeA", "TypeB", "TypeC"]

        clf = MorphologyClassifier(cell_types=cell_types, n_components=2)
        clf.fit(features, labels)

        test_features = np.random.rand(10, 12).astype(np.float32)
        log_likes = clf.log_likelihood(test_features)

        assert log_likes.shape == (10, 3)

    def test_handles_missing_class(self):
        """Should handle cell types with zero samples gracefully."""
        from CITEgeist.model.morphology_classifier import MorphologyClassifier

        features = np.random.rand(100, 12).astype(np.float32)
        labels = np.array([0] * 50 + [1] * 50)  # No class 2
        cell_types = ["TypeA", "TypeB", "TypeC"]

        clf = MorphologyClassifier(cell_types=cell_types, n_components=2)
        clf.fit(features, labels)

        test_features = np.random.rand(5, 12).astype(np.float32)
        log_likes = clf.log_likelihood(test_features)

        assert log_likes.shape == (5, 3)
        # Class 2 should have low (fallback) likelihood
        assert np.all(log_likes[:, 2] <= log_likes[:, :2].max())
