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


class TestConstrainedAssignment:
    """Test Hungarian assignment with count constraints."""

    def test_hungarian_respects_counts(self):
        """Assignments should match count constraints exactly."""
        from CITEgeist.model.constrained_assignment import hungarian_assign

        # 10 samples, 3 types, counts [4, 3, 3]
        log_likes = np.random.rand(10, 3)
        counts = np.array([4, 3, 3])

        assignments = hungarian_assign(log_likes, counts)

        assert len(assignments) == 10
        assert (assignments == 0).sum() == 4
        assert (assignments == 1).sum() == 3
        assert (assignments == 2).sum() == 3

    def test_hungarian_maximizes_likelihood(self):
        """Should prefer high-likelihood assignments."""
        from CITEgeist.model.constrained_assignment import hungarian_assign

        # Clear preference: sample 0 strongly prefers type 0
        log_likes = np.array([
            [10.0, 0.0, 0.0],  # Sample 0 strongly prefers type 0
            [0.0, 5.0, 0.0],  # Sample 1 prefers type 1
            [0.0, 0.0, 5.0],  # Sample 2 prefers type 2
        ])
        counts = np.array([1, 1, 1])

        assignments = hungarian_assign(log_likes, counts)

        assert assignments[0] == 0  # Should get preferred type
        assert assignments[1] == 1
        assert assignments[2] == 2

    def test_hungarian_handles_count_mismatch(self):
        """Should adjust when counts don't match n_samples."""
        from CITEgeist.model.constrained_assignment import hungarian_assign

        log_likes = np.random.rand(5, 3)
        counts = np.array([2, 2, 2])  # Sum = 6, but only 5 samples

        assignments = hungarian_assign(log_likes, counts)

        assert len(assignments) == 5
        assert all(0 <= a < 3 for a in assignments)

    def test_random_assign_respects_counts(self):
        """Random baseline should also respect counts."""
        from CITEgeist.model.constrained_assignment import random_assign

        counts = np.array([4, 3, 3])

        assignments = random_assign(counts, n_samples=10)

        assert len(assignments) == 10
        assert (assignments == 0).sum() == 4
        assert (assignments == 1).sum() == 3
        assert (assignments == 2).sum() == 3

    def test_proportions_to_counts_sums_correctly(self):
        """Counts should sum to n_total and approximate proportions."""
        from CITEgeist.model.constrained_assignment import proportions_to_counts

        proportions = np.array([0.4, 0.35, 0.25])
        counts = proportions_to_counts(proportions, 10)  # positional arg

        assert counts.sum() == 10
        assert np.allclose(counts / counts.sum(), proportions, atol=0.15)

    def test_proportions_to_counts_handles_zero(self):
        """Should handle n_total=0 gracefully."""
        from CITEgeist.model.constrained_assignment import proportions_to_counts

        proportions = np.array([0.5, 0.3, 0.2])
        counts = proportions_to_counts(proportions, 0)  # positional arg

        assert counts.sum() == 0
        assert len(counts) == 3
