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
