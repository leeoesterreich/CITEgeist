"""Tests for nucleus patch extraction with global normalization."""
import numpy as np
import pandas as pd
import pytest

from CITEgeist.model.patch_extraction import (
    extract_patch,
    compute_global_stats,
    extract_patch_with_size,
)


class TestComputeGlobalStats:
    """Test global stats computation."""

    def test_percentile_method(self):
        """Test percentile normalization stats."""
        # 2-channel image with known values
        image = np.array([
            [[0, 100, 200], [50, 150, 250]],  # Channel 0: range 0-250
            [[10, 20, 30], [40, 50, 60]],      # Channel 1: range 10-60
        ], dtype=np.float32)

        stats = compute_global_stats(image, norm_method="percentile")

        assert stats["method"] == "percentile"
        assert stats["p1"].shape == (2,)
        assert stats["p99"].shape == (2,)
        # 1st percentile should be near min, 99th near max
        assert stats["p1"][0] < stats["p99"][0]
        assert stats["p1"][1] < stats["p99"][1]

    def test_minmax_method(self):
        """Test minmax normalization stats."""
        image = np.array([
            [[0, 100, 200], [50, 150, 250]],
            [[10, 20, 30], [40, 50, 60]],
        ], dtype=np.float32)

        stats = compute_global_stats(image, norm_method="minmax")

        assert stats["method"] == "minmax"
        assert stats["min"][0] == 0
        assert stats["max"][0] == 250
        assert stats["min"][1] == 10
        assert stats["max"][1] == 60

    def test_invalid_method_raises(self):
        """Test that invalid method raises error."""
        image = np.zeros((2, 10, 10), dtype=np.float32)
        with pytest.raises(ValueError, match="Unknown norm_method"):
            compute_global_stats(image, norm_method="invalid")


class TestExtractPatch:
    """Test patch extraction with global normalization."""

    def test_requires_global_stats(self):
        """Test that missing global_stats raises error."""
        image = np.zeros((2, 100, 100), dtype=np.float32)
        bbox = (10, 10, 30, 30)

        with pytest.raises(ValueError, match="global_stats is required"):
            extract_patch(image, bbox, global_stats=None)

    def test_percentile_normalization_output_range(self):
        """Test that percentile normalization outputs [0, 1] range."""
        image = np.random.rand(2, 100, 100).astype(np.float32) * 1000
        bbox = (20, 20, 50, 50)
        stats = compute_global_stats(image, norm_method="percentile")

        patch = extract_patch(image, bbox, global_stats=stats)

        assert patch.shape == (2, 96, 96)
        assert patch.min() >= 0.0
        assert patch.max() <= 1.0

    def test_minmax_normalization_output_range(self):
        """Test that minmax normalization outputs [0, 1] range."""
        image = np.random.rand(2, 100, 100).astype(np.float32) * 1000
        bbox = (20, 20, 50, 50)
        stats = compute_global_stats(image, norm_method="minmax")

        patch = extract_patch(image, bbox, global_stats=stats)

        assert patch.shape == (2, 96, 96)
        assert patch.min() >= 0.0
        assert patch.max() <= 1.0

    def test_intensity_preserved_across_patches(self):
        """Test that intensity differences are preserved across patches."""
        # Create image with bright region (channel 0) and dim region
        image = np.zeros((2, 100, 100), dtype=np.float32)
        image[0, 10:30, 10:30] = 1000  # Bright patch
        image[0, 60:80, 60:80] = 100   # Dim patch
        image[1, :, :] = 500           # Uniform channel 1

        stats = compute_global_stats(image, norm_method="minmax")

        bright_patch = extract_patch(image, (10, 10, 30, 30), global_stats=stats)
        dim_patch = extract_patch(image, (60, 60, 80, 80), global_stats=stats)

        # Bright patch should have higher mean than dim patch
        assert bright_patch[0].mean() > dim_patch[0].mean()


class TestExtractPatchBasicFunctionality:
    """Test basic patch extraction functionality (with global stats)."""

    def test_extract_single_patch_shape(self):
        """Extracted patch should have correct shape."""
        # Fake 3-channel image (DAPI + 2 boundary markers)
        image = np.random.rand(3, 500, 500).astype(np.float32)
        stats = compute_global_stats(image, norm_method="percentile")

        # Nucleus bounding box (x_min, y_min, x_max, y_max)
        bbox = (100, 150, 120, 170)  # 20x20 nucleus

        patch = extract_patch(image, bbox, expansion=0.75, output_size=96, global_stats=stats)

        assert patch.shape == (3, 96, 96)
        assert patch.dtype == np.float32

    def test_extract_patch_with_expansion(self):
        """Expansion should capture context around nucleus."""
        # Create image with gradient to test expansion properly
        image = np.zeros((1, 200, 200), dtype=np.float32)
        # Nucleus region with high value
        image[0, 90:110, 90:110] = 1.0
        # Surrounding region with intermediate value (to create variance)
        image[0, 70:130, 70:130] = np.where(
            image[0, 70:130, 70:130] == 0, 0.3, image[0, 70:130, 70:130]
        )

        stats = compute_global_stats(image, norm_method="minmax")
        bbox = (90, 90, 110, 110)  # 20x20 nucleus

        # No expansion - just nucleus
        patch_no_exp = extract_patch(image, bbox, expansion=0.0, output_size=20, global_stats=stats)

        # With expansion - includes surrounding context
        patch_exp = extract_patch(image, bbox, expansion=1.0, output_size=60, global_stats=stats)

        # The expanded patch should capture more of the image area
        assert patch_exp.shape[1] == 60
        assert patch_no_exp.shape[1] == 20


class TestExtractPatchesForSpot:
    """Test batch extraction (with global stats)."""

    def test_extract_patches_for_spot(self):
        """Batch extraction should return correct shape."""
        from CITEgeist.model.patch_extraction import extract_patches_for_spot

        image = np.random.rand(3, 500, 500).astype(np.float32)
        stats = compute_global_stats(image, norm_method="percentile")

        nuclei_df = pd.DataFrame({
            'nucleus_id': [1, 2, 3],
            'bbox_x_min': [50, 150, 250],
            'bbox_y_min': [50, 150, 250],
            'bbox_x_max': [70, 170, 270],
            'bbox_y_max': [70, 170, 270],
        })

        patches = extract_patches_for_spot(
            image, nuclei_df, expansion=0.75, output_size=96, global_stats=stats
        )

        assert patches.shape == (3, 3, 96, 96)  # 3 nuclei, 3 channels, 96x96


class TestExtractPatchWithSize:
    """Test patch extraction with size features."""

    def test_returns_patch_and_size(self):
        """Test that function returns both patch and size features."""
        image = np.random.rand(2, 100, 100).astype(np.float32) * 1000
        bbox = (20, 20, 50, 50)  # 30x30 bbox
        stats = compute_global_stats(image, norm_method="percentile")

        patch, size_features = extract_patch_with_size(image, bbox, global_stats=stats)

        assert patch.shape == (2, 96, 96)
        assert size_features.shape == (3,)

    def test_size_features_are_log_transformed(self):
        """Test that size features are log1p transformed."""
        image = np.random.rand(2, 100, 100).astype(np.float32)
        bbox = (10, 10, 20, 30)  # w=10, h=20, area=200
        stats = compute_global_stats(image, norm_method="percentile")

        _, size_features = extract_patch_with_size(image, bbox, global_stats=stats)

        expected_log_w = np.log1p(10)
        expected_log_h = np.log1p(20)
        expected_log_area = np.log1p(200)

        np.testing.assert_almost_equal(size_features[0], expected_log_w, decimal=5)
        np.testing.assert_almost_equal(size_features[1], expected_log_h, decimal=5)
        np.testing.assert_almost_equal(size_features[2], expected_log_area, decimal=5)

    def test_different_bbox_sizes_give_different_features(self):
        """Test that different bbox sizes produce different features."""
        image = np.random.rand(2, 200, 200).astype(np.float32)
        stats = compute_global_stats(image, norm_method="percentile")

        small_bbox = (10, 10, 20, 20)  # 10x10
        large_bbox = (50, 50, 100, 100)  # 50x50

        _, small_size = extract_patch_with_size(image, small_bbox, global_stats=stats)
        _, large_size = extract_patch_with_size(image, large_bbox, global_stats=stats)

        # Large bbox should have larger size features
        assert large_size[0] > small_size[0]  # log_w
        assert large_size[1] > small_size[1]  # log_h
        assert large_size[2] > small_size[2]  # log_area
