"""Integration test for new preprocessing pipeline.

Tests the full pipeline flow:
1. compute_global_stats → global normalization stats
2. extract_patch_with_size → patches + size features
3. Output file validation for downstream VAE/prototype training
"""
import tempfile
from pathlib import Path
import numpy as np
import pandas as pd
import pytest

from CITEgeist.model.patch_extraction import (
    compute_global_stats,
    extract_patch,
    extract_patch_with_size,
)


class TestPreprocessingPipeline:
    """Integration tests for preprocessing pipeline."""

    def test_global_stats_savez_roundtrip(self):
        """Test that global stats can be saved and loaded correctly."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)

            # Create test image
            image = np.random.rand(2, 200, 200).astype(np.float32) * 65535

            # Compute and save stats
            stats = compute_global_stats(image, norm_method="percentile")
            np.savez(tmpdir / "global_stats.npz", **stats)

            # Load and verify
            loaded = dict(np.load(tmpdir / "global_stats.npz"))
            assert loaded["method"] == "percentile"
            np.testing.assert_array_equal(loaded["p1"], stats["p1"])
            np.testing.assert_array_equal(loaded["p99"], stats["p99"])

    def test_patches_and_sizes_output_format(self):
        """Test that patches and sizes are saved in correct format."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)

            # Create test image
            image = np.random.rand(2, 200, 200).astype(np.float32) * 65535
            stats = compute_global_stats(image, norm_method="percentile")

            # Extract patches for multiple nuclei
            bboxes = [
                (20, 20, 40, 40),   # 20x20
                (60, 60, 90, 90),   # 30x30
                (100, 100, 140, 140),  # 40x40
            ]

            patches = []
            sizes = []
            for bbox in bboxes:
                patch, size = extract_patch_with_size(image, bbox, global_stats=stats)
                patches.append(patch)
                sizes.append(size)

            patches_array = np.stack(patches, axis=0)
            sizes_array = np.stack(sizes, axis=0)

            # Save as pipeline would
            np.save(tmpdir / "spot_0_patches.npy", patches_array)
            np.save(tmpdir / "spot_0_sizes.npy", sizes_array)

            # Verify format
            loaded_patches = np.load(tmpdir / "spot_0_patches.npy")
            loaded_sizes = np.load(tmpdir / "spot_0_sizes.npy")

            assert loaded_patches.shape == (3, 2, 96, 96)
            assert loaded_sizes.shape == (3, 3)
            assert loaded_patches.dtype == np.float32
            assert loaded_sizes.dtype == np.float32

    def test_size_features_discriminate_nucleus_sizes(self):
        """Test that size features can discriminate different nucleus sizes."""
        image = np.random.rand(2, 300, 300).astype(np.float32)
        stats = compute_global_stats(image, norm_method="percentile")

        # Small, medium, large nuclei
        small_bbox = (10, 10, 20, 20)      # 10x10
        medium_bbox = (50, 50, 80, 80)     # 30x30
        large_bbox = (100, 100, 170, 170)  # 70x70

        _, small_size = extract_patch_with_size(image, small_bbox, global_stats=stats)
        _, medium_size = extract_patch_with_size(image, medium_bbox, global_stats=stats)
        _, large_size = extract_patch_with_size(image, large_bbox, global_stats=stats)

        # Size features should be ordered: small < medium < large
        assert small_size[2] < medium_size[2] < large_size[2]  # log_area

    def test_intensity_preserved_across_patches(self):
        """Test that global normalization preserves intensity differences."""
        # Create image with known bright and dim regions
        image = np.zeros((2, 200, 200), dtype=np.float32)
        image[0, 20:50, 20:50] = 10000   # Bright nucleus
        image[0, 100:130, 100:130] = 1000  # Dim nucleus
        image[1, :, :] = 5000  # Uniform channel 1

        stats = compute_global_stats(image, norm_method="minmax")

        bright_patch, _ = extract_patch_with_size(image, (20, 20, 50, 50), global_stats=stats)
        dim_patch, _ = extract_patch_with_size(image, (100, 100, 130, 130), global_stats=stats)

        # Bright patch should have higher mean intensity in channel 0
        assert bright_patch[0].mean() > dim_patch[0].mean()

    def test_minmax_vs_percentile_consistency(self):
        """Test that both normalization methods produce valid [0,1] output."""
        image = np.random.rand(2, 200, 200).astype(np.float32) * 65535
        bbox = (50, 50, 80, 80)

        stats_minmax = compute_global_stats(image, norm_method="minmax")
        stats_pctl = compute_global_stats(image, norm_method="percentile")

        patch_minmax, _ = extract_patch_with_size(image, bbox, global_stats=stats_minmax)
        patch_pctl, _ = extract_patch_with_size(image, bbox, global_stats=stats_pctl)

        # Both should be in [0, 1] range
        assert patch_minmax.min() >= 0.0
        assert patch_minmax.max() <= 1.0
        assert patch_pctl.min() >= 0.0
        assert patch_pctl.max() <= 1.0

    def test_full_pipeline_directory_structure(self):
        """Test complete output directory structure for downstream compatibility."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)

            # Create mock data
            image = np.random.rand(2, 200, 200).astype(np.float32) * 65535
            stats = compute_global_stats(image, norm_method="percentile")

            # Simulate what prepare_patches.py would create
            np.savez(tmpdir / "global_stats.npz", **stats)

            # Create patches for 2 spots
            for spot_id in ["spot_0", "spot_1"]:
                bboxes = [
                    (20, 20, 40, 40),
                    (60, 60, 80, 80),
                ]
                patches = []
                sizes = []
                for bbox in bboxes:
                    patch, size = extract_patch_with_size(image, bbox, global_stats=stats)
                    patches.append(patch)
                    sizes.append(size)

                np.save(tmpdir / f"spot_{spot_id}_patches.npy", np.stack(patches))
                np.save(tmpdir / f"spot_{spot_id}_sizes.npy", np.stack(sizes))
                np.save(tmpdir / f"spot_{spot_id}_nucleus_ids.npy", np.array([1, 2], dtype=np.int64))

            # Create extraction stats (as prepare_patches.py would)
            extraction_stats = {
                "successful_patches": 4,
                "failed_patches": 0,
                "empty_spots": 0,
                "norm_method": "percentile",
            }
            import json
            with open(tmpdir / "extraction_stats.json", "w") as f:
                json.dump(extraction_stats, f)

            # Verify all required files exist
            required_files = [
                "global_stats.npz",
                "extraction_stats.json",
            ]
            for f in required_files:
                assert (tmpdir / f).exists(), f"Missing required file: {f}"

            # Verify per-spot files
            for spot_id in ["spot_0", "spot_1"]:
                assert (tmpdir / f"spot_{spot_id}_patches.npy").exists()
                assert (tmpdir / f"spot_{spot_id}_sizes.npy").exists()
                assert (tmpdir / f"spot_{spot_id}_nucleus_ids.npy").exists()

    def test_size_features_dtype(self):
        """Test that size features are float32 for PyTorch compatibility."""
        image = np.random.rand(2, 100, 100).astype(np.float32)
        stats = compute_global_stats(image, norm_method="percentile")

        _, size = extract_patch_with_size(image, (20, 20, 40, 40), global_stats=stats)

        assert size.dtype == np.float32

    def test_edge_case_small_bbox(self):
        """Test handling of very small bounding boxes."""
        image = np.random.rand(2, 100, 100).astype(np.float32)
        stats = compute_global_stats(image, norm_method="percentile")

        # 5x5 bbox (very small)
        small_bbox = (10, 10, 15, 15)

        patch, size = extract_patch_with_size(image, small_bbox, global_stats=stats)

        # Should still produce valid output
        assert patch.shape == (2, 96, 96)
        assert size.shape == (3,)
        assert not np.isnan(patch).any()
        assert not np.isnan(size).any()
