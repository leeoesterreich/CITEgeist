# tests/benchmarking/test_patch_extractor.py
"""Tests for H&E patch extraction from WSI."""
import numpy as np
import pandas as pd
import pytest
from pathlib import Path


def test_extract_patch_correct_size():
    """Extracted patch should be 224×224×3."""
    from Benchmarking.morphology_audit.src.patch_extractor import extract_he_patch

    # Fake WSI: 1000×1000 RGB
    wsi = np.random.randint(0, 255, (1000, 1000, 3), dtype=np.uint8)
    patch = extract_he_patch(wsi, cx=500, cy=500, patch_size_px=224)
    assert patch.shape == (224, 224, 3)
    assert patch.dtype == np.uint8


def test_extract_patch_boundary_clamp():
    """Patch near boundary should be zero-padded, not crash."""
    from Benchmarking.morphology_audit.src.patch_extractor import extract_he_patch

    wsi = np.random.randint(0, 255, (200, 200, 3), dtype=np.uint8)
    # Center near edge — should not crash
    patch = extract_he_patch(wsi, cx=10, cy=10, patch_size_px=224)
    assert patch.shape == (224, 224, 3)


def test_extract_patches_batch():
    """Batch extraction should produce correct count and index."""
    from Benchmarking.morphology_audit.src.patch_extractor import extract_patches_batch

    wsi = np.random.randint(0, 255, (2000, 2000, 3), dtype=np.uint8)
    centroids = pd.DataFrame({
        "cell_x": [500, 1000, 1500],
        "cell_y": [500, 1000, 1500],
    })
    patches, valid_mask = extract_patches_batch(wsi, centroids, patch_size_px=224)
    assert len(patches) == 3
    assert all(p.shape == (224, 224, 3) for p in patches)


def test_save_patches_as_jpeg(tmp_path):
    """Patches should save as JPEG files with index CSV."""
    from Benchmarking.morphology_audit.src.patch_extractor import save_patches_jpeg

    patches = [np.random.randint(0, 255, (224, 224, 3), dtype=np.uint8) for _ in range(3)]
    cell_ids = ["cell_0", "cell_1", "cell_2"]
    spot_ids = [0, 0, 1]

    index_df = save_patches_jpeg(patches, cell_ids, spot_ids, tmp_path)
    assert len(index_df) == 3
    assert (tmp_path / "cell_0.jpg").exists()
    assert (tmp_path / "cell_1.jpg").exists()
    assert (tmp_path / "cell_2.jpg").exists()
    assert list(index_df.columns) == ["cell_id", "spot_id", "file_path"]
