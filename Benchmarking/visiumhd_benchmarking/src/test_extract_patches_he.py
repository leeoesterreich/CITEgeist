"""Tests for H&E patch extraction."""
import sys
from pathlib import Path
import numpy as np
import pytest

# Add src directory to path for imports
_src_dir = Path(__file__).parent
if str(_src_dir) not in sys.path:
    sys.path.insert(0, str(_src_dir))

from extract_patches_he import (
    extract_nucleus_patch,
    normalize_he_patch,
    extract_patches_for_spot,
)


def test_extract_nucleus_patch_basic():
    """Test basic patch extraction."""
    # Mock image and mask
    image = np.random.randint(0, 255, (500, 500, 3), dtype=np.uint8)
    mask = np.zeros((500, 500), dtype=np.int32)
    mask[100:120, 100:120] = 1  # Nucleus at (100-120, 100-120)

    patch = extract_nucleus_patch(
        image, mask, nucleus_id=1,
        output_size=224, expansion=0.75
    )

    assert patch.shape == (224, 224, 3)
    assert patch.dtype == np.uint8


def test_extract_nucleus_patch_edge():
    """Test patch extraction at image edge."""
    image = np.random.randint(0, 255, (100, 100, 3), dtype=np.uint8)
    mask = np.zeros((100, 100), dtype=np.int32)
    mask[0:10, 0:10] = 1  # Nucleus at corner

    patch = extract_nucleus_patch(
        image, mask, nucleus_id=1,
        output_size=224, expansion=0.75
    )

    assert patch.shape == (224, 224, 3)


def test_normalize_he_patch():
    """Test H&E normalization."""
    patch = np.random.randint(0, 255, (224, 224, 3), dtype=np.uint8)

    normalized = normalize_he_patch(patch)

    assert normalized.shape == (3, 224, 224)  # CHW format
    assert normalized.dtype == np.float32
    assert normalized.min() >= 0 and normalized.max() <= 1


def test_extract_patches_for_spot():
    """Test batch patch extraction for a spot."""
    image = np.random.randint(0, 255, (500, 500, 3), dtype=np.uint8)
    mask = np.zeros((500, 500), dtype=np.int32)
    mask[100:120, 100:120] = 1
    mask[200:220, 200:220] = 2
    mask[300:320, 300:320] = 3

    nucleus_ids = np.array([1, 2, 3])

    patches, valid_ids = extract_patches_for_spot(
        image, mask, nucleus_ids,
        output_size=224, expansion=0.75
    )

    assert patches.shape == (3, 3, 224, 224)  # (N, C, H, W)
    assert len(valid_ids) == 3


def test_extract_patches_for_spot_with_invalid():
    """Test batch extraction handles invalid nucleus IDs gracefully."""
    image = np.random.randint(0, 255, (500, 500, 3), dtype=np.uint8)
    mask = np.zeros((500, 500), dtype=np.int32)
    mask[100:120, 100:120] = 1
    mask[200:220, 200:220] = 2

    # Request 3 nuclei but only 2 exist
    nucleus_ids = np.array([1, 2, 99])

    patches, valid_ids = extract_patches_for_spot(
        image, mask, nucleus_ids,
        output_size=224, expansion=0.75
    )

    assert patches.shape == (2, 3, 224, 224)  # Only 2 valid
    assert len(valid_ids) == 2
    assert 99 not in valid_ids


def test_extract_patches_for_spot_empty():
    """Test batch extraction with no valid nuclei."""
    image = np.random.randint(0, 255, (500, 500, 3), dtype=np.uint8)
    mask = np.zeros((500, 500), dtype=np.int32)

    nucleus_ids = np.array([1, 2, 3])  # None exist

    patches, valid_ids = extract_patches_for_spot(
        image, mask, nucleus_ids,
        output_size=224, expansion=0.75
    )

    assert patches.shape == (0, 3, 224, 224)
    assert len(valid_ids) == 0


def test_extract_nucleus_patch_not_found():
    """Test error when nucleus ID not in mask."""
    image = np.random.randint(0, 255, (100, 100, 3), dtype=np.uint8)
    mask = np.zeros((100, 100), dtype=np.int32)

    with pytest.raises(ValueError, match="not found"):
        extract_nucleus_patch(image, mask, nucleus_id=1)
