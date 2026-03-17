"""Tests for StarDist H&E segmentation."""
import sys
from pathlib import Path
import numpy as np
import pytest

# Add src directory to path for imports
_src_dir = Path(__file__).parent
if str(_src_dir) not in sys.path:
    sys.path.insert(0, str(_src_dir))

from run_segmentation_he import (
    preprocess_he_for_stardist,
    segment_tile,
    stitch_masks,
    extract_centroids,
)


def test_preprocess_he_for_stardist():
    """Test H&E preprocessing for StarDist."""
    # Mock RGB H&E image
    he_image = np.random.randint(0, 255, (100, 100, 3), dtype=np.uint8)

    processed = preprocess_he_for_stardist(he_image)

    # StarDist expects RGB uint8
    assert processed.ndim == 3
    assert processed.shape[-1] == 3
    assert processed.dtype == np.uint8


def test_preprocess_he_grayscale_input():
    """Test preprocessing handles grayscale input by converting to RGB."""
    gray_image = np.random.randint(0, 255, (100, 100), dtype=np.uint8)

    processed = preprocess_he_for_stardist(gray_image)

    # Should be converted to 3-channel
    assert processed.ndim == 3
    assert processed.shape == (100, 100, 3)
    assert processed.dtype == np.uint8


def test_preprocess_he_preserves_rgb():
    """Test that preprocessing preserves RGB H&E images."""
    he_image = np.random.randint(0, 255, (100, 100, 3), dtype=np.uint8)

    processed = preprocess_he_for_stardist(he_image)

    # Should be unchanged for uint8 RGB input
    np.testing.assert_array_equal(processed, he_image)


def test_segment_tile():
    """Test single tile segmentation (mock, no actual StarDist)."""
    tile = np.random.randint(0, 255, (512, 512, 3), dtype=np.uint8)

    from unittest.mock import patch, MagicMock
    import pandas as pd

    mock_segmenter = MagicMock()
    mock_segmenter.segment.return_value = (
        np.zeros((512, 512), dtype=np.int32),  # masks
        pd.DataFrame(columns=['y_pixel', 'x_pixel', 'nucleus_id']),  # centroids
    )

    with patch('run_segmentation_he.get_segmenter', return_value=mock_segmenter):
        masks = segment_tile(tile)

    assert masks.shape == (512, 512)
    assert masks.dtype in [np.int32, np.int64, np.uint32]


def test_segment_tile_calls_stardist_correctly():
    """Test that segment_tile calls StarDist with correct parameters."""
    from unittest.mock import patch, MagicMock
    import pandas as pd

    tile = np.random.randint(0, 255, (256, 256, 3), dtype=np.uint8)

    mock_segmenter = MagicMock()
    mock_segmenter.segment.return_value = (
        np.zeros((256, 256), dtype=np.int32),
        pd.DataFrame(columns=['y_pixel', 'x_pixel', 'nucleus_id']),
    )

    with patch('run_segmentation_he.get_segmenter', return_value=mock_segmenter):
        segment_tile(tile, prob_thresh=0.5, nms_thresh=0.3)

    # Verify segment was called once
    mock_segmenter.segment.assert_called_once()

    # Check kwargs were passed through
    call_kwargs = mock_segmenter.segment.call_args[1]
    assert call_kwargs['prob_thresh'] == 0.5
    assert call_kwargs['nms_thresh'] == 0.3


def test_stitch_masks():
    """Test mask stitching from tiles."""
    # Two adjacent tiles with overlapping masks
    tile1 = np.array([[1, 1, 0], [1, 1, 0], [0, 0, 0]])
    tile2 = np.array([[0, 2, 2], [0, 2, 2], [0, 0, 0]])

    tiles = [(tile1, (0, 0)), (tile2, (0, 2))]  # (mask, (row_offset, col_offset))

    stitched = stitch_masks(tiles, output_shape=(3, 5), overlap=1)

    assert stitched.shape == (3, 5)
    # Should have 2 unique masks (excluding 0)
    assert len(np.unique(stitched)) - 1 <= 2


def test_stitch_masks_relabels():
    """Test that stitching relabels masks to avoid conflicts."""
    # Both tiles have label 1 - they should get relabeled to avoid conflicts
    tile1 = np.array([[1, 1], [1, 1]])
    tile2 = np.array([[1, 1], [1, 1]])

    tiles = [(tile1, (0, 0)), (tile2, (0, 2))]

    stitched = stitch_masks(tiles, output_shape=(2, 4), overlap=0)

    # Should have 2 distinct non-zero labels (relabeled from both having 1)
    unique_labels = np.unique(stitched)
    non_zero_labels = unique_labels[unique_labels > 0]
    assert len(non_zero_labels) == 2  # Should have labels 1 and 2


def test_stitch_masks_handles_overlap():
    """Test mask stitching handles overlapping regions."""
    # Create tiles with overlap region
    tile1 = np.array([
        [1, 1, 1, 0],
        [1, 1, 1, 0],
        [0, 0, 0, 0],
        [0, 0, 0, 0],
    ])
    tile2 = np.array([
        [0, 0, 2, 2],
        [0, 0, 2, 2],
        [0, 0, 0, 0],
        [0, 0, 0, 0],
    ])

    # 1 pixel overlap
    tiles = [(tile1, (0, 0)), (tile2, (0, 2))]

    stitched = stitch_masks(tiles, output_shape=(4, 6), overlap=2)

    assert stitched.shape == (4, 6)
    # Both masks should be present
    assert stitched.max() >= 2


def test_stitch_masks_empty_tiles():
    """Test stitching handles empty tiles."""
    tile1 = np.zeros((10, 10), dtype=np.int32)
    tile2 = np.array([[1, 1], [1, 1]], dtype=np.int32)

    tiles = [(tile1, (0, 0)), (tile2, (5, 5))]

    stitched = stitch_masks(tiles, output_shape=(10, 10), overlap=0)

    assert stitched.shape == (10, 10)
    assert stitched.max() == 1


def test_stitch_masks_large_output():
    """Test stitching works for larger output shapes."""
    tile_size = 64
    tiles = []

    # Create 2x2 grid of tiles
    for i in range(2):
        for j in range(2):
            mask = np.ones((tile_size, tile_size), dtype=np.int32) * (i * 2 + j + 1)
            tiles.append((mask, (i * tile_size, j * tile_size)))

    stitched = stitch_masks(tiles, output_shape=(128, 128), overlap=0)

    assert stitched.shape == (128, 128)
    # Should have 4 distinct labels (1-4) plus background
    unique = np.unique(stitched)
    assert len(unique) >= 4


class TestIntegration:
    """Integration tests (skip if StarDist unavailable or model weights not cached)."""

    @pytest.mark.slow
    @pytest.mark.integration
    def test_full_pipeline_small_image(self):
        """Test full segmentation pipeline on small synthetic image."""
        # Skip if StarDist is not available
        pytest.importorskip("stardist", reason="StarDist not installed")

        from run_segmentation_he import segment_tile

        # Create simple synthetic image with circular nuclei
        image = np.ones((128, 128, 3), dtype=np.uint8) * 200

        # Add some dark circles (nuclei)
        y, x = np.ogrid[:128, :128]
        for cy, cx in [(32, 32), (96, 96), (32, 96)]:
            circle = (y - cy) ** 2 + (x - cx) ** 2 <= 100
            image[circle] = 50

        # This will actually call StarDist
        masks = segment_tile(image, modality='he')

        assert masks.shape == (128, 128)
        assert masks.dtype in [np.int32, np.int64]


def test_extract_centroids():
    """Test centroid extraction from mask."""
    # Create mask with two nuclei at known positions
    mask = np.zeros((100, 100), dtype=np.int32)
    mask[10:20, 10:20] = 1  # Nucleus 1 centered around (15, 15)
    mask[70:80, 70:80] = 2  # Nucleus 2 centered around (75, 75)

    centroids = extract_centroids(mask)

    assert centroids.shape == (2, 2)  # 2 nuclei, (y, x) coords
    # Check centroids are approximately correct (center of each region)
    assert abs(centroids[0, 0] - 14.5) < 1  # y of nucleus 1
    assert abs(centroids[0, 1] - 14.5) < 1  # x of nucleus 1
    assert abs(centroids[1, 0] - 74.5) < 1  # y of nucleus 2
    assert abs(centroids[1, 1] - 74.5) < 1  # x of nucleus 2


def test_extract_centroids_empty():
    """Test centroid extraction with empty mask."""
    mask = np.zeros((50, 50), dtype=np.int32)

    centroids = extract_centroids(mask)

    assert centroids.shape == (0, 2)  # No nuclei
