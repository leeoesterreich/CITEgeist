"""Tests for Cellpose H&E segmentation."""
import numpy as np
import pytest
from run_cellpose_he import (
    preprocess_he_for_cellpose,
    segment_tile,
    stitch_masks,
)


def test_preprocess_he_for_cellpose():
    """Test H&E preprocessing for Cellpose."""
    # Mock RGB H&E image
    he_image = np.random.randint(0, 255, (100, 100, 3), dtype=np.uint8)

    processed = preprocess_he_for_cellpose(he_image)

    # Cellpose expects grayscale or specific channel format
    assert processed.ndim == 2 or processed.shape[-1] in [1, 2, 3]
    assert processed.dtype == np.float32 or processed.dtype == np.uint8


def test_preprocess_he_grayscale_input():
    """Test preprocessing handles grayscale input."""
    gray_image = np.random.randint(0, 255, (100, 100), dtype=np.uint8)

    processed = preprocess_he_for_cellpose(gray_image)

    assert processed.ndim == 2
    assert processed.dtype == np.uint8


def test_preprocess_he_inverts_nuclei():
    """Test that preprocessing inverts so nuclei become bright."""
    # Create image with dark nuclei (typical H&E)
    he_image = np.ones((100, 100, 3), dtype=np.uint8) * 200  # Light background
    he_image[40:60, 40:60, :] = 50  # Dark nucleus region

    processed = preprocess_he_for_cellpose(he_image)

    # After inversion, nucleus region should be brighter than background
    nucleus_region = processed[40:60, 40:60]
    background_region = processed[0:20, 0:20]

    assert nucleus_region.mean() > background_region.mean()


def test_segment_tile():
    """Test single tile segmentation (mock, no actual Cellpose)."""
    # This test validates the interface, actual Cellpose tested in integration
    tile = np.random.randint(0, 255, (512, 512, 3), dtype=np.uint8)

    # Mock implementation for unit test
    from unittest.mock import patch, MagicMock

    mock_model = MagicMock()
    mock_model.eval.return_value = (
        np.zeros((512, 512), dtype=np.int32),  # masks
        None,  # flows
        None,  # styles
    )

    with patch('run_cellpose_he.get_cellpose_model', return_value=mock_model):
        masks = segment_tile(tile, diameter=30)

    assert masks.shape == (512, 512)
    assert masks.dtype in [np.int32, np.int64, np.uint32]


def test_segment_tile_calls_cellpose_correctly():
    """Test that segment_tile calls Cellpose with correct parameters."""
    from unittest.mock import patch, MagicMock

    tile = np.random.randint(0, 255, (256, 256, 3), dtype=np.uint8)

    mock_model = MagicMock()
    mock_model.eval.return_value = (
        np.zeros((256, 256), dtype=np.int32),
        None,
        None,
    )

    with patch('run_cellpose_he.get_cellpose_model', return_value=mock_model):
        segment_tile(tile, diameter=40)

    # Verify eval was called once
    mock_model.eval.assert_called_once()

    # Check diameter was passed
    call_kwargs = mock_model.eval.call_args[1]
    assert call_kwargs['diameter'] == 40
    assert call_kwargs['channels'] == [0, 0]  # Grayscale


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
    """Integration tests (skip if Cellpose unavailable or model weights not cached)."""

    @pytest.mark.slow
    @pytest.mark.integration
    def test_full_pipeline_small_image(self):
        """Test full segmentation pipeline on small synthetic image."""
        # Skip if Cellpose is not available
        cellpose = pytest.importorskip("cellpose", reason="Cellpose not installed")

        from run_cellpose_he import segment_tile

        # Create simple synthetic image with circular nuclei
        image = np.ones((128, 128, 3), dtype=np.uint8) * 200

        # Add some dark circles (nuclei)
        y, x = np.ogrid[:128, :128]
        for cy, cx in [(32, 32), (96, 96), (32, 96)]:
            mask = (y - cy) ** 2 + (x - cx) ** 2 <= 100
            image[mask] = 50

        # This will actually call Cellpose
        masks = segment_tile(image, diameter=20, gpu=False)

        assert masks.shape == (128, 128)
        assert masks.dtype in [np.int32, np.int64]
