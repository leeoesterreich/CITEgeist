"""Tests for watershed cell segmentation."""
import numpy as np
import pytest
from CITEgeist.model.watershed_segmentation import watershed_from_nuclei


def test_watershed_single_nucleus():
    """Test watershed expansion from a single nucleus seed."""
    # 50x50 image with nucleus in center
    nucleus_mask = np.zeros((50, 50), dtype=np.int32)
    y, x = np.ogrid[:50, :50]
    nucleus_mask[((x - 25)**2 + (y - 25)**2) <= 25] = 1  # r=5

    # Boundary channel with ring around nucleus (simulating membrane)
    boundary_img = np.zeros((50, 50), dtype=np.float32)
    boundary_img[((x - 25)**2 + (y - 25)**2) > 64] = 0.0  # outside r=8
    boundary_img[((x - 25)**2 + (y - 25)**2) <= 64] = 1.0  # inside r=8
    boundary_img[((x - 25)**2 + (y - 25)**2) <= 25] = 0.5  # nucleus

    cell_mask = watershed_from_nuclei(nucleus_mask, boundary_img)

    # Cell mask should have same label (1) but larger area
    assert cell_mask.max() == 1
    assert (cell_mask == 1).sum() > (nucleus_mask == 1).sum()
    # Nucleus should be contained in cell
    assert np.all(cell_mask[nucleus_mask == 1] == 1)


def test_watershed_multiple_nuclei():
    """Test watershed with multiple nuclei - each gets separate cell."""
    nucleus_mask = np.zeros((100, 100), dtype=np.int32)
    y, x = np.ogrid[:100, :100]
    # Nucleus 1 at (25, 25)
    nucleus_mask[((x - 25)**2 + (y - 25)**2) <= 25] = 1
    # Nucleus 2 at (75, 75)
    nucleus_mask[((x - 75)**2 + (y - 75)**2) <= 25] = 2

    # Boundary: high signal between the two cells (around x=50)
    boundary_img = np.zeros((100, 100), dtype=np.float32)
    boundary_img[:, 45:55] = 1.0  # Vertical membrane between cells

    cell_mask = watershed_from_nuclei(nucleus_mask, boundary_img)

    # Both labels should exist
    assert set(np.unique(cell_mask)) == {1, 2}
    # Cell 1 should be mostly left half
    assert (cell_mask[:, :50] == 1).sum() > (cell_mask[:, :50] == 2).sum()
    # Cell 2 should be mostly right half
    assert (cell_mask[:, 50:] == 2).sum() > (cell_mask[:, 50:] == 1).sum()


def test_watershed_shape_mismatch():
    """Test that shape mismatch raises ValueError."""
    nucleus_mask = np.zeros((50, 50), dtype=np.int32)
    boundary_img = np.zeros((60, 60), dtype=np.float32)

    with pytest.raises(ValueError, match="Shape mismatch"):
        watershed_from_nuclei(nucleus_mask, boundary_img)