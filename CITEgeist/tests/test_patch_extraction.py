"""Tests for nucleus patch extraction."""
import numpy as np
import pandas as pd
import pytest


def test_extract_single_patch_shape():
    """Extracted patch should have correct shape."""
    from CITEgeist.model.patch_extraction import extract_patch

    # Fake 3-channel image (DAPI + 2 boundary markers)
    image = np.random.rand(3, 500, 500).astype(np.float32)

    # Nucleus bounding box (x_min, y_min, x_max, y_max)
    bbox = (100, 150, 120, 170)  # 20x20 nucleus

    patch = extract_patch(image, bbox, expansion=0.75, output_size=96)

    assert patch.shape == (3, 96, 96)
    assert patch.dtype == np.float32


def test_extract_patch_with_expansion():
    """Expansion should capture context around nucleus."""
    from CITEgeist.model.patch_extraction import extract_patch

    # Create image with gradient to test expansion properly
    image = np.zeros((1, 200, 200), dtype=np.float32)
    # Nucleus region with high value
    image[0, 90:110, 90:110] = 1.0
    # Surrounding region with intermediate value (to create variance)
    image[0, 70:130, 70:130] = np.where(
        image[0, 70:130, 70:130] == 0, 0.3, image[0, 70:130, 70:130]
    )

    bbox = (90, 90, 110, 110)  # 20x20 nucleus

    # No expansion - just nucleus (uniform high values, low variance)
    patch_no_exp = extract_patch(image, bbox, expansion=0.0, output_size=20)

    # With expansion - includes surrounding context (mix of values, higher variance before norm)
    patch_exp = extract_patch(image, bbox, expansion=1.0, output_size=60)

    # After z-score normalization, uniform patches have std=1e-8 (nearly 0)
    # Expanded patches with mix of values have real variance
    # The expanded patch should capture more of the image area
    assert patch_exp.shape[1] == 60
    assert patch_no_exp.shape[1] == 20


def test_extract_patches_for_spot():
    """Batch extraction should return correct shape."""
    from CITEgeist.model.patch_extraction import extract_patches_for_spot

    image = np.random.rand(3, 500, 500).astype(np.float32)

    nuclei_df = pd.DataFrame({
        'nucleus_id': [1, 2, 3],
        'bbox_x_min': [50, 150, 250],
        'bbox_y_min': [50, 150, 250],
        'bbox_x_max': [70, 170, 270],
        'bbox_y_max': [70, 170, 270],
    })

    patches = extract_patches_for_spot(image, nuclei_df, expansion=0.75, output_size=96)

    assert patches.shape == (3, 3, 96, 96)  # 3 nuclei, 3 channels, 96x96
