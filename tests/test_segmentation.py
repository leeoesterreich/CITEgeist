"""
Unit tests for segmentation utilities (non-Cellpose dependent pieces).
"""

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

import CITEgeist.model.morphology.segmentation as seg
from CITEgeist.model.morphology.segmentation import (
    assign_nuclei_centroids_to_spots,
    compute_spot_nuclei_counts_from_adata,
    detect_spot_diameter_pixels,
    get_resolution_image_and_scale,
    normalize_nuclei_counts_for_prior,
)


def _mock_visium_adata() -> AnnData:
    adata = AnnData(X=np.zeros((3, 2), dtype=float))
    adata.obs_names = pd.Index(["spot_0", "spot_1", "spot_2"])
    adata.obsm["spatial"] = np.array(
        [
            [100.0, 100.0],
            [200.0, 100.0],
            [300.0, 100.0],
        ]
    )
    adata.uns["spatial"] = {
        "library_1": {
            "images": {
                "hires": np.ones((50, 80, 3), dtype=np.uint8) * 127,
                "lowres": np.ones((10, 16, 3), dtype=np.uint8) * 127,
            },
            "scalefactors": {
                "tissue_hires_scalef": 0.5,
                "tissue_lowres_scalef": 0.1,
                "spot_diameter_fullres": 60.0,
            },
        }
    }
    return adata


@pytest.mark.unit
def test_assign_nuclei_centroids_to_spots_basic():
    spot_centers = np.array([[0.0, 0.0], [10.0, 0.0]])
    spot_names = pd.Index(["spot_a", "spot_b"])
    centroids = np.array(
        [
            [0.5, 0.4],   # spot_a
            [9.2, -0.3],  # spot_b
            [40.0, 40.0],  # outside radius
        ]
    )
    counts = assign_nuclei_centroids_to_spots(
        centroids_xy=centroids,
        spot_centers_xy=spot_centers,
        spot_radius_px=2.0,
        spot_names=spot_names,
    )
    assert counts.loc["spot_a"] == 1
    assert counts.loc["spot_b"] == 1
    assert int(counts.sum()) == 2


@pytest.mark.unit
def test_normalize_nuclei_counts_for_prior_clipping():
    raw = pd.Series([0, 10, 20, 100], index=["a", "b", "c", "d"], name="nuclei_count_raw")
    target = normalize_nuclei_counts_for_prior(raw, clip_min=0.25, clip_max=2.5)
    assert target.name == "nuclei_count_target"
    assert np.isclose(target.loc["c"], 20.0 / 15.0)
    assert target.min() >= 0.25
    assert target.max() <= 2.5


@pytest.mark.unit
def test_get_resolution_image_and_scale_modes():
    adata = _mock_visium_adata()
    img_low, scale_low = get_resolution_image_and_scale(adata, "lowres")
    img_hi, scale_hi = get_resolution_image_and_scale(adata, "hires")

    assert img_low.shape[:2] == (10, 16)
    assert img_hi.shape[:2] == (50, 80)
    assert np.isclose(scale_low, 0.1)
    assert np.isclose(scale_hi, 0.5)


@pytest.mark.unit
def test_get_resolution_image_and_scale_fullres_fallback():
    adata = _mock_visium_adata()
    img_full, scale_full = get_resolution_image_and_scale(adata, "fullres", max_fullres_side=500)
    # hires (50,80) with hires_scale 0.5 -> target (100,160)
    assert img_full.shape[:2] == (100, 160)
    assert np.isclose(scale_full, 1.0)


@pytest.mark.unit
def test_detect_spot_diameter_from_scalefactors():
    """Test that spot diameter is detected from scalefactors."""
    adata = _mock_visium_adata()
    diameter_px = detect_spot_diameter_pixels(adata)
    assert np.isclose(diameter_px, 60.0)  # From mock scalefactors


@pytest.mark.unit
def test_detect_spot_diameter_explicit_um():
    """Test explicit spot_diameter_um parameter."""
    adata = _mock_visium_adata()
    # Remove scalefactors to force explicit path
    del adata.uns["spatial"]["library_1"]["scalefactors"]["spot_diameter_fullres"]

    # 55 microns at 0.5 microns/pixel = 110 pixels
    diameter_px = detect_spot_diameter_pixels(
        adata,
        spot_diameter_um=55.0,
        pixel_size_um=0.5,
    )
    assert np.isclose(diameter_px, 110.0)


@pytest.mark.unit
def test_detect_spot_diameter_platform_based():
    """Test platform-based spot diameter detection."""
    adata = _mock_visium_adata()
    del adata.uns["spatial"]["library_1"]["scalefactors"]["spot_diameter_fullres"]
    adata.uns["platform"] = "visium"

    # Visium = 55 µm, at 0.2125 µm/px = ~259 px
    diameter_px = detect_spot_diameter_pixels(
        adata,
        pixel_size_um=0.2125,
    )
    assert np.isclose(diameter_px, 55.0 / 0.2125, rtol=0.01)


@pytest.mark.unit
def test_detect_spot_diameter_raises_without_info():
    """Test that error is raised when spot diameter cannot be determined."""
    adata = _mock_visium_adata()
    del adata.uns["spatial"]["library_1"]["scalefactors"]["spot_diameter_fullres"]

    with pytest.raises(ValueError, match="Could not auto-detect spot diameter"):
        detect_spot_diameter_pixels(adata)
