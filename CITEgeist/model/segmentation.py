"""
StarDist-based nuclei segmentation utilities for Visium-style datasets.

This module provides optional image segmentation support for generating
spot-level nuclei counts that can be used as a soft abundance prior during
cell proportion optimization.
"""

from __future__ import annotations

import logging
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Optional, Tuple

import numpy as np
import pandas as pd
from scipy.ndimage import center_of_mass
from scipy.spatial import cKDTree


@dataclass
class SegmentationResult:
    """Container for segmentation and mapping outputs."""

    masks: np.ndarray
    centroids_xy: np.ndarray
    nuclei_count_raw: pd.Series
    resolution_mode: str
    image_shape: Tuple[int, int]


# Known platform spot diameters in microns
PLATFORM_SPOT_DIAMETERS_UM = {
    "visium": 55.0,
    "visium_v1": 55.0,
    "visium_v2": 55.0,
    "visium_hd": 8.0,
    "visium_hd_2um": 2.0,
    "xenium_pseudovisium": 55.0,  # Simulated Visium from Xenium
}

# Known platform spot spacings in microns (center-to-center)
PLATFORM_SPOT_SPACINGS_UM = {
    "visium": 100.0,
    "visium_v1": 100.0,
    "visium_v2": 100.0,
    "visium_hd": 8.0,  # HD bins are contiguous
    "visium_hd_2um": 2.0,
}


class StarDistSegmenter:
    """Wrapper around StarDist 2D pretrained models for nuclei segmentation.

    Lazily imports StarDist/csbdeep on first use so the module can be imported
    without those optional dependencies installed.

    Args:
        modality: One of ``"he"`` (H&E brightfield) or ``"dapi"`` (fluorescence).
    """

    MODELS = {
        "he": "2D_versatile_he",
        "dapi": "2D_versatile_fluo",
    }

    def __init__(self, modality: str = "he"):
        if modality not in self.MODELS:
            raise ValueError(
                f"Unknown modality '{modality}'. Choose from {list(self.MODELS.keys())}."
            )
        self.modality = modality

        try:
            from stardist.models import StarDist2D  # noqa: F401
        except ImportError as exc:
            raise ImportError(
                "StarDist is required for segmentation. Install with "
                "`pip install stardist csbdeep` or `pip install \"citegeist[imaging]\"`."
            ) from exc

        model_name = self.MODELS[modality]
        logging.info("Loading StarDist model '%s' for modality '%s'", model_name, modality)
        self._model = StarDist2D.from_pretrained(model_name)

    def segment(
        self,
        image: np.ndarray,
        **kwargs,
    ) -> Tuple[np.ndarray, pd.DataFrame]:
        """Run StarDist prediction on an image.

        Args:
            image: Input image. For ``"he"`` modality this should be an RGB
                uint8 image; for ``"dapi"`` a single-channel grayscale or
                uint8 image.
            **kwargs: Extra keyword arguments forwarded to
                ``StarDist2D.predict_instances`` (e.g. ``prob_thresh``,
                ``nms_thresh``, ``scale``).

        Returns:
            masks: ``int32`` label image where each nucleus has a unique ID.
            centroids_df: DataFrame with columns ``[y_pixel, x_pixel, nucleus_id]``
                containing the centroid of each detected nucleus.
        """
        from csbdeep.utils import normalize as csb_normalize

        # Normalize image to [0, 1] as expected by StarDist
        if self.modality == "he":
            # H&E: normalize per-channel
            img_norm = csb_normalize(image, 1, 99.8, axis=(0, 1))
        else:
            # Fluorescence: normalize globally
            if image.ndim == 3:
                # Take first channel if multi-channel
                image = image[..., 0]
            img_norm = csb_normalize(image, 1, 99.8)

        masks, details = self._model.predict_instances(img_norm, **kwargs)
        masks = np.asarray(masks, dtype=np.int32)

        max_label = int(masks.max())
        if max_label <= 0:
            centroids_df = pd.DataFrame(
                columns=["y_pixel", "x_pixel", "nucleus_id"]
            )
        else:
            labels = np.arange(1, max_label + 1, dtype=np.int32)
            com = center_of_mass(
                np.ones_like(masks, dtype=np.float32),
                labels=masks,
                index=labels,
            )
            centroids_df = pd.DataFrame(
                {
                    "y_pixel": [float(c[0]) for c in com],
                    "x_pixel": [float(c[1]) for c in com],
                    "nucleus_id": labels,
                }
            )

        logging.info(
            "StarDist segmented %d nuclei from %s image (%s)",
            len(centroids_df),
            self.modality,
            image.shape[:2],
        )
        return masks, centroids_df


def run_nuclei_segmentation(
    image: np.ndarray,
    modality: str = "he",
    **kwargs,
) -> Tuple[np.ndarray, pd.DataFrame]:
    """Convenience function: segment nuclei with StarDist.

    Args:
        image: Input image (RGB uint8 for H&E, grayscale for DAPI).
        modality: ``"he"`` or ``"dapi"``.
        **kwargs: Forwarded to ``StarDistSegmenter.segment()``.

    Returns:
        masks: ``int32`` label image.
        centroids_df: DataFrame with ``[y_pixel, x_pixel, nucleus_id]``.
    """
    segmenter = StarDistSegmenter(modality=modality)
    return segmenter.segment(image, **kwargs)


def run_cellpose_nuclei_segmentation(
    image_rgb_uint8: np.ndarray,
    use_gpu: bool = False,
    diameter: Optional[float] = None,
    flow_threshold: float = 0.4,
    cellprob_threshold: float = 0.0,
    model=None,
    model_type: str = "nuclei",
) -> Tuple[np.ndarray, np.ndarray]:
    """Deprecated: use ``run_nuclei_segmentation(image, modality='dapi')`` instead.

    This wrapper calls StarDist under the hood and converts the output to the
    legacy ``(masks, centroids_xy)`` format expected by old callers.
    """
    warnings.warn(
        "run_cellpose_nuclei_segmentation() is deprecated. "
        "Use run_nuclei_segmentation(image, modality='dapi') instead.",
        DeprecationWarning,
        stacklevel=2,
    )
    masks, centroids_df = run_nuclei_segmentation(image_rgb_uint8, modality="dapi")

    # Convert centroids_df to legacy (x, y) array format
    if len(centroids_df) == 0:
        centroids_xy = np.zeros((0, 2), dtype=np.float64)
    else:
        centroids_xy = np.column_stack([
            centroids_df["x_pixel"].to_numpy(dtype=np.float64),
            centroids_df["y_pixel"].to_numpy(dtype=np.float64),
        ])

    return masks, centroids_xy


def detect_spot_diameter_pixels(
    adata,
    pixel_size_um: Optional[float] = None,
    spot_diameter_um: Optional[float] = None,
) -> float:
    """
    Auto-detect spot diameter in pixels for nuclei-to-spot assignment.

    Detection hierarchy:
    1. If spot_diameter_um is provided explicitly, use it
    2. If scalefactors contains spot_diameter_fullres, use it directly
    3. If platform metadata is available, use known platform geometry
    4. Raise error with guidance if auto-detection fails

    Args:
        adata: AnnData with spatial metadata
        pixel_size_um: Microns per pixel in the coordinate frame. Required if
            spot_diameter_um is provided or platform-based detection is used.
        spot_diameter_um: Explicit spot diameter in microns (overrides auto-detection)

    Returns:
        Spot diameter in pixels (in the fullres/coordinate frame)

    Raises:
        ValueError: If spot diameter cannot be determined
    """
    lib = _get_first_library_payload(adata)
    scalefactors = lib.get("scalefactors", {})

    # Priority 1: Explicit diameter in microns
    if spot_diameter_um is not None:
        if pixel_size_um is None or pixel_size_um <= 0:
            raise ValueError(
                "pixel_size_um must be provided when specifying spot_diameter_um. "
                "This is the microns-per-pixel scale of your coordinate system."
            )
        diameter_px = spot_diameter_um / pixel_size_um
        logging.info(
            "Using explicit spot diameter: %.1f um = %.1f pixels (pixel_size=%.4f um/px)",
            spot_diameter_um, diameter_px, pixel_size_um
        )
        return diameter_px

    # Priority 2: Existing scalefactors (standard Visium from SpaceRanger)
    if "spot_diameter_fullres" in scalefactors:
        diameter_px = float(scalefactors["spot_diameter_fullres"])
        if np.isfinite(diameter_px) and diameter_px > 0:
            logging.info(
                "Using spot_diameter_fullres from scalefactors: %.1f pixels", diameter_px
            )
            return diameter_px

    # Priority 3: Platform metadata
    platform = None
    if "platform" in adata.uns:
        platform = str(adata.uns["platform"]).lower().strip()
    elif "spatial" in adata.uns:
        # Check library key for platform hints
        lib_key = next(iter(adata.uns["spatial"].keys()), "").lower()
        for known_platform in PLATFORM_SPOT_DIAMETERS_UM:
            if known_platform.replace("_", "") in lib_key.replace("_", ""):
                platform = known_platform
                break

    if platform and platform in PLATFORM_SPOT_DIAMETERS_UM:
        if pixel_size_um is None or pixel_size_um <= 0:
            raise ValueError(
                f"Detected platform '{platform}' but pixel_size_um is required to convert "
                f"spot diameter from microns to pixels. Please provide pixel_size_um "
                f"(microns per pixel in your coordinate frame)."
            )
        diameter_um = PLATFORM_SPOT_DIAMETERS_UM[platform]
        diameter_px = diameter_um / pixel_size_um
        logging.info(
            "Detected platform '%s': spot diameter = %.1f um = %.1f pixels",
            platform, diameter_um, diameter_px
        )
        return diameter_px

    # Auto-detection failed - provide helpful error
    raise ValueError(
        "Could not auto-detect spot diameter. Please provide one of:\n"
        "  1. spot_diameter_um: Explicit spot diameter in microns (e.g., 55.0 for Visium)\n"
        "  2. Ensure adata.uns['spatial'][lib]['scalefactors']['spot_diameter_fullres'] exists\n"
        "  3. Set adata.uns['platform'] to a known platform: "
        f"{list(PLATFORM_SPOT_DIAMETERS_UM.keys())}\n\n"
        "Common spot diameters:\n"
        "  - Visium: 55 um diameter, 100 um center-to-center spacing\n"
        "  - Visium HD: 8 um bins (contiguous)\n"
        "  - Xenium pseudo-Visium: typically 55 um to match Visium geometry"
    )


def _require_spatial_uns(adata) -> Dict:
    if "spatial" not in adata.uns:
        raise ValueError("AnnData missing `uns['spatial']`; cannot access Visium images/scalefactors.")
    payload = adata.uns["spatial"]
    if not isinstance(payload, dict) or len(payload) == 0:
        raise ValueError("AnnData `uns['spatial']` is empty; cannot run image-based segmentation.")
    return payload


def _get_first_library_payload(adata) -> Dict:
    spatial_payload = _require_spatial_uns(adata)
    first_key = next(iter(spatial_payload.keys()))
    return spatial_payload[first_key]


def _prepare_rgb_uint8(image: np.ndarray) -> np.ndarray:
    """Ensure image is uint8 RGB for consistent input to segmentation models."""
    if image.ndim == 2:
        image = np.stack([image, image, image], axis=-1)
    if image.ndim != 3 or image.shape[-1] not in (3, 4):
        raise ValueError(f"Expected 2D grayscale or 3D RGB(A) image, got shape {image.shape}.")
    if image.shape[-1] == 4:
        image = image[..., :3]

    if image.dtype == np.uint8:
        return image

    image = np.asarray(image, dtype=np.float32)
    vmin = float(np.nanmin(image))
    vmax = float(np.nanmax(image))
    if not np.isfinite(vmin) or not np.isfinite(vmax):
        raise ValueError("Input image contains non-finite values.")
    if vmax <= vmin:
        return np.zeros_like(image, dtype=np.uint8)
    scaled = (image - vmin) / (vmax - vmin)
    return (255.0 * np.clip(scaled, 0.0, 1.0)).astype(np.uint8)


def get_resolution_image_and_scale(
    adata,
    resolution_mode: str = "hires",
    max_fullres_side: int = 9000,
) -> Tuple[np.ndarray, float]:
    """
    Return segmentation image and fullres->image coordinate scale.

    Returns:
        image_rgb_uint8: image used for segmentation
        fullres_to_image_scale: multiply fullres pixel coordinates by this
            value to map spot centers/radii into the returned image frame.
    """
    if resolution_mode not in {"lowres", "hires", "fullres"}:
        raise ValueError("resolution_mode must be one of {'lowres', 'hires', 'fullres'}.")

    lib = _get_first_library_payload(adata)
    images = lib.get("images", {})
    scalefactors = lib.get("scalefactors", {})

    if resolution_mode == "lowres":
        if "lowres" not in images:
            raise ValueError("No lowres image found in `uns['spatial'][lib]['images']['lowres']`.")
        scale = float(scalefactors.get("tissue_lowres_scalef", np.nan))
        if not np.isfinite(scale):
            raise ValueError("Missing/invalid `tissue_lowres_scalef` in scalefactors.")
        return _prepare_rgb_uint8(images["lowres"]), scale

    if resolution_mode == "hires":
        if "hires" not in images:
            raise ValueError("No hires image found in `uns['spatial'][lib]['images']['hires']`.")
        scale = float(scalefactors.get("tissue_hires_scalef", np.nan))
        if not np.isfinite(scale):
            raise ValueError("Missing/invalid `tissue_hires_scalef` in scalefactors.")
        return _prepare_rgb_uint8(images["hires"]), scale

    # fullres mode
    if "fullres" in images:
        return _prepare_rgb_uint8(images["fullres"]), 1.0

    # Fallback: approximate fullres by upscaling hires to fullres coordinate frame.
    if "hires" not in images:
        raise ValueError("Cannot construct fullres fallback: hires image not available.")

    hires = _prepare_rgb_uint8(images["hires"])
    hires_scale = float(scalefactors.get("tissue_hires_scalef", np.nan))
    if not np.isfinite(hires_scale) or hires_scale <= 0:
        raise ValueError("Missing/invalid `tissue_hires_scalef`; cannot build fullres fallback.")

    # Target dimensions in fullres frame.
    target_h = int(round(hires.shape[0] / hires_scale))
    target_w = int(round(hires.shape[1] / hires_scale))
    target_max = max(target_h, target_w)

    # Guard against extreme upscaling. Preserve coordinate frame with downscaling factor.
    if target_max > max_fullres_side:
        shrink = max_fullres_side / float(target_max)
        target_h = max(1, int(round(target_h * shrink)))
        target_w = max(1, int(round(target_w * shrink)))
        # Effective scale in image frame after shrink.
        effective_scale = shrink
    else:
        effective_scale = 1.0

    # Lightweight nearest-neighbor upscale using index mapping (avoids extra deps).
    y_idx = np.floor(np.linspace(0, hires.shape[0] - 1, target_h)).astype(int)
    x_idx = np.floor(np.linspace(0, hires.shape[1] - 1, target_w)).astype(int)
    fullres_approx = hires[y_idx][:, x_idx]
    logging.warning(
        "fullres image not present; using hires-upscaled fallback (%s -> %s).",
        hires.shape[:2],
        fullres_approx.shape[:2],
    )
    return fullres_approx, effective_scale


def assign_nuclei_centroids_to_spots(
    centroids_xy: np.ndarray,
    spot_centers_xy: np.ndarray,
    spot_radius_px: float,
    spot_names: pd.Index,
) -> pd.Series:
    """
    Assign each nucleus centroid to nearest spot within spot radius.
    """
    if spot_centers_xy.ndim != 2 or spot_centers_xy.shape[1] != 2:
        raise ValueError("spot_centers_xy must have shape (N, 2).")
    if len(spot_names) != spot_centers_xy.shape[0]:
        raise ValueError("spot_names length must match spot_centers_xy rows.")
    if spot_radius_px <= 0:
        raise ValueError("spot_radius_px must be positive.")

    counts = np.zeros(spot_centers_xy.shape[0], dtype=np.int64)
    if centroids_xy.size == 0:
        return pd.Series(counts, index=spot_names, name="nuclei_count_raw")

    # Filter out spots with NaN/Inf coordinates before building KDTree
    finite_mask = np.isfinite(spot_centers_xy).all(axis=1)
    if not finite_mask.all():
        n_bad = (~finite_mask).sum()
        logging.warning("Filtered %d spots with non-finite spatial coordinates", n_bad)
    finite_centers = spot_centers_xy[finite_mask]
    finite_indices = np.where(finite_mask)[0]

    if finite_centers.shape[0] == 0:
        logging.warning("No spots with finite coordinates; returning zero counts")
        return pd.Series(counts, index=spot_names, name="nuclei_count_raw")

    # Also filter NaN centroids
    finite_cent_mask = np.isfinite(centroids_xy).all(axis=1)
    clean_centroids = centroids_xy[finite_cent_mask]

    tree = cKDTree(finite_centers)
    dists, idxs = tree.query(clean_centroids, k=1, workers=-1)
    valid = dists <= float(spot_radius_px)
    if np.any(valid):
        # Map back from finite subset indices to original indices
        original_idxs = finite_indices[idxs[valid]]
        np.add.at(counts, original_idxs, 1)
    return pd.Series(counts, index=spot_names, name="nuclei_count_raw")


def normalize_nuclei_counts_for_prior(
    nuclei_count_raw: pd.Series,
    clip_min: float = 0.25,
    clip_max: float = 2.5,
) -> pd.Series:
    """
    Normalize raw nuclei counts into a soft abundance target around 1.0.
    """
    values = nuclei_count_raw.astype(float).to_numpy()
    median = float(np.median(values))
    if not np.isfinite(median) or median <= 1e-8:
        median = 1.0
    normalized = values / median
    normalized = np.clip(normalized, clip_min, clip_max)
    return pd.Series(normalized, index=nuclei_count_raw.index, name="nuclei_count_target")


def compute_spot_nuclei_counts_from_adata(
    adata,
    resolution_mode: str = "hires",
    max_fullres_side: int = 9000,
    spot_diameter_um: Optional[float] = None,
    pixel_size_um: Optional[float] = None,
    modality: str = "he",
    **stardist_kwargs,
) -> SegmentationResult:
    """
    End-to-end nuclei segmentation + centroid assignment for an AnnData Visium object.

    Uses StarDist to segment nuclei from the histology image, then assigns
    each detected nucleus to the nearest Visium spot within the spot radius.

    Spot diameter is determined by (in order of priority):
    1. spot_diameter_um parameter (explicit, requires pixel_size_um)
    2. scalefactors['spot_diameter_fullres'] from adata.uns['spatial']
    3. Platform-based detection from adata.uns['platform'] (requires pixel_size_um)

    Args:
        adata: AnnData with spatial coordinates and histology images.
        resolution_mode: Image resolution to use ('lowres', 'hires', 'fullres').
        max_fullres_side: Max image dimension for fullres fallback.
        spot_diameter_um: Explicit spot diameter in microns (e.g., 55.0 for Visium).
        pixel_size_um: Microns per pixel in coordinate frame (required for um-based params).
        modality: StarDist modality - ``"he"`` for H&E or ``"dapi"`` for fluorescence.
        **stardist_kwargs: Extra arguments forwarded to ``StarDistSegmenter.segment()``
            (e.g. ``prob_thresh``, ``nms_thresh``, ``scale``).

    Returns:
        SegmentationResult with nuclei counts per spot.
    """
    if "spatial" not in adata.obsm:
        raise ValueError("AnnData missing `obsm['spatial']`; cannot map nuclei to spots.")

    image, fullres_to_image_scale = get_resolution_image_and_scale(
        adata=adata,
        resolution_mode=resolution_mode,
        max_fullres_side=max_fullres_side,
    )

    # Auto-detect spot diameter using the detection hierarchy
    spot_diam_fullres = detect_spot_diameter_pixels(
        adata=adata,
        pixel_size_um=pixel_size_um,
        spot_diameter_um=spot_diameter_um,
    )

    spot_centers_fullres = np.asarray(adata.obsm["spatial"], dtype=np.float64)

    # StarDist processes the full image in one call (no patch mode needed)
    masks, centroids_df = run_nuclei_segmentation(
        image=image, modality=modality, **stardist_kwargs
    )

    # Convert centroids_df (y_pixel, x_pixel) to xy array for spot assignment
    if len(centroids_df) == 0:
        centroids_xy = np.zeros((0, 2), dtype=np.float64)
    else:
        centroids_xy = np.column_stack([
            centroids_df["x_pixel"].to_numpy(dtype=np.float64),
            centroids_df["y_pixel"].to_numpy(dtype=np.float64),
        ])

    spot_centers_xy = spot_centers_fullres * float(fullres_to_image_scale)
    spot_radius_px = (spot_diam_fullres * float(fullres_to_image_scale)) / 2.0
    raw_counts = assign_nuclei_centroids_to_spots(
        centroids_xy=centroids_xy,
        spot_centers_xy=spot_centers_xy,
        spot_radius_px=spot_radius_px,
        spot_names=adata.obs_names,
    )

    return SegmentationResult(
        masks=masks,
        centroids_xy=centroids_xy,
        nuclei_count_raw=raw_counts,
        resolution_mode=resolution_mode,
        image_shape=(int(image.shape[0]), int(image.shape[1])),
    )


def save_segmentation_artifacts(
    output_folder: str,
    sample_name: str,
    result: SegmentationResult,
    save_masks: bool = True,
) -> Dict[str, str]:
    """Persist segmentation artifacts and return output paths."""
    out = Path(output_folder)
    out.mkdir(parents=True, exist_ok=True)

    counts_path = out / f"{sample_name}_nuclei_counts_{result.resolution_mode}.csv"
    result.nuclei_count_raw.to_csv(counts_path, header=True)

    outputs: Dict[str, str] = {"nuclei_counts_csv": str(counts_path)}
    if save_masks:
        masks_path = out / f"{sample_name}_stardist_masks_{result.resolution_mode}.npy"
        np.save(masks_path, result.masks)
        outputs["stardist_masks_npy"] = str(masks_path)
    return outputs
