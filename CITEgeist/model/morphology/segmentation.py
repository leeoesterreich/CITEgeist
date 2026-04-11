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
from scipy.spatial import cKDTree


@dataclass
class SegmentationQC:
    """Quality-control metrics for a segmentation run."""

    n_raw: int
    n_after_area_filter: int
    n_removed_small: int
    n_removed_large: int
    pixel_size_um: Optional[float]
    min_area_um2: float
    max_area_um2: float
    area_median_um2: float
    area_p25_um2: float
    area_p75_um2: float
    density_per_mm2: float  # nuclei / tissue area in mm²
    density_flag: str  # "ok", "high", "low"

    def to_dict(self) -> dict:
        return {k: v for k, v in self.__dict__.items()}


# Expected nuclei density range (nuclei per mm²) for human tissue.
# Sparse stroma ~500/mm², dense lymphoid tissue ~5000/mm².
EXPECTED_DENSITY_RANGE = (500.0, 5000.0)
# Density flag thresholds (multiples of the upper/lower bound)
DENSITY_HIGH_FACTOR = 2.0
DENSITY_LOW_FACTOR = 0.5


@dataclass
class SegmentationResult:
    """Container for segmentation and mapping outputs."""

    masks: np.ndarray
    centroids_xy: np.ndarray
    nuclei_count_raw: pd.Series
    resolution_mode: str
    image_shape: Tuple[int, int]
    nucleus_spot_map: Optional[pd.DataFrame] = None
    qc: Optional[SegmentationQC] = None


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
            raise ValueError(f"Unknown modality '{modality}'. Choose from {list(self.MODELS.keys())}.")
        self.modality = modality

        try:
            from stardist.models import StarDist2D  # noqa: F401  # pylint: disable=import-outside-toplevel
        except ImportError as exc:
            raise ImportError(
                "StarDist is required for segmentation. Install with "
                '`pip install stardist csbdeep` or `pip install "citegeist[imaging]"`.'
            ) from exc

        model_name = self.MODELS[modality]
        logging.info("Loading StarDist model '%s' for modality '%s'", model_name, modality)
        self._model = StarDist2D.from_pretrained(model_name)

    # StarDist's 2D_versatile_fluo was trained at 0.25 um/pixel
    TRAINING_PIXEL_SIZE_UM = 0.25
    # Default prob_thresh — conservative default for quantitative work.
    # StarDist's own default (0.479) is optimized for visual recall, not
    # counting accuracy.  0.6 removes low-confidence debris.
    DEFAULT_PROB_THRESH = 0.6
    # Nuclear area bounds in µm² (biology constraint, tissue-agnostic).
    # Smallest human nucleus ~5 µm diameter (lymphocyte) → ~20 µm².
    # Largest plausible single nucleus ~25 µm diameter → ~500 µm².
    DEFAULT_MIN_AREA_UM2 = 20.0
    DEFAULT_MAX_AREA_UM2 = 500.0

    def segment(
        self,
        image: np.ndarray,
        pixel_size_um: Optional[float] = None,
        min_area_um2: Optional[float] = None,
        max_area_um2: Optional[float] = None,
        **kwargs,
    ) -> Tuple[np.ndarray, pd.DataFrame]:
        """Run StarDist prediction on an image with area-based post-filtering.

        Args:
            image: Input image. For ``"he"`` modality this should be an RGB
                uint8 image; for ``"dapi"`` a single-channel grayscale or
                uint8 image.
            pixel_size_um: Microns per pixel. When provided and ``scale`` is
                not in *kwargs*, automatically sets
                ``scale = pixel_size_um / 0.25`` to correct for the model's
                training resolution.  Also enables area filtering and QC.
            min_area_um2: Minimum nucleus area in µm² (default 20).
                Set to 0 to disable lower bound filtering.
            max_area_um2: Maximum nucleus area in µm² (default 500).
                Set to ``float('inf')`` to disable upper bound filtering.
            **kwargs: Extra keyword arguments forwarded to
                ``StarDist2D.predict_instances`` (e.g. ``prob_thresh``,
                ``nms_thresh``, ``scale``).

        Returns:
            masks: ``int32`` label image where each nucleus has a unique ID.
                Labels removed by area filtering are zeroed out.
            centroids_df: DataFrame with columns ``[y_pixel, x_pixel, nucleus_id]``
                containing the centroid of each detected nucleus (after filtering).
        """
        from csbdeep.utils import normalize as csb_normalize  # pylint: disable=import-outside-toplevel

        # Auto-compute scale from pixel size if not explicitly provided
        if pixel_size_um is not None and "scale" not in kwargs:
            kwargs["scale"] = pixel_size_um / self.TRAINING_PIXEL_SIZE_UM
            logging.info(
                "Auto scale=%.3f from pixel_size_um=%.4f",
                kwargs["scale"],
                pixel_size_um,
            )

        # Apply default prob_thresh if not explicitly provided
        if "prob_thresh" not in kwargs:
            kwargs["prob_thresh"] = self.DEFAULT_PROB_THRESH

        # Normalize image to [0, 1] as expected by StarDist
        if self.modality == "he":
            img_norm = csb_normalize(image, 1, 99.8, axis=(0, 1))
        else:
            if image.ndim == 3:
                image = image[..., 0]
            img_norm = csb_normalize(image, 1, 99.8)

        masks, details = self._model.predict_instances(img_norm, **kwargs)
        masks = np.asarray(masks, dtype=np.int32)

        # Use StarDist's native centroid output instead of center_of_mass
        # (center_of_mass OOMs on large images with >100K labels)
        points = details["points"]  # (N, 2) in (row, col) order
        n_raw = len(points)

        if n_raw == 0:
            centroids_df = pd.DataFrame(columns=["y_pixel", "x_pixel", "nucleus_id"])
            logging.info("StarDist segmented 0 nuclei from %s image", self.modality)
            return masks, centroids_df

        # --- Area-based post-filtering (requires pixel_size_um) ---
        if pixel_size_um is not None:
            px_area_per_um2 = pixel_size_um**2
            _min = self.DEFAULT_MIN_AREA_UM2 if min_area_um2 is None else min_area_um2
            _max = self.DEFAULT_MAX_AREA_UM2 if max_area_um2 is None else max_area_um2
            min_area_px = _min / px_area_per_um2
            max_area_px = _max / px_area_per_um2

            # Compute area per label via bincount (memory-efficient)
            areas_px = np.bincount(masks.ravel())  # index 0 = background
            # labels are 1..n_raw
            label_areas = areas_px[1 : n_raw + 1] if len(areas_px) > n_raw else areas_px[1:]

            keep_mask = (label_areas >= min_area_px) & (label_areas <= max_area_px)
            n_small = int(np.sum(label_areas < min_area_px))
            n_large = int(np.sum(label_areas > max_area_px))
            n_kept = int(keep_mask.sum())

            if n_kept < n_raw:
                # Zero out removed labels in the mask
                remove_labels = np.where(~keep_mask)[0] + 1  # 1-based
                if len(remove_labels) > 0:
                    remove_set = np.zeros(n_raw + 1, dtype=bool)
                    remove_set[remove_labels] = True
                    masks[remove_set[masks]] = 0

                # Filter points
                points = points[keep_mask]
                label_areas = label_areas[keep_mask]

            areas_um2 = label_areas * px_area_per_um2

            # QC: density estimation
            image_area_mm2 = (image.shape[0] * pixel_size_um * image.shape[1] * pixel_size_um) / 1e6
            density = n_kept / image_area_mm2 if image_area_mm2 > 0 else 0.0
            lo, hi = EXPECTED_DENSITY_RANGE
            if density > hi * DENSITY_HIGH_FACTOR:
                flag = "high"
            elif density < lo * DENSITY_LOW_FACTOR:
                flag = "low"
            else:
                flag = "ok"

            # Store QC on the class so callers can retrieve it
            self._last_qc = SegmentationQC(
                n_raw=n_raw,
                n_after_area_filter=n_kept,
                n_removed_small=n_small,
                n_removed_large=n_large,
                pixel_size_um=pixel_size_um,
                min_area_um2=_min,
                max_area_um2=_max,
                area_median_um2=float(np.median(areas_um2)) if len(areas_um2) > 0 else 0.0,
                area_p25_um2=float(np.percentile(areas_um2, 25)) if len(areas_um2) > 0 else 0.0,
                area_p75_um2=float(np.percentile(areas_um2, 75)) if len(areas_um2) > 0 else 0.0,
                density_per_mm2=density,
                density_flag=flag,
            )

            logging.info(
                "Area filter: %d → %d nuclei (-%d small < %.0f µm², -%d large > %.0f µm²)",
                n_raw,
                n_kept,
                n_small,
                _min,
                n_large,
                _max,
            )
            if flag != "ok":
                logging.warning(
                    "Nuclei density %.0f/mm² flagged as %s (expected %.0f–%.0f/mm²)",
                    density,
                    flag,
                    lo,
                    hi,
                )
        else:
            self._last_qc = None

        n_nuclei = len(points)
        centroids_df = pd.DataFrame(
            {
                "y_pixel": points[:, 0].astype(np.float64),
                "x_pixel": points[:, 1].astype(np.float64),
                "nucleus_id": np.arange(1, n_nuclei + 1, dtype=np.int32),
            }
        )

        logging.info(
            "StarDist segmented %d nuclei from %s image (%s)",
            n_nuclei,
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
    masks, centroids_df = segmenter.segment(image, **kwargs)
    # Stash QC on the function for callers that need it
    run_nuclei_segmentation._last_qc = getattr(segmenter, "_last_qc", None)  # type: ignore[attr-defined]
    return masks, centroids_df


def run_cellpose_nuclei_segmentation(
    image_rgb_uint8: np.ndarray,
    *,
    _use_gpu: bool = False,
    _diameter: Optional[float] = None,
    _flow_threshold: float = 0.4,
    _cellprob_threshold: float = 0.0,
    _model=None,
    _model_type: str = "nuclei",
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
        centroids_xy = np.column_stack(
            [
                centroids_df["x_pixel"].to_numpy(dtype=np.float64),
                centroids_df["y_pixel"].to_numpy(dtype=np.float64),
            ]
        )

    return masks, centroids_xy


def estimate_pixel_size_um(
    adata,
    spot_diameter_um: float = 55.0,
) -> Optional[float]:
    """Estimate microns-per-pixel from Visium scalefactors.

    For standard Visium data, SpaceRanger provides ``spot_diameter_fullres``
    (the spot diameter in fullres image pixels).  Since the physical spot
    diameter is known (55 µm for Visium), pixel size follows:

        pixel_size_um = spot_diameter_um / spot_diameter_fullres

    Args:
        adata: AnnData with ``uns['spatial'][lib]['scalefactors']``.
        spot_diameter_um: Physical spot diameter in microns (default 55.0
            for standard Visium).

    Returns:
        Estimated pixel size in µm/pixel, or ``None`` if scalefactors are
        missing or incomplete.
    """
    try:
        lib = _get_first_library_payload(adata)
    except ValueError:
        return None

    scalefactors = lib.get("scalefactors", {})
    diam_px = scalefactors.get("spot_diameter_fullres")
    if diam_px is None or not np.isfinite(float(diam_px)) or float(diam_px) <= 0:
        return None

    pixel_size = spot_diameter_um / float(diam_px)
    logging.info(
        "Estimated pixel_size_um=%.4f from spot_diameter_fullres=%.1f px " "and known spot diameter=%.1f µm",
        pixel_size,
        float(diam_px),
        spot_diameter_um,
    )
    return pixel_size


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
            spot_diameter_um,
            diameter_px,
            pixel_size_um,
        )
        return diameter_px

    # Priority 2: Existing scalefactors (standard Visium from SpaceRanger)
    if "spot_diameter_fullres" in scalefactors:
        diameter_px = float(scalefactors["spot_diameter_fullres"])
        if np.isfinite(diameter_px) and diameter_px > 0:
            logging.info("Using spot_diameter_fullres from scalefactors: %.1f pixels", diameter_px)
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
            "Detected platform '%s': spot diameter = %.1f um = %.1f pixels", platform, diameter_um, diameter_px
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
    *,
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

    # Auto-estimate pixel_size_um from scalefactors if not provided.
    # This enables the 3-layer defense (scale correction + area filtering)
    # without requiring the caller to know the image resolution.
    if pixel_size_um is None:
        pixel_size_um = estimate_pixel_size_um(adata)
        if pixel_size_um is not None:
            logging.info("Auto-estimated pixel_size_um=%.4f from scalefactors", pixel_size_um)

    image, fullres_to_image_scale = get_resolution_image_and_scale(
        adata=adata,
        resolution_mode=resolution_mode,
        max_fullres_side=max_fullres_side,
    )

    # Adjust pixel_size_um for the image resolution we're actually segmenting.
    # scalefactors give pixel size in the fullres frame; if we're segmenting
    # the hires image, pixels are larger by 1/scale.
    pixel_size_image = None
    if pixel_size_um is not None and fullres_to_image_scale > 0:
        pixel_size_image = pixel_size_um / fullres_to_image_scale
        logging.info(
            "Pixel size in %s image frame: %.4f µm/px (fullres: %.4f, scale: %.4f)",
            resolution_mode,
            pixel_size_image,
            pixel_size_um,
            fullres_to_image_scale,
        )

    # Auto-detect spot diameter using the detection hierarchy
    spot_diam_fullres = detect_spot_diameter_pixels(
        adata=adata,
        pixel_size_um=pixel_size_um,
        spot_diameter_um=spot_diameter_um,
    )

    spot_centers_fullres = np.asarray(adata.obsm["spatial"], dtype=np.float64)

    # StarDist processes the full image in one call (no patch mode needed).
    # Pass pixel_size in the image frame so StarDist auto-computes scale
    # correction and area filtering works at the correct resolution.
    masks, centroids_df = run_nuclei_segmentation(
        image=image, modality=modality, pixel_size_um=pixel_size_image, **stardist_kwargs
    )

    # Convert centroids_df (y_pixel, x_pixel) to xy array for spot assignment
    if len(centroids_df) == 0:
        centroids_xy = np.zeros((0, 2), dtype=np.float64)
    else:
        centroids_xy = np.column_stack(
            [
                centroids_df["x_pixel"].to_numpy(dtype=np.float64),
                centroids_df["y_pixel"].to_numpy(dtype=np.float64),
            ]
        )

    spot_centers_xy = spot_centers_fullres * float(fullres_to_image_scale)
    spot_radius_px = (spot_diam_fullres * float(fullres_to_image_scale)) / 2.0
    raw_counts = assign_nuclei_centroids_to_spots(
        centroids_xy=centroids_xy,
        spot_centers_xy=spot_centers_xy,
        spot_radius_px=spot_radius_px,
        spot_names=adata.obs_names,
    )

    # Build nucleus-spot mapping DataFrame
    nucleus_spot_map = None
    if centroids_xy.size > 0:
        spot_names_arr = np.asarray(adata.obs_names)
        finite_mask = np.isfinite(spot_centers_xy).all(axis=1)
        finite_centers = spot_centers_xy[finite_mask]
        finite_indices = np.where(finite_mask)[0]
        finite_cent_mask = np.isfinite(centroids_xy).all(axis=1)
        clean_centroids = centroids_xy[finite_cent_mask]
        clean_indices = np.where(finite_cent_mask)[0]

        if finite_centers.shape[0] > 0 and clean_centroids.shape[0] > 0:
            tree = cKDTree(finite_centers)
            dists, idxs = tree.query(clean_centroids, k=1, workers=-1)
            valid = dists <= float(spot_radius_px)
            assigned_nucleus_indices = clean_indices[valid]
            assigned_spot_original_indices = finite_indices[idxs[valid]]
            nucleus_spot_map = pd.DataFrame(
                {
                    "nucleus_id": assigned_nucleus_indices + 1,  # 1-based nucleus IDs
                    "spot_id": spot_names_arr[assigned_spot_original_indices],
                    "x_pixel": centroids_xy[assigned_nucleus_indices, 0],
                    "y_pixel": centroids_xy[assigned_nucleus_indices, 1],
                }
            )

    # Retrieve QC from the segmenter (stashed by run_nuclei_segmentation)
    qc = getattr(run_nuclei_segmentation, "_last_qc", None)

    return SegmentationResult(
        masks=masks,
        centroids_xy=centroids_xy,
        nuclei_count_raw=raw_counts,
        resolution_mode=resolution_mode,
        image_shape=(int(image.shape[0]), int(image.shape[1])),
        nucleus_spot_map=nucleus_spot_map,
        qc=qc,
    )


def compute_spot_nuclei_counts_patchwise(
    fullres_image: np.ndarray,
    spot_centers_fullres: np.ndarray,
    spot_names: pd.Index,
    pixel_size_um: float,
    *,
    spot_diameter_um: float = 55.0,
    patch_margin_um: float = 25.0,
    modality: str = "he",
    **stardist_kwargs,
) -> SegmentationResult:
    """Spot-patch-based nuclei segmentation for Visium on fullres images.

    Instead of processing the entire fullres image at once (which OOMs),
    this function crops a square patch around each Visium spot at native
    resolution, runs StarDist on each patch, and counts nuclei within the
    spot radius.  This preserves full resolution where it matters and keeps
    memory bounded per-patch.

    Patch size = ``spot_diameter_um + 2 * patch_margin_um`` (default 105 µm
    = ~132 px at 0.79 µm/px).  StarDist is loaded once and reused across
    all patches.

    Args:
        fullres_image: RGB uint8 fullres H&E image.
        spot_centers_fullres: (N, 2) spot centers in fullres pixel coords (x, y).
        spot_names: Index of spot barcodes, length N.
        pixel_size_um: Microns per pixel in the fullres image.
        spot_diameter_um: Physical spot diameter in microns (default 55.0).
        patch_margin_um: Extra margin around spot for context (default 25.0).
        modality: ``"he"`` or ``"dapi"``.
        **stardist_kwargs: Forwarded to ``StarDistSegmenter.segment()``
            (e.g. ``prob_thresh``, ``nms_thresh``).  ``scale`` is auto-set
            from ``pixel_size_um``.

    Returns:
        SegmentationResult with per-spot nuclei counts, centroids in fullres
        pixel coordinates, and a dummy masks array (empty — per-patch masks
        are not stitched).
    """
    segmenter = StarDistSegmenter(modality=modality)

    spot_radius_um = spot_diameter_um / 2.0
    patch_half_um = spot_radius_um + patch_margin_um
    patch_half_px = int(np.ceil(patch_half_um / pixel_size_um))
    spot_radius_px = spot_radius_um / pixel_size_um

    img_h, img_w = fullres_image.shape[:2]
    n_spots = len(spot_names)

    logging.info(
        "Patchwise segmentation: %d spots, patch=%d px (%.0f µm), " "spot_radius=%.1f px, pixel_size=%.4f µm/px",
        n_spots,
        2 * patch_half_px,
        2 * patch_half_um,
        spot_radius_px,
        pixel_size_um,
    )

    counts = np.zeros(n_spots, dtype=np.int64)
    all_centroids = []  # list of (x_fullres, y_fullres, spot_idx)
    total_raw = 0
    total_filtered = 0
    total_small = 0
    total_large = 0
    all_areas_um2 = []

    for i in range(n_spots):
        cx, cy = spot_centers_fullres[i]
        if not (np.isfinite(cx) and np.isfinite(cy)):
            continue

        cx_int, cy_int = int(round(cx)), int(round(cy))

        # Crop patch (clamp to image bounds)
        y0 = max(0, cy_int - patch_half_px)
        y1 = min(img_h, cy_int + patch_half_px)
        x0 = max(0, cx_int - patch_half_px)
        x1 = min(img_w, cx_int + patch_half_px)

        if y1 - y0 < 10 or x1 - x0 < 10:
            continue

        patch = fullres_image[y0:y1, x0:x1]

        # Skip near-blank patches (background/white space).
        # csbdeep.normalize divides by (vmax - vmin); near-constant patches
        # cause numerical instability and can crash TensorFlow.
        if patch.std() < 3.0:
            continue

        # Run StarDist on this patch
        _, centroids_df = segmenter.segment(patch, pixel_size_um=pixel_size_um, **stardist_kwargs)

        # Accumulate QC from segmenter
        patch_qc = getattr(segmenter, "_last_qc", None)
        if patch_qc is not None:
            total_raw += patch_qc.n_raw
            total_filtered += patch_qc.n_after_area_filter
            total_small += patch_qc.n_removed_small
            total_large += patch_qc.n_removed_large

        if len(centroids_df) == 0:
            continue

        # Convert patch-local centroids to fullres coordinates
        cent_x_full = centroids_df["x_pixel"].to_numpy() + x0
        cent_y_full = centroids_df["y_pixel"].to_numpy() + y0

        # Keep only nuclei within the spot radius
        dx = cent_x_full - cx
        dy = cent_y_full - cy
        dist = np.sqrt(dx**2 + dy**2)
        in_spot = dist <= spot_radius_px

        n_in = int(in_spot.sum())
        counts[i] = n_in

        for j in np.where(in_spot)[0]:
            all_centroids.append((cent_x_full[j], cent_y_full[j], i))

        # Collect area stats for QC
        if patch_qc is not None and patch_qc.area_median_um2 > 0:
            all_areas_um2.append(patch_qc.area_median_um2)

    # Build outputs
    counts_series: pd.Series = pd.Series(counts, index=spot_names, name="nuclei_count_raw")

    # Build nucleus-spot mapping
    nucleus_spot_map = None
    if all_centroids:
        cent_arr = np.array(all_centroids)
        spot_names_arr = np.asarray(spot_names)
        nucleus_spot_map = pd.DataFrame(
            {
                "nucleus_id": np.arange(1, len(cent_arr) + 1, dtype=np.int32),
                "spot_id": spot_names_arr[cent_arr[:, 2].astype(int)],
                "x_pixel": cent_arr[:, 0],
                "y_pixel": cent_arr[:, 1],
            }
        )

    # Build aggregate QC
    areas_arr = np.array(all_areas_um2) if all_areas_um2 else np.array([0.0])
    image_area_mm2 = (img_h * pixel_size_um * img_w * pixel_size_um) / 1e6
    total_nuclei = int(counts_series.sum())
    density = total_nuclei / image_area_mm2 if image_area_mm2 > 0 else 0.0
    lo, hi = EXPECTED_DENSITY_RANGE
    if density > hi * DENSITY_HIGH_FACTOR:
        flag = "high"
    elif density < lo * DENSITY_LOW_FACTOR:
        flag = "low"
    else:
        flag = "ok"

    qc = SegmentationQC(
        n_raw=total_raw,
        n_after_area_filter=total_filtered,
        n_removed_small=total_small,
        n_removed_large=total_large,
        pixel_size_um=pixel_size_um,
        min_area_um2=segmenter.DEFAULT_MIN_AREA_UM2,
        max_area_um2=segmenter.DEFAULT_MAX_AREA_UM2,
        area_median_um2=float(np.median(areas_arr)) if len(areas_arr) > 0 else 0.0,
        area_p25_um2=float(np.percentile(areas_arr, 25)) if len(areas_arr) > 0 else 0.0,
        area_p75_um2=float(np.percentile(areas_arr, 75)) if len(areas_arr) > 0 else 0.0,
        density_per_mm2=density,
        density_flag=flag,
    )

    logging.info(
        "Patchwise done: %d nuclei across %d spots (%.1f/spot avg), "
        "%d raw → %d filtered (-%d small, -%d large), density=%.0f/mm² [%s]",
        total_nuclei,
        n_spots,
        total_nuclei / max(1, n_spots),
        total_raw,
        total_filtered,
        total_small,
        total_large,
        density,
        flag,
    )

    # Centroids in fullres xy for downstream compatibility
    if all_centroids:
        cent_all = np.array(all_centroids)[:, :2]
    else:
        cent_all = np.zeros((0, 2), dtype=np.float64)

    return SegmentationResult(
        masks=np.zeros((1, 1), dtype=np.int32),  # dummy — patches not stitched
        centroids_xy=cent_all,
        nuclei_count_raw=counts_series,
        resolution_mode="fullres_patchwise",
        image_shape=(img_h, img_w),
        nucleus_spot_map=nucleus_spot_map,
        qc=qc,
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

    centroids_path = out / f"{sample_name}_centroids.csv"
    centroids_df = pd.DataFrame(
        {
            "nucleus_id": np.arange(1, len(result.centroids_xy) + 1),
            "x_pixel": result.centroids_xy[:, 0],
            "y_pixel": result.centroids_xy[:, 1],
        }
    )
    centroids_df.to_csv(centroids_path, index=False)
    outputs["centroids"] = str(centroids_path)

    if result.nucleus_spot_map is not None:
        map_path = out / f"{sample_name}_nucleus_spot_map.csv"
        result.nucleus_spot_map.to_csv(map_path, index=False)
        outputs["nucleus_spot_map"] = str(map_path)

    if result.qc is not None:
        import json  # pylint: disable=import-outside-toplevel

        qc_path = out / f"{sample_name}_segmentation_qc.json"
        with open(qc_path, "w") as fh:
            json.dump(result.qc.to_dict(), fh, indent=2)
        outputs["segmentation_qc"] = str(qc_path)

    return outputs
