"""
Cellpose-based nuclei segmentation utilities for Visium-style datasets.

This module provides optional image segmentation support for generating
spot-level nuclei counts that can be used as a soft abundance prior during
cell proportion optimization.
"""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

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
            "Using explicit spot diameter: %.1f µm = %.1f pixels (pixel_size=%.4f µm/px)",
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
            "Detected platform '%s': spot diameter = %.1f µm = %.1f pixels",
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
        "  - Visium: 55 µm diameter, 100 µm center-to-center spacing\n"
        "  - Visium HD: 8 µm bins (contiguous)\n"
        "  - Xenium pseudo-Visium: typically 55 µm to match Visium geometry"
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
    """Ensure image is uint8 RGB for robust Cellpose input."""
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


def run_cellpose_nuclei_segmentation(
    image_rgb_uint8: np.ndarray,
    use_gpu: bool = False,
    diameter: Optional[float] = None,
    flow_threshold: float = 0.4,
    cellprob_threshold: float = 0.0,
    model=None,
    model_type: str = "nuclei",
) -> Tuple[np.ndarray, np.ndarray]:
    """Run Cellpose nuclei model and return mask labels + centroid array."""
    if model is None:
        model = _build_cellpose_model(use_gpu=use_gpu, model_type=model_type)

    eval_out = model.eval(
        image_rgb_uint8,
        channels=[0, 0],
        diameter=diameter,
        flow_threshold=flow_threshold,
        cellprob_threshold=cellprob_threshold,
    )
    if isinstance(eval_out, tuple):
        masks = eval_out[0]
    else:
        masks = eval_out

    masks = np.asarray(masks)
    if masks.ndim != 2:
        raise ValueError(f"Cellpose masks should be 2D label image, got shape {masks.shape}.")

    max_label = int(masks.max())
    if max_label <= 0:
        centroids = np.zeros((0, 2), dtype=np.float64)
    else:
        # center_of_mass returns (row, col) -> convert to (x, y).
        labels = np.arange(1, max_label + 1, dtype=np.int32)
        com = center_of_mass(np.ones_like(masks, dtype=np.float32), labels=masks, index=labels)
        centroids = np.array([(float(c[1]), float(c[0])) for c in com], dtype=np.float64)

    return masks, centroids


def _build_cellpose_model(use_gpu: bool = False, model_type: str = "nuclei"):
    try:
        from cellpose import models
    except Exception as exc:
        raise ImportError(
            "Cellpose is required for segmentation. Install with `pip install \"citegeist[imaging]\"` "
            "or `pip install cellpose scikit-image`."
        ) from exc

    # Support both legacy Cellpose API (`Cellpose`) and newer API (`CellposeModel`).
    if hasattr(models, "Cellpose"):
        return models.Cellpose(gpu=use_gpu, model_type=model_type)
    if hasattr(models, "CellposeModel"):
        try:
            return models.CellposeModel(gpu=use_gpu, model_type=model_type)
        except TypeError:
            return models.CellposeModel(gpu=use_gpu, pretrained_model=model_type)
    raise AttributeError("cellpose.models has neither `Cellpose` nor `CellposeModel`; unsupported Cellpose version.")


def _count_nuclei_from_patch_centroids(
    centroids_xy: np.ndarray,
    local_spot_center_xy: Tuple[float, float],
    spot_radius_px: float,
) -> int:
    if centroids_xy.size == 0:
        return 0
    dx = centroids_xy[:, 0] - float(local_spot_center_xy[0])
    dy = centroids_xy[:, 1] - float(local_spot_center_xy[1])
    dist_sq = dx * dx + dy * dy
    return int(np.sum(dist_sq <= (float(spot_radius_px) ** 2)))


def _compute_patch_bounds(cx: float, cy: float, patch_radius_px: float, width: int, height: int) -> Tuple[int, int, int, int]:
    x0 = max(0, int(np.floor(cx - patch_radius_px)))
    x1 = min(width, int(np.ceil(cx + patch_radius_px)))
    y0 = max(0, int(np.floor(cy - patch_radius_px)))
    y1 = min(height, int(np.ceil(cy + patch_radius_px)))
    return x0, x1, y0, y1


def _segment_fullres_by_spot_patches(
    image_rgb_uint8: np.ndarray,
    spot_centers_fullres: np.ndarray,
    spot_radius_px: float,
    spot_names: pd.Index,
    use_gpu: bool = False,
    diameter: Optional[float] = None,
    flow_threshold: float = 0.4,
    cellprob_threshold: float = 0.0,
    patch_radius_multiplier: float = 1.5,
    patch_workers: int = 4,
) -> Tuple[pd.Series, np.ndarray]:
    """Segment fullres image by spot-local patches and count nuclei per spot."""
    if patch_radius_multiplier < 1.0:
        raise ValueError("patch_radius_multiplier must be >= 1.0.")
    if patch_workers < 1:
        raise ValueError("patch_workers must be >= 1.")

    h, w = int(image_rgb_uint8.shape[0]), int(image_rgb_uint8.shape[1])
    patch_radius_px = float(spot_radius_px) * float(patch_radius_multiplier)

    def _process_indices(indices: List[int]) -> Tuple[List[Tuple[int, int]], List[np.ndarray]]:
        model = _build_cellpose_model(use_gpu=use_gpu)
        local_counts: List[Tuple[int, int]] = []
        local_centroids_global: List[np.ndarray] = []
        for idx in indices:
            cx, cy = spot_centers_fullres[idx]
            x0, x1, y0, y1 = _compute_patch_bounds(cx, cy, patch_radius_px, w, h)
            if x1 <= x0 or y1 <= y0:
                local_counts.append((idx, 0))
                continue
            patch = image_rgb_uint8[y0:y1, x0:x1]
            if patch.size == 0:
                local_counts.append((idx, 0))
                continue
            _, centroids_local = run_cellpose_nuclei_segmentation(
                image_rgb_uint8=patch,
                use_gpu=use_gpu,
                diameter=diameter,
                flow_threshold=flow_threshold,
                cellprob_threshold=cellprob_threshold,
                model=model,
            )
            local_center = (float(cx - x0), float(cy - y0))
            count = _count_nuclei_from_patch_centroids(
                centroids_xy=centroids_local,
                local_spot_center_xy=local_center,
                spot_radius_px=spot_radius_px,
            )
            local_counts.append((idx, count))
            if centroids_local.size > 0:
                centroids_global = centroids_local.copy()
                centroids_global[:, 0] += float(x0)
                centroids_global[:, 1] += float(y0)
                local_centroids_global.append(centroids_global)
        return local_counts, local_centroids_global

    n_spots = spot_centers_fullres.shape[0]
    all_idx = list(range(n_spots))
    chunk_size = max(1, int(np.ceil(n_spots / float(patch_workers))))
    chunks = [all_idx[i : i + chunk_size] for i in range(0, n_spots, chunk_size)]

    counts = np.zeros(n_spots, dtype=np.int64)
    centroids_parts: List[np.ndarray] = []

    if patch_workers == 1:
        res_counts, res_centroids = _process_indices(chunks[0])
        for idx, c in res_counts:
            counts[idx] = c
        centroids_parts.extend(res_centroids)
    else:
        with ThreadPoolExecutor(max_workers=patch_workers) as ex:
            futures = [ex.submit(_process_indices, chunk) for chunk in chunks]
            for fut in futures:
                res_counts, res_centroids = fut.result()
                for idx, c in res_counts:
                    counts[idx] = c
                centroids_parts.extend(res_centroids)

    if len(centroids_parts) > 0:
        centroids_all = np.vstack(centroids_parts)
    else:
        centroids_all = np.zeros((0, 2), dtype=np.float64)
    return pd.Series(counts, index=spot_names, name="nuclei_count_raw"), centroids_all


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

    tree = cKDTree(spot_centers_xy)
    dists, idxs = tree.query(centroids_xy, k=1, workers=-1)
    valid = dists <= float(spot_radius_px)
    if np.any(valid):
        np.add.at(counts, idxs[valid], 1)
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
    use_gpu: bool = False,
    diameter: Optional[float] = None,
    flow_threshold: float = 0.4,
    cellprob_threshold: float = 0.0,
    max_fullres_side: int = 9000,
    fullres_patch_mode: bool = False,
    fullres_patch_radius_multiplier: float = 1.5,
    fullres_patch_workers: int = 4,
    spot_diameter_um: Optional[float] = None,
    pixel_size_um: Optional[float] = None,
) -> SegmentationResult:
    """
    End-to-end segmentation + centroid assignment for an AnnData Visium object.

    Spot diameter is determined by (in order of priority):
    1. spot_diameter_um parameter (explicit, requires pixel_size_um)
    2. scalefactors['spot_diameter_fullres'] from adata.uns['spatial']
    3. Platform-based detection from adata.uns['platform'] (requires pixel_size_um)

    Args:
        adata: AnnData with spatial coordinates and histology images
        resolution_mode: Image resolution to use ('lowres', 'hires', 'fullres')
        use_gpu: Whether to use GPU for Cellpose
        diameter: Cellpose diameter parameter (nucleus size, not spot size)
        flow_threshold: Cellpose flow threshold
        cellprob_threshold: Cellpose cell probability threshold
        max_fullres_side: Max image dimension for fullres fallback
        fullres_patch_mode: Use per-spot patch segmentation for fullres
        fullres_patch_radius_multiplier: Patch size as multiple of spot radius
        fullres_patch_workers: Number of parallel workers for patch mode
        spot_diameter_um: Explicit spot diameter in microns (e.g., 55.0 for Visium)
        pixel_size_um: Microns per pixel in coordinate frame (required for um-based params)

    Returns:
        SegmentationResult with nuclei counts per spot
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

    if resolution_mode == "fullres" and fullres_patch_mode and np.isclose(fullres_to_image_scale, 1.0):
        spot_radius_px = float(spot_diam_fullres) / 2.0
        raw_counts, centroids_xy = _segment_fullres_by_spot_patches(
            image_rgb_uint8=image,
            spot_centers_fullres=spot_centers_fullres,
            spot_radius_px=spot_radius_px,
            spot_names=adata.obs_names,
            use_gpu=use_gpu,
            diameter=diameter,
            flow_threshold=flow_threshold,
            cellprob_threshold=cellprob_threshold,
            patch_radius_multiplier=fullres_patch_radius_multiplier,
            patch_workers=fullres_patch_workers,
        )
        masks = np.zeros((1, 1), dtype=np.int32)
    else:
        masks, centroids_xy = run_cellpose_nuclei_segmentation(
            image_rgb_uint8=image,
            use_gpu=use_gpu,
            diameter=diameter,
            flow_threshold=flow_threshold,
            cellprob_threshold=cellprob_threshold,
        )
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
        masks_path = out / f"{sample_name}_cellpose_masks_{result.resolution_mode}.npy"
        np.save(masks_path, result.masks)
        outputs["cellpose_masks_npy"] = str(masks_path)
    return outputs
