"""Nucleus patch extraction for VAE training."""
import numpy as np
import pandas as pd
from typing import Dict, Any, Tuple
import cv2


def compute_global_stats(
    image: np.ndarray,
    norm_method: str = "percentile",
) -> Dict[str, Any]:
    """Compute normalization stats per channel for entire image.

    Args:
        image: (C, H, W) multi-channel image
        norm_method: "percentile" (1st/99th) or "minmax" (full range)

    Returns:
        Dict with normalization parameters
    """
    if norm_method == "percentile":
        p1 = np.percentile(image, 1, axis=(1, 2))
        p99 = np.percentile(image, 99, axis=(1, 2))
        return {"method": "percentile", "p1": p1, "p99": p99}
    elif norm_method == "minmax":
        cmin = image.min(axis=(1, 2))
        cmax = image.max(axis=(1, 2))
        return {"method": "minmax", "min": cmin, "max": cmax}
    else:
        raise ValueError(f"Unknown norm_method: {norm_method}")


def extract_patch(
    image: np.ndarray,
    bbox: Tuple[int, int, int, int],
    expansion: float = 0.75,
    output_size: int = 96,
    global_stats: Dict[str, Any] = None,
) -> np.ndarray:
    """Extract expanded patch around a nucleus with global normalization.

    Args:
        image: (C, H, W) multi-channel image
        bbox: (x_min, y_min, x_max, y_max) nucleus bounding box
        expansion: Fraction to expand bbox in each direction (0.75 = 75%)
        output_size: Final patch size after resize
        global_stats: REQUIRED - normalization stats from compute_global_stats()

    Returns:
        patch: (C, output_size, output_size) normalized patch in [0, 1] range

    Raises:
        ValueError: If global_stats is None
    """
    if global_stats is None:
        raise ValueError(
            "global_stats is required. Per-patch normalization is deprecated. "
            "Pass global_stats from compute_global_stats()."
        )

    x_min, y_min, x_max, y_max = bbox
    C, H, W = image.shape

    # Compute expansion
    w = x_max - x_min
    h = y_max - y_min

    exp_w = int(w * expansion)
    exp_h = int(h * expansion)

    # Expand bbox
    x_min_exp = max(0, x_min - exp_w)
    x_max_exp = min(W, x_max + exp_w)
    y_min_exp = max(0, y_min - exp_h)
    y_max_exp = min(H, y_max + exp_h)

    # Crop
    patch = image[:, y_min_exp:y_max_exp, x_min_exp:x_max_exp]

    # Resize each channel
    resized = np.zeros((C, output_size, output_size), dtype=np.float32)
    for c in range(C):
        resized[c] = cv2.resize(
            patch[c],
            (output_size, output_size),
            interpolation=cv2.INTER_LINEAR
        )

    # Apply global normalization
    method = global_stats["method"]
    for c in range(C):
        if method == "percentile":
            p1, p99 = global_stats["p1"][c], global_stats["p99"][c]
            clipped = np.clip(resized[c], p1, p99)
            resized[c] = (clipped - p1) / (p99 - p1 + 1e-8)
        elif method == "minmax":
            cmin, cmax = global_stats["min"][c], global_stats["max"][c]
            resized[c] = (resized[c] - cmin) / (cmax - cmin + 1e-8)

    return resized


def extract_patches_for_spot(
    image: np.ndarray,
    nuclei_df: pd.DataFrame,
    expansion: float = 0.75,
    output_size: int = 96,
    global_stats: Dict[str, Any] = None,
) -> np.ndarray:
    """
    Extract patches for all nuclei in a spot.

    Args:
        image: (C, H, W) multi-channel image
        nuclei_df: DataFrame with columns ['nucleus_id', 'bbox_x_min',
                   'bbox_y_min', 'bbox_x_max', 'bbox_y_max']
        expansion: Bbox expansion fraction
        output_size: Output patch size
        global_stats: REQUIRED - normalization stats from compute_global_stats()

    Returns:
        patches: (N, C, output_size, output_size) array of patches
    """
    patches = []

    for _, row in nuclei_df.iterrows():
        bbox = (
            int(row['bbox_x_min']),
            int(row['bbox_y_min']),
            int(row['bbox_x_max']),
            int(row['bbox_y_max']),
        )
        patch = extract_patch(image, bbox, expansion, output_size, global_stats)
        patches.append(patch)

    return np.stack(patches, axis=0)
