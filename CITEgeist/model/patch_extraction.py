"""Nucleus patch extraction for VAE training."""
import numpy as np
import pandas as pd
from typing import Tuple
import cv2


def extract_patch(
    image: np.ndarray,
    bbox: Tuple[int, int, int, int],
    expansion: float = 0.75,
    output_size: int = 96,
) -> np.ndarray:
    """
    Extract expanded patch around a nucleus.

    Args:
        image: (C, H, W) multi-channel image
        bbox: (x_min, y_min, x_max, y_max) nucleus bounding box
        expansion: Fraction to expand bbox in each direction (0.75 = 75%)
        output_size: Final patch size after resize

    Returns:
        patch: (C, output_size, output_size) normalized patch
    """
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

    # Normalize per channel (z-score)
    for c in range(C):
        mean = resized[c].mean()
        std = resized[c].std() + 1e-8
        resized[c] = (resized[c] - mean) / std

    return resized


def extract_patches_for_spot(
    image: np.ndarray,
    nuclei_df: pd.DataFrame,
    expansion: float = 0.75,
    output_size: int = 96,
) -> np.ndarray:
    """
    Extract patches for all nuclei in a spot.

    Args:
        image: (C, H, W) multi-channel image
        nuclei_df: DataFrame with columns ['nucleus_id', 'bbox_x_min',
                   'bbox_y_min', 'bbox_x_max', 'bbox_y_max']
        expansion: Bbox expansion fraction
        output_size: Output patch size

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
        patch = extract_patch(image, bbox, expansion, output_size)
        patches.append(patch)

    return np.stack(patches, axis=0)
