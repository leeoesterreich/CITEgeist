#!/usr/bin/env python
"""
Generate Cellpose-compatible synthetic images from scCube simulation data.

Supports two image modes:
- dapi: Grayscale nuclei on black background (for Cellpose 'nuclei' model)
- h_and_e: H&E-style with purple nuclei and pink cytoplasm (for Cellpose 'cyto2' model)
"""

import argparse
import logging
from pathlib import Path
from typing import Tuple

import numpy as np
import pandas as pd
from PIL import Image

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)

# Default paths
DEFAULT_INPUT_DIR = Path(
    "/ix1/alee/LO_LAB/Personal/Brent_Schlegel/bts76/Projects/CITEgeist/"
    "Wu_Visium/Simulations/scCube_12k/replicates"
)

# Color constants (from actual Visium H&E and Xenium DAPI analysis)
DAPI_BACKGROUND = (0, 0, 0)
DAPI_NUCLEUS_INTENSITY = 180  # Grayscale intensity

HE_BACKGROUND = (250, 250, 250)
HE_NUCLEUS_COLOR = (140, 90, 130)  # Purple/magenta (hematoxylin)
HE_CYTOPLASM_COLOR = (240, 220, 220)  # Pale pink (eosin)


def create_gaussian_kernel(sigma: float, size: int = None) -> np.ndarray:
    """
    Create a 2D Gaussian kernel for nucleus rendering.

    Args:
        sigma: Standard deviation of Gaussian
        size: Kernel size (default: 6*sigma, must be odd)

    Returns:
        2D numpy array with Gaussian values normalized to [0, 1]
    """
    if size is None:
        size = int(np.ceil(6 * sigma))
        if size % 2 == 0:
            size += 1

    center = size // 2
    y, x = np.ogrid[:size, :size]
    kernel = np.exp(-((x - center) ** 2 + (y - center) ** 2) / (2 * sigma ** 2))
    return kernel


def generate_dapi_image(
    cell_coords: np.ndarray,
    image_size: Tuple[int, int],
    sigma: float = 3.5,
    padding: int = 100,
) -> np.ndarray:
    """
    Generate DAPI-style grayscale image with Gaussian nuclei.

    Args:
        cell_coords: Nx2 array of (x, y) coordinates in original units
        image_size: (width, height) of output image
        sigma: Gaussian sigma for nuclei
        padding: Padding in pixels around coordinate extent

    Returns:
        RGB uint8 numpy array (grayscale stored as RGB)
    """
    width, height = image_size

    # Compute coordinate scaling
    x_min, x_max = cell_coords[:, 0].min(), cell_coords[:, 0].max()
    y_min, y_max = cell_coords[:, 1].min(), cell_coords[:, 1].max()

    scale_x = (width - 2 * padding) / (x_max - x_min + 1e-6)
    scale_y = (height - 2 * padding) / (y_max - y_min + 1e-6)

    # Create float accumulator for Gaussian blending
    img_float = np.zeros((height, width), dtype=np.float32)

    # Create Gaussian kernel once
    kernel = create_gaussian_kernel(sigma)
    k_size = kernel.shape[0]
    k_half = k_size // 2

    # Draw each nucleus
    for x_orig, y_orig in cell_coords:
        # Transform to image coordinates
        x_px = int(padding + (x_orig - x_min) * scale_x)
        y_px = int(padding + (y_orig - y_min) * scale_y)

        # Compute kernel bounds with clipping
        x_start = max(0, x_px - k_half)
        x_end = min(width, x_px + k_half + 1)
        y_start = max(0, y_px - k_half)
        y_end = min(height, y_px + k_half + 1)

        # Corresponding kernel region
        kx_start = x_start - (x_px - k_half)
        kx_end = k_size - ((x_px + k_half + 1) - x_end)
        ky_start = y_start - (y_px - k_half)
        ky_end = k_size - ((y_px + k_half + 1) - y_end)

        # Add Gaussian to image (max blending for overlapping nuclei)
        if x_end > x_start and y_end > y_start:
            img_float[y_start:y_end, x_start:x_end] = np.maximum(
                img_float[y_start:y_end, x_start:x_end],
                kernel[ky_start:ky_end, kx_start:kx_end] * DAPI_NUCLEUS_INTENSITY,
            )

    # Convert to uint8 RGB (grayscale replicated to 3 channels)
    img_gray = np.clip(img_float, 0, 255).astype(np.uint8)
    img_rgb = np.stack([img_gray, img_gray, img_gray], axis=-1)

    return img_rgb


def generate_he_image(
    cell_coords: np.ndarray,
    image_size: Tuple[int, int],
    nucleus_sigma: float = 3.5,
    cytoplasm_radius: int = 10,
    padding: int = 100,
) -> np.ndarray:
    """
    Generate H&E-style image with pink cytoplasm and purple nuclei.

    Args:
        cell_coords: Nx2 array of (x, y) coordinates in original units
        image_size: (width, height) of output image
        nucleus_sigma: Gaussian sigma for nuclei
        cytoplasm_radius: Radius of cytoplasm circle in pixels
        padding: Padding in pixels around coordinate extent

    Returns:
        RGB uint8 numpy array
    """
    width, height = image_size

    # Compute coordinate scaling
    x_min, x_max = cell_coords[:, 0].min(), cell_coords[:, 0].max()
    y_min, y_max = cell_coords[:, 1].min(), cell_coords[:, 1].max()

    scale_x = (width - 2 * padding) / (x_max - x_min + 1e-6)
    scale_y = (height - 2 * padding) / (y_max - y_min + 1e-6)

    # Start with white background
    img = np.full((height, width, 3), HE_BACKGROUND, dtype=np.uint8)

    # Transform all coordinates to pixel space
    coords_px = []
    for x_orig, y_orig in cell_coords:
        x_px = int(padding + (x_orig - x_min) * scale_x)
        y_px = int(padding + (y_orig - y_min) * scale_y)
        coords_px.append((x_px, y_px))

    # First pass: draw cytoplasm circles (pink)
    for x_px, y_px in coords_px:
        y_indices, x_indices = np.ogrid[
            max(0, y_px - cytoplasm_radius) : min(height, y_px + cytoplasm_radius + 1),
            max(0, x_px - cytoplasm_radius) : min(width, x_px + cytoplasm_radius + 1),
        ]
        # Create mask for circular cytoplasm
        local_y = y_indices - y_px
        local_x = x_indices - x_px
        mask = (local_x ** 2 + local_y ** 2) <= cytoplasm_radius ** 2

        # Apply cytoplasm color where mask is True
        y_start = max(0, y_px - cytoplasm_radius)
        y_end = min(height, y_px + cytoplasm_radius + 1)
        x_start = max(0, x_px - cytoplasm_radius)
        x_end = min(width, x_px + cytoplasm_radius + 1)

        for c in range(3):
            region = img[y_start:y_end, x_start:x_end, c]
            region[mask] = HE_CYTOPLASM_COLOR[c]

    # Second pass: draw nuclei (purple Gaussian)
    kernel = create_gaussian_kernel(nucleus_sigma)
    k_size = kernel.shape[0]
    k_half = k_size // 2

    for x_px, y_px in coords_px:
        x_start = max(0, x_px - k_half)
        x_end = min(width, x_px + k_half + 1)
        y_start = max(0, y_px - k_half)
        y_end = min(height, y_px + k_half + 1)

        kx_start = x_start - (x_px - k_half)
        kx_end = k_size - ((x_px + k_half + 1) - x_end)
        ky_start = y_start - (y_px - k_half)
        ky_end = k_size - ((y_px + k_half + 1) - y_end)

        if x_end > x_start and y_end > y_start:
            k_region = kernel[ky_start:ky_end, kx_start:kx_end]
            for c in range(3):
                bg = img[y_start:y_end, x_start:x_end, c].astype(np.float32)
                nucleus = np.full_like(bg, HE_NUCLEUS_COLOR[c])
                # Alpha blend: result = nucleus * alpha + background * (1 - alpha)
                blended = nucleus * k_region + bg * (1 - k_region)
                img[y_start:y_end, x_start:x_end, c] = np.clip(blended, 0, 255).astype(
                    np.uint8
                )

    return img