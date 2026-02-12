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


def load_cell_data(replicate_id: int, condition: str, input_dir: Path) -> pd.DataFrame:
    """
    Load cell-to-spot index data from scCube simulation.

    Args:
        replicate_id: Replicate number (0-4)
        condition: 'high_seg' or 'mixed'
        input_dir: Base path to scCube replicates directory

    Returns:
        DataFrame with Cell, Cell_type, point_x, point_y, spot columns
    """
    index_file = input_dir / condition / "ST_sim" / f"Wu_ST_{replicate_id}_index.csv"

    if not index_file.exists():
        raise FileNotFoundError(f"Index file not found: {index_file}")

    df = pd.read_csv(index_file, index_col=0)
    logger.info("Loaded %d cells from %s", len(df), index_file)

    # Log cell type distribution
    for ct, count in df["Cell_type"].value_counts().items():
        logger.info("  %s: %d", ct, count)

    return df


def compute_nuclei_counts(cell_df: pd.DataFrame) -> pd.Series:
    """
    Compute ground truth nuclei counts per spot.

    Args:
        cell_df: DataFrame with 'spot' column

    Returns:
        Series with spot names as index, cell counts as values
    """
    counts = cell_df.groupby("spot").size()
    counts.name = "nuclei_count"
    logger.info(
        "Nuclei counts: min=%d, max=%d, mean=%.1f, total=%d",
        counts.min(),
        counts.max(),
        counts.mean(),
        counts.sum(),
    )
    return counts


def generate_images_for_replicate(
    replicate_id: int,
    condition: str,
    input_dir: Path,
    output_dir: Path,
    modes: list,
    image_size: int,
    nucleus_sigma: float,
) -> None:
    """Generate images and nuclei counts for a single replicate."""
    logger.info("=== Replicate %d, Condition: %s ===", replicate_id, condition)

    # Load cell data
    cell_df = load_cell_data(replicate_id, condition, input_dir)
    cell_coords = cell_df[["point_x", "point_y"]].values

    # Compute and save nuclei counts
    nuclei_counts = compute_nuclei_counts(cell_df)
    counts_dir = output_dir / condition / "nuclei_counts"
    counts_dir.mkdir(parents=True, exist_ok=True)
    counts_path = counts_dir / f"Wu_rep_{replicate_id}_nuclei_counts.csv"
    nuclei_counts.to_csv(counts_path, header=True)
    logger.info("Saved nuclei counts: %s", counts_path)

    # Generate images for each mode
    for mode in modes:
        logger.info("Generating %s image...", mode)

        if mode == "dapi":
            img_array = generate_dapi_image(
                cell_coords,
                image_size=(image_size, image_size),
                sigma=nucleus_sigma,
            )
        elif mode == "h_and_e":
            img_array = generate_he_image(
                cell_coords,
                image_size=(image_size, image_size),
                nucleus_sigma=nucleus_sigma,
            )
        else:
            raise ValueError(f"Unknown mode: {mode}")

        # Save image
        img_dir = output_dir / condition / "images" / mode
        img_dir.mkdir(parents=True, exist_ok=True)
        img_path = img_dir / f"Wu_rep_{replicate_id}_cellpose.png"

        Image.fromarray(img_array).save(img_path)
        logger.info("Saved: %s (%dx%d)", img_path, img_array.shape[1], img_array.shape[0])


def main():
    parser = argparse.ArgumentParser(
        description="Generate Cellpose-compatible images from scCube simulation data"
    )
    parser.add_argument(
        "--replicate-id",
        type=int,
        required=True,
        help="Replicate ID (0-4)",
    )
    parser.add_argument(
        "--condition",
        choices=["high_seg", "mixed"],
        required=True,
        help="Simulation condition",
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=DEFAULT_INPUT_DIR,
        help="Path to scCube replicates directory",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Output directory for images and nuclei counts",
    )
    parser.add_argument(
        "--mode",
        choices=["dapi", "h_and_e", "both"],
        default="both",
        help="Image mode (default: both)",
    )
    parser.add_argument(
        "--image-size",
        type=int,
        default=8000,
        help="Image size in pixels (default: 8000)",
    )
    parser.add_argument(
        "--nucleus-sigma",
        type=float,
        default=3.5,
        help="Gaussian sigma for nuclei (default: 3.5)",
    )
    args = parser.parse_args()

    modes = ["dapi", "h_and_e"] if args.mode == "both" else [args.mode]

    generate_images_for_replicate(
        replicate_id=args.replicate_id,
        condition=args.condition,
        input_dir=args.input_dir,
        output_dir=args.output_dir,
        modes=modes,
        image_size=args.image_size,
        nucleus_sigma=args.nucleus_sigma,
    )

    logger.info("Done!")


if __name__ == "__main__":
    main()