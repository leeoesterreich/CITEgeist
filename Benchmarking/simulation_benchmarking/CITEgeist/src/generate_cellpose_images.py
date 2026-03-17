#!/usr/bin/env python
"""
Generate synthetic nuclei images from scCube simulation data.

Supports two image modes:
- dapi: Grayscale nuclei on black background (for nuclei segmentation models)
- h_and_e: H&E-style with purple nuclei and pink cytoplasm (for cytoplasm segmentation models)
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
# Increased intensity/contrast for better Cellpose detection
DAPI_BACKGROUND = (0, 0, 0)
DAPI_NUCLEUS_INTENSITY = 255  # Maximum intensity for best Cellpose detection

HE_BACKGROUND = (250, 250, 250)
HE_NUCLEUS_COLOR = (80, 40, 90)  # Darker purple for better contrast
HE_CYTOPLASM_COLOR = (240, 200, 210)  # Slightly more saturated pink


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


def create_dapi_kernel(radius: float, edge_sigma: float = 2.0, size: int = None) -> np.ndarray:
    """
    Create a filled disk with soft edge for realistic DAPI rendering.

    Real DAPI nuclei appear as bright filled disks with distinct edges,
    not smooth Gaussian falloffs. This kernel creates a flat-top disk
    with a soft Gaussian edge for more realistic appearance.

    Args:
        radius: Radius of the solid disk portion (pixels)
        edge_sigma: Gaussian sigma for soft edge falloff (pixels)
        size: Kernel size (default: 2*(radius + 3*edge_sigma) + 1)

    Returns:
        2D numpy array with filled disk and soft boundary, normalized to [0, 1]
    """
    if size is None:
        size = int(np.ceil(2 * (radius + 3 * edge_sigma)))
        if size % 2 == 0:
            size += 1

    center = size // 2
    y, x = np.ogrid[:size, :size]
    dist = np.sqrt((x - center) ** 2 + (y - center) ** 2)

    # Solid disk up to radius, then Gaussian falloff
    kernel = np.zeros((size, size), dtype=np.float32)
    inner = dist <= radius
    outer = dist > radius

    kernel[inner] = 1.0
    if edge_sigma > 0:
        kernel[outer] = np.exp(-((dist[outer] - radius) ** 2) / (2 * edge_sigma ** 2))

    return kernel


def generate_dapi_image(
    cell_coords: np.ndarray,
    image_size: Tuple[int, int],
    nucleus_radius: float = 15.0,
    edge_sigma: float = 3.0,
    padding: int = 100,
) -> np.ndarray:
    """
    Generate DAPI-style grayscale image with realistic filled-disk nuclei.

    Uses a flat-top disk kernel with soft edge to mimic real DAPI appearance,
    where nuclei appear as bright filled regions with distinct boundaries.

    Note: Nucleus size (radius=15, edge_sigma=3 -> ~36px total diameter) is
    calibrated for nuclei segmentation models (~30+ pixels).
    Smaller nuclei (e.g., radius=8) result in poor detection accuracy.

    Args:
        cell_coords: Nx2 array of (x, y) coordinates in original units
        image_size: (width, height) of output image
        nucleus_radius: Radius of the solid nucleus disk (pixels)
        edge_sigma: Gaussian sigma for soft edge falloff (pixels)
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

    # Create float accumulator for blending
    img_float = np.zeros((height, width), dtype=np.float32)

    # Create filled-disk kernel with soft edge (realistic DAPI appearance)
    kernel = create_dapi_kernel(nucleus_radius, edge_sigma)
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

        # Add kernel to image (max blending for overlapping nuclei)
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
    nucleus_sigma: float = 10.0,
    cytoplasm_radius: int = 25,
    padding: int = 100,
) -> np.ndarray:
    """
    Generate H&E-style image with pink cytoplasm and purple nuclei.

    Note: Nucleus sigma=10 gives FWHM ~24px and kernel extent ~60px,
    calibrated for cytoplasm segmentation models. Cytoplasm radius=25
    maintains realistic nucleus-to-cytoplasm ratio.

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
    nucleus_radius: float,
    edge_sigma: float,
    he_nucleus_sigma: float,
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
            # Use filled-disk kernel for realistic DAPI
            img_array = generate_dapi_image(
                cell_coords,
                image_size=(image_size, image_size),
                nucleus_radius=nucleus_radius,
                edge_sigma=edge_sigma,
            )
        elif mode == "h_and_e":
            img_array = generate_he_image(
                cell_coords,
                image_size=(image_size, image_size),
                nucleus_sigma=he_nucleus_sigma,
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
        description="Generate synthetic nuclei images from scCube simulation data"
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
        "--nucleus-radius",
        type=float,
        default=15.0,
        help="DAPI nucleus radius in pixels (default: 15.0, ~30px diameter solid disk)",
    )
    parser.add_argument(
        "--edge-sigma",
        type=float,
        default=3.0,
        help="DAPI nucleus edge softness sigma (default: 3.0)",
    )
    parser.add_argument(
        "--he-nucleus-sigma",
        type=float,
        default=10.0,
        help="H&E nucleus Gaussian sigma (default: 10.0, ~24px FWHM)",
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
        nucleus_radius=args.nucleus_radius,
        edge_sigma=args.edge_sigma,
        he_nucleus_sigma=args.he_nucleus_sigma,
    )

    logger.info("Done!")


if __name__ == "__main__":
    main()