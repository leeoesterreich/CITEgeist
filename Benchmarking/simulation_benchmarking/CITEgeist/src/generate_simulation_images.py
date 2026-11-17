#!/usr/bin/env python
"""
Generate StarDist-tuned synthetic nuclei images from scCube simulation data.

Supports two image modes:
- dapi: Grayscale nuclei on black background with per-nucleus size variation
  and Poisson background noise (tuned for StarDist segmentation)
- h_and_e: H&E-style with elliptical nuclei, Perlin-like texture background,
  and per-nucleus color variation

Key differences from generate_cellpose_images.py:
- DAPI: per-nucleus radius variation (truncated normal ±20%), Poisson background noise
- H&E: Perlin-like fractal noise background, elliptical nuclei with random orientation,
  per-nucleus hematoxylin (±15%) and eosin/cytoplasm (±10%) color variation
- Output filenames use *_stardist.png suffix
- Reproducible via --seed argument
"""

import sys as _sys
from pathlib import Path as _Path

_sys.path.insert(0, str(_Path(__file__).resolve().parents[3]))
import _paths as _cg  # noqa: E402

import argparse
import logging
from pathlib import Path
from typing import Tuple

import numpy as np
import pandas as pd
from PIL import Image
from scipy.ndimage import zoom
from scipy.stats import truncnorm

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)

# Default paths
DEFAULT_INPUT_DIR = Path(
    str(_cg.root("DATA") / "Wu_Visium" / "Simulations" / "scCube_12k" / "replicates")
)

# Color constants (from actual Visium H&E and Xenium DAPI analysis)
DAPI_BACKGROUND = (0, 0, 0)
DAPI_NUCLEUS_INTENSITY = 255  # Maximum intensity for best StarDist detection

HE_BACKGROUND = (250, 250, 250)
HE_NUCLEUS_COLOR = (80, 40, 90)  # Darker purple for better contrast
HE_CYTOPLASM_COLOR = (240, 200, 210)  # Slightly more saturated pink


# ---------------------------------------------------------------------------
# Kernel utilities
# ---------------------------------------------------------------------------



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
        kernel[outer] = np.exp(-((dist[outer] - radius) ** 2) / (2 * edge_sigma**2))

    return kernel


def create_elliptical_dapi_kernel(
    radius: float,
    axis_ratio: float,
    angle_deg: float,
    edge_sigma: float = 2.0,
) -> np.ndarray:
    """
    Create a filled elliptical disk with soft edge for DAPI rendering.

    Args:
        radius: Semi-major axis radius in pixels
        axis_ratio: Ratio of minor to major axis (0.7-1.0 gives mild ellipses)
        angle_deg: Rotation angle in degrees
        edge_sigma: Gaussian sigma for soft edge falloff (pixels)

    Returns:
        2D numpy array with filled ellipse and soft boundary, normalized to [0, 1]
    """
    size = int(np.ceil(2 * (radius + 3 * edge_sigma)))
    if size % 2 == 0:
        size += 1

    center = size // 2
    y, x = np.ogrid[:size, :size]
    # Shift to center
    yc = y - center
    xc = x - center

    # Rotate coordinates
    angle_rad = np.deg2rad(angle_deg)
    cos_a, sin_a = np.cos(angle_rad), np.sin(angle_rad)
    xr = cos_a * xc + sin_a * yc
    yr = -sin_a * xc + cos_a * yc

    # Ellipse distance: semi-major=radius, semi-minor=radius*axis_ratio
    semi_minor = radius * axis_ratio
    dist_ellipse = np.sqrt((xr / radius) ** 2 + (yr / semi_minor) ** 2) * radius

    kernel = np.zeros((size, size), dtype=np.float32)
    inner = dist_ellipse <= radius
    outer = dist_ellipse > radius

    kernel[inner] = 1.0
    if edge_sigma > 0:
        kernel[outer] = np.exp(-((dist_ellipse[outer] - radius) ** 2) / (2 * edge_sigma**2))

    return kernel


# ---------------------------------------------------------------------------
# Perlin-like background
# ---------------------------------------------------------------------------


def generate_perlin_background(
    width: int,
    height: int,
    scale: int = 50,
    octaves: int = 4,
    seed: int = 42,
) -> np.ndarray:
    """
    Generate Perlin-like (fractal) noise using multi-octave scipy zoom upscaling.

    Args:
        width: Image width in pixels
        height: Image height in pixels
        scale: Base spatial scale for noise (higher = coarser features)
        octaves: Number of frequency octaves to sum
        seed: Random seed for reproducibility

    Returns:
        2D float32 array with values in [0, 1]
    """
    rng = np.random.RandomState(seed)
    noise = np.zeros((height, width), dtype=np.float32)
    for octave in range(octaves):
        freq = 2**octave
        amp = 0.5**octave
        small_h = max(2, height // max(2, scale // freq))
        small_w = max(2, width // max(2, scale // freq))
        small = rng.randn(small_h, small_w).astype(np.float32)
        big = zoom(small, (height / small_h, width / small_w), order=1)
        noise += amp * big[:height, :width]
    noise = (noise - noise.min()) / (noise.max() - noise.min() + 1e-8)
    return noise


# ---------------------------------------------------------------------------
# Image generation
# ---------------------------------------------------------------------------


def generate_dapi_image(
    cell_coords: np.ndarray,
    image_size: Tuple[int, int],
    nucleus_radius: float = 15.0,
    edge_sigma: float = 3.0,
    padding: int = 100,
    radius_std: float = 0.2,
    noise_lambda: float = 3.0,
    rng_seed: int = 42,
) -> np.ndarray:
    """
    Generate DAPI-style grayscale image with per-nucleus size variation and
    Poisson background noise (StarDist-tuned).

    Each nucleus receives an independently varied radius drawn from a truncated
    normal distribution (±2 sigma clipped, scale=radius_std). After all nuclei
    are placed, Poisson background noise is blended in using np.maximum so that
    nuclei always dominate.

    Args:
        cell_coords: Nx2 array of (x, y) coordinates in original units
        image_size: (width, height) of output image
        nucleus_radius: Mean radius of the solid nucleus disk (pixels)
        edge_sigma: Gaussian sigma for soft edge falloff (pixels)
        padding: Padding in pixels around coordinate extent
        radius_std: Fractional standard deviation for per-nucleus radius variation
            (0.2 = ±20%, truncated at ±2 sigma)
        noise_lambda: Lambda parameter for Poisson background noise; values are
            clipped to [0, 15] before blending
        rng_seed: Random seed for reproducibility

    Returns:
        RGB uint8 numpy array (grayscale stored as RGB)
    """
    rng = np.random.default_rng(rng_seed)
    width, height = image_size

    # Compute coordinate scaling
    x_min, x_max = cell_coords[:, 0].min(), cell_coords[:, 0].max()
    y_min, y_max = cell_coords[:, 1].min(), cell_coords[:, 1].max()

    scale_x = (width - 2 * padding) / (x_max - x_min + 1e-6)
    scale_y = (height - 2 * padding) / (y_max - y_min + 1e-6)

    # Create float accumulator for blending
    img_float = np.zeros((height, width), dtype=np.float32)

    # Draw each nucleus with independently varied radius
    for x_orig, y_orig in cell_coords:
        # Sample per-nucleus radius variation from truncated normal
        variation = truncnorm.rvs(a=-2, b=2, scale=radius_std, random_state=rng.integers(0, 2**31))
        r = nucleus_radius * (1.0 + variation)
        r = max(r, 1.0)  # Guard against degenerate radius

        # Build per-nucleus kernel
        kernel = create_dapi_kernel(r, edge_sigma)
        k_size = kernel.shape[0]
        k_half = k_size // 2

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

        # Max blending so overlapping nuclei and brighter region wins
        if x_end > x_start and y_end > y_start:
            img_float[y_start:y_end, x_start:x_end] = np.maximum(
                img_float[y_start:y_end, x_start:x_end],
                kernel[ky_start:ky_end, kx_start:kx_end] * DAPI_NUCLEUS_INTENSITY,
            )

    # Add Poisson background noise; nuclei still dominate via np.maximum
    bg_noise = rng.poisson(lam=noise_lambda, size=(height, width)).astype(np.float32)
    bg_noise = bg_noise.clip(0, 15)
    img_float = np.maximum(img_float, bg_noise)

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
    axis_ratio_min: float = 0.7,
    intensity_std: float = 0.15,
    rng_seed: int = 42,
) -> np.ndarray:
    """
    Generate H&E-style image with Perlin-like background, elliptical nuclei,
    and per-nucleus color variation (StarDist-tuned).

    Enhancements over the Cellpose-tuned version:
    - Perlin-like fractal noise background instead of flat white
    - Elliptical nuclei with random orientation (uniform in [0, 180) degrees)
      and random axis ratio in [axis_ratio_min, 1.0]
    - Per-nucleus hematoxylin intensity variation: ±intensity_std (truncated normal)
    - Per-nucleus eosin/cytoplasm color variation: ±10% (truncated normal)

    Args:
        cell_coords: Nx2 array of (x, y) coordinates in original units
        image_size: (width, height) of output image
        nucleus_sigma: Gaussian sigma controlling nucleus kernel smoothness;
            also used as the mean semi-major radius for the H&E ellipse
        cytoplasm_radius: Radius of cytoplasm circle in pixels
        padding: Padding in pixels around coordinate extent
        axis_ratio_min: Minimum minor/major axis ratio (0.7 gives mild ellipses)
        intensity_std: Fractional std for per-nucleus hematoxylin intensity variation
        rng_seed: Random seed for reproducibility

    Returns:
        RGB uint8 numpy array
    """
    rng = np.random.default_rng(rng_seed)
    width, height = image_size

    # Compute coordinate scaling
    x_min, x_max = cell_coords[:, 0].min(), cell_coords[:, 0].max()
    y_min, y_max = cell_coords[:, 1].min(), cell_coords[:, 1].max()

    scale_x = (width - 2 * padding) / (x_max - x_min + 1e-6)
    scale_y = (height - 2 * padding) / (y_max - y_min + 1e-6)

    # --- Perlin-like background ---
    perlin = generate_perlin_background(width, height, scale=50, octaves=4, seed=rng_seed)
    # Blend into near-white background: small amplitude texture (±8 grey levels)
    bg_texture = (perlin * 16 - 8).astype(np.float32)  # [-8, +8] grey variation
    img = np.zeros((height, width, 3), dtype=np.float32)
    for c in range(3):
        img[:, :, c] = np.clip(HE_BACKGROUND[c] + bg_texture, 0, 255)

    # Transform all coordinates to pixel space
    coords_px = []
    for x_orig, y_orig in cell_coords:
        x_px = int(padding + (x_orig - x_min) * scale_x)
        y_px = int(padding + (y_orig - y_min) * scale_y)
        coords_px.append((x_px, y_px))

    n_cells = len(coords_px)

    # Sample per-nucleus parameters
    # Axis ratios: uniform in [axis_ratio_min, 1.0]
    axis_ratios = rng.uniform(axis_ratio_min, 1.0, size=n_cells)
    # Angles: uniform in [0, 180)
    angles = rng.uniform(0, 180, size=n_cells)
    # Hematoxylin intensity multipliers: truncated normal ±intensity_std
    hem_mults = np.array(
        [
            1.0 + truncnorm.rvs(a=-2, b=2, scale=intensity_std, random_state=int(rng.integers(0, 2**31)))
            for _ in range(n_cells)
        ]
    )
    # Cytoplasm color multipliers: ±10%
    cyto_mults = np.array(
        [1.0 + truncnorm.rvs(a=-2, b=2, scale=0.10, random_state=int(rng.integers(0, 2**31))) for _ in range(n_cells)]
    )

    # --- First pass: draw cytoplasm ellipses (eosin/pink) ---
    for i, (x_px, y_px) in enumerate(coords_px):
        cyto_color = tuple(int(np.clip(HE_CYTOPLASM_COLOR[c] * cyto_mults[i], 0, 255)) for c in range(3))

        y_indices, x_indices = np.ogrid[
            max(0, y_px - cytoplasm_radius) : min(height, y_px + cytoplasm_radius + 1),
            max(0, x_px - cytoplasm_radius) : min(width, x_px + cytoplasm_radius + 1),
        ]
        local_y = y_indices - y_px
        local_x = x_indices - x_px
        mask = (local_x**2 + local_y**2) <= cytoplasm_radius**2

        y_start = max(0, y_px - cytoplasm_radius)
        y_end = min(height, y_px + cytoplasm_radius + 1)
        x_start = max(0, x_px - cytoplasm_radius)
        x_end = min(width, x_px + cytoplasm_radius + 1)

        for c in range(3):
            region = img[y_start:y_end, x_start:x_end, c]
            region[mask] = cyto_color[c]

    # --- Second pass: draw elliptical nuclei (hematoxylin/purple) ---
    for i, (x_px, y_px) in enumerate(coords_px):
        hem_color = tuple(int(np.clip(HE_NUCLEUS_COLOR[c] * hem_mults[i], 0, 255)) for c in range(3))

        # Build per-nucleus elliptical Gaussian kernel
        kernel = create_elliptical_dapi_kernel(
            radius=nucleus_sigma,
            axis_ratio=axis_ratios[i],
            angle_deg=angles[i],
            edge_sigma=2.0,
        )
        k_size = kernel.shape[0]
        k_half = k_size // 2

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
                bg = img[y_start:y_end, x_start:x_end, c]
                nucleus_val = float(hem_color[c])
                blended = nucleus_val * k_region + bg * (1.0 - k_region)
                img[y_start:y_end, x_start:x_end, c] = blended

    return np.clip(img, 0, 255).astype(np.uint8)


# ---------------------------------------------------------------------------
# Data loading / replicate orchestration (kept identical to cellpose version)
# ---------------------------------------------------------------------------


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
    radius_std: float,
    noise_lambda: float,
    axis_ratio_min: float,
    intensity_std: float,
    seed: int,
) -> None:
    """Generate StarDist-tuned images and nuclei counts for a single replicate."""
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

    # Per-replicate seed offset for reproducible but distinct replicates
    rep_seed = seed + replicate_id

    # Generate images for each mode
    for mode in modes:
        logger.info("Generating %s image (StarDist-tuned, seed=%d)...", mode, rep_seed)

        if mode == "dapi":
            img_array = generate_dapi_image(
                cell_coords,
                image_size=(image_size, image_size),
                nucleus_radius=nucleus_radius,
                edge_sigma=edge_sigma,
                radius_std=radius_std,
                noise_lambda=noise_lambda,
                rng_seed=rep_seed,
            )
        elif mode == "h_and_e":
            img_array = generate_he_image(
                cell_coords,
                image_size=(image_size, image_size),
                nucleus_sigma=he_nucleus_sigma,
                axis_ratio_min=axis_ratio_min,
                intensity_std=intensity_std,
                rng_seed=rep_seed,
            )
        else:
            raise ValueError(f"Unknown mode: {mode}")

        # Save image with stardist suffix
        img_dir = output_dir / condition / "images" / mode
        img_dir.mkdir(parents=True, exist_ok=True)
        img_path = img_dir / f"Wu_rep_{replicate_id}_stardist.png"

        Image.fromarray(img_array).save(img_path)
        logger.info("Saved: %s (%dx%d)", img_path, img_array.shape[1], img_array.shape[0])


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser(
        description="Generate StarDist-tuned synthetic nuclei images from scCube simulation data"
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
        help="DAPI nucleus mean radius in pixels (default: 15.0, ~30px diameter solid disk)",
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
        help="H&E nucleus semi-major radius (default: 10.0)",
    )
    parser.add_argument(
        "--radius-std",
        type=float,
        default=0.2,
        help="Fractional std for DAPI per-nucleus radius variation (default: 0.2 = ±20%%)",
    )
    parser.add_argument(
        "--noise-lambda",
        type=float,
        default=3.0,
        help="Poisson background noise lambda for DAPI (default: 3.0)",
    )
    parser.add_argument(
        "--axis-ratio-min",
        type=float,
        default=0.7,
        help="Minimum minor/major axis ratio for H&E elliptical nuclei (default: 0.7)",
    )
    parser.add_argument(
        "--intensity-std",
        type=float,
        default=0.15,
        help="Fractional std for H&E hematoxylin intensity variation (default: 0.15 = ±15%%)",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Base random seed (per-replicate seed = seed + replicate_id, default: 42)",
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
        radius_std=args.radius_std,
        noise_lambda=args.noise_lambda,
        axis_ratio_min=args.axis_ratio_min,
        intensity_std=args.intensity_std,
        seed=args.seed,
    )

    logger.info("Done!")


if __name__ == "__main__":
    main()
