#!/usr/bin/env python
"""
Prepare nucleus patches for VAE training.

Extracts patches around detected nuclei from multi-channel images,
organizing them by spot for downstream Sinkhorn assignment.

Supports two modes:
1. Single image mode (--image): Load a single multi-channel image (PNG/TIFF)
2. OME-TIFF mode (--dapi + --boundary + --region): Load DAPI and boundary
   channels from separate OME-TIFF files and crop to a specific region

Usage (single image):
    python prepare_patches.py \
        --image morphology_mip.tiff \
        --mask cellpose_nuclei.npy \
        --nuclei_spot_map nuclei_spot_mapping.csv \
        --output_dir output/patches

Usage (OME-TIFF dual channel):
    python prepare_patches.py \
        --dapi ch0000_dapi.ome.tif \
        --boundary ch0001_atp1a1_cd45_e-cadherin.ome.tif \
        --region 0 \
        --mask cellpose_nuclei.npy \
        --nuclei_spot_map nuclei_spot_mapping.csv \
        --output_dir output/patches
"""
from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
import tifffile
from skimage.measure import regionprops
from tqdm import tqdm

# Add repository root to path for CITEgeist imports
REPO_ROOT = Path(__file__).parent.parent.parent.parent.parent
sys.path.insert(0, str(REPO_ROOT))

from CITEgeist.model.patch_extraction import (
    compute_global_stats,
    extract_patch_with_size,
)
from CITEgeist.model.morphology_features import extract_nucleus_features

# =============================================================================
# OME-TIFF REGION LOADING (for Xenium benchmark)
# =============================================================================

# Pixel size for Xenium morphology images
PIXEL_SIZE_UM = 0.2125

# Region boundaries in microns (from pseudo-Visium data)
# Format: (x_min, x_max, y_min, y_max)
REGION_BOUNDS_UM = {
    0: (29.01, 2279.01, 30.87, 5486.83),
    1: (2329.01, 4579.01, 30.87, 5400.23),
    2: (4629.01, 6829.01, 30.87, 5486.83),
    3: (6879.01, 9129.01, 30.87, 5660.04),
    4: (9179.01, 11429.01, 30.87, 5746.64),
}

PADDING_UM = 100.0  # Padding around region in microns


def micron_to_pixel(x_um: float, y_um: float) -> tuple:
    """Convert micron coordinates to pixel coordinates."""
    x_px = int(x_um / PIXEL_SIZE_UM)
    y_px = int(y_um / PIXEL_SIZE_UM)
    return x_px, y_px


def get_region_bounds_pixel(region_id: int) -> tuple:
    """Get pixel bounds for a region with padding."""
    x_min_um, x_max_um, y_min_um, y_max_um = REGION_BOUNDS_UM[region_id]

    # Add padding
    x_min_um -= PADDING_UM
    x_max_um += PADDING_UM
    y_min_um -= PADDING_UM
    y_max_um += PADDING_UM

    x_min_px, y_min_px = micron_to_pixel(x_min_um, y_min_um)
    x_max_px, y_max_px = micron_to_pixel(x_max_um, y_max_um)

    # Ensure non-negative
    x_min_px = max(0, x_min_px)
    y_min_px = max(0, y_min_px)

    return x_min_px, x_max_px, y_min_px, y_max_px


def load_ome_tiff_region(tiff_path: Path, region_id: int) -> np.ndarray:
    """
    Load a region from an OME-TIFF channel image.

    Args:
        tiff_path: Path to OME-TIFF file
        region_id: Region ID (0-4)

    Returns:
        2D numpy array of channel intensities (H, W)
    """
    x_min, x_max, y_min, y_max = get_region_bounds_pixel(region_id)

    logger.info(
        "Loading region %d from %s (y=[%d:%d], x=[%d:%d])",
        region_id, tiff_path.name, y_min, y_max, x_min, x_max
    )

    with tifffile.TiffFile(tiff_path) as tif:
        page = tif.pages[0]
        # OME-TIFF reads as (Y, X)
        region = page.asarray()[y_min:y_max, x_min:x_max]

    return region.astype(np.float32)


def load_dual_channel_image(
    dapi_path: Path,
    boundary_path: Path,
    region_id: int,
) -> np.ndarray:
    """
    Load DAPI and boundary channels and stack as 2-channel image.

    Args:
        dapi_path: Path to DAPI OME-TIFF
        boundary_path: Path to boundary OME-TIFF
        region_id: Region ID (0-4)

    Returns:
        Image array with shape (2, H, W), normalized to [0, 1]
    """
    logger.info("Loading dual-channel image (DAPI + boundary) for region %d", region_id)

    # Load both channels
    dapi = load_ome_tiff_region(dapi_path, region_id)
    boundary = load_ome_tiff_region(boundary_path, region_id)

    # Verify shapes match
    if dapi.shape != boundary.shape:
        raise ValueError(
            f"Channel shape mismatch: DAPI {dapi.shape} vs boundary {boundary.shape}"
        )

    # Normalize each channel to [0, 1] (uint16 data)
    dapi_norm = dapi / 65535.0
    boundary_norm = boundary / 65535.0

    # Stack as (2, H, W)
    image = np.stack([dapi_norm, boundary_norm], axis=0)

    logger.info(
        "Loaded dual-channel image: shape=%s, DAPI range=[%.3f, %.3f], "
        "boundary range=[%.3f, %.3f]",
        image.shape,
        dapi_norm.min(), dapi_norm.max(),
        boundary_norm.min(), boundary_norm.max(),
    )

    return image


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)


def load_image(image_path: Path) -> np.ndarray:
    """
    Load image and ensure (C, H, W) format.

    Args:
        image_path: Path to multi-channel TIFF or PNG image

    Returns:
        Image array with shape (C, H, W)
    """
    logger.info(f"Loading image from {image_path}")
    suffix = image_path.suffix.lower()

    if suffix in [".png", ".jpg", ".jpeg"]:
        # Use PIL for common image formats
        from PIL import Image
        # Disable decompression bomb check for large microscopy images
        Image.MAX_IMAGE_PIXELS = None
        img = Image.open(str(image_path))
        image = np.array(img)
        logger.info(f"Loaded PNG/JPEG image with shape {image.shape}")
    else:
        # Use tifffile for TIFF
        image = tifffile.imread(str(image_path))

    # Handle different input shapes
    if image.ndim == 2:
        # Single channel - add channel dimension
        image = image[np.newaxis, ...]
        logger.info(f"Added channel dimension: {image.shape}")
    elif image.ndim == 3:
        # Could be (C, H, W) or (H, W, C)
        # Assume (H, W, C) if last dim is small (common for RGB/multi-channel)
        if image.shape[-1] <= 4 and image.shape[-1] < image.shape[0]:
            image = np.moveaxis(image, -1, 0)
            logger.info(f"Transposed from (H, W, C) to (C, H, W): {image.shape}")
    elif image.ndim == 4:
        # Likely (1, C, H, W) or (T, H, W, C) - take first slice
        logger.warning(f"4D image detected with shape {image.shape}, taking first slice")
        image = image[0]
        if image.ndim == 3 and image.shape[-1] <= 4:
            image = np.moveaxis(image, -1, 0)

    # Convert to float32 for numerical stability
    if image.dtype == np.uint8:
        image = image.astype(np.float32) / 255.0
        logger.info("Converted uint8 to float32 [0,1]")
    elif image.dtype != np.float32:
        image = image.astype(np.float32)

    logger.info(f"Loaded image with shape {image.shape}, dtype {image.dtype}")
    return image


def load_mask(mask_path: Path) -> np.ndarray:
    """
    Load nucleus mask from .npy or .tiff file.

    Args:
        mask_path: Path to mask file

    Returns:
        2D labeled mask array
    """
    logger.info(f"Loading mask from {mask_path}")

    suffix = mask_path.suffix.lower()
    if suffix == ".npy":
        mask = np.load(str(mask_path))
    elif suffix in [".tiff", ".tif"]:
        mask = tifffile.imread(str(mask_path))
    else:
        raise ValueError(f"Unsupported mask format: {suffix}")

    # Ensure 2D
    if mask.ndim == 3:
        logger.warning(f"3D mask detected with shape {mask.shape}, taking first slice")
        mask = mask[0]

    logger.info(f"Loaded mask with shape {mask.shape}, {len(np.unique(mask))-1} nuclei")
    return mask


def add_bboxes_to_features(features_df: pd.DataFrame, mask: np.ndarray) -> pd.DataFrame:
    """
    Add bounding box coordinates to features DataFrame.

    Args:
        features_df: DataFrame with nucleus_id column
        mask: Labeled nucleus mask

    Returns:
        DataFrame with bbox_x_min, bbox_y_min, bbox_x_max, bbox_y_max columns
    """
    logger.info("Computing bounding boxes from regionprops...")

    # Get bounding boxes using regionprops
    # regionprops returns bbox as (min_row, min_col, max_row, max_col)
    # which is (y_min, x_min, y_max, x_max)
    props = regionprops(mask)

    bbox_data = {}
    for p in props:
        # Convert from (y_min, x_min, y_max, x_max) to (x_min, y_min, x_max, y_max)
        bbox_data[p.label] = (p.bbox[1], p.bbox[0], p.bbox[3], p.bbox[2])

    # Add columns to DataFrame
    features_df["bbox_x_min"] = features_df["nucleus_id"].map(
        lambda x: bbox_data.get(x, (0, 0, 0, 0))[0]
    )
    features_df["bbox_y_min"] = features_df["nucleus_id"].map(
        lambda x: bbox_data.get(x, (0, 0, 0, 0))[1]
    )
    features_df["bbox_x_max"] = features_df["nucleus_id"].map(
        lambda x: bbox_data.get(x, (0, 0, 0, 0))[2]
    )
    features_df["bbox_y_max"] = features_df["nucleus_id"].map(
        lambda x: bbox_data.get(x, (0, 0, 0, 0))[3]
    )

    return features_df


def prepare_patches(
    mask_path: str,
    nuclei_spot_map_path: str,
    output_dir: str,
    image_path: Optional[str] = None,
    dapi_path: Optional[str] = None,
    boundary_path: Optional[str] = None,
    region_id: Optional[int] = None,
    expansion: float = 0.75,
    patch_size: int = 96,
    min_bbox_size: int = 4,
    norm_method: str = "percentile",
) -> None:
    """
    Extract and save patches for all nuclei, organized by spot.

    Supports two modes:
    1. Single image mode: provide image_path
    2. OME-TIFF dual channel mode: provide dapi_path, boundary_path, region_id

    Args:
        mask_path: Path to Cellpose nucleus mask (.npy or .tiff)
        nuclei_spot_map_path: CSV mapping nucleus_id to spot_id
        output_dir: Output directory for patches
        image_path: Path to multi-channel TIFF/PNG image (mode 1)
        dapi_path: Path to DAPI OME-TIFF (mode 2)
        boundary_path: Path to boundary OME-TIFF (mode 2)
        region_id: Region ID 0-4 for OME-TIFF cropping (mode 2)
        expansion: Bbox expansion fraction (default 0.75)
        patch_size: Output patch size (default 96)
        min_bbox_size: Minimum bbox dimension to extract (default 4)
        norm_method: Normalization method ("percentile" or "minmax")
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load image (two modes)
    if dapi_path is not None and boundary_path is not None and region_id is not None:
        # Mode 2: Dual-channel OME-TIFF
        image = load_dual_channel_image(
            Path(dapi_path), Path(boundary_path), region_id
        )
    elif image_path is not None:
        # Mode 1: Single image
        image = load_image(Path(image_path))
    else:
        raise ValueError(
            "Must provide either --image OR (--dapi, --boundary, --region)"
        )

    # Compute and save global normalization stats
    logger.info(f"Computing global normalization stats (method={norm_method})...")
    global_stats = compute_global_stats(image, norm_method=norm_method)
    stats_path = output_dir / "global_stats.npz"
    np.savez(stats_path, **global_stats)
    logger.info(f"Saved global normalization stats to {stats_path}")

    # Load mask
    mask = load_mask(Path(mask_path))

    # Validate shapes match
    if image.shape[1:] != mask.shape:
        raise ValueError(
            f"Image spatial shape {image.shape[1:]} doesn't match mask shape {mask.shape}"
        )

    # Extract nucleus features
    logger.info("Extracting nucleus features...")
    features_df = extract_nucleus_features(mask)
    logger.info(f"Extracted features for {len(features_df)} nuclei")

    # Add bounding boxes
    features_df = add_bboxes_to_features(features_df, mask)

    # Load spot mapping
    logger.info(f"Loading spot mapping from {nuclei_spot_map_path}")
    spot_map = pd.read_csv(nuclei_spot_map_path)

    # Validate spot_map has required columns
    if "nucleus_id" not in spot_map.columns:
        raise ValueError("nuclei_spot_map CSV must have 'nucleus_id' column")
    if "spot_id" not in spot_map.columns:
        raise ValueError("nuclei_spot_map CSV must have 'spot_id' column")

    # Merge features with spot mapping
    features_df = features_df.merge(spot_map, on="nucleus_id", how="inner")
    logger.info(f"{len(features_df)} nuclei matched to spots")

    if len(features_df) == 0:
        logger.warning("No nuclei matched to spots! Check nucleus_id format.")
        return

    # Extract patches per spot
    spot_ids = features_df["spot_id"].unique()
    logger.info(f"Extracting patches for {len(spot_ids)} spots...")

    stats = {
        "total_spots": len(spot_ids),
        "total_nuclei": len(features_df),
        "successful_patches": 0,
        "failed_patches": 0,
        "empty_spots": 0,
        "norm_method": norm_method,
    }

    for spot_id in tqdm(spot_ids, desc="Processing spots"):
        spot_nuclei = features_df[features_df["spot_id"] == spot_id]

        patches = []
        size_features_list = []
        nucleus_ids = []

        for _, row in spot_nuclei.iterrows():
            bbox = (
                int(row["bbox_x_min"]),
                int(row["bbox_y_min"]),
                int(row["bbox_x_max"]),
                int(row["bbox_y_max"]),
            )

            # Skip very small bboxes
            bbox_w = bbox[2] - bbox[0]
            bbox_h = bbox[3] - bbox[1]
            if bbox_w < min_bbox_size or bbox_h < min_bbox_size:
                logger.debug(
                    f"Skipping nucleus {row['nucleus_id']} - bbox too small: {bbox_w}x{bbox_h}"
                )
                stats["failed_patches"] += 1
                continue

            try:
                patch, size_feats = extract_patch_with_size(
                    image, bbox, expansion, patch_size, global_stats
                )
                patches.append(patch)
                size_features_list.append(size_feats)
                nucleus_ids.append(row["nucleus_id"])
                stats["successful_patches"] += 1
            except Exception as e:
                logger.debug(
                    f"Failed to extract patch for nucleus {row['nucleus_id']}: {e}"
                )
                stats["failed_patches"] += 1
                continue

        if patches:
            # Stack patches: (N, C, H, W)
            patches_array = np.stack(patches, axis=0)
            sizes_array = np.stack(size_features_list, axis=0)

            np.save(output_dir / f"spot_{spot_id}_patches.npy", patches_array)
            np.save(output_dir / f"spot_{spot_id}_sizes.npy", sizes_array)

            # Save nucleus IDs for this spot (for later reference)
            np.save(
                output_dir / f"spot_{spot_id}_nucleus_ids.npy",
                np.array(nucleus_ids, dtype=np.int64),
            )
        else:
            stats["empty_spots"] += 1

    # Save features with spot assignments and bboxes
    features_path = output_dir / "nucleus_features.csv"
    features_df.to_csv(features_path, index=False)
    logger.info(f"Saved nucleus features to {features_path}")

    # Save summary stats
    stats_path = output_dir / "extraction_stats.json"
    with open(stats_path, "w") as f:
        json.dump(stats, f, indent=2)
    logger.info(f"Saved extraction stats to {stats_path}")

    # Validate required outputs exist
    required_base_files = ["global_stats.npz", "nucleus_features.csv"]
    for f in required_base_files:
        if not (output_dir / f).exists():
            raise RuntimeError(f"Missing required output: {f}")

    # Spot-check a few spots (validate patches have matching sizes)
    sample_spots = list(spot_ids)[:3]
    for spot_id in sample_spots:
        patches_file = output_dir / f"spot_{spot_id}_patches.npy"
        sizes_file = output_dir / f"spot_{spot_id}_sizes.npy"
        if patches_file.exists() and not sizes_file.exists():
            raise RuntimeError(f"Patches exist but sizes missing for spot {spot_id}")

    # Print summary
    logger.info("=" * 50)
    logger.info("Patch extraction complete:")
    logger.info(f"  Total spots: {stats['total_spots']}")
    logger.info(f"  Total nuclei: {stats['total_nuclei']}")
    logger.info(f"  Successful patches: {stats['successful_patches']}")
    logger.info(f"  Failed patches: {stats['failed_patches']}")
    logger.info(f"  Empty spots: {stats['empty_spots']}")
    logger.info(f"  Output directory: {output_dir}")


def main():
    parser = argparse.ArgumentParser(
        description="Prepare nucleus patches for VAE training",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example (single image):
    python prepare_patches.py \\
        --image morphology_mip.tiff \\
        --mask cellpose_nuclei.npy \\
        --nuclei_spot_map nuclei_spot_mapping.csv \\
        --output_dir output/patches

Example (OME-TIFF dual channel - RECOMMENDED):
    python prepare_patches.py \\
        --dapi ch0000_dapi.ome.tif \\
        --boundary ch0001_atp1a1_cd45_e-cadherin.ome.tif \\
        --region 0 \\
        --mask cellpose_nuclei.npy \\
        --nuclei_spot_map nuclei_spot_mapping.csv \\
        --output_dir output/patches
        """,
    )
    # Image input (mutually exclusive modes)
    parser.add_argument(
        "--image",
        help="Path to multi-channel TIFF/PNG image (mode 1)",
    )
    parser.add_argument(
        "--dapi",
        help="Path to DAPI OME-TIFF (mode 2, use with --boundary and --region)",
    )
    parser.add_argument(
        "--boundary",
        help="Path to boundary channel OME-TIFF (mode 2)",
    )
    parser.add_argument(
        "--region",
        type=int,
        choices=[0, 1, 2, 3, 4],
        help="Region ID for OME-TIFF cropping (mode 2)",
    )
    # Required arguments
    parser.add_argument(
        "--mask",
        required=True,
        help="Path to Cellpose nucleus mask (.npy or .tiff)",
    )
    parser.add_argument(
        "--nuclei_spot_map",
        required=True,
        help="CSV mapping nucleus_id to spot_id",
    )
    parser.add_argument(
        "--output_dir",
        required=True,
        help="Output directory for patches",
    )
    # Optional arguments
    parser.add_argument(
        "--expansion",
        type=float,
        default=0.75,
        help="Bbox expansion fraction (default: 0.75)",
    )
    parser.add_argument(
        "--patch_size",
        type=int,
        default=96,
        help="Output patch size (default: 96)",
    )
    parser.add_argument(
        "--min_bbox_size",
        type=int,
        default=4,
        help="Minimum bbox dimension to extract (default: 4)",
    )
    parser.add_argument(
        "--norm-method",
        type=str,
        choices=["minmax", "percentile"],
        default="percentile",
        help="Normalization method: minmax (dtype range) or percentile (1st/99th)",
    )

    args = parser.parse_args()

    # Validate arguments
    if args.image is None and (args.dapi is None or args.boundary is None or args.region is None):
        parser.error("Must provide either --image OR (--dapi, --boundary, --region)")

    prepare_patches(
        mask_path=args.mask,
        nuclei_spot_map_path=args.nuclei_spot_map,
        output_dir=args.output_dir,
        image_path=args.image,
        dapi_path=args.dapi,
        boundary_path=args.boundary,
        region_id=args.region,
        expansion=args.expansion,
        patch_size=args.patch_size,
        min_bbox_size=args.min_bbox_size,
        norm_method=args.norm_method,
    )


if __name__ == "__main__":
    main()
