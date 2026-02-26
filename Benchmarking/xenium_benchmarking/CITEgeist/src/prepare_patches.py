#!/usr/bin/env python
"""
Prepare nucleus patches for VAE training.

Extracts patches around detected nuclei from multi-channel images,
organizing them by spot for downstream Sinkhorn assignment.

Usage:
    python prepare_patches.py \
        --image morphology_mip.tiff \
        --mask cellpose_nuclei.npy \
        --nuclei_spot_map nuclei_spot_mapping.csv \
        --output_dir output/patches \
        --expansion 0.75 \
        --patch_size 96
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

from CITEgeist.model.patch_extraction import extract_patch
from CITEgeist.model.morphology_features import extract_nucleus_features

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
    image_path: str,
    mask_path: str,
    nuclei_spot_map_path: str,
    output_dir: str,
    expansion: float = 0.75,
    patch_size: int = 96,
    min_bbox_size: int = 4,
) -> None:
    """
    Extract and save patches for all nuclei, organized by spot.

    Args:
        image_path: Path to multi-channel TIFF image
        mask_path: Path to Cellpose nucleus mask (.npy or .tiff)
        nuclei_spot_map_path: CSV mapping nucleus_id to spot_id
        output_dir: Output directory for patches
        expansion: Bbox expansion fraction (default 0.75)
        patch_size: Output patch size (default 96)
        min_bbox_size: Minimum bbox dimension to extract (default 4)
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load image
    image = load_image(Path(image_path))

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
    }

    for spot_id in tqdm(spot_ids, desc="Processing spots"):
        spot_nuclei = features_df[features_df["spot_id"] == spot_id]

        patches = []
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
                patch = extract_patch(image, bbox, expansion, patch_size)
                patches.append(patch)
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
            np.save(output_dir / f"spot_{spot_id}_patches.npy", patches_array)

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
Example:
    python prepare_patches.py \\
        --image morphology_mip.tiff \\
        --mask cellpose_nuclei.npy \\
        --nuclei_spot_map nuclei_spot_mapping.csv \\
        --output_dir output/patches \\
        --expansion 0.75 \\
        --patch_size 96
        """,
    )
    parser.add_argument(
        "--image",
        required=True,
        help="Path to multi-channel TIFF image",
    )
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

    args = parser.parse_args()

    prepare_patches(
        image_path=args.image,
        mask_path=args.mask,
        nuclei_spot_map_path=args.nuclei_spot_map,
        output_dir=args.output_dir,
        expansion=args.expansion,
        patch_size=args.patch_size,
        min_bbox_size=args.min_bbox_size,
    )


if __name__ == "__main__":
    main()
