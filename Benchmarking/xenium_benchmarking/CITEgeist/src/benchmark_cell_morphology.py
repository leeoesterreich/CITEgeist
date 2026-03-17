#!/usr/bin/env python
"""
Benchmark cell morphology features for single-cell assignment.

Compares nuclear-only features vs cell-level features (via watershed)
for Module 3b cell type assignment accuracy.

This script:
1. Loads DAPI and boundary channels from Xenium RCC data
2. Runs StarDist on DAPI to get nuclear segmentation
3. Runs watershed using boundary channel to expand to cell boundaries
4. Extracts cell features (nuclear + cytoplasmic morphology)
5. Saves outputs for evaluation

Usage:
    python benchmark_cell_morphology.py --region 0 --output-dir ./output/cell_morphology
"""

import argparse
import json
import logging
import sys
import time
from pathlib import Path
from typing import Any, Dict, Optional, Tuple

import numpy as np
import pandas as pd
import tifffile

# Add CITEgeist to path
REPO_ROOT = Path(__file__).parent.parent.parent.parent.parent
sys.path.insert(0, str(REPO_ROOT))

from CITEgeist.model.segmentation import run_nuclei_segmentation
from CITEgeist.model.watershed_segmentation import watershed_from_nuclei
from CITEgeist.model.morphology_features import extract_cell_features

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)

# =============================================================================
# DEFAULT PATHS
# =============================================================================

XENIUM_DIR = Path(
    "/ix1/alee/LO_LAB/General/Public_Data/10x_PublicData/"
    "Xenium_RNA_Proteomic_RenalCellCarcinoma"
)
MORPHOLOGY_DIR = XENIUM_DIR / "morphology_focus"
IMAGE_DIR = REPO_ROOT / "Benchmarking/xenium_benchmarking/scResolve/images/morphology_hires"

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


# =============================================================================
# IMAGE LOADING FUNCTIONS
# =============================================================================


def micron_to_pixel(x_um: float, y_um: float) -> Tuple[int, int]:
    """Convert micron coordinates to pixel coordinates."""
    x_px = int(x_um / PIXEL_SIZE_UM)
    y_px = int(y_um / PIXEL_SIZE_UM)
    return x_px, y_px


def get_region_bounds_pixel(region_id: int) -> Tuple[int, int, int, int]:
    """
    Get pixel bounds for a region with padding.

    Returns:
        Tuple of (x_min_px, x_max_px, y_min_px, y_max_px)
    """
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


def load_channel_region(
    tiff_path: Path,
    region_id: int,
) -> np.ndarray:
    """
    Load a region from an OME-TIFF channel image.

    Args:
        tiff_path: Path to OME-TIFF file
        region_id: Region ID (0-4)

    Returns:
        2D numpy array of channel intensities
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


def load_dapi_region(region_id: int) -> np.ndarray:
    """Load DAPI channel for a region."""
    dapi_path = MORPHOLOGY_DIR / "ch0000_dapi.ome.tif"
    if not dapi_path.exists():
        raise FileNotFoundError(f"DAPI channel not found: {dapi_path}")
    return load_channel_region(dapi_path, region_id)


def load_boundary_region(region_id: int) -> np.ndarray:
    """Load boundary channel (ATP1A1/CD45/E-cadherin) for a region."""
    boundary_path = MORPHOLOGY_DIR / "ch0001_atp1a1_cd45_e-cadherin.ome.tif"
    if not boundary_path.exists():
        raise FileNotFoundError(f"Boundary channel not found: {boundary_path}")
    return load_channel_region(boundary_path, region_id)


def load_coord_info(region_id: int) -> Dict[str, Any]:
    """Load coordinate info JSON for a region."""
    region_name = f"Xenium_region_{region_id}"
    json_path = IMAGE_DIR / region_name / "coord_info.json"

    if not json_path.exists():
        # Generate coord_info from region bounds
        x_min, x_max, y_min, y_max = get_region_bounds_pixel(region_id)
        x_min_um, x_max_um, y_min_um, y_max_um = REGION_BOUNDS_UM[region_id]
        return {
            "region_id": region_id,
            "micron_bounds": {
                "x_min": x_min_um - PADDING_UM,
                "x_max": x_max_um + PADDING_UM,
                "y_min": y_min_um - PADDING_UM,
                "y_max": y_max_um + PADDING_UM,
            },
            "pixel_bounds": {
                "x_min": x_min,
                "x_max": x_max,
                "y_min": y_min,
                "y_max": y_max,
            },
            "pixel_size": PIXEL_SIZE_UM,
            "padding_microns": PADDING_UM,
        }

    with open(json_path) as f:
        return json.load(f)


def normalize_to_uint8(img: np.ndarray) -> np.ndarray:
    """Normalize image to uint8 range [0, 255]."""
    img_min = img.min()
    img_max = img.max()
    if img_max > img_min:
        img_norm = (img - img_min) / (img_max - img_min) * 255
    else:
        img_norm = np.zeros_like(img)
    return img_norm.astype(np.uint8)


# =============================================================================
# MAIN BENCHMARK FUNCTION
# =============================================================================


def run_benchmark(
    region_id: int,
    output_dir: Path,
) -> Dict[str, Any]:
    """
    Run cell morphology benchmark for a single region.

    Args:
        region_id: Region ID (0-4)
        output_dir: Output directory for results

    Returns:
        Dictionary with benchmark results
    """
    region_name = f"Xenium_region_{region_id}"
    logger.info("=" * 70)
    logger.info("CELL MORPHOLOGY BENCHMARK: %s", region_name)
    logger.info("=" * 70)

    result_dir = Path(output_dir) / region_name
    result_dir.mkdir(parents=True, exist_ok=True)

    timings = {}

    # =========================================================================
    # Step 1: Load DAPI and boundary channels
    # =========================================================================
    logger.info("-" * 60)
    logger.info("STEP 1: Loading image channels")
    logger.info("-" * 60)

    t0 = time.time()
    dapi_img = load_dapi_region(region_id)
    boundary_img = load_boundary_region(region_id)
    timings["load_images_sec"] = time.time() - t0

    logger.info("DAPI shape: %s, dtype: %s", dapi_img.shape, dapi_img.dtype)
    logger.info("Boundary shape: %s, dtype: %s", boundary_img.shape, boundary_img.dtype)

    if dapi_img.shape != boundary_img.shape:
        raise ValueError(
            f"Shape mismatch: DAPI {dapi_img.shape} vs boundary {boundary_img.shape}"
        )

    # =========================================================================
    # Step 2: StarDist nuclear segmentation
    # =========================================================================
    logger.info("-" * 60)
    logger.info("STEP 2: StarDist nuclear segmentation")
    logger.info("-" * 60)

    # Convert DAPI to uint8 for segmentation
    dapi_uint8 = normalize_to_uint8(dapi_img)

    logger.info("Running StarDist nuclei segmentation")

    t0 = time.time()
    nucleus_mask, centroids_df = run_nuclei_segmentation(dapi_uint8, modality="dapi")
    centroids = centroids_df[["x_pixel", "y_pixel"]].values
    timings["segmentation_sec"] = time.time() - t0
    n_nuclei = len(centroids)

    logger.info(
        "StarDist: %d nuclei detected in %.1fs (%.1f nuclei/sec)",
        n_nuclei, timings["segmentation_sec"],
        n_nuclei / timings["segmentation_sec"] if timings["segmentation_sec"] > 0 else 0,
    )

    # =========================================================================
    # Step 3: Watershed cell segmentation
    # =========================================================================
    logger.info("-" * 60)
    logger.info("STEP 3: Watershed cell segmentation")
    logger.info("-" * 60)

    t0 = time.time()
    cell_mask = watershed_from_nuclei(nucleus_mask, boundary_img, use_gradient=True)
    timings["watershed_sec"] = time.time() - t0

    n_cells = len(np.unique(cell_mask)) - 1  # Exclude background (0)
    logger.info(
        "Watershed: %d cells in %.1fs",
        n_cells, timings["watershed_sec"]
    )

    # =========================================================================
    # Step 4: Extract morphology features
    # =========================================================================
    logger.info("-" * 60)
    logger.info("STEP 4: Extracting morphology features")
    logger.info("-" * 60)

    t0 = time.time()
    features_df = extract_cell_features(nucleus_mask, cell_mask)
    timings["feature_extraction_sec"] = time.time() - t0

    logger.info(
        "Extracted features for %d cells (%d columns) in %.1fs",
        len(features_df), len(features_df.columns), timings["feature_extraction_sec"]
    )

    # Compute feature summary statistics
    feature_stats = {}
    numeric_cols = features_df.select_dtypes(include=[np.number]).columns
    for col in numeric_cols:
        feature_stats[col] = {
            "mean": float(features_df[col].mean()),
            "std": float(features_df[col].std()),
            "min": float(features_df[col].min()),
            "max": float(features_df[col].max()),
            "median": float(features_df[col].median()),
        }

    # =========================================================================
    # Step 5: Save outputs
    # =========================================================================
    logger.info("-" * 60)
    logger.info("STEP 5: Saving outputs")
    logger.info("-" * 60)

    # Save features
    features_path = result_dir / "cell_features.csv"
    features_df.to_csv(features_path, index=False)
    logger.info("Saved features to %s", features_path)

    # Save masks as compressed numpy arrays
    nucleus_mask_path = result_dir / "nucleus_mask.npz"
    np.savez_compressed(nucleus_mask_path, mask=nucleus_mask)
    logger.info("Saved nucleus mask to %s", nucleus_mask_path)

    cell_mask_path = result_dir / "cell_mask.npz"
    np.savez_compressed(cell_mask_path, mask=cell_mask)
    logger.info("Saved cell mask to %s", cell_mask_path)

    # Save centroids
    centroids_df = pd.DataFrame(centroids, columns=["x", "y"])
    centroids_path = result_dir / "nucleus_centroids.csv"
    centroids_df.to_csv(centroids_path, index=False)
    logger.info("Saved %d centroids to %s", len(centroids_df), centroids_path)

    # =========================================================================
    # Compile results
    # =========================================================================
    results = {
        "region": region_name,
        "region_id": region_id,
        "image_shape": list(dapi_img.shape),
        "n_nuclei": n_nuclei,
        "n_cells": n_cells,
        "n_cells_with_features": len(features_df),
        "timings": timings,
        "total_time_sec": sum(timings.values()),
        "feature_columns": list(features_df.columns),
        "feature_stats": feature_stats,
        "parameters": {
            "segmentation_backend": "stardist",
            "modality": "dapi",
            "use_gradient": True,
        },
    }

    results_path = result_dir / "benchmark_results.json"
    with open(results_path, "w") as f:
        json.dump(results, f, indent=2)
    logger.info("Saved results to %s", results_path)

    return results


def main():
    parser = argparse.ArgumentParser(
        description="Benchmark cell morphology features for single-cell assignment"
    )
    parser.add_argument(
        "--region", type=int, required=True, choices=[0, 1, 2, 3, 4],
        help="Region ID (0-4)"
    )
    parser.add_argument(
        "--output-dir", type=str,
        default=str(REPO_ROOT / "Benchmarking/xenium_benchmarking/CITEgeist/output/cell_morphology"),
        help="Output directory for results"
    )

    args = parser.parse_args()

    results = run_benchmark(
        region_id=args.region,
        output_dir=Path(args.output_dir),
    )

    logger.info("=" * 70)
    logger.info("BENCHMARK COMPLETE")
    logger.info("=" * 70)
    logger.info("Region: %s", results["region"])
    logger.info("Nuclei detected: %d", results["n_nuclei"])
    logger.info("Cells with features: %d", results["n_cells_with_features"])
    logger.info("Total time: %.1fs", results["total_time_sec"])
    logger.info("Results saved to: %s", args.output_dir)


if __name__ == "__main__":
    main()
