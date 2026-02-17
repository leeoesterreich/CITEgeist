#!/usr/bin/env python
"""
Discrete cell assignment benchmark with Cellpose nuclei segmentation.

This script runs the full CITEgeist discrete pipeline:
1. Load morphology image for a Xenium pseudo-Visium region
2. Run Cellpose segmentation to get nuclei counts per spot
3. Run discrete cell assignment (IQP) using nuclei counts as constraints
4. Run GEX deconvolution using discrete cell counts
5. Save outputs for evaluation

Usage:
    python benchmark_discrete_cellpose.py --region 0 --output-dir ./output/discrete
    python benchmark_discrete_cellpose.py --region 0 --use-gpu --cellpose-diameter 20
"""

import argparse
import json
import logging
import sys
import time
from pathlib import Path
from typing import Any, Dict, Optional

import cv2
import numpy as np
import pandas as pd
import scanpy as sc

# Add CITEgeist to path
REPO_ROOT = Path(__file__).parent.parent.parent.parent.parent
sys.path.insert(0, str(REPO_ROOT))

# Import shared benchmark constants
BENCHMARK_ROOT = Path(__file__).parent.parent.parent
sys.path.insert(0, str(BENCHMARK_ROOT))
from benchmark_constants import ACHIEVABLE_7_CELL_PROFILE_DICT

from CITEgeist.model.citegeist_model import CitegeistModel
from CITEgeist.model.segmentation import (
    SegmentationResult,
    assign_nuclei_centroids_to_spots,
    run_cellpose_nuclei_segmentation,
)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)


# =============================================================================
# DEFAULT PATHS
# =============================================================================

DATA_DIR = REPO_ROOT / "Benchmarking/xenium_pseudovisium/data_protein_gt"
IMAGE_DIR = REPO_ROOT / "Benchmarking/xenium_benchmarking/scResolve/images/morphology_hires"
N_REGIONS = 5  # Xenium_region_0 through Xenium_region_4


# =============================================================================
# IMAGE AND COORDINATE LOADING
# =============================================================================


def load_morphology_image(region_id: int, image_dir: Path) -> np.ndarray:
    """
    Load morphology PNG image for a region.

    Args:
        region_id: Xenium region ID (0-4)
        image_dir: Root directory containing region folders

    Returns:
        RGB uint8 numpy array
    """
    region_name = f"Xenium_region_{region_id}"
    image_path = image_dir / region_name / "morphology.png"

    if not image_path.exists():
        raise FileNotFoundError(f"Missing morphology image: {image_path}")

    logger.info("Loading morphology image: %s", image_path)
    bgr = cv2.imread(str(image_path), cv2.IMREAD_COLOR)
    if bgr is None:
        raise ValueError(f"Failed to load image: {image_path}")

    rgb = cv2.cvtColor(bgr, cv2.COLOR_BGR2RGB)
    logger.info("Image shape: %s, dtype: %s", rgb.shape, rgb.dtype)
    return rgb


def load_coord_info(region_id: int, image_dir: Path) -> Dict[str, Any]:
    """
    Load coordinate info JSON for a region.

    Args:
        region_id: Xenium region ID (0-4)
        image_dir: Root directory containing region folders

    Returns:
        Dictionary with pixel_size, micron_bounds, pixel_bounds, image_size
    """
    region_name = f"Xenium_region_{region_id}"
    coord_path = image_dir / region_name / "coord_info.json"

    if not coord_path.exists():
        raise FileNotFoundError(f"Missing coord_info.json: {coord_path}")

    with open(coord_path) as f:
        coord_info = json.load(f)

    logger.info(
        "Coord info: pixel_size=%.4f um/px, image_size=%s",
        coord_info["pixel_size"],
        coord_info["image_size"],
    )
    return coord_info


def convert_micron_to_pixel(
    spot_coords_micron: np.ndarray,
    coord_info: Dict[str, Any],
) -> np.ndarray:
    """
    Convert spot coordinates from microns to pixels in the morphology image frame.

    Uses pixel_bounds to handle image cropping and padding correctly. The pixel_bounds
    field stores where the image crop starts in the full Xenium image coordinate system.

    Args:
        spot_coords_micron: Array of shape (N, 2) with coordinates in microns
        coord_info: Coordinate info dictionary with pixel_size and pixel_bounds

    Returns:
        Array of shape (N, 2) with coordinates in pixels
    """
    pixel_size = float(coord_info["pixel_size"])

    # Get crop origin in full image pixel coordinates
    px_x_min = float(coord_info["pixel_bounds"]["x_min"])
    px_y_min = float(coord_info["pixel_bounds"]["y_min"])

    # Convert: first to full image pixels, then subtract crop origin
    # pixel_in_crop = micron / pixel_size - crop_origin
    spot_coords_pixel = np.empty_like(spot_coords_micron, dtype=np.float64)
    spot_coords_pixel[:, 0] = spot_coords_micron[:, 0] / pixel_size - px_x_min
    spot_coords_pixel[:, 1] = spot_coords_micron[:, 1] / pixel_size - px_y_min

    return spot_coords_pixel


# =============================================================================
# NUCLEI SEGMENTATION AND ASSIGNMENT
# =============================================================================


def run_cellpose_and_assign(
    image_rgb: np.ndarray,
    spot_coords_pixel: np.ndarray,
    spot_names: pd.Index,
    spot_diameter_um: float,
    pixel_size_um: float,
    use_gpu: bool = False,
    cellpose_diameter: Optional[float] = None,
) -> pd.Series:
    """
    Run Cellpose segmentation and assign nuclei to spots.

    Args:
        image_rgb: RGB morphology image
        spot_coords_pixel: Spot centers in pixel coordinates (N, 2)
        spot_names: Spot identifiers
        spot_diameter_um: Spot diameter in microns (e.g., 55.0 for Visium)
        pixel_size_um: Microns per pixel in image
        use_gpu: Whether to use GPU for Cellpose
        cellpose_diameter: Cellpose nucleus diameter (auto-detect if None)

    Returns:
        Series with nuclei counts per spot
    """
    logger.info("Running Cellpose nuclei segmentation...")
    t0 = time.time()

    masks, centroids_xy = run_cellpose_nuclei_segmentation(
        image_rgb_uint8=image_rgb,
        use_gpu=use_gpu,
        diameter=cellpose_diameter,
        flow_threshold=0.4,
        cellprob_threshold=0.0,
    )

    seg_time = time.time() - t0
    n_nuclei = int(masks.max()) if masks.size > 0 else 0
    logger.info(
        "Cellpose found %d nuclei in %.1fs (%.1f nuclei/sec)",
        n_nuclei,
        seg_time,
        n_nuclei / seg_time if seg_time > 0 else 0,
    )

    # Compute spot radius in pixels
    spot_radius_px = (spot_diameter_um / pixel_size_um) / 2.0
    logger.info(
        "Spot diameter: %.1f um = %.1f pixels (radius=%.1f px)",
        spot_diameter_um,
        spot_diameter_um / pixel_size_um,
        spot_radius_px,
    )

    # Assign nuclei centroids to spots
    nuclei_counts = assign_nuclei_centroids_to_spots(
        centroids_xy=centroids_xy,
        spot_centers_xy=spot_coords_pixel,
        spot_radius_px=spot_radius_px,
        spot_names=spot_names,
    )

    # Log distribution
    total_assigned = int(nuclei_counts.sum())
    pct_assigned = 100.0 * total_assigned / n_nuclei if n_nuclei > 0 else 0
    logger.info(
        "Assigned %d/%d nuclei to spots (%.1f%%)",
        total_assigned,
        n_nuclei,
        pct_assigned,
    )
    logger.info(
        "Nuclei per spot: min=%d, max=%d, mean=%.1f, median=%.1f",
        nuclei_counts.min(),
        nuclei_counts.max(),
        nuclei_counts.mean(),
        nuclei_counts.median(),
    )

    return nuclei_counts


# =============================================================================
# MAIN BENCHMARK FUNCTION
# =============================================================================


def run_discrete_benchmark(
    region_id: int,
    output_dir: Path,
    data_dir: Path = DATA_DIR,
    image_dir: Path = IMAGE_DIR,
    use_gpu: bool = False,
    cellpose_diameter: Optional[float] = None,
    max_em_iterations: int = 20,
    run_gex: bool = True,
    min_counts: int = 25,
    spot_diameter_um: float = 55.0,
    scale_mode: str = "per_marker",
) -> Dict[str, Any]:
    """
    Run full discrete CITEgeist benchmark for one region.

    Args:
        region_id: Xenium region ID (0-4)
        output_dir: Output directory for results
        data_dir: Directory containing h5ad_objects/
        image_dir: Directory containing morphology images
        use_gpu: Whether to use GPU for Cellpose
        cellpose_diameter: Cellpose nucleus diameter (auto if None)
        max_em_iterations: Maximum EM iterations for discrete assignment
        run_gex: Whether to run GEX deconvolution
        min_counts: Minimum counts filter for GEX preprocessing
        spot_diameter_um: Spot diameter in microns

    Returns:
        Dictionary with benchmark results and timing
    """
    sample_name = f"Xenium_region_{region_id}"
    result_dir = output_dir / sample_name
    result_dir.mkdir(parents=True, exist_ok=True)

    logger.info("=" * 70)
    logger.info("DISCRETE CELLPOSE BENCHMARK: %s", sample_name)
    logger.info("=" * 70)

    timings = {}

    # =========================================================================
    # Step 1: Load data
    # =========================================================================
    logger.info("-" * 60)
    logger.info("STEP 1: Loading data")
    logger.info("-" * 60)

    gex_path = data_dir / "h5ad_objects" / f"{sample_name}_GEX.h5ad"
    cite_path = data_dir / "h5ad_objects" / f"{sample_name}_CITE.h5ad"

    if not gex_path.exists() or not cite_path.exists():
        raise FileNotFoundError(f"Missing h5ad files for region {region_id}")

    gex_adata = sc.read_h5ad(gex_path)
    cite_adata = sc.read_h5ad(cite_path)

    logger.info("GEX: %d spots x %d genes", gex_adata.n_obs, gex_adata.n_vars)
    logger.info("CITE: %d spots x %d proteins", cite_adata.n_obs, cite_adata.n_vars)

    # =========================================================================
    # Step 2: Load morphology image and coordinate info
    # =========================================================================
    logger.info("-" * 60)
    logger.info("STEP 2: Loading morphology image")
    logger.info("-" * 60)

    t0 = time.time()
    image_rgb = load_morphology_image(region_id, image_dir)
    coord_info = load_coord_info(region_id, image_dir)
    timings["image_load_sec"] = time.time() - t0

    # Convert spot coordinates from microns to pixels
    spot_coords_micron = np.asarray(cite_adata.obsm["spatial"], dtype=np.float64)
    spot_coords_pixel = convert_micron_to_pixel(spot_coords_micron, coord_info)
    pixel_size_um = float(coord_info["pixel_size"])

    # =========================================================================
    # Step 3: Run Cellpose segmentation
    # =========================================================================
    logger.info("-" * 60)
    logger.info("STEP 3: Cellpose nuclei segmentation")
    logger.info("-" * 60)

    t0 = time.time()
    nuclei_counts = run_cellpose_and_assign(
        image_rgb=image_rgb,
        spot_coords_pixel=spot_coords_pixel,
        spot_names=cite_adata.obs_names,
        spot_diameter_um=spot_diameter_um,
        pixel_size_um=pixel_size_um,
        use_gpu=use_gpu,
        cellpose_diameter=cellpose_diameter,
    )
    timings["cellpose_sec"] = time.time() - t0

    # Save nuclei counts
    nuclei_counts_path = result_dir / f"{sample_name}_nuclei_counts.csv"
    nuclei_counts.to_csv(nuclei_counts_path, header=True)
    logger.info("Saved nuclei counts: %s", nuclei_counts_path)

    # =========================================================================
    # Step 4: Initialize CITEgeist model
    # =========================================================================
    logger.info("-" * 60)
    logger.info("STEP 4: Initialize CITEgeist model")
    logger.info("-" * 60)

    model = CitegeistModel(
        sample_name=sample_name,
        output_folder=str(result_dir),
        simulation=True,
        gene_expression_adata=gex_adata.copy(),
        antibody_capture_adata=cite_adata.copy(),
    )

    # Load achievable-7 cell profile dictionary
    model.load_cell_profile_dict(ACHIEVABLE_7_CELL_PROFILE_DICT)
    logger.info("Loaded ACHIEVABLE_7 profile with %d cell types", len(ACHIEVABLE_7_CELL_PROFILE_DICT))

    # =========================================================================
    # Step 5: Preprocess data
    # =========================================================================
    logger.info("-" * 60)
    logger.info("STEP 5: Preprocessing")
    logger.info("-" * 60)

    t0 = time.time()
    model.filter_gex(min_counts=min_counts)
    model.preprocess_gex(target_sum=10000)
    # Use discrete-specific preprocessing (preserves cellularity signal)
    model.preprocess_antibody_discrete(
        winsorize_lower=5,
        winsorize_upper=95,
        scale_mode=scale_mode,
    )
    timings["preprocess_sec"] = time.time() - t0

    # =========================================================================
    # Step 6: Run discrete cell assignment (IQP)
    # =========================================================================
    logger.info("-" * 60)
    logger.info("STEP 6: Discrete cell assignment (IQP)")
    logger.info("-" * 60)

    t0 = time.time()
    cell_counts_df = model.run_discrete_cell_assignment(
        nuclei_counts=nuclei_counts,
        max_em_iterations=max_em_iterations,
        beta_convergence_tol=1e-3,
        max_nuclei_cap=30,
        beta_min=0.1,
        beta_max=2.0,
        timeout_per_spot=60.0,
    )
    timings["discrete_assignment_sec"] = time.time() - t0

    # Log cell type distribution
    logger.info("Cell type distribution:")
    total_cells = cell_counts_df.values.sum()
    for ct in cell_counts_df.columns:
        ct_total = cell_counts_df[ct].sum()
        pct = 100 * ct_total / total_cells if total_cells > 0 else 0
        logger.info("  %s: %d cells (%.1f%%)", ct, ct_total, pct)

    # Save cell counts prediction for evaluation
    prediction_path = result_dir / f"{sample_name}_deconv_predictions.csv"
    # For evaluation compatibility, save as proportions
    row_sums = cell_counts_df.values.sum(axis=1, keepdims=True)
    row_sums = np.maximum(row_sums, 1)
    proportions_df = pd.DataFrame(
        cell_counts_df.values / row_sums,
        index=cell_counts_df.index,
        columns=cell_counts_df.columns,
    )
    proportions_df.to_csv(prediction_path)
    logger.info("Saved proportions: %s", prediction_path)

    # =========================================================================
    # Step 7: GEX deconvolution (optional)
    # =========================================================================
    if run_gex:
        logger.info("-" * 60)
        logger.info("STEP 7: Gene expression deconvolution")
        logger.info("-" * 60)

        t0 = time.time()
        try:
            model.run_cell_expression_pass1(
                radius=None,  # Auto-detect
                alpha=0.5,
                checkpoint_interval=100,
                output_dir=str(result_dir / "checkpoints"),
                rerun=True,
                use_discrete_mode=True,
                cell_counts=cell_counts_df,
            )
            timings["gex_deconv_sec"] = time.time() - t0
            logger.info("GEX deconvolution completed in %.1fs", timings["gex_deconv_sec"])
        except Exception as e:
            logger.error("GEX deconvolution failed: %s", e)
            timings["gex_deconv_sec"] = -1.0
    else:
        timings["gex_deconv_sec"] = None

    # =========================================================================
    # Step 8: Save summary
    # =========================================================================
    logger.info("-" * 60)
    logger.info("STEP 8: Saving summary")
    logger.info("-" * 60)

    results = {
        "region_id": region_id,
        "sample_name": sample_name,
        "mode": "discrete_cellpose",
        "n_spots": int(gex_adata.n_obs),
        "n_genes": int(gex_adata.n_vars),
        "n_proteins": int(cite_adata.n_vars),
        "n_cell_types": len(ACHIEVABLE_7_CELL_PROFILE_DICT),
        "cell_types": list(ACHIEVABLE_7_CELL_PROFILE_DICT.keys()),
        "nuclei_stats": {
            "total_nuclei": int(nuclei_counts.sum()),
            "mean_per_spot": float(nuclei_counts.mean()),
            "median_per_spot": float(nuclei_counts.median()),
            "min_per_spot": int(nuclei_counts.min()),
            "max_per_spot": int(nuclei_counts.max()),
        },
        "cellpose_params": {
            "use_gpu": use_gpu,
            "diameter": cellpose_diameter,
            "spot_diameter_um": spot_diameter_um,
            "pixel_size_um": pixel_size_um,
        },
        "discrete_params": {
            "max_em_iterations": max_em_iterations,
        },
        "timings": timings,
        "output_dir": str(result_dir),
    }

    summary_path = result_dir / "benchmark_summary.json"
    with open(summary_path, "w") as f:
        json.dump(results, f, indent=2)
    logger.info("Saved summary: %s", summary_path)

    logger.info("=" * 70)
    logger.info("BENCHMARK COMPLETE")
    logger.info("  Total nuclei: %d", results["nuclei_stats"]["total_nuclei"])
    logger.info("  Cellpose time: %.1fs", timings["cellpose_sec"])
    logger.info("  Discrete assignment time: %.1fs", timings["discrete_assignment_sec"])
    if timings["gex_deconv_sec"] is not None and timings["gex_deconv_sec"] > 0:
        logger.info("  GEX deconvolution time: %.1fs", timings["gex_deconv_sec"])
    logger.info("=" * 70)

    return results


# =============================================================================
# CLI
# =============================================================================


def main():
    parser = argparse.ArgumentParser(
        description="Run discrete cell assignment benchmark with Cellpose nuclei segmentation"
    )

    # Region selection
    parser.add_argument(
        "--region",
        type=int,
        required=True,
        help=f"Xenium region ID to process (0-{N_REGIONS - 1})",
    )

    # I/O directories
    parser.add_argument(
        "--output-dir",
        type=str,
        default=str(REPO_ROOT / "Benchmarking/xenium_benchmarking/CITEgeist/output/discrete"),
        help="Output directory for results",
    )
    parser.add_argument(
        "--data-dir",
        type=str,
        default=str(DATA_DIR),
        help="Input directory with h5ad_objects/",
    )
    parser.add_argument(
        "--image-dir",
        type=str,
        default=str(IMAGE_DIR),
        help="Directory with morphology images",
    )

    # Cellpose parameters
    parser.add_argument(
        "--use-gpu",
        action="store_true",
        help="Use GPU for Cellpose segmentation",
    )
    parser.add_argument(
        "--cellpose-diameter",
        type=float,
        default=None,
        help="Cellpose nucleus diameter in pixels (auto-detect if not specified)",
    )
    parser.add_argument(
        "--spot-diameter-um",
        type=float,
        default=55.0,
        help="Spot diameter in microns (default: 55.0 for Visium geometry)",
    )

    # Discrete assignment parameters
    parser.add_argument(
        "--max-em-iterations",
        type=int,
        default=20,
        help="Maximum EM iterations for discrete assignment",
    )

    # GEX parameters
    parser.add_argument(
        "--skip-gex",
        action="store_true",
        help="Skip gene expression deconvolution",
    )
    parser.add_argument(
        "--min-counts",
        type=int,
        default=25,
        help="Minimum counts filter for GEX preprocessing",
    )
    parser.add_argument(
        "--scale-mode",
        choices=["per_marker", "global", "none"],
        default="per_marker",
        help="Antibody scaling mode: per_marker (default), global, or none",
    )

    args = parser.parse_args()

    # Validate region
    if args.region < 0 or args.region >= N_REGIONS:
        parser.error(f"Invalid region ID: {args.region}. Must be 0-{N_REGIONS - 1}")

    # Run benchmark
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    results = run_discrete_benchmark(
        region_id=args.region,
        output_dir=output_dir,
        data_dir=Path(args.data_dir),
        image_dir=Path(args.image_dir),
        use_gpu=args.use_gpu,
        cellpose_diameter=args.cellpose_diameter,
        max_em_iterations=args.max_em_iterations,
        run_gex=not args.skip_gex,
        min_counts=args.min_counts,
        spot_diameter_um=args.spot_diameter_um,
        scale_mode=args.scale_mode,
    )

    # Print summary
    print("\n" + "=" * 60)
    print(f"DISCRETE BENCHMARK COMPLETE: {results['sample_name']}")
    print("=" * 60)
    print(f"  Spots: {results['n_spots']}")
    print(f"  Nuclei detected: {results['nuclei_stats']['total_nuclei']}")
    print(f"  Cell types: {results['n_cell_types']}")
    print(f"  Timings:")
    print(f"    Cellpose: {results['timings']['cellpose_sec']:.1f}s")
    print(f"    Discrete assignment: {results['timings']['discrete_assignment_sec']:.1f}s")
    if results['timings']['gex_deconv_sec'] is not None:
        if results['timings']['gex_deconv_sec'] > 0:
            print(f"    GEX deconvolution: {results['timings']['gex_deconv_sec']:.1f}s")
        else:
            print(f"    GEX deconvolution: FAILED")
    print(f"  Output: {results['output_dir']}")


if __name__ == "__main__":
    main()
