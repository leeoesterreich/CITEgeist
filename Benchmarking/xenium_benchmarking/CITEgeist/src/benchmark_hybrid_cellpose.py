#!/usr/bin/env python
"""
Hybrid cell assignment benchmark: continuous model + discretization via nuclei counts.

This approach achieves ~94% of continuous model performance while providing
discrete cell counts constrained by Cellpose nuclei segmentation.

Pipeline:
1. Load morphology image for a Xenium pseudo-Visium region
2. Run Cellpose segmentation to get nuclei counts per spot
3. Run continuous CITEgeist model (CLR preprocessing, QP optimization)
4. Discretize continuous proportions using nuclei counts (largest remainder method)
5. Run GEX deconvolution using discrete cell counts
6. Save outputs for evaluation

Usage:
    python benchmark_hybrid_cellpose.py --region 0 --output-dir ./output/hybrid
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
# IMAGE AND COORDINATE LOADING (same as discrete benchmark)
# =============================================================================


def load_morphology_image(region_id: int, image_dir: Path) -> np.ndarray:
    """Load morphology PNG image for a region."""
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
    """Load coordinate info JSON for a region."""
    region_name = f"Xenium_region_{region_id}"
    json_path = image_dir / region_name / "coord_info.json"

    if not json_path.exists():
        raise FileNotFoundError(f"Missing coord_info.json: {json_path}")

    with open(json_path) as f:
        return json.load(f)


def convert_micron_to_pixel(coords_micron: np.ndarray, coord_info: Dict) -> np.ndarray:
    """Convert micron coordinates to pixel coordinates."""
    pixel_size = float(coord_info["pixel_size"])
    # Handle nested structure: micron_bounds.x_min / y_min
    offset_x = float(coord_info["micron_bounds"]["x_min"])
    offset_y = float(coord_info["micron_bounds"]["y_min"])

    coords_pixel = np.zeros_like(coords_micron)
    coords_pixel[:, 0] = (coords_micron[:, 0] - offset_x) / pixel_size
    coords_pixel[:, 1] = (coords_micron[:, 1] - offset_y) / pixel_size

    return coords_pixel


def run_cellpose_and_assign(
    image_rgb: np.ndarray,
    spot_coords_pixel: np.ndarray,
    spot_names: pd.Index,
    spot_diameter_um: float,
    pixel_size_um: float,
    use_gpu: bool = False,
    cellpose_diameter: Optional[float] = None,
) -> pd.Series:
    """Run Cellpose segmentation and assign nuclei to spots."""
    logger.info("Running Cellpose segmentation (use_gpu=%s, diameter=%s)",
                use_gpu, cellpose_diameter)

    t0 = time.time()
    masks, centroids_xy = run_cellpose_nuclei_segmentation(
        image_rgb_uint8=image_rgb,
        use_gpu=use_gpu,
        diameter=cellpose_diameter,
        model_type="nuclei",
    )
    seg_time = time.time() - t0
    n_nuclei = len(centroids_xy)

    logger.info(
        "Segmented %d nuclei in %.1fs (%.1f nuclei/sec)",
        n_nuclei, seg_time, n_nuclei / seg_time if seg_time > 0 else 0,
    )

    # Compute spot radius in pixels
    spot_radius_px = (spot_diameter_um / pixel_size_um) / 2.0

    # Assign nuclei centroids to spots
    nuclei_counts = assign_nuclei_centroids_to_spots(
        centroids_xy=centroids_xy,
        spot_centers_xy=spot_coords_pixel,
        spot_radius_px=spot_radius_px,
        spot_names=spot_names,
    )

    total_assigned = int(nuclei_counts.sum())
    logger.info(
        "Assigned %d/%d nuclei to spots (%.1f%%)",
        total_assigned, n_nuclei, 100.0 * total_assigned / n_nuclei if n_nuclei > 0 else 0,
    )

    return nuclei_counts


# =============================================================================
# MAIN BENCHMARK FUNCTION
# =============================================================================


def run_hybrid_benchmark(
    region_id: int,
    output_dir: Path,
    data_dir: Path = DATA_DIR,
    image_dir: Path = IMAGE_DIR,
    use_gpu: bool = False,
    cellpose_diameter: Optional[float] = None,
    run_gex: bool = True,
    min_counts: int = 25,
    spot_diameter_um: float = 55.0,
    lambda_laplacian: float = 0.1,
    # NEW parameters for module enrichment and KL regularization
    use_module_enrichment: bool = False,
    module_weight: float = 0.5,
    use_kl_regularization: bool = False,
    kl_temperature: float = 0.3,
    lambda_kl: float = 0.1,
) -> Dict[str, Any]:
    """
    Run hybrid CITEgeist benchmark for one region.

    Hybrid approach:
    1. Run continuous model to get optimal proportions
    2. Discretize using nuclei counts from Cellpose
    3. Run GEX deconvolution with discrete counts

    Args:
        region_id: Xenium region ID (0-4)
        output_dir: Output directory for results
        data_dir: Directory containing h5ad_objects/
        image_dir: Directory containing morphology images
        use_gpu: Whether to use GPU for Cellpose
        cellpose_diameter: Cellpose nucleus diameter (auto if None)
        run_gex: Whether to run GEX deconvolution
        min_counts: Minimum counts filter for GEX preprocessing
        spot_diameter_um: Spot diameter in microns
        lambda_laplacian: Spatial smoothing weight for continuous model
        use_module_enrichment: Use module-aware enrichment based on anchor genes
        module_weight: Weight for module signal vs base enrichment (0-1)
        use_kl_regularization: Use KL-divergence regularization instead of L2
        kl_temperature: Softmax temperature for KL target distribution
        lambda_kl: KL penalty weight

    Returns:
        Dictionary with benchmark results and timing
    """
    sample_name = f"Xenium_region_{region_id}"
    result_dir = output_dir / sample_name
    result_dir.mkdir(parents=True, exist_ok=True)

    logger.info("=" * 70)
    logger.info("HYBRID BENCHMARK: %s", sample_name)
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
    # Step 2: Load morphology image and run Cellpose
    # =========================================================================
    logger.info("-" * 60)
    logger.info("STEP 2: Cellpose nuclei segmentation")
    logger.info("-" * 60)

    t0 = time.time()
    image_rgb = load_morphology_image(region_id, image_dir)
    coord_info = load_coord_info(region_id, image_dir)

    spot_coords_micron = np.asarray(cite_adata.obsm["spatial"], dtype=np.float64)
    spot_coords_pixel = convert_micron_to_pixel(spot_coords_micron, coord_info)
    pixel_size_um = float(coord_info["pixel_size"])

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

    # =========================================================================
    # Step 3: Initialize CITEgeist model with continuous preprocessing
    # =========================================================================
    logger.info("-" * 60)
    logger.info("STEP 3: Initialize model (continuous mode)")
    logger.info("-" * 60)

    model = CitegeistModel(
        sample_name=sample_name,
        output_folder=str(result_dir),
        simulation=True,
        gene_expression_adata=gex_adata.copy(),
        antibody_capture_adata=cite_adata.copy(),
    )

    model.load_cell_profile_dict(ACHIEVABLE_7_CELL_PROFILE_DICT)
    logger.info("Loaded ACHIEVABLE_7 profile with %d cell types", len(ACHIEVABLE_7_CELL_PROFILE_DICT))

    # Use CLR preprocessing (continuous model style)
    t0 = time.time()
    model.filter_gex(min_counts=min_counts)
    model.preprocess_gex(target_sum=10000)
    model.preprocess_antibody()  # CLR normalization
    timings["preprocess_sec"] = time.time() - t0

    # =========================================================================
    # Step 4: Run continuous model
    # =========================================================================
    logger.info("-" * 60)
    logger.info("STEP 4: Continuous cell proportion model")
    logger.info("-" * 60)

    t0 = time.time()
    global_props_df, finetuned_props_df = model.run_cell_proportion_model(
        max_workers=8,
        lambda_laplacian=lambda_laplacian,
        skip_finetuning=False,
    )
    timings["continuous_sec"] = time.time() - t0
    logger.info("Continuous model completed in %.1fs", timings["continuous_sec"])

    # =========================================================================
    # Step 5: Discretize using nuclei counts
    # =========================================================================
    logger.info("-" * 60)
    logger.info("STEP 5: Discretize proportions with nuclei counts")
    logger.info("-" * 60)

    t0 = time.time()
    cell_counts_df = model.discretize_proportions(finetuned_props_df, nuclei_counts)
    timings["discretize_sec"] = time.time() - t0

    # Verify total matches
    total_discretized = cell_counts_df.values.sum()
    total_nuclei = nuclei_counts.sum()
    logger.info(
        "Discretization: %d cells assigned (nuclei count: %d, match: %s)",
        total_discretized, total_nuclei,
        "YES" if total_discretized == total_nuclei else "NO"
    )

    # Log cell type distribution
    logger.info("Cell type distribution:")
    for ct in cell_counts_df.columns:
        ct_total = cell_counts_df[ct].sum()
        pct = 100 * ct_total / total_discretized if total_discretized > 0 else 0
        logger.info("  %s: %d cells (%.1f%%)", ct, ct_total, pct)

    # Save predictions (as proportions for evaluation compatibility)
    prediction_path = result_dir / f"{sample_name}_deconv_predictions.csv"
    row_sums = cell_counts_df.values.sum(axis=1, keepdims=True)
    row_sums = np.maximum(row_sums, 1)
    proportions_df = pd.DataFrame(
        cell_counts_df.values / row_sums,
        index=cell_counts_df.index,
        columns=cell_counts_df.columns,
    )
    proportions_df.to_csv(prediction_path)
    logger.info("Saved proportions: %s", prediction_path)

    # Also save raw cell counts
    cell_counts_path = result_dir / f"{sample_name}_cell_counts.csv"
    cell_counts_df.to_csv(cell_counts_path)

    # =========================================================================
    # Step 6: GEX deconvolution (optional)
    # =========================================================================
    if run_gex:
        logger.info("-" * 60)
        logger.info("STEP 6: Gene expression deconvolution")
        logger.info("-" * 60)

        t0 = time.time()
        try:
            model.run_cell_expression_pass1(
                radius=None,
                alpha=0.5,
                checkpoint_interval=100,
                output_dir=str(result_dir / "checkpoints"),
                rerun=True,
                use_discrete_mode=True,
                cell_counts=cell_counts_df,
                # NEW parameters
                use_module_enrichment=use_module_enrichment,
                module_weight=module_weight,
                use_kl_regularization=use_kl_regularization,
                kl_temperature=kl_temperature,
                lambda_kl=lambda_kl,
            )
            timings["gex_deconv_sec"] = time.time() - t0
            logger.info("GEX deconvolution completed in %.1fs", timings["gex_deconv_sec"])
        except Exception as e:
            logger.error("GEX deconvolution failed: %s", e)
            timings["gex_deconv_sec"] = -1.0
    else:
        timings["gex_deconv_sec"] = None

    # =========================================================================
    # Step 7: Save summary
    # =========================================================================
    logger.info("-" * 60)
    logger.info("STEP 7: Saving summary")
    logger.info("-" * 60)

    results = {
        "region_id": region_id,
        "sample_name": sample_name,
        "mode": "hybrid_cellpose",
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
        "continuous_params": {
            "lambda_laplacian": lambda_laplacian,
        },
        "timings": timings,
        "output_dir": str(result_dir),
    }

    summary_path = result_dir / "benchmark_summary.json"
    with open(summary_path, "w") as f:
        json.dump(results, f, indent=2)
    logger.info("Saved summary: %s", summary_path)

    logger.info("=" * 70)
    logger.info("BENCHMARK COMPLETE: %s", sample_name)
    logger.info("=" * 70)

    return results


# =============================================================================
# CLI
# =============================================================================


def main():
    parser = argparse.ArgumentParser(description="Hybrid CITEgeist benchmark with Cellpose")
    parser.add_argument("--region", type=int, required=True,
                        help="Xenium region ID (0-4)")
    parser.add_argument("--output-dir", type=str, required=True,
                        help="Output directory for results")
    parser.add_argument("--data-dir", type=str, default=str(DATA_DIR),
                        help="Directory containing h5ad_objects/")
    parser.add_argument("--image-dir", type=str, default=str(IMAGE_DIR),
                        help="Directory containing morphology images")
    parser.add_argument("--use-gpu", action="store_true", default=False,
                        help="Use GPU for Cellpose")
    parser.add_argument("--cellpose-diameter", type=float, default=None,
                        help="Cellpose nucleus diameter (auto if not set)")
    parser.add_argument("--no-gex", action="store_true", default=False,
                        help="Skip GEX deconvolution")
    parser.add_argument("--min-counts", type=int, default=25,
                        help="Minimum counts for GEX filtering")
    parser.add_argument("--spot-diameter-um", type=float, default=55.0,
                        help="Spot diameter in microns")
    parser.add_argument("--lambda-laplacian", type=float, default=0.1,
                        help="Spatial smoothing weight")
    # NEW: Module enrichment and KL regularization
    parser.add_argument("--use-module-enrichment", action="store_true", default=False,
                        help="Use module-aware enrichment based on anchor genes")
    parser.add_argument("--module-weight", type=float, default=0.5,
                        help="Weight for module signal vs base enrichment (0-1)")
    parser.add_argument("--use-kl-regularization", action="store_true", default=False,
                        help="Use KL-divergence regularization instead of L2")
    parser.add_argument("--kl-temperature", type=float, default=0.3,
                        help="Softmax temperature for KL target")
    parser.add_argument("--lambda-kl", type=float, default=0.1,
                        help="KL penalty weight")
    args = parser.parse_args()

    results = run_hybrid_benchmark(
        region_id=args.region,
        output_dir=Path(args.output_dir),
        data_dir=Path(args.data_dir),
        image_dir=Path(args.image_dir),
        use_gpu=args.use_gpu,
        cellpose_diameter=args.cellpose_diameter,
        run_gex=not args.no_gex,
        min_counts=args.min_counts,
        spot_diameter_um=args.spot_diameter_um,
        lambda_laplacian=args.lambda_laplacian,
        # NEW parameters
        use_module_enrichment=args.use_module_enrichment,
        module_weight=args.module_weight,
        use_kl_regularization=args.use_kl_regularization,
        kl_temperature=args.kl_temperature,
        lambda_kl=args.lambda_kl,
    )

    print(f"\nBenchmark complete. Results saved to: {results['output_dir']}")


if __name__ == "__main__":
    main()
