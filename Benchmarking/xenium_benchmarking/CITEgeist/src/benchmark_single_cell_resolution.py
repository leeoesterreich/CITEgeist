#!/usr/bin/env python
"""
Single-cell resolution benchmark using Module 3b nucleus assignment pipeline.

This script tests the new single-cell resolution features on Xenium pseudo-Visium data:
1. Loads existing CITEgeist hybrid results (proportions + deconvolved GEX)
2. Runs Cellpose on the morphology image to get nuclei masks (not just counts)
3. Maps nuclei centroids to spots
4. Runs the new run_nucleus_assignment() pipeline
5. Distributes GEX to individual cells
6. Creates single-cell AnnData output

This is a proof-of-concept test focused on verifying the pipeline runs end-to-end.

Usage:
    python benchmark_single_cell_resolution.py --region 0 --output-dir ./output/single_cell
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

from CITEgeist.model.segmentation import (
    run_cellpose_nuclei_segmentation,
)
from CITEgeist.model.module3b_nucleus_assignment import run_nucleus_assignment
from CITEgeist.model.cell_level_gex import distribute_gex_to_cells
from CITEgeist.model.single_cell_output import create_single_cell_adata

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
HYBRID_OUTPUT_DIR = REPO_ROOT / "Benchmarking/xenium_benchmarking/CITEgeist/output/hybrid_cellpose"
N_REGIONS = 5  # Xenium_region_0 through Xenium_region_4


# =============================================================================
# IMAGE AND COORDINATE LOADING
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
    offset_x = float(coord_info["micron_bounds"]["x_min"])
    offset_y = float(coord_info["micron_bounds"]["y_min"])

    coords_pixel = np.zeros_like(coords_micron)
    coords_pixel[:, 0] = (coords_micron[:, 0] - offset_x) / pixel_size
    coords_pixel[:, 1] = (coords_micron[:, 1] - offset_y) / pixel_size

    return coords_pixel


def map_nuclei_to_spots(
    centroids_xy: np.ndarray,
    spot_coords_pixel: np.ndarray,
    spot_names: pd.Index,
    spot_radius_px: float,
) -> pd.DataFrame:
    """
    Map each nucleus centroid to the nearest spot within spot radius.

    Returns DataFrame with columns: nucleus_id, spot_id
    """
    from scipy.spatial import cKDTree

    if centroids_xy.size == 0:
        return pd.DataFrame(columns=['nucleus_id', 'spot_id'])

    tree = cKDTree(spot_coords_pixel)
    dists, idxs = tree.query(centroids_xy, k=1, workers=-1)

    # Build mapping for nuclei within spot radius
    rows = []
    for i, (dist, spot_idx) in enumerate(zip(dists, idxs)):
        if dist <= spot_radius_px:
            rows.append({
                'nucleus_id': i + 1,  # 1-indexed to match mask labels
                'spot_id': spot_names[spot_idx],
            })

    return pd.DataFrame(rows)


# =============================================================================
# LOAD EXISTING RESULTS
# =============================================================================


def load_hybrid_results(region_id: int, hybrid_dir: Path) -> Dict[str, Any]:
    """Load existing hybrid benchmark results for a region."""
    sample_name = f"Xenium_region_{region_id}"
    result_dir = hybrid_dir / sample_name

    if not result_dir.exists():
        raise FileNotFoundError(f"Missing hybrid results: {result_dir}")

    # Load proportions (finetuned)
    prop_path = result_dir / f"{sample_name}_cell_prop_finetuned_results.csv"
    if not prop_path.exists():
        prop_path = result_dir / f"{sample_name}_deconv_predictions.csv"
    proportions_df = pd.read_csv(prop_path, index_col=0)

    # Load cell counts
    counts_path = result_dir / f"{sample_name}_cell_counts.csv"
    cell_counts_df = pd.read_csv(counts_path, index_col=0) if counts_path.exists() else None

    # Load nuclei counts
    nuclei_path = result_dir / f"{sample_name}_nuclei_counts.csv"
    nuclei_counts = pd.read_csv(nuclei_path, index_col=0).squeeze("columns") if nuclei_path.exists() else None

    # Load deconvolved GEX
    gex_path = result_dir / f"{sample_name}_gene_expression_pass1.parquet"
    if gex_path.exists():
        deconvolved_gex = pd.read_parquet(gex_path)
    else:
        deconvolved_gex = None

    # Load summary
    summary_path = result_dir / "benchmark_summary.json"
    with open(summary_path) as f:
        summary = json.load(f)

    return {
        'sample_name': sample_name,
        'proportions': proportions_df,
        'cell_counts': cell_counts_df,
        'nuclei_counts': nuclei_counts,
        'deconvolved_gex': deconvolved_gex,
        'summary': summary,
        'result_dir': result_dir,
    }


# =============================================================================
# MAIN BENCHMARK FUNCTION
# =============================================================================


def run_single_cell_benchmark(
    region_id: int,
    output_dir: Path,
    data_dir: Path = DATA_DIR,
    image_dir: Path = IMAGE_DIR,
    hybrid_dir: Path = HYBRID_OUTPUT_DIR,
    use_gpu: bool = False,
    cellpose_diameter: Optional[float] = None,
    spot_diameter_um: float = 55.0,
) -> Dict[str, Any]:
    """
    Run single-cell resolution pipeline for one region.

    Args:
        region_id: Xenium region ID (0-4)
        output_dir: Output directory for results
        data_dir: Directory containing h5ad_objects/
        image_dir: Directory containing morphology images
        hybrid_dir: Directory containing hybrid benchmark results
        use_gpu: Whether to use GPU for Cellpose
        cellpose_diameter: Cellpose nucleus diameter (auto if None)
        spot_diameter_um: Spot diameter in microns

    Returns:
        Dictionary with benchmark results and statistics
    """
    sample_name = f"Xenium_region_{region_id}"
    result_dir = output_dir / sample_name
    result_dir.mkdir(parents=True, exist_ok=True)

    logger.info("=" * 70)
    logger.info("SINGLE-CELL RESOLUTION BENCHMARK: %s", sample_name)
    logger.info("=" * 70)

    timings = {}

    # =========================================================================
    # Step 1: Load existing hybrid results
    # =========================================================================
    logger.info("-" * 60)
    logger.info("STEP 1: Loading existing hybrid results")
    logger.info("-" * 60)

    t0 = time.time()
    hybrid_results = load_hybrid_results(region_id, hybrid_dir)
    timings["load_hybrid_sec"] = time.time() - t0

    proportions_df = hybrid_results['proportions']
    deconvolved_gex = hybrid_results['deconvolved_gex']
    nuclei_counts = hybrid_results['nuclei_counts']

    logger.info("Loaded proportions: %d spots x %d cell types",
                proportions_df.shape[0], proportions_df.shape[1])

    if deconvolved_gex is not None:
        logger.info("Loaded deconvolved GEX: %d rows x %d genes",
                    deconvolved_gex.shape[0], deconvolved_gex.shape[1])
    else:
        logger.warning("No deconvolved GEX found - will skip GEX distribution")

    if nuclei_counts is not None:
        logger.info("Nuclei counts: total=%d, mean=%.1f/spot",
                    nuclei_counts.sum(), nuclei_counts.mean())

    # =========================================================================
    # Step 2: Load morphology image and run Cellpose (get mask, not just counts)
    # =========================================================================
    logger.info("-" * 60)
    logger.info("STEP 2: Cellpose segmentation (full mask)")
    logger.info("-" * 60)

    t0 = time.time()
    image_rgb = load_morphology_image(region_id, image_dir)
    coord_info = load_coord_info(region_id, image_dir)
    pixel_size_um = float(coord_info["pixel_size"])

    logger.info("Running Cellpose (use_gpu=%s, diameter=%s)", use_gpu, cellpose_diameter)
    masks, centroids_xy = run_cellpose_nuclei_segmentation(
        image_rgb_uint8=image_rgb,
        use_gpu=use_gpu,
        diameter=cellpose_diameter,
        model_type="nuclei",
    )
    timings["cellpose_sec"] = time.time() - t0

    n_nuclei = int(masks.max())
    logger.info("Segmented %d nuclei in %.1fs", n_nuclei, timings["cellpose_sec"])
    logger.info("Mask shape: %s, centroids shape: %s", masks.shape, centroids_xy.shape)

    # Save mask for debugging
    mask_path = result_dir / f"{sample_name}_cellpose_mask.npy"
    np.save(mask_path, masks)
    logger.info("Saved mask: %s", mask_path)

    # =========================================================================
    # Step 3: Map nuclei centroids to spots
    # =========================================================================
    logger.info("-" * 60)
    logger.info("STEP 3: Map nuclei to spots")
    logger.info("-" * 60)

    t0 = time.time()

    # Load spatial coordinates from original CITE data
    cite_path = data_dir / "h5ad_objects" / f"{sample_name}_CITE.h5ad"
    cite_adata = sc.read_h5ad(cite_path)

    spot_coords_micron = np.asarray(cite_adata.obsm["spatial"], dtype=np.float64)
    spot_coords_pixel = convert_micron_to_pixel(spot_coords_micron, coord_info)
    spot_radius_px = (spot_diameter_um / pixel_size_um) / 2.0

    nuclei_spot_map = map_nuclei_to_spots(
        centroids_xy=centroids_xy,
        spot_coords_pixel=spot_coords_pixel,
        spot_names=cite_adata.obs_names,
        spot_radius_px=spot_radius_px,
    )
    timings["map_nuclei_sec"] = time.time() - t0

    n_mapped = len(nuclei_spot_map)
    pct_mapped = 100.0 * n_mapped / n_nuclei if n_nuclei > 0 else 0
    logger.info("Mapped %d/%d nuclei to spots (%.1f%%)", n_mapped, n_nuclei, pct_mapped)

    # Save nuclei-spot mapping
    mapping_path = result_dir / f"{sample_name}_nuclei_spot_map.csv"
    nuclei_spot_map.to_csv(mapping_path, index=False)

    # =========================================================================
    # Step 4: Run Module 3b nucleus assignment
    # =========================================================================
    logger.info("-" * 60)
    logger.info("STEP 4: Module 3b nucleus assignment")
    logger.info("-" * 60)

    t0 = time.time()

    # Prepare proportions with spot_id column for pipeline
    proportions_with_spot = proportions_df.copy()
    proportions_with_spot['spot_id'] = proportions_with_spot.index

    # Get cell types from profile dict
    cell_types = list(ACHIEVABLE_7_CELL_PROFILE_DICT.keys())

    # Filter to only cell types present in proportions
    cell_types = [ct for ct in cell_types if ct in proportions_df.columns]

    # Compute nuclei counts per spot from mapping (not using pre-computed for accuracy)
    mapped_nuclei_counts = nuclei_spot_map.groupby('spot_id').size()
    # Reindex to match all spots (fill missing with 0)
    mapped_nuclei_counts = mapped_nuclei_counts.reindex(proportions_df.index, fill_value=0)

    logger.info("Cell types for assignment: %s", cell_types)
    logger.info("Spots with nuclei: %d/%d", (mapped_nuclei_counts > 0).sum(), len(mapped_nuclei_counts))

    try:
        assignment_result = run_nucleus_assignment(
            mask=masks,
            nuclei_spot_map=nuclei_spot_map,
            proportions=proportions_with_spot,
            nuclei_counts=mapped_nuclei_counts,
            cell_types=cell_types,
        )
        timings["assignment_sec"] = time.time() - t0

        n_assigned = len(assignment_result.assignments)
        logger.info("Assigned %d nuclei to cell types in %.1fs", n_assigned, timings["assignment_sec"])

        # Log cell type distribution
        ct_counts = pd.Series(assignment_result.assignments.values()).value_counts()
        logger.info("Cell type distribution:")
        for ct, count in ct_counts.items():
            pct = 100 * count / n_assigned if n_assigned > 0 else 0
            logger.info("  %s: %d (%.1f%%)", ct, count, pct)

    except Exception as e:
        logger.error("Nucleus assignment failed: %s", e)
        import traceback
        traceback.print_exc()
        assignment_result = None
        timings["assignment_sec"] = -1.0

    # =========================================================================
    # Step 5: Distribute GEX to individual cells
    # =========================================================================
    if assignment_result is not None and deconvolved_gex is not None:
        logger.info("-" * 60)
        logger.info("STEP 5: Distribute GEX to cells")
        logger.info("-" * 60)

        t0 = time.time()
        try:
            cell_gex = distribute_gex_to_cells(
                deconvolved_gex=deconvolved_gex,
                assignments=assignment_result.assignments,
                nucleus_spot_map=nuclei_spot_map,
            )
            timings["distribute_gex_sec"] = time.time() - t0

            n_cells_with_gex = (cell_gex.sum(axis=1) > 0).sum()
            logger.info("Distributed GEX to %d/%d cells in %.1fs",
                        n_cells_with_gex, len(cell_gex), timings["distribute_gex_sec"])

            # Basic GEX statistics
            total_counts = cell_gex.sum().sum()
            mean_per_cell = cell_gex.sum(axis=1).mean()
            logger.info("Total GEX counts: %.0f, mean per cell: %.1f", total_counts, mean_per_cell)

        except Exception as e:
            logger.error("GEX distribution failed: %s", e)
            import traceback
            traceback.print_exc()
            cell_gex = None
            timings["distribute_gex_sec"] = -1.0
    else:
        cell_gex = None
        timings["distribute_gex_sec"] = None

    # =========================================================================
    # Step 6: Create single-cell AnnData
    # =========================================================================
    if assignment_result is not None and cell_gex is not None:
        logger.info("-" * 60)
        logger.info("STEP 6: Create single-cell AnnData")
        logger.info("-" * 60)

        t0 = time.time()
        try:
            sc_adata = create_single_cell_adata(
                cell_gex=cell_gex,
                morphology_features=assignment_result.morphology_features,
                assignments=assignment_result.assignments,
                sample_name=sample_name,
                classifier=assignment_result.classifier,
            )
            timings["create_adata_sec"] = time.time() - t0

            logger.info("Created single-cell AnnData: %d cells x %d genes",
                        sc_adata.n_obs, sc_adata.n_vars)
            logger.info("Cell types: %s", sc_adata.obs['cell_type'].value_counts().to_dict())

            # Save AnnData
            adata_path = result_dir / f"{sample_name}_single_cell.h5ad"
            sc_adata.write_h5ad(adata_path)
            logger.info("Saved single-cell AnnData: %s", adata_path)

        except Exception as e:
            logger.error("AnnData creation failed: %s", e)
            import traceback
            traceback.print_exc()
            sc_adata = None
            timings["create_adata_sec"] = -1.0
    else:
        sc_adata = None
        timings["create_adata_sec"] = None

    # =========================================================================
    # Step 7: Save summary
    # =========================================================================
    logger.info("-" * 60)
    logger.info("STEP 7: Saving summary")
    logger.info("-" * 60)

    # Compute statistics
    stats = {
        "n_nuclei_segmented": n_nuclei,
        "n_nuclei_mapped_to_spots": n_mapped,
        "pct_nuclei_mapped": pct_mapped,
        "n_spots_total": len(proportions_df),
        "n_spots_with_nuclei": int((mapped_nuclei_counts > 0).sum()),
    }

    if assignment_result is not None:
        stats["n_nuclei_assigned"] = len(assignment_result.assignments)
        ct_counts = pd.Series(assignment_result.assignments.values()).value_counts()
        stats["cell_type_distribution"] = ct_counts.to_dict()

    if sc_adata is not None:
        stats["n_cells_in_adata"] = sc_adata.n_obs
        stats["n_genes_in_adata"] = sc_adata.n_vars
        stats["mean_gex_per_cell"] = float(sc_adata.X.sum(axis=1).mean())

    results = {
        "region_id": region_id,
        "sample_name": sample_name,
        "pipeline": "single_cell_resolution_module3b",
        "cell_types": cell_types,
        "cellpose_params": {
            "use_gpu": use_gpu,
            "diameter": cellpose_diameter,
            "spot_diameter_um": spot_diameter_um,
            "pixel_size_um": pixel_size_um,
        },
        "statistics": stats,
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
    parser = argparse.ArgumentParser(
        description="Single-cell resolution benchmark using Module 3b"
    )
    parser.add_argument("--region", type=int, required=True,
                        help="Xenium region ID (0-4)")
    parser.add_argument("--output-dir", type=str, required=True,
                        help="Output directory for results")
    parser.add_argument("--data-dir", type=str, default=str(DATA_DIR),
                        help="Directory containing h5ad_objects/")
    parser.add_argument("--image-dir", type=str, default=str(IMAGE_DIR),
                        help="Directory containing morphology images")
    parser.add_argument("--hybrid-dir", type=str, default=str(HYBRID_OUTPUT_DIR),
                        help="Directory containing hybrid benchmark results")
    parser.add_argument("--use-gpu", action="store_true", default=False,
                        help="Use GPU for Cellpose")
    parser.add_argument("--cellpose-diameter", type=float, default=None,
                        help="Cellpose nucleus diameter (auto if not set)")
    parser.add_argument("--spot-diameter-um", type=float, default=55.0,
                        help="Spot diameter in microns")
    args = parser.parse_args()

    results = run_single_cell_benchmark(
        region_id=args.region,
        output_dir=Path(args.output_dir),
        data_dir=Path(args.data_dir),
        image_dir=Path(args.image_dir),
        hybrid_dir=Path(args.hybrid_dir),
        use_gpu=args.use_gpu,
        cellpose_diameter=args.cellpose_diameter,
        spot_diameter_um=args.spot_diameter_um,
    )

    print(f"\nBenchmark complete. Results saved to: {results['output_dir']}")


if __name__ == "__main__":
    main()
