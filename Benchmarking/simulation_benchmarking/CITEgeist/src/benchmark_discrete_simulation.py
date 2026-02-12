#!/usr/bin/env python
"""
Discrete cell assignment benchmark on scCube simulation data.

Pipeline:
1. Load synthetic Cellpose-compatible image
2. Run Cellpose segmentation with appropriate model (nuclei or cyto2)
3. Compare predicted vs ground truth nuclei counts
4. Run discrete cell assignment (IQP)
5. Evaluate against ground truth proportions
6. Run GEX deconvolution and evaluate
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

# Add paths
REPO_ROOT = Path(__file__).parent.parent.parent.parent.parent
sys.path.insert(0, str(REPO_ROOT))

BENCHMARK_ROOT = Path(__file__).parent.parent.parent
sys.path.insert(0, str(BENCHMARK_ROOT))

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

# Default paths
DEFAULT_SCCUBE_DIR = Path(
    "/ix1/alee/LO_LAB/Personal/Brent_Schlegel/bts76/Projects/CITEgeist/"
    "Wu_Visium/Simulations/scCube_12k/replicates"
)
DEFAULT_IMAGE_DIR = REPO_ROOT / "Benchmarking/simulation_benchmarking/CITEgeist"
DEFAULT_H5AD_DIR = DEFAULT_SCCUBE_DIR  # h5ad_objects are under condition folders

# Simulation cell type profile (9 cell types in simulation)
SIMULATION_CELL_PROFILE_DICT = {
    "B-cells": ["B-cells_Protein_1", "B-cells_Protein_2"],
    "CAFs": ["CAFs_Protein_1", "CAFs_Protein_2"],
    "Cancer Epithelial": ["Cancer Epithelial_Protein_1", "Cancer Epithelial_Protein_2"],
    "Endothelial": ["Endothelial_Protein_1", "Endothelial_Protein_2"],
    "Myeloid": ["Myeloid_Protein_1", "Myeloid_Protein_2"],
    "Normal Epithelial": ["Normal Epithelial_Protein_1", "Normal Epithelial_Protein_2"],
    "PVL": ["PVL_Protein_1", "PVL_Protein_2"],
    "Plasmablasts": ["Plasmablasts_Protein_1", "Plasmablasts_Protein_2"],
    "T-cells": ["T-cells_Protein_1", "T-cells_Protein_2"],
}

# Model type mapping
MODE_TO_MODEL = {
    "dapi": "nuclei",
    "h_and_e": "cyto2",
}


def load_synthetic_image(
    replicate_id: int,
    condition: str,
    mode: str,
    image_dir: Path,
) -> np.ndarray:
    """Load synthetic Cellpose-compatible image."""
    image_path = (
        image_dir / condition / "images" / mode / f"Wu_rep_{replicate_id}_cellpose.png"
    )
    if not image_path.exists():
        raise FileNotFoundError(f"Image not found: {image_path}")

    logger.info("Loading image: %s", image_path)
    bgr = cv2.imread(str(image_path), cv2.IMREAD_COLOR)
    if bgr is None:
        raise ValueError(f"Failed to load image: {image_path}")

    rgb = cv2.cvtColor(bgr, cv2.COLOR_BGR2RGB)
    logger.info("Image shape: %s", rgb.shape)
    return rgb


def load_ground_truth_counts(
    replicate_id: int,
    condition: str,
    image_dir: Path,
) -> pd.Series:
    """Load ground truth nuclei counts per spot."""
    counts_path = (
        image_dir / condition / "nuclei_counts" / f"Wu_rep_{replicate_id}_nuclei_counts.csv"
    )
    if not counts_path.exists():
        raise FileNotFoundError(f"Nuclei counts not found: {counts_path}")

    counts = pd.read_csv(counts_path, index_col=0).squeeze()
    logger.info("Loaded ground truth counts: %d spots, total %d cells", len(counts), counts.sum())
    return counts


def load_ground_truth_proportions(
    replicate_id: int,
    condition: str,
    sccube_dir: Path,
) -> pd.DataFrame:
    """Load ground truth cell type proportions per spot."""
    prop_path = sccube_dir / condition / "ST_sim" / f"Wu_ST_{replicate_id}_prop.csv"
    if not prop_path.exists():
        raise FileNotFoundError(f"Proportions not found: {prop_path}")

    df = pd.read_csv(prop_path, index_col=0)
    # Drop coordinate columns if present
    df = df.drop(columns=["spot_x", "spot_y"], errors="ignore")
    logger.info("Loaded ground truth proportions: %d spots, %d cell types", len(df), len(df.columns))
    return df


def run_cellpose_segmentation(
    image: np.ndarray,
    mode: str,
    use_gpu: bool = False,
    diameter: Optional[float] = None,
) -> tuple:
    """Run Cellpose segmentation with mode-appropriate model."""
    model_type = MODE_TO_MODEL[mode]
    logger.info("Running Cellpose with model_type=%s", model_type)

    start = time.time()
    masks, centroids = run_cellpose_nuclei_segmentation(
        image_rgb_uint8=image,
        use_gpu=use_gpu,
        diameter=diameter,
        model_type=model_type,
    )
    elapsed = time.time() - start

    n_detected = int(masks.max()) if masks.size > 0 else 0
    logger.info("Detected %d nuclei in %.1fs", n_detected, elapsed)

    return masks, centroids, elapsed


def compute_spot_coordinates(sccube_dir: Path, condition: str, replicate_id: int) -> pd.DataFrame:
    """Load spot coordinates from ST_sim meta file."""
    meta_path = sccube_dir / condition / "ST_sim" / f"Wu_ST_{replicate_id}_meta.csv"
    df = pd.read_csv(meta_path, index_col=0)
    return df[["spot", "spot_x", "spot_y"]].drop_duplicates().set_index("spot")


def run_benchmark(
    replicate_id: int,
    condition: str,
    mode: str,
    image_dir: Path,
    sccube_dir: Path,
    output_dir: Path,
    use_gpu: bool = False,
    cellpose_diameter: Optional[float] = None,
    max_em_iterations: int = 20,
) -> Dict[str, Any]:
    """Run full discrete benchmark pipeline."""
    logger.info("=" * 60)
    logger.info("BENCHMARK: replicate=%d, condition=%s, mode=%s", replicate_id, condition, mode)
    logger.info("=" * 60)

    results = {"replicate_id": replicate_id, "condition": condition, "mode": mode}
    timings = {}

    # Step 1: Load image
    image = load_synthetic_image(replicate_id, condition, mode, image_dir)

    # Step 2: Run Cellpose
    masks, centroids, cellpose_time = run_cellpose_segmentation(
        image, mode, use_gpu=use_gpu, diameter=cellpose_diameter
    )
    timings["cellpose_sec"] = cellpose_time

    # Step 3: Load spot coordinates and assign nuclei
    spot_coords_df = compute_spot_coordinates(sccube_dir, condition, replicate_id)
    spot_centers = spot_coords_df[["spot_x", "spot_y"]].values

    # Scale spot coordinates to image space (same scaling as image generation)
    # Image padding=100, image_size=8000, coord range 0-50
    padding = 100
    image_size = image.shape[0]
    coord_range = 50.0
    scale = (image_size - 2 * padding) / coord_range
    spot_centers_px = padding + spot_centers * scale

    # Estimate spot radius (based on typical spot spacing)
    spot_radius_px = scale * 0.8  # ~80% of 1 coordinate unit

    pred_counts = assign_nuclei_centroids_to_spots(
        centroids_xy=centroids,
        spot_centers_xy=spot_centers_px,
        spot_radius_px=spot_radius_px,
        spot_names=spot_coords_df.index.tolist(),
    )

    # Step 4: Compare to ground truth counts
    gt_counts = load_ground_truth_counts(replicate_id, condition, image_dir)

    # Align indices
    common_spots = pred_counts.index.intersection(gt_counts.index)
    pred_aligned = pred_counts.loc[common_spots]
    gt_aligned = gt_counts.loc[common_spots]

    count_corr = np.corrcoef(pred_aligned.values, gt_aligned.values)[0, 1]
    count_rmse = np.sqrt(np.mean((pred_aligned.values - gt_aligned.values) ** 2))

    results["nuclei_count_correlation"] = float(count_corr)
    results["nuclei_count_rmse"] = float(count_rmse)
    results["total_predicted_nuclei"] = int(pred_counts.sum())
    results["total_gt_nuclei"] = int(gt_counts.sum())

    logger.info("Nuclei count correlation: %.3f", count_corr)
    logger.info("Nuclei count RMSE: %.2f", count_rmse)

    # Step 5: Load h5ad files and run discrete assignment
    h5ad_dir = sccube_dir / condition / "h5ad_objects"
    cite_path = h5ad_dir / f"Wu_rep_{replicate_id}_CITE.h5ad"
    gex_path = h5ad_dir / f"Wu_rep_{replicate_id}_GEX.h5ad"

    if not cite_path.exists() or not gex_path.exists():
        logger.warning("h5ad files not found, skipping discrete assignment")
        results["discrete_assignment_skipped"] = True
        return results

    logger.info("Loading h5ad files...")
    adata_cite = sc.read_h5ad(cite_path)
    adata_gex = sc.read_h5ad(gex_path)

    # Initialize model
    model = CitegeistModel(
        sample_name=f"Wu_rep_{replicate_id}",
        output_folder=str(output_dir / condition / mode / f"Wu_rep_{replicate_id}"),
        simulation=True,
        gene_expression_adata=adata_gex,
        antibody_capture_adata=adata_cite,
    )

    model.preprocess_antibody_discrete()
    model.preprocess_gex(target_sum=10000)
    model.load_cell_profile_dict(SIMULATION_CELL_PROFILE_DICT)

    # Run discrete assignment
    start = time.time()
    cell_counts = model.run_discrete_cell_assignment(
        nuclei_counts=pred_aligned,
        max_iterations=max_em_iterations,
        convergence_threshold=1e-4,
        max_nuclei_cap=30,
    )
    timings["discrete_assignment_sec"] = time.time() - start

    # Compute predicted proportions
    pred_props = cell_counts.div(cell_counts.sum(axis=1), axis=0).fillna(0)

    # Load and compare to ground truth proportions
    gt_props = load_ground_truth_proportions(replicate_id, condition, sccube_dir)

    # Align to common spots and cell types
    common_spots = pred_props.index.intersection(gt_props.index)
    common_types = pred_props.columns.intersection(gt_props.columns)

    pred_props_aligned = pred_props.loc[common_spots, common_types]
    gt_props_aligned = gt_props.loc[common_spots, common_types]

    # Compute metrics
    prop_corr = np.corrcoef(
        pred_props_aligned.values.flatten(), gt_props_aligned.values.flatten()
    )[0, 1]
    prop_rmse = np.sqrt(
        np.mean((pred_props_aligned.values - gt_props_aligned.values) ** 2)
    )

    results["proportion_correlation"] = float(prop_corr)
    results["proportion_rmse"] = float(prop_rmse)
    results["timings"] = timings

    logger.info("Proportion correlation: %.3f", prop_corr)
    logger.info("Proportion RMSE: %.4f", prop_rmse)

    # Save results
    result_dir = output_dir / condition / mode / f"Wu_rep_{replicate_id}"
    result_dir.mkdir(parents=True, exist_ok=True)

    with open(result_dir / "benchmark_results.json", "w") as f:
        json.dump(results, f, indent=2)

    pred_props.to_csv(result_dir / "predicted_proportions.csv")
    cell_counts.to_csv(result_dir / "cell_counts.csv")

    logger.info("Results saved to %s", result_dir)
    return results


def main():
    parser = argparse.ArgumentParser(
        description="Run discrete cell assignment benchmark on simulation data"
    )
    parser.add_argument("--replicate-id", type=int, required=True, help="Replicate ID (0-4)")
    parser.add_argument("--condition", choices=["high_seg", "mixed"], required=True)
    parser.add_argument("--mode", choices=["dapi", "h_and_e"], required=True)
    parser.add_argument("--image-dir", type=Path, default=DEFAULT_IMAGE_DIR)
    parser.add_argument("--sccube-dir", type=Path, default=DEFAULT_SCCUBE_DIR)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--use-gpu", action="store_true")
    parser.add_argument("--cellpose-diameter", type=float, default=None)
    parser.add_argument("--max-em-iterations", type=int, default=20)
    args = parser.parse_args()

    results = run_benchmark(
        replicate_id=args.replicate_id,
        condition=args.condition,
        mode=args.mode,
        image_dir=args.image_dir,
        sccube_dir=args.sccube_dir,
        output_dir=args.output_dir,
        use_gpu=args.use_gpu,
        cellpose_diameter=args.cellpose_diameter,
        max_em_iterations=args.max_em_iterations,
    )

    print("\n=== RESULTS ===")
    print(f"Nuclei count correlation: {results.get('nuclei_count_correlation', 'N/A'):.3f}")
    print(f"Proportion correlation: {results.get('proportion_correlation', 'N/A'):.3f}")


if __name__ == "__main__":
    main()