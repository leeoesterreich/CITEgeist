#!/usr/bin/env python
"""
Test hybrid approach: run continuous model, then discretize using nuclei counts.

This should give us near-0.90 quality with integer outputs.
"""

import argparse
import json
import logging
import sys
import time
from pathlib import Path

import cv2
import numpy as np
import pandas as pd
import scanpy as sc

# Add paths
REPO_ROOT = Path(__file__).parent.parent.parent.parent.parent
sys.path.insert(0, str(REPO_ROOT))

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

# Simulation cell type profile
SIMULATION_CELL_PROFILE_DICT = {
    "B-cells": {"Major": ["B-cells_Protein_1", "B-cells_Protein_2"]},
    "CAFs": {"Major": ["CAFs_Protein_1", "CAFs_Protein_2"]},
    "Cancer Epithelial": {"Major": ["Cancer Epithelial_Protein_1", "Cancer Epithelial_Protein_2"]},
    "Endothelial": {"Major": ["Endothelial_Protein_1", "Endothelial_Protein_2"]},
    "Myeloid": {"Major": ["Myeloid_Protein_1", "Myeloid_Protein_2"]},
    "Normal Epithelial": {"Major": ["Normal Epithelial_Protein_1", "Normal Epithelial_Protein_2"]},
    "PVL": {"Major": ["PVL_Protein_1", "PVL_Protein_2"]},
    "Plasmablasts": {"Major": ["Plasmablasts_Protein_1", "Plasmablasts_Protein_2"]},
    "T-cells": {"Major": ["T-cells_Protein_1", "T-cells_Protein_2"]},
}


def discretize_proportions(proportions_df: pd.DataFrame, nuclei_counts: pd.Series) -> pd.DataFrame:
    """
    Convert continuous proportions to integer cell counts using nuclei counts.

    Uses largest remainder method to ensure sum(counts) = nuclei_count for each spot.
    """
    # Align indices
    common_spots = proportions_df.index.intersection(nuclei_counts.index)
    props = proportions_df.loc[common_spots]
    nuclei = nuclei_counts.loc[common_spots]

    cell_counts = pd.DataFrame(index=common_spots, columns=props.columns, dtype=float)

    for spot in common_spots:
        N = int(nuclei[spot])
        p = props.loc[spot].values

        if N == 0:
            cell_counts.loc[spot] = 0
            continue

        # Initial floor allocation
        raw_counts = p * N
        floor_counts = np.floor(raw_counts).astype(int)

        # Remainders for largest remainder method
        remainders = raw_counts - floor_counts

        # How many more cells to allocate
        deficit = N - floor_counts.sum()

        # Allocate to largest remainders
        if deficit > 0:
            # Get indices sorted by remainder (descending)
            indices = np.argsort(-remainders)
            for i in range(int(deficit)):
                floor_counts[indices[i]] += 1

        cell_counts.loc[spot] = floor_counts

    return cell_counts.astype(int)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--replicate-id", type=int, required=True)
    parser.add_argument("--condition", type=str, default="mixed")
    parser.add_argument("--mode", type=str, default="dapi")
    parser.add_argument("--output-dir", type=str, required=True)
    args = parser.parse_args()

    # Paths
    sccube_dir = Path("/ix1/alee/LO_LAB/Personal/Brent_Schlegel/bts76/Projects/CITEgeist/Wu_Visium/Simulations/scCube_12k/replicates")
    image_dir = REPO_ROOT / "Benchmarking/simulation_benchmarking/CITEgeist"
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    results = {"replicate_id": args.replicate_id, "condition": args.condition, "mode": args.mode}

    logger.info("=" * 60)
    logger.info("HYBRID TEST: replicate=%d, condition=%s", args.replicate_id, args.condition)
    logger.info("=" * 60)

    # Step 1: Load image and run Cellpose
    image_path = image_dir / args.condition / "images" / args.mode / f"Wu_rep_{args.replicate_id}_cellpose.png"
    logger.info("Loading image: %s", image_path)
    bgr = cv2.imread(str(image_path), cv2.IMREAD_COLOR)
    rgb = cv2.cvtColor(bgr, cv2.COLOR_BGR2RGB)

    logger.info("Running Cellpose...")
    start = time.time()
    masks, centroids = run_cellpose_nuclei_segmentation(rgb, use_gpu=False, model_type="nuclei")
    cellpose_time = time.time() - start
    logger.info("Cellpose completed in %.1fs, detected %d nuclei", cellpose_time, len(centroids))

    # Step 2: Assign nuclei to spots
    # Load spot coordinates
    h5ad_dir = sccube_dir / args.condition / "h5ad_objects"
    cite_path = h5ad_dir / f"Wu_rep_{args.replicate_id}_CITE.h5ad"
    gex_path = h5ad_dir / f"Wu_rep_{args.replicate_id}_GEX.h5ad"

    adata_cite = sc.read_h5ad(cite_path)
    adata_gex = sc.read_h5ad(gex_path)

    # Get spot coordinates from adata
    spot_coords = adata_cite.obs[["spot_x", "spot_y"]].copy()

    # Load cell index for coordinate transform
    cell_index_path = sccube_dir / args.condition / "ST_sim" / f"Wu_ST_{args.replicate_id}_index.csv"
    cell_index_df = pd.read_csv(cell_index_path)

    # Compute transform parameters
    padding = 100
    image_size = rgb.shape[0]
    cell_xy = cell_index_df[["point_x", "point_y"]].to_numpy(dtype=float)
    x_min, x_max = cell_xy[:, 0].min(), cell_xy[:, 0].max()
    y_min, y_max = cell_xy[:, 1].min(), cell_xy[:, 1].max()
    scale_x = (image_size - 2 * padding) / (x_max - x_min + 1e-6)
    scale_y = (image_size - 2 * padding) / (y_max - y_min + 1e-6)

    # Transform spot coordinates to image space
    spot_centers_px = np.column_stack([
        padding + (spot_coords["spot_x"].values - x_min) * scale_x,
        padding + (spot_coords["spot_y"].values - y_min) * scale_y,
    ])

    # Estimate spot radius
    spot_map = spot_coords[["spot_x", "spot_y"]]
    labeled = spot_map.loc[cell_index_df["spot"], ["spot_x", "spot_y"]].to_numpy(dtype=float)
    dx = (cell_xy[:, 0] - labeled[:, 0]) * scale_x
    dy = (cell_xy[:, 1] - labeled[:, 1]) * scale_y
    dist_px = np.sqrt(dx * dx + dy * dy)
    spot_radius_px = float(np.quantile(dist_px, 0.99))

    # Assign nuclei to spots
    pred_counts = assign_nuclei_centroids_to_spots(
        centroids_xy=centroids,
        spot_centers_xy=spot_centers_px,
        spot_radius_px=spot_radius_px,
        spot_names=spot_coords.index.tolist(),
    )
    logger.info("Assigned %d nuclei to %d spots", pred_counts.sum(), len(pred_counts))

    # Step 3: Run continuous model
    logger.info("Running continuous model...")
    model = CitegeistModel(
        sample_name=f"Wu_rep_{args.replicate_id}",
        output_folder=str(output_dir / f"Wu_rep_{args.replicate_id}"),
        simulation=True,
        gene_expression_adata=adata_gex,
        antibody_capture_adata=adata_cite,
    )

    model.preprocess_antibody()  # CLR normalization
    model.preprocess_gex(target_sum=10000)
    model.load_cell_profile_dict(SIMULATION_CELL_PROFILE_DICT)

    start = time.time()
    global_props_df, finetuned_props_df = model.run_cell_proportion_model(
        max_workers=8,
        lambda_laplacian=0.1,
        skip_finetuning=False,
    )
    continuous_time = time.time() - start
    logger.info("Continuous model completed in %.1fs", continuous_time)

    # Step 4: Discretize using nuclei counts
    logger.info("Discretizing continuous proportions using nuclei counts...")
    cell_counts = discretize_proportions(finetuned_props_df, pred_counts)

    # Convert back to proportions for evaluation
    pred_props = cell_counts.div(cell_counts.sum(axis=1), axis=0).fillna(0)

    # Step 5: Evaluate against ground truth
    gt_props_path = sccube_dir / args.condition / "ST_sim" / f"Wu_ST_{args.replicate_id}_prop.csv"
    gt_props = pd.read_csv(gt_props_path, index_col=0)
    gt_props = gt_props.drop(columns=["spot_x", "spot_y"], errors="ignore")

    # Align
    common_spots = pred_props.index.intersection(gt_props.index)
    common_types = pred_props.columns.intersection(gt_props.columns)

    pred_aligned = pred_props.loc[common_spots, common_types]
    gt_aligned = gt_props.loc[common_spots, common_types]

    # Also evaluate the raw continuous proportions
    continuous_aligned = finetuned_props_df.loc[common_spots, common_types]

    # Compute metrics
    hybrid_corr = np.corrcoef(pred_aligned.values.flatten(), gt_aligned.values.flatten())[0, 1]
    hybrid_rmse = np.sqrt(np.mean((pred_aligned.values - gt_aligned.values) ** 2))

    continuous_corr = np.corrcoef(continuous_aligned.values.flatten(), gt_aligned.values.flatten())[0, 1]
    continuous_rmse = np.sqrt(np.mean((continuous_aligned.values - gt_aligned.values) ** 2))

    results["hybrid_proportion_correlation"] = float(hybrid_corr)
    results["hybrid_proportion_rmse"] = float(hybrid_rmse)
    results["continuous_proportion_correlation"] = float(continuous_corr)
    results["continuous_proportion_rmse"] = float(continuous_rmse)
    results["discretization_loss"] = float(continuous_corr - hybrid_corr)

    logger.info("=" * 60)
    logger.info("RESULTS:")
    logger.info("  Continuous prop_corr: %.4f", continuous_corr)
    logger.info("  Hybrid (discretized) prop_corr: %.4f", hybrid_corr)
    logger.info("  Discretization loss: %.4f", continuous_corr - hybrid_corr)
    logger.info("=" * 60)

    # Save results
    results_path = output_dir / f"Wu_rep_{args.replicate_id}" / "benchmark_results.json"
    results_path.parent.mkdir(parents=True, exist_ok=True)
    with open(results_path, "w") as f:
        json.dump(results, f, indent=2)

    logger.info("Results saved to %s", results_path)


if __name__ == "__main__":
    main()
