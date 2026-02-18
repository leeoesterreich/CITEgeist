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
# Format: {"cell_type": {"Major": [markers]}}
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


def load_cell_index(sccube_dir: Path, condition: str, replicate_id: int) -> pd.DataFrame:
    """Load cell-level simulation table with cell coordinates and spot labels."""
    index_path = sccube_dir / condition / "ST_sim" / f"Wu_ST_{replicate_id}_index.csv"
    return pd.read_csv(index_path, index_col=0)


def compute_image_space_spot_centers(
    spot_coords_df: pd.DataFrame,
    cell_index_df: pd.DataFrame,
    image_size: int,
    padding: int = 100,
) -> np.ndarray:
    """
    Map spot centers into image pixels using the same transform as image generation.

    The synthetic images are generated by scaling point_x/point_y using per-replicate
    min/max values from the cell coordinates, not a fixed [0, 50] range.
    """
    cell_xy = cell_index_df[["point_x", "point_y"]].to_numpy(dtype=float)
    x_min, x_max = cell_xy[:, 0].min(), cell_xy[:, 0].max()
    y_min, y_max = cell_xy[:, 1].min(), cell_xy[:, 1].max()

    scale_x = (image_size - 2 * padding) / (x_max - x_min + 1e-6)
    scale_y = (image_size - 2 * padding) / (y_max - y_min + 1e-6)

    spot_centers = spot_coords_df[["spot_x", "spot_y"]].to_numpy(dtype=float)
    return np.column_stack(
        [
            padding + (spot_centers[:, 0] - x_min) * scale_x,
            padding + (spot_centers[:, 1] - y_min) * scale_y,
        ]
    )


def estimate_spot_radius_px(
    spot_coords_df: pd.DataFrame,
    cell_index_df: pd.DataFrame,
    image_size: int,
    padding: int = 100,
    quantile: float = 0.99,
) -> float:
    """Estimate spot radius in pixels from empirical cell-to-labeled-spot distances."""
    cell_xy = cell_index_df[["point_x", "point_y"]].to_numpy(dtype=float)
    x_min, x_max = cell_xy[:, 0].min(), cell_xy[:, 0].max()
    y_min, y_max = cell_xy[:, 1].min(), cell_xy[:, 1].max()
    scale_x = (image_size - 2 * padding) / (x_max - x_min + 1e-6)
    scale_y = (image_size - 2 * padding) / (y_max - y_min + 1e-6)

    spot_map = spot_coords_df[["spot_x", "spot_y"]]
    labeled = spot_map.loc[cell_index_df["spot"], ["spot_x", "spot_y"]].to_numpy(dtype=float)
    dx = (cell_xy[:, 0] - labeled[:, 0]) * scale_x
    dy = (cell_xy[:, 1] - labeled[:, 1]) * scale_y
    dist_px = np.sqrt(dx * dx + dy * dy)

    return float(np.quantile(dist_px, quantile))


def compute_geometry_ground_truth_counts(
    spot_coords_df: pd.DataFrame,
    cell_index_df: pd.DataFrame,
    image_size: int,
    spot_radius_px: float,
    padding: int = 100,
) -> pd.Series:
    """Compute geometry-consistent GT nuclei counts by assigning true cell coordinates to spots."""
    cell_xy = cell_index_df[["point_x", "point_y"]].to_numpy(dtype=float)
    x_min, x_max = cell_xy[:, 0].min(), cell_xy[:, 0].max()
    y_min, y_max = cell_xy[:, 1].min(), cell_xy[:, 1].max()
    scale_x = (image_size - 2 * padding) / (x_max - x_min + 1e-6)
    scale_y = (image_size - 2 * padding) / (y_max - y_min + 1e-6)

    cell_centroids_px = np.column_stack(
        [
            padding + (cell_xy[:, 0] - x_min) * scale_x,
            padding + (cell_xy[:, 1] - y_min) * scale_y,
        ]
    )
    spot_centers_px = compute_image_space_spot_centers(
        spot_coords_df=spot_coords_df,
        cell_index_df=cell_index_df,
        image_size=image_size,
        padding=padding,
    )
    return assign_nuclei_centroids_to_spots(
        centroids_xy=cell_centroids_px,
        spot_centers_xy=spot_centers_px,
        spot_radius_px=spot_radius_px,
        spot_names=spot_coords_df.index.tolist(),
    )


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
    lambda_sparse: float = 0.0,
    scale_mode: str = "per_marker",
    lambda_prior: float = 0.0,
    global_solve: bool = True,
    global_time_limit: float = 300.0,
    global_mip_gap: float = 0.05,
    use_continuous_prior: bool = False,
    lambda_continuous: float = 50.0,
) -> Dict[str, Any]:
    """Run full discrete benchmark pipeline."""
    logger.info("=" * 60)
    logger.info("BENCHMARK: replicate=%d, condition=%s, mode=%s", replicate_id, condition, mode)
    logger.info("=" * 60)

    results = {"replicate_id": replicate_id, "condition": condition, "mode": mode, "lambda_sparse": lambda_sparse, "scale_mode": scale_mode, "lambda_prior": lambda_prior, "global_solve": global_solve, "global_time_limit": global_time_limit, "global_mip_gap": global_mip_gap, "use_continuous_prior": use_continuous_prior, "lambda_continuous": lambda_continuous}
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
    cell_index_df = load_cell_index(sccube_dir, condition, replicate_id)

    # Map spots to image pixels using the same transform as generate_cellpose_images.py
    padding = 100
    image_size = image.shape[0]
    spot_centers_px = compute_image_space_spot_centers(
        spot_coords_df=spot_coords_df,
        cell_index_df=cell_index_df,
        image_size=image_size,
        padding=padding,
    )
    # Use empirical radius from simulation geometry, not hard-coded 1.2 units.
    spot_radius_px = estimate_spot_radius_px(
        spot_coords_df=spot_coords_df,
        cell_index_df=cell_index_df,
        image_size=image_size,
        padding=padding,
        quantile=0.99,
    )

    pred_counts = assign_nuclei_centroids_to_spots(
        centroids_xy=centroids,
        spot_centers_xy=spot_centers_px,
        spot_radius_px=spot_radius_px,
        spot_names=spot_coords_df.index.tolist(),
    )

    # Step 4: Compare to geometry-consistent ground truth counts.
    gt_counts_geometry = compute_geometry_ground_truth_counts(
        spot_coords_df=spot_coords_df,
        cell_index_df=cell_index_df,
        image_size=image_size,
        spot_radius_px=spot_radius_px,
        padding=padding,
    )

    # Align indices
    common_spots = pred_counts.index.intersection(gt_counts_geometry.index)
    pred_aligned = pred_counts.loc[common_spots]
    gt_aligned = gt_counts_geometry.loc[common_spots]

    count_corr = np.corrcoef(pred_aligned.values, gt_aligned.values)[0, 1]
    count_rmse = np.sqrt(np.mean((pred_aligned.values - gt_aligned.values) ** 2))

    results["nuclei_count_correlation"] = float(count_corr)
    results["nuclei_count_rmse"] = float(count_rmse)
    results["total_predicted_nuclei"] = int(pred_counts.sum())
    results["total_gt_nuclei"] = int(gt_counts_geometry.sum())
    results["spot_radius_px"] = float(spot_radius_px)

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

    model.preprocess_antibody_discrete(scale_mode=scale_mode)
    model.preprocess_gex(target_sum=10000)
    model.load_cell_profile_dict(SIMULATION_CELL_PROFILE_DICT)

    # Estimate prior proportions from raw marker data if lambda_prior > 0
    prior_proportions = None
    if lambda_prior > 0:
        from CITEgeist.model.gurobi_impl import (
            estimate_prior_proportions_from_markers,
            map_antibodies_to_profiles_v2,
        )
        # Get raw marker data (before preprocessing) and assignment matrix
        raw_data, marker_names, assignment_matrix, cell_type_names = map_antibodies_to_profiles_v2(
            adata_cite, SIMULATION_CELL_PROFILE_DICT
        )
        prior_proportions = estimate_prior_proportions_from_markers(
            raw_marker_data=raw_data,
            assignment_matrix=assignment_matrix,
            cell_type_names=cell_type_names,
        )
        logger.info("Prior proportions: %s", dict(zip(cell_type_names, prior_proportions)))

    # Optional: Run continuous optimization first to get prior for discrete
    continuous_prior = None
    if use_continuous_prior:
        logger.info("Running continuous optimization first (hybrid mode)...")
        start_cont = time.time()

        # Run continuous proportion model (returns tuple: global, finetuned)
        global_props_df, finetuned_props_df = model.run_cell_proportion_model(
            max_workers=8,
            lambda_laplacian=0.1,
            skip_finetuning=False,
        )
        timings["continuous_sec"] = time.time() - start_cont
        logger.info("Continuous optimization completed in %.1fs", timings["continuous_sec"])

        # Use finetuned proportions as the continuous prior
        continuous_prior = finetuned_props_df.values
        logger.info("Continuous prior shape: %s", continuous_prior.shape)
        results["use_continuous_prior"] = True
        results["lambda_continuous"] = lambda_continuous
    else:
        results["use_continuous_prior"] = False

    # Run discrete assignment
    start = time.time()
    cell_counts = model.run_discrete_cell_assignment(
        nuclei_counts=pred_aligned,
        max_em_iterations=max_em_iterations,
        beta_convergence_tol=1e-4,
        max_nuclei_cap=30,
        lambda_sparse=lambda_sparse,
        prior_proportions=prior_proportions,
        lambda_prior=lambda_prior,
        global_solve=global_solve,
        global_time_limit=global_time_limit,
        global_mip_gap=global_mip_gap,
        continuous_prior=continuous_prior,
        lambda_continuous=lambda_continuous if use_continuous_prior else 0.0,
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

    logger.info("Proportion correlation: %.3f", prop_corr)
    logger.info("Proportion RMSE: %.4f", prop_rmse)

    # Step 6: Run GEX deconvolution and evaluate against ground truth
    logger.info("Running GEX deconvolution...")
    start = time.time()
    try:
        deconv_result = model.run_cell_expression_pass1(
            use_discrete_mode=True,
            cell_counts=cell_counts,
        )
        timings["gex_deconv_sec"] = time.time() - start
        logger.info("GEX deconvolution completed in %.1fs", timings["gex_deconv_sec"])

        # Extract spotwise profiles from result dict
        # Format: {spot_idx: np.ndarray of shape (n_cell_types, n_genes)}
        spotwise_profiles = deconv_result["spotwise_profiles"]
        n_spots, n_types, n_genes = deconv_result["dimensions"]

        # Get cell type names and gene names
        cell_type_names = list(SIMULATION_CELL_PROFILE_DICT.keys())
        gene_names = model.gene_expression_adata.var_names.tolist()
        spot_names = model.gene_expression_adata.obs_names.tolist()

        # Export predicted layers in format expected by benchmarking_gex.py
        # {cell_type}_layer.csv: spots as rows, genes as columns (same as GT format)
        result_dir = output_dir / condition / mode / f"Wu_rep_{replicate_id}"
        layers_dir = result_dir / "layers"
        layers_dir.mkdir(parents=True, exist_ok=True)

        for ct_idx, ct_name in enumerate(cell_type_names):
            # Extract this cell type's expression across all spots
            # Shape: (n_spots, n_genes) - spots as rows, genes as columns
            ct_expr = np.zeros((n_spots, n_genes))
            for spot_idx, profile in spotwise_profiles.items():
                ct_expr[spot_idx, :] = profile[ct_idx, :]

            # Create DataFrame with spots as rows, genes as columns (matches GT format)
            ct_df = pd.DataFrame(ct_expr, index=spot_names, columns=gene_names)
            layer_path = layers_dir / f"{ct_name}_layer.csv"
            ct_df.to_csv(layer_path)

        logger.info("Exported %d cell type layers to %s", len(cell_type_names), layers_dir)

        # Load ground truth GEX and compute metrics using existing framework
        gt_gex_dir = sccube_dir / condition / "ST_GEX_sim" / f"sample_{replicate_id}" / "layers"

        if gt_gex_dir.exists():
            # Import benchmarking function
            sys.path.insert(0, str(REPO_ROOT / "Benchmarking/simulation_benchmarking/src"))
            from benchmarking_gex import calculate_rmse

            logger.info("Comparing against ground truth GEX: %s", gt_gex_dir)
            gex_metrics = calculate_rmse(str(gt_gex_dir), str(layers_dir), normalize="range")

            results["gex_average_rmse"] = float(gex_metrics["average_rmse"])
            results["gex_median_rmse"] = float(gex_metrics["median_rmse"])
            results["gex_average_nrmse"] = float(gex_metrics["average_nrmse"])
            results["gex_median_nrmse"] = float(gex_metrics["median_nrmse"])
            results["gex_average_mae"] = float(gex_metrics["average_mae"])
            results["gex_median_mae"] = float(gex_metrics["median_mae"])
            results["gex_per_celltype"] = {
                ct: {k: float(v) for k, v in metrics.items()}
                for ct, metrics in gex_metrics["metrics_per_cell_type"].items()
            }

            logger.info("GEX Average RMSE: %.4f", gex_metrics["average_rmse"])
            logger.info("GEX Average NRMSE: %.4f", gex_metrics["average_nrmse"])
        else:
            logger.warning("Ground truth GEX not found at %s", gt_gex_dir)

    except Exception as e:
        logger.error("GEX deconvolution failed: %s", str(e))
        import traceback
        logger.error(traceback.format_exc())
        results["gex_error"] = str(e)

    results["timings"] = timings

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
    parser.add_argument("--lambda-sparse", type=float, default=0.0,
                        help="Sparsity penalty for IQP (default: 0.0)")
    parser.add_argument("--scale-mode", choices=["per_marker", "global", "none"],
                        default="per_marker",
                        help="Antibody scaling mode: per_marker (default), global, or none")
    parser.add_argument("--lambda-prior", type=float, default=0.0,
                        help="Prior regularization weight (default: 0.0). Uses marker-based prior.")
    parser.add_argument("--global-solve", action="store_true", default=True,
                        help="Use global IQP solver (default: True)")
    parser.add_argument("--no-global-solve", action="store_false", dest="global_solve",
                        help="Use per-spot IQP solver instead of global")
    parser.add_argument("--global-time-limit", type=float, default=300.0,
                        help="Time limit for global IQP solver in seconds (default: 300)")
    parser.add_argument("--global-mip-gap", type=float, default=0.05,
                        help="MIP gap tolerance for global solver (default: 0.05)")
    parser.add_argument("--use-continuous-prior", action="store_true", default=False,
                        help="Run continuous optimization first, use as prior for discrete (hybrid mode)")
    parser.add_argument("--lambda-continuous", type=float, default=50.0,
                        help="Regularization weight for continuous prior (default: 50.0)")
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
        lambda_sparse=args.lambda_sparse,
        scale_mode=args.scale_mode,
        lambda_prior=args.lambda_prior,
        global_solve=args.global_solve,
        global_time_limit=args.global_time_limit,
        global_mip_gap=args.global_mip_gap,
        use_continuous_prior=args.use_continuous_prior,
        lambda_continuous=args.lambda_continuous,
    )

    print("\n=== RESULTS ===")
    print(f"Nuclei count correlation: {results.get('nuclei_count_correlation', 'N/A'):.3f}")
    print(f"Proportion correlation: {results.get('proportion_correlation', 'N/A'):.3f}")


if __name__ == "__main__":
    main()
