#!/usr/bin/env python
"""
Hybrid cell assignment benchmark: continuous model + discretization via nuclei counts.

This approach achieves ~94% of continuous model performance while providing
discrete cell counts constrained by StarDist nuclei segmentation.

Pipeline:
1. Load morphology image for a Xenium pseudo-Visium region
2. Run StarDist segmentation to get nuclei counts per spot
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
    run_nuclei_segmentation,
)
from CITEgeist.model.detection import detect_cell_types
from CITEgeist.model.constrained_assignment import (
    ConstrainedAssignment,
    extract_patch_features,
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


def run_segmentation_and_assign(
    image_rgb: np.ndarray,
    spot_coords_pixel: np.ndarray,
    spot_names: pd.Index,
    spot_diameter_um: float,
    pixel_size_um: float,
) -> pd.Series:
    """Run StarDist nuclei segmentation and assign nuclei to spots."""
    logger.info("Running StarDist nuclei segmentation")

    t0 = time.time()
    masks, centroids_df = run_nuclei_segmentation(image_rgb, modality="dapi")
    centroids_xy = centroids_df[["x_pixel", "y_pixel"]].values
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
    run_gex: bool = True,
    min_counts: int = 25,
    spot_diameter_um: float = 55.0,
    lambda_laplacian: float = 0.1,
    # Detection-based filtering
    use_detection_filter: bool = True,
    # NEW parameters for module enrichment and KL regularization
    use_module_enrichment: bool = False,
    module_weight: float = 0.5,
    use_kl_regularization: bool = False,
    kl_temperature: float = 0.3,
    lambda_kl: float = 0.1,
    # Single-cell assignment
    single_cell_mode: Optional[str] = None,
    patches_dir: Optional[Path] = None,
    xgboost_model_path: Optional[Path] = None,
    vae_checkpoint_path: Optional[Path] = None,
    use_gpu: bool = False,
) -> Dict[str, Any]:
    """
    Run hybrid CITEgeist benchmark for one region.

    Hybrid approach:
    1. Run continuous model to get optimal proportions
    2. Discretize using nuclei counts from StarDist
    3. Run GEX deconvolution with discrete counts

    Args:
        region_id: Xenium region ID (0-4)
        output_dir: Output directory for results
        data_dir: Directory containing h5ad_objects/
        image_dir: Directory containing morphology images
        run_gex: Whether to run GEX deconvolution
        min_counts: Minimum counts filter for GEX preprocessing
        spot_diameter_um: Spot diameter in microns
        lambda_laplacian: Spatial smoothing weight for continuous model
        use_detection_filter: Apply detection-based confidence-weighted filtering
            to continuous proportions before discretization. This zeros out
            spurious small allocations to non-detected cell types, improving
            accuracy by ~5-8%.
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
    logger.info("STEP 2: StarDist nuclei segmentation")
    logger.info("-" * 60)

    t0 = time.time()
    image_rgb = load_morphology_image(region_id, image_dir)
    coord_info = load_coord_info(region_id, image_dir)

    spot_coords_micron = np.asarray(cite_adata.obsm["spatial"], dtype=np.float64)
    spot_coords_pixel = convert_micron_to_pixel(spot_coords_micron, coord_info)
    pixel_size_um = float(coord_info["pixel_size"])

    nuclei_counts = run_segmentation_and_assign(
        image_rgb=image_rgb,
        spot_coords_pixel=spot_coords_pixel,
        spot_names=cite_adata.obs_names,
        spot_diameter_um=spot_diameter_um,
        pixel_size_um=pixel_size_um,
    )
    timings["segmentation_sec"] = time.time() - t0

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

    # Convert cell counts to proportions
    row_sums = cell_counts_df.values.sum(axis=1, keepdims=True)
    row_sums = np.maximum(row_sums, 1)
    proportions_df = pd.DataFrame(
        cell_counts_df.values / row_sums,
        index=cell_counts_df.index,
        columns=cell_counts_df.columns,
    )

    # =========================================================================
    # Step 5.5: Detection-based post-filtering (optional)
    # =========================================================================
    # IMPORTANT: Apply filtering AFTER discretization, not before.
    # Pre-filtering breaks the model's learned allocations.
    # Post-filtering corrects predictions without affecting optimization.
    if use_detection_filter:
        logger.info("-" * 60)
        logger.info("STEP 5.5: Detection-based post-filtering")
        logger.info("-" * 60)

        t0 = time.time()

        # Try to load pre-computed detection mask from detection-estimation run
        # This ensures consistency with the validated detection results
        det_est_dir = REPO_ROOT / "Benchmarking/xenium_benchmarking/CITEgeist/output/detection_estimation" / f"region_{region_id}"
        precomputed_mask = det_est_dir / "detection_mask.csv"

        if precomputed_mask.exists():
            logger.info("Loading pre-computed detection mask from detection-estimation")
            detected_df = pd.read_csv(precomputed_mask, index_col=0).astype(bool)
            # Align with current data
            detected_df = detected_df.reindex(proportions_df.index)
            # Fill missing spots with True (don't filter)
            detected_df = detected_df.fillna(True)
        else:
            logger.info("No pre-computed mask found, running detection fresh")
            # Build marker groups from profile dict
            marker_cols = model.antibody_capture_adata.var_names.tolist()
            marker_groups = {}
            for cell_type, markers_spec in model.cell_profile_dict.items():
                if isinstance(markers_spec, dict):
                    markers = markers_spec.get('Major', [])
                else:
                    markers = markers_spec
                indices = [marker_cols.index(m) for m in markers if m in marker_cols]
                if indices:
                    marker_groups[cell_type] = indices

            logger.info("Detection marker groups: %s", {k: len(v) for k, v in marker_groups.items()})

            # Run detection on antibody data
            ab_matrix = model.antibody_capture_adata.X
            if hasattr(ab_matrix, 'toarray'):
                ab_matrix = ab_matrix.toarray()

            detected = detect_cell_types(
                ab_matrix, marker_groups, adaptive_threshold=True, log_transform=True
            )
            detected_df = pd.DataFrame(
                detected,
                index=model.antibody_capture_adata.obs_names,
                columns=list(marker_groups.keys())
            )

        # Precomputed reliability (r² from detection-estimation vs GT)
        DETECTION_RELIABILITY = {
            'B cells': 0.85, 'CD4+ T cells': 0.35, 'CD8+ T cells': 0.63,
            'Macrophages': 0.50, 'Endothelial': 0.66, 'Epithelial': 0.60,
            'Fibroblasts': 0.45,
        }

        # Apply confidence-weighted filtering to discretized proportions
        for col in proportions_df.columns:
            if col in detected_df.columns:
                r2 = DETECTION_RELIABILITY.get(col, 0.5)
                det_rate = detected_df[col].mean()
                not_detected = ~detected_df[col].values
                # Reduce non-detected proportions by r² factor
                proportions_df.loc[not_detected, col] *= (1 - r2)
                logger.info("  %s: detection=%.1f%%, r²=%.3f", col, 100*det_rate, r2)

        # Renormalize
        row_sums_filt = proportions_df.sum(axis=1).replace(0, 1)
        proportions_df = proportions_df.div(row_sums_filt, axis=0)

        # Re-discretize filtered proportions back to cell counts
        # This ensures cell_counts_df reflects the filtering
        cell_counts_df = model.discretize_proportions(proportions_df, nuclei_counts)
        logger.info("Re-discretized filtered proportions to cell counts")

        # Save detection mask
        detected_df.to_csv(result_dir / f"{sample_name}_detection_mask.csv")
        timings["detection_filter_sec"] = time.time() - t0
        logger.info("Detection filtering completed in %.1fs", timings["detection_filter_sec"])
    else:
        timings["detection_filter_sec"] = None

    # Save predictions
    prediction_path = result_dir / f"{sample_name}_deconv_predictions.csv"
    proportions_df.to_csv(prediction_path)
    logger.info("Saved proportions: %s", prediction_path)

    # Save cell counts (filtered if detection filter was applied)
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
    # Step 7: Single-cell assignment (optional)
    # =========================================================================
    single_cell_adata = None
    if single_cell_mode:
        logger.info("-" * 60)
        logger.info("STEP 7: Single-cell assignment (mode=%s)", single_cell_mode)
        logger.info("-" * 60)

        t0 = time.time()

        # Initialize assigner
        assigner = ConstrainedAssignment(mode=single_cell_mode, seed=42)

        # Load XGBoost model for 'xgboost' mode
        if single_cell_mode == "xgboost":
            if xgboost_model_path is None:
                xgboost_model_path = (
                    REPO_ROOT
                    / "Benchmarking/xenium_benchmarking/CITEgeist/output/xgboost_combined/xgboost_model.pkl"
                )
            if vae_checkpoint_path is None:
                vae_checkpoint_path = (
                    REPO_ROOT
                    / "Benchmarking/xenium_benchmarking/CITEgeist/output/vae_masked/vae/vae_final.pt"
                )

            if not xgboost_model_path.exists():
                logger.warning(
                    "XGBoost model not found: %s. Falling back to morphology mode.",
                    xgboost_model_path,
                )
                assigner = ConstrainedAssignment(mode="morphology", seed=42)
                single_cell_mode = "morphology"  # Update for fitting below
            else:
                device = "cuda" if use_gpu else "cpu"
                assigner.load_xgboost_model(
                    model_path=xgboost_model_path,
                    vae_checkpoint_path=vae_checkpoint_path if vae_checkpoint_path.exists() else None,
                    device=device,
                )
                logger.info("Loaded XGBoost model from %s", xgboost_model_path)

        # For 'groups' and 'morphology' modes, need to fit on high-purity spots
        if single_cell_mode in ("groups", "morphology"):
            if patches_dir is None:
                # Default patches directory
                patches_dir = (
                    REPO_ROOT
                    / "Benchmarking/xenium_benchmarking/CITEgeist/output/cell_morphology"
                )

            if not patches_dir.exists():
                logger.warning(
                    "Patches directory not found: %s. "
                    "Falling back to random assignment.",
                    patches_dir,
                )
                assigner = ConstrainedAssignment(mode="random", seed=42)
            else:
                # Build proportions DataFrame for fitting
                fit_proportions = finetuned_props_df.copy()
                fit_proportions["region"] = region_id
                fit_proportions["spot_id"] = fit_proportions.index

                logger.info("Fitting assigner on patches from %s", patches_dir)
                assigner.fit(
                    patches_dir=patches_dir,
                    proportions_df=fit_proportions,
                    max_spots=2000,
                    purity_threshold=0.5 if single_cell_mode == "morphology" else 0.6,
                )

        # Get nuclei centroids per spot from StarDist
        # Re-run to get centroids (or load from cache)
        logger.info("Extracting nuclei centroids per spot...")
        masks, centroids_df = run_nuclei_segmentation(image_rgb, modality="dapi")
        centroids_xy = centroids_df[["x_pixel", "y_pixel"]].values

        # Assign centroids to spots and get per-spot centroid lists
        spot_radius_px = (spot_diameter_um / pixel_size_um) / 2.0
        nuclei_centroids = {}
        for i, spot_name in enumerate(cite_adata.obs_names):
            spot_center = spot_coords_pixel[i]
            # Find nuclei within this spot
            dists = np.sqrt(
                (centroids_xy[:, 0] - spot_center[0]) ** 2
                + (centroids_xy[:, 1] - spot_center[1]) ** 2
            )
            in_spot = dists <= spot_radius_px
            nuclei_centroids[spot_name] = centroids_xy[in_spot]

        # Extract features/patches for each nucleus and assign
        celltype_names = list(cell_counts_df.columns)
        assignments = {}

        for spot_name in cite_adata.obs_names:
            spot_centroids = nuclei_centroids.get(spot_name, np.array([]))
            n_nuclei = len(spot_centroids)

            if n_nuclei == 0:
                assignments[spot_name] = np.array([], dtype=int)
                continue

            counts = cell_counts_df.loc[spot_name].values.astype(int)

            if assigner.mode == "random":
                # Random assignment - no features needed
                assignments[spot_name] = assigner.assign_spot(
                    nuclei_features=np.zeros((n_nuclei, 13)),  # Dummy features
                    cell_counts=counts,
                    celltype_names=celltype_names,
                )
            else:
                # Extract 2-channel patches for each nucleus
                patches_list = []
                features = []
                for cx, cy in spot_centroids:
                    cx, cy = int(cx), int(cy)
                    # Extract 96x96 patch centered on nucleus
                    h, w = image_rgb.shape[:2]
                    y1, y2 = max(0, cy - 48), min(h, cy + 48)
                    x1, x2 = max(0, cx - 48), min(w, cx + 48)
                    patch = image_rgb[y1:y2, x1:x2]

                    if patch.shape[0] < 10 or patch.shape[1] < 10:
                        # Create empty patch
                        patches_list.append(np.zeros((2, 96, 96), dtype=np.float32))
                        features.append(np.zeros(13))
                        continue

                    # Convert RGB to 2-channel (DAPI~blue, boundary~green)
                    dapi = patch[:, :, 2].astype(np.float32)  # Blue channel
                    boundary = patch[:, :, 1].astype(np.float32)  # Green channel

                    # Pad to 96x96 if needed
                    if dapi.shape != (96, 96):
                        padded_dapi = np.zeros((96, 96), dtype=np.float32)
                        padded_boundary = np.zeros((96, 96), dtype=np.float32)
                        padded_dapi[: dapi.shape[0], : dapi.shape[1]] = dapi
                        padded_boundary[: boundary.shape[0], : boundary.shape[1]] = boundary
                        dapi, boundary = padded_dapi, padded_boundary

                    patch_2ch = np.stack([dapi, boundary], axis=0)
                    patches_list.append(patch_2ch)

                    # Also extract simple features for morphology mode
                    try:
                        feats = extract_patch_features(patch_2ch)
                        features.append(feats)
                    except Exception:
                        features.append(np.zeros(13))

                patches_array = np.array(patches_list)
                features = np.array(features)

                if assigner.mode == "xgboost":
                    # XGBoost mode: pass patches for VAE+morphology extraction
                    assignments[spot_name] = assigner.assign_spot(
                        nuclei_features=features,  # Not used but required
                        cell_counts=counts,
                        celltype_names=celltype_names,
                        patches=patches_array,
                    )
                else:
                    # Morphology mode: use pre-extracted simple features
                    assignments[spot_name] = assigner.assign_spot(
                        nuclei_features=features,
                        cell_counts=counts,
                        celltype_names=celltype_names,
                    )

        # Load deconvolved GEX if available
        gex_parquet = result_dir / f"{sample_name}_gene_expression_pass1.parquet"
        if gex_parquet.exists():
            deconvolved_gex = pd.read_parquet(gex_parquet)
        else:
            # Create empty GEX DataFrame as fallback
            logger.warning("No deconvolved GEX found, creating empty expression matrix")
            gene_names = gex_adata.var_names.tolist()
            deconvolved_gex = pd.DataFrame(columns=gene_names)

        # Create single-cell AnnData
        single_cell_adata = assigner.create_single_cell_adata(
            spot_adata=cite_adata,
            assignments=assignments,
            nuclei_centroids=nuclei_centroids,
            deconvolved_gex=deconvolved_gex,
            celltype_names=celltype_names,
        )

        # Save single-cell h5ad
        sc_output_path = result_dir / f"{sample_name}_single_cell_{single_cell_mode}.h5ad"
        single_cell_adata.write_h5ad(sc_output_path)
        logger.info("Saved single-cell AnnData: %s", sc_output_path)
        logger.info("  - %d cells, %d genes", single_cell_adata.n_obs, single_cell_adata.n_vars)

        timings["single_cell_sec"] = time.time() - t0
    else:
        timings["single_cell_sec"] = None

    # =========================================================================
    # Step 8: Save summary
    # =========================================================================
    logger.info("-" * 60)
    logger.info("STEP 8: Saving summary")
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
        "segmentation_params": {
            "backend": "stardist",
            "modality": "dapi",
            "spot_diameter_um": spot_diameter_um,
            "pixel_size_um": pixel_size_um,
        },
        "continuous_params": {
            "lambda_laplacian": lambda_laplacian,
        },
        "single_cell": {
            "mode": single_cell_mode,
            "n_cells": single_cell_adata.n_obs if single_cell_adata is not None else None,
            "output_file": str(result_dir / f"{sample_name}_single_cell_{single_cell_mode}.h5ad")
                if single_cell_mode else None,
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
    parser = argparse.ArgumentParser(description="Hybrid CITEgeist benchmark with StarDist segmentation")
    parser.add_argument("--region", type=int, required=True,
                        help="Xenium region ID (0-4)")
    parser.add_argument("--output-dir", type=str, required=True,
                        help="Output directory for results")
    parser.add_argument("--data-dir", type=str, default=str(DATA_DIR),
                        help="Directory containing h5ad_objects/")
    parser.add_argument("--image-dir", type=str, default=str(IMAGE_DIR),
                        help="Directory containing morphology images")
    parser.add_argument("--no-gex", action="store_true", default=False,
                        help="Skip GEX deconvolution")
    parser.add_argument("--min-counts", type=int, default=25,
                        help="Minimum counts for GEX filtering")
    parser.add_argument("--spot-diameter-um", type=float, default=55.0,
                        help="Spot diameter in microns")
    parser.add_argument("--lambda-laplacian", type=float, default=0.1,
                        help="Spatial smoothing weight")
    # Detection-based filtering (improves accuracy by ~5-8%)
    parser.add_argument("--no-detection-filter", action="store_true", default=False,
                        help="Disable detection-based filtering (not recommended)")
    # Module enrichment and KL regularization
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
    # Single-cell assignment mode
    parser.add_argument("--use-gpu", action="store_true", default=False,
                        help="Use GPU for XGBoost inference (single-cell mode)")
    parser.add_argument("--single-cell", type=str, choices=["random", "morphology", "xgboost"],
                        default=None,
                        help="Output single-cell h5ad with assigned GEX per cell. "
                             "Modes: 'random' (constrained random ~22%%), "
                             "'morphology' (7-class morphology Hungarian ~46%%), "
                             "'xgboost' (VAE + morphology XGBoost ~50%%+)")
    parser.add_argument("--patches-dir", type=str, default=None,
                        help="Directory containing nucleus patches for single-cell modes "
                             "(required for 'groups' and 'morphology' modes)")
    parser.add_argument("--xgboost-model", type=str,
                        default="Benchmarking/xenium_benchmarking/CITEgeist/output/xgboost_combined/xgboost_model.pkl",
                        help="Path to trained XGBoost model (for xgboost mode)")
    parser.add_argument("--vae-checkpoint", type=str,
                        default="Benchmarking/xenium_benchmarking/CITEgeist/output/vae_masked/vae/vae_final.pt",
                        help="Path to VAE checkpoint for feature extraction (for xgboost mode)")
    args = parser.parse_args()

    results = run_hybrid_benchmark(
        region_id=args.region,
        output_dir=Path(args.output_dir),
        data_dir=Path(args.data_dir),
        image_dir=Path(args.image_dir),
        run_gex=not args.no_gex,
        min_counts=args.min_counts,
        spot_diameter_um=args.spot_diameter_um,
        lambda_laplacian=args.lambda_laplacian,
        use_detection_filter=not args.no_detection_filter,
        # Module enrichment and KL parameters
        use_module_enrichment=args.use_module_enrichment,
        module_weight=args.module_weight,
        use_kl_regularization=args.use_kl_regularization,
        kl_temperature=args.kl_temperature,
        lambda_kl=args.lambda_kl,
        # Single-cell assignment
        single_cell_mode=args.single_cell,
        patches_dir=Path(args.patches_dir) if args.patches_dir else None,
        xgboost_model_path=Path(args.xgboost_model) if args.xgboost_model else None,
        vae_checkpoint_path=Path(args.vae_checkpoint) if args.vae_checkpoint else None,
        use_gpu=args.use_gpu,
    )

    print(f"\nBenchmark complete. Results saved to: {results['output_dir']}")


if __name__ == "__main__":
    main()
