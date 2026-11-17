#!/usr/bin/env python
"""
Unified simulation benchmark: StarDist + continuous proportions + SC evaluation.

Pipeline:
1. Load/generate synthetic nuclei image
2. Run StarDist segmentation (or GT-ablation)
3. Run CITEgeist Module 3 continuous proportions
4. Discretize proportions to integer counts
5. Evaluate: spot-level proportions/GEX + single-cell accuracy
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
from scipy.spatial import cKDTree
from scipy.stats import pearsonr

# ---------------------------------------------------------------------------
# Repo / import paths
# ---------------------------------------------------------------------------
REPO_ROOT = Path(__file__).resolve().parent.parent.parent.parent.parent
sys.path.insert(0, str(REPO_ROOT))

# ---------------------------------------------------------------------------
# Default data directories
# ---------------------------------------------------------------------------
DEFAULT_SCCUBE_DIR = REPO_ROOT / "Benchmarking" / "simulation_benchmarking" / "replicates"
DEFAULT_IMAGE_DIR = REPO_ROOT / "Benchmarking" / "simulation_benchmarking" / "CITEgeist"

# ---------------------------------------------------------------------------
# Simulation cell type profile (9 cell types, 2 markers each)
# ---------------------------------------------------------------------------
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

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)


def discretize_proportions(proportions: np.ndarray, total: int) -> np.ndarray:
    """Largest-remainder method: convert proportions to integer counts summing to total."""
    prop_sum = proportions.sum()
    if prop_sum > 0:
        proportions = proportions / prop_sum
    scaled = proportions * total
    floors = np.floor(scaled).astype(int)
    remainders = scaled - floors
    deficit = total - floors.sum()
    if deficit > 0:
        indices = np.argsort(-remainders)
        for i in range(int(deficit)):
            floors[indices[i]] += 1
    elif deficit < 0:
        indices = np.argsort(remainders)
        for i in range(int(-deficit)):
            if floors[indices[i]] > 0:
                floors[indices[i]] -= 1
    return floors


def discretize_proportions_df(props_df: pd.DataFrame, nuclei_counts: pd.Series) -> pd.DataFrame:
    """Discretize proportions for all spots."""
    common = props_df.index.intersection(nuclei_counts.index)
    result = pd.DataFrame(0, index=common, columns=props_df.columns)
    for spot in common:
        result.loc[spot] = discretize_proportions(props_df.loc[spot].values, int(nuclei_counts[spot]))
    return result


def compute_hungarian_sc_accuracy(gt_labels: list, pred_counts: dict) -> float:
    """Best-case SC accuracy for one spot: sum of min(gt_count, pred_count) per type."""
    n = len(gt_labels)
    if n == 0:
        return 1.0
    gt_counts = {}
    for label in gt_labels:
        gt_counts[label] = gt_counts.get(label, 0) + 1
    types = set(gt_counts.keys()) | set(pred_counts.keys())
    correct = sum(min(gt_counts.get(t, 0), pred_counts.get(t, 0)) for t in types)
    return correct / n


def compute_hungarian_sc_accuracy_allspots(gt_index: pd.DataFrame, pred_counts_df: pd.DataFrame) -> float:
    """Aggregate best-case SC accuracy across all spots."""
    total_correct = 0
    total_cells = 0
    for spot in pred_counts_df.index:
        spot_gt = gt_index[gt_index["spot"] == spot]["Cell_type"].tolist()
        if len(spot_gt) == 0:
            continue
        spot_pred = pred_counts_df.loc[spot].to_dict()
        acc = compute_hungarian_sc_accuracy(spot_gt, spot_pred)
        total_correct += acc * len(spot_gt)
        total_cells += len(spot_gt)
    return total_correct / total_cells if total_cells > 0 else 0.0


def evaluate_segmentation_nn(
    gt_coords: np.ndarray,
    detected_coords: np.ndarray,
    threshold_px: float,
) -> dict:
    """Evaluate segmentation quality via nearest-neighbor matching."""
    if len(detected_coords) == 0:
        return {
            "detection_rate": 0.0,
            "false_positive_rate": 1.0,
            "centroid_rmse_px": float("inf"),
            "n_gt": len(gt_coords),
            "n_detected": 0,
        }
    if len(gt_coords) == 0:
        return {
            "detection_rate": 0.0,
            "false_positive_rate": 1.0 if len(detected_coords) > 0 else 0.0,
            "centroid_rmse_px": 0.0,
            "n_gt": 0,
            "n_detected": len(detected_coords),
        }

    tree = cKDTree(detected_coords)
    dists, indices = tree.query(gt_coords)
    matched_gt = dists <= threshold_px
    n_matched_gt = matched_gt.sum()

    gt_tree = cKDTree(gt_coords)
    det_dists, _ = gt_tree.query(detected_coords)
    matched_det = det_dists <= threshold_px

    detection_rate = n_matched_gt / len(gt_coords)
    false_positive_rate = 1.0 - matched_det.sum() / len(detected_coords)
    centroid_rmse = float(np.sqrt(np.mean(dists[matched_gt] ** 2))) if n_matched_gt > 0 else 0.0

    return {
        "detection_rate": detection_rate,
        "false_positive_rate": false_positive_rate,
        "centroid_rmse_px": centroid_rmse,
        "n_gt": len(gt_coords),
        "n_detected": len(detected_coords),
        "n_matched": int(n_matched_gt),
    }


# ===================================================================
# Data loading helpers
# ===================================================================


def load_ground_truth_proportions(replicate_id: int, condition: str, sccube_dir: Path) -> pd.DataFrame:
    """Load ground truth cell type proportions per spot."""
    prop_path = sccube_dir / condition / "ST_sim" / f"Wu_ST_{replicate_id}_prop.csv"
    if not prop_path.exists():
        raise FileNotFoundError(f"Proportions not found: {prop_path}")
    df = pd.read_csv(prop_path, index_col=0)
    df = df.drop(columns=["spot_x", "spot_y"], errors="ignore")
    logger.info(
        "Loaded ground truth proportions: %d spots, %d cell types",
        len(df),
        len(df.columns),
    )
    return df


def load_cell_index(sccube_dir: Path, condition: str, replicate_id: int) -> pd.DataFrame:
    """Load cell-level simulation table with cell coordinates and spot labels."""
    index_path = sccube_dir / condition / "ST_sim" / f"Wu_ST_{replicate_id}_index.csv"
    if not index_path.exists():
        raise FileNotFoundError(f"Cell index not found: {index_path}")
    return pd.read_csv(index_path, index_col=0)


def compute_spot_coordinates(sccube_dir: Path, condition: str, replicate_id: int) -> pd.DataFrame:
    """Load spot coordinates from ST_sim meta file."""
    meta_path = sccube_dir / condition / "ST_sim" / f"Wu_ST_{replicate_id}_meta.csv"
    if not meta_path.exists():
        raise FileNotFoundError(f"Meta file not found: {meta_path}")
    df = pd.read_csv(meta_path, index_col=0)
    return df[["spot", "spot_x", "spot_y"]].drop_duplicates().set_index("spot")


def compute_image_space_spot_centers(
    spot_coords_df: pd.DataFrame,
    cell_index_df: pd.DataFrame,
    image_size: int,
    padding: int = 100,
) -> np.ndarray:
    """Map spot centres to image pixel coordinates (same transform as image gen)."""
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
    """Estimate spot radius in pixels from empirical cell-to-spot distances."""
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


def compute_gt_coords_px(cell_index_df: pd.DataFrame, image_size: int, padding: int = 100) -> np.ndarray:
    """Transform GT cell positions to pixel space (x, y)."""
    cell_xy = cell_index_df[["point_x", "point_y"]].to_numpy(dtype=float)
    x_min, x_max = cell_xy[:, 0].min(), cell_xy[:, 0].max()
    y_min, y_max = cell_xy[:, 1].min(), cell_xy[:, 1].max()
    scale_x = (image_size - 2 * padding) / (x_max - x_min + 1e-6)
    scale_y = (image_size - 2 * padding) / (y_max - y_min + 1e-6)
    return np.column_stack(
        [
            padding + (cell_xy[:, 0] - x_min) * scale_x,
            padding + (cell_xy[:, 1] - y_min) * scale_y,
        ]
    )


def load_synthetic_image(replicate_id: int, condition: str, mode: str, image_dir: Path) -> np.ndarray:
    """Load a synthetic StarDist-compatible image."""
    import cv2

    image_path = image_dir / condition / "images" / mode / f"Wu_rep_{replicate_id}_stardist.png"
    if not image_path.exists():
        raise FileNotFoundError(f"Image not found: {image_path}")
    logger.info("Loading image: %s", image_path)
    bgr = cv2.imread(str(image_path), cv2.IMREAD_COLOR)
    if bgr is None:
        raise ValueError(f"Failed to load image: {image_path}")
    rgb = cv2.cvtColor(bgr, cv2.COLOR_BGR2RGB)
    logger.info("Image shape: %s", rgb.shape)
    return rgb


def save_json(data: Any, path: Path) -> None:
    """Write *data* as indented JSON."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        json.dump(data, fh, indent=2, default=str)
    logger.info("Saved JSON: %s", path)


# ===================================================================
# Pipeline step functions
# ===================================================================


def run_segmentation_step(image: np.ndarray, mode: str) -> Tuple[np.ndarray, pd.DataFrame]:
    """Run StarDist segmentation, returning (masks, centroids_df).

    centroids_df columns: [y_pixel, x_pixel, nucleus_id]
    """
    from CITEgeist.model.morphology.segmentation import run_nuclei_segmentation

    modality = "dapi" if mode == "dapi" else "he"
    logger.info("Running StarDist segmentation (modality=%s) ...", modality)
    start = time.time()
    masks, centroids_df = run_nuclei_segmentation(image, modality=modality)
    elapsed = time.time() - start
    logger.info("Segmentation complete: %d nuclei in %.1fs", len(centroids_df), elapsed)
    return masks, centroids_df


def run_gt_ablation(
    cell_index_df: pd.DataFrame,
) -> Tuple[pd.Series, pd.DataFrame]:
    """Use GT cell positions as 'perfect segmentation'.

    Returns:
        nuclei_counts: per-spot nuclei count Series
        centroids_df: DataFrame with [y_pixel, x_pixel, nucleus_id] (dummy coords)
    """
    nuclei_counts = cell_index_df.groupby("spot").size()
    nuclei_counts.name = "nuclei_count"
    centroids_df = pd.DataFrame(
        {
            "y_pixel": 0.0,
            "x_pixel": 0.0,
            "nucleus_id": range(len(cell_index_df)),
        }
    )
    return nuclei_counts, centroids_df


def run_proportion_estimation(
    replicate_id: int,
    condition: str,
    sccube_dir: Path,
    output_dir: Path,
    cellularity_enhanced: bool = False,
    nuclei_counts: Optional[pd.Series] = None,
) -> Tuple[Any, pd.DataFrame]:
    """Run CITEgeist Module 3 continuous proportion estimation.

    Returns:
        model: CitegeistModel instance (needed for GEX pass)
        proportions_df: DataFrame (spots x cell types) with finetuned proportions
    """
    import scanpy as sc

    from CITEgeist.model.citegeist_model import CitegeistModel

    h5ad_dir = sccube_dir / condition / "h5ad_objects"
    cite_path = h5ad_dir / f"Wu_rep_{replicate_id}_CITE.h5ad"
    gex_path = h5ad_dir / f"Wu_rep_{replicate_id}_GEX.h5ad"
    for p in (cite_path, gex_path):
        if not p.exists():
            raise FileNotFoundError(f"h5ad not found: {p}")

    logger.info("Loading h5ad files: %s, %s", cite_path, gex_path)
    adata_cite = sc.read_h5ad(cite_path)
    adata_gex = sc.read_h5ad(gex_path)

    model_output = str(output_dir / condition / f"Wu_rep_{replicate_id}" / "citegeist")
    model = CitegeistModel(
        sample_name=f"Wu_rep_{replicate_id}",
        output_folder=model_output,
        simulation=True,
        gene_expression_adata=adata_gex,
        antibody_capture_adata=adata_cite,
    )

    # Preprocessing
    if cellularity_enhanced and nuclei_counts is not None:
        logger.info("Using cellularity-enhanced preprocessing")
        model.preprocess_antibody_discrete(scale_mode="per_marker")
        # Scale by nuclei counts so high-cellularity spots have larger signal
        aligned = nuclei_counts.reindex(model.antibody_capture_adata.obs_names).fillna(1)
        mean_n = aligned.mean()
        if mean_n > 0:
            scale = (aligned / mean_n).values
            model.antibody_capture_adata.X = model.antibody_capture_adata.X * scale[:, np.newaxis]
    else:
        model.preprocess_antibody()

    model.preprocess_gex(target_sum=10000)

    # Set cell profile dict directly (simulation profiles are known)
    model.cell_profile_dict = SIMULATION_CELL_PROFILE_DICT

    logger.info("Running Module 3 continuous proportions (skip_finetuning=True) ...")
    start = time.time()
    _global_props, finetuned_props = model.run_cell_proportion_model(method="qp", max_workers=8)
    elapsed = time.time() - start
    logger.info("Proportion estimation done in %.1fs", elapsed)

    return model, finetuned_props


def evaluate_proportions(pred_props: pd.DataFrame, gt_props: pd.DataFrame) -> Dict[str, Any]:
    """Evaluate predicted vs GT proportions: Pearson r, RMSE, JSD."""
    from scipy.spatial.distance import jensenshannon

    common_spots = pred_props.index.intersection(gt_props.index)
    common_types = pred_props.columns.intersection(gt_props.columns)

    if len(common_spots) == 0 or len(common_types) == 0:
        logger.warning("No common spots/types for proportion evaluation")
        return {"error": "no_overlap"}

    pred = pred_props.loc[common_spots, common_types].values
    gt = gt_props.loc[common_spots, common_types].values

    # Flatten Pearson r
    flat_pred = pred.flatten()
    flat_gt = gt.flatten()
    r, p = pearsonr(flat_pred, flat_gt)

    # RMSE
    rmse = float(np.sqrt(np.mean((pred - gt) ** 2)))

    # Per-spot JSD (mean across spots)
    jsd_list = []
    for i in range(len(common_spots)):
        p_row = pred[i]
        g_row = gt[i]
        # Ensure valid distributions (non-negative, sum > 0)
        p_row = np.clip(p_row, 0, None)
        g_row = np.clip(g_row, 0, None)
        p_sum = p_row.sum()
        g_sum = g_row.sum()
        if p_sum > 0 and g_sum > 0:
            jsd_list.append(float(jensenshannon(p_row / p_sum, g_row / g_sum) ** 2))
    mean_jsd = float(np.mean(jsd_list)) if jsd_list else float("nan")

    # Per-cell-type Pearson r
    per_type_r = {}
    for ct in common_types:
        ct_pred = pred_props.loc[common_spots, ct].values
        ct_gt = gt_props.loc[common_spots, ct].values
        if ct_gt.std() > 0 and ct_pred.std() > 0:
            ct_r, _ = pearsonr(ct_pred, ct_gt)
            per_type_r[ct] = float(ct_r)
        else:
            per_type_r[ct] = float("nan")

    metrics = {
        "pearson_r": float(r),
        "pearson_p": float(p),
        "rmse": rmse,
        "mean_jsd": mean_jsd,
        "n_spots": len(common_spots),
        "n_cell_types": len(common_types),
        "per_type_r": per_type_r,
    }
    logger.info(
        "Proportions -- Pearson r=%.4f  RMSE=%.4f  JSD=%.4f",
        r,
        rmse,
        mean_jsd,
    )
    return metrics


def evaluate_gex_deconvolution(
    model: Any,
    cell_counts_df: pd.DataFrame,
    cell_index_df: pd.DataFrame,
    replicate_id: int,
    condition: str,
    sccube_dir: Path,
) -> Dict[str, Any]:
    """Run GEX deconvolution and evaluate against pre-computed GT layers.

    Uses the original evaluation approach (benchmarking_gex.py):
    - GT from ST_GEX_sim/sample_N/layers/{CellType}_GT.csv (raw counts)
    - Predicted from model output layers (raw counts)
    - log1p() applied to BOTH before computing metrics
    - Metrics: RMSE, NRMSE, MAE, Pearson r (all on log1p scale)
    """
    logger.info("Running GEX deconvolution (SACE) ...")
    start = time.time()

    # Build cell_assignments and cell_spot_map from GT cell index (simulation GT cells)
    # cell_index_df columns: Cell, Cell_type, point_x, point_y, spot_x, spot_y, spot
    cell_assignments = dict(zip(cell_index_df["Cell"], cell_index_df["Cell_type"]))
    cell_spot_map = pd.DataFrame(
        {
            "cell_id": cell_index_df["Cell"].values,
            "spot_barcode": cell_index_df["spot"].values,
            "x": cell_index_df["point_x"].values,
            "y": cell_index_df["point_y"].values,
        }
    )

    try:
        model.run_cell_expression_pass1(
            cell_assignments=cell_assignments,
            cell_spot_map=cell_spot_map,
        )
    except Exception as exc:
        logger.error("GEX deconvolution failed: %s", exc)
        return {"error": str(exc)}
    elapsed = time.time() - start
    logger.info("GEX deconvolution done in %.1fs", elapsed)

    cell_type_names = list(SIMULATION_CELL_PROFILE_DICT.keys())

    # Load pre-computed GT layers from ST_GEX_sim (matches the original evaluation)
    # sample_N maps to Wu_rep_N
    gt_layers_dir = sccube_dir / condition / "ST_GEX_sim" / f"sample_{replicate_id}" / "layers"
    if not gt_layers_dir.exists():
        logger.warning("GT GEX layers not found: %s", gt_layers_dir)
        return {"error": f"gt_layers_not_found: {gt_layers_dir}"}

    from sklearn.metrics import mean_absolute_error, mean_squared_error

    per_type_metrics = {}
    for ct_name in cell_type_names:
        # Load GT layer (spots × genes, raw counts)
        gt_path = gt_layers_dir / f"{ct_name}_GT.csv"
        if not gt_path.exists():
            logger.warning("GT layer not found: %s", gt_path)
            per_type_metrics[ct_name] = {"pearson_r": float("nan"), "rmse": float("nan")}
            continue

        gt_df = pd.read_csv(gt_path, index_col=0)  # spots × genes

        # Load predicted layer (spots × genes, raw counts)
        # Layers saved by CitegeistModel at {output_folder}/{sample}_pass1/layers/pass1/
        model_out = Path(model.output_folder)
        pred_path = model_out / f"Wu_rep_{replicate_id}_pass1" / "layers" / "pass1" / f"{ct_name}_layer_pass1.csv"
        if not pred_path.exists():
            # Try underscore variant
            ct_underscore = ct_name.replace(" ", "_")
            pred_path = (
                model_out / f"Wu_rep_{replicate_id}_pass1" / "layers" / "pass1" / f"{ct_underscore}_layer_pass1.csv"
            )
        if not pred_path.exists():
            logger.warning("Pred layer not found: %s", pred_path)
            per_type_metrics[ct_name] = {"pearson_r": float("nan"), "rmse": float("nan")}
            continue

        pred_df = pd.read_csv(pred_path, index_col=0)  # spots × genes

        # Align genes and spots
        common_genes = gt_df.columns.intersection(pred_df.columns)
        common_spots = gt_df.index.intersection(pred_df.index)
        if len(common_genes) < 5 or len(common_spots) < 3:
            per_type_metrics[ct_name] = {"pearson_r": float("nan"), "rmse": float("nan")}
            continue

        gt_aligned = gt_df.loc[common_spots, common_genes]
        pred_aligned = pred_df.loc[common_spots, common_genes]

        # Apply log1p to BOTH (the original approach — line 37 of benchmarking_gex.py)
        gt_log = np.log1p(gt_aligned.values)
        pred_log = np.log1p(pred_aligned.values)

        # Compute metrics on log1p scale
        rmse = float(np.sqrt(mean_squared_error(gt_log, pred_log)))
        mae = float(mean_absolute_error(gt_log, pred_log))
        range_gt = gt_log.max() - gt_log.min()
        nrmse = rmse / range_gt if range_gt > 0 else float("nan")

        # Pearson r on log1p-flattened
        gt_flat = gt_log.flatten()
        pred_flat = pred_log.flatten()
        if gt_flat.std() > 0 and pred_flat.std() > 0:
            r, _ = pearsonr(gt_flat, pred_flat)
        else:
            r = float("nan")

        per_type_metrics[ct_name] = {
            "pearson_r": float(r),
            "rmse": rmse,
            "nrmse": nrmse,
            "mae": mae,
            "n_spots": len(common_spots),
            "n_genes": len(common_genes),
        }
        logger.info("  %s: r=%.4f RMSE=%.4f NRMSE=%.4f MAE=%.4f", ct_name, r, rmse, nrmse, mae)

    # Aggregate across cell types
    valid_rs = [v["pearson_r"] for v in per_type_metrics.values() if not np.isnan(v.get("pearson_r", float("nan")))]
    valid_rmses = [v["rmse"] for v in per_type_metrics.values() if not np.isnan(v.get("rmse", float("nan")))]
    valid_maes = [v["mae"] for v in per_type_metrics.values() if not np.isnan(v.get("mae", float("nan")))]

    metrics = {
        "mean_gex_pearson_r": float(np.mean(valid_rs)) if valid_rs else float("nan"),
        "mean_gex_rmse": float(np.mean(valid_rmses)) if valid_rmses else float("nan"),
        "median_gex_rmse": float(np.median(valid_rmses)) if valid_rmses else float("nan"),
        "mean_gex_mae": float(np.mean(valid_maes)) if valid_maes else float("nan"),
        "per_type": per_type_metrics,
        "n_shared_genes": len(common_genes) if "common_genes" in dir() else 0,
        "gex_time_sec": elapsed,
    }
    logger.info(
        "GEX: mean r=%.4f  mean RMSE=%.4f  mean MAE=%.4f (%d types)",
        metrics["mean_gex_pearson_r"],
        metrics["mean_gex_rmse"],
        metrics["mean_gex_mae"],
        len(valid_rs),
    )
    return metrics


# ===================================================================
# Main orchestrator
# ===================================================================


def run_benchmark(
    replicate_id: int,
    condition: str,
    mode: str,
    output_dir: Path,
    image_dir: Path = DEFAULT_IMAGE_DIR,
    sccube_dir: Path = DEFAULT_SCCUBE_DIR,
    cellularity_enhanced: bool = False,
) -> Dict[str, Any]:
    """Run the full unified benchmark pipeline.

    Args:
        replicate_id: Replicate index (0-4).
        condition: Simulation condition (e.g. "high_seg", "mixed").
        mode: One of "dapi", "h_and_e", "gt_ablation".
        output_dir: Root directory for outputs.
        image_dir: Directory containing synthetic images.
        sccube_dir: Root of scCube replicate data.
        cellularity_enhanced: If True, scale antibody by nuclei counts.

    Returns:
        Dictionary of all benchmark metrics.
    """
    logger.info("=" * 70)
    logger.info("UNIFIED BENCHMARK: rep=%d  cond=%s  mode=%s", replicate_id, condition, mode)
    logger.info("=" * 70)

    results: Dict[str, Any] = {
        "replicate_id": replicate_id,
        "condition": condition,
        "mode": mode,
        "cellularity_enhanced": cellularity_enhanced,
    }
    timings: Dict[str, float] = {}

    # ------------------------------------------------------------------
    # 1.  Load ground truth
    # ------------------------------------------------------------------
    cell_index_df = load_cell_index(sccube_dir, condition, replicate_id)
    gt_props = load_ground_truth_proportions(replicate_id, condition, sccube_dir)
    spot_coords_df = compute_spot_coordinates(sccube_dir, condition, replicate_id)

    # ------------------------------------------------------------------
    # 2.  Segmentation (or GT ablation)
    # ------------------------------------------------------------------
    if mode == "gt_ablation":
        nuclei_counts, _centroids_df = run_gt_ablation(cell_index_df)
        results["segmentation"] = {"mode": "gt_ablation"}
    else:
        # Load image and run StarDist
        image = load_synthetic_image(replicate_id, condition, mode, image_dir)
        image_size = image.shape[0]

        seg_start = time.time()
        masks, centroids_df = run_segmentation_step(image, mode)
        timings["segmentation_sec"] = time.time() - seg_start

        # Assign detected centroids to spots
        # CRITICAL: extract (x, y) -- centroids_df columns are [y_pixel, x_pixel, ...]
        detected_xy = centroids_df[["x_pixel", "y_pixel"]].values

        spot_centers_px = compute_image_space_spot_centers(spot_coords_df, cell_index_df, image_size)
        spot_radius_px = estimate_spot_radius_px(spot_coords_df, cell_index_df, image_size)

        from CITEgeist.model.morphology.segmentation import assign_nuclei_centroids_to_spots

        nuclei_counts = assign_nuclei_centroids_to_spots(
            centroids_xy=detected_xy,
            spot_centers_xy=spot_centers_px,
            spot_radius_px=spot_radius_px,
            spot_names=spot_coords_df.index,
        )

        results["spot_radius_px"] = float(spot_radius_px)

        # Segmentation evaluation: NN matching against GT cell positions
        gt_coords_px = compute_gt_coords_px(cell_index_df, image_size)
        seg_metrics = evaluate_segmentation_nn(
            gt_coords=gt_coords_px,
            detected_coords=detected_xy,
            threshold_px=15.0,
        )
        results["segmentation"] = seg_metrics
        logger.info(
            "Segmentation: detection=%.2f  FP=%.2f  RMSE=%.1f px",
            seg_metrics["detection_rate"],
            seg_metrics["false_positive_rate"],
            seg_metrics["centroid_rmse_px"],
        )

        # Nuclei count accuracy (correlation with GT per-spot counts)
        gt_spot_counts = cell_index_df.groupby("spot").size()
        common_count_spots = nuclei_counts.index.intersection(gt_spot_counts.index)
        if len(common_count_spots) > 2:
            count_r, _ = pearsonr(
                nuclei_counts.loc[common_count_spots].values,
                gt_spot_counts.loc[common_count_spots].values,
            )
            count_rmse = float(
                np.sqrt(
                    np.mean(
                        (nuclei_counts.loc[common_count_spots].values - gt_spot_counts.loc[common_count_spots].values)
                        ** 2
                    )
                )
            )
            results["nuclei_count_correlation"] = float(count_r)
            results["nuclei_count_rmse"] = count_rmse

    results["total_nuclei_detected"] = int(nuclei_counts.sum())

    # Save nuclei counts
    out_rep_dir = output_dir / condition / mode / f"Wu_rep_{replicate_id}"
    out_rep_dir.mkdir(parents=True, exist_ok=True)
    nuclei_counts.to_csv(out_rep_dir / "nuclei_counts.csv")

    # ------------------------------------------------------------------
    # 3.  Continuous proportion estimation (Module 3)
    # ------------------------------------------------------------------
    prop_start = time.time()
    model, proportions_df = run_proportion_estimation(
        replicate_id=replicate_id,
        condition=condition,
        sccube_dir=sccube_dir,
        output_dir=output_dir,
        cellularity_enhanced=cellularity_enhanced,
        nuclei_counts=nuclei_counts,
    )
    timings["proportion_estimation_sec"] = time.time() - prop_start

    # Save continuous proportions
    proportions_df.to_csv(out_rep_dir / "cell_proportions.csv")

    # Evaluate proportions
    prop_metrics = evaluate_proportions(proportions_df, gt_props)
    results["proportion_metrics"] = prop_metrics

    # ------------------------------------------------------------------
    # 4.  Discretize proportions -> integer cell counts
    # ------------------------------------------------------------------
    cell_counts_df = discretize_proportions_df(proportions_df, nuclei_counts)
    cell_counts_df.to_csv(out_rep_dir / "cell_counts.csv")

    # ------------------------------------------------------------------
    # 5.  Hungarian single-cell accuracy
    # ------------------------------------------------------------------
    sc_accuracy = compute_hungarian_sc_accuracy_allspots(cell_index_df, cell_counts_df)
    results["hungarian_sc_accuracy"] = float(sc_accuracy)
    logger.info("Hungarian SC accuracy: %.4f", sc_accuracy)

    # ------------------------------------------------------------------
    # 6.  GEX deconvolution + evaluation
    # ------------------------------------------------------------------
    gex_start = time.time()
    gex_metrics = evaluate_gex_deconvolution(
        model=model,
        cell_counts_df=cell_counts_df,
        cell_index_df=cell_index_df,
        replicate_id=replicate_id,
        condition=condition,
        sccube_dir=sccube_dir,
    )
    timings["gex_deconv_sec"] = time.time() - gex_start
    results["gex_metrics"] = gex_metrics

    # ------------------------------------------------------------------
    # Wrap up
    # ------------------------------------------------------------------
    results["timings"] = timings
    save_json(results, out_rep_dir / "benchmark_results.json")
    logger.info("All results saved to %s", out_rep_dir)
    return results


# ===================================================================
# CLI entry point
# ===================================================================


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Unified simulation benchmark: StarDist + continuous proportions + SC evaluation"
    )
    parser.add_argument("--replicate-id", type=int, required=True, help="Replicate ID (0-4)")
    parser.add_argument(
        "--condition",
        type=str,
        required=True,
        help="Simulation condition (e.g. high_seg, mixed)",
    )
    parser.add_argument(
        "--mode",
        choices=["dapi", "h_and_e", "gt_ablation"],
        required=True,
        help="Segmentation mode or GT ablation",
    )
    parser.add_argument("--output-dir", type=Path, required=True, help="Root output directory")
    parser.add_argument("--image-dir", type=Path, default=DEFAULT_IMAGE_DIR, help="Image directory")
    parser.add_argument("--sccube-dir", type=Path, default=DEFAULT_SCCUBE_DIR, help="scCube data root")
    parser.add_argument(
        "--cellularity-enhanced",
        action="store_true",
        default=False,
        help="Scale antibody signal by nuclei counts",
    )

    args = parser.parse_args()

    results = run_benchmark(
        replicate_id=args.replicate_id,
        condition=args.condition,
        mode=args.mode,
        output_dir=args.output_dir,
        image_dir=args.image_dir,
        sccube_dir=args.sccube_dir,
        cellularity_enhanced=args.cellularity_enhanced,
    )

    print("\n=== RESULTS SUMMARY ===")
    pm = results.get("proportion_metrics", {})
    print(f"Proportion Pearson r: {pm.get('pearson_r', 'N/A')}")
    print(f"Proportion RMSE:      {pm.get('rmse', 'N/A')}")
    print(f"Hungarian SC acc:     {results.get('hungarian_sc_accuracy', 'N/A')}")
    gm = results.get("gex_metrics", {})
    print(f"GEX mean Pearson r:   {gm.get('mean_gex_pearson_r', 'N/A')}")


if __name__ == "__main__":
    main()
