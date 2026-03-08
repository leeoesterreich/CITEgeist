"""Benchmark Module 3b MIL single-cell assignment on Xenium pseudo-Visium.

Pipeline:
1. Load Module 3 hybrid proportions (pre-computed)
2. Load SimCLR backbone (pre-trained)
3. Extract DAPI+boundary embeddings via backbone
4. Train MIL head on proportions
5. Proportion-weighted Hungarian assignment
6. Evaluate: proportion r, single-cell accuracy
"""
import argparse
import json
import logging
import sys
from pathlib import Path
from typing import Dict

import numpy as np
import pandas as pd
import torch

# Add project root to path
PROJECT_ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(PROJECT_ROOT))

from CITEgeist.model.morphology_backbone import DAPIBackbone
from CITEgeist.model.single_cell_mil import SingleCellMIL, train_mil
from CITEgeist.model.hungarian_assignment import assign_nuclei_to_types
from CITEgeist.model.morphology_features import largest_remainder_discretize

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s')
logger = logging.getLogger(__name__)

# Xenium benchmark cell types (achievable-7)
CELL_TYPES = [
    "B cells", "CD4+ T cells", "CD8+ T cells",
    "Macrophages", "Endothelial", "Epithelial", "Fibroblasts",
]

# Paths
BENCHMARK_ROOT = PROJECT_ROOT / "Benchmarking" / "xenium_benchmarking" / "CITEgeist"
PSEUDOVISIUM_DIR = PROJECT_ROOT / "Benchmarking" / "xenium_pseudovisium"
DATA_DIR = PSEUDOVISIUM_DIR / "data_protein_gt"
HYBRID_OUTPUT_DIR = BENCHMARK_ROOT / "output" / "hybrid_cellpose"
PATCHES_DIR = BENCHMARK_ROOT / "output" / "patches_v2"
XENIUM_DATA_DIR = Path("/ix1/alee/LO_LAB/General/Public_Data/10x_PublicData/Xenium_RNA_Proteomic_RenalCellCarcinoma")

# Ground truth files
CELL_TYPES_PATH = DATA_DIR / "cell_type_assignments.csv"
CELL_TO_SPOT_PATH = DATA_DIR / "cell_to_spot_mapping.csv"
GT_PROPS_DIR = PROJECT_ROOT / "Benchmarking" / "xenium_benchmarking" / "ground_truth" / "proportions"

# Pixel size for Xenium morphology images (from experiment.xenium)
PIXEL_SIZE_UM = 0.2125

# Region boundaries in microns (from prepare_patches.py)
REGION_BOUNDS_UM = {
    0: (29.01, 2279.01, 30.87, 5486.83),
    1: (2329.01, 4579.01, 30.87, 5400.23),
    2: (4629.01, 6829.01, 30.87, 5486.83),
    3: (6879.01, 9129.01, 30.87, 5660.04),
    4: (9179.01, 11429.01, 30.87, 5746.64),
}
PADDING_UM = 100.0


def load_stage1_results(region_id: int) -> pd.DataFrame:
    """Load hybrid proportions from Module 3.

    Returns DataFrame with spot_id index and cell type columns.
    """
    region_name = f"Xenium_region_{region_id}"
    prop_path = HYBRID_OUTPUT_DIR / region_name / f"{region_name}_deconv_predictions.csv"
    if not prop_path.exists():
        raise FileNotFoundError(f"Stage 1 results not found: {prop_path}")
    props = pd.read_csv(prop_path, index_col=0)
    return props


def load_ground_truth(region_id: int):
    """Load ground truth proportions and single-cell assignments.

    Returns:
        gt_props: DataFrame with spot_id index and cell type columns
        gt_cells: DataFrame with columns [spot_id, cell_barcode, cell_type] or None
    """
    region_name = f"Xenium_region_{region_id}"

    # Spot-level proportions
    gt_props_path = GT_PROPS_DIR / f"{region_name}_prop.csv"
    gt_props = pd.read_csv(gt_props_path, index_col=0)

    # Single-cell assignments (global files, filter by region)
    gt_cells = None
    if CELL_TYPES_PATH.exists() and CELL_TO_SPOT_PATH.exists():
        cell_types_df = pd.read_csv(CELL_TYPES_PATH, index_col=0)
        mapping_df = pd.read_csv(CELL_TO_SPOT_PATH, index_col=0)

        # Filter to this region
        region_mapping = mapping_df[
            (mapping_df['region_id'] == region_id) &
            (mapping_df['spot_idx'] != -1)
        ].copy()

        if len(region_mapping) > 0:
            region_mapping['cell_barcode'] = region_mapping.index
            gt_cells = region_mapping.merge(cell_types_df, left_index=True, right_index=True)
            gt_cells = gt_cells[['spot_id', 'cell_barcode', 'cell_type']]
            logger.info("Loaded %d GT cells for region %d", len(gt_cells), region_id)

    return gt_props, gt_cells


def load_patches(region_id: int) -> dict:
    """Load pre-extracted nucleus patches per spot.

    Patch files: spot_spot_XXXX_patches.npy (where spot_id = "spot_XXXX")
    Nucleus IDs: spot_spot_XXXX_nucleus_ids.npy (Xenium cell barcodes)

    Returns:
        dict mapping spot_id -> dict with keys 'patches' (ndarray) and 'nucleus_ids' (ndarray)
    """
    region_dir = PATCHES_DIR / f"region_{region_id}"
    if not region_dir.exists():
        raise FileNotFoundError(f"Patches directory not found: {region_dir}")

    spot_data = {}
    for f in sorted(region_dir.glob("*_patches.npy")):
        # File: spot_spot_1032_patches.npy -> spot_id = "spot_1032"
        stem = f.stem.replace("_patches", "")  # "spot_spot_1032"
        spot_id = stem.replace("spot_", "", 1)  # "spot_1032"

        patches = np.load(f)
        if len(patches) == 0:
            continue

        # Load companion nucleus IDs
        nid_path = f.parent / f"{stem}_nucleus_ids.npy"
        nucleus_ids = np.load(nid_path) if nid_path.exists() else np.arange(len(patches))

        spot_data[spot_id] = {
            'patches': patches,
            'nucleus_ids': nucleus_ids,
        }

    return spot_data


def build_cellpose_to_xenium_map(
    region_id: int,
    cellpose_nucleus_ids: np.ndarray,
    max_distance_um: float = 5.0,
) -> Dict[int, str]:
    """Map Cellpose nucleus IDs to Xenium cell barcodes by spatial proximity.

    Cellpose centroids are in pixel coordinates of the cropped region image.
    Xenium cell centroids are in microns (global coordinate system).

    Args:
        region_id: Region index (0-4)
        cellpose_nucleus_ids: Array of Cellpose nucleus IDs from nucleus_features.csv
        max_distance_um: Maximum matching distance in microns

    Returns:
        Dict mapping Cellpose nucleus_id (int) -> Xenium barcode (str)
    """
    from scipy.spatial import cKDTree

    # Load Cellpose centroids
    features_path = PATCHES_DIR / f"region_{region_id}" / "nucleus_features.csv"
    cp_df = pd.read_csv(features_path)

    # Compute crop origin: region bounds with padding, clamped to 0
    x_min_um = max(0, REGION_BOUNDS_UM[region_id][0] - PADDING_UM)
    y_min_um = max(0, REGION_BOUNDS_UM[region_id][2] - PADDING_UM)

    # Convert Cellpose pixel coords to microns (global)
    cp_x = cp_df['centroid_x'].values * PIXEL_SIZE_UM + x_min_um * PIXEL_SIZE_UM
    cp_y = cp_df['centroid_y'].values * PIXEL_SIZE_UM + y_min_um * PIXEL_SIZE_UM

    # Wait — x_min_um is already in microns after clamping; the pixel origin
    # is int(x_min_um / PIXEL_SIZE_UM) pixels. So:
    # pixel_in_crop * PIXEL_SIZE = offset from crop origin in microns
    # global_micron = pixel_in_crop * PIXEL_SIZE + crop_origin_microns
    # crop_origin_microns = max(0, x_min_um - PADDING_UM) but x_min_um is already padded...
    # Actually the crop is: pixel bounds = micron_to_pixel(x_min_um - PADDING, y_min_um - PADDING)
    # clamped to 0. So crop_origin_px = max(0, int((x_min - pad) / pixel_size))
    # And global_micron = (crop_origin_px + pixel_in_crop) * pixel_size
    x_min_crop_um = REGION_BOUNDS_UM[region_id][0] - PADDING_UM
    y_min_crop_um = REGION_BOUNDS_UM[region_id][2] - PADDING_UM
    x_origin_px = max(0, int(x_min_crop_um / PIXEL_SIZE_UM))
    y_origin_px = max(0, int(y_min_crop_um / PIXEL_SIZE_UM))

    cp_global_x = (x_origin_px + cp_df['centroid_x'].values) * PIXEL_SIZE_UM
    cp_global_y = (y_origin_px + cp_df['centroid_y'].values) * PIXEL_SIZE_UM

    # Load ALL Xenium cells in spatial region (not just GT-mapped)
    xenium_cells = pd.read_parquet(XENIUM_DATA_DIR / "cells.parquet")
    x_lo = REGION_BOUNDS_UM[region_id][0] - PADDING_UM
    x_hi = REGION_BOUNDS_UM[region_id][1] + PADDING_UM
    y_lo = REGION_BOUNDS_UM[region_id][2] - PADDING_UM
    y_hi = REGION_BOUNDS_UM[region_id][3] + PADDING_UM
    region_xenium = xenium_cells[
        (xenium_cells['x_centroid'] >= max(0, x_lo)) & (xenium_cells['x_centroid'] <= x_hi) &
        (xenium_cells['y_centroid'] >= max(0, y_lo)) & (xenium_cells['y_centroid'] <= y_hi)
    ].copy()

    # Build KD-tree on Xenium coords, query with Cellpose
    xen_coords = np.column_stack([region_xenium['x_centroid'].values, region_xenium['y_centroid'].values])
    cp_coords = np.column_stack([cp_global_x, cp_global_y])
    tree = cKDTree(xen_coords)
    dists, idxs = tree.query(cp_coords, k=1)

    # Build mapping for matches within threshold
    cp_to_xenium = {}
    n_matched = 0
    for i, (d, idx) in enumerate(zip(dists, idxs)):
        if d <= max_distance_um:
            cp_nid = int(cp_df['nucleus_id'].iloc[i])
            xen_barcode = region_xenium['cell_id'].iloc[idx]
            cp_to_xenium[cp_nid] = xen_barcode
            n_matched += 1

    logger.info(
        "Spatial matching: %d/%d Cellpose nuclei matched to Xenium (%.1f%%, thresh=%.1fµm, "
        "median dist=%.2fµm)",
        n_matched, len(cp_df), n_matched / len(cp_df) * 100,
        max_distance_um, float(np.median(dists)),
    )
    return cp_to_xenium


def evaluate(
    assignments: dict,
    gt_props: pd.DataFrame,
    gt_cells: pd.DataFrame,
    nuclei_spot_map: pd.DataFrame,
    cell_types: list,
    cellpose_to_xenium: Dict[int, str] = None,
    gt_types_df: pd.DataFrame = None,
) -> dict:
    """Compute all evaluation metrics."""
    from scipy.stats import pearsonr
    from sklearn.metrics import mean_squared_error

    metrics = {}
    n_types = len(cell_types)

    # 1. Proportion correlation: re-aggregate assignments to proportions
    pred_spot_props = []
    gt_spot_props = []
    for spot_id in gt_props.index:
        spot_nuclei = nuclei_spot_map[nuclei_spot_map['spot_id'] == spot_id]
        n = len(spot_nuclei)
        if n == 0:
            continue
        type_counts = np.zeros(n_types)
        for nid in spot_nuclei['nucleus_id']:
            if nid in assignments:
                ct = assignments[nid]
                if ct in cell_types:
                    type_counts[cell_types.index(ct)] += 1
        pred_spot_props.append(type_counts / max(n, 1))
        gt_spot_props.append(gt_props.loc[spot_id, cell_types].values.astype(float))

    if pred_spot_props:
        pred_arr = np.array(pred_spot_props)
        gt_arr = np.array(gt_spot_props)
        pred_flat = pred_arr.flatten()
        gt_flat = gt_arr.flatten()
        r, p = pearsonr(pred_flat, gt_flat)
        metrics['proportion_pearson_r'] = float(r)
        metrics['proportion_pearson_p'] = float(p)
        metrics['proportion_rmse'] = float(np.sqrt(mean_squared_error(gt_arr, pred_arr)))

        # Per-type correlations
        per_type_corr = {}
        for i, ct in enumerate(cell_types):
            gt_col = gt_arr[:, i]
            pred_col = pred_arr[:, i]
            if gt_col.std() > 0 and pred_col.std() > 0:
                r_ct, p_ct = pearsonr(gt_col, pred_col)
                per_type_corr[ct] = {'pearson_r': float(r_ct), 'p_value': float(p_ct)}
            else:
                per_type_corr[ct] = {'pearson_r': 0.0, 'p_value': 1.0}
        metrics['per_type_correlations'] = per_type_corr

    # 2. Single-cell accuracy via spatial matching
    if cellpose_to_xenium is not None and gt_types_df is not None:
        correct, total = 0, 0
        per_type_correct = {ct: 0 for ct in cell_types}
        per_type_total = {ct: 0 for ct in cell_types}

        for cp_nid, pred_type in assignments.items():
            if cp_nid not in cellpose_to_xenium:
                continue
            xen_barcode = cellpose_to_xenium[cp_nid]
            if xen_barcode not in gt_types_df.index:
                continue
            gt_type = gt_types_df.loc[xen_barcode, 'cell_type']
            if gt_type not in cell_types:
                continue
            total += 1
            per_type_total[gt_type] = per_type_total.get(gt_type, 0) + 1
            if pred_type == gt_type:
                correct += 1
                per_type_correct[gt_type] = per_type_correct.get(gt_type, 0) + 1

        if total > 0:
            metrics['single_cell_accuracy'] = correct / total
            metrics['single_cell_total'] = total
            metrics['single_cell_correct'] = correct
            metrics['random_baseline'] = 1.0 / n_types
            metrics['per_type_accuracy'] = {
                ct: per_type_correct[ct] / max(per_type_total[ct], 1)
                for ct in cell_types if per_type_total[ct] > 0
            }
            metrics['per_type_total'] = {
                ct: per_type_total[ct]
                for ct in cell_types if per_type_total[ct] > 0
            }
            logger.info(
                "Single-cell accuracy: %d/%d = %.1f%% (random=%.1f%%)",
                correct, total, correct / total * 100, 100.0 / n_types,
            )

    return metrics


def main():
    parser = argparse.ArgumentParser(description="Benchmark MIL single-cell assignment")
    parser.add_argument("--region", type=int, required=True, help="Region index (0-4)")
    parser.add_argument("--simclr-checkpoint", type=str, required=True)
    parser.add_argument("--output-dir", type=str, required=True)
    parser.add_argument("--lambda-prior", type=float, default=1.0)
    parser.add_argument("--n-epochs", type=int, default=100)
    parser.add_argument("--device", type=str, default="cuda")
    parser.add_argument("--batch-size", type=int, default=64)
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    device = args.device if torch.cuda.is_available() else 'cpu'
    logger.info("=== MIL Benchmark: Region %d (device=%s) ===", args.region, device)

    # Load data
    stage1_props = load_stage1_results(args.region)
    gt_props, gt_cells = load_ground_truth(args.region)
    spot_data = load_patches(args.region)

    logger.info("Loaded %d spots with patches, %d stage1 spots",
                len(spot_data), len(stage1_props))

    # Initialize backbone
    backbone = DAPIBackbone(checkpoint=args.simclr_checkpoint, device=device)

    # Extract embeddings and build nuclei_spot_map with real Xenium cell IDs
    embeddings = {}
    nuclei_records = []
    for spot_id, data in spot_data.items():
        patches = data['patches']
        nucleus_ids = data['nucleus_ids']
        emb = backbone.extract_numpy(patches, batch_size=args.batch_size, device=device)
        embeddings[spot_id] = emb
        for nid in nucleus_ids:
            nuclei_records.append({'nucleus_id': int(nid), 'spot_id': spot_id})

    nuclei_spot_map = pd.DataFrame(nuclei_records)
    nuclei_counts = nuclei_spot_map.groupby('spot_id').size()

    logger.info("Extracted embeddings: %d spots, %d nuclei", len(embeddings), len(nuclei_spot_map))

    # Prepare proportions DataFrame for run_nucleus_assignment_mil
    # Needs 'spot_id' column + cell type columns
    props_for_mil = stage1_props[CELL_TYPES].copy()
    props_for_mil['spot_id'] = props_for_mil.index

    # Run MIL assignment
    from CITEgeist.model.module3b_nucleus_assignment import run_nucleus_assignment_mil

    result = run_nucleus_assignment_mil(
        embeddings=embeddings,
        nuclei_spot_map=nuclei_spot_map,
        proportions=props_for_mil,
        nuclei_counts=nuclei_counts,
        cell_types=CELL_TYPES,
        n_epochs=args.n_epochs,
        lambda_prior=args.lambda_prior,
        device=device,
    )

    logger.info("Assignment complete: %d nuclei assigned", len(result.assignments))

    # Build Cellpose → Xenium spatial mapping for single-cell accuracy
    all_cp_nids = nuclei_spot_map['nucleus_id'].values
    cellpose_to_xenium = build_cellpose_to_xenium_map(
        args.region, all_cp_nids, max_distance_um=5.0,
    )
    gt_types_df = pd.read_csv(CELL_TYPES_PATH, index_col=0)

    # Evaluate
    metrics = evaluate(
        result.assignments, gt_props, gt_cells,
        nuclei_spot_map, CELL_TYPES,
        cellpose_to_xenium=cellpose_to_xenium,
        gt_types_df=gt_types_df,
    )

    logger.info("Results: prop_r=%.4f, sc_acc=%.4f",
                metrics.get('proportion_pearson_r', -1),
                metrics.get('single_cell_accuracy', -1))

    # Save results
    region_output = output_dir / f"region_{args.region}"
    region_output.mkdir(parents=True, exist_ok=True)

    metrics_path = region_output / "metrics.json"
    with open(metrics_path, 'w') as f:
        json.dump(metrics, f, indent=2)

    assignments_path = region_output / "assignments.csv"
    pd.DataFrame([
        {'nucleus_id': nid, 'cell_type': ct}
        for nid, ct in result.assignments.items()
    ]).to_csv(assignments_path, index=False)

    logger.info("Saved results to %s", region_output)


if __name__ == "__main__":
    main()
