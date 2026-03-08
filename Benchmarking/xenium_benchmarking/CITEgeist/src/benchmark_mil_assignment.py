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

# Ground truth files
CELL_TYPES_PATH = DATA_DIR / "cell_type_assignments.csv"
CELL_TO_SPOT_PATH = DATA_DIR / "cell_to_spot_mapping.csv"
GT_PROPS_DIR = PROJECT_ROOT / "Benchmarking" / "xenium_benchmarking" / "ground_truth" / "proportions"


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


def evaluate(
    assignments: dict,
    gt_props: pd.DataFrame,
    gt_cells: pd.DataFrame,
    nuclei_spot_map: pd.DataFrame,
    cell_types: list,
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

    # 2. Per-spot accuracy against GT cell-level assignments
    # Note: Cellpose nucleus IDs != Xenium cell barcodes, so we compare
    # at spot level: GT proportions from single-cell assignments vs predicted
    if gt_cells is not None:
        gt_spot_props_sc = []
        pred_spot_props_sc = []
        for spot_id in gt_cells['spot_id'].unique():
            spot_gt = gt_cells[gt_cells['spot_id'] == spot_id]
            spot_nuclei = nuclei_spot_map[nuclei_spot_map['spot_id'] == spot_id]
            if len(spot_gt) == 0 or len(spot_nuclei) == 0:
                continue
            # GT proportions from single-cell assignments
            gt_counts = np.zeros(n_types)
            for ct in spot_gt['cell_type']:
                if ct in cell_types:
                    gt_counts[cell_types.index(ct)] += 1
            gt_prop_sc = gt_counts / gt_counts.sum() if gt_counts.sum() > 0 else gt_counts
            # Predicted proportions
            pred_counts = np.zeros(n_types)
            for nid in spot_nuclei['nucleus_id']:
                if nid in assignments:
                    ct = assignments[nid]
                    if ct in cell_types:
                        pred_counts[cell_types.index(ct)] += 1
            pred_prop_sc = pred_counts / pred_counts.sum() if pred_counts.sum() > 0 else pred_counts
            gt_spot_props_sc.append(gt_prop_sc)
            pred_spot_props_sc.append(pred_prop_sc)

        if gt_spot_props_sc:
            gt_arr_sc = np.array(gt_spot_props_sc)
            pred_arr_sc = np.array(pred_spot_props_sc)
            r_sc, p_sc = pearsonr(gt_arr_sc.flatten(), pred_arr_sc.flatten())
            metrics['sc_gt_proportion_r'] = float(r_sc)
            metrics['sc_gt_proportion_p'] = float(p_sc)
            metrics['sc_gt_spots_evaluated'] = len(gt_spot_props_sc)

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

    # Evaluate
    metrics = evaluate(
        result.assignments, gt_props, gt_cells,
        nuclei_spot_map, CELL_TYPES,
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
