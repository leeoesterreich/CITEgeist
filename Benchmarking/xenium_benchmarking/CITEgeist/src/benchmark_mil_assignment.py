"""Benchmark Module 3b MIL single-cell assignment on Xenium pseudo-Visium.

Pipeline:
1. Load Module 3 hybrid proportions (pre-computed)
2. Load SimCLR backbone (pre-trained)
3. Extract DAPI+boundary embeddings via backbone
4. Train MIL head on proportions
5. Proportion-weighted Hungarian assignment
6. Evaluate: proportion r, single-cell accuracy, GEX r
"""
import argparse
import json
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd

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

DATA_ROOT = Path(PROJECT_ROOT / "Benchmarking" / "xenium_pseudovisium" / "data_protein_gt")


def load_stage1_results(region_id: int) -> pd.DataFrame:
    """Load hybrid post-filter proportions from Module 3."""
    results_dir = PROJECT_ROOT / "Benchmarking" / "xenium_benchmarking" / "CITEgeist" / "output" / "hybrid_cellpose"
    region_name = f"Xenium_region_{region_id}"
    prop_path = results_dir / region_name / f"{region_name}_deconv_predictions.csv"
    if not prop_path.exists():
        raise FileNotFoundError(f"Stage 1 results not found: {prop_path}")
    return pd.read_csv(prop_path)


def load_ground_truth(region_id: int):
    """Load ground truth proportions and single-cell assignments."""
    region_name = f"Xenium_region_{region_id}"

    # Spot-level proportions
    gt_props_path = DATA_ROOT / "ground_truth" / f"{region_name}_proportions.csv"
    gt_props = pd.read_csv(gt_props_path)

    # Single-cell assignments
    gt_cells_path = DATA_ROOT / "ground_truth" / f"{region_name}_cell_assignments.csv"
    gt_cells = pd.read_csv(gt_cells_path) if gt_cells_path.exists() else None

    return gt_props, gt_cells


def load_patches(region_id: int) -> dict:
    """Load pre-extracted nucleus patches per spot."""
    patches_dir = PROJECT_ROOT / "Benchmarking" / "xenium_benchmarking" / "CITEgeist" / "output" / "patches"
    region_name = f"Xenium_region_{region_id}"
    region_dir = patches_dir / region_name

    spot_patches = {}
    for f in sorted(region_dir.glob("*_patches.npy")):
        spot_id = f.stem.replace("_patches", "")
        patches = np.load(f)
        if len(patches) > 0:
            spot_patches[spot_id] = patches

    return spot_patches


def evaluate(
    assignments: dict,
    proportions: pd.DataFrame,
    gt_props: pd.DataFrame,
    gt_cells: pd.DataFrame,
    nuclei_spot_map: pd.DataFrame,
    cell_types: list,
) -> dict:
    """Compute all evaluation metrics."""
    from scipy.stats import pearsonr

    metrics = {}

    # 1. Proportion correlation: re-aggregate assignments to proportions
    pred_spot_props = []
    for spot_id in gt_props['spot_id'].unique():
        spot_nuclei = nuclei_spot_map[nuclei_spot_map['spot_id'] == spot_id]
        n = len(spot_nuclei)
        if n == 0:
            continue
        type_counts = np.zeros(len(cell_types))
        for nid in spot_nuclei['nucleus_id']:
            if nid in assignments:
                ct = assignments[nid]
                if ct in cell_types:
                    type_counts[cell_types.index(ct)] += 1
        pred_spot_props.append(type_counts / max(n, 1))

    if pred_spot_props:
        pred_flat = np.array(pred_spot_props).flatten()
        gt_flat = gt_props[cell_types].values.flatten()
        min_len = min(len(pred_flat), len(gt_flat))
        r, p = pearsonr(pred_flat[:min_len], gt_flat[:min_len])
        metrics['proportion_pearson_r'] = float(r)
        metrics['proportion_pearson_p'] = float(p)

    # 2. Single-cell accuracy (if ground truth available)
    if gt_cells is not None:
        gt_map = dict(zip(gt_cells['nucleus_id'], gt_cells['cell_type']))
        correct, total = 0, 0
        per_type_correct = {ct: 0 for ct in cell_types}
        per_type_total = {ct: 0 for ct in cell_types}
        for nid, pred_type in assignments.items():
            if nid in gt_map:
                gt_type = gt_map[nid]
                total += 1
                if gt_type in per_type_total:
                    per_type_total[gt_type] += 1
                if pred_type == gt_type:
                    correct += 1
                    if gt_type in per_type_correct:
                        per_type_correct[gt_type] += 1

        metrics['single_cell_accuracy'] = correct / max(total, 1)
        metrics['single_cell_total'] = total
        metrics['per_type_accuracy'] = {
            ct: per_type_correct[ct] / max(per_type_total[ct], 1)
            for ct in cell_types
        }

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

    logger.info("=== MIL Benchmark: Region %d ===", args.region)

    # Load data
    stage1_props = load_stage1_results(args.region)
    gt_props, gt_cells = load_ground_truth(args.region)
    spot_patches = load_patches(args.region)

    logger.info("Loaded %d spots with patches", len(spot_patches))

    # Initialize backbone
    backbone = DAPIBackbone(checkpoint=args.simclr_checkpoint, device=args.device)

    # Extract embeddings
    embeddings = {}
    nuclei_records = []
    for spot_id, patches in spot_patches.items():
        emb = backbone.extract_numpy(patches, batch_size=args.batch_size, device=args.device)
        embeddings[spot_id] = emb
        for i in range(len(patches)):
            nuclei_records.append({'nucleus_id': f"{spot_id}_{i}", 'spot_id': spot_id})

    nuclei_spot_map = pd.DataFrame(nuclei_records)
    nuclei_counts = nuclei_spot_map.groupby('spot_id').size()

    logger.info("Extracted embeddings: %d spots, %d nuclei", len(embeddings), len(nuclei_spot_map))

    # Run MIL assignment
    from CITEgeist.model.module3b_nucleus_assignment import run_nucleus_assignment_mil

    result = run_nucleus_assignment_mil(
        embeddings=embeddings,
        nuclei_spot_map=nuclei_spot_map,
        proportions=stage1_props,
        nuclei_counts=nuclei_counts,
        cell_types=CELL_TYPES,
        n_epochs=args.n_epochs,
        lambda_prior=args.lambda_prior,
        device=args.device,
    )

    logger.info("Assignment complete: %d nuclei assigned", len(result.assignments))

    # Evaluate
    metrics = evaluate(
        result.assignments, stage1_props, gt_props, gt_cells,
        nuclei_spot_map, CELL_TYPES,
    )

    logger.info("Results: prop_r=%.4f, sc_acc=%.4f",
                metrics.get('proportion_pearson_r', -1),
                metrics.get('single_cell_accuracy', -1))

    # Save results
    metrics_path = output_dir / f"region_{args.region}_metrics.json"
    with open(metrics_path, 'w') as f:
        json.dump(metrics, f, indent=2)

    assignments_path = output_dir / f"region_{args.region}_assignments.csv"
    pd.DataFrame([
        {'nucleus_id': nid, 'cell_type': ct}
        for nid, ct in result.assignments.items()
    ]).to_csv(assignments_path, index=False)

    logger.info("Saved results to %s", output_dir)


if __name__ == "__main__":
    main()
