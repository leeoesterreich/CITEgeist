#!/usr/bin/env python -u
"""SACE GEX grid search runner — one (combo, region) per SLURM task.

Loads fixed proportions from Stage 1 best config, runs SACE with
grid params, evaluates per-cell GEX against Xenium GT, writes result JSON.
"""

import sys as _sys
from pathlib import Path as _Path

_sys.path.insert(0, str(_Path(__file__).resolve().parents[3]))
import _paths as _cg  # noqa: E402

import argparse
import json
import logging
import os
import sys
import tempfile
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import pearsonr

REPO_ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO_ROOT))

from CITEgeist.model.gex.sace_gex import run_sace

BENCHMARK_ROOT = REPO_ROOT / "Benchmarking" / "xenium_benchmarking"
DATA_DIR = REPO_ROOT / "Benchmarking" / "xenium_pseudovisium" / "data_protein_gt"
GT_GEX_DIR = BENCHMARK_ROOT / "ground_truth_singler" / "ground_truth_7type" / "ground_truth_gex"
XENIUM_RAW = Path(str(_cg.root("DATA") / "10x_PublicData" / "Xenium_RNA_Proteomic_RenalCellCarcinoma"))

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
)
logger = logging.getLogger(__name__)


def load_region_data(region_id: int):
    """Load spot GEX, cell-to-spot mapping, cell type assignments."""
    gex_path = DATA_DIR / "h5ad_objects" / f"Xenium_region_{region_id}_GEX.h5ad"
    spot_adata = sc.read_h5ad(gex_path)

    mapping = pd.read_csv(DATA_DIR / "cell_to_spot_mapping.csv", index_col=0)
    labels = pd.read_csv(DATA_DIR / "cell_type_assignments.csv", index_col=0)

    region_cells = mapping[(mapping["region_id"] == region_id) & (mapping["spot_idx"] != -1)].copy()
    region_cells = region_cells.join(labels, how="inner")
    region_cells = region_cells[region_cells["cell_type"] != "Unknown"]

    cell_assignments = dict(zip(region_cells.index.astype(str), region_cells["cell_type"]))
    # Remap global spot_idx to region-local indices matching spot_adata
    barcode_to_local = {bc: i for i, bc in enumerate(spot_adata.obs_names)}
    cell_spot_map = region_cells[["spot_id"]].copy()
    cell_spot_map.columns = ["spot_barcode"]
    cell_spot_map["spot_idx"] = cell_spot_map["spot_barcode"].map(barcode_to_local)
    cell_spot_map["cell_id"] = cell_spot_map.index.astype(str)
    # Drop cells whose spot isn't in spot_adata
    cell_spot_map = cell_spot_map.dropna(subset=["spot_idx"])
    cell_spot_map["spot_idx"] = cell_spot_map["spot_idx"].astype(int)
    # Update cell_assignments to match
    cell_assignments = {cid: cell_assignments[cid] for cid in cell_spot_map["cell_id"].values if cid in cell_assignments}
    if "x" in region_cells.columns:
        cell_spot_map["x"] = region_cells.loc[cell_spot_map.index, "x"].values
        cell_spot_map["y"] = region_cells.loc[cell_spot_map.index, "y"].values
    else:
        # Use spot coordinates as cell coordinates (cells inherit their spot's position)
        spot_coords = spot_adata.obsm.get("spatial", np.zeros((spot_adata.n_obs, 2)))
        cell_spot_map["x"] = spot_coords[cell_spot_map["spot_idx"].values, 0]
        cell_spot_map["y"] = spot_coords[cell_spot_map["spot_idx"].values, 1]

    return spot_adata, cell_assignments, cell_spot_map, region_cells


def load_xenium_gt(region_cells):
    """Load Xenium single-cell GEX as ground truth."""
    h5_path = XENIUM_RAW / "cell_feature_matrix.h5"
    gt_adata = sc.read_10x_h5(h5_path)
    gt_adata.var_names_make_unique()
    cell_ids = region_cells.index.astype(str).tolist()
    common = sorted(set(cell_ids) & set(gt_adata.obs_names))
    return gt_adata[common].copy()


def evaluate_percell(pred_adata, gt_adata, common_genes):
    """Per-gene Pearson r between predicted and GT per-cell expression."""
    common_cells = sorted(set(pred_adata.obs_names) & set(gt_adata.obs_names))
    if len(common_cells) == 0:
        return {"percell_mean_gene_r": 0.0, "percell_median_gene_r": 0.0, "n_cells": 0}

    pred_X = pred_adata[common_cells, common_genes].X
    gt_X = gt_adata[common_cells, common_genes].X
    if hasattr(pred_X, "toarray"):
        pred_X = pred_X.toarray()
    if hasattr(gt_X, "toarray"):
        gt_X = gt_X.toarray()

    gene_rs = []
    for g in range(pred_X.shape[1]):
        p, g_col = pred_X[:, g], gt_X[:, g]
        if np.std(p) > 0 and np.std(g_col) > 0:
            r, _ = pearsonr(p, g_col)
            if not np.isnan(r):
                gene_rs.append(r)

    return {
        "percell_mean_gene_r": float(np.mean(gene_rs)) if gene_rs else 0.0,
        "percell_median_gene_r": float(np.median(gene_rs)) if gene_rs else 0.0,
        "n_genes": len(gene_rs),
        "n_cells": len(common_cells),
    }


def evaluate_spot_layers(pred_adata, region_id, cell_assignments):
    """Aggregate per-cell predictions to spot-type layers, compute flat_r."""
    gt_dir = GT_GEX_DIR / f"Xenium_region_{region_id}"
    type_names = sorted(set(cell_assignments.values()))

    gt_filename_map = {
        "B cells": "B_cells_GT.csv",
        "CD4+ T cells": "CD4pos_T_cells_GT.csv",
        "CD8+ T cells": "CD8pos_T_cells_GT.csv",
        "Macrophages": "Macrophages_GT.csv",
        "Endothelial": "Endothelial_GT.csv",
        "Epithelial": "Epithelial_GT.csv",
        "Fibroblasts": "Fibroblasts_GT.csv",
    }

    all_gt_vals = []
    all_pred_vals = []
    per_type_r = {}

    for ct in type_names:
        gt_file = gt_filename_map.get(ct)
        if gt_file is None or not (gt_dir / gt_file).exists():
            continue

        gt_layer = pd.read_csv(gt_dir / gt_file, index_col=0)
        ct_cells = [c for c, t in cell_assignments.items() if t == ct and c in pred_adata.obs_names]
        if len(ct_cells) == 0:
            continue

        pred_ct = pred_adata[ct_cells].X
        if hasattr(pred_ct, "toarray"):
            pred_ct = pred_ct.toarray()
        pred_layer = pd.DataFrame(pred_ct, index=ct_cells, columns=pred_adata.var_names)

        common_genes = sorted(set(gt_layer.columns) & set(pred_layer.columns))
        if len(common_genes) == 0:
            continue

        # Aggregate per-cell predictions to spot level (sum cells per spot)
        common_spots = sorted(set(gt_layer.index) & set(pred_adata.obs["spot_barcode"].values))
        if len(common_spots) == 0:
            continue
        spot_pred = pd.DataFrame(0.0, index=common_spots, columns=common_genes)
        for cell_id in ct_cells:
            spot_bc = pred_adata.obs.loc[cell_id, "spot_barcode"]
            if spot_bc in spot_pred.index:
                cell_expr = pred_layer.loc[cell_id, common_genes].values
                spot_pred.loc[spot_bc] += cell_expr

        gt_vals = np.log1p(gt_layer.loc[common_spots, common_genes].values).flatten()
        pred_vals = np.log1p(spot_pred.values).flatten()

        mask = (gt_vals > 0) | (pred_vals > 0)
        if mask.sum() > 2:
            r, _ = pearsonr(gt_vals[mask], pred_vals[mask])
            per_type_r[ct] = float(r) if not np.isnan(r) else 0.0
            all_gt_vals.append(gt_vals[mask])
            all_pred_vals.append(pred_vals[mask])

    overall_r = 0.0
    if all_gt_vals:
        gt_flat = np.concatenate(all_gt_vals)
        pred_flat = np.concatenate(all_pred_vals)
        r, _ = pearsonr(gt_flat, pred_flat)
        overall_r = float(r) if not np.isnan(r) else 0.0

    return {
        "spot_layer_overall_r": overall_r,
        "spot_layer_per_type_r": per_type_r,
    }


def atomic_write_json(data: dict, path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_path = tempfile.mkstemp(dir=path.parent, suffix=".tmp")
    try:
        with os.fdopen(fd, "w") as f:
            json.dump(data, f, indent=2)
        os.rename(tmp_path, str(path))
    except Exception:
        os.unlink(tmp_path)
        raise


def run_one(combo_id: int, region_id: int, params: dict, props_dir: Path, output_dir: Path):
    logger.info("combo=%d region=%d params=%s", combo_id, region_id, params)

    spot_adata, cell_assignments, cell_spot_map, region_cells = load_region_data(region_id)
    gt_adata = load_xenium_gt(region_cells)

    props_path = props_dir / f"region_{region_id}_predictions.csv"
    props = pd.read_csv(props_path, index_col=0)

    spot_counts = spot_adata.X
    if hasattr(spot_counts, "toarray"):
        spot_counts = spot_counts.toarray()
    spot_coords = spot_adata.obsm.get("spatial", np.zeros((spot_adata.n_obs, 2)))

    bw = params["bandwidth"]
    if bw == "None" or bw is None or (isinstance(bw, float) and np.isnan(bw)):
        bw = None
    else:
        bw = float(bw)

    t0 = time.time()
    _, cell_adata, sace_meta = run_sace(
        spot_counts=spot_counts,
        proportions=props,
        cell_assignments=cell_assignments,
        cell_spot_map=cell_spot_map,
        spot_coords=spot_coords,
        gene_names=list(spot_adata.var_names),
        max_iter=int(params["max_iter"]),
        n_0=float(params["n_0"]),
        bandwidth=bw,
        damping_eta=float(params["damping_eta"]),
        eps=float(params["eps"]),
    )
    runtime = time.time() - t0

    common_genes = sorted(set(cell_adata.var_names) & set(gt_adata.var_names))
    percell_metrics = evaluate_percell(cell_adata, gt_adata, common_genes)
    layer_metrics = evaluate_spot_layers(cell_adata, region_id, cell_assignments)

    result = {
        "combo_id": combo_id,
        "region_id": region_id,
        "params": {k: (v if v != "None" else None) for k, v in params.items()},
        "metrics": {**percell_metrics, **layer_metrics},
        "runtime_seconds": round(runtime, 1),
    }

    result_dir = output_dir / str(combo_id)
    atomic_write_json(result, result_dir / f"region_{region_id}.json")
    logger.info(
        "combo=%d region=%d flat_r=%.4f percell_r=%.4f runtime=%.1fs",
        combo_id,
        region_id,
        layer_metrics["spot_layer_overall_r"],
        percell_metrics["percell_mean_gene_r"],
        runtime,
    )


def main():
    parser = argparse.ArgumentParser(description="SACE GEX grid search runner")
    parser.add_argument("--combo-id", type=int, required=True)
    parser.add_argument("--region-id", type=int, required=True)
    parser.add_argument(
        "--grid-file",
        type=Path,
        default=Path(__file__).parent / "combos_sace.csv",
    )
    parser.add_argument("--props-dir", type=Path, required=True)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(__file__).parent / "results" / "sace",
    )
    args = parser.parse_args()

    grid_df = pd.read_csv(args.grid_file)
    n_combos = len(grid_df)
    if not (0 <= args.combo_id < n_combos):
        logger.error("combo_id %d out of range [0, %d)", args.combo_id, n_combos)
        sys.exit(1)
    if not (0 <= args.region_id <= 4):
        logger.error("region_id %d out of range [0, 4]", args.region_id)
        sys.exit(1)

    row = grid_df.iloc[args.combo_id]
    params = {col: row[col] for col in grid_df.columns if col != "combo_id"}

    run_one(args.combo_id, args.region_id, params, args.props_dir, args.output_dir)


if __name__ == "__main__":
    main()
