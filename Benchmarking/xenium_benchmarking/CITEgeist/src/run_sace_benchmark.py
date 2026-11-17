#!/usr/bin/env python
"""SACE per-cell GEX benchmark on Xenium pseudo-Visium.

Runs the Spatially-Adaptive Compositional EM (SACE) algorithm on one
Xenium pseudo-Visium region and evaluates per-cell GEX predictions against
the Xenium single-cell ground truth.

Pipeline:
1. Load spot-level GEX h5ad for the region
2. Load cell-to-spot mapping + GT cell type assignments, filter to region
3. Build proportions from GT cell counts (non-Unknown only, renormalized)
4. Run SACE via CITEgeist.model.sace_gex.run_sace
5. Load Xenium single-cell GEX as ground truth
6. Evaluate per-gene Pearson r (overall and per-type)
7. Save results JSON

Usage:
    python run_sace_benchmark.py --region 0 --output_dir results_sace/
"""

import sys as _sys
from pathlib import Path as _Path

_sys.path.insert(0, str(_Path(__file__).resolve().parents[3]))
import _paths as _cg  # noqa: E402

import argparse
import json
import logging
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import pearsonr

# ---------------------------------------------------------------------------
# Paths (absolute, no dirname tricks)
# ---------------------------------------------------------------------------
REPO_ROOT = Path(str(_cg.REPO_ROOT))
sys.path.insert(0, str(REPO_ROOT))

DATA_DIR = REPO_ROOT / "Benchmarking" / "xenium_pseudovisium" / "data_protein_gt"
XENIUM_RAW = Path(str(_cg.root("DATA") / "10x_PublicData" / "Xenium_RNA_Proteomic_RenalCellCarcinoma"))

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
)
logger = logging.getLogger(__name__)


def load_region_data(region_id: int):
    """Load all data needed for SACE benchmark on one region.

    Returns:
        spot_adata: AnnData (n_spots, n_genes) with raw counts
        cell_assignments: Dict[cell_id_str -> type_name]
        cell_spot_map: DataFrame with columns cell_id, spot_barcode, spot_idx, x, y
        props: DataFrame (n_spots, n_types) renormalized to sum=1
        region_cells: DataFrame of cell metadata for this region
    """
    # 1. Spot-level GEX
    gex_path = DATA_DIR / "h5ad_objects" / f"Xenium_region_{region_id}_GEX.h5ad"
    spot_adata = sc.read_h5ad(gex_path)
    logger.info("Loaded spot GEX: %s", spot_adata.shape)

    # Build spot_id -> h5ad row index mapping
    spot_name_to_row = {name: idx for idx, name in enumerate(spot_adata.obs_names)}

    # 2. Cell-to-spot mapping + GT labels
    mapping = pd.read_csv(DATA_DIR / "cell_to_spot_mapping.csv", index_col=0)
    labels = pd.read_csv(DATA_DIR / "cell_type_assignments.csv", index_col=0)

    # Filter to this region, assigned spots, non-Unknown
    region_cells = mapping[(mapping["region_id"] == region_id) & (mapping["spot_idx"] != -1)].copy()
    region_cells = region_cells.join(labels, how="inner")
    region_cells = region_cells[region_cells["cell_type"] != "Unknown"]
    logger.info(
        "Region %d: %d non-Unknown cells in %d spots",
        region_id,
        len(region_cells),
        region_cells["spot_id"].nunique(),
    )

    # 3. Map spot_id -> h5ad row index (NOT the numeric spot_idx from the CSV)
    region_cells["h5ad_row"] = region_cells["spot_id"].map(spot_name_to_row)
    # Drop cells whose spot_id is not in the h5ad (should not happen)
    n_before = len(region_cells)
    region_cells = region_cells.dropna(subset=["h5ad_row"])
    region_cells["h5ad_row"] = region_cells["h5ad_row"].astype(int)
    if len(region_cells) < n_before:
        logger.warning("Dropped %d cells with unmapped spot_ids", n_before - len(region_cells))

    # 4. Build cell_assignments dict
    cell_assignments = dict(zip(region_cells.index.astype(str), region_cells["cell_type"]))

    # 5. Build cell_spot_map DataFrame
    #    spot_idx here must be the h5ad ROW index, not the CSV spot_idx
    cell_spot_map = pd.DataFrame(
        {
            "cell_id": region_cells.index.astype(str).values,
            "spot_barcode": region_cells["spot_id"].values,
            "spot_idx": region_cells["h5ad_row"].values,
            "x": region_cells["spot_x"].values.astype(float),
            "y": region_cells["spot_y"].values.astype(float),
        }
    )

    # 6. Build GT proportions from non-Unknown cell counts, renormalized to sum=1
    spot_names = list(spot_adata.obs_names)
    type_names = sorted(set(cell_assignments.values()))
    props = pd.DataFrame(0.0, index=spot_names, columns=type_names)
    for _, row in region_cells.iterrows():
        sid = row["spot_id"]
        ct = row["cell_type"]
        if sid in props.index and ct in props.columns:
            props.loc[sid, ct] += 1
    row_sums = props.sum(axis=1)
    row_sums = row_sums.replace(0, 1)
    props = props.div(row_sums, axis=0)

    return spot_adata, cell_assignments, cell_spot_map, props, region_cells


def load_xenium_singlecell_gex(cell_ids):
    """Load Xenium single-cell GEX for given cell IDs.

    Args:
        cell_ids: Iterable of cell ID strings to load.

    Returns:
        AnnData with raw single-cell counts for matched cells.
    """
    logger.info("Loading Xenium single-cell GEX from %s ...", XENIUM_RAW)
    adata = sc.read_10x_h5(XENIUM_RAW / "cell_feature_matrix.h5", gex_only=True)
    adata.var_names_make_unique()
    cell_id_set = set(cell_ids)
    mask = [name in cell_id_set for name in adata.obs_names]
    adata = adata[mask].copy()
    logger.info(
        "Loaded %d / %d requested single cells from Xenium",
        adata.n_obs,
        len(cell_id_set),
    )
    return adata


# ---------------------------------------------------------------------------
# Baseline methods
# ---------------------------------------------------------------------------


def baseline_equal_split(Y, cell_assignments, cell_spot_map, gene_names, type_names, props):
    """Split spot counts to cells proportionally by type, equal within type.

    For each cell of type t in spot s:
        x_c = (pi_{s,t} * y_s) / n_{s,t}

    Args:
        Y: (N, G) spot-level count matrix (dense).
        cell_assignments: Dict[cell_id -> type_name].
        cell_spot_map: DataFrame with cell_id, spot_barcode, spot_idx, x, y.
        gene_names: List of gene names.
        type_names: List of cell type names.
        props: DataFrame (N, T) proportions.

    Returns:
        AnnData (n_cells, G) with obs columns: cell_id, spot_barcode, cell_type.
    """
    import scanpy as sc

    N, G = Y.shape
    type_to_idx = {t: i for i, t in enumerate(type_names)}
    props_arr = props.values  # (N, T)

    # Count cells per (spot, type)
    n_st = np.zeros((N, len(type_names)), dtype=np.float64)
    for _, row in cell_spot_map.iterrows():
        cid = row["cell_id"]
        sidx = int(row["spot_idx"])
        ct = cell_assignments.get(cid)
        if ct is not None and ct in type_to_idx:
            n_st[sidx, type_to_idx[ct]] += 1

    # Allocate per cell
    cell_ids = []
    cell_types = []
    cell_spots = []
    cell_xs = []
    cell_ys = []
    cell_gex = np.zeros((len(cell_spot_map), G), dtype=np.float64)

    for i, (_, row) in enumerate(cell_spot_map.iterrows()):
        cid = row["cell_id"]
        sidx = int(row["spot_idx"])
        ct = cell_assignments.get(cid)
        if ct is None or ct not in type_to_idx:
            continue
        tidx = type_to_idx[ct]
        count = n_st[sidx, tidx]
        if count > 0:
            cell_gex[i, :] = (props_arr[sidx, tidx] * Y[sidx, :]) / count
        cell_ids.append(cid)
        cell_types.append(ct)
        cell_spots.append(row["spot_barcode"])
        cell_xs.append(row["x"])
        cell_ys.append(row["y"])

    cell_gex = cell_gex[: len(cell_ids), :]
    adata = sc.AnnData(
        X=cell_gex,
        obs=pd.DataFrame(
            {
                "cell_id": cell_ids,
                "spot_barcode": cell_spots,
                "cell_type": cell_types,
            }
        ),
        var=pd.DataFrame(index=gene_names),
    )
    adata.obs_names = cell_ids
    adata.obsm["spatial"] = np.column_stack([cell_xs, cell_ys])
    return adata


def baseline_type_reference(Y, cell_assignments, cell_spot_map, gene_names, type_names, props):
    """Solve for type-level reference profiles via least squares, then allocate.

    Solves Y = P @ beta via np.linalg.lstsq to get deconfounded type profiles.
    Allocates per (spot, type, gene) proportionally, then divides among cells.

    Args:
        Y: (N, G) spot-level count matrix (dense).
        cell_assignments: Dict[cell_id -> type_name].
        cell_spot_map: DataFrame with cell_id, spot_barcode, spot_idx, x, y.
        gene_names: List of gene names.
        type_names: List of cell type names.
        props: DataFrame (N, T) proportions.

    Returns:
        AnnData (n_cells, G) with obs columns: cell_id, spot_barcode, cell_type.
    """
    import scanpy as sc

    N, G = Y.shape
    T = len(type_names)
    type_to_idx = {t: i for i, t in enumerate(type_names)}
    P = props.values  # (N, T)

    # Solve Y = P @ beta for beta (T, G)
    beta, _, _, _ = np.linalg.lstsq(P, Y, rcond=None)
    beta = np.maximum(beta, 0.0)  # (T, G)

    # Spot-type expected expression: E[s,t,:] = props[s,t] * beta[t,:]
    # Normalize per spot so sum_t E[s,t,:] = Y[s,:], then multiply
    # weights[s,t,g] = props[s,t] * beta[t,g]
    # E[s,t,:] = weights[s,t,:] / sum_t(weights[s,t,:]) * Y[s,:]

    # Count cells per (spot, type)
    n_st = np.zeros((N, T), dtype=np.float64)
    for _, row in cell_spot_map.iterrows():
        cid = row["cell_id"]
        sidx = int(row["spot_idx"])
        ct = cell_assignments.get(cid)
        if ct is not None and ct in type_to_idx:
            n_st[sidx, type_to_idx[ct]] += 1

    # Precompute spot-type GEX allocation
    # spot_type_gex[s, t, :] = E[s,t,:] / n_st[s,t]
    # where E[s,t,:] = (props[s,t]*beta[t,:]) / sum_t'(props[s,t']*beta[t',:]) * Y[s,:]
    weights = P[:, :, None] * beta[None, :, :]  # (N, T, G)
    weight_sums = weights.sum(axis=1, keepdims=True)  # (N, 1, G)
    weight_sums = np.maximum(weight_sums, 1e-12)
    E = weights / weight_sums * Y[:, None, :]  # (N, T, G)

    cell_ids = []
    cell_types = []
    cell_spots = []
    cell_xs = []
    cell_ys = []
    cell_gex = np.zeros((len(cell_spot_map), G), dtype=np.float64)

    for i, (_, row) in enumerate(cell_spot_map.iterrows()):
        cid = row["cell_id"]
        sidx = int(row["spot_idx"])
        ct = cell_assignments.get(cid)
        if ct is None or ct not in type_to_idx:
            continue
        tidx = type_to_idx[ct]
        count = n_st[sidx, tidx]
        if count > 0:
            cell_gex[i, :] = E[sidx, tidx, :] / count
        cell_ids.append(cid)
        cell_types.append(ct)
        cell_spots.append(row["spot_barcode"])
        cell_xs.append(row["x"])
        cell_ys.append(row["y"])

    cell_gex = cell_gex[: len(cell_ids), :]
    adata = sc.AnnData(
        X=cell_gex,
        obs=pd.DataFrame(
            {
                "cell_id": cell_ids,
                "spot_barcode": cell_spots,
                "cell_type": cell_types,
            }
        ),
        var=pd.DataFrame(index=gene_names),
    )
    adata.obs_names = cell_ids
    adata.obsm["spatial"] = np.column_stack([cell_xs, cell_ys])
    return adata


def run_sace_iter0(Y, cell_assignments, cell_spot_map, gene_names, type_names, props, spot_coords):
    """Run SACE with max_iter=1 (single E-step, no M-step refinement).

    Args:
        Y: (N, G) spot-level count matrix (dense).
        cell_assignments: Dict[cell_id -> type_name].
        cell_spot_map: DataFrame with cell_id, spot_barcode, spot_idx, x, y.
        gene_names: List of gene names.
        type_names: List of cell type names.
        props: DataFrame (N, T) proportions.
        spot_coords: (N, 2) spatial coordinates.

    Returns:
        AnnData (n_cells, G) with obs columns: cell_id, spot_barcode, cell_type.
    """
    from CITEgeist.model.gex.sace_gex import run_sace

    _, cell_adata, _ = run_sace(
        spot_counts=Y,
        proportions=props,
        cell_assignments=cell_assignments,
        cell_spot_map=cell_spot_map,
        spot_coords=spot_coords,
        gene_names=gene_names,
        max_iter=1,
    )
    return cell_adata


def evaluate_per_cell(pred_adata, gt_adata, common_genes):
    """Per-gene Pearson r between predicted and GT per-cell expression.

    Args:
        pred_adata: AnnData from SACE (per-cell allocated counts).
        gt_adata: AnnData from Xenium raw single-cell.
        common_genes: List of gene names to evaluate.

    Returns:
        Dict with evaluation metrics.
    """
    common_cells = sorted(set(pred_adata.obs_names) & set(gt_adata.obs_names))
    if len(common_cells) == 0:
        logger.warning("No common cells found between prediction and GT!")
        return {"mean_gene_r": 0.0, "n_cells": 0, "n_genes": 0}

    pred_X = pred_adata[common_cells, common_genes].X
    gt_X = gt_adata[common_cells, common_genes].X
    if hasattr(pred_X, "toarray"):
        pred_X = pred_X.toarray()
    if hasattr(gt_X, "toarray"):
        gt_X = gt_X.toarray()

    # Overall per-gene r
    gene_rs = []
    gene_names_valid = []
    for g_idx, gene in enumerate(common_genes):
        p, g = pred_X[:, g_idx], gt_X[:, g_idx]
        if p.std() < 1e-10 or g.std() < 1e-10:
            continue
        r, _ = pearsonr(p, g)
        if np.isfinite(r):
            gene_rs.append(float(r))
            gene_names_valid.append(gene)

    # Per-type per-gene r
    type_rs = {}
    if "cell_type" in pred_adata.obs.columns:
        for ct in sorted(pred_adata.obs["cell_type"].unique()):
            ct_cells = [c for c in common_cells if pred_adata.obs.loc[c, "cell_type"] == ct]
            if len(ct_cells) < 10:
                continue
            ct_pred = pred_adata[ct_cells, common_genes].X
            ct_gt = gt_adata[ct_cells, common_genes].X
            if hasattr(ct_pred, "toarray"):
                ct_pred = ct_pred.toarray()
            if hasattr(ct_gt, "toarray"):
                ct_gt = ct_gt.toarray()
            ct_gene_rs = []
            for g_idx in range(len(common_genes)):
                p, g = ct_pred[:, g_idx], ct_gt[:, g_idx]
                if p.std() < 1e-10 or g.std() < 1e-10:
                    continue
                r, _ = pearsonr(p, g)
                if np.isfinite(r):
                    ct_gene_rs.append(float(r))
            if ct_gene_rs:
                type_rs[ct] = float(np.mean(ct_gene_rs))

    # Stratified by gene count level (high vs low, median split)
    gt_gene_totals = np.asarray(gt_X.sum(axis=0)).ravel()
    nonzero_totals = gt_gene_totals[gt_gene_totals > 0]
    median_total = float(np.median(nonzero_totals)) if len(nonzero_totals) > 0 else 0.0
    high_rs, low_rs = [], []
    for g_idx, gene in enumerate(common_genes):
        p, g = pred_X[:, g_idx], gt_X[:, g_idx]
        if p.std() < 1e-10 or g.std() < 1e-10:
            continue
        r, _ = pearsonr(p, g)
        if not np.isfinite(r):
            continue
        if gt_gene_totals[g_idx] >= median_total and gt_gene_totals[g_idx] > 0:
            high_rs.append(float(r))
        elif gt_gene_totals[g_idx] > 0:
            low_rs.append(float(r))

    return {
        "mean_gene_r": float(np.mean(gene_rs)) if gene_rs else 0.0,
        "median_gene_r": float(np.median(gene_rs)) if gene_rs else 0.0,
        "n_genes_evaluated": len(gene_rs),
        "n_cells_evaluated": len(common_cells),
        "per_type_mean_r": type_rs,
        "per_gene_r": dict(zip(gene_names_valid, gene_rs)),
        "high_count_genes_mean_r": float(np.mean(high_rs)) if high_rs else 0.0,
        "low_count_genes_mean_r": float(np.mean(low_rs)) if low_rs else 0.0,
        "n_high_count_genes": len(high_rs),
        "n_low_count_genes": len(low_rs),
    }


def evaluate_spot_type_layers(pred_adata, region_id):
    """Legacy spot-type layer evaluation (aggregated Pearson r).

    Aggregates per-cell predictions back to spot-type layers, then compares
    against ground truth CSV files using log1p + nonzero-masked Pearson r.

    Args:
        pred_adata: AnnData (n_cells, G) with obs columns cell_type, spot_barcode.
        region_id: Region ID for loading GT files.

    Returns:
        Dict with overall_r and per_type_r.
    """
    gt_dir = DATA_DIR / "ground_truth_gex" / f"Xenium_region_{region_id}"

    gt_filename_map = {
        "B cells": "B_cells_GT.csv",
        "T cells": None,  # merged from CD4 + CD8
        "Macrophages": "Macrophages_GT.csv",
        "Endothelial": "Endothelial_GT.csv",
        "Epithelial": "Epithelial_GT.csv",
        "Fibroblasts": "Fibroblasts_GT.csv",
    }

    # Aggregate per-cell predictions to spot-type layers
    pred_X = pred_adata.X
    if hasattr(pred_X, "toarray"):
        pred_X = pred_X.toarray()
    pred_genes = list(pred_adata.var_names)

    spot_type_sums = {}  # (spot, type) -> (G,) array
    for i in range(pred_adata.n_obs):
        spot = pred_adata.obs.iloc[i]["spot_barcode"]
        ct = pred_adata.obs.iloc[i]["cell_type"]
        # Merge T cell subtypes
        if "CD4" in ct or "CD8" in ct or ct in ("CD4pos_T_cells", "CD8pos_T_cells"):
            ct = "T cells"
        key = (spot, ct)
        if key not in spot_type_sums:
            spot_type_sums[key] = np.zeros(len(pred_genes), dtype=np.float64)
        spot_type_sums[key] += pred_X[i, :]

    type_rs = {}
    all_pv_parts = []
    all_gv_parts = []

    for type_name, gt_file in gt_filename_map.items():
        # Load GT
        if type_name == "T cells":
            # Sum CD4 + CD8 GT files
            cd4_path = gt_dir / "CD4pos_T_cells_GT.csv"
            cd8_path = gt_dir / "CD8pos_T_cells_GT.csv"
            if not cd4_path.exists() or not cd8_path.exists():
                logger.warning("Missing T cell GT for region %d", region_id)
                continue
            gt_cd4 = pd.read_csv(cd4_path, index_col=0).T  # spots x genes
            gt_cd8 = pd.read_csv(cd8_path, index_col=0).T
            # Align and sum
            common_spots = sorted(set(gt_cd4.index) & set(gt_cd8.index))
            common_gt_genes = sorted(set(gt_cd4.columns) & set(gt_cd8.columns))
            gt_df = gt_cd4.loc[common_spots, common_gt_genes] + gt_cd8.loc[common_spots, common_gt_genes]
        else:
            gt_path = gt_dir / gt_file
            if not gt_path.exists():
                logger.warning("Missing GT file %s", gt_path)
                continue
            gt_df = pd.read_csv(gt_path, index_col=0).T  # spots x genes

        # Build pred DataFrame for this type
        spots_in_pred = sorted(set(s for (s, t) in spot_type_sums if t == type_name))
        if not spots_in_pred:
            logger.warning("No predicted cells for type %s", type_name)
            continue

        pred_df = pd.DataFrame(
            np.stack([spot_type_sums[(s, type_name)] for s in spots_in_pred]),
            index=spots_in_pred,
            columns=pred_genes,
        )

        # Align spots and genes
        common_spots = sorted(set(pred_df.index) & set(gt_df.index))
        common_genes = sorted(set(pred_df.columns) & set(gt_df.columns))
        if len(common_spots) < 5 or len(common_genes) < 5:
            logger.warning(
                "Insufficient overlap for %s: %d spots, %d genes",
                type_name,
                len(common_spots),
                len(common_genes),
            )
            continue

        pred_vals = np.log1p(pred_df.loc[common_spots, common_genes].values)
        gt_vals = np.log1p(gt_df.loc[common_spots, common_genes].values)

        # Per-type: flatten, mask nonzero, compute r
        pv = pred_vals.ravel()
        gv = gt_vals.ravel()
        mask = (gv > 0) | (pv > 0)
        if mask.sum() < 10:
            continue
        r, _ = pearsonr(pv[mask], gv[mask])
        if np.isfinite(r):
            type_rs[type_name] = float(r)
            all_pv_parts.append(pv[mask])
            all_gv_parts.append(gv[mask])

    if all_pv_parts:
        pooled_pv = np.concatenate(all_pv_parts)
        pooled_gv = np.concatenate(all_gv_parts)
        overall_r = float(pearsonr(pooled_pv, pooled_gv)[0])
    else:
        overall_r = 0.0
    return {
        "overall_r": overall_r,
        "per_type_r": type_rs,
    }


def _collect_method_results(method_name, pred_adata, gt_adata, common_genes, region_id, include_per_gene=False):
    """Run both evaluation metrics on a method's predictions.

    Args:
        method_name: Name for logging.
        pred_adata: Predicted per-cell AnnData.
        gt_adata: Ground truth per-cell AnnData.
        common_genes: List of common gene names.
        region_id: Region ID for spot-type layer GT.
        include_per_gene: Whether to include per_gene_r dict.

    Returns:
        Dict with percell and spot_layer metrics.
    """
    percell = evaluate_per_cell(pred_adata, gt_adata, common_genes)
    spot_layer = evaluate_spot_type_layers(pred_adata, region_id)

    result = {
        "percell_mean_gene_r": percell["mean_gene_r"],
        "percell_median_gene_r": percell["median_gene_r"],
        "percell_n_genes": percell["n_genes_evaluated"],
        "percell_n_cells": percell["n_cells_evaluated"],
        "percell_per_type_r": percell["per_type_mean_r"],
        "spot_layer_overall_r": spot_layer["overall_r"],
        "spot_layer_per_type_r": spot_layer["per_type_r"],
    }
    if include_per_gene:
        result["per_gene_r"] = percell.get("per_gene_r", {})

    logger.info(
        "  %s: percell mean_r=%.4f, spot_layer r=%.4f",
        method_name,
        percell["mean_gene_r"],
        spot_layer["overall_r"],
    )
    return result


def main():
    parser = argparse.ArgumentParser(description="SACE per-cell GEX benchmark on Xenium pseudo-Visium")
    parser.add_argument("--region", type=int, required=True, help="Region ID (0-4)")
    parser.add_argument("--output_dir", type=str, default="results_sace", help="Output directory for results JSON")
    parser.add_argument("--max_iter", type=int, default=10, help="Max SACE EM iterations")
    parser.add_argument("--n_0", type=float, default=10.0, help="Shrinkage prior strength")
    parser.add_argument("--bandwidth", type=float, default=None, help="Kernel bandwidth (None=auto)")
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # --- Load data ---
    spot_adata, cell_assignments, cell_spot_map, props, region_cells = load_region_data(args.region)

    spot_X = spot_adata.X
    if hasattr(spot_X, "toarray"):
        spot_X = spot_X.toarray()
    spot_X = spot_X.astype(np.float64)

    gene_names = list(spot_adata.var_names)
    type_names = sorted(props.columns)
    spot_coords = np.column_stack([spot_adata.obs["x"].values, spot_adata.obs["y"].values])

    # --- Load GT single-cell GEX (once) ---
    gt_adata = load_xenium_singlecell_gex(list(cell_assignments.keys()))

    # --- Run all methods ---
    logger.info("=" * 60)
    logger.info("Running 4 methods on region %d", args.region)
    logger.info("=" * 60)

    # 1. Equal split baseline
    logger.info("--- Method 1/4: equal_split ---")
    t0 = time.time()
    equal_adata = baseline_equal_split(
        spot_X,
        cell_assignments,
        cell_spot_map,
        gene_names,
        type_names,
        props,
    )
    t_equal = time.time() - t0
    logger.info("equal_split: %d cells in %.1f s", equal_adata.n_obs, t_equal)

    # 2. Type reference baseline
    logger.info("--- Method 2/4: type_reference ---")
    t0 = time.time()
    typeref_adata = baseline_type_reference(
        spot_X,
        cell_assignments,
        cell_spot_map,
        gene_names,
        type_names,
        props,
    )
    t_typeref = time.time() - t0
    logger.info("type_reference: %d cells in %.1f s", typeref_adata.n_obs, t_typeref)

    # 3. SACE iter0 (single E-step)
    logger.info("--- Method 3/4: sace_iter0 ---")
    t0 = time.time()
    iter0_adata = run_sace_iter0(
        spot_X,
        cell_assignments,
        cell_spot_map,
        gene_names,
        type_names,
        props,
        spot_coords,
    )
    t_iter0 = time.time() - t0
    logger.info("sace_iter0: %d cells in %.1f s", iter0_adata.n_obs, t_iter0)

    # 4. SACE full
    logger.info("--- Method 4/4: sace_full ---")
    from CITEgeist.model.gex.sace_gex import run_sace

    t0 = time.time()
    spot_type_gex, sace_adata, diagnostics = run_sace(
        spot_counts=spot_X,
        proportions=props,
        cell_assignments=cell_assignments,
        cell_spot_map=cell_spot_map,
        spot_coords=spot_coords,
        gene_names=gene_names,
        max_iter=args.max_iter,
        n_0=args.n_0,
        bandwidth=args.bandwidth,
    )
    t_sace = time.time() - t0
    logger.info("sace_full: %d cells in %.1f s", sace_adata.n_obs, t_sace)

    # --- Evaluate all methods ---
    logger.info("=" * 60)
    logger.info("Evaluating all methods")
    logger.info("=" * 60)

    # Find common genes across all predictions and GT
    all_pred_genes = (
        set(equal_adata.var_names)
        & set(typeref_adata.var_names)
        & set(iter0_adata.var_names)
        & set(sace_adata.var_names)
        & set(gt_adata.var_names)
    )
    common_genes = sorted(all_pred_genes)
    logger.info("Common genes across all methods + GT: %d", len(common_genes))

    equal_results = _collect_method_results(
        "equal_split",
        equal_adata,
        gt_adata,
        common_genes,
        args.region,
    )
    typeref_results = _collect_method_results(
        "type_reference",
        typeref_adata,
        gt_adata,
        common_genes,
        args.region,
    )
    iter0_results = _collect_method_results(
        "sace_iter0",
        iter0_adata,
        gt_adata,
        common_genes,
        args.region,
    )
    sace_results = _collect_method_results(
        "sace_full",
        sace_adata,
        gt_adata,
        common_genes,
        args.region,
        include_per_gene=True,
    )

    # --- Build results JSON ---
    results = {
        "region": args.region,
        "n_common_genes": len(common_genes),
        "sace_params": {
            "max_iter": args.max_iter,
            "n_0": args.n_0,
            "bandwidth": args.bandwidth,
        },
        "sace_diagnostics": {
            "n_iterations": diagnostics.get("n_iterations"),
            "damping_activated": diagnostics.get("damping_activated"),
            "final_bandwidth": diagnostics.get("final_bandwidth"),
            "log_likelihood_trace": diagnostics.get("log_likelihood_trace"),
        },
        "sace_runtime_seconds": round(t_sace, 1),
        "equal_split": equal_results,
        "type_reference": typeref_results,
        "sace_iter0": iter0_results,
        "sace_full": sace_results,
    }

    # --- Save ---
    out_file = output_dir / f"sace_region_{args.region}.json"
    with open(out_file, "w") as f:
        json.dump(results, f, indent=2, default=str)

    # --- Summary table ---
    logger.info("=" * 60)
    logger.info(
        "SUMMARY - Region %d (%d genes, %d cells)", args.region, len(common_genes), sace_results["percell_n_cells"]
    )
    logger.info("%-16s %10s %10s %10s", "Method", "PC mean_r", "PC med_r", "SL r")
    logger.info("-" * 50)
    for name, res in [
        ("equal_split", equal_results),
        ("type_reference", typeref_results),
        ("sace_iter0", iter0_results),
        ("sace_full", sace_results),
    ]:
        logger.info(
            "%-16s %10.4f %10.4f %10.4f",
            name,
            res["percell_mean_gene_r"],
            res["percell_median_gene_r"],
            res["spot_layer_overall_r"],
        )
    logger.info("=" * 60)
    logger.info("Results saved to %s", out_file)


if __name__ == "__main__":
    main()
