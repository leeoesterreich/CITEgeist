#!/usr/bin/env python
"""Unified GEX benchmark: SACE vs QP vs NB on simulated data.

Runs 5 GEX deconvolution methods on one simulated replicate and evaluates
all through identical spot_layer_r pipeline. No per-cell GT exists in
simulation, so only spot-type level evaluation is performed.

Methods:
  1. equal_split      — divide spot counts equally within type (baseline)
  2. type_reference    — lstsq type profiles then equal split (baseline)
  3. sace_iter0        — SACE with max_iter=1 (production default)
  4. qp_regression     — load pre-existing QP layer CSVs
  5. nb_emission       — load pre-existing NB layer CSVs

Usage:
    python run_sace_simulation_benchmark.py --rep 0 --output_dir results_unified_gex/
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
# Paths
# ---------------------------------------------------------------------------
REPO_ROOT = Path(str(_cg.REPO_ROOT))
sys.path.insert(0, str(REPO_ROOT))
BENCH_DIR = REPO_ROOT / "Benchmarking" / "simulation_benchmarking" / "CITEgeist"
SIM_BASE = Path(
    str(_cg.root("DATA") / "Wu_Visium" / "Simulations" / "scCube_12k" / "replicates" / "high_seg")
)

# Pre-existing layer directories
QP_LAYER_BASE = BENCH_DIR / "results" / "high_seg"
NB_LAYER_BASE = BENCH_DIR / "results_nb" / "high_seg" / "gt_ablation"

# 9 cell types in simulation
SIM_TYPES = [
    "B-cells",
    "CAFs",
    "Cancer Epithelial",
    "Endothelial",
    "Myeloid",
    "Normal Epithelial",
    "PVL",
    "Plasmablasts",
    "T-cells",
]

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------


def load_sim_region_data(rep_id: int):
    """Load simulation data for one replicate.

    Args:
        rep_id: Replicate ID (0-4).

    Returns:
        spot_adata: AnnData (n_spots, n_genes) raw counts.
        cell_assignments: Dict[cell_id_str -> type_name].
        cell_spot_map: DataFrame with cell_id, spot_barcode, spot_idx, x, y.
        props: DataFrame (n_spots, n_types) GT proportions.
        spot_coords: (n_spots, 2) array.
    """
    # Spot-level GEX
    gex_path = SIM_BASE / "h5ad_objects" / f"Wu_rep_{rep_id}_GEX.h5ad"
    spot_adata = sc.read_h5ad(gex_path)
    logger.info("Loaded spot GEX: %s", spot_adata.shape)

    spot_name_to_row = {name: idx for idx, name in enumerate(spot_adata.obs_names)}

    # Cell index (true cell-to-type assignments)
    index_path = SIM_BASE / "ST_sim" / f"Wu_ST_{rep_id}_index.csv"
    index_df = pd.read_csv(index_path)
    logger.info("Loaded cell index: %d cells", len(index_df))

    # Build cell_assignments
    cell_assignments = dict(zip(index_df["Cell"].astype(str), index_df["Cell_type"]))

    # Map spot name to h5ad row index
    index_df["spot_idx"] = index_df["spot"].map(spot_name_to_row)
    # Drop cells whose spot is not in h5ad
    index_df = index_df.dropna(subset=["spot_idx"])
    index_df["spot_idx"] = index_df["spot_idx"].astype(int)

    cell_spot_map = pd.DataFrame(
        {
            "cell_id": index_df["Cell"].astype(str).values,
            "spot_barcode": index_df["spot"].values,
            "spot_idx": index_df["spot_idx"].values,
            "x": index_df["point_x"].values.astype(float),
            "y": index_df["point_y"].values.astype(float),
        }
    )

    # Update cell_assignments to only include mapped cells
    mapped_cells = set(cell_spot_map["cell_id"].values)
    cell_assignments = {k: v for k, v in cell_assignments.items() if k in mapped_cells}

    # GT proportions
    prop_path = SIM_BASE / "ST_sim" / f"Wu_ST_{rep_id}_prop.csv"
    prop_df = pd.read_csv(prop_path)
    # Build proportions DataFrame aligned to spot_adata
    type_cols = [c for c in prop_df.columns if c not in ("spot", "spot_x", "spot_y")]
    props = pd.DataFrame(0.0, index=spot_adata.obs_names, columns=type_cols)
    for _, row in prop_df.iterrows():
        sid = row["spot"]
        if sid in props.index:
            for tc in type_cols:
                props.loc[sid, tc] = row[tc]
    # Renormalize
    row_sums = props.sum(axis=1).replace(0, 1)
    props = props.div(row_sums, axis=0)

    spot_coords = np.column_stack(
        [
            spot_adata.obs["spot_x"].values.astype(float),
            spot_adata.obs["spot_y"].values.astype(float),
        ]
    )

    return spot_adata, cell_assignments, cell_spot_map, props, spot_coords


# ---------------------------------------------------------------------------
# Evaluation
# ---------------------------------------------------------------------------


def evaluate_sim_spot_layers(pred_layers, rep_id, type_names):
    """Evaluate spot-type layer predictions against simulation GT.

    Args:
        pred_layers: Dict[type_name -> DataFrame (spots x genes)] predicted layers.
        rep_id: Replicate ID for loading GT files.
        type_names: List of canonical type names.

    Returns:
        Dict with overall_r and per_type_r.
    """
    gt_dir = SIM_BASE / "ST_GEX_sim" / f"sample_{rep_id}" / "layers"

    type_rs = {}
    all_pv, all_gv = [], []

    for type_name in type_names:
        # GT file uses spaces in name
        gt_file = gt_dir / f"{type_name}_GT.csv"
        if not gt_file.exists():
            logger.warning("Missing GT: %s", gt_file)
            continue

        gt_df = pd.read_csv(gt_file, index_col=0)  # already spots x genes

        if type_name not in pred_layers:
            logger.warning("No prediction for type %s", type_name)
            continue
        pred_df = pred_layers[type_name]

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

        pv = pred_vals.ravel()
        gv = gt_vals.ravel()
        mask = (gv > 0) | (pv > 0)
        if mask.sum() < 10:
            continue
        r, _ = pearsonr(pv[mask], gv[mask])
        if np.isfinite(r):
            type_rs[type_name] = float(r)
            all_pv.append(pv[mask])
            all_gv.append(gv[mask])

    if all_pv:
        overall_r = float(pearsonr(np.concatenate(all_pv), np.concatenate(all_gv))[0])
    else:
        overall_r = 0.0

    return {"overall_r": overall_r, "per_type_r": type_rs}


def percell_to_spot_layers(pred_adata, spot_names, type_names):
    """Aggregate per-cell AnnData back to spot-type layers.

    Args:
        pred_adata: AnnData (n_cells, G) with obs cell_type, spot_barcode.
        spot_names: List of spot names (index for output DataFrames).
        type_names: List of canonical type names.

    Returns:
        Dict[type_name -> DataFrame (spots x genes)].
    """
    pred_X = pred_adata.X
    if hasattr(pred_X, "toarray"):
        pred_X = pred_X.toarray()
    gene_names = list(pred_adata.var_names)

    layers = {}
    for ct in type_names:
        ct_mask = pred_adata.obs["cell_type"].values == ct
        if ct_mask.sum() == 0:
            continue
        ct_adata = pred_adata[ct_mask]
        ct_X = ct_adata.X
        if hasattr(ct_X, "toarray"):
            ct_X = ct_X.toarray()

        # Sum per-cell values back to spots
        spot_sums = pd.DataFrame(0.0, index=spot_names, columns=gene_names)
        for i in range(ct_adata.n_obs):
            spot = ct_adata.obs.iloc[i]["spot_barcode"]
            if spot in spot_sums.index:
                spot_sums.loc[spot] += ct_X[i, :]
        layers[ct] = spot_sums

    return layers


# ---------------------------------------------------------------------------
# Layer loading
# ---------------------------------------------------------------------------


def load_existing_layers(layer_base_dir, rep_id, subfolder, type_names):
    """Load pre-existing spot-type layer CSVs.

    Args:
        layer_base_dir: Base directory (e.g., results/high_seg).
        rep_id: Replicate ID.
        subfolder: Subfolder under Wu_rep_{id}/ containing pass1 layers
                   (e.g., "citegeist" for QP, "citegeist_nb" for NB).
        type_names: List of canonical type names (with spaces/hyphens).

    Returns:
        Dict[type_name -> DataFrame (spots x genes)].
    """
    layer_dir = layer_base_dir / f"Wu_rep_{rep_id}" / subfolder / f"Wu_rep_{rep_id}_pass1" / "layers" / "pass1"

    layers = {}
    for type_name in type_names:
        # Convert spaces to underscores for filename
        file_stem = type_name.replace(" ", "_")
        csv_path = layer_dir / f"{file_stem}_layer_pass1.csv"
        if not csv_path.exists():
            logger.warning("Missing layer: %s", csv_path)
            continue
        # CSV is already spots x genes
        df = pd.read_csv(csv_path, index_col=0)
        layers[type_name] = df

    logger.info(
        "Loaded %d/%d type layers from %s",
        len(layers),
        len(type_names),
        layer_dir,
    )
    return layers


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main():
    """Run unified GEX benchmark for one simulation replicate."""
    parser = argparse.ArgumentParser(description="Unified GEX benchmark: SACE vs QP vs NB on simulation")
    parser.add_argument("--rep", type=int, required=True, help="Replicate ID (0-4)")
    parser.add_argument(
        "--output_dir", type=str, default="results_unified_gex", help="Output directory for results JSON"
    )
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    rep_id = args.rep

    # --- Load data ---
    spot_adata, cell_assignments, cell_spot_map, props, spot_coords = load_sim_region_data(rep_id)

    spot_X = spot_adata.X
    if hasattr(spot_X, "toarray"):
        spot_X = spot_X.toarray()
    spot_X = spot_X.astype(np.float64)

    gene_names = list(spot_adata.var_names)
    type_names = sorted(props.columns)
    spot_names = list(spot_adata.obs_names)
    N, G = spot_X.shape
    T = len(type_names)

    type_to_idx = {t: i for i, t in enumerate(type_names)}

    # Count cells per (spot, type) — shared across methods 1 and 2
    n_st = np.zeros((N, T), dtype=np.float64)
    for _, row in cell_spot_map.iterrows():
        cid = row["cell_id"]
        sidx = int(row["spot_idx"])
        ct = cell_assignments.get(cid)
        if ct is not None and ct in type_to_idx:
            n_st[sidx, type_to_idx[ct]] += 1

    # === Method 1: equal_split ===
    logger.info("--- Method 1/5: equal_split ---")
    props_arr = props.values
    t0 = time.time()
    cell_ids, cell_types_list, cell_spots_list = [], [], []
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
            cell_gex[i, :] = (props_arr[sidx, tidx] * spot_X[sidx, :]) / count
        cell_ids.append(cid)
        cell_types_list.append(ct)
        cell_spots_list.append(row["spot_barcode"])
    cell_gex = cell_gex[: len(cell_ids), :]
    equal_adata = sc.AnnData(
        X=cell_gex,
        obs=pd.DataFrame(
            {
                "cell_id": cell_ids,
                "spot_barcode": cell_spots_list,
                "cell_type": cell_types_list,
            }
        ),
        var=pd.DataFrame(index=gene_names),
    )
    equal_adata.obs_names = cell_ids
    t_equal = time.time() - t0
    logger.info("equal_split: %d cells in %.1f s", equal_adata.n_obs, t_equal)
    equal_layers = percell_to_spot_layers(equal_adata, spot_names, type_names)

    # === Method 2: type_reference ===
    logger.info("--- Method 2/5: type_reference ---")
    t0 = time.time()
    P = props.values
    beta, _, _, _ = np.linalg.lstsq(P, spot_X, rcond=None)
    beta = np.maximum(beta, 0.0)
    weights = P[:, :, None] * beta[None, :, :]
    weight_sums = np.maximum(weights.sum(axis=1, keepdims=True), 1e-12)
    E = weights / weight_sums * spot_X[:, None, :]

    cell_ids2, cell_types2, cell_spots2 = [], [], []
    cell_gex2 = np.zeros((len(cell_spot_map), G), dtype=np.float64)
    for i, (_, row) in enumerate(cell_spot_map.iterrows()):
        cid = row["cell_id"]
        sidx = int(row["spot_idx"])
        ct = cell_assignments.get(cid)
        if ct is None or ct not in type_to_idx:
            continue
        tidx = type_to_idx[ct]
        count = n_st[sidx, tidx]
        if count > 0:
            cell_gex2[i, :] = E[sidx, tidx, :] / count
        cell_ids2.append(cid)
        cell_types2.append(ct)
        cell_spots2.append(row["spot_barcode"])
    cell_gex2 = cell_gex2[: len(cell_ids2), :]
    typeref_adata = sc.AnnData(
        X=cell_gex2,
        obs=pd.DataFrame(
            {
                "cell_id": cell_ids2,
                "spot_barcode": cell_spots2,
                "cell_type": cell_types2,
            }
        ),
        var=pd.DataFrame(index=gene_names),
    )
    typeref_adata.obs_names = cell_ids2
    t_typeref = time.time() - t0
    logger.info("type_reference: %d cells in %.1f s", typeref_adata.n_obs, t_typeref)
    typeref_layers = percell_to_spot_layers(typeref_adata, spot_names, type_names)

    # === Method 3: sace_iter0 ===
    logger.info("--- Method 3/5: sace_iter0 ---")
    from CITEgeist.model.gex.sace_gex import run_sace

    t0 = time.time()
    _, sace_adata, diagnostics = run_sace(
        spot_counts=spot_X,
        proportions=props,
        cell_assignments=cell_assignments,
        cell_spot_map=cell_spot_map,
        spot_coords=spot_coords,
        gene_names=gene_names,
        max_iter=1,
    )
    t_sace = time.time() - t0
    logger.info("sace_iter0: %d cells in %.1f s", sace_adata.n_obs, t_sace)
    sace_layers = percell_to_spot_layers(sace_adata, spot_names, type_names)

    # === Method 4: qp_regression (load existing) ===
    logger.info("--- Method 4/5: qp_regression ---")
    qp_layers = load_existing_layers(QP_LAYER_BASE, rep_id, "citegeist", SIM_TYPES)

    # === Method 5: nb_emission (load existing) ===
    logger.info("--- Method 5/5: nb_emission ---")
    nb_layers = load_existing_layers(NB_LAYER_BASE, rep_id, "citegeist_nb", SIM_TYPES)

    # --- Evaluate all methods ---
    all_method_layers = {
        "equal_split": equal_layers,
        "type_reference": typeref_layers,
        "sace_iter0": sace_layers,
        "qp_regression": qp_layers,
        "nb_emission": nb_layers,
    }

    results = {"replicate": rep_id}
    for method_name, layers in all_method_layers.items():
        logger.info("Evaluating %s ...", method_name)
        eval_result = evaluate_sim_spot_layers(layers, rep_id, SIM_TYPES)
        results[method_name] = {
            "spot_layer_overall_r": eval_result["overall_r"],
            "spot_layer_per_type_r": eval_result["per_type_r"],
        }
        logger.info("  %s: spot_layer_r=%.4f", method_name, eval_result["overall_r"])

    # --- Save ---
    out_file = output_dir / f"simulation_unified_gex_rep_{rep_id}.json"
    with open(out_file, "w") as f:
        json.dump(results, f, indent=2, default=str)
    logger.info("Saved: %s", out_file)

    # --- Summary ---
    logger.info("=" * 60)
    logger.info("SUMMARY - Replicate %d", rep_id)
    logger.info("%-16s %12s", "Method", "spot_layer_r")
    logger.info("-" * 30)
    for m in all_method_layers:
        sl = results[m]["spot_layer_overall_r"]
        logger.info("%-16s %12.4f", m, sl)


if __name__ == "__main__":
    main()
