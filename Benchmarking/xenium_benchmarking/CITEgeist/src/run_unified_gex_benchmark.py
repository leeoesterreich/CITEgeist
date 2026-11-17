#!/usr/bin/env python
"""Unified GEX benchmark: SACE vs QP vs NB on Xenium pseudo-Visium.

Runs 6 GEX deconvolution methods on one Xenium pseudo-Visium region and
evaluates all through identical percell_r + spot_layer_r pipelines.

Methods:
  1. equal_split     -- divide spot counts equally within type (baseline)
  2. type_reference  -- lstsq type profiles then equal split (baseline)
  3. sace_iter0      -- SACE with max_iter=1 (production default)
  4. sace_full       -- SACE with max_iter=10
  5. qp_legacy       -- load pre-existing QP layer CSVs from output_singler_7type/
  6. nb_emission     -- load pre-existing NB layer CSVs from nb_vs_qp_gex/

All methods are evaluated at 6-type level (CD4+/CD8+ merged to T cells).

Usage:
    python run_unified_gex_benchmark.py --region 0 --output_dir results_unified_gex/
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
from scipy.stats import pearsonr  # noqa: F401 (used in imported functions)

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
REPO_ROOT = Path(str(_cg.REPO_ROOT))
sys.path.insert(0, str(REPO_ROOT))

# Add the src directory to sys.path so we can import run_sace_benchmark
_SRC_DIR = Path(__file__).parent
sys.path.insert(0, str(_SRC_DIR))

BENCH_DIR = REPO_ROOT / "Benchmarking" / "xenium_benchmarking" / "CITEgeist"

# Reuse everything from the SACE benchmark
from run_sace_benchmark import (  # noqa: E402
    load_region_data,
    load_xenium_singlecell_gex,
    baseline_equal_split,
    baseline_type_reference,
    run_sace_iter0,
    evaluate_per_cell,
    evaluate_spot_type_layers,
    _collect_method_results,
    DATA_DIR,
)

# Pre-existing layer directories
QP_LAYER_DIR = BENCH_DIR / "output_singler_7type"
NB_LAYER_DIR = BENCH_DIR / "_archive" / "output_experimental" / "nb_vs_qp_gex"

# 6-type system after merge
TYPES_6 = ["B cells", "T cells", "Macrophages", "Endothelial", "Epithelial", "Fibroblasts"]

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
)
logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# T-cell merge helper
# ---------------------------------------------------------------------------


def merge_tcells_in_adata(adata):
    """Remap CD4+/CD8+ T cell labels to 'T cells' in obs["cell_type"].

    Args:
        adata: AnnData with obs["cell_type"] column.

    Returns:
        Same AnnData with cell_type relabeled in-place.
    """
    ct = adata.obs["cell_type"].copy()
    mask = ct.str.contains("CD4|CD8", case=False, na=False)
    ct[mask] = "T cells"
    adata.obs["cell_type"] = ct
    return adata


# ---------------------------------------------------------------------------
# Layer CSV -> per-cell AnnData converter
# ---------------------------------------------------------------------------


def layers_to_percell_adata(layer_dir, region_id, cell_assignments, cell_spot_map, spot_adata, type_map):
    """Convert spot-type layer CSVs to per-cell AnnData.

    Loads pre-existing layer CSVs (genes x spots format), transposes to
    spots x genes, merges CD4+/CD8+ T cell layers by summing, then
    distributes to individual cells using equal-split-within-type.

    Args:
        layer_dir: Path containing {CellType}_layer_pass1.csv files.
        region_id: Region ID (for logging).
        cell_assignments: Dict[cell_id_str -> type_name] (already 6-type).
        cell_spot_map: DataFrame with cell_id, spot_barcode, spot_idx, x, y.
        spot_adata: AnnData for spot name/order mapping.
        type_map: Dict[csv_filename_stem -> canonical_type_name], e.g.
                  {"B_cells": "B cells", "T_cells": "T cells", ...}

    Returns:
        AnnData (n_cells, G) with obs columns: cell_id, spot_barcode, cell_type.
    """
    spot_names = list(spot_adata.obs_names)

    # Load all layer CSVs into a dict: canonical_type_name -> DataFrame (spots x genes)
    type_layers = {}
    gene_names = None
    for csv_stem, canon_name in type_map.items():
        csv_path = layer_dir / f"{csv_stem}_layer_pass1.csv"
        if not csv_path.exists():
            logger.warning("Missing layer file: %s", csv_path)
            continue
        # CSV is already spots x genes
        df = pd.read_csv(csv_path, index_col=0)
        # Align to spot_adata spot order (fill missing spots with 0)
        df = df.reindex(spot_names, fill_value=0.0)
        if gene_names is None:
            gene_names = list(df.columns)
        type_layers[canon_name] = df

    if not type_layers or gene_names is None:
        raise ValueError(f"No layer files loaded from {layer_dir} for region {region_id}")

    # Merge CD4+/CD8+ T cell layers if both exist (QP 7-type case)
    if "CD4+ T cells" in type_layers and "CD8+ T cells" in type_layers:
        cd4 = type_layers.pop("CD4+ T cells")
        cd8 = type_layers.pop("CD8+ T cells")
        common_genes = sorted(set(cd4.columns) & set(cd8.columns))
        type_layers["T cells"] = cd4[common_genes].add(cd8[common_genes], fill_value=0.0)
        gene_names = common_genes

    N = len(spot_names)
    G = len(gene_names)
    type_names = sorted(type_layers.keys())
    type_to_idx = {t: i for i, t in enumerate(type_names)}

    # Build layer array: (N, T, G)
    layer_arr = np.zeros((N, len(type_names), G), dtype=np.float64)
    for t_name, df in type_layers.items():
        tidx = type_to_idx[t_name]
        layer_arr[:, tidx, :] = df[gene_names].values

    # Count cells per (spot, type)
    n_st = np.zeros((N, len(type_names)), dtype=np.float64)
    for _, row in cell_spot_map.iterrows():
        cid = row["cell_id"]
        sidx = int(row["spot_idx"])
        ct = cell_assignments.get(cid)
        if ct is not None and ct in type_to_idx:
            n_st[sidx, type_to_idx[ct]] += 1

    # Allocate per cell: x_c = layer[s, t, :] / n_{s,t}
    cell_ids = []
    cell_types = []
    cell_spots = []
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
            cell_gex[i, :] = layer_arr[sidx, tidx, :] / count
        cell_ids.append(cid)
        cell_types.append(ct)
        cell_spots.append(row["spot_barcode"])

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
    logger.info("Loaded %d cells from layer CSVs in %s", len(cell_ids), layer_dir)
    return adata


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser(description="Unified GEX benchmark: SACE vs QP vs NB on Xenium pseudo-Visium")
    parser.add_argument("--region", type=int, required=True, help="Region ID (0-4)")
    parser.add_argument(
        "--output_dir", type=str, default="results_unified_gex", help="Output directory for results JSON"
    )
    parser.add_argument("--max_iter", type=int, default=10, help="Max SACE EM iterations for sace_full")
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    region_id = args.region

    # --- Load data ---
    logger.info("Loading region %d data ...", region_id)
    spot_adata, cell_assignments, cell_spot_map, props, region_cells = load_region_data(region_id)

    spot_X = spot_adata.X
    if hasattr(spot_X, "toarray"):
        spot_X = spot_X.toarray()
    spot_X = spot_X.astype(np.float64)

    gene_names = list(spot_adata.var_names)
    type_names = sorted(props.columns)
    spot_coords = np.column_stack([spot_adata.obs["x"].values, spot_adata.obs["y"].values])

    # --- Load GT single-cell GEX ---
    gt_adata = load_xenium_singlecell_gex(list(cell_assignments.keys()))

    # --- Build 6-type cell assignments (merge T cell subtypes) ---
    cell_assignments_6type = {}
    for cid, ct in cell_assignments.items():
        if "CD4" in ct or "CD8" in ct:
            cell_assignments_6type[cid] = "T cells"
        else:
            cell_assignments_6type[cid] = ct

    logger.info("=" * 70)
    logger.info("Running 6 methods on region %d", region_id)
    logger.info("=" * 70)

    # === Method 1: equal_split ===
    logger.info("--- Method 1/6: equal_split ---")
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
    merge_tcells_in_adata(equal_adata)
    logger.info("equal_split: %d cells in %.1f s", equal_adata.n_obs, t_equal)

    # === Method 2: type_reference ===
    logger.info("--- Method 2/6: type_reference ---")
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
    merge_tcells_in_adata(typeref_adata)
    logger.info("type_reference: %d cells in %.1f s", typeref_adata.n_obs, t_typeref)

    # === Method 3: sace_iter0 ===
    logger.info("--- Method 3/6: sace_iter0 ---")
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
    merge_tcells_in_adata(iter0_adata)
    logger.info("sace_iter0: %d cells in %.1f s", iter0_adata.n_obs, t_iter0)

    # === Method 4: sace_full ===
    logger.info("--- Method 4/6: sace_full ---")
    from CITEgeist.model.gex.sace_gex import run_sace  # noqa: E402

    t0 = time.time()
    _, sace_adata, diagnostics = run_sace(
        spot_counts=spot_X,
        proportions=props,
        cell_assignments=cell_assignments,
        cell_spot_map=cell_spot_map,
        spot_coords=spot_coords,
        gene_names=gene_names,
        max_iter=args.max_iter,
    )
    t_sace = time.time() - t0
    merge_tcells_in_adata(sace_adata)
    logger.info("sace_full: %d cells in %.1f s", sace_adata.n_obs, t_sace)

    # === Method 5: qp_legacy (load existing layers) ===
    logger.info("--- Method 5/6: qp_legacy ---")
    qp_layer_path = QP_LAYER_DIR / f"Xenium_region_{region_id}_pass1" / "layers" / "pass1"
    qp_type_map = {
        "B_cells": "B cells",
        "CD4+_T_cells": "CD4+ T cells",
        "CD8+_T_cells": "CD8+ T cells",
        "Macrophages": "Macrophages",
        "Endothelial": "Endothelial",
        "Epithelial": "Epithelial",
        "Fibroblasts": "Fibroblasts",
    }
    t0 = time.time()
    qp_adata = layers_to_percell_adata(
        qp_layer_path,
        region_id,
        cell_assignments_6type,
        cell_spot_map,
        spot_adata,
        qp_type_map,
    )
    t_qp = time.time() - t0
    logger.info("qp_legacy: %d cells in %.1f s", qp_adata.n_obs, t_qp)

    # === Method 6: nb_emission (load existing layers) ===
    logger.info("--- Method 6/6: nb_emission ---")
    nb_layer_path = (
        NB_LAYER_DIR / f"Xenium_region_{region_id}" / "gex_nb" / f"Xenium_region_{region_id}_pass1" / "layers" / "pass1"
    )
    nb_type_map = {
        "B_cells": "B cells",
        "T_cells": "T cells",
        "Macrophages": "Macrophages",
        "Endothelial": "Endothelial",
        "Epithelial": "Epithelial",
        "Fibroblasts": "Fibroblasts",
    }
    t0 = time.time()
    nb_adata = layers_to_percell_adata(
        nb_layer_path,
        region_id,
        cell_assignments_6type,
        cell_spot_map,
        spot_adata,
        nb_type_map,
    )
    t_nb = time.time() - t0
    logger.info("nb_emission: %d cells in %.1f s", nb_adata.n_obs, t_nb)

    # --- Evaluate all methods (merge GT T cells for percell) ---
    gt_adata_merged = gt_adata.copy()
    # GT adata from load_xenium_singlecell_gex() has no cell_type column;
    # add it from cell_assignments before merging T cells
    gt_adata_merged.obs["cell_type"] = gt_adata_merged.obs_names.map(lambda x: cell_assignments.get(x, "Unknown"))
    merge_tcells_in_adata(gt_adata_merged)

    # Common genes across ALL 6 methods + GT
    all_pred_genes = (
        set(equal_adata.var_names)
        & set(typeref_adata.var_names)
        & set(iter0_adata.var_names)
        & set(sace_adata.var_names)
        & set(qp_adata.var_names)
        & set(nb_adata.var_names)
        & set(gt_adata_merged.var_names)
    )
    common_genes = sorted(all_pred_genes)
    logger.info("Common genes across all 6 methods + GT: %d", len(common_genes))

    methods = {
        "equal_split": equal_adata,
        "type_reference": typeref_adata,
        "sace_iter0": iter0_adata,
        "sace_full": sace_adata,
        "qp_legacy": qp_adata,
        "nb_emission": nb_adata,
    }

    logger.info("=" * 70)
    logger.info("Evaluating all methods on region %d", region_id)
    logger.info("=" * 70)

    results = {
        "region": region_id,
        "n_common_genes": len(common_genes),
        "sace_params": {"max_iter": args.max_iter},
        "sace_diagnostics": {
            "n_iterations": diagnostics.get("n_iterations"),
            "damping_activated": diagnostics.get("damping_activated"),
            "final_bandwidth": diagnostics.get("final_bandwidth"),
            "log_likelihood_trace": diagnostics.get("log_likelihood_trace"),
        },
        "runtime_seconds": {
            "equal_split": round(t_equal, 1),
            "type_reference": round(t_typeref, 1),
            "sace_iter0": round(t_iter0, 1),
            "sace_full": round(t_sace, 1),
            "qp_legacy": round(t_qp, 1),
            "nb_emission": round(t_nb, 1),
        },
    }

    for method_name, pred_adata in methods.items():
        logger.info("Evaluating %s ...", method_name)
        results[method_name] = _collect_method_results(
            method_name,
            pred_adata,
            gt_adata_merged,
            common_genes,
            region_id,
        )

    # --- Save ---
    out_file = output_dir / f"xenium_unified_gex_region_{region_id}.json"
    with open(out_file, "w") as f:
        json.dump(results, f, indent=2, default=str)
    logger.info("Saved: %s", out_file)

    # --- Summary table ---
    logger.info("=" * 70)
    logger.info("SUMMARY - Region %d (%d genes)", region_id, len(common_genes))
    logger.info("%-16s %12s %12s", "Method", "percell_r", "spot_layer_r")
    logger.info("-" * 42)
    for method_name in methods:
        pc = results[method_name].get("percell_mean_gene_r", 0.0)
        sl = results[method_name].get("spot_layer_overall_r", 0.0)
        logger.info("%-16s %12.4f %12.4f", method_name, pc, sl)
    logger.info("=" * 70)


if __name__ == "__main__":
    main()
