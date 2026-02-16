"""
GEX deconvolution comparison across methods against protein-gated ground truth.

Evaluates per-cell-type gene expression layers from CITEgeist, Cell2Location,
Tangram, and scResolve. RCTD and Seurat only output proportions, so they are excluded.

Ground truth: per-cell-type summed gene expression at single-cell level,
aggregated to pseudo-Visium spots using protein-gated cell assignments.

scResolve: Uses the same hierarchical protein gating on scResolve multimodal-none
outputs to classify cells, then aggregates RNA per cell type per spot. This provides
a fair comparison using scResolve's joint RNA+protein super-resolution.
"""

import json
import logging
from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd
from scipy import stats

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

BASE_DIR = Path("/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking")
GT_DIR = Path("/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_pseudovisium/data_protein_gt/ground_truth_gex")

ACHIEVABLE_7_CELL_TYPES = [
    "B cells",
    "CD4+ T cells",
    "CD8+ T cells",
    "Macrophages",
    "Endothelial",
    "Epithelial",
    "Fibroblasts",
]

# Map cell type names to GT filenames
CELLTYPE_TO_FILENAME = {
    "B cells": "B_cells",
    "CD4+ T cells": "CD4pos_T_cells",
    "CD8+ T cells": "CD8pos_T_cells",
    "Macrophages": "Macrophages",
    "Endothelial": "Endothelial",
    "Epithelial": "Epithelial",
    "Fibroblasts": "Fibroblasts",
}


def load_gt_gex(region_id: int, cell_type: str) -> pd.DataFrame:
    """Load GT GEX for a cell type (genes x spots)."""
    ct_filename = CELLTYPE_TO_FILENAME.get(cell_type)
    if ct_filename is None:
        return None
    gt_path = GT_DIR / f"Xenium_region_{region_id}" / f"{ct_filename}_GT.csv"
    if not gt_path.exists():
        return None
    return pd.read_csv(gt_path, index_col=0)


def load_citegeist_gex(region_id: int) -> Dict[str, pd.DataFrame]:
    """Load CITEgeist GEX layers (spots x genes CSVs)."""
    # CITEgeist exports layers to {sample}_pass1/layers/pass1/
    sample_name = f"Xenium_region_{region_id}"
    # Updated path for asymmetric loss benchmark (output/manual)
    base_layers = BASE_DIR / "CITEgeist" / "output" / "manual" / f"{sample_name}_pass1" / "layers"

    # Check pass1 subdirectory first (standard export), then direct layers dir
    layers_dir = base_layers / "pass1"
    if not layers_dir.exists():
        layers_dir = base_layers

    if not layers_dir.exists():
        logger.warning(f"CITEgeist layers not found for region {region_id}: {layers_dir}")
        return {}

    ct_dfs = {}
    for ct in ACHIEVABLE_7_CELL_TYPES:
        # CITEgeist replaces spaces with _ but keeps + in filenames
        ct_file = ct.replace(" ", "_")
        for pattern in [
            f"{ct_file}_layer_pass1.csv",
            f"{ct_file}_layer.csv",
        ]:
            layer_file = layers_dir / pattern
            if layer_file.exists():
                ct_dfs[ct] = pd.read_csv(layer_file, index_col=0)
                break

    return ct_dfs


def load_citegeist_discrete_gex(region_id: int) -> Dict[str, pd.DataFrame]:
    """Load CITEgeist Discrete (Cellpose+IQP) GEX layers."""
    sample_name = f"Xenium_region_{region_id}"
    # v2 uses fixed coordinate conversion
    base_layers = BASE_DIR / "CITEgeist" / "output_discrete_cellpose_v2" / sample_name / f"{sample_name}_pass1" / "layers"

    # Check pass1 subdirectory first (standard export), then direct layers dir
    layers_dir = base_layers / "pass1"
    if not layers_dir.exists():
        layers_dir = base_layers

    if not layers_dir.exists():
        logger.warning(f"CITEgeist Discrete layers not found for region {region_id}: {layers_dir}")
        return {}

    ct_dfs = {}
    for ct in ACHIEVABLE_7_CELL_TYPES:
        # CITEgeist replaces spaces with _ but keeps + in filenames
        ct_file = ct.replace(" ", "_")
        for pattern in [
            f"{ct_file}_layer_pass1.csv",
            f"{ct_file}_layer.csv",
        ]:
            layer_file = layers_dir / pattern
            if layer_file.exists():
                ct_dfs[ct] = pd.read_csv(layer_file, index_col=0)
                break

    return ct_dfs


def load_competitor_gex(method: str, region_id: int) -> Dict[str, pd.DataFrame]:
    """Load C2L or Tangram GEX layers (spots x genes CSVs)."""
    layers_dir = BASE_DIR / method / "output_protein_gt" / f"Xenium_region_{region_id}" / "layers"

    if not layers_dir.exists():
        return {}

    ct_dfs = {}
    for ct in ACHIEVABLE_7_CELL_TYPES:
        layer_file = layers_dir / f"{ct}_layer.csv"
        if layer_file.exists():
            ct_dfs[ct] = pd.read_csv(layer_file, index_col=0)

    return ct_dfs


def load_scresolve_gex(region_id: int) -> Dict[str, pd.DataFrame]:
    """Load scResolve protein-gated GEX layers.

    scResolve outputs are generated by gate_and_aggregate_gex.py which applies
    the same hierarchical protein gating as create_protein_gt.py to scResolve's
    multimodal-none cells, then aggregates RNA per cell type per spot.

    Output format: (genes x spots) CSVs, same as GT format.
    """
    layers_dir = BASE_DIR / "scResolve" / "output_protein_gated" / f"Xenium_region_{region_id}" / "layers"

    if not layers_dir.exists():
        return {}

    ct_dfs = {}
    for ct in ACHIEVABLE_7_CELL_TYPES:
        # scResolve uses underscores but keeps + in filenames
        ct_file = ct.replace(" ", "_")
        layer_file = layers_dir / f"{ct_file}_layer.csv"
        if layer_file.exists():
            # scResolve outputs are (genes x spots), need to transpose to (spots x genes)
            # to match the expected input format for compute_gex_metrics
            df = pd.read_csv(layer_file, index_col=0)
            ct_dfs[ct] = df.T  # Transpose: (genes x spots) -> (spots x genes)

    return ct_dfs


def normalize_per_spot_cpm(vals: np.ndarray) -> np.ndarray:
    """Normalize each spot (column) to CPM across genes."""
    spot_totals = vals.sum(axis=0, keepdims=True)
    spot_totals[spot_totals == 0] = 1
    return vals / spot_totals * 1e6


def compute_gex_metrics(pred: pd.DataFrame, gt: pd.DataFrame) -> Dict[str, float]:
    """
    Compute GEX comparison metrics.

    pred: (spots x genes) — predictions
    gt: (genes x spots) — ground truth
    """
    # Transpose pred to (genes x spots) to match GT
    pred_t = pred.T

    # Case-insensitive gene matching
    pred_genes_lower = {g.lower(): g for g in pred_t.index}
    gt_genes_lower = {g.lower(): g for g in gt.index}
    common_genes_lower = set(pred_genes_lower.keys()) & set(gt_genes_lower.keys())

    common_spots = pred_t.columns.intersection(gt.columns)

    if len(common_genes_lower) == 0 or len(common_spots) == 0:
        return {"pearson_r": np.nan, "rmse": np.nan, "nrmse": np.nan, "mae": np.nan, "n_genes": 0, "n_spots": 0}

    pred_genes = [pred_genes_lower[g] for g in common_genes_lower]
    gt_genes = [gt_genes_lower[g] for g in common_genes_lower]

    pred_aligned = pred_t.loc[pred_genes, common_spots].values
    gt_aligned = gt.loc[gt_genes, common_spots].values

    # CPM normalization + log1p
    pred_log = np.log1p(normalize_per_spot_cpm(pred_aligned))
    gt_log = np.log1p(normalize_per_spot_cpm(gt_aligned))

    pred_flat = pred_log.flatten()
    gt_flat = gt_log.flatten()

    if np.std(pred_flat) == 0 or np.std(gt_flat) == 0:
        pearson_r = np.nan
    else:
        pearson_r, _ = stats.pearsonr(pred_flat, gt_flat)

    rmse = np.sqrt(np.mean((gt_flat - pred_flat) ** 2))
    mae = np.mean(np.abs(gt_flat - pred_flat))

    # NRMSE: normalize by range of ground truth
    gt_range = gt_flat.max() - gt_flat.min()
    nrmse = rmse / gt_range if gt_range > 0 else np.nan

    return {
        "pearson_r": float(pearson_r) if not np.isnan(pearson_r) else np.nan,
        "rmse": float(rmse),
        "nrmse": float(nrmse) if not np.isnan(nrmse) else np.nan,
        "mae": float(mae),
        "n_genes": len(common_genes_lower),
        "n_spots": len(common_spots),
    }


def main():
    print("=" * 80)
    print("GEX DECONVOLUTION COMPARISON - Protein-Gated Ground Truth")
    print("=" * 80)

    # Methods that produce GEX layers
    methods = {
        "CITEgeist": lambda rid: load_citegeist_gex(rid),
        "CITEgeist_Discrete": lambda rid: load_citegeist_discrete_gex(rid),
        "Cell2Location": lambda rid: load_competitor_gex("Cell2Location", rid),
        "Tangram": lambda rid: load_competitor_gex("Tangram", rid),
        "scResolve": lambda rid: load_scresolve_gex(rid),
    }

    all_results = {}

    for method_name, loader_fn in methods.items():
        method_metrics = []

        for region_id in range(5):
            gex_layers = loader_fn(region_id)

            if not gex_layers:
                logger.info(f"{method_name} region {region_id}: No GEX layers found")
                continue

            for ct in ACHIEVABLE_7_CELL_TYPES:
                if ct not in gex_layers:
                    continue

                gt = load_gt_gex(region_id, ct)
                if gt is None:
                    continue

                metrics = compute_gex_metrics(gex_layers[ct], gt)
                metrics["cell_type"] = ct
                metrics["region_id"] = region_id
                method_metrics.append(metrics)

        all_results[method_name] = method_metrics

    # Print summary
    print("\n--- Overall GEX Metrics ---")
    print(f"{'Method':<16} {'Pearson r':>10} {'RMSE':>10} {'NRMSE':>10} {'MAE':>10} {'Cell Types':>12} {'Regions':>8}")
    print("-" * 78)

    method_summaries = {}

    for method_name in ["CITEgeist", "CITEgeist_Discrete", "Cell2Location", "Tangram", "scResolve"]:
        metrics_list = all_results.get(method_name, [])
        if not metrics_list:
            print(f"{method_name:<16} {'N/A':>10} {'N/A':>10} {'N/A':>10} {'0':>12} {'0':>8}")
            continue

        df = pd.DataFrame(metrics_list)
        valid = df.dropna(subset=["pearson_r"])
        mean_r = valid["pearson_r"].mean()
        mean_rmse = valid["rmse"].mean()
        mean_nrmse = valid["nrmse"].mean()
        mean_mae = valid["mae"].mean()
        n_ct = valid["cell_type"].nunique()
        n_reg = valid["region_id"].nunique()

        print(f"{method_name:<16} {mean_r:>10.4f} {mean_rmse:>10.4f} {mean_nrmse:>10.4f} {mean_mae:>10.4f} {n_ct:>12} {n_reg:>8}")

        method_summaries[method_name] = {
            "mean_pearson_r": mean_r,
            "mean_rmse": mean_rmse,
            "mean_nrmse": mean_nrmse,
            "mean_mae": mean_mae,
            "n_cell_types": n_ct,
            "n_regions": n_reg,
        }

    # Per-cell-type table
    print("\n--- Per-Cell-Type GEX Pearson r (mean across regions) ---")
    sorted_methods = [m for m in ["CITEgeist", "Cell2Location", "Tangram", "scResolve"] if all_results.get(m)]

    header = f"{'Cell Type':<20}"
    for m in sorted_methods:
        header += f" {m:>16}"
    print(header)
    print("-" * (20 + 17 * len(sorted_methods)))

    for ct in ACHIEVABLE_7_CELL_TYPES:
        row = f"{ct:<20}"
        for m in sorted_methods:
            df = pd.DataFrame(all_results[m])
            ct_metrics = df[df["cell_type"] == ct]
            if len(ct_metrics) > 0:
                mean_r = ct_metrics["pearson_r"].dropna().mean()
                std_r = ct_metrics["pearson_r"].dropna().std()
                row += f" {mean_r:>7.3f}±{std_r:.3f}"
            else:
                row += f" {'N/A':>16}"
        print(row)

    # Per-region table
    print(f"\n--- Per-Region GEX Pearson r (mean across cell types) ---")
    header = f"{'Region':<10}"
    for m in sorted_methods:
        header += f" {m:>14}"
    print(header)
    print("-" * (10 + 15 * len(sorted_methods)))

    for region_id in range(5):
        row = f"Region {region_id:<3}"
        for m in sorted_methods:
            df = pd.DataFrame(all_results[m])
            reg = df[df["region_id"] == region_id]
            if len(reg) > 0:
                mean_r = reg["pearson_r"].dropna().mean()
                row += f" {mean_r:>14.4f}"
            else:
                row += f" {'N/A':>14}"
        print(row)

    # Save results
    results_dir = BASE_DIR / "evaluation" / "results" / "method_comparison"
    results_dir.mkdir(parents=True, exist_ok=True)

    def convert_numpy(obj):
        if isinstance(obj, (np.floating, np.float64)):
            return float(obj)
        elif isinstance(obj, (np.integer, np.int64)):
            return int(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, dict):
            return {k: convert_numpy(v) for k, v in obj.items()}
        elif isinstance(obj, list):
            return [convert_numpy(v) for v in obj]
        return obj

    with open(results_dir / "full_comparison_gex.json", "w") as f:
        json.dump(convert_numpy(all_results), f, indent=2)

    print(f"\nResults saved to: {results_dir / 'full_comparison_gex.json'}")
    print("=" * 80)


if __name__ == "__main__":
    main()
