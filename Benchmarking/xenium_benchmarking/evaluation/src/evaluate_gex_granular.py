"""
Evaluate GEX deconvolution across methods using granular 10-cell-type ground truth.

Compares CITEgeist, Cell2Location, and Tangram GEX predictions against
RNA-based per-cell-type gene expression ground truth.
"""

import json
import logging
from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import mean_squared_error, mean_absolute_error

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

# Constants
BASE_DIR = Path("/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking")
GT_DIR = BASE_DIR / "xenium_pseudovisium" / "data_granular_gt" / "ground_truth_gex"
BENCHMARK_DIR = BASE_DIR / "xenium_benchmarking"

# Cell type name mapping for file compatibility
CELLTYPE_TO_FILENAME = {
    "B cells": "B_cells",
    "CD8+ T cells": "CD8pos_T_cells",
    "Endothelial": "Endothelial",
    "Epithelial": "Epithelial",
    "Macrophages": "Macrophages",
    "Mixed Immune": "Mixed_Immune",
    "Myofibroblasts": "Myofibroblasts",
    "Proliferating T": "Proliferating_T",
    "Stromal": "Stromal",
    "Vascular Stromal": "Vascular_Stromal",
}

GRANULAR_CELL_TYPES = list(CELLTYPE_TO_FILENAME.keys())


def load_gt_gex(region_id: int, cell_type: str) -> pd.DataFrame:
    """Load ground truth GEX for a cell type in a region."""
    ct_filename = CELLTYPE_TO_FILENAME.get(cell_type, cell_type.replace(" ", "_").replace("+", "pos"))
    gt_path = GT_DIR / f"Xenium_region_{region_id}" / f"{ct_filename}_GT.csv"

    if not gt_path.exists():
        return None

    df = pd.read_csv(gt_path, index_col=0)
    return df


def load_citegeist_gex(region_id: int, profile_mapping: Dict[str, str]) -> Dict[str, pd.DataFrame]:
    """
    Load CITEgeist GEX and aggregate profiles to cell types using mapping.

    Args:
        region_id: Region index
        profile_mapping: Dict mapping Profile_X -> cell_type

    Returns:
        Dict mapping cell_type -> aggregated GEX DataFrame
    """
    layers_dir = (BENCHMARK_DIR / "CITEgeist" / "output_autodiscovery" /
                  f"Xenium_region_{region_id}_pass1" / "layers" / "pass1")

    if not layers_dir.exists():
        return {}

    # Load all profile layers and aggregate by cell type
    ct_dfs = {}

    for profile_name, cell_type in profile_mapping.items():
        if cell_type == "Unknown":
            continue

        layer_file = layers_dir / f"{profile_name}_layer_pass1.csv"
        if not layer_file.exists():
            continue

        df = pd.read_csv(layer_file, index_col=0)

        if cell_type not in ct_dfs:
            ct_dfs[cell_type] = df.copy()
        else:
            # Sum multiple profiles for same cell type
            ct_dfs[cell_type] = ct_dfs[cell_type].add(df, fill_value=0)

    return ct_dfs


def load_cell2location_gex(region_id: int) -> Dict[str, pd.DataFrame]:
    """Load Cell2Location GEX layers."""
    layers_dir = BENCHMARK_DIR / "Cell2Location" / "output_granular" / f"Xenium_region_{region_id}" / "layers"

    if not layers_dir.exists():
        return {}

    ct_dfs = {}
    for ct in GRANULAR_CELL_TYPES:
        layer_file = layers_dir / f"{ct}_layer.csv"
        if layer_file.exists():
            ct_dfs[ct] = pd.read_csv(layer_file, index_col=0)

    return ct_dfs


def load_tangram_gex(region_id: int) -> Dict[str, pd.DataFrame]:
    """Load Tangram GEX layers."""
    layers_dir = BENCHMARK_DIR / "Tangram" / "output_granular" / f"Xenium_region_{region_id}" / "layers"

    if not layers_dir.exists():
        return {}

    ct_dfs = {}
    for ct in GRANULAR_CELL_TYPES:
        layer_file = layers_dir / f"{ct}_layer.csv"
        if layer_file.exists():
            ct_dfs[ct] = pd.read_csv(layer_file, index_col=0)

    return ct_dfs


def _normalize_per_spot_cpm(vals):
    """Normalize each spot (column) to CPM across genes, handling zeros."""
    spot_totals = vals.sum(axis=0, keepdims=True)
    spot_totals[spot_totals == 0] = 1
    return vals / spot_totals * 1e6


def compute_gex_metrics(pred: pd.DataFrame, gt: pd.DataFrame) -> Dict[str, float]:
    """
    Compute GEX comparison metrics between prediction and ground truth.

    Both pred and GT are normalized per-spot to CPM before log1p and metric
    computation, removing library-size mismatches between methods.

    Args:
        pred: Predicted GEX (spots x genes) - needs transpose
        gt: Ground truth GEX (genes x spots)

    Returns:
        Dict with metrics
    """
    # Predictions are (spots x genes), GT is (genes x spots)
    # Transpose predictions to match GT format
    pred_t = pred.T  # Now (genes x spots)

    # Try case-insensitive gene matching (Tangram uses lowercase)
    pred_genes_lower = {g.lower(): g for g in pred_t.index}
    gt_genes_lower = {g.lower(): g for g in gt.index}

    common_genes_lower = set(pred_genes_lower.keys()).intersection(set(gt_genes_lower.keys()))

    if len(common_genes_lower) == 0:
        # Fallback to exact match
        common_genes = pred_t.index.intersection(gt.index)
        common_spots = pred_t.columns.intersection(gt.columns)

        if len(common_genes) == 0 or len(common_spots) == 0:
            return {"pearson_r": np.nan, "rmse": np.nan, "mae": np.nan, "n_genes": 0, "n_spots": 0}

        pred_aligned = pred_t.loc[common_genes, common_spots]
        gt_aligned = gt.loc[common_genes, common_spots]
    else:
        # Use case-insensitive matching
        common_spots = pred_t.columns.intersection(gt.columns)
        if len(common_spots) == 0:
            return {"pearson_r": np.nan, "rmse": np.nan, "mae": np.nan, "n_genes": 0, "n_spots": 0}

        # Build aligned DataFrames using mapped gene names
        pred_genes = [pred_genes_lower[g] for g in common_genes_lower]
        gt_genes = [gt_genes_lower[g] for g in common_genes_lower]

        pred_aligned = pred_t.loc[pred_genes, common_spots]
        gt_aligned = gt.loc[gt_genes, common_spots]

        # Rename pred_aligned index to match gt for proper alignment
        pred_aligned.index = gt_genes

    # Per-spot CPM normalization then log1p
    pred_log = np.log1p(_normalize_per_spot_cpm(pred_aligned.values))
    gt_log = np.log1p(_normalize_per_spot_cpm(gt_aligned.values))

    pred_flat = pred_log.flatten()
    gt_flat = gt_log.flatten()

    # Handle edge cases
    if np.std(pred_flat) == 0 or np.std(gt_flat) == 0:
        pearson_r = np.nan
    else:
        pearson_r, _ = stats.pearsonr(pred_flat, gt_flat)

    rmse = np.sqrt(mean_squared_error(gt_flat, pred_flat))
    mae = mean_absolute_error(gt_flat, pred_flat)

    n_genes = len(pred_aligned)  # Number of genes after alignment

    return {
        "pearson_r": pearson_r,
        "rmse": rmse,
        "mae": mae,
        "n_genes": n_genes,
        "n_spots": len(common_spots),
    }


def load_profile_mappings() -> Dict[int, Dict[str, str]]:
    """Load profile mappings from proportion evaluation results."""
    results_file = BENCHMARK_DIR / "xenium_benchmarking" / "evaluation" / "results_granular" / "full_results.json"

    if not results_file.exists():
        # Try alternate path
        results_file = BENCHMARK_DIR.parent / "xenium_benchmarking" / "evaluation" / "results_granular" / "full_results.json"

    if not results_file.exists():
        logger.warning("Could not find proportion evaluation results for profile mappings")
        return {}

    with open(results_file) as f:
        data = json.load(f)

    mappings = {}
    for region_data in data["methods"]["CITEgeist"]["regions"]:
        region_id = region_data["region_id"]
        if "profile_mapping" in region_data:
            mappings[region_id] = region_data["profile_mapping"]

    return mappings


def main():
    """Run GEX evaluation across all methods."""
    logger.info("=" * 60)
    logger.info("GEX Evaluation: Granular 10 Cell Type Ground Truth")
    logger.info("=" * 60)

    # Load profile mappings for CITEgeist
    profile_mappings = load_profile_mappings()
    logger.info(f"Loaded profile mappings for {len(profile_mappings)} regions")

    results = {
        "CITEgeist": [],
        "Cell2Location": [],
        "Tangram": [],
    }

    for region_id in range(5):
        logger.info(f"\nRegion {region_id}:")

        # Load method predictions
        citegeist_gex = load_citegeist_gex(region_id, profile_mappings.get(region_id, {}))
        cell2loc_gex = load_cell2location_gex(region_id)
        tangram_gex = load_tangram_gex(region_id)

        logger.info(f"  CITEgeist: {len(citegeist_gex)} cell types")
        logger.info(f"  Cell2Location: {len(cell2loc_gex)} cell types")
        logger.info(f"  Tangram: {len(tangram_gex)} cell types")

        for ct in GRANULAR_CELL_TYPES:
            gt = load_gt_gex(region_id, ct)
            if gt is None:
                continue

            # CITEgeist
            if ct in citegeist_gex:
                metrics = compute_gex_metrics(citegeist_gex[ct], gt)
                metrics["cell_type"] = ct
                metrics["region_id"] = region_id
                results["CITEgeist"].append(metrics)

            # Cell2Location
            if ct in cell2loc_gex:
                metrics = compute_gex_metrics(cell2loc_gex[ct], gt)
                metrics["cell_type"] = ct
                metrics["region_id"] = region_id
                results["Cell2Location"].append(metrics)

            # Tangram
            if ct in tangram_gex:
                metrics = compute_gex_metrics(tangram_gex[ct], gt)
                metrics["cell_type"] = ct
                metrics["region_id"] = region_id
                results["Tangram"].append(metrics)

    # Print summary
    logger.info("\n" + "=" * 60)
    logger.info("GEX EVALUATION SUMMARY")
    logger.info("=" * 60)

    summary_data = []

    for method, method_results in results.items():
        if not method_results:
            logger.info(f"\n{method}: No results")
            continue

        df = pd.DataFrame(method_results)

        # Overall metrics
        mean_r = df["pearson_r"].dropna().mean()
        std_r = df["pearson_r"].dropna().std()
        mean_rmse = df["rmse"].dropna().mean()
        mean_mae = df["mae"].dropna().mean()
        n_cell_types = df["cell_type"].nunique()

        logger.info(f"\n{method}:")
        logger.info(f"  Pearson r: {mean_r:.4f} +/- {std_r:.4f}")
        logger.info(f"  RMSE: {mean_rmse:.4f}")
        logger.info(f"  MAE: {mean_mae:.4f}")
        logger.info(f"  Cell types evaluated: {n_cell_types}")

        # Per-cell-type breakdown
        logger.info("\n  Per-cell-type Pearson r:")
        ct_means = df.groupby("cell_type")["pearson_r"].mean().sort_values(ascending=False)
        for ct, r in ct_means.items():
            logger.info(f"    {ct}: {r:.4f}")

        summary_data.append({
            "method": method,
            "mean_pearson_r": mean_r,
            "std_pearson_r": std_r,
            "mean_rmse": mean_rmse,
            "mean_mae": mean_mae,
            "n_cell_types": n_cell_types,
            "n_evaluations": len(method_results),
        })

    # Save results
    output_dir = BENCHMARK_DIR / "xenium_benchmarking" / "evaluation" / "results_granular"
    output_dir.mkdir(parents=True, exist_ok=True)

    # Detailed results
    for method, method_results in results.items():
        if method_results:
            df = pd.DataFrame(method_results)
            df.to_csv(output_dir / f"{method}_gex_metrics.csv", index=False)

    # Summary
    summary_df = pd.DataFrame(summary_data)
    summary_df.to_csv(output_dir / "gex_summary.csv", index=False)

    logger.info(f"\nResults saved to {output_dir}")
    logger.info("=" * 60)

    return results, summary_data


if __name__ == "__main__":
    main()
