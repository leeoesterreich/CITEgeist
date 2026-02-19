"""
Test scResolve performance on marker genes only.

scResolve paper claimed ~90% correlation on top 20 marker genes per cell type.
This script tests that claim by:
1. Computing top 20 DE genes per cell type from GT
2. Re-evaluating all methods on just those marker genes
"""

import json
import logging
from pathlib import Path
from typing import Dict, List, Set

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

BASE_DIR = Path("/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking")
GT_DIR = Path("/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_pseudovisium/data_protein_gt/ground_truth_gex")

ACHIEVABLE_7_CELL_TYPES = [
    "B cells", "CD4+ T cells", "CD8+ T cells", "Macrophages",
    "Endothelial", "Epithelial", "Fibroblasts",
]

CELLTYPE_TO_FILENAME = {
    "B cells": "B_cells",
    "CD4+ T cells": "CD4pos_T_cells",
    "CD8+ T cells": "CD8pos_T_cells",
    "Macrophages": "Macrophages",
    "Endothelial": "Endothelial",
    "Epithelial": "Epithelial",
    "Fibroblasts": "Fibroblasts",
}


def compute_marker_genes(n_top: int = 20) -> Dict[str, List[str]]:
    """
    Compute top marker genes for each cell type using DE analysis on GT data.

    For each cell type, compares that cell type vs all others across all spots
    using a simple fold-change + expression filter.
    """
    logger.info(f"Computing top {n_top} marker genes per cell type from GT...")

    marker_genes = {}

    # Load all GT data for region 0 (representative)
    region_id = 0

    # Collect expression per cell type
    ct_expr = {}
    for ct in ACHIEVABLE_7_CELL_TYPES:
        ct_filename = CELLTYPE_TO_FILENAME.get(ct)
        gt_path = GT_DIR / f"Xenium_region_{region_id}" / f"{ct_filename}_GT.csv"
        if gt_path.exists():
            df = pd.read_csv(gt_path, index_col=0)  # genes x spots
            ct_expr[ct] = df

    # Get common genes
    all_genes = None
    for ct, df in ct_expr.items():
        if all_genes is None:
            all_genes = set(df.index)
        else:
            all_genes = all_genes.intersection(df.index)

    all_genes = sorted(all_genes)
    logger.info(f"Found {len(all_genes)} common genes across cell types")

    # For each cell type, find genes with highest expression relative to others
    for target_ct in ACHIEVABLE_7_CELL_TYPES:
        if target_ct not in ct_expr:
            continue

        target_df = ct_expr[target_ct].loc[all_genes]
        target_mean = target_df.mean(axis=1)  # Mean expression per gene

        # Background: mean across other cell types
        other_means = []
        for ct, df in ct_expr.items():
            if ct != target_ct:
                other_means.append(df.loc[all_genes].mean(axis=1))

        background_mean = pd.concat(other_means, axis=1).mean(axis=1)

        # Compute fold change (add pseudocount)
        fc = np.log2((target_mean + 1) / (background_mean + 1))

        # Filter: require minimum expression in target
        min_expr_threshold = target_mean.quantile(0.5)  # Top 50% expressed
        fc_filtered = fc[target_mean > min_expr_threshold]

        # Get top N by fold change
        top_genes = fc_filtered.nlargest(n_top).index.tolist()
        marker_genes[target_ct] = top_genes

        logger.info(f"{target_ct}: top markers = {top_genes[:5]}... (FC range: {fc_filtered.nlargest(n_top).min():.2f} - {fc_filtered.nlargest(n_top).max():.2f})")

    return marker_genes


def load_gt_gex(region_id: int, cell_type: str) -> pd.DataFrame:
    """Load GT GEX for a cell type (genes x spots)."""
    ct_filename = CELLTYPE_TO_FILENAME.get(cell_type)
    if ct_filename is None:
        return None
    gt_path = GT_DIR / f"Xenium_region_{region_id}" / f"{ct_filename}_GT.csv"
    if not gt_path.exists():
        return None
    return pd.read_csv(gt_path, index_col=0)


def load_scresolve_gex(region_id: int) -> Dict[str, pd.DataFrame]:
    """Load scResolve protein-gated GEX layers."""
    layers_dir = BASE_DIR / "scResolve" / "output_protein_gated" / f"Xenium_region_{region_id}" / "layers"

    if not layers_dir.exists():
        return {}

    ct_dfs = {}
    for ct in ACHIEVABLE_7_CELL_TYPES:
        ct_file = ct.replace(" ", "_")
        layer_file = layers_dir / f"{ct_file}_layer.csv"
        if layer_file.exists():
            df = pd.read_csv(layer_file, index_col=0)
            ct_dfs[ct] = df.T  # Transpose: (genes x spots) -> (spots x genes)

    return ct_dfs


def load_citegeist_gex(region_id: int, variant: str = "CITEgeist_Discrete") -> Dict[str, pd.DataFrame]:
    """Load CITEgeist GEX layers."""
    CITEGEIST_VARIANTS = {
        "CITEgeist_Continuous": "output/manual",
        "CITEgeist_Hybrid": "output/hybrid_cellpose",
        "CITEgeist_Discrete": "output_discrete_cellpose_fixed",
    }

    variant_path = CITEGEIST_VARIANTS.get(variant)
    sample_name = f"Xenium_region_{region_id}"

    if variant == "CITEgeist_Continuous":
        base_layers = BASE_DIR / "CITEgeist" / variant_path / f"{sample_name}_pass1" / "layers"
    else:
        base_layers = BASE_DIR / "CITEgeist" / variant_path / sample_name / f"{sample_name}_pass1" / "layers"

    layers_dir = base_layers / "pass1"
    if not layers_dir.exists():
        layers_dir = base_layers

    if not layers_dir.exists():
        return {}

    ct_dfs = {}
    for ct in ACHIEVABLE_7_CELL_TYPES:
        ct_file = ct.replace(" ", "_")
        for pattern in [f"{ct_file}_layer_pass1.csv", f"{ct_file}_layer.csv"]:
            layer_file = layers_dir / pattern
            if layer_file.exists():
                ct_dfs[ct] = pd.read_csv(layer_file, index_col=0)
                break

    return ct_dfs


def normalize_per_spot_cpm(vals: np.ndarray) -> np.ndarray:
    """Normalize each spot (column) to CPM across genes."""
    spot_totals = vals.sum(axis=0, keepdims=True)
    spot_totals[spot_totals == 0] = 1
    return vals / spot_totals * 1e6


def compute_gex_metrics_marker_only(
    pred: pd.DataFrame,
    gt: pd.DataFrame,
    marker_genes: List[str]
) -> Dict[str, float]:
    """
    Compute GEX metrics restricted to marker genes only.
    """
    pred_t = pred.T  # (genes x spots)

    # Case-insensitive gene matching
    pred_genes_lower = {g.lower(): g for g in pred_t.index}
    gt_genes_lower = {g.lower(): g for g in gt.index}
    marker_genes_lower = {g.lower() for g in marker_genes}

    # Intersection: common genes that are also markers
    common_lower = set(pred_genes_lower.keys()) & set(gt_genes_lower.keys()) & marker_genes_lower
    common_spots = pred_t.columns.intersection(gt.columns)

    if len(common_lower) == 0 or len(common_spots) == 0:
        return {"pearson_r": np.nan, "n_marker_genes": 0, "n_spots": 0}

    pred_genes = [pred_genes_lower[g] for g in common_lower]
    gt_genes = [gt_genes_lower[g] for g in common_lower]

    pred_aligned = pred_t.loc[pred_genes, common_spots].values
    gt_aligned = gt.loc[gt_genes, common_spots].values

    # CPM + log1p
    pred_log = np.log1p(normalize_per_spot_cpm(pred_aligned))
    gt_log = np.log1p(normalize_per_spot_cpm(gt_aligned))

    pred_flat = pred_log.flatten()
    gt_flat = gt_log.flatten()

    if np.std(pred_flat) == 0 or np.std(gt_flat) == 0:
        pearson_r = np.nan
    else:
        pearson_r, _ = stats.pearsonr(pred_flat, gt_flat)

    return {
        "pearson_r": float(pearson_r) if not np.isnan(pearson_r) else np.nan,
        "n_marker_genes": len(common_lower),
        "n_spots": len(common_spots),
    }


def main():
    print("=" * 80)
    print("MARKER GENE ANALYSIS: Testing scResolve's claimed 90% correlation")
    print("=" * 80)

    # Step 1: Compute marker genes
    marker_genes = compute_marker_genes(n_top=20)

    # Save marker genes for reference
    results_dir = BASE_DIR / "evaluation" / "results" / "method_comparison"
    with open(results_dir / "marker_genes_top20.json", "w") as f:
        json.dump(marker_genes, f, indent=2)

    # Step 2: Evaluate methods on marker genes only
    methods = {
        "scResolve": lambda rid: load_scresolve_gex(rid),
        "CITEgeist_Discrete": lambda rid: load_citegeist_gex(rid, "CITEgeist_Discrete"),
        "CITEgeist_Continuous": lambda rid: load_citegeist_gex(rid, "CITEgeist_Continuous"),
    }

    all_results = {}

    for method_name, loader_fn in methods.items():
        method_metrics = []

        for region_id in range(5):
            gex_layers = loader_fn(region_id)

            if not gex_layers:
                continue

            for ct in ACHIEVABLE_7_CELL_TYPES:
                if ct not in gex_layers or ct not in marker_genes:
                    continue

                gt = load_gt_gex(region_id, ct)
                if gt is None:
                    continue

                metrics = compute_gex_metrics_marker_only(
                    gex_layers[ct], gt, marker_genes[ct]
                )
                metrics["cell_type"] = ct
                metrics["region_id"] = region_id
                method_metrics.append(metrics)

        all_results[method_name] = method_metrics

    # Print results
    print("\n" + "=" * 80)
    print("RESULTS: Pearson r on TOP 20 MARKER GENES per cell type")
    print("=" * 80)

    print(f"\n{'Method':<25} {'Mean r (markers)':>18} {'Mean r (all genes)':>20}")
    print("-" * 65)

    # Load full comparison for reference
    with open(results_dir / "full_comparison_gex.json") as f:
        full_results = json.load(f)["results"]

    for method_name in ["scResolve", "CITEgeist_Discrete", "CITEgeist_Continuous"]:
        marker_metrics = all_results.get(method_name, [])
        full_metrics = full_results.get(method_name, [])

        if marker_metrics:
            marker_r = np.nanmean([m["pearson_r"] for m in marker_metrics])
        else:
            marker_r = np.nan

        if full_metrics:
            full_r = np.nanmean([m["pearson_r"] for m in full_metrics])
        else:
            full_r = np.nan

        print(f"{method_name:<25} {marker_r:>18.4f} {full_r:>20.4f}")

    # Per cell-type breakdown
    print("\n--- Per-Cell-Type Marker Gene Correlation ---")
    print(f"{'Cell Type':<20} {'scResolve':>12} {'CG_Discrete':>12} {'CG_Continuous':>14}")
    print("-" * 60)

    for ct in ACHIEVABLE_7_CELL_TYPES:
        row = f"{ct:<20}"
        for method in ["scResolve", "CITEgeist_Discrete", "CITEgeist_Continuous"]:
            metrics = [m for m in all_results.get(method, []) if m["cell_type"] == ct]
            if metrics:
                mean_r = np.nanmean([m["pearson_r"] for m in metrics])
                row += f" {mean_r:>11.3f}"
            else:
                row += f" {'N/A':>11}"
        print(row)

    # Save results
    with open(results_dir / "marker_genes_comparison.json", "w") as f:
        json.dump(all_results, f, indent=2, default=float)

    print(f"\nResults saved to: {results_dir}")
    print("=" * 80)


if __name__ == "__main__":
    main()
