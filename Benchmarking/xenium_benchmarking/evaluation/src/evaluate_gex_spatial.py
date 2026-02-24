"""
Evaluate GEX deconvolution with spatial-aware metrics.

Measures spatial heterogeneity preservation, not just overall correlation:
- Per-spot correlation: Does each spot have correct gene profile?
- Per-gene correlation: Is spatial pattern of each gene preserved?
- Variance ratio: Is heterogeneity preserved or crushed to mean?
"""

import json
import logging
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
from scipy import stats

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

# Constants
BASE_DIR = Path("/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking")
GT_DIR = BASE_DIR / "xenium_pseudovisium" / "data_protein_gt_backup_20260209_200315" / "ground_truth_gex"
OUTPUT_DIR = BASE_DIR / "xenium_benchmarking" / "CITEgeist" / "output"

CELL_TYPES = [
    "B_cells", "CD4pos_T_cells", "CD8pos_T_cells", "Endothelial",
    "Epithelial", "Fibroblasts", "Macrophages"
]

# Marker genes for each cell type (canonical markers present in Xenium panel)
MARKER_GENES = {
    "B_cells": ["CD19", "MS4A1", "CD79A", "CD79B", "BANK1"],
    "CD4pos_T_cells": ["CD3E", "CD3D", "CD4", "IL7R", "CD40LG"],
    "CD8pos_T_cells": ["CD3E", "CD3D", "CD8A", "CD8B", "GZMB", "PRF1", "NKG7"],
    "Endothelial": ["PECAM1", "VWF", "CDH5", "CLDN5", "ESAM"],
    "Epithelial": ["EPCAM", "KRT19", "KRT8", "KRT18", "CDH1"],
    "Fibroblasts": ["COL1A1", "COL1A2", "DCN", "LUM", "PDGFRA", "FAP"],
    "Macrophages": ["CD68", "CD163", "CSF1R", "MARCO", "MRC1", "C1QA", "C1QB"],
}


def load_gt_gex(region_id: int, cell_type: str) -> pd.DataFrame:
    """Load ground truth GEX for a cell type. Returns (genes x spots)."""
    gt_path = GT_DIR / f"Xenium_region_{region_id}" / f"{cell_type}_GT.csv"
    if not gt_path.exists():
        return None
    return pd.read_csv(gt_path, index_col=0)


def load_citegeist_gex(region_id: int, output_subdir: str) -> Dict[str, pd.DataFrame]:
    """
    Load CITEgeist deconvolved GEX layers.

    Returns dict of cell_type -> DataFrame (spots x genes)
    """
    # Try multiple path patterns (different benchmark runs use different structures)
    possible_paths = [
        OUTPUT_DIR / output_subdir / f"Xenium_region_{region_id}_pass1" / "layers" / "pass1",
        OUTPUT_DIR / output_subdir / f"Xenium_region_{region_id}" / f"Xenium_region_{region_id}_pass1" / "layers" / "pass1",
    ]

    layers_dir = None
    for path in possible_paths:
        if path.exists():
            layers_dir = path
            break

    if layers_dir is None:
        logger.warning(f"Layers not found for region {region_id} in {output_subdir}")
        return {}

    ct_dfs = {}
    for layer_file in layers_dir.glob("*_layer_pass1.csv"):
        # Extract cell type from filename (e.g., "B_cells_layer_pass1.csv")
        ct_name = layer_file.stem.replace("_layer_pass1", "")
        # Normalize cell type names: "CD8+_T_cells" -> "CD8pos_T_cells"
        ct_name = ct_name.replace("+", "pos")
        df = pd.read_csv(layer_file, index_col=0)  # (spots x genes)
        ct_dfs[ct_name] = df

    return ct_dfs


def normalize_cpm(arr: np.ndarray, axis: int = 1) -> np.ndarray:
    """Normalize to CPM along specified axis (1=per-row/spot, 0=per-column/gene)."""
    totals = arr.sum(axis=axis, keepdims=True)
    totals = np.where(totals == 0, 1, totals)
    return arr / totals * 1e6


def compute_spatial_metrics(
    pred: pd.DataFrame,
    gt: pd.DataFrame,
    min_expression: float = 1.0
) -> Dict[str, float]:
    """
    Compute spatial-aware GEX metrics.

    Args:
        pred: Predicted GEX (spots x genes)
        gt: Ground truth GEX (genes x spots)
        min_expression: Minimum total expression to include gene/spot

    Returns:
        Dict with spatial metrics
    """
    # Transpose pred to match GT format (genes x spots)
    pred_t = pred.T

    # Align genes and spots
    common_genes = pred_t.index.intersection(gt.index)
    common_spots = pred_t.columns.intersection(gt.columns)

    if len(common_genes) < 10 or len(common_spots) < 10:
        return {
            "overall_pearson": np.nan,
            "per_spot_pearson_mean": np.nan,
            "per_spot_pearson_median": np.nan,
            "per_gene_pearson_mean": np.nan,
            "per_gene_pearson_median": np.nan,
            "variance_ratio": np.nan,
            "n_genes": len(common_genes),
            "n_spots": len(common_spots),
        }

    pred_aligned = pred_t.loc[common_genes, common_spots].values
    gt_aligned = gt.loc[common_genes, common_spots].values

    # Normalize per-spot to CPM, then log1p
    pred_norm = np.log1p(normalize_cpm(pred_aligned, axis=0))
    gt_norm = np.log1p(normalize_cpm(gt_aligned, axis=0))

    n_genes, n_spots = gt_norm.shape

    # 1. Overall Pearson (flattened - baseline)
    overall_r, _ = stats.pearsonr(pred_norm.flatten(), gt_norm.flatten())

    # 2. Per-spot correlation (across genes) - measures gene profile accuracy per location
    spot_correlations = []
    for s in range(n_spots):
        pred_spot = pred_norm[:, s]
        gt_spot = gt_norm[:, s]
        # Skip spots with no variation
        if np.std(pred_spot) > 0 and np.std(gt_spot) > 0:
            r, _ = stats.pearsonr(pred_spot, gt_spot)
            if not np.isnan(r):
                spot_correlations.append(r)

    per_spot_mean = np.mean(spot_correlations) if spot_correlations else np.nan
    per_spot_median = np.median(spot_correlations) if spot_correlations else np.nan

    # 3. Per-gene correlation (across spots) - measures spatial pattern preservation
    gene_correlations = []
    for g in range(n_genes):
        pred_gene = pred_norm[g, :]
        gt_gene = gt_norm[g, :]
        # Skip genes with no spatial variation in GT
        if np.std(gt_gene) > 0 and np.std(pred_gene) > 0:
            r, _ = stats.pearsonr(pred_gene, gt_gene)
            if not np.isnan(r):
                gene_correlations.append(r)

    per_gene_mean = np.mean(gene_correlations) if gene_correlations else np.nan
    per_gene_median = np.median(gene_correlations) if gene_correlations else np.nan

    # 4. Variance ratio - does method preserve heterogeneity?
    # Compute variance across spots for each gene, then average
    pred_var = np.var(pred_norm, axis=1).mean()  # Mean variance across genes
    gt_var = np.var(gt_norm, axis=1).mean()
    variance_ratio = pred_var / gt_var if gt_var > 0 else np.nan

    # 5. Per-gene variance correlation - does high-variance genes in GT also have high variance in pred?
    pred_gene_vars = np.var(pred_norm, axis=1)
    gt_gene_vars = np.var(gt_norm, axis=1)
    if np.std(pred_gene_vars) > 0 and np.std(gt_gene_vars) > 0:
        var_corr, _ = stats.pearsonr(pred_gene_vars, gt_gene_vars)
    else:
        var_corr = np.nan

    return {
        "overall_pearson": overall_r,
        "per_spot_pearson_mean": per_spot_mean,
        "per_spot_pearson_median": per_spot_median,
        "per_gene_pearson_mean": per_gene_mean,
        "per_gene_pearson_median": per_gene_median,
        "variance_ratio": variance_ratio,
        "variance_pattern_corr": var_corr,
        "n_genes": len(common_genes),
        "n_spots": len(common_spots),
        "n_spot_correlations": len(spot_correlations),
        "n_gene_correlations": len(gene_correlations),
    }


def compute_marker_gene_metrics(
    pred: pd.DataFrame,
    gt: pd.DataFrame,
    cell_type: str,
) -> Dict[str, any]:
    """
    Compute metrics specifically for marker genes of a cell type.

    Args:
        pred: Predicted GEX (spots x genes)
        gt: Ground truth GEX (genes x spots)
        cell_type: Cell type name to look up markers

    Returns:
        Dict with marker gene metrics
    """
    markers = MARKER_GENES.get(cell_type, [])
    if not markers:
        return {}

    # Transpose pred to match GT format (genes x spots)
    pred_t = pred.T

    # Find markers present in both pred and GT
    common_spots = pred_t.columns.intersection(gt.columns)
    if len(common_spots) < 10:
        return {}

    marker_results = []
    for marker in markers:
        # Case-insensitive gene matching
        pred_gene = None
        gt_gene = None

        for g in pred_t.index:
            if g.upper() == marker.upper():
                pred_gene = g
                break

        for g in gt.index:
            if g.upper() == marker.upper():
                gt_gene = g
                break

        if pred_gene is None or gt_gene is None:
            continue

        pred_vals = pred_t.loc[pred_gene, common_spots].values.astype(float)
        gt_vals = gt.loc[gt_gene, common_spots].values.astype(float)

        # Normalize to CPM per spot, then log1p
        # For single gene, just log1p the raw values
        pred_log = np.log1p(pred_vals)
        gt_log = np.log1p(gt_vals)

        # Skip if no variation
        if np.std(pred_log) == 0 or np.std(gt_log) == 0:
            marker_results.append({
                "marker": marker,
                "pearson_r": np.nan,
                "pred_mean": pred_vals.mean(),
                "gt_mean": gt_vals.mean(),
                "pred_std": pred_vals.std(),
                "gt_std": gt_vals.std(),
            })
            continue

        r, _ = stats.pearsonr(pred_log, gt_log)

        marker_results.append({
            "marker": marker,
            "pearson_r": r,
            "pred_mean": pred_vals.mean(),
            "gt_mean": gt_vals.mean(),
            "pred_std": pred_vals.std(),
            "gt_std": gt_vals.std(),
        })

    if not marker_results:
        return {}

    # Aggregate metrics
    valid_rs = [m["pearson_r"] for m in marker_results if not np.isnan(m["pearson_r"])]

    return {
        "marker_genes_found": len(marker_results),
        "marker_pearson_mean": np.mean(valid_rs) if valid_rs else np.nan,
        "marker_pearson_median": np.median(valid_rs) if valid_rs else np.nan,
        "marker_details": marker_results,
    }


def evaluate_method(output_subdir: str, method_name: str) -> pd.DataFrame:
    """Evaluate a method across all regions and cell types."""
    results = []

    for region_id in range(5):
        pred_layers = load_citegeist_gex(region_id, output_subdir)

        if not pred_layers:
            logger.warning(f"No predictions for region {region_id}")
            continue

        for ct in CELL_TYPES:
            gt = load_gt_gex(region_id, ct)
            if gt is None:
                continue

            # Try to find matching cell type in predictions
            pred = None
            # Direct match first
            if ct in pred_layers:
                pred = pred_layers[ct]
            else:
                # Fuzzy match: normalize both to compare
                ct_norm = ct.lower().replace("_", "").replace("+", "pos")
                for pred_ct, pred_df in pred_layers.items():
                    pred_ct_norm = pred_ct.lower().replace("_", "").replace("+", "pos")
                    if ct_norm == pred_ct_norm or ct_norm in pred_ct_norm or pred_ct_norm in ct_norm:
                        pred = pred_df
                        break

            if pred is None:
                continue

            metrics = compute_spatial_metrics(pred, gt)
            metrics["region_id"] = region_id
            metrics["cell_type"] = ct
            metrics["method"] = method_name

            # Add marker gene metrics
            marker_metrics = compute_marker_gene_metrics(pred, gt, ct)
            if marker_metrics:
                metrics["marker_genes_found"] = marker_metrics["marker_genes_found"]
                metrics["marker_pearson_mean"] = marker_metrics["marker_pearson_mean"]
                metrics["marker_pearson_median"] = marker_metrics["marker_pearson_median"]
                # Store detailed marker results for later analysis
                metrics["_marker_details"] = marker_metrics["marker_details"]

            results.append(metrics)

    return pd.DataFrame(results)


def main():
    """Run spatial GEX evaluation."""
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-subdir", type=str, required=True,
                        help="Subdirectory under CITEgeist/output/ to evaluate")
    parser.add_argument("--method-name", type=str, default="CITEgeist",
                        help="Name for this method variant")
    args = parser.parse_args()

    logger.info(f"Evaluating: {args.output_subdir}")
    logger.info("=" * 70)

    results = evaluate_method(args.output_subdir, args.method_name)

    if results.empty:
        logger.error("No results - check output directory")
        return

    # Summary statistics
    logger.info(f"\n{'='*70}")
    logger.info(f"SPATIAL GEX EVALUATION: {args.method_name}")
    logger.info(f"{'='*70}")

    # Overall metrics
    logger.info(f"\nOverall (N={len(results)} cell-type × region combinations):")
    for metric in ["overall_pearson", "per_spot_pearson_mean", "per_gene_pearson_mean", "variance_ratio"]:
        mean_val = results[metric].mean()
        std_val = results[metric].std()
        logger.info(f"  {metric}: {mean_val:.4f} ± {std_val:.4f}")

    # Per-cell-type breakdown
    logger.info(f"\nPer-cell-type per_gene_pearson_mean (spatial pattern preservation):")
    ct_means = results.groupby("cell_type")["per_gene_pearson_mean"].mean().sort_values(ascending=False)
    for ct, val in ct_means.items():
        logger.info(f"  {ct}: {val:.4f}")

    # Variance ratio interpretation
    mean_var_ratio = results["variance_ratio"].mean()
    logger.info(f"\nVariance ratio interpretation: {mean_var_ratio:.3f}")
    if mean_var_ratio < 0.5:
        logger.info("  ⚠️ SEVERE: Method is crushing spatial heterogeneity (variance < 50% of GT)")
    elif mean_var_ratio < 0.8:
        logger.info("  ⚠️ MODERATE: Some loss of spatial heterogeneity")
    elif mean_var_ratio > 1.2:
        logger.info("  ⚠️ Method is adding spurious variance")
    else:
        logger.info("  ✓ Good: Variance preserved (~80-120% of GT)")

    # Marker gene analysis
    logger.info(f"\n{'='*70}")
    logger.info("MARKER GENE ANALYSIS")
    logger.info(f"{'='*70}")

    if "marker_pearson_mean" in results.columns:
        # Overall marker gene performance
        marker_mean = results["marker_pearson_mean"].dropna().mean()
        marker_std = results["marker_pearson_mean"].dropna().std()
        logger.info(f"\nOverall marker gene spatial correlation: {marker_mean:.4f} ± {marker_std:.4f}")

        # Per-cell-type marker gene performance
        logger.info(f"\nPer-cell-type marker gene correlation:")
        ct_marker_means = results.groupby("cell_type")["marker_pearson_mean"].mean().sort_values(ascending=False)
        for ct, val in ct_marker_means.items():
            markers = MARKER_GENES.get(ct, [])
            logger.info(f"  {ct}: {val:.4f} (markers: {', '.join(markers[:3])}...)")

        # Detailed per-marker breakdown
        logger.info(f"\nDetailed marker gene correlations:")
        marker_details_all = []
        for _, row in results.iterrows():
            if "_marker_details" in row and row["_marker_details"]:
                for marker_info in row["_marker_details"]:
                    marker_details_all.append({
                        "region_id": row["region_id"],
                        "cell_type": row["cell_type"],
                        "marker": marker_info["marker"],
                        "pearson_r": marker_info["pearson_r"],
                        "pred_mean": marker_info["pred_mean"],
                        "gt_mean": marker_info["gt_mean"],
                    })

        if marker_details_all:
            marker_df = pd.DataFrame(marker_details_all)
            # Average per marker across regions
            marker_avg = marker_df.groupby(["cell_type", "marker"])["pearson_r"].mean()
            for (ct, marker), r in marker_avg.items():
                if not np.isnan(r):
                    logger.info(f"  {ct} / {marker}: {r:.4f}")

            # Save detailed marker results
            marker_output_path = OUTPUT_DIR / args.output_subdir / "marker_gene_details.csv"
            marker_df.to_csv(marker_output_path, index=False)
            logger.info(f"\nMarker details saved to: {marker_output_path}")
    else:
        logger.info("No marker gene data available")

    # Save results (drop internal marker details column)
    results_to_save = results.drop(columns=["_marker_details"], errors="ignore")
    output_path = OUTPUT_DIR / args.output_subdir / "spatial_gex_evaluation.csv"
    results_to_save.to_csv(output_path, index=False)
    logger.info(f"\nResults saved to: {output_path}")

    # Also save summary JSON
    summary = {
        "method": args.method_name,
        "output_subdir": args.output_subdir,
        "n_evaluations": len(results),
        "metrics": {
            "overall_pearson_mean": float(results["overall_pearson"].mean()),
            "per_spot_pearson_mean": float(results["per_spot_pearson_mean"].mean()),
            "per_gene_pearson_mean": float(results["per_gene_pearson_mean"].mean()),
            "variance_ratio_mean": float(results["variance_ratio"].mean()),
        },
        "per_celltype": ct_means.to_dict(),
    }

    # Add marker gene metrics to summary
    if "marker_pearson_mean" in results.columns:
        summary["marker_metrics"] = {
            "marker_pearson_mean": float(results["marker_pearson_mean"].dropna().mean()),
            "per_celltype_marker": ct_marker_means.to_dict() if "ct_marker_means" in dir() else {},
        }

    summary_path = OUTPUT_DIR / args.output_subdir / "spatial_gex_summary.json"
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)
    logger.info(f"Summary saved to: {summary_path}")


if __name__ == "__main__":
    main()
