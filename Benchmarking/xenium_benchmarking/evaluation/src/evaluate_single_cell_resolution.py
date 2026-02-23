"""
Evaluate single-cell resolution benchmark results.

This script:
1. Loads single-cell AnnData from Module 3b pipeline
2. Aggregates cell counts back to spot-level proportions
3. Compares against protein-gated ground truth proportions
4. Aggregates GEX back to spot level
5. Compares against ground truth GEX

Uses the same metrics as other benchmark evaluations for fair comparison.
"""

import argparse
import json
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats
from scipy.spatial.distance import jensenshannon

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Same cell types as other benchmarks
ACHIEVABLE_7_CELL_TYPES = [
    "B cells",
    "CD4+ T cells",
    "CD8+ T cells",
    "Macrophages",
    "Endothelial",
    "Epithelial",
    "Fibroblasts",
]


def load_single_cell_adata(region_id: int, sc_dir: str) -> Optional[sc.AnnData]:
    """Load single-cell AnnData for a region."""
    sc_path = Path(sc_dir) / f"Xenium_region_{region_id}" / f"Xenium_region_{region_id}_single_cell.h5ad"
    if not sc_path.exists():
        logger.warning(f"Single-cell AnnData not found: {sc_path}")
        return None
    return sc.read_h5ad(sc_path)


def aggregate_to_spot_proportions(adata: sc.AnnData, cell_types: List[str]) -> pd.DataFrame:
    """
    Aggregate single-cell assignments to spot-level proportions.

    Args:
        adata: Single-cell AnnData with 'spot_id' and 'cell_type' in obs
        cell_types: List of cell types to include

    Returns:
        DataFrame with spots as rows, cell types as columns, proportions as values
    """
    # Count cells per spot per type
    counts = adata.obs.groupby(['spot_id', 'cell_type']).size().unstack(fill_value=0)

    # Ensure all cell types are present
    for ct in cell_types:
        if ct not in counts.columns:
            counts[ct] = 0
    counts = counts[cell_types]

    # Convert to proportions
    row_sums = counts.sum(axis=1)
    proportions = counts.div(row_sums, axis=0).fillna(0)

    return proportions


def aggregate_to_spot_gex(adata: sc.AnnData) -> pd.DataFrame:
    """
    Aggregate single-cell GEX to spot level (total across all cell types).

    Args:
        adata: Single-cell AnnData with 'spot_id' in obs

    Returns:
        DataFrame with spots as rows, genes as columns
    """
    # Get expression matrix
    if hasattr(adata.X, 'toarray'):
        X = adata.X.toarray()
    else:
        X = adata.X

    # Create DataFrame
    gex_df = pd.DataFrame(X, index=adata.obs.index, columns=adata.var_names)
    gex_df['spot_id'] = adata.obs['spot_id'].values

    # Sum expression per spot
    spot_gex = gex_df.groupby('spot_id').sum()

    return spot_gex


def aggregate_to_spot_gex_by_celltype(adata: sc.AnnData, cell_types: List[str]) -> Dict[str, pd.DataFrame]:
    """
    Aggregate single-cell GEX to spot level, per cell type.

    Args:
        adata: Single-cell AnnData with 'spot_id' and 'cell_type' in obs
        cell_types: List of cell types

    Returns:
        Dict mapping cell_type -> DataFrame (spots x genes)
    """
    # Get expression matrix
    if hasattr(adata.X, 'toarray'):
        X = adata.X.toarray()
    else:
        X = adata.X

    result = {}
    for ct in cell_types:
        # Filter to cells of this type
        mask = adata.obs['cell_type'] == ct
        if mask.sum() == 0:
            continue

        # Create DataFrame for this cell type
        gex_df = pd.DataFrame(X[mask], index=adata.obs.index[mask], columns=adata.var_names)
        gex_df['spot_id'] = adata.obs.loc[mask, 'spot_id'].values

        # Sum expression per spot
        spot_gex = gex_df.groupby('spot_id').sum()
        result[ct] = spot_gex

    return result


def load_ground_truth_proportions(region_id: int, gt_dir: str) -> pd.DataFrame:
    """Load ground truth proportions for a region."""
    gt_path = Path(gt_dir) / "ground_truth" / f"Xenium_region_{region_id}_prop.csv"
    if not gt_path.exists():
        raise FileNotFoundError(f"Ground truth not found: {gt_path}")

    gt_df = pd.read_csv(gt_path, index_col=0)
    return gt_df


def load_ground_truth_gex_by_celltype(region_id: int, gt_dir: str) -> Dict[str, pd.DataFrame]:
    """
    Load ground truth GEX for a region, per cell type.

    Ground truth format: {cell_type}_GT.csv files with genes as rows, spots as columns.

    Returns:
        Dict mapping cell_type -> DataFrame (spots x genes)
    """
    gt_gex_dir = Path(gt_dir) / "ground_truth_gex" / f"Xenium_region_{region_id}"
    if not gt_gex_dir.exists():
        raise FileNotFoundError(f"Ground truth GEX dir not found: {gt_gex_dir}")

    # Map file names to standardized cell type names
    celltype_file_map = {
        "B cells": "B_cells_GT.csv",
        "CD4+ T cells": "CD4pos_T_cells_GT.csv",
        "CD8+ T cells": "CD8pos_T_cells_GT.csv",
        "Macrophages": "Macrophages_GT.csv",
        "Endothelial": "Endothelial_GT.csv",
        "Epithelial": "Epithelial_GT.csv",
        "Fibroblasts": "Fibroblasts_GT.csv",
    }

    result = {}
    for ct, filename in celltype_file_map.items():
        filepath = gt_gex_dir / filename
        if filepath.exists():
            # Load and transpose (genes x spots -> spots x genes)
            df = pd.read_csv(filepath, index_col=0).T
            result[ct] = df

    return result


def load_ground_truth_gex_total(region_id: int, gt_dir: str) -> pd.DataFrame:
    """Load total ground truth GEX (summed across all cell types)."""
    ct_gex = load_ground_truth_gex_by_celltype(region_id, gt_dir)

    if not ct_gex:
        raise FileNotFoundError(f"No ground truth GEX found for region {region_id}")

    # Sum across cell types
    total = None
    for ct, df in ct_gex.items():
        if total is None:
            total = df.copy()
        else:
            # Align indices before adding
            common_spots = total.index.intersection(df.index)
            common_genes = total.columns.intersection(df.columns)
            total = total.loc[common_spots, common_genes] + df.loc[common_spots, common_genes]

    return total


def calculate_proportion_metrics(
    gt_df: pd.DataFrame,
    pred_df: pd.DataFrame,
    cell_types: List[str],
) -> Dict[str, float]:
    """Calculate proportion evaluation metrics."""
    metrics = {}

    # Align indices
    common_spots = gt_df.index.intersection(pred_df.index)
    if len(common_spots) == 0:
        logger.warning("No common spots between GT and predictions")
        return {"error": "no_common_spots"}

    gt_aligned = gt_df.loc[common_spots]
    pred_aligned = pred_df.loc[common_spots]

    # Filter to common cell types
    common_types = [ct for ct in cell_types if ct in gt_aligned.columns and ct in pred_aligned.columns]

    gt_matrix = gt_aligned[common_types].values
    pred_matrix = pred_aligned[common_types].values

    # Per-cell-type metrics
    for i, ct in enumerate(common_types):
        gt_vals = gt_matrix[:, i]
        pred_vals = pred_matrix[:, i]

        r, p = stats.pearsonr(gt_vals, pred_vals)
        metrics[f"{ct}_pearson_r"] = r
        metrics[f"{ct}_rmse"] = np.sqrt(np.mean((gt_vals - pred_vals) ** 2))
        metrics[f"{ct}_mae"] = np.mean(np.abs(gt_vals - pred_vals))

    # Overall metrics (flattened)
    gt_flat = gt_matrix.flatten()
    pred_flat = pred_matrix.flatten()

    r_overall, _ = stats.pearsonr(gt_flat, pred_flat)
    metrics["overall_pearson_r"] = r_overall
    metrics["overall_rmse"] = np.sqrt(np.mean((gt_flat - pred_flat) ** 2))
    metrics["overall_mae"] = np.mean(np.abs(gt_flat - pred_flat))

    # Per-spot JSD
    jsd_values = []
    for i in range(len(gt_aligned)):
        gt_row = gt_matrix[i] + 1e-10
        pred_row = pred_matrix[i] + 1e-10
        gt_row = gt_row / gt_row.sum()
        pred_row = pred_row / pred_row.sum()
        jsd = jensenshannon(gt_row, pred_row)
        if not np.isnan(jsd):
            jsd_values.append(jsd)

    metrics["mean_jsd"] = np.mean(jsd_values) if jsd_values else np.nan
    metrics["n_spots"] = len(common_spots)

    return metrics


def calculate_gex_metrics(
    gt_gex: pd.DataFrame,
    pred_gex: pd.DataFrame,
) -> Dict[str, float]:
    """Calculate GEX evaluation metrics."""
    metrics = {}

    # Align spots and genes
    common_spots = gt_gex.index.intersection(pred_gex.index)
    common_genes = gt_gex.columns.intersection(pred_gex.columns)

    if len(common_spots) == 0 or len(common_genes) == 0:
        logger.warning(f"No common spots ({len(common_spots)}) or genes ({len(common_genes)})")
        return {"error": "no_common_data"}

    gt_aligned = gt_gex.loc[common_spots, common_genes]
    pred_aligned = pred_gex.loc[common_spots, common_genes]

    # Normalize to CPM for comparison
    gt_cpm = gt_aligned.div(gt_aligned.sum(axis=1), axis=0) * 1e6
    pred_cpm = pred_aligned.div(pred_aligned.sum(axis=1), axis=0) * 1e6

    # Fill NaNs from zero-sum rows
    gt_cpm = gt_cpm.fillna(0)
    pred_cpm = pred_cpm.fillna(0)

    # Overall correlation (flattened)
    gt_flat = gt_cpm.values.flatten()
    pred_flat = pred_cpm.values.flatten()

    # Remove zeros for correlation
    mask = (gt_flat > 0) | (pred_flat > 0)
    if mask.sum() > 10:
        r_overall, _ = stats.pearsonr(gt_flat[mask], pred_flat[mask])
        metrics["gex_pearson_r"] = r_overall
    else:
        metrics["gex_pearson_r"] = np.nan

    # RMSE on log-scale
    gt_log = np.log1p(gt_flat)
    pred_log = np.log1p(pred_flat)
    metrics["gex_log_rmse"] = np.sqrt(np.mean((gt_log - pred_log) ** 2))

    # Per-gene correlation (average)
    gene_corrs = []
    for gene in common_genes:
        gt_gene = gt_cpm[gene].values
        pred_gene = pred_cpm[gene].values
        if gt_gene.std() > 0 and pred_gene.std() > 0:
            r, _ = stats.pearsonr(gt_gene, pred_gene)
            if not np.isnan(r):
                gene_corrs.append(r)

    metrics["mean_gene_pearson_r"] = np.mean(gene_corrs) if gene_corrs else np.nan
    metrics["n_genes_evaluated"] = len(gene_corrs)
    metrics["n_genes_total"] = len(common_genes)
    metrics["n_spots"] = len(common_spots)

    return metrics


def evaluate_region(
    region_id: int,
    sc_dir: str,
    gt_dir: str,
    cell_types: List[str] = ACHIEVABLE_7_CELL_TYPES,
) -> Dict[str, any]:
    """Evaluate a single region."""
    logger.info(f"Evaluating region {region_id}")

    results = {
        "region_id": region_id,
        "status": "success",
    }

    # Load single-cell data
    adata = load_single_cell_adata(region_id, sc_dir)
    if adata is None:
        results["status"] = "failed"
        results["error"] = "single_cell_adata_not_found"
        return results

    logger.info(f"Loaded {adata.n_obs} cells from single-cell AnnData")
    results["n_cells"] = adata.n_obs

    # Aggregate to spot proportions
    pred_props = aggregate_to_spot_proportions(adata, cell_types)
    logger.info(f"Aggregated to {len(pred_props)} spots")

    # Load ground truth proportions
    try:
        gt_props = load_ground_truth_proportions(region_id, gt_dir)
        gt_props = gt_props[[ct for ct in cell_types if ct in gt_props.columns]]

        # Calculate proportion metrics
        prop_metrics = calculate_proportion_metrics(gt_props, pred_props, cell_types)
        results["proportion_metrics"] = prop_metrics
        logger.info(f"Proportion Pearson r: {prop_metrics.get('overall_pearson_r', 'N/A'):.4f}")
    except Exception as e:
        logger.error(f"Failed to evaluate proportions: {e}")
        results["proportion_metrics"] = {"error": str(e)}

    # Aggregate to spot GEX (total)
    try:
        pred_gex = aggregate_to_spot_gex(adata)
        logger.info(f"Aggregated total GEX: {pred_gex.shape}")

        # Load ground truth GEX (total)
        gt_gex = load_ground_truth_gex_total(region_id, gt_dir)

        # Calculate GEX metrics (total)
        gex_metrics = calculate_gex_metrics(gt_gex, pred_gex)
        results["gex_metrics_total"] = gex_metrics
        logger.info(f"Total GEX Pearson r: {gex_metrics.get('gex_pearson_r', 'N/A'):.4f}")
    except Exception as e:
        logger.error(f"Failed to evaluate total GEX: {e}")
        results["gex_metrics_total"] = {"error": str(e)}

    # Aggregate to spot GEX by cell type (more informative)
    try:
        pred_gex_by_ct = aggregate_to_spot_gex_by_celltype(adata, cell_types)
        gt_gex_by_ct = load_ground_truth_gex_by_celltype(region_id, gt_dir)

        # Calculate per-celltype GEX metrics
        ct_gex_metrics = {}
        ct_gex_corrs = []
        for ct in cell_types:
            if ct in pred_gex_by_ct and ct in gt_gex_by_ct:
                metrics = calculate_gex_metrics(gt_gex_by_ct[ct], pred_gex_by_ct[ct])
                ct_gex_metrics[ct] = metrics
                if "gex_pearson_r" in metrics and not np.isnan(metrics["gex_pearson_r"]):
                    ct_gex_corrs.append(metrics["gex_pearson_r"])
                    logger.info(f"  {ct} GEX Pearson r: {metrics['gex_pearson_r']:.4f}")

        results["gex_metrics_by_celltype"] = ct_gex_metrics
        if ct_gex_corrs:
            results["gex_metrics"] = {
                "mean_celltype_pearson_r": np.mean(ct_gex_corrs),
                "std_celltype_pearson_r": np.std(ct_gex_corrs),
                "n_celltypes_evaluated": len(ct_gex_corrs),
            }
            logger.info(f"Mean celltype GEX Pearson r: {np.mean(ct_gex_corrs):.4f}")
    except Exception as e:
        logger.error(f"Failed to evaluate celltype GEX: {e}")
        results["gex_metrics_by_celltype"] = {"error": str(e)}

    return results


def main():
    parser = argparse.ArgumentParser(description="Evaluate single-cell resolution benchmark")
    parser.add_argument("--sc-dir", type=str, required=True,
                       help="Directory containing single-cell resolution output")
    parser.add_argument("--gt-dir", type=str, required=True,
                       help="Directory containing ground truth data")
    parser.add_argument("--regions", type=str, default="0,1,2,4",
                       help="Comma-separated region IDs to evaluate")
    parser.add_argument("--output", type=str, required=True,
                       help="Output JSON file for results")

    args = parser.parse_args()

    regions = [int(r.strip()) for r in args.regions.split(",")]

    all_results = {
        "method": "CITEgeist_SingleCell",
        "regions": {},
        "aggregate": {},
    }

    # Evaluate each region
    prop_r_values = []
    gex_r_values = []

    for region_id in regions:
        results = evaluate_region(region_id, args.sc_dir, args.gt_dir)
        all_results["regions"][f"region_{region_id}"] = results

        if results["status"] == "success":
            if "proportion_metrics" in results and "overall_pearson_r" in results["proportion_metrics"]:
                prop_r_values.append(results["proportion_metrics"]["overall_pearson_r"])
            if "gex_metrics" in results and "mean_celltype_pearson_r" in results["gex_metrics"]:
                gex_r_values.append(results["gex_metrics"]["mean_celltype_pearson_r"])

    # Aggregate metrics
    if prop_r_values:
        all_results["aggregate"]["mean_proportion_pearson_r"] = np.mean(prop_r_values)
        all_results["aggregate"]["std_proportion_pearson_r"] = np.std(prop_r_values)

    if gex_r_values:
        all_results["aggregate"]["mean_gex_pearson_r"] = np.mean(gex_r_values)
        all_results["aggregate"]["std_gex_pearson_r"] = np.std(gex_r_values)

    all_results["aggregate"]["n_regions_evaluated"] = len(regions)
    all_results["aggregate"]["n_regions_successful"] = len(prop_r_values)

    # Save results
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, 'w') as f:
        json.dump(all_results, f, indent=2)

    logger.info(f"Results saved to {output_path}")

    # Print summary
    print("\n" + "="*60)
    print("SINGLE-CELL RESOLUTION BENCHMARK RESULTS")
    print("="*60)

    if prop_r_values:
        print(f"\nProportion Metrics ({len(prop_r_values)} regions):")
        print(f"  Mean Pearson r: {np.mean(prop_r_values):.4f} +/- {np.std(prop_r_values):.4f}")

    if gex_r_values:
        print(f"\nGEX Metrics ({len(gex_r_values)} regions):")
        print(f"  Mean Pearson r: {np.mean(gex_r_values):.4f} +/- {np.std(gex_r_values):.4f}")

    print("\n" + "="*60)


if __name__ == "__main__":
    main()
