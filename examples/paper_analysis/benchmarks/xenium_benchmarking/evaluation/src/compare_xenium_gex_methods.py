#!/usr/bin/env python
"""Cross-method GEX evaluation on Xenium pseudo-Visium (SingleR GT, 7 types).

Evaluates CITEgeist-SACE, Cell2Location, Tangram, and (optionally) scResolve
against the SingleR-based GEX ground truth.

Orientation detection:
- Checks whether row index values look like spot names (spot_*)
- If not, transposes (handles scResolve genes×spots format)

Cell type name mapping:
- GT files: B_cells_GT.csv, CD4pos_T_cells_GT.csv, etc.
- Prediction files: may use "B cells_layer.csv" or "B_cells_layer.csv"

Usage:
    python compare_xenium_gex_methods.py [--output_dir /path/]
    python compare_xenium_gex_methods.py --regions 0 1 2 3 4 --methods CITEgeist_SACE Cell2Location
"""

import argparse
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr

REPO_ROOT = Path("/path/to/CITEgeist_analysis")
sys.path.insert(0, str(REPO_ROOT))

BENCH_ROOT = REPO_ROOT / "benchmarks" / "xenium_benchmarking"

GT_ROOT = BENCH_ROOT / "ground_truth_singler" / "ground_truth_7type" / "ground_truth_gex"

ALL_METHODS = {
    "CITEgeist_SACE": BENCH_ROOT / "CITEgeist" / "results_sace_singler_7type",
    "Cell2Location": BENCH_ROOT / "Cell2Location" / "output_singler_7type",
    "Tangram": BENCH_ROOT / "Tangram" / "output_singler_7type",
    "scResolve": BENCH_ROOT / "scResolve" / "output_protein_gated",  # legacy
}

# GT filename stems (used as row keys)
CELL_TYPES = [
    "B_cells",
    "CD4pos_T_cells",
    "CD8pos_T_cells",
    "Endothelial",
    "Epithelial",
    "Fibroblasts",
    "Macrophages",
]

# GT stem → candidate prediction filename stems (methods use different naming conventions)
GT_TO_PRED_NAMES = {
    "B_cells": ["B cells", "B_cells"],
    "CD4pos_T_cells": ["CD4+ T cells", "CD4+_T_cells", "CD4pos_T_cells"],
    "CD8pos_T_cells": ["CD8+ T cells", "CD8+_T_cells", "CD8pos_T_cells"],
    "Endothelial": ["Endothelial"],
    "Epithelial": ["Epithelial"],
    "Fibroblasts": ["Fibroblasts"],
    "Macrophages": ["Macrophages"],
}

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger(__name__)


def load_csv_spots_genes(path: Path) -> pd.DataFrame:
    """Load a CSV and orient to spots x genes.

    Auto-detects orientation: if row index looks like spot_* names -> already correct.
    Otherwise transposes (genes x spots -> spots x genes).

    Args:
        path: path to CSV file

    Returns:
        DataFrame with spots as rows, genes as columns
    """
    df = pd.read_csv(path, index_col=0)
    # Check if index looks like spot names
    idx_sample = str(df.index[0]) if len(df.index) > 0 else ""
    if idx_sample.startswith("spot_") or idx_sample.startswith("AAAC"):
        return df  # already spots x genes
    # Check if columns look like spot names
    if len(df.columns) > 0 and str(df.columns[0]).startswith("spot_"):
        return df.T  # genes x spots -> transpose
    # Heuristic: if more rows than columns, likely genes x spots (Xenium has ~405 genes)
    if df.shape[0] > df.shape[1]:
        return df.T
    return df


def compute_metrics(pred: np.ndarray, gt: np.ndarray) -> dict:
    """Compute GEX evaluation metrics on log1p-transformed values.

    Args:
        pred: predicted values (spots x genes), raw counts
        gt: ground truth values (spots x genes), raw counts

    Returns:
        dict with flattened_r, spot_layer_r, rmse, nrmse, mae
    """
    pred = pred.astype(np.float64)
    gt = gt.astype(np.float64)

    # Rigor fix: uniform per-block GT-scale normalization for scale-fair RMSE.
    # Methods differ in their reference scale (e.g. Tangram uses raw counts 10-40x
    # higher than GT because its reference was not library-size normalised).
    # Pearson r is scale-invariant and unaffected; RMSE would otherwise be
    # confounded by scale differences across methods.
    # We rescale each (region, method, cell_type) prediction block to match
    # the GT block's total mass with a single global scalar, leaving GT unchanged.
    s_pred = float(pred.sum())
    s_gt = float(gt.sum())
    if s_pred > 0:
        pred = pred * (s_gt / s_pred)

    pred_log = np.log1p(pred)
    gt_log = np.log1p(gt)

    pv = pred_log.ravel()
    gv = gt_log.ravel()
    mask = (pv > 0) | (gv > 0)

    if mask.sum() < 10:
        return {"flattened_r": np.nan, "spot_layer_r": np.nan, "rmse": np.nan, "nrmse": np.nan, "mae": np.nan}

    # Flattened Pearson r (both-nonzero mask to avoid zero-inflation bias)
    flat_r = float(pearsonr(pv[mask], gv[mask])[0]) if mask.sum() > 1 else np.nan

    # Spot-layer r: per-gene Pearson r across spots, averaged
    per_gene_rs = []
    for g in range(pred_log.shape[1]):
        pg = pred_log[:, g]
        gg = gt_log[:, g]
        gmask = (pg > 0) | (gg > 0)
        if gmask.sum() < 3:
            continue
        r_g, _ = pearsonr(pg[gmask], gg[gmask])
        if np.isfinite(r_g):
            per_gene_rs.append(r_g)
    spot_layer_r = float(np.mean(per_gene_rs)) if per_gene_rs else np.nan

    # Error metrics on full (non-masked) log1p values
    diff = pred_log - gt_log
    rmse = float(np.sqrt(np.mean(diff**2)))
    gt_range = gt_log.max() - gt_log.min()
    nrmse = rmse / gt_range if gt_range > 0 else np.nan
    mae = float(np.mean(np.abs(diff)))

    return {
        "flattened_r": flat_r,
        "spot_layer_r": spot_layer_r,
        "rmse": rmse,
        "nrmse": nrmse,
        "mae": mae,
    }


def find_layer_file(layers_dir: Path, gt_stem: str) -> Path | None:
    """Find the prediction layer CSV for a given GT cell type stem.

    Tries all candidate names from GT_TO_PRED_NAMES, both with _layer.csv suffix.

    Args:
        layers_dir: directory containing {TypeName}_layer.csv files
        gt_stem: GT filename stem (e.g., "CD4pos_T_cells")

    Returns:
        Path if found, else None
    """
    for candidate in GT_TO_PRED_NAMES.get(gt_stem, [gt_stem]):
        path = layers_dir / f"{candidate}_layer.csv"
        if path.exists():
            return path
    return None


def evaluate_method_on_region(method_name: str, method_base_dir: Path, region_id: int) -> list:
    """Evaluate one method on one region.

    Args:
        method_name: display name for the method
        method_base_dir: base directory for the method's results
        region_id: 0-4

    Returns:
        List of result dicts, one per cell type.
    """
    layers_dir = method_base_dir / f"Xenium_region_{region_id}" / "layers"
    if not layers_dir.exists():
        logger.warning("  %s region %d: layers dir missing (%s)", method_name, region_id, layers_dir)
        return []

    results = []
    for gt_stem in CELL_TYPES:
        gt_path = GT_ROOT / f"Xenium_region_{region_id}" / f"{gt_stem}_GT.csv"
        if not gt_path.exists():
            logger.warning("  GT missing: %s", gt_path)
            continue

        layer_path = find_layer_file(layers_dir, gt_stem)
        if layer_path is None:
            logger.warning(
                "  %s region %d / %s: no layer file found in %s", method_name, region_id, gt_stem, layers_dir
            )
            continue

        # Load and orient to spots x genes
        gt_df = load_csv_spots_genes(gt_path)
        pred_df = load_csv_spots_genes(layer_path)

        # Normalize gene names to uppercase (avoid case mismatches and Excel date corruption)
        gt_df.columns = gt_df.columns.str.upper()
        pred_df.columns = pred_df.columns.str.upper()
        gt_df = gt_df.loc[:, ~gt_df.columns.duplicated()]
        pred_df = pred_df.loc[:, ~pred_df.columns.duplicated()]

        # Align on shared spots and genes
        common_spots = sorted(set(gt_df.index) & set(pred_df.index))
        common_genes = sorted(set(gt_df.columns) & set(pred_df.columns))
        n_shared_genes = len(common_genes)

        if len(common_spots) < 5 or n_shared_genes < 5:
            logger.warning(
                "  %s region %d / %s: insufficient overlap (%d spots, %d genes)",
                method_name,
                region_id,
                gt_stem,
                len(common_spots),
                n_shared_genes,
            )
            continue

        gt_vals = gt_df.loc[common_spots, common_genes].values.astype(np.float64)
        pred_vals = pred_df.loc[common_spots, common_genes].values.astype(np.float64)

        # Filter to GT-expressed genes (nonzero total before log1p)
        gt_gene_totals = gt_vals.sum(axis=0)  # sum across spots per gene
        expressed_mask = gt_gene_totals > 0
        n_total_genes = len(expressed_mask)
        n_expressed = int(expressed_mask.sum())
        gt_vals = gt_vals[:, expressed_mask]
        pred_vals = pred_vals[:, expressed_mask]
        common_genes = [g for g, m in zip(common_genes, expressed_mask) if m]
        n_shared_genes = n_expressed
        logger.info("  Expressed gene filter: %d/%d genes retained for %s", n_expressed, n_total_genes, gt_stem)

        if n_expressed < 5:
            logger.warning(
                "  %s region %d / %s: insufficient expressed genes (%d)",
                method_name,
                region_id,
                gt_stem,
                n_expressed,
            )
            continue

        metrics = compute_metrics(pred_vals, gt_vals)
        results.append(
            {
                "method": method_name,
                "region_id": region_id,
                "cell_type": gt_stem,
                "n_shared_spots": len(common_spots),
                "n_shared_genes": n_shared_genes,
                **metrics,
            }
        )
        logger.info(
            "  %s region %d / %-20s: flat_r=%.4f, spot_layer_r=%.4f, n_genes=%d",
            method_name,
            region_id,
            gt_stem,
            metrics["flattened_r"] if metrics["flattened_r"] is not None else -1,
            metrics["spot_layer_r"] if metrics["spot_layer_r"] is not None else -1,
            n_shared_genes,
        )

    return results


def main():
    parser = argparse.ArgumentParser(
        description="Cross-method GEX evaluation on Xenium pseudo-Visium (SingleR GT, 7 types)"
    )
    parser.add_argument(
        "--regions", type=int, nargs="+", default=list(range(5)), help="Region IDs to evaluate (default: 0 1 2 3 4)"
    )
    parser.add_argument(
        "--methods",
        type=str,
        nargs="+",
        default=list(ALL_METHODS.keys()),
        choices=list(ALL_METHODS.keys()),
        help="Methods to evaluate (default: all)",
    )
    parser.add_argument(
        "--output_dir", type=str, default=str(BENCH_ROOT / "evaluation"), help="Output directory for comparison CSVs"
    )
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Filter to requested methods
    methods_to_run = {k: v for k, v in ALL_METHODS.items() if k in args.methods}

    all_results = []

    for method_name, method_base_dir in methods_to_run.items():
        if not method_base_dir.exists():
            logger.warning("Method dir missing, skipping: %s (%s)", method_name, method_base_dir)
            continue

        logger.info("=" * 60)
        logger.info("Evaluating: %s", method_name)
        logger.info("=" * 60)

        for region_id in args.regions:
            region_results = evaluate_method_on_region(method_name, method_base_dir, region_id)
            all_results.extend(region_results)

    if not all_results:
        logger.error("No results produced. Check input paths.")
        return

    # Save per-type per-region CSV
    df = pd.DataFrame(all_results)
    out_csv = output_dir / "xenium_gex_comparison_7type.csv"
    df.to_csv(out_csv, index=False)
    logger.info("\nSaved per-region results to %s", out_csv)

    # Summary: mean across regions
    metric_cols = ["flattened_r", "spot_layer_r", "rmse", "nrmse", "mae"]
    summary_rows = []
    for method in df["method"].unique():
        for ct in df["cell_type"].unique():
            subset = df[(df["method"] == method) & (df["cell_type"] == ct)]
            if subset.empty:
                continue
            row = {"method": method, "cell_type": ct}
            for col in metric_cols:
                vals = subset[col].dropna()
                row[f"{col}_mean"] = round(float(vals.mean()), 4) if len(vals) > 0 else np.nan
                row[f"{col}_std"] = round(float(vals.std()), 4) if len(vals) > 1 else np.nan
            row["n_regions"] = len(subset)
            row["n_shared_genes_mean"] = round(float(subset["n_shared_genes"].mean()), 0)
            summary_rows.append(row)

    # Overall per method (across all types and regions)
    for method in df["method"].unique():
        subset = df[df["method"] == method]
        row = {"method": method, "cell_type": "ALL"}
        for col in metric_cols:
            vals = subset[col].dropna()
            row[f"{col}_mean"] = round(float(vals.mean()), 4) if len(vals) > 0 else np.nan
            row[f"{col}_std"] = round(float(vals.std()), 4) if len(vals) > 1 else np.nan
        row["n_regions"] = len(subset["region_id"].unique())
        row["n_shared_genes_mean"] = round(float(subset["n_shared_genes"].mean()), 0)
        summary_rows.append(row)

    summary_df = pd.DataFrame(summary_rows)
    out_summary = output_dir / "xenium_gex_comparison_7type_summary.csv"
    summary_df.to_csv(out_summary, index=False)
    logger.info("Saved summary to %s", out_summary)

    # Print summary table
    logger.info("\n%s", "=" * 70)
    logger.info("XENIUM GEX CROSS-METHOD SUMMARY (SingleR GT, 7 types)")
    logger.info("%s", "=" * 70)
    primary_methods = ["CITEgeist_SACE", "Cell2Location", "Tangram"]
    overall = summary_df[summary_df["cell_type"] == "ALL"]
    logger.info("%-20s %12s %12s %12s", "Method", "flat_r_mean", "spot_lr_mean", "n_genes_mean")
    logger.info("-" * 58)
    for method in primary_methods:
        row = overall[overall["method"] == method]
        if row.empty:
            continue
        logger.info(
            "%-20s %12.4f %12.4f %12.0f",
            method,
            row["flattened_r_mean"].values[0],
            row["spot_layer_r_mean"].values[0],
            row["n_shared_genes_mean"].values[0],
        )

    # scResolve as legacy
    scresolve_row = overall[overall["method"] == "scResolve"]
    if not scresolve_row.empty:
        logger.info(
            "%-20s %12.4f %12.4f %12.0f  [legacy]",
            "scResolve",
            scresolve_row["flattened_r_mean"].values[0],
            scresolve_row["spot_layer_r_mean"].values[0],
            scresolve_row["n_shared_genes_mean"].values[0],
        )


if __name__ == "__main__":
    main()
