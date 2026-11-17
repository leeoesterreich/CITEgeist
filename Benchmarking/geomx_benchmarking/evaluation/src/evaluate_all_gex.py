"""
Evaluate per-cell-type GEX deconvolution on pseudo-GeoMx benchmark.

Only CITEgeist, Cell2Location, and Tangram produce GEX layers.
RCTD, CARD, Seurat are proportions-only.

Metrics: per-type Pearson r on log1p-transformed counts, RMSE, overall flat r.
"""

import json
import logging
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger(__name__)

BASE_DIR = Path(__file__).resolve().parent.parent.parent
GT_DIR = BASE_DIR / "data" / "ground_truth_gex"

# Cell type name → filename prefix mapping
CT_TO_FILENAME = {
    "Epithelial": "Epithelial",
    "Fibroblasts": "Fibroblasts",
    "Endothelial": "Endothelial",
    "CD4+ T cells": "CD4pos_T_cells",
    "CD8+ T cells": "CD8pos_T_cells",
    "B cells": "B_cells",
    "Macrophages": "Macrophages",
}

GEX_METHODS = {
    "CITEgeist": BASE_DIR / "CITEgeist" / "output" / "pseudo_geomx" / "pipeline" / "layers",
    "Cell2Location": BASE_DIR / "Cell2Location" / "output" / "layers",
    "Tangram": BASE_DIR / "Tangram" / "output" / "layers",
}


def load_gt_gex(cell_type):
    """Load per-cell-type GEX ground truth."""
    fname = f"{CT_TO_FILENAME[cell_type]}_GT.csv"
    path = GT_DIR / fname
    if not path.exists():
        return None
    return pd.read_csv(path, index_col=0)


def load_method_gex(method, layers_dir, cell_type):
    """Load per-cell-type GEX prediction from a method's layers directory."""
    if not layers_dir.exists():
        return None

    # Try multiple naming conventions
    candidates = [
        f"{cell_type}_layer.csv",
        f"{CT_TO_FILENAME[cell_type]}_layer.csv",
    ]
    for fname in candidates:
        path = layers_dir / fname
        if path.exists():
            return pd.read_csv(path, index_col=0)
    return None


def evaluate_gex(gt_df, pred_df, cell_type, method):
    """Compute GEX evaluation metrics for one cell type."""
    # Align AOIs and genes
    common_aois = gt_df.index.intersection(pred_df.index)
    common_genes = gt_df.columns.intersection(pred_df.columns)
    if len(common_aois) == 0 or len(common_genes) == 0:
        return None

    gt = gt_df.loc[common_aois, common_genes].values.astype(float)
    pred = pred_df.loc[common_aois, common_genes].values.astype(float)

    # Filter to GT-expressed genes (nonzero total before log1p)
    gt_gene_totals = gt.sum(axis=0)  # sum across AOIs per gene
    expressed_mask = gt_gene_totals > 0
    n_total = len(expressed_mask)
    n_expressed = int(expressed_mask.sum())
    gt = gt[:, expressed_mask]
    pred = pred[:, expressed_mask]
    common_genes = [g for g, m in zip(common_genes, expressed_mask) if m]
    logger.info("  %s/%s: %d/%d genes expressed in GT", method, cell_type, n_expressed, n_total)

    if n_expressed == 0:
        return None

    # Log1p transform
    gt_log = np.log1p(gt)
    pred_log = np.log1p(pred)

    # Flattened Pearson r
    gt_flat = gt_log.flatten()
    pred_flat = pred_log.flatten()
    if np.std(gt_flat) == 0 or np.std(pred_flat) == 0:
        r = 0.0
    else:
        r, _ = stats.pearsonr(gt_flat, pred_flat)

    rmse = np.sqrt(np.mean((gt_log - pred_log) ** 2))

    return {
        "cell_type": cell_type,
        "method": method,
        "n_aois": len(common_aois),
        "n_genes": len(common_genes),
        "flat_r": round(r, 4),
        "rmse": round(rmse, 4),
    }


def main():
    output_dir = BASE_DIR / "evaluation" / "results"
    output_dir.mkdir(parents=True, exist_ok=True)

    all_results = []

    for method, layers_dir in GEX_METHODS.items():
        logger.info(f"\n=== {method} ===")
        if not layers_dir.exists():
            logger.warning(f"Layers dir not found: {layers_dir}")
            continue

        method_rs = []
        for ct in CT_TO_FILENAME.keys():
            gt_df = load_gt_gex(ct)
            if gt_df is None:
                logger.warning(f"  GT not found for {ct}")
                continue

            pred_df = load_method_gex(method, layers_dir, ct)
            if pred_df is None:
                logger.warning(f"  {method} layer not found for {ct}")
                continue

            result = evaluate_gex(gt_df, pred_df, ct, method)
            if result:
                all_results.append(result)
                method_rs.append(result["flat_r"])
                logger.info(f"  {ct}: flat_r={result['flat_r']}, rmse={result['rmse']}")

        if method_rs:
            logger.info(f"  {method} mean flat_r: {np.mean(method_rs):.4f}")

    # Save results
    if all_results:
        df = pd.DataFrame(all_results)
        csv_path = output_dir / "pseudo_geomx_all_methods_gex.csv"
        df.to_csv(csv_path, index=False)
        logger.info(f"\nWrote {csv_path}")

        # Summary per method
        summary = (
            df.groupby("method")
            .agg(
                mean_flat_r=("flat_r", "mean"),
                mean_rmse=("rmse", "mean"),
            )
            .sort_values("mean_flat_r", ascending=False)
        )

        summary_path = output_dir / "pseudo_geomx_all_methods_gex_summary.csv"
        summary.to_csv(summary_path)
        logger.info(f"Wrote {summary_path}")

        print("\n" + "=" * 60)
        print("Pseudo-GeoMx GEX Benchmark Results")
        print("=" * 60)
        print(summary.to_string())
    else:
        logger.warning("No GEX results to report")


if __name__ == "__main__":
    main()
