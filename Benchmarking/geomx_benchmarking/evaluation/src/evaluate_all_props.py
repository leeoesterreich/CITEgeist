"""
Evaluate all methods' cell type proportions on pseudo-GeoMx benchmark.

Metrics per cell type: Pearson r, RMSE, MAE
Overall: flattened Pearson r, mean JSD, RMSE, MAE

Reads predicted_proportions.csv from each method's output/ directory.
"""

import json
import logging
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from scipy.spatial.distance import jensenshannon

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger(__name__)

BASE_DIR = Path(__file__).resolve().parent.parent.parent  # Benchmarking/geomx_benchmarking/
GT_PATH = BASE_DIR / "data" / "aoi_gt_proportions.csv"

METHODS = {
    "CITEgeist": BASE_DIR / "CITEgeist" / "output" / "pseudo_geomx" / "pipeline" / "predicted_proportions.csv",
    "Cell2Location": BASE_DIR / "Cell2Location" / "output" / "predicted_proportions.csv",
    "RCTD": BASE_DIR / "RCTD" / "output" / "predicted_proportions.csv",
    "CARD": BASE_DIR / "CARD" / "output" / "predicted_proportions.csv",
    "Tangram": BASE_DIR / "Tangram" / "output" / "predicted_proportions.csv",
    "Seurat": BASE_DIR / "Seurat" / "output" / "predicted_proportions.csv",
}

CELL_TYPES = [
    "B cells",
    "CD4+ T cells",
    "CD8+ T cells",
    "Macrophages",
    "Endothelial",
    "Epithelial",
    "Fibroblasts",
]


def normalize_columns(df):
    """Normalize column names across methods.

    Handles: R's make.names (dots for spaces/+), C2L prefix (q05cell_abundance_w_sf_).
    """
    rename = {}
    for col in df.columns:
        clean = col
        clean = clean.replace("q05cell_abundance_w_sf_", "")
        clean = clean.replace(".", " ")
        clean = clean.replace("CD4  T cells", "CD4+ T cells")
        clean = clean.replace("CD8  T cells", "CD8+ T cells")
        rename[col] = clean.strip()
    return df.rename(columns=rename)


def evaluate_method(gt_df, pred_path, method_name):
    """Evaluate one method against GT."""
    if not pred_path.exists():
        logger.warning(f"{method_name}: {pred_path} not found, skipping")
        return None

    pred_df = pd.read_csv(pred_path, index_col=0)
    pred_df = normalize_columns(pred_df)

    # Align to common AOIs
    common = gt_df.index.intersection(pred_df.index)
    if len(common) == 0:
        logger.error(f"{method_name}: no common AOIs between GT and predictions")
        return None
    logger.info(f"{method_name}: {len(common)} common AOIs")

    gt = gt_df.loc[common]
    pred = pred_df.loc[common]

    results = {"method": method_name, "n_aois": len(common)}

    # Per-type metrics
    per_type = {}
    for ct in CELL_TYPES:
        if ct not in pred.columns:
            logger.warning(f"{method_name}: missing column '{ct}'")
            continue
        gt_vals = gt[ct].values
        pred_vals = pred[ct].values

        r, p = stats.pearsonr(gt_vals, pred_vals)
        rmse = np.sqrt(np.mean((gt_vals - pred_vals) ** 2))
        mae = np.mean(np.abs(gt_vals - pred_vals))
        per_type[ct] = {"r": round(r, 4), "p": round(p, 6), "rmse": round(rmse, 4), "mae": round(mae, 4)}

    results["per_type"] = per_type

    # Overall metrics
    common_types = [ct for ct in CELL_TYPES if ct in pred.columns and ct in gt.columns]
    gt_mat = gt[common_types].values
    pred_mat = pred[common_types].values

    # Flattened Pearson r
    r_flat, _ = stats.pearsonr(gt_mat.flatten(), pred_mat.flatten())
    results["overall_r"] = round(r_flat, 4)
    results["overall_rmse"] = round(np.sqrt(np.mean((gt_mat - pred_mat) ** 2)), 4)
    results["overall_mae"] = round(np.mean(np.abs(gt_mat - pred_mat)), 4)

    # Mean per-type r
    type_rs = [v["r"] for v in per_type.values()]
    results["mean_r"] = round(np.mean(type_rs), 4) if type_rs else None

    # JSD
    jsd_vals = []
    for i in range(len(gt)):
        gt_row = gt_mat[i] + 1e-10
        pred_row = pred_mat[i] + 1e-10
        gt_row /= gt_row.sum()
        pred_row /= pred_row.sum()
        jsd = jensenshannon(gt_row, pred_row)
        if not np.isnan(jsd):
            jsd_vals.append(jsd)
    results["mean_jsd"] = round(np.mean(jsd_vals), 4) if jsd_vals else None
    results["median_jsd"] = round(np.median(jsd_vals), 4) if jsd_vals else None

    return results


def main():
    output_dir = BASE_DIR / "evaluation" / "results"
    output_dir.mkdir(parents=True, exist_ok=True)

    gt_df = pd.read_csv(GT_PATH, index_col=0)
    logger.info(f"GT: {gt_df.shape} ({list(gt_df.columns)})")

    all_results = {}
    for method, pred_path in METHODS.items():
        result = evaluate_method(gt_df, pred_path, method)
        if result:
            all_results[method] = result
            logger.info(f"{method}: overall_r={result['overall_r']}, mean_r={result['mean_r']}")

    # Save JSON
    json_path = output_dir / "pseudo_geomx_all_methods_props.json"
    with open(json_path, "w") as f:
        json.dump(all_results, f, indent=2)
    logger.info(f"Wrote {json_path}")

    # Save summary CSV
    rows = []
    for method, r in all_results.items():
        row = {
            "method": method,
            "overall_r": r["overall_r"],
            "mean_r": r["mean_r"],
            "overall_rmse": r["overall_rmse"],
            "mean_jsd": r["mean_jsd"],
        }
        for ct, m in r["per_type"].items():
            row[f"{ct}_r"] = m["r"]
        rows.append(row)

    summary = pd.DataFrame(rows).sort_values("overall_r", ascending=False)
    csv_path = output_dir / "pseudo_geomx_all_methods_props_summary.csv"
    summary.to_csv(csv_path, index=False)
    logger.info(f"Wrote {csv_path}")

    # Print summary table
    print("\n" + "=" * 60)
    print("Pseudo-GeoMx Proportion Benchmark Results")
    print("=" * 60)
    print(summary[["method", "overall_r", "mean_r", "overall_rmse", "mean_jsd"]].to_string(index=False))


if __name__ == "__main__":
    main()
