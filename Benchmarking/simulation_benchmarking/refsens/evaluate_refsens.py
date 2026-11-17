"""
Evaluate reference-sensitivity benchmark results.

Loads predictions from all methods x conditions x ref_conditions,
computes Pearson r, RMSE, MAE, JSD per replicate, and compares
against the original adversarial_aggregated_metrics.csv.

Usage:
    python -u evaluate_refsens.py [--results-dir ./results] [--output-dir ./evaluation]
"""

import sys as _sys
from pathlib import Path as _Path

_sys.path.insert(0, str(_Path(__file__).resolve().parents[2]))
import _paths as _cg  # noqa: E402

import argparse
import logging
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial.distance import jensenshannon
from scipy.stats import pearsonr

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)

REFSENS_DIR = Path(__file__).resolve().parent
REPLICATES_BASE = (
    str(_cg.REPO_ROOT / "Benchmarking" / "simulation_benchmarking" / "replicates")
)
ORIGINAL_METRICS = (
    str(_cg.REPO_ROOT / "Benchmarking" / "simulation_benchmarking" / "Figures" / "adversarial_aggregated_metrics.csv")
)

# Method-specific prediction file patterns
METHOD_CONFIGS = {
    "RCTD": {
        "pred_pattern": "Wu_ST_{rep}_RCTD_deconv_predictions.csv",
        "subdir": "",  # predictions directly in condition/ref_cond/
    },
    "Seurat": {
        "pred_pattern": "Wu_ST_{rep}_Seurat_deconv_predictions.csv",
        "subdir": "",
    },
    "Tangram": {
        "pred_pattern": "cell_type_proportions.csv",
        "subdir": "tangram_map_{rep}",
    },
    "Cell2Location": {
        "pred_pattern": "cell2loc_deconv_predictions.csv",
        "subdir": "c2l_map_{rep}",
    },
}

CONDITIONS = ["high_seg", "mixed"]
REF_CONDITIONS = ["full", "downsampled"]
N_REPS = 5


def load_ground_truth(condition: str, rep: int) -> pd.DataFrame:
    """Load ground truth proportions for a replicate."""
    gt_path = os.path.join(
        REPLICATES_BASE, condition, "ST_sim", f"Wu_ST_{rep}_prop.csv"
    )
    gt = pd.read_csv(gt_path, index_col=0).sort_index()
    gt = gt.drop(columns=["spot_x", "spot_y"], errors="ignore")
    gt = gt.sort_index(axis=1)
    return gt


def load_prediction(method: str, condition: str, ref_cond: str, rep: int, results_dir: Path):
    """Load prediction file for a method/condition/ref/replicate."""
    config = METHOD_CONFIGS[method]
    pred_file = config["pred_pattern"].format(rep=rep)

    if config["subdir"]:
        subdir = config["subdir"].format(rep=rep)
        pred_dir = results_dir / method / condition / ref_cond / subdir
    else:
        pred_dir = results_dir / method / condition / ref_cond

    pred_path = pred_dir / pred_file

    # NFS negative-cache workaround: listdir sees files that stat() misses
    try:
        entries = os.listdir(pred_dir)
    except FileNotFoundError:
        return None
    if pred_file not in entries:
        return None

    try:
        pred = pd.read_csv(pred_path, index_col=0).sort_index()
    except (FileNotFoundError, OSError):
        # stat() negative cache — try reading via file descriptor
        try:
            fd = os.open(str(pred_dir), os.O_RDONLY | os.O_DIRECTORY)
            os.close(fd)
            pred = pd.read_csv(pred_path, index_col=0).sort_index()
        except (FileNotFoundError, OSError):
            logger.warning(f"NFS cache miss: {pred_path} visible in listdir but not readable")
            return None
    pred = pred.drop(columns=["spot_x", "spot_y", "spot"], errors="ignore")
    pred = pred.sort_index(axis=1)
    return pred


def compute_metrics(gt: pd.DataFrame, pred: pd.DataFrame) -> dict:
    """Compute benchmark metrics between ground truth and predictions."""
    # Align columns
    common_cols = sorted(set(gt.columns) & set(pred.columns))
    gt = gt[common_cols]
    pred = pred[common_cols]

    # Align rows
    common_idx = sorted(set(gt.index) & set(pred.index))
    gt = gt.loc[common_idx]
    pred = pred.loc[common_idx]

    gt_arr = gt.values
    pred_arr = pred.values

    # Pearson r (flattened)
    corr, _ = pearsonr(gt_arr.flatten(), pred_arr.flatten())

    # RMSE (overall)
    mse = np.mean((gt_arr - pred_arr) ** 2)
    rmse = np.sqrt(mse)

    # MAE (overall)
    mae = np.mean(np.abs(gt_arr - pred_arr))

    # JSD (per spot, then median)
    jsd_values = []
    for i in range(gt_arr.shape[0]):
        if pred_arr[i].sum() > 0:
            jsd = jensenshannon(gt_arr[i], pred_arr[i], base=2)
        else:
            jsd = 1.0
        jsd_values.append(jsd)
    median_jsd = np.median(jsd_values)

    return {
        "Corr": corr,
        "RMSE": rmse,
        "MAE": mae,
        "JSD": median_jsd,
    }


def main():
    parser = argparse.ArgumentParser(description="Evaluate refsens benchmark results")
    parser.add_argument(
        "--results-dir",
        type=str,
        default=str(REFSENS_DIR / "results"),
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default=str(REFSENS_DIR / "evaluation"),
    )
    args = parser.parse_args()

    results_dir = Path(args.results_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Collect results
    rows = []
    missing = []

    for method in METHOD_CONFIGS:
        for condition in CONDITIONS:
            pattern_label = "High Seg" if condition == "high_seg" else "Mixed"
            for ref_cond in REF_CONDITIONS:
                ref_label = "Full" if ref_cond == "full" else "Downsampled"
                for rep in range(N_REPS):
                    pred = load_prediction(method, condition, ref_cond, rep, results_dir)
                    if pred is None:
                        missing.append(f"{method}/{condition}/{ref_cond}/rep{rep}")
                        continue

                    gt = load_ground_truth(condition, rep)
                    metrics = compute_metrics(gt, pred)

                    rows.append({
                        "Method": method,
                        "Replicate": rep,
                        "Reference": ref_label,
                        "Pattern": pattern_label,
                        **metrics,
                    })

    if missing:
        logger.warning(f"Missing {len(missing)} prediction files:")
        for m in missing[:10]:
            logger.warning(f"  {m}")
        if len(missing) > 10:
            logger.warning(f"  ... and {len(missing) - 10} more")

    if not rows:
        logger.error("No predictions found. Nothing to evaluate.")
        sys.exit(1)

    results_df = pd.DataFrame(rows)

    # Save raw results
    raw_path = output_dir / "refsens_metrics.csv"
    results_df.to_csv(raw_path, index=False)
    logger.info(f"Saved raw metrics: {raw_path}")

    # Print summary
    logger.info("\n" + "=" * 70)
    logger.info("Summary: Mean Pearson r by Method x Reference")
    logger.info("=" * 70)

    summary = (
        results_df.groupby(["Method", "Reference"])["Corr"]
        .agg(["mean", "std", "count"])
        .round(4)
    )
    logger.info(f"\n{summary.to_string()}")

    # Compare against original adversarial metrics
    if os.path.exists(ORIGINAL_METRICS):
        logger.info("\n" + "=" * 70)
        logger.info("Comparison with original adversarial_aggregated_metrics.csv")
        logger.info("=" * 70)

        orig_df = pd.read_csv(ORIGINAL_METRICS, index_col=0)

        # Build comparison table
        comparison_rows = []
        for method in results_df["Method"].unique():
            for ref in results_df["Reference"].unique():
                new_sub = results_df[
                    (results_df["Method"] == method) & (results_df["Reference"] == ref)
                ]
                orig_sub = orig_df[
                    (orig_df["Method"] == method) & (orig_df["Reference"] == ref)
                ]

                if len(new_sub) > 0 and len(orig_sub) > 0:
                    comparison_rows.append({
                        "Method": method,
                        "Reference": ref,
                        "Original_Corr": orig_sub["Corr"].mean(),
                        "New_Corr": new_sub["Corr"].mean(),
                        "Delta_Corr": new_sub["Corr"].mean() - orig_sub["Corr"].mean(),
                        "Original_n": len(orig_sub),
                        "New_n": len(new_sub),
                    })
                elif len(new_sub) > 0:
                    comparison_rows.append({
                        "Method": method,
                        "Reference": ref,
                        "Original_Corr": np.nan,
                        "New_Corr": new_sub["Corr"].mean(),
                        "Delta_Corr": np.nan,
                        "Original_n": 0,
                        "New_n": len(new_sub),
                    })

        if comparison_rows:
            comp_df = pd.DataFrame(comparison_rows).round(4)
            comp_path = output_dir / "comparison_table.csv"
            comp_df.to_csv(comp_path, index=False)
            logger.info(f"\n{comp_df.to_string(index=False)}")
            logger.info(f"\nSaved comparison: {comp_path}")
    else:
        logger.warning(f"Original metrics not found at {ORIGINAL_METRICS}")

    logger.info("\nEvaluation complete.")


if __name__ == "__main__":
    main()
