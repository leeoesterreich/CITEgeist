# Recovered from 9c1c372b^ (body unchanged) for v40 baseline reproduction.
# Reference only — do not run.
#!/usr/bin/env python
"""
Gold standard QP benchmark for Xenium pseudo-Visium evaluation.

Pipeline: cuOPT QP + GMM detection gating + S/log(N+1) cellularity correction.
This is the production Module 3 proportion estimation (r=0.766).

Usage:
    python run_qp_only_benchmark.py [--output-dir /path/to/output]
"""

import sys as _sys
from pathlib import Path as _Path

_sys.path.insert(0, str(_Path(__file__).resolve().parents[2]))
import _paths as _cg  # noqa: E402

import argparse
import json
import logging
import sys
import time
import traceback
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import pearsonr

REPO_ROOT = Path(str(_cg.REPO_ROOT))
sys.path.insert(0, str(REPO_ROOT))

BENCHMARK_ROOT = REPO_ROOT / "Benchmarking" / "xenium_benchmarking"
sys.path.insert(0, str(BENCHMARK_ROOT))

from benchmark_constants import ACHIEVABLE_6_CELL_PROFILE_DICT
from CITEgeist.model.citegeist_model import CitegeistModel

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)

DATA_DIR = REPO_ROOT / "Benchmarking" / "xenium_pseudovisium" / "data_protein_gt"
GT_DIR = BENCHMARK_ROOT / "ground_truth_singler" / "ground_truth_6type" / "ground_truth"
DEFAULT_OUTPUT_DIR = BENCHMARK_ROOT / "CITEgeist" / "output_qp_singler_6type"

CELL_TYPES_6 = ["B cells", "T cells", "Macrophages", "Endothelial", "Epithelial", "Fibroblasts"]


def parse_args():
    parser = argparse.ArgumentParser(description="Run QP canonical 6-type benchmark on all Xenium regions")
    parser.add_argument(
        "--output-dir", type=str, default=str(DEFAULT_OUTPUT_DIR), help="Directory to write predictions and stats JSON"
    )
    return parser.parse_args()


def get_nuclei_counts(region_id: int) -> pd.Series:
    """Get per-spot nuclei counts from Xenium cell-to-spot mapping."""
    mapping = pd.read_csv(DATA_DIR / "cell_to_spot_mapping.csv", index_col=0)
    cell_types = pd.read_csv(DATA_DIR / "cell_type_assignments.csv", index_col=0)

    region = mapping[(mapping["region_id"] == region_id) & (mapping["spot_idx"] != -1)]
    known = region.merge(cell_types, left_index=True, right_index=True)
    known = known[known["cell_type"] != "Unknown"]
    return known.groupby("spot_id").size()


def evaluate_against_gt(props: pd.DataFrame, region_id: int) -> dict:
    """Compute per-type and overall Pearson r against SingleR 6-type GT."""
    gt_path = GT_DIR / f"Xenium_region_{region_id}_prop.csv"
    if not gt_path.exists():
        logger.warning("GT not found: %s", gt_path)
        return {}

    gt_df = pd.read_csv(gt_path, index_col=0)
    common = props.index.intersection(gt_df.index)
    if len(common) < 10:
        logger.warning("Too few common spots: %d", len(common))
        return {}

    per_type_r = {}
    all_pred, all_gt = [], []
    for ct in CELL_TYPES_6:
        if ct not in props.columns or ct not in gt_df.columns:
            continue
        p_vals = props.loc[common, ct].values.astype(float)
        g_vals = gt_df.loc[common, ct].values.astype(float)
        all_pred.extend(p_vals)
        all_gt.extend(g_vals)
        if p_vals.std() > 1e-6 and g_vals.std() > 1e-6:
            r, _ = pearsonr(p_vals, g_vals)
            per_type_r[ct] = round(float(r), 4)

    overall_r, _ = pearsonr(all_pred, all_gt)
    pred_arr = np.array(all_pred)
    gt_arr = np.array(all_gt)
    rmse = float(np.sqrt(np.mean((pred_arr - gt_arr) ** 2)))
    mae = float(np.mean(np.abs(pred_arr - gt_arr)))

    return {
        "n_spots": len(common),
        "overall_pearson_r": round(float(overall_r), 4),
        "rmse": round(rmse, 4),
        "mae": round(mae, 4),
        "per_type_r": per_type_r,
    }


def run_region(region_id: int, output_dir: Path) -> dict:
    sample_name = f"Xenium_region_{region_id}"
    logger.info("=" * 60)
    logger.info("Region %d: %s", region_id, sample_name)

    gex_path = DATA_DIR / "h5ad_objects" / f"{sample_name}_GEX.h5ad"
    cite_path = DATA_DIR / "h5ad_objects" / f"{sample_name}_CITE.h5ad"

    gex_adata = sc.read_h5ad(gex_path)
    cite_adata = sc.read_h5ad(cite_path)
    nuclei_counts = get_nuclei_counts(region_id)

    t0 = time.time()
    model = CitegeistModel(
        sample_name=sample_name,
        output_folder="/tmp/qp_only_benchmark",
        simulation=True,
        gene_expression_adata=gex_adata.copy(),
        antibody_capture_adata=cite_adata.copy(),
    )
    model.filter_gex(min_counts=25)
    model.preprocess_gex(target_sum=10000)
    model.preprocess_antibody()
    model.load_cell_profile_dict(ACHIEVABLE_6_CELL_PROFILE_DICT)
    # Gold standard: cuOPT QP + GMM gating (default on) + S/log(N+1) cellularity
    model.run_cell_proportion_model(
        method="qp",
        lambda_laplacian=0.0,
        nuclei_counts=nuclei_counts,
    )

    elapsed = time.time() - t0
    logger.info("QP done in %.1fs", elapsed)

    props = model.results["cell_prop"]
    for ct in CELL_TYPES_6:
        if ct not in props.columns:
            props[ct] = 0.0
    props = props[CELL_TYPES_6]

    # Renormalize
    row_sums = props.sum(axis=1)
    props = props.div(row_sums.replace(0, 1), axis=0)

    region_dir = output_dir / sample_name
    region_dir.mkdir(parents=True, exist_ok=True)
    out_path = region_dir / f"{sample_name}_deconv_predictions.csv"
    props.to_csv(out_path)
    logger.info("Saved: %s (%d spots x %d types)", out_path, *props.shape)

    # Evaluate against SingleR 6-type GT
    eval_metrics = evaluate_against_gt(props, region_id)
    stats = {"region_id": region_id, "time_s": round(elapsed, 1), **eval_metrics}

    stats_path = region_dir / "region_stats.json"
    with open(stats_path, "w") as f:
        json.dump(stats, f, indent=2)
    logger.info("Region %d: overall_r=%.4f", region_id, stats.get("overall_pearson_r", float("nan")))

    return stats


def main():
    args = parse_args()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info("Output dir: %s", output_dir)

    all_stats = []
    failed = []
    for rid in range(5):
        try:
            stats = run_region(rid, output_dir)
            all_stats.append(stats)
        except Exception as e:
            logger.error("Region %d FAILED: %s", rid, e)
            traceback.print_exc()
            failed.append(rid)

    # Write aggregate canonical summary
    if all_stats:
        rs = [s["overall_pearson_r"] for s in all_stats if "overall_pearson_r" in s]
        per_type_by_ct = {ct: [] for ct in CELL_TYPES_6}
        for s in all_stats:
            for ct, r in s.get("per_type_r", {}).items():
                per_type_by_ct[ct].append(r)

        summary = {
            "n_regions": len(all_stats),
            "overall_mean_r": round(float(np.mean(rs)), 4),
            "overall_std_r": round(float(np.std(rs)), 4),
            "per_region_r": {s["region_id"]: s.get("overall_pearson_r") for s in all_stats},
            "per_type_mean_r": {ct: round(float(np.mean(vals)), 4) for ct, vals in per_type_by_ct.items() if vals},
        }
        summary_path = output_dir / "canonical_summary_6type.json"
        with open(summary_path, "w") as f:
            json.dump(summary, f, indent=2)
        logger.info("Summary saved: %s", summary_path)
        logger.info("Overall mean r = %.4f ± %.4f", summary["overall_mean_r"], summary["overall_std_r"])

    if failed:
        logger.error("FAILED regions: %s — exiting with code 1", failed)
        sys.exit(1)

    logger.info("All %d regions complete. Output: %s", len(all_stats), output_dir)


if __name__ == "__main__":
    main()
