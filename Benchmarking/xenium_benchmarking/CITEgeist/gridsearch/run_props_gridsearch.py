#!/usr/bin/env python -u
"""Proportions grid search runner — one (combo, region) per SLURM task.

Loads a parameter combo from combos_props.csv, runs QP proportion
deconvolution on one Xenium region, evaluates against GT, and writes
a result JSON atomically.
"""

import argparse
import json
import logging
import os
import sys
import tempfile
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats
from scipy.spatial.distance import jensenshannon

REPO_ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO_ROOT))

BENCHMARK_ROOT = REPO_ROOT / "Benchmarking" / "xenium_benchmarking"
sys.path.insert(0, str(BENCHMARK_ROOT))
from benchmark_constants import ACHIEVABLE_7_CELL_PROFILE_DICT

from CITEgeist.model.citegeist_model import CitegeistModel

DATA_DIR = REPO_ROOT / "Benchmarking" / "xenium_pseudovisium" / "data_protein_gt"
GT_DIR = (
    REPO_ROOT / "Benchmarking" / "xenium_benchmarking" / "ground_truth_singler" / "ground_truth_7type" / "ground_truth"
)
CELL_TYPES = list(ACHIEVABLE_7_CELL_PROFILE_DICT.keys())

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
)
logger = logging.getLogger(__name__)


def load_region(region_id: int):
    """Load GEX and CITE h5ad for a region."""
    h5ad_dir = DATA_DIR / "h5ad_objects"
    gex = sc.read_h5ad(h5ad_dir / f"Xenium_region_{region_id}_GEX.h5ad")
    cite = sc.read_h5ad(h5ad_dir / f"Xenium_region_{region_id}_CITE.h5ad")
    return gex, cite


def load_nuclei_counts(nuclei_dir: Path, region_id: int, protein_adata):
    """Load StarDist nuclei counts and inject into protein_adata.obs."""
    path = nuclei_dir / f"stardist_counts_region_{region_id}.csv"
    counts_df = pd.read_csv(path, index_col=0)
    counts_series = counts_df["nuclei_count_raw"]
    max_count = counts_series.max()
    normalized = counts_series / max_count if max_count > 0 else counts_series
    protein_adata.obs["nuclei_count_target"] = normalized.reindex(protein_adata.obs_names, fill_value=0.0).values


def load_gt(region_id: int) -> pd.DataFrame:
    """Load ground truth proportions for a region (singler 7-type)."""
    gt_path = GT_DIR / f"Xenium_region_{region_id}_prop.csv"
    return pd.read_csv(gt_path, index_col=0)


def calculate_metrics(gt_df: pd.DataFrame, pred_df: pd.DataFrame) -> dict:
    """Compute evaluation metrics: Pearson r, RMSE, MAE, JSD."""
    common_spots = sorted(set(gt_df.index) & set(pred_df.index))
    common_types = [c for c in CELL_TYPES if c in gt_df.columns and c in pred_df.columns]

    gt = gt_df.loc[common_spots, common_types].values
    pred = pred_df.loc[common_spots, common_types].values

    if gt.max() > 1.5:
        row_sums = gt.sum(axis=1, keepdims=True)
        row_sums[row_sums == 0] = 1.0
        gt = gt / row_sums
    if pred.max() > 1.5:
        row_sums = pred.sum(axis=1, keepdims=True)
        row_sums[row_sums == 0] = 1.0
        pred = pred / row_sums

    metrics = {}

    per_type_r = {}
    for i, ct in enumerate(common_types):
        r, _ = stats.pearsonr(gt[:, i], pred[:, i])
        per_type_r[ct] = float(r)
    metrics["per_type_r"] = per_type_r

    r_overall, _ = stats.pearsonr(gt.flatten(), pred.flatten())
    metrics["overall_pearson_r"] = float(r_overall)
    metrics["overall_rmse"] = float(np.sqrt(np.mean((gt - pred) ** 2)))
    metrics["overall_mae"] = float(np.mean(np.abs(gt - pred)))

    jsd_vals = []
    for i in range(len(common_spots)):
        gt_row = gt[i] + 1e-10
        pred_row = pred[i] + 1e-10
        gt_row = gt_row / gt_row.sum()
        pred_row = pred_row / pred_row.sum()
        jsd = jensenshannon(gt_row, pred_row)
        if not np.isnan(jsd):
            jsd_vals.append(float(jsd))
    metrics["mean_jsd"] = float(np.mean(jsd_vals)) if jsd_vals else float("nan")

    return metrics


def atomic_write_json(data: dict, path: Path):
    """Write JSON atomically via tmp + rename."""
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_path = tempfile.mkstemp(dir=path.parent, suffix=".tmp")
    try:
        with os.fdopen(fd, "w") as f:
            json.dump(data, f, indent=2)
        os.rename(tmp_path, str(path))
    except Exception:
        os.unlink(tmp_path)
        raise


def run_one(combo_id: int, region_id: int, params: dict, nuclei_dir: Path, output_dir: Path):
    """Run QP proportions for one (combo, region) pair."""
    logger.info("combo=%d region=%d params=%s", combo_id, region_id, params)

    gex, cite = load_region(region_id)
    load_nuclei_counts(nuclei_dir, region_id, cite)
    gt_df = load_gt(region_id)

    sample_name = f"combo_{combo_id}_region_{region_id}"
    model = CitegeistModel(
        sample_name=sample_name,
        output_folder=str(output_dir / str(combo_id)),
        simulation=True,
        gene_expression_adata=gex,
        antibody_capture_adata=cite,
    )
    model.load_cell_profile_dict(ACHIEVABLE_7_CELL_PROFILE_DICT)
    model.filter_gex(min_counts=25)
    model.preprocess_gex(target_sum=10000)
    model.preprocess_antibody()

    use_nuclei = params["nuclei_prior_lambda"] > 0

    t0 = time.time()
    model.run_cell_proportion_model(
        method="qp",
        lambda_reg=params["lambda_reg"],
        alpha=params["alpha_elastic"],
        lambda_laplacian=params["lambda_laplacian"],
        detection_gate_ub=params["detection_gate_ub"],
        use_nuclei_prior=use_nuclei,
        nuclei_prior_lambda=params["nuclei_prior_lambda"],
        validation_warn_only=True,
        lambda_coverage=1.0,
    )
    runtime = time.time() - t0

    props_df = model.results["cell_prop"][CELL_TYPES]

    result_dir = output_dir / str(combo_id)
    result_dir.mkdir(parents=True, exist_ok=True)
    props_df.to_csv(result_dir / f"region_{region_id}_predictions.csv")

    metrics = calculate_metrics(gt_df, props_df)

    result = {
        "combo_id": combo_id,
        "region_id": region_id,
        "params": params,
        "metrics": metrics,
        "runtime_seconds": round(runtime, 1),
    }

    atomic_write_json(result, result_dir / f"region_{region_id}.json")
    logger.info(
        "combo=%d region=%d r=%.4f runtime=%.1fs",
        combo_id,
        region_id,
        metrics["overall_pearson_r"],
        runtime,
    )


def main():
    parser = argparse.ArgumentParser(description="Props grid search runner")
    parser.add_argument("--combo-id", type=int, required=True)
    parser.add_argument("--region-id", type=int, required=True)
    parser.add_argument(
        "--grid-file",
        type=Path,
        default=Path(__file__).parent / "combos_props.csv",
    )
    parser.add_argument(
        "--nuclei-dir",
        type=Path,
        default=REPO_ROOT / "Benchmarking" / "xenium_benchmarking" / "CITEgeist" / "output" / "stardist_counts",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(__file__).parent / "results" / "props",
    )
    args = parser.parse_args()

    grid_df = pd.read_csv(args.grid_file)
    n_combos = len(grid_df)
    if not (0 <= args.combo_id < n_combos):
        logger.error("combo_id %d out of range [0, %d)", args.combo_id, n_combos)
        sys.exit(1)
    if not (0 <= args.region_id <= 4):
        logger.error("region_id %d out of range [0, 4]", args.region_id)
        sys.exit(1)

    row = grid_df.iloc[args.combo_id]
    params = {
        "alpha_elastic": float(row["alpha_elastic"]),
        "lambda_reg": float(row["lambda_reg"]),
        "lambda_laplacian": float(row["lambda_laplacian"]),
        "detection_gate_ub": float(row["detection_gate_ub"]),
        "nuclei_prior_lambda": float(row["nuclei_prior_lambda"]),
    }

    run_one(args.combo_id, args.region_id, params, args.nuclei_dir, args.output_dir)


if __name__ == "__main__":
    main()
