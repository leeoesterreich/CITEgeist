#!/usr/bin/env python
"""
Regression benchmark: run cuOPT NB pipeline (with GMM gating, S/log(N+1) cellularity)
on simulation data and compare to stored baselines.

Usage:
    python benchmark_nb_regression.py --condition high_seg --reps 0 1 2 3 4
    python benchmark_nb_regression.py --condition high_seg --reps 0  # quick single-rep test
"""

import sys as _sys
from pathlib import Path as _Path

_sys.path.insert(0, str(_Path(__file__).resolve().parents[3]))
import _paths as _cg  # noqa: E402

import argparse
import json
import logging
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from scipy.spatial.distance import jensenshannon

REPO_ROOT = Path(__file__).resolve().parent.parent.parent.parent.parent
sys.path.insert(0, str(REPO_ROOT))

DEFAULT_SCCUBE_DIR = Path(
    str(_cg.root("DATA") / "Wu_Visium" / "Simulations" / "scCube_12k" / "replicates")
)

SIMULATION_CELL_PROFILE_DICT = {
    "B-cells": {"Major": ["B-cells_Protein_1", "B-cells_Protein_2"]},
    "CAFs": {"Major": ["CAFs_Protein_1", "CAFs_Protein_2"]},
    "Cancer Epithelial": {"Major": ["Cancer Epithelial_Protein_1", "Cancer Epithelial_Protein_2"]},
    "Endothelial": {"Major": ["Endothelial_Protein_1", "Endothelial_Protein_2"]},
    "Myeloid": {"Major": ["Myeloid_Protein_1", "Myeloid_Protein_2"]},
    "Normal Epithelial": {"Major": ["Normal Epithelial_Protein_1", "Normal Epithelial_Protein_2"]},
    "PVL": {"Major": ["PVL_Protein_1", "PVL_Protein_2"]},
    "Plasmablasts": {"Major": ["Plasmablasts_Protein_1", "Plasmablasts_Protein_2"]},
    "T-cells": {"Major": ["T-cells_Protein_1", "T-cells_Protein_2"]},
}

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

# Known baselines from previous NB runs (results_nb/high_seg/gt_ablation/)
BASELINES = {
    "high_seg": {
        0: {"pearson_r": 0.8657, "rmse": 0.1354, "sc_acc": 0.7324},
        1: {"pearson_r": 0.9333, "rmse": 0.1003, "sc_acc": None},
        2: {"pearson_r": 0.9408, "rmse": 0.0950, "sc_acc": None},
        3: {"pearson_r": 0.9266, "rmse": 0.1033, "sc_acc": None},
        4: {"pearson_r": 0.9249, "rmse": 0.1049, "sc_acc": None},
    },
    # QP baselines from results/summary_table.csv for comparison
    "high_seg_qp": {"pearson_r": 0.9536, "rmse": 0.1112, "sc_acc_gt_ablation": 0.801},
    "mixed_qp": {"pearson_r": 0.9186, "rmse": 0.0851, "sc_acc_gt_ablation": 0.8404},
}

# Regression tolerance: how much worse r can be before flagging
REGRESSION_TOLERANCE = 0.02  # 2% drop allowed


def load_ground_truth_proportions(replicate_id, condition, sccube_dir):
    prop_path = sccube_dir / condition / "ST_sim" / f"Wu_ST_{replicate_id}_prop.csv"
    df = pd.read_csv(prop_path, index_col=0)
    df = df.drop(columns=["spot_x", "spot_y"], errors="ignore")
    return df


def load_cell_index(sccube_dir, condition, replicate_id):
    index_path = sccube_dir / condition / "ST_sim" / f"Wu_ST_{replicate_id}_index.csv"
    return pd.read_csv(index_path, index_col=0)


def discretize_proportions(proportions, total):
    prop_sum = proportions.sum()
    if prop_sum > 0:
        proportions = proportions / prop_sum
    scaled = proportions * total
    floors = np.floor(scaled).astype(int)
    remainders = scaled - floors
    deficit = int(total - floors.sum())
    if deficit > 0:
        indices = np.argsort(-remainders)
        for i in range(min(deficit, len(indices))):
            floors[indices[i]] += 1
    elif deficit < 0:
        indices = np.argsort(remainders)
        for i in range(min(-deficit, len(indices))):
            if floors[indices[i]] > 0:
                floors[indices[i]] -= 1
    return floors


def hungarian_sc_accuracy(gt_index, pred_counts_df):
    total_correct = 0
    total_cells = 0
    for spot in pred_counts_df.index:
        spot_gt = gt_index[gt_index["spot"] == spot]["Cell_type"].tolist()
        if len(spot_gt) == 0:
            continue
        gt_counts = {}
        for label in spot_gt:
            gt_counts[label] = gt_counts.get(label, 0) + 1
        pred_counts = pred_counts_df.loc[spot].to_dict()
        types = set(gt_counts.keys()) | set(pred_counts.keys())
        correct = sum(min(gt_counts.get(t, 0), pred_counts.get(t, 0)) for t in types)
        total_correct += correct
        total_cells += len(spot_gt)
    return total_correct / total_cells if total_cells > 0 else 0.0


def evaluate_proportions(pred_props, gt_props):
    common_spots = pred_props.index.intersection(gt_props.index)
    common_types = pred_props.columns.intersection(gt_props.columns)
    pred = pred_props.loc[common_spots, common_types].values
    gt = gt_props.loc[common_spots, common_types].values

    r, p = pearsonr(pred.flatten(), gt.flatten())
    rmse = float(np.sqrt(np.mean((pred - gt) ** 2)))

    jsd_list = []
    for i in range(len(common_spots)):
        p_row = np.clip(pred[i], 0, None)
        g_row = np.clip(gt[i], 0, None)
        ps, gs = p_row.sum(), g_row.sum()
        if ps > 0 and gs > 0:
            jsd_list.append(float(jensenshannon(p_row / ps, g_row / gs) ** 2))

    per_type_r = {}
    for ct in common_types:
        ct_pred = pred_props.loc[common_spots, ct].values
        ct_gt = gt_props.loc[common_spots, ct].values
        if ct_gt.std() > 0 and ct_pred.std() > 0:
            ct_r, _ = pearsonr(ct_pred, ct_gt)
            per_type_r[ct] = float(ct_r)

    return {
        "pearson_r": float(r),
        "rmse": rmse,
        "mean_jsd": float(np.mean(jsd_list)) if jsd_list else float("nan"),
        "n_spots": len(common_spots),
        "n_cell_types": len(common_types),
        "per_type_r": per_type_r,
    }


def run_nb_benchmark_single(replicate_id, condition, sccube_dir, output_dir):
    """Run NB pipeline on one replicate with gt_ablation (perfect nuclei counts)."""
    import scanpy as sc
    from CITEgeist.model.citegeist_model import CitegeistModel

    logger.info("=" * 60)
    logger.info("NB REGRESSION BENCHMARK: rep=%d cond=%s", replicate_id, condition)
    logger.info("=" * 60)

    # Load data
    h5ad_dir = sccube_dir / condition / "h5ad_objects"
    adata_cite = sc.read_h5ad(h5ad_dir / f"Wu_rep_{replicate_id}_CITE.h5ad")
    adata_gex = sc.read_h5ad(h5ad_dir / f"Wu_rep_{replicate_id}_GEX.h5ad")
    gt_props = load_ground_truth_proportions(replicate_id, condition, sccube_dir)
    cell_index = load_cell_index(sccube_dir, condition, replicate_id)

    # GT nuclei counts (perfect segmentation ablation)
    nuclei_counts = cell_index.groupby("spot").size()
    nuclei_counts.name = "nuclei_count"

    # Setup model
    model_output = str(output_dir / condition / f"Wu_rep_{replicate_id}" / "citegeist")
    model = CitegeistModel(
        sample_name=f"Wu_rep_{replicate_id}",
        output_folder=model_output,
        simulation=True,
        gene_expression_adata=adata_gex,
        antibody_capture_adata=adata_cite,
    )

    # Preprocessing — NB needs preprocess_antibody() for raw_counts layer
    model.preprocess_antibody()
    model.preprocess_gex(target_sum=10000)

    # Simulation markers need explicit profile dict (unlike Xenium real data which
    # uses the hardcoded 18-marker panel via build_marker_config)
    model.cell_profile_dict = SIMULATION_CELL_PROFILE_DICT

    # Run NB proportion model with production defaults
    logger.info("Running NB proportions (GMM detection gating, cellularity) ...")
    start = time.time()
    # Run NB with production defaults (GMM detection, cellularity from GT nuclei)
    _global, finetuned = model.run_cell_proportion_model(
        method="nb",
        nb_device="cpu",
        nb_n_iter=100,
        nb_use_detection=True,
        nb_detection_threshold=0.5,
        nuclei_counts=nuclei_counts,
    )
    elapsed = time.time() - start
    logger.info("NB completed in %.1fs", elapsed)

    # Evaluate proportions
    prop_metrics = evaluate_proportions(finetuned, gt_props)
    logger.info(
        "Proportions: r=%.4f RMSE=%.4f JSD=%.4f",
        prop_metrics["pearson_r"],
        prop_metrics["rmse"],
        prop_metrics["mean_jsd"],
    )

    # Discretize + SC accuracy
    common = finetuned.index.intersection(nuclei_counts.index)
    cell_counts_df = pd.DataFrame(0, index=common, columns=finetuned.columns)
    for spot in common:
        cell_counts_df.loc[spot] = discretize_proportions(finetuned.loc[spot].values, int(nuclei_counts[spot]))
    sc_acc = hungarian_sc_accuracy(cell_index, cell_counts_df)
    logger.info("Hungarian SC accuracy: %.4f", sc_acc)

    result = {
        "method": "nb",
        "replicate_id": replicate_id,
        "condition": condition,
        "mode": "gt_ablation",
        "proportion_metrics": prop_metrics,
        "hungarian_sc_accuracy": sc_acc,
        "total_nuclei_detected": int(nuclei_counts.sum()),
        "timings": {"nb_deconv_sec": elapsed},
    }

    # Save
    out_dir = output_dir / condition / "gt_ablation" / f"Wu_rep_{replicate_id}"
    out_dir.mkdir(parents=True, exist_ok=True)
    with open(out_dir / "benchmark_results.json", "w") as f:
        json.dump(result, f, indent=2, default=str)
    finetuned.to_csv(out_dir / "cell_proportions.csv")

    return result


def check_regression(results, condition):
    """Compare results to baselines and report pass/fail."""
    baselines = BASELINES.get(condition, {})
    all_pass = True
    lines = []

    for result in results:
        rep = result["replicate_id"]
        r = result["proportion_metrics"]["pearson_r"]
        rmse = result["proportion_metrics"]["rmse"]
        sc = result["hungarian_sc_accuracy"]

        baseline = baselines.get(rep)
        if baseline:
            r_delta = r - baseline["pearson_r"]
            rmse_delta = rmse - baseline["rmse"]
            r_pass = r_delta >= -REGRESSION_TOLERANCE
            status = "PASS" if r_pass else "FAIL"
            if not r_pass:
                all_pass = False
            lines.append(
                f"  Rep {rep}: r={r:.4f} (baseline {baseline['pearson_r']:.4f}, "
                f"delta={r_delta:+.4f}) RMSE={rmse:.4f} SC={sc:.4f} [{status}]"
            )
        else:
            lines.append(f"  Rep {rep}: r={r:.4f} RMSE={rmse:.4f} SC={sc:.4f} [NO BASELINE]")

    # Aggregate
    rs = [r["proportion_metrics"]["pearson_r"] for r in results]
    rmses = [r["proportion_metrics"]["rmse"] for r in results]
    scs = [r["hungarian_sc_accuracy"] for r in results]

    print("\n" + "=" * 60)
    print(f"NB REGRESSION BENCHMARK: {condition} ({len(results)} reps)")
    print("=" * 60)
    print(f"  Mean r:    {np.mean(rs):.4f} +/- {np.std(rs):.4f}")
    print(f"  Mean RMSE: {np.mean(rmses):.4f} +/- {np.std(rmses):.4f}")
    print(f"  Mean SC:   {np.mean(scs):.4f} +/- {np.std(scs):.4f}")
    print()
    print("Per-replicate comparison:")
    for line in lines:
        print(line)
    print()

    qp_baseline = BASELINES.get(f"{condition}_qp", {})
    if qp_baseline:
        print(
            f"  Reference QP baseline: r={qp_baseline.get('pearson_r', 'N/A')} "
            f"RMSE={qp_baseline.get('rmse', 'N/A')}"
        )

    overall = "PASS" if all_pass else "REGRESSION DETECTED"
    print(f"\n  Overall: [{overall}]")
    print("=" * 60)
    return all_pass


def main():
    parser = argparse.ArgumentParser(description="NB regression benchmark on simulation data")
    parser.add_argument("--condition", default="high_seg", help="Simulation condition")
    parser.add_argument("--reps", type=int, nargs="+", default=[0, 1, 2, 3, 4])
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=REPO_ROOT / "Benchmarking" / "simulation_benchmarking" / "CITEgeist" / "results_nb_regression",
    )
    parser.add_argument("--sccube-dir", type=Path, default=DEFAULT_SCCUBE_DIR)
    args = parser.parse_args()

    results = []
    for rep in args.reps:
        result = run_nb_benchmark_single(rep, args.condition, args.sccube_dir, args.output_dir)
        results.append(result)

    passed = check_regression(results, args.condition)
    sys.exit(0 if passed else 1)


if __name__ == "__main__":
    main()
