"""
Test beta-weighted gating classifier vs QP argmax vs supervised logistic.

The QP learns per-marker betas (signal quality), then a simple gating
classifier uses those betas as weights with positive + negative evidence.
"""
import json
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).parent.parent.parent.parent.parent
sys.path.insert(0, str(REPO_ROOT / "CITEgeist"))
sys.path.insert(0, str(REPO_ROOT / "Benchmarking" / "xenium_pseudovisium" / "src"))
sys.path.insert(0, str(REPO_ROOT / "Benchmarking" / "xenium_benchmarking" / "CITEgeist" / "src"))
sys.path.insert(0, str(REPO_ROOT / "Benchmarking" / "xenium_benchmarking"))

from model.gurobi_impl import (
    optimize_cell_proportions_per_marker,
    map_antibodies_to_profiles_v2,
    beta_weighted_classification,
    supervised_cell_classification,
)
from load_xenium_singlecell import load_xenium_singlecell
from benchmark_constants import ACHIEVABLE_7_CELL_PROFILE_DICT

logger = logging.getLogger(__name__)
OUTPUT_DIR = REPO_ROOT / "Benchmarking" / "xenium_benchmarking" / "CITEgeist" / "output_cell_resolution"
GT_PATH = REPO_ROOT / "Benchmarking" / "xenium_pseudovisium" / "data_protein_gt" / "cell_type_assignments.csv"


def run_comparison(region_id: int = 2, max_cells: int = 5000):
    output_dir = OUTPUT_DIR / "diagnostics"
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info(f"Loading Xenium region {region_id} (max {max_cells} cells)")
    adata_gex, adata_protein = load_xenium_singlecell(
        region_id=region_id, max_cells=max_cells, seed=42
    )

    gt_df = pd.read_csv(GT_PATH, index_col=0)
    common_cells = sorted(set(adata_protein.obs_names) & set(gt_df.index))
    gt_labels = gt_df.loc[common_cells, "cell_type"].values
    coords = adata_protein.obsm["spatial"]

    from sklearn.metrics import accuracy_score, f1_score, classification_report

    # === Step 1: Run QP to learn betas ===
    logger.info("=" * 80)
    logger.info("STEP 1: QP optimization (learning betas)")
    logger.info("=" * 80)

    marker_data, mapped_markers, assignment_matrix, type_names = map_antibodies_to_profiles_v2(
        adata_protein, ACHIEVABLE_7_CELL_PROFILE_DICT
    )

    Y_qp, beta_values, beta_dict, alpha_values = optimize_cell_proportions_per_marker(
        marker_level_data=marker_data,
        marker_names=mapped_markers,
        assignment_matrix=assignment_matrix,
        cell_type_names=type_names,
        lambda_sparse=0.0,
        lambda_laplacian=0.01,
        coords=coords,
        laplacian_k=50,
        max_iterations=5,
    )

    logger.info("Learned betas:")
    for m_name, beta in beta_dict.items():
        logger.info(f"  {m_name:>12}: beta={beta:.4f}")

    # Helper to evaluate predictions
    def evaluate(Y, label):
        dominant = np.argmax(Y, axis=1)
        pred_all = [type_names[d] for d in dominant]
        cell_to_pred = dict(zip(adata_protein.obs_names, pred_all))
        pred_gt = np.array([cell_to_pred[c] for c in common_cells])

        gt_types = set(type_names)
        eval_mask = np.array([gt in gt_types for gt in gt_labels])
        gt_eval = gt_labels[eval_mask]
        pred_eval = pred_gt[eval_mask]

        acc = float(accuracy_score(gt_eval, pred_eval))
        f1_mac = float(f1_score(gt_eval, pred_eval, average="macro", zero_division=0))
        f1_wt = float(f1_score(gt_eval, pred_eval, average="weighted", zero_division=0))

        max_Y = Y.max(axis=1)
        one_hot = float(np.mean(max_Y > 0.9))

        logger.info(f"\n--- {label} ---")
        logger.info(f"  Accuracy: {acc:.3f}  F1_macro: {f1_mac:.3f}  F1_weighted: {f1_wt:.3f}")
        logger.info(f"  Max Y: mean={max_Y.mean():.3f} median={np.median(max_Y):.3f}")
        logger.info(f"  One-hot (>0.9): {one_hot*100:.1f}%")
        logger.info("\n" + classification_report(gt_eval, pred_eval, zero_division=0))

        return {
            "label": label,
            "accuracy": acc,
            "f1_macro": f1_mac,
            "f1_weighted": f1_wt,
            "max_Y_mean": float(max_Y.mean()),
            "one_hot_fraction": one_hot,
        }

    # === Evaluate all approaches ===
    results = []

    # 1. QP argmax baseline
    r_qp = evaluate(Y_qp, "QP argmax")
    results.append(r_qp)

    # 2. Beta-weighted gating
    logger.info("=" * 80)
    logger.info("Beta-weighted gating classifier")
    logger.info("=" * 80)
    Y_gate, scores_gate = beta_weighted_classification(
        marker_level_data=marker_data,
        marker_names=mapped_markers,
        assignment_matrix=assignment_matrix,
        cell_type_names=type_names,
        beta_values=beta_values,
    )
    r_gate = evaluate(Y_gate, "Beta-weighted gating")
    results.append(r_gate)

    # 3. Supervised logistic (p50, best from previous test)
    logger.info("=" * 80)
    logger.info("Supervised logistic (p50) on beta-gated labels")
    logger.info("=" * 80)
    Y_sup_qp, _ = supervised_cell_classification(
        marker_level_data=marker_data,
        Y_qp=Y_qp,
        cell_type_names=type_names,
        confidence_percentile=50.0,
        classifier_type="logistic",
    )
    r_sup_qp = evaluate(Y_sup_qp, "Supervised (QP labels)")
    results.append(r_sup_qp)

    # 4. Supervised logistic trained on beta-gated labels instead of QP labels
    logger.info("=" * 80)
    logger.info("Supervised logistic (p50) on beta-gated labels")
    logger.info("=" * 80)
    Y_sup_gate, _ = supervised_cell_classification(
        marker_level_data=marker_data,
        Y_qp=Y_gate,  # use gating labels instead of QP labels
        cell_type_names=type_names,
        confidence_percentile=50.0,
        classifier_type="logistic",
    )
    r_sup_gate = evaluate(Y_sup_gate, "Supervised (gating labels)")
    results.append(r_sup_gate)

    # === Summary ===
    logger.info("\n" + "=" * 80)
    logger.info("COMPARISON SUMMARY")
    logger.info("=" * 80)
    logger.info(f"{'Method':<30} {'Acc':>6} {'F1mac':>6} {'F1wt':>6} {'maxY':>6} {'1hot%':>6}")
    logger.info("-" * 80)
    for r in results:
        logger.info(f"{r['label']:<30} {r['accuracy']:6.3f} {r['f1_macro']:6.3f} "
                     f"{r['f1_weighted']:6.3f} {r['max_Y_mean']:6.3f} "
                     f"{r['one_hot_fraction']*100:5.1f}%")

    with open(output_dir / "gating_comparison.json", "w") as f:
        json.dump(results, f, indent=2, default=str)

    logger.info(f"\nSaved to {output_dir / 'gating_comparison.json'}")


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    run_comparison()
