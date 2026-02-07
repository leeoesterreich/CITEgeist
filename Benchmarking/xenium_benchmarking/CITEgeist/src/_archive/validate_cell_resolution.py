"""
Validate cell-resolution CITEgeist on Xenium single-cell data.

Runs Module 3 in resolution="cell" mode on Xenium RCC regions using
the curated Achievable-7 profiles, comparing against protein-gated
ground truth cell type labels.
"""
import argparse
import json
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

# Add paths
REPO_ROOT = Path(__file__).parent.parent.parent.parent.parent
sys.path.insert(0, str(REPO_ROOT / "CITEgeist"))
sys.path.insert(0, str(REPO_ROOT / "Benchmarking" / "xenium_pseudovisium" / "src"))
sys.path.insert(0, str(REPO_ROOT / "Benchmarking" / "xenium_benchmarking" / "CITEgeist" / "src"))
sys.path.insert(0, str(REPO_ROOT / "Benchmarking" / "xenium_benchmarking"))

from model.gurobi_impl import (
    optimize_cell_proportions_per_marker,
    map_antibodies_to_profiles_v2,
    estimate_true_expression_cell,
    beta_weighted_classification,
)
from load_xenium_singlecell import load_xenium_singlecell
from benchmark_constants import ACHIEVABLE_7_CELL_PROFILE_DICT

logger = logging.getLogger(__name__)
OUTPUT_DIR = REPO_ROOT / "Benchmarking" / "xenium_benchmarking" / "CITEgeist" / "output_cell_resolution"
GT_PATH = REPO_ROOT / "Benchmarking" / "xenium_pseudovisium" / "data_protein_gt" / "cell_type_assignments.csv"


def run_validation(region_id: int, max_cells: int = 10000):
    """Run cell-resolution validation for one region using Achievable-7 profiles."""
    output_dir = OUTPUT_DIR / f"region_{region_id}"
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load data
    logger.info(f"Loading Xenium region {region_id} (max {max_cells} cells)")
    adata_gex, adata_protein = load_xenium_singlecell(
        region_id=region_id, max_cells=max_cells, seed=42
    )
    logger.info(f"Loaded {adata_gex.shape[0]} cells")

    # Load ground truth (protein-gated cell type assignments)
    gt_df = pd.read_csv(GT_PATH, index_col=0)
    common_cells = sorted(set(adata_gex.obs_names) & set(gt_df.index))
    gt_labels = gt_df.loc[common_cells, "cell_type"].values
    logger.info(f"Ground truth available for {len(common_cells)} / {adata_gex.shape[0]} cells")

    X_protein = adata_protein.X
    if hasattr(X_protein, "toarray"):
        X_protein = X_protein.toarray()
    X_protein = np.asarray(X_protein, dtype=np.float64)
    coords = adata_protein.obsm["spatial"]

    # Save profiles used
    profiles_simple = {k: v["Major"] for k, v in ACHIEVABLE_7_CELL_PROFILE_DICT.items()}
    with open(output_dir / "profiles.json", "w") as f:
        json.dump(profiles_simple, f, indent=2)
    logger.info(f"Using Achievable-7 profiles: {list(ACHIEVABLE_7_CELL_PROFILE_DICT.keys())}")

    # === MODULE 3 PASS 1: Cell Type Soft Classification ===
    logger.info("=" * 60)
    logger.info("MODULE 3 PASS 1: Cell Type Classification (cell resolution, Gurobi QP)")

    marker_data, mapped_marker_names, assignment_matrix, type_names = map_antibodies_to_profiles_v2(
        adata_protein, ACHIEVABLE_7_CELL_PROFILE_DICT
    )

    Y_assignments, beta_values, beta_dict, alpha_values = optimize_cell_proportions_per_marker(
        marker_level_data=marker_data,
        marker_names=mapped_marker_names,
        assignment_matrix=assignment_matrix,
        cell_type_names=type_names,
        lambda_sparse=0.0,
        lambda_laplacian=0.01,
        coords=coords,
        laplacian_k=50,
        max_iterations=5,
    )

    from sklearn.metrics import accuracy_score, f1_score, classification_report

    def evaluate_assignments(Y, label):
        """Evaluate cell type assignments against ground truth."""
        dominant = np.argmax(Y, axis=1)
        pred_all = [type_names[dt] for dt in dominant]
        cell_to_pred = dict(zip(adata_protein.obs_names, pred_all))
        pred_labels = [cell_to_pred[c] for c in common_cells]

        gt_types_in_profiles = set(type_names)
        evaluable_mask = [gt in gt_types_in_profiles for gt in gt_labels]
        gt_eval = gt_labels[evaluable_mask]
        pred_eval = np.array(pred_labels)[evaluable_mask]

        acc = float(accuracy_score(gt_eval, pred_eval))
        f1_mac = float(f1_score(gt_eval, pred_eval, average="macro", zero_division=0))
        f1_wt = float(f1_score(gt_eval, pred_eval, average="weighted", zero_division=0))

        logger.info(f"\n--- {label} ---")
        logger.info(f"Accuracy: {acc:.3f}, F1 macro: {f1_mac:.3f}, F1 weighted: {f1_wt:.3f}")
        logger.info("\n" + classification_report(gt_eval, pred_eval, zero_division=0))

        return acc, f1_mac, f1_wt, dominant, int(sum(evaluable_mask))

    # Evaluate QP argmax
    qp_acc, qp_f1m, qp_f1w, dominant_qp, n_evaluable = evaluate_assignments(
        Y_assignments, "QP argmax"
    )

    # Beta-weighted gating classifier (uses QP-learned betas + positive/negative evidence)
    logger.info("=" * 60)
    logger.info("Beta-weighted gating classifier")
    Y_gate, _ = beta_weighted_classification(
        marker_level_data=marker_data,
        marker_names=mapped_marker_names,
        assignment_matrix=assignment_matrix,
        cell_type_names=type_names,
        beta_values=beta_values,
    )
    gate_acc, gate_f1m, gate_f1w, dominant_gate, _ = evaluate_assignments(
        Y_gate, "Beta-weighted gating"
    )

    # Use the better result for downstream (Pass 2)
    if gate_acc >= qp_acc:
        logger.info("Using beta-weighted gating for Pass 2 (better accuracy)")
        Y_final = Y_gate
        dominant_type = dominant_gate
        accuracy, f1_macro, f1_weighted = gate_acc, gate_f1m, gate_f1w
        method = "beta_weighted_gating"
    else:
        logger.info("Using QP argmax for Pass 2 (better accuracy)")
        Y_final = Y_assignments
        dominant_type = dominant_qp
        accuracy, f1_macro, f1_weighted = qp_acc, qp_f1m, qp_f1w
        method = "qp_argmax"

    results = {
        "region_id": region_id,
        "n_cells": int(adata_gex.shape[0]),
        "n_cells_with_gt": len(common_cells),
        "n_evaluable": n_evaluable,
        "n_profiles": len(type_names),
        "profiles": type_names,
        "method": method,
        "accuracy": accuracy,
        "f1_macro": f1_macro,
        "f1_weighted": f1_weighted,
        "qp_accuracy": qp_acc,
        "qp_f1_macro": qp_f1m,
        "gate_accuracy": gate_acc,
        "gate_f1_macro": gate_f1m,
        "max_Y_mean": float(Y_final.max(axis=1).mean()),
        "max_Y_median": float(np.median(Y_final.max(axis=1))),
        "doublet_fraction": float(np.mean(Y_final.max(axis=1) < 0.6)),
    }

    # Save assignments
    y_df = pd.DataFrame(Y_final, columns=type_names, index=adata_protein.obs_names)
    y_df.to_csv(output_dir / "cell_assignments.csv")

    # === MODULE 3 PASS 2: True Count Estimation ===
    logger.info("=" * 60)
    logger.info("MODULE 3 PASS 2: True Count Estimation (cell resolution)")

    X_gex = adata_gex.X
    if hasattr(X_gex, "toarray"):
        X_gex = X_gex.toarray()
    X_gex = np.asarray(X_gex, dtype=np.float64)

    # Build enrichment weights from Pass 1 assignments
    n_types = len(type_names)
    n_genes = X_gex.shape[1]
    enrichment = np.ones((n_types, n_genes)) * 0.1
    for ct_idx in range(n_types):
        ct_cells = np.where(dominant_type == ct_idx)[0]
        if len(ct_cells) > 0:
            ct_mean = X_gex[ct_cells].mean(axis=0)
            global_mean = X_gex.mean(axis=0) + 1e-10
            enrichment[ct_idx] = ct_mean / global_mean

    X_true = estimate_true_expression_cell(
        X_obs=X_gex,
        Y_assignments=Y_final,
        coords=coords,
        enrichment_weights=enrichment,
        library_slack=1.5,
        lambda_spatial=0.01,
        spatial_k=50,
    )

    # Evaluate dropout recovery
    obs_zeros = np.sum(X_gex == 0)
    true_zeros = np.sum(X_true == 0)
    results["dropout_recovery"] = {
        "obs_zero_fraction": float(obs_zeros / X_gex.size),
        "true_zero_fraction": float(true_zeros / X_true.size),
        "zeros_recovered": int(obs_zeros - true_zeros),
        "recovery_rate": float((obs_zeros - true_zeros) / max(obs_zeros, 1)),
    }

    # Save results
    with open(output_dir / "validation_results.json", "w") as f:
        json.dump(results, f, indent=2, default=str)

    logger.info(f"Results saved to {output_dir}")
    logger.info(f"Accuracy: {results['accuracy']:.3f}")
    logger.info(f"F1 (macro): {results['f1_macro']:.3f}")
    logger.info(f"Max Y mean: {results['max_Y_mean']:.3f}")
    logger.info(f"Doublet fraction: {results['doublet_fraction']:.3f}")
    logger.info(f"Dropout recovery rate: {results['dropout_recovery']['recovery_rate']:.3f}")

    return results


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Validate cell-resolution CITEgeist on Xenium single-cell data."
    )
    parser.add_argument("--region", type=int, default=0, help="Region ID (0-4)")
    parser.add_argument("--max-cells", type=int, default=10000, help="Max cells to load")
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    run_validation(args.region, args.max_cells)
