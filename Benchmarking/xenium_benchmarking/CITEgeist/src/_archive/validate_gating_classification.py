"""
Validate gating-based cell classification (cell_classification.py) on Xenium data.

Compares Module 3 gating classification against protein-gated ground truth
and the previous beta-weighted gating baseline.
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

from model.cell_classification import (
    GatingProfileSet,
    determine_thresholds,
    classify_cells_gating,
    compute_confidence,
    infer_negative_gates,
    generate_threshold_report,
)
from load_xenium_singlecell import load_xenium_singlecell
from benchmark_constants import ACHIEVABLE_7_CELL_PROFILE_DICT

logger = logging.getLogger(__name__)
OUTPUT_DIR = REPO_ROOT / "Benchmarking" / "xenium_benchmarking" / "CITEgeist" / "output_gating_validation"
GT_PATH = REPO_ROOT / "Benchmarking" / "xenium_pseudovisium" / "data_protein_gt" / "cell_type_assignments.csv"


def run_validation(region_id: int, max_cells: int = 10000):
    """Run gating-based classification validation for one region."""
    output_dir = OUTPUT_DIR / f"region_{region_id}"
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load data
    logger.info(f"Loading Xenium region {region_id} (max {max_cells} cells)")
    adata_gex, adata_protein = load_xenium_singlecell(
        region_id=region_id, max_cells=max_cells, seed=42
    )
    logger.info(f"Loaded {adata_gex.shape[0]} cells, {adata_protein.shape[1]} proteins")

    # Load ground truth
    gt_df = pd.read_csv(GT_PATH, index_col=0)
    common_cells = sorted(set(adata_protein.obs_names) & set(gt_df.index))
    gt_labels = gt_df.loc[common_cells, "cell_type"].values
    logger.info(f"Ground truth: {len(common_cells)} / {adata_protein.shape[0]} cells")

    # Prepare protein data
    X_protein = adata_protein.X
    if hasattr(X_protein, "toarray"):
        X_protein = X_protein.toarray()
    X_protein = np.asarray(X_protein, dtype=np.float64)
    marker_names = list(adata_protein.var_names)

    # === GATING CLASSIFICATION ===
    logger.info("=" * 60)
    logger.info("GATING-BASED CLASSIFICATION (cell_classification.py)")
    logger.info("=" * 60)

    # Build profile set from Achievable-7 dict
    profile_set = GatingProfileSet.from_flat_dict(ACHIEVABLE_7_CELL_PROFILE_DICT)
    type_names = profile_set.gating_order

    # Log profiles
    for pname in type_names:
        p = profile_set.profiles[pname]
        logger.info(f"  {pname}: Major={p.major_markers}, Minor={p.minor_markers}, priority={p.priority}")

    # Test with negative gates OFF (default, matches plan)
    logger.info("\n--- Thresholds (BIC-based GMM) ---")
    spatial_coords = adata_protein.obsm.get('spatial', None)
    thresholds = determine_thresholds(
        protein_data=X_protein,
        marker_names=marker_names,
        method="auto",
        spatial_coords=spatial_coords,
    )

    # Save threshold report
    report_df = generate_threshold_report(
        thresholds=thresholds,
        protein_data=X_protein,
        marker_names=marker_names,
        output_dir=str(output_dir / "threshold_report"),
    )

    # Classify without negative gates
    logger.info("\n--- Classification: No negative gates ---")
    result_no_neg = classify_cells_gating(
        protein_data=X_protein,
        marker_names=marker_names,
        profile_set=profile_set,
        thresholds=thresholds,
        use_negative_gates=False,
    )

    # Also test with negative gates ON for comparison
    logger.info("\n--- Classification: With inferred negative gates ---")
    infer_negative_gates(profile_set.profiles)
    for pname in type_names:
        p = profile_set.profiles[pname]
        if p.inferred_negatives:
            logger.info(f"  {pname}: inferred negatives = {p.inferred_negatives}")

    result_with_neg = classify_cells_gating(
        protein_data=X_protein,
        marker_names=marker_names,
        profile_set=profile_set,
        thresholds=thresholds,
        use_negative_gates=True,
    )

    # Compute confidence for both
    compute_confidence(X_protein, marker_names, result_no_neg, profile_set)
    compute_confidence(X_protein, marker_names, result_with_neg, profile_set)

    # === EVALUATION ===
    from sklearn.metrics import accuracy_score, f1_score, classification_report

    def evaluate(result, label):
        """Evaluate against GT."""
        cell_to_pred = {}
        for i, cell_id in enumerate(adata_protein.obs_names):
            cell_to_pred[cell_id] = result.cell_type_names[result.assignments[i]]

        pred_labels = [cell_to_pred[c] for c in common_cells]
        gt_types_in_profiles = set(type_names)
        evaluable = [(gt in gt_types_in_profiles) and (pred != "Unassigned")
                     for gt, pred in zip(gt_labels, pred_labels)]

        gt_eval = gt_labels[evaluable]
        pred_eval = np.array(pred_labels)[evaluable]

        n_unassigned_pred = sum(1 for p in pred_labels if p == "Unassigned")
        n_gt_not_in_profiles = sum(1 for gt in gt_labels if gt not in gt_types_in_profiles)

        acc = float(accuracy_score(gt_eval, pred_eval)) if len(gt_eval) > 0 else 0.0
        f1_mac = float(f1_score(gt_eval, pred_eval, average="macro", zero_division=0))
        f1_wt = float(f1_score(gt_eval, pred_eval, average="weighted", zero_division=0))

        logger.info(f"\n--- {label} ---")
        logger.info(f"Evaluable cells: {len(gt_eval)} / {len(common_cells)}")
        logger.info(f"  Unassigned (pred): {n_unassigned_pred}")
        logger.info(f"  GT not in profiles: {n_gt_not_in_profiles}")
        logger.info(f"Accuracy: {acc:.3f}, F1 macro: {f1_mac:.3f}, F1 weighted: {f1_wt:.3f}")
        logger.info("\n" + classification_report(gt_eval, pred_eval, zero_division=0))

        # Also compute accuracy including Unassigned as wrong
        evaluable_gt_only = [gt in gt_types_in_profiles for gt in gt_labels]
        gt_eval_all = gt_labels[evaluable_gt_only]
        pred_eval_all = np.array(pred_labels)[evaluable_gt_only]
        acc_incl_unassigned = float(accuracy_score(gt_eval_all, pred_eval_all)) if len(gt_eval_all) > 0 else 0.0
        logger.info(f"Accuracy (including Unassigned as wrong): {acc_incl_unassigned:.3f}")

        return {
            "accuracy": acc,
            "accuracy_incl_unassigned": acc_incl_unassigned,
            "f1_macro": f1_mac,
            "f1_weighted": f1_wt,
            "n_evaluable": len(gt_eval),
            "n_unassigned_pred": n_unassigned_pred,
            "pct_unassigned": 100.0 * n_unassigned_pred / len(common_cells),
        }

    results_no_neg = evaluate(result_no_neg, "Gating (no negative gates)")
    results_with_neg = evaluate(result_with_neg, "Gating (with negative gates)")

    # Per-type distribution
    logger.info("\n--- Cell type distribution ---")
    for label, result in [("No neg gates", result_no_neg), ("With neg gates", result_with_neg)]:
        logger.info(f"\n{label}:")
        for i, name in enumerate(result.cell_type_names):
            count = int((result.assignments == i).sum())
            pct = 100.0 * count / len(result.assignments)
            logger.info(f"  {name}: {count:,} ({pct:.1f}%)")

    # Doublet analysis
    n_doublets_no = int(result_no_neg.doublet_flags.sum())
    n_doublets_with = int(result_with_neg.doublet_flags.sum())
    logger.info(f"\nDoublets (no neg): {n_doublets_no} ({100*n_doublets_no/len(result_no_neg.assignments):.1f}%)")
    logger.info(f"Doublets (with neg): {n_doublets_with} ({100*n_doublets_with/len(result_with_neg.assignments):.1f}%)")

    # Confidence summary
    assigned_no = result_no_neg.assignments != result_no_neg.cell_type_names.index("Unassigned")
    assigned_with = result_with_neg.assignments != result_with_neg.cell_type_names.index("Unassigned")
    if assigned_no.any():
        logger.info(f"\nConfidence (no neg, assigned): mean={result_no_neg.confidence[assigned_no].mean():.3f}")
    if assigned_with.any():
        logger.info(f"Confidence (with neg, assigned): mean={result_with_neg.confidence[assigned_with].mean():.3f}")

    # Save results
    output = {
        "region_id": region_id,
        "n_cells": int(adata_protein.shape[0]),
        "n_cells_with_gt": len(common_cells),
        "n_proteins": int(adata_protein.shape[1]),
        "protein_names": marker_names,
        "profiles": type_names,
        "gating_no_neg": results_no_neg,
        "gating_with_neg": results_with_neg,
        "n_doublets_no_neg": n_doublets_no,
        "n_doublets_with_neg": n_doublets_with,
    }

    with open(output_dir / "validation_results.json", "w") as f:
        json.dump(output, f, indent=2, default=str)

    # Save thresholds
    thresholds.save(str(output_dir / "thresholds.json"))

    # Save assignments
    result_no_neg.to_dataframe(list(adata_protein.obs_names)).to_csv(
        output_dir / "assignments_no_neg.csv"
    )
    result_with_neg.to_dataframe(list(adata_protein.obs_names)).to_csv(
        output_dir / "assignments_with_neg.csv"
    )

    logger.info(f"\nResults saved to {output_dir}")
    logger.info(f"Best accuracy (no neg): {results_no_neg['accuracy']:.3f}")
    logger.info(f"Best accuracy (with neg): {results_with_neg['accuracy']:.3f}")

    return output


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Validate gating-based cell classification on Xenium data."
    )
    parser.add_argument("--region", type=int, default=2, help="Region ID (0-4)")
    parser.add_argument("--max-cells", type=int, default=10000, help="Max cells to load")
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    run_validation(args.region, args.max_cells)
