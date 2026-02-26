#!/usr/bin/env python3
"""
Evaluate VAE + Sinkhorn single-cell assignment against Xenium ground truth.

This script compares cell type predictions from the VAE + Sinkhorn module
against Xenium-derived ground truth labels.

Usage:
    python evaluate_vae_assignment.py \
        --predictions path/to/predictions.csv \
        --ground_truth path/to/ground_truth.csv \
        --output_dir path/to/output/

Inputs:
    --predictions: CSV with columns: nucleus_id, predicted_type, [confidence]
    --ground_truth: CSV with columns: nucleus_id, gt_type

Outputs:
    - evaluation_results.json: accuracy, n_samples, classification_report dict
    - confusion_matrix.csv: pandas DataFrame with labels as index/columns
    - confidence_calibration.csv (if confidence column exists)
"""

import argparse
import json
import logging
from pathlib import Path
from typing import Dict, Optional

import numpy as np
import pandas as pd
from sklearn.metrics import accuracy_score, classification_report, confusion_matrix

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)


def evaluate_assignment(
    predictions_file: str,
    ground_truth_file: str,
    output_dir: str,
) -> Dict:
    """
    Evaluate VAE + Sinkhorn assignment accuracy against ground truth.

    Args:
        predictions_file: CSV with nucleus_id, predicted_type, confidence columns
        ground_truth_file: CSV with nucleus_id, gt_type (from Xenium)
        output_dir: Directory for results

    Returns:
        Dict with evaluation results
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load data
    logger.info(f"Loading predictions from {predictions_file}")
    pred_df = pd.read_csv(predictions_file)

    logger.info(f"Loading ground truth from {ground_truth_file}")
    gt_df = pd.read_csv(ground_truth_file)

    # Validate required columns
    if "nucleus_id" not in pred_df.columns:
        raise ValueError("Predictions file must have 'nucleus_id' column")
    if "predicted_type" not in pred_df.columns:
        raise ValueError("Predictions file must have 'predicted_type' column")
    if "nucleus_id" not in gt_df.columns:
        raise ValueError("Ground truth file must have 'nucleus_id' column")
    if "gt_type" not in gt_df.columns:
        raise ValueError("Ground truth file must have 'gt_type' column")

    logger.info(f"Predictions: {len(pred_df)} cells")
    logger.info(f"Ground truth: {len(gt_df)} cells")

    # Merge on nucleus_id
    merged = pred_df.merge(gt_df, on="nucleus_id", how="inner")
    logger.info(f"Matched cells: {len(merged)}")

    if len(merged) == 0:
        logger.error("No matching nucleus_ids between predictions and ground truth")
        results = {
            "error": "no_matching_cells",
            "n_predictions": len(pred_df),
            "n_ground_truth": len(gt_df),
        }
        with open(output_dir / "evaluation_results.json", "w") as f:
            json.dump(results, f, indent=2)
        return results

    y_true = merged["gt_type"].values
    y_pred = merged["predicted_type"].values

    # Overall accuracy
    accuracy = accuracy_score(y_true, y_pred)
    logger.info(f"Overall accuracy: {accuracy:.4f}")

    # Per-class report
    report_dict = classification_report(y_true, y_pred, output_dict=True, zero_division=0)
    report_str = classification_report(y_true, y_pred, zero_division=0)
    print("\nClassification Report:")
    print(report_str)

    # Confusion matrix
    labels = sorted(set(y_true) | set(y_pred))
    cm = confusion_matrix(y_true, y_pred, labels=labels)
    cm_df = pd.DataFrame(cm, index=labels, columns=labels)

    # Save confusion matrix
    cm_df.to_csv(output_dir / "confusion_matrix.csv")
    logger.info(f"Saved confusion matrix to {output_dir / 'confusion_matrix.csv'}")

    # Build results dict
    results = {
        "accuracy": float(accuracy),
        "n_samples": len(merged),
        "n_predictions_total": len(pred_df),
        "n_ground_truth_total": len(gt_df),
        "n_unmatched_predictions": len(pred_df) - len(merged),
        "n_unmatched_ground_truth": len(gt_df) - len(merged),
        "classification_report": report_dict,
        "labels": labels,
    }

    # Per-class summary
    per_class_summary = {}
    for label in labels:
        if label in report_dict:
            per_class_summary[label] = {
                "precision": report_dict[label]["precision"],
                "recall": report_dict[label]["recall"],
                "f1-score": report_dict[label]["f1-score"],
                "support": report_dict[label]["support"],
            }
    results["per_class_summary"] = per_class_summary

    # Confidence calibration (if confidence column exists)
    if "confidence" in merged.columns:
        logger.info("Computing confidence calibration...")
        conf_calibration = compute_confidence_calibration(merged)
        results["confidence_calibration"] = conf_calibration

        # Save calibration data
        conf_df = pd.DataFrame(conf_calibration["bins"])
        conf_df.to_csv(output_dir / "confidence_calibration.csv", index=False)
        logger.info(f"Saved confidence calibration to {output_dir / 'confidence_calibration.csv'}")

        print("\nConfidence Calibration:")
        for bin_info in conf_calibration["bins"]:
            print(f"  {bin_info['bin']}: accuracy={bin_info['accuracy']:.4f} ({bin_info['n_samples']} samples)")
    else:
        logger.info("No confidence column found, skipping calibration")

    # Save results
    with open(output_dir / "evaluation_results.json", "w") as f:
        json.dump(results, f, indent=2, default=lambda x: float(x) if isinstance(x, np.floating) else x)
    logger.info(f"Saved evaluation results to {output_dir / 'evaluation_results.json'}")

    # Print summary
    print("\n" + "=" * 60)
    print("VAE ASSIGNMENT EVALUATION SUMMARY")
    print("=" * 60)
    print(f"Overall Accuracy: {accuracy:.4f}")
    print(f"Samples Evaluated: {len(merged)}")
    print(f"Number of Classes: {len(labels)}")
    print("=" * 60)

    return results


def compute_confidence_calibration(merged_df: pd.DataFrame) -> Dict:
    """
    Compute accuracy per confidence bin.

    Args:
        merged_df: DataFrame with 'confidence', 'predicted_type', 'gt_type' columns

    Returns:
        Dict with calibration data
    """
    bins = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
    bin_labels = ["0.0-0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1.0"]

    # Ensure confidence is numeric
    merged_df = merged_df.copy()
    merged_df["confidence"] = pd.to_numeric(merged_df["confidence"], errors="coerce")

    # Remove rows with NaN confidence
    valid_mask = ~merged_df["confidence"].isna()
    n_invalid = (~valid_mask).sum()
    if n_invalid > 0:
        logger.warning(f"Skipping {n_invalid} rows with invalid confidence values")
    merged_df = merged_df[valid_mask]

    # Bin the confidence values
    merged_df["conf_bin"] = pd.cut(
        merged_df["confidence"],
        bins=bins,
        labels=bin_labels,
        include_lowest=True,
    )

    # Compute accuracy per bin
    calibration_data = []
    for bin_label in bin_labels:
        bin_mask = merged_df["conf_bin"] == bin_label
        bin_df = merged_df[bin_mask]
        n_samples = len(bin_df)

        if n_samples > 0:
            correct = (bin_df["predicted_type"] == bin_df["gt_type"]).sum()
            acc = correct / n_samples
        else:
            acc = np.nan

        calibration_data.append({
            "bin": bin_label,
            "accuracy": float(acc) if not np.isnan(acc) else None,
            "n_samples": int(n_samples),
            "n_correct": int(correct) if n_samples > 0 else 0,
        })

    # Compute expected calibration error (ECE)
    total_samples = len(merged_df)
    ece = 0.0
    for bin_info in calibration_data:
        if bin_info["accuracy"] is not None and bin_info["n_samples"] > 0:
            # Expected confidence for this bin (midpoint)
            bin_start, bin_end = map(float, bin_info["bin"].split("-"))
            expected_conf = (bin_start + bin_end) / 2
            # Weight by proportion of samples in this bin
            weight = bin_info["n_samples"] / total_samples
            ece += weight * abs(bin_info["accuracy"] - expected_conf)

    return {
        "bins": calibration_data,
        "expected_calibration_error": float(ece),
        "n_total_samples": int(total_samples),
    }


def main():
    parser = argparse.ArgumentParser(
        description="Evaluate VAE + Sinkhorn single-cell assignment"
    )
    parser.add_argument(
        "--predictions",
        type=str,
        required=True,
        help="CSV with nucleus_id, predicted_type, [confidence] columns",
    )
    parser.add_argument(
        "--ground_truth",
        type=str,
        required=True,
        help="CSV with nucleus_id, gt_type columns",
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        required=True,
        help="Directory for evaluation results",
    )

    args = parser.parse_args()

    evaluate_assignment(
        predictions_file=args.predictions,
        ground_truth_file=args.ground_truth,
        output_dir=args.output_dir,
    )


if __name__ == "__main__":
    main()
