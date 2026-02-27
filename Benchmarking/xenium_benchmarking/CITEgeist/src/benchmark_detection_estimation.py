#!/usr/bin/env python
"""
Benchmark detection + estimation model on Xenium pseudo-Visium data.

This script runs the two-stage detection + estimation model:
- Stage 1: Multivariate GMM detection per cell type
- Stage 2: Masked IQP with learned alpha, beta, sigma^2

The goal is to compare against the existing continuous model, particularly
for B cells sparsity (current: 43% non-zero, target: <5% for rare types).
"""
import argparse
import json
import logging
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import pearsonr

# Add project root to path
REPO_ROOT = Path(__file__).parents[4]
sys.path.insert(0, str(REPO_ROOT))

# Add benchmark root for shared constants
BENCHMARK_ROOT = Path(__file__).parents[2]
sys.path.insert(0, str(BENCHMARK_ROOT))

from benchmark_constants import ACHIEVABLE_7_CELL_PROFILE_DICT

from CITEgeist.model.detection_estimation import solve_detection_estimation

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)


def build_profile_and_marker_groups(
    cell_profile_dict: dict,
    antibody_names: list,
) -> tuple:
    """
    Build profile matrix and marker_groups from cell_profile_dict.

    Args:
        cell_profile_dict: Dict mapping cell_type -> {"Major": [...], "Minor": [...]}
        antibody_names: List of antibody names in order

    Returns:
        profile: (n_types, n_markers) binary matrix
        marker_groups: Dict mapping cell_type -> marker indices
    """
    cell_types = list(cell_profile_dict.keys())
    n_types = len(cell_types)
    n_markers = len(antibody_names)

    # Build name to index mapping
    name_to_idx = {name: i for i, name in enumerate(antibody_names)}

    profile = np.zeros((n_types, n_markers))
    marker_groups = {}

    for k, cell_type in enumerate(cell_types):
        markers = cell_profile_dict[cell_type]["Major"]
        # Also include minor markers
        markers = markers + cell_profile_dict[cell_type].get("Minor", [])

        indices = []
        for marker in markers:
            if marker in name_to_idx:
                idx = name_to_idx[marker]
                profile[k, idx] = 1
                indices.append(idx)
            else:
                logger.warning(f"Marker {marker} not found in antibody panel")

        marker_groups[cell_type] = indices

    return profile, marker_groups


def compute_proportions(counts: np.ndarray) -> np.ndarray:
    """Convert integer counts to proportions."""
    sums = counts.sum(axis=1, keepdims=True)
    sums = np.maximum(sums, 1)  # Avoid division by zero
    return counts / sums


def evaluate_vs_ground_truth(
    pred_counts: np.ndarray,
    gt_df: pd.DataFrame,
    cell_types: list,
    detected: np.ndarray,
) -> dict:
    """
    Evaluate predicted counts against ground truth.

    Args:
        pred_counts: (n_spots, n_types) predicted integer counts
        gt_df: DataFrame with ground truth counts per cell type
        cell_types: List of cell type names
        detected: (n_spots, n_types) detection mask

    Returns:
        Dict with evaluation metrics
    """
    n_spots, n_types = pred_counts.shape

    # Compute proportions
    pred_props = compute_proportions(pred_counts)

    # Map predicted types to GT columns
    # GT columns are like "B cells", "CD4+ T cells", etc.
    gt_cols = [c for c in gt_df.columns if c in cell_types]

    if len(gt_cols) == 0:
        logger.warning("No matching GT columns found")
        return {}

    gt_values = gt_df[gt_cols].values

    # Check if GT is already proportions (values sum to ~1) or counts
    row_sums = gt_values.sum(axis=1)
    if np.allclose(row_sums, 1.0, atol=0.1):
        # Already proportions
        gt_props = gt_values
    else:
        # Counts - normalize
        gt_props = gt_values / np.maximum(row_sums[:, np.newaxis], 1e-8)

    # Overall correlation
    pred_flat = []
    gt_flat = []

    for k, cell_type in enumerate(cell_types):
        if cell_type in gt_cols:
            gt_idx = gt_cols.index(cell_type)
            pred_flat.extend(pred_props[:, k])
            gt_flat.extend(gt_props[:, gt_idx])

    overall_r, _ = pearsonr(pred_flat, gt_flat)

    # Per-type metrics
    per_type_metrics = {}
    for k, cell_type in enumerate(cell_types):
        if cell_type not in gt_cols:
            continue

        gt_idx = gt_cols.index(cell_type)

        # Proportion correlation
        r, _ = pearsonr(pred_props[:, k], gt_props[:, gt_idx])

        # Detection rate
        n_detected = detected[:, k].sum()
        detection_rate = n_detected / n_spots

        # GT presence rate (how often is this type actually present?)
        gt_present = (gt_props[:, gt_idx] > 0).sum() / n_spots

        # Sparsity analysis (key for B cells issue)
        pred_nonzero = (pred_counts[:, k] > 0).sum() / n_spots
        gt_nonzero = gt_present

        per_type_metrics[cell_type] = {
            "pearson_r": float(r),
            "detection_rate": float(detection_rate),
            "gt_presence_rate": float(gt_present),
            "pred_nonzero_rate": float(pred_nonzero),
            "gt_nonzero_rate": float(gt_nonzero),
        }

    return {
        "overall_pearson_r": float(overall_r),
        "per_type": per_type_metrics,
    }


def main():
    parser = argparse.ArgumentParser(
        description="Benchmark detection + estimation model on Xenium pseudo-Visium"
    )
    parser.add_argument("--region-id", type=int, default=0, help="Region ID (0-4)")
    parser.add_argument(
        "--input-dir",
        type=str,
        default=str(REPO_ROOT / "Benchmarking/xenium_pseudovisium/data_protein_gt"),
        help="Input directory with h5ad_objects/",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default=str(REPO_ROOT / "Benchmarking/xenium_benchmarking/CITEgeist/output/detection_estimation"),
        help="Output directory",
    )
    parser.add_argument("--max-iter", type=int, default=10, help="Max EM iterations")
    parser.add_argument(
        "--detection-threshold", type=float, default=0.5, help="GMM detection threshold"
    )

    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info(f"Running detection + estimation benchmark on region {args.region_id}")

    # Load data
    gex_path = input_dir / "h5ad_objects" / f"Xenium_region_{args.region_id}_GEX.h5ad"
    cite_path = input_dir / "h5ad_objects" / f"Xenium_region_{args.region_id}_CITE.h5ad"
    gt_path = input_dir / "ground_truth" / f"Xenium_region_{args.region_id}_prop.csv"

    for path in [gex_path, cite_path]:
        if not path.exists():
            logger.error(f"File not found: {path}")
            sys.exit(1)

    logger.info("Loading data...")
    gex_adata = sc.read_h5ad(gex_path)
    cite_adata = sc.read_h5ad(cite_path)

    # Load ground truth if available
    gt_df = None
    if gt_path.exists():
        gt_df = pd.read_csv(gt_path, index_col=0)
        logger.info(f"Loaded GT with columns: {list(gt_df.columns)}")

    # Extract data
    X = cite_adata.X
    if hasattr(X, "toarray"):
        X = X.toarray()
    X = np.asarray(X, dtype=float)

    nuclei_counts = cite_adata.obs["n_cells"].values.astype(int)
    antibody_names = list(cite_adata.var_names)

    logger.info(f"Data shape: {X.shape}, nuclei range: {nuclei_counts.min()}-{nuclei_counts.max()}")

    # Build profile and marker groups
    profile, marker_groups = build_profile_and_marker_groups(
        ACHIEVABLE_7_CELL_PROFILE_DICT, antibody_names
    )
    cell_types = list(ACHIEVABLE_7_CELL_PROFILE_DICT.keys())

    logger.info(f"Cell types: {cell_types}")
    logger.info(f"Marker groups: {marker_groups}")

    # Run detection + estimation
    logger.info("Running detection + estimation pipeline...")
    start_time = time.time()

    detected, counts, alpha, beta, sigma_sq = solve_detection_estimation(
        X=X,
        nuclei_counts=nuclei_counts,
        profile=profile,
        marker_groups=marker_groups,
        max_iter=args.max_iter,
        detection_threshold=args.detection_threshold,
    )

    runtime = time.time() - start_time
    logger.info(f"Completed in {runtime:.1f}s")

    # Log learned parameters
    logger.info(f"Learned alpha (baseline): {alpha}")
    logger.info(f"Learned beta (signal/cell): {beta}")
    logger.info(f"Learned sigma: {np.sqrt(sigma_sq)}")

    # Evaluate if GT available
    metrics = {}
    if gt_df is not None:
        metrics = evaluate_vs_ground_truth(counts, gt_df, cell_types, detected)
        logger.info(f"Overall Pearson r: {metrics['overall_pearson_r']:.4f}")

        for ct, m in metrics.get("per_type", {}).items():
            logger.info(
                f"  {ct}: r={m['pearson_r']:.3f}, "
                f"detection={m['detection_rate']:.1%}, "
                f"pred_nonzero={m['pred_nonzero_rate']:.1%}, "
                f"gt_nonzero={m['gt_nonzero_rate']:.1%}"
            )

    # Save results
    result_dir = output_dir / f"region_{args.region_id}"
    result_dir.mkdir(parents=True, exist_ok=True)

    # Save counts
    counts_df = pd.DataFrame(counts, columns=cell_types, index=cite_adata.obs_names)
    counts_df.to_csv(result_dir / "predicted_counts.csv")

    # Save detection mask
    detected_df = pd.DataFrame(detected, columns=cell_types, index=cite_adata.obs_names)
    detected_df.to_csv(result_dir / "detection_mask.csv")

    # Save summary
    results = {
        "region_id": args.region_id,
        "n_spots": X.shape[0],
        "n_markers": X.shape[1],
        "n_types": len(cell_types),
        "cell_types": cell_types,
        "runtime_sec": runtime,
        "parameters": {
            "max_iter": args.max_iter,
            "detection_threshold": args.detection_threshold,
        },
        "learned_alpha": alpha.tolist(),
        "learned_beta": beta.tolist(),
        "learned_sigma": np.sqrt(sigma_sq).tolist(),
        "metrics": metrics,
    }

    with open(result_dir / "results.json", "w") as f:
        json.dump(results, f, indent=2)

    logger.info(f"Results saved to {result_dir}")


if __name__ == "__main__":
    main()
