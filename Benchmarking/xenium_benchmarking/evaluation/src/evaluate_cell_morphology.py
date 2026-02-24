#!/usr/bin/env python
"""
Evaluate cell morphology features for cell type classification.

Compares:
1. Nuclear features only (current baseline)
2. Cell features (from watershed)
3. Combined nuclear + cell features

Reports accuracy, macro F1, and per-class metrics.

Usage:
    python evaluate_cell_morphology.py --region 0 --output-dir ./output/eval_cell_morph
"""

import argparse
import json
import logging
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, f1_score, classification_report
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler

# Add CITEgeist to path
REPO_ROOT = Path(__file__).parent.parent.parent.parent.parent
sys.path.insert(0, str(REPO_ROOT))

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)

# =============================================================================
# DEFAULT PATHS
# =============================================================================

XENIUM_DATA_DIR = Path(
    "/ix1/alee/LO_LAB/General/Public_Data/10x_PublicData/"
    "Xenium_RNA_Proteomic_RenalCellCarcinoma"
)
PSEUDOVISIUM_DIR = REPO_ROOT / "Benchmarking/xenium_pseudovisium"

# Pixel size for Xenium morphology images
PIXEL_SIZE_UM = 0.2125

# Region boundaries in microns (from pseudo-Visium data)
REGION_BOUNDS_UM = {
    0: (29.01, 2279.01, 30.87, 5486.83),
    1: (2329.01, 4579.01, 30.87, 5400.23),
    2: (4629.01, 6829.01, 30.87, 5486.83),
    3: (6879.01, 9129.01, 30.87, 5660.04),
    4: (9179.01, 11429.01, 30.87, 5746.64),
}
PADDING_UM = 100.0


# =============================================================================
# DATA LOADING FUNCTIONS
# =============================================================================


def load_ground_truth_cells(region_id: int) -> pd.DataFrame:
    """
    Load ground truth cell types and coordinates for a region.

    Loads the protein-gated cell type assignments and joins with
    Xenium cell coordinates. Filters to cells within the region.

    Args:
        region_id: Region ID (0-4)

    Returns:
        DataFrame with columns:
            - cell_id (index)
            - x_centroid, y_centroid: micron coordinates
            - cell_type: protein-gated ground truth type
    """
    # Load cell type assignments
    gt_path = PSEUDOVISIUM_DIR / "data_protein_gt" / "cell_type_assignments.csv"
    if not gt_path.exists():
        raise FileNotFoundError(f"Ground truth not found: {gt_path}")

    gt_df = pd.read_csv(gt_path)
    # Handle unnamed index column
    if gt_df.columns[0] in ("Unnamed: 0", ""):
        gt_df = gt_df.rename(columns={gt_df.columns[0]: "cell_id"})
    gt_df = gt_df.set_index("cell_id")

    # Load cell coordinates from Xenium
    cells_path = XENIUM_DATA_DIR / "cells.parquet"
    if not cells_path.exists():
        raise FileNotFoundError(f"Xenium cells not found: {cells_path}")

    cells_df = pd.read_parquet(cells_path)
    cells_df = cells_df.set_index("cell_id")

    # Join cell types with coordinates
    merged = gt_df.join(cells_df[["x_centroid", "y_centroid"]], how="inner")

    # Filter to cells in this region
    x_min, x_max, y_min, y_max = REGION_BOUNDS_UM[region_id]
    x_min -= PADDING_UM
    x_max += PADDING_UM
    y_min -= PADDING_UM
    y_max += PADDING_UM

    region_mask = (
        (merged["x_centroid"] >= x_min) &
        (merged["x_centroid"] <= x_max) &
        (merged["y_centroid"] >= y_min) &
        (merged["y_centroid"] <= y_max)
    )
    region_cells = merged[region_mask].copy()

    # Exclude "Unknown" cell types
    region_cells = region_cells[region_cells["cell_type"] != "Unknown"]

    logger.info(
        "Loaded %d ground truth cells for region %d",
        len(region_cells), region_id
    )

    return region_cells


def load_cell_features(features_path: Path) -> pd.DataFrame:
    """
    Load cell morphology features from benchmark output.

    Args:
        features_path: Path to cell_features.csv

    Returns:
        DataFrame with morphology features and cell coordinates
    """
    if not features_path.exists():
        raise FileNotFoundError(f"Features file not found: {features_path}")

    features_df = pd.read_csv(features_path)
    logger.info("Loaded %d cells with features", len(features_df))

    return features_df


def match_cells_to_gt(
    cell_features: pd.DataFrame,
    gt_df: pd.DataFrame,
    region_id: int,
    max_dist: float = 10.0,
) -> pd.DataFrame:
    """
    Match extracted cells to ground truth by spatial proximity.

    Converts Cellpose pixel coordinates to microns, then uses KDTree
    to find nearest ground truth cell within max_dist.

    Args:
        cell_features: DataFrame from benchmark with centroid_x, centroid_y in pixels
        gt_df: Ground truth DataFrame with x_centroid, y_centroid in microns
        region_id: Region ID for coordinate offset
        max_dist: Maximum matching distance in microns

    Returns:
        DataFrame with matched cells, including gt_cell_type column
    """
    # Convert Cellpose pixel coords to microns
    x_min_um, _, y_min_um, _ = REGION_BOUNDS_UM[region_id]
    x_offset = x_min_um - PADDING_UM
    y_offset = y_min_um - PADDING_UM

    cell_x_um = cell_features["centroid_x"].values * PIXEL_SIZE_UM + x_offset
    cell_y_um = cell_features["centroid_y"].values * PIXEL_SIZE_UM + y_offset

    # Build KDTree from GT coordinates
    gt_coords = gt_df[["x_centroid", "y_centroid"]].values
    tree = cKDTree(gt_coords)

    # Query for each cell
    cell_coords = np.column_stack([cell_x_um, cell_y_um])
    dists, indices = tree.query(cell_coords, k=1)

    # Filter by max distance
    valid_mask = dists <= max_dist

    matched = cell_features[valid_mask].copy()
    matched["gt_cell_type"] = gt_df.iloc[indices[valid_mask]]["cell_type"].values
    matched["match_dist"] = dists[valid_mask]
    matched["gt_cell_id"] = gt_df.index[indices[valid_mask]].values

    logger.info(
        "Matched %d/%d cells (%.1f%%) within %.1f microns",
        valid_mask.sum(), len(cell_features),
        100 * valid_mask.mean(), max_dist
    )

    return matched


# =============================================================================
# EVALUATION FUNCTIONS
# =============================================================================


def evaluate_feature_set(
    df: pd.DataFrame,
    feature_cols: List[str],
    target_col: str = "gt_cell_type",
    n_splits: int = 5,
) -> Dict[str, float]:
    """
    Evaluate classification accuracy using stratified cross-validation.

    Uses logistic regression as a simple classifier to measure how well
    morphology features can predict cell type.

    Args:
        df: DataFrame with features and target column
        feature_cols: List of feature column names to use
        target_col: Name of target column (ground truth cell type)
        n_splits: Number of cross-validation folds

    Returns:
        Dict with accuracy and F1 metrics (mean and std across folds)
    """
    # Check which columns exist
    available = [c for c in feature_cols if c in df.columns]
    if len(available) < len(feature_cols):
        logger.warning(
            "Feature set: only %d/%d columns available: %s",
            len(available), len(feature_cols), available
        )

    if len(available) == 0:
        return {"error": "no_features_available"}

    X = df[available].values
    y = df[target_col].values

    # Handle NaN values
    X = np.nan_to_num(X, nan=0.0, posinf=0.0, neginf=0.0)

    # Check for degenerate cases
    unique_classes = np.unique(y)
    if len(unique_classes) < 2:
        return {"error": "only_one_class"}

    # Standardize features
    scaler = StandardScaler()

    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=42)

    accuracies = []
    f1_scores = []

    for train_idx, test_idx in skf.split(X, y):
        X_train, X_test = X[train_idx], X[test_idx]
        y_train, y_test = y[train_idx], y[test_idx]

        X_train_scaled = scaler.fit_transform(X_train)
        X_test_scaled = scaler.transform(X_test)

        clf = LogisticRegression(
            multi_class="multinomial",
            solver="lbfgs",
            max_iter=1000,
            random_state=42,
        )
        clf.fit(X_train_scaled, y_train)

        y_pred = clf.predict(X_test_scaled)

        accuracies.append(accuracy_score(y_test, y_pred))
        f1_scores.append(f1_score(y_test, y_pred, average="macro"))

    return {
        "n_features": len(available),
        "features_used": available,
        "accuracy_mean": float(np.mean(accuracies)),
        "accuracy_std": float(np.std(accuracies)),
        "f1_macro_mean": float(np.mean(f1_scores)),
        "f1_macro_std": float(np.std(f1_scores)),
    }


def run_evaluation(
    region_id: int,
    cell_features_path: Path,
    output_dir: Path,
    max_dist: float = 10.0,
) -> Dict:
    """
    Run full evaluation comparing feature sets for a region.

    Args:
        region_id: Region ID (0-4)
        cell_features_path: Path to cell_features.csv from benchmark
        output_dir: Output directory for results

    Returns:
        Dict with evaluation results for all feature sets
    """
    region_name = f"Xenium_region_{region_id}"
    logger.info("=" * 60)
    logger.info("EVALUATING %s", region_name)
    logger.info("=" * 60)

    region_out = Path(output_dir) / region_name
    region_out.mkdir(parents=True, exist_ok=True)

    # Load cell features
    cell_features = load_cell_features(cell_features_path)

    # Load ground truth
    gt_df = load_ground_truth_cells(region_id)

    # Match cells to GT
    matched = match_cells_to_gt(cell_features, gt_df, region_id, max_dist)

    if len(matched) < 100:
        logger.warning("Too few matched cells (%d) for reliable evaluation", len(matched))
        return {
            "region": region_name,
            "status": "insufficient_data",
            "n_cells_matched": len(matched),
        }

    # Log cell type distribution
    ct_counts = matched["gt_cell_type"].value_counts()
    logger.info("Cell type distribution:")
    for ct, count in ct_counts.items():
        logger.info("  %s: %d (%.1f%%)", ct, count, 100 * count / len(matched))

    # Define feature sets
    nuclear_features = [
        "nucleus_area",
        "nucleus_circularity",
        "nucleus_eccentricity",
        "nucleus_solidity",
        "nucleus_aspect_ratio",
    ]

    cell_features_cols = [
        "cell_area",
        "cell_circularity",
        "cell_eccentricity",
        "cell_solidity",
        "cell_aspect_ratio",
    ]

    ratio_features = ["nc_ratio", "cytoplasm_area"]

    all_features = nuclear_features + cell_features_cols + ratio_features

    # Evaluate each feature set
    results = {
        "region": region_name,
        "region_id": region_id,
        "status": "success",
        "n_cells_matched": len(matched),
        "n_cell_types": matched["gt_cell_type"].nunique(),
        "cell_type_counts": ct_counts.to_dict(),
        "feature_sets": {},
    }

    feature_sets = {
        "nuclear_only": nuclear_features,
        "cell_only": cell_features_cols,
        "ratio_only": ratio_features,
        "nuclear_plus_ratio": nuclear_features + ratio_features,
        "all_features": all_features,
    }

    for set_name, cols in feature_sets.items():
        logger.info("-" * 40)
        logger.info("Evaluating %s...", set_name)

        metrics = evaluate_feature_set(matched, cols)
        results["feature_sets"][set_name] = metrics

        if "error" not in metrics:
            logger.info(
                "  Accuracy: %.3f +/- %.3f, F1: %.3f +/- %.3f",
                metrics["accuracy_mean"], metrics["accuracy_std"],
                metrics["f1_macro_mean"], metrics["f1_macro_std"]
            )

    # Save results
    with open(region_out / "evaluation_results.json", "w") as f:
        json.dump(results, f, indent=2)

    # Save matched cells for inspection
    matched.to_csv(region_out / "matched_cells.csv", index=False)

    logger.info("Results saved to %s", region_out)

    return results


def aggregate_results(all_results: List[Dict]) -> Dict:
    """
    Aggregate results across all regions.

    Args:
        all_results: List of per-region result dicts

    Returns:
        Dict with aggregated metrics (mean and std across regions)
    """
    successful = [r for r in all_results if r.get("status") == "success"]

    if not successful:
        return {"error": "no_successful_regions"}

    aggregated = {
        "n_regions": len(successful),
        "total_cells": sum(r["n_cells_matched"] for r in successful),
        "feature_sets": {},
    }

    # Get all feature set names
    set_names = list(successful[0]["feature_sets"].keys())

    for set_name in set_names:
        acc_values = []
        f1_values = []

        for r in successful:
            metrics = r["feature_sets"].get(set_name, {})
            if "accuracy_mean" in metrics:
                acc_values.append(metrics["accuracy_mean"])
            if "f1_macro_mean" in metrics:
                f1_values.append(metrics["f1_macro_mean"])

        if acc_values:
            aggregated["feature_sets"][set_name] = {
                "accuracy_mean": float(np.mean(acc_values)),
                "accuracy_std": float(np.std(acc_values)),
                "f1_macro_mean": float(np.mean(f1_values)) if f1_values else None,
                "f1_macro_std": float(np.std(f1_values)) if f1_values else None,
                "n_regions": len(acc_values),
            }

    return aggregated


def main():
    parser = argparse.ArgumentParser(
        description="Evaluate cell morphology features for cell type classification"
    )
    parser.add_argument(
        "--region", type=int, default=None,
        help="Region ID (0-4). If not specified, runs all regions."
    )
    parser.add_argument(
        "--features-path", type=str, default=None,
        help="Path to cell_features.csv from benchmark. "
             "If not specified, looks in default benchmark output directory."
    )
    parser.add_argument(
        "--features-dir", type=str, default=None,
        help="Directory containing per-region benchmark outputs. "
             "Used when --region is not specified."
    )
    parser.add_argument(
        "--output-dir", type=str,
        default=str(REPO_ROOT / "Benchmarking/xenium_benchmarking/evaluation/output/cell_morphology"),
        help="Output directory for evaluation results"
    )
    parser.add_argument(
        "--max-dist", type=float, default=10.0,
        help="Maximum distance in microns for cell matching (default: 10.0)"
    )

    args = parser.parse_args()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Default features directory
    default_features_dir = (
        REPO_ROOT / "Benchmarking/xenium_benchmarking/CITEgeist/output/cell_morphology"
    )

    all_results = []

    if args.region is not None:
        # Single region mode
        regions = [args.region]
        if args.features_path:
            features_paths = {args.region: Path(args.features_path)}
        else:
            features_dir = Path(args.features_dir) if args.features_dir else default_features_dir
            features_paths = {
                args.region: features_dir / f"Xenium_region_{args.region}" / "cell_features.csv"
            }
    else:
        # All regions mode
        regions = [0, 1, 2, 3, 4]
        features_dir = Path(args.features_dir) if args.features_dir else default_features_dir
        features_paths = {
            r: features_dir / f"Xenium_region_{r}" / "cell_features.csv"
            for r in regions
        }

    for region_id in regions:
        features_path = features_paths[region_id]
        if not features_path.exists():
            logger.warning("Features not found for region %d: %s", region_id, features_path)
            continue

        try:
            result = run_evaluation(
                region_id=region_id,
                cell_features_path=features_path,
                output_dir=output_dir,
                max_dist=args.max_dist,
            )
            all_results.append(result)
        except Exception as e:
            logger.error("Failed to evaluate region %d: %s", region_id, e)
            all_results.append({
                "region": f"Xenium_region_{region_id}",
                "status": "error",
                "error": str(e),
            })

    # Aggregate results if multiple regions
    if len(all_results) > 1:
        aggregated = aggregate_results(all_results)
        with open(output_dir / "aggregated_results.json", "w") as f:
            json.dump(aggregated, f, indent=2)

        # Print summary
        print("\n" + "=" * 60)
        print("CELL MORPHOLOGY FEATURE EVALUATION SUMMARY")
        print("=" * 60)
        print(f"\nRegions evaluated: {aggregated.get('n_regions', 0)}")
        print(f"Total cells: {aggregated.get('total_cells', 0)}")

        print("\nFeature Set Comparison:")
        print("-" * 60)
        for set_name, metrics in aggregated.get("feature_sets", {}).items():
            print(f"\n{set_name}:")
            print(f"  Accuracy: {metrics['accuracy_mean']:.3f} +/- {metrics['accuracy_std']:.3f}")
            if metrics.get("f1_macro_mean") is not None:
                print(f"  Macro F1: {metrics['f1_macro_mean']:.3f} +/- {metrics['f1_macro_std']:.3f}")

        print("\n" + "=" * 60)
    else:
        # Single region summary
        result = all_results[0]
        print("\n" + "=" * 60)
        print(f"EVALUATION RESULTS: {result.get('region', 'Unknown')}")
        print("=" * 60)

        if result.get("status") == "success":
            print(f"\nCells matched: {result['n_cells_matched']}")
            print(f"Cell types: {result['n_cell_types']}")

            print("\nFeature Set Comparison:")
            print("-" * 60)
            for set_name, metrics in result.get("feature_sets", {}).items():
                if "error" not in metrics:
                    print(f"\n{set_name}:")
                    print(f"  Accuracy: {metrics['accuracy_mean']:.3f} +/- {metrics['accuracy_std']:.3f}")
                    print(f"  Macro F1: {metrics['f1_macro_mean']:.3f} +/- {metrics['f1_macro_std']:.3f}")
        else:
            print(f"Status: {result.get('status', 'unknown')}")

        print("\n" + "=" * 60)

    logger.info("Evaluation complete. Results saved to %s", output_dir)


if __name__ == "__main__":
    main()
