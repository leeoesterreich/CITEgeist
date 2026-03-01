#!/usr/bin/env python
"""
Two-Stage Morphology Benchmark: Hybrid proportions + morphology-based single-cell assignment.

Stage 1: Hybrid CITEgeist (protein → proportions, r=0.74)
Stage 2: Morphology + Hungarian assignment (constrained by Stage 1 counts)

Evaluation: Compare per-cell assignments to original Xenium GT cell types.
"""
import argparse
import json
import logging
import sys
import time
from pathlib import Path
from typing import Dict, List, Tuple

import cv2
import numpy as np
import pandas as pd

# Setup paths
REPO_ROOT = Path(__file__).parent.parent.parent.parent.parent
sys.path.insert(0, str(REPO_ROOT))

from CITEgeist.model.stage2_pipeline import Stage2Pipeline
from CITEgeist.model.morphology_features import extract_patch
from CITEgeist.model.constrained_assignment import proportions_to_counts, random_assign

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

# Data paths
XENIUM_DIR = Path("/ix1/alee/LO_LAB/General/Public_Data/10x_PublicData/Xenium_RNA_Proteomic_RenalCellCarcinoma")
PSEUDOVISIUM_DIR = REPO_ROOT / "Benchmarking/xenium_pseudovisium/data_protein_gt"
IMAGE_DIR = REPO_ROOT / "Benchmarking/xenium_benchmarking/scResolve/images/morphology_hires"
HYBRID_OUTPUT_DIR = REPO_ROOT / "Benchmarking/xenium_benchmarking/CITEgeist/output/hybrid_detection_filter"

# Cell types (achievable-7)
CELL_TYPES = ["B cells", "CD4+ T cells", "CD8+ T cells", "Macrophages", "Endothelial", "Epithelial", "Fibroblasts"]
CELL_TYPE_TO_IDX = {ct: i for i, ct in enumerate(CELL_TYPES)}


def load_xenium_cells() -> pd.DataFrame:
    """Load Xenium cell coordinates and GT types."""
    # Load cell coordinates
    cells_df = pd.read_parquet(XENIUM_DIR / "cells.parquet")
    cells_df = cells_df.set_index("cell_id")

    # Load GT cell types
    gt_df = pd.read_csv(PSEUDOVISIUM_DIR / "cell_type_assignments.csv", index_col=0)

    # Merge
    cells_df = cells_df.join(gt_df, how="inner")

    # Load cell-to-spot mapping
    mapping_df = pd.read_csv(PSEUDOVISIUM_DIR / "cell_to_spot_mapping.csv", index_col=0)
    cells_df = cells_df.join(mapping_df[["spot_id", "region_id"]], how="inner")

    # Filter to cells assigned to spots
    cells_df = cells_df[cells_df["spot_id"].notna() & (cells_df["spot_id"] != "")]

    # Filter to known cell types
    cells_df = cells_df[cells_df["cell_type"].isin(CELL_TYPES)]

    logger.info(f"Loaded {len(cells_df)} Xenium cells with GT and spot assignments")
    return cells_df


def load_morphology_image(region_id: int) -> Tuple[np.ndarray, Dict]:
    """Load morphology image and coordinate info for a region."""
    region_name = f"Xenium_region_{region_id}"

    # Load image
    image_path = IMAGE_DIR / region_name / "morphology.png"
    bgr = cv2.imread(str(image_path), cv2.IMREAD_COLOR)
    rgb = cv2.cvtColor(bgr, cv2.COLOR_BGR2RGB)

    # Load coord info
    with open(IMAGE_DIR / region_name / "coord_info.json") as f:
        coord_info = json.load(f)

    return rgb, coord_info


def load_hybrid_proportions(region_id: int) -> pd.DataFrame:
    """Load Stage 1 hybrid proportions."""
    region_name = f"Xenium_region_{region_id}"
    props_path = HYBRID_OUTPUT_DIR / region_name / f"{region_name}_deconv_predictions.csv"

    if not props_path.exists():
        raise FileNotFoundError(f"No hybrid proportions found: {props_path}")

    return pd.read_csv(props_path, index_col=0)


def micron_to_pixel(x_um: float, y_um: float, coord_info: Dict) -> Tuple[float, float]:
    """Convert micron coordinates to pixel coordinates."""
    pixel_size = coord_info["pixel_size"]
    x_min = coord_info["micron_bounds"]["x_min"]
    y_min = coord_info["micron_bounds"]["y_min"]

    x_px = (x_um - x_min) / pixel_size
    y_px = (y_um - y_min) / pixel_size

    return x_px, y_px


def prepare_training_data(
    cells_df: pd.DataFrame,
    region_id: int,
    image: np.ndarray,
    coord_info: Dict,
    purity_threshold: float = 0.7,
    max_per_type: int = 500,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Prepare training data from high-purity spots.

    Args:
        cells_df: DataFrame with cell coords, GT types, spot assignments
        region_id: Region to use for training
        image: Morphology image (H, W, C)
        coord_info: Coordinate transform info
        purity_threshold: Minimum fraction of dominant type in spot
        max_per_type: Maximum samples per cell type

    Returns:
        patches: (N, 2, 64, 64) training patches
        labels: (N,) integer labels
    """
    region_cells = cells_df[cells_df["region_id"] == region_id]

    # Group by spot and find high-purity spots
    spot_groups = region_cells.groupby("spot_id")

    patches_list = []
    labels_list = []
    type_counts = {i: 0 for i in range(len(CELL_TYPES))}

    # Convert image to channel-first for extract_patch
    image_chw = np.transpose(image, (2, 0, 1))[:2]  # Take first 2 channels (DAPI, boundary)

    for spot_id, spot_cells in spot_groups:
        # Check purity
        type_counts_spot = spot_cells["cell_type"].value_counts()
        dominant_type = type_counts_spot.index[0]
        purity = type_counts_spot.iloc[0] / len(spot_cells)

        if purity < purity_threshold:
            continue

        dominant_idx = CELL_TYPE_TO_IDX[dominant_type]

        if type_counts[dominant_idx] >= max_per_type:
            continue

        # Extract patches for cells in this spot
        for _, cell in spot_cells.iterrows():
            if cell["cell_type"] != dominant_type:
                continue  # Only use dominant type cells from high-purity spots

            x_px, y_px = micron_to_pixel(cell["x_centroid"], cell["y_centroid"], coord_info)

            try:
                patch = extract_patch(image_chw, x_px, y_px, size=64)
                patches_list.append(patch)
                labels_list.append(dominant_idx)
                type_counts[dominant_idx] += 1
            except Exception:
                continue

    logger.info(f"Training samples per type: {type_counts}")

    return np.array(patches_list), np.array(labels_list)


def evaluate_region(
    pipeline: Stage2Pipeline,
    cells_df: pd.DataFrame,
    props_df: pd.DataFrame,
    region_id: int,
    image: np.ndarray,
    coord_info: Dict,
) -> Dict:
    """Evaluate Stage 2 on one region."""
    region_cells = cells_df[cells_df["region_id"] == region_id]

    # Convert image
    image_chw = np.transpose(image, (2, 0, 1))[:2]

    results = {
        "n_spots": 0,
        "n_cells": 0,
        "hungarian_correct": 0,
        "random_correct": 0,
        "per_type": {ct: {"correct": 0, "random": 0, "total": 0} for ct in CELL_TYPES},
    }

    # Group by spot
    spot_groups = region_cells.groupby("spot_id")

    for spot_id, spot_cells in spot_groups:
        if spot_id not in props_df.index:
            continue

        n_cells = len(spot_cells)
        if n_cells < 2:
            continue

        # Get Stage 1 proportions
        proportions = props_df.loc[spot_id, CELL_TYPES].values
        counts = proportions_to_counts(proportions, n_cells)

        # Extract patches
        patches = []
        gt_labels = []
        for _, cell in spot_cells.iterrows():
            x_px, y_px = micron_to_pixel(cell["x_centroid"], cell["y_centroid"], coord_info)
            try:
                patch = extract_patch(image_chw, x_px, y_px, size=64)
                patches.append(patch)
                gt_labels.append(CELL_TYPE_TO_IDX.get(cell["cell_type"], -1))
            except Exception:
                continue

        if len(patches) < 2:
            continue

        patches = np.array(patches)
        gt_labels = np.array(gt_labels)

        # Hungarian assignment
        hungarian_pred = pipeline.assign(patches, counts)

        # Random baseline
        random_pred = random_assign(counts, n_samples=len(patches))

        # Evaluate
        results["n_spots"] += 1
        results["n_cells"] += len(patches)
        results["hungarian_correct"] += (hungarian_pred == gt_labels).sum()
        results["random_correct"] += (random_pred == gt_labels).sum()

        # Per-type
        for i, ct in enumerate(CELL_TYPES):
            mask = gt_labels == i
            if mask.sum() > 0:
                results["per_type"][ct]["total"] += mask.sum()
                results["per_type"][ct]["correct"] += (hungarian_pred[mask] == i).sum()
                results["per_type"][ct]["random"] += (random_pred[mask] == i).sum()

    return results


def main():
    parser = argparse.ArgumentParser(description="Two-Stage Morphology Benchmark")
    parser.add_argument("--region", type=int, required=True, help="Region ID (0-4)")
    parser.add_argument("--output-dir", type=str, required=True, help="Output directory")
    parser.add_argument("--train-region", type=int, default=None,
                        help="Region to use for training (default: same as eval)")
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    train_region = args.train_region if args.train_region is not None else args.region
    eval_region = args.region

    logger.info("=" * 70)
    logger.info(f"TWO-STAGE MORPHOLOGY BENCHMARK: Region {eval_region}")
    logger.info(f"Training on region {train_region}")
    logger.info("=" * 70)

    # Load data
    logger.info("Loading Xenium cells...")
    cells_df = load_xenium_cells()

    logger.info(f"Loading morphology image for region {train_region}...")
    train_image, train_coord_info = load_morphology_image(train_region)

    # Prepare training data
    logger.info("Preparing training data from high-purity spots...")
    train_patches, train_labels = prepare_training_data(
        cells_df, train_region, train_image, train_coord_info,
        purity_threshold=0.7, max_per_type=500,
    )

    if len(train_patches) < 50:
        logger.error(f"Insufficient training data: {len(train_patches)} samples")
        return

    # Train pipeline
    logger.info("Training Stage 2 pipeline...")
    pipeline = Stage2Pipeline(cell_types=CELL_TYPES, n_gmm_components=2)
    pipeline.train(train_patches, train_labels)

    # Evaluate
    logger.info(f"Evaluating on region {eval_region}...")
    if eval_region != train_region:
        eval_image, eval_coord_info = load_morphology_image(eval_region)
    else:
        eval_image, eval_coord_info = train_image, train_coord_info

    props_df = load_hybrid_proportions(eval_region)

    results = evaluate_region(
        pipeline, cells_df, props_df, eval_region, eval_image, eval_coord_info
    )

    # Compute metrics
    hungarian_acc = results["hungarian_correct"] / results["n_cells"] if results["n_cells"] > 0 else 0
    random_acc = results["random_correct"] / results["n_cells"] if results["n_cells"] > 0 else 0

    logger.info("=" * 70)
    logger.info("RESULTS")
    logger.info("=" * 70)
    logger.info(f"Spots evaluated: {results['n_spots']}")
    logger.info(f"Cells evaluated: {results['n_cells']}")
    logger.info(f"Hungarian accuracy: {hungarian_acc:.4f}")
    logger.info(f"Random accuracy:    {random_acc:.4f}")
    logger.info(f"Improvement:        {hungarian_acc - random_acc:+.4f} ({(hungarian_acc/random_acc - 1)*100:+.1f}%)" if random_acc > 0 else "")

    logger.info("\nPer-type accuracy:")
    for ct in CELL_TYPES:
        stats = results["per_type"][ct]
        if stats["total"] > 0:
            h_acc = stats["correct"] / stats["total"]
            r_acc = stats["random"] / stats["total"]
            logger.info(f"  {ct}: Hungarian={h_acc:.3f}, Random={r_acc:.3f} (n={stats['total']})")

    # Save results
    output = {
        "region": eval_region,
        "train_region": train_region,
        "n_train_samples": len(train_patches),
        "results": results,
        "hungarian_accuracy": hungarian_acc,
        "random_accuracy": random_acc,
        "improvement": hungarian_acc - random_acc,
    }

    with open(output_dir / f"region_{eval_region}_results.json", "w") as f:
        json.dump(output, f, indent=2)

    logger.info(f"\nResults saved to {output_dir}")


if __name__ == "__main__":
    main()
