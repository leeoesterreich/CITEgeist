# evaluate_morphology_assignment.py
"""
Evaluate morphology-guided cell type assignment accuracy.

Compares Module 3b's morphology-based nucleus assignment against:
1. Ground truth (Protein-gated and RNA clustering)
2. Baseline methods (random, uniform, spot-proportion-only)
"""

import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Protein-gated ground truth cell types (7 types + Unknown)
PROTEIN_GT_CELL_TYPES = [
    "B cells",
    "CD4+ T cells",
    "CD8+ T cells",
    "Macrophages",
    "Endothelial",
    "Epithelial",
    "Fibroblasts",
]

# RNA clustering ground truth cell types (6 types, merged T cells)
RNA_GT_CELL_TYPES = [
    "B cells",
    "T cells",
    "Macrophages",
    "Endothelial",
    "Epithelial",
    "Fibroblasts",
]


def load_ground_truth(gt_type: str, gt_dir: Path) -> pd.DataFrame:
    """
    Load ground truth cell type assignments.

    Args:
        gt_type: "protein" or "rna"
        gt_dir: Path to xenium_pseudovisium directory

    Returns:
        DataFrame with cell_id index and cell_type column
    """
    gt_dir = Path(gt_dir)

    if gt_type == "protein":
        gt_path = gt_dir / "data_protein_gt" / "cell_type_assignments.csv"
    elif gt_type == "rna":
        gt_path = gt_dir / "data_rna_gt" / "cell_type_assignments.csv"
    else:
        raise ValueError(f"Unknown gt_type: {gt_type}. Use 'protein' or 'rna'.")

    if not gt_path.exists():
        raise FileNotFoundError(f"Ground truth file not found: {gt_path}")

    gt_df = pd.read_csv(gt_path)

    # Handle different index column formats
    # Protein GT has unnamed first column (index), RNA GT has explicit cell_id
    if "cell_id" in gt_df.columns:
        gt_df = gt_df.set_index("cell_id")
    elif gt_df.columns[0] == "Unnamed: 0" or gt_df.columns[0] == "":
        # First column is unnamed index
        gt_df = gt_df.set_index(gt_df.columns[0])
        gt_df.index.name = "cell_id"

    logger.info(f"Loaded {gt_type} GT: {len(gt_df)} cells, {gt_df['cell_type'].nunique()} types")

    return gt_df


def match_cellpose_to_gt(
    cellpose_coords: pd.DataFrame,
    gt_coords: pd.DataFrame,
    max_dist: float = 10.0,
) -> Dict[int, str]:
    """
    Match Cellpose nuclei to ground truth cells by spatial proximity.

    Uses a KD-tree for efficient nearest-neighbor lookup. Each Cellpose
    nucleus is matched to its closest ground truth cell if within max_dist.

    Args:
        cellpose_coords: DataFrame with nucleus_id, centroid_x, centroid_y
        gt_coords: DataFrame with x_centroid, y_centroid (index = cell_id)
        max_dist: Maximum distance in microns for a valid match

    Returns:
        Dict mapping nucleus_id -> gt_cell_id for matched nuclei
    """
    gt_xy = gt_coords[["x_centroid", "y_centroid"]].values
    gt_tree = cKDTree(gt_xy)
    gt_ids = gt_coords.index.tolist()

    cellpose_xy = cellpose_coords[["centroid_x", "centroid_y"]].values
    nucleus_ids = cellpose_coords["nucleus_id"].values

    distances, indices = gt_tree.query(cellpose_xy, k=1)

    matches = {}
    for i, (nid, dist, gt_idx) in enumerate(zip(nucleus_ids, distances, indices)):
        if dist <= max_dist:
            matches[int(nid)] = gt_ids[gt_idx]

    logger.info(f"Matched {len(matches)}/{len(nucleus_ids)} nuclei to GT (max_dist={max_dist}µm)")

    return matches
