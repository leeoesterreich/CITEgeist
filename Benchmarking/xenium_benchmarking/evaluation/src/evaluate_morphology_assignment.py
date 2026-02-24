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


# --- Baseline Assignment Methods ---

def run_baseline_random(
    original_assignments: pd.DataFrame,
    seed: int = 42,
) -> pd.DataFrame:
    """
    Random baseline: shuffle cell type assignments within each spot.

    This baseline preserves the cell type distribution per spot but
    randomizes which nucleus gets which cell type. Serves as a lower
    bound for morphology-guided assignment.

    Args:
        original_assignments: DataFrame with nucleus_id, spot_id, cell_type
        seed: Random seed for reproducibility

    Returns:
        DataFrame with same structure, shuffled cell_type within each spot
    """
    rng = np.random.default_rng(seed)
    result = original_assignments.copy()

    for spot_id in result["spot_id"].unique():
        mask = result["spot_id"] == spot_id
        types = result.loc[mask, "cell_type"].values.copy()
        rng.shuffle(types)
        result.loc[mask, "cell_type"] = types

    return result


def run_baseline_uniform(
    spot_props: pd.DataFrame,
    nuclei_per_spot: pd.Series,
    cell_types: List[str],
) -> pd.DataFrame:
    """
    Uniform baseline: all nuclei get equal probability, Hungarian assigns.

    This baseline ignores morphology completely and assigns all nuclei
    uniform probability across cell types. The Hungarian algorithm then
    assigns nuclei to satisfy the target counts from spot proportions.

    Args:
        spot_props: DataFrame indexed by spot_id, columns are cell types with proportions
        nuclei_per_spot: Series mapping spot_id -> nuclei count
        cell_types: List of cell type names (column order for proportions)

    Returns:
        DataFrame with nucleus_id, spot_id, cell_type
    """
    import sys
    sys.path.insert(0, "/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist")
    from CITEgeist.model.morphology_features import largest_remainder_discretize
    from CITEgeist.model.hungarian_assignment import assign_nuclei_to_types

    results = []
    nucleus_counter = 0

    for spot_id in spot_props.index:
        n_nuclei = int(nuclei_per_spot.get(spot_id, 0))
        if n_nuclei == 0:
            continue

        props = spot_props.loc[spot_id, cell_types].values.astype(float)
        counts = largest_remainder_discretize(props, n_nuclei)

        # Uniform probability: all nuclei have equal probability for all types
        probs = np.ones((n_nuclei, len(cell_types))) / len(cell_types)
        nucleus_ids = np.arange(nucleus_counter, nucleus_counter + n_nuclei)

        assignments = assign_nuclei_to_types(probs, counts, nucleus_ids)

        for nid, type_idx in assignments.items():
            results.append({
                "nucleus_id": nid,
                "spot_id": spot_id,
                "cell_type": cell_types[type_idx],
            })

        nucleus_counter += n_nuclei

    return pd.DataFrame(results)


def run_baseline_spot_proportion(
    nuclei_df: pd.DataFrame,
    spot_props: pd.DataFrame,
    nuclei_per_spot: pd.Series,
    cell_types: List[str],
) -> pd.DataFrame:
    """
    Spot-proportion baseline: use spot proportions as probability for all nuclei.

    This baseline assigns all nuclei within a spot the same probability vector
    (the spot's cell type proportions). Unlike morphology-guided assignment,
    it ignores individual nucleus features. The Hungarian algorithm then
    optimally assigns nuclei to cell types.

    Args:
        nuclei_df: DataFrame with nucleus_id, spot_id
        spot_props: DataFrame indexed by spot_id, columns are cell types with proportions
        nuclei_per_spot: Series mapping spot_id -> nuclei count
        cell_types: List of cell type names (column order for proportions)

    Returns:
        DataFrame with nucleus_id, spot_id, cell_type
    """
    import sys
    sys.path.insert(0, "/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist")
    from CITEgeist.model.morphology_features import largest_remainder_discretize
    from CITEgeist.model.hungarian_assignment import assign_nuclei_to_types

    results = []

    for spot_id in spot_props.index:
        spot_nuclei = nuclei_df[nuclei_df["spot_id"] == spot_id]
        n_nuclei = len(spot_nuclei)
        if n_nuclei == 0:
            continue

        props = spot_props.loc[spot_id, cell_types].values.astype(float)
        counts = largest_remainder_discretize(props, n_nuclei)

        # All nuclei in spot get the same probability vector (spot proportions)
        probs = np.tile(props, (n_nuclei, 1))
        nucleus_ids = spot_nuclei["nucleus_id"].values

        assignments = assign_nuclei_to_types(probs, counts, nucleus_ids)

        for nid, type_idx in assignments.items():
            results.append({
                "nucleus_id": int(nid),
                "spot_id": spot_id,
                "cell_type": cell_types[type_idx],
            })

    return pd.DataFrame(results)
