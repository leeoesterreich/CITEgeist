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

import matplotlib
matplotlib.use("Agg")  # Non-interactive backend for SLURM
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.spatial import cKDTree
from sklearn.metrics import precision_recall_fscore_support, confusion_matrix as sklearn_confusion_matrix

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


# --- Accuracy Metrics ---

def compute_accuracy_metrics(
    pred_labels: pd.Series,
    gt_labels: pd.Series,
    cell_types: List[str],
) -> Dict[str, float]:
    """
    Compute accuracy, precision, recall, F1 per cell type and overall.

    Args:
        pred_labels: Predicted cell type labels
        gt_labels: Ground truth cell type labels
        cell_types: List of cell type names

    Returns:
        Dict with all metrics including:
        - overall_accuracy: Fraction of correct predictions
        - n_cells: Number of cells evaluated
        - {celltype}_precision, {celltype}_recall, {celltype}_f1, {celltype}_support
        - macro_precision, macro_recall, macro_f1
    """
    common_idx = pred_labels.index.intersection(gt_labels.index)
    pred = pred_labels.loc[common_idx]
    gt = gt_labels.loc[common_idx]

    metrics = {}
    metrics["overall_accuracy"] = (pred == gt).mean()
    metrics["n_cells"] = len(common_idx)

    precision, recall, f1, support = precision_recall_fscore_support(
        gt, pred, labels=cell_types, average=None, zero_division=0
    )

    for i, ct in enumerate(cell_types):
        metrics[f"{ct}_precision"] = precision[i]
        metrics[f"{ct}_recall"] = recall[i]
        metrics[f"{ct}_f1"] = f1[i]
        metrics[f"{ct}_support"] = int(support[i])

    metrics["macro_precision"] = np.mean(precision)
    metrics["macro_recall"] = np.mean(recall)
    metrics["macro_f1"] = np.mean(f1)

    return metrics


def compute_confusion_matrix(
    pred_labels: pd.Series,
    gt_labels: pd.Series,
    cell_types: List[str],
) -> np.ndarray:
    """
    Compute confusion matrix (rows=actual/ground truth, cols=predicted).

    Args:
        pred_labels: Predicted cell type labels
        gt_labels: Ground truth cell type labels
        cell_types: List of cell type names (defines row/column order)

    Returns:
        Confusion matrix as numpy array, shape (n_types, n_types)
        Entry [i,j] = number of cells with true label i predicted as label j
    """
    common_idx = pred_labels.index.intersection(gt_labels.index)
    pred = pred_labels.loc[common_idx]
    gt = gt_labels.loc[common_idx]

    cm = sklearn_confusion_matrix(gt, pred, labels=cell_types)
    return cm


# --- Visualization Functions ---

def plot_confusion_matrix(
    cm: np.ndarray,
    cell_types: List[str],
    output_path: Path,
    title: str = "Confusion Matrix",
    normalize: bool = True,
) -> None:
    """
    Plot confusion matrix heatmap.

    Args:
        cm: Confusion matrix array (rows=actual, cols=predicted)
        cell_types: List of cell type names (row/column labels)
        output_path: Path to save the plot
        title: Plot title
        normalize: If True, normalize rows to sum to 1
    """
    if normalize:
        row_sums = cm.sum(axis=1, keepdims=True)
        row_sums[row_sums == 0] = 1  # Avoid division by zero
        cm_plot = cm.astype(float) / row_sums
        fmt = ".2f"
    else:
        cm_plot = cm
        fmt = "d"

    fig, ax = plt.subplots(figsize=(10, 8))
    sns.heatmap(
        cm_plot,
        annot=True,
        fmt=fmt,
        cmap="Blues",
        xticklabels=cell_types,
        yticklabels=cell_types,
        ax=ax,
    )
    ax.set_xlabel("Predicted")
    ax.set_ylabel("Actual")
    ax.set_title(title)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()

    logger.info(f"Saved confusion matrix: {output_path}")


def plot_baseline_comparison(
    results: Dict[str, Dict[str, float]],
    output_path: Path,
) -> None:
    """
    Plot bar chart comparing methods on accuracy and F1.

    Args:
        results: Dict mapping method name -> metrics dict with
                 'overall_accuracy' and 'macro_f1' keys
        output_path: Path to save the plot
    """
    methods = list(results.keys())
    accuracies = [results[m].get("overall_accuracy", 0) for m in methods]
    f1_scores = [results[m].get("macro_f1", 0) for m in methods]

    x = np.arange(len(methods))
    width = 0.35

    fig, ax = plt.subplots(figsize=(10, 6))
    bars1 = ax.bar(x - width / 2, accuracies, width, label="Accuracy", color="steelblue")
    bars2 = ax.bar(x + width / 2, f1_scores, width, label="Macro F1", color="darkorange")

    ax.set_xlabel("Method")
    ax.set_ylabel("Score")
    ax.set_title("Morphology Assignment: Baseline Comparison")
    ax.set_xticks(x)
    ax.set_xticklabels(methods, rotation=45, ha="right")
    ax.legend()
    ax.set_ylim(0, 1)

    # Add value labels on bars
    for bar in list(bars1) + list(bars2):
        height = bar.get_height()
        ax.annotate(
            f"{height:.2f}",
            xy=(bar.get_x() + bar.get_width() / 2, height),
            xytext=(0, 3),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=8,
        )

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()

    logger.info(f"Saved baseline comparison: {output_path}")


# --- Main Evaluation Functions ---

import argparse
import json
import scanpy as sc


def load_single_cell_adata(region_id: int, sc_dir: Path) -> Optional[sc.AnnData]:
    """
    Load single-cell AnnData for a region.

    Args:
        region_id: Region number (0, 1, 2, 3, or 4)
        sc_dir: Path to single_cell_resolution output directory

    Returns:
        AnnData object with single-cell assignments, or None if not found
    """
    sc_path = sc_dir / f"Xenium_region_{region_id}" / f"Xenium_region_{region_id}_single_cell.h5ad"
    if not sc_path.exists():
        logger.warning(f"Single-cell AnnData not found: {sc_path}")
        return None

    adata = sc.read_h5ad(sc_path)
    logger.info(f"Loaded single-cell AnnData: {adata.shape[0]} cells, {adata.obs['cell_type'].nunique()} types")
    return adata


def load_xenium_cell_coords(region_id: int, gt_dir: Path) -> pd.DataFrame:
    """
    Load Xenium cell coordinates for a specific region.

    Filters the cell_to_spot_mapping to only include cells in the specified
    region that are mapped to spots, then joins with the full cells.parquet
    to get exact centroid coordinates.

    Args:
        region_id: Region number (0, 1, 2, 3, or 4)
        gt_dir: Path to xenium_pseudovisium directory

    Returns:
        DataFrame with cell_id index and x_centroid, y_centroid columns
    """
    # Load cell-to-spot mapping to filter by region
    mapping_path = gt_dir / "data_protein_gt" / "cell_to_spot_mapping.csv"
    mapping_df = pd.read_csv(mapping_path, index_col=0)

    # Filter to cells in this region that are mapped to spots
    region_cells = mapping_df[
        (mapping_df["region_id"] == region_id) &
        (mapping_df["spot_idx"] != -1)
    ]
    region_cell_ids = set(region_cells.index)
    logger.info(f"Region {region_id}: {len(region_cell_ids)} GT cells mapped to spots")

    # Load full Xenium cell coordinates
    xenium_dir = Path("/ix1/alee/LO_LAB/General/Public_Data/10x_PublicData/Xenium_RNA_Proteomic_RenalCellCarcinoma")
    cells_path = xenium_dir / "cells.parquet"
    cells_df = pd.read_parquet(cells_path)
    cells_df = cells_df.set_index("cell_id")

    # Filter to region cells and select coordinates
    cells_df = cells_df.loc[cells_df.index.isin(region_cell_ids), ["x_centroid", "y_centroid"]]
    logger.info(f"Loaded {len(cells_df)} GT cell coordinates")

    return cells_df


def collapse_t_cells(labels: pd.Series) -> pd.Series:
    """
    Collapse CD4+ and CD8+ T cells to 'T cells' for RNA GT comparison.

    The RNA ground truth uses a merged T cell category, while protein GT
    distinguishes CD4+ and CD8+ T cells. This function harmonizes the
    predicted labels to match the RNA GT format.

    Args:
        labels: Series of cell type labels

    Returns:
        Series with T cell subtypes collapsed to 'T cells'
    """
    return labels.replace({"CD4+ T cells": "T cells", "CD8+ T cells": "T cells"})


def evaluate_region(
    region_id: int,
    sc_dir: Path,
    gt_dir: Path,
    output_dir: Path,
    max_dist: float = 10.0,
) -> Optional[Dict]:
    """
    Evaluate morphology assignment for a single region.

    This is the main evaluation function that:
    1. Loads single-cell AnnData from Module 3b
    2. Loads protein-gated ground truth
    3. Matches Cellpose nuclei to GT cells by spatial proximity
    4. Extracts morphology-guided predictions
    5. Runs all baseline methods
    6. Computes accuracy metrics for each method
    7. Generates confusion matrix visualizations

    Args:
        region_id: Region number (0, 1, 2, 3, or 4)
        sc_dir: Path to single_cell_resolution output directory
        gt_dir: Path to xenium_pseudovisium directory
        output_dir: Path to save results and plots
        max_dist: Maximum distance (microns) for cell matching

    Returns:
        Dict with comprehensive evaluation results, or None if evaluation fails
    """
    logger.info(f"=== Evaluating Region {region_id} ===")

    # Create region output directory
    region_out = output_dir / f"region_{region_id}"
    region_out.mkdir(parents=True, exist_ok=True)

    # Step 1: Load single-cell AnnData
    adata = load_single_cell_adata(region_id, sc_dir)
    if adata is None:
        return None

    # Step 2: Load ground truth
    try:
        gt_df = load_ground_truth("protein", gt_dir)
    except FileNotFoundError as e:
        logger.error(f"Failed to load ground truth: {e}")
        return None

    # Step 3: Load GT cell coordinates for spatial matching
    gt_coords = load_xenium_cell_coords(region_id, gt_dir)
    if len(gt_coords) == 0:
        logger.error(f"No GT cells found for region {region_id}")
        return None

    # Step 4: Prepare Cellpose coordinates for matching
    cellpose_coords = pd.DataFrame({
        "nucleus_id": adata.obs.index.astype(int).values,
        "centroid_x": adata.obs["x"].values,
        "centroid_y": adata.obs["y"].values,
    })

    # Step 5: Match Cellpose nuclei to GT cells
    matches = match_cellpose_to_gt(cellpose_coords, gt_coords, max_dist=max_dist)

    if len(matches) == 0:
        logger.error(f"No nuclei matched to GT cells (max_dist={max_dist})")
        return None

    logger.info(f"Matched {len(matches)}/{len(cellpose_coords)} nuclei")

    # Step 6: Build prediction and GT label series for matched nuclei
    matched_nucleus_ids = list(matches.keys())
    matched_gt_cell_ids = [matches[nid] for nid in matched_nucleus_ids]

    # Get morphology-guided predictions
    morphology_preds = adata.obs.loc[
        [str(nid) for nid in matched_nucleus_ids], "cell_type"
    ]
    morphology_preds.index = matched_nucleus_ids

    # Get ground truth labels
    gt_labels = gt_df.loc[matched_gt_cell_ids, "cell_type"]
    gt_labels.index = matched_nucleus_ids

    # Filter out "Unknown" cells from GT
    valid_mask = gt_labels != "Unknown"
    morphology_preds = morphology_preds[valid_mask]
    gt_labels = gt_labels[valid_mask]

    logger.info(f"Evaluating {len(gt_labels)} cells (excluding Unknown)")

    # Check if we have enough cells
    if len(gt_labels) < 10:
        logger.warning(f"Too few matched cells ({len(gt_labels)}) for evaluation")
        return None

    # Step 7: Prepare data for baselines
    # Get spot-level info from single-cell AnnData
    spot_ids = adata.obs.loc[[str(nid) for nid in matched_nucleus_ids], "spot_id"]
    spot_ids.index = matched_nucleus_ids
    spot_ids = spot_ids[valid_mask]

    nuclei_df = pd.DataFrame({
        "nucleus_id": morphology_preds.index,
        "spot_id": spot_ids.values,
        "cell_type": morphology_preds.values,
    })

    # Build spot proportions from morphology assignments
    spot_type_counts = nuclei_df.groupby(["spot_id", "cell_type"]).size().unstack(fill_value=0)
    spot_totals = spot_type_counts.sum(axis=1)
    spot_props = spot_type_counts.div(spot_totals, axis=0)

    # Ensure all cell types are present
    for ct in PROTEIN_GT_CELL_TYPES:
        if ct not in spot_props.columns:
            spot_props[ct] = 0.0
    spot_props = spot_props[PROTEIN_GT_CELL_TYPES]

    nuclei_per_spot = nuclei_df.groupby("spot_id").size()

    # Step 8: Run baselines
    logger.info("Running baseline methods...")

    # Baseline 1: Random shuffle within spot
    random_df = run_baseline_random(nuclei_df, seed=42)
    random_preds = pd.Series(
        random_df.set_index("nucleus_id")["cell_type"].values,
        index=random_df["nucleus_id"].values
    )
    random_preds = random_preds[random_preds.index.isin(morphology_preds.index)]

    # Baseline 2: Uniform probability (all nuclei equal)
    # Need to rebuild with consistent nucleus IDs
    uniform_df = run_baseline_uniform(spot_props, nuclei_per_spot, PROTEIN_GT_CELL_TYPES)

    # Baseline 3: Spot-proportion only
    spot_prop_df = run_baseline_spot_proportion(
        nuclei_df, spot_props, nuclei_per_spot, PROTEIN_GT_CELL_TYPES
    )
    spot_prop_preds = pd.Series(
        spot_prop_df.set_index("nucleus_id")["cell_type"].values,
        index=spot_prop_df["nucleus_id"].values
    )
    spot_prop_preds = spot_prop_preds[spot_prop_preds.index.isin(morphology_preds.index)]

    # Step 9: Compute metrics for each method
    logger.info("Computing accuracy metrics...")

    results = {
        "region_id": region_id,
        "n_matched_cells": len(gt_labels),
        "n_total_nuclei": len(cellpose_coords),
        "match_rate": len(gt_labels) / len(cellpose_coords),
        "methods": {},
    }

    # Morphology-guided method
    morph_metrics = compute_accuracy_metrics(morphology_preds, gt_labels, PROTEIN_GT_CELL_TYPES)
    results["methods"]["morphology"] = morph_metrics

    # Random baseline
    random_metrics = compute_accuracy_metrics(random_preds, gt_labels, PROTEIN_GT_CELL_TYPES)
    results["methods"]["random"] = random_metrics

    # Spot-proportion baseline
    spot_prop_metrics = compute_accuracy_metrics(spot_prop_preds, gt_labels, PROTEIN_GT_CELL_TYPES)
    results["methods"]["spot_proportion"] = spot_prop_metrics

    # Log summary
    logger.info(f"Morphology accuracy: {morph_metrics['overall_accuracy']:.3f}")
    logger.info(f"Random accuracy: {random_metrics['overall_accuracy']:.3f}")
    logger.info(f"Spot-prop accuracy: {spot_prop_metrics['overall_accuracy']:.3f}")

    # Step 10: Generate confusion matrices
    logger.info("Generating confusion matrices...")

    # Morphology confusion matrix
    cm_morph = compute_confusion_matrix(morphology_preds, gt_labels, PROTEIN_GT_CELL_TYPES)
    plot_confusion_matrix(
        cm_morph, PROTEIN_GT_CELL_TYPES,
        region_out / "confusion_matrix_morphology.png",
        title=f"Region {region_id}: Morphology-Guided Assignment"
    )

    # Random baseline confusion matrix
    cm_random = compute_confusion_matrix(random_preds, gt_labels, PROTEIN_GT_CELL_TYPES)
    plot_confusion_matrix(
        cm_random, PROTEIN_GT_CELL_TYPES,
        region_out / "confusion_matrix_random.png",
        title=f"Region {region_id}: Random Baseline"
    )

    # Spot-proportion confusion matrix
    cm_spot = compute_confusion_matrix(spot_prop_preds, gt_labels, PROTEIN_GT_CELL_TYPES)
    plot_confusion_matrix(
        cm_spot, PROTEIN_GT_CELL_TYPES,
        region_out / "confusion_matrix_spot_proportion.png",
        title=f"Region {region_id}: Spot-Proportion Baseline"
    )

    # Save raw confusion matrices
    np.savez(
        region_out / "confusion_matrices.npz",
        morphology=cm_morph,
        random=cm_random,
        spot_proportion=cm_spot,
        cell_types=PROTEIN_GT_CELL_TYPES,
    )

    logger.info(f"Region {region_id} evaluation complete")
    return results


def aggregate_results(region_results: List[Dict]) -> Dict:
    """
    Aggregate results across all regions.

    Args:
        region_results: List of per-region result dicts

    Returns:
        Dict with aggregated metrics (mean and std across regions)
    """
    methods = ["morphology", "random", "spot_proportion"]
    metrics = ["overall_accuracy", "macro_f1", "macro_precision", "macro_recall"]

    aggregated = {
        "n_regions": len(region_results),
        "total_cells": sum(r["n_matched_cells"] for r in region_results),
        "methods": {},
    }

    for method in methods:
        aggregated["methods"][method] = {}
        for metric in metrics:
            values = [r["methods"][method].get(metric, 0) for r in region_results]
            aggregated["methods"][method][f"{metric}_mean"] = float(np.mean(values))
            aggregated["methods"][method][f"{metric}_std"] = float(np.std(values))

    return aggregated


def main():
    """CLI entry point for morphology assignment evaluation."""
    parser = argparse.ArgumentParser(
        description="Evaluate morphology-guided cell type assignment accuracy",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Evaluate all regions
    python evaluate_morphology_assignment.py \\
        --sc-dir /path/to/single_cell_resolution \\
        --gt-dir /path/to/xenium_pseudovisium \\
        --output-dir /path/to/output

    # Evaluate specific regions
    python evaluate_morphology_assignment.py \\
        --sc-dir /path/to/single_cell_resolution \\
        --gt-dir /path/to/xenium_pseudovisium \\
        --output-dir /path/to/output \\
        --regions "0,1,2"
        """
    )
    parser.add_argument(
        "--sc-dir", type=str, required=True,
        help="Path to single_cell_resolution output directory"
    )
    parser.add_argument(
        "--gt-dir", type=str, required=True,
        help="Path to xenium_pseudovisium directory with ground truth"
    )
    parser.add_argument(
        "--regions", type=str, default="0,1,2,4",
        help="Comma-separated list of region IDs to evaluate (default: 0,1,2,4)"
    )
    parser.add_argument(
        "--output-dir", type=str, required=True,
        help="Directory to save evaluation results and plots"
    )
    parser.add_argument(
        "--max-dist", type=float, default=10.0,
        help="Maximum distance (microns) for Cellpose-to-GT matching (default: 10.0)"
    )

    args = parser.parse_args()

    # Parse paths
    sc_dir = Path(args.sc_dir)
    gt_dir = Path(args.gt_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Parse regions
    regions = [int(r.strip()) for r in args.regions.split(",")]
    logger.info(f"Evaluating regions: {regions}")

    # Evaluate each region
    all_results = []
    for region_id in regions:
        result = evaluate_region(
            region_id, sc_dir, gt_dir, output_dir, max_dist=args.max_dist
        )
        if result is not None:
            all_results.append(result)

    if len(all_results) == 0:
        logger.error("No regions were successfully evaluated")
        return 1

    # Aggregate results
    logger.info("=== Aggregating Results ===")
    aggregated = aggregate_results(all_results)

    # Save per-region results
    for result in all_results:
        region_out = output_dir / f"region_{result['region_id']}"
        with open(region_out / "metrics.json", "w") as f:
            json.dump(result, f, indent=2)

    # Save aggregated results
    with open(output_dir / "aggregated_results.json", "w") as f:
        json.dump(aggregated, f, indent=2)

    # Generate baseline comparison plot
    plot_baseline_comparison(
        aggregated["methods"],
        output_dir / "baseline_comparison.png"
    )

    # Print summary
    logger.info("=" * 50)
    logger.info("EVALUATION SUMMARY")
    logger.info("=" * 50)
    logger.info(f"Regions evaluated: {aggregated['n_regions']}")
    logger.info(f"Total cells: {aggregated['total_cells']}")
    logger.info("")
    for method, metrics in aggregated["methods"].items():
        logger.info(f"{method}:")
        logger.info(f"  Accuracy: {metrics['overall_accuracy_mean']:.3f} +/- {metrics['overall_accuracy_std']:.3f}")
        logger.info(f"  Macro F1: {metrics['macro_f1_mean']:.3f} +/- {metrics['macro_f1_std']:.3f}")

    logger.info(f"\nResults saved to: {output_dir}")
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
