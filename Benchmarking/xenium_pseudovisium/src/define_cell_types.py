"""
Define cell types based on protein marker expression.

This module provides cell type definitions and classification functions
for the Xenium kidney dataset based on the 27 available protein markers.
"""

import logging
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.mixture import GaussianMixture

logger = logging.getLogger(__name__)

# Cell profile dictionary for Xenium kidney dataset
# Organized by cell type with Major markers (required) and Minor markers (supportive)
XENIUM_CELL_PROFILE_DICT = {
    "CD4+ T cells": {
        "Major": ["CD3E", "CD4"],
        "Minor": ["CD45"],
    },
    "CD8+ T cells": {
        "Major": ["CD3E", "CD8A"],
        "Minor": ["CD45", "GranzymeB"],
    },
    "B cells": {
        "Major": ["CD20"],
        "Minor": ["CD45"],
    },
    "Plasma cells": {
        "Major": ["CD138"],
        "Minor": [],
    },
    "Macrophages": {
        "Major": ["CD68"],
        "Minor": ["CD163", "CD45", "HLA-DR"],
    },
    "Dendritic cells": {
        "Major": ["CD11c", "HLA-DR"],
        "Minor": ["CD45"],
    },
    "NK cells": {
        "Major": ["CD16", "GranzymeB"],
        "Minor": ["CD45"],
    },
    "Epithelial": {
        "Major": ["PanCK"],
        "Minor": ["E-Cadherin"],
    },
    "Endothelial": {
        "Major": ["CD31"],
        "Minor": [],
    },
    "Fibroblasts": {
        "Major": ["alphaSMA"],
        "Minor": ["Vimentin"],
    },
}

# Priority order for resolving conflicts (higher priority = assigned first)
CELL_TYPE_PRIORITY = {
    "CD4+ T cells": 10,
    "CD8+ T cells": 10,
    "B cells": 9,
    "Plasma cells": 8,
    "Macrophages": 7,
    "Dendritic cells": 6,
    "NK cells": 5,
    "Epithelial": 4,
    "Endothelial": 3,
    "Fibroblasts": 2,
}


def get_protein_thresholds_gmm(
    adata_protein: sc.AnnData,
    n_components: int = 2,
) -> Dict[str, float]:
    """
    Determine positive/negative thresholds for each protein using GMM.

    Fits a 2-component Gaussian Mixture Model to each protein's expression.
    The threshold is set at the intersection of the two Gaussians.

    Args:
        adata_protein: AnnData with protein expression
        n_components: Number of GMM components (default 2)

    Returns:
        Dict mapping protein name to threshold value
    """
    thresholds = {}

    for protein in adata_protein.var_names:
        # Get protein expression
        expr = np.asarray(adata_protein[:, protein].X.todense()).flatten()
        expr = expr[expr > 0]  # Focus on non-zero values

        if len(expr) < 100:
            # Not enough data, use median
            thresholds[protein] = np.median(expr) if len(expr) > 0 else 0
            continue

        # Fit GMM
        gmm = GaussianMixture(n_components=n_components, random_state=42)
        gmm.fit(expr.reshape(-1, 1))

        # Get means and find threshold between components
        means = gmm.means_.flatten()
        if len(means) == 2:
            # Threshold is midpoint between the two means
            thresholds[protein] = np.mean(means)
        else:
            thresholds[protein] = np.median(expr)

    return thresholds


def get_protein_thresholds_percentile(
    adata_protein: sc.AnnData,
    percentile: float = 75.0,
) -> Dict[str, float]:
    """
    Determine positive/negative thresholds using percentile.

    Args:
        adata_protein: AnnData with protein expression
        percentile: Percentile to use as threshold (default 75)

    Returns:
        Dict mapping protein name to threshold value
    """
    thresholds = {}

    for protein in adata_protein.var_names:
        expr = np.asarray(adata_protein[:, protein].X.todense()).flatten()
        expr_nonzero = expr[expr > 0]

        if len(expr_nonzero) > 0:
            thresholds[protein] = np.percentile(expr_nonzero, percentile)
        else:
            thresholds[protein] = 0

    return thresholds


def classify_cells_by_protein(
    adata_protein: sc.AnnData,
    cell_profile_dict: Dict[str, Dict[str, List[str]]],
    threshold_method: str = "percentile",
    percentile: float = 50.0,
    require_all_major: bool = True,
) -> pd.Series:
    """
    Classify each cell into a cell type based on protein expression.

    Args:
        adata_protein: AnnData with protein expression
        cell_profile_dict: Dict mapping cell types to marker lists
        threshold_method: "gmm" or "percentile"
        percentile: Percentile for thresholding (if method="percentile")
        require_all_major: If True, all Major markers must be positive

    Returns:
        Series with cell_id as index, cell_type as values
        ("Unassigned" for cells not matching any profile)
    """
    # Get thresholds
    if threshold_method == "gmm":
        thresholds = get_protein_thresholds_gmm(adata_protein)
    else:
        thresholds = get_protein_thresholds_percentile(adata_protein, percentile)

    logger.info(f"Using {threshold_method} thresholding (percentile={percentile})")

    # Get expression matrix
    if hasattr(adata_protein.X, "toarray"):
        X = adata_protein.X.toarray()
    else:
        X = adata_protein.X

    # Create binary matrix for each marker
    protein_names = list(adata_protein.var_names)
    is_positive = {}
    for protein in protein_names:
        idx = protein_names.index(protein)
        is_positive[protein] = X[:, idx] > thresholds.get(protein, 0)

    # Score each cell for each cell type
    n_cells = len(adata_protein)
    cell_scores = pd.DataFrame(0.0, index=adata_protein.obs_names, columns=list(cell_profile_dict.keys()))

    for cell_type, markers in cell_profile_dict.items():
        major_markers = markers.get("Major", [])
        minor_markers = markers.get("Minor", [])

        # Check which markers are available
        available_major = [m for m in major_markers if m in is_positive]
        available_minor = [m for m in minor_markers if m in is_positive]

        if not available_major:
            logger.warning(f"No major markers available for {cell_type}: {major_markers}")
            continue

        # Calculate major marker score
        if require_all_major:
            # All major markers must be positive
            major_score = np.ones(n_cells, dtype=bool)
            for marker in available_major:
                major_score = major_score & is_positive[marker]
            major_score = major_score.astype(float)
        else:
            # Average of major markers
            major_score = np.zeros(n_cells)
            for marker in available_major:
                major_score += is_positive[marker].astype(float)
            major_score /= len(available_major)

        # Add minor marker bonus
        minor_score = 0
        if available_minor:
            for marker in available_minor:
                minor_score += is_positive[marker].astype(float) * 0.1
            minor_score /= len(available_minor)

        cell_scores[cell_type] = major_score + minor_score

    # Assign cell types based on highest score
    # Add priority weighting for tie-breaking
    for cell_type in cell_scores.columns:
        priority = CELL_TYPE_PRIORITY.get(cell_type, 0)
        cell_scores[cell_type] += cell_scores[cell_type] * (priority * 0.001)

    # Get best cell type for each cell
    cell_types = cell_scores.idxmax(axis=1)
    max_scores = cell_scores.max(axis=1)

    # Mark cells with score 0 as Unassigned
    cell_types[max_scores == 0] = "Unassigned"

    # Summary
    type_counts = cell_types.value_counts()
    logger.info(f"Cell type distribution:\n{type_counts}")

    return cell_types


def get_cell_type_markers_available(
    adata_protein: sc.AnnData,
    cell_profile_dict: Dict[str, Dict[str, List[str]]],
) -> Dict[str, Dict[str, List[str]]]:
    """
    Check which markers from the profile dict are available in the data.

    Args:
        adata_protein: AnnData with protein expression
        cell_profile_dict: Cell type to marker mapping

    Returns:
        Dict showing available markers for each cell type
    """
    available = {}
    protein_names = set(adata_protein.var_names)

    for cell_type, markers in cell_profile_dict.items():
        available[cell_type] = {
            "Major": [m for m in markers.get("Major", []) if m in protein_names],
            "Minor": [m for m in markers.get("Minor", []) if m in protein_names],
            "Missing_Major": [m for m in markers.get("Major", []) if m not in protein_names],
            "Missing_Minor": [m for m in markers.get("Minor", []) if m not in protein_names],
        }

    return available


def summarize_protein_expression(
    adata_protein: sc.AnnData,
) -> pd.DataFrame:
    """
    Get summary statistics for each protein.

    Args:
        adata_protein: AnnData with protein expression

    Returns:
        DataFrame with expression statistics per protein
    """
    if hasattr(adata_protein.X, "toarray"):
        X = adata_protein.X.toarray()
    else:
        X = adata_protein.X

    summary = []
    for i, protein in enumerate(adata_protein.var_names):
        expr = X[:, i]
        nonzero = expr[expr > 0]
        summary.append({
            "protein": protein,
            "n_positive": len(nonzero),
            "pct_positive": 100 * len(nonzero) / len(expr),
            "mean_all": expr.mean(),
            "mean_positive": nonzero.mean() if len(nonzero) > 0 else 0,
            "median_positive": np.median(nonzero) if len(nonzero) > 0 else 0,
            "max": expr.max(),
        })

    return pd.DataFrame(summary).set_index("protein")


if __name__ == "__main__":
    import sys
    logging.basicConfig(level=logging.INFO)

    sys.path.insert(0, "Benchmarking/xenium_benchmarking/src")
    from load_xenium import load_xenium_data, split_gex_protein

    DATA_DIR = "/ix1/alee/LO_LAB/General/Public_Data/10x_PublicData/Xenium_RNA_Proteomic_RenalCellCarcinoma"

    print("Loading Xenium data...")
    adata = load_xenium_data(DATA_DIR)
    adata_gex, adata_protein = split_gex_protein(adata)

    print("\nProtein expression summary:")
    summary = summarize_protein_expression(adata_protein)
    print(summary.to_string())

    print("\nAvailable markers per cell type:")
    available = get_cell_type_markers_available(adata_protein, XENIUM_CELL_PROFILE_DICT)
    for ct, markers in available.items():
        print(f"{ct}:")
        print(f"  Major: {markers['Major']}")
        if markers['Missing_Major']:
            print(f"  Missing Major: {markers['Missing_Major']}")

    print("\nClassifying cells...")
    cell_types = classify_cells_by_protein(
        adata_protein,
        XENIUM_CELL_PROFILE_DICT,
        threshold_method="percentile",
        percentile=50.0,
    )

    print("\nCell type distribution:")
    print(cell_types.value_counts())
