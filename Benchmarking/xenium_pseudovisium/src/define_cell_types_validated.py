"""
Spatially-validated cell type definitions based on single-cell co-expression analysis.

These definitions are derived from bivariate Moran's I and single-cell co-expression
validation on 465,534 Xenium cells (2026-01-22 analysis).

Key changes from artificial GT:
1. Epithelial: PanCK only (E-Cadherin removed - r=0.065 with PanCK, doesn't co-express)
2. EMT cells: Vimentin + E-Cadherin (NEW - r=0.370, OR=3.7, spatially coherent 2.72x)
3. Fibroblasts: alphaSMA only (Vimentin goes to EMT or kept separate)
4. All immune profiles validated at single-cell level

Reference: docs/plans/2026-01-22-revised-gt-definitions-design.md
"""

import logging
from typing import Dict, List

import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.mixture import GaussianMixture

logger = logging.getLogger(__name__)

# ============================================================================
# SPATIALLY-VALIDATED CELL PROFILE DEFINITIONS
# ============================================================================
# These are the "best shot" definitions based on actual marker co-expression
# at the single-cell level. Each profile was validated with:
# - Pearson correlation > 0.3 for multi-marker profiles
# - Odds ratio > 2.0 for co-expression
# - Spatial coherence > 1.5x expected
# ============================================================================

VALIDATED_CELL_PROFILE_DICT = {
    # Immune populations (all validated with OR > 20)
    "B cells": {
        "Major": ["CD20", "CD45RA"],  # r=0.745, OR=488.8
        "Minor": [],
    },
    "T cells": {
        "Major": ["CD3E", "CD45RO"],  # r=0.636, OR=38.5
        "Minor": [],
    },
    "CD8+ T cells": {
        "Major": ["CD3E", "CD8A"],  # Subset of T cells
        "Minor": ["GranzymeB"],
    },
    "Exhausted T cells": {
        "Major": ["CD8A", "PD-1"],  # r=0.631, OR=93.6
        "Minor": [],
    },
    "Macrophages": {
        "Major": ["CD68", "CD16"],  # r=0.744, OR=22.0
        "Minor": [],
    },
    "M2 Macrophages": {
        "Major": ["CD68", "CD163"],  # Subset with M2 polarization
        "Minor": [],
    },
    # Structural populations
    "Endothelial": {
        "Major": ["CD31"],
        "Minor": [],
    },
    "Epithelial": {
        "Major": ["PanCK"],  # E-Cadherin REMOVED (r=0.065 with PanCK)
        "Minor": [],
    },
    "EMT cells": {
        "Major": ["Vimentin", "E-Cadherin"],  # r=0.370, OR=3.7, coherence=2.72x
        "Minor": [],
    },
    "Fibroblasts": {
        "Major": ["alphaSMA"],  # Vimentin removed (goes to EMT)
        "Minor": [],
    },
}

# Priority order (higher = assigned first when tied)
VALIDATED_CELL_TYPE_PRIORITY = {
    "B cells": 10,
    "CD8+ T cells": 10,
    "Exhausted T cells": 9,
    "T cells": 8,
    "Macrophages": 7,
    "M2 Macrophages": 6,
    "Endothelial": 5,
    "EMT cells": 4,  # NEW - captures transitional cells
    "Epithelial": 3,
    "Fibroblasts": 2,
}

# Simplified 8-type version for fair comparison with other methods
# (Merges T cell subtypes, macrophage subtypes)
VALIDATED_8TYPE_PROFILE_DICT = {
    "B cells": {
        "Major": ["CD20", "CD45RA"],
        "Minor": [],
    },
    "T cells": {
        "Major": ["CD3E"],
        "Minor": ["CD45RO", "CD8A", "CD4"],
    },
    "Macrophages": {
        "Major": ["CD68"],
        "Minor": ["CD16", "CD163"],
    },
    "Endothelial": {
        "Major": ["CD31"],
        "Minor": [],
    },
    "Epithelial": {
        "Major": ["PanCK"],
        "Minor": [],
    },
    "EMT cells": {
        "Major": ["Vimentin", "E-Cadherin"],
        "Minor": [],
    },
    "Fibroblasts": {
        "Major": ["alphaSMA"],
        "Minor": [],
    },
    "Stromal": {
        "Major": ["Vimentin"],  # Vimentin+ but NOT E-Cadherin+
        "Minor": [],
        "Negative": ["E-Cadherin", "alphaSMA"],  # Exclude EMT and Fibroblasts
    },
}

VALIDATED_8TYPE_PRIORITY = {
    "B cells": 10,
    "T cells": 9,
    "Macrophages": 8,
    "Endothelial": 7,
    "EMT cells": 6,  # Check EMT before pure stromal
    "Epithelial": 5,
    "Fibroblasts": 4,
    "Stromal": 3,
}


def get_protein_thresholds_gmm(
    adata_protein: sc.AnnData,
    n_components: int = 2,
) -> Dict[str, float]:
    """Determine positive/negative thresholds using GMM."""
    thresholds = {}

    for protein in adata_protein.var_names:
        expr = np.asarray(adata_protein[:, protein].X.todense()).flatten()
        expr = expr[expr > 0]

        if len(expr) < 100:
            thresholds[protein] = np.median(expr) if len(expr) > 0 else 0
            continue

        try:
            gmm = GaussianMixture(n_components=n_components, random_state=42)
            gmm.fit(expr.reshape(-1, 1))
            means = gmm.means_.flatten()
            thresholds[protein] = np.mean(means)
        except Exception:
            thresholds[protein] = np.median(expr)

    return thresholds


def classify_cells_validated(
    adata_protein: sc.AnnData,
    cell_profile_dict: Dict[str, Dict[str, List[str]]] = None,
    priority_dict: Dict[str, int] = None,
    threshold_method: str = "gmm",
) -> pd.Series:
    """
    Classify cells using spatially-validated profiles.

    Supports negative markers (cell must NOT express these).

    Args:
        adata_protein: AnnData with protein expression
        cell_profile_dict: Profile definitions (default: VALIDATED_8TYPE_PROFILE_DICT)
        priority_dict: Priority for tie-breaking (default: VALIDATED_8TYPE_PRIORITY)
        threshold_method: "gmm" or "percentile"

    Returns:
        Series mapping cell_id to cell_type
    """
    if cell_profile_dict is None:
        cell_profile_dict = VALIDATED_8TYPE_PROFILE_DICT
    if priority_dict is None:
        priority_dict = VALIDATED_8TYPE_PRIORITY

    # Get thresholds
    thresholds = get_protein_thresholds_gmm(adata_protein)

    # Get expression matrix
    if hasattr(adata_protein.X, "toarray"):
        X = adata_protein.X.toarray()
    else:
        X = adata_protein.X

    protein_names = list(adata_protein.var_names)
    n_cells = len(adata_protein)

    # Create binary positive matrix
    is_positive = {}
    for protein in protein_names:
        idx = protein_names.index(protein)
        is_positive[protein] = X[:, idx] > thresholds.get(protein, 0)

    # Score each cell for each type
    cell_scores = pd.DataFrame(0.0, index=adata_protein.obs_names, columns=list(cell_profile_dict.keys()))

    for cell_type, markers in cell_profile_dict.items():
        major_markers = markers.get("Major", [])
        minor_markers = markers.get("Minor", [])
        negative_markers = markers.get("Negative", [])

        available_major = [m for m in major_markers if m in is_positive]
        available_minor = [m for m in minor_markers if m in is_positive]
        available_negative = [m for m in negative_markers if m in is_positive]

        if not available_major:
            continue

        # All major markers must be positive
        major_score = np.ones(n_cells, dtype=bool)
        for marker in available_major:
            major_score = major_score & is_positive[marker]

        # Check negative markers (must NOT be positive)
        for marker in available_negative:
            major_score = major_score & ~is_positive[marker]

        major_score = major_score.astype(float)

        # Minor marker bonus
        minor_score = 0
        if available_minor:
            for marker in available_minor:
                minor_score += is_positive[marker].astype(float) * 0.1
            minor_score /= len(available_minor)

        cell_scores[cell_type] = major_score + minor_score

    # Add priority for tie-breaking
    for cell_type in cell_scores.columns:
        priority = priority_dict.get(cell_type, 0)
        cell_scores[cell_type] += cell_scores[cell_type] * (priority * 0.001)

    # Assign best type
    cell_types = cell_scores.idxmax(axis=1)
    max_scores = cell_scores.max(axis=1)
    cell_types[max_scores == 0] = "Unassigned"

    # Summary
    type_counts = cell_types.value_counts()
    logger.info(f"Validated cell type distribution:\n{type_counts}")

    return cell_types


# Mapping from validated 8-type to common evaluation types
VALIDATED_TO_COMMON_MAPPING = {
    "B cells": "B cells",
    "T cells": "T cells",
    "Macrophages": "Macrophages",
    "Endothelial": "Endothelial",
    "Epithelial": "Epithelial",
    "EMT cells": "EMT cells",
    "Fibroblasts": "Fibroblasts",
    "Stromal": "Stromal",
}
