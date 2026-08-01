"""
Shared constants for Xenium benchmarking.

Single source of truth for achievable-7 cell type definitions used by
both run_benchmark.py and evaluate_pipeline_stages.py.
"""

from typing import Dict, List, Set, Tuple

# =============================================================================
# ACHIEVABLE-7 CELL TYPE DEFINITIONS
# =============================================================================
#
# 7 cell types achievable with the 27-antibody Xenium panel.
# Ground truth: protein-gated single-cell classification (hierarchical gating).
#
# Profile rationale:
# - Fibroblasts: αSMA only (VIM excluded due to ubiquity)
# - CD4+ T cells: CD3E+ CD4+ CD8A-
# - CD8+ T cells: CD3E+ CD8A+

ACHIEVABLE_7_CELL_PROFILE_DICT: Dict[str, Dict[str, List[str]]] = {
    "B cells": {
        "Major": ["CD20"],
        "Minor": ["CD45RA"],
    },
    "CD4+ T cells": {
        "Major": ["CD3E", "CD4"],
        "Minor": ["CD45RO"],
    },
    "CD8+ T cells": {
        "Major": ["CD3E", "CD8A"],
        "Minor": ["GranzymeB"],
    },
    "Macrophages": {
        "Major": ["CD68", "CD163"],
        "Minor": ["CD16"],
    },
    "Endothelial": {
        "Major": ["CD31"],
        "Minor": [],
    },
    "Epithelial": {
        "Major": ["PanCK"],
        "Minor": [],
    },
    "Fibroblasts": {
        "Major": ["alphaSMA"],
        "Minor": [],
    },
}

# Set-based signatures for Jaccard matching of autodiscovered profiles
ACHIEVABLE_7_MARKER_SIGNATURES: Dict[str, Set[str]] = {
    "B cells": {"CD20", "CD45RA"},
    "CD4+ T cells": {"CD3E", "CD4", "CD45RO"},
    "CD8+ T cells": {"CD3E", "CD8A", "GranzymeB"},
    "Macrophages": {"CD68", "CD163", "CD16"},
    "Endothelial": {"CD31"},
    "Epithelial": {"PanCK"},
    "Fibroblasts": {"alphaSMA"},
}

# Primary/secondary format for profile scoring in staged evaluation
# Scoring formula: (2 * primary_overlap + secondary_overlap) / 3
ACHIEVABLE_7_GT_MARKERS: Dict[str, Dict[str, List[str]]] = {
    "B cells": {
        "primary": ["CD20"],
        "secondary": ["CD45RA"],
    },
    "CD4+ T cells": {
        "primary": ["CD4"],
        "secondary": ["CD3E", "CD45RO"],
    },
    "CD8+ T cells": {
        "primary": ["CD8A"],
        "secondary": ["CD3E", "GranzymeB"],
    },
    "Macrophages": {
        "primary": ["CD68", "CD163"],
        "secondary": ["CD16"],
    },
    "Endothelial": {
        "primary": ["CD31"],
        "secondary": [],
    },
    "Epithelial": {
        "primary": ["PanCK"],
        "secondary": [],
    },
    "Fibroblasts": {
        "primary": ["alphaSMA"],
        "secondary": [],
    },
}

# Protein GT → Achievable-7 mapping (identity — protein GT already has correct types)
GT_TO_ACHIEVABLE_7_MAPPING: Dict[str, str] = {
    "B cells": "B cells",
    "CD4+ T cells": "CD4+ T cells",
    "CD8+ T cells": "CD8+ T cells",
    "Macrophages": "Macrophages",
    "Endothelial": "Endothelial",
    "Epithelial": "Epithelial",
    "Fibroblasts": "Fibroblasts",
}

# GT → Achievable-6 mapping (merges CD4+/CD8+ T cells into "T cells")
GT_TO_ACHIEVABLE_6_MAPPING: Dict[str, str] = {
    "B cells": "B cells",
    "CD4+ T cells": "T cells",
    "CD8+ T cells": "T cells",
    "Macrophages": "Macrophages",
    "Endothelial": "Endothelial",
    "Epithelial": "Epithelial",
    "Fibroblasts": "Fibroblasts",
}

# Critical markers that MUST be flagged as interesting in Module 1
CRITICAL_MARKERS: List[str] = [
    "CD3E",
    "CD4",
    "CD8A",  # T cells
    "CD68",
    "CD163",  # Macrophages
    "CD20",  # B cells
    "PanCK",  # Epithelial
    "CD31",  # Endothelial
    "alphaSMA",  # Fibroblasts
]

# Expected colocalization pairs for Module 2a validation
EXPECTED_POSITIVE_PAIRS = [
    ("CD3E", "CD8A"),  # T cell markers
    ("CD68", "CD163"),  # Macrophage markers
    ("CD68", "HLA-DR"),  # Macrophage markers
    ("CD20", "CD45RA"),  # B cell markers
]

EXPECTED_NEGATIVE_PAIRS = [
    ("CD3E", "CD68"),  # T cells vs Macrophages
    ("CD20", "CD68"),  # B cells vs Macrophages
    ("PanCK", "CD20"),  # Epithelial vs B cells
]

# =============================================================================
# MARKER LINEAGE MAPPING (for discovery comparison evaluation)
# =============================================================================
#
# Maps each of the 27 Xenium protein markers to its canonical cell lineage.
# Used to evaluate whether discovered co-expression modules are biologically
# coherent (single-lineage) or mixed (cross-lineage).
#
# "Functional" markers (checkpoints, proliferation, broadly expressed) are
# excluded from coherence scoring since they can legitimately appear in
# any lineage.

MARKER_LINEAGE_MAP: Dict[str, str] = {
    # T cell lineage
    "CD3E": "T cell",
    "CD4": "T cell",
    "CD8A": "T cell",
    "CD45RO": "T cell",
    "GranzymeB": "T cell",
    # B cell lineage
    "CD20": "B cell",
    "CD45RA": "B cell",
    # Myeloid lineage
    "CD68": "Myeloid",
    "CD163": "Myeloid",
    "CD16": "Myeloid",
    "CD11c": "Myeloid",
    "HLA-DR": "Myeloid",
    # Stromal / Mesenchymal
    "alphaSMA": "Stromal",
    "Vimentin": "Stromal",
    # Epithelial
    "PanCK": "Epithelial",
    "E-Cadherin": "Epithelial",
    "Beta-catenin": "Epithelial",
    # Endothelial
    "CD31": "Endothelial",
    # Plasma cell
    "CD138": "Plasma cell",
}

# =============================================================================
# ACHIEVABLE-6 CELL TYPE DEFINITIONS (T cells merged)
# =============================================================================
#
# 6 cell types for fair comparison with RNA-based methods that cannot
# reliably distinguish CD4+ and CD8+ T cells.

ACHIEVABLE_6_CELL_PROFILE_DICT: Dict[str, Dict[str, List[str]]] = {
    "B cells": {
        "Major": ["CD20"],
        "Minor": ["CD45RA"],
    },
    "T cells": {
        "Major": ["CD3E"],
        "Minor": ["CD45RO"],
    },
    "Macrophages": {
        "Major": ["CD68", "CD163"],
        "Minor": ["CD16"],
    },
    "Endothelial": {
        "Major": ["CD31"],
        "Minor": [],
    },
    "Epithelial": {
        "Major": ["PanCK"],
        "Minor": [],
    },
    "Fibroblasts": {
        "Major": ["alphaSMA"],
        "Minor": [],
    },
}

# Markers excluded from coherence scoring (legitimately cross-lineage)
FUNCTIONAL_MARKERS: Set[str] = {
    "CD45",  # Pan-immune
    "PD-1",  # Checkpoint
    "PD-L1",  # Checkpoint
    "LAG-3",  # Checkpoint
    "VISTA",  # Checkpoint
    "Ki-67",  # Proliferation
    "PCNA",  # Proliferation
    "PTEN",  # Broadly expressed
}

# =============================================================================
# NB MARKER-TYPE TABLE (6-type, for graded marker ablation)
# =============================================================================
#
# Source: CITEgeist/model/_archive/nb_initialization.py MARKER_TYPE_TABLE
# Adapted from 7-type (CD4+T/CD8+T separate) to 6-type (merged T cells).
# "strong" -> 1.0 assignment, "soft" -> graded weight (default 0.1).

NB_MARKER_TYPE_TABLE_6TYPE: Dict[str, list] = {
    # --- Already Major in ACHIEVABLE_6 ---
    "CD20": [("B cells", "strong")],
    "CD3E": [("T cells", "strong")],
    "CD4": [("T cells", "strong")],
    "CD8A": [("T cells", "strong")],
    "CD68": [("Macrophages", "strong")],
    "CD163": [("Macrophages", "strong")],
    "CD31": [("Endothelial", "strong")],
    "PanCK": [("Epithelial", "strong")],
    "alphaSMA": [("Fibroblasts", "strong")],
    # --- New markers to test ---
    "CD16": [("Macrophages", "strong")],
    "HLA-DR": [("B cells", "soft"), ("Macrophages", "strong")],
    "CD11c": [("Macrophages", "soft")],
    "CD45": [("B cells", "strong"), ("T cells", "strong"), ("Macrophages", "strong")],
    "CD138": [("B cells", "soft")],
    "E-Cadherin": [("Macrophages", "soft")],
    "Vimentin": [("Macrophages", "soft"), ("Endothelial", "soft"), ("Fibroblasts", "strong")],
    "CD45RA": [("B cells", "strong")],
    "CD45RO": [("T cells", "strong")],
}

# The 9 markers NOT currently Major in ACHIEVABLE_6
ABLATION_MARKERS: List[str] = [
    "Vimentin",
    "HLA-DR",
    "CD11c",
    "CD45",
    "CD138",
    "E-Cadherin",
    "CD45RA",
    "CD45RO",
    "CD16",
]


def build_ablation_profile(
    base_profile: Dict[str, Dict[str, List[str]]],
    marker: str,
    config: str,
    nb_table: Dict[str, list] = None,
    soft_assign_weight: float = 0.1,
) -> Tuple[Dict[str, Dict[str, List]], float]:
    """Build a modified profile dict that adds one marker with the specified config.

    Args:
        base_profile: ACHIEVABLE_6_CELL_PROFILE_DICT (Major-only, copied internally).
        marker: Marker name to add (must be in nb_table).
        config: One of "binary", "soft_assign", "soft_both".
        nb_table: Marker-type table. Defaults to NB_MARKER_TYPE_TABLE_6TYPE.
        soft_assign_weight: Assignment weight for soft config (default 0.1).

    Returns:
        Tuple of (modified_profile_dict, marker_weight_scalar).
        marker_weight_scalar is 1.0 for binary/soft_assign, 0.3 for soft_both.
        Benchmark converts to dict {marker: scalar} for soft_both config.
    """
    import copy

    if nb_table is None:
        nb_table = NB_MARKER_TYPE_TABLE_6TYPE

    profile = copy.deepcopy(base_profile)
    associations = nb_table[marker]

    # Determine marker_weight for this marker
    mw_value = 0.3 if config == "soft_both" else 1.0

    for cell_type, strength in associations:
        if cell_type not in profile:
            continue

        if "Major" not in profile[cell_type]:
            profile[cell_type]["Major"] = []
        if "Soft" not in profile[cell_type]:
            profile[cell_type]["Soft"] = []

        if config == "binary":
            # Add as Major (binary 1.0)
            if marker not in profile[cell_type]["Major"]:
                profile[cell_type]["Major"].append(marker)
        else:
            # soft_assign or soft_both
            if strength == "strong":
                # Strong types get 1.0 (Major), not soft
                if marker not in profile[cell_type]["Major"]:
                    profile[cell_type]["Major"].append(marker)
            else:
                # Soft types get graded assignment
                profile[cell_type]["Soft"].append((marker, soft_assign_weight))

    return profile, mw_value
