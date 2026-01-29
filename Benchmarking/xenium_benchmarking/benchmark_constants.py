"""
Shared constants for Xenium benchmarking.

Single source of truth for achievable-7 cell type definitions used by
both run_benchmark.py and evaluate_pipeline_stages.py.
"""

from typing import Dict, List, Set

# =============================================================================
# ACHIEVABLE-7 CELL TYPE DEFINITIONS
# =============================================================================
#
# 7 cell types achievable with the 27-antibody Xenium panel.
# Collapsed from 10 granular RNA-based ground truth types.
#
# Collapse rationale:
# - Myofibroblasts + Stromal → Fibroblasts (both VIM+, αSMA overlaps)
# - Vascular Stromal → Endothelial (CD31+)
# - Proliferating T → CD8+ T cells (both CD3E+)
# - Mixed Immune → CD4+ T cells (HLA-DR+ T cells)

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
        "Major": ["alphaSMA", "Vimentin"],
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
    "Fibroblasts": {"alphaSMA", "Vimentin"},
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
        "primary": ["alphaSMA", "Vimentin"],
        "secondary": [],
    },
}

# 10 → 7 collapse mapping for ground truth evaluation
GT_TO_ACHIEVABLE_7_MAPPING: Dict[str, str] = {
    "B cells": "B cells",
    "Mixed Immune": "CD4+ T cells",
    "CD8+ T cells": "CD8+ T cells",
    "Proliferating T": "CD8+ T cells",
    "Macrophages": "Macrophages",
    "Endothelial": "Endothelial",
    "Vascular Stromal": "Endothelial",
    "Epithelial": "Epithelial",
    "Myofibroblasts": "Fibroblasts",
    "Stromal": "Fibroblasts",
}

# Critical markers that MUST be flagged as interesting in Module 1
CRITICAL_MARKERS: List[str] = [
    "CD3E", "CD4", "CD8A",  # T cells
    "CD68", "CD163",  # Macrophages
    "CD20",  # B cells
    "PanCK",  # Epithelial
    "CD31",  # Endothelial
    "alphaSMA", "Vimentin",  # Fibroblasts
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
