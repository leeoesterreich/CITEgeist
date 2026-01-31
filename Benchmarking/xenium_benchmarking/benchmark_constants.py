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

# Critical markers that MUST be flagged as interesting in Module 1
CRITICAL_MARKERS: List[str] = [
    "CD3E", "CD4", "CD8A",  # T cells
    "CD68", "CD163",  # Macrophages
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

# Markers excluded from coherence scoring (legitimately cross-lineage)
FUNCTIONAL_MARKERS: Set[str] = {
    "CD45",       # Pan-immune
    "PD-1",       # Checkpoint
    "PD-L1",      # Checkpoint
    "LAG-3",      # Checkpoint
    "VISTA",      # Checkpoint
    "Ki-67",      # Proliferation
    "PCNA",       # Proliferation
    "PTEN",       # Broadly expressed
}
