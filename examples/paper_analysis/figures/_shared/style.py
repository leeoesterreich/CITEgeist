#!/usr/bin/env python3
"""
Shared style configuration for CITEgeist manuscript figures.

This module provides consistent styling across all figures for Cell/Cell Reports
publication quality. All figure generation scripts should import and use this.

Usage:
    from _shared.style import apply_style, PALETTE, CELL_TYPE_COLORS, METHOD_COLORS
    apply_style()
"""

import os
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import font_manager

# Journal font profiles (pt).  These are the values the manuscript figures were
# generated with; "nature" is the default used throughout.
JOURNAL_FONTS = {
    "nature": {
        "panel_label_pt": 5,
        "title_pt": 8,
        "axis_label_pt": 7,
        "tick_label_pt": 6,
        "annotation_pt": 6.5,
        "data_label_pt": 6,
        "panel_label_weight": "bold",
        "max_cross_panel_delta_pt": 1.0,
    },
    "science": {
        "panel_label_pt": 5,
        "title_pt": 7,
        "axis_label_pt": 6.5,
        "tick_label_pt": 5.5,
        "annotation_pt": 6,
        "data_label_pt": 5.5,
        "panel_label_weight": "bold",
        "max_cross_panel_delta_pt": 1.0,
    },
    "cell": {
        "panel_label_pt": 6,
        "title_pt": 8,
        "axis_label_pt": 7,
        "tick_label_pt": 6,
        "annotation_pt": 6.5,
        "data_label_pt": 6,
        "panel_label_weight": "bold",
        "max_cross_panel_delta_pt": 1.0,
    },
}
DEFAULT_JOURNAL = "nature"


def resolve_font_profile(spec: dict) -> dict:
    """Resolve a font profile: journal defaults merged with explicit overrides."""
    base = dict(JOURNAL_FONTS.get(spec.get("journal", DEFAULT_JOURNAL), JOURNAL_FONTS[DEFAULT_JOURNAL]))
    base.update(spec.get("fonts", {}))
    return base


# The manuscript figures use Arial.  Point CITEGEIST_FIGURE_FONT at a local
# Arial.ttf to reproduce them exactly; otherwise matplotlib's sans-serif is used.
_ARIAL_TTF = Path(os.environ.get("CITEGEIST_FIGURE_FONT", "Arial.ttf"))
if _ARIAL_TTF.exists():
    font_manager.fontManager.addfont(str(_ARIAL_TTF))
    _FONT_FAMILY = "Arial"
else:
    _FONT_FAMILY = "sans-serif"

# Primary color palette for general use
PALETTE = {
    "primary": "#2171b5",  # Deep blue
    "secondary": "#6baed6",  # Light blue
    "accent1": "#fe9929",  # Orange
    "accent2": "#41ab5d",  # Green
    "neutral": "#636363",  # Gray
    "highlight": "#c51b8a",  # Magenta
    "background": "#f8f9fa",  # Light background
    "border": "#dee2e6",  # Border gray
}

# Consistent cell type colors across all figures
CELL_TYPE_COLORS = {
    # Major immune cells
    "T cells": "#c51b8a",
    "CD4+ T cells": "#c51b8a",
    "CD4_T_Cells": "#c51b8a",
    "CD8+ T cells": "#ae017e",
    "CD8_T_Cells": "#ae017e",
    "B cells": "#6baed6",
    "B_Cells": "#6baed6",
    "Macrophages": "#41ab5d",
    "Monocytes": "#78c679",
    # Stromal cells
    "Fibroblasts": "#fe9929",
    "Endothelial": "#636363",
    # Epithelial/Cancer
    "Epithelial": "#2171b5",
    "Cancer_Luminal": "#2171b5",
    "Cancer_Basal": "#08519c",
    # Other
    "Dendritic_Cells": "#fdd0a2",
    "Unknown": "#bdbdbd",
    "Unassigned": "#bdbdbd",
}

# Method colors for benchmarking comparisons
METHOD_COLORS = {
    "CITEgeist": "#c51b8a",  # Highlight color (our method)
    "Cell2Location": "#2171b5",  # Primary blue
    "RCTD": "#41ab5d",  # Green
    "Tangram": "#fe9929",  # Orange
    "Seurat": "#636363",  # Gray
    "CARD": "#fdd0a2",  # Light orange
}

# Display labels for methods in figures
METHOD_DISPLAY = {
    "CITEgeist": "CITEgeist",
    "Cell2Location": "Cell2Location",
    "RCTD": "RCTD",
    "Tangram": "Tangram",
    "Seurat": "Seurat",
    "CARD": "CARD",
}

# Module colors for pipeline diagrams
MODULE_COLORS = {
    1: "#2171b5",  # Blue - Marker Detection
    2: "#41ab5d",  # Green - Profile Discovery
    3: "#c51b8a",  # Magenta - Deconvolution
    4: "#fe9929",  # Orange - Program Discovery
    5: "#6baed6",  # Light blue - Integration
}


def apply_style(journal: str = "nature"):
    """
    Apply publication-quality style settings to matplotlib.

    Call this at the top of each figure generation script to ensure
    consistent styling across all figures.

    Font sizes are sourced from the journal font profile (default: "nature").
    Panels are saved at a fixed canvas size; bbox is not adjusted.

    Parameters
    ----------
    journal : str
        Journal profile to use for font sizes ("nature", "science", "cell").
        Defaults to "nature".
    """
    profile = resolve_font_profile({"journal": journal})
    plt.rcParams.update(
        {
            "font.family": _FONT_FAMILY,
            "font.size": profile["axis_label_pt"],
            "axes.labelsize": profile["axis_label_pt"],
            "axes.titlesize": profile["title_pt"],
            "xtick.labelsize": profile["tick_label_pt"],
            "ytick.labelsize": profile["tick_label_pt"],
            "legend.fontsize": profile["tick_label_pt"],
            "legend.title_fontsize": profile["axis_label_pt"],
            # Figure layout
            "figure.constrained_layout.use": True,
            "figure.dpi": 300,
            "savefig.dpi": 300,
            "savefig.bbox": None,
            "savefig.pad_inches": 0,
            "savefig.transparent": False,
            "savefig.facecolor": "white",
            # Axes styling
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.linewidth": 0.8,
            "axes.labelpad": 4,
            "axes.titlepad": 8,
            # Grid (off by default for clean look)
            "axes.grid": False,
            # Tick styling
            "xtick.major.width": 0.8,
            "ytick.major.width": 0.8,
            "xtick.major.size": 4,
            "ytick.major.size": 4,
            "xtick.direction": "out",
            "ytick.direction": "out",
            # Legend styling
            "legend.frameon": True,
            "legend.framealpha": 0.9,
            "legend.edgecolor": "none",
            "legend.borderpad": 0.4,
            # Line styling
            "lines.linewidth": 1.5,
            "lines.markersize": 6,
            # Patch styling
            "patch.linewidth": 1,
            # SVG output — <text> elements instead of <path> glyphs (required for font enforcement)
            "svg.fonttype": "none",
        }
    )


def get_cell_type_color(cell_type: str) -> str:
    """
    Get color for a cell type, with fallback to neutral gray.

    Args:
        cell_type: Cell type name (handles various naming conventions)

    Returns:
        Hex color string
    """
    # Direct lookup
    if cell_type in CELL_TYPE_COLORS:
        return CELL_TYPE_COLORS[cell_type]

    # Try with underscores replaced by spaces
    normalized = cell_type.replace("_", " ")
    if normalized in CELL_TYPE_COLORS:
        return CELL_TYPE_COLORS[normalized]

    # Try partial match
    for key, color in CELL_TYPE_COLORS.items():
        if key.lower() in cell_type.lower() or cell_type.lower() in key.lower():
            return color

    # Fallback to neutral
    return PALETTE["neutral"]


def get_method_color(method: str) -> str:
    """
    Get color for a deconvolution method, with fallback.

    Args:
        method: Method name

    Returns:
        Hex color string
    """
    # Handle suffixes like "_achievable_7"
    base_method = method.split("_")[0]

    if base_method in METHOD_COLORS:
        return METHOD_COLORS[base_method]
    if method in METHOD_COLORS:
        return METHOD_COLORS[method]

    return PALETTE["neutral"]


# Convenience function for creating consistent figure sizes
def get_figure_size(panels="2x2"):
    """
    Get standard figure sizes for common panel layouts.

    Args:
        panels: Layout descriptor ('1x1', '1x2', '2x1', '2x2', '3x1', etc.)

    Returns:
        tuple: (width, height) in inches
    """
    sizes = {
        "1x1": (5, 4),
        "1x2": (10, 4),
        "2x1": (6, 8),
        "2x2": (10, 8),
        "3x1": (5, 10),
        "1x3": (12, 4),
        "2x3": (12, 8),
        "3x2": (10, 10),
        "full_width": (7, 5),  # Cell Reports full width
        "half_width": (3.5, 3.5),  # Cell Reports half width
    }
    return sizes.get(panels, (10, 8))


if __name__ == "__main__":
    # Test the style module
    apply_style()
    print("Style applied successfully.")
    print(f"PALETTE keys: {list(PALETTE.keys())}")
    print(f"CELL_TYPE_COLORS: {len(CELL_TYPE_COLORS)} cell types")
    print(f"METHOD_COLORS: {list(METHOD_COLORS.keys())}")
