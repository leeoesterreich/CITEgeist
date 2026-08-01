#!/usr/bin/env python3
"""
Shared style configuration for CITEgeist manuscript figures.

This module provides consistent styling across all figures for Cell/Cell Reports
publication quality. All figure generation scripts should import and use this.

Usage:
    from figure_style import apply_style, PALETTE, CELL_TYPE_COLORS, METHOD_COLORS
    apply_style()
"""

import matplotlib.pyplot as plt

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
}

# Module colors for pipeline diagrams
MODULE_COLORS = {
    1: "#2171b5",  # Blue - Marker Detection
    2: "#41ab5d",  # Green - Profile Discovery
    3: "#c51b8a",  # Magenta - Deconvolution
    4: "#fe9929",  # Orange - Program Discovery
    5: "#6baed6",  # Light blue - Integration
}


def apply_style():
    """
    Apply publication-quality style settings to matplotlib.

    Call this at the top of each figure generation script to ensure
    consistent styling across all figures.

    Font sizes are set to minimum 10pt for readability in print.
    """
    plt.rcParams.update(
        {
            # Font settings (DejaVu Sans is universally available, similar to Arial)
            "font.family": "sans-serif",
            "font.size": 10,
            "axes.labelsize": 10,
            "axes.titlesize": 11,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "legend.fontsize": 9,
            "legend.title_fontsize": 10,
            # Figure layout
            "figure.constrained_layout.use": True,
            "figure.dpi": 300,
            "savefig.dpi": 300,
            "savefig.bbox": "tight",
            "savefig.pad_inches": 0.1,
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
