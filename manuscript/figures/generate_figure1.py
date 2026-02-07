#!/usr/bin/env python3
"""
Figure 1: CITEgeist Pipeline Overview

Panels:
  A: Module flow schematic (1 -> 2 -> 3 -> 4 -> 5) with labeled boxes
  B: Spatial operations highlight (Moran's I formula, colocalization concept, neighborhood graph)
  C: Resolution flexibility diagram (Visium spot vs single cell)

This is primarily a conceptual/schematic figure with text placeholders.
Final vector graphics should be created in a dedicated editor (Illustrator, Inkscape).
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Circle, Polygon
from matplotlib.gridspec import GridSpec
from pathlib import Path

# Import shared style
from figure_style import apply_style, PALETTE, CELL_TYPE_COLORS, MODULE_COLORS

# Apply publication style
apply_style()

# Paths
PROJECT_ROOT = Path(__file__).parent.parent.parent
OUTPUT_DIR = Path(__file__).parent / "output"

# Ensure output directory exists
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def draw_module_box(ax, x, y, width, height, module_num, title, subtitle, color):
    """Draw a module box with number, title, and subtitle."""
    # Main box
    box = FancyBboxPatch(
        (x, y), width, height,
        boxstyle="round,pad=0.02,rounding_size=0.03",
        facecolor=color, edgecolor='#2c3e50', linewidth=1.5, alpha=0.9
    )
    ax.add_patch(box)

    # Module number circle
    circle = Circle((x + 0.08, y + height - 0.08), 0.05,
                    facecolor='white', edgecolor=color, linewidth=2)
    ax.add_patch(circle)
    ax.text(x + 0.08, y + height - 0.08, str(module_num),
            ha='center', va='center', fontsize=11, fontweight='bold', color=color)

    # Title - increased font size
    ax.text(x + width/2, y + height - 0.14, title,
            ha='center', va='top', fontsize=10, fontweight='bold', color='white')

    # Subtitle - increased font size
    ax.text(x + width/2, y + 0.08, subtitle,
            ha='center', va='bottom', fontsize=8, color='white', style='italic')


def draw_arrow(ax, start, end, color='#2c3e50'):
    """Draw an arrow between two points."""
    arrow = FancyArrowPatch(
        start, end,
        arrowstyle='-|>', mutation_scale=15,
        color=color, linewidth=2
    )
    ax.add_patch(arrow)


def panel_a_pipeline_overview(ax):
    """Panel A: Module flow schematic (1 -> 2 -> 3 -> 4 -> 5)."""
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis('off')

    # Module specifications
    modules = [
        (1, "Marker\nInterest", "Kurtosis/GMM/Moran's I"),
        (2, "Profile\nDiscovery", "Colocalization + Clustering"),
        (3, "Deconvolution", "Proportions + GEX"),
        (4, "Program\nDiscovery", "NMF + Spatial Coherence"),
        (5, "Integration", "Cross-sample Alignment"),
    ]

    # Layout: 5 boxes in a row with arrows - tighter spacing
    box_width = 0.14
    box_height = 0.38
    start_x = 0.08
    spacing = 0.165
    y_pos = 0.32

    for i, (mod_num, title, subtitle) in enumerate(modules):
        x = start_x + i * spacing
        draw_module_box(ax, x, y_pos, box_width, box_height,
                       mod_num, title, subtitle, MODULE_COLORS[mod_num])

        # Draw arrow to next module
        if i < len(modules) - 1:
            arrow_start = (x + box_width + 0.005, y_pos + box_height/2)
            arrow_end = (x + spacing - 0.005, y_pos + box_height/2)
            draw_arrow(ax, arrow_start, arrow_end)

    # Input/Output labels - increased font size
    ax.text(0.03, y_pos + box_height/2, "Spatial\nTranscriptomics\n+ CITE-seq",
            ha='center', va='center', fontsize=9, style='italic',
            bbox=dict(boxstyle='round,pad=0.3', facecolor=PALETTE['background'],
                     edgecolor=PALETTE['border'], linewidth=1))

    ax.text(0.97, y_pos + box_height/2, "Spatial\nPrograms\n+ Relationships",
            ha='center', va='center', fontsize=9, style='italic',
            bbox=dict(boxstyle='round,pad=0.3', facecolor=PALETTE['background'],
                     edgecolor=PALETTE['border'], linewidth=1))

    # Title - panel label
    ax.text(0.02, 0.95, "A", fontsize=14, fontweight='bold', va='top')
    ax.text(0.5, 0.95, "CITEgeist Modular Pipeline",
            ha='center', va='top', fontsize=12, fontweight='bold')

    # Subtitle - increased font size
    ax.text(0.5, 0.18, "Protein-anchored deconvolution with automatic profile discovery",
            ha='center', va='top', fontsize=10, style='italic', color=PALETTE['neutral'])


def panel_b_spatial_operations(ax):
    """Panel B: Spatial operations highlight."""
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis('off')

    # Title
    ax.text(0.02, 0.95, "B", fontsize=14, fontweight='bold', va='top')
    ax.text(0.5, 0.95, "Spatial Statistics Foundation",
            ha='center', va='top', fontsize=12, fontweight='bold')

    # Section 1: Moran's I formula
    ax.text(0.17, 0.82, "Moran's I", ha='center', va='top', fontsize=10, fontweight='bold')
    morans_formula = r"$I = \frac{n}{\sum_{ij}w_{ij}} \frac{\sum_{ij}w_{ij}(x_i-\bar{x})(x_j-\bar{x})}{\sum_i(x_i-\bar{x})^2}$"
    ax.text(0.17, 0.68, morans_formula, ha='center', va='center', fontsize=11)
    ax.text(0.17, 0.52, "Spatial autocorrelation\n(marker clustering)",
            ha='center', va='top', fontsize=9, style='italic', color=PALETTE['neutral'])

    # Section 2: Colocalization concept
    ax.text(0.50, 0.82, "Colocalization", ha='center', va='top', fontsize=10, fontweight='bold')

    # Draw two overlapping circles (Venn-like) using CELL_TYPE_COLORS
    circle1 = Circle((0.45, 0.62), 0.08, facecolor=PALETTE['primary'], alpha=0.5, edgecolor='none')
    circle2 = Circle((0.55, 0.62), 0.08, facecolor=PALETTE['highlight'], alpha=0.5, edgecolor='none')
    ax.add_patch(circle1)
    ax.add_patch(circle2)
    ax.text(0.42, 0.62, "M1", ha='center', va='center', fontsize=9, fontweight='bold')
    ax.text(0.58, 0.62, "M2", ha='center', va='center', fontsize=9, fontweight='bold')

    ax.text(0.50, 0.47, "Bivariate Moran's I\nCo-occurrence scoring",
            ha='center', va='top', fontsize=9, style='italic', color=PALETTE['neutral'])

    # Section 3: Neighborhood graph
    ax.text(0.83, 0.82, "Spatial Graph", ha='center', va='top', fontsize=10, fontweight='bold')

    # Draw a simple neighborhood graph (hexagonal pattern)
    hex_centers = [
        (0.83, 0.64),
        (0.76, 0.60), (0.90, 0.60),
        (0.76, 0.52), (0.83, 0.52), (0.90, 0.52),
    ]

    for cx, cy in hex_centers:
        circle = Circle((cx, cy), 0.028, facecolor=PALETTE['secondary'],
                        edgecolor='#2c3e50', linewidth=0.8)
        ax.add_patch(circle)

    # Draw edges (simplified)
    edges = [
        ((0.83, 0.64), (0.76, 0.60)), ((0.83, 0.64), (0.90, 0.60)),
        ((0.83, 0.64), (0.83, 0.52)),
        ((0.76, 0.60), (0.76, 0.52)), ((0.76, 0.60), (0.83, 0.52)),
        ((0.90, 0.60), (0.90, 0.52)), ((0.90, 0.60), (0.83, 0.52)),
        ((0.76, 0.52), (0.83, 0.52)), ((0.83, 0.52), (0.90, 0.52)),
    ]
    for start, end in edges:
        ax.plot([start[0], end[0]], [start[1], end[1]],
                color='#2c3e50', linewidth=1, alpha=0.6)

    ax.text(0.83, 0.42, "k-NN / Radius\nNeighbor weighting",
            ha='center', va='top', fontsize=9, style='italic', color=PALETTE['neutral'])

    # Bottom annotation - tighter box
    ax.text(0.5, 0.22,
            "All modules leverage spatial context through neighbor-weighted statistics",
            ha='center', va='top', fontsize=10, style='italic',
            bbox=dict(boxstyle='round,pad=0.4', facecolor=PALETTE['background'],
                     edgecolor=PALETTE['border'], linewidth=1))


def panel_c_resolution_flexibility(ax):
    """Panel C: Resolution flexibility diagram (Visium spot vs single cell)."""
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis('off')

    # Title
    ax.text(0.02, 0.95, "C", fontsize=14, fontweight='bold', va='top')
    ax.text(0.5, 0.95, "Resolution-Agnostic Design",
            ha='center', va='top', fontsize=12, fontweight='bold')

    # Left side: Visium spot resolution
    ax.text(0.25, 0.80, "Spot Resolution", ha='center', va='top',
            fontsize=10, fontweight='bold', color=PALETTE['primary'])
    ax.text(0.25, 0.72, "(Visium, ~10-30 cells/spot)", ha='center', va='top',
            fontsize=9, style='italic', color=PALETTE['neutral'])

    # Draw Visium-like spots (large circles with mixed colors)
    spot_centers = [(0.15, 0.52), (0.25, 0.56), (0.35, 0.50), (0.20, 0.42)]
    spot_colors = [PALETTE['highlight'], PALETTE['primary'], PALETTE['accent2'], PALETTE['accent1']]
    np.random.seed(42)
    for (cx, cy), color in zip(spot_centers, spot_colors):
        # Outer spot
        circle = Circle((cx, cy), 0.06, facecolor=color, alpha=0.3,
                        edgecolor=color, linewidth=1.5)
        ax.add_patch(circle)
        # Inner cells (small dots)
        for _ in range(4):
            dx, dy = np.random.uniform(-0.03, 0.03, 2)
            cell = Circle((cx + dx, cy + dy), 0.012,
                         facecolor=np.random.choice(spot_colors),
                         alpha=0.8, edgecolor='none')
            ax.add_patch(cell)

    ax.text(0.25, 0.30, "Mixed cell types\nper spot", ha='center', va='top',
            fontsize=9, style='italic', color=PALETTE['neutral'])

    # Arrow in center
    ax.annotate('', xy=(0.60, 0.50), xytext=(0.40, 0.50),
                arrowprops=dict(arrowstyle='<->', color='#2c3e50', lw=2.5))
    ax.text(0.50, 0.58, "Same\nAlgorithm", ha='center', va='bottom',
            fontsize=10, fontweight='bold', color='#2c3e50')

    # Right side: Single-cell resolution
    ax.text(0.75, 0.80, "Cell Resolution", ha='center', va='top',
            fontsize=10, fontweight='bold', color=PALETTE['accent2'])
    ax.text(0.75, 0.72, "(Xenium, CosMx, MERFISH)", ha='center', va='top',
            fontsize=9, style='italic', color=PALETTE['neutral'])

    # Draw single cells (individual colored circles)
    np.random.seed(42)
    cell_colors = [PALETTE['highlight'], PALETTE['primary'], PALETTE['accent2'],
                   PALETTE['accent1'], PALETTE['neutral']]
    for _ in range(25):
        cx = np.random.uniform(0.62, 0.88)
        cy = np.random.uniform(0.38, 0.65)
        color = np.random.choice(cell_colors)
        cell = Circle((cx, cy), 0.018, facecolor=color, alpha=0.8,
                      edgecolor='white', linewidth=0.4)
        ax.add_patch(cell)

    ax.text(0.75, 0.30, "Direct cell type\nobservation", ha='center', va='top',
            fontsize=9, style='italic', color=PALETTE['neutral'])

    # Legend at bottom - larger font, better spacing
    legend_items = [
        ('T cells', CELL_TYPE_COLORS['T cells']),
        ('Macrophages', CELL_TYPE_COLORS['Macrophages']),
        ('Fibroblasts', CELL_TYPE_COLORS['Fibroblasts']),
        ('Epithelial', CELL_TYPE_COLORS['Epithelial']),
        ('Endothelial', CELL_TYPE_COLORS['Endothelial']),
    ]

    for i, (label, color) in enumerate(legend_items):
        x = 0.10 + i * 0.18
        circle = Circle((x, 0.12), 0.018, facecolor=color, edgecolor='none')
        ax.add_patch(circle)
        ax.text(x + 0.03, 0.12, label, ha='left', va='center', fontsize=9)


def generate_figure1():
    """Generate complete Figure 1."""
    print("Generating Figure 1: Pipeline Overview...")

    # Create figure with 3 panels (A on top spanning full width, B and C below)
    # Reduced figure size for tighter layout
    fig = plt.figure(figsize=(11, 9))
    gs = GridSpec(2, 2, figure=fig, height_ratios=[1, 1.1], hspace=0.20, wspace=0.12)

    # Panel A: Pipeline overview (top, full width)
    ax_a = fig.add_subplot(gs[0, :])
    panel_a_pipeline_overview(ax_a)

    # Panel B: Spatial operations (bottom left)
    ax_b = fig.add_subplot(gs[1, 0])
    panel_b_spatial_operations(ax_b)

    # Panel C: Resolution flexibility (bottom right)
    ax_c = fig.add_subplot(gs[1, 1])
    panel_c_resolution_flexibility(ax_c)

    # Save
    output_path = OUTPUT_DIR / "figure1_pipeline_overview.pdf"
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"Saved to {output_path}")

    # Also save PNG for quick preview
    png_path = OUTPUT_DIR / "figure1_pipeline_overview.png"
    plt.savefig(png_path, dpi=150, bbox_inches='tight', facecolor='white')
    print(f"Preview saved to {png_path}")

    # Save SVG for Illustrator
    svg_path = OUTPUT_DIR / "figure1_pipeline_overview.svg"
    plt.savefig(svg_path, format='svg', bbox_inches='tight', facecolor='white')
    print(f"SVG saved to {svg_path}")

    plt.close()

    # Print summary
    print("\n=== Figure 1 Summary ===")
    print("Panel A: CITEgeist modular pipeline (Modules 1-5)")
    print("  - Module 1: Marker Interest Detection (Kurtosis/GMM/Moran's I)")
    print("  - Module 2: Profile Discovery (Colocalization + Hierarchical Clustering)")
    print("  - Module 3: Deconvolution (Cell Proportions + Gene Expression)")
    print("  - Module 4: Program Discovery (NMF + Spatial Coherence)")
    print("  - Module 5: Cross-sample Integration (Alignment + Conserved Relationships)")
    print("\nPanel B: Spatial statistics foundation")
    print("  - Moran's I formula for spatial autocorrelation")
    print("  - Colocalization via bivariate statistics")
    print("  - Neighborhood graph construction")
    print("\nPanel C: Resolution flexibility")
    print("  - Visium spot resolution (mixed cells)")
    print("  - Single-cell resolution (Xenium/CosMx/MERFISH)")
    print("  - Same algorithm works at both scales")
    print("\nEnhancements applied:")
    print("  - Fonts increased to minimum 10pt")
    print("  - Consistent color palette from figure_style.py")
    print("  - Tightened layout with reduced spacing")


if __name__ == "__main__":
    generate_figure1()
