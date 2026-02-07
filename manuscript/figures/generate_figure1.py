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

# Paths
PROJECT_ROOT = Path(__file__).parent.parent.parent
OUTPUT_DIR = Path(__file__).parent / "output"

# Ensure output directory exists
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Style settings
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.size'] = 9
plt.rcParams['axes.linewidth'] = 0.5
plt.rcParams['xtick.major.width'] = 0.5
plt.rcParams['ytick.major.width'] = 0.5

# Module colors
MODULE_COLORS = {
    1: '#3498db',  # Blue - Marker Detection
    2: '#2ecc71',  # Green - Profile Discovery
    3: '#e74c3c',  # Red - Deconvolution
    4: '#9b59b6',  # Purple - Program Discovery
    5: '#f39c12',  # Orange - Integration
}


def draw_module_box(ax, x, y, width, height, module_num, title, subtitle, color):
    """Draw a module box with number, title, and subtitle."""
    # Main box
    box = FancyBboxPatch(
        (x, y), width, height,
        boxstyle="round,pad=0.02,rounding_size=0.03",
        facecolor=color, edgecolor='black', linewidth=1.5, alpha=0.85
    )
    ax.add_patch(box)

    # Module number circle
    circle = Circle((x + 0.08, y + height - 0.08), 0.05,
                    facecolor='white', edgecolor=color, linewidth=1.5)
    ax.add_patch(circle)
    ax.text(x + 0.08, y + height - 0.08, str(module_num),
            ha='center', va='center', fontsize=10, fontweight='bold', color=color)

    # Title
    ax.text(x + width/2, y + height - 0.15, title,
            ha='center', va='top', fontsize=9, fontweight='bold', color='white')

    # Subtitle
    ax.text(x + width/2, y + 0.08, subtitle,
            ha='center', va='bottom', fontsize=7, color='white', style='italic')


def draw_arrow(ax, start, end, color='black'):
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

    # Layout: 5 boxes in a row with arrows
    box_width = 0.15
    box_height = 0.35
    start_x = 0.05
    spacing = 0.17
    y_pos = 0.35

    for i, (mod_num, title, subtitle) in enumerate(modules):
        x = start_x + i * spacing
        draw_module_box(ax, x, y_pos, box_width, box_height,
                       mod_num, title, subtitle, MODULE_COLORS[mod_num])

        # Draw arrow to next module
        if i < len(modules) - 1:
            arrow_start = (x + box_width + 0.01, y_pos + box_height/2)
            arrow_end = (x + spacing - 0.01, y_pos + box_height/2)
            draw_arrow(ax, arrow_start, arrow_end, color='#2c3e50')

    # Input/Output labels
    ax.text(0.02, y_pos + box_height/2, "Spatial\nTranscriptomics\n+ CITE-seq",
            ha='center', va='center', fontsize=7, style='italic',
            bbox=dict(boxstyle='round', facecolor='#ecf0f1', edgecolor='gray'))

    ax.text(0.98, y_pos + box_height/2, "Spatial\nPrograms\n+ Relationships",
            ha='center', va='center', fontsize=7, style='italic',
            bbox=dict(boxstyle='round', facecolor='#ecf0f1', edgecolor='gray'))

    # Title
    ax.text(0.5, 0.95, "A. CITEgeist Modular Pipeline",
            ha='center', va='top', fontsize=12, fontweight='bold')

    # Subtitle
    ax.text(0.5, 0.15, "Protein-anchored deconvolution with automatic profile discovery",
            ha='center', va='top', fontsize=9, style='italic', color='#7f8c8d')


def panel_b_spatial_operations(ax):
    """Panel B: Spatial operations highlight."""
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis('off')

    # Title
    ax.text(0.5, 0.95, "B. Spatial Statistics Foundation",
            ha='center', va='top', fontsize=11, fontweight='bold')

    # Section 1: Moran's I formula
    ax.text(0.17, 0.80, "Moran's I", ha='center', va='top', fontsize=9, fontweight='bold')
    morans_formula = r"$I = \frac{n}{\sum_{ij}w_{ij}} \frac{\sum_{ij}w_{ij}(x_i-\bar{x})(x_j-\bar{x})}{\sum_i(x_i-\bar{x})^2}$"
    ax.text(0.17, 0.65, morans_formula, ha='center', va='center', fontsize=10)
    ax.text(0.17, 0.50, "Spatial autocorrelation\n(marker clustering)",
            ha='center', va='top', fontsize=7, style='italic', color='#7f8c8d')

    # Section 2: Colocalization concept
    ax.text(0.50, 0.80, "Colocalization", ha='center', va='top', fontsize=9, fontweight='bold')

    # Draw two overlapping circles (Venn-like)
    circle1 = Circle((0.45, 0.60), 0.08, facecolor='#3498db', alpha=0.5, edgecolor='none')
    circle2 = Circle((0.55, 0.60), 0.08, facecolor='#e74c3c', alpha=0.5, edgecolor='none')
    ax.add_patch(circle1)
    ax.add_patch(circle2)
    ax.text(0.42, 0.60, "M1", ha='center', va='center', fontsize=8, fontweight='bold')
    ax.text(0.58, 0.60, "M2", ha='center', va='center', fontsize=8, fontweight='bold')

    ax.text(0.50, 0.45, "Bivariate Moran's I\nCo-occurrence scoring",
            ha='center', va='top', fontsize=7, style='italic', color='#7f8c8d')

    # Section 3: Neighborhood graph
    ax.text(0.83, 0.80, "Spatial Graph", ha='center', va='top', fontsize=9, fontweight='bold')

    # Draw a simple neighborhood graph (hexagonal pattern)
    hex_centers = [
        (0.83, 0.62),
        (0.76, 0.58), (0.90, 0.58),
        (0.76, 0.50), (0.83, 0.50), (0.90, 0.50),
    ]

    for cx, cy in hex_centers:
        circle = Circle((cx, cy), 0.025, facecolor='#95a5a6', edgecolor='#2c3e50', linewidth=0.5)
        ax.add_patch(circle)

    # Draw edges (simplified)
    edges = [
        ((0.83, 0.62), (0.76, 0.58)), ((0.83, 0.62), (0.90, 0.58)),
        ((0.83, 0.62), (0.83, 0.50)),
        ((0.76, 0.58), (0.76, 0.50)), ((0.76, 0.58), (0.83, 0.50)),
        ((0.90, 0.58), (0.90, 0.50)), ((0.90, 0.58), (0.83, 0.50)),
        ((0.76, 0.50), (0.83, 0.50)), ((0.83, 0.50), (0.90, 0.50)),
    ]
    for start, end in edges:
        ax.plot([start[0], end[0]], [start[1], end[1]],
                color='#2c3e50', linewidth=0.8, alpha=0.6)

    ax.text(0.83, 0.42, "k-NN / Radius\nNeighbor weighting",
            ha='center', va='top', fontsize=7, style='italic', color='#7f8c8d')

    # Bottom annotation
    ax.text(0.5, 0.25,
            "All modules leverage spatial context through neighbor-weighted statistics",
            ha='center', va='top', fontsize=8, style='italic',
            bbox=dict(boxstyle='round', facecolor='#f8f9fa', edgecolor='#dee2e6'))


def panel_c_resolution_flexibility(ax):
    """Panel C: Resolution flexibility diagram (Visium spot vs single cell)."""
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis('off')

    # Title
    ax.text(0.5, 0.95, "C. Resolution-Agnostic Design",
            ha='center', va='top', fontsize=11, fontweight='bold')

    # Left side: Visium spot resolution
    ax.text(0.25, 0.78, "Spot Resolution", ha='center', va='top',
            fontsize=9, fontweight='bold', color='#2980b9')
    ax.text(0.25, 0.70, "(Visium, ~10-30 cells/spot)", ha='center', va='top',
            fontsize=7, style='italic', color='#7f8c8d')

    # Draw Visium-like spots (large circles with mixed colors)
    spot_centers = [(0.15, 0.50), (0.25, 0.55), (0.35, 0.48), (0.20, 0.40)]
    spot_colors = ['#e74c3c', '#3498db', '#2ecc71', '#9b59b6']
    for (cx, cy), color in zip(spot_centers, spot_colors):
        # Outer spot
        circle = Circle((cx, cy), 0.06, facecolor=color, alpha=0.3,
                        edgecolor=color, linewidth=1.5)
        ax.add_patch(circle)
        # Inner cells (small dots)
        for _ in range(4):
            dx, dy = np.random.uniform(-0.03, 0.03, 2)
            cell = Circle((cx + dx, cy + dy), 0.012,
                         facecolor=np.random.choice(['#e74c3c', '#3498db', '#2ecc71', '#f39c12']),
                         alpha=0.8, edgecolor='none')
            ax.add_patch(cell)

    ax.text(0.25, 0.28, "Mixed cell types\nper spot", ha='center', va='top',
            fontsize=7, style='italic', color='#7f8c8d')

    # Arrow in center
    ax.annotate('', xy=(0.60, 0.50), xytext=(0.40, 0.50),
                arrowprops=dict(arrowstyle='<->', color='#2c3e50', lw=2))
    ax.text(0.50, 0.58, "Same\nAlgorithm", ha='center', va='bottom',
            fontsize=8, fontweight='bold', color='#2c3e50')

    # Right side: Single-cell resolution
    ax.text(0.75, 0.78, "Cell Resolution", ha='center', va='top',
            fontsize=9, fontweight='bold', color='#27ae60')
    ax.text(0.75, 0.70, "(Xenium, CosMx, MERFISH)", ha='center', va='top',
            fontsize=7, style='italic', color='#7f8c8d')

    # Draw single cells (individual colored circles)
    np.random.seed(42)
    cell_types = ['#e74c3c', '#3498db', '#2ecc71', '#9b59b6', '#f39c12']
    for _ in range(25):
        cx = np.random.uniform(0.62, 0.88)
        cy = np.random.uniform(0.35, 0.62)
        color = np.random.choice(cell_types)
        cell = Circle((cx, cy), 0.018, facecolor=color, alpha=0.8,
                      edgecolor='white', linewidth=0.3)
        ax.add_patch(cell)

    ax.text(0.75, 0.28, "Direct cell type\nobservation", ha='center', va='top',
            fontsize=7, style='italic', color='#7f8c8d')

    # Legend at bottom
    legend_items = [
        ('T cells', '#e74c3c'),
        ('Macrophages', '#3498db'),
        ('Fibroblasts', '#2ecc71'),
        ('Epithelial', '#9b59b6'),
        ('Endothelial', '#f39c12'),
    ]

    for i, (label, color) in enumerate(legend_items):
        x = 0.12 + i * 0.18
        circle = Circle((x, 0.12), 0.015, facecolor=color, edgecolor='none')
        ax.add_patch(circle)
        ax.text(x + 0.025, 0.12, label, ha='left', va='center', fontsize=6)


def generate_figure1():
    """Generate complete Figure 1."""
    print("Generating Figure 1: Pipeline Overview...")

    # Create figure with 3 panels (A on top spanning full width, B and C below)
    fig = plt.figure(figsize=(12, 10))
    gs = GridSpec(2, 2, figure=fig, height_ratios=[1, 1.2], hspace=0.25, wspace=0.15)

    # Panel A: Pipeline overview (top, full width)
    ax_a = fig.add_subplot(gs[0, :])
    panel_a_pipeline_overview(ax_a)

    # Panel B: Spatial operations (bottom left)
    ax_b = fig.add_subplot(gs[1, 0])
    panel_b_spatial_operations(ax_b)

    # Panel C: Resolution flexibility (bottom right)
    ax_c = fig.add_subplot(gs[1, 1])
    panel_c_resolution_flexibility(ax_c)

    # Add figure label
    fig.text(0.02, 0.98, "Figure 1", fontsize=14, fontweight='bold', va='top')

    # Save
    output_path = OUTPUT_DIR / "figure1_pipeline_overview.pdf"
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"Saved to {output_path}")

    # Also save PNG for quick preview
    png_path = OUTPUT_DIR / "figure1_pipeline_overview.png"
    plt.savefig(png_path, dpi=150, bbox_inches='tight', facecolor='white')
    print(f"Preview saved to {png_path}")

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
    print("\nNote: This is a schematic figure. Final version should be")
    print("refined in vector graphics software (Illustrator, Inkscape).")


if __name__ == "__main__":
    generate_figure1()
