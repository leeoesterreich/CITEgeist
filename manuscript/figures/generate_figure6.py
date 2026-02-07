#!/usr/bin/env python3
"""
Figure 6: Interoperability & Output Formats

Enhanced version with consistent styling.

Panels:
  A: CITEgeist → External tools workflow
  B: Output formats table
"""

import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle
from matplotlib.gridspec import GridSpec
from pathlib import Path

# Import shared style
from figure_style import apply_style, PALETTE

# Apply publication style
apply_style()

OUTPUT_DIR = Path(__file__).parent / "output"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Tool-specific colors for interoperability diagram
TOOL_COLORS = {
    'citegeist': PALETTE['primary'],      # Deep blue
    'scanpy': PALETTE['accent2'],          # Green
    'pydeseq2': PALETTE['accent1'],        # Orange
    'gseapy': PALETTE['highlight'],        # Magenta
    'commot': '#e74c3c',                   # Red
}


def panel_a_workflow(ax):
    """Panel A: Interoperability workflow diagram."""
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis('off')

    # Panel label
    ax.text(0.02, 0.98, "A", fontsize=14, fontweight='bold', va='top')

    arrow_color = PALETTE['neutral']

    # CITEgeist outputs (left column)
    outputs = [
        ("Proportions", 0.72, TOOL_COLORS['citegeist']),
        ("GEX Layers", 0.52, TOOL_COLORS['citegeist']),
        ("Programs", 0.32, TOOL_COLORS['citegeist']),
        ("Relationships", 0.12, TOOL_COLORS['citegeist']),
    ]

    for label, y, color in outputs:
        rect = FancyBboxPatch((0.02, y), 0.18, 0.15, boxstyle="round,pad=0.01",
                              facecolor=color, edgecolor='#2980b9', linewidth=2, alpha=0.8)
        ax.add_patch(rect)
        ax.text(0.11, y + 0.075, label, ha='center', va='center', fontsize=10, fontweight='bold', color='white')

    ax.text(0.11, 0.92, "CITEgeist\nOutputs", ha='center', va='bottom', fontsize=11, fontweight='bold')

    # Arrows
    for y in [0.79, 0.59, 0.39, 0.19]:
        ax.annotate('', xy=(0.28, y), xytext=(0.21, y),
                    arrowprops=dict(arrowstyle='->', color=arrow_color, lw=1.5))

    # Format boxes (middle)
    ax.text(0.35, 0.67, "AnnData\nCSV", ha='center', va='center', fontsize=10, style='italic', color=PALETTE['neutral'])
    ax.text(0.35, 0.32, "Parquet\nJSON", ha='center', va='center', fontsize=10, style='italic', color=PALETTE['neutral'])

    # Arrows to tools
    for y in [0.79, 0.59, 0.39, 0.19]:
        ax.annotate('', xy=(0.48, y), xytext=(0.42, y),
                    arrowprops=dict(arrowstyle='->', color=arrow_color, lw=1.5))

    # External tools (right column)
    tools = [
        ("scanpy", 0.72, TOOL_COLORS['scanpy']),
        ("PyDESeq2", 0.52, TOOL_COLORS['pydeseq2']),
        ("GSEApy", 0.32, TOOL_COLORS['gseapy']),
        ("COMMOT", 0.12, TOOL_COLORS['commot']),
    ]

    for label, y, color in tools:
        rect = FancyBboxPatch((0.48, y), 0.18, 0.15, boxstyle="round,pad=0.01",
                              facecolor=color, edgecolor='black', linewidth=1.5, alpha=0.85)
        ax.add_patch(rect)
        ax.text(0.57, y + 0.075, label, ha='center', va='center', fontsize=10, fontweight='bold', color='white')

    ax.text(0.57, 0.92, "Downstream\nTools", ha='center', va='bottom', fontsize=11, fontweight='bold')

    # Use cases (far right)
    uses = [
        ("Clustering\nVisualization", 0.72),
        ("Differential\nExpression", 0.52),
        ("Pathway\nEnrichment", 0.32),
        ("Cell-Cell\nCommunication", 0.12),
    ]

    for label, y in uses:
        ax.text(0.78, y + 0.075, label, ha='left', va='center', fontsize=10, color=PALETTE['neutral'])


def panel_b_formats(ax):
    """Panel B: Output formats table."""
    ax.axis('off')

    # Panel label
    ax.text(0.02, 0.98, "B", fontsize=14, fontweight='bold', va='top')

    headers = ["Output", "Formats", "Contents", "Tools"]
    data = [
        ["Proportions", ".csv, .h5ad", "Cell fractions/spot", "scanpy, seaborn"],
        ["GEX Layers", ".h5ad, .parquet", "Per-cell-type expr", "PyDESeq2, GSEA"],
        ["Programs", ".csv, .h5ad", "W, H matrices", "NMF, scanpy"],
        ["Relationships", ".csv, .json", "Bivariate I", "networkx"],
        ["Coordinates", "AnnData.obsm", "x, y positions", "squidpy, COMMOT"],
    ]

    # Table
    table = ax.table(
        cellText=data,
        colLabels=headers,
        cellLoc='left',
        loc='center',
        bbox=[0.05, 0.15, 0.9, 0.75]
    )
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1, 2)

    # Style header with PALETTE primary color
    for j in range(len(headers)):
        table[(0, j)].set_facecolor(PALETTE['primary'])
        table[(0, j)].set_text_props(color='white', fontweight='bold')

    # Alternate row colors
    for i in range(1, len(data) + 1):
        for j in range(len(headers)):
            if i % 2 == 0:
                table[(i, j)].set_facecolor(PALETTE['background'])

    ax.text(0.5, 0.02, "Standard formats enable seamless integration with bioinformatics ecosystems",
            ha='center', va='bottom', fontsize=9, style='italic', color=PALETTE['neutral'])


def generate_figure6():
    """Generate complete Figure 6."""
    print("Generating Figure 6...")

    fig = plt.figure(figsize=(11, 7))
    gs = GridSpec(2, 1, figure=fig, height_ratios=[1.2, 1], hspace=0.25)

    # Panel A: Workflow
    ax_a = fig.add_subplot(gs[0])
    panel_a_workflow(ax_a)
    ax_a.set_title("Interoperability with External Tools", fontsize=12, fontweight='bold', loc='left', pad=10)

    # Panel B: Formats
    ax_b = fig.add_subplot(gs[1])
    panel_b_formats(ax_b)
    ax_b.set_title("Output Formats", fontsize=12, fontweight='bold', loc='left', pad=10)

    plt.tight_layout()

    output_path = OUTPUT_DIR / "figure6_interoperability.pdf"
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"Saved to {output_path}")

    png_path = OUTPUT_DIR / "figure6_interoperability.png"
    plt.savefig(png_path, dpi=150, bbox_inches='tight', facecolor='white')
    print(f"Preview saved to {png_path}")

    # Save SVG for Illustrator
    svg_path = OUTPUT_DIR / "figure6_interoperability.svg"
    plt.savefig(svg_path, format='svg', bbox_inches='tight', facecolor='white')
    print(f"SVG saved to {svg_path}")

    plt.close()

    print("\n=== Figure 6 Summary ===")
    print("Panels generated:")
    print("  A: Interoperability workflow (CITEgeist → External tools)")
    print("  B: Output formats table (5 output types)")

    print("\nEnhancements applied:")
    print("  - Consistent color palette from figure_style.py")
    print("  - Fonts increased to minimum 10pt")
    print("  - Panel labels in top-left corner")
    print("  - TOOL_COLORS for downstream tool distinction")


if __name__ == "__main__":
    generate_figure6()
