#!/usr/bin/env python3
"""
Figure 6: Interpretability & External Tool Compatibility

Panels:
  A: CITEgeist outputs -> external tools workflow diagram
     - Show: CITEgeist outputs (proportions, deconvolved GEX, programs) -> downstream tools
     - Tools: PyDESeq2, GSEApy/Enrichr, Harmony, COMMOT, scanpy (Leiden, rank_genes_groups)

  B: Midkine discovery pathway (condensed summary)
     - Key findings from the midkine/ER saturation mechanism
     - Show: CITEgeist discovery -> COMMOT cell-cell signaling -> Wet lab validation
     - Reference midkine analysis outputs

  C: Wet lab validation summary
     - MCF7 vs T47D differential response
     - ER saturation model concept
     - ChIP-seq binding validation mention

This is primarily a conceptual/summary figure with schematic diagrams.

Data sources:
  - Midkine mechanism: midkine/FINAL_MECHANISM_SUMMARY.md
  - Pipeline outputs: midkine/mdk_saturation_pipeline/outputs/
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Circle, Rectangle, Polygon
from matplotlib.gridspec import GridSpec
from matplotlib.collections import PatchCollection
from pathlib import Path

# Paths
PROJECT_ROOT = Path(__file__).parent.parent.parent
MIDKINE_DIR = PROJECT_ROOT / "midkine"
OUTPUT_DIR = Path(__file__).parent / "output"

# Ensure output directory exists
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Style settings
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.size'] = 8
plt.rcParams['axes.linewidth'] = 0.5
plt.rcParams['xtick.major.width'] = 0.5
plt.rcParams['ytick.major.width'] = 0.5

# Color palette
COLORS = {
    'citegeist_blue': '#3498db',
    'citegeist_dark': '#2980b9',
    'output_green': '#27ae60',
    'tool_orange': '#f39c12',
    'tool_purple': '#9b59b6',
    'tool_red': '#e74c3c',
    'arrow_gray': '#7f8c8d',
    'background_light': '#ecf0f1',
    'mcf7_color': '#e74c3c',
    't47d_color': '#3498db',
    'validation_green': '#27ae60',
}


def draw_rounded_box(ax, x, y, width, height, text, color, fontsize=8,
                      textcolor='white', subtitle=None, icon=None):
    """Draw a rounded box with text."""
    box = FancyBboxPatch(
        (x, y), width, height,
        boxstyle="round,pad=0.02,rounding_size=0.02",
        facecolor=color, edgecolor='black', linewidth=1, alpha=0.9
    )
    ax.add_patch(box)

    # Main text
    y_text = y + height/2 if subtitle is None else y + height*0.65
    ax.text(x + width/2, y_text, text,
            ha='center', va='center', fontsize=fontsize,
            fontweight='bold', color=textcolor)

    # Subtitle if provided
    if subtitle:
        ax.text(x + width/2, y + height*0.30, subtitle,
                ha='center', va='center', fontsize=fontsize-2,
                color=textcolor, style='italic')


def draw_arrow(ax, start, end, color='#7f8c8d', style='->', connectionstyle=None):
    """Draw an arrow between two points."""
    arrowprops = dict(
        arrowstyle=style,
        color=color,
        linewidth=1.5,
        shrinkA=3,
        shrinkB=3
    )
    if connectionstyle:
        arrowprops['connectionstyle'] = connectionstyle

    ax.annotate('', xy=end, xytext=start, arrowprops=arrowprops)


def panel_a_interoperability_workflow(ax):
    """Panel A: CITEgeist outputs to external tools workflow."""
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis('off')

    # Title
    ax.text(0.5, 0.97, "A. Interoperability with External Analysis Tools",
            ha='center', va='top', fontsize=11, fontweight='bold')

    # === Left side: CITEgeist Outputs ===
    ax.text(0.12, 0.85, "CITEgeist Outputs", ha='center', va='center',
            fontsize=9, fontweight='bold', color=COLORS['citegeist_dark'])

    # Output boxes
    outputs = [
        ("Cell Type\nProportions", 0.72, "Module 3"),
        ("Deconvolved\nGEX Layers", 0.54, "Module 3"),
        ("Spatial\nPrograms", 0.36, "Module 4"),
        ("Program\nRelationships", 0.18, "Module 5"),
    ]

    for text, y, module in outputs:
        draw_rounded_box(ax, 0.02, y, 0.20, 0.12, text,
                        COLORS['citegeist_blue'], fontsize=7, subtitle=module)

    # === Center: Arrow labels ===
    ax.text(0.35, 0.65, "AnnData\nLayers", ha='center', va='center',
            fontsize=7, style='italic', color=COLORS['arrow_gray'])
    ax.text(0.35, 0.35, "CSV/\nParquet", ha='center', va='center',
            fontsize=7, style='italic', color=COLORS['arrow_gray'])

    # === Right side: External Tools ===
    ax.text(0.72, 0.85, "Downstream Analysis", ha='center', va='center',
            fontsize=9, fontweight='bold', color=COLORS['tool_purple'])

    # Tool boxes - organized by analysis type
    tools = [
        # (name, y, color, description)
        ("PyDESeq2", 0.72, COLORS['tool_orange'], "Differential\nExpression"),
        ("GSEApy/\nEnrichr", 0.54, COLORS['tool_purple'], "Pathway\nEnrichment"),
        ("scanpy", 0.36, COLORS['tool_red'], "Clustering\n& Markers"),
        ("COMMOT", 0.18, COLORS['validation_green'], "Cell-Cell\nSignaling"),
    ]

    for name, y, color, desc in tools:
        draw_rounded_box(ax, 0.55, y, 0.18, 0.12, name, color, fontsize=7)
        ax.text(0.85, y + 0.06, desc, ha='left', va='center',
                fontsize=6, style='italic', color='#2c3e50')

    # Draw arrows from outputs to tools
    # Proportions -> scanpy (clustering by composition)
    draw_arrow(ax, (0.22, 0.78), (0.55, 0.42), color=COLORS['arrow_gray'])
    # GEX Layers -> PyDESeq2, GSEApy
    draw_arrow(ax, (0.22, 0.60), (0.55, 0.78), color=COLORS['arrow_gray'])
    draw_arrow(ax, (0.22, 0.60), (0.55, 0.60), color=COLORS['arrow_gray'])
    # Programs -> scanpy
    draw_arrow(ax, (0.22, 0.42), (0.55, 0.42), color=COLORS['arrow_gray'])
    # Relationships -> COMMOT
    draw_arrow(ax, (0.22, 0.24), (0.55, 0.24), color=COLORS['arrow_gray'])

    # Additional integration note
    ax.text(0.50, 0.02,
            "All outputs stored as AnnData layers or standard tabular formats for seamless integration",
            ha='center', va='bottom', fontsize=7, style='italic',
            bbox=dict(boxstyle='round', facecolor=COLORS['background_light'],
                     edgecolor='#bdc3c7', alpha=0.8))


def panel_b_midkine_discovery(ax):
    """Panel B: Midkine discovery pathway - condensed summary."""
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis('off')

    # Title
    ax.text(0.50, 0.95, "B. Discovery Pipeline: Midkine Secretion Mechanism",
            ha='center', va='top', fontsize=10, fontweight='bold')

    # === Discovery Flow (left to right) ===
    flow_y = 0.70

    # Step 1: CITEgeist Discovery
    draw_rounded_box(ax, 0.02, flow_y - 0.08, 0.22, 0.16,
                    "CITEgeist\nModule 4", COLORS['citegeist_blue'],
                    fontsize=8, subtitle="Program Discovery")

    # Arrow 1
    draw_arrow(ax, (0.24, flow_y), (0.30, flow_y), color=COLORS['arrow_gray'])

    # Step 2: COMMOT Analysis
    draw_rounded_box(ax, 0.30, flow_y - 0.08, 0.18, 0.16,
                    "COMMOT", COLORS['tool_purple'],
                    fontsize=8, subtitle="Signaling")

    # Arrow 2
    draw_arrow(ax, (0.48, flow_y), (0.54, flow_y), color=COLORS['arrow_gray'])

    # Step 3: MDK Hypothesis
    draw_rounded_box(ax, 0.54, flow_y - 0.08, 0.18, 0.16,
                    "MDK\nHypothesis", COLORS['tool_orange'],
                    fontsize=8, subtitle="Secretion")

    # Arrow 3
    draw_arrow(ax, (0.72, flow_y), (0.78, flow_y), color=COLORS['arrow_gray'])

    # Step 4: Validation
    draw_rounded_box(ax, 0.78, flow_y - 0.08, 0.20, 0.16,
                    "Wet Lab\nValidation", COLORS['validation_green'],
                    fontsize=8, subtitle="ChIP-seq")

    # === Key Finding Box ===
    finding_box = FancyBboxPatch(
        (0.05, 0.25), 0.90, 0.30,
        boxstyle="round,pad=0.02,rounding_size=0.02",
        facecolor='#f8f9fa', edgecolor='#2c3e50', linewidth=1.5
    )
    ax.add_patch(finding_box)

    ax.text(0.50, 0.52, "Key Discovery", ha='center', va='center',
            fontsize=9, fontweight='bold', color='#2c3e50')

    # Discovery text
    discovery_text = (
        "CITEgeist identified MDK (Midkine) as a spatially-variable gene\n"
        "in ER+ breast cancer. COMMOT analysis revealed cell-cell signaling\n"
        "patterns. Follow-up wet lab experiments discovered the ER saturation\n"
        "model explaining opposite MDK secretion in MCF7 vs T47D cell lines."
    )
    ax.text(0.50, 0.38, discovery_text, ha='center', va='center',
            fontsize=7.5, color='#2c3e50', linespacing=1.4)

    # Bottom note
    ax.text(0.50, 0.08,
            "Demonstrates CITEgeist's ability to generate testable biological hypotheses",
            ha='center', va='center', fontsize=7, style='italic', color='#7f8c8d')


def panel_c_validation_summary(ax):
    """Panel C: Wet lab validation summary - ER saturation model."""
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis('off')

    # Title
    ax.text(0.50, 0.95, "C. Validation: ER Saturation Model",
            ha='center', va='top', fontsize=10, fontweight='bold')

    # === MCF7 Box (Left) ===
    mcf7_box = FancyBboxPatch(
        (0.03, 0.35), 0.42, 0.52,
        boxstyle="round,pad=0.02,rounding_size=0.02",
        facecolor='#fdf2f2', edgecolor=COLORS['mcf7_color'], linewidth=2
    )
    ax.add_patch(mcf7_box)

    ax.text(0.24, 0.82, "MCF7 (ER-Saturated)", ha='center', va='center',
            fontsize=9, fontweight='bold', color=COLORS['mcf7_color'])

    mcf7_details = [
        "ESR1: 58.8 TPM (HIGH)",
        "ER Occupancy: 3.0%",
        "ER Sites: 12,472 peaks",
        "",
        "D538G Effect:",
        "  Loses binding sites",
        "  Chaperones UP",
        "  MDK secretion UP",
    ]

    for i, line in enumerate(mcf7_details):
        y_pos = 0.72 - i * 0.05
        weight = 'bold' if 'Effect' in line or line.startswith('  MDK') else 'normal'
        color = COLORS['mcf7_color'] if 'UP' in line else '#2c3e50'
        ax.text(0.06, y_pos, line, ha='left', va='center',
                fontsize=7, fontweight=weight, color=color)

    # === T47D Box (Right) ===
    t47d_box = FancyBboxPatch(
        (0.55, 0.35), 0.42, 0.52,
        boxstyle="round,pad=0.02,rounding_size=0.02",
        facecolor='#f0f7ff', edgecolor=COLORS['t47d_color'], linewidth=2
    )
    ax.add_patch(t47d_box)

    ax.text(0.76, 0.82, "T47D (ER-Unsaturated)", ha='center', va='center',
            fontsize=9, fontweight='bold', color=COLORS['t47d_color'])

    t47d_details = [
        "ESR1: 12.8 TPM (LOW)",
        "ER Occupancy: 0.9%",
        "ER Sites: 1,724 peaks",
        "",
        "D538G Effect:",
        "  Gains binding sites",
        "  Chaperones DOWN",
        "  MDK secretion DOWN",
    ]

    for i, line in enumerate(t47d_details):
        y_pos = 0.72 - i * 0.05
        weight = 'bold' if 'Effect' in line or line.startswith('  MDK') else 'normal'
        color = COLORS['t47d_color'] if 'DOWN' in line else '#2c3e50'
        ax.text(0.58, y_pos, line, ha='left', va='center',
                fontsize=7, fontweight=weight, color=color)

    # === VS Arrow between boxes ===
    ax.text(0.50, 0.60, "vs", ha='center', va='center',
            fontsize=10, fontweight='bold', color='#7f8c8d')

    # === Bottom: Key Insight ===
    insight_box = FancyBboxPatch(
        (0.05, 0.05), 0.90, 0.22,
        boxstyle="round,pad=0.02,rounding_size=0.02",
        facecolor=COLORS['background_light'], edgecolor='#7f8c8d', linewidth=1
    )
    ax.add_patch(insight_box)

    ax.text(0.50, 0.22, "Mechanism Validated by:", ha='center', va='center',
            fontsize=8, fontweight='bold', color='#2c3e50')

    validations = [
        "RNA-seq (GSE89888): Opposite chaperone regulation",
        "ChIP-seq (GSE125117): Opposite binding changes at HSP90B1",
        "ATAC-seq (GSE254216): Differential chromatin saturation",
    ]

    for i, val in enumerate(validations):
        ax.text(0.10, 0.15 - i * 0.04, f"  {val}", ha='left', va='center',
                fontsize=6.5, color='#2c3e50')


def panel_d_data_formats(ax):
    """Panel D: Data formats and file types output by CITEgeist."""
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis('off')

    # Title
    ax.text(0.50, 0.95, "D. Output Formats for Downstream Analysis",
            ha='center', va='top', fontsize=10, fontweight='bold')

    # Create table-like display
    headers = ["Output Type", "Format", "Key Contents", "Compatible Tools"]

    data = [
        ["Proportions", ".csv, AnnData", "Cell type fractions per spot", "scanpy, seaborn"],
        ["GEX Layers", ".h5ad, .parquet", "Per-cell-type expression", "PyDESeq2, GSEA"],
        ["Programs", ".csv, AnnData", "Gene loadings, H matrices", "scanpy, NMF tools"],
        ["Relationships", ".csv, .json", "Bivariate Moran's I", "networkx, Cytoscape"],
        ["Spatial Coords", "AnnData.obsm", "x, y coordinates", "squidpy, COMMOT"],
    ]

    # Draw header
    y_start = 0.78
    col_x = [0.05, 0.22, 0.42, 0.72]

    for i, (x, header) in enumerate(zip(col_x, headers)):
        ax.text(x, y_start, header, ha='left', va='center',
                fontsize=8, fontweight='bold', color='#2c3e50')

    # Horizontal line under header
    ax.plot([0.03, 0.97], [y_start - 0.04, y_start - 0.04],
            color='#2c3e50', linewidth=1)

    # Draw data rows
    for row_idx, row in enumerate(data):
        y = y_start - 0.12 - row_idx * 0.12

        # Alternating row background
        if row_idx % 2 == 0:
            bg = Rectangle((0.03, y - 0.05), 0.94, 0.10,
                           facecolor='#f8f9fa', edgecolor='none')
            ax.add_patch(bg)

        for col_idx, (x, cell) in enumerate(zip(col_x, row)):
            fontweight = 'bold' if col_idx == 0 else 'normal'
            color = COLORS['citegeist_dark'] if col_idx == 0 else '#2c3e50'
            ax.text(x, y, cell, ha='left', va='center',
                    fontsize=7, fontweight=fontweight, color=color)

    # Bottom note
    ax.text(0.50, 0.08,
            "Standard formats ensure compatibility with existing bioinformatics ecosystems",
            ha='center', va='center', fontsize=7, style='italic', color='#7f8c8d')


def generate_figure6():
    """Generate complete Figure 6."""
    print("Generating Figure 6: Interpretability & External Tool Compatibility...")

    # Create figure with 2x2 layout
    fig = plt.figure(figsize=(14, 12))
    gs = GridSpec(2, 2, figure=fig, height_ratios=[1, 1],
                  hspace=0.25, wspace=0.20)

    # Panel A: Interoperability workflow (top left)
    ax_a = fig.add_subplot(gs[0, 0])
    panel_a_interoperability_workflow(ax_a)

    # Panel B: Midkine discovery (top right)
    ax_b = fig.add_subplot(gs[0, 1])
    panel_b_midkine_discovery(ax_b)

    # Panel C: Validation summary (bottom left)
    ax_c = fig.add_subplot(gs[1, 0])
    panel_c_validation_summary(ax_c)

    # Panel D: Data formats (bottom right)
    ax_d = fig.add_subplot(gs[1, 1])
    panel_d_data_formats(ax_d)

    # Add figure label
    fig.text(0.02, 0.98, "Figure 6", fontsize=14, fontweight='bold', va='top')

    # Save PDF
    output_path = OUTPUT_DIR / "figure6_interoperability.pdf"
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"Saved to {output_path}")

    # Save PNG for preview
    png_path = OUTPUT_DIR / "figure6_interoperability.png"
    plt.savefig(png_path, dpi=150, bbox_inches='tight', facecolor='white')
    print(f"Preview saved to {png_path}")

    plt.close()

    # Print summary
    print("\n=== Figure 6 Summary ===")
    print("\nPanel A: Interoperability with External Analysis Tools")
    print("  - CITEgeist outputs: Proportions, GEX Layers, Programs, Relationships")
    print("  - Downstream tools: PyDESeq2, GSEApy/Enrichr, scanpy, COMMOT")
    print("  - Standard formats: AnnData, CSV, Parquet")

    print("\nPanel B: Midkine Discovery Pipeline")
    print("  - CITEgeist Module 4 -> COMMOT -> MDK Hypothesis -> Wet Lab Validation")
    print("  - Demonstrates hypothesis generation capability")

    print("\nPanel C: ER Saturation Model Validation")
    print("  - MCF7 (ER-saturated): Loses binding -> Chaperones UP -> MDK UP")
    print("  - T47D (ER-unsaturated): Gains binding -> Chaperones DOWN -> MDK DOWN")
    print("  - Validated by: RNA-seq (GSE89888), ChIP-seq (GSE125117), ATAC-seq (GSE254216)")

    print("\nPanel D: Output Formats for Downstream Analysis")
    print("  - Standard formats ensure ecosystem compatibility")
    print("  - Supports: scanpy, PyDESeq2, GSEA, networkx, squidpy")

    print("\nNote: This is primarily a conceptual/schematic figure.")
    print("Final version should be refined in vector graphics software if needed.")


if __name__ == "__main__":
    generate_figure6()
