#!/usr/bin/env python3
"""
Figure 1: CITEgeist Pipeline Overview

Assembles complete figure from pre-rendered SVG schematic PNGs.
- Panel A: CITEgeist Modular Pipeline
- Panel B: Spatial Statistics Foundation
- Panel C: Resolution-Agnostic Design
"""

import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.gridspec import GridSpec
from pathlib import Path

from figure_style import apply_style

apply_style()

# Paths
OUTPUT_DIR = Path(__file__).parent / "output"
SCHEMATIC_DIR = OUTPUT_DIR / "schematics" / "rendered"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def load_schematic(filename):
    """Load a rendered schematic PNG."""
    filepath = SCHEMATIC_DIR / filename
    if filepath.exists():
        return mpimg.imread(str(filepath))
    else:
        print(f"WARNING: {filename} not found")
        return None


def generate_figure1():
    """Generate Figure 1 by assembling schematic panels."""
    print("Loading rendered schematics...")

    panel_a = load_schematic("figure1_panel_a_pipeline.png")
    panel_b = load_schematic("figure1_panel_b_spatial_stats.png")
    panel_c = load_schematic("figure1_panel_c_resolution.png")

    # Create figure with 2 rows: Panel A on top (full width), B and C below
    fig = plt.figure(figsize=(12, 10))
    gs = GridSpec(2, 2, figure=fig, height_ratios=[1.2, 1], hspace=0.15, wspace=0.10)

    # Panel A: Full width top row
    ax_a = fig.add_subplot(gs[0, :])
    if panel_a is not None:
        ax_a.imshow(panel_a)
    ax_a.axis('off')
    ax_a.text(-0.02, 1.02, "A", fontsize=16, fontweight='bold', va='top',
              transform=ax_a.transAxes)

    # Panel B: Bottom left
    ax_b = fig.add_subplot(gs[1, 0])
    if panel_b is not None:
        ax_b.imshow(panel_b)
    ax_b.axis('off')
    ax_b.text(-0.02, 1.02, "B", fontsize=16, fontweight='bold', va='top',
              transform=ax_b.transAxes)

    # Panel C: Bottom right
    ax_c = fig.add_subplot(gs[1, 1])
    if panel_c is not None:
        ax_c.imshow(panel_c)
    ax_c.axis('off')
    ax_c.text(-0.02, 1.02, "C", fontsize=16, fontweight='bold', va='top',
              transform=ax_c.transAxes)

    plt.tight_layout()

    # Save outputs
    for fmt, dpi in [('pdf', 300), ('png', 150), ('svg', None)]:
        output_path = OUTPUT_DIR / f"figure1_pipeline_overview.{fmt}"
        if fmt == 'svg':
            plt.savefig(output_path, format='svg', bbox_inches='tight', facecolor='white')
        else:
            plt.savefig(output_path, dpi=dpi, bbox_inches='tight', facecolor='white')
        print(f"Saved: {output_path}")

    plt.close()

    print("\n" + "=" * 60)
    print("Figure 1: CITEgeist Pipeline Overview - COMPLETE")
    print("=" * 60)


if __name__ == "__main__":
    generate_figure1()
