#!/usr/bin/env python3
"""
Figure 2: Modules 1-2 Profile Discovery

Panels A, B, C: SCHEMATICS - need dedicated SVG files
Panel D: DATA TABLE - generated with matplotlib

Note: Figure 2 schematics are complex and not yet in svg_schematics.py.
For now, this script generates Panel D only.
Create A, B, C manually in Illustrator using consistent styling.
"""

import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from matplotlib.gridspec import GridSpec
from pathlib import Path

# Import shared style
from figure_style import apply_style, PALETTE, CELL_TYPE_COLORS

apply_style()

# Paths
PROJECT_ROOT = Path(__file__).parent.parent.parent
CELL_RES_DIR = PROJECT_ROOT / "Benchmarking/xenium_benchmarking/CITEgeist/output_cell_resolution"
OUTPUT_DIR = Path(__file__).parent / "output"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def load_cell_resolution_profiles(region_id=0):
    """Load discovered profiles from cell resolution analysis."""
    filepath = CELL_RES_DIR / f"region_{region_id}/profiles.json"
    if filepath.exists():
        with open(filepath, 'r') as f:
            return json.load(f)
    return None


def panel_d_profile_table(ax):
    """Panel D: Discovered profiles vs known markers table."""
    ax.text(0.02, 0.98, "D", fontsize=14, fontweight='bold', va='top')
    ax.axis('off')

    profiles = load_cell_resolution_profiles(0)

    # Known cell type markers (ground truth)
    known_markers = {
        'B cells': ['CD20'],
        'CD4+ T cells': ['CD3E', 'CD4'],
        'CD8+ T cells': ['CD3E', 'CD8A'],
        'Macrophages': ['CD68', 'CD163'],
        'Endothelial': ['CD31'],
        'Epithelial': ['PanCK'],
        'Fibroblasts': ['alphaSMA'],
    }

    cell_types = list(known_markers.keys())
    y_start = 0.88
    row_height = 0.10
    col_widths = [0.25, 0.35, 0.30]
    col_starts = [0.05, 0.32, 0.68]

    # Header
    headers = ['Cell Type', 'Known Markers', 'Discovered']
    for j, (header, col_start, col_width) in enumerate(zip(headers, col_starts, col_widths)):
        box = FancyBboxPatch((col_start, y_start), col_width - 0.02, row_height - 0.02,
                            boxstyle="round,pad=0.005", facecolor=PALETTE['primary'],
                            edgecolor='none', alpha=0.9)
        ax.add_patch(box)
        ax.text(col_start + col_width/2 - 0.01, y_start + row_height/2 - 0.01,
                header, ha='center', va='center', fontsize=10, fontweight='bold', color='white')

    # Data rows
    for i, cell_type in enumerate(cell_types):
        y = y_start - (i + 1) * row_height
        bg_color = PALETTE['background'] if i % 2 == 0 else 'white'

        for j, (col_start, col_width) in enumerate(zip(col_starts, col_widths)):
            box = FancyBboxPatch((col_start, y), col_width - 0.02, row_height - 0.02,
                                boxstyle="round,pad=0.005", facecolor=bg_color,
                                edgecolor=PALETTE['border'], linewidth=0.5)
            ax.add_patch(box)

        # Cell type name
        ax.text(col_starts[0] + col_widths[0]/2 - 0.01, y + row_height/2 - 0.01,
                cell_type, ha='center', va='center', fontsize=9, fontweight='bold')

        # Known markers
        known = known_markers[cell_type]
        ax.text(col_starts[1] + col_widths[1]/2 - 0.01, y + row_height/2 - 0.01,
                ', '.join(known), ha='center', va='center', fontsize=9)

        # Discovered markers
        if profiles:
            discovered = profiles.get(cell_type, [])
            if set(discovered) == set(known):
                status_color = PALETTE['accent2']
                status_symbol = "="
            elif set(known).issubset(set(discovered)):
                status_color = PALETTE['accent1']
                status_symbol = "+"
            else:
                status_color = PALETTE['neutral']
                status_symbol = "?"
            text_val = ', '.join(discovered) if discovered else '-'
            ax.text(col_starts[2] + col_widths[2]/2 - 0.01, y + row_height/2 - 0.01,
                    f"{status_symbol} {text_val}", ha='center', va='center', fontsize=9, color=status_color)
        else:
            ax.text(col_starts[2] + col_widths[2]/2 - 0.01, y + row_height/2 - 0.01,
                    "N/A", ha='center', va='center', fontsize=9, color=PALETTE['neutral'])

    # Legend
    ax.text(0.15, 0.08, "=", ha='center', va='center', fontsize=11, fontweight='bold', color=PALETTE['accent2'])
    ax.text(0.22, 0.08, "Exact match", ha='left', va='center', fontsize=9)
    ax.text(0.45, 0.08, "+", ha='center', va='center', fontsize=11, fontweight='bold', color=PALETTE['accent1'])
    ax.text(0.52, 0.08, "Superset", ha='left', va='center', fontsize=9)
    ax.text(0.75, 0.08, "?", ha='center', va='center', fontsize=11, fontweight='bold', color=PALETTE['neutral'])
    ax.text(0.82, 0.08, "Not found", ha='left', va='center', fontsize=9)


def generate_figure2():
    """Generate Figure 2 with placeholders for schematics."""
    print("Generating Figure 2...")

    fig = plt.figure(figsize=(12, 10))
    gs = GridSpec(2, 2, figure=fig, hspace=0.22, wspace=0.18)

    # Panel A: Placeholder
    ax_a = fig.add_subplot(gs[0, 0])
    ax_a.text(0.5, 0.5, "Panel A: Schematic\n\nMarker Interest Detection\n(Create in Illustrator)",
              ha='center', va='center', fontsize=10, style='italic',
              bbox=dict(boxstyle='round', facecolor='#f0f0f0', edgecolor='gray'))
    ax_a.set_xlim(0, 1)
    ax_a.set_ylim(0, 1)
    ax_a.axis('off')
    ax_a.set_title("Module 1: Marker Interest Detection", fontsize=12, fontweight='bold', loc='left')

    # Panel B: Placeholder
    ax_b = fig.add_subplot(gs[0, 1])
    ax_b.text(0.5, 0.5, "Panel B: Schematic\n\nProfile Discovery Workflow\n(Create in Illustrator)",
              ha='center', va='center', fontsize=10, style='italic',
              bbox=dict(boxstyle='round', facecolor='#f0f0f0', edgecolor='gray'))
    ax_b.set_xlim(0, 1)
    ax_b.set_ylim(0, 1)
    ax_b.axis('off')
    ax_b.set_title("Module 2: Profile Discovery Workflow", fontsize=12, fontweight='bold', loc='left')

    # Panel C: Placeholder
    ax_c = fig.add_subplot(gs[1, 0])
    ax_c.text(0.5, 0.5, "Panel C: Schematic\n\nXenium RCC Demonstration\n(Create in Illustrator)",
              ha='center', va='center', fontsize=10, style='italic',
              bbox=dict(boxstyle='round', facecolor='#f0f0f0', edgecolor='gray'))
    ax_c.set_xlim(0, 1)
    ax_c.set_ylim(0, 1)
    ax_c.axis('off')
    ax_c.set_title("Xenium RCC Single-Cell Demonstration", fontsize=12, fontweight='bold', loc='left')

    # Panel D: Profile comparison table (DATA)
    ax_d = fig.add_subplot(gs[1, 1])
    ax_d.set_xlim(0, 1)
    ax_d.set_ylim(0, 1)
    panel_d_profile_table(ax_d)
    ax_d.set_title("Discovered Profiles vs Known Markers", fontsize=12, fontweight='bold', loc='left')

    # Save
    for fmt, dpi in [('pdf', 300), ('png', 150), ('svg', None)]:
        output_path = OUTPUT_DIR / f"figure2_profile_discovery.{fmt}"
        if fmt == 'svg':
            plt.savefig(output_path, format='svg', bbox_inches='tight', facecolor='white')
        else:
            plt.savefig(output_path, dpi=dpi, bbox_inches='tight', facecolor='white')
        print(f"Saved: {output_path}")

    plt.close()

    print("\n" + "=" * 60)
    print("Figure 2: Profile Discovery")
    print("=" * 60)
    print("\nPanels A, B, C: SCHEMATICS - create in Illustrator")
    print("Panel D: DATA TABLE - generated above")
    print("\nFigure 2 schematics are complex. Recommended workflow:")
    print("  1. Use figure2_profile_discovery.svg as template")
    print("  2. Replace placeholders A, B, C with custom schematics")
    print("  3. Export final figure")


if __name__ == "__main__":
    generate_figure2()
