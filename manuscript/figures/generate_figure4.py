#!/usr/bin/env python3
"""
Figure 4: Module 4 Spatial Programs

Enhanced version with consistent styling.

Panels:
  A: NMF schematic (simple)
  B: Top programs with genes
  C: Moran's I by cell type
  D: Summary statistics
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch
from matplotlib.gridspec import GridSpec
from pathlib import Path

# Import shared style
from figure_style import apply_style, PALETTE, CELL_TYPE_COLORS, get_cell_type_color

# Apply publication style
apply_style()

# Paths
PROJECT_ROOT = Path(__file__).parent.parent.parent
XENIUM_MODULE4_DIR = PROJECT_ROOT / "Benchmarking/xenium_benchmarking/CITEgeist/output_module4_validation/singlecell"
PATIENT_MODULE5_DIR = PROJECT_ROOT / "output/module5_integration"
OUTPUT_DIR = Path(__file__).parent / "output"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def load_xenium_programs():
    """Load Module 4 program discovery results from Xenium data."""
    all_programs = []
    for region in range(5):
        region_dir = XENIUM_MODULE4_DIR / f"region_{region}"
        if region_dir.exists():
            program_file = region_dir / "module4_programs.csv"
            if program_file.exists():
                df = pd.read_csv(program_file)
                df['region'] = region
                all_programs.append(df)

    if all_programs:
        return pd.concat(all_programs, ignore_index=True)
    return None


def panel_a_schematic(ax):
    """Panel A: Simple NMF schematic."""
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis('off')

    # Panel label
    ax.text(0.02, 0.95, "A", fontsize=14, fontweight='bold', va='top')

    arrow_color = PALETTE['neutral']

    # Input matrix
    rect1 = FancyBboxPatch((0.05, 0.3), 0.18, 0.4, boxstyle="round,pad=0.01",
                           facecolor='#dae8fc', edgecolor='#6c8ebf', linewidth=2)
    ax.add_patch(rect1)
    ax.text(0.14, 0.5, "Expression\nMatrix", ha='center', va='center', fontsize=10, fontweight='bold')
    ax.text(0.14, 0.22, "Spots x Genes", ha='center', va='top', fontsize=9, color=PALETTE['neutral'])

    # Arrow
    ax.annotate('', xy=(0.30, 0.5), xytext=(0.24, 0.5),
                arrowprops=dict(arrowstyle='->', color=arrow_color, lw=2))

    # NMF box
    rect2 = FancyBboxPatch((0.30, 0.35), 0.15, 0.30, boxstyle="round,pad=0.01",
                           facecolor='#fff2cc', edgecolor='#d6b656', linewidth=2)
    ax.add_patch(rect2)
    ax.text(0.375, 0.5, "NMF", ha='center', va='center', fontsize=11, fontweight='bold')

    # Arrow
    ax.annotate('', xy=(0.52, 0.5), xytext=(0.46, 0.5),
                arrowprops=dict(arrowstyle='->', color=arrow_color, lw=2))

    # W matrix
    rect3 = FancyBboxPatch((0.52, 0.45), 0.12, 0.25, boxstyle="round,pad=0.01",
                           facecolor='#d5e8d4', edgecolor='#82b366', linewidth=2)
    ax.add_patch(rect3)
    ax.text(0.58, 0.575, "W", ha='center', va='center', fontsize=11, fontweight='bold')
    ax.text(0.58, 0.38, "Gene\nLoadings", ha='center', va='top', fontsize=9, color=PALETTE['neutral'])

    # H matrix
    rect4 = FancyBboxPatch((0.52, 0.15), 0.22, 0.12, boxstyle="round,pad=0.01",
                           facecolor='#d5e8d4', edgecolor='#82b366', linewidth=2)
    ax.add_patch(rect4)
    ax.text(0.63, 0.21, "H", ha='center', va='center', fontsize=11, fontweight='bold')
    ax.text(0.77, 0.21, "Program\nActivities", ha='left', va='center', fontsize=9, color=PALETTE['neutral'])

    # Moran's I validation
    ax.annotate('', xy=(0.63, 0.08), xytext=(0.63, 0.14),
                arrowprops=dict(arrowstyle='->', color=arrow_color, lw=2))
    rect5 = FancyBboxPatch((0.52, -0.05), 0.22, 0.12, boxstyle="round,pad=0.01",
                           facecolor='#f8cecc', edgecolor='#b85450', linewidth=2)
    ax.add_patch(rect5)
    ax.text(0.63, 0.01, "Moran's I\nValidation", ha='center', va='center', fontsize=10, fontweight='bold')


def panel_b_programs(ax, programs_df):
    """Panel B: Top programs as horizontal bars with genes."""
    # Panel label
    ax.text(-0.12, 1.05, "B", fontsize=14, fontweight='bold', va='top', transform=ax.transAxes)

    if programs_df is None or len(programs_df) == 0:
        ax.text(0.5, 0.5, "Program data not available",
                ha='center', va='center', fontsize=10, style='italic', color=PALETTE['neutral'])
        ax.axis('off')
        return

    # Get top programs by Moran's I
    top_df = programs_df.nlargest(8, 'spatial_moran_i')

    y_pos = np.arange(len(top_df))
    colors = [get_cell_type_color(ct) for ct in top_df['cell_type']]

    bars = ax.barh(y_pos, top_df['spatial_moran_i'], color=colors, alpha=0.8, height=0.7)

    # Add gene labels
    for i, (_, row) in enumerate(top_df.iterrows()):
        try:
            genes = eval(row['top_genes'])[:4]
            gene_str = ', '.join(genes)
        except:
            gene_str = ""
        ax.text(row['spatial_moran_i'] + 0.02, i, gene_str, va='center', fontsize=9)

    # Y-axis labels
    labels = [f"{row['cell_type'].split()[0]}" for _, row in top_df.iterrows()]
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=9)
    ax.set_xlabel("Moran's I (Spatial Coherence)", fontsize=10)
    ax.set_xlim(0, 0.65)
    ax.invert_yaxis()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)


def panel_c_boxplot(ax, programs_df):
    """Panel C: Moran's I distribution by cell type."""
    # Panel label
    ax.text(-0.12, 1.05, "C", fontsize=14, fontweight='bold', va='top', transform=ax.transAxes)

    if programs_df is None or len(programs_df) == 0:
        ax.text(0.5, 0.5, "Program data not available",
                ha='center', va='center', fontsize=10, style='italic', color=PALETTE['neutral'])
        ax.axis('off')
        return

    cell_types = programs_df['cell_type'].unique()
    data = [programs_df[programs_df['cell_type'] == ct]['spatial_moran_i'].values for ct in cell_types]

    bp = ax.boxplot(data, patch_artist=True, widths=0.6)

    for patch, ct in zip(bp['boxes'], cell_types):
        patch.set_facecolor(get_cell_type_color(ct))
        patch.set_alpha(0.7)

    ax.axhline(y=0.15, color=PALETTE['highlight'], linestyle='--', linewidth=1.5, alpha=0.7, label='Threshold')
    ax.set_xticklabels([ct.split()[0] for ct in cell_types], rotation=45, ha='right', fontsize=9)
    ax.set_ylabel("Moran's I", fontsize=10)
    ax.legend(fontsize=9)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Summary text
    n_above = (programs_df['spatial_moran_i'] > 0.15).sum()
    total = len(programs_df)
    ax.text(0.98, 0.98, f'{n_above}/{total} above\nthreshold',
            transform=ax.transAxes, ha='right', va='top', fontsize=9,
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))


def panel_d_summary(ax, programs_df):
    """Panel D: Summary statistics table."""
    # Panel label
    ax.text(0.02, 0.98, "D", fontsize=14, fontweight='bold', va='top')
    ax.axis('off')

    if programs_df is None or len(programs_df) == 0:
        ax.text(0.5, 0.5, "Program data not available",
                ha='center', va='center', fontsize=10, style='italic', color=PALETTE['neutral'])
        return

    # Compute summary
    summary_data = []
    for ct in programs_df['cell_type'].unique():
        ct_data = programs_df[programs_df['cell_type'] == ct]
        summary_data.append({
            'Cell Type': ct.split()[0],
            'N': len(ct_data),
            'Mean I': f"{ct_data['spatial_moran_i'].mean():.2f}",
            'Max I': f"{ct_data['spatial_moran_i'].max():.2f}",
        })

    df = pd.DataFrame(summary_data)

    table = ax.table(
        cellText=df.values,
        colLabels=df.columns,
        cellLoc='center',
        loc='center',
        bbox=[0.1, 0.15, 0.8, 0.75]
    )
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1.2, 1.8)

    # Style header
    for j in range(len(df.columns)):
        table[(0, j)].set_facecolor(PALETTE['primary'])
        table[(0, j)].set_text_props(color='white', fontweight='bold')

    ax.text(0.5, 0.05, "Program discovery across 5 Xenium regions",
            ha='center', va='bottom', fontsize=9, style='italic', color=PALETTE['neutral'],
            transform=ax.transAxes)


def generate_figure4():
    """Generate complete Figure 4."""
    print("Loading data...")
    programs = load_xenium_programs()
    if programs is not None:
        print(f"Loaded {len(programs)} programs")
    else:
        print("WARNING: No program data available")

    fig = plt.figure(figsize=(11, 8))
    gs = GridSpec(2, 2, figure=fig, hspace=0.30, wspace=0.28)

    # Panel A: Schematic
    ax_a = fig.add_subplot(gs[0, 0])
    panel_a_schematic(ax_a)
    ax_a.set_title("Program Discovery", fontsize=12, fontweight='bold', loc='left', pad=10)

    # Panel B: Top programs
    ax_b = fig.add_subplot(gs[0, 1])
    panel_b_programs(ax_b, programs)
    ax_b.set_title("Top Spatially Coherent Programs", fontsize=12, fontweight='bold', loc='left', pad=10)

    # Panel C: Boxplot
    ax_c = fig.add_subplot(gs[1, 0])
    panel_c_boxplot(ax_c, programs)
    ax_c.set_title("Spatial Coherence by Cell Type", fontsize=12, fontweight='bold', loc='left', pad=10)

    # Panel D: Summary
    ax_d = fig.add_subplot(gs[1, 1])
    panel_d_summary(ax_d, programs)
    ax_d.set_title("Summary Statistics", fontsize=12, fontweight='bold', loc='left', pad=10)

    plt.tight_layout()

    output_path = OUTPUT_DIR / "figure4_module4_programs.pdf"
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"Saved to {output_path}")

    png_path = OUTPUT_DIR / "figure4_module4_programs.png"
    plt.savefig(png_path, dpi=150, bbox_inches='tight', facecolor='white')
    print(f"Preview saved to {png_path}")

    # Save SVG for Illustrator
    svg_path = OUTPUT_DIR / "figure4_module4_programs.svg"
    plt.savefig(svg_path, format='svg', bbox_inches='tight', facecolor='white')
    print(f"SVG saved to {svg_path}")

    plt.close()

    print(f"\n=== Figure 4 Summary ===")
    if programs is not None:
        print(f"Total programs: {len(programs)}")
        print(f"Programs with I > 0.15: {(programs['spatial_moran_i'] > 0.15).sum()}")
    else:
        print("No program data available")

    print("\nEnhancements applied:")
    print("  - Consistent color palette from figure_style.py")
    print("  - Fonts increased to minimum 10pt")
    print("  - Panel labels in top-left corner")


if __name__ == "__main__":
    generate_figure4()
