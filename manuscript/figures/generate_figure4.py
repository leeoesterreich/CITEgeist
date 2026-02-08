#!/usr/bin/env python3
"""
Figure 4: Module 4 Spatial Programs

Panel A: SCHEMATIC - use output/schematics/figure4_panel_a_nmf.svg
Panels B, C, D: DATA - generated with matplotlib below
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.gridspec import GridSpec
from pathlib import Path

from figure_style import apply_style, PALETTE, get_cell_type_color

apply_style()

# Paths
PROJECT_ROOT = Path(__file__).parent.parent.parent
XENIUM_MODULE4_DIR = PROJECT_ROOT / "Benchmarking/xenium_benchmarking/CITEgeist/output_module4_validation/singlecell"
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


def panel_b_programs(ax, programs_df):
    """Panel B: Top programs as horizontal bars with genes."""
    ax.text(-0.12, 1.05, "B", fontsize=14, fontweight='bold', va='top', transform=ax.transAxes)

    if programs_df is None or len(programs_df) == 0:
        ax.text(0.5, 0.5, "Program data not available",
                ha='center', va='center', fontsize=10, style='italic', color=PALETTE['neutral'])
        ax.axis('off')
        return

    top_df = programs_df.nlargest(8, 'spatial_moran_i')
    y_pos = np.arange(len(top_df))
    colors = [get_cell_type_color(ct) for ct in top_df['cell_type']]

    bars = ax.barh(y_pos, top_df['spatial_moran_i'], color=colors, alpha=0.8, height=0.7)

    # Gene labels
    for i, (_, row) in enumerate(top_df.iterrows()):
        try:
            genes = eval(row['top_genes'])[:4]
            gene_str = ', '.join(genes)
        except:
            gene_str = ""
        ax.text(row['spatial_moran_i'] + 0.02, i, gene_str, va='center', fontsize=7, clip_on=False)

    labels = [f"{row['cell_type'].split()[0]}" for _, row in top_df.iterrows()]
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=9)
    ax.set_xlabel("Moran's I (Spatial Coherence)", fontsize=10)
    ax.set_xlim(0, 0.52)
    ax.invert_yaxis()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)


def panel_c_boxplot(ax, programs_df):
    """Panel C: Moran's I distribution by cell type."""
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

    # Abbreviated cell type labels
    abbrev_map = {
        'CD4+ T cells': 'CD4+', 'CD8+ T cells': 'CD8+', 'T cells': 'T',
        'Fibroblasts': 'Fibro', 'Macrophages': 'Mac', 'Endothelial': 'Endo',
        'Epithelial': 'Epi', 'B cells': 'B', 'NK cells': 'NK'
    }
    short_labels = [abbrev_map.get(ct, ct.split()[0][:6]) for ct in cell_types]
    ax.set_xticklabels(short_labels, rotation=45, ha='right', fontsize=9)
    ax.set_ylabel("Moran's I", fontsize=10)
    ax.legend(fontsize=9)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    n_above = (programs_df['spatial_moran_i'] > 0.15).sum()
    total = len(programs_df)
    ax.text(0.98, 0.98, f'{n_above}/{total} above\nthreshold',
            transform=ax.transAxes, ha='right', va='top', fontsize=9,
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))


def panel_d_summary(ax, programs_df):
    """Panel D: Summary statistics table."""
    ax.text(0.02, 0.98, "D", fontsize=14, fontweight='bold', va='top')
    ax.axis('off')

    if programs_df is None or len(programs_df) == 0:
        ax.text(0.5, 0.5, "Program data not available",
                ha='center', va='center', fontsize=10, style='italic', color=PALETTE['neutral'])
        return

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

    for j in range(len(df.columns)):
        table[(0, j)].set_facecolor(PALETTE['primary'])
        table[(0, j)].set_text_props(color='white', fontweight='bold')

    ax.text(0.5, 0.05, "Program discovery across 5 Xenium regions",
            ha='center', va='bottom', fontsize=9, style='italic', color=PALETTE['neutral'],
            transform=ax.transAxes)


def generate_figure4():
    """Generate Figure 4 by assembling schematic and data panels."""
    print("Loading data...")
    programs = load_xenium_programs()
    if programs is not None:
        print(f"Loaded {len(programs)} programs")
    else:
        print("WARNING: No program data available")

    # Load schematic
    panel_a_img = load_schematic("figure4_panel_a_nmf.png")

    fig = plt.figure(figsize=(11, 8))
    gs = GridSpec(2, 2, figure=fig, hspace=0.20, wspace=0.22)

    # Panel A: NMF Program Discovery schematic
    ax_a = fig.add_subplot(gs[0, 0])
    if panel_a_img is not None:
        ax_a.imshow(panel_a_img)
    ax_a.axis('off')
    ax_a.text(-0.02, 1.05, "A", fontsize=16, fontweight='bold', va='top', transform=ax_a.transAxes)

    # Panel B: Top programs (DATA)
    ax_b = fig.add_subplot(gs[0, 1])
    panel_b_programs(ax_b, programs)
    ax_b.set_title("Top Spatially Coherent Programs", fontsize=12, fontweight='bold', loc='left', pad=10)

    # Panel C: Boxplot (DATA)
    ax_c = fig.add_subplot(gs[1, 0])
    panel_c_boxplot(ax_c, programs)
    ax_c.set_title("Spatial Coherence by Cell Type", fontsize=12, fontweight='bold', loc='left', pad=10)

    # Panel D: Summary (DATA)
    ax_d = fig.add_subplot(gs[1, 1])
    panel_d_summary(ax_d, programs)
    ax_d.set_title("Summary Statistics", fontsize=12, fontweight='bold', loc='left', pad=10)

    plt.tight_layout()

    for fmt, dpi in [('pdf', 300), ('png', 150), ('svg', None)]:
        output_path = OUTPUT_DIR / f"figure4_module4_programs.{fmt}"
        if fmt == 'svg':
            plt.savefig(output_path, format='svg', bbox_inches='tight', facecolor='white')
        else:
            plt.savefig(output_path, dpi=dpi, bbox_inches='tight', facecolor='white')
        print(f"Saved: {output_path}")

    plt.close()

    print("\n" + "=" * 60)
    print("Figure 4: Module 4 Programs")
    print("=" * 60)
    print("\nPanel A: SCHEMATIC - use output/schematics/figure4_panel_a_nmf.svg")
    print("Panels B, C, D: DATA - generated above")


if __name__ == "__main__":
    generate_figure4()
