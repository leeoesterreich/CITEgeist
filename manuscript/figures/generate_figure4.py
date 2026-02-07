#!/usr/bin/env python3
"""
Figure 4: Module 4 Spatial Programs (Tour de Force #1)

Panels:
  A: NMF program discovery schematic (placeholder - to be created in vector editor)
  B: Example programs per cell type with top genes
  C: Moran's I validation (programs are spatially coherent)
  D: Xenium single-cell programs (resolution-agnostic proof)
  E: Bivariate relationships (co-localized vs exclusive programs)

Data sources:
  - Xenium singlecell: Benchmarking/xenium_benchmarking/CITEgeist/output_module4_validation/singlecell/
  - Patient Module 5: output/module5_integration/
"""

import os
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
import seaborn as sns
from pathlib import Path

# Paths
PROJECT_ROOT = Path(__file__).parent.parent.parent
XENIUM_MODULE4_DIR = PROJECT_ROOT / "Benchmarking/xenium_benchmarking/CITEgeist/output_module4_validation/singlecell"
PATIENT_MODULE5_DIR = PROJECT_ROOT / "output/module5_integration"
OUTPUT_DIR = Path(__file__).parent / "output"

# Ensure output directory exists
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Style settings
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.size'] = 8
plt.rcParams['axes.linewidth'] = 0.5
plt.rcParams['xtick.major.width'] = 0.5
plt.rcParams['ytick.major.width'] = 0.5

# Color palette for cell types
CELL_TYPE_COLORS = {
    'CD8+ T cells': '#E41A1C',
    'CD4+ T cells': '#377EB8',
    'Fibroblasts': '#4DAF4A',
    'Macrophages': '#984EA3',
    'Endothelial': '#FF7F00',
    'B cells': '#FFFF33',
    'Epithelial': '#A65628',
    # Patient data cell types
    'Cancer_Luminal': '#E41A1C',
    'Cancer_Basal': '#377EB8',
    'CD8_T_Cells': '#4DAF4A',
    'CD4_T_Cells': '#984EA3',
    'Macrophages': '#FF7F00',
    'Dendritic_Cells': '#FFFF33',
    'Fibroblasts': '#A65628',
    'Endothelial': '#F781BF',
    'B_Cells': '#999999',
    'Monocytes': '#66C2A5',
}


def load_xenium_programs():
    """Load all Xenium Module 4 program data across regions."""
    all_programs = []
    for region in range(5):
        region_dir = XENIUM_MODULE4_DIR / f"region_{region}"
        if region_dir.exists():
            df = pd.read_csv(region_dir / "module4_programs.csv")
            df['region'] = region
            all_programs.append(df)
    return pd.concat(all_programs, ignore_index=True)


def load_patient_relationships():
    """Load conserved relationships from patient Module 5."""
    return pd.read_csv(PATIENT_MODULE5_DIR / "module5_unified_conserved_relationships.csv")


def load_patient_aligned_programs():
    """Load aligned programs from patient Module 5."""
    return pd.read_csv(PATIENT_MODULE5_DIR / "module5_unified_aligned_programs.csv")


def panel_b_program_examples(ax, programs_df):
    """Panel B: Top programs with their genes as a heatmap-style display."""
    # Get top 2 programs per cell type by Moran's I
    top_programs = []
    for cell_type in programs_df['cell_type'].unique():
        ct_progs = programs_df[programs_df['cell_type'] == cell_type]
        top2 = ct_progs.nlargest(2, 'spatial_moran_i')
        top_programs.append(top2)

    top_df = pd.concat(top_programs)
    top_df = top_df.sort_values('spatial_moran_i', ascending=False).head(10)

    # Create text-based visualization
    y_pos = np.arange(len(top_df))
    colors = [CELL_TYPE_COLORS.get(ct, '#808080') for ct in top_df['cell_type']]

    ax.barh(y_pos, top_df['spatial_moran_i'], color=colors, alpha=0.7, height=0.6)

    # Add gene labels
    for i, (_, row) in enumerate(top_df.iterrows()):
        genes = eval(row['top_genes'])[:5]
        gene_str = ', '.join(genes)
        ax.text(row['spatial_moran_i'] + 0.02, i, gene_str, va='center', fontsize=6)

    ax.set_yticks(y_pos)
    ax.set_yticklabels([f"{row['cell_type']}\nR{row['region']}-P{row['program_id']}"
                        for _, row in top_df.iterrows()], fontsize=7)
    ax.set_xlabel("Moran's I", fontsize=8)
    ax.set_title("B. Top Spatially Coherent Programs", fontsize=10, fontweight='bold', loc='left')
    ax.set_xlim(0, 0.7)
    ax.invert_yaxis()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)


def panel_c_morans_i_validation(ax, programs_df):
    """Panel C: Moran's I distribution by cell type."""
    cell_types = programs_df['cell_type'].unique()

    # Box plot
    data_by_ct = [programs_df[programs_df['cell_type'] == ct]['spatial_moran_i'].values
                  for ct in cell_types]

    bp = ax.boxplot(data_by_ct, patch_artist=True, widths=0.6)

    # Color boxes
    for patch, ct in zip(bp['boxes'], cell_types):
        patch.set_facecolor(CELL_TYPE_COLORS.get(ct, '#808080'))
        patch.set_alpha(0.7)

    # Add threshold line
    ax.axhline(y=0.2, color='red', linestyle='--', linewidth=1, alpha=0.7, label='I=0.2 threshold')

    ax.set_xticklabels([ct.replace(' ', '\n') for ct in cell_types], fontsize=7, rotation=45, ha='right')
    ax.set_ylabel("Moran's I", fontsize=8)
    ax.set_title("C. Spatial Coherence by Cell Type", fontsize=10, fontweight='bold', loc='left')
    ax.legend(loc='upper right', fontsize=6)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Add summary stats
    n_above = (programs_df['spatial_moran_i'] > 0.2).sum()
    total = len(programs_df)
    ax.text(0.98, 0.02, f'{n_above}/{total} programs\nI > 0.2',
            transform=ax.transAxes, fontsize=7, ha='right', va='bottom')


def panel_d_xenium_summary(ax, programs_df):
    """Panel D: Xenium single-cell program summary statistics."""
    # Summary by region
    summary = programs_df.groupby('region').agg({
        'spatial_moran_i': ['mean', 'max', 'count'],
        'variance_explained': 'mean'
    }).round(3)
    summary.columns = ['Mean I', 'Max I', 'N Programs', 'Mean Var%']

    # Plot as table
    ax.axis('off')
    table = ax.table(
        cellText=summary.values,
        rowLabels=[f'Region {i}' for i in summary.index],
        colLabels=summary.columns,
        cellLoc='center',
        loc='center',
        bbox=[0, 0, 1, 0.8]
    )
    table.auto_set_font_size(False)
    table.set_fontsize(8)
    table.scale(1.2, 1.5)

    ax.set_title("D. Xenium Single-Cell Resolution\n(5 regions, 7 cell types, 175 programs)",
                 fontsize=10, fontweight='bold', loc='left')

    # Add note about resolution-agnostic
    ax.text(0.5, -0.1, "Same Module 4 workflow at single-cell resolution",
            transform=ax.transAxes, fontsize=7, ha='center', style='italic')


def panel_e_bivariate_relationships(ax, relationships_df):
    """Panel E: Bivariate relationships heatmap."""
    # Create pivot for relationship types
    rel_counts = relationships_df['relationship_type'].value_counts()

    # Pie chart of relationship types
    colors = {'co-localized': '#2ecc71', 'independent': '#95a5a6', 'exclusive': '#e74c3c'}
    wedges, texts, autotexts = ax.pie(
        rel_counts.values,
        labels=rel_counts.index,
        autopct='%1.0f%%',
        colors=[colors.get(r, '#808080') for r in rel_counts.index],
        startangle=90,
        explode=[0.05 if r != 'independent' else 0 for r in rel_counts.index]
    )

    for autotext in autotexts:
        autotext.set_fontsize(8)

    ax.set_title("E. Conserved Program Relationships\n(191 total across 14 patients)",
                 fontsize=10, fontweight='bold', loc='left')


def panel_e_bivariate_heatmap(ax, relationships_df, aligned_programs_df):
    """Panel E alternative: Heatmap of bivariate Moran's I."""
    # Get unique programs
    all_programs = set(relationships_df['program1_id'].unique()) | set(relationships_df['program2_id'].unique())
    all_programs = sorted(all_programs)[:20]  # Top 20 for readability

    # Create matrix
    n = len(all_programs)
    matrix = np.zeros((n, n))
    prog_to_idx = {p: i for i, p in enumerate(all_programs)}

    for _, row in relationships_df.iterrows():
        p1, p2 = row['program1_id'], row['program2_id']
        if p1 in prog_to_idx and p2 in prog_to_idx:
            i, j = prog_to_idx[p1], prog_to_idx[p2]
            matrix[i, j] = row['mean_bivariate_i']
            matrix[j, i] = row['mean_bivariate_i']

    # Plot heatmap
    im = ax.imshow(matrix, cmap='RdBu_r', vmin=-0.2, vmax=0.2, aspect='auto')

    # Labels - get cell types for programs
    prog_labels = []
    for p in all_programs:
        ct_row = aligned_programs_df[aligned_programs_df['program_id'] == p]
        if len(ct_row) > 0:
            ct = ct_row.iloc[0]['cell_type'][:3]  # First 3 chars
            prog_labels.append(f"{p[-3:]}({ct})")
        else:
            prog_labels.append(p[-3:])

    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xticklabels(prog_labels, fontsize=5, rotation=90)
    ax.set_yticklabels(prog_labels, fontsize=5)

    # Colorbar
    cbar = plt.colorbar(im, ax=ax, shrink=0.8)
    cbar.set_label("Bivariate Moran's I", fontsize=7)
    cbar.ax.tick_params(labelsize=6)

    ax.set_title("E. Program Spatial Relationships", fontsize=10, fontweight='bold', loc='left')


def generate_figure4():
    """Generate complete Figure 4."""
    print("Loading data...")

    # Load data
    xenium_programs = load_xenium_programs()
    patient_relationships = load_patient_relationships()
    patient_aligned = load_patient_aligned_programs()

    print(f"Loaded {len(xenium_programs)} Xenium programs")
    print(f"Loaded {len(patient_relationships)} patient relationships")

    # Create figure
    fig = plt.figure(figsize=(12, 10))
    gs = GridSpec(3, 2, figure=fig, height_ratios=[1.2, 1, 1], hspace=0.35, wspace=0.3)

    # Panel A placeholder
    ax_a = fig.add_subplot(gs[0, 0])
    ax_a.text(0.5, 0.5, "Panel A: NMF Schematic\n(Create in vector editor)",
              ha='center', va='center', fontsize=12, style='italic',
              bbox=dict(boxstyle='round', facecolor='#f0f0f0'))
    ax_a.set_xlim(0, 1)
    ax_a.set_ylim(0, 1)
    ax_a.axis('off')
    ax_a.set_title("A. Protein-Anchored Program Discovery", fontsize=10, fontweight='bold', loc='left')

    # Panel B: Program examples
    ax_b = fig.add_subplot(gs[0, 1])
    panel_b_program_examples(ax_b, xenium_programs)

    # Panel C: Moran's I validation
    ax_c = fig.add_subplot(gs[1, 0])
    panel_c_morans_i_validation(ax_c, xenium_programs)

    # Panel D: Xenium summary
    ax_d = fig.add_subplot(gs[1, 1])
    panel_d_xenium_summary(ax_d, xenium_programs)

    # Panel E: Bivariate relationships (use heatmap version)
    ax_e = fig.add_subplot(gs[2, :])
    panel_e_bivariate_heatmap(ax_e, patient_relationships, patient_aligned)

    # Save
    output_path = OUTPUT_DIR / "figure4_module4_programs.pdf"
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"Saved to {output_path}")

    # Also save PNG for quick preview
    png_path = OUTPUT_DIR / "figure4_module4_programs.png"
    plt.savefig(png_path, dpi=150, bbox_inches='tight', facecolor='white')
    print(f"Preview saved to {png_path}")

    plt.close()

    # Print summary statistics
    print("\n=== Figure 4 Summary Statistics ===")
    print(f"Total Xenium programs: {len(xenium_programs)}")
    print(f"Programs with Moran's I > 0.2: {(xenium_programs['spatial_moran_i'] > 0.2).sum()}")
    print(f"Max Moran's I: {xenium_programs['spatial_moran_i'].max():.3f}")
    print(f"Mean Moran's I: {xenium_programs['spatial_moran_i'].mean():.3f}")
    print(f"\nPatient relationships: {len(patient_relationships)}")
    print(f"  Co-localized: {(patient_relationships['relationship_type'] == 'co-localized').sum()}")
    print(f"  Exclusive: {(patient_relationships['relationship_type'] == 'exclusive').sum()}")
    print(f"  Independent: {(patient_relationships['relationship_type'] == 'independent').sum()}")


if __name__ == "__main__":
    generate_figure4()
