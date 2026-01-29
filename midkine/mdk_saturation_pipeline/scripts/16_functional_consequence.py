#!/usr/bin/env python
"""
16_functional_consequence.py

FUNCTIONAL CONSEQUENCE: Secretory Capacity, UPR Tone, and Stress Response
Diverge to Determine Protein Fate

Purpose:
Generates a 4-panel grant-quality figure showing how secretory capacity,
UPR activation, and ERAD stress response diverge between MCF7 and T47D
cell lines carrying D538G ESR1 mutations, and how these collectively
determine MDK protein fate (secretion vs degradation).

Panels:
  A - Module-level response (grouped bar chart: Secretory, UPR, ERAD)
  B - Three-module heatmap (ERAD, UPR, Secretory genes grouped)
  C - MDK mRNA vs secretory capacity scatter (4 conditions)
  D - Integrated evidence table

Datasets:
  GSE89888  - RNA-seq, MCF7/T47D WT vs D538G
  outputs/tables/mechanism_gene_signatures.csv (from script 13)

Outputs:
  outputs/figures/fig16_functional_consequence.pdf
  outputs/figures/fig16_functional_consequence.png
"""

import sys
from pathlib import Path
import pandas as pd
import numpy as np
import yaml
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import matplotlib.gridspec as gridspec
import seaborn as sns
import warnings
warnings.filterwarnings('ignore')

SCRIPT_DIR = Path(__file__).parent
PIPELINE_DIR = SCRIPT_DIR.parent
CONFIG_PATH = PIPELINE_DIR / "config.yaml"

with open(CONFIG_PATH) as f:
    config = yaml.safe_load(f)

DATA_DIR = PIPELINE_DIR / config['paths']['data_dir']
OUTPUT_DIR = PIPELINE_DIR / config['paths']['output_dir']
SIGNATURES = config['mechanism_discrimination']

# Sample groups for GSE89888
GROUPS = {
    'MCF7_WT': ['GSM2392606', 'GSM2392607', 'GSM2392608', 'GSM2392609'],
    'MCF7_D538G': ['GSM2392614', 'GSM2392615', 'GSM2392616', 'GSM2392617'],
    'T47D_WT': ['GSM2392582', 'GSM2392583', 'GSM2392584', 'GSM2392585'],
    'T47D_D538G': ['GSM2392590', 'GSM2392591', 'GSM2392592', 'GSM2392593'],
}

# Style constants
COLORS = {
    'MCF7': '#4682B4',
    'T47D': '#CD5C5C',
    'secretory': '#E8890C',
    'upr': '#2CA02C',
    'erad': '#7B2D8E',
    'hsp90b1': '#DAA520',
}
PANEL_LABEL_PROPS = dict(fontsize=14, fontweight='bold', va='top', ha='left')


def load_tpm():
    """Load GSE89888 TPM data with gene symbol mapping."""
    filepath = DATA_DIR / config['datasets']['GSE89888']
    df = pd.read_csv(filepath, sep='\t', compression='gzip', index_col=0)
    try:
        import mygene
        mg = mygene.MyGeneInfo()
        results = mg.querymany(df.index.astype(str).tolist(), scopes='entrezgene',
                              fields='symbol', species='human', returnall=True)
        id_to_sym = {str(h['query']): h['symbol'] for h in results['out'] if 'symbol' in h}
        df.index = df.index.astype(str).map(lambda x: id_to_sym.get(x, x))
    except Exception as e:
        print(f"Warning: Gene mapping failed: {e}")
    return df


def compute_module_log2fc(tpm, gene_list):
    """
    Compute per-gene log2FC (D538G/WT) for MCF7 and T47D, return arrays.

    For each gene: log2(mean D538G TPM / mean WT TPM).
    Returns dict with 'MCF7' and 'T47D' arrays of log2FC values.
    """
    mcf7_fcs = []
    t47d_fcs = []
    for gene in gene_list:
        if gene not in tpm.index:
            continue
        mcf7_wt = tpm.loc[gene, GROUPS['MCF7_WT']].values.astype(float)
        mcf7_d538g = tpm.loc[gene, GROUPS['MCF7_D538G']].values.astype(float)
        t47d_wt = tpm.loc[gene, GROUPS['T47D_WT']].values.astype(float)
        t47d_d538g = tpm.loc[gene, GROUPS['T47D_D538G']].values.astype(float)

        mcf7_wt_mean = np.mean(mcf7_wt)
        t47d_wt_mean = np.mean(t47d_wt)

        if mcf7_wt_mean > 0:
            mcf7_fcs.append(np.log2(np.mean(mcf7_d538g) / mcf7_wt_mean))
        if t47d_wt_mean > 0:
            t47d_fcs.append(np.log2(np.mean(t47d_d538g) / t47d_wt_mean))

    return {'MCF7': np.array(mcf7_fcs), 'T47D': np.array(t47d_fcs)}


def compute_condition_tpm(tpm, gene_list, condition_samples):
    """Compute mean TPM of gene_list genes for a given condition (sample list)."""
    gene_means = []
    for gene in gene_list:
        if gene not in tpm.index:
            continue
        vals = tpm.loc[gene, condition_samples].values.astype(float)
        gene_means.append(np.mean(vals))
    if len(gene_means) == 0:
        return 0.0
    return np.mean(gene_means)


# ============================================================================
# PANEL FUNCTIONS
# ============================================================================

def plot_panel_a(ax, tpm):
    """Panel A: Module-level response grouped bar chart."""
    modules = {
        'Secretory': SIGNATURES['secretory_chaperones'],
        'UPR': SIGNATURES['upr_genes'],
        'ERAD': SIGNATURES['erad_genes'],
    }
    module_colors = {
        'Secretory': COLORS['secretory'],
        'UPR': COLORS['upr'],
        'ERAD': COLORS['erad'],
    }

    module_names = list(modules.keys())
    mcf7_means = []
    mcf7_sems = []
    t47d_means = []
    t47d_sems = []

    for mod_name in module_names:
        fcs = compute_module_log2fc(tpm, modules[mod_name])
        mcf7_vals = fcs['MCF7']
        t47d_vals = fcs['T47D']
        mcf7_means.append(np.mean(mcf7_vals) if len(mcf7_vals) > 0 else 0)
        mcf7_sems.append(stats.sem(mcf7_vals) if len(mcf7_vals) > 1 else 0)
        t47d_means.append(np.mean(t47d_vals) if len(t47d_vals) > 0 else 0)
        t47d_sems.append(stats.sem(t47d_vals) if len(t47d_vals) > 1 else 0)

    x = np.arange(len(module_names))
    width = 0.35

    ax.bar(x - width / 2, mcf7_means, width, yerr=mcf7_sems,
           label='MCF7', color=COLORS['MCF7'], edgecolor='black',
           capsize=4, error_kw={'linewidth': 1.2})
    ax.bar(x + width / 2, t47d_means, width, yerr=t47d_sems,
           label='T47D', color=COLORS['T47D'], edgecolor='black',
           capsize=4, error_kw={'linewidth': 1.2})

    ax.axhline(0, color='black', linewidth=0.5)
    ax.set_xticks(x)

    # Module labels colored and bold
    tick_labels = []
    for mod_name in module_names:
        tick_labels.append(mod_name)
    ax.set_xticklabels(tick_labels, fontsize=11)
    for i, label in enumerate(ax.get_xticklabels()):
        label.set_color(module_colors[module_names[i]])
        label.set_fontweight('bold')

    ax.set_ylabel('Mean log2FC (D538G / WT)', fontsize=10)
    ax.legend(fontsize=9, framealpha=0.9)
    ax.set_title('Module-Level Response', fontsize=12, fontweight='bold')
    ax.text(-0.12, 1.05, 'A', transform=ax.transAxes, **PANEL_LABEL_PROPS)


def plot_panel_b(ax, sig_df):
    """Panel B: Three-module heatmap (grouped, NOT clustered)."""
    module_order = {
        'ERAD': SIGNATURES['erad_genes'],
        'UPR': SIGNATURES['upr_genes'],
        'Secretory': SIGNATURES['secretory_chaperones'],
    }
    module_color_map = {
        'ERAD': COLORS['erad'],
        'UPR': COLORS['upr'],
        'Secretory': COLORS['secretory'],
    }

    # Build ordered gene list and module assignments
    ordered_genes = []
    gene_module = {}
    for mod_name, gene_list in module_order.items():
        for gene in gene_list:
            if gene in sig_df['gene'].values:
                ordered_genes.append(gene)
                gene_module[gene] = mod_name

    if len(ordered_genes) == 0:
        ax.text(0.5, 0.5, 'No signature genes found',
                ha='center', va='center', transform=ax.transAxes)
        return

    # Build heatmap data (log2FC) and annotation data (FC)
    heat_rows = []
    annot_rows = []
    for gene in ordered_genes:
        row = sig_df[sig_df['gene'] == gene].iloc[0]
        heat_rows.append({
            'MCF7\nlog2FC': row['MCF7_log2FC'],
            'T47D\nlog2FC': row['T47D_log2FC'],
        })
        annot_rows.append({
            'MCF7\nlog2FC': f"{row['MCF7_FC']:.2f}",
            'T47D\nlog2FC': f"{row['T47D_FC']:.2f}",
        })

    heat_df = pd.DataFrame(heat_rows, index=ordered_genes)
    annot_df = pd.DataFrame(annot_rows, index=ordered_genes)

    # Find absolute max for symmetric color scale
    vmax = max(abs(heat_df.values.min()), abs(heat_df.values.max()), 0.5)

    sns.heatmap(heat_df, annot=annot_df, fmt='', cmap='RdBu_r', center=0,
                vmin=-vmax, vmax=vmax, ax=ax, linewidths=0.5, linecolor='white',
                cbar_kws={'label': 'log2FC', 'shrink': 0.7})

    # Add module color bar on left side
    # Determine row boundaries for each module
    current_row = 0
    module_boundaries = []
    for mod_name, gene_list in module_order.items():
        n_genes = sum(1 for g in gene_list if g in ordered_genes)
        if n_genes > 0:
            module_boundaries.append((mod_name, current_row, current_row + n_genes))
            current_row += n_genes

    # Draw colored rectangles as module indicators
    for mod_name, start_row, end_row in module_boundaries:
        rect = plt.Rectangle((-0.6, start_row), 0.4, end_row - start_row,
                             color=module_color_map[mod_name], clip_on=False,
                             transform=ax.get_yaxis_transform())
        ax.add_patch(rect)
        # Module label
        mid_y = (start_row + end_row) / 2
        ax.text(-0.4, mid_y, mod_name, ha='center', va='center', fontsize=7,
                fontweight='bold', color='white', rotation=90,
                transform=ax.get_yaxis_transform())

    # Add horizontal white lines between modules
    for mod_name, start_row, end_row in module_boundaries[:-1]:
        ax.axhline(y=end_row, color='white', linewidth=3)

    ax.set_title('Gene-Level Fold Changes by Module', fontsize=12, fontweight='bold')
    ax.text(-0.25, 1.05, 'B', transform=ax.transAxes, **PANEL_LABEL_PROPS)


def plot_panel_c(ax, tpm):
    """Panel C: MDK mRNA vs secretory capacity scatter (4 conditions, absolute TPM)."""
    sec_genes = SIGNATURES['secretory_chaperones']

    conditions = {
        'MCF7_WT': ('MCF7 WT', COLORS['MCF7'], 'o'),
        'MCF7_D538G': ('MCF7 D538G', COLORS['MCF7'], 's'),
        'T47D_WT': ('T47D WT', COLORS['T47D'], 'o'),
        'T47D_D538G': ('T47D D538G', COLORS['T47D'], 's'),
    }

    x_vals = []
    y_vals = []
    labels = []
    plot_colors = []
    markers = []

    for cond_key, (label, color, marker) in conditions.items():
        samples = GROUPS[cond_key]
        sec_mean = compute_condition_tpm(tpm, sec_genes, samples)

        if 'MDK' in tpm.index:
            mdk_mean = np.mean(tpm.loc['MDK', samples].values.astype(float))
        else:
            mdk_mean = 0.0

        x_vals.append(sec_mean)
        y_vals.append(mdk_mean)
        labels.append(label)
        plot_colors.append(color)
        markers.append(marker)

    # Plot each point
    for i in range(len(x_vals)):
        ax.scatter(x_vals[i], y_vals[i], c=plot_colors[i], marker=markers[i],
                   s=120, edgecolors='black', linewidth=1.2, zorder=3)
        ax.annotate(labels[i], (x_vals[i], y_vals[i]),
                    fontsize=8, textcoords='offset points', xytext=(8, 6))

    # Pearson correlation and trend line
    x_arr = np.array(x_vals)
    y_arr = np.array(y_vals)
    if len(x_arr) >= 3:
        r, p = stats.pearsonr(x_arr, y_arr)
        # Trend line
        slope, intercept = np.polyfit(x_arr, y_arr, 1)
        x_line = np.linspace(x_arr.min() * 0.9, x_arr.max() * 1.1, 50)
        y_line = slope * x_line + intercept
        ax.plot(x_line, y_line, '--', color='gray', alpha=0.7, linewidth=1.5)
        ax.text(0.05, 0.95, f'r = {r:.3f}\np = {p:.4f}',
                transform=ax.transAxes, fontsize=9, va='top',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    # Legend for markers
    legend_elements = [
        plt.Line2D([0], [0], marker='o', color='w', markeredgecolor='black',
                   markerfacecolor='gray', markersize=8, label='WT'),
        plt.Line2D([0], [0], marker='s', color='w', markeredgecolor='black',
                   markerfacecolor='gray', markersize=8, label='D538G'),
        Patch(facecolor=COLORS['MCF7'], edgecolor='black', label='MCF7'),
        Patch(facecolor=COLORS['T47D'], edgecolor='black', label='T47D'),
    ]
    ax.legend(handles=legend_elements, fontsize=8, loc='lower right', framealpha=0.9)

    ax.set_xlabel('Mean Secretory Chaperone TPM', fontsize=10)
    ax.set_ylabel('MDK TPM', fontsize=10)
    ax.set_title('MDK mRNA vs Secretory Capacity', fontsize=12, fontweight='bold')
    ax.text(-0.12, 1.05, 'C', transform=ax.transAxes, **PANEL_LABEL_PROPS)


def plot_panel_d(ax):
    """Panel D: Integrated evidence table."""
    ax.axis('off')

    col_labels = ['Context', 'Key Finding', 'Dataset', 'Statistic',
                  'Experimental\nPrediction']
    row_data = [
        ['Different\ncistrome',
         '168/632 opposite genes\nhave diff. ER binding',
         'GSE125117\n+ GSE89888',
         'p < 0.05\n(Fisher)',
         'CHX chase: T47D-D538G\nshorter MDK half-life'],
        ['Different\nproteostasis',
         '93% secretory co-regulation,\nUPR activation (3/5 UP)',
         'GSE89888',
         'p=0.0009,\np=0.04',
         'Pulse-chase: MCF7-D538G\nhigher fractional secretion'],
        ['Module\ncoordination',
         'MDK mRNA anti-correlates\nwith protein secretion',
         'GSE89888\n+ ELISA',
         'p=0.0008,\np=0.002',
         'Proteasome/lysosome\ninhibitors rescue T47D MDK'],
    ]

    # Context column colors
    context_colors = [COLORS['secretory'], COLORS['upr'], COLORS['erad']]

    table = ax.table(
        cellText=row_data,
        colLabels=col_labels,
        cellLoc='center',
        loc='center',
        colWidths=[0.15, 0.28, 0.15, 0.14, 0.28],
    )

    table.auto_set_font_size(False)
    table.set_fontsize(8)
    table.scale(1.0, 2.8)

    # Style header row
    for j in range(len(col_labels)):
        cell = table[0, j]
        cell.set_facecolor('#2C3E50')
        cell.set_text_props(color='white', fontweight='bold', fontsize=8)
        cell.set_edgecolor('white')

    # Style data rows
    for i in range(len(row_data)):
        # Context column colored
        cell_ctx = table[i + 1, 0]
        cell_ctx.set_facecolor(context_colors[i])
        cell_ctx.set_text_props(color='white', fontweight='bold')
        cell_ctx.set_edgecolor('white')

        # Other columns
        for j in range(1, len(col_labels)):
            cell = table[i + 1, j]
            cell.set_facecolor('#F8F9FA')
            cell.set_edgecolor('#DEE2E6')

    ax.set_title('Integrated Evidence for Divergent Protein Fate',
                 fontsize=12, fontweight='bold', pad=20)
    ax.text(-0.05, 1.05, 'D', transform=ax.transAxes, **PANEL_LABEL_PROPS)


# ============================================================================
# MAIN
# ============================================================================

def main():
    print("=" * 80)
    print("SCRIPT 16: FUNCTIONAL CONSEQUENCE FIGURE")
    print("  Secretory capacity, UPR tone, and stress response diverge")
    print("  to determine protein fate")
    print("=" * 80)

    # Ensure output directories exist
    (OUTPUT_DIR / "figures").mkdir(parents=True, exist_ok=True)

    # Load data
    print("\nLoading GSE89888 RNA-seq data...")
    tpm = load_tpm()
    all_samples = [s for g in GROUPS.values() for s in g if s in tpm.columns]
    tpm = tpm[tpm[all_samples].max(axis=1) > 0.1]
    print(f"Loaded {len(tpm)} genes")

    # Load mechanism gene signatures from script 13 output
    sig_path = OUTPUT_DIR / "tables" / "mechanism_gene_signatures.csv"
    print(f"\nLoading gene signatures from: {sig_path}")
    sig_df = pd.read_csv(sig_path)
    print(f"Loaded {len(sig_df)} signature genes")

    # Quick stats
    for mod_name, gene_list in [('Secretory', SIGNATURES['secretory_chaperones']),
                                 ('UPR', SIGNATURES['upr_genes']),
                                 ('ERAD', SIGNATURES['erad_genes'])]:
        found = [g for g in gene_list if g in sig_df['gene'].values]
        print(f"  {mod_name}: {len(found)}/{len(gene_list)} genes present")

    # Build figure
    print("\nGenerating fig16: Functional Consequence...")
    fig = plt.figure(figsize=(16, 12))
    gs = gridspec.GridSpec(2, 2, hspace=0.35, wspace=0.30,
                           left=0.08, right=0.95, top=0.93, bottom=0.05)

    # Panel A: Module-level response
    ax_a = fig.add_subplot(gs[0, 0])
    plot_panel_a(ax_a, tpm)

    # Panel B: Three-module heatmap
    ax_b = fig.add_subplot(gs[0, 1])
    plot_panel_b(ax_b, sig_df)

    # Panel C: MDK mRNA vs secretory capacity
    ax_c = fig.add_subplot(gs[1, 0])
    plot_panel_c(ax_c, tpm)

    # Panel D: Integrated evidence table
    ax_d = fig.add_subplot(gs[1, 1])
    plot_panel_d(ax_d)

    fig.suptitle('Functional Consequence: Divergent Proteostasis Determines MDK Fate',
                 fontsize=14, fontweight='bold', y=0.98)

    # Save
    fig_pdf = OUTPUT_DIR / "figures" / "fig16_functional_consequence.pdf"
    fig_png = OUTPUT_DIR / "figures" / "fig16_functional_consequence.png"
    fig.savefig(fig_pdf, bbox_inches='tight')
    fig.savefig(fig_png, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"\nSaved: {fig_pdf}")
    print(f"Saved: {fig_png}")

    print("\n" + "=" * 80)
    print("SCRIPT 16 COMPLETE")
    print("=" * 80)

    return 0


if __name__ == "__main__":
    sys.exit(main())
