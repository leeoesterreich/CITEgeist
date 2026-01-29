#!/usr/bin/env python
"""
15_proteostasis_setpoint.py

PROTEOSTASIS SETPOINT: MCF7 vs T47D Secretory/Stress States Under ESR1-D538G

Purpose:
Generate a 6-panel grant-quality figure showing that MCF7 and T47D occupy
distinct secretory/stress states under ESR1-D538G. The figure synthesizes
secretory chaperone co-regulation, UPR marker direction, ERAD gene divergence,
MDK mRNA anti-correlation, genome-wide opposite regulation, and chaperone
capacity into a cohesive visual narrative.

Datasets:
  GSE89888 - RNA-seq, MCF7/T47D WT vs D538G (primary)
  Pipeline outputs from scripts 10 and 13

Outputs:
  outputs/figures/fig15_proteostasis_setpoint.pdf
  outputs/figures/fig15_proteostasis_setpoint.png
  outputs/tables/proteostasis_module_summary.csv
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy import stats

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
import warnings
warnings.filterwarnings('ignore')

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

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

# ELISA MDK secretion fold changes (D538G/WT)
ELISA_MDK = {'MCF7': 1.83, 'T47D': 0.38}

# Style constants
COLORS = {
    'MCF7': '#4682B4',
    'T47D': '#CD5C5C',
    'secretory': '#E8890C',
    'upr': '#2CA02C',
    'erad': '#7B2D8E',
    'hsp90b1': '#DAA520',
    'same_as_secretion': '#E8890C',
    'inverse_repressor': '#7B2D8E',
}
PANEL_LABEL_PROPS = dict(fontsize=14, fontweight='bold', va='top', ha='left')

# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

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


# ---------------------------------------------------------------------------
# Panel functions
# ---------------------------------------------------------------------------

def panel_a_secretory_coregulation(ax, sig_df, tpm):
    """Panel A: Secretory module co-regulation grouped bar chart."""
    sec_genes = SIGNATURES['secretory_chaperones']
    genes_to_plot = sec_genes + ['HSP90B1']

    mcf7_log2fc = []
    t47d_log2fc = []
    labels = []

    for gene in genes_to_plot:
        row = sig_df[sig_df['gene'] == gene]
        if len(row) == 0:
            continue
        row = row.iloc[0]
        mcf7_log2fc.append(row['MCF7_log2FC'])
        t47d_log2fc.append(row['T47D_log2FC'])
        labels.append(gene)

    if not labels:
        ax.text(0.5, 0.5, 'No data', ha='center', va='center', transform=ax.transAxes)
        return

    x = np.arange(len(labels))
    width = 0.35

    bars_mcf7 = ax.bar(x - width / 2, mcf7_log2fc, width, label='MCF7',
                        color=COLORS['MCF7'], edgecolor='black', linewidth=0.5)
    bars_t47d = ax.bar(x + width / 2, t47d_log2fc, width, label='T47D',
                        color=COLORS['T47D'], edgecolor='black', linewidth=0.5)

    # Highlight HSP90B1 bars with gold
    for i, label in enumerate(labels):
        if label == 'HSP90B1':
            bars_mcf7[i].set_facecolor(COLORS['hsp90b1'])
            bars_t47d[i].set_facecolor(COLORS['hsp90b1'])

    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha='right', fontsize=8)
    ax.axhline(0, color='black', linewidth=0.5)
    ax.set_ylabel('log2(D538G / WT)')
    ax.legend(fontsize=8, loc='upper left')

    # Concordance annotation
    n_concordant = sum(1 for m, t in zip(mcf7_log2fc, t47d_log2fc)
                       if (m > 0) == (t > 0))
    n_total = len(labels)
    concordance_pct = int(round(100 * n_concordant / n_total)) if n_total > 0 else 0

    # Binomial test: concordance vs chance (50%)
    binom_p = stats.binomtest(n_concordant, n_total, 0.5, alternative='greater').pvalue

    ax.annotate(f"{concordance_pct}% concordant, p={binom_p:.4f}",
                xy=(0.98, 0.98), xycoords='axes fraction',
                ha='right', va='top', fontsize=8,
                bbox=dict(boxstyle='round,pad=0.3', facecolor='lightyellow', alpha=0.8))

    ax.text(-0.12, 1.08, 'A', transform=ax.transAxes, **PANEL_LABEL_PROPS)
    ax.set_title('Secretory Module Co-regulation', fontsize=10, fontweight='bold')


def panel_b_upr_direction(ax, sig_df):
    """Panel B: UPR marker direction bar chart (MCF7-D538G only)."""
    upr_genes = SIGNATURES['upr_genes']
    genes = []
    log2fcs = []
    bar_colors = []

    for gene in upr_genes:
        row = sig_df[sig_df['gene'] == gene]
        if len(row) == 0:
            continue
        row = row.iloc[0]
        fc = row['MCF7_log2FC']
        genes.append(gene)
        log2fcs.append(fc)
        if fc > 0.15:
            bar_colors.append('red')
        elif fc < -0.15:
            bar_colors.append('blue')
        else:
            bar_colors.append('gray')

    if not genes:
        ax.text(0.5, 0.5, 'No data', ha='center', va='center', transform=ax.transAxes)
        return

    x = np.arange(len(genes))
    ax.bar(x, log2fcs, color=bar_colors, edgecolor='black', linewidth=0.5)
    ax.set_xticks(x)
    ax.set_xticklabels(genes, rotation=45, ha='right', fontsize=8)
    ax.axhline(0, color='black', linewidth=0.5)
    ax.set_ylabel('log2(D538G / WT)')

    n_up = sum(1 for fc in log2fcs if fc > 0.15)
    n_down = sum(1 for fc in log2fcs if fc < -0.15)

    # One-sample t-test: mean log2FC vs 0
    t_stat, t_pval = stats.ttest_1samp(log2fcs, 0)

    ax.annotate(f"{n_up} UP, {n_down} DOWN; p={t_pval:.2f}",
                xy=(0.98, 0.98), xycoords='axes fraction',
                ha='right', va='top', fontsize=8,
                bbox=dict(boxstyle='round,pad=0.3', facecolor='lightyellow', alpha=0.8))

    ax.text(-0.12, 1.08, 'B', transform=ax.transAxes, **PANEL_LABEL_PROPS)
    ax.set_title('UPR Markers in MCF7-D538G', fontsize=10, fontweight='bold')


def panel_c_erad_heatmap(ax, sig_df):
    """Panel C: ERAD gene heatmap (8 genes x 2 cell lines)."""
    erad_genes = SIGNATURES['erad_genes']
    genes_present = []
    mcf7_log2fc = []
    t47d_log2fc = []
    mcf7_fc = []
    t47d_fc = []

    for gene in erad_genes:
        row = sig_df[sig_df['gene'] == gene]
        if len(row) == 0:
            continue
        row = row.iloc[0]
        genes_present.append(gene)
        mcf7_log2fc.append(row['MCF7_log2FC'])
        t47d_log2fc.append(row['T47D_log2FC'])
        mcf7_fc.append(row['MCF7_FC'])
        t47d_fc.append(row['T47D_FC'])

    if not genes_present:
        ax.text(0.5, 0.5, 'No data', ha='center', va='center', transform=ax.transAxes)
        return

    heatmap_data = pd.DataFrame({
        'MCF7 log2FC': mcf7_log2fc,
        'T47D log2FC': t47d_log2fc,
    }, index=genes_present)

    # Annotation: show FC (not log2FC)
    annot_data = pd.DataFrame({
        'MCF7 log2FC': [f"{fc:.2f}" for fc in mcf7_fc],
        'T47D log2FC': [f"{fc:.2f}" for fc in t47d_fc],
    }, index=genes_present)

    sns.heatmap(heatmap_data, annot=annot_data, fmt='', cmap='RdBu_r', center=0,
                ax=ax, cbar_kws={'label': 'log2FC', 'shrink': 0.8},
                linewidths=0.5, linecolor='white')

    ax.text(-0.12, 1.08, 'C', transform=ax.transAxes, **PANEL_LABEL_PROPS)
    ax.set_title('ERAD Gene Regulation', fontsize=10, fontweight='bold')


def panel_d_mdk_anticorrelation(ax, tpm):
    """Panel D: MDK mRNA anti-correlation with SEM error bars."""
    if 'MDK' not in tpm.index:
        ax.text(0.5, 0.5, 'MDK not found', ha='center', va='center', transform=ax.transAxes)
        return

    # Get replicate-level data
    mcf7_wt = tpm.loc['MDK', GROUPS['MCF7_WT']].values.astype(float)
    mcf7_d538g = tpm.loc['MDK', GROUPS['MCF7_D538G']].values.astype(float)
    t47d_wt = tpm.loc['MDK', GROUPS['T47D_WT']].values.astype(float)
    t47d_d538g = tpm.loc['MDK', GROUPS['T47D_D538G']].values.astype(float)

    # Replicate-level log2 ratios (each D538G replicate vs mean WT)
    mcf7_wt_mean = np.mean(mcf7_wt)
    t47d_wt_mean = np.mean(t47d_wt)

    mcf7_log2_ratios = np.log2(mcf7_d538g / mcf7_wt_mean) if mcf7_wt_mean > 0 else np.zeros(4)
    t47d_log2_ratios = np.log2(t47d_d538g / t47d_wt_mean) if t47d_wt_mean > 0 else np.zeros(4)

    mcf7_mean = np.mean(mcf7_log2_ratios)
    mcf7_sem = stats.sem(mcf7_log2_ratios)
    t47d_mean = np.mean(t47d_log2_ratios)
    t47d_sem = stats.sem(t47d_log2_ratios)

    x = np.arange(2)
    means = [mcf7_mean, t47d_mean]
    sems = [mcf7_sem, t47d_sem]
    colors = [COLORS['MCF7'], COLORS['T47D']]

    ax.bar(x, means, yerr=sems, color=colors, edgecolor='black', linewidth=0.5,
           capsize=5, error_kw={'linewidth': 1.5})
    ax.set_xticks(x)
    ax.set_xticklabels(['MCF7', 'T47D'], fontsize=10)
    ax.axhline(0, color='black', linewidth=0.5)
    ax.set_ylabel('MDK log2(D538G / WT)')

    # P-values from t-tests
    _, mcf7_pval = stats.ttest_ind(mcf7_d538g, mcf7_wt, equal_var=False)
    _, t47d_pval = stats.ttest_ind(t47d_d538g, t47d_wt, equal_var=False)

    # Secretion annotations
    y_offset_mcf7 = means[0] + sems[0] + 0.05 if means[0] > 0 else means[0] - sems[0] - 0.15
    y_offset_t47d = means[1] - sems[1] - 0.15 if means[1] < 0 else means[1] + sems[1] + 0.05

    ax.annotate(f"Secretion +83%\np={mcf7_pval:.4f}",
                xy=(0, means[0]), xytext=(0, y_offset_mcf7 + 0.1),
                ha='center', fontsize=7,
                bbox=dict(boxstyle='round,pad=0.2', facecolor='lightyellow', alpha=0.8))
    ax.annotate(f"Secretion -62%\np={t47d_pval:.3f}",
                xy=(1, means[1]), xytext=(1, y_offset_t47d - 0.1),
                ha='center', fontsize=7,
                bbox=dict(boxstyle='round,pad=0.2', facecolor='lightyellow', alpha=0.8))

    ax.text(-0.12, 1.08, 'D', transform=ax.transAxes, **PANEL_LABEL_PROPS)
    ax.set_title('MDK mRNA Anti-correlation', fontsize=10, fontweight='bold')


def panel_e_genome_wide_scatter(ax, ranking_df, sig_df):
    """Panel E: Genome-wide opposite regulation scatter (632 genes)."""
    if ranking_df is None or len(ranking_df) == 0:
        ax.text(0.5, 0.5, 'Ranking data not found', ha='center', va='center',
                transform=ax.transAxes)
        return

    # Classify genes by module
    sec_set = set(SIGNATURES['secretory_chaperones'])
    upr_set = set(SIGNATURES['upr_genes'])
    erad_set = set(SIGNATURES['erad_genes'])

    module_colors = []
    for _, row in ranking_df.iterrows():
        gene = row['gene']
        if gene in sec_set:
            module_colors.append(COLORS['secretory'])
        elif gene in upr_set:
            module_colors.append(COLORS['upr'])
        elif gene in erad_set:
            module_colors.append(COLORS['erad'])
        else:
            module_colors.append('lightgray')

    ax.scatter(ranking_df['MCF7_log2FC'], ranking_df['T47D_log2FC'],
               c=module_colors, s=15, alpha=0.6, edgecolors='none')

    # Highlight HSP90B1 with gold star
    hsp_row = ranking_df[ranking_df['gene'] == 'HSP90B1']
    if len(hsp_row) > 0:
        hsp = hsp_row.iloc[0]
        ax.scatter(hsp['MCF7_log2FC'], hsp['T47D_log2FC'],
                   marker='*', s=200, c=COLORS['hsp90b1'], edgecolors='black',
                   linewidths=0.5, zorder=5)
        ax.annotate('HSP90B1', (hsp['MCF7_log2FC'], hsp['T47D_log2FC']),
                    fontsize=7, fontweight='bold', xytext=(5, 5),
                    textcoords='offset points')

    # Highlight MDK with black star
    mdk_row = ranking_df[ranking_df['gene'] == 'MDK']
    if len(mdk_row) > 0:
        mdk = mdk_row.iloc[0]
        ax.scatter(mdk['MCF7_log2FC'], mdk['T47D_log2FC'],
                   marker='*', s=200, c='black', edgecolors='black',
                   linewidths=0.5, zorder=5)
        ax.annotate('MDK', (mdk['MCF7_log2FC'], mdk['T47D_log2FC']),
                    fontsize=7, fontweight='bold', xytext=(5, 5),
                    textcoords='offset points')

    ax.axhline(0, color='gray', linestyle='--', alpha=0.5)
    ax.axvline(0, color='gray', linestyle='--', alpha=0.5)
    ax.set_xlabel('MCF7 log2FC (D538G/WT)')
    ax.set_ylabel('T47D log2FC (D538G/WT)')

    # Legend
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor='lightgray',
               markersize=6, label='Other'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor=COLORS['secretory'],
               markersize=6, label='Secretory'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor=COLORS['upr'],
               markersize=6, label='UPR'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor=COLORS['erad'],
               markersize=6, label='ERAD'),
        Line2D([0], [0], marker='*', color='w', markerfacecolor=COLORS['hsp90b1'],
               markersize=10, label='HSP90B1'),
        Line2D([0], [0], marker='*', color='w', markerfacecolor='black',
               markersize=10, label='MDK'),
    ]
    ax.legend(handles=legend_elements, fontsize=6, loc='lower left', framealpha=0.8)

    ax.text(-0.12, 1.08, 'E', transform=ax.transAxes, **PANEL_LABEL_PROPS)
    ax.set_title(f'Genome-wide Opposite Regulation\n({len(ranking_df)} genes)',
                 fontsize=10, fontweight='bold')


def panel_f_chaperone_capacity(ax, tpm):
    """Panel F: Chaperone capacity bar chart (4 conditions)."""
    sec_genes = SIGNATURES['secretory_chaperones']
    present_genes = [g for g in sec_genes if g in tpm.index]

    if not present_genes:
        ax.text(0.5, 0.5, 'No data', ha='center', va='center', transform=ax.transAxes)
        return

    conditions = ['MCF7_WT', 'MCF7_D538G', 'T47D_WT', 'T47D_D538G']
    means = []
    sems = []

    for cond in conditions:
        gene_means = []
        for gene in present_genes:
            vals = tpm.loc[gene, GROUPS[cond]].values.astype(float)
            gene_means.append(np.mean(vals))
        means.append(np.mean(gene_means))
        sems.append(stats.sem(gene_means))

    x = np.arange(4)
    bar_colors = [COLORS['MCF7'], COLORS['MCF7'], COLORS['T47D'], COLORS['T47D']]
    alphas = [0.5, 1.0, 0.5, 1.0]

    for i in range(4):
        ax.bar(x[i], means[i], yerr=sems[i], color=bar_colors[i], alpha=alphas[i],
               edgecolor='black', linewidth=0.5, capsize=5, error_kw={'linewidth': 1.5})

    ax.set_xticks(x)
    ax.set_xticklabels(['MCF7\nWT', 'MCF7\nD538G', 'T47D\nWT', 'T47D\nD538G'], fontsize=8)
    ax.set_ylabel('Mean TPM (7 secretory chaperones)')

    ax.text(-0.12, 1.08, 'F', transform=ax.transAxes, **PANEL_LABEL_PROPS)
    ax.set_title('Chaperone Capacity', fontsize=10, fontweight='bold')


# ---------------------------------------------------------------------------
# Module summary table
# ---------------------------------------------------------------------------

def compute_module_summary(sig_df):
    """Compute per-module summary statistics."""
    modules = {
        'Secretory': SIGNATURES['secretory_chaperones'] + ['HSP90B1'],
        'UPR': SIGNATURES['upr_genes'],
        'ERAD': SIGNATURES['erad_genes'],
    }

    rows = []
    for module_name, gene_list in modules.items():
        present = sig_df[sig_df['gene'].isin(gene_list)]
        if len(present) == 0:
            continue
        rows.append({
            'module': module_name,
            'n_genes': len(present),
            'MCF7_mean_log2FC': present['MCF7_log2FC'].mean(),
            'MCF7_sem_log2FC': stats.sem(present['MCF7_log2FC'].values),
            'T47D_mean_log2FC': present['T47D_log2FC'].mean(),
            'T47D_sem_log2FC': stats.sem(present['T47D_log2FC'].values),
        })

    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("=" * 80)
    print("SCRIPT 15: PROTEOSTASIS SETPOINT FIGURE")
    print("  MCF7 and T47D occupy distinct secretory/stress states under ESR1-D538G")
    print("=" * 80)

    # Ensure output directories exist
    (OUTPUT_DIR / "figures").mkdir(parents=True, exist_ok=True)
    (OUTPUT_DIR / "tables").mkdir(parents=True, exist_ok=True)

    # Load data
    print("\nLoading GSE89888 RNA-seq data...")
    tpm = load_tpm()
    all_samples = [s for g in GROUPS.values() for s in g if s in tpm.columns]
    tpm = tpm[tpm[all_samples].max(axis=1) > 0.1]
    print(f"  Loaded {len(tpm)} genes")

    # Load mechanism gene signatures from script 13 output
    sig_path = OUTPUT_DIR / "tables" / "mechanism_gene_signatures.csv"
    print(f"\nLoading mechanism signatures: {sig_path}")
    sig_df = pd.read_csv(sig_path)
    print(f"  Loaded {len(sig_df)} signature genes")

    # Load unsupervised ranking from script 10 output
    ranking_path = OUTPUT_DIR / "tables" / "unsupervised_hsp90b1_ranking.csv"
    ranking_df = None
    if ranking_path.exists():
        print(f"Loading ranking data: {ranking_path}")
        ranking_df = pd.read_csv(ranking_path)
        print(f"  Loaded {len(ranking_df)} ranked genes")
    else:
        print(f"Warning: Ranking file not found: {ranking_path}")

    # -----------------------------------------------------------------------
    # Generate figure
    # -----------------------------------------------------------------------
    print("\nGenerating fig15: Proteostasis Setpoint...")

    fig = plt.figure(figsize=(20, 12))
    gs = gridspec.GridSpec(2, 3, hspace=0.40, wspace=0.35)

    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[0, 2])
    ax_d = fig.add_subplot(gs[1, 0])
    ax_e = fig.add_subplot(gs[1, 1])
    ax_f = fig.add_subplot(gs[1, 2])

    panel_a_secretory_coregulation(ax_a, sig_df, tpm)
    panel_b_upr_direction(ax_b, sig_df)
    panel_c_erad_heatmap(ax_c, sig_df)
    panel_d_mdk_anticorrelation(ax_d, tpm)
    panel_e_genome_wide_scatter(ax_e, ranking_df, sig_df)
    panel_f_chaperone_capacity(ax_f, tpm)

    fig.suptitle('Proteostasis Setpoint: MCF7 vs T47D Under ESR1-D538G',
                 fontsize=14, fontweight='bold', y=0.99)

    fig_pdf = OUTPUT_DIR / "figures" / "fig15_proteostasis_setpoint.pdf"
    fig_png = OUTPUT_DIR / "figures" / "fig15_proteostasis_setpoint.png"
    fig.savefig(fig_pdf, bbox_inches='tight')
    fig.savefig(fig_png, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {fig_pdf}")
    print(f"  Saved: {fig_png}")

    # -----------------------------------------------------------------------
    # Module summary table
    # -----------------------------------------------------------------------
    print("\nComputing module summary...")
    summary = compute_module_summary(sig_df)
    summary_path = OUTPUT_DIR / "tables" / "proteostasis_module_summary.csv"
    summary.to_csv(summary_path, index=False)
    print(f"  Saved: {summary_path}")
    print(summary.to_string(index=False))

    print("\n" + "=" * 80)
    print("SCRIPT 15 COMPLETE")
    print("=" * 80)

    return 0


if __name__ == "__main__":
    sys.exit(main())
