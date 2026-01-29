# Attractor State Evidence Figures Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Create three grant-quality figures (14-16) showing computational evidence that ESR1-D538G drives distinct attractor states in MCF7 vs T47D.

**Architecture:** Three standalone Python scripts following existing pipeline conventions (config.yaml, load_tpm with mygene, matplotlib Agg backend). Each script loads existing pipeline outputs (scripts 10, 11, 13) and GSE89888 data, computes panel-specific metrics, and generates a multi-panel figure. Shared patterns (load_tpm, load_peaks, GROUPS dict, color constants) are copy-pasted per script to maintain independence.

**Tech Stack:** Python 3.10, pandas, numpy, scipy, matplotlib, seaborn, mygene, gprofiler-official, yaml

---

## Shared Constants (copy into each script)

```python
# Standard pipeline preamble
import sys
import gzip
from pathlib import Path
import pandas as pd
import numpy as np
import yaml
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
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

GROUPS = {
    'MCF7_WT': ['GSM2392606', 'GSM2392607', 'GSM2392608', 'GSM2392609'],
    'MCF7_D538G': ['GSM2392614', 'GSM2392615', 'GSM2392616', 'GSM2392617'],
    'T47D_WT': ['GSM2392582', 'GSM2392583', 'GSM2392584', 'GSM2392585'],
    'T47D_D538G': ['GSM2392590', 'GSM2392591', 'GSM2392592', 'GSM2392593'],
}

# Grant-quality color palette
COLORS = {
    'MCF7': '#4682B4',       # steelblue
    'T47D': '#CD5C5C',       # indianred
    'secretory': '#E8890C',  # orange
    'upr': '#2CA02C',        # green
    'erad': '#7B2D8E',       # purple
    'hsp90b1': '#DAA520',    # goldenrod
    'same_as_secretion': '#E8890C',
    'inverse_repressor': '#7B2D8E',
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


def load_peaks(filepath):
    """Load BED peak file."""
    peaks = []
    with gzip.open(filepath, 'rt') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                try:
                    peaks.append({
                        'chr': parts[0],
                        'start': int(parts[1]),
                        'end': int(parts[2]),
                        'score': float(parts[4]) if len(parts) > 4 else 1.0,
                    })
                except (ValueError, IndexError):
                    continue
    return pd.DataFrame(peaks)
```

---

## Task 1: Script 14 — Cistrome Divergence Figure

**Files:**
- Create: `midkine/mdk_saturation_pipeline/scripts/14_cistrome_divergence.py`
- Read: `midkine/mdk_saturation_pipeline/outputs/tables/venn_expression_binding.csv`
- Read: `midkine/mdk_saturation_pipeline/outputs/tables/unsupervised_hsp90b1_ranking.csv`
- Output: `midkine/mdk_saturation_pipeline/outputs/figures/fig14_cistrome_divergence.{pdf,png}`
- Output: `midkine/mdk_saturation_pipeline/outputs/tables/cistrome_intersection_go_enrichment.csv`

**Step 1: Create script 14 with data loading and Panel A (Venn diagram)**

```python
#!/usr/bin/env python
"""
14_cistrome_divergence.py

ESR1-D538G engages distinct regulatory programs in MCF7 versus T47D.

Panels:
  A) Venn diagram: opposite expression ∩ differential ER binding
  B) ER binding change heatmap for 168 intersection genes
  C) ATAC-seq accessibility at intersection gene promoters
  D) GO functional enrichment by correlation type

Datasets: GSE89888, GSE125117, GSE254216, script 10/11 outputs
Outputs:  outputs/figures/fig14_cistrome_divergence.pdf/png
          outputs/tables/cistrome_intersection_go_enrichment.csv
"""

# [INSERT SHARED PREAMBLE HERE]
from matplotlib_venn import venn2

def load_venn_data():
    """Load script 11 output: 632 genes with binding data."""
    path = OUTPUT_DIR / "tables" / "venn_expression_binding.csv"
    if not path.exists():
        sys.exit(f"Missing: {path}. Run script 11 first.")
    return pd.read_csv(path)


def load_atac_peaks():
    """Load ATAC-seq summit files for MCF7 and T47D."""
    mcf7_path = DATA_DIR / config['datasets']['GSE254216_MCF7']
    t47d_path = DATA_DIR / config['datasets']['GSE254216_T47D']
    mcf7 = load_peaks(mcf7_path) if mcf7_path.exists() else pd.DataFrame()
    t47d = load_peaks(t47d_path) if t47d_path.exists() else pd.DataFrame()
    return mcf7, t47d


def count_peaks_in_region(peaks_df, chrom, center, window=10000):
    """Count peaks within ±window of a position."""
    if peaks_df is None or len(peaks_df) == 0:
        return 0
    chr_peaks = peaks_df[peaks_df['chr'] == chrom]
    overlapping = chr_peaks[
        (chr_peaks['start'] < center + window) &
        (chr_peaks['end'] > center - window)
    ]
    return len(overlapping)


def run_go_enrichment(gene_list, background_list):
    """Run GO enrichment using g:Profiler."""
    from gprofiler import GProfiler
    gp = GProfiler(return_dataframe=True)
    result = gp.profile(
        organism='hsapiens',
        query=gene_list,
        background=background_list,
        sources=['GO:BP'],
        significance_threshold_method='fdr',
        user_threshold=0.05,
        no_evidences=True,
    )
    return result


def main():
    print("=" * 80)
    print("SCRIPT 14: CISTROME DIVERGENCE FIGURE")
    print("=" * 80)

    # Load data
    venn_df = load_venn_data()
    intersection = venn_df[venn_df['has_differential_binding']].copy()
    no_binding = venn_df[~venn_df['has_differential_binding']].copy()

    print(f"Total opposite-expression genes: {len(venn_df)}")
    print(f"Intersection (+ diff binding): {len(intersection)}")

    same_genes = intersection[intersection['correlation_type'] == 'same_as_secretion']['gene'].tolist()
    inverse_genes = intersection[intersection['correlation_type'] == 'inverse_repressor']['gene'].tolist()
    all_background = venn_df['gene'].tolist()

    print(f"  same_as_secretion: {len(same_genes)}")
    print(f"  inverse_repressor: {len(inverse_genes)}")

    # ATAC-seq data
    print("\nLoading ATAC-seq data...")
    atac_mcf7, atac_t47d = load_atac_peaks()
    print(f"  MCF7 ATAC peaks: {len(atac_mcf7)}")
    print(f"  T47D ATAC peaks: {len(atac_t47d)}")

    # Compute ATAC signal at each intersection gene promoter
    atac_results = []
    for _, row in intersection.iterrows():
        if pd.isna(row.get('tss', np.nan)):
            continue
        mcf7_count = count_peaks_in_region(atac_mcf7, row['chr'], int(row['tss']))
        t47d_count = count_peaks_in_region(atac_t47d, row['chr'], int(row['tss']))
        atac_results.append({
            'gene': row['gene'],
            'mcf7_atac': mcf7_count,
            't47d_atac': t47d_count,
            'correlation_type': row['correlation_type'],
        })
    atac_df = pd.DataFrame(atac_results)
    print(f"  Computed ATAC for {len(atac_df)} genes")

    # GO enrichment
    print("\nRunning GO enrichment...")
    go_same = run_go_enrichment(same_genes, all_background)
    go_inverse = run_go_enrichment(inverse_genes, all_background)

    # Save GO results
    go_combined = pd.DataFrame()
    if len(go_same) > 0:
        go_same_top = go_same.head(20).copy()
        go_same_top['group'] = 'same_as_secretion'
        go_combined = pd.concat([go_combined, go_same_top])
    if len(go_inverse) > 0:
        go_inv_top = go_inverse.head(20).copy()
        go_inv_top['group'] = 'inverse_repressor'
        go_combined = pd.concat([go_combined, go_inv_top])
    go_combined.to_csv(OUTPUT_DIR / "tables" / "cistrome_intersection_go_enrichment.csv", index=False)
    print(f"Saved GO enrichment: {len(go_combined)} terms")

    # ---- FIGURE ----
    print("\nGenerating Figure 14...")
    fig = plt.figure(figsize=(16, 12))

    # Panel A: Venn diagram
    ax_a = fig.add_subplot(2, 2, 1)
    set_all = set(venn_df['gene'])
    set_diff_bind = set(intersection['gene'])
    venn2(
        [set_all - set_diff_bind, set_diff_bind],
        set_labels=('Opposite Expression\nOnly', 'Also Differential\nER Binding'),
        set_colors=(COLORS['MCF7'], COLORS['T47D']),
        alpha=0.6,
        ax=ax_a,
    )
    ax_a.set_title('Opposite Expression ∩ Differential ER Binding\n(GSE89888 + GSE125117)', fontsize=11)
    ax_a.text(-0.05, 1.05, 'A', transform=ax_a.transAxes, **PANEL_LABEL_PROPS)

    # Panel B: Binding change heatmap
    ax_b = fig.add_subplot(2, 2, 2)
    heatmap_genes = intersection.sort_values('total_effect', ascending=False)
    heatmap_data = heatmap_genes[['gene', 'MCF7_binding_diff', 'T47D_binding_diff']].set_index('gene')
    heatmap_data.columns = ['MCF7\n(D538G − WT)', 'T47D\n(D538G − WT)']

    # Clip for visualization
    vmax = np.percentile(np.abs(heatmap_data.values), 95)
    sns.heatmap(
        heatmap_data, cmap='RdBu_r', center=0, vmin=-vmax, vmax=vmax,
        ax=ax_b, cbar_kws={'label': 'ER Binding Change', 'shrink': 0.6},
        yticklabels=False,
    )
    # Add correlation_type color bar on right
    corr_colors = [COLORS.get(ct, 'gray') for ct in heatmap_genes['correlation_type']]
    for i, c in enumerate(corr_colors):
        ax_b.add_patch(plt.Rectangle((2.05, i), 0.15, 1, color=c, transform=ax_b.transData, clip_on=False))
    ax_b.set_title('ER Binding Changes at 168 Gene Promoters\n(GSE125117)', fontsize=11)
    ax_b.set_ylabel(f'{len(intersection)} genes (ranked by expression effect)')
    ax_b.text(-0.05, 1.05, 'B', transform=ax_b.transAxes, **PANEL_LABEL_PROPS)

    # Panel C: ATAC scatter
    ax_c = fig.add_subplot(2, 2, 3)
    if len(atac_df) > 0:
        for ct, color in [('same_as_secretion', COLORS['same_as_secretion']),
                          ('inverse_repressor', COLORS['inverse_repressor'])]:
            subset = atac_df[atac_df['correlation_type'] == ct]
            ax_c.scatter(subset['mcf7_atac'], subset['t47d_atac'],
                        c=color, alpha=0.6, s=30, label=ct.replace('_', ' ').title(),
                        edgecolors='black', linewidths=0.3)
        max_val = max(atac_df['mcf7_atac'].max(), atac_df['t47d_atac'].max()) + 2
        ax_c.plot([0, max_val], [0, max_val], 'k--', alpha=0.3, linewidth=1)
        ax_c.set_xlabel('MCF7 ATAC Peaks at Promoter (±10kb)')
        ax_c.set_ylabel('T47D ATAC Peaks at Promoter (±10kb)')
        ax_c.legend(fontsize=9, loc='upper left')
    ax_c.set_title('Chromatin Accessibility at Intersection Genes\n(GSE254216)', fontsize=11)
    ax_c.text(-0.05, 1.05, 'C', transform=ax_c.transAxes, **PANEL_LABEL_PROPS)

    # Panel D: GO enrichment dot plots
    ax_d = fig.add_subplot(2, 2, 4)
    go_plot_data = []
    for group_name, go_df, color in [
        ('Same as Secretion', go_same, COLORS['same_as_secretion']),
        ('Inverse Repressor', go_inverse, COLORS['inverse_repressor']),
    ]:
        if go_df is not None and len(go_df) > 0:
            top5 = go_df.head(5)
            for _, row in top5.iterrows():
                go_plot_data.append({
                    'term': row['name'][:40] + ('...' if len(row['name']) > 40 else ''),
                    'neg_log10_p': -np.log10(row['p_value']),
                    'intersection_size': row['intersection_size'],
                    'group': group_name,
                    'color': color,
                })

    if go_plot_data:
        go_plot = pd.DataFrame(go_plot_data)
        y_pos = np.arange(len(go_plot))
        bars = ax_d.barh(y_pos, go_plot['neg_log10_p'],
                        color=go_plot['color'], alpha=0.7, edgecolor='black', linewidth=0.5)
        ax_d.set_yticks(y_pos)
        ax_d.set_yticklabels(go_plot['term'], fontsize=8)
        ax_d.set_xlabel('-log10(FDR)')
        ax_d.axvline(-np.log10(0.05), color='gray', linestyle='--', alpha=0.5, label='FDR=0.05')

        # Add group labels
        from matplotlib.patches import Patch
        legend_elements = [
            Patch(facecolor=COLORS['same_as_secretion'], label='Same as Secretion'),
            Patch(facecolor=COLORS['inverse_repressor'], label='Inverse Repressor'),
        ]
        ax_d.legend(handles=legend_elements, fontsize=8, loc='lower right')
    else:
        ax_d.text(0.5, 0.5, 'No significant GO terms\n(FDR < 0.05)', ha='center', va='center',
                 transform=ax_d.transAxes, fontsize=12)
    ax_d.set_title('GO Biological Process Enrichment\n(Intersection Genes by Type)', fontsize=11)
    ax_d.text(-0.05, 1.05, 'D', transform=ax_d.transAxes, **PANEL_LABEL_PROPS)

    plt.tight_layout()
    fig.savefig(OUTPUT_DIR / "figures" / "fig14_cistrome_divergence.pdf", bbox_inches='tight')
    fig.savefig(OUTPUT_DIR / "figures" / "fig14_cistrome_divergence.png", dpi=300, bbox_inches='tight')
    plt.close()

    print(f"Saved: outputs/figures/fig14_cistrome_divergence.pdf/png")

    print("\n" + "=" * 80)
    print("SCRIPT 14 COMPLETE")
    print("=" * 80)
    return 0


if __name__ == "__main__":
    sys.exit(main())
```

**Step 2: Verify syntax**

Run: `python -c "import ast; ast.parse(open('scripts/14_cistrome_divergence.py').read()); print('OK')"`
Expected: `OK`

**Step 3: Commit**

```bash
git add midkine/mdk_saturation_pipeline/scripts/14_cistrome_divergence.py
git commit -m "feat: add script 14 - cistrome divergence figure"
```

---

## Task 2: Script 15 — Proteostasis Setpoint Figure

**Files:**
- Create: `midkine/mdk_saturation_pipeline/scripts/15_proteostasis_setpoint.py`
- Read: `midkine/mdk_saturation_pipeline/outputs/tables/mechanism_gene_signatures.csv`
- Read: `midkine/mdk_saturation_pipeline/outputs/tables/unsupervised_hsp90b1_ranking.csv`
- Output: `midkine/mdk_saturation_pipeline/outputs/figures/fig15_proteostasis_setpoint.{pdf,png}`
- Output: `midkine/mdk_saturation_pipeline/outputs/tables/proteostasis_module_summary.csv`

**Step 1: Create script 15**

```python
#!/usr/bin/env python
"""
15_proteostasis_setpoint.py

MCF7 and T47D occupy distinct secretory/stress states under ESR1-D538G.

Panels:
  A) Secretory module co-regulation dot plot (7 chaperones + HSP90B1)
  B) UPR marker direction bar chart (MCF7-D538G)
  C) ERAD gene heatmap (both cell lines)
  D) MDK mRNA anti-correlation with secretion
  E) Genome-wide opposite regulation scatter (632 genes)
  F) Chaperone capacity bar chart (4 conditions)

Datasets: GSE89888, script 10/13 outputs
Outputs:  outputs/figures/fig15_proteostasis_setpoint.pdf/png
          outputs/tables/proteostasis_module_summary.csv
"""

# [INSERT SHARED PREAMBLE HERE]

ELISA_MDK = {'MCF7': 1.83, 'T47D': 0.38}


def load_gene_signatures():
    """Load script 13 gene signature data."""
    path = OUTPUT_DIR / "tables" / "mechanism_gene_signatures.csv"
    if not path.exists():
        sys.exit(f"Missing: {path}. Run script 13 first.")
    return pd.read_csv(path)


def load_ranking_data():
    """Load script 10 genome-wide opposite-regulation ranking."""
    path = OUTPUT_DIR / "tables" / "unsupervised_hsp90b1_ranking.csv"
    if not path.exists():
        sys.exit(f"Missing: {path}. Run script 10 first.")
    return pd.read_csv(path)


def main():
    print("=" * 80)
    print("SCRIPT 15: PROTEOSTASIS SETPOINT FIGURE")
    print("=" * 80)

    # Load data
    sig_df = load_gene_signatures()
    ranking_df = load_ranking_data()
    tpm = load_tpm()

    print(f"Signature genes: {len(sig_df)}")
    print(f"Genome-wide opposite genes: {len(ranking_df)}")

    # Classify signature genes by module
    module_map = {}
    for gene in SIGNATURES['secretory_chaperones']:
        module_map[gene] = 'Secretory'
    for gene in SIGNATURES['upr_genes']:
        module_map[gene] = 'UPR'
    for gene in SIGNATURES['erad_genes']:
        module_map[gene] = 'ERAD'

    sig_df['module'] = sig_df['gene'].map(module_map).fillna('Other')

    # Compute module summary
    module_summary = []
    for module_name, genes in [('Secretory', SIGNATURES['secretory_chaperones']),
                                ('UPR', SIGNATURES['upr_genes']),
                                ('ERAD', SIGNATURES['erad_genes'])]:
        mod_data = sig_df[sig_df['gene'].isin(genes)]
        if len(mod_data) == 0:
            continue
        module_summary.append({
            'module': module_name,
            'n_genes': len(mod_data),
            'MCF7_mean_log2FC': mod_data['MCF7_log2FC'].mean(),
            'MCF7_sem_log2FC': mod_data['MCF7_log2FC'].sem(),
            'T47D_mean_log2FC': mod_data['T47D_log2FC'].mean(),
            'T47D_sem_log2FC': mod_data['T47D_log2FC'].sem(),
        })
    module_summary_df = pd.DataFrame(module_summary)
    module_summary_df.to_csv(OUTPUT_DIR / "tables" / "proteostasis_module_summary.csv", index=False)
    print(f"\nModule summary:")
    print(module_summary_df.to_string(index=False))

    # Get MDK data
    mdk_data = sig_df[sig_df['gene'] == 'MDK']
    if len(mdk_data) == 0:
        # MDK may not be in signature genes — compute from TPM
        if 'MDK' in tpm.index:
            mcf7_wt = tpm.loc['MDK', GROUPS['MCF7_WT']].values.astype(float)
            mcf7_d538g = tpm.loc['MDK', GROUPS['MCF7_D538G']].values.astype(float)
            t47d_wt = tpm.loc['MDK', GROUPS['T47D_WT']].values.astype(float)
            t47d_d538g = tpm.loc['MDK', GROUPS['T47D_D538G']].values.astype(float)
            mdk_mcf7_log2fc = np.log2(np.mean(mcf7_d538g) / np.mean(mcf7_wt))
            mdk_t47d_log2fc = np.log2(np.mean(t47d_d538g) / np.mean(t47d_wt))
            mdk_mcf7_sem = np.std(np.log2(mcf7_d538g / np.mean(mcf7_wt)), ddof=1) / 2
            mdk_t47d_sem = np.std(np.log2(t47d_d538g / np.mean(t47d_wt)), ddof=1) / 2
        else:
            mdk_mcf7_log2fc = mdk_t47d_log2fc = 0
            mdk_mcf7_sem = mdk_t47d_sem = 0
    else:
        mdk_row = mdk_data.iloc[0]
        mdk_mcf7_log2fc = mdk_row['MCF7_log2FC']
        mdk_t47d_log2fc = mdk_row['T47D_log2FC']
        # Compute SEM from replicates
        mcf7_wt = tpm.loc['MDK', GROUPS['MCF7_WT']].values.astype(float)
        mcf7_d538g = tpm.loc['MDK', GROUPS['MCF7_D538G']].values.astype(float)
        t47d_wt = tpm.loc['MDK', GROUPS['T47D_WT']].values.astype(float)
        t47d_d538g = tpm.loc['MDK', GROUPS['T47D_D538G']].values.astype(float)
        # SEM of log2(replicate/WT_mean)
        mdk_mcf7_sem = np.std(np.log2(mcf7_d538g / np.mean(mcf7_wt)), ddof=1) / 2
        mdk_t47d_sem = np.std(np.log2(t47d_d538g / np.mean(t47d_wt)), ddof=1) / 2

    # ---- FIGURE ----
    print("\nGenerating Figure 15...")
    fig = plt.figure(figsize=(20, 12))

    # Panel A: Secretory module co-regulation
    ax_a = fig.add_subplot(2, 3, 1)
    sec_genes = SIGNATURES['secretory_chaperones'] + ['HSP90B1']
    sec_data = sig_df[sig_df['gene'].isin(sec_genes)].copy()
    sec_data = sec_data.set_index('gene').reindex(sec_genes).reset_index()

    x = np.arange(len(sec_data))
    width = 0.35
    ax_a.bar(x - width/2, sec_data['MCF7_log2FC'], width, color=COLORS['MCF7'],
             label='MCF7', edgecolor='black', linewidth=0.5)
    ax_a.bar(x + width/2, sec_data['T47D_log2FC'], width, color=COLORS['T47D'],
             label='T47D', edgecolor='black', linewidth=0.5)

    # Highlight HSP90B1
    hsp_idx = list(sec_data['gene']).index('HSP90B1') if 'HSP90B1' in list(sec_data['gene']) else -1
    if hsp_idx >= 0:
        ax_a.bar(hsp_idx - width/2, sec_data.iloc[hsp_idx]['MCF7_log2FC'], width,
                color=COLORS['hsp90b1'], edgecolor='black', linewidth=0.5)
        ax_a.bar(hsp_idx + width/2, sec_data.iloc[hsp_idx]['T47D_log2FC'], width,
                color=COLORS['hsp90b1'], edgecolor='black', linewidth=0.5, alpha=0.6)

    ax_a.set_xticks(x)
    ax_a.set_xticklabels(sec_data['gene'], rotation=45, ha='right', fontsize=9)
    ax_a.axhline(0, color='black', linewidth=0.5)
    ax_a.set_ylabel('log2(D538G / WT)')
    ax_a.legend(fontsize=9)
    ax_a.set_title('Secretory Chaperone Co-regulation\n(93% concordant, p=0.0009)', fontsize=11)
    ax_a.text(-0.05, 1.05, 'A', transform=ax_a.transAxes, **PANEL_LABEL_PROPS)

    # Panel B: UPR marker direction
    ax_b = fig.add_subplot(2, 3, 2)
    upr_genes = SIGNATURES['upr_genes']
    upr_data = sig_df[sig_df['gene'].isin(upr_genes)].copy()
    upr_data = upr_data.set_index('gene').reindex(upr_genes).reset_index()

    bar_colors = []
    for _, row in upr_data.iterrows():
        if row['MCF7_log2FC'] > 0.15:
            bar_colors.append('#CD5C5C')  # red = stress activation
        elif row['MCF7_log2FC'] < -0.15:
            bar_colors.append('#4682B4')  # blue = stress relief
        else:
            bar_colors.append('#808080')  # gray = NC

    ax_b.bar(np.arange(len(upr_data)), upr_data['MCF7_log2FC'],
             color=bar_colors, edgecolor='black', linewidth=0.5)
    ax_b.set_xticks(np.arange(len(upr_data)))
    ax_b.set_xticklabels(upr_data['gene'], rotation=45, ha='right', fontsize=9)
    ax_b.axhline(0, color='black', linewidth=0.5)
    ax_b.set_ylabel('MCF7 log2(D538G / WT)')
    ax_b.set_title('UPR Markers in MCF7-D538G\n(3 UP, 0 DOWN; p=0.04)', fontsize=11)

    # Legend
    from matplotlib.patches import Patch
    ax_b.legend(handles=[
        Patch(facecolor='#CD5C5C', label='Stress activation (UP)'),
        Patch(facecolor='#808080', label='No change'),
        Patch(facecolor='#4682B4', label='Stress relief (DOWN)'),
    ], fontsize=8, loc='upper right')
    ax_b.text(-0.05, 1.05, 'B', transform=ax_b.transAxes, **PANEL_LABEL_PROPS)

    # Panel C: ERAD gene heatmap
    ax_c = fig.add_subplot(2, 3, 3)
    erad_genes = SIGNATURES['erad_genes']
    erad_data = sig_df[sig_df['gene'].isin(erad_genes)].copy()
    erad_data = erad_data.set_index('gene').reindex(erad_genes)

    hm = erad_data[['MCF7_log2FC', 'T47D_log2FC']].copy()
    hm.columns = ['MCF7', 'T47D']

    # Annotate with FC values
    annot = erad_data[['MCF7_FC', 'T47D_FC']].copy()
    annot.columns = ['MCF7', 'T47D']

    sns.heatmap(hm, cmap='RdBu_r', center=0, annot=annot, fmt='.2f',
                ax=ax_c, cbar_kws={'label': 'log2FC', 'shrink': 0.6},
                linewidths=0.5, linecolor='white')
    ax_c.set_title('ERAD Gene Response to D538G\n(GSE89888)', fontsize=11)
    ax_c.set_ylabel('')
    ax_c.text(-0.05, 1.05, 'C', transform=ax_c.transAxes, **PANEL_LABEL_PROPS)

    # Panel D: MDK mRNA anti-correlation
    ax_d = fig.add_subplot(2, 3, 4)
    x_pos = [0, 1]
    mdk_fcs = [mdk_mcf7_log2fc, mdk_t47d_log2fc]
    mdk_sems = [mdk_mcf7_sem, mdk_t47d_sem]
    bar_colors_mdk = [COLORS['MCF7'], COLORS['T47D']]

    bars = ax_d.bar(x_pos, mdk_fcs, yerr=mdk_sems, capsize=5,
                    color=bar_colors_mdk, edgecolor='black', linewidth=0.5, width=0.5)
    ax_d.set_xticks(x_pos)
    ax_d.set_xticklabels(['MCF7', 'T47D'], fontsize=11)
    ax_d.axhline(0, color='black', linewidth=0.5)
    ax_d.set_ylabel('MDK mRNA log2(D538G / WT)')

    # Annotate with secretion direction
    ax_d.annotate('Secretion +83%', xy=(0, mdk_fcs[0]), xytext=(0, mdk_fcs[0] - 0.15),
                 ha='center', fontsize=9, color='green', fontweight='bold')
    ax_d.annotate('Secretion −62%', xy=(1, mdk_fcs[1]), xytext=(1, mdk_fcs[1] + 0.1),
                 ha='center', fontsize=9, color='red', fontweight='bold')

    # p-values
    ax_d.text(0, mdk_fcs[0] + mdk_sems[0] + 0.05, 'p=0.0008', ha='center', fontsize=8)
    ax_d.text(1, mdk_fcs[1] + mdk_sems[1] + 0.05, 'p=0.002', ha='center', fontsize=8)

    ax_d.set_title('MDK mRNA Anti-correlates with Secretion\n(GSE89888 vs ELISA)', fontsize=11)
    ax_d.text(-0.05, 1.05, 'D', transform=ax_d.transAxes, **PANEL_LABEL_PROPS)

    # Panel E: Genome-wide opposite regulation scatter
    ax_e = fig.add_subplot(2, 3, 5)

    # Color by module membership
    ranking_df['module'] = ranking_df['gene'].map(module_map).fillna('Other')
    for mod_name, color, zorder in [
        ('Other', '#D3D3D3', 1),
        ('ERAD', COLORS['erad'], 3),
        ('UPR', COLORS['upr'], 3),
        ('Secretory', COLORS['secretory'], 3),
    ]:
        subset = ranking_df[ranking_df['module'] == mod_name]
        ax_e.scatter(subset['MCF7_log2FC'], subset['T47D_log2FC'],
                    c=color, s=15 if mod_name == 'Other' else 40,
                    alpha=0.4 if mod_name == 'Other' else 0.9,
                    edgecolors='black' if mod_name != 'Other' else 'none',
                    linewidths=0.3, zorder=zorder, label=mod_name)

    # Label HSP90B1 and MDK
    for label_gene in ['HSP90B1', 'MDK']:
        gene_row = ranking_df[ranking_df['gene'] == label_gene]
        if len(gene_row) > 0:
            row = gene_row.iloc[0]
            color = COLORS['hsp90b1'] if label_gene == 'HSP90B1' else 'black'
            ax_e.scatter(row['MCF7_log2FC'], row['T47D_log2FC'],
                        c=color, s=80, edgecolors='black', linewidths=1, zorder=5, marker='*')
            ax_e.annotate(label_gene, (row['MCF7_log2FC'], row['T47D_log2FC']),
                         fontsize=9, fontweight='bold',
                         xytext=(5, 5), textcoords='offset points')

    ax_e.axhline(0, color='black', linewidth=0.5, alpha=0.3)
    ax_e.axvline(0, color='black', linewidth=0.5, alpha=0.3)
    ax_e.set_xlabel('MCF7 log2(D538G / WT)')
    ax_e.set_ylabel('T47D log2(D538G / WT)')
    ax_e.legend(fontsize=8, loc='lower left')
    ax_e.set_title(f'Genome-wide Opposite Regulation\n({len(ranking_df)} genes, GSE89888)', fontsize=11)
    ax_e.text(-0.05, 1.05, 'E', transform=ax_e.transAxes, **PANEL_LABEL_PROPS)

    # Panel F: Chaperone capacity across conditions
    ax_f = fig.add_subplot(2, 3, 6)
    sec_chap_genes = SIGNATURES['secretory_chaperones']

    conditions = ['MCF7_WT', 'MCF7_D538G', 'T47D_WT', 'T47D_D538G']
    cond_labels = ['MCF7\nWT', 'MCF7\nD538G', 'T47D\nWT', 'T47D\nD538G']
    cond_colors = [COLORS['MCF7'], COLORS['MCF7'], COLORS['T47D'], COLORS['T47D']]
    cond_alphas = [0.5, 1.0, 0.5, 1.0]

    means = []
    sems = []
    for cond in conditions:
        vals = []
        for gene in sec_chap_genes:
            if gene in tpm.index:
                gene_mean = tpm.loc[gene, GROUPS[cond]].astype(float).mean()
                vals.append(gene_mean)
        means.append(np.mean(vals))
        sems.append(np.std(vals, ddof=1) / np.sqrt(len(vals)))

    bars = ax_f.bar(np.arange(4), means, yerr=sems, capsize=5,
                    color=cond_colors, edgecolor='black', linewidth=0.5, width=0.6)
    for bar, alpha in zip(bars, cond_alphas):
        bar.set_alpha(alpha)

    ax_f.set_xticks(np.arange(4))
    ax_f.set_xticklabels(cond_labels, fontsize=9)
    ax_f.set_ylabel('Mean Secretory Chaperone TPM')
    ax_f.set_title('Secretory Capacity Across Conditions\n(7 ER chaperones, GSE89888)', fontsize=11)
    ax_f.text(-0.05, 1.05, 'F', transform=ax_f.transAxes, **PANEL_LABEL_PROPS)

    plt.tight_layout()
    fig.savefig(OUTPUT_DIR / "figures" / "fig15_proteostasis_setpoint.pdf", bbox_inches='tight')
    fig.savefig(OUTPUT_DIR / "figures" / "fig15_proteostasis_setpoint.png", dpi=300, bbox_inches='tight')
    plt.close()

    print(f"\nSaved: outputs/figures/fig15_proteostasis_setpoint.pdf/png")

    print("\n" + "=" * 80)
    print("SCRIPT 15 COMPLETE")
    print("=" * 80)
    return 0


if __name__ == "__main__":
    sys.exit(main())
```

**Step 2: Verify syntax**

Run: `python -c "import ast; ast.parse(open('scripts/15_proteostasis_setpoint.py').read()); print('OK')"`
Expected: `OK`

**Step 3: Commit**

```bash
git add midkine/mdk_saturation_pipeline/scripts/15_proteostasis_setpoint.py
git commit -m "feat: add script 15 - proteostasis setpoint figure"
```

---

## Task 3: Script 16 — Functional Consequence Figure

**Files:**
- Create: `midkine/mdk_saturation_pipeline/scripts/16_functional_consequence.py`
- Read: `midkine/mdk_saturation_pipeline/outputs/tables/mechanism_gene_signatures.csv`
- Read: `midkine/mdk_saturation_pipeline/outputs/tables/mechanism_discrimination_scorecard.csv`
- Output: `midkine/mdk_saturation_pipeline/outputs/figures/fig16_functional_consequence.{pdf,png}`

**Step 1: Create script 16**

```python
#!/usr/bin/env python
"""
16_functional_consequence.py

Secretory capacity, UPR tone, and stress response diverge to determine protein fate.

Panels:
  A) Secretory module vs UPR trajectory (slope chart)
  B) Three-module heatmap (ERAD + UPR + Secretory, grouped)
  C) MDK mRNA vs secretory module direction (scatter, 4 conditions)
  D) Integrated evidence table

Datasets: GSE89888, script 13 outputs
Outputs:  outputs/figures/fig16_functional_consequence.pdf/png
"""

# [INSERT SHARED PREAMBLE HERE]


def load_gene_signatures():
    """Load script 13 gene signature data."""
    path = OUTPUT_DIR / "tables" / "mechanism_gene_signatures.csv"
    if not path.exists():
        sys.exit(f"Missing: {path}. Run script 13 first.")
    return pd.read_csv(path)


def main():
    print("=" * 80)
    print("SCRIPT 16: FUNCTIONAL CONSEQUENCE FIGURE")
    print("=" * 80)

    sig_df = load_gene_signatures()
    tpm = load_tpm()

    # Compute module-level means using absolute TPM across 4 conditions
    sec_genes = [g for g in SIGNATURES['secretory_chaperones'] if g in tpm.index]
    upr_genes = [g for g in SIGNATURES['upr_genes'] if g in tpm.index]
    erad_genes = [g for g in SIGNATURES['erad_genes'] if g in tpm.index]

    conditions = ['MCF7_WT', 'MCF7_D538G', 'T47D_WT', 'T47D_D538G']

    def module_mean_tpm(genes, condition):
        vals = [tpm.loc[g, GROUPS[condition]].astype(float).mean() for g in genes]
        return np.mean(vals)

    def module_log2fc(genes, cell_line):
        wt_key = f'{cell_line}_WT'
        d538g_key = f'{cell_line}_D538G'
        fcs = []
        for g in genes:
            wt = tpm.loc[g, GROUPS[wt_key]].astype(float).mean()
            d538g = tpm.loc[g, GROUPS[d538g_key]].astype(float).mean()
            if wt > 0:
                fcs.append(np.log2(d538g / wt))
        return np.mean(fcs), np.std(fcs, ddof=1) / np.sqrt(len(fcs)) if len(fcs) > 1 else 0

    sec_mcf7_fc, sec_mcf7_sem = module_log2fc(sec_genes, 'MCF7')
    sec_t47d_fc, sec_t47d_sem = module_log2fc(sec_genes, 'T47D')
    upr_mcf7_fc, upr_mcf7_sem = module_log2fc(upr_genes, 'MCF7')
    upr_t47d_fc, upr_t47d_sem = module_log2fc(upr_genes, 'T47D')
    erad_mcf7_fc, erad_mcf7_sem = module_log2fc(erad_genes, 'MCF7')
    erad_t47d_fc, erad_t47d_sem = module_log2fc(erad_genes, 'T47D')

    # MDK log2FC
    mdk_mcf7_fc = 0
    mdk_t47d_fc = 0
    if 'MDK' in tpm.index:
        mcf7_wt = tpm.loc['MDK', GROUPS['MCF7_WT']].astype(float).mean()
        mcf7_d538g = tpm.loc['MDK', GROUPS['MCF7_D538G']].astype(float).mean()
        t47d_wt = tpm.loc['MDK', GROUPS['T47D_WT']].astype(float).mean()
        t47d_d538g = tpm.loc['MDK', GROUPS['T47D_D538G']].astype(float).mean()
        mdk_mcf7_fc = np.log2(mcf7_d538g / mcf7_wt) if mcf7_wt > 0 else 0
        mdk_t47d_fc = np.log2(t47d_d538g / t47d_wt) if t47d_wt > 0 else 0

    print(f"Module log2FC (D538G/WT):")
    print(f"  Secretory: MCF7={sec_mcf7_fc:.3f}, T47D={sec_t47d_fc:.3f}")
    print(f"  UPR:       MCF7={upr_mcf7_fc:.3f}, T47D={upr_t47d_fc:.3f}")
    print(f"  ERAD:      MCF7={erad_mcf7_fc:.3f}, T47D={erad_t47d_fc:.3f}")
    print(f"  MDK mRNA:  MCF7={mdk_mcf7_fc:.3f}, T47D={mdk_t47d_fc:.3f}")

    # ---- FIGURE ----
    print("\nGenerating Figure 16...")
    fig = plt.figure(figsize=(16, 12))

    # Panel A: Secretory vs UPR trajectory (slope chart)
    ax_a = fig.add_subplot(2, 2, 1)

    modules = ['Secretory', 'UPR', 'ERAD']
    mcf7_vals = [sec_mcf7_fc, upr_mcf7_fc, erad_mcf7_fc]
    t47d_vals = [sec_t47d_fc, upr_t47d_fc, erad_t47d_fc]
    mcf7_sems = [sec_mcf7_sem, upr_mcf7_sem, erad_mcf7_sem]
    t47d_sems = [sec_t47d_sem, upr_t47d_sem, erad_t47d_sem]
    mod_colors = [COLORS['secretory'], COLORS['upr'], COLORS['erad']]

    x = np.arange(len(modules))
    width = 0.3
    bars1 = ax_a.bar(x - width/2, mcf7_vals, width, yerr=mcf7_sems, capsize=4,
                     color=COLORS['MCF7'], edgecolor='black', linewidth=0.5, label='MCF7')
    bars2 = ax_a.bar(x + width/2, t47d_vals, width, yerr=t47d_sems, capsize=4,
                     color=COLORS['T47D'], edgecolor='black', linewidth=0.5, label='T47D')

    # Color the module labels
    ax_a.set_xticks(x)
    labels = ax_a.set_xticklabels(modules, fontsize=11)
    for label, color in zip(labels, mod_colors):
        label.set_color(color)
        label.set_fontweight('bold')

    ax_a.axhline(0, color='black', linewidth=0.5)
    ax_a.set_ylabel('Mean log2(D538G / WT)')
    ax_a.legend(fontsize=10)
    ax_a.set_title('Module-Level Response to D538G\n(GSE89888)', fontsize=11)
    ax_a.text(-0.05, 1.05, 'A', transform=ax_a.transAxes, **PANEL_LABEL_PROPS)

    # Panel B: Three-module heatmap (grouped, not clustered)
    ax_b = fig.add_subplot(2, 2, 2)

    all_genes_ordered = erad_genes + upr_genes + sec_genes
    hm_data = []
    hm_genes = []
    hm_modules = []
    for gene in all_genes_ordered:
        row = sig_df[sig_df['gene'] == gene]
        if len(row) == 0:
            continue
        row = row.iloc[0]
        hm_data.append([row['MCF7_log2FC'], row['T47D_log2FC']])
        hm_genes.append(gene)
        if gene in SIGNATURES['erad_genes']:
            hm_modules.append('ERAD')
        elif gene in SIGNATURES['upr_genes']:
            hm_modules.append('UPR')
        else:
            hm_modules.append('Secretory')

    hm_df = pd.DataFrame(hm_data, index=hm_genes, columns=['MCF7', 'T47D'])

    # Annotate with FC values
    annot_df = pd.DataFrame(index=hm_genes, columns=['MCF7', 'T47D'])
    for gene in hm_genes:
        row = sig_df[sig_df['gene'] == gene].iloc[0]
        annot_df.loc[gene, 'MCF7'] = f"{row['MCF7_FC']:.2f}"
        annot_df.loc[gene, 'T47D'] = f"{row['T47D_FC']:.2f}"

    sns.heatmap(hm_df, cmap='RdBu_r', center=0, annot=annot_df, fmt='',
                ax=ax_b, cbar_kws={'label': 'log2FC', 'shrink': 0.6},
                linewidths=0.5, linecolor='white')

    # Add module color bar on left
    mod_color_map = {'ERAD': COLORS['erad'], 'UPR': COLORS['upr'], 'Secretory': COLORS['secretory']}
    for i, mod in enumerate(hm_modules):
        ax_b.add_patch(plt.Rectangle((-0.3, i), 0.25, 1, color=mod_color_map[mod],
                       transform=ax_b.transData, clip_on=False))

    # Add horizontal lines between modules
    erad_end = len([m for m in hm_modules if m == 'ERAD'])
    upr_end = erad_end + len([m for m in hm_modules if m == 'UPR'])
    ax_b.axhline(erad_end, color='white', linewidth=2)
    ax_b.axhline(upr_end, color='white', linewidth=2)

    ax_b.set_title('Three-Module Gene Heatmap\n(D538G / WT, GSE89888)', fontsize=11)
    ax_b.set_ylabel('')
    ax_b.text(-0.15, 1.05, 'B', transform=ax_b.transAxes, **PANEL_LABEL_PROPS)

    # Panel C: MDK mRNA vs secretory module (4 conditions, absolute TPM)
    ax_c = fig.add_subplot(2, 2, 3)

    # Use absolute TPM for 4 real data points
    sec_tpm_vals = []
    mdk_tpm_vals = []
    point_labels = []
    point_colors = []
    point_markers = []

    for cond, label, color, marker in [
        ('MCF7_WT', 'MCF7 WT', COLORS['MCF7'], 'o'),
        ('MCF7_D538G', 'MCF7 D538G', COLORS['MCF7'], 's'),
        ('T47D_WT', 'T47D WT', COLORS['T47D'], 'o'),
        ('T47D_D538G', 'T47D D538G', COLORS['T47D'], 's'),
    ]:
        sec_mean = module_mean_tpm(sec_genes, cond)
        mdk_mean = tpm.loc['MDK', GROUPS[cond]].astype(float).mean() if 'MDK' in tpm.index else 0
        sec_tpm_vals.append(sec_mean)
        mdk_tpm_vals.append(mdk_mean)
        point_labels.append(label)
        point_colors.append(color)
        point_markers.append(marker)

    for i in range(4):
        ax_c.scatter(sec_tpm_vals[i], mdk_tpm_vals[i], c=point_colors[i],
                    marker=point_markers[i], s=100, edgecolors='black', linewidths=1, zorder=3)
        ax_c.annotate(point_labels[i], (sec_tpm_vals[i], mdk_tpm_vals[i]),
                     fontsize=9, xytext=(8, 5), textcoords='offset points')

    # Trend line
    r, p = stats.pearsonr(sec_tpm_vals, mdk_tpm_vals)
    z = np.polyfit(sec_tpm_vals, mdk_tpm_vals, 1)
    poly = np.poly1d(z)
    x_line = np.linspace(min(sec_tpm_vals) * 0.9, max(sec_tpm_vals) * 1.1, 50)
    ax_c.plot(x_line, poly(x_line), '--', color='gray', alpha=0.5)
    ax_c.text(0.05, 0.95, f'r = {r:.2f}, p = {p:.3f}', transform=ax_c.transAxes,
             fontsize=10, va='top')

    # Legend for markers
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor='gray', markersize=8, label='WT'),
        Line2D([0], [0], marker='s', color='w', markerfacecolor='gray', markersize=8, label='D538G'),
    ]
    ax_c.legend(handles=legend_elements, fontsize=9, loc='lower right')

    ax_c.set_xlabel('Mean Secretory Chaperone TPM')
    ax_c.set_ylabel('MDK mRNA TPM')
    ax_c.set_title('MDK Expression vs Secretory Capacity\n(4 conditions, GSE89888)', fontsize=11)
    ax_c.text(-0.05, 1.05, 'C', transform=ax_c.transAxes, **PANEL_LABEL_PROPS)

    # Panel D: Integrated evidence table
    ax_d = fig.add_subplot(2, 2, 4)
    ax_d.axis('off')

    table_data = [
        ['Different\ncistrome',
         '168/632 opposite genes\nhave diff. ER binding',
         'GSE125117\n+ GSE89888',
         'p < 0.05\n(Fisher)',
         'CHX chase:\nT47D-D538G shorter\nMDK half-life'],
        ['Different\nproteostasis',
         '93% secretory co-regulation\nUPR activation (3/5 UP)',
         'GSE89888',
         'p = 0.0009\np = 0.04',
         'Pulse-chase:\nMCF7-D538G higher\nfractional secretion'],
        ['Module\ncoordination',
         'MDK mRNA anti-correlates\nwith protein secretion',
         'GSE89888\n+ ELISA',
         'p = 0.0008\np = 0.002',
         'Proteasome/lysosome\ninhibitors rescue\nT47D MDK'],
    ]

    col_labels = ['Context', 'Key Finding', 'Dataset', 'Statistic', 'Experimental\nPrediction']
    row_colors = [COLORS['secretory'], COLORS['upr'], COLORS['erad']]

    table = ax_d.table(
        cellText=table_data,
        colLabels=col_labels,
        cellLoc='center',
        loc='center',
        colWidths=[0.12, 0.28, 0.13, 0.13, 0.22],
    )
    table.auto_set_font_size(False)
    table.set_fontsize(8)
    table.scale(1.2, 2.0)

    # Color header row
    for j in range(5):
        table[0, j].set_facecolor('#2C3E50')
        table[0, j].set_text_props(color='white', fontweight='bold')

    # Color first column by "why"
    for i, color in enumerate(row_colors):
        table[i + 1, 0].set_facecolor(color)
        table[i + 1, 0].set_text_props(color='white', fontweight='bold')

    ax_d.set_title('Integrated Evidence Summary', fontsize=11, pad=20)
    ax_d.text(-0.05, 1.05, 'D', transform=ax_d.transAxes, **PANEL_LABEL_PROPS)

    plt.tight_layout()
    fig.savefig(OUTPUT_DIR / "figures" / "fig16_functional_consequence.pdf", bbox_inches='tight')
    fig.savefig(OUTPUT_DIR / "figures" / "fig16_functional_consequence.png", dpi=300, bbox_inches='tight')
    plt.close()

    print(f"\nSaved: outputs/figures/fig16_functional_consequence.pdf/png")

    print("\n" + "=" * 80)
    print("SCRIPT 16 COMPLETE")
    print("=" * 80)
    return 0


if __name__ == "__main__":
    sys.exit(main())
```

**Step 2: Verify syntax**

Run: `python -c "import ast; ast.parse(open('scripts/16_functional_consequence.py').read()); print('OK')"`
Expected: `OK`

**Step 3: Commit**

```bash
git add midkine/mdk_saturation_pipeline/scripts/16_functional_consequence.py
git commit -m "feat: add script 16 - functional consequence figure"
```

---

## Task 4: Pipeline Integration and SLURM Scripts

**Files:**
- Modify: `midkine/mdk_saturation_pipeline/run_pipeline.py` (line 31, add 3 scripts to SCRIPTS list)
- Create: `midkine/mdk_saturation_pipeline/slurm/run_scripts14_16.sh`

**Step 1: Update run_pipeline.py**

Add after `"13_discriminate_hsp90b1_mechanism.py"`:

```python
    "14_cistrome_divergence.py",
    "15_proteostasis_setpoint.py",
    "16_functional_consequence.py",
```

**Step 2: Create SLURM script**

```bash
#!/bin/bash
#SBATCH --job-name=fig14_16
#SBATCH --output=slurm_log/fig14_16_%j.out
#SBATCH --error=slurm_log/fig14_16_%j.err
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --cluster=htc
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alc376@pitt.edu

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/midkine/mdk_saturation_pipeline

module load gcc/8.2.0
source /ix1/alee/LO_LAB/Personal/Alexander_Chang/miniconda3/etc/profile.d/conda.sh
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

echo "Running Scripts 14-16: Attractor State Evidence Figures"
echo "======================================================="
date

python scripts/14_cistrome_divergence.py && \
python scripts/15_proteostasis_setpoint.py && \
python scripts/16_functional_consequence.py

echo ""
echo "All figure scripts finished"
date
```

**Step 3: Commit**

```bash
git add midkine/mdk_saturation_pipeline/run_pipeline.py \
        midkine/mdk_saturation_pipeline/slurm/run_scripts14_16.sh
git commit -m "feat: register scripts 14-16 in pipeline runner and add SLURM script"
```

---

## Task 5: Run Scripts and Verify Outputs

**Step 1: Submit via SLURM**

```bash
cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/midkine/mdk_saturation_pipeline
sbatch slurm/run_scripts14_16.sh
```

**Step 2: Monitor and verify**

Check SLURM output log for:
- "SCRIPT 14 COMPLETE"
- "SCRIPT 15 COMPLETE"
- "SCRIPT 16 COMPLETE"

Verify 8 output files exist:
- `outputs/figures/fig14_cistrome_divergence.{pdf,png}`
- `outputs/figures/fig15_proteostasis_setpoint.{pdf,png}`
- `outputs/figures/fig16_functional_consequence.{pdf,png}`
- `outputs/tables/cistrome_intersection_go_enrichment.csv`
- `outputs/tables/proteostasis_module_summary.csv`

**Step 3: Fix any runtime errors and re-submit if needed**

---

## Task 6: Final Review

**Step 1: Visual review of all 3 figures**

Read each PNG to verify:
- Panel labels (A, B, C, D...) are present and correctly placed
- Colors match the palette (MCF7=steelblue, T47D=indianred, etc.)
- Titles are accurate
- Statistical annotations are present
- No overlapping text or truncated labels

**Step 2: Commit any fixes**

```bash
git add -A
git commit -m "fix: figure adjustments from visual review"
```
