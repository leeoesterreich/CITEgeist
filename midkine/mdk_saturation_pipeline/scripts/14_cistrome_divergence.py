#!/usr/bin/env python
"""
14_cistrome_divergence.py

CISTROME DIVERGENCE: ESR1-D538G ENGAGES DISTINCT REGULATORY PROGRAMS
IN MCF7 vs T47D

Purpose:
Generates a 4-panel grant-quality figure demonstrating that the D538G
mutation activates fundamentally different transcriptional programs in
MCF7 vs T47D, by integrating expression changes (script 11 venn),
ER ChIP-seq binding, ATAC-seq accessibility, and GO enrichment.

Panels:
  A: Venn diagram - opposite-expression genes intersected with
     differential ER binding genes
  B: Heatmap of ER binding changes at intersection gene promoters
  C: ATAC-seq scatter - peak counts near intersection gene TSS
  D: GO enrichment of intersection gene subgroups

Datasets:
  outputs/tables/venn_expression_binding.csv (from script 11)
  outputs/tables/unsupervised_hsp90b1_ranking.csv (background genes)
  GSE254216 - ATAC-seq summit BED files

Outputs:
  outputs/figures/fig14_cistrome_divergence.pdf
  outputs/figures/fig14_cistrome_divergence.png
  outputs/tables/cistrome_intersection_go_enrichment.csv
"""

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
from matplotlib_venn import venn2
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

# Style constants
COLORS = {
    'MCF7': '#4682B4', 'T47D': '#CD5C5C',
    'secretory': '#E8890C', 'upr': '#2CA02C', 'erad': '#7B2D8E',
    'hsp90b1': '#DAA520',
    'same_as_secretion': '#E8890C', 'inverse_repressor': '#7B2D8E',
}
PANEL_LABEL_PROPS = dict(fontsize=14, fontweight='bold', va='top', ha='left')


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

def load_venn_data():
    """Load venn expression/binding data from script 11."""
    path = OUTPUT_DIR / "tables" / "venn_expression_binding.csv"
    df = pd.read_csv(path)
    print(f"Loaded venn data: {len(df)} genes")
    return df


def load_atac_peaks(bed_path):
    """Load ATAC-seq summit BED file (gzipped, 3-column: chr, start, end)."""
    peaks = []
    with gzip.open(bed_path, 'rt') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                try:
                    peaks.append((parts[0], int(parts[1]), int(parts[2])))
                except ValueError:
                    continue
    print(f"Loaded {len(peaks)} peaks from {bed_path.name}")
    return peaks


def count_peaks_near_tss(peaks, chrom, tss, window=10000):
    """Count peaks within +/-window of TSS on the given chromosome."""
    count = 0
    for p_chr, p_start, p_end in peaks:
        if p_chr != chrom:
            continue
        p_mid = (p_start + p_end) / 2
        if abs(p_mid - tss) <= window:
            count += 1
    return count


def run_go_enrichment(gene_list, background_list, source='GO:BP', top_n=5):
    """Run GO enrichment via gprofiler-official. Returns DataFrame."""
    try:
        from gprofiler import GProfiler
    except ImportError:
        print("WARNING: gprofiler-official not installed, skipping GO enrichment")
        return pd.DataFrame()

    gp = GProfiler(return_dataframe=True)
    try:
        result = gp.profile(
            organism='hsapiens',
            query=list(gene_list),
            background=list(background_list),
            sources=[source],
            significance_threshold_method='fdr',
            user_threshold=0.05,
            no_evidences=True,
        )
    except Exception as e:
        print(f"WARNING: GO enrichment failed: {e}")
        return pd.DataFrame()

    if result.empty:
        return pd.DataFrame()

    result = result.sort_values('p_value').head(top_n)
    return result[['native', 'name', 'p_value', 'intersection_size',
                    'term_size']].copy()


# ---------------------------------------------------------------------------
# Panel A: Venn diagram
# ---------------------------------------------------------------------------

def panel_a(ax, venn_df):
    """Venn diagram: opposite expression vs differential ER binding."""
    n_opposite = len(venn_df)  # full venn_df = opposite expression genes
    diff_binding = venn_df[venn_df['has_differential_binding'] == True]
    n_diff_binding = len(diff_binding)

    # Intersection is genes with BOTH opposite expression AND diff binding
    n_intersection = n_diff_binding  # diff binding is subset of venn_df
    n_only_expression = n_opposite - n_intersection
    # For the second circle, we use the same count since all diff_binding
    # genes are already in the opposite expression set
    n_only_binding = 0

    v = venn2(
        subsets=(n_only_expression, n_only_binding, n_intersection),
        set_labels=('Opposite\nExpression', 'Differential\nER Binding'),
        ax=ax,
    )

    # Style the patches
    if v.get_patch_by_id('10'):
        v.get_patch_by_id('10').set_color(COLORS['MCF7'])
        v.get_patch_by_id('10').set_alpha(0.4)
    if v.get_patch_by_id('01'):
        v.get_patch_by_id('01').set_color(COLORS['T47D'])
        v.get_patch_by_id('01').set_alpha(0.4)
    if v.get_patch_by_id('11'):
        v.get_patch_by_id('11').set_color(COLORS['hsp90b1'])
        v.get_patch_by_id('11').set_alpha(0.5)

    ax.set_title('Expression-Binding Intersection', fontsize=11,
                 fontweight='bold')
    ax.text(-0.05, 1.05, 'A', transform=ax.transAxes, **PANEL_LABEL_PROPS)


# ---------------------------------------------------------------------------
# Panel B: ER binding heatmap
# ---------------------------------------------------------------------------

def panel_b(ax, intersection_df):
    """Heatmap of ER binding changes at intersection gene promoters."""
    df = intersection_df.sort_values('total_effect', ascending=False).copy()

    # Prepare matrix
    mat = df[['MCF7_binding_diff', 'T47D_binding_diff']].values
    clip_val = np.percentile(np.abs(mat[np.isfinite(mat)]), 95)
    mat = np.clip(mat, -clip_val, clip_val)

    # Correlation type side colors
    ctype_colors = df['correlation_type'].map({
        'same_as_secretion': COLORS['same_as_secretion'],
        'inverse_repressor': COLORS['inverse_repressor'],
    }).fillna('#CCCCCC').values

    # Plot heatmap
    im = ax.imshow(mat, aspect='auto', cmap='RdBu_r',
                   vmin=-clip_val, vmax=clip_val,
                   interpolation='nearest')

    # Side color bar for correlation_type
    n_genes = len(df)
    for i, c in enumerate(ctype_colors):
        ax.add_patch(plt.Rectangle((-0.7, i - 0.5), 0.4, 1.0,
                                   color=c, clip_on=False))

    ax.set_xticks([0, 1])
    ax.set_xticklabels(['MCF7\nbinding diff', 'T47D\nbinding diff'],
                       fontsize=9)
    ax.set_ylabel(f'Genes (n={n_genes})', fontsize=10)
    ax.set_yticks([])
    ax.set_title('ER Binding Changes at\nIntersection Promoters',
                 fontsize=11, fontweight='bold')

    # Colorbar
    cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label('Binding difference', fontsize=9)

    ax.text(-0.15, 1.05, 'B', transform=ax.transAxes, **PANEL_LABEL_PROPS)


# ---------------------------------------------------------------------------
# Panel C: ATAC-seq scatter
# ---------------------------------------------------------------------------

def panel_c(ax, intersection_df, mcf7_peaks, t47d_peaks):
    """ATAC-seq peak count scatter for intersection genes."""
    mcf7_counts = []
    t47d_counts = []
    colors = []

    for _, row in intersection_df.iterrows():
        chrom = row['chr']
        tss = row['tss']
        if pd.isna(chrom) or pd.isna(tss):
            mcf7_counts.append(0)
            t47d_counts.append(0)
            colors.append('#CCCCCC')
            continue

        mc = count_peaks_near_tss(mcf7_peaks, str(chrom), float(tss))
        tc = count_peaks_near_tss(t47d_peaks, str(chrom), float(tss))
        mcf7_counts.append(mc)
        t47d_counts.append(tc)

        ctype = row.get('correlation_type', '')
        colors.append(COLORS.get(ctype, '#CCCCCC'))

    mcf7_counts = np.array(mcf7_counts)
    t47d_counts = np.array(t47d_counts)

    ax.scatter(mcf7_counts, t47d_counts, c=colors, alpha=0.6,
               edgecolors='k', linewidths=0.3, s=30)

    # Diagonal line
    max_val = max(mcf7_counts.max(), t47d_counts.max()) if len(mcf7_counts) > 0 else 10
    ax.plot([0, max_val], [0, max_val], 'k--', alpha=0.4, linewidth=1)

    ax.set_xlabel('MCF7 ATAC-seq peaks (within 10kb)', fontsize=10)
    ax.set_ylabel('T47D ATAC-seq peaks (within 10kb)', fontsize=10)
    ax.set_title('Chromatin Accessibility\nat Intersection Genes',
                 fontsize=11, fontweight='bold')

    # Legend
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker='o', color='w',
               markerfacecolor=COLORS['same_as_secretion'],
               markersize=8, label='Same as secretion'),
        Line2D([0], [0], marker='o', color='w',
               markerfacecolor=COLORS['inverse_repressor'],
               markersize=8, label='Inverse repressor'),
    ]
    ax.legend(handles=legend_elements, loc='upper left', fontsize=8,
              framealpha=0.8)

    ax.text(-0.05, 1.05, 'C', transform=ax.transAxes, **PANEL_LABEL_PROPS)


# ---------------------------------------------------------------------------
# Panel D: GO enrichment
# ---------------------------------------------------------------------------

def panel_d(ax, intersection_df, background_genes):
    """GO enrichment bar chart for intersection gene subgroups."""
    same_genes = intersection_df[
        intersection_df['correlation_type'] == 'same_as_secretion'
    ]['gene'].dropna().unique().tolist()

    inverse_genes = intersection_df[
        intersection_df['correlation_type'] == 'inverse_repressor'
    ]['gene'].dropna().unique().tolist()

    print(f"GO enrichment: {len(same_genes)} same_as_secretion, "
          f"{len(inverse_genes)} inverse_repressor genes")

    bg = list(background_genes)

    go_same = run_go_enrichment(same_genes, bg) if len(same_genes) >= 3 else pd.DataFrame()
    go_inverse = run_go_enrichment(inverse_genes, bg) if len(inverse_genes) >= 3 else pd.DataFrame()

    # Combine results for saving
    all_go = []
    if not go_same.empty:
        go_same = go_same.copy()
        go_same['group'] = 'same_as_secretion'
        all_go.append(go_same)
    if not go_inverse.empty:
        go_inverse = go_inverse.copy()
        go_inverse['group'] = 'inverse_repressor'
        all_go.append(go_inverse)

    if all_go:
        go_combined = pd.concat(all_go, ignore_index=True)
        go_save_path = OUTPUT_DIR / "tables" / "cistrome_intersection_go_enrichment.csv"
        go_combined.to_csv(go_save_path, index=False)
        print(f"Saved GO enrichment: {go_save_path}")
    else:
        go_combined = pd.DataFrame()
        # Save empty file
        go_save_path = OUTPUT_DIR / "tables" / "cistrome_intersection_go_enrichment.csv"
        go_combined.to_csv(go_save_path, index=False)
        print(f"No significant GO terms found; saved empty table.")

    # Plot
    if go_combined.empty:
        ax.text(0.5, 0.5, 'No significant\nGO terms found',
                ha='center', va='center', fontsize=12,
                transform=ax.transAxes, color='gray')
        ax.set_title('GO Enrichment (GO:BP)', fontsize=11, fontweight='bold')
        ax.text(-0.05, 1.05, 'D', transform=ax.transAxes, **PANEL_LABEL_PROPS)
        return

    # Build bar data
    bar_labels = []
    bar_values = []
    bar_colors = []

    for _, row in go_combined.iterrows():
        term_name = row['name']
        if len(term_name) > 40:
            term_name = term_name[:37] + '...'
        bar_labels.append(term_name)
        bar_values.append(-np.log10(row['p_value']))
        bar_colors.append(COLORS.get(row['group'], '#CCCCCC'))

    y_pos = np.arange(len(bar_labels))
    ax.barh(y_pos, bar_values, color=bar_colors, edgecolor='k',
            linewidth=0.5, height=0.7)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(bar_labels, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel('-log10(p-value)', fontsize=10)
    ax.set_title('GO Enrichment (GO:BP)', fontsize=11, fontweight='bold')

    # Legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor=COLORS['same_as_secretion'],
              label='Same as secretion'),
        Patch(facecolor=COLORS['inverse_repressor'],
              label='Inverse repressor'),
    ]
    ax.legend(handles=legend_elements, loc='lower right', fontsize=8,
              framealpha=0.8)

    ax.text(-0.05, 1.05, 'D', transform=ax.transAxes, **PANEL_LABEL_PROPS)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("=" * 70)
    print("Script 14: Cistrome Divergence Figure")
    print("=" * 70)

    # Ensure output directories exist
    (OUTPUT_DIR / "figures").mkdir(parents=True, exist_ok=True)
    (OUTPUT_DIR / "tables").mkdir(parents=True, exist_ok=True)

    # --- Load data ---
    venn_df = load_venn_data()

    # Intersection: opposite expression AND differential ER binding
    intersection_df = venn_df[venn_df['has_differential_binding'] == True].copy()
    print(f"Intersection genes (opposite expr + diff binding): "
          f"{len(intersection_df)}")

    # Load ATAC-seq peaks
    mcf7_atac_path = DATA_DIR / config['datasets']['GSE254216_MCF7']
    t47d_atac_path = DATA_DIR / config['datasets']['GSE254216_T47D']
    mcf7_peaks = load_atac_peaks(mcf7_atac_path)
    t47d_peaks = load_atac_peaks(t47d_atac_path)

    # Background gene list for GO enrichment
    bg_path = OUTPUT_DIR / "tables" / "unsupervised_hsp90b1_ranking.csv"
    if bg_path.exists():
        bg_df = pd.read_csv(bg_path)
        # Use the gene column if present
        bg_col = 'gene' if 'gene' in bg_df.columns else bg_df.columns[0]
        background_genes = bg_df[bg_col].dropna().unique().tolist()
        print(f"Background gene list: {len(background_genes)} genes")
    else:
        # Fall back to all genes in venn_df
        background_genes = venn_df['gene'].dropna().unique().tolist()
        print(f"Using venn_df genes as background: {len(background_genes)}")

    # --- Create figure ---
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle(
        'ESR1-D538G Engages Distinct Regulatory Programs in MCF7 vs T47D',
        fontsize=14, fontweight='bold', y=0.98,
    )

    panel_a(axes[0, 0], venn_df)
    panel_b(axes[0, 1], intersection_df)
    panel_c(axes[1, 0], intersection_df, mcf7_peaks, t47d_peaks)
    panel_d(axes[1, 1], intersection_df, background_genes)

    plt.tight_layout(rect=[0, 0, 1, 0.95])

    # --- Save ---
    fig_path_pdf = OUTPUT_DIR / "figures" / "fig14_cistrome_divergence.pdf"
    fig_path_png = OUTPUT_DIR / "figures" / "fig14_cistrome_divergence.png"
    fig.savefig(fig_path_pdf, dpi=300, bbox_inches='tight')
    fig.savefig(fig_path_png, dpi=300, bbox_inches='tight')
    plt.close(fig)

    print(f"\nSaved: {fig_path_pdf}")
    print(f"Saved: {fig_path_png}")
    print("=" * 70)
    print("Script 14 complete.")


if __name__ == '__main__':
    main()
