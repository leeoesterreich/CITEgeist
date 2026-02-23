#!/usr/bin/env python3
"""
Supplementary Figure S2: Differential Expression and Pathway Enrichment

Two-panel supplementary figure (1x2 layout, 12x5 inches) showing PyDESeq2 results
from pseudo-bulk differential expression analysis (responder vs progressor).

Panel A: Volcano plot (203 significant genes, 120 responder-up, 83 progressor-up)
Panel B: Pathway enrichment dot plot (MSigDB Hallmark 2020, split by direction)
"""

import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from pathlib import Path

from figure_style import apply_style, PALETTE

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
PROJECT_ROOT = Path(__file__).parent.parent.parent
PYDESEQ_DIR = PROJECT_ROOT / "examples/output_module5_pydeseq"
OUTPUT_DIR = Path(__file__).parent / "output"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Colors for responder vs progressor
RESPONDER_COLOR = '#2171b5'  # Blue
PROGRESSOR_COLOR = '#e74c3c'  # Red
NONSIG_COLOR = '#d4d4d4'  # Light gray for non-significant

# Significance thresholds
PADJ_THRESHOLD = 0.05
LOG2FC_THRESHOLD = 1.0  # For labeling top genes


# ---------------------------------------------------------------------------
# Data loading helpers
# ---------------------------------------------------------------------------

def load_de_results():
    """Load differential expression results from PyDESeq2."""
    path = PYDESEQ_DIR / "pseudobulk_de_results.csv"
    if path.exists():
        df = pd.read_csv(path, index_col=0)
        return df
    return None


def load_pathway_results():
    """Load pathway enrichment results for both directions."""
    resp_path = PYDESEQ_DIR / "pseudobulk_responder_up_MSigDB_Hallmark_2020.csv"
    prog_path = PYDESEQ_DIR / "pseudobulk_progressor_up_MSigDB_Hallmark_2020.csv"

    resp_df = pd.read_csv(resp_path) if resp_path.exists() else None
    prog_df = pd.read_csv(prog_path) if prog_path.exists() else None

    return resp_df, prog_df


def load_summary():
    """Load PyDESeq2 summary JSON."""
    path = PYDESEQ_DIR / "module5_pydeseq_summary.json"
    if path.exists():
        with open(path, 'r') as f:
            return json.load(f)
    return None


# ---------------------------------------------------------------------------
# Panel A: Volcano Plot
# ---------------------------------------------------------------------------

def panel_a_volcano(ax, de_df=None):
    """
    Volcano plot showing differential expression results.

    - X-axis: log2 fold change (positive = responder-up, negative = progressor-up)
    - Y-axis: -log10(adjusted p-value)
    - Colors: Blue (responder-up), Red (progressor-up), Gray (non-significant)
    - Labels: Top genes by significance and fold change
    """
    if de_df is None:
        de_df = load_de_results()

    if de_df is None:
        ax.text(0.5, 0.5, 'DE results not available', ha='center', va='center',
                transform=ax.transAxes)
        ax.set_title('A')
        return

    # Calculate -log10(padj)
    de_df = de_df.copy()
    de_df['neg_log10_padj'] = -np.log10(de_df['padj'].clip(lower=1e-100))

    # Classify genes
    de_df['category'] = 'nonsig'
    de_df.loc[(de_df['padj'] < PADJ_THRESHOLD) & (de_df['log2FoldChange'] > 0), 'category'] = 'responder_up'
    de_df.loc[(de_df['padj'] < PADJ_THRESHOLD) & (de_df['log2FoldChange'] < 0), 'category'] = 'progressor_up'

    # Plot non-significant first (background)
    nonsig = de_df[de_df['category'] == 'nonsig']
    ax.scatter(nonsig['log2FoldChange'], nonsig['neg_log10_padj'],
               c=NONSIG_COLOR, s=8, alpha=0.5, rasterized=True)

    # Plot responder-up
    resp = de_df[de_df['category'] == 'responder_up']
    ax.scatter(resp['log2FoldChange'], resp['neg_log10_padj'],
               c=RESPONDER_COLOR, s=15, alpha=0.7, label=f'Responder-up ({len(resp)})')

    # Plot progressor-up
    prog = de_df[de_df['category'] == 'progressor_up']
    ax.scatter(prog['log2FoldChange'], prog['neg_log10_padj'],
               c=PROGRESSOR_COLOR, s=15, alpha=0.7, label=f'Progressor-up ({len(prog)})')

    # Label top genes from each direction with manual offsets to avoid overlap
    # Separate highly significant genes (top of plot) from the clustered genes
    # Top responder-up genes (focus on most significant and well-separated)
    resp_gene_offsets = {
        'KLB': (5, 0),          # Most significant, far right, isolated
        'NEDD9': (5, 3),        # Second most significant
        'TMC5': (-40, 3),       # Third, offset left
        'GP2': (5, -8),         # Fifth, offset down
        'ZNF655': (-45, 0),     # Fourth, offset left
        'DACH1': (5, 5),        # Mentioned in legend
    }
    for gene, offset in resp_gene_offsets.items():
        if gene in de_df.index:
            row = de_df.loc[gene]
            if row['category'] == 'responder_up':
                ax.annotate(gene, (row['log2FoldChange'], row['neg_log10_padj']),
                           fontsize=7, ha='left' if offset[0] > 0 else 'right', va='center',
                           xytext=offset, textcoords='offset points')

    # Top progressor-up genes (manual offsets for each)
    prog_gene_offsets = {
        'VPREB3': (-5, 5),     # Most significant
        'COL4A4': (-5, -8),    # Second most significant, offset down
        'PLK1': (5, 3),        # Key proliferation gene
        'RAD51AP1': (-5, 0),   # DNA repair
        'CCNE1': (-5, 10),     # Cell cycle
        'FAM111B': (5, -3),    # Offset right
    }
    for gene, offset in prog_gene_offsets.items():
        if gene in de_df.index:
            row = de_df.loc[gene]
            if row['category'] == 'progressor_up':
                ax.annotate(gene, (row['log2FoldChange'], row['neg_log10_padj']),
                           fontsize=7, ha='right' if offset[0] < 0 else 'left', va='center',
                           xytext=offset, textcoords='offset points')

    # Add significance threshold line
    ax.axhline(-np.log10(PADJ_THRESHOLD), color='gray', linestyle='--', linewidth=0.8, alpha=0.5)
    ax.axvline(0, color='gray', linestyle='-', linewidth=0.5, alpha=0.3)

    # Styling
    ax.set_xlabel('log$_2$(Fold Change)')
    ax.set_ylabel('-log$_{10}$(adjusted p-value)')
    ax.set_title('A', fontweight='bold', loc='left')

    # Set axis limits to show data symmetrically
    xlim = max(abs(de_df['log2FoldChange'].min()), abs(de_df['log2FoldChange'].max()))
    ax.set_xlim(-xlim * 1.05, xlim * 1.05)

    # Legend
    ax.legend(loc='upper right', framealpha=0.9)

    # Add annotation for total genes
    n_sig = len(resp) + len(prog)
    ax.text(0.02, 0.98, f'203 significant genes\n(p$_{{adj}}$ < 0.05)',
            transform=ax.transAxes, fontsize=8, va='top', ha='left',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))


# ---------------------------------------------------------------------------
# Panel B: Pathway Enrichment Dot Plot
# ---------------------------------------------------------------------------

def panel_b_pathway_dotplot(ax, resp_df=None, prog_df=None):
    """
    Pathway enrichment dot plot split by direction.

    - Left side: Responder-up pathways (EMT, TNF-alpha, Myogenesis)
    - Right side: Progressor-up pathways (E2F Targets, G2-M Checkpoint)
    - Dot size: Gene overlap count
    - Dot color: -log10(adjusted p-value)
    """
    if resp_df is None or prog_df is None:
        resp_df, prog_df = load_pathway_results()

    if resp_df is None and prog_df is None:
        ax.text(0.5, 0.5, 'Pathway results not available', ha='center', va='center',
                transform=ax.transAxes)
        ax.set_title('B')
        return

    # Prepare data: combine both dataframes with direction label
    combined = []

    if resp_df is not None and len(resp_df) > 0:
        resp_df = resp_df.copy()
        resp_df['direction'] = 'Responder-up'
        resp_df['x_pos'] = -1  # Left side
        combined.append(resp_df)

    if prog_df is not None and len(prog_df) > 0:
        prog_df = prog_df.copy()
        prog_df['direction'] = 'Progressor-up'
        prog_df['x_pos'] = 1  # Right side
        combined.append(prog_df)

    if not combined:
        ax.text(0.5, 0.5, 'No significant pathways', ha='center', va='center',
                transform=ax.transAxes)
        return

    df = pd.concat(combined, ignore_index=True)

    # Parse overlap to get gene count
    df['n_genes'] = df['Overlap'].apply(lambda x: int(x.split('/')[0]))
    df['neg_log10_padj'] = -np.log10(df['Adjusted P-value'])

    # Clean up term names
    df['Term_clean'] = df['Term'].str.replace('_', ' ')

    # Sort by -log10 p-value within each direction
    df = df.sort_values(['direction', 'neg_log10_padj'], ascending=[True, False])

    # Assign y positions (0, 1, 2, ... for each direction)
    resp_data = df[df['direction'] == 'Responder-up'].reset_index(drop=True)
    prog_data = df[df['direction'] == 'Progressor-up'].reset_index(drop=True)

    # Plot responder-up pathways (left side, blue)
    for i, row in resp_data.iterrows():
        size = row['n_genes'] * 30  # Scale factor for dot size
        color_val = row['neg_log10_padj']
        ax.scatter(-0.5, i, s=size, c=RESPONDER_COLOR, alpha=0.8,
                   edgecolors='white', linewidths=0.5)
        ax.text(-1.0, i, row['Term_clean'], ha='right', va='center', fontsize=9)

    # Plot progressor-up pathways (right side, red)
    for i, row in prog_data.iterrows():
        size = row['n_genes'] * 30
        ax.scatter(0.5, i, s=size, c=PROGRESSOR_COLOR, alpha=0.8,
                   edgecolors='white', linewidths=0.5)
        ax.text(1.0, i, row['Term_clean'], ha='left', va='center', fontsize=9)

    # Add direction labels
    max_y = max(len(resp_data), len(prog_data)) - 1
    ax.text(-0.5, max_y + 0.8, 'Responder-up', ha='center', va='bottom',
            fontsize=10, fontweight='bold', color=RESPONDER_COLOR)
    ax.text(0.5, max_y + 0.8, 'Progressor-up', ha='center', va='bottom',
            fontsize=10, fontweight='bold', color=PROGRESSOR_COLOR)

    # Add center divider line
    ax.axvline(0, color='gray', linestyle='-', linewidth=0.5, alpha=0.3)

    # Add size legend
    sizes_legend = [3, 5, 7]
    for i, s in enumerate(sizes_legend):
        ax.scatter([], [], s=s*30, c='gray', alpha=0.7,
                   label=f'{s} genes', edgecolors='white')
    legend = ax.legend(title='Gene overlap', loc='lower right',
                       framealpha=0.9, fontsize=8)
    legend.get_title().set_fontsize(9)

    # Styling
    ax.set_xlim(-2.5, 2.5)
    ax.set_ylim(-0.5, max_y + 1.5)
    ax.set_title('B', fontweight='bold', loc='left')
    ax.set_xlabel('')
    ax.set_ylabel('')

    # Remove axes for cleaner look
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.set_xticks([])
    ax.set_yticks([])


# ---------------------------------------------------------------------------
# Main figure generation
# ---------------------------------------------------------------------------

def generate_supp_figure2():
    """Generate the complete Supplementary Figure S2."""
    apply_style()

    # Load data
    de_df = load_de_results()
    resp_df, prog_df = load_pathway_results()

    # Create figure with 1x2 layout
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Generate panels
    panel_a_volcano(axes[0], de_df)
    panel_b_pathway_dotplot(axes[1], resp_df, prog_df)

    # Adjust layout
    plt.tight_layout()

    # Save in multiple formats
    for fmt in ['png', 'pdf', 'svg']:
        outpath = OUTPUT_DIR / f"supp_figure2_de_pathway.{fmt}"
        fig.savefig(outpath, format=fmt, dpi=300 if fmt == 'png' else None,
                    facecolor='white', bbox_inches='tight')
        print(f"  Saved: {outpath}")

    plt.close(fig)

    return fig


if __name__ == "__main__":
    print("=" * 60)
    print("Generating Supplementary Figure S2: DE & Pathway Enrichment")
    print("=" * 60)
    generate_supp_figure2()
    print("\nDone.")
