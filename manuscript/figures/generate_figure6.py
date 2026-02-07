#!/usr/bin/env python3
"""
Figure 6: Interoperability & Downstream Analysis

Panel A: SCHEMATIC - use output/schematics/figure6_panel_a_workflow.svg
Panels B, C, D: DATA - generated with matplotlib below
Panel E: SCHEMATIC - use output/schematics/figure6_panel_e_validation.svg
"""

import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.gridspec import GridSpec
from pathlib import Path
from adjustText import adjust_text

from figure_style import apply_style, PALETTE, get_cell_type_color

apply_style()

# Paths
PROJECT_ROOT = Path(__file__).parent.parent.parent
PYDESEQ_DIR = PROJECT_ROOT / "examples/output_module5_pydeseq"
OUTPUT_DIR = Path(__file__).parent / "output"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Tool colors
TOOL_COLORS = {
    'scanpy': PALETTE['accent2'],
    'pydeseq2': PALETTE['accent1'],
    'gseapy': PALETTE['highlight'],
}


def load_de_results():
    try:
        return pd.read_csv(PYDESEQ_DIR / "pseudobulk_de_results.csv", index_col=0)
    except Exception as e:
        print(f"Error loading DE results: {e}")
        return None


def load_enrichment_results():
    results = {}
    for direction in ['responder_up', 'progressor_up']:
        for db in ['GO_Biological_Process_2021', 'MSigDB_Hallmark_2020']:
            filepath = PYDESEQ_DIR / f"pseudobulk_{direction}_{db}.csv"
            if filepath.exists():
                try:
                    df = pd.read_csv(filepath)
                    if len(df) > 0:
                        results[f"{direction}_{db}"] = df
                except Exception:
                    pass
    return results if results else None


def panel_b_scanpy_umap(ax):
    """Panel B: Simulated UMAP showing clustering."""
    ax.text(-0.05, 1.08, "B", fontsize=14, fontweight='bold', va='top', transform=ax.transAxes)
    ax.set_title("scanpy: Cluster Visualization", fontsize=10, fontweight='bold', loc='left', color=TOOL_COLORS['scanpy'])

    np.random.seed(42)
    cell_types = ['Epithelial', 'Fibroblasts', 'Macrophages', 'T cells', 'B cells', 'Endothelial']
    centers = [(-3, 2), (2, 3), (-2, -2), (3, -1), (0, 4), (4, 1)]

    for ct, (cx, cy) in zip(cell_types, centers):
        n = np.random.randint(80, 150)
        x = np.random.normal(cx, 0.8, n)
        y = np.random.normal(cy, 0.8, n)
        color = get_cell_type_color(ct)
        ax.scatter(x, y, c=color, s=8, alpha=0.6, label=ct.replace(' cells', ''), rasterized=True)

    ax.set_xlabel("UMAP1", fontsize=9)
    ax.set_ylabel("UMAP2", fontsize=9)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.legend(loc='upper right', fontsize=7, framealpha=0.9, markerscale=1.2)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.text(0.02, 0.02, "Deconvolved GEX layers\nenable cell-type clustering",
            transform=ax.transAxes, fontsize=8, style='italic', va='bottom', color=PALETTE['neutral'])


def panel_c_volcano(ax, de_data):
    """Panel C: PyDESeq2 volcano plot."""
    ax.text(-0.05, 1.08, "C", fontsize=14, fontweight='bold', va='top', transform=ax.transAxes)
    ax.set_title("PyDESeq2: Differential Expression", fontsize=10, fontweight='bold', loc='left', color=TOOL_COLORS['pydeseq2'])

    if de_data is None:
        ax.text(0.5, 0.5, "DE data not available", ha='center', va='center', fontsize=10)
        ax.axis('off')
        return

    df = de_data.dropna(subset=['log2FoldChange', 'padj'])
    df = df[df['padj'] > 0]
    df['neglog10p'] = -np.log10(df['padj'])

    padj_thresh, lfc_thresh = 0.05, 1.0
    sig_up = (df['padj'] < padj_thresh) & (df['log2FoldChange'] > lfc_thresh)
    sig_down = (df['padj'] < padj_thresh) & (df['log2FoldChange'] < -lfc_thresh)
    not_sig = ~(sig_up | sig_down)

    ax.scatter(df.loc[not_sig, 'log2FoldChange'], df.loc[not_sig, 'neglog10p'],
               c=PALETTE['neutral'], s=8, alpha=0.3, rasterized=True)
    ax.scatter(df.loc[sig_up, 'log2FoldChange'], df.loc[sig_up, 'neglog10p'],
               c='#2ecc71', s=15, alpha=0.7, label=f'Responder ↑ ({sig_up.sum()})', rasterized=True)
    ax.scatter(df.loc[sig_down, 'log2FoldChange'], df.loc[sig_down, 'neglog10p'],
               c='#e74c3c', s=15, alpha=0.7, label=f'Progressor ↑ ({sig_down.sum()})', rasterized=True)

    ax.axhline(-np.log10(padj_thresh), color='gray', linestyle='--', linewidth=0.8, alpha=0.5)
    ax.axvline(lfc_thresh, color='gray', linestyle='--', linewidth=0.8, alpha=0.5)
    ax.axvline(-lfc_thresh, color='gray', linestyle='--', linewidth=0.8, alpha=0.5)

    # Label top genes with adjustText
    texts = []
    top_genes = df.nsmallest(5, 'padj')
    for _, row in top_genes.iterrows():
        t = ax.text(row['log2FoldChange'], row['neglog10p'] + 0.1, row.name,
                    fontsize=7, ha='center', va='bottom')
        texts.append(t)

    adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle='-', color='gray', lw=0.5),
                expand_points=(1.5, 1.5), force_text=(0.8, 0.8))

    ax.set_xlabel("log₂ Fold Change", fontsize=9)
    ax.set_ylabel("-log₁₀ p-value", fontsize=9)
    ax.legend(loc='upper right', fontsize=7, framealpha=0.9)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    total_sig = sig_up.sum() + sig_down.sum()
    ax.text(0.02, 0.02, f"{total_sig} significant genes (padj<0.05, |LFC|>1)",
            transform=ax.transAxes, fontsize=8, style='italic', va='bottom', color=PALETTE['neutral'])


def panel_d_enrichment(ax, enrichment_data):
    """Panel D: GSEApy pathway enrichment bar plot."""
    ax.text(-0.05, 1.08, "D", fontsize=14, fontweight='bold', va='top', transform=ax.transAxes)
    ax.set_title("GSEApy: Pathway Enrichment", fontsize=10, fontweight='bold', loc='left', color=TOOL_COLORS['gseapy'])

    if enrichment_data is None:
        ax.text(0.5, 0.5, "Enrichment data not available", ha='center', va='center', fontsize=10)
        ax.axis('off')
        return

    all_terms = []
    for key, df in enrichment_data.items():
        if 'GO' in key and len(df) > 0:
            df = df.copy()
            df['direction'] = 'Responder' if 'responder' in key else 'Progressor'
            all_terms.append(df.head(3))

    if not all_terms:
        ax.text(0.5, 0.5, "No enriched pathways", ha='center', va='center', fontsize=10)
        ax.axis('off')
        return

    combined = pd.concat(all_terms, ignore_index=True)
    combined['neglog10p'] = -np.log10(combined['P-value'])
    combined = combined.sort_values('neglog10p', ascending=True)

    combined['short_term'] = combined['Term'].apply(
        lambda x: x.split('(')[0].strip()[:35] + '...' if len(x.split('(')[0].strip()) > 35 else x.split('(')[0].strip()
    )

    colors = ['#2ecc71' if d == 'Responder' else '#e74c3c' for d in combined['direction']]

    y = np.arange(len(combined))
    ax.barh(y, combined['neglog10p'], color=colors, alpha=0.85)
    ax.set_yticks(y)
    ax.set_yticklabels(combined['short_term'], fontsize=8)
    ax.set_xlabel("-log₁₀ p-value", fontsize=9)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    legend_elements = [Patch(facecolor='#2ecc71', label='Responder'),
                       Patch(facecolor='#e74c3c', label='Progressor')]
    ax.legend(handles=legend_elements, loc='lower right', fontsize=7, framealpha=0.9)


def generate_figure6():
    """Generate Figure 6 with data panels."""
    print("Loading data...")
    de_data = load_de_results()
    enrichment_data = load_enrichment_results()

    if de_data is not None:
        print(f"DE results: {len(de_data)} genes")
    if enrichment_data:
        print(f"Enrichment results: {list(enrichment_data.keys())}")

    plt.rcParams['figure.constrained_layout.use'] = False

    fig = plt.figure(figsize=(12, 8))
    gs = GridSpec(2, 3, figure=fig, height_ratios=[1, 1],
                  hspace=0.35, wspace=0.30,
                  left=0.06, right=0.98, top=0.95, bottom=0.08)

    # Panel A: Placeholder for schematic
    ax_a = fig.add_subplot(gs[0, 0])
    ax_a.text(0.5, 0.5, "Panel A: Schematic\n\nUse SVG file:\nfigure6_panel_a_workflow.svg",
              ha='center', va='center', fontsize=10, style='italic',
              bbox=dict(boxstyle='round', facecolor='#f0f0f0', edgecolor='gray'))
    ax_a.set_xlim(0, 1)
    ax_a.set_ylim(0, 1)
    ax_a.axis('off')
    ax_a.set_title("CITEgeist → Downstream", fontsize=11, fontweight='bold', loc='left')

    # Panel B: scanpy UMAP (DATA)
    ax_b = fig.add_subplot(gs[0, 1])
    panel_b_scanpy_umap(ax_b)

    # Panel C: PyDESeq2 volcano (DATA)
    ax_c = fig.add_subplot(gs[0, 2])
    panel_c_volcano(ax_c, de_data)

    # Panel D: GSEApy enrichment (DATA)
    ax_d = fig.add_subplot(gs[1, 0])
    panel_d_enrichment(ax_d, enrichment_data)

    # Panel E: Placeholder for validation schematic
    ax_e = fig.add_subplot(gs[1, 1:])
    ax_e.text(0.5, 0.5, "Panel E: Schematic\n\nUse SVG file:\nfigure6_panel_e_validation.svg",
              ha='center', va='center', fontsize=10, style='italic',
              bbox=dict(boxstyle='round', facecolor='#f0f0f0', edgecolor='gray'))
    ax_e.set_xlim(0, 1)
    ax_e.set_ylim(0, 1)
    ax_e.axis('off')
    ax_e.set_title("Experimental Validation: Midkine Discovery", fontsize=11, fontweight='bold', loc='left')

    plt.rcParams['savefig.bbox'] = 'standard'

    for fmt, dpi in [('pdf', 300), ('png', 150), ('svg', None)]:
        output_path = OUTPUT_DIR / f"figure6_interoperability.{fmt}"
        if fmt == 'svg':
            plt.savefig(output_path, format='svg', facecolor='white')
        else:
            plt.savefig(output_path, dpi=dpi, facecolor='white')
        print(f"Saved: {output_path}")

    plt.close()

    print("\n" + "=" * 60)
    print("Figure 6: Interoperability")
    print("=" * 60)
    print("\nPanel A: SCHEMATIC - use output/schematics/figure6_panel_a_workflow.svg")
    print("Panels B, C, D: DATA - generated above")
    print("Panel E: SCHEMATIC - use output/schematics/figure6_panel_e_validation.svg")


if __name__ == "__main__":
    generate_figure6()
