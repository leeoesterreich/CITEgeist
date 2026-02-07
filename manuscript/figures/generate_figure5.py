#!/usr/bin/env python3
"""
Figure 5: Module 5 Cross-Sample Integration

Panel A: SCHEMATIC - use output/schematics/figure5_panel_a_integration.svg
Panels B, C, D: DATA - generated with matplotlib below
"""

import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
from pathlib import Path
from adjustText import adjust_text

from figure_style import apply_style, PALETTE, get_cell_type_color

apply_style()

try:
    from sklearn.decomposition import PCA
    HAS_SKLEARN = True
except ImportError:
    HAS_SKLEARN = False

try:
    import umap
    HAS_UMAP = True
except ImportError:
    HAS_UMAP = False

# Paths
PROJECT_ROOT = Path(__file__).parent.parent.parent
MODULE5_DIR = PROJECT_ROOT / "output/module5_integration"
OUTPUT_DIR = Path(__file__).parent / "output"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Response colors
RESPONSE_COLORS = {
    'responder': '#2ecc71',
    'progressor': '#e74c3c',
}


def load_data():
    aligned = pd.read_csv(MODULE5_DIR / "module5_unified_aligned_programs.csv")
    embedding = np.load(MODULE5_DIR / "module5_unified_embedding.npy")
    metadata = pd.read_csv(MODULE5_DIR / "module5_unified_program_metadata.csv")
    with open(MODULE5_DIR / "module5_response_analysis.json", 'r') as f:
        response = json.load(f)
    return aligned, embedding, metadata, response


def panel_b_umap(ax, embedding, metadata):
    """Panel B: UMAP of programs."""
    ax.text(-0.12, 1.05, "B", fontsize=14, fontweight='bold', va='top', transform=ax.transAxes)

    n = min(embedding.shape[0], len(metadata))
    embedding = embedding[:n]
    metadata = metadata.iloc[:n]

    if HAS_UMAP and n > 15:
        try:
            reducer = umap.UMAP(n_neighbors=15, min_dist=0.1, random_state=42)
            coords = reducer.fit_transform(embedding)
        except:
            coords = embedding[:, :2] if embedding.shape[1] >= 2 else np.random.randn(n, 2)
    elif HAS_SKLEARN:
        pca = PCA(n_components=2, random_state=42)
        coords = pca.fit_transform(embedding)
    else:
        coords = embedding[:, :2]

    colors = [get_cell_type_color(ct) for ct in metadata['cell_type']]
    ax.scatter(coords[:, 0], coords[:, 1], c=colors, alpha=0.7, s=30, edgecolors='white', linewidths=0.3)

    ax.set_xlabel("UMAP 1", fontsize=10)
    ax.set_ylabel("UMAP 2", fontsize=10)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    present_ct = list(set(metadata['cell_type']))[:6]
    handles = [mpatches.Patch(color=get_cell_type_color(ct),
               label=ct.replace('_', ' ')[:12]) for ct in sorted(present_ct)]
    ax.legend(handles=handles, loc='upper right', fontsize=9, framealpha=0.9, ncol=1)

    ax.text(0.02, 0.02, f'n={len(metadata)} programs', transform=ax.transAxes, fontsize=9,
            color=PALETTE['neutral'])


def panel_c_response(ax, response_data):
    """Panel C: Response-associated programs bar chart."""
    ax.text(-0.12, 1.05, "C", fontsize=14, fontweight='bold', va='top', transform=ax.transAxes)

    resp = len(response_data.get('responder_enriched', []))
    prog = len(response_data.get('progressor_enriched', []))

    bars = ax.bar(['Responder\nenriched', 'Progressor\nenriched'], [resp, prog],
                  color=[RESPONSE_COLORS['responder'], RESPONSE_COLORS['progressor']], alpha=0.8, width=0.5)

    for bar, val in zip(bars, [resp, prog]):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.3,
                str(val), ha='center', va='bottom', fontsize=11, fontweight='bold')

    ax.set_ylabel("Number of Programs", fontsize=10)
    ax.set_ylim(0, max(resp, prog) * 1.35)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)


def panel_d_volcano(ax):
    """Panel D: Volcano plot."""
    ax.text(-0.12, 1.05, "D", fontsize=14, fontweight='bold', va='top', transform=ax.transAxes)

    deseq_file = PROJECT_ROOT / "examples/output_module5_pydeseq/pseudobulk_de_results.csv"
    if not deseq_file.exists():
        ax.text(0.5, 0.5, "DE results\nnot available", ha='center', va='center', fontsize=10,
                style='italic', color=PALETTE['neutral'])
        ax.axis('off')
        return

    de = pd.read_csv(deseq_file, index_col=0).dropna(subset=['padj', 'log2FoldChange'])
    de['neg_log10_padj'] = -np.log10(de['padj'].clip(lower=1e-50))

    de['cat'] = 'ns'
    de.loc[(de['padj'] < 0.05) & (de['log2FoldChange'] > 0), 'cat'] = 'resp'
    de.loc[(de['padj'] < 0.05) & (de['log2FoldChange'] < 0), 'cat'] = 'prog'

    colors = {'ns': PALETTE['border'], 'resp': RESPONSE_COLORS['responder'], 'prog': RESPONSE_COLORS['progressor']}

    for cat in ['ns', 'prog', 'resp']:
        sub = de[de['cat'] == cat]
        ax.scatter(sub['log2FoldChange'], sub['neg_log10_padj'],
                   c=colors[cat], alpha=0.5 if cat == 'ns' else 0.7,
                   s=10 if cat == 'ns' else 25, label=cat.title(), edgecolors='none')

    ax.axhline(-np.log10(0.05), color=PALETTE['neutral'], linestyle='--', lw=0.8, alpha=0.7)

    # Label top genes with adjustText
    texts = []
    for subset, n in [('prog', 5), ('resp', 3)]:
        top = de[de['cat'] == subset].nsmallest(n, 'padj')
        for _, row in top.iterrows():
            t = ax.text(row['log2FoldChange'], row['neg_log10_padj'], row.name,
                       fontsize=9, ha='center', va='bottom')
            texts.append(t)

    adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle='-', color='gray', lw=0.5),
                expand_points=(1.5, 1.5), force_text=(0.8, 0.8))

    ax.set_xlabel("log2 Fold Change", fontsize=10)
    ax.set_ylabel("-log10(p-adj)", fontsize=10)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    n_resp = (de['cat'] == 'resp').sum()
    n_prog = (de['cat'] == 'prog').sum()
    ax.text(0.98, 0.98, f"{n_resp} resp-enriched\n{n_prog} prog-enriched",
            transform=ax.transAxes, ha='right', va='top', fontsize=9,
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))


def generate_figure5():
    """Generate Figure 5 with data panels."""
    print("Loading data...")
    aligned, embedding, metadata, response = load_data()
    print(f"Programs: {len(aligned)}, Embedding: {embedding.shape}")

    fig = plt.figure(figsize=(11, 8))
    gs = GridSpec(2, 2, figure=fig, hspace=0.30, wspace=0.28)

    # Panel A: Placeholder for schematic
    ax_a = fig.add_subplot(gs[0, 0])
    ax_a.text(0.5, 0.5, "Panel A: Schematic\n\nUse SVG file:\nfigure5_panel_a_integration.svg",
              ha='center', va='center', fontsize=11, style='italic',
              bbox=dict(boxstyle='round', facecolor='#f0f0f0', edgecolor='gray'))
    ax_a.set_xlim(0, 1)
    ax_a.set_ylim(0, 1)
    ax_a.axis('off')
    ax_a.set_title("Cross-Sample Integration", fontsize=12, fontweight='bold', loc='left', pad=10)

    # Panel B: UMAP (DATA)
    ax_b = fig.add_subplot(gs[0, 1])
    panel_b_umap(ax_b, embedding, metadata)
    ax_b.set_title("Program Embedding", fontsize=12, fontweight='bold', loc='left', pad=10)

    # Panel C: Response (DATA)
    ax_c = fig.add_subplot(gs[1, 0])
    panel_c_response(ax_c, response)
    ax_c.set_title("Response-Associated Programs", fontsize=12, fontweight='bold', loc='left', pad=10)

    # Panel D: Volcano (DATA)
    ax_d = fig.add_subplot(gs[1, 1])
    panel_d_volcano(ax_d)
    ax_d.set_title("Differential Expression", fontsize=12, fontweight='bold', loc='left', pad=10)

    plt.tight_layout()

    for fmt, dpi in [('pdf', 300), ('png', 150), ('svg', None)]:
        output_path = OUTPUT_DIR / f"figure5_integration.{fmt}"
        if fmt == 'svg':
            plt.savefig(output_path, format='svg', bbox_inches='tight', facecolor='white')
        else:
            plt.savefig(output_path, dpi=dpi, bbox_inches='tight', facecolor='white')
        print(f"Saved: {output_path}")

    plt.close()

    print("\n" + "=" * 60)
    print("Figure 5: Cross-Sample Integration")
    print("=" * 60)
    print("\nPanel A: SCHEMATIC - use output/schematics/figure5_panel_a_integration.svg")
    print("Panels B, C, D: DATA - generated above")
    print(f"\nAligned programs: {len(aligned)}")
    print(f"Responder-enriched: {len(response.get('responder_enriched', []))}")
    print(f"Progressor-enriched: {len(response.get('progressor_enriched', []))}")


if __name__ == "__main__":
    generate_figure5()
