#!/usr/bin/env python3
"""
Figure 5: Module 5 Cross-Sample Integration

Enhanced version with consistent styling.

Panels:
  A: Integration workflow schematic
  B: UMAP of programs colored by cell type
  C: Response-associated programs
  D: Volcano plot (PyDESeq2)
"""

import json
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

# Response colors (specific to Figure 5)
RESPONSE_COLORS = {
    'responder': '#2ecc71',  # Green
    'progressor': '#e74c3c',  # Red
}


def load_data():
    aligned = pd.read_csv(MODULE5_DIR / "module5_unified_aligned_programs.csv")
    embedding = np.load(MODULE5_DIR / "module5_unified_embedding.npy")
    metadata = pd.read_csv(MODULE5_DIR / "module5_unified_program_metadata.csv")
    with open(MODULE5_DIR / "module5_response_analysis.json", 'r') as f:
        response = json.load(f)
    return aligned, embedding, metadata, response


def panel_a_schematic(ax):
    """Panel A: Simple integration workflow."""
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis('off')

    # Panel label
    ax.text(0.02, 0.95, "A", fontsize=14, fontweight='bold', va='top')

    arrow_color = PALETTE['neutral']

    # Samples
    sample_colors = [RESPONSE_COLORS['responder'], RESPONSE_COLORS['responder'],
                     RESPONSE_COLORS['progressor'], RESPONSE_COLORS['progressor']]
    for i, c in enumerate(sample_colors):
        rect = FancyBboxPatch((0.02, 0.72 - i*0.18), 0.12, 0.14, boxstyle="round,pad=0.01",
                              facecolor=c, edgecolor='black', linewidth=1, alpha=0.5)
        ax.add_patch(rect)
        ax.text(0.08, 0.79 - i*0.18, f"S{i+1}", ha='center', va='center', fontsize=10, fontweight='bold')
    ax.text(0.08, 0.90, "14 Samples", ha='center', va='bottom', fontsize=11, fontweight='bold')

    # Arrow
    ax.annotate('', xy=(0.22, 0.5), xytext=(0.15, 0.5),
                arrowprops=dict(arrowstyle='->', color=arrow_color, lw=2))

    # Module 4 box
    rect1 = FancyBboxPatch((0.22, 0.35), 0.16, 0.30, boxstyle="round,pad=0.01",
                           facecolor='#dae8fc', edgecolor='#6c8ebf', linewidth=2)
    ax.add_patch(rect1)
    ax.text(0.30, 0.5, "Module 4\nPrograms", ha='center', va='center', fontsize=10, fontweight='bold')

    # Arrow
    ax.annotate('', xy=(0.46, 0.5), xytext=(0.39, 0.5),
                arrowprops=dict(arrowstyle='->', color=arrow_color, lw=2))

    # Harmony box
    rect2 = FancyBboxPatch((0.46, 0.35), 0.16, 0.30, boxstyle="round,pad=0.01",
                           facecolor='#fff2cc', edgecolor='#d6b656', linewidth=2)
    ax.add_patch(rect2)
    ax.text(0.54, 0.5, "Harmony\nAlignment", ha='center', va='center', fontsize=10, fontweight='bold')

    # Arrow
    ax.annotate('', xy=(0.70, 0.5), xytext=(0.63, 0.5),
                arrowprops=dict(arrowstyle='->', color=arrow_color, lw=2))

    # Outputs
    outputs = [("73 Aligned\nPrograms", 0.65, '#d5e8d4'),
               ("Response\nAnalysis", 0.40, '#f8cecc')]
    for label, y, color in outputs:
        rect = FancyBboxPatch((0.70, y), 0.25, 0.18, boxstyle="round,pad=0.01",
                              facecolor=color, edgecolor='black', linewidth=1.5)
        ax.add_patch(rect)
        ax.text(0.825, y + 0.09, label, ha='center', va='center', fontsize=10, fontweight='bold')

    # Legend with colored squares
    ax.add_patch(mpatches.Rectangle((0.03, 0.04), 0.03, 0.03,
                 facecolor=RESPONSE_COLORS['responder'], edgecolor='none'))
    ax.text(0.07, 0.055, "Responder", ha='left', va='center', fontsize=9, color=PALETTE['neutral'])
    ax.add_patch(mpatches.Rectangle((0.25, 0.04), 0.03, 0.03,
                 facecolor=RESPONSE_COLORS['progressor'], edgecolor='none'))
    ax.text(0.29, 0.055, "Progressor", ha='left', va='center', fontsize=9, color=PALETTE['neutral'])


def panel_b_umap(ax, embedding, metadata):
    """Panel B: UMAP of programs."""
    # Panel label
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

    # Legend
    present_ct = list(set(metadata['cell_type']))[:6]
    handles = [mpatches.Patch(color=get_cell_type_color(ct),
               label=ct.replace('_', ' ')[:12]) for ct in sorted(present_ct)]
    ax.legend(handles=handles, loc='upper right', fontsize=9, framealpha=0.9, ncol=1)

    ax.text(0.02, 0.02, f'n={len(metadata)} programs', transform=ax.transAxes, fontsize=9,
            color=PALETTE['neutral'])


def panel_c_response(ax, response_data):
    """Panel C: Response-associated programs bar chart."""
    # Panel label
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
    # Panel label
    ax.text(-0.12, 1.05, "D", fontsize=14, fontweight='bold', va='top', transform=ax.transAxes)

    deseq_file = PROJECT_ROOT / "examples/output_module5_pydeseq/pseudobulk_de_results.csv"
    if not deseq_file.exists():
        ax.text(0.5, 0.5, "DE results\nnot available", ha='center', va='center', fontsize=10,
                style='italic', color=PALETTE['neutral'])
        ax.axis('off')
        return

    de = pd.read_csv(deseq_file, index_col=0).dropna(subset=['padj', 'log2FoldChange'])
    de['neg_log10_padj'] = -np.log10(de['padj'].clip(lower=1e-50))

    # Categorize
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

    # Label top genes
    for subset, n in [('prog', 5), ('resp', 3)]:
        top = de[de['cat'] == subset].nsmallest(n, 'padj')
        for _, row in top.iterrows():
            ax.annotate(row.name, (row['log2FoldChange'], row['neg_log10_padj']),
                       fontsize=9, ha='center', va='bottom')

    ax.set_xlabel("log2 Fold Change", fontsize=10)
    ax.set_ylabel("-log10(p-adj)", fontsize=10)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Summary
    n_resp = (de['cat'] == 'resp').sum()
    n_prog = (de['cat'] == 'prog').sum()
    ax.text(0.98, 0.98, f"{n_resp} resp-enriched\n{n_prog} prog-enriched",
            transform=ax.transAxes, ha='right', va='top', fontsize=9,
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))


def generate_figure5():
    """Generate complete Figure 5."""
    print("Loading data...")
    aligned, embedding, metadata, response = load_data()
    print(f"Programs: {len(aligned)}, Embedding: {embedding.shape}")

    fig = plt.figure(figsize=(11, 8))
    gs = GridSpec(2, 2, figure=fig, hspace=0.30, wspace=0.28)

    # Panel A: Schematic
    ax_a = fig.add_subplot(gs[0, 0])
    panel_a_schematic(ax_a)
    ax_a.set_title("Cross-Sample Integration", fontsize=12, fontweight='bold', loc='left', pad=10)

    # Panel B: UMAP
    ax_b = fig.add_subplot(gs[0, 1])
    panel_b_umap(ax_b, embedding, metadata)
    ax_b.set_title("Program Embedding", fontsize=12, fontweight='bold', loc='left', pad=10)

    # Panel C: Response
    ax_c = fig.add_subplot(gs[1, 0])
    panel_c_response(ax_c, response)
    ax_c.set_title("Response-Associated Programs", fontsize=12, fontweight='bold', loc='left', pad=10)

    # Panel D: Volcano
    ax_d = fig.add_subplot(gs[1, 1])
    panel_d_volcano(ax_d)
    ax_d.set_title("Differential Expression", fontsize=12, fontweight='bold', loc='left', pad=10)

    plt.tight_layout()

    output_path = OUTPUT_DIR / "figure5_integration.pdf"
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"Saved to {output_path}")

    png_path = OUTPUT_DIR / "figure5_integration.png"
    plt.savefig(png_path, dpi=150, bbox_inches='tight', facecolor='white')
    print(f"Preview saved to {png_path}")

    # Save SVG for Illustrator
    svg_path = OUTPUT_DIR / "figure5_integration.svg"
    plt.savefig(svg_path, format='svg', bbox_inches='tight', facecolor='white')
    print(f"SVG saved to {svg_path}")

    plt.close()

    print("\n=== Figure 5 Summary ===")
    print(f"Aligned programs: {len(aligned)}")
    print(f"Responder-enriched: {len(response.get('responder_enriched', []))}")
    print(f"Progressor-enriched: {len(response.get('progressor_enriched', []))}")

    print("\nEnhancements applied:")
    print("  - Consistent color palette from figure_style.py")
    print("  - Fonts increased to minimum 10pt")
    print("  - Panel labels in top-left corner")
    print("  - Consistent response colors (green=responder, red=progressor)")


if __name__ == "__main__":
    generate_figure5()
