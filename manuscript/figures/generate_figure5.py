#!/usr/bin/env python3
"""
Figure 5: Module 5 Cross-Sample Integration (Tour de Force #2)

Panels:
  A: 14-sample integration schematic (placeholder showing samples -> alignment -> unified space)
  B: UMAP/PCA of aligned programs colored by cell type and/or sample
  C: Responder vs Progressor enriched programs (bar chart showing 3 responder-enriched, 4 progressor-enriched)
  D: PyDESeq2 differential expression results (placeholder - will be filled when job completes)
  E: Conserved relationship network visualization

Data sources:
  - Module 5 outputs: output/module5_integration/
    - module5_unified_aligned_programs.csv (73 aligned programs with cell types, samples, top genes)
    - module5_unified_conserved_relationships.csv (191 relationships: co-localized, exclusive, independent)
    - module5_unified_embedding.npy (program embeddings for UMAP/PCA)
    - module5_unified_program_metadata.csv
    - module5_response_analysis.json (responder vs progressor program assignments)
    - module5_summary.json
"""

import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from pathlib import Path

try:
    import networkx as nx
    HAS_NETWORKX = True
except ImportError:
    HAS_NETWORKX = False
    print("Warning: networkx not available, network visualization will be simplified")

try:
    from sklearn.decomposition import PCA
    from sklearn.manifold import TSNE
    HAS_SKLEARN = True
except ImportError:
    HAS_SKLEARN = False
    print("Warning: sklearn not available, using raw embedding for visualization")

try:
    import umap
    HAS_UMAP = True
except ImportError:
    HAS_UMAP = False
    print("Warning: umap not available, will use PCA instead")

# Paths
PROJECT_ROOT = Path(__file__).parent.parent.parent
MODULE5_DIR = PROJECT_ROOT / "output/module5_integration"
OUTPUT_DIR = Path(__file__).parent / "output"

# Ensure output directory exists
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Style settings
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.size'] = 8
plt.rcParams['axes.linewidth'] = 0.5
plt.rcParams['xtick.major.width'] = 0.5
plt.rcParams['ytick.major.width'] = 0.5

# Color palette for cell types (consistent with patient data)
CELL_TYPE_COLORS = {
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

# Sample colors (for 14 samples, use a colormap)
SAMPLE_CMAP = plt.cm.tab20


def load_aligned_programs():
    """Load aligned programs from Module 5."""
    return pd.read_csv(MODULE5_DIR / "module5_unified_aligned_programs.csv")


def load_conserved_relationships():
    """Load conserved relationships from Module 5."""
    return pd.read_csv(MODULE5_DIR / "module5_unified_conserved_relationships.csv")


def load_response_analysis():
    """Load responder/progressor analysis."""
    with open(MODULE5_DIR / "module5_response_analysis.json", 'r') as f:
        return json.load(f)


def load_summary():
    """Load Module 5 summary."""
    with open(MODULE5_DIR / "module5_summary.json", 'r') as f:
        return json.load(f)


def load_embedding():
    """Load program embeddings."""
    return np.load(MODULE5_DIR / "module5_unified_embedding.npy")


def load_program_metadata():
    """Load program metadata for matching embedding rows to cell types."""
    return pd.read_csv(MODULE5_DIR / "module5_unified_program_metadata.csv")


def panel_a_schematic(ax):
    """Panel A: Integration schematic placeholder."""
    ax.text(0.5, 0.5, "Panel A: Cross-Sample Integration Schematic\n"
            "(Create in vector editor)\n\n"
            "Show:\n"
            "- 14 patient samples (P1-S1, P1-S2, ...)\n"
            "- Per-sample program discovery (Module 4)\n"
            "- Harmony alignment in latent space\n"
            "- Unified program identities across samples",
            ha='center', va='center', fontsize=9, style='italic',
            bbox=dict(boxstyle='round', facecolor='#f0f0f0', edgecolor='gray'))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis('off')
    ax.set_title("A. Module 5: Cross-Sample Integration", fontsize=10, fontweight='bold', loc='left')


def panel_b_embedding_visualization(ax, embedding, metadata):
    """Panel B: UMAP/PCA of aligned programs colored by cell type."""
    # Reduce to 2D for visualization
    if embedding.shape[0] != len(metadata):
        print(f"Warning: embedding ({embedding.shape[0]}) and metadata ({len(metadata)}) size mismatch")
        # Use minimum size
        n = min(embedding.shape[0], len(metadata))
        embedding = embedding[:n]
        metadata = metadata.iloc[:n]

    if HAS_UMAP and embedding.shape[0] > 15:
        try:
            reducer = umap.UMAP(n_neighbors=15, min_dist=0.1, random_state=42)
            coords_2d = reducer.fit_transform(embedding)
            method_name = "UMAP"
        except Exception as e:
            print(f"UMAP failed: {e}, falling back to PCA")
            if HAS_SKLEARN:
                pca = PCA(n_components=2, random_state=42)
                coords_2d = pca.fit_transform(embedding)
                method_name = "PCA"
            else:
                coords_2d = embedding[:, :2]
                method_name = "First 2 components"
    elif HAS_SKLEARN:
        pca = PCA(n_components=2, random_state=42)
        coords_2d = pca.fit_transform(embedding)
        method_name = "PCA"
    else:
        coords_2d = embedding[:, :2]
        method_name = "First 2 components"

    # Color by cell type
    cell_types = metadata['cell_type'].values
    unique_cell_types = list(set(cell_types))
    colors = [CELL_TYPE_COLORS.get(ct, '#808080') for ct in cell_types]

    # Plot
    scatter = ax.scatter(coords_2d[:, 0], coords_2d[:, 1],
                        c=colors, alpha=0.6, s=15, edgecolors='white', linewidths=0.3)

    ax.set_xlabel(f"{method_name} 1", fontsize=8)
    ax.set_ylabel(f"{method_name} 2", fontsize=8)
    ax.set_title("B. Program Embedding (colored by cell type)", fontsize=10, fontweight='bold', loc='left')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Legend - show only cell types present in data
    present_cell_types = [ct for ct in unique_cell_types if ct in CELL_TYPE_COLORS]
    handles = [mpatches.Patch(color=CELL_TYPE_COLORS.get(ct, '#808080'), label=ct.replace('_', ' '))
               for ct in sorted(present_cell_types)]
    ax.legend(handles=handles, loc='upper right', fontsize=5, framealpha=0.9,
              ncol=2, columnspacing=0.5, handletextpad=0.3)

    # Add annotation
    ax.text(0.02, 0.02, f'n={len(metadata)} programs\n14 samples',
            transform=ax.transAxes, fontsize=7, va='bottom')

    return coords_2d


def panel_c_response_enrichment(ax, response_data):
    """Panel C: Responder vs Progressor enriched programs."""
    responder_progs = response_data.get('responder_enriched', [])
    progressor_progs = response_data.get('progressor_enriched', [])

    # Create bar data
    categories = ['Responder-enriched', 'Progressor-enriched']
    counts = [len(responder_progs), len(progressor_progs)]
    colors = ['#2ecc71', '#e74c3c']  # Green for responder, red for progressor

    bars = ax.bar(categories, counts, color=colors, alpha=0.8, width=0.6)

    # Add count labels
    for bar, count in zip(bars, counts):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.1,
                str(count), ha='center', va='bottom', fontsize=10, fontweight='bold')

    ax.set_ylabel("Number of Programs", fontsize=8)
    ax.set_title("C. Response-Associated Programs", fontsize=10, fontweight='bold', loc='left')
    ax.set_ylim(0, max(counts) * 1.3)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Add details below bars
    y_offset = -0.15
    for i, (label, progs) in enumerate([('Responder', responder_progs), ('Progressor', progressor_progs)]):
        if progs:
            detail_text = '\n'.join([f"{p['cell_type']}: {', '.join(p['top_genes'][:3])}"
                                     for p in progs[:3]])
            ax.text(i, -0.5, detail_text,
                    ha='center', va='top', fontsize=6, transform=ax.get_xaxis_transform())

    # Adjust x-axis label position
    ax.tick_params(axis='x', pad=60)


def load_pydeseq2_results():
    """Load PyDESeq2 pseudo-bulk results."""
    deseq_file = PROJECT_ROOT / "examples/output_module5_pydeseq/pseudobulk_de_results.csv"
    if not deseq_file.exists():
        return None
    return pd.read_csv(deseq_file, index_col=0)


def panel_d_volcano_plot(ax):
    """Panel D: Volcano plot of PyDESeq2 results."""
    de_results = load_pydeseq2_results()

    if de_results is None:
        # Fallback to placeholder if results not available
        ax.text(0.5, 0.5, "Panel D: Differential Expression Analysis\n"
                "(PyDESeq2 results not found)",
                ha='center', va='center', fontsize=9, style='italic')
        ax.axis('off')
        ax.set_title("D. Differential Expression", fontsize=10, fontweight='bold', loc='left')
        return

    # Filter to valid results (non-NaN padj)
    de_valid = de_results.dropna(subset=['padj', 'log2FoldChange'])

    # Compute -log10(padj)
    de_valid = de_valid.copy()
    de_valid['neg_log10_padj'] = -np.log10(de_valid['padj'].clip(lower=1e-50))

    # Define significance
    sig_thresh = 0.05
    lfc_thresh = 1.0  # log2FC threshold for highlighting

    # Categorize genes
    de_valid['category'] = 'not_sig'
    de_valid.loc[(de_valid['padj'] < sig_thresh) & (de_valid['log2FoldChange'] > 0), 'category'] = 'responder_up'
    de_valid.loc[(de_valid['padj'] < sig_thresh) & (de_valid['log2FoldChange'] < 0), 'category'] = 'progressor_up'

    # Colors
    colors = {
        'not_sig': '#cccccc',
        'responder_up': '#2ecc71',  # Green
        'progressor_up': '#e74c3c'  # Red
    }

    # Plot all points
    for cat in ['not_sig', 'progressor_up', 'responder_up']:
        subset = de_valid[de_valid['category'] == cat]
        ax.scatter(subset['log2FoldChange'], subset['neg_log10_padj'],
                   c=colors[cat], alpha=0.5 if cat == 'not_sig' else 0.7,
                   s=8 if cat == 'not_sig' else 15, label=cat.replace('_', ' ').title(),
                   edgecolors='none')

    # Add significance threshold line
    ax.axhline(-np.log10(sig_thresh), color='gray', linestyle='--', linewidth=0.5, alpha=0.7)

    # Label top genes
    # Top 5 progressor-up genes (most significant negative log2FC)
    top_prog = de_valid[de_valid['category'] == 'progressor_up'].nsmallest(5, 'padj')
    for _, row in top_prog.iterrows():
        ax.annotate(row.name, (row['log2FoldChange'], row['neg_log10_padj']),
                    fontsize=6, ha='right', va='bottom',
                    xytext=(-3, 3), textcoords='offset points')

    # Top responder-up genes (most significant positive log2FC)
    top_resp = de_valid[de_valid['category'] == 'responder_up'].nsmallest(3, 'padj')
    for _, row in top_resp.iterrows():
        ax.annotate(row.name, (row['log2FoldChange'], row['neg_log10_padj']),
                    fontsize=6, ha='left', va='bottom',
                    xytext=(3, 3), textcoords='offset points')

    # Labels
    ax.set_xlabel("log2 Fold Change (Responder vs Progressor)", fontsize=8)
    ax.set_ylabel("-log10(adjusted p-value)", fontsize=8)
    ax.set_title("D. Differential Expression (Pseudo-bulk PyDESeq2)",
                 fontsize=10, fontweight='bold', loc='left')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Summary stats
    n_resp_up = (de_valid['category'] == 'responder_up').sum()
    n_prog_up = (de_valid['category'] == 'progressor_up').sum()
    n_total = len(de_valid)

    summary_text = f"n={n_total} genes\n{n_resp_up} responder↑\n{n_prog_up} progressor↑"
    ax.text(0.98, 0.98, summary_text, transform=ax.transAxes, fontsize=7,
            va='top', ha='right', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    # Legend
    ax.legend(loc='upper left', fontsize=6, framealpha=0.9, markerscale=1.5)


def panel_e_relationship_network(ax, relationships_df, aligned_programs_df):
    """Panel E: Conserved relationship network visualization."""
    if not HAS_NETWORKX:
        # Fallback: simple statistics display
        rel_counts = relationships_df['relationship_type'].value_counts()
        ax.text(0.5, 0.5, f"Conserved Relationships Summary\n\n"
                f"Co-localized: {rel_counts.get('co-localized', 0)}\n"
                f"Exclusive: {rel_counts.get('exclusive', 0)}\n"
                f"Independent: {rel_counts.get('independent', 0)}\n\n"
                f"Total: {len(relationships_df)} relationships",
                ha='center', va='center', fontsize=10)
        ax.axis('off')
        ax.set_title("E. Conserved Spatial Relationships", fontsize=10, fontweight='bold', loc='left')
        return

    # Build network from relationships
    G = nx.Graph()

    # Get program cell type mapping
    prog_to_celltype = dict(zip(aligned_programs_df['program_id'],
                                aligned_programs_df['cell_type']))

    # Filter to higher-conservation relationships for cleaner viz
    high_cons = relationships_df[relationships_df['n_samples'] >= 3].copy()

    if len(high_cons) < 10:
        high_cons = relationships_df.nlargest(30, 'n_samples')

    # Add nodes for programs that appear in relationships
    all_programs = set(high_cons['program1_id'].unique()) | set(high_cons['program2_id'].unique())

    for prog in all_programs:
        cell_type = prog_to_celltype.get(prog, 'Unknown')
        G.add_node(prog, cell_type=cell_type)

    # Add edges
    edge_colors = []
    edge_widths = []
    for _, row in high_cons.iterrows():
        p1, p2 = row['program1_id'], row['program2_id']
        rel_type = row['relationship_type']
        n_samples = row['n_samples']

        if rel_type == 'co-localized':
            color = '#2ecc71'  # Green
        elif rel_type == 'exclusive':
            color = '#e74c3c'  # Red
        else:
            color = '#95a5a6'  # Gray for independent

        G.add_edge(p1, p2, relationship=rel_type, weight=n_samples)
        edge_colors.append(color)
        edge_widths.append(0.5 + n_samples * 0.3)

    # Layout
    try:
        pos = nx.spring_layout(G, k=2, iterations=50, seed=42)
    except Exception:
        pos = nx.circular_layout(G)

    # Node colors by cell type
    node_colors = [CELL_TYPE_COLORS.get(G.nodes[n].get('cell_type', 'Unknown'), '#808080')
                   for n in G.nodes()]

    # Draw network
    nx.draw_networkx_edges(G, pos, ax=ax, edge_color=edge_colors,
                           width=edge_widths, alpha=0.6)
    nx.draw_networkx_nodes(G, pos, ax=ax, node_color=node_colors,
                           node_size=80, alpha=0.8, edgecolors='white', linewidths=0.5)

    # Add labels for highly connected nodes only
    degree = dict(G.degree())
    high_degree_nodes = {n: n[-3:] for n in G.nodes() if degree[n] >= 3}
    nx.draw_networkx_labels(G, pos, high_degree_nodes, ax=ax, font_size=5)

    ax.set_title("E. Conserved Spatial Relationships", fontsize=10, fontweight='bold', loc='left')
    ax.axis('off')

    # Legend for edge types
    legend_elements = [
        Line2D([0], [0], color='#2ecc71', linewidth=2, label='Co-localized'),
        Line2D([0], [0], color='#e74c3c', linewidth=2, label='Exclusive'),
        Line2D([0], [0], color='#95a5a6', linewidth=2, label='Independent'),
    ]
    ax.legend(handles=legend_elements, loc='lower right', fontsize=6, framealpha=0.9)

    # Add summary text
    rel_counts = relationships_df['relationship_type'].value_counts()
    summary = f"n={len(G.nodes())} programs, {len(G.edges())} relationships"
    ax.text(0.02, 0.02, summary, transform=ax.transAxes, fontsize=6, va='bottom')


def panel_e_relationship_summary(ax, relationships_df):
    """Panel E alternative: Relationship type distribution with conservation scores."""
    # Group by relationship type
    rel_summary = relationships_df.groupby('relationship_type').agg({
        'n_samples': ['count', 'mean', 'max'],
        'mean_bivariate_i': 'mean'
    }).round(3)
    rel_summary.columns = ['Count', 'Avg Samples', 'Max Samples', 'Mean Bivariate I']
    rel_summary = rel_summary.reset_index()

    # Create dual visualization
    # Left: pie chart of counts
    ax_pie = ax.inset_axes([0.0, 0.2, 0.4, 0.7])
    colors = {'co-localized': '#2ecc71', 'independent': '#95a5a6', 'exclusive': '#e74c3c'}
    pie_colors = [colors.get(r, '#808080') for r in rel_summary['relationship_type']]

    wedges, texts, autotexts = ax_pie.pie(
        rel_summary['Count'],
        labels=rel_summary['relationship_type'],
        autopct='%1.0f%%',
        colors=pie_colors,
        startangle=90,
        explode=[0.05 if r != 'independent' else 0 for r in rel_summary['relationship_type']]
    )
    for autotext in autotexts:
        autotext.set_fontsize(7)
    for text in texts:
        text.set_fontsize(7)
    ax_pie.set_title("Distribution", fontsize=8, y=1.05)

    # Right: bar chart of conservation
    ax_bar = ax.inset_axes([0.55, 0.2, 0.4, 0.6])
    x = range(len(rel_summary))
    bars = ax_bar.barh(x, rel_summary['Avg Samples'],
                       color=[colors.get(r, '#808080') for r in rel_summary['relationship_type']],
                       alpha=0.8)
    ax_bar.set_yticks(x)
    ax_bar.set_yticklabels([r.replace('-', '\n') for r in rel_summary['relationship_type']], fontsize=6)
    ax_bar.set_xlabel("Avg. Samples", fontsize=7)
    ax_bar.set_title("Conservation", fontsize=8)
    ax_bar.spines['top'].set_visible(False)
    ax_bar.spines['right'].set_visible(False)

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis('off')
    ax.set_title("E. Conserved Spatial Relationships (191 total, 14 samples)",
                 fontsize=10, fontweight='bold', loc='left')


def generate_figure5():
    """Generate complete Figure 5."""
    print("Loading Module 5 data...")

    # Load all data
    aligned_programs = load_aligned_programs()
    relationships = load_conserved_relationships()
    response_data = load_response_analysis()
    summary = load_summary()
    embedding = load_embedding()
    metadata = load_program_metadata()

    print(f"Loaded {len(aligned_programs)} aligned programs")
    print(f"Loaded {len(relationships)} conserved relationships")
    print(f"Embedding shape: {embedding.shape}")
    print(f"Metadata rows: {len(metadata)}")
    print(f"Response analysis: {summary.get('response_analysis', {})}")

    # Create figure
    fig = plt.figure(figsize=(12, 11))

    # Layout:
    # Row 1: Panel A (schematic) - spans left half, Panel B (embedding) - right half
    # Row 2: Panel C (response) - left, Panel D (DESeq placeholder) - right
    # Row 3: Panel E (network) - full width

    gs = GridSpec(3, 2, figure=fig, height_ratios=[1.2, 1.0, 1.2],
                  hspace=0.35, wspace=0.3)

    # Panel A: Schematic
    ax_a = fig.add_subplot(gs[0, 0])
    panel_a_schematic(ax_a)

    # Panel B: Embedding visualization
    ax_b = fig.add_subplot(gs[0, 1])
    panel_b_embedding_visualization(ax_b, embedding, metadata)

    # Panel C: Response enrichment
    ax_c = fig.add_subplot(gs[1, 0])
    panel_c_response_enrichment(ax_c, response_data)

    # Panel D: Volcano plot (PyDESeq2 results)
    ax_d = fig.add_subplot(gs[1, 1])
    panel_d_volcano_plot(ax_d)

    # Panel E: Relationship network (full width)
    ax_e = fig.add_subplot(gs[2, :])
    if HAS_NETWORKX:
        panel_e_relationship_network(ax_e, relationships, aligned_programs)
    else:
        panel_e_relationship_summary(ax_e, relationships)

    # Save
    output_path = OUTPUT_DIR / "figure5_integration.pdf"
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"Saved to {output_path}")

    # Also save PNG for quick preview
    png_path = OUTPUT_DIR / "figure5_integration.png"
    plt.savefig(png_path, dpi=150, bbox_inches='tight', facecolor='white')
    print(f"Preview saved to {png_path}")

    plt.close()

    # Print summary statistics
    print("\n=== Figure 5 Summary Statistics ===")
    print(f"\nModule 5 Cross-Sample Integration Summary:")
    print(f"  Samples: {summary.get('n_samples', 'N/A')}")
    print(f"  Aligned programs: {summary.get('n_aligned_programs', 'N/A')}")
    print(f"  Conserved relationships: {summary.get('n_conserved_relationships', 'N/A')}")
    print(f"  Highly conserved programs (>50% samples): {summary.get('n_highly_conserved_programs', 'N/A')}")

    print(f"\nRelationship Type Distribution:")
    rel_counts = relationships['relationship_type'].value_counts()
    for rel_type, count in rel_counts.items():
        pct = count / len(relationships) * 100
        print(f"  {rel_type}: {count} ({pct:.1f}%)")

    print(f"\nResponse-Associated Programs:")
    resp_analysis = summary.get('response_analysis', {})
    print(f"  Responder-enriched: {resp_analysis.get('n_responder_enriched', 'N/A')}")
    print(f"  Progressor-enriched: {resp_analysis.get('n_progressor_enriched', 'N/A')}")

    # Print details of response-enriched programs
    print("\n  Responder-enriched programs:")
    for prog in response_data.get('responder_enriched', []):
        genes = ', '.join(prog['top_genes'][:3])
        print(f"    - {prog['program_id']} ({prog['cell_type']}): {genes}")

    print("\n  Progressor-enriched programs:")
    for prog in response_data.get('progressor_enriched', []):
        genes = ', '.join(prog['top_genes'][:3])
        print(f"    - {prog['program_id']} ({prog['cell_type']}): {genes}")

    # Aligned program conservation
    print(f"\nProgram Conservation Across Samples:")
    conservation_bins = [
        (14, 14, '100% (all 14)'),
        (10, 13, '71-93%'),
        (5, 9, '36-64%'),
        (2, 4, '14-29%'),
        (1, 1, 'Single sample'),
    ]
    for low, high, label in conservation_bins:
        count = ((aligned_programs['n_samples'] >= low) &
                 (aligned_programs['n_samples'] <= high)).sum()
        print(f"  {label}: {count} programs")


if __name__ == "__main__":
    generate_figure5()
