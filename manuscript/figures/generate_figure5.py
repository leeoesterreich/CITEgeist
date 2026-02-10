#!/usr/bin/env python3
"""
Figure 5: Cross-Sample Integration and Response Biology (Section 2.5)

Dense 8-panel figure (4 rows x 2 columns, 12x16 inches) showing Module 3-5
cross-sample integration results and spatial response biology.

Panel A: Cell type proportion shift (biopsy -> surgical) grouped by response
Panel B: Program conservation grouped bar chart (programs per cell type, by response)
Panel C: Response enrichment dot plot (responder fraction per program)
Panel D: Conserved relationship network (65 nodes, 163 edges) with expanded layout
Panel E: Bivariate spatial co-localization (Responder P3-S2) with H&E backdrop
Panel F: Bivariate spatial exclusion (Progressor P1-S2) with H&E backdrop
Panel G: Moran's I spatial coherence validation (Module 4)
Panel H: Spatial proportion comparison (2x2 grid) with H&E backdrop
"""

import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
from matplotlib.colors import ListedColormap, Normalize, LinearSegmentedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable
from pathlib import Path
import warnings

import networkx as nx

# For H&E backdrop spatial plotting
import squidpy as sq
import scanpy as sc

from figure_style import apply_style, PALETTE, CELL_TYPE_COLORS, get_cell_type_color

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
PROJECT_ROOT = Path(__file__).parent.parent.parent
DATA_ROOT = Path(
    "/ix1/alee/LO_LAB/General/Lab_Data/"
    "20250210_CITEGeistPublicData_GEO_Alex/processed_files"
)
MODULE3_DIR = PROJECT_ROOT / "output/module3_unified"
MODULE4_DIR = PROJECT_ROOT / "output/module4_unified"
MODULE5_DIR = PROJECT_ROOT / "output/module5_integration"
OUTPUT_DIR = Path(__file__).parent / "output"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

SAMPLE_MAP = {
    'P1': {'S1': 'HCC22-088-P1-S1', 'S2': 'HCC22-088-P1-S2'},
    'P2': {'S1': 'HCC22-088-P2-S1', 'S2': 'HCC22-088-P2-S2'},
    'P3': {'S1': 'HCC22-088-P3-S1_A', 'S2': 'HCC22-088-P3-S2'},
    'P4': {'S1': 'HCC22-088-P4-S1', 'S2': 'HCC22-088-P4-S2_1i_rep'},
    'P5': {'S1': 'HCC22-088-P5-S1', 'S2': 'HCC22-088-P5-S2_F_rep'},
    'P6': {'S1': 'HCC22-088-P6-S1', 'S2': 'HCC22-088-P6-S2_D'},
}

RESPONSE_MAP = {
    'P1': 'progressor', 'P2': 'responder', 'P3': 'responder',
    'P4': 'progressor', 'P5': 'responder', 'P6': 'responder',
}

RESPONSE_COLORS = {'responder': '#2171b5', 'progressor': '#e74c3c'}
SPATIAL_POINT_SIZE = 1.6
SPATIAL_ALPHA_IMG = 0.55
SPATIAL_ALPHA_POINTS = 0.95

CELL_TYPES = [
    'Endothelial', 'Fibroblasts', 'B_Cells', 'Macrophages', 'Monocytes',
    'CD8_T_Cells', 'CD4_T_Cells', 'Cancer_Luminal', 'Cancer_Basal',
    'Dendritic_Cells',
]

# Ordered sample list (patient then timepoint)
SAMPLE_ORDER = []
for patient in ['P1', 'P2', 'P3', 'P4', 'P5', 'P6']:
    for tp in ['S1', 'S2']:
        SAMPLE_ORDER.append(SAMPLE_MAP[patient][tp])

# Patient from sample name helper
RESPONDER_SAMPLES = set()
PROGRESSOR_SAMPLES = set()
for patient, response in RESPONSE_MAP.items():
    if response == 'responder':
        for tp_sample in SAMPLE_MAP[patient].values():
            RESPONDER_SAMPLES.add(tp_sample)
    else:
        for tp_sample in SAMPLE_MAP[patient].values():
            PROGRESSOR_SAMPLES.add(tp_sample)

# Cell type abbreviation map used across panels
CELL_TYPE_ABBREV = {
    'Endothelial': 'Endo', 'Fibroblasts': 'Fibro', 'B_Cells': 'B',
    'Macrophages': 'Mac', 'Monocytes': 'Mono', 'CD8_T_Cells': 'CD8+T',
    'CD4_T_Cells': 'CD4+T', 'Cancer_Luminal': 'Ca Lum',
    'Cancer_Basal': 'Ca Bas', 'Dendritic_Cells': 'DC',
}


# ---------------------------------------------------------------------------
# Data loading helpers
# ---------------------------------------------------------------------------

def load_cell_proportions(sample_name):
    """Load finetuned cell proportions CSV for a given sample."""
    path = MODULE3_DIR / f"{sample_name}_cell_prop_finetuned_results.csv"
    if path.exists():
        return pd.read_csv(path, index_col=0)
    return None


def load_spatial_positions(sample_name):
    """Load tissue positions for a given sample, filtering in_tissue=1."""
    path = DATA_ROOT / sample_name / "outs/spatial/tissue_positions.csv"
    if path.exists():
        pos = pd.read_csv(path)
        return pos[pos['in_tissue'] == 1].set_index('barcode')
    return None


def load_aligned_programs():
    """Load Module 5 aligned programs CSV."""
    return pd.read_csv(MODULE5_DIR / "module5_unified_aligned_programs.csv")


def load_conserved_relationships():
    """Load Module 5 conserved relationships CSV."""
    return pd.read_csv(
        MODULE5_DIR / "module5_unified_conserved_relationships.csv"
    )


def load_response_analysis():
    """Load Module 5 response analysis JSON."""
    with open(MODULE5_DIR / "module5_response_analysis.json", 'r') as f:
        return json.load(f)


def load_module4_morans():
    """Load Moran's I values from all 12 Module 4 JSONs."""
    records = []
    for patient, timepoints in SAMPLE_MAP.items():
        for tp, sample in timepoints.items():
            path = MODULE4_DIR / f"{sample}_module4_discovery.json"
            if not path.exists():
                continue
            with open(path, 'r') as f:
                d = json.load(f)
            anchors = d.get('anchored', {}).get('anchors', {})
            for cell_type, info in anchors.items():
                for prog in info.get('programs', []):
                    records.append({
                        'sample': sample,
                        'patient': patient,
                        'cell_type': cell_type,
                        'spatial_moran_i': prog['spatial_moran_i'],
                    })
    return pd.DataFrame(records)


def compute_responder_frac(samples_str):
    """Given comma-separated sample names, compute fraction that are responders."""
    samples = [s.strip() for s in samples_str.split(',')]
    n_resp = sum(1 for s in samples if s in RESPONDER_SAMPLES)
    return n_resp / len(samples) if samples else 0.0


def get_sample_patient(sample_name):
    """Extract patient ID from sample name, e.g., HCC22-088-P3-S2 -> P3."""
    for patient, timepoints in SAMPLE_MAP.items():
        if sample_name in timepoints.values():
            return patient
    return None


def load_visium_with_image(sample_name):
    """Load Visium data with H&E image for sc.pl.spatial() plotting.

    Returns AnnData object with spatial coordinates and image loaded.
    """
    path = DATA_ROOT / sample_name / "outs"
    if not path.exists():
        return None

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        adata = sq.read.visium(
            path,
            counts_file='filtered_feature_bc_matrix.h5',
            load_images=True,
            gex_only=False
        )

    return adata


def plot_spatial_with_he(ax, adata, values, title="", cmap='viridis',
                          vmin=None, vmax=None, spot_size=1.0,
                          show_colorbar=True, colorbar_label=None):
    """Plot values on spatial coordinates with H&E backdrop using sc.pl.spatial.

    Args:
        ax: matplotlib axes
        adata: AnnData with spatial coordinates and image
        values: array of values to plot (same length as adata.n_obs)
        title: plot title
        cmap: colormap name
        vmin/vmax: color scale limits
        spot_size: multiplier for spot size
        show_colorbar: whether to show colorbar
        colorbar_label: label for colorbar

    Returns:
        The PathCollection for colorbar reference.
    """
    # Add values to adata.obs
    adata.obs['_plot_values'] = values

    if vmin is None:
        vmin = np.percentile(values[~np.isnan(values)], 1)
    if vmax is None:
        vmax = np.percentile(values[~np.isnan(values)], 99)

    # Use sc.pl.spatial for H&E backdrop
    sc.pl.spatial(
        adata,
        color='_plot_values',
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        spot_size=spot_size,
        alpha_img=0.8,
        ax=ax,
        show=False,
        colorbar_loc=None if not show_colorbar else 'right',
        title=title,
    )

    # Clean up temporary column
    del adata.obs['_plot_values']

    return ax


# ---------------------------------------------------------------------------
# Panel A: Cell Type Proportion Shift (Biopsy -> Surgical)
# ---------------------------------------------------------------------------

def panel_a_proportion_shift(ax):
    """Grouped bar chart of mean proportion change by response group."""
    deltas = []  # list of dicts: {patient, response, cell_type, delta}
    for patient, timepoints in SAMPLE_MAP.items():
        s1_name = timepoints['S1']
        s2_name = timepoints['S2']
        props_s1 = load_cell_proportions(s1_name)
        props_s2 = load_cell_proportions(s2_name)
        if props_s1 is None or props_s2 is None:
            continue
        mean_s1 = props_s1[CELL_TYPES].mean()
        mean_s2 = props_s2[CELL_TYPES].mean()
        delta = mean_s2 - mean_s1
        for ct in CELL_TYPES:
            deltas.append({
                'patient': patient,
                'response': RESPONSE_MAP[patient],
                'cell_type': ct,
                'delta': delta[ct],
            })

    df = pd.DataFrame(deltas)
    if df.empty:
        ax.text(0.5, 0.5, "Data not available", ha='center', va='center',
                fontsize=10, style='italic', color=PALETTE['neutral'])
        ax.axis('off')
        return

    x = np.arange(len(CELL_TYPES))
    bar_width = 0.35

    for i, (response, color) in enumerate([
        ('responder', RESPONSE_COLORS['responder']),
        ('progressor', RESPONSE_COLORS['progressor']),
    ]):
        sub = df[df['response'] == response]
        means = sub.groupby('cell_type')['delta'].mean().reindex(CELL_TYPES, fill_value=0)
        offset = (i - 0.5) * bar_width
        bars = ax.bar(x + offset, means, bar_width, color=color, alpha=0.7,
                      label=response.capitalize(), edgecolor='white', linewidth=0.5)

        # Overlay individual patient dots
        for ct_idx, ct in enumerate(CELL_TYPES):
            ct_sub = sub[sub['cell_type'] == ct]
            jitter = np.random.default_rng(42).uniform(-0.06, 0.06, len(ct_sub))
            ax.scatter(
                [ct_idx + offset] * len(ct_sub) + jitter,
                ct_sub['delta'].values,
                color='black', s=12, alpha=0.7, zorder=5, edgecolors='white',
                linewidths=0.3,
            )

    ax.axhline(0, color='gray', linewidth=0.5, linestyle='-')
    ax.set_xticks(x)
    ax.set_xticklabels([CELL_TYPE_ABBREV.get(ct, ct) for ct in CELL_TYPES],
                       rotation=45, ha='right', fontsize=9)
    ax.set_ylabel("Mean proportion change\n(Surgical - Biopsy)", fontsize=10)
    ax.legend(fontsize=9, loc='upper right', framealpha=0.9)

    n_resp = len([p for p in RESPONSE_MAP if RESPONSE_MAP[p] == 'responder'])
    n_prog = len([p for p in RESPONSE_MAP if RESPONSE_MAP[p] == 'progressor'])
    ax.text(0.02, 0.98, f'Resp n={n_resp}, Prog n={n_prog}',
            transform=ax.transAxes, ha='left', va='top', fontsize=9,
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))


# ---------------------------------------------------------------------------
# Panel B: Program Conservation Summary (Redesigned)
# ---------------------------------------------------------------------------

def panel_b_program_summary(ax):
    """Summary showing NMF gene programs by cell type with biological context.

    Programs = NMF-discovered gene expression patterns within each cell type layer.
    Each bar shows # programs found, with top genes annotated to explain biology.
    """
    aligned = load_aligned_programs()

    # Count programs per cell type
    ct_counts = aligned.groupby('cell_type').agg({
        'program_id': 'count',
        'conservation_score': 'mean',
        'n_samples': 'mean',
    }).rename(columns={'program_id': 'n_programs'})

    # Sort by number of programs
    ct_counts = ct_counts.sort_values('n_programs', ascending=True)

    # Create horizontal bar chart
    y_pos = np.arange(len(ct_counts))
    bars = ax.barh(y_pos, ct_counts['n_programs'],
                   color=[get_cell_type_color(ct) for ct in ct_counts.index],
                   edgecolor='white', linewidth=0.5, alpha=0.85)

    # Add clear count labels only (remove long per-row gene annotations).
    for i, (ct, row) in enumerate(ct_counts.iterrows()):
        n_progs = int(row['n_programs'])
        ax.text(row['n_programs'] + 0.25, i,
                f"{n_progs}",
                va='center', ha='left', fontsize=9, fontweight='bold')

        # Conservation indicator (star if highly conserved)
        if row['conservation_score'] > 0.5:
            ax.text(row['n_programs'] - 0.5, i, '*',
                    va='center', ha='right', fontsize=12,
                    fontweight='bold', color='white')

    # Y-axis labels with cell type names
    ax.set_yticks(y_pos)
    ax.set_yticklabels([CELL_TYPE_ABBREV.get(ct, ct) for ct in ct_counts.index],
                       fontsize=9)

    ax.set_xlabel("Aligned programs (count)", fontsize=10)

    # Add explanatory summary
    total_progs = len(aligned)
    highly_conserved = (aligned['conservation_score'] > 0.5).sum()
    ax.text(0.98, 0.98,
            f"Total programs: {total_progs}\n"
            f"Highly conserved (>6 samples): {highly_conserved}\n"
            f"Star (*) marks mean conservation > 0.5",
            transform=ax.transAxes, ha='right', va='top', fontsize=8,
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.9))
    ax.set_xlim(0, ct_counts['n_programs'].max() * 1.25)


# ---------------------------------------------------------------------------
# Panel C: Response Enrichment Dot Plot
# ---------------------------------------------------------------------------

def panel_c_response_dotplot(ax):
    """Dot plot: responder fraction vs program, sized by n_samples, colored by cell type."""
    aligned = load_aligned_programs()
    response_data = load_response_analysis()

    # Compute responder fraction for each program
    aligned['responder_frac'] = aligned['samples'].apply(compute_responder_frac)

    # Sort by cell type
    aligned = aligned.sort_values(
        ['cell_type', 'conservation_score'], ascending=[True, False]
    ).reset_index(drop=True)

    y_positions = np.arange(len(aligned))
    colors = [get_cell_type_color(ct) for ct in aligned['cell_type']]
    sizes = aligned['n_samples'].values * 8  # scale for visibility

    ax.scatter(aligned['responder_frac'], y_positions, c=colors, s=sizes,
               alpha=0.7, edgecolors='white', linewidths=0.3, rasterized=True)

    # Baseline: 8/12 = 0.667
    ax.axvline(8 / 12, color='gray', linestyle='--', linewidth=1.0, alpha=0.7,
               label='Expected (8/12)')

    # Highlight and label significant programs
    highlight_programs = {
        'aligned_004': ('Macrophages', 'resp'),
        'aligned_008': ('CD4_T_Cells', 'prog'),
        'aligned_002': ('Cancer_Luminal', 'prog'),
        'aligned_016': ('Monocytes', 'prog'),
    }

    texts = []
    for prog_id, (ct, direction) in highlight_programs.items():
        mask = aligned['program_id'] == prog_id
        if not mask.any():
            continue
        idx = mask.idxmax()
        row = aligned.loc[idx]
        y = y_positions[idx]
        x = row['responder_frac']

        # Emphasize with ring
        ax.scatter([x], [y], s=row['n_samples'] * 8 + 60,
                   edgecolors='black', facecolors='none', linewidths=1.5, zorder=6)

        # Get top genes from response_data
        gene_str = ""
        for entry_list_key in ['responder_enriched', 'progressor_enriched']:
            for entry in response_data.get(entry_list_key, []):
                if entry['program_id'] == prog_id:
                    gene_str = ", ".join(entry['top_genes'][:3])
                    break

        label = f"{prog_id}\n({ct.replace('_', ' ')})\n{gene_str}"
        t = ax.annotate(
            label, (x, y), fontsize=9, ha='left', va='center',
            xytext=(8, 0), textcoords='offset points',
            bbox=dict(boxstyle='round,pad=0.2', facecolor='white',
                      alpha=0.8, edgecolor='gray', linewidth=0.5),
        )
        texts.append(t)

    ax.set_xlabel("Responder fraction of samples", fontsize=10)
    ax.set_ylabel("Aligned programs (by cell type)", fontsize=10)
    ax.set_yticks([])
    ax.set_xlim(-0.05, 1.05)
    ax.legend(fontsize=9, loc='upper left')

    # Cell type sidebar
    ct_groups = aligned['cell_type'].values
    for i in range(len(ct_groups)):
        color = get_cell_type_color(ct_groups[i])
        ax.add_patch(plt.Rectangle((-0.08, i - 0.5), 0.03, 1,
                                   color=color, clip_on=False,
                                   transform=ax.get_yaxis_transform()))


# ---------------------------------------------------------------------------
# Panel D: Conserved Relationship Network (Expanded with Cell Type Labels)
# ---------------------------------------------------------------------------

def panel_d_network(ax):
    """Readable summary network of key conserved spatial relationships."""
    aligned = load_aligned_programs()
    relationships = load_conserved_relationships()

    G = nx.Graph()

    # Add nodes with attributes
    for _, row in aligned.iterrows():
        G.add_node(row['program_id'],
                   cell_type=row['cell_type'],
                   conservation_score=row['conservation_score'],
                   n_samples=row['n_samples'],
                   top_genes=row.get('top_genes', ''))

    # Add edges with attributes
    for _, row in relationships.iterrows():
        p1, p2 = row['program1_id'], row['program2_id']
        if p1 in G.nodes and p2 in G.nodes:
            G.add_edge(p1, p2,
                       mean_bivariate_i=row['mean_bivariate_i'],
                       relationship_type=row['relationship_type'])

    # Four key relationships highlighted in the manuscript narrative.
    biological_modules = {
        'stromal_immune': {
            'nodes': ['aligned_012', 'aligned_050'],
            'I': 0.358,
            'label': 'Stromal-Immune\nAxis',
            'cell_types': 'Fibroblast + CD4 T',
            'color': '#2ca02c',
        },
        'antigen_presentation': {
            'nodes': ['aligned_000', 'aligned_050'],
            'I': 0.258,
            'label': 'Antigen\nPresentation',
            'cell_types': 'DC + CD4 T',
            'color': '#1f77b4',
        },
        't_helper_diversity': {
            'nodes': ['aligned_008', 'aligned_050'],
            'I': 0.178,
            'label': 'T Helper\nDiversity',
            'cell_types': 'CD4 T subtypes',
            'color': '#9467bd',
        },
        'immune_exclusion': {
            'nodes': ['aligned_003', 'aligned_021'],
            'I': -0.194,
            'label': 'Immune\nExclusion',
            'cell_types': 'Cancer vs CD4 T',
            'color': '#d62728',
        },
    }

    # Build deterministic positions for only key nodes to avoid label pileups.
    module_centers = {
        'stromal_immune': (0.70, 0.70),
        'antigen_presentation': (0.32, 0.74),
        't_helper_diversity': (0.56, 0.42),
        'immune_exclusion': (0.30, 0.28),
    }
    node_positions = {}
    node_hits = {}
    offsets = [(-0.05, 0.00), (0.05, 0.00)]
    for name, module in biological_modules.items():
        cx, cy = module_centers[name]
        for idx, node in enumerate(module['nodes']):
            if node not in G.nodes:
                continue
            dx, dy = offsets[idx]
            node_positions.setdefault(node, np.array([0.0, 0.0]))
            node_hits[node] = node_hits.get(node, 0) + 1
            node_positions[node] += np.array([cx + dx, cy + dy])
    for node in list(node_positions.keys()):
        node_positions[node] = node_positions[node] / node_hits[node]

    # Draw module edges first.
    for module in biological_modules.values():
        n1, n2 = module['nodes']
        if n1 not in node_positions or n2 not in node_positions:
            continue
        x1, y1 = node_positions[n1]
        x2, y2 = node_positions[n2]
        linestyle = '-' if module['I'] > 0 else '--'
        ax.plot([x1, x2], [y1, y2], color=module['color'], linewidth=2.5,
                linestyle=linestyle, alpha=0.9, zorder=2)

    # Draw nodes with cell-type color and size by conservation.
    for node, (x, y) in node_positions.items():
        ct = G.nodes[node].get('cell_type', 'Unknown')
        node_color = get_cell_type_color(ct)
        size = 0.03 + 0.06 * float(G.nodes[node].get('conservation_score', 0.2))
        circle = plt.Circle((x, y), size, facecolor=node_color, edgecolor='black',
                            linewidth=1.0, alpha=0.95, zorder=3)
        ax.add_patch(circle)

    # Module callouts: spread initial positions, then repel with anchors.
    callout_pos = {
        'stromal_immune': (0.82, 0.88),
        'antigen_presentation': (0.16, 0.87),
        't_helper_diversity': (0.70, 0.11),
        'immune_exclusion': (0.15, 0.08),
    }
    callout_texts = []
    anchor_x = []
    anchor_y = []
    for name, module in biological_modules.items():
        txt_x, txt_y = callout_pos[name]
        sign = '+' if module['I'] > 0 else ''
        bg_color = '#e8f5e9' if module['I'] > 0 else '#ffebee'
        n1, n2 = module['nodes']
        if n1 in node_positions and n2 in node_positions:
            mid = (node_positions[n1] + node_positions[n2]) / 2.0
            anchor_x.append(float(mid[0]))
            anchor_y.append(float(mid[1]))
        else:
            anchor_x.append(txt_x)
            anchor_y.append(txt_y)
        callout_texts.append(
            ax.text(
                txt_x, txt_y,
                f"{module['label']}\n{module['cell_types']}\nI={sign}{module['I']:.2f}",
                fontsize=7, fontweight='bold', ha='center', va='center',
                bbox=dict(boxstyle='round,pad=0.30', facecolor=bg_color,
                          edgecolor=module['color'], linewidth=1.5, alpha=0.95),
                zorder=4,
            )
        )

    # Repel labels to reduce residual spacing collisions in dense layouts.
    try:
        from adjustText import adjust_text

        adjust_text(
            callout_texts,
            target_x=anchor_x,
            target_y=anchor_y,
            ax=ax,
            only_move={"text": "xy", "static": "xy", "explode": "xy", "pull": "xy"},
            force_text=0.5,
            force_static=0.2,
            force_pull=0.06,
            expand_text=(1.2, 1.35),
            arrowprops=dict(arrowstyle='-', color='#666666', lw=0.8, alpha=0.7),
            ensure_inside_axes=True,
        )
    except Exception:
        pass

    # Legend for edge types
    legend_handles = [
        plt.Line2D([0], [0], color='#2ca02c', linewidth=3,
                   label=f'Co-localized (n={(relationships["relationship_type"] == "co-localized").sum()})'),
        plt.Line2D([0], [0], color='#d62728', linewidth=3, linestyle='dashed',
                   label=f'Exclusive (n={(relationships["relationship_type"] == "exclusive").sum()})'),
    ]
    ax.legend(handles=legend_handles, fontsize=8, loc='lower right',
              framealpha=0.95, title='Spatial Relationships', title_fontsize=9)

    # Summary box with biological interpretation
    summary_text = (
        f'{G.number_of_nodes()} conserved programs\n'
        f'{(relationships["relationship_type"] == "co-localized").sum()} co-localized pairs\n'
        f'{(relationships["relationship_type"] == "exclusive").sum()} exclusive pairs\n'
        f'Node size = conservation score'
    )
    ax.text(0.02, 0.98, summary_text,
            transform=ax.transAxes, ha='left', va='top', fontsize=7,
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.9))

    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.0)
    ax.axis('off')


# ---------------------------------------------------------------------------
# Spatial plotting helper
# ---------------------------------------------------------------------------

def plot_spatial_proportion(ax, sample_name, cell_type, cmap='viridis',
                            vmin=None, vmax=None, show_colorbar=True,
                            colorbar_label=None):
    """Plot cell type proportion on tissue spatial coordinates.

    Returns the PathCollection (scatter) for optional colorbar sharing.
    """
    props = load_cell_proportions(sample_name)
    positions = load_spatial_positions(sample_name)

    if props is None or positions is None:
        ax.text(0.5, 0.5, "Data not available", ha='center', va='center',
                fontsize=9, style='italic', color=PALETTE['neutral'])
        ax.axis('off')
        return None

    # Merge on barcode index
    common = props.index.intersection(positions.index)
    props = props.loc[common]
    positions = positions.loc[common]

    x = positions['pxl_col_in_fullres'].values
    y = positions['pxl_row_in_fullres'].values
    values = props[cell_type].values

    if vmin is None:
        vmin = np.percentile(values, 1)
    if vmax is None:
        vmax = np.percentile(values, 99)

    sc = ax.scatter(x, y, c=values, cmap=cmap, s=4, alpha=0.8,
                    vmin=vmin, vmax=vmax, edgecolors='none', rasterized=True)

    ax.invert_yaxis()
    ax.set_aspect('equal')
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)

    if show_colorbar:
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="4%", pad=0.05)
        cbar = plt.colorbar(sc, cax=cax)
        if colorbar_label:
            cbar.set_label(colorbar_label, fontsize=8)
        cbar.ax.tick_params(labelsize=8)

    return sc


# ---------------------------------------------------------------------------
# Panel E: Co-localization (Responder P3-S2) with H&E Backdrop
# ---------------------------------------------------------------------------

def panel_e_colocalization(ax):
    """Bivariate spatial co-localization: Fibroblast + CD4 T cell on P3-S2 with H&E.

    Redesigned to:
    - Use sc.pl.spatial() with H&E tissue backdrop
    - Show clear cell type labels
    - Include bivariate scatter with trend line
    - Viridis colormap for clarity
    """
    sample = 'HCC22-088-P3-S2'
    ct1 = 'Fibroblasts'
    ct2 = 'CD4_T_Cells'

    # Load Visium data with H&E image
    adata = load_visium_with_image(sample)
    props = load_cell_proportions(sample)

    if adata is None or props is None:
        ax.text(0.5, 0.5, "Data not available", ha='center', va='center',
                fontsize=9, style='italic', color=PALETTE['neutral'])
        ax.axis('off')
        return

    # Align barcodes
    common = adata.obs_names.intersection(props.index)
    if len(common) == 0:
        ax.text(0.5, 0.5, "No matching barcodes", ha='center', va='center',
                fontsize=9, style='italic', color=PALETTE['neutral'])
        ax.axis('off')
        return

    adata = adata[common, :].copy()
    props = props.loc[common]

    v1 = props[ct1].values
    v2 = props[ct2].values

    ax.set_axis_off()

    # Keep maps on top, a dedicated legend band in the middle, and scatter below.
    left_ax = ax.inset_axes([0.0, 0.40, 0.48, 0.56])
    right_ax = ax.inset_axes([0.52, 0.40, 0.48, 0.56])
    scatter_ax = ax.inset_axes([0.08, 0.02, 0.84, 0.17])

    # Left: Fibroblasts with H&E - use Oranges colormap for fibroblasts
    fib_vmin = np.percentile(v1, 5)
    fib_vmax = np.percentile(v1, 95)
    adata.obs['_fibroblast_prop'] = v1
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sc.pl.spatial(
            adata, color='_fibroblast_prop', cmap='Oranges',
            vmin=fib_vmin, vmax=fib_vmax,
            size=SPATIAL_POINT_SIZE, alpha_img=SPATIAL_ALPHA_IMG,
            alpha=SPATIAL_ALPHA_POINTS, ax=left_ax, show=False,
            colorbar_loc=None, title='',
        )
    left_ax.set_title('Fibroblasts', fontsize=11, fontweight='bold', pad=4,
                      color='#d95f02')
    left_ax.set_xlabel("")
    left_ax.set_ylabel("")

    # Right: CD4 T Cells with H&E - use Purples colormap for T cells
    cd4_vmin = np.percentile(v2, 5)
    cd4_vmax = np.percentile(v2, 95)
    adata.obs['_cd4t_prop'] = v2
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sc.pl.spatial(
            adata, color='_cd4t_prop', cmap='Purples',
            vmin=cd4_vmin, vmax=cd4_vmax,
            size=SPATIAL_POINT_SIZE, alpha_img=SPATIAL_ALPHA_IMG,
            alpha=SPATIAL_ALPHA_POINTS, ax=right_ax, show=False,
            colorbar_loc=None, title='',
        )
    right_ax.set_title('CD4+ T Cells', fontsize=11, fontweight='bold', pad=4,
                       color='#7570b3')
    right_ax.set_xlabel("")
    right_ax.set_ylabel("")

    # Explicit color legends for both spatial maps.
    fib_cax = ax.inset_axes([0.06, 0.255, 0.34, 0.020])
    fib_sm = plt.cm.ScalarMappable(cmap='Oranges', norm=Normalize(vmin=fib_vmin, vmax=fib_vmax))
    fib_sm.set_array([])
    fib_cbar = plt.colorbar(fib_sm, cax=fib_cax, orientation='horizontal')
    fib_cbar.ax.tick_params(labelsize=6, pad=1)
    ax.text(0.23, 0.279, "Fibro proportion", transform=ax.transAxes,
            ha='center', va='bottom', fontsize=7)

    cd4_cax = ax.inset_axes([0.58, 0.255, 0.34, 0.020])
    cd4_sm = plt.cm.ScalarMappable(cmap='Purples', norm=Normalize(vmin=cd4_vmin, vmax=cd4_vmax))
    cd4_sm.set_array([])
    cd4_cbar = plt.colorbar(cd4_sm, cax=cd4_cax, orientation='horizontal')
    cd4_cbar.ax.tick_params(labelsize=6, pad=1)
    ax.text(0.75, 0.279, "CD4+ T proportion", transform=ax.transAxes,
            ha='center', va='bottom', fontsize=7)

    # Scatter plot: v1 vs v2 with regression line
    scatter_ax.scatter(v1, v2, s=3, alpha=0.4, c='#636363', edgecolors='none',
                       rasterized=True)

    # Add regression line
    from scipy import stats
    slope, intercept, r_value, p_value, std_err = stats.linregress(v1, v2)
    x_line = np.linspace(v1.min(), v1.max(), 100)
    y_line = slope * x_line + intercept
    scatter_ax.plot(x_line, y_line, 'g-', linewidth=2, alpha=0.8,
                    label=f'r={r_value:.2f}')

    scatter_ax.set_xlabel('Fibro prop', fontsize=8)
    scatter_ax.set_ylabel('CD4+T prop', fontsize=8)
    scatter_ax.tick_params(labelsize=7)
    scatter_ax.spines['top'].set_visible(False)
    scatter_ax.spines['right'].set_visible(False)
    scatter_ax.legend(fontsize=7, loc='upper right', framealpha=0.8)

    # Compact annotation block above scatter.
    ax.text(
        0.5, 0.205,
        "Co-localization Moran's I = 0.358",
        transform=ax.transAxes, ha='center', va='center', fontsize=8,
        fontweight='bold',
        bbox=dict(boxstyle='round,pad=0.25', facecolor='#e8f5e9', alpha=0.95,
                  edgecolor='#2ca02c', linewidth=1.3),
    )
    ax.text(
        0.5, 0.314,
        "Color encodes cell-type proportion per spot (5th-95th percentile scale)",
        transform=ax.transAxes, ha='center', va='center', fontsize=7,
        color=PALETTE['neutral'],
    )

    # Figure title already carries sample context; skip extra sample badge to avoid crowding.


# ---------------------------------------------------------------------------
# Panel F: Exclusion (Progressor P1-S2) with H&E Backdrop
# ---------------------------------------------------------------------------

def panel_f_exclusion(ax):
    """Bivariate spatial exclusion: Cancer Luminal vs CD4 T cell on P1-S2 with H&E.

    Redesigned to:
    - Use sc.pl.spatial() with H&E tissue backdrop
    - Show clear cell type labels
    - Include bivariate scatter with negative trend
    - Emphasize exclusion pattern
    """
    sample = 'HCC22-088-P1-S2'
    ct1 = 'Cancer_Luminal'
    ct2 = 'CD4_T_Cells'

    # Load Visium data with H&E image
    adata = load_visium_with_image(sample)
    props = load_cell_proportions(sample)

    if adata is None or props is None:
        ax.text(0.5, 0.5, "Data not available", ha='center', va='center',
                fontsize=9, style='italic', color=PALETTE['neutral'])
        ax.axis('off')
        return

    # Align barcodes
    common = adata.obs_names.intersection(props.index)
    if len(common) == 0:
        ax.text(0.5, 0.5, "No matching barcodes", ha='center', va='center',
                fontsize=9, style='italic', color=PALETTE['neutral'])
        ax.axis('off')
        return

    adata = adata[common, :].copy()
    props = props.loc[common]

    v1 = props[ct1].values
    v2 = props[ct2].values

    ax.set_axis_off()

    # Keep maps on top, a dedicated legend band in the middle, and scatter below.
    left_ax = ax.inset_axes([0.0, 0.40, 0.48, 0.56])
    right_ax = ax.inset_axes([0.52, 0.40, 0.48, 0.56])
    scatter_ax = ax.inset_axes([0.08, 0.02, 0.84, 0.17])

    # Left: Cancer Luminal with H&E - use Blues colormap for cancer
    cancer_vmin = np.percentile(v1, 5)
    cancer_vmax = np.percentile(v1, 95)
    adata.obs['_cancer_prop'] = v1
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sc.pl.spatial(
            adata, color='_cancer_prop', cmap='Blues',
            vmin=cancer_vmin, vmax=cancer_vmax,
            size=SPATIAL_POINT_SIZE, alpha_img=SPATIAL_ALPHA_IMG,
            alpha=SPATIAL_ALPHA_POINTS, ax=left_ax, show=False,
            colorbar_loc=None, title='',
        )
    left_ax.set_title('Cancer (Luminal)', fontsize=11, fontweight='bold', pad=4,
                      color='#2171b5')
    left_ax.set_xlabel("")
    left_ax.set_ylabel("")

    # Right: CD4 T Cells with H&E - use Purples colormap for T cells
    cd4_vmin = np.percentile(v2, 5)
    cd4_vmax = np.percentile(v2, 95)
    adata.obs['_cd4t_prop'] = v2
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sc.pl.spatial(
            adata, color='_cd4t_prop', cmap='Purples',
            vmin=cd4_vmin, vmax=cd4_vmax,
            size=SPATIAL_POINT_SIZE, alpha_img=SPATIAL_ALPHA_IMG,
            alpha=SPATIAL_ALPHA_POINTS, ax=right_ax, show=False,
            colorbar_loc=None, title='',
        )
    right_ax.set_title('CD4+ T Cells', fontsize=11, fontweight='bold', pad=4,
                       color='#7570b3')
    right_ax.set_xlabel("")
    right_ax.set_ylabel("")

    # Explicit color legends for both spatial maps.
    cancer_cax = ax.inset_axes([0.06, 0.255, 0.34, 0.020])
    cancer_sm = plt.cm.ScalarMappable(cmap='Blues', norm=Normalize(vmin=cancer_vmin, vmax=cancer_vmax))
    cancer_sm.set_array([])
    cancer_cbar = plt.colorbar(cancer_sm, cax=cancer_cax, orientation='horizontal')
    cancer_cbar.ax.tick_params(labelsize=6, pad=1)
    ax.text(0.23, 0.279, "Cancer proportion", transform=ax.transAxes,
            ha='center', va='bottom', fontsize=7)

    cd4_cax = ax.inset_axes([0.58, 0.255, 0.34, 0.020])
    cd4_sm = plt.cm.ScalarMappable(cmap='Purples', norm=Normalize(vmin=cd4_vmin, vmax=cd4_vmax))
    cd4_sm.set_array([])
    cd4_cbar = plt.colorbar(cd4_sm, cax=cd4_cax, orientation='horizontal')
    cd4_cbar.ax.tick_params(labelsize=6, pad=1)
    ax.text(0.75, 0.279, "CD4+ T proportion", transform=ax.transAxes,
            ha='center', va='bottom', fontsize=7)

    # Scatter plot: v1 vs v2 with regression line (negative slope expected)
    scatter_ax.scatter(v1, v2, s=3, alpha=0.4, c='#636363', edgecolors='none',
                       rasterized=True)

    # Add regression line
    from scipy import stats
    slope, intercept, r_value, p_value, std_err = stats.linregress(v1, v2)
    x_line = np.linspace(v1.min(), v1.max(), 100)
    y_line = slope * x_line + intercept
    scatter_ax.plot(x_line, y_line, 'r-', linewidth=2, alpha=0.8,
                    label=f'r={r_value:.2f}')

    scatter_ax.set_xlabel('Cancer prop', fontsize=8)
    scatter_ax.set_ylabel('CD4+T prop', fontsize=8)
    scatter_ax.tick_params(labelsize=7)
    scatter_ax.spines['top'].set_visible(False)
    scatter_ax.spines['right'].set_visible(False)
    scatter_ax.legend(fontsize=7, loc='upper right', framealpha=0.8)

    ax.text(
        0.5, 0.205,
        "Exclusion Moran's I = -0.194",
        transform=ax.transAxes, ha='center', va='center', fontsize=8,
        fontweight='bold',
        bbox=dict(boxstyle='round,pad=0.25', facecolor='#ffebee', alpha=0.95,
                  edgecolor='#d62728', linewidth=1.3),
    )
    ax.text(
        0.5, 0.314,
        "Color encodes cell-type proportion per spot (5th-95th percentile scale)",
        transform=ax.transAxes, ha='center', va='center', fontsize=7,
        color=PALETTE['neutral'],
    )

    # Figure title already carries sample context; skip extra sample badge to avoid crowding.


# ---------------------------------------------------------------------------
# Panel G: Moran's I Spatial Coherence Validation
# ---------------------------------------------------------------------------

def panel_g_morans_violin(ax):
    """Violin + strip plot of Moran's I values by cell type across all 12 samples."""
    morans_df = load_module4_morans()
    if morans_df.empty:
        ax.text(0.5, 0.5, "Module 4 data not available", ha='center',
                va='center', fontsize=10, style='italic',
                color=PALETTE['neutral'])
        ax.axis('off')
        return

    # Determine cell types present, sorted
    present_cts = sorted(morans_df['cell_type'].unique())

    data_by_ct = [morans_df[morans_df['cell_type'] == ct]['spatial_moran_i'].values
                  for ct in present_cts]

    parts = ax.violinplot(data_by_ct, positions=range(len(present_cts)),
                          showmeans=False, showmedians=True, showextrema=False)

    for i, (pc, ct) in enumerate(zip(parts['bodies'], present_cts)):
        pc.set_facecolor(get_cell_type_color(ct))
        pc.set_alpha(0.5)

    # Median line color
    parts['cmedians'].set_color('black')
    parts['cmedians'].set_linewidth(1.0)

    # Strip overlay
    for i, (vals, ct) in enumerate(zip(data_by_ct, present_cts)):
        jitter = np.random.default_rng(42).uniform(-0.15, 0.15, len(vals))
        ax.scatter(np.full(len(vals), i) + jitter, vals, s=6, alpha=0.5,
                   color=get_cell_type_color(ct), edgecolors='white',
                   linewidths=0.2, zorder=5)

    # Threshold line
    ax.axhline(0.15, color=PALETTE['highlight'], linestyle='--', linewidth=1.2,
               alpha=0.7, label='Threshold (I=0.15)')

    ax.set_xticks(range(len(present_cts)))
    ax.set_xticklabels([CELL_TYPE_ABBREV.get(ct, ct[:6]) for ct in present_cts],
                       rotation=45, ha='right', fontsize=9)
    ax.set_ylabel("Moran's I", fontsize=10)
    ax.legend(fontsize=9, loc='upper right')

    # Annotate n above threshold
    all_vals = morans_df['spatial_moran_i']
    n_above = (all_vals > 0.15).sum()
    total = len(all_vals)
    pct = n_above / total * 100 if total > 0 else 0
    ax.text(0.02, 0.98, f'{n_above}/{total} ({pct:.0f}%) above threshold',
            transform=ax.transAxes, ha='left', va='top', fontsize=9,
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))


# ---------------------------------------------------------------------------
# Panel H: Spatial Proportion Comparison (2x2 grid) with H&E Backdrop
# ---------------------------------------------------------------------------

def panel_h_spatial_comparison(ax):
    """2x2 spatial proportion comparison with H&E backdrop.

    Shows Fibroblast + Cancer_Luminal proportions for P3-S2 (responder) vs P1-S2 (progressor)
    to visually demonstrate the key biological difference.
    """
    ax.set_axis_off()

    # Layout: 2 rows x 2 cols
    # Top row: Fibroblasts, Bottom row: Cancer_Luminal
    # Left col: P3-S2 (responder), Right col: P1-S2 (progressor)
    # Using viridis for all spatial proportion plots (perceptually uniform)
    panels = [
        (0.02, 0.53, 0.45, 0.43, 'HCC22-088-P3-S2', 'Fibroblasts',
         'viridis', 'Fibro (Resp)', '#2171b5'),
        (0.53, 0.53, 0.45, 0.43, 'HCC22-088-P1-S2', 'Fibroblasts',
         'viridis', 'Fibro (Prog)', '#e74c3c'),
        (0.02, 0.05, 0.45, 0.43, 'HCC22-088-P3-S2', 'Cancer_Luminal',
         'viridis', 'Cancer (Resp)', '#2171b5'),
        (0.53, 0.05, 0.45, 0.43, 'HCC22-088-P1-S2', 'Cancer_Luminal',
         'viridis', 'Cancer (Prog)', '#e74c3c'),
    ]

    # Precompute global vmin/vmax for each cell type across both samples
    vranges = {}
    for ct in ['Fibroblasts', 'Cancer_Luminal']:
        all_vals = []
        for sample in ['HCC22-088-P3-S2', 'HCC22-088-P1-S2']:
            props = load_cell_proportions(sample)
            if props is not None and ct in props.columns:
                all_vals.extend(props[ct].values)
        if all_vals:
            vranges[ct] = (np.percentile(all_vals, 2), np.percentile(all_vals, 98))
        else:
            vranges[ct] = (0, 1)

    for (x0, y0, w, h, sample, ct, cmap, title, resp_color) in panels:
        sub_ax = ax.inset_axes([x0, y0, w, h])

        # Load Visium with H&E
        adata = load_visium_with_image(sample)
        props = load_cell_proportions(sample)

        if adata is None or props is None:
            sub_ax.text(0.5, 0.5, "N/A", ha='center', va='center',
                        fontsize=9, color=PALETTE['neutral'])
            sub_ax.axis('off')
            continue

        # Align barcodes
        common = adata.obs_names.intersection(props.index)
        if len(common) == 0:
            sub_ax.text(0.5, 0.5, "No data", ha='center', va='center',
                        fontsize=9, color=PALETTE['neutral'])
            sub_ax.axis('off')
            continue

        adata = adata[common, :].copy()
        props_common = props.loc[common]
        vals = props_common[ct].values

        # Add to adata and plot with H&E
        adata.obs['_prop'] = vals
        vmin, vmax = vranges.get(ct, (np.percentile(vals, 2), np.percentile(vals, 98)))

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            sc.pl.spatial(
                adata, color='_prop', cmap=cmap,
                vmin=np.percentile(vals, 5), vmax=np.percentile(vals, 95),
                size=SPATIAL_POINT_SIZE, alpha_img=SPATIAL_ALPHA_IMG,
                alpha=SPATIAL_ALPHA_POINTS, ax=sub_ax, show=False,
                colorbar_loc=None, title='',
            )
        sub_ax.set_xlabel("")
        sub_ax.set_ylabel("")

        sub_ax.set_title(title, fontsize=9, fontweight='bold', pad=2, color=resp_color)

        # Add mean proportion annotation
        mean_prop = np.mean(vals)
        sub_ax.text(0.98, 0.02, f'mean={mean_prop:.2f}', transform=sub_ax.transAxes,
                    ha='right', va='bottom', fontsize=8,
                    bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    # Add row labels
    ax.text(-0.04, 0.74, 'Fibroblasts', transform=ax.transAxes,
            ha='right', va='center', fontsize=9, fontweight='bold',
            rotation=90)
    ax.text(-0.04, 0.25, 'Cancer', transform=ax.transAxes,
            ha='right', va='center', fontsize=9, fontweight='bold',
            rotation=90)


# ---------------------------------------------------------------------------
# Main figure assembly
# ---------------------------------------------------------------------------

def generate_figure5():
    """Generate Figure 5: Cross-Sample Integration and Response Biology.

    Assembles 8 panels in 4 rows x 2 columns showing Module 3-5 results:
    cell type proportion shifts, program conservation, response enrichment,
    spatial relationships, bivariate co-localization/exclusion, Moran's I
    validation, and spatial proportion comparisons.
    """
    # Disable constrained_layout since we use manual GridSpec
    apply_style()
    plt.rcParams['figure.constrained_layout.use'] = False

    fig = plt.figure(figsize=(12, 16))

    # GridSpec: 4 rows x 2 cols with custom height ratios
    gs = gridspec.GridSpec(
        4, 2, figure=fig,
        height_ratios=[1.0, 1.0, 1.0, 1.0],
        hspace=0.30, wspace=0.30,
        left=0.06, right=0.94, top=0.96, bottom=0.04,
    )

    # --- Panel A (row 0, col 0) ---
    print("Generating Panel A: Cell type proportion shift...")
    ax_a = fig.add_subplot(gs[0, 0])
    ax_a.text(-0.12, 1.08, "A", fontsize=14, fontweight='bold', va='top',
              transform=ax_a.transAxes)
    panel_a_proportion_shift(ax_a)
    ax_a.set_title("Cell Type Proportion Change\n(Biopsy \u2192 Surgical)",
                    fontsize=10, fontweight='bold', loc='left', pad=8)

    # --- Panel B (row 0, col 1) ---
    print("Generating Panel B: Program conservation summary...")
    ax_b = fig.add_subplot(gs[0, 1])
    ax_b.text(-0.12, 1.08, "B", fontsize=14, fontweight='bold', va='top',
              transform=ax_b.transAxes)
    panel_b_program_summary(ax_b)
    ax_b.set_title("Aligned Programs by Cell Type",
                    fontsize=10, fontweight='bold', loc='left', pad=8)

    # --- Panel C (row 1, col 0) ---
    print("Generating Panel C: Response enrichment dot plot...")
    ax_c = fig.add_subplot(gs[1, 0])
    ax_c.text(-0.12, 1.08, "C", fontsize=14, fontweight='bold', va='top',
              transform=ax_c.transAxes)
    panel_c_response_dotplot(ax_c)
    ax_c.set_title("Response Enrichment",
                    fontsize=10, fontweight='bold', loc='left', pad=8)

    # --- Panel D (row 1, col 1) ---
    print("Generating Panel D: Conserved relationship network...")
    ax_d = fig.add_subplot(gs[1, 1])
    ax_d.text(-0.12, 1.08, "D", fontsize=14, fontweight='bold', va='top',
              transform=ax_d.transAxes)
    panel_d_network(ax_d)
    ax_d.set_title("Conserved Spatial Relationships",
                    fontsize=10, fontweight='bold', loc='left', pad=8)

    # --- Panel E (row 2, col 0) ---
    print("Generating Panel E: Co-localization (P3-S2)...")
    ax_e = fig.add_subplot(gs[2, 0])
    ax_e.text(-0.12, 1.08, "E", fontsize=14, fontweight='bold', va='top',
              transform=ax_e.transAxes)
    panel_e_colocalization(ax_e)
    ax_e.set_title("Co-localization (Responder P3-S2)",
                    fontsize=10, fontweight='bold', loc='left', pad=8)

    # --- Panel F (row 2, col 1) ---
    print("Generating Panel F: Exclusion (P1-S2)...")
    ax_f = fig.add_subplot(gs[2, 1])
    ax_f.text(-0.12, 1.08, "F", fontsize=14, fontweight='bold', va='top',
              transform=ax_f.transAxes)
    panel_f_exclusion(ax_f)
    ax_f.set_title("Exclusion (Progressor P1-S2)",
                    fontsize=10, fontweight='bold', loc='left', pad=8)

    # --- Panel G (row 3, col 0) ---
    print("Generating Panel G: Moran's I validation...")
    ax_g = fig.add_subplot(gs[3, 0])
    ax_g.text(-0.12, 1.08, "G", fontsize=14, fontweight='bold', va='top',
              transform=ax_g.transAxes)
    panel_g_morans_violin(ax_g)
    ax_g.set_title("Spatial Coherence Validation",
                    fontsize=10, fontweight='bold', loc='left', pad=8)

    # --- Panel H (row 3, col 1) ---
    print("Generating Panel H: Spatial proportion comparison...")
    ax_h = fig.add_subplot(gs[3, 1])
    ax_h.text(-0.12, 1.08, "H", fontsize=14, fontweight='bold', va='top',
              transform=ax_h.transAxes)
    panel_h_spatial_comparison(ax_h)
    ax_h.set_title("Spatial Proportion Comparison",
                    fontsize=10, fontweight='bold', loc='left', pad=8)

    # Save in multiple formats
    for fmt, dpi in [('pdf', 300), ('png', 150), ('svg', None)]:
        output_path = OUTPUT_DIR / f"figure5_full_pipeline.{fmt}"
        try:
            if fmt == 'svg':
                plt.savefig(output_path, format='svg', bbox_inches='tight',
                            facecolor='white')
            else:
                plt.savefig(output_path, dpi=dpi, bbox_inches='tight',
                            facecolor='white')
            print(f"Saved: {output_path}")
        except Exception as e:
            print(f"Warning: Could not save {fmt}: {e}")

    plt.close()

    # Summary
    print("\n" + "=" * 70)
    print("Figure 5: Cross-Sample Integration and Response Biology")
    print("=" * 70)
    print("\nPanel A: Cell type proportion shift (biopsy -> surgical)")
    print("Panel B: Aligned programs by cell type (bar chart)")
    print("Panel C: Response enrichment dot plot")
    print("Panel D: Conserved relationship network (expanded with labels)")
    print("Panel E: Co-localization - Fibroblasts & CD4 T Cells (P3-S2, H&E)")
    print("Panel F: Exclusion - Cancer Luminal & CD4 T Cells (P1-S2, H&E)")
    print("Panel G: Moran's I spatial coherence validation")
    print("Panel H: Spatial proportion comparison (2x2, H&E backdrop)")

    # Load summary statistics
    aligned = load_aligned_programs()
    relationships = load_conserved_relationships()
    response = load_response_analysis()
    morans_df = load_module4_morans()

    print(f"\nModule 5: {len(aligned)} aligned programs, "
          f"{len(relationships)} relationships")
    print(f"  Responder-enriched: "
          f"{len(response.get('responder_enriched', []))}")
    print(f"  Progressor-enriched: "
          f"{len(response.get('progressor_enriched', []))}")

    n_coloc = (relationships['relationship_type'] == 'co-localized').sum()
    n_excl = (relationships['relationship_type'] == 'exclusive').sum()
    print(f"  Co-localized: {n_coloc}, Exclusive: {n_excl}")

    if not morans_df.empty:
        n_above = (morans_df['spatial_moran_i'] > 0.15).sum()
        total = len(morans_df)
        print(f"\nModule 4: {n_above}/{total} ({n_above/total*100:.0f}%) "
              f"programs above I=0.15 threshold")


if __name__ == "__main__":
    generate_figure5()
