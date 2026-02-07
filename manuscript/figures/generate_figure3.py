#!/usr/bin/env python3
"""
Figure 3: Module 3 Deconvolution Benchmarking

Comprehensive benchmarking figure with BOTH proportion AND GEX benchmarking.

Panels:
  A: Two-pass deconvolution schematic
  B: Proportion benchmarking (3 metrics: Correlation, JSD, RMSE)
  C: GEX benchmarking (correlation + coverage by cell type)
  D: Predicted vs Ground Truth scatter
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
from figure_style import apply_style, PALETTE, CELL_TYPE_COLORS, METHOD_COLORS, get_method_color, get_cell_type_color

# Apply publication style
apply_style()

# Paths
PROJECT_ROOT = Path(__file__).parent.parent.parent
RESULTS_DIR = PROJECT_ROOT / "Benchmarking/xenium_benchmarking/evaluation/results/method_comparison"
GROUND_TRUTH_DIR = PROJECT_ROOT / "Benchmarking/xenium_pseudovisium/data/ground_truth"
CITEGEIST_OUTPUT_DIR = PROJECT_ROOT / "Benchmarking/xenium_benchmarking/CITEgeist/output_achievable_7"
GEX_RESULTS_FILE = PROJECT_ROOT / "Benchmarking/xenium_benchmarking/evaluation/results_gex_comparison_fair.json"
OUTPUT_DIR = Path(__file__).parent / "output"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def load_benchmark_summary():
    """Load proportion benchmark summary data."""
    try:
        df = pd.read_csv(RESULTS_DIR / "method_summary.csv", header=[0, 1], index_col=0)
        df.columns = ['_'.join(col).strip() for col in df.columns.values]
        return df
    except Exception as e:
        print(f"Error loading benchmark summary: {e}")
        return None


def load_gex_benchmark():
    """Load GEX benchmark data."""
    try:
        with open(GEX_RESULTS_FILE) as f:
            return json.load(f)
    except Exception as e:
        print(f"Error loading GEX benchmark: {e}")
        return None


def prepare_method_data(summary_df):
    """Prepare proportion method data for plotting."""
    method_mapping = {
        'CITEgeist_achievable_7': 'CITEgeist',
        'Cell2Location_achievable_7': 'Cell2Location',
        'RCTD_achievable_7': 'RCTD',
        'Tangram_achievable_7': 'Tangram',
        'Seurat_achievable_7': 'Seurat',
    }
    data = {'method': [], 'pearson_r_mean': [], 'pearson_r_std': [],
            'jsd_mean': [], 'jsd_std': [], 'rmse_mean': [], 'rmse_std': []}
    for idx_name, display_name in method_mapping.items():
        if idx_name in summary_df.index:
            data['method'].append(display_name)
            data['pearson_r_mean'].append(summary_df.loc[idx_name, 'pearson_r_mean'])
            data['pearson_r_std'].append(summary_df.loc[idx_name, 'pearson_r_std'])
            data['jsd_mean'].append(summary_df.loc[idx_name, 'jsd_mean'])
            data['jsd_std'].append(summary_df.loc[idx_name, 'jsd_std'])
            data['rmse_mean'].append(summary_df.loc[idx_name, 'rmse_mean'])
            data['rmse_std'].append(summary_df.loc[idx_name, 'rmse_std'])
    return pd.DataFrame(data)


def panel_a_schematic(ax):
    """Panel A: Clean two-pass deconvolution schematic."""
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis('off')

    # Panel label
    ax.text(-0.02, 1.02, "A", fontsize=14, fontweight='bold', va='top', transform=ax.transAxes)

    arrow_color = PALETTE['neutral']

    # Pass 1: Proportion estimation
    ax.text(0.5, 0.92, "Two-Pass Deconvolution", ha='center', fontsize=11, fontweight='bold')

    # Input: Mixed spot
    wedge_colors = [get_cell_type_color('T cells'), get_cell_type_color('Macrophages'),
                    get_cell_type_color('Fibroblasts'), get_cell_type_color('Epithelial')]
    for i, c in enumerate(wedge_colors):
        wedge = mpatches.Wedge((0.08, 0.55), 0.05, i*90, (i+1)*90,
                               facecolor=c, edgecolor='white', linewidth=1)
        ax.add_patch(wedge)
    ax.text(0.08, 0.68, "Mixed\nSpot", ha='center', va='bottom', fontsize=9)

    # Arrow to Pass 1
    ax.annotate('', xy=(0.20, 0.55), xytext=(0.14, 0.55),
                arrowprops=dict(arrowstyle='->', color=arrow_color, lw=1.5))

    # Pass 1 box
    box1 = FancyBboxPatch((0.20, 0.42), 0.22, 0.26, boxstyle="round,pad=0.02",
                          facecolor='#d5e8d4', edgecolor='#82b366', linewidth=2)
    ax.add_patch(box1)
    ax.text(0.31, 0.62, "Pass 1", ha='center', fontsize=10, fontweight='bold')
    ax.text(0.31, 0.52, "Protein → Proportions", ha='center', fontsize=9)
    ax.text(0.31, 0.45, "QP + Spatial smoothing", ha='center', fontsize=8, style='italic', color=PALETTE['neutral'])

    # Arrow to Pass 2
    ax.annotate('', xy=(0.50, 0.55), xytext=(0.43, 0.55),
                arrowprops=dict(arrowstyle='->', color=arrow_color, lw=1.5))

    # Pass 2 box
    box2 = FancyBboxPatch((0.50, 0.42), 0.22, 0.26, boxstyle="round,pad=0.02",
                          facecolor='#dae8fc', edgecolor='#6c8ebf', linewidth=2)
    ax.add_patch(box2)
    ax.text(0.61, 0.62, "Pass 2", ha='center', fontsize=10, fontweight='bold')
    ax.text(0.61, 0.52, "RNA → GEX Layers", ha='center', fontsize=9)
    ax.text(0.61, 0.45, "Per-cell-type expression", ha='center', fontsize=8, style='italic', color=PALETTE['neutral'])

    # Arrow to outputs
    ax.annotate('', xy=(0.80, 0.55), xytext=(0.73, 0.55),
                arrowprops=dict(arrowstyle='->', color=arrow_color, lw=1.5))

    # Outputs
    ax.text(0.90, 0.68, "Outputs", ha='center', fontsize=10, fontweight='bold')

    # Proportion bar
    bar_x, bar_y = 0.80, 0.56
    bar_w, bar_h = 0.18, 0.06
    props = [0.25, 0.30, 0.25, 0.20]
    cum = bar_x
    for p, c in zip(props, wedge_colors):
        w = p * bar_w
        rect = plt.Rectangle((cum, bar_y), w, bar_h, facecolor=c, edgecolor='white', linewidth=0.5)
        ax.add_patch(rect)
        cum += w
    ax.text(bar_x + bar_w/2, bar_y + bar_h + 0.02, "Proportions", ha='center', fontsize=8)

    # GEX layers representation
    for i, (c, label) in enumerate(zip(wedge_colors[:3], ['T', 'M', 'F'])):
        rect = plt.Rectangle((0.82 + i*0.05, 0.42), 0.04, 0.10, facecolor=c, edgecolor='white', linewidth=0.5, alpha=0.8)
        ax.add_patch(rect)
    ax.text(0.89, 0.40, "GEX Layers", ha='center', fontsize=8)

    # Key insight box
    ax.text(0.5, 0.18, "Key: Same-slide protein anchors deconvolution of both proportions and gene expression",
            ha='center', fontsize=9, style='italic', color=PALETTE['neutral'],
            bbox=dict(boxstyle='round,pad=0.3', facecolor=PALETTE['background'], edgecolor=PALETTE['border']))


def panel_b_proportion_benchmark(ax, method_data):
    """Panel B: Proportion benchmarking - grouped bar chart."""
    ax.text(-0.08, 1.08, "B", fontsize=14, fontweight='bold', va='top', transform=ax.transAxes)

    methods = method_data['method'].tolist()
    x = np.arange(len(methods))
    width = 0.25

    colors = [get_method_color(m) for m in methods]

    # Normalize metrics for comparison (higher = better for all)
    corr = method_data['pearson_r_mean'].values
    jsd_inv = 1 - method_data['jsd_mean'].values  # Invert so higher = better
    rmse_inv = 1 - method_data['rmse_mean'].values  # Invert so higher = better

    # Plot grouped bars
    bars1 = ax.bar(x - width, corr, width, label='Correlation', color=PALETTE['primary'], alpha=0.85)
    bars2 = ax.bar(x, jsd_inv, width, label='1 - JSD', color=PALETTE['accent1'], alpha=0.85)
    bars3 = ax.bar(x + width, rmse_inv, width, label='1 - RMSE', color=PALETTE['accent2'], alpha=0.85)

    # Highlight CITEgeist
    for bars in [bars1, bars2, bars3]:
        bars[0].set_edgecolor('black')
        bars[0].set_linewidth(1.5)

    ax.set_ylabel('Score (higher = better)', fontsize=10)
    ax.set_xticks(x)
    ax.set_xticklabels(methods, rotation=30, ha='right', fontsize=9)
    ax.set_ylim(0, 1.0)
    ax.legend(loc='upper right', fontsize=8, framealpha=0.9)
    ax.set_title('Proportion Benchmarking', fontsize=11, fontweight='bold', loc='left')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)


def panel_c_gex_benchmark(ax, gex_data):
    """Panel C: GEX benchmarking - correlation + coverage by cell type."""
    ax.text(-0.08, 1.08, "C", fontsize=14, fontweight='bold', va='top', transform=ax.transAxes)

    if gex_data is None:
        ax.text(0.5, 0.5, "GEX benchmark data not available", ha='center', va='center', fontsize=10)
        ax.axis('off')
        return

    # Extract cell types and metrics
    cell_types = [ct for ct in gex_data['CITEgeist'].keys() if not ct.startswith('_')]

    # Prepare data
    citegeist_corr = [gex_data['CITEgeist'][ct]['spot_pearson_r'] for ct in cell_types]
    scresolve_corr = [gex_data['scResolve'][ct]['spot_pearson_r'] for ct in cell_types]
    citegeist_cov = [gex_data['CITEgeist'][ct]['coverage'] for ct in cell_types]
    scresolve_cov = [gex_data['scResolve'][ct]['coverage'] for ct in cell_types]

    x = np.arange(len(cell_types))
    width = 0.35

    # Plot correlation bars
    bars1 = ax.bar(x - width/2, citegeist_corr, width, label='CITEgeist', color=METHOD_COLORS['CITEgeist'], alpha=0.85)
    bars2 = ax.bar(x + width/2, scresolve_corr, width, label='scResolve', color=PALETTE['neutral'], alpha=0.85)

    # Add coverage as text annotations
    for i, (cg_cov, sr_cov) in enumerate(zip(citegeist_cov, scresolve_cov)):
        ax.text(i - width/2, citegeist_corr[i] + 0.02, f'{cg_cov*100:.0f}%', ha='center', va='bottom', fontsize=7, color=METHOD_COLORS['CITEgeist'])
        ax.text(i + width/2, scresolve_corr[i] + 0.02, f'{sr_cov*100:.0f}%', ha='center', va='bottom', fontsize=7, color=PALETTE['neutral'])

    # Short cell type labels
    short_labels = [ct.replace(' cells', '').replace('+ ', '+') for ct in cell_types]
    ax.set_xticks(x)
    ax.set_xticklabels(short_labels, rotation=45, ha='right', fontsize=8)
    ax.set_ylabel('GEX Correlation (r)', fontsize=10)
    ax.set_ylim(0, 0.7)
    ax.legend(loc='upper right', fontsize=8, framealpha=0.9)
    ax.set_title('GEX Benchmarking (% = coverage)', fontsize=11, fontweight='bold', loc='left')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Overall summary
    cg_overall = gex_data['CITEgeist']['_overall']['mean_spot_pearson_r']
    sr_overall = gex_data['scResolve']['_overall']['mean_spot_pearson_r']
    cg_cov_overall = gex_data['CITEgeist']['_overall']['mean_coverage']
    sr_cov_overall = gex_data['scResolve']['_overall']['mean_coverage']

    summary_text = f"Overall: CITEgeist r={cg_overall:.2f} (100% cov) vs scResolve r={sr_overall:.2f} ({sr_cov_overall*100:.0f}% cov)"
    ax.text(0.5, -0.25, summary_text, ha='center', va='top', fontsize=9, style='italic',
            transform=ax.transAxes, color=PALETTE['neutral'])


def panel_d_scatter(ax, region_id=0):
    """Panel D: Predicted vs Ground Truth scatter plot."""
    ax.text(-0.08, 1.08, "D", fontsize=14, fontweight='bold', va='top', transform=ax.transAxes)

    try:
        gt = pd.read_csv(GROUND_TRUTH_DIR / f"Xenium_region_{region_id}_prop.csv", index_col=0)
        pred = pd.read_csv(CITEGEIST_OUTPUT_DIR / f"Xenium_region_{region_id}_cell_prop_finetuned_results.csv", index_col=0)

        common = list(set(gt.index) & set(pred.index))
        gt, pred = gt.loc[common], pred.loc[common]

        gt_cols = [c for c in gt.columns if c not in ['Unassigned', 'Unknown', 'n_cells', 'spot_x', 'spot_y']]
        pred_cols = [c for c in pred.columns if c not in ['Unassigned', 'Unknown']]
        common_ct = [c for c in gt_cols if c in pred_cols]

        # Plot each cell type with its color
        for ct in common_ct:
            color = get_cell_type_color(ct)
            ax.scatter(gt[ct].values, pred[ct].values, c=color, alpha=0.3, s=8,
                      label=ct.replace(' cells', '').replace('+ ', '+'), rasterized=True)

        # Diagonal line
        ax.plot([0, 1], [0, 1], 'k--', lw=1.5, alpha=0.5)

        # Calculate overall correlation
        all_gt = np.concatenate([gt[ct].values for ct in common_ct])
        all_pred = np.concatenate([pred[ct].values for ct in common_ct])
        corr = np.corrcoef(all_gt, all_pred)[0, 1]

        ax.text(0.05, 0.92, f'r = {corr:.3f}', transform=ax.transAxes, fontsize=11, fontweight='bold')

        ax.set_xlabel("Ground Truth", fontsize=10)
        ax.set_ylabel("CITEgeist", fontsize=10)
        ax.set_xlim(-0.02, 1.02)
        ax.set_ylim(-0.02, 1.02)
        ax.set_aspect('equal')
        ax.legend(loc='lower right', fontsize=7, framealpha=0.9, markerscale=1.5)
        ax.set_title('Predicted vs Ground Truth', fontsize=11, fontweight='bold', loc='left')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    except Exception as e:
        ax.text(0.5, 0.5, f"Data not available\n({str(e)[:40]}...)",
                ha='center', va='center', fontsize=10, style='italic', color=PALETTE['neutral'])
        ax.set_title('Predicted vs Ground Truth', fontsize=11, fontweight='bold', loc='left')


def generate_figure3():
    """Generate complete Figure 3."""
    print("Loading data...")
    summary_df = load_benchmark_summary()
    gex_data = load_gex_benchmark()

    if summary_df is None:
        print("ERROR: Could not load proportion benchmark summary")
        return

    method_data = prepare_method_data(summary_df)
    print(f"Proportion methods: {method_data['method'].tolist()}")

    if gex_data:
        print(f"GEX methods: {list(gex_data.keys())}")
    else:
        print("WARNING: GEX benchmark data not available")

    # Create figure with 2x2 layout
    fig = plt.figure(figsize=(12, 10))
    gs = GridSpec(2, 2, figure=fig, hspace=0.35, wspace=0.30)

    # Panel A: Schematic (top left)
    ax_a = fig.add_subplot(gs[0, 0])
    panel_a_schematic(ax_a)

    # Panel B: Proportion benchmark (top right)
    ax_b = fig.add_subplot(gs[0, 1])
    panel_b_proportion_benchmark(ax_b, method_data)

    # Panel C: GEX benchmark (bottom left)
    ax_c = fig.add_subplot(gs[1, 0])
    panel_c_gex_benchmark(ax_c, gex_data)

    # Panel D: Scatter (bottom right)
    ax_d = fig.add_subplot(gs[1, 1])
    panel_d_scatter(ax_d, region_id=0)

    plt.tight_layout()

    # Save outputs
    for fmt, dpi in [('pdf', 300), ('png', 150), ('svg', None)]:
        output_path = OUTPUT_DIR / f"figure3_benchmarking.{fmt}"
        if fmt == 'svg':
            plt.savefig(output_path, format='svg', bbox_inches='tight', facecolor='white')
        else:
            plt.savefig(output_path, dpi=dpi, bbox_inches='tight', facecolor='white')
        print(f"Saved: {output_path}")

    plt.close()

    # Print summary
    print("\n=== Figure 3 Summary ===")
    print("Panels:")
    print("  A: Two-pass deconvolution schematic")
    print("  B: Proportion benchmarking (5 methods, 3 metrics)")
    print("  C: GEX benchmarking (CITEgeist vs scResolve, 7 cell types)")
    print("  D: Predicted vs Ground Truth scatter")

    if gex_data:
        cg = gex_data['CITEgeist']['_overall']
        sr = gex_data['scResolve']['_overall']
        print(f"\nGEX Benchmark Summary:")
        print(f"  CITEgeist: r={cg['mean_spot_pearson_r']:.3f}, coverage={cg['mean_coverage']*100:.0f}%")
        print(f"  scResolve: r={sr['mean_spot_pearson_r']:.3f}, coverage={sr['mean_coverage']*100:.0f}%")


if __name__ == "__main__":
    generate_figure3()
