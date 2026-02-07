#!/usr/bin/env python3
"""
Figure 3: Module 3 Deconvolution Benchmarking

Panel A: SCHEMATIC - use output/schematics/figure3_panel_a_deconvolution.svg
Panels B, C, D: DATA - generated with matplotlib below

Workflow:
1. Run this script to generate data panels (B, C, D)
2. Run svg_schematics.py to generate Panel A SVG
3. Combine in Illustrator for final figure
"""

import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from pathlib import Path
from adjustText import adjust_text

# Import shared style
from figure_style import apply_style, PALETTE, METHOD_COLORS, get_method_color, get_cell_type_color

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


def panel_b_proportion_benchmark(ax, method_data):
    """Panel B: Proportion benchmarking - grouped bar chart."""
    ax.text(-0.08, 1.08, "B", fontsize=14, fontweight='bold', va='top', transform=ax.transAxes)

    methods = method_data['method'].tolist()
    x = np.arange(len(methods))
    width = 0.25

    # Normalize metrics for comparison (higher = better for all)
    corr = method_data['pearson_r_mean'].values
    jsd_inv = 1 - method_data['jsd_mean'].values
    rmse_inv = 1 - method_data['rmse_mean'].values

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

    cell_types = [ct for ct in gex_data['CITEgeist'].keys() if not ct.startswith('_')]
    citegeist_corr = [gex_data['CITEgeist'][ct]['spot_pearson_r'] for ct in cell_types]
    scresolve_corr = [gex_data['scResolve'][ct]['spot_pearson_r'] for ct in cell_types]
    citegeist_cov = [gex_data['CITEgeist'][ct]['coverage'] for ct in cell_types]
    scresolve_cov = [gex_data['scResolve'][ct]['coverage'] for ct in cell_types]

    x = np.arange(len(cell_types))
    width = 0.35

    bars1 = ax.bar(x - width/2, citegeist_corr, width, label='CITEgeist', color=METHOD_COLORS['CITEgeist'], alpha=0.85)
    bars2 = ax.bar(x + width/2, scresolve_corr, width, label='scResolve', color=PALETTE['neutral'], alpha=0.85)

    # Coverage annotations with adjustText
    texts = []
    for i, (cg_cov, sr_cov) in enumerate(zip(citegeist_cov, scresolve_cov)):
        t1 = ax.text(i - width/2, citegeist_corr[i] + 0.02, f'{cg_cov*100:.0f}%',
                     ha='center', va='bottom', fontsize=7, color=METHOD_COLORS['CITEgeist'])
        t2 = ax.text(i + width/2, scresolve_corr[i] + 0.02, f'{sr_cov*100:.0f}%',
                     ha='center', va='bottom', fontsize=7, color=PALETTE['neutral'])
        texts.extend([t1, t2])

    adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle='-', color='gray', lw=0.3),
                expand_points=(1.5, 1.8), force_text=(0.8, 1.5), force_points=(0.3, 0.5),
                only_move={'points': 'y', 'text': 'xy'})

    short_labels = [ct.replace(' cells', '').replace('+ ', '+') for ct in cell_types]
    ax.set_xticks(x)
    ax.set_xticklabels(short_labels, rotation=45, ha='right', fontsize=8)
    ax.set_ylabel('GEX Correlation (r)', fontsize=10)
    ax.set_ylim(0, 0.7)
    ax.legend(loc='upper right', fontsize=8, framealpha=0.9)
    ax.set_title('GEX Benchmarking (% = coverage)', fontsize=11, fontweight='bold', loc='left')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    cg_overall = gex_data['CITEgeist']['_overall']['mean_spot_pearson_r']
    sr_overall = gex_data['scResolve']['_overall']['mean_spot_pearson_r']
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

        for ct in common_ct:
            color = get_cell_type_color(ct)
            ax.scatter(gt[ct].values, pred[ct].values, c=color, alpha=0.3, s=8,
                      label=ct.replace(' cells', '').replace('+ ', '+'), rasterized=True)

        ax.plot([0, 1], [0, 1], 'k--', lw=1.5, alpha=0.5)

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
    """Generate Figure 3 data panels."""
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

    # Create figure with data panels only (Panel A is schematic SVG)
    fig = plt.figure(figsize=(12, 10))
    gs = GridSpec(2, 2, figure=fig, hspace=0.35, wspace=0.30)

    # Panel A: Placeholder for schematic
    ax_a = fig.add_subplot(gs[0, 0])
    ax_a.text(0.5, 0.5, "Panel A: Schematic\n\nUse SVG file:\nfigure3_panel_a_deconvolution.svg",
              ha='center', va='center', fontsize=11, style='italic',
              bbox=dict(boxstyle='round', facecolor='#f0f0f0', edgecolor='gray'))
    ax_a.set_xlim(0, 1)
    ax_a.set_ylim(0, 1)
    ax_a.axis('off')
    ax_a.set_title("Two-Pass Deconvolution", fontsize=12, fontweight='bold', loc='left', pad=10)

    # Panel B: Proportion benchmark (DATA)
    ax_b = fig.add_subplot(gs[0, 1])
    panel_b_proportion_benchmark(ax_b, method_data)

    # Panel C: GEX benchmark (DATA)
    ax_c = fig.add_subplot(gs[1, 0])
    panel_c_gex_benchmark(ax_c, gex_data)

    # Panel D: Scatter (DATA)
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
    print("\n" + "=" * 60)
    print("Figure 3: Deconvolution Benchmarking")
    print("=" * 60)
    print("\nPanel A: SCHEMATIC - use output/schematics/figure3_panel_a_deconvolution.svg")
    print("Panels B, C, D: DATA - generated above")
    print("\nTo create final figure:")
    print("  1. Open figure3_benchmarking.svg in Illustrator")
    print("  2. Replace Panel A placeholder with figure3_panel_a_deconvolution.svg")
    print("  3. Export as figure3_benchmarking_final.pdf")


if __name__ == "__main__":
    generate_figure3()
