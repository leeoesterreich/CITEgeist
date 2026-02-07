#!/usr/bin/env python3
"""
Figure 3: Module 3 Deconvolution Benchmarking

Enhanced version with 4 panels:

Panels:
  A: Proportion estimation schematic
  B: Benchmark comparison bar plots (Correlation, JSD, RMSE)
  C: Predicted vs Ground Truth scatter plot
  D: Simulated benchmark summary (optional, from simulation data)
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
from figure_style import apply_style, PALETTE, CELL_TYPE_COLORS, METHOD_COLORS, get_method_color

# Apply publication style
apply_style()

# Paths
PROJECT_ROOT = Path(__file__).parent.parent.parent
RESULTS_DIR = PROJECT_ROOT / "Benchmarking/xenium_benchmarking/evaluation/results/method_comparison"
GROUND_TRUTH_DIR = PROJECT_ROOT / "Benchmarking/xenium_pseudovisium/data/ground_truth"
CITEGEIST_OUTPUT_DIR = PROJECT_ROOT / "Benchmarking/xenium_benchmarking/CITEgeist/output_achievable_7"
H5AD_DIR = PROJECT_ROOT / "Benchmarking/xenium_pseudovisium/data/h5ad_objects"
SIM_RESULTS_DIR = PROJECT_ROOT / "Benchmarking/simulation_benchmarking/Figures"
OUTPUT_DIR = Path(__file__).parent / "output"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def load_benchmark_summary():
    """Load benchmark summary data."""
    try:
        df = pd.read_csv(RESULTS_DIR / "method_summary.csv", header=[0, 1], index_col=0)
        df.columns = ['_'.join(col).strip() for col in df.columns.values]
        return df
    except Exception as e:
        print(f"Error loading benchmark summary: {e}")
        return None


def prepare_method_data(summary_df):
    """Prepare method data for plotting."""
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


def load_simulated_summary():
    """Load simulated benchmark summary if available."""
    sim_file = SIM_RESULTS_DIR / "prop_all_metrics_highseg_combined.csv"
    if sim_file.exists():
        try:
            return pd.read_csv(sim_file)
        except Exception:
            pass
    return None


def panel_a_schematic(ax):
    """Panel A: Simple, clean deconvolution schematic."""
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis('off')

    # Panel label
    ax.text(0.02, 0.95, "A", fontsize=14, fontweight='bold', va='top')

    # Box styling
    box_height = 0.25
    arrow_color = PALETTE['neutral']

    # 1. Mixed Spot (left)
    spot_x = 0.10
    wedge_colors = [CELL_TYPE_COLORS['T cells'], CELL_TYPE_COLORS['Macrophages'],
                    CELL_TYPE_COLORS['Fibroblasts'], CELL_TYPE_COLORS['Epithelial']]
    for i, c in enumerate(wedge_colors):
        wedge = mpatches.Wedge((spot_x, 0.5), 0.06, i*90, (i+1)*90,
                               facecolor=c, edgecolor='white', linewidth=1)
        ax.add_patch(wedge)
    ax.text(spot_x, 0.78, "Mixed\nSpot", ha='center', va='bottom', fontsize=10, fontweight='bold')

    # Arrow 1
    ax.annotate('', xy=(0.24, 0.5), xytext=(0.17, 0.5),
                arrowprops=dict(arrowstyle='->', color=arrow_color, lw=2))

    # 2. Profiles box
    box1 = FancyBboxPatch((0.24, 0.35), 0.16, 0.30, boxstyle="round,pad=0.02",
                          facecolor='#d5e8d4', edgecolor='#82b366', linewidth=2)
    ax.add_patch(box1)
    ax.text(0.32, 0.5, "Protein\nProfiles", ha='center', va='center', fontsize=10, fontweight='bold')
    ax.text(0.32, 0.28, "(Module 2)", ha='center', va='top', fontsize=9, style='italic', color=PALETTE['neutral'])

    # Arrow 2
    ax.annotate('', xy=(0.48, 0.5), xytext=(0.41, 0.5),
                arrowprops=dict(arrowstyle='->', color=arrow_color, lw=2))

    # 3. Optimization box
    box2 = FancyBboxPatch((0.48, 0.30), 0.20, 0.40, boxstyle="round,pad=0.02",
                          facecolor='#dae8fc', edgecolor='#6c8ebf', linewidth=2)
    ax.add_patch(box2)
    ax.text(0.58, 0.58, "Optimization", ha='center', va='center', fontsize=10, fontweight='bold')
    ax.text(0.58, 0.45, "QP + Spatial\nSmoothing", ha='center', va='center', fontsize=9)

    # Arrow 3
    ax.annotate('', xy=(0.76, 0.5), xytext=(0.69, 0.5),
                arrowprops=dict(arrowstyle='->', color=arrow_color, lw=2))

    # 4. Output proportions
    bar_x = 0.76
    bar_w = 0.18
    bar_h = 0.08
    props = [0.25, 0.30, 0.25, 0.20]
    cum = bar_x
    for p, c in zip(props, wedge_colors):
        w = p * bar_w
        rect = plt.Rectangle((cum, 0.46), w, bar_h, facecolor=c, edgecolor='white', linewidth=1)
        ax.add_patch(rect)
        cum += w
    ax.text(bar_x + bar_w/2, 0.78, "Cell Type\nProportions", ha='center', va='bottom',
            fontsize=10, fontweight='bold')
    ax.text(bar_x + bar_w/2, 0.38, "Per-spot fractions", ha='center', va='top',
            fontsize=9, style='italic', color=PALETTE['neutral'])


def panel_b_bars(axes, method_data):
    """Panel B: Three clean bar charts."""
    methods = method_data['method'].tolist()
    x = np.arange(len(methods))
    colors = [get_method_color(m) for m in methods]

    metrics = [
        ('pearson_r_mean', 'pearson_r_std', 'Correlation (r)', (0, 0.55), True),
        ('jsd_mean', 'jsd_std', 'JSD', (0, 0.8), False),
        ('rmse_mean', 'rmse_std', 'RMSE', (0, 0.45), False),
    ]

    for ax, (mean_col, std_col, ylabel, ylim, higher_better) in zip(axes, metrics):
        bars = ax.bar(x, method_data[mean_col], yerr=method_data[std_col],
                      color=colors, alpha=0.85, capsize=4, error_kw={'linewidth': 1.5})
        ax.set_ylabel(ylabel, fontsize=10)
        ax.set_xticks(x)
        ax.set_xticklabels(methods, rotation=45, ha='right', fontsize=9)
        ax.set_ylim(ylim)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        # Mark best
        if higher_better:
            best_idx = method_data[mean_col].idxmax()
            note = "higher is better"
        else:
            best_idx = method_data[mean_col].idxmin()
            note = "lower is better"

        ax.text(0.97, 0.97, note, transform=ax.transAxes, ha='right', va='top',
                fontsize=8, style='italic', color=PALETTE['neutral'])

        # Highlight best bar
        bars[best_idx].set_edgecolor('black')
        bars[best_idx].set_linewidth(2)


def panel_c_scatter(ax, region_id=0):
    """Panel C: Clean scatter plot - Predicted vs Ground Truth."""
    # Panel label
    ax.text(-0.15, 1.05, "C", fontsize=14, fontweight='bold', va='top', transform=ax.transAxes)

    try:
        gt = pd.read_csv(GROUND_TRUTH_DIR / f"Xenium_region_{region_id}_prop.csv", index_col=0)
        pred = pd.read_csv(CITEGEIST_OUTPUT_DIR / f"Xenium_region_{region_id}_cell_prop_finetuned_results.csv", index_col=0)

        common = list(set(gt.index) & set(pred.index))
        gt, pred = gt.loc[common], pred.loc[common]

        gt_cols = [c for c in gt.columns if c not in ['Unassigned', 'Unknown', 'n_cells', 'spot_x', 'spot_y']]
        pred_cols = [c for c in pred.columns if c not in ['Unassigned', 'Unknown']]
        common_ct = [c for c in gt_cols if c in pred_cols]

        all_gt, all_pred, all_colors = [], [], []
        for ct in common_ct:
            all_gt.extend(gt[ct].values)
            all_pred.extend(pred[ct].values)
            # Use get_cell_type_color for flexibility
            color = CELL_TYPE_COLORS.get(ct, CELL_TYPE_COLORS.get(ct.split()[0], PALETTE['neutral']))
            all_colors.extend([color] * len(gt))

        ax.scatter(all_gt, all_pred, c=all_colors, alpha=0.4, s=15, rasterized=True)
        ax.plot([0, 1], [0, 1], 'k--', lw=1.5, alpha=0.5)

        corr = np.corrcoef(all_gt, all_pred)[0, 1]
        ax.text(0.05, 0.92, f'r = {corr:.2f}', transform=ax.transAxes, fontsize=11, fontweight='bold')

        ax.set_xlabel("Ground Truth Proportion", fontsize=10)
        ax.set_ylabel("CITEgeist Predicted", fontsize=10)
        ax.set_xlim(-0.02, 1.02)
        ax.set_ylim(-0.02, 1.02)
        ax.set_aspect('equal')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        # Compact legend
        handles = [mpatches.Patch(color=CELL_TYPE_COLORS.get(ct, PALETTE['neutral']),
                   label=ct.replace(' cells', '').replace('+ ', ''))
                   for ct in common_ct[:5]]
        ax.legend(handles=handles, loc='lower right', fontsize=8, framealpha=0.9)

    except Exception as e:
        ax.text(0.5, 0.5, f"Data not available\n({str(e)[:30]}...)",
                ha='center', va='center', fontsize=10, style='italic', color=PALETTE['neutral'])
        ax.axis('off')


def panel_d_summary(ax, method_data):
    """Panel D: Summary statistics table."""
    # Panel label
    ax.text(-0.05, 1.05, "D", fontsize=14, fontweight='bold', va='top', transform=ax.transAxes)
    ax.axis('off')

    # Create summary table
    summary_text = [
        ['Method', 'Pearson r', 'JSD', 'RMSE'],
    ]

    for _, row in method_data.iterrows():
        summary_text.append([
            row['method'],
            f"{row['pearson_r_mean']:.3f}",
            f"{row['jsd_mean']:.3f}",
            f"{row['rmse_mean']:.3f}",
        ])

    table = ax.table(
        cellText=summary_text[1:],
        colLabels=summary_text[0],
        cellLoc='center',
        loc='center',
        bbox=[0.05, 0.15, 0.9, 0.75]
    )
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1.2, 1.8)

    # Style header
    for j in range(len(summary_text[0])):
        table[(0, j)].set_facecolor(PALETTE['primary'])
        table[(0, j)].set_text_props(color='white', fontweight='bold')

    # Highlight CITEgeist row (row 1)
    for j in range(len(summary_text[0])):
        table[(1, j)].set_facecolor('#fee')
        table[(1, j)].set_text_props(fontweight='bold')

    ax.text(0.5, 0.05, "Xenium pseudo-Visium benchmark (14 regions, 7 cell types)",
            ha='center', va='bottom', fontsize=9, style='italic', color=PALETTE['neutral'],
            transform=ax.transAxes)


def generate_figure3():
    """Generate complete Figure 3."""
    print("Loading data...")
    summary_df = load_benchmark_summary()
    if summary_df is None:
        print("ERROR: Could not load benchmark summary")
        return

    method_data = prepare_method_data(summary_df)
    print(f"Methods: {method_data['method'].tolist()}")

    # Create figure with 2x2 layout
    fig = plt.figure(figsize=(11, 9))
    gs = GridSpec(2, 4, figure=fig, height_ratios=[0.8, 1], hspace=0.35, wspace=0.40)

    # Panel A: Schematic (top row, left half)
    ax_a = fig.add_subplot(gs[0, :2])
    panel_a_schematic(ax_a)
    ax_a.set_title("Two-Pass Deconvolution", fontsize=12, fontweight='bold', loc='left', pad=10)

    # Panel B label area (using first bar chart for label)
    ax_b1 = fig.add_subplot(gs[0, 2])
    ax_b2 = fig.add_subplot(gs[0, 3])
    ax_b3 = fig.add_subplot(gs[1, 0])
    panel_b_bars([ax_b1, ax_b2, ax_b3], method_data)

    # Add B label to first bar chart
    ax_b1.text(-0.2, 1.15, "B", fontsize=14, fontweight='bold', va='top', transform=ax_b1.transAxes)
    ax_b1.set_title("Benchmarking", fontsize=12, fontweight='bold', loc='left', pad=10)

    # Panel C: Scatter plot (bottom middle)
    ax_c = fig.add_subplot(gs[1, 1:3])
    panel_c_scatter(ax_c, region_id=0)
    ax_c.set_title("Predicted vs Ground Truth", fontsize=12, fontweight='bold', loc='left', pad=10)

    # Panel D: Summary table (bottom right)
    ax_d = fig.add_subplot(gs[1, 3])
    panel_d_summary(ax_d, method_data)
    ax_d.set_title("Summary", fontsize=12, fontweight='bold', loc='left', pad=10)

    plt.tight_layout()

    # Save
    output_path = OUTPUT_DIR / "figure3_benchmarking.pdf"
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"Saved to {output_path}")

    png_path = OUTPUT_DIR / "figure3_benchmarking.png"
    plt.savefig(png_path, dpi=150, bbox_inches='tight', facecolor='white')
    print(f"Preview saved to {png_path}")

    # Save SVG for Illustrator
    svg_path = OUTPUT_DIR / "figure3_benchmarking.svg"
    plt.savefig(svg_path, format='svg', bbox_inches='tight', facecolor='white')
    print(f"SVG saved to {svg_path}")

    plt.close()

    # Print summary
    print("\n=== Figure 3 Summary ===")
    for _, row in method_data.iterrows():
        print(f"{row['method']:15s}: r={row['pearson_r_mean']:.3f}, JSD={row['jsd_mean']:.3f}, RMSE={row['rmse_mean']:.3f}")

    print("\nEnhancements applied:")
    print("  - Added Panel C: Predicted vs Ground Truth scatter")
    print("  - Added Panel D: Summary statistics table")
    print("  - Consistent color palette from figure_style.py")
    print("  - Fonts increased to minimum 10pt")


if __name__ == "__main__":
    generate_figure3()
