#!/usr/bin/env python3
"""
Figure 3: Module 3 Deconvolution Benchmarking

Panels:
  A: Proportion estimation schematic (placeholder - to be created in vector editor)
  B: Benchmark comparison bar plots (JSD, RMSE, correlation)
     - CITEgeist vs Cell2Location, RCTD, Tangram, Seurat
  C: Spatial visualization of proportions (example region)

Data sources:
  - Benchmarking results: Benchmarking/xenium_benchmarking/evaluation/results/method_comparison/
  - Ground truth: Benchmarking/xenium_pseudovisium/data/ground_truth/
  - CITEgeist output: Benchmarking/xenium_benchmarking/CITEgeist/output_achievable_7/
  - Spatial data: Benchmarking/xenium_pseudovisium/data/h5ad_objects/
"""

import os
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
import seaborn as sns
from pathlib import Path

# Paths
PROJECT_ROOT = Path(__file__).parent.parent.parent
RESULTS_DIR = PROJECT_ROOT / "Benchmarking/xenium_benchmarking/evaluation/results/method_comparison"
GROUND_TRUTH_DIR = PROJECT_ROOT / "Benchmarking/xenium_pseudovisium/data/ground_truth"
CITEGEIST_OUTPUT_DIR = PROJECT_ROOT / "Benchmarking/xenium_benchmarking/CITEgeist/output_achievable_7"
H5AD_DIR = PROJECT_ROOT / "Benchmarking/xenium_pseudovisium/data/h5ad_objects"
OUTPUT_DIR = Path(__file__).parent / "output"

# Ensure output directory exists
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Style settings
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.size'] = 8
plt.rcParams['axes.linewidth'] = 0.5
plt.rcParams['xtick.major.width'] = 0.5
plt.rcParams['ytick.major.width'] = 0.5

# Method color palette - consistent across figures
METHOD_COLORS = {
    'CITEgeist': '#E41A1C',  # Red - primary method
    'Cell2Location': '#377EB8',  # Blue
    'RCTD': '#4DAF4A',  # Green
    'Tangram': '#984EA3',  # Purple
    'Seurat': '#FF7F00',  # Orange
}

# Cell type colors for spatial plots
CELL_TYPE_COLORS = {
    'B cells': '#E41A1C',
    'CD4+ T cells': '#377EB8',
    'CD8+ T cells': '#4DAF4A',
    'Macrophages': '#984EA3',
    'Endothelial': '#FF7F00',
    'Epithelial': '#FFFF33',
    'Fibroblasts': '#A65628',
}


def load_benchmark_summary():
    """Load method summary from CSV."""
    summary_path = RESULTS_DIR / "method_summary.csv"
    df = pd.read_csv(summary_path, header=[0, 1], index_col=0)

    # Flatten multi-index columns
    df.columns = ['_'.join(col).strip() for col in df.columns.values]
    return df


def load_comparison_table():
    """Load per-region comparison table."""
    table_path = RESULTS_DIR / "comparison_table.csv"
    df = pd.read_csv(table_path)
    return df


def load_full_results():
    """Load full results JSON for detailed per-cell-type metrics."""
    results_path = RESULTS_DIR / "full_results.json"
    with open(results_path, 'r') as f:
        return json.load(f)


def load_ground_truth_proportions(region_id=0):
    """Load ground truth proportions for a region."""
    gt_path = GROUND_TRUTH_DIR / f"Xenium_region_{region_id}_prop.csv"
    return pd.read_csv(gt_path, index_col=0)


def load_citegeist_proportions(region_id=0):
    """Load CITEgeist predicted proportions for a region."""
    pred_path = CITEGEIST_OUTPUT_DIR / f"Xenium_region_{region_id}_cell_prop_finetuned_results.csv"
    return pd.read_csv(pred_path, index_col=0)


def load_spatial_coordinates(region_id=0):
    """Load spatial coordinates from h5ad file."""
    try:
        import scanpy as sc
        h5ad_path = H5AD_DIR / f"Xenium_region_{region_id}_GEX.h5ad"
        adata = sc.read_h5ad(h5ad_path)
        coords = pd.DataFrame(
            adata.obsm['spatial'],
            index=adata.obs_names,
            columns=['x', 'y']
        )
        return coords
    except Exception as e:
        print(f"Warning: Could not load spatial coordinates: {e}")
        return None


def prepare_method_data(summary_df):
    """Prepare data for bar plots, extracting key metrics."""
    # Filter to main methods (achievable_7 suffix)
    methods = ['CITEgeist', 'Cell2Location', 'RCTD', 'Tangram', 'Seurat']
    method_mapping = {
        'CITEgeist_achievable_7': 'CITEgeist',
        'Cell2Location_achievable_7': 'Cell2Location',
        'RCTD_achievable_7': 'RCTD',
        'Tangram_achievable_7': 'Tangram',
        'Seurat_achievable_7': 'Seurat',
    }

    data = {
        'method': [],
        'pearson_r_mean': [],
        'pearson_r_std': [],
        'jsd_mean': [],
        'jsd_std': [],
        'rmse_mean': [],
        'rmse_std': [],
    }

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
    """Panel A: Placeholder for proportion estimation schematic."""
    ax.text(0.5, 0.5, "Panel A: Deconvolution Schematic\n(Create in vector editor)\n\n"
            "Show:\n- Spatial spots with mixed cell types\n- Reference profiles\n"
            "- Optimization (minimize reconstruction error)\n- Spatial regularization",
            ha='center', va='center', fontsize=9, style='italic',
            bbox=dict(boxstyle='round', facecolor='#f0f0f0', edgecolor='gray'))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis('off')
    ax.set_title("A. Module 3: Cell Proportion Estimation", fontsize=10, fontweight='bold', loc='left')


def panel_b_benchmark_bars(ax_pearson, ax_jsd, ax_rmse, method_data):
    """Panel B: Benchmark comparison bar plots."""
    methods = method_data['method'].tolist()
    x = np.arange(len(methods))
    width = 0.7

    colors = [METHOD_COLORS.get(m, '#808080') for m in methods]

    # Pearson correlation (higher is better)
    ax_pearson.bar(x, method_data['pearson_r_mean'], width,
                   yerr=method_data['pearson_r_std'],
                   color=colors, alpha=0.8, capsize=3, error_kw={'linewidth': 1})
    ax_pearson.set_ylabel("Pearson r", fontsize=8)
    ax_pearson.set_title("B. Proportion Accuracy", fontsize=10, fontweight='bold', loc='left')
    ax_pearson.set_xticks(x)
    ax_pearson.set_xticklabels(methods, rotation=45, ha='right', fontsize=7)
    ax_pearson.set_ylim(0, 0.6)
    ax_pearson.axhline(y=0, color='black', linewidth=0.5)
    ax_pearson.spines['top'].set_visible(False)
    ax_pearson.spines['right'].set_visible(False)

    # Add value labels
    for i, (v, std) in enumerate(zip(method_data['pearson_r_mean'], method_data['pearson_r_std'])):
        ax_pearson.text(i, v + std + 0.02, f'{v:.2f}', ha='center', va='bottom', fontsize=6)

    # JSD (lower is better)
    ax_jsd.bar(x, method_data['jsd_mean'], width,
               yerr=method_data['jsd_std'],
               color=colors, alpha=0.8, capsize=3, error_kw={'linewidth': 1})
    ax_jsd.set_ylabel("JSD", fontsize=8)
    ax_jsd.set_title("Jensen-Shannon Divergence", fontsize=9, loc='left')
    ax_jsd.set_xticks(x)
    ax_jsd.set_xticklabels(methods, rotation=45, ha='right', fontsize=7)
    ax_jsd.set_ylim(0, 0.8)
    ax_jsd.spines['top'].set_visible(False)
    ax_jsd.spines['right'].set_visible(False)

    # Add "lower is better" annotation
    ax_jsd.text(0.95, 0.95, r'$\downarrow$ better', transform=ax_jsd.transAxes,
                ha='right', va='top', fontsize=7, style='italic')

    # RMSE (lower is better)
    ax_rmse.bar(x, method_data['rmse_mean'], width,
                yerr=method_data['rmse_std'],
                color=colors, alpha=0.8, capsize=3, error_kw={'linewidth': 1})
    ax_rmse.set_ylabel("RMSE", fontsize=8)
    ax_rmse.set_title("Root Mean Square Error", fontsize=9, loc='left')
    ax_rmse.set_xticks(x)
    ax_rmse.set_xticklabels(methods, rotation=45, ha='right', fontsize=7)
    ax_rmse.set_ylim(0, 0.5)
    ax_rmse.spines['top'].set_visible(False)
    ax_rmse.spines['right'].set_visible(False)

    # Add "lower is better" annotation
    ax_rmse.text(0.95, 0.95, r'$\downarrow$ better', transform=ax_rmse.transAxes,
                 ha='right', va='top', fontsize=7, style='italic')


def panel_c_spatial_visualization(axes, region_id=0):
    """Panel C: Spatial visualization of proportions for example region."""
    # Load data
    gt = load_ground_truth_proportions(region_id)
    pred = load_citegeist_proportions(region_id)
    coords = load_spatial_coordinates(region_id)

    if coords is None:
        # Fallback: create placeholder
        for ax in axes:
            ax.text(0.5, 0.5, "Spatial data\nnot available",
                    ha='center', va='center', fontsize=9, style='italic')
            ax.axis('off')
        return

    # Get common spots
    common_spots = list(set(gt.index) & set(pred.index) & set(coords.index))
    common_spots.sort()

    gt = gt.loc[common_spots]
    pred = pred.loc[common_spots]
    coords = coords.loc[common_spots]

    # Select 3 cell types to visualize
    cell_types_to_show = ['Macrophages', 'Epithelial', 'Fibroblasts']

    # Find common cell types
    gt_celltypes = [c for c in gt.columns if c in cell_types_to_show]
    pred_celltypes = [c for c in pred.columns if c in cell_types_to_show]
    common_celltypes = [c for c in gt_celltypes if c in pred_celltypes][:3]

    if len(common_celltypes) < 3:
        # Use whatever cell types are available
        common_celltypes = [c for c in gt.columns if c in pred.columns and c not in ['Unassigned', 'Unknown']][:3]

    for i, (ax, ct) in enumerate(zip(axes, common_celltypes)):
        if ct in gt.columns and ct in pred.columns:
            # Plot ground truth as background scatter
            gt_vals = gt[ct].values
            pred_vals = pred[ct].values

            # Normalize for visualization
            vmin, vmax = 0, max(gt_vals.max(), pred_vals.max())

            # Plot predicted proportions
            sc = ax.scatter(coords['x'], coords['y'], c=pred_vals,
                           cmap='Reds', s=3, vmin=vmin, vmax=vmax, alpha=0.8)

            ax.set_aspect('equal')
            ax.set_title(ct, fontsize=8)
            ax.set_xticks([])
            ax.set_yticks([])
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['bottom'].set_visible(False)
            ax.spines['left'].set_visible(False)

            # Add colorbar for last plot
            if i == len(common_celltypes) - 1:
                cbar = plt.colorbar(sc, ax=ax, shrink=0.8)
                cbar.set_label('Proportion', fontsize=7)
                cbar.ax.tick_params(labelsize=6)


def panel_c_scatter_comparison(ax, region_id=0):
    """Alternative Panel C: Scatter plot of predicted vs ground truth."""
    # Load data
    gt = load_ground_truth_proportions(region_id)
    pred = load_citegeist_proportions(region_id)

    # Get common spots and cell types
    common_spots = list(set(gt.index) & set(pred.index))
    common_spots.sort()

    gt = gt.loc[common_spots]
    pred = pred.loc[common_spots]

    # Get common cell types (excluding Unassigned/Unknown)
    gt_celltypes = [c for c in gt.columns if c not in ['Unassigned', 'Unknown', 'n_cells', 'spot_x', 'spot_y']]
    pred_celltypes = [c for c in pred.columns if c not in ['Unassigned', 'Unknown']]
    common_celltypes = [c for c in gt_celltypes if c in pred_celltypes]

    # Collect all points for scatter
    all_gt = []
    all_pred = []
    all_colors = []

    for ct in common_celltypes:
        gt_vals = gt[ct].values
        pred_vals = pred[ct].values
        all_gt.extend(gt_vals)
        all_pred.extend(pred_vals)
        color = CELL_TYPE_COLORS.get(ct, '#808080')
        all_colors.extend([color] * len(gt_vals))

    # Plot scatter
    ax.scatter(all_gt, all_pred, c=all_colors, alpha=0.3, s=5, rasterized=True)

    # Add diagonal line
    ax.plot([0, 1], [0, 1], 'k--', linewidth=1, alpha=0.5)

    # Calculate correlation
    corr = np.corrcoef(all_gt, all_pred)[0, 1]

    ax.set_xlabel("Ground Truth Proportion", fontsize=8)
    ax.set_ylabel("Predicted Proportion", fontsize=8)
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    ax.set_aspect('equal')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Add correlation text
    ax.text(0.05, 0.95, f'r = {corr:.3f}', transform=ax.transAxes,
            fontsize=8, va='top', fontweight='bold')

    # Add legend
    handles = [mpatches.Patch(color=CELL_TYPE_COLORS.get(ct, '#808080'), label=ct)
               for ct in common_celltypes]
    ax.legend(handles=handles, loc='lower right', fontsize=6, framealpha=0.9)

    ax.set_title("C. Predicted vs Ground Truth (Region 0)", fontsize=10, fontweight='bold', loc='left')


def generate_figure3():
    """Generate complete Figure 3."""
    print("Loading data...")

    # Load benchmark data
    summary_df = load_benchmark_summary()
    method_data = prepare_method_data(summary_df)

    print(f"Methods loaded: {method_data['method'].tolist()}")
    print(f"Mean Pearson r: {dict(zip(method_data['method'], method_data['pearson_r_mean'].round(3)))}")

    # Create figure
    fig = plt.figure(figsize=(11, 9))

    # Layout:
    # Row 1: Panel A (schematic) - spans full width
    # Row 2: Panel B (3 bar plots side by side)
    # Row 3: Panel C (scatter comparison + spatial plots)

    gs = GridSpec(3, 4, figure=fig, height_ratios=[1, 1.2, 1.2],
                  hspace=0.4, wspace=0.35)

    # Panel A: Schematic (placeholder) - top row
    ax_a = fig.add_subplot(gs[0, :])
    panel_a_schematic(ax_a)

    # Panel B: Benchmark bar plots - middle row
    ax_pearson = fig.add_subplot(gs[1, 0:2])
    ax_jsd = fig.add_subplot(gs[1, 2])
    ax_rmse = fig.add_subplot(gs[1, 3])
    panel_b_benchmark_bars(ax_pearson, ax_jsd, ax_rmse, method_data)

    # Panel C: Scatter comparison - bottom left
    ax_scatter = fig.add_subplot(gs[2, 0:2])
    panel_c_scatter_comparison(ax_scatter, region_id=0)

    # Panel C (continued): Spatial visualization - bottom right
    ax_spatial1 = fig.add_subplot(gs[2, 2])
    ax_spatial2 = fig.add_subplot(gs[2, 3])

    # Create a mini figure for spatial plots
    try:
        # Load spatial data
        gt = load_ground_truth_proportions(0)
        pred = load_citegeist_proportions(0)
        coords = load_spatial_coordinates(0)

        if coords is not None:
            common_spots = list(set(gt.index) & set(pred.index) & set(coords.index))
            common_spots.sort()

            gt = gt.loc[common_spots]
            pred = pred.loc[common_spots]
            coords = coords.loc[common_spots]

            # Show 2 cell types: ground truth vs predicted
            ct = 'Macrophages'
            if ct in gt.columns and ct in pred.columns:
                gt_vals = gt[ct].values
                pred_vals = pred[ct].values
                vmax = max(gt_vals.max(), pred_vals.max())

                # Ground truth
                sc1 = ax_spatial1.scatter(coords['x'], coords['y'], c=gt_vals,
                                         cmap='Reds', s=3, vmin=0, vmax=vmax, alpha=0.8)
                ax_spatial1.set_title(f'{ct}\n(Ground Truth)', fontsize=8)
                ax_spatial1.set_aspect('equal')
                ax_spatial1.axis('off')

                # Predicted
                sc2 = ax_spatial2.scatter(coords['x'], coords['y'], c=pred_vals,
                                         cmap='Reds', s=3, vmin=0, vmax=vmax, alpha=0.8)
                ax_spatial2.set_title(f'{ct}\n(CITEgeist)', fontsize=8)
                ax_spatial2.set_aspect('equal')
                ax_spatial2.axis('off')

                # Add colorbar
                cbar = plt.colorbar(sc2, ax=ax_spatial2, shrink=0.6)
                cbar.set_label('Proportion', fontsize=7)
                cbar.ax.tick_params(labelsize=6)
        else:
            raise Exception("No spatial coordinates")

    except Exception as e:
        print(f"Warning: Could not create spatial plots: {e}")
        ax_spatial1.text(0.5, 0.5, "Spatial\nvisualization\n(requires\nscanpy)",
                        ha='center', va='center', fontsize=8, style='italic')
        ax_spatial1.axis('off')
        ax_spatial2.text(0.5, 0.5, "Predicted\nproportions",
                        ha='center', va='center', fontsize=8, style='italic')
        ax_spatial2.axis('off')

    # Save
    output_path = OUTPUT_DIR / "figure3_benchmarking.pdf"
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"Saved to {output_path}")

    # Also save PNG for quick preview
    png_path = OUTPUT_DIR / "figure3_benchmarking.png"
    plt.savefig(png_path, dpi=150, bbox_inches='tight', facecolor='white')
    print(f"Preview saved to {png_path}")

    plt.close()

    # Print summary statistics
    print("\n=== Figure 3 Summary Statistics ===")
    print("\nMethod Performance Summary (7 cell types, 5 regions):")
    print("-" * 60)
    for _, row in method_data.iterrows():
        print(f"{row['method']:15s}: r={row['pearson_r_mean']:.3f} +/- {row['pearson_r_std']:.3f}, "
              f"JSD={row['jsd_mean']:.3f}, RMSE={row['rmse_mean']:.3f}")

    # Highlight CITEgeist performance
    citegeist = method_data[method_data['method'] == 'CITEgeist'].iloc[0]
    best_r = method_data['pearson_r_mean'].max()
    best_jsd = method_data['jsd_mean'].min()
    best_rmse = method_data['rmse_mean'].min()

    print("\n" + "=" * 60)
    print("Key Findings:")
    if citegeist['pearson_r_mean'] == best_r:
        print(f"  - CITEgeist achieves BEST correlation (r={best_r:.3f})")
    if citegeist['jsd_mean'] == best_jsd:
        print(f"  - CITEgeist achieves BEST JSD ({best_jsd:.3f})")
    if citegeist['rmse_mean'] == best_rmse:
        print(f"  - CITEgeist achieves BEST RMSE ({best_rmse:.3f})")

    # Relative improvement
    cell2loc = method_data[method_data['method'] == 'Cell2Location'].iloc[0]
    r_improvement = ((citegeist['pearson_r_mean'] - cell2loc['pearson_r_mean']) /
                     cell2loc['pearson_r_mean'] * 100)
    jsd_improvement = ((cell2loc['jsd_mean'] - citegeist['jsd_mean']) /
                       cell2loc['jsd_mean'] * 100)
    print(f"  - CITEgeist vs Cell2Location: {r_improvement:.1f}% better correlation, "
          f"{jsd_improvement:.1f}% lower JSD")


if __name__ == "__main__":
    generate_figure3()
