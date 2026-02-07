"""
Generate publication-quality benchmark figures:
  1. benchmark_schematic.png - Flowchart of the Xenium benchmarking design
  2. benchmark_spatial_comparison.png - GT vs CITEgeist pie chart spatial maps (Region 4)
  3. benchmark_citegeist_vs_scresolve.png - GT | CITEgeist | scResolve spatial pies
     + performance metric bar charts (Region 0)
"""

import os
import numpy as np
import pandas as pd
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from matplotlib.collections import PatchCollection

# ── Paths ──────────────────────────────────────────────────────────────────────
BASE = os.path.dirname(os.path.abspath(__file__))
BENCH_DIR = os.path.dirname(BASE)  # xenium_benchmarking/
PSEUDO_DIR = os.path.join(os.path.dirname(BENCH_DIR), "xenium_pseudovisium")

GT_CSV = os.path.join(PSEUDO_DIR, "data_granular_gt", "ground_truth", "Xenium_region_4_prop.csv")
PRED_CSV = os.path.join(BENCH_DIR, "CITEgeist", "output_achievable_7",
                        "Xenium_region_4_cell_prop_global_results.csv")
H5AD_PATH = os.path.join(PSEUDO_DIR, "data", "h5ad_objects", "Xenium_region_4_GEX.h5ad")

OUTPUT_DIR = BASE  # figures/ directory

# ── Achievable-7 mapping ──────────────────────────────────────────────────────
GT_TO_ACHIEVABLE_7 = {
    "B cells": "B cells",
    "Mixed Immune": "CD4+ T cells",
    "CD8+ T cells": "CD8+ T cells",
    "Proliferating T": "CD8+ T cells",
    "Macrophages": "Macrophages",
    "Endothelial": "Endothelial",
    "Vascular Stromal": "Endothelial",
    "Epithelial": "Epithelial",
    "Myofibroblasts": "Fibroblasts",
    "Stromal": "Fibroblasts",
}

ACHIEVABLE_7 = [
    "B cells", "CD4+ T cells", "CD8+ T cells",
    "Macrophages", "Endothelial", "Epithelial", "Fibroblasts",
]

# Colorblind-friendly palette (7 distinct colors)
CELL_COLORS = {
    "B cells":       "#1b9e77",  # teal
    "CD4+ T cells":  "#d95f02",  # orange
    "CD8+ T cells":  "#7570b3",  # purple
    "Macrophages":   "#e7298a",  # pink
    "Endothelial":   "#66a61e",  # green
    "Epithelial":    "#e6ab02",  # gold
    "Fibroblasts":   "#a6761d",  # brown
}


# ── Data loading ──────────────────────────────────────────────────────────────

def load_spatial_coords(h5ad_path):
    """Load spot names and spatial coordinates from h5ad file."""
    with h5py.File(h5ad_path, 'r') as f:
        coords = f['obsm']['spatial'][:]
        names = [x.decode() if isinstance(x, bytes) else x for x in f['obs']['_index'][:]]
    return pd.DataFrame(coords, index=names, columns=['x', 'y'])


def collapse_gt_to_achievable_7(gt_df):
    """Collapse 10 granular cell types to achievable-7."""
    prop_cols = [c for c in gt_df.columns if c in GT_TO_ACHIEVABLE_7]
    collapsed = pd.DataFrame(0.0, index=gt_df.index, columns=ACHIEVABLE_7)
    for src, dst in GT_TO_ACHIEVABLE_7.items():
        if src in gt_df.columns:
            collapsed[dst] += gt_df[src]
    # Renormalize rows to sum to 1
    row_sums = collapsed.sum(axis=1)
    row_sums = row_sums.replace(0, 1)  # avoid div by zero
    collapsed = collapsed.div(row_sums, axis=0)
    return collapsed


def load_data():
    """Load GT, predictions, and coordinates for Region 4."""
    coords = load_spatial_coords(H5AD_PATH)

    # Ground truth
    gt_raw = pd.read_csv(GT_CSV, index_col=0)
    gt = collapse_gt_to_achievable_7(gt_raw)

    # CITEgeist predictions
    pred = pd.read_csv(PRED_CSV, index_col=0)
    # Drop "Unknown" column if present, renormalize
    if "Unknown" in pred.columns:
        pred = pred.drop(columns=["Unknown"])
    pred = pred[ACHIEVABLE_7]
    row_sums = pred.sum(axis=1)
    row_sums = row_sums.replace(0, 1)
    pred = pred.div(row_sums, axis=0)

    # Align all dataframes to common spots
    common = coords.index.intersection(gt.index).intersection(pred.index)
    print(f"Common spots: {len(common)}")
    return coords.loc[common], gt.loc[common], pred.loc[common]


# ── Figure 1: Benchmark Schematic ─────────────────────────────────────────────

def draw_schematic(output_path):
    """Draw horizontal flowchart of benchmark design."""
    fig, ax = plt.subplots(figsize=(14, 3.5))
    ax.set_xlim(-0.2, 10.8)
    ax.set_ylim(-1.2, 2.2)
    ax.axis('off')

    # Box positions (x_center, y_center)
    boxes = [
        (0.8, 0.7, "10x Xenium\nIn Situ", "#E8F4FD"),
        (3.0, 0.7, "Pseudo-Visium\nConversion", "#FFF3E0"),
        (5.2, 0.7, "Ground Truth\nDefinition", "#E8F5E9"),
        (7.4, 0.7, "Deconvolution\nMethods", "#F3E5F5"),
        (9.6, 0.7, "Evaluation\nMetrics", "#FFEBEE"),
    ]

    subtexts = [
        "Single-cell\nresolution RNA",
        "Aggregate into\n~55µm spots\n(1,414 spots)",
        "Achievable-7\ncell types from\nsingle-cell labels",
        "CITEgeist\nCell2Location\nRCTD\nTangram\nSeurat",
        "Pearson r\nRMSE, MAE\nJSD",
    ]

    box_width = 1.6
    box_height = 0.8

    for i, (cx, cy, title, color) in enumerate(boxes):
        rect = FancyBboxPatch(
            (cx - box_width/2, cy - box_height/2),
            box_width, box_height,
            boxstyle="round,pad=0.1",
            facecolor=color, edgecolor="#333333", linewidth=1.5
        )
        ax.add_patch(rect)
        ax.text(cx, cy, title, fontsize=11, fontweight='bold',
                ha='center', va='center', fontfamily='sans-serif')

        # Subtitle below
        ax.text(cx, cy - box_height/2 - 0.25, subtexts[i],
                fontsize=8.5, ha='center', va='top', color='#555555',
                fontfamily='sans-serif', linespacing=1.3)

        # Arrow to next box
        if i < len(boxes) - 1:
            next_cx = boxes[i + 1][0]
            ax.annotate(
                '', xy=(next_cx - box_width/2 - 0.05, cy),
                xytext=(cx + box_width/2 + 0.05, cy),
                arrowprops=dict(arrowstyle='->', color='#333333',
                                lw=2, connectionstyle='arc3,rad=0')
            )

    # Title
    ax.text(5.2, 1.9, "Xenium Spatial Deconvolution Benchmark",
            fontsize=15, fontweight='bold', ha='center', va='center',
            fontfamily='sans-serif')

    fig.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f"Saved: {output_path}")


# ── Figure 2: Spatial Pie Chart Comparison ────────────────────────────────────

def draw_pie_at(ax, x, y, proportions, colors, radius=30):
    """Draw a single pie chart at (x, y) on the axes."""
    # Filter out near-zero proportions to reduce wedge count
    mask = proportions > 0.01
    if mask.sum() == 0:
        mask[proportions.argmax()] = True
    props = proportions[mask]
    cols = [colors[i] for i in range(len(proportions)) if mask[i]]

    # Convert proportions to angles
    angles = props / props.sum() * 360.0
    start_angle = 90  # start from top

    wedges = []
    for angle, color in zip(angles, cols):
        theta1 = start_angle
        theta2 = start_angle + angle
        wedge = mpatches.Wedge(
            (x, y), radius, theta1, theta2,
            facecolor=color, edgecolor='none', linewidth=0
        )
        ax.add_patch(wedge)
        start_angle = theta2


def compute_pie_radius(x, y):
    """Compute pie radius from nearest-neighbor distances."""
    from scipy.spatial import distance
    sample_idx = np.random.choice(len(x), min(200, len(x)), replace=False)
    sample_coords = np.column_stack([x[sample_idx], y[sample_idx]])
    nn_dists = []
    for i in range(len(sample_coords)):
        d = distance.cdist([sample_coords[i]], sample_coords)[0]
        d[i] = np.inf
        nn_dists.append(d.min())
    return np.median(nn_dists) * 0.40


def draw_spatial_panel(ax, coords, proportions_df, title, radius, xlim, ylim):
    """Draw spatial pie chart map on given axes with fixed limits."""
    x = coords['x'].values
    y = coords['y'].values

    color_list = [CELL_COLORS[ct] for ct in ACHIEVABLE_7]
    prop_values = proportions_df[ACHIEVABLE_7].values

    for i in range(len(x)):
        draw_pie_at(ax, x[i], y[i], prop_values[i], color_list, radius=radius)

    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_aspect('equal')
    ax.set_title(title, fontsize=14, fontweight='bold', fontfamily='sans-serif', pad=10)
    ax.axis('off')


def draw_spatial_comparison(output_path, coords, gt, pred):
    """Draw side-by-side GT vs CITEgeist spatial pie chart maps."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 10))

    x = coords['x'].values
    y = coords['y'].values

    # Compute shared radius and axis limits
    radius = compute_pie_radius(x, y)
    pad = radius * 3
    xlim = (x.min() - pad, x.max() + pad)
    ylim = (y.min() - pad, y.max() + pad)

    print("Drawing ground truth panel...")
    draw_spatial_panel(ax1, coords, gt, "Ground Truth", radius, xlim, ylim)
    print("Drawing CITEgeist panel...")
    draw_spatial_panel(ax2, coords, pred, "CITEgeist (Global)", radius, xlim, ylim)

    # Shared legend
    legend_handles = [
        mpatches.Patch(facecolor=CELL_COLORS[ct], edgecolor='#333333',
                       linewidth=0.5, label=ct)
        for ct in ACHIEVABLE_7
    ]
    fig.legend(handles=legend_handles, loc='lower center', ncol=7,
               fontsize=11, frameon=True, edgecolor='#cccccc',
               bbox_to_anchor=(0.5, -0.01))

    fig.suptitle("Region 4: Ground Truth vs CITEgeist Predictions",
                 fontsize=16, fontweight='bold', fontfamily='sans-serif', y=0.98)

    plt.tight_layout(rect=[0, 0.03, 1, 0.96])
    fig.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f"Saved: {output_path}")


# ── Figure 3: CITEgeist vs scResolve (Region 0) ──────────────────────────────

# Region 0 paths
GT_CSV_R0 = os.path.join(PSEUDO_DIR, "data_granular_gt", "ground_truth", "Xenium_region_0_prop.csv")
PRED_CSV_R0 = os.path.join(BENCH_DIR, "CITEgeist", "output_achievable_7",
                           "Xenium_region_0_cell_prop_global_results.csv")
SCRESOLVE_CSV_R0 = os.path.join(BENCH_DIR, "scResolve", "output_morphology",
                                "Xenium_region_0", "Xenium_region_0_fair_proportions.csv")
H5AD_PATH_R0 = os.path.join(PSEUDO_DIR, "data", "h5ad_objects", "Xenium_region_0_GEX.h5ad")


def load_predictions(csv_path):
    """Load a prediction CSV and normalize to achievable-7."""
    pred = pd.read_csv(csv_path, index_col=0)
    if "Unknown" in pred.columns:
        pred = pred.drop(columns=["Unknown"])
    pred = pred[ACHIEVABLE_7]
    row_sums = pred.sum(axis=1).replace(0, 1)
    return pred.div(row_sums, axis=0)


def compute_metrics(gt, pred):
    """Compute per-cell-type and overall metrics between GT and prediction."""
    from scipy.spatial.distance import jensenshannon
    from scipy.stats import pearsonr

    metrics = {}

    # Per-cell-type Pearson r
    per_ct_r = {}
    for ct in ACHIEVABLE_7:
        r, _ = pearsonr(gt[ct].values, pred[ct].values)
        per_ct_r[ct] = r
    metrics['per_ct_pearson'] = per_ct_r

    # Overall metrics (flatten all cell types)
    gt_flat = gt[ACHIEVABLE_7].values.flatten()
    pred_flat = pred[ACHIEVABLE_7].values.flatten()
    r_overall, _ = pearsonr(gt_flat, pred_flat)
    metrics['pearson_r'] = r_overall

    # RMSE
    metrics['rmse'] = np.sqrt(np.mean((gt[ACHIEVABLE_7].values - pred[ACHIEVABLE_7].values) ** 2))

    # MAE
    metrics['mae'] = np.mean(np.abs(gt[ACHIEVABLE_7].values - pred[ACHIEVABLE_7].values))

    # Mean JSD per spot
    jsd_vals = []
    for i in range(len(gt)):
        g = gt[ACHIEVABLE_7].iloc[i].values
        p = pred[ACHIEVABLE_7].iloc[i].values
        # Ensure non-negative and sum > 0
        g = np.clip(g, 0, None)
        p = np.clip(p, 0, None)
        if g.sum() > 0 and p.sum() > 0:
            jsd_vals.append(jensenshannon(g, p) ** 2)  # squared JSD
    metrics['jsd'] = np.mean(jsd_vals)

    return metrics


def load_region0_data():
    """Load GT, CITEgeist, scResolve, and coordinates for Region 0."""
    coords = load_spatial_coords(H5AD_PATH_R0)

    gt_raw = pd.read_csv(GT_CSV_R0, index_col=0)
    gt = collapse_gt_to_achievable_7(gt_raw)

    citegeist = load_predictions(PRED_CSV_R0)
    scresolve = load_predictions(SCRESOLVE_CSV_R0)

    # Align to common spots
    common = coords.index.intersection(gt.index).intersection(
        citegeist.index).intersection(scresolve.index)
    print(f"Region 0 common spots: {len(common)}")
    return (coords.loc[common], gt.loc[common],
            citegeist.loc[common], scresolve.loc[common])


def draw_metric_bars(ax, cg_metrics, sr_metrics, metric_key, ylabel, title):
    """Draw grouped bar chart comparing CITEgeist vs scResolve for a metric."""
    if metric_key == 'per_ct_pearson':
        cell_types = ACHIEVABLE_7
        cg_vals = [cg_metrics['per_ct_pearson'][ct] for ct in cell_types]
        sr_vals = [sr_metrics['per_ct_pearson'][ct] for ct in cell_types]
        x = np.arange(len(cell_types))
        width = 0.35
        ax.bar(x - width/2, cg_vals, width, label='CITEgeist',
               color='#4A90D9', edgecolor='#333333', linewidth=0.5)
        ax.bar(x + width/2, sr_vals, width, label='scResolve',
               color='#E85D75', edgecolor='#333333', linewidth=0.5)
        ax.set_xticks(x)
        ax.set_xticklabels([ct.replace(' ', '\n') for ct in cell_types],
                           fontsize=8, fontfamily='sans-serif')
        ax.axhline(y=0, color='#999999', linewidth=0.5, linestyle='--')
    else:
        methods = ['CITEgeist', 'scResolve']
        vals = [cg_metrics[metric_key], sr_metrics[metric_key]]
        colors = ['#4A90D9', '#E85D75']
        bars = ax.bar(methods, vals, color=colors, edgecolor='#333333',
                      linewidth=0.5, width=0.5)
        # Add value labels on bars
        for bar, val in zip(bars, vals):
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.005,
                    f'{val:.3f}', ha='center', va='bottom', fontsize=9,
                    fontfamily='sans-serif')

    ax.set_ylabel(ylabel, fontsize=10, fontfamily='sans-serif')
    ax.set_title(title, fontsize=11, fontweight='bold', fontfamily='sans-serif')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)


def draw_citegeist_vs_scresolve(output_path, coords, gt, citegeist, scresolve):
    """Draw 3-panel spatial comparison + metric bar charts."""
    # Layout: top row = 3 spatial panels, bottom row = 4 metric charts
    fig = plt.figure(figsize=(22, 16))

    # GridSpec: 2 rows, top row taller
    gs = fig.add_gridspec(2, 3, height_ratios=[2.5, 1], hspace=0.25, wspace=0.15)

    # Spatial panels (top row)
    ax_gt = fig.add_subplot(gs[0, 0])
    ax_cg = fig.add_subplot(gs[0, 1])
    ax_sr = fig.add_subplot(gs[0, 2])

    x = coords['x'].values
    y = coords['y'].values
    radius = compute_pie_radius(x, y)
    pad = radius * 3
    xlim = (x.min() - pad, x.max() + pad)
    ylim = (y.min() - pad, y.max() + pad)

    print("Drawing GT panel (Region 0)...")
    draw_spatial_panel(ax_gt, coords, gt, "Ground Truth", radius, xlim, ylim)
    print("Drawing CITEgeist panel (Region 0)...")
    draw_spatial_panel(ax_cg, coords, citegeist, "CITEgeist", radius, xlim, ylim)
    print("Drawing scResolve panel (Region 0)...")
    draw_spatial_panel(ax_sr, coords, scresolve, "scResolve", radius, xlim, ylim)

    # Compute metrics
    print("Computing metrics...")
    cg_metrics = compute_metrics(gt, citegeist)
    sr_metrics = compute_metrics(gt, scresolve)

    print(f"CITEgeist  - Pearson r: {cg_metrics['pearson_r']:.3f}, "
          f"RMSE: {cg_metrics['rmse']:.3f}, MAE: {cg_metrics['mae']:.3f}, "
          f"JSD: {cg_metrics['jsd']:.3f}")
    print(f"scResolve  - Pearson r: {sr_metrics['pearson_r']:.3f}, "
          f"RMSE: {sr_metrics['rmse']:.3f}, MAE: {sr_metrics['mae']:.3f}, "
          f"JSD: {sr_metrics['jsd']:.3f}")

    # Metric bar charts (bottom row)
    # Split bottom row into 4 panels using nested gridspec
    gs_bottom = gs[1, :].subgridspec(1, 4, wspace=0.4)

    ax_pr = fig.add_subplot(gs_bottom[0])
    draw_metric_bars(ax_pr, cg_metrics, sr_metrics, 'pearson_r',
                     'Pearson r', 'Overall Pearson r')

    ax_rmse = fig.add_subplot(gs_bottom[1])
    draw_metric_bars(ax_rmse, cg_metrics, sr_metrics, 'rmse',
                     'RMSE', 'RMSE')

    ax_mae = fig.add_subplot(gs_bottom[2])
    draw_metric_bars(ax_mae, cg_metrics, sr_metrics, 'mae',
                     'MAE', 'MAE')

    ax_jsd = fig.add_subplot(gs_bottom[3])
    draw_metric_bars(ax_jsd, cg_metrics, sr_metrics, 'jsd',
                     'JSD', 'Jensen-Shannon Divergence')

    # Shared legend for spatial panels
    legend_handles = [
        mpatches.Patch(facecolor=CELL_COLORS[ct], edgecolor='#333333',
                       linewidth=0.5, label=ct)
        for ct in ACHIEVABLE_7
    ]
    fig.legend(handles=legend_handles, loc='upper center', ncol=7,
               fontsize=11, frameon=True, edgecolor='#cccccc',
               bbox_to_anchor=(0.5, 0.52))

    fig.suptitle("Region 0: CITEgeist vs scResolve",
                 fontsize=18, fontweight='bold', fontfamily='sans-serif', y=0.98)

    fig.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f"Saved: {output_path}")


# ── Main ──────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    # Figure 1: Schematic
    schematic_path = os.path.join(OUTPUT_DIR, "benchmark_schematic.png")
    draw_schematic(schematic_path)

    # Figure 2: Spatial comparison (Region 4)
    print("Loading Region 4 data...")
    coords, gt, pred = load_data()
    spatial_path = os.path.join(OUTPUT_DIR, "benchmark_spatial_comparison.png")
    draw_spatial_comparison(spatial_path, coords, gt, pred)

    # Figure 3: CITEgeist vs scResolve (Region 0)
    print("\nLoading Region 0 data...")
    coords_r0, gt_r0, cg_r0, sr_r0 = load_region0_data()
    vs_path = os.path.join(OUTPUT_DIR, "benchmark_citegeist_vs_scresolve.png")
    draw_citegeist_vs_scresolve(vs_path, coords_r0, gt_r0, cg_r0, sr_r0)

    print("\nAll figures done!")
