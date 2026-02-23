#!/usr/bin/env python3
"""
Generate Figure 3F spatial comparison panel with αSMA-only fibroblast GT.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
from pathlib import Path

plt.rcParams.update({
    'font.size': 10,
    'axes.titlesize': 11,
    'font.family': 'sans-serif',
    'figure.dpi': 150,
    'savefig.dpi': 300,
})

PROJECT_ROOT = Path(__file__).parent.parent.parent
GT_DIR = PROJECT_ROOT / "Benchmarking/xenium_pseudovisium/data_protein_gt/ground_truth"
H5AD_DIR = PROJECT_ROOT / "Benchmarking/xenium_pseudovisium/data_protein_gt/h5ad_objects"
PRED_DIR = PROJECT_ROOT / "Benchmarking/xenium_benchmarking/CITEgeist/output/manual"
OUTPUT_DIR = Path(__file__).parent / "output"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def load_data(region_id=0):
    """Load GT, predictions, and coordinates."""
    # Ground truth
    gt_path = GT_DIR / f"Xenium_region_{region_id}_prop.csv"
    gt = pd.read_csv(gt_path, index_col=0)

    # Predictions - check both finetuned and global (files are at top level, not in subdirs)
    pred_path = PRED_DIR / f"Xenium_region_{region_id}_cell_prop_finetuned_results.csv"
    if not pred_path.exists():
        pred_path = PRED_DIR / f"Xenium_region_{region_id}_cell_prop_global_results.csv"

    if not pred_path.exists():
        raise FileNotFoundError(f"No prediction file found for region {region_id} at {pred_path}")

    pred = pd.read_csv(pred_path, index_col=0)

    # Coordinates from h5ad file (obsm['spatial'])
    h5ad_path = H5AD_DIR / f"Xenium_region_{region_id}_CITE.h5ad"
    adata = sc.read_h5ad(h5ad_path)
    coords = pd.DataFrame(
        adata.obsm['spatial'],
        index=adata.obs_names,
        columns=['x', 'y']
    )

    return gt, pred, coords


def generate_spatial_panel(region_id=0):
    """Generate the spatial comparison panel."""
    print(f"Loading data for region {region_id}...")
    gt, pred, coords = load_data(region_id)

    # Align indices
    common_spots = sorted(set(gt.index) & set(pred.index))
    gt = gt.loc[common_spots]
    pred = pred.loc[common_spots]
    coords = coords.loc[common_spots]

    print(f"  Spots: {len(common_spots)}")
    print(f"  GT columns: {list(gt.columns)}")
    print(f"  Pred columns: {list(pred.columns)}")

    # Cell types to show - direct matches
    cell_types = ["B cells", "Macrophages", "Endothelial", "Epithelial", "Fibroblasts"]

    # Check for T cells (may need combining in GT)
    if "T cells" in pred.columns:
        cell_types.insert(1, "T cells")

    # Filter to available columns
    available_cts = []
    for ct in cell_types:
        if ct == "T cells":
            # Combine CD4+ and CD8+ T cells from GT
            if "CD4+ T cells" in gt.columns and "CD8+ T cells" in gt.columns:
                gt["T cells"] = gt["CD4+ T cells"] + gt["CD8+ T cells"]
                if "T cells" in pred.columns:
                    available_cts.append(ct)
        elif ct in gt.columns and ct in pred.columns:
            available_cts.append(ct)

    print(f"  Plotting cell types: {available_cts}")

    n_ct = len(available_cts)
    spot_x = coords["x"].values
    spot_y = coords["y"].values

    # Create figure: 2 rows (GT, Pred) x N columns (cell types)
    fig, axes = plt.subplots(2, n_ct, figsize=(2.5 * n_ct, 5))
    if n_ct == 1:
        axes = axes.reshape(2, 1)

    for i, ct in enumerate(available_cts):
        gt_vals = gt[ct].values
        pred_vals = pred[ct].values

        # Shared colorscale
        vmin = 0
        vmax = max(np.percentile(gt_vals, 98), np.percentile(pred_vals, 98))
        if vmax < 0.01:
            vmax = 0.1

        # Compute correlation for this cell type
        r = np.corrcoef(gt_vals, pred_vals)[0, 1]

        # Top row: Ground Truth
        ax_gt = axes[0, i]
        sc_gt = ax_gt.scatter(spot_x, spot_y, c=gt_vals, s=3, cmap="viridis",
                              vmin=vmin, vmax=vmax, alpha=0.8, rasterized=True)
        ax_gt.set_aspect("equal")
        ax_gt.set_xticks([])
        ax_gt.set_yticks([])

        # Title with cell type name
        short_name = ct.replace(" cells", "").replace("+ ", "+")
        ax_gt.set_title(f"{short_name}\n(r={r:.2f})", fontsize=10, fontweight='bold' if ct == "Fibroblasts" else 'normal')

        if i == 0:
            ax_gt.set_ylabel("Ground Truth", fontsize=10, fontweight='bold')

        # Highlight Fibroblasts with border
        if ct == "Fibroblasts":
            for spine in ax_gt.spines.values():
                spine.set_edgecolor('red')
                spine.set_linewidth(2)
                spine.set_visible(True)

        # Bottom row: Predicted
        ax_pred = axes[1, i]
        sc_pred = ax_pred.scatter(spot_x, spot_y, c=pred_vals, s=3, cmap="viridis",
                                   vmin=vmin, vmax=vmax, alpha=0.8, rasterized=True)
        ax_pred.set_aspect("equal")
        ax_pred.set_xticks([])
        ax_pred.set_yticks([])

        if i == 0:
            ax_pred.set_ylabel("CITEgeist", fontsize=10, fontweight='bold')

        # Highlight Fibroblasts with border
        if ct == "Fibroblasts":
            for spine in ax_pred.spines.values():
                spine.set_edgecolor('red')
                spine.set_linewidth(2)
                spine.set_visible(True)

    # Add colorbar
    cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])
    fig.colorbar(sc_pred, cax=cbar_ax, label="Proportion")

    fig.suptitle(f"Figure 3F: Spatial Comparison — Region {region_id} (αSMA-only Fibroblast GT)",
                 fontsize=12, fontweight='bold', y=1.02)

    plt.tight_layout(rect=[0, 0, 0.91, 1])

    # Save
    out_path = OUTPUT_DIR / f"figure3f_spatial_asma_fix_region{region_id}.png"
    plt.savefig(out_path, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"Saved to: {out_path}")

    plt.savefig(OUTPUT_DIR / f"figure3f_spatial_asma_fix_region{region_id}.pdf",
                bbox_inches='tight', facecolor='white')

    plt.close()

    # Print per-cell-type correlations
    print("\nPer-cell-type correlations:")
    for ct in available_cts:
        r = np.corrcoef(gt[ct].values, pred[ct].values)[0, 1]
        marker = " *** FIXED" if ct == "Fibroblasts" else ""
        print(f"  {ct}: r = {r:.3f}{marker}")


if __name__ == "__main__":
    generate_spatial_panel(region_id=0)
