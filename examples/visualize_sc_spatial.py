"""Visualize single-cell spatial outputs from CITEgeist 12-patient pipeline.

Overlays CITEgeist single-cell deconvolution results on H&E tissue images.

Generates:
1. Cell type spatial scatter on H&E (4 representative samples)
2. Marker gene expression on H&E (CD3E, CD68, EPCAM, VIM, MDK)
3. Summary stacked bar across all 12 samples
"""
import os
import sys
import json
import scanpy as sc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd
from PIL import Image

_version = os.environ.get("MORPHOLOGY_VERSION", "v3")
OUTPUT_BASE = f"/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/output/morphology_assignment_{_version}"
RAW_BASE = "/ix1/alee/LO_LAB/General/Lab_Data/20250210_CITEGeistPublicData_GEO_Alex/processed_files"
FIG_DIR = "/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/output/sc_spatial_figures"
os.makedirs(FIG_DIR, exist_ok=True)

# 4 representative samples: 2 biopsy, 2 surgical
SAMPLES = [
    "HCC22-088-P1-S2",        # biopsy
    "HCC22-088-P3-S1_A",      # biopsy
    "HCC22-088-P4-S2_1i_rep", # surgical
    "HCC22-088-P6-S2_D",      # surgical
]

MARKER_GENES = ["CD3E", "CD68", "EPCAM", "VIM", "MDK"]

TYPE_COLORS = {
    "Cancer_Basal": "#E41A1C",
    "Cancer_Luminal": "#FF7F00",
    "CD4_T_Cells": "#377EB8",
    "CD8_T_Cells": "#4DAF4A",
    "Macrophages": "#984EA3",
    "Dendritic_Cells": "#A65628",
    "B_Cells": "#F781BF",
    "Fibroblasts": "#999999",
    "Endothelial": "#66C2A5",
    "Monocytes": "#E7298A",
}


def load_he_image(sample_name):
    """Load H&E hires image and scale factor for a sample."""
    spatial_dir = os.path.join(RAW_BASE, sample_name, "outs", "spatial")
    img_path = os.path.join(spatial_dir, "tissue_hires_image.png")
    sf_path = os.path.join(spatial_dir, "scalefactors_json.json")

    if not os.path.exists(img_path):
        print(f"  WARNING: No H&E image for {sample_name}")
        return None, None

    img = np.array(Image.open(img_path))
    with open(sf_path) as f:
        sf = json.load(f)

    return img, sf["tissue_hires_scalef"]


def load_sample(sample_name):
    """Load single-cell h5ad and H&E image for a sample."""
    sc_path = os.path.join(OUTPUT_BASE, sample_name, f"{sample_name}_single_cell.h5ad")
    if not os.path.exists(sc_path):
        print(f"  WARNING: {sc_path} not found, skipping")
        return None, None, None

    adata = sc.read_h5ad(sc_path)
    img, scale = load_he_image(sample_name)

    print(f"  {sample_name}: {adata.n_obs} cells, {adata.n_vars} genes, "
          f"{adata.obs['cell_type'].nunique()} types"
          f"{', H&E loaded' if img is not None else ', NO H&E'}")

    return adata, img, scale


def plot_celltype_on_he(adatas, images, scales, sample_names):
    """Plot cell type assignments overlaid on H&E tissue images."""
    n = len(adatas)
    fig, axes = plt.subplots(1, n, figsize=(7 * n, 7), dpi=150)
    if n == 1:
        axes = [axes]

    for ax, adata, img, scale, name in zip(axes, adatas, images, scales, sample_names):
        # Show H&E background
        if img is not None:
            ax.imshow(img, alpha=0.8)

        types = adata.obs['cell_type']
        coords = adata.obsm['spatial']

        # Scale fullres coords to hires image space
        if scale is not None:
            coords_scaled = coords * scale
        else:
            coords_scaled = coords

        for ct in TYPE_COLORS:
            mask = (types == ct).values
            if mask.sum() == 0:
                continue
            ax.scatter(
                coords_scaled[mask, 0], coords_scaled[mask, 1],
                c=TYPE_COLORS[ct], label=ct,
                s=4, alpha=0.85, edgecolors='none', rasterized=True
            )

        short_name = name.replace("HCC22-088-", "")
        ax.set_title(f"{short_name}\n({adata.n_obs} cells)",
                     fontsize=12, fontweight='bold')
        ax.set_aspect('equal')
        ax.set_xticks([])
        ax.set_yticks([])

    # Shared legend
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='lower center', ncol=5,
               fontsize=10, markerscale=3, frameon=True,
               facecolor='white', edgecolor='gray',
               bbox_to_anchor=(0.5, -0.03))

    fig.suptitle("CITEgeist Single-Cell Deconvolution on H&E",
                 fontsize=15, fontweight='bold', y=1.02)
    plt.tight_layout()
    out = os.path.join(FIG_DIR, "celltype_on_he_4samples.png")
    fig.savefig(out, dpi=150, bbox_inches='tight', facecolor='white')
    print(f"  Saved: {out}")
    plt.close()


def plot_marker_genes_on_he(adatas, images, scales, sample_names):
    """Plot marker gene expression overlaid on H&E tissue images."""
    n_samples = len(adatas)
    n_genes = len(MARKER_GENES)

    fig, axes = plt.subplots(n_genes, n_samples,
                             figsize=(6 * n_samples, 5 * n_genes), dpi=150)

    for j, (adata, img, scale, name) in enumerate(
            zip(adatas, images, scales, sample_names)):
        # Log1p normalize for visualization
        adata_vis = adata.copy()
        sc.pp.normalize_total(adata_vis, target_sum=1e4)
        sc.pp.log1p(adata_vis)

        coords = adata_vis.obsm['spatial']
        coords_scaled = coords * scale if scale is not None else coords

        for i, gene in enumerate(MARKER_GENES):
            ax = axes[i, j] if n_samples > 1 else axes[i]

            # H&E background
            if img is not None:
                ax.imshow(img, alpha=0.6)

            if gene in adata_vis.var_names:
                expr = np.array(adata_vis[:, gene].X).flatten()
                vmax = (np.percentile(expr[expr > 0], 99)
                        if (expr > 0).sum() > 10 else max(expr.max(), 0.01))

                sc_plot = ax.scatter(
                    coords_scaled[:, 0], coords_scaled[:, 1],
                    c=expr, cmap='magma', s=3, alpha=0.85,
                    vmin=0, vmax=vmax, edgecolors='none', rasterized=True
                )
                plt.colorbar(sc_plot, ax=ax, shrink=0.5, pad=0.02)
            else:
                ax.text(0.5, 0.5, f"{gene}\nnot in\ngene set",
                        transform=ax.transAxes, ha='center', va='center',
                        fontsize=11, color='red', fontweight='bold')

            ax.set_aspect('equal')
            ax.set_xticks([])
            ax.set_yticks([])

            if j == 0:
                ax.set_ylabel(gene, fontsize=13, fontweight='bold')
            if i == 0:
                ax.set_title(name.replace("HCC22-088-", ""),
                             fontsize=12, fontweight='bold')

    fig.suptitle("Single-Cell Marker Gene Expression on H&E (log1p)",
                 fontsize=15, fontweight='bold', y=1.01)
    plt.tight_layout()
    out = os.path.join(FIG_DIR, "marker_genes_on_he_4samples.png")
    fig.savefig(out, dpi=150, bbox_inches='tight', facecolor='white')
    print(f"  Saved: {out}")
    plt.close()


def plot_summary_stats():
    """Plot cell count and type distribution summary across all 12 samples."""
    all_samples = sorted([d for d in os.listdir(OUTPUT_BASE)
                          if os.path.isdir(os.path.join(OUTPUT_BASE, d))
                          and d.startswith("HCC")])

    records = []
    for s in all_samples:
        path = os.path.join(OUTPUT_BASE, s, f"{s}_single_cell.h5ad")
        if not os.path.exists(path):
            continue
        ad = sc.read_h5ad(path)
        for ct, count in ad.obs['cell_type'].value_counts().items():
            records.append({
                "sample": s.replace("HCC22-088-", ""),
                "cell_type": ct, "count": count
            })

    df = pd.DataFrame(records)
    pivot = df.pivot_table(index='sample', columns='cell_type',
                           values='count', fill_value=0)

    fig, ax = plt.subplots(figsize=(14, 6), dpi=150)
    ordered_types = [t for t in TYPE_COLORS if t in pivot.columns]
    colors = [TYPE_COLORS[t] for t in ordered_types]
    pivot[ordered_types].plot(kind='bar', stacked=True, color=colors,
                             ax=ax, width=0.8)

    ax.set_ylabel("Number of Cells", fontsize=12)
    ax.set_xlabel("")
    ax.set_title("Single-Cell Counts by Type Across 12 Patient Samples",
                 fontsize=13, fontweight='bold')
    ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=9)
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()

    out = os.path.join(FIG_DIR, "cellcount_summary_12samples.png")
    fig.savefig(out, dpi=150, bbox_inches='tight', facecolor='white')
    print(f"  Saved: {out}")
    plt.close()

    print("\n=== SINGLE-CELL SUMMARY (ALL 12 SAMPLES) ===")
    totals = pivot[ordered_types].sum()
    print(f"Total cells: {int(totals.sum())} across {len(pivot)} samples")
    for t in ordered_types:
        print(f"  {t:20s}: {int(totals[t]):>6d}")


def main():
    print("Loading samples with H&E images...")
    adatas, images, scales, valid_names = [], [], [], []
    for s in SAMPLES:
        ad, img, scale = load_sample(s)
        if ad is not None:
            adatas.append(ad)
            images.append(img)
            scales.append(scale)
            valid_names.append(s)

    if not adatas:
        print("ERROR: No samples loaded!")
        sys.exit(1)

    print(f"\nLoaded {len(adatas)} samples, generating plots...\n")

    print("--- Cell type on H&E ---")
    plot_celltype_on_he(adatas, images, scales, valid_names)

    print("\n--- Marker genes on H&E ---")
    plot_marker_genes_on_he(adatas, images, scales, valid_names)

    print("\n--- Summary stats (all 12 samples) ---")
    plot_summary_stats()

    print(f"\nAll figures saved to: {FIG_DIR}")


if __name__ == "__main__":
    main()
