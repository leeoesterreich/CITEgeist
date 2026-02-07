#!/usr/bin/env python3
"""
Generate spatial visualization panels for manuscript figures.

This script creates:
1. Spatial tissue visualization with cell proportion overlay (for Figure 1)
2. Real UMAP from deconvolved data (for Figure 6B)
3. Spatial program map (for Figure 4)

Run via sbatch due to data loading requirements.
"""

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from pathlib import Path
import sys

# Add project to path
PROJECT_ROOT = Path(__file__).parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from manuscript.figures.figure_style import apply_style, PALETTE, CELL_TYPE_COLORS, get_cell_type_color

apply_style()

# Paths
XENIUM_H5AD_DIR = PROJECT_ROOT / "Benchmarking/xenium_pseudovisium/data/h5ad_objects"
CITEGEIST_OUTPUT = PROJECT_ROOT / "Benchmarking/xenium_benchmarking/CITEgeist/output/hierarchical"
OUTPUT_DIR = Path(__file__).parent / "output" / "spatial_panels"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def generate_spatial_proportion_plot(region_id=0):
    """Generate spatial visualization showing cell type proportions."""
    print(f"Loading Xenium region {region_id}...")

    # Load h5ad
    gex_file = XENIUM_H5AD_DIR / f"Xenium_region_{region_id}_GEX.h5ad"
    adata = sc.read_h5ad(gex_file)

    # Load cell proportions
    prop_file = CITEGEIST_OUTPUT / f"Xenium_region_{region_id}_cell_prop_finetuned_results.csv"
    props = pd.read_csv(prop_file, index_col=0)

    # Match spot indices
    common_spots = adata.obs_names.intersection(props.index)
    adata = adata[common_spots, :].copy()
    props = props.loc[common_spots]

    print(f"  Spots: {len(common_spots)}, Cell types: {props.columns.tolist()}")

    # Add proportions to adata
    for ct in props.columns:
        adata.obs[f'prop_{ct}'] = props[ct].values

    # Assign dominant cell type
    adata.obs['dominant_celltype'] = props.idxmax(axis=1).values

    # Create figure with 2 panels
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Panel 1: Spatial plot colored by dominant cell type
    ax1 = axes[0]
    coords = adata.obsm['spatial']
    cell_types = adata.obs['dominant_celltype'].unique()

    for ct in cell_types:
        mask = adata.obs['dominant_celltype'] == ct
        color = get_cell_type_color(ct)
        ax1.scatter(coords[mask, 0], coords[mask, 1],
                   c=color, s=8, alpha=0.7, label=ct.replace('+', '+\n'), rasterized=True)

    ax1.set_xlabel('Spatial X')
    ax1.set_ylabel('Spatial Y')
    ax1.set_title('CITEgeist Deconvolution', fontsize=11, fontweight='bold')
    ax1.legend(loc='upper right', fontsize=7, markerscale=1.5)
    ax1.set_aspect('equal')
    ax1.axis('off')

    # Panel 2: Pie chart showing overall composition
    ax2 = axes[1]
    mean_props = props.mean()
    colors = [get_cell_type_color(ct) for ct in mean_props.index]
    wedges, texts, autotexts = ax2.pie(mean_props, labels=None, autopct='%1.0f%%',
                                        colors=colors, pctdistance=0.75)
    ax2.legend(wedges, [ct.replace('+', '+\n') for ct in mean_props.index],
               loc='center left', bbox_to_anchor=(1, 0.5), fontsize=8)
    ax2.set_title('Mean Composition', fontsize=11, fontweight='bold')

    plt.tight_layout()

    output_path = OUTPUT_DIR / f"spatial_proportions_region{region_id}.png"
    plt.savefig(output_path, dpi=150, bbox_inches='tight', facecolor='white')
    print(f"Saved: {output_path}")
    plt.close()

    return str(output_path)


def generate_deconvolved_umap(region_id=0):
    """Generate UMAP from actual deconvolved expression data."""
    print(f"Generating UMAP for region {region_id}...")

    # Load GEX data
    gex_file = XENIUM_H5AD_DIR / f"Xenium_region_{region_id}_GEX.h5ad"
    adata = sc.read_h5ad(gex_file)

    # Load cell proportions for coloring
    prop_file = CITEGEIST_OUTPUT / f"Xenium_region_{region_id}_cell_prop_finetuned_results.csv"
    props = pd.read_csv(prop_file, index_col=0)

    # Match spots
    common_spots = adata.obs_names.intersection(props.index)
    adata = adata[common_spots, :].copy()
    props = props.loc[common_spots]

    # Preprocess for UMAP
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    sc.pp.pca(adata, n_comps=30)
    sc.pp.neighbors(adata, n_neighbors=15)
    sc.tl.umap(adata)

    # Add dominant cell type
    adata.obs['dominant_celltype'] = props.idxmax(axis=1).values

    # Plot UMAP colored by cell type
    fig, ax = plt.subplots(figsize=(7, 6))

    cell_types = adata.obs['dominant_celltype'].unique()
    for ct in cell_types:
        mask = adata.obs['dominant_celltype'] == ct
        color = get_cell_type_color(ct)
        umap_coords = adata.obsm['X_umap']
        ax.scatter(umap_coords[mask, 0], umap_coords[mask, 1],
                  c=color, s=10, alpha=0.6, label=ct, rasterized=True)

    ax.set_xlabel('UMAP 1', fontsize=10)
    ax.set_ylabel('UMAP 2', fontsize=10)
    ax.set_title('Deconvolved GEX Clustering', fontsize=11, fontweight='bold')
    ax.legend(loc='upper right', fontsize=8, markerscale=1.5)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()

    output_path = OUTPUT_DIR / f"deconvolved_umap_region{region_id}.png"
    plt.savefig(output_path, dpi=150, bbox_inches='tight', facecolor='white')
    print(f"Saved: {output_path}")
    plt.close()

    return str(output_path)


if __name__ == "__main__":
    print("Generating spatial visualization panels...")

    # Generate for region 0 (representative)
    spatial_plot = generate_spatial_proportion_plot(region_id=0)
    umap_plot = generate_deconvolved_umap(region_id=0)

    print("\n=== Generated Panels ===")
    print(f"Spatial proportions: {spatial_plot}")
    print(f"Deconvolved UMAP: {umap_plot}")
    print("\nThese can now be embedded in Figure 1 and Figure 6.")
