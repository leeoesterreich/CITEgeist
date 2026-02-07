#!/usr/bin/env python
"""
UMAP Validation of Xenium 10-Cell-Type Ground Truth

This script validates that the 10 granular cell type annotations derived from
RNA k-means clustering correctly capture distinct biological cell populations
by visualizing single-cell Xenium RNA expression via UMAP.

Outputs:
- UMAP colored by cell type
- UMAP colored by k-means cluster
- Marker gene expression panels
- Cell type composition charts
- Quantitative clustering metrics
- Confusion analysis

Author: Alexander Chang
Date: 2026-01-09
"""

import json
import logging
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from scipy import stats
from sklearn.metrics import (
    adjusted_rand_score,
    calinski_harabasz_score,
    davies_bouldin_score,
    silhouette_samples,
    silhouette_score,
)

# Add src to path for local imports
sys.path.insert(0, str(Path(__file__).parent))
from load_xenium import load_xenium_data, split_gex_protein
from rna_cell_types import (
    GRANULAR_CELLTYPE_MAP,
    XENIUM_KIDNEY_RNA_CLUSTER_MAP,
    load_rna_clusters,
)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)

# Constants
DATA_DIR = "/ix1/alee/LO_LAB/General/Public_Data/10x_PublicData/Xenium_RNA_Proteomic_RenalCellCarcinoma"
OUTPUT_DIR = Path(__file__).parent.parent / "figures"

# Comprehensive marker gene lists (will be filtered to available genes)
MARKER_GENES = {
    "T cells": ["CD3D", "CD3E", "CD3G", "CD4", "CD8A", "CD8B", "TRAC", "TRBC1", "TRBC2"],
    "B cells": ["MS4A1", "CD79A", "CD79B", "CD19", "PAX5", "BANK1", "BLK"],
    "Macrophages": ["CD68", "CD163", "CSF1R", "MARCO", "MSR1", "FCGR3A", "C1QA", "C1QB"],
    "Dendritic": ["ITGAX", "HLA-DRA", "HLA-DRB1", "CD1C", "CLEC9A", "FCER1A"],
    "NK cells": ["NCAM1", "NCR1", "KLRD1", "GNLY", "NKG7", "KLRB1", "KLRC1"],
    "Epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19", "CDH1", "MUC1", "TJP1"],
    "Endothelial": ["PECAM1", "VWF", "CDH5", "CD34", "KDR", "CLDN5", "ERG"],
    "Fibroblasts": ["COL1A1", "COL1A2", "DCN", "VIM", "FAP", "THY1", "PDGFRA"],
    "Myofibroblasts": ["ACTA2", "TAGLN", "MYH11", "CNN1", "DES"],
    "Proliferating": ["MKI67", "TOP2A", "PCNA", "MCM2", "CCNB1", "CDK1"],
}

# Cell type color palette (10 colors for 10 cell types)
CELLTYPE_COLORS = {
    "CD8+ T cells": "#E64B35",      # Red
    "Macrophages": "#4DBBD5",        # Cyan
    "Mixed Immune": "#00A087",       # Teal
    "Epithelial": "#3C5488",         # Blue
    "Myofibroblasts": "#F39B7F",     # Salmon
    "Stromal": "#8491B4",            # Gray-blue
    "Endothelial": "#91D1C2",        # Light teal
    "B cells": "#DC0000",            # Dark red
    "Proliferating T": "#7E6148",    # Brown
    "Vascular Stromal": "#B09C85",   # Tan
}


def preprocess_data(adata: sc.AnnData) -> sc.AnnData:
    """
    Standard single-cell RNA preprocessing pipeline.

    Args:
        adata: Raw AnnData object

    Returns:
        Preprocessed AnnData with UMAP coordinates
    """
    logger.info(f"Starting preprocessing with {adata.n_obs} cells and {adata.n_vars} genes")

    # Store raw counts
    adata.layers["counts"] = adata.X.copy()

    # Basic QC filtering
    logger.info("Filtering cells and genes...")
    sc.pp.filter_cells(adata, min_genes=50)
    sc.pp.filter_genes(adata, min_cells=10)
    logger.info(f"After filtering: {adata.n_obs} cells and {adata.n_vars} genes")

    # Normalization
    logger.info("Normalizing...")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # Store normalized counts for marker visualization
    adata.layers["normalized"] = adata.X.copy()

    # HVG selection (use all genes if panel is small)
    n_hvg = min(2000, adata.n_vars - 1)
    if n_hvg < adata.n_vars:
        logger.info(f"Selecting {n_hvg} highly variable genes...")
        sc.pp.highly_variable_genes(adata, n_top_genes=n_hvg, flavor="seurat_v3", layer="counts")
    else:
        logger.info("Using all genes (small panel)")
        adata.var["highly_variable"] = True

    # Scale for PCA
    logger.info("Scaling data...")
    sc.pp.scale(adata, max_value=10)

    # PCA
    n_pcs = min(50, adata.n_vars - 1, adata.n_obs - 1)
    logger.info(f"Computing PCA with {n_pcs} components...")
    sc.tl.pca(adata, n_comps=n_pcs)

    # Neighbors and UMAP
    n_neighbors = min(30, adata.n_obs - 1)
    n_pcs_for_neighbors = min(30, n_pcs)
    logger.info(f"Computing neighbors (n_neighbors={n_neighbors}, n_pcs={n_pcs_for_neighbors})...")
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs_for_neighbors)

    logger.info("Computing UMAP...")
    sc.tl.umap(adata)

    # Also compute Leiden clustering for comparison
    logger.info("Computing Leiden clustering...")
    sc.tl.leiden(adata, resolution=0.8)

    logger.info("Preprocessing complete")
    return adata


def add_cell_type_annotations(adata: sc.AnnData, clusters_df: pd.DataFrame) -> sc.AnnData:
    """
    Add cell type annotations based on k-means clusters.

    Args:
        adata: Preprocessed AnnData
        clusters_df: DataFrame with cluster assignments

    Returns:
        AnnData with cell type annotations
    """
    logger.info("Adding cell type annotations...")

    # Align cell IDs
    common_cells = adata.obs_names.intersection(clusters_df.index)
    logger.info(f"Found {len(common_cells)} cells with cluster assignments")

    # Subset to common cells
    adata = adata[list(common_cells)].copy()
    clusters_aligned = clusters_df.loc[common_cells]

    # Add cluster number
    adata.obs["kmeans_cluster"] = clusters_aligned["Cluster"].astype(str)

    # Add cell type name (granular, unsimplified)
    adata.obs["cell_type"] = clusters_aligned["Cluster"].map(XENIUM_KIDNEY_RNA_CLUSTER_MAP)

    # Log cell type distribution
    logger.info("Cell type distribution:")
    for ct, count in adata.obs["cell_type"].value_counts().items():
        logger.info(f"  {ct}: {count:,} cells ({100*count/adata.n_obs:.1f}%)")

    return adata


def get_available_markers(adata: sc.AnnData) -> Dict[str, List[str]]:
    """
    Filter marker genes to those available in the dataset.

    Args:
        adata: AnnData object

    Returns:
        Dict of cell type to available marker genes
    """
    available_genes = set(adata.var_names)
    available_markers = {}

    for celltype, markers in MARKER_GENES.items():
        found = [m for m in markers if m in available_genes]
        if found:
            available_markers[celltype] = found
            logger.info(f"{celltype}: {len(found)}/{len(markers)} markers available - {found}")
        else:
            logger.warning(f"{celltype}: No markers found in panel")

    return available_markers


def plot_umap_celltype(adata: sc.AnnData, output_path: Path) -> None:
    """Plot UMAP colored by cell type."""
    logger.info("Generating UMAP by cell type...")

    fig, ax = plt.subplots(figsize=(12, 10))

    # Plot each cell type with its color
    for ct in adata.obs["cell_type"].unique():
        mask = adata.obs["cell_type"] == ct
        color = CELLTYPE_COLORS.get(ct, "#999999")
        ax.scatter(
            adata.obsm["X_umap"][mask, 0],
            adata.obsm["X_umap"][mask, 1],
            c=color,
            label=f"{ct} ({mask.sum():,})",
            s=1,
            alpha=0.5,
            rasterized=True,
        )

    ax.set_xlabel("UMAP1", fontsize=12)
    ax.set_ylabel("UMAP2", fontsize=12)
    ax.set_title("Xenium Single-Cell RNA UMAP - Ground Truth Cell Types", fontsize=14)
    ax.legend(
        loc="center left",
        bbox_to_anchor=(1, 0.5),
        markerscale=5,
        frameon=False,
        fontsize=10,
    )

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()
    logger.info(f"Saved: {output_path}")


def plot_umap_clusters(adata: sc.AnnData, output_path: Path) -> None:
    """Plot UMAP colored by k-means cluster number."""
    logger.info("Generating UMAP by k-means cluster...")

    fig, ax = plt.subplots(figsize=(10, 8))

    scatter = ax.scatter(
        adata.obsm["X_umap"][:, 0],
        adata.obsm["X_umap"][:, 1],
        c=adata.obs["kmeans_cluster"].astype(int),
        cmap="tab10",
        s=1,
        alpha=0.5,
        rasterized=True,
    )

    ax.set_xlabel("UMAP1", fontsize=12)
    ax.set_ylabel("UMAP2", fontsize=12)
    ax.set_title("Xenium UMAP - K-means Clusters (k=10)", fontsize=14)
    plt.colorbar(scatter, ax=ax, label="Cluster")

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()
    logger.info(f"Saved: {output_path}")


def plot_marker_genes(adata: sc.AnnData, available_markers: Dict[str, List[str]], output_path: Path) -> None:
    """Plot comprehensive marker gene expression panels."""
    logger.info("Generating marker gene expression panels...")

    # Flatten markers and keep unique
    all_markers = []
    for markers in available_markers.values():
        all_markers.extend(markers)
    all_markers = list(dict.fromkeys(all_markers))  # Remove duplicates, preserve order

    if not all_markers:
        logger.warning("No markers available for plotting")
        return

    # Calculate grid size
    n_markers = len(all_markers)
    n_cols = 4
    n_rows = (n_markers + n_cols - 1) // n_cols

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(4 * n_cols, 3.5 * n_rows))
    axes = axes.flatten() if n_markers > 1 else [axes]

    # Use normalized layer for visualization
    adata_plot = adata.copy()
    adata_plot.X = adata_plot.layers["normalized"]

    for idx, marker in enumerate(all_markers):
        ax = axes[idx]

        if marker in adata_plot.var_names:
            expr = adata_plot[:, marker].X
            if hasattr(expr, "toarray"):
                expr = expr.toarray().flatten()
            else:
                expr = np.array(expr).flatten()

            scatter = ax.scatter(
                adata_plot.obsm["X_umap"][:, 0],
                adata_plot.obsm["X_umap"][:, 1],
                c=expr,
                cmap="viridis",
                s=0.5,
                alpha=0.5,
                rasterized=True,
            )
            plt.colorbar(scatter, ax=ax, shrink=0.8)

        ax.set_title(marker, fontsize=10)
        ax.set_xticks([])
        ax.set_yticks([])

    # Hide unused axes
    for idx in range(n_markers, len(axes)):
        axes[idx].set_visible(False)

    plt.suptitle("Marker Gene Expression on UMAP", fontsize=14, y=1.02)
    plt.tight_layout()
    plt.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close()
    logger.info(f"Saved: {output_path}")


def plot_composition(adata: sc.AnnData, output_path: Path) -> None:
    """Plot cell type composition charts."""
    logger.info("Generating composition charts...")

    counts = adata.obs["cell_type"].value_counts()
    colors = [CELLTYPE_COLORS.get(ct, "#999999") for ct in counts.index]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    # Bar chart
    bars = ax1.barh(counts.index, counts.values, color=colors)
    ax1.set_xlabel("Number of Cells", fontsize=12)
    ax1.set_title("Cell Type Counts", fontsize=14)

    # Add count labels
    for bar, count in zip(bars, counts.values):
        ax1.text(
            count + counts.max() * 0.01,
            bar.get_y() + bar.get_height() / 2,
            f"{count:,}",
            va="center",
            fontsize=9,
        )

    # Pie chart
    ax2.pie(
        counts.values,
        labels=[f"{ct}\n({100*v/counts.sum():.1f}%)" for ct, v in counts.items()],
        colors=colors,
        autopct="",
        startangle=90,
    )
    ax2.set_title("Cell Type Proportions", fontsize=14)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()
    logger.info(f"Saved: {output_path}")


def plot_silhouette_per_celltype(adata: sc.AnnData, silhouette_per_cell: np.ndarray, output_path: Path) -> None:
    """Plot silhouette scores per cell type."""
    logger.info("Generating silhouette plot...")

    df = pd.DataFrame({
        "cell_type": adata.obs["cell_type"].values,
        "silhouette": silhouette_per_cell,
    })

    # Order by mean silhouette
    order = df.groupby("cell_type")["silhouette"].mean().sort_values(ascending=False).index

    fig, ax = plt.subplots(figsize=(10, 6))

    sns.boxplot(
        data=df,
        x="cell_type",
        y="silhouette",
        order=order,
        palette=[CELLTYPE_COLORS.get(ct, "#999999") for ct in order],
        ax=ax,
    )

    ax.axhline(y=0, color="red", linestyle="--", alpha=0.5)
    ax.set_xlabel("Cell Type", fontsize=12)
    ax.set_ylabel("Silhouette Score", fontsize=12)
    ax.set_title("Silhouette Score by Cell Type", fontsize=14)
    plt.xticks(rotation=45, ha="right")

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()
    logger.info(f"Saved: {output_path}")


def plot_confusion_matrix(adata: sc.AnnData, output_path: Path) -> None:
    """Plot confusion matrix based on nearest neighbor analysis."""
    logger.info("Generating confusion matrix...")

    # Get connectivity matrix from neighbors
    if "connectivities" not in adata.obsp:
        logger.warning("No connectivity matrix found, skipping confusion matrix")
        return

    conn = adata.obsp["connectivities"]
    cell_types = adata.obs["cell_type"].values
    unique_types = sorted(adata.obs["cell_type"].unique())
    n_types = len(unique_types)

    # Build confusion matrix: for each cell, what % of neighbors are each type
    confusion = np.zeros((n_types, n_types))
    type_to_idx = {t: i for i, t in enumerate(unique_types)}

    # Sample cells for efficiency (full matrix too slow)
    n_sample = min(10000, adata.n_obs)
    sample_idx = np.random.choice(adata.n_obs, n_sample, replace=False)

    for i in sample_idx:
        my_type = type_to_idx[cell_types[i]]
        neighbors = conn[i].nonzero()[1]
        if len(neighbors) > 0:
            for j in neighbors:
                neighbor_type = type_to_idx[cell_types[j]]
                confusion[my_type, neighbor_type] += 1

    # Normalize rows to percentages
    row_sums = confusion.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1  # Avoid division by zero
    confusion = confusion / row_sums * 100

    fig, ax = plt.subplots(figsize=(10, 8))

    sns.heatmap(
        confusion,
        annot=True,
        fmt=".1f",
        cmap="Blues",
        xticklabels=unique_types,
        yticklabels=unique_types,
        ax=ax,
    )

    ax.set_xlabel("Neighbor Cell Type", fontsize=12)
    ax.set_ylabel("Cell Type", fontsize=12)
    ax.set_title("Nearest Neighbor Confusion Matrix (%)", fontsize=14)
    plt.xticks(rotation=45, ha="right")
    plt.yticks(rotation=0)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()
    logger.info(f"Saved: {output_path}")


def plot_marker_heatmap(adata: sc.AnnData, available_markers: Dict[str, List[str]], output_path: Path) -> None:
    """Plot mean marker expression heatmap per cell type."""
    logger.info("Generating marker expression heatmap...")

    # Get normalized expression
    X = adata.layers["normalized"]
    if hasattr(X, "toarray"):
        X = X.toarray()

    # Flatten markers
    all_markers = []
    for markers in available_markers.values():
        all_markers.extend(markers)
    all_markers = list(dict.fromkeys(all_markers))

    # Filter to available
    available = [m for m in all_markers if m in adata.var_names]
    if not available:
        logger.warning("No markers available for heatmap")
        return

    # Calculate mean expression per cell type
    cell_types = sorted(adata.obs["cell_type"].unique())
    expr_matrix = np.zeros((len(cell_types), len(available)))

    for i, ct in enumerate(cell_types):
        mask = adata.obs["cell_type"] == ct
        for j, marker in enumerate(available):
            gene_idx = list(adata.var_names).index(marker)
            expr_matrix[i, j] = X[mask, gene_idx].mean()

    # Z-score normalize columns
    expr_z = (expr_matrix - expr_matrix.mean(axis=0)) / (expr_matrix.std(axis=0) + 1e-6)

    fig, ax = plt.subplots(figsize=(max(12, len(available) * 0.5), 8))

    sns.heatmap(
        expr_z,
        xticklabels=available,
        yticklabels=cell_types,
        cmap="RdBu_r",
        center=0,
        ax=ax,
    )

    ax.set_xlabel("Marker Gene", fontsize=12)
    ax.set_ylabel("Cell Type", fontsize=12)
    ax.set_title("Marker Gene Expression (Z-score)", fontsize=14)
    plt.xticks(rotation=90, fontsize=8)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()
    logger.info(f"Saved: {output_path}")


def calculate_metrics(adata: sc.AnnData) -> Tuple[Dict, np.ndarray]:
    """
    Calculate comprehensive clustering metrics.

    Returns:
        Tuple of (metrics dict, per-cell silhouette scores)
    """
    logger.info("Calculating clustering metrics...")

    # Get PCA coordinates for metric calculation
    X_pca = adata.obsm["X_pca"]
    labels = adata.obs["cell_type"].astype("category").cat.codes.values
    leiden_labels = adata.obs["leiden"].astype("category").cat.codes.values

    metrics = {}

    # Global metrics
    logger.info("Computing silhouette score...")
    silhouette_per_cell = silhouette_samples(X_pca, labels)
    metrics["silhouette_score"] = float(np.mean(silhouette_per_cell))

    logger.info("Computing Calinski-Harabasz index...")
    metrics["calinski_harabasz"] = float(calinski_harabasz_score(X_pca, labels))

    logger.info("Computing Davies-Bouldin index...")
    metrics["davies_bouldin"] = float(davies_bouldin_score(X_pca, labels))

    logger.info("Computing Adjusted Rand Index vs Leiden...")
    metrics["ari_vs_leiden"] = float(adjusted_rand_score(labels, leiden_labels))

    # Per-cell-type silhouette
    metrics["silhouette_per_celltype"] = {}
    for ct in adata.obs["cell_type"].unique():
        mask = adata.obs["cell_type"] == ct
        metrics["silhouette_per_celltype"][ct] = float(silhouette_per_cell[mask].mean())

    # Cluster centroids and inter-cluster distances
    cell_types = sorted(adata.obs["cell_type"].unique())
    centroids = np.zeros((len(cell_types), X_pca.shape[1]))
    for i, ct in enumerate(cell_types):
        mask = adata.obs["cell_type"] == ct
        centroids[i] = X_pca[mask].mean(axis=0)

    # Compute pairwise distances
    from scipy.spatial.distance import cdist
    distances = cdist(centroids, centroids, metric="euclidean")

    metrics["inter_cluster_distances"] = {}
    for i, ct1 in enumerate(cell_types):
        for j, ct2 in enumerate(cell_types):
            if i < j:
                key = f"{ct1}_vs_{ct2}"
                metrics["inter_cluster_distances"][key] = float(distances[i, j])

    # Nearest neighbor purity
    logger.info("Computing nearest neighbor purity...")
    if "connectivities" in adata.obsp:
        conn = adata.obsp["connectivities"]
        cell_types_arr = adata.obs["cell_type"].values

        purity_scores = []
        n_sample = min(10000, adata.n_obs)
        sample_idx = np.random.choice(adata.n_obs, n_sample, replace=False)

        for i in sample_idx:
            neighbors = conn[i].nonzero()[1]
            if len(neighbors) > 0:
                my_type = cell_types_arr[i]
                same_type = sum(1 for j in neighbors if cell_types_arr[j] == my_type)
                purity_scores.append(same_type / len(neighbors))

        metrics["neighbor_purity"] = float(np.mean(purity_scores))

        # Per-cell-type purity
        metrics["neighbor_purity_per_celltype"] = {}
        for ct in adata.obs["cell_type"].unique():
            ct_indices = np.where(adata.obs["cell_type"] == ct)[0]
            ct_sample = np.random.choice(ct_indices, min(1000, len(ct_indices)), replace=False)

            ct_purity = []
            for i in ct_sample:
                neighbors = conn[i].nonzero()[1]
                if len(neighbors) > 0:
                    same_type = sum(1 for j in neighbors if cell_types_arr[j] == ct)
                    ct_purity.append(same_type / len(neighbors))

            if ct_purity:
                metrics["neighbor_purity_per_celltype"][ct] = float(np.mean(ct_purity))

    # Cell counts
    metrics["cell_counts"] = {
        ct: int(count) for ct, count in adata.obs["cell_type"].value_counts().items()
    }
    metrics["total_cells"] = int(adata.n_obs)

    logger.info("Metrics calculation complete")
    return metrics, silhouette_per_cell


def run_statistical_tests(adata: sc.AnnData, available_markers: Dict[str, List[str]]) -> Dict:
    """Run ANOVA and effect size calculations for marker genes."""
    logger.info("Running statistical tests...")

    results = {}

    # Get normalized expression
    X = adata.layers["normalized"]
    if hasattr(X, "toarray"):
        X = X.toarray()

    # Flatten markers
    all_markers = []
    for markers in available_markers.values():
        all_markers.extend(markers)
    all_markers = list(dict.fromkeys(all_markers))

    cell_types = sorted(adata.obs["cell_type"].unique())

    for marker in all_markers:
        if marker not in adata.var_names:
            continue

        gene_idx = list(adata.var_names).index(marker)
        expr = X[:, gene_idx]

        # Group by cell type
        groups = [expr[adata.obs["cell_type"] == ct] for ct in cell_types]

        # ANOVA
        try:
            f_stat, p_value = stats.f_oneway(*groups)
            results[marker] = {
                "f_statistic": float(f_stat) if not np.isnan(f_stat) else None,
                "p_value": float(p_value) if not np.isnan(p_value) else None,
            }
        except Exception as e:
            logger.warning(f"ANOVA failed for {marker}: {e}")
            results[marker] = {"f_statistic": None, "p_value": None}

        # Effect size (eta-squared)
        if results[marker]["p_value"] is not None:
            grand_mean = expr.mean()
            ss_between = sum(len(g) * (g.mean() - grand_mean) ** 2 for g in groups)
            ss_total = sum((expr - grand_mean) ** 2)
            if ss_total > 0:
                results[marker]["eta_squared"] = float(ss_between / ss_total)
            else:
                results[marker]["eta_squared"] = None

    logger.info(f"Statistical tests complete for {len(results)} markers")
    return results


def main():
    """Main validation pipeline."""
    logger.info("=" * 70)
    logger.info("UMAP VALIDATION OF XENIUM GROUND TRUTH")
    logger.info("=" * 70)

    # Create output directory
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    logger.info(f"Output directory: {OUTPUT_DIR}")

    # Step 1: Load Xenium data
    logger.info("\n--- Step 1: Loading Xenium data ---")
    adata = load_xenium_data(DATA_DIR)
    adata_gex, _ = split_gex_protein(adata)
    logger.info(f"Loaded {adata_gex.n_obs} cells with {adata_gex.n_vars} genes")

    # Step 2: Load k-means clusters
    logger.info("\n--- Step 2: Loading k-means cluster assignments ---")
    clusters_df = load_rna_clusters(DATA_DIR, n_clusters=10)

    # Step 3: Preprocess data
    logger.info("\n--- Step 3: Preprocessing data ---")
    adata_gex = preprocess_data(adata_gex)

    # Step 4: Add cell type annotations
    logger.info("\n--- Step 4: Adding cell type annotations ---")
    adata_gex = add_cell_type_annotations(adata_gex, clusters_df)

    # Step 5: Get available markers
    logger.info("\n--- Step 5: Identifying available marker genes ---")
    available_markers = get_available_markers(adata_gex)

    # Step 6: Calculate metrics
    logger.info("\n--- Step 6: Calculating clustering metrics ---")
    metrics, silhouette_per_cell = calculate_metrics(adata_gex)

    # Step 7: Run statistical tests
    logger.info("\n--- Step 7: Running statistical tests ---")
    stat_results = run_statistical_tests(adata_gex, available_markers)
    metrics["marker_statistics"] = stat_results

    # Step 8: Generate figures
    logger.info("\n--- Step 8: Generating figures ---")
    plot_umap_celltype(adata_gex, OUTPUT_DIR / "umap_celltype_validation.png")
    plot_umap_clusters(adata_gex, OUTPUT_DIR / "umap_kmeans_clusters.png")
    plot_marker_genes(adata_gex, available_markers, OUTPUT_DIR / "umap_marker_genes_comprehensive.png")
    plot_composition(adata_gex, OUTPUT_DIR / "celltype_composition.png")
    plot_silhouette_per_celltype(adata_gex, silhouette_per_cell, OUTPUT_DIR / "silhouette_per_celltype.png")
    plot_confusion_matrix(adata_gex, OUTPUT_DIR / "confusion_matrix.png")
    plot_marker_heatmap(adata_gex, available_markers, OUTPUT_DIR / "marker_expression_heatmap.png")

    # Step 9: Save metrics
    logger.info("\n--- Step 9: Saving metrics ---")
    metrics_path = OUTPUT_DIR / "validation_metrics.json"
    with open(metrics_path, "w") as f:
        json.dump(metrics, f, indent=2)
    logger.info(f"Saved metrics to: {metrics_path}")

    # Print summary
    logger.info("\n" + "=" * 70)
    logger.info("VALIDATION SUMMARY")
    logger.info("=" * 70)
    logger.info(f"Total cells analyzed: {metrics['total_cells']:,}")
    logger.info(f"Silhouette score: {metrics['silhouette_score']:.3f}")
    logger.info(f"Calinski-Harabasz index: {metrics['calinski_harabasz']:.1f}")
    logger.info(f"Davies-Bouldin index: {metrics['davies_bouldin']:.3f}")
    logger.info(f"ARI vs Leiden: {metrics['ari_vs_leiden']:.3f}")
    if "neighbor_purity" in metrics:
        logger.info(f"Neighbor purity: {metrics['neighbor_purity']:.3f}")

    logger.info("\nSilhouette per cell type:")
    for ct, score in sorted(metrics["silhouette_per_celltype"].items(), key=lambda x: -x[1]):
        logger.info(f"  {ct}: {score:.3f}")

    logger.info("\n" + "=" * 70)
    logger.info("VALIDATION COMPLETE")
    logger.info("=" * 70)
    logger.info(f"Figures saved to: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
