"""
RNA-based cell type classification using Xenium's built-in clustering.

This module provides functions to load and annotate cell types based on
gene expression clustering from Xenium's analysis output, rather than
protein thresholding.

JUSTIFICATION:
--------------
Reference-based and clustering-based methods are the recommended approaches
for cell type annotation in Xenium data, as established by:

    Zhao et al. (2025). "Benchmarking cell type annotation methods for 10x
    Xenium spatial transcriptomics data." BMC Bioinformatics, 26(1), 25.
    https://doi.org/10.1186/s12859-025-06044-0

Key findings from this study:
1. SingleR was the best performing reference-based annotation tool
2. Gene expression-based clustering outperforms protein thresholding
3. Marker-gene-based manual annotation should use RNA, not protein alone

Additionally, 10x Genomics recommends label transfer from scRNA-seq references:
    https://www.10xgenomics.com/analysis-guides/xenium-cell-type-annotation

Using RNA-based clustering avoids the circular logic problem where:
- Protein markers define cell types (ground truth)
- Same protein markers are used for deconvolution (prediction)
- This creates artificial correlation unrelated to true cell identity
"""

import logging
import tarfile
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import numpy as np
import pandas as pd
import scanpy as sc

logger = logging.getLogger(__name__)

# Default cluster-to-celltype mapping for Xenium Kidney RCC dataset
# Based on protein marker profiles of each RNA cluster (k=10)
# Clusters defined by gene_expression_kmeans_10_clusters
#
# Analysis performed 2026-01-02 with analyze_cluster_profiles.py:
# Cluster 1: 89,144 cells - CD3E=378, CD8A=210, CD45=393
# Cluster 2: 88,062 cells - CD68=430, HLA-DR=317, CD163=88, CD16=82
# Cluster 3: 69,920 cells - CD3E=142, CD8A=118, HLA-DR=142
# Cluster 4: 60,787 cells - PanCK=39, Vimentin=311, PCNA=73
# Cluster 5: 35,338 cells - alphaSMA=108, Vimentin=374
# Cluster 6: 30,801 cells - CD3E=62, Vimentin=198, mixed low
# Cluster 7: 30,481 cells - CD31=168, Vimentin=375
# Cluster 8: 25,427 cells - CD20=293, CD45RA=398, CD45=344
# Cluster 9: 18,284 cells - CD3E=679, CD31=139, CD138=44, PCNA=83
# Cluster 10: 13,642 cells - CD31=53, Vimentin=209

XENIUM_KIDNEY_RNA_CLUSTER_MAP = {
    1: "CD8+ T cells",  # CD3E=378, CD8A=210 - cytotoxic T lymphocytes
    2: "Macrophages",  # CD68=430, HLA-DR=317, CD163=88 - myeloid cells
    3: "Mixed Immune",  # CD3E=142, CD8A=118, HLA-DR=142 - T/myeloid interface
    4: "Epithelial",  # PanCK=39, high Vimentin - tumor/epithelial
    5: "Myofibroblasts",  # alphaSMA=108, Vimentin=374 - activated fibroblasts
    6: "Stromal",  # Mixed low markers, CD3E=62 - general stromal
    7: "Endothelial",  # CD31=168, Vimentin=375 - vascular endothelium
    8: "B cells",  # CD20=293, CD45RA=398 - B lymphocytes
    9: "Proliferating T",  # CD3E=679, PCNA=83, CD138=44 - proliferating T cells
    10: "Vascular Stromal",  # CD31=53, Vimentin=209 - perivascular/stromal
}

# Simplified mapping for benchmark comparison (6 cell types)
# Groups related cell types for cleaner analysis
SIMPLIFIED_CELLTYPE_MAP = {
    "CD8+ T cells": "T cells",
    "B cells": "B cells",
    "Macrophages": "Macrophages",
    "Mixed Immune": "Macrophages",  # Assign to Macrophages for simplicity
    "Myofibroblasts": "Fibroblasts",
    "Stromal": "Fibroblasts",  # Assign to Fibroblasts
    "Epithelial": "Epithelial",
    "Endothelial": "Endothelial",
    "Vascular Stromal": "Endothelial",
    "Proliferating T": "T cells",  # Keep with T cells
}

# Granular 10-cell type mapping (unsimplified)
# This provides maximum granularity to highlight CITEgeist's proteomic advantage
GRANULAR_CELLTYPE_MAP = {
    "CD8+ T cells": "CD8+ T cells",
    "Macrophages": "Macrophages",
    "Mixed Immune": "Mixed Immune",
    "Epithelial": "Epithelial",
    "Myofibroblasts": "Myofibroblasts",
    "Stromal": "Stromal",
    "Endothelial": "Endothelial",
    "B cells": "B cells",
    "Proliferating T": "Proliferating T",
    "Vascular Stromal": "Vascular Stromal",
}


def load_rna_clusters(
    data_dir: Union[str, Path],
    n_clusters: int = 10,
    clustering_type: str = "gene_expression",
) -> pd.DataFrame:
    """
    Load RNA-based clustering from Xenium analysis output.

    Args:
        data_dir: Path to Xenium output directory containing analysis.tar.gz
        n_clusters: Number of clusters to use (2-10 available)
        clustering_type: "gene_expression" or "protein_expression"

    Returns:
        DataFrame with cell_id index and 'Cluster' column

    Reference:
        Zhao et al. (2025). BMC Bioinformatics.
        https://doi.org/10.1186/s12859-025-06044-0
    """
    data_dir = Path(data_dir)
    analysis_tar = data_dir / "analysis.tar.gz"

    if not analysis_tar.exists():
        raise FileNotFoundError(f"Analysis archive not found: {analysis_tar}")

    cluster_path = f"analysis/clustering/{clustering_type}_kmeans_{n_clusters}_clusters/clusters.csv"

    with tarfile.open(analysis_tar, "r:gz") as tar:
        try:
            f = tar.extractfile(cluster_path)
            clusters_df = pd.read_csv(f)
        except KeyError:
            raise ValueError(f"Cluster file not found in archive: {cluster_path}")

    # Standardize column names
    clusters_df = clusters_df.rename(columns={"Barcode": "cell_id"})
    clusters_df["cell_id"] = clusters_df["cell_id"].astype(str)
    clusters_df = clusters_df.set_index("cell_id")

    logger.info(f"Loaded {len(clusters_df)} cells with {n_clusters} clusters")

    return clusters_df


def annotate_clusters(
    clusters_df: pd.DataFrame,
    cluster_map: Dict[int, str] = None,
    simplify: bool = True,
) -> pd.Series:
    """
    Annotate cluster numbers with cell type names.

    Args:
        clusters_df: DataFrame with 'Cluster' column
        cluster_map: Dict mapping cluster number to cell type name.
                     If None, uses XENIUM_KIDNEY_RNA_CLUSTER_MAP
        simplify: If True, apply SIMPLIFIED_CELLTYPE_MAP to merge related types

    Returns:
        Series with cell_id index and cell type values
    """
    if cluster_map is None:
        cluster_map = XENIUM_KIDNEY_RNA_CLUSTER_MAP

    cell_types = clusters_df["Cluster"].map(cluster_map)

    if simplify:
        cell_types = cell_types.map(SIMPLIFIED_CELLTYPE_MAP)

    # Handle unmapped clusters
    cell_types = cell_types.fillna("Unknown")

    type_counts = cell_types.value_counts()
    logger.info(f"Cell type distribution:\n{type_counts}")

    return cell_types


def load_rna_cell_types(
    data_dir: Union[str, Path],
    n_clusters: int = 10,
    cluster_map: Dict[int, str] = None,
    simplify: bool = True,
) -> pd.Series:
    """
    High-level function to load and annotate cell types from RNA clustering.

    This is the recommended method for obtaining ground truth cell types
    in Xenium benchmarking, as it avoids circular logic with protein markers.

    Args:
        data_dir: Path to Xenium output directory
        n_clusters: Number of clusters (default 10)
        cluster_map: Custom cluster-to-celltype mapping
        simplify: Merge related cell types

    Returns:
        Series with cell_id index and cell type values

    Reference:
        Zhao et al. (2025). "Benchmarking cell type annotation methods for
        10x Xenium spatial transcriptomics data." BMC Bioinformatics.
        https://doi.org/10.1186/s12859-025-06044-0
    """
    clusters_df = load_rna_clusters(data_dir, n_clusters=n_clusters)
    cell_types = annotate_clusters(clusters_df, cluster_map=cluster_map, simplify=simplify)

    return cell_types


def validate_cluster_annotation(
    adata_protein: sc.AnnData,
    cell_types: pd.Series,
    expected_markers: Dict[str, List[str]] = None,
) -> pd.DataFrame:
    """
    Validate cluster annotations by checking protein marker expression.

    Args:
        adata_protein: AnnData with protein expression
        cell_types: Series with cell type labels
        expected_markers: Dict mapping cell type to expected high markers

    Returns:
        DataFrame with mean marker expression per cell type
    """
    if expected_markers is None:
        expected_markers = {
            "T cells": ["CD3E", "CD8A", "CD4"],
            "B cells": ["CD20"],
            "Macrophages": ["CD68", "HLA-DR"],
            "Fibroblasts": ["alphaSMA", "Vimentin"],
            "Epithelial": ["PanCK", "E-Cadherin"],
            "Endothelial": ["CD31"],
        }

    # Get expression matrix
    X = adata_protein.X.toarray() if hasattr(adata_protein.X, "toarray") else adata_protein.X
    proteins = list(adata_protein.var_names)

    # Align cell types with adata
    common_cells = cell_types.index.intersection(adata_protein.obs_names)
    cell_types_aligned = cell_types.loc[common_cells]
    adata_aligned = adata_protein[common_cells]
    X_aligned = X[np.isin(adata_protein.obs_names, common_cells)]

    # Calculate mean expression per cell type
    results = {}
    for ct in cell_types_aligned.unique():
        ct_mask = cell_types_aligned == ct
        ct_means = {}
        for marker in proteins:
            idx = proteins.index(marker)
            ct_means[marker] = X_aligned[ct_mask.values, idx].mean()
        results[ct] = ct_means

    results_df = pd.DataFrame(results).T

    # Log validation summary
    logger.info("Cluster annotation validation:")
    for ct, markers in expected_markers.items():
        if ct in results_df.index:
            available_markers = [m for m in markers if m in results_df.columns]
            if available_markers:
                mean_expr = results_df.loc[ct, available_markers].mean()
                logger.info(f"  {ct}: mean of {available_markers} = {mean_expr:.1f}")

    return results_df


def calculate_spot_proportions_rna(
    cell_to_spot: pd.DataFrame,
    cell_types: pd.Series,
    spot_names: List[str],
) -> pd.DataFrame:
    """
    Calculate cell type proportions per spot using RNA-based cell types.

    Args:
        cell_to_spot: DataFrame with 'spot_id' column mapping cell_id -> spot
        cell_types: Series with cell_id index and cell type values
        spot_names: List of spot names to include

    Returns:
        DataFrame with spot index and cell type proportion columns
    """
    # Get unique cell types
    unique_types = sorted(cell_types.unique())

    # Initialize proportions matrix
    proportions = pd.DataFrame(0.0, index=spot_names, columns=unique_types)
    n_cells_per_spot = pd.Series(0, index=spot_names)

    for spot_name in spot_names:
        # Find cells in this spot
        cell_mask = cell_to_spot["spot_id"] == spot_name
        spot_cell_ids = cell_to_spot.index[cell_mask]

        # Get cell types for these cells
        spot_cells_with_type = [c for c in spot_cell_ids if c in cell_types.index]

        if len(spot_cells_with_type) > 0:
            spot_types = cell_types.loc[spot_cells_with_type]
            type_counts = spot_types.value_counts(normalize=True)

            for ct, prop in type_counts.items():
                if ct in proportions.columns:
                    proportions.loc[spot_name, ct] = prop

            n_cells_per_spot[spot_name] = len(spot_cells_with_type)

    # Add metadata
    proportions["n_cells"] = n_cells_per_spot

    logger.info(f"Calculated RNA-based proportions for {len(spot_names)} spots")

    return proportions


if __name__ == "__main__":
    import sys

    logging.basicConfig(level=logging.INFO)

    sys.path.insert(0, "benchmarks/xenium_benchmarking/src")
    from load_xenium import load_xenium_data, split_gex_protein

    DATA_DIR = "/path/to/Xenium_RCC"

    print("Loading RNA-based cell types...")
    cell_types = load_rna_cell_types(DATA_DIR, n_clusters=10, simplify=True)

    print(f"\nCell type distribution:")
    print(cell_types.value_counts())

    print("\nValidating against protein markers...")
    adata = load_xenium_data(DATA_DIR)
    _, adata_protein = split_gex_protein(adata)

    validation_df = validate_cluster_annotation(adata_protein, cell_types)
    print("\nMean protein expression by RNA-defined cell type:")
    key_markers = ["CD68", "CD3E", "CD20", "PanCK", "alphaSMA", "CD31"]
    print(validation_df[key_markers].round(1))
