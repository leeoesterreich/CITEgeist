"""
Process GSE156632 scRNA-seq reference data for Xenium benchmarking.

This script:
1. Loads all 12 CSV.gz files from GSE156632
2. Merges into single AnnData object with batch information
3. Performs QC filtering (min_genes, max_genes, mt%)
4. Normalizes (log1p, target_sum=1e4)
5. Integrates batches using Harmony
6. Clusters and annotates 6 cell types matching Xenium ground truth
7. Exports in formats for Cell2Location, Tangram, RCTD, Seurat

Target cell types (must match Xenium RNA-based ground truth):
- B cells
- T cells
- Macrophages
- Fibroblasts
- Epithelial
- Endothelial

Usage:
    python process_reference.py --input-dir /path/to/GSE156632 --output-dir /path/to/processed
"""

import os
import sys
import argparse
import logging
from pathlib import Path
from typing import Dict, List, Optional
import gc

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad

# Try to import harmony - will be available in the cell2location_env
try:
    import harmonypy
    HARMONY_AVAILABLE = True
except ImportError:
    HARMONY_AVAILABLE = False
    print("Warning: harmonypy not available. Will skip batch integration.")

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)


# Marker genes for cell type annotation
# These markers are used to identify the 6 cell types matching Xenium ground truth
MARKER_GENES = {
    "T cells": ["CD3D", "CD3E", "CD3G", "CD4", "CD8A", "CD8B", "TRAC", "TRBC1", "TRBC2"],
    "B cells": ["CD79A", "CD79B", "MS4A1", "CD19", "PAX5", "BANK1"],
    "Macrophages": ["CD68", "CD163", "CSF1R", "MARCO", "MSR1", "CD14", "FCGR3A"],
    "Fibroblasts": ["COL1A1", "COL1A2", "DCN", "VIM", "ACTA2", "FAP", "PDGFRA"],
    "Epithelial": ["EPCAM", "KRT18", "KRT8", "CA9", "PAX8", "CDH1", "KRT19"],
    "Endothelial": ["PECAM1", "VWF", "CDH5", "CD34", "KDR", "FLT1"],
}


def load_csv_files(input_dir: Path) -> ad.AnnData:
    """
    Load all GSE156632 CSV.gz files and merge into single AnnData.

    Args:
        input_dir: Directory containing CSV.gz files

    Returns:
        Merged AnnData object with batch information
    """
    csv_files = sorted(input_dir.glob("GSM*.csv.gz"))
    logger.info(f"Found {len(csv_files)} CSV files")

    adatas = []
    for csv_file in csv_files:
        logger.info(f"Loading {csv_file.name}...")

        # Extract sample name from filename
        # Format: GSM4735364_RCC1n.csv.gz -> RCC1n
        sample_name = csv_file.stem.replace(".csv", "").split("_")[1]

        # Determine if tumor or normal
        tissue_type = "normal" if sample_name.endswith("n") else "tumor"

        # Load CSV (genes as rows, cells as columns - need to transpose)
        # The CSV has Gene_ID, Symbol, then cell barcodes
        # We want to use Symbol (gene names) as var_names
        df = pd.read_csv(csv_file, index_col=0)
        logger.info(f"  Shape before processing: {df.shape}")

        # Extract gene symbols from the 'Symbol' column (first column after index)
        gene_symbols = df['Symbol'].values if 'Symbol' in df.columns else df.index.values

        # Drop the Symbol column if it exists (keep only count data)
        if 'Symbol' in df.columns:
            df = df.drop(columns=['Symbol'])

        logger.info(f"  Shape after dropping Symbol column: {df.shape}")

        # Transpose: cells should be rows, genes should be columns
        df = df.T
        logger.info(f"  Shape after transpose: {df.shape}")

        # Create AnnData
        adata = ad.AnnData(X=df.values.astype(np.float32))
        adata.obs_names = df.index.astype(str)
        adata.var_names = pd.Index(gene_symbols).astype(str)

        # Make var_names unique (some gene symbols may be duplicated)
        adata.var_names_make_unique()

        # Add metadata
        adata.obs["sample"] = sample_name
        adata.obs["tissue_type"] = tissue_type
        adata.obs["batch"] = sample_name  # For batch correction

        logger.info(f"  Created AnnData: {adata.shape[0]} cells, {adata.shape[1]} genes")
        adatas.append(adata)

        # Clean up
        del df
        gc.collect()

    # Merge all samples
    logger.info("Merging all samples...")
    adata = ad.concat(adatas, join="outer", fill_value=0)
    adata.obs_names_make_unique()

    logger.info(f"Merged AnnData: {adata.shape[0]} cells, {adata.shape[1]} genes")
    logger.info(f"Samples: {adata.obs['sample'].unique().tolist()}")

    # Clean up
    del adatas
    gc.collect()

    return adata


def qc_filtering(
    adata: ad.AnnData,
    min_genes: int = 200,
    max_genes: int = 5000,
    max_mt_pct: float = 20.0,
    min_cells: int = 3,
) -> ad.AnnData:
    """
    Quality control filtering.

    Args:
        adata: Input AnnData
        min_genes: Minimum genes per cell
        max_genes: Maximum genes per cell
        max_mt_pct: Maximum mitochondrial percentage
        min_cells: Minimum cells per gene

    Returns:
        Filtered AnnData
    """
    logger.info("Performing QC filtering...")
    n_cells_before = adata.n_obs
    n_genes_before = adata.n_vars

    # Calculate QC metrics
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
    )

    # Filter cells
    sc.pp.filter_cells(adata, min_genes=min_genes)
    logger.info(f"  After min_genes={min_genes}: {adata.n_obs} cells")

    # Filter by max genes
    adata = adata[adata.obs["n_genes_by_counts"] < max_genes, :].copy()
    logger.info(f"  After max_genes={max_genes}: {adata.n_obs} cells")

    # Filter by MT percentage
    adata = adata[adata.obs["pct_counts_mt"] < max_mt_pct, :].copy()
    logger.info(f"  After max_mt_pct={max_mt_pct}: {adata.n_obs} cells")

    # Filter genes
    sc.pp.filter_genes(adata, min_cells=min_cells)
    logger.info(f"  After min_cells={min_cells}: {adata.n_vars} genes")

    logger.info(f"QC complete: {n_cells_before} -> {adata.n_obs} cells, "
                f"{n_genes_before} -> {adata.n_vars} genes")

    return adata


def normalize_and_scale(adata: ad.AnnData, target_sum: float = 1e4) -> ad.AnnData:
    """
    Normalize, log-transform, and identify highly variable genes.

    Args:
        adata: QC-filtered AnnData
        target_sum: Target sum for normalization

    Returns:
        Normalized AnnData with raw counts in .raw
    """
    logger.info("Normalizing and scaling...")

    # Store raw counts
    adata.raw = adata.copy()

    # Normalize
    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)

    # Find highly variable genes
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=2000,
        batch_key="batch",
        flavor="seurat_v3",
        subset=False,
    )

    n_hvg = adata.var["highly_variable"].sum()
    logger.info(f"  Identified {n_hvg} highly variable genes")

    return adata


def run_harmony(adata: ad.AnnData, batch_key: str = "batch") -> ad.AnnData:
    """
    Run Harmony batch integration.

    Args:
        adata: Normalized AnnData
        batch_key: Column in obs for batch information

    Returns:
        AnnData with harmony-corrected embedding
    """
    if not HARMONY_AVAILABLE:
        logger.warning("Harmony not available, skipping batch integration")
        return adata

    logger.info("Running Harmony batch integration...")

    # PCA first
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=50, use_highly_variable=True)

    # Run Harmony
    ho = harmonypy.run_harmony(
        adata.obsm["X_pca"],
        adata.obs,
        batch_key,
        max_iter_harmony=20,
    )

    adata.obsm["X_pca_harmony"] = ho.Z_corr.T

    # Compute neighbors and UMAP on harmony-corrected embedding
    sc.pp.neighbors(adata, use_rep="X_pca_harmony", n_neighbors=15)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=0.8)

    logger.info(f"  Found {adata.obs['leiden'].nunique()} clusters")

    return adata


def annotate_cell_types(
    adata: ad.AnnData,
    marker_genes: Dict[str, List[str]] = MARKER_GENES,
) -> ad.AnnData:
    """
    Annotate cell types based on marker gene expression.

    Uses a scoring approach: for each cell, compute the mean expression
    of marker genes for each cell type. Assign the cell type with highest score.

    Args:
        adata: Clustered AnnData
        marker_genes: Dictionary mapping cell type to marker genes

    Returns:
        AnnData with cell_type annotation
    """
    logger.info("Annotating cell types...")

    # Get normalized expression matrix
    if adata.raw is not None:
        expr_df = pd.DataFrame(
            adata.raw.X.toarray() if hasattr(adata.raw.X, 'toarray') else adata.raw.X,
            index=adata.obs_names,
            columns=adata.raw.var_names,
        )
    else:
        expr_df = pd.DataFrame(
            adata.X.toarray() if hasattr(adata.X, 'toarray') else adata.X,
            index=adata.obs_names,
            columns=adata.var_names,
        )

    # Score each cell type
    scores = pd.DataFrame(index=adata.obs_names)

    for cell_type, markers in marker_genes.items():
        # Find markers present in data
        present_markers = [m for m in markers if m in expr_df.columns]
        logger.info(f"  {cell_type}: {len(present_markers)}/{len(markers)} markers present")

        if present_markers:
            scores[cell_type] = expr_df[present_markers].mean(axis=1)
        else:
            scores[cell_type] = 0.0

    # Assign cell type based on highest score
    adata.obs["cell_type"] = scores.idxmax(axis=1)

    # Add scores to obs
    for cell_type in marker_genes.keys():
        adata.obs[f"{cell_type}_score"] = scores[cell_type].values

    # Log cell type distribution
    ct_counts = adata.obs["cell_type"].value_counts()
    logger.info("Cell type distribution:")
    for ct, count in ct_counts.items():
        pct = 100 * count / adata.n_obs
        logger.info(f"  {ct}: {count} ({pct:.1f}%)")

    return adata


def export_for_cell2location(
    adata: ad.AnnData,
    output_dir: Path,
) -> None:
    """Export in format for Cell2Location reference training."""
    logger.info("Exporting for Cell2Location...")

    c2l_dir = output_dir / "cell2location"
    c2l_dir.mkdir(parents=True, exist_ok=True)

    # Cell2Location needs raw counts + cell type labels
    adata_export = adata.raw.to_adata() if adata.raw is not None else adata.copy()
    adata_export.obs["cell_type"] = adata.obs["cell_type"]
    adata_export.obs["sample"] = adata.obs["sample"]
    adata_export.obs["batch"] = adata.obs["batch"]

    adata_export.write_h5ad(c2l_dir / "reference.h5ad")
    logger.info(f"  Saved to {c2l_dir / 'reference.h5ad'}")


def export_for_tangram(
    adata: ad.AnnData,
    output_dir: Path,
) -> None:
    """Export in format for Tangram."""
    logger.info("Exporting for Tangram...")

    tg_dir = output_dir / "tangram"
    tg_dir.mkdir(parents=True, exist_ok=True)

    # Tangram needs normalized counts + cell type labels
    adata_export = adata.copy()
    adata_export.write_h5ad(tg_dir / "reference.h5ad")
    logger.info(f"  Saved to {tg_dir / 'reference.h5ad'}")


def export_for_rctd(
    adata: ad.AnnData,
    output_dir: Path,
) -> None:
    """Export in format for RCTD (R)."""
    logger.info("Exporting for RCTD...")

    rctd_dir = output_dir / "rctd"
    rctd_dir.mkdir(parents=True, exist_ok=True)

    # RCTD needs:
    # 1. counts matrix (genes x cells) as CSV
    # 2. cell types as CSV

    # Get raw counts
    if adata.raw is not None:
        counts = adata.raw.X.toarray() if hasattr(adata.raw.X, 'toarray') else adata.raw.X
        gene_names = adata.raw.var_names
    else:
        counts = adata.X.toarray() if hasattr(adata.X, 'toarray') else adata.X
        gene_names = adata.var_names

    # Save counts (genes x cells - transpose)
    counts_df = pd.DataFrame(
        counts.T,
        index=gene_names,
        columns=adata.obs_names,
    )
    counts_df.to_csv(rctd_dir / "counts.csv")
    logger.info(f"  Saved counts: {counts_df.shape}")

    # Save cell types
    cell_types = adata.obs[["cell_type"]].copy()
    cell_types.to_csv(rctd_dir / "cell_types.csv")
    logger.info(f"  Saved cell types: {cell_types.shape}")

    # Also save as h5ad for convenience
    adata_export = adata.raw.to_adata() if adata.raw is not None else adata.copy()
    adata_export.obs["cell_type"] = adata.obs["cell_type"]
    adata_export.write_h5ad(rctd_dir / "reference.h5ad")


def export_for_seurat(
    adata: ad.AnnData,
    output_dir: Path,
) -> None:
    """Export in format for Seurat."""
    logger.info("Exporting for Seurat...")

    seurat_dir = output_dir / "seurat"
    seurat_dir.mkdir(parents=True, exist_ok=True)

    # Seurat needs same format as RCTD
    # Get raw counts
    if adata.raw is not None:
        counts = adata.raw.X.toarray() if hasattr(adata.raw.X, 'toarray') else adata.raw.X
        gene_names = adata.raw.var_names
    else:
        counts = adata.X.toarray() if hasattr(adata.X, 'toarray') else adata.X
        gene_names = adata.var_names

    # Save counts (genes x cells - transpose)
    counts_df = pd.DataFrame(
        counts.T,
        index=gene_names,
        columns=adata.obs_names,
    )
    counts_df.to_csv(seurat_dir / "counts.csv")
    logger.info(f"  Saved counts: {counts_df.shape}")

    # Save cell types
    cell_types = adata.obs[["cell_type"]].copy()
    cell_types.to_csv(seurat_dir / "cell_types.csv")
    logger.info(f"  Saved cell types: {cell_types.shape}")

    # Also save as h5ad for convenience
    adata_export = adata.raw.to_adata() if adata.raw is not None else adata.copy()
    adata_export.obs["cell_type"] = adata.obs["cell_type"]
    adata_export.write_h5ad(seurat_dir / "reference.h5ad")


def main():
    parser = argparse.ArgumentParser(
        description="Process GSE156632 reference for Xenium benchmarking"
    )
    parser.add_argument(
        "--input-dir",
        type=str,
        required=True,
        help="Directory containing GSE156632 CSV.gz files",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        required=True,
        help="Output directory for processed data",
    )
    parser.add_argument(
        "--min-genes",
        type=int,
        default=200,
        help="Minimum genes per cell (default: 200)",
    )
    parser.add_argument(
        "--max-genes",
        type=int,
        default=5000,
        help="Maximum genes per cell (default: 5000)",
    )
    parser.add_argument(
        "--max-mt-pct",
        type=float,
        default=20.0,
        help="Maximum mitochondrial percentage (default: 20.0)",
    )
    parser.add_argument(
        "--skip-harmony",
        action="store_true",
        help="Skip Harmony batch integration",
    )

    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("=" * 60)
    logger.info("GSE156632 Reference Processing")
    logger.info("=" * 60)
    logger.info(f"Input directory: {input_dir}")
    logger.info(f"Output directory: {output_dir}")

    # Step 1: Load data
    logger.info("\n" + "=" * 60)
    logger.info("Step 1: Loading CSV files")
    logger.info("=" * 60)
    adata = load_csv_files(input_dir)

    # Step 2: QC filtering
    logger.info("\n" + "=" * 60)
    logger.info("Step 2: QC Filtering")
    logger.info("=" * 60)
    adata = qc_filtering(
        adata,
        min_genes=args.min_genes,
        max_genes=args.max_genes,
        max_mt_pct=args.max_mt_pct,
    )

    # Step 3: Normalize
    logger.info("\n" + "=" * 60)
    logger.info("Step 3: Normalization")
    logger.info("=" * 60)
    adata = normalize_and_scale(adata)

    # Step 4: Harmony batch integration
    if not args.skip_harmony:
        logger.info("\n" + "=" * 60)
        logger.info("Step 4: Harmony Batch Integration")
        logger.info("=" * 60)
        adata = run_harmony(adata)
    else:
        logger.info("\n" + "=" * 60)
        logger.info("Step 4: Skipping Harmony (--skip-harmony)")
        logger.info("=" * 60)
        # Still need to do PCA, neighbors, clustering
        sc.pp.scale(adata, max_value=10)
        sc.tl.pca(adata, n_comps=50, use_highly_variable=True)
        sc.pp.neighbors(adata, n_neighbors=15)
        sc.tl.umap(adata)
        sc.tl.leiden(adata, resolution=0.8)

    # Step 5: Cell type annotation
    logger.info("\n" + "=" * 60)
    logger.info("Step 5: Cell Type Annotation")
    logger.info("=" * 60)
    adata = annotate_cell_types(adata)

    # Step 6: Save combined dataset
    logger.info("\n" + "=" * 60)
    logger.info("Step 6: Saving Combined Dataset")
    logger.info("=" * 60)
    adata.write_h5ad(output_dir / "GSE156632_combined.h5ad")
    logger.info(f"Saved combined dataset to {output_dir / 'GSE156632_combined.h5ad'}")

    # Step 7: Export for each method
    logger.info("\n" + "=" * 60)
    logger.info("Step 7: Exporting for Deconvolution Methods")
    logger.info("=" * 60)

    export_for_cell2location(adata, output_dir)
    export_for_tangram(adata, output_dir)
    export_for_rctd(adata, output_dir)
    export_for_seurat(adata, output_dir)

    # Summary
    logger.info("\n" + "=" * 60)
    logger.info("Processing Complete!")
    logger.info("=" * 60)
    logger.info(f"Final dataset: {adata.n_obs} cells, {adata.n_vars} genes")
    logger.info(f"Cell types: {adata.obs['cell_type'].unique().tolist()}")
    logger.info(f"Output files:")
    logger.info(f"  - {output_dir / 'GSE156632_combined.h5ad'}")
    logger.info(f"  - {output_dir / 'cell2location' / 'reference.h5ad'}")
    logger.info(f"  - {output_dir / 'tangram' / 'reference.h5ad'}")
    logger.info(f"  - {output_dir / 'rctd' / 'counts.csv'}")
    logger.info(f"  - {output_dir / 'seurat' / 'counts.csv'}")


if __name__ == "__main__":
    main()
