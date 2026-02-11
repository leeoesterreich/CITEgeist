#!/usr/bin/env python
"""
Convert simulation reference h5ad to CSV format for CARD R scripts.

This avoids the need for SeuratDisk (which conflicts with other R packages).

Usage:
    python convert_reference_to_csv.py
"""

import logging
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)

# Paths
REF_H5AD = Path("/ix1/alee/LO_LAB/Personal/Brent_Schlegel/bts76/Projects/CITEgeist/Wu_Visium/Simulations/models/Wu_scRNA_ref_ERpos.h5ad")
OUTPUT_DIR = Path(__file__).parent / "reference_csv"


def main():
    OUTPUT_DIR.mkdir(exist_ok=True)

    logger.info(f"Loading reference: {REF_H5AD}")
    adata = sc.read_h5ad(REF_H5AD)
    logger.info(f"  Shape: {adata.shape}")
    logger.info(f"  Cell types: {adata.obs['celltype_major'].unique().tolist()}")

    # Get counts (use raw if available, otherwise use X)
    if adata.raw is not None:
        logger.info("Using raw counts")
        counts = adata.raw.X
        gene_names = adata.raw.var_names
    else:
        logger.info("Using X matrix")
        counts = adata.X
        gene_names = adata.var_names

    # Convert to dense if sparse
    if hasattr(counts, 'toarray'):
        counts = counts.toarray()

    # Check if data is log-transformed (if max is small, likely log)
    if counts.max() < 50:
        logger.info("Data appears log-transformed, converting to counts: exp(X) - 1")
        counts = np.round(np.expm1(counts))

    # Create counts DataFrame (genes x cells, as expected by CARD)
    counts_df = pd.DataFrame(
        counts.T,
        index=gene_names,
        columns=adata.obs_names
    )

    # Save counts
    counts_path = OUTPUT_DIR / "reference_counts.csv"
    counts_df.to_csv(counts_path)
    logger.info(f"Saved counts: {counts_df.shape} to {counts_path}")

    # Save cell types
    cell_types_df = pd.DataFrame({
        'cell_type': adata.obs['celltype_major']
    }, index=adata.obs_names)

    cell_types_path = OUTPUT_DIR / "reference_cell_types.csv"
    cell_types_df.to_csv(cell_types_path)
    logger.info(f"Saved cell types: {cell_types_df.shape} to {cell_types_path}")

    # Summary
    logger.info("\nReference conversion complete!")
    logger.info(f"  Genes: {len(gene_names)}")
    logger.info(f"  Cells: {len(adata.obs_names)}")
    logger.info(f"  Cell types: {adata.obs['celltype_major'].nunique()}")


if __name__ == "__main__":
    main()
