#!/usr/bin/env python
"""
Generate marker gene lists for CARD reference-free mode on simulation data.

Uses the simulation reference scRNA-seq data to identify differentially
expressed marker genes for each cell type.

Usage:
    python generate_marker_genes_simulation.py

Output:
    marker_genes_simulation.csv
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
REF_PATH = "/ix1/alee/LO_LAB/Personal/Brent_Schlegel/bts76/Projects/CITEgeist/Wu_Visium/Simulations/models/Wu_scRNA_ref_ERpos.h5ad"
OUTPUT_PATH = Path(__file__).parent / "marker_genes_simulation.csv"
N_MARKERS = 50


def find_marker_genes(adata, cell_type_col="celltype_major", n_markers=50):
    """
    Find marker genes for each cell type using Wilcoxon rank-sum test.
    """
    logger.info(f"Finding {n_markers} marker genes per cell type...")

    # Run differential expression
    sc.tl.rank_genes_groups(
        adata,
        groupby=cell_type_col,
        method="wilcoxon",
        n_genes=n_markers,
        use_raw=True if adata.raw is not None else False
    )

    # Extract results
    results = []
    for ct in adata.obs[cell_type_col].unique():
        de_results = sc.get.rank_genes_groups_df(adata, group=ct)
        top_genes = de_results.head(n_markers)

        for _, row in top_genes.iterrows():
            results.append({
                "cell_type": ct,
                "gene": row["names"],
                "score": row["scores"],
                "pval_adj": row["pvals_adj"]
            })

    return pd.DataFrame(results)


def main():
    logger.info(f"Loading reference: {REF_PATH}")
    adata = sc.read_h5ad(REF_PATH)
    logger.info(f"  Shape: {adata.shape}")
    logger.info(f"  Cell types: {adata.obs['celltype_major'].unique().tolist()}")

    # Find markers
    marker_df = find_marker_genes(
        adata,
        cell_type_col="celltype_major",
        n_markers=N_MARKERS
    )

    # Save
    marker_df.to_csv(OUTPUT_PATH, index=False)
    logger.info(f"Saved {len(marker_df)} markers to {OUTPUT_PATH}")

    # Summary
    for ct in marker_df["cell_type"].unique():
        n = len(marker_df[marker_df["cell_type"] == ct])
        logger.info(f"  {ct}: {n} markers")


if __name__ == "__main__":
    main()
