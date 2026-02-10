#!/usr/bin/env python
"""
Generate marker gene lists for CARD reference-free mode.

For the protein_gt benchmark, we use the achievable-7 cell type definitions
from benchmark_constants.py and find the corresponding marker GENES (not proteins)
in the reference scRNA-seq data.

This script identifies highly differentially expressed genes for each cell type
from the reference scRNA-seq dataset.

Usage:
    python generate_marker_genes.py --reference reference.h5ad --output marker_genes.csv --n-markers 50
"""

import argparse
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


def find_marker_genes(adata, cell_type_col="cell_type", n_markers=50):
    """
    Find marker genes for each cell type using Wilcoxon rank-sum test.

    Args:
        adata: AnnData with cell type annotations
        cell_type_col: Column name for cell types
        n_markers: Number of top markers per cell type

    Returns:
        DataFrame with columns: cell_type, gene, score
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
        # Get top genes for this cell type
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
    parser = argparse.ArgumentParser(
        description="Generate marker genes for CARD reference-free mode"
    )
    parser.add_argument(
        "--reference",
        type=str,
        required=True,
        help="Reference scRNA-seq h5ad file"
    )
    parser.add_argument(
        "--output",
        type=str,
        required=True,
        help="Output CSV file"
    )
    parser.add_argument(
        "--cell-type-col",
        type=str,
        default="cell_type",
        help="Cell type column name"
    )
    parser.add_argument(
        "--n-markers",
        type=int,
        default=50,
        help="Number of marker genes per cell type"
    )

    args = parser.parse_args()

    # Load reference
    logger.info(f"Loading reference: {args.reference}")
    adata = sc.read_h5ad(args.reference)
    logger.info(f"  Shape: {adata.shape}")
    logger.info(f"  Cell types: {adata.obs[args.cell_type_col].unique().tolist()}")

    # Find markers
    marker_df = find_marker_genes(
        adata,
        cell_type_col=args.cell_type_col,
        n_markers=args.n_markers
    )

    # Save
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    marker_df.to_csv(output_path, index=False)
    logger.info(f"Saved {len(marker_df)} markers to {output_path}")

    # Summary
    for ct in marker_df["cell_type"].unique():
        n = len(marker_df[marker_df["cell_type"] == ct])
        logger.info(f"  {ct}: {n} markers")


if __name__ == "__main__":
    main()
