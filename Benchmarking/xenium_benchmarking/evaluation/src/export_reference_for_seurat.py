#!/usr/bin/env python3
"""
Export GSE156632 reference h5ad to 10x-style MTX format for Seurat loading.

Avoids the SeuratDisk h5ad-to-h5seurat conversion that crashes with HDF5 errors.
Outputs: matrix.mtx.gz, barcodes.tsv.gz, features.tsv.gz, cell_types.csv
"""

import argparse
import gzip
import logging
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.io import mmwrite
from scipy.sparse import csc_matrix

logger = logging.getLogger(__name__)


def main(reference_h5ad, output_dir):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info(f'Loading reference: {reference_h5ad}')
    adata = sc.read_h5ad(reference_h5ad)
    logger.info(f'Reference: {adata.n_obs} cells x {adata.n_vars} genes')

    # Use raw counts if available
    if adata.raw is not None:
        logger.info('Using raw counts from adata.raw')
        X = adata.raw.X
        gene_names = list(adata.raw.var_names)
    else:
        X = adata.X
        gene_names = list(adata.var_names)

    # Seurat replaces underscores with dashes, which can create new duplicates.
    # Do the same replacement first, then make names unique.
    gene_names = [g.replace('_', '-') for g in gene_names]

    from collections import Counter
    name_counts = Counter(gene_names)
    dups = {n for n, c in name_counts.items() if c > 1}
    if dups:
        logger.info(f'Making {len(dups)} duplicate gene names unique after underscore replacement')
        seen = {}
        unique_names = []
        for name in gene_names:
            if name in dups:
                idx = seen.get(name, 0)
                seen[name] = idx + 1
                unique_names.append(f'{name}.{idx}' if idx > 0 else name)
            else:
                unique_names.append(name)
        gene_names = unique_names

    # Ensure sparse CSC for mmwrite
    if not hasattr(X, 'toarray'):
        from scipy.sparse import csc_matrix as csc
        X = csc(X)
    else:
        X = csc_matrix(X)

    cell_ids = list(adata.obs_names)

    # Write matrix.mtx.gz (genes x cells, transposed from anndata's cells x genes)
    mtx_path = output_dir / 'matrix.mtx'
    logger.info(f'Writing matrix ({X.shape[0]} cells x {X.shape[1]} genes, transposed for 10x format)...')
    mmwrite(str(mtx_path), X.T)  # 10x format is genes x cells

    # Gzip the mtx file
    with open(mtx_path, 'rb') as f_in:
        with gzip.open(str(mtx_path) + '.gz', 'wb') as f_out:
            f_out.write(f_in.read())
    mtx_path.unlink()
    logger.info('Wrote matrix.mtx.gz')

    # Write barcodes.tsv.gz
    barcodes_path = output_dir / 'barcodes.tsv.gz'
    with gzip.open(str(barcodes_path), 'wt') as f:
        for cell_id in cell_ids:
            f.write(f'{cell_id}\n')
    logger.info(f'Wrote barcodes.tsv.gz ({len(cell_ids)} cells)')

    # Write features.tsv.gz (gene_id \t gene_name \t feature_type)
    features_path = output_dir / 'features.tsv.gz'
    with gzip.open(str(features_path), 'wt') as f:
        for gene in gene_names:
            f.write(f'{gene}\t{gene}\tGene Expression\n')
    logger.info(f'Wrote features.tsv.gz ({len(gene_names)} genes)')

    # Write cell_types.csv
    ct_col = None
    for col in ['cell_type', 'celltype', 'CellType', 'cell_type_broad']:
        if col in adata.obs.columns:
            ct_col = col
            break
    if ct_col is None:
        logger.warning('No cell_type column found in obs. Available: %s', list(adata.obs.columns))
        ct_col = adata.obs.columns[0]
        logger.warning(f'Using first column: {ct_col}')

    ct_df = pd.DataFrame({
        'cell_id': cell_ids,
        'cell_type': adata.obs[ct_col].values
    })
    ct_df.to_csv(output_dir / 'cell_types.csv', index=False)
    logger.info(f'Wrote cell_types.csv using column "{ct_col}"')
    logger.info('Cell type distribution:')
    for ct, n in ct_df['cell_type'].value_counts().items():
        logger.info(f'  {ct}: {n}')

    logger.info('Done. Files ready for Seurat Read10X()')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--reference', type=str, required=True,
                        help='Path to reference h5ad file')
    parser.add_argument('--output_dir', type=str, required=True,
                        help='Output directory for 10x-format files')
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s - %(levelname)s - %(message)s')
    main(args.reference, args.output_dir)
