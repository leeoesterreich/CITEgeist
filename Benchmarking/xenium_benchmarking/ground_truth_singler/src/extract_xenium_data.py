#!/usr/bin/env python
"""
Extract cell-level GEX and protein data from Xenium cell_feature_matrix.h5.

Outputs:
- xenium_cell_gex.csv: cells x RNA genes (raw counts)
- xenium_cell_proteins.csv: cells x protein markers (raw counts)

Cell IDs are taken from cells.parquet to match cell_to_spot_mapping.csv format.
"""

import sys as _sys
from pathlib import Path as _Path

_sys.path.insert(0, str(_Path(__file__).resolve().parents[3]))
import _paths as _cg  # noqa: E402

import argparse
import logging
import h5py
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.sparse import csc_matrix

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(message)s")
logger = logging.getLogger(__name__)

XENIUM_DIR = str(_cg.root("DATA") / "10x_PublicData" / "Xenium_RNA_Proteomic_RenalCellCarcinoma")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load cell IDs from cells.parquet (these match cell_to_spot_mapping.csv index)
    cells_df = pd.read_parquet(Path(XENIUM_DIR) / "cells.parquet")
    cell_ids = cells_df["cell_id"].values
    logger.info(f"Loaded {len(cell_ids)} cell IDs from cells.parquet")

    # Load cell_feature_matrix.h5
    h5_path = Path(XENIUM_DIR) / "cell_feature_matrix.h5"
    with h5py.File(h5_path, "r") as f:
        # 10x h5 format: matrix group with data, indices, indptr, shape, features
        mat_group = f["matrix"]
        data = mat_group["data"][:]
        indices = mat_group["indices"][:]
        indptr = mat_group["indptr"][:]
        shape = mat_group["shape"][:]

        # Feature metadata
        feature_names = mat_group["features"]["name"][:].astype(str)
        feature_types = mat_group["features"]["feature_type"][:].astype(str)

    # Build sparse matrix (features x cells in 10x format)
    sparse_mat = csc_matrix((data, indices, indptr), shape=shape)
    logger.info(f"Loaded matrix: {sparse_mat.shape} (features x cells)")

    # Split by feature type
    rna_mask = np.array([ft == "Gene Expression" for ft in feature_types])
    protein_mask = np.array([ft != "Gene Expression" for ft in feature_types])

    rna_names = feature_names[rna_mask]
    protein_names = feature_names[protein_mask]

    logger.info(f"RNA features: {rna_mask.sum()}, Protein features: {protein_mask.sum()}")

    # Extract as dense (cells x features)
    rna_mat = sparse_mat[rna_mask, :].T.toarray()
    protein_mat = sparse_mat[protein_mask, :].T.toarray()

    # Create DataFrames with cell IDs as index
    gex_df = pd.DataFrame(rna_mat, index=cell_ids, columns=rna_names)
    protein_df = pd.DataFrame(protein_mat, index=cell_ids, columns=protein_names)

    # Save
    gex_path = output_dir / "xenium_cell_gex.csv"
    protein_path = output_dir / "xenium_cell_proteins.csv"
    gex_df.to_csv(gex_path)
    protein_df.to_csv(protein_path)
    logger.info(f"Saved GEX: {gex_df.shape} -> {gex_path}")
    logger.info(f"Saved Protein: {protein_df.shape} -> {protein_path}")


if __name__ == "__main__":
    main()
