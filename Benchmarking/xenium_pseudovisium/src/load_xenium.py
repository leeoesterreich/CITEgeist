"""
Load Xenium data into AnnData format.

This module provides functions to load Xenium single-cell data from
cell_feature_matrix.h5 and cells.parquet files into AnnData objects.
"""

import json
import logging
from pathlib import Path
from typing import Tuple, Dict, Optional

import h5py
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

logger = logging.getLogger(__name__)


def load_xenium_data(
    data_dir: str,
    min_counts: int = 0,
    min_cells: int = 0,
) -> sc.AnnData:
    """
    Load Xenium cell feature matrix and metadata into AnnData.

    Args:
        data_dir: Path to Xenium output directory containing
                  cell_feature_matrix.h5 and cells.parquet
        min_counts: Minimum total counts per cell to keep
        min_cells: Minimum cells per feature to keep

    Returns:
        AnnData with:
        - X: sparse matrix (cells x features)
        - obs: cell metadata including spatial coordinates
        - var: feature metadata with feature_types column
        - obsm['spatial']: (n_cells, 2) spatial coordinates
    """
    data_dir = Path(data_dir)
    h5_path = data_dir / "cell_feature_matrix.h5"
    parquet_path = data_dir / "cells.parquet"

    logger.info(f"Loading Xenium data from {data_dir}")

    # Load expression matrix from h5
    adata = _load_h5_matrix(h5_path)
    logger.info(f"Loaded matrix: {adata.shape[0]} cells x {adata.shape[1]} features")

    # Load cell metadata
    cells_df = pd.read_parquet(parquet_path)
    logger.info(f"Loaded cell metadata: {len(cells_df)} cells")

    # Align cell order
    adata.obs_names = adata.obs_names.astype(str)
    cells_df["cell_id"] = cells_df["cell_id"].astype(str)
    cells_df = cells_df.set_index("cell_id")

    # Ensure same cells in same order
    common_cells = adata.obs_names.intersection(cells_df.index)
    if len(common_cells) < len(adata.obs_names):
        logger.warning(
            f"Some cells in matrix not in metadata: "
            f"{len(adata.obs_names) - len(common_cells)} missing"
        )
    adata = adata[common_cells, :].copy()
    cells_df = cells_df.loc[common_cells]

    # Add cell metadata to obs
    for col in cells_df.columns:
        adata.obs[col] = cells_df[col].values

    # Set spatial coordinates
    adata.obsm["spatial"] = np.column_stack(
        [adata.obs["x_centroid"].values, adata.obs["y_centroid"].values]
    )

    # Basic filtering
    if min_counts > 0:
        sc.pp.filter_cells(adata, min_counts=min_counts)
    if min_cells > 0:
        sc.pp.filter_genes(adata, min_cells=min_cells)

    logger.info(f"Final shape after filtering: {adata.shape}")

    return adata


def _load_h5_matrix(h5_path: Path) -> sc.AnnData:
    """
    Load sparse matrix from Xenium cell_feature_matrix.h5.

    The h5 file has 10x-like structure:
    - matrix/data: non-zero values
    - matrix/indices: row indices
    - matrix/indptr: column pointers
    - matrix/barcodes: cell barcodes
    - matrix/features/name: feature names
    - matrix/features/feature_type: Gene Expression, Protein Expression, etc.
    """
    with h5py.File(h5_path, "r") as f:
        # Load sparse matrix components
        data = f["matrix/data"][:]
        indices = f["matrix/indices"][:]
        indptr = f["matrix/indptr"][:]
        shape = f["matrix/shape"][:]

        # Create sparse matrix (CSC format from 10x)
        matrix = sparse.csc_matrix((data, indices, indptr), shape=shape)
        # Transpose to cells x features
        matrix = matrix.T.tocsr()

        # Load barcodes (cell IDs)
        barcodes = [b.decode("utf-8") for b in f["matrix/barcodes"][:]]

        # Load feature info
        feature_names = [n.decode("utf-8") for n in f["matrix/features/name"][:]]
        feature_types = [t.decode("utf-8") for t in f["matrix/features/feature_type"][:]]
        feature_ids = [i.decode("utf-8") for i in f["matrix/features/id"][:]]

    # Create AnnData
    adata = sc.AnnData(X=matrix)
    adata.obs_names = pd.Index(barcodes)
    adata.var_names = pd.Index(feature_names)
    adata.var["feature_id"] = feature_ids
    adata.var["feature_types"] = feature_types

    return adata


def split_gex_protein(
    adata: sc.AnnData,
) -> Tuple[sc.AnnData, sc.AnnData]:
    """
    Split combined AnnData into separate GEX and protein objects.

    Args:
        adata: Combined AnnData with both gene expression and protein data

    Returns:
        Tuple of (adata_gex, adata_protein)
    """
    # Identify feature types
    gex_mask = adata.var["feature_types"] == "Gene Expression"
    protein_mask = adata.var["feature_types"] == "Protein Expression"

    logger.info(
        f"Splitting: {gex_mask.sum()} genes, {protein_mask.sum()} proteins"
    )

    # Split into separate objects
    adata_gex = adata[:, gex_mask].copy()
    adata_protein = adata[:, protein_mask].copy()

    # Update feature_types for CITEgeist compatibility
    # CITEgeist expects "Antibody Capture" not "Protein Expression"
    adata_protein.var["feature_types"] = "Antibody Capture"

    return adata_gex, adata_protein


def load_protein_panel(panel_json_path: str) -> Dict[str, Dict]:
    """
    Parse protein_panel.json to get protein metadata.

    Args:
        panel_json_path: Path to protein_panel.json

    Returns:
        Dict mapping protein short_name to metadata dict
    """
    with open(panel_json_path, "r") as f:
        panel = json.load(f)

    proteins = {}
    for target in panel["payload"]["panel"]["targets"]:
        short_name = target.get("short_name", "")
        if short_name and short_name != "n/a":
            proteins[short_name] = {
                "long_name": target.get("long_name", ""),
                "subpanel": target.get("from_subpanel", ""),
                "public_id": target.get("public_id", ""),
            }

    logger.info(f"Loaded {len(proteins)} proteins from panel")
    return proteins


def get_feature_summary(adata: sc.AnnData) -> pd.DataFrame:
    """
    Get summary of features by type.

    Args:
        adata: AnnData object

    Returns:
        DataFrame with feature type counts
    """
    return adata.var["feature_types"].value_counts().to_frame("count")


if __name__ == "__main__":
    # Test loading
    logging.basicConfig(level=logging.INFO)

    DATA_DIR = "/ix1/alee/LO_LAB/General/Public_Data/10x_PublicData/Xenium_RNA_Proteomic_RenalCellCarcinoma"

    print("Loading Xenium data...")
    adata = load_xenium_data(DATA_DIR)
    print(f"Loaded {adata.shape[0]} cells x {adata.shape[1]} features")

    print("\nFeature summary:")
    print(get_feature_summary(adata))

    print("\nSpatial coordinate ranges:")
    spatial = adata.obsm["spatial"]
    print(f"X: {spatial[:, 0].min():.1f} to {spatial[:, 0].max():.1f}")
    print(f"Y: {spatial[:, 1].min():.1f} to {spatial[:, 1].max():.1f}")

    print("\nSplitting into GEX and protein...")
    adata_gex, adata_protein = split_gex_protein(adata)
    print(f"GEX: {adata_gex.shape}")
    print(f"Protein: {adata_protein.shape}")
    print(f"Proteins: {list(adata_protein.var_names)}")
