"""
Convert h5ad files to CSV format for RCTD (R).

RCTD requires:
- Counts matrix (genes x cells/spots)
- Coordinates for spatial data
- Cell type annotations for reference

Usage:
    python convert_h5ad.py --input /path/to/file.h5ad --output-dir /path/to/output --mode spatial
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


def convert_spatial(adata, output_dir: Path, prefix: str = "spatial"):
    """
    Convert spatial h5ad to CSVs for RCTD.

    Args:
        adata: Spatial AnnData
        output_dir: Output directory
        prefix: Filename prefix
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info(f"Converting spatial data: {adata.shape}")

    # Get counts (genes x spots)
    if hasattr(adata.X, 'toarray'):
        counts = adata.X.toarray()
    else:
        counts = adata.X

    counts_df = pd.DataFrame(
        counts.T,  # Transpose to genes x spots
        index=adata.var_names,
        columns=adata.obs_names,
    )

    # Save counts
    counts_path = output_dir / f"{prefix}_counts.csv"
    counts_df.to_csv(counts_path)
    logger.info(f"  Saved counts: {counts_df.shape} to {counts_path}")

    # Get coordinates
    if 'spatial' in adata.obsm:
        coords = pd.DataFrame(
            adata.obsm['spatial'],
            index=adata.obs_names,
            columns=['x', 'y'],
        )
    else:
        # Create dummy coordinates
        logger.warning("No spatial coordinates found, creating dummy coordinates")
        n_spots = adata.n_obs
        side = int(np.ceil(np.sqrt(n_spots)))
        x = np.tile(np.arange(side), side)[:n_spots]
        y = np.repeat(np.arange(side), side)[:n_spots]
        coords = pd.DataFrame({'x': x, 'y': y}, index=adata.obs_names)

    coords_path = output_dir / f"{prefix}_coords.csv"
    coords.to_csv(coords_path)
    logger.info(f"  Saved coords: {coords.shape} to {coords_path}")

    return counts_path, coords_path


def convert_reference(adata, output_dir: Path, cell_type_col: str = "cell_type"):
    """
    Convert reference h5ad to CSVs for RCTD.

    Args:
        adata: Reference AnnData
        output_dir: Output directory
        cell_type_col: Column with cell type annotations
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info(f"Converting reference data: {adata.shape}")

    # Get raw counts if available
    if adata.raw is not None:
        counts_adata = adata.raw.to_adata()
        if hasattr(counts_adata.X, 'toarray'):
            counts = counts_adata.X.toarray()
        else:
            counts = counts_adata.X
        gene_names = counts_adata.var_names
    else:
        if hasattr(adata.X, 'toarray'):
            counts = adata.X.toarray()
        else:
            counts = adata.X
        gene_names = adata.var_names

    # Transpose to genes x cells
    counts_df = pd.DataFrame(
        counts.T,
        index=gene_names,
        columns=adata.obs_names,
    )

    # Save counts
    counts_path = output_dir / "reference_counts.csv"
    counts_df.to_csv(counts_path)
    logger.info(f"  Saved counts: {counts_df.shape} to {counts_path}")

    # Save cell types
    if cell_type_col in adata.obs.columns:
        cell_types = adata.obs[[cell_type_col]].copy()
        cell_types.columns = ['cell_type']
    else:
        raise ValueError(f"Cell type column '{cell_type_col}' not found in adata.obs")

    cell_types_path = output_dir / "reference_cell_types.csv"
    cell_types.to_csv(cell_types_path)
    logger.info(f"  Saved cell types: {cell_types.shape} to {cell_types_path}")

    return counts_path, cell_types_path


def main():
    parser = argparse.ArgumentParser(
        description="Convert h5ad to CSV for RCTD"
    )
    parser.add_argument(
        "--input",
        type=str,
        required=True,
        help="Input h5ad file path",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        required=True,
        help="Output directory for CSVs",
    )
    parser.add_argument(
        "--mode",
        type=str,
        choices=["spatial", "reference"],
        required=True,
        help="Conversion mode: spatial or reference",
    )
    parser.add_argument(
        "--prefix",
        type=str,
        default="",
        help="Filename prefix (for spatial mode)",
    )
    parser.add_argument(
        "--cell-type-col",
        type=str,
        default="cell_type",
        help="Column with cell type annotations (for reference mode)",
    )

    args = parser.parse_args()

    input_path = Path(args.input)
    output_dir = Path(args.output_dir)

    logger.info(f"Loading {input_path}...")
    adata = sc.read_h5ad(input_path)

    if args.mode == "spatial":
        prefix = args.prefix if args.prefix else input_path.stem
        convert_spatial(adata, output_dir, prefix)
    else:
        convert_reference(adata, output_dir, args.cell_type_col)

    logger.info("Conversion complete!")


if __name__ == "__main__":
    main()
