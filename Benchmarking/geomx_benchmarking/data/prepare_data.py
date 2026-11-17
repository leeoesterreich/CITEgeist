"""
Convert pseudo-GeoMx h5ad to CSV formats for R-based methods (RCTD, CARD, Seurat).

Writes:
  - aoi_counts.csv: genes × AOIs (R convention: rows=genes, cols=samples)
  - aoi_coords.csv: AOI spatial coordinates (aoi_id, x, y)
"""

import logging
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger(__name__)

DATA_DIR = Path(__file__).resolve().parent


def main():
    gex_path = DATA_DIR / "aoi_gex.h5ad"
    logger.info(f"Loading {gex_path}")
    adata = sc.read_h5ad(gex_path)
    logger.info(f"Shape: {adata.shape} (AOIs × genes)")

    X = adata.X.toarray() if hasattr(adata.X, "toarray") else np.asarray(adata.X)
    counts_df = pd.DataFrame(X.T, index=adata.var_names, columns=adata.obs_names)
    counts_path = DATA_DIR / "aoi_counts.csv"
    counts_df.to_csv(counts_path)
    logger.info(f"Wrote {counts_path}: {counts_df.shape[0]} genes × {counts_df.shape[1]} AOIs")

    meta = pd.read_csv(DATA_DIR / "aoi_metadata.csv", index_col=0)
    coords_df = meta[["x", "y"]].copy()
    coords_path = DATA_DIR / "aoi_coords.csv"
    coords_df.to_csv(coords_path)
    logger.info(f"Wrote {coords_path}: {len(coords_df)} AOIs")


if __name__ == "__main__":
    main()
