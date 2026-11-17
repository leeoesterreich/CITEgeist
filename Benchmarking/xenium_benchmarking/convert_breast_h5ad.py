"""
Convert breast pseudo-Visium h5ad files to CSV format for R methods.

Produces per-region:
  - {prefix}_counts.csv (genes x spots)
  - {prefix}_coords.csv (x, y per spot)
"""

import argparse
import logging
from pathlib import Path

import pandas as pd
import scanpy as sc
from scipy import sparse

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")


def convert_region(h5ad_path: Path, output_dir: Path, prefix: str = "Xenium") -> None:
    """Convert one region's h5ad to CSV format."""
    adata = sc.read_h5ad(h5ad_path)
    logger.info(f"Loaded {h5ad_path.name}: {adata.n_obs} spots x {adata.n_vars} genes")

    output_dir.mkdir(parents=True, exist_ok=True)

    X = adata.X.toarray() if sparse.issparse(adata.X) else adata.X
    counts_df = pd.DataFrame(X.T, index=adata.var_names, columns=adata.obs_names)
    counts_df.to_csv(output_dir / f"{prefix}_counts.csv")

    coords_df = pd.DataFrame(adata.obsm["spatial"], index=adata.obs_names, columns=["x", "y"])
    coords_df.to_csv(output_dir / f"{prefix}_coords.csv")

    logger.info(f"  Saved counts ({counts_df.shape}) and coords ({coords_df.shape})")


def main():
    parser = argparse.ArgumentParser(description="Convert breast h5ad to CSV for R methods")
    parser.add_argument("--input-dir", default="data_breast/h5ad_objects")
    parser.add_argument("--output-dir", default="data_breast/csv_for_r")
    parser.add_argument("--n-regions", type=int, default=5)
    parser.add_argument("--prefix", default="Xenium")
    args = parser.parse_args()

    bench_dir = Path(__file__).parent
    input_dir = bench_dir / args.input_dir
    output_dir = bench_dir / args.output_dir

    for i in range(args.n_regions):
        h5ad_path = input_dir / f"Xenium_region_{i}_GEX.h5ad"
        region_out = output_dir / f"Xenium_region_{i}"
        convert_region(h5ad_path, region_out, args.prefix)


if __name__ == "__main__":
    main()
