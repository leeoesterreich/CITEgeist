"""Convert h5ad reference to CSV for R methods (RCTD, Seurat).

Reads h5ad, writes genes×cells counts CSV. Only writes if CSV doesn't
already exist (idempotent). Runs in ~30s for 100k cells since it uses
sparse → dense conversion only once.

Usage:
    python convert_h5ad_for_r.py <ref_dir>
    # Reads <ref_dir>/reference.h5ad → writes <ref_dir>/counts.csv
"""
import sys
from pathlib import Path

import pandas as pd
import scanpy as sc


def main():
    ref_dir = Path(sys.argv[1])
    h5ad_path = ref_dir / "reference.h5ad"
    csv_path = ref_dir / "counts.csv"

    if csv_path.exists():
        print(f"counts.csv already exists at {csv_path}, skipping")
        return

    print(f"Converting {h5ad_path} -> {csv_path}")
    adata = sc.read_h5ad(str(h5ad_path))
    X = adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X
    counts = pd.DataFrame(X, index=adata.obs_names, columns=adata.var_names).T
    counts.to_csv(csv_path)
    print(f"Saved {csv_path} ({counts.shape}, {csv_path.stat().st_size / 1e6:.0f} MB)")


if __name__ == "__main__":
    main()
