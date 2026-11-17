#!/usr/bin/env python
"""Export reference.h5ad to a 10x MTX triplet for Read10X().

The dense X matrix (54850 cells × 22169 genes) is read in chunks to avoid
OOM. Output is genes × cells (transposed) as required by Seurat's Read10X.

Barcodes (ad.obs_names) are written as-is; they match cell_types.csv 1:1
(verified: 54850/54850 intersection, 0 mismatches).
"""
import sys as _sys
from pathlib import Path as _Path

_sys.path.insert(0, str(_Path(__file__).resolve().parents[3]))
import _paths as _cg  # noqa: E402

import h5py
import scipy.io
import scipy.sparse
import gzip
import numpy as np
from pathlib import Path

BASE = Path(
    str(_cg.REPO_ROOT / "Benchmarking" / "xenium_benchmarking" / "reference_data" / "GSE156632" / "processed" / "seurat")
)
H5AD = BASE / "reference.h5ad"
OUT = BASE / "mtx"
OUT.mkdir(parents=True, exist_ok=True)

CHUNK = 2000  # cells per chunk (safe for 64 G RAM)


def load_names(f, group_key):
    names = list(f[group_key]["_index"][:])
    if isinstance(names[0], bytes):
        names = [x.decode() for x in names]
    return names


def main():
    print(f"Opening {H5AD} ...")
    with h5py.File(H5AD, "r") as f:
        obs_names = load_names(f, "obs")
        var_names = load_names(f, "var")
        n_cells = len(obs_names)
        n_genes = len(var_names)
        print(f"Shape: {n_cells} cells × {n_genes} genes")

        # Build CSC genes×cells by reading the dense X in row (cell) chunks
        # and stacking into a COO then converting to CSC at the end.
        print("Building sparse matrix in chunks ...")
        row_chunks = []
        for start in range(0, n_cells, CHUNK):
            end = min(start + CHUNK, n_cells)
            block = f["X"][start:end, :]  # (chunk_cells, n_genes) dense
            row_chunks.append(scipy.sparse.csr_matrix(block.astype(np.float32)))
            if (start // CHUNK) % 10 == 0:
                print(f"  chunk {start}/{n_cells} ...")

        print("Stacking and transposing to genes×cells ...")
        cells_x_genes = scipy.sparse.vstack(row_chunks, format="csr")
        genes_x_cells = cells_x_genes.T.tocsc()
        print(f"Final matrix shape: {genes_x_cells.shape}  nnz={genes_x_cells.nnz}")

    # Write matrix.mtx.gz (scipy mmwrite writes in row-major COO order)
    mtx_path = OUT / "matrix.mtx.gz"
    print(f"Writing {mtx_path} ...")
    with gzip.open(mtx_path, "wb") as fh:
        scipy.io.mmwrite(fh, genes_x_cells)

    # Write barcodes.tsv.gz  (column names = cell barcodes)
    bc_path = OUT / "barcodes.tsv.gz"
    print(f"Writing {bc_path} ...")
    with gzip.open(bc_path, "wt") as fh:
        fh.write("\n".join(obs_names) + "\n")

    # Write features.tsv.gz  (row names = gene names)
    ft_path = OUT / "features.tsv.gz"
    print(f"Writing {ft_path} ...")
    with gzip.open(ft_path, "wt") as fh:
        for g in var_names:
            fh.write(f"{g}\t{g}\tGene Expression\n")

    print("Done.")
    print(f"Output: {OUT}")
    print(f"  matrix.mtx.gz   → {genes_x_cells.shape[0]} genes × {genes_x_cells.shape[1]} cells")
    print(f"  barcodes.tsv.gz → {len(obs_names)} barcodes")
    print(f"  features.tsv.gz → {len(var_names)} features")


if __name__ == "__main__":
    main()
