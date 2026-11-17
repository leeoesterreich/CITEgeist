"""
Prepare Wu et al. 2021 breast cancer scRNA-seq reference for benchmarking.

Loads MTX + metadata, remaps celltype_major labels to breast benchmark names,
exports in per-method formats under reference_data/Wu_2021/.

Label remapping happens HERE (once), so all downstream methods see
breast-standardized names (e.g., "Macrophages" not "Myeloid").
"""

import argparse
import logging
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from scipy.io import mmread

from breast_constants import WU_2021_DIR, WU_TO_BREAST_MAPPING

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")


def load_wu_2021() -> sc.AnnData:
    """Load Wu 2021 from MTX format."""
    logger.info(f"Loading Wu 2021 from {WU_2021_DIR}")

    matrix = mmread(str(WU_2021_DIR / "matrix.mtx"))
    barcodes = pd.read_csv(WU_2021_DIR / "barcodes.tsv", header=None, sep="\t")[0].values
    genes = pd.read_csv(WU_2021_DIR / "genes.tsv", header=None, sep="\t")
    gene_names = genes.iloc[:, 1].values if genes.shape[1] >= 2 else genes.iloc[:, 0].values

    adata = sc.AnnData(X=sparse.csr_matrix(matrix.T))
    adata.obs_names = pd.Index(barcodes)
    adata.var_names = pd.Index(gene_names)
    adata.var_names_make_unique()

    metadata = pd.read_csv(WU_2021_DIR / "metadata.csv", index_col=0)
    common = adata.obs_names.intersection(metadata.index)
    adata = adata[common].copy()
    adata.obs["celltype_major"] = metadata.loc[common, "celltype_major"].values

    logger.info(f"Loaded {adata.n_obs} cells x {adata.n_vars} genes")
    logger.info(f"Cell type distribution:\n{adata.obs['celltype_major'].value_counts().to_string()}")

    return adata


def remap_labels(adata: sc.AnnData) -> sc.AnnData:
    """Remap Wu celltype_major to breast benchmark names."""
    adata.obs["cell_type"] = adata.obs["celltype_major"].map(WU_TO_BREAST_MAPPING)

    unmapped = adata.obs["cell_type"].isna().sum()
    if unmapped > 0:
        logger.warning(f"{unmapped} cells with unmapped cell types — dropping")
        adata = adata[adata.obs["cell_type"].notna()].copy()

    logger.info(f"Remapped cell type distribution:\n{adata.obs['cell_type'].value_counts().to_string()}")
    return adata


def export_csv_format(adata: sc.AnnData, output_dir: Path, subdir: str, chunk_size: int = 500) -> None:
    """Export counts + cell_types CSVs for R-based methods.

    Writes the genes×cells counts matrix in chunks directly from the sparse
    matrix, converting each chunk to integer strings via vectorised numpy.
    ~50× faster than pd.DataFrame.to_csv() for a 30K×100K matrix.
    """
    import io

    out = output_dir / subdir
    out.mkdir(parents=True, exist_ok=True)

    n_genes = adata.n_vars
    csv_path = out / "counts.csv"
    logger.info(f"  Writing counts.csv ({n_genes} genes x {adata.n_obs} cells) ...")

    # genes × cells (R convention: rows=genes, cols=cells)
    X_T = adata.X.T.tocsr() if sparse.issparse(adata.X) else adata.X.T

    with open(csv_path, "w", buffering=8 * 1024 * 1024) as f:
        f.write("," + ",".join(adata.obs_names) + "\n")
        for start in range(0, n_genes, chunk_size):
            end = min(start + chunk_size, n_genes)
            chunk = X_T[start:end]
            if sparse.issparse(chunk):
                chunk = chunk.toarray()
            # scRNA-seq counts are integers — use %d via numpy for vectorised formatting
            chunk_int = chunk.astype(np.int32)
            buf = io.StringIO()
            np.savetxt(buf, chunk_int, fmt="%d", delimiter=",")
            lines = buf.getvalue().split("\n")
            gene_names = adata.var_names[start:end]
            for gene, line in zip(gene_names, lines):
                if line:
                    f.write(f"{gene},{line}\n")
            if end % 5000 < chunk_size or end == n_genes:
                logger.info(f"    {end}/{n_genes} genes written")

    logger.info(f"  Wrote counts.csv")

    ct_df = pd.DataFrame({"cell_type": adata.obs["cell_type"].values}, index=adata.obs_names)
    ct_df.to_csv(out / "cell_types.csv")
    logger.info(f"  Wrote cell_types.csv ({len(ct_df)} cells)")


def export_h5ad_format(adata: sc.AnnData, output_dir: Path, subdir: str) -> None:
    """Export AnnData h5ad for Python-based methods (Tangram, C2L)."""
    out = output_dir / subdir
    out.mkdir(parents=True, exist_ok=True)

    adata_out = adata.copy()
    adata_out.write(out / "reference.h5ad")
    logger.info(f"  Wrote reference.h5ad ({adata_out.n_obs} cells)")


def main():
    parser = argparse.ArgumentParser(description="Prepare Wu 2021 reference for breast benchmark")
    parser.add_argument(
        "--output-dir",
        default="reference_data/Wu_2021",
        help="Output directory relative to xenium_benchmarking/",
    )
    args = parser.parse_args()

    bench_dir = Path(__file__).parent
    output_dir = bench_dir / args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    adata = load_wu_2021()
    adata = remap_labels(adata)

    logger.info("Exporting CSV format (RCTD, CARD, Seurat)...")
    export_csv_format(adata, output_dir, "rctd")
    # Seurat and CARD use identical format — symlink to avoid 3x write of ~60GB CSV
    for method in ["seurat", "card"]:
        method_dir = output_dir / method
        method_dir.mkdir(parents=True, exist_ok=True)
        rctd_dir = output_dir / "rctd"
        for fname in ["counts.csv", "cell_types.csv"]:
            target = method_dir / fname
            if target.exists() or target.is_symlink():
                target.unlink()
            target.symlink_to(rctd_dir / fname)
        logger.info(f"  Symlinked {method}/ -> rctd/")

    logger.info("Exporting h5ad format (Tangram, C2L)...")
    export_h5ad_format(adata, output_dir, "tangram")
    export_h5ad_format(adata, output_dir, "cell2location")

    logger.info(f"All reference formats saved to {output_dir}")


if __name__ == "__main__":
    main()
