"""
Prepare ER+ filtered Wu et al. 2021 scRNA-seq reference for competitor GEX validation.

Loads the full Wu 2021 dataset (MTX format), filters to ER+ subtype patients,
remaps cell types to our breast benchmark names, normalizes, and saves as h5ad.

Output: output/competitor_gex/reference/wu_erpos_reference.h5ad
"""

import logging
from pathlib import Path

import pandas as pd
import scanpy as sc
from scipy import sparse
from scipy.io import mmread

from constants import ERPOS_CIDS, OUTPUT_ROOT, WU_DATA_DIR, WU_TO_BREAST_MAPPING

logger = logging.getLogger(__name__)
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
)


def load_wu_2021() -> sc.AnnData:
    """Load Wu 2021 from MTX format with metadata."""
    logger.info(f"Loading Wu 2021 from {WU_DATA_DIR}")

    matrix = mmread(str(WU_DATA_DIR / "matrix.mtx"))
    barcodes = pd.read_csv(WU_DATA_DIR / "barcodes.tsv", header=None, sep="\t")[0].values
    genes = pd.read_csv(WU_DATA_DIR / "genes.tsv", header=None, sep="\t")
    gene_names = genes.iloc[:, 1].values if genes.shape[1] >= 2 else genes.iloc[:, 0].values

    adata = sc.AnnData(X=sparse.csr_matrix(matrix.T))
    adata.obs_names = pd.Index(barcodes)
    adata.var_names = pd.Index(gene_names)
    adata.var_names_make_unique()

    # Join metadata
    metadata = pd.read_csv(WU_DATA_DIR / "metadata.csv", index_col=0)
    common = adata.obs_names.intersection(metadata.index)
    adata = adata[common].copy()
    adata.obs["celltype_major"] = metadata.loc[common, "celltype_major"].values
    adata.obs["subtype"] = metadata.loc[common, "subtype"].values

    # Extract patient ID from barcode (CID####_BARCODE format)
    adata.obs["patient_id"] = [idx.rsplit("_", 1)[0] for idx in adata.obs_names]

    logger.info(f"Loaded {adata.n_obs} cells x {adata.n_vars} genes")
    logger.info(f"Subtypes: {adata.obs['subtype'].value_counts().to_dict()}")
    logger.info(f"Patients: {adata.obs['patient_id'].nunique()}")

    return adata


def filter_erpos(adata: sc.AnnData) -> sc.AnnData:
    """Filter to ER+ subtype patients."""
    logger.info("Filtering to ER+ subtype ...")

    mask = adata.obs["subtype"] == "ER+"
    adata_erpos = adata[mask].copy()

    # Verify patient IDs match expected
    found_cids = sorted(adata_erpos.obs["patient_id"].unique())
    expected_cids = sorted(ERPOS_CIDS)
    logger.info(f"Expected ER+ patients: {expected_cids}")
    logger.info(f"Found ER+ patients:    {found_cids}")

    missing = set(expected_cids) - set(found_cids)
    extra = set(found_cids) - set(expected_cids)
    if missing:
        logger.warning(f"Missing expected patients: {missing}")
    if extra:
        logger.warning(f"Unexpected extra patients: {extra}")

    logger.info(f"ER+ subset: {adata_erpos.n_obs} cells from {len(found_cids)} patients")
    return adata_erpos


def remap_cell_types(adata: sc.AnnData) -> sc.AnnData:
    """Remap Wu celltype_major to breast benchmark names."""
    adata.obs["cell_type"] = adata.obs["celltype_major"].map(WU_TO_BREAST_MAPPING)

    unmapped = adata.obs["cell_type"].isna().sum()
    if unmapped > 0:
        unmapped_types = adata.obs.loc[adata.obs["cell_type"].isna(), "celltype_major"].value_counts()
        logger.warning(f"{unmapped} cells with unmapped cell types — dropping")
        logger.warning(f"Unmapped types:\n{unmapped_types.to_string()}")
        adata = adata[adata.obs["cell_type"].notna()].copy()

    logger.info(f"Remapped cell type distribution:\n{adata.obs['cell_type'].value_counts().to_string()}")
    return adata


def main():
    ref_dir = OUTPUT_ROOT / "reference"
    ref_dir.mkdir(parents=True, exist_ok=True)
    out_path = ref_dir / "wu_erpos_reference.h5ad"

    # Load full dataset
    adata = load_wu_2021()

    # Filter to ER+
    adata = filter_erpos(adata)

    # Remap cell types
    adata = remap_cell_types(adata)

    # Normalize (keep raw counts in .raw for methods that need them)
    logger.info("Normalizing: total-count + log1p ...")
    adata.raw = adata.copy()
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)

    # Summary stats
    logger.info("=" * 60)
    logger.info(f"Final reference: {adata.n_obs} cells x {adata.n_vars} genes")
    logger.info(f"Patients: {sorted(adata.obs['patient_id'].unique())}")
    logger.info(f"Cell types ({adata.obs['cell_type'].nunique()}):")
    for ct, count in adata.obs["cell_type"].value_counts().items():
        logger.info(f"  {ct}: {count}")
    logger.info("=" * 60)

    # Save
    adata.write(out_path)
    logger.info(f"Saved ER+ reference to {out_path}")


if __name__ == "__main__":
    main()
