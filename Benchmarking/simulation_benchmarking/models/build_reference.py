#!/usr/bin/env python
"""Build ER+ scRNA-seq reference from raw Wu et al. 2021 data.

Loads the full Wu_2021 dataset from the lab shared data folder,
filters to ER+ subtype (11 patients, ~38k cells), selects top 2000
highly variable genes, and writes the canonical reference h5ad.

Usage:
    python build_reference.py
    python build_reference.py --output models/Wu_scRNA_ref_ERpos.h5ad
"""

import sys as _sys
from pathlib import Path as _Path

_sys.path.insert(0, str(_Path(__file__).resolve().parents[2]))
import _paths as _cg  # noqa: E402

import argparse
import logging
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.io import mmread
from scipy.sparse import csc_matrix

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger(__name__)

WU_RAW = Path(str(_cg.root("DATA") / "BC-Datasets" / "Wu_2021" / "Wu_etal_2021_BRCA_scRNASeq"))
SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_OUTPUT = SCRIPT_DIR / "Wu_scRNA_ref_ERpos.h5ad"

EXPECTED_CELL_TYPES = sorted([
    "B-cells", "CAFs", "Cancer Epithelial", "Endothelial",
    "Myeloid", "Normal Epithelial", "PVL", "Plasmablasts", "T-cells",
])
N_HVG = 2000


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", type=str, default=str(DEFAULT_OUTPUT))
    parser.add_argument("--wu-path", type=str, default=str(WU_RAW))
    args = parser.parse_args()

    wu_path = Path(args.wu_path)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)

    logger.info("Loading raw Wu_2021 data from %s", wu_path)
    genes = pd.read_csv(wu_path / "genes.tsv", header=None, sep="\t")[0].values
    barcodes = pd.read_csv(wu_path / "barcodes.tsv", header=None, sep="\t")[0].values
    X = csc_matrix(mmread(wu_path / "matrix.mtx")).T  # genes×cells → cells×genes
    adata = sc.AnnData(X=X, obs=pd.DataFrame(index=barcodes), var=pd.DataFrame(index=genes))
    logger.info("Raw shape: %s", adata.shape)

    meta = pd.read_csv(wu_path / "metadata.csv", index_col=0)
    logger.info("Metadata: %d cells, columns: %s", len(meta), list(meta.columns))

    common = adata.obs_names.intersection(meta.index)
    logger.info("Common barcodes: %d / %d", len(common), adata.n_obs)
    adata = adata[common].copy()
    meta = meta.loc[common]

    for col in ["orig.ident", "subtype", "celltype_major", "celltype_minor"]:
        adata.obs[col] = meta[col].values
    adata.obs["patient_ID"] = meta["orig.ident"].values

    er_mask = adata.obs["subtype"] == "ER+"
    logger.info("ER+ cells: %d / %d", er_mask.sum(), adata.n_obs)
    adata = adata[er_mask].copy()

    actual_types = sorted(adata.obs["celltype_major"].unique().tolist())
    assert actual_types == EXPECTED_CELL_TYPES, (
        f"Cell type mismatch: expected {EXPECTED_CELL_TYPES}, got {actual_types}"
    )
    logger.info("Cell types: %s", dict(adata.obs["celltype_major"].value_counts()))

    logger.info("Selecting top %d HVGs...", N_HVG)
    sc.pp.highly_variable_genes(adata, n_top_genes=N_HVG, flavor="seurat_v3", span=0.3)
    adata = adata[:, adata.var["highly_variable"]].copy()
    adata.raw = adata.copy()
    adata.layers["raw"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    logger.info("Final shape: %s", adata.shape)

    assert adata.shape == (38241, N_HVG), f"Unexpected shape: {adata.shape}"
    assert len(adata.obs["patient_ID"].unique()) == 11

    logger.info("Writing to %s", output)
    adata.write_h5ad(output)
    logger.info("Done. File size: %.1f MB", output.stat().st_size / 1e6)


if __name__ == "__main__":
    main()
