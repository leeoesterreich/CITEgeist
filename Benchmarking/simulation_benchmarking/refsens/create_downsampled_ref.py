"""
Create patient-dropout reference for reference-sensitivity benchmarking.

Tests whether deconvolution methods degrade when the scRNA-seq reference
loses patient diversity. The "downsampled" reference drops HALF the patients
(not uniform cell subsampling — that showed no difference).

Strategy:
  - Full reference: all 26 Wu 2021 patients (~100k cells)
  - Dropout reference: 13 patients (half), selected to maximize cell-type
    composition shift vs full (adversarial dropout)
  - Seed: 42 for reproducibility

Usage:
    python -u create_downsampled_ref.py
"""

import sys as _sys
from pathlib import Path as _Path

_sys.path.insert(0, str(_Path(__file__).resolve().parents[2]))
import _paths as _cg  # noqa: E402

import logging
import os
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.io import mmread

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)

WU_MTX_DIR = str(_cg.root("DATA") / "BC-Datasets" / "Wu_2021" / "Wu_etal_2021_BRCA_scRNASeq")
REFSENS_DIR = Path(__file__).resolve().parent
OUTPUT_DIR = REFSENS_DIR / "reference"
SEED = 42


def load_wu_reference(wu_dir: str):
    """Load raw Wu 2021 scRNA-seq from Matrix Market files."""
    logger.info(f"Loading Wu reference from {wu_dir} ...")

    matrix = mmread(os.path.join(wu_dir, "matrix.mtx")).T.tocsr().astype("float32")
    genes_df = pd.read_csv(os.path.join(wu_dir, "genes.tsv"), header=None, sep="\t")
    barcodes_df = pd.read_csv(os.path.join(wu_dir, "barcodes.tsv"), header=None)
    metadata_df = pd.read_csv(os.path.join(wu_dir, "metadata.csv"), index_col=0)

    import anndata
    adata = anndata.AnnData(
        X=matrix,
        obs=pd.DataFrame(index=barcodes_df[0]),
        var=pd.DataFrame(index=genes_df[0]),
    )
    adata.obs.index.name = None
    adata.obs["cell_type"] = metadata_df["celltype_major"].values
    adata.obs["patient"] = metadata_df["orig.ident"].values

    logger.info(f"  Shape: {adata.shape}")
    n_patients = adata.obs["patient"].nunique()
    logger.info(f"  Patients: {n_patients}")
    logger.info(f"  Cell types: {adata.obs['cell_type'].value_counts().to_dict()}")
    return adata


def patient_dropout(adata, seed: int):
    """Drop half the patients to create a composition-shifted reference.

    Selects the dropout set to maximize cell-type composition divergence
    from the full reference (adversarial patient selection).
    """
    rng = np.random.default_rng(seed)
    patients = sorted(adata.obs["patient"].unique())
    n_keep = len(patients) // 2
    n_drop = len(patients) - n_keep

    logger.info(f"\nPatient dropout: keeping {n_keep}/{len(patients)} patients")

    # Compute per-patient cell-type composition
    ct_by_patient = pd.crosstab(adata.obs["patient"], adata.obs["cell_type"],
                                normalize="index")
    full_composition = adata.obs["cell_type"].value_counts(normalize=True).sort_index()

    # Greedy adversarial selection: iteratively drop the patient whose removal
    # maximizes composition shift (L1 distance from full)
    remaining = list(patients)
    dropped = []

    for _ in range(n_drop):
        best_patient = None
        best_shift = -1

        for p in remaining:
            test_remaining = [x for x in remaining if x != p]
            mask = adata.obs["patient"].isin(test_remaining)
            test_comp = adata.obs.loc[mask, "cell_type"].value_counts(
                normalize=True).reindex(full_composition.index, fill_value=0)
            shift = (test_comp - full_composition).abs().sum()
            if shift > best_shift:
                best_shift = shift
                best_patient = p

        remaining.remove(best_patient)
        dropped.append(best_patient)
        logger.info(f"  Dropped {best_patient} (L1 shift: {best_shift:.4f})")

    logger.info(f"\n  Kept patients ({n_keep}): {remaining}")
    logger.info(f"  Dropped patients ({n_drop}): {dropped}")

    mask = adata.obs["patient"].isin(remaining)
    adata_dropout = adata[mask].copy()

    # Report composition shift
    dropout_comp = adata_dropout.obs["cell_type"].value_counts(
        normalize=True).reindex(full_composition.index, fill_value=0)
    logger.info(f"\n  Composition comparison (full vs dropout):")
    for ct in full_composition.index:
        logger.info(f"    {ct:25s}: {full_composition[ct]:.4f} -> {dropout_comp[ct]:.4f} "
                     f"(delta={dropout_comp[ct] - full_composition[ct]:+.4f})")
    total_shift = (dropout_comp - full_composition).abs().sum()
    logger.info(f"  Total L1 composition shift: {total_shift:.4f}")
    logger.info(f"  Dropout reference: {adata_dropout.shape}")

    return adata_dropout, remaining, dropped


def save_reference(adata, output_dir: Path, label: str):
    """Save reference in h5ad + CSV formats."""
    out = output_dir / label
    out.mkdir(parents=True, exist_ok=True)

    adata.var.index.name = None
    adata.var.index = adata.var.index.astype(str)
    adata.obs.index.name = None
    adata.obs.index = adata.obs.index.astype(str)

    h5ad_path = out / "reference.h5ad"
    adata.write(h5ad_path)
    logger.info(f"  Saved h5ad: {h5ad_path} ({adata.shape})")

    # Cell types CSV only — R scripts use convert_h5ad_for_r.py for counts
    ct_df = adata.obs[["cell_type"]].copy()
    ct_path = out / "cell_types.csv"
    ct_df.to_csv(ct_path)
    logger.info(f"  Saved cell_types CSV: {ct_path} ({len(ct_df)} cells)")

    return out


def main():
    logger.info("=" * 60)
    logger.info("Create Patient-Dropout Reference for Refsens Benchmark")
    logger.info("=" * 60)

    adata_full = load_wu_reference(WU_MTX_DIR)

    logger.info("\nSaving full reference ...")
    save_reference(adata_full, OUTPUT_DIR, "full")

    logger.info("\nCreating patient-dropout reference ...")
    adata_dropout, kept, dropped = patient_dropout(adata_full, SEED)
    save_reference(adata_dropout, OUTPUT_DIR, "downsampled")

    # Save patient lists for documentation
    meta_path = OUTPUT_DIR / "dropout_metadata.json"
    import json
    with open(meta_path, "w") as f:
        json.dump({
            "strategy": "adversarial_patient_dropout",
            "seed": SEED,
            "n_patients_full": len(kept) + len(dropped),
            "n_patients_kept": len(kept),
            "n_patients_dropped": len(dropped),
            "kept_patients": kept,
            "dropped_patients": dropped,
            "full_n_cells": int(adata_full.shape[0]),
            "dropout_n_cells": int(adata_dropout.shape[0]),
        }, f, indent=2)
    logger.info(f"\nSaved dropout metadata: {meta_path}")

    logger.info("\n" + "=" * 60)
    logger.info("Summary")
    logger.info("=" * 60)
    logger.info(f"Full reference:    {adata_full.shape} ({len(kept) + len(dropped)} patients)")
    logger.info(f"Dropout reference: {adata_dropout.shape} ({len(kept)} patients)")
    logger.info(f"Output directory:  {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
