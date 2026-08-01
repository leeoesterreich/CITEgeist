"""
Train Cell2Location RegressionModel on the ER+ Wu 2021 reference.

Loads the normalized h5ad from Task 1, reverses log1p normalization to
recover raw counts (C2L expects raw counts), applies gene filtering,
trains the RegressionModel, and exports the inf_aver signatures.

Output:
    output/competitor_gex/cell2location/reference_model/model/  — saved C2L model
    output/competitor_gex/cell2location/reference_model/inf_aver.csv  — genes x cell_types
    output/competitor_gex/cell2location/reference_model/adata_ref_trained.h5ad
"""

import logging
import time
from pathlib import Path

import cell2location
import numpy as np
import pandas as pd
import scanpy as sc
from cell2location.models import RegressionModel
from cell2location.utils.filtering import filter_genes
from constants import OUTPUT_ROOT
from scipy import sparse

logger = logging.getLogger(__name__)
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
)


def main():
    ref_path = OUTPUT_ROOT / "reference" / "wu_erpos_reference.h5ad"
    out_dir = OUTPUT_ROOT / "cell2location" / "reference_model"
    model_dir = out_dir / "model"
    out_dir.mkdir(parents=True, exist_ok=True)
    model_dir.mkdir(parents=True, exist_ok=True)

    # ----------------------------------------------------------------
    # Load reference
    # ----------------------------------------------------------------
    logger.info(f"Loading ER+ reference from {ref_path}")
    adata_ref = sc.read_h5ad(ref_path)
    logger.info(f"Reference shape: {adata_ref.shape}")
    logger.info(f"Cell types: {adata_ref.obs['cell_type'].value_counts().to_dict()}")

    # ----------------------------------------------------------------
    # Use true raw counts for C2L (stored in .raw by prepare_erpos_reference)
    # ----------------------------------------------------------------
    if adata_ref.raw is not None:
        logger.info("Using .raw.X (true raw counts) for C2L training")
        adata_ref = adata_ref.raw.to_adata()
        # Re-attach cell_type from the processed object
        # (raw.to_adata() preserves obs)
    else:
        logger.warning("No .raw found — falling back to expm1 approximation")
        if sparse.issparse(adata_ref.X):
            X_dense = adata_ref.X.toarray()
        else:
            X_dense = np.array(adata_ref.X)
        X_raw = np.round(np.expm1(X_dense)).astype(np.float32)
        adata_ref.X = sparse.csr_matrix(X_raw)
        del X_dense, X_raw

    logger.info(f"  X range: [{adata_ref.X.min():.0f}, {adata_ref.X.max():.0f}]")

    # ----------------------------------------------------------------
    # Gene filtering (C2L recommended defaults)
    # ----------------------------------------------------------------
    logger.info("Filtering genes ...")
    selected = filter_genes(
        adata_ref,
        cell_count_cutoff=5,
        cell_percentage_cutoff2=0.03,
        nonz_mean_cutoff=1.12,
    )

    # filter_genes returns a boolean Series of selected genes
    adata_ref = adata_ref[:, selected].copy()
    logger.info(f"After gene filtering: {adata_ref.shape}")

    # ----------------------------------------------------------------
    # Add batch key (single-batch reference)
    # ----------------------------------------------------------------
    adata_ref.obs["sample"] = "batch_0"

    # ----------------------------------------------------------------
    # Setup + train RegressionModel
    # ----------------------------------------------------------------
    logger.info("Setting up RegressionModel ...")
    RegressionModel.setup_anndata(
        adata=adata_ref,
        batch_key="sample",
        labels_key="cell_type",
    )

    mod = RegressionModel(adata_ref)

    logger.info("Training RegressionModel (250 epochs, GPU) ...")
    t0 = time.time()
    mod.train(max_epochs=250, accelerator="gpu", batch_size=2500)
    train_time = time.time() - t0
    logger.info(f"Training completed in {train_time:.0f}s")

    # ----------------------------------------------------------------
    # Export posterior
    # ----------------------------------------------------------------
    logger.info("Exporting posterior (1000 samples) ...")
    adata_ref = mod.export_posterior(
        adata_ref,
        sample_kwargs={"num_samples": 1000, "batch_size": 2500, "accelerator": "gpu"},
    )

    # ----------------------------------------------------------------
    # Extract inf_aver (genes x cell_types)
    # ----------------------------------------------------------------
    logger.info("Extracting inf_aver signatures ...")
    if "means_per_cluster_mu_fg" in adata_ref.varm:
        inf_aver = adata_ref.varm["means_per_cluster_mu_fg"].copy()
    elif "q05_per_cluster_mu_fg" in adata_ref.varm:
        # Fallback: some C2L versions use q05
        inf_aver = adata_ref.varm["q05_per_cluster_mu_fg"].copy()
        logger.warning("Using q05_per_cluster_mu_fg instead of means_per_cluster_mu_fg")
    else:
        # Last resort: check .var columns
        factor_names = adata_ref.uns["mod"]["factor_names"]
        sig_cols = [f"means_per_cluster_mu_fg_{ct}" for ct in factor_names]
        inf_aver = adata_ref.var[sig_cols].copy()
        inf_aver.columns = factor_names
        logger.warning("Extracted inf_aver from .var columns")

    logger.info(f"inf_aver shape: {inf_aver.shape}")
    logger.info(f"inf_aver columns: {list(inf_aver.columns)}")

    # ----------------------------------------------------------------
    # Save outputs
    # ----------------------------------------------------------------
    # Save model
    mod.save(str(model_dir), overwrite=True)
    logger.info(f"Saved model to {model_dir}")

    # Save trained reference h5ad
    h5ad_path = out_dir / "adata_ref_trained.h5ad"
    adata_ref.write(h5ad_path)
    logger.info(f"Saved trained reference h5ad to {h5ad_path}")

    # Save inf_aver CSV
    inf_aver_path = out_dir / "inf_aver.csv"
    if isinstance(inf_aver, pd.DataFrame):
        inf_aver.to_csv(inf_aver_path)
    else:
        pd.DataFrame(inf_aver, index=adata_ref.var_names).to_csv(inf_aver_path)
    logger.info(f"Saved inf_aver to {inf_aver_path}")

    # Summary
    logger.info("=" * 60)
    logger.info(f"C2L reference training complete")
    logger.info(f"  Training time: {train_time:.0f}s")
    logger.info(f"  Model: {model_dir}")
    logger.info(f"  inf_aver: {inf_aver_path} ({inf_aver.shape})")
    logger.info("=" * 60)


if __name__ == "__main__":
    main()
