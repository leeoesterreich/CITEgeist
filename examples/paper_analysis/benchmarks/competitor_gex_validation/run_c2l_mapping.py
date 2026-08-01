"""
Cell2Location spatial mapping for a single patient Visium sample.

Loads pre-trained inf_aver signatures from the reference model, maps a
Visium sample via Cell2Location, and exports:
  - Normalized proportions CSV
  - Per-cell-type GEX layers as an h5ad with spatial coords

Usage:
    python -u run_c2l_mapping.py --sample HCC22-088-P1-S1
"""

import argparse
import logging
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

import cell2location

from constants import OUTPUT_ROOT, SPACERANGER_ROOT, SAMPLES

logger = logging.getLogger(__name__)
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
)


def load_visium_with_spatial(sample: str) -> sc.AnnData:
    """Load Visium h5 matrix and inject spatial coordinates from tissue_positions.csv.

    sq.read.visium(load_images=False) does NOT populate obsm["spatial"],
    so we load the matrix directly and inject coords manually.
    """
    h5_path = SPACERANGER_ROOT / sample / "outs" / "filtered_feature_bc_matrix.h5"
    spatial_path = SPACERANGER_ROOT / sample / "outs" / "spatial" / "tissue_positions.csv"

    logger.info(f"Loading h5 matrix from {h5_path}")
    adata = sc.read_10x_h5(str(h5_path))
    adata.var_names_make_unique()
    logger.info(f"  Loaded: {adata.shape}")

    # Inject spatial coordinates
    logger.info(f"Loading spatial coords from {spatial_path}")
    pos = pd.read_csv(spatial_path, index_col=0)

    # Subset to barcodes present in the filtered matrix
    common = adata.obs_names.intersection(pos.index)
    if len(common) < adata.n_obs:
        logger.warning(
            f"  {adata.n_obs - len(common)} barcodes missing from tissue_positions; "
            f"subsetting to {len(common)} common barcodes"
        )
        adata = adata[common].copy()

    pos = pos.loc[adata.obs_names]
    adata.obsm["spatial"] = pos[["pxl_col_in_fullres", "pxl_row_in_fullres"]].values.astype(np.float64)
    adata.obs["in_tissue"] = pos["in_tissue"].values
    adata.obs["array_row"] = pos["array_row"].values
    adata.obs["array_col"] = pos["array_col"].values

    # Filter to in-tissue spots only
    in_tissue_mask = adata.obs["in_tissue"] == 1
    if in_tissue_mask.sum() < adata.n_obs:
        logger.info(f"  Filtering to {in_tissue_mask.sum()} in-tissue spots (from {adata.n_obs})")
        adata = adata[in_tissue_mask].copy()

    return adata


def main():
    parser = argparse.ArgumentParser(description="C2L spatial mapping for one patient sample")
    parser.add_argument("--sample", type=str, required=True, choices=SAMPLES, help="Sample ID (e.g. HCC22-088-P1-S1)")
    args = parser.parse_args()

    sample = args.sample
    out_dir = OUTPUT_ROOT / "cell2location" / sample
    out_dir.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------
    # Load pre-trained inf_aver
    # ------------------------------------------------------------------
    inf_aver_path = OUTPUT_ROOT / "cell2location" / "reference_model" / "inf_aver.csv"
    logger.info(f"Loading inf_aver from {inf_aver_path}")
    inf_aver = pd.read_csv(inf_aver_path, index_col=0)
    logger.info(f"  inf_aver shape: {inf_aver.shape} (genes x cell_types)")
    logger.info(f"  Cell types: {list(inf_aver.columns)}")

    # ------------------------------------------------------------------
    # Load Visium sample
    # ------------------------------------------------------------------
    adata_vis = load_visium_with_spatial(sample)
    logger.info(f"Visium sample {sample}: {adata_vis.shape}")

    # Ensure raw integer counts (C2L expects raw counts)
    if sparse.issparse(adata_vis.X):
        x_max = adata_vis.X.max()
    else:
        x_max = np.max(adata_vis.X)
    logger.info(f"  X max value: {x_max:.1f}")

    # ------------------------------------------------------------------
    # Remove MT genes
    # ------------------------------------------------------------------
    adata_vis.var["SYMBOL"] = adata_vis.var_names
    adata_vis.var["MT_gene"] = [g.startswith("MT-") for g in adata_vis.var["SYMBOL"]]
    mt_count = adata_vis.var["MT_gene"].sum()
    logger.info(f"  Removing {mt_count} MT genes")
    adata_vis.obsm["MT"] = adata_vis[:, adata_vis.var["MT_gene"].values].X.toarray()
    adata_vis = adata_vis[:, ~adata_vis.var["MT_gene"].values].copy()

    # ------------------------------------------------------------------
    # Intersect genes between spatial and reference
    # ------------------------------------------------------------------
    intersect = np.intersect1d(adata_vis.var_names, inf_aver.index)
    logger.info(f"  Shared genes: {len(intersect)} " f"(spatial: {adata_vis.n_vars}, ref: {inf_aver.shape[0]})")
    adata_vis = adata_vis[:, intersect].copy()
    inf_aver_shared = inf_aver.loc[intersect, :].copy()

    # ------------------------------------------------------------------
    # Setup and train C2L spatial model
    # ------------------------------------------------------------------
    logger.info("Setting up Cell2location spatial model ...")
    cell2location.models.Cell2location.setup_anndata(adata=adata_vis)

    mod = cell2location.models.Cell2location(
        adata_vis,
        cell_state_df=inf_aver_shared,
        N_cells_per_location=5,
        detection_alpha=200,
    )

    logger.info("Training Cell2location (30000 epochs, GPU) ...")
    t0 = time.time()
    mod.train(max_epochs=30000, batch_size=None, train_size=1, accelerator="gpu")
    train_time = time.time() - t0
    logger.info(f"Training completed in {train_time:.0f}s")

    # ------------------------------------------------------------------
    # Export posterior
    # ------------------------------------------------------------------
    logger.info("Exporting posterior (3000 samples) ...")
    adata_vis = mod.export_posterior(
        adata_vis,
        sample_kwargs={
            "num_samples": 3000,
            "batch_size": mod.adata.n_obs,
            "accelerator": "gpu",
        },
    )

    # ------------------------------------------------------------------
    # Extract and save proportions
    # ------------------------------------------------------------------
    logger.info("Extracting proportions from q05_cell_abundance_w_sf ...")
    # NOTE: C2L .copy() on varm DataFrames returns zeros — use .values
    df_abund = adata_vis.obsm["q05_cell_abundance_w_sf"]
    if isinstance(df_abund, pd.DataFrame):
        abund_values = df_abund.values
        col_names = list(df_abund.columns)
    else:
        abund_values = np.array(df_abund)
        col_names = list(adata_vis.uns["mod"]["factor_names"])

    row_sums = abund_values.sum(axis=1, keepdims=True)
    row_sums = np.maximum(row_sums, 1e-10)  # avoid division by zero
    prop_values = abund_values / row_sums

    proportions = pd.DataFrame(
        prop_values,
        index=adata_vis.obs_names,
        columns=col_names,
    )
    proportions.index.name = "spot"

    prop_path = out_dir / "proportions.csv"
    proportions.to_csv(prop_path)
    logger.info(f"Saved proportions to {prop_path} ({proportions.shape})")

    # ------------------------------------------------------------------
    # Extract per-cell-type GEX layers
    # ------------------------------------------------------------------
    logger.info("Computing expected expression per cell type ...")
    expected_dict = mod.module.model.compute_expected_per_cell_type(mod.samples["post_sample_q05"], mod.adata_manager)

    # Build output h5ad with layers
    adata_out = sc.AnnData(
        X=adata_vis.X.copy() if sparse.issparse(adata_vis.X) else adata_vis.X.copy(),
        obs=adata_vis.obs[["array_row", "array_col"]].copy(),
        var=adata_vis.var[["SYMBOL"]].copy(),
    )
    adata_out.obsm["spatial"] = adata_vis.obsm["spatial"].copy()

    for i, ct_name in enumerate(mod.factor_names_):
        layer_data = expected_dict["mu"][i]
        if sparse.issparse(layer_data):
            layer_data = layer_data.toarray()
        adata_out.layers[ct_name] = np.array(layer_data, dtype=np.float32)
        logger.info(f"  Layer '{ct_name}': shape {layer_data.shape}, " f"mean={np.mean(layer_data):.4f}")

    gex_path = out_dir / "c2l_gex_layers.h5ad"
    adata_out.write(gex_path)
    logger.info(f"Saved GEX layers h5ad to {gex_path}")

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------
    logger.info("=" * 60)
    logger.info(f"C2L mapping complete for {sample}")
    logger.info(f"  Training time: {train_time:.0f}s")
    logger.info(f"  Spots: {adata_vis.n_obs}, Genes: {adata_vis.n_vars}")
    logger.info(f"  Cell types: {mod.factor_names_}")
    logger.info(f"  Proportions: {prop_path}")
    logger.info(f"  GEX layers: {gex_path}")
    logger.info("=" * 60)


if __name__ == "__main__":
    main()
