"""
Run Tangram spatial mapping for a single Visium sample.

Loads the ER+ Wu 2021 reference (log1p-normalized), loads a Visium sample
from SpaceRanger, finds marker genes, runs Tangram cell-to-space mapping,
extracts per-cell-type GEX decomposition, and saves results.

Output:
    output/competitor_gex/tangram/{sample}/proportions.csv
    output/competitor_gex/tangram/{sample}/gex_layers.h5ad
    output/competitor_gex/tangram/{sample}/mapping_matrix.npz
"""

import argparse
import logging
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import tangram as tg
from constants import CELL_TYPES, OUTPUT_ROOT, SPACERANGER_ROOT, WU_TO_BREAST_MAPPING

logger = logging.getLogger(__name__)
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
)


def load_visium_sample(sample: str) -> sc.AnnData:
    """Load a Visium sample from SpaceRanger output and inject spatial coords."""
    sample_dir = SPACERANGER_ROOT / sample / "outs"
    h5_path = sample_dir / "filtered_feature_bc_matrix.h5"

    logger.info(f"Loading Visium sample from {h5_path}")
    ad_sp = sc.read_10x_h5(str(h5_path))
    ad_sp.var_names_make_unique()

    # Inject spatial coordinates from tissue_positions.csv
    # sq.read.visium(load_images=False) does NOT populate obsm["spatial"]
    pos_path = sample_dir / "spatial" / "tissue_positions.csv"
    pos_df = pd.read_csv(pos_path, index_col=0)

    # Columns: in_tissue, array_row, array_col, pxl_row_in_fullres, pxl_col_in_fullres
    common = ad_sp.obs_names.intersection(pos_df.index)
    pos_df = pos_df.loc[common]
    ad_sp = ad_sp[common].copy()

    # Store spatial coords as (x, y) = (pxl_col, pxl_row)
    ad_sp.obsm["spatial"] = pos_df[["pxl_col_in_fullres", "pxl_row_in_fullres"]].values.astype(np.float64)

    # Filter to in-tissue spots
    in_tissue = pos_df["in_tissue"].values == 1
    ad_sp = ad_sp[in_tissue].copy()
    logger.info(f"Loaded {ad_sp.n_obs} in-tissue spots, {ad_sp.n_vars} genes")

    return ad_sp


def find_marker_genes(ad_sc: sc.AnnData, n_genes: int = 100) -> list:
    """Find marker genes per cell type using Wilcoxon rank-sum test."""
    logger.info(f"Finding top {n_genes} marker genes per cell type ...")
    sc.tl.rank_genes_groups(ad_sc, groupby="cell_type", method="wilcoxon", n_genes=n_genes, use_raw=False)
    markers_df = pd.DataFrame(ad_sc.uns["rank_genes_groups"]["names"]).iloc[:n_genes, :]
    markers = list(np.unique(markers_df.melt().value.values))
    logger.info(f"Found {len(markers)} unique marker genes")
    return markers


def decompose_gex_per_type(
    ad_map: sc.AnnData,
    ad_sc: sc.AnnData,
    ad_sp: sc.AnnData,
) -> dict:
    """
    Per-cell-type GEX decomposition from Tangram mapping matrix.

    For each cell type k with reference cells T_k:
        GEX_k(s) = sum_{c in T_k} M[c,s] * X_ref[c,:] / sum_{c in T_k} M[c,s]
    If denominator < 1e-6, set GEX to zero (type not present at that spot).

    Returns dict of cell_type -> (n_spots, n_common_genes) numpy arrays, plus
    the list of common gene names.
    """
    # Mapping matrix M: (reference_cells x spots)
    M = ad_map.X  # dense or sparse
    if sp.issparse(M):
        M = M.toarray()
    M = np.asarray(M, dtype=np.float64)

    # Reference expression (use the preprocessed genes from ad_sc)
    X_ref = ad_sc.X
    if sp.issparse(X_ref):
        X_ref = X_ref.toarray()
    X_ref = np.asarray(X_ref, dtype=np.float64)

    # Common genes between spatial and reference
    common_genes = ad_sp.var_names.intersection(ad_sc.var_names)
    sp_gene_idx = ad_sp.var_names.get_indexer(common_genes)
    sc_gene_idx = ad_sc.var_names.get_indexer(common_genes)

    X_ref_common = X_ref[:, sc_gene_idx]  # (n_ref_cells, n_common_genes)

    cell_types = sorted(ad_sc.obs["cell_type"].unique())
    gex_layers = {}

    for ct in cell_types:
        ct_mask = ad_sc.obs["cell_type"].values == ct
        M_k = M[ct_mask, :]  # (n_cells_k, n_spots)
        X_k = X_ref_common[ct_mask, :]  # (n_cells_k, n_common_genes)

        # Weighted sum: (n_spots, n_common_genes)
        denom = M_k.sum(axis=0)  # (n_spots,)
        numer = M_k.T @ X_k  # (n_spots, n_common_genes)

        # Avoid division by zero
        safe_denom = np.where(denom < 1e-6, 1.0, denom)
        gex_k = numer / safe_denom[:, None]
        # Zero out spots where type is absent
        gex_k[denom < 1e-6, :] = 0.0

        gex_layers[ct] = gex_k.astype(np.float32)
        logger.info(f"  {ct}: nonzero spots = {(denom >= 1e-6).sum()}/{len(denom)}")

    return gex_layers, list(common_genes)


def main():
    parser = argparse.ArgumentParser(description="Run Tangram mapping for a single Visium sample")
    parser.add_argument("--sample", required=True, help="Sample name (e.g. HCC22-088-P1-S1)")
    parser.add_argument("--device", default="cuda:0", help="Device for Tangram (default: cuda:0)")
    args = parser.parse_args()

    sample = args.sample
    device = args.device
    out_dir = OUTPUT_ROOT / "tangram" / sample
    out_dir.mkdir(parents=True, exist_ok=True)

    t_start = time.time()

    # ------------------------------------------------------------------
    # Load ER+ reference (already log1p-normalized)
    # ------------------------------------------------------------------
    ref_path = OUTPUT_ROOT / "reference" / "wu_erpos_reference.h5ad"
    logger.info(f"Loading ER+ reference from {ref_path}")
    ad_sc = sc.read_h5ad(ref_path)
    logger.info(f"Reference: {ad_sc.n_obs} cells, {ad_sc.n_vars} genes")
    logger.info(f"Cell types: {ad_sc.obs['cell_type'].value_counts().to_dict()}")

    # ------------------------------------------------------------------
    # Load and normalize Visium sample
    # ------------------------------------------------------------------
    ad_sp = load_visium_sample(sample)
    sc.pp.normalize_total(ad_sp)
    sc.pp.log1p(ad_sp)

    # ------------------------------------------------------------------
    # Find marker genes and preprocess for Tangram
    # ------------------------------------------------------------------
    markers = find_marker_genes(ad_sc, n_genes=100)
    tg.pp_adatas(ad_sc, ad_sp, genes=markers)

    # ------------------------------------------------------------------
    # Run Tangram mapping
    # ------------------------------------------------------------------
    logger.info(f"Running Tangram mapping (mode=cells, 500 epochs, device={device}) ...")
    t_map = time.time()
    ad_map = tg.map_cells_to_space(
        ad_sc,
        ad_sp,
        mode="cells",
        density_prior="rna_count_based",
        num_epochs=500,
        device=device,
    )
    logger.info(f"Mapping completed in {time.time() - t_map:.0f}s")

    # ------------------------------------------------------------------
    # Extract proportions
    # ------------------------------------------------------------------
    logger.info("Projecting cell type annotations ...")
    tg.project_cell_annotations(ad_map, ad_sp, annotation="cell_type")
    ct_pred = ad_sp.obsm["tangram_ct_pred"]  # DataFrame: spots x cell_types
    row_sums = ct_pred.sum(axis=1)
    proportions = ct_pred.div(row_sums, axis=0).fillna(0)

    prop_path = out_dir / "proportions.csv"
    proportions.to_csv(prop_path)
    logger.info(f"Saved proportions to {prop_path}")

    # ------------------------------------------------------------------
    # Per-cell-type GEX decomposition
    # ------------------------------------------------------------------
    logger.info("Computing per-cell-type GEX decomposition ...")
    gex_layers, common_genes = decompose_gex_per_type(ad_map, ad_sc, ad_sp)

    # tg.pp_adatas lowercases gene names — restore to uppercase for consistency
    common_genes = [g.upper() for g in common_genes]

    # Build h5ad with layers
    # Use spatial data's barcodes and the common gene set
    gex_adata = sc.AnnData(
        X=np.zeros((ad_sp.n_obs, len(common_genes)), dtype=np.float32),
        obs=pd.DataFrame(index=ad_sp.obs_names),
        var=pd.DataFrame(index=common_genes),
    )
    # Store spatial coords
    gex_adata.obsm["spatial"] = ad_sp.obsm["spatial"].copy()
    # Store proportions
    gex_adata.obsm["proportions"] = proportions.values.astype(np.float32)
    gex_adata.uns["proportion_columns"] = list(proportions.columns)

    for ct, gex_arr in gex_layers.items():
        gex_adata.layers[ct] = gex_arr

    # Compute reconstruction: sum_k prop_k * GEX_k
    logger.info("Computing weighted reconstruction for QC ...")
    recon = np.zeros((ad_sp.n_obs, len(common_genes)), dtype=np.float64)
    for ct in gex_layers:
        if ct in proportions.columns:
            w = proportions[ct].values[:, None]
            recon += w * gex_layers[ct]
    gex_adata.layers["reconstruction"] = recon.astype(np.float32)

    gex_path = out_dir / "gex_layers.h5ad"
    gex_adata.write(gex_path)
    logger.info(f"Saved GEX layers h5ad to {gex_path}")

    # QC: project genes via Tangram for comparison (non-critical)
    try:
        logger.info("Projecting genes via Tangram for QC comparison ...")
        tg.project_genes(ad_map, ad_sp)
    except Exception as e:
        logger.warning("tg.project_genes QC failed (non-critical): %s", e)

    # ------------------------------------------------------------------
    # Save mapping matrix for provenance
    # ------------------------------------------------------------------
    M = ad_map.X
    if not sp.issparse(M):
        M = sp.csr_matrix(M)
    npz_path = out_dir / "mapping_matrix.npz"
    sp.save_npz(str(npz_path), M)
    logger.info(f"Saved mapping matrix to {npz_path}")

    # ------------------------------------------------------------------
    # QC summary
    # ------------------------------------------------------------------
    elapsed = time.time() - t_start
    logger.info("=" * 60)
    logger.info(f"Tangram mapping complete for {sample}")
    logger.info(f"  Total time: {elapsed:.0f}s")
    logger.info(f"  Spots: {ad_sp.n_obs}")
    logger.info(f"  Common genes in GEX layers: {len(common_genes)}")
    logger.info(f"  Cell types: {sorted(gex_layers.keys())}")
    logger.info(f"  Output dir: {out_dir}")

    # Proportion summary
    logger.info("  Proportion summary (mean per type):")
    for ct in sorted(proportions.columns):
        logger.info(f"    {ct}: {proportions[ct].mean():.4f}")

    logger.info("=" * 60)


if __name__ == "__main__":
    main()
