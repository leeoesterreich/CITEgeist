#!/usr/bin/env python
"""
Task 7: COMMOT spatial communication for competitor GEX validation.

Runs COMMOT CellChat-based ligand-receptor communication analysis on
pseudo-single-cell data, extracting MDK-related pathway scores.

Usage:
    python -u run_commot.py --method cell2location --sample HCC22-088-P1-S1
"""

import argparse
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

sys.path.insert(0, str(Path(__file__).resolve().parent))
from constants import COMMOT_DISTANCE_THR, COMMOT_PATHWAYS, OUTPUT_ROOT

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
)
logger = logging.getLogger(__name__)


def main():
    parser = argparse.ArgumentParser(description="COMMOT MDK signaling analysis")
    parser.add_argument("--method", required=True, choices=["cell2location", "tangram"])
    parser.add_argument("--sample", required=True)
    args = parser.parse_args()

    method = args.method
    sample = args.sample

    # Paths
    sc_path = OUTPUT_ROOT / method / sample / f"{sample}_single_cell.h5ad"
    out_dir = OUTPUT_ROOT / method / "commot"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{sample}_commot_mdk.csv"

    if not sc_path.exists():
        logger.error("Single-cell h5ad not found: %s", sc_path)
        sys.exit(1)

    logger.info("=== COMMOT: method=%s sample=%s ===", method, sample)

    # Load pseudo-single-cell data
    adata = sc.read_h5ad(sc_path)
    adata.var_names_make_unique()
    logger.info("Loaded %d cells x %d genes", adata.n_obs, adata.n_vars)

    # Filter out low-UMI cells (empty pseudo-cells confuse COMMOT)
    total_umi = np.array(adata.X.sum(axis=1)).flatten()
    umi_mask = total_umi > 20
    n_removed = (~umi_mask).sum()
    if n_removed > 0:
        logger.info("Filtering %d/%d cells with <=20 UMI", n_removed, adata.n_obs)
        adata = adata[umi_mask].copy()
        logger.info("After filter: %d cells", adata.n_obs)

    # COMMOT expects "celltype" not "cell_type"
    adata.obs["celltype"] = adata.obs["cell_type"].values

    # Import COMMOT
    try:
        import commot as ct
    except ImportError:
        logger.error("COMMOT not available. Run in COMMOT conda environment.")
        sys.exit(1)

    # Load CellChat database
    logger.info("Loading CellChat L-R database...")
    df_cellchat = ct.pp.ligand_receptor_database(species="human", database="CellChat")

    # Pre-filter to only target ligands before expression filter (saves memory)
    # CellChat DB columns are integers: 0=ligand, 1=receptor, 2=pathway, 3=signaling_type
    target_ligands = {"MDK", "PTN", "MIF"}
    lig_col = df_cellchat.columns[0]
    df_cellchat_sub = df_cellchat[df_cellchat[lig_col].isin(target_ligands)]
    logger.info(
        "Pre-filtered to %d L-R pairs with target ligands %s (col=%s)", len(df_cellchat_sub), target_ligands, lig_col
    )

    # Filter L-R pairs to those expressed in this dataset
    df_filtered = ct.pp.filter_lr_database(df_cellchat_sub, adata, min_cell_pct=0.01)
    logger.info("Filtered to %d L-R pairs (from %d)", len(df_filtered), len(df_cellchat_sub))

    if len(df_filtered) == 0:
        logger.warning("No L-R pairs pass filter — saving empty results")
        empty_df = pd.DataFrame(columns=["spot_barcode", "cell_type"] + COMMOT_PATHWAYS + ["method", "sample"])
        empty_df.to_csv(out_path, index=False)
        return

    # Run spatial communication
    logger.info("Running COMMOT spatial communication (dis_thr=%d)...", COMMOT_DISTANCE_THR)
    ct.tl.spatial_communication(
        adata,
        database_name="cellchat",
        df_ligrec=df_filtered,
        dis_thr=COMMOT_DISTANCE_THR,
        heteromeric=True,
        pathway_sum=True,
    )

    # Extract MDK pathway scores from sum-sender DataFrame
    # COMMOT stores sender marginals in obsm['commot-cellchat-sum-sender'] as a DataFrame
    # with columns like 's-MDK-SDC4', 's-MK' (pathway sums), 's-total-total'
    logger.info("Extracting MDK pathway scores...")
    sender_key = "commot-cellchat-sum-sender"
    if sender_key not in adata.obsm:
        logger.error("No sender scores found in obsm. Keys: %s", list(adata.obsm.keys()))
        sys.exit(1)

    sender_df = adata.obsm[sender_key]
    logger.info("  Available sender columns: %s", list(sender_df.columns))

    result_df = pd.DataFrame(
        {
            "spot_barcode": (
                adata.obs["spot_barcode"].values if "spot_barcode" in adata.obs.columns else adata.obs_names
            ),
            "cell_type": adata.obs["celltype"].values,
        }
    )

    for pathway in COMMOT_PATHWAYS:
        if pathway in sender_df.columns:
            result_df[pathway] = sender_df[pathway].values
            vals = result_df[pathway]
            logger.info(
                "  %s: mean=%.4f, max=%.4f, nonzero=%d/%d",
                pathway,
                vals.mean(),
                vals.max(),
                (vals > 0).sum(),
                len(vals),
            )
        else:
            logger.warning("  %s: not found in sender columns", pathway)
            result_df[pathway] = np.nan

    result_df["method"] = method
    result_df["sample"] = sample

    # Log dominant sender cell type per pathway
    logger.info("Dominant sender cell types per pathway:")
    for pathway in COMMOT_PATHWAYS:
        if pathway in result_df.columns and result_df[pathway].notna().any():
            agg = result_df.groupby("cell_type")[pathway].mean()
            dominant = agg.idxmax()
            logger.info("  %s: %s (mean=%.4f)", pathway, dominant, agg[dominant])

    result_df.to_csv(out_path, index=False)
    logger.info("Saved %d rows to %s", len(result_df), out_path)


if __name__ == "__main__":
    main()
