#!/usr/bin/env python
"""Compute COMMOT MDK-SDC4 sender scores for Figure 4E.

This script runs in the COMMOT environment and saves the sender scores
to a CSV file that can be loaded by generate_figure4.py.

Usage:
    conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/COMMOT
    python compute_commot_scores.py
"""

import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq

# Add CITEgeist to path
PROJECT_ROOT = Path("/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist")
sys.path.insert(0, str(PROJECT_ROOT / "CITEgeist"))

# Import COMMOT
import commot as ct

# Constants
SAMPLE = "HCC22-088-P4-S2_1i_rep"
DATA_ROOT = Path("/ix1/alee/LO_LAB/General/Lab_Data/20250210_CITEGeistPublicData_GEO_Alex/processed_files")
OUTPUT_DIR = PROJECT_ROOT / "manuscript" / "figures" / "output"

# Cell type profiles (same as vignette 2)
CELL_PROFILES = {
    "Endothelial": ["CD31-1"],
    "Fibroblasts": ["aSMA-1"],
    "B_Cells": ["CD20-1"],
    "Macrophages": ["CD68-1", "CD163-1"],
    "Monocytes": ["CD14-1"],
    "CD8_T_Cells": ["CD3E-1", "CD8A-1"],
    "CD4_T_Cells": ["CD3E-1", "CD4-1"],
    "Cancer_Luminal": ["EPCAM-1", "panCK-1"],
    "Cancer_Basal": ["EPCAM-1", "panCK-1"],
    "Dendritic_Cells": ["CD11c-1"],
}


def main():
    print("=" * 60)
    print("Computing COMMOT MDK-SDC4 Sender Scores")
    print("=" * 60)

    # Load Visium data
    visium_path = DATA_ROOT / SAMPLE / "outs"
    print(f"\nLoading Visium data from {visium_path}...")
    adata = sq.read.visium(
        path=str(visium_path),
        counts_file="filtered_feature_bc_matrix.h5",
        load_images=True,
        gex_only=False,
    )
    print(f"  Loaded: {adata.shape}")

    # Load deconvolved gene expression
    gex_path = PROJECT_ROOT / "output" / "module3_unified" / f"{SAMPLE}_gene_expression_pass1.parquet"
    print(f"\nLoading deconvolved GEX from {gex_path}...")
    gex_df = pd.read_parquet(gex_path)
    print(f"  Loaded: {gex_df.shape}")

    # Load cell proportions
    prop_path = PROJECT_ROOT / "output" / "module3_unified" / f"{SAMPLE}_cell_prop_finetuned_results.csv"
    print(f"\nLoading cell proportions from {prop_path}...")
    prop_df = pd.read_csv(prop_path, index_col=0)
    print(f"  Loaded: {prop_df.shape}")

    # Import expand function
    from model.analysis_functions import expand_prop_gex_adata

    # Get cell types from gex_df
    cell_types = gex_df.index.str.split(":::").str[1].unique()
    print(f"\nCell types in gex_df: {list(cell_types)}")

    # Align barcodes
    common_barcodes = adata.obs_names.intersection(
        gex_df.index.str.split(":::").str[0].unique()
    )
    print(f"  Common barcodes: {len(common_barcodes)}")

    adata_sub = adata[common_barcodes].copy()

    # Attach deconvolved gene expression layers
    print("\nAttaching deconvolved GEX layers...")
    for ct_name in cell_types:
        ct_mask = gex_df.index.str.endswith(f":::{ct_name}")
        ct_df = gex_df.loc[ct_mask].copy()
        ct_df.index = ct_df.index.str.split(":::").str[0]
        ct_df = ct_df.loc[ct_df.index.intersection(common_barcodes)]

        common_genes = adata_sub.var_names.intersection(ct_df.columns)
        if len(common_genes) == 0:
            continue

        layer_name = ct_name.replace(" ", "_") + "_genes_pass1"
        adata_sub.layers[layer_name] = ct_df.loc[adata_sub.obs_names, common_genes].values
        print(f"  Added layer {layer_name}: {len(common_genes)} genes")

    # Add proportions to obs
    for col in prop_df.columns:
        if col in ["barcode"]:
            continue
        adata_sub.obs[col] = prop_df.loc[adata_sub.obs_names, col].values

    # Expand to pseudo-single-cell
    print("\nExpanding to pseudo-single-cell representation...")
    expanded_adata = expand_prop_gex_adata(adata_sub, celltype_profile_dict=CELL_PROFILES)
    print(f"  Expanded adata: {expanded_adata.shape}")

    # Get ligand-receptor database
    print("\nLoading CellChat ligand-receptor database...")
    df_ligrec = ct.pp.ligand_receptor_database(database="CellChat", species="human")
    df_cellchat_filtered = ct.pp.filter_lr_database(df_ligrec, expanded_adata, min_cell_pct=0.01)
    print(f"  Filtered L-R pairs: {len(df_cellchat_filtered)}")

    # Run COMMOT spatial communication
    print("\nRunning COMMOT spatial_communication...")
    ct.tl.spatial_communication(
        expanded_adata,
        database_name="cellchat",
        df_ligrec=df_cellchat_filtered,
        dis_thr=500,
        heteromeric=True,
        pathway_sum=True,
    )
    print("  COMMOT analysis complete")

    # Extract sender scores
    print("\nExtracting sender scores...")
    sender_signal = expanded_adata.obsm["commot-cellchat-sum-sender"]
    print(f"  Sender signal columns: {list(sender_signal.columns)}")

    # Filter to cancer cells only
    cancer_mask = expanded_adata.obs["cell_type"].isin(["Cancer_Luminal", "Cancer_Basal"])
    cancer_adata = expanded_adata[cancer_mask].copy()
    cancer_sender = cancer_adata.obsm["commot-cellchat-sum-sender"]
    print(f"  Cancer cells: {cancer_adata.shape[0]}")

    # Aggregate to spot level (mean across cancer cell types per spot)
    cancer_sender["barcode"] = cancer_adata.obs["barcode"].values
    spot_sender = cancer_sender.groupby("barcode").mean()
    print(f"  Aggregated to {len(spot_sender)} spots")

    # Save results
    output_file = OUTPUT_DIR / "commot_sender_scores_P4-S2.csv"
    spot_sender.to_csv(output_file)
    print(f"\nSaved sender scores to {output_file}")

    # Print summary statistics for key pathways
    print("\n" + "=" * 60)
    print("Sender Score Summary (cancer cells)")
    print("=" * 60)
    for col in ["s-MDK-SDC4", "s-MDK-NCL", "s-PTN-SDC4", "s-PTN-NCL", "s-MIF-CD74_CD44"]:
        if col in spot_sender.columns:
            print(f"  {col}: mean={spot_sender[col].mean():.3f}, max={spot_sender[col].max():.3f}")

    print("\nDone!")


if __name__ == "__main__":
    main()
