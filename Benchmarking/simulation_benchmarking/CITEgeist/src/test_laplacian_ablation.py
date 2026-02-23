#!/usr/bin/env python3
"""Test whether Laplacian spatial smoothing is the key difference maker."""

import sys
import os
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import pearsonr

# Add CITEgeist to path
sys.path.insert(0, "/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist")

from CITEgeist.model import CitegeistModel

# Simulation cell profile dict (same as benchmark - with nested Major keys)
SIMULATION_CELL_PROFILE_DICT = {
    "B-cells": {"Major": ["B-cells_Protein_1", "B-cells_Protein_2"]},
    "CAFs": {"Major": ["CAFs_Protein_1", "CAFs_Protein_2"]},
    "Cancer Epithelial": {"Major": ["Cancer Epithelial_Protein_1", "Cancer Epithelial_Protein_2"]},
    "Endothelial": {"Major": ["Endothelial_Protein_1", "Endothelial_Protein_2"]},
    "Myeloid": {"Major": ["Myeloid_Protein_1", "Myeloid_Protein_2"]},
    "Normal Epithelial": {"Major": ["Normal Epithelial_Protein_1", "Normal Epithelial_Protein_2"]},
    "PVL": {"Major": ["PVL_Protein_1", "PVL_Protein_2"]},
    "Plasmablasts": {"Major": ["Plasmablasts_Protein_1", "Plasmablasts_Protein_2"]},
    "T-cells": {"Major": ["T-cells_Protein_1", "T-cells_Protein_2"]},
}


def run_continuous_benchmark(replicate_id: int, condition: str, lambda_laplacian: float):
    """Run continuous model and compute proportion correlation."""

    # Data is in Brent's scCube simulation directory
    sccube_dir = "/ix1/alee/LO_LAB/Personal/Brent_Schlegel/bts76/Projects/CITEgeist/Wu_Visium/Simulations/scCube_12k/replicates"

    # Load data
    gex_path = f"{sccube_dir}/{condition}/h5ad_objects/Wu_rep_{replicate_id}_GEX.h5ad"
    cite_path = f"{sccube_dir}/{condition}/h5ad_objects/Wu_rep_{replicate_id}_CITE.h5ad"
    gt_path = f"{sccube_dir}/{condition}/ST_sim/Wu_ST_{replicate_id}_prop.csv"

    print(f"Loading replicate {replicate_id}, condition={condition}")
    adata_gex = sc.read_h5ad(gex_path)
    adata_cite = sc.read_h5ad(cite_path)
    gt_props = pd.read_csv(gt_path, index_col=0)
    # Drop coordinate columns if present
    gt_props = gt_props.drop(columns=["spot_x", "spot_y"], errors="ignore")

    # Initialize model
    model = CitegeistModel(
        sample_name=f"Wu_rep_{replicate_id}",
        output_folder=f"/tmp/laplacian_test_{replicate_id}_{lambda_laplacian}",
        simulation=True,
        gene_expression_adata=adata_gex,
        antibody_capture_adata=adata_cite,
    )

    # Preprocess (standard continuous mode)
    model.preprocess_antibody()
    model.preprocess_gex()
    model.load_cell_profile_dict(SIMULATION_CELL_PROFILE_DICT)

    # Run continuous proportion model
    print(f"Running continuous model with lambda_laplacian={lambda_laplacian}")
    global_props_df, finetuned_props_df = model.run_cell_proportion_model(
        max_iterations=20,
        lambda_laplacian=lambda_laplacian,  # Key parameter we're testing
    )

    # Align columns
    common_types = [c for c in gt_props.columns if c in finetuned_props_df.columns]
    pred = finetuned_props_df[common_types].values.flatten()
    gt = gt_props[common_types].values.flatten()

    # Compute correlation
    mask = ~(np.isnan(pred) | np.isnan(gt))
    corr, _ = pearsonr(pred[mask], gt[mask])

    return corr


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--replicate-id", type=int, default=0)
    parser.add_argument("--condition", type=str, default="mixed")
    parser.add_argument("--lambda-laplacian", type=float, default=0.1)
    args = parser.parse_args()

    corr = run_continuous_benchmark(
        args.replicate_id,
        args.condition,
        args.lambda_laplacian
    )

    print(f"\n=== RESULT ===")
    print(f"Condition: {args.condition}")
    print(f"Replicate: {args.replicate_id}")
    print(f"lambda_laplacian: {args.lambda_laplacian}")
    print(f"Proportion correlation: {corr:.4f}")
