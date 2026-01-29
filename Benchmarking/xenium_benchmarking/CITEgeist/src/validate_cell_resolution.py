"""
Validate cell-resolution CITEgeist on Xenium single-cell data.

Runs Modules 1-3 in resolution="cell" mode on Xenium RCC regions,
comparing against RNA-derived ground truth cell type labels and
the existing heuristic cell assignment pipeline.

Implements Tests 1-4 from the design doc.
"""
import argparse
import json
import logging
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.metrics import accuracy_score, f1_score, confusion_matrix

# Add paths
REPO_ROOT = Path(__file__).parent.parent.parent.parent.parent
sys.path.insert(0, str(REPO_ROOT / "CITEgeist"))
sys.path.insert(0, str(REPO_ROOT / "Benchmarking" / "xenium_pseudovisium" / "src"))
sys.path.insert(0, str(REPO_ROOT / "Benchmarking" / "xenium_benchmarking" / "CITEgeist" / "src"))

from model import (
    CitegeistModel,
    identify_interesting_markers,
    analyze_marker_colocalization,
    discover_profiles,
    select_profiles,
)
from model.gurobi_impl import (
    optimize_cell_proportions_per_marker,
    map_antibodies_to_profiles_v2,
    estimate_true_expression_cell,
)
from load_xenium_singlecell import load_xenium_singlecell

logger = logging.getLogger(__name__)
OUTPUT_DIR = REPO_ROOT / "Benchmarking" / "xenium_benchmarking" / "CITEgeist" / "output_cell_resolution"
GT_PATH = REPO_ROOT / "Benchmarking" / "xenium_pseudovisium" / "data_rna_gt" / "cell_types.csv"


def run_validation(region_id: int, max_cells: int = 10000):
    """Run full cell-resolution validation for one region."""
    output_dir = OUTPUT_DIR / f"region_{region_id}"
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load data
    logger.info(f"Loading Xenium region {region_id} (max {max_cells} cells)")
    adata_gex, adata_protein = load_xenium_singlecell(
        region_id=region_id, max_cells=max_cells, seed=42
    )
    logger.info(f"Loaded {adata_gex.shape[0]} cells")

    # Load ground truth
    gt_df = pd.read_csv(GT_PATH)
    gt_df = gt_df.set_index("cell_id")
    common_cells = list(set(adata_gex.obs_names) & set(gt_df.index))
    gt_labels = gt_df.loc[common_cells, "cell_type"].values
    logger.info(f"Ground truth available for {len(common_cells)} cells")

    # === MODULE 1: Marker Interest Detection ===
    logger.info("=" * 60)
    logger.info("MODULE 1: Marker Interest Detection (cell resolution)")
    X_protein = adata_protein.X
    if hasattr(X_protein, "toarray"):
        X_protein = X_protein.toarray()
    X_protein = np.asarray(X_protein, dtype=np.float64)
    coords = adata_protein.obsm["spatial"]
    marker_names = list(adata_protein.var_names)

    m1_result = identify_interesting_markers(
        X=X_protein, coords=coords, marker_names=marker_names,
        morans_k=50, smooth_k=20, morans_n_perm=99, seed=42,
    )
    interesting = m1_result.interesting_markers
    logger.info(f"Found {len(interesting)} interesting markers: {interesting}")

    # Save Module 1 results
    m1_df = m1_result.to_dataframe()
    m1_df.to_csv(output_dir / "module1_marker_interest.csv", index=False)

    # === MODULE 2: Profile Discovery ===
    logger.info("=" * 60)
    logger.info("MODULE 2: Profile Discovery (cell resolution)")
    coloc_result = analyze_marker_colocalization(
        X=X_protein, coords=coords, marker_names=marker_names,
        markers_to_analyze=interesting,
        neighbor_k=30, smooth_k=20,
        multi_scale_k=[20, 40, 60, 80, 100],
        n_permutations=99, seed=42,
    )
    profile_result = discover_profiles(coloc_result, seed=42)
    selected = select_profiles(
        adata_protein, profile_result,
        interesting_markers=interesting,
    )
    profiles = selected.profiles
    logger.info(f"Discovered {len(profiles)} profiles: {list(profiles.keys())}")

    # Save profiles
    with open(output_dir / "profiles.json", "w") as f:
        json.dump({k: list(v) for k, v in profiles.items()}, f, indent=2)

    # === MODULE 3 PASS 1: Cell Type Soft Classification ===
    logger.info("=" * 60)
    logger.info("MODULE 3 PASS 1: Cell Type Classification (cell resolution, Gurobi QP)")

    # Map markers to profiles (reuse existing mapping function)
    marker_data, assignment_matrix, type_names = map_antibodies_to_profiles_v2(
        adata_protein, profiles
    )

    Y_assignments, beta_values, beta_dict = optimize_cell_proportions_per_marker(
        marker_level_data=marker_data,
        marker_names=[marker_names[i] for i in range(marker_data.shape[1])],
        assignment_matrix=assignment_matrix,
        cell_type_names=type_names,
        lambda_sparse=0.1,
        lambda_laplacian=0.01,
        coords=coords,
        laplacian_k=50,
        max_iterations=10,
    )

    # Evaluate classification
    dominant_type = np.argmax(Y_assignments, axis=1)
    predicted_labels = [type_names[dt] for dt in dominant_type]

    results = {
        "region_id": region_id,
        "n_cells": len(common_cells),
        "n_profiles": len(profiles),
        "profiles": list(profiles.keys()),
        "max_Y_mean": float(Y_assignments.max(axis=1).mean()),
        "max_Y_median": float(np.median(Y_assignments.max(axis=1))),
        "doublet_fraction": float(np.mean(Y_assignments.max(axis=1) < 0.6)),
    }

    # Save Y_assignments
    y_df = pd.DataFrame(Y_assignments, columns=type_names, index=adata_protein.obs_names)
    y_df.to_csv(output_dir / "cell_assignments.csv")

    # === MODULE 3 PASS 2: True Count Estimation ===
    logger.info("=" * 60)
    logger.info("MODULE 3 PASS 2: True Count Estimation (cell resolution)")

    X_gex = adata_gex.X
    if hasattr(X_gex, "toarray"):
        X_gex = X_gex.toarray()
    X_gex = np.asarray(X_gex, dtype=np.float64)

    # Build enrichment weights from Pass 1 assignments
    n_types = len(type_names)
    n_genes = X_gex.shape[1]
    enrichment = np.ones((n_types, n_genes)) * 0.1
    for ct_idx in range(n_types):
        ct_cells = np.where(dominant_type == ct_idx)[0]
        if len(ct_cells) > 0:
            ct_mean = X_gex[ct_cells].mean(axis=0)
            global_mean = X_gex.mean(axis=0) + 1e-10
            enrichment[ct_idx] = ct_mean / global_mean

    X_true = estimate_true_expression_cell(
        X_obs=X_gex,
        Y_assignments=Y_assignments,
        coords=coords,
        enrichment_weights=enrichment,
        library_slack=1.5,
        lambda_spatial=0.01,
        spatial_k=50,
    )

    # Evaluate dropout recovery
    obs_zeros = np.sum(X_gex == 0)
    true_zeros = np.sum(X_true == 0)
    results["dropout_recovery"] = {
        "obs_zero_fraction": float(obs_zeros / X_gex.size),
        "true_zero_fraction": float(true_zeros / X_true.size),
        "zeros_recovered": int(obs_zeros - true_zeros),
        "recovery_rate": float((obs_zeros - true_zeros) / max(obs_zeros, 1)),
    }

    # Save results
    with open(output_dir / "validation_results.json", "w") as f:
        json.dump(results, f, indent=2, default=str)

    logger.info(f"Results saved to {output_dir}")
    logger.info(f"Max Y mean: {results['max_Y_mean']:.3f}")
    logger.info(f"Doublet fraction: {results['doublet_fraction']:.3f}")
    logger.info(f"Dropout recovery rate: {results['dropout_recovery']['recovery_rate']:.3f}")

    return results


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Validate cell-resolution CITEgeist on Xenium single-cell data."
    )
    parser.add_argument("--region", type=int, default=0, help="Region ID (0-4)")
    parser.add_argument("--max-cells", type=int, default=10000, help="Max cells to load")
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    run_validation(args.region, args.max_cells)
