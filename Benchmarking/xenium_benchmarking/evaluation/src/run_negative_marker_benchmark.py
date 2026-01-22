#!/usr/bin/env python3
"""
Benchmark CITEgeist with negative marker autodiscovery on Xenium data.

This script:
1. Loads Xenium pseudo-Visium data
2. Runs Module 2d to autodiscover profiles with negative markers
3. Runs Module 3 with negative marker optimization
4. Evaluates against validated GT

Usage:
    python run_negative_marker_benchmark.py --region 0 --output_dir /path/to/output
"""

import argparse
import json
import logging
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats

# Add paths
REPO_ROOT = Path(__file__).parent.parent.parent.parent.parent
sys.path.insert(0, str(REPO_ROOT / "Benchmarking/xenium_pseudovisium/src"))
sys.path.insert(0, str(REPO_ROOT / "CITEgeist"))

from load_xenium import load_xenium_data, split_gex_protein

# Import CITEgeist modules
from CITEgeist.model.spatial_colocalization import (
    analyze_marker_colocalization,
    discover_profiles,
    detect_residual_signal,
    resolve_hierarchical_profiles,
)
from CITEgeist.model.gurobi_impl import (
    optimize_cell_proportions_with_negatives,
    build_assignment_matrices_from_profiles,
)

logger = logging.getLogger(__name__)

# Paths
XENIUM_DATA_DIR = Path("/ix1/alee/LO_LAB/General/Public_Data/10x_PublicData/Xenium_RNA_Proteomic_RenalCellCarcinoma")
PSEUDOVISIUM_DIR = Path("/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_pseudovisium")
BENCHMARK_DIR = Path("/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking")

# Validated GT path
VALIDATED_GT_DIR = BENCHMARK_DIR / "evaluation" / "validated_benchmark"


def load_pseudovisium_region(region_id: int) -> Tuple[sc.AnnData, sc.AnnData]:
    """Load pseudo-Visium data for a region."""
    data_dir = PSEUDOVISIUM_DIR / "data"

    # Load GEX
    gex_path = data_dir / f"pseudovisium_region_{region_id}_gex.h5ad"
    adata_gex = sc.read_h5ad(gex_path)

    # Load protein
    protein_path = data_dir / f"pseudovisium_region_{region_id}_protein.h5ad"
    adata_protein = sc.read_h5ad(protein_path)

    return adata_gex, adata_protein


def load_validated_gt(region_id: int) -> pd.DataFrame:
    """Load validated ground truth for a region."""
    gt_path = VALIDATED_GT_DIR / f"validated_gt_region_{region_id}.csv"
    if not gt_path.exists():
        raise FileNotFoundError(f"Validated GT not found: {gt_path}")
    return pd.read_csv(gt_path, index_col=0)


def run_autodiscovery_with_negatives(
    adata_protein: sc.AnnData,
    fdr_alpha: float = 0.05,
    residual_threshold: float = 0.10,
) -> Dict[str, Dict[str, List[str]]]:
    """
    Run Module 2 autodiscovery with negative marker resolution.

    Returns:
        Profile dict with positive and negative markers
    """
    # Get protein data
    X = adata_protein.X
    if hasattr(X, 'toarray'):
        X = X.toarray()

    marker_names = list(adata_protein.var_names)
    coords = adata_protein.obsm['spatial']

    logger.info(f"Running autodiscovery on {X.shape[0]} spots, {X.shape[1]} markers")

    # Normalize to [0, 1] range
    X_norm = X.copy()
    for i in range(X_norm.shape[1]):
        col_min = X_norm[:, i].min()
        col_max = X_norm[:, i].max()
        if col_max > col_min:
            X_norm[:, i] = (X_norm[:, i] - col_min) / (col_max - col_min)

    # Module 2a: Colocalization analysis
    logger.info("Running Module 2a: Colocalization analysis...")
    coloc_result = analyze_marker_colocalization(
        X_norm, coords, marker_names,
        neighbor_k=6,
        n_permutations=500,
        seed=42,
    )

    # Module 2b: Profile discovery
    logger.info("Running Module 2b: Profile discovery...")
    profile_result = discover_profiles(
        coloc_result,
        fdr_alpha=fdr_alpha,
        pvalue_source="bivariate_morans",
        verbose=True,
    )

    logger.info(f"Discovered {len(profile_result.profiles)} initial profiles")
    for i, profile in enumerate(profile_result.profiles):
        logger.info(f"  Profile {i}: {profile}")

    # Module 2d: Hierarchical resolution with negative markers
    logger.info("Running Module 2d: Hierarchical profile resolution...")
    resolved_result = resolve_hierarchical_profiles(
        X_norm, marker_names, profile_result.profiles,
        residual_threshold=residual_threshold,
        verbose=True,
    )

    logger.info(f"Resolved to {resolved_result.n_resolved_profiles} profiles with negatives")

    # Convert to profile dict format
    profile_dict = resolved_result.to_profile_dict()

    return profile_dict


def run_citegeist_with_negatives(
    adata_protein: sc.AnnData,
    profile_dict: Dict[str, Dict[str, List[str]]],
    lambda_neg: float = 1.0,
) -> pd.DataFrame:
    """
    Run CITEgeist Module 3 with negative marker optimization.

    Returns:
        DataFrame of cell type proportions per spot
    """
    # Get protein data
    X = adata_protein.X
    if hasattr(X, 'toarray'):
        X = X.toarray()

    marker_names = list(adata_protein.var_names)
    coords = adata_protein.obsm['spatial']

    # Normalize
    X_norm = X.copy()
    for i in range(X_norm.shape[1]):
        col_min = X_norm[:, i].min()
        col_max = X_norm[:, i].max()
        if col_max > col_min:
            X_norm[:, i] = (X_norm[:, i] - col_min) / (col_max - col_min)

    # Build assignment matrices
    pos_assign, neg_assign, cell_type_names = build_assignment_matrices_from_profiles(
        marker_names, profile_dict
    )

    logger.info(f"Positive assignments: {pos_assign.sum():.0f}")
    logger.info(f"Negative assignments: {neg_assign.sum():.0f}")
    logger.info(f"Cell types: {cell_type_names}")

    # Run optimization
    logger.info("Running Module 3 optimization with negative markers...")
    Y, beta, marker_beta = optimize_cell_proportions_with_negatives(
        X_norm,
        marker_names,
        pos_assign,
        neg_assign,
        cell_type_names,
        lambda_neg=lambda_neg,
        lambda_reg=0.1,
        lambda_laplacian=0.1,
        coords=coords,
        max_iterations=30,
        warn_only=True,
    )

    # Create DataFrame
    proportions_df = pd.DataFrame(
        Y,
        index=adata_protein.obs_names,
        columns=cell_type_names,
    )

    return proportions_df


def evaluate_proportions(
    pred_df: pd.DataFrame,
    gt_df: pd.DataFrame,
) -> Dict:
    """Evaluate predicted proportions against ground truth."""
    # Align indices
    common_spots = pred_df.index.intersection(gt_df.index)

    if len(common_spots) == 0:
        raise ValueError("No common spots between prediction and GT")

    pred = pred_df.loc[common_spots]
    gt = gt_df.loc[common_spots]

    # Find common cell types
    common_types = [c for c in pred.columns if c in gt.columns]

    results = {
        "n_spots": len(common_spots),
        "n_types": len(common_types),
        "cell_types": common_types,
        "per_type": {},
    }

    # Per-type metrics
    all_pred = []
    all_gt = []

    for cell_type in common_types:
        p = pred[cell_type].values
        g = gt[cell_type].values

        all_pred.extend(p)
        all_gt.extend(g)

        r, pval = stats.pearsonr(p, g)
        rmse = np.sqrt(np.mean((p - g) ** 2))
        mae = np.mean(np.abs(p - g))

        results["per_type"][cell_type] = {
            "pearson_r": float(r),
            "pearson_p": float(pval),
            "rmse": float(rmse),
            "mae": float(mae),
            "mean_gt": float(g.mean()),
            "mean_pred": float(p.mean()),
        }

    # Overall metrics
    all_pred = np.array(all_pred)
    all_gt = np.array(all_gt)

    results["overall_pearson_r"], results["overall_pearson_p"] = stats.pearsonr(all_pred, all_gt)
    results["overall_rmse"] = float(np.sqrt(np.mean((all_pred - all_gt) ** 2)))
    results["overall_mae"] = float(np.mean(np.abs(all_pred - all_gt)))

    return results


def main():
    parser = argparse.ArgumentParser(description="Run negative marker benchmark")
    parser.add_argument("--region", type=int, required=True, help="Region ID (0-4)")
    parser.add_argument("--output_dir", type=str, required=True, help="Output directory")
    parser.add_argument("--lambda_neg", type=float, default=1.0, help="Negative marker penalty")
    parser.add_argument("--residual_threshold", type=float, default=0.10, help="Residual threshold")
    parser.add_argument("--fdr_alpha", type=float, default=0.05, help="FDR threshold")
    args = parser.parse_args()

    # Setup logging
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info(f"=" * 60)
    logger.info(f"Negative Marker Benchmark - Region {args.region}")
    logger.info(f"=" * 60)
    logger.info(f"Parameters:")
    logger.info(f"  lambda_neg: {args.lambda_neg}")
    logger.info(f"  residual_threshold: {args.residual_threshold}")
    logger.info(f"  fdr_alpha: {args.fdr_alpha}")

    # Load data
    logger.info("Loading pseudo-Visium data...")
    adata_gex, adata_protein = load_pseudovisium_region(args.region)
    logger.info(f"  GEX: {adata_gex.shape}")
    logger.info(f"  Protein: {adata_protein.shape}")

    # Load validated GT
    logger.info("Loading validated ground truth...")
    gt_df = load_validated_gt(args.region)
    logger.info(f"  GT shape: {gt_df.shape}")
    logger.info(f"  GT cell types: {list(gt_df.columns)}")

    # Run autodiscovery with negative markers
    logger.info("-" * 60)
    logger.info("Running Module 2 autodiscovery with negative markers...")
    profile_dict = run_autodiscovery_with_negatives(
        adata_protein,
        fdr_alpha=args.fdr_alpha,
        residual_threshold=args.residual_threshold,
    )

    # Save discovered profiles
    profiles_path = output_dir / f"discovered_profiles_region_{args.region}.json"
    with open(profiles_path, 'w') as f:
        json.dump(profile_dict, f, indent=2)
    logger.info(f"Saved discovered profiles to {profiles_path}")

    # Run CITEgeist with negative markers
    logger.info("-" * 60)
    logger.info("Running Module 3 with negative marker optimization...")
    pred_df = run_citegeist_with_negatives(
        adata_protein,
        profile_dict,
        lambda_neg=args.lambda_neg,
    )

    # Save predictions
    pred_path = output_dir / f"citegeist_negmarkers_region_{args.region}.csv"
    pred_df.to_csv(pred_path)
    logger.info(f"Saved predictions to {pred_path}")

    # Evaluate
    logger.info("-" * 60)
    logger.info("Evaluating against validated GT...")
    results = evaluate_proportions(pred_df, gt_df)

    # Print results
    logger.info(f"\nResults for Region {args.region}:")
    logger.info(f"  Overall Pearson r: {results['overall_pearson_r']:.3f}")
    logger.info(f"  Overall RMSE: {results['overall_rmse']:.3f}")
    logger.info(f"  Overall MAE: {results['overall_mae']:.3f}")
    logger.info(f"\n  Per-type metrics:")

    for cell_type, metrics in results["per_type"].items():
        logger.info(f"    {cell_type}:")
        logger.info(f"      r={metrics['pearson_r']:.3f}, "
                   f"mean_gt={metrics['mean_gt']:.3f}, "
                   f"mean_pred={metrics['mean_pred']:.3f}")

    # Save results
    results["region_id"] = args.region
    results["lambda_neg"] = args.lambda_neg
    results["residual_threshold"] = args.residual_threshold
    results["discovered_profiles"] = profile_dict

    results_path = output_dir / f"negmarker_benchmark_region_{args.region}.json"
    with open(results_path, 'w') as f:
        json.dump(results, f, indent=2, default=lambda x: float(x) if hasattr(x, 'item') else x)
    logger.info(f"Saved results to {results_path}")

    logger.info("=" * 60)
    logger.info("Benchmark complete!")
    logger.info("=" * 60)


if __name__ == "__main__":
    main()
