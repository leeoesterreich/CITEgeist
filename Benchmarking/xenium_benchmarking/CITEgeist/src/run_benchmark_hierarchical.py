#!/usr/bin/env python
"""
Mode C: Hierarchical Profile Discovery Benchmark

Runs the CITEgeist pipeline with HIERARCHICAL profile discovery:
1. Module 1: Identify interesting markers
2. Module 2a: Pairwise colocalization analysis
3. Module 2b: HIERARCHICAL profile discovery (learns shared markers like CD3)
4. Module 2c: Profile selection (if needed)
5. Module 3: Cell proportion optimization
6. Evaluation against RNA-based ground truth

This mode tests CITEgeist's ability to:
- Learn hierarchical marker relationships (e.g., CD3 shared by CD4+ and CD8+ T cells)
- Properly allocate shared markers across related cell types
- Achieve improved performance on data with marker hierarchies
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
from scipy.stats import pearsonr, spearmanr
from scipy.spatial.distance import jensenshannon

# Add CITEgeist to path
REPO_ROOT = Path(__file__).parent.parent.parent.parent.parent
sys.path.insert(0, str(REPO_ROOT))

from CITEgeist.model.marker_interest import identify_interesting_markers
from CITEgeist.model.spatial_colocalization import (
    analyze_marker_colocalization,
    discover_hierarchical_profiles,  # NEW: hierarchical discovery
    select_profiles,
)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)


# ============================================================================
# Profile-to-Cell-Type Mapping (same as autodiscovery)
# ============================================================================

CELLTYPE_MARKER_SIGNATURES = {
    "CD8+ T cells": {
        "primary": ["CD3E", "CD8A"],
        "secondary": ["CD45", "GranzymeB"],
    },
    "Macrophages": {
        "primary": ["CD68", "HLA-DR"],
        "secondary": ["CD163", "CD16"],
    },
    "Mixed Immune": {
        "primary": ["CD3E", "HLA-DR"],
        "secondary": ["CD8A"],
    },
    "Epithelial": {
        "primary": ["PanCK"],
        "secondary": ["E-Cadherin", "Vimentin"],
    },
    "Myofibroblasts": {
        "primary": ["alphaSMA", "Vimentin"],
        "secondary": [],
    },
    "Stromal": {
        "primary": ["Vimentin"],
        "secondary": [],
    },
    "Endothelial": {
        "primary": ["CD31"],
        "secondary": ["Vimentin"],
    },
    "B cells": {
        "primary": ["CD20"],
        "secondary": ["CD45RA", "CD45"],
    },
    "Proliferating T": {
        "primary": ["CD3E", "PCNA"],
        "secondary": ["Ki-67"],
    },
    "Vascular Stromal": {
        "primary": ["CD31", "Vimentin"],
        "secondary": [],
    },
    # Additional T cell subtypes for hierarchical detection
    "CD4+ T cells": {
        "primary": ["CD3E", "CD4"],
        "secondary": ["CD45"],
    },
    "T cells": {
        "primary": ["CD3E"],
        "secondary": ["CD45"],
    },
}


def score_profile_celltype_match(
    profile: List[str],
    celltype: str,
    signature: Dict[str, List[str]],
) -> float:
    """Score how well a profile matches a cell type signature."""
    profile_set = set(profile)
    primary = set(signature["primary"])
    secondary = set(signature["secondary"])

    primary_overlap = len(profile_set & primary)
    primary_score = primary_overlap / len(primary) if primary else 0.0

    secondary_overlap = len(profile_set & secondary)
    secondary_score = secondary_overlap / len(secondary) if secondary else 0.0

    combined_score = (2 * primary_score + secondary_score) / 3
    return combined_score


def map_profiles_to_celltypes(
    discovered_profiles: Dict[str, List[str]],
    ground_truth_celltypes: List[str],
) -> Tuple[Dict[str, str], Dict[str, float]]:
    """
    Map discovered profiles to ground truth cell types.

    Note: Takes Dict[str, List[str]] format from hierarchical discovery.
    """
    relevant_signatures = {
        ct: sig
        for ct, sig in CELLTYPE_MARKER_SIGNATURES.items()
        if ct in ground_truth_celltypes
    }

    profile_names = list(discovered_profiles.keys())
    celltype_list = list(relevant_signatures.keys())

    n_profiles = len(profile_names)
    n_celltypes = len(celltype_list)

    score_matrix = np.zeros((n_profiles, n_celltypes))

    for i, profile_name in enumerate(profile_names):
        profile_markers = discovered_profiles[profile_name]
        for j, celltype in enumerate(celltype_list):
            score_matrix[i, j] = score_profile_celltype_match(
                profile_markers, celltype, relevant_signatures[celltype]
            )

    # Greedy assignment
    profile_to_celltype = {}
    match_scores = {}
    used_celltypes = set()

    flat_indices = np.argsort(score_matrix.flatten())[::-1]

    for flat_idx in flat_indices:
        i, j = np.unravel_index(flat_idx, score_matrix.shape)
        profile_name = profile_names[i]
        if profile_name in profile_to_celltype or celltype_list[j] in used_celltypes:
            continue
        if score_matrix[i, j] > 0:
            profile_to_celltype[profile_name] = celltype_list[j]
            match_scores[profile_name] = score_matrix[i, j]
            used_celltypes.add(celltype_list[j])

    return profile_to_celltype, match_scores


def build_cell_profile_dict(
    discovered_profiles: Dict[str, List[str]],
    profile_to_celltype: Dict[str, str],
) -> Dict[str, Dict[str, List[str]]]:
    """Build a cell_profile_dict from discovered profiles."""
    cell_profile_dict = {}

    for profile_name, celltype in profile_to_celltype.items():
        profile_markers = discovered_profiles[profile_name]
        major_markers = profile_markers[:2] if len(profile_markers) >= 2 else profile_markers
        minor_markers = profile_markers[2:] if len(profile_markers) > 2 else []

        cell_profile_dict[celltype] = {
            "Major": list(major_markers),
            "Minor": list(minor_markers),
        }

    return cell_profile_dict


# ============================================================================
# Evaluation Metrics
# ============================================================================


def compute_deconvolution_metrics(
    predicted: pd.DataFrame,
    ground_truth: pd.DataFrame,
    matched_celltypes: List[str],
) -> Dict[str, float]:
    """Compute deconvolution accuracy metrics."""
    common_spots = predicted.index.intersection(ground_truth.index)
    pred_aligned = predicted.loc[common_spots, matched_celltypes]
    gt_aligned = ground_truth.loc[common_spots, matched_celltypes]

    pred_norm = pred_aligned.div(pred_aligned.sum(axis=1) + 1e-10, axis=0)
    gt_norm = gt_aligned.div(gt_aligned.sum(axis=1) + 1e-10, axis=0)

    # Per-spot JSD
    jsds = []
    for spot in common_spots:
        p = pred_norm.loc[spot].values + 1e-10
        q = gt_norm.loc[spot].values + 1e-10
        p /= p.sum()
        q /= q.sum()
        jsds.append(jensenshannon(p, q))

    mean_jsd = np.mean(jsds)

    # Per-celltype correlations
    pearson_by_ct = {}
    spearman_by_ct = {}
    for ct in matched_celltypes:
        pred_vals = pred_norm[ct].values
        gt_vals = gt_norm[ct].values
        if np.std(pred_vals) > 0 and np.std(gt_vals) > 0:
            pearson_by_ct[ct] = pearsonr(pred_vals, gt_vals)[0]
            spearman_by_ct[ct] = spearmanr(pred_vals, gt_vals)[0]

    mean_pearson = np.mean(list(pearson_by_ct.values())) if pearson_by_ct else 0.0
    mean_spearman = np.mean(list(spearman_by_ct.values())) if spearman_by_ct else 0.0

    diff = pred_norm.values - gt_norm.values
    rmse = np.sqrt(np.mean(diff**2))
    mae = np.mean(np.abs(diff))

    return {
        "mean_jsd": mean_jsd,
        "rmse": rmse,
        "mae": mae,
        "mean_pearson": mean_pearson,
        "mean_spearman": mean_spearman,
        "pearson_by_celltype": pearson_by_ct,
        "spearman_by_celltype": spearman_by_ct,
        "n_spots": len(common_spots),
        "n_matched_celltypes": len(matched_celltypes),
    }


# ============================================================================
# Main Benchmark
# ============================================================================


def run_hierarchical_benchmark(
    gex_adata: sc.AnnData,
    protein_adata: sc.AnnData,
    ground_truth: pd.DataFrame,
    region_id: int = 0,
    output_dir: str = ".",
    radius: float = 4.0,
    lambda_reg: float = 1.0,
    alpha_elastic: float = 0.7,
    max_y_change: float = 0.4,
    min_counts: int = 25,
    improvement_threshold: float = 0.05,
    max_depth: int = 5,
) -> Dict:
    """
    Run hierarchical profile discovery benchmark.

    Pipeline:
    1. Module 1: Identify interesting markers
    2. Module 2a: Colocalization analysis
    3. Module 2b: HIERARCHICAL profile discovery
    4. Map profiles to ground truth cell types
    5. Module 3: Cell proportion optimization
    6. Evaluate against ground truth
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    X_protein = (
        protein_adata.X.toarray()
        if hasattr(protein_adata.X, "toarray")
        else protein_adata.X
    )
    coords = protein_adata.obsm["spatial"]
    marker_names = list(protein_adata.var_names)

    logger.info("=" * 70)
    logger.info(f"HIERARCHICAL PROFILE DISCOVERY BENCHMARK: Region {region_id}")
    logger.info("=" * 70)
    logger.info(f"Spots: {X_protein.shape[0]}, Markers: {len(marker_names)}")

    # ========================================================================
    # Module 1: Identify Interesting Markers
    # ========================================================================
    logger.info("-" * 60)
    logger.info("MODULE 1: Identifying Interesting Markers")
    logger.info("-" * 60)

    interest_result = identify_interesting_markers(
        X=X_protein,
        coords=coords,
        marker_names=marker_names,
        kurtosis_threshold=2.0,
        morans_threshold=0.1,
        morans_k=8,
        morans_n_perm=199,
    )

    interesting_markers = interest_result.interesting_markers
    logger.info(f"Found {len(interesting_markers)} interesting markers: {interesting_markers}")

    if len(interesting_markers) < 3:
        logger.error("Not enough interesting markers for profile discovery")
        return {"error": "Not enough interesting markers", "n_interesting": len(interesting_markers)}

    # ========================================================================
    # Module 2a: Colocalization Analysis
    # ========================================================================
    logger.info("-" * 60)
    logger.info("MODULE 2a: Pairwise Colocalization Analysis")
    logger.info("-" * 60)

    coloc_result = analyze_marker_colocalization(
        X=X_protein,
        coords=coords,
        marker_names=marker_names,
        markers_to_analyze=interesting_markers,
        neighbor_k=6,
        n_permutations=199,
    )

    # ========================================================================
    # Module 2b: HIERARCHICAL Profile Discovery
    # ========================================================================
    logger.info("-" * 60)
    logger.info("MODULE 2b: HIERARCHICAL Profile Discovery")
    logger.info("-" * 60)

    hierarchical_result = discover_hierarchical_profiles(
        coloc_result=coloc_result,
        antibody_expression=X_protein,
        marker_names=marker_names,
        coords=coords,
        improvement_threshold=improvement_threshold,
        sharing_ratio=0.5,
        sharing_min_I=0.2,
        max_depth=max_depth,
        verbose=True,
    )

    logger.info(hierarchical_result.summary())
    logger.info(f"Tree depth: {hierarchical_result.tree.get_depth()}")
    logger.info(f"Shared markers: {hierarchical_result.shared_markers}")

    # Get flat profiles (Module 3 compatible)
    discovered_profiles = hierarchical_result.flat_profiles
    logger.info(f"Discovered {len(discovered_profiles)} flat profiles")
    for name, markers in discovered_profiles.items():
        logger.info(f"  {name}: {markers}")

    # ========================================================================
    # Map Profiles to Ground Truth Cell Types
    # ========================================================================
    logger.info("-" * 60)
    logger.info("Mapping Discovered Profiles to Cell Types")
    logger.info("-" * 60)

    gt_celltypes = [c for c in ground_truth.columns if c != "n_cells"]

    profile_to_celltype, match_scores = map_profiles_to_celltypes(
        discovered_profiles=discovered_profiles,
        ground_truth_celltypes=gt_celltypes,
    )

    for profile_name, celltype in profile_to_celltype.items():
        logger.info(
            f"  {profile_name} ({discovered_profiles[profile_name]}) -> {celltype} "
            f"(score: {match_scores[profile_name]:.2f})"
        )

    matched_celltypes = list(profile_to_celltype.values())
    unmatched_gt = [ct for ct in gt_celltypes if ct not in matched_celltypes]
    if unmatched_gt:
        logger.warning(f"Unmatched ground truth cell types: {unmatched_gt}")

    # Build cell profile dict
    cell_profile_dict = build_cell_profile_dict(discovered_profiles, profile_to_celltype)

    logger.info(f"Built cell_profile_dict with {len(cell_profile_dict)} cell types")

    # ========================================================================
    # Module 3: Cell Proportion Optimization
    # ========================================================================
    logger.info("-" * 60)
    logger.info("MODULE 3: Cell Proportion Optimization")
    logger.info("-" * 60)

    from CITEgeist.model.citegeist_model import CitegeistModel

    model = CitegeistModel(
        sample_name=f"hierarchical_region_{region_id}",
        output_folder=str(output_dir),
        simulation=True,
        gene_expression_adata=gex_adata,
        antibody_capture_adata=protein_adata,
    )

    model.load_cell_profile_dict(cell_profile_dict)
    model.filter_gex(min_counts=min_counts)
    model.preprocess_gex(target_sum=10000)
    model.preprocess_antibody()

    global_props, finetuned_props = model.run_cell_proportion_model(
        radius=radius,
        lambda_reg=lambda_reg,
        alpha=alpha_elastic,
        max_y_change=max_y_change,
        validation_warn_only=True,
    )

    if finetuned_props is not None:
        predicted_proportions = pd.DataFrame(
            finetuned_props,
            index=model.gene_expression_adata.obs_names,
            columns=list(cell_profile_dict.keys()),
        )
    else:
        predicted_proportions = pd.DataFrame(
            global_props,
            index=model.gene_expression_adata.obs_names,
            columns=list(cell_profile_dict.keys()),
        )

    # ========================================================================
    # Evaluation
    # ========================================================================
    logger.info("-" * 60)
    logger.info("EVALUATION: Comparing to Ground Truth")
    logger.info("-" * 60)

    metrics = compute_deconvolution_metrics(
        predicted=predicted_proportions,
        ground_truth=ground_truth,
        matched_celltypes=matched_celltypes,
    )

    logger.info(f"Mean JSD: {metrics['mean_jsd']:.4f}")
    logger.info(f"RMSE: {metrics['rmse']:.4f}")
    logger.info(f"MAE: {metrics['mae']:.4f}")
    logger.info(f"Mean Pearson: {metrics['mean_pearson']:.4f}")
    logger.info(f"Mean Spearman: {metrics['mean_spearman']:.4f}")

    # ========================================================================
    # Save Results
    # ========================================================================
    results = {
        "region_id": region_id,
        "mode": "hierarchical_discovery",
        "n_spots": X_protein.shape[0],
        "n_markers": len(marker_names),
        "n_interesting_markers": len(interesting_markers),
        "interesting_markers": interesting_markers,
        "n_discovered_profiles": len(discovered_profiles),
        "discovered_profiles": {k: list(v) for k, v in discovered_profiles.items()},
        "profile_to_celltype": profile_to_celltype,
        "match_scores": match_scores,
        "matched_celltypes": matched_celltypes,
        "unmatched_gt_celltypes": unmatched_gt,
        "cell_profile_dict": cell_profile_dict,
        # Hierarchical-specific results
        "tree_depth": hierarchical_result.tree.get_depth(),
        "shared_markers": {k: list(v) for k, v in hierarchical_result.shared_markers.items()},
        "reconstruction_error": hierarchical_result.reconstruction_error,
        "metrics": metrics,
        "parameters": {
            "radius": radius,
            "lambda_reg": lambda_reg,
            "alpha_elastic": alpha_elastic,
            "max_y_change": max_y_change,
            "min_counts": min_counts,
            "improvement_threshold": improvement_threshold,
            "max_depth": max_depth,
        },
    }

    # Save results
    with open(output_dir / f"region_{region_id}_hierarchical_results.json", "w") as f:
        def convert_numpy(obj):
            if isinstance(obj, np.ndarray):
                return obj.tolist()
            elif isinstance(obj, (np.int64, np.int32)):
                return int(obj)
            elif isinstance(obj, (np.float64, np.float32)):
                return float(obj)
            elif isinstance(obj, dict):
                return {k: convert_numpy(v) for k, v in obj.items()}
            elif isinstance(obj, list):
                return [convert_numpy(v) for v in obj]
            return obj

        json.dump(convert_numpy(results), f, indent=2)

    predicted_proportions.to_csv(
        output_dir / f"region_{region_id}_hierarchical_proportions.csv"
    )

    logger.info("=" * 70)
    logger.info("HIERARCHICAL BENCHMARK COMPLETE")
    logger.info("=" * 70)

    return results


def main():
    parser = argparse.ArgumentParser(
        description="Run CITEgeist hierarchical profile discovery benchmark (Mode C)"
    )
    parser.add_argument(
        "--region-id", type=int, default=0, help="Region ID to process"
    )
    parser.add_argument(
        "--input-dir",
        type=str,
        default=str(REPO_ROOT / "Benchmarking/xenium_pseudovisium/data_granular_gt"),
        help="Input directory with h5ad_objects/ and ground_truth/",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default=str(REPO_ROOT / "Benchmarking/xenium_benchmarking/CITEgeist/output_hierarchical"),
        help="Output directory for results",
    )
    parser.add_argument("--radius", type=float, default=4.0, help="Spatial neighbor radius")
    parser.add_argument("--lambda-reg", type=float, default=1.0, help="Regularization lambda")
    parser.add_argument("--alpha-elastic", type=float, default=0.7, help="Elastic net alpha")
    parser.add_argument("--max-y-change", type=float, default=0.4, help="Max Y change constraint")
    parser.add_argument("--min-counts", type=int, default=25, help="Min counts filter")
    parser.add_argument("--improvement-threshold", type=float, default=0.05, help="Reconstruction improvement threshold for tree cutting")
    parser.add_argument("--max-depth", type=int, default=5, help="Maximum tree depth")

    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)

    # Load data
    gex_path = input_dir / "h5ad_objects" / f"Xenium_region_{args.region_id}_GEX.h5ad"
    protein_path = input_dir / "h5ad_objects" / f"Xenium_region_{args.region_id}_CITE.h5ad"
    gt_path = input_dir / "ground_truth" / f"Xenium_region_{args.region_id}_prop.csv"

    for path in [gex_path, protein_path, gt_path]:
        if not path.exists():
            logger.error(f"File not found: {path}")
            sys.exit(1)

    logger.info(f"Loading data for region {args.region_id}...")
    gex_adata = sc.read_h5ad(gex_path)
    protein_adata = sc.read_h5ad(protein_path)
    ground_truth = pd.read_csv(gt_path, index_col=0)

    results = run_hierarchical_benchmark(
        gex_adata=gex_adata,
        protein_adata=protein_adata,
        ground_truth=ground_truth,
        region_id=args.region_id,
        output_dir=str(output_dir),
        radius=args.radius,
        lambda_reg=args.lambda_reg,
        alpha_elastic=args.alpha_elastic,
        max_y_change=args.max_y_change,
        min_counts=args.min_counts,
        improvement_threshold=args.improvement_threshold,
        max_depth=args.max_depth,
    )


if __name__ == "__main__":
    main()
