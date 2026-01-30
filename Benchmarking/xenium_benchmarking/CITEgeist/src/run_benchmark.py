#!/usr/bin/env python
"""
Unified CITEgeist benchmark runner for Xenium pseudo-Visium data.

Two modes:
- manual: Uses hand-crafted achievable-7 profiles (fair benchmark baseline)
- hierarchical: Auto-discovers profiles with hierarchy learning (flagship feature)

Usage:
    python run_benchmark.py --region-id 0 --mode manual       # Achievable-7 profiles
    python run_benchmark.py --region-id 0 --mode hierarchical # Auto-discover with hierarchy
"""

import argparse
import json
import logging
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.spatial.distance import jensenshannon
from scipy.stats import pearsonr, spearmanr

# Add CITEgeist to path
REPO_ROOT = Path(__file__).parent.parent.parent.parent.parent
sys.path.insert(0, str(REPO_ROOT))

# Import shared benchmark constants
BENCHMARK_ROOT = Path(__file__).parent.parent.parent
sys.path.insert(0, str(BENCHMARK_ROOT))
from benchmark_constants import (
    ACHIEVABLE_7_CELL_PROFILE_DICT,
    ACHIEVABLE_7_MARKER_SIGNATURES,
    GT_TO_ACHIEVABLE_7_MAPPING,
)

from CITEgeist.model.citegeist_model import CitegeistModel
from CITEgeist.model.marker_interest import identify_interesting_markers
from CITEgeist.model.spatial_colocalization import (
    analyze_marker_colocalization,
    discover_hierarchical_profiles,
    rescue_singletons,
    select_profiles,
)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)


# =============================================================================
# HIERARCHICAL PROFILE MAPPING
# =============================================================================


def compress_to_achievable_7(
    discovered_profiles: Dict[str, List[str]],
) -> Dict[str, str]:
    """
    Map N discovered profiles → 7 achievable types via Jaccard similarity.

    Multiple discovered profiles may map to same achievable type.
    This enables fair benchmarking while preserving discovery granularity.

    Args:
        discovered_profiles: Dict mapping profile names to marker lists

    Returns:
        Dict mapping profile names to achievable-7 cell type names
    """
    profile_mapping = {}

    for profile_name, markers in discovered_profiles.items():
        marker_set = set(markers)
        best_match = None
        best_jaccard = -1

        for celltype, signature in ACHIEVABLE_7_MARKER_SIGNATURES.items():
            intersection = len(marker_set & signature)
            union = len(marker_set | signature)
            jaccard = intersection / union if union > 0 else 0

            if jaccard > best_jaccard:
                best_jaccard = jaccard
                best_match = celltype

        profile_mapping[profile_name] = best_match if best_match else "Unknown"

    return profile_mapping


def build_cell_profile_dict_from_hierarchical(
    discovered_profiles: Dict[str, List[str]],
    profile_mapping: Dict[str, str],
) -> Dict[str, Dict[str, List[str]]]:
    """
    Convert discovered profiles to cell_profile_dict format for Module 3.

    Uses the achievable-7 cell type names from the mapping.
    """
    cell_profile_dict = {}

    for profile_name, celltype in profile_mapping.items():
        markers = discovered_profiles[profile_name]
        major = markers[:2] if len(markers) >= 2 else markers
        minor = markers[2:] if len(markers) > 2 else []

        # If multiple profiles map to same celltype, merge markers
        if celltype in cell_profile_dict:
            existing_major = set(cell_profile_dict[celltype]["Major"])
            existing_minor = set(cell_profile_dict[celltype]["Minor"])
            cell_profile_dict[celltype]["Major"] = list(existing_major | set(major))
            cell_profile_dict[celltype]["Minor"] = list(existing_minor | set(minor))
        else:
            cell_profile_dict[celltype] = {
                "Major": list(major),
                "Minor": list(minor),
            }

    return cell_profile_dict


# =============================================================================
# MAIN BENCHMARK FUNCTIONS
# =============================================================================


def run_manual_benchmark(
    region_id: int,
    gex_adata: sc.AnnData,
    protein_adata: sc.AnnData,
    output_dir: Path,
    radius: float = 4.0,
    lambda_reg: float = 1.0,
    alpha_elastic: float = 0.7,
    max_y_change: float = 0.4,
    min_counts: int = 25,
    run_gex: bool = False,
) -> Dict[str, Any]:
    """
    Run CITEgeist with manual achievable-7 profiles.

    This is the fair benchmark baseline for method comparison.
    """
    logger.info("=" * 70)
    logger.info(f"MANUAL MODE: Region {region_id} with ACHIEVABLE-7 profiles")
    logger.info("=" * 70)

    sample_name = f"Xenium_region_{region_id}"

    model = CitegeistModel(
        sample_name=sample_name,
        output_folder=str(output_dir),
        simulation=True,
        gene_expression_adata=gex_adata,
        antibody_capture_adata=protein_adata,
    )

    model.load_cell_profile_dict(ACHIEVABLE_7_CELL_PROFILE_DICT)

    logger.info("Preprocessing...")
    model.filter_gex(min_counts=min_counts)
    model.preprocess_gex(target_sum=10000)
    model.preprocess_antibody()

    logger.info("Running cell proportion optimization...")
    start_time = time.time()

    global_props, finetuned_props = model.run_cell_proportion_model(
        radius=radius,
        lambda_reg=lambda_reg,
        alpha=alpha_elastic,
        max_y_change=max_y_change,
        validation_warn_only=True,
    )

    prop_time = time.time() - start_time
    logger.info(f"Cell proportion optimization completed in {prop_time:.1f}s")

    # Save results
    result_dir = output_dir / sample_name
    result_dir.mkdir(parents=True, exist_ok=True)

    if finetuned_props is not None:
        props_df = pd.DataFrame(
            finetuned_props,
            index=model.antibody_capture_adata.obs_names,
            columns=list(ACHIEVABLE_7_CELL_PROFILE_DICT.keys()),
        )
        props_df.to_csv(result_dir / f"{sample_name}_deconv_predictions.csv")

    # Run GEX if requested
    gex_time = None
    if run_gex:
        logger.info("Running gene expression deconvolution...")
        gex_start = time.time()
        try:
            model.run_cell_expression_pass1(
                radius=radius,
                alpha=0.5,
                checkpoint_interval=100,
                output_dir=str(output_dir / "checkpoints"),
                rerun=True,
            )
            gex_time = time.time() - gex_start
            logger.info(f"GEX deconvolution completed in {gex_time:.1f}s")
        except Exception as e:
            logger.error(f"GEX deconvolution failed: {e}")
            gex_time = -1.0

    results = {
        "region_id": region_id,
        "mode": "manual",
        "n_spots": gex_adata.shape[0],
        "n_cell_types": len(ACHIEVABLE_7_CELL_PROFILE_DICT),
        "cell_types": list(ACHIEVABLE_7_CELL_PROFILE_DICT.keys()),
        "runtime_proportions_sec": prop_time,
        "runtime_gex_sec": gex_time,
        "output_dir": str(result_dir),
    }

    with open(result_dir / "run_stats.json", "w") as f:
        json.dump(results, f, indent=2)

    logger.info(f"Results saved to {result_dir}")
    return results


def run_hierarchical_benchmark(
    region_id: int,
    gex_adata: sc.AnnData,
    protein_adata: sc.AnnData,
    output_dir: Path,
    radius: float = 4.0,
    lambda_reg: float = 1.0,
    alpha_elastic: float = 0.7,
    max_y_change: float = 0.4,
    min_counts: int = 25,
    run_gex: bool = False,
    # Hierarchical-specific parameters
    fdr_alpha: float = 0.05,
    top_k: int = 3,
    improvement_threshold: float = 0.05,
    max_depth: int = 5,
) -> Dict[str, Any]:
    """
    Run CITEgeist with hierarchical profile autodiscovery.

    This is the flagship feature showcasing data-driven discovery.
    """
    logger.info("=" * 70)
    logger.info(f"HIERARCHICAL MODE: Region {region_id} with autodiscovery")
    logger.info("=" * 70)

    sample_name = f"Xenium_region_{region_id}"

    X_protein = (
        protein_adata.X.toarray()
        if hasattr(protein_adata.X, "toarray")
        else protein_adata.X
    )
    coords = protein_adata.obsm["spatial"]
    marker_names = list(protein_adata.var_names)

    logger.info(f"Spots: {X_protein.shape[0]}, Markers: {len(marker_names)}")

    # =========================================================================
    # Module 1: Identify Interesting Markers
    # =========================================================================
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
    logger.info(f"Found {len(interesting_markers)} interesting markers")

    if len(interesting_markers) < 3:
        logger.warning("Not enough interesting markers, using all markers")
        interesting_markers = marker_names

    # =========================================================================
    # Module 2a: Colocalization Analysis
    # =========================================================================
    logger.info("-" * 60)
    logger.info("MODULE 2a: Pairwise Colocalization Analysis")
    logger.info("-" * 60)

    coloc_result = analyze_marker_colocalization(
        X=X_protein,
        coords=coords,
        marker_names=marker_names,
        markers_to_analyze=interesting_markers,
        neighbor_k=6,
        n_permutations=999,
    )

    # =========================================================================
    # Module 2b: HIERARCHICAL Profile Discovery
    # =========================================================================
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
        fdr_alpha=fdr_alpha,
        top_k=top_k,
        verbose=True,
    )

    logger.info(hierarchical_result.summary())
    logger.info(f"Tree depth: {hierarchical_result.tree.get_depth()}")
    logger.info(f"Shared markers: {hierarchical_result.shared_markers}")

    discovered_profiles = hierarchical_result.flat_profiles
    logger.info(f"Discovered {len(discovered_profiles)} profiles:")
    for name, markers in discovered_profiles.items():
        logger.info(f"  {name}: {markers}")

    # =========================================================================
    # Singleton Rescue: Filter noise singletons by unique spatial coverage
    # =========================================================================
    if interest_result.signal_masks is not None and interest_result.signal_mask_marker_names is not None:
        logger.info("-" * 60)
        logger.info("SINGLETON RESCUE: Filtering by unique spatial coverage")
        logger.info("-" * 60)

        # Convert dict to list-of-lists for rescue function
        profile_names = list(discovered_profiles.keys())
        profile_marker_lists = [discovered_profiles[name] for name in profile_names]

        rescued_marker_lists = rescue_singletons(
            profiles=profile_marker_lists,
            signal_masks=interest_result.signal_masks,
            signal_mask_marker_names=interest_result.signal_mask_marker_names,
            min_unique_coverage=0.3,
            min_signal_fraction=0.05,
            verbose=True,
        )

        # Rebuild dict preserving names for kept profiles
        rescued_profiles = {}
        for markers in rescued_marker_lists:
            # Find original name
            matched = False
            for name, orig_markers in discovered_profiles.items():
                if orig_markers == markers and name not in rescued_profiles:
                    rescued_profiles[name] = markers
                    matched = True
                    break
            if not matched:
                # Fallback: generate name
                rescued_profiles[f"Profile_{len(rescued_profiles)}"] = markers

        logger.info(f"After rescue: {len(rescued_profiles)} profiles (was {len(discovered_profiles)})")
        discovered_profiles = rescued_profiles
    else:
        logger.warning("No GMM signal masks available, skipping singleton rescue")

    # =========================================================================
    # Map to Achievable-7 for Fair Comparison
    # =========================================================================
    logger.info("-" * 60)
    logger.info("Mapping profiles to achievable-7 cell types")
    logger.info("-" * 60)

    profile_mapping = compress_to_achievable_7(discovered_profiles)
    for profile, celltype in profile_mapping.items():
        logger.info(f"  {profile} → {celltype}")

    cell_profile_dict = build_cell_profile_dict_from_hierarchical(
        discovered_profiles, profile_mapping
    )

    logger.info(f"Mapped to {len(cell_profile_dict)} cell types for Module 3")

    # =========================================================================
    # Module 3: Cell Proportion Optimization
    # =========================================================================
    logger.info("-" * 60)
    logger.info("MODULE 3: Cell Proportion Optimization")
    logger.info("-" * 60)

    model = CitegeistModel(
        sample_name=sample_name,
        output_folder=str(output_dir),
        simulation=True,
        gene_expression_adata=gex_adata,
        antibody_capture_adata=protein_adata,
    )

    model.load_cell_profile_dict(cell_profile_dict)
    model.filter_gex(min_counts=min_counts)
    model.preprocess_gex(target_sum=10000)
    model.preprocess_antibody()

    start_time = time.time()
    global_props, finetuned_props = model.run_cell_proportion_model(
        radius=radius,
        lambda_reg=lambda_reg,
        alpha=alpha_elastic,
        max_y_change=max_y_change,
        validation_warn_only=True,
    )
    prop_time = time.time() - start_time
    logger.info(f"Cell proportion optimization completed in {prop_time:.1f}s")

    # Save results
    result_dir = output_dir / sample_name
    result_dir.mkdir(parents=True, exist_ok=True)

    if finetuned_props is not None:
        props_df = pd.DataFrame(
            finetuned_props,
            index=model.antibody_capture_adata.obs_names,
            columns=list(cell_profile_dict.keys()),
        )
        props_df.to_csv(result_dir / f"{sample_name}_deconv_predictions.csv")

    # Run GEX if requested
    gex_time = None
    if run_gex:
        logger.info("Running gene expression deconvolution...")
        gex_start = time.time()
        try:
            model.run_cell_expression_pass1(
                radius=radius,
                alpha=0.5,
                checkpoint_interval=100,
                output_dir=str(output_dir / "checkpoints"),
                rerun=True,
            )
            gex_time = time.time() - gex_start
            logger.info(f"GEX deconvolution completed in {gex_time:.1f}s")
        except Exception as e:
            logger.error(f"GEX deconvolution failed: {e}")
            gex_time = -1.0

    # =========================================================================
    # Save comprehensive results
    # =========================================================================
    results = {
        "region_id": region_id,
        "mode": "hierarchical",
        "n_spots": gex_adata.shape[0],
        "n_interesting_markers": len(interesting_markers),
        "interesting_markers": interesting_markers,
        # Hierarchical discovery results
        "n_discovered_profiles": len(discovered_profiles),
        "discovered_profiles": {k: list(v) for k, v in discovered_profiles.items()},
        "profile_mapping": profile_mapping,
        "tree_depth": hierarchical_result.tree.get_depth(),
        "shared_markers": {k: list(v) for k, v in hierarchical_result.shared_markers.items()},
        "reconstruction_error": hierarchical_result.reconstruction_error,
        # Final cell types used
        "n_cell_types": len(cell_profile_dict),
        "cell_types": list(cell_profile_dict.keys()),
        "cell_profile_dict": cell_profile_dict,
        # Timing
        "runtime_proportions_sec": prop_time,
        "runtime_gex_sec": gex_time,
        "output_dir": str(result_dir),
        # Parameters
        "parameters": {
            "fdr_alpha": fdr_alpha,
            "top_k": top_k,
            "improvement_threshold": improvement_threshold,
            "max_depth": max_depth,
        },
    }

    with open(result_dir / "run_stats.json", "w") as f:
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

    logger.info(f"Results saved to {result_dir}")
    return results


# =============================================================================
# CLI
# =============================================================================


def main():
    parser = argparse.ArgumentParser(
        description="Unified CITEgeist benchmark for Xenium pseudo-Visium data"
    )

    # Required arguments
    parser.add_argument(
        "--region-id",
        type=int,
        required=True,
        help="Region ID to process (0-13 for all regions)",
    )
    parser.add_argument(
        "--mode",
        choices=["manual", "hierarchical"],
        default="manual",
        help="Benchmark mode: 'manual' (achievable-7) or 'hierarchical' (autodiscovery)",
    )

    # I/O arguments
    parser.add_argument(
        "--input-dir",
        type=str,
        default=str(REPO_ROOT / "Benchmarking/xenium_pseudovisium/data_granular_gt"),
        help="Input directory with h5ad_objects/",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default=None,
        help="Output directory (auto-set to output/manual or output/hierarchical if not specified)",
    )

    # Optimization parameters
    parser.add_argument("--radius", type=float, default=4.0, help="Spatial neighbor radius")
    parser.add_argument("--lambda-reg", type=float, default=1.0, help="Regularization lambda")
    parser.add_argument("--alpha-elastic", type=float, default=0.7, help="Elastic net alpha")
    parser.add_argument("--max-y-change", type=float, default=0.4, help="Max Y change")
    parser.add_argument("--min-counts", type=int, default=25, help="Min counts filter")
    parser.add_argument("--run-gex", action="store_true", help="Run GEX deconvolution")
    parser.add_argument(
        "--no-unknown",
        action="store_true",
        help="Exclude unknown/unassigned profile from cell type dictionary",
    )

    # Hierarchical-specific parameters
    parser.add_argument(
        "--fdr-alpha",
        type=float,
        default=0.05,
        help="FDR threshold for colocalization (hierarchical mode)",
    )
    parser.add_argument(
        "--top-k",
        type=int,
        default=3,
        help="Mutual top-k for edge sparsification (hierarchical mode)",
    )
    parser.add_argument(
        "--improvement-threshold",
        type=float,
        default=0.05,
        help="Reconstruction improvement threshold for tree cutting (hierarchical mode)",
    )
    parser.add_argument(
        "--max-depth",
        type=int,
        default=5,
        help="Maximum tree depth (hierarchical mode)",
    )

    args = parser.parse_args()

    # Set output directory based on mode if not specified
    if args.output_dir is None:
        base_dir = REPO_ROOT / "Benchmarking/xenium_benchmarking/CITEgeist/output"
        args.output_dir = str(base_dir / args.mode)

    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load data
    gex_path = input_dir / "h5ad_objects" / f"Xenium_region_{args.region_id}_GEX.h5ad"
    protein_path = input_dir / "h5ad_objects" / f"Xenium_region_{args.region_id}_CITE.h5ad"

    for path in [gex_path, protein_path]:
        if not path.exists():
            logger.error(f"File not found: {path}")
            sys.exit(1)

    logger.info(f"Loading data for region {args.region_id}...")
    gex_adata = sc.read_h5ad(gex_path)
    protein_adata = sc.read_h5ad(protein_path)

    # Run benchmark
    if args.mode == "manual":
        results = run_manual_benchmark(
            region_id=args.region_id,
            gex_adata=gex_adata,
            protein_adata=protein_adata,
            output_dir=output_dir,
            radius=args.radius,
            lambda_reg=args.lambda_reg,
            alpha_elastic=args.alpha_elastic,
            max_y_change=args.max_y_change,
            min_counts=args.min_counts,
            run_gex=args.run_gex,
        )
    else:  # hierarchical
        results = run_hierarchical_benchmark(
            region_id=args.region_id,
            gex_adata=gex_adata,
            protein_adata=protein_adata,
            output_dir=output_dir,
            radius=args.radius,
            lambda_reg=args.lambda_reg,
            alpha_elastic=args.alpha_elastic,
            max_y_change=args.max_y_change,
            min_counts=args.min_counts,
            run_gex=args.run_gex,
            fdr_alpha=args.fdr_alpha,
            top_k=args.top_k,
            improvement_threshold=args.improvement_threshold,
            max_depth=args.max_depth,
        )

    print(f"\nCompleted region {args.region_id} ({args.mode} mode)")
    print(f"  Spots: {results['n_spots']}")
    print(f"  Cell types: {results['n_cell_types']}")
    print(f"  Runtime: {results['runtime_proportions_sec']:.1f}s")
    print(f"  Output: {results['output_dir']}")


if __name__ == "__main__":
    main()
