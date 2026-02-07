"""
Test specificity-weighted profile discovery on simulated and Xenium datasets.

This script compares profile discovery with and without specificity weighting
to evaluate the impact on handling promiscuous markers like CD45.

Usage:
    python tests/test_specificity_weighting.py --dataset simulated_high_seg
    python tests/test_specificity_weighting.py --dataset simulated_mixed
    python tests/test_specificity_weighting.py --dataset xenium
"""

import argparse
import json
import logging
import sys
from pathlib import Path
from typing import Dict, List, Tuple, Any

import numpy as np
import pandas as pd
import scanpy as sc

# Add repo root to path
REPO_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(REPO_ROOT))

from CITEgeist.model.marker_interest import identify_interesting_markers
from CITEgeist.model.spatial_colocalization import (
    analyze_marker_colocalization,
    discover_profiles,
    select_profiles,
)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)


def load_simulated_data(replicate: int = 0, variant: str = "high_seg") -> Tuple[np.ndarray, np.ndarray, List[str], pd.DataFrame]:
    """
    Load simulated Wu data.

    Args:
        replicate: Replicate index (0-4)
        variant: 'high_seg' or 'mixed'

    Returns:
        Tuple of (X, coords, marker_names, ground_truth_proportions)
    """
    base_path = REPO_ROOT / "replicates" / variant / "h5ad_objects"
    prop_path = REPO_ROOT / "replicates" / variant / "ST_sim"

    # Load CITE data
    cite_path = base_path / f"Wu_rep_{replicate}_CITE.h5ad"
    if not cite_path.exists():
        raise FileNotFoundError(f"Simulated CITE data not found: {cite_path}")

    adata = sc.read_h5ad(cite_path)
    X = adata.X if isinstance(adata.X, np.ndarray) else adata.X.toarray()
    coords = adata.obsm.get('spatial', adata.obs[['x', 'y']].values if 'x' in adata.obs else None)
    marker_names = list(adata.var_names)

    # Load ground truth proportions
    prop_file = prop_path / f"Wu_ST_{replicate}_prop.csv"
    if prop_file.exists():
        gt_props = pd.read_csv(prop_file, index_col=0)
    else:
        gt_props = None

    logger.info(f"Loaded simulated {variant} replicate {replicate}: {X.shape[0]} spots, {X.shape[1]} markers")

    return X, coords, marker_names, gt_props


def load_xenium_data(region_id: int = 0) -> Tuple[np.ndarray, np.ndarray, List[str], pd.DataFrame]:
    """
    Load Xenium pseudovisium data.

    Args:
        region_id: Region index (0-4)

    Returns:
        Tuple of (X, coords, marker_names, ground_truth_proportions)
    """
    base_path = REPO_ROOT / "Benchmarking" / "xenium_pseudovisium" / "data_granular_gt"
    h5ad_path = base_path / "h5ad_objects"

    # Load CITE data
    cite_path = h5ad_path / f"Xenium_region_{region_id}_CITE.h5ad"
    if not cite_path.exists():
        raise FileNotFoundError(f"Xenium CITE data not found: {cite_path}")

    adata = sc.read_h5ad(cite_path)
    X = adata.X if isinstance(adata.X, np.ndarray) else adata.X.toarray()
    coords = adata.obsm.get('spatial', adata.obs[['spot_x', 'spot_y']].values if 'spot_x' in adata.obs else None)
    marker_names = list(adata.var_names)

    # Load ground truth proportions
    gt_path = base_path / "ground_truth"
    prop_file = gt_path / f"Xenium_region_{region_id}_prop.csv"
    if prop_file.exists():
        gt_props = pd.read_csv(prop_file, index_col=0)
    else:
        gt_props = None

    logger.info(f"Loaded Xenium region {region_id}: {X.shape[0]} spots, {X.shape[1]} markers")

    return X, coords, marker_names, gt_props


def run_profile_discovery(
    X: np.ndarray,
    coords: np.ndarray,
    marker_names: List[str],
    use_specificity: bool = True,
    fdr_alpha: float = 0.05,
    top_k: int = 3,
) -> Dict[str, Any]:
    """
    Run the full profile discovery pipeline.

    Args:
        X: Expression matrix (n_spots, n_markers)
        coords: Spatial coordinates (n_spots, 2)
        marker_names: List of marker names
        use_specificity: Whether to use specificity weighting
        fdr_alpha: FDR threshold
        top_k: Number of top partners in mutual top-k

    Returns:
        Dictionary with discovery results and metrics
    """
    logger.info(f"Running profile discovery (specificity_weighted={use_specificity})...")

    # Module 1: Identify interesting markers
    interest_result = identify_interesting_markers(
        X, coords, marker_names,
        morans_threshold=0.1,
        morans_n_perm=199,
        verbose=True,
    )
    interesting_markers = interest_result.interesting_markers
    logger.info(f"Module 1: {len(interesting_markers)} interesting markers identified")

    # Module 2a: Pairwise colocalization
    coloc_result = analyze_marker_colocalization(
        X, coords, marker_names,
        markers_to_analyze=interesting_markers,
        n_permutations=199,
        verbose=True,
    )

    # Log specificity info
    if coloc_result.marker_specificity:
        spec = coloc_result.marker_specificity
        promiscuous = [(m, s) for m, s in spec.items() if s < 0.3]
        specific = [(m, s) for m, s in spec.items() if s > 0.6]
        logger.info(f"Marker specificity: {len(promiscuous)} promiscuous, {len(specific)} highly specific")
        if promiscuous:
            logger.info(f"  Promiscuous: {promiscuous[:5]}")
        if specific:
            logger.info(f"  Specific: {specific[:5]}")

    # Module 2b: Profile discovery
    # If NOT using specificity, temporarily clear it from the result
    if not use_specificity:
        original_specificity = coloc_result.marker_specificity
        coloc_result.marker_specificity = {}

    discovery_result = discover_profiles(
        coloc_result,
        fdr_alpha=fdr_alpha,
        top_k=top_k,
        verbose=True,
    )

    # Restore specificity if we cleared it
    if not use_specificity:
        coloc_result.marker_specificity = original_specificity

    logger.info(f"Module 2b: {len(discovery_result.profiles)} profiles discovered")
    logger.info(f"  Modularity: {discovery_result.modularity:.3f}")

    # Module 2c: Profile selection
    # Use 0.5% marginal gain (matching sbatch_CITEgeistArray.sh) and max 15 profiles
    # to prevent timeout with many singleton profiles
    selection_result = select_profiles(
        X, coords, marker_names,
        discovery_result.profiles,
        interesting_markers,
        coloc_result,
        min_marginal_gain=0.005,  # 0.5% matches working benchmarks
        max_profiles=15,           # Limit iterations with many singletons
        verbose=True,
    )

    logger.info(f"Module 2c: {selection_result.optimal_n} profiles selected")
    logger.info(f"  Variance explained: {selection_result.variance_explained[-1] if len(selection_result.variance_explained) > 0 else 0:.1%}")

    return {
        "use_specificity": use_specificity,
        "n_interesting_markers": len(interesting_markers),
        "interesting_markers": interesting_markers,
        "n_discovered_profiles": len(discovery_result.profiles),
        "discovered_profiles": discovery_result.profiles,
        "modularity": discovery_result.modularity,
        "n_significant_edges": discovery_result.n_significant_edges,
        "singletons": discovery_result.singletons,
        "n_selected_profiles": selection_result.optimal_n,
        "selected_profiles": selection_result.selected_profiles,
        "variance_explained": float(selection_result.variance_explained[-1]) if len(selection_result.variance_explained) > 0 else 0,
        "marginal_gains": selection_result.marginal_gains.tolist(),
        "stopping_reason": selection_result.stopping_reason,
        "marker_specificity": coloc_result.marker_specificity,
    }


def compare_results(with_spec: Dict, without_spec: Dict) -> Dict[str, Any]:
    """Compare results with and without specificity weighting."""
    comparison = {
        "n_discovered_profiles": {
            "with_specificity": with_spec["n_discovered_profiles"],
            "without_specificity": without_spec["n_discovered_profiles"],
            "difference": with_spec["n_discovered_profiles"] - without_spec["n_discovered_profiles"],
        },
        "n_selected_profiles": {
            "with_specificity": with_spec["n_selected_profiles"],
            "without_specificity": without_spec["n_selected_profiles"],
            "difference": with_spec["n_selected_profiles"] - without_spec["n_selected_profiles"],
        },
        "modularity": {
            "with_specificity": with_spec["modularity"],
            "without_specificity": without_spec["modularity"],
            "difference": with_spec["modularity"] - without_spec["modularity"],
        },
        "variance_explained": {
            "with_specificity": with_spec["variance_explained"],
            "without_specificity": without_spec["variance_explained"],
            "difference": with_spec["variance_explained"] - without_spec["variance_explained"],
        },
        "n_singletons": {
            "with_specificity": len(with_spec["singletons"]),
            "without_specificity": len(without_spec["singletons"]),
            "difference": len(with_spec["singletons"]) - len(without_spec["singletons"]),
        },
    }
    return comparison


def main():
    parser = argparse.ArgumentParser(description="Test specificity weighting on profile discovery")
    parser.add_argument(
        "--dataset",
        type=str,
        choices=["simulated_high_seg", "simulated_mixed", "xenium"],
        default="simulated_high_seg",
        help="Dataset to test on",
    )
    parser.add_argument(
        "--replicate",
        type=int,
        default=0,
        help="Replicate/region index (0-4)",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default=None,
        help="Output directory for results",
    )
    args = parser.parse_args()

    # Load data
    if args.dataset == "simulated_high_seg":
        X, coords, marker_names, gt_props = load_simulated_data(args.replicate, "high_seg")
        dataset_name = f"simulated_high_seg_rep{args.replicate}"
    elif args.dataset == "simulated_mixed":
        X, coords, marker_names, gt_props = load_simulated_data(args.replicate, "mixed")
        dataset_name = f"simulated_mixed_rep{args.replicate}"
    else:  # xenium
        X, coords, marker_names, gt_props = load_xenium_data(args.replicate)
        dataset_name = f"xenium_region{args.replicate}"

    logger.info(f"\n{'='*60}")
    logger.info(f"Testing specificity weighting on {dataset_name}")
    logger.info(f"{'='*60}")

    # Run with specificity weighting
    logger.info("\n--- WITH Specificity Weighting ---")
    results_with_spec = run_profile_discovery(
        X, coords, marker_names, use_specificity=True
    )

    # Run without specificity weighting
    logger.info("\n--- WITHOUT Specificity Weighting ---")
    results_without_spec = run_profile_discovery(
        X, coords, marker_names, use_specificity=False
    )

    # Compare results
    comparison = compare_results(results_with_spec, results_without_spec)

    logger.info("\n" + "="*60)
    logger.info("COMPARISON SUMMARY")
    logger.info("="*60)
    for metric, values in comparison.items():
        logger.info(f"{metric}:")
        logger.info(f"  With specificity:    {values['with_specificity']:.4f}" if isinstance(values['with_specificity'], float) else f"  With specificity:    {values['with_specificity']}")
        logger.info(f"  Without specificity: {values['without_specificity']:.4f}" if isinstance(values['without_specificity'], float) else f"  Without specificity: {values['without_specificity']}")
        logger.info(f"  Difference:          {values['difference']:+.4f}" if isinstance(values['difference'], float) else f"  Difference:          {values['difference']:+d}")

    # Show profile differences
    logger.info("\n--- Profile Comparison ---")
    logger.info(f"WITH specificity ({len(results_with_spec['selected_profiles'])} profiles):")
    for i, p in enumerate(results_with_spec['selected_profiles'][:10]):
        logger.info(f"  {i+1}. {p}")

    logger.info(f"\nWITHOUT specificity ({len(results_without_spec['selected_profiles'])} profiles):")
    for i, p in enumerate(results_without_spec['selected_profiles'][:10]):
        logger.info(f"  {i+1}. {p}")

    # Save results
    if args.output_dir:
        output_path = Path(args.output_dir)
    else:
        output_path = REPO_ROOT / "test_results" / "specificity_weighting"
    output_path.mkdir(parents=True, exist_ok=True)

    results = {
        "dataset": dataset_name,
        "with_specificity": results_with_spec,
        "without_specificity": results_without_spec,
        "comparison": comparison,
    }

    # Convert numpy arrays for JSON serialization
    def convert_for_json(obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, dict):
            return {k: convert_for_json(v) for k, v in obj.items()}
        elif isinstance(obj, list):
            return [convert_for_json(i) for i in obj]
        elif isinstance(obj, (np.int64, np.int32)):
            return int(obj)
        elif isinstance(obj, (np.float64, np.float32)):
            return float(obj)
        return obj

    results = convert_for_json(results)

    output_file = output_path / f"{dataset_name}_specificity_comparison.json"
    with open(output_file, "w") as f:
        json.dump(results, f, indent=2)
    logger.info(f"\nResults saved to: {output_file}")

    return results


if __name__ == "__main__":
    main()
