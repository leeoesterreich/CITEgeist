"""
Test adaptive colocalization enhancements on simulated datasets.

Compares profile discovery with different enhancement configurations:
- baseline: No enhancements
- expr_fallback: Expression correlation fallback for low spatial signal
- multiscale: Multi-scale neighborhood analysis
- both: Both enhancements combined

Usage:
    python tests/test_adaptive_benchmark.py --dataset high_seg --replicate 0 --config baseline
    python tests/test_adaptive_benchmark.py --dataset mixed --config expr_fallback
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
    get_marker_morans_dict,
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


def count_cell_types_in_profiles(profiles: List[List[str]], gt_cell_types: List[str] = None) -> Dict[str, Any]:
    """
    Count how many distinct cell types are represented in discovered profiles.

    For simulated data, markers are named like "B-cells_Protein_1", so we can
    extract the cell type from the marker name.

    Args:
        profiles: List of profiles, each containing marker names
        gt_cell_types: Ground truth cell types (optional)

    Returns:
        Dictionary with cell type counting metrics
    """
    # Extract cell types from marker names (format: "CellType_Protein_N")
    discovered_cell_types = set()
    for profile in profiles:
        for marker in profile:
            if "_Protein_" in marker:
                cell_type = marker.rsplit("_Protein_", 1)[0]
                discovered_cell_types.add(cell_type)

    # Ground truth cell types (9 total in Wu simulation)
    gt_types = gt_cell_types or [
        "B-cells", "CAFs", "Cancer_Epithelial", "Endothelial",
        "Myeloid", "Normal_Epithelial", "PVL", "Plasmablasts", "T-cells"
    ]

    found = discovered_cell_types.intersection(set(gt_types))

    return {
        "n_found": len(found),
        "n_total": len(gt_types),
        "found_types": sorted(list(found)),
        "missing_types": sorted(list(set(gt_types) - found)),
        "extra_types": sorted(list(discovered_cell_types - set(gt_types))),
    }


def run_profile_discovery(
    X: np.ndarray,
    coords: np.ndarray,
    marker_names: List[str],
    use_expr_fallback: bool = False,
    use_multiscale: bool = False,
    fdr_alpha: float = 0.05,
    top_k: int = 3,
) -> Dict[str, Any]:
    """
    Run the full profile discovery pipeline with optional enhancements.

    Args:
        X: Expression matrix (n_spots, n_markers)
        coords: Spatial coordinates (n_spots, 2)
        marker_names: List of marker names
        use_expr_fallback: Enable expression correlation fallback
        use_multiscale: Enable multi-scale neighborhood analysis
        fdr_alpha: FDR threshold
        top_k: Number of top partners in mutual top-k

    Returns:
        Dictionary with discovery results and metrics
    """
    config_str = f"expr_fallback={use_expr_fallback}, multiscale={use_multiscale}"
    logger.info(f"Running profile discovery ({config_str})...")

    # Module 1: Identify interesting markers
    interest_result = identify_interesting_markers(
        X, coords, marker_names,
        morans_threshold=0.1,
        morans_n_perm=199,
        verbose=True,
    )
    interesting_markers = interest_result.interesting_markers
    logger.info(f"Module 1: {len(interesting_markers)} interesting markers identified")

    # Get marker Moran's I for adaptive blending
    marker_morans_i = None
    if use_expr_fallback:
        marker_morans_i = get_marker_morans_dict(interest_result)
        avg_morans = np.mean(list(marker_morans_i.values()))
        logger.info(f"Average marker Moran's I: {avg_morans:.3f}")

    # Multi-scale parameters
    multi_scale_k = [6, 12, 24] if use_multiscale else None

    # Module 2a: Pairwise colocalization with enhancements
    coloc_result = analyze_marker_colocalization(
        X, coords, marker_names,
        markers_to_analyze=interesting_markers,
        marker_morans_i=marker_morans_i,  # For expression fallback
        multi_scale_k=multi_scale_k,      # For multi-scale
        n_permutations=199,
        verbose=True,
    )

    # Log enhancement info
    logger.info(f"Blend mode: {coloc_result.blend_mode}")
    logger.info(f"Blend weights: {coloc_result.blend_weights}")

    # Check for multi-scale info
    if use_multiscale and coloc_result.marker_pairs:
        scales_found = 0
        for pair in coloc_result.marker_pairs:
            if pair.bivariate_morans_per_scale:
                scales_found += 1
        logger.info(f"Multi-scale computed for {scales_found}/{len(coloc_result.marker_pairs)} pairs")

    # Module 2b: Profile discovery
    discovery_result = discover_profiles(
        coloc_result,
        fdr_alpha=fdr_alpha,
        top_k=top_k,
        verbose=True,
    )

    logger.info(f"Module 2b: {len(discovery_result.profiles)} profiles discovered")
    logger.info(f"  Modularity: {discovery_result.modularity:.3f}")

    # Module 2c: Profile selection
    selection_result = select_profiles(
        X, coords, marker_names,
        discovery_result.profiles,
        interesting_markers,
        coloc_result,
        min_marginal_gain=0.005,
        max_profiles=15,
        verbose=True,
    )

    logger.info(f"Module 2c: {selection_result.optimal_n} profiles selected")
    logger.info(f"  Variance explained: {selection_result.variance_explained[-1] if len(selection_result.variance_explained) > 0 else 0:.1%}")

    # Count cell types found
    cell_type_metrics = count_cell_types_in_profiles(selection_result.selected_profiles)
    logger.info(f"Cell types found: {cell_type_metrics['n_found']}/{cell_type_metrics['n_total']}")
    logger.info(f"  Found: {cell_type_metrics['found_types']}")
    logger.info(f"  Missing: {cell_type_metrics['missing_types']}")

    return {
        "config": {
            "use_expr_fallback": use_expr_fallback,
            "use_multiscale": use_multiscale,
        },
        "blend_mode": coloc_result.blend_mode,
        "blend_weights": coloc_result.blend_weights,
        "avg_marker_morans_i": coloc_result.avg_marker_morans_i,
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
        "cell_type_metrics": cell_type_metrics,
    }


def main():
    parser = argparse.ArgumentParser(description="Test adaptive colocalization enhancements")
    parser.add_argument(
        "--dataset",
        type=str,
        choices=["high_seg", "mixed"],
        default="high_seg",
        help="Dataset variant to test on",
    )
    parser.add_argument(
        "--replicate",
        type=int,
        default=0,
        help="Replicate index (0-4)",
    )
    parser.add_argument(
        "--config",
        type=str,
        choices=["baseline", "expr_fallback", "multiscale", "both"],
        default="baseline",
        help="Enhancement configuration",
    )
    parser.add_argument(
        "--use-expr-fallback",
        type=str,
        default="false",
        help="Enable expression fallback (true/false)",
    )
    parser.add_argument(
        "--use-multiscale",
        type=str,
        default="false",
        help="Enable multi-scale analysis (true/false)",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default=None,
        help="Output directory for results",
    )
    args = parser.parse_args()

    # Parse boolean flags
    use_expr_fallback = args.use_expr_fallback.lower() == "true"
    use_multiscale = args.use_multiscale.lower() == "true"

    # Load data
    X, coords, marker_names, gt_props = load_simulated_data(args.replicate, args.dataset)
    dataset_name = f"{args.dataset}_rep{args.replicate}"

    logger.info(f"\n{'='*60}")
    logger.info(f"Testing adaptive colocalization on {dataset_name}")
    logger.info(f"Config: {args.config}")
    logger.info(f"Expression fallback: {use_expr_fallback}")
    logger.info(f"Multi-scale: {use_multiscale}")
    logger.info(f"{'='*60}")

    # Run profile discovery
    results = run_profile_discovery(
        X, coords, marker_names,
        use_expr_fallback=use_expr_fallback,
        use_multiscale=use_multiscale,
    )

    # Add metadata
    results["dataset"] = args.dataset
    results["replicate"] = args.replicate
    results["config_name"] = args.config

    # Summary
    logger.info(f"\n{'='*60}")
    logger.info("SUMMARY")
    logger.info(f"{'='*60}")
    logger.info(f"Config: {args.config}")
    logger.info(f"Blend mode: {results['blend_mode']}")
    logger.info(f"Profiles discovered: {results['n_discovered_profiles']}")
    logger.info(f"Profiles selected: {results['n_selected_profiles']}")
    logger.info(f"Cell types found: {results['cell_type_metrics']['n_found']}/9")
    logger.info(f"Variance explained: {results['variance_explained']:.1%}")

    # Save results
    if args.output_dir:
        output_dir = Path(args.output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        output_file = output_dir / f"{dataset_name}_{args.config}.json"

        # Convert numpy arrays to lists for JSON serialization
        def convert_to_serializable(obj):
            if isinstance(obj, np.ndarray):
                return obj.tolist()
            if isinstance(obj, np.integer):
                return int(obj)
            if isinstance(obj, np.floating):
                return float(obj)
            if isinstance(obj, dict):
                return {k: convert_to_serializable(v) for k, v in obj.items()}
            if isinstance(obj, list):
                return [convert_to_serializable(v) for v in obj]
            return obj

        serializable_results = convert_to_serializable(results)

        with open(output_file, "w") as f:
            json.dump(serializable_results, f, indent=2)

        logger.info(f"Results saved to: {output_file}")

    return results


if __name__ == "__main__":
    main()
