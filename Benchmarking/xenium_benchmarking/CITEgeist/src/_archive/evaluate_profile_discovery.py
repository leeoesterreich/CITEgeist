#!/usr/bin/env python
"""
Mode A: Profile Discovery Evaluation

Tests whether CITEgeist's Module 1+2 correctly discovers and groups expected
marker combinations from Xenium protein data, independent of deconvolution.

Expected marker groupings from RNA cluster analysis (k=10):
- CD8+ T cells: CD3E, CD8A (possibly GranzymeB, CD45)
- Macrophages: CD68, HLA-DR, CD163, CD16
- Mixed Immune: CD3E, HLA-DR, CD8A (lower intensity)
- Epithelial: PanCK, Vimentin (high)
- Myofibroblasts: alphaSMA, Vimentin
- Stromal: Vimentin (low), mixed
- Endothelial: CD31, Vimentin (high)
- B cells: CD20, CD45RA, CD45
- Proliferating T: CD3E (very high), PCNA, Ki-67
- Vascular Stromal: CD31 (low), Vimentin

Metrics:
1. Profile Recovery Rate: How many expected groupings were discovered?
2. Marker Grouping Accuracy: Are expected markers grouped together?
3. Profile Purity: Are discovered profiles contaminated with unexpected markers?
4. Coverage: What fraction of expected markers appear in discovered profiles?

Reference:
    Zhao et al. (2025). "Benchmarking cell type annotation methods for 10x
    Xenium spatial transcriptomics data." BMC Bioinformatics, 26(1), 25.
    https://doi.org/10.1186/s12859-025-06044-0
"""

import argparse
import json
import logging
import sys
from pathlib import Path
from typing import Dict, List, Set, Tuple

import numpy as np
import pandas as pd
import scanpy as sc

# Add CITEgeist to path
REPO_ROOT = Path(__file__).parent.parent.parent.parent.parent
sys.path.insert(0, str(REPO_ROOT))

from CITEgeist.model.marker_interest import identify_interesting_markers
from CITEgeist.model.spatial_colocalization import (
    analyze_marker_colocalization,
    discover_profiles,
)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)


# ============================================================================
# Expected Marker Groupings from RNA Cluster Analysis
# ============================================================================

# Expected marker groupings based on Xenium Kidney RCC protein profiles per cluster
# These are the "ground truth" groupings we expect Module 2 to discover
EXPECTED_GROUPINGS = {
    "T_cell_core": {"CD3E", "CD8A"},  # Must group together
    "Macrophage_core": {"CD68", "HLA-DR"},  # Must group together
    "B_cell_core": {"CD20", "CD45RA"},  # Must group together
    "Endothelial_core": {"CD31"},  # May group with Vimentin
    "Myofibroblast_core": {"alphaSMA", "Vimentin"},  # Often group together
    "Proliferation": {"PCNA", "Ki-67"},  # If both present
}

# Markers that should NOT be grouped together (mutually exclusive cell types)
FORBIDDEN_GROUPINGS = [
    {"CD3E", "CD68"},  # T cells vs Macrophages
    {"CD3E", "CD20"},  # T cells vs B cells
    {"CD68", "CD20"},  # Macrophages vs B cells
    {"PanCK", "CD68"},  # Epithelial vs Macrophages
    {"PanCK", "CD3E"},  # Epithelial vs T cells
]

# Markers that characterize specific cell types (for coverage calculation)
CELL_TYPE_MARKERS = {
    "CD8+ T cells": ["CD3E", "CD8A", "CD45", "GranzymeB"],
    "Macrophages": ["CD68", "HLA-DR", "CD163", "CD16"],
    "B cells": ["CD20", "CD45RA", "CD45"],
    "Endothelial": ["CD31", "Vimentin"],
    "Myofibroblasts": ["alphaSMA", "Vimentin"],
    "Epithelial": ["PanCK", "E-Cadherin", "Vimentin"],
    "Proliferating": ["PCNA", "Ki-67"],
}


# ============================================================================
# Evaluation Metrics
# ============================================================================


def compute_grouping_accuracy(
    discovered_profiles: List[List[str]],
    expected_groupings: Dict[str, Set[str]],
    available_markers: Set[str],
) -> Dict[str, float]:
    """
    Compute how well discovered profiles match expected marker groupings.

    Returns:
        Dict with:
        - recovery_rate: Fraction of expected groupings that were discovered
        - grouping_precision: Of expected pairs found together, how many are correct?
        - grouping_recall: Of expected pairs, how many were found together?
    """
    # Filter expected groupings to only include available markers
    filtered_expected = {}
    for name, markers in expected_groupings.items():
        filtered = markers & available_markers
        if len(filtered) >= 2:  # Need at least 2 markers to form a grouping
            filtered_expected[name] = filtered

    if not filtered_expected:
        return {"recovery_rate": 0.0, "grouping_precision": 0.0, "grouping_recall": 0.0}

    # Build discovered pairs (which markers are in the same profile)
    discovered_pairs = set()
    for profile in discovered_profiles:
        for i, m1 in enumerate(profile):
            for m2 in profile[i + 1 :]:
                discovered_pairs.add(frozenset([m1, m2]))

    # Build expected pairs
    expected_pairs = set()
    for markers in filtered_expected.values():
        markers_list = list(markers)
        for i, m1 in enumerate(markers_list):
            for m2 in markers_list[i + 1 :]:
                expected_pairs.add(frozenset([m1, m2]))

    # Count recovered groupings (at least 2 markers from expected found together)
    recovered = 0
    for name, markers in filtered_expected.items():
        markers_list = list(markers)
        # Check if at least 2 markers are in the same profile
        for profile in discovered_profiles:
            overlap = set(profile) & markers
            if len(overlap) >= 2:
                recovered += 1
                break

    recovery_rate = recovered / len(filtered_expected) if filtered_expected else 0.0

    # Precision and recall for pairs
    true_positives = len(discovered_pairs & expected_pairs)
    precision = true_positives / len(discovered_pairs) if discovered_pairs else 0.0
    recall = true_positives / len(expected_pairs) if expected_pairs else 0.0

    return {
        "recovery_rate": recovery_rate,
        "grouping_precision": precision,
        "grouping_recall": recall,
        "n_expected_groupings": len(filtered_expected),
        "n_recovered": recovered,
        "n_expected_pairs": len(expected_pairs),
        "n_discovered_pairs": len(discovered_pairs),
        "n_true_positive_pairs": true_positives,
    }


def compute_profile_purity(
    discovered_profiles: List[List[str]],
    forbidden_groupings: List[Set[str]],
    available_markers: Set[str],
) -> Dict[str, float]:
    """
    Compute profile purity (absence of forbidden marker combinations).

    Returns:
        Dict with:
        - purity: Fraction of profiles without forbidden groupings
        - n_violations: Number of forbidden combinations found
    """
    # Filter forbidden groupings to only include available markers
    filtered_forbidden = []
    for forbidden in forbidden_groupings:
        filtered = forbidden & available_markers
        if len(filtered) >= 2:
            filtered_forbidden.append(filtered)

    if not filtered_forbidden:
        return {"purity": 1.0, "n_violations": 0}

    violations = 0
    for profile in discovered_profiles:
        profile_set = set(profile)
        for forbidden in filtered_forbidden:
            # Check if all forbidden markers are in this profile
            if forbidden <= profile_set:
                violations += 1

    purity = 1.0 - (violations / len(discovered_profiles)) if discovered_profiles else 1.0

    return {"purity": max(0.0, purity), "n_violations": violations}


def compute_marker_coverage(
    discovered_profiles: List[List[str]],
    cell_type_markers: Dict[str, List[str]],
    available_markers: Set[str],
) -> Dict[str, float]:
    """
    Compute what fraction of expected cell type markers appear in profiles.

    Returns:
        Dict with coverage statistics per cell type and overall.
    """
    # Get all markers in discovered profiles
    discovered_markers = set()
    for profile in discovered_profiles:
        discovered_markers.update(profile)

    coverage_by_type = {}
    total_expected = 0
    total_covered = 0

    for cell_type, markers in cell_type_markers.items():
        available = [m for m in markers if m in available_markers]
        if not available:
            continue

        covered = [m for m in available if m in discovered_markers]
        coverage = len(covered) / len(available)
        coverage_by_type[cell_type] = {
            "coverage": coverage,
            "n_available": len(available),
            "n_covered": len(covered),
            "covered_markers": covered,
        }

        total_expected += len(available)
        total_covered += len(covered)

    overall_coverage = total_covered / total_expected if total_expected > 0 else 0.0

    return {
        "overall_coverage": overall_coverage,
        "per_cell_type": coverage_by_type,
    }


def compute_profile_statistics(
    discovered_profiles: List[List[str]],
    n_interesting_markers: int,
) -> Dict[str, float]:
    """
    Compute basic statistics about discovered profiles.
    """
    if not discovered_profiles:
        return {
            "n_profiles": 0,
            "n_multi_marker_profiles": 0,
            "n_singletons": 0,
            "mean_profile_size": 0.0,
            "max_profile_size": 0,
        }

    sizes = [len(p) for p in discovered_profiles]
    return {
        "n_profiles": len(discovered_profiles),
        "n_multi_marker_profiles": sum(1 for s in sizes if s > 1),
        "n_singletons": sum(1 for s in sizes if s == 1),
        "mean_profile_size": np.mean(sizes),
        "max_profile_size": max(sizes),
        "n_interesting_markers_input": n_interesting_markers,
    }


# ============================================================================
# Main Evaluation
# ============================================================================


def run_profile_discovery_evaluation(
    protein_adata: sc.AnnData,
    region_id: int = 0,
    output_dir: str = ".",
    fdr_alpha: float = 0.05,
    top_k: int = 3,
    kurtosis_threshold: float = 2.0,
    morans_threshold: float = 0.1,
    n_permutations: int = 199,
) -> Dict:
    """
    Run Module 1+2 profile discovery evaluation on Xenium protein data.

    Args:
        protein_adata: AnnData with protein expression (.X) and spatial coords (.obsm['spatial'])
        region_id: Region identifier for output files
        output_dir: Directory for output files
        fdr_alpha: FDR threshold for colocalization significance
        top_k: Top-k partners for mutual sparsification
        kurtosis_threshold: Threshold for Module 1 kurtosis gating
        morans_threshold: Moran's I threshold for spatial autocorrelation
        n_permutations: Number of permutations for significance testing

    Returns:
        Dict with all evaluation metrics
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Extract data
    X = protein_adata.X.toarray() if hasattr(protein_adata.X, "toarray") else protein_adata.X
    coords = protein_adata.obsm["spatial"]
    marker_names = list(protein_adata.var_names)

    logger.info(f"Region {region_id}: {X.shape[0]} spots, {len(marker_names)} markers")
    logger.info(f"Available markers: {marker_names}")

    # ========================================================================
    # Module 1: Identify Interesting Markers
    # ========================================================================
    logger.info("=" * 60)
    logger.info("MODULE 1: Marker Interest Detection")
    logger.info("=" * 60)

    interest_result = identify_interesting_markers(
        X=X,
        coords=coords,
        marker_names=marker_names,
        kurtosis_threshold=kurtosis_threshold,
        morans_threshold=morans_threshold,
        morans_k=8,
        morans_n_perm=n_permutations,
    )

    interesting_markers = interest_result.interesting_markers
    boring_markers = interest_result.boring_markers

    logger.info(f"Interesting markers: {interesting_markers}")
    logger.info(f"Boring markers: {boring_markers}")

    # Save Module 1 results
    interest_df = interest_result.to_dataframe()
    interest_df.to_csv(output_dir / f"region_{region_id}_module1_interest.csv", index=False)

    # ========================================================================
    # Module 2a: Colocalization Analysis
    # ========================================================================
    logger.info("=" * 60)
    logger.info("MODULE 2a: Pairwise Colocalization Analysis")
    logger.info("=" * 60)

    if len(interesting_markers) < 2:
        logger.warning("Not enough interesting markers for colocalization analysis")
        return {"error": "Not enough interesting markers", "n_interesting": len(interesting_markers)}

    coloc_result = analyze_marker_colocalization(
        X=X,
        coords=coords,
        marker_names=marker_names,
        markers_to_analyze=interesting_markers,
        neighbor_k=6,
        n_permutations=n_permutations,
    )

    # Save colocalization results
    coloc_df = coloc_result.to_dataframe()
    coloc_df.to_csv(output_dir / f"region_{region_id}_module2a_colocalization.csv", index=False)

    # ========================================================================
    # Module 2b: Profile Discovery
    # ========================================================================
    logger.info("=" * 60)
    logger.info("MODULE 2b: Profile Discovery")
    logger.info("=" * 60)

    profile_result = discover_profiles(
        colocalization_result=coloc_result,
        fdr_alpha=fdr_alpha,
        top_k=top_k,
        use_triangle_assembly=False,
        pvalue_source="bivariate_morans",
    )

    discovered_profiles = profile_result.profiles
    logger.info(profile_result.summary())

    # Save profile results
    profile_df = profile_result.to_dataframe()
    profile_df.to_csv(output_dir / f"region_{region_id}_module2b_profiles.csv", index=False)

    # ========================================================================
    # Evaluation Metrics
    # ========================================================================
    logger.info("=" * 60)
    logger.info("EVALUATION: Comparing to Expected Groupings")
    logger.info("=" * 60)

    available_markers = set(marker_names)

    # Compute all metrics
    grouping_metrics = compute_grouping_accuracy(
        discovered_profiles=discovered_profiles,
        expected_groupings=EXPECTED_GROUPINGS,
        available_markers=available_markers,
    )

    purity_metrics = compute_profile_purity(
        discovered_profiles=discovered_profiles,
        forbidden_groupings=FORBIDDEN_GROUPINGS,
        available_markers=available_markers,
    )

    coverage_metrics = compute_marker_coverage(
        discovered_profiles=discovered_profiles,
        cell_type_markers=CELL_TYPE_MARKERS,
        available_markers=available_markers,
    )

    profile_stats = compute_profile_statistics(
        discovered_profiles=discovered_profiles,
        n_interesting_markers=len(interesting_markers),
    )

    # Compile results
    results = {
        "region_id": region_id,
        "n_spots": X.shape[0],
        "n_markers_total": len(marker_names),
        "n_interesting_markers": len(interesting_markers),
        "interesting_markers": interesting_markers,
        "n_discovered_profiles": len(discovered_profiles),
        "discovered_profiles": [list(p) for p in discovered_profiles],
        "modularity": profile_result.modularity,
        "grouping_metrics": grouping_metrics,
        "purity_metrics": purity_metrics,
        "coverage_metrics": coverage_metrics,
        "profile_statistics": profile_stats,
        "parameters": {
            "fdr_alpha": fdr_alpha,
            "top_k": top_k,
            "kurtosis_threshold": kurtosis_threshold,
            "morans_threshold": morans_threshold,
            "n_permutations": n_permutations,
        },
    }

    # Log summary
    logger.info(f"Grouping Recovery Rate: {grouping_metrics['recovery_rate']:.1%}")
    logger.info(f"Grouping Precision: {grouping_metrics['grouping_precision']:.1%}")
    logger.info(f"Grouping Recall: {grouping_metrics['grouping_recall']:.1%}")
    logger.info(f"Profile Purity: {purity_metrics['purity']:.1%}")
    logger.info(f"Marker Coverage: {coverage_metrics['overall_coverage']:.1%}")
    logger.info(f"Modularity: {profile_result.modularity:.3f}")

    # Save full results
    with open(output_dir / f"region_{region_id}_evaluation.json", "w") as f:
        # Convert numpy types for JSON serialization
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

    return results


def main():
    parser = argparse.ArgumentParser(
        description="Evaluate CITEgeist Module 1+2 profile discovery against expected groupings"
    )
    parser.add_argument(
        "--region-id", type=int, default=0, help="Region ID to process (default: 0)"
    )
    parser.add_argument(
        "--input-dir",
        type=str,
        default=str(REPO_ROOT / "Benchmarking/xenium_pseudovisium/data_granular_gt"),
        help="Input directory with h5ad_objects/",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default=str(REPO_ROOT / "Benchmarking/xenium_benchmarking/CITEgeist/output_profile_discovery"),
        help="Output directory for evaluation results",
    )
    parser.add_argument(
        "--fdr-alpha", type=float, default=0.05, help="FDR threshold (default: 0.05)"
    )
    parser.add_argument(
        "--top-k", type=int, default=3, help="Top-k partners for sparsification (default: 3)"
    )
    parser.add_argument(
        "--n-permutations",
        type=int,
        default=199,
        help="Number of permutations for significance testing (default: 199)",
    )

    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)

    # Load protein data
    protein_path = input_dir / "h5ad_objects" / f"Xenium_region_{args.region_id}_CITE.h5ad"
    if not protein_path.exists():
        logger.error(f"Protein data not found: {protein_path}")
        sys.exit(1)

    logger.info(f"Loading protein data from: {protein_path}")
    protein_adata = sc.read_h5ad(protein_path)

    results = run_profile_discovery_evaluation(
        protein_adata=protein_adata,
        region_id=args.region_id,
        output_dir=str(output_dir),
        fdr_alpha=args.fdr_alpha,
        top_k=args.top_k,
        n_permutations=args.n_permutations,
    )

    logger.info("=" * 60)
    logger.info("EVALUATION COMPLETE")
    logger.info("=" * 60)
    logger.info(f"Results saved to: {output_dir}")


if __name__ == "__main__":
    main()
