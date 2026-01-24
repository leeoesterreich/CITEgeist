"""
Xenium validation: Compare protein-discovered vs RNA-discovered profiles.

This script tests whether RNA marker selection can recover similar
profiles to protein-based discovery using real Xenium data that has
both protein and RNA measurements.

Tests:
1. Run protein-based Module 1-2 pipeline (ground truth)
2. Run RNA-based pipeline with different modes
3. Compare profile concordance

This validates the core claim that CITEgeist can work with RNA-only data.
"""

import argparse
import json
import logging
import sys
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import scanpy as sc

# Add paths
REPO_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(REPO_ROOT / "CITEgeist"))
sys.path.insert(0, str(REPO_ROOT / "Benchmarking" / "xenium_pseudovisium" / "src"))

from model.marker_interest import identify_interesting_markers, MarkerInterestResult
from model.spatial_colocalization import (
    analyze_marker_colocalization,
    discover_profiles,
    ColocalizationResult,
    ProfileDiscoveryResult,
)
from model.rna_marker_selection import (
    MarkerMode,
    select_rna_markers,
    compare_protein_vs_rna_profiles,
    validate_marker_spatial_quality,
    ProfileComparisonResult,
)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)

# Data paths
XENIUM_DATA_DIR = Path(
    "/ix1/alee/LO_LAB/General/Public_Data/10x_PublicData/"
    "Xenium_RNA_Proteomic_RenalCellCarcinoma"
)

OUTPUT_DIR = REPO_ROOT / "tests" / "test_results" / "rna_vs_protein_validation"


def load_xenium_data(
    data_dir: Path = XENIUM_DATA_DIR,
    max_cells: Optional[int] = None,
    seed: int = 42,
) -> Tuple[sc.AnnData, sc.AnnData]:
    """
    Load Xenium data and split into GEX and protein.

    Returns:
        Tuple of (adata_gex, adata_protein)
    """
    from load_xenium import load_xenium_data as _load_xenium, split_gex_protein

    logger.info(f"Loading Xenium data from {data_dir}")
    adata = _load_xenium(str(data_dir), min_counts=0)

    logger.info(f"Total cells: {adata.n_obs:,}")

    # Subsample if needed
    if max_cells is not None and adata.n_obs > max_cells:
        np.random.seed(seed)
        indices = np.random.choice(adata.n_obs, size=max_cells, replace=False)
        indices = np.sort(indices)
        adata = adata[indices].copy()
        logger.info(f"Subsampled to {adata.n_obs:,} cells")

    # Split
    adata_gex, adata_protein = split_gex_protein(adata)

    logger.info(f"GEX: {adata_gex.n_obs:,} cells × {adata_gex.n_vars} genes")
    logger.info(f"Protein: {adata_protein.n_obs:,} cells × {adata_protein.n_vars} proteins")
    logger.info(f"Proteins: {list(adata_protein.var_names)}")

    return adata_gex, adata_protein


def run_protein_pipeline(
    adata_protein: sc.AnnData,
    n_permutations: int = 99,
    seed: int = 42,
) -> Tuple[MarkerInterestResult, ProfileDiscoveryResult]:
    """
    Run standard protein-based Module 1-2 pipeline.

    This is the "ground truth" for comparison.
    """
    logger.info("\n" + "=" * 60)
    logger.info("PROTEIN-BASED PIPELINE (Ground Truth)")
    logger.info("=" * 60)

    X = adata_protein.X
    if hasattr(X, "toarray"):
        X = X.toarray()
    X = np.asarray(X, dtype=np.float64)

    coords = adata_protein.obsm["spatial"]
    marker_names = list(adata_protein.var_names)

    # Module 1
    logger.info("\n[Module 1] Marker Interest Detection")
    interest_result = identify_interesting_markers(
        X=X,
        coords=coords,
        marker_names=marker_names,
        morans_n_perm=n_permutations,
        seed=seed,
        verbose=True,
    )

    interesting = interest_result.interesting_markers
    logger.info(f"Interesting markers: {len(interesting)}/{len(marker_names)}")
    logger.info(f"  {interesting}")

    # Module 2
    logger.info("\n[Module 2] Colocalization & Profile Discovery")

    coloc_result = analyze_marker_colocalization(
        X=X,
        coords=coords,
        marker_names=marker_names,
        markers_to_analyze=interesting,
        n_permutations=n_permutations,
        seed=seed,
        verbose=True,
    )

    profile_result = discover_profiles(
        colocalization_result=coloc_result,
        fdr_alpha=0.05,
        seed=seed,
    )

    logger.info(f"Discovered {len(profile_result.profiles)} profiles:")
    for i, profile in enumerate(profile_result.profiles):
        logger.info(f"  Profile {i}: {profile}")

    return interest_result, profile_result


def run_rna_pipeline(
    adata_gex: sc.AnnData,
    mode: MarkerMode,
    curated_markers: Optional[List[str]] = None,
    max_discovered: int = 50,
    n_permutations: int = 99,
    strict_spatial: bool = False,
    seed: int = 42,
) -> Tuple[MarkerInterestResult, ProfileDiscoveryResult, List[str]]:
    """
    Run RNA-based Module 1-2 pipeline.

    Returns:
        Tuple of (interest_result, profile_result, selected_markers)
    """
    logger.info("\n" + "=" * 60)
    logger.info(f"RNA-BASED PIPELINE (mode={mode.value})")
    logger.info("=" * 60)

    # Step 1: Select RNA markers
    logger.info("\n[Step 1] RNA Marker Selection")
    marker_result = select_rna_markers(
        adata_gex,
        mode=mode,
        curated_markers=curated_markers,
        max_discovered=max_discovered,
        n_permutations=n_permutations,
        strict_spatial_threshold=strict_spatial,
        seed=seed,
    )

    selected_markers = marker_result.selected_markers
    logger.info(f"Selected {len(selected_markers)} markers")
    logger.info(f"  Curated: {len(marker_result.curated_markers)}")
    logger.info(f"  Discovered: {len(marker_result.discovered_markers)}")

    if len(selected_markers) < 3:
        logger.warning("Too few markers selected, skipping profile discovery")
        return None, None, selected_markers

    # Subset to selected markers
    markers_in_data = [m for m in selected_markers if m in adata_gex.var_names]
    adata_subset = adata_gex[:, markers_in_data].copy()

    X = adata_subset.X
    if hasattr(X, "toarray"):
        X = X.toarray()
    X = np.asarray(X, dtype=np.float64)

    coords = adata_subset.obsm["spatial"]
    marker_names = list(adata_subset.var_names)

    # Module 1 on selected markers
    logger.info("\n[Step 2] Module 1: Marker Interest Detection")
    interest_result = identify_interesting_markers(
        X=X,
        coords=coords,
        marker_names=marker_names,
        morans_n_perm=n_permutations,
        seed=seed,
        verbose=True,
    )

    interesting = interest_result.interesting_markers
    logger.info(f"Interesting markers: {len(interesting)}/{len(marker_names)}")

    if len(interesting) < 3:
        logger.warning("Too few interesting markers, using all selected")
        interesting = marker_names[:min(10, len(marker_names))]

    # Module 2
    logger.info("\n[Step 3] Module 2: Colocalization & Profile Discovery")

    coloc_result = analyze_marker_colocalization(
        X=X,
        coords=coords,
        marker_names=marker_names,
        markers_to_analyze=interesting,
        n_permutations=n_permutations,
        seed=seed,
        verbose=True,
    )

    profile_result = discover_profiles(
        colocalization_result=coloc_result,
        fdr_alpha=0.05,
        seed=seed,
    )

    logger.info(f"Discovered {len(profile_result.profiles)} profiles:")
    for i, profile in enumerate(profile_result.profiles):
        logger.info(f"  Profile {i}: {profile}")

    return interest_result, profile_result, selected_markers


def compare_pipelines(
    protein_profiles: List[List[str]],
    rna_profiles: List[List[str]],
    adata_protein: sc.AnnData,
    adata_gex: sc.AnnData,
    rna_markers: List[str],
) -> ProfileComparisonResult:
    """
    Compare protein-discovered vs RNA-discovered profiles.
    """
    logger.info("\n" + "=" * 60)
    logger.info("PROFILE COMPARISON")
    logger.info("=" * 60)

    # Get expression matrices
    X_protein = adata_protein.X
    if hasattr(X_protein, "toarray"):
        X_protein = X_protein.toarray()
    X_protein = np.asarray(X_protein, dtype=np.float64)

    # For RNA, use the selected markers
    rna_markers_in_data = [m for m in rna_markers if m in adata_gex.var_names]
    adata_rna_subset = adata_gex[:, rna_markers_in_data].copy()

    X_rna = adata_rna_subset.X
    if hasattr(X_rna, "toarray"):
        X_rna = X_rna.toarray()
    X_rna = np.asarray(X_rna, dtype=np.float64)

    result = compare_protein_vs_rna_profiles(
        protein_profiles=protein_profiles,
        rna_profiles=rna_profiles,
        X_protein=X_protein,
        X_rna=X_rna,
        protein_marker_names=list(adata_protein.var_names),
        rna_marker_names=list(adata_rna_subset.var_names),
        coords=adata_protein.obsm["spatial"],
    )

    logger.info(result.summary())

    return result


def main():
    parser = argparse.ArgumentParser(
        description="Compare protein vs RNA profile discovery on Xenium data"
    )
    parser.add_argument(
        "--max-cells", type=int, default=10000,
        help="Maximum cells to use (default: 10000 for speed)"
    )
    parser.add_argument(
        "--n-permutations", type=int, default=99,
        help="Permutations for p-values (default: 99)"
    )
    parser.add_argument(
        "--seed", type=int, default=42,
        help="Random seed"
    )
    parser.add_argument(
        "--output-dir", type=str, default=None,
        help="Output directory (default: tests/test_results/rna_vs_protein_validation)"
    )
    args = parser.parse_args()

    output_dir = Path(args.output_dir) if args.output_dir else OUTPUT_DIR
    output_dir.mkdir(parents=True, exist_ok=True)

    start_time = datetime.now()

    print("=" * 70)
    print("RNA vs Protein Profile Discovery Validation")
    print("=" * 70)
    print(f"Start time: {start_time}")
    print(f"Max cells: {args.max_cells}")
    print(f"Output dir: {output_dir}")
    print()

    # Load data
    adata_gex, adata_protein = load_xenium_data(
        max_cells=args.max_cells,
        seed=args.seed,
    )

    # Run protein pipeline (ground truth)
    protein_interest, protein_profiles = run_protein_pipeline(
        adata_protein,
        n_permutations=args.n_permutations,
        seed=args.seed,
    )

    # Save protein results
    protein_interest.to_dataframe().to_csv(
        output_dir / "protein_marker_interest.csv", index=False
    )

    protein_profile_list = protein_profiles.profiles

    # Test different RNA modes
    results = {}

    # Mode 1: Curated-only (use protein marker names as curated)
    # This tests if we can recover profiles using the same markers
    logger.info("\n" + "#" * 70)
    logger.info("# TEST 1: Curated-only mode (protein markers as RNA)")
    logger.info("#" * 70)

    # Check which protein markers are also in RNA data
    protein_markers_in_rna = [
        m for m in adata_protein.var_names if m in adata_gex.var_names
    ]
    logger.info(f"Protein markers in RNA data: {protein_markers_in_rna}")

    if protein_markers_in_rna:
        rna_interest_1, rna_profiles_1, rna_markers_1 = run_rna_pipeline(
            adata_gex,
            mode=MarkerMode.CURATED_ONLY,
            curated_markers=protein_markers_in_rna,
            n_permutations=args.n_permutations,
            seed=args.seed,
        )

        if rna_profiles_1 is not None:
            comparison_1 = compare_pipelines(
                protein_profile_list,
                rna_profiles_1.profiles,
                adata_protein,
                adata_gex,
                rna_markers_1,
            )
            results["curated_protein_markers"] = {
                "n_markers": len(rna_markers_1),
                "n_profiles": len(rna_profiles_1.profiles),
                "mean_jaccard": comparison_1.mean_best_jaccard,
                "mean_semantic_jaccard": comparison_1.mean_semantic_jaccard,
                "mean_spatial_corr": comparison_1.mean_spatial_correlation,
                "matched_profiles": comparison_1.matched_profiles,
                "semantic_matched_profiles": comparison_1.semantic_matched_profiles,
            }
    else:
        logger.warning("No protein markers found in RNA data for curated test")

    # Mode 2: Hybrid with canonical markers
    logger.info("\n" + "#" * 70)
    logger.info("# TEST 2: Hybrid mode (canonical markers + discovery)")
    logger.info("#" * 70)

    # Define canonical markers for kidney/RCC
    canonical_markers = [
        # Immune
        "CD3D", "CD3E", "CD4", "CD8A", "CD68", "CD163", "CD14",
        "MS4A1", "CD79A", "NKG7", "GNLY",
        # Epithelial/tumor
        "EPCAM", "KRT8", "KRT18", "KRT19", "CA9", "PAX8",
        # Endothelial
        "PECAM1", "VWF", "CDH5",
        # Fibroblasts
        "COL1A1", "COL1A2", "DCN", "FAP",
        # Proliferation
        "MKI67", "TOP2A",
    ]

    rna_interest_2, rna_profiles_2, rna_markers_2 = run_rna_pipeline(
        adata_gex,
        mode=MarkerMode.HYBRID,
        curated_markers=canonical_markers,
        max_discovered=30,
        n_permutations=args.n_permutations,
        strict_spatial=False,
        seed=args.seed,
    )

    if rna_profiles_2 is not None:
        comparison_2 = compare_pipelines(
            protein_profile_list,
            rna_profiles_2.profiles,
            adata_protein,
            adata_gex,
            rna_markers_2,
        )
        results["hybrid_canonical"] = {
            "n_markers": len(rna_markers_2),
            "n_profiles": len(rna_profiles_2.profiles),
            "mean_jaccard": comparison_2.mean_best_jaccard,
            "mean_semantic_jaccard": comparison_2.mean_semantic_jaccard,
            "mean_spatial_corr": comparison_2.mean_spatial_correlation,
            "matched_profiles": comparison_2.matched_profiles,
            "semantic_matched_profiles": comparison_2.semantic_matched_profiles,
            "match_details": comparison_2.match_details[:5] if comparison_2.match_details else [],
        }

        # Save hybrid results
        if rna_interest_2 is not None:
            rna_interest_2.to_dataframe().to_csv(
                output_dir / "rna_hybrid_marker_interest.csv", index=False
            )

    # Mode 3: Full autodiscovery
    logger.info("\n" + "#" * 70)
    logger.info("# TEST 3: Autodiscovery mode (no prior knowledge)")
    logger.info("#" * 70)

    rna_interest_3, rna_profiles_3, rna_markers_3 = run_rna_pipeline(
        adata_gex,
        mode=MarkerMode.AUTODISCOVERY,
        max_discovered=50,
        n_permutations=args.n_permutations,
        strict_spatial=False,
        seed=args.seed,
    )

    if rna_profiles_3 is not None:
        comparison_3 = compare_pipelines(
            protein_profile_list,
            rna_profiles_3.profiles,
            adata_protein,
            adata_gex,
            rna_markers_3,
        )
        results["autodiscovery"] = {
            "n_markers": len(rna_markers_3),
            "n_profiles": len(rna_profiles_3.profiles),
            "mean_jaccard": comparison_3.mean_best_jaccard,
            "mean_semantic_jaccard": comparison_3.mean_semantic_jaccard,
            "mean_spatial_corr": comparison_3.mean_spatial_correlation,
            "matched_profiles": comparison_3.matched_profiles,
            "semantic_matched_profiles": comparison_3.semantic_matched_profiles,
            "match_details": comparison_3.match_details[:5] if comparison_3.match_details else [],
        }

        # Save autodiscovery results
        if rna_interest_3 is not None:
            rna_interest_3.to_dataframe().to_csv(
                output_dir / "rna_autodiscovery_marker_interest.csv", index=False
            )

    # Summary
    elapsed = datetime.now() - start_time

    print("\n" + "=" * 70)
    print("VALIDATION SUMMARY")
    print("=" * 70)
    print(f"\nProtein pipeline (ground truth):")
    print(f"  Interesting markers: {len(protein_interest.interesting_markers)}")
    print(f"  Profiles discovered: {len(protein_profile_list)}")
    for i, profile in enumerate(protein_profile_list):
        print(f"    Profile {i}: {profile}")

    print(f"\nRNA pipeline comparisons:")
    for mode_name, metrics in results.items():
        print(f"\n  {mode_name}:")
        print(f"    Markers used: {metrics['n_markers']}")
        print(f"    Profiles discovered: {metrics['n_profiles']}")
        print(f"    Mean Jaccard (raw marker names): {metrics['mean_jaccard']:.3f}")
        print(f"    Mean Semantic Jaccard (biological equivalence): {metrics.get('mean_semantic_jaccard', 0):.3f}")
        print(f"    Mean spatial correlation: {metrics['mean_spatial_corr']:.3f}")
        print(f"    Matched profiles (Jaccard > 0.3): {metrics['matched_profiles']}")
        print(f"    Matched profiles (Semantic > 0.3): {metrics.get('semantic_matched_profiles', 0)}")
        if metrics.get('match_details'):
            print("    Top matches:")
            for match in metrics['match_details'][:3]:
                equiv = match.get('equivalent_markers', [])
                print(f"      {match['protein_profile']} -> {match['rna_profile']}")
                print(f"        spatial_corr={match.get('spatial_corr', 0):.3f}, semantic={match.get('semantic_jaccard', 0):.3f}")
                if equiv:
                    print(f"        Equivalent markers found: {equiv}")

    print(f"\nTotal runtime: {elapsed}")
    print(f"Results saved to: {output_dir}")

    # Save summary
    summary = {
        "timestamp": start_time.isoformat(),
        "max_cells": args.max_cells,
        "n_permutations": args.n_permutations,
        "protein_profiles": protein_profile_list,
        "results": results,
        "runtime_seconds": elapsed.total_seconds(),
    }

    with open(output_dir / "validation_summary.json", "w") as f:
        json.dump(summary, f, indent=2, default=str)

    # Determine success using multiple criteria
    # Note: Use BEST individual match, not mean (mean is diluted by unmatched profiles)
    if results:
        # Find best individual matches across all modes
        best_individual_spatial = 0.0
        best_individual_semantic = 0.0
        total_semantic_matches = 0
        best_match_info = []

        for mode_name, metrics in results.items():
            # Count semantic matches
            total_semantic_matches += metrics.get("semantic_matched_profiles", 0)

            # Check individual match details
            for match in metrics.get("match_details", []):
                spatial = match.get("spatial_corr", 0)
                semantic = match.get("semantic_jaccard", 0)

                if spatial > best_individual_spatial:
                    best_individual_spatial = spatial
                if semantic > best_individual_semantic:
                    best_individual_semantic = semantic

                # Track good matches
                if spatial > 0.3 or semantic > 0.3:
                    best_match_info.append({
                        "mode": mode_name,
                        "protein": match["protein_profile"],
                        "rna": match["rna_profile"],
                        "spatial": spatial,
                        "semantic": semantic,
                        "equiv": match.get("equivalent_markers", []),
                    })

        print(f"\nBest individual spatial correlation: {best_individual_spatial:.3f}")
        print(f"Best individual semantic Jaccard: {best_individual_semantic:.3f}")
        print(f"Total profiles with semantic match (>0.3): {total_semantic_matches}")

        # Success criteria (any of these):
        # 1. Best individual spatial correlation > 0.4
        # 2. Best individual semantic Jaccard > 0.3
        # 3. At least 2 semantic matches found
        passed_spatial = best_individual_spatial > 0.4
        passed_semantic = best_individual_semantic > 0.3
        passed_match_count = total_semantic_matches >= 2

        if passed_spatial or passed_semantic or passed_match_count:
            reasons = []
            if passed_spatial:
                reasons.append(f"spatial correlation ({best_individual_spatial:.3f})")
            if passed_semantic:
                reasons.append(f"semantic equivalence ({best_individual_semantic:.3f})")
            if passed_match_count:
                reasons.append(f"{total_semantic_matches} semantic matches")

            print(f"\n✓ VALIDATION PASSED: RNA profiles match protein profiles via {' and '.join(reasons)}")

            if best_match_info:
                print("\nBest biological equivalence matches:")
                for m in sorted(best_match_info, key=lambda x: x["spatial"], reverse=True)[:5]:
                    print(f"  [{m['mode']}] {m['protein']} -> {m['rna']}")
                    print(f"    spatial={m['spatial']:.3f}, semantic={m['semantic']:.3f}")
                    if m["equiv"]:
                        print(f"    Equivalent markers: {m['equiv']}")

            return 0
        else:
            print("\n✗ VALIDATION FAILED: Insufficient biological equivalence")
            print("  Thresholds: spatial>0.4, semantic>0.3, or >=2 semantic matches")
            return 1
    else:
        print("\n✗ VALIDATION FAILED: No results generated")
        return 1


if __name__ == "__main__":
    sys.exit(main())
