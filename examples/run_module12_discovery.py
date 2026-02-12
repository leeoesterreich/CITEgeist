#!/usr/bin/env python
"""
Module 1-2 Discovery Runner for Patient Samples.

Runs Modules 1-2 (marker interest detection and profile discovery) on
patient CITE-seq data and outputs discovered profiles to JSON for
comparison with curated profiles.

Usage:
    python run_module12_discovery.py --sample HCC22-088-P1-S1 --output-dir output/module12_discovery

    # With custom parameters:
    python run_module12_discovery.py --sample HCC22-088-P1-S1 \
        --top-k 3 --n-permutations 999 --output-dir output/module12_discovery
"""
import argparse
import json
import logging
import os
import sys
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import scipy.sparse
import squidpy as sq

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from CITEgeist.model import (
    identify_interesting_markers,
    analyze_marker_colocalization,
    select_profiles,
)
from CITEgeist.model.spatial_colocalization import discover_profiles_continuous

logger = logging.getLogger(__name__)
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)

# Patient data location
DATA_ROOT = Path("/ix1/alee/LO_LAB/General/Lab_Data/20250210_CITEGeistPublicData_GEO_Alex/processed_files")

# All 12 patient samples (deduplicated: P4-S2_1i_rep replaces P4-S2, P5-S2_F_rep replaces P5-S2)
ALL_SAMPLES = [
    "HCC22-088-P1-S1",
    "HCC22-088-P1-S2",
    "HCC22-088-P2-S1",
    "HCC22-088-P2-S2",
    "HCC22-088-P3-S1_A",
    "HCC22-088-P3-S2",
    "HCC22-088-P4-S1",
    "HCC22-088-P4-S2_1i_rep",
    "HCC22-088-P5-S1",
    "HCC22-088-P5-S2_F_rep",
    "HCC22-088-P6-S1",
    "HCC22-088-P6-S2_D",
]

# Original curated profiles for comparison
CURATED_PROFILES = {
    "Cancer Cells": ["EPCAM-1", "SDC1-1", "KRT5-1"],
    "Macrophages": ["CD68-1", "CD14-1"],
    "CD4 T Cells": ["CD3E-1", "CD4-1"],
    "CD8 T Cells": ["CD3E-1", "CD8A-1"],
    "B Cells": ["MS4A1-1", "CD19-1"],
    "Endothelial Cells": ["PECAM1-1"],
    "Fibroblasts": ["ACTA2-1"],
}


def load_patient_antibody_data(sample_name: str) -> tuple:
    """
    Load patient antibody data from SpaceRanger output.

    Returns:
        Tuple of (X_raw, coords, marker_names, n_spots)
    """
    sample_path = DATA_ROOT / sample_name / "outs"

    if not sample_path.exists():
        raise FileNotFoundError(f"Sample path not found: {sample_path}")

    logger.info(f"Loading data from {sample_path}")

    # Load combined data
    adata = sq.read.visium(
        str(sample_path),
        counts_file='filtered_feature_bc_matrix.h5',
        load_images=True,
        gex_only=False,
    )

    # Split to get antibody data only
    # Antibody features have feature_types == "Antibody Capture"
    is_antibody = adata.var['feature_types'] == 'Antibody Capture'
    antibody_adata = adata[:, is_antibody].copy()

    # Get raw counts (before any transformation)
    X_raw = antibody_adata.X
    if scipy.sparse.issparse(X_raw):
        X_raw = X_raw.toarray()

    # Get spatial coordinates
    coords = antibody_adata.obsm['spatial']

    # Filter spots with non-finite coordinates
    finite_mask = np.all(np.isfinite(coords), axis=1)
    if not np.all(finite_mask):
        n_invalid = (~finite_mask).sum()
        logger.warning(f"Filtering {n_invalid} spots with non-finite coordinates")
        coords = coords[finite_mask]
        X_raw = X_raw[finite_mask]

    marker_names = list(antibody_adata.var_names)
    n_spots = X_raw.shape[0]

    logger.info(f"Loaded {n_spots} spots, {len(marker_names)} markers")
    logger.info(f"Markers: {marker_names}")

    return X_raw, coords, marker_names, n_spots


def run_module12(
    X: np.ndarray,
    coords: np.ndarray,
    marker_names: List[str],
    top_k: int = 3,
    n_permutations: int = 999,
    fdr_threshold: float = 0.05,
    variance_target: float = 0.90,
    min_marginal_gain: float = 0.005,
    seed: int = 42,
) -> Dict:
    """
    Run Modules 1-2 on antibody expression data.

    Args:
        X: Raw antibody expression matrix (spots x markers)
        coords: Spatial coordinates (spots x 2)
        marker_names: List of marker names
        top_k: Mutual top-k for profile pair sparsification
        n_permutations: Number of permutations for significance testing
        fdr_threshold: FDR threshold for significant markers/pairs
        variance_target: Target fraction of spatial variance to explain
        min_marginal_gain: Minimum marginal variance gain per profile
        seed: Random seed for reproducibility

    Returns:
        Dict with module results and discovered profiles
    """
    results = {
        "parameters": {
            "top_k": top_k,
            "n_permutations": n_permutations,
            "fdr_threshold": fdr_threshold,
            "variance_target": variance_target,
            "min_marginal_gain": min_marginal_gain,
            "seed": seed,
        },
        "n_spots": X.shape[0],
        "n_markers": len(marker_names),
        "marker_names": marker_names,
    }

    # =========================================================================
    # Module 1: Marker Interest Detection
    # =========================================================================
    logger.info("Module 1: Identifying interesting markers...")

    m1_result = identify_interesting_markers(
        X=X,
        coords=coords,
        marker_names=marker_names,
        morans_k=8,
        morans_n_perm=min(n_permutations, 99),  # Use fewer perms for speed
        seed=seed,
        verbose=True,
    )

    interesting = m1_result.interesting_markers
    logger.info(f"Module 1: {len(interesting)}/{len(marker_names)} interesting markers")

    # Build score dictionaries from markers list
    marker_scores = {}
    for m in m1_result.markers:
        marker_scores[m.name] = {
            "kurtosis": float(m.kurtosis),
            "gmm_snr": float(m.gmm_snr),
            "morans_i": float(m.morans_i),
            "morans_pvalue": float(m.morans_i_pvalue),
            "passed_kurtosis": bool(m.passed_kurtosis),
            "passed_gmm": bool(m.passed_gmm),
            "passed_morans": bool(m.passed_morans),
        }

    results["module1"] = {
        "interesting_markers": interesting,
        "boring_markers": m1_result.boring_markers,
        "marker_scores": marker_scores,
    }

    if len(interesting) < 2:
        logger.warning("Fewer than 2 interesting markers - cannot discover profiles")
        results["module2a"] = {"error": "insufficient_markers"}
        results["module2b"] = {"error": "insufficient_markers"}
        results["module2c"] = {"error": "insufficient_markers"}
        results["discovered_profiles"] = []
        return results

    # =========================================================================
    # Module 2a: Marker Colocalization Analysis
    # =========================================================================
    logger.info("Module 2a: Analyzing marker colocalization...")

    coloc_result = analyze_marker_colocalization(
        X=X,
        coords=coords,
        marker_names=marker_names,
        markers_to_analyze=interesting,
        neighbor_k=8,
        n_permutations=n_permutations,
        seed=seed,
        verbose=True,
    )

    logger.info(f"Module 2a: {len(coloc_result.pairs)} marker pairs analyzed")

    # Store colocalization results
    results["module2a"] = {
        "n_pairs": len(coloc_result.pairs),
        "significant_pairs": [
            {
                "marker_a": p.marker_a,
                "marker_b": p.marker_b,
                "colocalization_score": float(p.colocalization_score),
                "bivariate_morans_pvalue": float(p.bivariate_morans_pvalue),
                "mutual_neighbor_enrichment": float(p.mutual_neighbor_enrichment),
            }
            for p in coloc_result.pairs
            if p.bivariate_morans_pvalue < fdr_threshold
        ],
    }

    # =========================================================================
    # Module 2b: Profile Discovery
    # =========================================================================
    logger.info(f"Module 2b: Discovering profiles (continuous, top_k={top_k})...")

    discovery_result = discover_profiles_continuous(
        colocalization_result=coloc_result,
        top_k=top_k,
        distance_metric="colocalization_score",
        seed=seed,
        verbose=True,
    )

    logger.info(f"Module 2b: {len(discovery_result.profiles)} candidate profiles")

    results["module2b"] = {
        "n_candidate_profiles": len(discovery_result.profiles),
        "candidate_profiles": [list(p) for p in discovery_result.profiles],
        "singletons": list(discovery_result.singletons) if hasattr(discovery_result, 'singletons') else [],
    }

    for i, profile in enumerate(discovery_result.profiles):
        logger.info(f"  Candidate {i+1}: {list(profile)}")

    if len(discovery_result.profiles) == 0:
        logger.warning("No profiles discovered")
        results["module2c"] = {"error": "no_profiles"}
        results["discovered_profiles"] = []
        return results

    # =========================================================================
    # Module 2c: Profile Selection by Spatial Variance
    # =========================================================================
    logger.info("Module 2c: Selecting profiles by spatial variance...")

    selection_result = select_profiles(
        X=X,
        coords=coords,
        marker_names=marker_names,
        profiles=discovery_result.profiles,
        interesting_markers=interesting,
        colocalization_result=coloc_result,
        min_spatial_explained=variance_target,
        min_protein_explained=variance_target,
        verbose=True,
    )

    selected_profiles = selection_result.selected_profiles
    n_selected = selection_result.optimal_n

    logger.info(f"Module 2c: Selected {n_selected} profiles")

    results["module2c"] = {
        "n_selected": int(n_selected),
        "variance_explained": [float(v) for v in selection_result.variance_explained],
        "marginal_gains": [float(g) for g in selection_result.marginal_gains] if hasattr(selection_result, 'marginal_gains') else [],
        "stopping_reason": str(selection_result.stopping_reason),
    }

    # Final discovered profiles
    discovered = []
    for i, profile in enumerate(selected_profiles):
        markers = list(profile) if not isinstance(profile, list) else profile
        discovered.append({
            "profile_id": i,
            "markers": markers,
            "n_markers": len(markers),
        })
        logger.info(f"  Selected {i+1}: {markers}")

    results["discovered_profiles"] = discovered

    return results


def main():
    parser = argparse.ArgumentParser(
        description="Run Module 1-2 discovery on patient samples"
    )
    parser.add_argument(
        "--sample",
        type=str,
        required=True,
        help=f"Sample name (one of {ALL_SAMPLES})",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="output/module12_discovery",
        help="Output directory for results",
    )
    parser.add_argument("--top-k", type=int, default=3, help="Mutual top-k for sparsification")
    parser.add_argument("--n-permutations", type=int, default=999, help="Permutations for significance")
    parser.add_argument("--fdr-threshold", type=float, default=0.05, help="FDR threshold")
    parser.add_argument("--variance-target", type=float, default=0.90, help="Target variance explained")
    parser.add_argument("--min-marginal-gain", type=float, default=0.005, help="Min marginal gain per profile")
    parser.add_argument("--seed", type=int, default=42, help="Random seed")

    args = parser.parse_args()

    # Validate sample name
    if args.sample not in ALL_SAMPLES:
        logger.error(f"Unknown sample: {args.sample}")
        logger.error(f"Valid samples: {ALL_SAMPLES}")
        sys.exit(1)

    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load data
    logger.info(f"Processing sample: {args.sample}")
    X, coords, marker_names, n_spots = load_patient_antibody_data(args.sample)

    # Run Module 1-2
    results = run_module12(
        X=X,
        coords=coords,
        marker_names=marker_names,
        top_k=args.top_k,
        n_permutations=args.n_permutations,
        fdr_threshold=args.fdr_threshold,
        variance_target=args.variance_target,
        min_marginal_gain=args.min_marginal_gain,
        seed=args.seed,
    )

    # Add metadata
    results["sample_name"] = args.sample
    results["timestamp"] = datetime.now().isoformat()
    results["curated_profiles"] = CURATED_PROFILES

    # Save results
    output_file = output_dir / f"{args.sample}_module12_discovery.json"
    with open(output_file, "w") as f:
        json.dump(results, f, indent=2)

    logger.info(f"Saved results to {output_file}")

    # Print summary
    print(f"\n{'='*60}")
    print(f"Module 1-2 Discovery Summary: {args.sample}")
    print(f"{'='*60}")
    print(f"Spots: {n_spots}")
    print(f"Markers: {len(marker_names)}")
    print(f"Interesting markers: {len(results['module1']['interesting_markers'])}")
    print(f"Discovered profiles: {len(results['discovered_profiles'])}")
    print()

    for p in results['discovered_profiles']:
        print(f"  Profile {p['profile_id']+1}: {p['markers']}")

    print(f"\nOutput: {output_file}")


if __name__ == "__main__":
    main()
