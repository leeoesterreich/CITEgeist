#!/usr/bin/env python3
"""
Quick diagnostic to test Module 2b with relaxed parameters.

Tests different top_k and min_score thresholds to see if we can get
more consistent profiles across regions.
"""

import json
import sys
from pathlib import Path

# Add CITEgeist to path
sys.path.insert(0, str(Path(__file__).parent.parent))

import numpy as np
import scanpy as sc
from collections import defaultdict

from CITEgeist.model.discovery.spatial_colocalization import (
    analyze_marker_colocalization,
    discover_profiles,
)

import pytest

# Paths (relative to repo root, resolved from this test file's location)
_REPO_ROOT = Path(__file__).resolve().parent.parent
XENIUM_BASE = _REPO_ROOT / "Benchmarking" / "xenium_pseudovisium"
OUTPUT_DIR = Path(__file__).resolve().parent / "test_results"
if not XENIUM_BASE.exists():
    pytest.skip("Xenium pseudovisium dataset not available", allow_module_level=True)
OUTPUT_DIR.mkdir(exist_ok=True)

# Ideal 7-type profiles for comparison
IDEAL_7 = {
    "B cells": {"CD20", "CD45RA"},
    "CD4+ T cells": {"CD3E", "CD4", "CD45RO"},
    "CD8+ T cells": {"CD3E", "CD8A"},
    "Macrophages": {"CD68", "CD163", "CD16"},
    "Endothelial": {"CD31"},
    "Epithelial": {"PanCK", "E-Cadherin"},
    "Fibroblasts": {"alphaSMA", "Vimentin"},
}


def load_region_data(region_id: int):
    """Load antibody data for a region."""
    adata_path = XENIUM_BASE / f"data/h5ad_objects/Xenium_region_{region_id}_CITE.h5ad"
    adata = sc.read_h5ad(adata_path)

    # CITE h5ad has antibody data in main matrix
    X = adata.X
    if hasattr(X, 'toarray'):
        X = X.toarray()
    marker_names = adata.var_names.tolist()
    coords = adata.obsm["spatial"]

    return X, coords, marker_names


def profile_to_celltype(profile_markers: set) -> str:
    """Match profile to best cell type."""
    best_match = "Unknown"
    best_score = 0

    for celltype, markers in IDEAL_7.items():
        overlap = len(profile_markers & markers)
        union = len(profile_markers | markers)
        jaccard = overlap / union if union > 0 else 0

        if jaccard > best_score:
            best_score = jaccard
            best_match = celltype

    return best_match if best_score > 0.3 else "Novel"


def run_parameters(top_k_values, min_score_values, regions=[0, 1, 2, 3, 4]):
    """Test different parameter combinations."""

    results = {}

    for top_k in top_k_values:
        for min_score in min_score_values:
            param_key = f"topk{top_k}_minscore{min_score}"
            print(f"\n{'='*60}")
            print(f"Testing: top_k={top_k}, min_score={min_score}")
            print(f"{'='*60}")

            all_profiles = defaultdict(list)  # profile_signature -> regions found
            region_results = {}

            for region_id in regions:
                print(f"\n  Region {region_id}...")

                try:
                    X, coords, marker_names = load_region_data(region_id)

                    # Run colocalization analysis
                    coloc_result = analyze_marker_colocalization(
                        X=X,
                        coords=coords,
                        marker_names=marker_names,
                        neighbor_k=8,
                        n_permutations=199,
                        seed=1234,
                        verbose=False,
                    )

                    # Run profile discovery with test parameters
                    profile_result = discover_profiles(
                        colocalization_result=coloc_result,
                        fdr_alpha=0.05,
                        top_k=top_k,
                        min_score=min_score,  # None = adaptive GMM
                        use_triangle_assembly=False,
                        pvalue_source="bivariate_morans",
                        verbose=False,
                    )

                    profiles = profile_result.profiles

                    # Track multi-marker profiles
                    multi_marker = []
                    for p in profiles:
                        if len(p) >= 2:
                            sig = tuple(sorted(p))
                            all_profiles[sig].append(region_id)
                            celltype = profile_to_celltype(set(p))
                            multi_marker.append((list(p), celltype))

                    region_results[region_id] = {
                        "n_profiles": len(profiles),
                        "n_multi_marker": len(multi_marker),
                        "multi_marker": multi_marker,
                    }

                    print(f"    {len(profiles)} profiles, {len(multi_marker)} multi-marker")
                    for markers, ct in multi_marker:
                        print(f"      {markers} -> {ct}")

                except Exception as e:
                    print(f"    ERROR: {e}")
                    region_results[region_id] = {"error": str(e)}

            # Analyze consistency
            consistent_profiles = {k: v for k, v in all_profiles.items() if len(v) >= 3}

            print(f"\n  CONSISTENT PROFILES (3+ regions):")
            for sig, regions_found in sorted(consistent_profiles.items(), key=lambda x: -len(x[1])):
                ct = profile_to_celltype(set(sig))
                print(f"    {list(sig)} -> {ct} (regions: {regions_found})")

            # Check ideal profile recovery
            ideal_found = {}
            for ct, markers in IDEAL_7.items():
                for sig in all_profiles.keys():
                    if markers.issubset(set(sig)) or set(sig) == markers:
                        ideal_found[ct] = all_profiles[sig]
                        break

            print(f"\n  IDEAL PROFILE RECOVERY:")
            for ct in IDEAL_7.keys():
                if ct in ideal_found:
                    print(f"    {ct}: Found in {ideal_found[ct]}")
                else:
                    print(f"    {ct}: NOT FOUND")

            results[param_key] = {
                "top_k": top_k,
                "min_score": min_score,
                "region_results": region_results,
                "consistent_profiles": {str(k): v for k, v in consistent_profiles.items()},
                "n_consistent": len(consistent_profiles),
                "ideal_recovery": {k: v for k, v in ideal_found.items()},
                "n_ideal_recovered": len(ideal_found),
            }

    return results


if __name__ == "__main__":
    # Test parameter combinations
    # Current default: top_k=6, min_score=None (adaptive)
    # Test more relaxed settings

    top_k_values = [8, 10, 15]  # More lenient (keep more edges)
    min_score_values = [None, 0.6, 0.5]  # None=adaptive, or fixed lower threshold

    print("Testing Module 2b with relaxed parameters...")
    print("Goal: Find settings that produce consistent profiles across regions")

    results = run_parameters(top_k_values, min_score_values)

    # Save results
    output_file = OUTPUT_DIR / "module2b_parameter_test.json"
    with open(output_file, "w") as f:
        json.dump(results, f, indent=2, default=str)

    print(f"\n\nResults saved to {output_file}")

    # Summary
    print("\n" + "="*60)
    print("SUMMARY: Best parameter combinations")
    print("="*60)

    best = max(results.values(), key=lambda x: (x["n_ideal_recovered"], x["n_consistent"]))
    print(f"Best: top_k={best['top_k']}, min_score={best['min_score']}")
    print(f"  Ideal profiles recovered: {best['n_ideal_recovered']}/7")
    print(f"  Consistent multi-marker: {best['n_consistent']}")
