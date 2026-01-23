#!/usr/bin/env python
"""
Verify that the connected component fix produces correct trees.

Test cases:
1. Simulated high_seg data: Should produce 9 separate cell type trees
2. Xenium data: Should produce reasonable lineage separation
"""

import sys
from pathlib import Path

import numpy as np
import scanpy as sc

# Add CITEgeist to path
REPO_ROOT = Path("/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/.worktrees/hierarchical-profiles")
sys.path.insert(0, str(REPO_ROOT))

from CITEgeist.model.spatial_colocalization import (
    analyze_marker_colocalization,
    discover_hierarchical_profiles,
)


def test_simulated_data():
    """Test on simulated high_seg data - should get 9 cell types."""
    print("=" * 70)
    print("TEST 1: Simulated High-Seg Data")
    print("=" * 70)

    # Load simulated data (in main repo, not worktree)
    sim_dir = Path("/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/replicates/high_seg/h5ad_objects")
    cite_path = sim_dir / "Wu_rep_1_CITE.h5ad"

    if not cite_path.exists():
        print(f"Simulated data not found at {cite_path}")
        return False

    adata = sc.read_h5ad(cite_path)
    X = adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X
    coords = adata.obsm.get("spatial", np.random.rand(X.shape[0], 2) * 100)
    marker_names = list(adata.var_names)

    print(f"Loaded {len(marker_names)} markers, {X.shape[0]} spots")

    # Filter to cell-type specific markers only
    specific_markers = [m for m in marker_names if not m.startswith("Nonspecific")]
    print(f"Specific markers ({len(specific_markers)}): {specific_markers}")

    # Run colocalization
    coloc_result = analyze_marker_colocalization(
        X=X,
        coords=coords,
        marker_names=marker_names,
        markers_to_analyze=specific_markers,
        neighbor_k=6,
        n_permutations=199,
    )

    # Run hierarchical discovery with FDR + top-k filtering
    result = discover_hierarchical_profiles(
        coloc_result=coloc_result,
        antibody_expression=X,
        marker_names=marker_names,
        coords=coords,
        fdr_alpha=0.05,
        top_k=3,
        verbose=True,
    )

    print(f"\nDiscovered {len(result.flat_profiles)} profiles:")
    for name, markers in result.flat_profiles.items():
        print(f"  {name}: {markers}")

    # Check: we should have at least 9 separate profiles for the 9 cell types
    # In high_seg, each cell type has 2 markers that should colocalize
    n_profiles = len(result.flat_profiles)
    print(f"\nExpected: ~9 profiles (9 cell types), Got: {n_profiles}")

    success = n_profiles >= 8  # Allow for some merging
    print(f"Test {'PASSED' if success else 'FAILED'}")
    return success


def test_xenium_data():
    """Test on Xenium data - should get reasonable lineage separation."""
    print("\n" + "=" * 70)
    print("TEST 2: Xenium Data")
    print("=" * 70)

    # Load Xenium data (in main repo, not worktree)
    xenium_dir = Path("/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_pseudovisium/data_granular_gt")
    cite_path = xenium_dir / "h5ad_objects" / "Xenium_region_0_CITE.h5ad"

    if not cite_path.exists():
        print(f"Xenium data not found at {cite_path}")
        return False

    adata = sc.read_h5ad(cite_path)
    X = adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X
    coords = adata.obsm.get("spatial", np.random.rand(X.shape[0], 2) * 100)
    marker_names = list(adata.var_names)

    print(f"Loaded {len(marker_names)} markers, {X.shape[0]} spots")
    print(f"Markers: {marker_names}")

    # Run colocalization
    coloc_result = analyze_marker_colocalization(
        X=X,
        coords=coords,
        marker_names=marker_names,
        markers_to_analyze=marker_names,
        neighbor_k=6,
        n_permutations=199,
    )

    # Run hierarchical discovery
    result = discover_hierarchical_profiles(
        coloc_result=coloc_result,
        antibody_expression=X,
        marker_names=marker_names,
        coords=coords,
        fdr_alpha=0.05,
        top_k=3,
        verbose=True,
    )

    print(f"\nDiscovered {len(result.flat_profiles)} profiles:")
    for name, markers in result.flat_profiles.items():
        print(f"  {name}: {markers}")

    # Check: critical markers should be in DIFFERENT profiles
    # - PanCK (epithelial) should be separate from CD68 (macrophage)
    # - CD31 (endothelial) should be separate from both
    panck_profile = None
    cd68_profile = None
    cd31_profile = None

    for name, markers in result.flat_profiles.items():
        if "PanCK" in markers:
            panck_profile = name
        if "CD68" in markers:
            cd68_profile = name
        if "CD31" in markers:
            cd31_profile = name

    print(f"\nPanCK in profile: {panck_profile}")
    print(f"CD68 in profile: {cd68_profile}")
    print(f"CD31 in profile: {cd31_profile}")

    # All three should be in different profiles
    profiles_are_different = (
        panck_profile != cd68_profile and
        panck_profile != cd31_profile and
        cd68_profile != cd31_profile
    )

    print(f"\nAll three in different profiles? {profiles_are_different}")

    if profiles_are_different:
        print("✓ SUCCESS: Critical lineage markers properly separated!")
    else:
        print("✗ FAILED: Markers incorrectly grouped together")

    return profiles_are_different


if __name__ == "__main__":
    results = []

    # Test 1: Simulated data
    try:
        results.append(("Simulated high_seg", test_simulated_data()))
    except Exception as e:
        print(f"Simulated test error: {e}")
        import traceback
        traceback.print_exc()
        results.append(("Simulated high_seg", False))

    # Test 2: Xenium data
    try:
        results.append(("Xenium region 0", test_xenium_data()))
    except Exception as e:
        print(f"Xenium test error: {e}")
        import traceback
        traceback.print_exc()
        results.append(("Xenium region 0", False))

    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    for name, passed in results:
        status = "✓ PASSED" if passed else "✗ FAILED"
        print(f"  {name}: {status}")

    all_passed = all(p for _, p in results)
    sys.exit(0 if all_passed else 1)
