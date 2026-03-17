#!/usr/bin/env python
"""
Verify the adaptive ratio threshold produces correct results.

Tests the new adaptive_ratio_threshold parameter which automatically
filters weak cross-celltype bridges by requiring edges to score at least
85% of each marker's best partnership score.
"""

import sys
from pathlib import Path
import numpy as np
import scanpy as sc

REPO_ROOT = Path("/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/.worktrees/hierarchical-profiles")
sys.path.insert(0, str(REPO_ROOT))

from CITEgeist.model.spatial_colocalization import (
    analyze_marker_colocalization,
    discover_hierarchical_profiles,
)


def get_celltype(marker_name: str) -> str:
    parts = marker_name.rsplit('_Protein_', 1)
    return parts[0] if len(parts) == 2 else marker_name


def test_dataset(name: str, cite_path: Path):
    """Test that the adaptive threshold produces pure components."""
    print(f"\n{'=' * 60}")
    print(f"TESTING: {name}")
    print(f"{'=' * 60}")

    if not cite_path.exists():
        print(f"  Data not found")
        return False

    adata = sc.read_h5ad(cite_path)
    X = adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X
    coords = adata.obsm.get("spatial", np.random.rand(X.shape[0], 2) * 100)
    marker_names = list(adata.var_names)
    markers_to_analyze = [m for m in marker_names if not m.startswith("Nonspecific")]

    print(f"Running colocalization on {len(markers_to_analyze)} markers...")
    coloc_result = analyze_marker_colocalization(
        X=X, coords=coords, marker_names=marker_names,
        markers_to_analyze=markers_to_analyze, neighbor_k=6, n_permutations=199,
    )

    # Test with adaptive threshold (default 0.85)
    print("\nTest 1: Adaptive ratio threshold = 0.85 (default)")
    result = discover_hierarchical_profiles(
        coloc_result=coloc_result,
        antibody_expression=X,
        marker_names=marker_names,
        coords=coords,
        adaptive_ratio_threshold=0.85,  # Default adaptive threshold
        min_edge_score=0.0,  # No fixed threshold
        verbose=True,
    )

    non_empty = {k: v for k, v in result.flat_profiles.items() if v}
    print(f"\nDiscovered {len(non_empty)} non-empty profiles:")

    all_pure = True
    for profile_name, markers in sorted(non_empty.items()):
        celltypes = set(get_celltype(m) for m in markers)
        is_pure = len(celltypes) == 1
        status = "PURE" if is_pure else f"MIXED: {celltypes}"
        print(f"  {profile_name}: {markers} - {status}")
        if not is_pure:
            all_pure = False

    # Also test with ratio=0.80 for comparison
    print("\nTest 2: Adaptive ratio threshold = 0.80 (more permissive)")
    result_80 = discover_hierarchical_profiles(
        coloc_result=coloc_result,
        antibody_expression=X,
        marker_names=marker_names,
        coords=coords,
        adaptive_ratio_threshold=0.80,
        min_edge_score=0.0,
        verbose=True,
    )

    non_empty_80 = {k: v for k, v in result_80.flat_profiles.items() if v}
    pure_80 = all(
        len(set(get_celltype(m) for m in markers)) == 1
        for markers in non_empty_80.values()
    )
    print(f"  Ratio 0.80: {len(non_empty_80)} profiles, all pure: {pure_80}")

    return all_pure


def main():
    datasets = [
        ("High-Seg Simulated",
         Path("/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/replicates/high_seg/h5ad_objects/Wu_rep_1_CITE.h5ad")),
        ("Mixed Simulated",
         Path("/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/replicates/mixed/h5ad_objects/Wu_rep_1_CITE.h5ad")),
    ]

    results = []
    for name, path in datasets:
        try:
            passed = test_dataset(name, path)
            results.append((name, passed))
        except Exception as e:
            print(f"\nError: {e}")
            import traceback
            traceback.print_exc()
            results.append((name, False))

    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    for name, passed in results:
        status = "PASSED" if passed else "FAILED"
        print(f"  {name}: {status}")

    all_passed = all(p for _, p in results)
    return 0 if all_passed else 1


if __name__ == "__main__":
    sys.exit(main())
