#!/usr/bin/env python
"""
Test automatic threshold detection on all datasets.

Verifies that auto-detection:
1. Finds appropriate thresholds per sample
2. Produces reasonable profiles
3. Matches or exceeds fixed threshold performance on simulated data
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
    """Extract cell type from marker name (for simulated data)."""
    parts = marker_name.rsplit('_Protein_', 1)
    return parts[0] if len(parts) == 2 else marker_name


def test_dataset(name: str, cite_path: Path, is_simulated: bool = True):
    """Test auto-detection on a dataset."""
    print(f"\n{'=' * 70}")
    print(f"AUTO-DETECTION TEST: {name}")
    print(f"{'=' * 70}")

    if not cite_path.exists():
        print(f"  Data not found: {cite_path}")
        return None

    adata = sc.read_h5ad(cite_path)
    X = adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X
    coords = adata.obsm.get("spatial", np.random.rand(X.shape[0], 2) * 100)
    marker_names = list(adata.var_names)

    if is_simulated:
        markers_to_analyze = [m for m in marker_names if not m.startswith("Nonspecific")]
    else:
        markers_to_analyze = marker_names

    print(f"Analyzing {len(markers_to_analyze)} markers...")

    coloc_result = analyze_marker_colocalization(
        X=X, coords=coords, marker_names=marker_names,
        markers_to_analyze=markers_to_analyze, neighbor_k=6, n_permutations=199,
    )

    # Test with auto-detection (None)
    print("\n--- Auto-detection mode ---")
    result_auto = discover_hierarchical_profiles(
        coloc_result=coloc_result,
        antibody_expression=X,
        marker_names=marker_names,
        coords=coords,
        adaptive_ratio_threshold=None,  # Auto-detect!
        min_edge_score=0.0,
        verbose=True,
    )

    non_empty_auto = {k: v for k, v in result_auto.flat_profiles.items() if v}
    print(f"\nAuto-detection result: {len(non_empty_auto)} profiles")

    if is_simulated:
        # Check purity for simulated data
        pure_auto = sum(1 for markers in non_empty_auto.values()
                       if len(set(get_celltype(m) for m in markers)) == 1)
        purity_auto = pure_auto / len(non_empty_auto) if non_empty_auto else 0
        print(f"Purity: {purity_auto:.1%} ({pure_auto}/{len(non_empty_auto)} pure)")

    # Compare with fixed 0.85 threshold
    print("\n--- Fixed threshold (0.85) for comparison ---")
    result_fixed = discover_hierarchical_profiles(
        coloc_result=coloc_result,
        antibody_expression=X,
        marker_names=marker_names,
        coords=coords,
        adaptive_ratio_threshold=0.85,
        min_edge_score=0.0,
        verbose=True,
    )

    non_empty_fixed = {k: v for k, v in result_fixed.flat_profiles.items() if v}
    print(f"\nFixed 0.85 result: {len(non_empty_fixed)} profiles")

    if is_simulated:
        pure_fixed = sum(1 for markers in non_empty_fixed.values()
                        if len(set(get_celltype(m) for m in markers)) == 1)
        purity_fixed = pure_fixed / len(non_empty_fixed) if non_empty_fixed else 0
        print(f"Purity: {purity_fixed:.1%} ({pure_fixed}/{len(non_empty_fixed)} pure)")

        # Check if auto-detection meets or exceeds fixed
        if purity_auto >= purity_fixed:
            print("\n✓ Auto-detection matches or exceeds fixed threshold!")
            return True
        else:
            print(f"\n✗ Auto-detection ({purity_auto:.1%}) < Fixed ({purity_fixed:.1%})")
            return False

    # For non-simulated, just show profiles
    print("\nProfiles discovered (auto):")
    for pname, markers in sorted(non_empty_auto.items()):
        if len(markers) > 1:
            print(f"  {pname}: {markers}")

    return True


def main():
    xenium_base = Path("/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_pseudovisium/data_granular_gt/h5ad_objects")

    datasets = [
        ("High-Seg Simulated",
         Path("/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/replicates/high_seg/h5ad_objects/Wu_rep_1_CITE.h5ad"),
         True),
        ("Mixed Simulated",
         Path("/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/replicates/mixed/h5ad_objects/Wu_rep_1_CITE.h5ad"),
         True),
        ("Xenium Region 0",
         xenium_base / "Xenium_region_0_CITE.h5ad",
         False),
        ("Xenium Region 3",
         xenium_base / "Xenium_region_3_CITE.h5ad",
         False),
    ]

    results = []
    for name, path, is_sim in datasets:
        try:
            passed = test_dataset(name, path, is_sim)
            results.append((name, passed))
        except Exception as e:
            print(f"\nError: {e}")
            import traceback
            traceback.print_exc()
            results.append((name, False))

    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    for name, passed in results:
        status = "PASSED" if passed else "FAILED" if passed is False else "N/A"
        print(f"  {name}: {status}")


if __name__ == "__main__":
    main()
