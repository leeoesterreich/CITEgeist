#!/usr/bin/env python
"""
Test Module 1-2c robustness across mixed and high_seg simulated data.
Settings should work without tuning across both domains.
"""
import os
import sys
import logging

os.chdir('/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist')
sys.path.insert(0, '.')

import scanpy as sc
import numpy as np

logging.basicConfig(level=logging.INFO, format='%(message)s')

from CITEgeist.model import (
    identify_interesting_markers,
    analyze_marker_colocalization,
    discover_profiles,
    select_profiles,
)

# Ground truth cell types
GT_CELL_TYPES = ['B-cells', 'CAFs', 'Cancer Epithelial', 'Endothelial',
                 'Myeloid', 'Normal Epithelial', 'PVL', 'Plasmablasts', 'T-cells']
GT_PAIRS = [(f'{ct}_Protein_1', f'{ct}_Protein_2') for ct in GT_CELL_TYPES]


def test_pipeline(dataset_name: str, h5ad_path: str):
    """Run full Module 1-2c pipeline and report metrics."""
    print(f"\n{'='*60}")
    print(f"TESTING: {dataset_name}")
    print('='*60)

    # Load data
    adata = sc.read_h5ad(h5ad_path)
    X = adata.X.toarray() if hasattr(adata.X, 'toarray') else adata.X
    coords = adata.obsm['spatial']
    marker_names = list(adata.var_names)
    print(f"Loaded: {X.shape[0]} spots, {X.shape[1]} markers")

    # ========================================
    # MODULE 1: Marker Interest Detection
    # Settings: All defaults (adaptive GMM thresholds)
    # ========================================
    print("\n--- Module 1: Marker Interest ---")
    result1 = identify_interesting_markers(
        X, coords, marker_names,
        morans_k=8,           # Default: 8 neighbors
        smooth_k=6,           # Default: 6 neighbors for smoothing
        morans_n_perm=199,    # Fast but sufficient
        gmm_snr_threshold=1.0,  # Default
        seed=1234,
        verbose=True,
    )

    # Count GT markers recovered
    gt_markers = [f'{ct}_Protein_{i}' for ct in GT_CELL_TYPES for i in [1, 2]]
    gt_recovered = [m for m in gt_markers if m in result1.interesting_markers]
    print(f"GT markers recovered: {len(gt_recovered)}/18")
    print(f"Total interesting: {len(result1.interesting_markers)}")

    # ========================================
    # MODULE 2a: Marker Colocalization
    # Settings: Use bivariate Moran's I with permutation test on smoothed data
    # ========================================
    print("\n--- Module 2a: Colocalization ---")
    coloc = analyze_marker_colocalization(
        X, coords, marker_names,
        markers_to_analyze=result1.interesting_markers,
        neighbor_k=8,          # Default: match Module 1
        smooth_k=6,            # Spatial smoothing before bivariate Moran's I
        n_permutations=999,    # Better p-value resolution
        seed=1234,
        verbose=True,
    )
    print(f"Colocalization pairs: {len(coloc.pairs)}")

    # Check GT pair scores and p-values
    print("\nGT pair bivariate Moran's I p-values:")
    gt_pair_info = []
    for a, b in GT_PAIRS:
        for p in coloc.pairs:
            if (p.marker_a == a and p.marker_b == b) or (p.marker_a == b and p.marker_b == a):
                gt_pair_info.append((a, b, p.bivariate_morans_pvalue, p.colocalization_score))
                break
        else:
            gt_pair_info.append((a, b, None, None))  # Not found

    for a, b, pval, score in gt_pair_info:
        status = 'FOUND' if pval is not None else 'MISSING'
        pval_str = f'{pval:.4f}' if pval is not None else 'N/A'
        score_str = f'{score:.3f}' if score is not None else 'N/A'
        print(f"  {a[:20]:20} - {b[:20]:20}: p={pval_str:>8}, score={score_str:>8} [{status}]")

    # ========================================
    # MODULE 2b: Profile Discovery
    # Settings: bivariate_morans FDR, top_k=3, no GMM score filter
    # ========================================
    print("\n--- Module 2b: Profile Discovery ---")
    profile_result = discover_profiles(
        coloc,
        fdr_alpha=0.05,           # Default FDR
        top_k=3,                  # Default: 3 top partners per marker
        pvalue_source='bivariate_morans',  # NEW: use bivariate Moran's I
        min_score=0.0,            # Disable GMM score filter (let FDR handle it)
        seed=1234,
        verbose=True,
    )

    # Count exact matches
    exact_matches = sum(1 for p in profile_result.profiles if set(p) in [set(pair) for pair in GT_PAIRS])
    multi_marker = sum(1 for p in profile_result.profiles if len(p) >= 2)
    print(f"\nDiscovered {len(profile_result.profiles)} profiles ({multi_marker} multi-marker)")
    print(f"GT exact matches: {exact_matches}/9")

    for i, p in enumerate(profile_result.profiles):
        match = 'EXACT' if set(p) in [set(pair) for pair in GT_PAIRS] else ''
        print(f"  {i+1}. {sorted(p)} {match}")

    # ========================================
    # MODULE 2c: Profile Selection
    # Settings: spatial variance-based
    # ========================================
    print("\n--- Module 2c: Profile Selection ---")
    selection = select_profiles(
        X, coords, marker_names,
        profiles=profile_result.profiles,
        interesting_markers=result1.interesting_markers,
        colocalization_result=coloc,
        min_spatial_explained=0.90,
        min_marginal_gain=0.005,
        verbose=True,
    )

    selected_matches = sum(1 for p in selection.selected_profiles if set(p) in [set(pair) for pair in GT_PAIRS])
    print(f"\nSelected {len(selection.selected_profiles)} profiles")

    # Handle variance_explained which may be scalar or array
    ve = selection.variance_explained
    if hasattr(ve, '__len__') and len(ve) > 0:
        var_explained = float(ve[-1])  # Get last cumulative value
    elif hasattr(ve, 'item'):
        var_explained = ve.item()  # numpy scalar
    else:
        var_explained = float(ve) if ve is not None else 0.0

    print(f"Spatial variance explained: {var_explained:.1%}")
    print(f"GT matches in selection: {selected_matches}/9")

    return {
        'dataset': dataset_name,
        'module1_gt_recovered': len(gt_recovered),
        'module1_total': len(result1.interesting_markers),
        'module2a_pairs': len(coloc.pairs),
        'module2b_profiles': len(profile_result.profiles),
        'module2b_exact_matches': exact_matches,
        'module2c_selected': len(selection.selected_profiles),
        'module2c_gt_matches': selected_matches,
        'spatial_variance': var_explained,  # Use scalar value
    }


if __name__ == '__main__':
    results = []

    # Test on mixed data
    results.append(test_pipeline(
        'MIXED (overlapping cell types)',
        'replicates/mixed/h5ad_objects/Wu_rep_1_CITE.h5ad'
    ))

    # Test on high_seg data
    results.append(test_pipeline(
        'HIGH_SEG (clear boundaries)',
        'replicates/high_seg/h5ad_objects/Wu_rep_1_CITE.h5ad'
    ))

    # Summary
    print("\n" + "="*60)
    print("SUMMARY: Settings that work across both domains")
    print("="*60)
    print("""
Module 1 (Marker Interest):
  - morans_k=8, smooth_k=6 (spatial smoothing before Moran's I)
  - Adaptive GMM thresholds for kurtosis AND Moran's I
  - OR logic: pass kurtosis OR morans, AND gmm_snr >= 1.0

Module 2a (Colocalization):
  - neighbor_k=8 (match Module 1)
  - smooth_k=6 (spatial smoothing before bivariate Moran's I)
  - n_permutations=999 (stable p-values for FDR)
  - Computes bivariate Moran's I on smoothed data with permutation test

Module 2b (Profile Discovery):
  - pvalue_source='bivariate_morans' (spatial cross-correlation)
  - fdr_alpha=0.05 (standard FDR threshold)
  - top_k=3 (mutual top-3 sparsification)
  - min_score=0.0 (no additional GMM score filter)

Module 2c (Profile Selection):
  - min_spatial_explained=0.90
  - min_marginal_gain=0.005
""")

    print("\nResults comparison:")
    print(f"{'Metric':<30} {'MIXED':>12} {'HIGH_SEG':>12}")
    print("-" * 55)
    for key in ['module1_gt_recovered', 'module2b_exact_matches', 'module2c_gt_matches', 'spatial_variance']:
        v1 = results[0][key]
        v2 = results[1][key]
        if key == 'spatial_variance':
            print(f"{key:<30} {v1:>11.1%} {v2:>11.1%}")
        else:
            print(f"{key:<30} {v1:>12} {v2:>12}")
