"""
Test script for Module 2c v2: Coverage-Based Profile Selection.

Compares the old variance-based selection with the new coverage-based selection
that prioritizes:
1. Rarity-weighted coverage (upweights low-abundance markers)
2. Spatial novelty (different spatial patterns)
3. Non-redundancy (avoids overlapping profiles)

Also tests the two-support attachment in Module 2b.
"""

import numpy as np
import h5py
import scipy.sparse as sp
import pandas as pd
import scanpy as sc
import logging
import sys

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(message)s')
logger = logging.getLogger(__name__)

# Add project root to path
sys.path.insert(0, '/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist')

from CITEgeist.model import (
    identify_interesting_markers,
    analyze_marker_colocalization,
    discover_profiles,
    select_profiles_by_reconstruction,
)
from CITEgeist.model.spatial_colocalization import (
    select_profiles_by_coverage,
    select_profiles_by_spatial_variance,
)

# Data paths
DATA_FOLDER = "/ix1/alee/LO_LAB/General/Lab_Data/20250210_CITEGeistPublicData_GEO_Alex/processed_files/"
SIMULATED_DATA = "/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/replicates/high_seg/h5ad_objects/Wu_rep_1_CITE.h5ad"


def load_real_patient_data(sample: str):
    """Load antibody data from a real patient sample."""
    h5_path = f'{DATA_FOLDER}/{sample}/outs/filtered_feature_bc_matrix.h5'
    spatial_path = f'{DATA_FOLDER}/{sample}/outs/spatial/tissue_positions.csv'

    with h5py.File(h5_path, 'r') as f:
        feature_names = f['matrix/features/name'][:].astype(str)
        feature_types = f['matrix/features/feature_type'][:].astype(str)
        data = f['matrix/data'][:]
        indices = f['matrix/indices'][:]
        indptr = f['matrix/indptr'][:]
        shape = f['matrix/shape'][:]
        barcodes = f['matrix/barcodes'][:].astype(str)
        X_sparse = sp.csc_matrix((data, indices, indptr), shape=shape).T.tocsr()

    antibody_mask = feature_types == 'Antibody Capture'
    antibody_X = X_sparse[:, antibody_mask].toarray()
    antibody_names = feature_names[antibody_mask]

    try:
        spatial_df = pd.read_csv(spatial_path, header=0)
        if 'barcode' in spatial_df.columns:
            spatial_df = spatial_df.set_index('barcode')
    except:
        spatial_df = pd.read_csv(spatial_path, header=None, index_col=0)
        spatial_df.columns = ['in_tissue', 'array_row', 'array_col', 'pxl_row_in_fullres', 'pxl_col_in_fullres']

    common = [b for b in barcodes if b in spatial_df.index]
    idx = [dict(zip(barcodes, range(len(barcodes))))[b] for b in common]

    if 'pxl_row_in_fullres' in spatial_df.columns:
        coords = spatial_df.loc[common, ['pxl_row_in_fullres', 'pxl_col_in_fullres']].values.astype(float)
    else:
        coords = spatial_df.iloc[:, 1:3].values.astype(float)

    return antibody_X[idx], coords, list(antibody_names)


def test_simulated_data():
    """Test coverage-based selection on simulated data."""
    logger.info("=" * 70)
    logger.info("TEST: SIMULATED DATA - Coverage vs Variance Selection")
    logger.info("=" * 70)

    # Load simulated data
    adata = sc.read_h5ad(SIMULATED_DATA)
    X = adata.X if not hasattr(adata.X, 'toarray') else adata.X.toarray()
    coords = adata.obsm['spatial']
    marker_names = list(adata.var_names)

    logger.info(f"Loaded: {X.shape[0]} spots, {len(marker_names)} markers")

    # Ground truth
    gt_cell_types = ['B-cells', 'CAFs', 'Cancer Epithelial', 'Endothelial',
                     'Myeloid', 'Normal Epithelial', 'PVL', 'Plasmablasts', 'T-cells']
    gt_profiles = {ct: set([f'{ct}_Protein_1', f'{ct}_Protein_2']) for ct in gt_cell_types}
    logger.info(f"Ground truth: {len(gt_cell_types)} cell types")

    # Module 1: Identify interesting markers
    logger.info("\nRunning Module 1 (marker interest)...")
    result1 = identify_interesting_markers(X, coords, marker_names, verbose=False)
    logger.info(f"Interesting markers: {len(result1.interesting_markers)}")

    # Module 2a: Colocalization analysis (bivariate Moran's I with permutation test)
    logger.info("\nRunning Module 2a (colocalization, bivariate Moran's I)...")
    coloc_result = analyze_marker_colocalization(
        X, coords, marker_names,
        markers_to_analyze=result1.interesting_markers,
        n_permutations=199,  # Perms for bivariate Moran's I p-value
        verbose=False
    )
    logger.info(f"Pairs analyzed: {len(coloc_result.pairs)}")

    # Module 2b: Profile discovery (with two-support attachment)
    logger.info("\nRunning Module 2b (profile discovery with two-support)...")
    profile_result = discover_profiles(
        coloc_result,
        fdr_alpha=0.05,
        top_k=3,
        use_triangle_assembly=True,  # Use triangle assembly to test two-support
        verbose=True
    )
    logger.info(f"Discovered {len(profile_result.profiles)} profiles")

    # ==========================================================================
    # Module 2c: COMPARE OLD vs NEW selection methods
    # ==========================================================================

    logger.info("\n" + "=" * 50)
    logger.info("OLD METHOD: Variance-based selection")
    logger.info("=" * 50)

    selection_old = select_profiles_by_reconstruction(
        X, marker_names, profile_result.profiles,
        colocalization_result=coloc_result,
        interesting_markers=result1.interesting_markers,
        verbose=True
    )

    # Check ground truth recovery
    exact_old = sum(1 for p in selection_old.selected_profiles if set(p) in gt_profiles.values())
    logger.info(f"\nGround truth recovery (old): {exact_old}/{len(gt_cell_types)}")

    logger.info("\n" + "=" * 50)
    logger.info("NEW METHOD: Coverage-based selection")
    logger.info("=" * 50)

    selection_new = select_profiles_by_coverage(
        X, coords, marker_names, profile_result.profiles,
        colocalization_result=coloc_result,
        interesting_markers=result1.interesting_markers,
        alpha=1.0,           # Coverage weight
        beta=0.5,            # Spatial novelty weight
        gamma=0.3,           # Redundancy penalty
        min_var_explained=0.90,   # Dual checkpoint: variance threshold
        min_coverage=0.90,        # Dual checkpoint: marker coverage threshold
        n_null_samples=50,        # Samples for null distribution
        stat_significance=0.05,   # P-value threshold for stopping
        verbose=True
    )

    # Check ground truth recovery
    exact_new = sum(1 for p in selection_new.selected_profiles if set(p) in gt_profiles.values())
    logger.info(f"\nGround truth recovery (new): {exact_new}/{len(gt_cell_types)}")

    logger.info("\n" + "=" * 50)
    logger.info("SPATIAL VARIANCE METHOD: Eigenvalue-based selection")
    logger.info("=" * 50)

    selection_spatial = select_profiles_by_spatial_variance(
        X, coords, marker_names, profile_result.profiles,
        colocalization_result=coloc_result,
        interesting_markers=result1.interesting_markers,
        alpha=1.0,                 # Spatial coverage weight
        beta=0.3,                  # Proportion smoothness weight
        gamma=0.2,                 # Redundancy penalty
        min_spatial_explained=0.90,
        min_marginal_gain=0.02,    # Lower threshold for more profiles
        verbose=True
    )

    # Check ground truth recovery
    exact_spatial = sum(1 for p in selection_spatial.selected_profiles if set(p) in gt_profiles.values())
    logger.info(f"\nGround truth recovery (spatial): {exact_spatial}/{len(gt_cell_types)}")

    # ==========================================================================
    # COMPARISON
    # ==========================================================================

    logger.info("\n" + "=" * 50)
    logger.info("COMPARISON: Old vs Coverage vs Spatial")
    logger.info("=" * 50)

    logger.info(f"\n{'Metric':<25} {'Variance':<15} {'Coverage':<15} {'Spatial':<15}")
    logger.info("-" * 70)
    logger.info(f"{'Profiles selected':<25} {str(selection_old.optimal_n):<15} {str(selection_new.optimal_n):<15} {str(selection_spatial.optimal_n):<15}")
    old_ve = f"{selection_old.variance_explained[-1]:.1%}"
    new_ve = f"{selection_new.variance_explained[-1]:.1%}"
    spatial_ve = f"{selection_spatial.variance_explained[-1]:.1%}"
    logger.info(f"{'Data variance explained':<25} {old_ve:<15} {new_ve:<15} {spatial_ve:<15}")
    spatial_sv = f"{selection_spatial.explained_spatial_variance[-1]:.1%}"
    logger.info(f"{'Spatial variance expl.':<25} {'N/A':<15} {'N/A':<15} {spatial_sv:<15}")
    old_gt = f"{exact_old}/{len(gt_cell_types)}"
    new_gt = f"{exact_new}/{len(gt_cell_types)}"
    spatial_gt = f"{exact_spatial}/{len(gt_cell_types)}"
    logger.info(f"{'Ground truth recovery':<25} {old_gt:<15} {new_gt:<15} {spatial_gt:<15}")
    smooth = f"{selection_spatial.proportion_smoothness[-1]:.3f}"
    logger.info(f"{'Proportion smoothness':<25} {'N/A':<15} {'N/A':<15} {smooth:<15}")
    logger.info(f"{'Stopping reason':<25} {'var_threshold':<15} {selection_new.stopping_reason:<15} {selection_spatial.stopping_reason:<15}")

    # Show which cell types are recovered by each method
    logger.info("\nGround truth profile recovery details:")
    for ct in gt_cell_types:
        gt_markers = gt_profiles[ct]
        old_match = "✓" if any(set(p) == gt_markers for p in selection_old.selected_profiles) else "✗"
        new_match = "✓" if any(set(p) == gt_markers for p in selection_new.selected_profiles) else "✗"
        spatial_match = "✓" if any(set(p) == gt_markers for p in selection_spatial.selected_profiles) else "✗"
        logger.info(f"  {ct:<25} Var: {old_match}  Cov: {new_match}  Spatial: {spatial_match}")

    return selection_old, selection_new, selection_spatial


def test_real_patient_data():
    """Test coverage-based selection on real patient data."""
    logger.info("\n" + "=" * 70)
    logger.info("TEST: REAL PATIENT DATA - Coverage vs Variance Selection")
    logger.info("=" * 70)

    sample = "HCC22-088-P1-S2"
    X, coords, marker_names = load_real_patient_data(sample)
    logger.info(f"Loaded: {X.shape[0]} spots, {len(marker_names)} markers")

    # Module 1: Identify interesting markers
    logger.info("\nRunning Module 1 (marker interest)...")
    result1 = identify_interesting_markers(X, coords, marker_names, verbose=False)
    logger.info(f"Interesting markers: {len(result1.interesting_markers)}")

    # Check for stromal markers
    stromal_markers = ['ACTA2-1', 'VIM-1', 'COL1A1-1']
    found_stromal = [m for m in stromal_markers if m in result1.interesting_markers]
    logger.info(f"Stromal markers in interesting set: {found_stromal}")

    # Module 2a: Colocalization analysis (bivariate Moran's I with permutation test)
    logger.info("\nRunning Module 2a (colocalization, bivariate Moran's I)...")
    coloc_result = analyze_marker_colocalization(
        X, coords, marker_names,
        markers_to_analyze=result1.interesting_markers,
        n_permutations=199,  # Perms for bivariate Moran's I p-value
        verbose=False
    )
    logger.info(f"Pairs analyzed: {len(coloc_result.pairs)}")

    # Module 2b: Profile discovery
    logger.info("\nRunning Module 2b (profile discovery with two-support)...")
    profile_result = discover_profiles(
        coloc_result,
        fdr_alpha=0.05,
        top_k=3,
        use_triangle_assembly=True,
        verbose=True
    )
    logger.info(f"Discovered {len(profile_result.profiles)} profiles")

    # ==========================================================================
    # Module 2c: COMPARE OLD vs NEW selection methods
    # ==========================================================================

    logger.info("\n" + "=" * 50)
    logger.info("OLD METHOD: Variance-based selection")
    logger.info("=" * 50)

    selection_old = select_profiles_by_reconstruction(
        X, marker_names, profile_result.profiles,
        colocalization_result=coloc_result,
        interesting_markers=result1.interesting_markers,
        verbose=True
    )

    logger.info("\n" + "=" * 50)
    logger.info("NEW METHOD: Coverage-based selection")
    logger.info("=" * 50)

    selection_new = select_profiles_by_coverage(
        X, coords, marker_names, profile_result.profiles,
        colocalization_result=coloc_result,
        interesting_markers=result1.interesting_markers,
        alpha=1.0,
        beta=0.5,
        gamma=0.3,
        min_var_explained=0.90,
        min_coverage=0.90,
        n_null_samples=50,
        stat_significance=0.05,
        verbose=True
    )

    logger.info("\n" + "=" * 50)
    logger.info("SPATIAL VARIANCE METHOD: Eigenvalue-based selection")
    logger.info("=" * 50)

    selection_spatial = select_profiles_by_spatial_variance(
        X, coords, marker_names, profile_result.profiles,
        colocalization_result=coloc_result,
        interesting_markers=result1.interesting_markers,
        alpha=1.0,
        beta=0.3,
        gamma=0.2,
        min_spatial_explained=0.90,
        min_marginal_gain=0.02,    # Lower threshold for more profiles
        verbose=True
    )

    # ==========================================================================
    # COMPARISON
    # ==========================================================================

    logger.info("\n" + "=" * 50)
    logger.info("COMPARISON: Old vs Coverage vs Spatial")
    logger.info("=" * 50)

    logger.info(f"\n{'Metric':<25} {'Variance':<15} {'Coverage':<15} {'Spatial':<15}")
    logger.info("-" * 70)
    logger.info(f"{'Profiles selected':<25} {str(selection_old.optimal_n):<15} {str(selection_new.optimal_n):<15} {str(selection_spatial.optimal_n):<15}")
    old_ve = f"{selection_old.variance_explained[-1]:.1%}"
    new_ve = f"{selection_new.variance_explained[-1]:.1%}"
    spatial_ve = f"{selection_spatial.variance_explained[-1]:.1%}"
    logger.info(f"{'Data variance explained':<25} {old_ve:<15} {new_ve:<15} {spatial_ve:<15}")
    spatial_sv = f"{selection_spatial.explained_spatial_variance[-1]:.1%}"
    logger.info(f"{'Spatial variance expl.':<25} {'N/A':<15} {'N/A':<15} {spatial_sv:<15}")
    smooth = f"{selection_spatial.proportion_smoothness[-1]:.3f}"
    logger.info(f"{'Proportion smoothness':<25} {'N/A':<15} {'N/A':<15} {smooth:<15}")
    logger.info(f"{'Stopping reason':<25} {'var_threshold':<15} {selection_new.stopping_reason:<15} {selection_spatial.stopping_reason:<15}")

    # Check which method picks up stromal markers earlier
    logger.info("\nStromal marker coverage analysis:")
    stromal_markers_covered_old = []
    stromal_markers_covered_new = []

    for m in stromal_markers:
        for i, p in enumerate(selection_old.selected_profiles):
            if m in p:
                stromal_markers_covered_old.append((m, i + 1))
                break
        for i, p in enumerate(selection_new.selected_profiles):
            if m in p:
                stromal_markers_covered_new.append((m, i + 1))
                break

    logger.info(f"  Old method covered: {[m for m, _ in stromal_markers_covered_old]}")
    logger.info(f"  New method covered: {[m for m, _ in stromal_markers_covered_new]}")

    # Biological interpretation
    logger.info("\nSelected profiles (new method) with interpretation:")
    for i, profile in enumerate(selection_new.selected_profiles):
        interp = ''
        ps = set(profile)
        # Strip -1 suffix for matching
        ps_clean = set(m.replace('-1', '') for m in ps)
        if 'CD3E' in ps_clean and 'CD8A' in ps_clean:
            interp = '-> CD8+ T cells'
        elif 'CD3E' in ps_clean:
            interp = '-> T cells'
        elif 'CD68' in ps_clean and 'ITGAX' in ps_clean:
            interp = '-> Macrophages/DC'
        elif 'CD68' in ps_clean:
            interp = '-> Myeloid'
        elif 'BCL2' in ps_clean and 'SDC1' in ps_clean:
            interp = '-> Plasma cells'
        elif 'EPCAM' in ps_clean:
            interp = '-> Epithelial'
        elif 'CD19' in ps_clean:
            interp = '-> B cells'
        elif 'ACTA2' in ps_clean or 'VIM' in ps_clean:
            interp = '-> Stromal/CAFs'
        logger.info(f"  {i+1}. {profile} {interp}")

    return selection_old, selection_new, selection_spatial


def main():
    """Run all coverage selection tests."""
    logger.info("Module 2c v2: Coverage-Based Profile Selection Test")
    logger.info("=" * 70)
    logger.info("Testing two-support attachment (Module 2b) and coverage-based selection (Module 2c)")
    logger.info("=" * 70)

    # Test on simulated data
    sim_old, sim_new, sim_spatial = test_simulated_data()

    # Test on real patient data
    real_old, real_new, real_spatial = test_real_patient_data()

    # Final summary
    logger.info("\n" + "=" * 70)
    logger.info("FINAL SUMMARY")
    logger.info("=" * 70)
    logger.info("\nSimulated Data:")
    logger.info(f"  Variance method: {sim_old.optimal_n} profiles, {sim_old.variance_explained[-1]:.1%} data variance")
    logger.info(f"  Coverage method: {sim_new.optimal_n} profiles, {sim_new.variance_explained[-1]:.1%} data variance, stopped by {sim_new.stopping_reason}")
    logger.info(f"  Spatial method:  {sim_spatial.optimal_n} profiles, {sim_spatial.explained_spatial_variance[-1]:.1%} spatial variance, stopped by {sim_spatial.stopping_reason}")

    logger.info("\nReal Patient Data:")
    logger.info(f"  Variance method: {real_old.optimal_n} profiles, {real_old.variance_explained[-1]:.1%} data variance")
    logger.info(f"  Coverage method: {real_new.optimal_n} profiles, {real_new.variance_explained[-1]:.1%} data variance, stopped by {real_new.stopping_reason}")
    logger.info(f"  Spatial method:  {real_spatial.optimal_n} profiles, {real_spatial.explained_spatial_variance[-1]:.1%} spatial variance, stopped by {real_spatial.stopping_reason}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
