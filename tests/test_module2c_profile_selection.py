"""
Test script for Module 2c: Reconstruction-Based Profile Selection.

This script tests the profile selection functionality that uses
reconstruction accuracy to determine optimal number of profiles.
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
    """Test Module 2c on simulated data with known ground truth."""
    logger.info("=" * 70)
    logger.info("TEST: SIMULATED DATA")
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
    logger.info(f"Ground truth: {len(gt_cell_types)} cell types")

    # Module 1: Identify interesting markers
    logger.info("\nRunning Module 1 (marker interest)...")
    result1 = identify_interesting_markers(X, coords, marker_names, verbose=False)
    logger.info(f"Interesting markers: {len(result1.interesting_markers)}")

    # Module 2a: Colocalization analysis
    logger.info("\nRunning Module 2a (colocalization)...")
    coloc_result = analyze_marker_colocalization(
        X, coords, marker_names,
        markers_to_analyze=result1.interesting_markers,
        n_permutations=5000,
        verbose=False
    )
    logger.info(f"Pairs analyzed: {len(coloc_result.pairs)}")

    # Module 2b: Profile discovery
    logger.info("\nRunning Module 2b (profile discovery, k=3)...")
    profile_result = discover_profiles(
        coloc_result,
        fdr_alpha=0.05,
        top_k=3,
        verbose=True
    )
    logger.info(f"Discovered {len(profile_result.profiles)} profiles")

    # Module 2c: Profile selection
    logger.info("\nRunning Module 2c (reconstruction-based selection)...")
    selection = select_profiles_by_reconstruction(
        X, marker_names, profile_result.profiles,
        colocalization_result=coloc_result,
        interesting_markers=result1.interesting_markers,  # Pass interesting markers
        verbose=True
    )

    logger.info("\n" + "=" * 50)
    logger.info("MODULE 2c RESULTS (SIMULATED):")
    logger.info("=" * 50)
    logger.info(selection.summary())

    # Show reconstruction curve
    logger.info("\nReconstruction curve:")
    logger.info(f"{'n':>3} {'RMSE':>10} {'VarExpl':>10}")
    logger.info("-" * 25)
    for i, (rmse, ve) in enumerate(zip(selection.reconstruction_errors, selection.variance_explained)):
        marker = " <-- elbow" if i + 1 == selection.optimal_n else ""
        logger.info(f"{i+1:>3} {rmse:>10.4f} {ve:>9.1%}{marker}")

    # Check ground truth recovery
    gt_profiles = {ct: set([f'{ct}_Protein_1', f'{ct}_Protein_2']) for ct in gt_cell_types}
    exact_matches = sum(1 for p in selection.selected_profiles if set(p) in gt_profiles.values())
    logger.info(f"\nGround truth recovery: {exact_matches}/{len(gt_cell_types)} exact matches")

    return selection


def test_real_patient_data():
    """Test Module 2c on real patient data."""
    logger.info("\n" + "=" * 70)
    logger.info("TEST: REAL PATIENT DATA (HCC22-088-P1-S2)")
    logger.info("=" * 70)

    sample = "HCC22-088-P1-S2"
    X, coords, marker_names = load_real_patient_data(sample)
    logger.info(f"Loaded: {X.shape[0]} spots, {len(marker_names)} markers")

    # Module 1: Identify interesting markers
    logger.info("\nRunning Module 1 (marker interest)...")
    result1 = identify_interesting_markers(X, coords, marker_names, verbose=False)
    logger.info(f"Interesting markers: {len(result1.interesting_markers)}")

    # Module 2a: Colocalization analysis
    logger.info("\nRunning Module 2a (colocalization)...")
    coloc_result = analyze_marker_colocalization(
        X, coords, marker_names,
        markers_to_analyze=result1.interesting_markers,
        n_permutations=5000,
        verbose=False
    )
    logger.info(f"Pairs analyzed: {len(coloc_result.pairs)}")

    # Module 2b: Profile discovery
    logger.info("\nRunning Module 2b (profile discovery, k=3)...")
    profile_result = discover_profiles(
        coloc_result,
        fdr_alpha=0.05,
        top_k=3,
        verbose=True
    )
    logger.info(f"Discovered {len(profile_result.profiles)} profiles")

    # Module 2c: Profile selection
    logger.info("\nRunning Module 2c (reconstruction-based selection)...")
    selection = select_profiles_by_reconstruction(
        X, marker_names, profile_result.profiles,
        colocalization_result=coloc_result,
        interesting_markers=result1.interesting_markers,  # Pass interesting markers
        verbose=True
    )

    logger.info("\n" + "=" * 50)
    logger.info("MODULE 2c RESULTS (REAL DATA):")
    logger.info("=" * 50)
    logger.info(selection.summary())

    # Show reconstruction curve
    logger.info("\nReconstruction curve:")
    logger.info(f"{'n':>3} {'RMSE':>10} {'VarExpl':>10}")
    logger.info("-" * 25)
    for i, (rmse, ve) in enumerate(zip(selection.reconstruction_errors, selection.variance_explained)):
        marker = " <-- elbow" if i + 1 == selection.optimal_n else ""
        logger.info(f"{i+1:>3} {rmse:>10.4f} {ve:>9.1%}{marker}")

    # Biological interpretation of selected profiles
    logger.info("\nSelected profiles with interpretation:")
    for i, profile in enumerate(selection.selected_profiles):
        interp = ''
        ps = set(profile)
        if 'CD3E' in ps and 'CD8A' in ps:
            interp = '-> CD8+ T cells'
        elif 'CD3E' in ps:
            interp = '-> T cells'
        elif 'CD68' in ps and 'ITGAX' in ps:
            interp = '-> Macrophages/DC'
        elif 'CD68' in ps:
            interp = '-> Myeloid'
        elif 'BCL2' in ps and 'SDC1' in ps:
            interp = '-> Plasma cells'
        elif 'EPCAM' in ps:
            interp = '-> Epithelial'
        elif 'CD19' in ps:
            interp = '-> B cells'
        logger.info(f"  {i+1}. {profile} {interp}")

    return selection


def main():
    """Run all Module 2c tests."""
    logger.info("Module 2c: Reconstruction-Based Profile Selection Test")
    logger.info("=" * 70)

    # Test on simulated data
    sim_result = test_simulated_data()

    # Test on real patient data
    real_result = test_real_patient_data()

    # Summary comparison
    logger.info("\n" + "=" * 70)
    logger.info("SUMMARY")
    logger.info("=" * 70)
    logger.info(f"Simulated: Selected {sim_result.optimal_n}/{len(sim_result.all_profiles_ranked)} profiles")
    logger.info(f"  Variance explained: {sim_result.variance_explained[sim_result.optimal_n - 1]:.1%}")
    logger.info(f"Real data: Selected {real_result.optimal_n}/{len(real_result.all_profiles_ranked)} profiles")
    logger.info(f"  Variance explained: {real_result.variance_explained[real_result.optimal_n - 1]:.1%}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
