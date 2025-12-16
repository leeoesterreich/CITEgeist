"""
Test script for Module 2: Profile Discovery from Spatial Colocalization.

This script tests the profile discovery functionality on both simulated
and real patient data, validating:
1. Simulated data: Should recover 9 ground truth cell type profiles
2. Real data: Should produce biologically sensible lineage groupings

Key changes being tested:
- GMM-BIC adaptive score threshold selection
- Gap-based lineage splitting for separate dendrograms
- Dynamic tree cutting within lineages
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

from CITEgeist.model import identify_interesting_markers, analyze_marker_colocalization, discover_profiles

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
    """Test profile discovery on simulated data with known ground truth."""
    logger.info("=" * 70)
    logger.info("TEST 1: SIMULATED DATA")
    logger.info("=" * 70)

    # Load simulated data
    adata = sc.read_h5ad(SIMULATED_DATA)
    X = adata.X if not hasattr(adata.X, 'toarray') else adata.X.toarray()
    coords = adata.obsm['spatial']
    marker_names = list(adata.var_names)

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
        n_permutations=199,
        verbose=False
    )
    logger.info(f"Pairs analyzed: {len(coloc_result.pairs)}")

    # Module 2b: Profile discovery
    logger.info("\nRunning Module 2b (profile discovery)...")
    profile_result = discover_profiles(coloc_result, verbose=True)

    # Validation against ground truth
    gt_cell_types = ['B-cells', 'CAFs', 'Cancer Epithelial', 'Endothelial',
                     'Myeloid', 'Normal Epithelial', 'PVL', 'Plasmablasts', 'T-cells']
    gt_profiles = {ct: set([f'{ct}_Protein_1', f'{ct}_Protein_2']) for ct in gt_cell_types}

    exact_matches = sum(1 for p in profile_result.profiles if set(p) in gt_profiles.values())

    logger.info(f"\n--- RESULTS ---")
    logger.info(f"Discovered {len(profile_result.profiles)} profiles:")
    for i, p in enumerate(profile_result.profiles):
        match = "EXACT" if set(p) in gt_profiles.values() else ""
        logger.info(f"  {i+1}. {sorted(p)} {match}")

    logger.info(f"\nExact matches: {exact_matches}/9")
    logger.info(f"Lineage dendrograms: {len(profile_result.lineage_dendrograms)}")
    logger.info(f"Modularity: {profile_result.modularity:.3f}")
    logger.info(f"Significant edges: {profile_result.n_significant_edges}")

    return exact_matches, profile_result


def test_real_patient_data():
    """Test profile discovery on real patient data."""
    logger.info("\n" + "=" * 70)
    logger.info("TEST 2: REAL PATIENT DATA (HCC22-088-P1-S2)")
    logger.info("=" * 70)

    sample = "HCC22-088-P1-S2"
    X, coords, marker_names = load_real_patient_data(sample)
    logger.info(f"Loaded: {X.shape[0]} spots, {len(marker_names)} markers")

    # Module 1: Identify interesting markers
    logger.info("\nRunning Module 1 (marker interest)...")
    result1 = identify_interesting_markers(X, coords, marker_names, verbose=False)
    logger.info(f"Interesting markers: {len(result1.interesting_markers)}")
    logger.info(f"  {list(result1.interesting_markers)}")

    # Module 2a: Colocalization analysis
    logger.info("\nRunning Module 2a (colocalization) with 99 permutations...")
    coloc_result = analyze_marker_colocalization(
        X, coords, marker_names,
        markers_to_analyze=result1.interesting_markers,
        n_permutations=99,
        verbose=False
    )
    logger.info(f"Pairs analyzed: {len(coloc_result.pairs)}")

    # Show top colocalized pairs
    logger.info("\nTop 10 colocalized pairs:")
    for pair in coloc_result.top_pairs(10):
        logger.info(f"  {pair.marker_a:12} - {pair.marker_b:12}: score={pair.colocalization_score:.3f}")

    # Module 2b: Profile discovery
    logger.info("\nRunning Module 2b (profile discovery)...")
    profile_result = discover_profiles(coloc_result, verbose=True)

    # Biological interpretation
    logger.info(f"\n--- DISCOVERED PROFILES ---")
    for i, profile in enumerate(profile_result.profiles):
        markers = ', '.join(sorted(profile))
        interp = ''
        ps = set(profile)
        if 'CD3E' in ps and 'CD8A' in ps:
            interp = '-> CD8+ T cells'
        elif 'CD3E' in ps and 'CD4' in ps:
            interp = '-> CD4+ T cells'
        elif 'CD3E' in ps:
            interp = '-> T cells'
        elif 'CD68' in ps and 'ITGAX' in ps:
            interp = '-> Macrophages/DC'
        elif 'CD68' in ps:
            interp = '-> Myeloid'
        elif 'ACTA2' in ps and 'VIM' in ps:
            interp = '-> Stromal/CAFs'
        elif 'BCL2' in ps and 'SDC1' in ps:
            interp = '-> Plasma cells'
        elif 'CD19' in ps:
            interp = '-> B cells'
        elif 'EPCAM' in ps and 'KRT5' in ps:
            interp = '-> Epithelial'
        elif 'EPCAM' in ps:
            interp = '-> Epithelial'
        elif 'CEACAM8' in ps and 'ITGAM' in ps:
            interp = '-> Neutrophils'
        logger.info(f"  {i+1}. ({len(profile)}) {markers} {interp}")

    logger.info(f"\n{profile_result.summary()}")
    logger.info(f"Lineage dendrograms: {len(profile_result.lineage_dendrograms)}")

    return profile_result


def main():
    """Run all tests."""
    logger.info("Module 2 Profile Discovery Test")
    logger.info("Testing gap-based lineage splitting and dynamic tree cutting")
    logger.info("=" * 70)

    # Test on simulated data
    exact_matches, sim_result = test_simulated_data()

    # Test on real data
    real_result = test_real_patient_data()

    # Summary
    logger.info("\n" + "=" * 70)
    logger.info("SUMMARY")
    logger.info("=" * 70)
    logger.info(f"Simulated data: {exact_matches}/9 exact matches")
    logger.info(f"Simulated modularity: {sim_result.modularity:.3f}")
    logger.info(f"Real data profiles: {len(real_result.profiles)}")
    logger.info(f"Real data lineages: {len(real_result.lineage_dendrograms)}")
    logger.info(f"Real data modularity: {real_result.modularity:.3f}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
