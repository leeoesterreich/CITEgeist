"""
Test script for Module 2: Profile Discovery from Spatial Colocalization.

This script tests the profile discovery functionality on both simulated
and real patient data with FDR correction + mutual top-k sparsification.

Key features being tested:
- BH-FDR correction (q < 0.05) for p-value filtering
- Mutual top-k sparsification to prevent hub marker collapse
- Triangle-first assembly vs hierarchical clustering
- Removed arbitrary 2-node threshold (fixes CD3E-CD8A split)
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

# Test configurations: compare hierarchical at different k values
TEST_CONFIGS = [
    {"top_k": 2, "fdr_alpha": 0.05, "use_triangle": False, "label": "hierarchical k=2"},
    {"top_k": 3, "fdr_alpha": 0.05, "use_triangle": False, "label": "hierarchical k=3"},
]

# Expected biological pairs for real data evaluation
EXPECTED_PAIRS = [
    ("EPCAM", "KRT5"),      # Epithelial
    ("CD3E", "CD8A"),       # T cells
    ("CD68", "ITGAX"),      # Myeloid
    ("SDC1", "BCL2"),       # Plasma cells
]


def count_recovered_pairs(profiles, expected_pairs=EXPECTED_PAIRS):
    """
    Count biologically expected marker pairs that were recovered in same profile.

    Args:
        profiles: List of profiles (each is a list of marker names)
        expected_pairs: List of (marker1, marker2) tuples that should be together

    Returns:
        Tuple of (recovered_count, total_expected, list of recovered pairs)
    """
    recovered = []
    for m1, m2 in expected_pairs:
        for profile in profiles:
            # Handle -1 suffix in real data marker names
            profile_set = set(profile)
            m1_match = m1 in profile_set or f"{m1}-1" in profile_set
            m2_match = m2 in profile_set or f"{m2}-1" in profile_set
            if m1_match and m2_match:
                recovered.append((m1, m2))
                break
    return len(recovered), len(expected_pairs), recovered


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


def test_simulated_data(coloc_result, config):
    """Test profile discovery on simulated data with known ground truth."""
    logger.info(f"\n--- Testing {config['label']} ---")

    # Module 2b: Profile discovery with config
    profile_result = discover_profiles(
        coloc_result,
        fdr_alpha=config['fdr_alpha'],
        top_k=config['top_k'],
        use_triangle_assembly=config.get('use_triangle', False),
        verbose=True
    )

    # Validation against ground truth
    gt_cell_types = ['B-cells', 'CAFs', 'Cancer Epithelial', 'Endothelial',
                     'Myeloid', 'Normal Epithelial', 'PVL', 'Plasmablasts', 'T-cells']
    gt_profiles = {ct: set([f'{ct}_Protein_1', f'{ct}_Protein_2']) for ct in gt_cell_types}

    exact_matches = sum(1 for p in profile_result.profiles if set(p) in gt_profiles.values())

    logger.info(f"Discovered {len(profile_result.profiles)} profiles:")
    for i, p in enumerate(profile_result.profiles):
        match = "EXACT" if set(p) in gt_profiles.values() else ""
        logger.info(f"  {i+1}. {sorted(p)} {match}")

    logger.info(f"Exact matches: {exact_matches}/9")
    logger.info(f"Modularity: {profile_result.modularity:.3f}")
    logger.info(f"Significant edges: {profile_result.n_significant_edges}")

    return {
        "config": config['label'],
        "exact_matches": exact_matches,
        "n_profiles": len(profile_result.profiles),
        "modularity": profile_result.modularity,
        "n_edges": profile_result.n_significant_edges,
        "result": profile_result
    }


def test_real_patient_data(coloc_result, config, marker_names):
    """Test profile discovery on real patient data."""
    logger.info(f"\n--- Testing {config['label']} ---")

    # Module 2b: Profile discovery with config
    profile_result = discover_profiles(
        coloc_result,
        fdr_alpha=config['fdr_alpha'],
        top_k=config['top_k'],
        use_triangle_assembly=config.get('use_triangle', False),
        verbose=True
    )

    # Count recovered pairs
    n_recovered, n_expected, recovered = count_recovered_pairs(profile_result.profiles)

    # Count multi-marker profiles
    n_multi = sum(1 for p in profile_result.profiles if len(p) >= 2)
    n_singleton = len(profile_result.profiles) - n_multi

    # Biological interpretation
    logger.info(f"Discovered {len(profile_result.profiles)} profiles ({n_multi} multi-marker, {n_singleton} singletons):")
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

    logger.info(f"Recovered pairs: {n_recovered}/{n_expected} {recovered}")
    logger.info(f"Modularity: {profile_result.modularity:.3f}")

    return {
        "config": config['label'],
        "n_profiles": len(profile_result.profiles),
        "n_multi": n_multi,
        "n_singleton": n_singleton,
        "recovered_pairs": n_recovered,
        "modularity": profile_result.modularity,
        "n_edges": profile_result.n_significant_edges,
        "result": profile_result
    }


def main():
    """Run all tests with parameter sweep."""
    logger.info("Module 2 Profile Discovery Test")
    logger.info("Testing FDR correction + mutual top-k sparsification")
    logger.info("=" * 70)

    # =================================================================
    # TEST 1: SIMULATED DATA
    # =================================================================
    logger.info("\n" + "=" * 70)
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

    # Module 2a: Colocalization analysis (run once)
    logger.info("\nRunning Module 2a (colocalization)...")
    sim_coloc = analyze_marker_colocalization(
        X, coords, marker_names,
        markers_to_analyze=result1.interesting_markers,
        n_permutations=199,
        verbose=False
    )
    logger.info(f"Pairs analyzed: {len(sim_coloc.pairs)}")

    # Test each config
    sim_results = []
    for config in TEST_CONFIGS:
        sim_results.append(test_simulated_data(sim_coloc, config))

    # =================================================================
    # TEST 2: REAL PATIENT DATA
    # =================================================================
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

    # Module 2a: Colocalization analysis (run once)
    logger.info("\nRunning Module 2a (colocalization) with 99 permutations...")
    real_coloc = analyze_marker_colocalization(
        X, coords, marker_names,
        markers_to_analyze=result1.interesting_markers,
        n_permutations=99,
        verbose=False
    )
    logger.info(f"Pairs analyzed: {len(real_coloc.pairs)}")

    # Show top colocalized pairs
    logger.info("\nTop 10 colocalized pairs:")
    for pair in real_coloc.top_pairs(10):
        logger.info(f"  {pair.marker_a:12} - {pair.marker_b:12}: score={pair.colocalization_score:.3f}")

    # Test each config
    real_results = []
    for config in TEST_CONFIGS:
        real_results.append(test_real_patient_data(real_coloc, config, marker_names))

    # =================================================================
    # SUMMARY
    # =================================================================
    logger.info("\n" + "=" * 70)
    logger.info("SUMMARY: PARAMETER SWEEP RESULTS")
    logger.info("=" * 70)

    logger.info("\nSIMULATED DATA:")
    logger.info(f"{'Config':<20} {'Exact/9':<10} {'Profiles':<10} {'Modularity':<12} {'Edges':<8}")
    logger.info("-" * 60)
    for r in sim_results:
        logger.info(f"{r['config']:<20} {r['exact_matches']}/9       {r['n_profiles']:<10} {r['modularity']:<12.3f} {r['n_edges']:<8}")

    logger.info("\nREAL PATIENT DATA:")
    logger.info(f"{'Config':<20} {'Recovered':<10} {'Profiles':<10} {'Multi':<8} {'Modularity':<12}")
    logger.info("-" * 70)
    for r in real_results:
        logger.info(f"{r['config']:<20} {r['recovered_pairs']}/4       {r['n_profiles']:<10} {r['n_multi']:<8} {r['modularity']:<12.3f}")

    # Determine best config for real data
    best_real = max(real_results, key=lambda x: (x['recovered_pairs'], x['modularity']))
    logger.info(f"\nBest real data config: {best_real['config']} (recovered {best_real['recovered_pairs']}/4 pairs)")

    return 0


if __name__ == "__main__":
    sys.exit(main())
