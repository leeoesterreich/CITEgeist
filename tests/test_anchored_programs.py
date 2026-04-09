#!/usr/bin/env python
"""
Test Module 4: Protein-Anchored Spatial Transcriptomic Program Discovery.

Tests both approaches:
1. Legacy: discover_anchored_programs (uses raw expression + contrastive)
2. Recommended: discover_programs_from_layers (uses deconvolved layers from Module 3)

Uses simulated data from scCubed to verify:
1. Programs are discovered for each cell type anchor
2. Cross-modal validation with proteins works
3. Spatial coherence is computed correctly
4. Deconvolved layer stacking works correctly
"""

import logging
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from CITEgeist.model import (
    # Legacy approach
    discover_anchored_programs,
    # New deconvolved layer approach
    discover_programs_from_layers,
    stack_deconvolved_layers,
    extract_celltype_expression,
    # Shared
    AnchoredProgramDiscoveryResult,
    store_results_in_adata,
    # Module 4b: Bivariate program relationships
    analyze_program_relationships,
    BivariateProgramResult,
)

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def test_anchored_programs_synthetic():
    """Test Module 4 on synthetic scCubed data."""

    # Load synthetic data
    data_dir = Path(__file__).parent.parent / "replicates" / "high_seg" / "h5ad_objects"

    gex_path = data_dir / "Wu_rep_0_GEX.h5ad"
    cite_path = data_dir / "Wu_rep_0_CITE.h5ad"

    if not gex_path.exists() or not cite_path.exists():
        logger.warning(f"Synthetic data not found at {data_dir}. Skipping test.")
        return None

    logger.info(f"Loading GEX from {gex_path}")
    adata_gex = sc.read_h5ad(gex_path)

    logger.info(f"Loading CITE from {cite_path}")
    adata_cite = sc.read_h5ad(cite_path)

    logger.info(f"GEX: {adata_gex.shape[0]} spots, {adata_gex.shape[1]} genes")
    logger.info(f"CITE: {adata_cite.shape[0]} spots, {adata_cite.shape[1]} proteins")

    # Check for ground truth cell type proportions
    # In simulated data, ground truth is typically in .obs
    gt_columns = [c for c in adata_gex.obs.columns if 'proportion' in c.lower() or 'celltype' in c.lower()]
    logger.info(f"Available obs columns: {list(adata_gex.obs.columns)[:10]}...")

    # Create mock cell type proportions (Module 3 output)
    # For synthetic data, we'll use the ground truth if available
    # or create random proportions for testing
    n_spots = adata_gex.shape[0]

    # Check for ground truth proportions
    if 'ground_truth_proportions' in adata_gex.uns:
        gt_props = adata_gex.uns['ground_truth_proportions']
        cell_type_proportions = pd.DataFrame(gt_props, index=adata_gex.obs_names)
        logger.info(f"Using ground truth proportions: {list(cell_type_proportions.columns)}")
    else:
        # Create mock proportions based on protein markers
        # The simulated data has cell-type specific proteins like "B-cells_Protein_1"
        protein_names = list(adata_cite.var_names)
        cell_types = set()
        for p in protein_names:
            if '_Protein_' in p:
                ct = p.split('_Protein_')[0]
                cell_types.add(ct)

        cell_types = sorted(list(cell_types))
        if not cell_types:
            cell_types = ['CellType_A', 'CellType_B', 'CellType_C']

        logger.info(f"Inferred cell types from protein names: {cell_types}")

        # Create proportions from protein expression
        cell_type_proportions = pd.DataFrame(index=adata_gex.obs_names)
        cite_data = adata_cite.X if not hasattr(adata_cite.X, 'toarray') else adata_cite.X.toarray()

        for ct in cell_types:
            # Find proteins for this cell type
            ct_proteins = [i for i, p in enumerate(protein_names) if p.startswith(ct + '_')]
            if ct_proteins:
                # Average expression of cell type proteins as proxy for proportion
                ct_expr = cite_data[:, ct_proteins].mean(axis=1)
                cell_type_proportions[ct] = ct_expr

        # Normalize to sum to 1
        row_sums = cell_type_proportions.sum(axis=1)
        row_sums[row_sums == 0] = 1
        cell_type_proportions = cell_type_proportions.div(row_sums, axis=0)

    logger.info(f"Cell type proportions shape: {cell_type_proportions.shape}")
    logger.info(f"Cell types: {list(cell_type_proportions.columns)}")

    # Run Module 4
    logger.info("Running Module 4: discover_anchored_programs()")

    result = discover_anchored_programs(
        adata=adata_gex,
        cell_type_proportions=cell_type_proportions,
        profile_discovery_result=None,  # No Module 2 result for this test
        protein_adata=adata_cite,
        K_programs=3,  # Discover 3 programs per cell type
        lambda_spatial=0.1,
        lambda_sparsity=0.01,
        min_proportion_threshold=0.05,  # Lower threshold for test
        validate_with_proteins=True,
        random_state=42,
    )

    # Verify results
    logger.info("=" * 50)
    logger.info("RESULTS")
    logger.info("=" * 50)
    logger.info(result.summary())

    # Check that we discovered programs for at least one anchor
    assert result.n_anchors > 0, "No anchors were processed"
    assert result.total_programs > 0, "No programs were discovered"

    # Check individual anchor results
    for anchor_name, anchor_result in result.results_by_anchor.items():
        logger.info(f"\n{anchor_name}:")
        logger.info(f"  Programs: {len(anchor_result.programs)}")
        logger.info(f"  Spots used: {anchor_result.n_spots_used}")
        logger.info(f"  Reconstruction error: {anchor_result.reconstruction_error:.2f}")

        # Check W and H matrices
        assert anchor_result.W.shape[0] == len(anchor_result.gene_names)
        assert anchor_result.W.shape[1] == len(anchor_result.programs)
        assert anchor_result.H.shape[0] == len(anchor_result.programs)
        assert anchor_result.H.shape[1] == n_spots

        # Check spatial Moran's I for each program
        for prog in anchor_result.programs:
            logger.info(f"    Program {prog.program_id}: Moran's I = {prog.spatial_moran_i:.3f} (p={prog.spatial_moran_pvalue:.3f})")
            logger.info(f"      Top genes: {prog.top_genes[:5]}")

    # Store results in AnnData
    store_results_in_adata(adata_gex, result)

    # Verify storage
    assert 'anchored_programs' in adata_gex.uns
    for anchor_name in result.results_by_anchor.keys():
        key = f'X_anchored_programs_{anchor_name}'
        assert key in adata_gex.obsm, f"Missing {key} in obsm"

    logger.info("\nTest PASSED: Module 4 (legacy contrastive) works correctly on synthetic data")

    return result


def test_deconvolved_layer_approach():
    """
    Test the new deconvolved layer-based approach (Module 4 v3).

    This test:
    1. Creates mock deconvolved layers (simulating Module 3 output)
    2. Tests stack_deconvolved_layers helper
    3. Tests discover_programs_from_layers
    """
    logger.info("=" * 60)
    logger.info("Testing DECONVOLVED LAYER approach (Module 4 v3)")
    logger.info("=" * 60)

    # Load synthetic data
    data_dir = Path(__file__).parent.parent / "replicates" / "high_seg" / "h5ad_objects"

    gex_path = data_dir / "Wu_rep_0_GEX.h5ad"
    cite_path = data_dir / "Wu_rep_0_CITE.h5ad"

    if not gex_path.exists() or not cite_path.exists():
        logger.warning(f"Synthetic data not found at {data_dir}. Skipping test.")
        return None

    logger.info(f"Loading GEX from {gex_path}")
    adata_gex = sc.read_h5ad(gex_path)

    logger.info(f"Loading CITE from {cite_path}")
    adata_cite = sc.read_h5ad(cite_path)

    logger.info(f"GEX: {adata_gex.shape[0]} spots, {adata_gex.shape[1]} genes")

    # Create mock deconvolved layers (simulating Module 3 output)
    # In real usage, these would come from optimize_gene_expression()
    protein_names = list(adata_cite.var_names)
    cell_types = set()
    for p in protein_names:
        if '_Protein_' in p:
            ct = p.split('_Protein_')[0]
            cell_types.add(ct)
    cell_types = sorted(list(cell_types))

    if not cell_types:
        cell_types = ['CellType_A', 'CellType_B', 'CellType_C']

    logger.info(f"Creating mock deconvolved layers for: {cell_types}")

    # Get raw expression
    if hasattr(adata_gex.X, 'toarray'):
        X_raw = adata_gex.X.toarray()
    else:
        X_raw = np.array(adata_gex.X)

    # Create mock cell type proportions based on protein expression
    cite_data = adata_cite.X if not hasattr(adata_cite.X, 'toarray') else adata_cite.X.toarray()
    cell_type_proportions = pd.DataFrame(index=adata_gex.obs_names)

    for ct in cell_types:
        ct_proteins = [i for i, p in enumerate(protein_names) if p.startswith(ct + '_')]
        if ct_proteins:
            ct_expr = cite_data[:, ct_proteins].mean(axis=1)
            cell_type_proportions[ct] = ct_expr

    # Normalize proportions
    row_sums = cell_type_proportions.sum(axis=1)
    row_sums[row_sums == 0] = 1
    cell_type_proportions = cell_type_proportions.div(row_sums, axis=0)

    # Create mock deconvolved layers: proportion-weighted expression
    # (In real Module 3, this comes from Gurobi optimization)
    for ct in cell_types:
        ct_props = cell_type_proportions[ct].values.reshape(-1, 1)
        # Deconvolved expression = raw * proportion (simplified mock)
        layer_data = X_raw * ct_props
        layer_name = f"{ct}_genes_pass1"
        adata_gex.layers[layer_name] = layer_data
        logger.info(f"  Created layer: {layer_name}")

    # Test 1: stack_deconvolved_layers
    logger.info("\n--- Test 1: stack_deconvolved_layers ---")
    stacked_adata = stack_deconvolved_layers(
        adata_gex,
        layer_pattern="_genes_pass1",
        coord_key="spatial",
    )

    n_spots = adata_gex.shape[0]
    n_celltypes = len(cell_types)
    expected_rows = n_spots * n_celltypes

    assert stacked_adata.shape[0] == expected_rows, \
        f"Expected {expected_rows} rows, got {stacked_adata.shape[0]}"
    assert 'original_spot' in stacked_adata.obs.columns
    assert 'cell_type' in stacked_adata.obs.columns
    assert 'spatial' in stacked_adata.obsm

    logger.info(f"  Stacked AnnData shape: {stacked_adata.shape}")
    logger.info(f"  Cell types: {stacked_adata.obs['cell_type'].unique().tolist()}")
    logger.info("  ✓ stack_deconvolved_layers works correctly")

    # Test 2: extract_celltype_expression
    logger.info("\n--- Test 2: extract_celltype_expression ---")
    for ct in cell_types[:2]:  # Test first 2 cell types
        X_ct, coords_ct = extract_celltype_expression(adata_gex, ct, "_genes_pass1")
        assert X_ct.shape == (n_spots, adata_gex.shape[1])
        assert coords_ct.shape == (n_spots, 2)
        logger.info(f"  {ct}: X shape = {X_ct.shape}, coords shape = {coords_ct.shape}")
    logger.info("  ✓ extract_celltype_expression works correctly")

    # Test 3: discover_programs_from_layers
    logger.info("\n--- Test 3: discover_programs_from_layers ---")
    result = discover_programs_from_layers(
        adata=adata_gex,
        layer_pattern="_genes_pass1",
        cell_type_proportions=cell_type_proportions,
        protein_adata=adata_cite,
        K_programs=3,
        lambda_spatial=0.1,
        lambda_sparsity=0.01,
        min_expression_threshold=0.0,
        validate_with_proteins=True,
        detect_subpopulations=True,
        n_subpopulations=3,
        random_state=42,
    )

    # Verify results
    assert result.n_anchors > 0, "No anchors were processed"
    assert result.total_programs > 0, "No programs were discovered"
    assert result.parameters.get('source') == 'deconvolved_layers', \
        "Result should indicate deconvolved_layers source"

    logger.info(f"\n  Results: {result.n_anchors} cell types, {result.total_programs} programs")

    for anchor_name, anchor_result in result.results_by_anchor.items():
        logger.info(f"\n  {anchor_name}:")
        logger.info(f"    Programs: {len(anchor_result.programs)}")
        logger.info(f"    Spots used: {anchor_result.n_spots_used}")

        # Check matrices
        assert anchor_result.W.shape[1] == len(anchor_result.programs)
        assert anchor_result.H.shape[0] == len(anchor_result.programs)

        # Check variance explained is reasonable (0-100%)
        for prog in anchor_result.programs:
            assert 0 <= prog.variance_explained <= 100, \
                f"Variance explained {prog.variance_explained} out of range"
            logger.info(
                f"    Program {prog.program_id}: "
                f"var_explained={prog.variance_explained:.1f}%, "
                f"Moran's I={prog.spatial_moran_i:.3f}"
            )

    # Store results in AnnData
    store_results_in_adata(adata_gex, result)
    assert 'anchored_programs' in adata_gex.uns
    logger.info("\n  ✓ Results stored in AnnData")

    logger.info("\n" + "=" * 60)
    logger.info("Test PASSED: Module 4 v3 (deconvolved layers) works correctly")
    logger.info("=" * 60)

    return result


def test_bivariate_program_relationships():
    """
    Test Module 4b: Bivariate Program Relationships.

    This test:
    1. Creates mock programs from deconvolved layers
    2. Runs analyze_program_relationships to find spatial cross-correlations
    3. Verifies output structure and basic statistics
    """
    logger.info("=" * 60)
    logger.info("Testing MODULE 4b: Bivariate Program Relationships")
    logger.info("=" * 60)

    # First, run Module 4 v3 to get programs
    logger.info("Step 1: Running Module 4 v3 to get programs...")

    # Load synthetic data
    data_dir = Path(__file__).parent.parent / "replicates" / "high_seg" / "h5ad_objects"

    gex_path = data_dir / "Wu_rep_0_GEX.h5ad"
    cite_path = data_dir / "Wu_rep_0_CITE.h5ad"

    if not gex_path.exists() or not cite_path.exists():
        logger.warning(f"Synthetic data not found at {data_dir}. Skipping test.")
        return None

    adata_gex = sc.read_h5ad(gex_path)
    adata_cite = sc.read_h5ad(cite_path)

    # Create mock deconvolved layers
    protein_names = list(adata_cite.var_names)
    cell_types = set()
    for p in protein_names:
        if '_Protein_' in p:
            ct = p.split('_Protein_')[0]
            cell_types.add(ct)
    cell_types = sorted(list(cell_types))

    if not cell_types:
        cell_types = ['CellType_A', 'CellType_B', 'CellType_C']

    # Get raw expression
    if hasattr(adata_gex.X, 'toarray'):
        X_raw = adata_gex.X.toarray()
    else:
        X_raw = np.array(adata_gex.X)

    # Create mock cell type proportions
    cite_data = adata_cite.X if not hasattr(adata_cite.X, 'toarray') else adata_cite.X.toarray()
    cell_type_proportions = pd.DataFrame(index=adata_gex.obs_names)

    for ct in cell_types:
        ct_proteins = [i for i, p in enumerate(protein_names) if p.startswith(ct + '_')]
        if ct_proteins:
            ct_expr = cite_data[:, ct_proteins].mean(axis=1)
            cell_type_proportions[ct] = ct_expr

    # Normalize
    row_sums = cell_type_proportions.sum(axis=1)
    row_sums[row_sums == 0] = 1
    cell_type_proportions = cell_type_proportions.div(row_sums, axis=0)

    # Create mock deconvolved layers
    for ct in cell_types:
        ct_props = cell_type_proportions[ct].values.reshape(-1, 1)
        layer_data = X_raw * ct_props
        layer_name = f"{ct}_genes_pass1"
        adata_gex.layers[layer_name] = layer_data

    # Run Module 4 v3
    result_m4 = discover_programs_from_layers(
        adata=adata_gex,
        layer_pattern="_genes_pass1",
        cell_type_proportions=cell_type_proportions,
        protein_adata=adata_cite,
        K_programs=3,
        lambda_spatial=0.1,
        lambda_sparsity=0.01,
        min_expression_threshold=0.0,
        validate_with_proteins=False,
        detect_subpopulations=False,
        random_state=42,
    )

    logger.info(f"  Module 4 discovered: {result_m4.n_anchors} cell types, {result_m4.total_programs} programs")

    # Step 2: Run Module 4b
    logger.info("\nStep 2: Running Module 4b (analyze_program_relationships)...")

    result_m4b = analyze_program_relationships(
        result=result_m4,
        adata=adata_gex,
        fdr_threshold=0.05,
        min_bivariate_i=0.1,
        n_permutations=99,  # Fewer permutations for faster test
        neighbor_k=8,
        include_within_anchor=True,
        random_state=42,
    )

    # Verify result structure
    assert isinstance(result_m4b, BivariateProgramResult), "Result should be BivariateProgramResult"
    assert result_m4b.n_pairs_tested > 0, "Should test at least one pair"
    assert len(result_m4b.all_pairs) == result_m4b.n_pairs_tested, "Pair count mismatch"

    # Print summary
    logger.info("\n" + result_m4b.summary())

    # Verify output structure
    logger.info("\nStep 3: Verifying output structure...")

    # Check that all pairs have valid relationship types
    valid_types = {"co-localized", "exclusive", "independent"}
    for pair in result_m4b.all_pairs:
        assert pair.relationship_type in valid_types, f"Invalid type: {pair.relationship_type}"
        assert -1 <= pair.bivariate_morans_i <= 1, f"Moran's I out of range: {pair.bivariate_morans_i}"
        assert 0 <= pair.bivariate_pvalue <= 1, f"P-value out of range: {pair.bivariate_pvalue}"

    logger.info(f"  ✓ All {len(result_m4b.all_pairs)} pairs have valid structure")

    # Check significant pairs
    logger.info(f"  ✓ Found {result_m4b.n_significant} significant pairs")

    # Test DataFrame export
    df = result_m4b.to_dataframe()
    assert len(df) == len(result_m4b.all_pairs), "DataFrame row count mismatch"
    expected_cols = ['anchor1', 'program1_id', 'anchor2', 'program2_id',
                     'bivariate_morans_i', 'bivariate_pvalue', 'relationship_type']
    for col in expected_cols:
        assert col in df.columns, f"Missing column: {col}"
    logger.info("  ✓ DataFrame export works correctly")

    logger.info("\n" + "=" * 60)
    logger.info("Test PASSED: Module 4b (bivariate relationships) works correctly")
    logger.info("=" * 60)

    return result_m4b


if __name__ == '__main__':
    # Test legacy approach
    result1 = test_anchored_programs_synthetic()
    if result1 is None:
        print("Legacy test skipped - synthetic data not found")
    else:
        print(f"\nLegacy: {result1.n_anchors} anchors, {result1.total_programs} programs")

    print("\n" + "-" * 60 + "\n")

    # Test new deconvolved layer approach
    result2 = test_deconvolved_layer_approach()
    if result2 is None:
        print("Deconvolved layer test skipped - synthetic data not found")
    else:
        print(f"\nDeconvolved layers: {result2.n_anchors} cell types, {result2.total_programs} programs")

    print("\n" + "-" * 60 + "\n")

    # Test Module 4b: Bivariate program relationships
    result3 = test_bivariate_program_relationships()
    if result3 is None:
        print("Module 4b test skipped - synthetic data not found")
    else:
        print(f"\nModule 4b: {result3.n_pairs_tested} pairs tested, {result3.n_significant} significant")
