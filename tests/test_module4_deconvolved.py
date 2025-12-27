#!/usr/bin/env python
"""
Test Module 4 v3: Deconvolved Layer-Based Program Discovery.

Tests on:
1. Simulated data (with mock layers or actual Module 3 output)
2. Real patient data (loading CSV layers from Module 3 output)

Usage:
    python tests/test_module4_deconvolved.py --mode simulated
    python tests/test_module4_deconvolved.py --mode real --sample HCC22-088-P1-S1
    python tests/test_module4_deconvolved.py --mode both
"""

import argparse
import logging
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from CITEgeist.model import (
    discover_programs_from_layers,
    stack_deconvolved_layers,
    extract_celltype_expression,
    store_results_in_adata,
    AnchoredProgramDiscoveryResult,
)

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def load_csv_layers_to_adata(
    adata: sc.AnnData,
    layers_dir: Path,
    layer_pattern: str = "_layer_pass1.csv",
) -> sc.AnnData:
    """
    Load CSV layer files from Module 3 output into AnnData layers.

    Args:
        adata: AnnData object to add layers to
        layers_dir: Directory containing CSV layer files
        layer_pattern: Pattern to match layer files

    Returns:
        AnnData with layers added
    """
    layer_files = list(layers_dir.glob(f"*{layer_pattern}"))

    if not layer_files:
        raise ValueError(f"No layer files found in {layers_dir} matching {layer_pattern}")

    logger.info(f"Found {len(layer_files)} layer files in {layers_dir}")

    for layer_file in layer_files:
        # Extract cell type name from filename
        # e.g., "Profile_1_layer_pass1.csv" -> "Profile_1"
        cell_type = layer_file.stem.replace("_layer_pass1", "")

        # Load CSV
        df = pd.read_csv(layer_file, index_col=0)

        # Align to adata - handle potential duplicates in gene names
        common_spots = list(adata.obs_names.intersection(df.index))

        if len(common_spots) == 0:
            logger.warning(f"No matching spots for {cell_type}, skipping")
            continue

        # Create layer matrix aligned to adata
        layer_matrix = np.zeros((adata.shape[0], adata.shape[1]))

        # Build gene name to index mapping (handle duplicates)
        adata_gene_to_idx = {}
        for idx, gene in enumerate(adata.var_names):
            if gene not in adata_gene_to_idx:
                adata_gene_to_idx[gene] = idx

        # Find common genes between adata and CSV
        common_genes = [g for g in df.columns if g in adata_gene_to_idx]

        if len(common_genes) == 0:
            logger.warning(f"No matching genes for {cell_type}, skipping")
            continue

        # Get indices
        spot_idx = np.array([adata.obs_names.get_loc(s) for s in common_spots])
        gene_idx = np.array([adata_gene_to_idx[g] for g in common_genes])

        # Fill in values using advanced indexing
        layer_data = df.loc[common_spots, common_genes].values
        for i, s_idx in enumerate(spot_idx):
            layer_matrix[s_idx, gene_idx] = layer_data[i, :]

        # Add as layer with standard naming
        layer_name = f"{cell_type}_genes_pass1"
        adata.layers[layer_name] = layer_matrix

        logger.info(f"  Added layer: {layer_name} ({len(common_spots)} spots, {len(common_genes)} genes)")

    return adata


def test_simulated_data(replicate: int = 0):
    """
    Test Module 4 v3 on simulated scCube data.

    Creates mock deconvolved layers since simulated data doesn't have Module 3 output.
    """
    logger.info("=" * 70)
    logger.info("TESTING ON SIMULATED DATA (with mock deconvolved layers)")
    logger.info("=" * 70)

    data_dir = Path(__file__).parent.parent / "replicates" / "high_seg" / "h5ad_objects"

    gex_path = data_dir / f"Wu_rep_{replicate}_GEX.h5ad"
    cite_path = data_dir / f"Wu_rep_{replicate}_CITE.h5ad"

    if not gex_path.exists() or not cite_path.exists():
        logger.error(f"Simulated data not found at {data_dir}")
        return None

    # Load data
    logger.info(f"Loading GEX from {gex_path}")
    adata_gex = sc.read_h5ad(gex_path)

    logger.info(f"Loading CITE from {cite_path}")
    adata_cite = sc.read_h5ad(cite_path)

    logger.info(f"GEX: {adata_gex.shape[0]} spots, {adata_gex.shape[1]} genes")
    logger.info(f"CITE: {adata_cite.shape[0]} spots, {adata_cite.shape[1]} proteins")

    # Create mock deconvolved layers based on protein expression
    protein_names = list(adata_cite.var_names)
    cell_types = sorted(set(
        p.split('_Protein_')[0] for p in protein_names if '_Protein_' in p
    ))

    if not cell_types:
        cell_types = ['CellType_A', 'CellType_B', 'CellType_C']

    logger.info(f"Creating mock layers for {len(cell_types)} cell types: {cell_types}")

    # Get raw expression
    X_raw = adata_gex.X.toarray() if hasattr(adata_gex.X, 'toarray') else np.array(adata_gex.X)

    # Create proportions from protein expression
    cite_data = adata_cite.X.toarray() if hasattr(adata_cite.X, 'toarray') else np.array(adata_cite.X)
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

    # Create mock deconvolved layers
    for ct in cell_types:
        ct_props = cell_type_proportions[ct].values.reshape(-1, 1)
        layer_data = X_raw * ct_props
        layer_name = f"{ct}_genes_pass1"
        adata_gex.layers[layer_name] = layer_data

    logger.info(f"Created {len(cell_types)} mock deconvolved layers")

    # Run Module 4 v3
    logger.info("\nRunning discover_programs_from_layers...")

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

    # Print results
    print_results(result, "SIMULATED")

    return result


def test_real_data(sample_name: str = "HCC22-088-P1-S1"):
    """
    Test Module 4 v3 on real patient data with Module 3 output.
    """
    logger.info("=" * 70)
    logger.info(f"TESTING ON REAL DATA: {sample_name}")
    logger.info("=" * 70)

    # Paths
    data_dir = Path("/ix1/alee/LO_LAB/General/Lab_Data/20250210_CITEGeistPublicData_GEO_Alex/processed_files")
    output_dir = Path(__file__).parent.parent / "output"

    sample_data_dir = data_dir / sample_name / "outs"
    layers_dir = output_dir / f"{sample_name}_pass1" / "layers" / "pass1"

    if not sample_data_dir.exists():
        logger.error(f"Sample data not found at {sample_data_dir}")
        return None

    if not layers_dir.exists():
        logger.error(f"Module 3 layers not found at {layers_dir}")
        logger.error("Run Module 3 first to generate deconvolved layers")
        return None

    # Load raw data using squidpy
    logger.info(f"Loading sample from {sample_data_dir}")
    adata = sq.read.visium(
        str(sample_data_dir),
        counts_file='filtered_feature_bc_matrix.h5',
        load_images=True,
        gex_only=False,
    )

    logger.info(f"Loaded: {adata.shape[0]} spots, {adata.shape[1]} features")

    # Split into GEX and antibody
    feature_types = adata.var['feature_types']
    gex_mask = feature_types == 'Gene Expression'
    ab_mask = feature_types == 'Antibody Capture'

    adata_gex = adata[:, gex_mask].copy()
    adata_cite = adata[:, ab_mask].copy()

    # Make gene names unique (required for layer loading)
    adata_gex.var_names_make_unique()

    logger.info(f"GEX: {adata_gex.shape[1]} genes")
    logger.info(f"Antibody: {adata_cite.shape[1]} proteins")

    # Load CSV layers into adata
    logger.info(f"\nLoading Module 3 layers from {layers_dir}")
    adata_gex = load_csv_layers_to_adata(adata_gex, layers_dir)

    # Check layers
    layer_names = [l for l in adata_gex.layers.keys() if '_genes_pass1' in l]
    logger.info(f"Loaded {len(layer_names)} layers: {layer_names[:5]}...")

    # Load cell type proportions if available
    prop_file = output_dir / f"{sample_name}_cell_prop_finetuned_results.csv"
    cell_type_proportions = None
    if prop_file.exists():
        cell_type_proportions_raw = pd.read_csv(prop_file, index_col=0)
        # Align to adata - create full-sized DataFrame with zeros for missing spots
        cell_type_proportions = pd.DataFrame(
            0.0,
            index=adata_gex.obs_names,
            columns=cell_type_proportions_raw.columns
        )
        common = adata_gex.obs_names.intersection(cell_type_proportions_raw.index)
        cell_type_proportions.loc[common, :] = cell_type_proportions_raw.loc[common, :]
        logger.info(f"Loaded cell type proportions: {cell_type_proportions.shape} ({len(common)} spots with data)")

    # Run Module 4 v3
    logger.info("\nRunning discover_programs_from_layers...")

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

    # Print results
    print_results(result, f"REAL ({sample_name})")

    # Save results
    output_file = output_dir / f"{sample_name}_module4_v3_programs.csv"
    result.to_dataframe().to_csv(output_file)
    logger.info(f"\nSaved results to {output_file}")

    return result


def print_results(result: AnchoredProgramDiscoveryResult, data_type: str):
    """Print summary of results."""
    print("\n" + "=" * 70)
    print(f"RESULTS: {data_type} DATA")
    print("=" * 70)
    print(f"Cell types analyzed: {result.n_anchors}")
    print(f"Total programs: {result.total_programs}")
    print(f"Source: {result.parameters.get('source', 'unknown')}")

    print("\nPer-cell-type summary:")
    print("-" * 70)

    for anchor_name, anchor_result in result.results_by_anchor.items():
        print(f"\n{anchor_name}:")
        print(f"  Spots used: {anchor_result.n_spots_used}")
        print(f"  Programs: {len(anchor_result.programs)}")

        for prog in anchor_result.programs:
            print(f"    Program {prog.program_id}: "
                  f"var_explained={prog.variance_explained:.1f}%, "
                  f"Moran's I={prog.spatial_moran_i:.3f} (p={prog.spatial_moran_pvalue:.3f})")
            print(f"      Top genes: {', '.join(prog.top_genes[:5])}")

        if anchor_result.subpopulations:
            print(f"  Subpopulations: {len(anchor_result.subpopulations)}")
            for sp in anchor_result.subpopulations:
                print(f"    {sp.location_label}: {sp.n_spots} spots, dominant program {sp.dominant_program}")

    print("\n" + "=" * 70)


def main():
    parser = argparse.ArgumentParser(description="Test Module 4 v3 on real and simulated data")
    parser.add_argument(
        "--mode",
        choices=["simulated", "real", "both"],
        default="both",
        help="Which data to test on"
    )
    parser.add_argument(
        "--sample",
        default="HCC22-088-P1-S1",
        help="Sample name for real data testing"
    )
    parser.add_argument(
        "--replicate",
        type=int,
        default=0,
        help="Replicate number for simulated data (0-4)"
    )

    args = parser.parse_args()

    results = {}

    if args.mode in ["simulated", "both"]:
        results['simulated'] = test_simulated_data(args.replicate)

    if args.mode in ["real", "both"]:
        results['real'] = test_real_data(args.sample)

    # Summary
    print("\n" + "=" * 70)
    print("FINAL SUMMARY")
    print("=" * 70)

    for data_type, result in results.items():
        if result is not None:
            print(f"{data_type}: {result.n_anchors} cell types, {result.total_programs} programs")

            # Check variance explained is reasonable
            all_var = []
            for anchor in result.results_by_anchor.values():
                for prog in anchor.programs:
                    all_var.append(prog.variance_explained)

            if all_var:
                print(f"  Variance explained range: {min(all_var):.1f}% - {max(all_var):.1f}%")
                if max(all_var) > 100:
                    print("  WARNING: Variance explained > 100% - possible bug")
                else:
                    print("  OK: Variance explained in valid range")
        else:
            print(f"{data_type}: FAILED or data not found")


if __name__ == "__main__":
    main()
