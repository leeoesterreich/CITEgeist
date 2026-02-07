#!/usr/bin/env python
"""
Test discover_joint_programs on real HCC22-088 patient data.

This script validates the joint discovery function works with actual
Module 3 deconvolved layers.
"""

import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq

# Add parent to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from CITEgeist.model import discover_joint_programs, JointDiscoveryResult

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
    """Load CSV layer files from Module 3 output into AnnData layers."""
    layer_files = list(layers_dir.glob(f"*{layer_pattern}"))

    if not layer_files:
        raise ValueError(f"No layer files found in {layers_dir} matching {layer_pattern}")

    logger.info(f"Found {len(layer_files)} layer files in {layers_dir}")

    for layer_file in layer_files:
        # Extract cell type name from filename
        cell_type = layer_file.stem.replace("_layer_pass1", "")

        # Clean up cell type name for layer key
        # Convert to simpler format for layer pattern matching
        layer_key = f"{cell_type}_genes_pass1"

        # Load CSV
        df = pd.read_csv(layer_file, index_col=0)

        # Align to adata
        common_spots = list(adata.obs_names.intersection(df.index))

        if len(common_spots) == 0:
            logger.warning(f"No matching spots for {cell_type}, skipping")
            continue

        # Create layer matrix aligned to adata
        layer_matrix = np.zeros((adata.shape[0], adata.shape[1]))

        # Build gene name to index mapping
        adata_gene_to_idx = {gene: idx for idx, gene in enumerate(adata.var_names)}

        # Find common genes
        common_genes = [g for g in df.columns if g in adata_gene_to_idx]

        if len(common_genes) == 0:
            logger.warning(f"No matching genes for {cell_type}, skipping")
            continue

        # Get indices
        spot_idx = np.array([list(adata.obs_names).index(s) for s in common_spots])
        gene_idx = np.array([adata_gene_to_idx[g] for g in common_genes])

        # Fill values
        csv_spot_idx = np.array([list(df.index).index(s) for s in common_spots])
        csv_gene_idx = np.array([list(df.columns).index(g) for g in common_genes])

        for i, (si, csv_si) in enumerate(zip(spot_idx, csv_spot_idx)):
            for j, (gi, csv_gi) in enumerate(zip(gene_idx, csv_gene_idx)):
                layer_matrix[si, gi] = df.iloc[csv_si, csv_gi]

        adata.layers[layer_key] = layer_matrix
        logger.info(f"  Loaded {cell_type}: {len(common_spots)} spots, {len(common_genes)} genes")

    return adata


def main():
    sample_name = "HCC22-088-P1-S1"

    # Paths
    base_dir = Path("/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeistNeilAnalysis/CITEgeist")
    data_dir = Path("/ix1/alee/LO_LAB/General/Lab_Data/20250210_CITEGeistPublicData_GEO_Alex/processed_files")
    output_dir = base_dir / "output"
    layers_dir = output_dir / f"{sample_name}_pass1" / "layers" / "pass1"

    logger.info(f"Testing discover_joint_programs on {sample_name}")

    # Load the sample data
    sample_path = data_dir / sample_name / "outs"
    logger.info(f"Loading sample from {sample_path}")

    adata = sq.read.visium(
        sample_path,
        counts_file='filtered_feature_bc_matrix.h5',
        load_images=True,
        gex_only=False
    )
    logger.info(f"Loaded adata: {adata.shape}")

    # Split by feature type (GEX only)
    feature_types = adata.var['feature_types']
    gex_mask = feature_types == 'Gene Expression'
    adata_gex = adata[:, gex_mask].copy()
    logger.info(f"GEX subset: {adata_gex.shape}")

    # Load cell proportions
    props_file = output_dir / f"{sample_name}_cell_prop_finetuned_results.csv"
    if not props_file.exists():
        props_file = output_dir / f"{sample_name}_cell_prop_global_results.csv"

    proportions = pd.read_csv(props_file, index_col=0)
    logger.info(f"Loaded proportions: {proportions.shape}")
    logger.info(f"Cell types: {list(proportions.columns)}")

    # Load deconvolved layers
    logger.info(f"Loading layers from {layers_dir}")
    adata_gex = load_csv_layers_to_adata(adata_gex, layers_dir)

    logger.info(f"Available layers: {list(adata_gex.layers.keys())}")

    # Run joint discovery
    logger.info("Running discover_joint_programs...")

    result = discover_joint_programs(
        adata=adata_gex,
        cell_type_proportions=proportions,
        K_programs=10,
        layer_pattern="_genes_pass1",
        random_state=42,
    )

    # Report results
    logger.info("\n" + "="*60)
    logger.info("RESULTS")
    logger.info("="*60)
    logger.info(result.summary())

    # Show programs by type
    logger.info("\nPrograms by type:")
    for prog in result.programs:
        logger.info(
            f"  Program {prog.program_id}: {prog.program_type} "
            f"(primary={prog.primary_cell_type}, "
            f"interaction_score={prog.interaction_score:.2f}, "
            f"Moran's I={prog.spatial_moran_i:.3f})"
        )
        logger.info(f"    Top genes: {', '.join(prog.top_genes[:5])}")

    # Save results
    output_file = output_dir / f"{sample_name}_joint_programs.csv"
    result.to_dataframe().to_csv(output_file, index=False)
    logger.info(f"\nSaved results to {output_file}")

    logger.info("\nTest completed successfully!")

    return result


if __name__ == "__main__":
    main()
