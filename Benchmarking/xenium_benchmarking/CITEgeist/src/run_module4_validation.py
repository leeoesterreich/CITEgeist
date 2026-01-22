#!/usr/bin/env python
"""
Module 4/5 Validation on Xenium Data.

Validates Module 4 (Anchored Program Discovery) and Module 5 (Cross-Sample
Integration) by running them on three data sources:
1. Single-cell Xenium RCC data (cell-level)
2. Manual deconvolved data (output_achievable_7)
3. Autodiscovery deconvolved data (output_autodiscovery)

This tests both functionality (code runs without errors) and enables
concordance validation (compare programs across resolutions).

Usage:
    python run_module4_validation.py --source singlecell --region-id 0
    python run_module4_validation.py --source manual --region-id 0
    python run_module4_validation.py --source autodiscovery --region-id 0
    python run_module4_validation.py --source all --all-regions
"""

import argparse
import json
import logging
import os
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import scanpy as sc

# Add paths
REPO_ROOT = Path(__file__).parent.parent.parent.parent.parent
sys.path.insert(0, str(REPO_ROOT / "CITEgeist"))
sys.path.insert(0, str(REPO_ROOT / "Benchmarking" / "xenium_pseudovisium" / "src"))
sys.path.insert(0, str(REPO_ROOT / "Benchmarking" / "xenium_benchmarking" / "CITEgeist" / "src"))

from model.anchored_program_discovery import (
    discover_programs_from_layers,
    analyze_program_relationships,
    AnchoredProgramDiscoveryResult,
    BivariateProgramResult,
)
from model.cross_sample_integration import (
    integrate_samples,
    IntegrationResult,
)

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Directories
BENCH_DIR = REPO_ROOT / "Benchmarking" / "xenium_benchmarking" / "CITEgeist"
PSEUDOVISIUM_DIR = REPO_ROOT / "Benchmarking" / "xenium_pseudovisium"
OUTPUT_DIR = BENCH_DIR / "output_module4_validation"

# Data source paths
DATA_SOURCES = {
    "manual": BENCH_DIR / "output_achievable_7",
    "autodiscovery": BENCH_DIR / "output_autodiscovery",
}


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
        cell_type = layer_file.stem.replace("_layer_pass1", "")
        df = pd.read_csv(layer_file, index_col=0)

        common_spots = list(adata.obs_names.intersection(df.index))
        if len(common_spots) == 0:
            logger.warning(f"No matching spots for {cell_type}, skipping")
            continue

        layer_matrix = np.zeros((adata.shape[0], adata.shape[1]))

        adata_gene_to_idx = {}
        for idx, gene in enumerate(adata.var_names):
            if gene not in adata_gene_to_idx:
                adata_gene_to_idx[gene] = idx

        common_genes = [g for g in df.columns if g in adata_gene_to_idx]
        if len(common_genes) == 0:
            logger.warning(f"No matching genes for {cell_type}, skipping")
            continue

        spot_idx = np.array([adata.obs_names.get_loc(s) for s in common_spots])
        gene_idx = np.array([adata_gene_to_idx[g] for g in common_genes])

        layer_data = df.loc[common_spots, common_genes].values
        for i, s_idx in enumerate(spot_idx):
            layer_matrix[s_idx, gene_idx] = layer_data[i, :]

        layer_name = f"{cell_type}_genes_pass1"
        adata.layers[layer_name] = layer_matrix
        logger.info(f"  Added layer: {layer_name} ({len(common_spots)} spots, {len(common_genes)} genes)")

    return adata


def load_deconvolved_data(
    source: str,
    region_id: int,
) -> Tuple[sc.AnnData, pd.DataFrame]:
    """
    Load deconvolved data from Module 3 output.

    Args:
        source: "manual" or "autodiscovery"
        region_id: Region ID (0-4)

    Returns:
        Tuple of (adata with layers, cell proportions DataFrame)
    """
    source_dir = DATA_SOURCES[source]

    # Load GEX data for gene names and spatial coordinates
    gex_path = PSEUDOVISIUM_DIR / "data_granular_gt" / "h5ad_objects" / f"Xenium_region_{region_id}_GEX.h5ad"
    adata = sc.read_h5ad(gex_path)
    logger.info(f"Loaded GEX: {adata.shape[0]} spots × {adata.shape[1]} genes")

    # Load cell proportions
    prop_path = source_dir / f"Xenium_region_{region_id}_cell_prop_finetuned_results.csv"
    if not prop_path.exists():
        prop_path = source_dir / f"Xenium_region_{region_id}_cell_prop_global_results.csv"

    proportions = pd.read_csv(prop_path, index_col=0)
    logger.info(f"Loaded proportions: {proportions.shape}")

    # Load deconvolved layers from CSV files
    layers_dir = source_dir / f"Xenium_region_{region_id}_pass1" / "layers" / "pass1"
    if not layers_dir.exists():
        # Try alternative structure
        layers_dir = source_dir / f"Xenium_region_{region_id}" / "layers" / "pass1"

    if layers_dir.exists():
        adata = load_csv_layers_to_adata(adata, layers_dir)
    else:
        logger.warning(f"Layer directory not found: {layers_dir}")
        # Try loading from parquet file
        parquet_path = source_dir / f"Xenium_region_{region_id}_gene_expression_pass1.parquet"
        if parquet_path.exists():
            logger.info(f"Loading from parquet: {parquet_path}")
            # Parquet contains flattened deconvolved expression - need to reshape
            gex_df = pd.read_parquet(parquet_path)
            logger.info(f"Parquet shape: {gex_df.shape}")

    return adata, proportions


def load_singlecell_data(
    region_id: int,
    max_cells: int = 50000,
) -> Tuple[sc.AnnData, sc.AnnData, pd.Series]:
    """
    Load single-cell Xenium data and assign cell types.

    Args:
        region_id: Region ID (0-4)
        max_cells: Maximum cells to load per region

    Returns:
        Tuple of (adata_gex, adata_protein, cell_type_assignments)
    """
    from load_xenium_singlecell import load_xenium_singlecell
    from define_cell_types import classify_cells_by_protein

    # Load single-cell data
    adata_gex, adata_protein = load_xenium_singlecell(
        region_id=region_id,
        max_cells=max_cells,
    )

    # Use ACHIEVABLE_7 cell types for consistency with deconvolved data
    SINGLECELL_PROFILE_DICT = {
        "B cells": {"Major": ["CD20"], "Minor": ["CD45RA"]},
        "CD4+ T cells": {"Major": ["CD3E", "CD4"], "Minor": ["CD45RO"]},
        "CD8+ T cells": {"Major": ["CD3E", "CD8A"], "Minor": ["GranzymeB"]},
        "Macrophages": {"Major": ["CD68", "CD163"], "Minor": ["CD16"]},
        "Endothelial": {"Major": ["CD31"], "Minor": []},
        "Epithelial": {"Major": ["PanCK"], "Minor": []},  # E-Cadherin removed
        "Fibroblasts": {"Major": ["alphaSMA", "Vimentin"], "Minor": []},
    }

    # Assign cell types based on protein expression
    cell_types = classify_cells_by_protein(
        adata_protein,
        SINGLECELL_PROFILE_DICT,
        threshold_method="gmm",
    )

    logger.info(f"Cell type assignments:")
    for ct, count in cell_types.value_counts().items():
        logger.info(f"  {ct}: {count} cells ({100*count/len(cell_types):.1f}%)")

    return adata_gex, adata_protein, cell_types


def create_singlecell_layers(
    adata_gex: sc.AnnData,
    cell_types: pd.Series,
) -> sc.AnnData:
    """
    Create pseudo-deconvolved layers from single-cell data.

    For single-cell data, each cell belongs to one cell type, so we create
    layers by assigning each cell's expression to its cell type layer.

    Args:
        adata_gex: Single-cell GEX data
        cell_types: Cell type assignments

    Returns:
        AnnData with cell-type layers
    """
    adata = adata_gex.copy()

    unique_types = [ct for ct in cell_types.unique() if ct != "Unassigned"]
    logger.info(f"Creating layers for {len(unique_types)} cell types")

    for ct in unique_types:
        mask = cell_types == ct
        layer_matrix = np.zeros(adata.X.shape)

        if hasattr(adata.X, 'toarray'):
            layer_matrix[mask, :] = adata.X[mask, :].toarray()
        else:
            layer_matrix[mask, :] = adata.X[mask, :]

        # Normalize to make comparable across cell types
        row_sums = layer_matrix.sum(axis=1, keepdims=True)
        row_sums[row_sums == 0] = 1
        layer_matrix = layer_matrix / row_sums * 10000  # CPM-like normalization

        layer_name = f"{ct}_genes_pass1"
        adata.layers[layer_name] = layer_matrix
        logger.info(f"  Added layer: {layer_name} ({mask.sum()} cells)")

    return adata


def run_module4(
    adata: sc.AnnData,
    proportions: Optional[pd.DataFrame] = None,
    K_programs: int = 5,
) -> AnchoredProgramDiscoveryResult:
    """
    Run Module 4 program discovery.

    Args:
        adata: AnnData with deconvolved layers
        proportions: Cell type proportions (optional)
        K_programs: Number of programs per cell type

    Returns:
        AnchoredProgramDiscoveryResult
    """
    logger.info("Running Module 4: Anchored Program Discovery")

    result = discover_programs_from_layers(
        adata=adata,
        layer_pattern="_genes_pass1",
        cell_type_proportions=proportions,
        K_programs=K_programs,
        validate_with_proteins=False,  # No protein data for validation
        detect_subpopulations=True,
        n_subpopulations=3,
    )

    logger.info(f"Discovered programs for {len(result.results)} cell types")
    for ct, ct_result in result.results.items():
        logger.info(f"  {ct}: {len(ct_result.programs)} programs")

    return result


def run_module4b(
    module4_result: AnchoredProgramDiscoveryResult,
    adata: sc.AnnData,
) -> BivariateProgramResult:
    """
    Run Module 4b bivariate relationship analysis.

    Args:
        module4_result: Output from Module 4
        adata: AnnData with spatial coordinates

    Returns:
        BivariateProgramResult
    """
    logger.info("Running Module 4b: Bivariate Program Relationships")

    result = analyze_program_relationships(
        module4_result=module4_result,
        adata=adata,
        n_permutations=99,
        fdr_alpha=0.05,
    )

    logger.info(f"Found {len(result.significant_pairs)} significant program pairs")

    return result


def run_module5(
    module4_results: Dict[str, AnchoredProgramDiscoveryResult],
    module4b_results: Optional[Dict[str, BivariateProgramResult]] = None,
) -> IntegrationResult:
    """
    Run Module 5 cross-region integration.

    Args:
        module4_results: Dict of Module 4 results keyed by region name
        module4b_results: Optional dict of Module 4b results

    Returns:
        IntegrationResult
    """
    logger.info(f"Running Module 5: Cross-Region Integration ({len(module4_results)} regions)")

    result = integrate_samples(
        module4_results=module4_results,
        module4b_results=module4b_results,
        n_components=30,
        similarity_threshold=0.7,
    )

    logger.info(f"Found {len(result.aligned_programs)} aligned programs")
    if result.conserved_relationships:
        logger.info(f"Found {len(result.conserved_relationships)} conserved relationships")

    return result


def save_results(
    source: str,
    region_id: int,
    module4_result: AnchoredProgramDiscoveryResult,
    module4b_result: Optional[BivariateProgramResult] = None,
    output_dir: Optional[Path] = None,
):
    """Save Module 4/4b results to disk."""
    if output_dir is None:
        output_dir = OUTPUT_DIR / source / f"region_{region_id}"

    output_dir.mkdir(parents=True, exist_ok=True)

    # Save Module 4 results as JSON
    module4_summary = {
        "source": source,
        "region_id": region_id,
        "n_cell_types": len(module4_result.results),
        "cell_types": list(module4_result.results.keys()),
        "programs_per_type": {
            ct: len(ct_result.programs)
            for ct, ct_result in module4_result.results.items()
        },
    }

    with open(output_dir / "module4_summary.json", "w") as f:
        json.dump(module4_summary, f, indent=2)

    # Save top genes per program
    top_genes_data = []
    for ct, ct_result in module4_result.results.items():
        for prog in ct_result.programs:
            top_genes_data.append({
                "cell_type": ct,
                "program_id": prog.program_id,
                "top_genes": prog.top_genes[:20],
                "variance_explained": prog.variance_explained,
                "spatial_moran_i": prog.spatial_moran_i,
            })

    pd.DataFrame(top_genes_data).to_csv(output_dir / "module4_programs.csv", index=False)

    # Save Module 4b results if available
    if module4b_result is not None:
        pairs_data = []
        for pair in module4b_result.all_pairs:
            pairs_data.append({
                "program_1": f"{pair.cell_type_1}_{pair.program_id_1}",
                "program_2": f"{pair.cell_type_2}_{pair.program_id_2}",
                "bivariate_morans_i": pair.bivariate_morans_i,
                "p_value": pair.p_value,
                "relationship_type": pair.relationship_type,
            })

        pd.DataFrame(pairs_data).to_csv(output_dir / "module4b_pairs.csv", index=False)

    logger.info(f"Results saved to {output_dir}")


def run_validation_for_source(
    source: str,
    region_ids: List[int],
    K_programs: int = 5,
    max_cells_singlecell: int = 50000,
) -> Dict[str, Dict]:
    """
    Run full Module 4/4b/5 validation for a data source.

    Args:
        source: "singlecell", "manual", or "autodiscovery"
        region_ids: List of region IDs to process
        K_programs: Number of programs per cell type
        max_cells_singlecell: Max cells for single-cell data

    Returns:
        Dict with Module 4 and 4b results per region, plus Module 5 integration
    """
    results = {
        "module4": {},
        "module4b": {},
        "module5": None,
    }

    for region_id in region_ids:
        logger.info(f"\n{'='*60}")
        logger.info(f"Processing {source} - Region {region_id}")
        logger.info(f"{'='*60}")

        try:
            if source == "singlecell":
                adata_gex, adata_protein, cell_types = load_singlecell_data(
                    region_id=region_id,
                    max_cells=max_cells_singlecell,
                )
                adata = create_singlecell_layers(adata_gex, cell_types)
                proportions = None
            else:
                adata, proportions = load_deconvolved_data(source, region_id)

            # Check if we have layers
            deconv_layers = [l for l in adata.layers.keys() if "_genes_pass1" in l]
            if not deconv_layers:
                logger.warning(f"No deconvolved layers found for {source} region {region_id}")
                continue

            # Run Module 4
            module4_result = run_module4(adata, proportions, K_programs)
            results["module4"][f"region_{region_id}"] = module4_result

            # Run Module 4b
            module4b_result = run_module4b(module4_result, adata)
            results["module4b"][f"region_{region_id}"] = module4b_result

            # Save results
            save_results(source, region_id, module4_result, module4b_result)

        except Exception as e:
            logger.error(f"Error processing {source} region {region_id}: {e}")
            import traceback
            traceback.print_exc()
            continue

    # Run Module 5 if we have multiple regions
    if len(results["module4"]) > 1:
        logger.info(f"\n{'='*60}")
        logger.info(f"Running Module 5 Integration for {source}")
        logger.info(f"{'='*60}")

        try:
            module5_result = run_module5(
                results["module4"],
                results["module4b"],
            )
            results["module5"] = module5_result

            # Save Module 5 results
            output_dir = OUTPUT_DIR / source
            output_dir.mkdir(parents=True, exist_ok=True)

            integration_summary = {
                "source": source,
                "n_regions": len(results["module4"]),
                "n_aligned_programs": len(module5_result.aligned_programs),
                "n_conserved_relationships": len(module5_result.conserved_relationships) if module5_result.conserved_relationships else 0,
            }

            with open(output_dir / "module5_integration.json", "w") as f:
                json.dump(integration_summary, f, indent=2)

        except Exception as e:
            logger.error(f"Error in Module 5 integration for {source}: {e}")
            import traceback
            traceback.print_exc()

    return results


def main():
    parser = argparse.ArgumentParser(description="Module 4/5 Validation on Xenium Data")
    parser.add_argument(
        "--source",
        type=str,
        choices=["singlecell", "manual", "autodiscovery", "all"],
        default="all",
        help="Data source to validate",
    )
    parser.add_argument(
        "--region-id",
        type=int,
        default=None,
        help="Single region to process (0-4)",
    )
    parser.add_argument(
        "--all-regions",
        action="store_true",
        help="Process all 5 regions",
    )
    parser.add_argument(
        "--K-programs",
        type=int,
        default=5,
        help="Number of programs per cell type",
    )
    parser.add_argument(
        "--max-cells",
        type=int,
        default=50000,
        help="Max cells for single-cell data",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default=None,
        help="Output directory",
    )

    args = parser.parse_args()

    # Determine regions to process
    if args.all_regions:
        region_ids = list(range(5))
    elif args.region_id is not None:
        region_ids = [args.region_id]
    else:
        region_ids = [0]  # Default to region 0

    # Override output directory if specified
    global OUTPUT_DIR
    if args.output_dir:
        OUTPUT_DIR = Path(args.output_dir)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Determine sources to process
    if args.source == "all":
        sources = ["singlecell", "manual", "autodiscovery"]
    else:
        sources = [args.source]

    # Run validation
    all_results = {}
    for source in sources:
        logger.info(f"\n{'#'*70}")
        logger.info(f"# VALIDATING SOURCE: {source.upper()}")
        logger.info(f"{'#'*70}")

        results = run_validation_for_source(
            source=source,
            region_ids=region_ids,
            K_programs=args.K_programs,
            max_cells_singlecell=args.max_cells,
        )
        all_results[source] = results

    # Save overall summary
    summary = {
        "sources": sources,
        "regions": region_ids,
        "K_programs": args.K_programs,
        "results": {
            source: {
                "n_regions_processed": len(results["module4"]),
                "module5_completed": results["module5"] is not None,
            }
            for source, results in all_results.items()
        },
    }

    with open(OUTPUT_DIR / "validation_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    logger.info(f"\n{'='*60}")
    logger.info("VALIDATION COMPLETE")
    logger.info(f"{'='*60}")
    logger.info(f"Results saved to: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
