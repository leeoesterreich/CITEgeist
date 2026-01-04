"""
Run CITEgeist deconvolution on Xenium pseudo-Visium data.

This module provides wrapper functions to run CITEgeist on the
pseudo-Visium datasets generated from Xenium single-cell data.
"""

import os
import sys
import argparse
import logging
import time
import json
from pathlib import Path
from typing import Dict, Any, Optional

import numpy as np
import pandas as pd
import scanpy as sc

# Add CITEgeist to path
# Path: Benchmarking/xenium_benchmarking/CITEgeist/src/run_benchmark.py -> 5 levels to repo root
REPO_ROOT = Path(__file__).parent.parent.parent.parent.parent
sys.path.insert(0, str(REPO_ROOT / "CITEgeist"))

from model.citegeist_model import CitegeistModel

logger = logging.getLogger(__name__)

# Cell profile dictionary for Xenium proteins (matching define_cell_types.py)
# Uses the same markers for deconvolution as for ground truth generation
XENIUM_CELL_PROFILE_DICT = {
    "CD4+ T cells": {
        "Major": ["CD3E", "CD4"],
        "Minor": ["CD45"],
    },
    "CD8+ T cells": {
        "Major": ["CD3E", "CD8A"],
        "Minor": ["CD45", "GranzymeB"],
    },
    "B cells": {
        "Major": ["CD20"],
        "Minor": ["CD45"],
    },
    "Plasma cells": {
        "Major": ["CD138"],
        "Minor": [],
    },
    "Macrophages": {
        "Major": ["CD68"],
        "Minor": ["CD163", "CD45", "HLA-DR"],
    },
    "Dendritic cells": {
        "Major": ["CD11c", "HLA-DR"],
        "Minor": ["CD45"],
    },
    "NK cells": {
        "Major": ["CD16", "GranzymeB"],
        "Minor": ["CD45"],
    },
    "Epithelial": {
        "Major": ["PanCK"],
        "Minor": ["E-Cadherin"],
    },
    "Endothelial": {
        "Major": ["CD31"],
        "Minor": [],
    },
    "Fibroblasts": {
        "Major": ["alphaSMA"],
        "Minor": ["Vimentin"],
    },
}

# RNA-based cell profile dictionary (matching RNA clustering ground truth)
# These 6 cell types match the simplified RNA-cluster-based ground truth.
# Using RNA-based clustering for ground truth avoids circular logic.
#
# Reference:
#   Zhao et al. (2025). "Benchmarking cell type annotation methods for 10x
#   Xenium spatial transcriptomics data." BMC Bioinformatics, 26(1), 25.
#   https://doi.org/10.1186/s12859-025-06044-0
RNA_CELL_PROFILE_DICT = {
    "B cells": {
        "Major": ["CD20"],
        "Minor": [],
    },
    "T cells": {
        "Major": ["CD3E", "CD8A"],
        "Minor": ["CD4"],
    },
    "Macrophages": {
        "Major": ["CD68", "HLA-DR"],
        "Minor": ["CD163"],
    },
    "Fibroblasts": {
        "Major": ["alphaSMA", "Vimentin"],
        "Minor": [],
    },
    "Epithelial": {
        "Major": ["PanCK", "E-Cadherin"],
        "Minor": [],
    },
    "Endothelial": {
        "Major": ["CD31"],
        "Minor": [],
    },
}

# Granular 10-cell type profile dictionary (unsimplified RNA clustering)
# These 10 cell types match the unsimplified RNA k-means clusters.
# This provides maximum granularity to highlight CITEgeist's proteomic advantage.
#
# Cluster analysis performed 2026-01-02 (see analyze_cluster_profiles.py):
# Cluster 1: CD8+ T cells     - CD3E=378, CD8A=210
# Cluster 2: Macrophages      - CD68=430, CD163=88
# Cluster 3: Mixed Immune     - CD3E=142, CD8A=118, HLA-DR=142
# Cluster 4: Epithelial       - PanCK=39, Vimentin=311
# Cluster 5: Myofibroblasts   - alphaSMA=108, Vimentin=374
# Cluster 6: Stromal          - Mixed low markers
# Cluster 7: Endothelial      - CD31=168
# Cluster 8: B cells          - CD20=293, CD45RA=398
# Cluster 9: Proliferating T  - CD3E=679, PCNA=83
# Cluster 10: Vascular Stromal - CD31=53, Vimentin=209
GRANULAR_CELL_PROFILE_DICT = {
    "CD8+ T cells": {
        "Major": ["CD3E", "CD8A"],
        "Minor": ["CD45", "GranzymeB"],
    },
    "Macrophages": {
        "Major": ["CD68", "HLA-DR"],
        "Minor": ["CD163", "CD16"],
    },
    "Mixed Immune": {
        "Major": ["CD3E", "HLA-DR"],
        "Minor": ["CD8A", "CD45"],
    },
    "Epithelial": {
        "Major": ["PanCK"],
        "Minor": ["E-Cadherin", "Vimentin"],
    },
    "Myofibroblasts": {
        "Major": ["alphaSMA", "Vimentin"],
        "Minor": [],
    },
    "Stromal": {
        "Major": ["Vimentin"],
        "Minor": ["CD3E"],  # Low but detectable
    },
    "Endothelial": {
        "Major": ["CD31"],
        "Minor": ["Vimentin"],
    },
    "B cells": {
        "Major": ["CD20"],
        "Minor": ["CD45", "CD45RA"],
    },
    "Proliferating T": {
        "Major": ["CD3E"],
        "Minor": ["PCNA", "Ki-67"],
    },
    "Vascular Stromal": {
        "Major": ["Vimentin"],
        "Minor": ["CD31"],
    },
}


def run_citegeist_on_region(
    region_id: int,
    input_dir: str,
    output_dir: str,
    cell_profile_dict: Optional[Dict] = None,
    radius: float = 4.0,
    lambda_reg: float = 1.0,
    alpha_elastic: float = 0.7,
    max_y_change: float = 0.4,
    lambda_laplacian: float = 0.1,
    laplacian_k: int = 8,
    min_counts: int = 25,
    nonzero_percentage: float = 0.01,
    mean_expression_threshold: float = 1.1,
    prefix: str = "Xenium",
    run_gex: bool = False,
) -> Dict[str, Any]:
    """
    Run CITEgeist deconvolution on a single region.

    Args:
        region_id: Region identifier (0-4)
        input_dir: Directory containing h5ad_objects/
        output_dir: Directory to save outputs
        cell_profile_dict: Cell type to marker mapping (default: XENIUM_CELL_PROFILE_DICT)
        radius: Neighbor detection radius for spatial prior
        lambda_reg: Regularization strength
        alpha_elastic: Elastic net mixing parameter
        max_y_change: Maximum change in Y during optimization
        lambda_laplacian: Laplacian smoothing strength
        laplacian_k: Number of neighbors for Laplacian
        min_counts: Minimum counts per spot
        nonzero_percentage: Min percentage of spots with expression
        mean_expression_threshold: Min mean expression in nonzero spots
        prefix: Filename prefix
        run_gex: Whether to run gene expression deconvolution (Pass 1)

    Returns:
        Dict with run statistics and output paths
    """
    if cell_profile_dict is None:
        cell_profile_dict = XENIUM_CELL_PROFILE_DICT

    input_dir = Path(input_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # File paths
    gex_path = input_dir / "h5ad_objects" / f"{prefix}_region_{region_id}_GEX.h5ad"
    protein_path = input_dir / "h5ad_objects" / f"{prefix}_region_{region_id}_CITE.h5ad"

    logger.info(f"Loading region {region_id} data...")
    logger.info(f"  GEX: {gex_path}")
    logger.info(f"  Protein: {protein_path}")

    # Load data
    adata_gex = sc.read_h5ad(gex_path)
    adata_protein = sc.read_h5ad(protein_path)

    logger.info(f"  GEX shape: {adata_gex.shape}")
    logger.info(f"  Protein shape: {adata_protein.shape}")

    # Initialize model in simulation mode
    sample_name = f"{prefix}_region_{region_id}"
    model = CitegeistModel(
        sample_name=sample_name,
        output_folder=str(output_dir),
        simulation=True,
        gene_expression_adata=adata_gex,
        antibody_capture_adata=adata_protein,
    )

    # Load cell profile dictionary
    model.load_cell_profile_dict(cell_profile_dict)

    # Preprocessing
    logger.info("Preprocessing gene expression...")
    model.filter_gex(
        nonzero_percentage=nonzero_percentage,
        mean_expression_threshold=mean_expression_threshold,
        min_counts=min_counts,
    )
    model.preprocess_gex(target_sum=10000)

    logger.info("Preprocessing antibody/protein...")
    model.preprocess_antibody()

    # Run cell proportion model
    logger.info("Running cell proportion optimization...")
    start_time = time.time()

    global_props, finetuned_props = model.run_cell_proportion_model(
        radius=radius,
        lambda_reg=lambda_reg,
        alpha=alpha_elastic,
        max_y_change=max_y_change,
        lambda_laplacian=lambda_laplacian,
        laplacian_k=laplacian_k,
        per_marker_beta=False,  # Disable per-marker beta (made results worse)
        validation_warn_only=True,  # Don't fail on high Unknown proportion
    )

    prop_time = time.time() - start_time
    logger.info(f"Cell proportion optimization completed in {prop_time:.1f}s")

    # Run gene expression deconvolution if requested
    gex_time = 0.0
    if run_gex:
        logger.info("Running gene expression deconvolution (Pass 1)...")
        gex_start = time.time()

        try:
            pass1_results = model.run_cell_expression_pass1(
                radius=radius,
                alpha=0.5,
                lambda_reg_gex=0.001,
                global_enrichment_weight=0.5,
                local_enrichment_weight=0.5,
                checkpoint_interval=100,
                output_dir=str(output_dir / "checkpoints"),
                rerun=True,
            )
            gex_time = time.time() - gex_start
            logger.info(f"Gene expression deconvolution completed in {gex_time:.1f}s")
        except Exception as e:
            logger.error(f"Gene expression deconvolution failed: {e}")
            gex_time = -1.0  # Mark as failed

    # Get output paths - results are saved by the model during run
    result_dir = output_dir / sample_name
    result_dir.mkdir(parents=True, exist_ok=True)
    prop_csv = result_dir / "cell_proportions.csv"

    # Save proportions in benchmarking-compatible format
    if finetuned_props is not None:
        props_df = pd.DataFrame(
            finetuned_props,
            index=model.antibody_capture_adata.obs_names,
            columns=list(cell_profile_dict.keys()) + ["Unknown"],
        )
        props_df.to_csv(result_dir / f"{sample_name}_deconv_predictions.csv")

    # Run statistics
    stats = {
        "region_id": region_id,
        "sample_name": sample_name,
        "n_spots": adata_gex.shape[0],
        "n_genes_filtered": model.gene_expression_adata.shape[1],
        "n_proteins": adata_protein.shape[1],
        "n_cell_types": len(cell_profile_dict) + 1,  # +1 for Unknown
        "runtime_proportions_sec": prop_time,
        "runtime_gex_sec": gex_time if run_gex else None,
        "output_dir": str(result_dir),
        "prop_csv": str(prop_csv),
        "gex_layers_dir": str(output_dir / f"{sample_name}_pass1" / "layers") if run_gex else None,
    }

    # Save stats
    with open(result_dir / "run_stats.json", "w") as f:
        json.dump(stats, f, indent=2)

    logger.info(f"Results saved to {result_dir}")

    return stats


def run_all_regions(
    input_dir: str,
    output_dir: str,
    n_regions: int = 5,
    **kwargs,
) -> Dict[str, Any]:
    """
    Run CITEgeist on all regions.

    Args:
        input_dir: Directory containing h5ad_objects/
        output_dir: Directory to save outputs
        n_regions: Number of regions
        **kwargs: Additional arguments for run_citegeist_on_region

    Returns:
        Dict with summary statistics
    """
    all_stats = []

    for region_id in range(n_regions):
        logger.info(f"\n{'='*60}")
        logger.info(f"Processing region {region_id}/{n_regions-1}")
        logger.info(f"{'='*60}")

        stats = run_citegeist_on_region(
            region_id=region_id,
            input_dir=input_dir,
            output_dir=output_dir,
            **kwargs,
        )
        all_stats.append(stats)

    # Summary
    total_runtime = sum(s["runtime_proportions_sec"] for s in all_stats)
    total_spots = sum(s["n_spots"] for s in all_stats)

    summary = {
        "n_regions": n_regions,
        "total_spots": total_spots,
        "total_runtime_sec": total_runtime,
        "mean_runtime_per_region": total_runtime / n_regions,
        "regions": all_stats,
    }

    # Save summary
    output_dir = Path(output_dir)
    with open(output_dir / "benchmark_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    logger.info(f"\n{'='*60}")
    logger.info("All regions completed!")
    logger.info(f"Total runtime: {total_runtime:.1f}s ({total_runtime/60:.1f} min)")
    logger.info(f"Total spots: {total_spots}")
    logger.info(f"{'='*60}")

    return summary


def main():
    parser = argparse.ArgumentParser(
        description="Run CITEgeist on Xenium pseudo-Visium data"
    )
    parser.add_argument(
        "--region-id",
        type=int,
        required=True,
        help="Region ID to process (0-4)",
    )
    parser.add_argument(
        "--input-dir",
        type=str,
        default="Benchmarking/xenium_pseudovisium/data",
        help="Input directory containing h5ad_objects/",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="Benchmarking/xenium_benchmarking/CITEgeist/output",
        help="Output directory for results",
    )
    parser.add_argument(
        "--radius",
        type=float,
        default=4.0,
        help="Neighbor detection radius",
    )
    parser.add_argument(
        "--lambda-reg",
        type=float,
        default=1.0,
        help="Regularization strength",
    )
    parser.add_argument(
        "--alpha-elastic",
        type=float,
        default=0.7,
        help="Elastic net mixing parameter",
    )
    parser.add_argument(
        "--max-y-change",
        type=float,
        default=0.4,
        help="Maximum change in Y during optimization",
    )
    parser.add_argument(
        "--min-counts",
        type=int,
        default=25,
        help="Minimum counts per spot",
    )
    parser.add_argument(
        "--use-rna-profiles",
        action="store_true",
        help="Use RNA-based cell profile dict (6 cell types matching simplified RNA clustering GT)",
    )
    parser.add_argument(
        "--use-granular-profiles",
        action="store_true",
        help="Use granular cell profile dict (10 cell types matching unsimplified RNA clustering GT)",
    )
    parser.add_argument(
        "--run-gex",
        action="store_true",
        help="Run gene expression deconvolution (Pass 1) after cell proportions",
    )

    args = parser.parse_args()

    # Setup logging
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    )

    # Select cell profile dict based on flags
    if args.use_granular_profiles:
        cell_profile_dict = GRANULAR_CELL_PROFILE_DICT
        logger.info("Using GRANULAR cell profile dict (10 cell types, unsimplified RNA clustering)")
    elif args.use_rna_profiles:
        cell_profile_dict = RNA_CELL_PROFILE_DICT
        logger.info("Using RNA-based cell profile dict (6 cell types, simplified)")
    else:
        cell_profile_dict = XENIUM_CELL_PROFILE_DICT
        logger.info("Using protein-based cell profile dict (10 cell types)")

    # Run
    stats = run_citegeist_on_region(
        region_id=args.region_id,
        input_dir=args.input_dir,
        output_dir=args.output_dir,
        cell_profile_dict=cell_profile_dict,
        radius=args.radius,
        lambda_reg=args.lambda_reg,
        alpha_elastic=args.alpha_elastic,
        max_y_change=args.max_y_change,
        min_counts=args.min_counts,
        run_gex=args.run_gex,
    )

    print(f"\nCompleted region {args.region_id}")
    print(f"  Spots: {stats['n_spots']}")
    print(f"  Runtime (proportions): {stats['runtime_proportions_sec']:.1f}s")
    if args.run_gex and stats.get('runtime_gex_sec'):
        print(f"  Runtime (GEX): {stats['runtime_gex_sec']:.1f}s")
        print(f"  GEX layers: {stats['gex_layers_dir']}")
    print(f"  Output: {stats['output_dir']}")


if __name__ == "__main__":
    main()
