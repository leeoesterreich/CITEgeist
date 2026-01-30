"""
Run Tangram deconvolution on Xenium pseudo-Visium data.

This script runs Tangram spatial mapping using the GSE156632 reference
to predict cell type proportions and per-cell-type gene expression.

Tangram maps single cells to spatial locations by learning a probabilistic
cell-to-spot assignment that maximizes gene expression correlation.

Usage:
    python run_benchmark.py --region-id 0 --ref-path /path/to/reference.h5ad --output-dir /path/to/output
"""

import os
import sys
import argparse
import logging
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

import tangram as tg

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)

# Target cell types (must match protein-gated ground truth)
TARGET_CELL_TYPES = [
    "B cells",
    "CD4+ T cells",
    "CD8+ T cells",
    "Macrophages",
    "Endothelial",
    "Epithelial",
    "Fibroblasts",
]


def find_marker_genes(adata_ref, labels_key: str = "cell_type", n_markers: int = 100):
    """
    Find marker genes for each cell type.

    Args:
        adata_ref: Reference AnnData with cell type labels
        labels_key: Column for cell type labels
        n_markers: Number of top markers per cell type

    Returns:
        List of marker gene names
    """
    logger.info(f"Finding {n_markers} marker genes per cell type...")

    # Copy to avoid modifying original
    adata = adata_ref.copy()

    # Rank genes
    sc.tl.rank_genes_groups(adata, groupby=labels_key, method='wilcoxon')

    # Extract top markers for each cell type
    markers = set()
    cell_types = adata.obs[labels_key].unique()

    for ct in cell_types:
        ct_markers = sc.get.rank_genes_groups_df(adata, group=ct)
        top_markers = ct_markers.head(n_markers)['names'].tolist()
        markers.update(top_markers)
        logger.info(f"  {ct}: {len(top_markers)} markers")

    marker_list = list(markers)
    logger.info(f"Total unique markers: {len(marker_list)}")

    return marker_list


def preprocess_for_tangram(adata_ref, adata_sp, markers, labels_key: str = "cell_type"):
    """
    Preprocess reference and spatial data for Tangram.

    Args:
        adata_ref: Reference AnnData
        adata_sp: Spatial AnnData
        markers: List of marker genes
        labels_key: Column for cell type labels

    Returns:
        Preprocessed reference and spatial AnnData
    """
    logger.info("Preprocessing data for Tangram...")

    # Filter to shared genes
    shared_genes = np.intersect1d(
        np.intersect1d(adata_ref.var_names, adata_sp.var_names),
        markers
    )
    logger.info(f"  Shared marker genes: {len(shared_genes)}")

    # Preprocess with Tangram
    tg.pp_adatas(adata_ref, adata_sp, genes=list(shared_genes))

    logger.info(f"  Reference after preprocessing: {adata_ref.shape}")
    logger.info(f"  Spatial after preprocessing: {adata_sp.shape}")

    return adata_ref, adata_sp


def run_tangram_mapping(
    adata_ref,
    adata_sp,
    labels_key: str = "cell_type",
    num_epochs: int = 500,
    device: str = "cuda:0",
):
    """
    Run Tangram cell-to-space mapping.

    Args:
        adata_ref: Preprocessed reference AnnData
        adata_sp: Preprocessed spatial AnnData
        labels_key: Column for cell type labels
        num_epochs: Training epochs
        device: Device for training

    Returns:
        Mapping AnnData
    """
    logger.info(f"Running Tangram mapping ({num_epochs} epochs)...")
    start_time = time.time()

    # Run mapping
    ad_map = tg.map_cells_to_space(
        adata_sc=adata_ref,
        adata_sp=adata_sp,
        device=device,
        num_epochs=num_epochs,
    )

    train_time = time.time() - start_time
    logger.info(f"Mapping completed in {train_time:.1f}s ({train_time/60:.1f} min)")

    return ad_map, train_time


def extract_proportions(ad_map, adata_ref, adata_sp, labels_key: str = "cell_type"):
    """
    Extract cell type proportions from Tangram mapping.

    Args:
        ad_map: Mapping AnnData from Tangram
        adata_ref: Reference AnnData
        adata_sp: Spatial AnnData
        labels_key: Column for cell type labels

    Returns:
        Proportions DataFrame
    """
    logger.info("Extracting cell type proportions...")

    # Project cell annotations to space
    tg.project_cell_annotations(ad_map, adata_sp, annotation=labels_key)

    # Get proportions (stored in adata_sp.obsm['tangram_ct_pred'])
    if 'tangram_ct_pred' in adata_sp.obsm:
        proportions = pd.DataFrame(
            adata_sp.obsm['tangram_ct_pred'],
            index=adata_sp.obs_names,
            columns=adata_sp.uns['tangram_ct_pred_names'] if 'tangram_ct_pred_names' in adata_sp.uns else None
        )
    else:
        # Alternative: compute from mapping matrix
        logger.info("  Computing proportions from mapping matrix...")
        cell_types = adata_ref.obs[labels_key]
        mapping_matrix = ad_map.X  # (n_cells, n_spots)

        # Get proportion of each cell type per spot
        proportions = pd.DataFrame(
            index=adata_sp.obs_names,
            columns=cell_types.unique()
        )

        for ct in cell_types.unique():
            ct_mask = (cell_types == ct).values
            ct_weights = mapping_matrix[ct_mask, :].sum(axis=0)
            proportions[ct] = ct_weights

        # Normalize
        proportions = proportions.div(proportions.sum(axis=1), axis=0)

    logger.info(f"  Proportions shape: {proportions.shape}")
    return proportions


def extract_gex_layers(ad_map, adata_ref, adata_sp, labels_key: str = "cell_type"):
    """
    Extract per-cell-type gene expression layers.

    Tangram projects gene expression by weighting reference cells
    by their mapping probability to each spot.

    Args:
        ad_map: Mapping AnnData from Tangram
        adata_ref: Reference AnnData
        adata_sp: Spatial AnnData
        labels_key: Column for cell type labels

    Returns:
        Dictionary of layer DataFrames
    """
    logger.info("Extracting per-cell-type gene expression layers...")

    # Get mapping matrix: (n_cells, n_spots)
    mapping_matrix = ad_map.X

    # Get cell types
    cell_types = adata_ref.obs[labels_key].values
    unique_cts = np.unique(cell_types)

    # Get gene expression matrix
    if hasattr(adata_ref.X, 'toarray'):
        ref_expr = adata_ref.X.toarray()
    else:
        ref_expr = adata_ref.X

    layer_dfs = {}

    for ct in unique_cts:
        logger.info(f"  Processing {ct}...")

        # Get mask for this cell type
        ct_mask = (cell_types == ct)
        n_cells_ct = ct_mask.sum()

        # Get mapping weights for this cell type
        ct_mapping = mapping_matrix[ct_mask, :]  # (n_cells_ct, n_spots)

        # Normalize mapping weights per spot
        ct_weights = ct_mapping / (ct_mapping.sum(axis=0, keepdims=True) + 1e-10)

        # Weighted average of gene expression
        ct_expr = ref_expr[ct_mask, :]  # (n_cells_ct, n_genes)
        layer_expr = ct_weights.T @ ct_expr  # (n_spots, n_genes)

        # Create DataFrame
        layer_df = pd.DataFrame(
            layer_expr,
            index=adata_sp.obs_names,
            columns=adata_ref.var_names,
        )
        layer_dfs[ct] = layer_df

    return layer_dfs


def main():
    parser = argparse.ArgumentParser(
        description="Run Tangram on Xenium pseudo-Visium data"
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
        default="../../xenium_pseudovisium/data_rna_gt",
        help="Input directory containing h5ad_objects/",
    )
    parser.add_argument(
        "--ref-path",
        type=str,
        required=True,
        help="Path to processed reference h5ad file",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="output_rna_gt",
        help="Output directory for results",
    )
    parser.add_argument(
        "--n-markers",
        type=int,
        default=100,
        help="Number of marker genes per cell type (default: 100)",
    )
    parser.add_argument(
        "--num-epochs",
        type=int,
        default=500,
        help="Training epochs (default: 500)",
    )
    parser.add_argument(
        "--no-gpu",
        action="store_true",
        help="Disable GPU usage",
    )
    parser.add_argument(
        "--prefix",
        type=str,
        default="Xenium",
        help="Filename prefix (default: Xenium)",
    )

    args = parser.parse_args()

    # Resolve paths
    script_dir = Path(__file__).parent.parent
    input_dir = Path(args.input_dir)
    if not input_dir.is_absolute():
        input_dir = script_dir / input_dir
    ref_path = Path(args.ref_path)
    if not ref_path.is_absolute():
        ref_path = script_dir / ref_path
    output_dir = Path(args.output_dir)
    if not output_dir.is_absolute():
        output_dir = script_dir / output_dir

    sample_name = f"{args.prefix}_region_{args.region_id}"
    device = "cpu" if args.no_gpu else "cuda:0"

    logger.info("=" * 60)
    logger.info(f"Tangram Deconvolution - Region {args.region_id}")
    logger.info("=" * 60)
    logger.info(f"Input: {input_dir}")
    logger.info(f"Reference: {ref_path}")
    logger.info(f"Output: {output_dir}")
    logger.info(f"Device: {device}")

    # Load reference
    logger.info("\nLoading reference data...")
    adata_ref = sc.read_h5ad(ref_path)
    logger.info(f"  Reference shape: {adata_ref.shape}")
    logger.info(f"  Cell types: {adata_ref.obs['cell_type'].value_counts().to_dict()}")

    # Load spatial data
    logger.info(f"\nLoading spatial data for region {args.region_id}...")
    gex_path = input_dir / "h5ad_objects" / f"{args.prefix}_region_{args.region_id}_GEX.h5ad"
    adata_sp = sc.read_h5ad(gex_path)
    logger.info(f"  Spatial shape: {adata_sp.shape}")

    # Find marker genes
    markers = find_marker_genes(adata_ref, n_markers=args.n_markers)

    # Preprocess
    adata_ref, adata_sp = preprocess_for_tangram(adata_ref, adata_sp, markers)

    # Run Tangram mapping
    ad_map, train_time = run_tangram_mapping(
        adata_ref, adata_sp,
        num_epochs=args.num_epochs,
        device=device,
    )

    # Create output directory
    region_output_dir = output_dir / sample_name
    region_output_dir.mkdir(parents=True, exist_ok=True)

    # Extract proportions
    proportions = extract_proportions(ad_map, adata_ref, adata_sp)

    # Save with standard naming convention for benchmarking
    prop_path = region_output_dir / f"{sample_name}_cell_prop_finetuned_results.csv"
    proportions.to_csv(prop_path)
    logger.info(f"Saved proportions to {prop_path}")

    # Also save with legacy naming for compatibility
    legacy_path = region_output_dir / f"{sample_name}_deconv_predictions.csv"
    proportions.to_csv(legacy_path)
    logger.info(f"Saved legacy proportions to {legacy_path}")

    # Extract GEX layers
    layer_dfs = extract_gex_layers(ad_map, adata_ref, adata_sp)
    layers_dir = region_output_dir / "layers"
    layers_dir.mkdir(parents=True, exist_ok=True)

    for ct, layer_df in layer_dfs.items():
        layer_path = layers_dir / f"{ct}_layer.csv"
        layer_df.to_csv(layer_path)
        logger.info(f"  Saved {ct}: {layer_df.shape}")

    # Save results summary
    stats = {
        "region_id": args.region_id,
        "sample_name": sample_name,
        "n_spots": adata_sp.n_obs,
        "n_genes": adata_sp.n_vars,
        "n_markers": len(markers),
        "n_cell_types": len(proportions.columns),
        "cell_types": proportions.columns.tolist(),
        "train_time_sec": train_time,
        "num_epochs": args.num_epochs,
    }

    with open(region_output_dir / "run_stats.json", "w") as f:
        json.dump(stats, f, indent=2)

    # Save mapping
    ad_map.write(region_output_dir / f"{sample_name}_mapping.h5ad")

    logger.info("\n" + "=" * 60)
    logger.info("Deconvolution Complete!")
    logger.info("=" * 60)
    logger.info(f"  Spots: {adata_sp.n_obs}")
    logger.info(f"  Cell types: {len(proportions.columns)}")
    logger.info(f"  Train time: {train_time:.1f}s ({train_time/60:.1f} min)")
    logger.info(f"  Proportions: {prop_path}")
    logger.info(f"  GEX layers: {layers_dir}")


if __name__ == "__main__":
    main()
