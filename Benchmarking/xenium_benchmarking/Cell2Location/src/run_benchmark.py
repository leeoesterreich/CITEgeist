"""
Run Cell2Location deconvolution on Xenium pseudo-Visium data.

This script runs Cell2Location spatial deconvolution using the trained
reference model from GSE156632 to predict cell type proportions and
per-cell-type gene expression layers.

Usage:
    python run_benchmark.py --region-id 0 --ref-dir /path/to/reference --output-dir /path/to/output
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

import cell2location
from cell2location.models import Cell2location, RegressionModel

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)

# Target cell types (must match Xenium RNA-based ground truth)
TARGET_CELL_TYPES = [
    "B cells",
    "T cells",
    "Macrophages",
    "Fibroblasts",
    "Epithelial",
    "Endothelial",
]


def load_reference_signatures(ref_dir: Path):
    """
    Load trained reference model and signatures.

    Args:
        ref_dir: Directory containing trained reference model

    Returns:
        Reference AnnData and signature DataFrame
    """
    logger.info(f"Loading reference from {ref_dir}...")

    # Load reference adata
    ref_path = ref_dir / "reference_trained.h5ad"
    adata_ref = sc.read_h5ad(ref_path)
    logger.info(f"  Reference shape: {adata_ref.shape}")

    # Load signatures
    sig_path = ref_dir / "cell_type_signatures.csv"
    if sig_path.exists():
        signatures = pd.read_csv(sig_path, index_col=0)
        logger.info(f"  Signatures shape: {signatures.shape}")
        logger.info(f"  Cell types: {signatures.columns.tolist()}")
    else:
        # Extract from adata
        factor_names = adata_ref.uns['mod']['factor_names']
        if 'means_per_cluster_mu_fg' in adata_ref.varm.keys():
            signature_cols = [f'means_per_cluster_mu_fg_{ct}' for ct in factor_names]
            signatures = adata_ref.varm['means_per_cluster_mu_fg'][signature_cols].copy()
        else:
            signature_cols = [f'means_per_cluster_mu_fg_{ct}' for ct in factor_names]
            signatures = adata_ref.var[signature_cols].copy()
        signatures.columns = factor_names
        logger.info(f"  Extracted signatures: {signatures.shape}")

    return adata_ref, signatures


def preprocess_spatial_data(adata_vis, signatures):
    """
    Preprocess spatial data for Cell2Location.

    Args:
        adata_vis: Spatial AnnData (gene expression)
        signatures: Cell type signature DataFrame

    Returns:
        Preprocessed AnnData and filtered signatures
    """
    logger.info("Preprocessing spatial data...")

    # Ensure gene symbols are set
    if 'SYMBOL' not in adata_vis.var.columns:
        adata_vis.var['SYMBOL'] = adata_vis.var_names

    # Handle mitochondrial genes
    adata_vis.var['MT_gene'] = adata_vis.var['SYMBOL'].str.startswith('MT-')
    n_mt = adata_vis.var['MT_gene'].sum()
    logger.info(f"  Mitochondrial genes: {n_mt}")

    if n_mt > 0:
        # Store MT counts separately
        adata_vis.obsm['MT'] = adata_vis[:, adata_vis.var['MT_gene'].values].X
        if hasattr(adata_vis.obsm['MT'], 'toarray'):
            adata_vis.obsm['MT'] = adata_vis.obsm['MT'].toarray()
        # Remove MT genes
        adata_vis = adata_vis[:, ~adata_vis.var['MT_gene'].values].copy()

    # Find intersection with signatures
    shared_genes = np.intersect1d(adata_vis.var_names, signatures.index)
    logger.info(f"  Shared genes: {len(shared_genes)} / {adata_vis.n_vars}")

    # Filter to shared genes
    adata_vis = adata_vis[:, shared_genes].copy()
    signatures_filtered = signatures.loc[shared_genes, :].copy()

    logger.info(f"  Final spatial shape: {adata_vis.shape}")
    logger.info(f"  Final signatures shape: {signatures_filtered.shape}")

    return adata_vis, signatures_filtered


def run_cell2location(
    adata_vis,
    signatures,
    output_dir: Path,
    n_cells_per_location: int = 5,
    detection_alpha: int = 200,
    max_epochs: int = 30000,
    use_gpu: bool = True,
):
    """
    Run Cell2Location deconvolution.

    Args:
        adata_vis: Preprocessed spatial AnnData
        signatures: Filtered cell type signatures
        output_dir: Output directory
        n_cells_per_location: Prior on cells per spot
        detection_alpha: Detection efficiency parameter
        max_epochs: Maximum training epochs
        use_gpu: Whether to use GPU

    Returns:
        Annotated AnnData with deconvolution results
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("Setting up Cell2Location model...")
    logger.info(f"  N cells per location: {n_cells_per_location}")
    logger.info(f"  Detection alpha: {detection_alpha}")
    logger.info(f"  Max epochs: {max_epochs}")

    # Setup anndata
    Cell2location.setup_anndata(adata=adata_vis)

    # Create model
    mod = Cell2location(
        adata_vis,
        cell_state_df=signatures,
        N_cells_per_location=n_cells_per_location,
        detection_alpha=detection_alpha,
    )

    # Train (use accelerator/devices instead of deprecated use_gpu)
    logger.info("Training Cell2Location model...")
    start_time = time.time()

    if use_gpu:
        mod.train(
            max_epochs=max_epochs,
            batch_size=None,
            train_size=1,
            accelerator="gpu",
        )
    else:
        mod.train(
            max_epochs=max_epochs,
            batch_size=None,
            train_size=1,
            accelerator="cpu",
        )

    train_time = time.time() - start_time
    logger.info(f"Training completed in {train_time:.1f}s ({train_time/60:.1f} min)")

    # Plot training history
    fig, ax = plt.subplots(1, 1, figsize=(8, 5))
    mod.plot_history(1000, ax=ax)
    plt.savefig(output_dir / "training_history.png", dpi=150, bbox_inches='tight')
    plt.close()

    # Export posterior (use accelerator instead of deprecated use_gpu)
    logger.info("Exporting posterior...")
    sample_kwargs = {
        'num_samples': 3000,
        'batch_size': mod.adata.n_obs,
    }
    if use_gpu:
        sample_kwargs['accelerator'] = 'gpu'
    adata_vis = mod.export_posterior(adata_vis, sample_kwargs=sample_kwargs)

    # Save model
    mod.save(str(output_dir / "model"), overwrite=True)

    return adata_vis, mod, train_time


def extract_proportions(adata_vis, output_dir: Path, sample_name: str):
    """
    Extract and save cell type proportions.

    Args:
        adata_vis: AnnData with Cell2Location results
        output_dir: Output directory
        sample_name: Sample/region name

    Returns:
        Proportions DataFrame
    """
    logger.info("Extracting cell type proportions...")

    # Get abundance at 5th percentile
    abundance = adata_vis.obsm['q05_cell_abundance_w_sf'].copy()

    # Normalize to proportions
    total = abundance.sum(axis=1)
    proportions = abundance.div(total, axis=0)

    # Clean column names (remove prefix if present)
    clean_cols = []
    for col in proportions.columns:
        if '_' in col:
            clean_cols.append(col.split('_')[-1])
        else:
            clean_cols.append(col)
    proportions.columns = clean_cols

    # Save in benchmarking-compatible format (finetuned naming convention)
    output_path = output_dir / f"{sample_name}_cell_prop_finetuned_results.csv"
    proportions.to_csv(output_path)
    logger.info(f"  Saved proportions to {output_path}")
    logger.info(f"  Shape: {proportions.shape}")

    return proportions


def extract_gex_layers(adata_vis, mod, output_dir: Path):
    """
    Extract per-cell-type gene expression layers.

    Args:
        adata_vis: AnnData with Cell2Location results
        mod: Trained Cell2Location model
        output_dir: Output directory

    Returns:
        Dictionary of layer DataFrames
    """
    logger.info("Extracting per-cell-type gene expression layers...")

    layers_dir = output_dir / "layers"
    layers_dir.mkdir(parents=True, exist_ok=True)

    # Compute expected expression per cell type
    expected_dict = mod.module.model.compute_expected_per_cell_type(
        mod.samples["post_sample_q05"],
        mod.adata_manager
    )

    # Add to AnnData layers and save
    layer_dfs = {}
    for i, cell_type in enumerate(mod.factor_names_):
        layer_data = expected_dict['mu'][i]

        # Convert to dense if sparse
        if hasattr(layer_data, 'toarray'):
            layer_data = layer_data.toarray()

        # Add to AnnData
        adata_vis.layers[cell_type] = layer_data

        # Create DataFrame
        df = pd.DataFrame(
            layer_data,
            index=adata_vis.obs_names,
            columns=adata_vis.var_names,
        )

        # Save to CSV
        layer_path = layers_dir / f"{cell_type}_layer.csv"
        df.to_csv(layer_path)
        layer_dfs[cell_type] = df

        logger.info(f"  Saved {cell_type}: {df.shape}")

    return layer_dfs


def main():
    parser = argparse.ArgumentParser(
        description="Run Cell2Location on Xenium pseudo-Visium data"
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
        "--ref-dir",
        type=str,
        required=True,
        help="Directory containing trained reference model",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="output_rna_gt",
        help="Output directory for results",
    )
    parser.add_argument(
        "--n-cells-per-location",
        type=int,
        default=5,
        help="Prior on cells per spot (default: 5)",
    )
    parser.add_argument(
        "--max-epochs",
        type=int,
        default=30000,
        help="Maximum training epochs (default: 30000)",
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
    ref_dir = Path(args.ref_dir)
    if not ref_dir.is_absolute():
        ref_dir = script_dir / ref_dir
    output_dir = Path(args.output_dir)
    if not output_dir.is_absolute():
        output_dir = script_dir / output_dir

    sample_name = f"{args.prefix}_region_{args.region_id}"

    logger.info("=" * 60)
    logger.info(f"Cell2Location Deconvolution - Region {args.region_id}")
    logger.info("=" * 60)
    logger.info(f"Input: {input_dir}")
    logger.info(f"Reference: {ref_dir}")
    logger.info(f"Output: {output_dir}")

    # Load reference
    adata_ref, signatures = load_reference_signatures(ref_dir)

    # Load spatial data
    logger.info(f"\nLoading spatial data for region {args.region_id}...")
    gex_path = input_dir / "h5ad_objects" / f"{args.prefix}_region_{args.region_id}_GEX.h5ad"
    adata_vis = sc.read_h5ad(gex_path)
    logger.info(f"  Spatial shape: {adata_vis.shape}")

    # Preprocess
    adata_vis, signatures_filtered = preprocess_spatial_data(adata_vis, signatures)

    # Run Cell2Location
    region_output_dir = output_dir / sample_name
    adata_vis, mod, train_time = run_cell2location(
        adata_vis,
        signatures_filtered,
        region_output_dir,
        n_cells_per_location=args.n_cells_per_location,
        max_epochs=args.max_epochs,
        use_gpu=not args.no_gpu,
    )

    # Extract proportions
    proportions = extract_proportions(adata_vis, region_output_dir, sample_name)

    # Extract GEX layers
    layer_dfs = extract_gex_layers(adata_vis, mod, region_output_dir)

    # Save results summary
    stats = {
        "region_id": args.region_id,
        "sample_name": sample_name,
        "n_spots": adata_vis.n_obs,
        "n_genes": adata_vis.n_vars,
        "n_cell_types": len(signatures_filtered.columns),
        "cell_types": signatures_filtered.columns.tolist(),
        "train_time_sec": train_time,
        "max_epochs": args.max_epochs,
        "n_cells_per_location": args.n_cells_per_location,
    }

    with open(region_output_dir / "run_stats.json", "w") as f:
        json.dump(stats, f, indent=2)

    # Save annotated AnnData
    adata_vis.write(region_output_dir / f"{sample_name}_annotated.h5ad")

    logger.info("\n" + "=" * 60)
    logger.info("Deconvolution Complete!")
    logger.info("=" * 60)
    logger.info(f"  Spots: {adata_vis.n_obs}")
    logger.info(f"  Cell types: {len(signatures_filtered.columns)}")
    logger.info(f"  Train time: {train_time:.1f}s ({train_time/60:.1f} min)")
    logger.info(f"  Proportions: {region_output_dir / f'{sample_name}_cell_prop_finetuned_results.csv'}")
    logger.info(f"  GEX layers: {region_output_dir / 'layers/'}")


if __name__ == "__main__":
    main()
