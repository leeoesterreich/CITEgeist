"""
Train Cell2Location reference model on scRNA-seq data.

Trains a RegressionModel on a reference scRNA-seq dataset to generate
cell-type-specific gene expression signatures for spatial deconvolution.
Shared across simulation and Xenium benchmarking pipelines.

Usage:
    python train_reference.py --ref-path /path/to/reference.h5ad --output-dir /path/to/output
    python train_reference.py --ref-path ref.h5ad --output-dir out/ --labels-key celltype_major
"""

import os
import sys
import argparse
import hashlib
import json
import logging
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

import cell2location
from cell2location.models import RegressionModel

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)


def filter_genes(
    adata,
    cell_count_cutoff: int = 5,
    cell_percentage_cutoff: float = 0.03,
    nonz_mean_cutoff: float = 1.12,
):
    """
    Filter genes based on expression criteria.

    Args:
        adata: AnnData object
        cell_count_cutoff: Minimum cells expressing gene
        cell_percentage_cutoff: Minimum percentage of cells
        nonz_mean_cutoff: Minimum mean expression in expressing cells

    Returns:
        Boolean mask of selected genes
    """
    # Calculate statistics
    cell_count = np.array((adata.X > 0).sum(axis=0)).flatten()
    cell_percentage = cell_count / adata.n_obs

    # Mean expression in expressing cells
    X = adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X
    nonz_mean = np.array([X[X[:, i] > 0, i].mean() if (X[:, i] > 0).any() else 0 for i in range(X.shape[1])])

    # Apply filters
    selected = (
        (cell_count >= cell_count_cutoff)
        & (cell_percentage >= cell_percentage_cutoff)
        & (nonz_mean >= nonz_mean_cutoff)
    )

    logger.info(f"Gene filtering: {selected.sum()} / {len(selected)} genes passed")
    return selected


def read_gene_list(path) -> list:
    """One gene per line; blanks stripped; duplicates and an empty file are errors."""
    genes = [ln.strip() for ln in Path(path).read_text().splitlines() if ln.strip()]
    if not genes:
        raise ValueError(f"gene list {path} is empty")
    dup = sorted({g for g in genes if genes.count(g) > 1})
    if dup:
        raise ValueError(f"gene list {path} has duplicates, e.g. {dup[:10]}")
    return genes


def subset_to_gene_list(adata, genes):
    """Keep exactly ``genes`` (in adata.var_names order). Raises ValueError naming absent genes."""
    present = set(adata.var_names)
    missing = [g for g in genes if g not in present]
    if missing:
        raise ValueError(f"{len(missing)} gene(s) not in the reference, e.g. {missing[:10]}")
    wanted = set(genes)
    out = adata[:, [g in wanted for g in adata.var_names]].copy()
    if set(out.var_names) != wanted:
        raise AssertionError("gene subset does not equal the gene list")
    return out


def train_reference_model(
    adata_ref,
    output_dir: Path,
    batch_key: str = "sample",
    labels_key: str = "cell_type",
    max_epochs: int = 250,
    use_gpu: bool = True,
):
    """
    Train Cell2Location RegressionModel on reference data.

    Args:
        adata_ref: Reference AnnData with cell type labels
        output_dir: Output directory for model and signatures
        batch_key: Column for batch effects
        labels_key: Column for cell type labels
        max_epochs: Maximum training epochs
        use_gpu: Whether to use GPU

    Returns:
        Trained model and annotated AnnData
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("Setting up RegressionModel...")
    logger.info(f"  Batch key: {batch_key}")
    logger.info(f"  Labels key: {labels_key}")
    logger.info(f"  Cell types: {adata_ref.obs[labels_key].unique().tolist()}")

    # Setup anndata for Cell2Location
    RegressionModel.setup_anndata(
        adata=adata_ref,
        batch_key=batch_key,
        labels_key=labels_key,
    )

    # Create model
    mod = RegressionModel(adata_ref)

    # Train (use accelerator instead of deprecated use_gpu)
    logger.info(f"Training RegressionModel for {max_epochs} epochs...")
    if use_gpu:
        mod.train(max_epochs=max_epochs, accelerator="gpu")
    else:
        mod.train(max_epochs=max_epochs, accelerator="cpu")

    # Plot training history
    logger.info("Saving training history...")
    fig, ax = plt.subplots(1, 1, figsize=(8, 5))
    mod.plot_history(ax=ax)
    plt.savefig(output_dir / "training_history.png", dpi=150, bbox_inches="tight")
    plt.close()

    # Export posterior (use accelerator instead of deprecated use_gpu)
    logger.info("Exporting posterior...")
    sample_kwargs = {
        "num_samples": 1000,
        "batch_size": 2500,
    }
    if use_gpu:
        sample_kwargs["accelerator"] = "gpu"
    adata_ref = mod.export_posterior(adata_ref, sample_kwargs=sample_kwargs)

    # Save model
    logger.info("Saving model...")
    mod.save(str(output_dir / "model"), overwrite=True)

    # Save reference with posterior
    adata_ref.write(output_dir / "reference_trained.h5ad")

    # Extract and save signatures
    logger.info("Extracting cell type signatures...")

    # Get factor names (cell types)
    factor_names = adata_ref.uns["mod"]["factor_names"]
    logger.info(f"  Factor names: {factor_names}")

    # Extract mean gene expression per cell type
    if "means_per_cluster_mu_fg" in adata_ref.varm.keys():
        signature_cols = [f"means_per_cluster_mu_fg_{ct}" for ct in factor_names]
        inf_aver = adata_ref.varm["means_per_cluster_mu_fg"][signature_cols].copy()
    else:
        signature_cols = [f"means_per_cluster_mu_fg_{ct}" for ct in factor_names]
        inf_aver = adata_ref.var[signature_cols].copy()

    inf_aver.columns = factor_names

    # Save signatures
    inf_aver.to_csv(output_dir / "cell_type_signatures.csv")
    logger.info(f"  Saved signatures: {inf_aver.shape}")

    return mod, adata_ref, inf_aver


def main():
    parser = argparse.ArgumentParser(description="Train Cell2Location reference model")
    parser.add_argument(
        "--ref-path",
        type=str,
        required=True,
        help="Path to processed reference h5ad file",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        required=True,
        help="Output directory for model and signatures",
    )
    parser.add_argument(
        "--batch-key",
        type=str,
        default="sample",
        help="Column for batch effects (default: sample)",
    )
    parser.add_argument(
        "--labels-key",
        type=str,
        default="cell_type",
        help="Column for cell type labels (default: cell_type)",
    )
    parser.add_argument(
        "--max-epochs",
        type=int,
        default=250,
        help="Maximum training epochs (default: 250)",
    )
    parser.add_argument(
        "--filter-genes",
        action="store_true",
        help="Apply gene filtering before training",
    )
    parser.add_argument(
        "--gene-list",
        type=str,
        default=None,
        help="Train on exactly these genes (one per line); mutually exclusive with --filter-genes",
    )
    parser.add_argument(
        "--drop-genes",
        type=str,
        default=None,
        help="Optional file (one gene per line) removed from the reference before any gene filtering. Off by default.",
    )
    parser.add_argument(
        "--no-gpu",
        action="store_true",
        help="Disable GPU usage",
    )

    args = parser.parse_args()
    if args.gene_list and args.filter_genes:
        parser.error("--gene-list and --filter-genes are mutually exclusive")

    ref_path = Path(args.ref_path)
    output_dir = Path(args.output_dir)

    logger.info("=" * 60)
    logger.info("Cell2Location Reference Training")
    logger.info("=" * 60)
    logger.info(f"Reference: {ref_path}")
    logger.info(f"Output: {output_dir}")

    # Load reference
    logger.info("\nLoading reference data...")
    adata_ref = sc.read_h5ad(ref_path)
    logger.info(f"  Shape: {adata_ref.shape}")
    logger.info(f"  Cell types: {adata_ref.obs[args.labels_key].value_counts().to_dict()}")

    drop_list = None
    if args.drop_genes:  # opt-in gene-removal control: applied before ANY gene selection
        sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
        from Benchmarking.xenium_benchmarking.rna_refinement.drop_genes import drop_genes, read_gene_file

        drop_list = read_gene_file(args.drop_genes)
        adata_ref = drop_genes(adata_ref, drop_list)
        logger.info(f"  After --drop-genes: {adata_ref.shape} {adata_ref.uns['drop_genes']['summary']}")

    # Filter genes if requested
    if args.filter_genes:
        logger.info("\nFiltering genes...")
        selected = filter_genes(adata_ref)
        adata_ref = adata_ref[:, selected].copy()
        logger.info(f"  After filtering: {adata_ref.shape}")

    gene_list = None
    if args.gene_list:
        gene_list = read_gene_list(args.gene_list)
        adata_ref = subset_to_gene_list(adata_ref, gene_list)
        logger.info(f"  Subset to gene list ({len(gene_list)} genes): {adata_ref.shape}")

    # Ensure raw counts (Cell2Location needs counts, not normalized)
    # If .raw exists, use it
    if adata_ref.raw is not None:
        logger.info("Using raw counts from adata.raw...")
        adata_ref = adata_ref.raw.to_adata()
        adata_ref.obs[args.labels_key] = adata_ref.obs[args.labels_key] if args.labels_key in adata_ref.obs else None
        if gene_list is not None:  # the .raw swap restores the full gene set; re-apply the list
            adata_ref = subset_to_gene_list(adata_ref, gene_list)
        if drop_list is not None:  # same reason: re-apply the removal
            adata_ref = drop_genes(adata_ref, drop_list)

    # Check for required columns
    if args.batch_key not in adata_ref.obs.columns:
        logger.warning(f"Batch key '{args.batch_key}' not found. Creating dummy batch...")
        adata_ref.obs[args.batch_key] = "batch_0"

    # Train model
    logger.info("\n" + "=" * 60)
    logger.info("Training Reference Model")
    logger.info("=" * 60)

    mod, adata_ref, signatures = train_reference_model(
        adata_ref,
        output_dir,
        batch_key=args.batch_key,
        labels_key=args.labels_key,
        max_epochs=args.max_epochs,
        use_gpu=not args.no_gpu,
    )

    logger.info("\n" + "=" * 60)
    logger.info("Training Complete!")
    logger.info("=" * 60)
    logger.info(f"Model saved to: {output_dir / 'model'}")
    logger.info(f"Reference saved to: {output_dir / 'reference_trained.h5ad'}")
    logger.info(f"Signatures saved to: {output_dir / 'cell_type_signatures.csv'}")
    logger.info(f"Signature shape: {signatures.shape}")

    if gene_list is not None:
        if set(signatures.index) != set(gene_list):
            raise AssertionError("trained signature genes differ from the gene list")
        (output_dir / "training_genes.json").write_text(
            json.dumps(
                {
                    "n_genes": len(gene_list),
                    "gene_list_sha256": hashlib.sha256(Path(args.gene_list).read_bytes()).hexdigest(),
                },
                indent=2,
            )
        )


if __name__ == "__main__":
    main()
