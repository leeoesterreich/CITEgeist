"""
Cell2Location reference-sensitivity benchmark.

Runs Cell2Location spatial deconvolution with either full or downsampled
reference on simulated Visium data. Uses the same hyperparameters as the
main simulation benchmark (N_cells_per_location=5, detection_alpha=200,
max_epochs=30000).

The reference model (RegressionModel) is trained from scratch for each
reference condition, since C2L requires its own trained signatures.
Reference training uses 250 epochs (same as train_reference.py).

Usage:
    python -u run_c2l_refsens.py \
        --replicate-id 0 \
        --condition high_seg \
        --ref-condition full \
        --ref-dir ./reference \
        --output-dir ./results/Cell2Location/high_seg/full
"""

import sys as _sys
from pathlib import Path as _Path

_sys.path.insert(0, str(_Path(__file__).resolve().parents[2]))
import _paths as _cg  # noqa: E402

import argparse
import json
import logging
import os
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

import cell2location
from cell2location.models import RegressionModel

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)

REPLICATES_BASE = (
    str(_cg.REPO_ROOT / "Benchmarking" / "simulation_benchmarking" / "replicates")
)


def train_or_load_reference(ref_path: str, output_dir: Path, use_gpu: bool = True):
    """
    Train C2L RegressionModel on reference, or load cached signatures.

    Returns inf_aver DataFrame (genes x cell_types).
    """
    sig_cache = output_dir / "cell_type_signatures.csv"
    model_cache = output_dir / "ref_model"

    if sig_cache.exists():
        logger.info(f"Loading cached signatures from {sig_cache}")
        inf_aver = pd.read_csv(sig_cache, index_col=0)
        return inf_aver

    # Load raw reference
    h5ad_path = os.path.join(ref_path, "reference.h5ad")
    logger.info(f"Loading reference from {h5ad_path}")
    adata_ref = sc.read_h5ad(h5ad_path)

    # C2L needs raw counts; our reference h5ad has raw counts
    logger.info(f"  Reference shape: {adata_ref.shape}")
    logger.info(f"  Cell types: {adata_ref.obs['cell_type'].value_counts().to_dict()}")

    # Add dummy batch key (single-batch reference)
    adata_ref.obs["sample"] = "batch_0"

    # Setup + train RegressionModel
    logger.info("Training RegressionModel (250 epochs) ...")
    RegressionModel.setup_anndata(
        adata=adata_ref,
        batch_key="sample",
        labels_key="cell_type",
    )
    mod = RegressionModel(adata_ref)

    accelerator = "gpu" if use_gpu else "cpu"
    mod.train(max_epochs=250, accelerator=accelerator)

    # Export posterior
    sample_kwargs = {"num_samples": 1000, "batch_size": 2500}
    if use_gpu:
        sample_kwargs["accelerator"] = "gpu"
    adata_ref = mod.export_posterior(adata_ref, sample_kwargs=sample_kwargs)

    # Extract signatures
    factor_names = adata_ref.uns["mod"]["factor_names"]
    if "means_per_cluster_mu_fg" in adata_ref.varm.keys():
        sig_cols = [f"means_per_cluster_mu_fg_{ct}" for ct in factor_names]
        inf_aver = adata_ref.varm["means_per_cluster_mu_fg"][sig_cols].copy()
    else:
        sig_cols = [f"means_per_cluster_mu_fg_{ct}" for ct in factor_names]
        inf_aver = adata_ref.var[sig_cols].copy()

    inf_aver.columns = factor_names

    # Cache
    output_dir.mkdir(parents=True, exist_ok=True)
    inf_aver.to_csv(sig_cache)
    mod.save(str(model_cache), overwrite=True)
    logger.info(f"Saved signatures: {sig_cache} ({inf_aver.shape})")

    return inf_aver


def run_c2l_mapping(
    adata_sp,
    inf_aver,
    output_dir: Path,
    replicate_id: int,
    use_gpu: bool = True,
):
    """Run Cell2Location spatial mapping on a single replicate."""

    # Remove MT genes
    adata_sp.var["SYMBOL"] = adata_sp.var.index
    adata_sp.var["MT_gene"] = [g.startswith("MT-") for g in adata_sp.var["SYMBOL"]]
    adata_sp.obsm["MT"] = adata_sp[:, adata_sp.var["MT_gene"].values].X.toarray()
    adata_sp = adata_sp[:, ~adata_sp.var["MT_gene"].values].copy()

    # Intersect genes
    intersect = np.intersect1d(adata_sp.var_names, inf_aver.index)
    adata_sp = adata_sp[:, intersect].copy()
    inf_aver_shared = inf_aver.loc[intersect, :].copy()

    logger.info(f"  Shared genes: {len(intersect)}")

    # Setup and train spatial model
    cell2location.models.Cell2location.setup_anndata(adata=adata_sp)

    mod = cell2location.models.Cell2location(
        adata_sp,
        cell_state_df=inf_aver_shared,
        N_cells_per_location=5,
        detection_alpha=200,
    )

    logger.info("Training Cell2location spatial model (30000 epochs) ...")
    start_time = time.time()
    accelerator = "gpu" if use_gpu else "cpu"
    mod.train(max_epochs=30000, batch_size=None, train_size=1, accelerator=accelerator)
    train_time = time.time() - start_time
    logger.info(f"Training completed in {train_time:.0f}s")

    # Export posterior
    sample_kwargs = {
        "num_samples": 3000,
        "batch_size": mod.adata.n_obs,
        "accelerator": accelerator,
    }
    adata_sp = mod.export_posterior(adata_sp, sample_kwargs=sample_kwargs)

    # Extract proportions from q05 abundances
    df = adata_sp.obsm["q05_cell_abundance_w_sf"]
    total_abundance = df.sum(axis=1)
    proportions = df.div(total_abundance, axis=0)
    proportions.columns = [col.split("_")[-1] for col in proportions.columns]
    proportions.reset_index(inplace=True)
    proportions.rename(columns={"index": "spot"}, inplace=True)

    return proportions, train_time


def main():
    parser = argparse.ArgumentParser(description="Cell2Location refsens benchmark")
    parser.add_argument("--replicate-id", type=int, required=True)
    parser.add_argument("--condition", type=str, required=True, choices=["high_seg", "mixed"])
    parser.add_argument("--ref-condition", type=str, required=True, choices=["full", "downsampled"])
    parser.add_argument(
        "--ref-dir",
        type=str,
        default=str(_cg.REPO_ROOT / "Benchmarking" / "simulation_benchmarking" / "refsens" / "reference"),
    )
    parser.add_argument("--output-dir", type=str, required=True)
    parser.add_argument("--no-gpu", action="store_true")
    args = parser.parse_args()

    use_gpu = not args.no_gpu
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    ref_path = os.path.join(args.ref_dir, args.ref_condition)

    logger.info("=" * 60)
    logger.info(
        f"Cell2Location Refsens: rep={args.replicate_id} "
        f"condition={args.condition} ref={args.ref_condition}"
    )
    logger.info("=" * 60)

    # Train / load reference signatures
    # Cache per ref_condition so training happens once, not per replicate
    ref_model_dir = Path(args.ref_dir) / args.ref_condition / "c2l_model"
    inf_aver = train_or_load_reference(ref_path, ref_model_dir, use_gpu=use_gpu)

    # Load spatial data
    spatial_path = os.path.join(
        REPLICATES_BASE,
        args.condition,
        "h5ad_objects",
        f"Wu_rep_{args.replicate_id}_GEX.h5ad",
    )
    logger.info(f"Loading spatial data: {spatial_path}")
    adata_sp = sc.read_h5ad(spatial_path)
    logger.info(f"  Spatial shape: {adata_sp.shape}")

    # Run mapping
    proportions, train_time = run_c2l_mapping(
        adata_sp, inf_aver, output_dir, args.replicate_id, use_gpu=use_gpu
    )

    # Save proportions
    prop_path = output_dir / "cell2loc_deconv_predictions.csv"
    proportions.to_csv(prop_path, index=False)
    logger.info(f"Saved proportions: {prop_path}")

    # Run stats
    stats = {
        "replicate_id": args.replicate_id,
        "condition": args.condition,
        "ref_condition": args.ref_condition,
        "train_time_sec": train_time,
    }
    with open(output_dir / "run_stats.json", "w") as f:
        json.dump(stats, f, indent=2)

    # Write .done marker
    (output_dir / "c2l.done").write_text(str(time.time()))

    logger.info("Done.")


if __name__ == "__main__":
    main()
