"""
Tangram reference-sensitivity benchmark.

Runs Tangram spatial deconvolution with either full or downsampled reference
on simulated Visium data. Uses the same hyperparameters as the main simulation
benchmark (100 markers/type, 500 epochs, rna_count_based density prior).

Usage:
    python -u run_tangram_refsens.py \
        --replicate-id 0 \
        --condition high_seg \
        --ref-condition full \
        --ref-dir ./reference \
        --output-dir ./results/Tangram/high_seg/full
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

try:
    import tangram as tg
except ImportError:
    tg = None

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)

REPLICATES_BASE = (
    str(_cg.REPO_ROOT / "Benchmarking" / "simulation_benchmarking" / "replicates")
)


def load_reference(ref_dir: str):
    """Load reference from h5ad saved by create_downsampled_ref.py."""
    h5ad_path = os.path.join(ref_dir, "reference.h5ad")
    logger.info(f"Loading reference from {h5ad_path} ...")
    adata = sc.read_h5ad(h5ad_path)

    # Preprocess (same as main Tangram benchmark)
    sc.pp.filter_cells(adata, min_genes=10)
    sc.pp.filter_genes(adata, min_cells=3)
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)

    logger.info(f"  Reference shape: {adata.shape}")
    logger.info(f"  Cell types: {adata.obs['cell_type'].value_counts().to_dict()}")
    return adata


def find_marker_genes(adata_ref, labels_key="cell_type", n_markers=100):
    """Find marker genes per cell type using Wilcoxon rank-sum test."""
    logger.info(f"Finding {n_markers} marker genes per cell type ...")
    adata = adata_ref.copy()
    sc.tl.rank_genes_groups(adata, groupby=labels_key, method="wilcoxon")

    markers = set()
    for ct in adata.obs[labels_key].unique():
        ct_markers = sc.get.rank_genes_groups_df(adata, group=ct)
        top = ct_markers.head(n_markers)["names"].tolist()
        markers.update(top)

    marker_list = list(markers)
    logger.info(f"  Total unique markers: {len(marker_list)}")
    return marker_list


def main():
    parser = argparse.ArgumentParser(description="Tangram refsens benchmark")
    parser.add_argument("--replicate-id", type=int, required=True)
    parser.add_argument("--condition", type=str, required=True, choices=["high_seg", "mixed"])
    parser.add_argument("--ref-condition", type=str, required=True, choices=["full", "downsampled"])
    parser.add_argument(
        "--ref-dir",
        type=str,
        default=str(_cg.REPO_ROOT / "Benchmarking" / "simulation_benchmarking" / "refsens" / "reference"),
    )
    parser.add_argument("--output-dir", type=str, required=True)
    parser.add_argument("--n-markers", type=int, default=100)
    parser.add_argument("--num-epochs", type=int, default=500)
    parser.add_argument("--no-gpu", action="store_true")
    args = parser.parse_args()

    device = "cpu" if args.no_gpu else "cuda:0"
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    ref_path = os.path.join(args.ref_dir, args.ref_condition)

    logger.info("=" * 60)
    logger.info(
        f"Tangram Refsens: rep={args.replicate_id} "
        f"condition={args.condition} ref={args.ref_condition}"
    )
    logger.info("=" * 60)

    # Load spatial data
    spatial_path = os.path.join(
        REPLICATES_BASE,
        args.condition,
        "h5ad_objects",
        f"Wu_rep_{args.replicate_id}_GEX.h5ad",
    )
    logger.info(f"Loading spatial data: {spatial_path}")
    adata_sp = sc.read_h5ad(spatial_path)
    sc.pp.normalize_total(adata_sp)
    sc.pp.log1p(adata_sp)
    logger.info(f"  Spatial shape: {adata_sp.shape}")

    # Load reference
    adata_ref = load_reference(ref_path)

    # Marker genes
    markers = find_marker_genes(adata_ref, n_markers=args.n_markers)

    # Preprocessing for Tangram
    shared_genes = np.intersect1d(
        np.intersect1d(adata_ref.var_names, adata_sp.var_names),
        markers,
    )
    logger.info(f"Shared marker genes: {len(shared_genes)}")

    if tg is not None:
        tg.pp_adatas(adata_ref, adata_sp, genes=list(shared_genes))

    # Run mapping
    logger.info(f"Running Tangram mapping ({args.num_epochs} epochs on {device}) ...")
    start_time = time.time()

    ad_map = tg.map_cells_to_space(
        adata_sc=adata_ref,
        adata_sp=adata_sp,
        mode="cells",
        density_prior="rna_count_based",
        num_epochs=args.num_epochs,
        device=device,
    )

    train_time = time.time() - start_time
    logger.info(f"Mapping completed in {train_time:.1f}s")

    # Extract proportions
    tg.project_cell_annotations(ad_map, adata_sp, annotation="cell_type")

    if "tangram_ct_pred" in adata_sp.obsm:
        proportions = pd.DataFrame(
            adata_sp.obsm["tangram_ct_pred"],
            index=adata_sp.obs_names,
        )
    else:
        cell_types = adata_ref.obs["cell_type"]
        mapping_matrix = ad_map.X
        proportions = pd.DataFrame(
            index=adata_sp.obs_names,
            columns=cell_types.unique(),
        )
        for ct in cell_types.unique():
            ct_mask = (cell_types == ct).values
            proportions[ct] = mapping_matrix[ct_mask, :].sum(axis=0)

    proportions = proportions.div(proportions.sum(axis=1), axis=0).fillna(0)

    # Save
    prop_path = output_dir / "cell_type_proportions.csv"
    proportions.to_csv(prop_path)
    logger.info(f"Saved proportions: {prop_path}")

    # Run stats
    stats = {
        "replicate_id": args.replicate_id,
        "condition": args.condition,
        "ref_condition": args.ref_condition,
        "n_spots": adata_sp.n_obs,
        "n_genes_shared": len(shared_genes),
        "n_ref_cells": adata_ref.n_obs,
        "train_time_sec": train_time,
    }
    with open(output_dir / "run_stats.json", "w") as f:
        json.dump(stats, f, indent=2)

    # Write .done marker
    (output_dir / "tangram.done").write_text(str(time.time()))

    logger.info("Done.")


if __name__ == "__main__":
    main()
