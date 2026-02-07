#!/usr/bin/env python
"""
Experiment: Add Vimentin to Epithelial profile in achievable-7 benchmark.

Hypothesis: Fibroblasts may be absorbing proportion from VIM+ epithelial cells
(EMT-positive) because only fibroblasts claim Vimentin. Adding VIM to the
Epithelial profile should let the optimizer attribute VIM signal to both
epithelial and fibroblast spots based on co-occurring markers (PanCK vs alphaSMA).

Changes from baseline achievable-7:
  - Epithelial: ["PanCK"] → ["PanCK", "Vimentin"]
  - Fibroblasts: ["alphaSMA", "Vimentin"] (unchanged)

Usage:
    python run_benchmark_vim_epithelial.py --region-id 0
    python run_benchmark_vim_epithelial.py --region-id 0 --run-gex
"""

import argparse
import json
import logging
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

# Add CITEgeist to path
REPO_ROOT = Path(__file__).parent.parent.parent.parent.parent
sys.path.insert(0, str(REPO_ROOT))

from CITEgeist.model.citegeist_model import CitegeistModel

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)

# =============================================================================
# MODIFIED ACHIEVABLE-7 PROFILES: Vimentin added to Epithelial
# =============================================================================

VIM_EPITHELIAL_CELL_PROFILE_DICT = {
    "B cells": {
        "Major": ["CD20"],
        "Minor": ["CD45RA"],
    },
    "CD4+ T cells": {
        "Major": ["CD3E", "CD4"],
        "Minor": ["CD45RO"],
    },
    "CD8+ T cells": {
        "Major": ["CD3E", "CD8A"],
        "Minor": ["GranzymeB"],
    },
    "Macrophages": {
        "Major": ["CD68", "CD163"],
        "Minor": ["CD16"],
    },
    "Endothelial": {
        "Major": ["CD31"],
        "Minor": [],
    },
    "Epithelial": {
        "Major": ["PanCK", "Vimentin"],  # <-- CHANGED: added Vimentin
        "Minor": [],
    },
    "Fibroblasts": {
        "Major": ["alphaSMA", "Vimentin"],
        "Minor": [],
    },
}


def run_vim_epithelial_benchmark(
    region_id: int,
    input_dir: Path,
    output_dir: Path,
    radius: float = 4.0,
    lambda_reg: float = 1.0,
    alpha_elastic: float = 0.7,
    max_y_change: float = 0.4,
    min_counts: int = 25,
    run_gex: bool = False,
) -> dict:
    """Run benchmark with VIM added to Epithelial profile."""

    logger.info("=" * 70)
    logger.info(f"VIM-EPITHELIAL EXPERIMENT: Region {region_id}")
    logger.info("Epithelial profile: PanCK + Vimentin (shared with Fibroblasts)")
    logger.info("=" * 70)

    # Load data
    gex_path = input_dir / "h5ad_objects" / f"Xenium_region_{region_id}_GEX.h5ad"
    protein_path = input_dir / "h5ad_objects" / f"Xenium_region_{region_id}_CITE.h5ad"

    for path in [gex_path, protein_path]:
        if not path.exists():
            logger.error(f"File not found: {path}")
            sys.exit(1)

    logger.info(f"Loading data for region {region_id}...")
    gex_adata = sc.read_h5ad(gex_path)
    protein_adata = sc.read_h5ad(protein_path)

    sample_name = f"Xenium_region_{region_id}"

    model = CitegeistModel(
        sample_name=sample_name,
        output_folder=str(output_dir),
        simulation=True,
        gene_expression_adata=gex_adata,
        antibody_capture_adata=protein_adata,
    )

    model.load_cell_profile_dict(VIM_EPITHELIAL_CELL_PROFILE_DICT)

    logger.info("Preprocessing...")
    model.filter_gex(min_counts=min_counts)
    model.preprocess_gex(target_sum=10000)
    model.preprocess_antibody()

    logger.info("Running cell proportion optimization...")
    start_time = time.time()

    global_props, finetuned_props = model.run_cell_proportion_model(
        radius=radius,
        lambda_reg=lambda_reg,
        alpha=alpha_elastic,
        max_y_change=max_y_change,
        validation_warn_only=True,
    )

    prop_time = time.time() - start_time
    logger.info(f"Cell proportion optimization completed in {prop_time:.1f}s")

    # Save results
    result_dir = output_dir / sample_name
    result_dir.mkdir(parents=True, exist_ok=True)

    if finetuned_props is not None:
        props_df = pd.DataFrame(
            finetuned_props,
            index=model.antibody_capture_adata.obs_names,
            columns=list(VIM_EPITHELIAL_CELL_PROFILE_DICT.keys()),
        )
        props_df.to_csv(result_dir / f"{sample_name}_deconv_predictions.csv")

    # Run GEX if requested
    gex_time = None
    if run_gex:
        logger.info("Running gene expression deconvolution...")
        gex_start = time.time()
        try:
            model.run_cell_expression_pass1(
                radius=radius,
                alpha=0.5,
                checkpoint_interval=100,
                output_dir=str(output_dir / "checkpoints"),
                rerun=True,
            )
            gex_time = time.time() - gex_start
            logger.info(f"GEX deconvolution completed in {gex_time:.1f}s")
        except Exception as e:
            logger.error(f"GEX deconvolution failed: {e}")
            gex_time = -1.0

    results = {
        "region_id": region_id,
        "mode": "vim_epithelial_experiment",
        "profile_change": "Added Vimentin to Epithelial Major markers",
        "cell_profile_dict": VIM_EPITHELIAL_CELL_PROFILE_DICT,
        "n_spots": gex_adata.shape[0],
        "n_cell_types": len(VIM_EPITHELIAL_CELL_PROFILE_DICT),
        "cell_types": list(VIM_EPITHELIAL_CELL_PROFILE_DICT.keys()),
        "runtime_proportions_sec": prop_time,
        "runtime_gex_sec": gex_time,
        "output_dir": str(result_dir),
    }

    with open(result_dir / "run_stats.json", "w") as f:
        json.dump(results, f, indent=2, default=str)

    logger.info(f"Results saved to {result_dir}")
    return results


def main():
    parser = argparse.ArgumentParser(
        description="VIM-Epithelial experiment: benchmark with Vimentin in Epithelial profile"
    )
    parser.add_argument("--region-id", type=int, required=True)
    parser.add_argument(
        "--input-dir",
        type=str,
        default=str(REPO_ROOT / "Benchmarking/xenium_pseudovisium/data_granular_gt"),
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default=str(
            REPO_ROOT / "Benchmarking/xenium_benchmarking/CITEgeist/output/vim_epithelial"
        ),
    )
    parser.add_argument("--radius", type=float, default=4.0)
    parser.add_argument("--lambda-reg", type=float, default=1.0)
    parser.add_argument("--alpha-elastic", type=float, default=0.7)
    parser.add_argument("--max-y-change", type=float, default=0.4)
    parser.add_argument("--min-counts", type=int, default=25)
    parser.add_argument("--run-gex", action="store_true")

    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    results = run_vim_epithelial_benchmark(
        region_id=args.region_id,
        input_dir=Path(args.input_dir),
        output_dir=output_dir,
        radius=args.radius,
        lambda_reg=args.lambda_reg,
        alpha_elastic=args.alpha_elastic,
        max_y_change=args.max_y_change,
        min_counts=args.min_counts,
        run_gex=args.run_gex,
    )

    print(f"\nCompleted region {args.region_id} (vim_epithelial experiment)")
    print(f"  Spots: {results['n_spots']}")
    print(f"  Cell types: {results['n_cell_types']}")
    print(f"  Runtime: {results['runtime_proportions_sec']:.1f}s")
    print(f"  Output: {results['output_dir']}")


if __name__ == "__main__":
    main()
