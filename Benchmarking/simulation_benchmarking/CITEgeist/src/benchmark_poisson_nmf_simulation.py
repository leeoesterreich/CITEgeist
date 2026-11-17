#!/usr/bin/env python
"""Poisson NMF GEX deconvolution benchmark on scCube simulation data.

Runs CITEgeist proportion estimation → Poisson NMF GEX allocation
in both profile_only and joint modes, then evaluates against GT.

Usage:
    python benchmark_poisson_nmf_simulation.py --replicate 0 --condition high_seg
    python benchmark_poisson_nmf_simulation.py --replicate 0 --condition high_seg --mode joint
"""

import sys as _sys
from pathlib import Path as _Path

_sys.path.insert(0, str(_Path(__file__).resolve().parents[3]))
import _paths as _cg  # noqa: E402

import argparse
import json
import logging
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

# ---------------------------------------------------------------------------
# Repo path setup
# ---------------------------------------------------------------------------
REPO_ROOT = Path(__file__).resolve().parent.parent.parent.parent.parent
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(REPO_ROOT / "Benchmarking" / "simulation_benchmarking" / "src"))

# ---------------------------------------------------------------------------
# Constants (shared with benchmark_sace_simulation.py)
# ---------------------------------------------------------------------------
DEFAULT_SCCUBE_DIR = Path(
    str(_cg.root("DATA") / "Wu_Visium" / "Simulations" / "scCube_12k" / "replicates")
)
DEFAULT_OUTPUT_BASE = REPO_ROOT / "Benchmarking" / "simulation_benchmarking" / "CITEgeist" / "results_poisson_nmf"

SIMULATION_CELL_PROFILE_DICT = {
    "B-cells": {"Major": ["B-cells_Protein_1", "B-cells_Protein_2"]},
    "CAFs": {"Major": ["CAFs_Protein_1", "CAFs_Protein_2"]},
    "Cancer Epithelial": {"Major": ["Cancer Epithelial_Protein_1", "Cancer Epithelial_Protein_2"]},
    "Endothelial": {"Major": ["Endothelial_Protein_1", "Endothelial_Protein_2"]},
    "Myeloid": {"Major": ["Myeloid_Protein_1", "Myeloid_Protein_2"]},
    "Normal Epithelial": {"Major": ["Normal Epithelial_Protein_1", "Normal Epithelial_Protein_2"]},
    "PVL": {"Major": ["PVL_Protein_1", "PVL_Protein_2"]},
    "Plasmablasts": {"Major": ["Plasmablasts_Protein_1", "Plasmablasts_Protein_2"]},
    "T-cells": {"Major": ["T-cells_Protein_1", "T-cells_Protein_2"]},
}

SIM_TYPES = list(SIMULATION_CELL_PROFILE_DICT.keys())
CELLS_PER_SPOT = 5

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Reuse helpers from SACE benchmark
# ---------------------------------------------------------------------------
from benchmark_sace_simulation import (
    largest_remainder_round,
    build_synthetic_cell_spot_map,
    spot_type_gex_to_layers,
    export_layers,
)

# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------


def run_poisson_nmf_benchmark(
    replicate_id: int,
    condition: str,
    mode: str,
    sccube_dir: Path,
    output_base: Path,
    max_workers: int = 8,
    cells_per_spot: int = CELLS_PER_SPOT,
    max_iter: int = 100,
    lambda_W: float = 0.01,
    lambda_H_init: float = 10.0,
    lambda_H_final: float = 0.1,
) -> dict:
    """Run Poisson NMF GEX benchmark for one replicate.

    Args:
        replicate_id: Replicate index (0-4).
        condition: "high_seg" or "mixed".
        mode: "profile_only" or "joint".
        sccube_dir: Root scCube replicates directory.
        output_base: Root output directory.
        max_workers: Parallel workers for proportion estimation.
        cells_per_spot: Synthetic cells per spot.
        max_iter: NMF iterations.
        lambda_W: L2 regularization for W.
        lambda_H_init: Initial KL regularization for H (joint mode).
        lambda_H_final: Final KL regularization for H (joint mode).

    Returns:
        Dict with evaluation metrics.
    """
    from CITEgeist.model.citegeist_model import CitegeistModel
    from CITEgeist.model.poisson_nmf_gex import run_poisson_nmf_gex

    # -----------------------------------------------------------------------
    # Load data (same as SACE benchmark)
    # -----------------------------------------------------------------------
    h5ad_dir = sccube_dir / condition / "h5ad_objects"
    gex_path = h5ad_dir / f"Wu_rep_{replicate_id}_GEX.h5ad"
    cite_path = h5ad_dir / f"Wu_rep_{replicate_id}_CITE.h5ad"

    logger.info("Loading GEX: %s", gex_path)
    adata_gex = sc.read_h5ad(gex_path)
    logger.info("Loading CITE: %s", cite_path)
    adata_cite = sc.read_h5ad(cite_path)

    # -----------------------------------------------------------------------
    # Run CitegeistModel for proportions (via QP)
    # -----------------------------------------------------------------------
    model_output = str(output_base / condition / f"Wu_rep_{replicate_id}" / mode / "citegeist")
    model = CitegeistModel(
        sample_name=f"Wu_rep_{replicate_id}",
        output_folder=model_output,
        simulation=True,
        gene_expression_adata=adata_gex,
        antibody_capture_adata=adata_cite,
    )

    logger.info("Preprocessing GEX (target_sum=10000) ...")
    model.preprocess_gex(target_sum=10000)
    logger.info("Preprocessing antibody ...")
    model.preprocess_antibody()
    model.cell_profile_dict = SIMULATION_CELL_PROFILE_DICT

    logger.info("Running proportion estimation ...")
    t0 = time.time()
    global_props_df, _ = model.run_cell_proportion_model(method="qp", max_workers=max_workers)
    elapsed_prop = time.time() - t0
    logger.info("Proportions done in %.1fs", elapsed_prop)

    # Clean proportions
    props_df = global_props_df.drop(columns=["recon_error"], errors="ignore")
    props_df = props_df.reindex(model.gene_expression_adata.obs_names)
    n_types = props_df.shape[1]
    props_df = props_df.fillna(1.0 / n_types)

    # -----------------------------------------------------------------------
    # Build synthetic cell assignments
    # -----------------------------------------------------------------------
    cell_assignments, cell_spot_map = build_synthetic_cell_spot_map(props_df, cells_per_spot=cells_per_spot)

    # -----------------------------------------------------------------------
    # Get raw counts for Poisson NMF
    # -----------------------------------------------------------------------
    raw_X = np.asarray(adata_gex.X, dtype=float)
    if hasattr(raw_X, "toarray"):
        raw_X = raw_X.toarray()

    # Get antibody data for marker-guided W init
    antibody_data = np.asarray(adata_cite.X, dtype=float)
    if hasattr(antibody_data, "toarray"):
        antibody_data = antibody_data.toarray()
    antibody_names = list(adata_cite.var_names)

    # -----------------------------------------------------------------------
    # Run Poisson NMF
    # -----------------------------------------------------------------------
    logger.info("Running Poisson NMF (mode=%s, max_iter=%d) ...", mode, max_iter)
    t0 = time.time()
    spot_type_gex, cell_adata, diagnostics = run_poisson_nmf_gex(
        spot_counts=raw_X,
        proportions=props_df,
        cell_assignments=cell_assignments,
        cell_spot_map=cell_spot_map,
        gene_names=list(adata_gex.var_names),
        mode=mode,
        max_iter=max_iter,
        lambda_W=lambda_W,
        lambda_H_init=lambda_H_init,
        lambda_H_final=lambda_H_final,
        antibody_data=antibody_data,
        antibody_names=antibody_names,
        cell_profile_dict=SIMULATION_CELL_PROFILE_DICT,
    )
    elapsed_nmf = time.time() - t0
    logger.info("Poisson NMF done in %.1fs (%d iterations)", elapsed_nmf, diagnostics.get("n_iterations", -1))

    # -----------------------------------------------------------------------
    # Convert to layers and export
    # -----------------------------------------------------------------------
    spot_names = list(adata_gex.obs_names)
    gene_names = list(adata_gex.var_names)
    type_names = list(props_df.columns)

    layers = spot_type_gex_to_layers(spot_type_gex, spot_names, gene_names, type_names)
    layers_dir = output_base / condition / f"Wu_rep_{replicate_id}" / mode / "layers"
    export_layers(layers, layers_dir)

    # -----------------------------------------------------------------------
    # Evaluate
    # -----------------------------------------------------------------------
    gt_dir = sccube_dir / condition / "ST_GEX_sim" / f"sample_{replicate_id}" / "layers"
    from benchmarking_gex import calculate_rmse

    metrics = calculate_rmse(
        ground_truth_dir=str(gt_dir),
        predictions_dir=str(layers_dir),
        normalize="range",
    )

    # Save metrics
    metrics_path = layers_dir / "poisson_nmf_metrics.json"
    with open(metrics_path, "w") as fh:
        json.dump(
            {
                "replicate": replicate_id,
                "condition": condition,
                "mode": mode,
                "elapsed_prop_s": elapsed_prop,
                "elapsed_nmf_s": elapsed_nmf,
                "max_iter": max_iter,
                "lambda_W": lambda_W,
                "lambda_H_init": lambda_H_init,
                "lambda_H_final": lambda_H_final,
                "n_iterations": diagnostics.get("n_iterations"),
                "final_kl": diagnostics.get("final_kl"),
                **{k: v for k, v in metrics.items() if k != "metrics_per_cell_type"},
                "metrics_per_cell_type": {
                    ct: {mk: float(mv) for mk, mv in cm.items()}
                    for ct, cm in metrics.get("metrics_per_cell_type", {}).items()
                },
            },
            fh,
            indent=2,
            default=str,
        )

    logger.info(
        "SUMMARY rep=%d cond=%s mode=%s | prop=%.1fs nmf=%.1fs | " "spot_layer_r=%.4f  pearson_r=%.4f  nrmse=%.4f",
        replicate_id,
        condition,
        mode,
        elapsed_prop,
        elapsed_nmf,
        metrics.get("mean_spot_layer_r", float("nan")),
        metrics.get("mean_pearson_r", float("nan")),
        metrics.get("average_nrmse", float("nan")),
    )

    return {
        "replicate": replicate_id,
        "condition": condition,
        "mode": mode,
        "elapsed_prop_s": elapsed_prop,
        "elapsed_nmf_s": elapsed_nmf,
        "evaluation": metrics,
        "diagnostics": diagnostics,
    }


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Poisson NMF GEX benchmark on scCube simulation data.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--replicate", type=int, default=0, help="Replicate ID (0-4)")
    parser.add_argument(
        "--condition",
        type=str,
        default="high_seg",
        choices=["high_seg", "mixed"],
        help="Simulation condition",
    )
    parser.add_argument(
        "--mode",
        type=str,
        default="profile_only",
        choices=["profile_only", "joint"],
        help="NMF mode",
    )
    parser.add_argument("--sccube_dir", type=Path, default=DEFAULT_SCCUBE_DIR)
    parser.add_argument("--output_dir", type=Path, default=DEFAULT_OUTPUT_BASE)
    parser.add_argument("--max_workers", type=int, default=8)
    parser.add_argument("--max_iter", type=int, default=100, help="NMF iterations")
    parser.add_argument("--lambda_W", type=float, default=0.01)
    parser.add_argument("--lambda_H_init", type=float, default=10.0)
    parser.add_argument("--lambda_H_final", type=float, default=0.1)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    logger.info(
        "Poisson NMF benchmark | rep=%d cond=%s mode=%s",
        args.replicate,
        args.condition,
        args.mode,
    )
    run_poisson_nmf_benchmark(
        replicate_id=args.replicate,
        condition=args.condition,
        mode=args.mode,
        sccube_dir=args.sccube_dir,
        output_base=args.output_dir,
        max_workers=args.max_workers,
        max_iter=args.max_iter,
        lambda_W=args.lambda_W,
        lambda_H_init=args.lambda_H_init,
        lambda_H_final=args.lambda_H_final,
    )


if __name__ == "__main__":
    main()
