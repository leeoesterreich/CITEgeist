#!/usr/bin/env python
"""
Module 4 Program Discovery Runner for Patient Samples.

Runs Module 4 (anchored + joint program discovery) on patient samples
using existing Module 3 deconvolved gene expression outputs.

Usage:
    python run_module4_discovery.py --sample HCC22-088-P1-S1 --output-dir output/module4_discovery

    # Run both anchored and joint discovery:
    python run_module4_discovery.py --sample HCC22-088-P1-S1 --mode both
"""
import argparse
import json
import logging
import os
import sys
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse
import squidpy as sq

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from CITEgeist.model.anchored_program_discovery import (
    discover_anchored_programs,
    discover_joint_programs,
)

logger = logging.getLogger(__name__)
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)

# Data locations
DATA_ROOT = Path("/ix1/alee/LO_LAB/General/Lab_Data/20250210_CITEGeistPublicData_GEO_Alex/processed_files")
# Use unified profile Module 3 outputs (10 cell types)
MODULE3_OUTPUT = Path("/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/output/module3_unified")

# All 12 patient samples (deduplicated: P4-S2_1i_rep replaces P4-S2, P5-S2_F_rep replaces P5-S2)
ALL_SAMPLES = [
    "HCC22-088-P1-S1",
    "HCC22-088-P1-S2",
    "HCC22-088-P2-S1",
    "HCC22-088-P2-S2",
    "HCC22-088-P3-S1_A",
    "HCC22-088-P3-S2",
    "HCC22-088-P4-S1",
    "HCC22-088-P4-S2_1i_rep",
    "HCC22-088-P5-S1",
    "HCC22-088-P5-S2_F_rep",
    "HCC22-088-P6-S1",
    "HCC22-088-P6-S2_D",
]

# Sample metadata
SAMPLE_METADATA = {
    "HCC22-088-P1-S1": {"patient": "P1", "timepoint": "S1", "response": "Progressor"},
    "HCC22-088-P1-S2": {"patient": "P1", "timepoint": "S2", "response": "Progressor"},
    "HCC22-088-P2-S1": {"patient": "P2", "timepoint": "S1", "response": "Responder"},
    "HCC22-088-P2-S2": {"patient": "P2", "timepoint": "S2", "response": "Responder"},
    "HCC22-088-P3-S1_A": {"patient": "P3", "timepoint": "S1", "response": "Responder"},
    "HCC22-088-P3-S2": {"patient": "P3", "timepoint": "S2", "response": "Responder"},
    "HCC22-088-P4-S1": {"patient": "P4", "timepoint": "S1", "response": "Progressor"},
    "HCC22-088-P4-S2_1i_rep": {"patient": "P4", "timepoint": "S2", "response": "Progressor"},
    "HCC22-088-P5-S1": {"patient": "P5", "timepoint": "S1", "response": "Responder"},
    "HCC22-088-P5-S2_F_rep": {"patient": "P5", "timepoint": "S2", "response": "Responder"},
    "HCC22-088-P6-S1": {"patient": "P6", "timepoint": "S1", "response": "Responder"},
    "HCC22-088-P6-S2_D": {"patient": "P6", "timepoint": "S2", "response": "Responder"},
}


def load_module3_outputs(sample_name: str) -> tuple:
    """
    Load Module 3 outputs for a sample.

    Returns:
        Tuple of (adata, cell_type_proportions)
    """
    # Load raw spatial data
    sample_path = DATA_ROOT / sample_name / "outs"
    logger.info(f"Loading spatial data from {sample_path}")

    adata = sq.read.visium(
        str(sample_path),
        counts_file='filtered_feature_bc_matrix.h5',
        load_images=True,
        gex_only=True,  # Only GEX for Module 4
    )

    # Make var names unique
    adata.var_names_make_unique()

    # Load cell type proportions
    prop_file = MODULE3_OUTPUT / f"{sample_name}_cell_prop_finetuned_results.csv"
    if not prop_file.exists():
        prop_file = MODULE3_OUTPUT / f"{sample_name}_cell_prop_global_results.csv"

    logger.info(f"Loading proportions from {prop_file}")
    proportions = pd.read_csv(prop_file, index_col=0)

    # Align indices
    common_spots = adata.obs_names.intersection(proportions.index)
    logger.info(f"Common spots: {len(common_spots)}")

    adata = adata[common_spots].copy()
    proportions = proportions.loc[common_spots]

    # Create weighted deconvolved layers for joint discovery
    # These are needed by stack_deconvolved_layers()
    raw_X = adata.X.toarray() if scipy.sparse.issparse(adata.X) else np.array(adata.X)
    for ct in proportions.columns:
        layer_name = f"{ct}_genes_pass1"
        ct_weights = proportions[ct].values[:, np.newaxis]
        adata.layers[layer_name] = raw_X * ct_weights
    logger.info(f"Created {len(proportions.columns)} weighted layers for joint discovery")

    logger.info(f"Loaded adata: {adata.shape[0]} spots, {adata.shape[1]} genes")
    logger.info(f"Cell types: {list(proportions.columns)}")

    return adata, proportions


def run_module4_anchored(
    adata: sc.AnnData,
    proportions: pd.DataFrame,
    K_programs: int = 5,
    top_n_genes: int = 50,
    seed: int = 42,
) -> Dict:
    """Run anchored program discovery (per cell type)."""
    logger.info("Running ANCHORED program discovery...")

    result = discover_anchored_programs(
        adata=adata,
        cell_type_proportions=proportions,
        K_programs=K_programs,
        top_n_genes=top_n_genes,
        random_state=seed,
        validate_with_proteins=False,  # No protein data for validation
    )

    # Convert to serializable format
    output = {
        "mode": "anchored",
        "n_cell_types": result.n_anchors,
        "total_programs": result.total_programs,
        "anchors": {},
    }

    for anchor_name, anchor in result.results_by_anchor.items():
        anchor_data = {
            "cell_type": anchor_name,
            "n_programs": len(anchor.programs),
            "n_spots_used": int(anchor.n_spots_used),
            "programs": [],
        }

        for prog in anchor.programs:
            anchor_data["programs"].append({
                "program_id": int(prog.program_id),
                "top_genes": prog.top_genes[:20],
                "variance_explained": float(prog.variance_explained),
                "spatial_moran_i": float(prog.spatial_moran_i),
                "spatial_moran_pvalue": float(prog.spatial_moran_pvalue),
                "mean_activity": float(prog.mean_activity),
            })

        output["anchors"][anchor_name] = anchor_data

    logger.info(f"Anchored discovery complete: {result.n_anchors} cell types, {result.total_programs} programs")
    return output, result


def run_module4_joint(
    adata: sc.AnnData,
    proportions: pd.DataFrame,
    K_programs: int = 10,
    top_n_genes: int = 50,
    seed: int = 42,
) -> Dict:
    """Run joint program discovery (all cell types together)."""
    logger.info("Running JOINT program discovery...")

    result = discover_joint_programs(
        adata=adata,
        cell_type_proportions=proportions,
        K_programs=K_programs,
        layer_pattern="_genes_pass1",
        top_n_genes=top_n_genes,
        random_state=seed,
    )

    # Convert to serializable format
    output = {
        "mode": "joint",
        "n_programs": len(result.programs),
        "n_spots": result.n_spots,
        "reconstruction_error": float(result.reconstruction_error),
        "programs": [],
    }

    for prog in result.programs:
        output["programs"].append({
            "program_id": int(prog.program_id),
            "top_genes": prog.top_genes[:20],
            "variance_explained": float(prog.variance_explained),
            "spatial_moran_i": float(prog.spatial_moran_i),
            "spatial_moran_pvalue": float(prog.spatial_moran_pvalue),
            "mean_activity": float(prog.mean_activity),
            "primary_cell_type": prog.primary_cell_type,
            "secondary_cell_type": prog.secondary_cell_type,
            "interaction_score": float(prog.interaction_score),
            "program_type": prog.program_type,
        })

    # Summary by program type
    type_counts = {}
    for prog in result.programs:
        ptype = prog.program_type
        type_counts[ptype] = type_counts.get(ptype, 0) + 1
    output["program_type_counts"] = type_counts

    logger.info(f"Joint discovery complete: {len(result.programs)} programs")
    logger.info(f"  Types: {type_counts}")

    return output, result


def main():
    parser = argparse.ArgumentParser(
        description="Run Module 4 program discovery on patient samples"
    )
    parser.add_argument(
        "--sample",
        type=str,
        required=True,
        help=f"Sample name (one of {ALL_SAMPLES})",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="output/module4_discovery",
        help="Output directory for results",
    )
    parser.add_argument(
        "--mode",
        type=str,
        choices=["anchored", "joint", "both"],
        default="both",
        help="Discovery mode: anchored, joint, or both",
    )
    parser.add_argument("--k-anchored", type=int, default=5, help="Programs per cell type (anchored)")
    parser.add_argument("--k-joint", type=int, default=10, help="Total programs (joint)")
    parser.add_argument("--top-n-genes", type=int, default=50, help="Top genes per program")
    parser.add_argument("--seed", type=int, default=42, help="Random seed")

    args = parser.parse_args()

    # Validate sample name
    if args.sample not in ALL_SAMPLES:
        logger.error(f"Unknown sample: {args.sample}")
        logger.error(f"Valid samples: {ALL_SAMPLES}")
        sys.exit(1)

    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load Module 3 outputs
    logger.info(f"Processing sample: {args.sample}")
    adata, proportions = load_module3_outputs(args.sample)

    results = {
        "sample_name": args.sample,
        "metadata": SAMPLE_METADATA.get(args.sample, {}),
        "timestamp": datetime.now().isoformat(),
        "parameters": {
            "mode": args.mode,
            "k_anchored": args.k_anchored,
            "k_joint": args.k_joint,
            "top_n_genes": args.top_n_genes,
            "seed": args.seed,
        },
        "n_spots": adata.shape[0],
        "n_genes": adata.shape[1],
        "cell_types": list(proportions.columns),
    }

    # Run discovery
    if args.mode in ["anchored", "both"]:
        anchored_output, anchored_result = run_module4_anchored(
            adata, proportions,
            K_programs=args.k_anchored,
            top_n_genes=args.top_n_genes,
            seed=args.seed,
        )
        results["anchored"] = anchored_output

        # Save H matrices for Module 5
        anchored_H = {ct: res.H for ct, res in anchored_result.results_by_anchor.items()}
        np.save(output_dir / f"{args.sample}_anchored_H.npy", anchored_H, allow_pickle=True)

    if args.mode in ["joint", "both"]:
        joint_output, joint_result = run_module4_joint(
            adata, proportions,
            K_programs=args.k_joint,
            top_n_genes=args.top_n_genes,
            seed=args.seed,
        )
        results["joint"] = joint_output

        # Save matrices for Module 5
        np.save(output_dir / f"{args.sample}_joint_W.npy", joint_result.W)
        np.save(output_dir / f"{args.sample}_joint_H.npy", joint_result.H)

    # Save results
    output_file = output_dir / f"{args.sample}_module4_discovery.json"
    with open(output_file, "w") as f:
        json.dump(results, f, indent=2)

    logger.info(f"Saved results to {output_file}")

    # Print summary
    print(f"\n{'='*60}")
    print(f"Module 4 Discovery Summary: {args.sample}")
    print(f"{'='*60}")
    print(f"Spots: {adata.shape[0]}")
    print(f"Genes: {adata.shape[1]}")
    print(f"Cell types: {len(proportions.columns)}")

    if "anchored" in results:
        print(f"\nAnchored Discovery:")
        for ct, data in results["anchored"]["anchors"].items():
            spatial_progs = sum(1 for p in data["programs"] if p["spatial_moran_i"] > 0.1)
            print(f"  {ct}: {data['n_programs']} programs ({spatial_progs} spatial)")

    if "joint" in results:
        print(f"\nJoint Discovery:")
        print(f"  Total programs: {results['joint']['n_programs']}")
        for ptype, count in results["joint"]["program_type_counts"].items():
            print(f"  {ptype}: {count}")

    print(f"\nOutput: {output_file}")


if __name__ == "__main__":
    main()
