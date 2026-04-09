#!/usr/bin/env python
"""
Module 4b Bivariate Program Relationship Analysis.

Computes spatial co-localization (bivariate Moran's I) between program pairs
within each sample. This enables identification of:
- Co-localized programs: Different cell type programs that peak together spatially
- Exclusive programs: Programs that avoid each other spatially

Usage:
    python run_module4b_bivariate.py --sample HCC22-088-P1-S1 --output-dir output/module4b_unified
"""
import argparse
import json
import logging
import sys
from datetime import datetime
from pathlib import Path
from typing import Dict, Optional

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse
import squidpy as sq

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from CITEgeist.model.programs.anchored_program_discovery import (
    AnchoredProgramDiscoveryResult,
    AnchoredProgramResult,
    BivariateProgramResult,
    SpatialProgram,
    analyze_program_relationships,
)

logger = logging.getLogger(__name__)
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)

# Data locations
DATA_ROOT = Path("/ix1/alee/LO_LAB/General/Lab_Data/20250210_CITEGeistPublicData_GEO_Alex/processed_files")
MODULE3_OUTPUT = Path("/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/output/module3_nb")
MODULE4_OUTPUT = Path("/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/output/module4_nb")

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


def load_module4_result(sample_name: str) -> Optional[AnchoredProgramDiscoveryResult]:
    """
    Load Module 4 results and reconstruct AnchoredProgramDiscoveryResult.

    Args:
        sample_name: Sample identifier

    Returns:
        AnchoredProgramDiscoveryResult or None if not found
    """
    json_file = MODULE4_OUTPUT / f"{sample_name}_module4_discovery.json"
    h_file = MODULE4_OUTPUT / f"{sample_name}_anchored_H.npy"

    if not json_file.exists():
        logger.warning(f"Module 4 JSON not found: {json_file}")
        return None

    with open(json_file) as f:
        data = json.load(f)

    if "anchored" not in data:
        logger.warning(f"No anchored results in {json_file}")
        return None

    # Load H matrices
    H_matrices = {}
    if h_file.exists():
        H_matrices = np.load(h_file, allow_pickle=True).item()
        logger.info(f"Loaded H matrices for {len(H_matrices)} cell types")

    # Convert JSON to AnchoredProgramResult objects
    results_by_anchor = {}

    for anchor_name, anchor_data in data["anchored"]["anchors"].items():
        programs = []

        for prog_data in anchor_data["programs"]:
            top_genes = prog_data["top_genes"]
            gene_loadings = {g: 1.0 / (i + 1) for i, g in enumerate(top_genes)}

            programs.append(SpatialProgram(
                program_id=prog_data["program_id"],
                top_genes=top_genes,
                gene_loadings=gene_loadings,
                variance_explained=prog_data["variance_explained"],
                spatial_moran_i=prog_data["spatial_moran_i"],
                spatial_moran_pvalue=prog_data["spatial_moran_pvalue"],
                mean_activity=prog_data["mean_activity"],
                active_spots_fraction=0.0,
            ))

        # Create W matrix from gene loadings
        all_genes = list(set(g for p in programs for g in p.top_genes))
        n_progs = len(programs)
        W = np.zeros((len(all_genes), n_progs))

        for k, prog in enumerate(programs):
            for g, loading in prog.gene_loadings.items():
                if g in all_genes:
                    W[all_genes.index(g), k] = loading

        # Get H matrix from saved file
        H = H_matrices.get(anchor_name, np.zeros((n_progs, 1)))

        results_by_anchor[anchor_name] = AnchoredProgramResult(
            anchor_name=anchor_name,
            anchor_proteins=[],
            programs=programs,
            W=W,
            H=H,
            gene_names=all_genes,
            protein_correlations=pd.DataFrame(),
            reconstruction_error=0.0,
            n_spots_used=anchor_data["n_spots_used"],
            parameters={"sample_name": sample_name},
        )

    return AnchoredProgramDiscoveryResult(
        results_by_anchor=results_by_anchor,
        n_anchors=len(results_by_anchor),
        total_programs=sum(len(r.programs) for r in results_by_anchor.values()),
        profile_discovery_result=None,
        parameters={"source_file": str(json_file)},
    )


def load_spatial_adata(sample_name: str) -> sc.AnnData:
    """
    Load spatial AnnData for a sample.

    Args:
        sample_name: Sample identifier

    Returns:
        AnnData with spatial coordinates
    """
    # Try loading from Module 3 output first (has deconvolved layers)
    module3_file = MODULE3_OUTPUT / f"{sample_name}_module3_results.h5ad"
    if module3_file.exists():
        logger.info(f"Loading from Module 3 output: {module3_file}")
        adata = sc.read_h5ad(module3_file)
        return adata

    # Fall back to raw data
    sample_path = DATA_ROOT / sample_name / "outs"
    logger.info(f"Loading from raw data: {sample_path}")

    adata = sq.read.visium(
        str(sample_path),
        counts_file='filtered_feature_bc_matrix.h5',
        load_images=True,
        gex_only=True,
    )

    # Filter NaN coordinates
    if 'spatial' in adata.obsm:
        coords = adata.obsm['spatial']
        valid_mask = np.all(np.isfinite(coords), axis=1)
        if (~valid_mask).sum() > 0:
            adata = adata[valid_mask].copy()

    return adata


def run_module4b(
    sample_name: str,
    output_dir: Path,
    n_permutations: int = 199,
    fdr_threshold: float = 0.05,
    min_bivariate_i: float = 0.1,
    neighbor_k: int = 8,
    seed: int = 42,
) -> Optional[BivariateProgramResult]:
    """
    Run Module 4b bivariate analysis for a sample.

    Args:
        sample_name: Sample identifier
        output_dir: Output directory
        n_permutations: Permutations for p-value
        fdr_threshold: FDR threshold
        min_bivariate_i: Minimum |I| for strong relationship
        neighbor_k: k-NN for spatial weights
        seed: Random seed

    Returns:
        BivariateProgramResult or None if failed
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load Module 4 results
    logger.info(f"Loading Module 4 results for {sample_name}...")
    module4_result = load_module4_result(sample_name)

    if module4_result is None:
        logger.error(f"Failed to load Module 4 results for {sample_name}")
        return None

    logger.info(f"Loaded {module4_result.n_anchors} anchors, {module4_result.total_programs} programs")

    # Load spatial data
    logger.info("Loading spatial data...")
    adata = load_spatial_adata(sample_name)
    logger.info(f"Loaded adata: {adata.shape[0]} spots")

    # Run bivariate analysis
    logger.info("Running bivariate analysis...")
    result = analyze_program_relationships(
        result=module4_result,
        adata=adata,
        fdr_threshold=fdr_threshold,
        min_bivariate_i=min_bivariate_i,
        n_permutations=n_permutations,
        neighbor_k=neighbor_k,
        include_within_anchor=True,
        random_state=seed,
    )

    # Save results
    output_file = output_dir / f"{sample_name}_module4b_relationships.csv"
    result.to_dataframe().to_csv(output_file, index=False)
    logger.info(f"Saved results to {output_file}")

    # Save summary JSON
    summary = {
        "sample_name": sample_name,
        "timestamp": datetime.now().isoformat(),
        "n_programs_total": result.n_programs_total,
        "n_pairs_tested": result.n_pairs_tested,
        "n_significant": result.n_significant,
        "fdr_threshold": fdr_threshold,
        "parameters": {
            "n_permutations": n_permutations,
            "min_bivariate_i": min_bivariate_i,
            "neighbor_k": neighbor_k,
        },
        "relationship_summary": {
            "co-localized": sum(1 for p in result.significant_pairs if p.relationship_type == "co-localized"),
            "exclusive": sum(1 for p in result.significant_pairs if p.relationship_type == "exclusive"),
            "other": sum(1 for p in result.significant_pairs if p.relationship_type not in ["co-localized", "exclusive"]),
        },
    }

    summary_file = output_dir / f"{sample_name}_module4b_summary.json"
    with open(summary_file, "w") as f:
        json.dump(summary, f, indent=2)

    return result


def main():
    parser = argparse.ArgumentParser(
        description="Run Module 4b bivariate program analysis"
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
        default="output/module4b_unified",
        help="Output directory for results",
    )
    parser.add_argument("--n-permutations", type=int, default=199, help="Permutations for p-value")
    parser.add_argument("--fdr-threshold", type=float, default=0.05, help="FDR threshold")
    parser.add_argument("--min-bivariate-i", type=float, default=0.1, help="Min |I| for strong")
    parser.add_argument("--neighbor-k", type=int, default=8, help="k-NN for spatial weights")
    parser.add_argument("--seed", type=int, default=42, help="Random seed")

    args = parser.parse_args()

    # Validate sample name
    if args.sample not in ALL_SAMPLES:
        logger.error(f"Unknown sample: {args.sample}")
        logger.error(f"Valid samples: {ALL_SAMPLES}")
        sys.exit(1)

    # Run Module 4b
    result = run_module4b(
        sample_name=args.sample,
        output_dir=Path(args.output_dir),
        n_permutations=args.n_permutations,
        fdr_threshold=args.fdr_threshold,
        min_bivariate_i=args.min_bivariate_i,
        neighbor_k=args.neighbor_k,
        seed=args.seed,
    )

    if result is None:
        print(f"Module 4b failed for {args.sample}")
        sys.exit(1)

    # Print summary
    print(f"\n{'='*60}")
    print(f"Module 4b Complete: {args.sample}")
    print(f"{'='*60}")
    print(f"Programs analyzed: {result.n_programs_total}")
    print(f"Pairs tested: {result.n_pairs_tested}")
    print(f"Significant pairs: {result.n_significant}")

    if result.significant_pairs:
        co_loc = sum(1 for p in result.significant_pairs if p.relationship_type == "co-localized")
        excl = sum(1 for p in result.significant_pairs if p.relationship_type == "exclusive")
        print(f"  Co-localized: {co_loc}")
        print(f"  Exclusive: {excl}")

        print(f"\nTop 5 relationships:")
        sorted_pairs = sorted(result.significant_pairs, key=lambda p: abs(p.bivariate_morans_i), reverse=True)
        for pair in sorted_pairs[:5]:
            print(f"  {pair.anchor1}_{pair.program1_id} <-> {pair.anchor2}_{pair.program2_id}: "
                  f"I={pair.bivariate_morans_i:.3f} ({pair.relationship_type})")

    print(f"\nOutput: {args.output_dir}")


if __name__ == "__main__":
    main()
