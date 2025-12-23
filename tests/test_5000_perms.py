"""
Test with 5000 permutations to see timing and FDR behavior.
"""

import numpy as np
import pandas as pd
import scanpy as sc
import logging
import sys
import time

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(message)s')
logger = logging.getLogger(__name__)

sys.path.insert(0, '/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist')

from CITEgeist.model import (
    identify_interesting_markers,
    analyze_marker_colocalization,
    discover_profiles,
    select_profiles_by_reconstruction,
)

SIMULATED_DATA = "/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/replicates/high_seg/h5ad_objects/Wu_rep_1_CITE.h5ad"


def main():
    logger.info("=" * 70)
    logger.info("TEST: 5000 PERMUTATIONS ON SIMULATED DATA")
    logger.info("=" * 70)

    # Load simulated data
    adata = sc.read_h5ad(SIMULATED_DATA)
    X = adata.X if not hasattr(adata.X, 'toarray') else adata.X.toarray()
    coords = adata.obsm['spatial']
    marker_names = list(adata.var_names)

    logger.info(f"Loaded: {X.shape[0]} spots, {len(marker_names)} markers")

    # Ground truth
    gt_cell_types = ['B-cells', 'CAFs', 'Cancer Epithelial', 'Endothelial',
                     'Myeloid', 'Normal Epithelial', 'PVL', 'Plasmablasts', 'T-cells']
    logger.info(f"Ground truth: {len(gt_cell_types)} cell types")

    # Module 1: Identify interesting markers
    logger.info("\nRunning Module 1 (marker interest)...")
    result1 = identify_interesting_markers(X, coords, marker_names, verbose=False)
    logger.info(f"Interesting markers: {len(result1.interesting_markers)}")

    # Module 2a: Colocalization analysis with 5000 permutations
    logger.info("\n" + "=" * 50)
    logger.info("Running Module 2a with 5000 permutations...")
    logger.info("=" * 50)

    start_time = time.time()
    coloc_result = analyze_marker_colocalization(
        X, coords, marker_names,
        markers_to_analyze=result1.interesting_markers,
        n_permutations=5000,
        seed=42,
        verbose=True
    )
    elapsed = time.time() - start_time

    logger.info(f"\n5000 permutations completed in {elapsed:.1f} seconds ({elapsed/60:.1f} minutes)")
    logger.info(f"Pairs analyzed: {len(coloc_result.pairs)}")

    # Check p-value distribution
    p_values = [p.neighbor_enrichment_pvalue for p in coloc_result.pairs]
    logger.info(f"\nP-value distribution:")
    logger.info(f"  Min p-value: {min(p_values):.6f}")
    logger.info(f"  Expected min with 5000 perms: {1/5001:.6f}")
    logger.info(f"  Pairs with p < 0.01: {sum(p < 0.01 for p in p_values)}")
    logger.info(f"  Pairs with p < 0.001: {sum(p < 0.001 for p in p_values)}")
    logger.info(f"  Pairs with p < 0.0005: {sum(p < 0.0005 for p in p_values)}")

    # Module 2b: Profile discovery
    logger.info("\nRunning Module 2b (profile discovery, k=3)...")
    profile_result = discover_profiles(
        coloc_result,
        fdr_alpha=0.05,
        top_k=3,
        verbose=True
    )
    logger.info(f"Discovered {len(profile_result.profiles)} profiles")

    # Check ground truth recovery
    gt_profiles = {ct: set([f'{ct}_Protein_1', f'{ct}_Protein_2']) for ct in gt_cell_types}
    recovered = 0
    logger.info("\nProfile recovery:")
    for i, profile in enumerate(profile_result.profiles):
        ps = set(profile)
        match = None
        for ct, gt_set in gt_profiles.items():
            if ps == gt_set:
                match = ct
                recovered += 1
                break
        status = f"-> {match}" if match else ""
        logger.info(f"  {i+1}. {profile} {status}")

    logger.info(f"\nGround truth recovery: {recovered}/{len(gt_cell_types)} exact matches")

    # Module 2c: Profile selection
    if len(profile_result.profiles) > 1:
        logger.info("\nRunning Module 2c (reconstruction-based selection)...")
        selection = select_profiles_by_reconstruction(
            X, marker_names, profile_result.profiles,
            colocalization_result=coloc_result,
            interesting_markers=result1.interesting_markers,  # FIXED: pass interesting markers
            verbose=True
        )
        logger.info(f"\nOptimal profiles: {selection.optimal_n}/{len(selection.all_profiles_ranked)}")
        logger.info(f"Variance explained at elbow: {selection.variance_explained[selection.optimal_n - 1]:.1%}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
