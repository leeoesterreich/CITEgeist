"""
Full CITEgeist pipeline on Xenium single-cell data.

This script demonstrates CITEgeist as a resolution-independent spatial
multi-omics integration toolbox by running Modules 1-4 on single-cell data.

Pipeline:
1. Module 1: Identify spatially-interesting protein markers
2. Module 2: Discover spatial proteomic profiles via colocalization
3. Cell Assignment: Assign cells to profiles (replaces deconvolution)
4. Module 4: Discover protein-anchored transcriptional programs

This demonstrates that:
- CITEgeist works on both spot-level and single-cell data
- Spatial proteomic patterns are discovered (not just cell types)
- Protein-anchored programs reveal heterogeneity within populations
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
from sklearn.decomposition import NMF
from scipy.spatial import cKDTree
from scipy.stats import pearsonr

# Add CITEgeist to path
REPO_ROOT = Path(__file__).parent.parent.parent.parent.parent
sys.path.insert(0, str(REPO_ROOT / "CITEgeist"))
sys.path.insert(0, str(REPO_ROOT / "Benchmarking" / "xenium_pseudovisium" / "src"))

from model.marker_interest import identify_interesting_markers, MarkerInterestResult
from model.spatial_colocalization import (
    analyze_marker_colocalization,
    discover_profiles,
    ColocalizationResult,
    ProfileDiscoveryResult,
)
from model.cell_assignment import assign_cells_to_profiles, CellAssignmentResult
from load_xenium_singlecell import load_xenium_singlecell, XENIUM_DATA_DIR

logger = logging.getLogger(__name__)

# Output directory
OUTPUT_DIR = REPO_ROOT / "Benchmarking" / "xenium_benchmarking" / "CITEgeist" / "output_singlecell"


def run_module1(
    adata_protein: sc.AnnData,
    morans_k: int = 20,
    smooth_k: int = 10,
    morans_n_perm: int = 99,
    seed: int = 42,
) -> MarkerInterestResult:
    """Run Module 1: Marker Interest Detection."""
    logger.info("=" * 60)
    logger.info("MODULE 1: Marker Interest Detection")
    logger.info("=" * 60)

    X = adata_protein.X
    if hasattr(X, "toarray"):
        X = X.toarray()
    X = np.asarray(X, dtype=np.float64)

    coords = adata_protein.obsm["spatial"]
    marker_names = list(adata_protein.var_names)

    result = identify_interesting_markers(
        X=X,
        coords=coords,
        marker_names=marker_names,
        morans_k=morans_k,
        smooth_k=smooth_k,
        morans_n_perm=morans_n_perm,
        seed=seed,
        verbose=True,
    )

    logger.info(f"Identified {len(result.interesting_markers)} interesting markers")
    logger.info(f"Top 5: {result.interesting_markers[:5]}")

    return result


def run_module2(
    adata_protein: sc.AnnData,
    interesting_markers: list,
    neighbor_k: int = 15,
    fdr_alpha: float = 0.05,
    seed: int = 42,
) -> ProfileDiscoveryResult:
    """Run Module 2: Spatial Colocalization and Profile Discovery."""
    logger.info("\n" + "=" * 60)
    logger.info("MODULE 2: Spatial Colocalization & Profile Discovery")
    logger.info("=" * 60)

    X = adata_protein.X
    if hasattr(X, "toarray"):
        X = X.toarray()
    X = np.asarray(X, dtype=np.float64)

    coords = adata_protein.obsm["spatial"]
    marker_names = list(adata_protein.var_names)

    # Module 2a: Colocalization analysis
    logger.info("Running colocalization analysis...")
    coloc_result = analyze_marker_colocalization(
        X=X,
        coords=coords,
        marker_names=marker_names,
        markers_to_analyze=interesting_markers,
        neighbor_k=neighbor_k,
        multi_scale_k=[neighbor_k // 2, neighbor_k, neighbor_k * 2],
        n_permutations=99,
        seed=seed,
    )
    logger.info(f"Analyzed {len(coloc_result.pairs)} marker pairs")

    # Module 2b: Profile discovery
    logger.info("Discovering profiles via hierarchical clustering...")
    profile_result = discover_profiles(
        colocalization_result=coloc_result,
        fdr_alpha=fdr_alpha,
        top_k=3,
        seed=seed,
    )
    logger.info(f"Discovered {len(profile_result.profiles)} profiles")

    for i, profile in enumerate(profile_result.profiles):
        logger.info(f"  Profile {i}: {profile}")

    return profile_result


def run_cell_assignment(
    adata_protein: sc.AnnData,
    profiles: list,
    method: str = "weighted_sum",
    spatial_smoothing: float = 0.1,
    k_neighbors: int = 15,
) -> CellAssignmentResult:
    """Run cell assignment (Module 3 adaptation for single cells)."""
    logger.info("\n" + "=" * 60)
    logger.info("CELL ASSIGNMENT (Module 3 Adaptation)")
    logger.info("=" * 60)

    X = adata_protein.X
    if hasattr(X, "toarray"):
        X = X.toarray()
    X = np.asarray(X, dtype=np.float64)

    coords = adata_protein.obsm["spatial"]
    marker_names = list(adata_protein.var_names)
    cell_ids = list(adata_protein.obs_names)

    # Convert profile list to dict
    profile_dict = {f"Profile_{i}": profile for i, profile in enumerate(profiles)}

    result = assign_cells_to_profiles(
        protein_data=X,
        coords=coords,
        profile_dict=profile_dict,
        marker_names=marker_names,
        cell_ids=cell_ids,
        method=method,
        spatial_smoothing=spatial_smoothing,
        k_neighbors=k_neighbors,
    )

    return result


def run_module4_singlecell(
    adata_gex: sc.AnnData,
    assignment_result: CellAssignmentResult,
    K_programs: int = 5,
    min_cells: int = 100,
    top_genes: int = 50,
    seed: int = 42,
) -> Dict:
    """
    Run Module 4: Protein-Anchored Program Discovery on single cells.

    For single-cell data, we have actual cell-type-specific expression
    (from cells assigned to each profile) rather than deconvolved layers.

    Args:
        adata_gex: Gene expression AnnData.
        assignment_result: Cell assignment results.
        K_programs: Number of programs per profile.
        min_cells: Minimum cells to analyze a profile.
        top_genes: Number of top genes to report per program.
        seed: Random seed.

    Returns:
        Dict with program discovery results per profile.
    """
    logger.info("\n" + "=" * 60)
    logger.info("MODULE 4: Protein-Anchored Program Discovery (Single-Cell)")
    logger.info("=" * 60)

    X_gex = adata_gex.X
    if hasattr(X_gex, "toarray"):
        X_gex = X_gex.toarray()
    X_gex = np.asarray(X_gex, dtype=np.float64)

    coords = adata_gex.obsm["spatial"]
    gene_names = list(adata_gex.var_names)

    results = {}

    for profile_idx, profile_name in enumerate(assignment_result.profile_names):
        if profile_name == "Unknown":
            continue

        # Get cells assigned to this profile (high confidence)
        cell_mask = (
            (assignment_result.hard_assignments == profile_idx) &
            (assignment_result.assignment_confidence >= 0.5)
        )
        n_cells = cell_mask.sum()

        if n_cells < min_cells:
            logger.info(f"Skipping {profile_name}: only {n_cells} cells (need {min_cells})")
            continue

        logger.info(f"\n{profile_name}: {n_cells} cells")

        # Extract GEX for these cells
        X_profile = X_gex[cell_mask, :]
        coords_profile = coords[cell_mask, :]

        # Filter low-expressed genes
        gene_expr = X_profile.sum(axis=0)
        active_genes = gene_expr > np.percentile(gene_expr, 50)
        X_profile_filtered = X_profile[:, active_genes]
        gene_names_filtered = [gene_names[i] for i in np.where(active_genes)[0]]

        # Normalize (log1p then scale)
        X_norm = np.log1p(X_profile_filtered)
        X_norm = X_norm / (X_norm.max(axis=0, keepdims=True) + 1e-10)

        # Run NMF
        logger.info(f"  Running NMF with K={K_programs} programs...")
        nmf = NMF(
            n_components=K_programs,
            init="nndsvda",
            random_state=seed,
            max_iter=500,
        )
        W = nmf.fit_transform(X_norm)  # (n_cells, K)
        H = nmf.components_  # (K, n_genes)

        # Compute spatial Moran's I for each program
        logger.info("  Computing spatial coherence (Moran's I)...")
        morans_i_values = []
        for k in range(K_programs):
            program_activity = W[:, k]
            mi = _compute_morans_i_simple(program_activity, coords_profile)
            morans_i_values.append(mi)

        # Extract top genes per program
        programs = []
        for k in range(K_programs):
            loadings = H[k, :]
            top_idx = np.argsort(loadings)[::-1][:top_genes]
            top_gene_list = [(gene_names_filtered[i], loadings[i]) for i in top_idx]

            programs.append({
                "program_id": k,
                "top_genes": [g[0] for g in top_gene_list[:10]],
                "top_genes_loadings": top_gene_list,
                "moran_i": morans_i_values[k],
                "mean_activity": W[:, k].mean(),
                "variance_explained": (W[:, k] ** 2).sum() / (W ** 2).sum(),
            })

            logger.info(f"  Program {k}: Moran's I={morans_i_values[k]:.3f}, "
                       f"top genes: {', '.join([g[0] for g in top_gene_list[:5]])}")

        results[profile_name] = {
            "n_cells": n_cells,
            "W": W,
            "H": H,
            "gene_names": gene_names_filtered,
            "programs": programs,
            "reconstruction_error": nmf.reconstruction_err_,
        }

    return results


def _compute_morans_i_simple(
    values: np.ndarray,
    coords: np.ndarray,
    k: int = 10,
) -> float:
    """Simple Moran's I computation."""
    n = len(values)
    if n < 10:
        return np.nan

    tree = cKDTree(coords)
    _, indices = tree.query(coords, k=min(k + 1, n))

    if indices.ndim == 1:
        return np.nan

    indices = indices[:, 1:]  # Remove self
    z = values - values.mean()
    denom = np.sum(z ** 2)

    if denom == 0:
        return np.nan

    num = 0.0
    S0 = 0.0
    for i in range(n):
        for j in indices[i]:
            num += z[i] * z[j]
            S0 += 1

    if S0 == 0:
        return np.nan

    I = (n / S0) * (num / denom)
    return float(I)


def main():
    parser = argparse.ArgumentParser(
        description="Run full CITEgeist pipeline on Xenium single-cell data"
    )
    parser.add_argument(
        "--region", type=int, default=0, help="Region ID (0-4)"
    )
    parser.add_argument(
        "--max-cells", type=int, default=None,
        help="Max cells to use (None for all)"
    )
    parser.add_argument(
        "--k-programs", type=int, default=5,
        help="Number of programs per profile (default: 5)"
    )
    parser.add_argument(
        "--seed", type=int, default=42,
        help="Random seed"
    )
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    )

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Load single-cell data
    logger.info("Loading Xenium single-cell data...")
    adata_gex, adata_protein = load_xenium_singlecell(
        region_id=args.region,
        max_cells=args.max_cells,
        seed=args.seed,
    )

    n_cells = adata_protein.shape[0]
    logger.info(f"Loaded {n_cells:,} cells")

    # Run pipeline
    start_time = datetime.now()

    # Module 1
    module1_result = run_module1(adata_protein, seed=args.seed)

    # Module 2
    module2_result = run_module2(
        adata_protein,
        interesting_markers=module1_result.interesting_markers,
        seed=args.seed,
    )

    # Cell Assignment
    assignment_result = run_cell_assignment(
        adata_protein,
        profiles=module2_result.profiles,
    )

    # Module 4
    module4_result = run_module4_singlecell(
        adata_gex,
        assignment_result,
        K_programs=args.k_programs,
        seed=args.seed,
    )

    elapsed = datetime.now() - start_time

    # Save results
    output_prefix = f"region_{args.region}"
    if args.max_cells:
        output_prefix += f"_maxcells_{args.max_cells}"

    # Module 1 results
    module1_result.to_dataframe().to_csv(
        OUTPUT_DIR / f"{output_prefix}_module1.csv", index=False
    )

    # Assignment results
    assignment_result.to_dataframe().to_csv(
        OUTPUT_DIR / f"{output_prefix}_cell_assignments.csv"
    )

    # Module 4 summary (convert numpy types to native Python for JSON)
    module4_summary = {
        "region_id": int(args.region),
        "n_cells": int(n_cells),
        "n_profiles_analyzed": int(len(module4_result)),
        "k_programs": int(args.k_programs),
        "profiles": {},
    }
    for profile_name, result in module4_result.items():
        module4_summary["profiles"][profile_name] = {
            "n_cells": int(result["n_cells"]),
            "programs": [
                {
                    "program_id": int(p["program_id"]),
                    "moran_i": float(p["moran_i"]) if not np.isnan(p["moran_i"]) else None,
                    "top_genes": p["top_genes"],
                    "variance_explained": float(p["variance_explained"]),
                }
                for p in result["programs"]
            ],
        }

    with open(OUTPUT_DIR / f"{output_prefix}_module4_summary.json", "w") as f:
        json.dump(module4_summary, f, indent=2)

    # Print summary
    print("\n" + "=" * 60)
    print("PIPELINE COMPLETE")
    print("=" * 60)
    print(f"Region: {args.region}")
    print(f"Cells analyzed: {n_cells:,}")
    print(f"Interesting markers: {len(module1_result.interesting_markers)}")
    print(f"Profiles discovered: {len(module2_result.profiles)}")
    print(f"Profiles with programs: {len(module4_result)}")
    print(f"Total runtime: {elapsed}")
    print(f"\nResults saved to: {OUTPUT_DIR}")

    # Detailed program summary
    print("\n" + "-" * 60)
    print("PROGRAM DISCOVERY SUMMARY")
    print("-" * 60)
    for profile_name, result in module4_result.items():
        print(f"\n{profile_name} ({result['n_cells']} cells):")
        for prog in result["programs"]:
            print(f"  Program {prog['program_id']}: "
                  f"Moran's I={prog['moran_i']:.3f}, "
                  f"top genes: {', '.join(prog['top_genes'][:3])}")


if __name__ == "__main__":
    main()
