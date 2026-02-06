#!/usr/bin/env python
"""
Run Module 4 on Xenium single-cell data using discovered profiles.

At single-cell resolution, assigns cells to profiles (no deconvolution needed),
then runs NMF to discover spatial transcriptomic programs within each cell type.

Usage:
    python run_singlecell_module4.py --mode full
    python run_singlecell_module4.py --mode quadrant --quadrant-id 0
"""

import argparse
import json
import logging
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.spatial import cKDTree
from sklearn.decomposition import NMF
from sklearn.preprocessing import normalize

# Add paths
REPO_ROOT = Path(__file__).parent.parent.parent.parent.parent
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(REPO_ROOT / "Benchmarking" / "xenium_benchmarking" / "CITEgeist" / "src"))

from load_xenium_singlecell import load_xenium_singlecell

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)

# Output directory
OUTPUT_BASE = REPO_ROOT / "Benchmarking" / "xenium_benchmarking" / "CITEgeist" / "output_singlecell_demonstration"


def get_quadrant_bounds(coords: np.ndarray, quadrant_id: int) -> Tuple[float, float, float, float]:
    """Get bounds for a spatial quadrant."""
    x_mid = (coords[:, 0].min() + coords[:, 0].max()) / 2
    y_mid = (coords[:, 1].min() + coords[:, 1].max()) / 2
    x_min, x_max = coords[:, 0].min(), coords[:, 0].max()
    y_min, y_max = coords[:, 1].min(), coords[:, 1].max()

    bounds = {
        0: (x_min, x_mid, y_min, y_mid),
        1: (x_mid, x_max, y_min, y_mid),
        2: (x_min, x_mid, y_mid, y_max),
        3: (x_mid, x_max, y_mid, y_max),
    }
    return bounds[quadrant_id]


def assign_cells_to_profiles(
    X_protein: np.ndarray,
    marker_names: List[str],
    profiles: Dict[str, Dict],
    method: str = "soft",
) -> Tuple[np.ndarray, np.ndarray, List[str]]:
    """
    Assign cells to discovered profiles based on protein expression.

    Args:
        X_protein: Protein expression matrix (n_cells x n_proteins)
        marker_names: Protein names
        profiles: Dict of profile_name -> {"markers": [list of markers]}
        method: "soft" for probability assignment, "hard" for argmax

    Returns:
        assignments: (n_cells,) array of profile indices (hard) or
                     (n_cells, n_profiles) probability matrix (soft)
        profile_names: List of profile names in order
    """
    profile_names = list(profiles.keys())
    n_cells = X_protein.shape[0]
    n_profiles = len(profile_names)

    # Compute profile scores for each cell
    scores = np.zeros((n_cells, n_profiles))

    for i, (name, profile) in enumerate(profiles.items()):
        markers = profile["markers"]
        marker_indices = [marker_names.index(m) for m in markers if m in marker_names]

        if marker_indices:
            # Average normalized expression of profile markers
            profile_expr = X_protein[:, marker_indices]
            # Normalize per marker, then average
            profile_expr_norm = profile_expr / (profile_expr.max(axis=0, keepdims=True) + 1e-10)
            scores[:, i] = profile_expr_norm.mean(axis=1)

    if method == "soft":
        # Softmax to get probabilities
        scores_exp = np.exp(scores - scores.max(axis=1, keepdims=True))
        assignments = scores_exp / scores_exp.sum(axis=1, keepdims=True)
    else:
        # Hard assignment
        assignments = scores.argmax(axis=1)

    return assignments, scores, profile_names


def compute_morans_i(
    values: np.ndarray,
    coords: np.ndarray,
    k: int = 15,
) -> float:
    """Compute Moran's I for spatial autocorrelation."""
    n = len(values)
    if n < k + 1:
        return 0.0

    # Build KNN graph
    tree = cKDTree(coords)
    _, indices = tree.query(coords, k=k + 1)

    # Compute Moran's I
    z = values - values.mean()
    z_sq_sum = (z ** 2).sum()

    if z_sq_sum == 0:
        return 0.0

    # Sum of cross-products for neighbors
    cross_sum = 0.0
    for i in range(n):
        for j in indices[i, 1:]:  # Skip self
            cross_sum += z[i] * z[j]

    # Total weight (k neighbors per cell)
    W = n * k

    morans_i = (n / W) * (cross_sum / z_sq_sum)
    return float(morans_i)


def run_nmf_for_celltype(
    adata_gex: sc.AnnData,
    cell_mask: np.ndarray,
    n_programs: int = 5,
    min_cells: int = 500,
) -> Optional[Dict]:
    """
    Run NMF on cells of a specific type to discover programs.

    Returns:
        Dict with W (gene loadings), H (cell scores), top genes per program, Moran's I
    """
    n_cells = cell_mask.sum()
    if n_cells < min_cells:
        logger.warning(f"Only {n_cells} cells, skipping (min: {min_cells})")
        return None

    # Extract expression for this cell type
    X = adata_gex[cell_mask, :].X
    if hasattr(X, "todense"):
        X = np.asarray(X.todense())

    # Filter low-variance genes
    gene_var = X.var(axis=0)
    top_genes_mask = gene_var > np.percentile(gene_var, 50)
    X_filtered = X[:, top_genes_mask]
    gene_names_filtered = adata_gex.var_names[top_genes_mask].tolist()

    # Normalize
    X_norm = X_filtered / (X_filtered.sum(axis=1, keepdims=True) + 1e-10) * 1e4
    X_norm = np.log1p(X_norm)

    # Make non-negative
    X_norm = X_norm - X_norm.min()

    # Run NMF
    logger.info(f"Running NMF with K={n_programs} on {n_cells} cells x {X_norm.shape[1]} genes")
    model = NMF(n_components=n_programs, init="nndsvd", random_state=42, max_iter=500)
    H = model.fit_transform(X_norm)  # (n_cells, n_programs)
    W = model.components_.T  # (n_genes, n_programs)

    # Get coordinates for this cell type
    coords = adata_gex[cell_mask, :].obsm["spatial"]

    # Compute top genes and Moran's I per program
    programs = []
    for prog_idx in range(n_programs):
        # Top genes
        gene_loadings = W[:, prog_idx]
        top_indices = np.argsort(gene_loadings)[::-1][:20]
        top_genes = [gene_names_filtered[i] for i in top_indices]
        top_loadings = gene_loadings[top_indices].tolist()

        # Moran's I for spatial coherence
        cell_scores = H[:, prog_idx]
        morans_i = compute_morans_i(cell_scores, coords, k=15)

        programs.append({
            "program_id": prog_idx,
            "top_genes": top_genes,
            "top_loadings": top_loadings,
            "morans_i": morans_i,
            "variance_explained": float((cell_scores ** 2).sum() / (H ** 2).sum()),
        })

    return {
        "n_cells": int(n_cells),
        "n_genes": int(X_norm.shape[1]),
        "n_programs": n_programs,
        "programs": programs,
        "W_shape": list(W.shape),
        "H_shape": list(H.shape),
        "reconstruction_error": float(model.reconstruction_err_),
    }


def main():
    parser = argparse.ArgumentParser(description="Run Module 4 on single-cell data")
    parser.add_argument(
        "--mode",
        choices=["full", "quadrant"],
        required=True,
        help="Run on full dataset or single quadrant",
    )
    parser.add_argument(
        "--quadrant-id",
        type=int,
        choices=[0, 1, 2, 3],
        help="Quadrant ID (required if mode=quadrant)",
    )
    parser.add_argument(
        "--n-programs",
        type=int,
        default=5,
        help="Number of NMF programs per cell type (default: 5)",
    )
    args = parser.parse_args()

    if args.mode == "quadrant" and args.quadrant_id is None:
        parser.error("--quadrant-id required when mode=quadrant")

    # Determine paths
    if args.mode == "full":
        output_subdir = "full"
    else:
        output_subdir = f"quadrants/Q{args.quadrant_id}"

    output_dir = OUTPUT_BASE / output_subdir

    # Load Module 2 profiles
    profiles_path = output_dir / "module2c_profiles_selected.json"
    if not profiles_path.exists():
        logger.error(f"Profiles not found at {profiles_path}. Run Module 1-2 first.")
        sys.exit(1)

    with open(profiles_path) as f:
        profiles = json.load(f)

    logger.info(f"Loaded {len(profiles)} profiles from Module 2")

    # Load data
    if args.mode == "full":
        adata_gex, adata_protein = load_xenium_singlecell()
    else:
        adata_full, _ = load_xenium_singlecell(max_cells=1000)
        bounds = get_quadrant_bounds(adata_full.obsm["spatial"], args.quadrant_id)
        adata_gex, adata_protein = load_xenium_singlecell(region_bounds=bounds)

    logger.info(f"Loaded {adata_gex.shape[0]:,} cells")

    # Extract arrays
    X_protein = np.asarray(adata_protein.X.todense() if hasattr(adata_protein.X, "todense") else adata_protein.X)
    marker_names = list(adata_protein.var_names)

    # Assign cells to profiles
    logger.info("Assigning cells to profiles...")
    assignments, scores, profile_names = assign_cells_to_profiles(
        X_protein=X_protein,
        marker_names=marker_names,
        profiles=profiles,
        method="hard",
    )

    # Save assignments
    assignments_df = pd.DataFrame({
        "cell_id": adata_gex.obs_names,
        "profile": [profile_names[i] for i in assignments],
        "profile_score": scores.max(axis=1),
    })
    assignments_df.to_csv(output_dir / "cell_assignments.csv", index=False)

    # Count cells per profile
    for i, name in enumerate(profile_names):
        n = (assignments == i).sum()
        logger.info(f"  {name}: {n:,} cells ({100*n/len(assignments):.1f}%)")

    # Create output directory for programs
    programs_dir = output_dir / "module4_programs"
    programs_dir.mkdir(exist_ok=True)

    # Run NMF for each cell type
    all_results = {}
    for i, profile_name in enumerate(profile_names):
        logger.info("=" * 70)
        logger.info(f"Processing: {profile_name}")
        logger.info("=" * 70)

        cell_mask = assignments == i
        result = run_nmf_for_celltype(
            adata_gex=adata_gex,
            cell_mask=cell_mask,
            n_programs=args.n_programs,
            min_cells=500,
        )

        if result is not None:
            all_results[profile_name] = result

            # Save per-cell-type results
            with open(programs_dir / f"{profile_name.replace(' ', '_')}_programs.json", "w") as f:
                json.dump(result, f, indent=2)

    # Save summary
    summary = {
        "mode": args.mode,
        "quadrant_id": args.quadrant_id,
        "n_cells_total": int(adata_gex.shape[0]),
        "n_profiles": len(profile_names),
        "profile_names": profile_names,
        "n_programs_per_type": args.n_programs,
        "results": {
            name: {
                "n_cells": r["n_cells"],
                "n_programs": r["n_programs"],
                "mean_morans_i": np.mean([p["morans_i"] for p in r["programs"]]),
            }
            for name, r in all_results.items()
        },
    }
    with open(output_dir / "module4_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    logger.info("=" * 70)
    logger.info("MODULE 4 COMPLETE")
    logger.info(f"Output saved to: {output_dir}")
    logger.info("=" * 70)


if __name__ == "__main__":
    main()
