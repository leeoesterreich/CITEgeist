#!/usr/bin/env python3
"""
Pan-Region Module 4: NMF Program Discovery with Consensus Profiles

Applies a unified consensus profile set across all 465K cells from the
Xenium RCC dataset to discover gene expression programs per cell type.

This enables cross-region comparison of transcriptomic programs within
consistent cell type definitions.
"""

import argparse
import json
import logging
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
from scipy.sparse import issparse
from sklearn.decomposition import NMF
from sklearn.neighbors import NearestNeighbors
import warnings

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Add parent paths for imports
import sys
sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent / "xenium_pseudovisium" / "src"))
from load_xenium_singlecell import load_xenium_singlecell

# Consensus profiles derived from Q0, Q2, Q3 analysis
CONSENSUS_PROFILES = {
    "Exhausted_CD8_T": ["LAG-3", "CD3E", "CD8A", "PD-1", "CD45RO", "CD45"],
    "Macrophages": ["CD68", "CD163", "CD16", "HLA-DR"],
    "B_cells": ["CD20", "CD45RA"],
    "Proliferating_Tumor": ["Ki-67", "PanCK", "PCNA"],
    "Stromal": ["CD31", "alphaSMA", "Vimentin", "PTEN"],
    "Epithelial": ["E-Cadherin", "Beta-catenin"],
    "Plasma_Thelper": ["CD138", "CD4"],
    "Dendritic": ["CD11c"],
    "Checkpoint_pos": ["PD-L1", "VISTA"],
    "Cytotoxic": ["GranzymeB"],
}


def assign_cells_to_profiles(
    X_protein: np.ndarray,
    marker_names: List[str],
    profiles: Dict[str, List[str]],
) -> Tuple[np.ndarray, Dict[str, int]]:
    """
    Assign each cell to the best-matching profile based on protein expression.

    Uses normalized mean expression of profile markers.

    Returns:
        assignments: Array of profile indices per cell
        profile_to_idx: Mapping of profile name to index
    """
    n_cells = X_protein.shape[0]
    profile_names = list(profiles.keys())
    profile_to_idx = {name: i for i, name in enumerate(profile_names)}

    # Compute profile scores for each cell
    scores = np.zeros((n_cells, len(profile_names)))

    for i, (profile_name, markers) in enumerate(profiles.items()):
        # Get indices of markers in this profile
        marker_indices = []
        for m in markers:
            if m in marker_names:
                marker_indices.append(marker_names.index(m))

        if len(marker_indices) == 0:
            logger.warning(f"No markers found for profile {profile_name}")
            continue

        # Mean expression of profile markers (normalized per marker)
        profile_expr = X_protein[:, marker_indices]

        # Z-score normalize each marker
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            profile_expr_norm = (profile_expr - np.mean(profile_expr, axis=0)) / (np.std(profile_expr, axis=0) + 1e-8)

        # Score = mean of z-scored markers
        scores[:, i] = np.mean(profile_expr_norm, axis=1)

    # Assign to highest-scoring profile
    assignments = np.argmax(scores, axis=1)

    return assignments, profile_to_idx


def compute_morans_i(values: np.ndarray, coords: np.ndarray, k: int = 15) -> float:
    """Compute Moran's I spatial autocorrelation."""
    n = len(values)
    if n < k + 1:
        return 0.0

    # Build neighbor graph
    nn = NearestNeighbors(n_neighbors=k + 1)
    nn.fit(coords)
    distances, indices = nn.kneighbors(coords)

    # Compute Moran's I
    mean_val = np.mean(values)
    numerator = 0.0
    denominator = np.sum((values - mean_val) ** 2)

    if denominator < 1e-10:
        return 0.0

    W = 0.0
    for i in range(n):
        for j in indices[i, 1:]:  # Skip self
            w_ij = 1.0 / (distances[i, list(indices[i]).index(j)] + 1e-8)
            numerator += w_ij * (values[i] - mean_val) * (values[j] - mean_val)
            W += w_ij

    if W < 1e-10:
        return 0.0

    morans_i = (n / W) * (numerator / denominator)
    return morans_i


def run_nmf_for_profile(
    X_gex: np.ndarray,
    coords: np.ndarray,
    gene_names: List[str],
    profile_name: str,
    n_programs: int = 5,
    min_cells: int = 1000,
) -> Dict:
    """
    Run NMF to discover gene expression programs for a cell type.

    Returns dict with program info including top genes and Moran's I.
    """
    n_cells = X_gex.shape[0]

    if n_cells < min_cells:
        logger.warning(f"{profile_name}: Only {n_cells} cells, skipping (min: {min_cells})")
        return None

    logger.info(f"{profile_name}: Running NMF with K={n_programs} on {n_cells:,} cells x {X_gex.shape[1]} genes")

    # Filter to variable genes
    gene_var = np.var(X_gex, axis=0)
    top_var_idx = np.argsort(gene_var)[-500:]  # Top 500 variable genes
    X_filtered = X_gex[:, top_var_idx]
    filtered_genes = [gene_names[i] for i in top_var_idx]

    # Ensure non-negative
    X_filtered = np.maximum(X_filtered, 0)

    # Run NMF
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        nmf = NMF(n_components=n_programs, init='nndsvda', random_state=42, max_iter=500)
        H = nmf.fit_transform(X_filtered)  # Cell x Program
        W = nmf.components_  # Program x Gene

    # Extract program info
    programs = []
    for k in range(n_programs):
        # Top genes for this program
        gene_loadings = W[k, :]
        top_gene_idx = np.argsort(gene_loadings)[-10:][::-1]
        top_genes = [filtered_genes[i] for i in top_gene_idx]
        top_loadings = gene_loadings[top_gene_idx].tolist()

        # Moran's I for program activity
        program_activity = H[:, k]
        morans_i = compute_morans_i(program_activity, coords, k=15)

        programs.append({
            "program_id": k,
            "top_genes": top_genes,
            "top_loadings": top_loadings,
            "morans_i": morans_i,
            "variance_explained": np.var(program_activity),
        })

    return {
        "profile_name": profile_name,
        "n_cells": n_cells,
        "n_programs": n_programs,
        "programs": programs,
        "mean_morans_i": np.mean([p["morans_i"] for p in programs]),
        "H_matrix_shape": H.shape,
    }


def main():
    parser = argparse.ArgumentParser(description="Pan-region Module 4 with consensus profiles")
    parser.add_argument("--output-dir", type=str,
                        default="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/.worktrees/singlecell-demonstration/Benchmarking/xenium_benchmarking/CITEgeist/output_singlecell_demonstration/panregion")
    parser.add_argument("--n-programs", type=int, default=5, help="Number of NMF programs per cell type")
    parser.add_argument("--min-cells", type=int, default=1000, help="Minimum cells per profile")
    parser.add_argument("--subsample", type=int, default=None, help="Subsample cells for testing")
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("=" * 70)
    logger.info("PAN-REGION MODULE 4: Consensus Profile Program Discovery")
    logger.info("=" * 70)

    # Save consensus profiles
    with open(output_dir / "consensus_profiles.json", "w") as f:
        json.dump(CONSENSUS_PROFILES, f, indent=2)
    logger.info(f"Using {len(CONSENSUS_PROFILES)} consensus profiles")

    # Load full dataset
    logger.info("Loading full Xenium dataset (all 465K cells)...")
    adata_gex, adata_protein = load_xenium_singlecell(
        max_cells=args.subsample,  # None = full dataset, or subsample for testing
        region_bounds=None,  # Full dataset
    )

    n_cells = adata_gex.shape[0]
    logger.info(f"Loaded {n_cells:,} cells × {adata_gex.shape[1]} genes × {adata_protein.shape[1]} proteins")

    # Get data matrices
    X_gex = adata_gex.X.toarray() if issparse(adata_gex.X) else adata_gex.X
    X_protein = adata_protein.X.toarray() if issparse(adata_protein.X) else adata_protein.X
    coords = adata_gex.obsm['spatial']
    gene_names = list(adata_gex.var_names)
    marker_names = list(adata_protein.var_names)

    # Log-normalize gene expression
    X_gex = np.log1p(X_gex)

    # Assign cells to consensus profiles
    logger.info("Assigning cells to consensus profiles...")
    assignments, profile_to_idx = assign_cells_to_profiles(
        X_protein, marker_names, CONSENSUS_PROFILES
    )

    # Count cells per profile
    idx_to_profile = {v: k for k, v in profile_to_idx.items()}
    profile_counts = {}
    for i in range(len(CONSENSUS_PROFILES)):
        count = np.sum(assignments == i)
        profile_name = idx_to_profile[i]
        profile_counts[profile_name] = count
        pct = 100 * count / n_cells
        logger.info(f"  {profile_name}: {count:,} cells ({pct:.1f}%)")

    # Save cell assignments
    cell_assignments_df = pd.DataFrame({
        "cell_id": adata_gex.obs_names,
        "x": coords[:, 0],
        "y": coords[:, 1],
        "profile_idx": assignments,
        "profile_name": [idx_to_profile[a] for a in assignments],
    })
    cell_assignments_df.to_csv(output_dir / "cell_assignments.csv", index=False)
    logger.info(f"Saved cell assignments to {output_dir / 'cell_assignments.csv'}")

    # Run NMF for each profile
    logger.info("=" * 70)
    logger.info("Running NMF Program Discovery per Cell Type")
    logger.info("=" * 70)

    all_results = {}
    programs_dir = output_dir / "programs"
    programs_dir.mkdir(exist_ok=True)

    for profile_name in CONSENSUS_PROFILES.keys():
        logger.info(f"\n{'=' * 70}")
        logger.info(f"Processing: {profile_name}")
        logger.info("=" * 70)

        # Get cells for this profile
        profile_idx = profile_to_idx[profile_name]
        cell_mask = assignments == profile_idx

        if np.sum(cell_mask) < args.min_cells:
            logger.warning(f"Skipping {profile_name}: only {np.sum(cell_mask)} cells")
            continue

        # Subset data
        X_profile = X_gex[cell_mask, :]
        coords_profile = coords[cell_mask, :]

        # Run NMF
        result = run_nmf_for_profile(
            X_profile, coords_profile, gene_names,
            profile_name, n_programs=args.n_programs, min_cells=args.min_cells
        )

        if result is not None:
            all_results[profile_name] = result

            # Save per-profile results
            with open(programs_dir / f"{profile_name}_programs.json", "w") as f:
                # Convert numpy types to python types
                result_serializable = {
                    k: (v.tolist() if isinstance(v, np.ndarray) else v)
                    for k, v in result.items()
                }
                result_serializable["programs"] = [
                    {k: (v.tolist() if isinstance(v, np.ndarray) else v) for k, v in p.items()}
                    for p in result["programs"]
                ]
                json.dump(result_serializable, f, indent=2)

    # Summary
    logger.info("\n" + "=" * 70)
    logger.info("PAN-REGION MODULE 4 COMPLETE")
    logger.info("=" * 70)

    summary = {
        "n_cells_total": int(n_cells),
        "n_profiles": len(CONSENSUS_PROFILES),
        "profiles_analyzed": len(all_results),
        "results": {
            name: {
                "n_cells": result["n_cells"],
                "n_programs": result["n_programs"],
                "mean_morans_i": result["mean_morans_i"],
            }
            for name, result in all_results.items()
        }
    }

    with open(output_dir / "panregion_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    logger.info(f"Results saved to: {output_dir}")
    logger.info(f"Profiles analyzed: {len(all_results)}/{len(CONSENSUS_PROFILES)}")

    # Print program highlights
    logger.info("\nTop Programs by Moran's I:")
    for profile_name, result in all_results.items():
        best_program = max(result["programs"], key=lambda p: p["morans_i"])
        logger.info(f"  {profile_name}: Program {best_program['program_id']} "
                   f"(I={best_program['morans_i']:.3f}) - {', '.join(best_program['top_genes'][:5])}")


if __name__ == "__main__":
    main()
