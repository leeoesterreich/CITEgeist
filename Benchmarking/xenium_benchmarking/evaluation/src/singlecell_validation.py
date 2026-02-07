#!/usr/bin/env python3
"""
Single-cell validation analysis for CITEgeist autodiscovery.

This script validates autodiscovered profiles at the single-cell level using
Xenium protein expression data. Key analyses:

1. Cell-level marker co-expression validation (e.g., EMT cells = Vimentin + E-Cadherin)
2. Compare autodiscovered profile assignments to artificial GT definitions
3. Spatial coherence of profile-assigned cells
4. Quantify how many cells are "misclassified" by artificial GT

Usage:
    python singlecell_validation.py --output_dir /path/to/output
"""

import argparse
import json
import logging
import sys
from pathlib import Path
from typing import Dict, List, Tuple, Optional

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.spatial.distance import cdist
from sklearn.mixture import GaussianMixture

# Add paths
REPO_ROOT = Path(__file__).parent.parent.parent.parent.parent
sys.path.insert(0, str(REPO_ROOT / "Benchmarking/xenium_pseudovisium/src"))

from load_xenium import load_xenium_data, split_gex_protein

logger = logging.getLogger(__name__)

# Xenium data location
XENIUM_DATA_DIR = Path("/ix1/alee/LO_LAB/General/Public_Data/10x_PublicData/Xenium_RNA_Proteomic_RenalCellCarcinoma")

# Autodiscovered profiles (spatially-validated from colocalization analysis)
AUTODISCOVERED_PROFILES = {
    "B cells": ["CD20", "CD45RA"],
    "Macrophages": ["CD68", "CD16"],
    "T cells": ["CD3E", "CD45RO"],
    "Exhausted T cells": ["CD8A", "PD-1"],
    "M2 Macrophages": ["CD163"],
    "Endothelial": ["CD31"],
    "Epithelial": ["PanCK"],
    "EMT cells": ["Vimentin", "E-Cadherin"],
}

# Artificial GT profiles (from manual definitions)
ARTIFICIAL_GT_PROFILES = {
    "B cells": ["CD20"],
    "CD4+ T cells": ["CD3E", "CD4"],
    "CD8+ T cells": ["CD3E", "CD8A"],
    "Macrophages": ["CD68"],
    "Endothelial": ["CD31"],
    "Epithelial": ["PanCK", "E-Cadherin"],
    "Fibroblasts": ["alphaSMA", "Vimentin"],
}

# Key validation pairs (marker1, marker2, expected_relationship)
VALIDATION_PAIRS = [
    ("Vimentin", "E-Cadherin", "co-express", "EMT cells"),
    ("PanCK", "E-Cadherin", "independent", "Epithelial definition"),
    ("alphaSMA", "Vimentin", "co-express", "Fibroblasts"),
    ("CD68", "CD16", "co-express", "Macrophages"),
    ("CD20", "CD45RA", "co-express", "B cells"),
    ("CD3E", "CD45RO", "co-express", "T cells"),
    ("CD8A", "PD-1", "co-express", "Exhausted T"),
]


def get_protein_thresholds(
    protein_matrix: np.ndarray,
    protein_names: List[str],
    method: str = "gmm"
) -> Dict[str, float]:
    """
    Compute positive/negative thresholds for each protein marker.

    Args:
        protein_matrix: (n_cells, n_proteins) expression matrix
        protein_names: List of protein names
        method: "gmm" for Gaussian mixture, "percentile" for 75th percentile

    Returns:
        Dict mapping protein name to threshold
    """
    thresholds = {}

    for i, protein in enumerate(protein_names):
        expr = protein_matrix[:, i]
        expr_nonzero = expr[expr > 0]

        if len(expr_nonzero) < 100:
            # Not enough signal, use simple percentile
            thresholds[protein] = np.percentile(expr[expr > 0], 75) if np.sum(expr > 0) > 0 else 0
            continue

        if method == "gmm":
            try:
                gmm = GaussianMixture(n_components=2, random_state=42)
                gmm.fit(expr_nonzero.reshape(-1, 1))

                # Threshold at intersection of two components
                means = gmm.means_.flatten()
                stds = np.sqrt(gmm.covariances_.flatten())

                # Find intersection point
                if means[0] < means[1]:
                    low_mean, high_mean = means[0], means[1]
                    low_std, high_std = stds[0], stds[1]
                else:
                    low_mean, high_mean = means[1], means[0]
                    low_std, high_std = stds[1], stds[0]

                # Use mean of the two component means as threshold
                thresholds[protein] = (low_mean + high_mean) / 2
            except Exception:
                thresholds[protein] = np.percentile(expr_nonzero, 75)
        else:
            thresholds[protein] = np.percentile(expr_nonzero, 75)

    return thresholds


def score_cells_by_profiles(
    protein_matrix: np.ndarray,
    protein_names: List[str],
    profiles: Dict[str, List[str]],
    thresholds: Dict[str, float],
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Score each cell for each profile based on marker expression.

    Args:
        protein_matrix: (n_cells, n_proteins) expression matrix
        protein_names: List of protein names
        profiles: Dict mapping profile name to list of marker names
        thresholds: Dict mapping protein name to positive threshold

    Returns:
        scores: (n_cells, n_profiles) continuous profile scores (mean z-score of markers)
        assignments: (n_cells,) hard assignment to best-scoring profile
    """
    n_cells = protein_matrix.shape[0]
    profile_names = list(profiles.keys())
    n_profiles = len(profile_names)

    # Z-score normalize each protein
    protein_zscore = np.zeros_like(protein_matrix)
    for i in range(protein_matrix.shape[1]):
        col = protein_matrix[:, i]
        if col.std() > 0:
            protein_zscore[:, i] = (col - col.mean()) / col.std()

    # Create protein name to index mapping
    name_to_idx = {name: i for i, name in enumerate(protein_names)}

    # Score each cell for each profile
    scores = np.zeros((n_cells, n_profiles))

    for p_idx, (profile_name, markers) in enumerate(profiles.items()):
        # Get indices of markers in this profile
        marker_indices = [name_to_idx[m] for m in markers if m in name_to_idx]

        if len(marker_indices) == 0:
            continue

        # Profile score = mean z-score of profile markers
        profile_scores = protein_zscore[:, marker_indices].mean(axis=1)
        scores[:, p_idx] = profile_scores

    # Hard assignment = argmax (with tie-breaking)
    assignments = np.argmax(scores, axis=1)

    return scores, assignments


def classify_cells_binary(
    protein_matrix: np.ndarray,
    protein_names: List[str],
    profiles: Dict[str, List[str]],
    thresholds: Dict[str, float],
) -> Dict[str, np.ndarray]:
    """
    Classify cells using binary positive/negative for each profile.

    A cell is positive for a profile if ALL major markers exceed threshold.

    Returns:
        Dict mapping profile name to boolean array of positive cells
    """
    name_to_idx = {name: i for i, name in enumerate(protein_names)}

    positive_cells = {}
    for profile_name, markers in profiles.items():
        # Get marker columns
        marker_indices = [name_to_idx[m] for m in markers if m in name_to_idx]

        if len(marker_indices) == 0:
            positive_cells[profile_name] = np.zeros(protein_matrix.shape[0], dtype=bool)
            continue

        # Cell is positive if ALL markers exceed threshold
        is_positive = np.ones(protein_matrix.shape[0], dtype=bool)
        for m_idx, marker in zip(marker_indices, [m for m in markers if m in name_to_idx]):
            threshold = thresholds.get(marker, 0)
            is_positive &= (protein_matrix[:, m_idx] > threshold)

        positive_cells[profile_name] = is_positive

    return positive_cells


def analyze_marker_coexpression(
    protein_matrix: np.ndarray,
    protein_names: List[str],
    marker1: str,
    marker2: str,
    thresholds: Dict[str, float],
) -> Dict:
    """
    Analyze co-expression of two markers at single-cell level.

    Returns:
        Dict with:
        - n_double_positive: cells positive for both
        - n_single_positive_1: cells positive for marker1 only
        - n_single_positive_2: cells positive for marker2 only
        - n_double_negative: cells negative for both
        - pearson_r: correlation coefficient
        - contingency: 2x2 contingency table
        - odds_ratio: association strength
    """
    name_to_idx = {name: i for i, name in enumerate(protein_names)}

    if marker1 not in name_to_idx or marker2 not in name_to_idx:
        return {"error": f"Marker not found: {marker1} or {marker2}"}

    expr1 = protein_matrix[:, name_to_idx[marker1]]
    expr2 = protein_matrix[:, name_to_idx[marker2]]

    thresh1 = thresholds.get(marker1, 0)
    thresh2 = thresholds.get(marker2, 0)

    pos1 = expr1 > thresh1
    pos2 = expr2 > thresh2

    # Contingency table
    double_pos = pos1 & pos2
    single_pos1 = pos1 & ~pos2
    single_pos2 = ~pos1 & pos2
    double_neg = ~pos1 & ~pos2

    n_double_pos = double_pos.sum()
    n_single_pos1 = single_pos1.sum()
    n_single_pos2 = single_pos2.sum()
    n_double_neg = double_neg.sum()

    # Pearson correlation (on log-transformed values for expression)
    expr1_log = np.log1p(expr1)
    expr2_log = np.log1p(expr2)
    pearson_r, pearson_p = stats.pearsonr(expr1_log, expr2_log)

    # Odds ratio
    contingency = np.array([[n_double_pos, n_single_pos2],
                           [n_single_pos1, n_double_neg]])

    # Add small constant to avoid division by zero
    odds_ratio = (n_double_pos * n_double_neg + 1) / (n_single_pos1 * n_single_pos2 + 1)

    return {
        "marker1": marker1,
        "marker2": marker2,
        "n_cells": len(expr1),
        "n_double_positive": int(n_double_pos),
        "n_single_positive_1": int(n_single_pos1),
        "n_single_positive_2": int(n_single_pos2),
        "n_double_negative": int(n_double_neg),
        "frac_double_positive": n_double_pos / len(expr1),
        "pearson_r": float(pearson_r),
        "pearson_p": float(pearson_p),
        "odds_ratio": float(odds_ratio),
        "contingency": contingency.tolist(),
    }


def compute_spatial_coherence(
    coords: np.ndarray,
    cell_mask: np.ndarray,
    n_neighbors: int = 10,
) -> float:
    """
    Compute spatial coherence (local enrichment) for a set of cells.

    Measures whether positive cells cluster together more than expected by chance.

    Args:
        coords: (n_cells, 2) spatial coordinates
        cell_mask: boolean array indicating cells of interest
        n_neighbors: number of neighbors to consider

    Returns:
        coherence_score: ratio of observed to expected neighbor fraction
    """
    if cell_mask.sum() < 10:
        return np.nan

    # Sample for efficiency
    n_sample = min(5000, coords.shape[0])
    sample_idx = np.random.choice(coords.shape[0], n_sample, replace=False)

    coords_sample = coords[sample_idx]
    mask_sample = cell_mask[sample_idx]

    # Compute pairwise distances
    dists = cdist(coords_sample, coords_sample)

    # For each cell, find k nearest neighbors
    observed_frac = []
    for i in range(n_sample):
        if not mask_sample[i]:
            continue

        # Get k nearest neighbors (excluding self)
        neighbor_idx = np.argsort(dists[i])[1:n_neighbors+1]
        neighbor_mask = mask_sample[neighbor_idx]
        observed_frac.append(neighbor_mask.mean())

    if len(observed_frac) == 0:
        return np.nan

    # Expected fraction (global rate)
    expected_frac = mask_sample.mean()

    # Coherence = observed / expected
    coherence = np.mean(observed_frac) / expected_frac if expected_frac > 0 else np.nan

    return coherence


def plot_marker_coexpression(
    protein_matrix: np.ndarray,
    protein_names: List[str],
    marker1: str,
    marker2: str,
    thresholds: Dict[str, float],
    output_path: Path,
    title: str = None,
):
    """
    Create scatter plot of marker co-expression with quadrant analysis.
    """
    name_to_idx = {name: i for i, name in enumerate(protein_names)}

    expr1 = protein_matrix[:, name_to_idx[marker1]]
    expr2 = protein_matrix[:, name_to_idx[marker2]]

    thresh1 = thresholds.get(marker1, 0)
    thresh2 = thresholds.get(marker2, 0)

    # Sample for plotting
    n_sample = min(50000, len(expr1))
    sample_idx = np.random.choice(len(expr1), n_sample, replace=False)

    expr1_sample = np.log1p(expr1[sample_idx])
    expr2_sample = np.log1p(expr2[sample_idx])

    # Create plot
    fig, ax = plt.subplots(figsize=(8, 8))

    # Scatter with density coloring
    ax.scatter(expr1_sample, expr2_sample, alpha=0.1, s=1, c='steelblue')

    # Add threshold lines
    ax.axvline(np.log1p(thresh1), color='red', linestyle='--', alpha=0.7, label=f'{marker1} threshold')
    ax.axhline(np.log1p(thresh2), color='green', linestyle='--', alpha=0.7, label=f'{marker2} threshold')

    # Compute quadrant counts
    pos1 = expr1 > thresh1
    pos2 = expr2 > thresh2
    n_total = len(expr1)

    q_labels = [
        f"Double+: {(pos1 & pos2).sum():,} ({100*(pos1 & pos2).sum()/n_total:.1f}%)",
        f"{marker1}+ only: {(pos1 & ~pos2).sum():,} ({100*(pos1 & ~pos2).sum()/n_total:.1f}%)",
        f"{marker2}+ only: {(~pos1 & pos2).sum():,} ({100*(~pos1 & pos2).sum()/n_total:.1f}%)",
        f"Double-: {(~pos1 & ~pos2).sum():,} ({100*(~pos1 & ~pos2).sum()/n_total:.1f}%)",
    ]

    # Add quadrant labels
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()

    ax.text(xlim[1]*0.95, ylim[1]*0.95, q_labels[0], ha='right', va='top', fontsize=10,
            bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.8))
    ax.text(xlim[1]*0.95, ylim[0]*0.05 + np.log1p(thresh2)*0.5, q_labels[1], ha='right', va='bottom', fontsize=10,
            bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))
    ax.text(xlim[0]*0.05 + np.log1p(thresh1)*0.5, ylim[1]*0.95, q_labels[2], ha='left', va='top', fontsize=10,
            bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))

    # Correlation
    r, p = stats.pearsonr(expr1_sample, expr2_sample)
    ax.set_xlabel(f'{marker1} (log1p)', fontsize=12)
    ax.set_ylabel(f'{marker2} (log1p)', fontsize=12)

    if title:
        ax.set_title(f'{title}\nPearson r = {r:.3f}', fontsize=14)
    else:
        ax.set_title(f'{marker1} vs {marker2} (n={n_total:,} cells)\nPearson r = {r:.3f}', fontsize=14)

    ax.legend(loc='lower right')

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()

    logger.info(f"Saved: {output_path}")


def plot_spatial_cells(
    coords: np.ndarray,
    cell_mask: np.ndarray,
    output_path: Path,
    title: str,
    n_sample: int = 50000,
):
    """
    Plot spatial distribution of cells matching a profile.
    """
    # Sample for plotting
    if coords.shape[0] > n_sample:
        sample_idx = np.random.choice(coords.shape[0], n_sample, replace=False)
        coords_sample = coords[sample_idx]
        mask_sample = cell_mask[sample_idx]
    else:
        coords_sample = coords
        mask_sample = cell_mask

    fig, ax = plt.subplots(figsize=(10, 10))

    # Plot all cells in gray
    ax.scatter(coords_sample[~mask_sample, 0], coords_sample[~mask_sample, 1],
               c='lightgray', s=0.5, alpha=0.3, label='Other cells')

    # Highlight positive cells
    ax.scatter(coords_sample[mask_sample, 0], coords_sample[mask_sample, 1],
               c='red', s=2, alpha=0.5, label=f'Positive ({mask_sample.sum():,} cells)')

    ax.set_xlabel('X coordinate')
    ax.set_ylabel('Y coordinate')
    ax.set_title(f'{title}\n{cell_mask.sum():,} / {len(cell_mask):,} cells ({100*cell_mask.mean():.1f}%)')
    ax.legend(loc='upper right')
    ax.set_aspect('equal')

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()

    logger.info(f"Saved: {output_path}")


def plot_profile_comparison(
    autodiscovery_positive: Dict[str, np.ndarray],
    artificial_positive: Dict[str, np.ndarray],
    output_path: Path,
):
    """
    Compare cell classifications between autodiscovered and artificial profiles.
    """
    # Find cells classified differently
    auto_profiles = list(autodiscovery_positive.keys())
    art_profiles = list(artificial_positive.keys())

    # Create comparison matrix
    data = []
    for auto_name, auto_mask in autodiscovery_positive.items():
        row = {"Autodiscovered": auto_name, "n_cells": int(auto_mask.sum())}
        for art_name, art_mask in artificial_positive.items():
            overlap = (auto_mask & art_mask).sum()
            row[f"Overlap_{art_name}"] = int(overlap)
        data.append(row)

    df = pd.DataFrame(data)

    # Create heatmap
    overlap_cols = [c for c in df.columns if c.startswith("Overlap_")]
    overlap_matrix = df[overlap_cols].values

    fig, ax = plt.subplots(figsize=(12, 8))

    # Normalize by row (autodiscovered profile)
    row_sums = overlap_matrix.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1
    overlap_norm = overlap_matrix / row_sums

    sns.heatmap(overlap_norm,
                xticklabels=[c.replace("Overlap_", "") for c in overlap_cols],
                yticklabels=df["Autodiscovered"].values,
                annot=True, fmt=".2f", cmap="YlOrRd",
                ax=ax)

    ax.set_xlabel("Artificial GT Profile")
    ax.set_ylabel("Autodiscovered Profile")
    ax.set_title("Profile Assignment Comparison\n(Fraction of autodiscovered cells matching artificial GT)")

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()

    logger.info(f"Saved: {output_path}")

    return df


def main(output_dir: Path):
    """Run single-cell validation analysis."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    figures_dir = output_dir / "figures"
    figures_dir.mkdir(exist_ok=True)

    logger.info("Loading Xenium single-cell data...")
    adata = load_xenium_data(XENIUM_DATA_DIR)

    logger.info("Splitting into GEX and protein...")
    adata_gex, adata_protein = split_gex_protein(adata)

    # Get protein matrix
    protein_matrix = adata_protein.X
    if hasattr(protein_matrix, 'toarray'):
        protein_matrix = protein_matrix.toarray()
    protein_names = list(adata_protein.var_names)
    coords = adata_protein.obsm["spatial"]

    logger.info(f"Protein matrix: {protein_matrix.shape}")
    logger.info(f"Proteins: {protein_names}")

    # Compute thresholds
    logger.info("Computing protein thresholds...")
    thresholds = get_protein_thresholds(protein_matrix, protein_names, method="gmm")

    # Save thresholds
    with open(output_dir / "protein_thresholds.json", "w") as f:
        json.dump(thresholds, f, indent=2)

    results = {
        "n_cells": int(protein_matrix.shape[0]),
        "n_proteins": int(protein_matrix.shape[1]),
        "proteins": protein_names,
        "thresholds": thresholds,
    }

    # =========================================================================
    # Analysis 1: Marker co-expression validation
    # =========================================================================
    logger.info("\n" + "="*60)
    logger.info("Analysis 1: Marker Co-expression Validation")
    logger.info("="*60)

    coexpression_results = []
    for marker1, marker2, expected, description in VALIDATION_PAIRS:
        logger.info(f"\nAnalyzing {marker1} vs {marker2} ({description})...")

        coexp = analyze_marker_coexpression(
            protein_matrix, protein_names, marker1, marker2, thresholds
        )
        coexp["expected"] = expected
        coexp["description"] = description
        coexpression_results.append(coexp)

        logger.info(f"  Double positive: {coexp['n_double_positive']:,} ({100*coexp['frac_double_positive']:.2f}%)")
        logger.info(f"  Pearson r: {coexp['pearson_r']:.3f}")
        logger.info(f"  Odds ratio: {coexp['odds_ratio']:.2f}")

        # Generate scatter plot
        plot_marker_coexpression(
            protein_matrix, protein_names, marker1, marker2, thresholds,
            figures_dir / f"coexpression_{marker1}_{marker2}.png",
            title=f"{description}: {marker1} vs {marker2}"
        )

    results["coexpression_analysis"] = coexpression_results

    # =========================================================================
    # Analysis 2: EMT cell validation (key finding)
    # =========================================================================
    logger.info("\n" + "="*60)
    logger.info("Analysis 2: EMT Cell Validation")
    logger.info("="*60)

    # Find EMT cells (Vimentin+ AND E-Cadherin+)
    vim_idx = protein_names.index("Vimentin")
    ecad_idx = protein_names.index("E-Cadherin")

    vim_pos = protein_matrix[:, vim_idx] > thresholds["Vimentin"]
    ecad_pos = protein_matrix[:, ecad_idx] > thresholds["E-Cadherin"]

    emt_cells = vim_pos & ecad_pos
    vim_only = vim_pos & ~ecad_pos
    ecad_only = ~vim_pos & ecad_pos

    emt_results = {
        "n_emt_cells": int(emt_cells.sum()),
        "frac_emt_cells": float(emt_cells.mean()),
        "n_vimentin_only": int(vim_only.sum()),
        "n_ecadherin_only": int(ecad_only.sum()),
        "spatial_coherence_emt": compute_spatial_coherence(coords, emt_cells),
        "spatial_coherence_vim_only": compute_spatial_coherence(coords, vim_only),
        "spatial_coherence_ecad_only": compute_spatial_coherence(coords, ecad_only),
    }

    logger.info(f"EMT cells (Vim+ ECad+): {emt_results['n_emt_cells']:,} ({100*emt_results['frac_emt_cells']:.2f}%)")
    logger.info(f"Vimentin only: {emt_results['n_vimentin_only']:,}")
    logger.info(f"E-Cadherin only: {emt_results['n_ecadherin_only']:,}")
    logger.info(f"EMT spatial coherence: {emt_results['spatial_coherence_emt']:.2f}")

    results["emt_validation"] = emt_results

    # Plot EMT cells spatially
    plot_spatial_cells(
        coords, emt_cells,
        figures_dir / "spatial_emt_cells.png",
        "EMT Cells (Vimentin+ E-Cadherin+)"
    )

    plot_spatial_cells(
        coords, vim_only,
        figures_dir / "spatial_vimentin_only.png",
        "Vimentin+ E-Cadherin- (Traditional 'Fibroblast')"
    )

    # =========================================================================
    # Analysis 3: Profile classification comparison
    # =========================================================================
    logger.info("\n" + "="*60)
    logger.info("Analysis 3: Profile Classification Comparison")
    logger.info("="*60)

    # Classify cells by autodiscovered profiles
    autodiscovery_positive = classify_cells_binary(
        protein_matrix, protein_names, AUTODISCOVERED_PROFILES, thresholds
    )

    # Classify cells by artificial GT profiles
    artificial_positive = classify_cells_binary(
        protein_matrix, protein_names, ARTIFICIAL_GT_PROFILES, thresholds
    )

    # Log counts
    logger.info("\nAutodiscovered profile cell counts:")
    for name, mask in autodiscovery_positive.items():
        coherence = compute_spatial_coherence(coords, mask)
        logger.info(f"  {name}: {mask.sum():,} cells, coherence={coherence:.2f}")

    logger.info("\nArtificial GT profile cell counts:")
    for name, mask in artificial_positive.items():
        coherence = compute_spatial_coherence(coords, mask)
        logger.info(f"  {name}: {mask.sum():,} cells, coherence={coherence:.2f}")

    # Compare classifications
    comparison_df = plot_profile_comparison(
        autodiscovery_positive, artificial_positive,
        figures_dir / "profile_comparison_heatmap.png"
    )
    comparison_df.to_csv(output_dir / "profile_comparison.csv", index=False)

    # =========================================================================
    # Analysis 4: Cells misclassified by artificial GT
    # =========================================================================
    logger.info("\n" + "="*60)
    logger.info("Analysis 4: Misclassification Analysis")
    logger.info("="*60)

    # Key finding: EMT cells would be split between Epithelial and Fibroblasts in artificial GT
    emt_in_fibroblasts = emt_cells & artificial_positive.get("Fibroblasts", np.zeros(len(emt_cells), dtype=bool))
    emt_in_epithelial = emt_cells & artificial_positive.get("Epithelial", np.zeros(len(emt_cells), dtype=bool))

    misclass_results = {
        "emt_classified_as_fibroblast": int(emt_in_fibroblasts.sum()),
        "emt_classified_as_epithelial": int(emt_in_epithelial.sum()),
        "total_emt": int(emt_cells.sum()),
        "pct_emt_as_fibroblast": float(100 * emt_in_fibroblasts.sum() / max(emt_cells.sum(), 1)),
        "pct_emt_as_epithelial": float(100 * emt_in_epithelial.sum() / max(emt_cells.sum(), 1)),
    }

    logger.info(f"EMT cells classified as Fibroblast: {misclass_results['emt_classified_as_fibroblast']:,} ({misclass_results['pct_emt_as_fibroblast']:.1f}%)")
    logger.info(f"EMT cells classified as Epithelial: {misclass_results['emt_classified_as_epithelial']:,} ({misclass_results['pct_emt_as_epithelial']:.1f}%)")

    results["misclassification_analysis"] = misclass_results

    # =========================================================================
    # Analysis 5: Soft assignment scoring
    # =========================================================================
    logger.info("\n" + "="*60)
    logger.info("Analysis 5: Soft Assignment Analysis")
    logger.info("="*60)

    # Score cells by autodiscovered profiles
    auto_scores, auto_assignments = score_cells_by_profiles(
        protein_matrix, protein_names, AUTODISCOVERED_PROFILES, thresholds
    )

    # Score cells by artificial profiles
    art_scores, art_assignments = score_cells_by_profiles(
        protein_matrix, protein_names, ARTIFICIAL_GT_PROFILES, thresholds
    )

    # Variance explained by each profile set
    # (Use total score variance as proxy)
    auto_var_explained = auto_scores.var(axis=0).sum()
    art_var_explained = art_scores.var(axis=0).sum()

    soft_results = {
        "autodiscovery_total_score_variance": float(auto_var_explained),
        "artificial_total_score_variance": float(art_var_explained),
        "variance_ratio": float(auto_var_explained / art_var_explained) if art_var_explained > 0 else np.nan,
    }

    logger.info(f"Autodiscovery score variance: {auto_var_explained:.2f}")
    logger.info(f"Artificial score variance: {art_var_explained:.2f}")
    logger.info(f"Ratio (auto/art): {soft_results['variance_ratio']:.2f}")

    results["soft_assignment_analysis"] = soft_results

    # Save all results
    with open(output_dir / "singlecell_validation_results.json", "w") as f:
        json.dump(results, f, indent=2, default=str)

    logger.info(f"\n{'='*60}")
    logger.info(f"Results saved to: {output_dir}")
    logger.info(f"{'='*60}")

    # Print summary
    print("\n" + "="*60)
    print("SINGLE-CELL VALIDATION SUMMARY")
    print("="*60)
    print(f"\nTotal cells analyzed: {results['n_cells']:,}")
    print(f"\nKey Finding - EMT Cells:")
    print(f"  Vimentin+ E-Cadherin+ cells: {emt_results['n_emt_cells']:,} ({100*emt_results['frac_emt_cells']:.2f}%)")
    print(f"  Spatial coherence: {emt_results['spatial_coherence_emt']:.2f}x expected")
    print(f"\nCo-expression Validation:")
    for coexp in coexpression_results:
        status = "VALIDATED" if coexp["odds_ratio"] > 2 else "WEAK"
        print(f"  {coexp['description']}: r={coexp['pearson_r']:.3f}, OR={coexp['odds_ratio']:.1f} [{status}]")
    print(f"\nMisclassification by Artificial GT:")
    print(f"  EMT cells → Fibroblast: {misclass_results['pct_emt_as_fibroblast']:.1f}%")
    print(f"  EMT cells → Epithelial: {misclass_results['pct_emt_as_epithelial']:.1f}%")
    print("="*60)

    return results


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Single-cell validation analysis")
    parser.add_argument(
        "--output_dir",
        type=str,
        default="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/evaluation/singlecell_validation",
        help="Output directory for results",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable verbose logging",
    )

    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO if args.verbose else logging.WARNING,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    )
    logger.setLevel(logging.INFO)

    main(args.output_dir)
