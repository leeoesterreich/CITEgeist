#!/usr/bin/env python3
"""
Spot-level validation: Compare autodiscovered profiles to actual cell counts.

For each pseudo-Visium spot, we know:
1. The actual cells within that spot (from cell_to_spot_mapping)
2. The protein expression of those cells
3. CITEgeist's predicted proportions

This script validates that autodiscovered profiles better capture
the actual cellular composition than artificial GT definitions.

Usage:
    python spot_level_validation.py --output_dir /path/to/output
"""

import argparse
import json
import logging
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from sklearn.mixture import GaussianMixture

# Add paths
REPO_ROOT = Path(__file__).parent.parent.parent.parent.parent
sys.path.insert(0, str(REPO_ROOT / "Benchmarking/xenium_pseudovisium/src"))

from load_xenium import load_xenium_data, split_gex_protein

logger = logging.getLogger(__name__)

# Paths
XENIUM_DATA_DIR = Path("/ix1/alee/LO_LAB/General/Public_Data/10x_PublicData/Xenium_RNA_Proteomic_RenalCellCarcinoma")
PSEUDOVISIUM_DIR = Path("/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_pseudovisium")
AUTODISCOVERY_DIR = Path("/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/CITEgeist/output_autodiscovery")

# Profile definitions
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

ARTIFICIAL_GT_PROFILES = {
    "B cells": ["CD20"],
    "CD4+ T cells": ["CD3E", "CD4"],
    "CD8+ T cells": ["CD3E", "CD8A"],
    "Macrophages": ["CD68"],
    "Endothelial": ["CD31"],
    "Epithelial": ["PanCK", "E-Cadherin"],
    "Fibroblasts": ["alphaSMA", "Vimentin"],
}


def load_cell_to_spot_mapping(region_id: int) -> pd.DataFrame:
    """Load mapping of cells to spots for a region."""
    mapping_path = PSEUDOVISIUM_DIR / "data" / "cell_to_spot_mapping.csv"
    df = pd.read_csv(mapping_path, index_col=0)

    # Filter to this region
    df = df[df["region_id"] == region_id].copy()

    # Filter out unmapped cells (spot_idx == -1)
    df = df[df["spot_idx"] != -1]

    return df


def get_protein_thresholds_gmm(protein_matrix: np.ndarray, protein_names: List[str]) -> Dict[str, float]:
    """Compute positive/negative thresholds using GMM."""
    thresholds = {}

    for i, protein in enumerate(protein_names):
        expr = protein_matrix[:, i]
        expr_nonzero = expr[expr > 0]

        if len(expr_nonzero) < 100:
            thresholds[protein] = np.percentile(expr[expr > 0], 75) if np.sum(expr > 0) > 0 else 0
            continue

        try:
            gmm = GaussianMixture(n_components=2, random_state=42)
            gmm.fit(expr_nonzero.reshape(-1, 1))
            means = gmm.means_.flatten()
            thresholds[protein] = (means[0] + means[1]) / 2
        except Exception:
            thresholds[protein] = np.percentile(expr_nonzero, 75)

    return thresholds


def classify_cells_to_profiles(
    protein_matrix: np.ndarray,
    protein_names: List[str],
    profiles: Dict[str, List[str]],
    thresholds: Dict[str, float],
) -> Dict[str, np.ndarray]:
    """Classify cells as positive/negative for each profile."""
    name_to_idx = {name: i for i, name in enumerate(protein_names)}

    positive_cells = {}
    for profile_name, markers in profiles.items():
        marker_indices = [name_to_idx[m] for m in markers if m in name_to_idx]

        if len(marker_indices) == 0:
            positive_cells[profile_name] = np.zeros(protein_matrix.shape[0], dtype=bool)
            continue

        is_positive = np.ones(protein_matrix.shape[0], dtype=bool)
        for m_idx, marker in zip(marker_indices, [m for m in markers if m in name_to_idx]):
            threshold = thresholds.get(marker, 0)
            is_positive &= (protein_matrix[:, m_idx] > threshold)

        positive_cells[profile_name] = is_positive

    return positive_cells


def compute_spot_proportions_from_cells(
    cell_to_spot: pd.DataFrame,
    cell_profiles: Dict[str, np.ndarray],
    cell_ids: np.ndarray,
) -> pd.DataFrame:
    """
    Compute actual cell type proportions per spot from cell-level classifications.

    Args:
        cell_to_spot: DataFrame with cell_id index and spot_id column
        cell_profiles: Dict mapping profile name to boolean array of positive cells
        cell_ids: Array of cell IDs matching profile arrays

    Returns:
        DataFrame with spot_id index and profile columns containing proportions
    """
    # Create cell_id to index mapping
    cell_id_to_idx = {cid: i for i, cid in enumerate(cell_ids)}

    # Get unique spots
    spots = cell_to_spot["spot_id"].unique()

    # Compute proportions for each spot
    results = []
    for spot_id in spots:
        spot_cells = cell_to_spot[cell_to_spot["spot_id"] == spot_id].index
        n_cells = len(spot_cells)

        if n_cells == 0:
            continue

        row = {"spot_id": spot_id, "n_cells": n_cells}

        for profile_name, pos_mask in cell_profiles.items():
            # Count positive cells in this spot
            n_positive = sum(
                1 for cell_id in spot_cells
                if cell_id in cell_id_to_idx and pos_mask[cell_id_to_idx[cell_id]]
            )
            row[profile_name] = n_positive / n_cells

        results.append(row)

    return pd.DataFrame(results).set_index("spot_id")


def load_citegeist_predictions(region_id: int) -> pd.DataFrame:
    """Load CITEgeist autodiscovery predictions for a region."""
    # Try different file patterns
    patterns = [
        AUTODISCOVERY_DIR / f"region_{region_id}_autodiscovery_proportions.csv",
        AUTODISCOVERY_DIR / f"Xenium_region_{region_id}_cell_prop_finetuned_results.csv",
        AUTODISCOVERY_DIR / f"autodiscovery_region_{region_id}_cell_prop_finetuned_results.csv",
    ]

    for path in patterns:
        if path.exists():
            return pd.read_csv(path, index_col=0)

    raise FileNotFoundError(f"Could not find predictions for region {region_id}")


def compare_proportions(
    actual: pd.DataFrame,
    predicted: pd.DataFrame,
    common_profiles: List[str],
) -> Dict:
    """
    Compare actual cell proportions to predicted proportions.

    Returns metrics including Pearson correlation per profile and overall.
    """
    # Align indices
    common_spots = actual.index.intersection(predicted.index)

    if len(common_spots) == 0:
        return {"error": "No common spots"}

    actual_aligned = actual.loc[common_spots, common_profiles]
    predicted_aligned = predicted.loc[common_spots, common_profiles]

    results = {
        "n_spots": len(common_spots),
        "profiles_compared": common_profiles,
        "per_profile": {},
    }

    # Per-profile correlation
    overall_actual = []
    overall_predicted = []

    for profile in common_profiles:
        if profile not in actual_aligned.columns or profile not in predicted_aligned.columns:
            continue

        a = actual_aligned[profile].values
        p = predicted_aligned[profile].values

        # Skip if no variance
        if a.std() == 0 or p.std() == 0:
            results["per_profile"][profile] = {
                "pearson_r": np.nan,
                "rmse": np.nan,
            }
            continue

        r, pval = stats.pearsonr(a, p)
        rmse = np.sqrt(np.mean((a - p) ** 2))

        results["per_profile"][profile] = {
            "pearson_r": float(r),
            "pearson_p": float(pval),
            "rmse": float(rmse),
            "mean_actual": float(a.mean()),
            "mean_predicted": float(p.mean()),
        }

        overall_actual.extend(a)
        overall_predicted.extend(p)

    # Overall correlation
    if len(overall_actual) > 0:
        r, pval = stats.pearsonr(overall_actual, overall_predicted)
        results["overall_pearson_r"] = float(r)
        results["overall_rmse"] = float(np.sqrt(np.mean((np.array(overall_actual) - np.array(overall_predicted)) ** 2)))

    return results


def main(output_dir: Path):
    """Run spot-level validation analysis."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    figures_dir = output_dir / "figures"
    figures_dir.mkdir(exist_ok=True)

    logger.info("Loading Xenium single-cell data...")
    adata = load_xenium_data(XENIUM_DATA_DIR)
    adata_gex, adata_protein = split_gex_protein(adata)

    protein_matrix = adata_protein.X
    if hasattr(protein_matrix, 'toarray'):
        protein_matrix = protein_matrix.toarray()
    protein_names = list(adata_protein.var_names)
    cell_ids = np.array(adata_protein.obs_names)

    logger.info(f"Loaded {len(cell_ids):,} cells with {len(protein_names)} proteins")

    # Compute thresholds
    logger.info("Computing protein thresholds...")
    thresholds = get_protein_thresholds_gmm(protein_matrix, protein_names)

    # Classify cells by both profile sets
    logger.info("Classifying cells by autodiscovered profiles...")
    auto_positive = classify_cells_to_profiles(
        protein_matrix, protein_names, AUTODISCOVERED_PROFILES, thresholds
    )

    logger.info("Classifying cells by artificial GT profiles...")
    art_positive = classify_cells_to_profiles(
        protein_matrix, protein_names, ARTIFICIAL_GT_PROFILES, thresholds
    )

    all_results = {
        "n_cells": len(cell_ids),
        "n_proteins": len(protein_names),
        "regions": {},
    }

    # Process each region
    for region_id in range(5):
        logger.info(f"\n{'='*60}")
        logger.info(f"Processing Region {region_id}")
        logger.info("="*60)

        try:
            # Load cell-to-spot mapping
            cell_to_spot = load_cell_to_spot_mapping(region_id)
            logger.info(f"  {len(cell_to_spot):,} cells mapped to spots")

            # Compute actual proportions from autodiscovered profiles
            actual_auto = compute_spot_proportions_from_cells(
                cell_to_spot, auto_positive, cell_ids
            )
            logger.info(f"  Computed autodiscovered proportions for {len(actual_auto)} spots")

            # Compute actual proportions from artificial profiles
            actual_art = compute_spot_proportions_from_cells(
                cell_to_spot, art_positive, cell_ids
            )

            # Try to load CITEgeist predictions
            try:
                predicted = load_citegeist_predictions(region_id)
                logger.info(f"  Loaded CITEgeist predictions: {predicted.shape}")

                # Find common profiles between actual and predicted
                common_auto = [p for p in AUTODISCOVERED_PROFILES.keys()
                              if p in actual_auto.columns and p in predicted.columns]

                if len(common_auto) > 0:
                    comparison = compare_proportions(actual_auto, predicted, common_auto)
                    logger.info(f"  Overall Pearson r (autodiscovered): {comparison.get('overall_pearson_r', 'N/A'):.3f}")

                    for profile, metrics in comparison.get("per_profile", {}).items():
                        if "pearson_r" in metrics and not np.isnan(metrics["pearson_r"]):
                            logger.info(f"    {profile}: r={metrics['pearson_r']:.3f}")
                else:
                    comparison = {"note": "No common profiles found"}

            except FileNotFoundError as e:
                logger.warning(f"  Could not load CITEgeist predictions: {e}")
                comparison = {"error": str(e)}

            # Save region results
            region_results = {
                "n_cells_mapped": len(cell_to_spot),
                "n_spots": len(actual_auto),
                "autodiscovered_comparison": comparison,
            }

            # Save actual proportions
            actual_auto.to_csv(output_dir / f"region_{region_id}_actual_autodiscovered.csv")
            actual_art.to_csv(output_dir / f"region_{region_id}_actual_artificial.csv")

            all_results["regions"][str(region_id)] = region_results

        except Exception as e:
            logger.error(f"  Error processing region {region_id}: {e}")
            all_results["regions"][str(region_id)] = {"error": str(e)}

    # Summary statistics
    logger.info("\n" + "="*60)
    logger.info("SPOT-LEVEL VALIDATION SUMMARY")
    logger.info("="*60)

    # Aggregate correlations across regions
    all_r_values = []
    for region_id, results in all_results["regions"].items():
        if "autodiscovered_comparison" in results:
            comp = results["autodiscovered_comparison"]
            if "overall_pearson_r" in comp:
                all_r_values.append(comp["overall_pearson_r"])
                logger.info(f"Region {region_id}: r = {comp['overall_pearson_r']:.3f}")

    if all_r_values:
        all_results["mean_overall_r"] = float(np.mean(all_r_values))
        all_results["std_overall_r"] = float(np.std(all_r_values))
        logger.info(f"\nMean r across regions: {all_results['mean_overall_r']:.3f} +/- {all_results['std_overall_r']:.3f}")

    # Save results
    with open(output_dir / "spot_level_validation_results.json", "w") as f:
        json.dump(all_results, f, indent=2, default=str)

    logger.info(f"\nResults saved to: {output_dir}")

    return all_results


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Spot-level validation analysis")
    parser.add_argument(
        "--output_dir",
        type=str,
        default="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/evaluation/spot_level_validation",
        help="Output directory",
    )
    parser.add_argument("--verbose", action="store_true")

    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
    )
    logger.setLevel(logging.INFO)

    main(args.output_dir)
