"""
Run Module 1 (Marker Interest Detection) on Xenium single-cell data.

This script demonstrates that CITEgeist's marker interest scoring works
at single-cell resolution, not just spot-level. It compares results
from single-cell and pseudo-Visium analyses to validate resolution
independence.

Key demonstration:
- Spatial autocorrelation (Moran's I) is computed per cell
- Kurtosis and GMM signal detection work identically
- Same markers should be identified as "interesting" at both resolutions
"""

import json
import logging
import sys
from pathlib import Path
from typing import Tuple, Dict, Optional
import argparse

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import spearmanr, pearsonr

# Add CITEgeist to path
REPO_ROOT = Path(__file__).parent.parent.parent.parent.parent
sys.path.insert(0, str(REPO_ROOT / "CITEgeist"))
sys.path.insert(0, str(REPO_ROOT / "Benchmarking" / "xenium_pseudovisium" / "src"))

from model.marker_interest import identify_interesting_markers, MarkerInterestResult
from load_xenium_singlecell import load_xenium_singlecell, XENIUM_DATA_DIR

logger = logging.getLogger(__name__)

# Output directory
OUTPUT_DIR = REPO_ROOT / "Benchmarking" / "xenium_benchmarking" / "CITEgeist" / "output_singlecell"


def run_module1_singlecell(
    region_id: int = 0,
    max_cells: Optional[int] = None,
    morans_k: int = 20,
    smooth_k: int = 10,
    morans_n_perm: int = 99,
    seed: int = 42,
) -> Tuple[MarkerInterestResult, sc.AnnData]:
    """
    Run Module 1 on Xenium single-cell data.

    Args:
        region_id: Which region to analyze (0-4).
        max_cells: Maximum cells to use (for faster testing).
        morans_k: Number of neighbors for Moran's I (increased for denser data).
        smooth_k: Number of neighbors for spatial smoothing.
        morans_n_perm: Permutations for Moran's I null distribution.
        seed: Random seed.

    Returns:
        Tuple of (MarkerInterestResult, adata_protein).
    """
    logger.info(f"Loading single-cell data for region {region_id}")
    adata_gex, adata_protein = load_xenium_singlecell(
        region_id=region_id,
        max_cells=max_cells,
        seed=seed,
    )

    n_cells = adata_protein.shape[0]
    n_proteins = adata_protein.shape[1]
    logger.info(f"Loaded {n_cells:,} cells × {n_proteins} proteins")

    # Get protein expression matrix
    X = adata_protein.X
    if hasattr(X, "toarray"):
        X = X.toarray()
    X = np.asarray(X, dtype=np.float64)

    # Get spatial coordinates
    coords = adata_protein.obsm["spatial"]
    marker_names = list(adata_protein.var_names)

    logger.info("Running Module 1 (marker interest detection) on single cells...")
    logger.info(f"  morans_k={morans_k}, smooth_k={smooth_k}, n_perm={morans_n_perm}")

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

    return result, adata_protein


def run_module1_spots(
    region_id: int = 0,
    morans_k: int = 8,
    smooth_k: int = 6,
    morans_n_perm: int = 199,
    seed: int = 42,
) -> Tuple[MarkerInterestResult, sc.AnnData]:
    """
    Run Module 1 on pseudo-Visium spots for comparison.

    Args:
        region_id: Which region to analyze (0-4).
        morans_k: Number of neighbors for Moran's I.
        smooth_k: Number of neighbors for spatial smoothing.
        morans_n_perm: Permutations for Moran's I null distribution.
        seed: Random seed.

    Returns:
        Tuple of (MarkerInterestResult, adata_protein).
    """
    # Load pseudo-Visium spot data
    spot_data_dir = REPO_ROOT / "Benchmarking" / "xenium_pseudovisium" / "data_granular_gt" / "h5ad_objects"
    protein_path = spot_data_dir / f"Xenium_region_{region_id}_CITE.h5ad"

    logger.info(f"Loading pseudo-Visium spot data from {protein_path}")
    adata_protein = sc.read_h5ad(protein_path)

    n_spots = adata_protein.shape[0]
    n_proteins = adata_protein.shape[1]
    logger.info(f"Loaded {n_spots:,} spots × {n_proteins} proteins")

    # Get protein expression matrix
    X = adata_protein.X
    if hasattr(X, "toarray"):
        X = X.toarray()
    X = np.asarray(X, dtype=np.float64)

    # Get spatial coordinates
    coords = adata_protein.obsm["spatial"]
    marker_names = list(adata_protein.var_names)

    logger.info("Running Module 1 (marker interest detection) on spots...")
    logger.info(f"  morans_k={morans_k}, smooth_k={smooth_k}, n_perm={morans_n_perm}")

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

    return result, adata_protein


def compare_results(
    cell_result: MarkerInterestResult,
    spot_result: MarkerInterestResult,
) -> Dict:
    """
    Compare Module 1 results from single-cell and spot analyses.

    Args:
        cell_result: Results from single-cell analysis.
        spot_result: Results from spot analysis.

    Returns:
        Dict with comparison metrics.
    """
    cell_df = cell_result.to_dataframe()
    spot_df = spot_result.to_dataframe()

    # Merge on marker name
    merged = cell_df.merge(
        spot_df,
        on="marker",
        suffixes=("_cell", "_spot"),
    )

    # Correlation of interest scores
    r_interest, p_interest = spearmanr(
        merged["interest_score_cell"],
        merged["interest_score_spot"],
    )

    # Correlation of Moran's I
    r_morans, p_morans = spearmanr(
        merged["morans_i_cell"],
        merged["morans_i_spot"],
    )

    # Correlation of kurtosis
    r_kurtosis, p_kurtosis = spearmanr(
        merged["kurtosis_cell"],
        merged["kurtosis_spot"],
    )

    # Jaccard similarity of top-10 markers
    cell_top10 = set(cell_df.head(10)["marker"])
    spot_top10 = set(spot_df.head(10)["marker"])
    jaccard_top10 = len(cell_top10 & spot_top10) / len(cell_top10 | spot_top10)

    # Jaccard of interesting markers
    cell_interesting = set(cell_result.interesting_markers)
    spot_interesting = set(spot_result.interesting_markers)
    jaccard_interesting = len(cell_interesting & spot_interesting) / max(
        len(cell_interesting | spot_interesting), 1
    )

    comparison = {
        "spearman_interest_score": {"r": r_interest, "p": p_interest},
        "spearman_morans_i": {"r": r_morans, "p": p_morans},
        "spearman_kurtosis": {"r": r_kurtosis, "p": p_kurtosis},
        "jaccard_top10": jaccard_top10,
        "jaccard_interesting": jaccard_interesting,
        "n_interesting_cell": len(cell_interesting),
        "n_interesting_spot": len(spot_interesting),
        "interesting_both": list(cell_interesting & spot_interesting),
        "interesting_cell_only": list(cell_interesting - spot_interesting),
        "interesting_spot_only": list(spot_interesting - cell_interesting),
    }

    return comparison


def main():
    parser = argparse.ArgumentParser(
        description="Run Module 1 on Xenium single-cell data and compare to spot-level"
    )
    parser.add_argument(
        "--region", type=int, default=0, help="Region ID (0-4)"
    )
    parser.add_argument(
        "--max-cells", type=int, default=None,
        help="Max cells to use (None for all)"
    )
    parser.add_argument(
        "--morans-k-cell", type=int, default=20,
        help="Moran's I neighbors for cells (default: 20)"
    )
    parser.add_argument(
        "--smooth-k-cell", type=int, default=10,
        help="Smoothing neighbors for cells (default: 10)"
    )
    parser.add_argument(
        "--morans-k-spot", type=int, default=8,
        help="Moran's I neighbors for spots (default: 8)"
    )
    parser.add_argument(
        "--smooth-k-spot", type=int, default=6,
        help="Smoothing neighbors for spots (default: 6)"
    )
    parser.add_argument(
        "--n-perm", type=int, default=99,
        help="Permutations for Moran's I (default: 99)"
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

    # Run Module 1 on single cells
    logger.info("=" * 60)
    logger.info("SINGLE-CELL ANALYSIS")
    logger.info("=" * 60)
    cell_result, adata_cells = run_module1_singlecell(
        region_id=args.region,
        max_cells=args.max_cells,
        morans_k=args.morans_k_cell,
        smooth_k=args.smooth_k_cell,
        morans_n_perm=args.n_perm,
        seed=args.seed,
    )

    # Run Module 1 on spots
    logger.info("\n" + "=" * 60)
    logger.info("SPOT-LEVEL ANALYSIS")
    logger.info("=" * 60)
    spot_result, adata_spots = run_module1_spots(
        region_id=args.region,
        morans_k=args.morans_k_spot,
        smooth_k=args.smooth_k_spot,
        morans_n_perm=args.n_perm,
        seed=args.seed,
    )

    # Compare results
    logger.info("\n" + "=" * 60)
    logger.info("RESOLUTION COMPARISON")
    logger.info("=" * 60)
    comparison = compare_results(cell_result, spot_result)

    # Print summary
    print("\n" + "=" * 60)
    print("RESOLUTION INDEPENDENCE VALIDATION")
    print("=" * 60)
    print(f"\nCorrelations between cell and spot level:")
    print(f"  Interest score: r={comparison['spearman_interest_score']['r']:.3f} "
          f"(p={comparison['spearman_interest_score']['p']:.2e})")
    print(f"  Moran's I:      r={comparison['spearman_morans_i']['r']:.3f} "
          f"(p={comparison['spearman_morans_i']['p']:.2e})")
    print(f"  Kurtosis:       r={comparison['spearman_kurtosis']['r']:.3f} "
          f"(p={comparison['spearman_kurtosis']['p']:.2e})")

    print(f"\nJaccard similarity:")
    print(f"  Top-10 markers:      {comparison['jaccard_top10']:.3f}")
    print(f"  Interesting markers: {comparison['jaccard_interesting']:.3f}")

    print(f"\nInteresting markers:")
    print(f"  Single-cell: {comparison['n_interesting_cell']} markers")
    print(f"  Spot-level:  {comparison['n_interesting_spot']} markers")
    print(f"  Both:        {comparison['interesting_both']}")
    print(f"  Cell only:   {comparison['interesting_cell_only']}")
    print(f"  Spot only:   {comparison['interesting_spot_only']}")

    # Validation thresholds
    passed_interest = comparison["spearman_interest_score"]["r"] > 0.7
    passed_morans = comparison["spearman_morans_i"]["r"] > 0.7
    passed_jaccard = comparison["jaccard_top10"] >= 0.5

    print("\n" + "=" * 60)
    print("VALIDATION STATUS")
    print("=" * 60)
    print(f"  Interest score correlation > 0.7: {'PASS' if passed_interest else 'FAIL'}")
    print(f"  Moran's I correlation > 0.7:      {'PASS' if passed_morans else 'FAIL'}")
    print(f"  Top-10 Jaccard >= 0.5:            {'PASS' if passed_jaccard else 'FAIL'}")

    overall = passed_interest and passed_morans and passed_jaccard
    print(f"\n  OVERALL: {'PASS - Resolution independent!' if overall else 'NEEDS INVESTIGATION'}")

    # Save results
    output_prefix = f"region_{args.region}"
    if args.max_cells:
        output_prefix += f"_maxcells_{args.max_cells}"

    cell_df = cell_result.to_dataframe()
    spot_df = spot_result.to_dataframe()
    cell_df.to_csv(OUTPUT_DIR / f"{output_prefix}_module1_singlecell.csv", index=False)
    spot_df.to_csv(OUTPUT_DIR / f"{output_prefix}_module1_spots.csv", index=False)

    with open(OUTPUT_DIR / f"{output_prefix}_module1_comparison.json", "w") as f:
        # Convert numpy types for JSON serialization
        comparison_json = {
            k: (
                {kk: float(vv) for kk, vv in v.items()} if isinstance(v, dict)
                else float(v) if isinstance(v, (np.floating, float))
                else int(v) if isinstance(v, (np.integer, int))
                else v
            )
            for k, v in comparison.items()
        }
        json.dump(comparison_json, f, indent=2)

    logger.info(f"\nResults saved to {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
