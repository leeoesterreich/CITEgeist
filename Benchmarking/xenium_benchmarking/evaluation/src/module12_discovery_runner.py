"""
Module 1-2 discovery runner for comparison experiment.

Runs CITEgeist Module 1 (marker interest) and Module 2 (profile discovery)
on Xenium protein data at spot or single-cell resolution. Optionally sweeps
top_k parameter for sensitivity analysis.

Usage:
    python module12_discovery_runner.py \
        --region 0 \
        --resolution-level spot \
        --output-dir results/discovery_comparison

    # With top_k sweep:
    python module12_discovery_runner.py \
        --region 0 \
        --resolution-level spot \
        --top-k-sweep \
        --output-dir results/discovery_comparison
"""
import argparse
import json
import logging
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(REPO_ROOT / "Benchmarking" / "xenium_benchmarking" / "CITEgeist" / "src"))

from CITEgeist.model.marker_interest import identify_interesting_markers
from CITEgeist.model.spatial_colocalization import (
    analyze_marker_colocalization,
    discover_profiles_continuous,
)
from leiden_baseline_comparison import load_data

logger = logging.getLogger(__name__)
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)


def run_module12(
    X: np.ndarray,
    coords: np.ndarray,
    marker_names: List[str],
    top_k: int = 3,
    neighbor_k: int = 15,
    seed: int = 42,
) -> Dict:
    """
    Run Module 1 and Module 2 on protein expression data.

    Returns:
        Dict with module1 results, module2 profiles, and parameters.
    """
    # Module 1: Marker interest detection
    logger.info("Running Module 1: Marker interest detection...")
    m1_result = identify_interesting_markers(
        X=X,
        coords=coords,
        marker_names=marker_names,
        morans_k=min(neighbor_k, 20),
        smooth_k=min(neighbor_k, 10),
        morans_n_perm=99,
        seed=seed,
    )

    interesting = m1_result.interesting_markers
    logger.info(f"Module 1: {len(interesting)}/{len(marker_names)} interesting markers")

    if len(interesting) < 2:
        logger.warning("Fewer than 2 interesting markers, returning empty profiles")
        return {
            "interesting_markers": interesting,
            "profiles": [],
            "singletons": [],
            "top_k": top_k,
        }

    # Module 2a: Colocalization analysis
    logger.info("Running Module 2a: Colocalization analysis...")
    coloc_result = analyze_marker_colocalization(
        X=X,
        coords=coords,
        marker_names=marker_names,
        markers_to_analyze=interesting,
        neighbor_k=neighbor_k,
        multi_scale_k=[neighbor_k // 2, neighbor_k, neighbor_k * 2],
        n_permutations=999,
        seed=seed,
    )
    logger.info(f"Module 2a: {len(coloc_result.pairs)} marker pairs analyzed")

    # Module 2b: Profile discovery (continuous)
    logger.info(f"Running Module 2b: Profile discovery (top_k={top_k})...")
    profile_result = discover_profiles_continuous(
        colocalization_result=coloc_result,
        top_k=top_k,
        distance_metric="colocalization_score",
        seed=seed,
    )

    profiles = [list(p) for p in profile_result.profiles]
    singletons = list(profile_result.singletons)

    logger.info(f"Module 2b: {len(profiles)} profiles, {len(singletons)} singletons")
    for i, p in enumerate(profiles):
        logger.info(f"  Profile {i}: {p}")
    if singletons:
        logger.info(f"  Singletons: {singletons}")

    return {
        "interesting_markers": interesting,
        "profiles": profiles,
        "singletons": singletons,
        "top_k": top_k,
    }


def main():
    parser = argparse.ArgumentParser(
        description="Run Module 1-2 for discovery comparison"
    )
    parser.add_argument("--region", type=int, required=True, help="Region ID (0-4)")
    parser.add_argument(
        "--resolution-level",
        choices=["spot", "cell"],
        required=True,
    )
    parser.add_argument("--top-k", type=int, default=3)
    parser.add_argument(
        "--top-k-sweep",
        action="store_true",
        help="Run top_k sensitivity sweep (2, 3, 4, 5)",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default=str(
            REPO_ROOT / "Benchmarking/xenium_benchmarking/evaluation/results/discovery_comparison"
        ),
    )
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    X, coords, marker_names = load_data(args.region, args.resolution_level)

    if args.top_k_sweep:
        sweep_results = {}
        for tk in [2, 3, 4, 5]:
            logger.info(f"\n{'='*60}\ntop_k = {tk}\n{'='*60}")
            result = run_module12(X, coords, marker_names, top_k=tk)
            sweep_results[str(tk)] = result

        output_file = (
            output_dir
            / f"module12_topk_sweep_region_{args.region}_{args.resolution_level}.json"
        )
        with open(output_file, "w") as f:
            json.dump(
                {
                    "region_id": args.region,
                    "resolution_level": args.resolution_level,
                    "n_observations": X.shape[0],
                    "n_markers": X.shape[1],
                    "marker_names": marker_names,
                    "sweep_results": sweep_results,
                },
                f,
                indent=2,
            )
        logger.info(f"Saved top_k sweep to {output_file}")
    else:
        result = run_module12(X, coords, marker_names, top_k=args.top_k)

        output_file = (
            output_dir
            / f"module12_region_{args.region}_{args.resolution_level}.json"
        )
        with open(output_file, "w") as f:
            json.dump(
                {
                    "region_id": args.region,
                    "resolution_level": args.resolution_level,
                    "n_observations": X.shape[0],
                    "n_markers": X.shape[1],
                    "marker_names": marker_names,
                    "result": result,
                },
                f,
                indent=2,
            )
        logger.info(f"Saved Module 1-2 results to {output_file}")


if __name__ == "__main__":
    main()
