#!/usr/bin/env python
"""
Radius sweep benchmark for CITEgeist Module 3.

Tests 0, 1, 2, 3 ring neighborhoods to determine optimal spatial context
for cell proportion deconvolution on Xenium pseudo-Visium data.

Radii calibrated for 100µm spot spacing:
- 50:  0 rings (center only, no neighbors)
- 105: 1 ring  (~6 neighbors)
- 205: 2 rings (~18 neighbors)
- 305: 3 rings (~36 neighbors)

Usage:
    python run_radius_sweep.py --region-id 0
    python run_radius_sweep.py --region-id 0 --run-gex  # Include GEX deconvolution
"""

import argparse
import json
import logging
import sys
import time
from pathlib import Path
from typing import Any, Dict, List

import scanpy as sc

# Add paths
REPO_ROOT = Path(__file__).parent.parent.parent.parent.parent
sys.path.insert(0, str(REPO_ROOT))

from Benchmarking.xenium_benchmarking.CITEgeist.src.run_benchmark import run_manual_benchmark

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)

# Radius sweep configuration
RADII = [50, 105, 205, 305]
RING_NAMES = {
    50: "0_rings",
    105: "1_ring",
    205: "2_rings",
    305: "3_rings",
}


def run_radius_sweep(
    region_id: int,
    input_dir: Path,
    output_base_dir: Path,
    lambda_reg: float = 1.0,
    alpha_elastic: float = 0.7,
    max_y_change: float = 0.4,
    min_counts: int = 25,
    run_gex: bool = False,
) -> Dict[str, Any]:
    """
    Run benchmark across all radii for a single region.

    Args:
        region_id: Region ID (0-4)
        input_dir: Directory containing h5ad_objects/
        output_base_dir: Base output directory for sweep results
        lambda_reg: Regularization lambda
        alpha_elastic: Elastic net alpha
        max_y_change: Max Y change constraint
        min_counts: Min counts filter
        run_gex: Whether to run GEX deconvolution

    Returns:
        Dict with results for each radius
    """
    # Load data once
    gex_path = input_dir / "h5ad_objects" / f"Xenium_region_{region_id}_GEX.h5ad"
    protein_path = input_dir / "h5ad_objects" / f"Xenium_region_{region_id}_CITE.h5ad"

    for path in [gex_path, protein_path]:
        if not path.exists():
            logger.error(f"File not found: {path}")
            sys.exit(1)

    logger.info(f"Loading data for region {region_id}...")
    gex_adata = sc.read_h5ad(gex_path)
    protein_adata = sc.read_h5ad(protein_path)

    sweep_results = {
        "region_id": region_id,
        "radii_tested": RADII,
        "results_by_radius": {},
    }

    for radius in RADII:
        ring_name = RING_NAMES[radius]
        logger.info("=" * 70)
        logger.info(f"RADIUS SWEEP: region={region_id}, radius={radius} ({ring_name})")
        logger.info("=" * 70)

        output_dir = output_base_dir / f"radius_{radius}"
        output_dir.mkdir(parents=True, exist_ok=True)

        start_time = time.time()

        try:
            # Need fresh copies since preprocessing modifies the data
            gex_copy = gex_adata.copy()
            protein_copy = protein_adata.copy()

            result = run_manual_benchmark(
                region_id=region_id,
                gex_adata=gex_copy,
                protein_adata=protein_copy,
                output_dir=output_dir,
                radius=radius,
                lambda_reg=lambda_reg,
                alpha_elastic=alpha_elastic,
                max_y_change=max_y_change,
                min_counts=min_counts,
                run_gex=run_gex,
                gex_radius=radius,  # Use same radius for GEX
            )

            result["ring_name"] = ring_name
            result["sweep_runtime_sec"] = time.time() - start_time
            sweep_results["results_by_radius"][radius] = result

            logger.info(f"Completed radius={radius} in {result['sweep_runtime_sec']:.1f}s")

        except Exception as e:
            logger.error(f"Failed for radius={radius}: {e}")
            sweep_results["results_by_radius"][radius] = {
                "error": str(e),
                "ring_name": ring_name,
            }

    # Save sweep summary
    summary_path = output_base_dir / f"sweep_summary_region_{region_id}.json"
    with open(summary_path, "w") as f:
        json.dump(sweep_results, f, indent=2, default=str)

    logger.info(f"Sweep summary saved to {summary_path}")
    return sweep_results


def main():
    parser = argparse.ArgumentParser(
        description="Run radius sweep benchmark for CITEgeist Module 3"
    )

    parser.add_argument(
        "--region-id",
        type=int,
        required=True,
        help="Region ID to process (0-4)",
    )
    parser.add_argument(
        "--input-dir",
        type=str,
        default=str(REPO_ROOT / "Benchmarking/xenium_pseudovisium/data_protein_gt"),
        help="Input directory with h5ad_objects/",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default=str(REPO_ROOT / "Benchmarking/xenium_benchmarking/CITEgeist/output/radius_sweep"),
        help="Output directory for sweep results",
    )

    # Optimization parameters (same defaults as run_benchmark.py)
    parser.add_argument("--lambda-reg", type=float, default=1.0)
    parser.add_argument("--alpha-elastic", type=float, default=0.7)
    parser.add_argument("--max-y-change", type=float, default=0.4)
    parser.add_argument("--min-counts", type=int, default=25)
    parser.add_argument("--run-gex", action="store_true", help="Run GEX deconvolution")

    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)

    results = run_radius_sweep(
        region_id=args.region_id,
        input_dir=input_dir,
        output_base_dir=output_dir,
        lambda_reg=args.lambda_reg,
        alpha_elastic=args.alpha_elastic,
        max_y_change=args.max_y_change,
        min_counts=args.min_counts,
        run_gex=args.run_gex,
    )

    print(f"\nCompleted radius sweep for region {args.region_id}")
    print(f"Radii tested: {RADII}")
    for radius, result in results["results_by_radius"].items():
        if "error" in result:
            print(f"  radius={radius}: FAILED - {result['error']}")
        else:
            print(f"  radius={radius} ({result['ring_name']}): {result['runtime_proportions_sec']:.1f}s")


if __name__ == "__main__":
    main()
