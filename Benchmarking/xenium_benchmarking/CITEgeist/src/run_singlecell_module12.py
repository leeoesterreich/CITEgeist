#!/usr/bin/env python
"""
Run Modules 1 and 2 on Xenium single-cell data.

Executes marker interest detection (Module 1) and profile discovery (Module 2)
at single-cell resolution. Supports full dataset or quadrant subsets.

Usage:
    python run_singlecell_module12.py --mode full
    python run_singlecell_module12.py --mode quadrant --quadrant-id 0
"""

import argparse
import json
import logging
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import scanpy as sc

# Add paths
REPO_ROOT = Path(__file__).parent.parent.parent.parent.parent
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(REPO_ROOT / "Benchmarking" / "xenium_benchmarking" / "CITEgeist" / "src"))

from CITEgeist.model.marker_interest import identify_interesting_markers
from CITEgeist.model.spatial_colocalization import (
    analyze_marker_colocalization,
    discover_profiles_continuous,
    select_profiles,
)
from load_xenium_singlecell import load_xenium_singlecell

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)

# Output directory
OUTPUT_BASE = REPO_ROOT / "Benchmarking" / "xenium_benchmarking" / "CITEgeist" / "output_singlecell_demonstration"


def get_quadrant_bounds(coords: np.ndarray, quadrant_id: int) -> Tuple[float, float, float, float]:
    """
    Get bounds for a spatial quadrant.

    Quadrants:
        Q0: bottom-left, Q1: bottom-right, Q2: top-left, Q3: top-right
    """
    x_mid = (coords[:, 0].min() + coords[:, 0].max()) / 2
    y_mid = (coords[:, 1].min() + coords[:, 1].max()) / 2
    x_min, x_max = coords[:, 0].min(), coords[:, 0].max()
    y_min, y_max = coords[:, 1].min(), coords[:, 1].max()

    bounds = {
        0: (x_min, x_mid, y_min, y_mid),  # bottom-left
        1: (x_mid, x_max, y_min, y_mid),  # bottom-right
        2: (x_min, x_mid, y_mid, y_max),  # top-left
        3: (x_mid, x_max, y_mid, y_max),  # top-right
    }
    return bounds[quadrant_id]


def load_data(
    mode: str,
    quadrant_id: Optional[int] = None,
) -> Tuple[sc.AnnData, sc.AnnData, str]:
    """
    Load Xenium data for full or quadrant mode.

    Returns:
        Tuple of (adata_gex, adata_protein, output_subdir)
    """
    if mode == "full":
        logger.info("Loading FULL dataset (all cells)")
        adata_gex, adata_protein = load_xenium_singlecell()
        output_subdir = "full"
    else:
        logger.info(f"Loading QUADRANT {quadrant_id}")
        # First load full to get bounds
        adata_full, _ = load_xenium_singlecell(max_cells=1000)  # Quick load for bounds
        full_coords = adata_full.obsm["spatial"]
        bounds = get_quadrant_bounds(full_coords, quadrant_id)

        # Now load with bounds
        adata_gex, adata_protein = load_xenium_singlecell(region_bounds=bounds)
        output_subdir = f"quadrants/Q{quadrant_id}"

    logger.info(f"Loaded {adata_gex.shape[0]:,} cells × {adata_gex.shape[1]} genes")
    logger.info(f"Proteins: {list(adata_protein.var_names)}")

    return adata_gex, adata_protein, output_subdir


def run_module1(
    X_protein: np.ndarray,
    coords: np.ndarray,
    marker_names: List[str],
    output_dir: Path,
) -> Tuple[Dict, Any]:
    """Run Module 1: Marker Interest Detection."""
    logger.info("=" * 70)
    logger.info("MODULE 1: Marker Interest Detection")
    logger.info("=" * 70)

    # At single-cell resolution, use more neighbors for Moran's I
    result = identify_interesting_markers(
        X=X_protein,
        coords=coords,
        marker_names=marker_names,
        kurtosis_threshold=2.0,
        morans_threshold=0.1,
        morans_k=15,  # More neighbors at single-cell resolution
        morans_n_perm=199,
        verbose=True,
    )

    # Build output
    output = {
        "n_markers_total": len(marker_names),
        "n_interesting": len(result.interesting_markers),
        "n_boring": len(result.boring_markers),
        "interesting_markers": result.interesting_markers,
        "boring_markers": result.boring_markers,
        "kurtosis_threshold": float(result.kurtosis_threshold),
        "morans_threshold": float(result.morans_threshold),
        "marker_details": [
            {
                "marker": m.name,
                "interest_score": float(m.interest_score),
                "kurtosis": float(m.kurtosis),
                "gmm_snr": float(m.gmm_snr),
                "morans_i": float(m.morans_i),
                "morans_i_pvalue": float(m.morans_i_pvalue),
                "passed_kurtosis": bool(m.passed_kurtosis),
                "passed_gmm": bool(m.passed_gmm),
                "passed_morans": bool(m.passed_morans),
            }
            for m in result.markers
        ],
    }

    # Save outputs
    with open(output_dir / "module1_marker_interest.json", "w") as f:
        json.dump(output, f, indent=2)

    df = result.to_dataframe()
    df.to_csv(output_dir / "module1_marker_interest.csv", index=False)

    logger.info(f"Interesting markers: {len(result.interesting_markers)}/{len(marker_names)}")
    logger.info(f"Markers: {result.interesting_markers}")

    return output, result


def run_module2(
    X_protein: np.ndarray,
    coords: np.ndarray,
    marker_names: List[str],
    interesting_markers: List[str],
    output_dir: Path,
) -> Dict:
    """Run Module 2: Profile Discovery."""
    logger.info("=" * 70)
    logger.info("MODULE 2: Profile Discovery")
    logger.info("=" * 70)

    # Filter to interesting markers
    marker_indices = [marker_names.index(m) for m in interesting_markers]
    X_filtered = X_protein[:, marker_indices]

    # Module 2a: Colocalization analysis
    logger.info("Module 2a: Colocalization Analysis")
    coloc_result = analyze_marker_colocalization(
        X=X_filtered,
        coords=coords,
        marker_names=interesting_markers,
        neighbor_k=15,  # More neighbors at single-cell resolution
        n_permutations=199,
    )

    # Save colocalization results
    coloc_df = coloc_result.to_dataframe()
    coloc_df.to_csv(output_dir / "module2a_colocalization.csv", index=False)

    # Module 2b: Profile discovery
    logger.info("Module 2b: Profile Discovery")
    profile_result = discover_profiles_continuous(
        colocalization_result=coloc_result,
    )

    # Save raw profiles (profiles is a List[List[str]])
    raw_profiles = {
        f"Profile_{i}": {"markers": list(markers)}
        for i, markers in enumerate(profile_result.profiles)
    }
    with open(output_dir / "module2b_profiles_raw.json", "w") as f:
        json.dump(raw_profiles, f, indent=2)

    # Module 2c: Profile selection
    logger.info("Module 2c: Profile Selection")
    # profiles is already a list of lists
    profiles_list = [list(markers) for markers in profile_result.profiles]
    selection_result = select_profiles(
        X=X_filtered,
        coords=coords,
        marker_names=interesting_markers,
        profiles=profiles_list,
        interesting_markers=interesting_markers,
        colocalization_result=coloc_result,
        max_profiles=10,
    )

    # Save selected profiles (convert list of lists to named dict)
    selected_profiles = selection_result.selected_profiles
    selected_output = {
        f"Profile_{i}": {"markers": list(markers)}
        for i, markers in enumerate(selected_profiles)
    }
    with open(output_dir / "module2c_profiles_selected.json", "w") as f:
        json.dump(selected_output, f, indent=2)

    logger.info(f"Discovered {len(selected_profiles)} profiles:")
    for i, markers in enumerate(selected_profiles):
        logger.info(f"  Profile_{i}: {list(markers)}")

    return {
        "n_profiles_raw": len(profile_result.profiles),
        "n_profiles_selected": len(selected_profiles),
        "profiles": selected_output,
    }


def main():
    parser = argparse.ArgumentParser(description="Run Modules 1-2 on single-cell data")
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
    args = parser.parse_args()

    if args.mode == "quadrant" and args.quadrant_id is None:
        parser.error("--quadrant-id required when mode=quadrant")

    # Load data
    adata_gex, adata_protein, output_subdir = load_data(args.mode, args.quadrant_id)

    # Setup output directory
    output_dir = OUTPUT_BASE / output_subdir
    output_dir.mkdir(parents=True, exist_ok=True)

    # Extract arrays
    X_protein = np.asarray(adata_protein.X.todense() if hasattr(adata_protein.X, "todense") else adata_protein.X)
    coords = adata_protein.obsm["spatial"]
    marker_names = list(adata_protein.var_names)

    # Save data summary
    summary = {
        "mode": args.mode,
        "quadrant_id": args.quadrant_id,
        "n_cells": int(adata_gex.shape[0]),
        "n_genes": int(adata_gex.shape[1]),
        "n_proteins": int(adata_protein.shape[1]),
        "protein_names": marker_names,
        "spatial_extent": {
            "x_min": float(coords[:, 0].min()),
            "x_max": float(coords[:, 0].max()),
            "y_min": float(coords[:, 1].min()),
            "y_max": float(coords[:, 1].max()),
        },
    }
    with open(output_dir / "data_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    # Run Module 1
    module1_output, module1_result = run_module1(
        X_protein=X_protein,
        coords=coords,
        marker_names=marker_names,
        output_dir=output_dir,
    )

    # Run Module 2
    module2_output = run_module2(
        X_protein=X_protein,
        coords=coords,
        marker_names=marker_names,
        interesting_markers=module1_result.interesting_markers,
        output_dir=output_dir,
    )

    # Save combined summary
    combined_summary = {
        "data": summary,
        "module1": module1_output,
        "module2": module2_output,
    }
    with open(output_dir / "module12_summary.json", "w") as f:
        json.dump(combined_summary, f, indent=2)

    logger.info("=" * 70)
    logger.info("MODULE 1-2 COMPLETE")
    logger.info(f"Output saved to: {output_dir}")
    logger.info("=" * 70)


if __name__ == "__main__":
    main()
