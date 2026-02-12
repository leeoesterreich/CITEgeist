#!/usr/bin/env python
"""
Module 3 Discrete Cell Assignment.

Alternative to continuous proportion estimation (run_module3_unified.py).
Uses Integer Quadratic Programming (IQP) to assign discrete cell counts
when nuclei segmentation data is available.

Advantages over continuous proportions:
- More interpretable: integer cell counts vs fractional proportions
- Better for downstream analyses requiring cell counts
- +10% Pearson correlation on proportion estimation (benchmarked)
- +10-16% Pearson correlation on GEX deconvolution (benchmarked)

Requirements:
- Nuclei counts per spot from Cellpose segmentation or Xenium cell mapping
- Gurobi license (academic free)

Usage:
    # With nuclei counts from Cellpose segmentation
    python run_module3_discrete.py --sample HCC22-088-P1-S1 \\
        --nuclei-file segmentation/nuclei_counts.csv \\
        --output-dir output/module3_discrete

    # For Xenium data with n_cells column in adata.obs
    python run_module3_discrete.py --sample xenium_region0 \\
        --use-obs-nuclei \\
        --output-dir output/module3_discrete
"""
import argparse
import json
import logging
import os
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq

# Add the parent directory to the system path
sys.path.insert(0, str(Path(__file__).parent.parent))

from CITEgeist.model.citegeist_model import CitegeistModel

logger = logging.getLogger(__name__)
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)

# Data location
DATA_ROOT = Path("/ix1/alee/LO_LAB/General/Lab_Data/20250210_CITEGeistPublicData_GEO_Alex/processed_files")

# Example unified profile (same as continuous method)
# Note: markers have -1 suffix as per real patient data convention
UNIFIED_PROFILE = {
    "Endothelial": {"Major": ["PECAM1-1"]},
    "Fibroblasts": {"Major": ["ACTA2-1"]},
    "B_Cells": {"Major": ["CD19-1"]},
    "Macrophages": {"Major": ["CD68-1", "CD163-1"]},
    "Monocytes": {"Major": ["CD14-1"]},
    "CD8_T_Cells": {"Major": ["CD3E-1", "CD8A-1"]},
    "CD4_T_Cells": {"Major": ["CD4-1"]},
    "Cancer_Luminal": {"Major": ["EPCAM-1"]},
}


def load_nuclei_counts(nuclei_file: str, adata_obs_names: pd.Index) -> pd.Series:
    """
    Load nuclei counts from CSV file.

    Expected format: CSV with spot names as index and integer counts.

    Args:
        nuclei_file: Path to nuclei counts CSV
        adata_obs_names: Spot names from AnnData to align

    Returns:
        pd.Series with spot names as index and nuclei counts as values
    """
    df = pd.read_csv(nuclei_file, index_col=0)
    if "nuclei_count" in df.columns:
        counts = df["nuclei_count"]
    elif df.shape[1] == 1:
        counts = df.iloc[:, 0]
    else:
        raise ValueError(f"Cannot determine nuclei count column in {nuclei_file}")

    # Align to adata spots
    common_spots = adata_obs_names.intersection(counts.index)
    if len(common_spots) < len(adata_obs_names):
        logger.warning(
            f"Only {len(common_spots)}/{len(adata_obs_names)} spots have nuclei counts"
        )
    return counts.loc[common_spots].astype(int)


def run_discrete_pipeline(
    sample_name: str,
    output_dir: str,
    nuclei_file: str = None,
    use_obs_nuclei: bool = False,
    cell_profile_dict: dict = None,
    min_counts: int = 100,
    max_em_iterations: int = 20,
    run_gex: bool = True,
) -> dict:
    """
    Run discrete cell assignment pipeline.

    Args:
        sample_name: Sample identifier (e.g., "HCC22-088-P1-S1")
        output_dir: Output directory for results
        nuclei_file: Path to nuclei counts CSV (optional if use_obs_nuclei)
        use_obs_nuclei: If True, use 'nuclei_count' or 'n_cells' from adata.obs
        cell_profile_dict: Cell type marker profile dictionary
        min_counts: Minimum counts filter for preprocessing
        max_em_iterations: Maximum EM iterations
        run_gex: Whether to run GEX deconvolution after cell assignment

    Returns:
        Dictionary with results and metrics
    """
    os.makedirs(output_dir, exist_ok=True)

    # Default to unified profile
    if cell_profile_dict is None:
        cell_profile_dict = UNIFIED_PROFILE

    # Load data
    sample_path = DATA_ROOT / sample_name / "outs"
    if not sample_path.exists():
        raise FileNotFoundError(f"Sample not found: {sample_path}")

    logger.info(f"Loading sample: {sample_name}")
    adata = sq.read.visium(
        sample_path,
        counts_file="filtered_feature_bc_matrix.h5",
        load_images=True,
        gex_only=False,  # Load both GEX and antibody
    )

    # Initialize model
    model = CitegeistModel(
        sample_name=sample_name,
        output_folder=output_dir,
        adata=adata,
    )

    # Split into GEX and antibody AnnData objects
    model.split_adata()

    # Preprocess
    model.filter_gex(min_counts=min_counts)
    model.copy_gex_to_protein_adata()
    model.preprocess_gex()

    # IMPORTANT: Use discrete preprocessing (preserves cellularity signal)
    model.preprocess_antibody_discrete()

    # Load cell profiles
    model.load_cell_profile_dict(cell_profile_dict)

    # Get nuclei counts
    if nuclei_file:
        nuclei_counts = load_nuclei_counts(
            nuclei_file, model.antibody_capture_adata.obs_names
        )
        # Add to adata.obs for the model to use
        model.antibody_capture_adata.obs["nuclei_count"] = nuclei_counts
    elif use_obs_nuclei:
        # Check for existing nuclei count column
        obs_cols = model.antibody_capture_adata.obs.columns
        if "nuclei_count" in obs_cols:
            logger.info("Using 'nuclei_count' from adata.obs")
        elif "n_cells" in obs_cols:
            logger.info("Using 'n_cells' from adata.obs as nuclei_count")
            model.antibody_capture_adata.obs["nuclei_count"] = (
                model.antibody_capture_adata.obs["n_cells"]
            )
        else:
            raise ValueError(
                "No nuclei count column found in adata.obs. "
                "Expected 'nuclei_count' or 'n_cells'."
            )
    else:
        raise ValueError(
            "Must provide either --nuclei-file or --use-obs-nuclei"
        )

    # Run discrete cell assignment (Phase 1 alternative)
    logger.info("Running discrete cell assignment (IQP + EM)...")
    start_time = datetime.now()

    cell_counts_df = model.run_discrete_cell_assignment(
        max_em_iterations=max_em_iterations,
    )

    phase1_time = (datetime.now() - start_time).total_seconds()
    logger.info(f"Phase 1 (discrete) completed in {phase1_time:.1f}s")

    # Log cell type distribution
    total_cells = cell_counts_df.values.sum()
    logger.info(f"Total cells assigned: {total_cells}")
    for ct in cell_counts_df.columns:
        ct_total = cell_counts_df[ct].sum()
        pct = 100 * ct_total / total_cells if total_cells > 0 else 0
        logger.info(f"  {ct}: {ct_total} cells ({pct:.1f}%)")

    results = {
        "sample_name": sample_name,
        "n_spots": len(cell_counts_df),
        "n_cell_types": len(cell_counts_df.columns),
        "total_cells": int(total_cells),
        "phase1_time_seconds": phase1_time,
        "cell_type_distribution": {
            ct: int(cell_counts_df[ct].sum()) for ct in cell_counts_df.columns
        },
    }

    # Run GEX deconvolution (Phase 2)
    if run_gex:
        logger.info("Running GEX deconvolution (Phase 2)...")
        start_time = datetime.now()

        model.run_cell_expression_pass1(
            use_discrete_mode=True,
            cell_counts=cell_counts_df,
        )

        phase2_time = (datetime.now() - start_time).total_seconds()
        logger.info(f"Phase 2 (GEX) completed in {phase2_time:.1f}s")
        results["phase2_time_seconds"] = phase2_time

    # Save results summary
    results_path = os.path.join(output_dir, f"{sample_name}_discrete_results.json")
    with open(results_path, "w") as f:
        json.dump(results, f, indent=2)
    logger.info(f"Saved results to {results_path}")

    return results


def main():
    parser = argparse.ArgumentParser(
        description="Run CITEgeist discrete cell assignment"
    )
    parser.add_argument(
        "--sample",
        type=str,
        required=True,
        help="Sample name (e.g., HCC22-088-P1-S1)",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="output/module3_discrete",
        help="Output directory",
    )
    parser.add_argument(
        "--nuclei-file",
        type=str,
        default=None,
        help="CSV file with nuclei counts per spot",
    )
    parser.add_argument(
        "--use-obs-nuclei",
        action="store_true",
        help="Use nuclei_count or n_cells from adata.obs",
    )
    parser.add_argument(
        "--min-counts",
        type=int,
        default=100,
        help="Minimum counts filter (default: 100 for biopsy, 25 for surgical)",
    )
    parser.add_argument(
        "--max-em-iter",
        type=int,
        default=20,
        help="Maximum EM iterations (default: 20)",
    )
    parser.add_argument(
        "--skip-gex",
        action="store_true",
        help="Skip GEX deconvolution (Phase 2)",
    )

    args = parser.parse_args()

    results = run_discrete_pipeline(
        sample_name=args.sample,
        output_dir=args.output_dir,
        nuclei_file=args.nuclei_file,
        use_obs_nuclei=args.use_obs_nuclei,
        min_counts=args.min_counts,
        max_em_iterations=args.max_em_iter,
        run_gex=not args.skip_gex,
    )

    logger.info("Pipeline completed successfully")
    return results


if __name__ == "__main__":
    main()
