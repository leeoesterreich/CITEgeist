#!/usr/bin/env python
"""
Module 3 Runner with Unified Profile.

Runs CITEgeist Module 3 (cell proportion + gene expression deconvolution)
using the unified profile derived from Module 1-2 discovery analysis.

Usage:
    python run_module3_unified.py --sample HCC22-088-P1-S1 --output-dir output/module3_unified
"""
import os
import sys
import argparse
import json
import logging
from datetime import datetime
from pathlib import Path

import numpy as np
import scanpy as sc
import pandas as pd
import squidpy as sq
import scipy.sparse

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

# All 14 patient samples
ALL_SAMPLES = [
    "HCC22-088-P1-S1",
    "HCC22-088-P1-S2",
    "HCC22-088-P2-S1",
    "HCC22-088-P2-S2",
    "HCC22-088-P3-S1_A",
    "HCC22-088-P3-S2",
    "HCC22-088-P4-S1",
    "HCC22-088-P4-S2",
    "HCC22-088-P4-S2_1i_rep",
    "HCC22-088-P5-S1",
    "HCC22-088-P5-S2",
    "HCC22-088-P5-S2_F_rep",
    "HCC22-088-P6-S1",
    "HCC22-088-P6-S2_D",
]

# Sample metadata (radius is auto-detected from spatial coordinates)
SAMPLE_METADATA = {
    "HCC22-088-P1-S1": {"type": "biopsy", "min_counts": 100},
    "HCC22-088-P1-S2": {"type": "surgical", "min_counts": 25},
    "HCC22-088-P2-S1": {"type": "biopsy", "min_counts": 100},
    "HCC22-088-P2-S2": {"type": "surgical", "min_counts": 25},
    "HCC22-088-P3-S1_A": {"type": "biopsy", "min_counts": 100},
    "HCC22-088-P3-S2": {"type": "surgical", "min_counts": 25},
    "HCC22-088-P4-S1": {"type": "biopsy", "min_counts": 100},
    "HCC22-088-P4-S2": {"type": "surgical", "min_counts": 25},
    "HCC22-088-P4-S2_1i_rep": {"type": "surgical", "min_counts": 25},
    "HCC22-088-P5-S1": {"type": "biopsy", "min_counts": 100},
    "HCC22-088-P5-S2": {"type": "surgical", "min_counts": 25},
    "HCC22-088-P5-S2_F_rep": {"type": "surgical", "min_counts": 25},
    "HCC22-088-P6-S1": {"type": "biopsy", "min_counts": 100},
    "HCC22-088-P6-S2_D": {"type": "surgical", "min_counts": 25},
}

# Unified profile from Module 1-2 discovery analysis (expanded version)
# Note: markers have -1 suffix as per real patient data convention
UNIFIED_PROFILE = {
    "Endothelial": {
        "Major": ["PECAM1-1"],
    },
    "Fibroblasts": {
        "Major": ["ACTA2-1"],
    },
    "B_Cells": {
        "Major": ["CD19-1"],
    },
    "Macrophages": {
        "Major": ["CD68-1", "CD163-1"],
    },
    "Monocytes": {
        "Major": ["CD14-1"],
    },
    "CD8_T_Cells": {
        "Major": ["CD3E-1", "CD8A-1"],
    },
    "CD4_T_Cells": {
        "Major": ["CD4-1"],
    },
    "Cancer_Luminal": {
        "Major": ["EPCAM-1"],
    },
    "Cancer_Basal": {
        "Major": ["KRT5-1", "SDC1-1"],
    },
    "Dendritic_Cells": {
        "Major": ["ITGAX-1", "HLA-DRA-1"],
    },
}


def run_module3(
    sample_name: str,
    output_dir: Path,
    profile: dict = None,
    skip_gex: bool = False,
) -> dict:
    """
    Run Module 3 (cell proportions + gene expression deconvolution).

    Args:
        sample_name: Sample identifier
        output_dir: Output directory
        profile: Cell profile dictionary (default: UNIFIED_PROFILE)
        skip_gex: If True, skip gene expression deconvolution (faster)

    Returns:
        Dict with result summary
    """
    if profile is None:
        profile = UNIFIED_PROFILE

    # Get sample-specific parameters (radius auto-detected from spatial coords)
    meta = SAMPLE_METADATA.get(sample_name, {"min_counts": 100})
    min_counts = meta["min_counts"]

    logger.info(f"Processing {sample_name} (radius=auto, min_counts={min_counts})")

    # Load data
    sample_path = DATA_ROOT / sample_name / "outs"
    if not sample_path.exists():
        raise FileNotFoundError(f"Sample path not found: {sample_path}")

    logger.info(f"Loading data from {sample_path}")
    adata = sq.read.visium(
        str(sample_path),
        counts_file='filtered_feature_bc_matrix.h5',
        load_images=True,
        gex_only=False,
    )

    # Filter spots with NaN spatial coordinates (common in surgical samples)
    if 'spatial' in adata.obsm:
        coords = adata.obsm['spatial']
        valid_mask = np.all(np.isfinite(coords), axis=1)
        n_invalid = (~valid_mask).sum()
        if n_invalid > 0:
            logger.warning(f"Filtering {n_invalid} spots with NaN/inf spatial coordinates")
            adata = adata[valid_mask].copy()
            logger.info(f"Remaining spots after coordinate filtering: {adata.shape[0]}")

    # Initialize model
    model = CitegeistModel(
        sample_name=sample_name,
        adata=adata,
        output_folder=str(output_dir),
    )

    # Split and preprocess
    model.split_adata()
    model.filter_gex(nonzero_percentage=0.01, mean_expression_threshold=1.1, min_counts=min_counts)
    model.copy_gex_to_protein_adata()
    model.preprocess_gex()
    model.preprocess_antibody()

    logger.info(f"After preprocessing: {model.gene_expression_adata.shape[0]} spots, "
                f"{model.gene_expression_adata.shape[1]} genes")

    # Load unified profile
    model.load_cell_profile_dict(profile)
    logger.info(f"Loaded unified profile with {len(profile)} cell types: {list(profile.keys())}")

    # Run cell proportion model (Pass 1) - radius auto-detected
    logger.info("Running cell proportion model...")
    global_props, finetuned_props = model.run_cell_proportion_model(
        validation_warn_only=True,  # Don't fail on rare cell types
    )

    # Append proportions to adata
    model.append_proportions_to_adata(key='global')
    model.append_proportions_to_adata(key='finetuned')

    result = {
        "sample_name": sample_name,
        "n_spots": model.gene_expression_adata.shape[0],
        "n_genes": model.gene_expression_adata.shape[1],
        "cell_types": list(profile.keys()),
        "n_cell_types": len(profile),
        "proportions_global": f"{sample_name}_cell_prop_global_results.csv",
        "proportions_finetuned": f"{sample_name}_cell_prop_finetuned_results.csv",
    }

    # Run gene expression deconvolution (Pass 2) unless skipped - radius auto-detected
    if not skip_gex:
        logger.info("Running gene expression deconvolution...")
        pass1_results = model.run_cell_expression_pass1(
            max_workers=None,
            checkpoint_interval=100,
            output_dir=str(output_dir / "checkpoints"),
            rerun=True,
        )

        model.append_gex_to_adata(pass_number=1)
        result["gex_pass1"] = f"{sample_name}_gene_expression_pass1.parquet"

    # Get final adata
    final_adata = model.get_adata()

    # Save adata with all results
    adata_path = output_dir / f"{sample_name}_module3_results.h5ad"
    final_adata.write_h5ad(adata_path)
    result["adata_path"] = str(adata_path)

    logger.info(f"Saved results to {output_dir}")

    return result


def main():
    parser = argparse.ArgumentParser(
        description="Run Module 3 with unified profile"
    )
    parser.add_argument(
        "--sample",
        type=str,
        required=True,
        help=f"Sample name (one of {ALL_SAMPLES})",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="output/module3_unified",
        help="Output directory for results",
    )
    parser.add_argument(
        "--skip-gex",
        action="store_true",
        help="Skip gene expression deconvolution (proportions only)",
    )
    parser.add_argument(
        "--profile-json",
        type=str,
        help="Optional path to custom profile JSON file",
    )

    args = parser.parse_args()

    # Validate sample name
    if args.sample not in ALL_SAMPLES:
        logger.error(f"Unknown sample: {args.sample}")
        logger.error(f"Valid samples: {ALL_SAMPLES}")
        sys.exit(1)

    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load custom profile if specified
    profile = UNIFIED_PROFILE
    if args.profile_json:
        logger.info(f"Loading custom profile from {args.profile_json}")
        with open(args.profile_json) as f:
            raw_profile = json.load(f)
        # Convert to CITEgeist format (add Major/Minor structure if missing)
        profile = {}
        for ct, markers in raw_profile.items():
            if isinstance(markers, list):
                profile[ct] = {"Major": markers}
            else:
                profile[ct] = markers

    # Run Module 3
    result = run_module3(
        sample_name=args.sample,
        output_dir=output_dir,
        profile=profile,
        skip_gex=args.skip_gex,
    )

    # Save metadata
    result["timestamp"] = datetime.now().isoformat()
    result["profile"] = {ct: list(v.get("Major", [])) for ct, v in profile.items()}

    meta_path = output_dir / f"{args.sample}_module3_meta.json"
    with open(meta_path, "w") as f:
        json.dump(result, f, indent=2)

    # Print summary
    print(f"\n{'='*60}")
    print(f"Module 3 Complete: {args.sample}")
    print(f"{'='*60}")
    print(f"Spots: {result['n_spots']}")
    print(f"Genes: {result['n_genes']}")
    print(f"Cell types: {result['n_cell_types']}")
    for ct in result['cell_types']:
        print(f"  - {ct}")
    print(f"\nOutput: {output_dir}")


if __name__ == "__main__":
    main()
