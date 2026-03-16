#!/usr/bin/env python
"""Step 1: Re-run Module 3 with 9-type profile for unified pipeline.

Usage:
    python run_unified_step1_module3.py --sample HCC22-088-P1-S1
    python run_unified_step1_module3.py --sample HCC22-088-P1-S1 --modality he
"""
import argparse
import logging
import os
import sys
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import scanpy as sc
import squidpy as sq
from model import CitegeistModel
from model.unified_config import CELL_PROFILES_NESTED, DATA_DIR, OUTPUT_BASE


def run_step1(sample_name, modality="he", xenium_gex=None, xenium_protein=None):
    output_dir = OUTPUT_BASE / sample_name / "module3"
    output_dir.mkdir(parents=True, exist_ok=True)

    marker_file = OUTPUT_BASE / sample_name / ".step1_complete"
    if marker_file.exists():
        logger.info(f"Step 1 already complete for {sample_name}, skipping")
        return

    if modality == "he":
        sample_path = DATA_DIR / sample_name / "outs"
        logger.info(f"Loading patient data from {sample_path}")
        adata = sq.read.visium(
            str(sample_path),
            counts_file="filtered_feature_bc_matrix.h5",
            load_images=True,
            gex_only=False,
        )
        model = CitegeistModel(
            sample_name=sample_name, adata=adata, output_folder=str(output_dir),
        )
        model.split_adata()
    elif modality == "dapi":
        logger.info(f"Loading Xenium data: {xenium_gex}, {xenium_protein}")
        adata_gex = sc.read_h5ad(xenium_gex)
        adata_cite = sc.read_h5ad(xenium_protein)
        model = CitegeistModel(
            sample_name=sample_name, output_folder=str(output_dir),
            simulation=True, gene_expression_adata=adata_gex,
            antibody_capture_adata=adata_cite,
        )
    else:
        raise ValueError(f"Unknown modality: {modality}")

    # Filter spots with NaN spatial coordinates (causes cKDTree and kneighbors_graph to fail)
    import numpy as np
    for adata_attr in ("gene_expression_adata", "antibody_capture_adata"):
        ad = getattr(model, adata_attr, None)
        if ad is not None and "spatial" in ad.obsm:
            finite_mask = np.isfinite(ad.obsm["spatial"]).all(axis=1)
            n_nan = (~finite_mask).sum()
            if n_nan > 0:
                logger.warning(f"Filtering {n_nan} spots with NaN spatial coords from {adata_attr}")
                setattr(model, adata_attr, ad[finite_mask].copy())

    model.preprocess_gex()
    model.preprocess_antibody()
    model.load_cell_profile_dict(CELL_PROFILES_NESTED)

    if modality == "he":
        nuclei_counts = model.compute_spot_nuclei_counts_cellpose(
            resolution_mode="hires", use_gpu=False, save_masks=True,
        )
        logger.info(f"Cellpose: {nuclei_counts.sum():.0f} total nuclei across {len(nuclei_counts)} spots")

    model.run_cell_proportion_model(validation_warn_only=True)
    model.run_cell_expression_pass1()

    # Save results h5ad for downstream steps (e.g., PC-MIL antibody signal)
    try:
        result_adata = model.get_results_adata()
        h5ad_path = output_dir / f"{sample_name}_module3_results.h5ad"
        result_adata.write_h5ad(str(h5ad_path))
        logger.info(f"Saved results h5ad to {h5ad_path}")
    except Exception as e:
        logger.warning(f"Could not save h5ad: {e}")

    marker_file.touch()
    logger.info(f"Step 1 complete for {sample_name}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Unified pipeline Step 1: Module 3")
    parser.add_argument("--sample", required=True)
    parser.add_argument("--modality", default="he", choices=["he", "dapi"])
    parser.add_argument("--xenium-gex", help="Xenium GEX h5ad path")
    parser.add_argument("--xenium-protein", help="Xenium protein h5ad path")
    args = parser.parse_args()
    run_step1(args.sample, args.modality, args.xenium_gex, args.xenium_protein)
