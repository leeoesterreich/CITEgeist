#!/usr/bin/env python
"""Phase 2: Module 3 deconvolution (cell-type proportion estimation) for patient samples.

Loads SpaceRanger data, runs CITEgeist Module 3 (Pass 1 proportions + Pass 2
gene expression deconvolution).  Two variants are supported:
  - baseline:     standard continuous proportion model
  - cellularity:  inject StarDist nuclei counts as a prior

Artifacts saved to --output-dir/<sample_name>/:
    <sample_name>_cell_prop_global_results.csv
    <sample_name>_cell_prop_finetuned_results.csv
    <sample_name>_gene_expression_pass1.parquet
    .phase2_complete

Usage:
    python run_patient_phase2.py --sample HCC22-088-P1-S1
    python run_patient_phase2.py --sample HCC22-088-P1-S1 --use-cellularity-prior
"""
import argparse
import logging
import os
import sys
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import numpy as np
import pandas as pd
import squidpy as sq
from model import CitegeistModel
from model.unified_config import CELL_PROFILES_NESTED

# Samples containing "S1" in name are core biopsies (higher quality / more counts).
# Samples with "S2" are surgical excisions (lower sequencing depth).
_BIOPSY_KEYWORDS = ("-S1",)


def _is_biopsy(sample_name):
    """Return True if the sample is a biopsy (S1 timepoint)."""
    return any(kw in sample_name for kw in _BIOPSY_KEYWORDS)


def run_phase2(sample_name, output_dir, seg_dir, data_dir, use_cellularity_prior=False):
    """Run Module 3 deconvolution for a single patient sample.

    Args:
        sample_name (str): Sample identifier, e.g. 'HCC22-088-P1-S1'.
        output_dir (str or Path): Root directory for phase2 outputs.
        seg_dir (str or Path): Root directory of phase1 outputs (for nuclei counts).
        data_dir (str or Path): Root of SpaceRanger processed_files directory.
        use_cellularity_prior (bool): If True, inject StarDist nuclei counts into
            the model as a cellularity prior.
    """
    output_dir = Path(output_dir)
    seg_dir = Path(seg_dir)
    data_dir = Path(data_dir)

    sample_out = output_dir / sample_name
    sample_out.mkdir(parents=True, exist_ok=True)

    marker_file = sample_out / ".phase2_complete"
    if marker_file.exists():
        logger.info("Phase 2 already complete for %s, skipping", sample_name)
        return

    # ----------------------------------------------------------------
    # Validate Phase 1 completion
    # ----------------------------------------------------------------
    phase1_marker = seg_dir / sample_name / ".phase1_complete"
    if not phase1_marker.exists():
        raise RuntimeError(f"Phase 1 not complete for {sample_name}. " f"Expected marker: {phase1_marker}")

    # ----------------------------------------------------------------
    # Load SpaceRanger data
    # ----------------------------------------------------------------
    sample_path = data_dir / sample_name / "outs"
    logger.info("Loading SpaceRanger data from %s", sample_path)
    adata = sq.read.visium(
        str(sample_path),
        counts_file="filtered_feature_bc_matrix.h5",
        load_images=True,
        gex_only=False,  # CRITICAL: must be False to load antibody capture data
    )

    # ----------------------------------------------------------------
    # Initialise CitegeistModel (real patient data, NOT simulation mode)
    # ----------------------------------------------------------------
    model = CitegeistModel(
        sample_name=sample_name,
        adata=adata,
        output_folder=str(sample_out),
    )
    model.split_adata()  # Splits combined adata into GEX + antibody AnnData

    # Filter NaN spatial coordinates (causes cKDTree / kneighbors_graph failures)
    for attr in ("gene_expression_adata", "antibody_capture_adata"):
        ad = getattr(model, attr, None)
        if ad is not None and "spatial" in ad.obsm:
            finite_mask = np.isfinite(ad.obsm["spatial"]).all(axis=1)
            n_nan = int((~finite_mask).sum())
            if n_nan > 0:
                logger.warning(
                    "Filtering %d spots with NaN spatial coords from %s",
                    n_nan,
                    attr,
                )
                setattr(model, attr, ad[finite_mask].copy())

    # min_counts: biopsy (S1) = 100, surgical (S2) = 25
    min_counts = 100 if _is_biopsy(sample_name) else 25
    logger.info("Preprocessing GEX with min_counts=%d", min_counts)
    model.filter_gex(min_counts=min_counts)
    model.copy_gex_to_protein_adata()  # align antibody spots to filtered GEX spots
    model.preprocess_gex()
    model.preprocess_antibody()

    # ----------------------------------------------------------------
    # Cell type profiles (patient data uses "-1" antibody suffix)
    # ----------------------------------------------------------------
    model.load_cell_profile_dict(CELL_PROFILES_NESTED)

    # ----------------------------------------------------------------
    # Optional: inject nuclei counts as cellularity prior
    # ----------------------------------------------------------------
    if use_cellularity_prior:
        nuclei_per_spot_path = seg_dir / sample_name / "segmentation" / "nuclei_per_spot.csv"
        if not nuclei_per_spot_path.exists():
            logger.warning(
                "Nuclei-per-spot file not found at %s; " "falling back to no cellularity prior",
                nuclei_per_spot_path,
            )
        else:
            nps_df = pd.read_csv(nuclei_per_spot_path)
            nuclei_series = nps_df.set_index("spot_barcode")["n_nuclei"]
            # Align to current GEX adata barcodes
            obs_barcodes = model.gene_expression_adata.obs_names
            nuclei_aligned = nuclei_series.reindex(obs_barcodes).fillna(0).astype(int)
            model.gene_expression_adata.obs["n_nuclei"] = nuclei_aligned.values
            logger.info(
                "Injected nuclei counts: total=%d across %d spots",
                int(nuclei_aligned.sum()),
                len(obs_barcodes),
            )

    # ----------------------------------------------------------------
    # Run Module 3 Pass 1 (proportions) + Pass 2 (GEX deconvolution)
    # ----------------------------------------------------------------
    logger.info("Running Module 3 proportion model for %s", sample_name)
    model.run_cell_proportion_model(
        validation_warn_only=True,
    )

    # NOTE: run_cell_expression_pass1() is intentionally omitted here.
    # Phase 2 only needs proportions for the SC assignment pipeline (Phase 4).
    # GEX deconvolution is skipped because it requires the same spot set as
    # the filtered GEX, which diverges from the proportion spot set and causes
    # a mismatch error ("Mismatch: gene_expression has N spots, cell_proportions
    # has M spots").

    marker_file.touch()
    logger.info("Phase 2 complete for %s", sample_name)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Patient pipeline Phase 2: Module 3 deconvolution")
    parser.add_argument("--sample", required=True, help="Sample name, e.g. HCC22-088-P1-S1")
    parser.add_argument(
        "--output-dir",
        default="output/patient_pipeline/phase2/baseline",
        help="Root directory for phase2 outputs",
    )
    parser.add_argument(
        "--seg-dir",
        default="output/patient_pipeline/phase1",
        help="Root directory of phase1 segmentation outputs",
    )
    parser.add_argument(
        "--data-dir",
        default="/ix1/alee/LO_LAB/General/Lab_Data/" "20250210_CITEGeistPublicData_GEO_Alex/processed_files",
        help="Root of SpaceRanger processed_files directory",
    )
    parser.add_argument(
        "--use-cellularity-prior",
        action="store_true",
        help="Inject StarDist nuclei counts as a cellularity prior",
    )
    args = parser.parse_args()
    run_phase2(
        args.sample,
        args.output_dir,
        args.seg_dir,
        args.data_dir,
        use_cellularity_prior=args.use_cellularity_prior,
    )
