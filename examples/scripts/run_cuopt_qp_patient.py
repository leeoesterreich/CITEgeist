"""Run cuOPT QP proportions for a single patient sample.

Called by sbatch_cuopt_qp_12patient.sh array job.

Usage:
  python run_cuopt_qp_patient.py --sample sample-P1-S1
"""

import argparse
import logging
import sys
from pathlib import Path

import squidpy as sq

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent.parent
sys.path.insert(0, str(REPO_ROOT))

# CITEgeist public data — download from GEO (accession in the manuscript) and set this path
PATIENT_DATA_ROOT = Path("/path/to/CITEgeist_public_data/processed_files")
OUTPUT_ROOT = REPO_ROOT / "output" / "module3_cuopt_qp"

SAMPLES = [
    "sample-P1-S1", "sample-P1-S2",
    "sample-P2-S1", "sample-P2-S2",
    "sample-P3-S1", "sample-P3-S2",
    "sample-P4-S1", "sample-P4-S2",
    "sample-P5-S1", "sample-P5-S2",
    "sample-P6-S1", "sample-P6-S2",
]

# Major markers only — minor markers (CD45RA/CD45RO/CD16) degrade QP r by 8.6%
# Format: {"CellType": {"Major": ["Marker1-1", ...], "Minor": [...]}, ...}
PROFILE_DICT = {
    "Endothelial": {"Major": ["PECAM1-1"]},
    "Fibroblasts": {"Major": ["ACTA2-1"]},
    "B_Cells": {"Major": ["CD19-1"]},
    "Macrophages": {"Major": ["CD68-1", "CD163-1"]},
    "Monocytes": {"Major": ["CD14-1"]},
    "CD8_T_Cells": {"Major": ["CD8A-1", "CD3E-1"]},
    "CD4_T_Cells": {"Major": ["CD4-1", "CD3E-1"]},
    "Cancer_Luminal": {"Major": ["EPCAM-1"]},
    "Cancer_Basal": {"Major": ["KRT5-1", "SDC1-1", "EPCAM-1"]},
    "Dendritic_Cells": {"Major": ["ITGAX-1", "HLA-DRA-1"]},
}

# min_counts: 100 for biopsy (S1), 25 for surgical (S2)
SAMPLE_METADATA = {s: {"min_counts": 100 if "S1" in s else 25} for s in SAMPLES}

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)


def run_qp_for_sample(sample_name: str):
    from CITEgeist.model import CitegeistModel

    patient_dir = PATIENT_DATA_ROOT / sample_name / "outs"
    output_dir = OUTPUT_ROOT / sample_name
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("Loading patient data: %s", patient_dir)
    adata = sq.read.visium(
        str(patient_dir),
        counts_file="filtered_feature_bc_matrix.h5",
        load_images=True,
        gex_only=False,
    )

    model = CitegeistModel(
        sample_name=sample_name,
        adata=adata,
        output_folder=str(output_dir),
    )

    logger.info("Splitting adata...")
    model.split_adata()

    # Filter spots with NaN spatial coordinates (common in surgical samples)
    import numpy as np
    for attr in ("gene_expression_adata", "antibody_capture_adata"):
        ad = getattr(model, attr)
        if ad is not None and "spatial" in ad.obsm:
            coords = ad.obsm["spatial"]
            valid = np.all(np.isfinite(coords), axis=1)
            n_invalid = (~valid).sum()
            if n_invalid > 0:
                logger.warning("Filtering %d spots with NaN spatial coords from %s", n_invalid, attr)
                setattr(model, attr, ad[valid].copy())

    meta = SAMPLE_METADATA[sample_name]
    logger.info("Filtering + preprocessing GEX (min_counts=%d)...", meta["min_counts"])
    model.filter_gex(min_counts=meta["min_counts"])
    model.preprocess_gex()

    logger.info("Preprocessing antibody...")
    model.preprocess_antibody()

    logger.info("Loading cell profile dict (major markers only)...")
    model.load_cell_profile_dict(PROFILE_DICT)

    logger.info("Running cuOPT QP proportions...")
    model.run_cell_proportion_model(
        method="qp",
        lambda_reg=1.0,
        lambda_laplacian=0.03,
        use_detection_gating=True,
        detection_gate_ub=0.01,
    )

    logger.info("QP complete. Writing .done marker.")
    (output_dir / ".done").touch()
    logger.info("Done: %s", sample_name)


def main():
    parser = argparse.ArgumentParser(description="cuOPT QP for a single patient")
    parser.add_argument("--sample", required=True, help="Sample name")
    args = parser.parse_args()

    assert args.sample in SAMPLES, f"Unknown sample: {args.sample}"
    run_qp_for_sample(args.sample)


if __name__ == "__main__":
    main()
