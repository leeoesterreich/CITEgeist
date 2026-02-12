#!/usr/bin/env python
"""
Discrete cell assignment benchmark on scCube simulation data.

Pipeline:
1. Load synthetic Cellpose-compatible image
2. Run Cellpose segmentation with appropriate model (nuclei or cyto2)
3. Compare predicted vs ground truth nuclei counts
4. Run discrete cell assignment (IQP)
5. Evaluate against ground truth proportions
6. Run GEX deconvolution and evaluate
"""

import argparse
import json
import logging
import sys
import time
from pathlib import Path
from typing import Any, Dict, Optional

import cv2
import numpy as np
import pandas as pd
import scanpy as sc

# Add paths
REPO_ROOT = Path(__file__).parent.parent.parent.parent.parent
sys.path.insert(0, str(REPO_ROOT))

BENCHMARK_ROOT = Path(__file__).parent.parent.parent
sys.path.insert(0, str(BENCHMARK_ROOT))

from CITEgeist.model.citegeist_model import CitegeistModel
from CITEgeist.model.segmentation import (
    assign_nuclei_centroids_to_spots,
    run_cellpose_nuclei_segmentation,
)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)

# Default paths
DEFAULT_SCCUBE_DIR = Path(
    "/ix1/alee/LO_LAB/Personal/Brent_Schlegel/bts76/Projects/CITEgeist/"
    "Wu_Visium/Simulations/scCube_12k/replicates"
)
DEFAULT_IMAGE_DIR = REPO_ROOT / "Benchmarking/simulation_benchmarking/CITEgeist"
DEFAULT_H5AD_DIR = DEFAULT_SCCUBE_DIR  # h5ad_objects are under condition folders

# Simulation cell type profile (9 cell types in simulation)
SIMULATION_CELL_PROFILE_DICT = {
    "B-cells": ["B-cells_Protein_1", "B-cells_Protein_2"],
    "CAFs": ["CAFs_Protein_1", "CAFs_Protein_2"],
    "Cancer Epithelial": ["Cancer Epithelial_Protein_1", "Cancer Epithelial_Protein_2"],
    "Endothelial": ["Endothelial_Protein_1", "Endothelial_Protein_2"],
    "Myeloid": ["Myeloid_Protein_1", "Myeloid_Protein_2"],
    "Normal Epithelial": ["Normal Epithelial_Protein_1", "Normal Epithelial_Protein_2"],
    "PVL": ["PVL_Protein_1", "PVL_Protein_2"],
    "Plasmablasts": ["Plasmablasts_Protein_1", "Plasmablasts_Protein_2"],
    "T-cells": ["T-cells_Protein_1", "T-cells_Protein_2"],
}

# Model type mapping
MODE_TO_MODEL = {
    "dapi": "nuclei",
    "h_and_e": "cyto2",
}


def load_synthetic_image(
    replicate_id: int,
    condition: str,
    mode: str,
    image_dir: Path,
) -> np.ndarray:
    """Load synthetic Cellpose-compatible image."""
    image_path = (
        image_dir / condition / "images" / mode / f"Wu_rep_{replicate_id}_cellpose.png"
    )
    if not image_path.exists():
        raise FileNotFoundError(f"Image not found: {image_path}")

    logger.info("Loading image: %s", image_path)
    bgr = cv2.imread(str(image_path), cv2.IMREAD_COLOR)
    if bgr is None:
        raise ValueError(f"Failed to load image: {image_path}")

    rgb = cv2.cvtColor(bgr, cv2.COLOR_BGR2RGB)
    logger.info("Image shape: %s", rgb.shape)
    return rgb


def load_ground_truth_counts(
    replicate_id: int,
    condition: str,
    image_dir: Path,
) -> pd.Series:
    """Load ground truth nuclei counts per spot."""
    counts_path = (
        image_dir / condition / "nuclei_counts" / f"Wu_rep_{replicate_id}_nuclei_counts.csv"
    )
    if not counts_path.exists():
        raise FileNotFoundError(f"Nuclei counts not found: {counts_path}")

    counts = pd.read_csv(counts_path, index_col=0).squeeze()
    logger.info("Loaded ground truth counts: %d spots, total %d cells", len(counts), counts.sum())
    return counts


def load_ground_truth_proportions(
    replicate_id: int,
    condition: str,
    sccube_dir: Path,
) -> pd.DataFrame:
    """Load ground truth cell type proportions per spot."""
    prop_path = sccube_dir / condition / "ST_sim" / f"Wu_ST_{replicate_id}_prop.csv"
    if not prop_path.exists():
        raise FileNotFoundError(f"Proportions not found: {prop_path}")

    df = pd.read_csv(prop_path, index_col=0)
    # Drop coordinate columns if present
    df = df.drop(columns=["spot_x", "spot_y"], errors="ignore")
    logger.info("Loaded ground truth proportions: %d spots, %d cell types", len(df), len(df.columns))
    return df


def run_cellpose_segmentation(
    image: np.ndarray,
    mode: str,
    use_gpu: bool = False,
    diameter: Optional[float] = None,
) -> tuple:
    """Run Cellpose segmentation with mode-appropriate model."""
    model_type = MODE_TO_MODEL[mode]
    logger.info("Running Cellpose with model_type=%s", model_type)

    start = time.time()
    masks, centroids = run_cellpose_nuclei_segmentation(
        image_rgb_uint8=image,
        use_gpu=use_gpu,
        diameter=diameter,
        model_type=model_type,
    )
    elapsed = time.time() - start

    n_detected = int(masks.max()) if masks.size > 0 else 0
    logger.info("Detected %d nuclei in %.1fs", n_detected, elapsed)

    return masks, centroids, elapsed