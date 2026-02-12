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
