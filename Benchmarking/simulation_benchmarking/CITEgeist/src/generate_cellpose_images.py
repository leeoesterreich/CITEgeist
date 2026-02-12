#!/usr/bin/env python
"""
Generate Cellpose-compatible synthetic images from scCube simulation data.

Supports two image modes:
- dapi: Grayscale nuclei on black background (for Cellpose 'nuclei' model)
- h_and_e: H&E-style with purple nuclei and pink cytoplasm (for Cellpose 'cyto2' model)
"""

import argparse
import logging
from pathlib import Path
from typing import Tuple

import numpy as np
import pandas as pd
from PIL import Image

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)

# Default paths
DEFAULT_INPUT_DIR = Path(
    "/ix1/alee/LO_LAB/Personal/Brent_Schlegel/bts76/Projects/CITEgeist/"
    "Wu_Visium/Simulations/scCube_12k/replicates"
)

# Color constants (from actual Visium H&E and Xenium DAPI analysis)
DAPI_BACKGROUND = (0, 0, 0)
DAPI_NUCLEUS_INTENSITY = 180  # Grayscale intensity

HE_BACKGROUND = (250, 250, 250)
HE_NUCLEUS_COLOR = (140, 90, 130)  # Purple/magenta (hematoxylin)
HE_CYTOPLASM_COLOR = (240, 220, 220)  # Pale pink (eosin)


def create_gaussian_kernel(sigma: float, size: int = None) -> np.ndarray:
    """
    Create a 2D Gaussian kernel for nucleus rendering.

    Args:
        sigma: Standard deviation of Gaussian
        size: Kernel size (default: 6*sigma, must be odd)

    Returns:
        2D numpy array with Gaussian values normalized to [0, 1]
    """
    if size is None:
        size = int(np.ceil(6 * sigma))
        if size % 2 == 0:
            size += 1

    center = size // 2
    y, x = np.ogrid[:size, :size]
    kernel = np.exp(-((x - center) ** 2 + (y - center) ** 2) / (2 * sigma ** 2))
    return kernel
