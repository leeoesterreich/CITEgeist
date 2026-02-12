# Cellpose Simulation Images Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Generate Cellpose-compatible synthetic morphology images from scCube simulation data and integrate with discrete cell assignment benchmarking.

**Architecture:** New image generation script creates DAPI and H&E style images from cell coordinates. New benchmark script adapts existing Xenium discrete pipeline for simulation data. Segmentation module extended to support configurable Cellpose model type.

**Tech Stack:** PIL/Pillow, NumPy, pandas, Cellpose, existing CITEgeist model

---

## Task 1: Create Directory Structure

**Files:**
- Create: `Benchmarking/simulation_benchmarking/CITEgeist/src/` (directory)
- Create: `Benchmarking/simulation_benchmarking/CITEgeist/slurm/` (directory)

**Step 1: Create src and slurm directories**

```bash
mkdir -p Benchmarking/simulation_benchmarking/CITEgeist/src
mkdir -p Benchmarking/simulation_benchmarking/CITEgeist/slurm
```

**Step 2: Verify directories exist**

```bash
ls -la Benchmarking/simulation_benchmarking/CITEgeist/
```

Expected: `high_seg/`, `mixed/`, `src/`, `slurm/`

**Step 3: Commit**

```bash
git add Benchmarking/simulation_benchmarking/CITEgeist/src/.gitkeep
git add Benchmarking/simulation_benchmarking/CITEgeist/slurm/.gitkeep
git commit -m "chore: add src and slurm directories for simulation discrete benchmark"
```

---

## Task 2: Implement Gaussian Nucleus Rendering Functions

**Files:**
- Create: `Benchmarking/simulation_benchmarking/CITEgeist/src/generate_cellpose_images.py`

**Step 1: Create script with imports and constants**

```python
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
```

**Step 2: Run syntax check**

```bash
python -m py_compile Benchmarking/simulation_benchmarking/CITEgeist/src/generate_cellpose_images.py
```

Expected: No output (success)

**Step 3: Commit**

```bash
git add Benchmarking/simulation_benchmarking/CITEgeist/src/generate_cellpose_images.py
git commit -m "feat(benchmark): add cellpose image generation script skeleton"
```

---

## Task 3: Implement DAPI Image Generation

**Files:**
- Modify: `Benchmarking/simulation_benchmarking/CITEgeist/src/generate_cellpose_images.py`

**Step 1: Add DAPI image generation function**

Append to the file:

```python
def generate_dapi_image(
    cell_coords: np.ndarray,
    image_size: Tuple[int, int],
    sigma: float = 3.5,
    padding: int = 100,
) -> np.ndarray:
    """
    Generate DAPI-style grayscale image with Gaussian nuclei.

    Args:
        cell_coords: Nx2 array of (x, y) coordinates in original units
        image_size: (width, height) of output image
        sigma: Gaussian sigma for nuclei
        padding: Padding in pixels around coordinate extent

    Returns:
        RGB uint8 numpy array (grayscale stored as RGB)
    """
    width, height = image_size

    # Compute coordinate scaling
    x_min, x_max = cell_coords[:, 0].min(), cell_coords[:, 0].max()
    y_min, y_max = cell_coords[:, 1].min(), cell_coords[:, 1].max()

    scale_x = (width - 2 * padding) / (x_max - x_min + 1e-6)
    scale_y = (height - 2 * padding) / (y_max - y_min + 1e-6)

    # Create float accumulator for Gaussian blending
    img_float = np.zeros((height, width), dtype=np.float32)

    # Create Gaussian kernel once
    kernel = create_gaussian_kernel(sigma)
    k_size = kernel.shape[0]
    k_half = k_size // 2

    # Draw each nucleus
    for x_orig, y_orig in cell_coords:
        # Transform to image coordinates
        x_px = int(padding + (x_orig - x_min) * scale_x)
        y_px = int(padding + (y_orig - y_min) * scale_y)

        # Compute kernel bounds with clipping
        x_start = max(0, x_px - k_half)
        x_end = min(width, x_px + k_half + 1)
        y_start = max(0, y_px - k_half)
        y_end = min(height, y_px + k_half + 1)

        # Corresponding kernel region
        kx_start = x_start - (x_px - k_half)
        kx_end = k_size - ((x_px + k_half + 1) - x_end)
        ky_start = y_start - (y_px - k_half)
        ky_end = k_size - ((y_px + k_half + 1) - y_end)

        # Add Gaussian to image (max blending for overlapping nuclei)
        if x_end > x_start and y_end > y_start:
            img_float[y_start:y_end, x_start:x_end] = np.maximum(
                img_float[y_start:y_end, x_start:x_end],
                kernel[ky_start:ky_end, kx_start:kx_end] * DAPI_NUCLEUS_INTENSITY,
            )

    # Convert to uint8 RGB (grayscale replicated to 3 channels)
    img_gray = np.clip(img_float, 0, 255).astype(np.uint8)
    img_rgb = np.stack([img_gray, img_gray, img_gray], axis=-1)

    return img_rgb
```

**Step 2: Run syntax check**

```bash
python -m py_compile Benchmarking/simulation_benchmarking/CITEgeist/src/generate_cellpose_images.py
```

**Step 3: Commit**

```bash
git add Benchmarking/simulation_benchmarking/CITEgeist/src/generate_cellpose_images.py
git commit -m "feat(benchmark): add DAPI-style Gaussian nucleus rendering"
```

---

## Task 4: Implement H&E Image Generation

**Files:**
- Modify: `Benchmarking/simulation_benchmarking/CITEgeist/src/generate_cellpose_images.py`

**Step 1: Add H&E image generation function**

Append to the file:

```python
def generate_he_image(
    cell_coords: np.ndarray,
    image_size: Tuple[int, int],
    nucleus_sigma: float = 3.5,
    cytoplasm_radius: int = 10,
    padding: int = 100,
) -> np.ndarray:
    """
    Generate H&E-style image with pink cytoplasm and purple nuclei.

    Args:
        cell_coords: Nx2 array of (x, y) coordinates in original units
        image_size: (width, height) of output image
        nucleus_sigma: Gaussian sigma for nuclei
        cytoplasm_radius: Radius of cytoplasm circle in pixels
        padding: Padding in pixels around coordinate extent

    Returns:
        RGB uint8 numpy array
    """
    width, height = image_size

    # Compute coordinate scaling
    x_min, x_max = cell_coords[:, 0].min(), cell_coords[:, 0].max()
    y_min, y_max = cell_coords[:, 1].min(), cell_coords[:, 1].max()

    scale_x = (width - 2 * padding) / (x_max - x_min + 1e-6)
    scale_y = (height - 2 * padding) / (y_max - y_min + 1e-6)

    # Start with white background
    img = np.full((height, width, 3), HE_BACKGROUND, dtype=np.uint8)

    # Transform all coordinates to pixel space
    coords_px = []
    for x_orig, y_orig in cell_coords:
        x_px = int(padding + (x_orig - x_min) * scale_x)
        y_px = int(padding + (y_orig - y_min) * scale_y)
        coords_px.append((x_px, y_px))

    # First pass: draw cytoplasm circles (pink)
    for x_px, y_px in coords_px:
        y_indices, x_indices = np.ogrid[
            max(0, y_px - cytoplasm_radius) : min(height, y_px + cytoplasm_radius + 1),
            max(0, x_px - cytoplasm_radius) : min(width, x_px + cytoplasm_radius + 1),
        ]
        # Create mask for circular cytoplasm
        local_y = y_indices - y_px
        local_x = x_indices - x_px
        mask = (local_x ** 2 + local_y ** 2) <= cytoplasm_radius ** 2

        # Apply cytoplasm color where mask is True
        y_start = max(0, y_px - cytoplasm_radius)
        y_end = min(height, y_px + cytoplasm_radius + 1)
        x_start = max(0, x_px - cytoplasm_radius)
        x_end = min(width, x_px + cytoplasm_radius + 1)

        for c in range(3):
            region = img[y_start:y_end, x_start:x_end, c]
            region[mask] = HE_CYTOPLASM_COLOR[c]

    # Second pass: draw nuclei (purple Gaussian)
    kernel = create_gaussian_kernel(nucleus_sigma)
    k_size = kernel.shape[0]
    k_half = k_size // 2

    for x_px, y_px in coords_px:
        x_start = max(0, x_px - k_half)
        x_end = min(width, x_px + k_half + 1)
        y_start = max(0, y_px - k_half)
        y_end = min(height, y_px + k_half + 1)

        kx_start = x_start - (x_px - k_half)
        kx_end = k_size - ((x_px + k_half + 1) - x_end)
        ky_start = y_start - (y_px - k_half)
        ky_end = k_size - ((y_px + k_half + 1) - y_end)

        if x_end > x_start and y_end > y_start:
            k_region = kernel[ky_start:ky_end, kx_start:kx_end]
            for c in range(3):
                bg = img[y_start:y_end, x_start:x_end, c].astype(np.float32)
                nucleus = np.full_like(bg, HE_NUCLEUS_COLOR[c])
                # Alpha blend: result = nucleus * alpha + background * (1 - alpha)
                blended = nucleus * k_region + bg * (1 - k_region)
                img[y_start:y_end, x_start:x_end, c] = np.clip(blended, 0, 255).astype(
                    np.uint8
                )

    return img
```

**Step 2: Run syntax check**

```bash
python -m py_compile Benchmarking/simulation_benchmarking/CITEgeist/src/generate_cellpose_images.py
```

**Step 3: Commit**

```bash
git add Benchmarking/simulation_benchmarking/CITEgeist/src/generate_cellpose_images.py
git commit -m "feat(benchmark): add H&E-style image generation with cytoplasm and nuclei"
```

---

## Task 5: Implement Data Loading and Nuclei Count Functions

**Files:**
- Modify: `Benchmarking/simulation_benchmarking/CITEgeist/src/generate_cellpose_images.py`

**Step 1: Add data loading and nuclei count functions**

Append to the file:

```python
def load_cell_data(replicate_id: int, condition: str, input_dir: Path) -> pd.DataFrame:
    """
    Load cell-to-spot index data from scCube simulation.

    Args:
        replicate_id: Replicate number (0-4)
        condition: 'high_seg' or 'mixed'
        input_dir: Base path to scCube replicates directory

    Returns:
        DataFrame with Cell, Cell_type, point_x, point_y, spot columns
    """
    index_file = input_dir / condition / "ST_sim" / f"Wu_ST_{replicate_id}_index.csv"

    if not index_file.exists():
        raise FileNotFoundError(f"Index file not found: {index_file}")

    df = pd.read_csv(index_file, index_col=0)
    logger.info("Loaded %d cells from %s", len(df), index_file)

    # Log cell type distribution
    for ct, count in df["Cell_type"].value_counts().items():
        logger.info("  %s: %d", ct, count)

    return df


def compute_nuclei_counts(cell_df: pd.DataFrame) -> pd.Series:
    """
    Compute ground truth nuclei counts per spot.

    Args:
        cell_df: DataFrame with 'spot' column

    Returns:
        Series with spot names as index, cell counts as values
    """
    counts = cell_df.groupby("spot").size()
    counts.name = "nuclei_count"
    logger.info(
        "Nuclei counts: min=%d, max=%d, mean=%.1f, total=%d",
        counts.min(),
        counts.max(),
        counts.mean(),
        counts.sum(),
    )
    return counts
```

**Step 2: Run syntax check**

```bash
python -m py_compile Benchmarking/simulation_benchmarking/CITEgeist/src/generate_cellpose_images.py
```

**Step 3: Commit**

```bash
git add Benchmarking/simulation_benchmarking/CITEgeist/src/generate_cellpose_images.py
git commit -m "feat(benchmark): add data loading and nuclei count computation"
```

---

## Task 6: Implement Main Function with CLI

**Files:**
- Modify: `Benchmarking/simulation_benchmarking/CITEgeist/src/generate_cellpose_images.py`

**Step 1: Add main function and argument parser**

Append to the file:

```python
def generate_images_for_replicate(
    replicate_id: int,
    condition: str,
    input_dir: Path,
    output_dir: Path,
    modes: list,
    image_size: int,
    nucleus_sigma: float,
) -> None:
    """Generate images and nuclei counts for a single replicate."""
    logger.info("=== Replicate %d, Condition: %s ===", replicate_id, condition)

    # Load cell data
    cell_df = load_cell_data(replicate_id, condition, input_dir)
    cell_coords = cell_df[["point_x", "point_y"]].values

    # Compute and save nuclei counts
    nuclei_counts = compute_nuclei_counts(cell_df)
    counts_dir = output_dir / condition / "nuclei_counts"
    counts_dir.mkdir(parents=True, exist_ok=True)
    counts_path = counts_dir / f"Wu_rep_{replicate_id}_nuclei_counts.csv"
    nuclei_counts.to_csv(counts_path, header=True)
    logger.info("Saved nuclei counts: %s", counts_path)

    # Generate images for each mode
    for mode in modes:
        logger.info("Generating %s image...", mode)

        if mode == "dapi":
            img_array = generate_dapi_image(
                cell_coords,
                image_size=(image_size, image_size),
                sigma=nucleus_sigma,
            )
        elif mode == "h_and_e":
            img_array = generate_he_image(
                cell_coords,
                image_size=(image_size, image_size),
                nucleus_sigma=nucleus_sigma,
            )
        else:
            raise ValueError(f"Unknown mode: {mode}")

        # Save image
        img_dir = output_dir / condition / "images" / mode
        img_dir.mkdir(parents=True, exist_ok=True)
        img_path = img_dir / f"Wu_rep_{replicate_id}_cellpose.png"

        Image.fromarray(img_array).save(img_path)
        logger.info("Saved: %s (%dx%d)", img_path, img_array.shape[1], img_array.shape[0])


def main():
    parser = argparse.ArgumentParser(
        description="Generate Cellpose-compatible images from scCube simulation data"
    )
    parser.add_argument(
        "--replicate-id",
        type=int,
        required=True,
        help="Replicate ID (0-4)",
    )
    parser.add_argument(
        "--condition",
        choices=["high_seg", "mixed"],
        required=True,
        help="Simulation condition",
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=DEFAULT_INPUT_DIR,
        help="Path to scCube replicates directory",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Output directory for images and nuclei counts",
    )
    parser.add_argument(
        "--mode",
        choices=["dapi", "h_and_e", "both"],
        default="both",
        help="Image mode (default: both)",
    )
    parser.add_argument(
        "--image-size",
        type=int,
        default=8000,
        help="Image size in pixels (default: 8000)",
    )
    parser.add_argument(
        "--nucleus-sigma",
        type=float,
        default=3.5,
        help="Gaussian sigma for nuclei (default: 3.5)",
    )
    args = parser.parse_args()

    modes = ["dapi", "h_and_e"] if args.mode == "both" else [args.mode]

    generate_images_for_replicate(
        replicate_id=args.replicate_id,
        condition=args.condition,
        input_dir=args.input_dir,
        output_dir=args.output_dir,
        modes=modes,
        image_size=args.image_size,
        nucleus_sigma=args.nucleus_sigma,
    )

    logger.info("Done!")


if __name__ == "__main__":
    main()
```

**Step 2: Run syntax check**

```bash
python -m py_compile Benchmarking/simulation_benchmarking/CITEgeist/src/generate_cellpose_images.py
```

**Step 3: Commit**

```bash
git add Benchmarking/simulation_benchmarking/CITEgeist/src/generate_cellpose_images.py
git commit -m "feat(benchmark): complete image generation script with CLI"
```

---

## Task 7: Create SLURM Script for Image Generation

**Files:**
- Create: `Benchmarking/simulation_benchmarking/CITEgeist/slurm/sbatch_generate_images.sh`

**Step 1: Create SLURM script**

```bash
#!/bin/bash
#SBATCH --job-name=gen_cellpose_img
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --array=0-4
#SBATCH --output=logs/gen_cellpose_img_%A_%a.out
#SBATCH --error=logs/gen_cellpose_img_%A_%a.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

# Generate Cellpose-compatible images for simulation benchmarking
# Usage: sbatch sbatch_generate_images.sh high_seg
#        sbatch sbatch_generate_images.sh mixed

CONDITION=${1:-high_seg}
REPLICATE_ID=$SLURM_ARRAY_TASK_ID

eval "$(conda shell.bash hook)"
conda activate ~/alc376_bgfs/envs/CITEgeist_env

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist

OUTPUT_DIR=Benchmarking/simulation_benchmarking/CITEgeist

python Benchmarking/simulation_benchmarking/CITEgeist/src/generate_cellpose_images.py \
    --replicate-id $REPLICATE_ID \
    --condition $CONDITION \
    --output-dir $OUTPUT_DIR \
    --mode both \
    --image-size 8000

echo "Completed replicate $REPLICATE_ID for condition $CONDITION"
```

**Step 2: Make executable and verify**

```bash
chmod +x Benchmarking/simulation_benchmarking/CITEgeist/slurm/sbatch_generate_images.sh
head -20 Benchmarking/simulation_benchmarking/CITEgeist/slurm/sbatch_generate_images.sh
```

**Step 3: Create logs directory**

```bash
mkdir -p Benchmarking/simulation_benchmarking/CITEgeist/slurm/logs
```

**Step 4: Commit**

```bash
git add Benchmarking/simulation_benchmarking/CITEgeist/slurm/sbatch_generate_images.sh
git commit -m "feat(slurm): add SLURM script for cellpose image generation"
```

---

## Task 8: Add model_type Parameter to Segmentation Module

**Files:**
- Modify: `CITEgeist/model/segmentation.py`

**Step 1: Update _build_cellpose_model to accept model_type**

Find and replace in `segmentation.py`:

```python
# OLD (around line 309-318):
def _build_cellpose_model(use_gpu: bool = False):
    ...
    if hasattr(models, "Cellpose"):
        return models.Cellpose(gpu=use_gpu, model_type="nuclei")
    if hasattr(models, "CellposeModel"):
        try:
            return models.CellposeModel(gpu=use_gpu, model_type="nuclei")
        except TypeError:
            return models.CellposeModel(gpu=use_gpu, pretrained_model="nuclei")

# NEW:
def _build_cellpose_model(use_gpu: bool = False, model_type: str = "nuclei"):
    """Build Cellpose model with specified type ('nuclei' or 'cyto2')."""
    from cellpose import models

    if hasattr(models, "Cellpose"):
        return models.Cellpose(gpu=use_gpu, model_type=model_type)
    if hasattr(models, "CellposeModel"):
        try:
            return models.CellposeModel(gpu=use_gpu, model_type=model_type)
        except TypeError:
            return models.CellposeModel(gpu=use_gpu, pretrained_model=model_type)
    raise AttributeError("cellpose.models has neither `Cellpose` nor `CellposeModel`")
```

**Step 2: Update run_cellpose_nuclei_segmentation signature**

Find and update the function signature (around line 261):

```python
# OLD:
def run_cellpose_nuclei_segmentation(
    image_rgb_uint8: np.ndarray,
    use_gpu: bool = False,
    diameter: Optional[float] = None,
    flow_threshold: float = 0.4,
    cellprob_threshold: float = 0.0,
    model=None,
) -> Tuple[np.ndarray, np.ndarray]:

# NEW:
def run_cellpose_nuclei_segmentation(
    image_rgb_uint8: np.ndarray,
    use_gpu: bool = False,
    diameter: Optional[float] = None,
    flow_threshold: float = 0.4,
    cellprob_threshold: float = 0.0,
    model=None,
    model_type: str = "nuclei",
) -> Tuple[np.ndarray, np.ndarray]:
```

And update the model creation line inside the function:

```python
# OLD:
    if model is None:
        model = _build_cellpose_model(use_gpu=use_gpu)

# NEW:
    if model is None:
        model = _build_cellpose_model(use_gpu=use_gpu, model_type=model_type)
```

**Step 3: Run syntax check**

```bash
python -m py_compile CITEgeist/model/segmentation.py
```

**Step 4: Commit**

```bash
git add CITEgeist/model/segmentation.py
git commit -m "feat(segmentation): add model_type parameter for cyto2 support"
```

---

## Task 9: Create Benchmark Script for Simulation Discrete Assignment

**Files:**
- Create: `Benchmarking/simulation_benchmarking/CITEgeist/src/benchmark_discrete_simulation.py`

**Step 1: Create script with imports and constants**

```python
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
```

**Step 2: Run syntax check**

```bash
python -m py_compile Benchmarking/simulation_benchmarking/CITEgeist/src/benchmark_discrete_simulation.py
```

**Step 3: Commit**

```bash
git add Benchmarking/simulation_benchmarking/CITEgeist/src/benchmark_discrete_simulation.py
git commit -m "feat(benchmark): add simulation discrete benchmark script skeleton"
```

---

## Task 10: Add Image Loading and Cellpose Functions to Benchmark Script

**Files:**
- Modify: `Benchmarking/simulation_benchmarking/CITEgeist/src/benchmark_discrete_simulation.py`

**Step 1: Add image loading and segmentation functions**

Append to the file:

```python
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
```

**Step 2: Run syntax check**

```bash
python -m py_compile Benchmarking/simulation_benchmarking/CITEgeist/src/benchmark_discrete_simulation.py
```

**Step 3: Commit**

```bash
git add Benchmarking/simulation_benchmarking/CITEgeist/src/benchmark_discrete_simulation.py
git commit -m "feat(benchmark): add image loading and cellpose functions"
```

---

## Task 11: Add Main Pipeline and CLI to Benchmark Script

**Files:**
- Modify: `Benchmarking/simulation_benchmarking/CITEgeist/src/benchmark_discrete_simulation.py`

**Step 1: Add main pipeline function**

Append to the file:

```python
def compute_spot_coordinates(sccube_dir: Path, condition: str, replicate_id: int) -> pd.DataFrame:
    """Load spot coordinates from ST_sim meta file."""
    meta_path = sccube_dir / condition / "ST_sim" / f"Wu_ST_{replicate_id}_meta.csv"
    df = pd.read_csv(meta_path, index_col=0)
    return df[["spot", "spot_x", "spot_y"]].drop_duplicates().set_index("spot")


def run_benchmark(
    replicate_id: int,
    condition: str,
    mode: str,
    image_dir: Path,
    sccube_dir: Path,
    output_dir: Path,
    use_gpu: bool = False,
    cellpose_diameter: Optional[float] = None,
    max_em_iterations: int = 20,
) -> Dict[str, Any]:
    """Run full discrete benchmark pipeline."""
    logger.info("=" * 60)
    logger.info("BENCHMARK: replicate=%d, condition=%s, mode=%s", replicate_id, condition, mode)
    logger.info("=" * 60)

    results = {"replicate_id": replicate_id, "condition": condition, "mode": mode}
    timings = {}

    # Step 1: Load image
    image = load_synthetic_image(replicate_id, condition, mode, image_dir)

    # Step 2: Run Cellpose
    masks, centroids, cellpose_time = run_cellpose_segmentation(
        image, mode, use_gpu=use_gpu, diameter=cellpose_diameter
    )
    timings["cellpose_sec"] = cellpose_time

    # Step 3: Load spot coordinates and assign nuclei
    spot_coords_df = compute_spot_coordinates(sccube_dir, condition, replicate_id)
    spot_centers = spot_coords_df[["spot_x", "spot_y"]].values

    # Scale spot coordinates to image space (same scaling as image generation)
    # Image padding=100, image_size=8000, coord range 0-50
    padding = 100
    image_size = image.shape[0]
    coord_range = 50.0
    scale = (image_size - 2 * padding) / coord_range
    spot_centers_px = padding + spot_centers * scale

    # Estimate spot radius (based on typical spot spacing)
    spot_radius_px = scale * 0.8  # ~80% of 1 coordinate unit

    pred_counts = assign_nuclei_centroids_to_spots(
        centroids_xy=centroids,
        spot_centers_xy=spot_centers_px,
        spot_radius_px=spot_radius_px,
        spot_names=spot_coords_df.index.tolist(),
    )

    # Step 4: Compare to ground truth counts
    gt_counts = load_ground_truth_counts(replicate_id, condition, image_dir)

    # Align indices
    common_spots = pred_counts.index.intersection(gt_counts.index)
    pred_aligned = pred_counts.loc[common_spots]
    gt_aligned = gt_counts.loc[common_spots]

    count_corr = np.corrcoef(pred_aligned.values, gt_aligned.values)[0, 1]
    count_rmse = np.sqrt(np.mean((pred_aligned.values - gt_aligned.values) ** 2))

    results["nuclei_count_correlation"] = float(count_corr)
    results["nuclei_count_rmse"] = float(count_rmse)
    results["total_predicted_nuclei"] = int(pred_counts.sum())
    results["total_gt_nuclei"] = int(gt_counts.sum())

    logger.info("Nuclei count correlation: %.3f", count_corr)
    logger.info("Nuclei count RMSE: %.2f", count_rmse)

    # Step 5: Load h5ad files and run discrete assignment
    h5ad_dir = sccube_dir / condition / "h5ad_objects"
    cite_path = h5ad_dir / f"Wu_rep_{replicate_id}_CITE.h5ad"
    gex_path = h5ad_dir / f"Wu_rep_{replicate_id}_GEX.h5ad"

    if not cite_path.exists() or not gex_path.exists():
        logger.warning("h5ad files not found, skipping discrete assignment")
        results["discrete_assignment_skipped"] = True
        return results

    logger.info("Loading h5ad files...")
    adata_cite = sc.read_h5ad(cite_path)
    adata_gex = sc.read_h5ad(gex_path)

    # Initialize model
    model = CitegeistModel(
        sample_name=f"Wu_rep_{replicate_id}",
        output_folder=str(output_dir / condition / mode / f"Wu_rep_{replicate_id}"),
        simulation=True,
        gene_expression_adata=adata_gex,
        antibody_capture_adata=adata_cite,
    )

    model.preprocess_antibody_discrete()
    model.preprocess_gex(target_sum=10000)
    model.load_cell_profile_dict(SIMULATION_CELL_PROFILE_DICT)

    # Run discrete assignment
    start = time.time()
    cell_counts = model.run_discrete_cell_assignment(
        nuclei_counts=pred_aligned,
        max_iterations=max_em_iterations,
        convergence_threshold=1e-4,
        max_nuclei_cap=30,
    )
    timings["discrete_assignment_sec"] = time.time() - start

    # Compute predicted proportions
    pred_props = cell_counts.div(cell_counts.sum(axis=1), axis=0).fillna(0)

    # Load and compare to ground truth proportions
    gt_props = load_ground_truth_proportions(replicate_id, condition, sccube_dir)

    # Align to common spots and cell types
    common_spots = pred_props.index.intersection(gt_props.index)
    common_types = pred_props.columns.intersection(gt_props.columns)

    pred_props_aligned = pred_props.loc[common_spots, common_types]
    gt_props_aligned = gt_props.loc[common_spots, common_types]

    # Compute metrics
    prop_corr = np.corrcoef(
        pred_props_aligned.values.flatten(), gt_props_aligned.values.flatten()
    )[0, 1]
    prop_rmse = np.sqrt(
        np.mean((pred_props_aligned.values - gt_props_aligned.values) ** 2)
    )

    results["proportion_correlation"] = float(prop_corr)
    results["proportion_rmse"] = float(prop_rmse)
    results["timings"] = timings

    logger.info("Proportion correlation: %.3f", prop_corr)
    logger.info("Proportion RMSE: %.4f", prop_rmse)

    # Save results
    result_dir = output_dir / condition / mode / f"Wu_rep_{replicate_id}"
    result_dir.mkdir(parents=True, exist_ok=True)

    with open(result_dir / "benchmark_results.json", "w") as f:
        json.dump(results, f, indent=2)

    pred_props.to_csv(result_dir / "predicted_proportions.csv")
    cell_counts.to_csv(result_dir / "cell_counts.csv")

    logger.info("Results saved to %s", result_dir)
    return results


def main():
    parser = argparse.ArgumentParser(
        description="Run discrete cell assignment benchmark on simulation data"
    )
    parser.add_argument("--replicate-id", type=int, required=True, help="Replicate ID (0-4)")
    parser.add_argument("--condition", choices=["high_seg", "mixed"], required=True)
    parser.add_argument("--mode", choices=["dapi", "h_and_e"], required=True)
    parser.add_argument("--image-dir", type=Path, default=DEFAULT_IMAGE_DIR)
    parser.add_argument("--sccube-dir", type=Path, default=DEFAULT_SCCUBE_DIR)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--use-gpu", action="store_true")
    parser.add_argument("--cellpose-diameter", type=float, default=None)
    parser.add_argument("--max-em-iterations", type=int, default=20)
    args = parser.parse_args()

    results = run_benchmark(
        replicate_id=args.replicate_id,
        condition=args.condition,
        mode=args.mode,
        image_dir=args.image_dir,
        sccube_dir=args.sccube_dir,
        output_dir=args.output_dir,
        use_gpu=args.use_gpu,
        cellpose_diameter=args.cellpose_diameter,
        max_em_iterations=args.max_em_iterations,
    )

    print("\n=== RESULTS ===")
    print(f"Nuclei count correlation: {results.get('nuclei_count_correlation', 'N/A'):.3f}")
    print(f"Proportion correlation: {results.get('proportion_correlation', 'N/A'):.3f}")


if __name__ == "__main__":
    main()
```

**Step 2: Run syntax check**

```bash
python -m py_compile Benchmarking/simulation_benchmarking/CITEgeist/src/benchmark_discrete_simulation.py
```

**Step 3: Commit**

```bash
git add Benchmarking/simulation_benchmarking/CITEgeist/src/benchmark_discrete_simulation.py
git commit -m "feat(benchmark): complete simulation discrete benchmark pipeline"
```

---

## Task 12: Create SLURM Script for Benchmark Pipeline

**Files:**
- Create: `Benchmarking/simulation_benchmarking/CITEgeist/slurm/sbatch_discrete_simulation.sh`

**Step 1: Create SLURM script**

```bash
#!/bin/bash
#SBATCH --job-name=sim_discrete
#SBATCH --cluster=gpu
#SBATCH --partition=l40s
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --array=0-4
#SBATCH --output=logs/sim_discrete_%A_%a.out
#SBATCH --error=logs/sim_discrete_%A_%a.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

# Run discrete cell assignment benchmark on simulation data
# Usage: sbatch sbatch_discrete_simulation.sh high_seg dapi
#        sbatch sbatch_discrete_simulation.sh mixed h_and_e

CONDITION=${1:-high_seg}
MODE=${2:-dapi}
REPLICATE_ID=$SLURM_ARRAY_TASK_ID

module load gurobi/12.0.3
eval "$(conda shell.bash hook)"
conda activate ~/alc376_bgfs/envs/CITEgeist_env

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist

OUTPUT_DIR=Benchmarking/simulation_benchmarking/CITEgeist/output_discrete

python Benchmarking/simulation_benchmarking/CITEgeist/src/benchmark_discrete_simulation.py \
    --replicate-id $REPLICATE_ID \
    --condition $CONDITION \
    --mode $MODE \
    --output-dir $OUTPUT_DIR \
    --use-gpu \
    --max-em-iterations 20

echo "Completed replicate $REPLICATE_ID, condition $CONDITION, mode $MODE"
```

**Step 2: Make executable**

```bash
chmod +x Benchmarking/simulation_benchmarking/CITEgeist/slurm/sbatch_discrete_simulation.sh
```

**Step 3: Commit**

```bash
git add Benchmarking/simulation_benchmarking/CITEgeist/slurm/sbatch_discrete_simulation.sh
git commit -m "feat(slurm): add SLURM script for simulation discrete benchmark"
```

---

## Task 13: Test Image Generation Locally

**Step 1: Run image generation for one replicate (dry run)**

```bash
cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist

eval "$(conda shell.bash hook)"
conda activate ~/alc376_bgfs/envs/CITEgeist_env

python Benchmarking/simulation_benchmarking/CITEgeist/src/generate_cellpose_images.py \
    --replicate-id 0 \
    --condition high_seg \
    --output-dir Benchmarking/simulation_benchmarking/CITEgeist \
    --mode both \
    --image-size 8000
```

Expected output:
- Creates `Benchmarking/simulation_benchmarking/CITEgeist/high_seg/images/dapi/Wu_rep_0_cellpose.png`
- Creates `Benchmarking/simulation_benchmarking/CITEgeist/high_seg/images/h_and_e/Wu_rep_0_cellpose.png`
- Creates `Benchmarking/simulation_benchmarking/CITEgeist/high_seg/nuclei_counts/Wu_rep_0_nuclei_counts.csv`

**Step 2: Verify outputs exist**

```bash
ls -la Benchmarking/simulation_benchmarking/CITEgeist/high_seg/images/dapi/
ls -la Benchmarking/simulation_benchmarking/CITEgeist/high_seg/images/h_and_e/
ls -la Benchmarking/simulation_benchmarking/CITEgeist/high_seg/nuclei_counts/
```

**Step 3: Visually inspect images**

Open the generated PNGs to verify:
- DAPI: White/gray Gaussian dots on black background
- H&E: Pink circles with purple centers on white background

---

## Task 14: Submit SLURM Jobs for Full Image Generation

**Step 1: Submit high_seg image generation**

```bash
cd Benchmarking/simulation_benchmarking/CITEgeist/slurm
sbatch sbatch_generate_images.sh high_seg
```

**Step 2: Submit mixed image generation**

```bash
sbatch sbatch_generate_images.sh mixed
```

**Step 3: Monitor jobs**

```bash
squeue -u $USER
```

---

## Task 15: Final Commit and Summary

**Step 1: Commit any remaining changes**

```bash
git status
git add -A
git commit -m "feat(benchmark): complete cellpose simulation image generation pipeline"
```

**Step 2: Push to remote**

```bash
git push origin dev
```

---

## Summary

| Task | Description | Files |
|------|-------------|-------|
| 1 | Create directories | `src/`, `slurm/` |
| 2-6 | Image generation script | `generate_cellpose_images.py` |
| 7 | SLURM for images | `sbatch_generate_images.sh` |
| 8 | Segmentation model_type | `segmentation.py` |
| 9-11 | Benchmark script | `benchmark_discrete_simulation.py` |
| 12 | SLURM for benchmark | `sbatch_discrete_simulation.sh` |
| 13-14 | Testing and execution | Local + SLURM |
| 15 | Final commit | Git |
