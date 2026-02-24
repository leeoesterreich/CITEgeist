# Cell Morphology Watershed Segmentation Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add watershed-based cell segmentation using boundary channel, extracting cell-level morphology features to improve Module 3b single-cell type assignment.

**Architecture:** DAPI nuclei (Cellpose) serve as watershed seeds; boundary channel gradient defines basin edges. Extract 12 features (nuclear + cell + ratio) instead of current 5 nuclear-only. Soft-label classifier uses expanded features for cell type prediction.

**Tech Stack:** scikit-image (watershed, regionprops), scipy (ndimage), numpy, existing Cellpose integration

---

## Task 1: Watershed Cell Segmentation Function

**Files:**
- Create: `CITEgeist/model/watershed_segmentation.py`
- Test: `CITEgeist/tests/test_watershed_segmentation.py`

### Step 1.1: Write the failing test for basic watershed

```python
# CITEgeist/tests/test_watershed_segmentation.py
"""Tests for watershed cell segmentation."""
import numpy as np
import pytest
from CITEgeist.model.watershed_segmentation import watershed_from_nuclei


def test_watershed_single_nucleus():
    """Test watershed expansion from a single nucleus seed."""
    # 50x50 image with nucleus in center
    nucleus_mask = np.zeros((50, 50), dtype=np.int32)
    y, x = np.ogrid[:50, :50]
    nucleus_mask[((x - 25)**2 + (y - 25)**2) <= 25] = 1  # r=5

    # Boundary channel with ring around nucleus (simulating membrane)
    boundary_img = np.zeros((50, 50), dtype=np.float32)
    boundary_img[((x - 25)**2 + (y - 25)**2) > 64] = 0.0  # outside r=8
    boundary_img[((x - 25)**2 + (y - 25)**2) <= 64] = 1.0  # inside r=8
    boundary_img[((x - 25)**2 + (y - 25)**2) <= 25] = 0.5  # nucleus

    cell_mask = watershed_from_nuclei(nucleus_mask, boundary_img)

    # Cell mask should have same label (1) but larger area
    assert cell_mask.max() == 1
    assert (cell_mask == 1).sum() > (nucleus_mask == 1).sum()
    # Nucleus should be contained in cell
    assert np.all(cell_mask[nucleus_mask == 1] == 1)
```

### Step 1.2: Run test to verify it fails

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest CITEgeist/tests/test_watershed_segmentation.py::test_watershed_single_nucleus -v`

Expected: FAIL with "ModuleNotFoundError: No module named 'CITEgeist.model.watershed_segmentation'"

### Step 1.3: Write minimal implementation

```python
# CITEgeist/model/watershed_segmentation.py
"""Watershed cell segmentation using nuclei as seeds."""
import numpy as np
from scipy import ndimage
from skimage.segmentation import watershed
from skimage.filters import sobel


def watershed_from_nuclei(
    nucleus_mask: np.ndarray,
    boundary_img: np.ndarray,
    use_gradient: bool = True,
) -> np.ndarray:
    """
    Expand nuclear masks to cell boundaries using watershed.

    Args:
        nucleus_mask: 2D labeled array where each nucleus has unique integer ID.
                      Background = 0.
        boundary_img: 2D float array of boundary channel (membrane stain).
                      Higher values = membrane signal.
        use_gradient: If True, use Sobel gradient of boundary_img as elevation.
                      If False, use inverted boundary_img directly.

    Returns:
        Cell mask with same labels as nucleus_mask, expanded to cell boundaries.
    """
    if nucleus_mask.shape != boundary_img.shape:
        raise ValueError(
            f"Shape mismatch: nucleus_mask {nucleus_mask.shape} vs "
            f"boundary_img {boundary_img.shape}"
        )

    # Compute elevation map for watershed
    if use_gradient:
        # Sobel gradient - edges become ridges (high values)
        elevation = sobel(boundary_img.astype(np.float64))
    else:
        # Invert so membrane (high) becomes ridge
        elevation = -boundary_img.astype(np.float64)

    # Generate markers from nucleus centroids
    # Each nucleus label becomes a seed
    markers = nucleus_mask.copy()

    # Run watershed
    cell_mask = watershed(elevation, markers=markers)

    return cell_mask.astype(np.int32)
```

### Step 1.4: Run test to verify it passes

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest CITEgeist/tests/test_watershed_segmentation.py::test_watershed_single_nucleus -v`

Expected: PASS

### Step 1.5: Add test for multiple nuclei

```python
# Add to test_watershed_segmentation.py

def test_watershed_multiple_nuclei():
    """Test watershed with multiple nuclei - each gets separate cell."""
    nucleus_mask = np.zeros((100, 100), dtype=np.int32)
    y, x = np.ogrid[:100, :100]
    # Nucleus 1 at (25, 25)
    nucleus_mask[((x - 25)**2 + (y - 25)**2) <= 25] = 1
    # Nucleus 2 at (75, 75)
    nucleus_mask[((x - 75)**2 + (y - 75)**2) <= 25] = 2

    # Boundary: high signal between the two cells (around x=50)
    boundary_img = np.zeros((100, 100), dtype=np.float32)
    boundary_img[:, 45:55] = 1.0  # Vertical membrane between cells

    cell_mask = watershed_from_nuclei(nucleus_mask, boundary_img)

    # Both labels should exist
    assert set(np.unique(cell_mask)) == {1, 2}
    # Cell 1 should be mostly left half
    assert (cell_mask[:, :50] == 1).sum() > (cell_mask[:, :50] == 2).sum()
    # Cell 2 should be mostly right half
    assert (cell_mask[:, 50:] == 2).sum() > (cell_mask[:, 50:] == 1).sum()
```

### Step 1.6: Run both tests

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest CITEgeist/tests/test_watershed_segmentation.py -v`

Expected: 2 tests PASS

### Step 1.7: Add test for shape mismatch error

```python
# Add to test_watershed_segmentation.py

def test_watershed_shape_mismatch():
    """Test that shape mismatch raises ValueError."""
    nucleus_mask = np.zeros((50, 50), dtype=np.int32)
    boundary_img = np.zeros((60, 60), dtype=np.float32)

    with pytest.raises(ValueError, match="Shape mismatch"):
        watershed_from_nuclei(nucleus_mask, boundary_img)
```

### Step 1.8: Run all watershed tests

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest CITEgeist/tests/test_watershed_segmentation.py -v`

Expected: 3 tests PASS

### Step 1.9: Commit

```bash
git add CITEgeist/model/watershed_segmentation.py CITEgeist/tests/test_watershed_segmentation.py
git commit -m "feat: add watershed cell segmentation from nuclei seeds"
```

---

## Task 2: Cell Morphology Feature Extraction

**Files:**
- Modify: `CITEgeist/model/morphology_features.py`
- Test: `CITEgeist/tests/test_morphology_features.py`

### Step 2.1: Write failing test for cell features

```python
# Add to CITEgeist/tests/test_morphology_features.py

def test_extract_cell_features_single_cell():
    """Test cell feature extraction from nucleus + cell masks."""
    from CITEgeist.model.morphology_features import extract_cell_features

    # Create nucleus (small circle) and cell (larger circle) masks
    nucleus_mask = np.zeros((50, 50), dtype=np.int32)
    cell_mask = np.zeros((50, 50), dtype=np.int32)
    y, x = np.ogrid[:50, :50]

    # Nucleus: r=5 at center
    nucleus_mask[((x - 25)**2 + (y - 25)**2) <= 25] = 1
    # Cell: r=15 at center (same label)
    cell_mask[((x - 25)**2 + (y - 25)**2) <= 225] = 1

    features_df = extract_cell_features(nucleus_mask, cell_mask)

    assert len(features_df) == 1
    assert 'cell_id' in features_df.columns
    # Nuclear features
    assert 'nucleus_area' in features_df.columns
    assert 'nucleus_circularity' in features_df.columns
    # Cell features
    assert 'cell_area' in features_df.columns
    assert 'cell_circularity' in features_df.columns
    # Ratio features
    assert 'nc_ratio' in features_df.columns
    assert 'cytoplasm_area' in features_df.columns

    # Cell area should be larger than nucleus
    row = features_df.iloc[0]
    assert row['cell_area'] > row['nucleus_area']
    # N:C ratio should be < 1
    assert 0 < row['nc_ratio'] < 1
    # Cytoplasm = cell - nucleus
    assert row['cytoplasm_area'] == row['cell_area'] - row['nucleus_area']
```

### Step 2.2: Run test to verify it fails

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest CITEgeist/tests/test_morphology_features.py::test_extract_cell_features_single_cell -v`

Expected: FAIL with "cannot import name 'extract_cell_features'"

### Step 2.3: Implement extract_cell_features

```python
# Add to CITEgeist/model/morphology_features.py after extract_nucleus_features

def extract_cell_features(
    nucleus_mask: np.ndarray,
    cell_mask: np.ndarray,
) -> pd.DataFrame:
    """
    Extract morphology features from paired nucleus and cell masks.

    Assumes nucleus_mask and cell_mask use the same integer labels,
    where each label corresponds to a nucleus-cell pair.

    Args:
        nucleus_mask: 2D labeled array of nuclei (from Cellpose)
        cell_mask: 2D labeled array of cells (from watershed)

    Returns:
        DataFrame with columns:
            - cell_id: label from masks
            - centroid_x, centroid_y: nucleus center coordinates
            Nuclear features (prefix 'nucleus_'):
            - nucleus_area, nucleus_circularity, nucleus_eccentricity,
              nucleus_solidity, nucleus_aspect_ratio
            Cell features (prefix 'cell_'):
            - cell_area, cell_circularity, cell_eccentricity,
              cell_solidity, cell_aspect_ratio
            Ratio features:
            - nc_ratio: nucleus_area / cell_area
            - cytoplasm_area: cell_area - nucleus_area
    """
    if nucleus_mask.shape != cell_mask.shape:
        raise ValueError(
            f"Shape mismatch: nucleus_mask {nucleus_mask.shape} vs "
            f"cell_mask {cell_mask.shape}"
        )

    # Get nucleus features
    nuc_df = extract_nucleus_features(nucleus_mask)
    if len(nuc_df) == 0:
        return pd.DataFrame()

    nuc_df = nuc_df.rename(columns={
        'nucleus_id': 'cell_id',
        'area': 'nucleus_area',
        'circularity': 'nucleus_circularity',
        'eccentricity': 'nucleus_eccentricity',
        'solidity': 'nucleus_solidity',
        'aspect_ratio': 'nucleus_aspect_ratio',
    })
    # Drop columns we'll rename differently
    nuc_df = nuc_df.drop(columns=[
        'perimeter', 'major_axis_length', 'minor_axis_length'
    ], errors='ignore')

    # Get cell features using regionprops
    if cell_mask.max() == 0:
        return pd.DataFrame()

    props = regionprops_table(
        cell_mask,
        properties=[
            'label',
            'area',
            'perimeter',
            'eccentricity',
            'solidity',
            'major_axis_length',
            'minor_axis_length',
        ]
    )
    cell_df = pd.DataFrame(props)
    cell_df = cell_df.rename(columns={'label': 'cell_id'})

    # Compute cell circularity
    perimeter = cell_df['perimeter'].replace(0, np.nan)
    cell_df['cell_circularity'] = (4 * np.pi * cell_df['area']) / (perimeter ** 2)
    cell_df['cell_circularity'] = cell_df['cell_circularity'].fillna(1.0).clip(0, 1)

    # Compute cell aspect ratio
    minor = cell_df['minor_axis_length'].replace(0, np.nan)
    cell_df['cell_aspect_ratio'] = cell_df['major_axis_length'] / minor
    cell_df['cell_aspect_ratio'] = cell_df['cell_aspect_ratio'].fillna(1.0)

    # Rename and select columns
    cell_df = cell_df.rename(columns={
        'area': 'cell_area',
        'eccentricity': 'cell_eccentricity',
        'solidity': 'cell_solidity',
    })
    cell_df = cell_df[['cell_id', 'cell_area', 'cell_circularity',
                       'cell_eccentricity', 'cell_solidity', 'cell_aspect_ratio']]

    # Merge nucleus and cell features
    merged = nuc_df.merge(cell_df, on='cell_id', how='inner')

    # Compute ratio features
    merged['nc_ratio'] = merged['nucleus_area'] / merged['cell_area'].replace(0, np.nan)
    merged['nc_ratio'] = merged['nc_ratio'].fillna(1.0).clip(0, 1)
    merged['cytoplasm_area'] = merged['cell_area'] - merged['nucleus_area']
    merged['cytoplasm_area'] = merged['cytoplasm_area'].clip(lower=0)

    return merged
```

### Step 2.4: Run test to verify it passes

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest CITEgeist/tests/test_morphology_features.py::test_extract_cell_features_single_cell -v`

Expected: PASS

### Step 2.5: Add test for multiple cells

```python
# Add to test_morphology_features.py

def test_extract_cell_features_multiple_cells():
    """Test cell features with multiple cells of different shapes."""
    from CITEgeist.model.morphology_features import extract_cell_features

    nucleus_mask = np.zeros((100, 100), dtype=np.int32)
    cell_mask = np.zeros((100, 100), dtype=np.int32)
    y, x = np.ogrid[:100, :100]

    # Cell 1: Small nucleus, large cell (low N:C ratio - like macrophage)
    nucleus_mask[((x - 25)**2 + (y - 25)**2) <= 16] = 1  # r=4
    cell_mask[((x - 25)**2 + (y - 25)**2) <= 225] = 1    # r=15

    # Cell 2: Large nucleus, small cell (high N:C ratio - like lymphocyte)
    nucleus_mask[((x - 75)**2 + (y - 75)**2) <= 36] = 2  # r=6
    cell_mask[((x - 75)**2 + (y - 75)**2) <= 64] = 2     # r=8

    features_df = extract_cell_features(nucleus_mask, cell_mask)

    assert len(features_df) == 2

    cell1 = features_df[features_df['cell_id'] == 1].iloc[0]
    cell2 = features_df[features_df['cell_id'] == 2].iloc[0]

    # Cell 1 should have lower N:C ratio (more cytoplasm)
    assert cell1['nc_ratio'] < cell2['nc_ratio']
    # Cell 1 should have larger cytoplasm area
    assert cell1['cytoplasm_area'] > cell2['cytoplasm_area']
```

### Step 2.6: Run all morphology tests

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest CITEgeist/tests/test_morphology_features.py -v`

Expected: All tests PASS (existing + 2 new)

### Step 2.7: Commit

```bash
git add CITEgeist/model/morphology_features.py CITEgeist/tests/test_morphology_features.py
git commit -m "feat: add cell-level morphology feature extraction"
```

---

## Task 3: Update Module 3b to Use Cell Features

**Files:**
- Modify: `CITEgeist/model/module3b_nucleus_assignment.py`
- Test: `CITEgeist/tests/test_module3b_nucleus_assignment.py`

### Step 3.1: Write failing test for cell-feature mode

```python
# Add to CITEgeist/tests/test_module3b_nucleus_assignment.py

def test_run_nucleus_assignment_with_cell_features():
    """Test assignment using cell morphology features."""
    from CITEgeist.model.module3b_nucleus_assignment import run_nucleus_assignment

    # Create synthetic nucleus and cell masks
    nucleus_mask = np.zeros((100, 100), dtype=np.int32)
    cell_mask = np.zeros((100, 100), dtype=np.int32)
    y, x = np.ogrid[:100, :100]

    # Two nuclei/cells
    nucleus_mask[((x - 30)**2 + (y - 30)**2) <= 25] = 1
    cell_mask[((x - 30)**2 + (y - 30)**2) <= 100] = 1
    nucleus_mask[((x - 70)**2 + (y - 70)**2) <= 25] = 2
    cell_mask[((x - 70)**2 + (y - 70)**2) <= 100] = 2

    # Both in same spot
    nuclei_spot_map = pd.DataFrame({
        'nucleus_id': [1, 2],
        'spot_id': ['spot_A', 'spot_A'],
    })

    proportions = pd.DataFrame({
        'spot_id': ['spot_A'],
        'TypeA': [0.5],
        'TypeB': [0.5],
    })

    nuclei_counts = pd.Series({'spot_A': 2})
    cell_types = ['TypeA', 'TypeB']

    result = run_nucleus_assignment(
        mask=nucleus_mask,
        nuclei_spot_map=nuclei_spot_map,
        proportions=proportions,
        nuclei_counts=nuclei_counts,
        cell_types=cell_types,
        cell_mask=cell_mask,  # NEW: pass cell mask for cell features
    )

    # Should have assignments for both nuclei
    assert len(result.assignments) == 2
    # Each should be assigned to one of the types
    assert set(result.assignments.values()) <= {'TypeA', 'TypeB'}
    # Morphology features should include cell features
    assert 'cell_area' in result.morphology_features.columns
    assert 'nc_ratio' in result.morphology_features.columns
```

### Step 3.2: Run test to verify it fails

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest CITEgeist/tests/test_module3b_nucleus_assignment.py::test_run_nucleus_assignment_with_cell_features -v`

Expected: FAIL (TypeError - unexpected keyword argument 'cell_mask')

### Step 3.3: Update run_nucleus_assignment to accept cell_mask

```python
# Replace run_nucleus_assignment in CITEgeist/model/module3b_nucleus_assignment.py

from .morphology_features import extract_nucleus_features, extract_cell_features, largest_remainder_discretize


def run_nucleus_assignment(
    mask: np.ndarray,
    nuclei_spot_map: pd.DataFrame,
    proportions: pd.DataFrame,
    nuclei_counts: pd.Series,
    cell_types: List[str],
    cell_mask: Optional[np.ndarray] = None,
) -> NucleusAssignmentResult:
    """
    Run full nucleus assignment pipeline.

    Args:
        mask: Cellpose label mask (H, W) with nucleus labels
        nuclei_spot_map: DataFrame with 'nucleus_id' and 'spot_id' columns
        proportions: DataFrame with 'spot_id' and one column per cell type
        nuclei_counts: Series mapping spot_id -> nuclei count
        cell_types: List of cell type names (column names in proportions)
        cell_mask: Optional cell mask from watershed (same labels as nucleus mask).
                   If provided, extracts cell-level features. If None, uses
                   nuclear features only.

    Returns:
        NucleusAssignmentResult with assignments and metadata
    """
    # Step 1: Extract morphology features
    if cell_mask is not None:
        # Use cell-level features (12 features)
        morph_df = extract_cell_features(mask, cell_mask)
        morph_df = morph_df.rename(columns={'cell_id': 'nucleus_id'})
        feature_cols = [
            'nucleus_area', 'nucleus_circularity', 'nucleus_eccentricity',
            'nucleus_solidity', 'nucleus_aspect_ratio',
            'cell_area', 'cell_circularity', 'cell_eccentricity',
            'cell_solidity', 'cell_aspect_ratio',
            'nc_ratio', 'cytoplasm_area',
        ]
    else:
        # Use nuclear features only (5 features)
        morph_df = extract_nucleus_features(mask)
        feature_cols = ['area', 'circularity', 'eccentricity', 'solidity', 'aspect_ratio']

    # Merge with spot assignments
    morph_df = morph_df.merge(nuclei_spot_map, on='nucleus_id', how='inner')

    # Step 2: Build training data (soft labels from spot proportions)
    prop_cols = cell_types
    spot_props = proportions.set_index('spot_id')[prop_cols]

    # Get soft labels for each nucleus
    y_soft = morph_df['spot_id'].map(lambda s: spot_props.loc[s].values).values
    y_soft = np.vstack(y_soft)  # (n_nuclei, n_types)

    # Feature matrix
    feature_cols = [c for c in feature_cols if c in morph_df.columns]
    X = morph_df[feature_cols].values

    # Step 3: Train classifier
    clf = SoftLabelClassifier(n_cell_types=len(cell_types))
    clf.fit(X, y_soft)

    # Step 4: Predict probabilities
    probs = clf.predict_proba(X)
    probs_df = pd.DataFrame(probs, columns=cell_types)
    probs_df['nucleus_id'] = morph_df['nucleus_id'].values
    probs_df['spot_id'] = morph_df['spot_id'].values

    # Step 5: Per-spot Hungarian assignment
    all_assignments = {}

    for spot_id in morph_df['spot_id'].unique():
        spot_mask = morph_df['spot_id'] == spot_id
        spot_nuclei = morph_df[spot_mask]
        spot_probs = probs[spot_mask]

        # Get proportions and nuclei count for this spot
        spot_proportions = spot_props.loc[spot_id].values
        n_nuclei = int(nuclei_counts.get(spot_id, len(spot_nuclei)))

        # Convert proportions to counts
        counts = largest_remainder_discretize(spot_proportions, n_nuclei)

        # Run Hungarian assignment
        nucleus_ids = spot_nuclei['nucleus_id'].values
        type_assignments = assign_nuclei_to_types(spot_probs, counts, nucleus_ids)

        # Convert type indices to names
        for nid, type_idx in type_assignments.items():
            all_assignments[nid] = cell_types[type_idx]

    return NucleusAssignmentResult(
        assignments=all_assignments,
        morphology_features=morph_df,
        classifier=clf,
        assignment_probs=probs_df,
    )
```

### Step 3.4: Update imports and add Optional

```python
# At top of module3b_nucleus_assignment.py, update imports:
from typing import Dict, List, Optional
```

### Step 3.5: Run test to verify it passes

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest CITEgeist/tests/test_module3b_nucleus_assignment.py::test_run_nucleus_assignment_with_cell_features -v`

Expected: PASS

### Step 3.6: Run all module3b tests

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest CITEgeist/tests/test_module3b_nucleus_assignment.py -v`

Expected: All tests PASS

### Step 3.7: Commit

```bash
git add CITEgeist/model/module3b_nucleus_assignment.py CITEgeist/tests/test_module3b_nucleus_assignment.py
git commit -m "feat: add cell_mask support to module3b for cell-level features"
```

---

## Task 4: Benchmark Script for Cell Morphology

**Files:**
- Create: `Benchmarking/xenium_benchmarking/CITEgeist/src/benchmark_cell_morphology.py`

### Step 4.1: Create benchmark script

```python
#!/usr/bin/env python
"""
Benchmark cell morphology features for single-cell assignment.

Compares nuclear-only features vs cell-level features (via watershed)
for Module 3b cell type assignment accuracy.

Usage:
    python benchmark_cell_morphology.py --region 0 --output-dir ./output/cell_morph
"""

import argparse
import json
import logging
import sys
import time
from pathlib import Path
from typing import Any, Dict, Optional, Tuple

import cv2
import numpy as np
import pandas as pd
import tifffile

# Add CITEgeist to path
REPO_ROOT = Path(__file__).parent.parent.parent.parent.parent
sys.path.insert(0, str(REPO_ROOT))

from CITEgeist.model.segmentation import run_cellpose_nuclei_segmentation
from CITEgeist.model.watershed_segmentation import watershed_from_nuclei
from CITEgeist.model.morphology_features import extract_cell_features

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)

# Paths
XENIUM_DIR = Path("/ix1/alee/LO_LAB/General/Public_Data/10x_PublicData/Xenium_RNA_Proteomic_RenalCellCarcinoma")
MORPHOLOGY_DIR = XENIUM_DIR / "morphology_focus"
PSEUDOVISIUM_DIR = REPO_ROOT / "Benchmarking/xenium_pseudovisium/data_protein_gt"


def load_boundary_channel(region_bbox: Dict[str, float]) -> np.ndarray:
    """
    Load boundary channel (ATP1A1/CD45/E-cadherin) for a region.

    Args:
        region_bbox: Dict with x_min, x_max, y_min, y_max in microns

    Returns:
        2D numpy array of boundary channel intensities
    """
    boundary_path = MORPHOLOGY_DIR / "ch0001_atp1a1_cd45_e-cadherin.ome.tif"

    # Read metadata for pixel size
    with tifffile.TiffFile(boundary_path) as tif:
        # Xenium pixel size is 0.2125 um/px
        pixel_size = 0.2125

        # Convert bbox to pixels
        x_min_px = int(region_bbox['x_min'] / pixel_size)
        x_max_px = int(region_bbox['x_max'] / pixel_size)
        y_min_px = int(region_bbox['y_min'] / pixel_size)
        y_max_px = int(region_bbox['y_max'] / pixel_size)

        # Read the region (assuming single Z-plane or MIP)
        # tifffile reads as (Y, X) for 2D
        page = tif.pages[0]
        region = page.asarray()[y_min_px:y_max_px, x_min_px:x_max_px]

    return region.astype(np.float32)


def run_benchmark(
    region_id: int,
    output_dir: Path,
    use_gpu: bool = False,
) -> Dict[str, Any]:
    """Run cell morphology benchmark for a single region."""

    region_name = f"Xenium_region_{region_id}"
    logger.info("=== Processing %s ===", region_name)

    output_dir = Path(output_dir) / region_name
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load region info
    region_dir = PSEUDOVISIUM_DIR / region_name
    coord_info_path = region_dir / "coord_info.json"

    if not coord_info_path.exists():
        raise FileNotFoundError(f"Missing coord_info.json: {coord_info_path}")

    with open(coord_info_path) as f:
        coord_info = json.load(f)

    # Load DAPI image for Cellpose
    dapi_path = MORPHOLOGY_DIR / "ch0000_dapi.ome.tif"
    logger.info("Loading DAPI and boundary channels...")

    pixel_size = 0.2125
    bbox = coord_info['micron_bounds']
    x_min_px = int(bbox['x_min'] / pixel_size)
    x_max_px = int(bbox['x_max'] / pixel_size)
    y_min_px = int(bbox['y_min'] / pixel_size)
    y_max_px = int(bbox['y_max'] / pixel_size)

    with tifffile.TiffFile(dapi_path) as tif:
        dapi_img = tif.pages[0].asarray()[y_min_px:y_max_px, x_min_px:x_max_px]

    # Load boundary channel
    boundary_img = load_boundary_channel(bbox)

    logger.info("Image shapes: DAPI=%s, boundary=%s", dapi_img.shape, boundary_img.shape)

    # Step 1: Cellpose nuclear segmentation
    logger.info("Running Cellpose nuclear segmentation...")
    t0 = time.time()

    # Convert to RGB for Cellpose (expects 3-channel)
    dapi_rgb = np.stack([dapi_img, dapi_img, dapi_img], axis=-1)
    dapi_rgb = (dapi_rgb / dapi_rgb.max() * 255).astype(np.uint8)

    nucleus_mask, centroids = run_cellpose_nuclei_segmentation(
        image_rgb_uint8=dapi_rgb,
        use_gpu=use_gpu,
        model_type="nuclei",
    )
    cellpose_time = time.time() - t0
    n_nuclei = len(centroids)
    logger.info("Cellpose: %d nuclei in %.1fs", n_nuclei, cellpose_time)

    # Step 2: Watershed cell segmentation
    logger.info("Running watershed cell segmentation...")
    t0 = time.time()
    cell_mask = watershed_from_nuclei(nucleus_mask, boundary_img)
    watershed_time = time.time() - t0
    logger.info("Watershed: completed in %.1fs", watershed_time)

    # Step 3: Extract features
    logger.info("Extracting morphology features...")
    features_df = extract_cell_features(nucleus_mask, cell_mask)
    logger.info("Extracted %d cells with %d features", len(features_df), len(features_df.columns))

    # Save outputs
    features_df.to_csv(output_dir / "cell_features.csv", index=False)
    np.save(output_dir / "nucleus_mask.npy", nucleus_mask)
    np.save(output_dir / "cell_mask.npy", cell_mask)

    results = {
        "region": region_name,
        "n_nuclei": n_nuclei,
        "n_cells_with_features": len(features_df),
        "cellpose_time_s": cellpose_time,
        "watershed_time_s": watershed_time,
        "feature_columns": list(features_df.columns),
    }

    with open(output_dir / "benchmark_results.json", "w") as f:
        json.dump(results, f, indent=2)

    return results


def main():
    parser = argparse.ArgumentParser(description="Benchmark cell morphology features")
    parser.add_argument("--region", type=int, required=True, help="Region ID (0-4)")
    parser.add_argument("--output-dir", type=str, default="./output/cell_morphology",
                        help="Output directory")
    parser.add_argument("--use-gpu", action="store_true", help="Use GPU for Cellpose")

    args = parser.parse_args()

    results = run_benchmark(
        region_id=args.region,
        output_dir=Path(args.output_dir),
        use_gpu=args.use_gpu,
    )

    logger.info("Results: %s", json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
```

### Step 4.2: Create SLURM submission script

```bash
#!/bin/bash
#SBATCH --job-name=cell_morph
#SBATCH --output=logs/cell_morph_%A_%a.out
#SBATCH --error=logs/cell_morph_%A_%a.err
#SBATCH --array=0-4
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --partition=htc
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

# Activate environment
source ~/.bashrc
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/CITEgeist/src

python benchmark_cell_morphology.py \
    --region ${SLURM_ARRAY_TASK_ID} \
    --output-dir ../output/cell_morphology
```

Save as: `Benchmarking/xenium_benchmarking/CITEgeist/slurm/sbatch_cell_morphology.sh`

### Step 4.3: Commit benchmark script

```bash
git add Benchmarking/xenium_benchmarking/CITEgeist/src/benchmark_cell_morphology.py
git add Benchmarking/xenium_benchmarking/CITEgeist/slurm/sbatch_cell_morphology.sh
git commit -m "feat: add cell morphology benchmark script"
```

---

## Task 5: Evaluation Script

**Files:**
- Create: `Benchmarking/xenium_benchmarking/evaluation/src/evaluate_cell_morphology.py`

### Step 5.1: Create evaluation script

This script compares nuclear-only vs cell-level features for classification accuracy.

```python
#!/usr/bin/env python
"""
Evaluate cell morphology features for cell type classification.

Compares:
1. Nuclear features only (current baseline)
2. Cell features (from watershed)
3. Combined nuclear + cell features

Reports accuracy, macro F1, and per-class metrics.

Usage:
    python evaluate_cell_morphology.py --region 0 --output-dir ./output/eval_cell_morph
"""

import argparse
import json
import logging
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import accuracy_score, f1_score, classification_report
from scipy.spatial import cKDTree

REPO_ROOT = Path(__file__).parent.parent.parent.parent.parent
sys.path.insert(0, str(REPO_ROOT))

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

XENIUM_DATA_DIR = Path("/ix1/alee/LO_LAB/General/Public_Data/10x_PublicData/Xenium_RNA_Proteomic_RenalCellCarcinoma")
PSEUDOVISIUM_DIR = REPO_ROOT / "Benchmarking/xenium_pseudovisium/data_protein_gt"


def load_ground_truth(region_id: int) -> pd.DataFrame:
    """Load protein-gated ground truth cell types."""
    region_name = f"Xenium_region_{region_id}"
    gt_path = PSEUDOVISIUM_DIR / region_name / "cells.csv"

    gt_df = pd.read_csv(gt_path)
    # Expected columns: cell_id, x_centroid, y_centroid, cell_type
    return gt_df


def match_cells_to_gt(
    cell_features: pd.DataFrame,
    gt_df: pd.DataFrame,
    max_dist: float = 10.0,
) -> pd.DataFrame:
    """Match extracted cells to ground truth by spatial proximity."""

    # Build KDTree from GT coordinates
    gt_coords = gt_df[['x_centroid', 'y_centroid']].values
    tree = cKDTree(gt_coords)

    # Query for each cell
    cell_coords = cell_features[['centroid_x', 'centroid_y']].values
    dists, indices = tree.query(cell_coords, k=1)

    # Filter by max distance
    valid_mask = dists <= max_dist

    matched = cell_features[valid_mask].copy()
    matched['gt_cell_type'] = gt_df.iloc[indices[valid_mask]]['cell_type'].values
    matched['match_dist'] = dists[valid_mask]

    logger.info("Matched %d/%d cells (%.1f%%) within %.1f pixels",
                valid_mask.sum(), len(cell_features),
                100 * valid_mask.mean(), max_dist)

    return matched


def evaluate_feature_set(
    df: pd.DataFrame,
    feature_cols: List[str],
    target_col: str = 'gt_cell_type',
    n_splits: int = 5,
) -> Dict[str, float]:
    """Evaluate classification accuracy using cross-validation."""

    from sklearn.model_selection import StratifiedKFold

    X = df[feature_cols].values
    y = df[target_col].values

    # Handle NaN values
    X = np.nan_to_num(X, nan=0.0, posinf=0.0, neginf=0.0)

    # Standardize features
    scaler = StandardScaler()

    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=42)

    accuracies = []
    f1_scores = []

    for train_idx, test_idx in skf.split(X, y):
        X_train, X_test = X[train_idx], X[test_idx]
        y_train, y_test = y[train_idx], y[test_idx]

        X_train_scaled = scaler.fit_transform(X_train)
        X_test_scaled = scaler.transform(X_test)

        clf = LogisticRegression(
            multi_class='multinomial',
            solver='lbfgs',
            max_iter=1000,
            random_state=42,
        )
        clf.fit(X_train_scaled, y_train)

        y_pred = clf.predict(X_test_scaled)

        accuracies.append(accuracy_score(y_test, y_pred))
        f1_scores.append(f1_score(y_test, y_pred, average='macro'))

    return {
        'accuracy_mean': np.mean(accuracies),
        'accuracy_std': np.std(accuracies),
        'f1_macro_mean': np.mean(f1_scores),
        'f1_macro_std': np.std(f1_scores),
    }


def run_evaluation(
    region_id: int,
    cell_features_path: Path,
    output_dir: Path,
) -> Dict:
    """Run full evaluation comparing feature sets."""

    region_name = f"Xenium_region_{region_id}"
    logger.info("=== Evaluating %s ===", region_name)

    output_dir = Path(output_dir) / region_name
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load cell features
    cell_features = pd.read_csv(cell_features_path)
    logger.info("Loaded %d cells with features", len(cell_features))

    # Load ground truth
    gt_df = load_ground_truth(region_id)

    # Match cells to GT
    matched = match_cells_to_gt(cell_features, gt_df)

    # Define feature sets
    nuclear_features = [
        'nucleus_area', 'nucleus_circularity', 'nucleus_eccentricity',
        'nucleus_solidity', 'nucleus_aspect_ratio',
    ]

    cell_features_cols = [
        'cell_area', 'cell_circularity', 'cell_eccentricity',
        'cell_solidity', 'cell_aspect_ratio',
    ]

    ratio_features = ['nc_ratio', 'cytoplasm_area']

    all_features = nuclear_features + cell_features_cols + ratio_features

    # Evaluate each feature set
    results = {
        'region': region_name,
        'n_cells_matched': len(matched),
        'n_cell_types': matched['gt_cell_type'].nunique(),
    }

    feature_sets = {
        'nuclear_only': nuclear_features,
        'cell_only': cell_features_cols,
        'ratio_only': ratio_features,
        'nuclear_plus_ratio': nuclear_features + ratio_features,
        'all_features': all_features,
    }

    for set_name, cols in feature_sets.items():
        # Check which columns exist
        available = [c for c in cols if c in matched.columns]
        if len(available) < len(cols):
            logger.warning("Feature set %s: only %d/%d columns available",
                          set_name, len(available), len(cols))

        if len(available) == 0:
            continue

        logger.info("Evaluating %s (%d features)...", set_name, len(available))
        metrics = evaluate_feature_set(matched, available)
        results[set_name] = {
            'features': available,
            **metrics,
        }
        logger.info("  Accuracy: %.3f +/- %.3f, F1: %.3f +/- %.3f",
                   metrics['accuracy_mean'], metrics['accuracy_std'],
                   metrics['f1_macro_mean'], metrics['f1_macro_std'])

    # Save results
    with open(output_dir / "evaluation_results.json", "w") as f:
        json.dump(results, f, indent=2)

    return results


def main():
    parser = argparse.ArgumentParser(description="Evaluate cell morphology features")
    parser.add_argument("--region", type=int, required=True, help="Region ID (0-4)")
    parser.add_argument("--features-path", type=str, required=True,
                        help="Path to cell_features.csv from benchmark")
    parser.add_argument("--output-dir", type=str, default="./output/eval_cell_morph",
                        help="Output directory")

    args = parser.parse_args()

    results = run_evaluation(
        region_id=args.region,
        cell_features_path=Path(args.features_path),
        output_dir=Path(args.output_dir),
    )

    logger.info("Evaluation complete. Results saved.")


if __name__ == "__main__":
    main()
```

### Step 5.2: Commit evaluation script

```bash
git add Benchmarking/xenium_benchmarking/evaluation/src/evaluate_cell_morphology.py
git commit -m "feat: add evaluation script for cell morphology features"
```

---

## Task 6: Update __init__.py Exports

**Files:**
- Modify: `CITEgeist/model/__init__.py`

### Step 6.1: Add new module exports

```python
# Add to CITEgeist/model/__init__.py

from .watershed_segmentation import watershed_from_nuclei
from .morphology_features import extract_cell_features
```

### Step 6.2: Commit

```bash
git add CITEgeist/model/__init__.py
git commit -m "feat: export watershed and cell feature functions from model package"
```

---

## Summary

| Task | Description | New Files | Modified Files |
|------|-------------|-----------|----------------|
| 1 | Watershed segmentation | `watershed_segmentation.py`, `test_watershed_segmentation.py` | - |
| 2 | Cell feature extraction | - | `morphology_features.py`, `test_morphology_features.py` |
| 3 | Module 3b cell_mask support | - | `module3b_nucleus_assignment.py`, `test_module3b_nucleus_assignment.py` |
| 4 | Benchmark script | `benchmark_cell_morphology.py`, `sbatch_cell_morphology.sh` | - |
| 5 | Evaluation script | `evaluate_cell_morphology.py` | - |
| 6 | Package exports | - | `__init__.py` |

**Total commits:** 6
**Estimated implementation time:** 2-3 hours

---

## Running the Full Pipeline

After implementation, run the full benchmark:

```bash
# 1. Submit benchmark jobs (extracts features + watershed for all regions)
cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/CITEgeist/slurm
sbatch sbatch_cell_morphology.sh

# 2. After completion, run evaluation for each region
cd ../evaluation/src
for i in 0 1 2 3 4; do
    python evaluate_cell_morphology.py \
        --region $i \
        --features-path ../../CITEgeist/output/cell_morphology/Xenium_region_$i/cell_features.csv \
        --output-dir ./output/eval_cell_morph
done

# 3. Compare results across regions
# Expected output: JSON files with accuracy/F1 for each feature set
```
