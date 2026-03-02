# H&E Morphology ViT+MIL Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Extend CITEgeist's morphology pipeline to support H&E images with pre-trained ViT + proportion-guided MIL for improved single-cell assignment accuracy.

**Architecture:** Pre-trained UNI ViT encoder (frozen) extracts 768-dim embeddings from 224x224 H&E nucleus patches. A proportion-guided MIL head aggregates nucleus embeddings to spot-level, trained with proportion MSE loss. Hungarian assignment converts attention weights to single-cell type predictions.

**Tech Stack:** PyTorch, timm (ViT), Cellpose, scipy (Hungarian), scikit-image

---

## Task 1: Create Visium HD Benchmarking Directory Structure

**Files:**
- Create: `Benchmarking/visiumhd_benchmarking/README.md`
- Create: `Benchmarking/visiumhd_benchmarking/src/__init__.py`
- Create: `Benchmarking/visiumhd_benchmarking/data/.gitkeep`

**Step 1: Create directory structure**

```bash
mkdir -p Benchmarking/visiumhd_benchmarking/{src,data,patches,models,results}
touch Benchmarking/visiumhd_benchmarking/src/__init__.py
touch Benchmarking/visiumhd_benchmarking/data/.gitkeep
```

**Step 2: Write README**

```markdown
# Visium HD H&E Morphology Benchmark

Benchmark morphology-based cell type assignment using pILC Visium HD data with H&E images.

## Data Sources
- pILC samples: TP08-2202, TP12-880, TP15-M509
- Location: /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/pILC_project/

## Pipeline
1. Create pseudo-Visium spots (55μm)
2. Run Cellpose segmentation on H&E
3. Extract 224x224 nucleus patches
4. Extract ViT features (UNI model)
5. Train MIL head with proportion supervision
6. Evaluate single-cell assignment accuracy
```

**Step 3: Commit**

```bash
git add Benchmarking/visiumhd_benchmarking/
git commit -m "feat: create Visium HD benchmarking directory structure"
```

---

## Task 2: Implement Pseudo-Visium Spot Creation

**Files:**
- Create: `Benchmarking/visiumhd_benchmarking/src/create_pseudo_visium.py`
- Test: `Benchmarking/visiumhd_benchmarking/src/test_create_pseudo_visium.py`

**Step 1: Write the failing test**

```python
# Benchmarking/visiumhd_benchmarking/src/test_create_pseudo_visium.py
"""Tests for pseudo-Visium spot creation."""
import numpy as np
import pandas as pd
import pytest
from create_pseudo_visium import (
    create_hex_grid,
    assign_cells_to_spots,
    compute_spot_proportions,
)


def test_create_hex_grid_basic():
    """Test hex grid generation covers bounding box."""
    bounds = (0, 1000, 0, 1000)  # x_min, x_max, y_min, y_max
    spacing_um = 100
    pixel_size = 0.5  # um/pixel

    spots = create_hex_grid(bounds, spacing_um, pixel_size)

    assert isinstance(spots, pd.DataFrame)
    assert 'spot_id' in spots.columns
    assert 'x' in spots.columns
    assert 'y' in spots.columns
    assert len(spots) > 0
    # Should have roughly (1000/100)^2 = 100 spots
    assert 50 < len(spots) < 200


def test_assign_cells_to_spots():
    """Test cell-to-spot assignment."""
    # Create mock cells
    cells = pd.DataFrame({
        'cell_id': [1, 2, 3, 4],
        'x': [50, 55, 150, 155],
        'y': [50, 55, 50, 55],
        'cell_type': ['A', 'A', 'B', 'B'],
    })

    # Create mock spots
    spots = pd.DataFrame({
        'spot_id': [0, 1],
        'x': [50, 150],
        'y': [50, 50],
    })

    mapping = assign_cells_to_spots(cells, spots, radius_um=30, pixel_size=1.0)

    assert 'spot_id' in mapping.columns
    # Cells 1,2 should be in spot 0; cells 3,4 in spot 1
    assert mapping[mapping['cell_id'] == 1]['spot_id'].iloc[0] == 0
    assert mapping[mapping['cell_id'] == 3]['spot_id'].iloc[0] == 1


def test_compute_spot_proportions():
    """Test proportion computation."""
    mapping = pd.DataFrame({
        'cell_id': [1, 2, 3, 4, 5],
        'spot_id': [0, 0, 0, 1, 1],
        'cell_type': ['A', 'A', 'B', 'B', 'B'],
    })

    proportions = compute_spot_proportions(mapping, min_cells=2)

    assert 0 in proportions.index
    assert 1 in proportions.index
    # Spot 0: 2A, 1B -> A=0.67, B=0.33
    assert abs(proportions.loc[0, 'A'] - 0.667) < 0.01
    assert abs(proportions.loc[0, 'B'] - 0.333) < 0.01
```

**Step 2: Run test to verify it fails**

```bash
cd Benchmarking/visiumhd_benchmarking/src
python -m pytest test_create_pseudo_visium.py -v
```

Expected: FAIL with "ModuleNotFoundError: No module named 'create_pseudo_visium'"

**Step 3: Write minimal implementation**

```python
# Benchmarking/visiumhd_benchmarking/src/create_pseudo_visium.py
"""Create pseudo-Visium spots from Visium HD single-cell data.

Aggregates single-cell data into Visium-like spots for benchmarking
morphology-based cell type assignment.
"""
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from typing import Tuple


def create_hex_grid(
    bounds: Tuple[float, float, float, float],
    spacing_um: float,
    pixel_size: float,
) -> pd.DataFrame:
    """Create hexagonal grid of pseudo-Visium spots.

    Args:
        bounds: (x_min, x_max, y_min, y_max) in pixels
        spacing_um: Center-to-center spacing in microns
        pixel_size: Microns per pixel

    Returns:
        DataFrame with spot_id, x, y columns
    """
    x_min, x_max, y_min, y_max = bounds
    spacing_px = spacing_um / pixel_size

    # Hex grid: offset every other row by half spacing
    spots = []
    spot_id = 0

    y = y_min + spacing_px / 2
    row = 0
    while y < y_max:
        x_offset = (spacing_px / 2) if row % 2 else 0
        x = x_min + spacing_px / 2 + x_offset
        while x < x_max:
            spots.append({'spot_id': spot_id, 'x': x, 'y': y})
            spot_id += 1
            x += spacing_px
        y += spacing_px * np.sqrt(3) / 2  # Hex row spacing
        row += 1

    return pd.DataFrame(spots)


def assign_cells_to_spots(
    cells: pd.DataFrame,
    spots: pd.DataFrame,
    radius_um: float,
    pixel_size: float,
) -> pd.DataFrame:
    """Assign cells to nearest spot within radius.

    Args:
        cells: DataFrame with cell_id, x, y, cell_type columns
        spots: DataFrame with spot_id, x, y columns
        radius_um: Maximum distance from spot center in microns
        pixel_size: Microns per pixel

    Returns:
        cells DataFrame with added spot_id column
    """
    radius_px = radius_um / pixel_size

    # Build KD-tree for spots
    spot_coords = spots[['x', 'y']].values
    tree = cKDTree(spot_coords)

    # Query nearest spot for each cell
    cell_coords = cells[['x', 'y']].values
    distances, indices = tree.query(cell_coords)

    # Assign to spot if within radius
    result = cells.copy()
    result['spot_id'] = -1  # Default: unassigned
    mask = distances <= radius_px
    result.loc[mask, 'spot_id'] = spots.iloc[indices[mask]]['spot_id'].values

    # Remove unassigned cells
    result = result[result['spot_id'] >= 0]

    return result


def compute_spot_proportions(
    mapping: pd.DataFrame,
    min_cells: int = 5,
) -> pd.DataFrame:
    """Compute cell type proportions per spot.

    Args:
        mapping: DataFrame with spot_id, cell_type columns
        min_cells: Minimum cells per spot (filter smaller spots)

    Returns:
        DataFrame with spot_id as index, cell types as columns, proportions as values
    """
    # Count cells per spot
    spot_counts = mapping.groupby('spot_id').size()
    valid_spots = spot_counts[spot_counts >= min_cells].index

    # Filter to valid spots
    filtered = mapping[mapping['spot_id'].isin(valid_spots)]

    # Compute proportions
    proportions = filtered.groupby(['spot_id', 'cell_type']).size().unstack(fill_value=0)
    proportions = proportions.div(proportions.sum(axis=1), axis=0)

    return proportions
```

**Step 4: Run test to verify it passes**

```bash
cd Benchmarking/visiumhd_benchmarking/src
python -m pytest test_create_pseudo_visium.py -v
```

Expected: PASS

**Step 5: Commit**

```bash
git add Benchmarking/visiumhd_benchmarking/src/create_pseudo_visium.py
git add Benchmarking/visiumhd_benchmarking/src/test_create_pseudo_visium.py
git commit -m "feat: implement pseudo-Visium spot creation for Visium HD"
```

---

## Task 3: Implement Cellpose H&E Segmentation Wrapper

**Files:**
- Create: `Benchmarking/visiumhd_benchmarking/src/run_cellpose_he.py`
- Test: `Benchmarking/visiumhd_benchmarking/src/test_run_cellpose_he.py`

**Step 1: Write the failing test**

```python
# Benchmarking/visiumhd_benchmarking/src/test_run_cellpose_he.py
"""Tests for Cellpose H&E segmentation."""
import numpy as np
import pytest
from run_cellpose_he import (
    preprocess_he_for_cellpose,
    segment_tile,
    stitch_masks,
)


def test_preprocess_he_for_cellpose():
    """Test H&E preprocessing for Cellpose."""
    # Mock RGB H&E image
    he_image = np.random.randint(0, 255, (100, 100, 3), dtype=np.uint8)

    processed = preprocess_he_for_cellpose(he_image)

    # Cellpose expects grayscale or specific channel format
    assert processed.ndim == 2 or processed.shape[-1] in [1, 2, 3]
    assert processed.dtype == np.float32 or processed.dtype == np.uint8


def test_segment_tile():
    """Test single tile segmentation (mock, no actual Cellpose)."""
    # This test validates the interface, actual Cellpose tested in integration
    tile = np.random.randint(0, 255, (512, 512, 3), dtype=np.uint8)

    # Mock implementation for unit test
    from unittest.mock import patch, MagicMock

    mock_model = MagicMock()
    mock_model.eval.return_value = (
        np.zeros((512, 512), dtype=np.int32),  # masks
        None,  # flows
        None,  # styles
    )

    with patch('run_cellpose_he.get_cellpose_model', return_value=mock_model):
        masks = segment_tile(tile, diameter=30)

    assert masks.shape == (512, 512)
    assert masks.dtype in [np.int32, np.int64, np.uint32]


def test_stitch_masks():
    """Test mask stitching from tiles."""
    # Two adjacent tiles with overlapping masks
    tile1 = np.array([[1, 1, 0], [1, 1, 0], [0, 0, 0]])
    tile2 = np.array([[0, 2, 2], [0, 2, 2], [0, 0, 0]])

    tiles = [(tile1, (0, 0)), (tile2, (0, 2))]  # (mask, (row_offset, col_offset))

    stitched = stitch_masks(tiles, output_shape=(3, 5), overlap=1)

    assert stitched.shape == (3, 5)
    # Should have 2 unique masks (excluding 0)
    assert len(np.unique(stitched)) - 1 <= 2
```

**Step 2: Run test to verify it fails**

```bash
python -m pytest test_run_cellpose_he.py -v
```

Expected: FAIL

**Step 3: Write minimal implementation**

```python
# Benchmarking/visiumhd_benchmarking/src/run_cellpose_he.py
"""Cellpose segmentation for H&E images.

Segments nuclei from H&E WSIs using Cellpose, processing in tiles
to handle large images.
"""
import numpy as np
from typing import List, Tuple, Optional
from pathlib import Path
import logging

logger = logging.getLogger(__name__)

# Lazy import Cellpose (heavy dependency)
_cellpose_model = None


def get_cellpose_model(model_type: str = 'nuclei', gpu: bool = True):
    """Get or create Cellpose model (singleton)."""
    global _cellpose_model
    if _cellpose_model is None:
        from cellpose import models
        _cellpose_model = models.Cellpose(model_type=model_type, gpu=gpu)
    return _cellpose_model


def preprocess_he_for_cellpose(image: np.ndarray) -> np.ndarray:
    """Preprocess H&E image for Cellpose nuclei detection.

    Cellpose nuclei model works best on images where nuclei are dark.
    For H&E, hematoxylin stains nuclei blue/purple.

    Args:
        image: (H, W, 3) RGB H&E image, uint8

    Returns:
        Preprocessed image for Cellpose (grayscale, nuclei dark)
    """
    # Convert to grayscale - nuclei appear darker
    if image.ndim == 3:
        # Use blue channel (hematoxylin) or luminance
        gray = 0.299 * image[..., 0] + 0.587 * image[..., 1] + 0.114 * image[..., 2]
        gray = gray.astype(np.float32)
    else:
        gray = image.astype(np.float32)

    # Invert so nuclei are bright (Cellpose expects bright objects)
    gray = 255 - gray

    # Normalize to [0, 255]
    gray = (gray - gray.min()) / (gray.max() - gray.min() + 1e-8) * 255

    return gray.astype(np.uint8)


def segment_tile(
    tile: np.ndarray,
    diameter: float = 30,
    model_type: str = 'nuclei',
    gpu: bool = True,
) -> np.ndarray:
    """Segment nuclei in a single tile.

    Args:
        tile: (H, W, 3) RGB tile
        diameter: Expected nucleus diameter in pixels
        model_type: Cellpose model type
        gpu: Use GPU if available

    Returns:
        (H, W) integer mask with unique ID per nucleus
    """
    model = get_cellpose_model(model_type, gpu)

    # Preprocess
    processed = preprocess_he_for_cellpose(tile)

    # Run Cellpose - channels=[0,0] for grayscale
    masks, flows, styles = model.eval(
        processed,
        diameter=diameter,
        channels=[0, 0],
        flow_threshold=0.4,
        cellprob_threshold=0.0,
    )

    return masks.astype(np.int32)


def stitch_masks(
    tiles: List[Tuple[np.ndarray, Tuple[int, int]]],
    output_shape: Tuple[int, int],
    overlap: int = 0,
) -> np.ndarray:
    """Stitch tile masks into full image mask.

    Args:
        tiles: List of (mask, (row_offset, col_offset)) tuples
        output_shape: (H, W) of full output
        overlap: Overlap between tiles in pixels

    Returns:
        Stitched mask with relabeled unique IDs
    """
    output = np.zeros(output_shape, dtype=np.int32)
    max_label = 0

    for mask, (row_off, col_off) in tiles:
        h, w = mask.shape

        # Relabel mask to avoid conflicts
        relabeled = mask.copy()
        relabeled[relabeled > 0] += max_label

        # Determine non-overlap region (center of tile)
        row_start = row_off + overlap // 2 if row_off > 0 else row_off
        col_start = col_off + overlap // 2 if col_off > 0 else col_off

        # Place in output
        row_end = min(row_off + h, output_shape[0])
        col_end = min(col_off + w, output_shape[1])

        mask_row_start = row_start - row_off
        mask_col_start = col_start - col_off

        region = output[row_start:row_end, col_start:col_end]
        mask_region = relabeled[mask_row_start:mask_row_start + (row_end - row_start),
                                mask_col_start:mask_col_start + (col_end - col_start)]

        # Only overwrite zeros
        region[region == 0] = mask_region[region == 0]

        max_label = output.max()

    return output


def segment_wsi(
    wsi_path: Path,
    output_path: Path,
    tile_size: int = 2048,
    overlap: int = 128,
    diameter: float = 30,
    gpu: bool = True,
) -> np.ndarray:
    """Segment entire WSI in tiles.

    Args:
        wsi_path: Path to WSI TIFF
        output_path: Path to save mask .npy
        tile_size: Tile size in pixels
        overlap: Overlap between tiles
        diameter: Expected nucleus diameter
        gpu: Use GPU

    Returns:
        Full segmentation mask
    """
    from PIL import Image
    from tqdm import tqdm

    Image.MAX_IMAGE_PIXELS = None

    logger.info(f"Loading WSI from {wsi_path}")
    wsi = np.array(Image.open(wsi_path))
    h, w = wsi.shape[:2]

    logger.info(f"WSI shape: {wsi.shape}")

    # Process tiles
    tiles = []
    step = tile_size - overlap

    positions = [(r, c) for r in range(0, h, step) for c in range(0, w, step)]

    for row, col in tqdm(positions, desc="Segmenting tiles"):
        row_end = min(row + tile_size, h)
        col_end = min(col + tile_size, w)

        tile = wsi[row:row_end, col:col_end]

        mask = segment_tile(tile, diameter=diameter, gpu=gpu)
        tiles.append((mask, (row, col)))

    # Stitch
    logger.info("Stitching masks")
    full_mask = stitch_masks(tiles, (h, w), overlap=overlap)

    # Save
    np.save(output_path, full_mask)
    logger.info(f"Saved mask to {output_path}, {full_mask.max()} nuclei detected")

    return full_mask
```

**Step 4: Run test to verify it passes**

```bash
python -m pytest test_run_cellpose_he.py -v
```

Expected: PASS

**Step 5: Commit**

```bash
git add Benchmarking/visiumhd_benchmarking/src/run_cellpose_he.py
git add Benchmarking/visiumhd_benchmarking/src/test_run_cellpose_he.py
git commit -m "feat: implement Cellpose H&E segmentation wrapper"
```

---

## Task 4: Implement H&E Patch Extraction (224x224)

**Files:**
- Create: `Benchmarking/visiumhd_benchmarking/src/extract_patches_he.py`
- Test: `Benchmarking/visiumhd_benchmarking/src/test_extract_patches_he.py`

**Step 1: Write the failing test**

```python
# Benchmarking/visiumhd_benchmarking/src/test_extract_patches_he.py
"""Tests for H&E patch extraction."""
import numpy as np
import pytest
from extract_patches_he import (
    extract_nucleus_patch,
    normalize_he_patch,
    extract_patches_for_spot,
)


def test_extract_nucleus_patch_basic():
    """Test basic patch extraction."""
    # Mock image and mask
    image = np.random.randint(0, 255, (500, 500, 3), dtype=np.uint8)
    mask = np.zeros((500, 500), dtype=np.int32)
    mask[100:120, 100:120] = 1  # Nucleus at (100-120, 100-120)

    patch = extract_nucleus_patch(
        image, mask, nucleus_id=1,
        output_size=224, expansion=0.75
    )

    assert patch.shape == (224, 224, 3)
    assert patch.dtype == np.uint8


def test_extract_nucleus_patch_edge():
    """Test patch extraction at image edge."""
    image = np.random.randint(0, 255, (100, 100, 3), dtype=np.uint8)
    mask = np.zeros((100, 100), dtype=np.int32)
    mask[0:10, 0:10] = 1  # Nucleus at corner

    patch = extract_nucleus_patch(
        image, mask, nucleus_id=1,
        output_size=224, expansion=0.75
    )

    assert patch.shape == (224, 224, 3)


def test_normalize_he_patch():
    """Test H&E normalization."""
    patch = np.random.randint(0, 255, (224, 224, 3), dtype=np.uint8)

    normalized = normalize_he_patch(patch)

    assert normalized.shape == (3, 224, 224)  # CHW format
    assert normalized.dtype == np.float32
    assert normalized.min() >= 0 and normalized.max() <= 1


def test_extract_patches_for_spot():
    """Test batch patch extraction for a spot."""
    image = np.random.randint(0, 255, (500, 500, 3), dtype=np.uint8)
    mask = np.zeros((500, 500), dtype=np.int32)
    mask[100:120, 100:120] = 1
    mask[200:220, 200:220] = 2
    mask[300:320, 300:320] = 3

    nucleus_ids = np.array([1, 2, 3])

    patches, valid_ids = extract_patches_for_spot(
        image, mask, nucleus_ids,
        output_size=224, expansion=0.75
    )

    assert patches.shape == (3, 3, 224, 224)  # (N, C, H, W)
    assert len(valid_ids) == 3
```

**Step 2: Run test to verify it fails**

```bash
python -m pytest test_extract_patches_he.py -v
```

Expected: FAIL

**Step 3: Write minimal implementation**

```python
# Benchmarking/visiumhd_benchmarking/src/extract_patches_he.py
"""Extract nucleus patches from H&E images.

Extracts 224x224 RGB patches around detected nuclei for ViT input.
"""
import numpy as np
from typing import Tuple, Optional
from skimage.measure import regionprops
from skimage.transform import resize


def extract_nucleus_patch(
    image: np.ndarray,
    mask: np.ndarray,
    nucleus_id: int,
    output_size: int = 224,
    expansion: float = 0.75,
) -> np.ndarray:
    """Extract a patch around a nucleus.

    Args:
        image: (H, W, 3) RGB image
        mask: (H, W) segmentation mask with nucleus IDs
        nucleus_id: ID of nucleus to extract
        output_size: Output patch size (square)
        expansion: Fraction to expand bounding box

    Returns:
        (output_size, output_size, 3) RGB patch
    """
    # Get nucleus bounding box
    nucleus_mask = mask == nucleus_id
    if not nucleus_mask.any():
        raise ValueError(f"Nucleus {nucleus_id} not found in mask")

    rows = np.where(nucleus_mask.any(axis=1))[0]
    cols = np.where(nucleus_mask.any(axis=0))[0]

    r_min, r_max = rows.min(), rows.max()
    c_min, c_max = cols.min(), cols.max()

    # Expand bounding box
    height = r_max - r_min
    width = c_max - c_min

    r_expand = int(height * expansion / 2)
    c_expand = int(width * expansion / 2)

    r_min = max(0, r_min - r_expand)
    r_max = min(image.shape[0], r_max + r_expand)
    c_min = max(0, c_min - c_expand)
    c_max = min(image.shape[1], c_max + c_expand)

    # Extract and resize
    crop = image[r_min:r_max, c_min:c_max]

    # Resize to output_size
    resized = resize(crop, (output_size, output_size),
                     preserve_range=True, anti_aliasing=True)

    return resized.astype(np.uint8)


def normalize_he_patch(patch: np.ndarray) -> np.ndarray:
    """Normalize H&E patch for ViT input.

    Args:
        patch: (H, W, 3) RGB patch, uint8

    Returns:
        (3, H, W) normalized float32 tensor
    """
    # Convert to float [0, 1]
    normalized = patch.astype(np.float32) / 255.0

    # HWC to CHW
    normalized = np.transpose(normalized, (2, 0, 1))

    return normalized


def extract_patches_for_spot(
    image: np.ndarray,
    mask: np.ndarray,
    nucleus_ids: np.ndarray,
    output_size: int = 224,
    expansion: float = 0.75,
) -> Tuple[np.ndarray, np.ndarray]:
    """Extract patches for all nuclei in a spot.

    Args:
        image: (H, W, 3) RGB image
        mask: (H, W) segmentation mask
        nucleus_ids: IDs of nuclei in this spot
        output_size: Output patch size
        expansion: Bounding box expansion

    Returns:
        patches: (N, 3, H, W) normalized patches
        valid_ids: IDs of successfully extracted nuclei
    """
    patches = []
    valid_ids = []

    for nid in nucleus_ids:
        try:
            patch = extract_nucleus_patch(
                image, mask, int(nid),
                output_size=output_size,
                expansion=expansion
            )
            normalized = normalize_he_patch(patch)
            patches.append(normalized)
            valid_ids.append(nid)
        except (ValueError, IndexError):
            continue

    if not patches:
        return np.empty((0, 3, output_size, output_size), dtype=np.float32), np.array([])

    return np.stack(patches), np.array(valid_ids)
```

**Step 4: Run test to verify it passes**

```bash
python -m pytest test_extract_patches_he.py -v
```

Expected: PASS

**Step 5: Commit**

```bash
git add Benchmarking/visiumhd_benchmarking/src/extract_patches_he.py
git add Benchmarking/visiumhd_benchmarking/src/test_extract_patches_he.py
git commit -m "feat: implement H&E patch extraction (224x224)"
```

---

## Task 5: Implement ViT Feature Extractor

**Files:**
- Create: `CITEgeist/model/vit_extractor.py`
- Test: `CITEgeist/tests/test_vit_extractor.py`

**Step 1: Write the failing test**

```python
# CITEgeist/tests/test_vit_extractor.py
"""Tests for ViT feature extraction."""
import numpy as np
import torch
import pytest


def test_vit_extractor_shape():
    """Test ViT output shape."""
    from CITEgeist.model.vit_extractor import ViTFeatureExtractor

    extractor = ViTFeatureExtractor(model_name='vit_small_patch16_224', pretrained=False)

    # Mock batch
    x = torch.randn(4, 3, 224, 224)

    with torch.no_grad():
        features = extractor(x)

    assert features.shape == (4, 384)  # vit_small has 384-dim output


def test_vit_extractor_numpy():
    """Test numpy input handling."""
    from CITEgeist.model.vit_extractor import ViTFeatureExtractor

    extractor = ViTFeatureExtractor(model_name='vit_small_patch16_224', pretrained=False)

    # Numpy input
    x = np.random.rand(2, 3, 224, 224).astype(np.float32)

    features = extractor.extract_numpy(x)

    assert isinstance(features, np.ndarray)
    assert features.shape == (2, 384)


def test_vit_extractor_imagenet_normalization():
    """Test ImageNet normalization is applied."""
    from CITEgeist.model.vit_extractor import ViTFeatureExtractor

    extractor = ViTFeatureExtractor(model_name='vit_small_patch16_224', pretrained=False)

    # Input in [0, 1]
    x = torch.ones(1, 3, 224, 224) * 0.5

    # After normalization, values should be different
    normalized = extractor.normalize(x)

    # ImageNet mean ~ 0.485, so (0.5 - 0.485) / 0.229 ~ 0.065
    assert not torch.allclose(normalized, x)
```

**Step 2: Run test to verify it fails**

```bash
cd CITEgeist
python -m pytest tests/test_vit_extractor.py -v
```

Expected: FAIL

**Step 3: Write minimal implementation**

```python
# CITEgeist/model/vit_extractor.py
"""Vision Transformer feature extractor for histopathology images.

Provides pre-trained ViT encoders for extracting features from H&E patches.
Supports multiple model variants including UNI foundation model.
"""
import torch
import torch.nn as nn
import numpy as np
from typing import Optional, Union
from pathlib import Path


class ViTFeatureExtractor(nn.Module):
    """ViT-based feature extractor for nucleus patches.

    Wraps timm ViT models with ImageNet normalization and optional
    pretrained weights from histopathology foundation models.

    Attributes:
        model: timm ViT model (classification head removed)
        embed_dim: Output embedding dimension
    """

    # ImageNet normalization constants
    IMAGENET_MEAN = (0.485, 0.456, 0.406)
    IMAGENET_STD = (0.229, 0.224, 0.225)

    def __init__(
        self,
        model_name: str = 'vit_large_patch16_224',
        pretrained: bool = True,
        weights_path: Optional[Path] = None,
        device: str = 'cuda' if torch.cuda.is_available() else 'cpu',
    ):
        """Initialize ViT feature extractor.

        Args:
            model_name: timm model name (e.g., 'vit_large_patch16_224')
            pretrained: Load ImageNet pretrained weights
            weights_path: Path to custom weights (e.g., UNI model)
            device: Device to run on
        """
        super().__init__()
        import timm

        # Create model without classification head
        self.model = timm.create_model(
            model_name,
            pretrained=pretrained,
            num_classes=0,  # Remove classification head
        )

        # Load custom weights if provided
        if weights_path is not None:
            state_dict = torch.load(weights_path, map_location='cpu')
            # Handle different state dict formats
            if 'model' in state_dict:
                state_dict = state_dict['model']
            self.model.load_state_dict(state_dict, strict=False)

        self.model = self.model.to(device)
        self.model.eval()

        self.device = device
        self.embed_dim = self.model.num_features

        # Register normalization as buffer
        self.register_buffer(
            'mean',
            torch.tensor(self.IMAGENET_MEAN).view(1, 3, 1, 1)
        )
        self.register_buffer(
            'std',
            torch.tensor(self.IMAGENET_STD).view(1, 3, 1, 1)
        )

    def normalize(self, x: torch.Tensor) -> torch.Tensor:
        """Apply ImageNet normalization.

        Args:
            x: (B, 3, H, W) input in [0, 1]

        Returns:
            Normalized tensor
        """
        return (x - self.mean) / self.std

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """Extract features from patches.

        Args:
            x: (B, 3, 224, 224) input patches in [0, 1]

        Returns:
            (B, embed_dim) feature vectors
        """
        x = self.normalize(x)
        return self.model(x)

    @torch.no_grad()
    def extract_numpy(
        self,
        patches: np.ndarray,
        batch_size: int = 32,
    ) -> np.ndarray:
        """Extract features from numpy array of patches.

        Args:
            patches: (N, 3, 224, 224) float32 in [0, 1]
            batch_size: Batch size for processing

        Returns:
            (N, embed_dim) feature array
        """
        self.eval()

        features = []
        n_patches = len(patches)

        for i in range(0, n_patches, batch_size):
            batch = patches[i:i + batch_size]
            batch_tensor = torch.from_numpy(batch).to(self.device)
            batch_features = self.forward(batch_tensor)
            features.append(batch_features.cpu().numpy())

        return np.concatenate(features, axis=0)


def load_uni_extractor(
    weights_path: Path,
    device: str = 'cuda',
) -> ViTFeatureExtractor:
    """Load UNI foundation model.

    UNI is a ViT-L/16 trained on 100M+ histopathology patches.

    Args:
        weights_path: Path to UNI weights
        device: Device to run on

    Returns:
        ViTFeatureExtractor with UNI weights
    """
    return ViTFeatureExtractor(
        model_name='vit_large_patch16_224',
        pretrained=False,
        weights_path=weights_path,
        device=device,
    )
```

**Step 4: Run test to verify it passes**

```bash
python -m pytest tests/test_vit_extractor.py -v
```

Expected: PASS

**Step 5: Commit**

```bash
git add CITEgeist/model/vit_extractor.py
git add CITEgeist/tests/test_vit_extractor.py
git commit -m "feat: implement ViT feature extractor with UNI support"
```

---

## Task 6: Implement Proportion-Guided MIL Module

**Files:**
- Create: `CITEgeist/model/proportion_mil.py`
- Test: `CITEgeist/tests/test_proportion_mil.py`

**Step 1: Write the failing test**

```python
# CITEgeist/tests/test_proportion_mil.py
"""Tests for proportion-guided MIL."""
import torch
import pytest


def test_proportion_mil_forward():
    """Test MIL forward pass."""
    from CITEgeist.model.proportion_mil import ProportionGuidedMIL

    model = ProportionGuidedMIL(input_dim=768, n_cell_types=9, hidden_dim=256)

    # Mock embeddings for 15 nuclei
    embeddings = torch.randn(15, 768)

    proportions, attention = model(embeddings)

    assert proportions.shape == (9,)
    assert attention.shape == (15, 9)
    assert torch.allclose(proportions.sum(), torch.tensor(1.0), atol=1e-5)
    assert torch.allclose(attention.sum(dim=0), torch.ones(9), atol=1e-5)


def test_proportion_mil_loss():
    """Test proportion loss computation."""
    from CITEgeist.model.proportion_mil import ProportionGuidedMIL, proportion_loss

    pred = torch.tensor([0.3, 0.5, 0.2])
    target = torch.tensor([0.25, 0.55, 0.2])

    loss = proportion_loss(pred, target)

    assert loss.ndim == 0  # Scalar
    assert loss >= 0


def test_proportion_mil_training_step():
    """Test full training step."""
    from CITEgeist.model.proportion_mil import ProportionGuidedMIL, proportion_loss

    model = ProportionGuidedMIL(input_dim=768, n_cell_types=5, hidden_dim=128)
    optimizer = torch.optim.Adam(model.parameters(), lr=1e-3)

    # Mock data
    embeddings = torch.randn(20, 768)
    gt_proportions = torch.tensor([0.2, 0.3, 0.25, 0.15, 0.1])

    # Training step
    model.train()
    pred_proportions, attention = model(embeddings)
    loss = proportion_loss(pred_proportions, gt_proportions)

    optimizer.zero_grad()
    loss.backward()
    optimizer.step()

    assert loss.item() > 0
```

**Step 2: Run test to verify it fails**

```bash
python -m pytest tests/test_proportion_mil.py -v
```

Expected: FAIL

**Step 3: Write minimal implementation**

```python
# CITEgeist/model/proportion_mil.py
"""Proportion-guided Multiple Instance Learning for cell type prediction.

Aggregates nucleus-level embeddings to spot-level proportions using
attention mechanism. Trained with proportion supervision.
"""
import torch
import torch.nn as nn
import torch.nn.functional as F
from typing import Tuple


class ProportionGuidedMIL(nn.Module):
    """MIL aggregator with proportion guidance.

    Uses gated attention to learn per-cell-type attention weights,
    enabling interpretable aggregation of nucleus features to
    spot-level cell type proportions.

    Attributes:
        n_cell_types: Number of cell types (K)
        attention_V: Value transformation for gated attention
        attention_U: Gate transformation for gated attention
        attention_W: Attention logit projection
    """

    def __init__(
        self,
        input_dim: int = 768,
        n_cell_types: int = 9,
        hidden_dim: int = 256,
        dropout: float = 0.1,
    ):
        """Initialize MIL module.

        Args:
            input_dim: Dimension of input embeddings
            n_cell_types: Number of cell types to predict
            hidden_dim: Hidden dimension in attention network
            dropout: Dropout rate
        """
        super().__init__()
        self.n_cell_types = n_cell_types

        # Gated attention mechanism
        self.attention_V = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.Tanh(),
            nn.Dropout(dropout),
        )
        self.attention_U = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.Sigmoid(),
            nn.Dropout(dropout),
        )
        self.attention_W = nn.Linear(hidden_dim, n_cell_types)

        # Optional: classifier for direct proportion prediction
        self.classifier = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, n_cell_types),
        )

    def forward(
        self,
        embeddings: torch.Tensor,
    ) -> Tuple[torch.Tensor, torch.Tensor]:
        """Forward pass.

        Args:
            embeddings: (N, input_dim) nucleus embeddings

        Returns:
            proportions: (K,) predicted cell type proportions
            attention: (N, K) attention weights per nucleus per type
        """
        N = embeddings.shape[0]

        # Compute gated attention
        A_V = self.attention_V(embeddings)  # (N, hidden)
        A_U = self.attention_U(embeddings)  # (N, hidden)
        A_gate = A_V * A_U  # (N, hidden)

        # Attention logits per cell type
        A_logits = self.attention_W(A_gate)  # (N, K)

        # Softmax over nuclei for each cell type
        attention = F.softmax(A_logits, dim=0)  # (N, K)

        # Weighted aggregation: (K, N) @ (N, D) -> (K, D)
        weighted_embeddings = torch.mm(attention.T, embeddings)  # (K, input_dim)

        # Classify each weighted embedding
        type_logits = self.classifier(weighted_embeddings)  # (K, K)

        # Diagonal gives prediction for each type's aggregated embedding
        # But simpler: just use attention weights directly as proportions
        proportions = attention.sum(dim=0) / N  # (K,)

        # Ensure sums to 1
        proportions = proportions / (proportions.sum() + 1e-8)

        return proportions, attention

    def get_nucleus_probabilities(
        self,
        attention: torch.Tensor,
    ) -> torch.Tensor:
        """Convert attention to per-nucleus type probabilities.

        Args:
            attention: (N, K) attention weights from forward()

        Returns:
            (N, K) probability of each type for each nucleus
        """
        # Normalize per nucleus (row-wise)
        return F.softmax(attention, dim=1)


def proportion_loss(
    pred: torch.Tensor,
    target: torch.Tensor,
    eps: float = 1e-8,
) -> torch.Tensor:
    """Compute proportion prediction loss.

    Combines MSE loss with KL divergence for better gradient flow.

    Args:
        pred: (K,) predicted proportions
        target: (K,) ground truth proportions
        eps: Small constant for numerical stability

    Returns:
        Scalar loss
    """
    # MSE loss
    mse = F.mse_loss(pred, target)

    # KL divergence (pred || target)
    kl = (pred * torch.log((pred + eps) / (target + eps))).sum()

    return mse + 0.1 * kl


def entropy_regularization(attention: torch.Tensor) -> torch.Tensor:
    """Compute entropy of attention weights for regularization.

    Higher entropy = more uniform attention = less confident.
    Can be used to encourage or discourage concentration.

    Args:
        attention: (N, K) attention weights

    Returns:
        Scalar entropy (mean across cell types)
    """
    eps = 1e-8
    # Entropy per cell type (column)
    entropy_per_type = -(attention * torch.log(attention + eps)).sum(dim=0)
    return entropy_per_type.mean()
```

**Step 4: Run test to verify it passes**

```bash
python -m pytest tests/test_proportion_mil.py -v
```

Expected: PASS

**Step 5: Commit**

```bash
git add CITEgeist/model/proportion_mil.py
git add CITEgeist/tests/test_proportion_mil.py
git commit -m "feat: implement proportion-guided MIL module"
```

---

## Task 7: Implement Training Pipeline

**Files:**
- Create: `Benchmarking/visiumhd_benchmarking/src/train_mil.py`
- Test: `Benchmarking/visiumhd_benchmarking/src/test_train_mil.py`

**Step 1: Write the failing test**

```python
# Benchmarking/visiumhd_benchmarking/src/test_train_mil.py
"""Tests for MIL training pipeline."""
import numpy as np
import torch
import pytest
from pathlib import Path
from train_mil import (
    SpotDataset,
    train_epoch,
    evaluate,
)


def test_spot_dataset():
    """Test spot dataset loading."""
    # Create mock data directory
    import tempfile
    import json

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)

        # Create mock spot files
        for i in range(3):
            spot_dir = tmpdir / f"spot_{i}"
            spot_dir.mkdir()

            # Mock embeddings and proportions
            np.save(spot_dir / "embeddings.npy", np.random.randn(10, 768).astype(np.float32))
            np.save(spot_dir / "proportions.npy", np.random.dirichlet(np.ones(5)).astype(np.float32))

        dataset = SpotDataset(tmpdir, n_cell_types=5)

        assert len(dataset) == 3

        embeddings, proportions = dataset[0]
        assert embeddings.shape[1] == 768
        assert proportions.shape == (5,)


def test_train_epoch():
    """Test single training epoch."""
    from CITEgeist.model.proportion_mil import ProportionGuidedMIL

    model = ProportionGuidedMIL(input_dim=768, n_cell_types=5, hidden_dim=64)
    optimizer = torch.optim.Adam(model.parameters(), lr=1e-3)

    # Mock dataloader
    mock_data = [
        (torch.randn(10, 768), torch.tensor([0.2, 0.3, 0.25, 0.15, 0.1])),
        (torch.randn(15, 768), torch.tensor([0.1, 0.4, 0.2, 0.2, 0.1])),
    ]

    loss = train_epoch(model, mock_data, optimizer)

    assert loss > 0


def test_evaluate():
    """Test evaluation."""
    from CITEgeist.model.proportion_mil import ProportionGuidedMIL

    model = ProportionGuidedMIL(input_dim=768, n_cell_types=5, hidden_dim=64)
    model.eval()

    mock_data = [
        (torch.randn(10, 768), torch.tensor([0.2, 0.3, 0.25, 0.15, 0.1])),
    ]

    metrics = evaluate(model, mock_data)

    assert 'loss' in metrics
    assert 'pearson_r' in metrics
```

**Step 2: Run test to verify it fails**

```bash
python -m pytest test_train_mil.py -v
```

Expected: FAIL

**Step 3: Write minimal implementation**

```python
# Benchmarking/visiumhd_benchmarking/src/train_mil.py
"""Training pipeline for proportion-guided MIL.

Trains MIL head on pre-extracted ViT embeddings with spot-level
proportion supervision.
"""
import sys
from pathlib import Path
import numpy as np
import torch
import torch.nn.functional as F
from torch.utils.data import Dataset, DataLoader
from typing import List, Tuple, Dict, Optional
import logging
from scipy.stats import pearsonr

# Add repo root
REPO_ROOT = Path(__file__).parent.parent.parent.parent
sys.path.insert(0, str(REPO_ROOT))

from CITEgeist.model.proportion_mil import ProportionGuidedMIL, proportion_loss, entropy_regularization

logger = logging.getLogger(__name__)


class SpotDataset(Dataset):
    """Dataset of spot embeddings and proportions.

    Expects directory structure:
        data_dir/
            spot_0/
                embeddings.npy  # (N, embed_dim)
                proportions.npy  # (K,)
            spot_1/
            ...
    """

    def __init__(
        self,
        data_dir: Path,
        n_cell_types: int,
        min_nuclei: int = 5,
    ):
        """Initialize dataset.

        Args:
            data_dir: Directory containing spot subdirectories
            n_cell_types: Number of cell types
            min_nuclei: Minimum nuclei per spot
        """
        self.data_dir = Path(data_dir)
        self.n_cell_types = n_cell_types

        # Find valid spots
        self.spots = []
        for spot_dir in sorted(self.data_dir.glob("spot_*")):
            emb_path = spot_dir / "embeddings.npy"
            prop_path = spot_dir / "proportions.npy"

            if emb_path.exists() and prop_path.exists():
                emb = np.load(emb_path)
                if len(emb) >= min_nuclei:
                    self.spots.append(spot_dir)

    def __len__(self) -> int:
        return len(self.spots)

    def __getitem__(self, idx: int) -> Tuple[torch.Tensor, torch.Tensor]:
        spot_dir = self.spots[idx]

        embeddings = np.load(spot_dir / "embeddings.npy")
        proportions = np.load(spot_dir / "proportions.npy")

        return (
            torch.from_numpy(embeddings).float(),
            torch.from_numpy(proportions).float(),
        )


def collate_spots(batch: List[Tuple[torch.Tensor, torch.Tensor]]):
    """Custom collate - each spot is processed independently."""
    return batch  # Return as list, not stacked


def train_epoch(
    model: ProportionGuidedMIL,
    data: List[Tuple[torch.Tensor, torch.Tensor]],
    optimizer: torch.optim.Optimizer,
    device: str = 'cpu',
    entropy_weight: float = 0.01,
) -> float:
    """Train for one epoch.

    Args:
        model: MIL model
        data: List of (embeddings, proportions) tuples
        optimizer: Optimizer
        device: Device
        entropy_weight: Weight for entropy regularization

    Returns:
        Average loss
    """
    model.train()
    total_loss = 0.0

    for embeddings, gt_proportions in data:
        embeddings = embeddings.to(device)
        gt_proportions = gt_proportions.to(device)

        # Forward
        pred_proportions, attention = model(embeddings)

        # Loss
        loss = proportion_loss(pred_proportions, gt_proportions)
        loss += entropy_weight * entropy_regularization(attention)

        # Backward
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()

        total_loss += loss.item()

    return total_loss / len(data)


@torch.no_grad()
def evaluate(
    model: ProportionGuidedMIL,
    data: List[Tuple[torch.Tensor, torch.Tensor]],
    device: str = 'cpu',
) -> Dict[str, float]:
    """Evaluate model.

    Args:
        model: MIL model
        data: List of (embeddings, proportions) tuples
        device: Device

    Returns:
        Dictionary of metrics
    """
    model.eval()

    total_loss = 0.0
    all_pred = []
    all_gt = []

    for embeddings, gt_proportions in data:
        embeddings = embeddings.to(device)
        gt_proportions = gt_proportions.to(device)

        pred_proportions, _ = model(embeddings)

        loss = proportion_loss(pred_proportions, gt_proportions)
        total_loss += loss.item()

        all_pred.append(pred_proportions.cpu().numpy())
        all_gt.append(gt_proportions.cpu().numpy())

    # Compute correlation
    all_pred = np.concatenate([p.flatten() for p in all_pred])
    all_gt = np.concatenate([g.flatten() for g in all_gt])

    r, p = pearsonr(all_pred, all_gt)

    return {
        'loss': total_loss / len(data),
        'pearson_r': r,
        'pearson_p': p,
    }


def train(
    model: ProportionGuidedMIL,
    train_data: List[Tuple[torch.Tensor, torch.Tensor]],
    val_data: List[Tuple[torch.Tensor, torch.Tensor]],
    n_epochs: int = 100,
    lr: float = 1e-4,
    device: str = 'cpu',
    save_path: Optional[Path] = None,
) -> Dict[str, List[float]]:
    """Full training loop.

    Args:
        model: MIL model
        train_data: Training data
        val_data: Validation data
        n_epochs: Number of epochs
        lr: Learning rate
        device: Device
        save_path: Path to save best model

    Returns:
        Training history
    """
    model = model.to(device)
    optimizer = torch.optim.AdamW(model.parameters(), lr=lr, weight_decay=1e-4)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, n_epochs)

    history = {'train_loss': [], 'val_loss': [], 'val_r': []}
    best_r = -1

    for epoch in range(n_epochs):
        train_loss = train_epoch(model, train_data, optimizer, device)
        val_metrics = evaluate(model, val_data, device)

        scheduler.step()

        history['train_loss'].append(train_loss)
        history['val_loss'].append(val_metrics['loss'])
        history['val_r'].append(val_metrics['pearson_r'])

        if val_metrics['pearson_r'] > best_r:
            best_r = val_metrics['pearson_r']
            if save_path:
                torch.save(model.state_dict(), save_path)

        if epoch % 10 == 0:
            logger.info(
                f"Epoch {epoch}: train_loss={train_loss:.4f}, "
                f"val_loss={val_metrics['loss']:.4f}, val_r={val_metrics['pearson_r']:.4f}"
            )

    return history
```

**Step 4: Run test to verify it passes**

```bash
python -m pytest test_train_mil.py -v
```

Expected: PASS

**Step 5: Commit**

```bash
git add Benchmarking/visiumhd_benchmarking/src/train_mil.py
git add Benchmarking/visiumhd_benchmarking/src/test_train_mil.py
git commit -m "feat: implement MIL training pipeline"
```

---

## Task 8: Implement Evaluation Pipeline

**Files:**
- Create: `Benchmarking/visiumhd_benchmarking/src/evaluate_single_cell.py`
- Test: `Benchmarking/visiumhd_benchmarking/src/test_evaluate_single_cell.py`

**Step 1: Write the failing test**

```python
# Benchmarking/visiumhd_benchmarking/src/test_evaluate_single_cell.py
"""Tests for single-cell evaluation."""
import numpy as np
import torch
import pytest
from evaluate_single_cell import (
    evaluate_spot_assignment,
    compute_accuracy_metrics,
)


def test_evaluate_spot_assignment():
    """Test single spot evaluation."""
    # Mock attention weights (5 nuclei, 3 types)
    attention = np.array([
        [0.8, 0.1, 0.1],  # Likely type 0
        [0.7, 0.2, 0.1],  # Likely type 0
        [0.1, 0.8, 0.1],  # Likely type 1
        [0.2, 0.1, 0.7],  # Likely type 2
        [0.1, 0.1, 0.8],  # Likely type 2
    ])

    gt_proportions = np.array([0.4, 0.2, 0.4])  # 2, 1, 2
    gt_types = np.array([0, 0, 1, 2, 2])

    results = evaluate_spot_assignment(
        attention, gt_proportions, gt_types,
        use_hungarian=True
    )

    assert 'accuracy' in results
    assert 'assignments' in results
    assert 0 <= results['accuracy'] <= 1


def test_compute_accuracy_metrics():
    """Test accuracy metric computation."""
    predictions = np.array([0, 0, 1, 2, 2])
    ground_truth = np.array([0, 0, 1, 2, 1])  # One mistake

    metrics = compute_accuracy_metrics(predictions, ground_truth, n_types=3)

    assert 'overall_accuracy' in metrics
    assert 'per_type_f1' in metrics
    assert metrics['overall_accuracy'] == 0.8  # 4/5
```

**Step 2: Run test to verify it fails**

```bash
python -m pytest test_evaluate_single_cell.py -v
```

Expected: FAIL

**Step 3: Write minimal implementation**

```python
# Benchmarking/visiumhd_benchmarking/src/evaluate_single_cell.py
"""Single-cell assignment evaluation.

Evaluates MIL attention weights against ground truth cell types
using Hungarian assignment for constrained evaluation.
"""
import sys
from pathlib import Path
import numpy as np
from typing import Dict, Optional
from sklearn.metrics import accuracy_score, f1_score, confusion_matrix

REPO_ROOT = Path(__file__).parent.parent.parent.parent
sys.path.insert(0, str(REPO_ROOT))

from CITEgeist.model.hungarian_assignment import assign_nuclei_to_types


def evaluate_spot_assignment(
    attention: np.ndarray,
    gt_proportions: np.ndarray,
    gt_types: np.ndarray,
    use_hungarian: bool = True,
) -> Dict:
    """Evaluate single-cell assignment for one spot.

    Args:
        attention: (N, K) attention weights from MIL
        gt_proportions: (K,) ground truth proportions
        gt_types: (N,) ground truth cell type indices
        use_hungarian: Use Hungarian assignment with count constraints

    Returns:
        Dictionary with accuracy and assignments
    """
    n_nuclei, n_types = attention.shape

    if use_hungarian:
        # Convert proportions to counts
        counts = np.round(gt_proportions * n_nuclei).astype(int)

        # Adjust rounding errors
        diff = n_nuclei - counts.sum()
        if diff > 0:
            counts[np.argmax(gt_proportions)] += diff
        elif diff < 0:
            counts[np.argmax(counts)] += diff

        # Run Hungarian
        nucleus_ids = np.arange(n_nuclei)
        assignments_dict = assign_nuclei_to_types(attention, counts, nucleus_ids)
        assignments = np.array([assignments_dict[i] for i in range(n_nuclei)])
    else:
        # Simple argmax
        assignments = np.argmax(attention, axis=1)

    # Compute accuracy
    accuracy = accuracy_score(gt_types, assignments)

    return {
        'accuracy': accuracy,
        'assignments': assignments,
        'gt_types': gt_types,
    }


def compute_accuracy_metrics(
    predictions: np.ndarray,
    ground_truth: np.ndarray,
    n_types: int,
    type_names: Optional[list] = None,
) -> Dict:
    """Compute comprehensive accuracy metrics.

    Args:
        predictions: (N,) predicted cell type indices
        ground_truth: (N,) ground truth cell type indices
        n_types: Number of cell types
        type_names: Optional names for cell types

    Returns:
        Dictionary of metrics
    """
    if type_names is None:
        type_names = [f"Type_{i}" for i in range(n_types)]

    overall_acc = accuracy_score(ground_truth, predictions)

    # Per-type F1
    f1_per_type = f1_score(
        ground_truth, predictions,
        labels=list(range(n_types)),
        average=None,
        zero_division=0,
    )

    # Confusion matrix
    cm = confusion_matrix(ground_truth, predictions, labels=list(range(n_types)))

    # Per-type accuracy (recall)
    per_type_acc = {}
    for i, name in enumerate(type_names):
        mask = ground_truth == i
        if mask.sum() > 0:
            per_type_acc[name] = accuracy_score(ground_truth[mask], predictions[mask])
        else:
            per_type_acc[name] = np.nan

    return {
        'overall_accuracy': overall_acc,
        'macro_f1': f1_score(ground_truth, predictions, average='macro', zero_division=0),
        'weighted_f1': f1_score(ground_truth, predictions, average='weighted', zero_division=0),
        'per_type_f1': dict(zip(type_names, f1_per_type)),
        'per_type_accuracy': per_type_acc,
        'confusion_matrix': cm,
    }


def evaluate_all_spots(
    model,
    data_loader,
    gt_types_per_spot: Dict[int, np.ndarray],
    device: str = 'cpu',
) -> Dict:
    """Evaluate across all spots.

    Args:
        model: Trained MIL model
        data_loader: DataLoader yielding (embeddings, proportions)
        gt_types_per_spot: Dict mapping spot_id to ground truth types
        device: Device

    Returns:
        Aggregated metrics
    """
    import torch

    model.eval()

    all_predictions = []
    all_gt = []
    spot_accuracies = []

    with torch.no_grad():
        for spot_id, (embeddings, gt_proportions) in enumerate(data_loader):
            embeddings = embeddings.to(device)

            _, attention = model(embeddings)
            attention_np = attention.cpu().numpy()

            gt_types = gt_types_per_spot.get(spot_id)
            if gt_types is None:
                continue

            results = evaluate_spot_assignment(
                attention_np,
                gt_proportions.numpy(),
                gt_types,
                use_hungarian=True,
            )

            all_predictions.extend(results['assignments'])
            all_gt.extend(gt_types)
            spot_accuracies.append(results['accuracy'])

    all_predictions = np.array(all_predictions)
    all_gt = np.array(all_gt)

    n_types = len(np.unique(all_gt))
    metrics = compute_accuracy_metrics(all_predictions, all_gt, n_types)
    metrics['mean_spot_accuracy'] = np.mean(spot_accuracies)

    return metrics
```

**Step 4: Run test to verify it passes**

```bash
python -m pytest test_evaluate_single_cell.py -v
```

Expected: PASS

**Step 5: Commit**

```bash
git add Benchmarking/visiumhd_benchmarking/src/evaluate_single_cell.py
git add Benchmarking/visiumhd_benchmarking/src/test_evaluate_single_cell.py
git commit -m "feat: implement single-cell evaluation pipeline"
```

---

## Task 9: Create Main Pipeline Script

**Files:**
- Create: `Benchmarking/visiumhd_benchmarking/src/run_benchmark.py`

**Step 1: Write the main script**

```python
# Benchmarking/visiumhd_benchmarking/src/run_benchmark.py
"""Main benchmark script for Visium HD H&E morphology.

End-to-end pipeline:
1. Create pseudo-Visium spots from Visium HD
2. Run Cellpose segmentation
3. Extract patches and ViT features
4. Train MIL with proportion supervision
5. Evaluate single-cell assignment
"""
import argparse
import logging
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import torch
from tqdm import tqdm

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

REPO_ROOT = Path(__file__).parent.parent.parent.parent
sys.path.insert(0, str(REPO_ROOT))


def create_pseudo_visium(
    sample_path: Path,
    output_dir: Path,
    spot_diameter_um: float = 55,
    spot_spacing_um: float = 100,
    pixel_size: float = 0.5,
    min_cells: int = 5,
):
    """Step 1: Create pseudo-Visium spots."""
    from create_pseudo_visium import create_hex_grid, assign_cells_to_spots, compute_spot_proportions
    import scanpy as sc

    logger.info(f"Loading sample from {sample_path}")
    adata = sc.read_h5ad(sample_path)

    # Get cell coordinates and types
    spatial = adata.obsm['spatial']
    if hasattr(spatial, 'values'):
        spatial = spatial.values

    cells = pd.DataFrame({
        'cell_id': range(len(adata)),
        'x': spatial[:, 0],
        'y': spatial[:, 1],
        'cell_type': adata.obs['cell_type_canonical'].values,
    })

    # Filter Unknown
    cells = cells[cells['cell_type'] != 'Unknown']
    logger.info(f"Filtered to {len(cells)} cells (excluded Unknown)")

    # Create grid
    bounds = (cells['x'].min(), cells['x'].max(), cells['y'].min(), cells['y'].max())
    spots = create_hex_grid(bounds, spot_spacing_um, pixel_size)
    logger.info(f"Created {len(spots)} spots")

    # Assign cells to spots
    mapping = assign_cells_to_spots(cells, spots, spot_diameter_um / 2, pixel_size)
    logger.info(f"Assigned {len(mapping)} cells to spots")

    # Compute proportions
    proportions = compute_spot_proportions(mapping, min_cells=min_cells)
    logger.info(f"Valid spots with >={min_cells} cells: {len(proportions)}")

    # Save
    output_dir.mkdir(parents=True, exist_ok=True)
    spots.to_csv(output_dir / "spots.csv", index=False)
    mapping.to_csv(output_dir / "cell_to_spot_mapping.csv", index=False)
    proportions.to_parquet(output_dir / "proportions.parquet")

    return spots, mapping, proportions


def run_cellpose(
    wsi_path: Path,
    output_path: Path,
    diameter: float = 30,
    gpu: bool = True,
):
    """Step 2: Run Cellpose segmentation."""
    from run_cellpose_he import segment_wsi

    logger.info(f"Running Cellpose on {wsi_path}")
    mask = segment_wsi(wsi_path, output_path, diameter=diameter, gpu=gpu)
    logger.info(f"Detected {mask.max()} nuclei")

    return mask


def extract_features(
    wsi_path: Path,
    mask_path: Path,
    mapping_path: Path,
    output_dir: Path,
    vit_weights: Path,
    batch_size: int = 32,
    device: str = 'cuda',
):
    """Step 3: Extract patches and ViT features."""
    from PIL import Image
    from extract_patches_he import extract_patches_for_spot
    from CITEgeist.model.vit_extractor import ViTFeatureExtractor

    Image.MAX_IMAGE_PIXELS = None

    logger.info("Loading data")
    wsi = np.array(Image.open(wsi_path))
    mask = np.load(mask_path)
    mapping = pd.read_csv(mapping_path)

    # Initialize ViT
    logger.info("Loading ViT model")
    vit = ViTFeatureExtractor(weights_path=vit_weights, device=device)

    # Process each spot
    output_dir.mkdir(parents=True, exist_ok=True)

    for spot_id, group in tqdm(mapping.groupby('spot_id'), desc="Extracting features"):
        spot_dir = output_dir / f"spot_{spot_id}"
        spot_dir.mkdir(exist_ok=True)

        nucleus_ids = group['cell_id'].values

        # Extract patches
        patches, valid_ids = extract_patches_for_spot(
            wsi, mask, nucleus_ids,
            output_size=224, expansion=0.75
        )

        if len(patches) == 0:
            continue

        # Extract ViT features
        embeddings = vit.extract_numpy(patches, batch_size=batch_size)

        # Save
        np.save(spot_dir / "embeddings.npy", embeddings)
        np.save(spot_dir / "nucleus_ids.npy", valid_ids)
        np.save(spot_dir / "gt_types.npy", group.set_index('cell_id').loc[valid_ids, 'cell_type'].values)


def train_and_evaluate(
    features_dir: Path,
    proportions_path: Path,
    output_dir: Path,
    n_cell_types: int,
    n_epochs: int = 100,
    device: str = 'cuda',
):
    """Step 4-5: Train MIL and evaluate."""
    from train_mil import SpotDataset, train, evaluate
    from evaluate_single_cell import evaluate_all_spots, compute_accuracy_metrics
    from CITEgeist.model.proportion_mil import ProportionGuidedMIL

    # Load proportions
    proportions = pd.read_parquet(proportions_path)

    # Create datasets
    dataset = SpotDataset(features_dir, n_cell_types=n_cell_types)

    # Split train/val
    n_train = int(0.8 * len(dataset))
    train_indices = list(range(n_train))
    val_indices = list(range(n_train, len(dataset)))

    train_data = [dataset[i] for i in train_indices]
    val_data = [dataset[i] for i in val_indices]

    logger.info(f"Train spots: {len(train_data)}, Val spots: {len(val_data)}")

    # Train
    model = ProportionGuidedMIL(input_dim=768, n_cell_types=n_cell_types)

    history = train(
        model, train_data, val_data,
        n_epochs=n_epochs,
        device=device,
        save_path=output_dir / "best_model.pt",
    )

    # Save history
    with open(output_dir / "training_history.json", 'w') as f:
        json.dump({k: [float(v) for v in vals] for k, vals in history.items()}, f)

    # Load best model and evaluate
    model.load_state_dict(torch.load(output_dir / "best_model.pt"))

    # Evaluate on validation
    val_metrics = evaluate(model, val_data, device=device)

    logger.info(f"Validation Pearson r: {val_metrics['pearson_r']:.4f}")

    # Save results
    with open(output_dir / "results.json", 'w') as f:
        json.dump(val_metrics, f, indent=2)

    return model, val_metrics


def main():
    parser = argparse.ArgumentParser(description="Visium HD H&E Morphology Benchmark")
    parser.add_argument("--sample", type=Path, required=True, help="Path to sample h5ad")
    parser.add_argument("--wsi", type=Path, required=True, help="Path to H&E WSI")
    parser.add_argument("--vit-weights", type=Path, required=True, help="Path to ViT weights")
    parser.add_argument("--output", type=Path, required=True, help="Output directory")
    parser.add_argument("--device", type=str, default="cuda", help="Device")
    parser.add_argument("--n-epochs", type=int, default=100, help="Training epochs")

    args = parser.parse_args()

    output_dir = args.output
    output_dir.mkdir(parents=True, exist_ok=True)

    # Step 1: Create pseudo-Visium
    logger.info("=== Step 1: Create pseudo-Visium ===")
    spots, mapping, proportions = create_pseudo_visium(
        args.sample,
        output_dir / "pseudo_visium",
    )

    # Step 2: Run Cellpose
    logger.info("=== Step 2: Run Cellpose ===")
    mask_path = output_dir / "cellpose_mask.npy"
    if not mask_path.exists():
        run_cellpose(args.wsi, mask_path)

    # Step 3: Extract features
    logger.info("=== Step 3: Extract ViT features ===")
    extract_features(
        args.wsi,
        mask_path,
        output_dir / "pseudo_visium" / "cell_to_spot_mapping.csv",
        output_dir / "features",
        args.vit_weights,
        device=args.device,
    )

    # Step 4-5: Train and evaluate
    logger.info("=== Step 4-5: Train MIL and evaluate ===")
    n_cell_types = len(proportions.columns)
    model, metrics = train_and_evaluate(
        output_dir / "features",
        output_dir / "pseudo_visium" / "proportions.parquet",
        output_dir / "model",
        n_cell_types=n_cell_types,
        n_epochs=args.n_epochs,
        device=args.device,
    )

    logger.info(f"=== Complete! Results saved to {output_dir} ===")


if __name__ == "__main__":
    main()
```

**Step 2: Commit**

```bash
git add Benchmarking/visiumhd_benchmarking/src/run_benchmark.py
git commit -m "feat: implement main benchmark pipeline script"
```

---

## Task 10: Create SLURM Submission Script

**Files:**
- Create: `Benchmarking/visiumhd_benchmarking/slurm/sbatch_benchmark.sh`

**Step 1: Write SLURM script**

```bash
#!/bin/bash
#SBATCH --job-name=visiumhd_benchmark
#SBATCH --output=logs/benchmark_%j.out
#SBATCH --error=logs/benchmark_%j.err
#SBATCH --time=48:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --partition=l40s
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

# Load modules
module load gurobi/12.0.3

# Activate environment
source activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

# Set paths
REPO_ROOT=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist
PILC_DIR=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/pILC_project
ENACT_DIR=/ihome/alee/alc376/alc376_bgfs/enact-pipeline/output_files

# Sample to process (override via environment variable)
SAMPLE=${SAMPLE:-TP08-2202}

# Map sample to paths
case $SAMPLE in
    TP08-2202)
        H5AD=$PILC_DIR/enact_adatas/TP08-2202_canonical_annotated.h5ad
        WSI=$ENACT_DIR/tp08-2202-pilc/tmap/wsi.tif
        ;;
    TP12-880)
        H5AD=$PILC_DIR/enact_adatas/TP12-880_canonical_annotated.h5ad
        WSI=$ENACT_DIR/tp12-880-pilc/tmap/wsi.tif
        ;;
    TP15-M509)
        H5AD=$PILC_DIR/enact_adatas/TP15-M509_canonical_annotated.h5ad
        WSI=$ENACT_DIR/tp15-m509-pilc/tmap/wsi.tif
        ;;
    *)
        echo "Unknown sample: $SAMPLE"
        exit 1
        ;;
esac

OUTPUT=$REPO_ROOT/Benchmarking/visiumhd_benchmarking/results/$SAMPLE

# ViT weights (UNI model - needs to be downloaded)
VIT_WEIGHTS=$REPO_ROOT/models/uni_weights.pth

# Create output directory
mkdir -p $OUTPUT
mkdir -p $REPO_ROOT/Benchmarking/visiumhd_benchmarking/logs

# Run benchmark
cd $REPO_ROOT/Benchmarking/visiumhd_benchmarking/src

python run_benchmark.py \
    --sample $H5AD \
    --wsi $WSI \
    --vit-weights $VIT_WEIGHTS \
    --output $OUTPUT \
    --device cuda \
    --n-epochs 100

echo "Benchmark complete for $SAMPLE"
```

**Step 2: Commit**

```bash
mkdir -p Benchmarking/visiumhd_benchmarking/slurm
chmod +x Benchmarking/visiumhd_benchmarking/slurm/sbatch_benchmark.sh
git add Benchmarking/visiumhd_benchmarking/slurm/sbatch_benchmark.sh
git commit -m "feat: add SLURM submission script for Visium HD benchmark"
```

---

## Task 11: Update Model __init__.py

**Files:**
- Modify: `CITEgeist/model/__init__.py`

**Step 1: Add new modules to exports**

```python
# Add to CITEgeist/model/__init__.py

from .vit_extractor import ViTFeatureExtractor, load_uni_extractor
from .proportion_mil import ProportionGuidedMIL, proportion_loss, entropy_regularization
```

**Step 2: Commit**

```bash
git add CITEgeist/model/__init__.py
git commit -m "feat: export ViT extractor and MIL modules"
```

---

## Summary

| Task | Component | Files | Status |
|------|-----------|-------|--------|
| 1 | Directory structure | `visiumhd_benchmarking/` | Pending |
| 2 | Pseudo-Visium creation | `create_pseudo_visium.py` | Pending |
| 3 | Cellpose H&E wrapper | `run_cellpose_he.py` | Pending |
| 4 | H&E patch extraction | `extract_patches_he.py` | Pending |
| 5 | ViT feature extractor | `vit_extractor.py` | Pending |
| 6 | Proportion-guided MIL | `proportion_mil.py` | Pending |
| 7 | Training pipeline | `train_mil.py` | Pending |
| 8 | Evaluation pipeline | `evaluate_single_cell.py` | Pending |
| 9 | Main benchmark script | `run_benchmark.py` | Pending |
| 10 | SLURM submission | `sbatch_benchmark.sh` | Pending |
| 11 | Module exports | `__init__.py` | Pending |

**Total estimated time:** 4-6 hours

**Dependencies to install:**
- `timm` (PyTorch Image Models)
- UNI weights download from Hugging Face

**Note:** UNI model weights require Hugging Face access. Alternative: use `vit_large_patch16_224` with ImageNet weights for initial testing.
