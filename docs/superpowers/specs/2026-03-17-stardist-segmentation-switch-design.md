# Design Spec: StarDist Segmentation Switch

**Date:** 2026-03-17
**Status:** Approved
**Branch:** `dev`

## Problem

Cellpose `nuclei` model is trained on fluorescence microscopy (bright nuclei on dark background). On H&E histology (dark nuclei on light background), it segments tissue structures and holes instead of nuclei. The unified PC-MIL pipeline (12 patient Visium samples) produced bad segmentation masks, leading to downstream attention collapse.

## Decision

Replace Cellpose with StarDist across the entire codebase. StarDist has purpose-built pretrained models for both H&E and fluorescence:

- **`2D_versatile_he`** — trained on H&E histopathology nuclei
- **`2D_versatile_fluo`** — trained on fluorescence (DAPI) nuclei

## Scope

### Files Modified (Code Changes)

| File | Change |
|------|--------|
| `CITEgeist/model/segmentation.py` | Replace `_build_cellpose_model()`, `run_cellpose_nuclei_segmentation()`, `_segment_fullres_by_spot_patches()`, `compute_spot_nuclei_counts_from_adata()`, `save_segmentation_artifacts()` with `StarDistSegmenter` class and `run_nuclei_segmentation()` convenience function. Update `save_segmentation_artifacts()` output filenames from `*_cellpose_*` → `*_stardist_*`. Update `_prepare_rgb_uint8()` docstring. |
| `CITEgeist/model/citegeist_model.py` | Rename `compute_spot_nuclei_counts_cellpose()` → `compute_spot_nuclei_counts()`. Update internals to use `StarDistSegmenter`. |
| `CITEgeist/examples/run_unified_step2_features.py` | Remove `_get_cellpose_model()`. Use `StarDistSegmenter(modality)`. Wire up DAPI path. Rename output dirs `cellpose/` → `segmentation/`. |
| `CITEgeist/examples/run_unified_step1_module3.py` | Update call from `compute_spot_nuclei_counts_cellpose()` → `compute_spot_nuclei_counts()`. |
| `CITEgeist/examples/run_unified_step3_pcmil.py` | Update path from `cellpose/` → `segmentation/` for centroids CSV. |
| `Benchmarking/visiumhd_benchmarking/src/run_cellpose_he.py` | Replace Cellpose with `StarDistSegmenter("he")`. Rename file to `run_segmentation_he.py`. |
| `Benchmarking/visiumhd_benchmarking/src/run_benchmark.py` | Replace `step2_run_cellpose()` with StarDist equivalent. Remove `--cellpose-diameter` and `--skip-cellpose` CLI args, replace with `--skip-segmentation`. Update mask filenames from `cellpose_mask.npy` → `stardist_mask.npy`. |
| `Benchmarking/visiumhd_benchmarking/src/test_run_cellpose_he.py` | Rename to `test_run_segmentation_he.py`. Update all imports and mock targets to reference new module name. |
| `Benchmarking/visiumhd_benchmarking/src/visualize_he_patches.py` | Update mask path from `cellpose_mask.npy` → `stardist_mask.npy`. Update labels. |
| `Benchmarking/visiumhd_benchmarking/src/__init__.py` | Update module docstring. |
| `Benchmarking/xenium_benchmarking/CITEgeist/src/prepare_patches.py` | Add optional `--run-segmentation` flag using `StarDistSegmenter("dapi")`. Default: load existing masks (unchanged behavior). |
| `Benchmarking/xenium_benchmarking/CITEgeist/src/benchmark_hybrid_cellpose.py` | Replace `run_cellpose_nuclei_segmentation` import with `run_nuclei_segmentation`. Update centroid handling for DataFrame return type. |
| `Benchmarking/xenium_benchmarking/CITEgeist/src/benchmark_multimodal_cellpose.py` | Same as above. |
| `Benchmarking/xenium_benchmarking/CITEgeist/src/benchmark_discrete_cellpose.py` | Replace Cellpose imports (`SegmentationResult`, `run_cellpose_nuclei_segmentation`). Update centroid handling. |
| `Benchmarking/xenium_benchmarking/CITEgeist/src/benchmark_single_cell_resolution.py` | Replace `run_cellpose_nuclei_segmentation` with `run_nuclei_segmentation`. Update `centroids_xy` to use DataFrame. |
| `Benchmarking/xenium_benchmarking/CITEgeist/src/benchmark_cell_morphology.py` | Replace `run_cellpose_nuclei_segmentation` import. |
| `Benchmarking/simulation_benchmarking/CITEgeist/src/generate_cellpose_images.py` | Update docstrings/comments referencing Cellpose models. |
| `Benchmarking/simulation_benchmarking/CITEgeist/src/test_discretize_continuous.py` | Update `run_cellpose_nuclei_segmentation` import → `run_nuclei_segmentation`. |
| `CITEgeist/tests/test_segmentation.py` | Update monkeypatched Cellpose functions to StarDist equivalents. |
| `examples/run_module3_unified.py` | Update `compute_spot_nuclei_counts_cellpose()` → `compute_spot_nuclei_counts()`. |
| `pyproject.toml` | Replace `cellpose>=3.0.0` with `stardist` in dependencies. |

### Backward Compatibility

The `run_nuclei_segmentation()` convenience function is NOT a true drop-in for `run_cellpose_nuclei_segmentation()` — the return type changes from `(masks, centroids_xy_array)` to `(masks, centroids_df)`. All callers must be updated.

To ease migration of Xenium benchmark scripts (which may be re-run for comparison), keep a thin deprecation wrapper:

```python
def run_cellpose_nuclei_segmentation(image, **kwargs):
    """DEPRECATED: Use run_nuclei_segmentation() instead."""
    import warnings
    warnings.warn("run_cellpose_nuclei_segmentation() is deprecated, use run_nuclei_segmentation()", DeprecationWarning)
    masks, centroids_df = run_nuclei_segmentation(image, modality="dapi", **kwargs)
    # Convert to legacy (x, y) array format
    centroids_xy = centroids_df[["x_pixel", "y_pixel"]].values
    return masks, centroids_xy
```

### TF GPU Fallback

If TensorFlow cannot access GPU (common on HPC login nodes or CUDA version mismatches):
- StarDist falls back to CPU automatically
- Set `TF_FORCE_GPU_ALLOW_GROWTH=true` in SLURM scripts
- If TF GPU init fails, log warning and continue on CPU

### Files with Comment/Docstring Updates Only

| File | Change |
|------|--------|
| `CITEgeist/model/gurobi_impl.py` | Update docstrings referencing "Cellpose" → "nuclei segmentation" |
| `CITEgeist/model/detection_estimation.py` | Update docstrings |
| `CITEgeist/model/morphology_features.py` | Update module docstring and `_prepare_rgb_uint8()` docstring |
| `CITEgeist/model/module3b_nucleus_assignment.py` | Update docstring |
| `CITEgeist/model/two_stage_pipeline.py` | Update comment |
| `Benchmarking/visiumhd_benchmarking/README.md` | Update references to `run_cellpose_he.py` → `run_segmentation_he.py`, `cellpose_masks.npy` → `stardist_mask.npy`, cellpose dependency → stardist |
| `Benchmarking/visiumhd_benchmarking/slurm/sbatch_benchmark*.sh` | Update `--skip-cellpose` → `--skip-segmentation` in comments |
| `Benchmarking/xenium_benchmarking/CITEgeist/src/unsupervised_kmeans_classifier.py` | Update comment referencing Cellpose pixels |
| `CLAUDE.md` | Update Module 3 Alternative description and `segmentation.py` references |

### Environment

- Install `stardist` into `CITEgeist_env` (brings `tensorflow` as dependency)
- Replace `cellpose>=3.0.0` with `stardist` in `pyproject.toml`

## Architecture

### StarDistSegmenter Class

Located in `CITEgeist/model/segmentation.py`:

```python
class StarDistSegmenter:
    """Nuclei segmentation via StarDist pretrained models.

    Args:
        modality: "he" for H&E histology, "dapi" for fluorescence.
    """

    MODELS = {
        "he": "2D_versatile_he",
        "dapi": "2D_versatile_fluo",
    }

    def __init__(self, modality: str = "he"):
        from stardist.models import StarDist2D
        model_name = self.MODELS[modality]
        self.model = StarDist2D.from_pretrained(model_name)
        self.modality = modality

    def segment(self, image: np.ndarray, **kwargs) -> tuple:
        """Segment nuclei and return (masks, centroids_df).

        Args:
            image: For H&E: (H, W, 3) uint8 RGB.
                   For DAPI: (H, W) single-channel float or uint16.

        Returns:
            masks: (H, W) int32 label array, 0 = background.
            centroids_df: DataFrame with columns [y_pixel, x_pixel, nucleus_id].
        """
        from csbdeep.utils import normalize

        # Normalize input
        if self.modality == "he":
            # StarDist H&E model expects normalized RGB
            img_norm = normalize(image, 1, 99.8, axis=(0, 1))
        else:
            img_norm = normalize(image, 1, 99.8)

        labels, details = self.model.predict_instances(img_norm, **kwargs)

        # Extract centroids
        from scipy.ndimage import center_of_mass
        nucleus_ids = np.unique(labels)
        nucleus_ids = nucleus_ids[nucleus_ids > 0]
        centroids = center_of_mass(labels > 0, labels, nucleus_ids)
        centroids_df = pd.DataFrame(centroids, columns=["y_pixel", "x_pixel"])
        centroids_df["nucleus_id"] = nucleus_ids.astype(int)

        return labels.astype(np.int32), centroids_df
```

### Centroid Coordinate Convention

**IMPORTANT:** The old `run_cellpose_nuclei_segmentation()` returned centroids as `(x, y)` numpy arrays. The new `StarDistSegmenter.segment()` returns a DataFrame with `[y_pixel, x_pixel]` columns.

All downstream consumers must be updated:
- `run_unified_step2_features.py` already swaps to `(x, y)` via `centroids_df[["x_pixel", "y_pixel"]].values` at line 149 — **no change needed**.
- `compute_spot_nuclei_counts_from_adata()` in `segmentation.py` — update to use DataFrame columns instead of `centroids_xy` array.
- `assign_nuclei_centroids_to_spots()` in `segmentation.py` — update to accept DataFrame or convert internally.

### Convenience Function

```python
def run_nuclei_segmentation(image, modality="he", **kwargs):
    """Drop-in replacement for run_cellpose_nuclei_segmentation().

    Returns:
        masks: (H, W) int32 label array.
        centroids_df: DataFrame with [y_pixel, x_pixel, nucleus_id].
    """
    segmenter = StarDistSegmenter(modality=modality)
    return segmenter.segment(image, **kwargs)
```

### CitegeistModel Integration

```python
# Old:
def compute_spot_nuclei_counts_cellpose(self, ...):
    ...
    masks, centroids_xy = run_cellpose_nuclei_segmentation(image, ...)

# New:
def compute_spot_nuclei_counts(self, ...):
    ...
    segmenter = StarDistSegmenter(modality="he")
    masks, centroids_df = segmenter.segment(image, ...)
```

### Unified Pipeline Step 2

```python
# Old:
cp_model = _get_cellpose_model("nuclei", torch.cuda.is_available())
masks, _, _, _ = cp_model.eval(fullres_crop, channels=[0, 0], diameter=None)

# New:
from model.segmentation import StarDistSegmenter
segmenter = StarDistSegmenter(modality=modality)  # "he" or "dapi"
masks, centroids_df = segmenter.segment(fullres_crop)
```

### DAPI Path (New)

Wire up the `else` branch in `load_image_and_segment()`:

```python
elif modality == "dapi":
    # Load Xenium OME-TIFF, extract DAPI channel only (channel 0)
    # Use StarDistSegmenter("dapi") for segmentation
    # Return masks, centroids, image, spatial_coords, barcodes, spot_radius
```

Xenium image loading extracted from `prepare_patches.py` into a shared helper.

### Parallelization Note

The existing `_segment_fullres_by_spot_patches()` uses multi-threaded Cellpose with `max_workers`. StarDist's TensorFlow backend does not support the same thread-per-worker pattern. Replace with:
- Single StarDist model instance, process image in one call (StarDist handles large images well)
- Or tile the image and process sequentially if memory is an issue

## Output Directory Rename

The `cellpose/` output directory in the unified pipeline is renamed to `segmentation/`:

```
output/unified_pipeline/{sample_name}/
├── module3/          # Step 1 (unchanged)
├── segmentation/     # Step 2 (was: cellpose/)
│   ├── nuclei_masks.npy
│   ├── nuclei_centroids.csv
│   └── nuclei_per_spot.csv
├── features/         # Step 2 (unchanged)
├── pcmil/            # Step 3 (unchanged)
└── validation/       # Step 4 (unchanged)
```

## Cache Invalidation

- New mask cache files: `*_stardist_masks_tissue.npy`
- Glob pattern in Step 2 updated from `*_cellpose_masks_tissue.npy` → `*_stardist_masks_tissue.npy`
- Old `*_cellpose_masks_tissue.npy` files are ignored (different glob pattern)
- `.step2_complete` markers must be cleared before re-running
- Old `cellpose/` output dirs are not deleted automatically; users can clean up manually

## GPU Memory Management

StarDist uses TensorFlow; ViT feature extraction uses PyTorch. Both run in the same Step 2 script sequentially.

**Critical:** TensorFlow pre-allocates all GPU memory by default. Must set memory growth BEFORE any TF import:

```python
import tensorflow as tf
tf.config.experimental.set_memory_growth(
    tf.config.list_physical_devices('GPU')[0], True
)
```

Or set environment variable in SLURM scripts:
```bash
export TF_FORCE_GPU_ALLOW_GROWTH=true
```

Sequence:
1. Set TF memory growth
2. StarDist segmentation (TF)
3. `tf.keras.backend.clear_session()` + `del segmenter`
4. ViT feature extraction (PyTorch)

## Error Handling

- **No detections**: Log warning, return empty masks/centroids. Existing Step 2 logic handles empty patches gracefully.
- **TF import failure**: Raise clear error with install instructions (`pip install stardist`).
- **Image format**: `StarDistSegmenter.segment()` validates input shape and dtype, converts if needed.

## Testing

- Verify StarDist `2D_versatile_he` produces reasonable nuclei masks on one patient H&E sample
- Verify `2D_versatile_fluo` works on Xenium DAPI channel
- Re-run unified pipeline Steps 2-4 for all 12 patient samples
- Compare nuclei counts to Cellpose baseline (expect more, smaller nuclei with StarDist on H&E)
- Run updated Visium HD benchmark test suite (`test_run_segmentation_he.py`)

## Not in Scope

- Retraining or fine-tuning StarDist models
- Changing patch size or ViT model
- Fixing PC-MIL attention collapse (separate issue)
