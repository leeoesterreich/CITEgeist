# StarDist Segmentation Switch Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace Cellpose with StarDist for nuclei segmentation across the entire CITEgeist codebase — `2D_versatile_he` for H&E, `2D_versatile_fluo` for DAPI.

**Architecture:** Single `StarDistSegmenter` class in `segmentation.py` with modality-based model selection. Deprecation wrapper preserves old API for gradual migration. TF memory growth required for GPU coexistence with PyTorch.

**Tech Stack:** StarDist, TensorFlow, csbdeep, scipy, pandas

**Spec:** `docs/superpowers/specs/2026-03-17-stardist-segmentation-switch-design.md`

**Worktree:** Execute in a git worktree branched from `dev` (feature/stardist-segmentation)

---

## Task 1: Environment Setup

**Files:**
- Modify: `pyproject.toml`

- [ ] **Step 1: Install stardist into CITEgeist_env**

```bash
source activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env
pip install stardist
```

- [ ] **Step 2: Verify stardist imports and pretrained models download**

```bash
python -c "
from stardist.models import StarDist2D
from csbdeep.utils import normalize
m = StarDist2D.from_pretrained('2D_versatile_he')
print('H&E model loaded:', m)
m2 = StarDist2D.from_pretrained('2D_versatile_fluo')
print('Fluo model loaded:', m2)
print('SUCCESS')
"
```

Expected: Both models download and load without error.

- [ ] **Step 3: Update pyproject.toml**

In `pyproject.toml`, find the `cellpose>=3.0.0` dependency and replace:

```python
# Old:
"cellpose>=3.0,<4.0",
# New:
"stardist",
```

- [ ] **Step 4: Commit**

```bash
git add pyproject.toml
git commit -m "feat: replace cellpose with stardist in dependencies"
```

---

## Task 2: Core StarDistSegmenter Class

**Files:**
- Modify: `CITEgeist/model/segmentation.py`
- Modify: `CITEgeist/tests/test_segmentation.py`

This is the core change. Add `StarDistSegmenter`, `run_nuclei_segmentation()`, and a deprecation wrapper for `run_cellpose_nuclei_segmentation()`. Keep existing helper functions (`assign_nuclei_centroids_to_spots`, `_prepare_rgb_uint8`, etc.) that don't reference Cellpose directly.

- [ ] **Step 1: Read current segmentation.py to understand all Cellpose-dependent functions**

Read the full file. Identify:
- `_build_cellpose_model()` — DELETE
- `run_cellpose_nuclei_segmentation()` — REPLACE with deprecation wrapper
- `_segment_fullres_by_spot_patches()` — DELETE (replaced by single StarDist call)
- `compute_spot_nuclei_counts_from_adata()` — UPDATE to use StarDistSegmenter
- `save_segmentation_artifacts()` — UPDATE filenames from `cellpose` → `stardist`
- `_prepare_rgb_uint8()` — KEEP, update docstring only

- [ ] **Step 2: Add StarDistSegmenter class at top of segmentation.py**

Add after imports, before any existing functions. See spec for full class code. Key points:
- Lazy import of `stardist.models.StarDist2D` inside `__init__`
- Lazy import of `csbdeep.utils.normalize` inside `segment()`
- Returns `(masks: np.ndarray, centroids_df: pd.DataFrame)` with columns `[y_pixel, x_pixel, nucleus_id]`
- H&E: `normalize(image, 1, 99.8, axis=(0, 1))` for RGB
- DAPI: `normalize(image, 1, 99.8)` for single-channel

- [ ] **Step 3: Add run_nuclei_segmentation() convenience function**

```python
def run_nuclei_segmentation(image, modality="he", **kwargs):
    """Segment nuclei using StarDist pretrained models.

    Args:
        image: Input image. H&E: (H, W, 3) uint8 RGB. DAPI: (H, W) float/uint16.
        modality: "he" or "dapi".

    Returns:
        masks: (H, W) int32 label array, 0 = background.
        centroids_df: DataFrame with columns [y_pixel, x_pixel, nucleus_id].
    """
    segmenter = StarDistSegmenter(modality=modality)
    return segmenter.segment(image, **kwargs)
```

- [ ] **Step 4: Replace _build_cellpose_model() and run_cellpose_nuclei_segmentation()**

Delete `_build_cellpose_model()` entirely. Replace `run_cellpose_nuclei_segmentation()` with deprecation wrapper:

```python
def run_cellpose_nuclei_segmentation(image, **kwargs):
    """DEPRECATED: Use run_nuclei_segmentation() instead."""
    import warnings
    warnings.warn(
        "run_cellpose_nuclei_segmentation() is deprecated, use run_nuclei_segmentation()",
        DeprecationWarning, stacklevel=2,
    )
    masks, centroids_df = run_nuclei_segmentation(image, modality="dapi", **kwargs)
    centroids_xy = centroids_df[["x_pixel", "y_pixel"]].values
    return masks, centroids_xy
```

- [ ] **Step 5: Delete _segment_fullres_by_spot_patches()**

This function used multi-threaded Cellpose. StarDist processes the full image in one call. Remove it entirely. If any code references it, update those callers to use `StarDistSegmenter.segment()` directly.

- [ ] **Step 6: Update compute_spot_nuclei_counts_from_adata()**

This function loads an image from AnnData and runs segmentation. Update:
- Replace `run_cellpose_nuclei_segmentation(image, ...)` with `run_nuclei_segmentation(image, modality="he")`
- Update centroid handling: old code used `centroids_xy` array, new code gets `centroids_df` DataFrame
- When calling `assign_nuclei_centroids_to_spots()`, convert DataFrame to `(x, y)` array: `centroids_df[["x_pixel", "y_pixel"]].values`

- [ ] **Step 7: Update save_segmentation_artifacts()**

Change mask output filename from `*_cellpose_masks_*` to `*_stardist_masks_*`:

```python
# Old:
masks_path = out / f"{sample_name}_cellpose_masks_{result.resolution_mode}.npy"
outputs["cellpose_masks_npy"] = str(masks_path)
# New:
masks_path = out / f"{sample_name}_stardist_masks_{result.resolution_mode}.npy"
outputs["stardist_masks_npy"] = str(masks_path)
```

- [ ] **Step 8: Update _prepare_rgb_uint8() docstring**

Change "for robust Cellpose input" to "for robust segmentation input".

- [ ] **Step 9: Update tests/test_segmentation.py (if it exists)**

Check if `CITEgeist/tests/test_segmentation.py` exists. If it does:
- Replace any monkeypatching of Cellpose functions (`_build_cellpose_model`, `run_cellpose_nuclei_segmentation`) with StarDist equivalents
- Keep all tests for `assign_nuclei_centroids_to_spots`, `normalize_nuclei_counts_for_prior`, etc.

If the file does NOT exist (it may have been `tests/test_segmentation.py` at the repo root instead), search for it:
```bash
find . -name "test_segmentation*" -not -path "*/output/*" -not -path "*/logs/*"
```
Update whichever file is found.

- [ ] **Step 10: Commit**

```bash
git add CITEgeist/model/segmentation.py CITEgeist/tests/test_segmentation.py
git commit -m "feat: replace Cellpose with StarDistSegmenter in core segmentation module"
```

---

## Task 3: CitegeistModel Integration

**Files:**
- Modify: `CITEgeist/model/citegeist_model.py`

- [ ] **Step 1: Read citegeist_model.py to find all Cellpose references**

Search for `compute_spot_nuclei_counts_cellpose`, `cellpose`, `Cellpose` in the file.

- [ ] **Step 2: Rename compute_spot_nuclei_counts_cellpose() → compute_spot_nuclei_counts()**

- Update the method name
- Update the import at the top of the file (from `segmentation` import list)
- Update the method body to use `StarDistSegmenter` or `run_nuclei_segmentation()`
- Update all internal references (e.g., `self.results["nuclei_counts"]` stays the same, but log messages referencing "Cellpose" should say "StarDist" or be generic)

- [ ] **Step 3: Update all callers of the renamed method within citegeist_model.py**

Search for `compute_spot_nuclei_counts_cellpose` within the file — update all call sites. The method `run_discrete_cell_assignment()` at ~line 1261 references it.

- [ ] **Step 4: Update docstrings referencing Cellpose**

Replace "Cellpose segmentation" with "nuclei segmentation" or "StarDist segmentation" in docstrings throughout the file (~lines 543, 569, 606, 730, 1085, 1094, 1191, 1261, 1987).

- [ ] **Step 5: Commit**

```bash
git add CITEgeist/model/citegeist_model.py
git commit -m "feat: rename compute_spot_nuclei_counts_cellpose to compute_spot_nuclei_counts"
```

---

## Task 4: Unified Pipeline Updates (Steps 1-3)

**Files:**
- Modify: `CITEgeist/examples/run_unified_step2_features.py`
- Modify: `CITEgeist/examples/run_unified_step1_module3.py`
- Modify: `CITEgeist/examples/run_unified_step3_pcmil.py`
- Modify: `CITEgeist/examples/sbatch_unified_step2.sh`

- [ ] **Step 1: Update run_unified_step2_features.py — remove Cellpose, add StarDist**

Major changes:
1. Remove `_get_cellpose_model()` function entirely
2. Add import: `from model.segmentation import StarDistSegmenter`
3. In `load_image_and_segment()`, replace the Cellpose block (lines ~105-119):

```python
# Old:
cached_crop_masks = list(module3_dir.glob("*_cellpose_masks_tissue.npy"))
...
cp_model = _get_cellpose_model("nuclei", torch.cuda.is_available())
masks, _, _, _ = cp_model.eval(fullres_crop, channels=[0, 0], diameter=None)
np.save(module3_dir / f"{sample_name}_cellpose_masks_tissue.npy", masks)

# New:
cached_crop_masks = list(module3_dir.glob("*_stardist_masks_tissue.npy"))
...
segmenter = StarDistSegmenter(modality="he")
masks, centroids_df = segmenter.segment(fullres_crop)
import tensorflow as tf
tf.keras.backend.clear_session()
del segmenter  # Free TF memory before PyTorch
np.save(module3_dir / f"{sample_name}_stardist_masks_tissue.npy", masks)
```

4. Rename output dir from `cellpose` → `segmentation`:
   - `cellpose_dir = OUTPUT_BASE / sample_name / "cellpose"` → `seg_dir = OUTPUT_BASE / sample_name / "segmentation"`
   - Update all references to `cellpose_dir` → `seg_dir`

5. Add TF GPU memory growth at the top of `load_image_and_segment()` (before any StarDist import):

```python
import os
os.environ.setdefault("TF_FORCE_GPU_ALLOW_GROWTH", "true")
```

6. Update docstrings and comments.

- [ ] **Step 2: Update run_unified_step1_module3.py**

Replace `compute_spot_nuclei_counts_cellpose` → `compute_spot_nuclei_counts`:
- Line ~75: `nuclei_counts = model.compute_spot_nuclei_counts_cellpose(...)` → `nuclei_counts = model.compute_spot_nuclei_counts(...)`
- Line ~78: Update log message from "Cellpose" to "Segmentation"

- [ ] **Step 3: Update run_unified_step3_pcmil.py**

Update centroid CSV path:
- Line ~72: `centroids = pd.read_csv(base / "cellpose" / "nuclei_centroids.csv")` → `centroids = pd.read_csv(base / "segmentation" / "nuclei_centroids.csv")`

- [ ] **Step 4: Update sbatch_unified_step2.sh**

Add TF memory growth env var and update echo message:

```bash
export TF_FORCE_GPU_ALLOW_GROWTH=true

echo "Running Step 2 (StarDist + ViT features) for ${SAMPLE}"
```

Also update `CITEgeist/examples/sbatch_unified_step2_rerun.sh` (same changes).

- [ ] **Step 5: Commit**

```bash
git add CITEgeist/examples/run_unified_step2_features.py \
       CITEgeist/examples/run_unified_step1_module3.py \
       CITEgeist/examples/run_unified_step3_pcmil.py \
       CITEgeist/examples/sbatch_unified_step2.sh \
       CITEgeist/examples/sbatch_unified_step2_rerun.sh
git commit -m "feat: switch unified pipeline from Cellpose to StarDist"
```

---

## Tasks 5, 6, 7: Benchmarking Updates (PARALLELIZABLE)

> Tasks 5, 6, and 7 are independent of each other and can be executed in parallel by separate subagents. They all depend on Tasks 1-3 being complete.

## Task 5: Visium HD Benchmarking Updates

**Files:**
- Rename: `Benchmarking/visiumhd_benchmarking/src/run_cellpose_he.py` → `run_segmentation_he.py`
- Modify: `Benchmarking/visiumhd_benchmarking/src/run_benchmark.py`
- Rename: `Benchmarking/visiumhd_benchmarking/src/test_run_cellpose_he.py` → `test_run_segmentation_he.py`
- Modify: `Benchmarking/visiumhd_benchmarking/src/visualize_he_patches.py`
- Modify: `Benchmarking/visiumhd_benchmarking/src/__init__.py`
- Modify: `Benchmarking/visiumhd_benchmarking/README.md`
- Modify: `Benchmarking/visiumhd_benchmarking/slurm/sbatch_benchmark*.sh`

- [ ] **Step 1: Read run_cellpose_he.py to understand current structure**

Understand all functions: `preprocess_he_for_cellpose`, `segment_tile`, `stitch_masks`, `extract_centroids`. These need StarDist equivalents.

- [ ] **Step 2: Create run_segmentation_he.py (rename + rewrite)**

```bash
git mv Benchmarking/visiumhd_benchmarking/src/run_cellpose_he.py \
       Benchmarking/visiumhd_benchmarking/src/run_segmentation_he.py
```

Then update the file:
- Replace all Cellpose imports with `from CITEgeist.model.segmentation import StarDistSegmenter`
- Replace `preprocess_he_for_cellpose` with StarDist normalization (or remove if StarDistSegmenter handles it internally)
- Replace `segment_tile` to use `StarDistSegmenter("he").segment(tile)`
- Keep `stitch_masks` and `extract_centroids` if they are model-agnostic
- Rename function names to remove "cellpose" references

- [ ] **Step 3: Update run_benchmark.py**

Read the file first. Then:
- Replace `step2_run_cellpose()` calls with StarDist equivalent
- Replace `--cellpose-diameter` and `--skip-cellpose` CLI args with `--skip-segmentation`
- Update mask filename from `cellpose_mask.npy` → `stardist_mask.npy`
- Update imports from `run_cellpose_he` → `run_segmentation_he`

- [ ] **Step 4: Update test file**

```bash
git mv Benchmarking/visiumhd_benchmarking/src/test_run_cellpose_he.py \
       Benchmarking/visiumhd_benchmarking/src/test_run_segmentation_he.py
```

Update all imports and mock targets to reference `run_segmentation_he` instead of `run_cellpose_he`.

- [ ] **Step 5: Update visualize_he_patches.py**

Replace `cellpose_mask.npy` → `stardist_mask.npy` and "Cellpose mask" labels.

- [ ] **Step 6: Update __init__.py and README.md**

- `__init__.py`: Update module docstring referencing `run_cellpose_he`
- `README.md`: Update file references and dependency lists

- [ ] **Step 7: Update examples/run_morphology_assignment.py**

This file imports from `run_cellpose_he` (the renamed file). Update:
```python
# Old:
from run_cellpose_he import segment_wsi, extract_centroids
# New:
from run_segmentation_he import segment_wsi, extract_centroids
```

- [ ] **Step 8: Update SLURM scripts**

In `Benchmarking/visiumhd_benchmarking/slurm/sbatch_benchmark*.sh`:
- Add `export TF_FORCE_GPU_ALLOW_GROWTH=true`
- Replace `--skip-cellpose` → `--skip-segmentation` in comments/args

- [ ] **Step 9: Commit**

```bash
git add Benchmarking/visiumhd_benchmarking/
git commit -m "feat: switch Visium HD benchmarking from Cellpose to StarDist"
```

---

## Task 6: Xenium Benchmarking Updates

**Files:**
- Modify: `Benchmarking/xenium_benchmarking/CITEgeist/src/benchmark_hybrid_cellpose.py`
- Modify: `Benchmarking/xenium_benchmarking/CITEgeist/src/benchmark_multimodal_cellpose.py`
- Modify: `Benchmarking/xenium_benchmarking/CITEgeist/src/benchmark_discrete_cellpose.py`
- Modify: `Benchmarking/xenium_benchmarking/CITEgeist/src/benchmark_single_cell_resolution.py`
- Modify: `Benchmarking/xenium_benchmarking/CITEgeist/src/benchmark_cell_morphology.py`
- Modify: `Benchmarking/xenium_benchmarking/CITEgeist/src/prepare_patches.py`
- Modify: `Benchmarking/xenium_benchmarking/CITEgeist/src/unsupervised_kmeans_classifier.py`

These scripts import `run_cellpose_nuclei_segmentation` from `CITEgeist.model.segmentation`. The deprecation wrapper will make them work without changes, but we should update them to use the new API.

- [ ] **Step 1: Update benchmark_hybrid_cellpose.py**

Replace:
```python
from CITEgeist.model.segmentation import run_cellpose_nuclei_segmentation
```
with:
```python
from CITEgeist.model.segmentation import run_nuclei_segmentation
```

Update call sites: `masks, centroids_xy = run_cellpose_nuclei_segmentation(...)` → use `run_nuclei_segmentation(image, modality="dapi")` and convert centroids_df to xy array.

- [ ] **Step 2: Repeat for benchmark_multimodal_cellpose.py**

Same pattern as Step 1.

- [ ] **Step 3: Update benchmark_discrete_cellpose.py**

Also imports `SegmentationResult` — check if this type is still needed or can be replaced.

- [ ] **Step 4: Update benchmark_single_cell_resolution.py**

Same import replacement. This file uses `centroids_xy` extensively — update variable names and array access.

- [ ] **Step 5: Update benchmark_cell_morphology.py**

Replace import of `run_cellpose_nuclei_segmentation`.

- [ ] **Step 6: Update prepare_patches.py**

Add optional `--run-segmentation` flag:
```python
parser.add_argument("--run-segmentation", action="store_true",
                    help="Run StarDist segmentation instead of loading pre-computed masks")
```

If flag is set, use `StarDistSegmenter("dapi")` on the DAPI channel. Otherwise load masks as before.

- [ ] **Step 7: Update benchmark_cellpose_resolution.py**

This file calls `compute_spot_nuclei_counts_cellpose()` directly and has `--cellpose-gpu` and `--cellpose-diameter` CLI args. Update:
- Replace import and call of `compute_spot_nuclei_counts_cellpose` → use `run_nuclei_segmentation` or the renamed model method
- Update CLI args to remove cellpose-specific flags

- [ ] **Step 8: Update unsupervised_kmeans_classifier.py comment**

Replace "Cellpose pixels" with "segmentation pixels" in comment.

- [ ] **Step 9: Commit**

```bash
git add Benchmarking/xenium_benchmarking/
git commit -m "feat: update Xenium benchmarks to use StarDist via run_nuclei_segmentation"
```

---

## Task 7: Simulation Benchmarking + Misc Updates

**Files:**
- Modify: `Benchmarking/simulation_benchmarking/CITEgeist/src/generate_cellpose_images.py`
- Modify: `Benchmarking/simulation_benchmarking/CITEgeist/src/test_discretize_continuous.py`
- Modify: `Benchmarking/simulation_benchmarking/CITEgeist/src/benchmark_discrete_simulation.py`
- Modify: `examples/run_module3_unified.py`
- Modify: `examples/run_module3_discrete.py`

- [ ] **Step 1: Update generate_cellpose_images.py**

This file generates synthetic images. Update docstrings/comments that reference "Cellpose 'nuclei' model" or "Cellpose 'cyto2' model" to be generic ("nuclei segmentation model").

- [ ] **Step 2: Update test_discretize_continuous.py**

Replace `run_cellpose_nuclei_segmentation` import with `run_nuclei_segmentation`.

- [ ] **Step 3: Update benchmark_discrete_simulation.py**

Replace `run_cellpose_nuclei_segmentation` import and calls. Update `centroids_xy` variable to use DataFrame-based centroids.

- [ ] **Step 4: Update examples/run_module3_unified.py**

Replace `compute_spot_nuclei_counts_cellpose()` → `compute_spot_nuclei_counts()`.

- [ ] **Step 5: Update examples/run_module3_discrete.py**

Update docstrings referencing "Cellpose segmentation".

- [ ] **Step 6: Commit**

```bash
git add Benchmarking/simulation_benchmarking/ examples/run_module3_unified.py examples/run_module3_discrete.py
git commit -m "feat: update simulation benchmarks and examples for StarDist"
```

---

## Task 8: Docstring/Comment Cleanup + CLAUDE.md

**Files:**
- Modify: `CITEgeist/model/gurobi_impl.py`
- Modify: `CITEgeist/model/detection_estimation.py`
- Modify: `CITEgeist/model/morphology_features.py`
- Modify: `CITEgeist/model/module3b_nucleus_assignment.py`
- Modify: `CITEgeist/model/two_stage_pipeline.py`
- Modify: `CLAUDE.md`

- [ ] **Step 1: Bulk update docstrings**

In each file, search for "Cellpose" and replace with "nuclei segmentation" or "StarDist" as appropriate. These are comment-only changes — no logic changes.

Use `replace_all` where safe:
- `"Cellpose segmentation"` → `"nuclei segmentation"`
- `"from Cellpose"` → `"from nuclei segmentation"`
- `"Cellpose label mask"` → `"nuclei label mask"`

- [ ] **Step 2: Update CLAUDE.md**

In the Module 3 Alternative section, replace references to Cellpose:
- "Cellpose segmentation or Xenium cell mapping" → "StarDist segmentation or Xenium cell mapping"
- Update `segmentation.py` description to mention StarDist

- [ ] **Step 3: Commit**

```bash
git add CITEgeist/model/gurobi_impl.py \
       CITEgeist/model/detection_estimation.py \
       CITEgeist/model/morphology_features.py \
       CITEgeist/model/module3b_nucleus_assignment.py \
       CITEgeist/model/two_stage_pipeline.py \
       CLAUDE.md
git commit -m "docs: replace Cellpose references with StarDist in docstrings and CLAUDE.md"
```

---

## Task 9: Smoke Test on Real Data

**Files:** None (validation only)

- [ ] **Step 1: Clear old Step 2 outputs for one sample**

```bash
cd CITEgeist
rm -f output/unified_pipeline/HCC22-088-P1-S1/.step2_complete
rm -rf output/unified_pipeline/HCC22-088-P1-S1/segmentation/
rm -rf output/unified_pipeline/HCC22-088-P1-S1/features/
rm -f output/unified_pipeline/HCC22-088-P1-S1/module3/*stardist_masks*
```

- [ ] **Step 2: Run Step 2 on one sample via SLURM**

Create a quick single-sample test script or run interactively:

```bash
sbatch --wrap="
export TF_FORCE_GPU_ALLOW_GROWTH=true
source activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env
cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/CITEgeist
python examples/run_unified_step2_features.py --sample HCC22-088-P1-S1 --modality he
" --job-name=test_stardist --time=01:00:00 --mem=64G --gres=gpu:1 \
  --cluster=gpu --partition=a100 \
  --output=output/unified_pipeline/logs/test_stardist_%j.out \
  --error=output/unified_pipeline/logs/test_stardist_%j.err \
  --mail-type=FAIL --mail-user=alc376@pitt.edu
```

- [ ] **Step 3: Verify results**

After job completes:
1. Check nuclei count: `wc -l output/unified_pipeline/HCC22-088-P1-S1/segmentation/nuclei_centroids.csv`
2. Check features shape: `python -c "import numpy as np; print(np.load('output/unified_pipeline/HCC22-088-P1-S1/features/vit_features.npy').shape)"`
3. Compare nuclei count to old Cellpose result (was 5,415 with Cellpose nuclei model — expect significantly more with StarDist H&E model)
4. Generate a demo overlay image to visually verify nuclei look correct

- [ ] **Step 4: If smoke test passes, commit spec as final**

```bash
git add docs/superpowers/specs/2026-03-17-stardist-segmentation-switch-design.md
git commit -m "docs: finalize StarDist segmentation switch spec"
```
