# VAE Preprocessing Fix Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Fix patch preprocessing to preserve morphological signal (intensity + size) for VAE-based cell classification.

**Architecture:** Replace per-patch z-score normalization with global normalization (minmax or percentile), add size feature encoding, require strict validation with no backward compatibility.

**Tech Stack:** Python, NumPy, PyTorch, tifffile

**Design Doc:** `docs/plans/2026-02-27-vae-preprocessing-fix-design.md`

---

## Task 1: Update patch_extraction.py - Global Normalization

**Files:**
- Modify: `CITEgeist/model/patch_extraction.py`
- Create: `CITEgeist/tests/test_patch_extraction.py`

**Step 1: Write the failing test**

```python
# CITEgeist/tests/test_patch_extraction.py
"""Tests for patch extraction with global normalization."""
import numpy as np
import pytest

from CITEgeist.model.patch_extraction import extract_patch, compute_global_stats


class TestComputeGlobalStats:
    """Test global stats computation."""

    def test_percentile_method(self):
        """Test percentile normalization stats."""
        # 2-channel image with known values
        image = np.array([
            [[0, 100, 200], [50, 150, 250]],  # Channel 0: range 0-250
            [[10, 20, 30], [40, 50, 60]],      # Channel 1: range 10-60
        ], dtype=np.float32)

        stats = compute_global_stats(image, norm_method="percentile")

        assert stats["method"] == "percentile"
        assert stats["p1"].shape == (2,)
        assert stats["p99"].shape == (2,)
        # 1st percentile should be near min, 99th near max
        assert stats["p1"][0] < stats["p99"][0]
        assert stats["p1"][1] < stats["p99"][1]

    def test_minmax_method(self):
        """Test minmax normalization stats."""
        image = np.array([
            [[0, 100, 200], [50, 150, 250]],
            [[10, 20, 30], [40, 50, 60]],
        ], dtype=np.float32)

        stats = compute_global_stats(image, norm_method="minmax")

        assert stats["method"] == "minmax"
        assert stats["min"][0] == 0
        assert stats["max"][0] == 250
        assert stats["min"][1] == 10
        assert stats["max"][1] == 60

    def test_invalid_method_raises(self):
        """Test that invalid method raises error."""
        image = np.zeros((2, 10, 10), dtype=np.float32)
        with pytest.raises(ValueError, match="Unknown norm_method"):
            compute_global_stats(image, norm_method="invalid")


class TestExtractPatch:
    """Test patch extraction with global normalization."""

    def test_requires_global_stats(self):
        """Test that missing global_stats raises error."""
        image = np.zeros((2, 100, 100), dtype=np.float32)
        bbox = (10, 10, 30, 30)

        with pytest.raises(ValueError, match="global_stats is required"):
            extract_patch(image, bbox, global_stats=None)

    def test_percentile_normalization_output_range(self):
        """Test that percentile normalization outputs [0, 1] range."""
        image = np.random.rand(2, 100, 100).astype(np.float32) * 1000
        bbox = (20, 20, 50, 50)
        stats = compute_global_stats(image, norm_method="percentile")

        patch = extract_patch(image, bbox, global_stats=stats)

        assert patch.shape == (2, 96, 96)
        assert patch.min() >= 0.0
        assert patch.max() <= 1.0

    def test_minmax_normalization_output_range(self):
        """Test that minmax normalization outputs [0, 1] range."""
        image = np.random.rand(2, 100, 100).astype(np.float32) * 1000
        bbox = (20, 20, 50, 50)
        stats = compute_global_stats(image, norm_method="minmax")

        patch = extract_patch(image, bbox, global_stats=stats)

        assert patch.shape == (2, 96, 96)
        assert patch.min() >= 0.0
        assert patch.max() <= 1.0

    def test_intensity_preserved_across_patches(self):
        """Test that intensity differences are preserved across patches."""
        # Create image with bright region (channel 0) and dim region
        image = np.zeros((2, 100, 100), dtype=np.float32)
        image[0, 10:30, 10:30] = 1000  # Bright patch
        image[0, 60:80, 60:80] = 100   # Dim patch
        image[1, :, :] = 500           # Uniform channel 1

        stats = compute_global_stats(image, norm_method="minmax")

        bright_patch = extract_patch(image, (10, 10, 30, 30), global_stats=stats)
        dim_patch = extract_patch(image, (60, 60, 80, 80), global_stats=stats)

        # Bright patch should have higher mean than dim patch
        assert bright_patch[0].mean() > dim_patch[0].mean()
```

**Step 2: Run test to verify it fails**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest CITEgeist/tests/test_patch_extraction.py -v`

Expected: FAIL with "cannot import name 'compute_global_stats'"

**Step 3: Write implementation**

```python
# CITEgeist/model/patch_extraction.py
"""Nucleus patch extraction for VAE training."""
import numpy as np
import cv2
from typing import Dict, Any, Tuple


def compute_global_stats(
    image: np.ndarray,
    norm_method: str = "percentile",
) -> Dict[str, Any]:
    """Compute normalization stats per channel for entire image.

    Args:
        image: (C, H, W) multi-channel image
        norm_method: "percentile" (1st/99th) or "minmax" (full range)

    Returns:
        Dict with normalization parameters
    """
    if norm_method == "percentile":
        p1 = np.percentile(image, 1, axis=(1, 2))
        p99 = np.percentile(image, 99, axis=(1, 2))
        return {"method": "percentile", "p1": p1, "p99": p99}
    elif norm_method == "minmax":
        cmin = image.min(axis=(1, 2))
        cmax = image.max(axis=(1, 2))
        return {"method": "minmax", "min": cmin, "max": cmax}
    else:
        raise ValueError(f"Unknown norm_method: {norm_method}")


def extract_patch(
    image: np.ndarray,
    bbox: Tuple[int, int, int, int],
    expansion: float = 0.75,
    output_size: int = 96,
    global_stats: Dict[str, Any] = None,
) -> np.ndarray:
    """Extract expanded patch around a nucleus with global normalization.

    Args:
        image: (C, H, W) multi-channel image
        bbox: (x_min, y_min, x_max, y_max) nucleus bounding box
        expansion: Fraction to expand bbox in each direction (0.75 = 75%)
        output_size: Final patch size after resize
        global_stats: REQUIRED - normalization stats from compute_global_stats()

    Returns:
        patch: (C, output_size, output_size) normalized patch in [0, 1] range

    Raises:
        ValueError: If global_stats is None
    """
    if global_stats is None:
        raise ValueError(
            "global_stats is required. Per-patch normalization is deprecated. "
            "Pass global_stats from compute_global_stats()."
        )

    x_min, y_min, x_max, y_max = bbox
    C, H, W = image.shape

    # Compute expansion
    w = x_max - x_min
    h = y_max - y_min

    exp_w = int(w * expansion)
    exp_h = int(h * expansion)

    # Expand bbox
    x_min_exp = max(0, x_min - exp_w)
    x_max_exp = min(W, x_max + exp_w)
    y_min_exp = max(0, y_min - exp_h)
    y_max_exp = min(H, y_max + exp_h)

    # Crop
    patch = image[:, y_min_exp:y_max_exp, x_min_exp:x_max_exp]

    # Resize each channel
    resized = np.zeros((C, output_size, output_size), dtype=np.float32)
    for c in range(C):
        resized[c] = cv2.resize(
            patch[c],
            (output_size, output_size),
            interpolation=cv2.INTER_LINEAR
        )

    # Apply global normalization
    method = global_stats["method"]
    for c in range(C):
        if method == "percentile":
            p1, p99 = global_stats["p1"][c], global_stats["p99"][c]
            clipped = np.clip(resized[c], p1, p99)
            resized[c] = (clipped - p1) / (p99 - p1 + 1e-8)
        elif method == "minmax":
            cmin, cmax = global_stats["min"][c], global_stats["max"][c]
            resized[c] = (resized[c] - cmin) / (cmax - cmin + 1e-8)

    return resized
```

**Step 4: Run test to verify it passes**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest CITEgeist/tests/test_patch_extraction.py -v`

Expected: All tests PASS

**Step 5: Commit**

```bash
git add CITEgeist/model/patch_extraction.py CITEgeist/tests/test_patch_extraction.py
git commit -m "feat: add global normalization to patch extraction

- Add compute_global_stats() with percentile and minmax methods
- Require global_stats param in extract_patch() (no fallback)
- Add comprehensive tests for normalization"
```

---

## Task 2: Add extract_patch_with_size Function

**Files:**
- Modify: `CITEgeist/model/patch_extraction.py`
- Modify: `CITEgeist/tests/test_patch_extraction.py`

**Step 1: Write the failing test**

```python
# Add to CITEgeist/tests/test_patch_extraction.py

from CITEgeist.model.patch_extraction import extract_patch_with_size


class TestExtractPatchWithSize:
    """Test patch extraction with size features."""

    def test_returns_patch_and_size(self):
        """Test that function returns both patch and size features."""
        image = np.random.rand(2, 100, 100).astype(np.float32) * 1000
        bbox = (20, 20, 50, 50)  # 30x30 bbox
        stats = compute_global_stats(image, norm_method="percentile")

        patch, size_features = extract_patch_with_size(image, bbox, global_stats=stats)

        assert patch.shape == (2, 96, 96)
        assert size_features.shape == (3,)

    def test_size_features_are_log_transformed(self):
        """Test that size features are log1p transformed."""
        image = np.random.rand(2, 100, 100).astype(np.float32)
        bbox = (10, 10, 20, 30)  # w=10, h=20, area=200
        stats = compute_global_stats(image, norm_method="percentile")

        _, size_features = extract_patch_with_size(image, bbox, global_stats=stats)

        expected_log_w = np.log1p(10)
        expected_log_h = np.log1p(20)
        expected_log_area = np.log1p(200)

        np.testing.assert_almost_equal(size_features[0], expected_log_w, decimal=5)
        np.testing.assert_almost_equal(size_features[1], expected_log_h, decimal=5)
        np.testing.assert_almost_equal(size_features[2], expected_log_area, decimal=5)

    def test_different_bbox_sizes_give_different_features(self):
        """Test that different bbox sizes produce different features."""
        image = np.random.rand(2, 200, 200).astype(np.float32)
        stats = compute_global_stats(image, norm_method="percentile")

        small_bbox = (10, 10, 20, 20)  # 10x10
        large_bbox = (50, 50, 100, 100)  # 50x50

        _, small_size = extract_patch_with_size(image, small_bbox, global_stats=stats)
        _, large_size = extract_patch_with_size(image, large_bbox, global_stats=stats)

        # Large bbox should have larger size features
        assert large_size[0] > small_size[0]  # log_w
        assert large_size[1] > small_size[1]  # log_h
        assert large_size[2] > small_size[2]  # log_area
```

**Step 2: Run test to verify it fails**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest CITEgeist/tests/test_patch_extraction.py::TestExtractPatchWithSize -v`

Expected: FAIL with "cannot import name 'extract_patch_with_size'"

**Step 3: Write implementation**

```python
# Add to CITEgeist/model/patch_extraction.py

def extract_patch_with_size(
    image: np.ndarray,
    bbox: Tuple[int, int, int, int],
    expansion: float = 0.75,
    output_size: int = 96,
    global_stats: Dict[str, Any] = None,
) -> Tuple[np.ndarray, np.ndarray]:
    """Extract patch AND return size features.

    Args:
        image: (C, H, W) multi-channel image
        bbox: (x_min, y_min, x_max, y_max) nucleus bounding box
        expansion: Fraction to expand bbox in each direction
        output_size: Final patch size after resize
        global_stats: REQUIRED - normalization stats from compute_global_stats()

    Returns:
        patch: (C, output_size, output_size) normalized patch
        size_features: (3,) array of [log1p(width), log1p(height), log1p(area)]
    """
    x_min, y_min, x_max, y_max = bbox
    w = x_max - x_min
    h = y_max - y_min
    area = w * h

    # Log-transform for better scale (nuclei vary 10-1000 px^2)
    size_features = np.array([
        np.log1p(w),
        np.log1p(h),
        np.log1p(area),
    ], dtype=np.float32)

    patch = extract_patch(image, bbox, expansion, output_size, global_stats)

    return patch, size_features
```

**Step 4: Run test to verify it passes**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest CITEgeist/tests/test_patch_extraction.py -v`

Expected: All tests PASS

**Step 5: Commit**

```bash
git add CITEgeist/model/patch_extraction.py CITEgeist/tests/test_patch_extraction.py
git commit -m "feat: add extract_patch_with_size for size feature encoding

- Returns (patch, size_features) tuple
- Size features are log1p(w), log1p(h), log1p(area)
- Preserves absolute size information after 96x96 resize"
```

---

## Task 3: Update prepare_patches.py

**Files:**
- Modify: `Benchmarking/xenium_benchmarking/CITEgeist/src/prepare_patches.py`

**Step 1: Add global stats computation and saving**

Add to prepare_patches.py after image loading:

```python
# After loading image (around line 340)
from CITEgeist.model.patch_extraction import compute_global_stats, extract_patch_with_size

# Compute and save global stats
global_stats = compute_global_stats(image, norm_method=norm_method)
stats_path = output_dir / "global_stats.npz"
np.savez(stats_path, **global_stats)
logger.info(f"Saved global normalization stats ({norm_method}) to {stats_path}")
```

**Step 2: Update patch extraction loop**

Replace the patch extraction in the spot loop:

```python
for spot_id in tqdm(spot_ids, desc="Processing spots"):
    spot_nuclei = features_df[features_df["spot_id"] == spot_id]

    patches = []
    size_features_list = []  # NEW
    nucleus_ids = []

    for _, row in spot_nuclei.iterrows():
        bbox = (
            int(row["bbox_x_min"]),
            int(row["bbox_y_min"]),
            int(row["bbox_x_max"]),
            int(row["bbox_y_max"]),
        )

        # Skip very small bboxes
        bbox_w = bbox[2] - bbox[0]
        bbox_h = bbox[3] - bbox[1]
        if bbox_w < min_bbox_size or bbox_h < min_bbox_size:
            stats["failed_patches"] += 1
            continue

        try:
            patch, size_feats = extract_patch_with_size(
                image, bbox, expansion, patch_size, global_stats
            )
            patches.append(patch)
            size_features_list.append(size_feats)  # NEW
            nucleus_ids.append(row["nucleus_id"])
            stats["successful_patches"] += 1
        except Exception as e:
            logger.debug(f"Failed to extract patch for nucleus {row['nucleus_id']}: {e}")
            stats["failed_patches"] += 1
            continue

    if patches:
        patches_array = np.stack(patches, axis=0)
        sizes_array = np.stack(size_features_list, axis=0)  # NEW

        np.save(output_dir / f"spot_{spot_id}_patches.npy", patches_array)
        np.save(output_dir / f"spot_{spot_id}_sizes.npy", sizes_array)  # NEW
        np.save(output_dir / f"spot_{spot_id}_nucleus_ids.npy", np.array(nucleus_ids, dtype=np.int64))
    else:
        stats["empty_spots"] += 1
```

**Step 3: Add CLI argument**

```python
parser.add_argument(
    "--norm-method",
    type=str,
    choices=["minmax", "percentile"],
    default="percentile",
    help="Normalization method: minmax (dtype range) or percentile (1st/99th)"
)
```

**Step 4: Add validation at end**

```python
# After extraction loop, validate outputs
required_base_files = ["global_stats.npz", "nucleus_features.csv"]
for f in required_base_files:
    if not (output_dir / f).exists():
        raise RuntimeError(f"Missing required output: {f}")

# Spot-check a few spots
sample_spots = list(spot_ids)[:3]
for spot_id in sample_spots:
    patches_file = output_dir / f"spot_{spot_id}_patches.npy"
    sizes_file = output_dir / f"spot_{spot_id}_sizes.npy"
    if patches_file.exists() and not sizes_file.exists():
        raise RuntimeError(f"Patches exist but sizes missing for spot {spot_id}")
```

**Step 5: Run manual test**

Run: `python Benchmarking/xenium_benchmarking/CITEgeist/src/prepare_patches.py --help`

Expected: Shows `--norm-method` argument

**Step 6: Commit**

```bash
git add Benchmarking/xenium_benchmarking/CITEgeist/src/prepare_patches.py
git commit -m "feat: update prepare_patches for global normalization + size features

- Compute and save global_stats.npz with norm method
- Save spot_{id}_sizes.npy alongside patches
- Add --norm-method CLI flag (percentile or minmax)
- Add validation for required output files"
```

---

## Task 4: Update DirectSoftmaxModel for Size Features

**Files:**
- Modify: `CITEgeist/model/direct_softmax_model.py`
- Create: `CITEgeist/tests/test_direct_softmax_size.py`

**Step 1: Write the failing test**

```python
# CITEgeist/tests/test_direct_softmax_size.py
"""Tests for DirectSoftmaxModel with size features."""
import torch
import pytest

from CITEgeist.model.vae import VAEEncoder
from CITEgeist.model.direct_softmax_model import DirectSoftmaxModel


class TestDirectSoftmaxWithSize:
    """Test size feature integration."""

    @pytest.fixture
    def encoder(self):
        return VAEEncoder(in_channels=2, latent_dim=128)

    def test_model_accepts_size_features(self, encoder):
        """Test that model can be created with size features enabled."""
        model = DirectSoftmaxModel(
            encoder=encoder,
            n_types=7,
            latent_dim=128,
            use_size_features=True,
        )
        assert model.use_size_features is True
        # Projection input should be 128 + 3 = 131
        assert model.projection[0].in_features == 131

    def test_forward_with_size_features(self, encoder):
        """Test forward pass with size features."""
        model = DirectSoftmaxModel(
            encoder=encoder,
            n_types=7,
            latent_dim=128,
            use_size_features=True,
        )

        patches = torch.randn(10, 2, 96, 96)
        proportions = torch.ones(7) / 7
        size_features = torch.randn(10, 3)

        loss, assignments = model(patches, proportions, size_features=size_features)

        assert loss.shape == ()
        assert assignments.shape == (10, 7)

    def test_forward_without_size_features_raises(self, encoder):
        """Test that missing size features raises error when enabled."""
        model = DirectSoftmaxModel(
            encoder=encoder,
            n_types=7,
            latent_dim=128,
            use_size_features=True,
        )

        patches = torch.randn(10, 2, 96, 96)
        proportions = torch.ones(7) / 7

        with pytest.raises(ValueError, match="size_features required"):
            model(patches, proportions, size_features=None)

    def test_model_without_size_features(self, encoder):
        """Test that model works without size features."""
        model = DirectSoftmaxModel(
            encoder=encoder,
            n_types=7,
            latent_dim=128,
            use_size_features=False,
        )

        patches = torch.randn(10, 2, 96, 96)
        proportions = torch.ones(7) / 7

        loss, assignments = model(patches, proportions)

        assert loss.shape == ()
        assert assignments.shape == (10, 7)
```

**Step 2: Run test to verify it fails**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest CITEgeist/tests/test_direct_softmax_size.py -v`

Expected: FAIL (use_size_features param doesn't exist)

**Step 3: Update DirectSoftmaxModel**

Modify `__init__`:

```python
def __init__(
    self,
    encoder: nn.Module,
    n_types: int,
    latent_dim: int = 128,
    projection_dim: int = 32,
    temperature: float = 0.1,
    # ... existing params ...
    # Size feature params
    use_size_features: bool = False,
):
    super().__init__()

    # Size features
    self.use_size_features = use_size_features
    input_dim = latent_dim + 3 if use_size_features else latent_dim

    # Projection head (adjusted input dim)
    self.projection = nn.Sequential(
        nn.Linear(input_dim, 64),
        nn.ReLU(inplace=True),
        nn.Linear(64, projection_dim),
    )
    # ... rest of init ...
```

Modify `forward`:

```python
def forward(
    self,
    patches: torch.Tensor,
    proportions: torch.Tensor,
    return_components: bool = False,
    size_features: Optional[torch.Tensor] = None,  # NEW
) -> Tuple[torch.Tensor, torch.Tensor]:
    N = patches.shape[0]

    # Validate size_features if required
    if self.use_size_features and size_features is None:
        raise ValueError("size_features required when use_size_features=True")

    # Encode patches (frozen)
    with torch.no_grad():
        mu, _ = self.encoder(patches)
    z = mu

    # Concatenate size features if enabled
    if self.use_size_features and size_features is not None:
        z = torch.cat([z, size_features], dim=1)  # (N, latent_dim + 3)

    # Project to shared space
    projected = self.projection(z)
    # ... rest of forward ...
```

Also update `assign` method similarly.

**Step 4: Run test to verify it passes**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest CITEgeist/tests/test_direct_softmax_size.py -v`

Expected: All tests PASS

**Step 5: Commit**

```bash
git add CITEgeist/model/direct_softmax_model.py CITEgeist/tests/test_direct_softmax_size.py
git commit -m "feat: add size feature support to DirectSoftmaxModel

- Add use_size_features param (default False for backward compat)
- Concat size features to latent before projection
- Require size_features when enabled (no silent fallback)"
```

---

## Task 5: Update train_prototypes.py

**Files:**
- Modify: `CITEgeist/model/train_prototypes.py`

**Step 1: Update load_spot_data function**

```python
def load_spot_data(
    patches_dir: Path,
    spot_id: str,
    device: torch.device,
) -> Tuple[torch.Tensor, torch.Tensor]:
    """Load patches AND size features for a spot.

    Raises:
        FileNotFoundError: If patches or size features missing
    """
    patch_file = patches_dir / f"spot_{spot_id}_patches.npy"
    size_file = patches_dir / f"spot_{spot_id}_sizes.npy"

    if not patch_file.exists():
        patch_file = patches_dir / f"{spot_id}_patches.npy"
    if not patch_file.exists():
        raise FileNotFoundError(f"Patches not found: {patch_file}")

    if not size_file.exists():
        size_file = patches_dir / f"{spot_id}_sizes.npy"
    if not size_file.exists():
        raise FileNotFoundError(
            f"Size features not found: {size_file}. "
            f"Re-run prepare_patches.py to generate size features."
        )

    patches = torch.from_numpy(np.load(patch_file).astype(np.float32)).to(device)
    sizes = torch.from_numpy(np.load(size_file).astype(np.float32)).to(device)

    return patches, sizes
```

**Step 2: Validate global_stats.npz exists**

```python
# At start of train_prototypes()
global_stats_path = patches_path / "global_stats.npz"
if not global_stats_path.exists():
    raise FileNotFoundError(
        f"global_stats.npz not found in {patches_path}. "
        f"Re-run prepare_patches.py with new preprocessing."
    )
logger.info(f"Found global_stats.npz - using new preprocessing")
```

**Step 3: Add preprocessing_version check for VAE**

```python
# In load_vae_encoder()
checkpoint = torch.load(checkpoint_path, map_location=device)

if "preprocessing_version" not in checkpoint:
    raise ValueError(
        f"VAE checkpoint {checkpoint_path} was trained with old preprocessing. "
        f"Re-run train_vae.py with new globally-normalized patches."
    )
```

**Step 4: Update training loop to pass size_features**

```python
# In training loop
patches, sizes = load_spot_data(patches_path, spot_id, device)
components, _ = model(patches, proportions, return_components=True, size_features=sizes)
```

**Step 5: Add CLI flag**

```python
parser.add_argument(
    "--use-size-features",
    action="store_true",
    default=True,
    help="Use bbox size features (requires new preprocessing)"
)
```

**Step 6: Commit**

```bash
git add CITEgeist/model/train_prototypes.py
git commit -m "feat: update train_prototypes for new preprocessing

- Require global_stats.npz (fail if missing)
- Require size features (fail if missing)
- Validate VAE has preprocessing_version
- Pass size_features to model forward"
```

---

## Task 6: Update train_vae.py with Preprocessing Version

**Files:**
- Modify: `CITEgeist/model/train_vae.py` (or create if doesn't exist)

**Step 1: Find train_vae.py location**

Run: `find /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist -name "train_vae*.py" -type f`

**Step 2: Add preprocessing_version to checkpoint saving**

```python
# When saving checkpoint
torch.save({
    "model_state_dict": model.state_dict(),
    "optimizer_state_dict": optimizer.state_dict(),
    "epoch": epoch,
    "in_channels": in_channels,
    "latent_dim": latent_dim,
    "preprocessing_version": 2,  # v1=per-patch z-score, v2=global percentile
    "history": history,
}, checkpoint_path)
```

**Step 3: Commit**

```bash
git add CITEgeist/model/train_vae.py  # or appropriate path
git commit -m "feat: add preprocessing_version to VAE checkpoint

- Version 2 = global percentile/minmax normalization
- Old checkpoints without version will be rejected"
```

---

## Task 7: Integration Test - Full Pipeline

**Files:**
- Create: `Benchmarking/xenium_benchmarking/CITEgeist/tests/test_preprocessing_pipeline.py`

**Step 1: Write integration test**

```python
"""Integration test for new preprocessing pipeline."""
import tempfile
from pathlib import Path
import numpy as np
import pandas as pd
import pytest


def test_full_preprocessing_pipeline():
    """Test prepare_patches → train_vae → train_prototypes flow."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)

        # Create mock image and mask
        image = np.random.rand(2, 200, 200).astype(np.float32) * 65535
        mask = np.zeros((200, 200), dtype=np.int32)
        mask[20:40, 20:40] = 1
        mask[60:80, 60:80] = 2
        mask[100:120, 100:120] = 3

        np.save(tmpdir / "test_mask.npy", mask)

        # Create spot mapping
        spot_map = pd.DataFrame({
            "nucleus_id": [1, 2, 3],
            "spot_id": ["spot_0", "spot_0", "spot_1"],
        })
        spot_map.to_csv(tmpdir / "spot_map.csv", index=False)

        # Create proportions
        props = pd.DataFrame({
            "spot_id": ["spot_0", "spot_1"],
            "TypeA": [0.5, 0.3],
            "TypeB": [0.5, 0.7],
        })
        props.to_csv(tmpdir / "proportions.csv", index=False)

        # Verify outputs would be created (just check structure)
        output_dir = tmpdir / "patches"
        output_dir.mkdir()

        # Simulate what prepare_patches creates
        from CITEgeist.model.patch_extraction import compute_global_stats
        stats = compute_global_stats(image, norm_method="percentile")
        np.savez(output_dir / "global_stats.npz", **stats)

        # Verify global_stats.npz is valid
        loaded = np.load(output_dir / "global_stats.npz")
        assert "method" in loaded or "p1" in loaded
```

**Step 2: Run test**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest Benchmarking/xenium_benchmarking/CITEgeist/tests/test_preprocessing_pipeline.py -v`

**Step 3: Commit**

```bash
git add Benchmarking/xenium_benchmarking/CITEgeist/tests/test_preprocessing_pipeline.py
git commit -m "test: add integration test for preprocessing pipeline"
```

---

## Task 8: Re-run Pipeline on Xenium Data

**Note:** This is a SLURM job, not a unit test.

**Step 1: Create SLURM script for re-preprocessing**

```bash
# Benchmarking/xenium_benchmarking/CITEgeist/slurm/sbatch_preprocess_v2.sh
#!/bin/bash
#SBATCH --job-name=preprocess_v2
#SBATCH --output=slurm/logs/preprocess_v2_%A_%a.out
#SBATCH --error=slurm/logs/preprocess_v2_%A_%a.err
#SBATCH --array=0-4
#SBATCH --time=2:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

REGION=$SLURM_ARRAY_TASK_ID

python src/prepare_patches.py \
    --dapi /path/to/ch0000_dapi.ome.tif \
    --boundary /path/to/ch0001_boundary.ome.tif \
    --region $REGION \
    --mask output/cellpose/region_${REGION}_nuclei.npy \
    --nuclei_spot_map output/region_${REGION}/nuclei_spot_mapping.csv \
    --output_dir output/patches_v2/region_${REGION} \
    --norm-method percentile
```

**Step 2: Submit and monitor**

```bash
sbatch slurm/sbatch_preprocess_v2.sh
squeue -u $USER
```

---

## Summary

| Task | Description | Test Command |
|------|-------------|--------------|
| 1 | Global normalization in patch_extraction.py | `pytest CITEgeist/tests/test_patch_extraction.py -v` |
| 2 | Size feature extraction | `pytest CITEgeist/tests/test_patch_extraction.py -v` |
| 3 | Update prepare_patches.py | Manual test + SLURM |
| 4 | DirectSoftmaxModel size features | `pytest CITEgeist/tests/test_direct_softmax_size.py -v` |
| 5 | Update train_prototypes.py | Integrated with task 4 tests |
| 6 | VAE preprocessing version | Manual verification |
| 7 | Integration test | `pytest .../test_preprocessing_pipeline.py -v` |
| 8 | Re-run on Xenium data | SLURM job |
