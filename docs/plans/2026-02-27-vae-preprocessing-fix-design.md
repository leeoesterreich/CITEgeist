# VAE Preprocessing Fix Design

**Date**: 2026-02-27
**Status**: Approved

## Problem Statement

The VAE-based single-cell classification pipeline has poor performance (~20% accuracy) because the patch preprocessing destroys morphological signal:

1. **Per-patch z-score normalization** (lines 54-58 in `patch_extraction.py`): Each patch is individually normalized to mean=0, std=1. This destroys absolute intensity differences between cell types (e.g., bright vs dim nuclei become identical).

2. **Fixed 96x96 resize**: All patches resized to same dimensions. A large macrophage and small T cell look identical after resize.

The model can only learn from within-patch shape patterns, not from actual morphological differences between cell types.

## Solution

### 1. Global Normalization (Two Methods to Compare)

Replace per-patch z-score with global normalization. Support two methods for ablation:

| Method | Description | When to Use |
|--------|-------------|-------------|
| `minmax` | Scale by dtype max (0-65535 → 0-1) | Preserves full dynamic range |
| `percentile` | Clip to 1st/99th percentile, then scale | Robust to outliers |

```python
# In prepare_patches.py - compute stats for both methods
def compute_global_stats(
    image: np.ndarray,
    norm_method: str = "percentile",  # "minmax" or "percentile"
) -> Dict[str, Any]:
    """Compute normalization stats per channel."""
    if norm_method == "percentile":
        p1 = np.percentile(image, 1, axis=(1, 2))   # (C,)
        p99 = np.percentile(image, 99, axis=(1, 2))  # (C,)
        return {"method": "percentile", "p1": p1, "p99": p99}
    elif norm_method == "minmax":
        cmin = image.min(axis=(1, 2))  # (C,)
        cmax = image.max(axis=(1, 2))  # (C,)
        return {"method": "minmax", "min": cmin, "max": cmax}
    else:
        raise ValueError(f"Unknown norm_method: {norm_method}")

# In patch_extraction.py - apply based on method
def extract_patch(
    image: np.ndarray,
    bbox: Tuple[int, int, int, int],
    expansion: float = 0.75,
    output_size: int = 96,
    global_stats: Dict[str, Any],  # REQUIRED, not Optional
) -> np.ndarray:
    ...
    if global_stats is None:
        raise ValueError(
            "global_stats is required. Per-patch normalization is deprecated."
        )

    method = global_stats["method"]
    for c in range(C):
        if method == "percentile":
            p1, p99 = global_stats["p1"][c], global_stats["p99"][c]
            clipped = np.clip(resized[c], p1, p99)
            resized[c] = (clipped - p1) / (p99 - p1 + 1e-8)
        elif method == "minmax":
            cmin, cmax = global_stats["min"][c], global_stats["max"][c]
            resized[c] = (resized[c] - cmin) / (cmax - cmin + 1e-8)
```

### CLI Flag for Normalization Method

```python
# In prepare_patches.py
parser.add_argument(
    "--norm-method",
    type=str,
    choices=["minmax", "percentile"],
    default="percentile",
    help="Normalization method: minmax (0-65535 range) or percentile (1st/99th)"
)
```

### Ablation Experiment

Run both and compare:

```bash
# Method 1: minmax
python prepare_patches.py --norm-method minmax --output_dir output/patches_minmax ...
python train_vae.py --patches-dir output/patches_minmax --output-dir output/vae_minmax

# Method 2: percentile
python prepare_patches.py --norm-method percentile --output_dir output/patches_pct ...
python train_vae.py --patches-dir output/patches_pct --output-dir output/vae_pct

# Compare results
python evaluate_direct_softmax.py --model output/vae_minmax/...
python evaluate_direct_softmax.py --model output/vae_pct/...
```

### 2. Size Feature Encoding

Encode original bbox dimensions as features concatenated to VAE latent:

```python
def extract_patch_with_size(
    image: np.ndarray,
    bbox: Tuple[int, int, int, int],
    ...
) -> Tuple[np.ndarray, np.ndarray]:
    """Extract patch AND return size features."""
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

### 3. Model Updates

Modify `DirectSoftmaxModel` to accept size features:

```python
def __init__(self, ..., use_size_features: bool = True):
    self.use_size_features = use_size_features
    input_dim = latent_dim + 3 if use_size_features else latent_dim

    self.projection = nn.Sequential(
        nn.Linear(input_dim, 64),  # 128 + 3 = 131 if using size
        nn.ReLU(inplace=True),
        nn.Linear(64, projection_dim),
    )

def forward(self, patches, proportions, size_features=None):
    with torch.no_grad():
        mu, _ = self.encoder(patches)

    if self.use_size_features and size_features is not None:
        z = torch.cat([mu, size_features], dim=1)  # (N, 131)
    else:
        z = mu

    projected = self.projection(z)
    ...
```

## Strict Validation (No Fallbacks)

Old preprocessing was broken - no backward compatibility:

### File Requirements

| File | Required | Error if Missing |
|------|----------|------------------|
| `global_stats.npz` | Yes | "Re-run prepare_patches.py" |
| `spot_{id}_sizes.npy` | Yes | "Re-run prepare_patches.py" |
| `global_stats` param | Yes | "Per-patch normalization deprecated" |

### Checkpoint Versioning

```python
# In train_vae.py - save version marker
torch.save({
    "model_state_dict": model.state_dict(),
    "preprocessing_version": 2,  # v1 = per-patch z-score, v2 = global percentile
    ...
}, checkpoint_path)

# In train_prototypes.py - reject old checkpoints
checkpoint = torch.load(vae_checkpoint)
if "preprocessing_version" not in checkpoint:
    raise ValueError(
        f"VAE checkpoint was trained with old preprocessing. "
        f"Re-run train_vae.py with new globally-normalized patches."
    )
```

## Files to Modify

| File | Change |
|------|--------|
| `CITEgeist/model/patch_extraction.py` | Add `global_stats` param (required), add `extract_patch_with_size()` |
| `Benchmarking/.../prepare_patches.py` | Add `compute_global_stats()`, save `global_stats.npz` and `*_sizes.npy` |
| `CITEgeist/model/direct_softmax_model.py` | Add `use_size_features` param, concat size to latent |
| `CITEgeist/model/train_prototypes.py` | Update `load_spot_data()`, add CLI flags, strict validation |
| `CITEgeist/model/train_vae.py` | Add `preprocessing_version` to checkpoint |

## New Output Files

```
output/patches/
├── global_stats.npz          # NEW: {p1: (C,), p99: (C,)}
├── spot_0_patches.npy        # (N, C, 96, 96) - now globally normalized
├── spot_0_sizes.npy          # NEW: (N, 3) - [log_w, log_h, log_area]
├── spot_0_nucleus_ids.npy
├── ...
└── nucleus_features.csv
```

## Mandatory Re-run Sequence

| Step | Script | Why Required |
|------|--------|--------------|
| 1 | `prepare_patches.py` | Generate globally-normalized patches + size features |
| 2 | `train_vae.py` | VAE must learn from correct intensity/contrast patterns |
| 3 | `train_prototypes.py` | Prototypes use new VAE encoder + size features |
| 4 | `evaluate_*.py` | Measure improvement |

**Why VAE Retraining is Mandatory:**
- Old VAE input: Each patch normalized to mean=0, std=1
- New VAE input: Patches scaled by region-wide 1st/99th percentile to [0,1]
- Old VAE expects ~N(0,1) input, new data has different distribution
- Without retrain: Reconstruction loss explodes, latents meaningless

## Expected Impact

| Issue | Before | After |
|-------|--------|-------|
| Intensity differences | Destroyed (per-patch z-score) | Preserved (global norm) |
| Size differences | Destroyed (96x96 resize) | Encoded (3 features) |
| Cross-cell comparisons | Impossible | Possible |
| Single-cell accuracy | ~20% | TBD (expect significant improvement) |

## Ablation: minmax vs percentile

| Method | Pros | Cons |
|--------|------|------|
| `minmax` | Full dynamic range preserved | Sensitive to outliers (hot pixels) |
| `percentile` | Robust to outliers | May clip real signal at extremes |

**Hypothesis**: `percentile` should work better for microscopy data which often has hot pixels and background noise. But `minmax` preserves true intensity ratios if data is clean.

Run both, compare single-cell accuracy and spot-level proportion correlation.
