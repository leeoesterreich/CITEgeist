# VICReg + Attention MIL Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Improve single-cell classification accuracy from 20.2% by adding discriminative VAE pretraining (VICReg) and attention-weighted MIL aggregation.

**Architecture:** Three-layer enhancement: (1) VICReg loss during VAE pretraining creates discriminative latents, (2) attention-weighted aggregation downweights ambiguous nuclei, (3) per-class attention heads let rare cell types focus on their characteristic morphology.

**Tech Stack:** PyTorch, torchvision (augmentations), existing VAE/DirectSoftmax infrastructure

---

## Task 1: Geometric Augmentation Module

**Files:**
- Create: `CITEgeist/model/augmentations.py`
- Test: `CITEgeist/tests/test_augmentations.py`

**Step 1: Write the failing test**

```python
# CITEgeist/tests/test_augmentations.py
"""Tests for geometric augmentation module."""
import torch
import pytest


def test_geometric_augment_preserves_shape():
    """Augmented patch should have same shape as input."""
    from CITEgeist.model.augmentations import GeometricAugment

    aug = GeometricAugment()
    x = torch.randn(2, 96, 96)  # 2-channel 96x96 patch
    x_aug = aug(x)

    assert x_aug.shape == x.shape


def test_geometric_augment_creates_different_views():
    """Two augmentations of same input should differ."""
    from CITEgeist.model.augmentations import GeometricAugment

    aug = GeometricAugment()
    x = torch.randn(2, 96, 96)

    torch.manual_seed(0)
    x_aug1 = aug(x)
    torch.manual_seed(1)
    x_aug2 = aug(x)

    # Views should be different (not identical)
    assert not torch.allclose(x_aug1, x_aug2)


def test_geometric_augment_batch():
    """Should work on batched input."""
    from CITEgeist.model.augmentations import GeometricAugment

    aug = GeometricAugment()
    x = torch.randn(8, 2, 96, 96)  # Batch of 8
    x_aug = aug(x)

    assert x_aug.shape == x.shape
```

**Step 2: Run test to verify it fails**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest CITEgeist/tests/test_augmentations.py -v`

Expected: FAIL with "No module named 'CITEgeist.model.augmentations'"

**Step 3: Write minimal implementation**

```python
# CITEgeist/model/augmentations.py
"""Geometric augmentations for nucleus patches.

Provides augmentations safe for microscopy data that preserve
intensity relationships while creating diverse views for VICReg.
"""
import torch
import torch.nn as nn
import torchvision.transforms.functional as TF
import random


class GeometricAugment(nn.Module):
    """Geometric-only augmentations for nucleus patches.

    Applies random combination of:
    - Rotation: 0, 90, 180, or 270 degrees
    - Horizontal flip (50% probability)
    - Vertical flip (50% probability)
    - Small translation: +/- 5 pixels

    These augmentations are safe for microscopy data as they
    don't alter intensity relationships between channels.
    """

    def __init__(self, max_translate: int = 5):
        """Initialize augmentation.

        Args:
            max_translate: Maximum translation in pixels.
        """
        super().__init__()
        self.max_translate = max_translate

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """Apply random geometric augmentations.

        Args:
            x: Input tensor of shape (C, H, W) or (B, C, H, W)

        Returns:
            Augmented tensor with same shape
        """
        # Handle both single image and batch
        if x.dim() == 3:
            return self._augment_single(x)
        elif x.dim() == 4:
            return torch.stack([self._augment_single(img) for img in x])
        else:
            raise ValueError(f"Expected 3D or 4D tensor, got {x.dim()}D")

    def _augment_single(self, x: torch.Tensor) -> torch.Tensor:
        """Augment a single image (C, H, W)."""
        # Random rotation (0, 90, 180, 270)
        k = random.randint(0, 3)
        x = torch.rot90(x, k, dims=(1, 2))

        # Random horizontal flip
        if random.random() > 0.5:
            x = torch.flip(x, dims=(2,))

        # Random vertical flip
        if random.random() > 0.5:
            x = torch.flip(x, dims=(1,))

        # Small random translation via padding and cropping
        if self.max_translate > 0:
            dx = random.randint(-self.max_translate, self.max_translate)
            dy = random.randint(-self.max_translate, self.max_translate)
            x = self._translate(x, dx, dy)

        return x

    def _translate(self, x: torch.Tensor, dx: int, dy: int) -> torch.Tensor:
        """Translate image by (dx, dy) pixels with zero padding."""
        C, H, W = x.shape
        result = torch.zeros_like(x)

        # Source and destination slices
        src_y = slice(max(0, -dy), min(H, H - dy))
        src_x = slice(max(0, -dx), min(W, W - dx))
        dst_y = slice(max(0, dy), min(H, H + dy))
        dst_x = slice(max(0, dx), min(W, W + dx))

        result[:, dst_y, dst_x] = x[:, src_y, src_x]
        return result
```

**Step 4: Run test to verify it passes**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest CITEgeist/tests/test_augmentations.py -v`

Expected: PASS (3 tests)

**Step 5: Commit**

```bash
git add CITEgeist/model/augmentations.py CITEgeist/tests/test_augmentations.py
git commit -m "feat: add geometric augmentation module for VICReg"
```

---

## Task 2: VICReg Loss Module

**Files:**
- Create: `CITEgeist/model/vicreg.py`
- Test: `CITEgeist/tests/test_vicreg.py`

**Step 1: Write the failing test**

```python
# CITEgeist/tests/test_vicreg.py
"""Tests for VICReg loss module."""
import torch
import pytest


def test_vicreg_loss_returns_components():
    """VICReg should return total loss and components."""
    from CITEgeist.model.vicreg import vicreg_loss

    z = torch.randn(32, 128)  # Batch of 32, 128-dim latents
    z_aug = torch.randn(32, 128)  # Augmented view

    loss, components = vicreg_loss(z, z_aug)

    assert isinstance(loss, torch.Tensor)
    assert loss.dim() == 0  # Scalar
    assert "invariance" in components
    assert "variance" in components
    assert "covariance" in components


def test_vicreg_invariance_zero_for_identical():
    """Invariance loss should be 0 for identical views."""
    from CITEgeist.model.vicreg import vicreg_loss

    z = torch.randn(32, 128)

    _, components = vicreg_loss(z, z)  # Same input twice

    assert components["invariance"].item() < 1e-6


def test_vicreg_variance_penalizes_collapse():
    """Variance loss should be high when all embeddings identical."""
    from CITEgeist.model.vicreg import vicreg_loss

    # Collapsed embeddings: all same vector
    z = torch.ones(32, 128)
    z_aug = torch.ones(32, 128) + torch.randn(32, 128) * 0.01

    _, components = vicreg_loss(z, z_aug)

    # Variance loss should be high (std < 1 in all dims)
    assert components["variance"].item() > 100  # 128 dims * ~1 each


def test_vicreg_covariance_penalizes_correlation():
    """Covariance loss should penalize correlated dimensions."""
    from CITEgeist.model.vicreg import vicreg_loss

    # Create embeddings with correlated dimensions
    base = torch.randn(32, 1)
    z = base.expand(32, 128).clone()  # All dims identical = highly correlated
    z_aug = z + torch.randn(32, 128) * 0.1

    _, components = vicreg_loss(z, z_aug)

    # Covariance loss should be high
    assert components["covariance"].item() > 10
```

**Step 2: Run test to verify it fails**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest CITEgeist/tests/test_vicreg.py -v`

Expected: FAIL with "No module named 'CITEgeist.model.vicreg'"

**Step 3: Write minimal implementation**

```python
# CITEgeist/model/vicreg.py
"""VICReg loss for self-supervised representation learning.

Implements Variance-Invariance-Covariance Regularization from
Bardes et al. 2021 (https://arxiv.org/abs/2105.04906).

The loss encourages:
- Invariance: Embeddings of augmented views should match
- Variance: Each dimension should have std >= 1 (prevents collapse)
- Covariance: Dimensions should be decorrelated (prevents redundancy)
"""
import torch
import torch.nn.functional as F
from typing import Tuple, Dict


def vicreg_loss(
    z: torch.Tensor,
    z_aug: torch.Tensor,
    invariance_weight: float = 25.0,
    variance_weight: float = 25.0,
    covariance_weight: float = 1.0,
    variance_target: float = 1.0,
) -> Tuple[torch.Tensor, Dict[str, torch.Tensor]]:
    """Compute VICReg loss between two views.

    Args:
        z: Embeddings from view 1, shape (B, D)
        z_aug: Embeddings from view 2 (augmented), shape (B, D)
        invariance_weight: Weight for invariance term (default 25.0)
        variance_weight: Weight for variance term (default 25.0)
        covariance_weight: Weight for covariance term (default 1.0)
        variance_target: Target std for each dimension (default 1.0)

    Returns:
        total_loss: Weighted sum of all components
        components: Dict with individual loss terms
    """
    batch_size, dim = z.shape

    # === Invariance Loss ===
    # MSE between embeddings of augmented views
    invariance_loss = F.mse_loss(z, z_aug)

    # === Variance Loss ===
    # Encourage std >= variance_target in each dimension
    # Hinge loss: penalize if std < target
    std_z = z.std(dim=0)
    std_z_aug = z_aug.std(dim=0)

    variance_loss = (
        F.relu(variance_target - std_z).mean() +
        F.relu(variance_target - std_z_aug).mean()
    )

    # === Covariance Loss ===
    # Penalize off-diagonal covariance (encourages decorrelated dims)
    z_centered = z - z.mean(dim=0)
    z_aug_centered = z_aug - z_aug.mean(dim=0)

    cov_z = (z_centered.T @ z_centered) / (batch_size - 1)
    cov_z_aug = (z_aug_centered.T @ z_aug_centered) / (batch_size - 1)

    # Zero out diagonal (we only penalize off-diagonal)
    off_diag_z = cov_z - torch.diag(cov_z.diag())
    off_diag_z_aug = cov_z_aug - torch.diag(cov_z_aug.diag())

    covariance_loss = (
        (off_diag_z ** 2).sum() / dim +
        (off_diag_z_aug ** 2).sum() / dim
    )

    # === Total Loss ===
    total_loss = (
        invariance_weight * invariance_loss +
        variance_weight * variance_loss +
        covariance_weight * covariance_loss
    )

    components = {
        "invariance": invariance_loss,
        "variance": variance_loss,
        "covariance": covariance_loss,
    }

    return total_loss, components
```

**Step 4: Run test to verify it passes**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest CITEgeist/tests/test_vicreg.py -v`

Expected: PASS (4 tests)

**Step 5: Commit**

```bash
git add CITEgeist/model/vicreg.py CITEgeist/tests/test_vicreg.py
git commit -m "feat: add VICReg loss module"
```

---

## Task 3: Integrate VICReg into VAE Training

**Files:**
- Modify: `CITEgeist/model/train_vae.py`
- Test: `CITEgeist/tests/test_train_vae_vicreg.py`

**Step 1: Write the failing test**

```python
# CITEgeist/tests/test_train_vae_vicreg.py
"""Tests for VICReg-enhanced VAE training."""
import torch
import pytest
import tempfile
import numpy as np
from pathlib import Path


def test_train_vae_with_vicreg_flag():
    """Training with --use-vicreg should include VICReg loss."""
    from CITEgeist.model.train_vae import train_vae

    # Create temporary patches
    with tempfile.TemporaryDirectory() as tmpdir:
        patches_dir = Path(tmpdir) / "patches"
        patches_dir.mkdir()
        output_dir = Path(tmpdir) / "output"

        # Create dummy patches
        patches = np.random.randn(100, 2, 96, 96).astype(np.float32)
        np.save(patches_dir / "test_patches.npy", patches)

        # Train with VICReg
        history = train_vae(
            patches_dir=str(patches_dir),
            output_dir=str(output_dir),
            in_channels=2,
            latent_dim=32,
            epochs=2,
            batch_size=16,
            device="cpu",
            use_vicreg=True,
        )

        # Should have VICReg loss components in history
        assert "vicreg_loss" in history
        assert "vicreg_invariance" in history
        assert "vicreg_variance" in history
        assert "vicreg_covariance" in history


def test_train_vae_without_vicreg_unchanged():
    """Training without --use-vicreg should work as before."""
    from CITEgeist.model.train_vae import train_vae

    with tempfile.TemporaryDirectory() as tmpdir:
        patches_dir = Path(tmpdir) / "patches"
        patches_dir.mkdir()
        output_dir = Path(tmpdir) / "output"

        patches = np.random.randn(100, 2, 96, 96).astype(np.float32)
        np.save(patches_dir / "test_patches.npy", patches)

        history = train_vae(
            patches_dir=str(patches_dir),
            output_dir=str(output_dir),
            in_channels=2,
            latent_dim=32,
            epochs=2,
            batch_size=16,
            device="cpu",
            use_vicreg=False,
        )

        # Should NOT have VICReg components
        assert "vicreg_loss" not in history
        assert "loss" in history
        assert "recon_loss" in history
```

**Step 2: Run test to verify it fails**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest CITEgeist/tests/test_train_vae_vicreg.py -v`

Expected: FAIL with "TypeError: train_vae() got an unexpected keyword argument 'use_vicreg'"

**Step 3: Modify train_vae.py to add VICReg support**

Add to imports (after line 28):
```python
try:
    from .augmentations import GeometricAugment
    from .vicreg import vicreg_loss
except ImportError:
    from augmentations import GeometricAugment
    from vicreg import vicreg_loss
```

Modify `train_vae` function signature (line 91-102) to add new parameters:
```python
def train_vae(
    patches_dir: str,
    output_dir: str,
    in_channels: int = 3,
    latent_dim: int = 128,
    batch_size: int = 64,
    epochs: int = 100,
    lr: float = 1e-4,
    beta: float = 0.5,
    device: str = "cuda",
    checkpoint_interval: int = 10,
    num_workers: int = 4,
    # VICReg parameters
    use_vicreg: bool = False,
    vicreg_weight: float = 1.0,
    vicreg_invariance: float = 25.0,
    vicreg_variance: float = 25.0,
    vicreg_covariance: float = 1.0,
) -> Dict[str, List[float]]:
```

Update history initialization (after line 151):
```python
    history = {
        "loss": [],
        "recon_loss": [],
        "kl_loss": [],
    }
    if use_vicreg:
        history.update({
            "vicreg_loss": [],
            "vicreg_invariance": [],
            "vicreg_variance": [],
            "vicreg_covariance": [],
        })
        augment = GeometricAugment()
        logger.info("VICReg enabled with geometric augmentations")
```

Modify training loop (replace lines 167-194):
```python
        pbar = tqdm(dataloader, desc=f"Epoch {epoch+1}/{epochs}")
        for batch in pbar:
            batch = batch.to(device)

            # Forward pass
            x_recon, mu, logvar = model(batch)

            # Compute VAE loss
            vae_loss, recon_loss, kl_loss = VAE.loss_function(
                batch, x_recon, mu, logvar, beta=beta
            )

            loss = vae_loss

            # Add VICReg loss if enabled
            if use_vicreg:
                # Create augmented view
                batch_aug = augment(batch)
                _, mu_aug, _ = model(batch_aug)

                # Compute VICReg on latents
                vic_loss, vic_components = vicreg_loss(
                    mu, mu_aug,
                    invariance_weight=vicreg_invariance,
                    variance_weight=vicreg_variance,
                    covariance_weight=vicreg_covariance,
                )
                loss = loss + vicreg_weight * vic_loss

                epoch_vicreg += vic_loss.item()
                epoch_vic_inv += vic_components["invariance"].item()
                epoch_vic_var += vic_components["variance"].item()
                epoch_vic_cov += vic_components["covariance"].item()

            # Backward pass
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()

            # Accumulate metrics
            epoch_loss += loss.item()
            epoch_recon += recon_loss.item()
            epoch_kl += kl_loss.item()
            n_batches += 1

            # Update progress bar
            postfix = {
                "loss": f"{loss.item():.4f}",
                "recon": f"{recon_loss.item():.4f}",
                "kl": f"{kl_loss.item():.4f}",
            }
            if use_vicreg:
                postfix["vic"] = f"{vic_loss.item():.4f}"
            pbar.set_postfix(postfix)
```

Initialize VICReg accumulators before the epoch loop (after line 159):
```python
        epoch_vicreg = 0.0
        epoch_vic_inv = 0.0
        epoch_vic_var = 0.0
        epoch_vic_cov = 0.0
```

Update history recording (after averaging, around line 200):
```python
        if use_vicreg:
            history["vicreg_loss"].append(epoch_vicreg / n_batches)
            history["vicreg_invariance"].append(epoch_vic_inv / n_batches)
            history["vicreg_variance"].append(epoch_vic_var / n_batches)
            history["vicreg_covariance"].append(epoch_vic_cov / n_batches)
```

Add CLI arguments (after line 293):
```python
    parser.add_argument(
        "--use-vicreg",
        action="store_true",
        help="Enable VICReg loss for discriminative pretraining"
    )
    parser.add_argument(
        "--vicreg-weight",
        type=float,
        default=1.0,
        help="Weight for VICReg loss relative to VAE loss"
    )
    parser.add_argument(
        "--vicreg-invariance",
        type=float,
        default=25.0,
        help="VICReg invariance term weight"
    )
    parser.add_argument(
        "--vicreg-variance",
        type=float,
        default=25.0,
        help="VICReg variance term weight"
    )
    parser.add_argument(
        "--vicreg-covariance",
        type=float,
        default=1.0,
        help="VICReg covariance term weight"
    )
```

Pass new args to train_vae() call (around line 316):
```python
        use_vicreg=args.use_vicreg,
        vicreg_weight=args.vicreg_weight,
        vicreg_invariance=args.vicreg_invariance,
        vicreg_variance=args.vicreg_variance,
        vicreg_covariance=args.vicreg_covariance,
```

**Step 4: Run test to verify it passes**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest CITEgeist/tests/test_train_vae_vicreg.py -v`

Expected: PASS (2 tests)

**Step 5: Commit**

```bash
git add CITEgeist/model/train_vae.py CITEgeist/tests/test_train_vae_vicreg.py
git commit -m "feat: integrate VICReg loss into VAE training"
```

---

## Task 4: Attention Aggregator Module

**Files:**
- Create: `CITEgeist/model/attention_aggregator.py`
- Test: `CITEgeist/tests/test_attention_aggregator.py`

**Step 1: Write the failing test**

```python
# CITEgeist/tests/test_attention_aggregator.py
"""Tests for attention-weighted MIL aggregation."""
import torch
import pytest


def test_attention_aggregator_output_shape():
    """Aggregated proportions should sum to 1."""
    from CITEgeist.model.attention_aggregator import AttentionAggregator

    agg = AttentionAggregator(embed_dim=32, n_types=7)

    embeddings = torch.randn(20, 32)  # 20 nuclei, 32-dim
    soft_assignments = torch.softmax(torch.randn(20, 7), dim=1)

    pred_props, entropy, attn_weights = agg(embeddings, soft_assignments)

    assert pred_props.shape == (7,)
    assert torch.allclose(pred_props.sum(), torch.tensor(1.0), atol=1e-5)
    assert attn_weights.shape == (20, 1)


def test_attention_weights_sum_to_one():
    """Attention weights should sum to 1 across nuclei."""
    from CITEgeist.model.attention_aggregator import AttentionAggregator

    agg = AttentionAggregator(embed_dim=32, n_types=7)

    embeddings = torch.randn(20, 32)
    soft_assignments = torch.softmax(torch.randn(20, 7), dim=1)

    _, _, attn_weights = agg(embeddings, soft_assignments)

    assert torch.allclose(attn_weights.sum(), torch.tensor(1.0), atol=1e-5)


def test_per_class_attention_output_shape():
    """Per-class attention should produce K separate attention heads."""
    from CITEgeist.model.attention_aggregator import PerClassAttentionAggregator

    agg = PerClassAttentionAggregator(embed_dim=32, n_types=7)

    embeddings = torch.randn(20, 32)
    soft_assignments = torch.softmax(torch.randn(20, 7), dim=1)

    pred_props, entropy, attn_weights_list = agg(embeddings, soft_assignments)

    assert pred_props.shape == (7,)
    assert len(attn_weights_list) == 7
    for w in attn_weights_list:
        assert w.shape == (20,)


def test_entropy_regularization():
    """Entropy should be higher for uniform attention, lower for peaked."""
    from CITEgeist.model.attention_aggregator import AttentionAggregator

    agg = AttentionAggregator(embed_dim=32, n_types=7)

    # Create embeddings that will produce different attention patterns
    # Uniform embeddings -> more uniform attention
    uniform_embed = torch.zeros(20, 32)
    soft_assignments = torch.softmax(torch.randn(20, 7), dim=1)

    _, entropy_uniform, _ = agg(uniform_embed, soft_assignments)

    # Varied embeddings -> potentially more peaked attention
    varied_embed = torch.randn(20, 32) * 10
    _, entropy_varied, _ = agg(varied_embed, soft_assignments)

    # Both should produce valid entropy values
    assert entropy_uniform.item() >= 0
    assert entropy_varied.item() >= 0
```

**Step 2: Run test to verify it fails**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest CITEgeist/tests/test_attention_aggregator.py -v`

Expected: FAIL with "No module named 'CITEgeist.model.attention_aggregator'"

**Step 3: Write minimal implementation**

```python
# CITEgeist/model/attention_aggregator.py
"""Attention-weighted MIL aggregation for single-cell assignment.

Provides attention mechanisms for aggregating per-nucleus predictions
to spot-level proportions, allowing the model to downweight ambiguous
or uninformative nuclei.

Two variants:
1. AttentionAggregator: Single attention head, shared across all classes
2. PerClassAttentionAggregator: K separate attention heads (MoE style)
"""
import torch
import torch.nn as nn
import torch.nn.functional as F
from typing import Tuple, List


class AttentionAggregator(nn.Module):
    """Attention-weighted aggregation for MIL.

    Learns to weight nuclei by their informativeness, then computes
    weighted average of soft assignments to get spot-level proportions.

    Attributes:
        attention: Gating network that produces attention weights
        n_types: Number of cell types
    """

    def __init__(
        self,
        embed_dim: int,
        n_types: int,
        hidden_dim: int = 64,
    ):
        """Initialize attention aggregator.

        Args:
            embed_dim: Dimension of nucleus embeddings
            n_types: Number of cell types (K)
            hidden_dim: Hidden dimension in attention network
        """
        super().__init__()
        self.n_types = n_types

        self.attention = nn.Sequential(
            nn.Linear(embed_dim, hidden_dim),
            nn.Tanh(),
            nn.Linear(hidden_dim, 1),
        )

    def forward(
        self,
        embeddings: torch.Tensor,
        soft_assignments: torch.Tensor,
    ) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        """Aggregate soft assignments with attention weighting.

        Args:
            embeddings: (N, D) nucleus embeddings
            soft_assignments: (N, K) per-nucleus soft type assignments

        Returns:
            pred_props: (K,) predicted spot-level proportions
            entropy: Scalar, entropy of attention weights (for regularization)
            attn_weights: (N, 1) attention weights per nucleus
        """
        N = embeddings.shape[0]

        # Compute attention logits and weights
        attn_logits = self.attention(embeddings)  # (N, 1)
        attn_weights = F.softmax(attn_logits, dim=0)  # (N, 1), sum to 1

        # Weighted aggregation of soft assignments
        weighted_assignments = attn_weights * soft_assignments  # (N, K)
        pred_props = weighted_assignments.sum(dim=0)  # (K,)

        # Entropy of attention weights (higher = more uniform)
        eps = 1e-8
        entropy = -(attn_weights * torch.log(attn_weights + eps)).sum()

        return pred_props, entropy, attn_weights


class PerClassAttentionAggregator(nn.Module):
    """Per-class attention aggregation (Mixture-of-Experts MIL).

    Each cell type has its own attention head, allowing the model to
    learn what features are important for each specific type. This helps
    rare classes that might be overwhelmed by majority types.

    Attributes:
        attentions: List of K attention networks (one per class)
        n_types: Number of cell types
    """

    def __init__(
        self,
        embed_dim: int,
        n_types: int,
        hidden_dim: int = 64,
    ):
        """Initialize per-class attention aggregator.

        Args:
            embed_dim: Dimension of nucleus embeddings
            n_types: Number of cell types (K)
            hidden_dim: Hidden dimension in each attention network
        """
        super().__init__()
        self.n_types = n_types

        self.attentions = nn.ModuleList([
            nn.Sequential(
                nn.Linear(embed_dim, hidden_dim),
                nn.Tanh(),
                nn.Linear(hidden_dim, 1),
            )
            for _ in range(n_types)
        ])

    def forward(
        self,
        embeddings: torch.Tensor,
        soft_assignments: torch.Tensor,
    ) -> Tuple[torch.Tensor, torch.Tensor, List[torch.Tensor]]:
        """Aggregate with per-class attention.

        Args:
            embeddings: (N, D) nucleus embeddings
            soft_assignments: (N, K) per-nucleus soft type assignments

        Returns:
            pred_props: (K,) predicted spot-level proportions
            mean_entropy: Scalar, mean entropy across all attention heads
            attn_weights_list: List of K tensors, each (N,) attention weights
        """
        N = embeddings.shape[0]
        pred_props = []
        entropies = []
        attn_weights_list = []

        for k in range(self.n_types):
            # Compute attention for this class
            attn_logits = self.attentions[k](embeddings)  # (N, 1)
            attn_weights = F.softmax(attn_logits, dim=0).squeeze(-1)  # (N,)

            # Weighted sum of this class's soft assignments
            weighted_k = (attn_weights * soft_assignments[:, k]).sum()
            pred_props.append(weighted_k)

            # Entropy for this head
            eps = 1e-8
            entropy_k = -(attn_weights * torch.log(attn_weights + eps)).sum()
            entropies.append(entropy_k)

            attn_weights_list.append(attn_weights)

        # Stack and normalize proportions
        pred_props = torch.stack(pred_props)  # (K,)
        pred_props = pred_props / (pred_props.sum() + 1e-8)  # Normalize to sum=1

        mean_entropy = torch.stack(entropies).mean()

        return pred_props, mean_entropy, attn_weights_list
```

**Step 4: Run test to verify it passes**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest CITEgeist/tests/test_attention_aggregator.py -v`

Expected: PASS (4 tests)

**Step 5: Commit**

```bash
git add CITEgeist/model/attention_aggregator.py CITEgeist/tests/test_attention_aggregator.py
git commit -m "feat: add attention-weighted MIL aggregation modules"
```

---

## Task 5: Integrate Attention into DirectSoftmax Model

**Files:**
- Modify: `CITEgeist/model/direct_softmax_model.py`
- Test: `CITEgeist/tests/test_direct_softmax_attention.py`

**Step 1: Write the failing test**

```python
# CITEgeist/tests/test_direct_softmax_attention.py
"""Tests for attention-enhanced DirectSoftmax model."""
import torch
import pytest


class MockEncoder:
    """Mock VAE encoder for testing."""
    def __init__(self):
        pass

    def __call__(self, x):
        B = x.shape[0]
        mu = torch.randn(B, 128)
        logvar = torch.zeros(B, 128)
        return mu, logvar

    def parameters(self):
        return iter([])

    def eval(self):
        pass

    def to(self, device):
        return self


def test_direct_softmax_with_attention():
    """DirectSoftmax with attention should produce valid proportions."""
    from CITEgeist.model.direct_softmax_model import DirectSoftmaxModel

    encoder = MockEncoder()
    model = DirectSoftmaxModel(
        encoder=encoder,
        n_types=7,
        latent_dim=128,
        use_attention=True,
    )

    patches = torch.randn(20, 2, 96, 96)
    proportions = torch.softmax(torch.randn(7), dim=0)

    loss, soft_assignments = model(patches, proportions)

    assert loss.dim() == 0  # Scalar
    assert soft_assignments.shape == (20, 7)


def test_direct_softmax_with_per_class_attention():
    """DirectSoftmax with per-class attention should work."""
    from CITEgeist.model.direct_softmax_model import DirectSoftmaxModel

    encoder = MockEncoder()
    model = DirectSoftmaxModel(
        encoder=encoder,
        n_types=7,
        latent_dim=128,
        use_attention=True,
        use_per_class_attention=True,
    )

    patches = torch.randn(20, 2, 96, 96)
    proportions = torch.softmax(torch.randn(7), dim=0)

    loss, soft_assignments = model(patches, proportions)

    assert loss.dim() == 0
    assert soft_assignments.shape == (20, 7)


def test_direct_softmax_attention_entropy_in_loss():
    """Entropy regularization should affect total loss."""
    from CITEgeist.model.direct_softmax_model import DirectSoftmaxModel

    encoder = MockEncoder()
    model = DirectSoftmaxModel(
        encoder=encoder,
        n_types=7,
        latent_dim=128,
        use_attention=True,
        attention_entropy_weight=0.1,
    )

    patches = torch.randn(20, 2, 96, 96)
    proportions = torch.softmax(torch.randn(7), dim=0)

    loss_components, _ = model(patches, proportions, return_components=True)

    assert "attention_entropy" in loss_components
```

**Step 2: Run test to verify it fails**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest CITEgeist/tests/test_direct_softmax_attention.py -v`

Expected: FAIL with "TypeError: DirectSoftmaxModel.__init__() got an unexpected keyword argument 'use_attention'"

**Step 3: Modify DirectSoftmaxModel to add attention support**

Add import at top of `direct_softmax_model.py` (after line 36):
```python
try:
    from .attention_aggregator import AttentionAggregator, PerClassAttentionAggregator
except ImportError:
    from attention_aggregator import AttentionAggregator, PerClassAttentionAggregator
```

Modify `__init__` signature (lines 53-64) to add new parameters:
```python
    def __init__(
        self,
        encoder: nn.Module,
        n_types: int,
        latent_dim: int = 128,
        projection_dim: int = 32,
        temperature: float = 0.1,
        repulsion_weight: float = 1.0,
        repulsion_margin: float = 0.5,
        entropy_weight: float = 0.1,
        use_cosine: bool = True,
        # Attention parameters
        use_attention: bool = False,
        use_per_class_attention: bool = False,
        attention_entropy_weight: float = 0.1,
    ):
```

Add attention initialization after line 100 (after self.use_cosine = use_cosine):
```python
        # Attention aggregation
        self.use_attention = use_attention
        self.use_per_class_attention = use_per_class_attention
        self.attention_entropy_weight = attention_entropy_weight

        if use_attention:
            if use_per_class_attention:
                self.aggregator = PerClassAttentionAggregator(
                    embed_dim=projection_dim,
                    n_types=n_types,
                )
            else:
                self.aggregator = AttentionAggregator(
                    embed_dim=projection_dim,
                    n_types=n_types,
                )
```

Modify forward method (lines 139-199) to use attention when enabled. Replace the aggregation section (around line 172-173):
```python
        # Softmax to get per-nucleus soft assignments
        soft_assignments = F.softmax(logits / self.temperature, dim=1)  # (N, K)

        # Aggregate to spot-level predicted proportions
        if self.use_attention:
            pred_props, attention_entropy, _ = self.aggregator(projected, soft_assignments)
        else:
            pred_props = soft_assignments.mean(dim=0)  # (K,)
            attention_entropy = torch.tensor(0.0, device=logits.device)
```

Update the loss computation (around line 194-198):
```python
        # Total loss
        loss = (
            kl_loss
            + self.repulsion_weight * repulsion_loss
            - self.entropy_weight * assignment_entropy  # Negative = encourage entropy
        )

        # Add attention entropy regularization if using attention
        if self.use_attention:
            loss = loss - self.attention_entropy_weight * attention_entropy
```

Update the return_components dict (around line 218-228):
```python
        if return_components:
            components = {
                "total": loss,
                "kl": kl_loss,
                "repulsion": repulsion_loss,
                "entropy": assignment_entropy,
                "prop_mae": prop_mae,
                "confidence": confidence,
                "proto_min_dist": proto_min_dist,
                "logit_std": logit_std,
            }
            if self.use_attention:
                components["attention_entropy"] = attention_entropy
            return components, soft_assignments
```

**Step 4: Run test to verify it passes**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest CITEgeist/tests/test_direct_softmax_attention.py -v`

Expected: PASS (3 tests)

**Step 5: Commit**

```bash
git add CITEgeist/model/direct_softmax_model.py CITEgeist/tests/test_direct_softmax_attention.py
git commit -m "feat: integrate attention aggregation into DirectSoftmax"
```

---

## Task 6: Update Training Script with Attention Flags

**Files:**
- Modify: `CITEgeist/model/train_prototypes.py`

**Step 1: Add CLI arguments for attention**

Add after line 165 (after the existing DirectSoftmax args):
```python
    # Attention aggregation parameters
    parser.add_argument(
        "--use-attention",
        action="store_true",
        help="Use attention-weighted MIL aggregation instead of mean pooling"
    )
    parser.add_argument(
        "--use-per-class-attention",
        action="store_true",
        help="Use per-class (MoE) attention heads"
    )
    parser.add_argument(
        "--attention-entropy-weight",
        type=float,
        default=0.1,
        help="Weight for attention entropy regularization"
    )
```

**Step 2: Pass attention args when creating DirectSoftmaxModel**

Find the DirectSoftmaxModel initialization in train_prototypes.py and add:
```python
            use_attention=args.use_attention,
            use_per_class_attention=args.use_per_class_attention,
            attention_entropy_weight=args.attention_entropy_weight,
```

**Step 3: Update logging to include attention metrics**

When logging training progress, add attention entropy if applicable.

**Step 4: Test manually**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m CITEgeist.model.train_prototypes --help | grep -i attention`

Expected: Should show the three new attention arguments

**Step 5: Commit**

```bash
git add CITEgeist/model/train_prototypes.py
git commit -m "feat: add attention CLI flags to prototype training"
```

---

## Task 7: Create SLURM Scripts for Training Pipeline

**Files:**
- Create: `Benchmarking/xenium_benchmarking/CITEgeist/slurm/sbatch_vae_vicreg.sh`
- Create: `Benchmarking/xenium_benchmarking/CITEgeist/slurm/sbatch_prototype_attention.sh`

**Step 1: Write VAE + VICReg training script**

```bash
#!/bin/bash
#SBATCH --job-name=vae_vicreg
#SBATCH --output=logs/vae_vicreg_%j.out
#SBATCH --error=logs/vae_vicreg_%j.err
#SBATCH --partition=l40s
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

# CITEgeist environment
source /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env/bin/activate

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist

OUTPUT_DIR="Benchmarking/xenium_benchmarking/CITEgeist/output/vae_vicreg"
mkdir -p "${OUTPUT_DIR}/vae"

python -m CITEgeist.model.train_vae \
    --patches-dir "Benchmarking/xenium_benchmarking/CITEgeist/output/vae_sinkhorn_2ch/spot_patches" \
    --output-dir "${OUTPUT_DIR}/vae" \
    --in-channels 2 \
    --latent-dim 128 \
    --epochs 100 \
    --batch-size 64 \
    --lr 1e-4 \
    --beta 0.5 \
    --use-vicreg \
    --vicreg-weight 1.0 \
    --vicreg-invariance 25.0 \
    --vicreg-variance 25.0 \
    --vicreg-covariance 1.0

echo "VAE + VICReg training complete"
```

**Step 2: Write prototype + attention training script**

```bash
#!/bin/bash
#SBATCH --job-name=proto_attn
#SBATCH --output=logs/proto_attn_%j.out
#SBATCH --error=logs/proto_attn_%j.err
#SBATCH --partition=l40s
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

# CITEgeist environment
source /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env/bin/activate

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist

OUTPUT_DIR="Benchmarking/xenium_benchmarking/CITEgeist/output/vae_vicreg"

python -m CITEgeist.model.train_prototypes \
    --vae-checkpoint "${OUTPUT_DIR}/vae/vae_final.pt" \
    --patches-dir "${OUTPUT_DIR}/../vae_sinkhorn_2ch/spot_patches" \
    --proportions-file "${OUTPUT_DIR}/../vae_sinkhorn_2ch/proportions_for_training.csv" \
    --output-dir "${OUTPUT_DIR}/prototypes_attention" \
    --n-types 7 \
    --use-direct-softmax \
    --softmax-temperature 0.1 \
    --repulsion-weight 1.0 \
    --repulsion-margin 0.5 \
    --entropy-weight 0.01 \
    --use-attention \
    --use-per-class-attention \
    --attention-entropy-weight 0.1 \
    --epochs 50 \
    --lr 1e-3

echo "Prototype + Attention training complete"
```

**Step 3: Commit**

```bash
git add Benchmarking/xenium_benchmarking/CITEgeist/slurm/sbatch_vae_vicreg.sh
git add Benchmarking/xenium_benchmarking/CITEgeist/slurm/sbatch_prototype_attention.sh
git commit -m "feat: add SLURM scripts for VICReg + attention training"
```

---

## Task 8: Run Full Pipeline and Evaluate

**Step 1: Submit VAE + VICReg training**

```bash
cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/CITEgeist/slurm
mkdir -p logs
sbatch sbatch_vae_vicreg.sh
```

**Step 2: After VAE completes, submit prototype + attention training**

```bash
sbatch sbatch_prototype_attention.sh
```

**Step 3: Run evaluation**

Use existing evaluation script with new checkpoint paths:
```bash
python -m CITEgeist.model.evaluate_direct_softmax \
    --vae-checkpoint output/vae_vicreg/vae/vae_final.pt \
    --prototype-checkpoint output/vae_vicreg/prototypes_attention/prototypes_final.pt \
    --output-dir output/vae_vicreg/eval
```

**Step 4: Compare results**

| Method | Single-cell Acc | Notes |
|--------|-----------------|-------|
| DirectSoftmax (baseline) | 20.2% | No VICReg, mean pooling |
| VICReg + Attention | TBD | This experiment |

---

## Summary

This plan implements three improvements to the VAE + DirectSoftmax pipeline:

1. **VICReg VAE pretraining** (Tasks 1-3): Adds variance/invariance/covariance regularization with geometric augmentations to learn discriminative latents

2. **Attention-weighted MIL** (Tasks 4-5): Replaces mean pooling with learned attention to downweight ambiguous nuclei

3. **Per-class attention** (Task 5): K separate attention heads for class-specific focus

All features are controlled by CLI flags for easy ablation.
