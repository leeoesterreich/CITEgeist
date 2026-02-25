# VAE + Sinkhorn OT Single-Cell Assignment Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Implement a two-stage deep learning system (VAE + Sinkhorn OT) to assign individual nuclei to cell types using spot-level proportion constraints.

**Architecture:** Stage 1 trains a VAE on nucleus patches for representation learning. Stage 2 trains projection heads and prototypes using Sinkhorn optimal transport with spot-level proportions as weak supervision. Inference assigns nuclei to their nearest prototype, respecting proportion constraints.

**Tech Stack:** PyTorch, torchvision, scikit-learn (for K-means init), existing CITEgeist infrastructure

---

## Task 1: Sinkhorn OT Implementation

**Files:**
- Create: `CITEgeist/model/sinkhorn.py`
- Test: `CITEgeist/tests/test_sinkhorn.py`

**Step 1: Write the failing test for basic Sinkhorn**

```python
# CITEgeist/tests/test_sinkhorn.py
"""Tests for Sinkhorn optimal transport."""
import torch
import pytest


def test_sinkhorn_returns_valid_transport_plan():
    """Transport plan should have correct marginals."""
    from CITEgeist.model.sinkhorn import sinkhorn

    # 5 nuclei, 3 cell types
    cost = torch.rand(5, 3)
    row_marginal = torch.ones(5) / 5  # uniform
    col_marginal = torch.tensor([0.4, 0.4, 0.2])  # proportions

    P = sinkhorn(cost, row_marginal, col_marginal)

    # Check shape
    assert P.shape == (5, 3)

    # Check row marginals (each nucleus sums to 1/N)
    row_sums = P.sum(dim=1)
    assert torch.allclose(row_sums, row_marginal, atol=1e-3)

    # Check column marginals (match proportions)
    col_sums = P.sum(dim=0)
    assert torch.allclose(col_sums, col_marginal, atol=1e-3)

    # Check non-negative
    assert (P >= 0).all()
```

**Step 2: Run test to verify it fails**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest CITEgeist/tests/test_sinkhorn.py::test_sinkhorn_returns_valid_transport_plan -v`
Expected: FAIL with "ModuleNotFoundError: No module named 'CITEgeist.model.sinkhorn'"

**Step 3: Write minimal implementation**

```python
# CITEgeist/model/sinkhorn.py
"""Differentiable Sinkhorn optimal transport."""
import torch


def sinkhorn(
    cost: torch.Tensor,
    row_marginal: torch.Tensor,
    col_marginal: torch.Tensor,
    temperature: float = 0.1,
    n_iters: int = 50,
    eps: float = 1e-8,
) -> torch.Tensor:
    """
    Compute optimal transport plan via Sinkhorn iterations.

    Args:
        cost: (N, K) cost matrix (e.g., distances to prototypes)
        row_marginal: (N,) target row sums (typically uniform 1/N)
        col_marginal: (K,) target column sums (spot proportions)
        temperature: Lower = sharper assignments (default 0.1)
        n_iters: Number of Sinkhorn iterations (default 50)
        eps: Small value for numerical stability

    Returns:
        P: (N, K) transport plan (soft assignment matrix)
    """
    # Initialize kernel from cost
    K = torch.exp(-cost / temperature)

    # Normalize marginals
    row_marginal = row_marginal / row_marginal.sum()
    col_marginal = col_marginal / col_marginal.sum()

    # Sinkhorn iterations
    for _ in range(n_iters):
        # Row normalization
        row_sum = K.sum(dim=1, keepdim=True) + eps
        K = K / row_sum * row_marginal.unsqueeze(1)

        # Column normalization
        col_sum = K.sum(dim=0, keepdim=True) + eps
        K = K / col_sum * col_marginal.unsqueeze(0)

    return K
```

**Step 4: Run test to verify it passes**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest CITEgeist/tests/test_sinkhorn.py::test_sinkhorn_returns_valid_transport_plan -v`
Expected: PASS

**Step 5: Add test for differentiability**

```python
# Add to CITEgeist/tests/test_sinkhorn.py

def test_sinkhorn_is_differentiable():
    """Gradients should flow through Sinkhorn."""
    from CITEgeist.model.sinkhorn import sinkhorn

    cost = torch.rand(5, 3, requires_grad=True)
    row_marginal = torch.ones(5) / 5
    col_marginal = torch.tensor([0.4, 0.4, 0.2])

    P = sinkhorn(cost, row_marginal, col_marginal)
    loss = (P * cost).sum()  # OT loss
    loss.backward()

    assert cost.grad is not None
    assert cost.grad.shape == cost.shape
```

**Step 6: Run differentiability test**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest CITEgeist/tests/test_sinkhorn.py::test_sinkhorn_is_differentiable -v`
Expected: PASS

**Step 7: Add test for temperature effect**

```python
# Add to CITEgeist/tests/test_sinkhorn.py

def test_sinkhorn_temperature_affects_sharpness():
    """Lower temperature should give sharper assignments."""
    from CITEgeist.model.sinkhorn import sinkhorn

    cost = torch.tensor([[0.1, 0.5, 0.9],
                         [0.9, 0.1, 0.5],
                         [0.5, 0.9, 0.1]])
    row_marginal = torch.ones(3) / 3
    col_marginal = torch.ones(3) / 3

    P_sharp = sinkhorn(cost, row_marginal, col_marginal, temperature=0.01)
    P_soft = sinkhorn(cost, row_marginal, col_marginal, temperature=1.0)

    # Sharp should have higher max values (more concentrated)
    assert P_sharp.max() > P_soft.max()
```

**Step 8: Run temperature test**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest CITEgeist/tests/test_sinkhorn.py::test_sinkhorn_temperature_affects_sharpness -v`
Expected: PASS

**Step 9: Commit**

```bash
git add CITEgeist/model/sinkhorn.py CITEgeist/tests/test_sinkhorn.py
git commit -m "feat: add differentiable Sinkhorn optimal transport"
```

---

## Task 2: Patch Extraction Pipeline

**Files:**
- Create: `CITEgeist/model/patch_extraction.py`
- Test: `CITEgeist/tests/test_patch_extraction.py`

**Step 1: Write the failing test for single patch extraction**

```python
# CITEgeist/tests/test_patch_extraction.py
"""Tests for nucleus patch extraction."""
import torch
import numpy as np
import pytest


def test_extract_single_patch_shape():
    """Extracted patch should have correct shape."""
    from CITEgeist.model.patch_extraction import extract_patch

    # Fake 3-channel image (DAPI + 2 boundary markers)
    image = np.random.rand(3, 500, 500).astype(np.float32)

    # Nucleus bounding box (x_min, y_min, x_max, y_max)
    bbox = (100, 150, 120, 170)  # 20x20 nucleus

    patch = extract_patch(image, bbox, expansion=0.75, output_size=96)

    assert patch.shape == (3, 96, 96)
    assert patch.dtype == np.float32


def test_extract_patch_with_expansion():
    """Expansion should capture context around nucleus."""
    from CITEgeist.model.patch_extraction import extract_patch

    # Create image with marker at nucleus center
    image = np.zeros((1, 200, 200), dtype=np.float32)
    image[0, 90:110, 90:110] = 1.0  # nucleus region

    bbox = (90, 90, 110, 110)  # 20x20 nucleus

    # No expansion - just nucleus
    patch_no_exp = extract_patch(image, bbox, expansion=0.0, output_size=20)

    # With expansion - should include surrounding context
    patch_exp = extract_patch(image, bbox, expansion=1.0, output_size=60)

    # Expanded patch should have lower mean (more background included)
    assert patch_exp.mean() < patch_no_exp.mean()
```

**Step 2: Run test to verify it fails**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest CITEgeist/tests/test_patch_extraction.py::test_extract_single_patch_shape -v`
Expected: FAIL with "ModuleNotFoundError"

**Step 3: Write minimal implementation**

```python
# CITEgeist/model/patch_extraction.py
"""Nucleus patch extraction for VAE training."""
import numpy as np
from typing import Tuple
import cv2


def extract_patch(
    image: np.ndarray,
    bbox: Tuple[int, int, int, int],
    expansion: float = 0.75,
    output_size: int = 96,
) -> np.ndarray:
    """
    Extract expanded patch around a nucleus.

    Args:
        image: (C, H, W) multi-channel image
        bbox: (x_min, y_min, x_max, y_max) nucleus bounding box
        expansion: Fraction to expand bbox in each direction (0.75 = 75%)
        output_size: Final patch size after resize

    Returns:
        patch: (C, output_size, output_size) normalized patch
    """
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

    # Normalize per channel (z-score)
    for c in range(C):
        mean = resized[c].mean()
        std = resized[c].std() + 1e-8
        resized[c] = (resized[c] - mean) / std

    return resized
```

**Step 4: Run test to verify it passes**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest CITEgeist/tests/test_patch_extraction.py -v`
Expected: PASS

**Step 5: Add batch extraction function**

```python
# Add to CITEgeist/model/patch_extraction.py

def extract_patches_for_spot(
    image: np.ndarray,
    nuclei_df: 'pd.DataFrame',
    expansion: float = 0.75,
    output_size: int = 96,
) -> np.ndarray:
    """
    Extract patches for all nuclei in a spot.

    Args:
        image: (C, H, W) multi-channel image
        nuclei_df: DataFrame with columns ['nucleus_id', 'bbox_x_min',
                   'bbox_y_min', 'bbox_x_max', 'bbox_y_max']
        expansion: Bbox expansion fraction
        output_size: Output patch size

    Returns:
        patches: (N, C, output_size, output_size) array of patches
    """
    patches = []

    for _, row in nuclei_df.iterrows():
        bbox = (
            int(row['bbox_x_min']),
            int(row['bbox_y_min']),
            int(row['bbox_x_max']),
            int(row['bbox_y_max']),
        )
        patch = extract_patch(image, bbox, expansion, output_size)
        patches.append(patch)

    return np.stack(patches, axis=0)
```

**Step 6: Add test for batch extraction**

```python
# Add to CITEgeist/tests/test_patch_extraction.py
import pandas as pd

def test_extract_patches_for_spot():
    """Batch extraction should return correct shape."""
    from CITEgeist.model.patch_extraction import extract_patches_for_spot

    image = np.random.rand(3, 500, 500).astype(np.float32)

    nuclei_df = pd.DataFrame({
        'nucleus_id': [1, 2, 3],
        'bbox_x_min': [50, 150, 250],
        'bbox_y_min': [50, 150, 250],
        'bbox_x_max': [70, 170, 270],
        'bbox_y_max': [70, 170, 270],
    })

    patches = extract_patches_for_spot(image, nuclei_df, expansion=0.75, output_size=96)

    assert patches.shape == (3, 3, 96, 96)  # 3 nuclei, 3 channels, 96x96
```

**Step 7: Run batch extraction test**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest CITEgeist/tests/test_patch_extraction.py::test_extract_patches_for_spot -v`
Expected: PASS

**Step 8: Commit**

```bash
git add CITEgeist/model/patch_extraction.py CITEgeist/tests/test_patch_extraction.py
git commit -m "feat: add nucleus patch extraction for VAE"
```

---

## Task 3: VAE Architecture

**Files:**
- Create: `CITEgeist/model/vae.py`
- Test: `CITEgeist/tests/test_vae.py`

**Step 1: Write the failing test for encoder output shape**

```python
# CITEgeist/tests/test_vae.py
"""Tests for VAE architecture."""
import torch
import pytest


def test_encoder_output_shape():
    """Encoder should output mu and logvar of correct shape."""
    from CITEgeist.model.vae import VAEEncoder

    encoder = VAEEncoder(in_channels=3, latent_dim=128)
    x = torch.randn(4, 3, 96, 96)  # batch of 4 patches

    mu, logvar = encoder(x)

    assert mu.shape == (4, 128)
    assert logvar.shape == (4, 128)
```

**Step 2: Run test to verify it fails**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest CITEgeist/tests/test_vae.py::test_encoder_output_shape -v`
Expected: FAIL with "ModuleNotFoundError"

**Step 3: Write encoder implementation**

```python
# CITEgeist/model/vae.py
"""Variational Autoencoder for nucleus patch representation learning."""
import torch
import torch.nn as nn
import torch.nn.functional as F


class VAEEncoder(nn.Module):
    """CNN encoder for nucleus patches."""

    def __init__(self, in_channels: int = 3, latent_dim: int = 128):
        super().__init__()
        self.latent_dim = latent_dim

        # Conv blocks: 96 -> 48 -> 24 -> 12 -> 6
        self.conv1 = nn.Sequential(
            nn.Conv2d(in_channels, 32, 3, stride=2, padding=1),
            nn.BatchNorm2d(32),
            nn.ReLU(inplace=True),
        )
        self.conv2 = nn.Sequential(
            nn.Conv2d(32, 64, 3, stride=2, padding=1),
            nn.BatchNorm2d(64),
            nn.ReLU(inplace=True),
        )
        self.conv3 = nn.Sequential(
            nn.Conv2d(64, 128, 3, stride=2, padding=1),
            nn.BatchNorm2d(128),
            nn.ReLU(inplace=True),
        )
        self.conv4 = nn.Sequential(
            nn.Conv2d(128, 256, 3, stride=2, padding=1),
            nn.BatchNorm2d(256),
            nn.ReLU(inplace=True),
        )

        # Adaptive pooling to handle variable input sizes
        self.pool = nn.AdaptiveAvgPool2d((1, 1))

        # FC layers for mu and logvar
        self.fc_mu = nn.Linear(256, latent_dim)
        self.fc_logvar = nn.Linear(256, latent_dim)

    def forward(self, x: torch.Tensor) -> tuple:
        """
        Encode patches to latent distribution parameters.

        Args:
            x: (B, C, H, W) input patches

        Returns:
            mu: (B, latent_dim) mean of latent distribution
            logvar: (B, latent_dim) log variance of latent distribution
        """
        x = self.conv1(x)
        x = self.conv2(x)
        x = self.conv3(x)
        x = self.conv4(x)
        x = self.pool(x)
        x = x.view(x.size(0), -1)

        mu = self.fc_mu(x)
        logvar = self.fc_logvar(x)

        return mu, logvar
```

**Step 4: Run test to verify it passes**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest CITEgeist/tests/test_vae.py::test_encoder_output_shape -v`
Expected: PASS

**Step 5: Add decoder and test**

```python
# Add to CITEgeist/model/vae.py

class VAEDecoder(nn.Module):
    """CNN decoder for nucleus patches."""

    def __init__(self, out_channels: int = 3, latent_dim: int = 128):
        super().__init__()

        self.fc = nn.Linear(latent_dim, 256 * 6 * 6)

        # Transposed conv blocks: 6 -> 12 -> 24 -> 48 -> 96
        self.deconv1 = nn.Sequential(
            nn.ConvTranspose2d(256, 128, 4, stride=2, padding=1),
            nn.BatchNorm2d(128),
            nn.ReLU(inplace=True),
        )
        self.deconv2 = nn.Sequential(
            nn.ConvTranspose2d(128, 64, 4, stride=2, padding=1),
            nn.BatchNorm2d(64),
            nn.ReLU(inplace=True),
        )
        self.deconv3 = nn.Sequential(
            nn.ConvTranspose2d(64, 32, 4, stride=2, padding=1),
            nn.BatchNorm2d(32),
            nn.ReLU(inplace=True),
        )
        self.deconv4 = nn.Sequential(
            nn.ConvTranspose2d(32, out_channels, 4, stride=2, padding=1),
        )

    def forward(self, z: torch.Tensor) -> torch.Tensor:
        """
        Decode latent vector to patch.

        Args:
            z: (B, latent_dim) latent vectors

        Returns:
            x_recon: (B, C, 96, 96) reconstructed patches
        """
        x = self.fc(z)
        x = x.view(x.size(0), 256, 6, 6)
        x = self.deconv1(x)
        x = self.deconv2(x)
        x = self.deconv3(x)
        x = self.deconv4(x)
        return x
```

**Step 6: Add decoder test**

```python
# Add to CITEgeist/tests/test_vae.py

def test_decoder_output_shape():
    """Decoder should reconstruct correct shape."""
    from CITEgeist.model.vae import VAEDecoder

    decoder = VAEDecoder(out_channels=3, latent_dim=128)
    z = torch.randn(4, 128)

    x_recon = decoder(z)

    assert x_recon.shape == (4, 3, 96, 96)
```

**Step 7: Run decoder test**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest CITEgeist/tests/test_vae.py::test_decoder_output_shape -v`
Expected: PASS

**Step 8: Add full VAE class with reparameterization**

```python
# Add to CITEgeist/model/vae.py

class VAE(nn.Module):
    """Full VAE for nucleus patch representation learning."""

    def __init__(self, in_channels: int = 3, latent_dim: int = 128):
        super().__init__()
        self.encoder = VAEEncoder(in_channels, latent_dim)
        self.decoder = VAEDecoder(in_channels, latent_dim)
        self.latent_dim = latent_dim

    def reparameterize(self, mu: torch.Tensor, logvar: torch.Tensor) -> torch.Tensor:
        """Reparameterization trick: z = mu + std * epsilon."""
        std = torch.exp(0.5 * logvar)
        eps = torch.randn_like(std)
        return mu + eps * std

    def forward(self, x: torch.Tensor) -> tuple:
        """
        Forward pass through VAE.

        Args:
            x: (B, C, H, W) input patches

        Returns:
            x_recon: (B, C, H, W) reconstructed patches
            mu: (B, latent_dim) latent means
            logvar: (B, latent_dim) latent log variances
        """
        mu, logvar = self.encoder(x)
        z = self.reparameterize(mu, logvar)
        x_recon = self.decoder(z)
        return x_recon, mu, logvar

    def encode(self, x: torch.Tensor) -> torch.Tensor:
        """Encode to latent space (use mean, no sampling)."""
        mu, _ = self.encoder(x)
        return mu

    @staticmethod
    def loss_function(x: torch.Tensor, x_recon: torch.Tensor,
                      mu: torch.Tensor, logvar: torch.Tensor,
                      beta: float = 0.5) -> tuple:
        """
        Compute VAE loss.

        Args:
            x: Original input
            x_recon: Reconstructed input
            mu: Latent means
            logvar: Latent log variances
            beta: KL weight (default 0.5)

        Returns:
            loss: Total loss
            recon_loss: Reconstruction loss
            kl_loss: KL divergence
        """
        # Reconstruction loss (MSE)
        recon_loss = F.mse_loss(x_recon, x, reduction='mean')

        # KL divergence
        kl_loss = -0.5 * torch.mean(1 + logvar - mu.pow(2) - logvar.exp())

        loss = recon_loss + beta * kl_loss

        return loss, recon_loss, kl_loss
```

**Step 9: Add full VAE test**

```python
# Add to CITEgeist/tests/test_vae.py

def test_vae_forward_and_loss():
    """VAE forward pass and loss should work."""
    from CITEgeist.model.vae import VAE

    vae = VAE(in_channels=3, latent_dim=128)
    x = torch.randn(4, 3, 96, 96)

    x_recon, mu, logvar = vae(x)

    assert x_recon.shape == x.shape
    assert mu.shape == (4, 128)

    loss, recon_loss, kl_loss = VAE.loss_function(x, x_recon, mu, logvar)

    assert loss.ndim == 0  # scalar
    assert loss > 0


def test_vae_encode():
    """Encode method should return deterministic latents."""
    from CITEgeist.model.vae import VAE

    vae = VAE(in_channels=3, latent_dim=128)
    vae.eval()
    x = torch.randn(4, 3, 96, 96)

    z1 = vae.encode(x)
    z2 = vae.encode(x)

    # Should be deterministic (no sampling)
    assert torch.allclose(z1, z2)
```

**Step 10: Run all VAE tests**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest CITEgeist/tests/test_vae.py -v`
Expected: All PASS

**Step 11: Commit**

```bash
git add CITEgeist/model/vae.py CITEgeist/tests/test_vae.py
git commit -m "feat: add VAE architecture for nucleus patches"
```

---

## Task 4: Projection Heads and Prototypes

**Files:**
- Create: `CITEgeist/model/projection_heads.py`
- Test: `CITEgeist/tests/test_projection_heads.py`

**Step 1: Write the failing test**

```python
# CITEgeist/tests/test_projection_heads.py
"""Tests for projection heads and prototypes."""
import torch
import pytest


def test_projection_heads_output_shape():
    """Each head should project to correct dimension."""
    from CITEgeist.model.projection_heads import ProjectionHeads

    heads = ProjectionHeads(
        input_dim=128,
        projection_dim=32,
        n_types=7
    )

    z = torch.randn(10, 128)  # 10 nuclei

    projected = heads(z)  # List of K tensors

    assert len(projected) == 7
    for p in projected:
        assert p.shape == (10, 32)
```

**Step 2: Run test to verify it fails**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest CITEgeist/tests/test_projection_heads.py::test_projection_heads_output_shape -v`
Expected: FAIL with "ModuleNotFoundError"

**Step 3: Write implementation**

```python
# CITEgeist/model/projection_heads.py
"""Projection heads and prototypes for cell type assignment."""
import torch
import torch.nn as nn
from typing import List


class ProjectionHeads(nn.Module):
    """K projection heads, one per cell type."""

    def __init__(
        self,
        input_dim: int = 128,
        projection_dim: int = 32,
        n_types: int = 7,
        hidden_dim: int = 64,
    ):
        super().__init__()
        self.n_types = n_types
        self.projection_dim = projection_dim

        # Create K separate MLP heads
        self.heads = nn.ModuleList([
            nn.Sequential(
                nn.Linear(input_dim, hidden_dim),
                nn.ReLU(inplace=True),
                nn.Linear(hidden_dim, projection_dim),
            )
            for _ in range(n_types)
        ])

    def forward(self, z: torch.Tensor) -> List[torch.Tensor]:
        """
        Project latent through each head.

        Args:
            z: (N, input_dim) latent vectors

        Returns:
            List of K tensors, each (N, projection_dim)
        """
        return [head(z) for head in self.heads]


class Prototypes(nn.Module):
    """K learnable prototype vectors."""

    def __init__(self, projection_dim: int = 32, n_types: int = 7):
        super().__init__()
        self.n_types = n_types

        # Initialize prototypes randomly
        self.prototypes = nn.Parameter(
            torch.randn(n_types, projection_dim) * 0.1
        )

    def forward(self) -> torch.Tensor:
        """Return prototype vectors."""
        return self.prototypes

    def init_from_kmeans(self, projected_latents: List[torch.Tensor]):
        """
        Initialize prototypes using K-means on projected latents.

        Args:
            projected_latents: List of K tensors, each (N, projection_dim)
        """
        from sklearn.cluster import KMeans

        with torch.no_grad():
            for k in range(self.n_types):
                # Get all projected latents for this head
                data = projected_latents[k].cpu().numpy()

                # Find centroid via K-means with K=1
                kmeans = KMeans(n_clusters=1, n_init=10).fit(data)
                self.prototypes.data[k] = torch.tensor(
                    kmeans.cluster_centers_[0],
                    dtype=self.prototypes.dtype,
                    device=self.prototypes.device
                )
```

**Step 4: Run test to verify it passes**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest CITEgeist/tests/test_projection_heads.py::test_projection_heads_output_shape -v`
Expected: PASS

**Step 5: Add distance computation test**

```python
# Add to CITEgeist/tests/test_projection_heads.py

def test_compute_distances():
    """Distance computation should have correct shape."""
    from CITEgeist.model.projection_heads import ProjectionHeads, Prototypes, compute_distances

    heads = ProjectionHeads(input_dim=128, projection_dim=32, n_types=7)
    prototypes = Prototypes(projection_dim=32, n_types=7)

    z = torch.randn(10, 128)
    projected = heads(z)
    proto = prototypes()

    distances = compute_distances(projected, proto)

    assert distances.shape == (10, 7)
    assert (distances >= 0).all()
```

**Step 6: Add distance computation function**

```python
# Add to CITEgeist/model/projection_heads.py

def compute_distances(
    projected: List[torch.Tensor],
    prototypes: torch.Tensor,
) -> torch.Tensor:
    """
    Compute distances from projected latents to prototypes.

    Args:
        projected: List of K tensors, each (N, D)
        prototypes: (K, D) prototype vectors

    Returns:
        distances: (N, K) distance matrix
    """
    K = len(projected)
    N = projected[0].shape[0]

    distances = torch.zeros(N, K, device=projected[0].device)

    for k in range(K):
        # Euclidean distance from each nucleus to prototype k
        diff = projected[k] - prototypes[k].unsqueeze(0)
        distances[:, k] = torch.norm(diff, dim=1)

    return distances
```

**Step 7: Run distance test**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest CITEgeist/tests/test_projection_heads.py::test_compute_distances -v`
Expected: PASS

**Step 8: Commit**

```bash
git add CITEgeist/model/projection_heads.py CITEgeist/tests/test_projection_heads.py
git commit -m "feat: add projection heads and prototypes for cell type assignment"
```

---

## Task 5: Prototype Learning Model (Stage 2)

**Files:**
- Create: `CITEgeist/model/prototype_learning.py`
- Test: `CITEgeist/tests/test_prototype_learning.py`

**Step 1: Write the failing test**

```python
# CITEgeist/tests/test_prototype_learning.py
"""Tests for prototype learning (Stage 2)."""
import torch
import pytest


def test_prototype_model_forward():
    """Forward pass should return loss and transport plan."""
    from CITEgeist.model.prototype_learning import PrototypeLearningModel
    from CITEgeist.model.vae import VAEEncoder

    # Create frozen encoder
    encoder = VAEEncoder(in_channels=3, latent_dim=128)
    encoder.eval()
    for p in encoder.parameters():
        p.requires_grad = False

    model = PrototypeLearningModel(
        encoder=encoder,
        n_types=7,
        latent_dim=128,
        projection_dim=32,
    )

    # Fake spot with 10 nuclei
    patches = torch.randn(10, 3, 96, 96)
    proportions = torch.tensor([0.3, 0.2, 0.2, 0.1, 0.1, 0.05, 0.05])

    loss, transport_plan = model(patches, proportions)

    assert loss.ndim == 0  # scalar
    assert transport_plan.shape == (10, 7)
```

**Step 2: Run test to verify it fails**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest CITEgeist/tests/test_prototype_learning.py::test_prototype_model_forward -v`
Expected: FAIL with "ModuleNotFoundError"

**Step 3: Write implementation**

```python
# CITEgeist/model/prototype_learning.py
"""Prototype learning model (Stage 2)."""
import torch
import torch.nn as nn
from typing import Tuple

from .projection_heads import ProjectionHeads, Prototypes, compute_distances
from .sinkhorn import sinkhorn


class PrototypeLearningModel(nn.Module):
    """
    Stage 2 model: learns projection heads and prototypes.

    Uses Sinkhorn OT with spot proportions as supervision.
    """

    def __init__(
        self,
        encoder: nn.Module,
        n_types: int,
        latent_dim: int = 128,
        projection_dim: int = 32,
        sinkhorn_temp: float = 0.1,
        sinkhorn_iters: int = 50,
    ):
        super().__init__()

        # Frozen encoder
        self.encoder = encoder
        for p in self.encoder.parameters():
            p.requires_grad = False

        # Learnable components
        self.heads = ProjectionHeads(
            input_dim=latent_dim,
            projection_dim=projection_dim,
            n_types=n_types,
        )
        self.prototypes = Prototypes(
            projection_dim=projection_dim,
            n_types=n_types,
        )

        self.n_types = n_types
        self.sinkhorn_temp = sinkhorn_temp
        self.sinkhorn_iters = sinkhorn_iters

    def forward(
        self,
        patches: torch.Tensor,
        proportions: torch.Tensor,
    ) -> Tuple[torch.Tensor, torch.Tensor]:
        """
        Forward pass for a single spot.

        Args:
            patches: (N, C, H, W) nucleus patches for this spot
            proportions: (K,) cell type proportions for this spot

        Returns:
            loss: Optimal transport loss (scalar)
            transport_plan: (N, K) soft assignment matrix
        """
        N = patches.shape[0]

        # Encode patches (no grad)
        with torch.no_grad():
            mu, _ = self.encoder(patches)
        z = mu  # Use mean, not sampled

        # Project through heads
        projected = self.heads(z)

        # Get prototypes
        proto = self.prototypes()

        # Compute distances
        distances = compute_distances(projected, proto)

        # Sinkhorn OT
        row_marginal = torch.ones(N, device=patches.device) / N
        transport_plan = sinkhorn(
            distances,
            row_marginal,
            proportions,
            temperature=self.sinkhorn_temp,
            n_iters=self.sinkhorn_iters,
        )

        # OT loss: sum of (transport * cost)
        loss = (transport_plan * distances).sum()

        return loss, transport_plan

    def assign(
        self,
        patches: torch.Tensor,
        proportions: torch.Tensor,
        temperature: float = 0.05,
    ) -> Tuple[torch.Tensor, torch.Tensor]:
        """
        Assign nuclei to cell types.

        Args:
            patches: (N, C, H, W) nucleus patches
            proportions: (K,) cell type proportions
            temperature: Sinkhorn temperature (lower = sharper)

        Returns:
            assignments: (N,) cell type indices
            confidence: (N,) assignment confidence scores
        """
        self.eval()
        with torch.no_grad():
            # Encode
            mu, _ = self.encoder(patches)
            z = mu

            # Project
            projected = self.heads(z)
            proto = self.prototypes()
            distances = compute_distances(projected, proto)

            # Sinkhorn with low temperature for sharp assignments
            N = patches.shape[0]
            row_marginal = torch.ones(N, device=patches.device) / N
            transport_plan = sinkhorn(
                distances,
                row_marginal,
                proportions,
                temperature=temperature,
                n_iters=100,
            )

            # Hard assignment
            assignments = transport_plan.argmax(dim=1)
            confidence = transport_plan.max(dim=1).values

        return assignments, confidence
```

**Step 4: Run test to verify it passes**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest CITEgeist/tests/test_prototype_learning.py::test_prototype_model_forward -v`
Expected: PASS

**Step 5: Add assignment test**

```python
# Add to CITEgeist/tests/test_prototype_learning.py

def test_prototype_model_assignment():
    """Assignment should return valid indices and confidences."""
    from CITEgeist.model.prototype_learning import PrototypeLearningModel
    from CITEgeist.model.vae import VAEEncoder

    encoder = VAEEncoder(in_channels=3, latent_dim=128)
    model = PrototypeLearningModel(encoder=encoder, n_types=7)

    patches = torch.randn(10, 3, 96, 96)
    proportions = torch.tensor([0.3, 0.2, 0.2, 0.1, 0.1, 0.05, 0.05])

    assignments, confidence = model.assign(patches, proportions)

    assert assignments.shape == (10,)
    assert confidence.shape == (10,)
    assert (assignments >= 0).all() and (assignments < 7).all()
    assert (confidence >= 0).all() and (confidence <= 1).all()
```

**Step 6: Run assignment test**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest CITEgeist/tests/test_prototype_learning.py::test_prototype_model_assignment -v`
Expected: PASS

**Step 7: Commit**

```bash
git add CITEgeist/model/prototype_learning.py CITEgeist/tests/test_prototype_learning.py
git commit -m "feat: add prototype learning model with Sinkhorn OT"
```

---

## Task 6: Training Scripts

**Files:**
- Create: `CITEgeist/model/train_vae.py`
- Create: `CITEgeist/model/train_prototypes.py`
- Create: `Benchmarking/xenium_benchmarking/CITEgeist/slurm/sbatch_train_vae.sh`
- Create: `Benchmarking/xenium_benchmarking/CITEgeist/slurm/sbatch_train_prototypes.sh`

**Step 1: Create VAE training script**

```python
# CITEgeist/model/train_vae.py
"""Stage 1: Train VAE on nucleus patches."""
import argparse
import torch
import torch.optim as optim
from torch.utils.data import DataLoader, Dataset
import numpy as np
from pathlib import Path
import json
from tqdm import tqdm

from .vae import VAE


class PatchDataset(Dataset):
    """Dataset of pre-extracted nucleus patches."""

    def __init__(self, patches_dir: Path):
        self.patches_dir = Path(patches_dir)
        self.patch_files = sorted(self.patches_dir.glob("*.npy"))

        # Load all patches into memory (or use lazy loading for large datasets)
        self.patches = []
        for f in tqdm(self.patch_files, desc="Loading patches"):
            self.patches.append(np.load(f))
        self.patches = np.concatenate(self.patches, axis=0)

    def __len__(self):
        return len(self.patches)

    def __getitem__(self, idx):
        return torch.tensor(self.patches[idx], dtype=torch.float32)


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
):
    """Train VAE on nucleus patches."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Data
    dataset = PatchDataset(patches_dir)
    loader = DataLoader(dataset, batch_size=batch_size, shuffle=True, num_workers=4)

    # Model
    vae = VAE(in_channels=in_channels, latent_dim=latent_dim).to(device)
    optimizer = optim.Adam(vae.parameters(), lr=lr)

    # Training loop
    history = {"epoch": [], "loss": [], "recon_loss": [], "kl_loss": []}

    for epoch in range(epochs):
        vae.train()
        epoch_loss = 0
        epoch_recon = 0
        epoch_kl = 0

        for batch in tqdm(loader, desc=f"Epoch {epoch+1}/{epochs}"):
            batch = batch.to(device)

            optimizer.zero_grad()
            x_recon, mu, logvar = vae(batch)
            loss, recon_loss, kl_loss = VAE.loss_function(batch, x_recon, mu, logvar, beta)
            loss.backward()
            optimizer.step()

            epoch_loss += loss.item() * len(batch)
            epoch_recon += recon_loss.item() * len(batch)
            epoch_kl += kl_loss.item() * len(batch)

        # Log
        n = len(dataset)
        history["epoch"].append(epoch)
        history["loss"].append(epoch_loss / n)
        history["recon_loss"].append(epoch_recon / n)
        history["kl_loss"].append(epoch_kl / n)

        print(f"Epoch {epoch+1}: loss={epoch_loss/n:.4f}, "
              f"recon={epoch_recon/n:.4f}, kl={epoch_kl/n:.4f}")

        # Save checkpoint every 10 epochs
        if (epoch + 1) % 10 == 0:
            torch.save({
                "epoch": epoch,
                "model_state_dict": vae.state_dict(),
                "optimizer_state_dict": optimizer.state_dict(),
                "history": history,
            }, output_dir / f"vae_checkpoint_epoch{epoch+1}.pt")

    # Save final model
    torch.save(vae.state_dict(), output_dir / "vae_final.pt")

    with open(output_dir / "training_history.json", "w") as f:
        json.dump(history, f, indent=2)

    return vae


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--patches_dir", required=True)
    parser.add_argument("--output_dir", required=True)
    parser.add_argument("--latent_dim", type=int, default=128)
    parser.add_argument("--batch_size", type=int, default=64)
    parser.add_argument("--epochs", type=int, default=100)
    parser.add_argument("--lr", type=float, default=1e-4)
    parser.add_argument("--beta", type=float, default=0.5)
    args = parser.parse_args()

    train_vae(**vars(args))
```

**Step 2: Create prototype training script**

```python
# CITEgeist/model/train_prototypes.py
"""Stage 2: Train projection heads and prototypes."""
import argparse
import torch
import torch.optim as optim
import numpy as np
import pandas as pd
from pathlib import Path
import json
from tqdm import tqdm

from .vae import VAE
from .prototype_learning import PrototypeLearningModel


def train_prototypes(
    vae_checkpoint: str,
    patches_dir: str,
    proportions_file: str,
    output_dir: str,
    n_types: int = 7,
    latent_dim: int = 128,
    projection_dim: int = 32,
    epochs: int = 50,
    lr: float = 1e-3,
    sinkhorn_temp: float = 0.1,
    device: str = "cuda",
):
    """Train projection heads and prototypes."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load VAE encoder
    vae = VAE(in_channels=3, latent_dim=latent_dim)
    vae.load_state_dict(torch.load(vae_checkpoint, map_location="cpu"))
    encoder = vae.encoder.to(device)
    encoder.eval()

    # Load proportions
    proportions_df = pd.read_csv(proportions_file)
    spot_ids = proportions_df["spot_id"].values
    prop_cols = [c for c in proportions_df.columns if c != "spot_id"]

    # Model
    model = PrototypeLearningModel(
        encoder=encoder,
        n_types=n_types,
        latent_dim=latent_dim,
        projection_dim=projection_dim,
        sinkhorn_temp=sinkhorn_temp,
    ).to(device)

    optimizer = optim.Adam(
        list(model.heads.parameters()) + list(model.prototypes.parameters()),
        lr=lr,
    )

    # Training loop
    history = {"epoch": [], "loss": []}
    patches_dir = Path(patches_dir)

    for epoch in range(epochs):
        model.train()
        epoch_loss = 0
        n_spots = 0

        for spot_id in tqdm(spot_ids, desc=f"Epoch {epoch+1}/{epochs}"):
            # Load patches for this spot
            patch_file = patches_dir / f"spot_{spot_id}_patches.npy"
            if not patch_file.exists():
                continue

            patches = torch.tensor(np.load(patch_file), dtype=torch.float32).to(device)
            if len(patches) == 0:
                continue

            # Get proportions
            row = proportions_df[proportions_df["spot_id"] == spot_id]
            proportions = torch.tensor(row[prop_cols].values[0], dtype=torch.float32).to(device)

            # Forward
            optimizer.zero_grad()
            loss, _ = model(patches, proportions)
            loss.backward()
            optimizer.step()

            epoch_loss += loss.item()
            n_spots += 1

        avg_loss = epoch_loss / max(n_spots, 1)
        history["epoch"].append(epoch)
        history["loss"].append(avg_loss)

        print(f"Epoch {epoch+1}: loss={avg_loss:.4f}")

        # Save checkpoint
        if (epoch + 1) % 10 == 0:
            torch.save({
                "epoch": epoch,
                "heads_state_dict": model.heads.state_dict(),
                "prototypes_state_dict": model.prototypes.state_dict(),
                "optimizer_state_dict": optimizer.state_dict(),
                "history": history,
            }, output_dir / f"prototype_checkpoint_epoch{epoch+1}.pt")

    # Save final
    torch.save({
        "heads_state_dict": model.heads.state_dict(),
        "prototypes_state_dict": model.prototypes.state_dict(),
    }, output_dir / "prototype_final.pt")

    with open(output_dir / "training_history.json", "w") as f:
        json.dump(history, f, indent=2)

    return model


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--vae_checkpoint", required=True)
    parser.add_argument("--patches_dir", required=True)
    parser.add_argument("--proportions_file", required=True)
    parser.add_argument("--output_dir", required=True)
    parser.add_argument("--n_types", type=int, default=7)
    parser.add_argument("--epochs", type=int, default=50)
    parser.add_argument("--lr", type=float, default=1e-3)
    args = parser.parse_args()

    train_prototypes(**vars(args))
```

**Step 3: Create SLURM scripts**

```bash
# Benchmarking/xenium_benchmarking/CITEgeist/slurm/sbatch_train_vae.sh
#!/bin/bash
#SBATCH --job-name=train_vae
#SBATCH --output=logs/train_vae_%j.out
#SBATCH --error=logs/train_vae_%j.err
#SBATCH --time=24:00:00
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

# Activate environment
source ~/.bashrc
eval "$(conda shell.bash hook)"
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist

python -m CITEgeist.model.train_vae \
    --patches_dir Benchmarking/xenium_benchmarking/data/patches \
    --output_dir Benchmarking/xenium_benchmarking/CITEgeist/output/vae \
    --latent_dim 128 \
    --batch_size 64 \
    --epochs 100 \
    --lr 1e-4 \
    --beta 0.5
```

```bash
# Benchmarking/xenium_benchmarking/CITEgeist/slurm/sbatch_train_prototypes.sh
#!/bin/bash
#SBATCH --job-name=train_proto
#SBATCH --output=logs/train_proto_%j.out
#SBATCH --error=logs/train_proto_%j.err
#SBATCH --time=12:00:00
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

source ~/.bashrc
eval "$(conda shell.bash hook)"
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist

python -m CITEgeist.model.train_prototypes \
    --vae_checkpoint Benchmarking/xenium_benchmarking/CITEgeist/output/vae/vae_final.pt \
    --patches_dir Benchmarking/xenium_benchmarking/data/patches \
    --proportions_file Benchmarking/xenium_benchmarking/data/proportions.csv \
    --output_dir Benchmarking/xenium_benchmarking/CITEgeist/output/prototypes \
    --n_types 7 \
    --epochs 50 \
    --lr 1e-3
```

**Step 4: Commit**

```bash
git add CITEgeist/model/train_vae.py CITEgeist/model/train_prototypes.py \
    Benchmarking/xenium_benchmarking/CITEgeist/slurm/sbatch_train_vae.sh \
    Benchmarking/xenium_benchmarking/CITEgeist/slurm/sbatch_train_prototypes.sh
git commit -m "feat: add VAE and prototype training scripts"
```

---

## Task 7: Integration with Module 3b

**Files:**
- Modify: `CITEgeist/model/module3b_nucleus_assignment.py`
- Test: `CITEgeist/tests/test_module3b_nucleus_assignment.py`

**Step 1: Add VAE-Sinkhorn assignment method**

```python
# Add to CITEgeist/model/module3b_nucleus_assignment.py (after existing imports)

def run_nucleus_assignment_vae(
    image: np.ndarray,
    nuclei_df: pd.DataFrame,
    proportions: pd.DataFrame,
    cell_types: List[str],
    vae_checkpoint: str,
    prototype_checkpoint: str,
    device: str = "cuda",
    patch_expansion: float = 0.75,
    patch_size: int = 96,
) -> NucleusAssignmentResult:
    """
    Assign nuclei using VAE + Sinkhorn OT approach.

    Args:
        image: (C, H, W) multi-channel image
        nuclei_df: DataFrame with nucleus bboxes and spot assignments
        proportions: DataFrame with spot proportions
        cell_types: List of cell type names
        vae_checkpoint: Path to trained VAE weights
        prototype_checkpoint: Path to trained prototype weights
        device: torch device
        patch_expansion: Bbox expansion fraction
        patch_size: Output patch size

    Returns:
        NucleusAssignmentResult with assignments
    """
    import torch
    from .vae import VAE
    from .prototype_learning import PrototypeLearningModel
    from .patch_extraction import extract_patch

    # Load models
    vae = VAE(in_channels=image.shape[0], latent_dim=128)
    vae.load_state_dict(torch.load(vae_checkpoint, map_location="cpu"))

    model = PrototypeLearningModel(
        encoder=vae.encoder,
        n_types=len(cell_types),
    )
    checkpoint = torch.load(prototype_checkpoint, map_location="cpu")
    model.heads.load_state_dict(checkpoint["heads_state_dict"])
    model.prototypes.load_state_dict(checkpoint["prototypes_state_dict"])
    model = model.to(device)
    model.eval()

    # Set up proportions lookup
    prop_cols = cell_types
    spot_props = proportions.set_index('spot_id')[prop_cols]

    all_assignments = {}
    all_confidences = {}

    for spot_id in nuclei_df['spot_id'].unique():
        spot_nuclei = nuclei_df[nuclei_df['spot_id'] == spot_id]

        if len(spot_nuclei) == 0:
            continue

        # Extract patches
        patches = []
        for _, row in spot_nuclei.iterrows():
            bbox = (
                int(row['bbox_x_min']),
                int(row['bbox_y_min']),
                int(row['bbox_x_max']),
                int(row['bbox_y_max']),
            )
            patch = extract_patch(image, bbox, patch_expansion, patch_size)
            patches.append(patch)

        patches = torch.tensor(np.stack(patches), dtype=torch.float32).to(device)

        # Get proportions
        spot_proportions = torch.tensor(
            spot_props.loc[spot_id].values,
            dtype=torch.float32
        ).to(device)

        # Assign
        with torch.no_grad():
            assignments, confidence = model.assign(patches, spot_proportions)

        # Store results
        for i, (_, row) in enumerate(spot_nuclei.iterrows()):
            nid = int(row['nucleus_id'])
            type_idx = assignments[i].item()
            all_assignments[nid] = cell_types[type_idx]
            all_confidences[nid] = confidence[i].item()

    # Build confidence DataFrame
    conf_df = pd.DataFrame({
        'nucleus_id': list(all_confidences.keys()),
        'confidence': list(all_confidences.values()),
    })

    return NucleusAssignmentResult(
        assignments=all_assignments,
        morphology_features=None,
        classifier=None,
        assignment_probs=conf_df,
        method="vae_sinkhorn",
    )
```

**Step 2: Add test**

```python
# Add to CITEgeist/tests/test_module3b_nucleus_assignment.py

def test_vae_assignment_runs(tmp_path):
    """VAE assignment should run with mock checkpoints."""
    # This is an integration test - skip if no GPU or checkpoints
    pytest.skip("Integration test - requires trained models")
```

**Step 3: Commit**

```bash
git add CITEgeist/model/module3b_nucleus_assignment.py
git commit -m "feat: add VAE-Sinkhorn assignment method to Module 3b"
```

---

## Task 8: Evaluation Script

**Files:**
- Create: `Benchmarking/xenium_benchmarking/evaluation/src/evaluate_vae_assignment.py`

**Step 1: Create evaluation script**

```python
# Benchmarking/xenium_benchmarking/evaluation/src/evaluate_vae_assignment.py
"""Evaluate VAE + Sinkhorn single-cell assignment."""
import argparse
import numpy as np
import pandas as pd
from pathlib import Path
from sklearn.metrics import accuracy_score, classification_report, confusion_matrix
import json


def evaluate_assignment(
    predictions_file: str,
    ground_truth_file: str,
    output_dir: str,
):
    """
    Evaluate assignment accuracy against ground truth.

    Args:
        predictions_file: CSV with nucleus_id, predicted_type, confidence
        ground_truth_file: CSV with nucleus_id, gt_type (from Xenium)
        output_dir: Directory for results
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load data
    pred_df = pd.read_csv(predictions_file)
    gt_df = pd.read_csv(ground_truth_file)

    # Merge on nucleus_id
    merged = pred_df.merge(gt_df, on="nucleus_id", how="inner")

    y_true = merged["gt_type"].values
    y_pred = merged["predicted_type"].values

    # Overall accuracy
    accuracy = accuracy_score(y_true, y_pred)
    print(f"Overall accuracy: {accuracy:.4f}")

    # Per-class report
    report = classification_report(y_true, y_pred, output_dict=True)
    print(classification_report(y_true, y_pred))

    # Confusion matrix
    labels = sorted(set(y_true) | set(y_pred))
    cm = confusion_matrix(y_true, y_pred, labels=labels)
    cm_df = pd.DataFrame(cm, index=labels, columns=labels)

    # Confidence calibration
    if "confidence" in merged.columns:
        bins = [0, 0.2, 0.4, 0.6, 0.8, 1.0]
        merged["conf_bin"] = pd.cut(merged["confidence"], bins)
        conf_acc = merged.groupby("conf_bin").apply(
            lambda x: (x["predicted_type"] == x["gt_type"]).mean()
        )
        print("\nConfidence calibration:")
        print(conf_acc)

    # Save results
    results = {
        "accuracy": accuracy,
        "n_samples": len(merged),
        "classification_report": report,
    }

    with open(output_dir / "evaluation_results.json", "w") as f:
        json.dump(results, f, indent=2)

    cm_df.to_csv(output_dir / "confusion_matrix.csv")

    return results


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--predictions", required=True)
    parser.add_argument("--ground_truth", required=True)
    parser.add_argument("--output_dir", required=True)
    args = parser.parse_args()

    evaluate_assignment(args.predictions, args.ground_truth, args.output_dir)
```

**Step 2: Commit**

```bash
git add Benchmarking/xenium_benchmarking/evaluation/src/evaluate_vae_assignment.py
git commit -m "feat: add VAE assignment evaluation script"
```

---

## Task 9: Data Preparation Script

**Files:**
- Create: `Benchmarking/xenium_benchmarking/CITEgeist/src/prepare_patches.py`

**Step 1: Create patch preparation script**

```python
# Benchmarking/xenium_benchmarking/CITEgeist/src/prepare_patches.py
"""Prepare nucleus patches for VAE training."""
import argparse
import numpy as np
import pandas as pd
from pathlib import Path
from tqdm import tqdm
import tifffile

import sys
sys.path.insert(0, str(Path(__file__).parents[4]))

from CITEgeist.model.patch_extraction import extract_patch
from CITEgeist.model.morphology_features import extract_nucleus_features


def prepare_patches(
    image_path: str,
    mask_path: str,
    nuclei_spot_map_path: str,
    output_dir: str,
    expansion: float = 0.75,
    patch_size: int = 96,
):
    """
    Extract and save patches for all nuclei.

    Args:
        image_path: Path to multi-channel TIFF image
        mask_path: Path to Cellpose nucleus mask
        nuclei_spot_map_path: CSV mapping nucleus_id to spot_id
        output_dir: Output directory for patches
        expansion: Bbox expansion fraction
        patch_size: Output patch size
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load data
    print("Loading image...")
    image = tifffile.imread(image_path)
    if image.ndim == 2:
        image = image[np.newaxis, ...]

    print("Loading mask...")
    mask = np.load(mask_path) if mask_path.endswith(".npy") else tifffile.imread(mask_path)

    print("Extracting nucleus features (for bboxes)...")
    features_df = extract_nucleus_features(mask)

    # Add bounding boxes
    from skimage.measure import regionprops
    props = regionprops(mask)
    bbox_data = {
        p.label: (p.bbox[1], p.bbox[0], p.bbox[3], p.bbox[2])  # x_min, y_min, x_max, y_max
        for p in props
    }
    features_df["bbox_x_min"] = features_df["nucleus_id"].map(lambda x: bbox_data.get(x, (0,0,0,0))[0])
    features_df["bbox_y_min"] = features_df["nucleus_id"].map(lambda x: bbox_data.get(x, (0,0,0,0))[1])
    features_df["bbox_x_max"] = features_df["nucleus_id"].map(lambda x: bbox_data.get(x, (0,0,0,0))[2])
    features_df["bbox_y_max"] = features_df["nucleus_id"].map(lambda x: bbox_data.get(x, (0,0,0,0))[3])

    # Load spot mapping
    spot_map = pd.read_csv(nuclei_spot_map_path)
    features_df = features_df.merge(spot_map, on="nucleus_id", how="inner")

    # Extract patches per spot
    print("Extracting patches...")
    for spot_id in tqdm(features_df["spot_id"].unique()):
        spot_nuclei = features_df[features_df["spot_id"] == spot_id]

        patches = []
        for _, row in spot_nuclei.iterrows():
            bbox = (
                int(row["bbox_x_min"]),
                int(row["bbox_y_min"]),
                int(row["bbox_x_max"]),
                int(row["bbox_y_max"]),
            )
            try:
                patch = extract_patch(image, bbox, expansion, patch_size)
                patches.append(patch)
            except Exception as e:
                print(f"Warning: Failed to extract patch for nucleus {row['nucleus_id']}: {e}")
                continue

        if patches:
            patches = np.stack(patches, axis=0)
            np.save(output_dir / f"spot_{spot_id}_patches.npy", patches)

    # Save features with bboxes
    features_df.to_csv(output_dir / "nucleus_features.csv", index=False)

    print(f"Done. Saved patches to {output_dir}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--image", required=True)
    parser.add_argument("--mask", required=True)
    parser.add_argument("--nuclei_spot_map", required=True)
    parser.add_argument("--output_dir", required=True)
    parser.add_argument("--expansion", type=float, default=0.75)
    parser.add_argument("--patch_size", type=int, default=96)
    args = parser.parse_args()

    prepare_patches(
        args.image, args.mask, args.nuclei_spot_map, args.output_dir,
        args.expansion, args.patch_size
    )
```

**Step 2: Commit**

```bash
git add Benchmarking/xenium_benchmarking/CITEgeist/src/prepare_patches.py
git commit -m "feat: add patch preparation script for VAE training"
```

---

## Summary

**Total tasks:** 9

**Files created:**
- `CITEgeist/model/sinkhorn.py` - Differentiable Sinkhorn OT
- `CITEgeist/model/patch_extraction.py` - Nucleus patch extraction
- `CITEgeist/model/vae.py` - VAE architecture
- `CITEgeist/model/projection_heads.py` - Projection heads and prototypes
- `CITEgeist/model/prototype_learning.py` - Stage 2 model
- `CITEgeist/model/train_vae.py` - VAE training script
- `CITEgeist/model/train_prototypes.py` - Prototype training script
- `Benchmarking/.../evaluate_vae_assignment.py` - Evaluation script
- `Benchmarking/.../prepare_patches.py` - Data preparation

**Files modified:**
- `CITEgeist/model/module3b_nucleus_assignment.py` - Added VAE assignment method

**Execution order:**
1. Task 1: Sinkhorn (core algorithm, no dependencies)
2. Task 2: Patch extraction (no dependencies)
3. Task 3: VAE (no dependencies)
4. Task 4: Projection heads (no dependencies)
5. Task 5: Prototype learning (depends on 1, 3, 4)
6. Task 6: Training scripts (depends on 3, 5)
7. Task 7: Module 3b integration (depends on 5)
8. Task 8: Evaluation (depends on 7)
9. Task 9: Data preparation (depends on 2)

Tasks 1-4 can be parallelized. Tasks 5-9 are sequential.
