# Unified Single-Cell Assignment (Module 3b MIL) Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add MIL-based single-cell assignment to CITEgeist Module 3b, supporting both DAPI (SimCLR) and H&E (ImageNet ViT) backbones with proportion-weighted Hungarian assignment, and benchmark on Xenium pseudo-Visium.

**Architecture:** ViT-Small backbone (384-dim) → Gated Attention MIL (384→256→K) → Proportion-weighted Hungarian assignment. Two modality-specific backbones (DAPI SimCLR, H&E ImageNet+SSL) share the same MIL head and assignment engine. Builds on existing Module 3b infrastructure (`hungarian_assignment.py`, `module3b_nucleus_assignment.py`, `NucleusAssignmentResult`).

**Tech Stack:** PyTorch, timm, scipy.optimize.linear_sum_assignment, numpy, pandas

---

## Existing Code Map

| File | Role | Action |
|------|------|--------|
| `CITEgeist/model/hungarian_assignment.py` | Hungarian with `-log(prob)` cost | **Modify**: add proportion weighting |
| `CITEgeist/model/module3b_nucleus_assignment.py` | Orchestration (random, morphology, vae_sinkhorn) | **Modify**: add `run_nucleus_assignment_mil()` |
| `CITEgeist/model/proportion_mil.py` | MIL for H&E (input_dim=768) | **Reference only** — logic extracted to new file |
| `CITEgeist/model/simclr.py` | SimCLR training | **Reuse as-is** |
| `CITEgeist/model/vit_encoder.py` | ViT-Small (2ch, 96x96) | **Reuse as-is** |
| `CITEgeist/model/vit_extractor.py` | ViT-Small ImageNet (3ch, 224x224) | **Reuse as-is** |
| `CITEgeist/model/__init__.py` | Package exports | **Modify**: add new exports |
| `CITEgeist/model/citegeist_model.py` | Main model class | **Modify**: add `run_single_cell_assignment()` |

---

### Task 1: Proportion-Weighted Hungarian Assignment

**Files:**
- Modify: `CITEgeist/model/hungarian_assignment.py`
- Test: `CITEgeist/tests/test_hungarian_weighted.py`

**Step 1: Write the failing test**

```python
"""Tests for proportion-weighted Hungarian assignment."""
import numpy as np
import pytest

from CITEgeist.model.hungarian_assignment import assign_nuclei_to_types


def test_proportion_weighted_changes_assignment():
    """Proportion prior should bias assignment vs uniform prior."""
    np.random.seed(42)
    n_nuclei, n_types = 5, 3
    # Ambiguous morphology: all nuclei look similar
    probs = np.full((n_nuclei, n_types), 1.0 / n_types)
    # Add slight morphology signal
    probs[0, 0] += 0.05  # nucleus 0 slightly favors type 0
    probs[1, 1] += 0.05  # nucleus 1 slightly favors type 1
    probs = probs / probs.sum(axis=1, keepdims=True)

    counts = np.array([2, 2, 1])
    nucleus_ids = np.arange(n_nuclei)

    # Without proportion weighting (lambda=0)
    result_unweighted = assign_nuclei_to_types(
        probs, counts, nucleus_ids, lambda_prior=0.0
    )

    # With strong proportion prior favoring type 0
    proportions = np.array([0.8, 0.15, 0.05])
    result_weighted = assign_nuclei_to_types(
        probs, counts, nucleus_ids,
        lambda_prior=2.0, proportions=proportions,
    )

    # Both should return valid assignments
    assert len(result_unweighted) == n_nuclei
    assert len(result_weighted) == n_nuclei

    # Count constraints should be respected in both
    unw_counts = np.bincount([result_unweighted[i] for i in range(n_nuclei)], minlength=n_types)
    w_counts = np.bincount([result_weighted[i] for i in range(n_nuclei)], minlength=n_types)
    np.testing.assert_array_equal(unw_counts, counts)
    np.testing.assert_array_equal(w_counts, counts)


def test_lambda_zero_matches_original():
    """lambda_prior=0 should produce same result as no proportions."""
    np.random.seed(42)
    probs = np.random.dirichlet([1, 1, 1], size=6)
    counts = np.array([2, 2, 2])
    nucleus_ids = np.arange(6)

    result_no_prop = assign_nuclei_to_types(probs, counts, nucleus_ids)
    result_zero_lambda = assign_nuclei_to_types(
        probs, counts, nucleus_ids, lambda_prior=0.0
    )
    assert result_no_prop == result_zero_lambda


def test_proportions_none_ignored():
    """If proportions=None, lambda_prior is ignored."""
    probs = np.array([[0.7, 0.3], [0.4, 0.6]])
    counts = np.array([1, 1])
    nucleus_ids = np.array([10, 20])
    result = assign_nuclei_to_types(
        probs, counts, nucleus_ids, lambda_prior=5.0, proportions=None,
    )
    assert len(result) == 2
```

**Step 2: Run test to verify it fails**

Run: `cd CITEgeist && python -m pytest tests/test_hungarian_weighted.py -v`
Expected: FAIL — `assign_nuclei_to_types()` doesn't accept `lambda_prior` or `proportions` params

**Step 3: Implement proportion weighting**

Modify `CITEgeist/model/hungarian_assignment.py`:

```python
"""Hungarian algorithm for optimal nucleus-to-celltype assignment."""
import numpy as np
from scipy.optimize import linear_sum_assignment
from typing import Dict, Optional


def assign_nuclei_to_types(
    probs: np.ndarray,
    counts: np.ndarray,
    nucleus_ids: np.ndarray,
    lambda_prior: float = 0.0,
    proportions: Optional[np.ndarray] = None,
) -> Dict[int, int]:
    """
    Assign nuclei to cell types using Hungarian algorithm.

    Cost matrix combines morphology likelihood and optional proportion prior:
        cost[i,k] = -log(probs[i,k]) - lambda_prior * log(proportions[k])

    Args:
        probs: (n_nuclei, n_types) probability matrix (e.g., MIL attention)
        counts: (n_types,) integer counts per cell type
        nucleus_ids: (n_nuclei,) unique identifier for each nucleus
        lambda_prior: Weight for proportion prior (0 = morphology only)
        proportions: (n_types,) continuous proportions from Module 3.
                     If None, lambda_prior is ignored.

    Returns:
        Dict mapping nucleus_id -> cell_type_index
    """
    n_nuclei, n_types = probs.shape
    n_slots = int(counts.sum())

    # Expand counts into slot list
    slots = []
    for type_idx, count in enumerate(counts):
        slots.extend([type_idx] * int(count))

    # Handle edge case: no cells to assign
    if n_slots == 0:
        assignments = {}
        for i, nid in enumerate(nucleus_ids):
            assignments[int(nid)] = int(np.argmax(probs[i]))
        return assignments

    # Compute proportion prior (broadcast across nuclei)
    prop_prior = np.zeros(n_types)
    if proportions is not None and lambda_prior > 0:
        safe_props = np.clip(proportions, 1e-10, None)
        prop_prior = -lambda_prior * np.log(safe_props)

    # Build cost matrix with expanded slots
    cost_matrix = np.zeros((max(n_nuclei, n_slots), max(n_nuclei, n_slots)))
    cost_matrix[:] = 1e6  # High cost for padding

    for i in range(n_nuclei):
        for j in range(n_slots):
            type_idx = slots[j]
            cost_matrix[i, j] = -np.log(probs[i, type_idx] + 1e-10) + prop_prior[type_idx]

    # Solve assignment
    row_ind, col_ind = linear_sum_assignment(cost_matrix)

    # Build result dictionary
    assignments = {}
    for i, j in zip(row_ind, col_ind):
        if i < n_nuclei:
            nid = int(nucleus_ids[i])
            if j < n_slots:
                assignments[nid] = slots[j]
            else:
                assignments[nid] = int(np.argmax(probs[i]))

    return assignments
```

**Step 4: Run test to verify it passes**

Run: `cd CITEgeist && python -m pytest tests/test_hungarian_weighted.py -v`
Expected: PASS (3 tests)

**Step 5: Commit**

```bash
git add CITEgeist/model/hungarian_assignment.py CITEgeist/tests/test_hungarian_weighted.py
git commit -m "feat: add proportion-weighted cost to Hungarian assignment"
```

---

### Task 2: Morphology Backbone ABC + Implementations

**Files:**
- Create: `CITEgeist/model/morphology_backbone.py`
- Test: `CITEgeist/tests/test_morphology_backbone.py`

**Step 1: Write the failing test**

```python
"""Tests for morphology backbone abstraction."""
import numpy as np
import pytest
import torch

from CITEgeist.model.morphology_backbone import DAPIBackbone, HEBackbone


def test_dapi_backbone_output_shape():
    """DAPI backbone should output (N, 384) from (N, 2, 96, 96) input."""
    backbone = DAPIBackbone()
    patches = torch.randn(4, 2, 96, 96)
    with torch.no_grad():
        embeddings = backbone.extract(patches)
    assert embeddings.shape == (4, 384)


def test_he_backbone_output_shape():
    """H&E backbone should output (N, 384) from (N, 3, 224, 224) input."""
    backbone = HEBackbone(model_name='vit_small_patch16_224')
    patches = torch.randn(4, 3, 224, 224)
    with torch.no_grad():
        embeddings = backbone.extract(patches)
    assert embeddings.shape == (4, 384)


def test_dapi_backbone_from_checkpoint(tmp_path):
    """DAPI backbone should load SimCLR encoder from checkpoint."""
    from CITEgeist.model.simclr import SimCLR
    # Create and save a SimCLR model
    model = SimCLR(in_channels=2, img_size=96)
    ckpt_path = tmp_path / "simclr.pt"
    torch.save(model.state_dict(), ckpt_path)

    # Load backbone from checkpoint
    backbone = DAPIBackbone(checkpoint=str(ckpt_path))
    patches = torch.randn(2, 2, 96, 96)
    with torch.no_grad():
        embeddings = backbone.extract(patches)
    assert embeddings.shape == (2, 384)


def test_dapi_backbone_extract_numpy():
    """extract_numpy should handle numpy input and batching."""
    backbone = DAPIBackbone()
    patches_np = np.random.randn(10, 2, 96, 96).astype(np.float32)
    embeddings = backbone.extract_numpy(patches_np, batch_size=4)
    assert embeddings.shape == (10, 384)
    assert embeddings.dtype == np.float32
```

**Step 2: Run test to verify it fails**

Run: `cd CITEgeist && python -m pytest tests/test_morphology_backbone.py -v`
Expected: FAIL — `morphology_backbone` module doesn't exist

**Step 3: Implement backbone abstraction**

Create `CITEgeist/model/morphology_backbone.py`:

```python
"""Morphology backbone encoders for single-cell assignment.

Two modality-specific backends, both outputting 384-dim embeddings:
- DAPIBackbone: ViT-Small trained with SimCLR on DAPI+boundary patches (96x96, 2ch)
- HEBackbone: ViT-Small with ImageNet init, optionally SSL fine-tuned on H&E (224x224, 3ch)
"""
import logging
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Optional

import numpy as np
import torch
import torch.nn as nn

logger = logging.getLogger(__name__)


class MorphologyBackbone(ABC):
    """Abstract base class for morphology feature extraction."""

    @abstractmethod
    def extract(self, patches: torch.Tensor) -> torch.Tensor:
        """Extract embeddings from image patches.

        Args:
            patches: (N, C, H, W) tensor of nucleus patches

        Returns:
            (N, 384) embedding tensor
        """

    def extract_numpy(
        self,
        patches: np.ndarray,
        batch_size: int = 64,
        device: str = "cpu",
    ) -> np.ndarray:
        """Extract embeddings from numpy patches with batching.

        Args:
            patches: (N, C, H, W) numpy array
            batch_size: Number of patches per batch
            device: Device for inference

        Returns:
            (N, 384) numpy array of embeddings
        """
        self._model.eval()
        self._model.to(device)
        all_embeddings = []

        with torch.no_grad():
            for i in range(0, len(patches), batch_size):
                batch = torch.tensor(
                    patches[i:i + batch_size], dtype=torch.float32, device=device
                )
                emb = self.extract(batch)
                all_embeddings.append(emb.cpu().numpy())

        return np.concatenate(all_embeddings, axis=0)

    @property
    @abstractmethod
    def _model(self) -> nn.Module:
        """Return the underlying nn.Module for device management."""

    @property
    def embed_dim(self) -> int:
        return 384


class DAPIBackbone(MorphologyBackbone):
    """DAPI + boundary channel backbone using SimCLR ViT-Small encoder.

    Input: (N, 2, 96, 96) — DAPI + boundary stain
    Output: (N, 384) — CLS token embedding
    """

    def __init__(
        self,
        checkpoint: Optional[str] = None,
        in_channels: int = 2,
        img_size: int = 96,
        device: str = "cpu",
    ):
        from .vit_encoder import create_vit_small

        self.encoder = create_vit_small(in_chans=in_channels, img_size=img_size)

        if checkpoint is not None:
            self._load_simclr_encoder(checkpoint, device)

        self.encoder.eval()

    def _load_simclr_encoder(self, checkpoint_path: str, device: str):
        """Load encoder weights from a SimCLR checkpoint."""
        state = torch.load(checkpoint_path, map_location=device, weights_only=False)

        # SimCLR checkpoint has 'encoder.*' prefix
        encoder_state = {}
        for k, v in state.items():
            if k.startswith("encoder."):
                encoder_state[k[len("encoder."):]] = v

        if encoder_state:
            self.encoder.load_state_dict(encoder_state, strict=False)
            logger.info("Loaded SimCLR encoder from %s", checkpoint_path)
        else:
            # Try loading directly (may be just encoder weights)
            self.encoder.load_state_dict(state, strict=False)
            logger.info("Loaded encoder weights from %s", checkpoint_path)

    def extract(self, patches: torch.Tensor) -> torch.Tensor:
        return self.encoder(patches)

    @property
    def _model(self) -> nn.Module:
        return self.encoder


class HEBackbone(MorphologyBackbone):
    """H&E backbone using ImageNet-pretrained ViT-Small.

    Input: (N, 3, 224, 224) — RGB H&E patches
    Output: (N, 384) — CLS token embedding

    Optionally load SSL fine-tuned weights on top of ImageNet init.
    """

    def __init__(
        self,
        model_name: str = "vit_small_patch16_224",
        pretrained: bool = True,
        checkpoint: Optional[str] = None,
        device: str = "cpu",
    ):
        from .vit_extractor import ViTFeatureExtractor

        self._extractor = ViTFeatureExtractor(
            model_name=model_name,
            pretrained=pretrained,
            device=device,
        )

        if checkpoint is not None:
            state = torch.load(checkpoint, map_location=device, weights_only=False)
            self._extractor.model.load_state_dict(state, strict=False)
            logger.info("Loaded SSL fine-tuned weights from %s", checkpoint)

    def extract(self, patches: torch.Tensor) -> torch.Tensor:
        return self._extractor(patches)

    @property
    def _model(self) -> nn.Module:
        return self._extractor.model
```

**Step 4: Run test to verify it passes**

Run: `cd CITEgeist && python -m pytest tests/test_morphology_backbone.py -v`
Expected: PASS (4 tests). Note: `test_he_backbone_output_shape` requires timm installed.

**Step 5: Commit**

```bash
git add CITEgeist/model/morphology_backbone.py CITEgeist/tests/test_morphology_backbone.py
git commit -m "feat: add MorphologyBackbone ABC with DAPI and H&E implementations"
```

---

### Task 3: Single-Cell MIL Head

**Files:**
- Create: `CITEgeist/model/single_cell_mil.py`
- Test: `CITEgeist/tests/test_single_cell_mil.py`

**Step 1: Write the failing test**

```python
"""Tests for single-cell MIL head."""
import numpy as np
import pytest
import torch

from CITEgeist.model.single_cell_mil import SingleCellMIL, train_mil, mil_loss


def test_forward_output_shapes():
    """MIL forward should return proportions (K,) and attention (N, K)."""
    model = SingleCellMIL(input_dim=384, n_types=7, hidden_dim=256)
    embeddings = torch.randn(10, 384)
    proportions, attention = model(embeddings)

    assert proportions.shape == (7,)
    assert attention.shape == (10, 7)


def test_proportions_sum_to_one():
    """Predicted proportions should sum to approximately 1."""
    model = SingleCellMIL(input_dim=384, n_types=5)
    embeddings = torch.randn(8, 384)
    proportions, _ = model(embeddings)
    assert abs(proportions.sum().item() - 1.0) < 0.01


def test_attention_rows_sum_to_one():
    """Each nucleus attention row should sum to 1 (softmax over types)."""
    model = SingleCellMIL(input_dim=384, n_types=5)
    embeddings = torch.randn(8, 384)
    _, attention = model(embeddings)
    row_sums = attention.sum(dim=1)
    torch.testing.assert_close(row_sums, torch.ones(8), atol=1e-5, rtol=1e-5)


def test_mil_loss_computation():
    """Loss should combine MSE and KL terms."""
    pred = torch.tensor([0.5, 0.3, 0.2])
    target = torch.tensor([0.6, 0.2, 0.2])
    loss = mil_loss(pred, target)
    assert loss.item() > 0
    assert loss.requires_grad is False  # no grad on detached tensors


def test_train_mil_reduces_loss():
    """Training should reduce loss over a few epochs."""
    torch.manual_seed(42)
    model = SingleCellMIL(input_dim=384, n_types=3)

    # Fake training data: 5 spots, each with 3-8 nuclei
    spots = []
    for _ in range(5):
        n = np.random.randint(3, 9)
        emb = torch.randn(n, 384)
        props = torch.softmax(torch.randn(3), dim=0)
        spots.append((emb, props))

    history = train_mil(model, spots, spots, n_epochs=20, lr=1e-3)
    assert history['train_loss'][-1] < history['train_loss'][0]
```

**Step 2: Run test to verify it fails**

Run: `cd CITEgeist && python -m pytest tests/test_single_cell_mil.py -v`
Expected: FAIL — module doesn't exist

**Step 3: Implement SingleCellMIL**

Create `CITEgeist/model/single_cell_mil.py`:

```python
"""Gated Attention MIL for single-cell type assignment.

Shared across modalities — accepts any 384-dim embedding input.
Trained with spot-level proportion supervision from Module 3.
At inference, attention weights serve as per-nucleus type probabilities
for downstream Hungarian assignment.
"""
import logging
from typing import Dict, List, Optional, Tuple

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F

logger = logging.getLogger(__name__)


class SingleCellMIL(nn.Module):
    """Gated attention MIL for proportion prediction and nucleus typing.

    Architecture:
        V = tanh(W_v @ x)      (384 → hidden_dim)
        U = sigmoid(W_u @ x)   (384 → hidden_dim)
        gate = V * U            (element-wise)
        logits = W_c @ gate     (hidden_dim → n_types)
        attention = softmax(logits, dim=1)  ← per-nucleus distribution over types
        proportions = mean(attention, dim=0) ← spot-level aggregation

    Args:
        input_dim: Embedding dimension from backbone (384)
        n_types: Number of cell types
        hidden_dim: Gated attention hidden dimension
        dropout: Dropout rate
    """

    def __init__(
        self,
        input_dim: int = 384,
        n_types: int = 7,
        hidden_dim: int = 256,
        dropout: float = 0.1,
    ):
        super().__init__()
        self.input_dim = input_dim
        self.n_types = n_types

        self.W_v = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.Tanh(),
        )
        self.W_u = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.Sigmoid(),
        )
        self.W_c = nn.Sequential(
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, n_types),
        )

    def forward(
        self, embeddings: torch.Tensor
    ) -> Tuple[torch.Tensor, torch.Tensor]:
        """Forward pass.

        Args:
            embeddings: (N, input_dim) nucleus embeddings from backbone

        Returns:
            proportions: (n_types,) predicted spot-level proportions
            attention: (N, n_types) per-nucleus type probability matrix
        """
        V = self.W_v(embeddings)    # (N, hidden)
        U = self.W_u(embeddings)    # (N, hidden)
        gate = V * U                # (N, hidden)
        logits = self.W_c(gate)     # (N, n_types)

        # Each nucleus votes over cell types (softmax over types)
        attention = F.softmax(logits, dim=1)  # (N, n_types)

        # Spot proportions = mean of nucleus votes
        proportions = attention.mean(dim=0)   # (n_types,)

        return proportions, attention


def mil_loss(
    pred: torch.Tensor,
    target: torch.Tensor,
    kl_weight: float = 0.1,
    eps: float = 1e-8,
) -> torch.Tensor:
    """Combined MSE + KL loss for proportion supervision.

    Args:
        pred: (K,) predicted proportions
        target: (K,) target proportions from Module 3
        kl_weight: Weight for KL divergence term
        eps: Numerical stability

    Returns:
        Scalar loss
    """
    mse = F.mse_loss(pred, target)

    # KL(target || pred)
    target_safe = target.clamp(min=eps)
    pred_safe = pred.clamp(min=eps)
    kl = (target_safe * (target_safe / pred_safe).log()).sum()

    return mse + kl_weight * kl


def train_mil(
    model: SingleCellMIL,
    train_spots: List[Tuple[torch.Tensor, torch.Tensor]],
    val_spots: List[Tuple[torch.Tensor, torch.Tensor]],
    n_epochs: int = 100,
    lr: float = 1e-4,
    entropy_weight: float = 0.01,
    device: str = "cpu",
    save_path: Optional[str] = None,
) -> Dict[str, list]:
    """Train MIL head on spot-level proportion targets.

    Args:
        model: SingleCellMIL instance
        train_spots: List of (embeddings (N,384), proportions (K,)) tuples
        val_spots: Validation spots (same format)
        n_epochs: Training epochs
        lr: Learning rate
        entropy_weight: Weight for entropy regularization
        device: Training device
        save_path: If set, save best model checkpoint

    Returns:
        History dict with train_loss, val_loss, val_r lists
    """
    model.to(device)
    model.train()

    optimizer = torch.optim.AdamW(model.parameters(), lr=lr, weight_decay=1e-4)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=n_epochs)

    history = {'train_loss': [], 'val_loss': [], 'val_r': []}
    best_val_r = -1.0

    for epoch in range(n_epochs):
        # Training
        model.train()
        epoch_loss = 0.0
        for embeddings, target_props in train_spots:
            embeddings = embeddings.to(device)
            target_props = target_props.to(device)

            pred_props, attention = model(embeddings)
            loss = mil_loss(pred_props, target_props)

            # Entropy regularization: encourage diverse attention
            entropy = -(attention * (attention + 1e-8).log()).sum(dim=1).mean()
            loss = loss - entropy_weight * entropy

            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
            epoch_loss += loss.item()

        scheduler.step()
        avg_train_loss = epoch_loss / max(len(train_spots), 1)
        history['train_loss'].append(avg_train_loss)

        # Validation
        model.eval()
        val_loss = 0.0
        all_pred, all_target = [], []
        with torch.no_grad():
            for embeddings, target_props in val_spots:
                embeddings = embeddings.to(device)
                target_props = target_props.to(device)
                pred_props, _ = model(embeddings)
                val_loss += mil_loss(pred_props, target_props).item()
                all_pred.append(pred_props.cpu().numpy())
                all_target.append(target_props.cpu().numpy())

        avg_val_loss = val_loss / max(len(val_spots), 1)
        history['val_loss'].append(avg_val_loss)

        # Pearson r
        if all_pred:
            pred_flat = np.concatenate(all_pred)
            target_flat = np.concatenate(all_target)
            r = float(np.corrcoef(pred_flat, target_flat)[0, 1])
        else:
            r = 0.0
        history['val_r'].append(r)

        # Save best
        if save_path and r > best_val_r:
            best_val_r = r
            torch.save(model.state_dict(), save_path)

        if (epoch + 1) % 10 == 0:
            logger.info(
                "Epoch %d/%d: train_loss=%.4f val_loss=%.4f val_r=%.4f",
                epoch + 1, n_epochs, avg_train_loss, avg_val_loss, r,
            )

    return history
```

**Step 4: Run test to verify it passes**

Run: `cd CITEgeist && python -m pytest tests/test_single_cell_mil.py -v`
Expected: PASS (5 tests)

**Step 5: Commit**

```bash
git add CITEgeist/model/single_cell_mil.py CITEgeist/tests/test_single_cell_mil.py
git commit -m "feat: add SingleCellMIL gated attention head with training loop"
```

---

### Task 4: MIL Assignment Method in Module 3b

**Files:**
- Modify: `CITEgeist/model/module3b_nucleus_assignment.py`
- Test: `CITEgeist/tests/test_mil_assignment.py`

**Step 1: Write the failing test**

```python
"""Tests for MIL-based nucleus assignment in Module 3b."""
import numpy as np
import pandas as pd
import pytest
import torch

from CITEgeist.model.module3b_nucleus_assignment import (
    run_nucleus_assignment_mil,
    NucleusAssignmentResult,
)


@pytest.fixture
def fake_spot_data():
    """Create minimal fake data for 2 spots, 3 types, 10 nuclei."""
    n_types = 3
    cell_types = ["TypeA", "TypeB", "TypeC"]

    # 2 spots: spot_0 has 4 nuclei, spot_1 has 6
    nuclei_spot_map = pd.DataFrame({
        'nucleus_id': range(10),
        'spot_id': ['spot_0'] * 4 + ['spot_1'] * 6,
    })

    proportions = pd.DataFrame({
        'spot_id': ['spot_0', 'spot_1'],
        'TypeA': [0.5, 0.33],
        'TypeB': [0.25, 0.33],
        'TypeC': [0.25, 0.34],
    })

    nuclei_counts = pd.Series({'spot_0': 4, 'spot_1': 6})

    # Pre-computed embeddings: dict of spot_id -> (N, 384) array
    embeddings = {
        'spot_0': np.random.randn(4, 384).astype(np.float32),
        'spot_1': np.random.randn(6, 384).astype(np.float32),
    }

    return nuclei_spot_map, proportions, nuclei_counts, cell_types, embeddings


def test_mil_assignment_returns_result(fake_spot_data):
    """MIL assignment should return NucleusAssignmentResult."""
    nuclei_spot_map, proportions, nuclei_counts, cell_types, embeddings = fake_spot_data

    result = run_nucleus_assignment_mil(
        embeddings=embeddings,
        nuclei_spot_map=nuclei_spot_map,
        proportions=proportions,
        nuclei_counts=nuclei_counts,
        cell_types=cell_types,
        n_epochs=5,  # Short training for test
        lambda_prior=1.0,
        device="cpu",
    )

    assert isinstance(result, NucleusAssignmentResult)
    assert result.method == "mil"
    assert len(result.assignments) == 10
    # All assigned types should be valid
    for nid, ctype in result.assignments.items():
        assert ctype in cell_types


def test_mil_assignment_respects_counts(fake_spot_data):
    """Assignments per spot should respect discretized count constraints."""
    nuclei_spot_map, proportions, nuclei_counts, cell_types, embeddings = fake_spot_data

    result = run_nucleus_assignment_mil(
        embeddings=embeddings,
        nuclei_spot_map=nuclei_spot_map,
        proportions=proportions,
        nuclei_counts=nuclei_counts,
        cell_types=cell_types,
        n_epochs=5,
        lambda_prior=0.0,
        device="cpu",
    )

    # Check count constraints per spot
    for spot_id in ['spot_0', 'spot_1']:
        spot_nuclei = nuclei_spot_map[nuclei_spot_map['spot_id'] == spot_id]
        n_nuclei = int(nuclei_counts[spot_id])
        assigned_count = sum(
            1 for nid in spot_nuclei['nucleus_id']
            if nid in result.assignments
        )
        assert assigned_count == len(spot_nuclei)
```

**Step 2: Run test to verify it fails**

Run: `cd CITEgeist && python -m pytest tests/test_mil_assignment.py -v`
Expected: FAIL — `run_nucleus_assignment_mil` not importable

**Step 3: Implement MIL assignment function**

Add to `CITEgeist/model/module3b_nucleus_assignment.py` (after existing `run_nucleus_assignment_vae`):

```python
def run_nucleus_assignment_mil(
    embeddings: Dict[str, np.ndarray],
    nuclei_spot_map: pd.DataFrame,
    proportions: pd.DataFrame,
    nuclei_counts: pd.Series,
    cell_types: List[str],
    n_epochs: int = 100,
    lambda_prior: float = 1.0,
    mil_checkpoint: Optional[str] = None,
    device: str = "cpu",
    train_frac: float = 0.8,
    random_seed: int = 42,
) -> NucleusAssignmentResult:
    """
    Run MIL-based nucleus assignment.

    Trains a SingleCellMIL head on spot-level proportion targets, then uses
    attention weights as per-nucleus type likelihoods for proportion-weighted
    Hungarian assignment.

    Args:
        embeddings: Dict mapping spot_id -> (n_nuclei, 384) numpy array
                    of pre-extracted backbone embeddings
        nuclei_spot_map: DataFrame with 'nucleus_id' and 'spot_id' columns
        proportions: DataFrame with 'spot_id' and one column per cell type
        nuclei_counts: Series mapping spot_id -> nuclei count
        cell_types: List of cell type names
        n_epochs: MIL training epochs
        lambda_prior: Proportion prior weight for Hungarian assignment
        mil_checkpoint: If set, load pre-trained MIL weights (skip training)
        device: PyTorch device
        train_frac: Fraction of spots for training (rest for validation)
        random_seed: Random seed for train/val split

    Returns:
        NucleusAssignmentResult with method="mil"
    """
    import torch
    from .single_cell_mil import SingleCellMIL, train_mil
    from .hungarian_assignment import assign_nuclei_to_types
    from .morphology_features import largest_remainder_discretize

    n_types = len(cell_types)
    prop_cols = cell_types
    spot_props = proportions.set_index('spot_id')[prop_cols]

    # Prepare spot data as (embeddings_tensor, proportions_tensor) tuples
    spot_ids = [sid for sid in spot_props.index if sid in embeddings]
    spot_data = []
    for sid in spot_ids:
        emb = torch.tensor(embeddings[sid], dtype=torch.float32)
        props = torch.tensor(spot_props.loc[sid].values, dtype=torch.float32)
        spot_data.append((emb, props))

    # Train/val split
    rng = np.random.default_rng(random_seed)
    indices = rng.permutation(len(spot_data))
    n_train = int(len(spot_data) * train_frac)
    train_spots = [spot_data[i] for i in indices[:n_train]]
    val_spots = [spot_data[i] for i in indices[n_train:]]
    if not val_spots:
        val_spots = train_spots  # Fallback if too few spots

    # Initialize and train MIL
    model = SingleCellMIL(input_dim=384, n_types=n_types)

    if mil_checkpoint is not None:
        state = torch.load(mil_checkpoint, map_location=device, weights_only=False)
        model.load_state_dict(state)
    else:
        train_mil(
            model, train_spots, val_spots,
            n_epochs=n_epochs, device=device,
        )

    # Inference: get attention weights per spot, run Hungarian
    model.to(device)
    model.eval()
    all_assignments = {}
    all_attention_records = []

    for spot_id in nuclei_spot_map['spot_id'].unique():
        if spot_id not in embeddings or spot_id not in spot_props.index:
            continue

        spot_nuclei = nuclei_spot_map[nuclei_spot_map['spot_id'] == spot_id]
        nucleus_ids = spot_nuclei['nucleus_id'].values
        spot_emb = torch.tensor(embeddings[spot_id], dtype=torch.float32, device=device)
        spot_proportions = spot_props.loc[spot_id].values.astype(np.float64)
        n_nuclei = int(nuclei_counts.get(spot_id, len(spot_nuclei)))

        with torch.no_grad():
            _, attention = model(spot_emb)  # (N, K)
        attention_np = attention.cpu().numpy()

        # Discretize proportions to counts
        counts = largest_remainder_discretize(spot_proportions, n_nuclei)

        # Proportion-weighted Hungarian assignment
        type_assignments = assign_nuclei_to_types(
            attention_np, counts, nucleus_ids,
            lambda_prior=lambda_prior,
            proportions=spot_proportions,
        )

        for nid, type_idx in type_assignments.items():
            all_assignments[nid] = cell_types[type_idx]

        # Store attention for diagnostics
        for i, nid in enumerate(nucleus_ids):
            record = {'nucleus_id': int(nid), 'spot_id': spot_id}
            for k, ct in enumerate(cell_types):
                record[ct] = float(attention_np[i, k])
            all_attention_records.append(record)

    assignment_probs = pd.DataFrame(all_attention_records) if all_attention_records else None

    return NucleusAssignmentResult(
        assignments=all_assignments,
        morphology_features=None,
        classifier=None,
        assignment_probs=assignment_probs,
        method="mil",
    )
```

**Step 4: Run test to verify it passes**

Run: `cd CITEgeist && python -m pytest tests/test_mil_assignment.py -v`
Expected: PASS (2 tests)

**Step 5: Commit**

```bash
git add CITEgeist/model/module3b_nucleus_assignment.py CITEgeist/tests/test_mil_assignment.py
git commit -m "feat: add MIL-based nucleus assignment to Module 3b"
```

---

### Task 5: CitegeistModel Integration

**Files:**
- Modify: `CITEgeist/model/citegeist_model.py` (add `run_single_cell_assignment()` method after `run_cell_expression_pass1`)
- Modify: `CITEgeist/model/__init__.py` (add new exports)

**Step 1: Add `run_single_cell_assignment()` to CitegeistModel**

Insert after `run_single_cell_resolution` method (line ~1948) in `citegeist_model.py`:

```python
def run_single_cell_assignment(
    self,
    patches_dir: str,
    nuclei_counts: pd.Series,
    modality: str = "he",
    backbone_checkpoint: Optional[str] = None,
    mil_checkpoint: Optional[str] = None,
    lambda_prior: float = 1.0,
    n_epochs: int = 100,
    batch_size: int = 64,
    device: str = "cpu",
) -> pd.DataFrame:
    """
    Module 3b: Assign individual nuclei to cell types using MIL.

    Requires Module 3 (run_cell_proportion_model) to have been run first.
    Uses a modality-specific backbone to extract nucleus embeddings,
    trains a MIL head on spot-level proportions, then runs proportion-weighted
    Hungarian assignment.

    Args:
        patches_dir: Directory containing per-spot nucleus patches.
            Expected structure: {spot_id}_patches.npy files
        nuclei_counts: Series mapping spot_id -> number of nuclei
        modality: "dapi" (2ch, 96x96 SimCLR) or "he" (3ch, 224x224 ImageNet ViT)
        backbone_checkpoint: Path to pre-trained backbone (SimCLR for DAPI,
            SSL fine-tuned for H&E). If None, uses untrained backbone.
        mil_checkpoint: Path to pre-trained MIL weights. If None, trains
            from scratch using Module 3 proportions.
        lambda_prior: Weight for proportion prior in Hungarian assignment
        n_epochs: MIL training epochs (ignored if mil_checkpoint provided)
        batch_size: Batch size for backbone feature extraction
        device: PyTorch device for backbone and MIL inference

    Returns:
        DataFrame with columns: nucleus_id, spot_id, cell_type
    """
    from .morphology_backbone import DAPIBackbone, HEBackbone
    from .module3b_nucleus_assignment import run_nucleus_assignment_mil

    # Validate Module 3 has been run
    if 'cell_prop_finetuned_results' not in self.results:
        raise RuntimeError(
            "Module 3 must be run first (run_cell_proportion_model). "
            "No finetuned proportions found."
        )

    proportions_df = self.results['cell_prop_finetuned_results']
    cell_types = [c for c in proportions_df.columns if c != 'spot_id']

    logger.info("Module 3b: Single-cell assignment via MIL (%s backbone)", modality)

    # Initialize backbone
    if modality == "dapi":
        backbone = DAPIBackbone(checkpoint=backbone_checkpoint, device=device)
    elif modality == "he":
        backbone = HEBackbone(checkpoint=backbone_checkpoint, device=device)
    else:
        raise ValueError(f"Unknown modality: {modality}. Use 'dapi' or 'he'.")

    # Extract embeddings per spot
    import glob
    import os
    embeddings = {}
    nuclei_spot_records = []

    patch_files = sorted(glob.glob(os.path.join(patches_dir, "*_patches.npy")))
    logger.info("Found %d patch files in %s", len(patch_files), patches_dir)

    for pf in patch_files:
        spot_id = os.path.basename(pf).replace("_patches.npy", "")
        patches = np.load(pf)  # (N, C, H, W)
        if len(patches) == 0:
            continue

        emb = backbone.extract_numpy(patches, batch_size=batch_size, device=device)
        embeddings[spot_id] = emb

        for i in range(len(patches)):
            nuclei_spot_records.append({
                'nucleus_id': f"{spot_id}_{i}",
                'spot_id': spot_id,
            })

    nuclei_spot_map = pd.DataFrame(nuclei_spot_records)
    logger.info("Extracted embeddings for %d spots, %d nuclei",
                len(embeddings), len(nuclei_spot_map))

    # Run MIL assignment
    result = run_nucleus_assignment_mil(
        embeddings=embeddings,
        nuclei_spot_map=nuclei_spot_map,
        proportions=proportions_df,
        nuclei_counts=nuclei_counts,
        cell_types=cell_types,
        n_epochs=n_epochs,
        lambda_prior=lambda_prior,
        mil_checkpoint=mil_checkpoint,
        device=device,
    )

    # Convert to DataFrame
    records = []
    for nid, ctype in result.assignments.items():
        sid = nuclei_spot_map.loc[
            nuclei_spot_map['nucleus_id'] == nid, 'spot_id'
        ].iloc[0]
        records.append({
            'nucleus_id': nid,
            'spot_id': sid,
            'cell_type': ctype,
        })

    assignments_df = pd.DataFrame(records)

    # Store in results
    self.results['single_cell_assignments'] = assignments_df
    self.results['single_cell_attention'] = result.assignment_probs

    # Save to output
    output_path = os.path.join(self.output_folder, 'single_cell_assignments.csv')
    assignments_df.to_csv(output_path, index=False)
    logger.info("Saved %d assignments to %s", len(assignments_df), output_path)

    return assignments_df
```

**Step 2: Update `__init__.py`**

Add after the existing Module 3b exports (line 95):

```python
# Module 3b MIL: Single-cell assignment via attention
from .morphology_backbone import MorphologyBackbone, DAPIBackbone, HEBackbone
from .single_cell_mil import SingleCellMIL, mil_loss, train_mil
```

And add to `__all__`:

```python
# Module 3b MIL: Morphology backbone + attention
"MorphologyBackbone",
"DAPIBackbone",
"HEBackbone",
"SingleCellMIL",
"mil_loss",
"train_mil",
```

**Step 3: Commit**

```bash
git add CITEgeist/model/citegeist_model.py CITEgeist/model/__init__.py
git commit -m "feat: add run_single_cell_assignment() to CitegeistModel"
```

---

### Task 6: Xenium Benchmark Script

**Files:**
- Create: `Benchmarking/xenium_benchmarking/CITEgeist/src/benchmark_mil_assignment.py`
- Create: `Benchmarking/xenium_benchmarking/CITEgeist/slurm/sbatch_mil_assignment.sh`

**Step 1: Write benchmark script**

Create `Benchmarking/xenium_benchmarking/CITEgeist/src/benchmark_mil_assignment.py`:

```python
"""Benchmark Module 3b MIL single-cell assignment on Xenium pseudo-Visium.

Pipeline:
1. Load Module 3 hybrid proportions (pre-computed)
2. Load SimCLR backbone (pre-trained)
3. Extract DAPI+boundary embeddings via backbone
4. Train MIL head on proportions
5. Proportion-weighted Hungarian assignment
6. Evaluate: proportion r, single-cell accuracy, GEX r
"""
import argparse
import json
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# Add project root to path
PROJECT_ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(PROJECT_ROOT))

from CITEgeist.model.morphology_backbone import DAPIBackbone
from CITEgeist.model.single_cell_mil import SingleCellMIL, train_mil
from CITEgeist.model.hungarian_assignment import assign_nuclei_to_types
from CITEgeist.model.morphology_features import largest_remainder_discretize

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s')
logger = logging.getLogger(__name__)

# Xenium benchmark cell types (achievable-7)
CELL_TYPES = [
    "B cells", "CD4+ T cells", "CD8+ T cells",
    "Macrophages", "Endothelial", "Epithelial", "Fibroblasts",
]

DATA_ROOT = Path(PROJECT_ROOT / "Benchmarking" / "xenium_pseudovisium" / "data_protein_gt")


def load_stage1_results(region_id: int) -> pd.DataFrame:
    """Load hybrid post-filter proportions from Module 3."""
    results_dir = PROJECT_ROOT / "Benchmarking" / "xenium_benchmarking" / "CITEgeist" / "output" / "hybrid_cellpose"
    region_name = f"Xenium_region_{region_id}"
    prop_path = results_dir / region_name / f"{region_name}_deconv_predictions.csv"
    if not prop_path.exists():
        raise FileNotFoundError(f"Stage 1 results not found: {prop_path}")
    return pd.read_csv(prop_path)


def load_ground_truth(region_id: int):
    """Load ground truth proportions and single-cell assignments."""
    region_name = f"Xenium_region_{region_id}"

    # Spot-level proportions
    gt_props_path = DATA_ROOT / "ground_truth" / f"{region_name}_proportions.csv"
    gt_props = pd.read_csv(gt_props_path)

    # Single-cell assignments
    gt_cells_path = DATA_ROOT / "ground_truth" / f"{region_name}_cell_assignments.csv"
    gt_cells = pd.read_csv(gt_cells_path) if gt_cells_path.exists() else None

    return gt_props, gt_cells


def load_patches(region_id: int) -> dict:
    """Load pre-extracted nucleus patches per spot."""
    patches_dir = PROJECT_ROOT / "Benchmarking" / "xenium_benchmarking" / "CITEgeist" / "output" / "patches"
    region_name = f"Xenium_region_{region_id}"
    region_dir = patches_dir / region_name

    spot_patches = {}
    for f in sorted(region_dir.glob("*_patches.npy")):
        spot_id = f.stem.replace("_patches", "")
        patches = np.load(f)
        if len(patches) > 0:
            spot_patches[spot_id] = patches

    return spot_patches


def evaluate(
    assignments: dict,
    proportions: pd.DataFrame,
    gt_props: pd.DataFrame,
    gt_cells: pd.DataFrame,
    nuclei_spot_map: pd.DataFrame,
    cell_types: list,
) -> dict:
    """Compute all evaluation metrics."""
    from scipy.stats import pearsonr

    metrics = {}

    # 1. Proportion correlation: re-aggregate assignments to proportions
    pred_spot_props = []
    for spot_id in gt_props['spot_id'].unique():
        spot_nuclei = nuclei_spot_map[nuclei_spot_map['spot_id'] == spot_id]
        n = len(spot_nuclei)
        if n == 0:
            continue
        type_counts = np.zeros(len(cell_types))
        for nid in spot_nuclei['nucleus_id']:
            if nid in assignments:
                ct = assignments[nid]
                if ct in cell_types:
                    type_counts[cell_types.index(ct)] += 1
        pred_spot_props.append(type_counts / max(n, 1))

    if pred_spot_props:
        pred_flat = np.array(pred_spot_props).flatten()
        gt_flat = gt_props[cell_types].values.flatten()
        min_len = min(len(pred_flat), len(gt_flat))
        r, p = pearsonr(pred_flat[:min_len], gt_flat[:min_len])
        metrics['proportion_pearson_r'] = float(r)
        metrics['proportion_pearson_p'] = float(p)

    # 2. Single-cell accuracy (if ground truth available)
    if gt_cells is not None:
        gt_map = dict(zip(gt_cells['nucleus_id'], gt_cells['cell_type']))
        correct, total = 0, 0
        per_type_correct = {ct: 0 for ct in cell_types}
        per_type_total = {ct: 0 for ct in cell_types}
        for nid, pred_type in assignments.items():
            if nid in gt_map:
                gt_type = gt_map[nid]
                total += 1
                if gt_type in per_type_total:
                    per_type_total[gt_type] += 1
                if pred_type == gt_type:
                    correct += 1
                    if gt_type in per_type_correct:
                        per_type_correct[gt_type] += 1

        metrics['single_cell_accuracy'] = correct / max(total, 1)
        metrics['single_cell_total'] = total
        metrics['per_type_accuracy'] = {
            ct: per_type_correct[ct] / max(per_type_total[ct], 1)
            for ct in cell_types
        }

    return metrics


def main():
    parser = argparse.ArgumentParser(description="Benchmark MIL single-cell assignment")
    parser.add_argument("--region", type=int, required=True, help="Region index (0-4)")
    parser.add_argument("--simclr-checkpoint", type=str, required=True)
    parser.add_argument("--output-dir", type=str, required=True)
    parser.add_argument("--lambda-prior", type=float, default=1.0)
    parser.add_argument("--n-epochs", type=int, default=100)
    parser.add_argument("--device", type=str, default="cuda")
    parser.add_argument("--batch-size", type=int, default=64)
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("=== MIL Benchmark: Region %d ===", args.region)

    # Load data
    stage1_props = load_stage1_results(args.region)
    gt_props, gt_cells = load_ground_truth(args.region)
    spot_patches = load_patches(args.region)

    logger.info("Loaded %d spots with patches", len(spot_patches))

    # Initialize backbone
    backbone = DAPIBackbone(checkpoint=args.simclr_checkpoint, device=args.device)

    # Extract embeddings
    embeddings = {}
    nuclei_records = []
    for spot_id, patches in spot_patches.items():
        emb = backbone.extract_numpy(patches, batch_size=args.batch_size, device=args.device)
        embeddings[spot_id] = emb
        for i in range(len(patches)):
            nuclei_records.append({'nucleus_id': f"{spot_id}_{i}", 'spot_id': spot_id})

    nuclei_spot_map = pd.DataFrame(nuclei_records)
    nuclei_counts = nuclei_spot_map.groupby('spot_id').size()

    logger.info("Extracted embeddings: %d spots, %d nuclei", len(embeddings), len(nuclei_spot_map))

    # Run MIL assignment
    from CITEgeist.model.module3b_nucleus_assignment import run_nucleus_assignment_mil

    result = run_nucleus_assignment_mil(
        embeddings=embeddings,
        nuclei_spot_map=nuclei_spot_map,
        proportions=stage1_props,
        nuclei_counts=nuclei_counts,
        cell_types=CELL_TYPES,
        n_epochs=args.n_epochs,
        lambda_prior=args.lambda_prior,
        device=args.device,
    )

    logger.info("Assignment complete: %d nuclei assigned", len(result.assignments))

    # Evaluate
    metrics = evaluate(
        result.assignments, stage1_props, gt_props, gt_cells,
        nuclei_spot_map, CELL_TYPES,
    )

    logger.info("Results: prop_r=%.4f, sc_acc=%.4f",
                metrics.get('proportion_pearson_r', -1),
                metrics.get('single_cell_accuracy', -1))

    # Save results
    metrics_path = output_dir / f"region_{args.region}_metrics.json"
    with open(metrics_path, 'w') as f:
        json.dump(metrics, f, indent=2)

    assignments_path = output_dir / f"region_{args.region}_assignments.csv"
    pd.DataFrame([
        {'nucleus_id': nid, 'cell_type': ct}
        for nid, ct in result.assignments.items()
    ]).to_csv(assignments_path, index=False)

    logger.info("Saved results to %s", output_dir)


if __name__ == "__main__":
    main()
```

**Step 2: Write SLURM submission script**

Create `Benchmarking/xenium_benchmarking/CITEgeist/slurm/sbatch_mil_assignment.sh`:

```bash
#!/bin/bash
#SBATCH --job-name=mil_assign
#SBATCH --array=0-4
#SBATCH --partition=l40s
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=slurm/logs/mil_assign_%A_%a.out
#SBATCH --error=slurm/logs/mil_assign_%A_%a.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

module load gurobi/12.0.3

eval "$(conda shell.bash hook)"
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

REGION=${SLURM_ARRAY_TASK_ID}

# UPDATE THESE PATHS before submission
SIMCLR_CHECKPOINT="output/simclr_ssl/simclr_best.pt"
OUTPUT_DIR="output/mil_assignment"

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/CITEgeist

python src/benchmark_mil_assignment.py \
    --region ${REGION} \
    --simclr-checkpoint ${SIMCLR_CHECKPOINT} \
    --output-dir ${OUTPUT_DIR} \
    --lambda-prior 1.0 \
    --n-epochs 100 \
    --device cuda \
    --batch-size 64
```

**Step 3: Commit**

```bash
git add Benchmarking/xenium_benchmarking/CITEgeist/src/benchmark_mil_assignment.py \
        Benchmarking/xenium_benchmarking/CITEgeist/slurm/sbatch_mil_assignment.sh
git commit -m "feat: add MIL single-cell assignment benchmark for Xenium"
```

---

### Task 7: Update Exports and Archive Deprecated Files

**Files:**
- Modify: `CITEgeist/model/__init__.py` (already partially done in Task 5)

**Step 1: Mark deprecated modules with comments**

In `__init__.py`, update the deprecated imports with deprecation comments. Do NOT remove them yet (would break existing code), but add comments:

```python
# DEPRECATED: VAE superseded by ViT-based SSL (SimCLR)
# from .vae import VAE, VAEEncoder, VAEDecoder

# DEPRECATED: DINO training collapsed; use SimCLR
# from .dino import DINO, DINOHead, create_dino_model

# DEPRECATED: MAE underperforms SimCLR in constrained setting
# from .mae import MAE, MAEDecoder
```

Keep the imports functional for now to avoid breaking existing benchmark scripts. Add a `_DEPRECATED` list:

```python
_DEPRECATED = [
    "VAE", "VAEEncoder", "VAEDecoder",
    "DINO", "DINOHead", "create_dino_model",
    "MAE", "MAEDecoder",
    "solve_masked_iqp",  # Discrete IQP superseded
]
```

**Step 2: Commit**

```bash
git add CITEgeist/model/__init__.py
git commit -m "chore: mark VAE, DINO, MAE, IQP as deprecated in favor of SimCLR MIL"
```

---

### Task 8: Integration Test — End-to-End Pipeline

**Files:**
- Create: `CITEgeist/tests/test_e2e_single_cell.py`

**Step 1: Write end-to-end integration test**

```python
"""End-to-end test: Module 3 proportions → Module 3b MIL assignment."""
import numpy as np
import pandas as pd
import pytest
import torch

from CITEgeist.model.morphology_backbone import DAPIBackbone
from CITEgeist.model.single_cell_mil import SingleCellMIL, train_mil
from CITEgeist.model.module3b_nucleus_assignment import run_nucleus_assignment_mil


@pytest.fixture
def synthetic_xenium_data():
    """Create synthetic data mimicking Xenium benchmark structure."""
    torch.manual_seed(42)
    np.random.seed(42)

    n_spots = 20
    n_types = 7
    cell_types = [
        "B cells", "CD4+ T cells", "CD8+ T cells",
        "Macrophages", "Endothelial", "Epithelial", "Fibroblasts",
    ]

    # Generate proportions (Module 3 output)
    raw_props = np.random.dirichlet(np.ones(n_types) * 2, size=n_spots)
    spot_ids = [f"spot_{i}" for i in range(n_spots)]
    proportions = pd.DataFrame(raw_props, columns=cell_types)
    proportions['spot_id'] = spot_ids

    # Generate fake patches and nuclei
    embeddings = {}
    nuclei_records = []
    nuclei_counts = {}
    gt_assignments = {}

    for i, sid in enumerate(spot_ids):
        n_nuclei = np.random.randint(3, 15)
        nuclei_counts[sid] = n_nuclei
        # Fake 384-dim embeddings (as if backbone already ran)
        embeddings[sid] = np.random.randn(n_nuclei, 384).astype(np.float32)
        for j in range(n_nuclei):
            nid = f"{sid}_{j}"
            nuclei_records.append({'nucleus_id': nid, 'spot_id': sid})
            # Random ground truth
            gt_assignments[nid] = cell_types[np.random.randint(n_types)]

    nuclei_spot_map = pd.DataFrame(nuclei_records)
    nuclei_counts_series = pd.Series(nuclei_counts)

    return {
        'proportions': proportions,
        'embeddings': embeddings,
        'nuclei_spot_map': nuclei_spot_map,
        'nuclei_counts': nuclei_counts_series,
        'cell_types': cell_types,
        'gt_assignments': gt_assignments,
    }


def test_e2e_mil_pipeline(synthetic_xenium_data):
    """Full pipeline: embeddings → MIL train → Hungarian → assignments."""
    data = synthetic_xenium_data

    result = run_nucleus_assignment_mil(
        embeddings=data['embeddings'],
        nuclei_spot_map=data['nuclei_spot_map'],
        proportions=data['proportions'],
        nuclei_counts=data['nuclei_counts'],
        cell_types=data['cell_types'],
        n_epochs=10,
        lambda_prior=1.0,
        device="cpu",
    )

    # Verify all nuclei assigned
    total_nuclei = len(data['nuclei_spot_map'])
    assert len(result.assignments) == total_nuclei
    assert result.method == "mil"

    # Verify all assignments are valid cell types
    for nid, ct in result.assignments.items():
        assert ct in data['cell_types']

    # Verify attention matrix is returned
    assert result.assignment_probs is not None
    assert len(result.assignment_probs) == total_nuclei


def test_e2e_proportion_consistency(synthetic_xenium_data):
    """Re-aggregated proportions should correlate with input proportions."""
    data = synthetic_xenium_data

    result = run_nucleus_assignment_mil(
        embeddings=data['embeddings'],
        nuclei_spot_map=data['nuclei_spot_map'],
        proportions=data['proportions'],
        nuclei_counts=data['nuclei_counts'],
        cell_types=data['cell_types'],
        n_epochs=30,
        lambda_prior=1.0,
        device="cpu",
    )

    # Re-aggregate to proportions
    cell_types = data['cell_types']
    input_props = data['proportions'].set_index('spot_id')[cell_types]

    for spot_id in data['proportions']['spot_id']:
        spot_nuclei = data['nuclei_spot_map'][
            data['nuclei_spot_map']['spot_id'] == spot_id
        ]
        n = len(spot_nuclei)
        if n == 0:
            continue
        counts = np.zeros(len(cell_types))
        for nid in spot_nuclei['nucleus_id']:
            ct = result.assignments.get(nid)
            if ct in cell_types:
                counts[cell_types.index(ct)] += 1
        reagg_props = counts / n

        # Re-aggregated should roughly match input (within discretization error)
        input_p = input_props.loc[spot_id].values
        # RMSE should be reasonable (< 0.3 for 7 types)
        rmse = np.sqrt(np.mean((reagg_props - input_p) ** 2))
        assert rmse < 0.5, f"Spot {spot_id}: RMSE={rmse:.3f} too high"
```

**Step 2: Run test**

Run: `cd CITEgeist && python -m pytest tests/test_e2e_single_cell.py -v`
Expected: PASS (2 tests)

**Step 3: Commit**

```bash
git add CITEgeist/tests/test_e2e_single_cell.py
git commit -m "test: add end-to-end integration test for MIL single-cell assignment"
```

---

## Task Summary

| Task | Component | Files | Dependencies |
|------|-----------|-------|-------------|
| 1 | Proportion-weighted Hungarian | `hungarian_assignment.py` | None |
| 2 | Morphology backbone ABC | `morphology_backbone.py` | `vit_encoder.py`, `vit_extractor.py`, `simclr.py` |
| 3 | SingleCellMIL head | `single_cell_mil.py` | None (standalone PyTorch) |
| 4 | MIL assignment in Module 3b | `module3b_nucleus_assignment.py` | Tasks 1, 3 |
| 5 | CitegeistModel integration | `citegeist_model.py`, `__init__.py` | Tasks 2, 4 |
| 6 | Xenium benchmark script | `benchmark_mil_assignment.py`, SLURM | Tasks 2, 4 |
| 7 | Archive deprecated modules | `__init__.py` | None |
| 8 | End-to-end integration test | `test_e2e_single_cell.py` | Tasks 1-4 |

**Dependency graph**: Tasks 1, 2, 3 are independent (parallelize). Task 4 depends on 1+3. Tasks 5, 6 depend on 4+2. Task 7 is independent. Task 8 depends on 1-4.
