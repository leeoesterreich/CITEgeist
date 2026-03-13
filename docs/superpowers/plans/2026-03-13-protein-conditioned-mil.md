# Protein-Conditioned MIL (PC-MIL) Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a learned model that fuses spot-level protein proportions with per-nucleus ViT image features via per-class gated attention to produce accurate single-cell type assignments, solving the MIL attention collapse problem.

**Architecture:** Frozen ViT-S extracts 384-dim features per nucleus patch. A projection head maps to 64-dim. Spot-level protein proportions are encoded to 32-dim and broadcast-concatenated to each nucleus (96-dim fused). Per-class gated attention heads produce (N, K) soft assignments. Multi-task loss: proportion MSE + protein reconstruction MSE + entropy + diversity regularization.

**Tech Stack:** PyTorch, timm (ViT-S), numpy, scipy (Hungarian), sklearn (GMM detection). Existing CITEgeist infrastructure for Module 3 proportions, detection, and constrained assignment.

**Spec:** `docs/superpowers/specs/2026-03-13-protein-conditioned-mil-design.md`

---

## File Structure

| File | Responsibility |
|------|---------------|
| **Create:** `CITEgeist/model/pc_mil.py` | PCMILModel nn.Module — architecture + forward pass |
| **Create:** `CITEgeist/model/pc_mil_training.py` | Training loop, loss functions, validation, early stopping |
| **Create:** `CITEgeist/model/pc_mil_inference.py` | Inference pipeline: detection masking, Hungarian, output formatting |
| **Create:** `CITEgeist/tests/test_pc_mil.py` | Unit tests for model, loss, inference |
| **Create:** `Benchmarking/xenium_benchmarking/CITEgeist/src/benchmark_pc_mil.py` | Xenium 5-fold CV benchmark harness |
| **Create:** `Benchmarking/xenium_benchmarking/CITEgeist/slurm/sbatch_pc_mil.sh` | SLURM job script |
| **Modify:** `CITEgeist/tests/test_asymmetric_loss.py:81,105,135,145` | Fix 4 unpack sites for 5-value return |

---

## Chunk 1: Bug Fix + Core Model

### Task 0: Fix test_asymmetric_loss.py API break

**Files:**
- Modify: `CITEgeist/tests/test_asymmetric_loss.py:81,105,135,145`

- [ ] **Step 1: Read current test file to confirm unpack sites**

Run: `grep -n "optimize_cell_proportions_per_marker" CITEgeist/tests/test_asymmetric_loss.py`

- [ ] **Step 2: Fix all 4 unpack sites**

At line 81, change:
```python
Y, beta, marker_beta_dict, alpha = optimize_cell_proportions_per_marker(
```
to:
```python
Y, beta, marker_beta_dict, alpha, _recon_error = optimize_cell_proportions_per_marker(
```

At line 105, change:
```python
Y_sym, _, _, _ = optimize_cell_proportions_per_marker(
```
to:
```python
Y_sym, _, _, _, _ = optimize_cell_proportions_per_marker(
```

At line 135, change:
```python
Y_sym, _, _, _ = optimize_cell_proportions_per_marker(
```
to:
```python
Y_sym, _, _, _, _ = optimize_cell_proportions_per_marker(
```

At line 145, change:
```python
Y_asym, _, _, _ = optimize_cell_proportions_per_marker(
```
to:
```python
Y_asym, _, _, _, _ = optimize_cell_proportions_per_marker(
```

- [ ] **Step 3: Verify fix compiles**

Run: `cd CITEgeist && python -c "from tests.test_asymmetric_loss import *; print('OK')"`
Expected: `OK` (no import errors)

- [ ] **Step 4: Commit**

```bash
git add CITEgeist/tests/test_asymmetric_loss.py
git commit -m "fix: update test_asymmetric_loss for 5-value return from optimize_cell_proportions_per_marker"
```

---

### Task 1: PCMILModel architecture

**Files:**
- Create: `CITEgeist/model/pc_mil.py`
- Create: `CITEgeist/tests/test_pc_mil.py`

- [ ] **Step 1: Write the failing test**

Create `CITEgeist/tests/test_pc_mil.py`:

```python
"""Tests for Protein-Conditioned MIL model."""
import torch
import torch.nn.functional as F
import numpy as np
import pytest


def test_pcmil_forward_shape():
    """PCMILModel forward pass produces correct output shapes."""
    from model.pc_mil import PCMILModel

    K = 7   # cell types
    M = 15  # protein markers
    N = 12  # nuclei in this spot

    model = PCMILModel(
        image_dim=384,
        n_types=K,
        n_markers=M,
        image_proj_dim=64,
        protein_context_dim=32,
        hidden_dim=128,
    )

    image_features = torch.randn(N, 384)
    protein_props = torch.rand(K)
    protein_props = protein_props / protein_props.sum()

    proportions, attention, reconstructed = model(image_features, protein_props)

    assert proportions.shape == (K,), f"Expected ({K},), got {proportions.shape}"
    assert attention.shape == (N, K), f"Expected ({N},{K}), got {attention.shape}"
    assert reconstructed.shape == (M,), f"Expected ({M},), got {reconstructed.shape}"

    # Proportions should sum to ~1
    assert abs(proportions.sum().item() - 1.0) < 1e-5

    # Each row of attention should sum to ~1 (softmax over types)
    row_sums = attention.sum(dim=1)
    assert torch.allclose(row_sums, torch.ones(N), atol=1e-5)


def test_pcmil_gradient_flow():
    """Gradients flow through all learnable components."""
    from model.pc_mil import PCMILModel

    model = PCMILModel(image_dim=384, n_types=7, n_markers=15)
    image_features = torch.randn(10, 384)
    protein_props = torch.rand(7)
    protein_props = protein_props / protein_props.sum()

    proportions, attention, reconstructed = model(image_features, protein_props)

    # Use proportion loss
    target = torch.rand(7)
    target = target / target.sum()
    loss = torch.nn.functional.mse_loss(proportions, target)
    loss.backward()

    # Check key parameter groups have gradients
    for name, param in model.named_parameters():
        if param.requires_grad:
            assert param.grad is not None, f"No gradient for {name}"
            assert not torch.all(param.grad == 0), f"Zero gradient for {name}"


def test_pcmil_profile_initialization():
    """Protein profile matrix can be initialized from cell_profile_dict."""
    from model.pc_mil import PCMILModel, build_profile_matrix

    cell_profile_dict = {
        "T_cells": ["CD3", "CD4"],
        "B_cells": ["CD19", "CD20"],
        "Macrophages": ["CD68"],
    }
    marker_names = ["CD3", "CD4", "CD19", "CD20", "CD68"]

    profile_matrix = build_profile_matrix(cell_profile_dict, marker_names)

    assert profile_matrix.shape == (3, 5)
    # T_cells should have 1.0 for CD3 (idx 0) and CD4 (idx 1)
    assert profile_matrix[0, 0] == 1.0
    assert profile_matrix[0, 1] == 1.0
    assert profile_matrix[0, 2] == 0.0
    # B_cells should have 1.0 for CD19 (idx 2) and CD20 (idx 3)
    assert profile_matrix[1, 2] == 1.0
    assert profile_matrix[1, 3] == 1.0

    # Use it to initialize model
    model = PCMILModel(
        image_dim=384,
        n_types=3,
        n_markers=5,
        init_profile_matrix=torch.tensor(profile_matrix, dtype=torch.float32),
    )
    assert torch.allclose(
        model.protein_profiles.data,
        torch.tensor(profile_matrix, dtype=torch.float32),
    )


def test_pcmil_variable_nuclei():
    """Model handles different numbers of nuclei per spot."""
    from model.pc_mil import PCMILModel

    model = PCMILModel(image_dim=384, n_types=7, n_markers=15)
    protein_props = torch.rand(7)
    protein_props = protein_props / protein_props.sum()

    for n in [1, 5, 20, 50]:
        feats = torch.randn(n, 384)
        props, attn, recon = model(feats, protein_props)
        assert props.shape == (7,)
        assert attn.shape == (n, 7)


def test_pcmil_softmax_dim():
    """CRITICAL: Verify softmax is over types (dim=1), not nuclei (dim=0)."""
    from model.pc_mil import PCMILModel

    model = PCMILModel(image_dim=384, n_types=3, n_markers=5)
    feats = torch.randn(10, 384)
    protein_props = torch.tensor([0.5, 0.3, 0.2])

    _, attention, _ = model(feats, protein_props)

    # Each nucleus row sums to 1 (softmax over types)
    row_sums = attention.sum(dim=1)
    assert torch.allclose(row_sums, torch.ones(10), atol=1e-5), \
        "Softmax must be over types (dim=1)! Row sums should be 1.0"

    # Column sums should NOT all be equal (would indicate dim=0 bug)
    col_sums = attention.sum(dim=0)
    assert not torch.allclose(col_sums, col_sums[0].expand(3), atol=0.1), \
        "Column sums are suspiciously uniform — possible dim=0 softmax bug"
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cd CITEgeist && python -m pytest tests/test_pc_mil.py -v 2>&1 | head -30`
Expected: FAIL with `ModuleNotFoundError: No module named 'model.pc_mil'`

- [ ] **Step 3: Write PCMILModel implementation**

Create `CITEgeist/model/pc_mil.py`:

```python
"""Protein-Conditioned Multiple Instance Learning (PC-MIL) model.

Fuses spot-level protein proportions with per-nucleus image features
via per-class gated attention to produce single-cell type assignments.

Architecture:
    image_features (N, 384) -> projection -> (N, 64)
    protein_props (K,) -> encoder -> (32,) -> broadcast to N nuclei
    concat -> (N, 96) -> per-class gated attention -> (N, K) logits
    softmax(dim=1) -> attention -> mean(dim=0) -> proportions (K,)
    proportions @ profile_matrix -> reconstructed protein (M,)

See: docs/superpowers/specs/2026-03-13-protein-conditioned-mil-design.md
"""
import torch
import torch.nn as nn
import torch.nn.functional as F
import numpy as np
from typing import Dict, List, Optional, Tuple


def build_profile_matrix(
    cell_profile_dict: Dict[str, List[str]],
    marker_names: List[str],
) -> np.ndarray:
    """Build (K, M) binary profile matrix from cell_profile_dict.

    Args:
        cell_profile_dict: {cell_type: [marker_name, ...]}
        marker_names: Ordered list of all marker names (defines column order)

    Returns:
        (K, M) numpy array where profile[k, m] = 1.0 if marker m belongs to type k
    """
    cell_types = list(cell_profile_dict.keys())
    K = len(cell_types)
    M = len(marker_names)
    marker_to_idx = {name: i for i, name in enumerate(marker_names)}

    profile = np.zeros((K, M), dtype=np.float32)
    for k, ct in enumerate(cell_types):
        for marker in cell_profile_dict[ct]:
            if marker in marker_to_idx:
                profile[k, marker_to_idx[marker]] = 1.0

    return profile


class PCMILModel(nn.Module):
    """Protein-Conditioned MIL for single-cell type assignment.

    Args:
        image_dim: Input dimension from frozen ViT-S (384)
        n_types: Number of cell types K
        n_markers: Number of protein markers M
        image_proj_dim: Image projection output dimension (64)
        protein_context_dim: Protein encoder output dimension (32)
        hidden_dim: Hidden dimension in gated attention (128)
        dropout: Dropout rate in attention heads
        init_profile_matrix: Optional (K, M) tensor to initialize protein profiles
    """

    def __init__(
        self,
        image_dim: int = 384,
        n_types: int = 7,
        n_markers: int = 15,
        image_proj_dim: int = 64,
        protein_context_dim: int = 32,
        hidden_dim: int = 128,
        dropout: float = 0.1,
        init_profile_matrix: Optional[torch.Tensor] = None,
    ):
        super().__init__()
        self.n_types = n_types
        self.n_markers = n_markers

        # Image projection: 384 -> 64
        self.image_projection = nn.Sequential(
            nn.Linear(image_dim, image_proj_dim),
            nn.LayerNorm(image_proj_dim),
        )

        # Protein encoder: K -> 32
        self.protein_encoder = nn.Sequential(
            nn.Linear(n_types, protein_context_dim),
            nn.ReLU(),
            nn.Linear(protein_context_dim, protein_context_dim),
        )

        # Fused dimension
        fused_dim = image_proj_dim + protein_context_dim

        # Per-class gated attention: K separate gate+score heads
        self.W_gate = nn.ModuleList([
            nn.Linear(fused_dim, hidden_dim) for _ in range(n_types)
        ])
        self.W_score = nn.ModuleList([
            nn.Linear(fused_dim, hidden_dim) for _ in range(n_types)
        ])
        self.W_out = nn.ModuleList([
            nn.Sequential(nn.Dropout(dropout), nn.Linear(hidden_dim, 1))
            for _ in range(n_types)
        ])

        # Protein reconstruction profile matrix: (K, M) learnable
        if init_profile_matrix is not None:
            self.protein_profiles = nn.Parameter(init_profile_matrix.clone())
        else:
            self.protein_profiles = nn.Parameter(torch.randn(n_types, n_markers) * 0.1)

        # Xavier init for projection and attention weights
        self._init_weights()

    def _init_weights(self):
        """Xavier initialization for linear layers."""
        for module in [self.image_projection, self.protein_encoder]:
            for m in module:
                if isinstance(m, nn.Linear):
                    nn.init.xavier_uniform_(m.weight)
                    nn.init.zeros_(m.bias)
        for k in range(self.n_types):
            nn.init.xavier_uniform_(self.W_gate[k].weight)
            nn.init.zeros_(self.W_gate[k].bias)
            nn.init.xavier_uniform_(self.W_score[k].weight)
            nn.init.zeros_(self.W_score[k].bias)
            for m in self.W_out[k]:
                if isinstance(m, nn.Linear):
                    nn.init.xavier_uniform_(m.weight)
                    nn.init.zeros_(m.bias)

    def forward_with_logits(
        self,
        image_features: torch.Tensor,
        protein_proportions: torch.Tensor,
    ) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor]:
        """Forward pass returning pre-softmax logits (for inference masking).

        Args:
            image_features: (N, image_dim) pre-extracted ViT features per nucleus
            protein_proportions: (K,) spot-level protein proportions

        Returns:
            logits: (N, K) pre-softmax logits
            attention: (N, K) per-nucleus type probability matrix
            proportions: (K,) predicted spot-level proportions
            reconstructed: (M,) reconstructed protein signal
        """
        N = image_features.shape[0]

        # Project image features
        image_emb = self.image_projection(image_features)  # (N, 64)

        # Encode protein context and broadcast to all nuclei
        protein_context = self.protein_encoder(protein_proportions)  # (32,)
        protein_context_broadcast = protein_context.unsqueeze(0).expand(N, -1)  # (N, 32)

        # Fuse: concat image + protein
        fused = torch.cat([image_emb, protein_context_broadcast], dim=1)  # (N, 96)

        # Per-class gated attention
        logits_list = []
        for k in range(self.n_types):
            gate_k = torch.sigmoid(self.W_gate[k](fused))    # (N, hidden)
            score_k = torch.tanh(self.W_score[k](fused))     # (N, hidden)
            gated = gate_k * score_k                          # (N, hidden)
            logit_k = self.W_out[k](gated)                    # (N, 1)
            logits_list.append(logit_k)

        logits = torch.cat(logits_list, dim=1)  # (N, K)

        # Softmax over types (dim=1) — CRITICAL: NOT dim=0
        attention = F.softmax(logits, dim=1)  # (N, K)

        # Spot-level proportions
        proportions = attention.mean(dim=0)  # (K,)

        # Protein reconstruction
        reconstructed = proportions @ self.protein_profiles  # (M,)

        return logits, attention, proportions, reconstructed

    def forward(
        self,
        image_features: torch.Tensor,
        protein_proportions: torch.Tensor,
    ) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        """Forward pass for a single spot.

        Args:
            image_features: (N, image_dim) pre-extracted ViT features per nucleus
            protein_proportions: (K,) spot-level protein proportions

        Returns:
            proportions: (K,) predicted spot-level proportions
            attention: (N, K) per-nucleus type probability matrix
            reconstructed: (M,) reconstructed protein signal
        """
        _, attention, proportions, reconstructed = self.forward_with_logits(
            image_features, protein_proportions,
        )
        return proportions, attention, reconstructed
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `cd CITEgeist && python -m pytest tests/test_pc_mil.py -v`
Expected: All 5 tests PASS

- [ ] **Step 5: Commit**

```bash
git add CITEgeist/model/pc_mil.py CITEgeist/tests/test_pc_mil.py
git commit -m "feat: add PCMILModel architecture with per-class gated attention"
```

---

## Chunk 2: Loss Functions + Training Pipeline

### Task 2: Loss functions

**Files:**
- Modify: `CITEgeist/tests/test_pc_mil.py` (add loss tests)
- Create: `CITEgeist/model/pc_mil_training.py`

- [ ] **Step 1: Write failing loss tests**

Append to `CITEgeist/tests/test_pc_mil.py`:

```python
def test_pcmil_loss_components():
    """All 4 loss components compute without errors."""
    from model.pc_mil_training import compute_pc_mil_loss

    K, M, N = 7, 15, 10
    pred_props = torch.rand(K)
    pred_props = pred_props / pred_props.sum()
    true_props = torch.rand(K)
    true_props = true_props / true_props.sum()
    attention = F.softmax(torch.randn(N, K), dim=1)
    reconstructed = torch.randn(M)
    observed_protein = torch.randn(M)

    loss, components = compute_pc_mil_loss(
        pred_proportions=pred_props,
        true_proportions=true_props,
        attention=attention,
        reconstructed_protein=reconstructed,
        observed_protein=observed_protein,
    )

    assert loss.ndim == 0, "Loss should be scalar"
    assert not torch.isnan(loss), "Loss should not be NaN"
    assert "proportion" in components
    assert "reconstruction" in components
    assert "entropy" in components
    assert "diversity" in components


def test_diversity_loss_penalizes_absent_types():
    """Diversity loss increases when active types have near-zero prediction."""
    from model.pc_mil_training import compute_pc_mil_loss

    K, M, N = 3, 5, 10
    true_props = torch.tensor([0.4, 0.4, 0.2])  # all types present
    observed = torch.randn(M)

    # Good prediction: all types represented
    attn_good = F.softmax(torch.randn(N, K), dim=1)
    props_good = attn_good.mean(dim=0)
    recon_good = torch.randn(M)
    _, comp_good = compute_pc_mil_loss(props_good, true_props, attn_good, recon_good, observed)

    # Bad prediction: one type collapsed to 0
    attn_bad = torch.zeros(N, K)
    attn_bad[:, 0] = 0.7
    attn_bad[:, 1] = 0.3  # type 2 is absent
    props_bad = attn_bad.mean(dim=0)
    recon_bad = torch.randn(M)
    _, comp_bad = compute_pc_mil_loss(props_bad, true_props, attn_bad, recon_bad, observed)

    assert comp_bad["diversity"] > comp_good["diversity"], \
        "Diversity loss should be higher when an active type is predicted as absent"


def test_entropy_loss_penalizes_uniform():
    """Entropy loss increases when assignments are uniform (high entropy)."""
    from model.pc_mil_training import compute_pc_mil_loss

    K, M, N = 5, 10, 10
    true_props = torch.ones(K) / K
    observed = torch.randn(M)

    # Confident: each nucleus assigned to one type
    attn_confident = torch.zeros(N, K)
    for i in range(N):
        attn_confident[i, i % K] = 1.0
    # Add small epsilon for log stability
    attn_confident = attn_confident.clamp(min=1e-6)
    attn_confident = attn_confident / attn_confident.sum(dim=1, keepdim=True)

    # Uniform: all nuclei spread evenly
    attn_uniform = torch.ones(N, K) / K

    _, comp_confident = compute_pc_mil_loss(
        attn_confident.mean(0), true_props, attn_confident, torch.randn(M), observed
    )
    _, comp_uniform = compute_pc_mil_loss(
        attn_uniform.mean(0), true_props, attn_uniform, torch.randn(M), observed
    )

    assert comp_uniform["entropy"] > comp_confident["entropy"], \
        "Entropy loss should be higher for uniform assignments"
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cd CITEgeist && python -m pytest tests/test_pc_mil.py::test_pcmil_loss_components -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'model.pc_mil_training'`

- [ ] **Step 3: Write loss functions and training loop**

Create `CITEgeist/model/pc_mil_training.py`:

```python
"""Training pipeline for Protein-Conditioned MIL.

Provides:
- compute_pc_mil_loss(): Multi-task loss with 4 components
- SpotDataset: PyTorch dataset for spot-level training data
- train_pc_mil(): Full training loop with early stopping
"""
import logging
from typing import Dict, List, Optional, Tuple

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import Dataset

from .pc_mil import PCMILModel

logger = logging.getLogger(__name__)


def compute_pc_mil_loss(
    pred_proportions: torch.Tensor,
    true_proportions: torch.Tensor,
    attention: torch.Tensor,
    reconstructed_protein: torch.Tensor,
    observed_protein: torch.Tensor,
    lambda_recon: float = 1.0,
    lambda_entropy: float = 0.1,
    lambda_diversity: float = 0.05,
    mask: Optional[torch.Tensor] = None,
) -> Tuple[torch.Tensor, Dict[str, float]]:
    """Compute multi-task PC-MIL loss.

    Args:
        pred_proportions: (K,) predicted spot proportions
        true_proportions: (K,) target proportions from Module 3
        attention: (N, K) per-nucleus attention weights
        reconstructed_protein: (M,) reconstructed protein signal
        observed_protein: (M,) observed CLR-normalized protein signal
        lambda_recon: Weight for reconstruction loss (default 1.0)
        lambda_entropy: Weight for entropy regularization (default 0.1)
        lambda_diversity: Weight for diversity loss (default 0.05)
        mask: Optional (N,) bool mask for padded nuclei (True = valid)

    Returns:
        total_loss: Scalar loss
        components: Dict of individual loss values for logging
    """
    eps = 1e-8

    # Apply padding mask if provided
    if mask is not None:
        attention = attention[mask]

    # 1. Proportion loss (MSE)
    l_prop = F.mse_loss(pred_proportions, true_proportions)

    # 2. Protein reconstruction loss (MSE)
    l_recon = F.mse_loss(reconstructed_protein, observed_protein)

    # 3. Entropy regularization (minimize = sharper assignments)
    per_nucleus_entropy = -(attention * torch.log(attention + eps)).sum(dim=1)
    l_entropy = per_nucleus_entropy.mean()

    # 4. Diversity loss (penalize absent active types)
    mean_attention = attention.mean(dim=0)  # (K,) = predicted proportions
    active_mask = (true_proportions > 0.01).float()
    # +0.01 floor preserves gradient (unlike clamp which zeros it)
    l_diversity = -(active_mask * torch.log(mean_attention + 0.01)).sum()

    # Combined
    total = l_prop + lambda_recon * l_recon + lambda_entropy * l_entropy + lambda_diversity * l_diversity

    components = {
        "proportion": l_prop.item(),
        "reconstruction": l_recon.item(),
        "entropy": l_entropy.item(),
        "diversity": l_diversity.item(),
        "total": total.item(),
    }

    return total, components


class SpotDataset(Dataset):
    """Dataset yielding per-spot training tuples.

    Each item contains:
    - image_features: (N_i, 384) pre-extracted ViT features
    - protein_props: (K,) Module 3 proportions (conditioning input)
    - protein_signal: (M,) CLR-normalized protein signal (reconstruction target)
    - true_props: (K,) ground truth proportions (supervision target)
    - n_nuclei: int, number of nuclei in this spot
    """

    def __init__(
        self,
        features_per_spot: List[np.ndarray],
        protein_props: np.ndarray,
        protein_signals: np.ndarray,
        true_props: np.ndarray,
        spot_weights: Optional[np.ndarray] = None,
    ):
        """Initialize SpotDataset.

        Args:
            features_per_spot: List of (N_i, 384) arrays, one per spot
            protein_props: (n_spots, K) Module 3 proportions
            protein_signals: (n_spots, M) CLR-normalized protein signals
            true_props: (n_spots, K) ground truth proportions
            spot_weights: Optional (n_spots,) inverse-frequency weights
        """
        self.features = features_per_spot
        self.protein_props = protein_props
        self.protein_signals = protein_signals
        self.true_props = true_props
        self.weights = spot_weights if spot_weights is not None else np.ones(len(features_per_spot))

    def __len__(self):
        return len(self.features)

    def __getitem__(self, idx):
        return {
            "image_features": torch.tensor(self.features[idx], dtype=torch.float32),
            "protein_props": torch.tensor(self.protein_props[idx], dtype=torch.float32),
            "protein_signal": torch.tensor(self.protein_signals[idx], dtype=torch.float32),
            "true_props": torch.tensor(self.true_props[idx], dtype=torch.float32),
            "weight": torch.tensor(self.weights[idx], dtype=torch.float32),
        }


def compute_inverse_frequency_weights(
    true_props: np.ndarray,
    min_weight: float = 1.0,
    max_weight: float = 10.0,
) -> np.ndarray:
    """Compute per-spot weights based on rarest type present.

    Args:
        true_props: (n_spots, K) ground truth proportions
        min_weight: Minimum weight
        max_weight: Maximum weight cap

    Returns:
        (n_spots,) weights
    """
    # For each spot, find the rarest type present (prop > 0.01)
    present = true_props > 0.01
    # Global frequency of each type across spots
    type_freq = present.mean(axis=0)  # (K,)
    type_freq = np.clip(type_freq, 0.01, None)  # avoid div by 0

    weights = np.ones(true_props.shape[0])
    for i in range(true_props.shape[0]):
        present_types = np.where(present[i])[0]
        if len(present_types) > 0:
            rarest_freq = type_freq[present_types].min()
            weights[i] = 1.0 / rarest_freq

    # Normalize and clip
    weights = weights / weights.mean()
    weights = np.clip(weights, min_weight, max_weight)
    return weights


def train_pc_mil(
    model: PCMILModel,
    train_dataset: SpotDataset,
    val_dataset: SpotDataset,
    n_epochs: int = 200,
    lr: float = 1e-3,
    weight_decay: float = 1e-4,
    lambda_recon: float = 1.0,
    lambda_entropy: float = 0.1,
    lambda_diversity: float = 0.05,
    patience: int = 20,
    grad_clip: float = 1.0,
    device: str = "cpu",
    save_path: Optional[str] = None,
) -> Dict[str, list]:
    """Train PC-MIL model.

    Args:
        model: PCMILModel instance
        train_dataset: Training SpotDataset
        val_dataset: Validation SpotDataset
        n_epochs: Max training epochs
        lr: Learning rate
        weight_decay: AdamW weight decay
        lambda_recon: Reconstruction loss weight
        lambda_entropy: Entropy loss weight
        lambda_diversity: Diversity loss weight
        patience: Early stopping patience (epochs without improvement)
        grad_clip: Max gradient norm
        device: 'cuda' or 'cpu'
        save_path: Path to save best model checkpoint

    Returns:
        History dict with per-epoch metrics
    """
    model.to(device)
    optimizer = torch.optim.AdamW(model.parameters(), lr=lr, weight_decay=weight_decay)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=n_epochs, eta_min=1e-5)

    history = {
        "train_loss": [], "val_loss": [], "val_r": [],
        "active_types": [], "mean_entropy": [],
    }
    best_val_r = -1.0
    epochs_no_improve = 0

    for epoch in range(n_epochs):
        # --- Training ---
        model.train()
        epoch_loss = 0.0
        n_train = len(train_dataset)

        # Stratified shuffle (simple: just random order)
        indices = np.random.permutation(n_train)

        for idx in indices:
            sample = train_dataset[idx]
            img_feats = sample["image_features"].to(device)
            prot_props = sample["protein_props"].to(device)
            prot_signal = sample["protein_signal"].to(device)
            true_props = sample["true_props"].to(device)
            weight = sample["weight"].to(device)

            # Skip spots with 0 nuclei
            if img_feats.shape[0] == 0:
                continue

            proportions, attention, reconstructed = model(img_feats, prot_props)

            loss, _ = compute_pc_mil_loss(
                proportions, true_props, attention, reconstructed, prot_signal,
                lambda_recon=lambda_recon,
                lambda_entropy=lambda_entropy,
                lambda_diversity=lambda_diversity,
            )
            loss = loss * weight

            optimizer.zero_grad()
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), grad_clip)
            optimizer.step()

            epoch_loss += loss.item()

        scheduler.step()
        history["train_loss"].append(epoch_loss / max(n_train, 1))

        # --- Validation ---
        model.eval()
        val_loss = 0.0
        all_pred, all_target = [], []
        active_count = 0
        total_entropy = 0.0
        n_val = len(val_dataset)

        with torch.no_grad():
            for idx in range(n_val):
                sample = val_dataset[idx]
                img_feats = sample["image_features"].to(device)
                prot_props = sample["protein_props"].to(device)
                prot_signal = sample["protein_signal"].to(device)
                true_props = sample["true_props"].to(device)

                if img_feats.shape[0] == 0:
                    continue

                proportions, attention, reconstructed = model(img_feats, prot_props)

                loss, comp = compute_pc_mil_loss(
                    proportions, true_props, attention, reconstructed, prot_signal,
                    lambda_recon=lambda_recon,
                    lambda_entropy=lambda_entropy,
                    lambda_diversity=lambda_diversity,
                )
                val_loss += loss.item()

                all_pred.append(proportions.cpu().numpy())
                all_target.append(true_props.cpu().numpy())

                # Collapse detection: count types with any nucleus having max attn > 0.3
                max_per_type = attention.max(dim=0).values
                active_count += (max_per_type > 0.3).sum().item()
                total_entropy += comp["entropy"]

        history["val_loss"].append(val_loss / max(n_val, 1))
        history["active_types"].append(active_count / max(n_val, 1))
        history["mean_entropy"].append(total_entropy / max(n_val, 1))

        # Pearson r
        if all_pred:
            pred_flat = np.concatenate(all_pred)
            target_flat = np.concatenate(all_target)
            if pred_flat.std() > 0 and target_flat.std() > 0:
                r = float(np.corrcoef(pred_flat, target_flat)[0, 1])
            else:
                r = 0.0
        else:
            r = 0.0
        history["val_r"].append(r)

        # Early stopping
        if r > best_val_r:
            best_val_r = r
            epochs_no_improve = 0
            if save_path:
                torch.save(model.state_dict(), save_path)
        else:
            epochs_no_improve += 1

        if (epoch + 1) % 10 == 0:
            logger.info(
                "Epoch %d/%d: train=%.4f val=%.4f r=%.4f active=%.1f entropy=%.4f",
                epoch + 1, n_epochs, history["train_loss"][-1],
                history["val_loss"][-1], r,
                history["active_types"][-1], history["mean_entropy"][-1],
            )

        if epochs_no_improve >= patience:
            logger.info("Early stopping at epoch %d (patience=%d)", epoch + 1, patience)
            break

    return history
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `cd CITEgeist && python -m pytest tests/test_pc_mil.py -v`
Expected: All 8 tests PASS

- [ ] **Step 5: Commit**

```bash
git add CITEgeist/model/pc_mil_training.py CITEgeist/tests/test_pc_mil.py
git commit -m "feat: add PC-MIL training pipeline with multi-task loss"
```

---

## Chunk 3: Inference Pipeline

### Task 3: Inference with detection masking and Hungarian assignment

**Files:**
- Create: `CITEgeist/model/pc_mil_inference.py`
- Modify: `CITEgeist/tests/test_pc_mil.py` (add inference tests)

- [ ] **Step 1: Write failing inference tests**

Append to `CITEgeist/tests/test_pc_mil.py`:

```python
def test_pcmil_inference_basic():
    """Inference produces per-nucleus assignments with correct format."""
    from model.pc_mil_inference import pc_mil_infer_spot

    from model.pc_mil import PCMILModel

    model = PCMILModel(image_dim=384, n_types=3, n_markers=5)
    model.eval()

    N = 8
    image_features = torch.randn(N, 384)
    protein_props = torch.tensor([0.5, 0.3, 0.2])
    detected = np.array([True, True, True])

    result = pc_mil_infer_spot(
        model=model,
        image_features=image_features,
        protein_proportions=protein_props,
        detected_types=detected,
        cell_type_names=["TypeA", "TypeB", "TypeC"],
    )

    assert len(result) == N
    assert "cell_type" in result.columns
    assert "confidence" in result.columns
    assert set(result["cell_type"].unique()).issubset({"TypeA", "TypeB", "TypeC"})


def test_pcmil_inference_detection_mask():
    """Undetected types get zero assignments."""
    from model.pc_mil_inference import pc_mil_infer_spot
    from model.pc_mil import PCMILModel

    model = PCMILModel(image_dim=384, n_types=3, n_markers=5)
    model.eval()

    N = 10
    image_features = torch.randn(N, 384)
    protein_props = torch.tensor([0.5, 0.3, 0.2])
    # Only type 0 detected
    detected = np.array([True, False, False])

    result = pc_mil_infer_spot(
        model=model,
        image_features=image_features,
        protein_proportions=protein_props,
        detected_types=detected,
        cell_type_names=["TypeA", "TypeB", "TypeC"],
    )

    # All assignments should be TypeA (only detected type)
    assert (result["cell_type"] == "TypeA").all()


def test_pcmil_inference_all_false_detection():
    """All-false detection mask falls back to no masking (no NaN)."""
    from model.pc_mil_inference import pc_mil_infer_spot
    from model.pc_mil import PCMILModel

    model = PCMILModel(image_dim=384, n_types=3, n_markers=5)
    model.eval()

    N = 5
    image_features = torch.randn(N, 384)
    protein_props = torch.tensor([0.5, 0.3, 0.2])
    detected = np.array([False, False, False])  # nothing detected

    result = pc_mil_infer_spot(
        model=model,
        image_features=image_features,
        protein_proportions=protein_props,
        detected_types=detected,
        cell_type_names=["TypeA", "TypeB", "TypeC"],
    )

    assert len(result) == N
    assert not result["confidence"].isna().any(), "Should not have NaN confidences"
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cd CITEgeist && python -m pytest tests/test_pc_mil.py::test_pcmil_inference_basic -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'model.pc_mil_inference'`

- [ ] **Step 3: Write inference module**

Create `CITEgeist/model/pc_mil_inference.py`:

```python
"""Inference pipeline for Protein-Conditioned MIL.

Provides per-spot inference with:
- Detection mask from GMM (model/detection.py)
- Hungarian assignment for discrete cell typing
- Output formatting with confidence scores
"""
import logging
from typing import Dict, List, Optional

import numpy as np
import pandas as pd
import torch
import torch.nn.functional as F

from .constrained_assignment import hungarian_assign
from .morphology_features import largest_remainder_discretize

logger = logging.getLogger(__name__)


def pc_mil_infer_spot(
    model: torch.nn.Module,
    image_features: torch.Tensor,
    protein_proportions: torch.Tensor,
    detected_types: np.ndarray,
    cell_type_names: List[str],
    nucleus_ids: Optional[List[str]] = None,
    barcode: Optional[str] = None,
) -> pd.DataFrame:
    """Run PC-MIL inference for a single spot.

    Args:
        model: Trained PCMILModel (in eval mode)
        image_features: (N, image_dim) pre-extracted ViT features
        protein_proportions: (K,) spot-level protein proportions
        detected_types: (K,) boolean mask from detection.detect_cell_types()
        cell_type_names: List of K cell type names (ordered)
        nucleus_ids: Optional list of N nucleus IDs
        barcode: Optional spot barcode

    Returns:
        DataFrame with columns: nucleus_id, barcode, cell_type, confidence, protein_score
    """
    N = image_features.shape[0]
    K = len(cell_type_names)
    device = next(model.parameters()).device

    if nucleus_ids is None:
        nucleus_ids = [f"nuc_{i:04d}" for i in range(N)]

    with torch.no_grad():
        img_feats = image_features.to(device)
        prot_props = protein_proportions.to(device)

        # Get pre-softmax logits via forward_with_logits
        logits, _, _, _ = model.forward_with_logits(img_feats, prot_props)

        # Apply detection mask (with all-false fallback)
        if detected_types.any():
            mask_tensor = torch.tensor(detected_types, dtype=torch.bool, device=device)
            logits[:, ~mask_tensor] = float("-inf")

        # Re-apply softmax with mask
        attention = F.softmax(logits, dim=1)  # (N, K)
        proportions = attention.mean(dim=0)   # (K,)

    # Convert to numpy
    attention_np = attention.cpu().numpy()  # (N, K)
    proportions_np = proportions.cpu().numpy()
    protein_props_np = protein_proportions.cpu().numpy() if isinstance(protein_proportions, torch.Tensor) else protein_proportions

    # Convert proportions to integer counts
    counts = largest_remainder_discretize(proportions_np, N)

    # Hungarian assignment using attention as log-likelihoods
    # hungarian_assign expects log_likes and internally negates for cost
    log_likes = np.log(attention_np + 1e-8)
    assignments = hungarian_assign(log_likes, counts)

    # Build output DataFrame
    rows = []
    for i in range(N):
        assigned_type = cell_type_names[assignments[i]]
        confidence = float(attention_np[i, assignments[i]])
        protein_score = float(protein_props_np[assignments[i]])

        rows.append({
            "nucleus_id": nucleus_ids[i],
            "barcode": barcode or "",
            "cell_type": assigned_type,
            "confidence": confidence,
            "protein_score": protein_score,
        })

    return pd.DataFrame(rows)
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `cd CITEgeist && python -m pytest tests/test_pc_mil.py -v`
Expected: All 11 tests PASS

- [ ] **Step 5: Commit**

```bash
git add CITEgeist/model/pc_mil_inference.py CITEgeist/tests/test_pc_mil.py
git commit -m "feat: add PC-MIL inference pipeline with detection masking and Hungarian assignment"
```

---

## Chunk 4: Xenium Benchmark Harness + SLURM

### Task 4: Benchmark script for Xenium 5-fold CV

**Files:**
- Create: `Benchmarking/xenium_benchmarking/CITEgeist/src/benchmark_pc_mil.py`
- Create: `Benchmarking/xenium_benchmarking/CITEgeist/slurm/sbatch_pc_mil.sh`

- [ ] **Step 1: Write the benchmark harness**

Create `Benchmarking/xenium_benchmarking/CITEgeist/src/benchmark_pc_mil.py`:

```python
#!/usr/bin/env python
"""
Benchmark PC-MIL on Xenium pseudo-Visium with 5-fold cross-validation.

Pipeline:
1. Load pre-extracted ViT-S features per nucleus
2. Load Module 3 hybrid proportions (conditioning input)
3. Load CLR-normalized protein signals (reconstruction target)
4. Load ground truth proportions + single-cell labels
5. 5-fold CV: train PC-MIL, evaluate proportion r + single-cell accuracy

Usage:
    python benchmark_pc_mil.py \
        --output-dir output/pc_mil \
        --n-epochs 200 \
        --device cuda
"""
import argparse
import json
import logging
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
import torch

# Setup paths
REPO_ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO_ROOT))

from CITEgeist.model.pc_mil import PCMILModel, build_profile_matrix
from CITEgeist.model.pc_mil_training import (
    SpotDataset, compute_inverse_frequency_weights, train_pc_mil,
)
from CITEgeist.model.pc_mil_inference import pc_mil_infer_spot
from CITEgeist.model.detection import detect_cell_types
from CITEgeist.model.constrained_assignment import random_assign

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[logging.StreamHandler(sys.stderr)],
)
logger = logging.getLogger(__name__)

# Data paths
PSEUDOVISIUM_DIR = REPO_ROOT / "Benchmarking/xenium_pseudovisium/data_protein_gt"
HYBRID_OUTPUT_DIR = REPO_ROOT / "Benchmarking/xenium_benchmarking/CITEgeist/output/hybrid_detection_filter"

# Cell types (achievable-7 for Xenium)
CELL_TYPES = [
    "B cells", "CD4+ T cells", "CD8+ T cells", "Macrophages",
    "Endothelial", "Epithelial", "Fibroblasts",
]

# 14 regions, 5 folds
FOLDS = [
    {"train": [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10], "val": [11, 12, 13]},
    {"train": [0, 1, 2, 3, 4, 5, 6, 7, 11, 12, 13], "val": [8, 9, 10]},
    {"train": [0, 1, 2, 3, 4, 8, 9, 10, 11, 12, 13], "val": [5, 6, 7]},
    {"train": [0, 1, 5, 6, 7, 8, 9, 10, 11, 12, 13], "val": [2, 3, 4]},
    {"train": [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13], "val": [0, 1]},
]


def load_region_data(region_id: int, features_dir: Path) -> Dict:
    """Load all data for one pseudo-Visium region.

    Returns dict with keys:
        features_per_spot: List of (N_i, 384) arrays per spot
        protein_props: (n_spots, K) Module 3 hybrid proportions
        protein_signals: (n_spots, M) CLR-normalized protein signals
        true_props: (n_spots, K) ground truth proportions
        cell_type_labels: List of (N_i,) string arrays per spot (GT labels)
        nucleus_ids: List of (N_i,) string arrays per spot
        spot_barcodes: (n_spots,) spot barcode strings
    """
    # This function needs to be adapted to actual data layout
    # Placeholder structure — will be filled during Task 6
    raise NotImplementedError(
        "Implement based on actual data layout in "
        f"{PSEUDOVISIUM_DIR} and {features_dir}"
    )


def evaluate_fold(
    model: PCMILModel,
    val_data: Dict,
    cell_type_names: List[str],
    device: str,
) -> Dict:
    """Evaluate model on validation fold.

    Returns dict with proportion_r, single_cell_accuracy, per_type_accuracy.
    Uses pc_mil_infer_spot for both proportion and single-cell evaluation
    to ensure consistency (detection masking applied uniformly).
    """
    model.eval()

    all_pred_props = []
    all_true_props = []
    correct = 0
    total = 0
    per_type = {ct: {"correct": 0, "total": 0} for ct in cell_type_names}

    for i in range(len(val_data["features_per_spot"])):
        feats = val_data["features_per_spot"][i]
        if len(feats) == 0:
            continue

        img_feats = torch.tensor(feats, dtype=torch.float32)
        prot_props = torch.tensor(val_data["protein_props"][i], dtype=torch.float32)

        # Detection mask
        detected = val_data.get("detected", np.ones(len(cell_type_names), dtype=bool))
        if detected.ndim == 2:
            detected = detected[i]

        result_df = pc_mil_infer_spot(
            model=model,
            image_features=img_feats,
            protein_proportions=prot_props,
            detected_types=detected,
            cell_type_names=cell_type_names,
        )

        # Proportion evaluation — derive from inference attention (consistent masking)
        # Count assigned types to get predicted proportions
        type_counts = result_df["cell_type"].value_counts()
        pred_p = np.zeros(len(cell_type_names))
        for k, ct in enumerate(cell_type_names):
            pred_p[k] = type_counts.get(ct, 0) / len(result_df)
        all_pred_props.append(pred_p)
        all_true_props.append(val_data["true_props"][i])

        # Single-cell accuracy (if GT labels available)
        if "cell_type_labels" in val_data:
            gt_labels = val_data["cell_type_labels"][i]
            pred_labels = result_df["cell_type"].values
            for j in range(len(gt_labels)):
                if j < len(pred_labels):
                    total += 1
                    gt = gt_labels[j]
                    pred = pred_labels[j]
                    if gt == pred:
                        correct += 1
                    if gt in per_type:
                        per_type[gt]["total"] += 1
                        if gt == pred:
                            per_type[gt]["correct"] += 1

    # Proportion Pearson r
    pred_flat = np.concatenate(all_pred_props)
    true_flat = np.concatenate(all_true_props)
    if pred_flat.std() > 0 and true_flat.std() > 0:
        prop_r = float(np.corrcoef(pred_flat, true_flat)[0, 1])
    else:
        prop_r = 0.0

    sc_accuracy = correct / total if total > 0 else 0.0
    random_accuracy = 1.0 / len(cell_type_names) if total > 0 else 0.0

    per_type_acc = {}
    for ct in cell_type_names:
        t = per_type[ct]["total"]
        c = per_type[ct]["correct"]
        per_type_acc[ct] = c / t if t > 0 else 0.0

    return {
        "proportion_r": prop_r,
        "single_cell_accuracy": sc_accuracy,
        "random_accuracy": random_accuracy,
        "per_type_accuracy": per_type_acc,
        "n_spots": len(all_pred_props),
        "n_cells": total,
    }


def main():
    parser = argparse.ArgumentParser(description="PC-MIL Xenium Benchmark")
    parser.add_argument("--output-dir", type=str, required=True)
    parser.add_argument("--features-dir", type=str, required=True,
                        help="Directory with pre-extracted ViT features per region")
    parser.add_argument("--n-epochs", type=int, default=200)
    parser.add_argument("--lr", type=float, default=1e-3)
    parser.add_argument("--lambda-recon", type=float, default=1.0)
    parser.add_argument("--lambda-entropy", type=float, default=0.1)
    parser.add_argument("--lambda-diversity", type=float, default=0.05)
    parser.add_argument("--patience", type=int, default=20)
    parser.add_argument("--device", type=str, default="cuda")
    parser.add_argument("--fold", type=int, default=None, help="Run single fold (0-4)")
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    features_dir = Path(args.features_dir)

    folds_to_run = [args.fold] if args.fold is not None else range(5)
    all_results = []

    for fold_idx in folds_to_run:
        fold = FOLDS[fold_idx]
        logger.info(f"=== Fold {fold_idx}: train={fold['train']}, val={fold['val']} ===")

        # Load data for train/val regions
        # NOTE: load_region_data must be implemented for actual data layout
        train_features, train_props, train_signals, train_true = [], [], [], []
        for rid in fold["train"]:
            data = load_region_data(rid, features_dir)
            train_features.extend(data["features_per_spot"])
            train_props.append(data["protein_props"])
            train_signals.append(data["protein_signals"])
            train_true.append(data["true_props"])

        train_props = np.concatenate(train_props)
        train_signals = np.concatenate(train_signals)
        train_true = np.concatenate(train_true)

        # Inverse frequency weights
        weights = compute_inverse_frequency_weights(train_true)

        train_ds = SpotDataset(train_features, train_props, train_signals, train_true, weights)

        # Val — collect per-spot lists, then concatenate arrays for SpotDataset
        val_features = []
        val_props_list, val_signals_list, val_true_list = [], [], []
        val_labels_list = []
        for rid in fold["val"]:
            data = load_region_data(rid, features_dir)
            val_features.extend(data["features_per_spot"])
            val_props_list.append(data["protein_props"])
            val_signals_list.append(data["protein_signals"])
            val_true_list.append(data["true_props"])
            if "cell_type_labels" in data:
                val_labels_list.extend(data["cell_type_labels"])

        val_props = np.concatenate(val_props_list)
        val_true = np.concatenate(val_true_list)
        val_signals = np.concatenate(val_signals_list)

        val_ds = SpotDataset(val_features, val_props, val_signals, val_true)

        # Build val_data dict for evaluate_fold (per-spot indexing)
        val_data_combined = {
            "features_per_spot": val_features,
            "protein_props": val_props,   # (n_spots, K) — indexable by spot
            "true_props": val_true,       # (n_spots, K)
        }
        if val_labels_list:
            val_data_combined["cell_type_labels"] = val_labels_list

        # Model
        K = len(CELL_TYPES)
        M = train_signals.shape[1]
        model = PCMILModel(
            image_dim=384, n_types=K, n_markers=M,
        )

        save_path = str(output_dir / f"pc_mil_fold{fold_idx}.pt")

        history = train_pc_mil(
            model=model,
            train_dataset=train_ds,
            val_dataset=val_ds,
            n_epochs=args.n_epochs,
            lr=args.lr,
            lambda_recon=args.lambda_recon,
            lambda_entropy=args.lambda_entropy,
            lambda_diversity=args.lambda_diversity,
            patience=args.patience,
            device=args.device,
            save_path=save_path,
        )

        # Evaluate best model
        model.load_state_dict(torch.load(save_path, weights_only=True))
        model.to(args.device)

        fold_results = evaluate_fold(
            model, val_data_combined, CELL_TYPES, args.device,
        )
        fold_results["fold"] = fold_idx
        fold_results["history"] = history
        all_results.append(fold_results)

        logger.info(f"Fold {fold_idx}: prop_r={fold_results['proportion_r']:.4f}, "
                     f"sc_acc={fold_results['single_cell_accuracy']:.4f}")

        # Save fold results
        with open(output_dir / f"fold{fold_idx}_results.json", "w") as f:
            json.dump({k: v for k, v in fold_results.items() if k != "history"},
                      f, indent=2, default=str)

    # Summary
    if all_results:
        mean_r = np.mean([r["proportion_r"] for r in all_results])
        mean_acc = np.mean([r["single_cell_accuracy"] for r in all_results])
        logger.info(f"=== SUMMARY: mean_prop_r={mean_r:.4f}, mean_sc_acc={mean_acc:.4f} ===")

        summary = {
            "mean_proportion_r": float(mean_r),
            "mean_single_cell_accuracy": float(mean_acc),
            "folds": [{k: v for k, v in r.items() if k != "history"} for r in all_results],
        }
        with open(output_dir / "summary.json", "w") as f:
            json.dump(summary, f, indent=2, default=str)


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Write SLURM script**

Create `Benchmarking/xenium_benchmarking/CITEgeist/slurm/sbatch_pc_mil.sh`:

```bash
#!/bin/bash
#SBATCH --job-name=pc_mil
#SBATCH --output=Benchmarking/xenium_benchmarking/CITEgeist/slurm/logs/pc_mil_%j.out
#SBATCH --error=Benchmarking/xenium_benchmarking/CITEgeist/slurm/logs/pc_mil_%j.err
#SBATCH --partition=l40s
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

set -euo pipefail

echo "=== PC-MIL Benchmark ==="
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: $(hostname)"
echo "GPU: $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null || echo 'N/A')"
date

# Activate environment
eval "$(conda shell.bash hook)"
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

# Run benchmark
cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist

python Benchmarking/xenium_benchmarking/CITEgeist/src/benchmark_pc_mil.py \
    --output-dir Benchmarking/xenium_benchmarking/CITEgeist/output/pc_mil \
    --features-dir Benchmarking/xenium_benchmarking/CITEgeist/output/vit_features \
    --n-epochs 200 \
    --device cuda \
    ${FOLD:+--fold $FOLD}

echo "=== Done ==="
date
```

- [ ] **Step 3: Verify scripts parse correctly**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -c "import ast; ast.parse(open('Benchmarking/xenium_benchmarking/CITEgeist/src/benchmark_pc_mil.py').read()); print('OK')"`
Expected: `OK`

- [ ] **Step 4: Commit**

```bash
git add Benchmarking/xenium_benchmarking/CITEgeist/src/benchmark_pc_mil.py \
       Benchmarking/xenium_benchmarking/CITEgeist/slurm/sbatch_pc_mil.sh
git commit -m "feat: add Xenium PC-MIL benchmark harness with 5-fold CV and SLURM script"
```

---

## Chunk 5: ViT Feature Extraction + Data Pipeline

### Task 5: Pre-extract ViT-S features for Xenium nucleus patches

This task adapts the existing ViT feature extraction from `vit_extractor.py` to process Xenium DAPI+boundary patches into cached 384-dim features. This is a prerequisite for running the benchmark.

**Files:**
- Create: `Benchmarking/xenium_benchmarking/CITEgeist/src/extract_vit_features.py`
- Create: `Benchmarking/xenium_benchmarking/CITEgeist/slurm/sbatch_extract_vit.sh`

- [ ] **Step 1: Write the feature extraction script**

Create `Benchmarking/xenium_benchmarking/CITEgeist/src/extract_vit_features.py`:

```python
#!/usr/bin/env python
"""
Extract ViT-S features from Xenium nucleus patches.

For each pseudo-Visium region:
1. Load pre-extracted patches (DAPI + boundary, 2-channel)
2. Convert to 3-channel: [DAPI, boundary, boundary]
3. Resize to 224x224
4. Run frozen ViT-S (ImageNet pretrained)
5. Save (N, 384) feature tensor + spot-nucleus mapping

Output structure per region:
    output_dir/region_{id}/
        vit_features.npy    # (N_total, 384)
        nucleus_ids.npy     # (N_total,) string array
        spot_mapping.csv    # nucleus_id -> spot_barcode mapping

Usage:
    python extract_vit_features.py \
        --patches-dir output/vae_masked/patches_combined \
        --output-dir output/vit_features \
        --device cuda
"""
import argparse
import logging
import sys
from pathlib import Path

import numpy as np
import torch

REPO_ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO_ROOT))

from CITEgeist.model.vit_extractor import ViTFeatureExtractor

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[logging.StreamHandler(sys.stderr)],
)
logger = logging.getLogger(__name__)


def convert_2ch_to_3ch(patches: np.ndarray) -> np.ndarray:
    """Convert 2-channel (DAPI+boundary) to 3-channel [DAPI, boundary, boundary].

    Args:
        patches: (N, 2, H, W) array

    Returns:
        (N, 3, H, W) array
    """
    dapi = patches[:, 0:1, :, :]     # (N, 1, H, W)
    boundary = patches[:, 1:2, :, :]  # (N, 1, H, W)
    return np.concatenate([dapi, boundary, boundary], axis=1)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--patches-dir", type=str, required=True)
    parser.add_argument("--output-dir", type=str, required=True)
    parser.add_argument("--vit-model", type=str, default="vit_small_patch16_224")
    parser.add_argument("--batch-size", type=int, default=64)
    parser.add_argument("--device", type=str, default="cuda")
    args = parser.parse_args()

    patches_dir = Path(args.patches_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Initialize ViT extractor
    logger.info(f"Loading ViT model: {args.vit_model}")
    vit = ViTFeatureExtractor(
        model_name=args.vit_model,
        pretrained=True,
        device=args.device,
    )
    logger.info(f"ViT embed_dim: {vit.embed_dim}")

    # Process each region
    for region_dir in sorted(patches_dir.iterdir()):
        if not region_dir.is_dir():
            continue

        region_name = region_dir.name
        logger.info(f"Processing region: {region_name}")

        # Load patches (expecting *_patches.npy files)
        patch_files = sorted(region_dir.glob("*_patches.npy"))
        if not patch_files:
            logger.warning(f"No patch files in {region_dir}, skipping")
            continue

        all_features = []
        all_ids = []

        for pf in patch_files:
            patches = np.load(pf)
            spot_id = pf.stem.replace("_patches", "")

            # Load nucleus IDs if available
            id_file = region_dir / f"{spot_id}_nucleus_ids.npy"
            if id_file.exists():
                nuc_ids = np.load(id_file, allow_pickle=True)
            else:
                nuc_ids = np.array([f"{spot_id}_nuc{i}" for i in range(len(patches))])

            if len(patches) == 0:
                continue

            # Convert 2ch -> 3ch if needed
            if patches.ndim == 4 and patches.shape[1] == 2:
                patches = convert_2ch_to_3ch(patches)

            # Extract features
            features = vit.extract_numpy(patches.astype(np.float32), batch_size=args.batch_size)
            all_features.append(features)
            all_ids.extend(nuc_ids)

        if not all_features:
            logger.warning(f"No features extracted for {region_name}")
            continue

        all_features = np.concatenate(all_features, axis=0)
        all_ids = np.array(all_ids)

        # Save
        out_region = output_dir / region_name
        out_region.mkdir(exist_ok=True)
        np.save(out_region / "vit_features.npy", all_features)
        np.save(out_region / "nucleus_ids.npy", all_ids)

        logger.info(f"  Saved {all_features.shape[0]} features ({all_features.shape[1]}-dim)")

    logger.info("Done.")


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Write SLURM script for extraction**

Create `Benchmarking/xenium_benchmarking/CITEgeist/slurm/sbatch_extract_vit.sh`:

```bash
#!/bin/bash
#SBATCH --job-name=extract_vit
#SBATCH --output=Benchmarking/xenium_benchmarking/CITEgeist/slurm/logs/extract_vit_%j.out
#SBATCH --error=Benchmarking/xenium_benchmarking/CITEgeist/slurm/logs/extract_vit_%j.err
#SBATCH --partition=l40s
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

set -euo pipefail

echo "=== ViT Feature Extraction ==="
date

eval "$(conda shell.bash hook)"
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist

python Benchmarking/xenium_benchmarking/CITEgeist/src/extract_vit_features.py \
    --patches-dir Benchmarking/xenium_benchmarking/CITEgeist/output/vae_masked/patches_combined \
    --output-dir Benchmarking/xenium_benchmarking/CITEgeist/output/vit_features \
    --device cuda

echo "=== Done ==="
date
```

- [ ] **Step 3: Verify script parses**

Run: `python -c "import ast; ast.parse(open('Benchmarking/xenium_benchmarking/CITEgeist/src/extract_vit_features.py').read()); print('OK')"`
Expected: `OK`

- [ ] **Step 4: Commit**

```bash
git add Benchmarking/xenium_benchmarking/CITEgeist/src/extract_vit_features.py \
       Benchmarking/xenium_benchmarking/CITEgeist/slurm/sbatch_extract_vit.sh
git commit -m "feat: add ViT-S feature extraction for Xenium nucleus patches"
```

---

## Chunk 6: Integration + End-to-End Validation

### Task 6: Implement `load_region_data` in benchmark harness

This task fills in the `load_region_data` function based on actual data layout. It needs to connect pre-extracted ViT features, Module 3 proportions, protein signals, and ground truth labels.

**Files:**
- Modify: `Benchmarking/xenium_benchmarking/CITEgeist/src/benchmark_pc_mil.py`

- [ ] **Step 1: Inspect actual data layout**

Run:
```bash
ls Benchmarking/xenium_pseudovisium/data_protein_gt/ | head -20
ls Benchmarking/xenium_benchmarking/CITEgeist/output/hybrid_detection_filter/ | head -20
```

Use output to implement `load_region_data()` with correct file paths, column names, and data formats.

- [ ] **Step 2: Implement `load_region_data` based on actual layout**

Replace the `raise NotImplementedError` in `benchmark_pc_mil.py` with the actual loading logic, matching column names and file formats from step 1.

- [ ] **Step 3: Dry-run test (CPU, 1 epoch, 1 fold)**

Run:
```bash
cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist
python Benchmarking/xenium_benchmarking/CITEgeist/src/benchmark_pc_mil.py \
    --output-dir /tmp/pc_mil_test \
    --features-dir Benchmarking/xenium_benchmarking/CITEgeist/output/vit_features \
    --n-epochs 1 --fold 0 --device cpu
```
Expected: Script runs to completion, produces `fold0_results.json`

- [ ] **Step 4: Commit**

```bash
git add Benchmarking/xenium_benchmarking/CITEgeist/src/benchmark_pc_mil.py
git commit -m "feat: implement data loading for PC-MIL Xenium benchmark"
```

### Task 7: Run full benchmark on GPU

- [ ] **Step 1: Extract ViT features (if not already cached)**

```bash
sbatch Benchmarking/xenium_benchmarking/CITEgeist/slurm/sbatch_extract_vit.sh
```

Wait for completion, verify output exists.

- [ ] **Step 2: Submit 5-fold benchmark**

```bash
sbatch Benchmarking/xenium_benchmarking/CITEgeist/slurm/sbatch_pc_mil.sh
```

- [ ] **Step 3: Check results**

After job completes:
```bash
cat Benchmarking/xenium_benchmarking/CITEgeist/output/pc_mil/summary.json
```

Verify:
- `mean_proportion_r >= 0.72` (at minimum matches baseline)
- `mean_single_cell_accuracy >= 0.50` (spec target, above 0.46 constrained Hungarian baseline)
- Per-type accuracy: no type below random baseline (1/K = ~14%)

- [ ] **Step 4: Commit results**

```bash
git add Benchmarking/xenium_benchmarking/CITEgeist/output/pc_mil/summary.json
git commit -m "results: PC-MIL Xenium benchmark 5-fold CV results"
```

---

## Task Dependencies

```
Task 0 (bug fix) ──────────────────────────┐
                                            │
Task 1 (PCMILModel) ──→ Task 2 (training) ──→ Task 3 (inference) ──→ Task 4 (benchmark)
                                                                           │
Task 5 (ViT extraction) ─────────────────────────────────→ Task 6 (integration)
                                                                           │
                                                                     Task 7 (run)
```

Tasks 0 and 5 can run in parallel with Tasks 1-3.
Task 4 depends on Tasks 1-3.
Task 6 depends on Tasks 4 and 5.
Task 7 depends on Task 6.
