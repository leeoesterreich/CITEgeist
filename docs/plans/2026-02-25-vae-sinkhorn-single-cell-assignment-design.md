# VAE + Sinkhorn OT Single-Cell Assignment Design

**Date:** 2026-02-25
**Status:** Approved
**Author:** Brainstorming session

## Problem Statement

CITEgeist produces accurate spot-level cell type proportions (r=0.726) but cannot determine which specific nucleus within a spot belongs to which cell type. The current morphology-based approach only achieves ~48% accuracy vs ~47% random baseline (+1% improvement), because nuclear/cell morphology features don't discriminate cell types in this tissue.

**Goal:** Develop a method to assign individual nuclei to cell types that significantly exceeds random assignment, using only Visium-level inputs (spot-aggregated molecular data + imaging).

**Constraints:**
- Input: Pseudo-Visium level data (spot-aggregated protein + GEX + DAPI imaging)
- Cannot use per-cell molecular data (not available in real Visium)
- Must work with cancer tissue (violates spatial smoothness assumptions)

**Target accuracy:** 60-70%+ (vs 14% random baseline for 7 cell types)

## Approach

Two-stage deep learning approach:
1. **Stage 1:** Train VAE on nucleus patches to learn compressed visual representations
2. **Stage 2:** Learn cell-type-specific projection heads and prototypes using spot-level proportion constraints via Sinkhorn optimal transport

### Why Previous Morphology Failed

Hand-crafted morphology features (area, circularity, eccentricity, N:C ratio, etc.) only captured ~5% of cell type signal. The VAE approach can learn richer features from raw pixels including texture, chromatin patterns, and subtle visual cues that aren't captured by simple geometric metrics.

### Why Sinkhorn OT

The key challenge is we don't have per-nucleus labels, only spot-level proportion constraints. Sinkhorn optimal transport:
- Finds assignments that respect proportion constraints
- Is fully differentiable (gradients flow to projection heads and prototypes)
- Naturally avoids degenerate solutions (predicting spot-average for every nucleus)

## Architecture

### Stage 1: VAE Pretraining (Unsupervised)

```
Input: Expanded nucleus patch (96×96 pixels)
       - Nucleus bounding box + 75% padding
       - Channels: DAPI + boundary markers (3-4 channels)
                    ↓
            ┌───────────────┐
            │  Encoder CNN  │  (ResNet-18 style)
            │  Conv blocks  │
            │  → FC layers  │
            └───────────────┘
                    ↓
            μ, log(σ²) ∈ ℝ^128
                    ↓
            Reparameterization: z = μ + σ × ε
                    ↓
            ┌───────────────┐
            │  Decoder CNN  │  (Transposed convs)
            └───────────────┘
                    ↓
            Reconstructed patch

Loss: L_VAE = L_recon + β × KL(q(z|x) || p(z))
      L_recon = MSE(x, x̂)
      β = 0.5 (tunable)
```

**Training:**
- Batch: Random sample of patches across all spots/regions
- No labels needed - pure reconstruction
- Train until reconstruction loss plateaus
- Output: Frozen encoder weights

### Stage 2: Prototype Learning (Weakly Supervised)

```
            z ∈ ℝ^128 (from frozen encoder)
                    ↓
    ┌─────────────────────────────────────┐
    │  K Projection Heads                 │
    │  head_k: MLP(128 → 64 → 32)         │
    │  Each head learns which features    │
    │  matter for its cell type           │
    └─────────────────────────────────────┘
                    ↓
            projected_k ∈ ℝ^32 for each type k
                    ↓
    ┌─────────────────────────────────────┐
    │  K Prototypes                       │
    │  pₖ ∈ ℝ^32 (learnable parameters)   │
    └─────────────────────────────────────┘
                    ↓
            Distance matrix:
            D[i,k] = ||head_k(zᵢ) - pₖ||²
                    ↓
    ┌─────────────────────────────────────┐
    │  Sinkhorn OT (per spot)             │
    │  Row marginals: uniform (1/N each)  │
    │  Col marginals: spot proportions    │
    │  Temperature: 0.1                   │
    └─────────────────────────────────────┘
                    ↓
            Transport plan P ∈ ℝ^(N×K)
                    ↓
            Loss = Σᵢⱼ P[i,j] × D[i,j]
```

**Training:**
- Process spots as units (variable N nuclei per spot)
- Supervision: CITEgeist continuous proportions (spot-level)
- Backprop updates: projection heads + prototypes
- Encoder remains frozen

### Sinkhorn Algorithm

```python
def sinkhorn(cost, row_marginal, col_marginal, temperature=0.1, n_iters=50):
    """
    Differentiable optimal transport via Sinkhorn iterations.

    Args:
        cost: (N, K) distance matrix
        row_marginal: (N,) uniform distribution
        col_marginal: (K,) spot proportions
        temperature: lower = sharper assignments
        n_iters: iteration count

    Returns:
        P: (N, K) transport plan (soft assignments)
    """
    # Initialize with scaled cost
    K_matrix = torch.exp(-cost / temperature)

    # Sinkhorn iterations (alternating row/col normalization)
    for _ in range(n_iters):
        K_matrix = K_matrix / K_matrix.sum(dim=1, keepdim=True) * row_marginal.unsqueeze(1)
        K_matrix = K_matrix / K_matrix.sum(dim=0, keepdim=True) * col_marginal.unsqueeze(0)

    return K_matrix
```

### Inference

```python
def assign_nuclei(patches, spot_proportions, encoder, heads, prototypes):
    # 1. Encode
    z = encoder(patches)  # (N, 128)

    # 2. Project through each head
    projected = [head_k(z) for head_k in heads]  # K × (N, 32)

    # 3. Compute distances to prototypes
    distances = torch.stack([
        torch.norm(projected[k] - prototypes[k], dim=1)
        for k in range(K)
    ], dim=1)  # (N, K)

    # 4. Sinkhorn assignment
    transport_plan = sinkhorn(distances, uniform(N), spot_proportions, temp=0.05)

    # 5. Hard assignment + confidence
    assignments = transport_plan.argmax(dim=1)
    confidence = transport_plan.max(dim=1).values

    return assignments, confidence
```

## Data Pipeline

### Patch Extraction

```python
def extract_patch(image, nucleus_bbox, expansion=0.75, output_size=96):
    """
    Extract expanded patch around nucleus.

    Args:
        image: (C, H, W) multi-channel image (DAPI + boundaries)
        nucleus_bbox: (x_min, y_min, x_max, y_max) from Cellpose
        expansion: fraction to expand bbox in each direction
        output_size: final patch size after resize
    """
    # Expand bounding box
    w, h = x_max - x_min, y_max - y_min
    x_min -= w * expansion
    x_max += w * expansion
    y_min -= h * expansion
    y_max += h * expansion

    # Crop and resize
    patch = image[:, y_min:y_max, x_min:x_max]
    patch = resize(patch, (output_size, output_size))

    # Normalize per channel
    patch = (patch - patch.mean(dim=(1,2), keepdim=True)) / patch.std(dim=(1,2), keepdim=True)

    return patch
```

### Data Split

| Split | Regions | Nuclei (approx) | Purpose |
|-------|---------|-----------------|---------|
| Train | 3-4 | ~70-90K | Model training |
| Val | 1 | ~20-25K | Hyperparameter tuning |
| Test | 1 | ~20-25K | Final evaluation |

## Hyperparameters

| Parameter | Value | Notes |
|-----------|-------|-------|
| Patch size | 96×96 | After expansion and resize |
| Bbox expansion | 75% | Captures cytoplasm and context |
| Image channels | 3-4 | DAPI + boundary markers |
| VAE latent dim | 128 | Shared representation |
| VAE β | 0.5 | KL weight, tune to avoid collapse |
| Projection dim | 32 | Per-type projected space |
| Projection head | MLP(128→64→32) | 2-layer with ReLU |
| Sinkhorn temp | 0.1 (train), 0.05 (inference) | Lower = sharper |
| Sinkhorn iters | 50 | Usually sufficient |
| Prototype init | K-means | On projected latents after head init |
| Learning rate | 1e-4 (VAE), 1e-3 (Stage 2) | Adam optimizer |
| Batch size | 64 patches (VAE), 16 spots (Stage 2) | |

## Evaluation

### Ground Truth

Two Xenium-derived ground truth annotations:

| GT Type | Source | Cell Types |
|---------|--------|------------|
| Protein-gated | Xenium antibody panel | 7 types + Unknown |
| RNA-clustered | Xenium transcript clustering | 6 types |

### Metrics

1. **Per-nucleus accuracy** - Overall and per-class
2. **Precision/Recall/F1** - Per cell type
3. **Confusion matrix** - Which types get confused
4. **Within-spot accuracy** - Distribution across spots
5. **Confidence calibration** - Are high-confidence predictions more accurate?

### Baselines

| Method | Expected Accuracy |
|--------|-------------------|
| Random assignment | ~14% (1/7 types) |
| Morphology features | ~48% |
| VAE + K-means (no Sinkhorn) | TBD |
| **This method** | Target 60-70%+ |

### Ablations

- Latent dimension: 64, 128, 256
- Projection dimension: 16, 32, 64
- Patch expansion: 50%, 75%, 100%
- Sinkhorn temperature: 0.05, 0.1, 0.2
- With/without VAE pretraining (train end-to-end)
- Prototype init: random vs K-means

## Visualizations

1. **Latent space UMAP** - Colored by predicted type and GT type
2. **Prototype locations** - Visualize K prototypes in projected space
3. **Per-spot accuracy histogram** - Distribution of within-spot accuracy
4. **Attention analysis** - Which latent features do heads use?
5. **Failure cases** - Examine misclassified nuclei

## Implementation Notes

### File Structure

```
CITEgeist/model/
├── vae_encoder.py          # VAE architecture
├── projection_heads.py     # K projection heads
├── sinkhorn.py            # Sinkhorn OT implementation
├── prototype_learning.py  # Stage 2 training loop
├── patch_extraction.py    # Data pipeline
└── single_cell_vae_assignment.py  # Main interface
```

### Dependencies

- PyTorch (existing)
- torchvision (for transforms)
- No new major dependencies

### Integration with Module 3b

This replaces the current morphology-based assignment in `module3b_nucleus_assignment.py`. The interface remains the same:

```python
# Current
assignments = assign_nuclei_random(nuclei, spot_proportions)

# New
assignments = assign_nuclei_vae_sinkhorn(
    nuclei_patches,
    spot_proportions,
    encoder, heads, prototypes
)
```

## Risks and Mitigations

| Risk | Mitigation |
|------|------------|
| VAE doesn't learn discriminative features | Ablate with supervised pretraining on dominated spots |
| Sinkhorn doesn't converge | Increase iterations, tune temperature |
| Overfitting to training regions | Cross-validate across regions |
| Cancer heterogeneity breaks prototypes | Per-patient or per-region prototype adaptation |
| Computational cost | Batch processing, GPU acceleration |

## Success Criteria

- [ ] VAE reconstruction visually reasonable
- [ ] Latent space shows structure (not random)
- [ ] Prototypes separate in projected space
- [ ] Test accuracy > 55% (meaningful improvement)
- [ ] Test accuracy > 60% (good)
- [ ] Test accuracy > 70% (excellent)

## Future Extensions

1. **H&E integration** - Add H&E channels when available for real Visium
2. **Multi-resolution patches** - Capture both local and context features
3. **Attention over nuclei** - Learn which nuclei are most informative per spot
4. **Domain adaptation** - Transfer to new tissue types
5. **Uncertainty quantification** - Bayesian projection heads
