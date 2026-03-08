# Unified Single-Cell Assignment Mode (Module 3b)

**Date:** 2026-03-08
**Status:** Approved
**Approach:** B (MIL-First)

## Summary

Integrate DAPI SimCLR (Xenium) and H&E ViT (Visium) morphology models with a shared Gated Attention MIL head and proportion-weighted Hungarian assignment into a unified Module 3b for CITEgeist. Archive discrete IQP mode (poor performance) and CNN VAE / DINO / MAE (SimCLR is best).

## Pipeline

```
Module 3 (Hybrid Protein QP)
    ↓
Continuous Proportions + Cellpose Nuclei Counts
    ↓
Discretize (low proportions → 0)
    ↓
Integer Cell Counts per Type
    ↓
Module 3b (Single-Cell Assignment)
    ├── Backbone: SimCLR ViT-S (DAPI) or ImageNet+SSL ViT-S (H&E)
    ├── MIL Head: Gated Attention (384 → 256 → K)
    ├── Assignment: Proportion-Weighted Hungarian
    └── Output: Nucleus → Cell Type labels
    ↓
GEX Pass 2 (per-cell deconvolution)
```

## Architecture

### Backbone Encoders (modality-specific, both output 384-dim)

| Modality | Backbone | Input | SSL Method | Pre-training |
|----------|----------|-------|------------|--------------|
| DAPI (Xenium) | ViT-Small | 2ch, 96×96 | SimCLR | From scratch on target dataset |
| H&E (Visium) | ViT-Small | 3ch, 224×224 | SimCLR/VICReg fine-tune | ImageNet init via timm |

Shared interface:
```python
class MorphologyBackbone(ABC):
    def extract(self, patches: np.ndarray) -> np.ndarray:
        """(N, C, H, W) → (N, 384) embeddings"""

class DAPIBackbone(MorphologyBackbone):   # ViT-S SimCLR, D=384
class HEBackbone(MorphologyBackbone):     # ViT-S ImageNet+SSL, D=384
```

Training is two-phase per dataset:
1. **Phase 1 (unsupervised):** Train/fine-tune backbone on all nucleus patches
2. **Phase 2 (proportion-supervised):** Freeze backbone, train MIL head using Module 3 proportions

### Gated Attention MIL (shared, modality-agnostic)

- Input: 384-dim nucleus embeddings (from either backbone)
- Architecture: V = tanh(W_v @ x), U = sigmoid(W_u @ x), gate = V * U
- Logits: W_c @ gate → attention = softmax(logits, dim=1)
- Spot proportions: mean(attention, dim=0)
- Training loss: MSE + KL against Module 3 proportions
- Output at inference: attention matrix (N × K) — per-nucleus type probabilities

### Proportion-Weighted Hungarian Assignment

Cost matrix combines morphology likelihood and protein-derived proportion prior:

```
cost[i,k] = -log(attention[i,k]) - λ * log(proportion[k])
```

| Term | Signal | Source |
|------|--------|--------|
| `-log(attention[i,k])` | Morphology: how much does nucleus i look like type k? | MIL |
| `-λ * log(proportion[k])` | Protein: how confident is Module 3 about type k? | Module 3 |

- Integer counts (from discretized proportions) enforce hard constraints
- `λ` controls proportion prior weight (tunable via grid search)
- Solved with `scipy.optimize.linear_sum_assignment`

## File Organization

### New files (`CITEgeist/model/`)

| File | Purpose |
|------|---------|
| `morphology_backbone.py` | ABC + `DAPIBackbone` + `HEBackbone` wrappers |
| `single_cell_mil.py` | Gated attention MIL (shared, 384→256→K) |
| `single_cell_assignment.py` | Proportion-weighted Hungarian engine |

### Modified files

| File | Change |
|------|--------|
| `citegeist_model.py` | Add `run_single_cell_assignment()` (Module 3b entry point) |
| `__init__.py` | Export new classes |

### Reused as-is

| File | Role |
|------|------|
| `vit_encoder.py` | ViT-Small backbone (DAPI SimCLR) |
| `vit_extractor.py` | ViT-Small + timm (H&E ImageNet init) |
| `simclr.py` | SimCLR training method |
| `proportion_mil.py` | Reference (logic migrates to `single_cell_mil.py`) |

### Archived (deprecated)

| File | Reason |
|------|--------|
| `masked_iqp.py` | Discrete IQP mode superseded by Module 3b |
| `vae.py` / `train_vae.py` | CNN VAE superseded by ViT-based SSL |
| `dino.py` | Training collapsed; SimCLR preferred |
| `mae.py` | SimCLR outperforms in constrained setting |

## CitegeistModel API

```python
# Module 3 (existing)
model.run_cell_proportion_model()  # hybrid protein QP → proportions

# Module 3b (new)
assignments_df = model.run_single_cell_assignment(
    patches_dir="path/to/nucleus_patches/",
    nuclei_counts=nuclei_series,          # from Cellpose
    modality="he",                        # or "dapi"
    backbone_checkpoint="path/to/ssl.pt", # pre-trained backbone
    mil_checkpoint=None,                  # None = train using Module 3 props
    lambda_prior=1.0,                     # proportion weighting
)
# Returns: DataFrame with nucleus_id → cell_type

# GEX Pass 2 (uses assignments)
model.run_cell_expression_pass1(
    use_single_cell=True,
    cell_assignments=assignments_df
)
```

## Xenium Benchmark Strategy

### Pipeline per region
1. Module 3 (Hybrid) — existing results (r=0.7428)
2. Cellpose segmentation — existing nuclei counts
3. SimCLR backbone — existing checkpoint
4. MIL training (NEW) — train on region's spots with Module 3 proportions
5. Proportion-weighted Hungarian (NEW) — nucleus labels
6. GEX Pass 2 (NEW) — per-cell deconvolution

### Evaluation metrics

| Metric | Measures | Baseline |
|--------|----------|----------|
| Spot-level proportion r | Assignment re-aggregated to proportions vs GT | 0.7428 |
| Single-cell accuracy | % nuclei correctly typed vs protein-gated GT | 57.5% (SimCLR constrained) |
| GEX Pearson r | Per-cell deconvolved expression quality | 0.3648 |

### Execution
- 5-region leave-one-out cross-validation
- MIL trained per-region; SimCLR backbone shared
- GPU partition for MIL training (~10-30 min/region)
- CPU partition for Hungarian + GEX deconvolution
- SLURM array job across 5 regions

## Key Hypothesis

MIL-learned attention weights provide better per-nucleus type likelihoods than Gaussian-fitted morphology features, leading to higher single-cell accuracy, which in turn improves GEX deconvolution quality.

## SSL Method Selection Rationale

| Method | Constrained Accuracy | Status |
|--------|---------------------|--------|
| **SimCLR** | **57.50%** | **Selected** |
| MAE | 55.79% | Archived |
| DINO | N/A (collapsed) | Archived |
| CNN VAE | ~46% (raw morph) | Archived |
