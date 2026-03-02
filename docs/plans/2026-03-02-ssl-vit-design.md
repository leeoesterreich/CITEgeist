# Self-Supervised ViT for Nucleus Morphology Embeddings

**Date:** 2026-03-02
**Status:** Approved
**Author:** Claude + Alex

## Overview

Replace the current VAE encoder (128-dim) with self-supervised Vision Transformer (ViT) embeddings trained via MAE and DINO. This is a drop-in replacement that keeps the existing morphology features (57-dim) and XGBoost classifier unchanged.

## Problem Statement

The current VAE embeddings have poor cell type discriminability:
- Silhouette score: **-0.088** (negative = worse than random)
- VAE-only XGBoost accuracy: ~25%
- Combined (VAE + morphology): **33.5%** (only +3.7% over morphology-only)

Self-supervised ViT methods (MAE, DINO) learn semantically meaningful representations and should produce embeddings with better cell type separation.

## Success Criteria

1. **Beat current accuracy:** Combined ViT + morphology > 33.5%
2. **Positive silhouette:** ViT embeddings have silhouette > 0 (current: -0.088)
3. **Target >40% accuracy:** Significant improvement to justify added complexity

## Design Decisions

| Decision | Choice | Rationale |
|----------|--------|-----------|
| Approach | Drop-in VAE replacement | Keep proven morphology features, minimize changes |
| Methods | Both MAE and DINO | Empirical comparison on this domain |
| Weights | Train from scratch | 2-channel DAPI+boundary differs from RGB ImageNet |
| Model size | ViT-Small (22M params) | Good capacity for 96x96 patches |
| Compute | L40S single GPU | Available on CRC, 48GB sufficient |

## Architecture

### ViT-Small Backbone

```
Input: (B, 2, 96, 96) - DAPI + boundary channels
Patch size: 16x16
Num patches: 36 (6x6 grid)

Encoder:
  - Patch embedding: Conv2d(2, 384, kernel=16, stride=16)
  - Position embedding: learnable (1 + 36 tokens)
  - Transformer: 12 layers, 6 heads, dim=384, mlp_ratio=4

Output: [CLS] token (384-dim) for downstream tasks
```

### MAE Head

```
Mask ratio: 0.75 (27 of 36 patches masked)
Decoder:
  - embed_dim: 192
  - depth: 4 layers
  - num_heads: 3

Loss: MSE reconstruction of masked patch pixels
```

### DINO Head

```
Projection head:
  - hidden_dim: 2048
  - bottleneck_dim: 256
  - out_dim: 4096

Teacher: EMA of student (momentum 0.996 → 1.0)
Centering: EMA of teacher outputs (prevents collapse)

Loss: Cross-entropy(student_softmax, teacher_softmax)
```

## Data Pipeline

### Input Data

- Source: `output/vae_masked/patches_combined/region_{0-4}/`
- Shape: (N, 2, 96, 96) per patch
- Total: ~50,000+ patches across 5 Xenium regions

### Augmentations

**MAE (minimal):**
- Random horizontal/vertical flip
- Random 90° rotation
- Per-channel normalization

**DINO (multi-crop):**
- 2 global crops (scale 0.4-1.0, size 96x96)
- 6 local crops (scale 0.05-0.4, size 48x48)
- Random flip, rotation, Gaussian blur, intensity jitter

## Training Configuration

### MAE

```yaml
epochs: 200
optimizer: AdamW
  lr: 1.5e-4
  weight_decay: 0.05
  betas: [0.9, 0.95]
scheduler: cosine with 10-epoch warmup
batch_size: 256
mixed_precision: fp16
estimated_time: 4-6 hours (L40S)
```

### DINO

```yaml
epochs: 300
optimizer: AdamW
  lr: 5e-4
  weight_decay: 0.04 → 0.4 (cosine)
  betas: [0.9, 0.999]
teacher_momentum: 0.996 → 1.0 (cosine)
teacher_temp: 0.04
student_temp: 0.1
scheduler: cosine with 10-epoch warmup
batch_size: 128 (multi-crop uses more memory)
mixed_precision: fp16
estimated_time: 8-12 hours (L40S)
```

## Integration

### Embedding Extraction

```python
# After training
encoder = load_trained_encoder("mae_final.pt")  # or dino_final.pt
encoder.eval()

with torch.no_grad():
    embeddings = encoder(patches)  # (N, 384)
```

### XGBoost Pipeline (unchanged structure)

```python
# Old
vae_embeddings = extract_vae_embeddings(...)     # (N, 128)
morph_features = extract_morphology_features(...)  # (N, 57)
combined = np.hstack([vae_embeddings, morph_features])  # (N, 185)

# New (drop-in replacement)
vit_embeddings = extract_vit_embeddings(...)     # (N, 384)
morph_features = extract_morphology_features(...)  # (N, 57)
combined = np.hstack([vit_embeddings, morph_features])  # (N, 441)
```

## Evaluation Plan

1. **Embedding quality:**
   - Silhouette score by cell type
   - t-SNE/UMAP visualization

2. **Classification:**
   - XGBoost accuracy (overall and per-class)
   - Ablation: ViT-only vs Morph-only vs Combined

3. **Comparison:**
   - MAE vs DINO vs VAE baseline
   - Feature importance analysis

## File Structure

### New Model Files

```
CITEgeist/model/
├── vit_encoder.py      # ViT-Small backbone
├── mae.py              # MAE pretext task
├── dino.py             # DINO pretext task
└── ssl_utils.py        # Augmentations, multi-crop
```

### New Training/Eval Scripts

```
Benchmarking/xenium_benchmarking/CITEgeist/src/
├── train_mae.py
├── train_dino.py
├── extract_ssl_embeddings.py
└── evaluate_ssl_embeddings.py
```

### New SLURM Scripts

```
Benchmarking/xenium_benchmarking/CITEgeist/slurm/
├── sbatch_train_mae.sh
├── sbatch_train_dino.sh
└── sbatch_evaluate_ssl.sh
```

### Output Directory

```
Benchmarking/xenium_benchmarking/CITEgeist/output/
├── mae_ssl/
│   ├── mae_final.pt
│   └── training_history.json
├── dino_ssl/
│   ├── dino_final.pt
│   └── training_history.json
└── ssl_evaluation/
    ├── silhouette_scores.json
    ├── comparison_results.json
    └── tsne_plots/
```

## Implementation Order

1. **Phase 1: Infrastructure**
   - `vit_encoder.py` - ViT-Small backbone
   - `ssl_utils.py` - augmentations and data loading

2. **Phase 2: MAE**
   - `mae.py` - MAE head and training logic
   - `train_mae.py` - training script
   - `sbatch_train_mae.sh` - SLURM submission

3. **Phase 3: DINO**
   - `dino.py` - DINO head with teacher-student
   - `train_dino.py` - training script
   - `sbatch_train_dino.sh` - SLURM submission

4. **Phase 4: Evaluation**
   - `extract_ssl_embeddings.py`
   - `evaluate_ssl_embeddings.py`
   - Run comparison against VAE baseline

## References

- MAE: He et al., "Masked Autoencoders Are Scalable Vision Learners" (2021)
- DINO: Caron et al., "Emerging Properties in Self-Supervised Vision Transformers" (2021)
- ViT: Dosovitskiy et al., "An Image is Worth 16x16 Words" (2020)
