# H&E Morphology Pipeline with ViT + MIL Design

**Date:** 2026-03-02
**Status:** Approved
**Author:** Claude + User

## Overview

Extend CITEgeist's morphology-based cell assignment system to support H&E images (3-channel RGB) and improve accuracy via pre-trained Vision Transformer (ViT) features with Multiple Instance Learning (MIL) aggregation.

### Goals

1. **Wire H&E support** into the existing morphology pipeline
2. **Improve accuracy** via ViT + proportion-guided MIL
3. **Compare feature importance** between H&E (3ch) and DAPI+boundary (2ch)
4. **Create pseudo-Visium benchmark** from Visium HD pILC data

## Architecture

```
┌─────────────────────────────────────────────────────────────────────┐
│                    H&E Morphology Pipeline                          │
├─────────────────────────────────────────────────────────────────────┤
│                                                                     │
│  1. DATA PREPARATION                                                │
│     Visium HD H&E → Cellpose segmentation → Pseudo-Visium spots     │
│     (55μm grid)                                                     │
│                                                                     │
│  2. PATCH EXTRACTION                                                │
│     Per nucleus: crop 224×224 RGB patch, normalize                  │
│     Output: spot_{id}_patches.npy (N, 3, 224, 224)                  │
│                                                                     │
│  3. FEATURE EXTRACTION                                              │
│     Branch A: Pre-trained ViT (UNI) → 768-dim embeddings            │
│     Branch B: VAE encoder → 128-dim + morphology 57-dim (baseline)  │
│                                                                     │
│  4. MIL AGGREGATION (proportion-guided)                             │
│     [N nucleus embeddings] → Attention weights → Spot proportions   │
│     Loss: MSE(pred_props, gt_props) + KL regularization             │
│                                                                     │
│  5. SINGLE-CELL ASSIGNMENT                                          │
│     Attention weights → Per-nucleus type probabilities              │
│     + Hungarian assignment with spot counts                         │
│                                                                     │
└─────────────────────────────────────────────────────────────────────┘
```

## Design Decisions

| Decision | Choice | Rationale |
|----------|--------|-----------|
| Patch size | 224×224 | Match ViT input size |
| ViT model | UNI (MahmoodLab) | Trained on 100M+ histopathology patches |
| Spot diameter | 55μm | Standard Visium size |
| Unknown cells | Exclude | Cleaner ground truth |
| Segmentation | Cellpose re-segmentation | Consistency with Xenium workflow |
| Training signal | Proportion MSE | Soft guidance via MIL |

## Data Sources

### pILC Visium HD (H&E benchmark)

| Sample | Cells | Key Types |
|--------|-------|-----------|
| TP08-2202 | 464,680 | Fibroblast (170K), Plasma (45K), Macrophage (37K), Cancer (32K) |
| TP12-880 | 120,577 | Fibroblast (58K), Cancer (16K), Macrophage (3K) |
| TP15-M509 | 440,535 | Fibroblast (175K), Macrophage (58K), Cancer (45K), T_Cell (38K) |

**Total:** ~1M cells with ScType ground truth (excluding Unknown)

**H&E Images:** ~27K × 27K pixels RGB TIFF per sample

**Location:** `/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/pILC_project/`

### Xenium (DAPI+boundary baseline)

Existing benchmark at `Benchmarking/xenium_benchmarking/`

## Component Designs

### 1. Pseudo-Visium Creation

```python
SPOT_DIAMETER_UM = 55
SPOT_SPACING_UM = 100  # center-to-center (standard Visium hex grid)

# Grid generation
spots = create_hex_grid(wsi_bounds, spacing=SPOT_SPACING_UM)

# Cell-to-spot assignment
for cell in cells:
    spot_id = find_nearest_spot(cell.centroid, spots, max_dist=SPOT_DIAMETER_UM/2)

# Ground truth proportions per spot
spot_proportions = cells.groupby('spot_id')['cell_type'].value_counts(normalize=True)
```

**Filtering:**
- Exclude "Unknown" cells
- Minimum 5 cells per spot
- Merge rare types (<1%) into "Other" or exclude

### 2. Cellpose Re-segmentation

```python
from cellpose import models

model = models.Cellpose(model_type='nuclei')

# Process in tiles (WSI too large for single pass)
for tile in tiles:
    masks, flows, styles = model.eval(tile, diameter=30, channels=[0,0])
```

### 3. Feature Extraction

**Branch A: Pre-trained ViT**
```python
from timm import create_model

class ViTFeatureExtractor:
    def __init__(self):
        self.model = create_model('vit_large_patch16_224', pretrained=False, num_classes=0)
        self.model.load_state_dict(torch.load('uni_weights.pth'))
        self.model.eval()

    def extract(self, patches):  # (N, 3, 224, 224)
        with torch.no_grad():
            return self.model(patches)  # (N, 768)
```

**Branch B: Multi-Channel VAE**
```python
class MultiChannelVAE(nn.Module):
    def __init__(self, in_channels=3, latent_dim=128):
        self.encoder = VAEEncoder(in_channels, latent_dim)
        self.decoder = VAEDecoder(in_channels, latent_dim)
```

**H&E Morphology Features (~65 dims):**
- Shape features (12 dims) - existing
- Stain deconvolution: hematoxylin_mean, eosin_mean, he_ratio
- Per-channel stats: R/G/B mean, std (6 dims)
- Texture features

### 4. MIL Aggregation

```python
class ProportionGuidedMIL(nn.Module):
    def __init__(self, input_dim=768, n_cell_types=9, hidden_dim=256):
        # Gated attention mechanism
        self.attention_V = nn.Sequential(nn.Linear(input_dim, hidden_dim), nn.Tanh())
        self.attention_U = nn.Sequential(nn.Linear(input_dim, hidden_dim), nn.Sigmoid())
        self.attention_W = nn.Linear(hidden_dim, n_cell_types)

        # Proportion prediction head
        self.classifier = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, n_cell_types),
            nn.Softmax(dim=-1)
        )

    def forward(self, embeddings):  # (N_nuclei, 768)
        A_V = self.attention_V(embeddings)
        A_U = self.attention_U(embeddings)
        A = self.attention_W(A_V * A_U)
        A = F.softmax(A, dim=0)  # Normalize across nuclei

        spot_repr = torch.mm(A.T, embeddings)
        proportions = self.classifier(spot_repr.mean(dim=0))

        return proportions, A
```

**Training:**
```python
loss_prop = F.mse_loss(pred_proportions, gt_proportions)
loss_entropy = -torch.mean(torch.sum(attention * torch.log(attention + 1e-8), dim=0))
loss = loss_prop + 0.1 * loss_entropy
```

**Strategy:**
- Stage 1: Train MIL head with frozen ViT
- Stage 2: Optional fine-tuning of ViT last layers
- Augmentation: Random rotation, flip, color jitter

### 5. Single-Cell Assignment

```python
from scipy.optimize import linear_sum_assignment

def assign_nuclei_to_types(probs, spot_proportions, n_nuclei):
    # Convert proportions to integer counts
    counts = np.round(spot_proportions * n_nuclei).astype(int)

    # Build cost matrix
    slots = []
    for type_idx, count in enumerate(counts):
        slots.extend([type_idx] * count)

    cost_matrix = np.zeros((n_nuclei, n_nuclei))
    for i in range(n_nuclei):
        for j, type_idx in enumerate(slots):
            cost_matrix[i, j] = -np.log(probs[i, type_idx] + 1e-8)

    row_ind, col_ind = linear_sum_assignment(cost_matrix)
    return {i: slots[col_ind[i]] for i in range(n_nuclei)}
```

### 6. Feature Importance Comparison

**Channel Ablation:**
```python
def channel_ablation(model, patches, gt):
    results = {}
    for drop in ['R', 'G', 'B', 'H', 'E']:
        ablated = zero_channel(patches, drop)
        acc = evaluate(model, ablated, gt)
        results[f'drop_{drop}'] = acc
    return results
```

**Attention-Morphology Correlation:**
```python
def analyze_attention_by_morphology(attention, morphology_features):
    correlations = {}
    for cell_type in range(n_types):
        for feat_name, feat_values in morphology_features.items():
            r, p = spearmanr(attention[:, cell_type], feat_values)
            correlations[(cell_type, feat_name)] = {'r': r, 'p': p}
    return pd.DataFrame(correlations)
```

## Output Structure

```
Benchmarking/visiumhd_benchmarking/
├── data/
│   ├── TP08-2202/
│   │   ├── pseudo_visium_spots.csv
│   │   ├── cell_to_spot_mapping.csv
│   │   ├── ground_truth_proportions.parquet
│   │   └── cellpose_masks.npy
│   ├── TP12-880/
│   └── TP15-M509/
├── patches/
│   └── {sample}/spot_{id}_patches.npy
├── models/
│   ├── mil_head.pt
│   └── vae_baseline.pt
├── results/
│   ├── feature_comparison/
│   │   ├── channel_ablation.csv
│   │   ├── attention_morphology_corr.csv
│   │   └── summary_report.md
│   └── evaluation/
│       ├── proportion_correlation.csv
│       ├── single_cell_accuracy.csv
│       └── per_celltype_f1.csv
└── src/
    ├── create_pseudo_visium.py
    ├── run_cellpose.py
    ├── extract_patches_he.py
    ├── vit_feature_extractor.py
    ├── train_mil.py
    ├── evaluate.py
    └── compare_features.py
```

## Evaluation Metrics

| Metric | Description |
|--------|-------------|
| **Proportion Pearson r** | Correlation between predicted and GT spot proportions |
| **Global accuracy** | Single-cell accuracy without constraints (argmax) |
| **Constrained accuracy (oracle)** | Hungarian with GT proportions |
| **Constrained accuracy (predicted)** | Hungarian with predicted proportions |
| **Per-type F1** | F1 score for each cell type |

## Expected Outcomes

| Method | Estimated Global Acc | Estimated Constrained Acc |
|--------|---------------------|---------------------------|
| Current (VAE + XGBoost, 2ch) | 33.5% | 46% |
| VAE baseline (3ch H&E) | 35-40% | 50-55% |
| **ViT + MIL (3ch H&E)** | **45-55%** | **60-70%** |

## Key Questions to Answer

1. Does H&E's cytoplasm info help distinguish immune subtypes (T vs B vs NK)?
2. Does DAPI's nuclear clarity help for cancer vs stromal discrimination?
3. Which cell types are hardest in each modality?
4. Are the learned features complementary or redundant?

## Dependencies

**New packages:**
- `timm` - PyTorch Image Models (for ViT)
- UNI weights download (requires Hugging Face access)

**Existing:**
- `cellpose` - Already in CITEgeist_env
- `torch`, `scipy`, `numpy`, `pandas` - Already available

## References

- UNI: Chen et al. "A General-Purpose Self-Supervised Model for Computational Pathology" (2024)
- MIL for histopathology: Ilse et al. "Attention-based Deep Multiple Instance Learning" (2018)
- Stain deconvolution: Ruifrok & Johnston (2001)
