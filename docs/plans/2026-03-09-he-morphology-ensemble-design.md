# H&E Morphology Ensemble Pipeline for Patient Cohort

**Date**: 2026-03-09
**Status**: Approved
**Goal**: Apply H&E morphology-based single-cell assignment to 12 patient samples, ensemble with protein proportions, and validate qualitatively (no ground truth).

## Motivation

CITEgeist Module 3 produces spot-level cell type proportions from protein signal. Adding H&E morphology (via MIL) provides a complementary modality that can:
1. Assign individual nuclei to cell types (sub-spot resolution)
2. Improve proportion estimates via confidence-weighted ensemble
3. Feed improved proportions into GEX deconvolution (Pass 2)

A previous version had a known failure mode: assigning collagen/stroma regions to T cells. This pipeline must avoid such obvious errors.

## Pipeline Overview

```
STAGE 1: Segmentation
  H&E (cytassist_image.tiff) → Cellpose → nuclei mask → nuclei-to-spot assignment

STAGE 2: Feature Extraction
  224×224 patches per nucleus → ViT-Small (ImageNet) → 384-dim embeddings

STAGE 3: MIL Training
  Supervised by Module 3 GLOBAL proportions (finetuning OFF)
  ProportionGuidedMIL, ~100 epochs, one model per sample

STAGE 4: Ensemble Proportions
  Module 3 protein proportions (weighted by 1/recon_error)
  + MIL predicted proportions (weighted by 1/entropy)
  → confidence-weighted ensemble per spot

STAGE 5a: GEX Pass 2
  Ensemble proportions → deconvolved gene expression

STAGE 5b: Cell Assignment
  MIL attention weights → cost matrix → Hungarian assignment

STAGE 6: Visualization & Sanity Checks
  • Spatial overlay (cell types on H&E)
  • Attention heatmaps per cell type
  • Proportion correlation (ensemble vs protein-only)
```

## Design Decisions

| Decision | Choice | Rationale |
|----------|--------|-----------|
| Supervision | Module 3 global proportions | Finetuned proportions have inconsistent performance |
| Training scope | One model per sample | Avoids batch effects, enables per-sample debugging |
| ViT backbone | vit_small_patch16_224 (ImageNet) | UNI foundation model showed no improvement in benchmarks |
| Epochs | 100 (default) | Training plateaus here; more causes overfitting |
| Finetuning | OFF by default | User-reported inconsistent performance |
| Ensemble weighting | Per-spot confidence | Protein dominates where protein is reliable; morphology contributes where it's confident |

## Confidence-Weighted Ensemble (Stage 4)

```python
# Per-spot weights
w_protein = 1 / (recon_error + eps)   # Module 3 reconstruction error
w_mil = 1 / (entropy + eps)           # MIL prediction entropy (Shannon)

# Normalize
w_protein_norm = w_protein / (w_protein + w_mil)
w_mil_norm = w_mil / (w_protein + w_mil)

# Ensemble
prop_ensemble = w_protein_norm * prop_module3 + w_mil_norm * prop_mil
```

- **Module 3 reconstruction error**: How well proportions explain observed protein signal (already computed during optimization)
- **MIL entropy**: Shannon entropy of predicted proportion vector. Uniform = uncertain, peaked = confident.

## Files to Create

| File | Purpose |
|------|---------|
| `examples/run_morphology_assignment.py` | Main patient pipeline script (all 6 stages) |
| `examples/slurm/sbatch_morphology.sh` | SLURM array job for 12 samples |
| `model/ensemble_proportions.py` | Confidence-weighted ensemble logic |

## Existing Code to Reuse

| Component | Source | Notes |
|-----------|--------|-------|
| Cellpose segmentation | `model/segmentation.py` | `run_cellpose_nuclei_segmentation()` |
| Patch extraction | `model/patch_extraction.py` | `extract_patch()` with global normalization |
| ViT feature extraction | `model/vit_extractor.py` | `ViTFeatureExtractor` (384-dim) |
| ProportionGuidedMIL | `model/proportion_mil.py` | Gated attention + proportion loss |
| Hungarian assignment | `model/hungarian_assignment.py` | Count-constrained optimal matching |
| GEX Pass 2 | `model/gurobi_impl.py` | `optimize_gene_expression()` |

## Parameters

| Parameter | Default | Notes |
|-----------|---------|-------|
| `epochs` | 100 | Plateau point from Visium HD benchmarks |
| `vit_model` | `vit_small_patch16_224` | ImageNet pretrained, 384-dim output |
| `finetune_proportions` | False | Off by default |
| `cellpose_model` | `nuclei` | For H&E nucleus detection |
| `min_nuclei_per_spot` | 1 | Skip spots with no detected nuclei |
| `patch_size` | 224 | Standard ViT input size |
| `batch_size` | 256 | For ViT forward pass |
| `learning_rate` | 1e-3 | From benchmark tuning |

## Output Structure

```
output/morphology_assignment/
├── HCC22-088-P1-S1/
│   ├── nuclei_mask.npy              # Cellpose segmentation
│   ├── embeddings/                  # Per-spot ViT embeddings
│   ├── mil_checkpoint.pt            # Trained MIL model
│   ├── cell_assignments.csv         # nucleus_id, spot, assigned_type
│   ├── ensemble_proportions.csv     # Blended protein+MIL proportions
│   ├── ensemble_gex_pass2.parquet   # GEX deconvolved with ensemble props
│   ├── confidence_weights.csv       # Per-spot w_protein, w_mil
│   ├── spatial_overlay.png          # Cell types overlaid on H&E
│   └── attention_heatmaps/          # Per-cell-type attention maps
├── HCC22-088-P1-S2/
│   └── ...
└── summary/
    ├── proportion_correlation.csv
    └── ensemble_vs_protein_comparison.csv
```

## Qualitative Validation (No Ground Truth)

### Spatial Overlay Checks
- T cells should cluster in immune-infiltrated regions, NOT stroma/collagen
- Fibroblasts should be in stromal/connective tissue regions
- Cancer cells should align with tumor nests visible in H&E
- Endothelial cells should trace vascular structures

### Attention Heatmap Checks
- Fibroblast attention: spindle-shaped nuclei in stroma
- Lymphocyte attention: small round nuclei in immune infiltrates
- Cancer attention: large pleomorphic nuclei in epithelial regions
- Endothelial attention: elongated nuclei along vessel walls

### Quantitative Diagnostics
- MIL-predicted proportions vs Module 3 proportions (Pearson r per sample)
- Confidence weight distribution: how often does morphology dominate vs protein?
- Per-cell-type entropy: which types does MIL predict confidently?

## SLURM Strategy

- Array job: 12 tasks (one per sample)
- GPU partition (L40S or A100)
- Resources: 1 GPU, 8 CPUs, 64GB RAM, 4h time limit
- Stages 1-2 are GPU-intensive (Cellpose + ViT), Stages 3-6 are lighter

## 12 Patient Samples

| Sample | Type | Response |
|--------|------|----------|
| HCC22-088-P1-S1 | Biopsy | Progressor |
| HCC22-088-P1-S2 | Surgical | Progressor |
| HCC22-088-P2-S1 | Biopsy | Responder |
| HCC22-088-P2-S2 | Surgical | Responder |
| HCC22-088-P3-S1_A | Biopsy | Responder |
| HCC22-088-P3-S2 | Surgical | Responder |
| HCC22-088-P4-S1 | Biopsy | Progressor |
| HCC22-088-P4-S2_1i_rep | Surgical | Progressor |
| HCC22-088-P5-S1 | Biopsy | Responder |
| HCC22-088-P5-S2_F_rep | Surgical | Responder |
| HCC22-088-P6-S1 | Biopsy | Responder |
| HCC22-088-P6-S2_D | Surgical | Responder |

## Data Paths

- Patient data: `/ix1/alee/LO_LAB/General/Lab_Data/20250210_CITEGeistPublicData_GEO_Alex/processed_files/`
- H&E images: `{patient_dir}/outs/spatial/cytassist_image.tiff`
- Module 3 outputs: `output/module3_unified/HCC22-088-{SAMPLE}/`
- Module 3 proportions: `*_cell_prop_global_results.csv`
