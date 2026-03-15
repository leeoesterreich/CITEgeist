# Unified PC-MIL Pipeline: Patient + Xenium Single-Cell Assignment

**Date**: 2026-03-14
**Status**: Design approved

## Overview

A unified pipeline that runs CITEgeist Module 3 (deconvolution) through PC-MIL (single-cell assignment) on both Xenium pseudo-Visium regions (DAPI + boundary) and the 12 patient Visium CITE-seq samples (H&E). The same codebase, architecture, and supervision strategy handles both imaging modalities.

## Design Decisions

| Decision | Choice | Rationale |
|----------|--------|-----------|
| ViT backbone | ImageNet ViT-S/16 (frozen) | UNI showed no improvement on Visium HD; accessible without HuggingFace auth |
| Supervision | Module 3 finetuned proportions | Self-contained — no external ground truth needed for training |
| Training scope | Per-sample | Avoids mixing biology across patients/timepoints |
| Cancer types | Single "Cancer" (EPCAM-1) | No need for Luminal/Basal split at this level |
| Unknown cell type | Disabled | Simplifies downstream assignment |
| Count constraints | None at inference | Global classification outperformed constrained assignment in benchmarks |
| Evaluation | Single-cell marker gene validation | No ground truth for patient data; checks biological consistency |

## Cell Type Profile (9 types, K=9)

```python
CELL_PROFILES = {
    "Cancer": {"Major": ["EPCAM-1"]},
    "Macrophages": {"Major": ["CD68-1", "CD163-1"]},
    "CD8_T_Cells": {"Major": ["CD3E-1", "CD8A-1"]},
    "CD4_T_Cells": {"Major": ["CD4-1"]},
    "B_Cells": {"Major": ["CD19-1"]},
    "Endothelial": {"Major": ["PECAM1-1"]},
    "Fibroblasts": {"Major": ["ACTA2-1"]},
    "Monocytes": {"Major": ["CD14-1"]},
    "Dendritic_Cells": {"Major": ["ITGAX-1", "HLA-DRA-1"]},
}
```

Note: The codebase requires nested dict format with `"Major"` keys (validated by `validate_cell_profile_dict()`).

## Pipeline Steps

### Step 1: Module 3 Re-run (Proportions + GEX Deconvolution)

**Input**: Raw SpaceRanger output (patient) or Xenium pseudo-Visium AnnData
**Output**: 9-type proportions + deconvolved GEX parquet

- Cell profile dict: 9 types (Cancer merged, no Luminal/Basal split)
- Unknown cell type: disabled by omitting "Unknown" from the profile dict (the `unknown_threshold` parameter only applies when "Unknown" is a profile key)
- Cellpose nuclei counts computed within this step via `compute_spot_nuclei_counts_cellpose()`
- Nuclei counts injected as prior to guide proportion estimation
- Cellpose model: `nuclei` for DAPI (Xenium), `cyto2` for H&E (patient)
- **Cellpose masks are cached** and reused by Step 2 (no duplicate segmentation run)

### Step 2: Cellpose Segmentation + Patch Extraction

**Input**: Tissue images (DAPI+boundary or H&E RGB)
**Output**: Nucleus masks, centroids, 224x224 patches

**Cellpose configuration:**
- Xenium: 2-channel input (DAPI + boundary), `model_type='nuclei'`, `channels=[1, 2]`
- Patient H&E: 3-channel RGB input, `model_type='cyto2'`, `channels=[0, 0]`

**Patch extraction:**
- 224x224 crop centered on each nucleus centroid
- Xenium: 2-channel patches zero-padded to 3 channels for ViT
- H&E: 3-channel patches fed directly to ViT
- Each nucleus assigned to nearest Visium spot (within spot radius); nuclei outside spots discarded

### Step 3: ViT Feature Extraction

**Input**: 224x224x3 patches
**Output**: 384-dim feature vectors per nucleus

- Frozen ImageNet ViT-S/16
- GPU-accelerated batch inference (~30x faster than CPU)
- Output: `vit_features.npy` (N_nuclei x 384) + `nucleus_ids.npy`

### Step 4: PC-MIL Training + Inference (Per-Sample)

**Input**: ViT features + Module 3 proportions + spot-level antibody signal
**Output**: Per-nucleus cell type assignments

#### Architecture
- Image projection: 384 → 64 (Linear + ReLU)
- Protein encoder: K → 32 (Linear + ReLU)
- Fusion: concat → 96-dim per nucleus
- Per-class gated attention: K=9 separate heads (gate + score networks on 96-dim)
- Output: softmax over types (dim=1) → per-nucleus soft assignments
- Aggregation: mean across nuclei → predicted spot proportions

#### Loss (multi-task)
- Proportion MSE: predicted vs Module 3 finetuned proportions
- Protein reconstruction MSE: proportions @ learnable profile matrix vs observed antibody signal
- Entropy regularization: prevent uniform assignments
- Diversity loss: prevent ignoring present types

#### Training
- All spots from one sample as training data
- No train/val split (evaluated post-hoc via marker genes)
- Early stopping on training loss plateau, with `max_epochs=200` hard cap (Visium HD showed overfitting beyond ~100 epochs)
- Attention hidden dim: 128

#### Inference
- Per-nucleus: argmax of soft assignment → discrete cell type label
- No count constraints enforced

### Step 5: Marker Gene Validation (Single-Cell Level)

**Input**: Per-nucleus assignments + Module 3 deconvolved GEX parquet
**Output**: Per-nucleus marker scores + summary metrics

For each nucleus assigned type T in spot S:
1. Pull deconvolved GEX row `S:::T` from Module 3 parquet
2. Extract expression for known marker genes of type T
3. Extract expression for marker genes of other types (negative controls)
4. Compute marker score: assigned-type markers vs non-assigned-type markers

#### RNA Marker Validation Dictionary

| Cell Type | Expected Markers |
|-----------|-----------------|
| Cancer | EPCAM |
| Macrophages | CD68, CD163, CSF1R |
| CD8_T_Cells | CD8A, CD8B, GZMB |
| CD4_T_Cells | CD4, IL7R, CCR7 |
| B_Cells | CD19, MS4A1, CD79A |
| Fibroblasts | COL1A1, DCN, VIM |
| Endothelial | PECAM1, VWF, CDH5 |
| Monocytes | CD14, FCGR3A, S100A8 |
| Dendritic_Cells | ITGAX, HLA-DRA, CLEC10A |

#### Summary Metrics
- Per-type: median marker score, fraction of nuclei where assigned-type markers > other-type markers
- Overall: weighted average across types
- Visualization: violin plots of marker scores by assigned type

## Output Structure

```
output/unified_pipeline/{sample_name}/
├── module3/
│   ├── cell_prop_finetuned_results.csv
│   ├── gene_expression_pass1.parquet
│   ├── module3_meta.json
│   └── module3_results.h5ad
├── cellpose/
│   ├── nuclei_masks.npy
│   ├── nuclei_centroids.csv
│   └── nuclei_per_spot.csv
├── features/
│   ├── vit_features.npy
│   ├── nucleus_ids.npy
│   └── patches.npy
├── pcmil/
│   ├── model_weights.pt
│   ├── assignments.csv
│   └── training_log.json
└── validation/
    ├── marker_gene_scores.csv
    └── validation_summary.json
```

## SLURM Execution

Four sequential array jobs with marker files between steps:

| Step | Job Type | Resources | Time | Dependencies |
|------|----------|-----------|------|-------------|
| 1. Module 3 | Array (12 patient + 5/14 Xenium) | 64GB RAM, CPU, Gurobi | ~6h | None |
| 2. Patches + ViT | Array (same) | 32GB RAM, 1 GPU | ~30min | Step 1 complete |
| 3. PC-MIL train + infer | Array (same) | 32GB RAM, 1 GPU | ~1-2h | Step 2 complete |
| 4. Validation | Array (same) | 16GB RAM, CPU | ~15min | Step 3 complete |

Marker files (e.g., `output/unified_pipeline/{sample}/.step1_complete`) used to gate subsequent steps since cross-cluster SLURM dependencies don't work on this HPC.

## Xenium-Specific Evaluation

On Xenium pseudo-Visium regions, ground truth from protein gating is available for additional evaluation beyond marker genes:
- Per-nucleus accuracy vs Cellpose-Xenium spatial matches
- Proportion Pearson r vs ground truth proportions
- Confusion matrix by cell type

This is evaluation only — ground truth is never used during training.

## Cohort

### Patient Samples (12)
- **Responders (8):** P2-S1, P2-S2, P3-S1_A, P3-S2, P5-S1, P5-S2_F_rep, P6-S1, P6-S2_D
- **Progressors (4):** P1-S1, P1-S2, P4-S1, P4-S2_1i_rep

### Xenium Pseudo-Visium Regions
- 14 regions (all available); use 5 for quick iteration during development
- DAPI + boundary morphology
- Ground truth available for held-out evaluation

## Key Code Changes Required

1. **Unified PC-MIL class**: Merge Xenium and Visium HD PC-MIL into single `UnifiedPCMIL`
2. **Module 3 runner**: New script with 9-type profile, unknown disabled, nuclei count injection
3. **Cellpose wrapper**: Modality-aware (DAPI vs H&E) with shared output format
4. **Patch extraction**: Handle 2-channel (Xenium) and 3-channel (H&E) inputs
5. **Validation module**: Marker gene scoring at single-cell level
6. **SLURM scripts**: 4 array job templates with marker file gating

## Code Gaps (Current State → Target State)

These are gaps between the existing codebase and what this spec requires. The implementation plan must address each.

### Critical

1. **Profile dict schema mismatch in PC-MIL**: Module 3 uses nested `{"Major": [...]}` dicts, but `build_profile_matrix()` in `pc_mil.py` expects flat `{type: [markers]}`. Passing nested dicts will iterate over `"Major"` instead of marker names, producing broken profile rows. **Fix**: Add an adapter that flattens nested profile dicts for PC-MIL's profile matrix builder.

2. **Inference mode hardcoded to Hungarian**: Spec says argmax with no count constraints, but `pc_mil_infer_spot()` in `pc_mil_inference.py` always applies largest-remainder discretization + Hungarian constrained assignment. **Fix**: Add an `inference_mode` parameter (`"argmax_global"` vs `"hungarian_constrained"`) and default to `"argmax_global"` for this pipeline.

3. **Cellpose API not modality-aware**: Spec requires `model_type='nuclei'` + `channels=[1,2]` for Xenium and `model_type='cyto2'` + `channels=[0,0]` for H&E, but `compute_spot_nuclei_counts_cellpose()` in `citegeist_model.py` does not expose `model_type` or `channels` — it's hardcoded to `channels=[0,0]` on RGB. **Fix**: Add `model_type` and `channels` parameters to the Cellpose API, or create a new modality-aware wrapper.

### Warnings

4. **Trainer requires val_dataset**: Current `pc_mil_training.py` requires a `val_dataset` and early-stops on `val_r`. Spec says no train/val split. **Fix**: Make `val_dataset` optional; when absent, use training loss plateau + `max_epochs=200` hard cap for stopping.

5. **Loss hyperparameter drift**: Code defaults to `lambda_diversity=0.5` (strong regularization). Spec doesn't pin values. **Fix**: Document and pin loss weights in the pipeline config to ensure reproducibility.

6. **Feature extractor output format**: Spec expects unified `vit_features.npy` (all nuclei in one array), but current `extract_vit_features.py` writes per-spot files. **Fix**: Write a unified array with matching `nucleus_ids.npy` index, or adapt PC-MIL data loader to read per-spot files.

7. **Core API not wired to PC-MIL**: `run_single_cell_resolution()` in `citegeist_model.py` dispatches to older constrained assignment/MIL paths, not PC-MIL. **Fix**: Add PC-MIL as a dispatch option, or bypass the core API with a standalone orchestrator script.

## Error Handling

- **Zero nuclei in spot**: Skip spot during PC-MIL training; proportion for that spot comes from Module 3 only
- **Too few spots for training**: Warn if sample has <100 spots; proceed but flag in validation summary
- **Module 3 validation failure**: Log warning and continue (set `validation_warn_only=True`)
- **Corrupt patches / ViT errors**: Skip nucleus, log warning, exclude from assignment
