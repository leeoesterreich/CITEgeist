# Protein-Conditioned MIL (PC-MIL) Design Spec

## Problem Statement

CITEgeist achieves strong spot-level cell type proportion estimation (r=0.74 on Xenium pseudo-Visium) using protein deconvolution, but extending this to **single-cell resolution** remains unsolved. The current MIL approaches suffer from **attention collapse** — the model predicts only 2-3 dominant cell types and ignores rare types, particularly on imbalanced patient data.

**Root cause**: Existing MIL heads see image features in isolation and must independently rediscover cell type information that the protein channel already provides. Evidence: constrained Hungarian assignment (which uses protein-derived counts as hard constraints) achieves 46% single-cell accuracy vs 25% for unconstrained morphology classification — proving that protein-as-constraint dramatically reduces problem difficulty.

**Goal**: A learned model that fuses spot-level protein proportions with per-nucleus image features to produce accurate single-cell type assignments and downstream gene expression deconvolution.

## Design Decisions

| Decision | Choice | Rationale |
|----------|--------|-----------|
| Evaluation target | Xenium pseudo-Visium (validate), patient data (deploy) | If Xenium works but patient fails, pseudo-Visium construction is wrong |
| Modalities | Protein (spot) + morphology (nucleus) → type + GEX | Core problem definition |
| Protein-image coupling | Soft constraint (protein as input feature) | Most learnable; model discovers optimal weighting |
| Training objective | Multi-task (proportion + protein reconstruction) | Denser gradients; self-consistency loop |
| Image backbone | Frozen ViT-S (ImageNet pretrained) for both modalities | Proven (r=0.65); DAPI=ch0+boundary+boundary, H&E=RGB natively |
| Design complexity | Full design upfront | Existing components available; avoid dead-end iterations |
| Fractional assignment | Detection mask at inference; soft assignments during training | Prevents noise from sub-threshold types without losing gradient signal |

## Architecture

### Protein-Conditioned MIL (PC-MIL)

```
+----------------------------------------------------------+
|                     Per Spot                               |
|                                                           |
|  Protein proportions (K,)                                 |
|       |                                                   |
|       v                                                   |
|  ProteinEncoder: Linear(K, 32) -> ReLU -> Linear(32, 32)  |
|       |                                                   |
|       v                                                   |
|  protein_context (32,)  ---- broadcast to all N nuclei    |
|                                                           |
|  +-----------------------------------------------------+  |
|  |              Per Nucleus (xN)                        |  |
|  |                                                      |  |
|  |  Image patch (3, 224, 224)                           |  |
|  |       |                                              |  |
|  |  Frozen ViT-S -> 384-dim                             |  |
|  |       |                                              |  |
|  |  ImageProjection: Linear(384, 64) -> LayerNorm       |  |
|  |       |                                              |  |
|  |  image_emb (64,)                                     |  |
|  |       |                                              |  |
|  |  concat(image_emb, protein_context) -> (96,)         |  |
|  |       |                                              |  |
|  |  +------------------------------------------------+  |  |
|  |  |  Per-Class Gated Attention (xK heads)          |  |  |
|  |  |                                                 |  |  |
|  |  |  gate_k = sigmoid(W_gate_k @ fused + b_gate)   |  |  |
|  |  |  score_k = tanh(W_score_k @ fused + b)         |  |  |
|  |  |  logit_k = gate_k * score_k                    |  |  |
|  |  +------------------------------------------------+  |  |
|  |       |                                              |  |
|  |  logits (K,) per nucleus                             |  |
|  +-----------------------------------------------------+  |
|                                                           |
|  Stack logits -> (N, K)                                   |
|       |                                                   |
|  softmax(dim=1) -> attention (N, K)                       |
|       |                                                   |
|  +-- mean(dim=0) -> predicted_proportions (K,)            |
|  |                                                        |
|  +-- predicted_proportions @ protein_profile_matrix ->     |
|       reconstructed protein (M,) per spot                 |
+----------------------------------------------------------+
```

### Key Dimensions

- N = nuclei per spot (variable, ~5-30)
- K = cell types (10 for Xenium benchmark)
- M = protein markers (~15-20)
- Image embedding: 64-dim
- Protein context: 32-dim
- Fused: 96-dim

### Protein Reconstruction Branch

Per-nucleus attention weights (N, K) are averaged to proportions (K,), then multiplied by a learnable **protein profile matrix** (K, M) to reconstruct the spot's protein signal. This creates a self-consistency loop: wrong assignments -> bad reconstruction -> gradient correction.

The protein profile matrix is initialized from Module 2's `cell_profile_dict`. Since `cell_profile_dict` is membership-based (`{cell_type: [marker_list]}`), initialization constructs a binary (K, M) matrix where `profile[k, m] = 1.0` if marker `m` belongs to type `k`, else `0.0`. Marker and cell type ordering must be explicitly aligned with the training data columns. The matrix is then allowed to fine-tune during training to learn continuous profile values.

### Input Channels

- **Xenium (DAPI + boundary)**: 3 channels = [DAPI, boundary, boundary]
- **Patient H&E**: 3 channels = [R, G, B] natively

Single frozen ViT-S handles both; only the front-end stacking differs.

## Loss Functions

### Primary Losses

**1. Proportion loss (L_prop)**
```python
L_prop = MSE(predicted_proportions, true_proportions)
```
MSE over KL divergence — KL blows up when predicted proportions approach 0 for present types, destabilizing early training.

**2. Protein reconstruction loss (L_recon)**
```python
# protein_profiles: (K, M) learnable matrix, initialized from cell_profile_dict
reconstructed = predicted_proportions @ protein_profiles  # (M,)
L_recon = MSE(reconstructed, observed_protein_signal)
```
Self-consistency loop: if assignments are wrong, protein reconstruction fails, providing dense gradient signal. The `observed_protein_signal` is CLR-normalized, consistent with Module 3 preprocessing (`preprocess_antibody()`).

### Regularization Losses

**3. Assignment entropy loss (L_entropy)**
```python
per_nucleus_entropy = -sum(attention * log(attention + eps), dim=1)  # (N,)
L_entropy = per_nucleus_entropy.mean()
```
Prevents uniform-assignment collapse. Each nucleus should commit to 1-2 types.

**4. Diversity loss (L_diversity)**
```python
mean_attention = attention.mean(dim=0)  # (K,) = predicted proportions
active_mask = (true_proportions > 0.01).float()
L_diversity = -(active_mask * log(mean_attention + 0.01)).sum()
```
Prevents dominant-type collapse. If a type is present in the spot, at least some nuclei must be assigned to it. The `+ 0.01` floor prevents log(0) while preserving gradient flow — unlike `clamp(min=0.01)` which would zero the gradient exactly when recovery is needed most. Combined with `lambda_diversity = 0.05`, this keeps the regularization gentle.

### Combined Loss
```python
L = L_prop + lambda_recon * L_recon + lambda_entropy * L_entropy + lambda_diversity * L_diversity
```

**Default weights:**
- lambda_recon = 1.0 (equal to proportion loss)
- lambda_entropy = 0.1 (gentle push toward confidence)
- lambda_diversity = 0.05 (light touch, prevent total collapse)

## Training Pipeline

### Data Flow

1. **Pre-extract features (offline, once):**
   - Run frozen ViT-S on all nucleus patches -> save (N_total, 384) tensor
   - Store spot->nucleus mapping (which nuclei belong to which spot)
   - Store spot-level protein proportions (from Module 3 continuous/hybrid) and raw protein signals

2. **Training (online):**
   - DataLoader yields per-spot batches:
     `{image_features: (N_i, 384), protein_props: (K,), protein_signal: (M,), true_props: (K,)}`
   - Variable N_i per spot handled via collate with padding + mask

### Training Parameters

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| Optimizer | AdamW | Standard for transformer-adjacent architectures |
| Learning rate | 1e-3 | Only training small head (~50K params) |
| Weight decay | 1e-4 | Light regularization |
| Batch size | 64 spots | ~300-1000 nuclei per batch |
| Epochs | 200 | With early stopping (patience=20 on val proportion r) |
| LR schedule | Cosine annealing to 1e-5 | Smooth decay |
| Gradient clipping | max_norm=1.0 | Prevent spikes from reconstruction loss |

### Class Imbalance Handling

1. **Inverse-frequency spot weighting**: Spots where rare types are present get upweighted. Weight = 1/frequency of rarest type in that spot.
2. **Diversity loss**: Directly prevents ignoring present types.
3. **Stratified batching**: Each batch sampled to ensure representation of all K types.

### Initialization

- **Protein profile matrix**: From Module 2's `cell_profile_dict`
- **Per-class gate/score weights**: Xavier initialization
- **Image projection**: Xavier initialization

## Fractional Assignment Handling

**Problem**: A spot might predict 0.05 proportion of macrophages (= 0.5 cells at 10 nuclei). Assigning even one nucleus to a sub-threshold type adds noise.

**Solution (two levels):**

- **Training time**: Soft assignments are fine. 0.05 probability contributes proportionally to loss. No issue.
- **Inference time**: Use the existing `detection.py` GMM module as a hard gate:

```python
detected = detection_gmm(protein_signal)  # binary (K,) per spot
masked_logits = logits.clone()
masked_logits[:, ~detected] = -inf  # prevent assignment to undetected types
attention = softmax(masked_logits, dim=1)
```

During training, the detection mask is **NOT applied** — all types receive gradients via soft assignments. The mask is only used at inference time. This avoids train/inference distribution mismatch while keeping gradient flow to all type heads. The protein reconstruction loss naturally provides signal about which types are absent (reconstruction fails if absent types get assigned nuclei).

## Inference Pipeline

```
Per spot:
  1. Detection gate: detected = detect_cell_types(X, marker_groups, ...)  # returns (n_spots, n_types) bool
  2. Forward pass: image_features + protein_context -> fused -> per-class gated attention -> logits
  3. Mask undetected (with all-false fallback):
     if detected[spot].any():
         logits[:, ~detected[spot]] = -inf
     else:
         pass  # no masking if detection finds nothing — prevents NaN softmax
  4. Softmax: attention = softmax(logits, dim=1)
  5. Proportions: proportions = attention.mean(dim=0)
  6. Single-cell assignment:
     counts = largest_remainder(proportions, n_nuclei)
     assignments = hungarian(attention, counts)  # constrained_assignment.py negates internally
  7. Aggregate to spot counts: spot_counts = assignments.value_counts() per type
  8. GEX deconvolution: run_cell_expression_pass1(use_discrete_mode=True, cell_counts=spot_counts)
```

### Output Format

```
nucleus_id | barcode | cell_type | confidence | image_score | protein_score
nuc_001    | AACAC.. | CD4_T     | 0.82       | 0.71        | 0.89
nuc_002    | AACAC.. | Fibroblast| 0.94       | 0.91        | 0.95
nuc_003    | AACAC.. | Macrophage| 0.45       | 0.38        | 0.67
```

- **confidence**: attention weight for assigned type (from softmax output)
- **image_score**: run a second forward pass with zero protein context; the attention weight for the assigned type measures image-only contribution
- **protein_score**: `protein_proportions[assigned_type]` — the spot-level protein proportion for that type, measuring how much protein alone supports this assignment

## Pipeline Integration

PC-MIL slots after Module 3 Pass 1 (continuous/hybrid proportions), before Pass 2 (GEX deconvolution):

```
Module 1 -> Module 2 -> Module 3 Pass 1 (continuous proportions, r=0.74)
                            |
                            +-- Current: Hungarian with morphology (46% sc accuracy)
                            |
                            +-- New: PC-MIL (proportions + image -> refined assignments)
                                      |
                                      v
                                Module 3 Pass 2 (GEX deconvolution)
                                      |
                                      v
                                Module 4 (program discovery)
```

PC-MIL **consumes** the continuous/hybrid proportion output as its protein conditioning signal. No changes to Modules 1, 2, 4, or 5.

### Integration Method

One new method on CitegeistModel:

```python
def run_pc_mil_assignment(self,
                          patch_features_path: str,
                          nucleus_spot_mapping: pd.DataFrame,
                          n_epochs: int = 200,
                          device: str = 'cuda') -> pd.DataFrame:
    """
    Protein-Conditioned MIL assignment.
    Uses Module 3 continuous proportions + nucleus image features
    to produce per-nucleus cell type assignments.

    Returns DataFrame: nucleus_id, barcode, cell_type, confidence
    """
```

## Evaluation & Success Criteria

### Metrics

**Spot-level (proportion quality):**
- Pearson r vs ground truth
- Baseline: continuous/hybrid = 0.74
- Target: >= 0.72 (at minimum equal — we're refining assignments, not degrading proportions)

**Single-cell (assignment quality):**
- Accuracy vs Xenium ground truth (protein-gated cell types)
- Baseline: constrained Hungarian = 46%
- Target: >= 50%
- Per-type accuracy breakdown (if any type drops below random, collapse is occurring)

**Collapse detection (monitored during training):**
- Active types per epoch (types where any nucleus has max attention > 0.3)
- If active types < K, collapse is happening
- Entropy of mean attention across all spots

### Ablation Experiments (on Xenium)

| Experiment | What it tests |
|-----------|---------------|
| PC-MIL (full) | Complete model |
| No protein conditioning | Image-only MIL — isolates protein's contribution |
| No reconstruction loss | Proportion loss only — tests if reconstruction prevents collapse |
| No diversity loss | Tests if dominant-type collapse returns |
| No detection mask | Tests if fractional assignment noise hurts |
| Frozen vs learnable protein profiles | Tests if profile drift helps or hurts |

### Edge Cases

- **Spots with zero nuclei**: Skip entirely — no assignment possible. Output empty rows.
- **Detection false negatives for rare types**: The all-false fallback (no masking) handles this. Additionally, the diversity loss during training penalizes ignoring present types.
- **Padded nuclei in variable-length batches**: Mask padded positions in the loss computation — padded nuclei do not contribute to proportion mean or entropy/diversity losses.
- **Cell type/marker ordering**: Enforce strict ordering alignment at data loading time via explicit column name matching between `protein_props`, `true_props`, detection mask, and profile matrix.

### Failure Modes to Watch

1. **Protein shortcut**: Model ignores image, echoes protein proportions -> check "no protein conditioning" ablation
2. **Attention collapse**: All nuclei get ~uniform assignments -> monitor active types + entropy
3. **Xenium->patient gap**: Works on Xenium, fails on patient -> implies pseudo-Visium construction issue

### Validation Strategy

- **Xenium**: 5-fold cross-validation over 14 pseudo-Visium regions. Train on 11, validate on 3.
- **Patient**: Train/val split by patient. Report proportion r only (no single-cell ground truth).

## File Organization

### New Files

```
CITEgeist/model/
  pc_mil.py                    # PCMILModel class (architecture + forward pass)
  pc_mil_training.py           # Training loop, loss computation, validation
  pc_mil_inference.py          # Inference pipeline, detection masking, Hungarian output

Benchmarking/xenium_benchmarking/CITEgeist/
  src/benchmark_pc_mil.py      # Xenium benchmark harness (5-fold CV)
  slurm/sbatch_pc_mil.sh       # SLURM submission script
```

### Reused Existing Code

| Component | Existing file | Usage |
|-----------|--------------|-------|
| ViT-S feature extraction | `visiumhd_benchmarking/src/run_benchmark.py` | Extract and cache patch features |
| Detection GMM | `model/detection.py` | Inference-time type masking |
| Hungarian assignment | `model/constrained_assignment.py` | Convert attention -> discrete assignments |
| Protein profile dict | `model/citegeist_model.py` | Initialize protein profile matrix |
| Module 3 Pass 2 | `model/gurobi_impl.py` | GEX deconvolution from assignments |
| Patch extraction | `xenium_benchmarking/CITEgeist/src/prepare_patches.py` | Nucleus patch preparation |

### Relationship to Existing MIL Files

PC-MIL supersedes `single_cell_mil.py`, `direct_softmax_model.py`, and `proportion_mil.py` for production use. Those files are retained for ablation/historical reference but are not part of the active pipeline going forward.

### Changes to Existing Files

- **`citegeist_model.py`**: Add `run_pc_mil_assignment()` method as new entry point (additive, no existing methods modified)
- **`test_asymmetric_loss.py`**: Fix 4 unpack sites for API break (bug fix, see below)

## Bug Fix: Critical #1 (alongside this work)

**Issue**: `optimize_cell_proportions_per_marker()` in `gurobi_impl.py:941` now returns 5 values (`..., recon_error`), but `test_asymmetric_loss.py` still unpacks 4 values at lines 81, 105, 135, 145.

**Fix**: Update all 4 unpack sites in `test_asymmetric_loss.py` to accept the 5th return value:
```python
# Before:
Y_values, beta_new, marker_beta_dict, alpha_values = optimize_cell_proportions_per_marker(...)
# After:
Y_values, beta_new, marker_beta_dict, alpha_values, recon_error = optimize_cell_proportions_per_marker(...)
```

## Appendix: Why Not Other Approaches

### Cross-Attention Fusion (Approach 2)
Elegant but overkill — ~15 protein markers is too small a set for meaningful cross-attention. The attention over such a small set would be noisy. Concatenation achieves the same information flow with fewer parameters.

### Prototype-Anchored with Protein Prior (Approach 3)
Mathematically clean (Bayesian) but rigid. The log-prior + log-likelihood fusion assumes modality independence, which doesn't hold — fibroblast morphology correlates with its protein signature. Learned fusion captures these non-linear interactions.

### Hard Constraint (Current Hungarian)
No learning, no gradient flow. Simple but rigid — if protein counts are wrong, image can't correct them. PC-MIL's soft conditioning allows the model to learn when to trust each modality.
