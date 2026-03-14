# PC-MIL Xenium Benchmark: Lessons Learned

**Date:** 2026-03-14
**Context:** Implemented Protein-Conditioned MIL (PC-MIL) for single-cell type assignment on Xenium pseudo-Visium data with DAPI+boundary imaging.

## Summary

PC-MIL achieves **excellent proportion estimation** (r=0.756, exceeding the 0.743 Hybrid baseline) but **single-cell accuracy is fundamentally limited** (~21-24%) by the DAPI+boundary input signal. This is ~1.5-1.7x random (14.3%), confirming that nucleus morphology from fluorescence imaging alone cannot reliably distinguish cell types.

## Key Results

| Metric | Value | Baseline |
|--------|-------|----------|
| Proportion Pearson r | **0.756** | Hybrid: 0.743 |
| Single-cell accuracy | 21.5% (constrained) / 24.3% (global) | Random: 14.3% |
| Active cell types | 4-5 / 7 | Collapsed: 1 |

## What Worked

1. **Z-score normalization of protein signals** — Protein MSE was 127M× larger than proportion MSE. Z-scoring brought reconstruction loss to the same scale as proportion loss, preventing attention collapse.

2. **Reconstruction warmup (20 epochs)** — Letting proportion + diversity losses establish multi-type attention before reconstruction loss kicks in prevents early collapse.

3. **Hungarian assignment loss during training** — Using GT proportions → counts → optimal assignment → CE pseudo-labels improved proportion r from 0.65 to 0.74. The EM-style loop (attention → assignment → loss → better attention) works for proportions.

4. **Per-class gated attention** — K separate gate+score heads with softmax(dim=1) works correctly. The dim=0 vs dim=1 softmax bug is a known pitfall (caused r=0.017 in earlier MIL work).

## What Didn't Work

1. **Detection mask during training** — Applying detection mask (setting undetected type logits to -inf) during training caused partial attention collapse (active types dropped to 1.1-1.4). The -inf values propagate through diversity loss. Detection mask should ONLY be applied at inference.

2. **Morphology feature branch (3-dim sizes)** — Adding width/height/area features to the model had negligible impact. These features don't discriminate between cell types.

3. **Protein-context dropout (30%)** — Randomly zeroing protein input during training didn't force meaningful image feature learning. The model simply performed worse on those steps without developing compensating image representations.

4. **Temperature-controlled Hungarian pseudo-labels** — Annealing temperature from 2.0→0.5 had no measurable impact on either proportion r or SC accuracy.

5. **Ensemble (PC-MIL + morphology Gaussians)** — Combining PC-MIL attention log-probs with Xenium cell morphology Gaussian log-likelihoods at any weight (0-100%) gave ~21% SC accuracy. Both signals are at the constrained-random floor.

6. **Confidence-first greedy assignment** — Assigning high-confidence nuclei first, then using Hungarian for the rest, did not improve over pure Hungarian. The model's confidence reflects spot composition knowledge, not per-nucleus identity.

## Fundamental Insights

### Why proportion r is high but SC accuracy is low
The model learns "how many of each type in this spot" (proportion estimation) but NOT "which nucleus is which type" (single-cell identity). These are decoupled — you can get perfect proportions with random per-cell assignment.

### The constrained assignment paradox
- **Global classification** (no count constraints): 24.3% accuracy
- **Constrained Hungarian** (using count constraints): 21.2% accuracy
- Constraints HURT because they force some cells into wrong types even when the global classifier would have been right. The count constraints only help when the per-cell signal is strong enough to benefit from disambiguation.

### DAPI+boundary imaging ceiling
- Cellpose nucleus features (8-dim from masks): 17.9% with perfect counts
- Xenium cell features (12-dim with cytoplasm): 23.3% with perfect counts
- The boundary channel shows ALL cell membranes (not just the target cell's), making it non-discriminative at the per-nucleus level
- ViT features from DAPI+boundary: ~21% (essentially constrained random)

### H&E is fundamentally richer
The Visium HD H&E pilot achieved **56% SC accuracy** with the same MIL architecture because:
- Hematoxylin stains nuclei with cell-type-specific chromatin patterns
- Eosin stains cytoplasm with type-specific morphology
- RGB color provides 3 truly independent information channels
- Cell boundaries are visible in the tissue context

## Architecture Summary (for future reference)

```
image_features (N, 384) → projection → (N, 64)
protein_props (K,) → encoder → (32,) → broadcast to N nuclei
[optional] morph_features (N, 3) → projection → (N, 16)
concat → (N, 96-112) → per-class gated attention → (N, K) logits
softmax(dim=1) → attention → mean(dim=0) → proportions (K,)
proportions @ profile_matrix → reconstructed protein (M,)
```

Loss: proportion MSE + z-scored reconstruction MSE + entropy + diversity + Hungarian CE

## Files Created

| File | Purpose |
|------|---------|
| `CITEgeist/model/pc_mil.py` | PCMILModel architecture |
| `CITEgeist/model/pc_mil_training.py` | Training loop with 5-component loss |
| `CITEgeist/model/pc_mil_inference.py` | Inference with detection mask + ensemble |
| `CITEgeist/tests/test_pc_mil.py` | 11 unit tests |
| `Benchmarking/.../benchmark_pc_mil.py` | 5-fold CV benchmark harness |
| `Benchmarking/.../extract_vit_features.py` | ViT-S feature extraction |

## Next Steps

1. **Apply PC-MIL to Visium HD H&E data** — The architecture is proven; H&E provides the richer signal needed for SC accuracy. Expect significant improvement over the 56% MIL baseline due to Hungarian + protein conditioning.

2. **Keep single-cell assignments for both modalities** — Even at ~24% on Xenium DAPI, it's still 1.7× random and provides meaningful spatial information when aggregated.

3. **Consider Cellpose cyto model** — Running `model_type='cyto'` on H&E (hematoxylin=nucleus, eosin=cytoplasm) would give full cell segmentation with nc_ratio and cytoplasm features, potentially boosting SC accuracy further.

4. **Global classification + selective constraints** — Use global argmax for all cells, but apply count constraints only when spot-level proportions are confident and the global prediction is uncertain. This avoids the constrained assignment paradox.
