# Design: Spatial GMM Component Characterization for Thresholding

## Date: 2026-02-02

## Problem

PanCK has a 3-component GMM (background mean=25, dim mean=115, bright mean=368) with threshold at 59.4. Epithelial cells (median expression ~20) are absorbed into the background component. The sub-threshold Moran's I approach (v1) detected spatial signal below threshold but refitting a 2-comp GMM on sub-threshold values still placed the boundary too high (30.3), yielding Epithelial recall stuck at 0.66.

The fundamental issue: BIC picks 3 components, but the lowest component is a mixture of true noise AND dim epithelial cells. We need to split it further until we isolate a pure noise component.

## Approach: Per-Component Spatial Characterization

Instead of treating GMM components as "background/dim/bright" by expression level, classify each as **noise vs signal** by its spatial character:

- **Noise**: spatially random (flat field / fuzzy TV screen). Moran's I ~ 0.
- **Signal**: spatially structured (puncta, regions of dominance). Moran's I >> 0.

Then threshold at the bottom of the lowest signal component's bell curve.

## Algorithm

For each marker (when `spatial_coords` is available):

1. **Fit GMMs** with k=1..5 components. Store all fitted models.
2. **Start with BIC-selected k**. Sort components by mean ascending.
3. **Walk components from lowest mean upward**, computing Moran's I on posterior-weighted expression for each:
   - `z_k(i) = P(component k | x_i) * x_i` for all cells
   - Compute Moran's I on `z_k` using spatial KNN (k=15, binary weights)
   - 199 permutations for p-value
   - **No spatial structure** (I < I_min or p > 0.05) -> label as "noise", continue up
   - **Has spatial structure** (I > I_min and p < 0.05) -> this is the **lowest signal component**. Stop walking.
4. **Set threshold**: expression value where P(lowest signal component | x) = 0.1. This captures the component's full bell curve down to its meaningful tail.
5. **If the lowest component already has spatial structure** (no clean noise floor found at this k):
   - Increment k by 1, refit GMM, repeat from step 3
   - This splits the lowest component, hopefully separating noise from dim signal
6. **Stopping conditions**:
   - Found a noise/signal boundary -> done, use P=0.1 threshold
   - Reached k=5 and lowest component still has spatial structure -> fall back to standard BIC posterior crossing threshold
   - All components are noise (no spatial structure at any level) -> marker is unimodal noise. Use BIC fallback threshold (90th percentile of fitted Gaussian)
   - k=1 selected by BIC -> skip spatial characterization, use existing unimodal logic

## Threshold Computation Detail

For the lowest signal component (say component j with mean mu_j):
- Evaluate GMM posteriors on a grid from 0 to mu_j (1000 points)
- Find x where P(component j | x) = 0.1
- This is the threshold: below here, even if a cell belongs to component j, we can't distinguish it from noise

## Expected Behavior

### PanCK (current: 3-comp, means 25/115/368, threshold 59.4)
- k=3: Component 1 (mean=25) -> check Moran's I. Likely has spatial structure (contains dim epithelial). Increment to k=4.
- k=4: Might split into (15, 40, 115, 368). Component 1 (mean=15) -> noise (no spatial structure). Component 2 (mean=40) -> signal (spatially structured dim epithelial).
- Threshold = P(comp2 | x) = 0.1, landing around ~20.
- Result: captures dim epithelial cells, Epithelial recall should improve.

### CD3E (current: 3-comp, means 21/313/967, threshold 60.5)
- k=3: Component 1 (mean=21) -> likely noise (T cells are punctate, background is random). Component 2 (mean=313) -> signal.
- Threshold = P(comp2 | x) = 0.1. Given the large gap between 21 and 313, this might land around 60-80.
- Result: similar to current threshold. Unchanged or slightly adjusted.

### CD68 (current: 3-comp, means 28/171/613, threshold 73.5)
- Similar to CD3E: large gap between background and signal. Threshold likely similar.

## Implementation

### Only file changed: `CITEgeist/model/cell_classification.py`

### New function: `_classify_components_by_spatial_structure()`

```python
def _classify_components_by_spatial_structure(
    expr: NDArray[np.floating],
    gmm: GaussianMixture,
    spatial_coords: NDArray[np.floating],
    spatial_k: int = 15,
    n_perm: int = 199,
    morans_i_min: float = 0.05,
    pvalue_max: float = 0.05,
    seed: int = 42,
) -> Tuple[List[bool], List[float], List[float]]:
    """
    Classify each GMM component as noise (False) or signal (True)
    by computing Moran's I on posterior-weighted expression.

    Returns:
        (is_signal, morans_i_values, pvalues) - per component
    """
```

### New function: `_find_spatial_threshold()`

```python
def _find_spatial_threshold(
    expr: NDArray[np.floating],
    spatial_coords: NDArray[np.floating],
    max_components: int = 5,
    spatial_k: int = 15,
    n_perm: int = 199,
    morans_i_min: float = 0.05,
    pvalue_max: float = 0.05,
    posterior_cutoff: float = 0.1,
    seed: int = 42,
) -> Optional[MarkerThreshold]:
    """
    Iterative spatial GMM thresholding.
    Increments k until a clean noise/signal boundary is found.
    """
```

### Modified: `_compute_subthreshold_morans_i()`
- **Replaced** by per-component Moran's I inside `_classify_components_by_spatial_structure()`
- The old function is removed

### Modified: `_fit_gmm_bic()`
- Extend max_components parameter range to 5
- Return all fitted models (not just BIC-best) for spatial characterization

### Modified: `determine_thresholds()`
- When spatial_coords available: use `_find_spatial_threshold()` as primary method
- Falls back to standard BIC threshold on failure

### Modified: `_determine_thresholds_per_cluster()`
- Remove sub-threshold Moran's I second-chance block (replaced by new primary method)
- Keep per-cluster adaptive logic (CV + gain) as independent mechanism

### Modified: `AdaptiveThresholdInfo`
- Replace `subthreshold_morans_i` / `subthreshold_pvalue` / `original_threshold` with:
  - `component_morans_i: Optional[Dict[int, float]]` - per-component Moran's I
  - `component_is_signal: Optional[Dict[int, bool]]` - noise/signal classification
  - `noise_signal_boundary_k: Optional[int]` - k at which boundary was found
  - `spatial_threshold_method: Optional[str]` - "spatial_gmm", "bic_fallback", etc.

### Modified: `MarkerThreshold`
- New field: `spatial_component_info: Optional[Dict]` - component-level spatial characterization

### Serialization
- All new fields Optional with None defaults -> backward compatible
- `ThresholdSet.save/load` extended for new fields

### Report
- `generate_threshold_report()` adds per-component Moran's I to summary CSV
- Diagnostic plots show component-by-component spatial classification

## Verification

Submit all 5 regions (`sbatch --array=0-4`):
- PanCK should get a lower threshold (~20) via spatial GMM
- Epithelial recall should improve from 0.66
- CD3E, CD68 should remain similar (large gap between noise and signal components)
- Overall accuracy should stay stable or improve
- No markers should regress significantly vs current BIC-only approach
