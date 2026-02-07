# Module 2 Autodiscovery Fixes Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Fix Module 1/2 autodiscovery evaluation bugs and add singleton rescue to improve profile quality on Xenium benchmark.

**Architecture:** Three independent fixes: (1) increase permutation count for better FDR resolution, (2) create shared benchmark constants with achievable-7 GT replacing old 10-type GT, (3) add singleton rescue based on unique spatial coverage using Module 1 GMM signal masks. Multi-scale colocalization is already enabled by default — no change needed there.

**Tech Stack:** Python, numpy, scipy, sklearn (GaussianMixture), existing CITEgeist model infrastructure.

---

## Important Context

- `run_benchmark.py` uses `discover_hierarchical_profiles()` (line 3623 of spatial_colocalization.py)
- `evaluate_pipeline_stages.py` uses `discover_profiles()` (line 2831 of spatial_colocalization.py)
- These are DIFFERENT functions. Both need the permutation fix. Only the evaluation script needs the GT fix.
- Multi-scale colocalization (`multi_scale_k=[6, 12, 24, 48, 64]`) is already the default and already running. It was NOT disabled. Stromal markers fail to group in some regions due to genuine biological variability, not a missing parameter.
- The `_fit_gmm_per_marker()` function (marker_interest.py:115) fits per-marker GMMs but only returns aggregate `snr_values` and `signal_fractions` — not per-spot signal assignments needed for Fix 3.

---

### Task 1: Shared Benchmark Constants Module

**Files:**
- Create: `Benchmarking/xenium_benchmarking/benchmark_constants.py`

**Step 1: Create the shared constants module**

```python
"""
Shared constants for Xenium benchmarking.

Single source of truth for achievable-7 cell type definitions used by
both run_benchmark.py and evaluate_pipeline_stages.py.
"""

from typing import Dict, List, Set

# =============================================================================
# ACHIEVABLE-7 CELL TYPE DEFINITIONS
# =============================================================================
#
# 7 cell types achievable with the 27-antibody Xenium panel.
# Collapsed from 10 granular RNA-based ground truth types.
#
# Collapse rationale:
# - Myofibroblasts + Stromal → Fibroblasts (both VIM+, αSMA overlaps)
# - Vascular Stromal → Endothelial (CD31+)
# - Proliferating T → CD8+ T cells (both CD3E+)
# - Mixed Immune → CD4+ T cells (HLA-DR+ T cells)

ACHIEVABLE_7_CELL_PROFILE_DICT: Dict[str, Dict[str, List[str]]] = {
    "B cells": {
        "Major": ["CD20"],
        "Minor": ["CD45RA"],
    },
    "CD4+ T cells": {
        "Major": ["CD3E", "CD4"],
        "Minor": ["CD45RO"],
    },
    "CD8+ T cells": {
        "Major": ["CD3E", "CD8A"],
        "Minor": ["GranzymeB"],
    },
    "Macrophages": {
        "Major": ["CD68", "CD163"],
        "Minor": ["CD16"],
    },
    "Endothelial": {
        "Major": ["CD31"],
        "Minor": [],
    },
    "Epithelial": {
        "Major": ["PanCK"],
        "Minor": [],
    },
    "Fibroblasts": {
        "Major": ["alphaSMA", "Vimentin"],
        "Minor": [],
    },
}

# Set-based signatures for Jaccard matching of autodiscovered profiles
ACHIEVABLE_7_MARKER_SIGNATURES: Dict[str, Set[str]] = {
    "B cells": {"CD20", "CD45RA"},
    "CD4+ T cells": {"CD3E", "CD4", "CD45RO"},
    "CD8+ T cells": {"CD3E", "CD8A", "GranzymeB"},
    "Macrophages": {"CD68", "CD163", "CD16"},
    "Endothelial": {"CD31"},
    "Epithelial": {"PanCK"},
    "Fibroblasts": {"alphaSMA", "Vimentin"},
}

# Primary/secondary format for profile scoring in staged evaluation
# Scoring formula: (2 * primary_overlap + secondary_overlap) / 3
ACHIEVABLE_7_GT_MARKERS: Dict[str, Dict[str, List[str]]] = {
    "B cells": {
        "primary": ["CD20"],
        "secondary": ["CD45RA"],
    },
    "CD4+ T cells": {
        "primary": ["CD3E", "CD4"],
        "secondary": ["CD45RO"],
    },
    "CD8+ T cells": {
        "primary": ["CD3E", "CD8A"],
        "secondary": ["GranzymeB"],
    },
    "Macrophages": {
        "primary": ["CD68", "CD163"],
        "secondary": ["CD16"],
    },
    "Endothelial": {
        "primary": ["CD31"],
        "secondary": [],
    },
    "Epithelial": {
        "primary": ["PanCK"],
        "secondary": [],
    },
    "Fibroblasts": {
        "primary": ["alphaSMA", "Vimentin"],
        "secondary": [],
    },
}

# 10 → 7 collapse mapping for ground truth evaluation
GT_TO_ACHIEVABLE_7_MAPPING: Dict[str, str] = {
    "B cells": "B cells",
    "Mixed Immune": "CD4+ T cells",
    "CD8+ T cells": "CD8+ T cells",
    "Proliferating T": "CD8+ T cells",
    "Macrophages": "Macrophages",
    "Endothelial": "Endothelial",
    "Vascular Stromal": "Endothelial",
    "Epithelial": "Epithelial",
    "Myofibroblasts": "Fibroblasts",
    "Stromal": "Fibroblasts",
}

# Critical markers that MUST be flagged as interesting in Module 1
CRITICAL_MARKERS: List[str] = [
    "CD3E", "CD4", "CD8A",  # T cells
    "CD68", "CD163",  # Macrophages
    "CD20",  # B cells
    "PanCK",  # Epithelial
    "CD31",  # Endothelial
    "alphaSMA", "Vimentin",  # Fibroblasts
]

# Expected colocalization pairs for Module 2a validation
EXPECTED_POSITIVE_PAIRS = [
    ("CD3E", "CD8A"),  # T cell markers
    ("CD68", "CD163"),  # Macrophage markers
    ("CD68", "HLA-DR"),  # Macrophage markers
    ("CD20", "CD45RA"),  # B cell markers
]

EXPECTED_NEGATIVE_PAIRS = [
    ("CD3E", "CD68"),  # T cells vs Macrophages
    ("CD20", "CD68"),  # B cells vs Macrophages
    ("PanCK", "CD20"),  # Epithelial vs B cells
]
```

**Step 2: Commit**

```bash
git add Benchmarking/xenium_benchmarking/benchmark_constants.py
git commit -m "feat: add shared benchmark constants module with achievable-7 GT"
```

---

### Task 2: Update evaluate_pipeline_stages.py to Use Shared Constants

**Files:**
- Modify: `Benchmarking/xenium_benchmarking/CITEgeist/src/evaluate_pipeline_stages.py`

**Step 1: Replace imports and constants**

Replace lines 54-133 (the old `GT_CELLTYPE_MARKERS`, `CRITICAL_MARKERS`, `EXPECTED_POSITIVE_PAIRS`, `EXPECTED_NEGATIVE_PAIRS`, `ORACLE_PROFILES` definitions) with imports from the shared module:

```python
# Import shared benchmark constants
BENCHMARK_ROOT = Path(__file__).parent.parent.parent.parent
sys.path.insert(0, str(BENCHMARK_ROOT))
from benchmark_constants import (
    ACHIEVABLE_7_GT_MARKERS,
    ACHIEVABLE_7_MARKER_SIGNATURES,
    CRITICAL_MARKERS,
    EXPECTED_POSITIVE_PAIRS,
    EXPECTED_NEGATIVE_PAIRS,
    GT_TO_ACHIEVABLE_7_MAPPING,
)
```

Keep `ORACLE_PROFILES` locally since it's specific to oracle evaluation and uses the 10-type granularity for that purpose.

**Step 2: Replace GT_CELLTYPE_MARKERS references**

There are 3 locations where `GT_CELLTYPE_MARKERS` is used:

1. **Line 455** (stage 2b profile scoring loop):
   ```python
   # OLD:
   for ct, markers in GT_CELLTYPE_MARKERS.items():
   # NEW:
   for ct, markers in ACHIEVABLE_7_GT_MARKERS.items():
   ```

2. **Line 626** (stage 2c profile scoring loop):
   ```python
   # OLD:
   for ct, markers in GT_CELLTYPE_MARKERS.items():
   # NEW:
   for ct, markers in ACHIEVABLE_7_GT_MARKERS.items():
   ```

3. **Lines 507, 510, 651, 652** (GT coverage calculations):
   ```python
   # OLD:
   gt_coverage = len(matched_gt_types) / len(GT_CELLTYPE_MARKERS)
   missing_gt = [ct for ct in GT_CELLTYPE_MARKERS if ct not in matched_gt_types]
   # NEW:
   gt_coverage = len(matched_gt_types) / len(ACHIEVABLE_7_GT_MARKERS)
   missing_gt = [ct for ct in ACHIEVABLE_7_GT_MARKERS if ct not in matched_gt_types]
   ```

**Step 3: Commit**

```bash
git add Benchmarking/xenium_benchmarking/CITEgeist/src/evaluate_pipeline_stages.py
git commit -m "fix: use achievable-7 GT instead of old 10-type GT in staged evaluation"
```

---

### Task 3: Update run_benchmark.py to Use Shared Constants

**Files:**
- Modify: `Benchmarking/xenium_benchmarking/CITEgeist/src/run_benchmark.py`

**Step 1: Replace local constant definitions with imports**

Replace lines 48-114 (the `ACHIEVABLE_7_CELL_PROFILE_DICT`, `GT_TO_ACHIEVABLE_7_MAPPING`, `ACHIEVABLE_7_MARKER_SIGNATURES` definitions) with imports:

```python
# Import shared benchmark constants
BENCHMARK_ROOT = Path(__file__).parent.parent.parent
sys.path.insert(0, str(BENCHMARK_ROOT))
from benchmark_constants import (
    ACHIEVABLE_7_CELL_PROFILE_DICT,
    ACHIEVABLE_7_MARKER_SIGNATURES,
    GT_TO_ACHIEVABLE_7_MAPPING,
)
```

Keep the `compress_to_achievable_7()` and `build_cell_profile_dict_from_hierarchical()` functions locally — they use the imported constants but contain logic specific to the benchmark runner.

**Step 2: Verify references still work**

The functions `compress_to_achievable_7()` (line 122) and `build_cell_profile_dict_from_hierarchical()` (line 158) reference `ACHIEVABLE_7_MARKER_SIGNATURES` — confirm these names match the imports.

**Step 3: Commit**

```bash
git add Benchmarking/xenium_benchmarking/CITEgeist/src/run_benchmark.py
git commit -m "refactor: use shared benchmark constants in run_benchmark.py"
```

---

### Task 4: Increase Default Permutation Count

**Files:**
- Modify: `CITEgeist/model/spatial_colocalization.py:753`
- Modify: `Benchmarking/xenium_benchmarking/CITEgeist/src/evaluate_pipeline_stages.py:296`
- Modify: `Benchmarking/xenium_benchmarking/CITEgeist/src/run_benchmark.py:371`

**Step 1: Update default in spatial_colocalization.py**

At line 753, change:
```python
# OLD:
n_permutations: int = 199,
# NEW:
n_permutations: int = 999,
```

**Step 2: Update explicit calls in evaluate_pipeline_stages.py**

At line 296, change:
```python
# OLD:
n_permutations=199,
# NEW:
n_permutations=999,
```

**Step 3: Update explicit calls in run_benchmark.py**

At line 371, change:
```python
# OLD:
n_permutations=199,
# NEW:
n_permutations=999,
```

**Step 4: Commit**

```bash
git add CITEgeist/model/spatial_colocalization.py
git add Benchmarking/xenium_benchmarking/CITEgeist/src/evaluate_pipeline_stages.py
git add Benchmarking/xenium_benchmarking/CITEgeist/src/run_benchmark.py
git commit -m "feat: increase default permutations from 199 to 999 for better FDR resolution"
```

---

### Task 5: Add Per-Spot GMM Signal Masks to Module 1

**Files:**
- Modify: `CITEgeist/model/marker_interest.py:25-38` (MarkerInterest dataclass)
- Modify: `CITEgeist/model/marker_interest.py:115-174` (_fit_gmm_per_marker function)
- Modify: `CITEgeist/model/marker_interest.py:610-637` (identify_interesting_markers result building)

**Step 1: Update _fit_gmm_per_marker to return per-spot signal masks**

At `marker_interest.py:115`, the function currently returns `(snr_values, signal_fractions)`. Change it to also return per-spot signal masks:

```python
def _fit_gmm_per_marker(
    X: NDArray[np.floating],
    seed: int,
) -> Tuple[NDArray[np.floating], NDArray[np.floating], NDArray[np.bool_]]:
    """
    Fit 2-component GMM to each marker to separate signal from noise.

    Args:
        X: Expression matrix (n_spots, n_markers).
        seed: Random seed for GMM initialization.

    Returns:
        Tuple of:
        - snr_values: SNR = (mu_signal - mu_background) / sigma_background (n_markers,)
        - signal_fractions: Fraction of spots in signal component (n_markers,)
        - signal_masks: Boolean mask of signal spots per marker (n_spots, n_markers)
    """
    n_spots, n_markers = X.shape
    snr_values = np.zeros(n_markers)
    signal_fractions = np.zeros(n_markers)
    signal_masks = np.zeros((n_spots, n_markers), dtype=bool)

    for m in range(n_markers):
        values = X[:, m].reshape(-1, 1)

        if np.std(values) < 1e-10:
            snr_values[m] = 0.0
            signal_fractions[m] = 0.0
            continue

        try:
            gmm = GaussianMixture(
                n_components=2,
                covariance_type="full",
                random_state=seed,
                n_init=3,
                max_iter=100,
            )
            gmm.fit(values)

            means = gmm.means_.flatten()
            stds = np.sqrt(gmm.covariances_.flatten())
            weights = gmm.weights_

            if means[0] > means[1]:
                signal_idx, bg_idx = 0, 1
            else:
                signal_idx, bg_idx = 1, 0

            mu_signal = means[signal_idx]
            mu_bg = means[bg_idx]
            sigma_bg = max(stds[bg_idx], 1e-6)

            snr_values[m] = (mu_signal - mu_bg) / sigma_bg
            signal_fractions[m] = weights[signal_idx]

            # Per-spot signal assignment
            labels = gmm.predict(values)
            signal_masks[:, m] = (labels == signal_idx)

        except Exception as e:
            logging.debug(f"GMM fitting failed for marker {m}: {e}")
            snr_values[m] = 0.0
            signal_fractions[m] = 0.0

    return snr_values, signal_fractions, signal_masks
```

**Step 2: Update identify_interesting_markers to store signal_masks**

At line 574, the call currently is:
```python
snr_values, signal_fractions = _fit_gmm_per_marker(X, seed)
```

Change to:
```python
snr_values, signal_fractions, signal_masks = _fit_gmm_per_marker(X, seed)
```

**Step 3: Add signal_masks to MarkerInterestResult**

At lines 40-48, add a new field:

```python
@dataclass
class MarkerInterestResult:
    """Results container for marker interest analysis."""

    markers: List[MarkerInterest]
    kurtosis_threshold: float
    morans_threshold: float
    morans_k: int
    morans_alpha: float
    signal_masks: Optional[NDArray[np.bool_]] = None  # (n_spots, n_markers) per Module 1 order
    signal_mask_marker_names: Optional[List[str]] = None  # Marker names corresponding to signal_masks columns
```

**Step 4: Set signal_masks in the result construction**

At line 642, where `MarkerInterestResult` is created, add the signal_masks:

```python
result = MarkerInterestResult(
    markers=markers,
    kurtosis_threshold=kurtosis_thresh_learned,
    morans_threshold=morans_thresh_learned,
    morans_k=morans_k,
    morans_alpha=0.05,
    signal_masks=signal_masks,
    signal_mask_marker_names=list(marker_names),
)
```

**Step 5: Commit**

```bash
git add CITEgeist/model/marker_interest.py
git commit -m "feat: add per-spot GMM signal masks to MarkerInterestResult"
```

---

### Task 6: Add Singleton Rescue to Module 2c

**Files:**
- Modify: `CITEgeist/model/spatial_colocalization.py` (add rescue function, integrate into select_profiles or as standalone)

**Step 1: Add the rescue function**

Add this function before `select_profiles()` (around line 3395):

```python
def rescue_singletons(
    profiles: List[List[str]],
    signal_masks: NDArray[np.bool_],
    signal_mask_marker_names: List[str],
    min_unique_coverage: float = 0.3,
    min_signal_fraction: float = 0.05,
    verbose: bool = False,
) -> List[List[str]]:
    """
    Filter singletons by unique spatial coverage.

    Multi-marker profiles (2+ markers) are always kept. Singletons are kept
    only if they cover unique spatial territory not explained by multi-marker
    profiles.

    Uses Module 1 GMM signal masks to define "on" spots per marker (spots
    classified as signal vs background by 2-component GMM).

    Args:
        profiles: Candidate profiles from Module 2b.
        signal_masks: Boolean mask (n_spots, n_markers) from Module 1 GMM.
        signal_mask_marker_names: Marker names for signal_masks columns.
        min_unique_coverage: Min fraction of singleton's signal spots that
            are outside multi-marker profile territory (default: 0.3).
        min_signal_fraction: Min fraction of total spots classified as signal
            for this marker (default: 0.05).
        verbose: Log rescue decisions.

    Returns:
        Filtered list of profiles with noise singletons removed.
    """
    n_spots = signal_masks.shape[0]

    # Build marker name → column index mapping
    name_to_idx = {name: i for i, name in enumerate(signal_mask_marker_names)}

    # Separate multi-marker profiles from singletons
    multi_marker = [p for p in profiles if len(p) >= 2]
    singletons = [p for p in profiles if len(p) == 1]

    if not singletons:
        return profiles

    # Compute explained territory: union of all signal spots from multi-marker profiles
    explained = np.zeros(n_spots, dtype=bool)
    for profile in multi_marker:
        for marker in profile:
            if marker in name_to_idx:
                explained |= signal_masks[:, name_to_idx[marker]]

    # Evaluate each singleton
    rescued = list(multi_marker)  # Always keep multi-marker profiles
    dropped = []

    for profile in singletons:
        marker = profile[0]
        if marker not in name_to_idx:
            if verbose:
                logger.info(f"  Singleton '{marker}': not in signal masks, dropping")
            dropped.append(marker)
            continue

        idx = name_to_idx[marker]
        signal_spots = signal_masks[:, idx]
        signal_fraction = np.sum(signal_spots) / n_spots

        # Check minimum area
        if signal_fraction < min_signal_fraction:
            if verbose:
                logger.info(
                    f"  Singleton '{marker}': signal_fraction={signal_fraction:.3f} "
                    f"< {min_signal_fraction}, dropping"
                )
            dropped.append(marker)
            continue

        # Check unique coverage
        unique_spots = signal_spots & ~explained
        unique_coverage = np.sum(unique_spots) / max(np.sum(signal_spots), 1)

        if unique_coverage >= min_unique_coverage:
            rescued.append(profile)
            if verbose:
                logger.info(
                    f"  Singleton '{marker}': RESCUED "
                    f"(unique_coverage={unique_coverage:.3f}, "
                    f"signal_fraction={signal_fraction:.3f})"
                )
        else:
            dropped.append(marker)
            if verbose:
                logger.info(
                    f"  Singleton '{marker}': DROPPED "
                    f"(unique_coverage={unique_coverage:.3f} "
                    f"< {min_unique_coverage})"
                )

    if verbose:
        logger.info(
            f"Singleton rescue: {len(rescued) - len(multi_marker)} rescued, "
            f"{len(dropped)} dropped out of {len(singletons)} singletons"
        )

    return rescued
```

**Step 2: Export the function**

Add `rescue_singletons` to the module's imports in `CITEgeist/model/__init__.py` if it has an explicit `__all__` or import list. Check line 31 area of `__init__.py`.

**Step 3: Commit**

```bash
git add CITEgeist/model/spatial_colocalization.py
git commit -m "feat: add singleton rescue based on unique spatial coverage"
```

---

### Task 7: Integrate Singleton Rescue into Staged Evaluation

**Files:**
- Modify: `Benchmarking/xenium_benchmarking/CITEgeist/src/evaluate_pipeline_stages.py`

**Step 1: Import rescue_singletons**

At line 36-40, add to the imports:
```python
from CITEgeist.model.spatial_colocalization import (
    analyze_marker_colocalization,
    discover_profiles,
    select_profiles,
    rescue_singletons,
)
```

**Step 2: Thread Module 1 result through**

The staged evaluation runs Module 1 first, then Module 2a, then 2b, then 2c. Find where the Module 1 `MarkerInterestResult` is stored and ensure it's available when calling Module 2c.

The `evaluate_stage1()` function returns the result. The main orchestration function needs to pass `marker_interest_result` to the stage 2c evaluation function.

**Step 3: Call rescue_singletons between Module 2b and 2c**

After `discover_profiles()` returns (around line 444) and before `select_profiles()` is called (around line 602), add:

```python
# Rescue singletons using Module 1 GMM signal masks
if marker_interest_result.signal_masks is not None:
    profiles_before = len(discovered_profiles)
    discovered_profiles = rescue_singletons(
        profiles=discovered_profiles,
        signal_masks=marker_interest_result.signal_masks,
        signal_mask_marker_names=marker_interest_result.signal_mask_marker_names,
        min_unique_coverage=0.3,
        min_signal_fraction=0.05,
        verbose=True,
    )
    logger.info(f"Singleton rescue: {profiles_before} -> {len(discovered_profiles)} profiles")
```

Note: The exact integration depends on how the staged evaluation passes data between stages. The `marker_interest_result` from stage 1 needs to be available in the stage 2b/2c evaluation function. Check the main orchestration function to see how to thread this through.

**Step 4: Commit**

```bash
git add Benchmarking/xenium_benchmarking/CITEgeist/src/evaluate_pipeline_stages.py
git commit -m "feat: integrate singleton rescue into staged evaluation pipeline"
```

---

### Task 8: Integrate Singleton Rescue into run_benchmark.py

**Files:**
- Modify: `Benchmarking/xenium_benchmarking/CITEgeist/src/run_benchmark.py`

**Step 1: Import rescue_singletons**

Add to imports at line 34-38:
```python
from CITEgeist.model.spatial_colocalization import (
    analyze_marker_colocalization,
    discover_hierarchical_profiles,
    select_profiles,
    rescue_singletons,
)
```

**Step 2: Call rescue after hierarchical profile discovery**

After `hierarchical_result.flat_profiles` is obtained (around line 399), convert to list format and apply rescue:

```python
discovered_profiles = hierarchical_result.flat_profiles
# Convert dict to list of lists for rescue_singletons
profile_lists = [markers for markers in discovered_profiles.values()]

if marker_interest_result.signal_masks is not None:
    profile_lists = rescue_singletons(
        profiles=profile_lists,
        signal_masks=marker_interest_result.signal_masks,
        signal_mask_marker_names=marker_interest_result.signal_mask_marker_names,
        verbose=True,
    )
    # Rebuild dict from rescued profiles
    discovered_profiles = {
        f"Profile_{i}": markers for i, markers in enumerate(profile_lists)
    }
```

Note: `marker_interest_result` is the return value of `identify_interesting_markers()` called earlier in the benchmark script. Ensure it's in scope at this point.

**Step 3: Commit**

```bash
git add Benchmarking/xenium_benchmarking/CITEgeist/src/run_benchmark.py
git commit -m "feat: integrate singleton rescue into benchmark runner"
```

---

### Task 9: Update __init__.py Exports

**Files:**
- Modify: `CITEgeist/model/__init__.py`

**Step 1: Add rescue_singletons to exports**

Find the import of `spatial_colocalization` functions (around line 31) and add `rescue_singletons`:

```python
from .spatial_colocalization import (
    ...,
    rescue_singletons,
)
```

Also add to `__all__` if it exists.

**Step 2: Commit**

```bash
git add CITEgeist/model/__init__.py
git commit -m "chore: export rescue_singletons from model package"
```

---

### Task 10: Validation - Re-run Staged Evaluation

**Files:**
- Modify (if needed): `Benchmarking/xenium_benchmarking/CITEgeist/slurm/run_staged_evaluation.sh`

**Step 1: Verify SLURM script has mail directives**

Read `run_staged_evaluation.sh` and ensure it has:
```bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alc376@pitt.edu
```

**Step 2: Submit the job**

```bash
sbatch Benchmarking/xenium_benchmarking/CITEgeist/slurm/run_staged_evaluation.sh
```

**Step 3: After completion, compare results**

Check the output files in `output_staged_evaluation/` for:
- Singleton count decreased from ~70% to <40%
- GT coverage improved (old scores misattributed CD4+ T cells)
- PanCK and Vimentin retained as singletons
- PTEN and E-Cadherin dropped
- Macrophage and B cell profiles unchanged (no regression)

**Step 4: Commit results summary**

```bash
git add docs/plans/2026-01-29-module2-autodiscovery-fixes-design.md
git commit -m "docs: add autodiscovery fixes design document"
```
