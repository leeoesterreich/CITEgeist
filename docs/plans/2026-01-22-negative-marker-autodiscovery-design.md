# Negative Marker Autodiscovery Design

**Date:** 2026-01-22
**Status:** Proposed
**Author:** Claude + Alex

## Problem Statement

CITEgeist's current profile discovery (Module 2) and deconvolution (Module 3) only support positive markers. This causes problems when cell types share markers but need to be distinguished by marker absence:

**Example from Xenium benchmark:**
- EMT cells: Vimentin+ E-Cadherin+ (both high)
- Stromal: Vimentin+ E-Cadherin- (Vimentin high, E-Cadherin low)

**Current result:**
- EMT overpredicted: 35-45% vs 1.6-6.2% actual (10-25x)
- Stromal underpredicted: 1.4% vs 6-9% actual (5x under)

The optimizer cannot distinguish "Vimentin alone" from "Vimentin + E-Cadherin" without negative marker handling.

## Objective

Add unsupervised negative marker discovery to Module 2 and negative marker handling to Module 3, enabling flow cytometry-style gating logic (positive AND negative markers) derived purely from spatial data.

## Design

### Module 2b Extension: Residual Signal Detection

After discovering positive co-localization clusters, detect markers with "residual signal" - markers that appear in contexts beyond their discovered profile.

```python
def detect_residual_signal(
    marker_data: np.ndarray,           # (N_spots, M_markers)
    profiles: Dict[str, List[str]],    # discovered positive profiles
    marker_names: List[str],
    residual_threshold: float = 0.05,  # min fraction of spots with residual
) -> Dict[str, List[str]]:
    """
    Detect markers with significant signal outside their primary profile.

    Returns:
        Dict mapping marker -> list of profiles where it has residual signal
    """
    residual_markers = {}

    for profile_name, profile_markers in profiles.items():
        for marker in profile_markers:
            m_idx = marker_names.index(marker)
            other_markers = [m for m in profile_markers if m != marker]

            if not other_markers:
                continue  # single-marker profile, no residual possible

            # Profile signature: spots where OTHER markers are high
            other_indices = [marker_names.index(m) for m in other_markers]
            other_high = np.all(marker_data[:, other_indices] > threshold, axis=1)

            # Marker high spots
            marker_high = marker_data[:, m_idx] > threshold

            # Residual: marker high BUT other markers low
            residual_spots = marker_high & ~other_high
            residual_fraction = residual_spots.sum() / marker_high.sum()

            if residual_fraction > residual_threshold:
                if marker not in residual_markers:
                    residual_markers[marker] = []
                residual_markers[marker].append(profile_name)

    return residual_markers
```

### Module 2c Extension: Hierarchical Profile Resolution

For markers with residual signal, create new profiles with negative markers to distinguish from the original profile.

```python
def resolve_hierarchical_profiles(
    profiles: Dict[str, List[str]],           # positive-only profiles
    residual_markers: Dict[str, List[str]],   # markers with residual
    marker_data: np.ndarray,
    marker_names: List[str],
) -> Dict[str, Dict[str, List[str]]]:
    """
    Create profiles with positive AND negative markers.

    Returns:
        Dict[profile_name, {"positive": [...], "negative": [...]}]
    """
    resolved_profiles = {}

    for profile_name, profile_markers in profiles.items():
        resolved_profiles[profile_name] = {
            "positive": profile_markers.copy(),
            "negative": []
        }

    for marker, source_profiles in residual_markers.items():
        for source_profile in source_profiles:
            # Create new profile for "marker alone" pattern
            other_markers = [m for m in profiles[source_profile] if m != marker]

            new_profile_name = f"{marker}_only"  # or biologically meaningful name
            resolved_profiles[new_profile_name] = {
                "positive": [marker],
                "negative": other_markers  # exclude the co-localizing markers
            }

    return resolved_profiles
```

### Module 3: Optimization with Negative Markers

**New assignment structure:**

```python
@dataclass
class ProfileAssignment:
    """Assignment of markers to cell types with positive/negative distinction."""
    positive_matrix: np.ndarray  # (M_markers, T_types) - positive assignments
    negative_matrix: np.ndarray  # (M_markers, T_types) - negative assignments
    cell_type_names: List[str]
    marker_names: List[str]
```

**Updated objective function:**

```python
def optimize_cell_proportions_with_negatives(
    marker_data: np.ndarray,           # (N, M)
    assignment: ProfileAssignment,
    lambda_reg: float = 0.01,
    lambda_neg: float = 1.0,           # NEW: negative marker penalty weight
    **kwargs
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Optimize cell proportions with positive and negative marker constraints.

    Objective:
        minimize Σ_i Σ_t Σ_m [
            positive_assignment[m,t] * (S[i,m] - β[m,t] * Y[i,t])²
            + λ_neg * negative_assignment[m,t] * (S[i,m] * Y[i,t])²
        ] + regularization

    The negative term (S[i,m] * Y[i,t])² penalizes high Y when negative
    marker signal is high, but has no effect when marker signal is low.
    """
    N, M = marker_data.shape
    T = len(assignment.cell_type_names)

    model = gp.Model("CellProportions_WithNegatives")

    # Variables
    Y = model.addVars(N, T, lb=0, ub=1, vtype=GRB.CONTINUOUS, name="Y")
    beta = model.addVars(M, T, lb=0.1, ub=2.0, vtype=GRB.CONTINUOUS, name="beta")

    # Objective terms
    error_terms = []

    for i in range(N):
        for t in range(T):
            for m in range(M):
                S_im = marker_data[i, m]

                # Positive marker term
                if assignment.positive_matrix[m, t] > 0:
                    error_terms.append(
                        (S_im - beta[m,t] * Y[i,t]) * (S_im - beta[m,t] * Y[i,t])
                    )

                # Negative marker term
                if assignment.negative_matrix[m, t] > 0:
                    error_terms.append(
                        lambda_neg * S_im * S_im * Y[i,t] * Y[i,t]
                    )

    # Sum constraint
    for i in range(N):
        model.addConstr(gp.quicksum(Y[i,t] for t in range(T)) >= 0.9)
        model.addConstr(gp.quicksum(Y[i,t] for t in range(T)) <= 1.1)

    model.setObjective(gp.quicksum(error_terms), GRB.MINIMIZE)
    model.optimize()

    # Extract results
    Y_values = np.array([[Y[i,t].X for t in range(T)] for i in range(N)])
    beta_values = np.array([[beta[m,t].X for t in range(T)] for m in range(M)])

    return Y_values, beta_values
```

### Mathematical Justification

**For a spot with HIGH Vimentin, LOW E-Cadherin:**

EMT (positive: [Vim, Ecad]):
```
(HIGH - β × Y_EMT)² + (LOW - β × Y_EMT)²
→ Vimentin wants Y_EMT high
→ E-Cadherin wants Y_EMT ≈ 0
→ Compromise: Y_EMT moderate/low
```

Stromal (positive: [Vim], negative: [Ecad]):
```
(HIGH - β × Y_stromal)² + λ_neg × (LOW × Y_stromal)²
→ Vimentin wants Y_stromal high
→ Negative term: LOW × Y_stromal ≈ 0 (minimal penalty)
→ Result: Y_stromal can be HIGH
```

**For a spot with HIGH Vimentin, HIGH E-Cadherin:**

EMT: Both terms want Y_EMT high → Y_EMT high ✓

Stromal:
```
(HIGH - β × Y_stromal)² + λ_neg × (HIGH × Y_stromal)²
→ Vimentin wants Y_stromal high
→ Negative term: BIG penalty if Y_stromal is high
→ Result: Y_stromal pushed DOWN ✓
```

## Output Format

Discovered profiles will use an extended dictionary format:

```python
discovered_profiles = {
    "B cells": {
        "positive": ["CD20", "CD45RA"],
        "negative": []
    },
    "T cells": {
        "positive": ["CD3E"],
        "negative": []
    },
    "EMT cells": {
        "positive": ["Vimentin", "E-Cadherin"],
        "negative": []
    },
    "Stromal": {
        "positive": ["Vimentin"],
        "negative": ["E-Cadherin", "alphaSMA"]
    },
    "Fibroblasts": {
        "positive": ["Vimentin", "alphaSMA"],
        "negative": ["E-Cadherin"]
    },
}
```

## Implementation Plan

### Phase 1: Module 2 Extensions
1. Add `detect_residual_signal()` to `spatial_colocalization.py`
2. Add `resolve_hierarchical_profiles()` to create positive/negative profiles
3. Update `ProfileDiscoveryResult` dataclass to include negative markers
4. Add unit tests with synthetic data

### Phase 2: Module 3 Extensions
1. Create `ProfileAssignment` dataclass with positive/negative matrices
2. Add `optimize_cell_proportions_with_negatives()` to `gurobi_impl.py`
3. Update `CitegeistModel` to use new optimization when negative markers present
4. Add `lambda_neg` parameter with sensible default (1.0)
5. Update validation to handle profiles with same positives but different negatives

### Phase 3: Integration & Testing
1. End-to-end test on Xenium benchmark data
2. Compare results: positive-only vs positive+negative profiles
3. Validate EMT/Stromal separation improves
4. Tune `lambda_neg` if needed

## Success Criteria

1. Module 2 automatically discovers that Vimentin appears in two contexts (with/without E-Cadherin)
2. Module 2 outputs Stromal profile with negative E-Cadherin marker
3. Module 3 correctly separates EMT from Stromal in benchmark
4. EMT prediction drops from 35-45% to ~5% (closer to 1.6-6.2% GT)
5. Stromal prediction increases from 1.4% to ~7% (closer to 6-9% GT)

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `residual_threshold` | 0.05 | Min fraction of marker-high spots with residual signal |
| `lambda_neg` | 1.0 | Weight for negative marker penalty in optimization |
| `negative_marker_pvalue` | 0.05 | Significance threshold for residual detection |

## Risks and Mitigations

1. **Over-fragmentation**: Too many profiles created from residual detection
   - Mitigation: Require minimum residual fraction (5%) and statistical significance

2. **Lambda_neg sensitivity**: Results may vary with penalty weight
   - Mitigation: Grid search on validation data, provide sensible default

3. **Circular dependencies**: Profile A negative for marker M, Profile B positive for M
   - Mitigation: This is expected and correct (EMT+ Ecad, Stromal- Ecad)

## References

- Flow cytometry gating strategies for mesenchymal populations
- CITEgeist Module 2 implementation: `CITEgeist/model/spatial_colocalization.py`
- CITEgeist Module 3 implementation: `CITEgeist/model/gurobi_impl.py`
- Xenium benchmark results: `Benchmarking/xenium_benchmarking/evaluation/validated_benchmark/`
