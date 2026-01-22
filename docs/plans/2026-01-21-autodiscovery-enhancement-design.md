# Autodiscovery Enhancement: Correlation Fallback + Reconstruction Merging

**Date:** 2026-01-21
**Status:** Design approved, ready for implementation
**Branch:** `hierarchical_approach`

## Problem Statement

Module 2 autodiscovery over-fragments into single-marker profiles instead of discovering biologically meaningful multi-marker cell types.

**Current behavior (Region 0):**
- 23 discovered profiles, mostly single-marker
- alphaSMA and Vimentin in separate profiles (should be Fibroblasts)
- PanCK and E-Cadherin in separate profiles (should be Epithelial)
- CD3E and CD8A in separate profiles (should be CD8+ T cells)

**Root cause:** Markers that should cluster together lack significant bivariate Moran's I edges, so they end up in separate connected components and can never be merged by the current algorithm.

**Target:** Autodiscovery should produce profiles that map cleanly to the 7 achievable cell types without using biological priors.

## Solution: Two-Stage Enhancement

### Stage 1: Expression Correlation Fallback (Module 2a)

Add secondary linkage based on expression correlation to connect markers that are co-expressed but lack strong spatial colocalization signal.

**Algorithm:**

```python
def compute_correlation_edges(
    adata: AnnData,
    markers: List[str],
    correlation_threshold: float = 0.5,
    co_occurrence_threshold: float = 0.20,
) -> List[Tuple[str, str, float]]:
    """
    Find marker pairs with high expression correlation.

    Returns edges for pairs where:
    1. Pearson r >= correlation_threshold
    2. Both markers expressed (>0) in >= co_occurrence_threshold of spots
    """
    correlation_edges = []

    for i, m1 in enumerate(markers):
        for m2 in markers[i+1:]:
            expr1 = adata[:, m1].X.toarray().flatten()
            expr2 = adata[:, m2].X.toarray().flatten()

            # Check co-occurrence
            both_expressed = (expr1 > 0) & (expr2 > 0)
            co_occurrence_rate = both_expressed.sum() / len(expr1)

            if co_occurrence_rate >= co_occurrence_threshold:
                # Compute correlation only where both expressed
                r, _ = pearsonr(expr1[both_expressed], expr2[both_expressed])

                if r >= correlation_threshold:
                    correlation_edges.append((m1, m2, r))

    return correlation_edges
```

**Integration:**
- Moran's I significant edges: weight = 1.0 (primary)
- Correlation edges: weight = 0.5 (secondary)
- Both edge types contribute to connected components
- During hierarchical clustering, use combined distance matrix

### Stage 2: Reconstruction-Guided Merging (Module 2d)

After initial discovery, iteratively merge profile pairs if reconstruction improves or doesn't significantly degrade.

**Algorithm:**

```python
def merge_profiles_by_reconstruction(
    adata: AnnData,
    profiles: List[List[str]],
    merge_tolerance: float = 0.05,
    max_iterations: int = 10,
) -> List[List[str]]:
    """
    Iteratively merge profile pairs if reconstruction improves.

    Args:
        adata: AnnData with antibody expression
        profiles: Initial profiles from Module 2c
        merge_tolerance: Accept merge if error increases by at most this fraction
        max_iterations: Maximum merge rounds

    Returns:
        Merged profiles
    """
    current_profiles = [list(p) for p in profiles]

    for iteration in range(max_iterations):
        best_merge = None
        best_improvement = -float('inf')

        # Try all pairwise merges
        for i in range(len(current_profiles)):
            for j in range(i + 1, len(current_profiles)):
                p1, p2 = current_profiles[i], current_profiles[j]

                # Compute reconstruction error: separate vs merged
                error_separate = compute_reconstruction_error(adata, [p1, p2])
                error_merged = compute_reconstruction_error(adata, [p1 + p2])

                # Improvement = how much better (or not worse) merged is
                improvement = error_separate - error_merged

                # Accept if merged is better OR within tolerance
                if error_merged <= error_separate * (1 + merge_tolerance):
                    if improvement > best_improvement:
                        best_improvement = improvement
                        best_merge = (i, j)

        if best_merge is None:
            break  # No beneficial merges found

        # Apply best merge
        i, j = best_merge
        merged = current_profiles[i] + current_profiles[j]
        current_profiles = [p for k, p in enumerate(current_profiles) if k not in (i, j)]
        current_profiles.append(merged)

        logging.info(f"Iteration {iteration}: Merged profiles {i}+{j} -> {merged}")

    return current_profiles


def compute_reconstruction_error(
    adata: AnnData,
    profiles: List[List[str]],
) -> float:
    """
    Compute how well profiles reconstruct the antibody expression.

    Uses non-negative least squares to find best coefficients,
    returns mean squared error.
    """
    # Build profile matrix: each column is mean expression of profile markers
    profile_matrix = []
    for profile in profiles:
        profile_expr = adata[:, profile].X.toarray().mean(axis=1)
        profile_matrix.append(profile_expr)

    P = np.column_stack(profile_matrix)  # (n_spots, n_profiles)

    # For each marker, find best reconstruction from profiles
    total_error = 0
    all_markers = set(m for p in profiles for m in p)

    for marker in all_markers:
        y = adata[:, marker].X.toarray().flatten()
        # Non-negative least squares
        coeffs, _ = nnls(P, y)
        y_pred = P @ coeffs
        total_error += np.mean((y - y_pred) ** 2)

    return total_error / len(all_markers)
```

## Pipeline Integration

```python
def run_autodiscovery_enhanced(adata, ...):
    # Module 1: Marker interest (unchanged)
    interesting_markers = identify_interesting_markers(adata)

    # Module 2a: Colocalization + correlation edges (ENHANCED)
    coloc_result = analyze_marker_colocalization(
        adata,
        interesting_markers,
        include_correlation_edges=True,  # NEW
        correlation_threshold=0.5,        # NEW
        co_occurrence_threshold=0.20,     # NEW
    )

    # Module 2b: Profile discovery (uses enhanced graph)
    discovery_result = discover_profiles(coloc_result, ...)

    # Module 2c: Profile selection (unchanged)
    selected_profiles = select_profiles(discovery_result, adata, ...)

    # Module 2d: Reconstruction-guided merging (NEW)
    final_profiles = merge_profiles_by_reconstruction(
        adata,
        selected_profiles,
        merge_tolerance=0.05,
        max_iterations=10,
    )

    return final_profiles
```

## Parameters

| Parameter | Default | Location | Rationale |
|-----------|---------|----------|-----------|
| `include_correlation_edges` | True | Module 2a | Enable correlation fallback |
| `correlation_threshold` | 0.5 | Module 2a | Moderate correlation, not too strict |
| `co_occurrence_threshold` | 0.20 | Module 2a | At least 20% of spots have both markers |
| `merge_tolerance` | 0.05 | Module 2d | Allow merges that increase error by up to 5% |
| `max_iterations` | 10 | Module 2d | Prevent runaway merging |

## Expected Outcomes

**Profile consolidation:**

| Profile | Before (fragmented) | After (expected) |
|---------|---------------------|------------------|
| Fibroblasts | alphaSMA, Vimentin separate | alphaSMA + Vimentin merged |
| Epithelial | PanCK, E-Cadherin separate | PanCK + E-Cadherin merged |
| CD8+ T cells | CD8A alone | CD3E + CD8A merged |
| Macrophages | CD68+CD16, CD163 separate | CD68 + CD163 + CD16 merged |
| B cells | CD20 + CD45RA | CD20 + CD45RA (unchanged) |
| Endothelial | CD31 | CD31 (unchanged) |
| CD4+ T cells | CD3E + CD45RO | CD3E + CD4 + CD45RO merged |

**Success criteria:**
1. Autodiscovery produces ≤10 profiles (down from 23)
2. At least 5 of 7 achievable types are discoverable
3. JSD gap between autodiscovery and manual profiles < 0.10

## File Changes Required

### 1. `CITEgeist/model/spatial_colocalization.py`

- Add `compute_correlation_edges()` function
- Modify `analyze_marker_colocalization()` to include correlation edges
- Add `merge_profiles_by_reconstruction()` function (Module 2d)
- Add `compute_reconstruction_error()` helper function

### 2. `CITEgeist/model/citegeist_model.py`

- Update autodiscovery flow to call Module 2d after Module 2c

### 3. `Benchmarking/xenium_benchmarking/CITEgeist/src/run_benchmark.py`

- Add new parameters for enhanced autodiscovery
- Update `--use-autodiscovery` to use enhanced pipeline

### 4. Tests

- Add `tests/test_correlation_edges.py`
- Add `tests/test_reconstruction_merging.py`

## Implementation Checklist

- [ ] Add `compute_correlation_edges()` to spatial_colocalization.py
- [ ] Modify `analyze_marker_colocalization()` to accept correlation edge parameters
- [ ] Update `_find_connected_components()` to handle weighted edges
- [ ] Add `compute_reconstruction_error()` helper
- [ ] Add `merge_profiles_by_reconstruction()` (Module 2d)
- [ ] Update `discover_profiles()` to return intermediate results for Module 2d
- [ ] Update run_benchmark.py with new parameters
- [ ] Add unit tests
- [ ] Run on Xenium region 0 and verify profile consolidation
- [ ] Benchmark against 7-type achievable GT

## Risks & Mitigations

| Risk | Likelihood | Mitigation |
|------|------------|------------|
| Over-merging (too few profiles) | Medium | Tolerance parameter, max iterations limit |
| Correlation edges add noise | Low | Co-occurrence threshold filters spurious correlations |
| Slow reconstruction computation | Low | NNLS is fast; cache profile matrices |
| No improvement | Low | Fallback to current algorithm if metrics don't improve |
