# Hierarchical Profile Discovery for Module 2

**Date:** 2026-01-22
**Status:** Design approved, ready for implementation
**Branch:** `hierarchical_approach`

## Problem Statement

Current Module 2 autodiscovery produces fragmented single-marker profiles because it treats markers as belonging to exactly one cell type. In reality, markers form hierarchies:

- CD3 is shared by all T cells (CD4+, CD8+, Tregs)
- Vimentin is shared by Fibroblasts, Myofibroblasts, EMT cells
- CD45 spans all immune lineages

**Current behavior:** Over-fragmentation into 20+ single-marker profiles, with related markers (CD3, CD4, CD8) in separate profiles.

**Target:** Autodiscovery should learn hierarchical marker relationships and produce profiles where shared markers (CD3) appear in multiple related profiles (CD4+ T cells, CD8+ T cells).

## Solution: Hierarchical NMF-based Profile Discovery

A two-phase algorithm that:
1. Learns a tree structure from spatial colocalization patterns
2. Assigns markers to appropriate hierarchy levels (shared vs specific)
3. Flattens to standard profiles with proper weight allocation for Module 3

**Scope:** Module 2 enhancement only. Module 3 deconvolution unchanged - it receives better-structured flat profiles.

## Algorithm Overview

### Phase 1: Structure Learning (Hierarchical Clustering)

**Input:** Interesting markers from Module 1, spatial colocalization matrix from Module 2a

**Step 1: Build colocalization distance matrix**
```
D[i,j] = 1 - normalized_bivariate_morans_I(marker_i, marker_j)
```
High Moran's I (spatially colocalized) → low distance → same branch

**Step 2: Hierarchical clustering**
- Agglomerative clustering with Ward's linkage on distance matrix D
- Produces full dendrogram from individual markers up to root

**Step 3: Reconstruction-guided tree cutting**
```python
def cut_tree(node, markers, threshold=0.05):
    if node is leaf:
        return [markers]

    left_markers, right_markers = split by clustering

    # Test if split improves reconstruction
    error_merged = reconstruction_error(markers)
    error_split = reconstruction_error(left) + reconstruction_error(right)

    improvement = (error_merged - error_split) / error_merged

    if improvement > threshold:
        # Split is worthwhile, recurse
        return cut_tree(left) + cut_tree(right)
    else:
        # Stop here, keep as single node
        return [markers]
```

**Output:** Tree structure where each node contains a set of markers, with variable depth per branch.

### Phase 2: Weight Learning (NMF at Each Level)

**Input:** Tree structure from Phase 1, marker expression matrix X (spots × markers)

**Step 1: NMF at each internal node**

For each non-leaf node with k children:
```
X_node = X[:, node_markers]  # Expression of markers in this subtree
W, H = NMF(X_node, n_components=k)

# W: (n_spots × k) - spot loadings per child branch
# H: (k × n_markers) - marker weights per child branch
```

The NMF rank k equals the number of child branches (determined by clustering).

**Step 2: Propagate weights down the tree**

Each marker accumulates weights as we traverse from root to leaves:
```
marker_weight[leaf][marker] = product(H values along path from root to leaf)
```

For multi-branch markers (e.g., CD45 in both T-cell and Myeloid branches):
- Marker appears in multiple paths
- Gets separate weights in each branch
- These propagate independently to respective leaves

**Step 3: Expression-weighted allocation for shared markers**

When a marker M belongs to parent node P with children C1, C2:
```
weight_C1 = H[C1, M] / (H[C1, M] + H[C2, M])
weight_C2 = H[C2, M] / (H[C1, M] + H[C2, M])
```

Higher NMF loading → marker contributes more to that child's final profile.

**Output:** For each leaf node (final cell type), a weighted marker profile.

### Phase 3: Flattening to Standard Profiles

**Input:** Hierarchical tree with NMF weights at each level

**Goal:** Convert to flat `Dict[str, List[str]]` format for Module 3

**Step 1: Identify leaf nodes as final cell types**
```python
cell_types = [leaf.name for leaf in tree.leaves()]
# e.g., ["CD4+ T cells", "CD8+ T cells", "Fibroblasts", "Myofibroblasts", ...]
```

**Step 2: Collect markers with weights for each leaf**
```python
for each leaf:
    profile[leaf] = {}
    for each marker in leaf's subtree (including ancestors):
        weight = propagated_weight(marker, leaf)
        if weight > min_weight_threshold:  # e.g., 0.05
            profile[leaf][marker] = weight
```

**Step 3: Convert weighted profiles to binary marker lists**
```python
profile_dict[cell_type] = [m for m, w in profile.items() if w > threshold]
```

Binary output maintains compatibility with existing Module 3.

## Mathematical Formulation

**Notation:**
- X ∈ ℝ^(n×m): Antibody expression matrix (n spots, m markers)
- T: Tree structure with nodes N, leaves L ⊂ N
- For node p with children C(p) = {c₁, c₂, ...}

**NMF at each internal node p:**

Given markers M(p) in node p's subtree:
```
X_p = X[:, M(p)]  ∈ ℝ^(n × |M(p)|)

Minimize: ||X_p - W_p H_p||²_F + λ||H_p||₁

Subject to: W_p ≥ 0, H_p ≥ 0
Where:
  W_p ∈ ℝ^(n × |C(p)|)       # spot loadings for each child
  H_p ∈ ℝ^(|C(p)| × |M(p)|)  # marker weights for each child
```

**Reconstruction error for tree cutting:**
```
E(node) = ||X_node - W_node H_node||²_F / ||X_node||²_F
```

Split node if:
```
[E(parent) - (E(left) + E(right))] / E(parent) > threshold
```

**Weight propagation for marker j to leaf l:**
```
w(j, l) = ∏_{p ∈ path(root, l)} H_p[child_index(p), j] / Σ_k H_p[k, j]
```

Normalizing at each level ensures weights sum appropriately.

**Final profile extraction:**
```
Profile(l) = {j : w(j, l) > min_weight_threshold}
```

## Key Design Decisions

| Decision | Choice | Rationale |
|----------|--------|-----------|
| Hierarchy scope | Module 2 only | Module 3 unchanged, receives better flat profiles |
| Depth selection | Data-driven (reconstruction + parsimony) | Automatically adapts to data complexity |
| Structure learning | Spatial colocalization-guided | Leverages CITEgeist's core strength |
| NMF rank | Constrained by cluster count | Consistency between structure and factorization |
| Shared marker allocation | Expression-weighted | Most biologically accurate |
| Tree cutting | Reconstruction-guided, per-branch | Allows asymmetric depth |
| Orphan markers | Multi-branch membership | Handles pan-lineage markers (CD45) |
| Output format | Binary profiles | Backward compatible with Module 3 |

## Integration with Existing Pipeline

**Current Module 2 pipeline:**
```
Module 2a: analyze_marker_colocalization() → ColocalizationResult
Module 2b: discover_profiles() → ProfileDiscoveryResult
Module 2c: select_profiles_by_reconstruction() → Selected profiles
```

**Enhanced pipeline:**
```
Module 2a: analyze_marker_colocalization() → ColocalizationResult  [unchanged]
Module 2b: discover_hierarchical_profiles() → HierarchicalProfileResult  [NEW]
Module 2c: select_profiles_by_reconstruction() → Selected profiles  [minor update]
```

**New function signature:**
```python
def discover_hierarchical_profiles(
    coloc_result: ColocalizationResult,
    antibody_adata: AnnData,
    improvement_threshold: float = 0.05,
    min_marker_weight: float = 0.05,
    max_depth: int = 5,  # safety limit
) -> HierarchicalProfileResult:
    """
    Phase 1: Build tree from colocalization distances
    Phase 2: NMF weights at each level
    Phase 3: Flatten with expression-weighted allocation
    """
```

**New dataclass:**
```python
@dataclass
class HierarchicalProfileResult:
    tree: ProfileTree  # Full hierarchy for inspection
    flat_profiles: Dict[str, List[str]]  # For Module 3
    depth_per_branch: Dict[str, int]  # How deep each lineage went
    shared_markers: Dict[str, List[str]]  # Which markers are shared
```

## Validation & Success Criteria

### Test 1: Xenium Benchmarking (Hierarchy Detection)

| Lineage | Expected Hierarchy | Success Criterion |
|---------|-------------------|-------------------|
| T cells | Parent: T cells (CD3) → Children: CD4+ (CD4), CD8+ (CD8) | CD3 appears in both CD4+ and CD8+ flat profiles |
| Stromal | Parent: Mesenchymal (Vimentin) → Children: Fibroblasts (FAP?), Myofibroblasts (αSMA) | Vimentin shared, αSMA/FAP distinguish subtypes |
| Myeloid | Parent: Myeloid (CD68?) → Children: Macrophages, etc. | Proper subdivision if markers support it |
| Epithelial | Standalone based on PanCK | PanCK anchors its own lineage (E-cadherin not required) |

**Key validation for T cells:** CD4+ T cells profile should contain [CD3, CD4] and CD8+ T cells profile should contain [CD3, CD8] as the top markers after flattening. This has been previously validated to work well with manual profiles in Module 3.

**Key validation for stromal:** The hierarchy should effectively parse stromal compartments (Fibroblasts, Myofibroblasts, CAFs, EMT cells) which share Vimentin but differ in other markers (αSMA, FAP).

**Metrics:**
- Autodiscovery JSD gap vs manual profiles < 0.10
- Pearson r improvement on CD4+ T cells (currently negative)
- Stromal subtypes correctly separated

### Test 2: Simulated Benchmarking (Adversarial - No Hierarchy)

The simulated data has 9 cell types with 2 specific markers each (18 total markers). No markers are shared across cell types - this is an adversarial test to ensure the algorithm doesn't create spurious hierarchy.

| Aspect | Expected Behavior |
|--------|------------------|
| Detected depth | 1 (flat) |
| Profile count | 9 (matches ground truth) |
| Markers per profile | 2 (no spurious additions) |
| Performance | No regression vs current Module 2 |

The reconstruction-guided cutting should automatically detect that no hierarchy is needed:
- Splitting any node won't improve reconstruction (markers are already specific)
- Algorithm stops at level 1, outputs flat profiles

### Test 3: Edge Cases

| Case | Expected Behavior |
|------|------------------|
| Single-marker profiles | Remain as leaves |
| Pan-immune markers (CD45) | Appear in multiple lineage branches |
| Markers with no colocalization | Isolated nodes, single-marker profiles |

## Implementation Plan

### File Changes

| File | Changes |
|------|---------|
| `CITEgeist/model/spatial_colocalization.py` | Add `discover_hierarchical_profiles()`, `HierarchicalProfileResult`, `ProfileTree` classes |
| `CITEgeist/model/citegeist_model.py` | Update autodiscovery flow to use new function |
| `tests/test_hierarchical_profiles.py` | New test file for hierarchy discovery |
| `Benchmarking/xenium_benchmarking/CITEgeist/src/run_benchmark.py` | Add `--hierarchical-discovery` flag |

### New Dependencies

None - NMF from sklearn is already available in the environment.

### Estimated Scope

- ~300 lines for core algorithm in spatial_colocalization.py
- ~100 lines for dataclasses and tree structure
- ~150 lines for tests

## Parameters

| Parameter | Default | Location | Rationale |
|-----------|---------|----------|-----------|
| `improvement_threshold` | 0.05 | Phase 1 | 5% reconstruction improvement required to split |
| `min_marker_weight` | 0.05 | Phase 3 | Minimum weight to include marker in flat profile |
| `max_depth` | 5 | Phase 1 | Safety limit to prevent runaway recursion |
| `nmf_lambda` | 0.1 | Phase 2 | L1 sparsity penalty for NMF |

## Risks & Mitigations

| Risk | Likelihood | Mitigation |
|------|------------|------------|
| Over-splitting (too many levels) | Low | Reconstruction threshold, max_depth limit |
| Under-splitting (misses hierarchy) | Medium | Tune improvement_threshold, validate on Xenium |
| NMF instability | Low | Use sklearn's robust implementation, multiple restarts |
| Slow on large marker sets | Low | NMF is fast; hierarchy limits problem size at each node |
| Regression on simulated data | Low | Adversarial test explicitly validates flat detection |

## Implementation Checklist

- [ ] Add `ProfileTree` class to represent hierarchy
- [ ] Add `HierarchicalProfileResult` dataclass
- [ ] Implement Phase 1: `_build_colocalization_tree()` with reconstruction-guided cutting
- [ ] Implement Phase 2: `_compute_nmf_weights()` at each tree level
- [ ] Implement Phase 3: `_flatten_to_profiles()` with weight propagation
- [ ] Add main entry point `discover_hierarchical_profiles()`
- [ ] Update `citegeist_model.py` to use new function in autodiscovery flow
- [ ] Add unit tests for hierarchy detection (Xenium-like data)
- [ ] Add unit tests for flat detection (simulated-like data)
- [ ] Add integration test with full pipeline
- [ ] Run Xenium benchmark with hierarchical discovery
- [ ] Run simulated benchmark to verify no regression
- [ ] Validate stromal compartment parsing

---

**Last Updated:** 2026-01-22
**Author:** Claude + Alex
