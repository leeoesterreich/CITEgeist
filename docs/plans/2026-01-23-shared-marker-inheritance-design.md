# Shared Marker Inheritance in Hierarchical Clustering

**Date:** 2026-01-23
**Status:** Design Complete
**Branch:** hierarchical_approach

## Problem Statement

The current hierarchical clustering assigns all markers to leaf nodes, leaving internal nodes with empty marker lists. This prevents shared lineage markers (e.g., CD3 for all T cells) from being inherited by descendant profiles.

**Current behavior:**
```
T_cells (markers=[])           ← empty!
├── CD4_T (markers=['CD3', 'CD4'])
└── CD8_T (markers=['CD3', 'CD8'])
```

**Desired behavior:**
```
T_cells (markers=['CD3'])      ← shared marker
├── CD4_T (markers=['CD4'])    ← specific marker
└── CD8_T (markers=['CD8'])    ← specific marker
```

When flattening, CD4+ profile inherits CD3 → `['CD3', 'CD4']`

## Solution: Bivariate Moran's I at Split Time

At each tree split, classify markers as **shared** (assign to internal node) or **specific** (pass to child) based on spatial correlation with both children.

### Algorithm

```
At each internal node with markers M and children L, R:

1. Compute representative signals:
   - left_signal = mean expression of markers in L subtree
   - right_signal = mean expression of markers in R subtree

2. For each marker m in M:
   - I_left = bivariate_morans_i(m, left_signal)
   - I_right = bivariate_morans_i(m, right_signal)

3. Classify marker:
   - ratio = min(I_left, I_right) / max(I_left, I_right)
   - min_I = min(I_left, I_right)

   If ratio > 0.5 AND min_I > 0.2:
       → SHARED: assign to this internal node
   Else if I_left > I_right:
       → SPECIFIC: pass to left child
   Else:
       → SPECIFIC: pass to right child

4. Recurse:
   - Left child receives: left_specific_markers
   - Right child receives: right_specific_markers
   - This node stores: shared_markers
```

### Threshold Logic

Combined threshold ensures markers are genuinely shared:
- **Ratio > 0.5**: Both correlations are similar (not dominated by one side)
- **Min > 0.2**: Both correlations are meaningful (not just two weak values)

### Separate Trees per Lineage

Disconnected components from Module 2a's colocalization graph become separate trees. Immune markers that don't colocalize with stromal markers naturally form separate hierarchies. The shared marker logic applies **within** each tree.

## Implementation Details

### New Helper Function

```python
def _compute_bivariate_morans_i_fast(
    marker_expr: np.ndarray,       # (n_spots,) expression of one marker
    reference_signal: np.ndarray,  # (n_spots,) mean expression of child group
    spatial_weights: np.ndarray,   # (n_spots, n_spots) precomputed spatial weights
) -> float:
    """Fast bivariate Moran's I for marker vs reference signal."""
```

### Modified `_recursive_tree_cut` Signature

```python
def _recursive_tree_cut(
    scipy_node,
    X: NDArray[np.floating],
    marker_names: List[str],
    marker_indices: List[int],      # NEW: which columns this node owns
    spatial_weights: np.ndarray,    # NEW: precomputed W matrix
    coords: np.ndarray,             # NEW: spatial coordinates
    improvement_threshold: float,
    sharing_ratio: float,           # NEW: default 0.5
    sharing_min_I: float,           # NEW: default 0.2
    current_depth: int,
    max_depth: int,
    parent_id: Optional[str],
) -> ProfileTreeNode:
```

### Core Classification Logic

```python
# Inside _recursive_tree_cut, after deciding to split:

# Get marker indices for left and right subtrees (from clustering)
left_indices = _get_leaf_indices(scipy_node.get_left())
right_indices = _get_leaf_indices(scipy_node.get_right())

# Compute representative signals (mean expression)
left_signal = X[:, left_indices].mean(axis=1)   # (n_spots,)
right_signal = X[:, right_indices].mean(axis=1) # (n_spots,)

# Classify each marker this node is responsible for
shared_markers = []
shared_indices = []
left_specific_markers = []
left_specific_indices = []
right_specific_markers = []
right_specific_indices = []

for idx in marker_indices:
    marker_expr = X[:, idx]
    I_left = _compute_bivariate_morans_i_fast(marker_expr, left_signal, spatial_weights)
    I_right = _compute_bivariate_morans_i_fast(marker_expr, right_signal, spatial_weights)

    max_I = max(I_left, I_right)
    min_I = min(I_left, I_right)
    ratio = min_I / max_I if max_I > 1e-10 else 0.0

    if ratio > sharing_ratio and min_I > sharing_min_I:
        # Shared: stays at this internal node
        shared_markers.append(marker_names[idx])
        shared_indices.append(idx)
    elif I_left >= I_right:
        # Specific to left
        left_specific_markers.append(marker_names[idx])
        left_specific_indices.append(idx)
    else:
        # Specific to right
        right_specific_markers.append(marker_names[idx])
        right_specific_indices.append(idx)

# Create node with shared markers (not empty!)
return ProfileTreeNode(
    node_id=node_id,
    markers=shared_markers,  # Internal nodes now have markers
    children=[left_child, right_child],
    parent_id=parent_id,
    depth=current_depth,
)
```

### Flattening with Inheritance

The existing `_collect_path_markers()` function walks from root to leaf, collecting `node.markers` at each step. With internal nodes now having markers, this naturally produces inherited profiles:

```python
def _flatten_tree_to_profiles(tree, ...):
    flat_profiles = {}
    marker_to_profiles = defaultdict(list)

    for leaf in tree.get_leaves():
        # Collect ALL markers from root to this leaf
        path_markers = _collect_path_markers(leaf, tree.root)
        flat_profiles[leaf.node_id] = list(path_markers)

        for marker in path_markers:
            marker_to_profiles[marker].append(leaf.node_id)

    # Shared markers = those appearing in multiple profiles
    shared_markers = {
        m: profiles for m, profiles in marker_to_profiles.items()
        if len(profiles) > 1
    }

    return flat_profiles, shared_markers
```

## Files to Modify

| File | Changes |
|------|---------|
| `spatial_colocalization.py` | Add `_compute_bivariate_morans_i_fast()`, modify `_recursive_tree_cut()`, modify `_build_hierarchical_tree()`, cleanup `_flatten_tree_to_profiles()` |

## Changes Summary

1. **New helper function:** `_compute_bivariate_morans_i_fast()`
2. **Modified `_recursive_tree_cut()`:** Add marker classification logic, internal nodes store shared markers
3. **Modified `_build_hierarchical_tree()`:** Precompute spatial weights, pass new parameters
4. **Modified `discover_hierarchical_profiles()`:** Accept `coords` parameter, pass threshold parameters
5. **Cleanup `_flatten_tree_to_profiles()`:** Remove failed NMF redistribution code, use simple path collection

## Estimated Scope

~150 lines of new/modified code

## Testing Strategy

1. **Unit test:** Manually constructed tree with known shared/specific markers
2. **Integration test:** Simulated T-cell hierarchy (CD3 shared, CD4/CD8 specific)
3. **Benchmark:** Re-run Xenium and simulated benchmarks to verify distinct profiles with proper inheritance

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `sharing_ratio` | 0.5 | Min ratio of bivariate Moran's I values for marker to be shared |
| `sharing_min_I` | 0.2 | Min bivariate Moran's I with both children for marker to be shared |
