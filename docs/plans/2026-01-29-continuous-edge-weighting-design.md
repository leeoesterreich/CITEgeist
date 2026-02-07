# Continuous Edge Weighting for Profile Discovery - Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Replace the binary FDR edge filter in Module 2b profile discovery with continuous colocalization-score-weighted hierarchical clustering over all markers, eliminating sensitivity to permutation count.

**Architecture:** Instead of filtering edges by permutation-based FDR (binary gate), build a complete distance matrix from all marker pairs and cluster directly. The dendrogram naturally separates noise (high distance) from signal (low distance) via gap-based cutting. Singletons emerge as markers with uniformly high distance to all others, detected by the existing `_split_dendrogram_by_gaps` and `_dynamic_tree_cut` logic. The mutual top-k sparsification is retained as an optional noise-reduction step on the distance matrix (setting non-top-k distances to max), but FDR is removed entirely from the clustering path.

**Tech Stack:** Python 3.10, scipy (linkage, fcluster), numpy, existing CITEgeist spatial_colocalization module.

---

## Context for Implementer

### What we're changing and why

Module 2b discovers cell type marker profiles by:
1. Computing bivariate Moran's I between all marker pairs (Module 2a)
2. Filtering pairs by permutation-based BH-FDR correction → binary edge set
3. Building a graph from significant edges → connected components
4. Hierarchical clustering within each component → profiles

The problem: step 2 is sensitive to the permutation count. Changing from 199→999 permutations shifts which edges survive FDR, which changes connected components, which changes profiles. Marginal pairs (p≈0.05) flip in/out, causing CD4/CD8 T cells to merge or separate depending on permutation count.

The fix: skip step 2 entirely. Build the distance matrix from ALL pairs (step 3 already uses continuous distances internally). Let the dendrogram's gap-based cutting handle noise separation — that's what dendrograms are designed for.

### Two discovery functions affected

There are two profile discovery entry points:
1. **`discover_profiles()`** (line 2831) — used by `evaluate_pipeline_stages.py`
2. **`discover_hierarchical_profiles()`** (line 3758) — used by `run_benchmark.py`

Both start with FDR → top-k → connected components. Both need the same change.

### What stays the same

- Module 2a (`analyze_marker_colocalization()`) — still computes all pair statistics including permutation p-values (these are still useful for reporting/diagnostics)
- `_build_distance_matrix()` (line 2362) — uses `colocalization_score`, not p-values
- `_build_colocalization_distance_matrix()` (line 1137) — uses `bivariate_morans_i`, not p-values
- `_dynamic_tree_cut()` (line 2565) — gap-based cutting, unchanged
- `_split_dendrogram_by_gaps()` (line 2465) — gap-based lineage splitting, unchanged
- Singleton rescue (downstream) — unchanged

### Key files

- `CITEgeist/model/spatial_colocalization.py` — core changes
- `Benchmarking/xenium_benchmarking/CITEgeist/src/evaluate_pipeline_stages.py` — caller updates
- `Benchmarking/xenium_benchmarking/CITEgeist/src/run_benchmark.py` — caller updates
- `tests/test_module2_profile_discovery.py` — test updates

---

## Task 1: Add `discover_profiles_continuous()` function

This is the core change. We add a new function that replaces the FDR→components→per-component-clustering pipeline with a single all-markers dendrogram approach.

**Files:**
- Modify: `CITEgeist/model/spatial_colocalization.py` (add new function after `discover_profiles` at ~line 3128)

**Step 1: Write the function**

Add this function after the existing `discover_profiles()` (after line 3128):

```python
def discover_profiles_continuous(
    colocalization_result: ColocalizationResult,
    *,
    top_k: int = 5,
    distance_metric: str = "colocalization_score",
    seed: int = 1234,
    verbose: bool = True,
) -> ProfileDiscoveryResult:
    """
    Discover cell type profiles using continuous edge weighting (no FDR gate).

    Instead of binary FDR filtering, builds a complete distance matrix from all
    marker pairs and uses hierarchical clustering with gap-based cutting to
    separate lineages and profiles. This eliminates sensitivity to permutation
    count since no p-value thresholds are used.

    Algorithm:
    1. Build complete distance matrix from colocalization scores (all pairs)
    2. Optional: mutual top-k masking (set non-top-k distances to max)
    3. Hierarchical clustering (average linkage) over all markers
    4. Gap-based lineage splitting (large gaps = separate lineages)
    5. Dynamic tree cutting within each lineage to extract profiles
    6. Markers with uniformly high distance become singletons naturally

    Args:
        colocalization_result: Result from analyze_marker_colocalization().
        top_k: Number of top partners to keep per marker in mutual top-k
            distance masking (default: 5). Set to 0 to disable.
            Unlike FDR-based filtering, this only sets non-top-k distances
            to maximum (1.0) rather than removing edges entirely.
        distance_metric: Which metric to use for distance:
            - 'colocalization_score' (default): Combined score (0.3*spearman +
              0.3*cosine + 0.4*bivariate_morans), normalized to [0,1].
              Distance = 1 - score.
            - 'bivariate_morans': Uses bivariate Moran's I only.
              Distance = 1 - normalized_I.
        seed: Random seed for reproducibility (default: 1234).
        verbose: Log progress information (default: True).

    Returns:
        ProfileDiscoveryResult with discovered profiles, dendrograms, and metrics.
    """
    pairs = colocalization_result.pairs
    all_markers = colocalization_result.marker_names
    marker_specificity = colocalization_result.marker_specificity
    n_markers = len(all_markers)

    if verbose:
        logging.info(f"Discovering profiles (continuous weighting) from {len(pairs)} marker pairs...")
        logging.info(f"Distance metric: {distance_metric}, top_k: {top_k}")

    if n_markers < 2:
        return ProfileDiscoveryResult(
            profiles=[[m] for m in all_markers],
            lineage_dendrograms={},
            singletons=list(all_markers),
            modularity=0.0,
            n_significant_edges=0,
            alpha=0.0,
        )

    # Step 1: Build complete distance matrix from ALL pairs
    marker_to_idx = {m: i for i, m in enumerate(all_markers)}

    if distance_metric == "bivariate_morans":
        # Use bivariate Moran's I with min-max normalization
        morans_values = [p.bivariate_morans_i for p in pairs]
        if len(morans_values) > 0:
            min_i = min(morans_values)
            max_i = max(morans_values)
            range_i = max_i - min_i if max_i > min_i else 1.0
        else:
            min_i, range_i = 0.0, 1.0

        dist_matrix = np.ones((n_markers, n_markers))
        np.fill_diagonal(dist_matrix, 0.0)

        for pair in pairs:
            i = marker_to_idx[pair.marker_a]
            j = marker_to_idx[pair.marker_b]
            normalized = (pair.bivariate_morans_i - min_i) / range_i
            d = 1.0 - normalized
            dist_matrix[i, j] = d
            dist_matrix[j, i] = d
    else:
        # Use colocalization_score (already in [0, 1])
        dist_matrix = np.ones((n_markers, n_markers))
        np.fill_diagonal(dist_matrix, 0.0)

        for pair in pairs:
            i = marker_to_idx[pair.marker_a]
            j = marker_to_idx[pair.marker_b]
            d = 1.0 - pair.colocalization_score
            dist_matrix[i, j] = d
            dist_matrix[j, i] = d

    if verbose:
        # Report distance statistics
        upper_tri = dist_matrix[np.triu_indices(n_markers, k=1)]
        logging.info(
            f"Distance matrix: min={upper_tri.min():.3f}, "
            f"median={np.median(upper_tri):.3f}, max={upper_tri.max():.3f}"
        )

    # Step 2: Optional mutual top-k distance masking
    # Instead of removing edges, set non-top-k distances to max (1.0)
    # This preserves the continuous nature while reducing noise
    if top_k > 0 and len(pairs) > 0:
        # For each marker, find its top-k closest partners
        top_k_mask = np.zeros((n_markers, n_markers), dtype=bool)

        for i in range(n_markers):
            # Get distances to all other markers (excluding self)
            dists = dist_matrix[i, :].copy()
            dists[i] = np.inf  # Exclude self

            # Apply specificity weighting if available
            if marker_specificity:
                marker_name = all_markers[i]
                s_i = marker_specificity.get(marker_name, 1.0)
                for j in range(n_markers):
                    if j != i:
                        s_j = marker_specificity.get(all_markers[j], 1.0)
                        # Lower specificity → higher effective distance
                        dists[j] = dists[j] / (np.sqrt(s_i * s_j) + 1e-10)

            # Find top-k closest
            k_actual = min(top_k, n_markers - 1)
            top_k_indices = np.argsort(dists)[:k_actual]
            for j in top_k_indices:
                top_k_mask[i, j] = True

        # Mutual: both must rank each other in top-k
        mutual_mask = top_k_mask & top_k_mask.T

        # Set non-mutual-top-k distances to max
        n_masked = 0
        for i in range(n_markers):
            for j in range(i + 1, n_markers):
                if not mutual_mask[i, j]:
                    dist_matrix[i, j] = 1.0
                    dist_matrix[j, i] = 1.0
                    n_masked += 1

        n_total_pairs = n_markers * (n_markers - 1) // 2
        n_kept = n_total_pairs - n_masked
        if verbose:
            logging.info(
                f"Mutual top-{top_k} masking: {n_kept}/{n_total_pairs} pairs kept "
                f"(specificity-weighted: {marker_specificity is not None})"
            )

    # Step 3: Hierarchical clustering over all markers
    condensed_dist = squareform(dist_matrix)
    Z = linkage(condensed_dist, method='average')

    if verbose:
        merge_distances = Z[:, 2]
        logging.info(
            f"Linkage: {len(Z)} merges, "
            f"distances [{merge_distances.min():.3f} ... {merge_distances.max():.3f}]"
        )

    # Step 4: Gap-based lineage splitting
    lineages = _split_dendrogram_by_gaps(Z, all_markers, seed=seed)

    if verbose:
        logging.info(f"Gap analysis: {len(lineages)} lineages")
        for i, lin in enumerate(lineages):
            logging.info(f"  Lineage {i + 1} ({len(lin)} markers): {lin}")

    # Step 5: Dynamic tree cutting within each lineage
    profiles = []
    lineage_dendrograms = {}
    lineage_idx = 0

    for lineage_markers in lineages:
        n_lineage = len(lineage_markers)

        if n_lineage == 1:
            profiles.append(lineage_markers)
            continue

        if n_lineage == 2:
            # Check distance between the pair
            i = marker_to_idx[lineage_markers[0]]
            j = marker_to_idx[lineage_markers[1]]
            if dist_matrix[i, j] < 0.7:  # Reasonable colocalization
                profiles.append(lineage_markers)
            else:
                # Too distant — separate singletons
                profiles.append([lineage_markers[0]])
                profiles.append([lineage_markers[1]])
            continue

        # Build lineage-specific dendrogram
        lin_marker_indices = [marker_to_idx[m] for m in lineage_markers]
        lin_dist = dist_matrix[np.ix_(lin_marker_indices, lin_marker_indices)]
        lin_condensed = squareform(lin_dist)
        lin_Z = linkage(lin_condensed, method='average')

        # Store dendrogram
        lineage_dendrograms[lineage_idx] = LineageDendrogram(
            markers=lineage_markers,
            linkage_matrix=lin_Z,
            distance_matrix=lin_dist,
        )
        lineage_idx += 1

        # Dynamic tree cutting
        cluster_labels = _dynamic_tree_cut(lin_Z, n_lineage)

        cluster_to_markers = defaultdict(list)
        for i, label in enumerate(cluster_labels):
            cluster_to_markers[label].append(lineage_markers[i])

        for cluster_markers in cluster_to_markers.values():
            profiles.append(cluster_markers)

    # Identify singletons (profiles with 1 marker)
    singletons = [p[0] for p in profiles if len(p) == 1]

    # Compute modularity using all pairs (not just FDR-filtered)
    # Use pairs with colocalization_score > 0.5 as "significant" for modularity calc
    significant_pairs = [p for p in pairs if p.colocalization_score > 0.5]
    modularity = _compute_modularity(profiles, pairs, significant_pairs)

    # Count edges with meaningful colocalization for reporting
    n_meaningful_edges = sum(1 for p in pairs if p.colocalization_score > 0.5)

    if verbose:
        logging.info(f"Discovered {len(profiles)} profiles, modularity = {modularity:.3f}")

    result = ProfileDiscoveryResult(
        profiles=profiles,
        lineage_dendrograms=lineage_dendrograms,
        singletons=singletons,
        modularity=modularity,
        n_significant_edges=n_meaningful_edges,
        alpha=0.0,  # No FDR threshold used
    )

    if verbose:
        logging.info(result.summary())

    return result
```

**Step 2: Add export to `__init__.py`**

In `CITEgeist/model/__init__.py`, add `discover_profiles_continuous` to the imports from `spatial_colocalization` and to `__all__`.

**Step 3: Commit**

```bash
git add CITEgeist/model/spatial_colocalization.py CITEgeist/model/__init__.py
git commit -m "feat: add discover_profiles_continuous() for permutation-insensitive profile discovery"
```

---

## Task 2: Add `discover_hierarchical_profiles_continuous()` function

The benchmark runner uses `discover_hierarchical_profiles()` which has additional logic for shared marker detection and reconstruction-guided tree cutting. We need an equivalent continuous version.

**Files:**
- Modify: `CITEgeist/model/spatial_colocalization.py` (add new function after `discover_hierarchical_profiles` at ~line 3960)

**Step 1: Write the function**

The key change is replacing the FDR→top-k→adaptive-ratio→connected-components pipeline (lines 3823-3903) with a continuous distance matrix approach. The rest of the function (tree building with shared marker detection) stays the same.

```python
def discover_hierarchical_profiles_continuous(
    coloc_result: ColocalizationResult,
    antibody_expression: NDArray[np.floating],
    marker_names: List[str],
    coords: NDArray[np.floating],
    improvement_threshold: float = 0.05,
    sharing_ratio: float = 0.5,
    sharing_min_I: float = 0.2,
    max_depth: int = 5,
    neighbor_k: int = 6,
    top_k: int = 3,
    distance_metric: str = "colocalization_score",
    verbose: bool = True,
) -> HierarchicalProfileResult:
    """
    Discover hierarchical profiles using continuous edge weighting (no FDR gate).

    Same as discover_hierarchical_profiles() but replaces the FDR→top-k→
    adaptive-ratio→connected-components pipeline with continuous distance
    matrix clustering. This eliminates sensitivity to permutation count.

    Algorithm:
    1. Build complete distance matrix from all colocalization scores
    2. Optional mutual top-k distance masking
    3. Hierarchical clustering over all markers
    4. Gap-based splitting into lineages (replaces connected components)
    5. Per-lineage reconstruction-guided tree building with shared marker detection

    Args:
        coloc_result: Result from analyze_marker_colocalization()
        antibody_expression: Expression matrix (n_spots, n_markers)
        marker_names: Names for each column in expression matrix
        coords: Spatial coordinates (n_spots, 2)
        improvement_threshold: Min reconstruction improvement to split (default 5%)
        sharing_ratio: Min ratio of bivariate Moran's I for shared classification
        sharing_min_I: Min bivariate Moran's I with both children for shared
        max_depth: Maximum tree depth (default 5)
        neighbor_k: Number of neighbors for spatial weights (default 6)
        top_k: Mutual top-k for distance masking (default 3). Set to 0 to disable.
        distance_metric: 'colocalization_score' (default) or 'bivariate_morans'
        verbose: Log progress (default True)

    Returns:
        HierarchicalProfileResult with tree structure and flat profiles.
    """
    if verbose:
        logger.info("Starting hierarchical profile discovery (continuous weighting)...")

    pairs = coloc_result.pairs
    all_markers = coloc_result.marker_names
    n_markers = len(all_markers)
    marker_to_idx_all = {m: i for i, m in enumerate(all_markers)}

    # Phase 1: Build complete distance matrix and find lineages
    if verbose:
        logger.info("Phase 1: Building continuous distance matrix...")

    # Build distance matrix
    if distance_metric == "bivariate_morans":
        morans_values = [p.bivariate_morans_i for p in pairs]
        min_i = min(morans_values) if morans_values else 0.0
        max_i = max(morans_values) if morans_values else 1.0
        range_i = max_i - min_i if max_i > min_i else 1.0

        dist_matrix = np.ones((n_markers, n_markers))
        np.fill_diagonal(dist_matrix, 0.0)
        for pair in pairs:
            i = marker_to_idx_all[pair.marker_a]
            j = marker_to_idx_all[pair.marker_b]
            normalized = (pair.bivariate_morans_i - min_i) / range_i
            d = 1.0 - normalized
            dist_matrix[i, j] = d
            dist_matrix[j, i] = d
    else:
        dist_matrix = np.ones((n_markers, n_markers))
        np.fill_diagonal(dist_matrix, 0.0)
        for pair in pairs:
            i = marker_to_idx_all[pair.marker_a]
            j = marker_to_idx_all[pair.marker_b]
            d = 1.0 - pair.colocalization_score
            dist_matrix[i, j] = d
            dist_matrix[j, i] = d

    # Optional mutual top-k masking
    if top_k > 0:
        specificity = coloc_result.marker_specificity
        top_k_mask = np.zeros((n_markers, n_markers), dtype=bool)

        for i in range(n_markers):
            dists = dist_matrix[i, :].copy()
            dists[i] = np.inf
            if specificity:
                s_i = specificity.get(all_markers[i], 1.0)
                for j in range(n_markers):
                    if j != i:
                        s_j = specificity.get(all_markers[j], 1.0)
                        dists[j] = dists[j] / (np.sqrt(s_i * s_j) + 1e-10)
            k_actual = min(top_k, n_markers - 1)
            top_k_indices = np.argsort(dists)[:k_actual]
            for j in top_k_indices:
                top_k_mask[i, j] = True

        mutual_mask = top_k_mask & top_k_mask.T
        for i in range(n_markers):
            for j in range(i + 1, n_markers):
                if not mutual_mask[i, j]:
                    dist_matrix[i, j] = 1.0
                    dist_matrix[j, i] = 1.0

        if verbose:
            n_total = n_markers * (n_markers - 1) // 2
            n_kept = np.sum(mutual_mask[np.triu_indices(n_markers, k=1)])
            logger.info(f"Mutual top-{top_k} masking: {n_kept}/{n_total} pairs kept")

    # Hierarchical clustering over all markers
    condensed_dist = squareform(dist_matrix)
    Z = linkage(condensed_dist, method='average')

    # Gap-based lineage splitting (replaces connected components)
    lineages = _split_dendrogram_by_gaps(Z, all_markers, seed=1234)

    if verbose:
        logger.info(f"Gap analysis: {len(lineages)} lineages")

    # Separate singletons: markers in 1-marker lineages with high min-distance
    singletons = []
    multi_marker_lineages = []
    for lineage in lineages:
        if len(lineage) == 1:
            marker = lineage[0]
            idx = marker_to_idx_all[marker]
            # Check if this marker is truly isolated (high distance to everything)
            other_dists = [dist_matrix[idx, j] for j in range(n_markers) if j != idx]
            min_dist = min(other_dists) if other_dists else 1.0
            if min_dist > 0.7:  # Isolated marker
                singletons.append(marker)
            else:
                multi_marker_lineages.append(lineage)
        else:
            multi_marker_lineages.append(lineage)

    if verbose:
        logger.info(f"Found {len(multi_marker_lineages)} multi-marker lineages, "
                    f"{len(singletons)} singletons")

    # Phase 2: Build hierarchical tree for each lineage
    # (Same as discover_hierarchical_profiles from here on)
    all_flat_profiles: Dict[str, List[str]] = {}
    all_shared_markers: Dict[str, List[str]] = defaultdict(list)
    all_trees: List[ProfileTree] = []
    profile_idx = 0

    for comp_idx, lineage_markers in enumerate(multi_marker_lineages):
        comp_markers = sorted(lineage_markers)
        n_comp = len(comp_markers)

        if verbose:
            logger.info(f"Processing lineage {comp_idx + 1} ({n_comp} markers): {comp_markers}")

        if n_comp == 1:
            profile_name = f"profile_{profile_idx}"
            all_flat_profiles[profile_name] = comp_markers
            profile_idx += 1
            continue

        if n_comp == 2:
            profile_name = f"profile_{profile_idx}"
            all_flat_profiles[profile_name] = comp_markers
            profile_idx += 1
            continue

        # Build distance matrix for this lineage only
        comp_coloc = _filter_coloc_result_for_markers(coloc_result, comp_markers)
        D, ordered_markers = _build_colocalization_distance_matrix(comp_coloc)

        # Build hierarchical clustering
        condensed = squareform(D)
        condensed = np.maximum(condensed, 0)
        lin_Z = linkage(condensed, method='ward')

        # Get marker indices
        marker_to_idx = {m: i for i, m in enumerate(marker_names)}
        comp_marker_indices = [marker_to_idx[m] for m in ordered_markers]

        # Build tree with reconstruction-guided cutting
        tree = _build_hierarchical_tree(
            X=antibody_expression,
            marker_names=marker_names,
            linkage_matrix=lin_Z,
            coords=coords,
            improvement_threshold=improvement_threshold,
            sharing_ratio=sharing_ratio,
            sharing_min_I=sharing_min_I,
            max_depth=max_depth,
            neighbor_k=neighbor_k,
            coloc_result=comp_coloc,
            ordered_markers=ordered_markers,
            marker_indices=comp_marker_indices,
            verbose=verbose,
        )

        all_trees.append(tree)

        # Flatten tree to profiles
        flat = tree.to_flat_profiles()
        for name, markers_list in flat.items():
            profile_name = f"profile_{profile_idx}"
            all_flat_profiles[profile_name] = markers_list
            profile_idx += 1

        # Track shared markers
        shared = tree.get_shared_markers()
        for marker_name, shared_with in shared.items():
            all_shared_markers[marker_name].extend(shared_with)

    # Add singletons
    for singleton in singletons:
        profile_name = f"profile_{profile_idx}"
        all_flat_profiles[profile_name] = [singleton]
        profile_idx += 1

    # Build combined tree (for visualization)
    if len(all_trees) == 1:
        combined_tree = all_trees[0]
    elif len(all_trees) > 1:
        combined_tree = ProfileTree(
            markers=list({m for t in all_trees for m in t.markers}),
            left=all_trees[0] if len(all_trees) > 0 else None,
            right=all_trees[1] if len(all_trees) > 1 else None,
            shared_markers=[],
            split_distance=1.0,
            reconstruction_improvement=0.0,
        )
    else:
        combined_tree = ProfileTree(
            markers=singletons,
            left=None,
            right=None,
            shared_markers=[],
            split_distance=0.0,
            reconstruction_improvement=0.0,
        )

    result = HierarchicalProfileResult(
        tree=combined_tree,
        flat_profiles=all_flat_profiles,
        shared_markers=dict(all_shared_markers),
        n_profiles=len(all_flat_profiles),
        n_singletons=len(singletons),
    )

    if verbose:
        logger.info(result.summary())

    return result
```

**Step 2: Add export to `__init__.py`**

Add `discover_hierarchical_profiles_continuous` to the imports and `__all__` in `CITEgeist/model/__init__.py`.

**Step 3: Commit**

```bash
git add CITEgeist/model/spatial_colocalization.py CITEgeist/model/__init__.py
git commit -m "feat: add discover_hierarchical_profiles_continuous() for benchmark runner"
```

---

## Task 3: Wire up `evaluate_pipeline_stages.py` to use continuous discovery

**Files:**
- Modify: `Benchmarking/xenium_benchmarking/CITEgeist/src/evaluate_pipeline_stages.py`

**Step 1: Update import**

At the top where `discover_profiles` is imported (around line 28), add `discover_profiles_continuous`:

```python
from CITEgeist.model import (
    identify_interesting_markers,
    analyze_marker_colocalization,
    discover_profiles,
    discover_profiles_continuous,  # ADD THIS
    rescue_singletons,
)
```

**Step 2: Replace `discover_profiles` call with `discover_profiles_continuous`**

In `run_stage2b()` (around line 382), change:

```python
    # BEFORE:
    profile_result = discover_profiles(
        colocalization_result=coloc_result,
        fdr_alpha=0.05,
        top_k=6,
        use_triangle_assembly=False,
        pvalue_source="bivariate_morans",
        verbose=True,
    )
```

To:

```python
    # AFTER:
    profile_result = discover_profiles_continuous(
        colocalization_result=coloc_result,
        top_k=6,
        distance_metric="colocalization_score",
        verbose=True,
    )
```

Note: `fdr_alpha`, `use_triangle_assembly`, and `pvalue_source` parameters are no longer needed.

**Step 3: Commit**

```bash
git add Benchmarking/xenium_benchmarking/CITEgeist/src/evaluate_pipeline_stages.py
git commit -m "feat: use continuous profile discovery in staged evaluation"
```

---

## Task 4: Wire up `run_benchmark.py` to use continuous discovery

**Files:**
- Modify: `Benchmarking/xenium_benchmarking/CITEgeist/src/run_benchmark.py`

**Step 1: Update import**

Add `discover_hierarchical_profiles_continuous` to imports:

```python
from CITEgeist.model import (
    identify_interesting_markers,
    analyze_marker_colocalization,
    discover_hierarchical_profiles,
    discover_hierarchical_profiles_continuous,  # ADD THIS
    rescue_singletons,
)
```

**Step 2: Replace `discover_hierarchical_profiles` call**

In the Module 2b section (around line 320), change:

```python
    # BEFORE:
    hierarchical_result = discover_hierarchical_profiles(
        coloc_result=coloc_result,
        antibody_expression=X_protein,
        marker_names=marker_names,
        coords=coords,
        improvement_threshold=improvement_threshold,
        sharing_ratio=0.5,
        sharing_min_I=0.2,
        max_depth=max_depth,
        fdr_alpha=fdr_alpha,
        top_k=top_k,
        verbose=True,
    )
```

To:

```python
    # AFTER:
    hierarchical_result = discover_hierarchical_profiles_continuous(
        coloc_result=coloc_result,
        antibody_expression=X_protein,
        marker_names=marker_names,
        coords=coords,
        improvement_threshold=improvement_threshold,
        sharing_ratio=0.5,
        sharing_min_I=0.2,
        max_depth=max_depth,
        top_k=top_k,
        distance_metric="colocalization_score",
        verbose=True,
    )
```

Note: `fdr_alpha` parameter is no longer passed.

**Step 3: Commit**

```bash
git add Benchmarking/xenium_benchmarking/CITEgeist/src/run_benchmark.py
git commit -m "feat: use continuous hierarchical discovery in benchmark runner"
```

---

## Task 5: Validate on Xenium benchmarking data

Submit the staged evaluation SLURM job to verify continuous discovery produces equivalent or better results than the FDR-gated version.

**Files:**
- Run: `Benchmarking/xenium_benchmarking/CITEgeist/slurm/run_staged_evaluation.sh`

**Step 1: Submit the job**

```bash
sbatch Benchmarking/xenium_benchmarking/CITEgeist/slurm/run_staged_evaluation.sh
```

**Step 2: Monitor and verify**

Spawn a monitoring agent to watch the SLURM job. When complete, check:

1. All 5 regions complete without errors
2. Stage 2b GT coverage is >= 5/7 (71%) in all regions (at least matching previous)
3. CD8+ T cells appear in separate profiles in >= 3/5 regions
4. Fibroblasts (alphaSMA) appear in profiles in >= 4/5 regions
5. Auto-discovery outperforms oracle in >= 3/5 regions

**Step 3: Compare to previous results (job 7855800)**

Previous results for reference:

| Region | Stage 2b GT | Stage 2c GT | Auto Pearson | Oracle Pearson |
|--------|------------|------------|-------------|---------------|
| 0 | 7/7 | 6/7 | 0.551 | 0.426 |
| 1 | 7/7 | 6/7 | 0.395 | 0.393 |
| 2 | 5/7 | 5/7 | 0.540 | 0.501 |
| 3 | 7/7 | 6/7 | 0.555 | 0.488 |
| 4 | 6/7 | 6/7 | 0.531 | 0.480 |

The continuous version should produce consistent results regardless of permutation count, with comparable or better GT coverage and deconvolution accuracy.

**Step 4: Commit results**

If validation passes, no code changes needed. If issues found, iterate on Tasks 1-4.

---

## Key Design Decisions

### Why keep `top_k` masking?

The mutual top-k masking serves a different purpose than FDR: it prevents "hub markers" (like CD45, which colocalize with everything) from pulling unrelated lineages into the same component. Without top-k, CD45 would bridge immune and epithelial markers into one giant dendrogram. Top-k masking is deterministic — it only depends on the colocalization scores, not permutations.

### Why `distance_metric="colocalization_score"` as default?

The `colocalization_score` combines Spearman (0.3), cosine similarity (0.3), and bivariate Moran's I (0.4). Using the combined score is more robust than Moran's I alone because it captures both non-spatial correlation and spatial structure. The combined score is already in [0,1] and doesn't need min-max normalization.

### Why keep `_split_dendrogram_by_gaps` for lineage separation?

Connected components were previously responsible for separating unrelated lineages (immune vs epithelial). Without FDR, all markers are in one dendrogram. `_split_dendrogram_by_gaps` handles this — it looks for unusually large gaps in merge distances, which naturally occur between biologically unrelated marker groups. The gap between immune markers (Moran's I ≈ 0.3-0.7) and epithelial markers (Moran's I ≈ 0.0) creates a clear distance gap.

### What about the singleton isolation threshold (0.7)?

For 2-marker lineages in `discover_profiles_continuous`, we use `dist < 0.7` to decide if they should be a pair or two singletons. This corresponds to `colocalization_score > 0.3`, which is a low bar — any detectable colocalization keeps them together. For singleton detection in `discover_hierarchical_profiles_continuous`, `min_dist > 0.7` identifies truly isolated markers (no meaningful colocalization with any other marker).
