# Post-Hoc Hierarchy Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Replace the tree-first hierarchical profile discovery with flat-first discovery + post-hoc hierarchy construction, preserving profile quality while adding organizational tree structure.

**Architecture:** Stage 1 calls `discover_profiles_continuous()` (proven 86-100% GT coverage) to get flat profiles. Stage 2 clusters those profiles by marker-set Jaccard similarity into a `ProfileTree` for reporting. The function signature and return type of `discover_hierarchical_profiles_continuous()` are unchanged, so no downstream code needs modification.

**Tech Stack:** scipy (linkage, to_tree), numpy, existing CITEgeist dataclasses

---

### Task 1: Add `_build_posthoc_hierarchy()` helper

**Files:**
- Modify: `CITEgeist/model/spatial_colocalization.py` (add after `_filter_coloc_result_for_markers` at ~line 4690)

**Step 1: Write the helper function**

Add the following function:

```python
def _build_posthoc_hierarchy(
    flat_profiles: Dict[str, List[str]],
) -> Tuple[ProfileTree, Dict[str, List[str]]]:
    """
    Build a post-hoc hierarchy over pre-computed flat profiles.

    Clusters profiles by Jaccard similarity of their marker sets, producing
    a ProfileTree for visualization/reporting.  The flat profiles themselves
    are NOT modified — this is purely organizational.

    Args:
        flat_profiles: Dict mapping profile_name -> list of marker names.

    Returns:
        Tuple of (ProfileTree, shared_markers_dict).
        shared_markers_dict maps marker -> list of profile names that contain it.
    """
    from scipy.cluster.hierarchy import linkage, to_tree
    from scipy.spatial.distance import squareform

    profile_names = list(flat_profiles.keys())
    profile_sets = [set(markers) for markers in flat_profiles.values()]
    n = len(profile_names)

    # Identify shared markers (appear in 2+ profiles)
    marker_to_profiles: Dict[str, List[str]] = defaultdict(list)
    for name, markers in flat_profiles.items():
        for m in markers:
            marker_to_profiles[m].append(name)
    shared_markers = {m: ps for m, ps in marker_to_profiles.items() if len(ps) > 1}

    # Trivial cases
    if n <= 1:
        if n == 1:
            node = ProfileTreeNode(
                node_id=profile_names[0],
                markers=flat_profiles[profile_names[0]],
                children=[], parent_id=None, depth=0,
            )
        else:
            node = ProfileTreeNode(
                node_id="root", markers=[], children=[],
                parent_id=None, depth=0,
            )
        return ProfileTree(root=node, n_levels=1), shared_markers

    # Pairwise Jaccard distance
    D = np.ones((n, n))
    np.fill_diagonal(D, 0.0)
    for i in range(n):
        for j in range(i + 1, n):
            intersection = len(profile_sets[i] & profile_sets[j])
            union = len(profile_sets[i] | profile_sets[j])
            jaccard = intersection / union if union > 0 else 0.0
            D[i, j] = 1.0 - jaccard
            D[j, i] = 1.0 - jaccard

    # Hierarchical clustering
    condensed = squareform(D)
    Z = linkage(condensed, method='average')
    scipy_root = to_tree(Z)

    # Convert scipy tree to ProfileTree
    def _convert(node, depth=0, parent_id=None):
        if node.is_leaf():
            idx = node.id
            return ProfileTreeNode(
                node_id=profile_names[idx],
                markers=flat_profiles[profile_names[idx]],
                children=[], parent_id=parent_id, depth=depth,
            )
        node_id = f"internal_{node.id}"
        left = _convert(node.get_left(), depth + 1, node_id)
        right = _convert(node.get_right(), depth + 1, node_id)
        return ProfileTreeNode(
            node_id=node_id, markers=[], children=[left, right],
            parent_id=parent_id, depth=depth,
        )

    root = _convert(scipy_root)
    n_levels = _compute_tree_depth(root)
    return ProfileTree(root=root, n_levels=n_levels), shared_markers
```

**Step 2: Commit**

```bash
git add CITEgeist/model/spatial_colocalization.py
git commit -m "feat: add _build_posthoc_hierarchy() for flat-first approach"
```

---

### Task 2: Rewrite `discover_hierarchical_profiles_continuous()`

**Files:**
- Modify: `CITEgeist/model/spatial_colocalization.py` (replace function body at lines 4322-4668)

**Step 1: Replace the function**

Replace the entire body of `discover_hierarchical_profiles_continuous()` with:

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
    Discover hierarchical profiles using flat-first, hierarchy-second approach.

    Stage 1: Runs discover_profiles_continuous() to get high-quality flat
    profiles (proven 86-100% GT coverage on Xenium benchmarks).

    Stage 2: Builds a post-hoc hierarchy over the flat profiles by clustering
    them by Jaccard similarity of marker sets. The tree is purely
    organizational — it does not modify the profiles.

    The function signature is unchanged from the tree-first version so that
    callers (run_benchmark.py) require no modification.

    Args:
        coloc_result: Result from analyze_marker_colocalization()
        antibody_expression: Expression matrix (n_spots, n_markers)
        marker_names: Names for each column in expression matrix
        coords: Spatial coordinates (n_spots, 2) [unused, kept for API compat]
        improvement_threshold: Unused, kept for API compatibility.
        sharing_ratio: Unused, kept for API compatibility.
        sharing_min_I: Unused, kept for API compatibility.
        max_depth: Unused, kept for API compatibility.
        neighbor_k: Unused, kept for API compatibility.
        top_k: Number of top partners for mutual top-k masking (default 3).
        distance_metric: Which metric for distance (default 'colocalization_score').
        verbose: Log progress (default True)

    Returns:
        HierarchicalProfileResult with flat profiles and post-hoc tree.
    """
    if verbose:
        logger.info("Hierarchical profile discovery (flat-first, hierarchy post-hoc)...")

    # ── Stage 1: Flat profile discovery ───────────────────────────────────
    flat_result = discover_profiles_continuous(
        colocalization_result=coloc_result,
        top_k=top_k,
        distance_metric=distance_metric,
        verbose=verbose,
    )

    # Convert to dict format expected by HierarchicalProfileResult
    flat_profiles: Dict[str, List[str]] = {}
    for i, profile_markers in enumerate(flat_result.profiles):
        flat_profiles[f"profile_{i}"] = list(profile_markers)
    for i, singleton in enumerate(flat_result.singletons):
        flat_profiles[f"singleton_{i}"] = [singleton]

    if verbose:
        logger.info(f"Stage 1: {len(flat_result.profiles)} multi-marker profiles, "
                     f"{len(flat_result.singletons)} singletons")

    # ── Stage 2: Post-hoc hierarchy ───────────────────────────────────────
    tree, shared_markers = _build_posthoc_hierarchy(flat_profiles)

    if verbose:
        logger.info(f"Stage 2: tree depth {tree.n_levels}, "
                     f"{len(shared_markers)} shared markers")

    # Compute reconstruction error
    final_error = _compute_final_reconstruction_error(
        antibody_expression, marker_names, list(flat_profiles.values())
    )

    if verbose:
        logger.info(f"Discovered {len(flat_profiles)} profiles, "
                     f"reconstruction error: {final_error:.4f}")

    return HierarchicalProfileResult(
        tree=tree,
        flat_profiles=flat_profiles,
        depth_per_branch={},
        shared_markers=shared_markers,
        reconstruction_error=final_error,
    )
```

**Step 2: Commit**

```bash
git add CITEgeist/model/spatial_colocalization.py
git commit -m "feat: rewrite hierarchical discovery as flat-first + posthoc hierarchy"
```

---

### Task 3: Update `run_benchmark.py` caller

**Files:**
- Modify: `Benchmarking/xenium_benchmarking/CITEgeist/src/run_benchmark.py` (lines 320-332)

**Step 1: Simplify the call**

The old call passes `improvement_threshold`, `sharing_ratio`, `sharing_min_I`, `max_depth`, `top_k`. These are kept for API compatibility but unused. Simplify to only pass meaningful params:

```python
    hierarchical_result = discover_hierarchical_profiles_continuous(
        coloc_result=coloc_result,
        antibody_expression=X_protein,
        marker_names=marker_names,
        coords=coords,
        top_k=top_k,
        distance_metric="colocalization_score",
        verbose=True,
    )
```

Also remove the now-unused `improvement_threshold` and `max_depth` from `run_hierarchical_benchmark()` function signature (lines 246-250) and the argparse parameters that feed them.

**Step 2: Commit**

```bash
git add Benchmarking/xenium_benchmarking/CITEgeist/src/run_benchmark.py
git commit -m "refactor: simplify hierarchical benchmark caller for flat-first approach"
```

---

### Task 4: Validate on Xenium benchmarking data

**Files:**
- Run: `Benchmarking/xenium_benchmarking/CITEgeist/slurm/run_benchmark.sh`

**Step 1: Submit SLURM job**

```bash
cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/CITEgeist/slurm
MODE=hierarchical sbatch run_benchmark.sh
```

**Step 2: Monitor and verify**

Expected results per region:
- 9-15 profiles (matching staged eval quality)
- 86-100% GT type coverage (5-7 of 7 types)
- Tree depth and shared markers logged correctly
- Deconvolution completes for all 5 regions

Check logs with:
```bash
grep -aE "Stage [12]|Discovered|After rescue|→" slurm_log/citegeist_bench_*.err
```

**Step 3: Commit validation results**

```bash
git commit -m "validate: flat-first hierarchical discovery on Xenium benchmarks"
```
