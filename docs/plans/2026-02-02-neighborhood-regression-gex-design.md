# Neighborhood Regression GEX Deconvolution Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Replace the current enrichment-heuristic GEX deconvolution with a joint neighborhood regression that solves `E[n,k] ≈ Σ_j P[n,j] · G[j,k]` across all neighborhood spots in Gurobi, estimating cell-type-specific expression rates directly.

**Architecture:** The new formulation replaces the per-center-spot QP with a per-neighborhood non-negative least squares (NNLS) regression. For each center spot's neighborhood (~7 spots), we define variables `G[j,k]` representing the expression rate of gene k in cell type j. The objective minimizes reconstruction error `Σ_n Σ_k (E[n,k] - Σ_j P[n,j] · G[j,k])²` across all neighborhood spots simultaneously. The center spot's deconvolved expression is then `X[center,j,k] = P[center,j] · G[j,k]` scaled to conserve observed counts. This replaces all enrichment computation, inverse-frequency weighting, and target heuristics.

**Tech Stack:** gurobipy (QP), numpy, scipy.sparse, scanpy (AnnData)

---

### Task 1: Write the new neighborhood regression solver

**Files:**
- Modify: `CITEgeist/model/gurobi_impl.py:2072-2309` (replace `deconvolute_spot_with_neighbors_with_prior`)

**Step 1: Write the new function**

Replace the body of `deconvolute_spot_with_neighbors_with_prior()` with the neighborhood regression formulation. Keep the same function signature for backward compatibility (unused parameters like `local_enrichment_weight`, `global_enrichment_weight` become no-ops).

The new function does:

```python
def deconvolute_spot_with_neighbors_with_prior(
    spot_idx, adata, cell_type_numbers_array, radius,
    global_prior=None, lambda_prior_weight=0.0,
    local_enrichment_weight=0.5, global_enrichment_weight=0.5,
    continuous_relaxation=True, lambda_gex_reg=0.01,
):
    """
    Deconvolve gene expression using neighborhood regression.

    Solves: min Σ_n Σ_k (E[n,k] - Σ_j P[n,j] · G[j,k])²
    where G[j,k] >= 0 are cell-type-specific expression rates.

    Center spot allocation: X[j,k] = P[center,j] · G[j,k] scaled to
    conserve observed counts E[center,k].
    """
    import gurobipy as gp
    from gurobipy import GRB
    model = None
    try:
        # 1. Get neighborhood
        neighborhood_indices = get_neighbors_with_fixed_radius(
            spot_idx, adata, radius=int(radius), include_center=True
        )
        if not neighborhood_indices:
            return None
        neighborhood_indices = np.array(
            [int(idx) for idx in neighborhood_indices
             if isinstance(idx, (int, np.integer))], dtype=int
        )

        # 2. Extract data
        expr_data = adata.X
        if scipy.sparse.issparse(expr_data):
            expr_data = expr_data.toarray()
        elif not isinstance(expr_data, np.ndarray):
            expr_data = np.array(expr_data)

        K = len(neighborhood_indices)  # number of spots in neighborhood
        T = cell_type_numbers_array.shape[1]  # cell types
        M = expr_data.shape[1]  # genes

        # Neighborhood expression (K x M) and proportions (K x T)
        E = expr_data[neighborhood_indices, :]
        P = cell_type_numbers_array[neighborhood_indices, :]

        # Center spot data (index 0 in neighborhood)
        center_expr = E[0, :]
        center_props = P[0, :]

        # 3. Build Gurobi model
        model = gp.Model(f"gex_neighborhood_{spot_idx}")
        model.setParam("OutputFlag", 0)
        model.setParam("Threads", 1)
        model.setParam("TimeLimit", 30)

        # 4. Variables: G[j,k] >= 0 — expression rate per cell type per gene
        G = {}
        for j in range(T):
            for k in range(M):
                G[j, k] = model.addVar(
                    vtype=GRB.CONTINUOUS, lb=0,
                    name=f"G_{j}_{k}"
                )

        # 5. Objective: min Σ_n Σ_k (E[n,k] - Σ_j P[n,j] · G[j,k])²
        obj_terms = []
        for n in range(K):
            for k in range(M):
                # Reconstruction: Σ_j P[n,j] * G[j,k]
                residual_terms = []
                for j in range(T):
                    if P[n, j] > 0:
                        residual_terms.append(P[n, j] * G[j, k])

                if residual_terms:
                    reconstruction = gp.quicksum(residual_terms)
                else:
                    reconstruction = 0

                observed = float(E[n, k])
                # (observed - reconstruction)^2
                # = observed^2 - 2*observed*reconstruction + reconstruction^2
                # observed^2 is constant, skip
                obj_terms.append(
                    reconstruction * reconstruction
                    - 2 * observed * reconstruction
                )

        # L2 regularization on G
        if lambda_gex_reg > 0:
            for j in range(T):
                for k in range(M):
                    obj_terms.append(lambda_gex_reg * G[j, k] * G[j, k])

        model.setObjective(gp.quicksum(obj_terms), GRB.MINIMIZE)

        # 6. Solve
        model.optimize()

        if model.status == GRB.OPTIMAL:
            # 7. Extract G values and compute center spot allocation
            result = np.zeros((T, M))
            for k in range(M):
                observed_total = float(center_expr[k])
                if observed_total <= 0:
                    continue

                # Raw allocation: P[center,j] * G[j,k]
                raw = np.zeros(T)
                for j in range(T):
                    raw[j] = center_props[j] * G[j, k].X

                # Scale to conserve observed counts
                raw_sum = np.sum(raw)
                if raw_sum > 0:
                    result[:, k] = (raw / raw_sum) * observed_total
                else:
                    # Fallback: uniform allocation
                    result[:, k] = observed_total / T

            return result
        else:
            logging.error(
                f"No feasible solution for spot {spot_idx}, "
                f"status={model.status}"
            )
            return None

    except Exception as e:
        logging.error(f"Error deconvolving spot {spot_idx}: {e}")
        logging.error(traceback.format_exc())
        return None
    finally:
        if model:
            del model
        gc.collect()
```

Key differences from current implementation:
- **Variables**: `G[j,k]` (expression rates) instead of `X[j,k]` (count allocations)
- **No enrichment computation** — removed entirely
- **All neighborhood spots** contribute to the objective (K equations per gene)
- **Count conservation** is post-hoc scaling, not a constraint (the regression finds rates, we scale to counts)
- **No `normalized_weights`**, no `celltype_frequencies`, no `compute_expression_aware_enrichment()`

**Step 2: Run the benchmark on one region to verify it works**

```bash
# Delete stale checkpoints first
cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/CITEgeist/output_protein_gt
rm -f Xenium_region_0_gene_expression_pass1.parquet
rm -rf Xenium_region_0_pass1/ checkpoints/

# Submit single region
sbatch --array=0 /path/to/run_protein_gt_benchmark.sh
```

Expected: Job completes without errors, produces `Xenium_region_0_gene_expression_pass1.parquet`.

**Step 3: Commit**

```bash
git add CITEgeist/model/gurobi_impl.py
git commit -m "feat: replace enrichment heuristic with neighborhood regression QP for GEX deconvolution"
```

---

### Task 2: Run full 5-region benchmark and comparison

**Files:**
- No code changes — benchmark infrastructure already exists

**Step 1: Delete all stale GEX outputs**

```bash
cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/CITEgeist/output_protein_gt
rm -f Xenium_region_*_gene_expression_pass1.parquet
rm -rf Xenium_region_*_pass1/ checkpoints/
```

**Step 2: Submit full benchmark**

```bash
sbatch /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/CITEgeist/slurm/run_protein_gt_benchmark.sh
```

**Step 3: After benchmark completes, submit comparison**

```bash
sbatch /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/evaluation/slurm/run_full_comparison_both.sh
```

**Step 4: Read comparison output and compare to baseline**

Baseline (current QP + enrichment):
- GEX Pearson r: 0.4081
- GEX RMSE: 4.4975
- GEX MAE: 2.7446
- Proportions r: 0.6643 (should be unchanged)

Look at: `evaluation/slurm/slurm_log/compare_both_*.out`

---

### Task 3: Tune lambda_gex_reg if needed

If results improve, sweep the L2 weight to find optimal:

**Files:**
- Modify: `Benchmarking/xenium_benchmarking/CITEgeist/src/sweep_gex_lambda.py` (update to match new function signature if needed)

**Step 1: Run lambda sweep on region 0**

Test lambdas: 0.0, 0.001, 0.01, 0.1, 1.0

**Step 2: Set optimal default in three places**

- `gurobi_impl.py` function signature default
- `gurobi_impl.py` `optimize_gene_expression()` signature default
- `citegeist_model.py` `run_cell_expression_pass1()` signature default

---

## Design Notes

### Why this should work better

The current approach computes a **heuristic enrichment score** by splitting spots into "high expression" vs "low expression" groups and comparing their cell type proportions. This is a weak statistical signal — it uses a median split, averages over groups, and produces a single scalar per cell-type-per-gene.

The new approach uses the **actual expression values and proportions** from all neighborhood spots as direct constraints in the QP. For 7 spots with 7 cell types, each gene has 7 equations and 7 unknowns — a well-determined NNLS system. The solver finds the rates G[j,k] that best explain the observed expression pattern across the whole neighborhood, then allocates the center spot's counts accordingly.

### Model size

Current: T×M variables, M constraints per spot (e.g., 7×400 = 2,800 vars, 400 constraints)
New: T×M variables, 0 hard constraints per spot (e.g., 7×400 = 2,800 vars, QP only)

The number of variables is the same. The objective has K×M quadratic terms instead of M, so ~7× more terms, but these are simple and Gurobi handles QPs efficiently. Runtime should be comparable.

### What's removed

- `compute_expression_aware_enrichment()` — no longer called
- `celltype_frequencies`, `inverse_frequency_weights`, `normalized_weights` — no longer computed
- `gene_specific_enrichment` array — no longer needed
- `center_props` as target input — replaced by regression
- The entire enrichment blending (`local_enrichment_weight`, `global_enrichment_weight`) — parameters kept in signature for compatibility but unused

### Proportions are unchanged

This only affects Pass 2 (GEX). Pass 1 (proportions) is a completely separate code path and will produce identical results.
