# Global IQP for Discrete Cell Assignment

**Date:** 2026-02-17
**Status:** Approved
**Author:** Claude + Alex

## Problem Statement

The discrete cell assignment model underperforms on the mixed simulation dataset compared to the continuous model:

| Model | Mixed Prop Corr | High_seg Prop Corr |
|-------|-----------------|---------------------|
| Continuous | ~0.70+ | ~0.76 |
| Discrete (current) | ~0.43 | ~0.76 |

**Root cause:** The discrete model's E-step solves each spot independently via per-spot IQP, while the continuous model solves all spots jointly in a single global optimization. This means:

- Continuous: Information flows between spots through shared beta in global objective
- Discrete: Each spot is an island; locally optimal but globally inconsistent assignments

When markers overlap between cell types (mixed dataset), per-spot IQP can't distinguish cell types. The continuous model uses global context to resolve ambiguity.

**Key insight:** We do NOT assume cells cluster spatially (no Laplacian). We assume **marker behavior per cell type is globally consistent** across all spots.

## Success Criteria

- Mixed dataset proportion correlation improves from ~0.43 to ≥0.65
- High_seg and Xenium datasets show no regression (≥0.76 and ≥0.66 respectively)
- Runtime remains reasonable (<30 min for 5000 spots)

## Design

### Core Concept

Replace per-spot IQP loop with a single global IQP that solves for all cell counts across all spots simultaneously.

**Current (per-spot):**
```
For i = 1 to N (independently):
    minimize: ||X[i] - α - β ⊙ (c[i] @ profile)||²
    s.t.: Σₜ c[i,t] = Nᵢ, c[i,t] ∈ Z≥0
```

**Proposed (global):**
```
minimize: Σᵢ ||X[i] - α - β ⊙ (c[i] @ profile)||²

s.t.: Σₜ c[i,t] = Nᵢ  ∀i    (per-spot nuclei constraint)
      c[i,t] ∈ Z≥0    ∀i,t  (integer cell counts)
```

### Output Semantics

`c_values[i, t]` = number of cells of type `t` in spot `i` (nuclei labeling)

Example for spot i with 10 nuclei:
```
c_values[i, :] = [3, 2, 0, 1, 4, 0, 0, 0, 0]
                  │  │     │  │
                  │  │     │  └─ 4 Cancer_Epithelial cells
                  │  │     └─ 1 Endothelial cell
                  │  └─ 2 CAFs
                  └─ 3 T-cells

Sum = 10 = nuclei_counts[i] ✓
```

Proportions are derived post-hoc for GT comparison only: `prop[i,t] = c[i,t] / N[i]`

### New Function: `solve_discrete_cell_counts_global()`

```python
def solve_discrete_cell_counts_global(
    marker_level_data: np.ndarray,      # (N, M) antibody signal
    marker_names: List[str],
    assignment_matrix: np.ndarray,       # (M, T) marker-to-celltype
    cell_type_names: List[str],
    nuclei_counts: np.ndarray,           # (N,) integer counts per spot
    beta_values: np.ndarray,             # (M,) marker scaling factors
    alpha_values: np.ndarray,            # (M,) marker baselines
    time_limit: float = 300.0,           # 5 min default
    mip_gap: float = 0.05,               # Accept 5% suboptimal
    prev_c_values: Optional[np.ndarray] = None,  # Warm-start
) -> np.ndarray:  # (N, T) integer cell counts
```

### Gurobi Model

**Variables:**
```python
c = model.addVars(N, T, vtype=GRB.INTEGER, lb=0, name="c")
for i in range(N):
    for t in range(T):
        c[i, t].ub = int(nuclei_counts[i])
```

**Constraints:**
```python
for i in range(N):
    model.addConstr(
        gp.quicksum(c[i, t] for t in range(T)) == int(nuclei_counts[i]),
        name=f"nuclei_{i}"
    )
```

**Objective:**
```python
error_terms = []
for i in range(N):
    for m in range(M):
        pred_im = gp.quicksum(c[i, t] * profile[t, m] for t in range(T))
        residual = X[i, m] - alpha[m] - beta[m] * pred_im
        error_terms.append(residual * residual)

model.setObjective(gp.quicksum(error_terms), GRB.MINIMIZE)
```

**Solver settings:**
```python
model.setParam("TimeLimit", time_limit)
model.setParam("MIPGap", mip_gap)
model.setParam("OutputFlag", 1)
model.setParam("Threads", 0)  # Use all cores
```

### Warm-Start Strategy

```python
# From previous EM iteration (preferred)
if prev_c_values is not None:
    for i in range(N):
        for t in range(T):
            c[i, t].Start = int(prev_c_values[i, t])
# Or uniform distribution
else:
    for i in range(N):
        N_i = int(nuclei_counts[i])
        base = N_i // T
        remainder = N_i % T
        for t in range(T):
            c[i, t].Start = base + (1 if t < remainder else 0)
```

### Solution Extraction

```python
model.optimize()

if model.status in [GRB.OPTIMAL, GRB.TIME_LIMIT, GRB.MIP_GAP]:
    c_values = np.zeros((N, T), dtype=np.int64)
    for i in range(N):
        for t in range(T):
            c_values[i, t] = int(round(c[i, t].X))

    logging.info(f"Global IQP: status={model.status}, "
                 f"gap={model.MIPGap:.2%}, time={model.Runtime:.1f}s")
else:
    raise RuntimeError(f"Global IQP failed: status={model.status}")
```

### EM Integration

Modified `optimize_discrete_cell_assignment_em()` E-step:

```python
# Current (per-spot):
for i in range(N):
    c_values[i] = solve_discrete_cell_counts(spot_i_data, ...)

# New (global):
c_values = solve_discrete_cell_counts_global(all_spot_data, ...)
```

M-step unchanged (already global beta/alpha update via OLS).

### API Changes

**`optimize_discrete_cell_assignment_em()`:**
```python
def optimize_discrete_cell_assignment_em(
    ...
    # New params
    global_solve: bool = True,
    global_time_limit: float = 300.0,
    global_mip_gap: float = 0.05,
) -> Tuple[np.ndarray, np.ndarray, Dict[str, float], np.ndarray]:
```

**`CitegeistModel.run_discrete_cell_assignment()`:**
```python
def run_discrete_cell_assignment(
    self,
    ...
    global_solve: bool = True,
    global_time_limit: float = 300.0,
    global_mip_gap: float = 0.05,
) -> pd.DataFrame:
```

**Backward compatibility:** `global_solve=False` falls back to per-spot behavior.

## Testing

### Unit Tests

`tests/test_discrete_global_solve.py`:

1. `test_global_solve_basic()` - Valid cell counts summing to nuclei
2. `test_global_vs_perspot_consistency()` - High correlation on easy data
3. `test_global_solve_time_limit()` - Respects time limit, returns feasible solution

### Benchmark Validation

| Dataset | Baseline (per-spot) | Target (global) | Purpose |
|---------|---------------------|-----------------|---------|
| high_seg | ~0.76 prop corr | ≥0.76 | No regression |
| mixed | ~0.43 prop corr | ≥0.65 | Primary improvement |
| xenium | ~0.66 prop corr | ≥0.66 | Real data validation |

## File Changes

| File | Changes |
|------|---------|
| `CITEgeist/model/gurobi_impl.py` | Add `solve_discrete_cell_counts_global()` (~100 lines) |
| `CITEgeist/model/gurobi_impl.py` | Modify `optimize_discrete_cell_assignment_em()` E-step |
| `CITEgeist/model/citegeist_model.py` | Add `global_solve`, `global_time_limit`, `global_mip_gap` params |
| `tests/test_discrete_global_solve.py` | New test file |
| `Benchmarking/.../benchmark_discrete_simulation.py` | Add `--global-solve` CLI flag |
| `Benchmarking/.../slurm/sbatch_discrete_global_*.sh` | New SLURM scripts |

## Risks & Mitigations

| Risk | Mitigation |
|------|------------|
| 45K integer variables too slow | Time limit + MIP gap; fallback to per-spot |
| Memory usage | Gurobi handles sparse models; use 64GB node |
| Regression on high_seg/xenium | Run all 3 datasets before merging |

## Summary

Replace per-spot IQP with global IQP to enforce globally consistent marker-celltype behavior. Output is integer cell counts (nuclei labels) per spot, with proportions derived only for evaluation. Expected improvement: mixed prop corr from ~0.43 to ≥0.65.
