# CLR Preprocessing Test for Discrete Model - Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Test whether using CLR preprocessing (from continuous model) with discrete IQP improves proportion correlation from ~0.43 to ≥0.65 on mixed dataset.

**Architecture:** Add `--use-clr-preprocessing` flag to the discrete benchmark script that calls `model.preprocess_antibody()` (CLR) instead of `model.preprocess_antibody_discrete()` (per-marker scaling). Run 5-replicate benchmark on mixed dataset and compare results.

**Tech Stack:** Python 3.10, Gurobi 12.0.3, SLURM, argparse

---

## Task 1: Add CLI Flag and Function Parameter

**Files:**
- Modify: `Benchmarking/simulation_benchmarking/CITEgeist/src/benchmark_discrete_simulation.py:266-273` (function signature)
- Modify: `Benchmarking/simulation_benchmarking/CITEgeist/src/benchmark_discrete_simulation.py:279` (results dict)
- Modify: `Benchmarking/simulation_benchmarking/CITEgeist/src/benchmark_discrete_simulation.py:586-590` (argparse)
- Modify: `Benchmarking/simulation_benchmarking/CITEgeist/src/benchmark_discrete_simulation.py:592-610` (function call)

### Step 1.1: Add parameter to function signature

Edit line 272-273 to add `use_clr_preprocessing` parameter:

```python
    lambda_continuous: float = 50.0,
    use_clr_preprocessing: bool = False,
) -> Dict[str, Any]:
```

### Step 1.2: Add to results dict

Edit line 279 to include the new parameter in logging:

```python
    results = {"replicate_id": replicate_id, "condition": condition, "mode": mode, "lambda_sparse": lambda_sparse, "scale_mode": scale_mode, "lambda_prior": lambda_prior, "global_solve": global_solve, "global_time_limit": global_time_limit, "global_mip_gap": global_mip_gap, "use_continuous_prior": use_continuous_prior, "lambda_continuous": lambda_continuous, "use_clr_preprocessing": use_clr_preprocessing}
```

### Step 1.3: Add argparse argument

After line 589, add:

```python
    parser.add_argument("--use-clr-preprocessing", action="store_true", default=False,
                        help="Use continuous CLR preprocessing instead of discrete per-marker scaling")
```

### Step 1.4: Pass to function call

Edit line 609-610 to add the new argument:

```python
        lambda_continuous=args.lambda_continuous,
        use_clr_preprocessing=args.use_clr_preprocessing,
    )
```

### Step 1.5: Test CLI parsing

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python Benchmarking/simulation_benchmarking/CITEgeist/src/benchmark_discrete_simulation.py --help | grep -A1 clr`

Expected: See `--use-clr-preprocessing` in help output

### Step 1.6: Commit

```bash
git add Benchmarking/simulation_benchmarking/CITEgeist/src/benchmark_discrete_simulation.py
git commit -m "feat(benchmark): add --use-clr-preprocessing parameter to discrete benchmark

Adds CLI flag and function parameter. Does not change behavior yet.

Co-Authored-By: Claude Opus 4.5 <noreply@anthropic.com>"
```

---

## Task 2: Add Preprocessing Branch

**Files:**
- Modify: `Benchmarking/simulation_benchmarking/CITEgeist/src/benchmark_discrete_simulation.py:369` (preprocessing call)

### Step 2.1: Replace preprocessing line

Replace line 369:

```python
    model.preprocess_antibody_discrete(scale_mode=scale_mode)
```

With:

```python
    if use_clr_preprocessing:
        logger.info("Using CLR preprocessing (continuous model style)")
        model.preprocess_antibody()  # CLR normalization
    else:
        model.preprocess_antibody_discrete(scale_mode=scale_mode)
```

### Step 2.2: Verify syntax

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m py_compile Benchmarking/simulation_benchmarking/CITEgeist/src/benchmark_discrete_simulation.py && echo "Syntax OK"`

Expected: `Syntax OK`

### Step 2.3: Commit

```bash
git add Benchmarking/simulation_benchmarking/CITEgeist/src/benchmark_discrete_simulation.py
git commit -m "feat(benchmark): implement CLR preprocessing branch for discrete benchmark

When --use-clr-preprocessing is passed, uses model.preprocess_antibody()
instead of model.preprocess_antibody_discrete().

Co-Authored-By: Claude Opus 4.5 <noreply@anthropic.com>"
```

---

## Task 3: Create SLURM Script

**Files:**
- Create: `Benchmarking/simulation_benchmarking/CITEgeist/slurm/sbatch_discrete_clr_test.sh`

### Step 3.1: Write SLURM script

```bash
#!/bin/bash
#SBATCH --job-name=discrete_clr_test
#SBATCH --output=logs/discrete_clr_test_%A_%a.out
#SBATCH --error=logs/discrete_clr_test_%A_%a.err
#SBATCH --array=0-4
#SBATCH --time=02:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --partition=htc
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

module load gurobi/12.0.3

eval "$(conda shell.bash hook)"
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist

python Benchmarking/simulation_benchmarking/CITEgeist/src/benchmark_discrete_simulation.py \
    --replicate-id ${SLURM_ARRAY_TASK_ID} \
    --condition mixed \
    --mode dapi \
    --output-dir Benchmarking/simulation_benchmarking/CITEgeist/output_discrete_clr/mixed/dapi \
    --use-clr-preprocessing \
    --global-solve \
    --global-time-limit 600 \
    --global-mip-gap 0.05 \
    --max-em-iterations 10
```

### Step 3.2: Create output directories

Run: `mkdir -p Benchmarking/simulation_benchmarking/CITEgeist/output_discrete_clr/mixed/dapi`

### Step 3.3: Commit

```bash
git add Benchmarking/simulation_benchmarking/CITEgeist/slurm/sbatch_discrete_clr_test.sh
git commit -m "chore(benchmark): add SLURM script for CLR preprocessing test

Tests discrete IQP with CLR preprocessing on mixed dataset (5 replicates).

Co-Authored-By: Claude Opus 4.5 <noreply@anthropic.com>"
```

---

## Task 4: Submit Benchmark

### Step 4.1: Submit job

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/simulation_benchmarking/CITEgeist/slurm && sbatch sbatch_discrete_clr_test.sh`

Expected: `Submitted batch job XXXXXXX`

### Step 4.2: Monitor job

Run: `squeue -u $USER --format="%.10i %.9P %.30j %.8u %.2t %.10M %.6D %R"`

Expected: See `discrete_clr_test` jobs in queue or running

---

## Task 5: Validate Results (After Jobs Complete)

### Step 5.1: Check job completion

Run: `ls -la Benchmarking/simulation_benchmarking/CITEgeist/output_discrete_clr/mixed/dapi/Wu_rep_*/benchmark_results.json | wc -l`

Expected: `5` (all replicates complete)

### Step 5.2: Extract CLR results

Run:
```bash
echo "=== CLR Preprocessing ===" && for f in Benchmarking/simulation_benchmarking/CITEgeist/output_discrete_clr/mixed/dapi/Wu_rep_*/benchmark_results.json; do python3 -c "import json; d=json.load(open('$f')); print(f'{f}: prop_corr={d[\"proportion_correlation\"]:.3f}')"; done
```

### Step 5.3: Compare to baseline

Run:
```bash
echo "=== Per-marker baseline ===" && for f in Benchmarking/simulation_benchmarking/CITEgeist/output_discrete_global/mixed/dapi/Wu_rep_*/benchmark_results.json; do python3 -c "import json; d=json.load(open('$f')); print(f'{f}: prop_corr={d[\"proportion_correlation\"]:.3f}')"; done
```

### Step 5.4: Update MEMORY.md

Based on results, update `/ihome/alee/alc376/.claude/projects/-ix1-alee-LO-LAB-Personal-Alexander-Chang-alc376-CITEgeist/memory/MEMORY.md` with findings under "## Discrete Model Investigation".

---

## Success Criteria

| Outcome | CLR prop_corr | Interpretation | Next Step |
|---------|--------------|----------------|-----------|
| **Success** | ≥0.65 | Preprocessing is root cause | Add `scale_mode='clr'` to discrete API |
| **Partial** | 0.50-0.64 | CLR helps but not enough | Investigate IQP objective formulation |
| **Failure** | ~0.43 | Preprocessing NOT the cause | Investigate beta/alpha fitting |

---

## Summary

| Task | Description | Estimated Time |
|------|-------------|----------------|
| 1 | Add CLI flag and parameter | 5 min |
| 2 | Add preprocessing branch | 3 min |
| 3 | Create SLURM script | 3 min |
| 4 | Submit benchmark | 2 min |
| 5 | Validate results | 5 min (after ~30-60 min job runtime) |

**Total implementation time:** ~15 min
**Total wall time (including job):** ~45-75 min
