# Multi-Variant CITEgeist Xenium Comparison Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Extend Xenium benchmarking to compare all three CITEgeist variants (Continuous, Hybrid, Discrete) with proper documentation of preprocessing differences.

**Architecture:** Modify existing comparison scripts to add CITEgeist variants to the METHODS dict, add documentation constants, and refactor the GEX loader to handle all variants via a parameter.

**Tech Stack:** Python, pandas, numpy, scipy, JSON

---

## Task 1: Verify All Prediction Files Exist

**Files:**
- Verify: `Benchmarking/xenium_benchmarking/CITEgeist/output/manual/Xenium_region_{0-4}/`
- Verify: `Benchmarking/xenium_benchmarking/CITEgeist/output/hybrid_cellpose/Xenium_region_{0-4}/`
- Verify: `Benchmarking/xenium_benchmarking/CITEgeist/output_discrete_cellpose_fixed/Xenium_region_{0-4}/`

**Step 1: Run verification command**

```bash
for variant in "output/manual" "output/hybrid_cellpose" "output_discrete_cellpose_fixed"; do
  echo "=== $variant ==="
  for i in 0 1 2 3 4; do
    pred_file="Benchmarking/xenium_benchmarking/CITEgeist/${variant}/Xenium_region_${i}/Xenium_region_${i}_deconv_predictions.csv"
    if [ -f "$pred_file" ]; then
      echo "  Region $i: OK"
    else
      echo "  Region $i: MISSING"
    fi
  done
done
```

Expected: All 15 files (5 regions x 3 variants) exist.

**Step 2: Verify GEX layers exist**

```bash
for variant in "output/manual" "output/hybrid_cellpose" "output_discrete_cellpose_fixed"; do
  echo "=== $variant GEX ==="
  for i in 0 1 2 3 4; do
    layers_dir="Benchmarking/xenium_benchmarking/CITEgeist/${variant}/Xenium_region_${i}/Xenium_region_${i}_pass1/layers"
    if [ -d "$layers_dir" ]; then
      n_files=$(ls -1 "$layers_dir"/*.csv 2>/dev/null | wc -l)
      echo "  Region $i: $n_files layer files"
    else
      echo "  Region $i: NO LAYERS DIR"
    fi
  done
done
```

Expected: Each region has 7 layer files (one per cell type).

---

## Task 2: Update compare_all_methods.py - Add Constants

**Files:**
- Modify: `Benchmarking/xenium_benchmarking/evaluation/src/compare_all_methods.py:35-45`

**Step 1: Add METHOD_NOTES and CITEGEIST_COMMON_PARAMS after imports**

Add these constants after line 33 (after RESULTS_DIR definition):

```python
# =============================================================================
# METHOD DOCUMENTATION
# =============================================================================

METHOD_NOTES = {
    "CITEgeist_Continuous": "CLR antibody normalization, QP optimization",
    "CITEgeist_Hybrid": "CLR antibody normalization, QP optimization -> discretize via nuclei counts",
    "CITEgeist_Discrete": "Per-marker scaling (no CLR), IQP optimization with nuclei constraints",
}

CITEGEIST_COMMON_PARAMS = {
    "data_source": "data_protein_gt (protein-gated achievable-7)",
    "min_counts": 25,
    "gex_preprocessing": "filter_gex + preprocess_gex(target_sum=10000)",
    "cell_profile_dict": "ACHIEVABLE_7_CELL_PROFILE_DICT (alphaSMA-only fibroblasts)",
    "n_regions": 5,
}
```

**Step 2: Update METHODS dict**

Replace lines 35-45 with:

```python
METHODS = {
    # CITEgeist variants
    "CITEgeist_Continuous": BASE_DIR / "CITEgeist" / "output" / "manual",
    "CITEgeist_Hybrid": BASE_DIR / "CITEgeist" / "output" / "hybrid_cellpose",
    "CITEgeist_Discrete": BASE_DIR / "CITEgeist" / "output_discrete_cellpose_fixed",
    # Other methods
    "Cell2Location": BASE_DIR / "Cell2Location" / "output_protein_gt",
    "Tangram": BASE_DIR / "Tangram" / "output_protein_gt",
    "RCTD": BASE_DIR / "RCTD" / "output_protein_gt",
    "Seurat": BASE_DIR / "Seurat" / "output_protein_gt",
    # CARD: VIM excluded - VIM is a pan-mesenchymal marker expressed in RCC cells,
    # not fibroblast-specific in kidney tissue. Including it caused severe bias (93% fibroblasts).
    "CARD": BASE_DIR / "CARD" / "output_protein_gt_novim",
}
```

**Step 3: Commit**

```bash
git add Benchmarking/xenium_benchmarking/evaluation/src/compare_all_methods.py
git commit -m "feat(benchmark): add CITEgeist variant constants to compare_all_methods"
```

---

## Task 3: Update compare_all_methods.py - Add Documentation Output

**Files:**
- Modify: `Benchmarking/xenium_benchmarking/evaluation/src/compare_all_methods.py:86-92`

**Step 1: Add method notes printing after the header**

After line 91 (after printing the header), add:

```python
    # Print method notes
    print("\n--- Method Notes ---")
    for method, note in METHOD_NOTES.items():
        print(f"  {method}: {note}")

    print("\n--- Common CITEgeist Parameters ---")
    for key, val in CITEGEIST_COMMON_PARAMS.items():
        print(f"  {key}: {val}")
```

**Step 2: Widen column formatting for longer method names**

Update line 95 to use wider columns (20 chars instead of 16):

```python
    print(f"{'Method':<24} {'Pearson r':>10} {'RMSE':>10} {'MAE':>10} {'JSD':>10}")
    print("-" * 64)
```

Update line 105-114 to use 24-char method column:

```python
    for method in sorted_methods:
        s = all_results[method]["summary"]
        jsd_str = f"{s['overall_mean_jsd']:.4f}" if "overall_mean_jsd" in s else "N/A"
        print(
            f"{method:<24} "
            f"{s['overall_mean_pearson_r']:>10.4f} "
            f"{s['overall_mean_rmse']:>10.4f} "
            f"{s['overall_mean_mae']:>10.4f} "
            f"{jsd_str:>10}"
        )
```

**Step 3: Add method notes to JSON output**

Before line 209 (before writing JSON), add method_notes to the output:

```python
    output_data = {
        "method_notes": METHOD_NOTES,
        "citegeist_common_params": CITEGEIST_COMMON_PARAMS,
        "results": convert_numpy(all_results),
    }

    with open(combined_output, "w") as f:
        json.dump(output_data, f, indent=2)
```

**Step 4: Commit**

```bash
git add Benchmarking/xenium_benchmarking/evaluation/src/compare_all_methods.py
git commit -m "feat(benchmark): add documentation output to compare_all_methods"
```

---

## Task 4: Update compare_all_methods_gex.py - Add Unified Loader

**Files:**
- Modify: `Benchmarking/xenium_benchmarking/evaluation/src/compare_all_methods_gex.py:39-124`

**Step 1: Add CITEGEIST_VARIANTS dict after CELLTYPE_TO_FILENAME**

Add after line 49:

```python
# =============================================================================
# CITEGEIST VARIANT PATHS
# =============================================================================

CITEGEIST_VARIANTS = {
    "CITEgeist_Continuous": "output/manual",
    "CITEgeist_Hybrid": "output/hybrid_cellpose",
    "CITEgeist_Discrete": "output_discrete_cellpose_fixed",
}

METHOD_NOTES = {
    "CITEgeist_Continuous": "CLR antibody normalization, QP optimization",
    "CITEgeist_Hybrid": "CLR antibody normalization, QP optimization -> discretize via nuclei counts",
    "CITEgeist_Discrete": "Per-marker scaling (no CLR), IQP optimization with nuclei constraints",
}
```

**Step 2: Replace load_citegeist_gex with unified version**

Replace lines 63-92 with:

```python
def load_citegeist_gex(region_id: int, variant: str = "CITEgeist_Continuous") -> Dict[str, pd.DataFrame]:
    """Load CITEgeist GEX layers for any variant.

    Args:
        region_id: Xenium region ID (0-4)
        variant: One of CITEGEIST_VARIANTS keys

    Returns:
        Dict mapping cell type names to (spots x genes) DataFrames
    """
    variant_path = CITEGEIST_VARIANTS.get(variant)
    if variant_path is None:
        raise ValueError(f"Unknown CITEgeist variant: {variant}. Valid: {list(CITEGEIST_VARIANTS.keys())}")

    sample_name = f"Xenium_region_{region_id}"
    base_layers = BASE_DIR / "CITEgeist" / variant_path / sample_name / f"{sample_name}_pass1" / "layers"

    # Check pass1 subdirectory first (standard export), then direct layers dir
    layers_dir = base_layers / "pass1"
    if not layers_dir.exists():
        layers_dir = base_layers

    if not layers_dir.exists():
        logger.warning(f"{variant} layers not found for region {region_id}: {layers_dir}")
        return {}

    ct_dfs = {}
    for ct in ACHIEVABLE_7_CELL_TYPES:
        # CITEgeist replaces spaces with _ but keeps + in filenames
        ct_file = ct.replace(" ", "_")
        for pattern in [
            f"{ct_file}_layer_pass1.csv",
            f"{ct_file}_layer.csv",
        ]:
            layer_file = layers_dir / pattern
            if layer_file.exists():
                ct_dfs[ct] = pd.read_csv(layer_file, index_col=0)
                break

    return ct_dfs
```

**Step 3: Remove load_citegeist_discrete_gex function**

Delete lines 95-123 (the old `load_citegeist_discrete_gex` function).

**Step 4: Commit**

```bash
git add Benchmarking/xenium_benchmarking/evaluation/src/compare_all_methods_gex.py
git commit -m "feat(benchmark): add unified CITEgeist GEX loader with variant parameter"
```

---

## Task 5: Update compare_all_methods_gex.py - Update Main Loop

**Files:**
- Modify: `Benchmarking/xenium_benchmarking/evaluation/src/compare_all_methods_gex.py:232-244`

**Step 1: Update methods dict to use all CITEgeist variants**

Replace lines 238-244 with:

```python
    # Methods that produce GEX layers
    methods = {
        # CITEgeist variants
        "CITEgeist_Continuous": lambda rid: load_citegeist_gex(rid, "CITEgeist_Continuous"),
        "CITEgeist_Hybrid": lambda rid: load_citegeist_gex(rid, "CITEgeist_Hybrid"),
        "CITEgeist_Discrete": lambda rid: load_citegeist_gex(rid, "CITEgeist_Discrete"),
        # Other methods
        "Cell2Location": lambda rid: load_competitor_gex("Cell2Location", rid),
        "Tangram": lambda rid: load_competitor_gex("Tangram", rid),
        "scResolve": lambda rid: load_scresolve_gex(rid),
    }
```

**Step 2: Update summary print order**

Replace line 280 with:

```python
    for method_name in ["CITEgeist_Continuous", "CITEgeist_Hybrid", "CITEgeist_Discrete", "Cell2Location", "Tangram", "scResolve"]:
```

**Step 3: Update per-cell-type table methods**

Replace line 308 with:

```python
    sorted_methods = [m for m in ["CITEgeist_Continuous", "CITEgeist_Hybrid", "CITEgeist_Discrete", "Cell2Location", "Tangram", "scResolve"] if all_results.get(m)]
```

**Step 4: Add method notes to JSON output**

Before line 366 (before writing JSON), wrap results with metadata:

```python
    output_data = {
        "method_notes": METHOD_NOTES,
        "results": convert_numpy(all_results),
    }

    with open(results_dir / "full_comparison_gex.json", "w") as f:
        json.dump(output_data, f, indent=2)
```

**Step 5: Commit**

```bash
git add Benchmarking/xenium_benchmarking/evaluation/src/compare_all_methods_gex.py
git commit -m "feat(benchmark): update GEX comparison to include all CITEgeist variants"
```

---

## Task 6: Create SLURM Submission Scripts

**Files:**
- Create: `Benchmarking/xenium_benchmarking/evaluation/src/sbatch_compare_multivariant.sh`

**Step 1: Create SLURM script**

```bash
#!/bin/bash
#SBATCH --job-name=compare_multivariant
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=logs/compare_multivariant_%j.out
#SBATCH --error=logs/compare_multivariant_%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

module load gurobi/12.0.3
eval "$(conda shell.bash hook)"
conda activate ~/alc376_bgfs/envs/CITEgeist_env

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist

# Run proportion comparison
echo "=== Running Proportion Comparison ==="
python Benchmarking/xenium_benchmarking/evaluation/src/compare_all_methods.py

# Run GEX comparison
echo "=== Running GEX Comparison ==="
python Benchmarking/xenium_benchmarking/evaluation/src/compare_all_methods_gex.py

echo "=== Complete ==="
```

**Step 2: Ensure logs directory exists**

```bash
mkdir -p Benchmarking/xenium_benchmarking/evaluation/src/logs
```

**Step 3: Commit**

```bash
git add Benchmarking/xenium_benchmarking/evaluation/src/sbatch_compare_multivariant.sh
git commit -m "feat(benchmark): add SLURM script for multi-variant comparison"
```

---

## Task 7: Run Comparison and Verify Results

**Step 1: Submit SLURM job**

```bash
cd Benchmarking/xenium_benchmarking/evaluation/src
sbatch sbatch_compare_multivariant.sh
```

**Step 2: Monitor job**

```bash
squeue -u alc376 --cluster=htc
```

**Step 3: Verify output files exist**

```bash
ls -la Benchmarking/xenium_benchmarking/evaluation/results/method_comparison/
```

Expected files:
- `full_comparison.json` (proportion results with all 8 methods)
- `full_comparison_gex.json` (GEX results with 6 methods)
- Individual method result files

**Step 4: Verify JSON structure includes method notes**

```bash
head -30 Benchmarking/xenium_benchmarking/evaluation/results/method_comparison/full_comparison.json
```

Expected: JSON should start with `method_notes` and `citegeist_common_params` keys.

---

## Task 8: Final Commit

**Step 1: Review changes**

```bash
git status
git diff --stat HEAD~5
```

**Step 2: Tag completion**

```bash
git tag -a v0.1.0-multivariant-comparison -m "Multi-variant CITEgeist Xenium comparison complete"
```

---

## Summary

| Task | Description | Files |
|------|-------------|-------|
| 1 | Verify prediction files | Shell verification |
| 2 | Add constants | `compare_all_methods.py` |
| 3 | Add documentation output | `compare_all_methods.py` |
| 4 | Add unified GEX loader | `compare_all_methods_gex.py` |
| 5 | Update GEX main loop | `compare_all_methods_gex.py` |
| 6 | Create SLURM script | `sbatch_compare_multivariant.sh` |
| 7 | Run and verify | SLURM job |
| 8 | Final commit | Git tag |
