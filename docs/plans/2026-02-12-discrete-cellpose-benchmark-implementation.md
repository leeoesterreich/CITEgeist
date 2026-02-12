# Discrete Cellpose Benchmark Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Merge discrete cell assignment feature with Cellpose fixes, then benchmark on Xenium pseudo-Visium with full Cellpose → IQP → GEX pipeline.

**Architecture:** Two-phase approach: (1) Git merge with conflict resolution in segmentation.py, (2) New benchmark script that runs Cellpose on morphology images, feeds nuclei counts to discrete IQP solver, runs GEX deconvolution, and outputs for evaluation.

**Tech Stack:** Python 3.10, Gurobi (IQP), Cellpose (nuclei segmentation), pandas/numpy, SLURM (GPU + CPU partitions)

---

## Task 1: Backup Existing Results

**Files:**
- Copy: `Benchmarking/xenium_benchmarking/evaluation/results/method_comparison/full_comparison_gex.json`

**Step 1: Create backup of continuous baseline results**

```bash
cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist
cp Benchmarking/xenium_benchmarking/evaluation/results/method_comparison/full_comparison_gex.json \
   Benchmarking/xenium_benchmarking/evaluation/results/method_comparison/full_comparison_gex_continuous_baseline.json
cp Benchmarking/xenium_benchmarking/evaluation/results/method_comparison/full_comparison.json \
   Benchmarking/xenium_benchmarking/evaluation/results/method_comparison/full_comparison_continuous_baseline.json
```

**Step 2: Git track comparison files**

```bash
git add Benchmarking/xenium_benchmarking/evaluation/results/method_comparison/*.json
git commit -m "results(benchmark): backup continuous baseline and track comparison JSONs"
```

---

## Task 2: Merge Dev into Feature Branch

**Files:**
- Modify: `.worktrees/feat-discrete-cell-assignment/CITEgeist/model/segmentation.py`

**Step 1: Navigate to feature worktree and merge dev**

```bash
cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/.worktrees/feat-discrete-cell-assignment
git fetch origin
git merge dev
```

Expected: Conflict in `CITEgeist/model/segmentation.py`

**Step 2: Resolve segmentation.py conflict**

The conflict is between:
- **Dev branch additions** (lines ~34-145): `PLATFORM_SPOT_DIAMETERS_UM`, `PLATFORM_SPOT_SPACINGS_UM` constants, and `detect_spot_diameter_pixels()` function
- **Feature branch**: No modifications to this area (changes are in `gurobi_impl.py` and `citegeist_model.py`)

Resolution: Accept BOTH changes (they're additive, non-overlapping). The dev version has all the content from the feature branch plus the new auto-detection features.

```bash
# After reviewing the conflict, accept the merged version
git checkout --theirs CITEgeist/model/segmentation.py
git add CITEgeist/model/segmentation.py
```

**Step 3: Complete the merge**

```bash
git commit -m "Merge branch 'dev' into feat/discrete-cell-assignment

Incorporate Cellpose auto-detection features:
- PLATFORM_SPOT_DIAMETERS_UM constants
- detect_spot_diameter_pixels() function
- CARD benchmarking updates"
```

---

## Task 3: Verify Merge with Unit Tests

**Files:**
- Test: `tests/test_discrete_assignment.py`

**Step 1: Load environment and run discrete assignment tests**

```bash
cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/.worktrees/feat-discrete-cell-assignment
module load gurobi/12.0.3
eval "$(conda shell.bash hook)" && conda activate ~/alc376_bgfs/envs/CITEgeist_env
pytest tests/test_discrete_assignment.py -v
```

Expected: All tests pass

**Step 2: Run segmentation tests to verify Cellpose features intact**

```bash
pytest tests/test_segmentation.py -v
```

Expected: All tests pass

---

## Task 4: Merge Feature Branch Back to Dev

**Files:**
- None (git operation)

**Step 1: Return to main repo and merge**

```bash
cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist
git merge feat/discrete-cell-assignment -m "$(cat <<'EOF'
Merge branch 'feat/discrete-cell-assignment' into dev

Add discrete cell assignment for Module 3:
- IQP solver for integer cell counts (solve_discrete_cell_counts)
- EM algorithm for beta optimization (optimize_discrete_cell_assignment_em)
- CitegeistModel.run_discrete_cell_assignment() entry point
- Discrete mode for run_cell_expression_pass1()

Benchmarked on Xenium pseudo-Visium:
- +10% Pearson correlation on proportions
- +10-16% on GEX deconvolution

Co-Authored-By: Claude Opus 4.5 <noreply@anthropic.com>
EOF
"
```

**Step 2: Verify merge completed cleanly**

```bash
git log --oneline -5
git status
```

---

## Task 5: Add NRMSE to GEX Evaluation

**Files:**
- Modify: `Benchmarking/xenium_benchmarking/evaluation/src/compare_all_methods_gex.py:184-193`

**Step 1: Write test for NRMSE metric**

Create temporary test file to verify NRMSE calculation:

```python
# Quick verification in Python REPL
import numpy as np
gt = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
pred = np.array([1.1, 2.1, 2.9, 4.2, 4.8])
rmse = np.sqrt(np.mean((gt - pred) ** 2))
nrmse = rmse / (gt.max() - gt.min())
print(f"RMSE: {rmse:.4f}, NRMSE: {nrmse:.4f}")
# Expected: RMSE ~0.15, NRMSE ~0.038
```

**Step 2: Add NRMSE to compute_gex_metrics function**

In `compare_all_methods_gex.py`, modify the `compute_gex_metrics` function (around line 184):

```python
    rmse = np.sqrt(np.mean((gt_flat - pred_flat) ** 2))
    mae = np.mean(np.abs(gt_flat - pred_flat))

    # NRMSE: normalize by range of ground truth
    gt_range = gt_flat.max() - gt_flat.min()
    nrmse = rmse / gt_range if gt_range > 0 else np.nan

    return {
        "pearson_r": float(pearson_r) if not np.isnan(pearson_r) else np.nan,
        "rmse": float(rmse),
        "nrmse": float(nrmse) if not np.isnan(nrmse) else np.nan,
        "mae": float(mae),
        "n_genes": len(common_genes_lower),
        "n_spots": len(common_spots),
    }
```

**Step 3: Update summary printing to include NRMSE**

Find the summary print section (around line 238) and add NRMSE column:

```python
    print(f"{'Method':<16} {'Pearson r':>10} {'RMSE':>10} {'NRMSE':>10} {'MAE':>10} {'Cell Types':>12} {'Regions':>8}")
    print("-" * 76)
```

And update the per-method printing (around line 250-257):

```python
        mean_nrmse = valid["nrmse"].mean()
        print(f"{method_name:<16} {mean_r:>10.4f} {mean_rmse:>10.4f} {mean_nrmse:>10.4f} {mean_mae:>10.4f} {n_ct:>12} {n_reg:>8}")

        method_summaries[method_name] = {
            "mean_pearson_r": mean_r,
            "mean_rmse": mean_rmse,
            "mean_nrmse": mean_nrmse,
            "mean_mae": mean_mae,
            ...
        }
```

**Step 4: Commit the NRMSE addition**

```bash
git add Benchmarking/xenium_benchmarking/evaluation/src/compare_all_methods_gex.py
git commit -m "feat(eval): add NRMSE metric to GEX comparison

NRMSE normalizes RMSE by ground truth range for scale-independent
comparison across cell types and methods."
```

---

## Task 6: Create Discrete Cellpose Benchmark Script

**Files:**
- Create: `Benchmarking/xenium_benchmarking/CITEgeist/src/benchmark_discrete_cellpose.py`

**Step 1: Create the benchmark script**

```python
#!/usr/bin/env python
"""
Benchmark discrete cell assignment with Cellpose nuclei detection.

Full pipeline: Morphology image → Cellpose → nuclei counts → IQP → GEX deconvolution

Usage:
    python benchmark_discrete_cellpose.py --region 0 --output-dir output_discrete_cellpose

    # All regions
    for r in 0 1 2 3 4; do
        python benchmark_discrete_cellpose.py --region $r --output-dir output_discrete_cellpose
    done
"""
from __future__ import annotations

import argparse
import json
import logging
import sys
import time
from pathlib import Path
from typing import Dict, Optional

import cv2
import numpy as np
import pandas as pd
import scanpy as sc

REPO_ROOT = Path(__file__).parent.parent.parent.parent.parent
sys.path.insert(0, str(REPO_ROOT))

from CITEgeist.model.citegeist_model import CitegeistModel
from CITEgeist.model.segmentation import (
    run_cellpose_nuclei_segmentation,
    assign_nuclei_centroids_to_spots,
    _prepare_rgb_uint8,
)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)

# Paths
DATA_DIR = REPO_ROOT / "Benchmarking" / "xenium_pseudovisium" / "data_protein_gt"
IMAGE_DIR = REPO_ROOT / "Benchmarking" / "xenium_benchmarking" / "scResolve" / "images" / "morphology_hires"
COORD_DIR = IMAGE_DIR  # coord_info.json is in same directory as hires images

# Achievable-7 profile for Xenium pseudo-Visium benchmark
ACHIEVABLE_7_PROFILE = {
    "B cells": {"Major": ["CD19"]},
    "CD4+ T cells": {"Major": ["CD4"]},
    "CD8+ T cells": {"Major": ["CD8A"]},
    "Macrophages": {"Major": ["CD68"]},
    "Endothelial": {"Major": ["PECAM1"]},
    "Epithelial": {"Major": ["EPCAM"]},
    "Fibroblasts": {"Major": ["ACTA2"]},
}


def load_morphology_image(region_id: int) -> np.ndarray:
    """Load hires morphology image for region."""
    image_path = IMAGE_DIR / f"Xenium_region_{region_id}" / "morphology.png"
    if not image_path.exists():
        raise FileNotFoundError(f"Morphology image not found: {image_path}")
    bgr = cv2.imread(str(image_path), cv2.IMREAD_COLOR)
    if bgr is None:
        raise ValueError(f"Failed to load image: {image_path}")
    rgb = cv2.cvtColor(bgr, cv2.COLOR_BGR2RGB)
    return _prepare_rgb_uint8(rgb)


def load_coord_info(region_id: int) -> Dict:
    """Load coordinate info for region."""
    coord_path = COORD_DIR / f"Xenium_region_{region_id}" / "coord_info.json"
    if not coord_path.exists():
        raise FileNotFoundError(f"Coord info not found: {coord_path}")
    with open(coord_path) as f:
        return json.load(f)


def run_cellpose_on_region(
    region_id: int,
    use_gpu: bool = False,
    diameter: Optional[float] = None,
) -> pd.Series:
    """
    Run Cellpose on region morphology and return nuclei counts per spot.

    Args:
        region_id: Xenium region ID (0-4)
        use_gpu: Use GPU for Cellpose
        diameter: Cellpose diameter param (auto-detect if None)

    Returns:
        Series with spot names as index, nuclei counts as values
    """
    logger.info(f"Loading morphology image for region {region_id}")
    image = load_morphology_image(region_id)
    coord = load_coord_info(region_id)

    # Get spot info
    adata_gex = sc.read_h5ad(DATA_DIR / "h5ad_objects" / f"Xenium_region_{region_id}_GEX.h5ad")
    spot_names = adata_gex.obs_names

    # Spot coordinates are in microns; convert to pixels using pixel_size
    pixel_size_um = float(coord["pixel_size"])  # microns per pixel
    spot_diameter_um = 55.0  # Visium standard
    spot_radius_px = (spot_diameter_um / pixel_size_um) / 2.0

    # Get spot centers in pixel coordinates
    # Coordinates in h5ad are in microns relative to crop origin
    spatial_coords = adata_gex.obsm["spatial"]  # microns
    spot_centers_px = spatial_coords / pixel_size_um  # convert to pixels

    logger.info(f"Running Cellpose segmentation (GPU={use_gpu})")
    t0 = time.time()
    masks, centroids = run_cellpose_nuclei_segmentation(
        image_rgb_uint8=image,
        use_gpu=use_gpu,
        diameter=diameter,
    )
    logger.info(f"Cellpose found {len(centroids)} nuclei in {time.time()-t0:.1f}s")

    # Assign nuclei to spots
    nuclei_counts = assign_nuclei_centroids_to_spots(
        centroids_xy=centroids,
        spot_centers_xy=spot_centers_px,
        spot_radius_px=spot_radius_px,
        spot_names=spot_names,
    )

    logger.info(f"Nuclei per spot: mean={nuclei_counts.mean():.1f}, max={nuclei_counts.max()}")
    return nuclei_counts


def run_discrete_benchmark(
    region_id: int,
    output_dir: Path,
    use_gpu: bool = False,
    cellpose_diameter: Optional[float] = None,
    max_em_iterations: int = 20,
) -> Dict:
    """
    Run full discrete benchmark pipeline for one region.

    Returns dict with timing and basic stats.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    sample_name = f"Xenium_region_{region_id}"
    logger.info(f"=== Starting discrete benchmark for {sample_name} ===")

    # Step 1: Cellpose segmentation
    t0 = time.time()
    nuclei_counts = run_cellpose_on_region(
        region_id=region_id,
        use_gpu=use_gpu,
        diameter=cellpose_diameter,
    )
    cellpose_time = time.time() - t0

    # Save nuclei counts
    nuclei_path = output_dir / f"{sample_name}_nuclei_counts.csv"
    nuclei_counts.to_csv(nuclei_path)

    # Step 2: Load data and initialize model
    adata_gex = sc.read_h5ad(DATA_DIR / "h5ad_objects" / f"{sample_name}_GEX.h5ad")
    adata_cite = sc.read_h5ad(DATA_DIR / "h5ad_objects" / f"{sample_name}_CITE.h5ad")

    model = CitegeistModel(
        sample_name=sample_name,
        output_folder=str(output_dir),
        simulation=True,
        gene_expression_adata=adata_gex,
        antibody_capture_adata=adata_cite,
    )

    # Step 3: Preprocess
    model.preprocess_gex(min_counts=25)
    model.preprocess_antibody_discrete()  # CLR without per-spot normalization
    model.load_cell_profile_dict(ACHIEVABLE_7_PROFILE)

    # Step 4: Discrete cell assignment
    t1 = time.time()
    cell_counts_df = model.run_discrete_cell_assignment(
        nuclei_counts=nuclei_counts,
        max_em_iterations=max_em_iterations,
    )
    discrete_time = time.time() - t1

    # Step 5: GEX deconvolution
    t2 = time.time()
    model.run_cell_expression_pass1(
        use_discrete_mode=True,
        cell_counts=cell_counts_df,
    )
    gex_time = time.time() - t2

    # Save summary
    results = {
        "region_id": region_id,
        "n_spots": len(nuclei_counts),
        "total_nuclei": int(nuclei_counts.sum()),
        "mean_nuclei_per_spot": float(nuclei_counts.mean()),
        "cellpose_time_s": cellpose_time,
        "discrete_time_s": discrete_time,
        "gex_time_s": gex_time,
        "total_time_s": cellpose_time + discrete_time + gex_time,
    }

    with open(output_dir / f"{sample_name}_benchmark_summary.json", "w") as f:
        json.dump(results, f, indent=2)

    logger.info(f"=== Completed {sample_name} in {results['total_time_s']:.1f}s ===")
    return results


def main():
    parser = argparse.ArgumentParser(description="Discrete Cellpose benchmark")
    parser.add_argument("--region", type=int, required=True, help="Region ID (0-4)")
    parser.add_argument("--output-dir", type=Path, default=Path("output_discrete_cellpose"))
    parser.add_argument("--use-gpu", action="store_true", help="Use GPU for Cellpose")
    parser.add_argument("--cellpose-diameter", type=float, default=None)
    parser.add_argument("--max-em-iterations", type=int, default=20)
    args = parser.parse_args()

    run_discrete_benchmark(
        region_id=args.region,
        output_dir=args.output_dir,
        use_gpu=args.use_gpu,
        cellpose_diameter=args.cellpose_diameter,
        max_em_iterations=args.max_em_iterations,
    )


if __name__ == "__main__":
    main()
```

**Step 2: Commit the benchmark script**

```bash
git add Benchmarking/xenium_benchmarking/CITEgeist/src/benchmark_discrete_cellpose.py
git commit -m "feat(benchmark): add discrete cellpose benchmark script

Full pipeline: Cellpose → nuclei counts → IQP discrete assignment → GEX
Outputs compatible with existing evaluation framework."
```

---

## Task 7: Create SLURM Submission Scripts

**Files:**
- Create: `Benchmarking/xenium_benchmarking/CITEgeist/slurm/sbatch_discrete_cellpose_gpu.sh`
- Create: `Benchmarking/xenium_benchmarking/CITEgeist/slurm/sbatch_discrete_cellpose_cpu.sh`

**Step 1: Create GPU script for Cellpose**

```bash
#!/bin/bash
#SBATCH --job-name=cellpose_discrete
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --array=0-4
#SBATCH --output=logs/cellpose_discrete_%A_%a.out
#SBATCH --error=logs/cellpose_discrete_%A_%a.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

module load gurobi/12.0.3
eval "$(conda shell.bash hook)"
conda activate ~/alc376_bgfs/envs/CITEgeist_env

REGION_ID=$SLURM_ARRAY_TASK_ID
OUTPUT_DIR=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/CITEgeist/output_discrete_cellpose

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist

python Benchmarking/xenium_benchmarking/CITEgeist/src/benchmark_discrete_cellpose.py \
    --region $REGION_ID \
    --output-dir $OUTPUT_DIR \
    --use-gpu \
    --max-em-iterations 20

echo "Completed region $REGION_ID"
```

**Step 2: Create CPU-only fallback script**

```bash
#!/bin/bash
#SBATCH --job-name=cellpose_discrete_cpu
#SBATCH --partition=htc
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --array=0-4
#SBATCH --output=logs/cellpose_discrete_cpu_%A_%a.out
#SBATCH --error=logs/cellpose_discrete_cpu_%A_%a.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

module load gurobi/12.0.3
eval "$(conda shell.bash hook)"
conda activate ~/alc376_bgfs/envs/CITEgeist_env

REGION_ID=$SLURM_ARRAY_TASK_ID
OUTPUT_DIR=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/CITEgeist/output_discrete_cellpose

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist

python Benchmarking/xenium_benchmarking/CITEgeist/src/benchmark_discrete_cellpose.py \
    --region $REGION_ID \
    --output-dir $OUTPUT_DIR \
    --max-em-iterations 20

echo "Completed region $REGION_ID"
```

**Step 3: Create logs directory and commit**

```bash
mkdir -p Benchmarking/xenium_benchmarking/CITEgeist/slurm/logs
git add Benchmarking/xenium_benchmarking/CITEgeist/slurm/sbatch_discrete_cellpose_gpu.sh
git add Benchmarking/xenium_benchmarking/CITEgeist/slurm/sbatch_discrete_cellpose_cpu.sh
git commit -m "feat(slurm): add discrete cellpose benchmark submission scripts

GPU version for fast Cellpose, CPU fallback for queue availability."
```

---

## Task 8: Run Benchmarks

**Step 1: Submit GPU jobs**

```bash
cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/CITEgeist/slurm
sbatch sbatch_discrete_cellpose_gpu.sh
```

**Step 2: Monitor progress**

```bash
squeue -u alc376
# Check logs
tail -f logs/cellpose_discrete_*.out
```

**Step 3: Verify outputs exist for all regions**

```bash
ls -la ../output_discrete_cellpose/
# Should see Xenium_region_0 through Xenium_region_4 directories
```

---

## Task 9: Run Evaluation

**Step 1: Update compare_all_methods_gex.py to include discrete results**

Add loader function for discrete outputs (they use same format as continuous):

```python
def load_citegeist_discrete_gex(region_id: int) -> Dict[str, pd.DataFrame]:
    """Load CITEgeist discrete GEX layers."""
    sample_name = f"Xenium_region_{region_id}"
    layers_dir = BASE_DIR / "CITEgeist" / "output_discrete_cellpose" / f"{sample_name}_pass1" / "layers" / "pass1"

    if not layers_dir.exists():
        return {}

    ct_dfs = {}
    for ct in ACHIEVABLE_7_CELL_TYPES:
        ct_file = ct.replace(" ", "_")
        for pattern in [f"{ct_file}_layer_pass1.csv", f"{ct_file}_layer.csv"]:
            layer_file = layers_dir / pattern
            if layer_file.exists():
                ct_dfs[ct] = pd.read_csv(layer_file, index_col=0)
                break
    return ct_dfs
```

**Step 2: Run evaluation**

```bash
cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/evaluation/src
python compare_all_methods_gex.py
```

**Step 3: Commit updated results**

```bash
cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist
git add Benchmarking/xenium_benchmarking/evaluation/results/method_comparison/*.json
git commit -m "results(benchmark): add CITEgeist discrete cellpose results

Metrics include Pearson, RMSE, NRMSE, MAE for achievable-7 cell types."
```

---

## Task 10: Final Verification and Summary

**Step 1: Compare discrete vs continuous**

```python
import json
with open("full_comparison_gex.json") as f:
    discrete = json.load(f)
with open("full_comparison_gex_continuous_baseline.json") as f:
    continuous = json.load(f)

print("CITEgeist Continuous:", continuous.get("CITEgeist", {}).get("mean_pearson_r"))
print("CITEgeist Discrete:", discrete.get("CITEgeist-Discrete", {}).get("mean_pearson_r"))
```

**Step 2: Verify all commits are in place**

```bash
git log --oneline -10
```

Expected commits:
1. docs: add discrete cellpose benchmark design
2. results(benchmark): backup continuous baseline and track comparison JSONs
3. Merge branch 'feat/discrete-cell-assignment' into dev
4. feat(eval): add NRMSE metric to GEX comparison
5. feat(benchmark): add discrete cellpose benchmark script
6. feat(slurm): add discrete cellpose benchmark submission scripts
7. results(benchmark): add CITEgeist discrete cellpose results

---

## Summary

| Task | Description | Files |
|------|-------------|-------|
| 1 | Backup existing results | JSON backups |
| 2 | Merge dev into feature branch | segmentation.py conflict |
| 3 | Verify with unit tests | test_discrete_assignment.py |
| 4 | Merge feature back to dev | Git merge |
| 5 | Add NRMSE metric | compare_all_methods_gex.py |
| 6 | Create benchmark script | benchmark_discrete_cellpose.py |
| 7 | Create SLURM scripts | sbatch scripts |
| 8 | Run benchmarks | SLURM submission |
| 9 | Run evaluation | compare_all_methods_gex.py |
| 10 | Final verification | Git log, JSON comparison |
