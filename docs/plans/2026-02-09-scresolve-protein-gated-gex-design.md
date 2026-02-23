# scResolve Protein-Gated GEX Benchmarking Design

**Date:** 2026-02-09
**Status:** Ready for implementation
**Goal:** Add scResolve to GEX deconvolution benchmark using protein-gated cell type assignment

---

## Background

### Current GEX Benchmarking State

| Method | Mean Pearson r | Coverage | Notes |
|--------|---------------|----------|-------|
| CITEgeist | 0.41 | 100% | Pass 2 deconvolved layers |
| Cell2Location | 0.40 | 100% | Posterior mean GEX |
| Tangram | NaN | - | Output format issue |
| **scResolve** | **Missing** | - | Not in comparison |

### Why scResolve is Missing

scResolve outputs per-cell molecule assignments, not per-cell-type aggregates. Previous approaches used Hungarian matching to transfer GT cell types to scResolve cells, which gives scResolve "perfect" cell type information and doesn't test its actual capability.

### Proposed Solution

Apply the **same protein-based hierarchical gating** used for ground truth (`create_protein_gt.py`) to scResolve's protein outputs. This tests whether scResolve's protein recovery is good enough for accurate cell type calling + GEX aggregation.

### Why This Should Work

From `PROTEIN_RECOVERY_REPORT.md`, scResolve multimodal-none achieves:
- Per-cell protein Pearson r: 0.65 (RNA-matched) to 0.82 (direct-matched)
- Per-cell protein cosine sim: 0.72 to 0.85

These strong correlations indicate per-cell protein profiles are well-preserved, so hierarchical gating should work effectively.

---

## Data Flow

```
scResolve multimodal-none output
    │
    ├── cells_exp_loc.h5ad
    │   ├── RNA features (369 genes)
    │   └── Protein_ features (27 markers)
    │
    ▼
┌─────────────────────────────────────┐
│ Step 1: Extract & Threshold Protein │
│ - Strip Protein_ prefix             │
│ - Calculate percentile thresholds   │
│   on scResolve's distribution       │
│   (handles 3.4x inflation)          │
└─────────────────────────────────────┘
    │
    ▼
┌─────────────────────────────────────┐
│ Step 2: Hierarchical Gating         │
│ - Same logic as create_protein_gt.py│
│ - CD20+ → B cells                   │
│ - CD3E+CD4+CD8A- → CD4+ T cells     │
│ - CD3E+CD8A+ → CD8+ T cells         │
│ - CD68+CD3E- → Macrophages          │
│ - CD31+CD68-CD3E- → Endothelial     │
│ - PanCK+|E-Cad high → Epithelial    │
│ - alphaSMA/Vim → Fibroblasts        │
└─────────────────────────────────────┘
    │
    ▼
┌─────────────────────────────────────┐
│ Step 3: Assign Cells to Spots       │
│ - KDTree nearest-neighbor on coords │
│ - Map scResolve cells → spot IDs    │
└─────────────────────────────────────┘
    │
    ▼
┌─────────────────────────────────────┐
│ Step 4: Aggregate RNA by Cell Type  │
│ - Sum RNA per cell type per spot    │
│ - Output: 7 CSVs (genes × spots)    │
└─────────────────────────────────────┘
    │
    ▼
┌─────────────────────────────────────┐
│ Step 5: Compare to GT GEX           │
│ - Same metrics as other methods     │
│ - Pearson r, RMSE, MAE per cell type│
└─────────────────────────────────────┘
```

---

## Gating Logic

Identical to `create_protein_gt.py` (lines 50-145):

### Thresholds (percentiles of non-zero values)

| Marker | Percentile | Rationale |
|--------|------------|-----------|
| CD3E | 50 | Standard |
| CD4 | 50 | Standard |
| CD8A | 50 | Standard |
| CD20 | 25 | Lower for sparse marker |
| CD68 | 50 | Standard |
| CD31 | 50 | Standard |
| PanCK | 25 | Lower for sparse marker |
| E-Cadherin | 90 | High for specificity |
| alphaSMA | 75 | Higher for specificity |
| Vimentin | 50 | Standard |

### Hierarchical Gates (order matters)

1. **B cells**: CD20+
2. **CD4+ T cells**: CD3E+ AND CD4+ AND CD8A- AND Unknown
3. **CD8+ T cells**: CD3E+ AND CD8A+ AND Unknown
4. **Macrophages**: CD68+ AND CD3E- AND Unknown
5. **Endothelial**: CD31+ AND CD68- AND CD3E- AND Unknown
6. **Epithelial**: (PanCK+ OR E-Cadherin high) AND Unknown
7. **Fibroblasts**: alphaSMA high AND NOT CD31+ AND CD68- AND CD3E- AND Unknown
8. **Fibroblasts**: Vimentin+ AND Unknown

---

## Implementation

### Task 1: Create gate_and_aggregate_gex.py

**File:** `Benchmarking/xenium_benchmarking/scResolve/src/gate_and_aggregate_gex.py`

```python
#!/usr/bin/env python
"""
Gate scResolve cells by protein markers and aggregate GEX per cell type.

Applies same hierarchical gating as create_protein_gt.py to scResolve's
multimodal output, then aggregates RNA per cell type per spot for
comparison with other deconvolution methods.

Usage:
    python gate_and_aggregate_gex.py --region 0 --output-dir /path/to/output
"""

import argparse
import json
import logging
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from scipy.spatial import cKDTree

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

# Paths
SCRESOLVE_DIR = Path("/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/"
                     "Benchmarking/xenium_benchmarking/scResolve")
PSEUDO_DIR = Path("/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/"
                  "Benchmarking/xenium_pseudovisium")

CELL_TYPE_ORDER = [
    "B cells", "CD4+ T cells", "CD8+ T cells", "Macrophages",
    "Endothelial", "Epithelial", "Fibroblasts",
]

GATING_PERCENTILES = {
    'CD3E': 50, 'CD4': 50, 'CD8A': 50, 'CD20': 25, 'CD68': 50,
    'CD31': 50, 'PanCK': 25, 'E-Cadherin': 90, 'alphaSMA': 75, 'Vimentin': 50,
}


def get_threshold(expr_df: pd.DataFrame, marker: str, percentile: float) -> float:
    """Get threshold as percentile of non-zero values."""
    if marker not in expr_df.columns:
        logger.warning(f"Marker {marker} not found in data")
        return 0
    vals = expr_df[marker]
    nonzero = vals[vals > 0]
    if len(nonzero) > 0:
        return np.percentile(nonzero, percentile)
    return 0


def classify_cells_by_protein(expr_df: pd.DataFrame) -> pd.Series:
    """Hierarchical protein gating - same logic as create_protein_gt.py."""
    cell_ids = expr_df.index
    cell_types = pd.Series('Unknown', index=cell_ids)

    # Calculate thresholds on THIS distribution
    thresholds = {m: get_threshold(expr_df, m, p) for m, p in GATING_PERCENTILES.items()}

    logger.info("Gating thresholds (on scResolve distribution):")
    logger.info(f"  CD3E: {thresholds['CD3E']:.1f}, CD4: {thresholds['CD4']:.1f}, "
                f"CD8A: {thresholds['CD8A']:.1f}")
    logger.info(f"  CD20: {thresholds['CD20']:.1f}, CD68: {thresholds['CD68']:.1f}, "
                f"CD31: {thresholds['CD31']:.1f}")
    logger.info(f"  PanCK: {thresholds['PanCK']:.1f}, E-Cad: {thresholds['E-Cadherin']:.1f}")
    logger.info(f"  alphaSMA: {thresholds['alphaSMA']:.1f}, Vimentin: {thresholds['Vimentin']:.1f}")

    # Gate 1: B cells (CD20+)
    b_cells = expr_df['CD20'] > thresholds['CD20']
    cell_types[b_cells] = 'B cells'
    logger.info(f"1. B cells (CD20+): {b_cells.sum()} ({b_cells.mean()*100:.1f}%)")

    # Gate 2: CD4+ T cells (CD3E+ CD4+ CD8A-)
    t_cell_base = expr_df['CD3E'] > thresholds['CD3E']
    cd4_pos = expr_df['CD4'] > thresholds['CD4']
    cd8_neg = expr_df['CD8A'] < thresholds['CD8A']
    cd4_tcells = t_cell_base & cd4_pos & cd8_neg & (cell_types == 'Unknown')
    cell_types[cd4_tcells] = 'CD4+ T cells'
    logger.info(f"2. CD4+ T cells: {cd4_tcells.sum()} ({cd4_tcells.mean()*100:.1f}%)")

    # Gate 3: CD8+ T cells (CD3E+ CD8A+)
    cd8_pos = expr_df['CD8A'] > thresholds['CD8A']
    cd8_tcells = t_cell_base & cd8_pos & (cell_types == 'Unknown')
    cell_types[cd8_tcells] = 'CD8+ T cells'
    logger.info(f"3. CD8+ T cells: {cd8_tcells.sum()} ({cd8_tcells.mean()*100:.1f}%)")

    # Gate 4: Macrophages (CD68+ CD3E-)
    cd68_pos = expr_df['CD68'] > thresholds['CD68']
    cd3e_neg = expr_df['CD3E'] < thresholds['CD3E']
    macrophages = cd68_pos & cd3e_neg & (cell_types == 'Unknown')
    cell_types[macrophages] = 'Macrophages'
    logger.info(f"4. Macrophages: {macrophages.sum()} ({macrophages.mean()*100:.1f}%)")

    # Gate 5: Endothelial (CD31+ CD68- CD3E-)
    cd31_pos = expr_df['CD31'] > thresholds['CD31']
    cd68_neg = expr_df['CD68'] < thresholds['CD68']
    endothelial = cd31_pos & cd68_neg & cd3e_neg & (cell_types == 'Unknown')
    cell_types[endothelial] = 'Endothelial'
    logger.info(f"5. Endothelial: {endothelial.sum()} ({endothelial.mean()*100:.1f}%)")

    # Gate 6: Epithelial (PanCK+ or E-Cadherin high)
    panck_pos = expr_df['PanCK'] > thresholds['PanCK']
    ecad_high = expr_df['E-Cadherin'] > thresholds['E-Cadherin']
    epithelial = (panck_pos | ecad_high) & (cell_types == 'Unknown')
    cell_types[epithelial] = 'Epithelial'
    logger.info(f"6. Epithelial: {epithelial.sum()} ({epithelial.mean()*100:.1f}%)")

    # Gate 7: Fibroblasts (alphaSMA high)
    asma_high = expr_df['alphaSMA'] > thresholds['alphaSMA']
    fibroblasts_asma = asma_high & ~cd31_pos & cd68_neg & cd3e_neg & (cell_types == 'Unknown')
    cell_types[fibroblasts_asma] = 'Fibroblasts'
    logger.info(f"7. Fibroblasts (alphaSMA): {fibroblasts_asma.sum()} ({fibroblasts_asma.mean()*100:.1f}%)")

    # Gate 8: Fibroblasts (Vimentin+ remaining)
    vim_pos = expr_df['Vimentin'] > thresholds['Vimentin']
    fibroblasts_vim = vim_pos & (cell_types == 'Unknown')
    cell_types[fibroblasts_vim] = 'Fibroblasts'
    logger.info(f"8. Fibroblasts (Vimentin): {fibroblasts_vim.sum()} ({fibroblasts_vim.mean()*100:.1f}%)")

    logger.info(f"Unknown: {(cell_types == 'Unknown').sum()} ({(cell_types == 'Unknown').mean()*100:.1f}%)")

    return cell_types


def load_scresolve_cells(region: int):
    """Load scResolve multimodal-none output."""
    region_name = f"Xenium_region_{region}"
    cells_path = (SCRESOLVE_DIR / "output_multimodal_v2" /
                  f"{region_name}_multimodal_none" / "segmentation" / "cells_exp_loc.h5ad")

    if not cells_path.exists():
        raise FileNotFoundError(f"scResolve cells not found: {cells_path}")

    adata = sc.read_h5ad(cells_path)
    logger.info(f"Loaded {adata.n_obs} cells, {adata.n_vars} features")

    return adata


def extract_protein_and_rna(adata):
    """Split features into protein and RNA."""
    all_names = np.array(adata.var_names)
    protein_mask = np.array([v.startswith("Protein_") for v in all_names])
    rna_mask = ~protein_mask

    X = adata.X.toarray() if sparse.issparse(adata.X) else adata.X.copy()

    # Protein DataFrame (strip Protein_ prefix)
    protein_names = [v.replace("Protein_", "") for v in all_names[protein_mask]]
    protein_df = pd.DataFrame(X[:, protein_mask], index=adata.obs_names, columns=protein_names)

    # RNA DataFrame
    rna_names = all_names[rna_mask]
    rna_df = pd.DataFrame(X[:, rna_mask], index=adata.obs_names, columns=rna_names)

    logger.info(f"  Protein features: {len(protein_names)}")
    logger.info(f"  RNA features: {len(rna_names)}")

    return protein_df, rna_df


def assign_cells_to_spots(adata, spot_coords: np.ndarray, spot_ids: np.ndarray):
    """Assign scResolve cells to nearest spots using KDTree."""
    cell_coords = adata.obs[['x', 'y']].values

    spot_tree = cKDTree(spot_coords)
    _, spot_indices = spot_tree.query(cell_coords)

    cell_spot_ids = spot_ids[spot_indices]

    return pd.Series(cell_spot_ids, index=adata.obs_names)


def aggregate_rna_by_celltype(
    rna_df: pd.DataFrame,
    cell_types: pd.Series,
    cell_spot_ids: pd.Series,
    spot_ids: np.ndarray,
) -> dict:
    """Aggregate RNA per cell type per spot. Returns dict of DataFrames (genes x spots)."""

    gex_layers = {}

    for ct in CELL_TYPE_ORDER:
        # Get cells of this type
        ct_mask = cell_types == ct
        ct_cells = cell_types[ct_mask].index

        if len(ct_cells) == 0:
            logger.warning(f"  No cells for {ct}")
            continue

        # Initialize spot x gene matrix
        ct_gex = pd.DataFrame(0.0, index=spot_ids, columns=rna_df.columns)

        # Sum RNA for each spot
        for cell_id in ct_cells:
            spot_id = cell_spot_ids[cell_id]
            ct_gex.loc[spot_id] += rna_df.loc[cell_id]

        # Transpose to genes x spots (matching GT format)
        gex_layers[ct] = ct_gex.T

        n_active_spots = (ct_gex.sum(axis=1) > 0).sum()
        logger.info(f"  {ct}: {len(ct_cells)} cells, {n_active_spots} active spots")

    return gex_layers


def main():
    parser = argparse.ArgumentParser(description="Gate scResolve cells and aggregate GEX")
    parser.add_argument("--region", type=int, required=True, help="Region ID (0-4)")
    parser.add_argument("--output-dir", type=str, default=None,
                        help="Output directory for layers")
    args = parser.parse_args()

    region = args.region
    region_name = f"Xenium_region_{region}"

    # Output directory
    if args.output_dir:
        output_dir = Path(args.output_dir)
    else:
        output_dir = SCRESOLVE_DIR / "output_protein_gated" / region_name / "layers"
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("=" * 60)
    logger.info(f"scResolve Protein-Gated GEX Aggregation - {region_name}")
    logger.info("=" * 60)

    # Step 1: Load scResolve cells
    logger.info("\n[Step 1] Loading scResolve cells...")
    adata = load_scresolve_cells(region)

    # Step 2: Extract protein and RNA
    logger.info("\n[Step 2] Extracting protein and RNA features...")
    protein_df, rna_df = extract_protein_and_rna(adata)

    # Step 3: Gate cells by protein
    logger.info("\n[Step 3] Gating cells by protein markers...")
    cell_types = classify_cells_by_protein(protein_df)

    # Step 4: Load spot coordinates
    logger.info("\n[Step 4] Loading spot coordinates...")
    spots_adata = sc.read_h5ad(PSEUDO_DIR / "data_protein_gt/h5ad_objects" / f"{region_name}_GEX.h5ad")
    spot_coords = spots_adata.obsm['spatial']
    spot_ids = spots_adata.obs_names.values
    logger.info(f"  {len(spot_ids)} spots")

    # Step 5: Assign cells to spots
    logger.info("\n[Step 5] Assigning cells to spots...")
    cell_spot_ids = assign_cells_to_spots(adata, spot_coords, spot_ids)

    # Step 6: Aggregate RNA by cell type
    logger.info("\n[Step 6] Aggregating RNA by cell type...")
    gex_layers = aggregate_rna_by_celltype(rna_df, cell_types, cell_spot_ids, spot_ids)

    # Step 7: Save layers
    logger.info("\n[Step 7] Saving layers...")
    for ct, df in gex_layers.items():
        ct_filename = ct.replace(" ", "_").replace("+", "+")
        out_path = output_dir / f"{ct_filename}_layer.csv"
        df.to_csv(out_path)
        logger.info(f"  Saved: {out_path}")

    # Save summary
    summary = {
        "region": region,
        "n_cells": len(cell_types),
        "n_spots": len(spot_ids),
        "cell_type_counts": cell_types.value_counts().to_dict(),
        "n_rna_genes": len(rna_df.columns),
        "n_protein_markers": len(protein_df.columns),
    }
    with open(output_dir / "summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    logger.info("\n" + "=" * 60)
    logger.info("Done!")
    logger.info(f"Output: {output_dir}")
    logger.info("=" * 60)


if __name__ == "__main__":
    main()
```

---

### Task 2: Create SLURM submission script

**File:** `Benchmarking/xenium_benchmarking/scResolve/slurm/gate_and_aggregate_gex.sh`

```bash
#!/bin/bash
#SBATCH --job-name=scresolve_gate
#SBATCH --output=slurm_log/gate_gex_%A_%a.out
#SBATCH --error=slurm_log/gate_gex_%A_%a.err
#SBATCH --array=0-4
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

set -e

REPO="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist"
export PATH="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env/bin:$PATH"
source activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

REGION=${SLURM_ARRAY_TASK_ID}

cd "${REPO}/Benchmarking/xenium_benchmarking/scResolve"
mkdir -p slurm_log

python src/gate_and_aggregate_gex.py --region ${REGION}
```

---

### Task 3: Update compare_all_methods_gex.py

**File:** `Benchmarking/xenium_benchmarking/evaluation/src/compare_all_methods_gex.py`

Add scResolve loader and include in methods dict:

```python
def load_scresolve_gex(region_id: int) -> Dict[str, pd.DataFrame]:
    """Load scResolve protein-gated GEX layers."""
    layers_dir = BASE_DIR / "scResolve" / "output_protein_gated" / f"Xenium_region_{region_id}" / "layers"

    if not layers_dir.exists():
        return {}

    ct_dfs = {}
    for ct in ACHIEVABLE_7_CELL_TYPES:
        ct_filename = ct.replace(" ", "_").replace("+", "+")
        layer_file = layers_dir / f"{ct_filename}_layer.csv"
        if layer_file.exists():
            ct_dfs[ct] = pd.read_csv(layer_file, index_col=0)

    return ct_dfs


# Update methods dict (around line 169):
methods = {
    "CITEgeist": lambda rid: load_citegeist_gex(rid),
    "Cell2Location": lambda rid: load_competitor_gex("Cell2Location", rid),
    "Tangram": lambda rid: load_competitor_gex("Tangram", rid),
    "scResolve": lambda rid: load_scresolve_gex(rid),
}
```

---

## Expected Output

### Directory Structure

```
Benchmarking/xenium_benchmarking/scResolve/
├── output_protein_gated/
│   ├── Xenium_region_0/
│   │   └── layers/
│   │       ├── B_cells_layer.csv
│   │       ├── CD4+_T_cells_layer.csv
│   │       ├── CD8+_T_cells_layer.csv
│   │       ├── Macrophages_layer.csv
│   │       ├── Endothelial_layer.csv
│   │       ├── Epithelial_layer.csv
│   │       ├── Fibroblasts_layer.csv
│   │       └── summary.json
│   ├── Xenium_region_1/
│   │   └── layers/...
│   └── ...
```

### Expected Results

Based on protein recovery metrics (r=0.65-0.82):

| Cell Type | Expected r | Reasoning |
|-----------|-----------|-----------|
| Macrophages | 0.55-0.60 | CD68 profile preserved |
| Fibroblasts | 0.50-0.55 | alphaSMA/Vim profiles preserved |
| CD8+ T cells | 0.45-0.50 | CD3E+CD8A pattern preserved |
| Epithelial | 0.40-0.50 | PanCK/E-Cad distinctive |
| Endothelial | 0.40-0.45 | CD31 profile intact |
| CD4+ T cells | 0.25-0.35 | CD4 vs CD8 discrimination trickier |
| B cells | 0.15-0.25 | CD20 sparse |
| **Mean** | **0.40-0.45** | **Competitive with CITEgeist/C2L** |

---

## Execution Order

1. **Create script:** Write `gate_and_aggregate_gex.py`
2. **Create SLURM:** Write `gate_and_aggregate_gex.sh`
3. **Submit jobs:** `sbatch slurm/gate_and_aggregate_gex.sh`
4. **Wait for completion:** Monitor with `squeue -u alc376`
5. **Update comparison:** Add scResolve loader to `compare_all_methods_gex.py`
6. **Run comparison:** `python evaluation/src/compare_all_methods_gex.py`
7. **Analyze results:** Compare scResolve to CITEgeist/C2L

---

## Validation Checklist

- [ ] Gating thresholds match create_protein_gt.py percentiles
- [ ] Gating logic matches create_protein_gt.py hierarchy
- [ ] Output format matches CITEgeist layers (genes x spots CSV)
- [ ] All 5 regions processed successfully
- [ ] Cell type distribution is reasonable (not all Unknown)
- [ ] scResolve appears in comparison results

---

## References

- `Benchmarking/xenium_pseudovisium/src/create_protein_gt.py` - GT gating logic
- `Benchmarking/xenium_benchmarking/scResolve/comparison_results/PROTEIN_RECOVERY_REPORT.md` - Protein recovery analysis
- `Benchmarking/xenium_benchmarking/evaluation/src/compare_all_methods_gex.py` - GEX comparison framework
