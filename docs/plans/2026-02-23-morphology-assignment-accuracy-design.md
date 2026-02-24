# Morphology Assignment Accuracy Assessment Design

**Date:** 2026-02-23
**Status:** Approved
**Goal:** Evaluate whether Cellpose nuclei morphology features improve cell type assignment accuracy in Module 3b

## Overview

The single-cell resolution pipeline (Module 3b) uses morphology features from Cellpose segmentation to guide per-nucleus cell type assignment. This assessment will determine:
1. Whether morphology-guided assignment outperforms baseline methods
2. Cell-level accuracy against ground truth labels

## Architecture

**New file:** `Benchmarking/xenium_benchmarking/evaluation/src/evaluate_morphology_assignment.py`

**Inputs:**
- Single-cell AnnData from Module 3b (`*_single_cell.h5ad`)
- Protein-gated ground truth (`data_protein_gt/cell_type_assignments.csv`)
- RNA clustering ground truth (`data_rna_gt/cell_type_assignments.csv`)

**Outputs:**
- `morphology_evaluation_results.json` — all metrics
- `confusion_matrix_protein_gt.png` — visualization
- `confusion_matrix_rna_gt.png`
- `baseline_comparison_table.csv` — side-by-side comparison
- `per_celltype_f1_barplot.png` — which types benefit most

**Flow:**
```
Load single-cell AnnData
    ↓
Load ground truth (Protein GT + RNA GT)
    ↓
Match Cellpose nuclei to GT cells (KDTree, <10µm)
    ↓
For each assignment method:
    ├─ Morphology-guided (existing)
    ├─ Random baseline (10 seeds)
    ├─ Uniform probability baseline
    └─ Spot-proportion-only baseline
    ↓
Compute metrics (accuracy, P/R/F1, confusion matrix)
    ↓
Output comparison report
```

## Ground Truth Sources

### Protein-Gated GT (7 cell types)
- Source: `data_protein_gt/cell_type_assignments.csv`
- Hierarchical gating on protein markers
- Cell types: B cells, CD4+ T cells, CD8+ T cells, Macrophages, Endothelial, Epithelial, Fibroblasts
- Has "Unknown" cells (excluded from accuracy calculation)

**Gating logic:**
| Gate | Logic | Cell Type |
|------|-------|-----------|
| 1 | CD20+ | B cells |
| 2 | CD3E+ CD4+ CD8A- | CD4+ T cells |
| 3 | CD3E+ CD8A+ | CD8+ T cells |
| 4 | CD68+ CD3E- | Macrophages |
| 5 | CD31+ CD68- CD3E- | Endothelial |
| 6 | PanCK+ or E-Cadherin high | Epithelial |
| 7 | alphaSMA high | Fibroblasts |

### RNA Clustering GT (6 cell types)
- Source: `data_rna_gt/cell_type_assignments.csv`
- Xenium k-means clustering on gene expression
- Cell types: B cells, T cells (combined), Macrophages, Endothelial, Epithelial, Fibroblasts
- ALL cells assigned (no Unknown)

**Evaluation approach:**
- Run accuracy metrics against BOTH ground truths
- For Protein GT: evaluate on 7 cell types directly
- For RNA GT: collapse CD4+/CD8+ T cell predictions → "T cells" before comparison
- Report concordance between the two GTs as sanity check

**Matching to Cellpose nuclei:**
- Both GTs have cell_ids from Xenium native segmentation
- Match to Cellpose nuclei by spatial proximity (KDTree nearest neighbor)
- Threshold: 10µm max distance
- Report match rate as diagnostic

## Baseline Assignment Methods

### Morphology-Guided (existing)
- Uses trained soft-label classifier on 7 morphology features (area, circularity, eccentricity, solidity, major_axis, minor_axis, aspect_ratio)
- Hungarian assignment respects per-spot count constraints
- **This is what we're evaluating**

### Random Baseline
- For each spot: randomly shuffle cell type assignments among nuclei
- Respects count constraints (same number of each type as morphology-guided)
- Run 10 random seeds, report mean ± std
- **Tests:** Does morphology beat pure chance?

### Uniform Probability Baseline
- All nuclei get equal probability (1/K) for each cell type
- Hungarian algorithm picks assignments based on count constraints alone
- **Tests:** Does the classifier add value beyond Hungarian matching?

### Spot-Proportion-Only Baseline
- Each nucleus gets spot's proportion vector as its probability
- No morphology features used — just spot context
- Hungarian assignment as usual
- **Tests:** Does morphology add anything beyond knowing "this spot is 60% macrophages"?

## Evaluation Metrics

### Per-Cell Accuracy
- Overall accuracy: % of cells correctly assigned
- Computed separately for each GT (Protein, RNA)

### Per-Cell-Type Metrics
For each cell type:
- **Precision**: TP / (TP + FP) — of cells we called type X, how many are correct?
- **Recall**: TP / (TP + FN) — of actual type X cells, how many did we find?
- **F1**: harmonic mean of precision and recall
- Macro-averaged F1 across all types

### Confusion Matrix
- Full K×K matrix (predicted vs actual)
- Normalized by row (actual) to show recall per type
- Highlights which cell types get confused with which
- Separate matrices for Protein GT and RNA GT

### Summary Statistics
- **Morphology lift**: (morphology accuracy - best baseline accuracy) / best baseline accuracy
- **Per-region breakdown**: metrics computed per Xenium region (0,1,2,3,4)
- **Aggregate**: mean ± std across regions

## Implementation Structure

```python
# Core functions
def load_ground_truth(gt_type: str, region_id: int) -> pd.DataFrame
def match_cellpose_to_gt(cellpose_coords, gt_coords, max_dist=10) -> Dict[int, str]
def run_baseline_random(spot_assignments, n_seeds=10) -> List[Dict]
def run_baseline_uniform(spot_props, nuclei_per_spot) -> Dict
def run_baseline_spot_proportion(spot_props, nuclei_per_spot) -> Dict
def compute_accuracy_metrics(pred_labels, gt_labels) -> Dict
def compute_confusion_matrix(pred_labels, gt_labels, cell_types) -> np.ndarray
def plot_confusion_matrix(cm, cell_types, output_path)
def plot_baseline_comparison(results, output_path)

# Main entry point
def evaluate_region(region_id, sc_adata_path, gt_dir) -> Dict
def main(args)
```

**CLI Interface:**
```bash
python evaluate_morphology_assignment.py \
    --sc-dir /path/to/single_cell_outputs \
    --gt-dir /path/to/xenium_pseudovisium \
    --regions 0,1,2,3,4 \
    --output-dir /path/to/results
```

**SLURM wrapper:** `sbatch_evaluate_morphology.sh` (10 min, 32GB)

## Expected Outcomes

**Hypothesis:** Morphology features provide modest improvement over spot-proportion-only baseline for certain cell types (e.g., larger Macrophages vs smaller lymphocytes) but limited improvement for morphologically similar types.

**Success criteria:**
- Morphology-guided accuracy > all three baselines
- At least one cell type shows >5% F1 improvement from morphology
- Results consistent across Protein GT and RNA GT

## Future Extensions (out of scope)

- Feature importance analysis (which morphology features matter most)
- Alternative classifiers (random forest, neural network)
- Cross-region generalization testing
