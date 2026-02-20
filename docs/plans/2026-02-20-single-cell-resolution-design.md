# CITEgeist Single-Cell Resolution Design

**Date:** 2026-02-20
**Status:** Approved
**Author:** Design session with Claude

## Overview

Extend CITEgeist from spot-level deconvolution to true single-cell resolution by assigning individual nuclei (detected via Cellpose) to cell types using a combination of:
1. Antibody-based proportion estimation (existing Module 3)
2. Nuclear morphology features (new)
3. Optimal assignment via Hungarian algorithm (new)

## Motivation

Current CITEgeist outputs:
- Cell type **proportions** per spot (e.g., "40% macrophage, 35% fibroblast, 25% cancer")
- Deconvolved GEX as **compartment layers** (e.g., `barcode:::macrophage`)

This is spot-level resolution (~55µm diameter, ~5-15 cells per spot).

Proposed extension:
- Assign each **individual nucleus** a cell type identity
- Output true **single-cell expression matrices**
- Enable downstream single-cell analyses on spatial data

## Architecture: Module 3 Refactor

Module 3 becomes three composable sub-modules:

```
┌─────────────────────────────────────────────────────────────────────┐
│  MODULE 3a: Proportion Estimation                                   │
│                                                                      │
│  Input:  Antibody data, cell profiles, (nuclei counts optional)     │
│  Output: Cell type proportions per spot [0-1, sum to 1]             │
│  Method: Continuous optimization (existing QP/EM)                   │
│  File:   module3a_proportion.py                                     │
└─────────────────────────────────────────────────────────────────────┘
                               ↓
┌─────────────────────────────────────────────────────────────────────┐
│  MODULE 3b: Per-Nucleus Assignment (NEW)                            │
│                                                                      │
│  Input:  Proportions, Cellpose masks, nuclei count                  │
│  Output: Cell type identity for each individual nucleus             │
│  Method: Soft-label classifier + Hungarian matching                 │
│  File:   module3b_nucleus_assignment.py                             │
└─────────────────────────────────────────────────────────────────────┘
                               ↓
┌─────────────────────────────────────────────────────────────────────┐
│  MODULE 3c: GEX Deconvolution                                       │
│                                                                      │
│  Input:  GEX data + assignments (spot-level OR cell-level)          │
│  Output: Deconvolved expression                                     │
│          - Spot mode: barcode:::cell_type layers (existing)         │
│          - Cell mode: per-nucleus expression matrix (new)           │
│  File:   module3c_gex_deconvolution.py                              │
└─────────────────────────────────────────────────────────────────────┘
```

### Backwards Compatibility

| Use case | Sub-modules | Output |
|----------|-------------|--------|
| Proportions only | 3a | Spot-level proportions |
| Current pipeline | 3a → 3c | Compartment layers |
| Single-cell resolution | 3a → 3b → 3c | Per-cell AnnData |

### Why Continuous for 3a?

Benchmark results (Xenium pseudo-Visium):

| Method | Proportion Pearson r |
|--------|---------------------|
| Continuous | 0.911 |
| Hybrid (continuous → discretize) | 0.852 |
| Discrete IQP | 0.435-0.632 |

Continuous optimization finds better proportions. Module 3b handles discretization via morphology-informed assignment, which is more intelligent than simple largest-remainder rounding.

## Module 3b: Per-Nucleus Assignment

### Step 1: Morphology Feature Extraction

Extract features from Cellpose label masks using `skimage.measure.regionprops`:

| Feature | Formula | Biological Relevance |
|---------|---------|---------------------|
| Area | pixel count | Cell size (cancer often larger) |
| Perimeter | boundary length | Complexity |
| Circularity | 4π × area / perimeter² | Round (lymphocytes) vs elongated (fibroblasts) |
| Eccentricity | 1 - (minor/major)² | Elongation |
| Solidity | area / convex_hull_area | Irregular boundaries (cancer) |
| Major axis | fitted ellipse | Size in longest dimension |
| Minor axis | fitted ellipse | Size in shortest dimension |
| Aspect ratio | major / minor | Shape |

**Output:** DataFrame with `[nucleus_id, spot_id, x, y, area, circularity, ...]`

### Step 2: Soft-Label Classifier Training

**Key insight:** We don't know individual nucleus labels, but we know spot-level proportions.

**Training data construction:**
```python
# For each nucleus n in spot s:
X[n] = [area, circularity, eccentricity, ...]  # 8 features
y[n] = proportions[s]  # soft label, e.g., [0.6, 0.3, 0.1, 0, ...]
```

**Why this works:**
- Spots dominated by macrophages will have mostly macrophage-like nuclei
- The classifier learns morphology↔type correlations across thousands of nuclei
- Individual labels are noisy, but aggregate signal emerges statistically

**Model:** Multinomial logistic regression (sklearn)
- Simple, interpretable, fast
- L2 regularization
- Outputs calibrated probabilities
- Can inspect coefficients (e.g., "larger area → higher cancer probability")

```python
from sklearn.linear_model import LogisticRegression

clf = LogisticRegression(multi_class='multinomial', max_iter=1000)
# Train with sample weights derived from soft labels
clf.fit(X_expanded, y_expanded, sample_weight=weights)
```

### Step 3: Proportion to Count Conversion

```python
# 3a outputs proportions (sum to 1)
proportions = [0.40, 0.35, 0.25, 0, 0, 0, 0]  # 7 cell types

# Cellpose outputs nuclei count
nuclei_count = 5

# Convert to counts
raw_counts = proportions * nuclei_count  # [2.0, 1.75, 1.25, 0, 0, 0, 0]

# Discretize using largest remainder method
discrete_counts = largest_remainder(raw_counts)  # [2, 2, 1, 0, 0, 0, 0]
```

### Step 4: Hungarian Assignment

**Per-spot procedure:**

```python
from scipy.optimize import linear_sum_assignment

def assign_nuclei_to_types(nuclei_features, discrete_counts, classifier):
    # Get probabilities from classifier
    probs = classifier.predict_proba(nuclei_features)  # (N_nuclei, K_types)

    # Expand counts to slots
    # [2 mac, 2 fib, 1 cancer] → ['mac', 'mac', 'fib', 'fib', 'cancer']
    slots = []
    slot_types = []
    for t, count in enumerate(discrete_counts):
        slots.extend([t] * count)

    # Build cost matrix: -log(prob) for numerical stability
    n_nuclei = len(nuclei_features)
    n_slots = len(slots)
    cost_matrix = np.zeros((n_nuclei, n_slots))
    for i in range(n_nuclei):
        for j, t in enumerate(slots):
            cost_matrix[i, j] = -np.log(probs[i, t] + 1e-10)

    # Solve assignment
    row_ind, col_ind = linear_sum_assignment(cost_matrix)

    # Return assignments
    assignments = {nucleus_id[row_ind[i]]: cell_types[slots[col_ind[i]]]
                   for i in range(len(row_ind))}
    return assignments
```

### Edge Cases

**Few nuclei:**
- Proportions = [0.4, 0.35, 0.25], Cellpose = 1 nucleus
- Counts = [0.4, 0.35, 0.25] → [1, 0, 0]
- Single nucleus assigned to highest proportion type

**Rounding artifacts:**
- Proportions = [0.33, 0.33, 0.34], Cellpose = 2 nuclei
- Largest remainder handles deterministically

## Module 3c: GEX Distribution

### Spot-Level Mode (Existing)

Output: `barcode:::cell_type` layers in parquet format

### Cell-Level Mode (New)

```python
for spot_id in spots:
    for cell_type in cell_types:
        # Get nuclei assigned to this type
        nuclei = assignments[(spot_id, cell_type)]

        if len(nuclei) == 0:
            continue

        # Get deconvolved expression for compartment
        layer_expr = deconvolved_gex[f"{spot_id}:::{cell_type}"]

        # Equal split among cells
        per_cell_expr = layer_expr / len(nuclei)

        for nucleus in nuclei:
            cell_expression[nucleus] = per_cell_expr
```

## Output Format

### Single-Cell AnnData

```python
adata_singlecell = AnnData(
    X = cell_expression_matrix,  # (N_cells, N_genes)
    obs = pd.DataFrame({
        'cell_id': [...],           # unique nucleus ID
        'spot_id': [...],           # parent spot barcode
        'cell_type': [...],         # assigned type
        'x': [...],                 # centroid x (Cellpose)
        'y': [...],                 # centroid y (Cellpose)
        'area': [...],              # morphology features
        'circularity': [...],
        'eccentricity': [...],
        'solidity': [...],
        'major_axis': [...],
        'minor_axis': [...],
        'aspect_ratio': [...],
    }),
    var = gene_metadata,
    uns = {
        'spatial': {...},
        'morphology_classifier': clf,
        'assignment_method': 'hungarian',
        'source_sample': sample_name,
    }
)
```

**File:** `{sample}_single_cell.h5ad`

## File Structure

```
CITEgeist/model/
├── module3a_proportion.py        # Refactored from gurobi_impl.py
├── module3b_nucleus_assignment.py    # NEW
├── module3c_gex_deconvolution.py     # Refactored from existing
├── morphology_features.py        # NEW: feature extraction
└── citegeist_model.py            # Orchestrates sub-modules
```

## Known Limitations

1. **Morphologically similar subtypes:** Cell types with identical nuclear morphology (e.g., CD4+ vs CD8+ T cells) will be assigned based on counts but not spatially resolved within a spot. The antibody data distinguishes them at the spot level, but we cannot determine which specific nucleus is CD4 vs CD8.

2. **Equal GEX split:** All cells of the same type in a spot get equal expression. No cell-to-cell variation within type within spot.

3. **Cellpose accuracy:** Assignment quality depends on segmentation quality. Under-segmentation (merged nuclei) or over-segmentation (split nuclei) will propagate to assignments.

## Dependencies

- `scikit-image`: `regionprops` for morphology features
- `scikit-learn`: `LogisticRegression` for classifier
- `scipy`: `linear_sum_assignment` for Hungarian algorithm
- Existing: `cellpose`, `gurobi`, `anndata`

## Success Criteria

1. Pipeline runs end-to-end on Xenium benchmark data
2. Output AnnData passes validation (correct dimensions, no NaN)
3. Cell type assignments show spatial coherence when visualized
4. Morphology classifier shows interpretable feature importances
5. Downstream single-cell analyses (clustering, DE) produce sensible results

## Future Extensions (Out of Scope)

- Morphology + intensity features (texture within nucleus)
- Stochastic GEX variation (sample from learned variance)
- Iterative classifier refinement
- Integration with Xenium cell boundaries (when available)
