# Cell Morphology Feature Extraction via Watershed Segmentation

**Date:** 2026-02-24
**Status:** Approved
**Goal:** Improve Module 3b single-cell type assignment by extracting cell-level morphology features (not just nuclear)

## Background

### Problem Statement

CITEgeist Module 3b assigns individual nuclei to cell types using:
1. Cellpose nuclear segmentation on DAPI
2. Soft-label classifier trained on spot proportions
3. Hungarian assignment respecting spot-level cell counts

**Current results (4 regions, 56,741 cells):**

| Method | Accuracy | Macro F1 |
|--------|----------|----------|
| Morphology-guided | 28.4% | 23.5% |
| Random baseline | 28.5% | 23.6% |
| Spot-proportion | 28.6% | 23.6% |

Nuclear morphology features (area, circularity, eccentricity, solidity, aspect_ratio) provide **zero discriminative signal** - accuracy equals random baseline.

### Hypothesis

Cell-level morphology features may be more discriminative than nuclear features alone:
- Macrophages: large cells, low N:C ratio
- T-cells: small cells, high N:C ratio
- Epithelial: polygonal cells, variable N:C ratio
- Fibroblasts: elongated cells

### Future Goal

Develop a general method that works on H&E images:
- Hematoxylin stains nuclei (analogous to DAPI)
- Eosin stains cytoplasm (analogous to boundary channel)

## Data Available

### Xenium RCC Dataset

Location: `/ix1/alee/LO_LAB/General/Public_Data/10x_PublicData/Xenium_RNA_Proteomic_RenalCellCarcinoma/`

**Relevant channels:**
- `ch0000_dapi.ome.tif` - Nuclear stain (for Cellpose)
- `ch0001_atp1a1_cd45_e-cadherin.ome.tif` - Boundary/membrane stain (for watershed)

**Note:** Boundary channel is a combination stain:
- ATP1A1: Ubiquitous membrane marker (Na+/K+ ATPase)
- CD45: Leukocyte-specific membrane marker
- E-cadherin: Epithelial cell-cell junction marker

This means signal intensity varies by cell type, but edges should be visible for most cells.

## Design

### Pipeline Overview

```
DAPI Image ──► Cellpose ──► Nuclear masks + centroids
                                  │
Boundary Channel ──► Gradient ────┼──► Watershed ──► Cell masks
                                  │        ▲
                           (nuclei as seeds)
                                  │
                                  ▼
                    ┌─────────────────────────┐
                    │   Feature Extraction    │
                    │   - Nuclear features    │
                    │   - Cell features       │
                    │   - Ratio features      │
                    └─────────────────────────┘
                                  │
                                  ▼
                    ┌─────────────────────────┐
                    │   Soft-label Classifier │
                    │   + Hungarian Assignment│
                    └─────────────────────────┘
```

### Step 1: Nuclear Segmentation (Existing)

Use existing Cellpose pipeline on DAPI:
- Model: `nuclei`
- Output: Label mask where each nucleus has unique ID
- Extract centroids for each nucleus

### Step 2: Watershed Cell Segmentation (New)

**Input:**
- Nuclear mask from Cellpose (seeds)
- Boundary channel image (gradient source)

**Algorithm:**
1. Load boundary channel image
2. Compute gradient magnitude (Sobel filter)
3. Use nuclear centroids as watershed seeds (markers)
4. Run `skimage.segmentation.watershed` with gradient as elevation map
5. Output: Cell mask where each cell has same ID as its nucleus

**Edge cases:**
- Nuclei at image boundary: Allow watershed to extend to image edge
- Touching nuclei: Each gets separate watershed basin (guaranteed by seeding)
- Dim boundary regions: Watershed will find nearest edges (may over-segment)

### Step 3: Feature Extraction (Extended)

**Current nuclear features (5):**
- `area`: Nuclear area in pixels
- `circularity`: 4π × area / perimeter²
- `eccentricity`: Ellipse eccentricity (0=circle, 1=line)
- `solidity`: area / convex_hull_area
- `aspect_ratio`: major_axis / minor_axis

**New cell features (5):**
- `cell_area`: Total cell area in pixels
- `cell_circularity`: 4π × cell_area / cell_perimeter²
- `cell_eccentricity`: Cell ellipse eccentricity
- `cell_solidity`: cell_area / cell_convex_hull_area
- `cell_aspect_ratio`: Cell major_axis / minor_axis

**New ratio features (2):**
- `nc_ratio`: nucleus_area / cell_area (N:C ratio)
- `cytoplasm_area`: cell_area - nucleus_area

**Total: 12 features** (vs current 5)

### Step 4: Classification (Modified)

Modify `SoftLabelClassifier` to use expanded feature set:
- Input: 12 features per nucleus/cell
- Training: Soft labels from spot proportions (unchanged)
- Output: Cell type probabilities per cell

Hungarian assignment unchanged - uses probability matrix + spot counts.

### Implementation Decisions

| Decision | Choice | Rationale |
|----------|--------|-----------|
| Watershed library | `skimage.segmentation.watershed` | Standard, well-documented, handles seeds properly |
| Gradient computation | Sobel magnitude | Robust edge detection, works on 16-bit images |
| Marker generation | Binary dilation of nuclear centroids | Ensures each nucleus seeds exactly one basin |
| Feature normalization | Per-spot z-score | Accounts for local tissue/staining variation |
| Missing cells | Exclude from training/eval | Cells without valid watershed = no features |

### File Structure

New/modified files:
```
CITEgeist/model/
├── watershed_segmentation.py     # NEW: Watershed cell segmentation
├── cell_morphology_features.py   # NEW: Cell-level feature extraction
├── morphology_features.py        # MODIFY: Add cell feature support
├── soft_label_classifier.py      # MODIFY: Handle 12 features
└── module3b_nucleus_assignment.py # MODIFY: Use cell features

Benchmarking/xenium_benchmarking/
├── CITEgeist/src/
│   └── benchmark_cell_morphology.py  # NEW: Benchmark script
└── evaluation/src/
    └── evaluate_cell_morphology.py   # NEW: Evaluation vs baselines
```

## Success Criteria

**Primary metric:** Classification accuracy on Xenium ground truth

| Outcome | Interpretation |
|---------|----------------|
| Accuracy > 35% | Cell features provide meaningful signal; continue development |
| Accuracy 30-35% | Marginal improvement; may not justify complexity |
| Accuracy ~28% | Cell morphology also not discriminative; accept uncertainty |

**Secondary metrics:**
- Per-cell-type F1 scores (some types may be more distinguishable)
- Feature importance analysis (which features matter most)

## Risks and Mitigations

| Risk | Mitigation |
|------|------------|
| Boundary channel gives poor segmentation | Fall back to Voronoi from nuclei; report both |
| Variable staining across regions | Per-region feature normalization |
| Cell features also not discriminative | Accept that single-cell assignment has inherent uncertainty given only morphology |
| Doesn't transfer to H&E | Test on H&E data before claiming generalization |

## Out of Scope

- Using boundary channel intensity as a feature (would not transfer to H&E)
- Deep learning approaches (e.g., training a CNN on cell patches)
- Multi-scale features (e.g., neighborhood context)

These could be future extensions if basic morphology shows promise.

## References

- Current Module 3b implementation: `CITEgeist/model/module3b_nucleus_assignment.py`
- Morphology evaluation results: Job 8092351, 8094734
- Xenium data: `/ix1/alee/LO_LAB/General/Public_Data/10x_PublicData/Xenium_RNA_Proteomic_RenalCellCarcinoma/`
