# Single-Cell Resolution Demonstration Design

**Date:** 2026-02-06
**Purpose:** Demonstrate CITEgeist Modules 1, 2, and 4 at single-cell resolution on Xenium RCC data to future-proof the paper as spatial technology improves.

## Goals

1. **Validate against known biology** - Recover expected cell type profiles and marker relationships from RCC literature
2. **Discover novel findings** - Identify spatial programs and cell-cell relationships not obvious from standard clustering

## Data

- **Source:** Xenium RCC dataset (10x public data)
- **Scale:** 465,534 cells × 430 genes × 27 proteins
- **Spatial extent:** ~67 mm² tissue area

## Analysis Tracks

| Track | Data | Jobs | Purpose |
|-------|------|------|---------|
| **A (Full)** | 465K cells | 1 sbatch | Comprehensive results |
| **B (Quadrants)** | ~116K cells × 4 | 4 sbatch | Faster iteration, regional analysis |

Quadrants defined by X/Y midpoint:
- Q0: bottom-left
- Q1: bottom-right
- Q2: top-left
- Q3: top-right

## Pipeline Stages

### Stage 0: Data Preparation

**Script:** `run_singlecell_module12.py` (initial section)

**Processing:**
1. Load full single-cell data via existing `load_xenium_singlecell()`
2. Split into 4 spatial quadrants by X/Y midpoint
3. ~116K cells per quadrant

**Outputs:**
```
output_singlecell_demonstration/
├── data_summary.json          # Cell counts per quadrant
└── quadrants/
    ├── Q0/
    ├── Q1/
    ├── Q2/
    └── Q3/
```

### Stage 1-2: Module 1 & 2 (Marker Interest + Profile Discovery)

**Script:** `run_singlecell_module12.py`

**Module 1 (Marker Interest Detection):**
- Input: Protein expression matrix (n_cells × 27 proteins)
- Compute kurtosis, GMM SNR, Moran's I for each protein
- At single-cell resolution, Moran's I uses k=15 neighbors
- Output: Ranked interesting markers

**Module 2 (Profile Discovery):**
- Input: Interesting markers + spatial coordinates
- 2a: Pairwise colocalization (Jaccard, Pearson, bivariate Moran's I)
- 2b: Graph-based profile discovery via hierarchical clustering
- 2c: Select optimal profiles by reconstruction accuracy
- Output: Discovered cell type profiles

**Outputs per run:**
```
{output_dir}/
├── module1_marker_interest.json
├── module1_marker_interest.csv
├── module2a_colocalization.csv
├── module2b_profiles_raw.json
├── module2c_profiles_selected.json
└── module2_dendrogram.png
```

**Expected profiles based on literature:**
- T cells: CD3E + CD4 or CD8A
- Macrophages: CD68 + CD163
- B cells: CD20
- Epithelial: PanCK ± E-cadherin
- Stromal: Vimentin + αSMA
- Endothelial: CD31

### Stage 3: Module 4 (Program Discovery)

**Script:** `run_singlecell_module4.py`

**Key difference from spot-level:** No deconvolution needed. At single-cell resolution, assign each cell to a profile, then run NMF on RNA within each cell type.

**Processing:**
1. Load discovered profiles from Module 2c output
2. Assign cells to profiles via protein expression similarity
3. For each cell type anchor:
   - Extract cells assigned to that type
   - Run NMF on their gene expression (K=5 programs)
   - Compute Moran's I for spatial coherence
   - Identify top genes per program
4. Run Module 4b: Bivariate relationships between programs

**Outputs:**
```
{output_dir}/
├── cell_assignments.csv
├── module4_programs/
│   ├── {cell_type}_programs.json
│   ├── {cell_type}_gene_loadings.csv
│   └── {cell_type}_cell_scores.csv
├── module4_summary.json
└── module4b_relationships.json
```

**Parameters:**
- `n_programs`: 5 (default)
- `assignment_method`: "soft" or "hard"
- `min_cells_per_type`: 500

### Stage 4: Evaluation & Validation

**Script:** `evaluate_singlecell.py`

**Part A: Profile Validation (Modules 1-2)**

Compare discovered profiles against expected biology:
```python
EXPECTED_PROFILES = {
    "T_helper": {"required": ["CD3E", "CD4"], "optional": ["CD45RO"]},
    "T_cytotoxic": {"required": ["CD3E", "CD8A"], "optional": ["GranzymeB"]},
    "Macrophage_M2": {"required": ["CD68", "CD163"], "optional": ["HLA-DR"]},
    "B_cell": {"required": ["CD20"], "optional": ["CD45RA"]},
    "Epithelial": {"required": ["PanCK"], "optional": ["E-Cadherin", "Beta-catenin"]},
    "Stromal_CAF": {"required": ["Vimentin", "alphaSMA"]},
    "Endothelial": {"required": ["CD31"], "optional": ["Vimentin"]},
}
```

Metrics: precision/recall/F1 for marker assignments

**Part B: Program Validation (Module 4)**

Gene set enrichment against:
- MSigDB Hallmark gene sets
- KEGG pathways
- Published ccRCC signatures (hypoxia, angiogenesis, immune infiltration)

**Outputs:**
```
{output_dir}/evaluation/
├── profile_validation.json
├── program_gsea_results.csv
├── validation_summary.md
└── discovery_catalog.md
```

### Stage 5: Figures & Final Outputs

**Script:** `generate_figures.py`

**Figures:**
1. `fig_profiles.png` - Profile discovery overview (heatmap + comparison)
2. `fig_spatial_celltypes.png` - Full tissue with cells colored by type
3. `fig_programs.png` - Top genes per program, Moran's I coherence
4. `fig_spatial_programs.png` - Selected programs plotted spatially
5. `fig_validation.png` - Profile recovery matrix + GSEA dot plot

**Final Deliverables:**
```
output_singlecell_demonstration/
├── figures/
│   ├── fig_profiles.png
│   ├── fig_spatial_celltypes.png
│   ├── fig_programs.png
│   ├── fig_spatial_programs.png
│   └── fig_validation.png
├── validation_report.json
├── discovery_catalog.md
└── paper_supplement.md
```

## Known Biology to Validate

### Cell Type Composition (from literature)
- **T cells** (~40% of immune cells): CD8+ (14.8%), CD4+, with exhausted subsets
- **Macrophages** (~30%): M2-polarized (CD163+, CD68+) predominate
- **B cells** (~3%): CD20+ cells, often in tertiary lymphoid structures
- **Epithelial/Tumor** (~26%): PanCK+, E-cadherin often lost in ccRCC
- **Stromal**: αSMA+ CAFs, Vimentin+ mesenchymal cells
- **Endothelial**: CD31+ vasculature

### Expected Spatial Patterns
- **Tertiary Lymphoid Structures (TLS)**: CD20+ B cells in center, surrounded by CD3+ T cells
- **Immune phenotypes**: Infiltrated vs excluded vs desert regions
- **M2 macrophages near tumor**: Immunosuppressive niche
- **Exhausted CD8+ T cells**: Enriched in metabolically stressed tumor areas
- **CAFs at tumor-stroma interface**: αSMA+ fibroblasts

### Marker Co-expression Expectations
- CD3E + CD4 or CD8A (T cell subsets)
- CD68 + CD163 (M2 macrophages)
- PanCK + E-cadherin (epithelial, though E-cad often lost in ccRCC)
- Vimentin + αSMA (activated fibroblasts)

## Literature Sources

- [Metabolic heterogeneity in ccRCC - J Transl Med 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC10900752/)
- [Spatial transcriptomics in metastatic RCC 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC12496100/)
- [TLS impact on immunotherapy in ccRCC](https://pmc.ncbi.nlm.nih.gov/articles/PMC10784869/)
- [Spatial heterogeneity in ccRCC prognosis](https://pmc.ncbi.nlm.nih.gov/articles/PMC10360235/)
- [Comprehensive immunoprofiles of RCC subtypes](https://pmc.ncbi.nlm.nih.gov/articles/PMC7139472/)

## File Structure

```
Benchmarking/xenium_benchmarking/CITEgeist/
├── src/
│   ├── run_singlecell_module12.py     # NEW: Module 1-2 at single-cell
│   ├── run_singlecell_module4.py      # NEW: Module 4 using discovered profiles
│   ├── evaluate_singlecell.py         # NEW: Validation + GSEA
│   └── generate_figures.py            # NEW: Publication figures
├── slurm/
│   ├── run_singlecell_full.sh         # NEW: Full 465K cell run
│   ├── run_singlecell_quadrants.sh    # NEW: Per-quadrant runs (array job)
│   └── run_evaluation.sh              # NEW: Evaluation job
└── output_singlecell_demonstration/   # NEW: All outputs
    ├── full/                          # Full dataset results
    ├── quadrants/                     # Per-quadrant results
    │   ├── Q0/
    │   ├── Q1/
    │   ├── Q2/
    │   └── Q3/
    ├── evaluation/                    # Validation reports
    └── figures/                       # Publication figures
```

## Reused Components

- `load_xenium_singlecell.py` - Data loading (exists)
- `CITEgeist/model/marker_interest.py` - Module 1 (exists)
- `CITEgeist/model/spatial_colocalization.py` - Module 2 (exists)
- `CITEgeist/model/anchored_program_discovery.py` - Module 4 (exists)

## Implementation Notes

- Module 1-2 parameters may need tuning for single-cell resolution (k neighbors for Moran's I)
- Cell assignment in Module 4 replaces deconvolution (single-cell = no mixing)
- GSEA requires MSigDB gene sets (download or use gseapy)
- Figures should be publication-quality (300 DPI, proper fonts)
