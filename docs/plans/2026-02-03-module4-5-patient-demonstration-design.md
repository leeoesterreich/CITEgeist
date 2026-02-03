# Module 4-5 Patient Sample Demonstration Design

**Date:** 2026-02-03
**Status:** Draft
**Goal:** Demonstrate utility of Modules 4 and 5 using 6 paired patient samples (biopsy vs post-endocrine-therapy surgical specimens)

---

## Overview

Use the HCC22-088 patient cohort to demonstrate the full CITEgeist pipeline (Modules 1-5), with focus on:
1. **Treatment response** - Programs that change between biopsy (S1) and surgical (S2)
2. **Conserved biology** - Programs consistent across all 6 patients
3. **Patient heterogeneity** - Differential trajectories between responders and progressors

## Patient Cohort

| Patient | Biopsy (S1) | Surgical (S2) | Response Status |
|---------|-------------|---------------|-----------------|
| P1 | HCC22-088-P1-S1 | HCC22-088-P1-S2 | **Progressor** |
| P2 | HCC22-088-P2-S1 | HCC22-088-P2-S2 | Responder |
| P3 | HCC22-088-P3-S1_A | HCC22-088-P3-S2 | Responder |
| P4 | HCC22-088-P4-S1 | HCC22-088-P4-S2, P4-S2_1i_rep | **Progressor** |
| P5 | HCC22-088-P5-S1 | HCC22-088-P5-S2, P5-S2_F_rep | Responder |
| P6 | HCC22-088-P6-S1 | HCC22-088-P6-S2_D | Responder |

**Treatment:** Neoadjuvant endocrine therapy (ER+ breast cancer)
**Existing outputs:** Module 3 complete (cell proportions + gene expression deconvolution)
**Data location:** `/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeistNeilAnalysis/CITEgeist/output/`

---

## Analysis Architecture

```
Module 1-2 (per sample)     Profile Curation     Module 3 (existing)
─────────────────────       ────────────────     ───────────────────
Auto-discover profiles  →   Manual refinement →  Cell proportions +
from all 14 samples         with expert input    Deconvolved GEX
                                    ↓
                            Hybrid Profile
                                    ↓
        ┌───────────────────────────┴───────────────────────────┐
        ↓                                                       ↓
Module 4 Anchored (existing)                    Module 4 Joint (new)
────────────────────────────                    ────────────────────
Per cell-type NMF                               All cell types together
Cell-type-specific programs                     Cross-cell-type programs
Module 4b: cross-anchor relationships           Post-hoc cell type assignment
        │                                                       │
        └───────────────────────────┬───────────────────────────┘
                                    ↓
                        Module 5 Integration
                        ────────────────────
                        Harmony-style alignment
                        Meta-program catalog
                        Conservation scoring
                                    ↓
                          Analysis Layer
                          ──────────────
                        Treatment response (S1→S2)
                        Differential trajectories
                        Responder vs Progressor
```

---

## Profile Discovery Workflow (Modules 1-2)

### Step 1: Run Modules 1-2 on all samples independently

Per sample:
- Module 1: Identify interesting markers (kurtosis, GMM, Moran's I)
- Module 2a: Compute marker colocalization
- Module 2b: Discover profiles via clustering
- Module 2c: Select optimal profile count

### Step 2: Compare discovered vs original profiles

| Original Profile | Discovered Profile(s) | Action |
|------------------|----------------------|--------|
| Cancer (EPCAM+) | Matches | Keep |
| CD163+ Macrophage | Splits into 2 clusters? | Evaluate biological meaning |
| Missing | New marker combination | Add if validated |

### Step 3: Manual curation session

Present:
- Side-by-side comparison of original vs discovered markers per cell type
- Marker colocalization heatmaps
- Recommendations based on consistency across samples

### Step 4: Finalized hybrid profile

Curated profile becomes input for downstream analysis.

---

## Module 4: Hybrid Discovery Strategy

### Mode A: Anchored Discovery (existing implementation)

```python
discover_anchored_programs(adata, cell_type_proportions, ...)
```

- Runs NMF per cell type anchor independently
- Contrastive NMF finds anchor-specific programs
- Module 4b computes cross-anchor bivariate Moran's I

**Strengths:** Cell-type-specific programs, clear biological assignment
**Captures:** "Macrophage inflammatory program", "Cancer proliferation program"

### Mode B: Joint Discovery (new implementation needed)

```python
discover_joint_programs(adata, cell_type_proportions, ...)
```

- Stack all deconvolved layers using `stack_deconvolved_layers()`
- Single NMF on combined matrix
- Post-hoc cell type assignment based on proportion correlations

**Strengths:** Programs can naturally span cell types
**Captures:** "Tumor-immune interface program", "Hypoxic niche program"

### Program Cell Type Assignment (for joint mode)

For each discovered program, compute:

1. **Cell type enrichment scores** - correlate program activity (H matrix) with cell type proportions across spots
2. **Gene source attribution** - fraction of top-loaded genes that are canonical markers for each cell type

| Label Type | Example | Criteria |
|------------|---------|----------|
| Single cell type | "Macrophage_inflammatory" | >70% enrichment in one cell type |
| Interaction | "Cancer-Tcell_interface" | Two cell types each >25% |
| Microenvironment | "TME_hypoxia" | Broad signal, no dominant cell type |

### Comparison Analysis

| Question | How to answer |
|----------|---------------|
| Do joint programs recapitulate anchored ones? | Correlate W matrices |
| What does joint find that anchored misses? | Programs with balanced multi-cell-type loadings |
| What does anchored+4b find that joint misses? | Spatially adjacent but expression-distinct programs |
| Which predicts treatment response better? | Compare responder/progressor discrimination |

---

## Module 5: Cross-Sample Integration

Uses existing Harmony-style implementation in `cross_sample_integration.py`.

### Integration approach

1. **Build program similarity matrix** - cosine similarity of W matrices (gene loadings)
2. **Harmony batch correction** - align programs across samples
3. **Cluster into meta-programs** - hierarchical clustering on similarity
4. **Conservation scoring:**

| Conservation Level | Definition | Interpretation |
|--------------------|------------|----------------|
| Universal (6/6 patients) | All patients, both timepoints | Core biology |
| Treatment-associated | Mainly S1 or mainly S2 | Treatment response |
| Response-specific | Mainly responders or progressors | Resistance/sensitivity |
| Patient-specific | 1-2 patients only | Individual heterogeneity |

---

## Analysis Layer

### Analysis 1: Treatment Effect (S1 → S2)

For each meta-program:
- Compute mean activity in S1 (biopsy) vs S2 (surgical)
- Paired test across 6 patients
- Classify as: **induced**, **suppressed**, or **stable**

Look for all types of changes:
- **Activity shifts** - same programs, different intensity
- **Program emergence/disappearance** - programs only in S1 or S2
- **Composition shifts** - same concept, different gene loadings

### Analysis 2: Conserved Biology

Programs found in ≥5/6 patients regardless of timepoint:
- Represent core tissue biology
- High conservation validates method finding real signal
- Report top genes and pathway enrichment

### Analysis 3: Differential Trajectories

Focus: How programs change differently in responders vs progressors.

For each meta-program, compute trajectory = (S2 activity - S1 activity):

| Program | Responder trajectory (n=4) | Progressor trajectory (n=2) | Interpretation |
|---------|---------------------------|----------------------------|----------------|
| Proliferation | -0.35 ± 0.08 | -0.05 ± 0.12 | Progressors maintain proliferation |
| Immune activation | +0.25 ± 0.10 | +0.05 ± 0.08 | Responders activate immune |

Statistical test: Compare trajectory distributions between groups.

---

## Modular Component Library

**Philosophy:** Lego blocks, not a fixed pipeline. Each component is standalone and can be rearranged as needed.

| Component | Input | Output | Reusable for |
|-----------|-------|--------|--------------|
| `run_module12_discovery.py` | Raw antibody data | Discovered profiles JSON | Any sample |
| `compare_profiles.py` | Discovered + curated profiles | Comparison table, heatmap | Profile refinement |
| `run_module4_anchored.py` | Module 3 output | Per-anchor programs | Any sample |
| `run_module4_joint.py` | Module 3 output | Joint programs | Any sample |
| `run_module5_integration.py` | Multi-sample Module 4 | Aligned meta-programs | Any sample set |
| `analyze_treatment_response.py` | Module 5 + sample metadata | S1 vs S2 comparisons | Paired designs |
| `analyze_trajectories.py` | Module 5 + response labels | Trajectory divergence | Case/control designs |
| `plot_program_spatial.py` | Module 4 + coords | Spatial activity maps | Any visualization |
| `plot_program_heatmap.py` | Module 5 results | Program × sample heatmap | Any summary |

### Design Principles

1. **Flat file outputs** - CSVs and PNGs that can be loaded anywhere
2. **Minimal dependencies between components** - each reads files, not in-memory objects
3. **Consistent naming conventions** - `{sample}_{module}_{output_type}.{ext}`
4. **Metadata sidecar files** - JSON with parameters used for reproducibility

---

## Figure Components (Lego Blocks)

### Available for Main Figures
- Method overview schematic
- Anchored vs Joint program comparison (Venn, correlations)
- Spatial maps of program activity
- Treatment response heatmap (programs × samples)
- Volcano plot (S1→S2 change vs significance)
- Trajectory divergence plot (responder vs progressor arrows)

### Available for Supplementary
- Module 1-2 discovered vs curated profiles
- Conservation score distributions
- Per-patient program activity
- Pathway enrichment tables

### Tutorial/Vignette Components
- Data loading example
- Module 1-2 + curation workflow
- Module 4 both modes
- Module 5 integration
- Analysis examples

---

## Implementation Order

### Phase 1: Foundation
- [ ] Implement `discover_joint_programs()` in anchored_program_discovery.py
- [ ] Add cell type assignment logic for joint programs
- [ ] Test on one sample to validate

### Phase 2: Profile Discovery
- [ ] Run Modules 1-2 on all 14 samples
- [ ] Generate comparison with original curated profiles
- [ ] Manual curation session → finalized hybrid profile

### Phase 3: Module 4 Execution
- [ ] Run anchored discovery on all samples with new profile
- [ ] Run joint discovery on all samples
- [ ] Compare outputs between approaches

### Phase 4: Module 5 Integration
- [ ] Integrate anchored programs across samples
- [ ] Integrate joint programs across samples
- [ ] Build meta-program catalog

### Phase 5: Analysis Components
- [ ] Treatment response analysis (S1 vs S2)
- [ ] Trajectory analysis (responders vs progressors)
- [ ] Conservation scoring

### Phase 6: Figures & Polish
- [ ] Generate figure components
- [ ] Assemble tutorial notebook
- [ ] Document component library

---

## Dependencies

### New Implementation Required
- `discover_joint_programs()` function
- Cell type assignment for joint programs
- Comparison utilities (anchored vs joint)

### Existing Infrastructure (ready to use)
- `stack_deconvolved_layers()` - stacks cell-type layers
- `unstack_program_results()` - splits results back
- `integrate_samples()` - Harmony-style integration (Module 5)
- Module 3 outputs already computed for all samples

---

## Notes

- Keep components modular for flexibility in final paper architecture
- Other elements floating in project may be incorporated later
- Focus on reusable blocks rather than rigid pipeline
