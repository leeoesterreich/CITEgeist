# CITEgeist Module 4-5 Patient Demonstration Results

**Date**: 2026-02-05
**Status**: Complete

## Overview

Full CITEgeist pipeline execution on 14 HCC patient samples (6 patients, pre/post treatment timepoints) using a unified 10-cell-type profile derived from Module 1-2 discovery analysis.

## Pipeline Execution Summary

| Module | Description | Status | Key Output |
|--------|-------------|--------|------------|
| **1-2** | Profile Discovery | ✅ Complete | Unified 10-cell-type profile |
| **3** | Deconvolution | ✅ Complete | Cell proportions + GEX layers |
| **4** | Anchored + Joint Discovery | ✅ Complete | 590 programs (14 samples) |
| **4b** | Bivariate Relationships | ✅ Complete | 878 significant pairs |
| **5** | Cross-Sample Integration | ✅ Complete | 73 aligned programs, 191 relationships |

## Unified Cell Type Profile (10 types)

Derived from Module 1-2 marker colocalization analysis across all 14 samples:

| Cell Type | Markers |
|-----------|---------|
| Endothelial | PECAM1-1 |
| Fibroblasts | ACTA2-1 |
| B_Cells | CD19-1 |
| Macrophages | CD68-1, CD163-1 |
| Monocytes | CD14-1 |
| CD8_T_Cells | CD3E-1, CD8A-1 |
| CD4_T_Cells | CD4-1 |
| Cancer_Luminal | EPCAM-1 |
| Cancer_Basal | KRT5-1, SDC1-1 |
| Dendritic_Cells | ITGAX-1, HLA-DRA-1 |

## Sample Metadata

| Sample | Patient | Timepoint | Response |
|--------|---------|-----------|----------|
| HCC22-088-P1-S1 | P1 | Pre-treatment | Progressor |
| HCC22-088-P1-S2 | P1 | Post-treatment | Progressor |
| HCC22-088-P2-S1 | P2 | Pre-treatment | Responder |
| HCC22-088-P2-S2 | P2 | Post-treatment | Responder |
| HCC22-088-P3-S1_A | P3 | Pre-treatment | Responder |
| HCC22-088-P3-S2 | P3 | Post-treatment | Responder |
| HCC22-088-P4-S1 | P4 | Pre-treatment | Progressor |
| HCC22-088-P4-S2 | P4 | Post-treatment | Progressor |
| HCC22-088-P4-S2_1i_rep | P4 | Post-treatment (rep) | Progressor |
| HCC22-088-P5-S1 | P5 | Pre-treatment | Responder |
| HCC22-088-P5-S2 | P5 | Post-treatment | Responder |
| HCC22-088-P5-S2_F_rep | P5 | Post-treatment (rep) | Responder |
| HCC22-088-P6-S1 | P6 | Pre-treatment | Responder |
| HCC22-088-P6-S2_D | P6 | Post-treatment | Responder |

**Patient Distribution**: 4 Responders (P2, P3, P5, P6), 2 Progressors (P1, P4)

## Module 4 Results: Program Discovery

### Summary Statistics

| Metric | Value |
|--------|-------|
| Total samples | 14 |
| Total anchored programs | 590 |
| Total joint programs | 140 |
| Programs per sample | 35-50 |
| Cell types per sample | 3-10 |

### Per-Sample Breakdown

| Sample | Anchors | Anchored Programs | Joint Programs |
|--------|---------|-------------------|----------------|
| HCC22-088-P1-S1 | 8 | 40 | 10 |
| HCC22-088-P1-S2 | 9 | 45 | 10 |
| HCC22-088-P2-S1 | 9 | 45 | 10 |
| HCC22-088-P2-S2 | 10 | 50 | 10 |
| HCC22-088-P3-S1_A | 8 | 40 | 10 |
| HCC22-088-P3-S2 | 7 | 35 | 10 |
| HCC22-088-P4-S1 | 9 | 45 | 10 |
| HCC22-088-P4-S2 | 10 | 50 | 10 |
| HCC22-088-P4-S2_1i_rep | 10 | 50 | 10 |
| HCC22-088-P5-S1 | 8 | 40 | 10 |
| HCC22-088-P5-S2 | 10 | 50 | 10 |
| HCC22-088-P5-S2_F_rep | 8 | 40 | 10 |
| HCC22-088-P6-S1 | 3 | 15 | 10 |
| HCC22-088-P6-S2_D | 9 | 45 | 10 |

### Joint Program Type Distribution

Across all 14 samples:
- **Interaction programs**: 57% (cross-cell-type spatial relationships)
- **Microenvironment programs**: 26% (multi-cell-type niches)
- **Single cell type programs**: 16% (cell-type-specific signatures)

## Module 4b Results: Bivariate Relationships

### Summary Statistics

| Metric | Value |
|--------|-------|
| Total samples | 14 |
| Total pairs tested | ~1000 per sample |
| Total significant pairs | 878 |
| Samples with significant pairs | 9/14 |

### Per-Sample Results

| Sample | Programs | Pairs Tested | Significant | Co-localized | Exclusive |
|--------|----------|--------------|-------------|--------------|-----------|
| HCC22-088-P1-S1 | 40 | 780 | 82 | - | - |
| HCC22-088-P1-S2 | 45 | 990 | 113 | - | - |
| HCC22-088-P2-S1 | 45 | 990 | 151 | - | - |
| HCC22-088-P2-S2 | 50 | 1225 | 150 | - | - |
| HCC22-088-P3-S1_A | 40 | 780 | 108 | - | - |
| HCC22-088-P3-S2 | 35 | 595 | 79 | - | - |
| HCC22-088-P4-S1 | 45 | 990 | 0 | 0 | 0 |
| HCC22-088-P4-S2 | 50 | 1225 | 75 | - | - |
| HCC22-088-P4-S2_1i_rep | 50 | 1225 | 0 | 0 | 0 |
| HCC22-088-P5-S1 | 40 | 780 | 0 | 0 | 0 |
| HCC22-088-P5-S2 | 50 | 1225 | 0 | 0 | 0 |
| HCC22-088-P5-S2_F_rep | 40 | 780 | 0 | 0 | 0 |
| HCC22-088-P6-S1 | 15 | 105 | 22 | - | - |
| HCC22-088-P6-S2_D | 45 | 990 | 98 | - | - |

**Note**: Samples with 0 significant pairs (P4-S1, P4-S2_1i_rep, P5-S1, P5-S2, P5-S2_F_rep) passed stringent FDR correction but showed no significant spatial co-localization above threshold.

## Module 5 Results: Cross-Sample Integration

### Integration Statistics

| Metric | Value |
|--------|-------|
| Samples integrated | 14 |
| Input programs | 590 |
| Genes aligned | 2,824 |
| PCA components | 30 |
| PCA variance explained | 55.5% |
| Harmony converged | Yes |
| Harmony iterations | 8 |
| Aligned program groups | 73 |
| Highly conserved (>50%) | 5 |
| Conserved relationships | 191 |

### Conserved Relationship Types

| Type | Count | Description |
|------|-------|-------------|
| **Co-localized** | 26 | Programs that spatially cluster together |
| **Exclusive** | 6 | Programs that avoid each other spatially |
| **Independent** | 159 | No consistent spatial relationship |

### Top Conserved Programs

Programs present in the majority of patient samples, representing core biological processes:

| Program ID | Cell Type | Conservation | Samples | Top Genes |
|------------|-----------|--------------|---------|-----------|
| aligned_000 | Cancer_Luminal | 100% | 14/14 | ZNF550, C1orf159, NEDD4L, SPDEF |
| aligned_036 | CD8_T_Cells | 100% | 14/14 | SAMD11, MT-CYB, MT-ND6, MT-ND5 |
| aligned_009 | CD4_T_Cells | 86% | 12/14 | IGKC, IGHA1, JCHAIN, IGHM |
| aligned_060 | CD4_T_Cells | 71% | 10/14 | COL3A1, COL1A1, COL1A2, SPARC |
| aligned_004 | Dendritic_Cells | 64% | 9/14 | SPP1, CTSK, APOE, FTL |

### Top Conserved Co-localized Relationships

Program pairs that consistently appear together spatially across patients:

| Program 1 | Program 2 | Mean I | Samples | Interpretation |
|-----------|-----------|--------|---------|----------------|
| aligned_036 (CD8_T) | aligned_045 (CD8_T) | 0.152 | 4 | CD8 T cell subpopulation niches |
| aligned_004 (DC) | aligned_019 (CD4_T) | 0.196 | 4 | Dendritic-T cell interaction zones |
| aligned_004 (DC) | aligned_027 (DC) | 0.162 | 4 | Dendritic cell aggregation |
| aligned_000 (Cancer) | aligned_028 (Cancer_Basal) | 0.119 | 4 | Cancer phenotype co-occurrence |
| aligned_013 (Mac) | aligned_060 (CD4_T) | 0.101 | 4 | Macrophage-T cell niches |

### Top Conserved Exclusive Relationships

Program pairs that consistently avoid each other spatially:

| Program 1 | Program 2 | Mean I | Samples | Interpretation |
|-----------|-----------|--------|---------|----------------|
| aligned_000 (Cancer) | aligned_001 (CD4_T) | -0.145 | 4 | Tumor-immune exclusion |

### Response-Associated Programs

#### Responder-Enriched Programs (n=3)

Programs more frequent in treatment responders:

| Program ID | Cell Type | Responder | Progressor | Top Genes | Pathway |
|------------|-----------|-----------|------------|-----------|---------|
| aligned_021 | Macrophages | 67% | 20% | FABP4, HBA2, CFD | Lipid metabolism |
| aligned_013 | Macrophages | 44% | 20% | VIM, S100A6, C3, CD74 | Immune activation |
| aligned_025 | Cancer_Luminal | 44% | 20% | MGP, FOS, MT-CO3 | Stress response |

#### Progressor-Enriched Programs (n=4)

Programs more frequent in treatment non-responders:

| Program ID | Cell Type | Responder | Progressor | Top Genes | Pathway |
|------------|-----------|-----------|------------|-----------|---------|
| aligned_027 | Dendritic_Cells | 22% | 60% | FTL, FN1, TIMP1, APOE | ECM remodeling |
| aligned_028 | Cancer_Basal | 22% | 60% | MT-CO3, MT-ND4, MT-CO2 | Metabolic shift |
| aligned_005 | Cancer_Luminal | 11% | 40% | MGP, HBA2, LTF | Unknown |
| aligned_043 | Endothelial | 11% | 40% | KLHL17, MT-ND1 | Unknown |

## Biological Interpretation

### Key Findings

1. **Universal tumor signatures**: Cancer_Luminal (aligned_000) and CD8_T_Cells (aligned_036) programs are 100% conserved across all patients, representing core tumor biology independent of treatment response.

2. **Immune contexture in responders**: Responder-enriched macrophage programs show immune activation (C3, CD74) and lipid metabolism (FABP4, CFD), suggesting active anti-tumor immunity.

3. **ECM remodeling in progressors**: Progressor-enriched programs show ECM remodeling signatures (FN1, TIMP1) and metabolic adaptations (mitochondrial genes), potentially indicating immune evasion mechanisms.

4. **Stromal involvement**: The CD4_T_Cells programs showing collagen genes (COL1A1, COL3A1) likely represent stromal-associated T cell niches rather than pure T cell signatures.

5. **Spatial co-localization patterns**:
   - Dendritic cell-T cell interaction zones are conserved across patients
   - CD8 T cell subpopulations form spatial niches
   - Tumor-immune exclusion is observed (Cancer vs CD4_T programs are exclusive)

6. **Variable spatial structure**: Some samples (P4-S1, P5 samples) show no significant bivariate relationships, suggesting these tumors have more diffuse/mixed spatial organization.

### Treatment Response Biomarkers

Based on the integration analysis:

| Category | Biomarker Programs | Genes | Significance |
|----------|-------------------|-------|--------------|
| **Responder Signature** | Macrophage immune activation | C3, CD74, FABP4, CFD | Complement activation, lipid metabolism |
| **Progressor Signature** | DC ECM remodeling | FN1, TIMP1, APOE | Matrix remodeling, potential immune exclusion |
| **Progressor Signature** | Cancer metabolic shift | MT-CO3, MT-ND4 | Mitochondrial metabolism adaptation |

## Output Files

### Module 3 (Deconvolution)
- Location: `output/module3_unified/`
- Per-sample files:
  - `{sample}_cell_prop_global_results.csv`
  - `{sample}_cell_prop_finetuned_results.csv`
  - `{sample}_module3_results.h5ad`

### Module 4 (Program Discovery)
- Location: `output/module4_unified/`
- Per-sample files:
  - `{sample}_module4_discovery.json`
  - `{sample}_anchored_H.npy`
  - `{sample}_joint_W.npy`
  - `{sample}_joint_H.npy`

### Module 4b (Bivariate Relationships)
- Location: `output/module4b_unified/`
- Per-sample files:
  - `{sample}_module4b_relationships.csv`
  - `{sample}_module4b_summary.json`

### Module 5 (Integration)
- Location: `output/module5_integration/`
- Files:
  - `module5_unified_aligned_programs.csv` - 73 aligned programs
  - `module5_unified_conserved_relationships.csv` - 191 conserved relationships
  - `module5_unified_embedding.npy` - Harmony integration embedding
  - `module5_unified_program_metadata.csv` - Program metadata
  - `module5_unified_similarity_network.graphml` - Program network
  - `module5_response_analysis.json` - Responder vs Progressor patterns
  - `module5_summary.json` - Integration summary

## Technical Notes

### Computational Resources

| Module | Time per Sample | Memory | Total Time |
|--------|-----------------|--------|------------|
| Module 3 | ~30-60 min | 64 GB | ~10 hours |
| Module 4 | ~5-10 min | 64 GB | ~2 hours |
| Module 4b | ~5-30 min | 32 GB | ~3 hours |
| Module 5 | ~1 min | 64 GB | ~1 min |

### SLURM Scripts

- `examples/sbatch_module3_unified.sh` - Module 3 deconvolution
- `examples/sbatch_module4_unified.sh` - Module 4 program discovery
- `examples/sbatch_module4b_unified.sh` - Module 4b bivariate analysis
- `examples/sbatch_module5_unified.sh` - Module 5 integration

---

**Last Updated**: 2026-02-05 21:30 EST
**Pipeline Status**: Complete
