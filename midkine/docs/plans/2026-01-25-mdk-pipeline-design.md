# MDK Saturation Pipeline Design

**Date:** 2026-01-25
**Status:** Approved
**Purpose:** Consolidate MDK/D538G analysis into a unified, reproducible pipeline

## Overview

This pipeline proves the ER saturation model explaining why ESR1-D538G causes opposite MDK secretion effects in MCF7 vs T47D cells. It integrates 7 datasets (6 GEO + 1 spatial) into a narrative that builds from observation to mechanism to validation.

## Directory Structure

```
midkine/mdk_saturation_pipeline/
├── scripts/
│   ├── 00_download_data.py                  # Optional: fetch GEO files
│   ├── 01_summarize_spatial_finding.py      # Vignette 4
│   ├── 02_analyze_chaperone_expression.py   # GSE89888
│   ├── 03_analyze_er_binding_changes.py     # GSE125117
│   ├── 04_quantify_saturation.py            # GSE254216 + GSE72249
│   ├── 05_analyze_foxa1_perturbations.py    # GSE254218 + GSE75329
│   ├── 06_cross_validate.py                 # Integrates 01-05
│   └── 07_generate_report.py                # Final summary
├── data/
│   ├── manifest.yaml                        # Dataset documentation
│   └── geo/                                 # GEO files (gitignored)
├── outputs/
│   ├── figures/                             # Publication-ready PNGs/PDFs
│   ├── tables/                              # CSVs for each analysis
│   └── report.md                            # Auto-generated summary
├── run_pipeline.py                          # Runs scripts 01-07 in order
├── config.yaml                              # Paths & parameters
└── README.md                                # Setup & usage instructions
```

## Datasets

| Dataset | Type | Description | Used In |
|---------|------|-------------|---------|
| Vignette 4 | Spatial | CITEgeist MDK program analysis | 01 |
| GSE89888 | RNA-seq | D538G vs WT expression (MCF7 + T47D) | 02 |
| GSE125117 | ER ChIP-seq | D538G vs WT binding changes | 03 |
| GSE254216 | ATAC-seq | Chromatin accessibility (saturation) | 04 |
| GSE72249 | FOXA1/ER ChIP-seq | FOXA1 binding at chaperone loci | 04 |
| GSE254218 | RNA-seq | FOXA1 knockdown in T47D | 05 |
| GSE75329 | RNA-seq | FOXA1/ER knockdown/OE in MCF7 | 05 |

## Script Responsibilities

### 01_summarize_spatial_finding.py
- **Question:** What did CITEgeist find about MDK?
- **Datasets:** Vignette 4 outputs
- **Outputs:** `spatial_summary.csv`, observation text

### 02_analyze_chaperone_expression.py
- **Question:** Do chaperones show opposite regulation in MCF7 vs T47D?
- **Datasets:** GSE89888
- **Outputs:** `chaperone_expression.csv`, 2-way ANOVA stats, heatmap

### 03_analyze_er_binding_changes.py
- **Question:** Does ER binding change oppositely at chaperone loci?
- **Datasets:** GSE125117
- **Outputs:** `binding_changes.csv`, peak counts at HSP90B1/HSPA5/CALR

### 04_quantify_saturation.py
- **Question:** Why do MCF7 and T47D respond differently?
- **Datasets:** GSE254216 (ATAC) + GSE72249 (FOXA1 ChIP)
- **Outputs:** `saturation_metrics.csv`, ER/ATAC occupancy, FOXA1 signal

### 05_analyze_foxa1_perturbations.py
- **Question:** Does perturbing FOXA1/ER confirm the model?
- **Datasets:** GSE254218 + GSE75329
- **Outputs:** `perturbation_results.csv`, KD/OE effects on chaperones

### 06_cross_validate.py
- **Question:** Do all lines of evidence agree?
- **Datasets:** Outputs from 01-05
- **Outputs:** `cross_validation.csv`, concordance matrix

### 07_generate_report.py
- **Question:** Compile everything into narrative
- **Datasets:** All outputs
- **Outputs:** `report.md`, assembled figures

## Cross-Validation Logic

Script 06 performs concordance checks across independent datasets:

| Check | Datasets Compared | Expected Agreement |
|-------|-------------------|-------------------|
| Expression ↔ Binding | GSE89888 vs GSE125117 | Chaperones UP in MCF7 = ER binding LOST |
| Binding ↔ Saturation | GSE125117 vs GSE254216 | MCF7 loses peaks because saturated (3%) |
| Saturation ↔ FOXA1 | GSE254216 vs GSE72249 | T47D has more FOXA1 at chaperones |
| Perturbation ↔ Model | GSE254218/GSE75329 vs all | FOXA1-KD derepresses chaperones |

Output format:
```
check,dataset_a,dataset_b,metric,expected,observed,concordant
expr_binding,GSE89888,GSE125117,HSP90B1_direction,opposite,opposite,TRUE
```

## Configuration

**config.yaml:**
```yaml
# Paths
data_dir: "data/geo"
spatial_dir: "../CITEgeist/examples/output_vignette4_mdk"
output_dir: "outputs"

# Dataset files
datasets:
  GSE89888: "GSE89888_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz"
  GSE125117_MCF7_WT: "GSM3563751_MCF7_WT_E2_peaks.bed.gz"
  GSE125117_MCF7_D538G: "GSM3563757_MCF7_D538G_E2_peaks.bed.gz"
  GSE125117_T47D_WT: "GSM3563760_T47D_WT_E2_peaks.bed.gz"
  GSE125117_T47D_D538G: "GSM3563766_T47D_D538G_E2_peaks.bed.gz"
  GSE254216_MCF7_ATAC: "GSM8036144_MCF7_summits.bed.gz"
  GSE254216_T47D_ATAC: "GSM8036152_T47D_summits.bed.gz"
  GSE254218: "GSE254218_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz"
  GSE75329: "GSE75329_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz"
  GSE72249_MCF7_FOXA1: "GSM1858624_GH917_FoxA1_unt_MCF7_ChIP_hg19.bedgraph.gz"
  GSE72249_T47D_FOXA1: "GSM1858654_GH1070_T47D_FOXA1_Unt_ChIP_hg19.bedgraph.gz"

# Parameters
significance_threshold: 0.05
fold_change_threshold: 1.2

# Flags
rerun_spatial: false
```

## Output Report Structure

**outputs/report.md:**
```markdown
# MDK Secretion Mechanism: Evidence Report
Generated: [timestamp]

## 1. Observation
[From 01: MCF7-D538G MDK UP, T47D-D538G MDK DOWN]

## 2. Chaperone Expression (GSE89888)
[Table + heatmap]

## 3. ER Binding Changes (GSE125117)
[Binding delta table]

## 4. Saturation Model (GSE254216 + GSE72249)
[Occupancy metrics]

## 5. Causal Validation (GSE254218 + GSE75329)
[Perturbation results]

## 6. Cross-Validation Summary
[X/Y checks passed]

## 7. Conclusion
[Auto-generated]
```

## Figure Outputs

| Figure | Content | Format |
|--------|---------|--------|
| fig1_spatial_observation.pdf | MDK programs from CITEgeist | PDF + PNG |
| fig2_chaperone_heatmap.pdf | Expression heatmap | PDF + PNG |
| fig3_binding_changes.pdf | ER peak counts at chaperone loci | PDF + PNG |
| fig4_saturation_model.pdf | ATAC vs ER occupancy | PDF + PNG |
| fig5_perturbation_validation.pdf | FOXA1 KD/OE effects | PDF + PNG |

## Error Handling

Each script validates inputs and fails early:
```python
def validate_inputs():
    required = [config['datasets']['GSE89888'], ...]
    missing = [f for f in required if not Path(f).exists()]
    if missing:
        sys.exit(f"Missing required files:\n" + "\n".join(missing))
```

Pipeline runner provides clear failure messages:
```
❌ Pipeline failed at 03_analyze_er_binding_changes.py
   Missing: data/geo/GSM3563751_MCF7_WT_E2_peaks.bed.gz
   Run: python scripts/00_download_data.py to fetch GEO files
```

## Usage

```bash
# Full pipeline
python run_pipeline.py

# Individual script
python scripts/02_analyze_chaperone_expression.py

# With spatial rerun
python run_pipeline.py --rerun-spatial
```

## Implementation Notes

- Each script runs independently for testing/debugging
- Scripts read from `data/` or previous outputs in `outputs/tables/`
- Scripts write to `outputs/tables/` and `outputs/figures/`
- Return exit code 0 on success for pipeline orchestration
- Existing analysis scripts in `midkine/` can be refactored into this structure
