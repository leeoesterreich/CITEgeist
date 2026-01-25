# MDK Saturation Pipeline

Unified analysis pipeline proving the ER saturation model for MDK secretion in D538G breast cancer cells.

## Overview

This pipeline integrates 7 datasets to explain why ESR1-D538G causes **opposite** effects on MDK secretion:
- **MCF7-D538G:** MDK secretion UP
- **T47D-D538G:** MDK secretion DOWN

## Quick Start

```bash
# Run full pipeline
python run_pipeline.py

# Run specific script
python scripts/02_analyze_chaperone_expression.py

# Start from specific step
python run_pipeline.py --from 3
```

## Pipeline Structure

| Script | Question | Datasets |
|--------|----------|----------|
| 01 | What did CITEgeist find? | Vignette 4 |
| 02 | Opposite chaperone regulation? | GSE89888 |
| 03 | Opposite ER binding? | GSE125117 |
| 04 | Why different responses? | GSE254216, GSE72249 |
| 05 | Does perturbation confirm? | GSE254218, GSE75329 |
| 06 | Do all datasets agree? | Outputs 01-05 |
| 07 | Generate report | All outputs |

## Datasets

See `data/manifest.yaml` for full documentation.

## Outputs

- `outputs/tables/` - CSV files from each analysis
- `outputs/figures/` - Publication-ready figures (PDF + PNG)
- `outputs/report.md` - Auto-generated summary report

## Requirements

- Python 3.10
- pandas, numpy, scipy, statsmodels
- matplotlib, seaborn
- pyyaml
- mygene (optional, for gene symbol mapping)

## Configuration

Edit `config.yaml` to modify paths or parameters.
