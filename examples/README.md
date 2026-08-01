# CITEgeist Example Scripts

End-to-end runners for the CITEgeist spatial multi-omics pipeline on Visium + CITE-seq
patient samples. Each script is a thin command-line wrapper around the `CITEgeist` package
modules; run them in the order below to reproduce the per-sample manuscript analysis.

> **Paths are hardcoded.** Input-data and per-module output locations are module-level
> constants near the top of every script (`DATA_ROOT` / `PATIENT_DATA_ROOT`, `MODULE3_OUTPUT`,
> `MODULE4_OUTPUT`, `MODULE4B_OUTPUT`). Edit them to point at your local copy before running.
> Every script ships with a `/path/to/CITEgeist_public_data` placeholder.

## Data

All 12 patient samples (`sample-P1-S1` … `sample-P6-S2`) are standard 10x Visium SpaceRanger
output directories, each laid out as `<DATA_ROOT>/<sample>/outs/` containing
`filtered_feature_bc_matrix.h5`, `spatial/`, and (for single-cell assignment) the H&E
`tissue_fullres_image.tif`. CITE-seq antibody capture is in the same SpaceRanger matrix
(loaded with `gex_only=False`); antibody names carry the `-1` suffix (e.g. `CD68-1`).

The processed data are publicly available at NCBI GEO under accession
[GSE289326](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE289326). After downloading, set `DATA_ROOT` /
`PATIENT_DATA_ROOT` in each script to your local
`.../CITEgeist_public_data/processed_files` directory.

## Requirements

- A Python environment with CITEgeist installed (`pip install citegeist`, or from source —
  `git clone` the repo, then `pip install -e .`; the repo `CITEgeist_env.yml` reproduces the
  full analysis environment).
- **A GPU node with the NVIDIA cuOPT library** for the Module 3 QP solver
  (`run_module3_unified.py`, `run_cuopt_qp_patient.py`). cuOPT has no CPU fallback.
- StarDist (H&E nuclei segmentation) for `run_single_cell_assignment.py` Stage 1
  (GPU recommended).

## Scripts and run order

| # | Script | Module(s) | Consumes | Produces | GPU |
|---|--------|-----------|----------|----------|-----|
| 1 | `run_module12_discovery.py` | M1 + M2 | raw CITE-seq antibody counts (`DATA_ROOT`) | discovered cell-type profiles (JSON) | no |
| 2 | `run_module3_unified.py` | M3 (QP proportions + SACE GEX) | `DATA_ROOT` + M1–2 profile | proportions + deconvolved GEX layers | **yes** |
| 2′ | `run_cuopt_qp_patient.py` | M3 (QP proportions only) | `DATA_ROOT` (single sample) | proportions (`output/module3_cuopt_qp`) | **yes** |
| 3 | `run_single_cell_assignment.py` | M3-post | M3 proportions + StarDist nuclei | per-cell type assignments + single-cell AnnData | Stage 1 (StarDist) |
| 4 | `run_module4_discovery.py` | M4 | `MODULE3_OUTPUT` deconvolved GEX | anchored NMF gene programs | no |
| 5 | `run_module4b_bivariate.py` | M4b | `MODULE3_OUTPUT` + `MODULE4_OUTPUT` | bivariate Moran's I between program pairs | no |
| 6 | `run_module5_integration.py` | M5 | `MODULE4_OUTPUT` + `MODULE4B_OUTPUT` | conserved cross-sample programs | no |

`run_cuopt_qp_patient.py` (2′) is a standalone single-patient alternative to the QP half of
`run_module3_unified.py`; use one or the other to produce the Module 3 proportions that
steps 3–5 consume.

## Example invocations

```bash
# 1. Marker interest + profile discovery
python run_module12_discovery.py --sample sample-P1-S1 --output-dir output/module12_discovery

# 2. Deconvolution (proportions + GEX) — GPU node
python run_module3_unified.py --sample sample-P1-S1 --output-dir output/module3_unified

# 3. Per-cell assignment (StarDist -> per-nucleus type assignment -> SACE per-cell GEX)
python run_single_cell_assignment.py --sample sample-P1-S1 --stages 1,5,6

# 4. Anchored gene-program discovery
python run_module4_discovery.py --sample sample-P1-S1 --mode both --output-dir output/module4_discovery

# 5. Bivariate program relationships
python run_module4b_bivariate.py --sample sample-P1-S1 --output-dir output/module4b_unified

# 6. Cross-sample integration (all 12 samples)
python run_module5_integration.py --output-dir output/module5_integration
```

See each script's module-level docstring for the full argument list and per-sample metadata.
