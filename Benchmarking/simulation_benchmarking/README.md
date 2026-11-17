# Simulation Benchmarking

Multi-method cell-type deconvolution benchmark on simulated spatial transcriptomics data (Fig 2 in manuscript).

## Conditions

Two simulation conditions from scCube (see `../simulation_generation/`):
- **high_seg**: Highly segmented -- distinct spatial regions with relatively homogeneous cell type compositions
- **mixed**: Mixed composition -- heterogeneous cell type mixing across tissue

3 replicates per condition.

## Methods

| Directory | Method | Language |
|-----------|--------|----------|
| `CARD/` | CARD | R |
| `Cell2Location/` | Cell2Location | Python |
| `CITEgeist/` | CITEgeist | Python |
| `RCTD/` | RCTD | R |
| `Seurat/` | Seurat (SPOTlight) | R |
| `Tangram/` | Tangram | Python |
| `scResolve/` | scResolve | Python |

## Shared Code

- `src/` -- Benchmarking metrics and method wrappers (see `src/README.md`)
- `evaluation/` -- Cross-method evaluation and comparison
- `Figures/` -- R visualization scripts for manuscript figures
- `Runtime/` -- Runtime comparison across methods

## Running

Each method directory contains per-condition scripts (e.g., `high_seg/`, `mixed/`). Submit via SLURM:

```bash
sbatch <Method>/high_seg/<script>.sh
sbatch <Method>/mixed/<script>.sh
```

After all methods complete, evaluate:
```bash
python src/compare_prop_methods.py
python src/compare_gex_methods.py
```

## Archived

`_archive/` -- deprecated method variants and exploratory benchmarks.
