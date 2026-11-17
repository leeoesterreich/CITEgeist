# Seurat Simulation Benchmark

Seurat (SPOTlight-based) deconvolution on simulated spatial transcriptomics data.

## Directory Structure

```
Seurat/
├── high_seg/              # Highly segmented condition
│   ├── Seurat_hs.R        # Deconvolution script
│   ├── Seurat_hs.sh       # SLURM submission
│   ├── h5seurat/           # Converted input data
│   └── output/             # Results
└── mixed/                 # Mixed composition condition
    ├── Seurat_mx.R
    ├── Seurat_mx.sh
    ├── h5seurat/
    └── output/
```

## Running

```bash
sbatch high_seg/Seurat_hs.sh
sbatch mixed/Seurat_mx.sh
```

Results are evaluated by `../src/seurat_bench_wrapper.py`.
