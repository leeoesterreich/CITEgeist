# Xenium Benchmarking

Multi-method benchmarking of spatial transcriptomics deconvolution on Xenium breast cancer pseudo-Visium data with SingleR ground truth.

## Ground Truth

SingleR label transfer from GSE156632 scRNA-seq reference (generated 2026-03-17):

- **6-type**: B cells, T cells (merged), Macrophages, Endothelial, Epithelial, Fibroblasts
- **7-type**: B cells, CD4+ T cells, CD8+ T cells, Macrophages, Endothelial, Epithelial, Fibroblasts

5 tissue regions, 7,054 total spots (55 um pseudo-Visium diameter).

Ground truth files: `ground_truth_singler/ground_truth_{6,7}type/ground_truth/Xenium_region_{0-4}_prop.csv`

## Results Summary

See `evaluation/results/canonical_results_all.json` for full per-region results.

| Method | 6-type r | 7-type r |
|--------|----------|----------|
| CITEgeist NB | 0.764 | 0.726 |
| RCTD | 0.686 | 0.570 |
| Cell2Location | 0.636 | 0.578 |
| CITEgeist QP | 0.513 | 0.509 |
| Tangram | 0.480 | 0.420 |

**Note:** "CITEgeist NB" refers to the NB-emission model variant used in this benchmark. The current production solver is cuOPT QP (Module 3). See `CITEgeist/README.md` for the current pipeline.

## Directory Structure

```
xenium_benchmarking/
├── ground_truth_singler/          # Canonical SingleR GT
├── evaluation/
│   ├── results/canonical_results_all.json
│   ├── src/
│   └── consolidate_results.py
├── CITEgeist/
│   ├── output_nb_singler_{6,7}type/   # NB model (canonical)
│   ├── output_singler_{6,7}type/      # QP baseline
│   ├── output_hybrid_{6,7}type/       # Hybrid pipeline
│   ├── src/
│   ├── slurm/
│   └── _archive/
├── Cell2Location/output_singler_{6,7}type/
├── RCTD/output_singler_{6,7}type/
├── Seurat/output_singler_{6,7}type/
├── Tangram/output_singler_{6,7}type/
└── reference_data/
```

## Reproducing

```bash
sbatch CITEgeist/slurm/sbatch_nb_emission.sh         # NB 6-type
sbatch CITEgeist/slurm/sbatch_nb_7type_benchmark.sh   # NB 7-type
python evaluation/consolidate_results.py               # Consolidate
```

## Archives

Deprecated outputs in `_archive/` subdirectories. Not used in manuscript.
