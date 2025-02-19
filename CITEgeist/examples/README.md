# CITEgeist Examples

This directory contains example scripts and configuration files that demonstrate how to use CITEgeist for spatial transcriptomics analysis.

## Contents

### `compute_sample.py`
A comprehensive example script that demonstrates how to:
- Load and preprocess Visium spatial transcriptomics data
- Apply cell type deconvolution using predefined cell profiles
- Run the CITEgeist model for spatial gene expression analysis
- Generate visualizations and save results

Usage:
```bash
python compute_sample.py --path /path/to/visium/data --radius 400 --min_counts 100 \
    --output_folder output --gurobi_license /path/to/gurobi.lic
```

### `sbatch_sample.sh`
A SLURM batch script template for running CITEgeist analysis on a high-performance computing cluster. This script:
- Sets up the necessary environment and dependencies
- Processes multiple samples in parallel using SLURM array jobs
- Handles resource allocation and job management

To use this script:
1. Modify the paths at the top of the script to match your installation:
   ```bash
   CITEGEIST_FOLDER="/path/to/local/install/CITEgeist/CITEgeist"
   GUROBI_LICENSE_FILE="/path/to/your/gurobi.lic"
   ```
2. Adjust the SLURM parameters as needed for your computing environment
3. Submit the job using:
   ```bash
   sbatch sbatch_sample.sh
   ```

### Interactive Vignettes

The examples directory includes three Jupyter notebook vignettes that demonstrate real-world applications of CITEgeist:

1. `vignette_1_biopsy_heterogeneity.ipynb`
   - Analysis of tumor heterogeneity in breast cancer biopsies
   - Demonstrates spatial profiling of different cell populations
   - Shows how to interpret CITEgeist results in a clinical context

2. `vignette_2_surgical_d538g.ipynb`
   - Focused analysis of ESR1 D538G mutant breast cancer samples
   - Explores the spatial organization of immune cells around tumor regions
   - Includes pathway analysis and mutation-specific patterns

3. `vignette_3_responder_macrophages.ipynb`
   - Investigation of macrophage populations in treatment response
   - Shows how to analyze specific immune cell subtypes
   - Demonstrates correlation of spatial patterns with clinical outcomes

These notebooks provide step-by-step walkthroughs with detailed explanations and visualizations. They serve as both tutorials and templates for your own analyses.

## Requirements
- Python environment with required packages (scanpy, squidpy, numpy, pandas)
- Gurobi optimizer with valid license
- CITEgeist package installed in your environment

For detailed information about CITEgeist and its methodology, please refer to the main README in the root directory. 
