#!/bin/bash
#SBATCH --job-name=install_music
#SBATCH --output=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/environments/install_music_%j.out
#SBATCH --error=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/environments/install_music_%j.err
#SBATCH --time=01:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

# Install MuSiC R package (required for CARD reference-based deconvolution)

set -e

echo "=============================================="
echo "Installing MuSiC R package"
echo "=============================================="
echo "Start time: $(date)"
echo "Node: $(hostname)"

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CARD_env

echo "R: $(which R)"
echo ""

# Install MuSiC from GitHub
echo "Installing MuSiC dependencies and MuSiC..."
Rscript -e '
options(repos = c(CRAN = "https://cloud.r-project.org"))
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
# Install all Bioconductor dependencies including TOAST (required by MuSiC)
BiocManager::install(c("Biobase", "SingleCellExperiment", "SummarizedExperiment", "TOAST"), update = FALSE, ask = FALSE)
devtools::install_github("xuranw/MuSiC", upgrade = "never")
'

# Verify installation
echo ""
echo "Verifying MuSiC installation..."
Rscript -e 'library(MuSiC); cat("MuSiC version:", as.character(packageVersion("MuSiC")), "\n")'

echo ""
echo "=============================================="
echo "MuSiC installation complete!"
echo "=============================================="
echo "End time: $(date)"
