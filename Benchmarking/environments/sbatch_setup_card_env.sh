#!/bin/bash
#SBATCH --job-name=setup_card_env
#SBATCH --output=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/environments/setup_card_env_%j.out
#SBATCH --error=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/environments/setup_card_env_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

# Setup CARD conda environment
# This installs R, dependencies, and CARD from GitHub

set -e

echo "=============================================="
echo "Setting up CARD environment"
echo "=============================================="
echo "Start time: $(date)"
echo "Node: $(hostname)"
echo ""

SCRIPT_DIR="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/environments"
ENV_PATH="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CARD_env"

cd "${SCRIPT_DIR}"

# Initialize conda
source "$(conda info --base)/etc/profile.d/conda.sh"

# Create conda environment with explicit prefix
echo "Creating conda environment at ${ENV_PATH}..."
conda env create -f card_env.yml -p "${ENV_PATH}" || {
    echo "Environment may already exist. Trying to update..."
    conda env update -f card_env.yml -p "${ENV_PATH}"
}

# Activate environment
echo ""
echo "Activating CARD_env..."
conda activate "${ENV_PATH}"

echo "R: $(which R)"
echo "R version: $(R --version | head -1)"

# Install MuSiC from GitHub (CARD dependency)
echo ""
echo "Installing MuSiC from GitHub..."
Rscript -e 'if (!requireNamespace("MuSiC", quietly = TRUE)) devtools::install_github("xuranw/MuSiC")'

# Install CARD from GitHub
echo ""
echo "Installing CARD from GitHub..."
Rscript -e 'devtools::install_github("YingMa0107/CARD")'

# Verify installation
echo ""
echo "Verifying installation..."
Rscript -e 'library(CARD); cat("CARD version:", as.character(packageVersion("CARD")), "\n")'

echo ""
echo "=============================================="
echo "CARD environment setup complete!"
echo "=============================================="
echo "End time: $(date)"
echo ""
echo "Environment location: ${ENV_PATH}"
echo "To activate: conda activate ${ENV_PATH}"
