#!/bin/bash
# Setup script for CARD benchmarking environment
#
# Usage:
#   ./setup_card_env.sh
#
# This script:
# 1. Creates the conda environment from card_env.yml
# 2. Installs CARD from GitHub
# 3. Installs MuSiC (dependency for reference-based mode)

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "============================================"
echo "Setting up CARD environment"
echo "============================================"

# Create conda environment
echo "Creating conda environment from card_env.yml..."
conda env create -f "${SCRIPT_DIR}/card_env.yml" || {
    echo "Environment may already exist. Trying to update..."
    conda env update -f "${SCRIPT_DIR}/card_env.yml"
}

# Activate environment
echo "Activating CARD_env..."
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate CARD_env

# Install MuSiC from GitHub (required dependency)
echo "Installing MuSiC from GitHub..."
Rscript -e 'if (!requireNamespace("MuSiC", quietly = TRUE)) devtools::install_github("xuranw/MuSiC")'

# Install CARD from GitHub
echo "Installing CARD from GitHub..."
Rscript -e 'devtools::install_github("YingMa0107/CARD")'

# Verify installation
echo ""
echo "Verifying installation..."
Rscript -e 'library(CARD); cat("CARD version:", as.character(packageVersion("CARD")), "\n")'

echo ""
echo "============================================"
echo "CARD environment setup complete!"
echo "============================================"
echo ""
echo "To activate: conda activate CARD_env"
