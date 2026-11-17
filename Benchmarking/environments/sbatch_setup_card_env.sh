#!/bin/bash
#SBATCH --job-name=setup_card_env
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --cluster=htc
#SBATCH --partition=htc

# Setup CARD conda environment
# This installs R, dependencies, and CARD from GitHub

_cg_find() {
    if [ -n "${CITEGEIST_PATHS_SH:-}" ]; then echo "${CITEGEIST_PATHS_SH}"; return 0; fi
    local p d
    p="$(cd "$(dirname "${BASH_SOURCE[0]}")" 2>/dev/null && pwd)/../_paths.sh"
    if [ -f "$p" ]; then echo "$p"; return 0; fi
    d="${SLURM_SUBMIT_DIR:-}"
    while [ -n "$d" ] && [ "$d" != "/" ] && [ "$d" != "." ]; do
        if [ -f "$d/Benchmarking/_paths.sh" ]; then echo "$d/Benchmarking/_paths.sh"; return 0; fi
        d="$(dirname "$d")"
    done
    echo "cannot find Benchmarking/_paths.sh; set CITEGEIST_PATHS_SH" >&2
    return 1
}
_cg_paths="$(_cg_find)" || exit 1
source "${_cg_paths}"
unset -f _cg_find
unset _cg_paths
set -e

echo "=============================================="
echo "Setting up CARD environment"
echo "=============================================="
echo "Start time: $(date)"
echo "Node: $(hostname)"
echo ""

SCRIPT_DIR="${CITEGEIST_REPO_ROOT}/Benchmarking/environments"
ENV_PATH="${CITEGEIST_ENVS_ROOT:?CITEGEIST_ENVS_ROOT is not set and conda was not found; set it to your conda envs directory}/CARD_env"

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
