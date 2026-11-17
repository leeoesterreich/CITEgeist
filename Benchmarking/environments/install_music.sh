#!/bin/bash
#SBATCH --job-name=install_music
#SBATCH --time=01:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --cluster=htc
#SBATCH --partition=htc

# Install MuSiC R package (required for CARD reference-based deconvolution)

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
echo "Installing MuSiC R package"
echo "=============================================="
echo "Start time: $(date)"
echo "Node: $(hostname)"

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate ${CITEGEIST_ENVS_ROOT:?CITEGEIST_ENVS_ROOT is not set and conda was not found; set it to your conda envs directory}/CARD_env

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
