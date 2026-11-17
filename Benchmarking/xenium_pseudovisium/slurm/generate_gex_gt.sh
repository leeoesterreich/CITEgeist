#!/bin/bash
#SBATCH --job-name=gen_gex_gt          # Job name
#SBATCH --output=slurm_log/%x_%j.out   # Standard output
#SBATCH --error=slurm_log/%x_%j.err    # Standard error
#SBATCH --ntasks=1                     # Number of tasks
#SBATCH --cpus-per-task=4              # CPU cores per task
#SBATCH --mem=128G                     # Memory (large for 465K cells)
#SBATCH --time=04:00:00                # Time limit
#SBATCH --cluster=htc                  # Cluster
#SBATCH --partition=htc                # Partition

# ============================================================================
# Generate GEX Ground Truth for Xenium Benchmarking (RNA-based)
# ============================================================================
# Run this BEFORE the main benchmarking pipeline
# Creates per-cell-type gene expression ground truth files
#
# Uses RNA-based cell type classification (recommended) to avoid circular
# logic with protein markers.
#
# Reference:
#   Zhao et al. (2025). "Benchmarking cell type annotation methods for 10x
#   Xenium spatial transcriptomics data." BMC Bioinformatics, 26(1), 25.
#   https://doi.org/10.1186/s12859-025-06044-0
# ============================================================================

_cg_find() {
    if [ -n "${CITEGEIST_PATHS_SH:-}" ]; then echo "${CITEGEIST_PATHS_SH}"; return 0; fi
    local p d
    p="$(cd "$(dirname "${BASH_SOURCE[0]}")" 2>/dev/null && pwd)/../../_paths.sh"
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
REPO_ROOT="${CITEGEIST_REPO_ROOT}"
PSEUDOVISIUM_DIR="${REPO_ROOT}/Benchmarking/xenium_pseudovisium"
OUTPUT_DIR="${PSEUDOVISIUM_DIR}/data_rna_gt"

mkdir -p "${PSEUDOVISIUM_DIR}/slurm/slurm_log"

echo "============================================"
echo "Generate GEX Ground Truth (RNA-based)"
echo "============================================"
echo "Start time: $(date)"
echo "Output directory: ${OUTPUT_DIR}"
echo ""

# Activate environment
export PATH="${CITEGEIST_ENVS_ROOT:?CITEGEIST_ENVS_ROOT is not set and conda was not found; set it to your conda envs directory}/CITEgeist_env/bin:$PATH"
source activate ${CITEGEIST_ENVS_ROOT:?CITEGEIST_ENVS_ROOT is not set and conda was not found; set it to your conda envs directory}/CITEgeist_env

cd "${REPO_ROOT}"

# Run the GEX ground truth generation script with RNA-based cell types (default)
python "${PSEUDOVISIUM_DIR}/src/generate_gex_ground_truth.py" \
    --use-rna-types \
    --output-dir "${OUTPUT_DIR}"

echo ""
echo "============================================"
echo "GEX Ground Truth Generation Complete"
echo "============================================"
echo "End time: $(date)"

# List generated files
echo ""
echo "Generated files:"
ls -la "${OUTPUT_DIR}/ground_truth_gex/" 2>/dev/null || echo "No files generated yet"
for region_dir in "${OUTPUT_DIR}/ground_truth_gex/Xenium_region_"*; do
    if [ -d "$region_dir" ]; then
        echo ""
        echo "$(basename $region_dir):"
        ls -la "$region_dir"
    fi
done
