#!/bin/bash
#SBATCH --job-name=card_setup_markers
#SBATCH --time=01:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=4
#SBATCH --cluster=htc
#SBATCH --partition=htc

# Generate marker genes for CARD reference-free mode
# Run this once before running run_all_regions_reffree.sh
#
# Uses the reference scRNA-seq data to identify top marker genes per cell type.

_cg_find() {
    if [ -n "${CITEGEIST_PATHS_SH:-}" ]; then echo "${CITEGEIST_PATHS_SH}"; return 0; fi
    local p d
    p="$(cd "$(dirname "${BASH_SOURCE[0]}")" 2>/dev/null && pwd)/../../../_paths.sh"
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

# Directories
BASE_DIR="${CITEGEIST_REPO_ROOT}/Benchmarking/xenium_benchmarking"
CARD_DIR="${BASE_DIR}/CARD"
SRC_DIR="${CARD_DIR}/src"
REF_FILE="${BASE_DIR}/reference_data/GSE156632/processed_protein7/rctd/reference.h5ad"
OUTPUT_FILE="${CARD_DIR}/marker_genes_protein7.csv"

echo "=============================================="
echo "CARD Marker Gene Generation"
echo "=============================================="
echo "Start time: $(date)"
echo ""
echo "Reference: ${REF_FILE}"
echo "Output: ${OUTPUT_FILE}"
echo ""

# Check if reference exists
if [ ! -f "${REF_FILE}" ]; then
    echo "ERROR: Reference h5ad not found at ${REF_FILE}"
    exit 1
fi

# Environment
CITEGEIST_ENV="${CITEGEIST_ENVS_ROOT:?CITEGEIST_ENVS_ROOT is not set and conda was not found; set it to your conda envs directory}/CITEgeist_env"

# Activate CITEgeist_env (has scanpy)
echo "Activating CITEgeist_env..."
export PATH="${CITEGEIST_ENV}/bin:$PATH"
source activate "${CITEGEIST_ENV}"

# Generate markers
python "${SRC_DIR}/generate_marker_genes.py" \
    --reference "${REF_FILE}" \
    --output "${OUTPUT_FILE}" \
    --cell-type-col "cell_type" \
    --n-markers 50

echo ""
echo "=============================================="
echo "Marker Generation Complete!"
echo "=============================================="
echo "End time: $(date)"
echo ""
echo "Generated marker file:"
ls -la "${OUTPUT_FILE}"
echo ""
echo "Preview:"
head -20 "${OUTPUT_FILE}"
