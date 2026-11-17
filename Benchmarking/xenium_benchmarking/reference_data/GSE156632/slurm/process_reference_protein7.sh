#!/bin/bash
#SBATCH --job-name=process_ref_p7
#SBATCH --time=12:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=16
#SBATCH --cluster=htc
#SBATCH --partition=htc

# Process GSE156632 reference data with 7 cell types matching protein GT
# Output: processed_protein7/ with cell2location, tangram, rctd, seurat subdirs

_cg_find() {
    if [ -n "${CITEGEIST_PATHS_SH:-}" ]; then echo "${CITEGEIST_PATHS_SH}"; return 0; fi
    local p d
    p="$(cd "$(dirname "${BASH_SOURCE[0]}")" 2>/dev/null && pwd)/../../../../_paths.sh"
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

BASE_DIR="${CITEGEIST_REPO_ROOT}/Benchmarking/xenium_benchmarking/reference_data/GSE156632"
SRC_DIR="${BASE_DIR}/src"
INPUT_DIR="${BASE_DIR}"
OUTPUT_DIR="${BASE_DIR}/processed_protein7"

mkdir -p "${OUTPUT_DIR}"

echo "=============================================="
echo "Processing GSE156632 Reference (7 Protein-GT Types)"
echo "=============================================="
echo "Start time: $(date)"
echo "Output directory: ${OUTPUT_DIR}"

# Environment
ENV_BASE="${CITEGEIST_ENVS_ROOT:?CITEGEIST_ENVS_ROOT is not set and conda was not found; set it to your conda envs directory}"
CELL2LOC_ENV="${ENV_BASE}/cell2location_env"
CITEGEIST_ENV="${ENV_BASE}/CITEgeist_env"

if [ -d "${CELL2LOC_ENV}" ]; then
    echo "Activating cell2location_env..."
    export PATH="${CELL2LOC_ENV}/bin:$PATH"
    source activate "${CELL2LOC_ENV}"
elif [ -d "${CITEGEIST_ENV}" ]; then
    echo "Using CITEgeist_env..."
    export PATH="${CITEGEIST_ENV}/bin:$PATH"
    source activate "${CITEGEIST_ENV}"
else
    echo "ERROR: No suitable conda environment found"
    exit 1
fi

echo "Python: $(which python)"
echo ""

python "${SRC_DIR}/process_reference.py" \
    --input-dir "${INPUT_DIR}" \
    --output-dir "${OUTPUT_DIR}" \
    --min-genes 200 \
    --max-genes 5000 \
    --max-mt-pct 20.0

echo ""
echo "=============================================="
echo "Processing Complete!"
echo "=============================================="
echo "End time: $(date)"
echo ""
echo "Output:"
ls -la "${OUTPUT_DIR}/"
for subdir in cell2location tangram rctd seurat; do
    echo "${subdir}:"
    ls -la "${OUTPUT_DIR}/${subdir}/" 2>/dev/null || echo "  (not created)"
done
