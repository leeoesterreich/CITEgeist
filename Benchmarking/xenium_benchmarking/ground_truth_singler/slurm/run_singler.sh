#!/bin/bash
#SBATCH --job-name=singler_xenium
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=06:00:00
#SBATCH --cluster=htc
#SBATCH --partition=htc

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
set -euo pipefail

BASE_DIR="${CITEGEIST_REPO_ROOT}/Benchmarking/xenium_benchmarking"
GT_DIR="${BASE_DIR}/ground_truth_singler"
LABELS_DIR="${GT_DIR}/singler_labels"
REF_MTX_DIR="${BASE_DIR}/reference_data/GSE156632/processed/seurat/mtx"
REF_LABELS="${BASE_DIR}/reference_data/GSE156632/processed/rctd/cell_types.csv"

# Step 1: Extract Xenium cell-level data (Python)
echo "=== Step 1: Extracting Xenium cell-level data ==="
PY_ENV="${CITEGEIST_ENVS_ROOT:?CITEGEIST_ENVS_ROOT is not set and conda was not found; set it to your conda envs directory}/CITEgeist_env"
export PATH="${PY_ENV}/bin:$PATH"
eval "$(conda shell.bash hook)"
set +u  # conda activate.d binutils script references unbound ADDR2LINE under `set -u`
conda activate "${PY_ENV}"
set -u

python "${GT_DIR}/src/extract_xenium_data.py" \
    --output-dir "${LABELS_DIR}"

echo "=== Step 1 complete: GEX and protein CSVs extracted ==="

# Step 2: Run SingleR (R)
echo "=== Step 2: Running SingleR label transfer ==="
R_ENV="${CITEGEIST_ENVS_ROOT:?CITEGEIST_ENVS_ROOT is not set and conda was not found; set it to your conda envs directory}/r_deconv_env"

"${R_ENV}/bin/Rscript" "${GT_DIR}/src/run_singler.R" \
    --xenium-counts "${LABELS_DIR}/xenium_cell_gex.csv" \
    --ref-mtx-dir "${REF_MTX_DIR}" \
    --ref-labels "${REF_LABELS}" \
    --output-dir "${LABELS_DIR}" \
    --ncpus 16

echo "=== SingleR complete ==="
