#!/bin/bash
#SBATCH --job-name=c2l_train_p7
#SBATCH --time=12:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1

# Train Cell2Location reference model on GSE156632 with 7 protein-GT types

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
echo "race namespace=${GPU_RACE_NAMESPACE:-none} partition=${GPU_RACE_PARTITION:-none}"

module load cuda/11.8.0
echo "Loaded CUDA: $(nvcc --version 2>/dev/null | grep release || echo 'CUDA loaded')"

BASE_DIR="${CITEGEIST_REPO_ROOT}/Benchmarking/xenium_benchmarking"
C2L_DIR="${BASE_DIR}/Cell2Location"
SRC_DIR="${C2L_DIR}/src"
REF_DATA="${REF_DATA:-${BASE_DIR}/reference_data/GSE156632/processed_protein7/cell2location/reference.h5ad}"
OUTPUT_DIR="${OUTPUT_DIR:-${C2L_DIR}/reference_model_protein7}"

echo "=============================================="
echo "Cell2Location Reference Training (7 Protein-GT Types)"
echo "=============================================="
echo "Start time: $(date)"
echo "Reference data: ${REF_DATA}"
echo "Output directory: ${OUTPUT_DIR}"

if [ ! -f "${REF_DATA}" ]; then
    echo "ERROR: Reference data not found at ${REF_DATA}"
    echo "Please run process_reference_protein7.sh first"
    exit 1
fi

# Use Cell2Location environment Python directly (no conda activation needed)
PYTHON="$(cg_env_python cell2location_env)"

echo "Python: ${PYTHON}"
${PYTHON} -c "import torch; print(f'CUDA available: {torch.cuda.is_available()}'); print(f'GPU: {torch.cuda.get_device_name(0) if torch.cuda.is_available() else None}')"
echo ""

${PYTHON} "${SRC_DIR}/train_reference.py" \
    --ref-path "${REF_DATA}" \
    --output-dir "${OUTPUT_DIR}" \
    --batch-key "sample" \
    --labels-key "cell_type" \
    --max-epochs 250 \
    --filter-genes

echo ""
echo "=============================================="
echo "Training Complete!"
echo "=============================================="
echo "End time: $(date)"
ls -la "${OUTPUT_DIR}/"
