#!/bin/bash
#SBATCH --job-name=citegeist_preprocess
#SBATCH --cluster=gpu
#SBATCH --partition=l40s
#SBATCH --gres=gpu:1
#SBATCH --mem=128G
#SBATCH --time=4:00:00
#SBATCH --output=slurm/logs/preprocess_%j.out
#SBATCH --error=slurm/logs/preprocess_%j.err

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
set -eo pipefail

eval "$(conda shell.bash hook)"
conda activate ${CITEGEIST_ENVS_ROOT:?CITEGEIST_ENVS_ROOT is not set and conda was not found; set it to your conda envs directory}/CITEgeist_env
set -u  # Re-enable nounset after conda activation

# Paths
XENIUM_MORPH="${CITEGEIST_DATA_ROOT}/10x_PublicData/Xenium_RNA_Proteomic_RenalCellCarcinoma/morphology_focus"
DAPI="${XENIUM_MORPH}/ch0000_dapi.ome.tif"
BOUNDARY="${XENIUM_MORPH}/ch0001_atp1a1_cd45_e-cadherin.ome.tif"

SCRIPT_DIR="${CITEGEIST_REPO_ROOT}/Benchmarking/xenium_benchmarking/CITEgeist"
INPUT_DIR="${CITEGEIST_REPO_ROOT}/Benchmarking/xenium_pseudovisium/data_protein_gt"
OUTPUT_DIR="${SCRIPT_DIR}/output/preprocessing"

echo "Starting preprocessing at $(date)"
python "${SCRIPT_DIR}/src/run_preprocess.py" \
    --dapi "${DAPI}" \
    --boundary "${BOUNDARY}" \
    --input-dir "${INPUT_DIR}" \
    --output-dir "${OUTPUT_DIR}" \
    --skip-segmentation \
    --skip-patches \
    --epochs 300 \
    --batch-size 256

echo "Preprocessing complete at $(date)"
