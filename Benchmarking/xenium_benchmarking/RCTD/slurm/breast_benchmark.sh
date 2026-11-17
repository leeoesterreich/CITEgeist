#!/bin/bash
#SBATCH --job-name=rctd_breast
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --array=0-4
#SBATCH --output=slurm/logs/rctd_breast_%A_%a.log

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
module purge
module load gcc/12.2.0 r/4.3.1

REPO_ROOT="${CITEGEIST_REPO_ROOT}"
BASE_DIR="${REPO_ROOT}/Benchmarking/xenium_benchmarking"
REGION_ID=${SLURM_ARRAY_TASK_ID}

REF_COUNTS="${BASE_DIR}/reference_data/Wu_2021/rctd/counts.csv"
REF_CELLTYPES="${BASE_DIR}/reference_data/Wu_2021/rctd/cell_types.csv"
SPATIAL_COUNTS="${BASE_DIR}/data_breast/csv_for_r/Xenium_region_${REGION_ID}/Xenium_counts.csv"
SPATIAL_COORDS="${BASE_DIR}/data_breast/csv_for_r/Xenium_region_${REGION_ID}/Xenium_coords.csv"
OUTPUT_DIR="${BASE_DIR}/RCTD/output_breast"

source activate ${CITEGEIST_ENVS_ROOT:?CITEGEIST_ENVS_ROOT is not set and conda was not found; set it to your conda envs directory}/r_deconv_env

Rscript "${BASE_DIR}/RCTD/src/run_benchmark.R" \
    --region-id "${REGION_ID}" \
    --ref-counts "${REF_COUNTS}" \
    --ref-cell-types "${REF_CELLTYPES}" \
    --spatial-counts "${SPATIAL_COUNTS}" \
    --spatial-coords "${SPATIAL_COORDS}" \
    --output-dir "${OUTPUT_DIR}" \
    --prefix "Xenium" \
    --max-cores 8

echo "RCTD breast region ${REGION_ID} complete"
