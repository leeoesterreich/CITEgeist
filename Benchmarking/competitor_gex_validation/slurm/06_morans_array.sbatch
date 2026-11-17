#!/bin/bash
#SBATCH --job-name=morans_cgv
#SBATCH --output=Benchmarking/competitor_gex_validation/slurm/logs/06_morans_%A_%a.out
#SBATCH --error=Benchmarking/competitor_gex_validation/slurm/logs/06_morans_%A_%a.err
#SBATCH --array=0-23
#SBATCH --partition=htc
#SBATCH --cluster=htc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=4:00:00

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
set -eo pipefail

REPO_ROOT="${CITEGEIST_REPO_ROOT}"
cd "${REPO_ROOT}"

# Map array index to method + sample
if [ ${SLURM_ARRAY_TASK_ID} -lt 12 ]; then
    METHOD="cell2location"
    SAMPLE_IDX=${SLURM_ARRAY_TASK_ID}
else
    METHOD="tangram"
    SAMPLE_IDX=$((SLURM_ARRAY_TASK_ID - 12))
fi
SAMPLE=$(sed -n "$((SAMPLE_IDX + 1))p" Benchmarking/competitor_gex_validation/slurm/samples.txt)

echo "=== Bivariate Moran's I ==="
echo "Job ID: ${SLURM_JOB_ID}, Array: ${SLURM_ARRAY_TASK_ID}"
echo "Method: ${METHOD}, Sample: ${SAMPLE}"
echo "Node: $(hostname)"
echo "Start: $(date)"

source activate ${CITEGEIST_ENVS_ROOT:?CITEGEIST_ENVS_ROOT is not set and conda was not found; set it to your conda envs directory}/CITEgeist_env

"$(cg_env_python CITEgeist_env)" -u Benchmarking/competitor_gex_validation/run_bivariate_morans.py \
    --method "${METHOD}" \
    --sample "${SAMPLE}"

echo "=== Done: $(date) ==="
