#!/bin/bash
#SBATCH --job-name=v2_refit
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --array=0-11
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
# v2 breast chain step 3 (stage B1): the frozen RNA refinement of the fresh p0 on the rebuilt nuclei, ALL p0 spots, one specimen per array task.
# Env: CODE_DIR (clean detached worktree at HEAD), P0_DIR, NUCLEI_DIR (stamped B2 dir <root>/v2b2), MODULE3_DIR (NEW; receives the module3 layout + stamps).
# Submitted by .commands/submit-v2-breast-chain.sh refit
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
set -euo pipefail
PY="$(cg_env_python CITEgeist_env)"
SAMPLES=(HCC22-088-P1-S1 HCC22-088-P1-S2 HCC22-088-P2-S1 HCC22-088-P2-S2 HCC22-088-P3-S1_A HCC22-088-P3-S2
         HCC22-088-P4-S1 HCC22-088-P4-S2_1i_rep HCC22-088-P5-S1 HCC22-088-P5-S2_F_rep HCC22-088-P6-S1 HCC22-088-P6-S2_D)
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}
export HDF5_USE_FILE_LOCKING=FALSE MPLCONFIGDIR=/tmp/mpl_$SLURM_JOB_ID NUMBA_CACHE_DIR=/tmp/numba_$SLURM_JOB_ID
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4 PYTHONDONTWRITEBYTECODE=1
cd "${CODE_DIR:?CODE_DIR required}"
echo "CPU job; v2 refit $SAMPLE code=$(git --no-optional-locks rev-parse HEAD) start=$(date -Is)"
"$PY" -u Benchmarking/midkine_remediation/src/rna_refined_props.py refit --specimen "$SAMPLE" \
  --p0-dir "${P0_DIR:?P0_DIR required}" --nuclei-dir "${NUCLEI_DIR:?NUCLEI_DIR required}" --module3-dir "${MODULE3_DIR:?MODULE3_DIR required}"
echo "done end=$(date -Is)"
