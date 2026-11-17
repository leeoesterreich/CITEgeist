#!/bin/bash
#SBATCH --job-name=v2_p0
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --array=0-11
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=Benchmarking/midkine_remediation/slurm/logs/v2_p0_%A_%a.log
# v2 breast chain step 1: the FRESH protein start p0 (closed-form marker share; stamped), one specimen per array task, CPU, no cuOPT.
# Env: CODE_DIR (clean detached worktree at HEAD), P0_DIR (NEW abs dir = <module3 dir>/p0).  Submitted by .commands/submit-v2-breast-chain.sh p0
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
SAMPLES=(HCC22-088-P1-S1 HCC22-088-P1-S2 HCC22-088-P2-S1 HCC22-088-P2-S2 HCC22-088-P3-S1_A HCC22-088-P3-S2
         HCC22-088-P4-S1 HCC22-088-P4-S2_1i_rep HCC22-088-P5-S1 HCC22-088-P5-S2_F_rep HCC22-088-P6-S1 HCC22-088-P6-S2_D)
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}
cd "${CODE_DIR:?CODE_DIR required}"
export HDF5_USE_FILE_LOCKING=FALSE MPLCONFIGDIR=/tmp/mpl_$SLURM_JOB_ID NUMBA_CACHE_DIR=/tmp/numba_$SLURM_JOB_ID OMP_NUM_THREADS=8
PY="$(cg_env_python CITEgeist_env)"
echo "v2 p0 ${SAMPLE} code=$(git rev-parse HEAD) start=$(date -Is)"
"$PY" -u Benchmarking/midkine_remediation/src/v2_inputs.py p0 --specimen "$SAMPLE" --p0-dir "${P0_DIR:?P0_DIR required}"
echo "done ${SAMPLE} end=$(date -Is)"
