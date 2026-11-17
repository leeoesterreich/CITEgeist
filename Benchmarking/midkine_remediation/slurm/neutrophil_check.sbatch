#!/bin/bash
#SBATCH --job-name=neutro_check
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=Benchmarking/midkine_remediation/slurm/logs/neutrophil_check_%j.log
# EVIDENCE ONLY (Neutrophils dropped from v41; needs the superseded 8-type proportions): Neutrophils plausibility check on the marker-share proportions. Env: OUT (json path), optional MODULE3_DIR.
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
cd ${CITEGEIST_REPO_ROOT}
export HDF5_USE_FILE_LOCKING=FALSE MPLCONFIGDIR=/tmp/mpl_$SLURM_JOB_ID NUMBA_CACHE_DIR=/tmp/numba_$SLURM_JOB_ID OMP_NUM_THREADS=4
PY="$(cg_env_python CITEgeist_env)"
"$PY" -u Benchmarking/midkine_remediation/src/neutrophil_check.py \
  --module3-dir "${MODULE3_DIR:-output/module3_marker_share_8type}" \
  --data-root ${CITEGEIST_DATA_ROOT}/GSE289326/processed_files \
  --out "${OUT:?OUT required}"
echo "done end=$(date -Is)"
