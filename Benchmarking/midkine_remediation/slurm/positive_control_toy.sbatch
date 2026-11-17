#!/bin/bash
#SBATCH --job-name=a12_pc_toy
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=Benchmarking/midkine_remediation/slurm/logs/pctoy_%j.out
#SBATCH --error=Benchmarking/midkine_remediation/slurm/logs/pctoy_%j.err
# Positive-control COMMOT rerun on the SYNTHETIC toy specimen only (envs/COMMOT; one `import commot`, 100-220 s on NFS). Not a real-data run.
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
export MPLCONFIGDIR=/tmp HDF5_USE_FILE_LOCKING=FALSE NUMBA_CACHE_DIR=/tmp
cd ${CITEGEIST_REPO_ROOT}
OUT=${PC_TOY_OUT:-Benchmarking/midkine_remediation/synthetic/positive_control_toy.json}
PY="$(cg_env_python COMMOT)"
"$PY" -u Benchmarking/midkine_remediation/src/positive_control.py --toy --out "$OUT"
echo "wrote $OUT"
