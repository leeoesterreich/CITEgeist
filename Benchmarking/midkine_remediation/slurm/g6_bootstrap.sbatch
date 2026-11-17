#!/bin/bash
#SBATCH --job-name=a12_boot
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=06:00:00
#SBATCH --array=0-9
#SBATCH --output=Benchmarking/midkine_remediation/slurm/logs/boot_%A_%a.out
#SBATCH --error=Benchmarking/midkine_remediation/slurm/logs/boot_%A_%a.err
# A12.3 descriptive path: block-bootstrap interval coverage. Task = (case, tile side). Optional env: BOOT_N, BOOT_B, BOOT_OUT.
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
export MPLCONFIGDIR=/tmp OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd ${CITEGEIST_REPO_ROOT}
CASES=(null_autocorrelated null_autocorrelated_stress null_noisy_nst_kshare null_depth_confounded INJ)
SIDES=(6 12)
C=${CASES[$((SLURM_ARRAY_TASK_ID % 5))]}
S=${SIDES[$((SLURM_ARRAY_TASK_ID / 5))]}
ARGS=(--side "$S" --procs "${SLURM_CPUS_PER_TASK}")
if [ "$C" = INJ ]; then ARGS+=(--case null_autocorrelated --inj); else ARGS+=(--case "$C"); fi
# record the env that shapes the run in the log (run_bootstrap.py also writes it into the output json under "settings")
echo "BOOT_N=${BOOT_N:-unset} BOOT_B=${BOOT_B:-unset} BOOT_OUT=${BOOT_OUT:-unset} SLURM_JOB_ID=${SLURM_JOB_ID:-} SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID:-}"
[ -n "${BOOT_N:-}" ] && ARGS+=(--n "$BOOT_N")
[ -n "${BOOT_B:-}" ] && ARGS+=(--B "$BOOT_B")
[ -n "${BOOT_OUT:-}" ] && ARGS+=(--out "$BOOT_OUT")
"$(cg_env_python CITEgeist_env)" -u Benchmarking/midkine_remediation/src/run_bootstrap.py "${ARGS[@]}"
