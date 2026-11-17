#!/bin/bash
#SBATCH --job-name=a12_g6fix
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=04:00:00
#SBATCH --array=0-19
#SBATCH --output=Benchmarking/midkine_remediation/slurm/logs/g6fix_%A_%a.out
#SBATCH --error=Benchmarking/midkine_remediation/slurm/logs/g6fix_%A_%a.err
# A12.3 G6 fix round: 4 variants x (4 null cases + injection). Rule: g6_fix_round_selection_rule.md. Optional env: G6F_N, G6F_OUT.
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
VARIANTS=(v_b3_nocov v_b3_cov v_b6_nocov v_b6_cov)
CASES=(null_autocorrelated null_autocorrelated_stress null_noisy_nst_kshare null_depth_confounded INJ)
V=${VARIANTS[$((SLURM_ARRAY_TASK_ID / 5))]}
C=${CASES[$((SLURM_ARRAY_TASK_ID % 5))]}
ARGS=(--variant "$V" --procs "${SLURM_CPUS_PER_TASK}")
if [ "$C" = INJ ]; then ARGS+=(--inj); else ARGS+=(--case "$C"); fi
[ -n "${G6F_N:-}" ] && ARGS+=(--n "$G6F_N")
[ -n "${G6F_OUT:-}" ] && ARGS+=(--out "$G6F_OUT")
echo "variant=$V case=$C node=$(hostname) start=$(date)"
"$(cg_env_python CITEgeist_env)" -u Benchmarking/midkine_remediation/src/run_g6fix.py "${ARGS[@]}"
echo "done $(date)"
