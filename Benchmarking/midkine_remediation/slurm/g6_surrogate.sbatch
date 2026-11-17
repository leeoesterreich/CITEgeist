#!/bin/bash
#SBATCH --job-name=a12_g6sur
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=04:00:00
#SBATCH --array=0-23
#SBATCH --output=Benchmarking/midkine_remediation/slurm/logs/g6sur_%A_%a.out
#SBATCH --error=Benchmarking/midkine_remediation/slurm/logs/g6sur_%A_%a.err
# A12.3 EXPLORATORY (not protocol): candidate surrogate-label null (sur_*) and 12-pitch tiles (b12_*) for the hurdle.
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
RUNS=(sur_null_autocorrelated sur_null_autocorrelated_stress sur_null_noisy_nst sur_null_depth_confounded inj_sur \
      b12_null_autocorrelated b12_null_autocorrelated_stress b12_null_noisy_nst b12_null_depth_confounded \
      b3_null_autocorrelated b3_null_autocorrelated_stress b3_null_noisy_nst b3_null_depth_confounded \
      b2_null_autocorrelated b2_null_noisy_nst b2_null_depth_confounded b4_null_autocorrelated b4_null_noisy_nst b4_null_depth_confounded \
      b3m_drop_log_n_st_eff b3m_drop_depth b3m_break_blocks)
RUN=${RUNS[$SLURM_ARRAY_TASK_ID]}
echo "run=$RUN node=$(hostname) start=$(date)"
ARGS=(--run "$RUN" --procs "${SLURM_CPUS_PER_TASK}")
[ -n "${G6_N:-}" ] && ARGS+=(--n "$G6_N")
[ -n "${G6_OUT:-}" ] && ARGS+=(--out "$G6_OUT")
"$(cg_env_python CITEgeist_env)" -u Benchmarking/midkine_remediation/src/run_g6.py "${ARGS[@]}"
echo "done $(date)"
