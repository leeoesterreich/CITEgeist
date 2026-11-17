#!/bin/bash
#SBATCH --job-name=a12_prodcovHO
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=06:00:00
#SBATCH --array=0-3
#SBATCH --output=Benchmarking/midkine_remediation/slurm/logs/prodcovHO_%A_%a.out
#SBATCH --error=Benchmarking/midkine_remediation/slurm/logs/prodcovHO_%A_%a.err
# A12 rev 8 held-out prospective run: same four cases, seeds 70000000..70000499, same criterion. Submit ONLY after the criterion commit exists: A12_CRITERION_COMMIT=<hash> sbatch Benchmarking/midkine_remediation/slurm/prod_coverage_heldout.sbatch
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
: "${A12_CRITERION_COMMIT:?set A12_CRITERION_COMMIT}"
export MPLCONFIGDIR=/tmp OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd ${CITEGEIST_REPO_ROOT}
git merge-base --is-ancestor "$A12_CRITERION_COMMIT" HEAD
CASES=(null_autocorrelated null_autocorrelated_stress null_depth_confounded null_noisy_nst_kshare)
C=${CASES[$SLURM_ARRAY_TASK_ID]}
"$(cg_env_python CITEgeist_env)" -u Benchmarking/midkine_remediation/src/run_prod_coverage.py --heldout --case "$C" --procs "${SLURM_CPUS_PER_TASK}"
