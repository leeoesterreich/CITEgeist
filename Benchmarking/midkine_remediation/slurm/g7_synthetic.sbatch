#!/bin/bash
#SBATCH --job-name=a12_g7
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=08:00:00
#SBATCH --array=0-3
#SBATCH --output=Benchmarking/midkine_remediation/slurm/logs/g7_%A_%a.out
#SBATCH --error=Benchmarking/midkine_remediation/slurm/logs/g7_%A_%a.err
# A12.3 G7 synthetic Moran gates; task 0 = definition C grid, 1 = definition A, 2 = runtime/memory at n~5000.
# Optional env: G7_N (pairs per setting), G7_OUT (output dir), G7_LATTICES (comma list, task 0 only).
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
export MPLCONFIGDIR=/tmp
cd ${CITEGEIST_REPO_ROOT}
RUNS=(c a bench cg)
RUN=${RUNS[$SLURM_ARRAY_TASK_ID]}
echo "run=$RUN node=$(hostname) start=$(date)"
ARGS=(--run "$RUN" --procs "${SLURM_CPUS_PER_TASK}")
[ -n "${G7_N:-}" ] && ARGS+=(--n "$G7_N")
[ -n "${G7_OUT:-}" ] && ARGS+=(--out "$G7_OUT")
[ -n "${G7_LATTICES:-}" ] && ARGS+=(--lattices "$G7_LATTICES")
if [ "$RUN" = "c" ]; then export OMP_NUM_THREADS=16; else export OMP_NUM_THREADS=16; fi
"$(cg_env_python CITEgeist_env)" -u Benchmarking/midkine_remediation/src/run_g7.py "${ARGS[@]}"
echo "done $(date)"
