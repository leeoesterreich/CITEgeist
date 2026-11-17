#!/bin/bash
#SBATCH --job-name=a12_tests
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=04:00:00
#SBATCH --output=Benchmarking/midkine_remediation/slurm/logs/a12tests_%j.out
#SBATCH --error=Benchmarking/midkine_remediation/slurm/logs/a12tests_%j.err
# A12 tests on htc (pass --mutations to also run the mutation harness; env A12_PYTEST_ARGS narrows the pytest selection, env A12_MUT_ARGS e.g. '--only NAME ...' narrows the mutations). The COMMOT toy tests shell out to envs/COMMOT.
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
export MPLCONFIGDIR=/tmp HDF5_USE_FILE_LOCKING=FALSE
cd ${CITEGEIST_REPO_ROOT}
PY="$(cg_env_python CITEgeist_env)"
"$PY" -u -m pytest ${A12_PYTEST_ARGS:-tests/midkine_remediation} -q -p no:cacheprovider --no-cov -rxXs --durations=10
if [ "${1:-}" = "--mutations" ]; then "$PY" -u tests/midkine_remediation/run_mutations.py ${A12_MUT_ARGS:-}; fi
