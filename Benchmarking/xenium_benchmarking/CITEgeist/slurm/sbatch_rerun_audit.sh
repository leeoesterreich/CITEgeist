#!/bin/bash
#SBATCH --job-name=nb_audit_rerun
#SBATCH --output=slurm_%j.log
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G

_cg_find() {
    if [ -n "${CITEGEIST_PATHS_SH:-}" ]; then echo "${CITEGEIST_PATHS_SH}"; return 0; fi
    local p d
    p="$(cd "$(dirname "${BASH_SOURCE[0]}")" 2>/dev/null && pwd)/../../../_paths.sh"
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

module load gurobi/12.0.3

source activate ${CITEGEIST_ENVS_ROOT:?CITEGEIST_ENVS_ROOT is not set and conda was not found; set it to your conda envs directory}/CITEgeist_env

set -u

REPO="${CITEGEIST_REPO_ROOT}"
cd "$REPO"

echo "=== Re-running NB canonical benchmark (CD68 fix) ==="
echo "Start: $(date)"

# Run all 5 regions with fresh predictions
python Benchmarking/xenium_benchmarking/CITEgeist/src/run_nb_canonical_benchmark.py \
    --output-dir Benchmarking/xenium_benchmarking/CITEgeist/output_nb_singler_6type \
    --n-iter 100 \
    --n-starts 3

echo "=== Benchmark complete, running audit ==="

# Run comprehensive audit
python Benchmarking/xenium_benchmarking/evaluation/src/audit_qp_pipeline.py

echo "=== Done: $(date) ==="
