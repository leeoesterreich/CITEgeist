#!/bin/bash
#SBATCH --job-name=unified_gex_sim
#SBATCH --output=unified_gex_sim_%j.log
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G

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

eval "$(conda shell.bash hook)"
conda activate ${CITEGEIST_ENVS_ROOT:?CITEGEIST_ENVS_ROOT is not set and conda was not found; set it to your conda envs directory}/CITEgeist_env

set -u

SCRIPT_DIR="${CITEGEIST_REPO_ROOT}/Benchmarking/simulation_benchmarking/CITEgeist/src"
OUTPUT_DIR="${CITEGEIST_REPO_ROOT}/Benchmarking/simulation_benchmarking/CITEgeist/results_unified_gex"
AGG_SCRIPT="${CITEGEIST_REPO_ROOT}/Benchmarking/xenium_benchmarking/CITEgeist/src/aggregate_unified_results.py"

cd "$SCRIPT_DIR"

for rep in 0 1 2 3 4; do
    echo "=== Running replicate $rep ==="
    python run_sace_simulation_benchmark.py --rep "$rep" --output_dir "$OUTPUT_DIR"
    echo "=== Replicate $rep done ==="
done

echo "=== Aggregating ==="
python "$AGG_SCRIPT" \
    --input_dir "$OUTPUT_DIR" \
    --prefix simulation_unified_gex_rep \
    --n 5

echo "=== ALL DONE ==="
