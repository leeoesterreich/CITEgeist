#!/bin/bash
#SBATCH --job-name=eval_gex
#SBATCH --time=02:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=2
#SBATCH --cluster=htc
#SBATCH --partition=htc

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

eval "$(conda shell.bash hook 2>/dev/null)"
conda activate ${CITEGEIST_ENVS_ROOT:?CITEGEIST_ENVS_ROOT is not set and conda was not found; set it to your conda envs directory}/CITEgeist_env
set -u

cd ${CITEGEIST_REPO_ROOT}

SRC_DIR="Benchmarking/xenium_benchmarking/evaluation/src"
OUT_DIR="Benchmarking/xenium_benchmarking/evaluation"

echo "=== Xenium GEX Method Comparison ==="
echo "Start: $(date)"

python -u "${SRC_DIR}/compare_xenium_gex_methods.py" \
    --output_dir "${OUT_DIR}"

# Mirror outputs to the figure's read path (Fig 4D reads results/gex_rerun_v8/).
FIG_DIR="${OUT_DIR}/results/gex_rerun_v8"
mkdir -p "${FIG_DIR}"
cp -f "${OUT_DIR}/xenium_gex_comparison_7type_summary.csv" "${FIG_DIR}/"
cp -f "${OUT_DIR}/xenium_gex_comparison_7type.csv" "${FIG_DIR}/"
echo "Mirrored summary+detail CSVs to ${FIG_DIR}"

echo "DONE: $(date)"
