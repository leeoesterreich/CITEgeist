#!/bin/bash
#SBATCH --job-name=geomx_tangram
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:1
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=Benchmarking/geomx_benchmarking/Tangram/slurm/tangram_%j.out
# GPU_RACE_READONLY_OK=1

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
module load cuda/11.8.0

ENV="${CITEGEIST_ENVS_ROOT:?CITEGEIST_ENVS_ROOT is not set and conda was not found; set it to your conda envs directory}/tangram_env"
REPO="${CITEGEIST_REPO_ROOT}"
REF="$REPO/Benchmarking/xenium_benchmarking/reference_data/GSE156632/processed/tangram/reference.h5ad"

source activate "$ENV"

python -u "$REPO/Benchmarking/geomx_benchmarking/Tangram/src/run_pseudo_geomx.py" \
    --input-h5ad "$REPO/Benchmarking/geomx_benchmarking/data/aoi_gex.h5ad" \
    --ref-path "$REF" \
    --output-dir "$REPO/Benchmarking/geomx_benchmarking/Tangram/output" \
    --n-markers 100 \
    --num-epochs 1000

echo "Done: $(date)"
