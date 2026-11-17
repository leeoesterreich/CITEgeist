#!/bin/bash
#SBATCH --job-name=geomx_pseudo
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --gres=gpu:1
#SBATCH --output=geomx_benchmarking/CITEgeist/slurm/pseudo_geomx_%j.out
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
module load cuda/12.6

REPO_ROOT="${CITEGEIST_REPO_ROOT}"
ENV="${CITEGEIST_ENVS_ROOT:?CITEGEIST_ENVS_ROOT is not set and conda was not found; set it to your conda envs directory}/CITEgeist_env"
OUTPUT="${REPO_ROOT}/geomx_benchmarking/CITEgeist/output/pseudo_geomx"

source activate ${ENV}
cd ${REPO_ROOT}

echo "=== Step 1: Create pseudo-GeoMx AOIs ==="
python -u geomx_benchmarking/CITEgeist/src/create_pseudo_geomx.py \
    --output-dir ${OUTPUT}/data \
    --n-regions 5 \
    --aois-per-segment 10 \
    --aoi-diameter 400.0 \
    --min-separation 400.0 \
    --min-cells 10

echo "=== Step 2: Run CITEgeist pipeline (known profiles) ==="
python -u geomx_benchmarking/CITEgeist/src/run_pseudo_geomx.py \
    --gex ${OUTPUT}/data/aoi_gex.h5ad \
    --protein ${OUTPUT}/data/aoi_protein.h5ad \
    --output-dir ${OUTPUT}/pipeline \
    --skip-m1m2 \
    --profile-dict ${REPO_ROOT}/geomx_benchmarking/CITEgeist/src/xenium_7type_profiles.json

echo "=== Step 3: Evaluate ==="
python -u geomx_benchmarking/CITEgeist/src/evaluate_pseudo_geomx.py \
    --pred-dir ${OUTPUT}/pipeline \
    --gt-dir ${OUTPUT}/data \
    --output ${REPO_ROOT}/geomx_benchmarking/CITEgeist/evaluation/pseudo_geomx_results.json

echo "=== Done ==="
