#!/bin/bash
#SBATCH --job-name=geomx_hnscc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --gres=gpu:1
#SBATCH --output=geomx_benchmarking/CITEgeist/slurm/geomx_hnscc_%j.out
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
export PYTHONUNBUFFERED=1
export HDF5_USE_FILE_LOCKING=FALSE
export NUMBA_CACHE_DIR=/tmp/numba_cache_${SLURM_JOB_ID}
export MPLBACKEND=Agg
export MPLCONFIGDIR=/tmp/mpl_${SLURM_JOB_ID}
mkdir -p /tmp/mpl_${SLURM_JOB_ID}
cp -r ~/.cache/matplotlib/* /tmp/mpl_${SLURM_JOB_ID}/ 2>/dev/null || true

module load cuda/12.6

REPO_ROOT="${CITEGEIST_REPO_ROOT}"
ENV="${CITEGEIST_ENVS_ROOT:?CITEGEIST_ENVS_ROOT is not set and conda was not found; set it to your conda envs directory}/CITEgeist_env"
DATA="${REPO_ROOT}/geomx_benchmarking/data/hnscc"
OUTPUT="${REPO_ROOT}/geomx_benchmarking/CITEgeist/output/geomx_hnscc"

source activate ${ENV}
cd ${REPO_ROOT}

echo "=== Step 1: Run CITEgeist pipeline ==="

python -u geomx_benchmarking/CITEgeist/src/run_geomx_hnscc.py \
    --gex ${DATA}/geomx_gex.h5ad \
    --protein ${DATA}/geomx_protein.h5ad \
    --output-dir ${OUTPUT} \
    --segment-col segment \
    --n-markers 20 \
    --holdout-frac 0.5

echo "=== Step 2: Evaluate ==="
python -u geomx_benchmarking/CITEgeist/src/evaluate_geomx_hnscc.py \
    --pipeline-dir ${OUTPUT} \
    --protein-adata ${DATA}/geomx_protein.h5ad \
    --sace-gex ${OUTPUT}/sace_cell_gex.h5ad \
    --holdout-file ${OUTPUT}/holdout_aois.csv \
    --segment-col segment \
    --output ${REPO_ROOT}/geomx_benchmarking/CITEgeist/evaluation/geomx_hnscc_results.json

echo "=== Done ==="
