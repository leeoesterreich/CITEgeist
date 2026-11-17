#!/bin/bash
#SBATCH --job-name=cg_no_unk
#SBATCH --output=slurm_log/cg_no_unk_%j.out
#SBATCH --error=slurm_log/cg_no_unk_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --cluster=htc
#SBATCH --partition=htc

# ============================================================================
# CITEgeist Benchmark: Achievable-7 WITHOUT Unknown profile
# ============================================================================
# Region 0 only (for GEX comparison figure)

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
REPO_ROOT="${CITEGEIST_REPO_ROOT}"
PSEUDOVISIUM_DIR="${REPO_ROOT}/Benchmarking/xenium_pseudovisium"
BENCH_DIR="${REPO_ROOT}/Benchmarking/xenium_benchmarking"
CITEGEIST_DIR="${BENCH_DIR}/CITEgeist"
INPUT_DIR="${PSEUDOVISIUM_DIR}/data_granular_gt"
OUTPUT_DIR="${CITEGEIST_DIR}/output_achievable_7_no_unknown"

mkdir -p "${OUTPUT_DIR}"

# Environment
export PATH="${CITEGEIST_ENVS_ROOT:?CITEGEIST_ENVS_ROOT is not set and conda was not found; set it to your conda envs directory}/CITEgeist_env/bin:$PATH"
source activate ${CITEGEIST_ENVS_ROOT:?CITEGEIST_ENVS_ROOT is not set and conda was not found; set it to your conda envs directory}/CITEgeist_env
module load gurobi/12.0.3

cd "${REPO_ROOT}"

echo "Running CITEgeist region 0: no Unknown, radius=150, NO Laplacian smoothing..."
python "${CITEGEIST_DIR}/src/run_benchmark.py" \
    --region-id 0 \
    --mode manual \
    --input-dir "${INPUT_DIR}" \
    --output-dir "${OUTPUT_DIR}" \
    --radius 250.0 \
    --lambda-reg 1.0 \
    --alpha-elastic 0.7 \
    --max-y-change 0.4 \
    --min-counts 25 \
    --lambda-laplacian 0.0 \
    --lambda-reg-gex 0.01 \
    --run-gex \
    --no-unknown

echo "Done!"
