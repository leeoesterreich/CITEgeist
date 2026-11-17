#!/bin/bash
#SBATCH --job-name=protein_gt_gex
#SBATCH --time=02:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --cluster=htc
#SBATCH --partition=htc

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
set -e

echo "=============================================="
echo "Generate GEX Ground Truth (Protein-Gated)"
echo "=============================================="
echo "Start time: $(date)"

CITEGEIST_ENV="${CITEGEIST_ENVS_ROOT:?CITEGEIST_ENVS_ROOT is not set and conda was not found; set it to your conda envs directory}/CITEgeist_env"
export PATH="${CITEGEIST_ENV}/bin:$PATH"
source activate "${CITEGEIST_ENV}"

SRC_DIR="${CITEGEIST_REPO_ROOT}/Benchmarking/xenium_pseudovisium/src"

python "${SRC_DIR}/generate_protein_gt_gex.py"

echo ""
echo "Verifying output..."
ls -la ${CITEGEIST_REPO_ROOT}/Benchmarking/xenium_pseudovisium/data_protein_gt/ground_truth_gex/
for i in 0 1 2 3 4; do
    echo "Region ${i}:"
    ls ${CITEGEIST_REPO_ROOT}/Benchmarking/xenium_pseudovisium/data_protein_gt/ground_truth_gex/Xenium_region_${i}/
done

echo ""
echo "End time: $(date)"
