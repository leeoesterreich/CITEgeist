#!/bin/bash
#SBATCH --job-name=sv_singler
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=06:00:00
#SBATCH --cluster=htc
#SBATCH --partition=htc
# Usage: VARIANT={control|A|B} sbatch run_singler_variant.sbatch
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
set -euo pipefail
XB="${CITEGEIST_REPO_ROOT}/Benchmarking/xenium_benchmarking"
SV="${XB}/reference_audit/singler_variants"
"${CITEGEIST_ENVS_ROOT:?CITEGEIST_ENVS_ROOT is not set and conda was not found; set it to your conda envs directory}/r_deconv_env/bin/Rscript" "${SV}/run_singler_variant.R" \
    --xenium-counts "${XB}/ground_truth_singler/singler_labels/xenium_cell_gex.csv" \
    --ref-mtx-dir "${XB}/reference_data/GSE156632/processed/seurat/mtx" \
    --ref-labels "${SV}/labels/ref_labels_${VARIANT}.csv" \
    --output-dir "${SV}/singler_out/${VARIANT}" \
    --ncpus 16
