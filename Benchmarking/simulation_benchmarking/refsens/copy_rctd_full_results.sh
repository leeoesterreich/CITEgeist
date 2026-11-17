#!/bin/bash
#SBATCH --job-name=copy_rctd
#SBATCH --time=00:10:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
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
cd ${CITEGEIST_REPO_ROOT}/Benchmarking/simulation_benchmarking/refsens

for cond in high_seg mixed; do
    for f in results/RCTD/${cond}/full/Wu_ST_*_RCTD_deconv_predictions.csv; do
        if [ -f "$f" ]; then
            dir=$(dirname "$f")
            base=$(basename "$f" .csv)
            cp "$f" "${dir}/${base}_v2.csv" && echo "Copied: ${dir}/${base}_v2.csv"
        fi
    done
done
echo "DONE"
