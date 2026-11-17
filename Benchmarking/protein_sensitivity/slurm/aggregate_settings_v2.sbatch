#!/bin/bash
#SBATCH --job-name=settings_v2_report
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=Benchmarking/protein_sensitivity/logs/%x_%j.log
# A7b round 2: run the fail-closed settings_v2 report (gates, OAT, joint LORO) into a fresh output dir.
# Report lives in the main repo, so SC.REPO (default ref_root) is the main repo: no --ref-root flag needed.
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
set -euo pipefail
cd ${CITEGEIST_REPO_ROOT}
PY="$(cg_env_python CITEgeist_env)"
OUT=${1:-Benchmarking/protein_sensitivity/results/settings_v2_report}
$PY -u Benchmarking/protein_sensitivity/settings_v2_report.py --out "$OUT"
