#!/bin/bash
#SBATCH --job-name=gex_enrich_diag
#SBATCH --output=Benchmarking/simulation_benchmarking/CITEgeist/slurm/logs/gex_enrichment_diagnostic_%j.log
#SBATCH --time=03:00:00
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

REPO_ROOT="${CITEGEIST_REPO_ROOT}"
SCRIPT="$REPO_ROOT/Benchmarking/simulation_benchmarking/CITEgeist/src/benchmark_gex_enrichment_diagnostic.py"
OUTPUT_DIR="$REPO_ROOT/Benchmarking/simulation_benchmarking/CITEgeist/results_enrichment_diagnostic"

mkdir -p "$OUTPUT_DIR"
mkdir -p "$REPO_ROOT/Benchmarking/simulation_benchmarking/CITEgeist/slurm/logs"

CONDITION="${CONDITION:-mixed}"

echo "=== GEX Enrichment Diagnostic ==="
echo "Condition: $CONDITION"
echo "Output: $OUTPUT_DIR"
echo ""

python "$SCRIPT" --condition "$CONDITION" --output_dir "$OUTPUT_DIR"

echo ""
echo "=== Done ==="
