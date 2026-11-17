#!/bin/bash
#SBATCH --job-name=analyze_clusters
#SBATCH --output=slurm_log/analyze_clusters_%j.out
#SBATCH --error=slurm_log/analyze_clusters_%j.err
#SBATCH --time=01:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --cluster=htc
#SBATCH --partition=htc

# Analyze RNA cluster protein profiles for cell type annotation

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
cd ${CITEGEIST_REPO_ROOT}/Benchmarking/xenium_pseudovisium

# Set up environment
export PATH="${CITEGEIST_ENVS_ROOT:?CITEGEIST_ENVS_ROOT is not set and conda was not found; set it to your conda envs directory}/CITEgeist_env/bin:$PATH"
PYTHON="$(cg_env_python CITEgeist_env)"

echo "Starting cluster profile analysis..."
echo "Job ID: $SLURM_JOB_ID"
echo "Time: $(date)"
echo "Python: $PYTHON"

$PYTHON analyze_cluster_profiles.py

echo "Completed at: $(date)"
