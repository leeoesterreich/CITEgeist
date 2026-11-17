#!/bin/bash
#SBATCH --job-name=umap_validation
#SBATCH --partition=htc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=6:00:00

# UMAP Validation of Xenium 10-Cell-Type Ground Truth
# This script validates the RNA k-means clustering ground truth by
# generating UMAP visualizations and computing clustering metrics.
#
# Expected runtime: 2-4 hours for 465K cells
# Output: figures/ directory with PNGs and JSON metrics

# Use absolute paths
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
PROJECT_DIR="${CITEGEIST_REPO_ROOT}/Benchmarking/xenium_pseudovisium"
CONDA_ENV="${CITEGEIST_ENVS_ROOT:?CITEGEIST_ENVS_ROOT is not set and conda was not found; set it to your conda envs directory}/CITEgeist_env"

echo "=============================================="
echo "UMAP VALIDATION OF XENIUM GROUND TRUTH"
echo "=============================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Memory: $SLURM_MEM_PER_NODE MB"
echo "Start time: $(date)"
echo "=============================================="

# Change to project directory
cd "$PROJECT_DIR"
echo "Working directory: $(pwd)"

# Activate conda environment using full path
source ${HOME}/.bashrc
source $(conda info --base)/etc/profile.d/conda.sh
conda activate "$CONDA_ENV"

# Verify Python and key packages
echo ""
echo "Environment check:"
which python
python --version
python -c "import scanpy; print(f'scanpy: {scanpy.__version__}')"
python -c "import sklearn; print(f'sklearn: {sklearn.__version__}')"
echo ""

# Run validation script with absolute path
echo "Running UMAP validation..."
python "${PROJECT_DIR}/src/validate_ground_truth_umap.py"

# Check exit status
if [ $? -eq 0 ]; then
    echo ""
    echo "=============================================="
    echo "VALIDATION COMPLETE"
    echo "=============================================="
    echo "Output files in: figures/"
    ls -la figures/
else
    echo ""
    echo "=============================================="
    echo "VALIDATION FAILED"
    echo "=============================================="
    exit 1
fi

echo ""
echo "End time: $(date)"
