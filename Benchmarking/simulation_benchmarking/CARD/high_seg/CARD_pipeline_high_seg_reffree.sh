#!/bin/bash
#SBATCH --job-name=CARD_hs_rf
#SBATCH --output=./logs/CARD_hs_reffree_%a.out
#SBATCH --error=./logs/CARD_hs_reffree_%a.err
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --array=0-4
#SBATCH --cluster=htc
#SBATCH --partition=htc

# Run CARD deconvolution on high_seg simulation (REFERENCE-FREE mode)
# Array job: processes replicates 0-4
#
# Prerequisites: marker_genes_simulation.csv must exist in parent directory

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
set -e

# Record start time
START_TIMESTAMP=$(date +%s)
START_TIME=$(date +'%Y-%m-%d %H:%M:%S')
echo "[$START_TIME] Job started for replicate $SLURM_ARRAY_TASK_ID (ref-free)"

# Directories
SCRIPT_DIR="${CITEGEIST_REPO_ROOT}/Benchmarking/simulation_benchmarking/CARD/high_seg"
OUTPUT_DIR="${SCRIPT_DIR}"

# Create log directory
mkdir -p "${SCRIPT_DIR}/logs"

cd "${SCRIPT_DIR}"

# Environment
CARD_ENV="${CITEGEIST_ENVS_ROOT:?CITEGEIST_ENVS_ROOT is not set and conda was not found; set it to your conda envs directory}/CARD_env"

# Activate CARD environment
echo "Activating CARD_env..."
export PATH="${CARD_ENV}/bin:$PATH"
source activate "${CARD_ENV}"

echo "R: $(which R)"
echo "R version: $(R --version | head -1)"

# Set replicate index from array ID
REPLICATE_INDEX=$SLURM_ARRAY_TASK_ID

# Run CARD (reference-free mode)
Rscript --vanilla CARD_pipeline_high_seg.R \
    --replicates $REPLICATE_INDEX \
    --output_dir $OUTPUT_DIR \
    --mode reference_free

# Record end time
END_TIMESTAMP=$(date +%s)
END_TIME=$(date +'%Y-%m-%d %H:%M:%S')
RUNTIME=$((END_TIMESTAMP - START_TIMESTAMP))
RUNTIME_MINUTES=$(echo "scale=2; $RUNTIME / 60" | bc)

echo "[$END_TIME] Job completed for replicate $SLURM_ARRAY_TASK_ID (ref-free)"
echo "CARD_REFFREE_TOTAL_RUNTIME: Replicate $SLURM_ARRAY_TASK_ID took $RUNTIME seconds ($RUNTIME_MINUTES minutes)"
