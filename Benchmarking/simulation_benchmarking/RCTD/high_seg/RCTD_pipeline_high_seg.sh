#!/bin/bash
#SBATCH --job-name=RCTD_hs
#SBATCH --output=./logs/RCTD_hs_%a.out
#SBATCH --error=./logs/RCTD_hs_%a.err
#SBATCH --time=00:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --array=0-4

# Record the start time
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
START_TIMESTAMP=$(date +%s)
START_TIME=$(date +'%Y-%m-%d %H:%M:%S')
echo "[$START_TIME] Job started for replicate $SLURM_ARRAY_TASK_ID" | tee -a ./logs/RCTD_hs_runtime.log

cd ${CITEGEIST_REPO_ROOT}/Benchmarking/simulation_benchmarking/RCTD/high_seg

# Activate R environment with spacexr
eval "$(conda shell.bash hook 2>/dev/null)"
conda activate ${CITEGEIST_ENVS_ROOT:?CITEGEIST_ENVS_ROOT is not set and conda was not found; set it to your conda envs directory}/r_deconv_env

# Set replicate index from job array ID
REPLICATE_INDEX=$SLURM_ARRAY_TASK_ID

# Define output directory
OUTPUT_DIR="${CITEGEIST_REPO_ROOT}/Benchmarking/simulation_benchmarking/RCTD/high_seg"

# Run the R script for the specific replicate
Rscript --vanilla RCTD_pipeline_high_seg.R --replicates $REPLICATE_INDEX --output_dir $OUTPUT_DIR

# Record the end time
END_TIMESTAMP=$(date +%s)
END_TIME=$(date +'%Y-%m-%d %H:%M:%S')

# Calculate total runtime
RUNTIME=$((END_TIMESTAMP - START_TIMESTAMP))
RUNTIME_MINUTES=$(echo "scale=2; $RUNTIME / 60" | bc)

echo "[$END_TIME] Job completed for replicate $SLURM_ARRAY_TASK_ID" | tee -a ./logs/RCTD_hs_runtime.log
echo "RCTD_TOTAL_RUNTIME: Replicate $SLURM_ARRAY_TASK_ID took $RUNTIME seconds ($RUNTIME_MINUTES minutes)" | tee -a ./logs/RCTD_hs_runtime.log
