#! /usr/bin/bash
#SBATCH --job-name=Seurat_mx_deconv
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --array=0-4  # Adjust this based on the number of replicates

# Load R environment
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
eval "$(conda shell.bash hook 2>/dev/null)"
conda activate ${CITEGEIST_ENVS_ROOT:?CITEGEIST_ENVS_ROOT is not set and conda was not found; set it to your conda envs directory}/r_deconv_env

SCRIPT_DIR="${CITEGEIST_REPO_ROOT}/Benchmarking/simulation_benchmarking/Seurat/mixed"
cd "$SCRIPT_DIR"

# Define input directory
INPUT_DIR="${CITEGEIST_REPO_ROOT}/Benchmarking/simulation_benchmarking/replicates/mixed/h5ad_objects"
OUTPUT_H5SEURAT_DIR="${CITEGEIST_REPO_ROOT}/Benchmarking/simulation_benchmarking/Seurat/mixed/h5seurat"
OUTPUT_DIR="${CITEGEIST_REPO_ROOT}/Benchmarking/simulation_benchmarking/Seurat/mixed/output"

mkdir -pv "$OUTPUT_H5SEURAT_DIR"

# Get list of h5ad files
FILES=($(ls $INPUT_DIR/*GEX.h5ad))

# Ensure we do not go out of bounds
if [ "$SLURM_ARRAY_TASK_ID" -ge "${#FILES[@]}" ]; then
    echo "Error: SLURM_ARRAY_TASK_ID $SLURM_ARRAY_TASK_ID is out of range for available files."
    exit 1
fi

# Get the specific file for this array job
FILE=${FILES[$SLURM_ARRAY_TASK_ID]}

# Extract the replicate name from the file
BASENAME=$(basename "$FILE" _GEX.h5ad)

# Capture start time
START_TIME=$(date +%s)
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Running Seurat deconvolution for replicate: $BASENAME"

# Run the Seurat R script for this replicate
Rscript --vanilla Seurat_mx.R "$FILE" "$OUTPUT_DIR" "$BASENAME"

# Capture end time
END_TIME=$(date +%s)
RUNTIME=$((END_TIME - START_TIME))
RUNTIME_MINUTES=$(echo "scale=2; $RUNTIME / 60" | bc)

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Completed Seurat deconvolution for replicate: $BASENAME"
echo "SEURAT_TOTAL_RUNTIME: $BASENAME took $RUNTIME seconds ($RUNTIME_MINUTES minutes)"

# Move converted h5Seurat files (if any exist)
mv "$INPUT_DIR"/${BASENAME}*.h5seurat "$OUTPUT_H5SEURAT_DIR" 2>/dev/null || true
