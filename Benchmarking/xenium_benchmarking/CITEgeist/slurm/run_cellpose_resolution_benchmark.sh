#!/bin/bash
#SBATCH --job-name=cg_cellpose_res
#SBATCH --output=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/CITEgeist/slurm/slurm_log/cellpose_res_%A_%a.out
#SBATCH --error=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/CITEgeist/slurm/slurm_log/cellpose_res_%A_%a.err
#SBATCH --array=0-4
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=24:00:00
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alc376@pitt.edu

# Benchmark Cellpose nuclei-count fidelity across lowres/hires/fullres.
# Default array covers Xenium pseudo-Visium regions 0-4.

set -eo pipefail

# Use submission directory so SLURM spool-script location does not break paths.
REPO_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
SRC_SCRIPT="${REPO_ROOT}/Benchmarking/xenium_benchmarking/CITEgeist/src/benchmark_cellpose_resolution.py"
INPUT_DIR="${REPO_ROOT}/Benchmarking/xenium_pseudovisium/data_protein_gt/h5ad_objects"
MAPPING_CSV="${REPO_ROOT}/Benchmarking/xenium_pseudovisium/data_protein_gt/cell_to_spot_mapping.csv"
FALLBACK_DATA_ROOT="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_pseudovisium/data_protein_gt"
if [[ ! -d "${INPUT_DIR}" ]]; then
  INPUT_DIR="${FALLBACK_DATA_ROOT}/h5ad_objects"
fi
if [[ ! -f "${MAPPING_CSV}" ]]; then
  MAPPING_CSV="${FALLBACK_DATA_ROOT}/cell_to_spot_mapping.csv"
fi
OUT_DIR="${REPO_ROOT}/Benchmarking/xenium_benchmarking/CITEgeist/output_cellpose_resolution"
LOG_DIR="${REPO_ROOT}/Benchmarking/xenium_benchmarking/CITEgeist/slurm/slurm_log"
LOWRES_IMAGE_ROOT="${REPO_ROOT}/Benchmarking/xenium_benchmarking/scResolve/images/morphology"
HIRES_IMAGE_ROOT="${REPO_ROOT}/Benchmarking/xenium_benchmarking/scResolve/images/morphology_hires"
COORD_INFO_ROOT="${REPO_ROOT}/Benchmarking/xenium_benchmarking/scResolve/images/morphology_hires"
if [[ ! -d "${HIRES_IMAGE_ROOT}" ]]; then
  LOWRES_IMAGE_ROOT="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/scResolve/images/morphology"
  HIRES_IMAGE_ROOT="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/scResolve/images/morphology_hires"
  COORD_INFO_ROOT="${HIRES_IMAGE_ROOT}"
fi
FULLRES_IMAGE_PATH=${FULLRES_IMAGE_PATH:-"/ix1/alee/LO_LAB/General/Public_Data/10x_PublicData/Xenium_RNA_Proteomic_RenalCellCarcinoma/morphology.ome.tif"}

mkdir -p "${LOG_DIR}" "${OUT_DIR}"

eval "$(conda shell.bash hook)"
conda activate /ihome/alee/alc376/alc376_bgfs/envs/CITEgeist_env
module load gurobi/12.0.3
export CELLPOSE_LOCAL_MODELS_PATH="${CELLPOSE_LOCAL_MODELS_PATH:-/scratch/slurm-${SLURM_JOB_ID}/cellpose_models}"
mkdir -p "${CELLPOSE_LOCAL_MODELS_PATH}"

REGION_ID=${SLURM_ARRAY_TASK_ID}
RESOLUTIONS=${RESOLUTIONS:-"lowres hires fullres"}
CELLPOSE_GPU=${CELLPOSE_GPU:-0}
CELLPOSE_DIAMETER=${CELLPOSE_DIAMETER:-""}
MAX_FULLRES_SIDE=${MAX_FULLRES_SIDE:-9000}
FULLRES_PATCH_RADIUS_MULTIPLIER=${FULLRES_PATCH_RADIUS_MULTIPLIER:-1.5}
FULLRES_PATCH_WORKERS=${FULLRES_PATCH_WORKERS:-4}
# Standard Visium spot diameter in microns
SPOT_DIAMETER_UM=${SPOT_DIAMETER_UM:-55.0}

GPU_FLAG=""
if [[ "${CELLPOSE_GPU}" == "1" ]]; then
  GPU_FLAG="--cellpose-gpu"
fi

DIAM_FLAG=""
if [[ -n "${CELLPOSE_DIAMETER}" ]]; then
  DIAM_FLAG="--cellpose-diameter ${CELLPOSE_DIAMETER}"
fi

cd "${REPO_ROOT}"

echo "===================================================="
echo "Cellpose resolution benchmark"
echo "Region: ${REGION_ID}"
echo "Resolutions: ${RESOLUTIONS}"
echo "Spot diameter: ${SPOT_DIAMETER_UM} µm"
echo "Output: ${OUT_DIR}"
echo "Start: $(date)"
echo "===================================================="

python "${SRC_SCRIPT}" \
  --region-id "${REGION_ID}" \
  --h5ad-dir "${INPUT_DIR}" \
  --mapping-csv "${MAPPING_CSV}" \
  --output-dir "${OUT_DIR}" \
  --lowres-image-root "${LOWRES_IMAGE_ROOT}" \
  --hires-image-root "${HIRES_IMAGE_ROOT}" \
  --coord-info-root "${COORD_INFO_ROOT}" \
  --fullres-image-path "${FULLRES_IMAGE_PATH}" \
  --resolutions ${RESOLUTIONS} \
  --spot-diameter-um "${SPOT_DIAMETER_UM}" \
  --fullres-patch-radius-multiplier "${FULLRES_PATCH_RADIUS_MULTIPLIER}" \
  --fullres-patch-workers "${FULLRES_PATCH_WORKERS}" \
  ${GPU_FLAG} \
  ${DIAM_FLAG} \
  --max-fullres-side "${MAX_FULLRES_SIDE}"

echo "Done: $(date)"
