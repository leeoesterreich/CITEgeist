#!/bin/bash
#SBATCH --job-name=cg_cellpose_gpu
#SBATCH --output=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/CITEgeist/slurm/slurm_log/cellpose_gpu_%A_%a.out
#SBATCH --error=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/CITEgeist/slurm/slurm_log/cellpose_gpu_%A_%a.err
#SBATCH --array=0-4
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --gres=gpu:1
#SBATCH --time=4:00:00
#SBATCH --cluster=gpu
#SBATCH --partition=l40s
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

# Benchmark Cellpose nuclei-count fidelity with GPU acceleration.
# Expected runtime: ~10-30 min per region (vs 8+ hours on CPU)

set -eo pipefail

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

# Load CUDA BEFORE conda to ensure libraries are found
module load cuda/12.8.0
module load gurobi/12.0.3

eval "$(conda shell.bash hook)"
conda activate /ihome/alee/alc376/alc376_bgfs/envs/CITEgeist_env

# Verify GPU is accessible
echo "CUDA_VISIBLE_DEVICES: ${CUDA_VISIBLE_DEVICES:-not set}"
python -c "import torch; print(f'PyTorch CUDA available: {torch.cuda.is_available()}')"
export CELLPOSE_LOCAL_MODELS_PATH="${CELLPOSE_LOCAL_MODELS_PATH:-/scratch/slurm-${SLURM_JOB_ID}/cellpose_models}"
mkdir -p "${CELLPOSE_LOCAL_MODELS_PATH}"

REGION_ID=${SLURM_ARRAY_TASK_ID}
RESOLUTIONS=${RESOLUTIONS:-"lowres hires fullres"}
MAX_FULLRES_SIDE=${MAX_FULLRES_SIDE:-9000}
FULLRES_PATCH_RADIUS_MULTIPLIER=${FULLRES_PATCH_RADIUS_MULTIPLIER:-1.5}
FULLRES_PATCH_WORKERS=${FULLRES_PATCH_WORKERS:-4}
# Standard Visium spot diameter in microns
SPOT_DIAMETER_UM=${SPOT_DIAMETER_UM:-55.0}

cd "${REPO_ROOT}"

echo "===================================================="
echo "Cellpose resolution benchmark (GPU)"
echo "Region: ${REGION_ID}"
echo "Resolutions: ${RESOLUTIONS}"
echo "Spot diameter: ${SPOT_DIAMETER_UM} µm"
echo "GPU: $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null || echo 'N/A')"
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
  --cellpose-gpu \
  --max-fullres-side "${MAX_FULLRES_SIDE}"

echo "Done: $(date)"
