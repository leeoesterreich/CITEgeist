#!/bin/bash
#SBATCH --job-name=preprocess_v2
#SBATCH --output=slurm/logs/preprocess_v2_%A_%a.out
#SBATCH --error=slurm/logs/preprocess_v2_%A_%a.err
#SBATCH --array=0-4
#SBATCH --time=2:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

# Re-preprocess patches with global normalization (v2)
#
# This script re-runs patch preparation with:
# 1. Global percentile normalization (instead of per-patch z-score)
# 2. Size feature extraction (log1p(w, h, area))
# 3. global_stats.npz and spot_*_sizes.npy outputs
#
# Output directory: output/patches_v2/region_{0-4}

set -e  # Exit on error

# Environment setup
source ~/.bashrc
eval "$(conda shell.bash hook)"
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist

REGION=$SLURM_ARRAY_TASK_ID

# Paths
BASE_DIR="Benchmarking/xenium_benchmarking"
XENIUM_MORPHOLOGY="/ix1/alee/LO_LAB/General/Public_Data/10x_PublicData/Xenium_RNA_Proteomic_RenalCellCarcinoma/morphology_focus"
DAPI_PATH="${XENIUM_MORPHOLOGY}/ch0000_dapi.ome.tif"
BOUNDARY_PATH="${XENIUM_MORPHOLOGY}/ch0001_atp1a1_cd45_e-cadherin.ome.tif"
SC_OUTPUT="${BASE_DIR}/CITEgeist/output/single_cell_resolution"
OUTPUT_BASE="${BASE_DIR}/CITEgeist/output/patches_v2"

# Create output directory
mkdir -p ${OUTPUT_BASE}/region_${REGION}
mkdir -p ${BASE_DIR}/CITEgeist/slurm/logs

MASK_PATH="${SC_OUTPUT}/Xenium_region_${REGION}/Xenium_region_${REGION}_cellpose_mask.npy"
SPOT_MAP="${SC_OUTPUT}/Xenium_region_${REGION}/Xenium_region_${REGION}_nuclei_spot_map.csv"
REGION_OUTPUT="${OUTPUT_BASE}/region_${REGION}"

echo "============================================"
echo "Re-preprocessing Region ${REGION} with v2 pipeline"
echo "============================================"
echo "DAPI: ${DAPI_PATH}"
echo "Boundary: ${BOUNDARY_PATH}"
echo "Mask: ${MASK_PATH}"
echo "Spot map: ${SPOT_MAP}"
echo "Output: ${REGION_OUTPUT}"
echo "Normalization: percentile (1st/99th)"
echo "============================================"

# Check if mask exists
if [ ! -f "${MASK_PATH}" ]; then
    echo "ERROR: Mask not found at ${MASK_PATH}"
    echo "Run single_cell_resolution pipeline first."
    exit 1
fi

if [ ! -f "${SPOT_MAP}" ]; then
    echo "ERROR: Spot map not found at ${SPOT_MAP}"
    echo "Run single_cell_resolution pipeline first."
    exit 1
fi

# Run patch preparation with new preprocessing
python -m Benchmarking.xenium_benchmarking.CITEgeist.src.prepare_patches \
    --dapi "${DAPI_PATH}" \
    --boundary "${BOUNDARY_PATH}" \
    --region ${REGION} \
    --mask "${MASK_PATH}" \
    --nuclei_spot_map "${SPOT_MAP}" \
    --output_dir "${REGION_OUTPUT}" \
    --expansion 0.75 \
    --patch_size 96 \
    --norm-method percentile

echo ""
echo "============================================"
echo "Verifying outputs..."
echo "============================================"

# Verify outputs
if [ ! -f "${REGION_OUTPUT}/global_stats.npz" ]; then
    echo "ERROR: global_stats.npz not created"
    exit 1
fi

# Count patch and size files
N_PATCHES=$(ls ${REGION_OUTPUT}/spot_*_patches.npy 2>/dev/null | wc -l)
N_SIZES=$(ls ${REGION_OUTPUT}/spot_*_sizes.npy 2>/dev/null | wc -l)

echo "Patch files: ${N_PATCHES}"
echo "Size files: ${N_SIZES}"

if [ "${N_PATCHES}" -ne "${N_SIZES}" ]; then
    echo "WARNING: Mismatch between patches (${N_PATCHES}) and sizes (${N_SIZES})"
fi

echo ""
echo "============================================"
echo "Region ${REGION} preprocessing complete!"
echo "============================================"
