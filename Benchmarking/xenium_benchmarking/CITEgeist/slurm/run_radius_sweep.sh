#!/bin/bash
#SBATCH --job-name=cg_radius_sweep
#SBATCH --output=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/CITEgeist/slurm/slurm_log/radius_sweep_%A_%a.out
#SBATCH --error=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/CITEgeist/slurm/slurm_log/radius_sweep_%A_%a.err
#SBATCH --array=0-4
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alc376@pitt.edu

# Radius sweep benchmark for CITEgeist Module 3
# Tests radii [50, 105, 205, 305] corresponding to 0, 1, 2, 3 rings
# Array job: regions 0-4 (each task runs all 4 radii)

set -e
module load gurobi/12.0.3

REPO="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist"
export PATH="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env/bin:$PATH"
source activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

REGION_ID=${SLURM_ARRAY_TASK_ID}

echo "=============================================="
echo "CITEgeist Radius Sweep - Region ${REGION_ID}"
echo "Testing radii: 50 (0 rings), 105 (1 ring), 205 (2 rings), 305 (3 rings)"
echo "=============================================="
echo "Start time: $(date)"

python "${REPO}/Benchmarking/xenium_benchmarking/CITEgeist/src/run_radius_sweep.py" \
    --region-id ${REGION_ID} \
    --input-dir "${REPO}/Benchmarking/xenium_pseudovisium/data_protein_gt" \
    --output-dir "${REPO}/Benchmarking/xenium_benchmarking/CITEgeist/output/radius_sweep"

echo ""
echo "Region ${REGION_ID} complete. End time: $(date)"
