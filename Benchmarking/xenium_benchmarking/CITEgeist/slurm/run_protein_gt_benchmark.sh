#!/bin/bash
#SBATCH --job-name=cg_protein_gt
#SBATCH --output=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/CITEgeist/slurm/slurm_log/protein_gt_%A_%a.out
#SBATCH --error=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/CITEgeist/slurm/slurm_log/protein_gt_%A_%a.err
#SBATCH --array=0-4
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alc376@pitt.edu

# Run CITEgeist manual benchmark against protein-gated ground truth
# Array job: regions 0-4

set -e
module load gurobi/12.0.3

REPO="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist"
export PATH="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env/bin:$PATH"
source activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

REGION_ID=${SLURM_ARRAY_TASK_ID}

echo "=============================================="
echo "CITEgeist Protein GT Benchmark - Region ${REGION_ID}"
echo "=============================================="
echo "Start time: $(date)"

python "${REPO}/Benchmarking/xenium_benchmarking/CITEgeist/src/run_benchmark.py" \
    --region-id ${REGION_ID} \
    --mode manual \
    --input-dir "${REPO}/Benchmarking/xenium_pseudovisium/data_protein_gt" \
    --output-dir "${REPO}/Benchmarking/xenium_benchmarking/CITEgeist/output_protein_gt" \
    --run-gex

echo ""
echo "Region ${REGION_ID} complete. End time: $(date)"
