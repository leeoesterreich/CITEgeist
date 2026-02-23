#!/bin/bash
#SBATCH --job-name=eval_asma_fix
#SBATCH --output=eval_asma_fix_%j.out
#SBATCH --error=eval_asma_fix_%j.err
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

# Evaluate CITEgeist with αSMA-only fibroblast ground truth

export PATH="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env/bin:$PATH"
source activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist

EVAL_DIR="Benchmarking/xenium_benchmarking/evaluation"
GT_DIR="Benchmarking/xenium_pseudovisium/data_protein_gt"
PRED_DIR="Benchmarking/xenium_benchmarking/CITEgeist/output/manual"

echo "=== Evaluating CITEgeist with αSMA-only Fibroblast GT ==="
echo "GT dir: ${GT_DIR}"
echo "Pred dir: ${PRED_DIR}"

python ${EVAL_DIR}/src/evaluate_benchmark.py \
    --gt-dir ${GT_DIR} \
    --pred-dir ${PRED_DIR} \
    --n-regions 5 \
    --output ${EVAL_DIR}/results_asma_only_gt.json

echo "=== Done ==="
cat ${EVAL_DIR}/results_asma_only_gt.json
