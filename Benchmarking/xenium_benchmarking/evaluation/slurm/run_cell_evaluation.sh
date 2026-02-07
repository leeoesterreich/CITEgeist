#!/bin/bash
#SBATCH --job-name=cell_classification_eval
#SBATCH --output=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/evaluation/slurm/slurm_log/cell_classification_eval_%j.out
#SBATCH --error=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/evaluation/slurm/slurm_log/cell_classification_eval_%j.err
#SBATCH --time=03:00:00
#SBATCH --mem=96G
#SBATCH --cpus-per-task=4
#SBATCH --partition=smp
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alc376@pitt.edu

echo "Starting cell classification evaluation at $(date)"

# Activate conda environment
source /ihome/alee/alc376/.bashrc
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist

EVAL_BASE=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/evaluation/results

# Run three-tier evaluation
# NOTE: Requires Module 3, ScType, and Seurat results to already exist.
# Run run_sctype_annotation.sh and run_seurat_label_transfer.sh first.
python Benchmarking/xenium_benchmarking/evaluation/src/evaluate_cell_classification.py \
    --module3_dir ${EVAL_BASE}/module3_classification \
    --sctype_dir ${EVAL_BASE}/sctype_annotation \
    --seurat_dir ${EVAL_BASE}/seurat_label_transfer \
    --output_dir ${EVAL_BASE}/cell_classification_evaluation \
    --verbose

echo "Finished at $(date)"
