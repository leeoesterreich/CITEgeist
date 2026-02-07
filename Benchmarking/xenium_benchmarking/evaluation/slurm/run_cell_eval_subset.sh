#!/bin/bash
#SBATCH --job-name=cell_eval_subset
#SBATCH --output=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/evaluation/slurm/slurm_log/cell_eval_subset_%j.out
#SBATCH --error=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/evaluation/slurm/slurm_log/cell_eval_subset_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=96G
#SBATCH --cpus-per-task=4
#SBATCH --partition=smp
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alc376@pitt.edu

source /ihome/alee/alc376/.bashrc
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist

EVAL_BASE=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/evaluation/results
GATING_DIR=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/CITEgeist/output_gating_validation

python Benchmarking/xenium_benchmarking/evaluation/src/run_cell_evaluation_subset.py \
    --gating_output_dir ${GATING_DIR} \
    --sctype_dir ${EVAL_BASE}/sctype_annotation \
    --seurat_dir ${EVAL_BASE}/seurat_label_transfer \
    --output_dir ${EVAL_BASE}/cell_classification_evaluation_both \
    --verbose

