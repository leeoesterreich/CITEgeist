#!/bin/bash
#SBATCH --job-name=seurat_label_transfer
#SBATCH --output=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/evaluation/slurm/slurm_log/seurat_label_transfer_%j.out
#SBATCH --error=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/evaluation/slurm/slurm_log/seurat_label_transfer_%j.err
#SBATCH --time=04:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=8
#SBATCH --partition=smp
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alc376@pitt.edu

echo "Starting Seurat label transfer at $(date)"

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist

REF_H5AD=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/reference_data/GSE156632/processed/seurat/reference.h5ad
MTX_DIR=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/reference_data/GSE156632/processed/seurat/mtx
OUTPUT_DIR=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/evaluation/results/seurat_label_transfer

# Step 1: Export reference to 10x MTX format (Python)
echo "=== Step 1: Exporting reference to MTX format ==="
source /ihome/alee/alc376/.bashrc
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

python Benchmarking/xenium_benchmarking/evaluation/src/export_reference_for_seurat.py \
    --reference ${REF_H5AD} \
    --output_dir ${MTX_DIR}

if [ $? -ne 0 ]; then
    echo "ERROR: MTX export failed"
    exit 1
fi
echo "MTX export complete"

# Step 2: Run Seurat label transfer (R)
echo "=== Step 2: Running Seurat label transfer ==="
module load gcc/12.2.0
module load r/4.4.0

Rscript Benchmarking/xenium_benchmarking/evaluation/src/seurat_label_transfer.R \
    --reference ${MTX_DIR} \
    --cell_types ${MTX_DIR}/cell_types.csv \
    --query_dir /ix1/alee/LO_LAB/General/Public_Data/10x_PublicData/Xenium_RNA_Proteomic_RenalCellCarcinoma \
    --output_dir ${OUTPUT_DIR} \
    --n_dims 30

echo "Finished at $(date)"
