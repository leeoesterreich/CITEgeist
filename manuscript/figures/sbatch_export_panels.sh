#!/bin/bash
#SBATCH --job-name=export_panels
#SBATCH --output=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/manuscript/figures/slurm_log/export_panels_%j.out
#SBATCH --error=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/manuscript/figures/slurm_log/export_panels_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alc376@pitt.edu

set -e

REPO="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist"
FIG_DIR="${REPO}/manuscript/figures"

export PATH="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env/bin:$PATH"
source activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

cd "${REPO}"

echo "=============================================="
echo "Exporting individual panels"
echo "Start time: $(date)"
echo "=============================================="

python "${FIG_DIR}/export_panels.py"

echo ""
echo "=============================================="
echo "Panel export complete"
echo "End time: $(date)"
echo "=============================================="
