#!/bin/bash
#SBATCH --job-name=gen_figures
#SBATCH --output=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/manuscript/figures/slurm_log/generate_all_%j.out
#SBATCH --error=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/manuscript/figures/slurm_log/generate_all_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alc376@pitt.edu

# Generate all manuscript figures

set -e

REPO="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist"
FIG_DIR="${REPO}/manuscript/figures"

export PATH="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env/bin:$PATH"
source activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

cd "${REPO}"

echo "=============================================="
echo "Generating Manuscript Figures"
echo "Start time: $(date)"
echo "=============================================="

for i in 1 2 3 4 5 6; do
    echo ""
    echo "--- Figure ${i} ---"
    python "${FIG_DIR}/generate_figure${i}.py"
    echo "Figure ${i} complete"
done

echo ""
echo "=============================================="
echo "All figures generated"
echo "End time: $(date)"
echo "=============================================="
