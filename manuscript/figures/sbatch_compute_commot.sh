#!/bin/bash
#SBATCH --job-name=commot_fig4
#SBATCH --output=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/manuscript/figures/slurm_log/commot_fig4_%j.out
#SBATCH --error=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/manuscript/figures/slurm_log/commot_fig4_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

# Compute COMMOT MDK-SDC4 sender scores for Figure 4E
# Uses the COMMOT environment (has compatible numpy)

set -e

REPO="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist"
FIG_DIR="${REPO}/manuscript/figures"
COMMOT_ENV="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/COMMOT"

export PATH="${COMMOT_ENV}/bin:$PATH"
source activate "${COMMOT_ENV}"

cd "${REPO}"

echo "=============================================="
echo "Computing COMMOT Sender Scores for Figure 4E"
echo "Environment: ${COMMOT_ENV}"
echo "Start time: $(date)"
echo "=============================================="

python "${FIG_DIR}/compute_commot_scores.py"

echo ""
echo "=============================================="
echo "COMMOT computation complete"
echo "End time: $(date)"
echo "=============================================="
