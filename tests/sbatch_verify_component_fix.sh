#!/bin/bash
#SBATCH --job-name=verify_component_fix
#SBATCH --output=slurm_log/verify_component_fix_%j.out
#SBATCH --error=slurm_log/verify_component_fix_%j.err
#SBATCH --time=1:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alc376@pitt.edu

echo "========================================"
echo "Verifying Connected Component Fix"
echo "========================================"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: ${SLURMD_NODENAME}"
echo "Start time: $(date)"
echo ""

# Load modules
module purge
module load gcc/12.2.0
module load gurobi/12.0.3
module load anaconda3/2023.09-0-python_3.11.5
source activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

PYTHON="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env/bin/python"
SCRIPT_DIR="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/.worktrees/hierarchical-profiles/tests"

${PYTHON} "${SCRIPT_DIR}/verify_component_fix.py"

exit_code=$?

echo ""
echo "========================================"
echo "Verification Complete"
echo "Exit code: ${exit_code}"
echo "End time: $(date)"
echo "========================================"

exit ${exit_code}
