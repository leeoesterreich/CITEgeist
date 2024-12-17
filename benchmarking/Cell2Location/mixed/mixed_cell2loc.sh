#!/bin/bash
#SBATCH --job-name=mx_cell2location  # Job name
#SBATCH --output=logs/%x_%a.out   # Standard output (%x = job name, %j = job ID)
#SBATCH --error=logs/%x_%a.err    # Standard error (%x = job name, %j = job ID)
#SBATCH --ntasks=1                # Number of tasks (1 per job)
#SBATCH --cpus-per-task=2         # Number of CPU cores per task
#SBATCH --mem=32G                 # Memory per node
#SBATCH --time=24:00:00           # Time limit (hh:mm:ss)
#SBATCH -M gpu                    # Partition/queue name
#SBATCH --gres=gpu:1              # Request one GPU
#SBATCH --array=0-4               # Array of jobs (one per replicate)

# load Miniconda environment
export PATH=/bgfs/alee/LO_LAB/Personal/Brent_Schlegel/miniconda3/bin:$PATH
source activate /bgfs/alee/LO_LAB/Personal/Brent_Schlegel/envs/CITEgeist

cd /bgfs/alee/LO_LAB/Personal/Brent_Schlegel/Projects/Wu_Visium/Simulations/scCube_12k/cell2location/mixed

# replicate names
REPLICATES=(
    "Wu_rep_0_GEX.h5ad"
    "Wu_rep_1_GEX.h5ad"
    "Wu_rep_2_GEX.h5ad"
    "Wu_rep_3_GEX.h5ad"
    "Wu_rep_4_GEX.h5ad"
)

# run py script with the replicate name
python mixed_cell2loc.py "${REPLICATES[$SLURM_ARRAY_TASK_ID]}"

