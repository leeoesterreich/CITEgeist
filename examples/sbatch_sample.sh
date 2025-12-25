#!/bin/bash
#SBATCH --job-name=CITEgeist_sample
#SBATCH --output=benchmarking_logs/CITEgeist_sample_%A_%a.log
#SBATCH --error=benchmarking_logs/CITEgeist_sample_%A_%a.log
#SBATCH --time=72:00:00
#SBATCH --mem=64G
##SBATCH --mail-type=ALL
##SBATCH --mail-user=alc376@pitt.edu
#SBATCH --cluster=htc
#SBATCH --cpus-per-task=16
#SBATCH --partition=HTC
#SBATCH --array=1-14

# Activate conda environment
source /ix1/alee/LO_LAB/Personal/Alexander_Chang/miniconda3/bin/activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env/
echo "Activated conda environment"


# CHANGE THIS FOR REPRODUCIBILITY ON YOUR OWN COMPUTER
DATA_FOLDER="/ix1/alee/LO_LAB/General/Lab_Data/20250210_CITEGeistPublicData_GEO_Alex/processed_files"

# Profile discovery mode: "manual" or "auto"
# Usage: PROFILE_MODE=auto sbatch examples/sbatch_sample.sh
PROFILE_MODE="${PROFILE_MODE:-manual}"

# =============================================================================
# AUTO-PROFILE DISCOVERY PARAMETERS
# Uses spatial colocalization pipeline (Modules 1 + 2a + 2b + 2c)
#
# IMPORTANT: Modules 1-2c use RAW antibody data (before CLR transformation)
# This is handled automatically by compute_sample.py
#
# Module 1: Marker Interest Detection (adaptive GMM thresholds)
# Module 2a: Marker Colocalization (bivariate Moran's I)
# Module 2b: Profile Discovery (FDR + hierarchical clustering)
# Module 2c: Profile Selection (spatial variance optimization)
# =============================================================================
TOP_K=3                    # Mutual top-k for profile pairing (default: 3, tested value)
DISCOVERY_SEED=1234        # Random seed for reproducibility
N_PERMUTATIONS=999         # Permutations for significance testing (999 for stable p-values)
FDR_THRESHOLD=0.05         # FDR threshold for significant markers/pairs
VARIANCE_TARGET=0.90       # Target fraction of spatial variance to explain
MIN_MARGINAL_GAIN=0.005    # Minimum marginal variance gain (0.5%) - lower = more profiles

# Load Gurobi module
module load gurobi/12.0.3
echo "Loaded gurobi module"

# Change to working directory
cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist
echo "Changed to working directory"

# Get the path for this array job
path_to_visium_folder=$(sed -n "${SLURM_ARRAY_TASK_ID}p" examples/sample_paths.txt)


echo "Processing sample: $path_to_visium_folder"
echo "Profile mode: $PROFILE_MODE"

# Build command
CMD="python examples/compute_sample.py --path \"$DATA_FOLDER/$path_to_visium_folder\" --output_folder output --radius 400 --min_counts 100"

# Add auto-profiles flags if in auto mode
if [ "$PROFILE_MODE" == "auto" ]; then
    CMD="$CMD --auto-profiles"
    CMD="$CMD --top-k $TOP_K"
    CMD="$CMD --discovery-seed $DISCOVERY_SEED"
    CMD="$CMD --n-permutations $N_PERMUTATIONS"
    CMD="$CMD --fdr-threshold $FDR_THRESHOLD"
    # Module 2c: Spatial variance-based profile selection
    CMD="$CMD --variance-target $VARIANCE_TARGET"
    CMD="$CMD --min-marginal-gain $MIN_MARGINAL_GAIN"
    echo "Using auto profile discovery (spatial colocalization pipeline)"
    echo "  - top_k=$TOP_K"
    echo "  - discovery_seed=$DISCOVERY_SEED"
    echo "  - n_permutations=$N_PERMUTATIONS"
    echo "  - fdr_threshold=$FDR_THRESHOLD"
    echo "  - variance_target=$VARIANCE_TARGET"
    echo "  - min_marginal_gain=$MIN_MARGINAL_GAIN"
else
    echo "Using manual cell profiles"
fi

echo "Executing: $CMD"
eval $CMD
