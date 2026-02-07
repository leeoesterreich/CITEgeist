#!/bin/bash
#SBATCH --job-name=cg_no_unk
#SBATCH --output=slurm_log/cg_no_unk_%j.out
#SBATCH --error=slurm_log/cg_no_unk_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alc376@pitt.edu

# ============================================================================
# CITEgeist Benchmark: Achievable-7 WITHOUT Unknown profile
# ============================================================================
# Region 0 only (for GEX comparison figure)

REPO_ROOT="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist"
PSEUDOVISIUM_DIR="${REPO_ROOT}/Benchmarking/xenium_pseudovisium"
BENCH_DIR="${REPO_ROOT}/Benchmarking/xenium_benchmarking"
CITEGEIST_DIR="${BENCH_DIR}/CITEgeist"
INPUT_DIR="${PSEUDOVISIUM_DIR}/data_granular_gt"
OUTPUT_DIR="${CITEGEIST_DIR}/output_achievable_7_no_unknown"

mkdir -p "${OUTPUT_DIR}"

# Environment
export PATH="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env/bin:$PATH"
source activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env
module load gurobi/12.0.3

cd "${REPO_ROOT}"

echo "Running CITEgeist region 0: no Unknown, radius=150, NO Laplacian smoothing..."
python "${CITEGEIST_DIR}/src/run_benchmark.py" \
    --region-id 0 \
    --mode manual \
    --input-dir "${INPUT_DIR}" \
    --output-dir "${OUTPUT_DIR}" \
    --radius 250.0 \
    --lambda-reg 1.0 \
    --alpha-elastic 0.7 \
    --max-y-change 0.4 \
    --min-counts 25 \
    --lambda-laplacian 0.0 \
    --lambda-reg-gex 0.01 \
    --run-gex \
    --no-unknown

echo "Done!"
