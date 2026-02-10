#!/bin/bash
#SBATCH --job-name=gen_supp_figs
#SBATCH --output=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/manuscript/figures/slurm_log/gen_supp_figs_%j.out
#SBATCH --error=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/manuscript/figures/slurm_log/gen_supp_figs_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

set -e

REPO="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist"
FIG_DIR="${REPO}/manuscript/figures"

export PATH="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env/bin:$PATH"
source activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

cd "${REPO}"

export CITEGEIST_DATA_ROOT="${REPO}"
export CITEGEIST_OUTPUT_ROOT="${REPO}/output"
export CITEGEIST_LICENSE_FILE="${REPO}/LICENSE"

echo "=============================================="
echo "Generating Supplementary Figures"
echo "Start time: $(date)"
echo "=============================================="
echo "DEPRECATION NOTICE: canonical entrypoint is now"
echo "  python -m repro.cli.run_figures --set v5_supp --config repro/config/example.paths.yaml"

python -m repro.cli.run_figures --set v5_supp

echo ""
echo "=============================================="
echo "Supplementary figure generation complete"
echo "End time: $(date)"
echo "=============================================="
