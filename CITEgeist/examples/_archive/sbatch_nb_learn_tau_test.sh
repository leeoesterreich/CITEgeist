#!/bin/bash
#SBATCH --job-name=nb_tau_test
#SBATCH --time=02:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --output=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/CITEgeist/examples/logs/nb_tau_test_%j.out
#SBATCH --error=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/CITEgeist/examples/logs/nb_tau_test_%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

set -eo pipefail

module load gurobi/12.0.3
eval "$(conda shell.bash hook 2>/dev/null)"
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

set -u

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/CITEgeist

python examples/run_nb_learn_tau_test.py
