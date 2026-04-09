#!/bin/bash
#SBATCH --job-name=sc_spatial_viz
#SBATCH --output=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/output/sc_spatial_figures/viz_%j.log
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

source activate /ihome/alee/alc376/alc376_bgfs/envs/CITEgeist_env

mkdir -p /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/output/sc_spatial_figures

# Use v3 outputs (StarDist patchwise)
export MORPHOLOGY_VERSION=v3
python /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/examples/visualize_sc_spatial.py
