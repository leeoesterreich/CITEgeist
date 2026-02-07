#!/bin/bash
#SBATCH --job-name=CITEgeist_fixed
#SBATCH --output=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/logs/citegeist_fixed_%A_%a.log
#SBATCH --error=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/logs/citegeist_fixed_%A_%a.log
#SBATCH --array=0-4
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --cluster=htc

source ~/.bashrc
module load gurobi/12.0.3

export PATH="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env/bin:$PATH"

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist

python << PYEOF
import sys
import os
import scanpy as sc
import pandas as pd
import numpy as np

sys.path.insert(0, "/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/CITEgeist")
from model.citegeist_model import CitegeistModel

region_id = int(os.environ.get("SLURM_ARRAY_TASK_ID", 0))
print(f"Processing region {region_id}")

# Paths
data_dir = "Benchmarking/xenium_benchmarking/data"
gex_path = f"{data_dir}/h5ad_objects/Xenium_region_{region_id}_GEX.h5ad"
cite_path = f"{data_dir}/h5ad_objects/Xenium_region_{region_id}_CITE.h5ad"
output_dir = f"Benchmarking/xenium_benchmarking/CITEgeist/fixed_region_{region_id}"

os.makedirs(output_dir, exist_ok=True)

print(f"Loading data...")
adata_gex = sc.read_h5ad(gex_path)
adata_cite = sc.read_h5ad(cite_path)

print(f"GEX: {adata_gex.shape}, CITE: {adata_cite.shape}")

# Cell profile dict for Xenium proteins
cell_profile_dict = {
    "CD4+ T cells": {"Major": ["CD3E", "CD4"], "Minor": ["CD45"]},
    "CD8+ T cells": {"Major": ["CD3E", "CD8A"], "Minor": ["CD45", "GranzymeB"]},
    "B cells": {"Major": ["CD20"], "Minor": ["CD45"]},
    "Plasma cells": {"Major": ["CD138"], "Minor": []},
    "Macrophages": {"Major": ["CD68"], "Minor": ["CD163", "CD45", "HLA-DR"]},
    "Dendritic cells": {"Major": ["CD11c", "HLA-DR"], "Minor": ["CD45"]},
    "NK cells": {"Major": ["CD16", "GranzymeB"], "Minor": ["CD45"]},
    "Epithelial": {"Major": ["PanCK"], "Minor": ["E-Cadherin"]},
    "Endothelial": {"Major": ["CD31"], "Minor": []},
    "Fibroblasts": {"Major": ["alphaSMA"], "Minor": ["Vimentin"]},
}

# Initialize model
model = CitegeistModel(
    sample_name=f"Xenium_region_{region_id}",
    output_folder=output_dir,
    simulation=True,
    gene_expression_adata=adata_gex,
    antibody_capture_adata=adata_cite,
)

model.load_cell_profile_dict(cell_profile_dict)
model.preprocess_gex(target_sum=1e4)
model.preprocess_antibody()

print("Running cell proportion model...")
model.run_cell_proportion_model(
    radius=4,
    lambda_reg=0.01,
    alpha=0.5,
    per_marker_beta=True,
    checkpoint_interval=50,
    skip_finetuning=False,
    validation_warn_only=True,  # Allow low-proportion cell types
)

print(f"Done with region {region_id}")
PYEOF
