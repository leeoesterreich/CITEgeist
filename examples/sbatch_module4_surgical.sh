#!/bin/bash
#SBATCH --job-name=Module4_Surgical
#SBATCH --output=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/output/logs/module4_surgical_%A_%a.log
#SBATCH --error=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/output/logs/module4_surgical_%A_%a.err
#SBATCH --time=4:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --array=0-4
#SBATCH --partition=htc

# Surgical samples that need Module 4 analysis
SAMPLES=(
    "HCC22-088-P1-S2"
    "HCC22-088-P2-S2"
    "HCC22-088-P3-S2"
    "HCC22-088-P4-S2"
    "HCC22-088-P6-S2_D"
)

SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

echo "=========================================="
echo "Processing surgical sample: $SAMPLE"
echo "Task ID: $SLURM_ARRAY_TASK_ID"
echo "Date: $(date)"
echo "=========================================="

# Use direct Python path (conda activate doesn't work on compute nodes)
PYTHON=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env/bin/python
export PYTHONPATH=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist:$PYTHONPATH

# Create log directory
mkdir -p /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/output/logs

# Run Module 4 for this sample
$PYTHON << EOF
import sys
import logging
from pathlib import Path

sys.path.insert(0, "/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist")

import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq

from CITEgeist.model import (
    discover_programs_from_layers,
    analyze_program_relationships,
)

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

SAMPLE = "$SAMPLE"
DATA_DIR = Path("/ix1/alee/LO_LAB/General/Lab_Data/20250210_CITEGeistPublicData_GEO_Alex/processed_files")
OUTPUT_DIR = Path("/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/output")

logger.info(f"Processing {SAMPLE}")

# Check if already done
m4_file = OUTPUT_DIR / f"{SAMPLE}_module4_programs.csv"
if m4_file.exists():
    logger.info(f"Module 4 output already exists for {SAMPLE}, skipping")
    sys.exit(0)

# Load data
sample_dir = DATA_DIR / SAMPLE / "outs"
if not sample_dir.exists():
    logger.error(f"Sample directory not found: {sample_dir}")
    sys.exit(1)

logger.info(f"Loading from {sample_dir}")
adata = sq.read.visium(
    str(sample_dir),
    counts_file='filtered_feature_bc_matrix.h5',
    load_images=True,
    gex_only=False,
)

# Split GEX and antibody
feature_types = adata.var['feature_types']
gex_mask = feature_types == 'Gene Expression'
ab_mask = feature_types == 'Antibody Capture'

adata_gex = adata[:, gex_mask].copy()
adata_cite = adata[:, ab_mask].copy()
adata_gex.var_names_make_unique()

logger.info(f"GEX: {adata_gex.shape}, Antibody: {adata_cite.shape}")

# Load deconvolved layers if available
layers_dir = OUTPUT_DIR / f"{SAMPLE}_pass1" / "layers" / "pass1"
n_layers = 0

if layers_dir.exists():
    layer_files = list(layers_dir.glob("*_layer_pass1.csv"))

    for layer_file in layer_files:
        cell_type = layer_file.stem.replace("_layer_pass1", "")
        df = pd.read_csv(layer_file, index_col=0)

        common_spots = list(adata_gex.obs_names.intersection(df.index))
        if len(common_spots) == 0:
            continue

        layer_matrix = np.zeros((adata_gex.shape[0], adata_gex.shape[1]))
        adata_gene_to_idx = {g: i for i, g in enumerate(adata_gex.var_names)}
        common_genes = [g for g in df.columns if g in adata_gene_to_idx]

        if len(common_genes) == 0:
            continue

        spot_idx = np.array([adata_gex.obs_names.get_loc(s) for s in common_spots])
        gene_idx = np.array([adata_gene_to_idx[g] for g in common_genes])

        layer_data = df.loc[common_spots, common_genes].values
        for i, s_idx in enumerate(spot_idx):
            layer_matrix[s_idx, gene_idx] = layer_data[i, :]

        layer_name = f"{cell_type}_genes_pass1"
        adata_gex.layers[layer_name] = layer_matrix
        n_layers += 1

logger.info(f"Loaded {n_layers} deconvolved layers")

if n_layers == 0:
    logger.warning("No deconvolved layers found - using raw expression with cell type proportions")

    # Load cell type proportions
    prop_file = OUTPUT_DIR / f"{SAMPLE}_cell_prop_finetuned_results.csv"
    if not prop_file.exists():
        logger.error(f"No proportions file found: {prop_file}")
        sys.exit(1)

    props = pd.read_csv(prop_file, index_col=0)

    # Create mock layers from proportions
    X_raw = adata_gex.X.toarray() if hasattr(adata_gex.X, 'toarray') else np.array(adata_gex.X)

    # Align proportions to adata
    common = adata_gex.obs_names.intersection(props.index)
    props_aligned = pd.DataFrame(0.0, index=adata_gex.obs_names, columns=props.columns)
    props_aligned.loc[common, :] = props.loc[common, :]

    for ct in props.columns:
        ct_props = props_aligned[ct].values.reshape(-1, 1)
        layer_data = X_raw * ct_props
        layer_name = f"{ct}_genes_pass1"
        adata_gex.layers[layer_name] = layer_data
        n_layers += 1

    logger.info(f"Created {n_layers} mock deconvolved layers from proportions")

# Load cell type proportions for Module 4
prop_file = OUTPUT_DIR / f"{SAMPLE}_cell_prop_finetuned_results.csv"
cell_type_proportions = None
if prop_file.exists():
    props_raw = pd.read_csv(prop_file, index_col=0)
    cell_type_proportions = pd.DataFrame(0.0, index=adata_gex.obs_names, columns=props_raw.columns)
    common = adata_gex.obs_names.intersection(props_raw.index)
    cell_type_proportions.loc[common, :] = props_raw.loc[common, :]
    logger.info(f"Loaded proportions: {cell_type_proportions.shape}")

# Run Module 4
logger.info("Running Module 4: Program Discovery")
result = discover_programs_from_layers(
    adata=adata_gex,
    layer_pattern="_genes_pass1",
    cell_type_proportions=cell_type_proportions,
    protein_adata=adata_cite,
    K_programs=3,
    lambda_spatial=0.1,
    lambda_sparsity=0.01,
    validate_with_proteins=True,
    detect_subpopulations=True,
    n_subpopulations=3,
    random_state=42,
)

# Save Module 4 results
result.to_dataframe().to_csv(m4_file)
logger.info(f"Saved Module 4 results to {m4_file}")

# Run Module 4b
logger.info("Running Module 4b: Bivariate Relationships")
m4b_result = analyze_program_relationships(
    result=result,
    adata=adata_gex,
    neighbor_k=6,
    n_permutations=199,
    fdr_threshold=0.1,
)

m4b_file = OUTPUT_DIR / f"{SAMPLE}_module4b_relationships.csv"
m4b_result.to_dataframe().to_csv(m4b_file)
logger.info(f"Saved Module 4b results to {m4b_file}")

logger.info(f"Complete! {SAMPLE} processed successfully")

EOF

echo "Job complete for $SAMPLE"
