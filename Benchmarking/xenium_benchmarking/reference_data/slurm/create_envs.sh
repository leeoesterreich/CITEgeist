#!/bin/bash
#SBATCH --job-name=create_envs
#SBATCH --output=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/reference_data/slurm/slurm_log/create_envs_%j.out
#SBATCH --error=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/reference_data/slurm/slurm_log/create_envs_%j.err
#SBATCH --time=04:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --cluster=htc
#SBATCH --partition=htc

# Create conda environments for Xenium benchmarking deconvolution methods
# Uses Pitt CRC approach: -p flag for prefix path, source activate for activation
# Reference: https://crc-pages.pitt.edu/user-manual/applications/python/

set -e

# Set up directories
YAML_DIR="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/reference_data/envs"
ENV_BASE="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs"

# Create env base directory
mkdir -p "${ENV_BASE}"

echo "=============================================="
echo "Creating conda environments for deconvolution"
echo "=============================================="
echo "Start time: $(date)"
echo "YAML directory: ${YAML_DIR}"
echo "Environment base: ${ENV_BASE}"
echo ""

# Initialize conda from user's miniconda3 installation
# (module load anaconda/2020.11 is no longer available on this cluster)
MINICONDA_PATH="/ix1/alee/LO_LAB/Personal/Alexander_Chang/miniconda3"
export PATH="${MINICONDA_PATH}/bin:$PATH"

# Initialize conda for shell
eval "$(${MINICONDA_PATH}/bin/conda shell.bash hook)"

echo "Conda version: $(conda --version)"
echo "Conda path: $(which conda)"
echo ""

# Function to create environment with -p prefix
create_env_with_prefix() {
    local env_name=$1
    local yaml_file=$2
    local env_path="${ENV_BASE}/${env_name}"

    echo "----------------------------------------------"
    echo "Creating environment: ${env_name}"
    echo "Path: ${env_path}"
    echo "YAML file: ${yaml_file}"
    echo "----------------------------------------------"

    # Remove existing environment if it exists
    if [ -d "${env_path}" ]; then
        echo "Environment directory exists. Removing..."
        rm -rf "${env_path}"
    fi

    # Create environment with -p prefix
    echo "Creating ${env_name} at ${env_path}..."
    conda env create -f "${yaml_file}" -p "${env_path}"

    if [ $? -eq 0 ]; then
        echo "SUCCESS: ${env_name} created at ${env_path}"
    else
        echo "ERROR: Failed to create ${env_name}"
        return 1
    fi
    echo ""
}

# Create Cell2Location environment
create_env_with_prefix "cell2location_env" "${YAML_DIR}/cell2location_env.yml"

# Create Tangram environment
create_env_with_prefix "tangram_env" "${YAML_DIR}/tangram_env.yml"

# Create R deconvolution environment (RCTD + Seurat)
create_env_with_prefix "r_deconv_env" "${YAML_DIR}/r_deconv_env.yml"

echo "=============================================="
echo "Installing spacexr (RCTD) package in R env"
echo "=============================================="

# Activate R environment using source activate with -p path
source activate "${ENV_BASE}/r_deconv_env"

# Install spacexr (RCTD) from GitHub
Rscript -e "
if (!requireNamespace('remotes', quietly = TRUE)) {
    install.packages('remotes', repos='https://cloud.r-project.org')
}
remotes::install_github('dmcable/spacexr', build_vignettes=FALSE, upgrade='never')
"

# Install SeuratDisk for h5seurat support
Rscript -e "
if (!requireNamespace('remotes', quietly = TRUE)) {
    install.packages('remotes', repos='https://cloud.r-project.org')
}
remotes::install_github('mojaveazure/seurat-disk', upgrade='never')
"

conda deactivate

echo ""
echo "=============================================="
echo "Environment creation complete!"
echo "=============================================="
echo "End time: $(date)"
echo ""
echo "Created environments at ${ENV_BASE}/:"
ls -la "${ENV_BASE}/"
echo ""
echo "To activate (use source activate with -p path):"
echo "  source activate ${ENV_BASE}/cell2location_env"
echo "  source activate ${ENV_BASE}/tangram_env"
echo "  source activate ${ENV_BASE}/r_deconv_env"
