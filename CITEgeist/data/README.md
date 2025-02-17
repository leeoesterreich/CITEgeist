# Reproducibility Guide

## For Reviewers

### 1. Download the Code

- Download the code from Figshare (see the Manuscript Submission for private access)
- Unzip the downloaded file to your preferred location

### 2. Download and Prepare the Data

1. Download the data from GEO (see Data Availability section in the manuscript for the GEO link and Private Access Token)
2. Run the following script to reorganize the data into SpaceRanger output folders, and strip unique identifiers required by GEO:

```bash
# untar the raw files
mkdir -pv ./GEO_data
tar -xvf GEO_data_RAW.tar -C ./GEO_data

# run the py preprocessing script
python3 ./organize_spatial_data.py --folder GEO_data
```

*Note*: When prompted, type 'yes' to confirm.

### 3. Set Up the Environment

1. Install dependencies using the provided environment file:
```bash
conda env create -f CITEgeist_env.yml
```

2. Activate the environment and set up a Jupyter kernel
```bash
conda activate CITEgeist_env
```

### 4. Obtain Gurobi License

CITEgeist requires a Gurobi license (free for academic use):

1. Sign up for a Gurobi [academic license](https://www.gurobi.com/downloads/end-user-license-agreement-academic/)
2. Follow the instructions to download and install your license
3. Update the license file path in the notebooks to match your local license location

### 5. Run the Analysis ([See the "Examples" folder for details](../examples/README.md))

You can either:

A. Run the notebooks directly:
- Update data paths in the notebooks to match your local directory structure
- Expected runtime: ~2 hours on a standard computer (16 threads, 32GB RAM)

B. Use SLURM distribution:
- Use the provided `examples/sbatch_sample.sh` script for distributed computing

## System Requirements

- RAM: 32GB (minimum)
- CPU: 16 threads (recommended)
- Storage: Sufficient space for the GEO dataset
- Operating System: Linux/Unix recommended (Windows users may need additional configuration)
