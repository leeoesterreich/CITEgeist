#!/usr/bin/env Rscript

library(optparse)
library(Seurat)
library(SeuratDisk)
library(tidyverse)
library(spacexr)
library(SingleCellExperiment)

# Command-line arguments parsing
option_list <- list(
  make_option(c("-r", "--replicates"), type = "character", default = "0",
              help = "Single replicate index to process", metavar = "character"),
  make_option(c("-o", "--output_dir"), type = "character", default = ".",
              help = "Output directory for results", metavar = "character")
)

opt <- parse_args(OptionParser(option_list = option_list))

# Parse replicate index
rep_index <- as.numeric(opt$replicates)
output_dir <- opt$output_dir

# Load and preprocess reference data (done only once)
cat("Loading and preprocessing single-cell reference data...\n")
sc_ref_path <- "/bgfs/alee/LO_LAB/Personal/Brent_Schlegel/Projects/Wu_Visium/Simulations/models/Wu_scRNA_ref_ERpos.h5seurat"
sc_ref <- LoadH5Seurat(sc_ref_path)
sc_data <- sc_ref@assays$RNA@counts
sc_data <- round(exp(sc_data) - 1)  # Convert logcounts to counts
sc_meta <- sc_ref@meta.data
cell_types <- as.factor(sc_meta$celltype_major)
names(cell_types) <- sc_meta$X
reference <- Reference(sc_data, cell_types, require_int = FALSE)

# Define replicate name
replicate_name <- paste0("Wu_ST_", rep_index)

# Define file paths for current replicate
st_data_path <- paste0("/bgfs/alee/LO_LAB/Personal/Brent_Schlegel/Projects/Wu_Visium/Simulations/scCube_12k/replicates/mixed/ST_sim/", replicate_name, "_data.csv")
st_meta_path <- paste0("/bgfs/alee/LO_LAB/Personal/Brent_Schlegel/Projects/Wu_Visium/Simulations/scCube_12k/replicates/mixed/ST_sim/", replicate_name, "_meta.csv")

# Load spatial transcriptomics data
st_data <- read.csv(st_data_path, row.names = 1)
st_meta <- read.csv(st_meta_path)
st_data <- round(st_data)  # Ensure counts are integers

# Extract coordinates
coords <- data.frame(x = st_meta$spot_x, y = st_meta$spot_y)
rownames(coords) <- colnames(st_data)

# Create SpatialRNA object
puck <- SpatialRNA(coords, st_data, require_int = FALSE)

# Run RCTD pipeline
myRCTD <- create.RCTD(puck, reference)
myRCTD <- run.RCTD(myRCTD, doublet_mode = 'full')
results <- myRCTD@results

# Normalize weights to sum to 1
norm_weights <- normalize_weights(results$weights)
norm_weights <- as.matrix(norm_weights)

# Define output path
output_file <- paste0(output_dir, "/", replicate_name, "_RCTD_deconv_predictions.csv")

# Save normalized weights to CSV
write.csv(norm_weights, output_file)

# Print progress
cat("Completed RCTD pipeline for replicate:", replicate_name, "\n")

