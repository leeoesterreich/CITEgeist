#!/usr/bin/env Rscript
#
# Run CARD deconvolution on mixed simulation data.
#
# CARD incorporates spatial correlation in cell type composition estimation.
# This script supports both reference and reference-free modes.
#
# Usage:
#   Rscript CARD_pipeline_mixed.R --replicates 0 --output_dir . --mode reference

library(optparse)
library(CARD)
library(Matrix)

# Command-line arguments
option_list <- list(
  make_option(c("-r", "--replicates"), type="character", default="0",
              help="Single replicate index to process", metavar="character"),
  make_option(c("-o", "--output_dir"), type="character", default=".",
              help="Output directory for results", metavar="character"),
  make_option(c("-m", "--mode"), type="character", default="reference",
              help="Mode: 'reference' or 'reference_free'", metavar="character")
)

opt <- parse_args(OptionParser(option_list=option_list))

rep_index <- as.numeric(opt$replicates)
output_dir <- opt$output_dir
mode <- opt$mode

cat("============================================\n")
cat(sprintf("CARD Deconvolution - mixed Replicate %d (%s mode)\n", rep_index, mode))
cat("============================================\n")
cat(sprintf("Start time: %s\n", Sys.time()))

# Base paths
BASE_PATH <- "/ix1/alee/LO_LAB/Personal/Brent_Schlegel/bts76/Projects/CITEgeist/Wu_Visium/Simulations"
CARD_DIR <- "/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/simulation_benchmarking/CARD"
REF_COUNTS_PATH <- file.path(CARD_DIR, "reference_csv/reference_counts.csv")
REF_CELLTYPES_PATH <- file.path(CARD_DIR, "reference_csv/reference_cell_types.csv")
DATA_PATH <- file.path(BASE_PATH, "scCube_12k/replicates/mixed/ST_sim")

replicate_name <- paste0("Wu_ST_", rep_index)

# Load spatial data
cat("\nLoading spatial data...\n")
st_data_path <- file.path(DATA_PATH, paste0(replicate_name, "_data.csv"))
st_meta_path <- file.path(DATA_PATH, paste0(replicate_name, "_meta.csv"))

st_data <- as.matrix(read.csv(st_data_path, row.names=1, check.names=FALSE))
st_meta <- read.csv(st_meta_path)

cat(sprintf("  Spatial data: %d genes x %d spots\n", nrow(st_data), ncol(st_data)))

# Prepare spatial location
spatial_location <- data.frame(
  x = st_meta$spot_x,
  y = st_meta$spot_y,
  row.names = colnames(st_data)
)

start_time <- Sys.time()

if (mode == "reference") {
  # =========================================
  # Reference-based CARD
  # =========================================
  cat("\nLoading reference data from CSV...\n")

  # Load pre-converted reference CSV files
  sc_data <- as.matrix(read.csv(REF_COUNTS_PATH, row.names=1, check.names=FALSE))
  sc_celltypes <- read.csv(REF_CELLTYPES_PATH, row.names=1, check.names=FALSE)

  cat(sprintf("  Reference: %d genes x %d cells\n", nrow(sc_data), ncol(sc_data)))
  cat(sprintf("  Cell types: %s\n", paste(unique(sc_celltypes$cell_type), collapse=", ")))

  # Create metadata for CARD
  sc_meta_card <- data.frame(
    cellID = rownames(sc_celltypes),
    cellType = sc_celltypes$cell_type,
    sampleInfo = "sample1",
    row.names = rownames(sc_celltypes)
  )

  # Create CARD object
  cat("\nCreating CARD object...\n")
  CARD_obj <- createCARDObject(
    sc_count = sc_data,
    sc_meta = sc_meta_card,
    spatial_count = st_data,
    spatial_location = spatial_location,
    ct.varname = "cellType",
    ct.select = unique(sc_meta_card$cellType),
    sample.varname = "sampleInfo",
    minCountGene = 100,
    minCountSpot = 5
  )

  # Run deconvolution
  cat("\nRunning CARD deconvolution...\n")
  CARD_obj <- CARD_deconvolution(CARD_object = CARD_obj)

  # Extract proportions
  proportions <- as.data.frame(CARD_obj@Proportion_CARD)

} else if (mode == "reference_free") {
  # =========================================
  # Reference-free CARD
  # =========================================
  cat("\nRunning CARD in reference-free mode...\n")

  # Load marker genes (need to generate these first)
  marker_file <- file.path(dirname(output_dir), "marker_genes_simulation.csv")

  if (!file.exists(marker_file)) {
    stop(sprintf("Marker genes file not found: %s\nPlease generate markers first.", marker_file))
  }

  marker_df <- read.csv(marker_file, check.names=FALSE)
  markerList <- split(marker_df$gene, marker_df$cell_type)

  cat(sprintf("  Loaded markers for %d cell types\n", length(markerList)))

  # Create CARDfree object
  cat("\nCreating CARDfree object...\n")
  CARDfree_obj <- createCARDfreeObject(
    markerList = markerList,
    spatial_count = st_data,
    spatial_location = spatial_location,
    minCountGene = 100,
    minCountSpot = 5
  )

  # Run reference-free deconvolution
  cat("\nRunning CARD reference-free deconvolution...\n")
  CARDfree_obj <- CARD_refFree(CARDfree_obj)

  # Extract proportions
  proportions <- as.data.frame(CARDfree_obj@Proportion_CARD)

} else {
  stop(sprintf("Unknown mode: %s", mode))
}

end_time <- Sys.time()
runtime <- as.numeric(difftime(end_time, start_time, units="secs"))

cat(sprintf("\nCARD completed in %.1f seconds (%.1f min)\n", runtime, runtime/60))
cat(sprintf("  Proportions: %d spots x %d cell types\n", nrow(proportions), ncol(proportions)))

# Save results
dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)

if (mode == "reference") {
  output_file <- file.path(output_dir, paste0(replicate_name, "_CARD_deconv_predictions.csv"))
} else {
  output_file <- file.path(output_dir, paste0(replicate_name, "_CARD_reffree_deconv_predictions.csv"))
}

write.csv(proportions, output_file, row.names=TRUE)
cat(sprintf("\nSaved predictions to: %s\n", output_file))

cat("\n============================================\n")
cat(sprintf("Completed CARD pipeline for replicate: %s\n", replicate_name))
cat(sprintf("End time: %s\n", Sys.time()))
