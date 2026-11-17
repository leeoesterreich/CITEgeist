#!/usr/bin/env Rscript

.cg_paths_r <- local({
  cand <- Sys.getenv("CITEGEIST_PATHS_R", unset = "")
  if (!nzchar(cand)) {
    f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
    if (length(f)) cand <- file.path(dirname(normalizePath(f[1])), "..", "..", "..", "_paths.R")
  }
  if (!nzchar(cand) || !file.exists(cand)) {
    d <- getwd()
    repeat {
      cand <- file.path(d, "Benchmarking", "_paths.R")
      if (file.exists(cand)) break
      if (dirname(d) == d) {
        stop("cannot find Benchmarking/_paths.R; set CITEGEIST_PATHS_R", call. = FALSE)
      }
      d <- dirname(d)
    }
  }
  cand
})
source(.cg_paths_r)
rm(.cg_paths_r)
library(optparse)
library(Seurat)
library(SeuratDisk)
library(spacexr)


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
sc_ref_path <- file.path(cg_root("REPO"), "Benchmarking", "simulation_benchmarking", "models", "Wu_scRNA_ref_ERpos.Rds")
sc_ref <- readRDS(sc_ref_path)
sc_data <- sc_ref@assays$RNA@counts  # Already raw integer counts
sc_meta <- sc_ref@meta.data
cell_types <- as.factor(sc_meta$celltype_major)
names(cell_types) <- sc_meta$X
reference <- Reference(sc_data, cell_types, require_int = FALSE)

# Define replicate name
replicate_name <- paste0("Wu_ST_", rep_index)

# Define file paths for current replicate
h5ad_dir <- file.path(cg_root("REPO"), "Benchmarking", "simulation_benchmarking", "replicates", "mixed", "h5ad_objects")
rep_name <- paste0("Wu_rep_", rep_index)
st_data_path <- file.path(h5ad_dir, paste0(rep_name, "_GEX_counts.csv"))
st_meta_path <- paste0(file.path(cg_root("REPO"), "Benchmarking", "simulation_benchmarking", "replicates", "mixed", "ST_sim"), "/")
st_meta_path <- paste0(st_meta_path, replicate_name, "_meta.csv")

# Load spatial transcriptomics data (corrected pseudo-counts from rebuild_sim_h5ad.py)
st_data <- read.csv(st_data_path, row.names = 1, check.names = FALSE)
st_meta <- read.csv(st_meta_path)

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

