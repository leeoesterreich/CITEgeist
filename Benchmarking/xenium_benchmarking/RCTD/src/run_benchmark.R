#!/usr/bin/env Rscript
#
# Run RCTD deconvolution on Xenium pseudo-Visium data.
#
# RCTD (Robust Cell Type Decomposition) uses a probabilistic model
# to estimate cell type proportions in spatial transcriptomics data.
#
# Note: RCTD only outputs cell type proportions, not per-cell-type
# gene expression layers.
#
# Usage:
#   Rscript run_benchmark.R --region-id 0 --ref-counts ref_counts.csv \
#     --ref-cell-types ref_cell_types.csv --spatial-counts spatial_counts.csv \
#     --spatial-coords spatial_coords.csv --output-dir output/

library(optparse)
library(spacexr)
library(Matrix)

# Parse command line arguments
option_list <- list(
  make_option(c("--region-id"), type="integer", default=0,
              help="Region ID (0-4)"),
  make_option(c("--ref-counts"), type="character", default=NULL,
              help="Reference counts CSV (genes x cells)"),
  make_option(c("--ref-cell-types"), type="character", default=NULL,
              help="Reference cell types CSV"),
  make_option(c("--spatial-counts"), type="character", default=NULL,
              help="Spatial counts CSV (genes x spots)"),
  make_option(c("--spatial-coords"), type="character", default=NULL,
              help="Spatial coordinates CSV"),
  make_option(c("--output-dir"), type="character", default="output_rna_gt",
              help="Output directory"),
  make_option(c("--prefix"), type="character", default="Xenium",
              help="Filename prefix"),
  make_option(c("--max-cores"), type="integer", default=8,
              help="Maximum cores for parallel processing")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

cat("============================================\n")
cat(sprintf("RCTD Deconvolution - Region %d\n", opt$`region-id`))
cat("============================================\n")
cat(sprintf("Start time: %s\n", Sys.time()))
cat(sprintf("Reference counts: %s\n", opt$`ref-counts`))
cat(sprintf("Reference cell types: %s\n", opt$`ref-cell-types`))
cat(sprintf("Spatial counts: %s\n", opt$`spatial-counts`))
cat(sprintf("Spatial coords: %s\n", opt$`spatial-coords`))
cat(sprintf("Output directory: %s\n", opt$`output-dir`))
cat("\n")

# Load reference data
cat("Loading reference data...\n")
ref_counts <- as.matrix(read.csv(opt$`ref-counts`, row.names=1, check.names=FALSE))
ref_cell_types <- read.csv(opt$`ref-cell-types`, row.names=1, check.names=FALSE)

cat(sprintf("  Reference counts: %d genes x %d cells\n", nrow(ref_counts), ncol(ref_counts)))
cat(sprintf("  Cell types: %s\n", paste(unique(ref_cell_types$cell_type), collapse=", ")))

# Create Reference object
cat("Creating Reference object...\n")
cell_types <- factor(ref_cell_types$cell_type)
names(cell_types) <- rownames(ref_cell_types)

reference <- Reference(
  counts = ref_counts,
  cell_types = cell_types,
  require_int = FALSE  # Allow non-integer counts
)

# Load spatial data
cat("\nLoading spatial data...\n")
spatial_counts <- as.matrix(read.csv(opt$`spatial-counts`, row.names=1, check.names=FALSE))
spatial_coords <- read.csv(opt$`spatial-coords`, row.names=1, check.names=FALSE)

cat(sprintf("  Spatial counts: %d genes x %d spots\n", nrow(spatial_counts), ncol(spatial_counts)))
cat(sprintf("  Coordinates: %d spots\n", nrow(spatial_coords)))

# Create SpatialRNA object (puck)
cat("Creating SpatialRNA object...\n")
coords <- as.data.frame(spatial_coords)
colnames(coords) <- c("x", "y")

puck <- SpatialRNA(
  coords = coords,
  counts = spatial_counts,
  require_int = FALSE
)

# Run RCTD
cat("\nRunning RCTD deconvolution...\n")
start_time <- Sys.time()

myRCTD <- create.RCTD(
  puck,
  reference,
  max_cores = opt$`max-cores`,
  CELL_MIN_INSTANCE = 5  # Minimum cells per cell type
)

myRCTD <- run.RCTD(
  myRCTD,
  doublet_mode = "full"  # Full mode for mixed cell types
)

end_time <- Sys.time()
runtime <- as.numeric(difftime(end_time, start_time, units="secs"))
cat(sprintf("RCTD completed in %.1f seconds (%.1f min)\n", runtime, runtime/60))

# Extract proportions
cat("\nExtracting cell type proportions...\n")
results <- myRCTD@results

# Normalize weights using spacexr function and convert to matrix
norm_weights <- normalize_weights(results$weights)
proportions <- as.data.frame(as.matrix(norm_weights))

cat(sprintf("  Proportions shape: %d spots x %d cell types\n",
            nrow(proportions), ncol(proportions)))

# Create output directory
sample_name <- sprintf("%s_region_%d", opt$prefix, opt$`region-id`)
output_dir <- file.path(opt$`output-dir`, sample_name)
dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)

# Save proportions (using standardized output filename)
output_file <- file.path(output_dir, sprintf("%s_cell_prop_finetuned_results.csv", sample_name))
write.csv(proportions, output_file, row.names=TRUE)
cat(sprintf("Saved proportions to: %s\n", output_file))

# Save run statistics
stats <- list(
  region_id = opt$`region-id`,
  sample_name = sample_name,
  n_spots = nrow(proportions),
  n_cell_types = ncol(proportions),
  cell_types = colnames(proportions),
  runtime_sec = runtime,
  method = "RCTD",
  doublet_mode = "full"
)

stats_file <- file.path(output_dir, "run_stats.json")
cat(jsonlite::toJSON(stats, auto_unbox=TRUE, pretty=TRUE), file=stats_file)
cat(sprintf("Saved stats to: %s\n", stats_file))

# Save RCTD object for further analysis
rctd_file <- file.path(output_dir, sprintf("%s_RCTD.rds", sample_name))
saveRDS(myRCTD, rctd_file)
cat(sprintf("Saved RCTD object to: %s\n", rctd_file))

cat("\n============================================\n")
cat("Deconvolution Complete!\n")
cat("============================================\n")
cat(sprintf("  Spots: %d\n", nrow(proportions)))
cat(sprintf("  Cell types: %d\n", ncol(proportions)))
cat(sprintf("  Runtime: %.1f seconds\n", runtime))
cat(sprintf("  Output: %s\n", output_dir))
cat(sprintf("End time: %s\n", Sys.time()))
