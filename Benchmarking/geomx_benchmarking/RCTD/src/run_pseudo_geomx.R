#!/usr/bin/env Rscript
# Run RCTD deconvolution on pseudo-GeoMx AOI data.

library(optparse)
library(spacexr)
library(Matrix)
library(data.table)

option_list <- list(
  make_option(c("--ref-counts"), type="character",
              default="reference_data/processed/rctd/counts.csv"),
  make_option(c("--ref-cell-types"), type="character",
              default="reference_data/processed/rctd/cell_types.csv"),
  make_option(c("--spatial-counts"), type="character",
              default="data/aoi_counts.csv"),
  make_option(c("--spatial-coords"), type="character",
              default="data/aoi_coords.csv"),
  make_option(c("--output-dir"), type="character",
              default="RCTD/output"),
  make_option(c("--max-cores"), type="integer", default=8)
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

cat("============================================\n")
cat("RCTD Deconvolution - Pseudo-GeoMx\n")
cat("============================================\n")
cat(sprintf("Start time: %s\n", Sys.time()))

base_dir <- getwd()
resolve <- function(p) normalizePath(file.path(base_dir, p), mustWork=TRUE)

ref_counts_path <- resolve(opt$`ref-counts`)
ref_ct_path <- resolve(opt$`ref-cell-types`)
spatial_counts_path <- resolve(opt$`spatial-counts`)
spatial_coords_path <- resolve(opt$`spatial-coords`)
output_dir <- file.path(base_dir, opt$`output-dir`)

dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)

cat("\nLoading reference data...\n")
ref_dt <- fread(ref_counts_path, header=TRUE)
rn <- ref_dt[[1]]
ref_dt[, 1 := NULL]
ref_counts <- as.matrix(ref_dt)
rownames(ref_counts) <- rn
ref_cell_types <- read.csv(ref_ct_path, row.names=1, check.names=FALSE)
cat(sprintf("  Reference: %d genes x %d cells\n", nrow(ref_counts), ncol(ref_counts)))

cell_types <- factor(ref_cell_types$cell_type)
names(cell_types) <- rownames(ref_cell_types)
reference <- Reference(counts=ref_counts, cell_types=cell_types, require_int=FALSE)

cat("\nLoading spatial data...\n")
spatial_counts <- as.matrix(read.csv(spatial_counts_path, row.names=1, check.names=FALSE))
spatial_coords <- read.csv(spatial_coords_path, row.names=1, check.names=FALSE)
cat(sprintf("  Spatial: %d genes x %d AOIs\n", nrow(spatial_counts), ncol(spatial_counts)))

coords <- as.data.frame(spatial_coords)
colnames(coords) <- c("x", "y")
puck <- SpatialRNA(coords=coords, counts=spatial_counts, require_int=FALSE)

cat("\nRunning RCTD...\n")
start_time <- Sys.time()
myRCTD <- create.RCTD(puck, reference, max_cores=opt$`max-cores`, CELL_MIN_INSTANCE=5)
myRCTD@config$max_cores <- 1
myRCTD <- run.RCTD(myRCTD, doublet_mode="full")
elapsed <- as.numeric(difftime(Sys.time(), start_time, units="mins"))
cat(sprintf("RCTD completed in %.1f minutes\n", elapsed))

cat("\nExtracting proportions...\n")
# In full mode, weights are a matrix at @results$weights (not per-pixel S4 objects)
weights <- as.data.frame(as.matrix(myRCTD@results$weights))
weights <- sweep(weights, 1, rowSums(weights), "/")

prop_path <- file.path(output_dir, "predicted_proportions.csv")
write.csv(weights, prop_path)
cat(sprintf("Wrote: %s (%d x %d)\n", prop_path, nrow(weights), ncol(weights)))

saveRDS(myRCTD, file.path(output_dir, "rctd_object.rds"))
cat(sprintf("\nDone: %s\n", Sys.time()))
