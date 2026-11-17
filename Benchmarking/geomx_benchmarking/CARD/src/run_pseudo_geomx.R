#!/usr/bin/env Rscript
# Run CARD deconvolution on pseudo-GeoMx AOI data.

library(optparse)
library(CARD)
library(Matrix)
library(data.table)

option_list <- list(
  make_option(c("--ref-counts"), type="character",
              default="reference_data/processed/rctd/counts.csv"),
  make_option(c("--ref-cell-types"), type="character",
              default="reference_data/processed/rctd/cell_types_with_sample.csv"),
  make_option(c("--spatial-counts"), type="character",
              default="data/aoi_counts.csv"),
  make_option(c("--spatial-coords"), type="character",
              default="data/aoi_coords.csv"),
  make_option(c("--output-dir"), type="character",
              default="CARD/output")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

cat("============================================\n")
cat("CARD Deconvolution - Pseudo-GeoMx\n")
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
ref_meta <- read.csv(ref_ct_path, row.names=1, check.names=FALSE)
dup_genes <- duplicated(rownames(ref_counts))
if (any(dup_genes)) {
    cat(sprintf("  Removing %d duplicate gene names from reference\n", sum(dup_genes)))
    ref_counts <- ref_counts[!dup_genes, ]
}
cat(sprintf("  Reference: %d genes x %d cells\n", nrow(ref_counts), ncol(ref_counts)))

cat("\nLoading spatial data...\n")
spatial_counts <- as.matrix(read.csv(spatial_counts_path, row.names=1, check.names=FALSE))
spatial_coords <- read.csv(spatial_coords_path, row.names=1, check.names=FALSE)
cat(sprintf("  Spatial: %d genes x %d AOIs\n", nrow(spatial_counts), ncol(spatial_counts)))

coords <- as.data.frame(spatial_coords)
colnames(coords) <- c("x", "y")

cat("\nCreating CARD object...\n")
sc_meta <- data.frame(
    cellID = rownames(ref_meta),
    cellType = ref_meta$cell_type,
    sampleInfo = if ("sample" %in% colnames(ref_meta)) ref_meta$sample else "sample1",
    row.names = rownames(ref_meta)
)

CARD_obj <- createCARDObject(
    sc_count = ref_counts,
    sc_meta = sc_meta,
    spatial_count = spatial_counts,
    spatial_location = coords,
    ct.varname = "cellType",
    ct.select = unique(sc_meta$cellType),
    sample.varname = "sampleInfo",
    minCountGene = 100,
    minCountSpot = 5
)

cat("\nRunning CARD deconvolution...\n")
start_time <- Sys.time()
CARD_obj <- CARD_deconvolution(CARD_obj)
elapsed <- as.numeric(difftime(Sys.time(), start_time, units="mins"))
cat(sprintf("CARD completed in %.1f minutes\n", elapsed))

proportions <- as.data.frame(CARD_obj@Proportion_CARD)
prop_path <- file.path(output_dir, "predicted_proportions.csv")
write.csv(proportions, prop_path)
cat(sprintf("Wrote: %s (%d x %d)\n", prop_path, nrow(proportions), ncol(proportions)))

saveRDS(CARD_obj, file.path(output_dir, "card_object.rds"))
cat(sprintf("\nDone: %s\n", Sys.time()))
