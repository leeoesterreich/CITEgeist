#!/usr/bin/env Rscript
#
# Test CARD with VIM excluded to see if performance improves
#

library(CARD)
library(Matrix)

cat("============================================\n")
cat("CARD Test - Excluding VIM\n")
cat("============================================\n")

# Paths
REF_DIR <- "/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/reference_data/GSE156632/processed_protein7/rctd"
SPATIAL_DIR <- "/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/CARD/output_protein_gt/temp_csvs"

# Load data
cat("\nLoading data...\n")
spatial_counts <- as.matrix(read.csv(file.path(SPATIAL_DIR, "Xenium_region_0_counts.csv"), row.names=1, check.names=FALSE))
spatial_coords <- read.csv(file.path(SPATIAL_DIR, "Xenium_region_0_coords.csv"), row.names=1, check.names=FALSE)
ref_counts <- as.matrix(read.csv(file.path(REF_DIR, "counts.csv"), row.names=1, check.names=FALSE))
ref_cell_types <- read.csv(file.path(REF_DIR, "cell_types.csv"), row.names=1, check.names=FALSE)

cat(sprintf("  Spatial: %d genes x %d spots\n", nrow(spatial_counts), ncol(spatial_counts)))
cat(sprintf("  Reference: %d genes x %d cells\n", nrow(ref_counts), ncol(ref_counts)))

# EXCLUDE VIM from both spatial and reference
cat("\nExcluding VIM...\n")
vim_in_spatial <- "VIM" %in% rownames(spatial_counts)
vim_in_ref <- "VIM" %in% rownames(ref_counts)
cat(sprintf("  VIM in spatial: %s\n", vim_in_spatial))
cat(sprintf("  VIM in reference: %s\n", vim_in_ref))

if (vim_in_spatial) {
  spatial_counts <- spatial_counts[rownames(spatial_counts) != "VIM", ]
}
if (vim_in_ref) {
  ref_counts <- ref_counts[rownames(ref_counts) != "VIM", ]
}

cat(sprintf("  After exclusion - Spatial: %d genes\n", nrow(spatial_counts)))
cat(sprintf("  After exclusion - Reference: %d genes\n", nrow(ref_counts)))

# Prepare data
spatial_location <- as.data.frame(spatial_coords)
colnames(spatial_location) <- c("x", "y")

sc_meta <- data.frame(
  cellID = rownames(ref_cell_types),
  cellType = ref_cell_types$cell_type,
  sampleInfo = "sample1",
  row.names = rownames(ref_cell_types)
)

# Run CARD
cat("\nCreating CARD object (without VIM)...\n")
CARD_obj <- createCARDObject(
  sc_count = ref_counts,
  sc_meta = sc_meta,
  spatial_count = spatial_counts,
  spatial_location = spatial_location,
  ct.varname = "cellType",
  ct.select = unique(sc_meta$cellType),
  sample.varname = "sampleInfo",
  minCountGene = 100,
  minCountSpot = 5
)

cat("\nRunning CARD deconvolution...\n")
CARD_obj <- CARD_deconvolution(CARD_object = CARD_obj)

# Extract and summarize results
proportions <- as.data.frame(CARD_obj@Proportion_CARD)

cat("\n============================================\n")
cat("Results (VIM excluded)\n")
cat("============================================\n")
cat(sprintf("Proportions: %d spots x %d cell types\n", nrow(proportions), ncol(proportions)))

cat("\nMean proportions per cell type:\n")
means <- colMeans(proportions)
for (ct in names(sort(means, decreasing=TRUE))) {
  cat(sprintf("  %s: %.1f%%\n", ct, means[ct] * 100))
}

# Save results
output_dir <- "/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/CARD/output_protein_gt_novim"
dir.create(file.path(output_dir, "Xenium_region_0"), recursive=TRUE, showWarnings=FALSE)
output_file <- file.path(output_dir, "Xenium_region_0", "Xenium_region_0_cell_prop_finetuned_results.csv")
write.csv(proportions, output_file, row.names=TRUE)
cat(sprintf("\nSaved to: %s\n", output_file))
