#!/usr/bin/env Rscript
#
# Run CARD deconvolution on Xenium pseudo-Visium data.
#
# CARD (Conditional Autoregressive-based Deconvolution) uses a probabilistic
# model that incorporates spatial correlation to estimate cell type proportions.
#
# Supports two modes:
#   - reference: Uses scRNA-seq reference (standard CARD)
#   - reference_free: Uses marker gene lists only (CARD-Free)
#
# Usage:
#   Rscript run_benchmark.R --mode reference --region-id 0 \
#     --ref-counts ref_counts.csv --ref-cell-types ref_cell_types.csv \
#     --spatial-counts spatial_counts.csv --spatial-coords spatial_coords.csv \
#     --output-dir output/
#
#   Rscript run_benchmark.R --mode reference_free --region-id 0 \
#     --marker-genes marker_genes.csv \
#     --spatial-counts spatial_counts.csv --spatial-coords spatial_coords.csv \
#     --output-dir output/

library(optparse)
library(CARD)
library(Matrix)

# Parse command line arguments
option_list <- list(
  make_option(c("--mode"), type="character", default="reference",
              help="Deconvolution mode: 'reference' or 'reference_free'"),
  make_option(c("--region-id"), type="integer", default=0,
              help="Region ID (0-4)"),
  make_option(c("--ref-counts"), type="character", default=NULL,
              help="Reference counts CSV (genes x cells) - for reference mode"),
  make_option(c("--ref-cell-types"), type="character", default=NULL,
              help="Reference cell types CSV - for reference mode"),
  make_option(c("--marker-genes"), type="character", default=NULL,
              help="Marker genes CSV (cell_type, gene columns) - for reference_free mode"),
  make_option(c("--spatial-counts"), type="character", default=NULL,
              help="Spatial counts CSV (genes x spots)"),
  make_option(c("--spatial-coords"), type="character", default=NULL,
              help="Spatial coordinates CSV"),
  make_option(c("--output-dir"), type="character", default="output",
              help="Output directory"),
  make_option(c("--prefix"), type="character", default="Xenium",
              help="Filename prefix"),
  make_option(c("--min-count-gene"), type="integer", default=100,
              help="Minimum counts per gene across spots"),
  make_option(c("--min-count-spot"), type="integer", default=5,
              help="Minimum counts per spot")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

cat("============================================\n")
cat(sprintf("CARD Deconvolution - Region %d (%s mode)\n", opt$`region-id`, opt$mode))
cat("============================================\n")
cat(sprintf("Start time: %s\n", Sys.time()))
cat(sprintf("Mode: %s\n", opt$mode))
cat(sprintf("Spatial counts: %s\n", opt$`spatial-counts`))
cat(sprintf("Spatial coords: %s\n", opt$`spatial-coords`))
cat(sprintf("Output directory: %s\n", opt$`output-dir`))
cat("\n")

# Load spatial data
cat("Loading spatial data...\n")
spatial_counts <- as.matrix(read.csv(opt$`spatial-counts`, row.names=1, check.names=FALSE))
spatial_coords <- read.csv(opt$`spatial-coords`, row.names=1, check.names=FALSE)

cat(sprintf("  Spatial counts: %d genes x %d spots\n", nrow(spatial_counts), ncol(spatial_counts)))
cat(sprintf("  Coordinates: %d spots\n", nrow(spatial_coords)))

# Prepare spatial location data frame
spatial_location <- as.data.frame(spatial_coords)
colnames(spatial_location) <- c("x", "y")

# Run deconvolution based on mode
start_time <- Sys.time()

if (opt$mode == "reference") {
  # =========================================
  # Reference-based CARD
  # =========================================
  cat("\nRunning CARD with scRNA-seq reference...\n")

  if (is.null(opt$`ref-counts`) || is.null(opt$`ref-cell-types`)) {
    stop("Reference mode requires --ref-counts and --ref-cell-types")
  }

  cat(sprintf("Reference counts: %s\n", opt$`ref-counts`))
  cat(sprintf("Reference cell types: %s\n", opt$`ref-cell-types`))

  # Load reference data
  cat("\nLoading reference data...\n")
  ref_counts <- as.matrix(read.csv(opt$`ref-counts`, row.names=1, check.names=FALSE))
  ref_cell_types <- read.csv(opt$`ref-cell-types`, row.names=1, check.names=FALSE)

  cat(sprintf("  Reference counts: %d genes x %d cells\n", nrow(ref_counts), ncol(ref_counts)))
  cat(sprintf("  Cell types: %s\n", paste(unique(ref_cell_types$cell_type), collapse=", ")))

  # Create metadata for CARD
  sc_meta <- data.frame(
    cellID = rownames(ref_cell_types),
    cellType = ref_cell_types$cell_type,
    sampleInfo = "sample1",  # Single sample for simplicity
    row.names = rownames(ref_cell_types)
  )

  # Create CARD object
  cat("\nCreating CARD object...\n")
  CARD_obj <- createCARDObject(
    sc_count = ref_counts,
    sc_meta = sc_meta,
    spatial_count = spatial_counts,
    spatial_location = spatial_location,
    ct.varname = "cellType",
    ct.select = unique(sc_meta$cellType),
    sample.varname = "sampleInfo",
    minCountGene = opt$`min-count-gene`,
    minCountSpot = opt$`min-count-spot`
  )

  # Run deconvolution
  cat("\nRunning CARD deconvolution...\n")
  CARD_obj <- CARD_deconvolution(CARD_object = CARD_obj)

  # Extract proportions
  proportions <- as.data.frame(CARD_obj@Proportion_CARD)

} else if (opt$mode == "reference_free") {
  # =========================================
  # Reference-free CARD (marker-based)
  # =========================================
  cat("\nRunning CARD in reference-free mode...\n")

  if (is.null(opt$`marker-genes`)) {
    stop("Reference-free mode requires --marker-genes")
  }

  cat(sprintf("Marker genes: %s\n", opt$`marker-genes`))

  # Load marker genes
  marker_df <- read.csv(opt$`marker-genes`, check.names=FALSE)
  cat(sprintf("  Loaded markers for %d cell types\n", length(unique(marker_df$cell_type))))

  # Convert to list format expected by CARD
  # Format: list(CellType1 = c("gene1", "gene2"), CellType2 = c("gene3", "gene4"), ...)
  markerList <- split(marker_df$gene, marker_df$cell_type)

  for (ct in names(markerList)) {
    cat(sprintf("    %s: %d markers\n", ct, length(markerList[[ct]])))
  }

  # Create CARD-free object
  cat("\nCreating CARDfree object...\n")
  CARDfree_obj <- createCARDfreeObject(
    markerList = markerList,
    spatial_count = spatial_counts,
    spatial_location = spatial_location,
    minCountGene = opt$`min-count-gene`,
    minCountSpot = opt$`min-count-spot`
  )

  # Run reference-free deconvolution
  cat("\nRunning CARD reference-free deconvolution...\n")
  CARDfree_obj <- CARD_refFree(CARDfree_obj)

  # Extract proportions
  proportions <- as.data.frame(CARDfree_obj@Proportion_CARD)

} else {
  stop(sprintf("Unknown mode: %s. Use 'reference' or 'reference_free'", opt$mode))
}

end_time <- Sys.time()
runtime <- as.numeric(difftime(end_time, start_time, units="secs"))
cat(sprintf("\nCARD completed in %.1f seconds (%.1f min)\n", runtime, runtime/60))

cat(sprintf("  Proportions shape: %d spots x %d cell types\n",
            nrow(proportions), ncol(proportions)))

# Create output directory
sample_name <- sprintf("%s_region_%d", opt$prefix, opt$`region-id`)
output_subdir <- file.path(opt$`output-dir`, sample_name)
dir.create(output_subdir, recursive=TRUE, showWarnings=FALSE)

# Save proportions (standardized output filename matching other methods)
output_file <- file.path(output_subdir, sprintf("%s_cell_prop_finetuned_results.csv", sample_name))
write.csv(proportions, output_file, row.names=TRUE)
cat(sprintf("Saved proportions to: %s\n", output_file))

# Save run statistics
stats <- list(
  region_id = opt$`region-id`,
  sample_name = sample_name,
  mode = opt$mode,
  n_spots = nrow(proportions),
  n_cell_types = ncol(proportions),
  cell_types = colnames(proportions),
  runtime_sec = runtime,
  method = ifelse(opt$mode == "reference", "CARD", "CARD (ref-free)")
)

stats_file <- file.path(output_subdir, "run_stats.json")
cat(jsonlite::toJSON(stats, auto_unbox=TRUE, pretty=TRUE), file=stats_file)
cat(sprintf("Saved stats to: %s\n", stats_file))

# Save CARD object for further analysis
if (opt$mode == "reference") {
  card_file <- file.path(output_subdir, sprintf("%s_CARD.rds", sample_name))
  saveRDS(CARD_obj, card_file)
} else {
  card_file <- file.path(output_subdir, sprintf("%s_CARDfree.rds", sample_name))
  saveRDS(CARDfree_obj, card_file)
}
cat(sprintf("Saved CARD object to: %s\n", card_file))

cat("\n============================================\n")
cat("Deconvolution Complete!\n")
cat("============================================\n")
cat(sprintf("  Mode: %s\n", opt$mode))
cat(sprintf("  Spots: %d\n", nrow(proportions)))
cat(sprintf("  Cell types: %d\n", ncol(proportions)))
cat(sprintf("  Runtime: %.1f seconds\n", runtime))
cat(sprintf("  Output: %s\n", output_subdir))
cat(sprintf("End time: %s\n", Sys.time()))
