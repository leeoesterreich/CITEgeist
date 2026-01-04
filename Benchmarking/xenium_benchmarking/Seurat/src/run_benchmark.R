#!/usr/bin/env Rscript
#
# Run Seurat TransferData deconvolution on Xenium pseudo-Visium data.
#
# Seurat's TransferData uses label transfer to estimate cell type
# proportions in spatial transcriptomics data. It finds anchors between
# reference and query, then transfers labels with prediction scores.
#
# Note: Seurat only outputs cell type proportions (prediction scores),
# not per-cell-type gene expression layers.
#
# Usage:
#   Rscript run_benchmark.R --region-id 0 --ref-counts ref_counts.csv \
#     --ref-cell-types ref_cell_types.csv --spatial-counts spatial_counts.csv \
#     --output-dir output/

library(optparse)
library(Seurat)
library(dplyr)

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
  make_option(c("--output-dir"), type="character", default="output_rna_gt",
              help="Output directory"),
  make_option(c("--prefix"), type="character", default="Xenium",
              help="Filename prefix"),
  make_option(c("--n-features"), type="integer", default=2000,
              help="Number of variable features")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

cat("============================================\n")
cat(sprintf("Seurat TransferData - Region %d\n", opt$`region-id`))
cat("============================================\n")
cat(sprintf("Start time: %s\n", Sys.time()))
cat(sprintf("Reference counts: %s\n", opt$`ref-counts`))
cat(sprintf("Reference cell types: %s\n", opt$`ref-cell-types`))
cat(sprintf("Spatial counts: %s\n", opt$`spatial-counts`))
cat(sprintf("Output directory: %s\n", opt$`output-dir`))
cat("\n")

# Load reference data
cat("Loading reference data...\n")
ref_counts <- as.matrix(read.csv(opt$`ref-counts`, row.names=1, check.names=FALSE))
ref_cell_types <- read.csv(opt$`ref-cell-types`, row.names=1, check.names=FALSE)

cat(sprintf("  Reference counts: %d genes x %d cells\n", nrow(ref_counts), ncol(ref_counts)))
cat(sprintf("  Cell types: %s\n", paste(unique(ref_cell_types$cell_type), collapse=", ")))

# Create Reference Seurat object
cat("Creating Reference Seurat object...\n")
ref_seurat <- CreateSeuratObject(
  counts = ref_counts,
  project = "Reference"
)
ref_seurat$cell_type <- ref_cell_types[colnames(ref_seurat), "cell_type"]

# Normalize and find variable features
cat("  Normalizing reference...\n")
ref_seurat <- NormalizeData(ref_seurat, verbose = FALSE)

# Handle FindVariableFeatures with error handling for constant expression genes
cat("  Finding variable features...\n")
tryCatch({
  ref_seurat <- FindVariableFeatures(ref_seurat, selection.method = "vst",
                                      nfeatures = opt$`n-features`, verbose = FALSE)
}, error = function(e) {
  cat("    VST failed, trying mean.var.plot method...\n")
  ref_seurat <<- FindVariableFeatures(ref_seurat, selection.method = "mean.var.plot",
                                       nfeatures = opt$`n-features`, verbose = FALSE)
})

ref_seurat <- ScaleData(ref_seurat, verbose = FALSE)
ref_seurat <- RunPCA(ref_seurat, verbose = FALSE)

# Load spatial data
cat("\nLoading spatial data...\n")
spatial_counts <- as.matrix(read.csv(opt$`spatial-counts`, row.names=1, check.names=FALSE))
cat(sprintf("  Spatial counts: %d genes x %d spots\n", nrow(spatial_counts), ncol(spatial_counts)))

# Create Spatial Seurat object
cat("Creating Spatial Seurat object...\n")
spatial_seurat <- CreateSeuratObject(
  counts = spatial_counts,
  project = "Spatial"
)

# Normalize
cat("  Normalizing spatial data...\n")
spatial_seurat <- NormalizeData(spatial_seurat, verbose = FALSE)

# Run TransferData
cat("\nRunning Seurat TransferData...\n")
start_time <- Sys.time()

# Find transfer anchors
cat("  Finding transfer anchors...\n")
anchors <- FindTransferAnchors(
  reference = ref_seurat,
  query = spatial_seurat,
  normalization.method = "LogNormalize",
  reference.assay = "RNA",
  query.assay = "RNA",
  dims = 1:30,
  verbose = FALSE
)

# Transfer labels
cat("  Transferring labels...\n")
predictions <- TransferData(
  anchorset = anchors,
  refdata = ref_seurat$cell_type,
  dims = 1:30,
  weight.reduction = "pcaproject",
  verbose = FALSE
)

end_time <- Sys.time()
runtime <- as.numeric(difftime(end_time, start_time, units="secs"))
cat(sprintf("TransferData completed in %.1f seconds (%.1f min)\n", runtime, runtime/60))

# Add predictions to spatial object
spatial_seurat <- AddMetaData(spatial_seurat, predictions)

# Extract proportions (prediction scores)
cat("\nExtracting cell type proportions...\n")

# Get prediction score columns (prediction.score.*)
score_cols <- grep("^prediction.score\\.", colnames(predictions), value = TRUE)
score_cols <- score_cols[score_cols != "prediction.score.max"]

proportions <- predictions[, score_cols]
colnames(proportions) <- gsub("prediction.score.", "", colnames(proportions))

# Ensure proportions sum to 1
row_sums <- rowSums(proportions)
proportions <- proportions / row_sums

cat(sprintf("  Proportions shape: %d spots x %d cell types\n",
            nrow(proportions), ncol(proportions)))

# Create output directory
sample_name <- sprintf("%s_region_%d", opt$prefix, opt$`region-id`)
output_dir <- file.path(opt$`output-dir`, sample_name)
dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)

# Save proportions with the expected filename format
output_file <- file.path(output_dir, sprintf("%s_cell_prop_finetuned_results.csv", sample_name))
write.csv(proportions, output_file, row.names=TRUE)
cat(sprintf("Saved proportions to: %s\n", output_file))

# Also save a copy with the standard name for evaluation compatibility
output_file_alt <- file.path(output_dir, sprintf("%s_deconv_predictions.csv", sample_name))
write.csv(proportions, output_file_alt, row.names=TRUE)
cat(sprintf("Saved proportions (alt) to: %s\n", output_file_alt))

# Save predicted labels
labels_file <- file.path(output_dir, sprintf("%s_predicted_labels.csv", sample_name))
labels_df <- data.frame(
  spot = rownames(predictions),
  predicted_cell_type = predictions$predicted.id,
  prediction_score = predictions$prediction.score.max
)
write.csv(labels_df, labels_file, row.names=FALSE)
cat(sprintf("Saved predicted labels to: %s\n", labels_file))

# Save run statistics
stats <- list(
  region_id = opt$`region-id`,
  sample_name = sample_name,
  n_spots = nrow(proportions),
  n_cell_types = ncol(proportions),
  cell_types = colnames(proportions),
  runtime_sec = runtime,
  method = "Seurat_TransferData",
  n_features = opt$`n-features`,
  n_anchors = nrow(anchors@anchors)
)

stats_file <- file.path(output_dir, "run_stats.json")
cat(jsonlite::toJSON(stats, auto_unbox=TRUE, pretty=TRUE), file=stats_file)
cat(sprintf("Saved stats to: %s\n", stats_file))

# Save Seurat object for further analysis
seurat_file <- file.path(output_dir, sprintf("%s_Seurat.rds", sample_name))
saveRDS(spatial_seurat, seurat_file)
cat(sprintf("Saved Seurat object to: %s\n", seurat_file))

cat("\n============================================\n")
cat("Deconvolution Complete!\n")
cat("============================================\n")
cat(sprintf("  Spots: %d\n", nrow(proportions)))
cat(sprintf("  Cell types: %d\n", ncol(proportions)))
cat(sprintf("  Anchors: %d\n", nrow(anchors@anchors)))
cat(sprintf("  Runtime: %.1f seconds\n", runtime))
cat(sprintf("  Output: %s\n", output_dir))
cat(sprintf("End time: %s\n", Sys.time()))
