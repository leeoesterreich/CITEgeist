#!/usr/bin/env Rscript
# Run Seurat label transfer on pseudo-GeoMx AOI data.

library(optparse)
library(Seurat)

option_list <- list(
  make_option(c("--ref-counts"), type="character",
              default="reference_data/processed/seurat/counts.csv"),
  make_option(c("--ref-cell-types"), type="character",
              default="reference_data/processed/seurat/cell_types.csv"),
  make_option(c("--spatial-counts"), type="character",
              default="data/aoi_counts.csv"),
  make_option(c("--output-dir"), type="character",
              default="Seurat/output")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

cat("============================================\n")
cat("Seurat Label Transfer - Pseudo-GeoMx\n")
cat("============================================\n")
cat(sprintf("Start time: %s\n", Sys.time()))

base_dir <- getwd()
resolve <- function(p) normalizePath(file.path(base_dir, p), mustWork=TRUE)

ref_counts_path <- resolve(opt$`ref-counts`)
ref_ct_path <- resolve(opt$`ref-cell-types`)
spatial_counts_path <- resolve(opt$`spatial-counts`)
output_dir <- file.path(base_dir, opt$`output-dir`)

dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)

cat("\nLoading reference...\n")
ref_counts <- as.matrix(read.csv(ref_counts_path, row.names=1, check.names=FALSE))
ref_meta <- read.csv(ref_ct_path, row.names=1, check.names=FALSE)
cat(sprintf("  Reference: %d genes x %d cells\n", nrow(ref_counts), ncol(ref_counts)))

# Seurat replaces underscores with dashes in feature names — do this pre-dedup
# to catch cases where GENE_X and GENE-X become duplicates after replacement
rownames(ref_counts) <- gsub("_", "-", rownames(ref_counts))
dup_genes <- duplicated(rownames(ref_counts))
if (any(dup_genes)) {
    cat(sprintf("  Removing %d duplicate gene names\n", sum(dup_genes)))
    ref_counts <- ref_counts[!dup_genes, ]
}

ref_obj <- CreateSeuratObject(counts=ref_counts, meta.data=ref_meta)
ref_obj <- NormalizeData(ref_obj)
ref_obj <- FindVariableFeatures(ref_obj, nfeatures=2000)
ref_obj <- ScaleData(ref_obj)
ref_obj <- RunPCA(ref_obj, npcs=30)

cat("\nLoading spatial data...\n")
spatial_counts <- as.matrix(read.csv(spatial_counts_path, row.names=1, check.names=FALSE))
cat(sprintf("  Spatial: %d genes x %d AOIs\n", nrow(spatial_counts), ncol(spatial_counts)))

spatial_obj <- CreateSeuratObject(counts=spatial_counts)
spatial_obj <- NormalizeData(spatial_obj)

cat("\nFinding transfer anchors...\n")
start_time <- Sys.time()

anchors <- FindTransferAnchors(
    reference=ref_obj,
    query=spatial_obj,
    dims=1:30,
    reference.reduction="pca"
)

predictions <- TransferData(
    anchorset=anchors,
    refdata=ref_obj$cell_type,
    dims=1:30,
    k.weight=20
)

elapsed <- as.numeric(difftime(Sys.time(), start_time, units="mins"))
cat(sprintf("Transfer completed in %.1f minutes\n", elapsed))

cat("\nExtracting proportions...\n")
score_cols <- grep("^prediction.score.", colnames(predictions), value=TRUE)
score_cols <- score_cols[score_cols != "prediction.score.max"]

proportions <- predictions[, score_cols]
colnames(proportions) <- gsub("prediction.score.", "", colnames(proportions))
colnames(proportions) <- gsub("\\.", " ", colnames(proportions))
colnames(proportions) <- gsub("CD4  T cells", "CD4+ T cells", colnames(proportions))
colnames(proportions) <- gsub("CD8  T cells", "CD8+ T cells", colnames(proportions))

prop_path <- file.path(output_dir, "predicted_proportions.csv")
write.csv(proportions, prop_path)
cat(sprintf("Wrote: %s (%d x %d)\n", prop_path, nrow(proportions), ncol(proportions)))

cat(sprintf("\nDone: %s\n", Sys.time()))
