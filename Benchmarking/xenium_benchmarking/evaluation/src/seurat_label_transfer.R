#!/usr/bin/env Rscript
# Seurat v5 label transfer for Xenium single-cell RNA data.
#
# Transfers cell type labels from GSE156632 scRNA-seq reference (RCC)
# to Xenium single-cell gene expression (405 genes) using
# FindTransferAnchors() + TransferData().
#
# Usage:
#   Rscript seurat_label_transfer.R \
#     --reference /path/to/GSE156632_combined.h5ad \
#     --query_dir /path/to/xenium_data \
#     --output_dir /path/to/output

suppressPackageStartupMessages({
    library(Seurat)
    library(Matrix)
    library(optparse)
    library(jsonlite)
})

# Parse command line arguments
option_list <- list(
    make_option(c("--reference"), type = "character",
        default = "/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/reference_data/GSE156632/processed/seurat/mtx",
        help = "Path to reference 10x MTX directory (with matrix.mtx.gz, barcodes.tsv.gz, features.tsv.gz)"),
    make_option(c("--cell_types"), type = "character",
        default = "/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/reference_data/GSE156632/processed/seurat/mtx/cell_types.csv",
        help = "Path to cell_types.csv with cell_id and cell_type columns"),
    make_option(c("--query_dir"), type = "character",
        default = "/ix1/alee/LO_LAB/General/Public_Data/10x_PublicData/Xenium_RNA_Proteomic_RenalCellCarcinoma",
        help = "Path to Xenium data directory"),
    make_option(c("--output_dir"), type = "character",
        default = "/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/evaluation/results/seurat_label_transfer",
        help = "Output directory"),
    make_option(c("--n_dims"), type = "integer", default = 30,
        help = "Number of PCA dimensions for transfer [default: 30]"),
    make_option(c("--n_features"), type = "integer", default = 2000,
        help = "Number of variable features [default: 2000]")
)

opt <- parse_args(OptionParser(option_list = option_list))

cat("=== Seurat Label Transfer for Xenium ===\n")
cat(sprintf("Reference: %s\n", opt$reference))
cat(sprintf("Query: %s\n", opt$query_dir))
cat(sprintf("Output: %s\n", opt$output_dir))

dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================================
# 1. Load reference scRNA-seq data from 10x MTX format
# ============================================================================
cat("\n--- Loading reference scRNA-seq data (10x MTX format) ---\n")
cat(sprintf("Loading from: %s\n", opt$reference))

ref_counts <- Read10X(data.dir = opt$reference)
ref <- CreateSeuratObject(counts = ref_counts, project = "GSE156632_ref")
cat(sprintf("Reference: %d cells x %d genes\n", ncol(ref), nrow(ref)))

# Load cell type annotations from CSV
cat(sprintf("Loading cell types from: %s\n", opt$cell_types))
ct_df <- read.csv(opt$cell_types, stringsAsFactors = FALSE)
rownames(ct_df) <- ct_df$cell_id

# Match to reference cells
common <- intersect(colnames(ref), ct_df$cell_id)
cat(sprintf("Cells with annotations: %d / %d\n", length(common), ncol(ref)))
ref <- ref[, common]
ref$cell_type <- ct_df[common, "cell_type"]
ref_labels <- "cell_type"

cat(sprintf("Using reference labels from: %s\n", ref_labels))
cat("Reference cell types:\n")
print(table(ref@meta.data[[ref_labels]]))

# Normalize reference
cat("\nNormalizing reference...\n")
ref <- NormalizeData(ref, verbose = FALSE)
ref <- FindVariableFeatures(ref, nfeatures = opt$n_features, verbose = FALSE)
ref <- ScaleData(ref, verbose = FALSE)
ref <- RunPCA(ref, npcs = opt$n_dims, verbose = FALSE)

# ============================================================================
# 2. Load Xenium query data
# ============================================================================
cat("\n--- Loading Xenium query data ---\n")

# Load from 10x h5 format
h5_path <- file.path(opt$query_dir, "cell_feature_matrix.h5")
if (!file.exists(h5_path)) {
    stop(sprintf("Cell feature matrix not found: %s", h5_path))
}

xenium <- Read10X_h5(h5_path)

# Xenium h5 may contain both Gene Expression and Protein Expression
# Extract just gene expression
if (is.list(xenium)) {
    if ("Gene Expression" %in% names(xenium)) {
        gex_matrix <- xenium[["Gene Expression"]]
    } else {
        gex_matrix <- xenium[[1]]
    }
} else {
    gex_matrix <- xenium
}

# Create Seurat object for query
query <- CreateSeuratObject(counts = gex_matrix, project = "Xenium_query")
cat(sprintf("Query: %d cells x %d genes\n", ncol(query), nrow(query)))

# Load spatial coordinates from cells.parquet
parquet_path <- file.path(opt$query_dir, "cells.parquet")
if (file.exists(parquet_path)) {
    cat("Loading spatial coordinates from cells.parquet...\n")
    cells_df <- arrow::read_parquet(parquet_path)
    cells_df$cell_id <- as.character(cells_df$cell_id)

    # Match cell order
    common_cells <- intersect(colnames(query), cells_df$cell_id)
    if (length(common_cells) > 0) {
        query <- query[, common_cells]
        cells_df <- cells_df[match(common_cells, cells_df$cell_id), ]
        query@meta.data$x_centroid <- cells_df$x_centroid
        query@meta.data$y_centroid <- cells_df$y_centroid
    }
}

# Normalize query
cat("Normalizing query...\n")
query <- NormalizeData(query, verbose = FALSE)
query <- FindVariableFeatures(query, nfeatures = opt$n_features, verbose = FALSE)

# ============================================================================
# 3. Find transfer anchors and transfer labels
# ============================================================================
cat("\n--- Finding transfer anchors ---\n")

# Find common genes
common_genes <- intersect(rownames(ref), rownames(query))
cat(sprintf("Common genes between reference and query: %d\n", length(common_genes)))

if (length(common_genes) < 50) {
    stop("Too few common genes for label transfer")
}

# Find anchors
anchors <- FindTransferAnchors(
    reference = ref,
    query = query,
    dims = 1:opt$n_dims,
    features = common_genes,
    reduction = "pcaproject",
    verbose = TRUE
)

cat(sprintf("Found %d transfer anchors\n", nrow(anchors@anchors)))

# Transfer labels
cat("\n--- Transferring labels ---\n")
predictions <- TransferData(
    anchorset = anchors,
    refdata = ref@meta.data[[ref_labels]],
    dims = 1:opt$n_dims,
    verbose = TRUE
)

# Add predictions to query metadata
query <- AddMetaData(query, metadata = predictions)

cat("\nPredicted cell types:\n")
print(table(query$predicted.id))

# ============================================================================
# 4. Compute prediction confidence metrics
# ============================================================================
cat("\n--- Computing confidence metrics ---\n")

# Prediction scores are in predictions columns
pred_scores <- predictions[, grep("^prediction.score", colnames(predictions))]
max_scores <- unname(apply(pred_scores, 1, max))

query$prediction_score <- max_scores
query$high_confidence <- max_scores > 0.5

cat(sprintf("High-confidence predictions (>0.5): %d/%d (%.1f%%)\n",
    sum(query$high_confidence),
    ncol(query),
    100 * mean(query$high_confidence)))

# ============================================================================
# 5. Map to broad lineages
# ============================================================================
lineage_map <- c(
    "T cells" = "Immune",
    "CD4+ T cells" = "Immune",
    "CD8+ T cells" = "Immune",
    "Treg" = "Immune",
    "NK cells" = "Immune",
    "B cells" = "Immune",
    "Plasma cells" = "Immune",
    "Macrophages" = "Immune",
    "M1 Macrophages" = "Immune",
    "M2 Macrophages" = "Immune",
    "Monocytes" = "Immune",
    "Dendritic cells" = "Immune",
    "DC" = "Immune",
    "pDC" = "Immune",
    "Mast cells" = "Immune",
    "Epithelial" = "Epithelial",
    "Tumor" = "Epithelial",
    "ccRCC" = "Epithelial",
    "pRCC" = "Epithelial",
    "Proximal tubule" = "Epithelial",
    "Distal tubule" = "Epithelial",
    "Collecting duct" = "Epithelial",
    "Endothelial" = "Stromal",
    "Fibroblasts" = "Stromal",
    "Fibroblast" = "Stromal",
    "Myofibroblast" = "Stromal",
    "Pericytes" = "Stromal",
    "Smooth muscle" = "Stromal"
)

# Map predicted types to lineages (unname to avoid Seurat metadata name conflicts)
predicted_lineage <- unname(sapply(as.character(query$predicted.id), function(x) {
    if (x %in% names(lineage_map)) {
        return(unname(lineage_map[x]))
    }
    return("Other")
}))
query$predicted_lineage <- predicted_lineage

cat("\nPredicted lineages:\n")
print(table(query$predicted_lineage))

# ============================================================================
# 6. Save results
# ============================================================================
cat("\n--- Saving results ---\n")

# Save per-cell predictions (build from metadata to avoid indexing issues)
results_df <- data.frame(
    cell_id = colnames(query),
    predicted_type = as.character(query$predicted.id),
    prediction_score = as.numeric(query$prediction_score),
    predicted_lineage = as.character(query$predicted_lineage),
    high_confidence = as.logical(query$high_confidence),
    stringsAsFactors = FALSE,
    row.names = NULL
)

# Add spatial coordinates if available
if ("x_centroid" %in% colnames(query@meta.data)) {
    results_df$x_centroid <- query$x_centroid
    results_df$y_centroid <- query$y_centroid
}

write.csv(results_df, file.path(opt$output_dir, "seurat_label_transfer.csv"),
          row.names = FALSE)

# Save summary
type_counts <- as.list(table(query$predicted.id))
lineage_counts <- as.list(table(query$predicted_lineage))

summary_data <- list(
    n_cells = ncol(query),
    n_common_genes = length(common_genes),
    n_anchors = nrow(anchors@anchors),
    reference_labels_column = ref_labels,
    n_reference_types = length(unique(ref@meta.data[[ref_labels]])),
    type_counts = type_counts,
    lineage_counts = lineage_counts,
    mean_prediction_score = mean(query$prediction_score),
    high_confidence_rate = mean(query$high_confidence)
)

write_json(summary_data, file.path(opt$output_dir, "seurat_summary.json"),
           pretty = TRUE, auto_unbox = TRUE)

cat("\n=== SEURAT LABEL TRANSFER COMPLETE ===\n")
cat(sprintf("Results saved to: %s\n", opt$output_dir))
cat(sprintf("Total cells: %d\n", ncol(query)))
cat(sprintf("Common genes: %d\n", length(common_genes)))
cat(sprintf("Mean prediction score: %.3f\n", mean(query$prediction_score)))
cat(sprintf("High-confidence rate: %.1f%%\n", 100 * mean(query$high_confidence)))
