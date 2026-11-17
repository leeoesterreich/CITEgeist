#!/usr/bin/env Rscript
# Convert Wu_scRNA_ref_ERpos.h5ad to h5seurat and Rds formats.
#
# SeuratDisk::Convert writes a valid h5seurat, but LoadH5Seurat fails on
# scanpy categorical columns ("Missing required datasets 'levels' and 'values'").
# Workaround: read sparse matrix + metadata from h5seurat via hdf5r directly,
# decoding categorical groups (categories + codes) manually.
#
# Usage:
#   Rscript convert_reference_formats.R
#   Rscript convert_reference_formats.R /path/to/Wu_scRNA_ref_ERpos.h5ad

suppressMessages(library(Seurat))
suppressMessages(library(SeuratDisk))
suppressMessages(library(Matrix))
suppressMessages(library(hdf5r))

args <- commandArgs(trailingOnly = TRUE)
h5ad_path <- if (length(args) >= 1) args[1] else file.path(getwd(), "Benchmarking/simulation_benchmarking/models/Wu_scRNA_ref_ERpos.h5ad")

cat(sprintf("[%s] Converting: %s\n", Sys.time(), h5ad_path))

if (!file.exists(h5ad_path)) {
  stop(sprintf("Input h5ad not found: %s", h5ad_path))
}

# Step 1: h5ad -> h5seurat
h5seurat_path <- sub("\\.h5ad$", ".h5seurat", h5ad_path)
cat(sprintf("[%s] Converting h5ad -> h5seurat: %s\n", Sys.time(), h5seurat_path))
Convert(h5ad_path, dest = "h5seurat", overwrite = TRUE)
cat(sprintf("[%s] h5seurat written: %.1f MB\n", Sys.time(), file.size(h5seurat_path) / 1e6))

# Step 2: Build Seurat object from h5seurat via hdf5r
cat(sprintf("[%s] Building Seurat object from h5seurat via hdf5r...\n", Sys.time()))
h5 <- H5File$new(h5seurat_path, mode = "r")

# Read sparse counts matrix
assay_group <- h5[["assays/RNA"]]
if ("counts" %in% names(assay_group)) {
  data_group <- assay_group[["counts"]]
} else {
  data_group <- assay_group[["data"]]
}
i_vals <- data_group[["indices"]]$read() + 1L
p_vals <- data_group[["indptr"]]$read()
x_vals <- data_group[["data"]]$read()
features <- h5[["assays/RNA/features"]]$read()
barcodes <- h5[["cell.names"]]$read()

counts <- sparseMatrix(i = i_vals, p = p_vals, x = x_vals,
                       dims = c(length(features), length(barcodes)),
                       dimnames = list(features, barcodes))

# Read metadata — handle both plain datasets and categorical groups
meta_group <- h5[["meta.data"]]
meta_df <- data.frame(row.names = barcodes)
for (col in names(meta_group)) {
  item <- meta_group[[col]]
  if (inherits(item, "H5Group")) {
    # Categorical: has "categories" (levels) and "codes" (integer indices)
    if ("categories" %in% names(item) && "codes" %in% names(item)) {
      cats <- item[["categories"]]$read()
      codes <- item[["codes"]]$read()
      meta_df[[col]] <- cats[codes + 1L]
      cat(sprintf("  Read categorical %s: %d levels\n", col, length(cats)))
    } else if ("levels" %in% names(item) && "values" %in% names(item)) {
      lvls <- item[["levels"]]$read()
      vals <- item[["values"]]$read()
      meta_df[[col]] <- lvls[vals + 1L]
      cat(sprintf("  Read factor %s: %d levels\n", col, length(lvls)))
    } else {
      cat(sprintf("  Skipping group %s (unknown structure: %s)\n", col, paste(names(item), collapse=", ")))
    }
  } else {
    meta_df[[col]] <- item$read()
    cat(sprintf("  Read dataset %s\n", col))
  }
}

h5$close_all()

# Create Seurat object
obj <- CreateSeuratObject(counts = counts, meta.data = meta_df)
cat(sprintf("[%s] Seurat object: %d cells x %d features\n", Sys.time(), ncol(obj), nrow(obj)))

# Save as Rds
rds_path <- sub("\\.h5ad$", ".Rds", h5ad_path)
cat(sprintf("[%s] Saving Rds: %s\n", Sys.time(), rds_path))
saveRDS(obj, rds_path)
cat(sprintf("[%s] Rds written: %.1f MB\n", Sys.time(), file.size(rds_path) / 1e6))

# Validate
cat(sprintf("[%s] Cell types: %s\n", Sys.time(),
    paste(sort(unique(obj@meta.data$celltype_major)), collapse = ", ")))
cat(sprintf("[%s] Done.\n", Sys.time()))
