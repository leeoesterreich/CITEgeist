#!/usr/bin/env Rscript
# SingleR label transfer: GSE156632 -> Xenium cells
#
# Usage: Rscript run_singler.R \
#   --xenium-counts <cells x genes CSV> \
#   --ref-mtx-dir <dir with matrix.mtx.gz, barcodes.tsv.gz, features.tsv.gz> \
#   --ref-labels <cell_types.csv> \
#   --output-dir <path> \
#   --ncpus 16

library(optparse)
library(SingleR)
library(BiocParallel)
library(Matrix)
library(Seurat)
library(data.table)

option_list <- list(
  make_option("--xenium-counts", type="character", help="Xenium GEX CSV (cells x genes, raw counts)"),
  make_option("--ref-mtx-dir", type="character", help="Dir with MTX reference (matrix.mtx.gz etc)"),
  make_option("--ref-labels", type="character", help="Reference cell_types.csv"),
  make_option("--output-dir", type="character", default="singler_labels"),
  make_option("--ncpus", type="integer", default=16)
)

opt <- parse_args(OptionParser(option_list=option_list))

# --- Load reference (MTX format) ---
cat("Loading reference from MTX...\n")
ref_mat <- Read10X(opt$`ref-mtx-dir`)  # Returns genes x cells sparse matrix
ref_labels_df <- read.csv(opt$`ref-labels`, row.names=1)
keep <- colnames(ref_mat) %in% rownames(ref_labels_df)
ref_mat <- ref_mat[, keep]
ref_labels <- ref_labels_df$cell_type[match(colnames(ref_mat), rownames(ref_labels_df))]
cat(sprintf("Reference: %d genes x %d cells, %d labels\n",
            nrow(ref_mat), ncol(ref_mat), length(ref_labels)))
cat("Reference label distribution:\n")
print(table(ref_labels))

# Log-normalize reference
cat("Log-normalizing reference...\n")
ref_norm <- LogNormalize(ref_mat)  # Seurat's log-normalize (log1p of CPM to 1e4)

# --- Load Xenium GEX ---
cat("Loading Xenium GEX counts...\n")
xenium_dt <- fread(opt$`xenium-counts`, header=TRUE)
cell_ids <- xenium_dt[[1]]  # First column is cell_id (row names)
xenium_dt[, 1 := NULL]      # Remove cell_id column
xenium_mat <- as.matrix(xenium_dt)
rownames(xenium_mat) <- cell_ids
cat(sprintf("Xenium: %d cells x %d genes\n", nrow(xenium_mat), ncol(xenium_mat)))

# Find common genes
common_genes <- intersect(colnames(xenium_mat), rownames(ref_norm))
cat(sprintf("Common genes: %d\n", length(common_genes)))

if (length(common_genes) < 50) {
  stop("Too few common genes for reliable annotation!")
}

# Log-normalize Xenium
cat("Log-normalizing Xenium...\n")
xenium_sparse <- Matrix(t(xenium_mat[, common_genes]), sparse=TRUE)  # genes x cells
xenium_norm <- LogNormalize(xenium_sparse)

# Subset reference to common genes
ref_norm_sub <- ref_norm[common_genes, ]

# --- Run SingleR ---
cat("Running SingleR...\n")
cat(sprintf("Using %d CPUs\n", opt$ncpus))
results <- SingleR(
  test = xenium_norm,       # genes x cells (log-normalized)
  ref = ref_norm_sub,       # genes x cells (log-normalized)
  labels = ref_labels,
  BPPARAM = MulticoreParam(opt$ncpus)
)

# --- Extract and save results ---
out_df <- data.frame(
  cell_id = cell_ids,
  singler_label = results$labels,
  singler_score = apply(results$scores, 1, max),
  delta_next = results$delta.next,
  pruned = is.na(results$pruned.labels),
  row.names = NULL
)

cat("Label distribution:\n")
print(table(out_df$singler_label))
cat(sprintf("Pruned cells: %d (%.1f%%)\n", sum(out_df$pruned), 100*mean(out_df$pruned)))

# Per-type median delta.next
cat("Per-type median delta.next:\n")
for (ct in unique(out_df$singler_label)) {
  med <- median(out_df$delta_next[out_df$singler_label == ct], na.rm=TRUE)
  cat(sprintf("  %s: %.4f\n", ct, med))
}

dir.create(opt$`output-dir`, recursive=TRUE, showWarnings=FALSE)
write.csv(out_df, file.path(opt$`output-dir`, "singler_all_cells.csv"), row.names=FALSE)
cat("Done.\n")
