#!/usr/bin/env Rscript
# Cluster-level SingleR naming: one label per level-1 Leiden cluster (SingleR(clusters=)).
# Loading, LogNormalize and gene intersection are copied from ground_truth_singler/src/run_singler.R.
#
# Usage: Rscript name_clusters.R --xenium-counts <cells x genes CSV> --clusters-csv <level1_clusters.csv> \
#   --ref-mtx-dir <dir> --ref-labels <cell_types.csv> --output-dir <path> --ncpus 16
# Writes cluster_names.csv (cluster,label,delta_next,n_cells,score_<label>...) and singler_versions.json.

library(optparse)
library(SingleR)
library(BiocParallel)
library(Matrix)
library(Seurat)
library(data.table)
library(jsonlite)

option_list <- list(
  make_option("--xenium-counts", type="character", help="Xenium GEX CSV (cells x genes, raw counts)"),
  make_option("--clusters-csv", type="character", help="level1_clusters.csv (cell_id, cluster)"),
  make_option("--ref-mtx-dir", type="character", help="Dir with MTX reference (matrix.mtx.gz etc)"),
  make_option("--ref-labels", type="character", help="Reference cell_types.csv (6 coarse labels)"),
  make_option("--output-dir", type="character", default="cluster_names"),
  make_option("--ncpus", type="integer", default=16)
)
opt <- parse_args(OptionParser(option_list=option_list))

# --- Reference ---
ref_mat <- Read10X(opt$`ref-mtx-dir`)
ref_labels_df <- read.csv(opt$`ref-labels`, row.names=1)
ref_labels <- ref_labels_df$cell_type[match(colnames(ref_mat), rownames(ref_labels_df))]
cat(sprintf("Reference: %d genes x %d cells\n", nrow(ref_mat), ncol(ref_mat)))
print(table(ref_labels))
ref_norm <- LogNormalize(ref_mat)

# --- Xenium ---
xenium_dt <- fread(opt$`xenium-counts`, header=TRUE)
cell_ids <- as.character(xenium_dt[[1]])
xenium_dt[, 1 := NULL]
xenium_mat <- as.matrix(xenium_dt)
rownames(xenium_mat) <- cell_ids

# --- Clusters: subset test cells to those in the clusters CSV ---
cl <- read.csv(opt$`clusters-csv`, colClasses="character")
cl_ids <- as.character(cl[[1]])
cl_vec <- cl$cluster
if (anyDuplicated(cl_ids) > 0) stop("Duplicated cell ids in clusters CSV")
if (any(is.na(cl_vec) | trimws(cl_vec) == "")) stop("NA or empty cluster value in clusters CSV")
keep <- intersect(cl_ids, cell_ids)
if (length(keep) < length(cl_ids)) {
  stop(sprintf("%d clusters-CSV cells are missing from the counts; every cell must inherit a label",
               length(cl_ids) - length(keep)))
}
cat(sprintf("Cells in counts: %d; in clusters CSV: %d; kept: %d\n", length(cell_ids), length(cl_ids), length(keep)))
clusters <- cl_vec[match(keep, cl_ids)]
xenium_mat <- xenium_mat[keep, , drop=FALSE]

common_genes <- intersect(colnames(xenium_mat), rownames(ref_norm))
cat(sprintf("Common genes: %d\n", length(common_genes)))
if (length(common_genes) < 50) stop("Too few common genes for reliable annotation!")

xenium_sparse <- Matrix(t(xenium_mat[, common_genes]), sparse=TRUE)
xenium_norm <- LogNormalize(xenium_sparse)
ref_norm_sub <- ref_norm[common_genes, ]

# --- SingleR at cluster level ---
res <- SingleR(
  test = xenium_norm,
  ref = ref_norm_sub,
  labels = ref_labels,
  clusters = clusters,
  BPPARAM = MulticoreParam(opt$ncpus)
)

# delta_next is SingleR's fine-tuned delta.next; label may differ from the argmax of the score_* columns.
scores <- as.data.frame(res$scores)
colnames(scores) <- paste0("score_", colnames(res$scores))
out <- data.frame(
  cluster = rownames(res),
  label = res$labels,
  delta_next = res$delta.next,
  n_cells = as.integer(table(clusters)[rownames(res)]),
  scores,
  check.names = FALSE,
  row.names = NULL
)
print(out[, c("cluster", "label", "delta_next", "n_cells")])

dir.create(opt$`output-dir`, recursive=TRUE, showWarnings=FALSE)
write.csv(out, file.path(opt$`output-dir`, "cluster_names.csv"), row.names=FALSE)
write_json(
  list(SingleR = as.character(packageVersion("SingleR")),
       R = paste(R.version$major, R.version$minor, sep="."),
       Seurat = as.character(packageVersion("Seurat"))),
  file.path(opt$`output-dir`, "singler_versions.json"), auto_unbox=TRUE, pretty=TRUE
)
cat("Done.\n")
