#!/usr/bin/env Rscript
#
# Seurat reference-sensitivity benchmark.
#
# Runs Seurat label transfer deconvolution with either full or downsampled
# reference on simulated Visium data. Uses the same hyperparameters as the
# main simulation benchmark (pcaproject reduction, LogNormalize).
#
# Unlike the main benchmark which loads h5seurat, this script loads from
# CSV (counts + cell_types) to avoid NFS h5seurat issues and to use the
# same reference format across all methods.
#
# Usage:
#   Rscript --vanilla run_seurat_refsens.R \
#     --replicate 0 \
#     --condition high_seg \
#     --ref-condition full \
#     --ref-dir ./reference \
#     --spatial-h5ad-dir /path/to/h5ad_objects \
#     --output-dir ./results/Seurat

.cg_paths_r <- local({
  cand <- Sys.getenv("CITEGEIST_PATHS_R", unset = "")
  if (!nzchar(cand)) {
    f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
    if (length(f)) cand <- file.path(dirname(normalizePath(f[1])), "..", "..", "_paths.R")
  }
  if (!nzchar(cand) || !file.exists(cand)) {
    d <- getwd()
    repeat {
      cand <- file.path(d, "Benchmarking", "_paths.R")
      if (file.exists(cand)) break
      if (dirname(d) == d) {
        stop("cannot find Benchmarking/_paths.R; set CITEGEIST_PATHS_R", call. = FALSE)
      }
      d <- dirname(d)
    }
  }
  cand
})
source(.cg_paths_r)
rm(.cg_paths_r)
library(optparse)
suppressMessages(library(Seurat))

option_list <- list(
  make_option(c("-r", "--replicate"), type = "integer", default = 0,
              help = "Replicate index (0-4)", metavar = "INT"),
  make_option(c("-c", "--condition"), type = "character", default = "high_seg",
              help = "Simulation condition: high_seg or mixed", metavar = "STR"),
  make_option(c("--ref-condition"), type = "character", default = "full",
              help = "Reference condition: full or downsampled", metavar = "STR"),
  make_option(c("--ref-dir"), type = "character",
              default = file.path(cg_root("REPO"), "Benchmarking", "simulation_benchmarking", "refsens", "reference"),
              help = "Directory containing full/ and downsampled/ subdirs", metavar = "PATH"),
  make_option(c("--spatial-h5ad-dir"), type = "character",
              default = file.path(cg_root("REPO"), "Benchmarking", "simulation_benchmarking", "replicates"),
              help = "Base directory containing condition/h5ad_objects/", metavar = "PATH"),
  make_option(c("-o", "--output-dir"), type = "character", default = ".",
              help = "Output directory for predictions", metavar = "PATH")
)

opt <- parse_args(OptionParser(option_list = option_list))

rep_index   <- opt$replicate
condition   <- opt$condition
ref_cond    <- opt[["ref-condition"]]
ref_dir     <- opt[["ref-dir"]]
spatial_dir <- opt[["spatial-h5ad-dir"]]
output_dir  <- opt[["output-dir"]]

cat(sprintf("=== Seurat Refsens: rep=%d, condition=%s, ref=%s ===\n",
            rep_index, condition, ref_cond))

# ── Load reference from CSV ─────────────────────────────────────────────────
ref_path <- file.path(ref_dir, ref_cond)
cat(sprintf("Loading reference from %s ...\n", ref_path))

counts_path <- file.path(ref_path, "counts.csv")
ct_path     <- file.path(ref_path, "cell_types.csv")

# counts.csv is genes x cells
sc_counts <- read.csv(counts_path, row.names = 1, check.names = FALSE)
sc_counts <- as.matrix(sc_counts)

ct_df <- read.csv(ct_path, row.names = 1)

cat(sprintf("  Reference: %d genes x %d cells\n", nrow(sc_counts), ncol(sc_counts)))

# Create Seurat reference object
ref_obj <- CreateSeuratObject(counts = sc_counts, project = "Wu_ref")
ref_obj$celltype <- ct_df[colnames(ref_obj), "cell_type"]
ref_obj <- NormalizeData(ref_obj, normalization.method = "LogNormalize")
ref_obj <- FindVariableFeatures(ref_obj)

cat(sprintf("  Cell types: %s\n", paste(unique(ref_obj$celltype), collapse = ", ")))

# ── Load spatial query data ─────────────────────────────────────────────────
# Load from the simulation CSV files (same as RCTD approach for consistency)
replicate_name <- paste0("Wu_ST_", rep_index)
st_sim_dir <- file.path(spatial_dir, condition, "ST_sim")
st_data_path <- file.path(st_sim_dir, paste0(replicate_name, "_data.csv"))

cat(sprintf("Loading spatial data: %s\n", st_data_path))

st_data <- read.csv(st_data_path, row.names = 1, check.names = FALSE)
st_data <- as.matrix(st_data)

# Create Seurat query object (spots as columns already in st_data)
query_obj <- CreateSeuratObject(counts = st_data, project = "spatial_query")
query_obj <- NormalizeData(query_obj, normalization.method = "LogNormalize")
query_obj <- ScaleData(query_obj)

cat(sprintf("  Query: %d genes x %d spots\n", nrow(st_data), ncol(st_data)))

# ── Run Seurat label transfer ───────────────────────────────────────────────
cat("Finding transfer anchors (pcaproject, LogNormalize) ...\n")

anchors <- FindTransferAnchors(
  reference = ref_obj,
  query = query_obj,
  scale = TRUE,
  reduction = "pcaproject",
  normalization.method = "LogNormalize"
)

predictions.assay <- TransferData(
  anchorset = anchors,
  refdata = ref_obj$celltype,
  prediction.assay = TRUE
)

# Extract and normalize predictions
seurat_predictions <- as.data.frame(t(predictions.assay@data))

# Get the cell type columns (exclude "max" prediction row)
ct_cols <- setdiff(colnames(seurat_predictions), "max")
seurat_predictions <- seurat_predictions[, ct_cols]

# Normalize to sum to 1 per spot
row_sums <- rowSums(seurat_predictions)
seurat_predictions <- seurat_predictions / row_sums

# ── Save ────────────────────────────────────────────────────────────────────
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_file <- file.path(
  output_dir,
  paste0(replicate_name, "_Seurat_deconv_predictions.csv")
)
write.csv(seurat_predictions, output_file, row.names = TRUE)
cat(sprintf("Saved predictions: %s\n", output_file))

# Write .done marker
done_file <- file.path(
  output_dir,
  paste0(replicate_name, "_Seurat.done")
)
writeLines(as.character(Sys.time()), done_file)

cat(sprintf("Completed Seurat refsens for %s (ref=%s)\n", replicate_name, ref_cond))
