#!/usr/bin/env Rscript
#
# RCTD reference-sensitivity benchmark.
#
# Runs RCTD deconvolution with either full or downsampled reference on
# simulated Visium data. Uses the same hyperparameters as the main
# simulation benchmark (doublet_mode='full', require_int=FALSE).
#
# Usage:
#   Rscript --vanilla run_rctd_refsens.R \
#     --replicate 0 \
#     --condition high_seg \
#     --ref-condition full \
#     --ref-dir ./reference \
#     --output-dir ./results/RCTD
#
# IMPORTANT: RCTD@config$max_cores is set to 1 to avoid PSOCK port 11235
# collision when running as array jobs.

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
library(spacexr)

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
  make_option(c("-o", "--output-dir"), type = "character", default = ".",
              help = "Output directory for predictions", metavar = "PATH")
)

opt <- parse_args(OptionParser(option_list = option_list))

rep_index   <- opt$replicate
condition   <- opt$condition
ref_cond    <- opt[["ref-condition"]]
ref_dir     <- opt[["ref-dir"]]
output_dir  <- opt[["output-dir"]]

cat(sprintf("=== RCTD Refsens: rep=%d, condition=%s, ref=%s ===\n",
            rep_index, condition, ref_cond))

# ── Load reference ──────────────────────────────────────────────────────────
ref_path <- file.path(ref_dir, ref_cond)
cat(sprintf("Loading reference from %s ...\n", ref_path))

counts_path <- file.path(ref_path, "counts.csv")
ct_path     <- file.path(ref_path, "cell_types.csv")

# counts.csv is genes x cells (rows=genes, cols=cells)
sc_data <- read.csv(counts_path, row.names = 1, check.names = FALSE)
sc_data <- as.matrix(sc_data)
sc_data <- round(sc_data)  # Ensure integers

ct_df <- read.csv(ct_path, row.names = 1)
cell_types <- as.factor(ct_df$cell_type)
names(cell_types) <- rownames(ct_df)

cat(sprintf("  Reference: %d genes x %d cells\n", nrow(sc_data), ncol(sc_data)))
cat(sprintf("  Cell types: %s\n", paste(levels(cell_types), collapse = ", ")))

reference <- Reference(sc_data, cell_types, require_int = FALSE)

# ── Load spatial data ───────────────────────────────────────────────────────
replicate_name <- paste0("Wu_ST_", rep_index)
st_sim_dir <- file.path(
  file.path(cg_root("REPO"), "Benchmarking", "simulation_benchmarking", "replicates"),
  condition, "ST_sim"
)
st_data_path <- file.path(st_sim_dir, paste0(replicate_name, "_data.csv"))
st_meta_path <- file.path(st_sim_dir, paste0(replicate_name, "_meta.csv"))

cat(sprintf("Loading spatial data: %s\n", st_data_path))

st_data <- read.csv(st_data_path, row.names = 1)
st_meta <- read.csv(st_meta_path)
st_data <- round(st_data)  # Ensure integers

coords <- data.frame(x = st_meta$spot_x, y = st_meta$spot_y)
rownames(coords) <- colnames(st_data)

puck <- SpatialRNA(coords, st_data, require_int = FALSE)

# ── Run RCTD ────────────────────────────────────────────────────────────────
cat("Creating RCTD object ...\n")
myRCTD <- create.RCTD(puck, reference)

# CRITICAL: set max_cores to 1 to avoid PSOCK port collision in array jobs
myRCTD@config$max_cores <- 1

cat("Running RCTD (doublet_mode='full', max_cores=1) ...\n")
myRCTD <- run.RCTD(myRCTD, doublet_mode = "full")
results <- myRCTD@results

# Normalize weights
norm_weights <- normalize_weights(results$weights)
norm_weights <- as.matrix(norm_weights)

# ── Save ────────────────────────────────────────────────────────────────────
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_file <- file.path(
  output_dir,
  paste0(replicate_name, "_RCTD_deconv_predictions.csv")
)
write.csv(norm_weights, output_file)
cat(sprintf("Saved predictions: %s\n", output_file))

# Write .done marker
done_file <- file.path(
  output_dir,
  paste0(replicate_name, "_RCTD.done")
)
writeLines(as.character(Sys.time()), done_file)

cat(sprintf("Completed RCTD refsens for %s (ref=%s)\n", replicate_name, ref_cond))
