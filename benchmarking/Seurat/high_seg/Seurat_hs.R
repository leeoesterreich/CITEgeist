#!/usr/bin/env Rscript

# Capture command-line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop("Usage: Rscript seurat_deconvolution.R <file_path> <output_dir> <rep_name>")
}

file_path  <- args[1]
output_dir <- args[2]
rep_name   <- args[3]

# Suppress messages while loading libraries
suppressMessages(library(Seurat))
suppressMessages(library(SeuratDisk))
suppressMessages(library(ggplot2))

# Convert h5ad to h5Seurat
converted_file <- sub("\\.h5ad$", ".h5seurat", file_path)
cat(sprintf("[%s] Converting h5ad to h5Seurat: %s -> %s\n", Sys.time(), file_path, converted_file))
Convert(file_path, dest = "h5Seurat", overwrite = TRUE, verbose = TRUE)

# Load converted Seurat object
cat(sprintf("[%s] Loading converted Seurat object: %s\n", Sys.time(), converted_file))
query_obj <- LoadH5Seurat(converted_file)

# Load reference Seurat object
ref_file  <- "/bgfs/alee/LO_LAB/Personal/Brent_Schlegel/Projects/Wu_Visium/Simulations/models/Wu_scRNA_ref_ERpos.h5seurat"
cat(sprintf("[%s] Loading reference Seurat object: %s\n", Sys.time(), ref_file))
Wu_scRNA_ref <- LoadH5Seurat(ref_file)
Wu_scRNA_ref <- FindVariableFeatures(Wu_scRNA_ref)

# Function to perform Seurat-based deconvolution
run_spatial_deconvolution <- function(query_obj, ref_obj, rep_name, output_dir) {
  cat(sprintf("[%s] Running spatial deconvolution for replicate: %s\n", Sys.time(), rep_name))
  
  # Ensure output directory exists
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  # Find PCA-based transfer anchors
  anchors <- FindTransferAnchors(
    reference = ref_obj,
    query = query_obj,
    scale = TRUE,
    reduction = "pcaproject",
    normalization.method = "LogNormalize"
  )
  
  # Predict query cell-type labels from reference
  predictions.assay <- TransferData(
    anchorset = anchors,
    refdata = ref_obj$celltype,
    prediction.assay = TRUE
  )
  
  # Save predictions to CSV
  seurat_predictions <- as.data.frame(t(predictions.assay@data))
  seurat_predictions <- seurat_predictions[, c(
    "B-cells", "CAFs", "Cancer Epithelial", "Endothelial", "Myeloid",
    "Normal Epithelial", "PVL", "Plasmablasts", "T-cells"
  )]
  
  output_file <- file.path(output_dir, paste0(rep_name, "_Seurat_deconv_predictions.csv"))
  write.csv(seurat_predictions, output_file, row.names = TRUE)
  cat(sprintf("[%s] Saved deconvolution predictions: %s\n", Sys.time(), output_file))
  
  # Add predictions to the Seurat object
  query_obj[["predictions"]] <- predictions.assay
  query_obj@meta.data$prediction <- apply(seurat_predictions, 1, function(row) {
    colnames(seurat_predictions)[which.max(row)]
  })
  
  # Save updated Seurat object
  updated_file <- file.path(output_dir, paste0(rep_name, "_updated_query.h5Seurat"))
  SaveH5Seurat(query_obj, filename = updated_file, overwrite = TRUE)
  cat(sprintf("[%s] Saved updated Seurat object: %s\n", Sys.time(), updated_file))
  
  # Generate and save visualization
  plot_file <- file.path(output_dir, paste0(rep_name, "_Seurat_Prediction.pdf"))
  pdf(plot_file)
  print(DimPlot(query_obj, group.by = "prediction", reduction = "spatial") +
          ggtitle("Seurat Prediction", subtitle = rep_name))
  dev.off()
  cat(sprintf("[%s] Saved visualization: %s\n", Sys.time(), plot_file))
  
  # Display prediction summary
  cat(sprintf("[%s] Prediction summary for %s:\n", Sys.time(), rep_name))
  print(table(query_obj$prediction))
}

# Preprocess query object
query_obj@assays$RNA$data <- log1p(query_obj@assays$RNA$data)
query_obj <- ScaleData(query_obj)

# Run the deconvolution
run_spatial_deconvolution(query_obj, Wu_scRNA_ref, rep_name, output_dir)

cat(sprintf("[%s] Seurat deconvolution completed for %s\n", Sys.time(), rep_name))
