# Load required libraries
library(ggplot2)
library(dplyr)
library(stringr)
library(ggsignif)
library(multcompView)

# Define function to parse log files
parse_runtime <- function(file, method) {
  data <- readLines(file)
  df <- data.frame(
    Method = method,
    Replicate = as.numeric(str_extract(data, "(?<=Replicate )\\d+")),
    Runtime_Seconds = as.numeric(str_extract(data, "(?<=took )\\d+\\.?\\d*"))
  )
  return(df)
}

# Set file paths
files <- list(
  "CITEgeist" = "CITEgeist_high_seg_runtime.log",
  "RCTD" = "RCTD_high_seg_runtime.log",
  "Seurat" = "Seurat_high_seg_runtime.log",
  "Cell2Location" = "cell2location_high_seg_runtime.log",
  "Tangram" = "tangram_high_seg_runtime.log"
)

# Parse all files
runtime_data <- bind_rows(lapply(names(files), function(method) {
  parse_runtime(files[[method]], method)
}))

# Ensure 'Method' is a factor
runtime_data$Method <- factor(runtime_data$Method, levels = names(files))

# Run ANOVA
anova_res <- aov(Runtime_Seconds ~ Method, data = runtime_data)
anova_summary <- summary(anova_res)[[1]]
anova_Fstat <- anova_summary$`F value`[1]
anova_pval <- anova_summary$`Pr(>F)`[1]

# Perform Tukey's HSD post-hoc test
tukey_res <- TukeyHSD(anova_res)
tukey_df <- as.data.frame(tukey_res$Method)

# Extract only comparisons involving CITEgeist
tukey_df$Comparison <- rownames(tukey_df)
tukey_df <- tukey_df[grepl("CITEgeist", tukey_df$Comparison), ]

# Map p-values to significance levels
tukey_df$Significance <- cut(
  tukey_df$`p adj`,
  breaks = c(0, 0.001, 0.01, 0.05, 1),
  labels = c("***", "**", "*", "NS"),  # NS = Not Significant
  right = FALSE
)

# Remove non-significant comparisons
tukey_df <- tukey_df[tukey_df$Significance != "NS", ]

# Ensure proper formatting of comparison labels for ggsignif
tukey_df$Comparison <- gsub("-", " vs. ", tukey_df$Comparison)
tukey_df$groups <- strsplit(tukey_df$Comparison, " vs. ")

# Determine max y value for scaling
y_max <- max(runtime_data$Runtime_Seconds, na.rm = TRUE) * 1.2  # Scale up for annotations

# Plot runtimes
p <- ggplot(runtime_data, aes(x = Method, y = Runtime_Seconds, fill = Method)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 3, alpha = 0.6) +
  labs(
    title = "Comparison of Deconvolution Method Runtimes",
    subtitle = sprintf("Highly Segmented; ANOVA: F = %.2f, p = %.3g", anova_Fstat, anova_pval),
    x = "Method",
    y = "Runtime (seconds)"
  ) + 
  scale_fill_manual(values = c(
    "CITEgeist" = "green",
    "Cell2Location" = "blue",
    "Seurat" = "red",
    "RCTD" = "orange",
    "Tangram" = "purple"
  )) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    panel.grid.major.y = element_line(color = "gray80"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  coord_cartesian(ylim = c(0, y_max))  # Adjust y-axis dynamically

# Add Tukey's significant comparisons only for CITEgeist
if (nrow(tukey_df) > 0) {
  p <- p + geom_signif(
    comparisons = tukey_df$groups,
    annotations = tukey_df$Significance,
    y_position = seq(y_max * 0.9, y_max, length.out = nrow(tukey_df)), 
    tip_length = 0.02
  )
}

# Save plot
ggsave("hs_runtime_comparison.png", plot = p, width = 8, height = 5, dpi = 300)

# Display plot
print(p)
