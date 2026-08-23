#!/usr/bin/env Rscript
# =============================================================================
# STEP 19B6: SAMPLE CORRELATION
# =============================================================================

suppressPackageStartupMessages({
  library(edgeR)
})

base_dir <- "~/Downloads/ESCC Neoadjuvant"
qc_dir <- file.path(base_dir, "06_tables/cross_dataset/epithelial_pseudobulk/QC")
fig_dir <- file.path(base_dir, "05_figures/cross_dataset/epithelial_pseudobulk/QC")

cat("Loading logCPM matrices and filter vectors...\n")

logcpm_221561 <- readRDS(file.path(qc_dir, "GSE221561_logCPM_QC.rds"))
logcpm_197677 <- readRDS(file.path(qc_dir, "GSE197677_logCPM_QC.rds"))
logcpm_160269 <- readRDS(file.path(qc_dir, "GSE160269_logCPM_QC.rds"))

keep_221561 <- readRDS(file.path(qc_dir, "GSE221561_filter_keep.rds"))
keep_197677 <- readRDS(file.path(qc_dir, "GSE197677_filter_keep.rds"))
keep_160269 <- readRDS(file.path(qc_dir, "GSE160269_filter_keep.rds"))

cat("Loaded.\n\n")

# ---- Correlation function ----

run_correlation <- function(logcpm, keep, dataset_name, fig_path, csv_path) {
  cat(sprintf("=== %s SAMPLE CORRELATION ===\n\n", dataset_name))
  
  # Filter to genes passing filter
  logcpm_filtered <- logcpm[keep, ]
  cat(sprintf("  Genes used: %d\n", nrow(logcpm_filtered)))
  
  # Calculate Spearman correlation
  cor_mat <- cor(logcpm_filtered, method = "spearman", use = "pairwise.complete.obs")
  
  cat(sprintf("  Correlation matrix: %d x %d\n", nrow(cor_mat), ncol(cor_mat)))
  cat(sprintf("  Off-diagonal range: [%.3f, %.3f]\n",
              min(cor_mat[lower.tri(cor_mat)]),
              max(cor_mat[lower.tri(cor_mat)])))
  
  # Save CSV
  write.csv(cor_mat, csv_path)
  cat(sprintf("  Saved: %s\n", basename(csv_path)))
  
  # Plot heatmap
  png(fig_path, width = max(7, ncol(cor_mat) * 0.6), height = max(6, nrow(cor_mat) * 0.5),
      units = "in", res = 150)
  
  # Order by hierarchical clustering for visualization
  hc <- hclust(as.dist(1 - cor_mat))
  cor_ordered <- cor_mat[hc$order, hc$order]
  
  # Heatmap with color scale
  n <- nrow(cor_ordered)
  image(1:n, 1:n, cor_ordered,
        col = colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100),
        zlim = c(min(cor_ordered), 1),
        xlab = "", ylab = "",
        main = sprintf("%s Sample-Sample Spearman Correlation", dataset_name),
        cex.main = 1.2,
        axes = FALSE)
  
  axis(1, at = 1:n, labels = colnames(cor_ordered), las = 2, cex.axis = 0.7)
  axis(2, at = 1:n, labels = rownames(cor_ordered), las = 2, cex.axis = 0.7)
  
  # Add correlation values
  for (i in 1:n) {
    for (j in 1:n) {
      text(i, j, sprintf("%.2f", cor_ordered[i, j]), cex = 0.5)
    }
  }
  
  dev.off()
  cat(sprintf("  Saved: %s\n\n", basename(fig_path)))
  
  return(cor_mat)
}

# ---- GSE221561 correlation ----

cor_221561 <- run_correlation(logcpm_221561, keep_221561, "GSE221561",
                              file.path(fig_dir, "STEP19B_GSE221561_sample_correlation_heatmap.png"),
                              file.path(qc_dir, "STEP19B_GSE221561_sample_correlation.csv"))

# ---- GSE197677 correlation ----

cor_197677 <- run_correlation(logcpm_197677, keep_197677, "GSE197677",
                              file.path(fig_dir, "STEP19B_GSE197677_sample_correlation_heatmap.png"),
                              file.path(qc_dir, "STEP19B_GSE197677_sample_correlation.csv"))

# ---- GSE160269 correlation ----

cor_160269 <- run_correlation(logcpm_160269, keep_160269, "GSE160269",
                              file.path(fig_dir, "STEP19B_GSE160269_sample_correlation_heatmap.png"),
                              file.path(qc_dir, "STEP19B_GSE160269_sample_correlation.csv"))

cat("=== ALL SAMPLE CORRELATION ANALYSES COMPLETE ===\n")
