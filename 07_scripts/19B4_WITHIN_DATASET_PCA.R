#!/usr/bin/env Rscript
# =============================================================================
# STEP 19B4: WITHIN-DATASET PCA
# =============================================================================

suppressPackageStartupMessages({
  library(edgeR)
})

base_dir <- "~/Downloads/ESCC Neoadjuvant"
qc_dir <- file.path(base_dir, "06_tables/cross_dataset/epithelial_pseudobulk/QC")
fig_dir <- file.path(base_dir, "05_figures/cross_dataset/epithelial_pseudobulk/QC")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

cat("Loading data...\n")

dge_221561 <- readRDS(file.path(qc_dir, "GSE221561_DGEList_TMM.rds"))
dge_197677 <- readRDS(file.path(qc_dir, "GSE197677_DGEList_TMM.rds"))
dge_160269 <- readRDS(file.path(qc_dir, "GSE160269_DGEList_TMM.rds"))

logcpm_221561 <- readRDS(file.path(qc_dir, "GSE221561_logCPM_QC.rds"))
logcpm_197677 <- readRDS(file.path(qc_dir, "GSE197677_logCPM_QC.rds"))
logcpm_160269 <- readRDS(file.path(qc_dir, "GSE160269_logCPM_QC.rds"))

keep_221561 <- readRDS(file.path(qc_dir, "GSE221561_filter_keep.rds"))
keep_197677 <- readRDS(file.path(qc_dir, "GSE197677_filter_keep.rds"))
keep_160269 <- readRDS(file.path(qc_dir, "GSE160269_filter_keep.rds"))

cat("Loaded.\n\n")

# ---- Helper function for PCA ----

run_pca_and_plot <- function(logcpm, keep, meta, dataset_name, treatment_label, color_values, fig_path, coord_path) {
  cat(sprintf("=== %s WITHIN-DATASET PCA ===\n\n", dataset_name))
  
  # Filter to genes passing filter
  logcpm_filtered <- logcpm[keep, ]
  cat(sprintf("  Genes used: %d\n", nrow(logcpm_filtered)))
  
  # Select top 2000 most variable genes (or all if fewer)
  gene_vars <- apply(logcpm_filtered, 1, var)
  n_select <- min(2000, nrow(logcpm_filtered))
  top_genes <- names(sort(gene_vars, decreasing = TRUE))[1:n_select]
  logcpm_top <- logcpm_filtered[top_genes, ]
  cat(sprintf("  Genes for PCA: %d (top %d by variance)\n", nrow(logcpm_top), n_select))
  
  # PCA (center = TRUE, scale = FALSE)
  pca_res <- prcomp(t(logcpm_top), center = TRUE, scale. = FALSE)
  
  # Variance explained
  var_explained <- summary(pca_res)$importance[2, ] * 100
  cat(sprintf("  PC1 variance: %.1f%%\n", var_explained[1]))
  cat(sprintf("  PC2 variance: %.1f%%\n", var_explained[2]))
  cat(sprintf("  PC3 variance: %.1f%%\n\n", var_explained[3]))
  
  # Build coordinates table
  coords <- data.frame(
    dataset = dataset_name,
    sample = rownames(pca_res$x),
    treatment = meta$treatment_standardized[match(rownames(pca_res$x), meta$sample)],
    PC1 = pca_res$x[, 1],
    PC2 = pca_res$x[, 2],
    PC3 = pca_res$x[, 3],
    PC4 = pca_res$x[, 4],
    PC5 = pca_res$x[, 5],
    epithelial_cells = meta$epithelial_cells[match(rownames(pca_res$x), meta$sample)],
    library_size = meta$lib.size[match(rownames(pca_res$x), meta$sample)],
    genes_detected = colSums(logcpm > 0)[match(rownames(pca_res$x), names(colSums(logcpm > 0)))],
    stringsAsFactors = FALSE
  )
  
  write.csv(coords, coord_path, row.names = FALSE)
  cat(sprintf("  Saved: %s\n", basename(coord_path)))
  
  # Plot
  png(fig_path, width = 7, height = 5, units = "in", res = 150)
  
  treatment <- coords$treatment
  unique_treat <- unique(treatment)
  
  # Default colors if not provided
  if (is.null(color_values)) {
    cols <- RColorBrewer::brewer.pal(max(3, length(unique_treat)), "Set2")[1:length(unique_treat)]
    names(cols) <- unique_treat
  } else {
    cols <- color_values
  }
  
  plot(coords$PC1, coords$PC2,
       col = cols[treatment], pch = 19, cex = 1.5,
       xlab = sprintf("PC1 (%.1f%% variance)", var_explained[1]),
       ylab = sprintf("PC2 (%.1f%% variance)", var_explained[2]),
       main = sprintf("%s Within-Dataset PCA", dataset_name),
       cex.main = 1.2)
  
  # Add sample labels
  text(coords$PC1, coords$PC2, labels = coords$sample, cex = 0.6, pos = 3, offset = 0.3)
  
  # Legend
  legend("topright", legend = names(cols), col = cols, pch = 19, cex = 0.8,
         title = treatment_label)
  
  dev.off()
  cat(sprintf("  Saved: %s\n\n", basename(fig_path)))
  
  return(coords)
}

# ---- GSE221561 PCA ----

coords_221561 <- run_pca_and_plot(
  logcpm = logcpm_221561,
  keep = keep_221561,
  meta = dge_221561$samples,
  dataset_name = "GSE221561",
  treatment_label = "Treatment",
  color_values = c("Neoadjuvant_treated" = "#E41A1C", "Surgery_alone" = "#377EB8"),
  fig_path = file.path(fig_dir, "STEP19B_GSE221561_PCA.png"),
  coord_path = file.path(qc_dir, "STEP19B_GSE221561_PCA_coordinates.csv")
)

# ---- GSE197677 PCA ----

coords_197677 <- run_pca_and_plot(
  logcpm = logcpm_197677,
  keep = keep_197677,
  meta = dge_197677$samples,
  dataset_name = "GSE197677",
  treatment_label = "Treatment",
  color_values = c("Neoadjuvant_chemotherapy" = "#E41A1C", "No_neoadjuvant_chemotherapy" = "#377EB8"),
  fig_path = file.path(fig_dir, "STEP19B_GSE197677_PCA.png"),
  coord_path = file.path(qc_dir, "STEP19B_GSE197677_PCA_coordinates.csv")
)

# ---- GSE160269 PCA ----

coords_160269 <- run_pca_and_plot(
  logcpm = logcpm_160269,
  keep = keep_160269,
  meta = dge_160269$samples,
  dataset_name = "GSE160269",
  treatment_label = "Treatment Group",
  color_values = NULL,  # use default
  fig_path = file.path(fig_dir, "STEP19B_GSE160269_PCA.png"),
  coord_path = file.path(qc_dir, "STEP19B_GSE160269_PCA_coordinates.csv")
)

cat("=== ALL WITHIN-DATASET PCA COMPLETE ===\n")
