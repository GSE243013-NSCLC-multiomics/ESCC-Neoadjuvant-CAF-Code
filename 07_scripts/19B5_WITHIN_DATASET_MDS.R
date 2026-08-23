#!/usr/bin/env Rscript
# =============================================================================
# STEP 19B5: WITHIN-DATASET MDS
# =============================================================================

suppressPackageStartupMessages({
  library(edgeR)
})

base_dir <- "~/Downloads/ESCC Neoadjuvant"
qc_dir <- file.path(base_dir, "06_tables/cross_dataset/epithelial_pseudobulk/QC")
fig_dir <- file.path(base_dir, "05_figures/cross_dataset/epithelial_pseudobulk/QC")

cat("Loading DGEList objects and filter vectors...\n")

dge_221561 <- readRDS(file.path(qc_dir, "GSE221561_DGEList_TMM.rds"))
dge_197677 <- readRDS(file.path(qc_dir, "GSE197677_DGEList_TMM.rds"))
dge_160269 <- readRDS(file.path(qc_dir, "GSE160269_DGEList_TMM.rds"))

keep_221561 <- readRDS(file.path(qc_dir, "GSE221561_filter_keep.rds"))
keep_197677 <- readRDS(file.path(qc_dir, "GSE197677_filter_keep.rds"))
keep_160269 <- readRDS(file.path(qc_dir, "GSE160269_filter_keep.rds"))

cat("Loaded.\n\n")

# ---- MDS function ----

run_mds_plot <- function(dge, keep, dataset_name, treatment_label, color_values, fig_path) {
  cat(sprintf("=== %s MDS PLOT ===\n\n", dataset_name))
  
  # Subset to filtered genes
  dge_filtered <- dge[keep, ]
  
  # Calculate MDS
  mds <- plotMDS(dge_filtered, plot = FALSE)
  
  cat(sprintf("  MDS calculated for %d genes\n", nrow(dge_filtered$counts)))
  cat(sprintf("  Top 2 dimensions explained: %.1f%% and %.1f%% of variation\n",
              mds$var.explained[1] * 100, mds$var.explained[2] * 100))
  
  # Get treatment labels
  treatment <- dge$samples$treatment_standardized
  unique_treat <- unique(treatment)
  
  if (is.null(color_values)) {
    cols <- RColorBrewer::brewer.pal(max(3, length(unique_treat)), "Set2")[1:length(unique_treat)]
    names(cols) <- unique_treat
  } else {
    cols <- color_values
  }
  
  # Plot
  png(fig_path, width = 7, height = 5, units = "in", res = 150)
  
  plot(mds$x, mds$y,
       col = cols[treatment], pch = 19, cex = 1.5,
       xlab = sprintf("Leading logFC dimension 1 (%.1f%%)", mds$var.explained[1] * 100),
       ylab = sprintf("Leading logFC dimension 2 (%.1f%%)", mds$var.explained[2] * 100),
       main = sprintf("%s MDS Plot", dataset_name),
       cex.main = 1.2)
  
  text(mds$x, mds$y, labels = colnames(dge$counts), cex = 0.6, pos = 3, offset = 0.3)
  
  legend("topright", legend = names(cols), col = cols, pch = 19, cex = 0.8,
         title = treatment_label)
  
  dev.off()
  cat(sprintf("  Saved: %s\n\n", basename(fig_path)))
}

# ---- GSE221561 MDS ----

run_mds_plot(dge_221561, keep_221561, "GSE221561", "Treatment",
             c("Neoadjuvant_treated" = "#E41A1C", "Surgery_alone" = "#377EB8"),
             file.path(fig_dir, "STEP19B_GSE221561_MDS.png"))

# ---- GSE197677 MDS ----

run_mds_plot(dge_197677, keep_197677, "GSE197677", "Treatment",
             c("Neoadjuvant_chemotherapy" = "#E41A1C", "No_neoadjuvant_chemotherapy" = "#377EB8"),
             file.path(fig_dir, "STEP19B_GSE197677_MDS.png"))

# ---- GSE160269 MDS ----

run_mds_plot(dge_160269, keep_160269, "GSE160269", "Treatment Group",
             NULL,
             file.path(fig_dir, "STEP19B_GSE160269_MDS.png"))

cat("=== ALL MDS PLOTS COMPLETE ===\n")
