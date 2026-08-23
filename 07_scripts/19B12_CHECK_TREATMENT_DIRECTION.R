#!/usr/bin/env Rscript
# =============================================================================
# STEP 19B12: CHECK TREATMENT DIRECTION WITHOUT TESTING
# =============================================================================

suppressPackageStartupMessages({
  library(edgeR)
})

base_dir <- "~/Downloads/ESCC Neoadjuvant"
qc_dir <- file.path(base_dir, "06_tables/cross_dataset/epithelial_pseudobulk/QC")

cat("Loading PCA coordinates...\n")

pca_221561 <- read.csv(file.path(qc_dir, "STEP19B_GSE221561_PCA_coordinates.csv"), stringsAsFactors = FALSE)
pca_197677 <- read.csv(file.path(qc_dir, "STEP19B_GSE197677_PCA_coordinates.csv"), stringsAsFactors = FALSE)

cat("Loaded.\n\n")

# ---- Calculate group centroids ----

cat("=== PCA GROUP CENTROIDS (DESCRIPTIVE ONLY) ===\n\n")

calc_centroids <- function(pca_coords, dataset) {
  groups <- unique(pca_coords$treatment)
  
  results <- list()
  for (g in groups) {
    sub <- pca_coords[pca_coords$treatment == g, ]
    results[[g]] <- data.frame(
      dataset = dataset,
      treatment = g,
      n_samples = nrow(sub),
      PC1_mean = mean(sub$PC1),
      PC1_sd = sd(sub$PC1),
      PC2_mean = mean(sub$PC2),
      PC2_sd = sd(sub$PC2),
      PC3_mean = mean(sub$PC3),
      PC3_sd = sd(sub$PC3),
      stringsAsFactors = FALSE
    )
  }
  
  res <- do.call(rbind, results)
  
  # Print
  cat(sprintf("--- %s ---\n", dataset))
  for (g in groups) {
    r <- res[res$treatment == g, ]
    cat(sprintf("  %s (n=%d):\n", g, r$n_samples))
    cat(sprintf("    PC1: %.1f (+/- %.1f)\n", r$PC1_mean, r$PC1_sd))
    cat(sprintf("    PC2: %.1f (+/- %.1f)\n", r$PC2_mean, r$PC2_sd))
  }
  cat("\n")
  
  return(res)
}

centroids_221561 <- calc_centroids(pca_221561, "GSE221561")
centroids_197677 <- calc_centroids(pca_197677, "GSE197677")

# ---- Combine and save ----

centroids_all <- rbind(centroids_221561, centroids_197677)
centroids_all$annotation <- "DESCRIPTIVE_ONLY"

write.csv(centroids_all, file.path(qc_dir, "STEP19B_PCA_group_centroids_descriptive.csv"), row.names = FALSE)
cat("Saved: STEP19B_PCA_group_centroids_descriptive.csv\n")

cat("\nNOTE: These are descriptive summaries only.\n")
cat("Do NOT interpret as treatment effects.\n")

cat("\nDone.\n")
