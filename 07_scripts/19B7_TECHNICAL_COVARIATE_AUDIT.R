#!/usr/bin/env Rscript
# =============================================================================
# STEP 19B7: TECHNICAL COVARIATE AUDIT
# =============================================================================

suppressPackageStartupMessages({
  library(edgeR)
})

base_dir <- "~/Downloads/ESCC Neoadjuvant"
qc_dir <- file.path(base_dir, "06_tables/cross_dataset/epithelial_pseudobulk/QC")

cat("Loading data...\n")

dge_221561 <- readRDS(file.path(qc_dir, "GSE221561_DGEList_TMM.rds"))
dge_197677 <- readRDS(file.path(qc_dir, "GSE197677_DGEList_TMM.rds"))
dge_160269 <- readRDS(file.path(qc_dir, "GSE160269_DGEList_TMM.rds"))

logcpm_221561 <- readRDS(file.path(qc_dir, "GSE221561_logCPM_QC.rds"))
logcpm_197677 <- readRDS(file.path(qc_dir, "GSE197677_logCPM_QC.rds"))
logcpm_160269 <- readRDS(file.path(qc_dir, "GSE160269_logCPM_QC.rds"))

pca_221561 <- read.csv(file.path(qc_dir, "STEP19B_GSE221561_PCA_coordinates.csv"), stringsAsFactors = FALSE)
pca_197677 <- read.csv(file.path(qc_dir, "STEP19B_GSE197677_PCA_coordinates.csv"), stringsAsFactors = FALSE)
pca_160269 <- read.csv(file.path(qc_dir, "STEP19B_GSE160269_PCA_coordinates.csv"), stringsAsFactors = FALSE)

cat("Loaded.\n\n")

# ---- Helper function ----

interpret_rho <- function(rho) {
  abs_rho <- abs(rho)
  if (abs_rho < 0.4) return("WEAK")
  if (abs_rho < 0.7) return("MODERATE")
  return("STRONG")
}

# ---- Build technical covariate table for each dataset ----

build_tech_table <- function(dge, logcpm, pca_coords, dataset) {
  samples <- dge$samples$sample
  n <- length(samples)
  
  data.frame(
    dataset = dataset,
    sample = samples,
    treatment = dge$samples$treatment_standardized,
    epithelial_cells = dge$samples$epithelial_cells,
    library_size = dge$samples$lib.size,
    effective_library_size = dge$samples$lib.size * dge$samples$norm.factors,
    genes_detected = colSums(logcpm > 0)[samples],
    TMM_factor = dge$samples$norm.factors,
    PC1 = pca_coords$PC1[match(samples, pca_coords$sample)],
    PC2 = pca_coords$PC2[match(samples, pca_coords$sample)],
    PC3 = pca_coords$PC3[match(samples, pca_coords$sample)],
    stringsAsFactors = FALSE
  )
}

tech_221561 <- build_tech_table(dge_221561, logcpm_221561, pca_221561, "GSE221561")
tech_197677 <- build_tech_table(dge_197677, logcpm_197677, pca_197677, "GSE197677")
tech_160269 <- build_tech_table(dge_160269, logcpm_160269, pca_160269, "GSE160269")

# ---- Calculate Spearman correlations ----

cat("=== TECHNICAL COVARIATE CORRELATIONS ===\n\n")

calc_correlations <- function(tech, dataset) {
  vars <- c("epithelial_cells", "library_size", "genes_detected")
  pcs <- c("PC1", "PC2", "PC3")
  
  results <- list()
  
  for (v in vars) {
    for (pc in pcs) {
      rho <- cor(tech[[v]], tech[[pc]], method = "spearman")
      results[[length(results) + 1]] <- data.frame(
        dataset = dataset,
        variable_x = v,
        variable_y = pc,
        spearman_rho = round(rho, 4),
        interpretation = interpret_rho(rho),
        stringsAsFactors = FALSE
      )
    }
  }
  
  # Also: epithelial_cells vs library_size, genes_detected
  rho_ec_lib <- cor(tech$epithelial_cells, tech$library_size, method = "spearman")
  results[[length(results) + 1]] <- data.frame(
    dataset = dataset,
    variable_x = "epithelial_cells",
    variable_y = "library_size",
    spearman_rho = round(rho_ec_lib, 4),
    interpretation = interpret_rho(rho_ec_lib),
    stringsAsFactors = FALSE
  )
  
  rho_ec_genes <- cor(tech$epithelial_cells, tech$genes_detected, method = "spearman")
  results[[length(results) + 1]] <- data.frame(
    dataset = dataset,
    variable_x = "epithelial_cells",
    variable_y = "genes_detected",
    spearman_rho = round(rho_ec_genes, 4),
    interpretation = interpret_rho(rho_ec_genes),
    stringsAsFactors = FALSE
  )
  
  do.call(rbind, results)
}

cor_221561 <- calc_correlations(tech_221561, "GSE221561")
cor_197677 <- calc_correlations(tech_197677, "GSE197677")
cor_160269 <- calc_correlations(tech_160269, "GSE160269")

cor_all <- rbind(cor_221561, cor_197677, cor_160269)

# Print notable correlations
cat("Notable correlations (|rho| >= 0.4):\n\n")
notable <- cor_all[abs(cor_all$spearman_rho) >= 0.4, ]
if (nrow(notable) > 0) {
  for (i in 1:nrow(notable)) {
    cat(sprintf("  %s: %s vs %s: rho=%.3f (%s)\n",
                notable$dataset[i], notable$variable_x[i], notable$variable_y[i],
                notable$spearman_rho[i], notable$interpretation[i]))
  }
} else {
  cat("  None found.\n")
}

cat("\n")

# Save
write.csv(cor_all, file.path(qc_dir, "STEP19B_technical_covariate_correlations.csv"), row.names = FALSE)
cat("Saved: STEP19B_technical_covariate_correlations.csv\n")

cat("\nDone.\n")
