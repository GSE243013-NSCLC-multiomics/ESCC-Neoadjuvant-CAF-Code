#!/usr/bin/env Rscript
# =============================================================================
# STEP 19C12: LEAVE-ONE-SAMPLE-OUT SENSITIVITY
# =============================================================================

suppressPackageStartupMessages({
  library(edgeR)
})

base_dir <- "~/Downloads/ESCC Neoadjuvant"
de_dir <- file.path(base_dir, "06_tables/cross_dataset/epithelial_pseudobulk/GSE197677_DE")

cat("Loading data...\n")

raw_counts <- readRDS(file.path(base_dir, "06_tables/cross_dataset/epithelial_pseudobulk/counts/GSE197677_epithelial_pseudobulk_raw_counts.rds"))
meta <- read.csv(file.path(base_dir, "06_tables/cross_dataset/epithelial_pseudobulk/counts/GSE197677_epithelial_pseudobulk_sample_metadata.csv"), stringsAsFactors = FALSE)
res_full <- read.csv(file.path(de_dir, "STEP19C_GSE197677_NACT_vs_nNACT_ALL_GENES.csv"), row.names = 1, check.names = FALSE)

cat("Loaded.\n\n")

# Top 100 genes from full model
n_top <- min(100, nrow(res_full))
top100_genes <- rownames(res_full)[1:n_top]
cat(sprintf("Top %d genes from full model.\n\n", n_top))

# ---- Run LOSO ----

cat("=== LEAVE-ONE-OUT ANALYSIS ===\n\n")

samples <- meta$sample
n_samples <- length(samples)

# Store results: for each gene, track direction across LOSO runs
loso_results <- data.frame(
  gene = top100_genes,
  full_logFC = res_full$logFC[1:n_top],
  full_FDR = res_full$FDR[1:n_top],
  runs_tested = 0,
  same_direction_runs = 0,
  LOSO_logFC_sum = 0,
  LOSO_logFC_sq_sum = 0,
  LOSO_logFC_min = Inf,
  LOSO_logFC_max = -Inf,
  stringsAsFactors = FALSE
)
rownames(loso_results) <- top100_genes

for (i in 1:n_samples) {
  leave_out <- samples[i]
  cat(sprintf("LOSO run %d/%d: leaving out %s\n", i, n_samples, leave_out))
  
  # Subset data
  counts_lo <- raw_counts[, -i, drop = FALSE]
  meta_lo <- meta[-i, ]
  
  # Skip if only one group has samples
  if (length(unique(meta_lo$treatment_standardized)) < 2) {
    cat("  WARNING: Only one treatment group after leaving out sample. Skipping.\n")
    next
  }
  
  # Build DGEList
  dge_lo <- DGEList(counts = counts_lo, samples = meta_lo)
  
  # Filter
  treatment_lo <- factor(meta_lo$treatment_standardized,
                         levels = c("No_neoadjuvant_chemotherapy", "Neoadjuvant_chemotherapy"))
  design_lo <- model.matrix(~ treatment_lo)
  rownames(design_lo) <- meta_lo$sample
  
  # Find the NACT coefficient name (might vary if factor levels drop)
  nact_coef <- grep("Neoadjuvant_chemotherapy", colnames(design_lo), value = TRUE)
  if (length(nact_coef) == 0) {
    cat("  WARNING: NACT coefficient not found. Skipping.\n")
    next
  }
  nact_coef <- nact_coef[1]  # take first match
  
  keep_lo <- filterByExpr(dge_lo, design = design_lo)
  dge_lo <- dge_lo[keep_lo, , keep.lib.sizes = FALSE]
  dge_lo <- calcNormFactors(dge_lo, method = "TMM")
  
  # Dispersion
  dge_lo <- estimateDisp(dge_lo, design_lo, robust = TRUE)
  
  # Fit
  fit_lo <- glmQLFit(dge_lo, design_lo, robust = TRUE)
  qlf_lo <- glmQLFTest(fit_lo, coef = nact_coef)
  
  # Extract results for top 100 genes
  res_lo <- topTags(qlf_lo, n = Inf, sort.by = "none")$table
  
  for (gene in top100_genes) {
    if (gene %in% rownames(res_lo)) {
      loso_results[gene, "runs_tested"] <- loso_results[gene, "runs_tested"] + 1
      
      logfc_lo <- res_lo[gene, "logFC"]
      
      # Check direction
      same_dir <- sign(logfc_lo) == sign(loso_results[gene, "full_logFC"])
      if (same_dir) {
        loso_results[gene, "same_direction_runs"] <- loso_results[gene, "same_direction_runs"] + 1
      }
      
      # Accumulate statistics
      loso_results[gene, "LOSO_logFC_sum"] <- loso_results[gene, "LOSO_logFC_sum"] + logfc_lo
      loso_results[gene, "LOSO_logFC_sq_sum"] <- loso_results[gene, "LOSO_logFC_sq_sum"] + logfc_lo^2
      loso_results[gene, "LOSO_logFC_min"] <- min(loso_results[gene, "LOSO_logFC_min"], logfc_lo)
      loso_results[gene, "LOSO_logFC_max"] <- max(loso_results[gene, "LOSO_logFC_max"], logfc_lo)
    }
  }
}

cat("\n")

# ---- Calculate summary statistics ----

cat("=== CALCULATING SUMMARY STATISTICS ===\n\n")

# Direction consistency
loso_results$direction_consistency_fraction <- loso_results$same_direction_runs / loso_results$runs_tested

# Median logFC across LOSO runs
loso_results$LOSO_median_logFC <- loso_results$LOSO_logFC_sum / loso_results$runs_tested

# Replace Inf/Min/Max with NA if no runs
loso_results$LOSO_logFC_min[loso_results$LOSO_logFC_min == Inf] <- NA
loso_results$LOSO_logFC_max[loso_results$LOSO_logFC_max == -Inf] <- NA

# Classify stability
loso_results$stability <- NA
loso_results$stability[loso_results$direction_consistency_fraction >= 0.9] <- "ROBUST_DIRECTION"
loso_results$stability[loso_results$direction_consistency_fraction >= 0.7 &
                       loso_results$direction_consistency_fraction < 0.9] <- "MODERATE_DIRECTION"
loso_results$stability[loso_results$direction_consistency_fraction < 0.7] <- "UNSTABLE_DIRECTION"

# Print summary
cat("Stability classification:\n")
print(table(loso_results$stability))

cat("\n")

# ---- Save ----

output <- loso_results[, c("gene", "full_logFC", "full_FDR", "runs_tested",
                            "same_direction_runs", "direction_consistency_fraction",
                            "LOSO_median_logFC", "LOSO_logFC_min", "LOSO_logFC_max", "stability")]

write.csv(output, file.path(de_dir, "STEP19C_GSE197677_LOSO_top100_stability.csv"), row.names = FALSE)
cat("Saved: STEP19C_GSE197677_LOSO_top100_stability.csv\n\n")

# Print top 10 stable genes
cat("Top 10 genes by direction consistency:\n")
top_stable <- output[order(-output$direction_consistency_fraction), ][1:10, ]
for (i in 1:nrow(top_stable)) {
  cat(sprintf("  %s: consistency=%.2f, median_logFC=%.3f, stability=%s\n",
              top_stable$gene[i], top_stable$direction_consistency_fraction[i],
              top_stable$LOSO_median_logFC[i], top_stable$stability[i]))
}

cat("\n=== LOSO SENSITIVITY ANALYSIS COMPLETE ===\n")
