#!/usr/bin/env Rscript
# =============================================================================
# STEP 19B9: TREATMENT-vs-TECHNICAL BALANCE AUDIT
# =============================================================================

suppressPackageStartupMessages({
  library(edgeR)
})

base_dir <- "~/Downloads/ESCC Neoadjuvant"
qc_dir <- file.path(base_dir, "06_tables/cross_dataset/epithelial_pseudobulk/QC")

cat("Loading data...\n")

dge_221561 <- readRDS(file.path(qc_dir, "GSE221561_DGEList_TMM.rds"))
dge_197677 <- readRDS(file.path(qc_dir, "GSE197677_DGEList_TMM.rds"))

cat("Loaded.\n\n")

# ---- Helper function ----

summarize_group <- function(dge, dataset) {
  samples <- dge$samples
  groups <- unique(samples$treatment_standardized)
  
  results <- list()
  for (g in groups) {
    sub <- samples[samples$treatment_standardized == g, ]
    results[[g]] <- data.frame(
      dataset = dataset,
      treatment = g,
      n_samples = nrow(sub),
      median_epithelial_cells = median(sub$epithelial_cells),
      range_epithelial_cells = paste0(range(sub$epithelial_cells), collapse = "-"),
      median_library_size = median(sub$lib.size),
      range_library_size = paste0(format(range(sub$lib.size), big.mark = ","), collapse = "-"),
      median_genes_detected = NA,  # will fill from logcpm
      range_genes_detected = NA,
      median_TMM_factor = median(sub$norm.factors),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, results)
}

# ---- GSE197677 ----

cat("=== GSE197677 TREATMENT BALANCE ===\n\n")

logcpm_197677 <- readRDS(file.path(qc_dir, "GSE197677_logCPM_QC.rds"))
samples_197677 <- dge_197677$samples$sample

# Add genes_detected
for (s in samples_197677) {
  dge_197677$samples$genes_detected[dge_197677$samples$sample == s] <- sum(logcpm_197677[, s] > 0)
}

groups_197677 <- unique(dge_197677$samples$treatment_standardized)
balance_197677 <- list()
for (g in groups_197677) {
  sub <- dge_197677$samples[dge_197677$samples$treatment_standardized == g, ]
  balance_197677[[g]] <- data.frame(
    dataset = "GSE197677",
    treatment = g,
    n_samples = nrow(sub),
    median_epithelial_cells = median(sub$epithelial_cells),
    range_epithelial_cells = paste0(range(sub$epithelial_cells), collapse = "-"),
    median_library_size = median(sub$lib.size),
    range_library_size = paste0(format(range(sub$lib.size), big.mark = ","), collapse = "-"),
    median_genes_detected = median(sub$genes_detected),
    range_genes_detected = paste0(range(sub$genes_detected), collapse = "-"),
    median_TMM_factor = median(sub$norm.factors),
    stringsAsFactors = FALSE
  )
  cat(sprintf("  %s (n=%d): epi_cells=%s, lib_size=%s, genes=%s\n",
              g, nrow(sub),
              paste0(range(sub$epithelial_cells), collapse = "-"),
              paste0(format(range(sub$lib.size), big.mark = ","), collapse = "-"),
              paste0(range(sub$genes_detected), collapse = "-")))
}

cat("\n")

# ---- GSE221561 ----

cat("=== GSE221561 TREATMENT BALANCE ===\n\n")

logcpm_221561 <- readRDS(file.path(qc_dir, "GSE221561_logCPM_QC.rds"))
samples_221561 <- dge_221561$samples$sample

for (s in samples_221561) {
  dge_221561$samples$genes_detected[dge_221561$samples$sample == s] <- sum(logcpm_221561[, s] > 0)
}

groups_221561 <- unique(dge_221561$samples$treatment_standardized)
balance_221561 <- list()
for (g in groups_221561) {
  sub <- dge_221561$samples[dge_221561$samples$treatment_standardized == g, ]
  n <- nrow(sub)
  
  if (n == 2) {
    cat(sprintf("  NOTE: SMALL_REFERENCE_GROUP for %s (n=%d)\n", g, n))
  }
  
  balance_221561[[g]] <- data.frame(
    dataset = "GSE221561",
    treatment = g,
    n_samples = n,
    median_epithelial_cells = median(sub$epithelial_cells),
    range_epithelial_cells = paste0(range(sub$epithelial_cells), collapse = "-"),
    median_library_size = median(sub$lib.size),
    range_library_size = paste0(format(range(sub$lib.size), big.mark = ","), collapse = "-"),
    median_genes_detected = median(sub$genes_detected),
    range_genes_detected = paste0(range(sub$genes_detected), collapse = "-"),
    median_TMM_factor = median(sub$norm.factors),
    stringsAsFactors = FALSE
  )
  cat(sprintf("  %s (n=%d): epi_cells=%s, lib_size=%s, genes=%s\n",
              g, n,
              paste0(range(sub$epithelial_cells), collapse = "-"),
              paste0(format(range(sub$lib.size), big.mark = ","), collapse = "-"),
              paste0(range(sub$genes_detected), collapse = "-")))
}

cat("\n")

# ---- Combine and save ----

balance_all <- rbind(do.call(rbind, balance_197677), do.call(rbind, balance_221561))

write.csv(balance_all, file.path(qc_dir, "STEP19B_group_technical_balance.csv"), row.names = FALSE)
cat("Saved: STEP19B_group_technical_balance.csv\n")

cat("\nDone.\n")
