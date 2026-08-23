#!/usr/bin/env Rscript
# =============================================================================
# STEP 19C0: VERIFY GSE197677 DE INPUTS
# =============================================================================

suppressPackageStartupMessages({
  library(edgeR)
})

base_dir <- "~/Downloads/ESCC Neoadjuvant"

cat("Loading inputs...\n")

raw_counts <- readRDS(file.path(base_dir, "06_tables/cross_dataset/epithelial_pseudobulk/counts/GSE197677_epithelial_pseudobulk_raw_counts.rds"))
meta <- read.csv(file.path(base_dir, "06_tables/cross_dataset/epithelial_pseudobulk/counts/GSE197677_epithelial_pseudobulk_sample_metadata.csv"), stringsAsFactors = FALSE)
dge_tmm <- readRDS(file.path(base_dir, "06_tables/cross_dataset/epithelial_pseudobulk/QC/GSE197677_DGEList_TMM.rds"))
design_meta <- read.csv(file.path(base_dir, "06_tables/cross_dataset/epithelial_pseudobulk/STEP19A_GSE197677_design_matrix.csv"), row.names = 1, check.names = FALSE)
outlier <- read.csv(file.path(base_dir, "06_tables/cross_dataset/epithelial_pseudobulk/QC/STEP19B_sample_outlier_audit.csv"), stringsAsFactors = FALSE)
filter_audit <- read.csv(file.path(base_dir, "06_tables/cross_dataset/epithelial_pseudobulk/QC/STEP19B_expression_filter_audit.csv"), stringsAsFactors = FALSE)

cat("Loaded.\n\n")

errors <- c()

# ---- Sample count ----

cat("=== SAMPLE COUNT CHECKS ===\n\n")

n_total <- ncol(raw_counts)
n_nact <- sum(meta$treatment_standardized == "Neoadjuvant_chemotherapy")
n_nnact <- sum(meta$treatment_standardized == "No_neoadjuvant_chemotherapy")

cat(sprintf("Total samples: %d\n", n_total))
cat(sprintf("NACT: %d\n", n_nact))
cat(sprintf("nNACT: %d\n", n_nnact))

if (n_total != 10) errors <- c(errors, paste("Expected 10 samples, got", n_total))
if (n_nact != 6) errors <- c(errors, paste("Expected 6 NACT, got", n_nact))
if (n_nnact != 4) errors <- c(errors, paste("Expected 4 nNACT, got", n_nnact))

cat("\n")

# ---- Duplicate check ----

cat("=== DUPLICATE SAMPLE ID CHECK ===\n\n")

dup_ids <- colnames(raw_counts)[duplicated(colnames(raw_counts))]
if (length(dup_ids) > 0) {
  errors <- c(errors, paste("Duplicated sample IDs:", paste(dup_ids, collapse = ", ")))
  cat("DUPLICATED:", paste(dup_ids, collapse = ", "), "\n")
} else {
  cat("No duplicated sample IDs.\n")
}

cat("\n")

# ---- QC status check ----

cat("=== QC STATUS CHECK ===\n\n")

outlier_197677 <- outlier[outlier$dataset == "GSE197677", ]
for (i in 1:nrow(outlier_197677)) {
  cat(sprintf("  %s: %s\n", outlier_197677$sample[i], outlier_197677$QC_status[i]))
}

n_multi <- sum(outlier_197677$QC_status == "QC_REVIEW_MULTI_FLAG")
n_low <- sum(outlier_197677$QC_status == "QC_LOW_CELL_PREDEFINED")
cat(sprintf("\nMulti-flag: %d, Low-cell: %d\n", n_multi, n_low))
if (n_multi > 0) errors <- c(errors, "Samples with QC_REVIEW_MULTI_FLAG exist")
if (n_low > 0) errors <- c(errors, "Low-cell samples found in GSE197677")

cat("\n")

# ---- Count integrity ----

cat("=== COUNT INTEGRITY CHECKS ===\n\n")

nas <- sum(is.na(raw_counts))
negs <- sum(raw_counts < 0, na.rm = TRUE)
is_int <- all(raw_counts == floor(raw_counts), na.rm = TRUE)

cat(sprintf("NA counts: %d\n", nas))
cat(sprintf("Negative counts: %d\n", negs))
cat(sprintf("Integer-like: %s\n", is_int))

if (nas > 0) errors <- c(errors, "Has NA counts")
if (negs > 0) errors <- c(errors, "Has negative counts")
if (!is_int) errors <- c(errors, "Counts not integer-like")

cat("\n")

# ---- Design matrix check ----

cat("=== DESIGN MATRIX CHECK ===\n\n")

# Build proper model matrix
treatment_factor <- factor(meta$treatment_standardized, levels = c("No_neoadjuvant_chemotherapy", "Neoadjuvant_chemotherapy"))
design <- model.matrix(~ treatment_factor)
rownames(design) <- meta$sample

cat("Design matrix columns:\n")
print(colnames(design))
cat("\nDesign matrix rank:", qr(design)$rank, "\n")
cat("Expected rank: 2\n")

if (qr(design)$rank != 2) errors <- c(errors, paste("Design matrix rank", qr(design)$rank, "!= 2"))

cat("\n")

# ---- Summary ----

if (length(errors) > 0) {
  cat("========================================\n")
  cat("VERIFICATION FAILED\n")
  cat("========================================\n\n")
  for (e in errors) cat("  ERROR:", e, "\n")
  stop("Verification failed.")
} else {
  cat("========================================\n")
  cat("VERIFICATION PASSED\n")
  cat("========================================\n")
  
  cat("\nCONTRAST DIRECTION:\n")
  cat("  Positive logFC = Higher in NACT (Neoadjuvant_chemotherapy)\n")
  cat("  Negative logFC = Higher in nNACT (No_neoadjuvant_chemotherapy)\n")
  cat("  Reference level: No_neoadjuvant_chemotherapy\n")
  cat("  coef = treatmentNACT\n")
  
  # Save verification
  verification <- list(
    status = "PASS",
    n_total = n_total,
    n_nact = n_nact,
    n_nnact = n_nnact,
    design_columns = colnames(design),
    contrast_direction = "Positive logFC = Higher in NACT"
  )
  saveRDS(verification, file.path(base_dir, "06_tables/cross_dataset/epithelial_pseudobulk/GSE197677_DE/STEP19C0_verification.rds"))
}

cat("\nDone.\n")
