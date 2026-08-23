#!/usr/bin/env Rscript
# =============================================================================
# STEP 19B0: VERIFY STEP 19A INPUTS
# =============================================================================

suppressPackageStartupMessages({
  library(edgeR)
})

base_dir <- "~/Downloads/ESCC Neoadjuvant"

# ---- Load Step 19A inputs ----

cat("Loading pseudobulk counts and metadata...\n")

counts_221561 <- readRDS(file.path(base_dir, "06_tables/cross_dataset/epithelial_pseudobulk/counts/GSE221561_epithelial_pseudobulk_raw_counts.rds"))
counts_197677 <- readRDS(file.path(base_dir, "06_tables/cross_dataset/epithelial_pseudobulk/counts/GSE197677_epithelial_pseudobulk_raw_counts.rds"))
counts_160269 <- readRDS(file.path(base_dir, "06_tables/cross_dataset/epithelial_pseudobulk/counts/GSE160269_epithelial_pseudobulk_raw_counts.rds"))

meta_221561 <- read.csv(file.path(base_dir, "06_tables/cross_dataset/epithelial_pseudobulk/counts/GSE221561_epithelial_pseudobulk_sample_metadata.csv"), stringsAsFactors = FALSE)
meta_197677 <- read.csv(file.path(base_dir, "06_tables/cross_dataset/epithelial_pseudobulk/counts/GSE197677_epithelial_pseudobulk_sample_metadata.csv"), stringsAsFactors = FALSE)
meta_160269 <- read.csv(file.path(base_dir, "06_tables/cross_dataset/epithelial_pseudobulk/counts/GSE160269_epithelial_pseudobulk_sample_metadata.csv"), stringsAsFactors = FALSE)

common_genes <- scan(file.path(base_dir, "06_tables/cross_dataset/epithelial_pseudobulk/STEP19A_three_dataset_common_gene_symbols.txt"), character(), quiet = TRUE)
qc_table <- read.csv(file.path(base_dir, "06_tables/cross_dataset/epithelial_pseudobulk/STEP19A_pseudobulk_sample_QC.csv"), stringsAsFactors = FALSE)

cat("Loaded.\n\n")

# ---- Dimension checks ----

errors <- c()

cat("=== DIMENSION CHECKS ===\n\n")

cat("GSE221561: dim =", paste(dim(counts_221561), collapse = " x "), "\n")
if (ncol(counts_221561) != 9) errors <- c(errors, "GSE221561: expected 9 columns, got {ncol(counts_221561)}")

cat("GSE197677: dim =", paste(dim(counts_197677), collapse = " x "), "\n")
if (ncol(counts_197677) != 10) errors <- c(errors, "GSE197677: expected 10 tumor columns, got {ncol(counts_197677)}")

cat("GSE160269: dim =", paste(dim(counts_160269), collapse = " x "), "\n")
if (ncol(counts_160269) != 60) errors <- c(errors, "GSE160269: expected 60 tumor columns, got {ncol(counts_160269)}")

cat("\n")

# ---- Duplicated sample ID checks ----

cat("=== DUPLICATED SAMPLE ID CHECKS ===\n\n")

dup_221561 <- colnames(counts_221561)[duplicated(colnames(counts_221561))]
dup_197677 <- colnames(counts_197677)[duplicated(colnames(counts_197677))]
dup_160269 <- colnames(counts_160269)[duplicated(colnames(counts_160269))]

if (length(dup_221561) > 0) {
  errors <- c(errors, paste("GSE221561 duplicated IDs:", paste(dup_221561, collapse = ", ")))
  cat("GSE221561 DUPLICATED:", paste(dup_221561, collapse = ", "), "\n")
} else {
  cat("GSE221561: 0 duplicated sample IDs\n")
}

if (length(dup_197677) > 0) {
  errors <- c(errors, paste("GSE197677 duplicated IDs:", paste(dup_197677, collapse = ", ")))
  cat("GSE197677 DUPLICATED:", paste(dup_197677, collapse = ", "), "\n")
} else {
  cat("GSE197677: 0 duplicated sample IDs\n")
}

if (length(dup_160269) > 0) {
  errors <- c(errors, paste("GSE160269 duplicated IDs:", paste(dup_160269, collapse = ", ")))
  cat("GSE160269 DUPLICATED:", paste(dup_160269, collapse = ", "), "\n")
} else {
  cat("GSE160269: 0 duplicated sample IDs\n")
}

cat("\n")

# ---- Count integrity checks ----

cat("=== COUNT INTEGRITY CHECKS ===\n\n")

check_counts <- function(mat, name) {
  nas <- sum(is.na(mat))
  negs <- sum(mat < 0, na.rm = TRUE)
  is_int <- all(mat == floor(mat), na.rm = TRUE)
  cat(sprintf("%s: NA=%d, negative=%d, integer_like=%s\n", name, nas, negs, is_int))
  list(nas = nas, negs = negs, is_int = is_int)
}

c1 <- check_counts(counts_221561, "GSE221561")
c2 <- check_counts(counts_197677, "GSE197677")
c3 <- check_counts(counts_160269, "GSE160269")

if (c1$nas > 0) errors <- c(errors, "GSE221561: has NA counts")
if (c1$negs > 0) errors <- c(errors, "GSE221561: has negative counts")
if (!c1$is_int) errors <- c(errors, "GSE221561: counts not integer-like")

if (c2$nas > 0) errors <- c(errors, "GSE197677: has NA counts")
if (c2$negs > 0) errors <- c(errors, "GSE197677: has negative counts")
if (!c2$is_int) errors <- c(errors, "GSE197677: counts not integer-like")

if (c3$nas > 0) errors <- c(errors, "GSE160269: has NA counts")
if (c3$negs > 0) errors <- c(errors, "GSE160269: has negative counts")
if (!c3$is_int) errors <- c(errors, "GSE160269: counts not integer-like")

cat("\n")

# ---- Column name to metadata match ----

cat("=== COLUMN-METADATA MATCH CHECKS ===\n\n")

check_match <- function(count_mat, meta, name) {
  count_ids <- colnames(count_mat)
  meta_ids <- meta$sample
  match_count <- sum(count_ids %in% meta_ids)
  cat(sprintf("%s: count_cols=%d, meta_rows=%d, matched=%d\n", name, length(count_ids), length(meta_ids), match_count))
  unmatched <- setdiff(count_ids, meta_ids)
  if (length(unmatched) > 0) cat(sprintf("  Unmatched in metadata: %s\n", paste(unmatched, collapse = ", ")))
  unmatched_meta <- setdiff(meta_ids, count_ids)
  if (length(unmatched_meta) > 0) cat(sprintf("  Unmatched in counts: %s\n", paste(unmatched_meta, collapse = ", ")))
  list(match_count = match_count, unmatched = unmatched, unmatched_meta = unmatched_meta)
}

m1 <- check_match(counts_221561, meta_221561, "GSE221561")
m2 <- check_match(counts_197677, meta_197677, "GSE197677")
m3 <- check_match(counts_160269, meta_160269, "GSE160269")

if (length(m1$unmatched) > 0) errors <- c(errors, paste("GSE221561 unmatched:", paste(m1$unmatched, collapse = ", ")))
if (length(m2$unmatched) > 0) errors <- c(errors, paste("GSE197677 unmatched:", paste(m2$unmatched, collapse = ", ")))
if (length(m3$unmatched) > 0) errors <- c(errors, paste("GSE160269 unmatched:", paste(m3$unmatched, collapse = ", ")))

cat("\n")

# ---- Common genes check ----

cat("=== COMMON GENES CHECK ===\n\n")

common_in_221561 <- sum(common_genes %in% rownames(counts_221561))
common_in_197677 <- sum(common_genes %in% rownames(counts_197677))
common_in_160269 <- sum(common_genes %in% rownames(counts_160269))

cat(sprintf("Three-dataset common genes: %d\n", length(common_genes)))
cat(sprintf("  Found in GSE221561: %d\n", common_in_221561))
cat(sprintf("  Found in GSE197677: %d\n", common_in_197677))
cat(sprintf("  Found in GSE160269: %d\n", common_in_160269))

if (common_in_221561 != length(common_genes)) errors <- c(errors, "GSE221561 missing some common genes")
if (common_in_197677 != length(common_genes)) errors <- c(errors, "GSE197677 missing some common genes")
if (common_in_160269 != length(common_genes)) errors <- c(errors, "GSE160269 missing some common genes")

cat("\n")

# ---- Define primary analysis sample sets ----

cat("=== PRIMARY ANALYSIS SAMPLE SETS ===\n\n")

primary_221561 <- colnames(counts_221561)
cat(sprintf("GSE221561 primary: %d samples\n", length(primary_221561)))

primary_197677 <- colnames(counts_197677)
cat(sprintf("GSE197677 primary: %d samples\n", length(primary_197677)))

primary_160269 <- meta_160269$sample[meta_160269$epithelial_cells >= 30]
excluded_160269 <- meta_160269$sample[meta_160269$epithelial_cells < 30]

cat(sprintf("GSE160269 primary: %d samples (epi_cells >= 30)\n", length(primary_160269)))
cat(sprintf("GSE160269 excluded: %d samples (kept in inventory, not in primary analysis)\n", length(excluded_160269)))
if (length(excluded_160269) > 0) {
  cat("  Excluded samples:", paste(excluded_160269, collapse = ", "), "\n")
  # Mark PRIMARY_ANALYSIS = FALSE in qc_table
  qc_table$PRIMARY_ANALYSIS <- ifelse(qc_table$sample %in% excluded_160269, FALSE, TRUE)
  qc_table$exclusion_reason <- ifelse(qc_table$sample %in% excluded_160269, "LOW_EPITHELIAL_CELL_COUNT", "")
}

cat("\n")

# ---- Summary ----

if (length(errors) > 0) {
  cat("========================================\n")
  cat("VERIFICATION FAILED\n")
  cat("========================================\n\n")
  for (e in errors) cat("  ERROR:", e, "\n")
  cat("\nSTOP.\n")
  saveRDS(list(status = "FAIL", errors = errors), file.path(base_dir, "06_tables/cross_dataset/epithelial_pseudobulk/QC/STEP19B0_input_verification.rds"))
  stop("Verification failed. See errors above.")
} else {
  cat("========================================\n")
  cat("VERIFICATION PASSED\n")
  cat("========================================\n")
  cat("\nAll dimension, integrity, and match checks passed.\n")
  
  # Save verification results
  verification <- list(
    status = "PASS",
    GSE221561_cols = ncol(counts_221561),
    GSE197677_cols = ncol(counts_197677),
    GSE160269_cols = ncol(counts_160269),
    GSE160269_primary = length(primary_160269),
    GSE160269_excluded = excluded_160269,
    common_genes = length(common_genes),
    primary_221561 = primary_221561,
    primary_197677 = primary_197677,
    primary_160269 = primary_160269
  )
  saveRDS(verification, file.path(base_dir, "06_tables/cross_dataset/epithelial_pseudobulk/QC/STEP19B0_input_verification.rds"))
  
  # Save updated QC table with exclusion flags
  write.csv(qc_table, file.path(base_dir, "06_tables/cross_dataset/epithelial_pseudobulk/STEP19A_pseudobulk_sample_QC.csv"), row.names = FALSE)
}

cat("\nDone.\n")
