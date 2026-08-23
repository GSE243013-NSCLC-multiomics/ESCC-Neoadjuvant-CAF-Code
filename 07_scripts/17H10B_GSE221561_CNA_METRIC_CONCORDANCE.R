# ============================================================
# ESCC Neoadjuvant Project
# STEP 17H10B
#
# METRIC CONCORDANCE / DISCORDANCE AUDIT
#
# Pure audit. NO InferCNA, NO findMalignant, NO cell filtering.
# ============================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
})

cat("\n============================================================\n")
cat("STEP 17H10B: METRIC CONCORDANCE / DISCORDANCE AUDIT\n")
cat("============================================================\n\n")

# ============================================================
# 1. LOAD LANDSCAPE
# ============================================================

landscape_file <- "06_tables/cross_dataset/malignant_epithelial/GSE221561_CNV_strategy_audit/STEP17H10_sample_level_CNA_landscape.csv"

if (!file.exists(landscape_file)) stop("Missing: ", landscape_file)

dt <- fread(landscape_file)
cat("Samples: ", nrow(dt), "\n", sep = "")

# ============================================================
# 2. RANK MATRIX
# ============================================================

cat("\n===== 2. METRIC RANK MATRIX =====\n\n")

metrics <- c("CNA_signal_ratio", "tumor_median_cnaCor", "correlation_separation",
             "tumor_above_normal95", "coherent_CNA_arms")
metric_labels <- c("signal_ratio", "cnaCor", "corr_separation", "above95", "coherent_arms")

rank_matrix <- data.table(sample = dt$sample, treatment = dt$treatment)

for (j in seq_along(metrics)) {
  rank_matrix[, (metric_labels[j]) := rank(-dt[[metrics[j]]])]
}

cat("Rank matrix:\n")
print(rank_matrix, nrows = 20)

fwrite(rank_matrix, file.path("06_tables/cross_dataset/malignant_epithelial/GSE221561_CNV_strategy_audit",
                               "STEP17H10_CNA_metric_rank_matrix.csv"), bom = TRUE)
cat("Saved: STEP17H10_CNA_metric_rank_matrix.csv\n")

# ============================================================
# 3. PER-SAMPLE RANK STATISTICS
# ============================================================

cat("\n===== 3. PER-SAMPLE RANK STATISTICS =====\n\n")

rank_stats <- data.table(
  sample = dt$sample,
  treatment = dt$treatment,
  rank_median = apply(rank_matrix[, ..metric_labels], 1, median, na.rm = TRUE),
  rank_min = apply(rank_matrix[, ..metric_labels], 1, min, na.rm = TRUE),
  rank_max = apply(rank_matrix[, ..metric_labels], 1, max, na.rm = TRUE),
  rank_iqr = apply(rank_matrix[, ..metric_labels], 1, IQR, na.rm = TRUE)
)

cat("Rank statistics:\n")
print(rank_stats, nrows = 20)

# ============================================================
# 4. CONCORDANCE STATUS
# ============================================================

cat("\n===== 4. CONCORDANCE CLASSIFICATION =====\n\n")

# Classification based on rank consistency:
# HIGH_CONCORDANCE: IQR <= 1.5 AND range <= 3
# INTERMEDIATE_CONCORDANCE: IQR <= 2.5 AND range <= 5
# DISCORDANT_REVIEW: everything else

rank_stats[, concordance_status := fcase(
  rank_iqr <= 1.5 & (rank_max - rank_min) <= 3, "HIGH_CONCORDANCE",
  rank_iqr <= 2.5 & (rank_max - rank_min) <= 5, "INTERMEDIATE_CONCORDANCE",
  default = "DISCORDANT_REVIEW"
)]

cat("Concordance classification:\n\n")
for (i in seq_len(nrow(rank_stats))) {
  cat(sprintf("  %-10s  median_rank=%.1f  IQR=%.1f  range=%d-%d  %s\n",
              rank_stats$sample[i], rank_stats$rank_median[i],
              rank_stats$rank_iqr[i], rank_stats$rank_min[i], rank_stats$rank_max[i],
              rank_stats$concordance_status[i]))
}

# ============================================================
# 5. SPECIAL REVIEW: S0520T AND S2487T
# ============================================================

cat("\n===== 5. SPECIAL REVIEW =====\n\n")

cat("S0520T:\n")
s0520 <- rank_stats[sample == "S0520T"]
s0520_ranks <- as.numeric(rank_matrix[sample == "S0520T", ..metric_labels])
cat("  Ranks: ", paste(s0520_ranks, collapse = ", "), "\n", sep = "")
cat("  Concordance: ", s0520$concordance_status, "\n", sep = "")

cat("\nS2487T:\n")
s2487 <- rank_stats[sample == "S2487T"]
s2487_ranks <- as.numeric(rank_matrix[sample == "S2487T", ..metric_labels])
cat("  Ranks: ", paste(s2487_ranks, collapse = ", "), "\n", sep = "")
cat("  Concordance: ", s2487$concordance_status, "\n", sep = "")

# ============================================================
# 6. SAVE
# ============================================================

# Merge concordance into landscape
dt_merged <- merge(dt, rank_stats[, .(sample, rank_median, rank_min, rank_max, rank_iqr, concordance_status)],
                   by = "sample", all.x = TRUE)

fwrite(dt_merged, file.path("06_tables/cross_dataset/malignant_epithelial/GSE221561_CNV_strategy_audit",
                             "STEP17H10_sample_level_CNA_landscape.csv"), bom = TRUE)

cat("\nSTEP 17H10B COMPLETE\n")
