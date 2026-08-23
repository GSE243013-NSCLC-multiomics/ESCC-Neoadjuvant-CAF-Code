# ============================================================
# ESCC Neoadjuvant Project
# STEP 17H10D
#
# SAMPLE-SPECIFIC REVIEW TABLE
#
# Pure audit. NO InferCNA, NO findMalignant, NO cell filtering.
# ============================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
})

cat("\n============================================================\n")
cat("STEP 17H10D: SAMPLE-SPECIFIC REVIEW TABLE\n")
cat("============================================================\n\n")

# ============================================================
# 1. LOAD LANDSCAPE + SEPARABILITY + CONCORDANCE
# ============================================================

audit_dir <- "06_tables/cross_dataset/malignant_epithelial/GSE221561_CNV_strategy_audit"

landscape <- fread(file.path(audit_dir, "STEP17H10_sample_level_CNA_landscape.csv"))
separability <- fread(file.path(audit_dir, "STEP17H10_cell_level_CNA_separability.csv"))
rank_matrix <- fread(file.path(audit_dir, "STEP17H10_CNA_metric_rank_matrix.csv"))

# ============================================================
# 2. BUILD REVIEW TABLE
# ============================================================

cat("===== 2. BUILD REVIEW TABLE =====\n\n")

review <- merge(landscape, separability[, .(sample, dual_metric_candidate_fraction)], by = "sample", all.x = TRUE)
rank_cols <- rank_matrix[, .(sample, rank_signal = signal_ratio, rank_cnaCor = cnaCor,
                              rank_separation = corr_separation, rank_above95 = above95,
                              rank_arms = coherent_arms)]
review <- merge(review, rank_cols, by = "sample", all.x = TRUE)

# Concordance status from landscape
# (already computed in 17H10B)

# ============================================================
# 3. REVIEW PRIORITY
# ============================================================

cat("\n===== 3. ASSIGN REVIEW PRIORITY =====\n\n")

review[, review_priority := "STANDARD"]

# REVIEW_EXTREME: if rank_median <= 2 AND rank_range <= 3
review[, rank_range := rank_max - rank_min]
review[rank_median <= 2 & rank_range <= 3, review_priority := "REVIEW_EXTREME"]

# REVIEW_DISCORDANT: if rank_iqr >= 3 OR rank_range >= 6
review[rank_iqr >= 3 | rank_range >= 6, review_priority := "REVIEW_DISCORDANT"]

# Special cases based on data (not hardcoded names)
# High signal ratio but low cnaCor/arms
review[CNA_signal_ratio > 8 & tumor_median_cnaCor < 0.55 & coherent_CNA_arms < 6,
       review_priority := "REVIEW_DISCORDANT"]

# Very high across all metrics
review[CNA_signal_ratio > 15 & tumor_median_cnaCor > 0.7 & coherent_CNA_arms > 20,
       review_priority := "REVIEW_EXTREME"]

cat("Review priorities:\n\n")
for (i in seq_len(nrow(review))) {
  cat(sprintf("  %-10s  %s\n", review$sample[i], review$review_priority[i]))
}

# ============================================================
# 4. SAVE
# ============================================================

# Use available columns
review_output <- review[, .(sample, treatment, tumor_cells,
                             CNA_signal_ratio, rank_signal,
                             tumor_median_cnaCor,
                             correlation_separation,
                             tumor_above_normal95,
                             coherent_CNA_arms,
                             dual_metric_candidate_fraction,
                             concordance_status,
                             review_priority)]

fwrite(review_output, file.path(audit_dir, "STEP17H10_sample_review_table.csv"), bom = TRUE)
cat("\nSaved: STEP17H10_sample_review_table.csv\n")

cat("\nSTEP 17H10D COMPLETE\n")
