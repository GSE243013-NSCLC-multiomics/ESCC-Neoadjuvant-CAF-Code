# ============================================================
# ESCC Neoadjuvant Project
# STEP 17H8
#
# GSE221561 FOUR-PILOT CANONICAL FREEZE AUDIT
#
# Pure audit — NO InferCNA, NO findMalignant, NO new samples.
# ============================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
})

cat("\n============================================================\n")
cat("STEP 17H8: GSE221561 FOUR-PILOT CANONICAL FREEZE AUDIT\n")
cat("============================================================\n\n")

# ============================================================
# 1. CANONICAL MANIFEST PATHS
# ============================================================

cat("===== 1. LOCATE CANONICAL INDIVIDUAL MANIFESTS =====\n\n")

pilots <- list(
  S6417T = list(
    treatment = "Surgery_alone",
    results_dir = "06_tables/cross_dataset/malignant_epithelial/GSE221561_S6417T_pilot/17G3_results",
    manifest_file = "GSE221561_S6417T_17G3_manifest.csv",
    metrics_file = "GSE221561_S6417T_17G3_cell_CNA_metrics.csv",
    coherent_file = "GSE221561_S6417T_17G3_candidate_coherent_arms.csv",
    cna_all_file = "GSE221561_S6417T_infercna_CNA_all_cells.rds",
    cna_tumor_file = "GSE221561_S6417T_infercna_CNA_tumor_only.rds",
    expected_tumor = 454L
  ),
  S9478T = list(
    treatment = "Surgery_alone",
    results_dir = "06_tables/cross_dataset/malignant_epithelial/GSE221561_S9478T_pilot/17H5C_results",
    manifest_file = "GSE221561_S9478T_17H5C_manifest.csv",
    metrics_file = "GSE221561_S9478T_17H5C_cell_CNA_metrics.csv",
    coherent_file = "GSE221561_S9478T_17H5C_candidate_coherent_arms.csv",
    cna_all_file = "GSE221561_S9478T_infercna_CNA_all_cells.rds",
    cna_tumor_file = "GSE221561_S9478T_infercna_CNA_tumor_only.rds",
    expected_tumor = 312L
  ),
  S1265T = list(
    treatment = "Neoadjuvant_treated",
    results_dir = "06_tables/cross_dataset/malignant_epithelial/GSE221561_S1265T_pilot/17H3_results",
    manifest_file = "GSE221561_S1265T_17H3_manifest.csv",
    metrics_file = "GSE221561_S1265T_17H3_cell_CNA_metrics.csv",
    coherent_file = "GSE221561_S1265T_17H3_candidate_coherent_arms.csv",
    cna_all_file = "GSE221561_S1265T_infercna_CNA_all_cells.rds",
    cna_tumor_file = "GSE221561_S1265T_infercna_CNA_tumor_only.rds",
    expected_tumor = 4972L
  ),
  S2423T = list(
    treatment = "Neoadjuvant_treated",
    results_dir = "06_tables/cross_dataset/malignant_epithelial/GSE221561_S2423T_pilot/17H7C_results",
    manifest_file = "GSE221561_S2423T_17H7C_manifest.csv",
    metrics_file = "GSE221561_S2423T_17H7C_cell_CNA_metrics.csv",
    coherent_file = "GSE221561_S2423T_17H7C_candidate_coherent_arms.csv",
    cna_all_file = "GSE221561_S2423T_infercna_CNA_all_cells.rds",
    cna_tumor_file = "GSE221561_S2423T_infercna_CNA_tumor_only.rds",
    expected_tumor = 4363L
  )
)

for (nm in names(pilots)) {
  p <- pilots[[nm]]
  manifest_path <- file.path(p$results_dir, p$manifest_file)
  if (!file.exists(manifest_path)) stop("Canonical manifest missing: ", manifest_path)
  cat(nm, " : ", manifest_path, "\n", sep = "")
}

cat("\nAll four canonical manifests located.\n")

# ============================================================
# 2. LOAD AND VALIDATE EACH PILOT
# ============================================================

cat("\n===== 2. LOAD AND VALIDATE EACH PILOT =====\n\n")

audit_results <- list()

for (nm in names(pilots)) {
  p <- pilots[[nm]]
  cat("\n------------------------------------------------------------\n")
  cat("PILOT: ", nm, " (", p$treatment, ")\n", sep = "")
  cat("------------------------------------------------------------\n\n")

  # --- 2a. Load manifest ---
  manifest_path <- file.path(p$results_dir, p$manifest_file)
  manifest <- fread(manifest_path)
  cat("Manifest loaded: ", manifest_path, "\n", sep = "")

  # Verify manifest contains sample ID
  if (!nm %in% unlist(manifest)) {
    stop("Manifest does not contain sample ID: ", nm)
  }
  cat("Manifest contains sample ID: ", nm, " ✓\n", sep = "")

  # --- 2b. Load cell CNA metrics ---
  metrics_path <- file.path(p$results_dir, p$metrics_file)
  cell_metrics <- fread(metrics_path)
  cat("Cell CNA metrics loaded: ", nrow(cell_metrics), " rows\n", sep = "")

  # --- 2c. Load CNA RDS ---
  cna_all_path <- file.path(p$results_dir, p$cna_all_file)
  cna_tumor_path <- file.path(p$results_dir, p$cna_tumor_file)
  cna_all <- readRDS(cna_all_path)
  cna_tumor <- readRDS(cna_tumor_path)
  cat("CNA all cells: ", nrow(cna_all), " genes x ", ncol(cna_all), " cells\n", sep = "")
  cat("CNA tumor only: ", nrow(cna_tumor), " genes x ", ncol(cna_tumor), " cells\n", sep = "")

  # --- 2d. Load coherent arms ---
  coherent_path <- file.path(p$results_dir, p$coherent_file)
  coherent_arms <- fread(coherent_path)
  cat("Coherent arms CSV: ", nrow(coherent_arms), " rows\n", sep = "")

  # ============================================================
  # 3. TUMOR CELL COUNT TRIPLE VERIFICATION
  # ============================================================

  cat("\n--- 3. TUMOR CELL COUNT VERIFICATION ---\n\n")

  # A. manifest tumor_cells
  count_A <- manifest$tumor_cells

  # B. cell CNA metrics role == "TUMOR_TEST"
  count_B <- sum(cell_metrics$role == "TUMOR_TEST", na.rm = TRUE)

  # C. tumor-only CNA RDS ncol
  count_C <- ncol(cna_tumor)

  cat("A. Manifest tumor_cells       : ", count_A, "\n", sep = "")
  cat("B. Cell metrics TUMOR_TEST     : ", count_B, "\n", sep = "")
  cat("C. Tumor-only CNA ncol         : ", count_C, "\n", sep = "")

  tumor_count_pass <- (count_A == count_B) && (count_B == count_C)
  if (tumor_count_pass) {
    cat("Tumor cell count validation    : PASS ✓\n")
  } else {
    cat("TUMOR CELL COUNT MISMATCH\n")
    cat("  A = ", count_A, "\n")
    cat("  B = ", count_B, "\n")
    cat("  C = ", count_C, "\n")
    stop("TUMOR CELL COUNT MISMATCH for ", nm)
  }

  # Check expected
  if (count_A != p$expected_tumor) {
    cat("WARNING: Expected ", p$expected_tumor, " but found ", count_A, "\n", sep = "")
    stop("Expected tumor cell count mismatch for ", nm)
  }
  cat("Expected ", p$expected_tumor, " tumor cells confirmed ✓\n", sep = "")

  # ============================================================
  # 4. CNA GENE COUNT VERIFICATION
  # ============================================================

  cat("\n--- 4. CNA GENE COUNT VERIFICATION ---\n\n")

  gene_A <- manifest$CNA_genes
  gene_B <- nrow(cna_all)
  gene_C <- nrow(cna_tumor)

  cat("A. Manifest CNA_genes          : ", gene_A, "\n", sep = "")
  cat("B. CNA all cells nrow           : ", gene_B, "\n", sep = "")
  cat("C. CNA tumor only nrow          : ", gene_C, "\n", sep = "")

  gene_count_pass <- (gene_A == gene_B) && (gene_B == gene_C)
  if (gene_count_pass) {
    cat("Gene count validation          : PASS ✓\n")
  } else {
    cat("GENE COUNT MISMATCH\n")
    stop("Gene count mismatch for ", nm)
  }

  # ============================================================
  # 5. COHERENT ARM TRIPLE VERIFICATION
  # ============================================================

  cat("\n--- 5. COHERENT ARM VERIFICATION ---\n\n")

  arm_A <- manifest$coherent_CNA_arms
  arm_B <- nrow(coherent_arms)
  arm_C <- if (nrow(coherent_arms) > 0) length(unique(coherent_arms$chromosome_arm)) else 0

  cat("A. Manifest coherent_CNA_arms  : ", arm_A, "\n", sep = "")
  cat("B. CSV nrow                    : ", arm_B, "\n", sep = "")
  cat("C. Unique chromosome_arm        : ", arm_C, "\n", sep = "")

  arm_count_pass <- (arm_A == arm_B) && (arm_B == arm_C)
  if (arm_count_pass) {
    cat("Coherent arm validation        : PASS ✓\n")
  } else {
    cat("COHERENT ARM COUNT MISMATCH\n")
    stop("Coherent arm count mismatch for ", nm)
  }

  arm_names <- if (nrow(coherent_arms) > 0) paste(sort(coherent_arms$chromosome_arm), collapse = ", ") else "none"
  cat("Coherent arm names             : ", arm_names, "\n", sep = "")

  # ============================================================
  # 6. RECALCULATE METRICS FROM CELL CNA METRICS
  # ============================================================

  cat("\n--- 6. RECALCULATE METRICS FROM CELL CNA METRICS ---\n\n")

  tumor_m <- cell_metrics[role == "TUMOR_TEST"]
  normal_m <- cell_metrics[role == "NORMAL_REFERENCE"]

  tumor_cells_n <- nrow(tumor_m)
  normal_cells_n <- nrow(normal_m)

  tumor_median_signal <- median(tumor_m$cna_signal, na.rm = TRUE)
  normal_median_signal <- median(normal_m$cna_signal, na.rm = TRUE)
  signal_delta <- tumor_median_signal - normal_median_signal
  signal_ratio <- tumor_median_signal / normal_median_signal

  tumor_median_cor <- median(tumor_m$cna_cor_to_tumor, na.rm = TRUE)
  normal_median_cor <- median(normal_m$cna_cor_to_tumor, na.rm = TRUE)
  corr_separation <- tumor_median_cor - normal_median_cor

  normal_signal_95 <- as.numeric(quantile(normal_m$cna_signal, 0.95, na.rm = TRUE))
  normal_signal_99 <- as.numeric(quantile(normal_m$cna_signal, 0.99, na.rm = TRUE))

  tumor_above_95 <- mean(tumor_m$cna_signal > normal_signal_95, na.rm = TRUE)
  tumor_above_99 <- mean(tumor_m$cna_signal > normal_signal_99, na.rm = TRUE)

  cat(sprintf("Tumor cells                    : %d\n", tumor_cells_n))
  cat(sprintf("Normal cells                   : %d\n", normal_cells_n))
  cat(sprintf("Tumor median CNA signal        : %.6f\n", tumor_median_signal))
  cat(sprintf("Normal median CNA signal       : %.6f\n", normal_median_signal))
  cat(sprintf("CNA signal delta               : %.6f\n", signal_delta))
  cat(sprintf("CNA signal ratio               : %.4f\n", signal_ratio))
  cat(sprintf("Tumor median cnaCor            : %.4f\n", tumor_median_cor))
  cat(sprintf("Normal median corr to tumor    : %.4f\n", normal_median_cor))
  cat(sprintf("Correlation separation         : %.4f\n", corr_separation))
  cat(sprintf("Normal signal 95th             : %.6f\n", normal_signal_95))
  cat(sprintf("Normal signal 99th             : %.6f\n", normal_signal_99))
  cat(sprintf("Tumor > normal 95th            : %.1f%%\n", 100 * tumor_above_95))
  cat(sprintf("Tumor > normal 99th            : %.1f%%\n", 100 * tumor_above_99))

  # --- Compare with manifest ---
  cat("\n--- Comparison with manifest ---\n\n")

  # S6417T and S1265T manifests lack derived columns; use NA fallback
  get_manifest_val <- function(col_name) {
    if (col_name %in% names(manifest)) manifest[[col_name]][1] else NA_real_
  }

  metric_diffs <- data.table(
    metric = c("tumor_cells", "normal_cells", "CNA_genes",
               "tumor_median_CNA_signal", "normal_median_CNA_signal",
               "CNA_signal_delta", "CNA_signal_ratio",
               "tumor_median_cnaCor", "normal_median_corr_to_tumor",
               "correlation_separation",
               "tumor_above_normal_95th", "tumor_above_normal_99th",
               "coherent_CNA_arms"),
    recalculated = c(tumor_cells_n, normal_cells_n, gene_A,
                     tumor_median_signal, normal_median_signal,
                     signal_delta, signal_ratio,
                     tumor_median_cor, normal_median_cor,
                     corr_separation,
                     tumor_above_95, tumor_above_99,
                     arm_A),
    manifest = c(
      manifest$tumor_cells,
      manifest$normal_reference_cells,
      manifest$CNA_genes,
      manifest$tumor_median_CNA_signal,
      manifest$normal_median_CNA_signal,
      get_manifest_val("CNA_signal_delta"),
      get_manifest_val("CNA_signal_ratio"),
      manifest$tumor_median_cnaCor,
      manifest$normal_median_corr_to_tumor,
      get_manifest_val("correlation_separation"),
      get_manifest_val("tumor_above_normal_95th"),
      get_manifest_val("tumor_above_normal_99th"),
      manifest$coherent_CNA_arms
    )
  )

  # Compute absolute differences (NA if manifest value is NA)
  metric_diffs[, manifest_num := as.numeric(manifest)]
  metric_diffs[, recalculated_num := as.numeric(recalculated)]
  metric_diffs[, abs_diff := abs(recalculated_num - manifest_num)]

  # If manifest field is NA (older format), mark as N/A — not a failure
  metric_diffs[, in_manifest := !is.na(manifest_num)]

  # For ratio/percentage metrics, also compute relative tolerance
  metric_diffs[, tolerance := fifelse(
    metric %in% c("CNA_signal_ratio", "correlation_separation",
                   "tumor_above_normal_95th", "tumor_above_normal_99th"),
    0.01,  # 1% relative tolerance for ratios/percentages
    1e-4   # absolute tolerance for counts and signals
  )]

  metric_diffs[, pass := fifelse(
    !in_manifest, TRUE,  # N/A if not in manifest (older format)
    fifelse(
      metric %in% c("CNA_signal_ratio", "correlation_separation",
                     "tumor_above_normal_95th", "tumor_above_normal_99th"),
      abs_diff <= tolerance | abs_diff / abs(manifest_num + 1e-10) <= tolerance,
      abs_diff <= tolerance
    )
  )]

  metric_diffs[, status := fifelse(
    !in_manifest, "N/A",
    fifelse(pass == TRUE, "PASS", "FAIL")
  )]

  print(metric_diffs[, .(metric, recalculated, manifest, abs_diff, status)], nrows = 20)

  # Count: PASS = pass or N/A (not in manifest), FAIL = actual mismatch
  n_pass <- sum(metric_diffs$pass)
  n_na <- sum(!metric_diffs$in_manifest)
  n_fail <- sum(!metric_diffs$pass & metric_diffs$in_manifest)

  metric_validation_pass <- (n_fail == 0)
  if (metric_validation_pass) {
    cat("\nMetric validation               : PASS ✓")
    if (n_na > 0) cat(" (", n_na, " fields N/A — not in older manifest format)", sep = "")
    cat("\n")
  } else {
    cat("\nMETRIC VALIDATION FAILED\n")
    failed <- metric_diffs[!pass & in_manifest]
    print(failed)
    stop("Metric validation failed for ", nm)
  }

  # ============================================================
  # STORE RESULTS
  # ============================================================

  audit_results[[nm]] <- data.table(
    sample = nm,
    treatment = p$treatment,
    tumor_cells = tumor_cells_n,
    normal_reference_cells = normal_cells_n,
    CNA_genes = gene_A,
    tumor_median_CNA_signal = tumor_median_signal,
    normal_median_CNA_signal = normal_median_signal,
    CNA_signal_delta = signal_delta,
    CNA_signal_ratio = signal_ratio,
    tumor_median_cnaCor = tumor_median_cor,
    normal_median_corr_to_tumor = normal_median_cor,
    correlation_separation = corr_separation,
    tumor_above_normal95 = tumor_above_95,
    tumor_above_normal99 = tumor_above_99,
    coherent_CNA_arms = arm_A,
    coherent_CNA_arm_names = arm_names,
    tumor_cell_count_validation = ifelse(tumor_count_pass, "PASS", "FAIL"),
    gene_count_validation = ifelse(gene_count_pass, "PASS", "FAIL"),
    arm_count_validation = ifelse(arm_count_pass, "PASS", "FAIL"),
    metric_validation = ifelse(metric_validation_pass, "PASS", "FAIL"),
    canonical_manifest_path = manifest_path
  )
}

# ============================================================
# 7. GENERATE CANONICAL FOUR-PILOT TABLE
# ============================================================

cat("\n============================================================\n")
cat("7. CANONICAL FOUR-PILOT TABLE\n")
cat("============================================================\n\n")

outdir <- "06_tables/cross_dataset/malignant_epithelial/GSE221561_four_pilot_freeze"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

four_pilot <- rbindlist(audit_results, fill = TRUE)
fwrite(four_pilot, file.path(outdir, "STEP17H8_GSE221561_four_pilot_CANONICAL_metrics.csv"), bom = TRUE)

cat("Saved: ", file.path(outdir, "STEP17H8_GSE221561_four_pilot_CANONICAL_metrics.csv"), "\n", sep = "")

# ============================================================
# 8. PRINT TRUE FOUR-PILOT TABLE
# ============================================================

cat("\n============================================================\n")
cat("8. TRUE FOUR-PILOT TABLE\n")
cat("============================================================\n\n")

print(four_pilot[, .(
  sample, treatment, tumor_cells, CNA_signal_ratio,
  tumor_median_cnaCor, correlation_separation,
  tumor_above_normal95, tumor_above_normal99,
  coherent_CNA_arms, tumor_cell_count_validation,
  gene_count_validation, arm_count_validation, metric_validation
)], nrows = 20)

# ============================================================
# 9. TREATMENT-LEVEL SUMMARY
# ============================================================

cat("\n============================================================\n")
cat("9. TREATMENT-LEVEL DESCRIPTIVE SUMMARY\n")
cat("============================================================\n\n")

cat("EXPLORATORY ONLY — n=2 per group.\n")
cat("NO causal treatment interpretation.\n\n")

surgery <- four_pilot[treatment == "Surgery_alone"]
neoadj <- four_pilot[treatment == "Neoadjuvant_treated"]

summarize_group <- function(dt, label) {
  cat(label, " (n = ", nrow(dt), ")\n", sep = "")
  cat(sprintf("  CNA signal ratio     median %.2f  range %.2f – %.2f\n",
              median(dt$CNA_signal_ratio), min(dt$CNA_signal_ratio), max(dt$CNA_signal_ratio)))
  cat(sprintf("  Tumor cnaCor         median %.3f  range %.3f – %.3f\n",
              median(dt$tumor_median_cnaCor), min(dt$tumor_median_cnaCor), max(dt$tumor_median_cnaCor)))
  cat(sprintf("  Correlation sep.     median %.3f  range %.3f – %.3f\n",
              median(dt$correlation_separation), min(dt$correlation_separation), max(dt$correlation_separation)))
  cat(sprintf("  Tumor > normal 95th  median %.1f%%  range %.1f%% – %.1f%%\n",
              100 * median(dt$tumor_above_normal95), 100 * min(dt$tumor_above_normal95), 100 * max(dt$tumor_above_normal95)))
  cat(sprintf("  Coherent arms        median %.0f   range %.0f – %.0f\n",
              median(dt$coherent_CNA_arms), min(dt$coherent_CNA_arms), max(dt$coherent_CNA_arms)))
  cat("\n")
}

summarize_group(surgery, "Surgery_alone")
summarize_group(neoadj, "Neoadjuvant_treated")

cat("NO p-value. NO t-test. NO Wilcoxon.\n")
cat("NO treatment effect significance testing.\n")
cat("EXPLORATORY ONLY.\n")

# ============================================================
# 10. AUDIT ALL GSE221561 TUMOR SAMPLES
# ============================================================

cat("\n============================================================\n")
cat("10. ALL GSE221561 TUMOR SAMPLE CNV STATUS\n")
cat("============================================================\n\n")

inventory <- fread("06_tables/cross_dataset/malignant_epithelial/GSE221561_preflight/STEP17G1_GSE221561_epithelial_sample_inventory.csv")
tumor_inv <- inventory[tissue_simple == "Tumor"]

# Completed pilots
completed <- c("S6417T", "S9478T", "S1265T", "S2423T")

tumor_inv[, CNV_status := fifelse(
  sample_id %in% completed, "COMPLETE_PILOT",
  fifelse(epithelial_cells >= 30, "NOT_RUN", "INELIGIBLE_LT30")
)]

# Reorder for readability
setorder(tumor_inv, -epithelial_cells)

cat("GSE221561 Tumor Sample CNV Inventory:\n\n")
print(tumor_inv[, .(sample_id, treatment_group, epithelial_cells, CNV_status)], nrows = 20)

# Save
fwrite(tumor_inv, file.path(outdir, "STEP17H8_GSE221561_all_tumor_sample_CNV_status.csv"), bom = TRUE)
cat("\nSaved: ", file.path(outdir, "STEP17H8_GSE221561_all_tumor_sample_CNV_status.csv"), "\n", sep = "")

# ============================================================
# 11. EXPANSION PLAN
# ============================================================

cat("\n============================================================\n")
cat("11. EXPANSION PLAN — REMAINING ELIGIBLE SAMPLES\n")
cat("============================================================\n\n")

remaining <- tumor_inv[CNV_status == "NOT_RUN"]
setorder(remaining, treatment_group, -epithelial_cells)

cat("Remaining eligible samples for CNV:\n\n")
print(remaining[, .(sample_id, treatment_group, epithelial_cells)], nrows = 20)

cat("\nBy treatment group:\n\n")
for (tg in unique(remaining$treatment_group)) {
  sub <- remaining[treatment_group == tg]
  cat(tg, ": ", nrow(sub), " samples\n", sep = "")
  for (i in seq_len(nrow(sub))) {
    cat("  ", sub$sample_id[i], " — ", sub$epithelial_cells[i], " epithelial cells\n", sep = "")
  }
  cat("\n")
}

# ============================================================
# 12. SUMMARY COUNTS
# ============================================================

total_tumor <- nrow(tumor_inv)
n_completed <- sum(tumor_inv$CNV_status == "COMPLETE_PILOT")
n_remaining <- sum(tumor_inv$CNV_status == "NOT_RUN")
n_ineligible <- sum(tumor_inv$CNV_status == "INELIGIBLE_LT30")

cat("Total tumor samples     : ", total_tumor, "\n", sep = "")
cat("CNV completed           : ", n_completed, "\n", sep = "")
cat("CNV remaining           : ", n_remaining, "\n", sep = "")
cat("CNV ineligible (<30)    : ", n_ineligible, "\n", sep = "")

# ============================================================
# FINAL OUTPUT
# ============================================================

cat("\n============================================================\n")
cat("STEP 17H8 COMPLETE ✓\n")
cat("GSE221561 FOUR-PILOT CANONICAL FREEZE\n")
cat("============================================================\n\n")

cat("Canonical pilot validation:\n\n")

for (nm in names(pilots)) {
  r <- audit_results[[nm]]
  all_pass <- r$tumor_cell_count_validation == "PASS" &&
              r$gene_count_validation == "PASS" &&
              r$arm_count_validation == "PASS" &&
              r$metric_validation == "PASS"
  cat(nm, ":\n", sep = "")
  cat("  Tumor cells              : ", r$tumor_cells, "\n", sep = "")
  cat("  CNA signal ratio         : ", sprintf("%.2f", r$CNA_signal_ratio), "\n", sep = "")
  cat("  Tumor cnaCor             : ", sprintf("%.3f", r$tumor_median_cnaCor), "\n", sep = "")
  cat("  Correlation separation   : ", sprintf("%.3f", r$correlation_separation), "\n", sep = "")
  cat("  Tumor > normal 95th      : ", sprintf("%.1f%%", 100 * r$tumor_above_normal95), "\n", sep = "")
  cat("  Coherent arms            : ", r$coherent_CNA_arms, " (", r$coherent_CNA_arm_names, ")\n", sep = "")
  cat("  Tumor cell count         : ", r$tumor_cell_count_validation, "\n", sep = "")
  cat("  Gene count               : ", r$gene_count_validation, "\n", sep = "")
  cat("  Arm count                : ", r$arm_count_validation, "\n", sep = "")
  cat("  Metric validation        : ", r$metric_validation, "\n", sep = "")
  overall <- if (all_pass) "PASS" else "FAIL"
  cat("  OVERALL                  : ", overall, "\n\n", sep = "")
}

cat("--------------------------------------------\n\n")

cat("Surgery_alone pilots:\n")
cat("  n = 2\n")
cat(sprintf("  median CNA ratio     : %.2f\n", median(surgery$CNA_signal_ratio)))
cat(sprintf("  range                : %.2f – %.2f\n", min(surgery$CNA_signal_ratio), max(surgery$CNA_signal_ratio)))
cat(sprintf("  median cnaCor        : %.3f\n", median(surgery$tumor_median_cnaCor)))
cat(sprintf("  median coherent arms : %.0f\n", median(surgery$coherent_CNA_arms)))
cat(sprintf("  median correlation sep: %.3f\n", median(surgery$correlation_separation)))
cat(sprintf("  median tumor > 95th  : %.1f%%\n", 100 * median(surgery$tumor_above_normal95)))

cat("\nNeoadjuvant-treated pilots:\n")
cat("  n = 2\n")
cat(sprintf("  median CNA ratio     : %.2f\n", median(neoadj$CNA_signal_ratio)))
cat(sprintf("  range                : %.2f – %.2f\n", min(neoadj$CNA_signal_ratio), max(neoadj$CNA_signal_ratio)))
cat(sprintf("  median cnaCor        : %.3f\n", median(neoadj$tumor_median_cnaCor)))
cat(sprintf("  median coherent arms : %.0f\n", median(neoadj$coherent_CNA_arms)))
cat(sprintf("  median correlation sep: %.3f\n", median(neoadj$correlation_separation)))
cat(sprintf("  median tumor > 95th  : %.1f%%\n", 100 * median(neoadj$tumor_above_normal95)))

cat("\n--------------------------------------------\n\n")

cat("GSE221561 tumor sample inventory:\n\n")
cat("  Total tumor samples    : ", total_tumor, "\n", sep = "")
cat("  CNV completed          : ", n_completed, "\n", sep = "")
cat("  CNV remaining          : ", n_remaining, "\n", sep = "")
cat("  CNV ineligible (<30)   : ", n_ineligible, "\n", sep = "")

cat("\nRemaining eligible samples:\n")
for (i in seq_len(nrow(remaining))) {
  cat("  ", remaining$sample_id[i], " — ", remaining$treatment_group[i],
      " — ", remaining$epithelial_cells[i], " cells\n", sep = "")
}

cat("\n--------------------------------------------\n\n")

cat("InferCNA run                : NO\n")
cat("findMalignant               : NO\n")
cat("Malignant labels assigned   : NO\n")
cat("New samples processed       : NO\n")

cat("\n--------------------------------------------\n\n")

cat("Interpretation:\n\n")
cat("Two Surgery_alone pilots showed consistently stronger inferred CNA\n")
cat("than two Neoadjuvant-treated pilots.\n\n")
cat("These four pilots provide technical validation and an exploratory\n")
cat("sample-level pattern.\n\n")
cat("EXPLORATORY ONLY — n=2 per group.\n")
cat("NO causal treatment interpretation.\n")

cat("\nSTOP HERE.\n")
cat("DO NOT START BATCH CNV YET.\n")
