# ============================================================
# ESCC Neoadjuvant Project
# STEP 17H9C
#
# Sequential controller for remaining GSE221561 CNV samples.
#
# Reads STEP17H9_batch_plan.csv.
# Runs 17H9B worker for each READY_TO_RUN sample, one at a time.
# Crash/restart safe: skips already-complete samples.
#
# Usage:
#   Rscript 17H9C_RUN_REMAINING_GSE221561_CNV_BATCH.R
# ============================================================

options(stringsAsFactors = FALSE, timeout = 7200)

suppressPackageStartupMessages({
  library(data.table)
})

cat("\n============================================================\n")
cat("STEP 17H9C: REMAINING GSE221561 CNV BATCH CONTROLLER\n")
cat("============================================================\n\n")

# ============================================================
# 1. READ BATCH PLAN
# ============================================================

plan_file <- "06_tables/cross_dataset/malignant_epithelial/GSE221561_CNV_batch/STEP17H9_batch_plan.csv"

if (!file.exists(plan_file)) stop("Missing batch plan: ", plan_file)

plan <- fread(plan_file)
ready <- plan[planned_status == "READY_TO_RUN"]
setorder(ready, run_order)

cat("Total plan rows     : ", nrow(plan), "\n", sep = "")
cat("READY_TO_RUN        : ", nrow(ready), "\n", sep = "")
cat("INELIGIBLE_LT30     : ", sum(plan$planned_status == "INELIGIBLE_LT30"), "\n", sep = "")

if (nrow(ready) == 0) {
  cat("No samples to run.\n")
  stop("Nothing to do.")
}

cat("\nRun order:\n")
for (i in seq_len(nrow(ready))) {
  cat(sprintf("  %d. %s — %d cells\n", ready$run_order[i], ready$sample[i], ready$epithelial_cells[i]))
}

# ============================================================
# 2. LOG SETUP
# ============================================================

log_dir <- "08_logs"
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

controller_log <- file.path(log_dir, "17H9_REMAINING_GSE221561_CNV_BATCH.log")
cat(paste0("\n=== CONTROLLER START: ", Sys.time(), " ===\n"), file = controller_log)

# ============================================================
# 3. SEQUENTIAL EXECUTION
# ============================================================

worker_script <- "07_scripts/17H9B_RUN_ONE_GSE221561_CNV_SAMPLE.R"

results_base <- "06_tables/cross_dataset/malignant_epithelial"

completed_count <- 0L
failed_count <- 0L
skipped_count <- 0L

for (i in seq_len(nrow(ready))) {
  sid <- ready$sample[i]
  cat("\n============================================================\n")
  cat("SAMPLE ", i, " of ", nrow(ready), ": ", sid, "\n", sep = "")
  cat("============================================================\n\n")

  # --- Crash/restart safety: check if already complete ---
  manifest_path <- file.path(
    results_base,
    paste0("GSE221561_", sid, "_pilot"),
    "17H9_results",
    paste0("GSE221561_", sid, "_17H9_manifest.csv")
  )

  if (file.exists(manifest_path)) {
    existing <- fread(manifest_path)
    if (nrow(existing) > 0 && existing$validation_status[1] == "PASS") {
      cat("SKIP_ALREADY_COMPLETE: ", sid, "\n", sep = "")
      cat(paste0("SKIP: ", sid, " at ", Sys.time(), "\n"), file = controller_log, append = TRUE)
      skipped_count <- skipped_count + 1L

      # Print summary of skipped
      cat("============================================\n")
      cat("17H9 SAMPLE COMPLETE (SKIPPED — ALREADY DONE)\n")
      cat("Sample: ", sid, "\n", sep = "")
      cat("Tumor cells: ", existing$tumor_cells, "\n", sep = "")
      cat("Signal ratio: ", sprintf("%.2f", existing$CNA_signal_ratio), "\n", sep = "")
      cat("Tumor cnaCor: ", sprintf("%.3f", existing$tumor_median_cnaCor), "\n", sep = "")
      cat("Coherent arms: ", existing$coherent_CNA_arms, "\n", sep = "")
      cat("Validation: ", existing$validation_status, "\n", sep = "")
      cat("============================================\n")
      next
    } else if (nrow(existing) > 0 && existing$validation_status[1] != "PASS") {
      stop("Existing manifest for ", sid, " has validation_status != PASS. STOPPING.")
    }
  }

  # --- Run worker ---
  sample_log <- file.path(log_dir, paste0("17H9_", sid, "_CNV.log"))
  cat(paste0("=== START: ", sid, " at ", Sys.time(), " ===\n"), file = sample_log)

  cmd <- paste("Rscript", worker_script, sid, ">", sample_log, "2>&1")
  cat("Running: ", cmd, "\n\n", sep = "")

  exit_code <- system(cmd)

  if (exit_code != 0) {
    cat("\n*** WORKER FAILED FOR ", sid, " (exit code ", exit_code, ") ***\n", sep = "")
    cat(paste0("FAIL: ", sid, " at ", Sys.time(), " exit=", exit_code, "\n"), file = controller_log, append = TRUE)
    failed_count <- failed_count + 1L
    stop("WORKER FAILED FOR ", sid, ". STOPPING ENTIRE BATCH.")
  }

  # --- Verify manifest was created ---
  if (!file.exists(manifest_path)) {
    stop("Worker completed but manifest not found: ", manifest_path)
  }

  result <- fread(manifest_path)
  if (result$validation_status[1] != "PASS") {
    stop("Worker completed but validation_status != PASS for ", sid)
  }

  cat(paste0("PASS: ", sid, " at ", Sys.time(), "\n"), file = controller_log, append = TRUE)
  completed_count <- completed_count + 1L

  # --- Print completion summary ---
  cat("\n============================================\n")
  cat("17H9 SAMPLE COMPLETE\n")
  cat("Sample: ", sid, "\n", sep = "")
  cat("Tumor cells: ", result$tumor_cells, "\n", sep = "")
  cat("CNA genes: ", result$CNA_genes, "\n", sep = "")
  cat("Signal ratio: ", sprintf("%.2f", result$CNA_signal_ratio), "\n", sep = "")
  cat("Tumor cnaCor: ", sprintf("%.3f", result$tumor_median_cnaCor), "\n", sep = "")
  cat("Correlation separation: ", sprintf("%.3f", result$correlation_separation), "\n", sep = "")
  cat("Tumor > normal 95th: ", sprintf("%.1f%%", 100 * result$tumor_above_normal95), "\n", sep = "")
  cat("Coherent arms: ", result$coherent_CNA_arms, "\n", sep = "")
  cat("Validation: ", result$validation_status, "\n", sep = "")
  cat("findMalignant: NO\n")
  cat("Malignant labels: NO\n")
  cat("============================================\n")

  # --- Cleanup between samples ---
  gc(verbose = FALSE)
}

# ============================================================
# 4. FINAL CONTROLLER SUMMARY
# ============================================================

cat("\n============================================================\n")
cat("STEP 17H9C CONTROLLER FINISHED\n")
cat("============================================================\n\n")

cat("Completed  : ", completed_count, "\n", sep = "")
cat("Skipped    : ", skipped_count, "\n", sep = "")
cat("Failed     : ", failed_count, "\n", sep = "")
cat("Total ready: ", nrow(ready), "\n", sep = "")

cat(paste0("=== CONTROLLER END: ", Sys.time(), " ===\n"), file = controller_log, append = TRUE)

if (failed_count > 0) {
  stop("BATCH COMPLETED WITH FAILURES.")
}

cat("\nAll READY_TO_RUN samples completed.\n")
cat("Now run STEP 17H9D to freeze all 9 tumor samples.\n")
