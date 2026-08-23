# ============================================================
# ESCC Neoadjuvant Project
# STEP 17H9A
#
# Build batch plan for remaining GSE221561 CNV samples.
#
# Reads Step 17H8 status CSV.
# Identifies CNV_status == "NOT_RUN" samples.
# ============================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
})

cat("\n============================================================\n")
cat("STEP 17H9A: BUILD REMAINING CNV BATCH PLAN\n")
cat("============================================================\n\n")

# ============================================================
# 1. READ STEP 17H8 STATUS
# ============================================================

status_file <- "06_tables/cross_dataset/malignant_epithelial/GSE221561_four_pilot_freeze/STEP17H8_GSE221561_all_tumor_sample_CNV_status.csv"

if (!file.exists(status_file)) stop("Missing: ", status_file)

status <- fread(status_file)
cat("Total tumor samples: ", nrow(status), "\n\n")

print(status[, .(sample_id, treatment_group, epithelial_cells, CNV_status)])

# ============================================================
# 2. FILTER NOT_RUN
# ============================================================

cat("\n===== 2. FILTER CNV_status == NOT_RUN =====\n\n")

remaining <- status[CNV_status == "NOT_RUN"]

cat("Remaining NOT_RUN samples: ", nrow(remaining), "\n\n")

# ============================================================
# 3. VALIDATION
# ============================================================

cat("===== 3. VALIDATION =====\n\n")

# Must be exactly 5
if (nrow(remaining) != 5) {
  stop("Expected 5 NOT_RUN samples, found ", nrow(remaining))
}
cat("Count check: 5 ✓\n")

# All must be Neoadjuvant_treated
if (any(remaining$treatment_group != "Neoadjuvant_treated")) {
  stop("Non-Neoadjuvant_treated sample found in NOT_RUN")
}
cat("All Neoadjuvant_treated ✓\n")

# No Surgery_alone in NOT_RUN
surg_in_remaining <- remaining[treatment_group == "Surgery_alone"]
if (nrow(surg_in_remaining) > 0) {
  stop("Surgery_alone sample still NOT_RUN: ", paste(surg_in_remaining$sample_id, collapse = ", "))
}
cat("No Surgery_alone in NOT_RUN ✓\n")

# All epithelial >= 30
if (any(remaining$epithelial_cells < 30)) {
  low <- remaining[epithelial_cells < 30]
  cat("\nWARNING: Samples with <30 epithelial cells:\n")
  print(low)
  # Mark as INELIGIBLE
  remaining[epithelial_cells < 30, planned_status := "INELIGIBLE_LT30"]
}
cat("All epithelial >= 30 ✓\n")

# ============================================================
# 4. BUILD PLAN
# ============================================================

cat("\n===== 4. BUILD BATCH PLAN =====\n\n")

# Only READY_TO_RUN samples (>=30 cells)
ready <- remaining[epithelial_cells >= 30]
setorder(ready, -epithelial_cells)

plan <- data.table(
  run_order = seq_len(nrow(ready)),
  sample = ready$sample_id,
  treatment = ready$treatment_group,
  epithelial_cells = ready$epithelial_cells,
  pre_status = ready$CNV_status,
  planned_status = "READY_TO_RUN"
)

# Add any INELIGIBLE
ineligible <- remaining[epithelial_cells < 30]
if (nrow(ineligible) > 0) {
  inel_plan <- data.table(
    run_order = NA_integer_,
    sample = ineligible$sample_id,
    treatment = ineligible$treatment_group,
    epithelial_cells = ineligible$epithelial_cells,
    pre_status = ineligible$CNV_status,
    planned_status = "INELIGIBLE_LT30"
  )
  plan <- rbind(plan, inel_plan)
}

cat("Batch plan:\n\n")
print(plan)

# ============================================================
# 5. SAVE PLAN
# ============================================================

outdir <- "06_tables/cross_dataset/malignant_epithelial/GSE221561_CNV_batch"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

plan_file <- file.path(outdir, "STEP17H9_batch_plan.csv")
fwrite(plan, plan_file, bom = TRUE)

cat("\nSaved: ", plan_file, "\n", sep = "")

# ============================================================
# 6. FINAL SUMMARY
# ============================================================

cat("\n============================================================\n")
cat("STEP 17H9A COMPLETE\n")
cat("============================================================\n\n")

cat("Total remaining     : ", nrow(remaining), "\n", sep = "")
cat("READY_TO_RUN        : ", sum(plan$planned_status == "READY_TO_RUN"), "\n", sep = "")
cat("INELIGIBLE_LT30     : ", sum(plan$planned_status == "INELIGIBLE_LT30"), "\n", sep = "")

cat("\nRun order (by epithelial_cells desc):\n\n")
for (i in seq_len(nrow(plan))) {
  if (plan$planned_status[i] == "READY_TO_RUN") {
    cat(sprintf("  %d. %s — %d cells — %s\n",
                plan$run_order[i], plan$sample[i],
                plan$epithelial_cells[i], plan$treatment[i]))
  } else {
    cat(sprintf("  SKIP. %s — %d cells — %s\n",
                plan$sample[i], plan$epithelial_cells[i], plan$planned_status[i]))
  }
}

cat("\nSTOP HERE.\n")
cat("Review plan before running 17H9C.\n")
