#!/usr/bin/env Rscript
# ==============================================================================
# 18D3: RUN COPYKAT ON P39T
# ==============================================================================

library(copykat)
library(data.table)

cat("============================================================\n")
cat("STEP 18D3: COPYKAT P39T\n")
cat("============================================================\n\n")

OUTDIR <- "06_tables/cross_dataset/malignant_epithelial/GSE160269_copykat_validation/P39T"
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

BASEDIR <- "06_tables/cross_dataset/malignant_epithelial/GSE160269_copykat_validation"

# ============================================================
# 1. LOAD DATA
# ============================================================

cat("1. Loading P39T raw counts...\n")

raw_counts <- readRDS(file.path(BASEDIR, "P39T_raw_counts.rds"))
annot <- fread(file.path(BASEDIR, "P39T_cell_annotations.csv"))

cat("   Raw count matrix:", nrow(raw_counts), "genes x", ncol(raw_counts), "cells\n")
cat("   Tumor cells:", sum(annot$role == "TUMOR"), "\n")
cat("   Normal cells:", sum(annot$role == "NORMAL"), "\n\n")

# ============================================================
# 2. PREPARE COPYKAT INPUT
# ============================================================

cat("2. Preparing CopyKAT input...\n")

# CopyKAT expects genes as rows, cells as columns
# Gene symbols as rownames
# Raw UMI counts (NOT normalized)

# Check that rownames are gene symbols
cat("   First 10 gene symbols:", head(rownames(raw_counts), 10), "\n")
cat("   Gene count:", nrow(raw_counts), "\n\n")

# ============================================================
# 3. RUN COPYKAT
# ============================================================

cat("3. Running CopyKAT (this may take several minutes)...\n\n")

start_time <- Sys.time()

# CopyKAT with default parameters
copykat_result <- tryCatch({
  copykat(
    raw.data = raw_counts,
    id.type = "S",
    cell.line = "no",
    ngene.chr = 5,
    LOW.DR = 0.05,
    UP.DR = 0.1,
    win.size = 25,
    norm.cell.names = colnames(raw_counts)[annot$role == "NORMAL"],
    sam.name = "GSE160269_P39T",
    n.cores = 1,
    output.seg = FALSE
  )
}, error = function(e) {
  cat("\n============================================\n")
  cat("COPYKAT ERROR\n")
  cat("============================================\n")
  cat(conditionMessage(e), "\n")
  stop("Step 18D3 stopped during CopyKAT.")
})

end_time <- Sys.time()
runtime_min <- as.numeric(difftime(end_time, start_time, units = "mins"))
cat(sprintf("\nCopyKAT runtime: %.2f min\n", runtime_min))

# ============================================================
# 4. EXTRACT RESULTS
# ============================================================

cat("\n4. Extracting CopyKAT results...\n")

# Get prediction results
predictions <- copykat_result$prediction

if (is.null(predictions)) {
  cat("ERROR: No predictions returned\n")
  stop("CopyKAT returned no predictions")
}

cat("   Prediction table dimensions:", nrow(predictions), "x", ncol(predictions), "\n")
cat("   Prediction columns:", colnames(predictions), "\n\n")

# Count predictions
pred_table <- table(predictions$copykat.pred)
cat("Prediction counts:\n")
for (pred_name in names(pred_table)) {
  cat(sprintf("   %s: %d\n", pred_name, pred_table[pred_name]))
}

# ============================================================
# 5. SAVE RESULTS
# ============================================================

cat("\n5. Saving results...\n")

# Save full CopyKAT object
saveRDS(copykat_result, file.path(OUTDIR, "P39T_copykat_object.rds"))
cat("   Saved: P39T_copykat_object.rds\n")

# Save predictions
fwrite(as.data.table(predictions, keep.rownames = "cell_barcode"), 
       file.path(OUTDIR, "P39T_copykat_predictions.csv"))
cat("   Saved: P39T_copykat_predictions.csv\n")

# Save CNA matrix if available
if (!is.null(copykat_result$CNV.matrix)) {
  saveRDS(copykat_result$CNV.matrix, file.path(OUTDIR, "P39T_copykat_CNA_matrix.rds"))
  cat("   Saved: P39T_copykat_CNA_matrix.rds\n")
}

# Create summary
summary_dt <- data.table(
  sample = "P39T",
  total_cells = nrow(predictions),
  aneuploid = sum(predictions$copykat.pred == "aneuploid", na.rm = TRUE),
  diploid = sum(predictions$copykat.pred == "diploid", na.rm = TRUE),
  not_defined = sum(predictions$copykat.pred == "not.defined", na.rm = TRUE),
  aneuploid_fraction = mean(predictions$copykat.pred == "aneuploid", na.rm = TRUE),
  diploid_fraction = mean(predictions$copykat.pred == "diploid", na.rm = TRUE),
  not_defined_fraction = mean(predictions$copykat.pred == "not.defined", na.rm = TRUE),
  runtime_minutes = runtime_min
)
fwrite(summary_dt, file.path(OUTDIR, "P39T_copykat_summary.csv"))
cat("   Saved: P39T_copykat_summary.csv\n")

# ============================================================
# 6. SUMMARY
# ============================================================

cat("\n============================================================\n")
cat("STEP 18D3 COPYKAT P39T COMPLETE\n")
cat("============================================================\n\n")

cat("Total cells:", summary_dt$total_cells, "\n")
cat("Aneuploid:", summary_dt$aneuploid, sprintf("(%.1f%%)", summary_dt$aneuploid_fraction * 100), "\n")
cat("Diploid:", summary_dt$diploid, sprintf("(%.1f%%)", summary_dt$diploid_fraction * 100), "\n")
cat("Not.defined:", summary_dt$not_defined, sprintf("(%.1f%%)", summary_dt$not_defined_fraction * 100), "\n")
cat("Runtime:", sprintf("%.2f min", runtime_min), "\n\n")

cat("Files saved to:", OUTDIR, "\n")
