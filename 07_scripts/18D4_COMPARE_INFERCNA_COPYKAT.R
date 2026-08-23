#!/usr/bin/env Rscript
# ==============================================================================
# 18D4: ORTHOGONAL COMPARISON - InferCNA vs CopyKAT
# ==============================================================================

library(data.table)

cat("============================================================\n")
cat("STEP 18D4: ORTHOGONAL CNV COMPARISON\n")
cat("============================================================\n\n")

OUTDIR <- "06_tables/cross_dataset/malignant_epithelial/GSE160269_copykat_validation"

# ============================================================
# 1. LOAD INFERCNA RESULTS
# ============================================================

cat("1. Loading InferCNA results...\n")

# Load sensitivity metrics
infercna_metrics <- fread("06_tables/cross_dataset/malignant_epithelial/GSE160269_reference_sensitivity/GSE160269_ALL_SENSITIVITY_METRICS.csv")

# Get FULL_REFERENCE results for both tumors
infercna_full <- infercna_metrics[config == "FULL_REFERENCE"]

cat("   InferCNA results loaded\n\n")

# ============================================================
# 2. LOAD COPYKAT RESULTS
# ============================================================

cat("2. Loading CopyKAT results...\n")

copykat_results <- list()

for (tumor in c("P16T", "P39T")) {
  pred_file <- file.path(OUTDIR, tumor, paste0(tumor, "_copykat_predictions.csv"))
  annot_file <- file.path(OUTDIR, paste0(tumor, "_cell_annotations.csv"))
  
  if (!file.exists(pred_file)) {
    cat(sprintf("   WARNING: %s predictions not found\n", tumor))
    next
  }
  
  predictions <- fread(pred_file)
  annotations <- fread(annot_file)
  
  # CopyKAT may modify cell barcodes (replace "-" with ".")
  # Normalize barcodes for matching
  normalize_barcodes <- function(x) gsub("-", ".", x)
  predictions[, cell_barcode_norm := normalize_barcodes(cell_barcode)]
  annotations[, cell_barcode_norm := normalize_barcodes(cell_barcode)]
  
  # Count tumor cells only
  tumor_cells <- annotations[role == "TUMOR", cell_barcode_norm]
  tumor_preds <- predictions[cell_barcode_norm %in% tumor_cells]
  
  # Count normal cells
  normal_cells <- annotations[role == "NORMAL", cell_barcode_norm]
  normal_preds <- predictions[cell_barcode_norm %in% normal_cells]
  
  copykat_results[[tumor]] <- list(
    predictions = predictions,
    annotations = annotations,
    tumor_cells = tumor_cells,
    normal_cells = normal_cells,
    tumor_preds = tumor_preds,
    normal_preds = normal_preds
  )
  
  cat(sprintf("   %s: %d tumor cells, %d normal cells\n", 
              tumor, length(tumor_cells), length(normal_cells)))
}

cat("\n")

# ============================================================
# 3. COMPUTE METRICS FOR EACH TUMOR
# ============================================================

cat("3. Computing orthogonal comparison metrics...\n\n")

comparison_results <- list()

for (tumor in c("P16T", "P39T")) {
  cat(sprintf("--- %s ---\n", tumor))
  
  # InferCNA metrics
  inf <- infercna_full[tumor_id == tumor]
  
  if (nrow(inf) == 0) {
    cat(sprintf("   WARNING: No InferCNA results for %s\n", tumor))
    next
  }
  
  # CopyKAT metrics
  ck <- copykat_results[[tumor]]
  
  if (is.null(ck)) {
    cat(sprintf("   WARNING: No CopyKAT results for %s\n", tumor))
    next
  }
  
  # Tumor cell predictions
  tumor_preds <- ck$tumor_preds$copykat.pred
  n_aneuploid <- sum(tumor_preds == "aneuploid", na.rm = TRUE)
  n_diploid <- sum(tumor_preds == "diploid", na.rm = TRUE)
  n_not_defined <- sum(tumor_preds == "not.defined", na.rm = TRUE)
  n_total <- length(tumor_preds)
  
  aneuploid_frac <- n_aneuploid / n_total
  not_defined_frac <- n_not_defined / n_total
  
  # Normal cell predictions (QC)
  normal_preds <- ck$normal_preds$copykat.pred
  n_normal_aneuploid <- sum(normal_preds == "aneuploid", na.rm = TRUE)
  n_normal_total <- length(normal_preds)
  normal_aneuploid_frac <- n_normal_aneuploid / n_normal_total
  
  cat(sprintf("   InferCNA:\n"))
  cat(sprintf("     Ratio: %.3f\n", inf$cna_signal_ratio))
  cat(sprintf("     cnaCor: %.3f\n", inf$tumor_cnaCor))
  cat(sprintf("     Coherent arms: %d\n", inf$coherent_arm_count))
  cat(sprintf("     Reference status: %s\n", "REFERENCE_STABLE_WEAK"))
  
  cat(sprintf("   CopyKAT:\n"))
  cat(sprintf("     Aneuploid: %d / %d (%.1f%%)\n", n_aneuploid, n_total, aneuploid_frac * 100))
  cat(sprintf("     Diploid: %d / %d (%.1f%%)\n", n_diploid, n_total, (n_diploid/n_total) * 100))
  cat(sprintf("     Not.defined: %d / %d (%.1f%%)\n", n_not_defined, n_total, not_defined_frac * 100))
  cat(sprintf("     Known-normal aneuploid: %d / %d (%.1f%%)\n", 
              n_normal_aneuploid, n_normal_total, normal_aneuploid_frac * 100))
  
  # Classification
  if (normal_aneuploid_frac > 0.20) {
    orthogonal_status <- "COPYKAT_REFERENCE_QC_FAIL"
    cat(sprintf("   *** COPYKAT REFERENCE QC FAIL: %.1f%% known-normal aneuploid ***\n", 
                normal_aneuploid_frac * 100))
  } else if (inf$cna_signal_ratio < 3.0 && inf$coherent_arm_count == 0 && 
             aneuploid_frac < 0.30) {
    orthogonal_status <- "ORTHOGONAL_WEAK_CONCORDANT"
  } else if (inf$cna_signal_ratio < 3.0 && aneuploid_frac >= 0.50) {
    orthogonal_status <- "ORTHOGONAL_STRONG_DISCORDANT"
  } else if (not_defined_frac > 0.50) {
    orthogonal_status <- "COPYKAT_INCONCLUSIVE"
  } else {
    orthogonal_status <- "MIXED_EVIDENCE"
  }
  
  cat(sprintf("   Orthogonal status: %s\n\n", orthogonal_status))
  
  comparison_results[[tumor]] <- data.table(
    sample = tumor,
    tumor_cells = n_total,
    infercna_ratio = inf$cna_signal_ratio,
    infercna_cnaCor = inf$tumor_cnaCor,
    infercna_coherent_arms = inf$coherent_arm_count,
    infercna_reference_status = "REFERENCE_STABLE_WEAK",
    copykat_aneuploid = n_aneuploid,
    copykat_diploid = n_diploid,
    copykat_not_defined = n_not_defined,
    copykat_aneuploid_candidate_fraction = aneuploid_frac,
    copykat_not_defined_fraction = not_defined_frac,
    known_normal_aneuploid_fraction = normal_aneuploid_frac,
    orthogonal_status = orthogonal_status
  )
}

# ============================================================
# 4. SAVE FINAL COMPARISON TABLE
# ============================================================

cat("4. Saving final comparison table...\n")

final_dt <- rbindlist(comparison_results)
fwrite(final_dt, file.path(OUTDIR, "STEP18D_P16T_P39T_orthogonal_CNV_validation.csv"))
cat("   Saved: STEP18D_P16T_P39T_orthogonal_CNV_validation.csv\n\n")

# ============================================================
# 5. FINAL OUTPUT
# ============================================================

cat("\n============================================\n")
cat("STEP 18D COMPLETE\n")
cat("GSE160269 ORTHOGONAL CNV VALIDATION\n")
cat("============================================\n\n")

for (tumor in c("P16T", "P39T")) {
  row <- final_dt[sample == tumor]
  
  cat(sprintf("%s\n", tumor))
  cat(sprintf("InferCNA ratio              : %.3f\n", row$infercna_ratio))
  cat(sprintf("InferCNA cnaCor             : %.3f\n", row$infercna_cnaCor))
  cat(sprintf("InferCNA coherent arms      : %d\n", row$infercna_coherent_arms))
  cat(sprintf("InferCNA reference status   : %s\n", row$infercna_reference_status))
  cat(sprintf("CopyKAT aneuploid           : %d\n", row$copykat_aneuploid))
  cat(sprintf("CopyKAT diploid             : %d\n", row$copykat_diploid))
  cat(sprintf("CopyKAT not-defined         : %d\n", row$copykat_not_defined))
  cat(sprintf("Aneuploid candidate %%       : %.1f%%\n", row$copykat_aneuploid_candidate_fraction * 100))
  cat(sprintf("Known-normal aneuploid %%    : %.1f%%\n", row$known_normal_aneuploid_fraction * 100))
  cat(sprintf("Orthogonal status           : %s\n", row$orthogonal_status))
  cat("--------------------------------------------\n\n")
}

cat("InferCNA rerun              : NO\n")
cat("findMalignant               : NO\n")
cat("Malignant labels assigned   : NO\n")
cat("Other tumor samples run     : NO\n")
cat("Cells removed               : NO\n\n")

cat("STOP HERE.\n")
cat("DO NOT RUN THE OTHER 58 GSE160269 TUMORS.\n")
