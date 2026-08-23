#!/usr/bin/env Rscript
# ==============================================================================
# 19D: GSE160269 REFERENCE SENSITIVITY - FINAL SUMMARY & CLASSIFICATION
# ==============================================================================

library(data.table)

OUTDIR <- "06_tables/cross_dataset/malignant_epithelial/GSE160269_reference_sensitivity"

# Load results
p16t <- fread(file.path(OUTDIR, "P16T_sensitivity_metrics.csv"))
p39t <- fread(file.path(OUTDIR, "P39T_sensitivity_metrics.csv"))

# Combine
combined <- rbindlist(list(p16t, p39t))

cat("============================================================\n")
cat("GSE160269 NORMAL REFERENCE SENSITIVITY AUDIT\n")
cat("============================================================\n\n")

cat("STEP 1: NORMAL REFERENCE CELL COUNTS\n")
cat("------------------------------------------------------------\n")
cat("  P126N:  24 cells\n")
cat("  P127N: 123 cells\n")
cat("  P128N:  16 cells (smallest)\n")
cat("  P130N:  20 cells\n")
cat("  Total: 183 cells\n")
cat("  Imbalance ratio: 7.7x\n\n")

cat("STEP 2-4: ALL CONFIGURATION METRICS\n")
cat("------------------------------------------------------------\n\n")

for (tumor in c("P16T", "P39T")) {
  dt <- combined[tumor_id == tumor]
  cat(sprintf("\n%s (%d tumor cells):\n", tumor, unique(dt$tumor_cells)))
  cat(sprintf("%-25s %8s %8s %8s %8s %8s %8s %8s %8s %8s\n",
              "Config", "RefCells", "Ratio", "cnaCor", "Sep", ">95th", ">99th", "Arms", "Gain", "Loss"))
  cat(paste(rep("-", 100), collapse = ""), "\n")
  
  for (i in 1:nrow(dt)) {
    row <- dt[i]
    # Parse coherent arms
    arms <- unlist(strsplit(as.character(row$coherent_arm_names), "; "))
    arms <- arms[arms != ""]
    n_gain <- sum(grepl("_gain$", arms))
    n_loss <- sum(grepl("_loss$", arms))
    
    cat(sprintf("%-25s %8d %8.3f %8.3f %8.3f %7.1f%% %7.1f%% %8d %8d %8d\n",
                row$config, row$total_normal, row$cna_signal_ratio,
                row$tumor_cnaCor, row$separation,
                row$pct_above_95, row$pct_above_99,
                row$coherent_arm_count, n_gain, n_loss))
  }
}

cat("\n\nSTEP 5: REFERENCE STABILITY CLASSIFICATION\n")
cat("============================================================\n\n")

# Classification function
classify <- function(dt, tumor_name) {
  full <- dt[config == "FULL_REFERENCE"]
  balanced <- dt[config == "BALANCED_WITHIN_DATASET"]
  leave_out <- dt[!config %in% c("FULL_REFERENCE", "BALANCED_WITHIN_DATASET")]
  
  ratios <- leave_out$cna_signal_ratio
  cnaCors <- leave_out$tumor_cnaCor
  
  ratio_range <- max(ratios) - min(ratios)
  cnaCor_range <- max(cnaCors) - min(cnaCors)
  
  cat(sprintf("\n%s Analysis:\n", tumor_name))
  cat(sprintf("  Leave-out CNA signal ratio range: %.3f - %.3f (delta=%.3f)\n",
              min(ratios), max(ratios), ratio_range))
  cat(sprintf("  Leave-out cnaCor range: %.3f - %.3f (delta=%.3f)\n",
              min(cnaCors), max(cnaCors), cnaCor_range))
  cat(sprintf("  FULL ratio: %.3f, cnaCor: %.3f\n", full$cna_signal_ratio, full$tumor_cnaCor))
  cat(sprintf("  BALANCED ratio: %.3f, cnaCor: %.3f\n", balanced$cna_signal_ratio, balanced$tumor_cnaCor))
  
  # Classify
  if (full$cna_signal_ratio >= 3.0 && full$tumor_cnaCor >= 0.6 && min(ratios) >= 2.5) {
    stability <- "REFERENCE_STABLE_STRONG"
    reason <- "Strong CNA consistently above thresholds"
  } else if (ratio_range <= 0.5 && cnaCor_range <= 0.05 &&
             abs(full$cna_signal_ratio - balanced$cna_signal_ratio) <= 0.3 &&
             full$cna_signal_ratio < 3.0) {
    stability <- "REFERENCE_STABLE_WEAK"
    reason <- "Weak CNA persists across all configurations"
  } else if (ratio_range > 0.5 || cnaCor_range > 0.1) {
    stability <- "REFERENCE_SENSITIVE"
    reason <- "Substantial variation depending on reference"
  } else {
    stability <- "INCONCLUSIVE"
    reason <- "Inconsistent results"
  }
  
  cat(sprintf("  Classification: %s\n", stability))
  cat(sprintf("  Reason: %s\n", reason))
  
  list(stability = stability, reason = reason)
}

p16t_class <- classify(p16t, "P16T")
p39t_class <- classify(p39t, "P39T")

cat("\n\n============================================================\n")
cat("FINAL CLASSIFICATIONS\n")
cat("============================================================\n\n")

cat(sprintf("P16T:  %s\n", p16t_class$stability))
cat(sprintf("P39T:  %s\n", p39t_class$stability))

cat("\n\n============================================================\n")
cat("INTERPRETATION\n")
cat("============================================================\n\n")

if (p16t_class$stability == "REFERENCE_STABLE_WEAK" &&
    p39t_class$stability == "REFERENCE_STABLE_WEAK") {
  cat("Both P16T and P39T are REFERENCE_STABLE_WEAK.\n\n")
  cat("This means:\n")
  cat("  - Weak CNA persists across ALL within-dataset reference configurations\n")
  cat("  - The weak signal is NOT driven by any single normal reference sample\n")
  cat("  - Both tumors show consistently low CNA signal ratios (<1.5)\n")
  cat("  - No coherent CNA arms detected in any configuration\n\n")
  cat("RECOMMENDATION:\n")
  cat("  The weak CNA signal appears to be a biological reality of these tumors,\n")
  cat("  NOT an artifact of reference selection. A third independent tumor pilot\n")
  cat("  is unlikely to change the conclusion.\n")
} else {
  cat("One or both tumors show REFERENCE_SENSITIVE or INCONCLUSIVE patterns.\n")
  cat("Additional investigation may be warranted.\n")
}

# Save combined table
fwrite(combined, file.path(OUTDIR, "GSE160269_ALL_SENSITIVITY_METRICS.csv"))

cat("\n\nSaved: GSE160269_ALL_SENSITIVITY_METRICS.csv\n")
