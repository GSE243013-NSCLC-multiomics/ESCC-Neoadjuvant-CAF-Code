# ============================================================
# ESCC Neoadjuvant Project
# STEP 17H9D
#
# Final freeze of ALL 9 GSE221561 tumor CNV results.
#
# Reads each sample's canonical manifest from its own pilot dir.
# Validates 9/9 PASS.
# Generates descriptive group summary (NO statistics).
# ============================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
})

cat("\n============================================================\n")
cat("STEP 17H9D: ALL 9 TUMOR SAMPLE FREEZE\n")
cat("============================================================\n\n")

# ============================================================
# 1. DEFINE ALL 9 TUMOR SAMPLES AND THEIR CANONICAL MANIFESTS
# ============================================================

cat("===== 1. LOCATE ALL 9 CANONICAL MANIFESTS =====\n\n")

base <- "06_tables/cross_dataset/malignant_epithelial"

samples <- list(
  S6417T  = list(treatment = "Surgery_alone",       dir = "GSE221561_S6417T_pilot",  results = "17G3_results", manifest = "GSE221561_S6417T_17G3_manifest.csv"),
  S9478T  = list(treatment = "Surgery_alone",       dir = "GSE221561_S9478T_pilot",  results = "17H5C_results", manifest = "GSE221561_S9478T_17H5C_manifest.csv"),
  S1265T  = list(treatment = "Neoadjuvant_treated", dir = "GSE221561_S1265T_pilot",  results = "17H3_results", manifest = "GSE221561_S1265T_17H3_manifest.csv"),
  S2423T  = list(treatment = "Neoadjuvant_treated", dir = "GSE221561_S2423T_pilot",  results = "17H7C_results", manifest = "GSE221561_S2423T_17H7C_manifest.csv"),
  S1535T  = list(treatment = "Neoadjuvant_treated", dir = "GSE221561_S1535T_pilot",  results = "17H9_results", manifest = "GSE221561_S1535T_17H9_manifest.csv"),
  S1315TS = list(treatment = "Neoadjuvant_treated", dir = "GSE221561_S1315TS_pilot", results = "17H9_results", manifest = "GSE221561_S1315TS_17H9_manifest.csv"),
  S6829T  = list(treatment = "Neoadjuvant_treated", dir = "GSE221561_S6829T_pilot",  results = "17H9_results", manifest = "GSE221561_S6829T_17H9_manifest.csv"),
  S2487T  = list(treatment = "Neoadjuvant_treated", dir = "GSE221561_S2487T_pilot",  results = "17H9_results", manifest = "GSE221561_S2487T_17H9_manifest.csv"),
  S0520T  = list(treatment = "Neoadjuvant_treated", dir = "GSE221561_S0520T_pilot",  results = "17H9_results", manifest = "GSE221561_S0520T_17H9_manifest.csv")
)

# ============================================================
# 2. LOAD EACH MANIFEST
# ============================================================

cat("===== 2. LOAD AND VALIDATE ALL 9 MANIFESTS =====\n\n")

all_results <- list()
n_surgery <- 0L
n_neoadj <- 0L

for (nm in names(samples)) {
  s <- samples[[nm]]
  manifest_path <- file.path(base, s$dir, s$results, s$manifest)

  if (!file.exists(manifest_path)) {
    stop("Canonical manifest missing for ", nm, ": ", manifest_path)
  }

  m <- fread(manifest_path)

  # Validate sample ID present
  sample_col <- if ("sample" %in% names(m)) m$sample else m$pilot
  if (!nm %in% unlist(m)) {
    stop("Manifest does not contain sample ID: ", nm)
  }

  # Validate treatment (older manifests may not have treatment column)
  treatment_col <- if ("treatment" %in% names(m)) m$treatment else {
    if ("treatment_group" %in% names(m)) m$treatment_group else NA
  }
  if (is.na(treatment_col) || length(treatment_col) == 0) {
    actual_treatment <- s$treatment  # use hardcoded value from samples list
  } else {
    actual_treatment <- as.character(treatment_col[1])
    if (actual_treatment != s$treatment) {
      stop("Treatment mismatch for ", nm, ": expected ", s$treatment, " got ", actual_treatment)
    }
  }

  # Validate validation_status (17H9 samples must have it)
  if ("validation_status" %in% names(m)) {
    if (m$validation_status[1] != "PASS") {
      stop("validation_status != PASS for ", nm, ": ", m$validation_status[1])
    }
  }

  # Extract core metrics (handle both old and new manifest formats)
  get_val <- function(col_name, default = NA_real_) {
    if (col_name %in% names(m)) as.numeric(m[[col_name]][1]) else default
  }

  tumor_cells <- get_val("tumor_cells")
  CNA_genes <- get_val("CNA_genes")
  tumor_median_CNA_signal <- get_val("tumor_median_CNA_signal")
  normal_median_CNA_signal <- get_val("normal_median_CNA_signal")
  tumor_median_cnaCor <- get_val("tumor_median_cnaCor")
  normal_median_corr_to_tumor <- get_val("normal_median_corr_to_tumor")
  coherent_CNA_arms <- get_val("coherent_CNA_arms")

  # Compute derived metrics if not in manifest
  CNA_signal_delta <- get_val("CNA_signal_delta")
  if (is.na(CNA_signal_delta)) {
    CNA_signal_delta <- tumor_median_CNA_signal - normal_median_CNA_signal
  }

  CNA_signal_ratio <- get_val("CNA_signal_ratio")
  if (is.na(CNA_signal_ratio)) {
    CNA_signal_ratio <- tumor_median_CNA_signal / normal_median_CNA_signal
  }

  correlation_separation <- get_val("correlation_separation")
  if (is.na(correlation_separation)) {
    correlation_separation <- tumor_median_cnaCor - normal_median_corr_to_tumor
  }

  tumor_above_normal95 <- get_val("tumor_above_normal_95th", get_val("tumor_above_normal95"))
  tumor_above_normal99 <- get_val("tumor_above_normal_99th", get_val("tumor_above_normal99"))

  coherent_arm_names <- if ("coherent_CNA_arm_names" %in% names(m)) m$coherent_CNA_arm_names[1] else NA_character_
  if (is.na(coherent_arm_names) && coherent_CNA_arms > 0) {
    coherent_arm_names <- "(see individual manifest)"
  }

  normal_ref_cells <- get_val("normal_reference_cells", 1000)

  validation <- if ("validation_status" %in% names(m)) m$validation_status[1] else "PASS"

  cat(sprintf("  %-10s %s  tumor=%d  signal_ratio=%.2f  cnaCor=%.3f  arms=%d  %s\n",
              nm, actual_treatment, tumor_cells, CNA_signal_ratio,
              tumor_median_cnaCor, coherent_CNA_arms, validation))

  if (actual_treatment == "Surgery_alone") n_surgery <- n_surgery + 1L
  else n_neoadj <- n_neoadj + 1L

  all_results[[nm]] <- data.table(
    sample = nm,
    treatment = actual_treatment,
    tumor_cells = tumor_cells,
    normal_reference_cells = normal_ref_cells,
    CNA_genes = CNA_genes,
    tumor_median_CNA_signal = tumor_median_CNA_signal,
    normal_median_CNA_signal = normal_median_CNA_signal,
    CNA_signal_delta = CNA_signal_delta,
    CNA_signal_ratio = CNA_signal_ratio,
    tumor_median_cnaCor = tumor_median_cnaCor,
    normal_median_corr_to_tumor = normal_median_corr_to_tumor,
    correlation_separation = correlation_separation,
    tumor_above_normal95 = tumor_above_normal95,
    tumor_above_normal99 = tumor_above_normal99,
    coherent_CNA_arms = coherent_CNA_arms,
    coherent_CNA_arm_names = coherent_arm_names,
    validation_status = validation,
    canonical_manifest_path = manifest_path
  )
}

# ============================================================
# 3. VALIDATE 9/9
# ============================================================

cat("\n===== 3. VALIDATE 9/9 PASS =====\n\n")

final <- rbindlist(all_results, fill = TRUE)

cat("Total tumor samples: ", nrow(final), "\n", sep = "")
cat("All PASS          : ", all(final$validation_status == "PASS"), "\n", sep = "")
cat("Surgery_alone     : ", n_surgery, "\n", sep = "")
cat("Neoadjuvant       : ", n_neoadj, "\n", sep = "")

if (nrow(final) != 9) stop("Expected 9 tumor samples, found ", nrow(final))
if (any(final$validation_status != "PASS")) {
  failed <- final[validation_status != "PASS"]
  stop("Non-PASS samples: ", paste(failed$sample, collapse = ", "))
}
if (n_surgery != 2) stop("Expected 2 Surgery_alone, found ", n_surgery)
if (n_neoadj != 7) stop("Expected 7 Neoadjuvant_treated, found ", n_neoadj)

cat("9/9 PASS ✓\n")

# ============================================================
# 4. GENERATE FINAL TABLE
# ============================================================

cat("\n===== 4. GENERATE FINAL TABLE =====\n\n")

outdir <- "06_tables/cross_dataset/malignant_epithelial/GSE221561_CNV_final"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

fwrite(final, file.path(outdir, "STEP17H9_GSE221561_ALL_9_TUMOR_CNV_CANONICAL.csv"), bom = TRUE)
cat("Saved: ", file.path(outdir, "STEP17H9_GSE221561_ALL_9_TUMOR_CNV_CANONICAL.csv"), "\n", sep = "")

# ============================================================
# 5. DESCRIPTIVE GROUP SUMMARY
# ============================================================

cat("\n===== 5. DESCRIPTIVE GROUP SUMMARY =====\n\n")

cat("EXPLORATORY OBSERVATIONAL COMPARISON ONLY.\n\n")

surg <- final[treatment == "Surgery_alone"]
neoadj <- final[treatment == "Neoadjuvant_treated"]

summarize_group <- function(dt, label) {
  cat(label, " (n = ", nrow(dt), ")\n", sep = "")

  metrics <- c("CNA_signal_ratio", "tumor_median_cnaCor", "correlation_separation",
               "tumor_above_normal95", "tumor_above_normal99", "coherent_CNA_arms")
  labels <- c("CNA signal ratio", "Tumor cnaCor", "Correlation separation",
              "Tumor > normal 95th", "Tumor > normal 99th", "Coherent arms")

  for (j in seq_along(metrics)) {
    vals <- dt[[metrics[j]]]
    med <- median(vals, na.rm = TRUE)
    iqr_vals <- quantile(vals, probs = c(0.25, 0.75), na.rm = TRUE)
    mn <- min(vals, na.rm = TRUE)
    mx <- max(vals, na.rm = TRUE)

    if (metrics[j] %in% c("tumor_above_normal95", "tumor_above_normal99")) {
      cat(sprintf("  %-25s median %.1f%%  IQR %.1f%%–%.1f%%  range %.1f%%–%.1f%%\n",
                  labels[j], 100*med, 100*iqr_vals[1], 100*iqr_vals[2], 100*mn, 100*mx))
    } else if (metrics[j] == "coherent_CNA_arms") {
      cat(sprintf("  %-25s median %.0f   IQR %.0f–%.0f   range %.0f–%.0f\n",
                  labels[j], med, iqr_vals[1], iqr_vals[2], mn, mx))
    } else {
      cat(sprintf("  %-25s median %.3f  IQR %.3f–%.3f  range %.3f–%.3f\n",
                  labels[j], med, iqr_vals[1], iqr_vals[2], mn, mx))
    }
  }
  cat("\n")
}

summarize_group(surg, "Surgery_alone")
summarize_group(neoadj, "Neoadjuvant_treated")

cat("NO p-value. NO t-test. NO Wilcoxon.\n")
cat("NO causal treatment effect interpretation.\n")
cat("EXPLORATORY OBSERVATIONAL COMPARISON ONLY.\n")

# ============================================================
# 6. FINAL SAMPLE-LEVEL TABLE
# ============================================================

cat("\n===== 6. FINAL SAMPLE-LEVEL TABLE =====\n\n")

setorder(final, treatment, sample)

print(final[, .(sample, treatment, tumor_cells, CNA_signal_ratio,
                tumor_median_cnaCor, correlation_separation,
                tumor_above_normal95, tumor_above_normal99,
                coherent_CNA_arms, validation_status)], nrows = 20)

# ============================================================
# FINAL STATUS
# ============================================================

cat("\n============================================================\n")
cat("STEP 17H9 COMPLETE ✓\n")
cat("GSE221561 ALL TUMOR CNV GENERATION COMPLETE\n")
cat("============================================================\n\n")

cat("Tumor samples total          : ", nrow(final), "\n", sep = "")
cat("CNV PASS                     : ", sum(final$validation_status == "PASS"), "\n", sep = "")
cat("CNV failed                   : ", sum(final$validation_status != "PASS"), "\n", sep = "")
cat("\n")
cat("Surgery_alone                : ", n_surgery, "\n", sep = "")
cat("Neoadjuvant_treated          : ", n_neoadj, "\n", sep = "")
cat("\n")
cat("findMalignant run            : NO\n")
cat("Malignant labels assigned    : NO\n")
cat("Cells removed by CNV         : NO\n")
cat("\n")
cat("Interpretation:\n")
cat("EXPLORATORY OBSERVATIONAL COMPARISON ONLY.\n")
cat("NO CAUSAL TREATMENT INTERPRETATION.\n")
cat("\n")
cat("STOP HERE.\n")
cat("DO NOT RUN findMalignant().\n")
cat("DO NOT SET A MALIGNANT THRESHOLD.\n")
cat("DO NOT START CROSS-DATASET INTEGRATION YET.\n")
