# ============================================================
# ESCC Neoadjuvant Project
# STEP 17H10E
#
# MALIGNANT-CALLING STRATEGY AUDIT + S0520T/S2487T TARGETED
#
# Pure audit. NO InferCNA, NO findMalignant, NO cell filtering.
# DO NOT CALL MALIGNANT CELLS. Only compare candidate strategies.
# ============================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
})

cat("\n============================================================\n")
cat("STEP 17H10E: MALIGNANT-CALLING STRATEGY AUDIT\n")
cat("============================================================\n\n")

# ============================================================
# 1. LOAD CELL-LEVEL DATA
# ============================================================

base <- "06_tables/cross_dataset/malignant_epithelial"
audit_dir <- file.path(base, "GSE221561_CNV_strategy_audit")

sample_info <- list(
  S6417T  = list(dir = "GSE221561_S6417T_pilot",  res = "17G3_results", file = "GSE221561_S6417T_17G3_cell_CNA_metrics.csv"),
  S9478T  = list(dir = "GSE221561_S9478T_pilot",  res = "17H5C_results", file = "GSE221561_S9478T_17H5C_cell_CNA_metrics.csv"),
  S1265T  = list(dir = "GSE221561_S1265T_pilot",  res = "17H3_results", file = "GSE221561_S1265T_17H3_cell_CNA_metrics.csv"),
  S2423T  = list(dir = "GSE221561_S2423T_pilot",  res = "17H7C_results", file = "GSE221561_S2423T_17H7C_cell_CNA_metrics.csv"),
  S1535T  = list(dir = "GSE221561_S1535T_pilot",  res = "17H9_results", file = "GSE221561_S1535T_17H9_cell_CNA_metrics.csv"),
  S1315TS = list(dir = "GSE221561_S1315TS_pilot", res = "17H9_results", file = "GSE221561_S1315TS_17H9_cell_CNA_metrics.csv"),
  S6829T  = list(dir = "GSE221561_S6829T_pilot",  res = "17H9_results", file = "GSE221561_S6829T_17H9_cell_CNA_metrics.csv"),
  S2487T  = list(dir = "GSE221561_S2487T_pilot",  res = "17H9_results", file = "GSE221561_S2487T_17H9_cell_CNA_metrics.csv"),
  S0520T  = list(dir = "GSE221561_S0520T_pilot",  res = "17H9_results", file = "GSE221561_S0520T_17H9_cell_CNA_metrics.csv")
)

landscape <- fread(file.path(audit_dir, "STEP17H10_sample_level_CNA_landscape.csv"))

# ============================================================
# 2. FOUR CANDIDATE STRATEGIES
# ============================================================

cat("===== 2. FOUR CANDIDATE STRATEGIES =====\n\n")

strategy_results <- list()

for (nm in names(sample_info)) {
  si <- sample_info[[nm]]
  cm <- fread(file.path(base, si$dir, si$res, si$file))

  tumor <- cm[role == "TUMOR_TEST"]
  normal <- cm[role == "NORMAL_REFERENCE"]

  n_sig_95 <- as.numeric(quantile(normal$cna_signal, 0.95, na.rm = TRUE))
  n_sig_99 <- as.numeric(quantile(normal$cna_signal, 0.99, na.rm = TRUE))
  n_cor_95 <- as.numeric(quantile(normal$cna_cor_to_tumor, 0.95, na.rm = TRUE))

  # Strategy A: signal > normal Q95
  frac_A <- mean(tumor$cna_signal > n_sig_95, na.rm = TRUE)

  # Strategy B: cnaCor > normal-to-tumor Q95
  frac_B <- mean(tumor$cna_cor_to_tumor > n_cor_95, na.rm = TRUE)

  # Strategy C: signal > normal Q95 AND cnaCor > normal Q95
  frac_C <- mean(tumor$cna_signal > n_sig_95 & tumor$cna_cor_to_tumor > n_cor_95, na.rm = TRUE)

  # Strategy D: signal > normal Q99 AND cnaCor > normal Q95
  frac_D <- mean(tumor$cna_signal > n_sig_99 & tumor$cna_cor_to_tumor > n_cor_95, na.rm = TRUE)

  strategy_results[[nm]] <- data.table(
    sample = nm,
    treatment = landscape[sample == nm, treatment],
    tumor_cells = nrow(tumor),
    strategy_A_signal_Q95 = frac_A,
    strategy_B_corr_Q95 = frac_B,
    strategy_C_dual_Q95 = frac_C,
    strategy_D_signal_Q99_corr_Q95 = frac_D,
    normal_signal_Q95 = n_sig_95,
    normal_signal_Q99 = n_sig_99,
    normal_corr_Q95 = n_cor_95
  )

  cat(sprintf("  %-10s  A=%.1f%%  B=%.1f%%  C=%.1f%%  D=%.1f%%\n",
              nm, 100*frac_A, 100*frac_B, 100*frac_C, 100*frac_D))
}

strat_dt <- rbindlist(strategy_results)

fwrite(strat_dt, file.path(audit_dir, "STEP17H10_candidate_strategy_comparison.csv"), bom = TRUE)
cat("\nSaved: STEP17H10_candidate_strategy_comparison.csv\n")

# ============================================================
# 3. STRATEGY STABILITY ASSESSMENT
# ============================================================

cat("\n===== 3. STRATEGY STABILITY ASSESSMENT =====\n\n")

strat_metrics <- data.table(
  strategy = c("A_signal_Q95", "B_corr_Q95", "C_dual_Q95", "D_signal_Q99_corr_Q95"),
  median_frac = c(median(strat_dt$strategy_A_signal_Q95),
                  median(strat_dt$strategy_B_corr_Q95),
                  median(strat_dt$strategy_C_dual_Q95),
                  median(strat_dt$strategy_D_signal_Q99_corr_Q95)),
  iqr_frac = c(IQR(strat_dt$strategy_A_signal_Q95),
               IQR(strat_dt$strategy_B_corr_Q95),
               IQR(strat_dt$strategy_C_dual_Q95),
               IQR(strat_dt$strategy_D_signal_Q99_corr_Q95)),
  range_frac = c(diff(range(strat_dt$strategy_A_signal_Q95)),
                 diff(range(strat_dt$strategy_B_corr_Q95)),
                 diff(range(strat_dt$strategy_C_dual_Q95)),
                 diff(range(strat_dt$strategy_D_signal_Q99_corr_Q95)))
)

# Check if any strategy has near-0 in weak samples and reasonable range in strong
# A stable strategy should have: low in weak, moderate-high in strong, not saturating
strat_metrics[, dynamic_range := range_frac]

cat("Strategy stability:\n\n")
print(strat_metrics)

# Most stable: smallest IQR that still has meaningful dynamic range
# But also check if CNA-weak samples have near-0
weak_samples <- strat_dt[treatment == "Neoadjuvant_treated" & tumor_cells > 1000]
strong_samples <- strat_dt[treatment == "Surgery_alone"]

strat_metrics[, weak_median := c(
  median(weak_samples$strategy_A_signal_Q95),
  median(weak_samples$strategy_B_corr_Q95),
  median(weak_samples$strategy_C_dual_Q95),
  median(weak_samples$strategy_D_signal_Q99_corr_Q95)
)]

strat_metrics[, strong_median := c(
  median(strong_samples$strategy_A_signal_Q95),
  median(strong_samples$strategy_B_corr_Q95),
  median(strong_samples$strategy_C_dual_Q95),
  median(strong_samples$strategy_D_signal_Q99_corr_Q95)
)]

cat("\nWeak vs strong group medians:\n")
print(strat_metrics[, .(strategy, weak_median, strong_median, iqr_frac)])

# Determine most stable
# Criteria: reasonable separation between weak and strong, low IQR
strat_metrics[, stability_score := (strong_median - weak_median) / (iqr_frac + 0.01)]
best <- strat_metrics[which.max(stability_score)]

most_stable <- best$strategy[1]
cat("\nMOST_STABLE_CANDIDATE_STRATEGY: ", most_stable, "\n", sep = "")
cat("Rationale: highest strong-weak separation relative to IQR.\n")

# ============================================================
# 4. S0520T TARGETED AUDIT
# ============================================================

cat("\n===== 4. S0520T TARGETED AUDIT =====\n\n")

s0520_cm <- fread(file.path(base, "GSE221561_S0520T_pilot", "17H9_results", "GSE221561_S0520T_17H9_cell_CNA_metrics.csv"))
s0520_tumor <- s0520_cm[role == "TUMOR_TEST"]
s0520_normal <- s0520_cm[role == "NORMAL_REFERENCE"]

s0520_audit <- data.table(
  metric = c("tumor_cells", "normal_cells",
             "tumor_signal_median", "tumor_signal_Q25", "tumor_signal_Q75", "tumor_signal_Q90", "tumor_signal_Q95", "tumor_signal_Q99",
             "normal_signal_median", "normal_signal_Q95", "normal_signal_Q99",
             "tumor_cnaCor_median", "tumor_cnaCor_Q25", "tumor_cnaCor_Q75", "tumor_cnaCor_Q90", "tumor_cnaCor_Q95",
             "normal_corr_median", "normal_corr_Q95",
             "frac_signal_gt_normalQ95", "frac_signal_gt_normalQ99",
             "frac_corr_gt_normalQ95",
             "dual_metric_candidate_fraction"),
  value = c(nrow(s0520_tumor), nrow(s0520_normal),
            median(s0520_tumor$cna_signal), quantile(s0520_tumor$cna_signal, 0.25), quantile(s0520_tumor$cna_signal, 0.75),
            quantile(s0520_tumor$cna_signal, 0.90), quantile(s0520_tumor$cna_signal, 0.95), quantile(s0520_tumor$cna_signal, 0.99),
            median(s0520_normal$cna_signal), quantile(s0520_normal$cna_signal, 0.95), quantile(s0520_normal$cna_signal, 0.99),
            median(s0520_tumor$cna_cor_to_tumor), quantile(s0520_tumor$cna_cor_to_tumor, 0.25),
            quantile(s0520_tumor$cna_cor_to_tumor, 0.75), quantile(s0520_tumor$cna_cor_to_tumor, 0.90),
            quantile(s0520_tumor$cna_cor_to_tumor, 0.95),
            median(s0520_normal$cna_cor_to_tumor), quantile(s0520_normal$cna_cor_to_tumor, 0.95),
            mean(s0520_tumor$cna_signal > quantile(s0520_normal$cna_signal, 0.95)),
            mean(s0520_tumor$cna_signal > quantile(s0520_normal$cna_signal, 0.99)),
            mean(s0520_tumor$cna_cor_to_tumor > quantile(s0520_normal$cna_cor_to_tumor, 0.95)),
            mean(s0520_tumor$cna_signal > quantile(s0520_normal$cna_signal, 0.95) &
                 s0520_tumor$cna_cor_to_tumor > quantile(s0520_normal$cna_cor_to_tumor, 0.95))
  )
)

fwrite(s0520_audit, file.path(audit_dir, "STEP17H10_S0520T_targeted_audit.csv"), bom = TRUE)

# Interpretation
s0520_dual <- s0520_audit[metric == "dual_metric_candidate_fraction", value]
s0520_tumor_q25 <- s0520_audit[metric == "tumor_signal_Q25", value]
s0520_tumor_median <- s0520_audit[metric == "tumor_signal_median", value]

cat("S0520T audit results:\n")
print(s0520_audit, nrows = 25)

# Check if broad cellular signal or subpopulation-driven
# If Q25 is already well above normal Q95, it's broad
n_q95 <- s0520_audit[metric == "normal_signal_Q95", value]

if (s0520_tumor_q25 > n_q95) {
  s0520_interp <- "BROAD_CELLULAR_SIGNAL"
} else if (s0520_dual > 0.5) {
  s0520_interp <- "BROAD_CELLULAR_SIGNAL"
} else {
  s0520_interp <- "SUBPOPULATION_DRIVEN"
}

cat("\nS0520T_STRONG_SIGNAL_INTERPRETATION: ", s0520_interp, "\n", sep = "")

# ============================================================
# 5. S2487T DISCORDANCE AUDIT
# ============================================================

cat("\n===== 5. S2487T DISCORDANCE AUDIT =====\n\n")

s2487_cm <- fread(file.path(base, "GSE221561_S2487T_pilot", "17H9_results", "GSE221561_S2487T_17H9_cell_CNA_metrics.csv"))
s2487_tumor <- s2487_cm[role == "TUMOR_TEST"]
s2487_normal <- s2487_cm[role == "NORMAL_REFERENCE"]

n_sig_95 <- quantile(s2487_normal$cna_signal, 0.95, na.rm = TRUE)
n_cor_95 <- quantile(s2487_normal$cna_cor_to_tumor, 0.95, na.rm = TRUE)

s2487_discordance <- data.table(
  metric = c("tumor_cells", "normal_cells",
             "tumor_signal_median", "tumor_signal_Q25", "tumor_signal_Q75",
             "normal_signal_Q95",
             "tumor_cnaCor_median", "tumor_cnaCor_Q25", "tumor_cnaCor_Q75",
             "normal_corr_Q95",
             "frac_signal_gt_normalQ95", "frac_corr_gt_normalQ95",
             "dual_metric_candidate_fraction",
             "signal_ratio_from_landscape"),
  value = c(nrow(s2487_tumor), nrow(s2487_normal),
            median(s2487_tumor$cna_signal), quantile(s2487_tumor$cna_signal, 0.25),
            quantile(s2487_tumor$cna_signal, 0.75), n_sig_95,
            median(s2487_tumor$cna_cor_to_tumor), quantile(s2487_tumor$cna_cor_to_tumor, 0.25),
            quantile(s2487_tumor$cna_cor_to_tumor, 0.75), n_cor_95,
            mean(s2487_tumor$cna_signal > n_sig_95),
            mean(s2487_tumor$cna_cor_to_tumor > n_cor_95),
            mean(s2487_tumor$cna_signal > n_sig_95 & s2487_tumor$cna_cor_to_tumor > n_cor_95),
            landscape[sample == "S2487T", CNA_signal_ratio]
  )
)

fwrite(s2487_discordance, file.path(audit_dir, "STEP17H10_S2487T_discordance_audit.csv"), bom = TRUE)

cat("S2487T discordance audit:\n")
print(s2487_discordance, nrows = 20)

# Interpretation
s2487_sig_frac <- s2487_discordance[metric == "frac_signal_gt_normalQ95", value]
s2487_cor_frac <- s2487_discordance[metric == "frac_corr_gt_normalQ95", value]

if (s2487_sig_frac > 0.5 & s2487_cor_frac < 0.3) {
  s2487_pattern <- "HIGH_SIGNAL_LOW_COHERENCE"
} else if (s2487_sig_frac > 0.3 & s2487_cor_frac > 0.3) {
  s2487_pattern <- "SUBPOPULATION_HETEROGENEITY"
} else {
  s2487_pattern <- "UNRESOLVED"
}

cat("\nS2487T_DISCORDANCE_PATTERN: ", s2487_pattern, "\n", sep = "")

# ============================================================
# 6. TREATMENT GROUP SUMMARY
# ============================================================

cat("\n===== 6. TREATMENT GROUP SUMMARY =====\n\n")

cat("EXPLORATORY OBSERVATIONAL DESCRIPTION ONLY.\n\n")

surg <- landscape[treatment == "Surgery_alone"]
neoadj <- landscape[treatment == "Neoadjuvant_treated"]

summarize <- function(dt, label) {
  cat(label, " (n = ", nrow(dt), ")\n", sep = "")
  for (m in c("CNA_signal_ratio", "tumor_median_cnaCor", "correlation_separation",
              "tumor_above_normal95", "coherent_CNA_arms")) {
    vals <- dt[[m]]
    med <- median(vals, na.rm = TRUE)
    q <- quantile(vals, probs = c(0.25, 0.75), na.rm = TRUE)
    mn <- min(vals, na.rm = TRUE)
    mx <- max(vals, na.rm = TRUE)
    cat(sprintf("  %-25s median %.3f  IQR %.3f–%.3f  range %.3f–%.3f\n",
                m, med, q[1], q[2], mn, mx))
  }
  cat("\n")
}

summarize(surg, "Surgery_alone")
summarize(neoadj, "Neoadjuvant_treated")

cat("Substantial within-group heterogeneity is present.\n")
cat("NO causal treatment interpretation.\n")

# ============================================================
# 7. SAVE STRATEGY RECOMMENDATION
# ============================================================

rec_text <- paste0(
  "STEP 17H10 STRATEGY RECOMMENDATION\n",
  "===================================\n\n",
  "Most stable candidate strategy: ", most_stable, "\n",
  "Rationale: highest strong-weak separation relative to IQR across 9 samples.\n\n",
  "S0520T interpretation: ", s0520_interp, "\n",
  "S2487T discordance: ", s2487_pattern, "\n\n",
  "IMPORTANT: No malignant labels have been assigned.\n",
  "All strategies remain CANDIDATE FRACTIONS only.\n"
)

writeLines(rec_text, file.path(audit_dir, "STEP17H10_strategy_recommendation.txt"))
cat("\nSaved: STEP17H10_strategy_recommendation.txt\n")

cat("\nSTEP 17H10E COMPLETE\n")
