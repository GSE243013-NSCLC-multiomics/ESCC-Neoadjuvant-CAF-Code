# ============================================================
# ESCC Neoadjuvant Project
# STEP 17H10C
#
# CELL-LEVEL SEPARABILITY AUDIT + DISTRIBUTION FIGURES
#
# Pure audit. NO InferCNA, NO findMalignant, NO cell filtering.
# ============================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
})

cat("\n============================================================\n")
cat("STEP 17H10C: CELL-LEVEL SEPARABILITY AUDIT\n")
cat("============================================================\n\n")

# ============================================================
# 1. LOAD ALL CELL-LEVEL METRICS
# ============================================================

cat("===== 1. LOAD CELL-LEVEL METRICS =====\n\n")

base <- "06_tables/cross_dataset/malignant_epithelial"

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

all_cell_data <- list()
separability <- list()

for (nm in names(sample_info)) {
  si <- sample_info[[nm]]
  fpath <- file.path(base, si$dir, si$res, si$file)
  cm <- fread(fpath)
  cm[, sample := nm]
  all_cell_data[[nm]] <- cm

  tumor <- cm[role == "TUMOR_TEST"]
  normal <- cm[role == "NORMAL_REFERENCE"]

  # Quantiles
  t_sig_q <- quantile(tumor$cna_signal, probs = c(0.25, 0.5, 0.75, 0.90, 0.95, 0.99), na.rm = TRUE)
  n_sig_q <- quantile(normal$cna_signal, probs = c(0.25, 0.5, 0.75, 0.90, 0.95, 0.99), na.rm = TRUE)
  t_cor_q <- quantile(tumor$cna_cor_to_tumor, probs = c(0.25, 0.5, 0.75, 0.90, 0.95), na.rm = TRUE)
  n_cor_q <- quantile(normal$cna_cor_to_tumor, probs = c(0.25, 0.5, 0.75, 0.90, 0.95), na.rm = TRUE)

  # Separability fractions
  n_sig_95 <- n_sig_q["95%"]
  n_sig_99 <- n_sig_q["99%"]
  n_cor_95 <- n_cor_q["95%"]

  frac_sig_gt_n95 <- mean(tumor$cna_signal > n_sig_95, na.rm = TRUE)
  frac_sig_gt_n99 <- mean(tumor$cna_signal > n_sig_99, na.rm = TRUE)
  frac_cor_gt_n95 <- mean(tumor$cna_cor_to_tumor > n_cor_95, na.rm = TRUE)
  frac_dual <- mean(tumor$cna_signal > n_sig_95 & tumor$cna_cor_to_tumor > n_cor_95, na.rm = TRUE)

  separability[[nm]] <- data.table(
    sample = nm,
    tumor_cells = nrow(tumor),
    normal_cells = nrow(normal),
    tumor_signal_median = t_sig_q["50%"],
    tumor_signal_Q25 = t_sig_q["25%"],
    tumor_signal_Q75 = t_sig_q["75%"],
    tumor_signal_Q90 = t_sig_q["90%"],
    tumor_signal_Q95 = t_sig_q["95%"],
    tumor_signal_Q99 = t_sig_q["99%"],
    normal_signal_median = n_sig_q["50%"],
    normal_signal_Q95 = n_sig_95,
    normal_signal_Q99 = n_sig_99,
    tumor_cnaCor_median = t_cor_q["50%"],
    tumor_cnaCor_Q25 = t_cor_q["25%"],
    tumor_cnaCor_Q75 = t_cor_q["75%"],
    tumor_cnaCor_Q90 = t_cor_q["90%"],
    tumor_cnaCor_Q95 = t_cor_q["95%"],
    normal_corr_median = n_cor_q["50%"],
    normal_corr_Q95 = n_cor_95,
    frac_tumor_signal_gt_normal_Q95 = frac_sig_gt_n95,
    frac_tumor_signal_gt_normal_Q99 = frac_sig_gt_n99,
    frac_tumor_corr_gt_normal_Q95 = frac_cor_gt_n95,
    dual_metric_candidate_fraction = frac_dual
  )

  cat(sprintf("  %-10s  tumor=%d  normal=%d  dual=%.1f%%\n",
              nm, nrow(tumor), nrow(normal), 100 * frac_dual))
}

sep_dt <- rbindlist(separability)

fwrite(sep_dt, file.path("06_tables/cross_dataset/malignant_epithelial/GSE221561_CNV_strategy_audit",
                          "STEP17H10_cell_level_CNA_separability.csv"), bom = TRUE)
cat("\nSaved: STEP17H10_cell_level_CNA_separability.csv\n")

# ============================================================
# 2. FIGURES
# ============================================================

cat("\n===== 2. GENERATE FIGURES =====\n\n")

figdir <- "05_figures/GSE221561/malignant_epithelial/strategy_audit"
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

# Read canonical landscape for treatment info
landscape <- fread("06_tables/cross_dataset/malignant_epithelial/GSE221561_CNV_strategy_audit/STEP17H10_sample_level_CNA_landscape.csv")

# Combine all cell data
all_cells <- rbindlist(all_cell_data, use.names = TRUE)

# Add treatment info
all_cells <- merge(all_cells, landscape[, .(sample, treatment)], by = "sample", all.x = TRUE)

# --- Figure 1: CNA signal ratio by sample ---
png(file.path(figdir, "STEP17H10_all9_CNA_signal_ratio_by_sample.png"), width = 3600, height = 2400, res = 300)
par(mar = c(7, 5, 4, 2))
boxplot(CNA_signal_ratio ~ sample, data = landscape,
        col = ifelse(landscape$treatment == "Surgery_alone", "steelblue", "coral"),
        main = "CNA Signal Ratio by Sample",
        ylab = "CNA Signal Ratio (Tumor / Normal median)",
        xlab = "", las = 2, cex.axis = 0.8)
legend("topright", legend = c("Surgery_alone", "Neoadjuvant_treated"),
       fill = c("steelblue", "coral"), cex = 0.8, bty = "n")
dev.off()

# --- Figure 2: cnaCor by sample ---
tumor_cells <- all_cells[role == "TUMOR_TEST"]
png(file.path(figdir, "STEP17H10_all9_cnaCor_by_sample.png"), width = 3600, height = 2400, res = 300)
par(mar = c(7, 5, 4, 2))
boxplot(cna_cor_to_tumor ~ sample, data = tumor_cells,
        col = ifelse(landscape$treatment[match(unique(tumor_cells$sample), landscape$sample)] == "Surgery_alone", "steelblue", "coral"),
        main = "Tumor cnaCor by Sample",
        ylab = "cnaCor (correlation to tumor mean profile)",
        xlab = "", las = 2, cex.axis = 0.8)
dev.off()

# --- Figure 3: Coherent arms by sample ---
png(file.path(figdir, "STEP17H10_all9_coherent_arms_by_sample.png"), width = 3600, height = 2400, res = 300)
par(mar = c(7, 5, 4, 2))
barplot(landscape$coherent_CNA_arms, names.arg = landscape$sample,
        col = ifelse(landscape$treatment == "Surgery_alone", "steelblue", "coral"),
        main = "Coherent CNA Arms by Sample",
        ylab = "Number of Coherent Arms", las = 2, cex.axis = 0.8)
dev.off()

# --- Figure 4: Tumor cells vs CNA signal ratio ---
png(file.path(figdir, "STEP17H10_tumor_cells_vs_CNA_signal_ratio.png"), width = 3200, height = 2800, res = 300)
par(mar = c(5, 5, 4, 2))
plot(landscape$tumor_cells, landscape$CNA_signal_ratio,
     col = ifelse(landscape$treatment == "Surgery_alone", "steelblue", "coral"),
     pch = 19, cex = 1.5,
     xlab = "Tumor Cell Count",
     ylab = "CNA Signal Ratio",
     main = "Tumor Cell Count vs CNA Signal Ratio")
text(landscape$tumor_cells, landscape$CNA_signal_ratio, labels = landscape$sample,
     pos = 3, cex = 0.7, offset = 0.3)
abline(lm(CNA_signal_ratio ~ tumor_cells, data = landscape), lty = 2, col = "grey50")
dev.off()

# --- Figure 5: Signal ratio vs cnaCor ---
png(file.path(figdir, "STEP17H10_signal_ratio_vs_cnaCor.png"), width = 3200, height = 2800, res = 300)
par(mar = c(5, 5, 4, 2))
plot(landscape$CNA_signal_ratio, landscape$tumor_median_cnaCor,
     col = ifelse(landscape$treatment == "Surgery_alone", "steelblue", "coral"),
     pch = 19, cex = 1.5,
     xlab = "CNA Signal Ratio",
     ylab = "Tumor median cnaCor",
     main = "CNA Signal Ratio vs cnaCor")
text(landscape$CNA_signal_ratio, landscape$tumor_median_cnaCor, labels = landscape$sample,
     pos = 3, cex = 0.7, offset = 0.3)
dev.off()

# --- Figure 6: cnaCor vs coherent arms ---
png(file.path(figdir, "STEP17H10_cnaCor_vs_coherent_arms.png"), width = 3200, height = 2800, res = 300)
par(mar = c(5, 5, 4, 2))
plot(landscape$tumor_median_cnaCor, landscape$coherent_CNA_arms,
     col = ifelse(landscape$treatment == "Surgery_alone", "steelblue", "coral"),
     pch = 19, cex = 1.5,
     xlab = "Tumor median cnaCor",
     ylab = "Coherent CNA Arms",
     main = "cnaCor vs Coherent Arms")
text(landscape$tumor_median_cnaCor, landscape$coherent_CNA_arms, labels = landscape$sample,
     pos = 3, cex = 0.7, offset = 0.3)
dev.off()

# --- Figure 7: Metric rank heatmap ---
rank_matrix <- fread("06_tables/cross_dataset/malignant_epithelial/GSE221561_CNV_strategy_audit/STEP17H10_CNA_metric_rank_matrix.csv")

rank_mat <- as.matrix(rank_matrix[, .(signal_ratio, cnaCor, corr_separation, above95, coherent_arms)])
rownames(rank_mat) <- rank_matrix$sample

png(file.path(figdir, "STEP17H10_metric_rank_heatmap.png"), width = 3200, height = 2400, res = 300)
par(mar = c(5, 8, 4, 2))
image(t(rank_mat), col = hcl.colors(9, "YlOrRd", rev = TRUE), axes = FALSE,
      main = "Metric Rank Heatmap (1=strongest)")
axis(2, at = seq(0, 1, length.out = nrow(rank_mat)), labels = rownames(rank_mat), las = 2, cex.axis = 0.8)
axis(1, at = seq(0, 1, length.out = ncol(rank_mat)), labels = colnames(rank_mat), las = 2, cex.axis = 0.8)
for (i in seq_len(nrow(rank_mat))) {
  for (j in seq_len(ncol(rank_mat))) {
    text((j - 0.5) / ncol(rank_mat), (i - 0.5) / nrow(rank_mat),
         labels = rank_mat[i, j], cex = 0.9, font = 2)
  }
}
dev.off()

# --- Figure 8: Dual metric fraction ---
png(file.path(figdir, "STEP17H10_cell_level_dual_metric_fraction.png"), width = 3600, height = 2400, res = 300)
par(mar = c(7, 5, 4, 2))
barplot(sep_dt$dual_metric_candidate_fraction * 100, names.arg = sep_dt$sample,
        col = ifelse(landscape$treatment[match(sep_dt$sample, landscape$sample)] == "Surgery_alone", "steelblue", "coral"),
        main = "Dual-Metric Candidate Fraction by Sample",
        ylab = "Candidate Fraction (%)", las = 2, cex.axis = 0.8,
        ylim = c(0, 100))
dev.off()

cat("8 figures saved to: ", figdir, "\n", sep = "")

cat("\nSTEP 17H10C COMPLETE\n")
