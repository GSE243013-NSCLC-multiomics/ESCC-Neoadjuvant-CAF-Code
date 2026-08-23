# ============================================================
# ESCC Neoadjuvant Project
# STEP 17H10A
#
# GSE221561 ALL-9-SAMPLE CNA LANDSCAPE
#
# Pure audit. NO InferCNA, NO findMalignant, NO cell filtering.
# ============================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
})

cat("\n============================================================\n")
cat("STEP 17H10A: GSE221561 CNA SAMPLE-LEVEL LANDSCAPE\n")
cat("============================================================\n\n")

# ============================================================
# 1. LOAD CANONICAL TABLE
# ============================================================

canonical_file <- "06_tables/cross_dataset/malignant_epithelial/GSE221561_CNV_final/STEP17H9_GSE221561_ALL_9_TUMOR_CNV_CANONICAL.csv"

if (!file.exists(canonical_file)) stop("Missing: ", canonical_file)

dt <- fread(canonical_file)

cat("Samples loaded: ", nrow(dt), "\n", sep = "")
if (nrow(dt) != 9) stop("Expected 9 samples, found ", nrow(dt))
if (any(dt$validation_status != "PASS")) stop("Non-PASS samples found")
cat("9/9 PASS ✓\n")

# ============================================================
# 2. RECOMPUTE MISSING METRICS FROM CELL-LEVEL DATA
# ============================================================

cat("\n===== 2. RECOMPUTE METRICS FROM CELL-LEVEL DATA =====\n\n")

# For each sample, read its canonical cell CNA metrics and recompute
base <- "06_tables/cross_dataset/malignant_epithelial"

sample_dirs <- list(
  S6417T  = list(dir = "GSE221561_S6417T_pilot",  results = "17G3_results", file = "GSE221561_S6417T_17G3_cell_CNA_metrics.csv"),
  S9478T  = list(dir = "GSE221561_S9478T_pilot",  results = "17H5C_results", file = "GSE221561_S9478T_17H5C_cell_CNA_metrics.csv"),
  S1265T  = list(dir = "GSE221561_S1265T_pilot",  results = "17H3_results", file = "GSE221561_S1265T_17H3_cell_CNA_metrics.csv"),
  S2423T  = list(dir = "GSE221561_S2423T_pilot",  results = "17H7C_results", file = "GSE221561_S2423T_17H7C_cell_CNA_metrics.csv"),
  S1535T  = list(dir = "GSE221561_S1535T_pilot",  results = "17H9_results", file = "GSE221561_S1535T_17H9_cell_CNA_metrics.csv"),
  S1315TS = list(dir = "GSE221561_S1315TS_pilot", results = "17H9_results", file = "GSE221561_S1315TS_17H9_cell_CNA_metrics.csv"),
  S6829T  = list(dir = "GSE221561_S6829T_pilot",  results = "17H9_results", file = "GSE221561_S6829T_17H9_cell_CNA_metrics.csv"),
  S2487T  = list(dir = "GSE221561_S2487T_pilot",  results = "17H9_results", file = "GSE221561_S2487T_17H9_cell_CNA_metrics.csv"),
  S0520T  = list(dir = "GSE221561_S0520T_pilot",  results = "17H9_results", file = "GSE221561_S0520T_17H9_cell_CNA_metrics.csv")
)

recomputed <- list()

for (nm in names(sample_dirs)) {
  sd <- sample_dirs[[nm]]
  metrics_path <- file.path(base, sd$dir, sd$results, sd$file)
  if (!file.exists(metrics_path)) {
    stop("Cell metrics missing for ", nm, ": ", metrics_path)
  }
  cm <- fread(metrics_path)
  tumor_m <- cm[role == "TUMOR_TEST"]
  normal_m <- cm[role == "NORMAL_REFERENCE"]

  normal_signal_95 <- as.numeric(quantile(normal_m$cna_signal, 0.95, na.rm = TRUE))

  recomputed[[nm]] <- data.table(
    sample = nm,
    tumor_cells = nrow(tumor_m),
    normal_cells = nrow(normal_m),
    tumor_median_CNA_signal = median(tumor_m$cna_signal, na.rm = TRUE),
    normal_median_CNA_signal = median(normal_m$cna_signal, na.rm = TRUE),
    tumor_median_cnaCor = median(tumor_m$cna_cor_to_tumor, na.rm = TRUE),
    normal_median_corr_to_tumor = median(normal_m$cna_cor_to_tumor, na.rm = TRUE),
    tumor_above_normal95 = mean(tumor_m$cna_signal > normal_signal_95, na.rm = TRUE)
  )
}

rec <- rbindlist(recomputed)

# Merge with canonical
dt <- merge(dt, rec[, .(sample, tumor_median_CNA_signal_rc = tumor_median_CNA_signal,
                         normal_median_CNA_signal_rc = normal_median_CNA_signal,
                         tumor_median_cnaCor_rc = tumor_median_cnaCor,
                         normal_median_corr_to_tumor_rc = normal_median_corr_to_tumor,
                         tumor_above_normal95_rc = tumor_above_normal95)],
            by = "sample", all.x = TRUE)

# Use recomputed values where canonical has NA
dt[, tumor_above_normal95 := fifelse(is.na(tumor_above_normal95), tumor_above_normal95_rc, tumor_above_normal95)]
dt[, tumor_above_normal99 := fifelse(is.na(tumor_above_normal99), NA_real_, tumor_above_normal99)]

# ============================================================
# 3. RANKS
# ============================================================

cat("\n===== 3. COMPUTE RANKS =====\n\n")

dt[, rank_signal_ratio := rank(-CNA_signal_ratio)]
dt[, rank_cnaCor := rank(-tumor_median_cnaCor)]
dt[, rank_corr_separation := rank(-correlation_separation)]
dt[, rank_above95 := rank(-tumor_above_normal95)]
dt[, rank_coherent_arms := rank(-coherent_CNA_arms)]

cat("Ranks computed.\n")

# ============================================================
# 4. SAVE LANDSCAPE TABLE
# ============================================================

outdir <- "06_tables/cross_dataset/malignant_epithelial/GSE221561_CNV_strategy_audit"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

landscape <- dt[, .(sample, treatment, tumor_cells,
                     CNA_signal_ratio, tumor_median_CNA_signal, normal_median_CNA_signal,
                     tumor_median_cnaCor, normal_median_corr_to_tumor, correlation_separation,
                     tumor_above_normal95, tumor_above_normal99, coherent_CNA_arms,
                     rank_signal_ratio, rank_cnaCor, rank_corr_separation,
                     rank_above95, rank_coherent_arms)]

fwrite(landscape, file.path(outdir, "STEP17H10_sample_level_CNA_landscape.csv"), bom = TRUE)
cat("Saved: STEP17H10_sample_level_CNA_landscape.csv\n")

# ============================================================
# 5. CELL-COUNT SENSITIVITY ANALYSIS
# ============================================================

cat("\n===== 5. CELL-COUNT SENSITIVITY ANALYSIS =====\n\n")

# Spearman rho (no p-value)
rho_ratio <- cor(dt$tumor_cells, dt$CNA_signal_ratio, method = "spearman", use = "complete.obs")
rho_cnaCor <- cor(dt$tumor_cells, dt$tumor_median_cnaCor, method = "spearman", use = "complete.obs")
rho_arms <- cor(dt$tumor_cells, dt$coherent_CNA_arms, method = "spearman", use = "complete.obs")
rho_above95 <- cor(dt$tumor_cells, dt$tumor_above_normal95, method = "spearman", use = "complete.obs")

cat(sprintf("rho(tumor_cells, CNA_signal_ratio)   = %.3f\n", rho_ratio))
cat(sprintf("rho(tumor_cells, cnaCor)             = %.3f\n", rho_cnaCor))
cat(sprintf("rho(tumor_cells, coherent_arms)      = %.3f\n", rho_arms))
cat(sprintf("rho(tumor_cells, above95)            = %.3f\n", rho_above95))
cat("\nDESCRIPTIVE ONLY — n=9.\n")

sensitivity <- data.table(
  comparison = c("tumor_cells vs CNA_signal_ratio",
                 "tumor_cells vs tumor_median_cnaCor",
                 "tumor_cells vs coherent_CNA_arms",
                 "tumor_cells vs tumor_above_normal95"),
  spearman_rho = c(rho_ratio, rho_cnaCor, rho_arms, rho_above95),
  note = rep("DESCRIPTIVE ONLY — n=9", 4)
)

fwrite(sensitivity, file.path(outdir, "STEP17H10_tumor_cellcount_sensitivity.csv"), bom = TRUE)
cat("Saved: STEP17H10_tumor_cellcount_sensitivity.csv\n")

# ============================================================
# 6. PRINT LANDSCAPE
# ============================================================

cat("\n===== 6. SAMPLE-LEVEL LANDSCAPE =====\n\n")

print(landscape[order(rank_signal_ratio), .(
  sample, treatment, tumor_cells, CNA_signal_ratio, tumor_median_cnaCor,
  correlation_separation, tumor_above_normal95, coherent_CNA_arms,
  rank_signal_ratio, rank_cnaCor, rank_coherent_arms
)], nrows = 20)

cat("\nSTEP 17H10A COMPLETE\n")
