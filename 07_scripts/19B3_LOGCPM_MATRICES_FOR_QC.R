#!/usr/bin/env Rscript
# =============================================================================
# STEP 19B3: LOGCPM MATRICES FOR QC ONLY
# =============================================================================

suppressPackageStartupMessages({
  library(edgeR)
})

base_dir <- "~/Downloads/ESCC Neoadjuvant"
qc_dir <- file.path(base_dir, "06_tables/cross_dataset/epithelial_pseudobulk/QC")

cat("Loading DGEList objects...\n")

dge_221561 <- readRDS(file.path(qc_dir, "GSE221561_DGEList_TMM.rds"))
dge_197677 <- readRDS(file.path(qc_dir, "GSE197677_DGEList_TMM.rds"))
dge_160269 <- readRDS(file.path(qc_dir, "GSE160269_DGEList_TMM.rds"))

cat("Loaded.\n\n")

# ---- Calculate logCPM ----

cat("=== CALCULATING logCPM (TMM-normalized, QC ONLY) ===\n\n")

logcpm_221561 <- cpm(dge_221561, log = TRUE, prior.count = 2)
cat("GSE221561 logCPM: ", paste(dim(logcpm_221561), collapse = " x "), "\n")

logcpm_197677 <- cpm(dge_197677, log = TRUE, prior.count = 2)
cat("GSE197677 logCPM: ", paste(dim(logcpm_197677), collapse = " x "), "\n")

logcpm_160269 <- cpm(dge_160269, log = TRUE, prior.count = 2)
cat("GSE160269 logCPM: ", paste(dim(logcpm_160269), collapse = " x "), "\n")

cat("\n")

# ---- Save logCPM matrices ----

cat("=== SAVING LOGCPM MATRICES ===\n\n")

saveRDS(logcpm_221561, file.path(qc_dir, "GSE221561_logCPM_QC.rds"))
cat("Saved: GSE221561_logCPM_QC.rds\n")

saveRDS(logcpm_197677, file.path(qc_dir, "GSE197677_logCPM_QC.rds"))
cat("Saved: GSE197677_logCPM_QC.rds\n")

saveRDS(logcpm_160269, file.path(qc_dir, "GSE160269_logCPM_QC.rds"))
cat("Saved: GSE160269_logCPM_QC.rds\n")

cat("\n")
cat("NOTE: These logCPM matrices are QC_ONLY_NOT_DE_INPUT.\n")
cat("Do NOT use them for differential expression analysis.\n")

cat("\nDone.\n")
