#!/usr/bin/env Rscript
# =============================================================================
# STEPS 19C1-19C8: GSE197677 PRIMARY DE ANALYSIS
# =============================================================================

suppressPackageStartupMessages({
  library(edgeR)
})

base_dir <- "~/Downloads/ESCC Neoadjuvant"
de_dir <- file.path(base_dir, "06_tables/cross_dataset/epithelial_pseudobulk/GSE197677_DE")
dir.create(de_dir, showWarnings = FALSE, recursive = TRUE)

cat("Loading inputs...\n")

raw_counts <- readRDS(file.path(base_dir, "06_tables/cross_dataset/epithelial_pseudobulk/counts/GSE197677_epithelial_pseudobulk_raw_counts.rds"))
meta <- read.csv(file.path(base_dir, "06_tables/cross_dataset/epithelial_pseudobulk/counts/GSE197677_epithelial_pseudobulk_sample_metadata.csv"), stringsAsFactors = FALSE)
outlier <- read.csv(file.path(base_dir, "06_tables/cross_dataset/epithelial_pseudobulk/QC/STEP19B_sample_outlier_audit.csv"), stringsAsFactors = FALSE)

cat("Loaded.\n\n")

# =============================================================================
# STEP 19C1: FREEZE SAMPLE TABLE
# =============================================================================

cat("=== STEP 19C1: FREEZE SAMPLE TABLE ===\n\n")

outlier_197677 <- outlier[outlier$dataset == "GSE197677", ]

sample_table <- data.frame(
  sample = meta$sample,
  patient = meta$patient,
  treatment = meta$treatment_standardized,
  epithelial_cells = meta$epithelial_cells,
  library_size = colSums(raw_counts)[meta$sample],
  TMM_factor = NA,  # will fill after TMM
  QC_status = outlier_197677$QC_status[match(meta$sample, outlier_197677$sample)],
  included_in_DE = TRUE,
  stringsAsFactors = FALSE
)

cat("NACT samples:\n")
cat(paste("  ", sample_table$sample[sample_table$treatment == "Neoadjuvant_chemotherapy"], collapse = "\n"), "\n\n")

cat("nNACT samples:\n")
cat(paste("  ", sample_table$sample[sample_table$treatment == "No_neoadjuvant_chemotherapy"], collapse = "\n"), "\n\n")

write.csv(sample_table, file.path(de_dir, "STEP19C_GSE197677_DE_SAMPLE_TABLE.csv"), row.names = FALSE)
cat("Saved: STEP19C_GSE197677_DE_SAMPLE_TABLE.csv\n\n")

# =============================================================================
# STEP 19C2: REBUILD FINAL DGE OBJECT
# =============================================================================

cat("=== STEP 19C2: REBUILD FINAL DGE OBJECT ===\n\n")

dge <- DGEList(counts = raw_counts, samples = meta)

# Treatment factor: reference = nNACT, treatment = NACT
treatment <- factor(meta$treatment_standardized,
                    levels = c("No_neoadjuvant_chemotherapy", "Neoadjuvant_chemotherapy"))

# Design matrix
design <- model.matrix(~ treatment)

cat("Design matrix columns:\n")
print(colnames(design))
cat("\nCoefficient for NACT vs nNACT: 'treatmentNeoadjuvant_chemotherapy'\n")
cat("Positive logFC = Higher in NACT\n")
cat("Negative logFC = Higher in nNACT\n\n")

# =============================================================================
# STEP 19C3: EXPRESSION FILTER
# =============================================================================

cat("=== STEP 19C3: EXPRESSION FILTER ===\n\n")

keep <- filterByExpr(dge, design = design)
cat(sprintf("Genes kept: %d / %d\n", sum(keep), nrow(dge)))

# Compare with Step 19B
filter_19b <- read.csv(file.path(base_dir, "06_tables/cross_dataset/epithelial_pseudobulk/QC/STEP19B_expression_filter_audit.csv"), stringsAsFactors = FALSE)
genes_19b <- filter_19b$genes_keep[filter_19b$dataset == "GSE197677"]
cat(sprintf("Step 19B filter: %d genes\n", genes_19b))
cat(sprintf("Difference: %d genes\n", sum(keep) - genes_19b))

if (abs(sum(keep) - genes_19b) > 100) {
  cat("WARNING: Large difference from Step 19B. Review before proceeding.\n")
}

# Subset and recalculate
dge_f <- dge[keep, , keep.lib.sizes = FALSE]
dge_f <- calcNormFactors(dge_f, method = "TMM")

# Update sample table with TMM factors
sample_table$TMM_factor <- dge_f$samples$norm.factors
write.csv(sample_table, file.path(de_dir, "STEP19C_GSE197677_DE_SAMPLE_TABLE.csv"), row.names = FALSE)

cat("\n")

# =============================================================================
# STEP 19C4: DISPERSION ESTIMATION
# =============================================================================

cat("=== STEP 19C4: DISPERSION ESTIMATION ===\n\n")

dge_f <- estimateDisp(dge_f, design, robust = TRUE)

common_disp <- dge_f$common.dispersion
bcv <- sqrt(common_disp)

cat(sprintf("Common dispersion: %.6f\n", common_disp))
cat(sprintf("BCV (biological coefficient of variation): %.4f\n", bcv))
cat("\nNOTE: Do not interpret BCV as treatment effect.\n\n")

saveRDS(dge_f, file.path(de_dir, "GSE197677_DGE_filtered_TMM_dispersion.rds"))
cat("Saved: GSE197677_DGE_filtered_TMM_dispersion.rds\n\n")

# =============================================================================
# STEP 19C5: QUASI-LIKELIHOOD MODEL
# =============================================================================

cat("=== STEP 19C5: QUASI-LIKELIHOOD MODEL ===\n\n")

fit <- glmQLFit(dge_f, design, robust = TRUE)

# Identify the NACT coefficient
nact_coef <- "treatmentNeoadjuvant_chemotherapy"
cat(sprintf("Coefficient for testing: %s\n", nact_coef))
cat("This tests: NACT - nNACT\n")
cat("Positive logFC = Higher in NACT\n\n")

qlf <- glmQLFTest(fit, coef = nact_coef)

cat("QL F-test complete.\n\n")

# =============================================================================
# STEP 19C6: FULL DE TABLE
# =============================================================================

cat("=== STEP 19C6: FULL DE TABLE ===\n\n")

res <- topTags(qlf, n = Inf, sort.by = "PValue")$table

# Add direction column
res$direction <- ifelse(res$logFC > 0, "Higher_in_NACT",
                 ifelse(res$logFC < 0, "Higher_in_nNACT", "No_direction"))

write.csv(res, file.path(de_dir, "STEP19C_GSE197677_NACT_vs_nNACT_ALL_GENES.csv"), row.names = TRUE)
cat(sprintf("Full DE table: %d genes\n", nrow(res)))
cat("Saved: STEP19C_GSE197677_NACT_vs_nNACT_ALL_GENES.csv\n\n")

# =============================================================================
# STEP 19C7: SIGNIFICANCE TIERS
# =============================================================================

cat("=== STEP 19C7: SIGNIFICANCE TIERS ===\n\n")

log2fc_threshold <- log2(1.5)
cat(sprintf("Effect size threshold: log2(1.5) = %.4f\n\n", log2fc_threshold))

# Tier A: FDR < 0.05
tier_a <- res[res$FDR < 0.05, ]
cat(sprintf("TIER A (FDR < 0.05): %d genes\n", nrow(tier_a)))
cat(sprintf("  Higher in NACT: %d\n", sum(tier_a$direction == "Higher_in_NACT")))
cat(sprintf("  Higher in nNACT: %d\n", sum(tier_a$direction == "Higher_in_nNACT")))

# Tier B: FDR < 0.05 AND abs(logFC) >= log2(1.5)
tier_b <- res[res$FDR < 0.05 & abs(res$logFC) >= log2fc_threshold, ]
cat(sprintf("\nTIER B (FDR < 0.05 + |FC| >= 1.5x): %d genes\n", nrow(tier_b)))
cat(sprintf("  Higher in NACT: %d\n", sum(tier_b$direction == "Higher_in_NACT")))
cat(sprintf("  Higher in nNACT: %d\n", sum(tier_b$direction == "Higher_in_nNACT")))

# Tier C: PValue < 0.05 AND abs(logFC) >= log2(1.5) AND FDR >= 0.05
tier_c <- res[res$PValue < 0.05 & abs(res$logFC) >= log2fc_threshold & res$FDR >= 0.05, ]
cat(sprintf("\nTIER C (Nominal exploratory): %d genes\n", nrow(tier_c)))
cat("  NOTE: NOT statistically significant. For cross-dataset directional replication only.\n")

# Save tier tables
write.csv(tier_a, file.path(de_dir, "STEP19C_GSE197677_FDR05.csv"), row.names = TRUE)
write.csv(tier_b, file.path(de_dir, "STEP19C_GSE197677_FDR05_log2FC1p5.csv"), row.names = TRUE)
write.csv(tier_c, file.path(de_dir, "STEP19C_GSE197677_nominal_directional_candidates.csv"), row.names = TRUE)

cat("\nSaved tier tables.\n\n")

# =============================================================================
# STEP 19C8: DE SUMMARY
# =============================================================================

cat("=== STEP 19C8: DE SUMMARY ===\n\n")

summary_stats <- data.frame(
  metric = c("genes_tested",
             "FDR05_total", "FDR05_higher_NACT", "FDR05_higher_nNACT",
             "FDR05_effect_total", "FDR05_effect_higher_NACT", "FDR05_effect_higher_nNACT",
             "nominal_effect_total",
             "median_abs_logFC", "max_abs_logFC", "min_FDR"),
  value = c(nrow(res),
            nrow(tier_a), sum(tier_a$direction == "Higher_in_NACT"), sum(tier_a$direction == "Higher_in_nNACT"),
            nrow(tier_b), sum(tier_b$direction == "Higher_in_NACT"), sum(tier_b$direction == "Higher_in_nNACT"),
            nrow(tier_c),
            round(median(abs(res$logFC)), 4),
            round(max(abs(res$logFC)), 4),
            round(min(res$FDR), 6)),
  stringsAsFactors = FALSE
)

write.csv(summary_stats, file.path(de_dir, "STEP19C_GSE197677_DE_SUMMARY.csv"), row.names = FALSE)
cat("Saved: STEP19C_GSE197677_DE_SUMMARY.csv\n\n")

# Print summary
for (i in 1:nrow(summary_stats)) {
  cat(sprintf("  %s: %s\n", summary_stats$metric[i], summary_stats$value[i]))
}

cat("\nSteps 19C1-19C8 complete.\n")
