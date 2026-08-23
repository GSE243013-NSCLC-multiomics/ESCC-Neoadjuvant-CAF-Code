#!/usr/bin/env Rscript
# =============================================================================
# STEPS 19D0-19D3: FREEZE DISCOVERY + VERIFY GSE221561 + GSE221561 DE + HARMONIZE
# =============================================================================

suppressPackageStartupMessages({
  library(edgeR)
})

base_dir <- "~/Downloads/ESCC Neoadjuvant"
de_dir <- file.path(base_dir, "06_tables/cross_dataset/epithelial_pseudobulk/GSE197677_DE")
repl_dir <- file.path(base_dir, "06_tables/cross_dataset/epithelial_pseudobulk/replication")
dir.create(repl_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# STEP 19D0: FREEZE GSE197677 DISCOVERY
# =============================================================================

cat("=== STEP 19D0: FREEZE GSE197677 DISCOVERY ===\n\n")

res_197677 <- read.csv(file.path(de_dir, "STEP19C_GSE197677_NACT_vs_nNACT_ALL_GENES.csv"), row.names = 1, check.names = FALSE)
loso_197677 <- read.csv(file.path(de_dir, "STEP19C_GSE197677_LOSO_top100_stability.csv"), stringsAsFactors = FALSE)
nominal_197677 <- read.csv(file.path(de_dir, "STEP19C_GSE197677_nominal_directional_candidates.csv"), row.names = 1, check.names = FALSE)

cat(sprintf("GSE197677 genes tested: %d\n", nrow(res_197677)))
cat(sprintf("FDR < 0.05: %d\n", sum(res_197677$FDR < 0.05)))
cat(sprintf("Lowest FDR: %.4f\n", min(res_197677$FDR)))
cat(sprintf("Nominal candidates: %d\n", nrow(nominal_197677)))
cat(sprintf("LOSO top-100 ROBUST: %d\n\n", sum(loso_197677$stability == "ROBUST_DIRECTION")))

# Rank genes: by PValue ascending, then abs(logFC) descending
res_197677$abs_logFC <- abs(res_197677$logFC)
res_197677_sorted <- res_197677[order(res_197677$PValue, -res_197677$abs_logFC), ]
res_197677_sorted$discovery_rank <- 1:nrow(res_197677_sorted)

# Build frozen freeze table
freeze <- data.frame(
  gene = rownames(res_197677_sorted),
  discovery_logFC = res_197677_sorted$logFC,
  discovery_PValue = res_197677_sorted$PValue,
  discovery_FDR = res_197677_sorted$FDR,
  discovery_direction = res_197677_sorted$direction,
  discovery_rank = res_197677_sorted$discovery_rank,
  LOSO_status = loso_197677$stability[match(rownames(res_197677_sorted), loso_197677$gene)],
  LOSO_direction_consistency = loso_197677$direction_consistency_fraction[match(rownames(res_197677_sorted), loso_197677$gene)],
  stringsAsFactors = FALSE
)

# Top 100 for primary replication
top100_freeze <- freeze[1:min(100, nrow(freeze)), ]

write.csv(freeze, file.path(repl_dir, "STEP19D_GSE197677_DISCOVERY_RANK_FREEZE.csv"), row.names = FALSE)
cat("Saved: STEP19D_GSE197677_DISCOVERY_RANK_FREEZE.csv\n\n")

# =============================================================================
# STEP 19D1: VERIFY GSE221561 INPUT
# =============================================================================

cat("=== STEP 19D1: VERIFY GSE221561 INPUT ===\n\n")

counts_221561 <- readRDS(file.path(base_dir, "06_tables/cross_dataset/epithelial_pseudobulk/counts/GSE221561_epithelial_pseudobulk_raw_counts.rds"))
meta_221561 <- read.csv(file.path(base_dir, "06_tables/cross_dataset/epithelial_pseudobulk/counts/GSE221561_epithelial_pseudobulk_sample_metadata.csv"), stringsAsFactors = FALSE)

cat(sprintf("Samples: %d\n", ncol(counts_221561)))
cat(sprintf("Neoadjuvant_treated: %d\n", sum(meta_221561$treatment_standardized == "Neoadjuvant_treated")))
cat(sprintf("Surgery_alone: %d\n", sum(meta_221561$treatment_standardized == "Surgery_alone")))

cat("\nNACT samples:\n")
cat(paste("  ", meta_221561$sample[meta_221561$treatment_standardized == "Neoadjuvant_treated"], collapse = "\n"), "\n")
cat("\nSurgery samples:\n")
cat(paste("  ", meta_221561$sample[meta_221561$treatment_standardized == "Surgery_alone"], collapse = "\n"), "\n")

# Checks
errors <- c()
if (ncol(counts_221561) != 9) errors <- c(errors, "Expected 9 samples")
if (sum(meta_221561$treatment_standardized == "Neoadjuvant_treated") != 7) errors <- c(errors, "Expected 7 Neoadjuvant_treated")
if (sum(meta_221561$treatment_standardized == "Surgery_alone") != 2) errors <- c(errors, "Expected 2 Surgery_alone")
if (any(duplicated(colnames(counts_221561)))) errors <- c(errors, "Duplicated sample IDs")
if (any(is.na(counts_221561))) errors <- c(errors, "NA counts")
if (any(counts_221561 < 0)) errors <- c(errors, "Negative counts")
if (!all(counts_221561 == floor(counts_221561))) errors <- c(errors, "Non-integer counts")

if (length(errors) > 0) {
  for (e in errors) cat("ERROR:", e, "\n")
  stop("Verification failed.")
}

cat("\nVerification PASSED.\n\n")

# =============================================================================
# STEP 19D2: GSE221561 edgeR QL MODEL
# =============================================================================

cat("=== STEP 19D2: GSE221561 edgeR QL MODEL ===\n\n")

dge_221561 <- DGEList(counts = counts_221561, samples = meta_221561)

# Treatment factor: reference = Surgery_alone, treatment = Neoadjuvant_treated
treatment_221561 <- factor(meta_221561$treatment_standardized,
                           levels = c("Surgery_alone", "Neoadjuvant_treated"))

design_221561 <- model.matrix(~ treatment_221561)
rownames(design_221561) <- meta_221561$sample

cat("Design matrix columns:\n")
print(colnames(design_221561))
cat("\nPositive logFC = Higher in Neoadjuvant_treated\n")
cat("Negative logFC = Higher in Surgery_alone\n\n")

# Filter
keep_221561 <- filterByExpr(dge_221561, design = design_221561)
cat(sprintf("Genes kept: %d / %d\n", sum(keep_221561), nrow(dge_221561)))

# Subset and normalize
dge_f_221561 <- dge_221561[keep_221561, , keep.lib.sizes = FALSE]
dge_f_221561 <- calcNormFactors(dge_f_221561, method = "TMM")

# Dispersion
dge_f_221561 <- estimateDisp(dge_f_221561, design_221561, robust = TRUE)
cat(sprintf("Common dispersion: %.6f\n", dge_f_221561$common.dispersion))
cat(sprintf("BCV: %.4f\n\n", sqrt(dge_f_221561$common.dispersion)))

# Fit
nact_coef_221561 <- grep("Neoadjuvant_treated", colnames(design_221561), value = TRUE)
cat(sprintf("Coefficient: %s\n", nact_coef_221561))

fit_221561 <- glmQLFit(dge_f_221561, design_221561, robust = TRUE)
qlf_221561 <- glmQLFTest(fit_221561, coef = nact_coef_221561)

# Extract results
res_221561 <- topTags(qlf_221561, n = Inf, sort.by = "PValue")$table
res_221561$direction <- ifelse(res_221561$logFC > 0, "Higher_in_Neoadjuvant_treated",
                        ifelse(res_221561$logFC < 0, "Higher_in_Surgery_alone", "No_direction"))

write.csv(res_221561, file.path(repl_dir, "STEP19D_GSE221561_NEOADJUVANT_vs_SURGERY_ALL_GENES.csv"), row.names = TRUE)
cat(sprintf("GSE221561 genes tested: %d\n", nrow(res_221561)))
cat(sprintf("FDR < 0.05: %d\n", sum(res_221561$FDR < 0.05)))
cat("Saved: STEP19D_GSE221561_NEOADJUVANT_vs_SURGERY_ALL_GENES.csv\n\n")

# =============================================================================
# STEP 19D3: GENE HARMONIZATION
# =============================================================================

cat("=== STEP 19D3: GENE HARMONIZATION ===\n\n")

common_genes_197677 <- rownames(res_197677)
common_genes_221561 <- rownames(res_221561)
common_tested <- intersect(common_genes_197677, common_genes_221561)

cat(sprintf("GSE197677 tested: %d\n", length(common_genes_197677)))
cat(sprintf("GSE221561 tested: %d\n", length(common_genes_221561)))
cat(sprintf("Common tested: %d\n\n", length(common_tested)))

harmonize <- data.frame(
  gene = common_tested,
  GSE197677_tested = TRUE,
  GSE221561_tested = TRUE,
  common_tested = TRUE,
  stringsAsFactors = FALSE
)

write.csv(harmonize, file.path(repl_dir, "STEP19D_common_tested_genes.csv"), row.names = FALSE)
cat("Saved: STEP19D_common_tested_genes.csv\n\n")

cat("Steps 19D0-19D3 complete.\n")
