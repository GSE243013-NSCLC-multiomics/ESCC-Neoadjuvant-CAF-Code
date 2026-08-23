#!/usr/bin/env Rscript
# =============================================================================
# STEP 19D-FIX: CONTRAST ORIENTATION AUDIT + CANONICAL REPLICATION REBUILD
# =============================================================================

suppressPackageStartupMessages({
  library(edgeR)
})

base_dir <- "~/Downloads/ESCC Neoadjuvant"
repl_dir <- file.path(base_dir, "06_tables/cross_dataset/epithelial_pseudobulk/replication")
de_dir <- file.path(base_dir, "06_tables/cross_dataset/epithelial_pseudobulk/GSE197677_DE")

# Check for global 'file' object issue
if (exists("file", envir = .GlobalEnv, inherits = FALSE)) {
  obj <- get("file", envir = .GlobalEnv)
  if (!is.function(obj)) {
    rm(file, envir = .GlobalEnv)
    cat("Removed non-function 'file' from global environment.\n")
  }
}

# =============================================================================
# 19D-FIX0: FREEZE CURRENT FILES
# =============================================================================

cat("=== 19D-FIX0: FREEZE CURRENT FILES ===\n\n")

old_repl <- read.csv(file.path(repl_dir, "STEP19D_GSE221561_NEOADJUVANT_vs_SURGERY_ALL_GENES.csv"),
                     row.names = 1, check.names = FALSE)

utils::write.csv(old_repl, file.path(repl_dir, "STEP19D_GSE221561_REPLICATION_DE_PRE_ORIENTATION_FIX.csv"))
stopifnot(file.exists(file.path(repl_dir, "STEP19D_GSE221561_REPLICATION_DE_PRE_ORIENTATION_FIX.csv")))
cat("Saved: STEP19D_GSE221561_REPLICATION_DE_PRE_ORIENTATION_FIX.csv\n\n")

# =============================================================================
# 19D-FIX1: AUDIT GSE221561 MODEL ORIENTATION
# =============================================================================

cat("=== 19D-FIX1: AUDIT GSE221561 MODEL ORIENTATION ===\n\n")

raw_counts <- readRDS(file.path(base_dir, "06_tables/cross_dataset/epithelial_pseudobulk/counts/GSE221561_epithelial_pseudobulk_raw_counts.rds"))
meta <- read.csv(file.path(base_dir, "06_tables/cross_dataset/epithelial_pseudobulk/counts/GSE221561_epithelial_pseudobulk_sample_metadata.csv"), stringsAsFactors = FALSE)

cat("Unique treatment values:\n")
print(unique(meta$treatment_standardized))
cat("\nTreatment table:\n")
print(table(meta$treatment_standardized))

# Create factor with explicit levels
treatment <- factor(meta$treatment_standardized,
                    levels = c("Surgery_alone", "Neoadjuvant_treated"))
cat("\nFactor levels:\n")
print(levels(treatment))

# No-intercept design
design <- model.matrix(~ 0 + treatment)
colnames(design) <- c("Surgery_alone", "Neoadjuvant_treated")
rownames(design) <- meta$sample

cat("\nDesign matrix:\n")
print(design)
cat("\nDesign colnames:\n")
print(colnames(design))

# Canonical contrast
contrast_vec <- limma::makeContrasts(
  Neoadjuvant_treated - Surgery_alone,
  levels = design
)

cat("\nCanonical contrast vector:\n")
print(contrast_vec)

cat("\n")
cat("CANONICAL POSITIVE LOGFC = HIGHER IN Neoadjuvant_treated\n")
cat("CANONICAL NEGATIVE LOGFC = HIGHER IN Surgery_alone\n\n")

# =============================================================================
# 19D-FIX2: RERUN GSE221561 edgeR MODEL
# =============================================================================

cat("=== 19D-FIX2: RERUN GSE221561 edgeR MODEL ===\n\n")

dge <- DGEList(counts = raw_counts, samples = meta)

keep <- filterByExpr(dge, design = design)
cat(sprintf("Genes kept: %d / %d\n", sum(keep), nrow(dge)))

dge_f <- dge[keep, , keep.lib.sizes = FALSE]
dge_f <- calcNormFactors(dge_f, method = "TMM")
dge_f <- estimateDisp(dge_f, design, robust = TRUE)

cat(sprintf("Common dispersion: %.6f\n", dge_f$common.dispersion))
cat(sprintf("BCV: %.4f\n\n", sqrt(dge_f$common.dispersion)))

fit <- glmQLFit(dge_f, design, robust = TRUE)
qlf <- glmQLFTest(fit, contrast = contrast_vec)

res_canonical <- topTags(qlf, n = Inf, sort.by = "PValue")$table
res_canonical$direction <- ifelse(res_canonical$logFC > 0, "Higher_in_Neoadjuvant_treated",
                           ifelse(res_canonical$logFC < 0, "Higher_in_Surgery_alone", "No_direction"))

canonical_path <- file.path(repl_dir, "STEP19D_GSE221561_NEOADJUVANT_MINUS_SURGERY_CANONICAL.csv")
utils::write.csv(res_canonical, canonical_path)
stopifnot(file.exists(canonical_path))
cat(sprintf("Canonical DE genes: %d\n", nrow(res_canonical)))
cat(sprintf("FDR < 0.05: %d\n\n", sum(res_canonical$FDR < 0.05)))

# =============================================================================
# 19D-FIX3: DIRECT SIGN-INVERSION DIAGNOSTIC
# =============================================================================

cat("=== 19D-FIX3: SIGN-INVERSION DIAGNOSTIC ===\n\n")

res_197677 <- read.csv(file.path(de_dir, "STEP19C_GSE197677_NACT_vs_nNACT_ALL_GENES.csv"),
                       row.names = 1, check.names = FALSE)

common_genes <- intersect(rownames(res_197677), rownames(res_canonical))
cat(sprintf("Common tested genes: %d\n\n", length(common_genes)))

# A. Old orientation
rho_old <- cor(res_197677$logFC[common_genes], old_repl$logFC[common_genes], method = "spearman")
same_old <- sum(sign(res_197677$logFC[common_genes]) == sign(old_repl$logFC[common_genes]))
cat(sprintf("A. Discovery vs OLD:     rho = %.4f, same_sign = %d / %d = %.3f\n",
            rho_old, same_old, length(common_genes), same_old/length(common_genes)))

# B. Flipped old
rho_flipped <- cor(res_197677$logFC[common_genes], -old_repl$logFC[common_genes], method = "spearman")
same_flipped <- sum(sign(res_197677$logFC[common_genes]) == sign(-old_repl$logFC[common_genes]))
cat(sprintf("B. Discovery vs -OLD:    rho = %.4f, same_sign = %d / %d = %.3f\n",
            rho_flipped, same_flipped, length(common_genes), same_flipped/length(common_genes)))

# C. New canonical
rho_new <- cor(res_197677$logFC[common_genes], res_canonical$logFC[common_genes], method = "spearman")
same_new <- sum(sign(res_197677$logFC[common_genes]) == sign(res_canonical$logFC[common_genes]))
cat(sprintf("C. Discovery vs CANONICAL: rho = %.4f, same_sign = %d / %d = %.3f\n\n",
            rho_new, same_new, length(common_genes), same_new/length(common_genes)))

# =============================================================================
# 19D-FIX4: TOP-100 SIGN TEST
# =============================================================================

cat("=== 19D-FIX4: TOP-100 SIGN TEST ===\n\n")

freeze <- read.csv(file.path(repl_dir, "STEP19D_GSE197677_DISCOVERY_RANK_FREEZE.csv"), stringsAsFactors = FALSE)
top100 <- freeze[1:min(100, nrow(freeze)), ]

# Match to all three orientations
match_old <- match(top100$gene, rownames(old_repl))
match_new <- match(top100$gene, rownames(res_canonical))

top100$old_logFC <- old_repl$logFC[match_old]
top100$new_logFC <- res_canonical$logFC[match_new]
top100$tested_old <- !is.na(top100$old_logFC)
top100$tested_new <- !is.na(top100$new_logFC)

for (k in c(25, 50, 100)) {
  sub <- top100[1:k, ]
  tested_idx <- sub$tested_new
  n_t <- sum(tested_idx)

  n_same_old <- sum(sub$discovery_direction[tested_idx] == old_repl$direction[match(sub$gene[tested_idx], rownames(old_repl))])
  n_same_flipped <- sum(sub$discovery_direction[tested_idx] == ifelse(old_repl$logFC[match(sub$gene[tested_idx], rownames(old_repl))] > 0,
                                                                      "Higher_in_Surgery_alone",
                                                                      "Higher_in_Neoadjuvant_treated"))
  n_same_new <- sum(sub$discovery_direction[tested_idx] == res_canonical$direction[match(sub$gene[tested_idx], rownames(res_canonical))])

  cat(sprintf("Top %d (tested=%d):\n", k, n_t))
  cat(sprintf("  OLD same dir:          %d / %d = %.3f\n", n_same_old, n_t, n_same_old/n_t))
  cat(sprintf("  -OLD same dir:         %d / %d = %.3f\n", n_same_flipped, n_t, n_same_flipped/n_t))
  cat(sprintf("  CANONICAL same dir:    %d / %d = %.3f\n\n", n_same_new, n_t, n_same_new/n_t))
}

# =============================================================================
# 19D-FIX5: CLASSIFY THE PROBLEM
# =============================================================================

cat("=== 19D-FIX5: CLASSIFY THE PROBLEM ===\n\n")

# Use top100 results
top100_tested <- top100[top100$tested_new, ]
n_tested_100 <- nrow(top100_tested)

# OLD concordance
n_same_old_100 <- sum(top100_tested$discovery_direction == old_repl$direction[match(top100_tested$gene, rownames(old_repl))])
conc_old <- n_same_old_100 / n_tested_100

# Flipped concordance
n_same_flip_100 <- sum(top100_tested$discovery_direction == ifelse(old_repl$logFC[match(top100_tested$gene, rownames(old_repl))] > 0,
                                                                    "Higher_in_Surgery_alone",
                                                                    "Higher_in_Neoadjuvant_treated"))
conc_flip <- n_same_flip_100 / n_tested_100

# Canonical concordance
n_same_new_100 <- sum(top100_tested$discovery_direction == res_canonical$direction[match(top100_tested$gene, rownames(res_canonical))])
conc_new <- n_same_new_100 / n_tested_100

cat(sprintf("OLD concordance:     %d / %d = %.3f\n", n_same_old_100, n_tested_100, conc_old))
cat(sprintf("-OLD concordance:    %d / %d = %.3f\n", n_same_flip_100, n_tested_100, conc_flip))
cat(sprintf("CANONICAL concordance: %d / %d = %.3f\n\n", n_same_new_100, n_tested_100, conc_new))

# Classification
if (conc_old <= 0.10 && conc_flip >= 0.90 && conc_new >= 0.90) {
  status <- "GLOBAL_SIGN_INVERSION_CONFIRMED"
  prev_valid <- "NO"
} else if (conc_old <= 0.10 && conc_flip >= 0.70) {
  status <- "PARTIAL_ORIENTATION_PROBLEM"
  prev_valid <- "NO"
} else if (conc_new < 0.50) {
  status <- "TRUE_DIRECTIONAL_DISCORDANCE"
  prev_valid <- "YES (but reflects biological discordance)"
} else {
  status <- "UNRESOLVED"
  prev_valid <- "UNCLEAR"
}

cat(sprintf("Orientation audit status: %s\n", status))
cat(sprintf("Previous 0/99 biological interpretation valid: %s\n\n", prev_valid))

# =============================================================================
# 19D-FIX6: CHECK SPECIFIC GENES
# =============================================================================

cat("=== 19D-FIX6: CHECK SPECIFIC GENES ===\n\n")

for (gene in c("RCN3", "SDC2", "SENP5")) {
  disc_fc <- res_197677$logFC[gene]
  old_fc <- old_repl$logFC[gene]
  new_fc <- res_canonical$logFC[gene]

  cat(sprintf("%s:\n", gene))
  cat(sprintf("  Discovery logFC:           %.4f\n", disc_fc))
  cat(sprintf("  OLD replication logFC:      %.4f\n", old_fc))
  cat(sprintf("  -OLD replication logFC:     %.4f\n", -old_fc))
  cat(sprintf("  CANONICAL replication logFC: %.4f\n", new_fc))
  cat(sprintf("  Same direction (CANONICAL): %s\n\n",
              ifelse(sign(disc_fc) == sign(new_fc), "YES", "NO")))
}

# =============================================================================
# 19D-FIX8: PROVENANCE NOTE
# =============================================================================

cat("=== 19D-FIX8: PROVENANCE NOTE ===\n\n")

provenance <- sprintf("# Step 19D Contrast Orientation Audit

## Date
%s

## Background

Step 19D initially reported 0/99 direction concordance between GSE197677
discovery and GSE221561 replication top-100 genes. This was flagged as
suspicious for a global contrast sign inversion.

## Contrast Definitions

1. Step 19C (GSE197677): positive logFC = Higher in NACT
   - Design: model.matrix(~ treatment), reference = nNACT
   - Coefficient: treatmentNeoadjuvant_chemotherapy

2. Step 19D original (GSE221561): positive logFC = Higher in Neoadjuvant_treated
   - Design: model.matrix(~ treatment), reference = Surgery_alone
   - Coefficient: treatment_221561Neoadjuvant_treated

3. Step 19D canonical (GSE221561): positive logFC = Higher in Neoadjuvant_treated
   - Design: model.matrix(~ 0 + treatment), no intercept
   - Contrast: Neoadjuvant_treated - Surgery_alone

## Pre-fix Results

- Top100 concordance (OLD): %d / %d = %.3f
- Top100 concordance (-OLD): %d / %d = %.3f
- Top100 concordance (CANONICAL): %d / %d = %.3f

## Genome-wide Spearman

- Discovery vs OLD: %.4f
- Discovery vs -OLD: %.4f
- Discovery vs CANONICAL: %.4f

## Classification

Status: %s
Previous 0/99 interpretation valid: %s

## Note

The old coefficient-based contrast and the new no-intercept contrast
should yield identical results if properly specified. A discrepancy
indicates the original model had an orientation issue.

%s
",
Sys.Date(),
n_same_old_100, n_tested_100, conc_old,
n_same_flip_100, n_tested_100, conc_flip,
n_same_new_100, n_tested_100, conc_new,
rho_old, rho_flipped, rho_new,
status, prev_valid,
if (status == "GLOBAL_SIGN_INVERSION_CONFIRMED") {
  "The previous 0/99 result is INVALID_DUE_TO_CONTRAST_ORIENTATION and must not be interpreted as biological discordance."
} else {
  "Further investigation required."
})

writeLines(provenance, file.path(repl_dir, "STEP19D_CONTRAST_ORIENTATION_AUDIT.md"))
cat("Saved: STEP19D_CONTRAST_ORIENTATION_AUDIT.md\n\n")

# =============================================================================
# FINAL OUTPUT
# =============================================================================

cat("\n")
cat("============================================\n")
cat("STEP 19D CONTRAST AUDIT COMPLETE\n")
cat("============================================\n\n")

cat("GSE197677 positive logFC:\nHigher in NACT\n\n")

cat("GSE221561 intended positive logFC:\nHigher in Neoadjuvant_treated\n\n")

cat("Old GSE221561 coefficient/contrast:\n")
cat("  model.matrix(~ treatment), coef = treatmentNeoadjuvant_treated\n\n")

cat("Canonical contrast:\n")
cat("  Neoadjuvant_treated - Surgery_alone (no-intercept design)\n\n")

cat("--------------------------------------------\n\n")

cat("Genome-wide Spearman:\n")
cat(sprintf("  Discovery vs OLD:      rho = %.4f\n", rho_old))
cat(sprintf("  Discovery vs -OLD:     rho = %.4f\n", rho_flipped))
cat(sprintf("  Discovery vs CANONICAL: rho = %.4f\n\n", rho_new))

cat("Frozen Top100:\n")
cat(sprintf("  Tested: %d\n", n_tested_100))
cat(sprintf("  OLD same direction:      %d / %d = %.3f\n", n_same_old_100, n_tested_100, conc_old))
cat(sprintf("  -OLD same direction:     %d / %d = %.3f\n", n_same_flip_100, n_tested_100, conc_flip))
cat(sprintf("  CANONICAL same direction: %d / %d = %.3f\n\n", n_same_new_100, n_tested_100, conc_new))

cat("RCN3:\n")
cat(sprintf("  Discovery logFC:            %.4f\n", res_197677$logFC["RCN3"]))
cat(sprintf("  OLD replication logFC:       %.4f\n", old_repl$logFC["RCN3"]))
cat(sprintf("  CANONICAL replication logFC: %.4f\n\n", res_canonical$logFC["RCN3"]))

cat("SDC2:\n")
cat(sprintf("  Discovery logFC:            %.4f\n", res_197677$logFC["SDC2"]))
cat(sprintf("  OLD replication logFC:       %.4f\n", old_repl$logFC["SDC2"]))
cat(sprintf("  CANONICAL replication logFC: %.4f\n\n", res_canonical$logFC["SDC2"]))

cat("SENP5:\n")
cat(sprintf("  Discovery logFC:            %.4f\n", res_197677$logFC["SENP5"]))
cat(sprintf("  OLD replication logFC:       %.4f\n", old_repl$logFC["SENP5"]))
cat(sprintf("  CANONICAL replication logFC: %.4f\n\n", res_canonical$logFC["SENP5"]))

cat("--------------------------------------------\n\n")
cat(sprintf("Orientation audit status: %s\n", status))
cat(sprintf("Previous 0/99 biological interpretation valid: %s\n\n", prev_valid))

cat("STOP HERE.\n")
cat("DO NOT CONTINUE STEP 19D8 AUTOMATICALLY.\n")
