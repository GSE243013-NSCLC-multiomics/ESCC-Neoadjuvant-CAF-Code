#!/usr/bin/env Rscript
# =============================================================================
# STEPS 19C13-19C14: MODEL LIMITATIONS + FREEZE DISCOVERY
# =============================================================================

base_dir <- "~/Downloads/ESCC Neoadjuvant"
de_dir <- file.path(base_dir, "06_tables/cross_dataset/epithelial_pseudobulk/GSE197677_DE")

# ---- STEP 19C13: MODEL LIMITATIONS NOTE ----

cat("=== STEP 19C13: MODEL LIMITATIONS NOTE ===\n\n")

limitations_text <- "# GSE197677 Model Limitations

## Statistical Design

- **Statistical unit**: Sample/patient (pseudobulk)
- **Comparison**: NACT vs nNACT tumor-specimen epithelial pseudobulk
- **Sample size**: NACT = 6, nNACT = 4 (total n = 10)
- **Cell-level replication**: Not used; all inference is at the sample level

## Covariates

- No additional clinical covariates were added because none were previously validated for inclusion in the frozen design
- The model contains only the treatment group indicator

## Batch Correction

- No cross-dataset batch correction was used
- This analysis is GSE197677 internal only

## Interpretation

- Results represent **association with treatment-group status**, not causal treatment effects
- The epithelial compartment is **not equivalent to CNV-confirmed malignant cells**
- All epithelial cells from tumor tissue were included regardless of CNV status

## Limitations

1. Small sample size (n = 6 vs 4) limits statistical power
2. No validation cohort within this dataset
3. Potential confounders (e.g., tumor stage, patient age) not adjusted
4. Biological coefficient of variation (BCV) reflects both true biological variation and unmeasured technical factors
5. FDR threshold controls false discovery rate but does not guarantee biological significance

## Reporting

- This analysis follows the EXPLORATORY OBSERVATIONAL DESCRIPTION ONLY framework
- No p-values are presented in the main text; FDR is reported for transparency
- Direction of effect is reported without causal attribution
"

writeLines(limitations_text, file.path(de_dir, "STEP19C_GSE197677_MODEL_LIMITATIONS.md"))
cat("Saved: STEP19C_GSE197677_MODEL_LIMITATIONS.md\n\n")

# ---- STEP 19C14: FREEZE DISCOVERY RESULT ----

cat("=== STEP 19C14: FREEZE DISCOVERY RESULT ===\n\n")

# Load summary data
res <- read.csv(file.path(de_dir, "STEP19C_GSE197677_NACT_vs_nNACT_ALL_GENES.csv"), row.names = 1, check.names = FALSE)

# Count genes
n_tested <- nrow(res)
n_fdr05 <- sum(res$FDR < 0.05)
n_fdr05_nact <- sum(res$FDR < 0.05 & res$logFC > 0)
n_fdr05_nnact <- sum(res$FDR < 0.05 & res$logFC < 0)

log2fc_thresh <- log2(1.5)
n_fdr05_effect <- sum(res$FDR < 0.05 & abs(res$logFC) >= log2fc_thresh)
n_fdr05_effect_nact <- sum(res$FDR < 0.05 & res$logFC > 0 & abs(res$logFC) >= log2fc_thresh)
n_fdr05_effect_nnact <- sum(res$FDR < 0.05 & res$logFC < 0 & abs(res$logFC) >= log2fc_thresh)

freeze <- data.frame(
  dataset = "GSE197677",
  contrast = "NACT_vs_nNACT",
  positive_logFC_definition = "Higher_in_NACT",
  samples_group1 = 6,
  samples_group2 = 4,
  genes_tested = n_tested,
  FDR05_genes = n_fdr05,
  FDR05_effect_genes = n_fdr05_effect,
  analysis_status = "PRIMARY_DISCOVERY_COMPLETE",
  statistical_unit = "SAMPLE_PATIENT",
  compartment = "TUMOR_TISSUE_EPITHELIAL",
  replication_dataset = "GSE221561",
  stringsAsFactors = FALSE
)

write.csv(freeze, file.path(de_dir, "STEP19C_GSE197677_DISCOVERY_FREEZE.csv"), row.names = FALSE)
cat("Saved: STEP19C_GSE197677_DISCOVERY_FREEZE.csv\n\n")

# ---- Print final output ----

cat("\n")
cat("============================================\n")
cat("STEP 19C COMPLETE\n")
cat("GSE197677 PRIMARY PSEUDOBULK DE\n")
cat("============================================\n\n")

cat("Comparison:\nNACT vs nNACT\n\n")

cat("Positive logFC:\nHigher in NACT\n\n")

cat("Samples:\n")
cat(sprintf("NACT                        : 6\n"))
cat(sprintf("nNACT                       : 4\n\n"))

# Load dispersion
dge_f <- readRDS(file.path(de_dir, "GSE197677_DGE_filtered_TMM_dispersion.rds"))
cat(sprintf("Genes tested                : %d\n", n_tested))
cat(sprintf("Common dispersion           : %.6f\n", dge_f$common.dispersion))
cat(sprintf("BCV                         : %.4f\n\n", sqrt(dge_f$common.dispersion)))

cat(sprintf("FDR < 0.05                  : %d\n", n_fdr05))
cat(sprintf("  Higher in NACT            : %d\n", n_fdr05_nact))
cat(sprintf("  Higher in nNACT           : %d\n\n", n_fdr05_nnact))

cat(sprintf("FDR < 0.05 + |FC| >= 1.5x  : %d\n", n_fdr05_effect))
cat(sprintf("  Higher in NACT            : %d\n", n_fdr05_effect_nact))
cat(sprintf("  Higher in nNACT           : %d\n\n", n_fdr05_effect_nnact))

cat(sprintf("Lowest FDR                  : %.6f\n", min(res$FDR)))
cat(sprintf("Largest |logFC|             : %.4f\n\n", max(abs(res$logFC))))

# LOSO
loso <- read.csv(file.path(de_dir, "STEP19C_GSE197677_LOSO_top100_stability.csv"), stringsAsFactors = FALSE)
cat("LOSO sensitivity:\n")
cat(sprintf("Top-100 robust direction    : %d\n", sum(loso$stability == "ROBUST_DIRECTION")))
cat(sprintf("Top-100 moderate direction  : %d\n", sum(loso$stability == "MODERATE_DIRECTION")))
cat(sprintf("Top-100 unstable direction  : %d\n\n", sum(loso$stability == "UNSTABLE_DIRECTION")))

# Top 10 genes
cat("--------------------------------------------\n\n")
cat("Top 10 genes by FDR:\n\n")

top10 <- head(res, 10)
for (i in 1:nrow(top10)) {
  gene <- rownames(top10)[i]
  loso_status <- loso$stability[loso$gene == gene]
  if (length(loso_status) == 0) loso_status <- "NA"
  cat(sprintf("%s\n  logFC: %.4f\n  FDR: %.6f\n  Direction: %s\n  LOSO: %s\n\n",
              gene, top10$logFC[i], top10$FDR[i], top10$direction[i], loso_status))
}

cat("--------------------------------------------\n\n")
cat("Statistical unit            : SAMPLE / PATIENT\n")
cat("Compartment                 : TUMOR-TISSUE EPITHELIAL\n")
cat("Malignant-specific claim    : NO\n\n")

cat("Samples removed             : NO\n")
cat("Cells removed               : NO\n")
cat("Cross-dataset integration   : NO\n")
cat("Pathway enrichment          : NO\n")
cat("Replication DE              : NO\n\n")

cat("STOP HERE.\n")
cat("DO NOT START STEP 19D AUTOMATICALLY.\n")
