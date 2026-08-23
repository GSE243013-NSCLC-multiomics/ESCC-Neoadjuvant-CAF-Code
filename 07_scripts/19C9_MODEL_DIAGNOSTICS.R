#!/usr/bin/env Rscript
# =============================================================================
# STEP 19C9: MODEL DIAGNOSTICS
# =============================================================================

suppressPackageStartupMessages({
  library(edgeR)
})

base_dir <- "~/Downloads/ESCC Neoadjuvant"
de_dir <- file.path(base_dir, "06_tables/cross_dataset/epithelial_pseudobulk/GSE197677_DE")
fig_dir <- file.path(base_dir, "05_figures/cross_dataset/epithelial_pseudobulk/GSE197677_DE")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

cat("Loading DGE object and fit...\n")

dge_f <- readRDS(file.path(de_dir, "GSE197677_DGE_filtered_TMM_dispersion.rds"))
res <- read.csv(file.path(de_dir, "STEP19C_GSE197677_NACT_vs_nNACT_ALL_GENES.csv"), row.names = 1, check.names = FALSE)

# Rebuild design and fit
meta <- read.csv(file.path(base_dir, "06_tables/cross_dataset/epithelial_pseudobulk/counts/GSE197677_epithelial_pseudobulk_sample_metadata.csv"), stringsAsFactors = FALSE)
treatment <- factor(meta$treatment_standardized, levels = c("No_neoadjuvant_chemotherapy", "Neoadjuvant_chemotherapy"))
design <- model.matrix(~ treatment)
rownames(design) <- meta$sample

fit <- glmQLFit(dge_f, design, robust = TRUE)
qlf <- glmQLFTest(fit, coef = "treatmentNeoadjuvant_chemotherapy")

cat("Loaded.\n\n")

log2fc_threshold <- log2(1.5)

# ---- 1. BCV Plot ----

cat("=== BCV PLOT ===\n\n")

png(file.path(fig_dir, "STEP19C_GSE197677_BCV_plot.png"), width = 7, height = 5, units = "in", res = 150)
plotBCV(dge_f, main = "GSE197677 Biological Coefficient of Variation")
dev.off()
cat("Saved: STEP19C_GSE197677_BCV_plot.png\n\n")

# ---- 2. QL Dispersion Plot ----

cat("=== QL DISPERSION PLOT ===\n\n")

png(file.path(fig_dir, "STEP19C_GSE197677_QL_dispersion_plot.png"), width = 7, height = 5, units = "in", res = 150)
plotQLDisp(fit, main = "GSE197677 Quasi-Likelihood Dispersion")
dev.off()
cat("Saved: STEP19C_GSE197677_QL_dispersion_plot.png\n\n")

# ---- 3. MD Plot ----

cat("=== MD PLOT ===\n\n")

png(file.path(fig_dir, "STEP19C_GSE197677_MD_plot.png"), width = 7, height = 5, units = "in", res = 150)
plotMD(qlf, main = "GSE197677 Mean-Difference Plot",
       hl.col = c("red", "blue"),
       legend = "topright")
abline(h = c(log2fc_threshold, -log2fc_threshold), lty = 2, col = "grey50")
dev.off()
cat("Saved: STEP19C_GSE197677_MD_plot.png\n\n")

# ---- 4. Volcano Plot ----

cat("=== VOLCANO PLOT ===\n\n")

# Define categories
fdr_threshold <- 0.05
volcano_category <- rep("Not significant", nrow(res))
volcano_category[res$FDR < fdr_threshold & abs(res$logFC) >= log2fc_threshold] <- "FDR < 0.05 + |FC| >= 1.5x"
volcano_category[res$FDR < fdr_threshold & abs(res$logFC) < log2fc_threshold] <- "FDR < 0.05 only"
volcano_category[res$FDR >= fdr_threshold & res$PValue < 0.05 & abs(res$logFC) >= log2fc_threshold] <- "Nominal (P < 0.05 + |FC| >= 1.5x)"

res$volcano_category <- volcano_category

# Order for plotting
cat("Volcano categories:\n")
print(table(volcano_category))

png(file.path(fig_dir, "STEP19C_GSE197677_volcano.png"), width = 8, height = 6, units = "in", res = 150)

# Color palette
cols <- c("FDR < 0.05 + |FC| >= 1.5x" = "#E41A1C",
          "FDR < 0.05 only" = "#FF7F00",
          "Nominal (P < 0.05 + |FC| >= 1.5x)" = "#377EB8",
          "Not significant" = "grey70")

# Plot
plot(res$logFC, -log10(res$FDR),
     col = cols[volcano_category],
     pch = 16, cex = 0.6,
     xlab = "log2 Fold Change (NACT vs nNACT)",
     ylab = "-log10(FDR)",
     main = "GSE197677 NACT vs nNACT Volcano Plot",
     cex.main = 1.2)

# Threshold lines
abline(v = c(log2fc_threshold, -log2fc_threshold), lty = 2, col = "grey50")
abline(h = -log10(fdr_threshold), lty = 2, col = "grey50")

# Legend
legend("topright", legend = names(cols), col = cols, pch = 16, cex = 0.7)

dev.off()
cat("Saved: STEP19C_GSE197677_volcano.png\n\n")

cat("=== MODEL DIAGNOSTICS COMPLETE ===\n")
