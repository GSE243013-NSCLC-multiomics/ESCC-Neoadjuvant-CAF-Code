#!/usr/bin/env Rscript
# =============================================================================
# STEP 19B10: CROSS-DATASET UNCORRECTED PCA
# =============================================================================

suppressPackageStartupMessages({
  library(edgeR)
})

base_dir <- "~/Downloads/ESCC Neoadjuvant"
qc_dir <- file.path(base_dir, "06_tables/cross_dataset/epithelial_pseudobulk/QC")
fig_dir <- file.path(base_dir, "05_figures/cross_dataset/epithelial_pseudobulk/QC")

cat("Loading data...\n")

logcpm_221561 <- readRDS(file.path(qc_dir, "GSE221561_logCPM_QC.rds"))
logcpm_197677 <- readRDS(file.path(qc_dir, "GSE197677_logCPM_QC.rds"))
logcpm_160269 <- readRDS(file.path(qc_dir, "GSE160269_logCPM_QC.rds"))

common_genes <- scan(file.path(base_dir, "06_tables/cross_dataset/epithelial_pseudobulk/STEP19A_three_dataset_common_gene_symbols.txt"), character(), quiet = TRUE)

meta_221561 <- read.csv(file.path(base_dir, "06_tables/cross_dataset/epithelial_pseudobulk/counts/GSE221561_epithelial_pseudobulk_sample_metadata.csv"), stringsAsFactors = FALSE)
meta_197677 <- read.csv(file.path(base_dir, "06_tables/cross_dataset/epithelial_pseudobulk/counts/GSE197677_epithelial_pseudobulk_sample_metadata.csv"), stringsAsFactors = FALSE)
meta_160269 <- read.csv(file.path(base_dir, "06_tables/cross_dataset/epithelial_pseudobulk/counts/GSE160269_epithelial_pseudobulk_sample_metadata.csv"), stringsAsFactors = FALSE)

cat("Loaded.\n\n")

# ---- Subset to common genes ----

cat("=== SUBSETTING TO COMMON GENES ===\n\n")

# Intersect common genes with actual matrix rownames
common_221561 <- intersect(common_genes, rownames(logcpm_221561))
common_197677 <- intersect(common_genes, rownames(logcpm_197677))
common_160269 <- intersect(common_genes, rownames(logcpm_160269))

# Use only genes present in all three
common_all <- Reduce(intersect, list(common_221561, common_197677, common_160269))
cat(sprintf("Common genes in all three datasets: %d\n", length(common_all)))

# Subset logCPM matrices
logcpm_221561_common <- logcpm_221561[common_all, ]
logcpm_197677_common <- logcpm_197677[common_all, ]
logcpm_160269_common <- logcpm_160269[common_all, ]

# ---- Subset to primary eligible samples ----

cat("\n=== SUBSETTING TO PRIMARY SAMPLES ===\n\n")

# GSE160269: exclude P63T and P87T
primary_160269 <- meta_160269$sample[meta_160269$epithelial_cells >= 30]
logcpm_160269_primary <- logcpm_160269_common[, primary_160269]

cat(sprintf("GSE221561: %d samples\n", ncol(logcpm_221561_common)))
cat(sprintf("GSE197677: %d samples\n", ncol(logcpm_197677_common)))
cat(sprintf("GSE160269: %d samples (primary)\n", ncol(logcpm_160269_primary)))
cat(sprintf("Total: %d samples\n", ncol(logcpm_221561_common) + ncol(logcpm_197677_common) + ncol(logcpm_160269_primary)))

# ---- Ensure gene order identical ----

cat("\n=== ALIGNING GENE ORDER ===\n\n")

# All should already have same genes, but verify
stopifnot(all(rownames(logcpm_221561_common) == rownames(logcpm_197677_common)))
stopifnot(all(rownames(logcpm_221561_common) == rownames(logcpm_160269_primary)))
cat("Gene order verified.\n")

# ---- Combine ----

cat("\n=== COMBINING MATRICES ===\n\n")

combined <- cbind(logcpm_221561_common, logcpm_197677_common, logcpm_160269_primary)
cat(sprintf("Combined matrix: %d genes x %d samples\n", nrow(combined), ncol(combined)))

# Build metadata
samples_221561 <- data.frame(
  sample = colnames(logcpm_221561_common),
  dataset = "GSE221561",
  treatment = meta_221561$treatment_standardized[match(colnames(logcpm_221561_common), meta_221561$sample)],
  stringsAsFactors = FALSE
)

samples_197677 <- data.frame(
  sample = colnames(logcpm_197677_common),
  dataset = "GSE197677",
  treatment = meta_197677$treatment_standardized[match(colnames(logcpm_197677_common), meta_197677$sample)],
  stringsAsFactors = FALSE
)

samples_160269 <- data.frame(
  sample = colnames(logcpm_160269_primary),
  dataset = "GSE160269",
  treatment = "Untreated_baseline",
  stringsAsFactors = FALSE
)

sample_info <- rbind(samples_221561, samples_197677, samples_160269)

# ---- Select top variable genes ----

cat("\n=== SELECTING VARIABLE GENES ===\n\n")

gene_vars <- apply(combined, 1, var)
n_select <- min(2000, nrow(combined))
top_genes <- names(sort(gene_vars, decreasing = TRUE))[1:n_select]
combined_top <- combined[top_genes, ]
cat(sprintf("Genes for PCA: %d (top %d by variance)\n", nrow(combined_top), n_select))

# ---- PCA ----

cat("\n=== RUNNING PCA ===\n\n")

pca_res <- prcomp(t(combined_top), center = TRUE, scale. = FALSE)
var_explained <- summary(pca_res)$importance[2, ] * 100

cat(sprintf("PC1 variance: %.1f%%\n", var_explained[1]))
cat(sprintf("PC2 variance: %.1f%%\n", var_explained[2]))

# ---- Plot ----

cat("\n=== GENERATING PLOTS ===\n\n")

# Determine dataset separation strength
dataset_labels <- sample_info$dataset
pc1_vals <- pca_res$x[, 1]

# Check if datasets separate strongly
dataset_means <- tapply(pc1_vals, dataset_labels, mean)
dataset_range <- diff(range(dataset_means))
within_range <- max(tapply(pc1_vals, dataset_labels, function(x) diff(range(x))))
separation_ratio <- dataset_range / within_range

cat(sprintf("Dataset separation ratio: %.2f\n", separation_ratio))
if (separation_ratio > 2) {
  structure_strength <- "STRONG"
} else if (separation_ratio > 1) {
  structure_strength <- "MODERATE"
} else {
  structure_strength <- "WEAK"
}
cat(sprintf("Dataset structure: %s\n", structure_strength))

# Color palette
dataset_colors <- c("GSE221561" = "#E41A1C", "GSE197677" = "#377EB8", "GSE160269" = "#4DAF4A")
treatment_pch <- c("Neoadjuvant_treated" = 19, "No_neoadjuvant_chemotherapy" = 17,
                    "Surgery_alone" = 15, "Untreated_baseline" = 18)

# Uncorrected PCA
png(file.path(fig_dir, "STEP19B_cross_dataset_UNCORRECTED_PCA.png"),
    width = 9, height = 6, units = "in", res = 150)

plot(pca_res$x[, 1], pca_res$x[, 2],
     col = dataset_colors[dataset_labels],
     pch = treatment_pch[sample_info$treatment],
     cex = 1.5,
     xlab = sprintf("PC1 (%.1f%% variance)", var_explained[1]),
     ylab = sprintf("PC2 (%.1f%% variance)", var_explained[2]),
     main = "Uncorrected Cross-Dataset PCA - Diagnostic Only",
     cex.main = 1.2)

# Legend
legend("topright",
       legend = c(names(dataset_colors), "", names(treatment_pch)),
       col = c(dataset_colors, NA, rep("black", length(treatment_pch))),
       pch = c(rep(19, 3), NA, treatment_pch),
       cex = 0.7,
       title = "Dataset / Treatment")

dev.off()
cat("Saved: STEP19B_cross_dataset_UNCORRECTED_PCA.png\n")

# Labeled version
png(file.path(fig_dir, "STEP19B_cross_dataset_UNCORRECTED_PCA_labeled.png"),
    width = 12, height = 8, units = "in", res = 150)

plot(pca_res$x[, 1], pca_res$x[, 2],
     col = dataset_colors[dataset_labels],
     pch = treatment_pch[sample_info$treatment],
     cex = 1.5,
     xlab = sprintf("PC1 (%.1f%% variance)", var_explained[1]),
     ylab = sprintf("PC2 (%.1f%% variance)", var_explained[2]),
     main = "Uncorrected Cross-Dataset PCA - Diagnostic Only",
     cex.main = 1.2)

text(pca_res$x[, 1], pca_res$x[, 2],
     labels = sample_info$sample, cex = 0.5, pos = 3, offset = 0.3)

legend("topright",
       legend = c(names(dataset_colors), "", names(treatment_pch)),
       col = c(dataset_colors, NA, rep("black", length(treatment_pch))),
       pch = c(rep(19, 3), NA, treatment_pch),
       cex = 0.7,
       title = "Dataset / Treatment")

dev.off()
cat("Saved: STEP19B_cross_dataset_UNCORRECTED_PCA_labeled.png\n")

cat("\nDone.\n")
