#!/usr/bin/env Rscript
# =============================================================================
# STEP 19B11: DATASET-CENTERED DIAGNOSTIC PCA
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

# ---- Subset to common genes and primary samples ----

cat("=== SUBSETTING TO COMMON GENES AND PRIMARY SAMPLES ===\n\n")

common_all <- Reduce(intersect, list(
  intersect(common_genes, rownames(logcpm_221561)),
  intersect(common_genes, rownames(logcpm_197677)),
  intersect(common_genes, rownames(logcpm_160269))
))
cat(sprintf("Common genes: %d\n", length(common_all)))

primary_160269 <- meta_160269$sample[meta_160269$epithelial_cells >= 30]

logcpm_221561_common <- logcpm_221561[common_all, ]
logcpm_197677_common <- logcpm_197677[common_all, ]
logcpm_160269_common <- logcpm_160269[common_all, primary_160269]

cat(sprintf("Samples: GSE221561=%d, GSE197677=%d, GSE160269=%d\n",
            ncol(logcpm_221561_common), ncol(logcpm_197677_common), ncol(logcpm_160269_common)))

# ---- Dataset centering ----

cat("\n=== DATASET CENTERING ===\n\n")
cat("For each gene within each dataset: subtract dataset-specific mean.\n")
cat("Do NOT divide by SD.\n\n")

center_dataset <- function(mat) {
  # Subtract row means (gene-wise centering within dataset)
  mat_centered <- mat - rowMeans(mat)
  mat_centered
}

centered_221561 <- center_dataset(logcpm_221561_common)
centered_197677 <- center_dataset(logcpm_197677_common)
centered_160269 <- center_dataset(logcpm_160269_common)

cat("Centering complete.\n")

# ---- Combine ----

combined_centered <- cbind(centered_221561, centered_197677, centered_160269)
cat(sprintf("Combined centered matrix: %d genes x %d samples\n", nrow(combined_centered), ncol(combined_centered)))

# Also build uncorrected combined for coordinate comparison
combined_uncorrected <- cbind(logcpm_221561_common, logcpm_197677_common, logcpm_160269_common)

# Sample info
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
  sample = colnames(logcpm_160269_common),
  dataset = "GSE160269",
  treatment = "Untreated_baseline",
  stringsAsFactors = FALSE
)

sample_info <- rbind(samples_221561, samples_197677, samples_160269)

# ---- Select top variable genes ----

cat("\n=== SELECTING VARIABLE GENES ===\n\n")

gene_vars <- apply(combined_centered, 1, var)
n_select <- min(2000, nrow(combined_centered))
top_genes <- names(sort(gene_vars, decreasing = TRUE))[1:n_select]
combined_top <- combined_centered[top_genes, ]
cat(sprintf("Genes for PCA: %d\n", nrow(combined_top)))

# ---- PCA on centered data ----

cat("\n=== RUNNING PCA (CENTERED) ===\n\n")

pca_centered <- prcomp(t(combined_top), center = TRUE, scale. = FALSE)
var_centered <- summary(pca_centered)$importance[2, ] * 100

cat(sprintf("PC1 variance: %.1f%%\n", var_centered[1]))
cat(sprintf("PC2 variance: %.1f%%\n", var_centered[2]))

# ---- PCA on uncorrected data for comparison ----

cat("\n=== RUNNING PCA (UNCORRECTED) ===\n\n")

uncorrected_top <- combined_uncorrected[top_genes, ]
pca_uncorrected <- prcomp(t(uncorrected_top), center = TRUE, scale. = FALSE)
var_uncorrected <- summary(pca_uncorrected)$importance[2, ] * 100

cat(sprintf("PC1 variance: %.1f%%\n", var_uncorrected[1]))
cat(sprintf("PC2 variance: %.1f%%\n", var_uncorrected[2]))

# ---- Plot ----

cat("\n=== GENERATING PLOTS ===\n\n")

dataset_colors <- c("GSE221561" = "#E41A1C", "GSE197677" = "#377EB8", "GSE160269" = "#4DAF4A")
treatment_pch <- c("Neoadjuvant_treated" = 19, "No_neoadjuvant_chemotherapy" = 17,
                    "Surgery_alone" = 15, "Untreated_baseline" = 18)

dataset_labels <- sample_info$dataset

png(file.path(fig_dir, "STEP19B_cross_dataset_DATASET_CENTERED_PCA.png"),
    width = 9, height = 6, units = "in", res = 150)

plot(pca_centered$x[, 1], pca_centered$x[, 2],
     col = dataset_colors[dataset_labels],
     pch = treatment_pch[sample_info$treatment],
     cex = 1.5,
     xlab = sprintf("PC1 (%.1f%% variance)", var_centered[1]),
     ylab = sprintf("PC2 (%.1f%% variance)", var_centered[2]),
     main = "Dataset-Centered PCA - Visualization Only",
     cex.main = 1.2)

legend("topright",
       legend = c(names(dataset_colors), "", names(treatment_pch)),
       col = c(dataset_colors, NA, rep("black", length(treatment_pch))),
       pch = c(rep(19, 3), NA, treatment_pch),
       cex = 0.7,
       title = "Dataset / Treatment")

dev.off()
cat("Saved: STEP19B_cross_dataset_DATASET_CENTERED_PCA.png\n")

# ---- Save coordinates ----

cat("\n=== SAVING COORDINATES ===\n\n")

coords <- data.frame(
  sample = sample_info$sample,
  dataset = sample_info$dataset,
  treatment = sample_info$treatment,
  PC1_uncorrected = pca_uncorrected$x[, 1],
  PC2_uncorrected = pca_uncorrected$x[, 2],
  PC1_centered = pca_centered$x[, 1],
  PC2_centered = pca_centered$x[, 2],
  stringsAsFactors = FALSE
)

write.csv(coords, file.path(qc_dir, "STEP19B_cross_dataset_PCA_coordinates.csv"), row.names = FALSE)
cat("Saved: STEP19B_cross_dataset_PCA_coordinates.csv\n")

# ---- Assess residual structure ----

pc1_centered <- pca_centered$x[, 1]
dataset_means_centered <- tapply(pc1_centered, dataset_labels, mean)
dataset_range_centered <- diff(range(dataset_means_centered))
within_range_centered <- max(tapply(pc1_centered, dataset_labels, function(x) diff(range(x))))
separation_ratio_centered <- dataset_range_centered / within_range_centered

cat(sprintf("\nDataset-centered separation ratio: %.2f\n", separation_ratio_centered))
if (separation_ratio_centered > 2) {
  residual_strength <- "STRONG"
} else if (separation_ratio_centered > 1) {
  residual_strength <- "MODERATE"
} else {
  residual_strength <- "WEAK"
}
cat(sprintf("Residual structure: %s\n", residual_strength))

cat("\nDone.\n")
