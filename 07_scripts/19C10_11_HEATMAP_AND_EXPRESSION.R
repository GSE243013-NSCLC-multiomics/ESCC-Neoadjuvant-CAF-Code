#!/usr/bin/env Rscript
# =============================================================================
# STEPS 19C10-19C11: TOP GENE HEATMAP AND SAMPLE-LEVEL EXPRESSION
# =============================================================================

suppressPackageStartupMessages({
  library(edgeR)
})

base_dir <- "~/Downloads/ESCC Neoadjuvant"
de_dir <- file.path(base_dir, "06_tables/cross_dataset/epithelial_pseudobulk/GSE197677_DE")
fig_dir <- file.path(base_dir, "05_figures/cross_dataset/epithelial_pseudobulk/GSE197677_DE")

cat("Loading data...\n")

dge_f <- readRDS(file.path(de_dir, "GSE197677_DGE_filtered_TMM_dispersion.rds"))
res <- read.csv(file.path(de_dir, "STEP19C_GSE197677_NACT_vs_nNACT_ALL_GENES.csv"), row.names = 1, check.names = FALSE)
meta <- read.csv(file.path(base_dir, "06_tables/cross_dataset/epithelial_pseudobulk/counts/GSE197677_epithelial_pseudobulk_sample_metadata.csv"), stringsAsFactors = FALSE)

cat("Loaded.\n\n")

# ---- Calculate logCPM for visualization ----

logcpm <- cpm(dge_f, log = TRUE, prior.count = 2)

# ---- STEP 19C10: TOP-GENE HEATMAP ----

cat("=== STEP 19C10: TOP-GENE HEATMAP ===\n\n")

# Select top 30 genes by FDR
n_genes <- min(30, nrow(res))
top_genes <- rownames(res)[1:n_genes]

# Check if all FDR < 0.05
all_sig <- all(res$FDR[1:n_genes] < 0.05)
if (!all_sig) {
  cat("NOTE: Not all top 30 genes are FDR < 0.05. Figure labeled accordingly.\n\n")
}

# Subset logCPM to top genes
logcpm_top <- logcpm[top_genes, ]

# Z-score each gene across samples (visualization only)
logcpm_z <- t(scale(t(logcpm_top)))

# Order samples by treatment
sample_order <- meta$sample[order(meta$treatment_standardized)]
logcpm_z <- logcpm_z[, sample_order]

# Treatment annotation
treatment_ann <- meta$treatment_standardized[match(sample_order, meta$sample)]
treatment_colors <- ifelse(treatment_ann == "Neoadjuvant_chemotherapy", "#E41A1C", "#377EB8")

# Plot
png(file.path(fig_dir, "STEP19C_GSE197677_top30_DE_heatmap.png"),
    width = 10, height = 8, units = "in", res = 150)

# Use base R heatmap-compatible image
# Transpose so samples are rows, genes are columns for standard heatmap layout
heatmap_mat <- t(logcpm_z)  # samples x genes

# Simple heatmap with annotations
par(mar = c(8, 6, 4, 2))
image(t(heatmap_mat),
      col = colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100),
      axes = FALSE,
      main = sprintf("Top %d Genes by FDR (Z-scored)", n_genes),
      cex.main = 1.1)

# X-axis: gene names
axis(1, at = seq(0, 1, length.out = n_genes), labels = rownames(logcpm_z),
     las = 2, cex.axis = 0.6, tick = FALSE)

# Y-axis: sample names with treatment color
axis(2, at = seq(0, 1, length.out = length(sample_order)),
     labels = paste0(sample_order, " (",
                     ifelse(treatment_ann == "Neoadjuvant_chemotherapy", "NACT", "nNACT"), ")"),
     las = 2, cex.axis = 0.7, tick = FALSE)

# Legend
legend("topright", legend = c("NACT", "nNACT"), col = c("#E41A1C", "#377EB8"),
       pch = 15, cex = 0.8, title = "Treatment")

if (!all_sig) {
  mtext("Top-ranked genes - not all statistically significant", side = 3, line = -0.5, cex = 0.8, col = "grey40")
}

dev.off()
cat("Saved: STEP19C_GSE197677_top30_DE_heatmap.png\n\n")

# ---- STEP 19C11: SAMPLE-LEVEL EXPRESSION CHECK ----

cat("=== STEP 19C11: SAMPLE-LEVEL EXPRESSION CHECK ===\n\n")

# Top 20 genes
n_top20 <- min(20, nrow(res))
top20_genes <- rownames(res)[1:n_top20]

# Build long-format table
long_expr <- data.frame()
for (gene in top20_genes) {
  for (s in meta$sample) {
    long_expr <- rbind(long_expr, data.frame(
      gene = gene,
      sample = s,
      treatment = meta$treatment_standardized[meta$sample == s],
      logCPM = logcpm[gene, s],
      stringsAsFactors = FALSE
    ))
  }
}

write.csv(long_expr, file.path(de_dir, "STEP19C_GSE197677_top20_gene_sample_expression.csv"), row.names = FALSE)
cat("Saved: STEP19C_GSE197677_top20_gene_sample_expression.csv\n\n")

# Faceted plot for top 12 genes
n_plot <- min(12, length(top20_genes))
plot_genes <- top20_genes[1:n_plot]

png(file.path(fig_dir, "STEP19C_GSE197677_top12_gene_expression.png"),
    width = 12, height = 10, units = "in", res = 150)

par(mfrow = c(4, 3), mar = c(3, 3, 2, 1))

for (gene in plot_genes) {
  gene_data <- long_expr[long_expr$gene == gene, ]
  
  # Box plot
  box_data_nact <- gene_data$logCPM[gene_data$treatment == "Neoadjuvant_chemotherapy"]
  box_data_nnact <- gene_data$logCPM[gene_data$treatment == "No_neoadjuvant_chemotherapy"]
  
  boxplot(list(NACT = box_data_nact, nNACT = box_data_nnact),
          main = gene, cex.main = 0.9,
          col = c("#E41A1C", "#377EB8"),
          ylab = "logCPM", cex.axis = 0.8)
  
  # Add individual points
  points(jitter(rep(1, length(box_data_nact))), box_data_nact, pch = 16, cex = 0.8, col = "#E41A1C")
  points(jitter(rep(2, length(box_data_nnact))), box_data_nnact, pch = 16, cex = 0.8, col = "#377EB8")
}

dev.off()
cat("Saved: STEP19C_GSE197677_top12_gene_expression.png\n\n")

cat("=== STEPS 19C10-19C11 COMPLETE ===\n")
