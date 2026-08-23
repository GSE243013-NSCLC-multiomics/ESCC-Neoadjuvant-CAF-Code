#!/usr/bin/env Rscript
# =============================================================================
# STEP 19B2: EXPRESSION FILTER AUDIT
# =============================================================================

suppressPackageStartupMessages({
  library(edgeR)
})

base_dir <- "~/Downloads/ESCC Neoadjuvant"
qc_dir <- file.path(base_dir, "06_tables/cross_dataset/epithelial_pseudobulk/QC")

cat("Loading DGEList objects and design matrices...\n")

dge_221561 <- readRDS(file.path(qc_dir, "GSE221561_DGEList_TMM.rds"))
dge_197677 <- readRDS(file.path(qc_dir, "GSE197677_DGEList_TMM.rds"))
dge_160269 <- readRDS(file.path(qc_dir, "GSE160269_DGEList_TMM.rds"))

# These are metadata CSVs; build proper model matrices for filterByExpr
meta_design_221561 <- read.csv(file.path(base_dir, "06_tables/cross_dataset/epithelial_pseudobulk/STEP19A_GSE221561_design_matrix.csv"), row.names = 1, check.names = FALSE)
meta_design_197677 <- read.csv(file.path(base_dir, "06_tables/cross_dataset/epithelial_pseudobulk/STEP19A_GSE197677_design_matrix.csv"), row.names = 1, check.names = FALSE)

# Build proper numeric design matrices
design_221561 <- model.matrix(~ treatment_standardized, data = meta_design_221561)
rownames(design_221561) <- rownames(meta_design_221561)

design_197677 <- model.matrix(~ treatment_standardized, data = meta_design_197677)
rownames(design_197677) <- rownames(meta_design_197677)

cat("Design matrices built.\n")

cat("Loaded.\n\n")

# ---- GSE197677 DE candidate filter ----

cat("=== GSE197677: DE_CANDIDATE_FILTER ===\n\n")

# filterByExpr uses the design matrix to determine minimum group size
# For pseudobulk with design: min.sample.size = min(colSums(design))
keep_197677 <- filterByExpr(dge_197677, design = design_197677)
genes_total_197677 <- nrow(dge_197677)
genes_keep_197677 <- sum(keep_197677)
genes_removed_197677 <- genes_total_197677 - genes_keep_197677

cat(sprintf("  genes_total:    %d\n", genes_total_197677))
cat(sprintf("  genes_keep:     %d\n", genes_keep_197677))
cat(sprintf("  genes_removed:  %d\n", genes_removed_197677))
cat(sprintf("  pct_keep:       %.1f%%\n", 100 * genes_keep_197677 / genes_total_197677))

cat("\n")

# ---- GSE221561 DE candidate filter ----

cat("=== GSE221561: DE_CANDIDATE_FILTER ===\n\n")

keep_221561 <- filterByExpr(dge_221561, design = design_221561)
genes_total_221561 <- nrow(dge_221561)
genes_keep_221561 <- sum(keep_221561)
genes_removed_221561 <- genes_total_221561 - genes_keep_221561

cat(sprintf("  genes_total:    %d\n", genes_total_221561))
cat(sprintf("  genes_keep:     %d\n", genes_keep_221561))
cat(sprintf("  genes_removed:  %d\n", genes_removed_221561))
cat(sprintf("  pct_keep:       %.1f%%\n", 100 * genes_keep_221561 / genes_total_221561))

cat("\n")

# ---- GSE160269 QC expression filter ----

cat("=== GSE160269: QC_EXPRESSION_FILTER ===\n\n")
cat("  No treatment contrast. Using CPM >= 1 in >= 6 of 58 primary samples.\n\n")

# Identify primary samples (epithelial_cells >= 30)
primary_meta <- dge_160269$samples[dge_160269$samples$epithelial_cells >= 30, ]
primary_samples <- primary_meta$sample
cat(sprintf("  Primary samples: %d\n", length(primary_samples)))

# Calculate CPM on TMM-normalized data
logcpm_full <- cpm(dge_160269, log = TRUE, prior.count = 2)
# Subset to primary samples
logcpm_primary <- logcpm_full[, primary_samples]

# CPM >= 1 means logCPM >= 0
keep_160269 <- rowSums(logcpm_primary >= 0) >= 6  # at least 6 of 58 primary samples

genes_total_160269 <- nrow(dge_160269)
genes_keep_160269 <- sum(keep_160269)
genes_removed_160269 <- genes_total_160269 - genes_keep_160269

cat(sprintf("  genes_total:    %d\n", genes_total_160269))
cat(sprintf("  genes_keep:     %d\n", genes_keep_160269))
cat(sprintf("  genes_removed:  %d\n", genes_removed_160269))
cat(sprintf("  pct_keep:       %.1f%%\n", 100 * genes_keep_160269 / genes_total_160269))
cat("  NOTE: This filter is for QC visualization only, NOT for DE.\n")

cat("\n")

# ---- Save filter audit ----

cat("=== SAVING FILTER AUDIT ===\n\n")

filter_audit <- data.frame(
  dataset = c("GSE197677", "GSE221561", "GSE160269"),
  filter_type = c("DE_CANDIDATE_FILTER", "DE_CANDIDATE_FILTER", "QC_EXPRESSION_FILTER"),
  genes_total = c(genes_total_197677, genes_total_221561, genes_total_160269),
  genes_keep = c(genes_keep_197677, genes_keep_221561, genes_keep_160269),
  genes_removed = c(genes_removed_197677, genes_removed_221561, genes_removed_160269),
  filter_definition = c(
    "edgeR::filterByExpr with design matrix (min.count=10, min.total.count=15, large.n=10, min.prop=0.7)",
    "edgeR::filterByExpr with design matrix (min.count=10, min.total.count=15, large.n=10, min.prop=0.7)",
    "CPM >= 1 in at least 6 of 58 primary samples (~10%); QC visualization only"
  ),
  stringsAsFactors = FALSE
)

write.csv(filter_audit, file.path(qc_dir, "STEP19B_expression_filter_audit.csv"), row.names = FALSE)
cat("Saved: STEP19B_expression_filter_audit.csv\n")

# ---- Save filter vectors for downstream steps ----

saveRDS(keep_197677, file.path(qc_dir, "GSE197677_filter_keep.rds"))
saveRDS(keep_221561, file.path(qc_dir, "GSE221561_filter_keep.rds"))
saveRDS(keep_160269, file.path(qc_dir, "GSE160269_filter_keep.rds"))
cat("Saved filter vectors for downstream use.\n")

cat("\nDone.\n")
