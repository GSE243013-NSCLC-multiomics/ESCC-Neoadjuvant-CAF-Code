#!/usr/bin/env Rscript
# =============================================================================
# STEP 19B1: BUILD DGEList OBJECTS AND TMM NORMALIZATION
# =============================================================================

suppressPackageStartupMessages({
  library(edgeR)
})

base_dir <- "~/Downloads/ESCC Neoadjuvant"
qc_dir <- file.path(base_dir, "06_tables/cross_dataset/epithelial_pseudobulk/QC")

cat("Loading pseudobulk counts and metadata...\n")

counts_221561 <- readRDS(file.path(base_dir, "06_tables/cross_dataset/epithelial_pseudobulk/counts/GSE221561_epithelial_pseudobulk_raw_counts.rds"))
counts_197677 <- readRDS(file.path(base_dir, "06_tables/cross_dataset/epithelial_pseudobulk/counts/GSE197677_epithelial_pseudobulk_raw_counts.rds"))
counts_160269 <- readRDS(file.path(base_dir, "06_tables/cross_dataset/epithelial_pseudobulk/counts/GSE160269_epithelial_pseudobulk_raw_counts.rds"))

meta_221561 <- read.csv(file.path(base_dir, "06_tables/cross_dataset/epithelial_pseudobulk/counts/GSE221561_epithelial_pseudobulk_sample_metadata.csv"), stringsAsFactors = FALSE)
meta_197677 <- read.csv(file.path(base_dir, "06_tables/cross_dataset/epithelial_pseudobulk/counts/GSE197677_epithelial_pseudobulk_sample_metadata.csv"), stringsAsFactors = FALSE)
meta_160269 <- read.csv(file.path(base_dir, "06_tables/cross_dataset/epithelial_pseudobulk/counts/GSE160269_epithelial_pseudobulk_sample_metadata.csv"), stringsAsFactors = FALSE)

cat("Loaded.\n\n")

# ---- Build DGEList objects ----

cat("=== BUILDING DGEList OBJECTS ===\n\n")

dge_221561 <- DGEList(counts = counts_221561, samples = meta_221561, genes = data.frame(gene = rownames(counts_221561), stringsAsFactors = FALSE))
cat("GSE221561 DGEList: ", nrow(dge_221561$counts), " genes x ", ncol(dge_221561$counts), " samples\n")

dge_197677 <- DGEList(counts = counts_197677, samples = meta_197677, genes = data.frame(gene = rownames(counts_197677), stringsAsFactors = FALSE))
cat("GSE197677 DGEList: ", nrow(dge_197677$counts), " genes x ", ncol(dge_197677$counts), " samples\n")

dge_160269 <- DGEList(counts = counts_160269, samples = meta_160269, genes = data.frame(gene = rownames(counts_160269), stringsAsFactors = FALSE))
cat("GSE160269 DGEList: ", nrow(dge_160269$counts), " genes x ", ncol(dge_160269$counts), " samples\n")

cat("\n")

# ---- TMM normalization (each dataset independently) ----

cat("=== TMM NORMALIZATION ===\n\n")

dge_221561 <- calcNormFactors(dge_221561, method = "TMM")
cat("GSE221561 TMM norm factors: range = [", min(dge_221561$samples$norm.factors), ", ", max(dge_221561$samples$norm.factors), "]\n")

dge_197677 <- calcNormFactors(dge_197677, method = "TMM")
cat("GSE197677 TMM norm factors: range = [", min(dge_197677$samples$norm.factors), ", ", max(dge_197677$samples$norm.factors), "]\n")

dge_160269 <- calcNormFactors(dge_160269, method = "TMM")
cat("GSE160269 TMM norm factors: range = [", min(dge_160269$samples$norm.factors), ", ", max(dge_160269$samples$norm.factors), "]\n")

cat("\n")

# ---- Save DGEList objects ----

cat("=== SAVING DGEList OBJECTS ===\n\n")

saveRDS(dge_221561, file.path(qc_dir, "GSE221561_DGEList_TMM.rds"))
cat("Saved: GSE221561_DGEList_TMM.rds\n")

saveRDS(dge_197677, file.path(qc_dir, "GSE197677_DGEList_TMM.rds"))
cat("Saved: GSE197677_DGEList_TMM.rds\n")

saveRDS(dge_160269, file.path(qc_dir, "GSE160269_DGEList_TMM.rds"))
cat("Saved: GSE160269_DGEList_TMM.rds\n")

cat("\n")

# ---- Build normalization factor table ----

cat("=== BUILDING NORMALIZATION FACTOR TABLE ===\n\n")

build_norm_table <- function(dge, dataset) {
  data.frame(
    dataset = dataset,
    sample = dge$samples$sample,
    library_size = dge$samples$lib.size,
    norm_factor = dge$samples$norm.factors,
    effective_library_size = dge$samples$lib.size * dge$samples$norm.factors,
    epithelial_cells = dge$samples$epithelial_cells,
    treatment = dge$samples$treatment_standardized,
    stringsAsFactors = FALSE
  )
}

norm_221561 <- build_norm_table(dge_221561, "GSE221561")
norm_197677 <- build_norm_table(dge_197677, "GSE197677")
norm_160269 <- build_norm_table(dge_160269, "GSE160269")

norm_all <- rbind(norm_221561, norm_197677, norm_160269)

write.csv(norm_all, file.path(qc_dir, "STEP19B_TMM_normalization_factors.csv"), row.names = FALSE)
cat("Saved: STEP19B_TMM_normalization_factors.csv\n")

# Print summary
cat("\n=== NORMALIZATION SUMMARY ===\n\n")
for (ds in c("GSE221561", "GSE197677", "GSE160269")) {
  sub <- norm_all[norm_all$dataset == ds, ]
  cat(sprintf("%s: library_size range = [%s, %s], norm_factor range = [%s, %s]\n",
              ds,
              format(min(sub$library_size), big.mark = ","),
              format(max(sub$library_size), big.mark = ","),
              round(min(sub$norm_factor), 4),
              round(max(sub$norm_factor), 4)))
}

cat("\nDone.\n")
