#!/usr/bin/env Rscript
# =============================================================================
# STEP 19E: CROSS-DATASET HALLMARK PRERANKED PATHWAY ANALYSIS
# =============================================================================

suppressPackageStartupMessages({
  library(fgsea)
  library(msigdbr)
  library(data.table)
  library(ggplot2)
})

base_dir <- "~/Downloads/ESCC Neoadjuvant"
de_dir <- file.path(base_dir, "06_tables/cross_dataset/epithelial_pseudobulk/GSE197677_DE")
repl_dir <- file.path(base_dir, "06_tables/cross_dataset/epithelial_pseudobulk/replication")
pw_dir <- file.path(base_dir, "06_tables/cross_dataset/pathway_analysis")
fig_dir <- file.path(base_dir, "05_figures/cross_dataset/pathway_analysis/Hallmark")
dir.create(pw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 19E0: PACKAGE + HALLMARK RESOURCE PREFLIGHT
# =============================================================================

cat("=== 19E0: PACKAGE + HALLMARK RESOURCE PREFLIGHT ===\n\n")

cat("Package versions:\n")
cat(sprintf("  fgsea:      %s\n", as.character(packageVersion("fgsea"))))
cat(sprintf("  msigdbr:    %s\n", as.character(packageVersion("msigdbr"))))
cat(sprintf("  data.table: %s\n\n", as.character(packageVersion("data.table"))))

writeLines(c(
  sprintf("fgsea: %s", as.character(packageVersion("fgsea"))),
  sprintf("msigdbr: %s", as.character(packageVersion("msigdbr"))),
  sprintf("data.table: %s", as.character(packageVersion("data.table"))),
  sprintf("R: %s", R.version.string),
  sprintf("Date: %s", Sys.Date())
), file.path(pw_dir, "STEP19E_PACKAGE_VERSIONS.txt"))

cat("Loading Hallmark gene sets from msigdbr cache...\n")

cache_file <- file.path(Sys.getenv("ESCC_CAF_PROJECT_ROOT", unset = "."), "cache", "msigdb.2026.1.Hs.H.rds")
hallmark_df <- readRDS(cache_file)

cat(sprintf("Loaded from cache: %d rows\n", nrow(hallmark_df)))

pathway_col <- "gs_name"
gene_col <- "db_gene_symbol"
hallmark_pathways <- split(hallmark_df[[gene_col]], hallmark_df[[pathway_col]])
hallmark_pathways <- lapply(hallmark_pathways, unique)

cat(sprintf("Hallmark pathways retrieved: %d\n", length(hallmark_pathways)))

hallmark_df_frozen <- data.frame(
  pathway = hallmark_df[[pathway_col]],
  gene_symbol = hallmark_df[[gene_col]],
  stringsAsFactors = FALSE
)
utils::write.csv(hallmark_df_frozen, file.path(pw_dir, "STEP19E_HALLMARK_GENESETS_FROZEN.csv"), row.names = FALSE)
cat("Saved: STEP19E_HALLMARK_GENESETS_FROZEN.csv\n\n")

# =============================================================================
# 19E1: BUILD COMMON GENE UNIVERSE
# =============================================================================

cat("=== 19E1: BUILD COMMON GENE UNIVERSE ===\n\n")

disc_all <- read.csv(file.path(de_dir, "STEP19C_GSE197677_NACT_vs_nNACT_ALL_GENES.csv"), row.names = 1, check.names = FALSE)
repl_all <- read.csv(file.path(repl_dir, "STEP19D_GSE221561_NEOADJUVANT_MINUS_SURGERY_CANONICAL.csv"), row.names = 1, check.names = FALSE)

cat(sprintf("GSE197677 genes: %d\n", nrow(disc_all)))
cat(sprintf("GSE221561 genes: %d\n", nrow(repl_all)))

common_genes <- intersect(rownames(disc_all), rownames(repl_all))
cat(sprintf("Common genes: %d\n", length(common_genes)))

if (abs(length(common_genes) - 9033) > 50) {
  cat("WARNING: Common gene count differs substantially from 9033\n")
}

common_universe <- data.frame(
  gene = common_genes,
  GSE197677_logFC = disc_all[common_genes, "logFC"],
  GSE197677_F = disc_all[common_genes, "F"],
  GSE197677_PValue = disc_all[common_genes, "PValue"],
  GSE197677_FDR = disc_all[common_genes, "FDR"],
  GSE221561_logFC = repl_all[common_genes, "logFC"],
  GSE221561_F = repl_all[common_genes, "F"],
  GSE221561_PValue = repl_all[common_genes, "PValue"],
  GSE221561_FDR = repl_all[common_genes, "FDR"],
  stringsAsFactors = FALSE
)

utils::write.csv(common_universe, file.path(pw_dir, "STEP19E_COMMON_GENE_UNIVERSE.csv"), row.names = FALSE)
cat("Saved: STEP19E_COMMON_GENE_UNIVERSE.csv\n\n")

# =============================================================================
# 19E2: BUILD PRERANKED STATISTICS
# =============================================================================

cat("=== 19E2: BUILD PRERANKED STATISTICS ===\n\n")

disc_rank <- sign(disc_all[common_genes, "logFC"]) * sqrt(disc_all[common_genes, "F"])
names(disc_rank) <- common_genes

repl_rank <- sign(repl_all[common_genes, "logFC"]) * sqrt(repl_all[common_genes, "F"])
names(repl_rank) <- common_genes

disc_valid <- !is.na(disc_rank) & is.finite(disc_rank)
repl_valid <- !is.na(repl_rank) & is.finite(repl_rank)

cat(sprintf("GSE197677: %d valid / %d total\n", sum(disc_valid), length(disc_rank)))
cat(sprintf("GSE221561: %d valid / %d total\n\n", sum(repl_valid), length(repl_rank)))

cat("GSE197677 sign distribution:\n")
cat(sprintf("  Positive (Higher in NACT): %d\n", sum(disc_rank[disc_valid] > 0)))
cat(sprintf("  Negative (Higher in nNACT): %d\n", sum(disc_rank[disc_valid] < 0)))
cat(sprintf("  Zero: %d\n\n", sum(disc_rank[disc_valid] == 0)))

cat("GSE221561 sign distribution:\n")
cat(sprintf("  Positive (Higher in Neoadjuvant_treated): %d\n", sum(repl_rank[repl_valid] > 0)))
cat(sprintf("  Negative (Higher in Surgery_alone): %d\n", sum(repl_rank[repl_valid] < 0)))
cat(sprintf("  Zero: %d\n\n", sum(repl_rank[repl_valid] == 0)))

disc_rank_df <- data.frame(
  gene = names(disc_rank), rank_statistic = disc_rank,
  logFC = disc_all[names(disc_rank), "logFC"],
  F_stat = disc_all[names(disc_rank), "F"],
  PValue = disc_all[names(disc_rank), "PValue"],
  FDR = disc_all[names(disc_rank), "FDR"],
  stringsAsFactors = FALSE
)
utils::write.csv(disc_rank_df, file.path(pw_dir, "STEP19E_GSE197677_COMMON_GENE_RANKS.csv"), row.names = FALSE)

repl_rank_df <- data.frame(
  gene = names(repl_rank), rank_statistic = repl_rank,
  logFC = repl_all[names(repl_rank), "logFC"],
  F_stat = repl_all[names(repl_rank), "F"],
  PValue = repl_all[names(repl_rank), "PValue"],
  FDR = repl_all[names(repl_rank), "FDR"],
  stringsAsFactors = FALSE
)
utils::write.csv(repl_rank_df, file.path(pw_dir, "STEP19E_GSE221561_COMMON_GENE_RANKS.csv"), row.names = FALSE)

cat("Saved: STEP19E_GSE197677_COMMON_GENE_RANKS.csv\n")
cat("Saved: STEP19E_GSE221561_COMMON_GENE_RANKS.csv\n\n")

# =============================================================================
# 19E3: RANK QUALITY AUDIT
# =============================================================================

cat("=== 19E3: RANK QUALITY AUDIT ===\n\n")

for (ds in c("GSE197677", "GSE221561")) {
  r <- if (ds == "GSE197677") disc_rank[disc_valid] else repl_rank[repl_valid]
  cat(sprintf("%s:\n", ds))
  cat(sprintf("  Total: %d, Pos: %d, Neg: %d, Zero: %d, Dup: %d\n",
              length(r), sum(r > 0), sum(r < 0), sum(r == 0), sum(duplicated(names(r)))))
  cat(sprintf("  Min=%.2f, Q1=%.2f, Median=%.2f, Q3=%.2f, Max=%.2f\n\n",
              min(r), quantile(r, 0.25), median(r), quantile(r, 0.75), max(r)))
}

png(file.path(fig_dir, "STEP19E_GSE197677_rank_distribution.png"), width = 800, height = 600, res = 150)
hist(disc_rank[disc_valid], breaks = 50, col = "steelblue", main = "GSE197677 Rank Distribution",
     xlab = "Rank Statistic (sign(logFC) * sqrt(F))", ylab = "Frequency")
abline(v = 0, col = "red", lty = 2)
dev.off()

png(file.path(fig_dir, "STEP19E_GSE221561_rank_distribution.png"), width = 800, height = 600, res = 150)
hist(repl_rank[repl_valid], breaks = 50, col = "tomato", main = "GSE221561 Rank Distribution",
     xlab = "Rank Statistic (sign(logFC) * sqrt(F))", ylab = "Frequency")
abline(v = 0, col = "red", lty = 2)
dev.off()

cat("Saved rank distribution figures.\n\n")

# =============================================================================
# 19E4: HALLMARK FILTER
# =============================================================================

cat("=== 19E4: HALLMARK FILTER ===\n\n")

common_set <- sort(unique(common_genes))
hallmark_common <- lapply(hallmark_pathways, function(genes) intersect(genes, common_set))

minSize <- 15
maxSize <- 500

hallmark_filtered <- hallmark_common[sapply(hallmark_common, function(genes) {
  length(genes) >= minSize & length(genes) <= maxSize
})]

cat(sprintf("Hallmark pathways before filter: %d\n", length(hallmark_common)))
cat(sprintf("Hallmark pathways after filter (size %d-%d): %d\n\n", minSize, maxSize, length(hallmark_filtered)))

size_audit <- data.frame(
  pathway = names(hallmark_filtered),
  effective_size = sapply(hallmark_filtered, length),
  stringsAsFactors = FALSE
)
size_audit <- size_audit[order(-size_audit$effective_size), ]

utils::write.csv(size_audit, file.path(pw_dir, "STEP19E_HALLMARK_EFFECTIVE_SIZE_AUDIT.csv"), row.names = FALSE)
cat("Saved: STEP19E_HALLMARK_EFFECTIVE_SIZE_AUDIT.csv\n\n")

# =============================================================================
# 19E5: GSE197677 HALLMARK FGSEA
# =============================================================================

cat("=== 19E5: GSE197677 HALLMARK FGSEA ===\n\n")

disc_rank_fgsea <- disc_rank[disc_valid]
disc_rank_fgsea <- disc_rank_fgsea[!is.na(disc_rank_fgsea) & is.finite(disc_rank_fgsea)]

cat(sprintf("Running fgseaMultilevel: %d genes, %d pathways...\n", length(disc_rank_fgsea), length(hallmark_filtered)))

fgsea_disc <- fgseaMultilevel(
  pathways = hallmark_filtered, stats = disc_rank_fgsea,
  minSize = minSize, maxSize = maxSize, eps = 0
)

cat(sprintf("Completed. Pathways tested: %d\n\n", nrow(fgsea_disc)))

fgsea_disc$direction <- ifelse(fgsea_disc$NES > 0, "Higher_in_NACT", "Higher_in_nNACT")
fgsea_disc <- fgsea_disc[order(fgsea_disc$padj), ]

fgsea_disc_out <- data.frame(
  pathway = fgsea_disc$pathway, pval = fgsea_disc$pval, padj = fgsea_disc$padj,
  ES = fgsea_disc$ES, NES = fgsea_disc$NES, size = fgsea_disc$size,
  leadingEdge = sapply(fgsea_disc$leadingEdge, function(x) paste(x, collapse = ",")),
  direction = fgsea_disc$direction, stringsAsFactors = FALSE
)
utils::write.csv(fgsea_disc_out, file.path(pw_dir, "STEP19E_GSE197677_HALLMARK_FGSEA_ALL.csv"), row.names = FALSE)
cat("Saved: STEP19E_GSE197677_HALLMARK_FGSEA_ALL.csv\n")

fdr05_disc <- fgsea_disc[fgsea_disc$padj < 0.05, ]
fdr10_disc <- fgsea_disc[fgsea_disc$padj < 0.10, ]

cat(sprintf("FDR < 0.05: %d pathways\n", nrow(fdr05_disc)))
cat(sprintf("FDR < 0.10: %d pathways\n\n", nrow(fdr10_disc)))

if (nrow(fdr05_disc) > 0) {
  cat("FDR < 0.05 pathways:\n")
  for (i in 1:nrow(fdr05_disc)) {
    cat(sprintf("  %s: NES=%.3f, padj=%.4f, %s\n", fdr05_disc$pathway[i], fdr05_disc$NES[i], fdr05_disc$padj[i], fdr05_disc$direction[i]))
  }
  cat("\n")
  fdr05_out <- data.frame(pathway=fdr05_disc$pathway, pval=fdr05_disc$pval, padj=fdr05_disc$padj, ES=fdr05_disc$ES, NES=fdr05_disc$NES, size=fdr05_disc$size, direction=fdr05_disc$direction, stringsAsFactors=FALSE)
  utils::write.csv(fdr05_out, file.path(pw_dir, "STEP19E_GSE197677_HALLMARK_FDR05.csv"), row.names = FALSE)
}

if (nrow(fdr10_disc) > 0) {
  fdr10_out <- data.frame(pathway=fdr10_disc$pathway, pval=fdr10_disc$pval, padj=fdr10_disc$padj, ES=fdr10_disc$ES, NES=fdr10_disc$NES, size=fdr10_disc$size, direction=fdr10_disc$direction, stringsAsFactors=FALSE)
  utils::write.csv(fdr10_out, file.path(pw_dir, "STEP19E_GSE197677_HALLMARK_FDR10.csv"), row.names = FALSE)
}

disc_rank_freeze <- data.frame(
  pathway = fgsea_disc$pathway, discovery_rank = 1:nrow(fgsea_disc),
  NES = fgsea_disc$NES, pval = fgsea_disc$pval, padj = fgsea_disc$padj,
  direction = fgsea_disc$direction, size = fgsea_disc$size,
  leadingEdge = sapply(fgsea_disc$leadingEdge, function(x) paste(x, collapse = ",")),
  stringsAsFactors = FALSE
)
utils::write.csv(disc_rank_freeze, file.path(pw_dir, "STEP19E_GSE197677_HALLMARK_DISCOVERY_RANK_FREEZE.csv"), row.names = FALSE)
cat("Saved: STEP19E_GSE197677_HALLMARK_DISCOVERY_RANK_FREEZE.csv\n\n")

# =============================================================================
# 19E6: GSE221561 HALLMARK FGSEA
# =============================================================================

cat("=== 19E6: GSE221561 HALLMARK FGSEA ===\n\n")

repl_rank_fgsea <- repl_rank[repl_valid]
repl_rank_fgsea <- repl_rank_fgsea[!is.na(repl_rank_fgsea) & is.finite(repl_rank_fgsea)]

cat(sprintf("Running fgseaMultilevel: %d genes, %d pathways...\n", length(repl_rank_fgsea), length(hallmark_filtered)))

fgsea_repl <- fgseaMultilevel(
  pathways = hallmark_filtered, stats = repl_rank_fgsea,
  minSize = minSize, maxSize = maxSize, eps = 0
)

cat(sprintf("Completed. Pathways tested: %d\n\n", nrow(fgsea_repl)))

fgsea_repl$direction <- ifelse(fgsea_repl$NES > 0, "Higher_in_Neoadjuvant_treated", "Higher_in_Surgery_alone")
fgsea_repl <- fgsea_repl[order(fgsea_repl$padj), ]

fgsea_repl_out <- data.frame(
  pathway = fgsea_repl$pathway, pval = fgsea_repl$pval, padj = fgsea_repl$padj,
  ES = fgsea_repl$ES, NES = fgsea_repl$NES, size = fgsea_repl$size,
  leadingEdge = sapply(fgsea_repl$leadingEdge, function(x) paste(x, collapse = ",")),
  direction = fgsea_repl$direction, stringsAsFactors = FALSE
)
utils::write.csv(fgsea_repl_out, file.path(pw_dir, "STEP19E_GSE221561_HALLMARK_FGSEA_ALL.csv"), row.names = FALSE)
cat("Saved: STEP19E_GSE221561_HALLMARK_FGSEA_ALL.csv\n")

cat(sprintf("GSE221561 FDR < 0.05: %d pathways\n", sum(fgsea_repl$padj < 0.05)))
cat(sprintf("GSE221561 FDR < 0.10: %d pathways\n\n", sum(fgsea_repl$padj < 0.10)))

# =============================================================================
# 19E7: CROSS-DATASET HALLMARK COMPARISON
# =============================================================================

cat("=== 19E7: CROSS-DATASET HALLMARK COMPARISON ===\n\n")

# Merge by pathway name
disc_dt <- data.table(fgsea_disc_out[, c("pathway", "NES", "pval", "padj", "direction", "size", "leadingEdge")])
repl_dt <- data.table(fgsea_repl_out[, c("pathway", "NES", "pval", "padj", "direction", "leadingEdge")])

setnames(disc_dt, c("NES", "pval", "padj", "direction", "leadingEdge"),
         c("GSE197677_NES", "GSE197677_pval", "GSE197677_padj", "GSE197677_direction", "GSE197677_leadingEdge"))
setnames(repl_dt, c("NES", "pval", "padj", "direction", "leadingEdge"),
         c("GSE221561_NES", "GSE221561_pval", "GSE221561_padj", "GSE221561_direction", "GSE221561_leadingEdge"))

merged <- merge(disc_dt, repl_dt, by = "pathway", all = TRUE)

# Use numeric sign for NES direction (NOT string comparison)
merged$same_NES_direction <- sign(merged$GSE197677_NES) == sign(merged$GSE221561_NES)
merged$NES_delta <- merged$GSE197677_NES - merged$GSE221561_NES

# Add discovery rank
rank_lookup <- setNames(disc_rank_freeze$discovery_rank, disc_rank_freeze$pathway)
merged$discovery_rank <- rank_lookup[merged$pathway]

# Classify evidence
merged$pathway_evidence_class <- "EXPLORATORY_ONLY"
for (i in 1:nrow(merged)) {
  d_padj <- merged$GSE197677_padj[i]
  r_padj <- merged$GSE221561_padj[i]
  r_pval <- merged$GSE221561_pval[i]
  same_dir <- merged$same_NES_direction[i]

  if (!is.na(d_padj) && d_padj < 0.05 && !is.na(r_padj) && r_padj < 0.05 && same_dir) {
    merged$pathway_evidence_class[i] <- "CROSS_DATASET_FDR_CONSISTENT"
  } else if (!is.na(d_padj) && d_padj < 0.05 && !is.na(r_pval) && r_pval < 0.05 && same_dir) {
    merged$pathway_evidence_class[i] <- "DISCOVERY_FDR_REPLICATION_NOMINAL"
  } else if (!is.na(d_padj) && d_padj < 0.05 && same_dir) {
    merged$pathway_evidence_class[i] <- "DISCOVERY_FDR_SAME_DIRECTION_ONLY"
  } else if (!is.na(d_padj) && d_padj < 0.05 && !same_dir) {
    merged$pathway_evidence_class[i] <- "OPPOSITE_DIRECTION"
  }
}

merged <- merged[order(merged$discovery_rank), ]

merged_out <- data.frame(
  pathway = merged$pathway,
  discovery_rank = merged$discovery_rank,
  GSE197677_NES = merged$GSE197677_NES,
  GSE197677_pval = merged$GSE197677_pval,
  GSE197677_padj = merged$GSE197677_padj,
  GSE221561_NES = merged$GSE221561_NES,
  GSE221561_pval = merged$GSE221561_pval,
  GSE221561_padj = merged$GSE221561_padj,
  same_NES_direction = merged$same_NES_direction,
  NES_delta = merged$NES_delta,
  discovery_significant_FDR05 = merged$GSE197677_padj < 0.05,
  replication_nominal_P05 = merged$GSE221561_pval < 0.05,
  replication_significant_FDR05 = merged$GSE221561_padj < 0.05,
  pathway_evidence_class = merged$pathway_evidence_class,
  stringsAsFactors = FALSE
)

utils::write.csv(merged_out, file.path(pw_dir, "STEP19E_HALLMARK_CROSS_DATASET_COMPARISON.csv"), row.names = FALSE)
cat("Saved: STEP19E_HALLMARK_CROSS_DATASET_COMPARISON.csv\n\n")

# Evidence class summary
cat("Pathway evidence class summary:\n")
for (cls in unique(merged_out$pathway_evidence_class)) {
  cat(sprintf("  %s: %d\n", cls, sum(merged_out$pathway_evidence_class == cls)))
}
cat("\n")

# =============================================================================
# 19E9: GLOBAL PATHWAY CONCORDANCE
# =============================================================================

cat("=== 19E9: GLOBAL PATHWAY CONCORDANCE ===\n\n")

both_tested <- !is.na(merged$GSE197677_NES) & !is.na(merged$GSE221561_NES)
n_both <- sum(both_tested)
n_same <- sum(merged$same_NES_direction[both_tested], na.rm = TRUE)
frac_same <- n_same / n_both

p_binom <- binom.test(n_same, n_both, 0.5, alternative = "two.sided")$p.value
rho_nes <- cor(merged$GSE197677_NES[both_tested], merged$GSE221561_NES[both_tested], method = "spearman")

cat(sprintf("Hallmarks compared: %d\n", n_both))
cat(sprintf("Same NES direction: %d / %d = %.3f\n", n_same, n_both, frac_same))
cat(sprintf("Binomial P (two-sided): %.4f\n", p_binom))
cat(sprintf("NES Spearman rho: %.4f\n\n", rho_nes))

# =============================================================================
# 19E10: FROZEN TOP-RANK PATHWAYS
# =============================================================================

cat("=== 19E10: FROZEN TOP-RANK PATHWAYS ===\n\n")

for (k in c(5, 10, 20)) {
  top_k <- disc_rank_freeze$pathway[1:min(k, nrow(disc_rank_freeze))]
  sub <- merged[merged$pathway %in% top_k & both_tested, ]
  n_sub <- nrow(sub)
  n_same_sub <- sum(sub$same_NES_direction)
  frac_sub <- n_same_sub / n_sub
  p_sub <- if (n_sub >= 5) binom.test(n_same_sub, n_sub, 0.5, alternative = "greater")$p.value else NA

  cat(sprintf("Top %d (tested=%d): same_dir=%d, frac=%.3f", k, n_sub, n_same_sub, frac_sub))
  if (!is.na(p_sub)) cat(sprintf(", P=%.4f", p_sub))
  cat("\n")
}
cat("\n")

# =============================================================================
# 19E11: LEADING EDGE AUDIT
# =============================================================================

cat("=== 19E11: LEADING EDGE AUDIT ===\n\n")

# Pathways of interest: padj < 0.10 in discovery OR top 10
disc_of_interest <- fgsea_disc$pathway[fgsea_disc$padj < 0.10 | seq_len(nrow(fgsea_disc)) <= 10]
disc_of_interest <- unique(disc_of_interest)

le_overlap <- data.frame()

for (pw in disc_of_interest) {
  d_le_str <- merged$GSE197677_leadingEdge[merged$pathway == pw]
  r_le_str <- merged$GSE221561_leadingEdge[merged$pathway == pw]
  same_dir <- merged$same_NES_direction[merged$pathway == pw]

  if (is.na(d_le_str) || is.na(r_le_str)) next

  d_le <- unlist(strsplit(d_le_str, ","))
  r_le <- unlist(strsplit(r_le_str, ","))

  overlap_genes <- intersect(d_le, r_le)
  union_genes <- union(d_le, r_le)
  jaccard <- length(overlap_genes) / length(union_genes)

  le_overlap <- rbind(le_overlap, data.frame(
    pathway = pw,
    discovery_leading_edge_n = length(d_le),
    replication_leading_edge_n = length(r_le),
    overlap_n = length(overlap_genes),
    union_n = length(union_genes),
    jaccard = jaccard,
    same_NES_direction = same_dir,
    stringsAsFactors = FALSE
  ))
}

utils::write.csv(le_overlap, file.path(pw_dir, "STEP19E_HALLMARK_LEADING_EDGE_OVERLAP.csv"), row.names = FALSE)
cat(sprintf("Leading-edge audit: %d pathways analyzed\n", nrow(le_overlap)))
cat("Saved: STEP19E_HALLMARK_LEADING_EDGE_OVERLAP.csv\n\n")

# =============================================================================
# 19E12: GENE EVIDENCE INTEGRATION
# =============================================================================

cat("=== 19E12: GENE EVIDENCE INTEGRATION ===\n\n")

gene_evidence <- read.csv(file.path(repl_dir, "STEP19D_CROSS_DATASET_RANKED_GENE_EVIDENCE_FINAL.csv"),
                          stringsAsFactors = FALSE)

le_gene_evidence <- data.frame()

for (pw in disc_of_interest) {
  d_le_str <- merged$GSE197677_leadingEdge[merged$pathway == pw]
  if (is.na(d_le_str)) next

  d_le <- unlist(strsplit(d_le_str, ","))
  same_dir <- merged$same_NES_direction[merged$pathway == pw]

  for (gene in d_le) {
    ge_row <- gene_evidence[gene_evidence$gene == gene, ]
    if (nrow(ge_row) == 0) next

    le_gene_evidence <- rbind(le_gene_evidence, data.frame(
      pathway = pw,
      gene = gene,
      GSE197677_logFC = ge_row$GSE197677_logFC,
      GSE197677_FDR = ge_row$GSE197677_FDR,
      GSE197677_rank = ge_row$GSE197677_rank,
      GSE197677_LOSO = ge_row$GSE197677_LOSO_status,
      GSE221561_logFC = ge_row$GSE221561_logFC,
      GSE221561_PValue = ge_row$GSE221561_PValue,
      same_direction = ge_row$same_direction,
      GSE221561_LOSO = ge_row$GSE221561_LOSO_status,
      directional_evidence_class = ge_row$directional_evidence_class,
      same_NES_direction = same_dir,
      stringsAsFactors = FALSE
    ))
  }
}

utils::write.csv(le_gene_evidence, file.path(pw_dir, "STEP19E_HALLMARK_LEADING_EDGE_GENE_EVIDENCE.csv"), row.names = FALSE)
cat(sprintf("Leading-edge gene evidence: %d rows\n", nrow(le_gene_evidence)))
cat("Saved: STEP19E_HALLMARK_LEADING_EDGE_GENE_EVIDENCE.csv\n\n")

# =============================================================================
# 19E13: FIGURES
# =============================================================================

cat("=== 19E13: FIGURES ===\n\n")

# 1. NES concordance scatter
plot_data <- merged[both_tested, ]
p1 <- ggplot(plot_data, aes(x = GSE197677_NES, y = GSE221561_NES)) +
  geom_point(aes(color = same_NES_direction), alpha = 0.7, size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_smooth(method = "lm", se = FALSE, color = "red", linewidth = 1) +
  scale_color_manual(values = c("TRUE" = "steelblue", "FALSE" = "tomato"),
                     labels = c("TRUE" = "Same direction", "FALSE" = "Opposite direction")) +
  labs(title = "Hallmark Pathway NES Concordance",
       subtitle = "GSE197677 discovery vs GSE221561 replication",
       x = "GSE197677 NES", y = "GSE221561 NES", color = "") +
  theme_minimal()
ggsave(file.path(fig_dir, "STEP19E_HALLMARK_NES_CONCORDANCE.png"), p1, width = 8, height = 7, dpi = 300)
cat("Saved: STEP19E_HALLMARK_NES_CONCORDANCE.png\n")

# 2. GSE197677 dotplot (top 15)
top15 <- head(fgsea_disc, 15)
p2 <- ggplot(top15, aes(x = NES, y = reorder(pathway, NES))) +
  geom_point(aes(size = size, color = padj)) +
  scale_color_gradient(low = "red", high = "blue") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  labs(title = "GSE197677 Hallmark Pathway Enrichment",
       subtitle = "NACT vs nNACT (discovery)",
       x = "NES", y = "", color = "padj", size = "Gene set size") +
  theme_minimal()
ggsave(file.path(fig_dir, "STEP19E_GSE197677_HALLMARK_DOTPLOT.png"), p2, width = 10, height = 8, dpi = 300)
cat("Saved: STEP19E_GSE197677_HALLMARK_DOTPLOT.png\n")

# 3. GSE221561 dotplot (top 15)
top15_repl <- head(fgsea_repl, 15)
p3 <- ggplot(top15_repl, aes(x = NES, y = reorder(pathway, NES))) +
  geom_point(aes(size = size, color = padj)) +
  scale_color_gradient(low = "red", high = "blue") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  labs(title = "GSE221561 Hallmark Pathway Enrichment",
       subtitle = "Neoadjuvant_treated vs Surgery_alone (replication; Surgery_alone n=2)",
       x = "NES", y = "", color = "padj", size = "Gene set size") +
  theme_minimal()
ggsave(file.path(fig_dir, "STEP19E_GSE221561_HALLMARK_DOTPLOT.png"), p3, width = 10, height = 8, dpi = 300)
cat("Saved: STEP19E_GSE221561_HALLMARK_DOTPLOT.png\n")

# 4. Heatmap (top 20)
top20_pw <- disc_rank_freeze$pathway[1:min(20, nrow(disc_rank_freeze))]
heatmap_data <- merged[merged$pathway %in% top20_pw, c("pathway", "GSE197677_NES", "GSE221561_NES")]
heatmap_data <- heatmap_data[order(heatmap_data$GSE197677_NES), ]

heatmap_matrix <- as.matrix(heatmap_data[, c("GSE197677_NES", "GSE221561_NES")])
rownames(heatmap_matrix) <- heatmap_data$pathway

png(file.path(fig_dir, "STEP19E_CROSS_DATASET_HALLMARK_HEATMAP.png"), width = 600, height = 800, res = 150)
par(mar = c(2, 12, 2, 2))
image(t(heatmap_matrix), col = colorRampPalette(c("tomato", "white", "steelblue"))(100),
      axes = FALSE, main = "Hallmark NES Heatmap (Top 20)")
axis(2, at = seq(0, 1, length.out = nrow(heatmap_matrix)), labels = rownames(heatmap_matrix), las = 2, cex.axis = 0.7)
axis(1, at = c(0, 1), labels = colnames(heatmap_matrix), cex.axis = 0.8)
dev.off()
cat("Saved: STEP19E_CROSS_DATASET_HALLMARK_HEATMAP.png\n\n")

# =============================================================================
# 19E14: FINAL PATHWAY DECISION
# =============================================================================

cat("=== 19E14: FINAL PATHWAY DECISION ===\n\n")

n_tested <- length(hallmark_filtered)
n_fdr05_disc <- sum(fgsea_disc$padj < 0.05, na.rm = TRUE)
n_fdr10_disc <- sum(fgsea_disc$padj < 0.10, na.rm = TRUE)
n_fdr05_repl <- sum(fgsea_repl$padj < 0.05, na.rm = TRUE)

n_cross_fdr <- sum(merged_out$pathway_evidence_class == "CROSS_DATASET_FDR_CONSISTENT", na.rm = TRUE)
n_disc_fdr_nominal <- sum(merged_out$pathway_evidence_class == "DISCOVERY_FDR_REPLICATION_NOMINAL", na.rm = TRUE)
n_disc_fdr_dir <- sum(merged_out$pathway_evidence_class == "DISCOVERY_FDR_SAME_DIRECTION_ONLY", na.rm = TRUE)

# Determine status
if (n_cross_fdr > 0) {
  overall_status <- "FORMAL_PATHWAY_DISCOVERY_WITH_CROSS_DATASET_SUPPORT"
} else if (n_fdr05_disc > 0 && (n_disc_fdr_nominal + n_disc_fdr_dir > 0)) {
  overall_status <- "FORMAL_PATHWAY_DISCOVERY_WITHOUT_REPLICATION_SUPPORT"
} else if (n_fdr05_disc == 0 && frac_same > 0.55) {
  overall_status <- "EXPLORATORY_PATHWAY_CONCORDANCE_ONLY"
} else {
  overall_status <- "NO_CLEAR_PATHWAY_SIGNAL"
}

cat(sprintf("Overall status: %s\n\n", overall_status))

# Decision table
decision <- data.frame(
  metric = c("n_hallmarks_tested", "GSE197677_FDR05", "GSE197677_FDR10", "GSE221561_FDR05",
             "discovery_FDR05_same_direction_replication", "discovery_FDR05_replication_P05",
             "both_FDR05_same_direction", "global_same_direction_fraction",
             "global_binomial_P", "NES_spearman_rho", "overall_status"),
  value = as.character(c(n_tested, n_fdr05_disc, n_fdr10_disc, n_fdr05_repl,
                         n_disc_fdr_dir + n_cross_fdr, n_disc_fdr_nominal + n_cross_fdr,
                         n_cross_fdr, round(frac_same, 4), round(p_binom, 4),
                         round(rho_nes, 4), overall_status)),
  stringsAsFactors = FALSE
)

utils::write.csv(decision, file.path(pw_dir, "STEP19E_HALLMARK_PATHWAY_DECISION.csv"), row.names = FALSE)
cat("Saved: STEP19E_HALLMARK_PATHWAY_DECISION.csv\n\n")

# =============================================================================
# 19E15: INTERPRETATION FREEZE
# =============================================================================

cat("=== 19E15: INTERPRETATION FREEZE ===\n\n")

interpretation <- sprintf("# Step 19E Hallmark Pathway Interpretation

## Date
%s

## Method
- Preranked gene-level analysis using fgseaMultilevel
- Ranking statistic: sign(logFC) * sqrt(F)
- Gene universe: %d common genes tested in both datasets
- Pathway database: MSigDB Hallmark (50 gene sets)
- minSize = 15, maxSize = 500, eps = 0

## Key Findings

### GSE197677 Discovery
- FDR < 0.05 pathways: %d
- FDR < 0.10 pathways: %d

### GSE221561 Replication
- FDR < 0.05 pathways: %d
- Note: Surgery_alone n=2 (small reference group)

### Cross-Dataset Concordance
- Hallmarks compared: %d
- Same NES direction: %d/%d (%.1f%%)
- Binomial P: %.4f
- NES Spearman rho: %.4f

### Pathway Evidence Classes
- CROSS_DATASET_FDR_CONSISTENT: %d
- DISCOVERY_FDR_REPLICATION_NOMINAL: %d
- DISCOVERY_FDR_SAME_DIRECTION_ONLY: %d
- OPPOSITE_DIRECTION: %d
- EXPLORATORY_ONLY: %d

### Overall Status
**%s**

## Important Caveats

1. Analysis uses preranked genes, not a nominal DEG cutoff
2. Same %d-gene universe used in both datasets
3. Datasets were analyzed independently
4. Samples were not pooled
5. Positive NES corresponds to treatment-associated direction within each dataset
6. Clinical contrasts are related but non-identical (NACT vs nNACT; Neoadjuvant_treated vs Surgery_alone)
7. GSE221561 Surgery_alone n=2
8. GSE197677 had zero single-gene FDR<0.05 DEGs
9. Pathway enrichment does not make individual genes significant
10. No causal treatment claim
11. Epithelial compartment is not equivalent to CNV-confirmed malignant epithelium
",
Sys.Date(),
n_both,
n_fdr05_disc, n_fdr10_disc, n_fdr05_repl,
n_both, n_same, n_both, frac_same * 100,
p_binom, rho_nes,
n_cross_fdr, n_disc_fdr_nominal, n_disc_fdr_dir,
sum(merged_out$pathway_evidence_class == "OPPOSITE_DIRECTION", na.rm = TRUE),
sum(merged_out$pathway_evidence_class == "EXPLORATORY_ONLY", na.rm = TRUE),
overall_status,
n_both
)

writeLines(interpretation, file.path(pw_dir, "STEP19E_HALLMARK_INTERPRETATION.md"))
cat("Saved: STEP19E_HALLMARK_INTERPRETATION.md\n\n")

# =============================================================================
# 19E16: FINAL FREEZE
# =============================================================================

cat("=== 19E16: FINAL FREEZE ===\n\n")

freeze_final <- data.frame(
  pathway = merged_out$pathway,
  discovery_rank = merged_out$discovery_rank,
  GSE197677_NES = merged_out$GSE197677_NES,
  GSE197677_pval = merged_out$GSE197677_pval,
  GSE197677_padj = merged_out$GSE197677_padj,
  GSE221561_NES = merged_out$GSE221561_NES,
  GSE221561_pval = merged_out$GSE221561_pval,
  GSE221561_padj = merged_out$GSE221561_padj,
  same_direction = merged_out$same_NES_direction,
  leading_edge_jaccard = NA_real_,
  pathway_evidence_class = merged_out$pathway_evidence_class,
  primary_collection = "HALLMARK",
  stringsAsFactors = FALSE
)

# Fill Jaccard from le_overlap
for (i in 1:nrow(freeze_final)) {
  pw <- freeze_final$pathway[i]
  le_row <- le_overlap[le_overlap$pathway == pw, ]
  if (nrow(le_row) > 0) freeze_final$leading_edge_jaccard[i] <- le_row$jaccard
}

utils::write.csv(freeze_final, file.path(pw_dir, "STEP19E_HALLMARK_CROSS_DATASET_FREEZE.csv"), row.names = FALSE)
cat("Saved: STEP19E_HALLMARK_CROSS_DATASET_FREEZE.csv\n\n")

# =============================================================================
# FINAL OUTPUT
# =============================================================================

cat("\n")
cat("============================================\n")
cat("STEP 19E COMPLETE\n")
cat("HALLMARK PRERANKED PATHWAY ANALYSIS\n")
cat("============================================\n\n")

cat(sprintf("Packages:\n"))
cat(sprintf("  fgsea:      %s\n", as.character(packageVersion("fgsea"))))
cat(sprintf("  msigdbr:    %s\n\n", as.character(packageVersion("msigdbr"))))

cat(sprintf("Common gene universe         : %d\n", n_both))
cat(sprintf("Hallmark pathways tested     : %d\n\n", n_tested))

cat("--------------------------------------------\n\n")
cat("GSE197677 DISCOVERY\n\n")
cat(sprintf("FDR < 0.05 pathways          : %d\n", n_fdr05_disc))
cat(sprintf("FDR < 0.10 pathways          : %d\n\n", n_fdr10_disc))

cat("Top 5 discovery pathways:\n\n")
cat("Pathway | NES | P | FDR | Direction\n")
cat("--------|-----|---|-----|----------\n")
for (i in 1:min(5, nrow(fgsea_disc))) {
  cat(sprintf("%s | %.3f | %.4f | %.4f | %s\n",
              fgsea_disc$pathway[i], fgsea_disc$NES[i], fgsea_disc$pval[i],
              fgsea_disc$padj[i], fgsea_disc$direction[i]))
}

cat("\n--------------------------------------------\n\n")
cat("GSE221561 REPLICATION\n\n")
cat(sprintf("FDR < 0.05 pathways          : %d\n\n", n_fdr05_repl))

cat("--------------------------------------------\n\n")
cat("CROSS-DATASET\n\n")
cat(sprintf("Discovery FDR05 pathways     : %d\n", n_fdr05_disc))
cat(sprintf("Same NES direction           : %d\n", n_disc_fdr_dir + n_cross_fdr))
cat(sprintf("Replication P < 0.05         : %d\n", n_disc_fdr_nominal + n_cross_fdr))
cat(sprintf("Both FDR < 0.05              : %d\n\n", n_cross_fdr))

cat("All Hallmarks:\n")
cat(sprintf("  Same direction             : %d / %d = %.3f\n", n_same, n_both, frac_same))
cat(sprintf("  Binomial P                 : %.4f\n", p_binom))
cat(sprintf("  NES Spearman rho           : %.4f\n\n", rho_nes))

# Top10
top10_pw <- disc_rank_freeze$pathway[1:min(10, nrow(disc_rank_freeze))]
sub10 <- merged[merged$pathway %in% top10_pw & both_tested, ]
cat("Top10 frozen discovery:\n")
cat(sprintf("  Same direction             : %d / %d = %.3f\n\n",
            sum(sub10$same_NES_direction), nrow(sub10), sum(sub10$same_NES_direction)/nrow(sub10)))

cat("--------------------------------------------\n\n")

# Top cross-dataset supported pathways
cat("Top cross-dataset supported pathways:\n\n")
supported <- merged_out[merged_out$pathway_evidence_class %in% c("CROSS_DATASET_FDR_CONSISTENT", "DISCOVERY_FDR_REPLICATION_NOMINAL"), ]
if (nrow(supported) > 0) {
  for (i in 1:min(10, nrow(supported))) {
    cat(sprintf("  %s\n", supported$pathway[i]))
    cat(sprintf("    Discovery NES/FDR: %.3f / %.4f\n", supported$GSE197677_NES[i], supported$GSE197677_padj[i]))
    cat(sprintf("    Replication NES/P/FDR: %.3f / %.4f / %.4f\n", supported$GSE221561_NES[i], supported$GSE221561_pval[i], supported$GSE221561_padj[i]))
    cat(sprintf("    Same direction: %s\n", supported$same_NES_direction[i]))
    cat(sprintf("    Evidence class: %s\n\n", supported$pathway_evidence_class[i]))
  }
} else {
  cat("  None\n\n")
}

cat("--------------------------------------------\n\n")
cat(sprintf("Overall pathway status:\n%s\n\n", overall_status))
cat("Single-gene replicated DEGs claimed: NO\n")
cat("Counts pooled: NO\n")
cat("Datasets integrated: NO\n")
cat("GSE160269 included: NO\n")
cat("GO/KEGG/Reactome run: NO\n\n")
cat("STOP HERE.\n")
cat("DO NOT START STEP 19F AUTOMATICALLY.\n")
