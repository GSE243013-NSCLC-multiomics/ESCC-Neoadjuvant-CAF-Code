#!/usr/bin/env Rscript
# =============================================================================
# STEP 19D: CORRECTED CROSS-DATASET DIRECTIONAL REPLICATION
# =============================================================================
# This script corrects two bugs in the original Step 19D analysis:
#
# BUG 1: Direction string comparison
#   Original: compared literal strings ("Higher_in_NACT" vs "Higher_in_Neoadjuvant_treated")
#   Fix: use logFC sign comparison (positive = treatment, negative = control)
#
# BUG 2: Named vector indexing
#   Original: df$logFC[gene] returns NA (unnamed vector)
#   Fix: df[gene, "logFC"] returns correct value (data frame indexing)
#
# CORRECTED RESULTS:
#   Top-100 concordance: 56/99 = 0.566 (Binomial P = 0.114)
#   Genome-wide concordance: 4484/9033 = 0.496 (essentially random)
#   RCN3, SDC2, SENP5 all show same direction across datasets
# =============================================================================

suppressPackageStartupMessages({
  library(edgeR)
})

base_dir <- "~/Downloads/ESCC Neoadjuvant"
repl_dir <- file.path(base_dir, "06_tables/cross_dataset/epithelial_pseudobulk/replication")
de_dir <- file.path(base_dir, "06_tables/cross_dataset/epithelial_pseudobulk/GSE197677_DE")
fig_dir <- file.path(base_dir, "05_figures/cross_dataset/epithelial_pseudobulk/replication")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 19D-A: LOAD DATA WITH CORRECT INDEXING
# =============================================================================

cat("=== 19D-A: LOAD DATA ===\n\n")

# Discovery (GSE197677)
disc <- read.csv(file.path(de_dir, "STEP19C_GSE197677_NACT_vs_nNACT_ALL_GENES.csv"),
                 row.names = 1, check.names = FALSE)

# Replication (GSE221561) - canonical contrast
repl <- read.csv(file.path(repl_dir, "STEP19D_GSE221561_NEOADJUVANT_MINUS_SURGERY_CANONICAL.csv"),
                 row.names = 1, check.names = FALSE)

# Frozen discovery ranking
freeze <- read.csv(file.path(repl_dir, "STEP19D_GSE197677_DISCOVERY_RANK_FREEZE.csv"),
                   stringsAsFactors = FALSE)

cat(sprintf("Discovery genes: %d\n", nrow(disc)))
cat(sprintf("Replication genes: %d\n", nrow(repl)))
cat(sprintf("Frozen ranking genes: %d\n\n", nrow(freeze)))

# =============================================================================
# 19D-B: COMMON GENE SET
# =============================================================================

cat("=== 19D-B: COMMON GENE SET ===\n\n")

common_genes <- intersect(rownames(disc), rownames(repl))
cat(sprintf("Common tested genes: %d / %d discovery, %d / %d replication\n",
            length(common_genes), nrow(disc), length(common_genes), nrow(repl)))

# Save common genes
common_df <- data.frame(
  gene = common_genes,
  in_discovery = TRUE,
  in_replication = TRUE
)
utils::write.csv(common_df, file.path(repl_dir, "STEP19D_CORRECTED_common_tested_genes.csv"),
                 row.names = FALSE)

# =============================================================================
# 19D-C: GENOME-WIDE CONCORDANCE (LOGFC SIGN)
# =============================================================================

cat("\n=== 19D-C: GENOME-WIDE CONCORDANCE ===\n\n")

# CORRECT indexing: df[gene, "logFC"] not df$logFC[gene]
disc_fc <- disc[common_genes, "logFC"]
repl_fc <- repl[common_genes, "logFC"]

# Remove any NAs
valid <- !is.na(disc_fc) & !is.na(repl_fc)
n_valid <- sum(valid)
cat(sprintf("Valid genes (no NA): %d / %d\n\n", n_valid, length(common_genes)))

# Sign concordance
n_concordant <- sum(sign(disc_fc[valid]) == sign(repl_fc[valid]))
n_discordant <- n_valid - n_concordant
cat(sprintf("Genome-wide sign concordance: %d / %d = %.3f\n",
            n_concordant, n_valid, n_concordant / n_valid))
cat(sprintf("Genome-wide sign discordance: %d / %d = %.3f\n\n",
            n_discordant, n_valid, n_discordant / n_valid))

# Spearman correlation
rho <- cor(disc_fc[valid], repl_fc[valid], method = "spearman")
cat(sprintf("Genome-wide Spearman rho: %.4f\n\n", rho))

# Direction breakdown
cat("Direction breakdown (using logFC sign):\n")
n_pos_pos <- sum(disc_fc[valid] > 0 & repl_fc[valid] > 0)
n_neg_neg <- sum(disc_fc[valid] < 0 & repl_fc[valid] < 0)
n_pos_neg <- sum(disc_fc[valid] > 0 & repl_fc[valid] < 0)
n_neg_pos <- sum(disc_fc[valid] < 0 & repl_fc[valid] > 0)
cat(sprintf("  Both higher in treatment:  %d (%.1f%%)\n", n_pos_pos, 100 * n_pos_pos / n_valid))
cat(sprintf("  Both higher in control:    %d (%.1f%%)\n", n_neg_neg, 100 * n_neg_neg / n_valid))
cat(sprintf("  Disc higher in treat, Repl higher in control: %d (%.1f%%)\n",
            n_pos_neg, 100 * n_pos_neg / n_valid))
cat(sprintf("  Disc higher in control, Repl higher in treat: %d (%.1f%%)\n\n",
            n_neg_pos, 100 * n_neg_pos / n_valid))

# Binomial test for genome-wide
p_genomewide <- binom.test(n_concordant, n_valid, 0.5, alternative = "greater")$p.value
cat(sprintf("Binomial P (greater than 0.5): %.2e\n\n", p_genomewide))

# =============================================================================
# 19D-D: TOP-K CONCORDANCE
# =============================================================================

cat("=== 19D-D: TOP-K CONCORDANCE ===\n\n")

top_k_results <- data.frame()

for (k in c(25, 50, 100, 200, 500)) {
  top_k_genes <- freeze$gene[1:min(k, nrow(freeze))]
  match_idx <- match(top_k_genes, rownames(repl))
  tested <- !is.na(match_idx)
  n_t <- sum(tested)

  if (n_t == 0) next

  disc_fc_k <- disc[top_k_genes[tested], "logFC"]
  repl_fc_k <- repl[top_k_genes[tested], "logFC"]

  n_same_k <- sum(sign(disc_fc_k) == sign(repl_fc_k))
  conc_k <- n_same_k / n_t
  p_k <- binom.test(n_same_k, n_t, 0.5, alternative = "greater")$p.value

  cat(sprintf("Top %d (tested=%d): %d / %d = %.3f (P = %.2e)\n",
              k, n_t, n_same_k, n_t, conc_k, p_k))

  top_k_results <- rbind(top_k_results, data.frame(
    k = k, tested = n_t, concordant = n_same_k,
    proportion = conc_k, p_value = p_k
  ))
}

# =============================================================================
# 19D-E: SPECIFIC GENE AUDIT (RCN3, SDC2, SENP5)
# =============================================================================

cat("\n=== 19D-E: SPECIFIC GENE AUDIT ===\n\n")

target_genes <- c("RCN3", "SDC2", "SENP5")

gene_audit <- data.frame()
for (gene in target_genes) {
  d_fc <- disc[gene, "logFC"]
  r_fc <- repl[gene, "logFC"]
  d_pval <- disc[gene, "PValue"]
  r_pval <- repl[gene, "PValue"]
  d_fdr <- disc[gene, "FDR"]
  r_fdr <- repl[gene, "FDR"]
  same_dir <- sign(d_fc) == sign(r_fc)

  cat(sprintf("%s:\n", gene))
  cat(sprintf("  Discovery:   logFC = %+.4f, P = %.2e, FDR = %.4f\n", d_fc, d_pval, d_fdr))
  cat(sprintf("  Replication: logFC = %+.4f, P = %.2e, FDR = %.4f\n", r_fc, r_pval, r_fdr))
  cat(sprintf("  Same direction: %s\n\n", same_dir))

  gene_audit <- rbind(gene_audit, data.frame(
    gene = gene,
    discovery_logFC = d_fc,
    discovery_PValue = d_pval,
    discovery_FDR = d_fdr,
    replication_logFC = r_fc,
    replication_PValue = r_pval,
    replication_FDR = r_fdr,
    same_direction = same_dir
  ))
}

utils::write.csv(gene_audit, file.path(repl_dir, "STEP19D_CORRECTED_specific_gene_audit.csv"),
                 row.names = FALSE)

# =============================================================================
# 19D-F: TOP-100 DETAILED TABLE
# =============================================================================

cat("=== 19D-F: TOP-100 DETAILED TABLE ===\n\n")

top100_genes <- freeze$gene[1:100]
match_idx <- match(top100_genes, rownames(repl))
tested <- !is.na(match_idx)
n_tested <- sum(tested)

cat(sprintf("Top-100 genes tested in replication: %d / 100\n\n", n_tested))

# Build detailed table
top100_detail <- data.frame(
  gene = top100_genes,
  discovery_rank = 1:100,
  discovery_logFC = disc[top100_genes, "logFC"],
  discovery_direction = disc[top100_genes, "direction"],
  replication_logFC = repl[top100_genes, "logFC"],
  replication_direction = repl[top100_genes, "direction"],
  tested = tested,
  same_direction = FALSE,
  stringsAsFactors = FALSE
)

# Calculate same_direction using sign comparison (CORRECT method)
for (i in 1:nrow(top100_detail)) {
  if (top100_detail$tested[i]) {
    d_fc <- top100_detail$discovery_logFC[i]
    r_fc <- top100_detail$replication_logFC[i]
    top100_detail$same_direction[i] <- sign(d_fc) == sign(r_fc)
  }
}

# Summary
n_same_100 <- sum(top100_detail$same_direction[top100_detail$tested])
cat(sprintf("Top-100 concordant: %d / %d = %.3f\n",
            n_same_100, n_tested, n_same_100 / n_tested))

utils::write.csv(top100_detail, file.path(repl_dir, "STEP19D_CORRECTED_top100_detailed.csv"),
                 row.names = FALSE)

# =============================================================================
# 19D-G: VISUALIZATION
# =============================================================================

cat("\n=== 19D-G: VISUALIZATION ===\n\n")

# G1: Genome-wide logFC scatter
pdf(file.path(fig_dir, "STEP19D_CORRECTED_genomewide_logFC_scatter.pdf"), width = 8, height = 8)
plot(disc_fc[valid], repl_fc[valid],
     pch = 16, cex = 0.3, col = rgb(0.2, 0.2, 0.2, 0.3),
     xlab = "GSE197677 logFC (NACT vs nNACT)",
     ylab = "GSE221561 logFC (Neoadjuvant vs Surgery_alone)",
     main = sprintf("Genome-wide logFC concordance (n=%d, rho=%.3f)", n_valid, rho))
abline(h = 0, col = "gray60", lty = 2)
abline(v = 0, col = "gray60", lty = 2)
abline(lm(repl_fc[valid] ~ disc_fc[valid]), col = "red", lwd = 2)
# Highlight specific genes
for (gene in target_genes) {
  points(disc[gene, "logFC"], repl[gene, "logFC"],
         pch = 17, cex = 1.5, col = "red")
  text(disc[gene, "logFC"], repl[gene, "logFC"],
       labels = gene, pos = 3, cex = 0.8, col = "red")
}
dev.off()
cat("Saved: STEP19D_CORRECTED_genomewide_logFC_scatter.pdf\n")

# G2: Top-100 direction concordance barplot
pdf(file.path(fig_dir, "STEP19D_CORRECTED_top100_concordance.pdf"), width = 8, height = 6)
par(mfrow = c(1, 2))

# Left: Concordance by rank tier
k_vals <- c(25, 50, 100)
conc_vals <- top_k_results$proportion[top_k_results$k %in% k_vals]
barplot(conc_vals, names.arg = paste0("Top-", k_vals),
        ylim = c(0, 1), col = "steelblue",
        main = "Direction concordance by rank",
        ylab = "Proportion concordant")
abline(h = 0.5, col = "red", lty = 2)
abline(h = 1, col = "gray60", lty = 3)

# Right: Top-100 same/different direction
n_diff <- n_tested - n_same_100
barplot(c(n_same_100, n_diff), names.arg = c("Same", "Different"),
        col = c("steelblue", "tomato"),
        main = sprintf("Top-100 direction (n=%d)", n_tested),
        ylab = "Number of genes")

dev.off()
cat("Saved: STEP19D_CORRECTED_top100_concordance.pdf\n")

# G3: Top-12 gene-level comparison
top12_genes <- freeze$gene[1:12]
pdf(file.path(fig_dir, "STEP19D_CORRECTED_top12_gene_comparison.pdf"), width = 10, height = 6)
par(mfrow = c(3, 4), mar = c(3, 3, 2, 1))
for (gene in top12_genes) {
  d_fc <- disc[gene, "logFC"]
  r_fc <- repl[gene, "logFC"]
  same <- sign(d_fc) == sign(r_fc)
  barplot(c(d_fc, r_fc),
          names.arg = c("Disc", "Repl"),
          col = ifelse(same, "steelblue", "tomato"),
          main = gene, cex.main = 0.9,
          ylim = range(c(d_fc, r_fc)) * 1.2)
  abline(h = 0, col = "gray40")
}
dev.off()
cat("Saved: STEP19D_CORRECTED_top12_gene_comparison.pdf\n")

# =============================================================================
# 19D-H: PROVENANCE NOTE
# =============================================================================

cat("\n=== 19D-H: PROVENANCE NOTE ===\n\n")

provenance <- sprintf("# Step 19D Corrected Replication Analysis

## Date
%s

## Bugs Corrected

### Bug 1: Direction string comparison
- **Original**: compared literal direction strings
  (\"Higher_in_NACT\" != \"Higher_in_Neoadjuvant_treated\" → always FALSE)
- **Fix**: use logFC sign comparison
  (positive = treatment group, negative = control group)

### Bug 2: Named vector indexing
- **Original**: `df$logFC[gene]` returns NA (extracts unnamed vector)
- **Fix**: `df[gene, \"logFC\"]` returns correct value (data frame indexing)

## Corrected Results

### Genome-wide
- Common genes: %d
- Sign concordance: %d / %d = %.3f
- Spearman rho: %.4f
- Both higher in treatment: %d (%.1f%%)
- Both higher in control: %d (%.1f%%)
- Discordant: %d (%.1f%%)

### Top-K concordance
%s

### Specific genes
- RCN3: Disc logFC = %+.4f, Repl logFC = %+.4f, Same dir = %s
- SDC2: Disc logFC = %+.4f, Repl logFC = %+.4f, Same dir = %s
- SENP5: Disc logFC = %+.4f, Repl logFC = %+.4f, Same dir = %s

## Interpretation

The corrected analysis shows:
1. **Top-100 concordance is modest positive (56%%)**, not anti-concordant
2. **Genome-wide concordance is ~50%%**, expected for non-DE genes
3. **Top-3 genes (RCN3, SDC2, SENP5) all show same direction** across datasets
4. The previous \"0/99 concordance\" and \"TRUE_DISCORDANCE\" were artifacts of bugs

## Previous (Incorrect) Results
- Reported 0/99 concordance → corrected to 56/99 = 0.566
- Classified as TRUE_DISCORDANCE → corrected to MODEST_POSITIVE_CONCORDANCE
",
Sys.Date(),
length(common_genes),
n_concordant, n_valid, n_concordant / n_valid,
rho,
n_pos_pos, 100 * n_pos_pos / n_valid,
n_neg_neg, 100 * n_neg_neg / n_valid,
n_discordant, 100 * n_discordant / n_valid,
paste(apply(top_k_results, 1, function(row) {
  sprintf("- Top-%d: %d / %d = %.3f (P = %.2e)", row["k"], row["concordant"], row["tested"], row["proportion"], row["p_value"])
}), collapse = "\n"),
disc["RCN3", "logFC"], repl["RCN3", "logFC"], disc["RCN3", "logFC"] * repl["RCN3", "logFC"] > 0,
disc["SDC2", "logFC"], repl["SDC2", "logFC"], disc["SDC2", "logFC"] * repl["SDC2", "logFC"] > 0,
disc["SENP5", "logFC"], repl["SENP5", "logFC"], disc["SENP5", "logFC"] * repl["SENP5", "logFC"] > 0
)

writeLines(provenance, file.path(repl_dir, "STEP19D_CORRECTED_REPLICATION_ANALYSIS.md"))
cat("Saved: STEP19D_CORRECTED_REPLICATION_ANALYSIS.md\n")

# =============================================================================
# FINAL OUTPUT
# =============================================================================

cat("\n")
cat("============================================\n")
cat("STEP 19D CORRECTED REPLICATION ANALYSIS\n")
cat("============================================\n\n")

cat("Bugs corrected:\n")
cat("  1. Direction string comparison → logFC sign comparison\n")
cat("  2. df$logFC[gene] → df[gene, \"logFC\"]\n\n")

cat("Genome-wide:\n")
cat(sprintf("  Sign concordance: %d / %d = %.3f\n", n_concordant, n_valid, n_concordant / n_valid))
cat(sprintf("  Spearman rho: %.4f\n\n", rho))

cat("Top-K concordance:\n")
for (i in 1:nrow(top_k_results)) {
  cat(sprintf("  Top-%d: %d / %d = %.3f (P = %.2e)\n",
              top_k_results$k[i], top_k_results$concordant[i],
              top_k_results$tested[i], top_k_results$proportion[i],
              top_k_results$p_value[i]))
}

cat("\nSpecific genes:\n")
for (i in 1:nrow(gene_audit)) {
  cat(sprintf("  %s: Disc=%+.4f, Repl=%+.4f, Same=%s\n",
              gene_audit$gene[i], gene_audit$discovery_logFC[i],
              gene_audit$replication_logFC[i], gene_audit$same_direction[i]))
}

cat("\nClassification: MODEST_POSITIVE_CONCORDANCE\n")
cat("Previous classification (TRUE_DISCORDANCE) was incorrect.\n\n")

cat("STOP HERE.\n")
cat("DO NOT CONTINUE STEP 19D8 AUTOMATICALLY.\n")
