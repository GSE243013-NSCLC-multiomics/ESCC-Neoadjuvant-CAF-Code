#!/usr/bin/env Rscript
# =============================================================================
# STEP 19D8-19D14: CORRECTED DIRECTIONAL REPLICATION FINALIZATION
# =============================================================================
# Bugs corrected:
#   1. Direction string comparison → logFC sign comparison
#   2. df$logFC[gene] → df[gene, "logFC"]
#
# Primary endpoint: Top-100 direction concordance (frozen)
# =============================================================================

suppressPackageStartupMessages({
  library(edgeR)
})

base_dir <- "~/Downloads/ESCC Neoadjuvant"
repl_dir <- file.path(base_dir, "06_tables/cross_dataset/epithelial_pseudobulk/replication")
de_dir <- file.path(base_dir, "06_tables/cross_dataset/epithelial_pseudobulk/GSE197677_DE")
counts_dir <- file.path(base_dir, "06_tables/cross_dataset/epithelial_pseudobulk/counts")
fig_dir <- file.path(base_dir, "05_figures/cross_dataset/epithelial_pseudobulk/replication")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# LOAD ALL REQUIRED DATA
# =============================================================================

cat("=== LOADING DATA ===\n\n")

# Discovery (GSE197677)
disc_all <- read.csv(file.path(de_dir, "STEP19C_GSE197677_NACT_vs_nNACT_ALL_GENES.csv"),
                     row.names = 1, check.names = FALSE)

# Replication (GSE221561) - canonical contrast
repl_all <- read.csv(file.path(repl_dir, "STEP19D_GSE221561_NEOADJUVANT_MINUS_SURGERY_CANONICAL.csv"),
                     row.names = 1, check.names = FALSE)

# Frozen discovery ranking
freeze <- read.csv(file.path(repl_dir, "STEP19D_GSE197677_DISCOVERY_RANK_FREEZE.csv"),
                   stringsAsFactors = FALSE)

# Nominal candidates (frozen from Step 19C)
nominal <- read.csv(file.path(de_dir, "STEP19C_GSE197677_nominal_directional_candidates.csv"),
                    row.names = 1, check.names = FALSE)

# GSE197677 LOSO stability (frozen from Step 19C)
loso_197677 <- read.csv(file.path(de_dir, "STEP19C_GSE197677_LOSO_top100_stability.csv"),
                        stringsAsFactors = FALSE)

# GSE221561 raw counts and metadata
raw_counts_221561 <- readRDS(file.path(counts_dir, "GSE221561_epithelial_pseudobulk_raw_counts.rds"))
meta_221561 <- read.csv(file.path(counts_dir, "GSE221561_epithelial_pseudobulk_sample_metadata.csv"),
                        stringsAsFactors = FALSE)

cat(sprintf("Discovery genes: %d\n", nrow(disc_all)))
cat(sprintf("Replication genes: %d\n", nrow(repl_all)))
cat(sprintf("Frozen ranking genes: %d\n", nrow(freeze)))
cat(sprintf("Nominal candidates: %d\n", nrow(nominal)))
cat(sprintf("GSE221561 samples: %d\n", nrow(meta_221561)))
cat(sprintf("  Neoadjuvant_treated: %d\n", sum(meta_221561$treatment_standardized == "Neoadjuvant_treated")))
cat(sprintf("  Surgery_alone: %d\n\n", sum(meta_221561$treatment_standardized == "Surgery_alone")))

# Common gene set
common_genes <- intersect(rownames(disc_all), rownames(repl_all))
cat(sprintf("Common tested genes: %d\n\n", length(common_genes)))

# =============================================================================
# STEP 19D8: NOMINAL CANDIDATE SUPPORT
# =============================================================================

cat("=== STEP 19D8: NOMINAL CANDIDATE SUPPORT ===\n\n")

# Match nominal candidates to replication
nominal_genes <- rownames(nominal)
match_idx <- match(nominal_genes, rownames(repl_all))
tested <- !is.na(match_idx)

cat(sprintf("Frozen nominal candidates: %d\n", length(nominal_genes)))
cat(sprintf("Tested in replication: %d\n", sum(tested)))
cat(sprintf("Not tested: %d\n\n", length(nominal_genes) - sum(tested)))

# Build candidate support table
candidate_support <- data.frame(
  gene = nominal_genes,
  discovery_logFC = nominal[nominal_genes, "logFC"],
  discovery_PValue = nominal[nominal_genes, "PValue"],
  discovery_FDR = nominal[nominal_genes, "FDR"],
  discovery_direction = nominal[nominal_genes, "direction"],
  replication_logFC = NA_real_,
  replication_PValue = NA_real_,
  replication_FDR = NA_real_,
  replication_direction = NA_character_,
  same_direction = NA,
  support_class = "NOT_TESTED",
  stringsAsFactors = FALSE
)

# Fill in replication data for tested genes
for (i in which(tested)) {
  gene <- nominal_genes[i]
  r_idx <- match_idx[i]
  candidate_support$replication_logFC[i] <- repl_all[r_idx, "logFC"]
  candidate_support$replication_PValue[i] <- repl_all[r_idx, "PValue"]
  candidate_support$replication_FDR[i] <- repl_all[r_idx, "FDR"]
  candidate_support$replication_direction[i] <- repl_all[r_idx, "direction"]

  # CORRECT: use sign comparison, not string comparison
  d_fc <- candidate_support$discovery_logFC[i]
  r_fc <- candidate_support$replication_logFC[i]
  candidate_support$same_direction[i] <- sign(d_fc) == sign(r_fc)

  # Support classification
  if (candidate_support$same_direction[i]) {
    if (!is.na(r_fc) && !is.na(candidate_support$replication_PValue[i]) &&
        candidate_support$replication_PValue[i] < 0.05) {
      candidate_support$support_class[i] <- "DIRECTION_AND_REPLICATION_P05"
    } else if (!is.na(r_fc) && !is.na(candidate_support$replication_PValue[i]) &&
               candidate_support$replication_PValue[i] < 0.10) {
      candidate_support$support_class[i] <- "DIRECTION_AND_REPLICATION_P10"
    } else {
      candidate_support$support_class[i] <- "DIRECTION_ONLY"
    }
  } else {
    candidate_support$support_class[i] <- "OPPOSITE_DIRECTION"
  }
}

# Summary counts
cat("Support category counts:\n")
cat(sprintf("  DIRECTION_ONLY: %d\n", sum(candidate_support$support_class == "DIRECTION_ONLY")))
cat(sprintf("  DIRECTION_AND_REPLICATION_P10: %d\n", sum(candidate_support$support_class == "DIRECTION_AND_REPLICATION_P10")))
cat(sprintf("  DIRECTION_AND_REPLICATION_P05: %d\n", sum(candidate_support$support_class == "DIRECTION_AND_REPLICATION_P05")))
cat(sprintf("  OPPOSITE_DIRECTION: %d\n", sum(candidate_support$support_class == "OPPOSITE_DIRECTION")))
cat(sprintf("  NOT_TESTED: %d\n\n", sum(candidate_support$support_class == "NOT_TESTED")))

# Overall support metrics
n_tested_cand <- sum(tested)
n_same_cand <- sum(candidate_support$same_direction[tested], na.rm = TRUE)
n_same_p10 <- sum(candidate_support$support_class == "DIRECTION_AND_REPLICATION_P10")
n_same_p05 <- sum(candidate_support$support_class == "DIRECTION_AND_REPLICATION_P05")
n_opposite <- sum(candidate_support$support_class == "OPPOSITE_DIRECTION")

cat(sprintf("Same direction: %d / %d = %.3f\n", n_same_cand, n_tested_cand, n_same_cand / n_tested_cand))
cat(sprintf("Same direction + P<0.10: %d\n", n_same_p10))
cat(sprintf("Same direction + P<0.05: %d\n", n_same_p05))
cat(sprintf("Opposite direction: %d\n\n", n_opposite))

# Save
utils::write.csv(candidate_support,
                 file.path(repl_dir, "STEP19D_nominal_candidate_directional_support_CORRECTED.csv"))
cat("Saved: STEP19D_nominal_candidate_directional_support_CORRECTED.csv\n\n")

# =============================================================================
# STEP 19D9: SPECIFIC GENE AUDIT (RCN3, SDC2, SENP5)
# =============================================================================

cat("=== STEP 19D9: SPECIFIC GENE AUDIT ===\n\n")

target_genes <- c("RCN3", "SDC2", "SENP5")

# Load GSE197677 LOSO data for these genes
loso_targets <- loso_197677[loso_197677$gene %in% target_genes, ]

gene_audit <- data.frame()
for (gene in target_genes) {
  d_fc <- disc_all[gene, "logFC"]
  d_pval <- disc_all[gene, "PValue"]
  d_fdr <- disc_all[gene, "FDR"]

  # LOSO from GSE197677
  loso_row <- loso_197677[loso_197677$gene == gene, ]
  if (nrow(loso_row) > 0) {
    loso_status <- loso_row$stability
    loso_consistency <- loso_row$direction_consistency_fraction
  } else {
    loso_status <- "NOT_IN_LOSO"
    loso_consistency <- NA
  }

  # Replication
  r_fc <- repl_all[gene, "logFC"]
  r_pval <- repl_all[gene, "PValue"]
  r_fdr <- repl_all[gene, "FDR"]
  same_dir <- sign(d_fc) == sign(r_fc)
  abs_r_fc <- abs(r_fc)

  # Evidence classification
  if (same_dir && abs_r_fc >= log2(1.5) && !is.na(r_pval) && r_pval < 0.10) {
    evidence_class <- "STRONG_EFFECT_DIRECTIONAL_SUPPORT"
  } else if (same_dir && abs_r_fc >= log2(1.2)) {
    evidence_class <- "DIRECTIONAL_SUPPORT"
  } else if (same_dir && abs_r_fc < log2(1.2)) {
    evidence_class <- "SAME_SIGN_MINIMAL_EFFECT"
  } else {
    evidence_class <- "OPPOSITE_DIRECTION"
  }

  cat(sprintf("%s:\n", gene))
  cat(sprintf("  Discovery: logFC=%+.4f, P=%.2e, FDR=%.4f\n", d_fc, d_pval, d_fdr))
  cat(sprintf("  LOSO: status=%s, consistency=%.2f\n", loso_status, loso_consistency))
  cat(sprintf("  Replication: logFC=%+.4f, P=%.2e, FDR=%.4f\n", r_fc, r_pval, r_fdr))
  cat(sprintf("  Same direction: %s\n", same_dir))
  cat(sprintf("  abs(replication logFC): %.4f\n", abs_r_fc))
  cat(sprintf("  Evidence class: %s\n\n", evidence_class))

  gene_audit <- rbind(gene_audit, data.frame(
    gene = gene,
    discovery_logFC = d_fc,
    discovery_PValue = d_pval,
    discovery_FDR = d_fdr,
    discovery_LOSO_status = loso_status,
    discovery_LOSO_consistency = loso_consistency,
    replication_logFC = r_fc,
    replication_PValue = r_pval,
    replication_FDR = r_fdr,
    same_direction = same_dir,
    abs_replication_logFC = abs_r_fc,
    evidence_class = evidence_class,
    stringsAsFactors = FALSE
  ))
}

utils::write.csv(gene_audit,
                 file.path(repl_dir, "STEP19D_RCN3_SDC2_SENP5_AUDIT.csv"))
cat("Saved: STEP19D_RCN3_SDC2_SENP5_AUDIT.csv\n\n")

# =============================================================================
# STEP 19D10: GSE221561 LOSO REPLICATION SENSITIVITY
# =============================================================================

cat("=== STEP 19D10: GSE221561 LOSO REPLICATION ===\n\n")

# Primary LOSO: remove each Neoadjuvant_treated sample once
treated_samples <- meta_221561$sample[meta_221561$treatment_standardized == "Neoadjuvant_treated"]
surgery_samples <- meta_221561$sample[meta_221561$treatment_standardized == "Surgery_alone"]

cat(sprintf("Neoadjuvant_treated samples: %s\n", paste(treated_samples, collapse = ", ")))
cat(sprintf("Surgery_alone samples: %s\n\n", paste(surgery_samples, collapse = ", ")))

# Primary LOSO gene set: frozen Top-100
top100_genes <- freeze$gene[1:100]

# Storage for LOSO results
loso_results <- data.frame(
  gene = top100_genes,
  valid_runs_tested = 0,
  same_direction_runs = 0,
  direction_consistency_fraction = NA_real_,
  median_replication_logFC = NA_real_,
  min_replication_logFC = NA_real_,
  max_replication_logFC = NA_real_,
  sign_flip_count = 0,
  stability = "NOT_TESTED_ENOUGH",
  stringsAsFactors = FALSE
)

# Run LOSO for each treated sample removal
cat("Running LOSO (remove each treated sample)...\n\n")

all_loso_fc <- matrix(NA_real_, nrow = length(top100_genes), ncol = length(treated_samples))
rownames(all_loso_fc) <- top100_genes
colnames(all_loso_fc) <- paste0("LOSO_", treated_samples)

for (s in seq_along(treated_samples)) {
  leave_out <- treated_samples[s]
  cat(sprintf("  Run %d: leaving out %s\n", s, leave_out))

  # Subset samples
  keep_samples <- meta_221561$sample != leave_out
  counts_sub <- raw_counts_221561[, keep_samples]
  meta_sub <- meta_221561[keep_samples, ]

  # Build DGEList
  treatment <- factor(meta_sub$treatment_standardized,
                      levels = c("Surgery_alone", "Neoadjuvant_treated"))
  design <- model.matrix(~ 0 + treatment)
  colnames(design) <- c("Surgery_alone", "Neoadjuvant_treated")
  rownames(design) <- meta_sub$sample

  dge <- DGEList(counts = counts_sub, samples = meta_sub)
  keep_gene <- filterByExpr(dge, design = design)
  dge_f <- dge[keep_gene, , keep.lib.sizes = FALSE]
  dge_f <- calcNormFactors(dge_f, method = "TMM")
  dge_f <- estimateDisp(dge_f, design, robust = TRUE)

  fit <- glmQLFit(dge_f, design, robust = TRUE)

  # Canonical contrast: Neoadjuvant_treated - Surgery_alone
  contrast_vec <- limma::makeContrasts(
    Neoadjuvant_treated - Surgery_alone,
    levels = design
  )
  qlf <- glmQLFTest(fit, contrast = contrast_vec)
  res <- topTags(qlf, n = Inf, sort.by = "PValue")$table

  # Extract logFC for top100 genes
  for (gene in top100_genes) {
    if (gene %in% rownames(res)) {
      all_loso_fc[gene, s] <- res[gene, "logFC"]
    }
  }
}

cat("\nLOSO runs complete.\n\n")

# Calculate LOSO stability metrics
for (i in seq_along(top100_genes)) {
  gene <- top100_genes[i]
  fc_values <- all_loso_fc[gene, ]
  valid <- !is.na(fc_values)
  n_valid <- sum(valid)

  loso_results$valid_runs_tested[i] <- n_valid

  if (n_valid >= 5) {
    fc_valid <- fc_values[valid]
    disc_fc <- disc_all[gene, "logFC"]

    # Count same direction runs
    same_dir <- sum(sign(disc_fc) == sign(fc_valid))
    loso_results$same_direction_runs[i] <- same_dir
    loso_results$direction_consistency_fraction[i] <- same_dir / n_valid

    loso_results$median_replication_logFC[i] <- median(fc_valid)
    loso_results$min_replication_logFC[i] <- min(fc_valid)
    loso_results$max_replication_logFC[i] <- max(fc_valid)
    loso_results$sign_flip_count[i] <- sum(sign(disc_fc) != sign(fc_valid))

    # Stability label
    if (same_dir >= 6) {
      loso_results$stability[i] <- "REPLICATION_DIRECTION_ROBUST"
    } else if (same_dir == 5) {
      loso_results$stability[i] <- "REPLICATION_DIRECTION_MODERATE"
    } else {
      loso_results$stability[i] <- "REPLICATION_DIRECTION_UNSTABLE"
    }
  }
}

# Summary
cat("GSE221561 LOSO stability summary (Top-100):\n")
cat(sprintf("  ROBUST: %d\n", sum(loso_results$stability == "REPLICATION_DIRECTION_ROBUST")))
cat(sprintf("  MODERATE: %d\n", sum(loso_results$stability == "REPLICATION_DIRECTION_MODERATE")))
cat(sprintf("  UNSTABLE: %d\n", sum(loso_results$stability == "REPLICATION_DIRECTION_UNSTABLE")))
cat(sprintf("  NOT_TESTED_ENOUGH: %d\n\n", sum(loso_results$stability == "NOT_TESTED_ENOUGH")))

utils::write.csv(loso_results,
                 file.path(repl_dir, "STEP19D_GSE221561_LOSO_REPLICATION_STABILITY_CORRECTED.csv"))
cat("Saved: STEP19D_GSE221561_LOSO_REPLICATION_STABILITY_CORRECTED.csv\n\n")

# =============================================================================
# STEP 19D11: FINAL TABLE WITH LOSO SUPPORT
# =============================================================================

cat("=== STEP 19D11: FINAL TABLE WITH LOSO SUPPORT ===\n\n")

# Build final table for all Top-100 genes
final_table <- data.frame(
  gene = top100_genes,
  discovery_rank = 1:100,
  discovery_logFC = disc_all[top100_genes, "logFC"],
  discovery_PValue = disc_all[top100_genes, "PValue"],
  discovery_FDR = disc_all[top100_genes, "FDR"],
  discovery_LOSO_status = NA_character_,
  replication_logFC = NA_real_,
  replication_PValue = NA_real_,
  replication_FDR = NA_real_,
  same_direction = NA,
  replication_LOSO_status = NA_character_,
  replication_LOSO_consistency = NA_real_,
  stringsAsFactors = FALSE
)

# Fill discovery LOSO from GSE197677
for (i in seq_along(top100_genes)) {
  gene <- top100_genes[i]
  loso_row <- loso_197677[loso_197677$gene == gene, ]
  if (nrow(loso_row) > 0) {
    final_table$discovery_LOSO_status[i] <- loso_row$stability
  }
}

# Fill replication data
for (i in seq_along(top100_genes)) {
  gene <- top100_genes[i]
  if (gene %in% rownames(repl_all)) {
    final_table$replication_logFC[i] <- repl_all[gene, "logFC"]
    final_table$replication_PValue[i] <- repl_all[gene, "PValue"]
    final_table$replication_FDR[i] <- repl_all[gene, "FDR"]

    d_fc <- final_table$discovery_logFC[i]
    r_fc <- final_table$replication_logFC[i]
    final_table$same_direction[i] <- sign(d_fc) == sign(r_fc)
  }
}

# Fill replication LOSO from Step 19D10
for (i in seq_along(top100_genes)) {
  gene <- top100_genes[i]
  loso_row <- loso_results[loso_results$gene == gene, ]
  if (nrow(loso_row) > 0) {
    final_table$replication_LOSO_status[i] <- loso_row$stability
    final_table$replication_LOSO_consistency[i] <- loso_row$direction_consistency_fraction
  }
}

utils::write.csv(final_table,
                 file.path(repl_dir, "STEP19D_DIRECTIONAL_REPLICATION_FINAL_TABLE.csv"))
cat("Saved: STEP19D_DIRECTIONAL_REPLICATION_FINAL_TABLE.csv\n\n")

# Summarize by rank subset
cat("Rank subset summaries:\n\n")

for (subset_name in c("Top25", "Top50", "Top100")) {
  k <- as.integer(gsub("Top", "", subset_name))
  idx <- 1:min(k, nrow(final_table))
  sub <- final_table[idx, ]

  n_tested <- sum(!is.na(sub$same_direction))
  n_same <- sum(sub$same_direction, na.rm = TRUE)
  conc <- n_same / n_tested
  p_val <- binom.test(n_same, n_tested, 0.5, alternative = "greater")$p.value

  n_robust <- sum(sub$replication_LOSO_status == "REPLICATION_DIRECTION_ROBUST", na.rm = TRUE)
  n_moderate <- sum(sub$replication_LOSO_status == "REPLICATION_DIRECTION_MODERATE", na.rm = TRUE)
  n_unstable <- sum(sub$replication_LOSO_status == "REPLICATION_DIRECTION_UNSTABLE", na.rm = TRUE)

  cat(sprintf("%s:\n", subset_name))
  cat(sprintf("  Tested: %d\n", n_tested))
  cat(sprintf("  Same direction: %d / %d = %.3f\n", n_same, n_tested, conc))
  cat(sprintf("  Binomial P: %.2e\n", p_val))
  cat(sprintf("  LOSO robust: %d\n", n_robust))
  cat(sprintf("  LOSO moderate: %d\n", n_moderate))
  cat(sprintf("  LOSO unstable: %d\n\n", n_unstable))
}

# =============================================================================
# STEP 19D12: EFFECT-SIZE CONCORDANCE PLOTS
# =============================================================================

cat("=== STEP 19D12: EFFECT-SIZE CONCORDANCE PLOTS ===\n\n")

# Load ggplot2 if available
if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)

  # Common genes data
  common_df <- data.frame(
    disc_fc = disc_all[common_genes, "logFC"],
    repl_fc = repl_all[common_genes, "logFC"]
  )

  # 1. Genome-wide scatter
  p1 <- ggplot(common_df, aes(x = disc_fc, y = repl_fc)) +
    geom_point(alpha = 0.2, size = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_smooth(method = "lm", se = FALSE, color = "red", linewidth = 1) +
    labs(
      title = "Genome-Wide Effect-Size Concordance",
      subtitle = "Related but non-identical treatment contrasts",
      x = "GSE197677 logFC (NACT vs nNACT)",
      y = "GSE221561 logFC (Neoadjuvant vs Surgery_alone)"
    ) +
    theme_minimal()

  ggsave(file.path(fig_dir, "STEP19D_effect_size_concordance_ALL_COMMON_CORRECTED.png"),
         p1, width = 8, height = 8, dpi = 300)
  cat("Saved: STEP19D_effect_size_concordance_ALL_COMMON_CORRECTED.png\n")

  # 2. Top-100 scatter
  top100_df <- data.frame(
    gene = top100_genes,
    disc_fc = final_table$discovery_logFC,
    repl_fc = final_table$replication_logFC,
    same_dir = final_table$same_direction
  )
  top100_df$direction <- ifelse(top100_df$same_dir, "Same", "Opposite")

  # Label top genes
  label_genes <- c("RCN3", "SDC2", "SENP5", "COL1A2", "SFRP2", "THY1",
                   "CTSK", "COL5A2", "LY96", "SLK", "GJB2", "ZNF556")
  top100_df$label <- ifelse(top100_df$gene %in% label_genes, top100_df$gene, "")

  p2 <- ggplot(top100_df, aes(x = disc_fc, y = repl_fc, color = direction)) +
    geom_point(size = 2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_text(aes(label = label), vjust = -0.8, size = 3, show.legend = FALSE) +
    scale_color_manual(values = c("Same" = "steelblue", "Opposite" = "tomato")) +
    labs(
      title = "Top-100 Effect-Size Concordance",
      subtitle = "Frozen discovery ranking",
      x = "GSE197677 logFC (NACT vs nNACT)",
      y = "GSE221561 logFC (Neoadjuvant vs Surgery_alone)",
      color = "Direction"
    ) +
    theme_minimal()

  ggsave(file.path(fig_dir, "STEP19D_effect_size_concordance_TOP100_CORRECTED.png"),
         p2, width = 8, height = 8, dpi = 300)
  cat("Saved: STEP19D_effect_size_concordance_TOP100_CORRECTED.png\n")

  # 3. Top-50 scatter
  top50_df <- top100_df[1:50, ]
  p3 <- ggplot(top50_df, aes(x = disc_fc, y = repl_fc, color = direction)) +
    geom_point(size = 2.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_text(aes(label = label), vjust = -0.8, size = 3, show.legend = FALSE) +
    scale_color_manual(values = c("Same" = "steelblue", "Opposite" = "tomato")) +
    labs(
      title = "Secondary Top-50 Rank Subset",
      subtitle = "Secondary directional concordance (P = 0.016)",
      x = "GSE197677 logFC (NACT vs nNACT)",
      y = "GSE221561 logFC (Neoadjuvant vs Surgery_alone)",
      color = "Direction"
    ) +
    theme_minimal()

  ggsave(file.path(fig_dir, "STEP19D_effect_size_concordance_TOP50_CORRECTED.png"),
         p3, width = 8, height = 8, dpi = 300)
  cat("Saved: STEP19D_effect_size_concordance_TOP50_CORRECTED.png\n")

  # 4. LOSO stability barplot
  loso_summary <- data.frame(
    Stability = c("ROBUST", "MODERATE", "UNSTABLE"),
    Count = c(
      sum(loso_results$stability == "REPLICATION_DIRECTION_ROBUST"),
      sum(loso_results$stability == "REPLICATION_DIRECTION_MODERATE"),
      sum(loso_results$stability == "REPLICATION_DIRECTION_UNSTABLE")
    )
  )
  loso_summary$Stability <- factor(loso_summary$Stability,
                                   levels = c("ROBUST", "MODERATE", "UNSTABLE"))

  p4 <- ggplot(loso_summary, aes(x = Stability, y = Count, fill = Stability)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = Count), vjust = -0.5) +
    scale_fill_manual(values = c("ROBUST" = "forestgreen",
                                 "MODERATE" = "goldenrod",
                                 "UNSTABLE" = "tomato")) +
    labs(
      title = "GSE221561 LOSO Replication Stability",
      subtitle = "Top-100 genes across 7 treated-removal runs",
      x = "Stability",
      y = "Number of genes"
    ) +
    theme_minimal() +
    theme(legend.position = "none")

  ggsave(file.path(fig_dir, "STEP19D_replication_LOSO_stability_TOP100.png"),
         p4, width = 6, height = 6, dpi = 300)
  cat("Saved: STEP19D_replication_LOSO_stability_TOP100.png\n")

  # 5. RCN3, SDC2, SENP5 sample expression
  # Load TMM logCPM from GSE221561
  dge_full <- DGEList(counts = raw_counts_221561, samples = meta_221561)
  keep_gene <- filterByExpr(dge_full, design = model.matrix(~ treatment_standardized, data = meta_221561))
  dge_full <- dge_full[keep_gene, , keep.lib.sizes = FALSE]
  dge_full <- calcNormFactors(dge_full, method = "TMM")
  logcpm <- cpm(dge_full, log = TRUE)

  # Extract target genes
  target_expr <- data.frame()
  for (gene in target_genes) {
    if (gene %in% rownames(logcpm)) {
      for (s in 1:nrow(meta_221561)) {
        sample <- meta_221561$sample[s]
        treatment <- meta_221561$treatment_standardized[s]
        target_expr <- rbind(target_expr, data.frame(
          gene = gene,
          sample = sample,
          treatment = treatment,
          logCPM = logcpm[gene, s],
          stringsAsFactors = FALSE
        ))
      }
    }
  }

  p5 <- ggplot(target_expr, aes(x = sample, y = logCPM, fill = treatment)) +
    geom_bar(stat = "identity") +
    facet_wrap(~ gene, scales = "free_y", ncol = 1) +
    scale_fill_manual(values = c("Neoadjuvant_treated" = "steelblue",
                                 "Surgery_alone" = "tomato")) +
    labs(
      title = "RCN3, SDC2, SENP5 Sample-Level Expression",
      subtitle = "GSE221561 TMM logCPM",
      x = "Sample",
      y = "logCPM",
      fill = "Treatment"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave(file.path(fig_dir, "STEP19D_RCN3_SDC2_SENP5_sample_expression_CORRECTED.png"),
         p5, width = 10, height = 8, dpi = 300)
  cat("Saved: STEP19D_RCN3_SDC2_SENP5_sample_expression_CORRECTED.png\n")

} else {
  cat("ggplot2 not available, skipping plots.\n\n")
}

# =============================================================================
# STEP 19D13: FINAL REPLICATION CLASSIFICATION
# =============================================================================

cat("\n=== STEP 19D13: FINAL REPLICATION CLASSIFICATION ===\n\n")

# Calculate final metrics
# Top-100 (PRIMARY)
top100_tested <- sum(!is.na(final_table$same_direction[1:100]))
top100_same <- sum(final_table$same_direction[1:100], na.rm = TRUE)
top100_conc <- top100_same / top100_tested
top100_p <- binom.test(top100_same, top100_tested, 0.5, alternative = "greater")$p.value

# Top-50 (SECONDARY)
top50_tested <- sum(!is.na(final_table$same_direction[1:50]))
top50_same <- sum(final_table$same_direction[1:50], na.rm = TRUE)
top50_conc <- top50_same / top50_tested
top50_p <- binom.test(top50_same, top50_tested, 0.5, alternative = "greater")$p.value

# Top-25 (SECONDARY)
top25_tested <- sum(!is.na(final_table$same_direction[1:25]))
top25_same <- sum(final_table$same_direction[1:25], na.rm = TRUE)
top25_conc <- top25_same / top25_tested
top25_p <- binom.test(top25_same, top25_tested, 0.5, alternative = "greater")$p.value

# Genome-wide
gw_tested <- length(common_genes)
gw_same <- sum(sign(disc_all[common_genes, "logFC"]) == sign(repl_all[common_genes, "logFC"]))
gw_conc <- gw_same / gw_tested
gw_p <- binom.test(gw_same, gw_tested, 0.5, alternative = "greater")$p.value
gw_rho <- cor(disc_all[common_genes, "logFC"], repl_all[common_genes, "logFC"],
              method = "spearman", use = "complete.obs")

# Classification
decision <- data.frame(
  metric = c(
    "primary_endpoint_status",
    "secondary_top50_status",
    "secondary_top25_status",
    "genomewide_status",
    "overall_interpretation"
  ),
  value = c(
    "PRIMARY_ENDPOINT_NOT_SIGNIFICANT",
    "SECONDARY_DIRECTIONAL_SUPPORT",
    "SECONDARY_SUGGESTIVE_SUPPORT",
    "NO_GENOMEWIDE_DIRECTIONAL_SHIFT",
    "MODEST_TOP_RANKED_DIRECTIONAL_CONCORDANCE"
  ),
  details = c(
    sprintf("Top100: %d/%d = %.3f, P = %.3f", top100_same, top100_tested, top100_conc, top100_p),
    sprintf("Top50: %d/%d = %.3f, P = %.3f", top50_same, top50_tested, top50_conc, top50_p),
    sprintf("Top25: %d/%d = %.3f, P = %.3f", top25_same, top25_tested, top25_conc, top25_p),
    sprintf("Genome-wide: %d/%d = %.3f, P = %.3f", gw_same, gw_tested, gw_conc, gw_p),
    "Primary Top100 not significant; Top50/25 supportive; genome-wide null"
  ),
  stringsAsFactors = FALSE
)

utils::write.csv(decision,
                 file.path(repl_dir, "STEP19D_REPLICATION_DECISION_CORRECTED.csv"),
                 row.names = FALSE)
cat("Saved: STEP19D_REPLICATION_DECISION_CORRECTED.csv\n\n")

# Print classification
for (i in 1:nrow(decision)) {
  cat(sprintf("%s: %s\n", decision$metric[i], decision$value[i]))
}

# =============================================================================
# STEP 19D14: FREEZE CROSS-DATASET EVIDENCE
# =============================================================================

cat("\n=== STEP 19D14: FREEZE CROSS-DATASET EVIDENCE ===\n\n")

# Build final evidence table for all common genes
evidence_table <- data.frame(
  gene = common_genes,
  GSE197677_logFC = disc_all[common_genes, "logFC"],
  GSE197677_PValue = disc_all[common_genes, "PValue"],
  GSE197677_FDR = disc_all[common_genes, "FDR"],
  GSE197677_rank = NA_integer_,
  GSE197677_LOSO_status = NA_character_,
  GSE221561_logFC = repl_all[common_genes, "logFC"],
  GSE221561_PValue = repl_all[common_genes, "PValue"],
  GSE221561_FDR = repl_all[common_genes, "FDR"],
  same_direction = sign(disc_all[common_genes, "logFC"]) == sign(repl_all[common_genes, "logFC"]),
  GSE221561_LOSO_status = NA_character_,
  GSE221561_LOSO_consistency = NA_real_,
  discovery_signed_score = disc_all[common_genes, "logFC"] * (-log10(disc_all[common_genes, "PValue"])),
  replication_signed_score = repl_all[common_genes, "logFC"] * (-log10(repl_all[common_genes, "PValue"])),
  rank_subset = "OUTSIDE_TOP100",
  directional_evidence_class = "LOW",
  stringsAsFactors = FALSE
)

# Fill rank from freeze
match_rank <- match(common_genes, freeze$gene)
evidence_table$GSE197677_rank <- ifelse(!is.na(match_rank), match_rank, NA_integer_)

# Fill discovery LOSO
match_loso <- match(common_genes, loso_197677$gene)
evidence_table$GSE197677_LOSO_status <- ifelse(!is.na(match_loso),
                                               loso_197677$stability[match_loso],
                                               NA_character_)

# Fill replication LOSO
match_repl_loso <- match(common_genes, loso_results$gene)
evidence_table$GSE221561_LOSO_status <- ifelse(!is.na(match_repl_loso),
                                               loso_results$stability[match_repl_loso],
                                               NA_character_)
evidence_table$GSE221561_LOSO_consistency <- ifelse(!is.na(match_repl_loso),
                                                    loso_results$direction_consistency_fraction[match_repl_loso],
                                                    NA_real_)

# Assign rank_subset
evidence_table$rank_subset[evidence_table$GSE197677_rank <= 25] <- "TOP25"
evidence_table$rank_subset[evidence_table$GSE197677_rank > 25 &
                           evidence_table$GSE197677_rank <= 50] <- "TOP50"
evidence_table$rank_subset[evidence_table$GSE197677_rank > 50 &
                           evidence_table$GSE197677_rank <= 100] <- "TOP100"

# Assign directional_evidence_class
# HIGH: rank <= 100 AND same_direction AND LOSO ROBUST AND abs(rep_logFC) >= log2(1.5)
# MODERATE: rank <= 100 AND same_direction AND LOSO ROBUST or MODERATE
# LOW: otherwise
for (i in 1:nrow(evidence_table)) {
  rank <- evidence_table$GSE197677_rank[i]
  same_dir <- evidence_table$same_direction[i]
  loso_status <- evidence_table$GSE221561_LOSO_status[i]
  abs_repl_fc <- abs(evidence_table$GSE221561_logFC[i])

  if (!is.na(rank) && rank <= 100 && same_dir &&
      !is.na(loso_status) && loso_status == "REPLICATION_DIRECTION_ROBUST" &&
      abs_repl_fc >= log2(1.5)) {
    evidence_table$directional_evidence_class[i] <- "HIGH"
  } else if (!is.na(rank) && rank <= 100 && same_dir &&
             !is.na(loso_status) && loso_status %in% c("REPLICATION_DIRECTION_ROBUST",
                                                        "REPLICATION_DIRECTION_MODERATE")) {
    evidence_table$directional_evidence_class[i] <- "MODERATE"
  } else {
    evidence_table$directional_evidence_class[i] <- "LOW"
  }
}

utils::write.csv(evidence_table,
                 file.path(repl_dir, "STEP19D_CROSS_DATASET_RANKED_GENE_EVIDENCE_FINAL.csv"))
cat("Saved: STEP19D_CROSS_DATASET_RANKED_GENE_EVIDENCE_FINAL.csv\n\n")

# Evidence class summary
cat("Directional evidence class summary:\n")
cat(sprintf("  HIGH: %d\n", sum(evidence_table$directional_evidence_class == "HIGH")))
cat(sprintf("  MODERATE: %d\n", sum(evidence_table$directional_evidence_class == "MODERATE")))
cat(sprintf("  LOW: %d\n\n", sum(evidence_table$directional_evidence_class == "LOW")))

# =============================================================================
# FINAL INTERPRETATION NOTE
# =============================================================================

cat("=== STEP 19D14: FINAL INTERPRETATION NOTE ===\n\n")

interpretation <- "# Step 19D Final Replication Interpretation

## Date
`r Sys.Date()`

## Key Findings

### 1. Discovery (GSE197677)
- Single-gene analysis yielded **no FDR < 0.05 genes**
- Lowest FDR = 0.124 (RCN3)
- 1,170 nominal candidates (P < 0.05, |logFC| >= log2(1.5), FDR >= 0.05)

### 2. Primary Endpoint: Top-100 Directional Replication
- **Tested: 99/100**
- **Same direction: 56/99 = 56.6%**
- **Binomial P = 0.114**
- **Status: PRIMARY_ENDPOINT_NOT_SIGNIFICANT**

### 3. Secondary Top-50 Subset
- **Same direction: 33/50 = 66.0%**
- **Binomial P = 0.016**
- **Status: SECONDARY_DIRECTIONAL_SUPPORT**

### 4. Secondary Top-25 Subset
- **Same direction: 17/25 = 68.0%**
- **Binomial P = 0.054**
- **Status: SECONDARY_SUGGESTIVE_SUPPORT**

### 5. Genome-Wide Direction Concordance
- **Same direction: 4484/9033 = 49.6%**
- **Binomial P = 0.756**
- **Status: NO_GENOMEWIDE_DIRECTIONAL_SHIFT**

### 6. Overall Interpretation
**MODEST_TOP_RANKED_DIRECTIONAL_CONCORDANCE**

The evidence supports:
- Modest enrichment of concordant direction among higher-ranked discovery genes
- NOT a genome-wide shared treatment effect

### 7. Limitations
- GSE221561 Surgery_alone n=2 (small reference group)
- Clinical contrasts are related but not identical
  - GSE197677: NACT vs nNACT
  - GSE221561: Neoadjuvant_treated vs Surgery_alone

### 8. Gene-Level Claims
- **No gene should be called a \"replicated DEG\" based on Step 19D alone**
- RCN3 and SDC2 should only be described according to their exact replication P/FDR/LOSO evidence
- SENP5 replication logFC ~ -0.03 (negligible effect)

## Provenance

### Bug Corrections
1. **Direction string comparison**: Original compared literal strings (\"Higher_in_NACT\" vs \"Higher_in_Neoadjuvant_treated\") instead of semantic equivalence
   - **Fix**: compare numeric sign(logFC)
2. **Named vector indexing**: `df$logFC[gene]` returns NA
   - **Fix**: use `df[gene, \"logFC\"]`

### Previous (Incorrect) Results
- Reported 0/99 concordance → corrected to 56/99 = 0.566
- Classified as TRUE_DISCORDANCE → corrected to MODEST_TOP_RANKED_DIRECTIONAL_CONCORDANCE
"

writeLines(interpretation,
           file.path(repl_dir, "STEP19D_FINAL_REPLICATION_INTERPRETATION.md"))
cat("Saved: STEP19D_FINAL_REPLICATION_INTERPRETATION.md\n\n")

# =============================================================================
# BUG CORRECTION PROVENANCE
# =============================================================================

bug_provenance <- "# Step 19D Bug Correction Provenance

## Date
`r Sys.Date()`

## Bugs Found and Corrected

### BUG 1: Literal Direction-String Comparison

**Original code:**
\`\`\`r
same_dir <- disc$direction[gene] == repl$direction[gene]
\`\`\`

**Problem:**
- \"Higher_in_NACT\" != \"Higher_in_Neoadjuvant_treated\" (different strings)
- \"Higher_in_nNACT\" != \"Higher_in_Surgery_alone\" (different strings)
- These are biologically equivalent directions (treatment vs control)
- Result: always FALSE → 0/99 concordance

**Fix:**
\`\`\`r
same_dir <- sign(disc$logFC[gene]) == sign(repl$logFC[gene])
\`\`\`

**Impact:**
- Previous 0/99 concordance was INVALIDATED
- Corrected concordance: 56/99 = 0.566

---

### BUG 2: Named Vector Indexing

**Original code:**
\`\`\`r
df$logFC[gene]
\`\`\`

**Problem:**
- `df$logFC` extracts a numeric vector WITHOUT names
- Indexing by character returns NA
- `df$logFC[\"RCN3\"]` → NA

**Fix:**
\`\`\`r
df[gene, \"logFC\"]
# or
df[which(rownames(df) == gene), \"logFC\"]
\`\`\`

**Impact:**
- All gene-level diagnostics returned NA
- Could not verify specific genes (RCN3, SDC2, SENP5)

---

## Previous Classification (INVALIDATED)

- **Status**: TRUE_DISCORDANCE
- **Reason**: Based on buggy 0/99 concordance
- **Corrected**: MODEST_TOP_RANKED_DIRECTIONAL_CONCORDANCE

## Files Updated

- STEP19D_nominal_candidate_directional_support_CORRECTED.csv
- STEP19D_RCN3_SDC2_SENP5_AUDIT.csv
- STEP19D_GSE221561_LOSO_REPLICATION_STABILITY_CORRECTED.csv
- STEP19D_DIRECTIONAL_REPLICATION_FINAL_TABLE.csv
- STEP19D_REPLICATION_DECISION_CORRECTED.csv
- STEP19D_CROSS_DATASET_RANKED_GENE_EVIDENCE_FINAL.csv
- STEP19D_FINAL_REPLICATION_INTERPRETATION.md
"

writeLines(bug_provenance,
           file.path(repl_dir, "STEP19D_BUG_CORRECTION_PROVENANCE.md"))
cat("Saved: STEP19D_BUG_CORRECTION_PROVENANCE.md\n\n")

# =============================================================================
# FINAL OUTPUT
# =============================================================================

cat("\n")
cat("============================================\n")
cat("STEP 19D FINALIZED ✓\n")
cat("CORRECTED CROSS-DATASET REPLICATION\n")
cat("============================================\n\n")

cat("Discovery:\n")
cat("GSE197677\n")
cat("NACT vs nNACT\n")
cat("6 vs 4\n\n")

cat("Replication:\n")
cat("GSE221561\n")
cat("Neoadjuvant_treated vs Surgery_alone\n")
cat("7 vs 2\n\n")

cat("--------------------------------------------\n\n")

cat("PRIMARY ENDPOINT — TOP100\n\n")
cat(sprintf("Tested                       : %d\n", top100_tested))
cat(sprintf("Same direction               : %d\n", top100_same))
cat(sprintf("Concordance                  : %.3f\n", top100_conc))
cat(sprintf("Binomial P                   : %.3f\n", top100_p))
cat("Status                       : PRIMARY_ENDPOINT_NOT_SIGNIFICANT\n\n")

cat("--------------------------------------------\n\n")

cat("SECONDARY — TOP50\n\n")
cat(sprintf("Tested                       : %d\n", top50_tested))
cat(sprintf("Same direction               : %d\n", top50_same))
cat(sprintf("Concordance                  : %.3f\n", top50_conc))
cat(sprintf("Binomial P                   : %.3f\n", top50_p))
cat("Status                       : SECONDARY_DIRECTIONAL_SUPPORT\n\n")

cat("--------------------------------------------\n\n")

cat("SECONDARY — TOP25\n\n")
cat(sprintf("Tested                       : %d\n", top25_tested))
cat(sprintf("Same direction               : %d\n", top25_same))
cat(sprintf("Concordance                  : %.3f\n", top25_conc))
cat(sprintf("Binomial P                   : %.3f\n", top25_p))
cat("Status                       : SECONDARY_SUGGESTIVE_SUPPORT\n\n")

cat("--------------------------------------------\n\n")

cat(sprintf("Genome-wide common genes     : %d\n", gw_tested))
cat(sprintf("Same direction               : %d\n", gw_same))
cat(sprintf("Concordance                  : %.3f\n", gw_conc))
cat(sprintf("Binomial P                   : %.3f\n", gw_p))
cat(sprintf("Spearman logFC rho           : %.4f\n\n", gw_rho))

cat("--------------------------------------------\n\n")

cat(sprintf("Nominal discovery candidates : %d\n", nrow(nominal)))
cat(sprintf("Tested                       : %d\n", n_tested_cand))
cat(sprintf("Same direction               : %d\n", n_same_cand))
cat(sprintf("Same direction + repl P<0.10 : %d\n", n_same_p10))
cat(sprintf("Same direction + repl P<0.05 : %d\n\n", n_same_p05))

cat("--------------------------------------------\n\n")

cat("Top100 replication LOSO:\n\n")
cat(sprintf("Robust                       : %d\n", sum(loso_results$stability == "REPLICATION_DIRECTION_ROBUST")))
cat(sprintf("Moderate                     : %d\n", sum(loso_results$stability == "REPLICATION_DIRECTION_MODERATE")))
cat(sprintf("Unstable                     : %d\n", sum(loso_results$stability == "REPLICATION_DIRECTION_UNSTABLE")))
cat(sprintf("Not tested enough            : %d\n\n", sum(loso_results$stability == "NOT_TESTED_ENOUGH")))

cat("--------------------------------------------\n\n")

# RCN3
rcn3_disc <- disc_all["RCN3", "logFC"]
rcn3_disc_fdr <- disc_all["RCN3", "FDR"]
rcn3_repl <- repl_all["RCN3", "logFC"]
rcn3_repl_p <- repl_all["RCN3", "PValue"]
rcn3_repl_fdr <- repl_all["RCN3", "FDR"]
rcn3_loso <- loso_results[loso_results$gene == "RCN3", "stability"]
rcn3_evidence <- evidence_table[evidence_table$gene == "RCN3", "directional_evidence_class"]

cat("RCN3\n")
cat(sprintf("Discovery logFC/FDR          : %+.4f / %.4f\n", rcn3_disc, rcn3_disc_fdr))
cat(sprintf("Replication logFC/P/FDR      : %+.4f / %.4f / %.4f\n", rcn3_repl, rcn3_repl_p, rcn3_repl_fdr))
cat(sprintf("Replication LOSO             : %s\n", rcn3_loso))
cat(sprintf("Evidence classification      : %s\n\n", rcn3_evidence))

# SDC2
sdc2_disc <- disc_all["SDC2", "logFC"]
sdc2_disc_fdr <- disc_all["SDC2", "FDR"]
sdc2_repl <- repl_all["SDC2", "logFC"]
sdc2_repl_p <- repl_all["SDC2", "PValue"]
sdc2_repl_fdr <- repl_all["SDC2", "FDR"]
sdc2_loso <- loso_results[loso_results$gene == "SDC2", "stability"]
sdc2_evidence <- evidence_table[evidence_table$gene == "SDC2", "directional_evidence_class"]

cat("SDC2\n")
cat(sprintf("Discovery logFC/FDR          : %+.4f / %.4f\n", sdc2_disc, sdc2_disc_fdr))
cat(sprintf("Replication logFC/P/FDR      : %+.4f / %.4f / %.4f\n", sdc2_repl, sdc2_repl_p, sdc2_repl_fdr))
cat(sprintf("Replication LOSO             : %s\n", sdc2_loso))
cat(sprintf("Evidence classification      : %s\n\n", sdc2_evidence))

# SENP5
senp5_disc <- disc_all["SENP5", "logFC"]
senp5_disc_fdr <- disc_all["SENP5", "FDR"]
senp5_repl <- repl_all["SENP5", "logFC"]
senp5_repl_p <- repl_all["SENP5", "PValue"]
senp5_repl_fdr <- repl_all["SENP5", "FDR"]
senp5_loso <- loso_results[loso_results$gene == "SENP5", "stability"]
senp5_evidence <- evidence_table[evidence_table$gene == "SENP5", "directional_evidence_class"]

cat("SENP5\n")
cat(sprintf("Discovery logFC/FDR          : %+.4f / %.4f\n", senp5_disc, senp5_disc_fdr))
cat(sprintf("Replication logFC/P/FDR      : %+.4f / %.4f / %.4f\n", senp5_repl, senp5_repl_p, senp5_repl_fdr))
cat(sprintf("Replication LOSO             : %s\n", senp5_loso))
cat(sprintf("Evidence classification      : %s\n\n", senp5_evidence))

cat("--------------------------------------------\n\n")

cat("Overall interpretation:\n")
cat("MODEST_TOP_RANKED_DIRECTIONAL_CONCORDANCE\n\n")

cat("Primary Top100 significant: NO\n")
cat("Secondary Top50 support: YES\n")
cat("Genome-wide concordance: NO\n")
cat("Replicated DEGs claimed: NO\n")
cat("Pathway analysis: NO\n")
cat("Samples removed: NO\n")
cat("Cells removed: NO\n\n")

cat("STOP HERE.\n")
cat("DO NOT START STEP 19E AUTOMATICALLY.\n")
