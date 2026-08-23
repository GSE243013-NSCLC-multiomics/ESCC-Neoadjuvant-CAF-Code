#!/usr/bin/env Rscript
# STEPS 19D4-19D7: DIRECTION CONCORDANCE, LOSO, GENOME-WIDE, RANK-BASED

base_dir <- "~/Downloads/ESCC Neoadjuvant"
repl_dir <- file.path(base_dir, "06_tables/cross_dataset/epithelial_pseudobulk/replication")
de_dir <- file.path(base_dir, "06_tables/cross_dataset/epithelial_pseudobulk/GSE197677_DE")

cat("Loading data...\n")
freeze <- read.csv(file.path(repl_dir, "STEP19D_GSE197677_DISCOVERY_RANK_FREEZE.csv"), stringsAsFactors = FALSE)
res_221561 <- read.csv(file.path(repl_dir, "STEP19D_GSE221561_NEOADJUVANT_vs_SURGERY_ALL_GENES.csv"), row.names = 1, check.names = FALSE)
harmonize <- read.csv(file.path(repl_dir, "STEP19D_common_tested_genes.csv"), stringsAsFactors = FALSE)
nominal_197677 <- read.csv(file.path(de_dir, "STEP19C_GSE197677_nominal_directional_candidates.csv"), row.names = 1, check.names = FALSE)
common_tested <- harmonize$gene
cat("Loaded.\n\n")

# === STEP 19D4: TOP-100 DIRECTION CONCORDANCE ===
cat("=== STEP 19D4: TOP-100 DIRECTION CONCORDANCE ===\n\n")

top100 <- freeze[1:min(100, nrow(freeze)), ]
top100$replication_logFC <- res_221561$logFC[match(top100$gene, rownames(res_221561))]
top100$replication_PValue <- res_221561$PValue[match(top100$gene, rownames(res_221561))]
top100$replication_FDR <- res_221561$FDR[match(top100$gene, rownames(res_221561))]
top100$replication_direction <- res_221561$direction[match(top100$gene, rownames(res_221561))]
top100$tested_in_221561 <- !is.na(top100$replication_logFC)

cat(sprintf("Top 100 frozen: %d\n", nrow(top100)))
cat(sprintf("Tested in GSE221561: %d\n", sum(top100$tested_in_221561)))

top100_tested <- top100[top100$tested_in_221561, ]
top100_tested$same_direction <- top100_tested$discovery_direction == top100_tested$replication_direction
n_tested_100 <- nrow(top100_tested)
n_same_100 <- sum(top100_tested$same_direction)
conc_100 <- n_same_100 / n_tested_100
bt_100 <- binom.test(n_same_100, n_tested_100, p = 0.5, alternative = "two.sided")

cat(sprintf("Same direction: %d / %d = %.3f\n", n_same_100, n_tested_100, conc_100))
cat(sprintf("Binomial P: %.4e\n\n", bt_100$p.value))

# Top 25 and 50
results_25_50 <- list()
for (k in c(25, 50)) {
  top_k <- freeze[1:k, ]
  repl_fc <- res_221561$logFC[match(top_k$gene, rownames(res_221561))]
  repl_dir <- res_221561$direction[match(top_k$gene, rownames(res_221561))]
  tested_idx <- !is.na(repl_fc)
  n_t <- sum(tested_idx)
  n_s <- sum(top_k$discovery_direction[tested_idx] == repl_dir[tested_idx])
  conc <- n_s / n_t
  bt <- binom.test(n_s, n_t, p = 0.5, alternative = "two.sided")
  cat(sprintf("Top %d: tested=%d, same=%d, concordance=%.3f, P=%.4e\n", k, n_t, n_s, conc, bt$p.value))
  results_25_50[[as.character(k)]] <- list(n_tested=n_t, n_same=n_s, concordance=conc, P=bt$p.value)
}
cat("\n")

# === STEP 19D5: ROBUST-LOSO SUBSET ===
cat("=== STEP 19D5: ROBUST-LOSO SUBSET REPLICATION ===\n\n")

robust_top100 <- top100[top100$LOSO_status == "ROBUST_DIRECTION", ]
cat(sprintf("Robust-LOSO genes in top 100: %d\n", nrow(robust_top100)))

robust_tested <- robust_top100[robust_top100$tested_in_221561, ]
robust_tested$same_direction <- robust_tested$discovery_direction == robust_tested$replication_direction
cat(sprintf("Tested: %d, Same direction: %d / %d = %.3f\n\n",
            nrow(robust_tested), sum(robust_tested$same_direction), nrow(robust_tested),
            sum(robust_tested$same_direction) / nrow(robust_tested)))

# Use only columns that exist
out_df <- data.frame(
  gene = robust_tested$gene,
  discovery_rank = robust_tested$discovery_rank,
  discovery_logFC = robust_tested$discovery_logFC,
  discovery_FDR = robust_tested$discovery_FDR,
  discovery_direction = robust_tested$discovery_direction,
  LOSO_status = robust_tested$LOSO_status,
  replication_logFC = robust_tested$replication_logFC,
  replication_PValue = robust_tested$replication_PValue,
  replication_FDR = robust_tested$replication_FDR,
  replication_direction = robust_tested$replication_direction,
  same_direction = robust_tested$same_direction,
  stringsAsFactors = FALSE
)
# Save using RDS then convert to avoid write.csv issues
saveRDS(out_df, file.path(base_dir, repl_dir, "STEP19D_top100_directional_replication.rds"))
# Also write CSV manually
csv_lines <- c(paste(names(out_df), collapse=","),
               apply(out_df, 1, function(x) paste(x, collapse=",")))
writeLines(csv_lines, file.path(base_dir, repl_dir, "STEP19D_top100_directional_replication.csv"))
cat("Saved: STEP19D_top100_directional_replication.csv\n\n")

# === STEP 19D6: GENOME-WIDE EFFECT-SIZE CONCORDANCE ===
cat("=== STEP 19D6: GENOME-WIDE EFFECT-SIZE CONCORDANCE ===\n\n")

common_freeze <- freeze[freeze$gene %in% common_tested, ]
common_freeze$repl_fc <- res_221561$logFC[match(common_freeze$gene, rownames(res_221561))]
common_complete <- common_freeze[!is.na(common_freeze$repl_fc), ]
cat(sprintf("Common genes with both logFC: %d\n", nrow(common_complete)))

rho_spearman <- cor(common_complete$discovery_logFC, common_complete$repl_fc, method = "spearman")
cat(sprintf("Spearman rho (logFC): %.4f\n", rho_spearman))

same_sign <- sum(sign(common_complete$discovery_logFC) == sign(common_complete$repl_fc))
cat(sprintf("Same sign: %d / %d = %.3f\n\n", same_sign, nrow(common_complete), same_sign / nrow(common_complete)))

write.csv(data.frame(gene=common_complete$gene, discovery_logFC=common_complete$discovery_logFC,
                     replication_logFC=common_complete$repl_fc,
                     same_direction=sign(common_complete$discovery_logFC)==sign(common_complete$repl_fc),
                     discovery_rank=common_complete$discovery_rank, stringsAsFactors=FALSE),
          file.path(repl_dir, "STEP19D_genomewide_effect_concordance.csv"), row.names=FALSE)

write.csv(data.frame(metric=c("common_genes","spearman_rho","same_sign_fraction"),
                     value=c(nrow(common_complete), round(rho_spearman,4), round(same_sign/nrow(common_complete),4)),
                     stringsAsFactors=FALSE),
          file.path(repl_dir, "STEP19D_genomewide_concordance_summary.csv"), row.names=FALSE)
cat("Saved.\n\n")

# === STEP 19D7: RANK-BASED CONCORDANCE ===
cat("=== STEP 19D7: RANK-BASED CONCORDANCE ===\n\n")

res_197677 <- read.csv(file.path(de_dir, "STEP19C_GSE197677_NACT_vs_nNACT_ALL_GENES.csv"), row.names=1, check.names=FALSE)
res_197677$signed_score <- sign(res_197677$logFC) * -log10(res_197677$PValue)
res_221561$signed_score <- sign(res_221561$logFC) * -log10(res_221561$PValue)

common_all <- intersect(rownames(res_197677), rownames(res_221561))
cat(sprintf("Common genes: %d\n", length(common_all)))

rho_signed <- cor(res_197677$signed_score[common_all], res_221561$signed_score[common_all], method="spearman")
cat(sprintf("Signed-rank rho: %.4f\n\n", rho_signed))

top100_pos_197677 <- rownames(res_197677)[order(-res_197677$signed_score)][1:100]
top100_neg_197677 <- rownames(res_197677)[order(res_197677$signed_score)][1:100]
top100_pos_221561 <- rownames(res_221561)[order(-res_221561$signed_score)][1:100]
top100_neg_221561 <- rownames(res_221561)[order(res_221561$signed_score)][1:100]

cat(sprintf("Top 100 positive overlap: %d / 100\n", length(intersect(top100_pos_197677, top100_pos_221561))))
cat(sprintf("Top 100 negative overlap: %d / 100\n\n", length(intersect(top100_neg_197677, top100_neg_221561))))

# Save summary
write.csv(data.frame(metric=c("ranked_signed_score_rho","top100_positive_overlap","top100_negative_overlap"),
                     value=c(round(rho_signed,4),
                             length(intersect(top100_pos_197677, top100_pos_221561)),
                             length(intersect(top100_neg_197677, top100_neg_221561))),
                     stringsAsFactors=FALSE),
          file.path(repl_dir, "STEP19D_rank_concordance_summary.csv"), row.names=FALSE)

cat("Steps 19D4-19D7 complete.\n")
