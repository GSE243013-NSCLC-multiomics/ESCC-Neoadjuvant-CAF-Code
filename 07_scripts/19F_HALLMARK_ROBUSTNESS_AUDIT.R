#!/usr/bin/env Rscript
# STEP 19F: CROSS-DATASET HALLMARK ROBUSTNESS, REDUNDANCY, AND LEADING-EDGE AUDIT

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(fgsea)
  library(data.table)
  library(ggplot2)
})

base_dir <- "~/Downloads/ESCC Neoadjuvant"
de_dir <- file.path(base_dir, "06_tables/cross_dataset/epithelial_pseudobulk/GSE197677_DE")
repl_dir <- file.path(base_dir, "06_tables/cross_dataset/epithelial_pseudobulk/replication")
pw_dir <- file.path(base_dir, "06_tables/cross_dataset/pathway_analysis")
counts_dir <- file.path(base_dir, "06_tables/cross_dataset/epithelial_pseudobulk/counts")
fig_dir <- file.path(base_dir, "05_figures/cross_dataset/pathway_analysis/Step19F")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

cat("=== 19F0: FREEZE INPUT PROVENANCE ===\n\n")

freeze <- read.csv(file.path(pw_dir, "STEP19E_HALLMARK_CROSS_DATASET_FREEZE.csv"), stringsAsFactors = FALSE)
supported_pathways <- freeze$pathway[freeze$pathway_evidence_class == "CROSS_DATASET_FDR_CONSISTENT"]
cat(sprintf("Frozen supported pathways: %d\n", length(supported_pathways)))
for (pw in supported_pathways) cat(sprintf("  - %s\n", pw))

disc_all <- read.csv(file.path(de_dir, "STEP19C_GSE197677_NACT_vs_nNACT_ALL_GENES.csv"), row.names = 1, check.names = FALSE)
repl_all <- read.csv(file.path(repl_dir, "STEP19D_GSE221561_NEOADJUVANT_MINUS_SURGERY_CANONICAL.csv"), row.names = 1, check.names = FALSE)
gene_evidence <- read.csv(file.path(repl_dir, "STEP19D_CROSS_DATASET_RANKED_GENE_EVIDENCE_FINAL.csv"), stringsAsFactors = FALSE)
disc_ranks <- read.csv(file.path(pw_dir, "STEP19E_GSE197677_COMMON_GENE_RANKS.csv"), stringsAsFactors = FALSE)
repl_ranks <- read.csv(file.path(pw_dir, "STEP19E_GSE221561_COMMON_GENE_RANKS.csv"), stringsAsFactors = FALSE)
common_genes <- intersect(rownames(disc_all), rownames(repl_all))

raw_counts_disc <- readRDS(file.path(counts_dir, "GSE197677_epithelial_pseudobulk_raw_counts.rds"))
meta_disc <- read.csv(file.path(counts_dir, "GSE197677_epithelial_pseudobulk_sample_metadata.csv"), stringsAsFactors = FALSE)
raw_counts_repl <- readRDS(file.path(counts_dir, "GSE221561_epithelial_pseudobulk_raw_counts.rds"))
meta_repl <- read.csv(file.path(counts_dir, "GSE221561_epithelial_pseudobulk_sample_metadata.csv"), stringsAsFactors = FALSE)

cache_file <- file.path(Sys.getenv("ESCC_CAF_PROJECT_ROOT", unset = "."), "cache", "msigdb.2026.1.Hs.H.rds")
hallmark_df <- readRDS(cache_file)
hallmark_pathways <- split(hallmark_df$db_gene_symbol, hallmark_df$gs_name)
hallmark_pathways <- lapply(hallmark_pathways, unique)

cat(sprintf("Common genes: %d\n\n", length(common_genes)))

# === 19F1: FDR OVERLAP AUDIT ===
cat("=== 19F1: FDR-SIGNIFICANT PATHWAY OVERLAP AUDIT ===\n\n")
fgsea_disc <- read.csv(file.path(pw_dir, "STEP19E_GSE197677_HALLMARK_FGSEA_ALL.csv"), stringsAsFactors = FALSE)
fgsea_repl <- read.csv(file.path(pw_dir, "STEP19E_GSE221561_HALLMARK_FGSEA_ALL.csv"), stringsAsFactors = FALSE)

disc_fdr05 <- fgsea_disc$pathway[fgsea_disc$padj < 0.05]
repl_fdr05 <- fgsea_repl$pathway[fgsea_repl$padj < 0.05]
both_fdr05 <- intersect(disc_fdr05, repl_fdr05)
both_fdr05_same <- character()
both_fdr05_opp <- character()

for (pw in both_fdr05) {
  d_nes <- fgsea_disc$NES[fgsea_disc$pathway == pw]
  r_nes <- fgsea_repl$NES[fgsea_repl$pathway == pw]
  if (sign(d_nes) == sign(r_nes)) both_fdr05_same <- c(both_fdr05_same, pw) else both_fdr05_opp <- c(both_fdr05_opp, pw)
}

cat(sprintf("Discovery FDR05: %d\nReplication FDR05: %d\nBoth FDR05: %d\nBoth same dir: %d\nBoth opposite: %d\n\n",
            length(disc_fdr05), length(repl_fdr05), length(both_fdr05), length(both_fdr05_same), length(both_fdr05_opp)))

pop <- nrow(fgsea_disc)
p_hyper <- phyper(length(both_fdr05) - 1, length(disc_fdr05), pop - length(disc_fdr05), length(repl_fdr05), lower.tail = FALSE)
cat(sprintf("Hypergeom P: %.4e\n", p_hyper))

if (length(both_fdr05) > 0) {
  p_dir <- binom.test(length(both_fdr05_same), length(both_fdr05), 0.5, alternative = "greater")$p.value
  cat(sprintf("Same dir among both: %d/%d = %.3f (P=%.4f)\n\n", length(both_fdr05_same), length(both_fdr05), length(both_fdr05_same)/length(both_fdr05), p_dir))
}

# === 19F2: LEADING-EDGE REDUNDANCY ===
cat("=== 19F2: LEADING-EDGE REDUNDANCY AUDIT ===\n\n")
disc_le <- list()
repl_le <- list()

for (pw in supported_pathways) {
  d_le_str <- fgsea_disc$leadingEdge[fgsea_disc$pathway == pw]
  r_le_str <- fgsea_repl$leadingEdge[fgsea_repl$pathway == pw]
  disc_le[[pw]] <- if (!is.na(d_le_str) && nchar(as.character(d_le_str)) > 0) unlist(strsplit(as.character(d_le_str), ",")) else character(0)
  repl_le[[pw]] <- if (!is.na(r_le_str) && nchar(as.character(r_le_str)) > 0) unlist(strsplit(as.character(r_le_str), ",")) else character(0)
}

calc_jaccard <- function(s1, s2) { i <- length(intersect(s1, s2)); u <- length(union(s1, s2)); if (u == 0) return(0); i/u }

n_pw <- length(supported_pathways)
disc_jaccard <- matrix(0, n_pw, n_pw); rownames(disc_jaccard) <- colnames(disc_jaccard) <- supported_pathways
repl_jaccard <- matrix(0, n_pw, n_pw); rownames(repl_jaccard) <- colnames(repl_jaccard) <- supported_pathways
for (i in 1:n_pw) for (j in 1:n_pw) {
  disc_jaccard[i,j] <- calc_jaccard(disc_le[[supported_pathways[i]]], disc_le[[supported_pathways[j]]])
  repl_jaccard[i,j] <- calc_jaccard(repl_le[[supported_pathways[i]]], repl_le[[supported_pathways[j]]])
}

cat("High overlap pairs (Jaccard >= 0.25):\n")
for (i in 1:(n_pw-1)) for (j in (i+1):n_pw) {
  if (disc_jaccard[i,j] >= 0.25 || repl_jaccard[i,j] >= 0.25)
    cat(sprintf("  %s vs %s: Disc=%.3f, Repl=%.3f\n", supported_pathways[i], supported_pathways[j], disc_jaccard[i,j], repl_jaccard[i,j]))
}
cat("\n")

# === 19F3: CROSS-DATASET LE REPRODUCIBILITY ===
cat("=== 19F3: CROSS-DATASET LE REPRODUCIBILITY ===\n\n")
le_repro <- data.frame()
for (pw in supported_pathways) {
  d_le <- disc_le[[pw]]; r_le <- repl_le[[pw]]
  overlap <- intersect(d_le, r_le); union_g <- union(d_le, r_le)
  d_nes <- fgsea_disc$NES[fgsea_disc$pathway == pw]; r_nes <- fgsea_repl$NES[fgsea_repl$pathway == pw]
  le_repro <- rbind(le_repro, data.frame(pathway=pw, disc_LE_n=length(d_le), repl_LE_n=length(r_le),
    overlap_n=length(overlap), union_n=length(union_g), jaccard=length(overlap)/length(union_g),
    disc_NES=d_nes, repl_NES=r_nes, same_dir=sign(d_nes)==sign(r_nes),
    frac_d_to_r=if(length(d_le)>0) length(overlap)/length(d_le) else 0,
    frac_r_to_d=if(length(r_le)>0) length(overlap)/length(r_le) else 0, stringsAsFactors=FALSE))
}
for (i in 1:nrow(le_repro)) cat(sprintf("  %s: Disc=%d, Repl=%d, Overlap=%d, J=%.3f\n", le_repro$pathway[i], le_repro$disc_LE_n[i], le_repro$repl_LE_n[i], le_repro$overlap_n[i], le_repro$jaccard[i]))
cat("\n")

# === 19F4: CORE LEADING-EDGE GENES ===
cat("=== 19F4: CORE RECURRENT LEADING-EDGE GENES ===\n\n")
all_le_genes <- unique(c(unlist(disc_le), unlist(repl_le)))
core_genes <- data.frame(gene=all_le_genes, n_disc_LE=0, n_repl_LE=0, n_both_LE=0, stringsAsFactors=FALSE)
for (pw in supported_pathways) {
  for (g in disc_le[[pw]]) if (g %in% core_genes$gene) core_genes$n_disc_LE[core_genes$gene==g] <- core_genes$n_disc_LE[core_genes$gene==g]+1
  for (g in repl_le[[pw]]) if (g %in% core_genes$gene) core_genes$n_repl_LE[core_genes$gene==g] <- core_genes$n_repl_LE[core_genes$gene==g]+1
  for (g in intersect(disc_le[[pw]], repl_le[[pw]])) if (g %in% core_genes$gene) core_genes$n_both_LE[core_genes$gene==g] <- core_genes$n_both_LE[core_genes$gene==g]+1
}

for (i in 1:nrow(core_genes)) {
  ge <- gene_evidence[gene_evidence$gene==core_genes$gene[i],]
  if (nrow(ge)>0) { core_genes$disc_logFC[i] <- ge$GSE197677_logFC; core_genes$repl_logFC[i] <- ge$GSE221561_logFC; core_genes$same_dir[i] <- ge$same_direction; core_genes$evidence[i] <- ge$directional_evidence_class }
}

core_genes <- core_genes[order(-core_genes$n_both_LE, -core_genes$n_disc_LE),]
cat("Top core genes (both LE >= 1):\n")
for (i in 1:min(20, sum(core_genes$n_both_LE>=1))) {
  g <- core_genes[i,]
  cat(sprintf("  %s: Disc=%d, Repl=%d, Both=%d\n", g$gene, g$n_disc_LE, g$n_repl_LE, g$n_both_LE))
}

cat("\nTarget gene presence:\n")
for (t in c("RCN3","SDC2","SENP5")) {
  d_pws <- sapply(supported_pathways, function(pw) t %in% disc_le[[pw]])
  r_pws <- sapply(supported_pathways, function(pw) t %in% repl_le[[pw]])
  cat(sprintf("%s: Disc LE in %d pathways (%s), Repl LE in %d (%s)\n", t, sum(d_pws), paste(names(d_pws)[d_pws], collapse=", "), sum(r_pws), paste(names(r_pws)[r_pws], collapse=", ")))
}
cat("\n")

# === 19F5: GSE197677 PATHWAY LOSO ===
cat("=== 19F5: GSE197677 PATHWAY LOSO ===\n\n")
disc_rank_vec <- setNames(disc_ranks$rank_statistic, disc_ranks$gene)
disc_rank_vec <- disc_rank_vec[!is.na(disc_rank_vec) & is.finite(disc_rank_vec)]

all_samples <- meta_disc$sample
loso_results_disc <- list()

for (s in seq_along(all_samples)) {
  leave_out <- all_samples[s]
  cat(sprintf("  Run %d/%d: leaving out %s\n", s, length(all_samples), leave_out))
  keep <- meta_disc$sample != leave_out
  counts_sub <- raw_counts_disc[, keep]; meta_sub <- meta_disc[keep,]
  treatment <- factor(meta_sub$treatment_standardized, levels=c("No_neoadjuvant_chemotherapy","Neoadjuvant_chemotherapy"))
  design <- model.matrix(~treatment); rownames(design) <- meta_sub$sample
  dge <- DGEList(counts=counts_sub, samples=meta_sub)
  keep_gene <- filterByExpr(dge, design=design)
  dge_f <- dge[keep_gene,,keep.lib.sizes=FALSE]
  dge_f <- calcNormFactors(dge_f, method="TMM")
  dge_f <- estimateDisp(dge_f, design, robust=TRUE)
  fit <- glmQLFit(dge_f, design, robust=TRUE)
  qlf <- glmQLFTest(fit, coef=2)
  res <- topTags(qlf, n=Inf, sort.by="PValue")$table
  loso_rank <- sign(res$logFC)*sqrt(res$F); names(loso_rank) <- rownames(res)
  loso_rank <- loso_rank[!is.na(loso_rank) & is.finite(loso_rank)]
  loso_results_disc[[leave_out]] <- fgseaMultilevel(pathways=hallmark_pathways, stats=loso_rank, minSize=15, maxSize=500, eps=0)
}

disc_loso_summary <- data.frame()
for (pw in supported_pathways) {
  full_nes <- fgsea_disc$NES[fgsea_disc$pathway==pw]
  nes_vec <- numeric(); ss <- 0; f5 <- 0; p5 <- 0; nv <- 0
  for (sn in names(loso_results_disc)) {
    lr <- loso_results_disc[[sn]][loso_results_disc[[sn]]$pathway==pw,]
    if (nrow(lr)>0 && !is.na(lr$NES)) { nv<-nv+1; nes_vec<-c(nes_vec,lr$NES)
      if (sign(lr$NES)==sign(full_nes)) ss<-ss+1
      if (!is.na(lr$padj) && lr$padj<0.05) f5<-f5+1
      if (!is.na(lr$pval) && lr$pval<0.05) p5<-p5+1 }
  }
  sf <- ss/nv; st <- if(sf>=0.90) "DIRECTION_ROBUST" else if(sf>=0.80) "DIRECTION_MODERATE" else "DIRECTION_UNSTABLE"
  disc_loso_summary <- rbind(disc_loso_summary, data.frame(pathway=pw, full_NES=full_nes, valid_runs=nv,
    same_sign=ss, same_sign_frac=sf, FDR05_runs=f5, P05_runs=p5,
    median_NES=median(nes_vec), min_NES=min(nes_vec), max_NES=max(nes_vec), status=st, stringsAsFactors=FALSE))
}

for (i in 1:nrow(disc_loso_summary)) cat(sprintf("  %s: %d/%d (%.0f%%) %s\n", disc_loso_summary$pathway[i], disc_loso_summary$same_sign[i], disc_loso_summary$valid_runs[i], disc_loso_summary$same_sign_frac[i]*100, disc_loso_summary$status[i]))
cat("\n")

# === 19F6: GSE221561 PATHWAY LOSO ===
cat("=== 19F6: GSE221561 PATHWAY LOSO ===\n\n")
treated_samples <- meta_repl$sample[meta_repl$treatment_standardized=="Neoadjuvant_treated"]
loso_results_repl <- list()

for (s in seq_along(treated_samples)) {
  leave_out <- treated_samples[s]
  cat(sprintf("  Run %d/%d: leaving out %s\n", s, length(treated_samples), leave_out))
  keep <- meta_repl$sample != leave_out
  counts_sub <- raw_counts_repl[, keep]; meta_sub <- meta_repl[keep,]
  treatment <- factor(meta_sub$treatment_standardized, levels=c("Surgery_alone","Neoadjuvant_treated"))
  design <- model.matrix(~0+treatment); colnames(design) <- c("Surgery_alone","Neoadjuvant_treated"); rownames(design) <- meta_sub$sample
  contrast_vec <- limma::makeContrasts(Neoadjuvant_treated-Surgery_alone, levels=design)
  dge <- DGEList(counts=counts_sub, samples=meta_sub)
  keep_gene <- filterByExpr(dge, design=design)
  dge_f <- dge[keep_gene,,keep.lib.sizes=FALSE]
  dge_f <- calcNormFactors(dge_f, method="TMM")
  dge_f <- estimateDisp(dge_f, design, robust=TRUE)
  fit <- glmQLFit(dge_f, design, robust=TRUE)
  qlf <- glmQLFTest(fit, contrast=contrast_vec)
  res <- topTags(qlf, n=Inf, sort.by="PValue")$table
  loso_rank <- sign(res$logFC)*sqrt(res$F); names(loso_rank) <- rownames(res)
  loso_rank <- loso_rank[!is.na(loso_rank) & is.finite(loso_rank)]
  loso_results_repl[[leave_out]] <- fgseaMultilevel(pathways=hallmark_pathways, stats=loso_rank, minSize=15, maxSize=500, eps=0)
}

repl_loso_summary <- data.frame()
for (pw in supported_pathways) {
  full_nes <- fgsea_repl$NES[fgsea_repl$pathway==pw]
  nes_vec <- numeric(); ss <- 0; f5 <- 0; p5 <- 0; nv <- 0
  for (sn in names(loso_results_repl)) {
    lr <- loso_results_repl[[sn]][loso_results_repl[[sn]]$pathway==pw,]
    if (nrow(lr)>0 && !is.na(lr$NES)) { nv<-nv+1; nes_vec<-c(nes_vec,lr$NES)
      if (sign(lr$NES)==sign(full_nes)) ss<-ss+1
      if (!is.na(lr$padj) && lr$padj<0.05) f5<-f5+1
      if (!is.na(lr$pval) && lr$pval<0.05) p5<-p5+1 }
  }
  sf <- ss/nv; st <- if(sf>=0.90) "DIRECTION_ROBUST" else if(sf>=0.80) "DIRECTION_MODERATE" else "DIRECTION_UNSTABLE"
  repl_loso_summary <- rbind(repl_loso_summary, data.frame(pathway=pw, full_NES=full_nes, valid_runs=nv,
    same_sign=ss, same_sign_frac=sf, FDR05_runs=f5, P05_runs=p5,
    median_NES=median(nes_vec), min_NES=min(nes_vec), max_NES=max(nes_vec), status=st, stringsAsFactors=FALSE))
}

for (i in 1:nrow(repl_loso_summary)) cat(sprintf("  %s: %d/%d (%.0f%%) %s\n", repl_loso_summary$pathway[i], repl_loso_summary$same_sign[i], repl_loso_summary$valid_runs[i], repl_loso_summary$same_sign_frac[i]*100, repl_loso_summary$status[i]))
cat("\n")

# === 19F7: CROSS-DATASET ROBUSTNESS ===
cat("=== 19F7: CROSS-DATASET ROBUSTNESS CLASSIFICATION ===\n\n")
robustness_final <- data.frame()
for (pw in supported_pathways) {
  dr <- disc_loso_summary[disc_loso_summary$pathway==pw,]; rr <- repl_loso_summary[repl_loso_summary$pathway==pw,]
  lr <- le_repro[le_repro$pathway==pw,]
  df <- fgsea_disc$padj[fgsea_disc$pathway==pw]; rf <- fgsea_repl$padj[fgsea_repl$pathway==pw]
  dn <- fgsea_disc$NES[fgsea_disc$pathway==pw]; rn <- fgsea_repl$NES[fgsea_repl$pathway==pw]
  sd <- sign(dn)==sign(rn)
  fc <- if(df<0.05 && rf<0.05 && sd && dr$same_sign_frac>=0.90 && rr$same_sign_frac>=0.90) "ROBUST_CROSS_DATASET_PATHWAY" else "SUPPORTED_BUT_SAMPLE_SENSITIVE"
  robustness_final <- rbind(robustness_final, data.frame(pathway=pw, disc_NES=dn, disc_FDR=df, repl_NES=rn, repl_FDR=rf,
    same_dir=sd, disc_LOSO_frac=dr$same_sign_frac, disc_LOSO_status=dr$status,
    repl_LOSO_frac=rr$same_sign_frac, repl_LOSO_status=rr$status, LE_Jaccard=lr$jaccard,
    disc_LE_n=lr$disc_LE_n, repl_LE_n=lr$repl_LE_n, LE_overlap=lr$overlap_n, final_class=fc, stringsAsFactors=FALSE))
}

n_robust <- sum(robustness_final$final_class=="ROBUST_CROSS_DATASET_PATHWAY")
n_sens <- sum(robustness_final$final_class=="SUPPORTED_BUT_SAMPLE_SENSITIVE")
cat(sprintf("ROBUST: %d/7\nSUPPORTED_BUT_SENSITIVE: %d/7\n\n", n_robust, n_sens))
for (i in 1:nrow(robustness_final)) cat(sprintf("  %s: %s\n", robustness_final$pathway[i], robustness_final$final_class[i]))

# === 19F8: GENE-DIRECTION AUDIT ===
cat("\n=== 19F8: GENE-DIRECTION SUPPORT ===\n\n")
gene_dir_audit <- data.frame()
for (pw in supported_pathways) {
  pw_genes <- intersect(hallmark_pathways[[pw]], common_genes)
  n_test <- length(pw_genes); n_same <- sum(sign(disc_all[pw_genes,"logFC"])==sign(repl_all[pw_genes,"logFC"]), na.rm=TRUE)
  d_le_c <- intersect(disc_le[[pw]], common_genes)
  d_le_same <- sum(sign(disc_all[d_le_c,"logFC"])==sign(repl_all[d_le_c,"logFC"]), na.rm=TRUE)
  bl <- intersect(disc_le[[pw]], repl_le[[pw]]); bl_c <- intersect(bl, common_genes)
  bl_same <- if(length(bl_c)>0) sum(sign(disc_all[bl_c,"logFC"])==sign(repl_all[bl_c,"logFC"]), na.rm=TRUE) else 0
  gene_dir_audit <- rbind(gene_dir_audit, data.frame(pathway=pw, n_tested=n_test, n_same=n_same, frac_same=n_same/n_test,
    disc_LE_n=length(d_le_c), disc_LE_same=d_le_same, disc_LE_frac=d_le_same/length(d_le_c),
    both_LE_n=length(bl_c), both_LE_same=bl_same, both_LE_frac=if(length(bl_c)>0) bl_same/length(bl_c) else NA, stringsAsFactors=FALSE))
}
for (i in 1:nrow(gene_dir_audit)) cat(sprintf("  %s: %d/%d (%.3f), Disc_LE=%d (%.3f), Both_LE=%d (%.3f)\n",
  gene_dir_audit$pathway[i], gene_dir_audit$n_same[i], gene_dir_audit$n_tested[i], gene_dir_audit$frac_same[i],
  gene_dir_audit$disc_LE_n[i], gene_dir_audit$disc_LE_frac[i], gene_dir_audit$both_LE_n[i], gene_dir_audit$both_LE_frac[i]))

# === 19F9: FIGURES ===
cat("\n=== 19F9: FIGURES ===\n\n")

# NES comparison
nes_comp <- data.frame(pathway=supported_pathways, disc_NES=robustness_final$disc_NES, repl_NES=robustness_final$repl_NES)
p1 <- ggplot(nes_comp, aes(x=disc_NES, y=repl_NES)) + geom_point(size=3) + geom_text(aes(label=pathway), vjust=-0.8, size=3) +
  geom_hline(yintercept=0, linetype="dashed") + geom_vline(xintercept=0, linetype="dashed") +
  labs(title="Seven Supported Hallmarks: NES Comparison", x="GSE197677 NES", y="GSE221561 NES") + theme_minimal()
ggsave(file.path(fig_dir, "STEP19F_SEVEN_PATHWAY_NES_COMPARISON.png"), p1, width=8, height=7, dpi=300)

# LOSO stability
loso_pd <- data.frame()
for (pw in supported_pathways) {
  loso_pd <- rbind(loso_pd, data.frame(pathway=pw, dataset="GSE197677", frac=disc_loso_summary$same_sign_frac[disc_loso_summary$pathway==pw], stringsAsFactors=FALSE))
  loso_pd <- rbind(loso_pd, data.frame(pathway=pw, dataset="GSE221561", frac=repl_loso_summary$same_sign_frac[repl_loso_summary$pathway==pw], stringsAsFactors=FALSE))
}
p2 <- ggplot(loso_pd, aes(x=pathway, y=frac, fill=dataset)) + geom_bar(stat="identity", position="dodge") +
  geom_hline(yintercept=0.9, linetype="dashed", color="red") + labs(title="LOSO Direction Stability", x="", y="Same-sign fraction") +
  scale_fill_manual(values=c("GSE197677"="steelblue","GSE221561"="tomato")) + theme_minimal() + theme(axis.text.x=element_text(angle=45,hjust=1))
ggsave(file.path(fig_dir, "STEP19F_SEVEN_PATHWAY_LOSO_STABILITY.png"), p2, width=10, height=6, dpi=300)

# LE Jaccard
p3 <- ggplot(le_repro, aes(x=reorder(pathway,jaccard), y=jaccard)) + geom_bar(stat="identity", fill="steelblue") +
  geom_hline(yintercept=0.25, linetype="dashed", color="red") + labs(title="Leading-Edge Jaccard", x="", y="Jaccard") +
  theme_minimal() + theme(axis.text.x=element_text(angle=45,hjust=1))
ggsave(file.path(fig_dir, "STEP19F_LEADING_EDGE_CROSS_DATASET_JACCARD.png"), p3, width=8, height=6, dpi=300)

cat("Figures saved.\n\n")

# === 19F10: GUARDRAILS ===
cat("=== 19F10: GUARDRAILS ===\n\n")
guardrails <- "# Step 19F Biological Interpretation Guardrails\n\n1. The 7 pathways represent pathway-specific cross-dataset support.\n2. Global Hallmark NES concordance was null: 23/49 same direction, P=0.775, rho=0.024.\n3. Do not describe the two cohorts as globally transcriptionally concordant.\n4. All 7 supported pathways had positive NES in both comparisons.\n5. Interpret only as higher pathway-ranked expression in neoadjuvant-associated group.\n6. Do NOT automatically write chemotherapy activates pathway X.\n7. Epithelial compartment not restricted to CNV-confirmed malignant cells.\n8. Hallmark gene sets not statistically independent. Pathway redundancy must be reported.\n9. Leading-edge genes are candidate contributors, not validated drivers.\n"
writeLines(guardrails, file.path(pw_dir, "STEP19F_BIOLOGICAL_INTERPRETATION_GUARDRAILS.md"))

# === 19F11: FINAL FREEZE ===
cat("=== 19F11: FINAL FREEZE ===\n\n")
final_freeze <- robustness_final[,c("pathway","disc_NES","disc_FDR","repl_NES","repl_FDR",
  "disc_LOSO_frac","repl_LOSO_frac","disc_LOSO_status","repl_LOSO_status","LE_Jaccard","final_class")]
names(final_freeze) <- c("pathway","discovery_NES","discovery_FDR","replication_NES","replication_FDR",
  "discovery_LOSO_direction_fraction","replication_LOSO_direction_fraction",
  "discovery_LOSO_status","replication_LOSO_status","leading_edge_Jaccard","final_robustness_class")
final_freeze$gene_direction_fraction <- gene_dir_audit$frac_same[match(final_freeze$pathway, gene_dir_audit$pathway)]
utils::write.csv(final_freeze, file.path(pw_dir, "STEP19F_SEVEN_HALLMARKS_FINAL_FREEZE.csv"), row.names=FALSE)
cat("Saved: STEP19F_SEVEN_HALLMARKS_FINAL_FREEZE.csv\n\n")

# Overall status
if (n_robust >= 4) {
  overall <- "PATHWAY_SPECIFIC_ROBUSTNESS_CONFIRMED"
} else if (n_robust >= 1) {
  overall <- "PATHWAY_SUPPORT_PRESENT_BUT_SAMPLE_SENSITIVE"
} else {
  overall <- "PATHWAY_SUPPORT_REDUCED_AFTER_LOSO"
}

cat("============================================\n")
cat("STEP 19F COMPLETE\n")
cat("SUPPORTED HALLMARK ROBUSTNESS AUDIT\n")
cat("============================================\n\n")
cat(sprintf("Frozen supported pathways: %d\n\n", length(supported_pathways)))
cat("--- Significance-set overlap ---\n")
cat(sprintf("Discovery FDR05: %d\n", length(disc_fdr05)))
cat(sprintf("Replication FDR05: %d\n", length(repl_fdr05)))
cat(sprintf("Both FDR05: %d\n", length(both_fdr05)))
cat(sprintf("Both same dir: %d\n", length(both_fdr05_same)))
cat(sprintf("Both opposite: %d\n", length(both_fdr05_opp)))
cat(sprintf("Hypergeom P: %.4e\n\n", p_hyper))
cat("--- Pathway robustness ---\n")
for (i in 1:nrow(robustness_final)) {
  cat(sprintf("%s: Disc NES/FDR=%.3f/%.4f, Repl NES/FDR=%.3f/%.4f, Disc_LOSO=%.0f%%, Repl_LOSO=%.0f%%, LE_J=%.3f, %s\n",
    robustness_final$pathway[i], robustness_final$disc_NES[i], robustness_final$disc_FDR[i],
    robustness_final$repl_NES[i], robustness_final$repl_FDR[i],
    robustness_final$disc_LOSO_frac[i]*100, robustness_final$repl_LOSO_frac[i]*100,
    robustness_final$LE_Jaccard[i], robustness_final$final_class[i]))
}
cat(sprintf("\nROBUST: %d/7\nSUPPORTED_SENSITIVE: %d/7\n", n_robust, n_sens))
cat(sprintf("\nOverall: %s\n", overall))
cat("\nGO: NO\nKEGG: NO\nReactome: NO\nSamples removed: NO\nCells removed: NO\n\n")
cat("STOP HERE.\nDO NOT START STEP 19G AUTOMATICALLY.\n")
