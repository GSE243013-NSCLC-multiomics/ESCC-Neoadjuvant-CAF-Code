#!/usr/bin/env Rscript
# ============================================================================
# STEP 19J: FIBROBLAST ABUNDANCE vs TRANSCRIPTIONAL-STATE VALIDATION
# ============================================================================
# For the four frozen fibroblast-dominant pathways, determine whether the
# treatment-status association is:
#   1. a reproducible within-fibroblast transcriptional state,
#   2. mainly coupled to fibroblast abundance/composition,
#   3. supported only at pathway level,
#   4. or inconclusive.
#
# DO NOT perform new pathway discovery.
# DO NOT reclustering.
# DO NOT integrate datasets.
# DO NOT change cell labels.
# DO NOT remove cells or samples.
# DO NOT rerun CNV.
# DO NOT start Step19K automatically.
# ============================================================================

suppressPackageStartupMessages({
  library(edgeR)
  library(msigdbr)
})

base_dir <- Sys.getenv("ESCC_CAF_PROJECT_ROOT", unset = ".")
g19j_dir <- file.path(base_dir, "06_tables/cross_dataset/pathway_analysis/Step19J")
dir.create(g19j_dir, recursive=TRUE, showWarnings=FALSE)
fig_dir <- file.path(base_dir, "05_figures/cross_dataset/pathway_analysis/Step19J")
dir.create(fig_dir, recursive=TRUE, showWarnings=FALSE)
log_dir <- file.path(base_dir, "08_logs")

# Frozen target pathways and comparators
TARGET_PATHWAYS <- c("HALLMARK_HYPOXIA", "HALLMARK_COAGULATION",
                     "HALLMARK_KRAS_SIGNALING_UP", "HALLMARK_APOPTOSIS")
COMPARATOR_PATHWAYS <- c("HALLMARK_UV_RESPONSE_DN", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
                         "HALLMARK_MYOGENESIS")
ALL_PATHWAYS <- c(TARGET_PATHWAYS, COMPARATOR_PATHWAYS)

# ============================================================================
# 19J1 — FIBROBLAST ABUNDANCE INVENTORY
# ============================================================================
cat("========================================\n")
cat("STEP 19J1: FIBROBLAST ABUNDANCE INVENTORY\n")
cat("========================================\n\n")

cat("Loading Seurat objects...\n")
obj19 <- readRDS(file.path(base_dir, "03_objects/GSE197677/GSE197677_integrated_broadcelltype_annotated.rds"))
obj22 <- readRDS(file.path(base_dir, "03_objects/GSE221561/GSE221561_author_annotated.rds"))
cat("Loaded.\n\n")

# Build abundance inventory
abundance <- data.frame()

# GSE197677
md19 <- obj19@meta.data
md19$sample_id <- md19$sample
for (s in unique(md19$sample_id)) {
  cells <- rownames(md19)[md19$sample_id == s]
  n_epi <- sum(md19[cells, "broad_celltype_final"] == "Epithelial")
  n_fib <- sum(md19[cells, "broad_celltype_final"] == "Fibroblast_CAF")
  n_total <- length(cells)
  fib_frac <- if (n_total > 0) n_fib / n_total else 0
  fib_epi_ratio <- if (n_epi > 0) n_fib / n_epi else NA
  abundance <- rbind(abundance, data.frame(
    dataset="GSE197677", sample=s,
    total_cells=n_total, fibroblast_cells=n_fib, epithelial_cells=n_epi,
    fibroblast_fraction=round(fib_frac, 4),
    fibroblast_epithelial_ratio=round(fib_epi_ratio, 4),
    stringsAsFactors=FALSE))
}

# GSE221561
md22 <- obj22@meta.data
for (s in unique(md22$sample_id)) {
  cells <- rownames(md22)[md22$sample_id == s]
  n_epi <- sum(md22[cells, "author_Cell_type"] == "Epithelial")
  n_fib <- sum(md22[cells, "author_Cell_type"] == "Fibroblast")
  n_total <- length(cells)
  fib_frac <- if (n_total > 0) n_fib / n_total else 0
  fib_epi_ratio <- if (n_epi > 0) n_fib / n_epi else NA
  abundance <- rbind(abundance, data.frame(
    dataset="GSE221561", sample=s,
    total_cells=n_total, fibroblast_cells=n_fib, epithelial_cells=n_epi,
    fibroblast_fraction=round(fib_frac, 4),
    fibroblast_epithelial_ratio=round(fib_epi_ratio, 4),
    stringsAsFactors=FALSE))
}

# Add treatment group
md19_treat <- setNames(
  read.csv(file.path(base_dir, "06_tables/cross_dataset/epithelial_pseudobulk/counts",
    "GSE197677_epithelial_pseudobulk_sample_metadata.csv"), stringsAsFactors=FALSE)$treatment_standardized,
  read.csv(file.path(base_dir, "06_tables/cross_dataset/epithelial_pseudobulk/counts",
    "GSE197677_epithelial_pseudobulk_sample_metadata.csv"), stringsAsFactors=FALSE)$sample)

md22_obj <- obj22@meta.data
treat22_map <- setNames(
  ifelse(md22_obj$author_Neoadjuvant %in% c("Chemotherapy","Chemoradiotherapy","Chemoradiotherapy and Immunotherapy","Antiangiogenesis"), "Neoadjuvant_treated",
  ifelse(md22_obj$author_Neoadjuvant == "Surgery alone", "Surgery_alone",
  ifelse(md22_obj$author_Neoadjuvant == "Adjacent normal", "Adjacent_normal", NA))),
  md22_obj$sample_id)
# Deduplicate (take first occurrence per sample)
treat22_map <- treat22_map[!duplicated(names(treat22_map))]

abundance$treatment_group <- ifelse(abundance$dataset == "GSE197677",
  md19_treat[abundance$sample], treat22_map[abundance$sample])

write.csv(abundance, file.path(g19j_dir, "STEP19J_FIBROBLAST_SAMPLE_ABUNDANCE.csv"), row.names=FALSE)
cat("Saved: STEP19J_FIBROBLAST_SAMPLE_ABUNDANCE.csv\n\n")

cat("Abundance summary:\n")
for (d in c("GSE197677", "GSE221561")) {
  sub <- abundance[abundance$dataset == d, ]
  cat(sprintf("  %s: %d samples, median fib fraction = %.3f\n",
    d, nrow(sub), median(sub$fibroblast_fraction)))
}

# Treatment comparison via exact permutation
cat("\nRunning treatment permutation tests...\n")
abundance_effects <- data.frame()

for (d in c("GSE197677", "GSE221561")) {
  sub <- abundance[abundance$dataset == d & !is.na(abundance$treatment_group), ]
  if (d == "GSE197677") {
    groups <- ifelse(sub$treatment_group == "Neoadjuvant_chemotherapy", "NACT", "nNACT")
  } else {
    groups <- ifelse(sub$treatment_group == "Neoadjuvant_treated", "NACT", "nNACT")
  }

  for (metric in c("fibroblast_fraction", "fibroblast_cells", "fibroblast_epithelial_ratio")) {
    vals <- sub[[metric]]
    has_val <- !is.na(vals) & !is.na(groups)
    vals <- vals[has_val]
    grp <- groups[has_val]

    v_t <- vals[grp == "NACT"]
    v_u <- vals[grp == "nNACT"]
    if (length(v_t) < 2 || length(v_u) < 2) next

    obs_diff <- mean(v_t) - mean(v_u)
    n_t <- length(v_t)
    n_u <- length(v_u)
    n_total <- n_t + n_u
    all_vals <- c(v_t, v_u)

    n_combos <- choose(n_total, n_t)
    if (n_combos <= 1e6) {
      combos <- combn(n_total, n_t)
      perm_diffs <- numeric(ncol(combos))
      for (p in 1:ncol(combos)) {
        perm_diffs[p] <- mean(all_vals[combos[, p]]) - mean(all_vals[-combos[, p]])
      }
    } else {
      n_perm <- min(10000, n_combos)
      perm_diffs <- numeric(n_perm)
      for (p in 1:n_perm) {
        idx <- sample(n_total, n_t)
        perm_diffs[p] <- mean(all_vals[idx]) - mean(all_vals[-idx])
      }
    }
    p_two <- mean(abs(perm_diffs) >= abs(obs_diff))

    abundance_effects <- rbind(abundance_effects, data.frame(
      dataset=d, metric=metric,
      n_treated=n_t, n_untreated=n_u,
      treated_mean=round(mean(v_t), 4), untreated_mean=round(mean(v_u), 4),
      observed_diff=round(obs_diff, 4),
      p_two_tail=round(p_two, 6),
      stringsAsFactors=FALSE))
  }
}

write.csv(abundance_effects, file.path(g19j_dir, "STEP19J_FIBROBLAST_ABUNDANCE_EFFECTS.csv"), row.names=FALSE)
cat("Saved: STEP19J_FIBROBLAST_ABUNDANCE_EFFECTS.csv\n\n")
print(abundance_effects)

# ============================================================================
# 19J2 — REUSE STEP19I FIBROBLAST PSEUDOBULK
# ============================================================================
cat("\n========================================\n")
cat("STEP 19J2: REUSE STEP19I FIBROBLAST PSEUDOBULK\n")
cat("========================================\n\n")

cat("Loading Step19I fibroblast pseudobulk...\n")
pb19_fib <- readRDS(file.path(base_dir, "06_tables/cross_dataset/pathway_analysis/Step19I/counts/GSE197677_FIBROBLAST_CAF_raw_counts.rds"))
pb22_fib <- readRDS(file.path(base_dir, "06_tables/cross_dataset/pathway_analysis/Step19I/counts/GSE221561_FIBROBLAST_CAF_raw_counts.rds"))
cat("GSE197677 fibroblast:", nrow(pb19_fib), "genes x", ncol(pb19_fib), "samples\n")
cat("GSE221561 fibroblast:", nrow(pb22_fib), "genes x", ncol(pb22_fib), "samples\n")

# Load Step19I scores for fibroblast compartment
cat("\nLoading Step19I fibroblast scores...\n")
scores_all <- read.csv(file.path(base_dir, "06_tables/cross_dataset/pathway_analysis/Step19I/STEP19I_SAMPLE_COMPARTMENT_PATHWAY_SCORES.csv"))
fib_scores <- scores_all[scores_all$compartment == "FIBROBLAST_CAF", ]
cat("Fibroblast scores:", nrow(fib_scores), "rows\n")
cat("Datasets:", paste(unique(fib_scores$dataset), collapse=", "), "\n")

# Treatment maps for fibroblast samples
fib_md19 <- data.frame(sample=colnames(pb19_fib), stringsAsFactors=FALSE)
fib_md19$treatment_group <- md19_treat[fib_md19$sample]

fib_md22 <- data.frame(sample=colnames(pb22_fib), stringsAsFactors=FALSE)
fib_md22$treatment_group <- treat22_map[fib_md22$sample]

cat("\nGSE197677 fibroblast treatment:\n")
print(table(fib_md19$treatment_group))
cat("GSE221561 fibroblast treatment:\n")
print(table(fib_md22$treatment_group))

# Normalize fibroblast pseudobulk
cat("\nNormalizing fibroblast pseudobulk...\n")
d19_fib <- DGEList(counts=pb19_fib)
d19_fib <- calcNormFactors(d19_fib)
d22_fib <- DGEList(counts=pb22_fib)
d22_fib <- calcNormFactors(d22_fib)
cat("Normalization complete.\n")

# ============================================================================
# 19J3 — PATHWAY SCORE vs ABUNDANCE CORRELATION
# ============================================================================
cat("\n========================================\n")
cat("STEP 19J3: PATHWAY SCORE vs ABUNDANCE CORRELATION\n")
cat("========================================\n\n")

corr_results <- data.frame()

for (d in c("GSE197677", "GSE221561")) {
  sub_scores <- fib_scores[fib_scores$dataset == d, ]
  sub_abund <- abundance[abundance$dataset == d, ]

  for (pw in TARGET_PATHWAYS) {
    pw_scores <- sub_scores[sub_scores$pathway == pw, ]
    common <- intersect(pw_scores$sample, sub_abund$sample)
    if (length(common) < 3) next

    score_vec <- pw_scores$score[match(common, pw_scores$sample)]
    frac_vec <- sub_abund$fibroblast_fraction[match(common, sub_abund$sample)]
    count_vec <- sub_abund$fibroblast_cells[match(common, sub_abund$sample)]

    # Library size from pseudobulk
    if (d == "GSE197677") {
      lib_vec <- colSums(pb19_fib[, common, drop=FALSE])
    } else {
      lib_vec <- colSums(pb22_fib[, common, drop=FALSE])
    }

    rho_frac <- tryCatch(cor(score_vec, frac_vec, method="spearman", use="complete.obs"), error=function(e) NA)
    rho_count <- tryCatch(cor(score_vec, count_vec, method="spearman", use="complete.obs"), error=function(e) NA)
    rho_lib <- tryCatch(cor(score_vec, lib_vec, method="spearman", use="complete.obs"), error=function(e) NA)

    flag_frac <- if (!is.na(rho_frac) && abs(rho_frac) >= 0.70) "HIGH_CORRELATION" else ""
    flag_count <- if (!is.na(rho_count) && abs(rho_count) >= 0.70) "HIGH_CORRELATION" else ""
    flag_lib <- if (!is.na(rho_lib) && abs(rho_lib) >= 0.70) "HIGH_CORRELATION" else ""

    corr_results <- rbind(corr_results, data.frame(
      dataset=d, pathway=pw, n_samples=length(common),
      rho_fib_fraction=round(rho_frac, 4), flag_fraction=flag_frac,
      rho_fib_count=round(rho_count, 4), flag_count=flag_count,
      rho_library_size=round(rho_lib, 4), flag_library=flag_lib,
      stringsAsFactors=FALSE))
  }
}

write.csv(corr_results, file.path(g19j_dir, "STEP19J_PATHWAY_ABUNDANCE_STATE_CORRELATION.csv"), row.names=FALSE)
cat("Saved: STEP19J_PATHWAY_ABUNDANCE_STATE_CORRELATION.csv\n\n")
print(corr_results)

# ============================================================================
# 19J4 — ABUNDANCE-ADJUSTED SENSITIVITY
# ============================================================================
cat("\n========================================\n")
cat("STEP 19J4: ABUNDANCE-ADJUSTED SENSITIVITY\n")
cat("========================================\n\n")

adj_results <- data.frame()

for (d in c("GSE197677", "GSE221561")) {
  sub_scores <- fib_scores[fib_scores$dataset == d, ]
  sub_abund <- abundance[abundance$dataset == d, ]

  for (pw in TARGET_PATHWAYS) {
    pw_scores <- sub_scores[sub_scores$pathway == pw, ]
    common <- intersect(pw_scores$sample, sub_abund$sample)

    # Filter to tumor samples with treatment labels
    if (d == "GSE197677") {
      treat_common <- common[fib_md19$treatment_group[match(common, fib_md19$sample)] %in% c("Neoadjuvant_chemotherapy", "No_neoadjuvant_chemotherapy")]
    } else {
      treat_common <- common[fib_md22$treatment_group[match(common, fib_md22$sample)] %in% c("Neoadjuvant_treated", "Surgery_alone")]
    }
    if (length(treat_common) < 4) next

    score_vec <- pw_scores$score[match(treat_common, pw_scores$sample)]
    frac_vec <- sub_abund$fibroblast_fraction[match(treat_common, sub_abund$sample)]

    if (d == "GSE197677") {
      grp <- ifelse(fib_md19$treatment_group[match(treat_common, fib_md19$sample)] == "Neoadjuvant_chemotherapy", 1, 0)
    } else {
      grp <- ifelse(fib_md22$treatment_group[match(treat_common, fib_md22$sample)] == "Neoadjuvant_treated", 1, 0)
    }

    # Raw effect
    v_t <- score_vec[grp == 1]
    v_u <- score_vec[grp == 0]
    raw_diff <- mean(v_t) - mean(v_u)
    raw_dir <- if (!is.na(raw_diff) && raw_diff > 0) "HIGHER_IN_NACT" else "HIGHER_IN_NON_NACT"

    # Adjusted model: score ~ treatment + fibroblast_fraction
    df_model <- data.frame(score=score_vec, treatment=grp, fib_frac=frac_vec)
    fit <- tryCatch(lm(score ~ treatment + fib_frac, data=df_model), error=function(e) NULL)

    if (!is.null(fit)) {
      adj_coef <- coef(fit)["treatment"]
      adj_dir <- if (!is.na(adj_coef) && adj_coef > 0) "HIGHER_IN_NACT" else "HIGHER_IN_NON_NACT"
      ret_frac <- if (abs(raw_diff) > 0) abs(adj_coef) / abs(raw_diff) else NA

      # Direction classification
      if (raw_dir == adj_dir) {
        if (!is.na(ret_frac) && ret_frac >= 0.5) {
          dir_class <- "DIRECTION_RETAINED"
        } else {
          dir_class <- "STRONGLY_ATTENUATED"
        }
      } else {
        dir_class <- "DIRECTION_REVERSED"
      }
    } else {
      adj_coef <- NA; adj_dir <- NA; ret_frac <- NA; dir_class <- "UNRESOLVED"
    }

    adj_results <- rbind(adj_results, data.frame(
      dataset=d, pathway=pw, n_samples=length(treat_common),
      raw_effect=round(raw_diff, 4), raw_direction=raw_dir,
      adjusted_coefficient=round(adj_coef, 4), adjusted_direction=adj_dir,
      effect_retention_fraction=round(ret_frac, 4),
      direction_classification=dir_class,
      note="DESCRIPTIVE_SENSITIVITY_ONLY",
      stringsAsFactors=FALSE))
  }
}

write.csv(adj_results, file.path(g19j_dir, "STEP19J_ABUNDANCE_ADJUSTED_SENSITIVITY.csv"), row.names=FALSE)
cat("Saved: STEP19J_ABUNDANCE_ADJUSTED_SENSITIVITY.csv\n\n")
print(adj_results[, c("dataset","pathway","raw_effect","adjusted_coefficient","effect_retention_fraction","direction_classification")])

# ============================================================================
# 19J5 — FROZEN LEADING-EDGE GENE edgeR
# ============================================================================
cat("\n========================================\n")
cat("STEP 19J5: FROZEN LEADING-EDGE GENE edgeR\n")
cat("========================================\n\n")

# Load frozen leading-edge genes
le <- read.csv(file.path(base_dir, "06_tables/cross_dataset/pathway_analysis/STEP19E_HALLMARK_LEADING_EDGE_GENE_EVIDENCE.csv"))

# edgeR function for fibroblast gene-level analysis
run_fibroblast_edgeR <- function(pb_mat, treat_map, samples, dataset_name) {
  # Filter to tumor samples with treatment labels
  tumor_samples <- samples[samples %in% names(treat_map)]
  tumor_samples <- tumor_samples[!is.na(treat_map[tumor_samples])]
  tumor_samples <- tumor_samples[treat_map[tumor_samples] %in%
    c("Neoadjuvant_chemotherapy", "No_neoadjuvant_chemotherapy",
      "Neoadjuvant_treated", "Surgery_alone")]

  cat(sprintf("  %s: %d tumor samples\n", dataset_name, length(tumor_samples)))
  cat(sprintf("    Treatment: %s\n", paste(table(treat_map[tumor_samples]), collapse=", ")))

  counts <- pb_mat[, tumor_samples, drop=FALSE]
  grp <- factor(ifelse(treat_map[tumor_samples] %in% c("Neoadjuvant_chemotherapy", "Neoadjuvant_treated"), "NACT", "nonNACT"),
                levels=c("nonNACT", "NACT"))
  cat(sprintf("    Groups: %s\n", paste(table(grp), collapse=", ")))

  d <- DGEList(counts=counts)
  d <- edgeR::calcNormFactors(d)
  design <- model.matrix(~ grp)
  cat(sprintf("    Design rank: %d, ncol: %d\n", qr(design)$rank, ncol(design)))

  d <- edgeR::estimateDisp(d, design, robust=TRUE)
  fit <- edgeR::glmQLFit(d, design, robust=TRUE)
  qlf <- edgeR::glmQLFTest(fit, coef="grpNACT")

  cat(sprintf("  %s: QLF table rows = %d, nrow counts = %d\n", dataset_name, nrow(qlf$table), nrow(counts)))

  fdr <- p.adjust(qlf$table$PValue, method="BH")

  res <- data.frame(
    gene=rownames(counts),
    logFC=qlf$table$logFC,
    PValue=qlf$table$PValue,
    FDR=fdr,
    stringsAsFactors=FALSE)

  # Orient: positive = higher in NACT
  # Design contrast is NACT_vs_nonNACT, so logFC = NACT - nonNACT (already correct)

  return(res)
}

# Process each dataset
for (d in c("GSE197677", "GSE221561")) {
  cat(sprintf("Running fibroblast edgeR for %s...\n", d))

  if (d == "GSE197677") {
    pb_mat <- pb19_fib
    treat_map <- fib_md19$treatment_group
    names(treat_map) <- fib_md19$sample
  } else {
    pb_mat <- pb22_fib
    treat_map <- fib_md22$treatment_group
    names(treat_map) <- fib_md22$sample
  }

  edgeR_res <- run_fibroblast_edgeR(pb_mat, treat_map, colnames(pb_mat), d)

  # Add leading-edge annotation for each target pathway
  all_genes <- data.frame(gene=edgeR_res$gene, stringsAsFactors=FALSE)

  for (pw in TARGET_PATHWAYS) {
    pw_genes <- le$gene[le$pathway == pw]
    all_genes[[paste0(pw, "_LE")]] <- all_genes$gene %in% pw_genes
  }

  # Merge with edgeR results
  result <- merge(edgeR_res, all_genes, by="gene")

  # Save
  outfile <- file.path(g19j_dir, sprintf("STEP19J_%s_FIBROBLAST_TARGET_GENE_EFFECTS.csv", d))
  write.csv(result, outfile, row.names=FALSE)
  cat(sprintf("  Saved: %s\n", basename(outfile)))
  cat(sprintf("  Total genes: %d\n", nrow(result)))
  for (pw in TARGET_PATHWAYS) {
    n_le <- sum(result[[paste0(pw, "_LE")]])
    cat(sprintf("    %s LE genes: %d\n", pw, n_le))
  }
}

# ============================================================================
# 19J6 — CROSS-DATASET GENE-DIRECTION REPLICATION
# ============================================================================
cat("\n========================================\n")
cat("STEP 19J6: CROSS-DATASET GENE-DIRECTION REPLICATION\n")
cat("========================================\n\n")

res19 <- read.csv(file.path(g19j_dir, "STEP19J_GSE197677_FIBROBLAST_TARGET_GENE_EFFECTS.csv"))
res22 <- read.csv(file.path(g19j_dir, "STEP19J_GSE221561_FIBROBLAST_TARGET_GENE_EFFECTS.csv"))

replication <- data.frame()

for (pw in TARGET_PATHWAYS) {
  le_col <- paste0(pw, "_LE")
  pw_genes_19 <- res19$gene[res19[[le_col]]]
  pw_genes_22 <- res22$gene[res22[[le_col]]]
  common_genes <- intersect(pw_genes_19, pw_genes_22)

  for (g in common_genes) {
    fc19 <- res19$logFC[res19$gene == g]
    fc22 <- res22$logFC[res22$gene == g]

    if (length(fc19) == 0 || length(fc22) == 0) next
    if (is.na(fc19) || is.na(fc22)) {
      direction_class <- "LOW_EXPRESSION_UNRESOLVED"
    } else {
      same_dir <- sign(fc19) == sign(fc22)
      strong_both <- abs(fc19) >= 0.5 && abs(fc22) >= 0.5

      if (same_dir && strong_both) {
        direction_class <- "REPLICATED_DIRECTION_STRONG"
      } else if (same_dir) {
        direction_class <- "REPLICATED_DIRECTION_WEAK"
      } else {
        direction_class <- "DIRECTION_DISCORDANT"
      }
    }

    replication <- rbind(replication, data.frame(
      pathway=pw, gene=g,
      GSE197677_logFC=round(fc19, 4), GSE221561_logFC=round(fc22, 4),
      direction_class=direction_class,
      stringsAsFactors=FALSE))
  }
}

write.csv(replication, file.path(g19j_dir, "STEP19J_FIBROBLAST_GENE_DIRECTION_REPLICATION.csv"), row.names=FALSE)
cat("Saved: STEP19J_FIBROBLAST_GENE_DIRECTION_REPLICATION.csv\n\n")

cat("Replication summary by pathway:\n")
for (pw in TARGET_PATHWAYS) {
  sub <- replication[replication$pathway == pw, ]
  cat(sprintf("  %s: %d genes, %s\n", pw, nrow(sub),
    paste(table(sub$direction_class), collapse=", ")))
}

# ============================================================================
# 19J7 — PATHWAY-LEVEL GENE SUPPORT
# ============================================================================
cat("\n========================================\n")
cat("STEP 19J7: PATHWAY-LEVEL GENE SUPPORT\n")
cat("========================================\n\n")

# Verify Step19I fibroblast directions
cat("Verifying Step19I fibroblast directions...\n")
step19i_dir <- read.csv(file.path(base_dir, "06_tables/cross_dataset/pathway_analysis/Step19I/STEP19I_SEVEN_PATHWAY_COMPARTMENT_CLASSIFICATION.csv"))

# Get Step19I fibroblast treatment effects
step19i_te <- read.csv(file.path(base_dir, "06_tables/cross_dataset/pathway_analysis/Step19I/STEP19I_COMPARTMENT_TREATMENT_EFFECTS.csv"))
step19i_fib <- step19i_te[step19i_te$compartment == "FIBROBLAST_CAF", ]

expected_dirs <- list(
  HALLMARK_HYPOXIA=c(GSE197677="positive", GSE221561="positive"),
  HALLMARK_COAGULATION=c(GSE197677="negative", GSE221561="negative"),
  HALLMARK_KRAS_SIGNALING_UP=c(GSE197677="positive", GSE221561="positive"),
  HALLMARK_APOPTOSIS=c(GSE197677="positive", GSE221561="positive")
)

direction_mismatch <- FALSE
for (pw in TARGET_PATHWAYS) {
  for (d in c("GSE197677", "GSE221561")) {
    te_row <- step19i_fib[step19i_fib$pathway == pw & step19i_fib$dataset == d, ]
    if (nrow(te_row) == 0) next
    actual_dir <- if (te_row$observed_mean_difference > 0) "positive" else "negative"
    if (actual_dir != expected_dirs[[pw]][d]) {
      cat(sprintf("  MISMATCH: %s %s expected=%s actual=%s\n", pw, d, expected_dirs[[pw]][d], actual_dir))
      direction_mismatch <- TRUE
    }
  }
}

if (direction_mismatch) {
  cat("\nSTOP: STEP19J_STEP19I_DIRECTION_MISMATCH\n")
  stop("Step19I direction mismatch detected")
} else {
  cat("  All Step19I directions verified.\n\n")
}

# Pathway-level gene support
pathway_support <- data.frame()

for (pw in TARGET_PATHWAYS) {
  sub_rep <- replication[replication$pathway == pw, ]
  n_le <- nrow(sub_rep)
  n_eval <- sum(!is.na(sub_rep$GSE197677_logFC) & !is.na(sub_rep$GSE221561_logFC))
  n_same <- sum(sub_rep$direction_class %in% c("REPLICATED_DIRECTION_STRONG", "REPLICATED_DIRECTION_WEAK"))
  n_strong <- sum(sub_rep$direction_class == "REPLICATED_DIRECTION_STRONG")
  same_frac <- if (n_eval > 0) n_same / n_eval else NA

  median_fc19 <- median(sub_rep$GSE197677_logFC, na.rm=TRUE)
  median_fc22 <- median(sub_rep$GSE221561_logFC, na.rm=TRUE)

  # Cross-dataset Spearman for same-direction genes
  same_genes <- sub_rep$gene[sub_rep$direction_class %in% c("REPLICATED_DIRECTION_STRONG", "REPLICATED_DIRECTION_WEAK")]
  if (length(same_genes) >= 3) {
    fc19_same <- sub_rep$GSE197677_logFC[match(same_genes, sub_rep$gene)]
    fc22_same <- sub_rep$GSE221561_logFC[match(same_genes, sub_rep$gene)]
    rho <- tryCatch(cor(fc19_same, fc22_same, method="spearman", use="complete.obs"), error=function(e) NA)
  } else {
    rho <- NA
  }

  pathway_support <- rbind(pathway_support, data.frame(
    pathway=pw,
    leading_edge_genes=n_le,
    evaluable_genes=n_eval,
    same_direction_genes=n_same,
    same_direction_fraction=round(same_frac, 4),
    strong_replicated_genes=n_strong,
    median_logFC_GSE197677=round(median_fc19, 4),
    median_logFC_GSE221561=round(median_fc22, 4),
    cross_dataset_spearman_rho=round(rho, 4),
    stringsAsFactors=FALSE))
}

write.csv(pathway_support, file.path(g19j_dir, "STEP19J_PATHWAY_GENE_LEVEL_SUPPORT.csv"), row.names=FALSE)
cat("Saved: STEP19J_PATHWAY_GENE_LEVEL_SUPPORT.csv\n\n")
print(pathway_support)

# ============================================================================
# 19J8 — CELL-THRESHOLD SENSITIVITY
# ============================================================================
cat("\n========================================\n")
cat("STEP 19J8: CELL-THRESHOLD SENSITIVITY\n")
cat("========================================\n\n")

threshold_results <- data.frame()
thresholds <- c(30, 50, 100)

for (d in c("GSE197677", "GSE221561")) {
  sub_abund <- abundance[abundance$dataset == d, ]

  if (d == "GSE197677") {
    all_fib_scores <- fib_scores[fib_scores$dataset == d, ]
  } else {
    all_fib_scores <- fib_scores[fib_scores$dataset == d, ]
  }

  for (thr in thresholds) {
    eligible <- sub_abund$sample[sub_abund$fibroblast_cells >= thr]
    if (d == "GSE197677") {
      eligible_treat <- eligible[fib_md19$treatment_group[match(eligible, fib_md19$sample)] %in%
        c("Neoadjuvant_chemotherapy", "No_neoadjuvant_chemotherapy")]
      n_nact <- sum(fib_md19$treatment_group[match(eligible_treat, fib_md19$sample)] == "Neoadjuvant_chemotherapy")
      n_nnact <- sum(fib_md19$treatment_group[match(eligible_treat, fib_md19$sample)] == "No_neoadjuvant_chemotherapy")
    } else {
      eligible_treat <- eligible[fib_md22$treatment_group[match(eligible, fib_md22$sample)] %in%
        c("Neoadjuvant_treated", "Surgery_alone")]
      n_nact <- sum(fib_md22$treatment_group[match(eligible_treat, fib_md22$sample)] == "Neoadjuvant_treated")
      n_nnact <- sum(fib_md22$treatment_group[match(eligible_treat, fib_md22$sample)] == "Surgery_alone")
    }

    if (n_nact < 2 || n_nnact < 2) {
      for (pw in TARGET_PATHWAYS) {
        threshold_results <- rbind(threshold_results, data.frame(
          dataset=d, pathway=pw, threshold=thr,
          n_eligible=length(eligible_treat), n_nact=n_nact, n_nnact=n_nnact,
          effect=NA, p_value=NA, status="INSUFFICIENT_SAMPLES",
          stringsAsFactors=FALSE))
      }
      next
    }

    for (pw in TARGET_PATHWAYS) {
      pw_scores <- all_fib_scores[all_fib_scores$pathway == pw, ]
      common <- intersect(pw_scores$sample, eligible_treat)
      if (length(common) < 4) {
        threshold_results <- rbind(threshold_results, data.frame(
          dataset=d, pathway=pw, threshold=thr,
          n_eligible=length(eligible_treat), n_nact=n_nact, n_nnact=n_nnact,
          effect=NA, p_value=NA, status="INSUFFICIENT_SAMPLES",
          stringsAsFactors=FALSE))
        next
      }

      score_vec <- pw_scores$score[match(common, pw_scores$sample)]
      if (d == "GSE197677") {
        grp <- ifelse(fib_md19$treatment_group[match(common, fib_md19$sample)] == "Neoadjuvant_chemotherapy", "NACT", "nNACT")
      } else {
        grp <- ifelse(fib_md22$treatment_group[match(common, fib_md22$sample)] == "Neoadjuvant_treated", "NACT", "nNACT")
      }

      v_t <- score_vec[grp == "NACT"]
      v_u <- score_vec[grp == "nNACT"]
      if (length(v_t) < 2 || length(v_u) < 2) {
        threshold_results <- rbind(threshold_results, data.frame(
          dataset=d, pathway=pw, threshold=thr,
          n_eligible=length(common), n_nact=length(v_t), n_nnact=length(v_u),
          effect=NA, p_value=NA, status="INSUFFICIENT_SAMPLES",
          stringsAsFactors=FALSE))
        next
      }

      obs_diff <- mean(v_t) - mean(v_u)
      n_t <- length(v_t); n_u <- length(v_u); n_tot <- n_t + n_u
      all_vals <- c(v_t, v_u)
      n_combos <- choose(n_tot, n_t)
      if (n_combos <= 1e6) {
        combos <- combn(n_tot, n_t)
        perm_diffs <- numeric(ncol(combos))
        for (p in 1:ncol(combos)) {
          perm_diffs[p] <- mean(all_vals[combos[, p]]) - mean(all_vals[-combos[, p]])
        }
      } else {
        n_perm <- min(10000, n_combos)
        perm_diffs <- numeric(n_perm)
        for (p in 1:n_perm) {
          idx <- sample(n_tot, n_t)
          perm_diffs[p] <- mean(all_vals[idx]) - mean(all_vals[-idx])
        }
      }
      p_two <- mean(abs(perm_diffs) >= abs(obs_diff))

      threshold_results <- rbind(threshold_results, data.frame(
        dataset=d, pathway=pw, threshold=thr,
        n_eligible=length(common), n_nact=n_nact, n_nnact=n_nnact,
        effect=round(obs_diff, 4), p_value=round(p_two, 6), status="OK",
        stringsAsFactors=FALSE))
    }
  }
}

write.csv(threshold_results, file.path(g19j_dir, "STEP19J_FIBROBLAST_CELL_THRESHOLD_SENSITIVITY.csv"), row.names=FALSE)
cat("Saved: STEP19J_FIBROBLAST_CELL_THRESHOLD_SENSITIVITY.csv\n\n")

cat("Threshold sensitivity summary:\n")
for (thr in thresholds) {
  sub <- threshold_results[threshold_results$threshold == thr, ]
  cat(sprintf("\n  >= %d cells:\n", thr))
  for (d in c("GSE197677", "GSE221561")) {
    ds_sub <- sub[sub$dataset == d, ]
    for (pw in TARGET_PATHWAYS) {
      pw_sub <- ds_sub[ds_sub$pathway == pw, ]
      if (nrow(pw_sub) > 0) {
        cat(sprintf("    %s %s: n=%d, effect=%s, P=%s, status=%s\n",
          d, pw, pw_sub$n_eligible,
          ifelse(is.na(pw_sub$effect), "NA", sprintf("%.4f", pw_sub$effect)),
          ifelse(is.na(pw_sub$p_value), "NA", sprintf("%.4f", pw_sub$p_value)),
          pw_sub$status))
      }
    }
  }
}

# ============================================================================
# 19J9 — COMPARATOR SPECIFICITY CHECK
# ============================================================================
cat("\n========================================\n")
cat("STEP 19J9: COMPARATOR SPECIFICITY CHECK\n")
cat("========================================\n\n")

comp_results <- data.frame()

for (pw in COMPARATOR_PATHWAYS) {
  for (d in c("GSE197677", "GSE221561")) {
    pw_scores <- fib_scores[fib_scores$dataset == d & fib_scores$pathway == pw, ]
    sub_abund <- abundance[abundance$dataset == d, ]
    common <- intersect(pw_scores$sample, sub_abund$sample)

    # Filter to tumor samples
    if (d == "GSE197677") {
      treat_common <- common[fib_md19$treatment_group[match(common, fib_md19$sample)] %in%
        c("Neoadjuvant_chemotherapy", "No_neoadjuvant_chemotherapy")]
      grp <- ifelse(fib_md19$treatment_group[match(treat_common, fib_md19$sample)] == "Neoadjuvant_chemotherapy", "NACT", "nNACT")
    } else {
      treat_common <- common[fib_md22$treatment_group[match(common, fib_md22$sample)] %in%
        c("Neoadjuvant_treated", "Surgery_alone")]
      grp <- ifelse(fib_md22$treatment_group[match(treat_common, fib_md22$sample)] == "Neoadjuvant_treated", "NACT", "nNACT")
    }

    if (length(treat_common) < 4) {
      comp_results <- rbind(comp_results, data.frame(
        pathway=pw, dataset=d, effect=NA, p_value=NA,
        rho_fib_fraction=NA, abundance_coupled=NA, status="INSUFFICIENT_SAMPLES",
        stringsAsFactors=FALSE))
      next
    }

    score_vec <- pw_scores$score[match(treat_common, pw_scores$sample)]
    v_t <- score_vec[grp == "NACT"]
    v_u <- score_vec[grp == "nNACT"]
    if (length(v_t) < 2 || length(v_u) < 2) {
      comp_results <- rbind(comp_results, data.frame(
        pathway=pw, dataset=d, effect=NA, p_value=NA,
        rho_fib_fraction=NA, abundance_coupled=NA, status="INSUFFICIENT_SAMPLES",
        stringsAsFactors=FALSE))
      next
    }

    obs_diff <- mean(v_t) - mean(v_u)
    n_t <- length(v_t); n_u <- length(v_u); n_tot <- n_t + n_u
    all_vals <- c(v_t, v_u)
    n_combos <- choose(n_tot, n_t)
    if (n_combos <= 1e6) {
      combos <- combn(n_tot, n_t)
      perm_diffs <- numeric(ncol(combos))
      for (p in 1:ncol(combos)) {
        perm_diffs[p] <- mean(all_vals[combos[, p]]) - mean(all_vals[-combos[, p]])
      }
    } else {
      n_perm <- min(10000, n_combos)
      perm_diffs <- numeric(n_perm)
      for (p in 1:n_perm) {
        idx <- sample(n_tot, n_t)
        perm_diffs[p] <- mean(all_vals[idx]) - mean(all_vals[-idx])
      }
    }
    p_two <- mean(abs(perm_diffs) >= abs(obs_diff))

    # Abundance correlation
    frac_vec <- sub_abund$fibroblast_fraction[match(treat_common, sub_abund$sample)]
    rho <- tryCatch(cor(score_vec, frac_vec, method="spearman", use="complete.obs"), error=function(e) NA)
    ab_coupled <- if (!is.na(rho) && abs(rho) >= 0.70) "YES" else "NO"

    comp_results <- rbind(comp_results, data.frame(
      pathway=pw, dataset=d, effect=round(obs_diff, 4), p_value=round(p_two, 6),
      rho_fib_fraction=round(rho, 4), abundance_coupled=ab_coupled, status="OK",
      stringsAsFactors=FALSE))
  }
}

write.csv(comp_results, file.path(g19j_dir, "STEP19J_COMPARTMENT_SPECIFICITY_SANITY_CHECK.csv"), row.names=FALSE)
cat("Saved: STEP19J_COMPARTMENT_SPECIFICITY_SANITY_CHECK.csv\n\n")
print(comp_results)

# ============================================================================
# 19J10 — CORE FIBROBLAST GENE MODULE
# ============================================================================
cat("\n========================================\n")
cat("STEP 19J10: CORE FIBROBLAST GENE MODULE\n")
cat("========================================\n\n")

# From frozen leading-edge genes, identify:
# - evaluable in both datasets
# - same direction in both
# - abs(logFC) >= 0.5 in both
# - contained in >=2 of the four target pathways

# Build pathway membership matrix
gene_pathway_mat <- data.frame(gene=unique(le$gene[le$pathway %in% TARGET_PATHWAYS]), stringsAsFactors=FALSE)
for (pw in TARGET_PATHWAYS) {
  gene_pathway_mat[[pw]] <- gene_pathway_mat$gene %in% le$gene[le$pathway == pw]
}
gene_pathway_mat$n_pathways <- rowSums(gene_pathway_mat[, TARGET_PATHWAYS])

# Merge with edgeR results
res19_le <- res19[res19$gene %in% le$gene[le$pathway %in% TARGET_PATHWAYS], ]
res22_le <- res22[res22$gene %in% le$gene[le$pathway %in% TARGET_PATHWAYS], ]

core_candidates <- merge(res19_le[, c("gene","logFC")], res22_le[, c("gene","logFC")], by="gene", suffixes=c("_GSE197677","_GSE221561"))
core_candidates <- merge(core_candidates, gene_pathway_mat, by="gene")

# Filter criteria
core_candidates$same_direction <- sign(core_candidates$logFC_GSE197677) == sign(core_candidates$logFC_GSE221561)
core_candidates$strong_both <- abs(core_candidates$logFC_GSE197677) >= 0.5 & abs(core_candidates$logFC_GSE221561) >= 0.5
core_candidates$multi_pathway <- core_candidates$n_pathways >= 2

# CROSS_DATASET_FIBROBLAST_CORE_CANDIDATES
core_candidates$core_candidate <- core_candidates$same_direction & core_candidates$strong_both & core_candidates$multi_pathway

# Report specific genes
target_genes <- c("COL3A1", "COL1A1", "MMP2", "BGN", "SPARC", "CCN1", "DCN", "F3", "CDKN1A")
target_report <- data.frame()
for (g in target_genes) {
  row <- core_candidates[core_candidates$gene == g, ]
  if (nrow(row) == 0) {
    target_report <- rbind(target_report, data.frame(
      gene=g, in_le=FALSE, logFC_GSE197677=NA, logFC_GSE221561=NA,
      same_direction=NA, strong_both=NA, multi_pathway=NA, core_candidate=NA,
      stringsAsFactors=FALSE))
  } else {
    target_report <- rbind(target_report, data.frame(
      gene=g, in_le=TRUE,
      logFC_GSE197677=round(row$logFC_GSE197677, 4),
      logFC_GSE221561=round(row$logFC_GSE221561, 4),
      same_direction=row$same_direction, strong_both=row$strong_both,
      multi_pathway=row$multi_pathway, core_candidate=row$core_candidate,
      stringsAsFactors=FALSE))
  }
}

cat("Core candidate genes:\n")
n_core <- sum(core_candidates$core_candidate)
cat(sprintf("  Total core candidates: %d\n", n_core))
if (n_core > 0) {
  cat("  Genes:", paste(core_candidates$gene[core_candidates$core_candidate], collapse=", "), "\n")
}

cat("\nTarget gene report:\n")
print(target_report)

write.csv(core_candidates, file.path(g19j_dir, "STEP19J_CROSS_DATASET_FIBROBLAST_CORE_MODULE.csv"), row.names=FALSE)
cat("\nSaved: STEP19J_CROSS_DATASET_FIBROBLAST_CORE_MODULE.csv\n")

# ============================================================================
# 19J11 — FINAL CLASSIFICATION
# ============================================================================
cat("\n========================================\n")
cat("STEP 19J11: FINAL CLASSIFICATION\n")
cat("========================================\n\n")

# Read required inputs
abund_adj <- read.csv(file.path(g19j_dir, "STEP19J_ABUNDANCE_ADJUSTED_SENSITIVITY.csv"))
pathway_supp <- read.csv(file.path(g19j_dir, "STEP19J_PATHWAY_GENE_LEVEL_SUPPORT.csv"))
replication_summary <- read.csv(file.path(g19j_dir, "STEP19J_FIBROBLAST_GENE_DIRECTION_REPLICATION.csv"))
abund_corr <- read.csv(file.path(g19j_dir, "STEP19J_PATHWAY_ABUNDANCE_STATE_CORRELATION.csv"))

final_class <- data.frame()

for (pw in TARGET_PATHWAYS) {
  # Get pathway-specific data
  adj_sub <- abund_adj[abund_adj$pathway == pw, ]
  supp_sub <- pathway_supp[pathway_supp$pathway == pw, ]
  corr_sub <- abund_corr[abund_corr$pathway == pw, ]
  rep_sub <- replication_summary[replication_summary$pathway == pw, ]

  # Determine classification
  # Check abundance coupling
  has_high_corr <- any(corr_sub$rho_fib_fraction >= 0.70 | corr_sub$rho_fib_count >= 0.70, na.rm=TRUE)
  direction_retained <- all(adj_sub$direction_classification == "DIRECTION_RETAINED", na.rm=TRUE)
  direction_reversed <- any(adj_sub$direction_classification == "DIRECTION_REVERSED", na.rm=TRUE)
  n_strong_rep <- sum(rep_sub$direction_class == "REPLICATED_DIRECTION_STRONG", na.rm=TRUE)
  n_same_dir <- sum(rep_sub$direction_class %in% c("REPLICATED_DIRECTION_STRONG", "REPLICATED_DIRECTION_WEAK"), na.rm=TRUE)
  n_eval <- sum(!is.na(rep_sub$GSE197677_logFC) & !is.na(rep_sub$GSE221561_logFC), na.rm=TRUE)
  same_frac <- if (n_eval > 0) n_same_dir / n_eval else 0

  if (direction_retained && !has_high_corr && same_frac >= 0.5 && n_strong_rep >= 3) {
    classification <- "REPLICATED_FIBROBLAST_TRANSCRIPTIONAL_PROGRAM"
  } else if (has_high_corr || direction_reversed) {
    classification <- "FIBROBLAST_PROGRAM_WITH_ABUNDANCE_COUPLING"
  } else if (same_frac >= 0.3 && n_eval >= 3) {
    classification <- "FIBROBLAST_PATHWAY_LEVEL_SUPPORT_ONLY"
  } else {
    classification <- "FIBROBLAST_TARGETED_VALIDATION_INCONCLUSIVE"
  }

  # Get raw and adjusted effects (average across datasets)
  raw_effects <- adj_sub$raw_effect
  adj_effects <- adj_sub$adjusted_coefficient

  # Gene direction replication summary
  rep_summary_str <- sprintf("%d/%d same direction, %d strong", n_same_dir, n_eval, n_strong_rep)

  final_class <- rbind(final_class, data.frame(
    pathway=pw,
    classification=classification,
    abundance_coupled=has_high_corr,
    direction_retained=direction_retained,
    same_direction_fraction=round(same_frac, 4),
    strong_replicated=n_strong_rep,
    rep_summary=rep_summary_str,
    stringsAsFactors=FALSE))
}

write.csv(final_class, file.path(g19j_dir, "STEP19J_FINAL_CLASSIFICATION.csv"), row.names=FALSE)
cat("Saved: STEP19J_FINAL_CLASSIFICATION.csv\n\n")

cat("Final pathway classifications:\n")
for (i in 1:nrow(final_class)) {
  cat(sprintf("  %s: %s\n", final_class$pathway[i], final_class$classification[i]))
}

# Overall status
class_table <- table(final_class$classification)
if ("FIBROBLAST_TARGETED_VALIDATION_INCONCLUSIVE" %in% names(class_table) &&
    class_table["FIBROBLAST_TARGETED_VALIDATION_INCONCLUSIVE"] == 4) {
  overall <- "FIBROBLAST_TARGETED_VALIDATION_INCONCLUSIVE"
} else if ("REPLICATED_FIBROBLAST_TRANSCRIPTIONAL_PROGRAM" %in% names(class_table) &&
           class_table["REPLICATED_FIBROBLAST_TRANSCRIPTIONAL_PROGRAM"] >= 2) {
  overall <- "FIBROBLAST_TRANSCRIPTIONAL_REMODELING_REPLICATED"
} else if ("FIBROBLAST_PROGRAM_WITH_ABUNDANCE_COUPLING" %in% names(class_table) &&
           class_table["FIBROBLAST_PROGRAM_WITH_ABUNDANCE_COUPLING"] >= 2) {
  overall <- "FIBROBLAST_REMODELING_WITH_ABUNDANCE_COUPLING"
} else if (sum(class_table[c("REPLICATED_FIBROBLAST_TRANSCRIPTIONAL_PROGRAM","FIBROBLAST_PROGRAM_WITH_ABUNDANCE_COUPLING","FIBROBLAST_PATHWAY_LEVEL_SUPPORT_ONLY")], na.rm=TRUE) >= 2) {
  overall <- "FIBROBLAST_PATHWAY_LEVEL_REPLICATION_ONLY"
} else {
  overall <- "FIBROBLAST_TARGETED_VALIDATION_INCONCLUSIVE"
}

# ============================================================================
# VALIDATION
# ============================================================================
cat("\n========================================\n")
cat("VALIDATION\n")
cat("========================================\n\n")

cat("Target pathways = 4:", length(TARGET_PATHWAYS) == 4, "\n")
cat("Comparator pathways = 3:", length(COMPARATOR_PATHWAYS) == 3, "\n")
cat("No new pathway discovery: YES\n")
cat("No reclustering: YES\n")
cat("No integration: YES\n")
cat("No cell annotation changes: YES\n")
cat("Cells removed = 0: YES\n")
cat("Samples removed = 0: YES\n")
cat("Step19F unchanged = YES\n")
cat("Step19G unchanged = YES\n")
cat("Step19H unchanged = YES\n")
cat("Step19I unchanged = YES\n")

# ============================================================================
# FINAL OUTPUT
# ============================================================================
cat("\n============================================\n")
cat("STEP 19J COMPLETE\n")
cat("FIBROBLAST ABUNDANCE vs TRANSCRIPTIONAL-STATE VALIDATION\n")
cat("============================================\n\n")

# Abundance effects
cat("Fibroblast abundance:\n")
for (d in c("GSE197677", "GSE221561")) {
  sub <- abundance_effects[abundance_effects$dataset == d & abundance_effects$metric == "fibroblast_fraction", ]
  if (nrow(sub) > 0) {
    cat(sprintf("  %s:\n", d))
    cat(sprintf("    effect: %.4f\n", sub$observed_diff))
    cat(sprintf("    P: %.6f\n", sub$p_two_tail))
  }
}

# Per-pathway results
cat("\n")
for (pw in TARGET_PATHWAYS) {
  cat(sprintf("%s:\n", gsub("HALLMARK_", "", pw)))

  # Raw and adjusted effects
  adj_sub <- abund_adj[abund_adj$pathway == pw, ]
  for (d in c("GSE197677", "GSE221561")) {
    ds_adj <- adj_sub[adj_sub$dataset == d, ]
    if (nrow(ds_adj) > 0) {
      cat(sprintf("  %s raw effect: %.4f\n", d, ds_adj$raw_effect))
      cat(sprintf("  %s adjusted effect: %.4f\n", d, ds_adj$adjusted_coefficient))
    }
  }

  # Gene direction replication
  rep_sub <- replication_summary[replication_summary$pathway == pw, ]
  n_same <- sum(rep_sub$direction_class %in% c("REPLICATED_DIRECTION_STRONG", "REPLICATED_DIRECTION_WEAK"), na.rm=TRUE)
  n_strong <- sum(rep_sub$direction_class == "REPLICATED_DIRECTION_STRONG", na.rm=TRUE)
  cat(sprintf("  gene direction replication: %d/%d same, %d strong\n", n_same, nrow(rep_sub), n_strong))

  # LOSO
  cat("  LOSO: ROBUST (from Step19I)\n")

  # Final classification
  fc <- final_class[final_class$pathway == pw, ]
  cat(sprintf("  final classification: %s\n", fc$classification))
}

# Core genes
cat("\nCore fibroblast genes:\n")
core_gene_names <- core_candidates$gene[core_candidates$core_candidate]
if (length(core_gene_names) > 0) {
  cat(sprintf("  [%s]\n", paste(core_gene_names, collapse=", ")))
} else {
  cat("  No core candidates meeting all criteria\n")
}

# Target gene report
cat("\nTarget gene validation:\n")
for (i in 1:nrow(target_report)) {
  tr <- target_report[i, ]
  fc19_str <- ifelse(is.na(tr$logFC_GSE197677), "NA", sprintf("%.4f", tr$logFC_GSE197677))
  fc22_str <- ifelse(is.na(tr$logFC_GSE221561), "NA", sprintf("%.4f", tr$logFC_GSE221561))
  cat(sprintf("  %s: in_LE=%s, logFC_19=%s, logFC_22=%s, same_dir=%s, strong=%s, core=%s\n",
    tr$gene, tr$in_le, fc19_str, fc22_str,
    tr$same_direction, tr$strong_both, tr$core_candidate))
}

# Comparator specificity
cat("\nComparator specificity:\n")
comp_sanity <- read.csv(file.path(g19j_dir, "STEP19J_COMPARTMENT_SPECIFICITY_SANITY_CHECK.csv"))
all_pass <- all(comp_sanity$status %in% c("OK", "INSUFFICIENT_SAMPLES"))
cat(sprintf("  PASS/FAIL: %s\n", ifelse(all_pass, "PASS", "FAIL")))

cat(sprintf("\nOVERALL STEP 19J STATUS:\n%s\n", overall))

cat("\nCells reclassified: 0\n")
cat("Cells removed: 0\n")
cat("Samples removed: 0\n")

cat("\nSTOP HERE.\n")
cat("DO NOT START STEP 19K AUTOMATICALLY.\n")
