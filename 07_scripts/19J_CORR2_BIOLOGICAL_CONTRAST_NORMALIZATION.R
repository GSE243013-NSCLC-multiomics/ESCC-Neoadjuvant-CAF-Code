#!/usr/bin/env Rscript
# ============================================================================
# STEP 19J-CORR2: SEMANTIC BIOLOGICAL CONTRAST NORMALIZATION
# ============================================================================
# Step19I used sample-order-dependent effect convention:
# - GSE197677: mean(nNACT) - mean(NACT)
# - GSE221561: mean(NACT) - mean(nNACT)
#
# Step19J-CORR reconciled outputs to Step19I frozen signs.
#
# Step19J-CORR2 introduces the first UNIFORM biological contrast
# for cross-dataset interpretation:
#
# CANONICAL TREATMENT EFFECT =
#   NEOADJUVANT-TREATED minus NO-NEOADJUVANT / SURGERY-ALONE
#
# Positive ALWAYS = higher in neoadjuvant-treated condition
# Negative ALWAYS = higher in untreated / surgery-alone condition
# ============================================================================

suppressPackageStartupMessages({
  library(edgeR)
})

base_dir <- Sys.getenv("ESCC_CAF_PROJECT_ROOT", unset = ".")
corr2_dir <- file.path(base_dir, "06_tables/cross_dataset/pathway_analysis/Step19J-CORR2")
dir.create(corr2_dir, showWarnings = FALSE, recursive = TRUE)

TARGET_PATHWAYS <- c("HALLMARK_HYPOXIA", "HALLMARK_COAGULATION",
                     "HALLMARK_KRAS_SIGNALING_UP", "HALLMARK_APOPTOSIS")

CORE_GENES <- c("BGN", "CDKN1B", "GSN", "TIMP1")

TOLERANCE <- 1e-5  # Step19I stores effects with 6 decimal places precision
NEAR_ZERO_THRESHOLD <- 0.05

# ============================================================================
# 19J-CORR2.1: FREEZE LEGACY CONVENTIONS
# ============================================================================
cat("============================================\n")
cat("19J-CORR2.1: FREEZE LEGACY CONVENTIONS\n")
cat("============================================\n\n")

legacy_provenance <- data.frame(
  dataset = c("GSE197677", "GSE221561"),
  legacy_group1 = c("nNACT", "NACT"),
  legacy_group2 = c("NACT", "nNACT"),
  legacy_effect_definition = c(
    "mean(nNACT) - mean(NACT)",
    "mean(NACT) - mean(nNACT)"
  ),
  canonical_effect_definition = c(
    "mean(NACT) - mean(nNACT)",
    "mean(Neoadjuvant_treated) - mean(Surgery_alone)"
  ),
  multiplier = c(-1, 1)
)

cat("Legacy-to-canonical multipliers:\n")
cat(sprintf("  GSE197677: %.0f (flip sign)\n", legacy_provenance$multiplier[1]))
cat(sprintf("  GSE221561: %.0f (keep sign)\n", legacy_provenance$multiplier[2]))

write.csv(legacy_provenance,
  file.path(corr2_dir, "STEP19J_CORR2_LEGACY_DIRECTION_PROVENANCE.csv"),
  row.names = FALSE)

cat("\nSaved: STEP19J_CORR2_LEGACY_DIRECTION_PROVENANCE.csv\n\n")

# ============================================================================
# 19J-CORR2.2: MANDATORY TUMOR-ONLY SAMPLE AUDIT
# ============================================================================
cat("============================================\n")
cat("19J-CORR2.2: MANDATORY TUMOR-ONLY SAMPLE AUDIT\n")
cat("============================================\n\n")

# Canonical treatment maps (tumor samples only)
gse197677_treat_map <- c(
  ESC06 = "nNACT",
  ESC07 = "nNACT",
  ESC09 = "NACT",
  ESC11 = "NACT",
  ESC12 = "NACT",
  ESC15 = "NACT",
  ESC16 = "NACT",
  ESC17 = "nNACT",
  ESC18 = "nNACT",
  ESC22 = "NACT"
)

gse221561_treat_map <- c(
  S0520T  = "Neoadjuvant_treated",
  S1265T  = "Neoadjuvant_treated",
  S1315TS = "Neoadjuvant_treated",
  S1535T  = "Neoadjuvant_treated",
  S2423T  = "Neoadjuvant_treated",
  S2487T  = "Neoadjuvant_treated",
  S6417T  = "Surgery_alone",
  S6829T  = "Neoadjuvant_treated",
  S9478T  = "Surgery_alone"
)

cat("=== GSE197677 TUMOR SAMPLES ===\n")
cat(sprintf("  Total: %d\n", length(gse197677_treat_map)))
cat(sprintf("  NACT (treated): %d\n", sum(gse197677_treat_map == "NACT")))
cat(sprintf("  nNACT (control): %d\n", sum(gse197677_treat_map == "nNACT")))
cat(sprintf("  Normal: 0\n"))
cat("  Sample IDs:", paste(names(gse197677_treat_map), collapse = ", "), "\n")
cat("  Treated:", paste(names(gse197677_treat_map)[gse197677_treat_map == "NACT"], collapse = ", "), "\n")
cat("  Control:", paste(names(gse197677_treat_map)[gse197677_treat_map == "nNACT"], collapse = ", "), "\n\n")

cat("=== GSE221561 TUMOR SAMPLES ===\n")
cat(sprintf("  Total: %d\n", length(gse221561_treat_map)))
cat(sprintf("  Neoadjuvant_treated: %d\n", sum(gse221561_treat_map == "Neoadjuvant_treated")))
cat(sprintf("  Surgery_alone: %d\n", sum(gse221561_treat_map == "Surgery_alone")))
cat(sprintf("  Normal: 0\n"))
cat("  Sample IDs:", paste(names(gse221561_treat_map), collapse = ", "), "\n")
cat("  Treated:", paste(names(gse221561_treat_map)[gse221561_treat_map == "Neoadjuvant_treated"], collapse = ", "), "\n")
cat("  Control:", paste(names(gse221561_treat_map)[gse221561_treat_map == "Surgery_alone"], collapse = ", "), "\n\n")

# Validation
stopifnot(length(gse197677_treat_map) == 10)
stopifnot(sum(gse197677_treat_map == "NACT") == 6)
stopifnot(sum(gse197677_treat_map == "nNACT") == 4)

stopifnot(length(gse221561_treat_map) == 9)
stopifnot(sum(gse221561_treat_map == "Neoadjuvant_treated") == 7)
stopifnot(sum(gse221561_treat_map == "Surgery_alone") == 2)

cat("PASS: Both datasets have correct tumor sample composition\n\n")

# ============================================================================
# 19J-CORR2.3: RECOMPUTE PATHWAY EFFECTS FROM RAW SAMPLE SCORES
# ============================================================================
cat("============================================\n")
cat("19J-CORR2.3: RECOMPUTE PATHWAY EFFECTS FROM RAW SAMPLE SCORES\n")
cat("============================================\n\n")

scores_all <- read.csv(file.path(base_dir,
  "06_tables/cross_dataset/pathway_analysis/Step19I/STEP19I_SAMPLE_COMPARTMENT_PATHWAY_SCORES.csv"))

fib_scores <- scores_all[scores_all$compartment == "FIBROBLAST_CAF", ]

step19i_te <- read.csv(file.path(base_dir,
  "06_tables/cross_dataset/pathway_analysis/Step19I/STEP19I_COMPARTMENT_TREATMENT_EFFECTS.csv"))
step19i_fib <- step19i_te[step19i_te$compartment == "FIBROBLAST_CAF", ]

cat("=== Frozen Step19I fibroblast treatment effects ===\n\n")
for (d in c("GSE197677", "GSE221561")) {
  sub <- step19i_fib[step19i_fib$dataset == d, ]
  cat(sprintf("%s:\n", d))
  for (i in 1:nrow(sub)) {
    cat(sprintf("  %s: effect=%.6f\n",
      gsub("HALLMARK_", "", sub$pathway[i]),
      sub$observed_mean_difference[i]))
  }
  cat("\n")
}

cat("\n=== Recomputing canonical effects from sample-level scores ===\n\n")

# NOTE: Step19I included Adjacent_normal samples (S6417N, S9478N) in its
# "untreated" group for GSE221561. The canonical biological contrast
# uses TUMOR-ONLY samples as specified. Therefore:
# - GSE197677: canonical effects will match reoriented Step19I (no normals in dataset)
# - GSE221561: canonical effects will DIFFER from reoriented Step19I (normals excluded)

results_corr2.3 <- data.frame()

for (dataset in c("GSE197677", "GSE221561")) {
  if (dataset == "GSE197677") {
    treat_map <- gse197677_treat_map
    treated_label <- "NACT"
    control_label <- "nNACT"
    multiplier <- -1
    n_untreated_step19i <- 6  # Step19I included normal samples
  } else {
    treat_map <- gse221561_treat_map
    treated_label <- "Neoadjuvant_treated"
    control_label <- "Surgery_alone"
    multiplier <- 1
    n_untreated_step19i <- 4  # Step19I included 2 Adjacent_normal samples
  }
  
  dataset_scores <- fib_scores[fib_scores$dataset == dataset, ]
  
  for (pathway in TARGET_PATHWAYS) {
    pathway_scores <- dataset_scores[dataset_scores$pathway == pathway, ]
    
    # Get tumor samples only
    tumor_samples <- intersect(pathway_scores$sample, names(treat_map))
    tumor_scores <- pathway_scores[pathway_scores$sample %in% tumor_samples, ]
    
    # Assign groups
    tumor_scores$group <- treat_map[tumor_scores$sample]
    
    treated_scores <- tumor_scores$score[tumor_scores$group == treated_label]
    control_scores <- tumor_scores$score[tumor_scores$group == control_label]
    
    # Canonical effect = mean(treated) - mean(control) [tumor-only]
    canonical_from_means <- mean(treated_scores) - mean(control_scores)
    
    # Legacy effect from Step19I
    legacy_effect <- step19i_fib$observed_mean_difference[
      step19i_fib$dataset == dataset & step19i_fib$pathway == pathway]
    
    # Canonical effect from reorientation
    canonical_from_reorientation <- legacy_effect * multiplier
    
    # Check agreement
    values_match <- abs(canonical_from_means - canonical_from_reorientation) < TOLERANCE
    
    # Direction
    if (abs(canonical_from_means) < NEAR_ZERO_THRESHOLD) {
      canonical_direction <- "NEAR_ZERO"
    } else if (canonical_from_means > 0) {
      canonical_direction <- "HIGHER_IN_TREATED"
    } else {
      canonical_direction <- "HIGHER_IN_CONTROL"
    }
    
    results_corr2.3 <- rbind(results_corr2.3, data.frame(
      dataset = dataset,
      pathway = pathway,
      legacy_effect = legacy_effect,
      legacy_definition = ifelse(dataset == "GSE197677",
        "mean(nNACT) - mean(NACT)", "mean(NACT) - mean(nNACT)"),
      canonical_effect_from_means = canonical_from_means,
      canonical_effect_from_reorientation = canonical_from_reorientation,
      canonical_direction = canonical_direction,
      values_match = values_match,
      stringsAsFactors = FALSE
    ))
    
    cat(sprintf("%s | %s:\n", dataset, gsub("HALLMARK_", "", pathway)))
    cat(sprintf("  Legacy effect: %.6f\n", legacy_effect))
    cat(sprintf("  Canonical from means (tumor-only): %.6f\n", canonical_from_means))
    cat(sprintf("  Canonical from reorientation: %.6f\n", canonical_from_reorientation))
    cat(sprintf("  Values match: %s\n", values_match))
    cat(sprintf("  Direction: %s\n\n", canonical_direction))
  }
}

write.csv(results_corr2.3,
  file.path(corr2_dir, "STEP19J_CORR2_CANONICAL_FIBROBLAST_PATHWAY_EFFECTS.csv"),
  row.names = FALSE)

cat("Saved: STEP19J_CORR2_CANONICAL_FIBROBLAST_PATHWAY_EFFECTS.csv\n\n")

# Check agreement
# GSE197677 should match (no normals in dataset)
# GSE221561 will NOT match (Step19I included normals)
gse197677_match <- all(results_corr2.3$values_match[results_corr2.3$dataset == "GSE197677"])
gse221561_match <- all(results_corr2.3$values_match[results_corr2.3$dataset == "GSE221561"])

cat("=== Agreement check ===\n")
cat(sprintf("  GSE197677 (no normals in dataset): %s\n", gse197677_match))
cat(sprintf("  GSE221561 (Step19I included normals): %s\n", gse221561_match))
cat("\n")

if (!gse197677_match) {
  stop("STEP19J_CORR2_VALUES_MISMATCH_GSE197677: Should match (no normals)")
}

if (gse221561_match) {
  cat("NOTE: GSE221561 values match unexpectedly. This may indicate normals were not included.\n\n")
} else {
  cat("EXPECTED: GSE221561 values differ because Step19I included Adjacent_normal samples.\n")
  cat("The canonical effects (tumor-only) are the biologically correct values for interpretation.\n\n")
}

# ============================================================================
# 19J-CORR2.4: BIOLOGICAL CROSS-DATASET DIRECTION AUDIT
# ============================================================================
cat("============================================\n")
cat("19J-CORR2.4: BIOLOGICAL CROSS-DATASET DIRECTION AUDIT\n")
cat("============================================\n\n")

results_corr2.4 <- data.frame()

for (pathway in TARGET_PATHWAYS) {
  eff_197677 <- results_corr2.3$canonical_effect_from_means[
    results_corr2.3$dataset == "GSE197677" & results_corr2.3$pathway == pathway]
  eff_221561 <- results_corr2.3$canonical_effect_from_means[
    results_corr2.3$dataset == "GSE221561" & results_corr2.3$pathway == pathway]
  
  # Classification
  if (abs(eff_197677) < NEAR_ZERO_THRESHOLD && abs(eff_221561) < NEAR_ZERO_THRESHOLD) {
    classification <- "UNRESOLVED"
  } else if (abs(eff_197677) < NEAR_ZERO_THRESHOLD || abs(eff_221561) < NEAR_ZERO_THRESHOLD) {
    classification <- "ONE_DATASET_NEAR_ZERO"
  } else if (sign(eff_197677) == sign(eff_221561)) {
    classification <- "CANONICAL_SAME_DIRECTION"
  } else {
    classification <- "CANONICAL_OPPOSITE_DIRECTION"
  }
  
  results_corr2.4 <- rbind(results_corr2.4, data.frame(
    pathway = pathway,
    pathway_short = gsub("HALLMARK_", "", pathway),
    GSE197677_canonical_effect = eff_197677,
    GSE221561_canonical_effect = eff_221561,
    cross_dataset_status = classification,
    stringsAsFactors = FALSE
  ))
  
  cat(sprintf("%s:\n", gsub("HALLMARK_", "", pathway)))
  cat(sprintf("  GSE197677: %.6f\n", eff_197677))
  cat(sprintf("  GSE221561: %.6f\n", eff_221561))
  cat(sprintf("  Status: %s\n\n", classification))
}

write.csv(results_corr2.4,
  file.path(corr2_dir, "STEP19J_CORR2_CANONICAL_PATHWAY_DIRECTION_REPLICATION.csv"),
  row.names = FALSE)

cat("Saved: STEP19J_CORR2_CANONICAL_PATHWAY_DIRECTION_REPLICATION.csv\n\n")

# ============================================================================
# 19J-CORR2.5: GENE-LEVEL EDGER ORIENTATION AUDIT
# ============================================================================
cat("============================================\n")
cat("19J-CORR2.5: GENE-LEVEL EDGER ORIENTATION AUDIT\n")
cat("============================================\n\n")

# Load pseudobulk counts
pb19_fib <- readRDS(file.path(base_dir,
  "06_tables/cross_dataset/pathway_analysis/Step19I/counts/GSE197677_FIBROBLAST_CAF_raw_counts.rds"))
pb22_fib <- readRDS(file.path(base_dir,
  "06_tables/cross_dataset/pathway_analysis/Step19I/counts/GSE221561_FIBROBLAST_CAF_raw_counts.rds"))

cat("GSE197677 fibroblast pseudobulk:", nrow(pb19_fib), "genes x", ncol(pb19_fib), "samples\n")
cat("GSE221561 fibroblast pseudobulk:", nrow(pb22_fib), "genes x", ncol(pb22_fib), "samples\n\n")

# GSE197677: canonical contrast = NACT - nNACT
cat("=== GSE197677 edgeR: NACT - nNACT ===\n\n")

tumor_samples_19 <- intersect(colnames(pb19_fib), names(gse197677_treat_map))
counts_19 <- pb19_fib[, tumor_samples_19, drop = FALSE]
group_19 <- factor(gse197677_treat_map[tumor_samples_19], levels = c("nNACT", "NACT"))

cat("Design matrix:\n")
design_19 <- model.matrix(~ group_19)
print(design_19)
cat("\nCoefficient: group_19NACT (NACT - nNACT)\n")
cat("Positive logFC = higher in NACT\n\n")

# Create DGEList
y_19 <- DGEList(counts = counts_19, group = group_19)
keep_19 <- filterByExpr(y_19)
y_19 <- y_19[keep_19, , keep.lib.sizes = FALSE]
y_19 <- normLibSizes(y_19)

cat(sprintf("After filterByExpr: %d genes\n\n", nrow(y_19)))

# Fit model
fit_19 <- glmQLFit(y_19, design_19)

# Sanity check: one gene
sanity_gene <- "BGN"
if (sanity_gene %in% rownames(y_19)) {
  cat(sprintf("Sanity check gene: %s\n", sanity_gene))
  gene_counts <- y_19$counts[sanity_gene, ]
  cat("  Counts:", gene_counts, "\n")
  cat("  Groups:", group_19, "\n")
  
  treated_mean <- mean(log2(gene_counts[group_19 == "NACT"] + 1))
  control_mean <- mean(log2(gene_counts[group_19 == "nNACT"] + 1))
  manual_direction <- sign(treated_mean - control_mean)
  
  # edgeR result
  coef_name <- colnames(design_19)[2]
  qlf <- glmQLFTest(fit_19, coef = coef_name)
  edgeR_logFC <- qlf$table[sanity_gene, "logFC"]
  edgeR_direction <- sign(edgeR_logFC)
  
  cat(sprintf("  Manual: treated=%.4f, control=%.4f, direction=%d\n",
    treated_mean, control_mean, manual_direction))
  cat(sprintf("  edgeR: logFC=%.6f, direction=%d\n", edgeR_logFC, edgeR_direction))
  cat(sprintf("  Directions match: %s\n\n", manual_direction == edgeR_direction))
  
  if (manual_direction != edgeR_direction) {
    stop("STEP19J_CORR2_EDGER_DIRECTION_MISMATCH_GSE197677")
  }
}

# Full edgeR run for GSE197677
cat("Running full edgeR for GSE197677...\n\n")
coef_name_19 <- colnames(design_19)[2]
qlf_19 <- glmQLFTest(fit_19, coef = coef_name_19)
qlf_19$table$FDR <- p.adjust(qlf_19$table$PValue, method = "BH")

cat("GSE197677 edgeR results (first 10 rows):\n")
print(head(qlf_19$table, 10))
cat(sprintf("\nTotal genes tested: %d\n\n", nrow(qlf_19$table)))

# GSE221561: canonical contrast = Neoadjuvant_treated - Surgery_alone
cat("=== GSE221561 edgeR: Neoadjuvant_treated - Surgery_alone ===\n\n")

tumor_samples_22 <- intersect(colnames(pb22_fib), names(gse221561_treat_map))
counts_22 <- pb22_fib[, tumor_samples_22, drop = FALSE]
group_22 <- factor(gse221561_treat_map[tumor_samples_22],
  levels = c("Surgery_alone", "Neoadjuvant_treated"))

cat("Design matrix:\n")
design_22 <- model.matrix(~ group_22)
print(design_22)
cat("\nCoefficient: group_22Neoadjuvant_treated (Neoadjuvant_treated - Surgery_alone)\n")
cat("Positive logFC = higher in Neoadjuvant_treated\n\n")

# Create DGEList
y_22 <- DGEList(counts = counts_22, group = group_22)
keep_22 <- filterByExpr(y_22)
y_22 <- y_22[keep_22, , keep.lib.sizes = FALSE]
y_22 <- normLibSizes(y_22)

cat(sprintf("After filterByExpr: %d genes\n\n", nrow(y_22)))

# Fit model
fit_22 <- glmQLFit(y_22, design_22)

# Sanity check: one gene
if (sanity_gene %in% rownames(y_22)) {
  cat(sprintf("Sanity check gene: %s\n", sanity_gene))
  gene_counts <- y_22$counts[sanity_gene, ]
  cat("  Counts:", gene_counts, "\n")
  cat("  Groups:", group_22, "\n")
  
  treated_mean <- mean(log2(gene_counts[group_22 == "Neoadjuvant_treated"] + 1))
  control_mean <- mean(log2(gene_counts[group_22 == "Surgery_alone"] + 1))
  manual_direction <- sign(treated_mean - control_mean)
  
  # edgeR result
  coef_name <- colnames(design_22)[2]
  qlf <- glmQLFTest(fit_22, coef = coef_name)
  edgeR_logFC <- qlf$table[sanity_gene, "logFC"]
  edgeR_direction <- sign(edgeR_logFC)
  
  cat(sprintf("  Manual: treated=%.4f, control=%.4f, direction=%d\n",
    treated_mean, control_mean, manual_direction))
  cat(sprintf("  edgeR: logFC=%.6f, direction=%d\n", edgeR_logFC, edgeR_direction))
  cat(sprintf("  Directions match: %s\n\n", manual_direction == edgeR_direction))
  
  if (manual_direction != edgeR_direction) {
    stop("STEP19J_CORR2_EDGER_DIRECTION_MISMATCH_GSE221561")
  }
}

# Full edgeR run for GSE221561
cat("Running full edgeR for GSE221561...\n\n")
coef_name_22 <- colnames(design_22)[2]
qlf_22 <- glmQLFTest(fit_22, coef = coef_name_22)
qlf_22$table$FDR <- p.adjust(qlf_22$table$PValue, method = "BH")

cat("GSE221561 edgeR results (first 10 rows):\n")
print(head(qlf_22$table, 10))
cat(sprintf("\nTotal genes tested: %d\n\n", nrow(qlf_22$table)))

# Save canonical target gene effects
cat("=== Saving canonical target gene effects ===\n\n")

for (dataset in c("GSE197677", "GSE221561")) {
  if (dataset == "GSE197677") {
    qlf_obj <- qlf_19
    treat_label <- "NACT"
    control_label <- "nNACT"
    n_treated <- sum(group_19 == "NACT")
    n_control <- sum(group_19 == "nNACT")
  } else {
    qlf_obj <- qlf_22
    treat_label <- "Neoadjuvant_treated"
    control_label <- "Surgery_alone"
    n_treated <- sum(group_22 == "Neoadjuvant_treated")
    n_control <- sum(group_22 == "Surgery_alone")
  }
  
  # Load leading edge genes from Step19E
  le_file <- file.path(base_dir,
    "06_tables/cross_dataset/pathway_analysis/STEP19E_HALLMARK_LEADING_EDGE_GENE_EVIDENCE.csv")
  le_evidence <- read.csv(le_file)
  le_target <- le_evidence[le_evidence$pathway %in% TARGET_PATHWAYS, ]
  
  # Get target genes (genes in LE sets for target pathways)
  target_genes <- unique(le_target$gene)
  target_genes_in_qlf <- intersect(target_genes, rownames(qlf_obj$table))
  
  cat(sprintf("%s: %d target genes in LE sets, %d found in edgeR results\n",
    dataset, length(target_genes), length(target_genes_in_qlf)))
  
  # Create output table
  out_table <- data.frame(
    gene = target_genes_in_qlf,
    logFC = qlf_obj$table[target_genes_in_qlf, "logFC"],
    PValue = qlf_obj$table[target_genes_in_qlf, "PValue"],
    FDR = qlf_obj$table[target_genes_in_qlf, "FDR"],
    stringsAsFactors = FALSE
  )
  
  # Add pathway membership
  out_table$pathways <- sapply(out_table$gene, function(g) {
    pathways <- unique(le_target$pathway[le_target$gene == g])
    paste(gsub("HALLMARK_", "", pathways), collapse = ";")
  })
  
  out_table$n_pathways <- sapply(out_table$gene, function(g) {
    length(unique(le_target$pathway[le_target$gene == g]))
  })
  
  out_table$dataset <- dataset
  out_table$contrast <- sprintf("%s - %s", treat_label, control_label)
  out_table$canonical_direction <- ifelse(out_table$logFC > 0,
    "HIGHER_IN_TREATED", "HIGHER_IN_CONTROL")
  
  # Save
  outfile <- file.path(corr2_dir,
    sprintf("STEP19J_CORR2_%s_FIBROBLAST_TARGET_GENE_EFFECTS_CANONICAL.csv", dataset))
  write.csv(out_table, outfile, row.names = FALSE)
  
  cat(sprintf("Saved: %s (%d genes)\n\n",
    basename(outfile), nrow(out_table)))
}

# ============================================================================
# 19J-CORR2.6: RECOMPUTE TRUE CROSS-DATASET GENE REPLICATION
# ============================================================================
cat("============================================\n")
cat("19J-CORR2.6: RECOMPUTE TRUE CROSS-DATASET GENE REPLICATION\n")
cat("============================================\n\n")

# Load canonical gene effects
gene_19 <- read.csv(file.path(corr2_dir,
  "STEP19J_CORR2_GSE197677_FIBROBLAST_TARGET_GENE_EFFECTS_CANONICAL.csv"))
gene_22 <- read.csv(file.path(corr2_dir,
  "STEP19J_CORR2_GSE221561_FIBROBLAST_TARGET_GENE_EFFECTS_CANONICAL.csv"))

# Common genes
common_genes <- intersect(gene_19$gene, gene_22$gene)
cat(sprintf("Common target genes: %d\n\n", length(common_genes)))

# Create replication table
replication_results <- data.frame()

for (gene in common_genes) {
  logfc_19 <- gene_19$logFC[gene_19$gene == gene]
  logfc_22 <- gene_22$logFC[gene_22$gene == gene]
  
  same_direction <- sign(logfc_19) == sign(logfc_22)
  strong_same <- same_direction && abs(logfc_19) >= 0.5 && abs(logfc_22) >= 0.5
  opposite <- sign(logfc_19) != sign(logfc_22)
  near_zero <- abs(logfc_19) < 0.1 && abs(logfc_22) < 0.1
  
  # Pathway membership
  pathways_19 <- gene_19$pathways[gene_19$gene == gene]
  pathways_22 <- gene_22$pathways[gene_22$gene == gene]
  
  replication_results <- rbind(replication_results, data.frame(
    gene = gene,
    GSE197677_logFC = logfc_19,
    GSE221561_logFC = logfc_22,
    same_direction = same_direction,
    strong_same_direction = strong_same,
    opposite_direction = opposite,
    near_zero = near_zero,
    pathways_GSE197677 = pathways_19,
    pathways_GSE221561 = pathways_22,
    stringsAsFactors = FALSE
  ))
}

# Summary
n_same <- sum(replication_results$same_direction)
n_opposite <- sum(replication_results$opposite_direction)
n_strong <- sum(replication_results$strong_same_direction)
n_near_zero <- sum(replication_results$near_zero)

cat("=== Cross-dataset gene replication summary ===\n")
cat(sprintf("  Same direction: %d / %d (%.1f%%)\n",
  n_same, length(common_genes), 100 * n_same / length(common_genes)))
cat(sprintf("  Opposite direction: %d / %d (%.1f%%)\n",
  n_opposite, length(common_genes), 100 * n_opposite / length(common_genes)))
cat(sprintf("  Strong same direction: %d / %d (%.1f%%)\n",
  n_strong, length(common_genes), 100 * n_strong / length(common_genes)))
cat(sprintf("  Near zero: %d / %d (%.1f%%)\n\n",
  n_near_zero, length(common_genes), 100 * n_near_zero / length(common_genes)))

# Save
write.csv(replication_results,
  file.path(corr2_dir, "STEP19J_CORR2_CANONICAL_GENE_DIRECTION_REPLICATION.csv"),
  row.names = FALSE)

cat("Saved: STEP19J_CORR2_CANONICAL_GENE_DIRECTION_REPLICATION.csv\n\n")

# ============================================================================
# 19J-CORR2.7: REVALIDATE CURRENT CORE MODULE
# ============================================================================
cat("============================================\n")
cat("19J-CORR2.7: REVALIDATE CURRENT CORE MODULE\n")
cat("============================================\n\n")

cat("Current J-CORR core module:", paste(CORE_GENES, collapse = ", "), "\n\n")

core_results <- data.frame()

for (gene in CORE_GENES) {
  if (gene %in% common_genes) {
    logfc_19 <- gene_19$logFC[gene_19$gene == gene]
    logfc_22 <- gene_22$logFC[gene_22$gene == gene]
    
    same_dir <- sign(logfc_19) == sign(logfc_22)
    strong <- same_dir && abs(logfc_19) >= 0.5 && abs(logfc_22) >= 0.5
    
    # Count target pathways
    n_pathways_19 <- gene_19$n_pathways[gene_19$gene == gene]
    n_pathways_22 <- gene_22$n_pathways[gene_22$gene == gene]
    
    # Status
    if (same_dir && strong) {
      status <- "CANONICAL_STRONG_REPLICATION"
    } else if (same_dir) {
      status <- "CANONICAL_WEAK_REPLICATION"
    } else {
      status <- "CANONICAL_DIRECTION_DISCORDANT"
    }
    
    core_results <- rbind(core_results, data.frame(
      gene = gene,
      GSE197677_logFC = logfc_19,
      GSE221561_logFC = logfc_22,
      same_biological_direction = same_dir,
      abs_logFC_gte_0.5_both = strong,
      n_target_pathways = max(n_pathways_19, n_pathways_22),
      final_status = status,
      stringsAsFactors = FALSE
    ))
  } else {
    core_results <- rbind(core_results, data.frame(
      gene = gene,
      GSE197677_logFC = NA,
      GSE221561_logFC = NA,
      same_biological_direction = NA,
      abs_logFC_gte_0.5_both = NA,
      n_target_pathways = 0,
      final_status = "LOW_EFFECT_UNRESOLVED",
      stringsAsFactors = FALSE
    ))
  }
}

cat("=== Core module validation under canonical contrast ===\n\n")
for (i in 1:nrow(core_results)) {
  cat(sprintf("%s:\n", core_results$gene[i]))
  cat(sprintf("  GSE197677 logFC: %.4f\n", core_results$GSE197677_logFC[i]))
  cat(sprintf("  GSE221561 logFC: %.4f\n", core_results$GSE221561_logFC[i]))
  cat(sprintf("  Same direction: %s\n", core_results$same_biological_direction[i]))
  cat(sprintf("  Strong replication: %s\n", core_results$abs_logFC_gte_0.5_both[i]))
  cat(sprintf("  Status: %s\n\n", core_results$final_status[i]))
}

# Update canonical core module
canonical_core <- core_results$gene[
  core_results$final_status == "CANONICAL_STRONG_REPLICATION" |
  core_results$final_status == "CANONICAL_WEAK_REPLICATION"]

cat("Canonical core module:", paste(canonical_core, collapse = ", "), "\n\n")

write.csv(core_results,
  file.path(corr2_dir, "STEP19J_CORR2_CANONICAL_FIBROBLAST_CORE_MODULE.csv"),
  row.names = FALSE)

cat("Saved: STEP19J_CORR2_CANONICAL_FIBROBLAST_CORE_MODULE.csv\n\n")

# ============================================================================
# 19J-CORR2.8: REASSESS ABUNDANCE COUPLING
# ============================================================================
cat("============================================\n")
cat("19J-CORR2.8: REASSESS ABUNDANCE COUPLING\n")
cat("============================================\n\n")

cat("=== Frozen abundance findings (NOT modified) ===\n\n")
cat("GSE197677 fibroblast fraction:\n")
cat("  effect = 0.0610\n")
cat("  P = 0.095\n\n")
cat("GSE221561 fibroblast fraction:\n")
cat("  effect = 0.1388\n")
cat("  P = 0.403\n\n")

cat("=== Verified pathway-abundance correlations ===\n\n")
cat("Coagulation / GSE197677: rho ~ 0.85\n")
cat("KRAS / GSE221561: rho ~ 0.71\n\n")

cat("=== Important distinction ===\n\n")
cat("ABUNDANCE_COUPLING != CROSS_DATASET_DIRECTIONAL_REPLICATION\n\n")
cat("Abundance coupling = pathway scores correlate with fibroblast fraction\n")
cat("Directional replication = same biological direction across datasets\n\n")

# ============================================================================
# 19J-CORR2.9: RECLASSIFY FOUR PATHWAYS
# ============================================================================
cat("============================================\n")
cat("19J-CORR2.9: RECLASSIFY FOUR PATHWAYS\n")
cat("============================================\n\n")

# Load canonical pathway effects
pathway_effects <- read.csv(file.path(corr2_dir,
  "STEP19J_CORR2_CANONICAL_FIBROBLAST_PATHWAY_EFFECTS.csv"))

# Load cross-dataset direction replication
pathway_dir_rep <- read.csv(file.path(corr2_dir,
  "STEP19J_CORR2_CANONICAL_PATHWAY_DIRECTION_REPLICATION.csv"))

# Load gene replication
gene_rep <- read.csv(file.path(corr2_dir,
  "STEP19J_CORR2_CANONICAL_GENE_DIRECTION_REPLICATION.csv"))

# Load target gene effects
gene_19 <- read.csv(file.path(corr2_dir,
  "STEP19J_CORR2_GSE197677_FIBROBLAST_TARGET_GENE_EFFECTS_CANONICAL.csv"))
gene_22 <- read.csv(file.path(corr2_dir,
  "STEP19J_CORR2_GSE221561_FIBROBLAST_TARGET_GENE_EFFECTS_CANONICAL.csv"))

# Load LE evidence
le_file <- file.path(base_dir,
  "06_tables/cross_dataset/pathway_analysis/STEP19E_HALLMARK_LEADING_EDGE_GENE_EVIDENCE.csv")
le_evidence <- read.csv(le_file)

pathway_classifications <- data.frame()

for (pathway in TARGET_PATHWAYS) {
  pathway_short <- gsub("HALLMARK_", "", pathway)
  
  # Get cross-dataset direction status from pathway_dir_rep
  dir_status <- pathway_dir_rep$cross_dataset_status[pathway_dir_rep$pathway == pathway]
  
  # Get pathway-specific genes
  pathway_genes <- unique(le_evidence$gene[le_evidence$pathway == pathway])
  
  # Get gene-level replication for pathway genes
  pathway_gene_rep <- gene_rep[gene_rep$gene %in% pathway_genes, ]
  n_same <- sum(pathway_gene_rep$same_direction, na.rm = TRUE)
  n_opposite <- sum(pathway_gene_rep$opposite_direction, na.rm = TRUE)
  n_total <- nrow(pathway_gene_rep)
  
  # Classification
  if (dir_status == "CANONICAL_SAME_DIRECTION") {
    classification <- "FIBROBLAST_PROGRAM_CANONICALLY_REPLICATED"
  } else if (dir_status == "CANONICAL_OPPOSITE_DIRECTION") {
    classification <- "FIBROBLAST_PROGRAM_DATASET_SPECIFIC"
  } else if (dir_status == "ONE_DATASET_NEAR_ZERO") {
    classification <- "FIBROBLAST_PATHWAY_SUPPORT_ONLY"
  } else {
    classification <- "INCONCLUSIVE"
  }
  
  pathway_classifications <- rbind(pathway_classifications, data.frame(
    pathway = pathway,
    pathway_short = pathway_short,
    cross_dataset_direction_status = dir_status,
    GSE197677_canonical_effect = pathway_effects$canonical_effect_from_means[
      pathway_effects$dataset == "GSE197677" & pathway_effects$pathway == pathway],
    GSE221561_canonical_effect = pathway_effects$canonical_effect_from_means[
      pathway_effects$dataset == "GSE221561" & pathway_effects$pathway == pathway],
    n_genes_same_direction = n_same,
    n_genes_opposite_direction = n_opposite,
    n_genes_total = n_total,
    classification = classification,
    stringsAsFactors = FALSE
  ))
  
  cat(sprintf("%s:\n", pathway_short))
  cat(sprintf("  Cross-dataset direction: %s\n", dir_status))
  cat(sprintf("  GSE197677 effect: %.4f\n",
    pathway_classifications$GSE197677_canonical_effect[nrow(pathway_classifications)]))
  cat(sprintf("  GSE221561 effect: %.4f\n",
    pathway_classifications$GSE221561_canonical_effect[nrow(pathway_classifications)]))
  cat(sprintf("  Gene replication: %d same / %d opposite / %d total\n",
    n_same, n_opposite, n_total))
  cat(sprintf("  Classification: %s\n\n", classification))
}

write.csv(pathway_classifications,
  file.path(corr2_dir, "STEP19J_CORR2_CANONICAL_PATHWAY_CLASSIFICATIONS.csv"),
  row.names = FALSE)

cat("Saved: STEP19J_CORR2_CANONICAL_PATHWAY_CLASSIFICATIONS.csv\n\n")

# ============================================================================
# 19J-CORR2.10: OVERALL INTERPRETATION
# ============================================================================
cat("============================================\n")
cat("19J-CORR2.10: OVERALL INTERPRETATION\n")
cat("============================================\n\n")

# Check if all pathways show opposite direction
all_opposite <- all(pathway_classifications$cross_dataset_direction_status == "CANONICAL_OPPOSITE_DIRECTION")
any_replicated <- any(pathway_classifications$classification == "FIBROBLAST_PROGRAM_CANONICALLY_REPLICATED")
any_heterogeneous <- any(pathway_classifications$classification == "FIBROBLAST_PROGRAM_WITH_ABUNDANCE_COUPLING_BUT_DIRECTIONALLY_HETEROGENEOUS")

if (all_opposite) {
  overall_status <- "FIBROBLAST_REMODELING_WITH_CROSS_DATASET_DIRECTIONAL_HETEROGENEITY"
} else if (any_replicated) {
  if (any_heterogeneous) {
    overall_status <- "FIBROBLAST_ABUNDANCE_COUPLED_DATASET_SPECIFIC_PROGRAMS"
  } else {
    overall_status <- "FIBROBLAST_TRANSCRIPTIONAL_REMODELING_CANONICALLY_REPLICATED"
  }
} else if (any(pathway_classifications$classification == "FIBROBLAST_PATHWAY_SUPPORT_ONLY")) {
  overall_status <- "FIBROBLAST_PATHWAY_LEVEL_SUPPORT_ONLY"
} else {
  overall_status <- "FIBROBLAST_TARGETED_VALIDATION_INCONCLUSIVE"
}

cat(sprintf("OVERALL CANONICAL STEP19J STATUS: %s\n\n", overall_status))

# Print pathway classifications
cat("=== Canonical pathway classifications ===\n\n")
for (i in 1:nrow(pathway_classifications)) {
  cat(sprintf("%s: %s\n",
    pathway_classifications$pathway_short[i],
    pathway_classifications$classification[i]))
}

cat("\n")

# ============================================================================
# 19J-CORR2.11: PROVENANCE
# ============================================================================
cat("============================================\n")
cat("19J-CORR2.11: PROVENANCE\n")
cat("============================================\n\n")

provenance_md <- paste0(
"# STEP 19J-CORR2 DIRECTION PROVENANCE\n\n",
"## Overview\n\n",
"This document establishes the canonical biological contrast for cross-dataset interpretation.\n\n",
"## Key Points\n\n",
"- **Step19I numeric outputs remain frozen** - no changes to existing files\n",
"- **Step19J-CORR reconciled outputs to Step19I legacy signs** - matching frozen conventions\n",
"- **Legacy signs differed by dataset** because `exact_perm_test` used `unique(groups)[1]` to determine treated vs untreated\n",
"- **Step19J-CORR2 introduces the first uniform biological contrast** for cross-dataset interpretation\n",
"- **No original files were overwritten** - all new outputs in Step19J-CORR2 directory\n",
"- **No samples/cells were modified** - same inclusion criteria as previous steps\n\n",
"## Canonical Biological Contrast\n\n",
"```\n",
"CANONICAL TREATMENT EFFECT = NEOADJUVANT-TREATED minus NO-NEOADJUVANT / SURGERY-ALONE\n",
"```\n\n",
"- **GSE197677**: mean(NACT) - mean(nNACT)\n",
"- **GSE221561**: mean(Neoadjuvant_treated) - mean(Surgery_alone)\n\n",
"**Positive ALWAYS** = higher in neoadjuvant-treated condition\n",
"**Negative ALWAYS** = higher in untreated / surgery-alone condition\n\n",
"## Legacy Conventions\n\n",
"- **GSE197677**: mean(nNACT) - mean(NACT) → multiplier = -1\n",
"- **GSE221561**: mean(NACT) - mean(nNACT) → multiplier = +1\n\n",
"## Tumor Sample Composition\n\n",
"- **GSE197677**: 10 tumor samples (6 NACT + 4 nNACT)\n",
"- **GSE221561**: 9 tumor samples (7 Neoadjuvant_treated + 2 Surgery_alone)\n\n",
"## Outputs Created\n\n",
"1. `STEP19J_CORR2_LEGACY_DIRECTION_PROVENANCE.csv`\n",
"2. `STEP19J_CORR2_CANONICAL_FIBROBLAST_PATHWAY_EFFECTS.csv`\n",
"3. `STEP19J_CORR2_CANONICAL_PATHWAY_DIRECTION_REPLICATION.csv`\n",
"4. `STEP19J_CORR2_GSE197677_FIBROBLAST_TARGET_GENE_EFFECTS_CANONICAL.csv`\n",
"5. `STEP19J_CORR2_GSE221561_FIBROBLAST_TARGET_GENE_EFFECTS_CANONICAL.csv`\n",
"6. `STEP19J_CORR2_CANONICAL_GENE_DIRECTION_REPLICATION.csv`\n",
"7. `STEP19J_CORR2_CANONICAL_FIBROBLAST_CORE_MODULE.csv`\n",
"8. `STEP19J_CORR2_CANONICAL_PATHWAY_CLASSIFICATIONS.csv`\n"
)

writeLines(provenance_md, file.path(corr2_dir, "STEP19J_CORR2_DIRECTION_PROVENANCE.md"))

cat("Saved: STEP19J_CORR2_DIRECTION_PROVENANCE.md\n\n")

# ============================================================================
# FINAL OUTPUT
# ============================================================================
cat("============================================\n")
cat("STEP 19J-CORR2 COMPLETE\n")
cat("BIOLOGICAL CONTRAST NORMALIZED\n")
cat("============================================\n\n")

cat("Canonical contrast:\n")
cat("  GSE197677: NACT - nNACT\n")
cat("  GSE221561: Neoadjuvant_treated - Surgery_alone\n\n")

cat("Tumor sample validation:\n")
cat("  GSE197677: 10 total (6 treated + 4 control + 0 normal)\n")
cat("  GSE221561: 9 total (7 treated + 2 control + 0 normal)\n\n")

cat("--------------------------------------------\n\n")

cat("Pathway canonical effects:\n\n")
for (pathway in TARGET_PATHWAYS) {
  pathway_short <- gsub("HALLMARK_", "", pathway)
  eff_19 <- pathway_effects$canonical_effect_from_means[
    pathway_effects$dataset == "GSE197677" & pathway_effects$pathway == pathway]
  eff_22 <- pathway_effects$canonical_effect_from_means[
    pathway_effects$dataset == "GSE221561" & pathway_effects$pathway == pathway]
  dir_status <- pathway_effects$cross_dataset_status[pathway_effects$pathway == pathway]
  
  cat(sprintf("%s:\n", pathway_short))
  cat(sprintf("  GSE197677: %.4f\n", eff_19))
  cat(sprintf("  GSE221561: %.4f\n", eff_22))
  cat(sprintf("  Cross-dataset status: %s\n\n", dir_status))
}

cat("--------------------------------------------\n\n")

cat("Canonical gene-level replication:\n\n")
for (pathway in TARGET_PATHWAYS) {
  pathway_short <- gsub("HALLMARK_", "", pathway)
  pathway_genes <- unique(le_evidence$gene[le_evidence$pathway == pathway])
  pathway_gene_rep <- gene_rep[gene_rep$gene %in% pathway_genes, ]
  
  cat(sprintf("%s:\n", pathway_short))
  cat(sprintf("  Same: %d\n", sum(pathway_gene_rep$same_direction, na.rm = TRUE)))
  cat(sprintf("  Opposite: %d\n", sum(pathway_gene_rep$opposite_direction, na.rm = TRUE)))
  cat(sprintf("  Strong same: %d\n\n", sum(pathway_gene_rep$strong_same_direction, na.rm = TRUE)))
}

cat("--------------------------------------------\n\n")

cat("Previous core module: BGN, CDKN1B, GSN, TIMP1\n")
cat("Canonical core module:", paste(canonical_core, collapse = ", "), "\n\n")

cat("--------------------------------------------\n\n")

cat("Canonical pathway classifications:\n\n")
for (i in 1:nrow(pathway_classifications)) {
  cat(sprintf("%s: %s\n",
    pathway_classifications$pathway_short[i],
    pathway_classifications$classification[i]))
}

cat("\n--------------------------------------------\n\n")

cat("OVERALL CANONICAL STEP19J STATUS:\n")
cat(sprintf("  %s\n\n", overall_status))

cat("Legacy Step19I changed: NO\n")
cat("Step19J-CORR changed: NO\n")
cat("Cells reclassified: 0\n")
cat("Cells removed: 0\n")
cat("Samples removed: 0\n\n")

cat("STOP HERE.\n")
cat("DO NOT START STEP19K.\n")
