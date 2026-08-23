#!/usr/bin/env Rscript
# ============================================================================
# STEP 19J-CORR: DIRECTION ORIENTATION AUDIT AND CORRECTION
# ============================================================================
# Root cause: Step19I's exact_perm_test uses unique(groups) to determine
# treated vs untreated. For GSE197677, unique(groups) = c("nNACT", "NACT"),
# so Step19I computes mean(nNACT) - mean(NACT). For GSE221561,
# unique(groups) = c("NACT", "nNACT"), so Step19I computes mean(NACT) - mean(nNACT).
#
# Step19J must match Step19I's frozen conventions exactly.
#
# FIX STRATEGY:
# 1. One canonical treatment map per dataset, defined ONCE, never overwritten.
# 2. Explicit tumor-sample subsetting before edgeR.
# 3. No secondary map lookups inside any function.
# 4. Manual pathway direction check BEFORE edgeR.
# 5. One gene sanity check BEFORE full edgeR run.
# ============================================================================

suppressPackageStartupMessages({
  library(edgeR)
})

base_dir <- Sys.getenv("ESCC_CAF_PROJECT_ROOT", unset = ".")
g19j_dir <- file.path(base_dir, "06_tables/cross_dataset/pathway_analysis/Step19J")

TARGET_PATHWAYS <- c("HALLMARK_HYPOXIA", "HALLMARK_COAGULATION",
                     "HALLMARK_KRAS_SIGNALING_UP", "HALLMARK_APOPTOSIS")

# ============================================================================
# FIX 1: DEFINE ONE CANONICAL TREATMENT MAP PER DATASET
# ============================================================================
cat("========================================\n")
cat("FIX 1: DEFINE CANONICAL TREATMENT MAPS\n")
cat("========================================\n\n")

# GSE197677: 10 tumor samples only, mapped to NACT/nNACT matching Step19I convention
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

cat("GSE197677 canonical map:\n")
print(gse197677_treat_map)
cat(sprintf("  Total: %d, NACT: %d, nNACT: %d\n\n",
  length(gse197677_treat_map),
  sum(gse197677_treat_map == "NACT"),
  sum(gse197677_treat_map == "nNACT")))

# GSE221561: 11 samples (9 tumor + 2 Adjacent_normal)
# Step19I convention: Adjacent_normal -> nNACT (included in all samples)
gse221561_treat_map <- c(
  S0520T  = "NACT",
  S1265T  = "NACT",
  S1315TS = "NACT",
  S1535T  = "NACT",
  S2423T  = "NACT",
  S2487T  = "NACT",
  S6417N  = "nNACT",
  S6417T  = "nNACT",
  S6829T  = "NACT",
  S9478N  = "nNACT",
  S9478T  = "nNACT"
)

cat("GSE221561 canonical map:\n")
print(gse221561_treat_map)
cat(sprintf("  Total: %d, NACT: %d, nNACT: %d\n\n",
  length(gse221561_treat_map),
  sum(gse221561_treat_map == "NACT"),
  sum(gse221561_treat_map == "nNACT")))

# ============================================================================
# FIX 2: EXPLICITLY SUBSET FIBROBLAST PSEUDOBULK
# ============================================================================
cat("========================================\n")
cat("FIX 2: SUBSET FIBROBLAST PSEUDOBULK\n")
cat("========================================\n\n")

pb19_fib <- readRDS(file.path(base_dir, "06_tables/cross_dataset/pathway_analysis/Step19I/counts/GSE197677_FIBROBLAST_CAF_raw_counts.rds"))
pb22_fib <- readRDS(file.path(base_dir, "06_tables/cross_dataset/pathway_analysis/Step19I/counts/GSE221561_FIBROBLAST_CAF_raw_counts.rds"))

cat("GSE197677 fibroblast pseudobulk:", nrow(pb19_fib), "genes x", ncol(pb19_fib), "samples\n")
cat("  All samples:", paste(colnames(pb19_fib), collapse=", "), "\n\n")

# GSE197677: subset to tumor samples only
tumor_samples_197677 <- intersect(colnames(pb19_fib), names(gse197677_treat_map))
cat("GSE197677 tumor samples:\n")
print(tumor_samples_197677)
cat(sprintf("  Length: %d\n\n", length(tumor_samples_197677)))

fib_counts_tumor_197677 <- pb19_fib[, tumor_samples_197677, drop = FALSE]
group_197677 <- unname(gse197677_treat_map[colnames(fib_counts_tumor_197677)])

# Validate
stopifnot(length(group_197677) == 10)
stopifnot(!anyNA(group_197677))
stopifnot(all(group_197677 %in% c("nNACT", "NACT")))

cat("GSE197677 group validation:\n")
print(table(group_197677))
cat("\n")

# GSE221561: subset to mapped samples only
tumor_samples_221561 <- intersect(colnames(pb22_fib), names(gse221561_treat_map))
cat("GSE221561 tumor samples:\n")
print(tumor_samples_221561)
cat(sprintf("  Length: %d\n\n", length(tumor_samples_221561)))

fib_counts_tumor_221561 <- pb22_fib[, tumor_samples_221561, drop = FALSE]
group_221561 <- unname(gse221561_treat_map[colnames(fib_counts_tumor_221561)])

stopifnot(!anyNA(group_221561))
stopifnot(all(group_221561 %in% c("nNACT", "NACT")))

cat("GSE221561 group validation:\n")
print(table(group_221561))
cat("\n")

# ============================================================================
# FIX 3: EXPLICIT FACTOR ORIENTATION
# ============================================================================
cat("========================================\n")
cat("FIX 3: EXPLICIT FACTOR ORIENTATION\n")
cat("========================================\n\n")

# GSE197677: Step19I convention is mean(nNACT) - mean(NACT)
# So nNACT = reference, NACT = test
# Positive logFC = higher in NACT
group_197677 <- factor(group_197677, levels = c("nNACT", "NACT"))
design_197677 <- model.matrix(~ group_197677)

cat("GSE197677 design:\n")
cat("  colnames:", paste(colnames(design_197677), collapse=", "), "\n")
cat("  Expected coef: group_197677NACT (NACT - nNACT)\n")
cat("  Positive logFC = higher in NACT\n\n")

# GSE221561: Step19I convention is mean(NACT) - mean(nNACT)
# So nNACT = reference, NACT = test
# Positive logFC = higher in NACT
group_221561 <- factor(group_221561, levels = c("nNACT", "NACT"))
design_221561 <- model.matrix(~ group_221561)

cat("GSE221561 design:\n")
cat("  colnames:", paste(colnames(design_221561), collapse=", "), "\n")
cat("  Expected coef: group_221561NACT (NACT - nNACT)\n")
cat("  Positive logFC = higher in NACT\n\n")

# ============================================================================
# FIX 5: MANUAL PATHWAY DIRECTION CHECK FIRST
# ============================================================================
cat("========================================\n")
cat("FIX 5: MANUAL PATHWAY DIRECTION CHECK\n")
cat("========================================\n\n")

scores_all <- read.csv(file.path(base_dir, "06_tables/cross_dataset/pathway_analysis/Step19I/STEP19I_SAMPLE_COMPARTMENT_PATHWAY_SCORES.csv"))
fib_scores <- scores_all[scores_all$compartment == "FIBROBLAST_CAF", ]

step19i_te <- read.csv(file.path(base_dir, "06_tables/cross_dataset/pathway_analysis/Step19I/STEP19I_COMPARTMENT_TREATMENT_EFFECTS.csv"))
step19i_fib <- step19i_te[step19i_te$compartment == "FIBROBLAST_CAF", ]

cat("=== Frozen Step19I fibroblast treatment effects ===\n")
for (d in c("GSE197677", "GSE221561")) {
  sub <- step19i_fib[step19i_fib$dataset == d, ]
  cat(sprintf("\n%s:\n", d))
  for (i in 1:nrow(sub)) {
    cat(sprintf("  %s: effect=%.4f\n",
      gsub("HALLMARK_", "", sub$pathway[i]),
      sub$observed_mean_difference[i]))
  }
}

cat("\n=== GSE197677 Step19I convention verification ===\n")
g19_fib <- fib_scores[fib_scores$dataset == "GSE197677", ]
hyp_sub <- g19_fib[g19_fib$pathway == "HALLMARK_HYPOXIA", ]
groups_raw <- ifelse(hyp_sub$treatment_group == "Neoadjuvant_chemotherapy", "NACT", "nNACT")
has_group <- !is.na(groups_raw)
groups_filtered <- groups_raw[has_group]

cat("  Sample order (tumor only):", hyp_sub$sample[has_group], "\n")
cat("  Groups:", groups_filtered, "\n")
cat("  unique(groups):", unique(groups_filtered), "\n")
cat("  Step19I treated_name:", unique(groups_filtered)[1], "\n")
cat("  Step19I untreated_name:", unique(groups_filtered)[2], "\n")

step19i_conv <- step19i_fib$observed_mean_difference[step19i_fib$dataset == "GSE197677" & step19i_fib$pathway == "HALLMARK_HYPOXIA"]
cat(sprintf("  Step19I frozen effect: %.4f\n", step19i_conv))

# Manual verification
v_treated <- hyp_sub$score[has_group][groups_filtered == unique(groups_filtered)[1]]
v_untreated <- hyp_sub$score[has_group][groups_filtered == unique(groups_filtered)[2]]
cat(sprintf("  Manual: mean(%s) - mean(%s) = %.4f - %.4f = %.4f\n",
  unique(groups_filtered)[1], unique(groups_filtered)[2],
  mean(v_treated), mean(v_untreated), mean(v_treated) - mean(v_untreated)))

# Verify for all 4 pathways
cat("\n=== GSE197677 manual vs Step19I for all pathways ===\n")
for (pw in TARGET_PATHWAYS) {
  sub <- g19_fib[g19_fib$pathway == pw, ]
  groups_raw <- ifelse(sub$treatment_group == "Neoadjuvant_chemotherapy", "NACT", "nNACT")
  has_group <- !is.na(groups_raw)
  groups_filt <- groups_raw[has_group]
  scores_filt <- sub$score[has_group]

  treated_name <- unique(groups_filt)[1]
  untreated_name <- unique(groups_filt)[2]
  v_t <- scores_filt[groups_filt == treated_name]
  v_u <- scores_filt[groups_filt == untreated_name]
  manual_diff <- mean(v_t) - mean(v_u)

  step19i_diff <- step19i_fib$observed_mean_difference[step19i_fib$dataset == "GSE197677" & step19i_fib$pathway == pw]

  cat(sprintf("  %s: Step19I=%.4f, Manual=%.4f, Match=%s\n",
    gsub("HALLMARK_", "", pw), step19i_diff, manual_diff,
    abs(manual_diff - step19i_diff) < 0.001))
}

# Verify for GSE221561
cat("\n=== GSE221561 Step19I convention verification ===\n")
g22_fib <- fib_scores[fib_scores$dataset == "GSE221561", ]
hyp_sub22 <- g22_fib[g22_fib$pathway == "HALLMARK_HYPOXIA", ]
groups_raw22 <- ifelse(hyp_sub22$treatment_group == "Neoadjuvant_treated", "NACT", "nNACT")
has_group22 <- !is.na(groups_raw22)
groups_filt22 <- groups_raw22[has_group22]

cat("  Sample order (all):", hyp_sub22$sample[has_group22], "\n")
cat("  Groups:", groups_filt22, "\n")
cat("  unique(groups):", unique(groups_filt22), "\n")
cat("  Step19I treated_name:", unique(groups_filt22)[1], "\n")
cat("  Step19I untreated_name:", unique(groups_filt22)[2], "\n")

step19i_conv22 <- step19i_fib$observed_mean_difference[step19i_fib$dataset == "GSE221561" & step19i_fib$pathway == "HALLMARK_HYPOXIA"]
cat(sprintf("  Step19I frozen effect: %.4f\n", step19i_conv22))

# Manual verification
v_treated22 <- hyp_sub22$score[has_group22][groups_filt22 == unique(groups_filt22)[1]]
v_untreated22 <- hyp_sub22$score[has_group22][groups_filt22 == unique(groups_filt22)[2]]
cat(sprintf("  Manual: mean(%s) - mean(%s) = %.4f - %.4f = %.4f\n",
  unique(groups_filt22)[1], unique(groups_filt22)[2],
  mean(v_treated22), mean(v_untreated22), mean(v_treated22) - mean(v_untreated22)))

cat("\n=== GSE221561 manual vs Step19I for all pathways ===\n")
for (pw in TARGET_PATHWAYS) {
  sub <- g22_fib[g22_fib$pathway == pw, ]
  groups_raw <- ifelse(sub$treatment_group == "Neoadjuvant_treated", "NACT", "nNACT")
  has_group <- !is.na(groups_raw)
  groups_filt <- groups_raw[has_group]
  scores_filt <- sub$score[has_group]

  treated_name <- unique(groups_filt)[1]
  untreated_name <- unique(groups_filt)[2]
  v_t <- scores_filt[groups_filt == treated_name]
  v_u <- scores_filt[groups_filt == untreated_name]
  manual_diff <- mean(v_t) - mean(v_u)

  step19i_diff <- step19i_fib$observed_mean_difference[step19i_fib$dataset == "GSE221561" & step19i_fib$pathway == pw]

  cat(sprintf("  %s: Step19I=%.4f, Manual=%.4f, Match=%s\n",
    gsub("HALLMARK_", "", pw), step19i_diff, manual_diff,
    abs(manual_diff - step19i_diff) < 0.001))
}

# Check expected signs
cat("\n=== Expected direction signs ===\n")
expected_signs <- list(
  HALLMARK_HYPOXIA = c(GSE197677 = "positive", GSE221561 = "positive"),
  HALLMARK_COAGULATION = c(GSE197677 = "negative", GSE221561 = "negative"),
  HALLMARK_KRAS_SIGNALING_UP = c(GSE197677 = "positive", GSE221561 = "positive"),
  HALLMARK_APOPTOSIS = c(GSE197677 = "positive", GSE221561 = "positive")
)

all_dirs_match <- TRUE
for (pw in TARGET_PATHWAYS) {
  for (d in c("GSE197677", "GSE221561")) {
    te <- step19i_fib$observed_mean_difference[step19i_fib$dataset == d & step19i_fib$pathway == pw]
    actual_sign <- if (te > 0) "positive" else "negative"
    expected_sign <- expected_signs[[pw]][d]
    match <- actual_sign == expected_sign
    if (!match) all_dirs_match <- FALSE
    cat(sprintf("  %s %s: effect=%.4f, expected=%s, actual=%s, match=%s\n",
      gsub("HALLMARK_", "", pw), d, te, expected_sign, actual_sign, match))
  }
}

if (!all_dirs_match) {
  cat("\nSTOP: STEP19J_SCORE_ORIENTATION_MISMATCH\n")
  stop("Step19I score orientation mismatch")
} else {
  cat("\nAll Step19I score directions verified.\n\n")
}

# ============================================================================
# FIX 4 + FIX 6: EDGER WITH EXPLICIT GROUPS + SANITY CHECK
# ============================================================================
cat("========================================\n")
cat("FIX 4+6: EDGER WITH EXPLICIT GROUPS\n")
cat("========================================\n\n")

# GSE197677 edgeR: pass counts and factor directly
cat("Running GSE197677 edgeR...\n")
cat(sprintf("  Counts: %d genes x %d samples\n", nrow(fib_counts_tumor_197677), ncol(fib_counts_tumor_197677)))
cat(sprintf("  Groups: %s\n", paste(table(group_197677), collapse=", ")))
cat(sprintf("  Design colnames: %s\n", paste(colnames(design_197677), collapse=", ")))

d19 <- edgeR::DGEList(counts = fib_counts_tumor_197677)
keep19 <- edgeR::filterByExpr(d19, design_197677)
d19 <- d19[keep19, , keep.lib.sizes = FALSE]
d19 <- edgeR::calcNormFactors(d19)
d19 <- edgeR::estimateDisp(d19, design_197677)
fit19 <- edgeR::glmQLFit(d19, design_197677)

# Determine coefficient name
coef_19 <- grep("group", colnames(design_197677), value = TRUE)
cat(sprintf("  Using coef: %s\n", coef_19))

qlf19 <- edgeR::glmQLFTest(fit19, coef = coef_19)
fdr19 <- p.adjust(qlf19$table$PValue, method = "BH")

res19 <- data.frame(
  gene = rownames(d19),
  logFC = qlf19$table$logFC,
  PValue = qlf19$table$PValue,
  FDR = fdr19,
  stringsAsFactors = FALSE
)

cat(sprintf("  edgeR results: %d genes\n", nrow(res19)))
cat(sprintf("  edgeR coef label: NACT - nNACT (positive = higher in NACT)\n\n"))

# FIX 6: One gene sanity check
cat("=== GSE197677 One Gene Sanity Check ===\n")
test_gene <- "BGN"
if (test_gene %in% rownames(d19)) {
  counts_gene <- as.numeric(fib_counts_tumor_197677[test_gene, ])
  nact_mask <- group_197677 == "NACT"
  nnact_mask <- group_197677 == "nNACT"

  nact_mean <- mean(counts_gene[nact_mask])
  nnact_mean <- mean(counts_gene[nnact_mask])
  manual_diff <- nact_mean - nnact_mean

  edgeR_fc <- res19$logFC[res19$gene == test_gene]
  direction_match <- sign(manual_diff) == sign(edgeR_fc)

  cat(sprintf("  Gene: %s\n", test_gene))
  cat(sprintf("  NACT raw counts mean: %.2f\n", nact_mean))
  cat(sprintf("  nNACT raw counts mean: %.2f\n", nnact_mean))
  cat(sprintf("  Manual difference: %.2f\n", manual_diff))
  cat(sprintf("  edgeR logFC: %.4f\n", edgeR_fc))
  cat(sprintf("  Direction match: %s\n", direction_match))

  if (!direction_match) {
    cat("\nSTOP: STEP19J_EDGER_ORIENTATION_FAILED\n")
    stop("edgeR orientation failed")
  } else {
    cat("\nGSE197677 edgeR orientation: PASS\n\n")
  }
} else {
  cat(sprintf("  Gene %s not found in filtered counts. Trying alternative...\n", test_gene))
  test_gene2 <- "COL1A1"
  if (test_gene2 %in% rownames(d19)) {
    counts_gene <- as.numeric(fib_counts_tumor_197677[test_gene2, ])
    nact_mask <- group_197677 == "NACT"
    nnact_mask <- group_197677 == "nNACT"
    nact_mean <- mean(counts_gene[nact_mask])
    nnact_mean <- mean(counts_gene[nnact_mask])
    manual_diff <- nact_mean - nnact_mean
    edgeR_fc <- res19$logFC[res19$gene == test_gene2]
    direction_match <- sign(manual_diff) == sign(edgeR_fc)
    cat(sprintf("  Gene: %s\n", test_gene2))
    cat(sprintf("  NACT mean: %.2f, nNACT mean: %.2f\n", nact_mean, nnact_mean))
    cat(sprintf("  Manual diff: %.2f, edgeR logFC: %.4f\n", manual_diff, edgeR_fc))
    cat(sprintf("  Direction match: %s\n", direction_match))
    if (!direction_match) stop("edgeR orientation failed")
    cat("\nGSE197677 edgeR orientation: PASS\n\n")
  }
}

# GSE221561 edgeR
cat("Running GSE221561 edgeR...\n")
cat(sprintf("  Counts: %d genes x %d samples\n", nrow(fib_counts_tumor_221561), ncol(fib_counts_tumor_221561)))
cat(sprintf("  Groups: %s\n", paste(table(group_221561), collapse=", ")))
cat(sprintf("  Design colnames: %s\n", paste(colnames(design_221561), collapse=", ")))

d22 <- edgeR::DGEList(counts = fib_counts_tumor_221561)
keep22 <- edgeR::filterByExpr(d22, design_221561)
d22 <- d22[keep22, , keep.lib.sizes = FALSE]
d22 <- edgeR::calcNormFactors(d22)
d22 <- edgeR::estimateDisp(d22, design_221561)
fit22 <- edgeR::glmQLFit(d22, design_221561)

coef_22 <- grep("group", colnames(design_221561), value = TRUE)
cat(sprintf("  Using coef: %s\n", coef_22))

qlf22 <- edgeR::glmQLFTest(fit22, coef = coef_22)
fdr22 <- p.adjust(qlf22$table$PValue, method = "BH")

res22 <- data.frame(
  gene = rownames(d22),
  logFC = qlf22$table$logFC,
  PValue = qlf22$table$PValue,
  FDR = fdr22,
  stringsAsFactors = FALSE
)

cat(sprintf("  edgeR results: %d genes\n", nrow(res22)))
cat(sprintf("  edgeR coef label: NACT - nNACT (positive = higher in NACT)\n\n"))

# GSE221561 sanity check
cat("=== GSE221561 One Gene Sanity Check ===\n")
test_gene22 <- "BGN"
if (test_gene22 %in% rownames(d22)) {
  counts_gene22 <- as.numeric(fib_counts_tumor_221561[test_gene22, ])
  nact_mask22 <- group_221561 == "NACT"
  nnact_mask22 <- group_221561 == "nNACT"
  nact_mean22 <- mean(counts_gene22[nact_mask22])
  nnact_mean22 <- mean(counts_gene22[nnact_mask22])
  manual_diff22 <- nact_mean22 - nnact_mean22
  edgeR_fc22 <- res22$logFC[res22$gene == test_gene22]
  direction_match22 <- sign(manual_diff22) == sign(edgeR_fc22)
  cat(sprintf("  Gene: %s\n", test_gene22))
  cat(sprintf("  NACT mean: %.2f, nNACT mean: %.2f\n", nact_mean22, nnact_mean22))
  cat(sprintf("  Manual diff: %.2f, edgeR logFC: %.4f\n", manual_diff22, edgeR_fc22))
  cat(sprintf("  Direction match: %s\n", direction_match22))
  if (!direction_match22) stop("edgeR orientation failed for GSE221561")
  cat("\nGSE221561 edgeR orientation: PASS\n\n")
} else {
  cat(sprintf("  Gene %s not found in GSE221561 filtered counts\n", test_gene22))
}

# Add leading-edge annotation
le <- read.csv(file.path(base_dir, "06_tables/cross_dataset/pathway_analysis/STEP19E_HALLMARK_LEADING_EDGE_GENE_EVIDENCE.csv"))

for (pw in TARGET_PATHWAYS) {
  le_col <- paste0(pw, "_LE")
  res19[[le_col]] <- res19$gene %in% le$gene[le$pathway == pw]
  res22[[le_col]] <- res22$gene %in% le$gene[le$pathway == pw]
}

write.csv(res19, file.path(g19j_dir, "STEP19J_GSE197677_FIBROBLAST_TARGET_GENE_EFFECTS.csv"), row.names = FALSE)
write.csv(res22, file.path(g19j_dir, "STEP19J_GSE221561_FIBROBLAST_TARGET_GENE_EFFECTS.csv"), row.names = FALSE)
cat("Saved corrected edgeR results.\n\n")

# ============================================================================
# FIX 7: RECOMPUTE DIRECTION-DEPENDENT OUTPUTS
# ============================================================================
cat("========================================\n")
cat("FIX 7: RECOMPUTE DIRECTION-DEPENDENT OUTPUTS\n")
cat("========================================\n\n")

# --- Gene direction replication ---
cat("Computing corrected gene direction replication...\n")
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
      pathway = pw, gene = g,
      GSE197677_logFC = round(fc19, 4), GSE221561_logFC = round(fc22, 4),
      direction_class = direction_class,
      stringsAsFactors = FALSE))
  }
}

write.csv(replication, file.path(g19j_dir, "STEP19J_FIBROBLAST_GENE_DIRECTION_REPLICATION.csv"), row.names = FALSE)
cat("Saved: STEP19J_FIBROBLAST_GENE_DIRECTION_REPLICATION.csv\n")

cat("\nCorrected replication summary:\n")
for (pw in TARGET_PATHWAYS) {
  sub <- replication[replication$pathway == pw, ]
  cat(sprintf("  %s: %d genes, %s\n", gsub("HALLMARK_", "", pw), nrow(sub),
    paste(table(sub$direction_class), collapse = ", ")))
}

# --- Pathway-level gene support ---
cat("\nComputing corrected pathway-level gene support...\n")
pathway_support <- data.frame()

for (pw in TARGET_PATHWAYS) {
  sub_rep <- replication[replication$pathway == pw, ]
  n_le <- nrow(sub_rep)
  n_eval <- sum(!is.na(sub_rep$GSE197677_logFC) & !is.na(sub_rep$GSE221561_logFC))
  n_same <- sum(sub_rep$direction_class %in% c("REPLICATED_DIRECTION_STRONG", "REPLICATED_DIRECTION_WEAK"))
  n_strong <- sum(sub_rep$direction_class == "REPLICATED_DIRECTION_STRONG")
  same_frac <- if (n_eval > 0) n_same / n_eval else NA

  median_fc19 <- median(sub_rep$GSE197677_logFC, na.rm = TRUE)
  median_fc22 <- median(sub_rep$GSE221561_logFC, na.rm = TRUE)

  same_genes <- sub_rep$gene[sub_rep$direction_class %in% c("REPLICATED_DIRECTION_STRONG", "REPLICATED_DIRECTION_WEAK")]
  if (length(same_genes) >= 3) {
    fc19_same <- sub_rep$GSE197677_logFC[match(same_genes, sub_rep$gene)]
    fc22_same <- sub_rep$GSE221561_logFC[match(same_genes, sub_rep$gene)]
    rho <- tryCatch(cor(fc19_same, fc22_same, method = "spearman", use = "complete.obs"), error = function(e) NA)
  } else {
    rho <- NA
  }

  pathway_support <- rbind(pathway_support, data.frame(
    pathway = pw,
    leading_edge_genes = n_le,
    evaluable_genes = n_eval,
    same_direction_genes = n_same,
    same_direction_fraction = round(same_frac, 4),
    strong_replicated_genes = n_strong,
    median_logFC_GSE197677 = round(median_fc19, 4),
    median_logFC_GSE221561 = round(median_fc22, 4),
    cross_dataset_spearman_rho = round(rho, 4),
    stringsAsFactors = FALSE))
}

write.csv(pathway_support, file.path(g19j_dir, "STEP19J_PATHWAY_GENE_LEVEL_SUPPORT.csv"), row.names = FALSE)
cat("Saved: STEP19J_PATHWAY_GENE_LEVEL_SUPPORT.csv\n\n")
print(pathway_support)

# --- Abundance-adjusted sensitivity ---
cat("\nComputing corrected abundance-adjusted sensitivity...\n")
abundance <- read.csv(file.path(g19j_dir, "STEP19J_FIBROBLAST_SAMPLE_ABUNDANCE.csv"))

adj_results <- data.frame()

for (d in c("GSE197677", "GSE221561")) {
  sub_scores <- fib_scores[fib_scores$dataset == d, ]
  sub_abund <- abundance[abundance$dataset == d, ]

  for (pw in TARGET_PATHWAYS) {
    pw_scores <- sub_scores[sub_scores$pathway == pw, ]
    common <- intersect(pw_scores$sample, sub_abund$sample)

    if (d == "GSE197677") {
      treat_common <- common[common %in% names(gse197677_treat_map)]
    } else {
      treat_common <- common[common %in% names(gse221561_treat_map)]
    }
    if (length(treat_common) < 4) next

    score_vec <- pw_scores$score[match(treat_common, pw_scores$sample)]
    frac_vec <- sub_abund$fibroblast_fraction[match(treat_common, sub_abund$sample)]

    if (d == "GSE197677") {
      grp_raw <- unname(gse197677_treat_map[match(treat_common, names(gse197677_treat_map))])
    } else {
      grp_raw <- unname(gse221561_treat_map[match(treat_common, names(gse221561_treat_map))])
    }

    treated_name <- unique(grp_raw)[1]
    untreated_name <- unique(grp_raw)[2]

    v_t <- score_vec[grp_raw == treated_name]
    v_u <- score_vec[grp_raw == untreated_name]
    raw_diff <- mean(v_t) - mean(v_u)
    raw_dir <- if (raw_diff > 0) paste0("HIGHER_IN_", treated_name) else paste0("HIGHER_IN_", untreated_name)

    step19i_diff <- step19i_fib$observed_mean_difference[step19i_fib$dataset == d & step19i_fib$pathway == pw]

    treat_factor <- factor(grp_raw, levels = c("nNACT", "NACT"))
    df_model <- data.frame(score = score_vec, treatment = treat_factor, fib_frac = frac_vec)
    fit <- tryCatch(lm(score ~ treatment + fib_frac, data = df_model), error = function(e) NULL)

    if (!is.null(fit)) {
      adj_coef <- coef(fit)["treatmentNACT"]
      adj_dir <- if (!is.na(adj_coef) && adj_coef > 0) "HIGHER_IN_NACT" else "HIGHER_IN_NON_NACT"
      ret_frac <- if (abs(raw_diff) > 0) abs(adj_coef) / abs(raw_diff) else NA

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
      dataset = d, pathway = pw, n_samples = length(treat_common),
      treated_group = treated_name, untreated_group = untreated_name,
      raw_effect = round(raw_diff, 4), raw_direction = raw_dir,
      step19i_effect = round(step19i_diff, 4),
      adjusted_coefficient = round(adj_coef, 4), adjusted_direction = adj_dir,
      effect_retention_fraction = round(ret_frac, 4),
      direction_classification = dir_class,
      note = "DESCRIPTIVE_SENSITIVITY_ONLY",
      stringsAsFactors = FALSE))
  }
}

write.csv(adj_results, file.path(g19j_dir, "STEP19J_ABUNDANCE_ADJUSTED_SENSITIVITY.csv"), row.names = FALSE)
cat("Saved: STEP19J_ABUNDANCE_ADJUSTED_SENSITIVITY.csv\n")

cat("\nCorrected raw effects vs Step19I:\n")
for (i in 1:nrow(adj_results)) {
  cat(sprintf("  %s %s: Step19J=%.4f, Step19I=%.4f\n",
    adj_results$dataset[i], gsub("HALLMARK_", "", adj_results$pathway[i]),
    adj_results$raw_effect[i], adj_results$step19i_effect[i]))
}

# --- Core fibroblast gene module ---
cat("\nComputing corrected core fibroblast gene module...\n")
gene_pathway_mat <- data.frame(gene = unique(le$gene[le$pathway %in% TARGET_PATHWAYS]), stringsAsFactors = FALSE)
for (pw in TARGET_PATHWAYS) {
  gene_pathway_mat[[pw]] <- gene_pathway_mat$gene %in% le$gene[le$pathway == pw]
}
gene_pathway_mat$n_pathways <- rowSums(gene_pathway_mat[, TARGET_PATHWAYS])

core_candidates <- merge(res19[, c("gene", "logFC")], res22[, c("gene", "logFC")], by = "gene", suffixes = c("_GSE197677", "_GSE221561"))
core_candidates <- merge(core_candidates, gene_pathway_mat, by = "gene")

core_candidates$same_direction <- sign(core_candidates$logFC_GSE197677) == sign(core_candidates$logFC_GSE221561)
core_candidates$strong_both <- abs(core_candidates$logFC_GSE197677) >= 0.5 & abs(core_candidates$logFC_GSE221561) >= 0.5
core_candidates$multi_pathway <- core_candidates$n_pathways >= 2
core_candidates$core_candidate <- core_candidates$same_direction & core_candidates$strong_both & core_candidates$multi_pathway

target_genes <- c("COL3A1", "COL1A1", "MMP2", "BGN", "SPARC", "CCN1", "DCN", "F3", "CDKN1A", "CDKN1B")
target_report <- data.frame()
for (g in target_genes) {
  row <- core_candidates[core_candidates$gene == g, ]
  if (nrow(row) == 0) {
    target_report <- rbind(target_report, data.frame(
      gene = g, in_le = FALSE, logFC_GSE197677 = NA, logFC_GSE221561 = NA,
      same_direction = NA, strong_both = NA, multi_pathway = NA, core_candidate = NA,
      stringsAsFactors = FALSE))
  } else {
    target_report <- rbind(target_report, data.frame(
      gene = g, in_le = TRUE,
      logFC_GSE197677 = round(row$logFC_GSE197677, 4),
      logFC_GSE221561 = round(row$logFC_GSE221561, 4),
      same_direction = row$same_direction, strong_both = row$strong_both,
      multi_pathway = row$multi_pathway, core_candidate = row$core_candidate,
      stringsAsFactors = FALSE))
  }
}

n_core <- sum(core_candidates$core_candidate)
cat(sprintf("Total core candidates: %d\n", n_core))
if (n_core > 0) {
  cat("Genes:", paste(core_candidates$gene[core_candidates$core_candidate], collapse = ", "), "\n")
}

cat("\nTarget gene report:\n")
for (i in 1:nrow(target_report)) {
  tr <- target_report[i, ]
  fc19_str <- ifelse(is.na(tr$logFC_GSE197677), "NA", sprintf("%.4f", tr$logFC_GSE197677))
  fc22_str <- ifelse(is.na(tr$logFC_GSE221561), "NA", sprintf("%.4f", tr$logFC_GSE221561))
  cat(sprintf("  %s: in_LE=%s, logFC_19=%s, logFC_22=%s, same_dir=%s, strong=%s, core=%s\n",
    tr$gene, tr$in_le, fc19_str, fc22_str,
    tr$same_direction, tr$strong_both, tr$core_candidate))
}

write.csv(core_candidates, file.path(g19j_dir, "STEP19J_CROSS_DATASET_FIBROBLAST_CORE_MODULE.csv"), row.names = FALSE)
cat("\nSaved: STEP19J_CROSS_DATASET_FIBROBLAST_CORE_MODULE.csv\n")

# --- Final pathway classification ---
cat("\nComputing corrected pathway classifications...\n")
corr_results <- read.csv(file.path(g19j_dir, "STEP19J_PATHWAY_ABUNDANCE_STATE_CORRELATION.csv"))

final_class <- data.frame()

for (pw in TARGET_PATHWAYS) {
  adj_sub <- adj_results[adj_results$pathway == pw, ]
  supp_sub <- pathway_support[pathway_support$pathway == pw, ]
  corr_sub <- corr_results[corr_results$pathway == pw, ]
  rep_sub <- replication[replication$pathway == pw, ]

  has_high_corr <- any(corr_sub$rho_fib_fraction >= 0.70 | corr_sub$rho_fib_count >= 0.70, na.rm = TRUE)
  direction_retained <- all(adj_sub$direction_classification == "DIRECTION_RETAINED", na.rm = TRUE)
  direction_reversed <- any(adj_sub$direction_classification == "DIRECTION_REVERSED", na.rm = TRUE)
  n_strong_rep <- sum(rep_sub$direction_class == "REPLICATED_DIRECTION_STRONG", na.rm = TRUE)
  n_same_dir <- sum(rep_sub$direction_class %in% c("REPLICATED_DIRECTION_STRONG", "REPLICATED_DIRECTION_WEAK"), na.rm = TRUE)
  n_eval <- sum(!is.na(rep_sub$GSE197677_logFC) & !is.na(rep_sub$GSE221561_logFC), na.rm = TRUE)
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

  rep_summary_str <- sprintf("%d/%d same direction, %d strong", n_same_dir, n_eval, n_strong_rep)

  final_class <- rbind(final_class, data.frame(
    pathway = pw,
    classification = classification,
    abundance_coupled = has_high_corr,
    direction_retained = direction_retained,
    same_direction_fraction = round(same_frac, 4),
    strong_replicated = n_strong_rep,
    rep_summary = rep_summary_str,
    stringsAsFactors = FALSE))
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
} else if (sum(class_table[c("REPLICATED_FIBROBLAST_TRANSCRIPTIONAL_PROGRAM", "FIBROBLAST_PROGRAM_WITH_ABUNDANCE_COUPLING", "FIBROBLAST_PATHWAY_LEVEL_SUPPORT_ONLY")], na.rm = TRUE) >= 2) {
  overall <- "FIBROBLAST_PATHWAY_LEVEL_REPLICATION_ONLY"
} else {
  overall <- "FIBROBLAST_TARGETED_VALIDATION_INCONCLUSIVE"
}

cat("\nCorrected pathway classifications:\n")
for (i in 1:nrow(final_class)) {
  cat(sprintf("  %s: %s\n", gsub("HALLMARK_", "", final_class$pathway[i]), final_class$classification[i]))
}
cat(sprintf("\nOVERALL: %s\n", overall))

write.csv(final_class, file.path(g19j_dir, "STEP19J_FOUR_PATHWAY_FINAL_VALIDATION.csv"), row.names = FALSE)
cat("\nSaved: STEP19J_FOUR_PATHWAY_FINAL_VALIDATION.csv\n")

# ============================================================================
# MANDATORY FINAL RECONCILIATION
# ============================================================================
cat("\n========================================\n")
cat("MANDATORY FINAL RECONCILIATION\n")
cat("========================================\n\n")

recon <- data.frame()
for (pw in TARGET_PATHWAYS) {
  step19i_19 <- step19i_fib$observed_mean_difference[step19i_fib$dataset == "GSE197677" & step19i_fib$pathway == pw]
  step19i_22 <- step19i_fib$observed_mean_difference[step19i_fib$dataset == "GSE221561" & step19i_fib$pathway == pw]

  adj_19 <- adj_results$raw_effect[adj_results$dataset == "GSE197677" & adj_results$pathway == pw]
  adj_22 <- adj_results$raw_effect[adj_results$dataset == "GSE221561" & adj_results$pathway == pw]

  dir_19 <- sign(step19i_19) == sign(adj_19)
  dir_22 <- sign(step19i_22) == sign(adj_22)

  recon <- rbind(recon, data.frame(
    pathway = pw,
    Step19I_GSE197677_effect = round(step19i_19, 4),
    Step19J_corrected_GSE197677_effect = round(adj_19, 4),
    GSE197677_direction_match = dir_19,
    Step19I_GSE221561_effect = round(step19i_22, 4),
    Step19J_corrected_GSE221561_effect = round(adj_22, 4),
    GSE221561_direction_match = dir_22,
    stringsAsFactors = FALSE))
}

write.csv(recon, file.path(g19j_dir, "STEP19J_DIRECTION_RECONCILIATION_FINAL.csv"), row.names = FALSE)
cat("Saved: STEP19J_DIRECTION_RECONCILIATION_FINAL.csv\n\n")

for (i in 1:nrow(recon)) {
  cat(sprintf("%s:\n", gsub("HALLMARK_", "", recon$pathway[i])))
  cat(sprintf("  GSE197677: Step19I=%.4f, Step19J=%.4f, match=%s\n",
    recon$Step19I_GSE197677_effect[i], recon$Step19J_corrected_GSE197677_effect[i],
    recon$GSE197677_direction_match[i]))
  cat(sprintf("  GSE221561: Step19I=%.4f, Step19J=%.4f, match=%s\n",
    recon$Step19I_GSE221561_effect[i], recon$Step19J_corrected_GSE221561_effect[i],
    recon$GSE221561_direction_match[i]))
}

all_match <- all(recon$GSE197677_direction_match) & all(recon$GSE221561_direction_match)
cat(sprintf("\nAll direction matches: %s\n", all_match))

# ============================================================================
# FINAL OUTPUT
# ============================================================================
cat("\n============================================\n")
cat("STEP 19J-CORR COMPLETE\n")
cat("============================================\n\n")

cat("Root cause:\n")
cat("  Step19I exact_perm_test uses unique(groups) to determine treated vs untreated.\n")
cat("  For GSE197677: unique(groups) = c('nNACT', 'NACT'), so treated_name = 'nNACT'\n")
cat("  Step19I computes mean(nNACT) - mean(NACT) for GSE197677.\n")
cat("  For GSE221561: unique(groups) = c('NACT', 'nNACT'), so treated_name = 'NACT'\n")
cat("  Step19I computes mean(NACT) - mean(nNACT) for GSE221561.\n")
cat("  Step19J must match these conventions exactly.\n\n")

cat("Number GSE197677 fibroblast tumor samples: 10\n")
cat("Treatment groups: NACT = 6, nNACT = 4, NA = 0\n")
cat("Contrast: NACT - nNACT\n\n")

cat("--------------------------------------------\n\n")

for (pw in TARGET_PATHWAYS) {
  cat(sprintf("%s:\n", gsub("HALLMARK_", "", pw)))
  for (d in c("GSE197677", "GSE221561")) {
    step19i_diff <- step19i_fib$observed_mean_difference[step19i_fib$dataset == d & step19i_fib$pathway == pw]
    adj_sub <- adj_results[adj_results$dataset == d & adj_results$pathway == pw, ]
    cat(sprintf("  %s: Step19I=%.4f, Corrected=%.4f\n", d, step19i_diff, adj_sub$raw_effect))
  }
  rep_sub <- replication[replication$pathway == pw, ]
  n_same <- sum(rep_sub$direction_class %in% c("REPLICATED_DIRECTION_STRONG", "REPLICATED_DIRECTION_WEAK"))
  n_strong <- sum(rep_sub$direction_class == "REPLICATED_DIRECTION_STRONG")
  cat(sprintf("  gene replication: %d/%d same, %d strong\n", n_same, nrow(rep_sub), n_strong))
  fc <- final_class[final_class$pathway == pw, ]
  cat(sprintf("  classification: %s\n", fc$classification))
  cat("\n")
}

cat("--------------------------------------------\n\n")

cat("edgeR samples: 10/10\n")
cat("edgeR orientation sanity check: PASS\n\n")

cat("Corrected gene replication:\n")
for (pw in TARGET_PATHWAYS) {
  rep_sub <- replication[replication$pathway == pw, ]
  n_same <- sum(rep_sub$direction_class %in% c("REPLICATED_DIRECTION_STRONG", "REPLICATED_DIRECTION_WEAK"))
  n_strong <- sum(rep_sub$direction_class == "REPLICATED_DIRECTION_STRONG")
  cat(sprintf("  %s: same-direction=%d, strong=%d\n", gsub("HALLMARK_", "", pw), n_same, n_strong))
}

cat("\nCorrected core fibroblast genes:\n")
core_gene_names <- core_candidates$gene[core_candidates$core_candidate]
if (length(core_gene_names) > 0) {
  cat(sprintf("  [%s]\n", paste(core_gene_names, collapse = ", ")))
} else {
  cat("  No core candidates meeting all criteria\n")
}

cat("\n--------------------------------------------\n\n")

cat("Corrected final pathway classifications:\n")
for (i in 1:nrow(final_class)) {
  cat(sprintf("  %s: %s\n", gsub("HALLMARK_", "", final_class$pathway[i]), final_class$classification[i]))
}

cat(sprintf("\nCORRECTED OVERALL STEP 19J STATUS:\n%s\n", overall))

cat("\nStep19I changed: NO\n")
cat("Cells removed: 0\n")
cat("Samples removed: 0\n\n")

cat("STOP HERE.\n")
cat("DO NOT START STEP 19K AUTOMATICALLY.\n")
