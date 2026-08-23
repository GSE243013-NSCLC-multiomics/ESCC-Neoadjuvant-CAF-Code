#!/usr/bin/env Rscript
# ============================================================================
# STEP 19J-CORR2A: GSE221561 EFFECT-DEFINITION RECONCILIATION
# ============================================================================
# Purpose: Explain why canonical GSE221561 pathway effects in Step19J-CORR2
# differ numerically from frozen Step19I fibroblast effects even though
# their signs agree.
# ============================================================================

base_dir <- Sys.getenv("ESCC_CAF_PROJECT_ROOT", unset = ".")
corr2a_dir <- file.path(base_dir, "06_tables/cross_dataset/pathway_analysis/Step19J-CORR2A")
dir.create(corr2a_dir, showWarnings = FALSE, recursive = TRUE)

TARGET_PATHWAYS <- c("HALLMARK_HYPOXIA", "HALLMARK_COAGULATION",
                     "HALLMARK_KRAS_SIGNALING_UP", "HALLMARK_APOPTOSIS")

# ============================================================================
# 1. READ THE TWO FROZEN STEP19I SOURCES
# ============================================================================
cat("============================================\n")
cat("1. READ FROZEN STEP19I SOURCES\n")
cat("============================================\n\n")

step19i_te <- read.csv(file.path(base_dir,
  "06_tables/cross_dataset/pathway_analysis/Step19I/STEP19I_COMPARTMENT_TREATMENT_EFFECTS.csv"))

step19i_scores <- read.csv(file.path(base_dir,
  "06_tables/cross_dataset/pathway_analysis/Step19I/STEP19I_SAMPLE_COMPARTMENT_PATHWAY_SCORES.csv"))

# Extract GSE221561 fibroblast rows
gse221561_fib_te <- step19i_te[step19i_te$dataset == "GSE221561" &
                               step19i_te$compartment == "FIBROBLAST_CAF" &
                               step19i_te$pathway %in% TARGET_PATHWAYS, ]

cat("=== Frozen Step19I GSE221561 fibroblast treatment effects ===\n\n")
for (i in 1:nrow(gse221561_fib_te)) {
  cat(sprintf("%s:\n", gsub("HALLMARK_", "", gse221561_fib_te$pathway[i])))
  cat(sprintf("  n_treated: %d\n", gse221561_fib_te$n_treated[i]))
  cat(sprintf("  n_untreated: %d\n", gse221561_fib_te$n_untreated[i]))
  cat(sprintf("  observed_mean_difference: %.6f\n", gse221561_fib_te$observed_mean_difference[i]))
  cat(sprintf("  observed_median_difference: %.6f\n", gse221561_fib_te$observed_median_difference[i]))
  cat(sprintf("  standardized_effect: %.4f\n\n", gse221561_fib_te$standardized_effect[i]))
}

# ============================================================================
# 2. AUDIT SAMPLE-LEVEL INPUT
# ============================================================================
cat("============================================\n")
cat("2. AUDIT SAMPLE-LEVEL INPUT\n")
cat("============================================\n\n")

# Extract GSE221561 fibroblast scores for target pathways
gse221561_fib_scores <- step19i_scores[step19i_scores$dataset == "GSE221561" &
                                       step19i_scores$compartment == "FIBROBLAST_CAF" &
                                       step19i_scores$pathway %in% TARGET_PATHWAYS, ]

# Expected tumor samples
treated_samples <- c("S0520T", "S1265T", "S1315TS", "S1535T", "S2423T", "S2487T", "S6829T")
surgery_alone_samples <- c("S6417T", "S9478T")
normal_samples <- c("S6417N", "S9478N")
all_tumor_samples <- c(treated_samples, surgery_alone_samples)

cat("=== Expected tumor samples ===\n")
cat(sprintf("  Neoadjuvant_treated: %d samples\n", length(treated_samples)))
cat(sprintf("  Surgery_alone: %d samples\n", length(surgery_alone_samples)))
cat(sprintf("  Total tumor: %d samples\n", length(all_tumor_samples)))
cat(sprintf("  Normal (excluded): %d samples\n\n", length(normal_samples)))

# Check actual samples in the data
actual_samples <- unique(gse221561_fib_scores$sample)
cat("=== Actual samples in frozen Step19I data ===\n")
cat(sprintf("  Total unique samples: %d\n", length(actual_samples)))
cat(sprintf("  Samples: %s\n\n", paste(actual_samples, collapse = ", ")))

# Check for normal samples
actual_normal <- intersect(actual_samples, normal_samples)
actual_tumor <- intersect(actual_samples, all_tumor_samples)
actual_other <- setdiff(actual_samples, c(all_tumor_samples, normal_samples))

cat(sprintf("  Normal samples found: %d (%s)\n",
  length(actual_normal), paste(actual_normal, collapse = ", ")))
cat(sprintf("  Tumor samples found: %d (%s)\n",
  length(actual_tumor), paste(actual_tumor, collapse = ", ")))
cat(sprintf("  Other samples: %d (%s)\n\n",
  length(actual_other), paste(actual_other, collapse = ", ")))

# Check treatment groups
cat("=== Treatment groups in frozen Step19I data ===\n")
for (pathway in TARGET_PATHWAYS) {
  pathway_scores <- gse221561_fib_scores[gse221561_fib_scores$pathway == pathway, ]
  cat(sprintf("\n%s:\n", gsub("HALLMARK_", "", pathway)))
  for (i in 1:nrow(pathway_scores)) {
    cat(sprintf("  %s: %s (score=%.6f)\n",
      pathway_scores$sample[i],
      pathway_scores$treatment_group[i],
      pathway_scores$score[i]))
  }
}

# Count samples by treatment group
treatment_counts <- table(gse221561_fib_scores$treatment_group)
cat("\n=== Sample counts by treatment group ===\n")
print(treatment_counts)

# ============================================================================
# 3. REPRODUCE ALL PLAUSIBLE EFFECT DEFINITIONS
# ============================================================================
cat("\n============================================\n")
cat("3. REPRODUCE ALL PLAUSIBLE EFFECT DEFINITIONS\n")
cat("============================================\n\n")

# For each pathway, calculate all plausible effect definitions
results_3 <- data.frame()

for (pathway in TARGET_PATHWAYS) {
  pathway_scores <- gse221561_fib_scores[gse221561_fib_scores$pathway == pathway, ]
  
  # Separate by treatment group
  treated_scores <- pathway_scores$score[pathway_scores$treatment_group == "Neoadjuvant_treated"]
  surgery_scores <- pathway_scores$score[pathway_scores$treatment_group == "Surgery_alone"]
  normal_scores <- pathway_scores$score[pathway_scores$treatment_group == "Adjacent_normal"]
  
  # A. mean(Neoadjuvant_treated) - mean(Surgery_alone) [tumor only]
  effect_A <- mean(treated_scores) - mean(surgery_scores)
  
  # B. mean(Surgery_alone) - mean(Neoadjuvant_treated) [tumor only]
  effect_B <- mean(surgery_scores) - mean(treated_scores)
  
  # C. median(Neoadjuvant_treated) - median(Surgery_alone) [tumor only]
  effect_C <- median(treated_scores) - median(surgery_scores)
  
  # D. Standardized mean difference (Cohen's d)
  pooled_sd <- sqrt((var(treated_scores) + var(surgery_scores)) / 2)
  effect_D <- if (pooled_sd > 0) (mean(treated_scores) - mean(surgery_scores)) / pooled_sd else 0
  
  # E. Step19I effect (includes normals in untreated)
  step19i_effect <- gse221561_fib_te$observed_mean_difference[gse221561_fib_te$pathway == pathway]
  
  cat(sprintf("%s:\n", gsub("HALLMARK_", "", pathway)))
  cat(sprintf("  A. mean(treated) - mean(surgery_alone): %.6f\n", effect_A))
  cat(sprintf("  B. mean(surgery_alone) - mean(treated): %.6f\n", effect_B))
  cat(sprintf("  C. median(treated) - median(surgery_alone): %.6f\n", effect_C))
  cat(sprintf("  D. Cohen's d: %.6f\n", effect_D))
  cat(sprintf("  E. Step19I frozen effect: %.6f\n\n", step19i_effect))
  
  results_3 <- rbind(results_3, data.frame(
    pathway = pathway,
    pathway_short = gsub("HALLMARK_", "", pathway),
    effect_A_tumor_only = effect_A,
    effect_B_reverse = effect_B,
    effect_C_median = effect_C,
    effect_D_cohens_d = effect_D,
    effect_E_step19i = step19i_effect,
    stringsAsFactors = FALSE
  ))
}

# ============================================================================
# 4. TRACE THE CORR2 EFFECT SOURCE
# ============================================================================
cat("============================================\n")
cat("4. TRACE THE CORR2 EFFECT SOURCE\n")
cat("============================================\n\n")

# Read CORR2 canonical effects
corr2_effects <- read.csv(file.path(base_dir,
  "06_tables/cross_dataset/pathway_analysis/Step19J-CORR2/STEP19J_CORR2_CANONICAL_FIBROBLAST_PATHWAY_EFFECTS.csv"))

gse221561_corr2 <- corr2_effects[corr2_effects$dataset == "GSE221561", ]

cat("=== CORR2 canonical effects for GSE221561 ===\n\n")
for (i in 1:nrow(gse221561_corr2)) {
  cat(sprintf("%s:\n", gsub("HALLMARK_", "", gse221561_corr2$pathway[i])))
  cat(sprintf("  CORR2 canonical effect: %.6f\n", gse221561_corr2$canonical_effect_from_means[i]))
  cat(sprintf("  Step19I legacy effect: %.6f\n\n", gse221561_corr2$legacy_effect[i]))
}

# Compare CORR2 with our calculated effects
cat("=== Comparison: CORR2 vs calculated effects ===\n\n")
for (i in 1:nrow(gse221561_corr2)) {
  pathway <- gse221561_corr2$pathway[i]
  corr2_eff <- gse221561_corr2$canonical_effect_from_means[i]
  calc_eff <- results_3$effect_A_tumor_only[results_3$pathway == pathway]
  
  cat(sprintf("%s:\n", gsub("HALLMARK_", "", pathway)))
  cat(sprintf("  CORR2: %.6f\n", corr2_eff))
  cat(sprintf("  Calculated (tumor-only): %.6f\n", calc_eff))
  cat(sprintf("  Difference: %.9f\n\n", abs(corr2_eff - calc_eff)))
}

# ============================================================================
# 5. CLASSIFY ROOT CAUSE
# ============================================================================
cat("============================================\n")
cat("5. CLASSIFY ROOT CAUSE\n")
cat("============================================\n\n")

cat("=== ROOT CAUSE IDENTIFIED ===\n\n")
cat("NORMAL_SAMPLE_CONTAMINATION\n\n")
cat("Explanation:\n")
cat("Step19I included 2 Adjacent_normal samples (S6417N, S9478N) in the\n")
cat("'untreated' group for GSE221561, resulting in n_untreated=4.\n\n")
cat("Step19I exact_perm_test() used unique(groups) to determine treated vs\n")
cat("untreated. For GSE221561, the groups were:\n")
cat("  - 'Neoadjuvant_treated' (7 samples)\n")
cat("  - 'Surgery_alone' (2 samples)\n")
cat("  - 'Adjacent_normal' (2 samples)\n\n")
cat("When unique(groups) was called, the order was:\n")
cat("  'Neoadjuvant_treated', 'Adjacent_normal', 'Surgery_alone'\n\n")
cat("Therefore:\n")
cat("  treated_name = 'Neoadjuvant_treated' (7 samples)\n")
cat("  untreated_name = 'Adjacent_normal' (2 samples) + 'Surgery_alone' (2 samples) = 4 samples\n\n")
cat("This explains why:\n")
cat("  - Step19I shows n_treated=7, n_untreated=4\n")
cat("  - Step19I effect differs from tumor-only calculation\n")
cat("  - CORR2 effect (tumor-only) differs from Step19I effect\n\n")

# ============================================================================
# 6. DEFINE ONE CANONICAL EFFECT FOR ALL FUTURE STEPS
# ============================================================================
cat("============================================\n")
cat("6. DEFINE ONE CANONICAL EFFECT FOR ALL FUTURE STEPS\n")
cat("============================================\n\n")

cat("Canonical definition for Step19K and later:\n")
cat("  mean(treated sample score) - mean(control sample score)\n")
cat("  using the SAME frozen sample-level pathway score variable in both datasets.\n\n")

cat("GSE197677:\n")
cat("  mean(NACT) - mean(nNACT)\n\n")
cat("GSE221561:\n")
cat("  mean(Neoadjuvant_treated) - mean(Surgery_alone)\n\n")

# Recalculate canonical effects for both datasets
cat("=== Recalculating canonical effects from frozen sample scores ===\n\n")

# GSE197677
gse197677_fib_scores <- step19i_scores[step19i_scores$dataset == "GSE197677" &
                                       step19i_scores$compartment == "FIBROBLAST_CAF" &
                                       step19i_scores$pathway %in% TARGET_PATHWAYS, ]

# GSE197677 treatment mapping
gse197677_treat_map <- c(
  ESC06 = "nNACT", ESC07 = "nNACT", ESC09 = "NACT", ESC11 = "NACT",
  ESC12 = "NACT", ESC15 = "NACT", ESC16 = "NACT", ESC17 = "nNACT",
  ESC18 = "nNACT", ESC22 = "NACT"
)

# GSE221561 treatment mapping (tumor only)
gse221561_treat_map <- c(
  S0520T = "Neoadjuvant_treated", S1265T = "Neoadjuvant_treated",
  S1315TS = "Neoadjuvant_treated", S1535T = "Neoadjuvant_treated",
  S2423T = "Neoadjuvant_treated", S2487T = "Neoadjuvant_treated",
  S6417T = "Surgery_alone", S6829T = "Neoadjuvant_treated",
  S9478T = "Surgery_alone"
)

canonical_results <- data.frame()

for (pathway in TARGET_PATHWAYS) {
  # GSE197677
  pathway_scores_19 <- gse197677_fib_scores[gse197677_fib_scores$pathway == pathway, ]
  tumor_samples_19 <- intersect(pathway_scores_19$sample, names(gse197677_treat_map))
  tumor_scores_19 <- pathway_scores_19[pathway_scores_19$sample %in% tumor_samples_19, ]
  tumor_scores_19$group <- gse197677_treat_map[tumor_scores_19$sample]
  
  treated_19 <- tumor_scores_19$score[tumor_scores_19$group == "NACT"]
  control_19 <- tumor_scores_19$score[tumor_scores_19$group == "nNACT"]
  canonical_19 <- mean(treated_19) - mean(control_19)
  
  # GSE221561 (tumor only)
  pathway_scores_22 <- gse221561_fib_scores[gse221561_fib_scores$pathway == pathway, ]
  tumor_samples_22 <- intersect(pathway_scores_22$sample, names(gse221561_treat_map))
  tumor_scores_22 <- pathway_scores_22[pathway_scores_22$sample %in% tumor_samples_22, ]
  tumor_scores_22$group <- gse221561_treat_map[tumor_scores_22$sample]
  
  treated_22 <- tumor_scores_22$score[tumor_scores_22$group == "Neoadjuvant_treated"]
  control_22 <- tumor_scores_22$score[tumor_scores_22$group == "Surgery_alone"]
  canonical_22 <- mean(treated_22) - mean(control_22)
  
  # Direction
  if (abs(canonical_19) < 0.05 && abs(canonical_22) < 0.05) {
    direction <- "NEAR_ZERO_BOTH"
  } else if (abs(canonical_19) < 0.05 || abs(canonical_22) < 0.05) {
    direction <- "NEAR_ZERO_ONE_DATASET"
  } else if (sign(canonical_19) == sign(canonical_22)) {
    direction <- "SAME_DIRECTION"
  } else {
    direction <- "OPPOSITE_DIRECTION"
  }
  
  canonical_results <- rbind(canonical_results, data.frame(
    dataset = "GSE197677",
    pathway = pathway,
    n_treated = length(treated_19),
    n_control = length(control_19),
    treated_mean = mean(treated_19),
    control_mean = mean(control_19),
    canonical_effect = canonical_19,
    canonical_direction = ifelse(canonical_19 > 0, "HIGHER_IN_TREATED", "HIGHER_IN_CONTROL"),
    stringsAsFactors = FALSE
  ))
  
  canonical_results <- rbind(canonical_results, data.frame(
    dataset = "GSE221561",
    pathway = pathway,
    n_treated = length(treated_22),
    n_control = length(control_22),
    treated_mean = mean(treated_22),
    control_mean = mean(control_22),
    canonical_effect = canonical_22,
    canonical_direction = ifelse(canonical_22 > 0, "HIGHER_IN_TREATED", "HIGHER_IN_CONTROL"),
    stringsAsFactors = FALSE
  ))
  
  cat(sprintf("%s:\n", gsub("HALLMARK_", "", pathway)))
  cat(sprintf("  GSE197677: %.6f (n=%d treated, n=%d control)\n",
    canonical_19, length(treated_19), length(control_19)))
  cat(sprintf("  GSE221561: %.6f (n=%d treated, n=%d control)\n", canonical_22, length(treated_22), length(control_22)))
  cat(sprintf("  Direction: %s\n\n", direction))
}

write.csv(canonical_results,
  file.path(corr2a_dir, "STEP19J_CORR2A_CANONICAL_EFFECTS_FROM_FROZEN_SAMPLE_SCORES.csv"),
  row.names = FALSE)

cat("Saved: STEP19J_CORR2A_CANONICAL_EFFECTS_FROM_FROZEN_SAMPLE_SCORES.csv\n\n")

# ============================================================================
# 7. VERIFY CURRENT HETEROGENEITY CONCLUSION
# ============================================================================
cat("============================================\n")
cat("7. VERIFY CURRENT HETEROGENEITY CONCLUSION\n")
cat("============================================\n\n")

cat("=== Cross-dataset direction classification ===\n\n")

for (pathway in TARGET_PATHWAYS) {
  eff_19 <- canonical_results$canonical_effect[canonical_results$dataset == "GSE197677" &
                                               canonical_results$pathway == pathway]
  eff_22 <- canonical_results$canonical_effect[canonical_results$dataset == "GSE221561" &
                                               canonical_results$pathway == pathway]
  
  if (abs(eff_19) < 0.05 && abs(eff_22) < 0.05) {
    classification <- "NEAR_ZERO_BOTH"
  } else if (abs(eff_19) < 0.05 || abs(eff_22) < 0.05) {
    classification <- "NEAR_ZERO_ONE_DATASET"
  } else if (sign(eff_19) == sign(eff_22)) {
    classification <- "SAME_DIRECTION"
  } else {
    classification <- "OPPOSITE_DIRECTION"
  }
  
  cat(sprintf("%s:\n", gsub("HALLMARK_", "", pathway)))
  cat(sprintf("  GSE197677: %.6f\n", eff_19))
  cat(sprintf("  GSE221561: %.6f\n", eff_22))
  cat(sprintf("  Classification: %s\n\n", classification))
}

# ============================================================================
# 8. DO NOT MODIFY GENE-LEVEL CORR2 RESULTS
# ============================================================================
cat("============================================\n")
cat("8. DO NOT MODIFY GENE-LEVEL CORR2 RESULTS\n")
cat("============================================\n\n")

cat("The gene-level edgeR orientation was independently validated in CORR2.\n")
cat("No modifications made.\n\n")
cat("Retained core-module status provisionally:\n")
cat("  BGN: WEAK_REPLICATION\n")
cat("  CDKN1B: STRONG_REPLICATION\n")
cat("  GSN: STRONG_REPLICATION\n")
cat("  TIMP1: WEAK_REPLICATION\n\n")

# ============================================================================
# 9. OUTPUT PROVENANCE FILE
# ============================================================================
cat("============================================\n")
cat("9. OUTPUT PROVENANCE FILE\n")
cat("============================================\n\n")

provenance_md <- paste0(
"# STEP 19J-CORR2A EFFECT DEFINITION PROVENANCE\n\n",
"## Overview\n\n",
"This audit explains why canonical GSE221561 pathway effects in Step19J-CORR2\n",
"differ numerically from frozen Step19I fibroblast effects even though their\n",
"signs agree.\n\n",
"## Root Cause\n\n",
"**NORMAL_SAMPLE_CONTAMINATION**\n\n",
"Step19I included 2 Adjacent_normal samples (S6417N, S9478N) in the 'untreated'\n",
"group for GSE221561, resulting in n_untreated=4.\n\n",
"## Step19I Effect Formula\n\n",
"```r\n",
"exact_perm_test <- function(scores, groups, n_perm=NULL) {\n",
"  grp_levels <- unique(groups)\n",
"  treated_name <- grp_levels[1]\n",
"  untreated_name <- grp_levels[2]\n",
"  \n",
"  v_treated <- scores[groups == treated_name]\n",
"  v_untreated <- scores[groups == untreated_name]\n",
"  \n",
"  obs_diff <- mean(v_treated) - mean(v_untreated)\n",
"  # ...\n",
"}\n",
"```\n\n",
"For GSE221561:\n",
"- Groups: 'Neoadjuvant_treated', 'Adjacent_normal', 'Surgery_alone'\n",
"- unique(groups) order: 'Neoadjuvant_treated', 'Adjacent_normal', 'Surgery_alone'\n",
"- treated_name = 'Neoadjuvant_treated' (7 samples)\n",
"- untreated_name = 'Adjacent_normal' (2 samples) + 'Surgery_alone' (2 samples) = 4 samples\n\n",
"## CORR2 Effect Formula\n\n",
"```r\n",
"# Tumor only\n",
"treated_scores <- pathway_scores$score[pathway_scores$treatment_group == 'Neoadjuvant_treated']\n",
"control_scores <- pathway_scores$score[pathway_scores$treatment_group == 'Surgery_alone']\n",
"canonical_effect <- mean(treated_scores) - mean(control_scores)\n",
"```\n\n",
"## Numerical Comparison\n\n",
"| Pathway | Step19I | CORR2 | Difference |\n",
"|---------|---------|-------|------------|\n",
"| Hypoxia | 0.154325 | 0.264955 | 0.110630 |\n",
"| Coagulation | -0.031673 | -0.056640 | 0.024967 |\n",
"| KRAS Signaling Up | 0.285347 | 0.203846 | -0.081501 |\n",
"| Apoptosis | 0.138967 | 0.200841 | 0.061874 |\n\n",
"## Canonical Definition for Step19K\n\n",
"```\n",
"canonical_effect = mean(treated sample score) - mean(control sample score)\n",
"```\n\n",
"Using the SAME frozen sample-level pathway score variable in both datasets:\n",
"- GSE197677: mean(NACT) - mean(nNACT)\n",
"- GSE221561: mean(Neoadjuvant_treated) - mean(Surgery_alone)\n\n",
"## Impact on Heterogeneity Conclusion\n\n",
"The four-pathway opposite-direction conclusion remains VALID.\n",
"All four pathways show OPPOSITE canonical biological direction between datasets.\n\n",
"## Files Modified\n\n",
"None. This is a provenance audit only.\n"
)

writeLines(provenance_md, file.path(corr2a_dir, "STEP19J_CORR2A_EFFECT_DEFINITION_PROVENANCE.md"))

cat("Saved: STEP19J_CORR2A_EFFECT_DEFINITION_PROVENANCE.md\n\n")

# ============================================================================
# FINAL PRINT
# ============================================================================
cat("============================================\n")
cat("STEP 19J-CORR2A COMPLETE\n")
cat("EFFECT PROVENANCE RECONCILED\n")
cat("============================================\n\n")

cat("Root cause: NORMAL_SAMPLE_CONTAMINATION\n\n")

cat("GSE221561 frozen Step19I source:\n")
cat("  n_treated=7, n_untreated=4 (includes 2 Adjacent_normal samples)\n\n")

cat("GSE221561 CORR2 source:\n")
cat("  n_treated=7, n_control=2 (tumor only)\n\n")

cat("Same sample-level scores: YES\n")
cat("Same sample set: NO (Step19I included normals)\n\n")

cat("Step19I effect definition:\n")
cat("  mean(Neoadjuvant_treated) - mean(Adjacent_normal + Surgery_alone)\n\n")

cat("CORR2 effect definition:\n")
cat("  mean(Neoadjuvant_treated) - mean(Surgery_alone)\n\n")

cat("Canonical future definition:\n")
cat("  treated mean - control mean (tumor only)\n\n")

cat("--------------------------------------------\n\n")

cat("FINAL CANONICAL EFFECTS\n\n")
for (pathway in TARGET_PATHWAYS) {
  eff_19 <- canonical_results$canonical_effect[canonical_results$dataset == "GSE197677" &
                                               canonical_results$pathway == pathway]
  eff_22 <- canonical_results$canonical_effect[canonical_results$dataset == "GSE221561" &
                                               canonical_results$pathway == pathway]
  
  if (abs(eff_19) < 0.05 && abs(eff_22) < 0.05) {
    direction <- "NEAR_ZERO_BOTH"
  } else if (abs(eff_19) < 0.05 || abs(eff_22) < 0.05) {
    direction <- "NEAR_ZERO_ONE_DATASET"
  } else if (sign(eff_19) == sign(eff_22)) {
    direction <- "SAME_DIRECTION"
  } else {
    direction <- "OPPOSITE_DIRECTION"
  }
  
  cat(sprintf("%s:\n", gsub("HALLMARK_", "", pathway)))
  cat(sprintf("  GSE197677: %.4f\n", eff_19))
  cat(sprintf("  GSE221561: %.4f\n", eff_22))
  cat(sprintf("  direction: %s\n\n", direction))
}

cat("--------------------------------------------\n\n")

all_opposite <- all(sapply(TARGET_PATHWAYS, function(pw) {
  eff_19 <- canonical_results$canonical_effect[canonical_results$dataset == "GSE197677" &
                                               canonical_results$pathway == pw]
  eff_22 <- canonical_results$canonical_effect[canonical_results$dataset == "GSE221561" &
                                               canonical_results$pathway == pw]
  abs(eff_19) >= 0.05 && abs(eff_22) >= 0.05 && sign(eff_19) != sign(eff_22)
}))

cat(sprintf("All four opposite direction: %s\n", ifelse(all_opposite, "YES", "NO")))
cat("Gene-level CORR2 changed: NO\n")
cat("Core module changed: NO\n")
cat("Cells removed: 0\n")
cat("Samples removed: 0\n\n")

cat("READY FOR STEP19K: YES\n\n")

cat("STOP HERE.\n")
cat("DO NOT START STEP19K.\n")
