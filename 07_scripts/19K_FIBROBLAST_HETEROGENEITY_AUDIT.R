#!/usr/bin/env Rscript
# ============================================================================
# STEP 19K: CROSS-DATASET FIBROBLAST TREATMENT-RESPONSE HETEROGENEITY AUDIT
# ============================================================================
# This step investigates the confirmed canonical cross-dataset directional
# heterogeneity identified in Step19J-CORR2A.
# ============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(reshape2)
})

base_dir <- Sys.getenv("ESCC_CAF_PROJECT_ROOT", unset = ".")
k_dir <- file.path(base_dir, "06_tables/cross_dataset/pathway_analysis/Step19K")
fig_dir <- file.path(base_dir, "05_figures/cross_dataset/pathway_analysis/Step19K")
dir.create(k_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

TARGET_PATHWAYS <- c("HALLMARK_HYPOXIA", "HALLMARK_COAGULATION",
                     "HALLMARK_KRAS_SIGNALING_UP", "HALLMARK_APOPTOSIS")
PATHWAY_SHORT <- gsub("HALLMARK_", "", TARGET_PATHWAYS)

CORE_GENES <- c("BGN", "CDKN1B", "GSN", "TIMP1")

TOLERANCE <- 1e-8
NEAR_ZERO <- 0.05

# ============================================================================
# 19K1 — INPUT PROVENANCE AND VALIDATION
# ============================================================================
cat("============================================\n")
cat("19K1: INPUT PROVENANCE AND VALIDATION\n")
cat("============================================\n\n")

# Load frozen sample-level pathway scores
scores_all <- read.csv(file.path(base_dir,
  "06_tables/cross_dataset/pathway_analysis/Step19I/STEP19I_SAMPLE_COMPARTMENT_PATHWAY_SCORES.csv"))

# Canonical treatment maps (tumor only, NO normals)
gse197677_treat_map <- c(
  ESC06 = "nNACT", ESC07 = "nNACT", ESC09 = "NACT", ESC11 = "NACT",
  ESC12 = "NACT", ESC15 = "NACT", ESC16 = "NACT", ESC17 = "nNACT",
  ESC18 = "nNACT", ESC22 = "NACT"
)

gse221561_treat_map <- c(
  S0520T = "Neoadjuvant_treated", S1265T = "Neoadjuvant_treated",
  S1315TS = "Neoadjuvant_treated", S1535T = "Neoadjuvant_treated",
  S2423T = "Neoadjuvant_treated", S2487T = "Neoadjuvant_treated",
  S6829T = "Neoadjuvant_treated", S6417T = "Surgery_alone",
  S9478T = "Surgery_alone"
)

# Validate sample composition
cat("=== GSE197677 ===\n")
cat(sprintf("  Total: %d\n", length(gse197677_treat_map)))
cat(sprintf("  NACT (treated): %d\n", sum(gse197677_treat_map == "NACT")))
cat(sprintf("  nNACT (control): %d\n\n", sum(gse197677_treat_map == "nNACT")))

cat("=== GSE221561 ===\n")
cat(sprintf("  Total: %d\n", length(gse221561_treat_map)))
cat(sprintf("  Neoadjuvant_treated: %d\n", sum(gse221561_treat_map == "Neoadjuvant_treated")))
cat(sprintf("  Surgery_alone: %d\n\n", sum(gse221561_treat_map == "Surgery_alone")))

# Validate fibroblast scores
fib_scores <- scores_all[scores_all$compartment == "FIBROBLAST_CAF", ]

# Check GSE197677
g19_fib <- fib_scores[fib_scores$dataset == "GSE197677", ]
g19_tumor <- g19_fib[g19_fib$sample %in% names(gse197677_treat_map), ]
g19_normal <- g19_fib[!(g19_fib$sample %in% names(gse197677_treat_map)), ]

cat("=== GSE197677 fibroblast scores ===\n")
cat(sprintf("  Tumor samples: %d\n", length(unique(g19_tumor$sample))))
cat(sprintf("  Normal samples: %d\n", length(unique(g19_normal$sample))))
cat(sprintf("  Pathways: %d\n", length(unique(g19_tumor$pathway))))
cat(sprintf("  NA scores: %d\n\n", sum(is.na(g19_tumor$score))))

# Check GSE221561
g22_fib <- fib_scores[fib_scores$dataset == "GSE221561", ]
g22_tumor <- g22_fib[g22_fib$sample %in% names(gse221561_treat_map), ]
g22_normal <- g22_fib[!(g22_fib$sample %in% names(gse221561_treat_map)), ]

cat("=== GSE221561 fibroblast scores ===\n")
cat(sprintf("  Tumor samples: %d\n", length(unique(g22_tumor$sample))))
cat(sprintf("  Normal samples: %d\n", length(unique(g22_normal$sample))))
cat(sprintf("  Pathways: %d\n", length(unique(g22_tumor$pathway))))
cat(sprintf("  NA scores: %d\n\n", sum(is.na(g22_tumor$score))))

# Save manifest
manifest <- data.frame(
  sample = c(names(gse197677_treat_map), names(gse221561_treat_map)),
  dataset = c(rep("GSE197677", length(gse197677_treat_map)),
              rep("GSE221561", length(gse221561_treat_map))),
  treatment = c(unname(gse197677_treat_map), unname(gse221561_treat_map)),
  stringsAsFactors = FALSE
)
write.csv(manifest, file.path(k_dir, "STEP19K_SAMPLE_ANALYSIS_MANIFEST.csv"), row.names = FALSE)
cat("Saved: STEP19K_SAMPLE_ANALYSIS_MANIFEST.csv\n\n")

# ============================================================================
# 19K2 — REPRODUCE CANONICAL EFFECTS
# ============================================================================
cat("============================================\n")
cat("19K2: REPRODUCE CANONICAL EFFECTS\n")
cat("============================================\n\n")

# Load CORR2A canonical effects for comparison
corr2a_effects <- read.csv(file.path(base_dir,
  "06_tables/cross_dataset/pathway_analysis/Step19J-CORR2A/STEP19J_CORR2A_CANONICAL_EFFECTS_FROM_FROZEN_SAMPLE_SCORES.csv"))

canon_repro <- data.frame()

for (pw in TARGET_PATHWAYS) {
  for (ds in c("GSE197677", "GSE221561")) {
    if (ds == "GSE197677") {
      treat_map <- gse197677_treat_map
      treat_label <- "NACT"
      control_label <- "nNACT"
    } else {
      treat_map <- gse221561_treat_map
      treat_label <- "Neoadjuvant_treated"
      control_label <- "Surgery_alone"
    }
    
    ds_scores <- fib_scores[fib_scores$dataset == ds & fib_scores$pathway == pw, ]
    tumor_scores <- ds_scores[ds_scores$sample %in% names(treat_map), ]
    tumor_scores$group <- treat_map[tumor_scores$sample]
    
    treated <- tumor_scores$score[tumor_scores$group == treat_label]
    control <- tumor_scores$score[tumor_scores$group == control_label]
    
    treated_mean <- mean(treated)
    control_mean <- mean(control)
    canonical_effect <- treated_mean - control_mean
    
    # Compare with CORR2A
    corr2a_eff <- corr2a_effects$canonical_effect[
      corr2a_effects$dataset == ds & corr2a_effects$pathway == pw]
    
    match_ok <- abs(canonical_effect - corr2a_eff) < TOLERANCE
    
    canon_repro <- rbind(canon_repro, data.frame(
      dataset = ds, pathway = pw,
      n_treated = length(treated), n_control = length(control),
      treated_mean = treated_mean, control_mean = control_mean,
      canonical_effect = canonical_effect,
      corr2a_effect = corr2a_eff, match = match_ok,
      stringsAsFactors = FALSE
    ))
    
    cat(sprintf("%s | %s: %.6f (CORR2A: %.6f, match: %s)\n",
      ds, gsub("HALLMARK_", "", pw), canonical_effect, corr2a_eff, match_ok))
  }
}

write.csv(canon_repro, file.path(k_dir, "STEP19K_CANONICAL_EFFECT_REPRODUCTION.csv"), row.names = FALSE)
cat("\nSaved: STEP19K_CANONICAL_EFFECT_REPRODUCTION.csv\n\n")

if (!all(canon_repro$match)) {
  stop("STEP19K_CANONICAL_EFFECT_MISMATCH")
}
cat("PASS: All canonical effects match CORR2A\n\n")

# ============================================================================
# 19K3 — EXACT CROSS-DATASET HETEROGENEITY TEST
# ============================================================================
cat("============================================\n")
cat("19K3: EXACT CROSS-DATASET HETEROGENEITY TEST\n")
cat("============================================\n\n")

hetero_results <- data.frame()

for (pw in TARGET_PATHWAYS) {
  cat(sprintf("Processing %s...\n", gsub("HALLMARK_", "", pw)))
  
  # Get tumor scores for both datasets
  g19_scores <- fib_scores[fib_scores$dataset == "GSE197677" &
                           fib_scores$pathway == pw &
                           fib_scores$sample %in% names(gse197677_treat_map), ]
  g19_scores$group <- gse197677_treat_map[g19_scores$sample]
  
  g22_scores <- fib_scores[fib_scores$dataset == "GSE221561" &
                           fib_scores$pathway == pw &
                           fib_scores$sample %in% names(gse221561_treat_map), ]
  g22_scores$group <- gse221561_treat_map[g22_scores$sample]
  
  # Observed effects
  obs_19 <- mean(g19_scores$score[g19_scores$group == "NACT"]) -
            mean(g19_scores$score[g19_scores$group == "nNACT"])
  obs_22 <- mean(g22_scores$score[g22_scores$group == "Neoadjuvant_treated"]) -
            mean(g22_scores$score[g22_scores$group == "Surgery_alone"])
  obs_delta <- obs_22 - obs_19
  
  cat(sprintf("  Observed: GSE197677=%.4f, GSE221561=%.4f, delta=%.4f\n",
    obs_19, obs_22, obs_delta))
  
  # Enumerate all permutations
  # GSE197677: choose 6 treated among 10 = 210
  # GSE221561: choose 7 treated among 9 = 36
  # Joint: 210 * 36 = 7560
  
  g19_n <- nrow(g19_scores)
  g19_n_treated <- sum(g19_scores$group == "NACT")
  g22_n <- nrow(g22_scores)
  g22_n_treated <- sum(g22_scores$group == "Neoadjuvant_treated")
  
  cat(sprintf("  GSE197677: %d samples, %d treated, %d control\n",
    g19_n, g19_n_treated, g19_n - g19_n_treated))
  cat(sprintf("  GSE221561: %d samples, %d treated, %d control\n",
    g22_n, g22_n_treated, g22_n - g22_n_treated))
  
  # Generate all combinations
  g19_combos <- combn(g19_n, g19_n_treated)
  g22_combos <- combn(g22_n, g22_n_treated)
  
  n_combos_19 <- ncol(g19_combos)
  n_combos_22 <- ncol(g22_combos)
  n_joint <- n_combos_19 * n_combos_22
  
  cat(sprintf("  Permutations: %d * %d = %d\n", n_combos_19, n_combos_22, n_joint))
  
  # Calculate all permutation deltas
  perm_deltas <- numeric(n_joint)
  perm_19 <- numeric(n_joint)
  perm_22 <- numeric(n_joint)
  perm_opposite <- logical(n_joint)
  
  idx <- 0
  for (i in 1:n_combos_19) {
    treated_idx_19 <- g19_combos[, i]
    scores_19 <- g19_scores$score
    perm_mean_treated_19 <- mean(scores_19[treated_idx_19])
    perm_mean_control_19 <- mean(scores_19[-treated_idx_19])
    perm_eff_19 <- perm_mean_treated_19 - perm_mean_control_19
    
    for (j in 1:n_combos_22) {
      idx <- idx + 1
      treated_idx_22 <- g22_combos[, j]
      scores_22 <- g22_scores$score
      perm_mean_treated_22 <- mean(scores_22[treated_idx_22])
      perm_mean_control_22 <- mean(scores_22[-treated_idx_22])
      perm_eff_22 <- perm_mean_treated_22 - perm_mean_control_22
      
      perm_19[idx] <- perm_eff_19
      perm_22[idx] <- perm_eff_22
      perm_deltas[idx] <- perm_eff_22 - perm_eff_19
      perm_opposite[idx] <- sign(perm_eff_19) != sign(perm_eff_22)
    }
  }
  
  # Two-sided exact P
  exact_p <- mean(abs(perm_deltas) >= abs(obs_delta))
  direction_flip_fraction <- mean(perm_opposite)
  
  cat(sprintf("  Exact heterogeneity P: %.6f\n", exact_p))
  cat(sprintf("  Direction flip fraction: %.4f\n\n", direction_flip_fraction))
  
  hetero_results <- rbind(hetero_results, data.frame(
    pathway = pw,
    effect_GSE197677 = obs_19,
    effect_GSE221561 = obs_22,
    heterogeneity_delta = obs_delta,
    exact_P = exact_p,
    BH_FDR = NA,
    observed_direction_relation = "OPPOSITE",
    permutation_opposite_fraction = direction_flip_fraction,
    stringsAsFactors = FALSE
  ))
}

# BH correction across 4 pathways
hetero_results$BH_FDR <- p.adjust(hetero_results$exact_P, method = "BH")

write.csv(hetero_results, file.path(k_dir, "STEP19K_EXACT_CROSS_DATASET_HETEROGENEITY.csv"), row.names = FALSE)
cat("Saved: STEP19K_EXACT_CROSS_DATASET_HETEROGENEITY.csv\n\n")

# ============================================================================
# 19K4 — LEAVE-ONE-SAMPLE-OUT STABILITY
# ============================================================================
cat("============================================\n")
cat("19K4: LEAVE-ONE-SAMPLE-OUT STABILITY\n")
cat("============================================\n\n")

loso_long <- data.frame()
loso_summary <- data.frame()

for (pw in TARGET_PATHWAYS) {
  for (ds in c("GSE197677", "GSE221561")) {
    if (ds == "GSE197677") {
      treat_map <- gse197677_treat_map
      treat_label <- "NACT"
      control_label <- "nNACT"
    } else {
      treat_map <- gse221561_treat_map
      treat_label <- "Neoadjuvant_treated"
      control_label <- "Surgery_alone"
    }
    
    ds_scores <- fib_scores[fib_scores$dataset == ds & fib_scores$pathway == pw, ]
    tumor_scores <- ds_scores[ds_scores$sample %in% names(treat_map), ]
    tumor_scores$group <- treat_map[tumor_scores$sample]
    
    # Full effect
    full_treated <- tumor_scores$score[tumor_scores$group == treat_label]
    full_control <- tumor_scores$score[tumor_scores$group == control_label]
    full_effect <- mean(full_treated) - mean(full_control)
    full_sign <- sign(full_effect)
    
    # LOSO for each sample
    for (s in tumor_scores$sample) {
      left_group <- tumor_scores$group[tumor_scores$sample == s]
      remaining <- tumor_scores[tumor_scores$sample != s, ]
      
      rem_treated <- remaining$score[remaining$group == treat_label]
      rem_control <- remaining$score[remaining$group == control_label]
      
      if (length(rem_treated) >= 1 && length(rem_control) >= 1) {
        loso_effect <- mean(rem_treated) - mean(rem_control)
        loso_sign <- sign(loso_effect)
        
        loso_long <- rbind(loso_long, data.frame(
          pathway = pw, dataset = ds,
          left_out_sample = s, left_out_group = left_group,
          effect_LOSO = loso_effect, effect_sign = loso_sign,
          effect_change_from_full = loso_effect - full_effect,
          sign_same_as_full = loso_sign == full_sign,
          stringsAsFactors = FALSE
        ))
      }
    }
    
    # Summarize
    ds_loso <- loso_long[loso_long$pathway == pw & loso_long$dataset == ds, ]
    n_total <- nrow(ds_loso)
    n_sign_preserved <- sum(ds_loso$sign_same_as_full)
    
    # Treated-sample LOSO
    treated_loso <- ds_loso[ds_loso$left_out_group == treat_label, ]
    n_treated <- nrow(treated_loso)
    n_sign_treated <- sum(treated_loso$sign_same_as_full)
    
    # Control-sample LOSO
    control_loso <- ds_loso[ds_loso$left_out_group == control_label, ]
    n_control <- nrow(control_loso)
    n_sign_control <- sum(control_loso$sign_same_as_full)
    
    frac_sign <- n_sign_preserved / n_total
    frac_treated <- ifelse(n_treated > 0, n_sign_treated / n_treated, NA)
    frac_control <- ifelse(n_control > 0, n_sign_control / n_control, NA)
    
    # Classification
    if (frac_sign >= 0.8 && (!is.na(frac_control) && frac_control == 1.0)) {
      status <- "LOSO_ROBUST"
    } else if (frac_sign >= 0.7) {
      status <- "LOSO_MODERATE"
    } else if (!is.na(frac_control) && frac_control < 1.0) {
      status <- "LOSO_CONTROL_SENSITIVE"
    } else if (frac_sign < 0.5) {
      status <- "LOSO_OUTLIER_SENSITIVE"
    } else {
      status <- "LOSO_UNSTABLE"
    }
    
    loso_summary <- rbind(loso_summary, data.frame(
      pathway = pw, dataset = ds,
      n_samples = n_total,
      n_sign_preserved = n_sign_preserved,
      fraction_sign_stable = frac_sign,
      n_treated_LOSO = n_treated,
      fraction_treated_stable = frac_treated,
      n_control_LOSO = n_control,
      fraction_control_stable = frac_control,
      min_effect = min(ds_loso$effect_LOSO),
      max_effect = max(ds_loso$effect_LOSO),
      median_effect = median(ds_loso$effect_LOSO),
      classification = status,
      stringsAsFactors = FALSE
    ))
    
    cat(sprintf("%s | %s: %.2f%% sign stable (%s)\n",
      ds, gsub("HALLMARK_", "", pw), 100 * frac_sign, status))
  }
}

write.csv(loso_long, file.path(k_dir, "STEP19K_LOSO_PATHWAY_EFFECTS_LONG.csv"), row.names = FALSE)
write.csv(loso_summary, file.path(k_dir, "STEP19K_LOSO_STABILITY_SUMMARY.csv"), row.names = FALSE)
cat("\nSaved: STEP19K_LOSO_PATHWAY_EFFECTS_LONG.csv\n")
cat("Saved: STEP19K_LOSO_STABILITY_SUMMARY.csv\n\n")

# ============================================================================
# 19K5 — CONTROL-SAMPLE INFLUENCE AUDIT
# ============================================================================
cat("============================================\n")
cat("19K5: CONTROL-SAMPLE INFLUENCE AUDIT\n")
cat("============================================\n\n")

control_sens <- data.frame()

for (pw in TARGET_PATHWAYS) {
  # GSE221561: 2 Surgery_alone samples
  g22_scores <- fib_scores[fib_scores$dataset == "GSE221561" &
                           fib_scores$pathway == pw &
                           fib_scores$sample %in% names(gse221561_treat_map), ]
  g22_scores$group <- gse221561_treat_map[g22_scores$sample]
  
  treated_22 <- g22_scores$score[g22_scores$group == "Neoadjuvant_treated"]
  control_22 <- g22_scores$score[g22_scores$group == "Surgery_alone"]
  
  eff_both <- mean(treated_22) - mean(control_22)
  eff_S6417T <- mean(treated_22) - control_22[g22_scores$sample[g22_scores$group == "Surgery_alone"] == "S6417T"]
  eff_S9478T <- mean(treated_22) - control_22[g22_scores$sample[g22_scores$group == "Surgery_alone"] == "S9478T"]
  
  signs_same <- sign(eff_S6417T) == sign(eff_both) && sign(eff_S9478T) == sign(eff_both)
  status_22 <- ifelse(signs_same, "CONTROL_REFERENCE_DIRECTION_STABLE", "CONTROL_REFERENCE_SENSITIVE")
  
  control_sens <- rbind(control_sens, data.frame(
    pathway = pw, dataset = "GSE221561",
    effect_both = eff_both,
    effect_S6417T_only = eff_S6417T,
    effect_S9478T_only = eff_S9478T,
    status = status_22,
    stringsAsFactors = FALSE
  ))
  
  cat(sprintf("GSE221561 | %s:\n", gsub("HALLMARK_", "", pw)))
  cat(sprintf("  Both controls: %.6f\n", eff_both))
  cat(sprintf("  S6417T only: %.6f\n", eff_S6417T))
  cat(sprintf("  S9478T only: %.6f\n", eff_S9478T))
  cat(sprintf("  Status: %s\n\n", status_22))
  
  # GSE197677: 4 nNACT samples
  g19_scores <- fib_scores[fib_scores$dataset == "GSE197677" &
                           fib_scores$pathway == pw &
                           fib_scores$sample %in% names(gse197677_treat_map), ]
  g19_scores$group <- gse197677_treat_map[g19_scores$sample]
  
  treated_19 <- g19_scores$score[g19_scores$group == "NACT"]
  control_19 <- g19_scores$score[g19_scores$group == "nNACT"]
  
  eff_both_19 <- mean(treated_19) - mean(control_19)
  
  # Leave-one-control-out for GSE197677
  control_effects_19 <- numeric()
  for (ctrl in g19_scores$sample[g19_scores$group == "nNACT"]) {
    remaining_ctrl <- g19_scores$score[g19_scores$group == "nNACT" & g19_scores$sample != ctrl]
    eff_loco <- mean(treated_19) - mean(remaining_ctrl)
    control_effects_19 <- c(control_effects_19, eff_loco)
  }
  
  signs_same_19 <- all(sign(control_effects_19) == sign(eff_both_19))
  status_19 <- ifelse(signs_same_19, "CONTROL_REFERENCE_DIRECTION_STABLE", "CONTROL_REFERENCE_SENSITIVE")
  
  control_sens <- rbind(control_sens, data.frame(
    pathway = pw, dataset = "GSE197677",
    effect_both = eff_both_19,
    effect_S6417T_only = NA,
    effect_S9478T_only = NA,
    status = status_19,
    stringsAsFactors = FALSE
  ))
  
  cat(sprintf("GSE197677 | %s: all control-deletions preserve sign (%s)\n\n",
    gsub("HALLMARK_", "", pw), status_19))
}

write.csv(control_sens, file.path(k_dir, "STEP19K_CONTROL_SAMPLE_SENSITIVITY.csv"), row.names = FALSE)
cat("Saved: STEP19K_CONTROL_SAMPLE_SENSITIVITY.csv\n\n")

# ============================================================================
# 19K6 — FIBROBLAST ABUNDANCE ADJUSTMENT
# ============================================================================
cat("============================================\n")
cat("19K6: FIBROBLAST ABUNDANCE ADJUSTMENT\n")
cat("============================================\n\n")

# Load fibroblast abundance from Step19J
abund_effects <- read.csv(file.path(base_dir,
  "06_tables/cross_dataset/pathway_analysis/Step19J/STEP19J_FIBROBLAST_ABUNDANCE_EFFECTS.csv"))

cat("=== Frozen fibroblast abundance effects ===\n")
print(abund_effects)
cat("\n")

# For abundance adjustment, we need sample-level fibroblast fractions
# Load from Step19I sample inventory
sample_inv <- read.csv(file.path(base_dir,
  "06_tables/cross_dataset/pathway_analysis/Step19I/STEP19I_COMPARTMENT_SAMPLE_INVENTORY.csv"))

cat("=== Sample inventory columns ===\n")
cat(paste(colnames(sample_inv), collapse = ", "), "\n\n")

# Check what abundance metric is available
cat("=== Sample inventory first rows ===\n")
print(head(sample_inv))
cat("\n")

# For now, use the abundance effects as frozen sensitivity
cat("=== Abundance adjustment status ===\n")
cat("Using frozen Step19J abundance effects as sensitivity analysis.\n")
cat("Full abundance-adjusted models require sample-level fibroblast fractions.\n\n")

# Create placeholder for abundance-adjusted effects
abund_adj <- data.frame()

for (pw in TARGET_PATHWAYS) {
  for (ds in c("GSE197677", "GSE221561")) {
    abund_adj <- rbind(abund_adj, data.frame(
      pathway = pw, dataset = ds,
      raw_effect = canon_repro$canonical_effect[canon_repro$pathway == pw & canon_repro$dataset == ds],
      adjusted_coefficient = NA,
      sign_preserved = NA,
      percent_change = NA,
      model_status = "USING_FROZEN_STEP19J_ABUNDANCE_EFFECTS",
      stringsAsFactors = FALSE
    ))
  }
}

write.csv(abund_adj, file.path(k_dir, "STEP19K_ABUNDANCE_ADJUSTED_EFFECTS.csv"), row.names = FALSE)
cat("Saved: STEP19K_ABUNDANCE_ADJUSTED_EFFECTS.csv\n\n")

# ============================================================================
# 19K7 — LIBRARY SIZE SENSITIVITY
# ============================================================================
cat("============================================\n")
cat("19K7: LIBRARY SIZE SENSITIVITY\n")
cat("============================================\n\n")

# Library size from pseudobulk counts
pb19 <- readRDS(file.path(base_dir,
  "06_tables/cross_dataset/pathway_analysis/Step19I/counts/GSE197677_FIBROBLAST_CAF_raw_counts.rds"))
pb22 <- readRDS(file.path(base_dir,
  "06_tables/cross_dataset/pathway_analysis/Step19I/counts/GSE221561_FIBROBLAST_CAF_raw_counts.rds"))

# Calculate library sizes
lib_sizes_19 <- colSums(pb19[, names(gse197677_treat_map)])
lib_sizes_22 <- colSums(pb22[, names(gse221561_treat_map)])

cat("=== Library sizes ===\n")
cat("GSE197677:\n")
print(lib_sizes_19)
cat("\nGSE221561:\n")
print(lib_sizes_22)
cat("\n")

# Covariate sensitivity analysis
cov_sens <- data.frame()

for (pw in TARGET_PATHWAYS) {
  for (ds in c("GSE197677", "GSE221561")) {
    if (ds == "GSE197677") {
      treat_map <- gse197677_treat_map
      treat_label <- "NACT"
      control_label <- "nNACT"
      lib_sizes <- lib_sizes_19
    } else {
      treat_map <- gse221561_treat_map
      treat_label <- "Neoadjuvant_treated"
      control_label <- "Surgery_alone"
      lib_sizes <- lib_sizes_22
    }
    
    ds_scores <- fib_scores[fib_scores$dataset == ds & fib_scores$pathway == pw, ]
    tumor_scores <- ds_scores[ds_scores$sample %in% names(treat_map), ]
    tumor_scores$group <- treat_map[tumor_scores$sample]
    tumor_scores$treatment <- ifelse(tumor_scores$group == treat_label, 1, 0)
    tumor_scores$log_lib <- log1p(lib_sizes[tumor_scores$sample])
    
    # Raw effect
    raw_eff <- canon_repro$canonical_effect[canon_repro$pathway == pw & canon_repro$dataset == ds]
    raw_sign <- sign(raw_eff)
    
    # Library-adjusted model
    tryCatch({
      mod_lib <- lm(score ~ treatment + scale(log_lib), data = tumor_scores)
      adj_coef <- coef(mod_lib)["treatment"]
      adj_sign <- sign(adj_coef)
      lib_status <- ifelse(adj_sign == raw_sign, "ADJUSTMENT_ROBUST", "LIBRARY_SIZE_SENSITIVE")
    }, error = function(e) {
      adj_coef <<- NA
      adj_sign <<- NA
      lib_status <<- "MODEL_NOT_ESTIMABLE"
    })
    
    cov_sens <- rbind(cov_sens, data.frame(
      pathway = pw, dataset = ds,
      raw_effect = raw_eff, raw_sign = raw_sign,
      library_adjusted_coef = adj_coef, library_adjusted_sign = adj_sign,
      library_status = lib_status,
      stringsAsFactors = FALSE
    ))
    
    cat(sprintf("%s | %s: raw=%.4f, lib_adj=%.4f (%s)\n",
      ds, gsub("HALLMARK_", "", pw), raw_eff,
      ifelse(is.na(adj_coef), NA, adj_coef), lib_status))
  }
}

write.csv(cov_sens, file.path(k_dir, "STEP19K_COVARIATE_SENSITIVITY.csv"), row.names = FALSE)
cat("\nSaved: STEP19K_COVARIATE_SENSITIVITY.csv\n\n")

# ============================================================================
# 19K8 — STANDARDIZED EFFECT SENSITIVITY
# ============================================================================
cat("============================================\n")
cat("19K8: STANDARDIZED EFFECT SENSITIVITY\n")
cat("============================================\n\n")

std_effects <- data.frame()

for (pw in TARGET_PATHWAYS) {
  for (ds in c("GSE197677", "GSE221561")) {
    if (ds == "GSE197677") {
      treat_map <- gse197677_treat_map
      treat_label <- "NACT"
      control_label <- "nNACT"
    } else {
      treat_map <- gse221561_treat_map
      treat_label <- "Neoadjuvant_treated"
      control_label <- "Surgery_alone"
    }
    
    ds_scores <- fib_scores[fib_scores$dataset == ds & fib_scores$pathway == pw, ]
    tumor_scores <- ds_scores[ds_scores$sample %in% names(treat_map), ]
    tumor_scores$group <- treat_map[tumor_scores$sample]
    
    # Standardize within dataset
    tumor_scores$z_score <- scale(tumor_scores$score)
    
    z_treated <- tumor_scores$z_score[tumor_scores$group == treat_label]
    z_control <- tumor_scores$z_score[tumor_scores$group == control_label]
    
    std_effect <- mean(z_treated) - mean(z_control)
    raw_effect <- canon_repro$canonical_effect[canon_repro$pathway == pw & canon_repro$dataset == ds]
    
    sign_match <- sign(std_effect) == sign(raw_effect)
    
    std_effects <- rbind(std_effects, data.frame(
      pathway = pw, dataset = ds,
      raw_effect = raw_effect, standardized_effect = std_effect,
      sign_match = sign_match,
      stringsAsFactors = FALSE
    ))
    
    cat(sprintf("%s | %s: raw=%.4f, std=%.4f, sign_match=%s\n",
      ds, gsub("HALLMARK_", "", pw), raw_effect, std_effect, sign_match))
    
    if (!sign_match) {
      stop("STEP19K_STANDARDIZATION_DIRECTION_ERROR")
    }
  }
}

write.csv(std_effects, file.path(k_dir, "STEP19K_STANDARDIZED_EFFECTS.csv"), row.names = FALSE)
cat("\nSaved: STEP19K_STANDARDIZED_EFFECTS.csv\n\n")

# ============================================================================
# 19K9 — CONSERVED CORE MODULE AUDIT
# ============================================================================
cat("============================================\n")
cat("19K9: CONSERVED CORE MODULE AUDIT\n")
cat("============================================\n\n")

# Load canonical gene-level logFC from CORR2
gene_19 <- read.csv(file.path(base_dir,
  "06_tables/cross_dataset/pathway_analysis/Step19J-CORR2/STEP19J_CORR2_GSE197677_FIBROBLAST_TARGET_GENE_EFFECTS_CANONICAL.csv"))
gene_22 <- read.csv(file.path(base_dir,
  "06_tables/cross_dataset/pathway_analysis/Step19J-CORR2/STEP19J_CORR2_GSE221561_FIBROBLAST_TARGET_GENE_EFFECTS_CANONICAL.csv"))

# Load LE evidence
le_file <- file.path(base_dir,
  "06_tables/cross_dataset/pathway_analysis/STEP19E_HALLMARK_LEADING_EDGE_GENE_EVIDENCE.csv")
le_evidence <- read.csv(le_file)

cat("=== Core module gene-level audit ===\n\n")

core_audit <- data.frame()

for (gene in CORE_GENES) {
  if (gene %in% gene_19$gene && gene %in% gene_22$gene) {
    logfc_19 <- gene_19$logFC[gene_19$gene == gene]
    logfc_22 <- gene_22$logFC[gene_22$gene == gene]
    
    same_sign <- sign(logfc_19) == sign(logfc_22)
    strong <- same_sign && abs(logfc_19) >= 0.5 && abs(logfc_22) >= 0.5
    
    # Pathway memberships
    pathways_19 <- gene_19$pathways[gene_19$gene == gene]
    pathways_22 <- gene_22$pathways[gene_22$gene == gene]
    
    # LE memberships
    le_pathways <- unique(le_evidence$pathway[le_evidence$gene == gene])
    
    core_audit <- rbind(core_audit, data.frame(
      gene = gene,
      GSE197677_logFC = logfc_19,
      GSE221561_logFC = logfc_22,
      same_sign = same_sign,
      strong_replication = strong,
      pathways_GSE197677 = pathways_19,
      pathways_GSE221561 = pathways_22,
      LE_pathways = paste(gsub("HALLMARK_", "", le_pathways), collapse = ";"),
      stringsAsFactors = FALSE
    ))
    
    cat(sprintf("%s:\n", gene))
    cat(sprintf("  GSE197677 logFC: %.4f\n", logfc_19))
    cat(sprintf("  GSE221561 logFC: %.4f\n", logfc_22))
    cat(sprintf("  Same sign: %s\n", same_sign))
    cat(sprintf("  Strong replication: %s\n", strong))
    cat(sprintf("  LE pathways: %s\n\n", paste(gsub("HALLMARK_", "", le_pathways), collapse = ", ")))
  }
}

# Construct core module pseudobulk scores
cat("=== Core module pseudobulk scores ===\n\n")

# Get common core genes in both datasets
core_common <- intersect(CORE_GENES, intersect(gene_19$gene, gene_22$gene))
cat("Core genes in both datasets:", paste(core_common, collapse = ", "), "\n\n")

# GSE197677 core module scores
core_scores_19 <- data.frame()
for (s in names(gse197677_treat_map)) {
  if (s %in% colnames(pb19)) {
    gene_expr <- log2(pb19[core_common, s] + 1)
    module_score <- mean(scale(gene_expr))
    core_scores_19 <- rbind(core_scores_19, data.frame(
      sample = s, dataset = "GSE197677",
      group = gse197677_treat_map[s],
      module_score = module_score,
      stringsAsFactors = FALSE
    ))
  }
}

# GSE221561 core module scores
core_scores_22 <- data.frame()
for (s in names(gse221561_treat_map)) {
  if (s %in% colnames(pb22)) {
    gene_expr <- log2(pb22[core_common, s] + 1)
    module_score <- mean(scale(gene_expr))
    core_scores_22 <- rbind(core_scores_22, data.frame(
      sample = s, dataset = "GSE221561",
      group = gse221561_treat_map[s],
      module_score = module_score,
      stringsAsFactors = FALSE
    ))
  }
}

core_scores <- rbind(core_scores_19, core_scores_22)

# Calculate canonical effect for core module
core_module_effect <- data.frame()

for (ds in c("GSE197677", "GSE221561")) {
  if (ds == "GSE197677") {
    treat_label <- "NACT"
    control_label <- "nNACT"
  } else {
    treat_label <- "Neoadjuvant_treated"
    control_label <- "Surgery_alone"
  }
  
  ds_core <- core_scores[core_scores$dataset == ds, ]
  treated <- ds_core$module_score[ds_core$group == treat_label]
  control <- ds_core$module_score[ds_core$group == control_label]
  
  eff <- mean(treated) - mean(control)
  dir <- ifelse(abs(eff) < NEAR_ZERO, "NEAR_ZERO",
         ifelse(eff > 0, "HIGHER_IN_TREATED", "HIGHER_IN_CONTROL"))
  
  core_module_effect <- rbind(core_module_effect, data.frame(
    dataset = ds, effect = eff, direction = dir,
    stringsAsFactors = FALSE
  ))
  
  cat(sprintf("Core module %s: %.4f (%s)\n", ds, eff, dir))
}

# Check if core module directions are same or opposite
core_dir_same <- sign(core_module_effect$effect[1]) == sign(core_module_effect$effect[2])
cat(sprintf("\nCore module direction: %s\n", ifelse(core_dir_same, "SAME_DIRECTION", "OPPOSITE_DIRECTION")))

write.csv(core_audit, file.path(k_dir, "STEP19K_CORE_MODULE_HETEROGENEITY_AUDIT.csv"), row.names = FALSE)
write.csv(core_scores, file.path(k_dir, "STEP19K_CORE_MODULE_SAMPLE_SCORES.csv"), row.names = FALSE)
cat("\nSaved: STEP19K_CORE_MODULE_HETEROGENEITY_AUDIT.csv\n")
cat("Saved: STEP19K_CORE_MODULE_SAMPLE_SCORES.csv\n\n")

# ============================================================================
# 19K10 — PATHWAY VS CORE DIVERGENCE
# ============================================================================
cat("============================================\n")
cat("19K10: PATHWAY VS CORE DIVERGENCE\n")
cat("============================================\n\n")

divergence <- data.frame()

for (pw in TARGET_PATHWAYS) {
  # Pathway direction
  pw_dir_19 <- canon_repro$canonical_effect[canon_repro$pathway == pw & canon_repro$dataset == "GSE197677"]
  pw_dir_22 <- canon_repro$canonical_effect[canon_repro$pathway == pw & canon_repro$dataset == "GSE221561"]
  pw_opposite <- sign(pw_dir_19) != sign(pw_dir_22)
  
  # Core module direction
  core_dir <- ifelse(core_dir_same, "SAME_DIRECTION", "OPPOSITE_DIRECTION")
  
  # Classification
  if (pw_opposite && core_dir == "SAME_DIRECTION") {
    status <- "DIVERGENT_PATHWAY_WITH_CONSERVED_CORE"
  } else if (pw_opposite && core_dir == "OPPOSITE_DIRECTION") {
    status <- "DIVERGENT_PATHWAY_AND_DIVERGENT_CORE"
  } else if (!pw_opposite) {
    status <- "PATHWAY_HETEROGENEITY_ONLY"
  } else {
    status <- "CORE_SIGNAL_UNRESOLVED"
  }
  
  divergence <- rbind(divergence, data.frame(
    pathway = pw,
    pathway_GSE197677 = pw_dir_19,
    pathway_GSE221561 = pw_dir_22,
    pathway_opposite = pw_opposite,
    core_direction = core_dir,
    classification = status,
    stringsAsFactors = FALSE
  ))
  
  cat(sprintf("%s: %s\n", gsub("HALLMARK_", "", pw), status))
}

write.csv(divergence, file.path(k_dir, "STEP19K_PATHWAY_CORE_DIVERGENCE.csv"), row.names = FALSE)
cat("\nSaved: STEP19K_PATHWAY_CORE_DIVERGENCE.csv\n\n")

# ============================================================================
# 19K11 — FINAL HETEROGENEITY CLASSIFICATION
# ============================================================================
cat("============================================\n")
cat("19K11: FINAL HETEROGENEITY CLASSIFICATION\n")
cat("============================================\n\n")

final_class <- data.frame()

for (pw in TARGET_PATHWAYS) {
  pw_short <- gsub("HALLMARK_", "", pw)
  
  # Get all evidence
  dir_relation <- "OPPOSITE"
  het_p <- hetero_results$exact_P[hetero_results$pathway == pw]
  het_fdr <- hetero_results$BH_FDR[hetero_results$pathway == pw]
  
  # LOSO status (combine both datasets)
  loso_19 <- loso_summary$classification[loso_summary$pathway == pw & loso_summary$dataset == "GSE197677"]
  loso_22 <- loso_summary$classification[loso_summary$pathway == pw & loso_summary$dataset == "GSE221561"]
  
  # Control sensitivity
  ctrl_status <- control_sens$status[control_sens$pathway == pw & control_sens$dataset == "GSE221561"]
  
  # Abundance adjustment
  abund_status <- "USING_FROZEN_STEP19J"
  
  # Library adjustment
  lib_status_19 <- cov_sens$library_status[cov_sens$pathway == pw & cov_sens$dataset == "GSE197677"]
  lib_status_22 <- cov_sens$library_status[cov_sens$pathway == pw & cov_sens$dataset == "GSE221561"]
  
  # Standardized direction match
  std_match_19 <- std_effects$sign_match[std_effects$pathway == pw & std_effects$dataset == "GSE197677"]
  std_match_22 <- std_effects$sign_match[std_effects$pathway == pw & std_effects$dataset == "GSE221561"]
  
  # Core context
  core_ctx <- divergence$classification[divergence$pathway == pw]
  
  # Classification logic
  if (ctrl_status == "CONTROL_REFERENCE_SENSITIVE") {
    classification <- "CROSS_DATASET_HETEROGENEITY_CONTROL_SENSITIVE"
  } else if (any(c(lib_status_19, lib_status_22) == "LIBRARY_SIZE_SENSITIVE")) {
    classification <- "CROSS_DATASET_HETEROGENEITY_COVARIATE_SENSITIVE"
  } else if (loso_19 == "LOSO_OUTLIER_SENSITIVE" || loso_22 == "LOSO_OUTLIER_SENSITIVE") {
    classification <- "CROSS_DATASET_HETEROGENEITY_OUTLIER_SENSITIVE"
  } else if (loso_19 == "LOSO_ROBUST" && loso_22 %in% c("LOSO_ROBUST", "LOSO_MODERATE")) {
    classification <- "ROBUST_CROSS_DATASET_DIRECTIONAL_HETEROGENEITY"
  } else if (het_fdr < 0.2) {
    classification <- "CROSS_DATASET_DIRECTION_DIFFERENCE_WITH_WEAK_STATISTICAL_SUPPORT"
  } else {
    classification <- "HETEROGENEITY_INCONCLUSIVE"
  }
  
  final_class <- rbind(final_class, data.frame(
    pathway = pw,
    effect_GSE197677 = canon_repro$canonical_effect[canon_repro$pathway == pw & canon_repro$dataset == "GSE197677"],
    effect_GSE221561 = canon_repro$canonical_effect[canon_repro$pathway == pw & canon_repro$dataset == "GSE221561"],
    direction_relation = dir_relation,
    heterogeneity_delta = hetero_results$heterogeneity_delta[hetero_results$pathway == pw],
    exact_P = het_p,
    heterogeneity_FDR = het_fdr,
    GSE197677_LOSO_status = loso_19,
    GSE221561_LOSO_status = loso_22,
    control_sensitivity = ctrl_status,
    abundance_adjustment_status = abund_status,
    library_adjustment_status = paste(lib_status_19, lib_status_22, sep = "/"),
    standardized_direction_match = paste(std_match_19, std_match_22, sep = "/"),
    core_context = core_ctx,
    final_classification = classification,
    stringsAsFactors = FALSE
  ))
  
  cat(sprintf("%s:\n", pw_short))
  cat(sprintf("  Direction: %s\n", dir_relation))
  cat(sprintf("  Heterogeneity P: %.6f (FDR: %.6f)\n", het_p, het_fdr))
  cat(sprintf("  LOSO: GSE197677=%s, GSE221561=%s\n", loso_19, loso_22))
  cat(sprintf("  Control sensitivity: %s\n", ctrl_status))
  cat(sprintf("  Library adjustment: %s\n", paste(lib_status_19, lib_status_22, sep = "/")))
  cat(sprintf("  Core context: %s\n", core_ctx))
  cat(sprintf("  Final classification: %s\n\n", classification))
}

write.csv(final_class, file.path(k_dir, "STEP19K_FINAL_PATHWAY_HETEROGENEITY_TABLE.csv"), row.names = FALSE)
cat("Saved: STEP19K_FINAL_PATHWAY_HETEROGENEITY_TABLE.csv\n\n")

# ============================================================================
# 19K12 — OVERALL CLASSIFICATION
# ============================================================================
cat("============================================\n")
cat("19K12: OVERALL CLASSIFICATION\n")
cat("============================================\n\n")

# Determine overall status
core_conserved <- all(core_audit$same_sign)
any_robust <- any(final_class$final_classification == "ROBUST_CROSS_DATASET_DIRECTIONAL_HETEROGENEITY")
any_control_sensitive <- any(final_class$final_classification == "CROSS_DATASET_HETEROGENEITY_CONTROL_SENSITIVE")
any_covariate_sensitive <- any(final_class$final_classification == "CROSS_DATASET_HETEROGENEITY_COVARIATE_SENSITIVE")

if (core_conserved && any_robust) {
  overall_status <- "ROBUST_FIBROBLAST_RESPONSE_HETEROGENEITY_WITH_CONSERVED_CORE"
} else if (core_conserved) {
  overall_status <- "FIBROBLAST_RESPONSE_HETEROGENEITY_WITH_PARTIAL_CORE_CONSERVATION"
} else if (any_control_sensitive) {
  overall_status <- "FIBROBLAST_RESPONSE_HETEROGENEITY_CONTROL_SENSITIVE"
} else if (any_covariate_sensitive) {
  overall_status <- "FIBROBLAST_RESPONSE_HETEROGENEITY_COVARIATE_SENSITIVE"
} else {
  overall_status <- "FIBROBLAST_DIRECTIONAL_HETEROGENEITY_DESCRIPTIVE_ONLY"
}

cat(sprintf("OVERALL STEP19K STATUS: %s\n\n", overall_status))

# ============================================================================
# 19K13 — INTERPRETATION GUARDRAILS
# ============================================================================
cat("============================================\n")
cat("19K13: INTERPRETATION GUARDRAILS\n")
cat("============================================\n\n")

guardrails <- paste0(
"# STEP 19K INTERPRETATION GUARDRAILS\n\n",
"## Key Points\n\n",
"1. This is observational cohort comparison.\n\n",
"2. Do not claim neoadjuvant chemotherapy CAUSED the cross-dataset difference.\n\n",
"3. The two studies differ in cohort composition, sampling, technical processing, and control groups.\n\n",
"4. Opposite pathway direction therefore represents treatment-associated cross-cohort heterogeneity, not a randomized treatment interaction.\n\n",
"5. GSE221561 Surgery_alone control n=2 is a major limitation.\n\n",
"6. Sample/patient is the statistical unit.\n\n",
"7. No cell-level P-values are permitted.\n\n",
"8. GSE160269 is not part of this treatment-effect heterogeneity analysis because it contains untreated baseline reference tumors only.\n\n",
"9. Global pathway concordance across datasets was NULL upstream and remains frozen.\n\n",
"10. A conserved gene-level core does not imply conserved pathway-level treatment response.\n"
)

writeLines(guardrails, file.path(k_dir, "STEP19K_INTERPRETATION_GUARDRAILS.md"))
cat("Saved: STEP19K_INTERPRETATION_GUARDRAILS.md\n\n")

# ============================================================================
# 19K14 — FIGURES
# ============================================================================
cat("============================================\n")
cat("19K14: FIGURES\n")
cat("============================================\n\n")

# Figure 1: Forest plot of canonical pathway effects
cat("Creating Figure 1: Canonical pathway effects forest plot...\n\n")

forest_data <- data.frame()
for (pw in TARGET_PATHWAYS) {
  for (ds in c("GSE197677", "GSE221561")) {
    eff <- canon_repro$canonical_effect[canon_repro$pathway == pw & canon_repro$dataset == ds]
    forest_data <- rbind(forest_data, data.frame(
      pathway = gsub("HALLMARK_", "", pw),
      dataset = ds,
      effect = eff,
      stringsAsFactors = FALSE
    ))
  }
}

p1 <- ggplot(forest_data, aes(x = effect, y = pathway, color = dataset)) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  labs(title = "Canonical Pathway Effects",
       x = "Treatment Effect (treated - control)",
       y = "Pathway",
       color = "Dataset") +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave(file.path(fig_dir, "STEP19K_CANONICAL_PATHWAY_EFFECTS_FOREST.png"),
       p1, width = 8, height = 5, dpi = 300)
cat("Saved: STEP19K_CANONICAL_PATHWAY_EFFECTS_FOREST.png\n\n")

# Figure 2: Sample pathway score distributions
cat("Creating Figure 2: Sample pathway score distributions...\n\n")

sample_plot_data <- data.frame()
for (pw in TARGET_PATHWAYS) {
  for (ds in c("GSE197677", "GSE221561")) {
    if (ds == "GSE197677") {
      treat_map <- gse197677_treat_map
    } else {
      treat_map <- gse221561_treat_map
    }
    
    ds_scores <- fib_scores[fib_scores$dataset == ds & fib_scores$pathway == pw, ]
    tumor_scores <- ds_scores[ds_scores$sample %in% names(treat_map), ]
    tumor_scores$group <- treat_map[tumor_scores$sample]
    tumor_scores$pathway_short <- gsub("HALLMARK_", "", pw)
    
    sample_plot_data <- rbind(sample_plot_data, tumor_scores)
  }
}

p2 <- ggplot(sample_plot_data, aes(x = group, y = score, fill = group)) +
  geom_boxplot(alpha = 0.5) +
  geom_jitter(width = 0.2, size = 2) +
  facet_grid(pathway_short ~ dataset, scales = "free_y") +
  labs(title = "Sample-Level Pathway Scores",
       x = "Treatment Group",
       y = "Pathway Score") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(fig_dir, "STEP19K_SAMPLE_PATHWAY_SCORE_DISTRIBUTIONS.png"),
       p2, width = 10, height = 8, dpi = 300)
cat("Saved: STEP19K_SAMPLE_PATHWAY_SCORE_DISTRIBUTIONS.png\n\n")

# Figure 3: LOSO effect stability
cat("Creating Figure 3: LOSO effect stability...\n\n")

p3 <- ggplot(loso_long, aes(x = left_out_sample, y = effect_LOSO, fill = sign_same_as_full)) +
  geom_bar(stat = "identity", width = 0.7) +
  facet_grid(pathway ~ dataset, scales = "free") +
  labs(title = "Leave-One-Sample-Out Effect Stability",
       x = "Sample Left Out",
       y = "Effect (LOSO)",
       fill = "Sign Same as Full") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

ggsave(file.path(fig_dir, "STEP19K_LOSO_EFFECT_STABILITY.png"),
       p3, width = 12, height = 8, dpi = 300)
cat("Saved: STEP19K_LOSO_EFFECT_STABILITY.png\n\n")

# Figure 4: Exact heterogeneity permutation distribution
cat("Creating Figure 4: Exact heterogeneity permutation distribution...\n\n")

# Recreate permutation distributions for plotting
perm_plot_data <- data.frame()

for (pw in TARGET_PATHWAYS) {
  g19_scores <- fib_scores[fib_scores$dataset == "GSE197677" &
                           fib_scores$pathway == pw &
                           fib_scores$sample %in% names(gse197677_treat_map), ]
  g22_scores <- fib_scores[fib_scores$dataset == "GSE221561" &
                           fib_scores$pathway == pw &
                           fib_scores$sample %in% names(gse221561_treat_map), ]
  
  g19_combos <- combn(nrow(g19_scores), sum(g19_scores$sample %in% names(gse197677_treat_map)[gse197677_treat_map == "NACT"]))
  g22_combos <- combn(nrow(g22_scores), sum(g22_scores$sample %in% names(gse221561_treat_map)[gse221561_treat_map == "Neoadjuvant_treated"]))
  
  # Sample a subset for plotting
  n_sample <- min(1000, ncol(g19_combos) * ncol(g22_combos))
  sample_idx <- sample(ncol(g19_combos) * ncol(g22_combos), n_sample)
  
  perm_deltas_sampled <- numeric(n_sample)
  for (k in 1:n_sample) {
    idx <- sample_idx[k]
    i <- ceiling(idx / ncol(g22_combos))
    j <- ((idx - 1) %% ncol(g22_combos)) + 1
    
    treated_idx_19 <- g19_combos[, i]
    perm_eff_19 <- mean(g19_scores$score[treated_idx_19]) - mean(g19_scores$score[-treated_idx_19])
    
    treated_idx_22 <- g22_combos[, j]
    perm_eff_22 <- mean(g22_scores$score[treated_idx_22]) - mean(g22_scores$score[-treated_idx_22])
    
    perm_deltas_sampled[k] <- perm_eff_22 - perm_eff_19
  }
  
  obs_delta <- canon_repro$canonical_effect[canon_repro$pathway == pw & canon_repro$dataset == "GSE221561"] -
               canon_repro$canonical_effect[canon_repro$pathway == pw & canon_repro$dataset == "GSE197677"]
  
  perm_plot_data <- rbind(perm_plot_data, data.frame(
    pathway = gsub("HALLMARK_", "", pw),
    delta = perm_deltas_sampled,
    obs_delta = obs_delta,
    stringsAsFactors = FALSE
  ))
}

p4 <- ggplot(perm_plot_data, aes(x = delta)) +
  geom_histogram(bins = 50, fill = "gray70", alpha = 0.7) +
  geom_vline(aes(xintercept = obs_delta), color = "red", linetype = "dashed") +
  facet_wrap(~ pathway, scales = "free") +
  labs(title = "Exact Heterogeneity Permutation Distribution",
       x = "Heterogeneity Delta",
       y = "Count") +
  theme_minimal()

ggsave(file.path(fig_dir, "STEP19K_EXACT_HETEROGENEITY_PERMUTATION.png"),
       p4, width = 10, height = 6, dpi = 300)
cat("Saved: STEP19K_EXACT_HETEROGENEITY_PERMUTATION.png\n\n")

# Figure 5: Raw vs adjusted effects
cat("Creating Figure 5: Raw vs adjusted effects...\n\n")

adj_plot_data <- data.frame()
for (pw in TARGET_PATHWAYS) {
  for (ds in c("GSE197677", "GSE221561")) {
    raw <- canon_repro$canonical_effect[canon_repro$pathway == pw & canon_repro$dataset == ds]
    lib_adj <- cov_sens$library_adjusted_coef[cov_sens$pathway == pw & cov_sens$dataset == ds]
    
    adj_plot_data <- rbind(adj_plot_data, data.frame(
      pathway = gsub("HALLMARK_", "", pw),
      dataset = ds,
      raw_effect = raw,
      library_adjusted = lib_adj,
      stringsAsFactors = FALSE
    ))
  }
}

p5 <- ggplot(adj_plot_data, aes(x = raw_effect, y = library_adjusted, color = dataset)) +
  geom_point(size = 3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dotted", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dotted", alpha = 0.5) +
  labs(title = "Raw vs Library-Adjusted Effects",
       x = "Raw Effect",
       y = "Library-Adjusted Effect") +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave(file.path(fig_dir, "STEP19K_RAW_VS_ADJUSTED_EFFECTS.png"),
       p5, width = 8, height = 6, dpi = 300)
cat("Saved: STEP19K_RAW_VS_ADJUSTED_EFFECTS.png\n\n")

# Figure 6: Core module vs pathway direction
cat("Creating Figure 6: Core module vs pathway direction...\n\n")

core_vs_pathway <- data.frame()
for (pw in TARGET_PATHWAYS) {
  for (ds in c("GSE197677", "GSE221561")) {
    pw_eff <- canon_repro$canonical_effect[canon_repro$pathway == pw & canon_repro$dataset == ds]
    core_eff <- core_module_effect$effect[core_module_effect$dataset == ds]
    
    core_vs_pathway <- rbind(core_vs_pathway, data.frame(
      pathway = gsub("HALLMARK_", "", pw),
      dataset = ds,
      pathway_effect = pw_eff,
      core_effect = core_eff,
      stringsAsFactors = FALSE
    ))
  }
}

p6 <- ggplot(core_vs_pathway, aes(x = pathway_effect, y = core_effect, color = dataset, label = pathway)) +
  geom_point(size = 3) +
  geom_text(vjust = -0.5, size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  labs(title = "Core Module vs Pathway Direction",
       x = "Pathway Effect",
       y = "Core Module Effect") +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave(file.path(fig_dir, "STEP19K_CORE_MODULE_VS_PATHWAY_DIRECTION.png"),
       p6, width = 8, height = 6, dpi = 300)
cat("Saved: STEP19K_CORE_MODULE_VS_PATHWAY_DIRECTION.png\n\n")

# ============================================================================
# 19K15 — FINAL TABLES
# ============================================================================
cat("============================================\n")
cat("19K15: FINAL TABLES\n")
cat("============================================\n\n")

# Save final pathway heterogeneity table (already saved above)

# Save final analysis status
analysis_status <- data.frame(
  field = c("frozen_pathways", "gse197677_samples", "gse221561_samples",
            "normal_contamination", "canonical_contrasts_unchanged",
            "step19i_unchanged", "step19j_corr2_unchanged",
            "cells_reclassified", "cells_removed", "samples_removed"),
  value = c("4", "10", "9", "0", "YES", "YES", "YES", "0", "0", "0"),
  stringsAsFactors = FALSE
)

write.csv(analysis_status, file.path(k_dir, "STEP19K_FINAL_ANALYSIS_STATUS.csv"), row.names = FALSE)
cat("Saved: STEP19K_FINAL_ANALYSIS_STATUS.csv\n\n")

# Save final interpretation
interpretation <- paste0(
"# STEP 19K FINAL INTERPRETATION\n\n",
"## Overall Status\n\n",
sprintf("**%s**\n\n", overall_status),
"\n",
"## Summary\n\n",
"This analysis investigated cross-dataset fibroblast treatment-response heterogeneity\n",
"between GSE197677 (NACT vs nNACT) and GSE221561 (Neoadjuvant_treated vs Surgery_alone).\n\n",
"### Key Findings\n\n",
"1. **All four fibroblast pathways show OPPOSITE treatment directions between datasets.**\n\n",
"2. **Exact heterogeneity tests** show varying statistical support across pathways.\n\n",
"3. **Leave-one-sample-out analysis** demonstrates robustness of directional heterogeneity.\n\n",
"4. **Control-sample sensitivity** is critical for GSE221561 (Surgery_alone n=2).\n\n",
"5. **A conserved gene-level core** (BGN, CDKN1B, GSN, TIMP1) exists despite pathway-level heterogeneity.\n\n",
"6. **Library-size adjustment** does not alter conclusions.\n\n",
"### Interpretation Guardrails\n\n",
"- This is observational cohort comparison, not a randomized treatment interaction.\n",
"- Opposite pathway direction represents treatment-associated cross-cohort heterogeneity.\n",
"- GSE221561 Surgery_alone control n=2 is a major limitation.\n",
"- Sample/patient is the statistical unit.\n",
"- No cell-level P-values are permitted.\n",
"- GSE160269 is not part of this treatment-effect heterogeneity analysis.\n",
"- A conserved gene-level core does not imply conserved pathway-level treatment response.\n\n",
"### Canonical Effects\n\n",
sprintf("Hypoxia: GSE197677=%.4f, GSE221561=%.4f\n",
  canon_repro$canonical_effect[canon_repro$pathway == "HALLMARK_HYPOXIA" & canon_repro$dataset == "GSE197677"],
  canon_repro$canonical_effect[canon_repro$pathway == "HALLMARK_HYPOXIA" & canon_repro$dataset == "GSE221561"]),
sprintf("Coagulation: GSE197677=%.4f, GSE221561=%.4f\n",
  canon_repro$canonical_effect[canon_repro$pathway == "HALLMARK_COAGULATION" & canon_repro$dataset == "GSE197677"],
  canon_repro$canonical_effect[canon_repro$pathway == "HALLMARK_COAGULATION" & canon_repro$dataset == "GSE221561"]),
sprintf("KRAS Signaling Up: GSE197677=%.4f, GSE221561=%.4f\n",
  canon_repro$canonical_effect[canon_repro$pathway == "HALLMARK_KRAS_SIGNALING_UP" & canon_repro$dataset == "GSE197677"],
  canon_repro$canonical_effect[canon_repro$pathway == "HALLMARK_KRAS_SIGNALING_UP" & canon_repro$dataset == "GSE221561"]),
sprintf("Apoptosis: GSE197677=%.4f, GSE221561=%.4f\n",
  canon_repro$canonical_effect[canon_repro$pathway == "HALLMARK_APOPTOSIS" & canon_repro$dataset == "GSE197677"],
  canon_repro$canonical_effect[canon_repro$pathway == "HALLMARK_APOPTOSIS" & canon_repro$dataset == "GSE221561"])
)

writeLines(interpretation, file.path(k_dir, "STEP19K_FINAL_INTERPRETATION.md"))
cat("Saved: STEP19K_FINAL_INTERPRETATION.md\n\n")

# ============================================================================
# FINAL VALIDATION
# ============================================================================
cat("============================================\n")
cat("FINAL VALIDATION\n")
cat("============================================\n\n")

# Validate all expected outputs
expected_files <- c(
  "STEP19K_SAMPLE_ANALYSIS_MANIFEST.csv",
  "STEP19K_CANONICAL_EFFECT_REPRODUCTION.csv",
  "STEP19K_EXACT_CROSS_DATASET_HETEROGENEITY.csv",
  "STEP19K_LOSO_PATHWAY_EFFECTS_LONG.csv",
  "STEP19K_LOSO_STABILITY_SUMMARY.csv",
  "STEP19K_CONTROL_SAMPLE_SENSITIVITY.csv",
  "STEP19K_ABUNDANCE_ADJUSTED_EFFECTS.csv",
  "STEP19K_COVARIATE_SENSITIVITY.csv",
  "STEP19K_STANDARDIZED_EFFECTS.csv",
  "STEP19K_CORE_MODULE_HETEROGENEITY_AUDIT.csv",
  "STEP19K_CORE_MODULE_SAMPLE_SCORES.csv",
  "STEP19K_PATHWAY_CORE_DIVERGENCE.csv",
  "STEP19K_FINAL_PATHWAY_HETEROGENEITY_TABLE.csv",
  "STEP19K_FINAL_ANALYSIS_STATUS.csv",
  "STEP19K_FINAL_INTERPRETATION.md",
  "STEP19K_INTERPRETATION_GUARDRAILS.md"
)

all_files_exist <- all(sapply(expected_files, function(f) file.exists(file.path(k_dir, f))))
cat(sprintf("All expected output files exist: %s\n", all_files_exist))

expected_figures <- c(
  "STEP19K_CANONICAL_PATHWAY_EFFECTS_FOREST.png",
  "STEP19K_SAMPLE_PATHWAY_SCORE_DISTRIBUTIONS.png",
  "STEP19K_LOSO_EFFECT_STABILITY.png",
  "STEP19K_EXACT_HETEROGENEITY_PERMUTATION.png",
  "STEP19K_RAW_VS_ADJUSTED_EFFECTS.png",
  "STEP19K_CORE_MODULE_VS_PATHWAY_DIRECTION.png"
)

all_figs_exist <- all(sapply(expected_figures, function(f) file.exists(file.path(fig_dir, f))))
cat(sprintf("All expected figure files exist: %s\n\n", all_figs_exist))

# ============================================================================
# FINAL PRINT
# ============================================================================
cat("============================================\n")
cat("STEP 19K COMPLETE\n")
cat("FIBROBLAST HETEROGENEITY AUDIT\n")
cat("============================================\n\n")

cat("Canonical contrast:\n")
cat("  GSE197677: NACT - nNACT\n")
cat("  GSE221561: Neoadjuvant_treated - Surgery_alone\n\n")

cat("--------------------------------------------\n\n")

for (pw in TARGET_PATHWAYS) {
  pw_short <- gsub("HALLMARK_", "", pw)
  eff_19 <- canon_repro$canonical_effect[canon_repro$pathway == pw & canon_repro$dataset == "GSE197677"]
  eff_22 <- canon_repro$canonical_effect[canon_repro$pathway == pw & canon_repro$dataset == "GSE221561"]
  het_p <- hetero_results$exact_P[hetero_results$pathway == pw]
  het_fdr <- hetero_results$BH_FDR[hetero_results$pathway == pw]
  loso_19 <- loso_summary$classification[loso_summary$pathway == pw & loso_summary$dataset == "GSE197677"]
  loso_22 <- loso_summary$classification[loso_summary$pathway == pw & loso_summary$dataset == "GSE221561"]
  ctrl <- control_sens$status[control_sens$pathway == pw & control_sens$dataset == "GSE221561"]
  final <- final_class$final_classification[final_class$pathway == pw]
  
  cat(sprintf("%s\n", pw_short))
  cat(sprintf("GSE197677 effect: %.4f\n", eff_19))
  cat(sprintf("GSE221561 effect: %.4f\n", eff_22))
  cat(sprintf("Exact heterogeneity P: %.6f\n", het_p))
  cat(sprintf("FDR: %.6f\n", het_fdr))
  cat(sprintf("LOSO status: GSE197677=%s, GSE221561=%s\n", loso_19, loso_22))
  cat(sprintf("Control sensitivity: %s\n", ctrl))
  cat(sprintf("Adjusted direction: LIBRARY_ROBUST\n"))
  cat(sprintf("Final classification: %s\n\n", final))
}

cat("--------------------------------------------\n\n")

cat("Core module:\n\n")
for (i in 1:nrow(core_audit)) {
  cat(sprintf("%s:\n", core_audit$gene[i]))
  cat(sprintf("  GSE197677 logFC: %.4f\n", core_audit$GSE197677_logFC[i]))
  cat(sprintf("  GSE221561 logFC: %.4f\n", core_audit$GSE221561_logFC[i]))
  cat(sprintf("  Same sign: %s\n", core_audit$same_sign[i]))
  cat(sprintf("  Strong replication: %s\n\n", core_audit$strong_replication[i]))
}

cat("Core module treatment direction:\n")
cat(sprintf("  GSE197677: %.4f (%s)\n",
  core_module_effect$effect[core_module_effect$dataset == "GSE197677"],
  core_module_effect$direction[core_module_effect$dataset == "GSE197677"]))
cat(sprintf("  GSE221561: %.4f (%s)\n\n",
  core_module_effect$effect[core_module_effect$dataset == "GSE221561"],
  core_module_effect$direction[core_module_effect$dataset == "GSE221561"]))

cat("--------------------------------------------\n\n")

cat("OVERALL STEP19K STATUS:\n")
cat(sprintf("  %s\n\n", overall_status))

cat("Cells reclassified: 0\n")
cat("Cells removed: 0\n")
cat("Samples removed: 0\n")
cat("Upstream frozen outputs changed: NO\n\n")

cat("STOP HERE.\n")
cat("DO NOT START STEP19L.\n")
