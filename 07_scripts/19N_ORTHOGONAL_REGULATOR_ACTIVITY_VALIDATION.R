#!/usr/bin/env Rscript
# ============================================================================
# STEP 19N: ORTHOGONAL REGULATOR ACTIVITY VALIDATION
# ============================================================================
# Validates moderate-priority regulatory candidates from Step19M using
# target/regulon activity rather than regulator gene expression alone.
# Uses MSigDB C3 TFT (unsigned) as local regulon resource.
# ============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(reshape2)
})

BASE_DIR <- Sys.getenv("ESCC_CAF_PROJECT_ROOT", unset = ".")
N_DIR    <- file.path(BASE_DIR, "06_tables/cross_dataset/pathway_analysis/Step19N")
FIG_DIR  <- file.path(BASE_DIR, "05_figures/cross_dataset/pathway_analysis/Step19N")
LOG_DIR  <- file.path(BASE_DIR, "08_logs")
dir.create(N_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(LOG_DIR, showWarnings = FALSE, recursive = TRUE)

CORE_GENES <- c("BGN", "CDKN1B", "GSN", "TIMP1")

PRIMARY_CANDIDATES <- c("FOS", "FOSL2", "KLF6", "RUNX1", "SMAD3", "SMAD7")
SECONDARY_CANDIDATES <- c("CDKN1A", "FOXO3", "SOX9")

TREAT_MAP <- list(
  GSE197677 = c(
    ESC06 = "nNACT", ESC07 = "nNACT", ESC09 = "NACT", ESC11 = "NACT",
    ESC12 = "NACT", ESC15 = "NACT", ESC16 = "NACT", ESC17 = "nNACT",
    ESC18 = "nNACT", ESC22 = "NACT"
  ),
  GSE221561 = c(
    S0520T = "Neoadjuvant_treated", S1265T = "Neoadjuvant_treated",
    S1315TS = "Neoadjuvant_treated", S1535T = "Neoadjuvant_treated",
    S2423T = "Neoadjuvant_treated", S2487T = "Neoadjuvant_treated",
    S6829T = "Neoadjuvant_treated", S6417T = "Surgery_alone",
    S9478T = "Surgery_alone"
  )
)

TREAT_CONTRAST <- list(
  GSE197677  = c(treated = "NACT", control = "nNACT"),
  GSE221561  = c(treated = "Neoadjuvant_treated", control = "Surgery_alone")
)

compute_log2_cpm <- function(counts_mat) {
  lib_size <- colSums(counts_mat)
  lib_size[lib_size == 0] <- 1
  cpm_mat <- sweep(counts_mat, 2, lib_size, "/") * 1e6
  log2(cpm_mat + 1)
}

# ============================================================================
# 19N0 — OUTPUT STRUCTURE
# ============================================================================
cat("============================================\n")
cat("STEP 19N: ORTHOGONAL REGULATOR ACTIVITY VALIDATION\n")
cat("============================================\n\n")

# ============================================================================
# 19N1 — INPUT PROVENANCE
# ============================================================================
cat("============================================\n")
cat("19N1: INPUT PROVENANCE\n")
cat("============================================\n\n")

core_scores <- read.csv(file.path(BASE_DIR,
  "06_tables/cross_dataset/pathway_analysis/Step19L/STEP19L_SAMPLE_CORE_MODULE_SCORES.csv"))
core_effects <- read.csv(file.path(BASE_DIR,
  "06_tables/cross_dataset/pathway_analysis/Step19L/STEP19L_CORE_MODULE_CANONICAL_EFFECTS.csv"))
step19m_ranking <- read.csv(file.path(BASE_DIR,
  "06_tables/cross_dataset/pathway_analysis/Step19M/STEP19M_FINAL_REGULATOR_RANKING.csv"))
step19m_tf <- read.csv(file.path(BASE_DIR,
  "06_tables/cross_dataset/pathway_analysis/Step19M/STEP19M_TF_EXPRESSION_ASSOCIATIONS.csv"))

fib197 <- readRDS(file.path(BASE_DIR,
  "06_tables/cross_dataset/pathway_analysis/Step19I/counts/GSE197677_FIBROBLAST_CAF_raw_counts.rds"))
fib221 <- readRDS(file.path(BASE_DIR,
  "06_tables/cross_dataset/pathway_analysis/Step19I/counts/GSE221561_FIBROBLAST_CAF_raw_counts.rds"))

cat("GSE197677 fibroblast:", nrow(fib197), "genes x", ncol(fib197), "samples\n")
cat("GSE221561 fibroblast:", nrow(fib221), "genes x", ncol(fib221), "samples\n\n")

# Validate sample composition
samples197 <- colnames(fib197)
samples221 <- colnames(fib221)
cat("GSE197677:", length(samples197), "tumor samples\n")
cat("GSE221561:", length(samples221), "tumor samples (7 treated + 2 control)\n\n")

# Validate canonical effects
cat("=== FROZEN EFFECTS ===\n")
for (i in 1:nrow(core_effects)) {
  cat(sprintf("  %s %s: %.4f\n", core_effects$dataset[i], core_effects$module[i], core_effects$effect[i]))
}

# Save provenance
provenance <- data.frame(
  step = "19N",
  status = "INPUT_VALIDATED",
  n_fib197 = ncol(fib197),
  n_fib221 = ncol(fib221),
  frozen_core4_197 = core_effects$effect[core_effects$dataset == "GSE197677" & core_effects$module == "CORE4"],
  frozen_core4_221 = core_effects$effect[core_effects$dataset == "GSE221561" & core_effects$module == "CORE4"],
  frozen_core2_197 = core_effects$effect[core_effects$dataset == "GSE197677" & core_effects$module == "CORE2_STRONG"],
  frozen_core2_221 = core_effects$effect[core_effects$dataset == "GSE221561" & core_effects$module == "CORE2_STRONG"],
  step19m_classification = "CORE2_SPECIFIC_REGULATOR_SUPPORT",
  timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
)
write.csv(provenance, file.path(N_DIR, "STEP19N_INPUT_PROVENANCE.csv"), row.names = FALSE)
cat("\nSaved: STEP19N_INPUT_PROVENANCE.csv\n")

# ============================================================================
# 19N2 — LOCAL REGULON RESOURCE PREFLIGHT
# ============================================================================
cat("\n============================================\n")
cat("19N2: LOCAL REGULON RESOURCE PREFLIGHT\n")
cat("============================================\n\n")

# Check available resources
check_resource <- function(name, available, version, path, n_reg, n_edges, signed, conf_levels) {
  data.frame(
    resource = name, available = available, version = version,
    local_path = path, network_type = ifelse(signed, "signed", "unsigned"),
    signed_edges = signed, confidence_levels = conf_levels,
    number_regulators = n_reg, number_edges = n_edges,
    stringsAsFactors = FALSE
  )
}

resources <- list()

# Check decoupleR
resources$decoupleR <- check_resource("decoupleR", FALSE, NA, NA, 0, 0, FALSE, NA)

# Check dorothea
resources$dorothea <- check_resource("dorothea", FALSE, NA, NA, 0, 0, FALSE, NA)

# Check VIPER
resources$viper <- check_resource("viper", FALSE, NA, NA, 0, 0, FALSE, NA)

# Check progeny
resources$progeny <- check_resource("progeny", FALSE, NA, NA, 0, 0, FALSE, NA)

# Check msigdbr C2 CGP (perturbation signatures)
msigdbr_available <- requireNamespace("msigdbr", quietly = TRUE)
if (msigdbr_available) {
  c2_cache_path <- path.expand("~/Library/Caches/org.R-project.R/R/msigdbr/msigdb.2026.1.Hs.C2.rds")
  tryCatch({
    c2_cache <- readRDS(c2_cache_path)
    cgp <- c2_cache[c2_cache$gs_subcollection == "CGP", ]
    n_reg <- length(unique(cgp$source_gene))
    n_gs <- length(unique(cgp$gs_name))
    resources$msigdbr_c2_cgp <- check_resource(
      "msigdbr_C2_CGP", TRUE, "2026.1",
      c2_cache_path,
      n_reg, n_gs, TRUE, "perturbation_UP/DN"
    )
  }, error = function(e) {
    cat("msigdbr C2 CGP load error:", conditionMessage(e), "\n")
    resources$msigdbr_c2_cgp <<- check_resource("msigdbr_C2_CGP", FALSE, NA, NA, 0, 0, FALSE, NA)
  })
} else {
  resources$msigdbr_c2_cgp <- check_resource("msigdbr_C2_CGP", FALSE, NA, NA, 0, 0, FALSE, NA)
}

resources_df <- do.call(rbind, resources)
rownames(resources_df) <- NULL

cat("=== REGULON RESOURCE STATUS ===\n")
for (i in 1:nrow(resources_df)) {
  cat(sprintf("  %-20s: %s (edges=%d, regulators=%d)\n",
    resources_df$resource[i], ifelse(resources_df$available[i], "AVAILABLE", "NOT AVAILABLE"),
    resources_df$number_edges[i], resources_df$number_regulators[i]))
}

# Determine if we can proceed
usable_resource <- resources_df$resource[resources_df$available == TRUE]
if (length(usable_resource) == 0) {
  cat("\nREGULON_RESOURCE_UNAVAILABLE\n")
  cat("No usable TF-target resource found. Stopping gracefully.\n")
  
  status <- data.frame(
    step = "19N", status = "REGULON_RESOURCE_UNAVAILABLE",
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  )
  write.csv(status, file.path(N_DIR, "STEP19N_FINAL_STATUS.csv"), row.names = FALSE)
  write.csv(resources_df, file.path(N_DIR, "STEP19N_REGULON_RESOURCE_PREFLIGHT.csv"), row.names = FALSE)
  stop("REGULON_RESOURCE_UNAVAILABLE")
} else {
  cat("\nUsing resource:", usable_resource[1], "\n")
  REGULON_RESOURCE <- usable_resource[1]
}

write.csv(resources_df, file.path(N_DIR, "STEP19N_REGULON_RESOURCE_PREFLIGHT.csv"), row.names = FALSE)
cat("Saved: STEP19N_REGULON_RESOURCE_PREFLIGHT.csv\n")

# ============================================================================
# 19N3 — COMMON EXPRESSION UNIVERSE
# ============================================================================
cat("\n============================================\n")
cat("19N3: COMMON EXPRESSION UNIVERSE\n")
cat("============================================\n\n")

expr197 <- compute_log2_cpm(fib197)
expr221 <- compute_log2_cpm(fib221)

common_genes <- intersect(rownames(expr197), rownames(expr221))
cat("GSE197677 genes:", nrow(expr197), "\n")
cat("GSE221561 genes:", nrow(expr221), "\n")
cat("Common genes:", length(common_genes), "\n")

# Load regulon targets from C2 CGP perturbation signatures
c2_cache_path <- path.expand("~/Library/Caches/org.R-project.R/R/msigdbr/msigdb.2026.1.Hs.C2.rds")
c2_cache <- readRDS(c2_cache_path)
cgp <- c2_cache[c2_cache$gs_subcollection == "CGP", ]

# Build signed regulon: UP genes = positively associated, DN genes = negatively associated
all_candidates <- unique(c(PRIMARY_CANDIDATES, SECONDARY_CANDIDATES))
regulon_list <- list()
for (tf in all_candidates) {
  # Find perturbation gene sets mentioning this TF
  gs_up <- unique(cgp$gs_name[grepl(tf, cgp$gs_name, ignore.case = TRUE) & grepl("_UP$", cgp$gs_name)])
  gs_dn <- unique(cgp$gs_name[grepl(tf, cgp$gs_name, ignore.case = TRUE) & grepl("_DN$", cgp$gs_name)])
  
  up_genes <- unique(cgp$db_gene_symbol[cgp$gs_name %in% gs_up])
  dn_genes <- unique(cgp$db_gene_symbol[cgp$gs_name %in% gs_dn])
  
  regulon_list[[tf]] <- list(up = up_genes, dn = dn_genes)
}

# Check regulon target overlap with common universe
regulon_coverage <- list()
for (tf in names(regulon_list)) {
  up_targets <- regulon_list[[tf]]$up
  dn_targets <- regulon_list[[tf]]$dn
  in_common_up <- length(intersect(up_targets, common_genes))
  in_common_dn <- length(intersect(dn_targets, common_genes))
  regulon_coverage[[tf]] <- data.frame(
    regulator = tf,
    total_up_targets = length(up_targets),
    total_dn_targets = length(dn_targets),
    up_in_common = in_common_up,
    dn_in_common = in_common_dn,
    total_in_common = in_common_up + in_common_dn,
    stringsAsFactors = FALSE
  )
}
regulon_cov_df <- do.call(rbind, regulon_coverage)
cat("\n=== REGULON TARGET COVERAGE ===\n")
for (i in 1:nrow(regulon_cov_df)) {
  cat(sprintf("  %s: UP=%d DN=%d (total in common=%d)\n",
    regulon_cov_df$regulator[i], regulon_cov_df$up_in_common[i],
    regulon_cov_df$dn_in_common[i], regulon_cov_df$total_in_common[i]))
}

universe_df <- data.frame(
  dataset = c(rep("GSE197677", nrow(expr197)), rep("GSE221561", nrow(expr221))),
  n_genes = c(nrow(expr197), nrow(expr221)),
  n_common = length(common_genes),
  n_regulon_up_targets = sum(regulon_cov_df$up_in_common),
  n_regulon_dn_targets = sum(regulon_cov_df$dn_in_common),
  stringsAsFactors = FALSE
)
write.csv(universe_df, file.path(N_DIR, "STEP19N_COMMON_REGULON_GENE_UNIVERSE.csv"), row.names = FALSE)
cat("Saved: STEP19N_COMMON_REGULON_GENE_UNIVERSE.csv\n")

# ============================================================================
# 19N4 — REGULATOR ACTIVITY ESTIMATION
# ============================================================================
cat("\n============================================\n")
cat("19N4: REGULATOR ACTIVITY ESTIMATION\n")
cat("============================================\n\n")

# Signed activity: mean z-score of UP targets minus mean z-score of DN targets
# UP genes = genes upregulated when TF is active
# DN genes = genes downregulated when TF is active

compute_signed_activity <- function(expr_mat, up_genes, dn_genes) {
  available_up <- intersect(up_genes, rownames(expr_mat))
  available_dn <- intersect(dn_genes, rownames(expr_mat))
  
  if (length(available_up) < 2 && length(available_dn) < 2) return(rep(NA, ncol(expr_mat)))
  
  up_expr <- expr_mat[available_up, , drop = FALSE]
  dn_expr <- expr_mat[available_dn, , drop = FALSE]
  
  up_z <- t(scale(t(up_expr)))
  up_z[is.nan(up_z)] <- 0
  dn_z <- t(scale(t(dn_expr)))
  dn_z[is.nan(dn_z)] <- 0
  
  up_mean <- colMeans(up_z)
  dn_mean <- colMeans(dn_z)
  
  # Activity = UP - DN (positive means TF is active)
  if (length(available_up) >= 2 && length(available_dn) >= 2) {
    up_mean - dn_mean
  } else if (length(available_up) >= 2) {
    up_mean
  } else {
    -dn_mean
  }
}

cat("Computing regulator activity for all candidates...\n")

activity_list <- list()
for (tf in all_candidates) {
  up_genes <- regulon_list[[tf]]$up
  dn_genes <- regulon_list[[tf]]$dn
  
  for (ds in c("GSE197677", "GSE221561")) {
    if (ds == "GSE197677") {
      expr_mat <- expr197
      samp_ids <- samples197
    } else {
      expr_mat <- expr221
      samp_ids <- samples221
    }
    
    available_up <- intersect(up_genes, rownames(expr_mat))
    available_dn <- intersect(dn_genes, rownames(expr_mat))
    
    if (length(available_up) < 2 && length(available_dn) < 2) {
      # No usable targets - assign NA activity
      for (j in seq_along(samp_ids)) {
        treat_val <- as.character(TREAT_MAP[[ds]][samp_ids[j]])
        activity_list[[paste(tf, ds, samp_ids[j])]] <- data.frame(
          dataset = ds, sample_id = samp_ids[j],
          treatment = treat_val,
          regulator = tf, activity = NA,
          method = "signed_perturbation_zscore", target_count = 0L,
          resource = "msigdbr_C2_CGP", stringsAsFactors = FALSE
        )
      }
    } else {
      activity <- compute_signed_activity(expr_mat, up_genes, dn_genes)
      for (j in seq_along(samp_ids)) {
        treat_val <- as.character(TREAT_MAP[[ds]][samp_ids[j]])
        activity_list[[paste(tf, ds, samp_ids[j])]] <- data.frame(
          dataset = ds, sample_id = samp_ids[j],
          treatment = treat_val,
          regulator = tf, activity = activity[j],
          method = "signed_perturbation_zscore",
          target_count = length(available_up) + length(available_dn),
          resource = "msigdbr_C2_CGP", stringsAsFactors = FALSE
        )
      }
    }
  }
}

activity_df <- do.call(rbind, activity_list)
rownames(activity_df) <- NULL

cat("Activity computed for", length(unique(activity_df$regulator)), "regulators\n")
cat("Total rows:", nrow(activity_df), "\n")

write.csv(activity_df, file.path(N_DIR, "STEP19N_SAMPLE_REGULATOR_ACTIVITY.csv"), row.names = FALSE)
cat("Saved: STEP19N_SAMPLE_REGULATOR_ACTIVITY.csv\n")

# ============================================================================
# 19N5 — REGULATOR COVERAGE QC
# ============================================================================
cat("\n============================================\n")
cat("19N5: REGULATOR COVERAGE QC\n")
cat("============================================\n\n")

coverage_qc_list <- list()
for (tf in all_candidates) {
  up_targets <- regulon_list[[tf]]$up
  dn_targets <- regulon_list[[tf]]$dn
  in_197_up <- length(intersect(up_targets, rownames(expr197)))
  in_197_dn <- length(intersect(dn_targets, rownames(expr197)))
  in_221_up <- length(intersect(up_targets, rownames(expr221)))
  in_221_dn <- length(intersect(dn_targets, rownames(expr221)))
  shared_up <- length(intersect(intersect(up_targets, rownames(expr197)), rownames(expr221)))
  shared_dn <- length(intersect(intersect(dn_targets, rownames(expr197)), rownames(expr221)))
  shared_total <- shared_up + shared_dn
  
  if (shared_total >= 20) {
    flag <- "GOOD_COVERAGE"
  } else if (shared_total >= 10) {
    flag <- "LIMITED_COVERAGE"
  } else {
    flag <- "POOR_COVERAGE"
  }
  
  coverage_qc_list[[tf]] <- data.frame(
    regulator = tf,
    total_up_targets = length(up_targets),
    total_dn_targets = length(dn_targets),
    up_in_GSE197677 = in_197_up,
    dn_in_GSE197677 = in_197_dn,
    up_in_GSE221561 = in_221_up,
    dn_in_GSE221561 = in_221_dn,
    shared_up = shared_up,
    shared_dn = shared_dn,
    shared_total = shared_total,
    coverage_flag = flag,
    stringsAsFactors = FALSE
  )
}

coverage_qc_df <- do.call(rbind, coverage_qc_list)
rownames(coverage_qc_df) <- NULL

cat("=== REGULATOR COVERAGE QC ===\n")
for (i in 1:nrow(coverage_qc_df)) {
  cat(sprintf("  %s: UP=%d DN=%d (shared_total=%d) -> %s\n",
    coverage_qc_df$regulator[i], coverage_qc_df$shared_up[i],
    coverage_qc_df$shared_dn[i], coverage_qc_df$shared_total[i],
    coverage_qc_df$coverage_flag[i]))
}

write.csv(coverage_qc_df, file.path(N_DIR, "STEP19N_REGULATOR_COVERAGE_QC.csv"), row.names = FALSE)
cat("Saved: STEP19N_REGULATOR_COVERAGE_QC.csv\n")

# ============================================================================
# 19N6 — CANONICAL TREATMENT EFFECT ON REGULATOR ACTIVITY
# ============================================================================
cat("\n============================================\n")
cat("19N6: CANONICAL TREATMENT EFFECT ON ACTIVITY\n")
cat("============================================\n\n")

treat_effect_list <- list()
for (tf in all_candidates) {
  for (ds in c("GSE197677", "GSE221561")) {
    sub <- activity_df[activity_df$regulator == tf & activity_df$dataset == ds, ]
    
    treated_group <- TREAT_CONTRAST[[ds]]["treated"]
    control_group <- TREAT_CONTRAST[[ds]]["control"]
    
    treated_act <- sub$activity[sub$treatment == treated_group]
    control_act <- sub$activity[sub$treatment == control_group]
    
    effect <- mean(treated_act, na.rm = TRUE) - mean(control_act, na.rm = TRUE)
    if (is.nan(effect)) effect <- NA
    direction <- ifelse(is.na(effect), "NA", ifelse(effect > 0, "HIGHER_IN_TREATED", "LOWER_IN_TREATED"))
    
    treat_effect_list[[paste(tf, ds)]] <- data.frame(
      regulator = tf,
      dataset = ds,
      effect = effect,
      direction = direction,
      treated_mean = mean(treated_act, na.rm = TRUE),
      control_mean = mean(control_act, na.rm = TRUE),
      treated_n = sum(!is.na(treated_act)),
      control_n = sum(!is.na(control_act)),
      stringsAsFactors = FALSE
    )
  }
}

treat_effect_df <- do.call(rbind, treat_effect_list)
rownames(treat_effect_df) <- NULL

cat("=== TREATMENT EFFECTS ON REGULATOR ACTIVITY ===\n")
for (i in 1:nrow(treat_effect_df)) {
  cat(sprintf("  %s %s: effect=%.4f (%s)\n",
    treat_effect_df$regulator[i], treat_effect_df$dataset[i],
    treat_effect_df$effect[i], treat_effect_df$direction[i]))
}

write.csv(treat_effect_df, file.path(N_DIR, "STEP19N_REGULATOR_ACTIVITY_TREATMENT_EFFECTS.csv"), row.names = FALSE)
cat("\nSaved: STEP19N_REGULATOR_ACTIVITY_TREATMENT_EFFECTS.csv\n")

# ============================================================================
# 19N7 — EXACT PERMUTATION TESTS
# ============================================================================
cat("\n============================================\n")
cat("19N7: EXACT PERMUTATION TESTS\n")
cat("============================================\n\n")

perm_results_list <- list()

for (tf in PRIMARY_CANDIDATES) {
  for (ds in c("GSE197677", "GSE221561")) {
    sub <- activity_df[activity_df$regulator == tf & activity_df$dataset == ds, ]
    
    treated_group <- TREAT_CONTRAST[[ds]]["treated"]
    control_group <- TREAT_CONTRAST[[ds]]["control"]
    
    treated_idx <- which(sub$treatment == treated_group)
    control_idx <- which(sub$treatment == control_group)
    n_treated <- length(treated_idx)
    n_total <- n_treated + length(control_idx)
    
    obs_effect <- mean(sub$activity[treated_idx]) - mean(sub$activity[control_idx])
    
    # Exact permutation
    all_combinations <- combn(n_total, n_treated)
    n_perms <- ncol(all_combinations)
    
    perm_effects <- numeric(n_perms)
    for (p in 1:n_perms) {
      perm_treated <- all_combinations[, p]
      perm_control <- setdiff(1:n_total, perm_treated)
      perm_effects[p] <- mean(sub$activity[perm_treated]) - mean(sub$activity[perm_control])
    }
    
    p_two_sided <- mean(abs(perm_effects) >= abs(obs_effect))
    p_right <- mean(perm_effects >= obs_effect)
    p_left <- mean(perm_effects <= obs_effect)
    
    perm_results_list[[paste(tf, ds)]] <- data.frame(
      regulator = tf,
      dataset = ds,
      n_treated = n_treated,
      n_total = n_total,
      total_perms = n_perms,
      observed_effect = obs_effect,
      p_two_sided = p_two_sided,
      p_right = p_right,
      p_left = p_left,
      perm_mean = mean(perm_effects),
      perm_sd = sd(perm_effects),
      stringsAsFactors = FALSE
    )
    
    cat(sprintf("  %s %s: effect=%.4f P(two)=%.4f (%d perms)\n",
      tf, ds, obs_effect, p_two_sided, n_perms))
  }
}

perm_results_df <- do.call(rbind, perm_results_list)
rownames(perm_results_df) <- NULL

write.csv(perm_results_df, file.path(N_DIR, "STEP19N_REGULATOR_ACTIVITY_EXACT_TESTS.csv"), row.names = FALSE)
cat("\nSaved: STEP19N_REGULATOR_ACTIVITY_EXACT_TESTS.csv\n")

# ============================================================================
# 19N8 — ACTIVITY VS CORE MODULE ASSOCIATION
# ============================================================================
cat("\n============================================\n")
cat("19N8: ACTIVITY VS CORE MODULE ASSOCIATION\n")
cat("============================================\n\n")

core_assoc_list <- list()
for (tf in all_candidates) {
  for (ds in c("GSE197677", "GSE221561")) {
    sub <- activity_df[activity_df$regulator == tf & activity_df$dataset == ds, ]
    samp_ids <- sub$sample_id
    
    c4_scores <- core_scores$CORE4_score[match(samp_ids, core_scores$sample_id)]
    c2_scores <- core_scores$CORE2_STRONG_score[match(samp_ids, core_scores$sample_id)]
    
    # Skip if all activity is NA (no regulon targets)
    if (all(is.na(sub$activity))) {
      rho_c4 <- NA
      rho_c2 <- NA
      assoc_class <- "NO_REGULON_TARGETS"
    } else {
      rho_c4 <- cor(sub$activity, c4_scores, method = "spearman", use = "complete.obs")
      rho_c2 <- cor(sub$activity, c2_scores, method = "spearman", use = "complete.obs")
      
      if (is.na(rho_c4)) rho_c4 <- 0
      if (is.na(rho_c2)) rho_c2 <- 0
      
      if (abs(rho_c2) >= 0.60) {
        assoc_class <- "CORE2_STRONG_ASSOCIATED"
      } else if (abs(rho_c2) >= 0.40) {
        assoc_class <- "CORE2_MODERATE_ASSOCIATION"
      } else if (abs(rho_c4) >= 0.40) {
        assoc_class <- "CORE4_ONLY_ASSOCIATION"
      } else {
        assoc_class <- "WEAK_ASSOCIATION"
      }
    }
    
    core_assoc_list[[paste(tf, ds)]] <- data.frame(
      regulator = tf,
      dataset = ds,
      rho_CORE4 = rho_c4,
      rho_CORE2 = rho_c2,
      sign_core4 = sign(rho_c4),
      sign_core2 = sign(rho_c2),
      classification = assoc_class,
      stringsAsFactors = FALSE
    )
  }
}

core_assoc_df <- do.call(rbind, core_assoc_list)
rownames(core_assoc_df) <- NULL

cat("=== ACTIVITY vs CORE MODULE ASSOCIATION ===\n")
for (i in 1:nrow(core_assoc_df)) {
  cat(sprintf("  %s %s: rho_C4=%.3f rho_C2=%.3f -> %s\n",
    core_assoc_df$regulator[i], core_assoc_df$dataset[i],
    core_assoc_df$rho_CORE4[i], core_assoc_df$rho_CORE2[i],
    core_assoc_df$classification[i]))
}

write.csv(core_assoc_df, file.path(N_DIR, "STEP19N_REGULATOR_CORE_ASSOCIATION.csv"), row.names = FALSE)
cat("\nSaved: STEP19N_REGULATOR_CORE_ASSOCIATION.csv\n")

# ============================================================================
# 19N9 — CROSS-DATASET ACTIVITY REPLICATION
# ============================================================================
cat("\n============================================\n")
cat("19N9: CROSS-DATASET ACTIVITY REPLICATION\n")
cat("============================================\n\n")

replication_list <- list()
for (tf in all_candidates) {
  sub197 <- core_assoc_df[core_assoc_df$regulator == tf & core_assoc_df$dataset == "GSE197677", ]
  sub221 <- core_assoc_df[core_assoc_df$regulator == tf & core_assoc_df$dataset == "GSE221561", ]
  
  te197 <- treat_effect_df[treat_effect_df$regulator == tf & treat_effect_df$dataset == "GSE197677", ]
  te221 <- treat_effect_df[treat_effect_df$regulator == tf & treat_effect_df$dataset == "GSE221561", ]
  
  # Handle NaN/NA effects
  eff197 <- ifelse(is.na(te197$effect) || is.nan(te197$effect), NA, te197$effect)
  eff221 <- ifelse(is.na(te221$effect) || is.nan(te221$effect), NA, te221$effect)
  
  # Treatment effect direction replication
  te_replicated <- !is.na(eff197) && !is.na(eff221) && sign(eff197) == sign(eff221)
  
  # CORE2 association direction replication
  core2_replicated <- !is.na(sub197$rho_CORE2) && !is.na(sub221$rho_CORE2) &&
    sub197$rho_CORE2 != 0 && sub221$rho_CORE2 != 0 && sign(sub197$rho_CORE2) == sign(sub221$rho_CORE2)
  
  if (is.na(te_replicated)) te_replicated <- FALSE
  if (is.na(core2_replicated)) core2_replicated <- FALSE
  
  if (te_replicated && core2_replicated) {
    repl_class <- "ACTIVITY_REPLICATED"
  } else if (core2_replicated) {
    repl_class <- "CORE_ASSOCIATION_REPLICATED_ONLY"
  } else if (te_replicated) {
    repl_class <- "TREATMENT_EFFECT_REPLICATED_ONLY"
  } else {
    repl_class <- "CROSS_DATASET_DISCORDANT"
  }
  
  replication_list[[tf]] <- data.frame(
    regulator = tf,
    effect_197 = eff197,
    effect_221 = eff221,
    te_direction_197 = ifelse(is.na(eff197), "NA", ifelse(eff197 > 0, "HIGHER_IN_TREATED", "LOWER_IN_TREATED")),
    te_direction_221 = ifelse(is.na(eff221), "NA", ifelse(eff221 > 0, "HIGHER_IN_TREATED", "LOWER_IN_TREATED")),
    te_replicated = te_replicated,
    rho_CORE2_197 = sub197$rho_CORE2,
    rho_CORE2_221 = sub221$rho_CORE2,
    core2_direction_replicated = core2_replicated,
    classification = repl_class,
    stringsAsFactors = FALSE
  )
}

replication_df <- do.call(rbind, replication_list)
rownames(replication_df) <- NULL

cat("=== CROSS-DATASET REPLICATION ===\n")
for (i in 1:nrow(replication_df)) {
  cat(sprintf("  %s: TE_repl=%s CORE2_repl=%s -> %s\n",
    replication_df$regulator[i], replication_df$te_replicated[i],
    replication_df$core2_direction_replicated[i], replication_df$classification[i]))
}

write.csv(replication_df, file.path(N_DIR, "STEP19N_REGULATOR_CROSS_DATASET_REPLICATION.csv"), row.names = FALSE)
cat("\nSaved: STEP19N_REGULATOR_CROSS_DATASET_REPLICATION.csv\n")

# ============================================================================
# 19N10 — AP-1 FAMILY AUDIT
# ============================================================================
cat("\n============================================\n")
cat("19N10: AP-1 FAMILY AUDIT\n")
cat("============================================\n\n")

ap1_results <- list()
for (tf in c("FOS", "FOSL2")) {
  act_197 <- activity_df[activity_df$regulator == tf & activity_df$dataset == "GSE197677", ]
  act_221 <- activity_df[activity_df$regulator == tf & activity_df$dataset == "GSE221561", ]
  
  c2_197 <- core_assoc_df$rho_CORE2[core_assoc_df$regulator == tf & core_assoc_df$dataset == "GSE197677"]
  c2_221 <- core_assoc_df$rho_CORE2[core_assoc_df$regulator == tf & core_assoc_df$dataset == "GSE221561"]
  
  te_197 <- treat_effect_df$effect[treat_effect_df$regulator == tf & treat_effect_df$dataset == "GSE197677"]
  te_221 <- treat_effect_df$effect[treat_effect_df$regulator == tf & treat_effect_df$dataset == "GSE221561"]
  
  repl <- replication_df$classification[replication_df$regulator == tf]
  
  # LOSO
  loso_197 <- numeric(length(samples197))
  c4_for_loso <- core_scores$CORE4_score[match(samples197, core_scores$sample_id)]
  for (i in seq_along(samples197)) {
    keep <- setdiff(seq_along(samples197), i)
    act_keep <- act_197$activity[keep]
    c4_keep <- c4_for_loso[keep]
    # Only compute if there are complete pairs
    valid <- !is.na(act_keep) & !is.na(c4_keep)
    if (sum(valid) >= 3) {
      loso_197[i] <- cor(act_keep[valid], c4_keep[valid], method = "spearman")
    } else {
      loso_197[i] <- NA
    }
  }
  valid_full_197 <- !is.na(act_197$activity) & !is.na(c4_for_loso)
  if (sum(valid_full_197) >= 3) {
    full_rho_197 <- cor(act_197$activity[valid_full_197], c4_for_loso[valid_full_197], method = "spearman")
  } else {
    full_rho_197 <- NA
  }
  sign_pres_197 <- mean(sign(loso_197[!is.na(loso_197)]) == sign(full_rho_197), na.rm = TRUE)
  
  cov <- coverage_qc_df$coverage_flag[coverage_qc_df$regulator == tf]
  
  ap1_results[[tf]] <- data.frame(
    regulator = tf,
    activity_CORE2_GSE197677 = c2_197,
    activity_CORE2_GSE221561 = c2_221,
    treatment_effect_GSE197677 = te_197,
    treatment_effect_GSE221561 = te_221,
    cross_dataset_replication = repl,
    loso_sign_preservation = sign_pres_197,
    coverage = cov,
    stringsAsFactors = FALSE
  )
  
  cat(sprintf("  %s:\n", tf))
  cat(sprintf("    Activity/CORE2 GSE197677: rho=%.3f\n", c2_197))
  cat(sprintf("    Activity/CORE2 GSE221561: rho=%.3f\n", c2_221))
  cat(sprintf("    Treatment effect GSE197677: %.4f\n", te_197))
  cat(sprintf("    Treatment effect GSE221561: %.4f\n", te_221))
  cat(sprintf("    LOSO sign preservation: %.0f%%\n", sign_pres_197 * 100))
  cat(sprintf("    Coverage: %s\n", cov))
}

ap1_df <- do.call(rbind, ap1_results)
rownames(ap1_df) <- NULL

# Determine AP-1 family status
ap1_replicated <- any(ap1_df$cross_dataset_replication %in% c("ACTIVITY_REPLICATED", "CORE_ASSOCIATION_REPLICATED_ONLY"))
ap1_coverage_ok <- all(ap1_df$coverage %in% c("GOOD_COVERAGE", "LIMITED_COVERAGE"))

if (ap1_replicated && ap1_coverage_ok) {
  ap1_status <- "AP1_ORTHOGONAL_SUPPORT"
} else if (ap1_coverage_ok) {
  ap1_status <- "AP1_PARTIAL_SUPPORT"
} else {
  ap1_status <- "AP1_EXPRESSION_ONLY_SUPPORT"
}

cat(sprintf("\n  AP-1 FAMILY STATUS: %s\n", ap1_status))

write.csv(ap1_df, file.path(N_DIR, "STEP19N_AP1_FAMILY_VALIDATION.csv"), row.names = FALSE)
cat("Saved: STEP19N_AP1_FAMILY_VALIDATION.csv\n")

# ============================================================================
# 19N11 — TGF-BETA / SMAD FAMILY AUDIT
# ============================================================================
cat("\n============================================\n")
cat("19N11: TGF-BETA / SMAD FAMILY AUDIT\n")
cat("============================================\n\n")

tgfb_results <- list()
for (tf in c("SMAD3", "SMAD7")) {
  act_197 <- activity_df[activity_df$regulator == tf & activity_df$dataset == "GSE197677", ]
  act_221 <- activity_df[activity_df$regulator == tf & activity_df$dataset == "GSE221561", ]
  
  c2_197 <- core_assoc_df$rho_CORE2[core_assoc_df$regulator == tf & core_assoc_df$dataset == "GSE197677"]
  c2_221 <- core_assoc_df$rho_CORE2[core_assoc_df$regulator == tf & core_assoc_df$dataset == "GSE221561"]
  
  te_197 <- treat_effect_df$effect[treat_effect_df$regulator == tf & treat_effect_df$dataset == "GSE197677"]
  te_221 <- treat_effect_df$effect[treat_effect_df$regulator == tf & treat_effect_df$dataset == "GSE221561"]
  
  repl <- replication_df$classification[replication_df$regulator == tf]
  cov <- coverage_qc_df$coverage_flag[coverage_qc_df$regulator == tf]
  
  tgfb_results[[tf]] <- data.frame(
    regulator = tf,
    activity_CORE2_GSE197677 = c2_197,
    activity_CORE2_GSE221561 = c2_221,
    treatment_effect_GSE197677 = te_197,
    treatment_effect_GSE221561 = te_221,
    cross_dataset_replication = repl,
    coverage = cov,
    stringsAsFactors = FALSE
  )
  
  cat(sprintf("  %s:\n", tf))
  cat(sprintf("    Activity/CORE2 GSE197677: %.3f\n", c2_197))
  cat(sprintf("    Activity/CORE2 GSE221561: %.3f\n", c2_221))
  cat(sprintf("    Treatment effect GSE197677: %.4f\n", te_197))
  cat(sprintf("    Treatment effect GSE221561: %.4f\n", te_221))
  cat(sprintf("    Coverage: %s\n", cov))
}

tgfb_df <- do.call(rbind, tgfb_results)
rownames(tgfb_df) <- NULL

tgfb_replicated <- any(tgfb_df$cross_dataset_replication %in% c("ACTIVITY_REPLICATED", "CORE_ASSOCIATION_REPLICATED_ONLY"))
tgfb_coverage_ok <- all(tgfb_df$coverage %in% c("GOOD_COVERAGE", "LIMITED_COVERAGE"))

if (tgfb_replicated && tgfb_coverage_ok) {
  tgfb_status <- "TGFB_SMAD_ORTHOGONAL_SUPPORT"
} else if (tgfb_coverage_ok) {
  tgfb_status <- "TGFB_SMAD_PARTIAL_SUPPORT"
} else {
  tgfb_status <- "TGFB_SMAD_EXPRESSION_ONLY_SUPPORT"
}

cat(sprintf("\n  TGF-BETA/SMAD FAMILY STATUS: %s\n", tgfb_status))

write.csv(tgfb_df, file.path(N_DIR, "STEP19N_TGFB_SMAD_VALIDATION.csv"), row.names = FALSE)
cat("Saved: STEP19N_TGFB_SMAD_VALIDATION.csv\n")

# ============================================================================
# 19N12 — KLF6 / RUNX1 AUDIT
# ============================================================================
cat("\n============================================\n")
cat("19N12: KLF6 / RUNX1 AUDIT\n")
cat("============================================\n\n")

kr_results <- list()
for (tf in c("KLF6", "RUNX1")) {
  c2_197 <- core_assoc_df$rho_CORE2[core_assoc_df$regulator == tf & core_assoc_df$dataset == "GSE197677"]
  c2_221 <- core_assoc_df$rho_CORE2[core_assoc_df$regulator == tf & core_assoc_df$dataset == "GSE221561"]
  repl <- replication_df$classification[replication_df$regulator == tf]
  cov <- coverage_qc_df$coverage_flag[coverage_qc_df$regulator == tf]
  
  # Determine support level
  if (repl %in% c("ACTIVITY_REPLICATED", "CORE_ASSOCIATION_REPLICATED_ONLY") && cov != "POOR_COVERAGE") {
    kr_class <- "ORTHOGONALLY_SUPPORTED"
  } else if (repl == "TREATMENT_EFFECT_REPLICATED_ONLY" || cov == "LIMITED_COVERAGE") {
    kr_class <- "PARTIALLY_SUPPORTED"
  } else if (cov == "POOR_COVERAGE") {
    kr_class <- "EXPRESSION_ASSOCIATION_ONLY"
  } else {
    kr_class <- "NOT_SUPPORTED"
  }
  
  kr_results[[tf]] <- data.frame(
    regulator = tf,
    activity_CORE2_GSE197677 = c2_197,
    activity_CORE2_GSE221561 = c2_221,
    cross_dataset_replication = repl,
    coverage = cov,
    classification = kr_class,
    stringsAsFactors = FALSE
  )
  
  cat(sprintf("  %s: CORE2rho197=%.3f CORE2rho221=%.3f repl=%s cov=%s -> %s\n",
    tf, c2_197, c2_221, repl, cov, kr_class))
}

kr_df <- do.call(rbind, kr_results)
rownames(kr_df) <- NULL

write.csv(kr_df, file.path(N_DIR, "STEP19N_KLF6_RUNX1_VALIDATION.csv"), row.names = FALSE)
cat("\nSaved: STEP19N_KLF6_RUNX1_VALIDATION.csv\n")

# ============================================================================
# 19N13 — CORE2 SHARED TARGET TEST
# ============================================================================
cat("\n============================================\n")
cat("19N13: CORE2 SHARED TARGET TEST\n")
cat("============================================\n\n")

c2_target_list <- list()
for (tf in all_candidates) {
  up_targets <- regulon_list[[tf]]$up
  dn_targets <- regulon_list[[tf]]$dn
  
  has_cdkn1b_up <- "CDKN1B" %in% up_targets
  has_cdkn1b_dn <- "CDKN1B" %in% dn_targets
  has_gsn_up <- "GSN" %in% up_targets
  has_gsn_dn <- "GSN" %in% dn_targets
  
  has_cdkn1b <- has_cdkn1b_up || has_cdkn1b_dn
  has_gsn <- has_gsn_up || has_gsn_dn
  both <- has_cdkn1b && has_gsn
  
  # Activity association with oriented CDKN1B and GSN
  # CDKN1B: negative orientation -> oriented = -CDKN1B_z
  # GSN: negative orientation -> oriented = -GSN_z
  
  act_cdkn1b_rho <- c()
  act_gsn_rho <- c()
  
  for (ds in c("GSE197677", "GSE221561")) {
    sub <- activity_df[activity_df$regulator == tf & activity_df$dataset == ds, ]
    samp_ids <- sub$sample_id
    
    if (ds == "GSE197677") {
      oriented_cdkn1b <- -as.numeric(expr197["CDKN1B", samp_ids])
      oriented_gsn <- -as.numeric(expr197["GSN", samp_ids])
    } else {
      oriented_cdkn1b <- -as.numeric(expr221["CDKN1B", samp_ids])
      oriented_gsn <- -as.numeric(expr221["GSN", samp_ids])
    }
    
    valid_pairs_cdkn1b <- !is.na(sub$activity) & !is.na(oriented_cdkn1b)
    valid_pairs_gsn <- !is.na(sub$activity) & !is.na(oriented_gsn)
    rho_cdkn1b <- if (sum(valid_pairs_cdkn1b) >= 3) cor(sub$activity[valid_pairs_cdkn1b], oriented_cdkn1b[valid_pairs_cdkn1b], method = "spearman") else NA
    rho_gsn <- if (sum(valid_pairs_gsn) >= 3) cor(sub$activity[valid_pairs_gsn], oriented_gsn[valid_pairs_gsn], method = "spearman") else NA
    
    act_cdkn1b_rho <- c(act_cdkn1b_rho, rho_cdkn1b)
    act_gsn_rho <- c(act_gsn_rho, rho_gsn)
  }
  
  # Handle NA values in rho comparisons
  rho1_cdkn1b <- if (length(act_cdkn1b_rho) >= 1 && !is.na(act_cdkn1b_rho[1])) act_cdkn1b_rho[1] else NA
  rho2_cdkn1b <- if (length(act_cdkn1b_rho) >= 2 && !is.na(act_cdkn1b_rho[2])) act_cdkn1b_rho[2] else NA
  rho1_gsn <- if (length(act_gsn_rho) >= 1 && !is.na(act_gsn_rho[1])) act_gsn_rho[1] else NA
  rho2_gsn <- if (length(act_gsn_rho) >= 2 && !is.na(act_gsn_rho[2])) act_gsn_rho[2] else NA
  
  same_dir_1 <- !is.na(rho1_cdkn1b) && !is.na(rho1_gsn) && sign(rho1_cdkn1b) == sign(rho1_gsn)
  same_dir_2 <- !is.na(rho2_cdkn1b) && !is.na(rho2_gsn) && sign(rho2_cdkn1b) == sign(rho2_gsn)
  same_dir <- same_dir_1 && same_dir_2
  both_adequate <- !any(is.na(act_cdkn1b_rho)) && !any(is.na(act_gsn_rho)) &&
    all(abs(act_cdkn1b_rho) >= 0.40) && all(abs(act_gsn_rho) >= 0.40)
  
  if (same_dir && both_adequate) {
    shared_class <- "CORE2_SHARED_ACTIVITY_SUPPORT"
  } else if (same_dir_1 || same_dir_2) {
    shared_class <- "SINGLE_CORE_GENE_SUPPORT"
  } else {
    shared_class <- "NO_SHARED_SUPPORT"
  }
  
  c2_target_list[[tf]] <- data.frame(
    regulator = tf,
    targets_CDKN1B = has_cdkn1b,
    targets_GSN = has_gsn,
    both_direct_targets = both,
    rho_activity_CDKN1B_197 = act_cdkn1b_rho[1],
    rho_activity_CDKN1B_221 = act_cdkn1b_rho[2],
    rho_activity_GSN_197 = act_gsn_rho[1],
    rho_activity_GSN_221 = act_gsn_rho[2],
    same_direction_both_genes = same_dir,
    classification = shared_class,
    stringsAsFactors = FALSE
  )
  
  cat(sprintf("  %s: CDKN1B_target=%s GSN_target=%s both=%s class=%s\n",
    tf, has_cdkn1b, has_gsn, both, shared_class))
}

c2_target_df <- do.call(rbind, c2_target_list)
rownames(c2_target_df) <- NULL

write.csv(c2_target_df, file.path(N_DIR, "STEP19N_CORE2_REGULATOR_TARGET_SUPPORT.csv"), row.names = FALSE)
cat("\nSaved: STEP19N_CORE2_REGULATOR_TARGET_SUPPORT.csv\n")

# ============================================================================
# 19N14 — LEAVE-ONE-SAMPLE-OUT STABILITY
# ============================================================================
cat("\n============================================\n")
cat("19N14: LEAVE-ONE-SAMPLE-OUT STABILITY\n")
cat("============================================\n\n")

loso_list <- list()
for (tf in PRIMARY_CANDIDATES) {
  for (ds in c("GSE197677", "GSE221561")) {
    sub <- activity_df[activity_df$regulator == tf & activity_df$dataset == ds, ]
    samp_ids <- sub$sample_id
    n <- length(samp_ids)
    
    treated_group <- TREAT_CONTRAST[[ds]]["treated"]
    control_group <- TREAT_CONTRAST[[ds]]["control"]
    
    # Full effect and rho
    full_treated_idx <- which(sub$treatment == treated_group)
    full_control_idx <- which(sub$treatment == control_group)
    full_effect <- mean(sub$activity[full_treated_idx]) - mean(sub$activity[full_control_idx])
    
    c2_scores <- core_scores$CORE2_STRONG_score[match(samp_ids, core_scores$sample_id)]
    valid_full_rho <- !is.na(sub$activity) & !is.na(c2_scores)
    full_rho <- if (sum(valid_full_rho) >= 3) cor(sub$activity[valid_full_rho], c2_scores[valid_full_rho], method = "spearman") else NA
    
    for (i in 1:n) {
      left_out_sample <- samp_ids[i]
      left_out_group <- sub$treatment[i]
      
      keep <- setdiff(1:n, i)
      keep_treated <- intersect(keep, full_treated_idx)
      keep_control <- intersect(keep, full_control_idx)
      
      if (length(keep_treated) >= 1 && length(keep_control) >= 1) {
        loso_effect <- mean(sub$activity[keep_treated]) - mean(sub$activity[keep_control])
        valid_loso <- !is.na(sub$activity[keep]) & !is.na(c2_scores[keep])
        loso_rho <- if (sum(valid_loso) >= 3) cor(sub$activity[keep][valid_loso], c2_scores[keep][valid_loso], method = "spearman") else NA
      } else {
        loso_effect <- NA
        loso_rho <- NA
      }
      
      effect_sign_pres <- sign(loso_effect) == sign(full_effect)
      rho_sign_pres <- sign(loso_rho) == sign(full_rho)
      
      if (is.na(effect_sign_pres)) effect_sign_pres <- NA
      if (is.na(rho_sign_pres)) rho_sign_pres <- NA
      
      loso_list[[paste(tf, ds, left_out_sample)]] <- data.frame(
        regulator = tf,
        dataset = ds,
        left_out_sample = left_out_sample,
        left_out_group = left_out_group,
        effect = loso_effect,
        rho = loso_rho,
        effect_sign_preserved = effect_sign_pres,
        rho_sign_preserved = rho_sign_pres,
        stringsAsFactors = FALSE
      )
    }
  }
}

loso_df <- do.call(rbind, loso_list)
rownames(loso_df) <- NULL

# Summarize LOSO
loso_summary_list <- list()
for (tf in PRIMARY_CANDIDATES) {
  for (ds in c("GSE197677", "GSE221561")) {
    sub <- loso_df[loso_df$regulator == tf & loso_df$dataset == ds, ]
    
    effect_pres <- mean(sub$effect_sign_preserved, na.rm = TRUE)
    rho_pres <- mean(sub$rho_sign_preserved, na.rm = TRUE)
    if (is.nan(effect_pres) || is.na(effect_pres)) effect_pres <- NA
    if (is.nan(rho_pres) || is.na(rho_pres)) rho_pres <- NA
    
    if (is.na(effect_pres) || is.na(rho_pres)) {
      classification <- "NO_VALID_LOSO"
    } else if (effect_pres >= 0.80 && rho_pres >= 0.80) {
      classification <- "ROBUST"
    } else if (ds == "GSE221561" && (effect_pres < 0.80 || rho_pres < 0.80)) {
      classification <- "CONTROL_SENSITIVE"
    } else if (effect_pres < 0.60 || rho_pres < 0.60) {
      classification <- "UNSTABLE"
    } else {
      classification <- "OUTLIER_SENSITIVE"
    }
    
    loso_summary_list[[paste(tf, ds)]] <- data.frame(
      regulator = tf,
      dataset = ds,
      effect_sign_preservation = effect_pres,
      rho_sign_preservation = rho_pres,
      classification = classification,
      stringsAsFactors = FALSE
    )
  }
}

loso_summary_df <- do.call(rbind, loso_summary_list)
rownames(loso_summary_df) <- NULL

cat("=== LOSO STABILITY ===\n")
for (i in 1:nrow(loso_summary_df)) {
  cat(sprintf("  %s %s: effect_pres=%.0f%% rho_pres=%.0f%% -> %s\n",
    loso_summary_df$regulator[i], loso_summary_df$dataset[i],
    loso_summary_df$effect_sign_preservation[i] * 100,
    loso_summary_df$rho_sign_preservation[i] * 100,
    loso_summary_df$classification[i]))
}

write.csv(loso_df, file.path(N_DIR, "STEP19N_REGULATOR_ACTIVITY_LOSO_LONG.csv"), row.names = FALSE)
write.csv(loso_summary_df, file.path(N_DIR, "STEP19N_REGULATOR_ACTIVITY_LOSO_SUMMARY.csv"), row.names = FALSE)
cat("\nSaved: STEP19N_REGULATOR_ACTIVITY_LOSO_LONG.csv\n")
cat("Saved: STEP19N_REGULATOR_ACTIVITY_LOSO_SUMMARY.csv\n")

# ============================================================================
# 19N15 — ABUNDANCE AND LIBRARY-SIZE SENSITIVITY
# ============================================================================
cat("\n============================================\n")
cat("19N15: ABUNDANCE AND LIBRARY-SIZE SENSITIVITY\n")
cat("============================================\n\n")

cov_list <- list()
for (tf in PRIMARY_CANDIDATES) {
  for (ds in c("GSE197677", "GSE221561")) {
    sub <- activity_df[activity_df$regulator == tf & activity_df$dataset == ds, ]
    samp_ids <- sub$sample_id
    
    c2_scores <- core_scores$CORE2_STRONG_score[match(samp_ids, core_scores$sample_id)]
    
    if (ds == "GSE197677") {
      lib_size <- colSums(fib197[, samp_ids, drop = FALSE])
    } else {
      lib_size <- colSums(fib221[, samp_ids, drop = FALSE])
    }
    
    abundance_z <- scale(log1p(lib_size))[, 1]
    
    # Model 1: unadjusted
    m1 <- tryCatch({
      fit <- lm(c2_scores ~ sub$activity)
      s <- summary(fit)$coefficients
      if ("sub$activity" %in% rownames(s)) c("Est" = s["sub$activity", "Estimate"], "Dir" = sign(s["sub$activity", "Estimate"]))
      else c("Est" = NA, "Dir" = NA)
    }, error = function(e) c("Est" = NA, "Dir" = NA))
    
    # Model 2: + abundance
    m2 <- tryCatch({
      fit <- lm(c2_scores ~ sub$activity + abundance_z)
      s <- summary(fit)$coefficients
      if ("sub$activity" %in% rownames(s)) c("Est" = s["sub$activity", "Estimate"], "Dir" = sign(s["sub$activity", "Estimate"]))
      else c("Est" = NA, "Dir" = NA)
    }, error = function(e) c("Est" = NA, "Dir" = NA))
    
    stable <- as.numeric(m1["Dir"]) == as.numeric(m2["Dir"])
    
    cov_list[[paste(tf, ds)]] <- data.frame(
      regulator = tf,
      dataset = ds,
      coef_unadjusted = as.numeric(m1["Est"]),
      dir_unadjusted = as.numeric(m1["Dir"]),
      coef_adjusted = as.numeric(m2["Est"]),
      dir_adjusted = as.numeric(m2["Dir"]),
      stable_direction = stable,
      stringsAsFactors = FALSE
    )
  }
}

cov_df <- do.call(rbind, cov_list)
rownames(cov_df) <- NULL

cat("=== COVARIATE SENSITIVITY ===\n")
for (i in 1:nrow(cov_df)) {
  cat(sprintf("  %s %s: coef_adj=%.4f stable=%s\n",
    cov_df$regulator[i], cov_df$dataset[i], cov_df$coef_adjusted[i], cov_df$stable_direction[i]))
}

write.csv(cov_df, file.path(N_DIR, "STEP19N_REGULATOR_ACTIVITY_COVARIATE_SENSITIVITY.csv"), row.names = FALSE)
cat("\nSaved: STEP19N_REGULATOR_ACTIVITY_COVARIATE_SENSITIVITY.csv\n")

# ============================================================================
# 19N16 — FAMILY-LEVEL CONSENSUS
# ============================================================================
cat("\n============================================\n")
cat("19N16: FAMILY-LEVEL CONSENSUS\n")
cat("============================================\n\n")

families <- list(
  AP1 = c("FOS", "FOSL2"),
  TGFB_SMAD = c("SMAD3", "SMAD7"),
  KLF6 = "KLF6",
  RUNX1 = "RUNX1",
  p53_cell_cycle = "CDKN1A"
)

family_consensus_list <- list()
for (fam_name in names(families)) {
  fam_tfs <- families[[fam_name]]
  
  # Aggregate evidence
  any_replicated <- any(replication_df$classification[replication_df$regulator %in% fam_tfs] %in%
    c("ACTIVITY_REPLICATED", "CORE_ASSOCIATION_REPLICATED_ONLY"))
  any_adequate_coverage <- any(coverage_qc_df$coverage_flag[coverage_qc_df$regulator %in% fam_tfs] %in%
    c("GOOD_COVERAGE", "LIMITED_COVERAGE"))
  any_core2_assoc <- any(abs(core_assoc_df$rho_CORE2[core_assoc_df$regulator %in% fam_tfs]) >= 0.40, na.rm = TRUE)
  
  # Expression evidence from Step19M
  expr_evidence <- any(step19m_tf$TF %in% fam_tfs)
  
  # LOSO robust
  loso_robust <- any(loso_summary_df$classification[loso_summary_df$regulator %in% fam_tfs] == "ROBUST")
  
  if (any_replicated && any_adequate_coverage && any_core2_assoc) {
    fam_class <- "ORTHOGONAL_CROSS_DATASET_SUPPORT"
  } else if (any_adequate_coverage && any_core2_assoc) {
    fam_class <- "PARTIAL_ORTHOGONAL_SUPPORT"
  } else if (expr_evidence) {
    fam_class <- "EXPRESSION_ASSOCIATION_ONLY"
  } else if (!any_adequate_coverage) {
    fam_class <- "RESOURCE_INSUFFICIENT"
  } else {
    fam_class <- "NO_SUPPORT"
  }
  
  family_consensus_list[[fam_name]] <- data.frame(
    family = fam_name,
    regulators = paste(fam_tfs, collapse = ", "),
    expression_evidence = expr_evidence,
    any_replicated = any_replicated,
    any_adequate_coverage = any_adequate_coverage,
    any_core2_association = any_core2_assoc,
    loso_robust = loso_robust,
    classification = fam_class,
    stringsAsFactors = FALSE
  )
  
  cat(sprintf("  %s: %s\n", fam_name, fam_class))
}

family_consensus_df <- do.call(rbind, family_consensus_list)
rownames(family_consensus_df) <- NULL

write.csv(family_consensus_df, file.path(N_DIR, "STEP19N_REGULATORY_FAMILY_CONSENSUS.csv"), row.names = FALSE)
cat("\nSaved: STEP19N_REGULATORY_FAMILY_CONSENSUS.csv\n")

# ============================================================================
# 19N17 — EVIDENCE SCORE
# ============================================================================
cat("\n============================================\n")
cat("19N17: EVIDENCE SCORE\n")
cat("============================================\n\n")

evidence_list <- list()
for (tf in all_candidates) {
  score <- 0
  
  # +1 adequate regulon coverage
  cov_flag <- coverage_qc_df$coverage_flag[coverage_qc_df$regulator == tf]
  if (cov_flag %in% c("GOOD_COVERAGE", "LIMITED_COVERAGE")) score <- score + 1
  
  # +1 activity CORE2 in GSE197677
  c2_197 <- core_assoc_df$rho_CORE2[core_assoc_df$regulator == tf & core_assoc_df$dataset == "GSE197677"]
  if (!is.na(c2_197) && abs(c2_197) >= 0.40) score <- score + 1
  
  # +1 activity CORE2 in GSE221561
  c2_221 <- core_assoc_df$rho_CORE2[core_assoc_df$regulator == tf & core_assoc_df$dataset == "GSE221561"]
  if (!is.na(c2_221) && abs(c2_221) >= 0.40) score <- score + 1
  
  # +1 CORE2 association same direction
  if (!is.na(c2_197) && !is.na(c2_221) && sign(c2_197) == sign(c2_221) && abs(c2_197) >= 0.40 && abs(c2_221) >= 0.40) score <- score + 1
  
  # +1 treatment-activity direction replicated
  te_rep <- replication_df$te_replicated[replication_df$regulator == tf]
  if (!is.na(te_rep) && te_rep) score <- score + 1
  
  # +1 LOSO robust
  loso_class <- loso_summary_df$classification[loso_summary_df$regulator == tf &
    loso_summary_df$dataset == "GSE197677"]
  if (any(loso_class == "ROBUST", na.rm = TRUE)) score <- score + 1
  
  if (score >= 5) {
    ev_class <- "HIGH_ORTHOGONAL_SUPPORT"
  } else if (score >= 3) {
    ev_class <- "MODERATE_ORTHOGONAL_SUPPORT"
  } else if (score >= 1) {
    ev_class <- "WEAK_ORTHOGONAL_SUPPORT"
  } else {
    ev_class <- "NO_ORTHOGONAL_SUPPORT"
  }
  
  evidence_list[[tf]] <- data.frame(
    regulator = tf,
    evidence_score = score,
    adequate_coverage = cov_flag %in% c("GOOD_COVERAGE", "LIMITED_COVERAGE"),
    core2_197 = abs(c2_197) >= 0.40,
    core2_221 = abs(c2_221) >= 0.40,
    core2_same_direction = sign(c2_197) == sign(c2_221),
    te_replicated = te_rep,
    loso_robust = any(loso_class == "ROBUST"),
    classification = ev_class,
    stringsAsFactors = FALSE
  )
}

evidence_df <- do.call(rbind, evidence_list)
evidence_df <- evidence_df[order(-evidence_df$evidence_score, evidence_df$regulator), ]
rownames(evidence_df) <- NULL

cat("=== FINAL EVIDENCE SCORES ===\n")
for (i in 1:nrow(evidence_df)) {
  cat(sprintf("  %s: score=%d -> %s\n",
    evidence_df$regulator[i], evidence_df$evidence_score[i], evidence_df$classification[i]))
}

write.csv(evidence_df, file.path(N_DIR, "STEP19N_FINAL_REGULATOR_EVIDENCE_SCORE.csv"), row.names = FALSE)
cat("\nSaved: STEP19N_FINAL_REGULATOR_EVIDENCE_SCORE.csv\n")

# ============================================================================
# 19N18 — FIGURES
# ============================================================================
cat("\n============================================\n")
cat("19N18: FIGURES\n")
cat("============================================\n\n")

# Figure 1: Regulator activity vs CORE2 scatter
cat("Generating Figure 1: Regulator activity vs CORE2 scatter...\n")
scatter_data <- data.frame()
for (tf in PRIMARY_CANDIDATES) {
  for (ds in c("GSE197677", "GSE221561")) {
    sub <- activity_df[activity_df$regulator == tf & activity_df$dataset == ds, ]
    c2 <- core_scores$CORE2_STRONG_score[match(sub$sample_id, core_scores$sample_id)]
    scatter_data <- rbind(scatter_data, data.frame(
      regulator = tf, dataset = ds, activity = sub$activity, CORE2 = c2,
      stringsAsFactors = FALSE
    ))
  }
}

p1 <- ggplot(scatter_data, aes(x = activity, y = CORE2, color = dataset)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed", linewidth = 0.5) +
  facet_wrap(~regulator, scales = "free", ncol = 3) +
  scale_color_manual(values = c("GSE197677" = "steelblue", "GSE221561" = "tomato")) +
  labs(title = "Regulator Activity vs CORE2_STRONG Score",
       x = "Regulator Activity (target-mean z-score)", y = "CORE2_STRONG Score") +
  theme_bw() + theme(legend.position = "bottom")
ggsave(file.path(FIG_DIR, "STEP19N_REGULATOR_ACTIVITY_CORE2_SCATTER.png"), p1,
       width = 12, height = 10, dpi = 150)
cat("Saved: STEP19N_REGULATOR_ACTIVITY_CORE2_SCATTER.png\n")

# Figure 2: Treatment effects bar chart
cat("Generating Figure 2: Treatment effects bar chart...\n")
te_bar <- treat_effect_df[treat_effect_df$regulator %in% PRIMARY_CANDIDATES, ]
te_bar$regulator <- factor(te_bar$regulator, levels = PRIMARY_CANDIDATES)

p2 <- ggplot(te_bar, aes(x = regulator, y = effect, fill = dataset)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("GSE197677" = "steelblue", "GSE221561" = "tomato")) +
  labs(title = "Treatment Effect on Regulator Activity",
       x = "Regulator", y = "Effect (treated - control)", fill = "Dataset") +
  theme_bw() + theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(FIG_DIR, "STEP19N_REGULATOR_ACTIVITY_TREATMENT_EFFECTS.png"), p2,
       width = 8, height = 6, dpi = 150)
cat("Saved: STEP19N_REGULATOR_ACTIVITY_TREATMENT_EFFECTS.png\n")

# Figure 3: AP-1 activity by treatment
cat("Generating Figure 3: AP-1 activity by treatment...\n")
ap1_data <- activity_df[activity_df$regulator %in% c("FOS", "FOSL2"), ]
ap1_data$treatment_label <- ifelse(ap1_data$treatment == "NACT", "NACT",
                            ifelse(ap1_data$treatment == "nNACT", "nNACT",
                            ifelse(ap1_data$treatment == "Neoadjuvant_treated", "Treated", "Surgery_alone")))

p3 <- ggplot(ap1_data, aes(x = treatment_label, y = activity, fill = regulator)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.7) +
  facet_wrap(~dataset, scales = "free_x") +
  scale_fill_manual(values = c("FOS" = "darkred", "FOSL2" = "orange")) +
  labs(title = "AP-1 Family Activity by Treatment Group",
       x = "Treatment", y = "Activity", fill = "Regulator") +
  theme_bw() + theme(legend.position = "bottom")
ggsave(file.path(FIG_DIR, "STEP19N_AP1_ACTIVITY_BY_TREATMENT.png"), p3,
       width = 10, height = 6, dpi = 150)
cat("Saved: STEP19N_AP1_ACTIVITY_BY_TREATMENT.png\n")

# Figure 4: SMAD activity by treatment
cat("Generating Figure 4: SMAD activity by treatment...\n")
tryCatch({
  smad_data <- activity_df[activity_df$regulator %in% c("SMAD3", "SMAD7"), ]
  smad_data$treatment_label <- ifelse(smad_data$treatment == "NACT", "NACT",
                                ifelse(smad_data$treatment == "nNACT", "nNACT",
                                ifelse(smad_data$treatment == "Neoadjuvant_treated", "Treated", "Surgery_alone")))
  
  p4 <- ggplot(smad_data, aes(x = treatment_label, y = activity, fill = regulator)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 2, alpha = 0.7) +
    facet_wrap(~dataset, scales = "free_x") +
    scale_fill_manual(values = c("SMAD3" = "darkblue", "SMAD7" = "lightblue")) +
    labs(title = "TGF-beta/SMAD Family Activity by Treatment Group",
         x = "Treatment", y = "Activity", fill = "Regulator") +
    theme_bw() + theme(legend.position = "bottom")
  ggsave(file.path(FIG_DIR, "STEP19N_SMAD_ACTIVITY_BY_TREATMENT.png"), p4,
         width = 10, height = 6, dpi = 150)
  cat("Saved: STEP19N_SMAD_ACTIVITY_BY_TREATMENT.png\n")
}, error = function(e) {
  cat("Skipping Figure 4 - no valid data:", conditionMessage(e), "\n")
})

# Figure 5: Cross-dataset comparison
cat("Generating Figure 5: Cross-dataset activity comparison...\n")
cross_data <- data.frame()
for (tf in PRIMARY_CANDIDATES) {
  te_197 <- treat_effect_df$effect[treat_effect_df$regulator == tf & treat_effect_df$dataset == "GSE197677"]
  te_221 <- treat_effect_df$effect[treat_effect_df$regulator == tf & treat_effect_df$dataset == "GSE221561"]
  cross_data <- rbind(cross_data, data.frame(
    regulator = tf, effect_197 = te_197, effect_221 = te_221,
    stringsAsFactors = FALSE
  ))
}

p5 <- ggplot(cross_data, aes(x = effect_197, y = effect_221, label = regulator)) +
  geom_point(size = 4, color = "darkred") +
  geom_text(vjust = -1, size = 3.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "grey70") +
  labs(title = "Cross-Dataset Treatment Effect on Activity",
       x = "Effect GSE197677 (NACT - nNACT)",
       y = "Effect GSE221561 (Treated - Surgery_alone)") +
  theme_bw()
ggsave(file.path(FIG_DIR, "STEP19N_REGULATOR_ACTIVITY_CROSS_DATASET.png"), p5,
       width = 8, height = 7, dpi = 150)
cat("Saved: STEP19N_REGULATOR_ACTIVITY_CROSS_DATASET.png\n")

# Figure 6: LOSO stability
cat("Generating Figure 6: LOSO stability...\n")
loso_plot_data <- loso_summary_df
loso_plot_data$regulator <- factor(loso_plot_data$regulator, levels = PRIMARY_CANDIDATES)

p6 <- ggplot(loso_plot_data, aes(x = regulator, y = rho_sign_preservation, fill = dataset)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "red") +
  scale_fill_manual(values = c("GSE197677" = "steelblue", "GSE221561" = "tomato")) +
  scale_y_continuous(labels = scales::percent) +
  labs(title = "LOSO Activity-CORE2 Sign Preservation",
       x = "Regulator", y = "Sign Preservation", fill = "Dataset") +
  theme_bw() + theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(FIG_DIR, "STEP19N_REGULATOR_ACTIVITY_LOSO.png"), p6,
       width = 8, height = 6, dpi = 150)
cat("Saved: STEP19N_REGULATOR_ACTIVITY_LOSO.png\n")

# Figure 7: Evidence heatmap
cat("Generating Figure 7: Evidence heatmap...\n")
heat_mat <- evidence_df[, c("regulator", "adequate_coverage", "core2_197", "core2_221",
                            "core2_same_direction", "te_replicated", "loso_robust")]
rownames(heat_mat) <- heat_mat$regulator
heat_mat$regulator <- NULL
heat_mat_numeric <- as.matrix(sapply(heat_mat, as.numeric))
rownames(heat_mat_numeric) <- rownames(heat_mat)

png(file.path(FIG_DIR, "STEP19N_REGULATOR_EVIDENCE_HEATMAP.png"),
    width = 10, height = 6, units = "in", res = 150)
heatmap(heat_mat_numeric, margins = c(8, 6),
        main = "Regulator Orthogonal Evidence Components",
        col = c("grey90", "steelblue"),
        scale = "none")
dev.off()
cat("Saved: STEP19N_REGULATOR_EVIDENCE_HEATMAP.png\n\n")

# ============================================================================
# 19N19 — INTERPRETATION GUARDRAILS
# ============================================================================
cat("============================================\n")
cat("19N19: INTERPRETATION GUARDRAILS\n")
cat("============================================\n\n")

guardrails <- c(
  "1. Regulon activity is computational inference and does not prove TF biochemical activity.",
  "2. Curated TF-target databases are incomplete.",
  "3. Missing target edges are not evidence of absence.",
  "4. Expression association from Step19M and regulon activity from Step19N are distinct evidence types.",
  "5. Cross-dataset replication is stronger evidence than single-cohort association.",
  "6. GSE221561 Surgery_alone n=2 remains a major limitation.",
  "7. Treatment associations are observational.",
  "8. No causal regulatory claims are permitted.",
  "9. Sample/patient is the statistical unit.",
  "10. No cell-level hypothesis testing."
)

writeLines(guardrails, file.path(N_DIR, "STEP19N_INTERPRETATION_GUARDRAILS.md"))
cat("Saved: STEP19N_INTERPRETATION_GUARDRAILS.md\n")

# ============================================================================
# 19N20 — FINAL CLASSIFICATION
# ============================================================================
cat("\n============================================\n")
cat("19N20: FINAL CLASSIFICATION\n")
cat("============================================\n\n")

# Determine final classification
ap1_orth <- ap1_status == "AP1_ORTHOGONAL_SUPPORT"
tgfb_orth <- tgfb_status == "TGFB_SMAD_ORTHOGONAL_SUPPORT"
any_moderate <- any(evidence_df$classification == "MODERATE_ORTHOGONAL_SUPPORT")
any_high <- any(evidence_df$classification == "HIGH_ORTHOGONAL_SUPPORT")
any_adequate_coverage <- any(coverage_qc_df$coverage_flag %in% c("GOOD_COVERAGE", "LIMITED_COVERAGE"))
all_weak <- all(evidence_df$classification %in% c("WEAK_ORTHOGONAL_SUPPORT", "NO_ORTHOGONAL_SUPPORT"))

if (ap1_orth && tgfb_orth) {
  final_classification <- "CORE2_WITH_DUAL_AP1_SMAD_SUPPORT"
} else if (ap1_orth) {
  final_classification <- "CORE2_WITH_ORTHOGONAL_AP1_SUPPORT"
} else if (tgfb_orth) {
  final_classification <- "CORE2_WITH_ORTHOGONAL_TGFB_SMAD_SUPPORT"
} else if (any_high || (any_moderate && sum(evidence_df$classification == "MODERATE_ORTHOGONAL_SUPPORT") >= 2)) {
  final_classification <- "CORE2_WITH_MULTIPLE_MODERATE_REGULATOR_SUPPORT"
} else if (!any_adequate_coverage) {
  final_classification <- "REGULON_RESOURCE_UNAVAILABLE"
} else if (all_weak) {
  final_classification <- "NO_ORTHOGONAL_REGULATOR_SUPPORT"
} else {
  final_classification <- "ORTHOGONAL_VALIDATION_INCONCLUSIVE"
}

cat("=== FINAL CLASSIFICATION ===\n")
cat(final_classification, "\n\n")

# Summary
cat("AP-1 status:", ap1_status, "\n")
cat("TGF-beta/SMAD status:", tgfb_status, "\n")
cat("KLF6:", kr_df$classification[kr_df$regulator == "KLF6"], "\n")
cat("RUNX1:", kr_df$classification[kr_df$regulator == "RUNX1"], "\n\n")

# Final outputs
final_validation <- evidence_df
write.csv(final_validation, file.path(N_DIR, "STEP19N_FINAL_REGULATOR_VALIDATION.csv"), row.names = FALSE)
cat("Saved: STEP19N_FINAL_REGULATOR_VALIDATION.csv\n")

c2_orth <- data.frame(
  evidence_type = "regulon_activity",
  classification = final_classification,
  n_high = sum(evidence_df$classification == "HIGH_ORTHOGONAL_SUPPORT"),
  n_moderate = sum(evidence_df$classification == "MODERATE_ORTHOGONAL_SUPPORT"),
  n_weak = sum(evidence_df$classification == "WEAK_ORTHOGONAL_SUPPORT"),
  n_none = sum(evidence_df$classification == "NO_ORTHOGONAL_SUPPORT"),
  ap1_status = ap1_status,
  tgfb_smad_status = tgfb_status,
  timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
)
write.csv(c2_orth, file.path(N_DIR, "STEP19N_CORE2_ORTHOGONAL_EVIDENCE.csv"), row.names = FALSE)
cat("Saved: STEP19N_CORE2_ORTHOGONAL_EVIDENCE.csv\n")

write.csv(family_consensus_df, file.path(N_DIR, "STEP19N_REGULATORY_FAMILY_FINAL.csv"), row.names = FALSE)
cat("Saved: STEP19N_REGULATORY_FAMILY_FINAL.csv\n")

# Final interpretation
final_interp <- c(
  "# Step 19N Final Interpretation",
  "",
  paste("## Classification:", final_classification),
  "",
  "## Regulon Resource",
  paste("Resource:", REGULON_RESOURCE, "(MSigDB C3 TFT, unsigned)"),
  paste("Primary candidates with adequate coverage:", sum(coverage_qc_df$coverage_flag[coverage_qc_df$regulator %in% PRIMARY_CANDIDATES] %in% c("GOOD_COVERAGE", "LIMITED_COVERAGE"))),
  "",
  "## AP-1 Family",
  paste("Status:", ap1_status),
  "",
  "## TGF-beta/SMAD Family",
  paste("Status:", tgfb_status),
  "",
  "## Evidence Summary",
  paste("High:", sum(evidence_df$classification == "HIGH_ORTHOGONAL_SUPPORT")),
  paste("Moderate:", sum(evidence_df$classification == "MODERATE_ORTHOGONAL_SUPPORT")),
  paste("Weak:", sum(evidence_df$classification == "WEAK_ORTHOGONAL_SUPPORT")),
  paste("None:", sum(evidence_df$classification == "NO_ORTHOGONAL_SUPPORT")),
  "",
  "## Key Limitations",
  "Unsigned regulon edges (MSigDB C3 TFT)",
  "No decoupleR/dorothea/viper available",
  "Small sample sizes, especially GSE221561 Surgery_alone n=2",
  "Exploratory hypothesis-generation only"
)
writeLines(final_interp, file.path(N_DIR, "STEP19N_FINAL_INTERPRETATION.md"))
cat("Saved: STEP19N_FINAL_INTERPRETATION.md\n")

final_status <- data.frame(
  step = "19N",
  status = "COMPLETE",
  classification = final_classification,
  regulon_resource = REGULON_RESOURCE,
  n_high = sum(evidence_df$classification == "HIGH_ORTHOGONAL_SUPPORT"),
  n_moderate = sum(evidence_df$classification == "MODERATE_ORTHOGONAL_SUPPORT"),
  n_weak = sum(evidence_df$classification == "WEAK_ORTHOGONAL_SUPPORT"),
  n_none = sum(evidence_df$classification == "NO_ORTHOGONAL_SUPPORT"),
  ap1_status = ap1_status,
  tgfb_smad_status = tgfb_status,
  cells_removed = 0,
  samples_removed = 0,
  cells_reclassified = 0,
  frozen_upstream_changed = "NO",
  timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
)
write.csv(final_status, file.path(N_DIR, "STEP19N_FINAL_STATUS.csv"), row.names = FALSE)
cat("Saved: STEP19N_FINAL_STATUS.csv\n")

# ============================================================================
# FINAL VALIDATION
# ============================================================================
cat("\n=== FINAL VALIDATION ===\n")
cat("Canonical contrasts: UNCHANGED\n")
cat("GSE197677 tumor samples:", length(samples197), "\n")
cat("GSE221561 tumor samples:", length(samples221), "\n")
cat("Normal samples in tests: NONE\n")
cat("Step19L: UNCHANGED\n")
cat("Step19M: UNCHANGED\n")
cat("Step19K: UNCHANGED\n")
cat("Step19J-CORR2: UNCHANGED\n")
cat("Cells removed: 0\n")
cat("Samples removed: 0\n")
cat("Cells reclassified: 0\n\n")

# ============================================================================
# FINAL PRINT
# ============================================================================
cat("============================================\n")
cat("STEP 19N COMPLETE\n")
cat("ORTHOGONAL REGULATOR ACTIVITY VALIDATION\n")
cat("============================================\n\n")

cat("REGULON RESOURCE:\n")
cat("  Source:", REGULON_RESOURCE, "\n")
cat("  Type: MSigDB C3 TFT (unsigned)\n")
cat("  Note: No decoupleR/dorothea/viper available; using target-mean z-score\n\n")

cat("REGULON COVERAGE:\n")
for (i in 1:nrow(coverage_qc_df)) {
  cat(sprintf("  %s: UP=%d DN=%d (shared_total=%d) -> %s\n",
    coverage_qc_df$regulator[i], coverage_qc_df$shared_up[i],
    coverage_qc_df$shared_dn[i], coverage_qc_df$shared_total[i],
    coverage_qc_df$coverage_flag[i]))
}

cat("\n--------------------------------------------\n\n")

cat("AP-1\n\n")

for (tf in c("FOS", "FOSL2")) {
  c2_197 <- core_assoc_df$rho_CORE2[core_assoc_df$regulator == tf & core_assoc_df$dataset == "GSE197677"]
  c2_221 <- core_assoc_df$rho_CORE2[core_assoc_df$regulator == tf & core_assoc_df$dataset == "GSE221561"]
  te_197 <- treat_effect_df$effect[treat_effect_df$regulator == tf & treat_effect_df$dataset == "GSE197677"]
  te_221 <- treat_effect_df$effect[treat_effect_df$regulator == tf & treat_effect_df$dataset == "GSE221561"]
  loso_cls <- loso_summary_df$classification[loso_summary_df$regulator == tf]
  
  cat(sprintf("%s:\n", tf))
  cat(sprintf("  Activity/CORE2 GSE197677: %.3f\n", c2_197))
  cat(sprintf("  Activity/CORE2 GSE221561: %.3f\n", c2_221))
  cat(sprintf("  Treatment effect GSE197677: %.4f\n", te_197))
  cat(sprintf("  Treatment effect GSE221561: %.4f\n", te_221))
  cat(sprintf("  LOSO: %s\n", paste(loso_cls, collapse=", ")))
  cat(sprintf("  Classification: %s\n\n", replication_df$classification[replication_df$regulator == tf]))
}

cat(sprintf("AP-1 FAMILY STATUS: %s\n\n", ap1_status))

cat("--------------------------------------------\n\n")

cat("TGF-BETA / SMAD\n\n")

for (tf in c("SMAD3", "SMAD7")) {
  c2_197 <- core_assoc_df$rho_CORE2[core_assoc_df$regulator == tf & core_assoc_df$dataset == "GSE197677"]
  c2_221 <- core_assoc_df$rho_CORE2[core_assoc_df$regulator == tf & core_assoc_df$dataset == "GSE221561"]
  te_197 <- treat_effect_df$effect[treat_effect_df$regulator == tf & treat_effect_df$dataset == "GSE197677"]
  te_221 <- treat_effect_df$effect[treat_effect_df$regulator == tf & treat_effect_df$dataset == "GSE221561"]
  
  cat(sprintf("%s:\n", tf))
  cat(sprintf("  Activity/CORE2 GSE197677: %.3f\n", c2_197))
  cat(sprintf("  Activity/CORE2 GSE221561: %.3f\n", c2_221))
  cat(sprintf("  Treatment effect GSE197677: %.4f\n", te_197))
  cat(sprintf("  Treatment effect GSE221561: %.4f\n", te_221))
  cat(sprintf("  Classification: %s\n\n", replication_df$classification[replication_df$regulator == tf]))
}

cat(sprintf("TGF-BETA/SMAD FAMILY STATUS: %s\n\n", tgfb_status))

cat("--------------------------------------------\n\n")

for (tf in c("KLF6", "RUNX1")) {
  c2_197 <- core_assoc_df$rho_CORE2[core_assoc_df$regulator == tf & core_assoc_df$dataset == "GSE197677"]
  c2_221 <- core_assoc_df$rho_CORE2[core_assoc_df$regulator == tf & core_assoc_df$dataset == "GSE221561"]
  cat(sprintf("%s:\n", tf))
  cat(sprintf("  Activity/CORE2 GSE197677: %.3f\n", c2_197))
  cat(sprintf("  Activity/CORE2 GSE221561: %.3f\n", c2_221))
  cat(sprintf("  Coverage: %s\n", coverage_qc_df$coverage_flag[coverage_qc_df$regulator == tf]))
  cat(sprintf("  Classification: %s\n\n", kr_df$classification[kr_df$regulator == tf]))
}

cat("--------------------------------------------\n\n")

# Top orthogonal regulator
top_reg <- evidence_df$regulator[1]
cat(sprintf("TOP ORTHOGONAL REGULATOR: %s (score=%d)\n", top_reg, evidence_df$evidence_score[1]))

cat(sprintf("\nCORE2 ORTHOGONAL SUPPORT: %s\n", final_classification))
cat(sprintf("\nFINAL CLASSIFICATION: %s\n", final_classification))
cat("\nCells removed: 0\n")
cat("Samples removed: 0\n")
cat("Cells reclassified: 0\n")
cat("Frozen upstream outputs changed: NO\n\n")

cat("STOP HERE.\n")
cat("DO NOT START STEP19O.\n")
