#!/usr/bin/env Rscript
# ============================================================================
# STEP 19L: CONSERVED FIBROBLAST CORE MODULE VALIDATION
# ============================================================================
# Validates BGN(+), CDKN1B(-), GSN(-), TIMP1(+) as a conserved fibroblast
# core module across GSE197677 and GSE221561 canonical contrasts.
# ============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
})

BASE_DIR <- Sys.getenv("ESCC_CAF_PROJECT_ROOT", unset = ".")
L_DIR    <- file.path(BASE_DIR, "06_tables/cross_dataset/pathway_analysis/Step19L")
FIG_DIR  <- file.path(BASE_DIR, "05_figures/cross_dataset/pathway_analysis/Step19L")
LOG_DIR  <- file.path(BASE_DIR, "08_logs")
dir.create(L_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(LOG_DIR, showWarnings = FALSE, recursive = TRUE)

CORE_GENES <- c("BGN", "CDKN1B", "GSN", "TIMP1")

TARGET_PATHWAYS <- c("HALLMARK_HYPOXIA", "HALLMARK_COAGULATION",
                     "HALLMARK_KRAS_SIGNALING_UP", "HALLMARK_APOPTOSIS")
PATHWAY_SHORT <- gsub("HALLMARK_", "", TARGET_PATHWAYS)

FROZEN_LOGFC <- list(
  GSE197677  = c(BGN = 1.3665, CDKN1B = -0.8393, GSN = -0.7420, TIMP1 = 0.6126),
  GSE221561  = c(BGN = 0.0108, CDKN1B = -0.6750, GSN = -0.5219, TIMP1 = 0.0255)
)

FROZEN_SIGN <- c(BGN = +1, CDKN1B = -1, GSN = -1, TIMP1 = +1)

TREATMENT_MAP <- list(
  GSE197677 = c(
    ESC06  = "nNACT", ESC07  = "nNACT", ESC09  = "NACT",
    ESC11  = "NACT",  ESC12  = "NACT",  ESC15  = "NACT",
    ESC16  = "NACT",  ESC17  = "nNACT", ESC18  = "nNACT",
    ESC22  = "NACT"
  ),
  GSE221561 = c(
    S0520T  = "Neoadjuvant_treated",  S1265T   = "Neoadjuvant_treated",
    S1315TS = "Neoadjuvant_treated",  S1535T   = "Neoadjuvant_treated",
    S2423T  = "Neoadjuvant_treated",  S2487T   = "Neoadjuvant_treated",
    S6829T  = "Neoadjuvant_treated",  S6417T   = "Surgery_alone",
    S9478T  = "Surgery_alone"
  )
)

TREATMENT_CONTRAST <- list(
  GSE197677  = c(treated = "NACT",      control = "nNACT"),
  GSE221561  = c(treated = "Neoadjuvant_treated", control = "Surgery_alone")
)

compute_log2_cpm <- function(counts_mat) {
  lib_size <- colSums(counts_mat)
  lib_size[lib_size == 0] <- 1
  cpm_mat <- sweep(counts_mat, 2, lib_size, "/") * 1e6
  log2(cpm_mat + 1)
}

# ============================================================================
# 19L1: INPUT PROVENANCE
# ============================================================================
cat("============================================\n")
cat("19L1: INPUT PROVENANCE\n")
cat("============================================\n\n")

fib197_file <- file.path(BASE_DIR,
  "06_tables/cross_dataset/pathway_analysis/Step19I/counts/GSE197677_FIBROBLAST_CAF_raw_counts.rds")
fib221_file <- file.path(BASE_DIR,
  "06_tables/cross_dataset/pathway_analysis/Step19I/counts/GSE221561_FIBROBLAST_CAF_raw_counts.rds")
corr197_file <- file.path(BASE_DIR,
  "06_tables/cross_dataset/pathway_analysis/Step19J-CORR2/STEP19J_CORR2_GSE197677_FIBROBLAST_TARGET_GENE_EFFECTS_CANONICAL.csv")
corr221_file <- file.path(BASE_DIR,
  "06_tables/cross_dataset/pathway_analysis/Step19J-CORR2/STEP19J_CORR2_GSE221561_FIBROBLAST_TARGET_GENE_EFFECTS_CANONICAL.csv")
inventory_file <- file.path(BASE_DIR,
  "06_tables/cross_dataset/pathway_analysis/Step19I/STEP19I_COMPARTMENT_SAMPLE_INVENTORY.csv")
pathway_scores_file <- file.path(BASE_DIR,
  "06_tables/cross_dataset/pathway_analysis/Step19I/STEP19I_SAMPLE_COMPARTMENT_PATHWAY_SCORES.csv")
epi197_file <- file.path(BASE_DIR,
  "06_tables/cross_dataset/pathway_analysis/Step19I/counts/GSE197677_EPITHELIAL_raw_counts.rds")
epi221_file <- file.path(BASE_DIR,
  "06_tables/cross_dataset/pathway_analysis/Step19I/counts/GSE221561_EPITHELIAL_raw_counts.rds")
gse160269_file <- file.path(BASE_DIR,
  "03_objects/GSE160269/GSE160269_annotated_merged_raw.rds")

stopifnot(file.exists(fib197_file), file.exists(fib221_file))
stopifnot(file.exists(corr197_file), file.exists(corr221_file))
stopifnot(file.exists(inventory_file), file.exists(pathway_scores_file))

fib197 <- readRDS(fib197_file)
fib221 <- readRDS(fib221_file)
corr197 <- read.csv(corr197_file)
corr221 <- read.csv(corr221_file)
inventory <- read.csv(inventory_file)
pathway_scores <- read.csv(pathway_scores_file)

cat("GSE197677 fibroblast:", nrow(fib197), "genes x", ncol(fib197), "samples\n")
cat("GSE221561 fibroblast:", nrow(fib221), "genes x", ncol(fib221), "samples\n")

missing_197 <- setdiff(CORE_GENES, rownames(fib197))
missing_221 <- setdiff(CORE_GENES, rownames(fib221))
if (length(missing_197) > 0) stop("Missing core genes in GSE197677: ", paste(missing_197, collapse=", "))
if (length(missing_221) > 0) stop("Missing core genes in GSE221561: ", paste(missing_221, collapse=", "))
cat("All 4 core genes present in both datasets.\n\n")

prov <- data.frame(
  file = c("fib197", "fib221", "corr197", "corr221", "inventory", "pathway_scores"),
  path = c(fib197_file, fib221_file, corr197_file, corr221_file, inventory_file, pathway_scores_file),
  exists = TRUE,
  rows = c(nrow(fib197), nrow(fib221), nrow(corr197), nrow(corr221), nrow(inventory), nrow(pathway_scores)),
  cols = c(ncol(fib197), ncol(fib221), ncol(corr197), ncol(corr221), ncol(inventory), ncol(pathway_scores))
)
write.csv(prov, file.path(L_DIR, "STEP19L_INPUT_PROVENANCE.csv"), row.names = FALSE)

# ============================================================================
# 19L2: VERIFY CANONICAL GENE EFFECTS
# ============================================================================
cat("============================================\n")
cat("19L2: VERIFY CANONICAL GENE EFFECTS\n")
cat("============================================\n\n")

verify_signs <- function(corr_df, dataset_label) {
  sub <- corr_df[corr_df$gene %in% CORE_GENES, ]
  sub <- sub[match(CORE_GENES, sub$gene), ]
  expected_sign <- FROZEN_SIGN[CORE_GENES]
  sign_match <- (sign(sub$logFC) == expected_sign)
  log2fc_match <- all(abs(sub$logFC - FROZEN_LOGFC[[dataset_label]][CORE_GENES]) < 0.01)
  result <- data.frame(
    gene = CORE_GENES,
    frozen_logFC = FROZEN_LOGFC[[dataset_label]][CORE_GENES],
    corr2_logFC = sub$logFC,
    PValue = sub$PValue,
    FDR = sub$FDR,
    canonical_direction = sub$canonical_direction,
    expected_sign = expected_sign,
    actual_logFC_sign = sign(sub$logFC),
    sign_match = sign_match,
    log2fc_frozen_match = log2fc_match,
    dataset = dataset_label,
    stringsAsFactors = FALSE
  )
  mismatches <- sum(!sign_match)
  if (mismatches > 0) {
    cat("STOP: Sign mismatch in", dataset_label, "for:", paste(sub$gene[!sign_match], collapse=", "), "\n")
    stop("Canonical sign mismatch in ", dataset_label)
  }
  cat(dataset_label, ": all", length(CORE_GENES), "core gene signs MATCH.\n")
  cat("  BGN =", round(sub$logFC[sub$gene == "BGN"], 4),
      " CDKN1B =", round(sub$logFC[sub$gene == "CDKN1B"], 4),
      " GSN =", round(sub$logFC[sub$gene == "GSN"], 4),
      " TIMP1 =", round(sub$logFC[sub$gene == "TIMP1"], 4), "\n")
  return(result)
}

v197 <- verify_signs(corr197, "GSE197677")
v221 <- verify_signs(corr221, "GSE221561")

cat("\nIndependent verification via log2 CPM from pseudobulk:\n")
for (ds_name in c("GSE197677", "GSE221561")) {
  counts <- if (ds_name == "GSE197677") fib197 else fib221
  tm <- TREATMENT_MAP[[ds_name]]
  tumor_ids <- names(tm)
  tumor_ids <- intersect(tumor_ids, colnames(counts))
  treated_ids <- names(tm[tm == TREATMENT_CONTRAST[[ds_name]]["treated"]])
  control_ids <- names(tm[tm == TREATMENT_CONTRAST[[ds_name]]["control"]])
  treated_ids <- intersect(treated_ids, tumor_ids)
  control_ids <- intersect(control_ids, tumor_ids)
  l2c <- compute_log2_cpm(counts)
  for (g in CORE_GENES) {
    t_mean <- mean(l2c[g, treated_ids])
    c_mean <- mean(l2c[g, control_ids])
    diff <- t_mean - c_mean
    cat("  ", ds_name, g, ": treated_mean =", round(t_mean, 3),
        " control_mean =", round(c_mean, 3),
        " diff =", round(diff, 3),
        " sign =", ifelse(diff > 0, "+", "-"), "\n")
  }
}

verify_combined <- rbind(v197, v221)
write.csv(verify_combined, file.path(L_DIR, "STEP19L_CORE_GENE_EFFECT_VERIFICATION.csv"), row.names = FALSE)
cat("\n")

# ============================================================================
# 19L3: BUILD SAMPLE-LEVEL CORE MODULE SCORE
# ============================================================================
cat("============================================\n")
cat("19L3: BUILD SAMPLE-LEVEL CORE MODULE SCORE\n")
cat("============================================\n\n")

build_core_scores <- function(counts_mat, dataset_label) {
  tm <- TREATMENT_MAP[[dataset_label]]
  tumor_ids <- names(tm)
  tumor_ids <- intersect(tumor_ids, colnames(counts_mat))
  l2c <- compute_log2_cpm(counts_mat)
  l2c_tumor <- l2c[, tumor_ids, drop = FALSE]
  gene_means <- rowMeans(l2c_tumor)
  gene_sds <- apply(l2c_tumor, 1, sd)
  gene_sds[gene_sds < 1e-10] <- 1e-10
  z_mat <- sweep(l2c_tumor, 1, gene_means) / gene_sds

  results <- data.frame(
    dataset = dataset_label,
    sample_id = tumor_ids,
    treatment = tm[tumor_ids],
    BGN_z = z_mat["BGN", ],
    CDKN1B_z = z_mat["CDKN1B", ],
    GSN_z = z_mat["GSN", ],
    TIMP1_z = z_mat["TIMP1", ],
    stringsAsFactors = FALSE
  )
  results$CORE4_score <- (results$BGN_z + (-1) * results$CDKN1B_z +
                          (-1) * results$GSN_z + results$TIMP1_z) / 4
  results$CORE2_STRONG_score <- ((-1) * results$CDKN1B_z + (-1) * results$GSN_z) / 2
  return(results)
}

scores_197 <- build_core_scores(fib197, "GSE197677")
scores_221 <- build_core_scores(fib221, "GSE221561")

all_scores <- rbind(scores_197, scores_221)
write.csv(all_scores, file.path(L_DIR, "STEP19L_SAMPLE_CORE_MODULE_SCORES.csv"), row.names = FALSE)

cat("GSE197677:", nrow(scores_197), "tumor samples scored\n")
cat("  Treated (NACT):", sum(scores_197$treatment == "NACT"), "\n")
cat("  Control (nNACT):", sum(scores_197$treatment == "nNACT"), "\n")
cat("  CORE4 mean:", round(mean(scores_197$CORE4_score), 4), "\n")
cat("  CORE2_STRONG mean:", round(mean(scores_197$CORE2_STRONG_score), 4), "\n\n")

cat("GSE221561:", nrow(scores_221), "tumor samples scored\n")
cat("  Treated (Neoadjuvant_treated):", sum(scores_221$treatment == "Neoadjuvant_treated"), "\n")
cat("  Control (Surgery_alone):", sum(scores_221$treatment == "Surgery_alone"), "\n")
cat("  CORE4 mean:", round(mean(scores_221$CORE4_score), 4), "\n")
cat("  CORE2_STRONG mean:", round(mean(scores_221$CORE2_STRONG_score), 4), "\n\n")

# ============================================================================
# 19L4: CANONICAL MODULE TREATMENT EFFECT
# ============================================================================
cat("============================================\n")
cat("19L4: CANONICAL MODULE TREATMENT EFFECT\n")
cat("============================================\n\n")

compute_effect <- function(scores_df, dataset_label) {
  tc <- TREATMENT_CONTRAST[[dataset_label]]
  t_scores <- scores_df$CORE4_score[scores_df$treatment == tc["treated"]]
  c_scores <- scores_df$CORE4_score[scores_df$treatment == tc["control"]]
  t2_scores <- scores_df$CORE2_STRONG_score[scores_df$treatment == tc["treated"]]
  c2_scores <- scores_df$CORE2_STRONG_score[scores_df$treatment == tc["control"]]
  data.frame(
    dataset = dataset_label,
    module = c("CORE4", "CORE2_STRONG"),
    treated_mean = c(mean(t_scores), mean(t2_scores)),
    control_mean = c(mean(c_scores), mean(c2_scores)),
    treated_n = c(length(t_scores), length(t_scores)),
    control_n = c(length(c_scores), length(c2_scores)),
    effect = c(mean(t_scores) - mean(c_scores), mean(t2_scores) - mean(c2_scores)),
    treated_sd = c(sd(t_scores), sd(t2_scores)),
    control_sd = c(sd(c_scores), sd(c2_scores)),
    stringsAsFactors = FALSE
  )
}

eff197 <- compute_effect(scores_197, "GSE197677")
eff221 <- compute_effect(scores_221, "GSE221561")
effects <- rbind(eff197, eff221)

effects$core4_direction <- ifelse(effects$module == "CORE4",
  ifelse(effects$effect > 0, "TREATED_HIGH", "CONTROL_HIGH"), NA)
effects$core4_conserved <- ifelse(effects$module == "CORE4",
  effects$effect > 0, NA)
effects$conserved_classification <- NA
core4_pos <- effects[effects$module == "CORE4", ]
if (all(core4_pos$effect > 0)) {
  effects$conserved_classification[effects$module == "CORE4"] <- "CONSERVED_POSITIVE"
} else {
  effects$conserved_classification[effects$module == "CORE4"] <- "NOT_CONSERVED"
}

cat("CORE4 Effects:\n")
for (i in 1:nrow(effects)) {
  if (effects$module[i] == "CORE4") {
    cat("  ", effects$dataset[i], ": effect =", round(effects$effect[i], 4),
        " (treated", round(effects$treated_mean[i], 4),
        " vs control", round(effects$control_mean[i], 4), ")\n")
  }
}
cat("Conserved classification:", effects$conserved_classification[effects$module == "CORE4"][1], "\n\n")

write.csv(effects, file.path(L_DIR, "STEP19L_CORE_MODULE_CANONICAL_EFFECTS.csv"), row.names = FALSE)

# ============================================================================
# 19L5: EXACT WITHIN-DATASET PERMUTATION TESTS
# ============================================================================
cat("============================================\n")
cat("19L5: EXACT WITHIN-DATASET PERMUTATION TESTS\n")
cat("============================================\n\n")

exact_perm_test <- function(scores_df, dataset_label) {
  tc <- TREATMENT_CONTRAST[[dataset_label]]
  treated_ids <- scores_df$sample_id[scores_df$treatment == tc["treated"]]
  control_ids <- scores_df$sample_id[scores_df$treatment == tc["control"]]
  n_treated <- length(treated_ids)
  n_control <- length(control_ids)
  n_total <- n_treated + n_control

  all_perms <- combn(n_total, n_treated)

  results <- data.frame()
  for (mod_name in c("CORE4_score", "CORE2_STRONG_score")) {
    observed_treated_mean <- mean(scores_df[[mod_name]][scores_df$treatment == tc["treated"]])
    observed_control_mean <- mean(scores_df[[mod_name]][scores_df$treatment == tc["control"]])
    observed_diff <- observed_treated_mean - observed_control_mean

    perm_diffs <- numeric(ncol(all_perms))
    all_values <- scores_df[[mod_name]]
    for (p in 1:ncol(all_perms)) {
      perm_treated <- all_values[all_perms[, p]]
      perm_control <- all_values[setdiff(1:n_total, all_perms[, p])]
      perm_diffs[p] <- mean(perm_treated) - mean(perm_control)
    }

    p_right <- mean(perm_diffs >= observed_diff)
    p_left <- mean(perm_diffs <= observed_diff)
    p_two <- 2 * min(p_right, p_left)

    module_label <- gsub("_score", "", mod_name)
    results <- rbind(results, data.frame(
      dataset = dataset_label,
      module = module_label,
      n_treated = n_treated,
      n_control = n_control,
      total_perms = ncol(all_perms),
      observed_treated_mean = round(observed_treated_mean, 4),
      observed_control_mean = round(observed_control_mean, 4),
      observed_diff = round(observed_diff, 4),
      p_right = round(p_right, 4),
      p_left = round(p_left, 4),
      p_two_sided = round(p_two, 4),
      perm_mean = round(mean(perm_diffs), 4),
      perm_sd = round(sd(perm_diffs), 4),
      stringsAsFactors = FALSE
    ))
  }
  return(results)
}

perm_197 <- exact_perm_test(scores_197, "GSE197677")
perm_221 <- exact_perm_test(scores_221, "GSE221561")
perm_results <- rbind(perm_197, perm_221)
write.csv(perm_results, file.path(L_DIR, "STEP19L_CORE_MODULE_EXACT_TESTS.csv"), row.names = FALSE)

cat("GSE197677 exact permutation (C(10,6) =", choose(10, 6), "perms):\n")
cat("  CORE4: observed_diff =", perm_197$observed_diff[perm_197$module == "CORE4"],
    " P(two) =", perm_197$p_two_sided[perm_197$module == "CORE4"], "\n")
cat("  CORE2_STRONG: observed_diff =", perm_197$observed_diff[perm_197$module == "CORE2_STRONG"],
    " P(two) =", perm_197$p_two_sided[perm_197$module == "CORE2_STRONG"], "\n\n")

cat("GSE221561 exact permutation (C(9,7) =", choose(9, 7), "perms):\n")
cat("  CORE4: observed_diff =", perm_221$observed_diff[perm_221$module == "CORE4"],
    " P(two) =", perm_221$p_two_sided[perm_221$module == "CORE4"], "\n")
cat("  CORE2_STRONG: observed_diff =", perm_221$observed_diff[perm_221$module == "CORE2_STRONG"],
    " P(two) =", perm_221$p_two_sided[perm_221$module == "CORE2_STRONG"], "\n\n")

# ============================================================================
# 19L6: LEAVE-ONE-SAMPLE-OUT ROBUSTNESS
# ============================================================================
cat("============================================\n")
cat("19L6: LEAVE-ONE-SAMPLE-OUT ROBUSTNESS\n")
cat("============================================\n\n")

run_loso <- function(scores_df, dataset_label) {
  tc <- TREATMENT_CONTRAST[[dataset_label]]
  loso_long <- data.frame()

  for (leave_out in scores_df$sample_id) {
    reduced <- scores_df[scores_df$sample_id != leave_out, ]
    treated_remaining <- reduced$treatment == tc["treated"]
    control_remaining <- reduced$treatment == tc["control"]

    if (sum(treated_remaining) == 0 || sum(control_remaining) == 0) next

    for (mod in c("CORE4_score", "CORE2_STRONG_score")) {
      t_mean <- mean(reduced[[mod]][treated_remaining])
      c_mean <- mean(reduced[[mod]][control_remaining])
      diff <- t_mean - c_mean
      direction <- ifelse(diff > 0, "+", "-")
      module_label <- gsub("_score", "", mod)
      loso_long <- rbind(loso_long, data.frame(
        dataset = dataset_label,
        module = module_label,
        left_out_sample = leave_out,
        left_out_treatment = scores_df$treatment[scores_df$sample_id == leave_out],
        n_treated = sum(treated_remaining),
        n_control = sum(control_remaining),
        effect_diff = round(diff, 4),
        direction = direction,
        stringsAsFactors = FALSE
      ))
    }
  }

  full_diffs <- c()
  for (mod in c("CORE4_score", "CORE2_STRONG_score")) {
    t_mean <- mean(scores_df[[mod]][scores_df$treatment == tc["treated"]])
    c_mean <- mean(scores_df[[mod]][scores_df$treatment == tc["control"]])
    full_diffs <- c(full_diffs, t_mean - c_mean)
  }
  full_direction <- ifelse(full_diffs > 0, "+", "-")
  names(full_direction) <- c("CORE4", "CORE2_STRONG")

  return(list(long = loso_long, full_direction = full_direction))
}

loso_197 <- run_loso(scores_197, "GSE197677")
loso_221 <- run_loso(scores_221, "GSE221561")
loso_all <- rbind(loso_197$long, loso_221$long)

full_dir_map <- c(
  GSE197677_CORE4 = loso_197$full_direction["CORE4"],
  GSE197677_CORE2_STRONG = loso_197$full_direction["CORE2_STRONG"],
  GSE221561_CORE4 = loso_221$full_direction["CORE4"],
  GSE221561_CORE2_STRONG = loso_221$full_direction["CORE2_STRONG"]
)

loso_summary <- data.frame()
for (ds in c("GSE197677", "GSE221561")) {
  for (mod in c("CORE4", "CORE2_STRONG")) {
    sub <- loso_all[loso_all$dataset == ds & loso_all$module == mod, ]
    total <- nrow(sub)
    if (total == 0) {
      loso_summary <- rbind(loso_summary, data.frame(
        dataset = ds, module = mod, n_total = 0, n_stable = 0,
        frac_stable = NA, n_control_left_out = 0, n_control_stable = 0,
        control_all_stable = NA, classification = "MODULE_UNSTABLE",
        note = "no LOSO data", stringsAsFactors = FALSE
      ))
      next
    }
    ref_dir <- unname(as.character(full_dir_map[paste0(ds, "_", mod)]))
    stable <- sum(sub$direction == ref_dir, na.rm = TRUE)
    control_sub <- sub[sub$left_out_treatment == as.character(unname(TREATMENT_CONTRAST[[ds]]["control"])), ]
    control_stable <- sum(control_sub$direction == ref_dir, na.rm = TRUE)
    frac_stable <- stable / total
    control_all_stable <- isTRUE((nrow(control_sub) == 0) || (control_stable == nrow(control_sub)))

    if (ds == "GSE221561" && mod == "CORE2_STRONG") {
      classification <- "MODULE_CONTROL_SENSITIVE"
      note <- "GSE221561 control n=2 sensitivity only"
    } else if (!is.na(frac_stable) && frac_stable >= 0.8 && isTRUE(control_all_stable)) {
      classification <- "LOSO_ROBUST"
      note <- ""
    } else if (isTRUE(control_all_stable)) {
      classification <- "MODULE_OUTLIER_SENSITIVE"
      note <- ""
    } else {
      classification <- "MODULE_UNSTABLE"
      note <- ""
    }

    loso_summary <- rbind(loso_summary, data.frame(
      dataset = ds,
      module = mod,
      n_total = total,
      n_stable = stable,
      frac_stable = round(frac_stable, 3),
      n_control_left_out = nrow(control_sub),
      n_control_stable = control_stable,
      control_all_stable = control_all_stable,
      classification = classification,
      note = note,
      stringsAsFactors = FALSE
    ))
  }
}

write.csv(loso_all, file.path(L_DIR, "STEP19L_CORE_MODULE_LOSO_LONG.csv"), row.names = FALSE)
write.csv(loso_summary, file.path(L_DIR, "STEP19L_CORE_MODULE_LOSO_SUMMARY.csv"), row.names = FALSE)

cat("LOSO Results:\n")
for (i in 1:nrow(loso_summary)) {
  cat("  ", loso_summary$dataset[i], loso_summary$module[i], ":",
      loso_summary$n_stable[i], "/", loso_summary$n_total[i], "stable (",
      round(loso_summary$frac_stable[i] * 100, 1), "%) ->",
      loso_summary$classification[i], "\n")
}
cat("\n")

# ============================================================================
# 19L7: LEAVE-ONE-GENE-OUT
# ============================================================================
cat("============================================\n")
cat("19L7: LEAVE-ONE-GENE-OUT\n")
cat("============================================\n\n")

build_loGO_scores <- function(counts_mat, dataset_label) {
  tm <- TREATMENT_MAP[[dataset_label]]
  tumor_ids <- names(tm)
  tumor_ids <- intersect(tumor_ids, colnames(counts_mat))
  l2c <- compute_log2_cpm(counts_mat)
  l2c_tumor <- l2c[, tumor_ids, drop = FALSE]
  gene_means <- rowMeans(l2c_tumor)
  gene_sds <- apply(l2c_tumor, 1, sd)
  gene_sds[gene_sds < 1e-10] <- 1e-10
  z_mat <- sweep(l2c_tumor, 1, gene_means) / gene_sds

  orients <- c(BGN = +1, CDKN1B = -1, GSN = -1, TIMP1 = +1)

  results <- data.frame()
  for (leave_gene in c(CORE_GENES, "NONE")) {
    remaining <- setdiff(CORE_GENES, leave_gene)
    if (leave_gene == "NONE") remaining <- CORE_GENES

    for (i in seq_along(tumor_ids)) {
      sid <- tumor_ids[i]
      mod_score <- mean(orients[remaining] * z_mat[remaining, sid])
      results <- rbind(results, data.frame(
        dataset = dataset_label,
        sample_id = sid,
        treatment = tm[sid],
        left_out_gene = leave_gene,
        n_genes = length(remaining),
        module_score = round(mod_score, 4),
        stringsAsFactors = FALSE
      ))
    }
  }
  return(results)
}

loGO_197 <- build_loGO_scores(fib197, "GSE197677")
loGO_221 <- build_loGO_scores(fib221, "GSE221561")
loGO_all <- rbind(loGO_197, loGO_221)

loGO_effects <- data.frame()
for (ds in c("GSE197677", "GSE221561")) {
  tc <- TREATMENT_CONTRAST[[ds]]
  sub <- loGO_all[loGO_all$dataset == ds, ]
  for (lg in unique(sub$left_out_gene)) {
    lg_sub <- sub[sub$left_out_gene == lg, ]
    t_mean <- mean(lg_sub$module_score[lg_sub$treatment == tc["treated"]])
    c_mean <- mean(lg_sub$module_score[lg_sub$treatment == tc["control"]])
    diff <- t_mean - c_mean
    loGO_effects <- rbind(loGO_effects, data.frame(
      dataset = ds,
      left_out_gene = lg,
      n_genes = lg_sub$n_genes[1],
      treated_mean = round(t_mean, 4),
      control_mean = round(c_mean, 4),
      effect = round(diff, 4),
      direction = ifelse(diff > 0, "+", "-"),
      stringsAsFactors = FALSE
    ))
  }
}

full_directions <- loGO_effects$direction[loGO_effects$left_out_gene == "NONE"]
names(full_directions) <- c("GSE197677", "GSE221561")

loGO_effects$gene_importance <- NA
for (i in 1:nrow(loGO_effects)) {
  if (loGO_effects$left_out_gene[i] == "NONE") {
    loGO_effects$gene_importance[i] <- "FULL_MODULE"
  } else if (loGO_effects$direction[i] == full_directions[loGO_effects$dataset[i]]) {
    if (abs(loGO_effects$effect[i]) < abs(loGO_effects$effect[loGO_effects$left_out_gene == "NONE" & loGO_effects$dataset == loGO_effects$dataset[i]]) * 0.5) {
      loGO_effects$gene_importance[i] <- "IMPORTANT_FOR_MAGNITUDE"
    } else {
      loGO_effects$gene_importance[i] <- "NONESSENTIAL_FOR_DIRECTION"
    }
  } else {
    loGO_effects$gene_importance[i] <- "DIRECTION_DEPENDENT_ON_GENE"
  }
}

write.csv(loGO_all, file.path(L_DIR, "STEP19L_LEAVE_ONE_GENE_OUT.csv"), row.names = FALSE)

cat("Leave-one-gene-out effects:\n")
for (i in 1:nrow(loGO_effects)) {
  cat("  ", loGO_effects$dataset[i], "- left_out:", loGO_effects$left_out_gene[i],
      " effect =", loGO_effects$effect[i], " dir:", loGO_effects$direction[i],
      " ->", loGO_effects$gene_importance[i], "\n")
}
cat("\n")

# ============================================================================
# 19L8: STRONG VS WEAK CORE MODULE
# ============================================================================
cat("============================================\n")
cat("19L8: STRONG VS WEAK CORE MODULE\n")
cat("============================================\n\n")

all_scores$CORE2_WEAK_score <- (all_scores$BGN_z + all_scores$TIMP1_z) / 2

strong_weak <- data.frame()
for (ds in c("GSE197677", "GSE221561")) {
  tc <- TREATMENT_CONTRAST[[ds]]
  sub <- all_scores[all_scores$dataset == ds, ]
  for (mod in c("CORE2_STRONG_score", "CORE2_WEAK_score")) {
    t_mean <- mean(sub[[mod]][sub$treatment == tc["treated"]])
    c_mean <- mean(sub[[mod]][sub$treatment == tc["control"]])
    diff <- t_mean - c_mean
    mod_label <- gsub("_score", "", mod)
    strong_weak <- rbind(strong_weak, data.frame(
      dataset = ds,
      module = mod_label,
      treated_mean = round(t_mean, 4),
      control_mean = round(c_mean, 4),
      effect = round(diff, 4),
      direction = ifelse(diff > 0, "+", "-"),
      stringsAsFactors = FALSE
    ))
  }
}

cat("Strong (CDKN1B+GSN) vs Weak (BGN+TIMP1):\n")
for (i in 1:nrow(strong_weak)) {
  cat("  ", strong_weak$dataset[i], strong_weak$module[i], ":",
      "effect =", strong_weak$effect[i], "dir:", strong_weak$direction[i], "\n")
}
cat("\n")

write.csv(strong_weak, file.path(L_DIR, "STEP19L_STRONG_VS_WEAK_CORE_MODULE.csv"), row.names = FALSE)

# ============================================================================
# 19L9: ABUNDANCE ADJUSTMENT
# ============================================================================
cat("============================================\n")
cat("19L9: ABUNDANCE ADJUSTMENT\n")
cat("============================================\n\n")

abund_results <- data.frame()
for (ds in c("GSE197677", "GSE221561")) {
  tc <- TREATMENT_CONTRAST[[ds]]
  tm <- TREATMENT_MAP[[ds]]
  sub <- all_scores[all_scores$dataset == ds, ]

  inv_sub <- inventory[inventory$dataset == ds & inventory$sample %in% sub$sample_id, ]
  sub <- merge(sub, inv_sub[, c("sample", "fibroblast_cells", "epithelial_cells")], by.x = "sample_id", by.y = "sample", all.x = TRUE)
  sub$fib_fraction <- sub$fibroblast_cells / (sub$epithelial_cells + sub$fibroblast_cells)

  for (mod in c("CORE4_score", "CORE2_STRONG_score")) {
    mod_label <- gsub("_score", "", mod)

    model_res <- tryCatch({
      f1 <- reformulate("treatment", response = mod)
      m1 <- lm(f1, data = sub)
      coef_trt1 <- coef(m1)[2]
      p_trt1 <- summary(m1)$coefficients[2, 4]

      f2 <- reformulate(c("treatment", "scale(fib_fraction)"), response = mod)
      m2 <- lm(f2, data = sub)
      coef_trt2 <- coef(m2)[2]
      p_trt2 <- summary(m2)$coefficients[2, 4]

      list(unadjusted_effect = coef_trt1, unadjusted_p = p_trt1,
           adjusted_effect = coef_trt2, adjusted_p = p_trt2,
           fib_fraction_median = round(median(sub$fib_fraction, na.rm = TRUE), 4))
    }, error = function(e) {
      cat("  Warning: model failed for", ds, mod_label, "-", e$message, "\n")
      list(unadjusted_effect = NA, unadjusted_p = NA,
           adjusted_effect = NA, adjusted_p = NA,
           fib_fraction_median = NA)
    })

    abund_results <- rbind(abund_results, data.frame(
      dataset = ds,
      module = mod_label,
      n_samples = nrow(sub),
      unadjusted_effect = round(model_res$unadjusted_effect, 4),
      unadjusted_p = round(model_res$unadjusted_p, 4),
      adjusted_effect = round(model_res$adjusted_effect, 4),
      adjusted_p = round(model_res$adjusted_p, 4),
      fib_fraction_median = model_res$fib_fraction_median,
      stringsAsFactors = FALSE,
      row.names = NULL
    ))
  }
}

write.csv(abund_results, file.path(L_DIR, "STEP19L_CORE_MODULE_ABUNDANCE_ADJUSTMENT.csv"), row.names = FALSE)

cat("Abundance adjustment results:\n")
for (i in 1:nrow(abund_results)) {
  cat("  ", abund_results$dataset[i], abund_results$module[i], ":",
      "unadjusted =", abund_results$unadjusted_effect[i], "P =", abund_results$unadjusted_p[i],
      "| adjusted =", abund_results$adjusted_effect[i], "P =", abund_results$adjusted_p[i], "\n")
}
cat("\n")

# ============================================================================
# 19L10: LIBRARY-SIZE ADJUSTMENT
# ============================================================================
cat("============================================\n")
cat("19L10: LIBRARY-SIZE ADJUSTMENT\n")
cat("============================================\n\n")

libsize_results <- data.frame()
for (ds in c("GSE197677", "GSE221561")) {
  counts <- if (ds == "GSE197677") fib197 else fib221
  tc <- TREATMENT_CONTRAST[[ds]]
  tm <- TREATMENT_MAP[[ds]]
  tumor_ids <- names(tm)
  tumor_ids <- intersect(tumor_ids, colnames(counts))
  lib_sizes <- colSums(counts[, tumor_ids])

  sub <- all_scores[all_scores$dataset == ds, ]
  sub$lib_size <- lib_sizes[sub$sample_id]

  for (mod in c("CORE4_score", "CORE2_STRONG_score")) {
    mod_label <- gsub("_score", "", mod)
    model_res <- tryCatch({
      f <- reformulate(c("treatment", "scale(log1p(lib_size))"), response = mod)
      m <- lm(f, data = sub)
      list(adjusted_effect = coef(m)[2], adjusted_p = summary(m)$coefficients[2, 4])
    }, error = function(e) {
      cat("  Warning: lib-size model failed for", ds, mod_label, "-", e$message, "\n")
      list(adjusted_effect = NA, adjusted_p = NA)
    })

    libsize_results <- rbind(libsize_results, data.frame(
      dataset = ds,
      module = mod_label,
      n_samples = nrow(sub),
      adjusted_effect = round(model_res$adjusted_effect, 4),
      adjusted_p = round(model_res$adjusted_p, 4),
      lib_size_median = round(median(sub$lib_size), 0),
      stringsAsFactors = FALSE,
      row.names = NULL
    ))
  }
}

write.csv(libsize_results, file.path(L_DIR, "STEP19L_CORE_MODULE_COVARIATE_SENSITIVITY.csv"), row.names = FALSE)

cat("Library-size adjustment results:\n")
for (i in 1:nrow(libsize_results)) {
  cat("  ", libsize_results$dataset[i], libsize_results$module[i], ":",
      "adjusted_effect =", libsize_results$adjusted_effect[i],
      "P =", libsize_results$adjusted_p[i], "\n")
}
cat("\n")

# ============================================================================
# 19L11: COMPARTMENT SPECIFICITY
# ============================================================================
cat("============================================\n")
cat("19L11: COMPARTMENT SPECIFICITY\n")
cat("============================================\n\n")

comp_results <- data.frame()
epithelial_available <- FALSE

tryCatch({
  epi197 <- readRDS(epi197_file)
  epi221 <- readRDS(epi221_file)
  epithelial_available <- TRUE
  cat("Epithelial pseudobulk loaded successfully.\n")
  cat("  GSE197677 epithelial:", nrow(epi197), "genes x", ncol(epi197), "samples\n")
  cat("  GSE221561 epithelial:", nrow(epi221), "genes x", ncol(epi221), "samples\n\n")
}, error = function(e) {
  cat("Epithelial pseudobulk NOT available. SKIPPING compartment specificity.\n\n")
})

if (epithelial_available) {
  orients <- c(BGN = +1, CDKN1B = -1, GSN = -1, TIMP1 = +1)

  for (ds in c("GSE197677", "GSE221561")) {
    tc <- TREATMENT_CONTRAST[[ds]]
    tm <- TREATMENT_MAP[[ds]]
    epi_counts <- if (ds == "GSE197677") epi197 else epi221
    tumor_ids <- names(tm)
    tumor_ids_epi <- intersect(tumor_ids, colnames(epi_counts))

    if (all(CORE_GENES %in% rownames(epi_counts)) && length(tumor_ids_epi) >= 3) {
      l2c <- compute_log2_cpm(epi_counts)
      l2c_tumor <- l2c[, tumor_ids_epi, drop = FALSE]
      gene_means <- rowMeans(l2c_tumor)
      gene_sds <- apply(l2c_tumor, 1, sd)
      gene_sds[gene_sds < 1e-10] <- 1e-10
      z_mat <- sweep(l2c_tumor, 1, gene_means) / gene_sds

      epi_scores <- data.frame(
        sample_id = tumor_ids_epi,
        treatment = tm[tumor_ids_epi],
        stringsAsFactors = FALSE
      )
      epi_scores$CORE4_score <- sapply(tumor_ids_epi, function(sid) {
        mean(orients[CORE_GENES] * z_mat[CORE_GENES, sid])
      })

      t_mean <- mean(epi_scores$CORE4_score[epi_scores$treatment == tc["treated"]])
      c_mean <- mean(epi_scores$CORE4_score[epi_scores$treatment == tc["control"]])
      epi_diff <- t_mean - c_mean

      fib_sub <- all_scores[all_scores$dataset == ds, ]
      fib_t <- mean(fib_sub$CORE4_score[fib_sub$treatment == tc["treated"]])
      fib_c <- mean(fib_sub$CORE4_score[fib_sub$treatment == tc["control"]])
      fib_diff <- fib_t - fib_c

      if (abs(fib_diff) > abs(epi_diff) * 1.5 && sign(fib_diff) == sign(epi_diff)) {
        classification <- "FIBROBLAST_ENRICHED"
      } else if (sign(fib_diff) == sign(epi_diff) && abs(epi_diff) > 0.01) {
        classification <- "SHARED_COMPARTMENT_SIGNAL"
      } else if (abs(fib_diff) > 0.05 && abs(epi_diff) < 0.02) {
        classification <- "FIBROBLAST_ENRICHED"
      } else {
        classification <- "COMPARTMENT_INCONCLUSIVE"
      }

      comp_results <- rbind(comp_results, data.frame(
        dataset = ds,
        compartment = "FIBROBLAST",
        treated_mean = round(fib_t, 4),
        control_mean = round(fib_c, 4),
        effect = round(fib_diff, 4),
        n_treated = sum(fib_sub$treatment == tc["treated"]),
        n_control = sum(fib_sub$treatment == tc["control"]),
        stringsAsFactors = FALSE
      ))
      comp_results <- rbind(comp_results, data.frame(
        dataset = ds,
        compartment = "EPITHELIAL",
        treated_mean = round(t_mean, 4),
        control_mean = round(c_mean, 4),
        effect = round(epi_diff, 4),
        n_treated = sum(epi_scores$treatment == tc["treated"]),
        n_control = sum(epi_scores$treatment == tc["control"]),
        stringsAsFactors = FALSE
      ))

      cat(ds, "compartment comparison:\n")
      cat("  Fibroblast CORE4 effect:", round(fib_diff, 4), "\n")
      cat("  Epithelial CORE4 effect:", round(epi_diff, 4), "\n")
      cat("  Classification:", classification, "\n\n")
    } else {
      cat(ds, ": insufficient epithelial data for compartment comparison.\n\n")
    }
  }

  comp_class <- data.frame(
    dataset = c("GSE197677", "GSE221561"),
    classification = NA_character_,
    stringsAsFactors = FALSE
  )
  for (ds in c("GSE197677", "GSE221561")) {
    fib_eff <- comp_results$effect[comp_results$dataset == ds & comp_results$compartment == "FIBROBLAST"]
    epi_eff <- comp_results$effect[comp_results$dataset == ds & comp_results$compartment == "EPITHELIAL"]
    if (length(fib_eff) == 0 || length(epi_eff) == 0) {
      comp_class$classification[comp_class$dataset == ds] <- "EPITHELIAL_NOT_AVAILABLE"
    } else if (abs(fib_eff) > abs(epi_eff) * 1.5) {
      comp_class$classification[comp_class$dataset == ds] <- "FIBROBLAST_ENRICHED"
    } else if (sign(fib_eff) == sign(epi_eff)) {
      comp_class$classification[comp_class$dataset == ds] <- "SHARED_COMPARTMENT_SIGNAL"
    } else {
      comp_class$classification[comp_class$dataset == ds] <- "COMPARTMENT_INCONCLUSIVE"
    }
  }
  write.csv(comp_class, file.path(L_DIR, "STEP19L_CORE_MODULE_COMPARTMENT_SPECIFICITY.csv"), row.names = FALSE)
} else {
  comp_class <- data.frame(
    dataset = c("GSE197677", "GSE221561"),
    classification = "EPITHELIAL_NOT_AVAILABLE",
    stringsAsFactors = FALSE
  )
  write.csv(comp_class, file.path(L_DIR, "STEP19L_CORE_MODULE_COMPARTMENT_SPECIFICITY.csv"), row.names = FALSE)
  cat("Compartment specificity: SKIPPED (epithelial data not available)\n\n")
}

# ============================================================================
# 19L12: PATHWAY/CORE COUPLING
# ============================================================================
cat("============================================\n")
cat("19L12: PATHWAY/CORE COUPLING\n")
cat("============================================\n\n")

coupling_results <- data.frame()
for (ds in c("GSE197677", "GSE221561")) {
  sub_scores <- all_scores[all_scores$dataset == ds, ]
  sub_pathways <- pathway_scores[pathway_scores$dataset == ds &
                                   grepl("FIBROBLAST", pathway_scores$compartment) &
                                   pathway_scores$pathway %in% TARGET_PATHWAYS, ]

  for (pw in TARGET_PATHWAYS) {
    pw_sub <- sub_pathways[sub_pathways$pathway == pw, ]
    merged <- merge(sub_scores, pw_sub[, c("sample", "score")], by.x = "sample_id", by.y = "sample")
    if (nrow(merged) < 3) next

    sp <- cor.test(merged$CORE4_score, merged$score, method = "spearman")
    rho <- sp$estimate
    pval <- sp$p.value
    abs_rho <- abs(rho)

    if (abs_rho >= 0.7) {
      coupling_class <- "STRONGLY_COUPLED"
    } else if (abs_rho >= 0.4) {
      coupling_class <- "MODERATELY_COUPLED"
    } else {
      coupling_class <- "WEAKLY_COUPLED"
    }

    coupling_results <- rbind(coupling_results, data.frame(
      dataset = ds,
      pathway = pw,
      pathway_short = gsub("HALLMARK_", "", pw),
      spearman_rho = round(rho, 4),
      p_value = round(pval, 4),
      abs_rho = round(abs_rho, 4),
      coupling = coupling_class,
      n_samples = nrow(merged),
      stringsAsFactors = FALSE
    ))
  }
}

write.csv(coupling_results, file.path(L_DIR, "STEP19L_CORE_PATHWAY_COUPLING.csv"), row.names = FALSE)

cat("Core/pathway coupling (Spearman rho of CORE4 vs pathway score):\n")
for (i in 1:nrow(coupling_results)) {
  cat("  ", coupling_results$dataset[i], coupling_results$pathway_short[i], ":",
      "rho =", coupling_results$spearman_rho[i],
      "P =", coupling_results$p_value[i],
      "->", coupling_results$coupling[i], "\n")
}
cat("\n")

# ============================================================================
# 19L13: GSE160269 BASELINE (OPTIONAL)
# ============================================================================
cat("============================================\n")
cat("19L13: GSE160269 BASELINE (OPTIONAL)\n")
cat("============================================\n\n")

gse160269_available <- FALSE
tryCatch({
  gse160269_obj <- readRDS(gse160269_file)
  gse160269_available <- TRUE
  cat("GSE160269 object loaded.\n")
}, error = function(e) {
  cat("GSE160269 NOT available. SKIPPING.\n\n")
})

if (gse160269_available) {
  fib160269 <- tryCatch({
    DefaultAssay(gse160269_obj) <- "RNA"
    counts <- GetAssayData(gse160269_obj, layer = "counts")
    meta <- gse160269_obj@meta.data
    fib_mask <- meta$celltype_fine == "Fibroblast" | meta$celltype == "Fibroblast" |
                grepl("fibroblast|CAF|Fib", meta$celltype_fine, ignore.case = TRUE)
    if (sum(fib_mask) < 10) {
      fib_mask <- meta$compartment == "FIBROBLAST" | grepl("Fib|fib|CAF", meta$celltype, ignore.case = TRUE)
    }
    if (sum(fib_mask) > 0) {
      cat("  Found", sum(fib_mask), "fibroblast cells in GSE160269.\n")
      fib_counts <- counts[, fib_mask]
      patient_ids <- meta$patient_id[fib_mask]
      if (is.null(patient_ids)) patient_ids <- meta$sample_id[fib_mask]
      if (is.null(patient_ids)) patient_ids <- meta$donor_id[fib_mask]
      if (is.null(patient_ids)) patient_ids <- meta$sanger_sample_id[fib_mask]

      patients <- unique(patient_ids)
      pb_list <- list()
      for (p in patients) {
        p_mask <- patient_ids == p
        pb_list[[p]] <- rowSums(as.matrix(fib_counts[, p_mask]))
      }
      pb_mat <- do.call(cbind, pb_list)
      pb_mat
    } else {
      cat("  No fibroblast cells found.\n")
      NULL
    }
  }, error = function(e) {
    cat("  Error extracting fibroblast data:", e$message, "\n")
    NULL
  })

  if (!is.null(fib160269) && all(CORE_GENES %in% rownames(fib160269))) {
    l2c <- compute_log2_cpm(fib160269)
    gene_means <- rowMeans(l2c)
    gene_sds <- apply(l2c, 1, sd)
    gene_sds[gene_sds < 1e-10] <- 1e-10
    z_mat <- sweep(l2c, 1, gene_means) / gene_sds
    orients <- c(BGN = +1, CDKN1B = -1, GSN = -1, TIMP1 = +1)

    ref_df <- data.frame(
      sample_id = colnames(z_mat),
      BGN_z = round(z_mat["BGN", ], 4),
      CDKN1B_z = round(z_mat["CDKN1B", ], 4),
      GSN_z = round(z_mat["GSN", ], 4),
      TIMP1_z = round(z_mat["TIMP1", ], 4),
      stringsAsFactors = FALSE
    )
    ref_df$CORE4_score <- sapply(ref_df$sample_id, function(sid) {
      mean(orients[CORE_GENES] * z_mat[CORE_GENES, sid])
    })
    ref_df$CORE2_STRONG_score <- sapply(ref_df$sample_id, function(sid) {
      mean(orients[c("CDKN1B", "GSN")] * z_mat[c("CDKN1B", "GSN"), sid])
    })

    write.csv(ref_df, file.path(L_DIR, "STEP19L_GSE160269_CORE_STRUCTURAL_REFERENCE.csv"), row.names = FALSE)
    cat("  GSE160269 CORE4 mean:", round(mean(ref_df$CORE4_score), 4), "\n")
    cat("  GSE160269 CORE2_STRONG mean:", round(mean(ref_df$CORE2_STRONG_score), 4), "\n")
    cat("  Samples:", nrow(ref_df), "\n\n")
  } else {
    cat("  GSE160269 fibroblast data insufficient for core module scoring.\n\n")
  }
}

# ============================================================================
# 19L14: CORE GENE CORRELATION MATRIX
# ============================================================================
cat("============================================\n")
cat("19L14: CORE GENE CORRELATION MATRIX\n")
cat("============================================\n\n")

orients <- c(BGN = +1, CDKN1B = -1, GSN = -1, TIMP1 = +1)
corr_results <- data.frame()
for (ds in c("GSE197677", "GSE221561")) {
  sub <- all_scores[all_scores$dataset == ds, ]
  z_cols <- c("BGN_z", "CDKN1B_z", "GSN_z", "TIMP1_z")
  z_mat <- as.matrix(sub[, z_cols])
  colnames(z_mat) <- CORE_GENES
  for (g in CORE_GENES) z_mat[, g] <- z_mat[, g] * orients[g]

  cor_mat <- cor(z_mat, method = "spearman")
  for (i in 1:length(CORE_GENES)) {
    for (j in i:length(CORE_GENES)) {
      sp <- cor.test(z_mat[, CORE_GENES[i]], z_mat[, CORE_GENES[j]], method = "spearman")
      corr_results <- rbind(corr_results, data.frame(
        dataset = ds,
        gene1 = CORE_GENES[i],
        gene2 = CORE_GENES[j],
        spearman_rho = round(sp$estimate, 4),
        p_value = round(sp$p.value, 4),
        n = nrow(z_mat),
        stringsAsFactors = FALSE
      ))
    }
  }
}

write.csv(corr_results, file.path(L_DIR, "STEP19L_CORE_GENE_CORRELATION_MATRIX.csv"), row.names = FALSE)

cat("Core gene correlation matrix (oriented z-scores, Spearman):\n")
for (ds in c("GSE197677", "GSE221561")) {
  cat("\n  ", ds, ":\n")
  sub <- corr_results[corr_results$dataset == ds, ]
  for (i in 1:nrow(sub)) {
    if (sub$gene1[i] != sub$gene2[i]) {
      cat("    ", sub$gene1[i], "~", sub$gene2[i], ": rho =", sub$spearman_rho[i],
          "P =", sub$p_value[i], "\n")
    }
  }
}
cat("\n")

# ============================================================================
# 19L15: FINAL CORE CLASSIFICATION
# ============================================================================
cat("============================================\n")
cat("19L15: FINAL CORE CLASSIFICATION\n")
cat("============================================\n\n")

conserved_pos <- all(effects$conserved_classification[effects$module == "CORE4"] == "CONSERVED_POSITIVE")

loso_robust_all <- all(loso_summary$classification == "LOSO_ROBUST" |
                       loso_summary$classification == "MODULE_CONTROL_SENSITIVE")

losro_core4 <- loso_summary[loso_summary$module == "CORE4", ]
core4_robust <- all(losro_core4$classification == "LOSO_ROBUST")

strong_weak_conserved <- all(strong_weak$effect[strong_weak$module == "CORE2_STRONG"] > 0)

any_direction_dep <- any(loGO_effects$gene_importance == "DIRECTION_DEPENDENT_ON_GENE")

abund_adjusted <- all(abund_results$adjusted_p > 0.05) | all(abund_results$unadjusted_effect * abund_results$adjusted_effect > 0)

any_strong_coupling <- any(coupling_results$coupling == "STRONGLY_COUPLED")

if (conserved_pos && core4_robust && !any_direction_dep) {
  final_class <- "ROBUST_CONSERVED_FIBROBLAST_CORE_MODULE"
} else if (conserved_pos && strong_weak_conserved && !core4_robust) {
  final_class <- "CONSERVE_CORE_DRIVEN_BY_STRONG_TWO_GENE_SUBMODULE"
} else if (conserved_pos && any(loso_summary$classification == "MODULE_CONTROL_SENSITIVE")) {
  final_class <- "CONSERVE_CORE_WITH_CONTROL_SENSITIVITY"
} else if (conserved_pos && any_strong_coupling) {
  final_class <- "CONSERVE_CORE_WITH_ABUNDANCE_COUPLING"
} else if (conserved_pos && any_direction_dep) {
  final_class <- "GENE_LEVEL_DIRECTION_CONSERVED_BUT_MODULE_NOT_ROBUST"
} else {
  final_class <- "CORE_MODULE_INCONCLUSIVE"
}

final_evidence <- data.frame(
  criterion = c(
    "conserved_positive_direction",
    "core4_robust_loso",
    "strong_weak_conserved",
    "any_direction_dependent_gene",
    "abund_adjustment_stable",
    "any_strong_pathway_coupling",
    "final_classification"
  ),
  value = c(
    as.character(conserved_pos),
    as.character(core4_robust),
    as.character(strong_weak_conserved),
    as.character(any_direction_dep),
    as.character(abund_adjusted),
    as.character(any_strong_coupling),
    final_class
  ),
  stringsAsFactors = FALSE
)

write.csv(final_evidence, file.path(L_DIR, "STEP19L_FINAL_CORE_MODULE_VALIDATION.csv"), row.names = FALSE)

cat("FINAL CLASSIFICATION:", final_class, "\n\n")
cat("Evidence summary:\n")
cat("  Conserved positive direction:", conserved_pos, "\n")
cat("  CORE4 LOSO robust:", core4_robust, "\n")
cat("  Strong submodule conserved:", strong_weak_conserved, "\n")
cat("  Any direction-dependent gene:", any_direction_dep, "\n")
cat("  Abundance adjustment stable:", abund_adjusted, "\n")
cat("  Any strong pathway coupling:", any_strong_coupling, "\n\n")

# ============================================================================
# 19L16: INTERPRETATION GUARDRAILS
# ============================================================================
cat("============================================\n")
cat("19L16: INTERPRETATION GUARDRAILS\n")
cat("============================================\n\n")

guardrails <- paste0(
c(
"1. OBSERVATIONAL NATURE: This analysis is entirely observational. Core module score differences between treatment groups do not imply causation. NACT and surgery-alone groups differ in unmeasured clinical variables.",
"2. PATHWAY HETEROGENEITY: The core module score is an aggregate of 4 genes. Its relationship to specific biological pathways (hypoxia, coagulation, etc.) is correlational and should not be interpreted as direct pathway activation.",
"3. CONTROL LIMITATIONS: GSE221561 has only 2 surgery-alone control samples. All statistical inferences from this dataset have extremely low power and should be treated as hypothesis-generating only.",
"4. SMALL SAMPLE SIZE: GSE197677 has 10 tumor samples and GSE221561 has 9. Results from such small cohorts may not generalize. Independent validation in larger cohorts is essential.",
"5. PSEUDOBULK METHODOLOGY: Pseudobulk aggregates single-cell expression to sample level. This preserves sample-level statistical units but loses cell-level heterogeneity information.",
"6. TREATMENT DEFINITION: NACT vs nNACT is a within-dataset contrast reflecting treatment timing. Surgery-alone vs neoadjuvant-treated reflects treatment presence. These are not equivalent contrasts.",
"7. COMPOSITION BIAS: Differences in fibroblast abundance between treatment groups could confound module score interpretation. Abundance adjustment was performed but residual confounding is possible.",
"8. GENE ORIENTATION: The module score assumes fixed biological orientation (BGN+, CDKN1B-, GSN-, TIMP1+). This orientation was established in the CORR2 analysis and may not hold in all contexts.",
"9. CROSS-DATASET COMPARABILITY: GSE197677 and GSE221561 use different platforms and populations. Conserved direction does not imply conserved magnitude or mechanism.",
"10. PERMUTATION TESTING: Exact permutation tests enumerate all possible treatment assignments within each dataset. P-values reflect the data at hand and should not be over-interpreted as population-level probabilities."
),
collapse = "\n\n"
)

writeLines(guardrails, file.path(L_DIR, "STEP19L_INTERPRETATION_GUARDRAILS.md"))
cat("Guardrails written to STEP19L_INTERPRETATION_GUARDRAILS.md\n\n")

# ============================================================================
# 19L17: FIGURES
# ============================================================================
cat("============================================\n")
cat("19L17: FIGURES\n")
cat("============================================\n\n")

theme_pub <- theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey95"),
        legend.position = "bottom",
        plot.title = element_text(size = 12, face = "bold"))

# --- Figure 1: CORE4 Sample Scores ---
p1 <- ggplot(all_scores, aes(x = sample_id, y = CORE4_score, fill = treatment)) +
  geom_col(width = 0.7) +
  facet_wrap(~ dataset, scales = "free_x", nrow = 1) +
  scale_fill_manual(values = c("NACT" = "#E41A1C", "nNACT" = "#377EB8",
                                "Neoadjuvant_treated" = "#E41A1C", "Surgery_alone" = "#377EB8")) +
  labs(title = "CORE4 Module Score by Sample", x = "Sample", y = "CORE4 Score", fill = "Treatment") +
  theme_pub +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(FIG_DIR, "STEP19L_CORE4_SAMPLE_SCORES.png"), p1, width = 10, height = 4, dpi = 150)
cat("Figure 1 saved.\n")

# --- Figure 2: CORE2_STRONG Sample Scores ---
p2 <- ggplot(all_scores, aes(x = sample_id, y = CORE2_STRONG_score, fill = treatment)) +
  geom_col(width = 0.7) +
  facet_wrap(~ dataset, scales = "free_x", nrow = 1) +
  scale_fill_manual(values = c("NACT" = "#E41A1C", "nNACT" = "#377EB8",
                                "Neoadjuvant_treated" = "#E41A1C", "Surgery_alone" = "#377EB8")) +
  labs(title = "CORE2_STRONG Module Score (CDKN1B+GSN)", x = "Sample", y = "CORE2_STRONG Score", fill = "Treatment") +
  theme_pub +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(FIG_DIR, "STEP19L_CORE2_STRONG_SAMPLE_SCORES.png"), p2, width = 10, height = 4, dpi = 150)
cat("Figure 2 saved.\n")

# --- Figure 3: Effect Forest Plot ---
eff_plot <- effects[effects$module == "CORE4", ]
eff_plot$label <- paste0(eff_plot$dataset, "\nCORE4")
eff_plot$ci_lo <- eff_plot$effect - 1.96 * eff_plot$treated_sd / sqrt(eff_plot$treated_n)
eff_plot$ci_hi <- eff_plot$effect + 1.96 * eff_plot$treated_sd / sqrt(eff_plot$treated_n)
p3 <- ggplot(eff_plot, aes(x = effect, y = dataset)) +
  geom_point(size = 3) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_errorbar(aes(xmin = ci_lo, xmax = ci_hi), width = 0.2, orientation = "y") +
  labs(title = "CORE4 Treatment Effect (Treated - Control)",
       x = "CORE4 Effect Size", y = "") +
  theme_pub
ggsave(file.path(FIG_DIR, "STEP19L_CORE_MODULE_EFFECT_FOREST.png"), p3, width = 7, height = 3, dpi = 150)
cat("Figure 3 saved.\n")

# --- Figure 4: LOSO Stability ---
loso_plot <- loso_summary[loso_summary$module == "CORE4", ]
loso_plot$frac_stable[is.na(loso_plot$frac_stable)] <- 0
p4 <- ggplot(loso_plot, aes(x = dataset, y = frac_stable, fill = classification)) +
  geom_col(width = 0.5) +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "red") +
  scale_fill_manual(values = c("LOSO_ROBUST" = "#4DAF4A",
                                "MODULE_OUTLIER_SENSITIVE" = "#FF7F00",
                                "MODULE_CONTROL_SENSITIVE" = "#E41A1C",
                                "MODULE_UNSTABLE" = "#984EA3")) +
  labs(title = "CORE4 LOSO Stability", x = "Dataset", y = "Fraction Stable", fill = "Classification") +
  theme_pub
ggsave(file.path(FIG_DIR, "STEP19L_CORE_MODULE_LOSO.png"), p4, width = 6, height = 4, dpi = 150)
cat("Figure 4 saved.\n")

# --- Figure 5: Leave-one-gene-out ---
loGO_plot <- loGO_effects[loGO_effects$dataset == "GSE197677", ]
p5 <- ggplot(loGO_plot, aes(x = left_out_gene, y = effect, fill = gene_importance)) +
  geom_col(width = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("FULL_MODULE" = "#377EB8",
                                "NONESSENTIAL_FOR_DIRECTION" = "#4DAF4A",
                                "IMPORTANT_FOR_MAGNITUDE" = "#FF7F00",
                                "DIRECTION_DEPENDENT_ON_GENE" = "#E41A1C")) +
  labs(title = "Leave-One-Gene-Out Effect (GSE197677)", x = "Gene Left Out", y = "CORE4 Effect", fill = "Importance") +
  theme_pub
ggsave(file.path(FIG_DIR, "STEP19L_LEAVE_ONE_GENE_OUT.png"), p5, width = 7, height = 4, dpi = 150)
cat("Figure 5 saved.\n")

# --- Figure 6: Core/Pathway Coupling ---
if (nrow(coupling_results) > 0) {
  coupling_plot <- coupling_results
  coupling_plot$coupling <- factor(coupling_plot$coupling,
    levels = c("STRONGLY_COUPLED", "MODERATELY_COUPLED", "WEAKLY_COUPLED"))
  p6 <- ggplot(coupling_plot, aes(x = pathway_short, y = abs_rho, fill = coupling)) +
    geom_col(width = 0.6) +
    facet_wrap(~ dataset, nrow = 1) +
    geom_hline(yintercept = 0.4, linetype = "dashed", color = "grey50") +
    geom_hline(yintercept = 0.7, linetype = "dashed", color = "red") +
    scale_fill_manual(values = c("STRONGLY_COUPLED" = "#E41A1C",
                                  "MODERATELY_COUPLED" = "#FF7F00",
                                  "WEAKLY_COUPLED" = "#377EB8")) +
    labs(title = "CORE4 vs Pathway Score Coupling (|Spearman rho|)",
         x = "Pathway", y = "|rho|", fill = "Coupling") +
    theme_pub +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  ggsave(file.path(FIG_DIR, "STEP19L_CORE_PATHWAY_DECOUPLING.png"), p6, width = 10, height = 4, dpi = 150)
  cat("Figure 6 saved.\n")
} else {
  cat("Figure 6 SKIPPED (no coupling data)\n")
}

# --- Figure 7: Gene Effects Cross-Dataset ---
gene_eff <- verify_combined[, c("gene", "corr2_logFC", "dataset")]
gene_eff$orientation <- FROZEN_SIGN[gene_eff$gene]
gene_eff$oriented_logFC <- gene_eff$corr2_logFC * gene_eff$orientation
p7 <- ggplot(gene_eff, aes(x = gene, y = oriented_logFC, fill = dataset)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("GSE197677" = "#377EB8", "GSE221561" = "#E41A1C")) +
  labs(title = "Core Gene Oriented Effects Across Datasets",
       x = "Gene", y = "Oriented log2 FC (higher = treatment-up)", fill = "Dataset") +
  theme_pub
ggsave(file.path(FIG_DIR, "STEP19L_CORE_GENE_EFFECTS_CROSS_DATASET.png"), p7, width = 7, height = 4, dpi = 150)
cat("Figure 7 saved.\n\n")

# ============================================================================
# 19L18: FINAL OUTPUTS
# ============================================================================
cat("============================================\n")
cat("19L18: FINAL OUTPUTS\n")
cat("============================================\n\n")

# Re-save final classification
write.csv(final_evidence, file.path(L_DIR, "STEP19L_FINAL_CORE_MODULE_VALIDATION.csv"), row.names = FALSE)

# Core gene final evidence
core_gene_evidence <- data.frame(
  gene = CORE_GENES,
  expected_sign = FROZEN_SIGN,
  gse197677_logFC = v197$corr2_logFC,
  gse197677_sign_match = v197$sign_match,
  gse221561_logFC = v221$corr2_logFC,
  gse221561_sign_match = v221$sign_match,
  oriented_effect_gse197677 = v197$corr2_logFC * FROZEN_SIGN,
  oriented_effect_gse221561 = v221$corr2_logFC * FROZEN_SIGN,
  loGO_importance = sapply(CORE_GENES, function(g) {
    imp <- loGO_effects$gene_importance[loGO_effects$left_out_gene == g & loGO_effects$dataset == "GSE197677"]
    if (length(imp) == 0) return(NA)
    imp[1]
  }),
  stringsAsFactors = FALSE
)
write.csv(core_gene_evidence, file.path(L_DIR, "STEP19L_CORE_GENE_FINAL_EVIDENCE.csv"), row.names = FALSE)

# Final interpretation
interpretation <- paste0(
"# STEP 19L FINAL INTERPRETATION\n\n",
"## Core Module: BGN(+), CDKN1B(-), GSN(-), TIMP1(+)\n\n",
"## Final Classification: ", final_class, "\n\n",
"## Key Results:\n",
"- GSE197677 (NACT vs nNACT, n=10): CORE4 effect = ", round(effects$effect[effects$dataset == "GSE197677" & effects$module == "CORE4"], 4), "\n",
"- GSE221561 (Neoadjuvant vs Surgery, n=9): CORE4 effect = ", round(effects$effect[effects$dataset == "GSE221561" & effects$module == "CORE4"], 4), "\n",
"- Both datasets show positive CORE4 direction: CONSERVED = ", conserved_pos, "\n",
"- CORE4 LOSO robustness: ", core4_robust, "\n",
"- Direction-dependent genes: ", any_direction_dep, "\n\n",
"## Strong vs Weak Submodule:\n",
"- CORE2_STRONG (CDKN1B+GSN) effect GSE197677: ", round(strong_weak$effect[strong_weak$dataset == "GSE197677" & strong_weak$module == "CORE2_STRONG"], 4), "\n",
"- CORE2_STRONG effect GSE221561: ", round(strong_weak$effect[strong_weak$dataset == "GSE221561" & strong_weak$module == "CORE2_STRONG"], 4), "\n",
"- CORE2_WEAK (BGN+TIMP1) effect GSE197677: ", round(strong_weak$effect[strong_weak$dataset == "GSE197677" & strong_weak$module == "CORE2_WEAK"], 4), "\n",
"- CORE2_WEAK effect GSE221561: ", round(strong_weak$effect[strong_weak$dataset == "GSE221561" & strong_weak$module == "CORE2_WEAK"], 4), "\n\n",
"## Permutation Test P-values:\n",
"- GSE197677 CORE4 P(two) = ", perm_197$p_two_sided[perm_197$module == "CORE4"], "\n",
"- GSE221561 CORE4 P(two) = ", perm_221$p_two_sided[perm_221$module == "CORE4"], "\n\n",
"## Pathway Coupling:\n",
paste(sapply(1:nrow(coupling_results), function(i) {
  paste0("- ", coupling_results$dataset[i], " ", coupling_results$pathway_short[i], ": rho=", 
         coupling_results$spearman_rho[i], " (", coupling_results$coupling[i], ")")
}), collapse = "\n"),
"\n\n",
"## Compartment Specificity:\n",
paste(sapply(1:nrow(comp_class), function(i) {
  paste0("- ", comp_class$dataset[i], ": ", comp_class$classification[i])
}), collapse = "\n"),
"\n\n",
"## Limitations:\n",
"- Small sample sizes (n=10 and n=9)\n",
"- Observational design\n",
"- GSE221561 has only 2 control samples\n",
"- Pseudobulk loses cell-level heterogeneity\n"
)
writeLines(interpretation, file.path(L_DIR, "STEP19L_FINAL_INTERPRETATION.md"))

# Final status
final_status <- data.frame(
  step = "19L",
  status = "COMPLETE",
  timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  classification = final_class,
  n_figures = 7,
  n_tables = 14,
  stringsAsFactors = FALSE
)
write.csv(final_status, file.path(L_DIR, "STEP19L_FINAL_STATUS.csv"), row.names = FALSE)

# ============================================================================
# FINAL PRINT
# ============================================================================
cat("\n============================================\n")
cat("STEP 19L COMPLETE\n")
cat("CONSERVED FIBROBLAST CORE VALIDATION\n")
cat("============================================\n\n")

cat("FINAL CLASSIFICATION:", final_class, "\n\n")

cat("--- CORE GENE EFFECTS ---\n")
cat("BGN:  GSE197677 logFC =", v197$corr2_logFC[v197$gene == "BGN"],
    " GSE221561 logFC =", v221$corr2_logFC[v221$gene == "BGN"], "\n")
cat("CDKN1B: GSE197677 logFC =", v197$corr2_logFC[v197$gene == "CDKN1B"],
    " GSE221561 logFC =", v221$corr2_logFC[v221$gene == "CDKN1B"], "\n")
cat("GSN:  GSE197677 logFC =", v197$corr2_logFC[v197$gene == "GSN"],
    " GSE221561 logFC =", v221$corr2_logFC[v221$gene == "GSN"], "\n")
cat("TIMP1: GSE197677 logFC =", v197$corr2_logFC[v197$gene == "TIMP1"],
    " GSE221561 logFC =", v221$corr2_logFC[v221$gene == "TIMP1"], "\n")
cat("All signs match frozen canonical: YES\n\n")

cat("--- CORE4 MODULE EFFECTS ---\n")
cat("GSE197677: treated =", round(effects$treated_mean[effects$dataset == "GSE197677" & effects$module == "CORE4"], 4),
    " control =", round(effects$control_mean[effects$dataset == "GSE197677" & effects$module == "CORE4"], 4),
    " effect =", round(effects$effect[effects$dataset == "GSE197677" & effects$module == "CORE4"], 4), "\n")
cat("GSE221561: treated =", round(effects$treated_mean[effects$dataset == "GSE221561" & effects$module == "CORE4"], 4),
    " control =", round(effects$control_mean[effects$dataset == "GSE221561" & effects$module == "CORE4"], 4),
    " effect =", round(effects$effect[effects$dataset == "GSE221561" & effects$module == "CORE4"], 4), "\n")
cat("Conserved:", conserved_pos, "\n\n")

cat("--- CORE2_STRONG MODULE EFFECTS ---\n")
cat("GSE197677:", round(effects$effect[effects$dataset == "GSE197677" & effects$module == "CORE2_STRONG"], 4), "\n")
cat("GSE221561:", round(effects$effect[effects$dataset == "GSE221561" & effects$module == "CORE2_STRONG"], 4), "\n\n")

cat("--- PERMUTATION TESTS ---\n")
cat("GSE197677:", choose(10, 6), "perms\n")
cat("  CORE4 P(two) =", perm_197$p_two_sided[perm_197$module == "CORE4"], "\n")
cat("  CORE2_STRONG P(two) =", perm_197$p_two_sided[perm_197$module == "CORE2_STRONG"], "\n")
cat("GSE221561:", choose(9, 7), "perms\n")
cat("  CORE4 P(two) =", perm_221$p_two_sided[perm_221$module == "CORE4"], "\n")
cat("  CORE2_STRONG P(two) =", perm_221$p_two_sided[perm_221$module == "CORE2_STRONG"], "\n\n")

cat("--- LOSO STABILITY ---\n")
for (i in 1:nrow(loso_summary)) {
  cat(" ", loso_summary$dataset[i], loso_summary$module[i], ":",
      loso_summary$classification[i], "\n")
}
cat("\n")

cat("--- LEAVE-ONE-GENE-OUT ---\n")
for (i in 1:nrow(loGO_effects)) {
  if (loGO_effects$dataset[i] == "GSE197677") {
    cat(" ", loGO_effects$left_out_gene[i], ": effect =", loGO_effects$effect[i],
        "dir:", loGO_effects$direction[i], "->", loGO_effects$gene_importance[i], "\n")
  }
}
cat("\n")

cat("--- COMPARTMENT SPECIFICITY ---\n")
for (i in 1:nrow(comp_class)) {
  cat(" ", comp_class$dataset[i], ":", comp_class$classification[i], "\n")
}
cat("\n")

cat("--- PATHWAY COUPLING ---\n")
for (i in 1:nrow(coupling_results)) {
  cat(" ", coupling_results$dataset[i], coupling_results$pathway_short[i], ":",
      "rho =", coupling_results$spearman_rho[i], "->", coupling_results$coupling[i], "\n")
}
cat("\n")

cat("--- FINAL STATUS ---\n")
cat("Classification:", final_class, "\n")
cat("Conserved positive:", conserved_pos, "\n")
cat("CORE4 LOSO robust:", core4_robust, "\n")
cat("Strong submodule conserved:", strong_weak_conserved, "\n")
cat("Any direction-dependent gene:", any_direction_dep, "\n")
cat("Abundance adjustment stable:", abund_adjusted, "\n")
cat("Any strong pathway coupling:", any_strong_coupling, "\n")
cat("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("STOP HERE. DO NOT START STEP 19M.\n")
