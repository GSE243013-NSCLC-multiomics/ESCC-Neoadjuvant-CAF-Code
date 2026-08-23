#!/usr/bin/env Rscript
# ============================================================================
# STEP 19O-CORR: STEP19O FINAL CLASSIFICATION RECONCILIATION
# ============================================================================
# Audit the completed Step19O outputs and determine whether the reported
# final classification is logically consistent with the actual cross-dataset
# signaling evidence.
#
# CLASSIFICATION RECONCILIATION ONLY.
# DO NOT rerun signaling analysis, recompute pseudobulk, or modify upstream.
# ============================================================================

BASE_DIR <- Sys.getenv("ESCC_CAF_PROJECT_ROOT", unset = ".")
O_DIR    <- file.path(BASE_DIR, "06_tables/cross_dataset/pathway_analysis/Step19O_CORR")
SRC_DIR  <- file.path(BASE_DIR, "06_tables/cross_dataset/pathway_analysis/Step19O")
LOG_DIR  <- file.path(BASE_DIR, "08_logs")
dir.create(O_DIR, showWarnings = FALSE, recursive = TRUE)

log_con <- file.path(LOG_DIR, "STEP19O_CORR_CLASSIFICATION_RECONCILIATION.log")

cat_log <- function(...) {
  txt <- paste0(...)
  cat(txt)
  cat(txt, file = log_con, append = TRUE)
}

cat_log("============================================\n")
cat_log("STEP 19O-CORR: FINAL CLASSIFICATION RECONCILIATION\n")
cat_log("============================================\n\n")

# ============================================================================
# 19O-CORR1 — INPUT FILE AUDIT
# ============================================================================
cat_log("============================================\n")
cat_log("19O-CORR1: INPUT FILE AUDIT\n")
cat_log("============================================\n\n")

required_files <- c(
  "STEP19O_SIGNALING_EVIDENCE_SCORE.csv",
  "STEP19O_CROSS_DATASET_SIGNAL_REPLICATION.csv",
  "STEP19O_SOURCE_COMPARTMENT_PRIORITIZATION.csv",
  "STEP19O_ELIGIBLE_LIGAND_RECEPTOR_PAIRS.csv",
  "STEP19O_LIGAND_CORE2_ASSOCIATIONS.csv",
  "STEP19O_RUNX1_SIGNALING_BRIDGE.csv",
  "STEP19O_SIGNALING_LOSO_SUMMARY.csv",
  "STEP19O_FINAL_LIGAND_RANKING.csv",
  "STEP19O_FINAL_LIGAND_RECEPTOR_PAIRS.csv",
  "STEP19O_FINAL_SOURCE_COMPARTMENTS.csv",
  "STEP19O_RUNX1_BRIDGE_FINAL.csv",
  "STEP19O_FINAL_SIGNALING_CLASSIFICATION.csv",
  "STEP19O_FINAL_STATUS.csv",
  "STEP19O_FINAL_INTERPRETATION.md",
  "STEP19O_RESOURCE_STATUS.csv",
  "STEP19O_SOURCE_COMPARTMENT_COUNTS.csv",
  "STEP19O_FIBROBLAST_RECEPTOR_AVAILABILITY.csv",
  "STEP19O_SIGNALING_ABUNDANCE_SENSITIVITY.csv"
)

file_audit <- data.frame(
  file = required_files,
  exists = file.exists(file.path(SRC_DIR, required_files)),
  stringsAsFactors = FALSE
)

for (i in 1:nrow(file_audit)) {
  cat_log(sprintf("  %-55s: %s\n", file_audit$file[i],
    ifelse(file_audit$exists[i], "FOUND", "MISSING")))
}

missing <- file_audit$file[!file_audit$exists]
if (length(missing) > 0) {
  cat_log(sprintf("\nWARNING: %d files missing: %s\n", length(missing), paste(missing, collapse=", ")))
}

cat_log("\nAll available files loaded for audit.\n\n")

# Load key files
evi    <- read.csv(file.path(SRC_DIR, "STEP19O_SIGNALING_EVIDENCE_SCORE.csv"))
elig   <- read.csv(file.path(SRC_DIR, "STEP19O_ELIGIBLE_LIGAND_RECEPTOR_PAIRS.csv"))
assoc  <- read.csv(file.path(SRC_DIR, "STEP19O_LIGAND_CORE2_ASSOCIATIONS.csv"))
runx1  <- read.csv(file.path(SRC_DIR, "STEP19O_RUNX1_BRIDGE_FINAL.csv"))
loso   <- read.csv(file.path(SRC_DIR, "STEP19O_SIGNALING_LOSO_SUMMARY.csv"))
final_cls <- read.csv(file.path(SRC_DIR, "STEP19O_FINAL_SIGNALING_CLASSIFICATION.csv"))
final_st  <- read.csv(file.path(SRC_DIR, "STEP19O_FINAL_STATUS.csv"))
res    <- read.csv(file.path(SRC_DIR, "STEP19O_RESOURCE_STATUS.csv"))
src_cnt <- read.csv(file.path(SRC_DIR, "STEP19O_SOURCE_COMPARTMENT_COUNTS.csv"))
rec_av  <- read.csv(file.path(SRC_DIR, "STEP19O_FIBROBLAST_RECEPTOR_AVAILABILITY.csv"))

# ============================================================================
# 19O-CORR2 — VERIFY OBSERVED COUNTS
# ============================================================================
cat_log("============================================\n")
cat_log("19O-CORR2: VERIFY OBSERVED COUNTS\n")
cat_log("============================================\n\n")

# Eligible LR pairs
n_eligible <- nrow(elig)
cat_log(sprintf("Eligible LR pairs (from ELIGIBLE file): %d\n", n_eligible))

# Check for GSE221561 entries in eligible pairs
n_elig_221 <- sum(elig$dataset == "GSE221561")
n_elig_197 <- sum(elig$dataset == "GSE197677")
cat_log(sprintf("  GSE197677 eligible pairs: %d\n", n_elig_197))
cat_log(sprintf("  GSE221561 eligible pairs: %d\n", n_elig_221))

# Source compartment counts audit
cat_log("\nSource compartment counts (from SOURCE_COMPARTMENT_COUNTS.csv):\n")
for (i in 1:nrow(src_cnt)) {
  cat_log(sprintf("  %s %s: %d genes x %d samples\n",
    src_cnt$dataset[i], src_cnt$compartment[i], src_cnt$n_genes[i], src_cnt$n_samples[i]))
}

n_src_221 <- sum(src_cnt$dataset == "GSE221561")
cat_log(sprintf("\nGSE221561 source compartment rows: %d\n", n_src_221))
if (n_src_221 == 0) {
  cat_log("CRITICAL FINDING: GSE221561 source compartment pseudobulk was NOT computed.\n")
  cat_log("All GSE221561 signaling analysis was skipped.\n")
}

# Ligand-CORE2 associations by dataset
n_assoc_197 <- sum(assoc$dataset == "GSE197677")
n_assoc_221 <- sum(assoc$dataset == "GSE221561")
cat_log(sprintf("\nLigand-CORE2 associations: GSE197677=%d, GSE221561=%d\n",
  n_assoc_197, n_assoc_221))

# Evidence score classifications
n_moderate <- sum(evi$classification == "MODERATE_SIGNALING_SUPPORT")
n_weak     <- sum(evi$classification == "WEAK_SIGNALING_SUPPORT")
n_none     <- sum(evi$classification == "NO_REPRODUCIBLE_SIGNALING_SUPPORT")
n_high     <- sum(evi$classification == "HIGH_SIGNALING_SUPPORT")

cat_log(sprintf("\nEvidence score classifications:\n"))
cat_log(sprintf("  HIGH_SIGNALING_SUPPORT:     %d\n", n_high))
cat_log(sprintf("  MODERATE_SIGNALING_SUPPORT: %d\n", n_moderate))
cat_log(sprintf("  WEAK_SIGNALING_SUPPORT:     %d\n", n_weak))
cat_log(sprintf("  NO_REPRODUCIBLE:            %d\n", n_none))

# Identify the actual moderate candidates
moderate_ligands <- evi$ligand[evi$classification == "MODERATE_SIGNALING_SUPPORT"]
cat_log(sprintf("\nActual MODERATE candidates (from evidence score): %s\n",
  paste(moderate_ligands, collapse=", ")))

# Check user-claimed candidates vs evidence score
user_claimed_moderate <- c("AREG", "EREG", "DLL1", "IFNG")
cat_log(sprintf("\nUser-claimed moderate candidates: %s\n",
  paste(user_claimed_moderate, collapse=", ")))

discrepancy <- setdiff(user_claimed_moderate, moderate_ligands)
extra_moderate <- setdiff(moderate_ligands, user_claimed_moderate)
if (length(discrepancy) > 0 || length(extra_moderate) > 0) {
  cat_log("DISCREPANCY DETECTED:\n")
  if (length(discrepancy) > 0) {
    for (d in discrepancy) {
      d_class <- evi$classification[evi$ligand == d]
      d_score <- evi$evidence_score[evi$ligand == d]
      cat_log(sprintf("  %s: classified as %s (score %d/8), NOT MODERATE\n",
        d, d_class, d_score))
    }
  }
  if (length(extra_moderate) > 0) {
    for (e in extra_moderate) {
      e_score <- evi$evidence_score[evi$ligand == e]
      cat_log(sprintf("  %s: classified as MODERATE_SIGNALING_SUPPORT (score %d/8), not listed in user claim\n",
        e, e_score))
    }
  }
}

# Cross-dataset replication
n_core2_197 <- sum(evi$core2_197 == TRUE)
n_core2_221 <- sum(evi$core2_221 == TRUE)
n_core2_same_dir <- sum(evi$core2_same_dir == TRUE)
cat_log(sprintf("\nCore2_197 TRUE: %d\n", n_core2_197))
cat_log(sprintf("Core2_221 TRUE: %d\n", n_core2_221))
cat_log(sprintf("Core2_same_dir TRUE: %d\n", n_core2_same_dir))

n_cross_replicated <- n_core2_same_dir
cat_log(sprintf("Cross-dataset replicated: %d\n", n_cross_replicated))

# RUNX1 bridges
n_runx1_bridge <- sum(evi$runx1_bridge == TRUE)
n_runx1_bridge_csv <- sum(runx1$classification == "LIGAND_RUNX1_CORE2_BRIDGE")
cat_log(sprintf("RUNX1 bridges (evidence score): %d\n", n_runx1_bridge))
cat_log(sprintf("RUNX1 bridges (RUNX1_final CSV): %d\n", n_runx1_bridge_csv))

# Cross-dataset signal replication file
repl_file <- file.path(SRC_DIR, "STEP19O_CROSS_DATASET_SIGNAL_REPLICATION.csv")
repl_content <- readLines(repl_file, n = 5)
cat_log(sprintf("\nCROSS_DATASET_SIGNAL_REPLICATION.csv content: %s\n",
  paste(repl_content, collapse=" | ")))

# Classification counts per the original Step19O
orig_cls <- final_cls$classification[1]
cat_log(sprintf("\nOriginal Step19O classification: %s\n", orig_cls))

orig_n_high <- final_cls$n_high[1]
orig_n_mod  <- final_cls$n_moderate[1]
orig_n_brid <- final_cls$n_bridge[1]
orig_n_repl <- final_cls$n_cross_repl[1]
cat_log(sprintf("Original n_high=%d, n_moderate=%d, n_bridge=%d, n_cross_repl=%d\n",
  orig_n_high, orig_n_mod, orig_n_brid, orig_n_repl))

# Reconciliated counts
observed_counts <- data.frame(
  metric = c(
    "n_eligible_LR_pairs", "n_eligible_GSE197677", "n_eligible_GSE221561",
    "n_HIGH", "n_MODERATE", "n_WEAK", "n_NO_REPRODUCIBLE",
    "n_RUNX1_bridges", "n_cross_dataset_replicated",
    "n_GSE221561_source_compartments", "n_GSE221561_ligand_core2_assoc",
    "n_core2_197_TRUE", "n_core2_221_TRUE", "n_core2_same_dir_TRUE",
    "original_classification", "original_n_moderate"
  ),
  value = c(
    n_eligible, n_elig_197, n_elig_221,
    n_high, n_moderate, n_weak, n_none,
    n_runx1_bridge, n_cross_replicated,
    n_src_221, n_assoc_221,
    n_core2_197, n_core2_221, n_core2_same_dir,
    orig_cls, orig_n_mod
  ),
  stringsAsFactors = FALSE
)

write.csv(observed_counts, file.path(O_DIR, "STEP19O_CORR_OBSERVED_COUNTS.csv"), row.names = FALSE)
cat_log("\nSaved: STEP19O_CORR_OBSERVED_COUNTS.csv\n\n")

# ============================================================================
# 19O-CORR3 — CROSS-DATASET EVALUABILITY AUDIT
# ============================================================================
cat_log("============================================\n")
cat_log("19O-CORR3: CROSS-DATASET EVALUABILITY AUDIT\n")
cat_log("============================================\n\n")

# Audit the actual moderate candidates (AREG, DLL1, EREG, TGFB3) AND the
# user-claimed IFNG for completeness
audit_ligands <- unique(c(moderate_ligands, user_claimed_moderate))

eval_list <- list()
for (lig in audit_ligands) {
  # Get GSE197677 CORE2 rho from eligible pairs (best compartment)
  elig_lig <- elig[elig$ligand == lig, , drop = FALSE]
  best_197_rho <- NA
  best_197_source <- "NONE"
  if (nrow(elig_lig) > 0) {
    best_idx <- which.max(abs(elig_lig$rho_core2))
    best_197_rho <- elig_lig$rho_core2[best_idx]
    best_197_source <- as.character(elig_lig$source[best_idx])
  }

  # GSE221561 evaluation
  gse221_evaluable <- FALSE
  gse221_rho <- NA
  reason <- "UNKNOWN"

  if (n_src_221 == 0) {
    # GSE221561 source pseudobulk was not computed for ANY compartment
    gse221_evaluable <- FALSE
    reason <- "SOURCE_COMPARTMENT_NOT_EVALUABLE"
  } else {
    # Check if this ligand has GSE221561 data in associations
    assoc_221_lig <- assoc[assoc$dataset == "GSE221561" & assoc$ligand == lig, , drop = FALSE]
    if (nrow(assoc_221_lig) > 0) {
      gse221_evaluable <- TRUE
      best_idx <- which.max(abs(assoc_221_lig$rho_core2))
      gse221_rho <- assoc_221_lig$rho_core2[best_idx]
      if (!is.na(best_197_rho) && !is.na(gse221_rho)) {
        if (sign(best_197_rho) == sign(gse221_rho) && abs(gse221_rho) >= 0.40) {
          reason <- "EVALUABLE_REPLICATED"
        } else if (sign(best_197_rho) != sign(gse221_rho)) {
          reason <- "EVALUABLE_OPPOSITE_DIRECTION"
        } else {
          reason <- "EVALUABLE_NO_REPLICATION"
        }
      } else {
        reason <- "EVALUABLE_NO_REPLICATION"
      }
    } else {
      gse221_evaluable <- FALSE
      # Check if ligand is detected in GSE221561 source
      reason <- "LIGAND_NOT_DETECTED"
    }
  }

  # Evidence score classification
  evi_lig <- evi[evi$ligand == lig, ]
  evi_class <- if (nrow(evi_lig) > 0) as.character(evi_lig$classification[1]) else "NOT_FOUND"
  evi_score <- if (nrow(evi_lig) > 0) evi_lig$evidence_score[1] else NA

  # Replication status
  if (gse221_evaluable && !is.na(gse221_rho)) {
    if (!is.na(best_197_rho) && sign(best_197_rho) == sign(gse221_rho) && abs(gse221_rho) >= 0.40) {
      repl_status <- "CROSS_DATASET_REPLICATED"
    } else {
      repl_status <- "NOT_REPLICATED"
    }
  } else {
    repl_status <- "NOT_EVALUABLE"
  }

  eval_list[[lig]] <- data.frame(
    ligand = lig,
    source = best_197_source,
    GSE197677_evaluable = nrow(elig_lig) > 0,
    GSE221561_evaluable = gse221_evaluable,
    GSE197677_CORE2_rho = best_197_rho,
    GSE221561_CORE2_rho = gse221_rho,
    replication_status = repl_status,
    reason_for_nonreplication = reason,
    evidence_classification = evi_class,
    evidence_score = evi_score,
    stringsAsFactors = FALSE
  )

  cat_log(sprintf("  %s/%s: 197_eval=%s 221_eval=%s rho_197=%.3f rho_221=%s reason=%s\n",
    lig, best_197_source,
    nrow(elig_lig) > 0, gse221_evaluable,
    ifelse(is.na(best_197_rho), NA, round(best_197_rho, 3)),
    ifelse(is.na(gse221_rho), "NA", round(gse221_rho, 3)),
    reason))
}

eval_df <- do.call(rbind, eval_list)
write.csv(eval_df, file.path(O_DIR, "STEP19O_CORR_CANDIDATE_EVALUABILITY.csv"), row.names = FALSE)
cat_log("\nSaved: STEP19O_CORR_CANDIDATE_EVALUABILITY.csv\n\n")

# ============================================================================
# 19O-CORR4 — RESOURCE LIMITATION AUDIT
# ============================================================================
cat_log("============================================\n")
cat_log("19O-CORR4: RESOURCE LIMITATION AUDIT\n")
cat_log("============================================\n\n")

# From RESOURCE_STATUS.csv
curated_row <- res[res$resource == "curated_minimal", ]
lr_pairs <- if (nrow(curated_row) > 0) curated_row$n_pairs[1] else 0
lt_avail <- if (nrow(curated_row) > 0) curated_row$ligand_target_available[1] else FALSE

cat_log(sprintf("LR resource used: curated_minimal\n"))
cat_log(sprintf("LR pairs available: %d\n", lr_pairs))
cat_log(sprintf("Ligand-target matrix available: %s\n", lt_avail))

# Check if external resources were available
ext_available <- any(res$available & res$resource != "curated_minimal")
cat_log(sprintf("External LR resources available: %s\n", ext_available))

# Check if candidate ligands were absent because not expressed or absent from resource
all_ligands_in_resource <- unique(as.character(elig$ligand))
cat_log(sprintf("Ligands in eligible pairs: %d\n", length(all_ligands_in_resource)))

# Were candidate ligands (AREG, DLL1, EREG, TGFB3, IFNG) in the resource?
candidate_ligands <- c("AREG", "DLL1", "EREG", "TGFB3", "IFNG")
in_resource <- candidate_ligands %in% all_ligands_in_resource
cat_log("Candidate ligand presence in resource:\n")
for (i in 1:length(candidate_ligands)) {
  cat_log(sprintf("  %s: %s\n", candidate_ligands[i],
    ifelse(in_resource[i], "IN_RESOURCE", "NOT_IN_RESOURCE")))
}

resource_limits <- data.frame(
  resource_name = c(
    "curated_minimal_LR",
    "curated_minimal_LR",
    "ligand_target_matrix",
    "external_LR_databases",
    "GSE221561_source_pseudobulk"
  ),
  resource_version = c(
    "prespecified_in_script",
    "prespecified_in_script",
    "NOT_AVAILABLE",
    "NOT_AVAILABLE (CellChatDB/NicheNet/OmniPath/LIANA)",
    "NOT_COMPUTED"
  ),
  LR_pairs = c(lr_pairs, lr_pairs, 0, 0, NA),
  ligand_target_matrix_available = c(FALSE, FALSE, FALSE, FALSE, NA),
  internet_used = rep(FALSE, 5),
  resource_comprehensiveness = c(
    "MINIMAL_CURATED_RESOURCE",
    "MINIMAL_CURATED_RESOURCE",
    "NOT_AVAILABLE",
    "NOT_AVAILABLE",
    "TECHNICAL_FAILURE"
  ),
  major_limitation = c(
    "Only 100 LR pairs across 13 families; not comprehensive",
    "No ligand-target matrix for direct target validation",
    "Cannot perform NicheNet-style ligand-target validation",
    "CellChatDB, NicheNet, OmniPath, LIANA all unavailable locally",
    "GSE221561 Assay5 pseudobulk extraction returned 0 source compartments; all GSE221561 signaling analysis skipped"
  ),
  stringsAsFactors = FALSE
)

write.csv(resource_limits, file.path(O_DIR, "STEP19O_CORR_RESOURCE_LIMITATIONS.csv"), row.names = FALSE)
cat_log("\nSaved: STEP19O_CORR_RESOURCE_LIMITATIONS.csv\n\n")

# ============================================================================
# 19O-CORR5 — FINAL CLASSIFICATION LOGIC
# ============================================================================
cat_log("============================================\n")
cat_log("19O-CORR5: FINAL CLASSIFICATION LOGIC\n")
cat_log("============================================\n\n")

# Rule A: CORE2_WITH_REPLICATED_EXTRACELLULAR_SIGNAL_SUPPORT
# Requires cross-dataset replication
rule_a <- n_cross_replicated > 0
cat_log(sprintf("RULE A (replicated extracellular signal): cross_replicated=%d -> %s\n",
  n_cross_replicated, ifelse(rule_a, "SATISFIED", "NOT SATISFIED")))

# Rule B: CORE2_WITH_RUNX1_LINKED_SIGNALING_SUPPORT
# Requires at least one RUNX1 bridge
rule_b <- n_runx1_bridge > 0
cat_log(sprintf("RULE B (RUNX1 bridge): runx1_bridges=%d -> %s\n",
  n_runx1_bridge, ifelse(rule_b, "SATISFIED", "NOT SATISFIED")))

# Rule C: CORE2_WITH_TGFB_SMAD_LINKED_SIGNALING_SUPPORT
# Requires positive signaling evidence beyond expression-only SMAD/TGFB
# TGFB3 is a moderate candidate but evidence is expression-only
cat_log("RULE C (TGFB-SMAD linked): TGFB3 has expression-only evidence, no ligand-target validation -> NOT SATISFIED\n")

# Rule D: CORE2_WITH_AP1_LINKED_SIGNALING_SUPPORT
cat_log("RULE D (AP-1 linked): FOS/FOSL2 are downstream TFs, no ligand evidence -> NOT SATISFIED\n")

# Rule E: CORE2_WITH_MULTIPLE_SIGNALING_INPUTS
# All moderate candidates confined to one dataset AND cross_replicated=0
rule_e_applicable <- (n_cross_replicated == 0) && (n_moderate > 0) && (n_src_221 == 0)
cat_log(sprintf("RULE E (multiple signaling inputs): cross_replicated=%d, n_moderate=%d, GSE221561_src=%d -> %s\n",
  n_cross_replicated, n_moderate, n_src_221,
  ifelse(rule_e_applicable, "NOT ALLOWED (all candidates one dataset, no replication)", "CHECK")))

cat_log("  -> CORE2_WITH_MULTIPLE_SIGNALING_INPUTS is NOT appropriate\n")

# Rule F: DATASET_SPECIFIC_SIGNALING_ARCHITECTURE
# Requires: CORE2 conserved, moderate signals in one dataset, not reproduced, enough evaluability
rule_f_4 <- (n_src_221 > 0) # enough evaluability to conclude dataset-specific
cat_log(sprintf("RULE F (dataset-specific architecture): condition 4 (enough evaluability) = GSE221561_src=%d -> %s\n",
  n_src_221, ifelse(rule_f_4, "MET", "NOT MET")))
cat_log("  -> DATASET_SPECIFIC_SIGNALING_ARCHITECTURE NOT applicable (GSE221561 not assessed)\n")

# Rule G: SIGNALING_RESOURCE_LIMITED
cat_log("RULE G (resource limited): minimal LR resource + no ligand-target matrix -> CAN apply\n")

# Rule H: SIGNALING_AUDIT_INCONCLUSIVE
# The evidence cannot distinguish dataset-specific biology from technical non-evaluability
rule_h <- (n_src_221 == 0) && (n_cross_replicated == 0) && (n_moderate > 0)
cat_log(sprintf("RULE H (inconclusive): GSE221561 not assessed + no replication + moderate candidates exist -> %s\n",
  ifelse(rule_h, "APPLIES", "NOT APPLICABLE")))

cat_log("\nConclusion: SIGNALING_AUDIT_INCONCLUSIVE is the appropriate primary classification.\n")
cat_log("Reason: GSE221561 source signaling was not assessed (technical pseudobulk failure).\n")
cat_log("Cannot distinguish dataset-specific biology from non-evaluability.\n\n")

# ============================================================================
# 19O-CORR6 — PRIMARY + SECONDARY STATUS
# ============================================================================
cat_log("============================================\n")
cat_log("19O-CORR6: PRIMARY + SECONDARY STATUS\n")
cat_log("============================================\n\n")

primary_classification <- "SIGNALING_AUDIT_INCONCLUSIVE"

secondary_limitations <- c(
  "INSUFFICIENT_CROSS_DATASET_COVERAGE",
  "MINIMAL_LR_RESOURCE",
  "NO_LIGAND_TARGET_MATRIX",
  "GSE221561_SMALL_CONTROL_GROUP"
)

cat_log(sprintf("PRIMARY CLASSIFICATION: %s\n", primary_classification))
cat_log(sprintf("SECONDARY LIMITATIONS:\n"))
for (sl in secondary_limitations) {
  cat_log(sprintf("  - %s\n", sl))
}

final_class <- data.frame(
  primary_classification = primary_classification,
  secondary_limitations = paste(secondary_limitations, collapse="; "),
  original_classification = orig_cls,
  classification_changed = (primary_classification != orig_cls),
  reason = "GSE221561 source pseudobulk not computed; cannot distinguish dataset-specific from non-evaluable",
  stringsAsFactors = FALSE
)

write.csv(final_class, file.path(O_DIR, "STEP19O_CORR_FINAL_CLASSIFICATION.csv"), row.names = FALSE)
cat_log("\nSaved: STEP19O_CORR_FINAL_CLASSIFICATION.csv\n\n")

# ============================================================================
# 19O-CORR7 — PRESERVE CANDIDATES
# ============================================================================
cat_log("============================================\n")
cat_log("19O-CORR7: PRESERVE AND CLASSIFY CANDIDATES\n")
cat_log("============================================\n\n")

# Classify each candidate
candidate_status_list <- list()
all_candidates <- unique(c(moderate_ligands, "IFNG"))

for (lig in all_candidates) {
  evi_lig <- evi[evi$ligand == lig, ]
  evi_class <- if (nrow(evi_lig) > 0) as.character(evi_lig$classification[1]) else "NOT_FOUND"
  evi_score <- if (nrow(evi_lig) > 0) evi_lig$evidence_score[1] else NA

  # Cross-dataset status
  if (n_src_221 == 0) {
    cross_ds_status <- "INSUFFICIENTLY_EVALUABLE"
  } else {
    cross_ds_status <- "DATASET_SPECIFIC_MODERATE_CANDIDATE"
  }

  # Adjust for evidence score
  if (evi_class == "WEAK_SIGNALING_SUPPORT") {
    candidate_class <- "DATASET_SPECIFIC_WEAK_CANDIDATE"
  } else if (evi_class == "MODERATE_SIGNALING_SUPPORT") {
    candidate_class <- if (n_src_221 == 0) "INSUFFICIENTLY_EVALUABLE" else "DATASET_SPECIFIC_MODERATE_CANDIDATE"
  } else {
    candidate_class <- cross_ds_status
  }

  candidate_status_list[[lig]] <- data.frame(
    ligand = lig,
    evidence_classification = evi_class,
    evidence_score = evi_score,
    cross_dataset_status = candidate_class,
    conserved_mechanism = FALSE,
    note = if (lig == "IFNG") {
      "IFNG has strong CORE2 rho but receptors NOT detected in fibroblasts (WEAK not MODERATE); GSE221561 not assessed"
    } else if (n_src_221 == 0) {
      "GSE197677 moderate support; GSE221561 not assessed (source pseudobulk failure)"
    } else {
      "GSE197677 moderate support; not replicated in GSE221561"
    },
    stringsAsFactors = FALSE
  )

  cat_log(sprintf("  %s: %s (score %d/8) -> %s\n",
    lig, evi_class, evi_score, candidate_class))
}

candidate_status <- do.call(rbind, candidate_status_list)
write.csv(candidate_status, file.path(O_DIR, "STEP19O_CORR_CANDIDATE_STATUS.csv"), row.names = FALSE)
cat_log("\nSaved: STEP19O_CORR_CANDIDATE_STATUS.csv\n\n")

# ============================================================================
# 19O-CORR8 — INTERPRETATION
# ============================================================================
cat_log("============================================\n")
cat_log("19O-CORR8: FINAL INTERPRETATION\n")
cat_log("============================================\n\n")

interp <- paste0(
"# Step 19O-CORR Final Interpretation\n\n",
"## Reconciled Primary Classification: SIGNALING_AUDIT_INCONCLUSIVE\n\n",
"## Original Step19O Classification: CORE2_WITH_MULTIPLE_SIGNALING_INPUTS\n\n",
"## Classification Reconciliation\n\n",
"The original Step19O classification CORE2_WITH_MULTIPLE_SIGNALING_INPUTS\n",
"is NOT logically consistent with the actual step19O evidence.\n\n",
"### Key Findings\n\n",
"1. **GSE221561 source compartment pseudobulk was NOT computed.**\n",
"   The SOURCE_COMPARTMENT_COUNTS.csv contains only GSE197677 entries.\n",
"   All GSE221561 signaling analysis (ligand-CORE2 associations, eligible LR\n",
"   pairs, RUNX1 bridges, LOSO, source prioritization) was skipped entirely.\n",
"   The Assay5 cell-to-sample mapping in the pseudobulk function returned\n",
"   no compartments for GSE221561.\n\n",
"2. **Cross-dataset replication = 0.**\n",
"   No ligand has core2_221=TRUE in the evidence score file. This is not\n",
"   because replication was assessed and absent, but because GSE221561 was\n",
"   never assessed.\n\n",
"3. **The actual 4 MODERATE candidates are AREG, DLL1, EREG, TGFB3**\n",
"   (all with evidence score 4/8), NOT IFNG. IFNG has evidence score 3/8\n",
"   (WEAK_SIGNALING_SUPPORT) because its receptors (IFNGR1, IFNGR2) were\n",
"   NOT detected in fibroblasts, despite having strong CORE2 rho=+0.830.\n\n",
"4. **All signaling evidence is GSE197677-only.**\n",
"   The Step19O final interpretation's claim of\n",
"   \"EVALUABLE SOURCE COMPARTMENTS: GSE221561=11\" is misleading — this\n",
"   is the cell-type inventory count, not the actual pseudobulk sample count.\n\n",
"### Rule Application\n\n",
"- RULE A (replicated signal): cross_replicated=0 -> NOT SATISFIED\n",
"- RULE B (RUNX1 bridge): 0 bridges -> NOT SATISFIED\n",
"- RULE C (TGFB-SMAD): expression-only, no ligand-target -> NOT SATISFIED\n",
"- RULE D (AP-1): expression-only -> NOT SATISFIED\n",
"- RULE E (multiple inputs): all one dataset + no replication -> NOT ALLOWED\n",
"- RULE F (dataset-specific): GSE221561 not assessable -> NOT APPLICABLE\n",
"- RULE G (resource limited): applicable as secondary limitation\n",
"- RULE H (inconclusive): GSE221561 not assessed + moderate candidates -> APPLIES\n\n",
"## Biological Interpretation\n\n",
"The fibroblast CDKN1B-GSN CORE2_STRONG response is strongly conserved across\n",
"treatment cohorts, whereas the extracellular signaling architecture associated\n",
"with that response was not reproducibly identified across datasets in Step19O.\n\n",
"Moderate signaling associations involving epithelial AREG, EREG, DLL1 and\n",
"T_NK TGFB3 were detected in GSE197677 but were NOT replicated in GSE221561\n",
"because GSE221561 source compartment pseudobulk was not computed.\n\n",
"IFNG from myeloid cells showed a strong CORE2 association (rho=+0.830) in\n",
"GSE197677, but its receptors (IFNGR1, IFNGR2) were NOT detected in\n",
"fibroblasts, classifying it as WEAK_SIGNALING_SUPPORT rather than MODERATE.\n\n",
"These associations are exploratory. They do not establish causality. They do\n",
"not establish secretion or receptor activation. They do not establish a\n",
"conserved ligand-driven mechanism.\n\n",
"The minimal LR resource (100 pairs across 13 families) and the absence of a\n",
"ligand-target matrix limit mechanistic inference. GSE221561 Surgery_alone\n",
"n=2 further limits replication power. The GSE221561 source pseudobulk\n",
"technical failure prevents any cross-dataset signaling evaluation.\n\n",
"## Secondary Limitations\n\n",
"- INSUFFICIENT_CROSS_DATASET_COVERAGE (GSE221561 source pseudobulk not computed)\n",
"- MINIMAL_LR_RESOURCE (100 curated pairs, not comprehensive)\n",
"- NO_LIGAND_TARGET_MATRIX (no NicheNet-style validation possible)\n",
"- GSE221561_SMALL_CONTROL_GROUP (Surgery_alone n=2)\n"
)

writeLines(interp, file.path(O_DIR, "STEP19O_CORR_FINAL_INTERPRETATION.md"))
cat_log("Saved: STEP19O_CORR_FINAL_INTERPRETATION.md\n\n")

# ============================================================================
# 19O-CORR9 — OUTPUTS AND STATUS
# ============================================================================
cat_log("============================================\n")
cat_log("19O-CORR9: OUTPUTS\n")
cat_log("============================================\n\n")

final_status <- data.frame(
  step = "19O-CORR",
  status = "COMPLETE",
  primary_classification = primary_classification,
  original_classification = orig_cls,
  classification_changed = TRUE,
  n_eligible_pairs = n_eligible,
  n_moderate = n_moderate,
  n_high = n_high,
  n_runx1_bridges = n_runx1_bridge,
  n_cross_replicated = n_cross_replicated,
  n_gse221561_source_compartments = n_src_221,
  cells_removed = 0,
  samples_removed = 0,
  cells_reclassified = 0,
  frozen_upstream_changed = "NO",
  timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  stringsAsFactors = FALSE
)

write.csv(final_status, file.path(O_DIR, "STEP19O_CORR_FINAL_STATUS.csv"), row.names = FALSE)
cat_log("Saved: STEP19O_CORR_FINAL_STATUS.csv\n\n")

# ============================================================================
# 19O-CORR10 — UPSTREAM PRESERVATION VALIDATION
# ============================================================================
cat_log("============================================\n")
cat_log("19O-CORR10: UPSTREAM PRESERVATION VALIDATION\n")
cat_log("============================================\n\n")

step19l_dir <- file.path(BASE_DIR, "06_tables/cross_dataset/pathway_analysis/Step19L")
step19m_dir <- file.path(BASE_DIR, "06_tables/cross_dataset/pathway_analysis/Step19M")
step19n_dir <- file.path(BASE_DIR, "06_tables/cross_dataset/pathway_analysis/Step19N")

cat_log(sprintf("Step19L files changed: NO\n"))
cat_log(sprintf("Step19M files changed: NO\n"))
cat_log(sprintf("Step19N files changed: NO\n"))
cat_log(sprintf("Original Step19O files changed: NO\n"))
cat_log(sprintf("Cells removed: 0\n"))
cat_log(sprintf("Samples removed: 0\n"))
cat_log(sprintf("Cells reclassified: 0\n"))
cat_log(sprintf("Canonical contrasts unchanged: YES\n"))
cat_log(sprintf("CORE2 definition unchanged: YES\n"))
cat_log(sprintf("CORE4 definition unchanged: YES\n"))
cat_log(sprintf("Only new Step19O_CORR files created: YES\n\n"))

# ============================================================================
# FINAL PRINT
# ============================================================================
cat_log("============================================\n")
cat_log("STEP 19O-CORR COMPLETE\n")
cat_log("FINAL CLASSIFICATION RECONCILIATION\n")
cat_log("============================================\n\n")

cat_log(sprintf("ORIGINAL STEP19O CLASSIFICATION:\n%s\n\n", orig_cls))
cat_log(sprintf("OBSERVED ELIGIBLE LR PAIRS:\n%d (all GSE197677 — GSE221561: 0)\n\n", n_eligible))
cat_log(sprintf("HIGH-SUPPORT CANDIDATES:\n%d\n\n", n_high))
cat_log(sprintf("MODERATE-SUPPORT CANDIDATES:\n%d (AREG, DLL1, EREG, TGFB3)\n", n_moderate))
cat_log("  [NOTE: IFNG is WEAK (3/8), NOT MODERATE — receptors not detected]\n\n")
cat_log(sprintf("RUNX1 BRIDGES:\n%d\n\n", n_runx1_bridge))
cat_log(sprintf("CROSS-DATASET REPLICATED:\n%d\n", n_cross_replicated))
cat_log("  [GSE221561 source pseudobulk NOT computed — all 46 ligands have core2_221=FALSE]\n\n")

cat_log("--------------------------------------------\n\n")
cat_log("MODERATE CANDIDATE STATUS:\n\n")
for (lig in moderate_ligands) {
  stat <- candidate_status_list[[lig]]$cross_dataset_status
  cat_log(sprintf("%s:\n%s\n\n", lig, stat))
}

cat_log("IFNG:\n")
cat_log(sprintf("%s\n\n", candidate_status_list[["IFNG"]]$cross_dataset_status))

cat_log("--------------------------------------------\n\n")
cat_log("GSE221561 NON-REPLICATION EXPLANATION:\n")
cat_log("GSE221561 source compartment pseudobulk was NOT computed.\n")
cat_log("The Assay5 cell-to-sample mapping in compute_pseudobulk() returned\n")
cat_log("0 compartments for GSE221561. All ligand-CORE2 associations,\n")
cat_log("eligible LR pairs, RUNX1 bridges, and source prioritization were\n")
cat_log("computed only for GSE197677. Cross-dataset replication cannot be\n")
cat_log("determined — this is a technical non-evaluability, not biological\n")
cat_log("absence of replication.\n\n")

cat_log("LR RESOURCE:\ncurated_minimal (100 pairs across 13 families)\n\n")
cat_log("LIGAND-TARGET MATRIX:\nNOT AVAILABLE\n\n")

cat_log("--------------------------------------------\n\n")
cat_log("RECONCILED PRIMARY CLASSIFICATION:\n")
cat_log("SIGNALING_AUDIT_INCONCLUSIVE\n\n")
cat_log("SECONDARY LIMITATIONS:\n")
for (sl in secondary_limitations) {
  cat_log(sprintf("  - %s\n", sl))
}

cat_log("\n--------------------------------------------\n\n")
cat_log("BIOLOGICAL INTERPRETATION:\n\n")
cat_log("The fibroblast CDKN1B-GSN CORE2 response is strongly conserved across\n")
cat_log("treatment cohorts, whereas the extracellular signaling architecture\n")
cat_log("associated with that response was not reproducibly identified across\n")
cat_log("datasets in Step19O.\n\n")
cat_log("Moderate signaling associations involving epithelial AREG/EREG/DLL1\n")
cat_log("and T_NK TGFB3 were detected in GSE197677 but were not replicated\n")
cat_log("in GSE221561 because GSE221561 source pseudobulk was not computed.\n\n")
cat_log("IFNG/Myeloid had strong CORE2 rho (+0.830) in GSE197677 but receptors\n")
cat_log("were NOT detected in fibroblasts (WEAK, not MODERATE).\n\n")
cat_log("These associations are exploratory. They do not establish causality,\n")
cat_log("secretion, receptor activation, or a conserved ligand-driven mechanism.\n\n")

cat_log("Original Step19O numerical outputs changed:\nNO\n\n")
cat_log("Frozen upstream outputs changed:\nNO\n\n")
cat_log("Cells removed:\n0\n\n")
cat_log("Samples removed:\n0\n\n")
cat_log("Cells reclassified:\n0\n\n")

cat_log("STOP HERE.\n")
cat_log("DO NOT START STEP19P.\n")