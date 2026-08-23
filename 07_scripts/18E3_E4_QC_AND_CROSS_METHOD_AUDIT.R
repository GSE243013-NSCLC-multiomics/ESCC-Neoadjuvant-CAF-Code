#!/usr/bin/env Rscript
# ==============================================================================
# 18E3 + 18E4: PER-CONFIGURATION QC AND CROSS-METHOD AUDIT
# ==============================================================================

library(data.table)
library(ggplot2)
library(patchwork)

cat("============================================================\n")
cat("STEP 18E3 + 18E4: QC AND CROSS-METHOD AUDIT\n")
cat("============================================================\n\n")

OUTDIR <- "06_tables/cross_dataset/malignant_epithelial/GSE160269_copykat_validation/reference_audit"
FIGDIR <- "05_figures/GSE160269/malignant_epithelial/reference_audit"
dir.create(FIGDIR, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 1. LOAD ALL COPYKAT RESULTS
# ============================================================

cat("1. Loading CopyKAT results...\n\n")

configs <- list(
  list(tumor = "P16T", strategy = "CURRENT_18D", 
       pred_file = "06_tables/cross_dataset/malignant_epithelial/GSE160269_copykat_validation/P16T/P16T_copykat_predictions.csv",
       obj_file = "06_tables/cross_dataset/malignant_epithelial/GSE160269_copykat_validation/P16T/P16T_copykat_object.rds"),
  list(tumor = "P16T", strategy = "SAME_SAMPLE_NON_EPITHELIAL",
       pred_file = "06_tables/cross_dataset/malignant_epithelial/GSE160269_copykat_validation/reference_audit/P16T_SAME_SAMPLE_NON_EPITHELIAL/P16T_SAME_SAMPLE_NON_EPITHELIAL_predictions.csv",
       obj_file = "06_tables/cross_dataset/malignant_epithelial/GSE160269_copykat_validation/reference_audit/P16T_SAME_SAMPLE_NON_EPITHELIAL/P16T_SAME_SAMPLE_NON_EPITHELIAL_copykat_object.rds"),
  list(tumor = "P16T", strategy = "ADJACENT_NORMAL_EPITHELIAL",
       pred_file = "06_tables/cross_dataset/malignant_epithelial/GSE160269_copykat_validation/reference_audit/P16T_ADJACENT_NORMAL_EPITHELIAL/P16T_ADJACENT_NORMAL_EPITHELIAL_predictions.csv",
       obj_file = "06_tables/cross_dataset/malignant_epithelial/GSE160269_copykat_validation/reference_audit/P16T_ADJACENT_NORMAL_EPITHELIAL/P16T_ADJACENT_NORMAL_EPITHELIAL_copykat_object.rds"),
  list(tumor = "P39T", strategy = "CURRENT_18D",
       pred_file = "06_tables/cross_dataset/malignant_epithelial/GSE160269_copykat_validation/P39T/P39T_copykat_predictions.csv",
       obj_file = "06_tables/cross_dataset/malignant_epithelial/GSE160269_copykat_validation/P39T/P39T_copykat_object.rds"),
  list(tumor = "P39T", strategy = "SAME_SAMPLE_NON_EPITHELIAL",
       pred_file = "06_tables/cross_dataset/malignant_epithelial/GSE160269_copykat_validation/reference_audit/P39T_SAME_SAMPLE_NON_EPITHELIAL/P39T_SAME_SAMPLE_NON_EPITHELIAL_predictions.csv",
       obj_file = "06_tables/cross_dataset/malignant_epithelial/GSE160269_copykat_validation/reference_audit/P39T_SAME_SAMPLE_NON_EPITHELIAL/P39T_SAME_SAMPLE_NON_EPITHELIAL_copykat_object.rds"),
  list(tumor = "P39T", strategy = "ADJACENT_NORMAL_EPITHELIAL",
       pred_file = "06_tables/cross_dataset/malignant_epithelial/GSE160269_copykat_validation/reference_audit/P39T_ADJACENT_NORMAL_EPITHELIAL/P39T_ADJACENT_NORMAL_EPITHELIAL_predictions.csv",
       obj_file = "06_tables/cross_dataset/malignant_epithelial/GSE160269_copykat_validation/reference_audit/P39T_ADJACENT_NORMAL_EPITHELIAL/P39T_ADJACENT_NORMAL_EPITHELIAL_copykat_object.rds")
)

# Load annotations
p16t_annot <- fread("06_tables/cross_dataset/malignant_epithelial/GSE160269_copykat_validation/P16T_cell_annotations.csv")
p39t_annot <- fread("06_tables/cross_dataset/malignant_epithelial/GSE160269_copykat_validation/P39T_cell_annotations.csv")

# Load InferCNA results
infercna_metrics <- fread("06_tables/cross_dataset/malignant_epithelial/GSE160269_reference_sensitivity/GSE160269_ALL_SENSITIVITY_METRICS.csv")
infercna_full <- infercna_metrics[config == "BALANCED_WITHIN_DATASET"]

# ============================================================
# 2. COMPUTE METRICS FOR EACH CONFIGURATION
# ============================================================

cat("2. Computing per-configuration metrics...\n\n")

results <- list()

for (cfg in configs) {
  cat(sprintf("--- %s x %s ---\n", cfg$tumor, cfg$strategy))
  
  # Load predictions
  pred <- fread(cfg$pred_file)
  
  # Load annotations
  annot <- if (cfg$tumor == "P16T") p16t_annot else p39t_annot
  
  # Normalize barcodes for matching
  normalize_barcodes <- function(x) gsub("-", ".", x)
  pred[, cell_barcode_norm := normalize_barcodes(cell_barcode)]
  annot[, cell_barcode_norm := normalize_barcodes(cell_barcode)]
  
  # Separate tumor and normal
  tumor_cells <- annot[role == "TUMOR", cell_barcode_norm]
  normal_cells <- annot[role == "NORMAL", cell_barcode_norm]
  
  tumor_preds <- pred[cell_barcode_norm %in% tumor_cells]
  normal_preds <- pred[cell_barcode_norm %in% normal_cells]
  
  # Compute fractions
  n_tumor <- nrow(tumor_preds)
  n_normal <- nrow(normal_preds)
  
  tumor_aneuploid <- sum(tumor_preds$copykat.pred == "aneuploid", na.rm = TRUE)
  tumor_diploid <- sum(tumor_preds$copykat.pred == "diploid", na.rm = TRUE)
  tumor_not_defined <- sum(tumor_preds$copykat.pred == "not.defined", na.rm = TRUE)
  
  normal_aneuploid <- sum(normal_preds$copykat.pred == "aneuploid", na.rm = TRUE)
  normal_diploid <- sum(normal_preds$copykat.pred == "diploid", na.rm = TRUE)
  normal_not_defined <- sum(normal_preds$copykat.pred == "not.defined", na.rm = TRUE)
  
  tumor_aneuploid_frac <- tumor_aneuploid / n_tumor
  tumor_diploid_frac <- tumor_diploid / n_tumor
  tumor_not_defined_frac <- tumor_not_defined / n_tumor
  
  normal_aneuploid_frac <- if (n_normal > 0) normal_aneuploid / n_normal else NA
  normal_not_defined_frac <- if (n_normal > 0) normal_not_defined / n_normal else NA
  
  cat(sprintf("   Tumor: %d cells\n", n_tumor))
  cat(sprintf("     Aneuploid: %d (%.1f%%)\n", tumor_aneuploid, tumor_aneuploid_frac * 100))
  cat(sprintf("     Diploid: %d (%.1f%%)\n", tumor_diploid, tumor_diploid_frac * 100))
  cat(sprintf("     Not.defined: %d (%.1f%%)\n", tumor_not_defined, tumor_not_defined_frac * 100))
  cat(sprintf("   Normal: %d cells\n", n_normal))
  cat(sprintf("     Aneuploid: %d (%.1f%%)\n", normal_aneuploid, normal_aneuploid_frac * 100))
  
  # Load CNV matrix for profile analysis
  obj_exists <- file.exists(cfg$obj_file)
  profile_separation <- NA
  tumor_profile_cor <- NA
  normal_profile_cor <- NA
  chr_pattern <- "UNRESOLVED"
  
  if (obj_exists) {
    cat("   Loading CopyKAT object...\n")
    copykat_obj <- readRDS(cfg$obj_file)
    # CNA matrix is in CNAmat (data.frame: genes x cells)
    cnv_mat <- as.matrix(copykat_obj$CNAmat)
    rm(copykat_obj); gc(verbose = FALSE)
    
    # Compute profiles
    tumor_in_cnv <- intersect(tumor_cells, colnames(cnv_mat))
    normal_in_cnv <- intersect(normal_cells, colnames(cnv_mat))
    
    if (length(tumor_in_cnv) > 10) {
      tumor_cnv <- cnv_mat[, tumor_in_cnv, drop = FALSE]
      tumor_mean <- rowMeans(tumor_cnv, na.rm = TRUE)
      
      # Tumor cell correlation to tumor mean
      tumor_cors <- sapply(tumor_in_cnv, function(cell) {
        cor(cnv_mat[, cell], tumor_mean, use = "pairwise.complete.obs")
      })
      tumor_profile_cor <- median(tumor_cors, na.rm = TRUE)
      
      if (length(normal_in_cnv) > 10) {
        normal_cnv <- cnv_mat[, normal_in_cnv, drop = FALSE]
        normal_mean <- rowMeans(normal_cnv, na.rm = TRUE)
        profile_separation <- mean(abs(tumor_mean - normal_mean), na.rm = TRUE)
        normal_cors <- sapply(normal_in_cnv, function(cell) {
          cor(cnv_mat[, cell], tumor_mean, use = "pairwise.complete.obs")
        })
        normal_profile_cor <- median(normal_cors, na.rm = TRUE)
      } else {
        # No normal cells in CNV matrix (reference-only config)
        profile_separation <- NA
        normal_profile_cor <- NA
      }
      
      # Chromosome-scale pattern
      chr_var <- var(tumor_mean, na.rm = TRUE)
      if (chr_var > 0.01 && !is.na(profile_separation) && profile_separation > 0.1) {
        chr_pattern <- "BROAD_COHERENT_CNA"
      } else if (chr_var > 0.005) {
        chr_pattern <- "WEAK_DIFFUSE_SHIFT"
      } else {
        chr_pattern <- "MINIMAL_CNA"
      }
      
    cat(sprintf("   Profile separation: %s\n", ifelse(is.na(profile_separation), "NA (no normal cells)", sprintf("%.4f", profile_separation))))
    cat(sprintf("   Tumor profile cor: %.4f\n", tumor_profile_cor))
    cat(sprintf("   Normal profile cor: %s\n", ifelse(is.na(normal_profile_cor), "NA (no normal cells)", sprintf("%.4f", normal_profile_cor))))
    cat(sprintf("   Chromosome pattern: %s\n", chr_pattern))
    } else {
      cat("   Too few tumor cells in CNV matrix for profile analysis\n")
    }
  } else {
    cat("   CopyKAT object not found, skipping CNV profile analysis\n")
  }
  
  results[[paste0(cfg$tumor, "_", cfg$strategy)]] <- data.table(
    sample = cfg$tumor,
    reference_strategy = cfg$strategy,
    tumor_cells = n_tumor,
    normal_cells = n_normal,
    tumor_aneuploid = tumor_aneuploid,
    tumor_diploid = tumor_diploid,
    tumor_not_defined = tumor_not_defined,
    tumor_aneuploid_candidate_fraction = tumor_aneuploid_frac,
    tumor_diploid_fraction = tumor_diploid_frac,
    tumor_not_defined_fraction = tumor_not_defined_frac,
    known_normal_aneuploid = normal_aneuploid,
    known_normal_aneuploid_fraction = normal_aneuploid_frac,
    known_normal_not_defined_fraction = normal_not_defined_frac,
    profile_separation = profile_separation,
    tumor_profile_correlation = tumor_profile_cor,
    normal_profile_correlation = normal_profile_cor,
    chromosome_scale_pattern = chr_pattern
  )
  
  cat("\n")
}

# ============================================================
# 3. REFERENCE ROBUSTNESS CLASSIFICATION
# ============================================================

cat("3. Classifying reference robustness...\n\n")

classify_robustness <- function(tumor_name, dt) {
  tumor_dt <- dt[sample == tumor_name]
  
  aneuploid_fracs <- tumor_dt$tumor_aneuploid_candidate_fraction
  normal_aneuploid_fracs <- tumor_dt$known_normal_aneuploid_fraction
  patterns <- tumor_dt$chromosome_scale_pattern
  
  # Check if aneuploid fraction is consistently high
  high_aneuploid <- all(aneuploid_fracs > 0.80, na.rm = TRUE)
  
  # Check if normal aneuploid is low
  low_normal <- all(normal_aneuploid_fracs < 0.10, na.rm = TRUE)
  
  # Check if chromosome pattern persists
  coherent_pattern <- sum(patterns == "BROAD_COHERENT_CNA") >= 2
  
  # Check if profile separation persists
  separations <- tumor_dt$profile_separation
  good_separation <- all(separations > 0.05, na.rm = TRUE)
  
  # Classification
  if (high_aneuploid && low_normal && (coherent_pattern || good_separation)) {
    return("COPYKAT_REFERENCE_ROBUST_STRONG")
  } else if (length(unique(round(aneuploid_fracs, 2))) > 2 || 
             diff(range(aneuploid_fracs, na.rm = TRUE)) > 0.30) {
    return("COPYKAT_REFERENCE_SENSITIVE")
  } else if (all(aneuploid_fracs < 0.50, na.rm = TRUE)) {
    return("COPYKAT_REFERENCE_ROBUST_WEAK")
  } else {
    return("COPYKAT_INCONCLUSIVE")
  }
}

final_dt <- rbindlist(results)

for (tumor in c("P16T", "P39T")) {
  robustness <- classify_robustness(tumor, final_dt)
  final_dt[sample == tumor, reference_robustness_status := robustness]
  cat(sprintf("%s: %s\n", tumor, robustness))
}

# ============================================================
# 4. CROSS-METHOD AUDIT
# ============================================================

cat("\n4. Cross-method audit...\n\n")

# Create cross-method table
cross_method <- list()

for (tumor in c("P16T", "P39T")) {
  # InferCNA (frozen)
  inf <- infercna_full[tumor_id == tumor]
  
  # CopyKAT (average across strategies)
  ck_dt <- final_dt[sample == tumor]
  avg_aneuploid <- mean(ck_dt$tumor_aneuploid_candidate_fraction, na.rm = TRUE)
  avg_normal_aneuploid <- mean(ck_dt$known_normal_aneuploid_fraction, na.rm = TRUE)
  robustness <- unique(ck_dt$reference_robustness_status)
  
  # Interpretation
  if (inf$cna_signal_ratio < 3.0 && inf$coherent_arm_count == 0 && avg_aneuploid > 0.80) {
    interpretation <- "METHOD_DISCORDANCE: InferCNA weak but CopyKAT strong"
  } else if (inf$cna_signal_ratio < 3.0 && avg_aneuploid < 0.50) {
    interpretation <- "CONCORDANT_WEAK: Both methods show weak CNA"
  } else {
    interpretation <- "MIXED_EVIDENCE"
  }
  
  cross_method[[tumor]] <- data.table(
    sample = tumor,
    method = "InferCNA",
    reference_strategy = "FULL_REFERENCE",
    tumor_cells = inf$tumor_cells,
    normal_cells = 183,
    aneuploid_candidate_fraction = NA,
    known_normal_aneuploid_fraction = NA,
    not_defined_fraction = NA,
    profile_separation = NA,
    tumor_profile_correlation = NA,
    chromosome_scale_pattern = "0_coherent_arms",
    reference_robustness_status = paste0("cnaCor=", round(inf$tumor_cnaCor, 3), "; ratio=", round(inf$cna_signal_ratio, 1))
  )
  
  cross_method[[paste0(tumor, "_CopyKAT")]] <- data.table(
    sample = tumor,
    method = "CopyKAT",
    reference_strategy = "MULTIPLE",
    tumor_cells = mean(ck_dt$tumor_cells),
    normal_cells = mean(ck_dt$normal_cells),
    aneuploid_candidate_fraction = avg_aneuploid,
    known_normal_aneuploid_fraction = avg_normal_aneuploid,
    not_defined_fraction = mean(ck_dt$tumor_not_defined_fraction, na.rm = TRUE),
    profile_separation = mean(ck_dt$profile_separation, na.rm = TRUE),
    tumor_profile_correlation = mean(ck_dt$tumor_profile_correlation, na.rm = TRUE),
    chromosome_scale_pattern = paste(unique(ck_dt$chromosome_scale_pattern), collapse = "/"),
    reference_robustness_status = robustness
  )
}

cross_dt <- rbindlist(cross_method, fill = TRUE)
fwrite(cross_dt, file.path(OUTDIR, "STEP18E_P16T_P39T_CROSS_METHOD_REFERENCE_AUDIT.csv"))
cat("Saved: STEP18E_P16T_P39T_CROSS_METHOD_REFERENCE_AUDIT.csv\n")

# ============================================================
# 5. GENERATE FIGURES
# ============================================================

cat("\n5. Generating figures...\n\n")

# Figure 1: Aneuploid fraction by reference
p1 <- ggplot(final_dt, aes(x = reference_strategy, y = tumor_aneuploid_candidate_fraction, fill = reference_strategy)) +
  geom_col() +
  facet_wrap(~sample, scales = "free_x") +
  coord_flip() +
  labs(title = "CopyKAT Aneuploid Candidate Fraction by Reference Strategy",
       x = "Reference Strategy", y = "Aneuploid Candidate Fraction") +
  theme_minimal() +
  theme(legend.position = "none")

ggsave(file.path(FIGDIR, "STEP18E_CopyKAT_aneuploid_fraction_by_reference.png"), p1, width = 10, height = 6, dpi = 300)

# Figure 2: Known-normal false positive
p2 <- ggplot(final_dt, aes(x = reference_strategy, y = known_normal_aneuploid_fraction, fill = reference_strategy)) +
  geom_col() +
  facet_wrap(~sample, scales = "free_x") +
  coord_flip() +
  labs(title = "Known-Normal Aneuploid False Positive Rate",
       x = "Reference Strategy", y = "Known-Normal Aneuploid Fraction") +
  theme_minimal() +
  theme(legend.position = "none")

ggsave(file.path(FIGDIR, "STEP18E_known_normal_false_positive_by_reference.png"), p2, width = 10, height = 6, dpi = 300)

cat("Figures saved to:", FIGDIR, "\n")

# ============================================================
# 6. FINAL OUTPUT
# ============================================================

cat("\n============================================\n")
cat("STEP 18E COMPLETE\n")
cat("GSE160269 COPYKAT REFERENCE AUDIT\n")
cat("============================================\n\n")

for (tumor in c("P16T", "P39T")) {
  cat(sprintf("%s\n", tumor))
  
  for (strategy in c("CURRENT_18D", "SAME_SAMPLE_NON_EPITHELIAL", "ADJACENT_NORMAL_EPITHELIAL")) {
    row <- final_dt[sample == tumor & reference_strategy == strategy]
    if (nrow(row) > 0) {
      cat(sprintf("%s:\n", strategy))
      cat(sprintf("Tumor aneuploid candidate %% : %.1f%%\n", row$tumor_aneuploid_candidate_fraction * 100))
      cat(sprintf("Known-normal aneuploid %%    : %.1f%%\n", row$known_normal_aneuploid_fraction * 100))
    }
  }
  
  robustness <- unique(final_dt[sample == tumor, reference_robustness_status])
  cat(sprintf("Chromosome-scale CNA pattern: %s\n", paste(unique(final_dt[sample == tumor, chromosome_scale_pattern]), collapse = "/")))
  cat(sprintf("Reference robustness status : %s\n", robustness))
  
  # Cross-method
  inf_row <- infercna_full[tumor_id == tumor]
  inf_ratio <- inf_row$cna_signal_ratio
  inf_cnaCor <- inf_row$tumor_cnaCor
  inf_arms <- inf_row$coherent_arm_count
  ck_status <- robustness
  
  # Interpretation based on InferCNA metrics + CopyKAT results
  if (inf_ratio < 3.0 && inf_arms == 0 && grepl("SENSITIVE", ck_status)) {
    interp <- "METHOD_DISCORDANCE: InferCNA weak (ratio<3, 0 arms) but CopyKAT reference-sensitive"
  } else if (inf_ratio < 3.0 && inf_arms == 0) {
    interp <- "CONCORDANT_WEAK: InferCNA weak, CopyKAT variable"
  } else {
    interp <- "MIXED_EVIDENCE"
  }
  
  cat(sprintf("\nCross-method status:\n"))
  cat(sprintf("InferCNA ratio              : %.1f\n", inf_ratio))
  cat(sprintf("InferCNA cnaCor             : %.3f\n", inf_cnaCor))
  cat(sprintf("InferCNA coherent arms      : %d\n", inf_arms))
  cat(sprintf("CopyKAT robustness          : %s\n", ck_status))
  cat(sprintf("Interpretation              : %s\n", interp))
  
  cat("--------------------------------------------\n\n")
}

cat("CopyKAT reruns              : P16T/P39T only\n")
cat("InferCNA rerun              : NO\n")
cat("findMalignant               : NO\n")
cat("Malignant labels assigned   : NO\n")
cat("Other 58 tumors processed   : NO\n")
cat("Cells removed               : NO\n\n")
cat("STOP HERE.\n")
