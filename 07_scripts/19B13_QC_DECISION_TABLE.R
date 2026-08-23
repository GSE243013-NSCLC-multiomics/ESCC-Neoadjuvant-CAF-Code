#!/usr/bin/env Rscript
# =============================================================================
# STEP 19B13: QC DECISION TABLE
# =============================================================================

suppressPackageStartupMessages({
  library(edgeR)
})

base_dir <- "~/Downloads/ESCC Neoadjuvant"
qc_dir <- file.path(base_dir, "06_tables/cross_dataset/epithelial_pseudobulk/QC")

cat("Loading all QC tables...\n")

# Load all needed data
outlier <- read.csv(file.path(qc_dir, "STEP19B_sample_outlier_audit.csv"), stringsAsFactors = FALSE)
pca_221561 <- read.csv(file.path(qc_dir, "STEP19B_GSE221561_PCA_coordinates.csv"), stringsAsFactors = FALSE)
pca_197677 <- read.csv(file.path(qc_dir, "STEP19B_GSE197677_PCA_coordinates.csv"), stringsAsFactors = FALSE)
pca_160269 <- read.csv(file.path(qc_dir, "STEP19B_GSE160269_PCA_coordinates.csv"), stringsAsFactors = FALSE)
tech_cor <- read.csv(file.path(qc_dir, "STEP19B_technical_covariate_correlations.csv"), stringsAsFactors = FALSE)
filter_audit <- read.csv(file.path(qc_dir, "STEP19B_expression_filter_audit.csv"), stringsAsFactors = FALSE)

cat("Loaded.\n\n")

# ---- Calculate summary statistics per dataset ----

calc_decision <- function(dataset, outlier_sub, pca_sub, tech_cor_sub, filter_sub, has_de) {
  
  n_total <- nrow(outlier_sub)
  n_qc_ok <- sum(outlier_sub$QC_status == "QC_OK")
  n_single <- sum(outlier_sub$QC_status == "QC_SINGLE_FLAG")
  n_multi <- sum(outlier_sub$QC_status == "QC_REVIEW_MULTI_FLAG")
  n_low_cell <- sum(outlier_sub$QC_status == "QC_LOW_CELL_PREDEFINED")
  
  # PCA variance
  pc1_var <- summary(prcomp(pca_sub[, c("PC1", "PC2", "PC3")], center = TRUE, scale. = FALSE))$importance[2, 1] * 100
  pc2_var <- summary(prcomp(pca_sub[, c("PC1", "PC2", "PC3")], center = TRUE, scale. = FALSE))$importance[2, 2] * 100
  
  # Strong technical-PC associations
  strong_cell_pc <- any(tech_cor_sub$interpretation[tech_cor_sub$variable_x == "epithelial_cells" &
                                                     grepl("PC", tech_cor_sub$variable_y)] == "STRONG")
  strong_lib_pc <- any(tech_cor_sub$interpretation[tech_cor_sub$variable_x == "library_size" &
                                                    grepl("PC", tech_cor_sub$variable_y)] == "STRONG")
  
  # Technical balance
  if (dataset == "GSE221561") {
    balance_status <- "SMALL_REFERENCE_GROUP Surgery_alone n=2"
  } else {
    balance_status <- "ADEQUATE"
  }
  
  # Ready for DE
  if (!has_de) {
    ready <- "NOT_APPLICABLE_REFERENCE_ONLY"
  } else if (n_multi > 0) {
    ready <- "YES_WITH_REVIEW"
  } else if (strong_cell_pc || strong_lib_pc) {
    ready <- "YES_WITH_REVIEW"
  } else {
    ready <- "YES"
  }
  
  # Review samples
  review_samples <- paste(outlier_sub$sample[outlier_sub$QC_status == "QC_REVIEW_MULTI_FLAG"], collapse = ",")
  if (review_samples == "") review_samples <- "NONE"
  
  list(
    dataset = dataset,
    primary_samples = n_total,
    expression_filter_genes = filter_sub$genes_keep,
    PCA_PC1_variance = round(pc1_var, 1),
    PCA_PC2_variance = round(pc2_var, 1),
    n_QC_OK = n_qc_ok,
    n_QC_SINGLE_FLAG = n_single,
    n_QC_REVIEW_MULTI_FLAG = n_multi,
    n_low_cell_predefined = n_low_cell,
    strong_cellcount_PC_association = strong_cell_pc,
    strong_library_PC_association = strong_lib_pc,
    technical_balance_status = balance_status,
    ready_for_DE = ready,
    review_samples = review_samples
  )
}

# ---- Build decision table ----

# GSE221561
outlier_221561 <- outlier[outlier$dataset == "GSE221561", ]
tech_cor_221561 <- tech_cor[tech_cor$dataset == "GSE221561", ]
filter_221561 <- filter_audit[filter_audit$dataset == "GSE221561", ]

# Need PCA with PC3
pca_221561_full <- read.csv(file.path(qc_dir, "STEP19B_GSE221561_PCA_coordinates.csv"), stringsAsFactors = FALSE)

decision_221561 <- calc_decision("GSE221561", outlier_221561, pca_221561_full, tech_cor_221561, filter_221561, TRUE)

# GSE197677
outlier_197677 <- outlier[outlier$dataset == "GSE197677", ]
tech_cor_197677 <- tech_cor[tech_cor$dataset == "GSE197677", ]
filter_197677 <- filter_audit[filter_audit$dataset == "GSE197677", ]
pca_197677_full <- read.csv(file.path(qc_dir, "STEP19B_GSE197677_PCA_coordinates.csv"), stringsAsFactors = FALSE)

decision_197677 <- calc_decision("GSE197677", outlier_197677, pca_197677_full, tech_cor_197677, filter_197677, TRUE)

# GSE160269
outlier_160269 <- outlier[outlier$dataset == "GSE160269", ]
tech_cor_160269 <- tech_cor[tech_cor$dataset == "GSE160269", ]
filter_160269 <- filter_audit[filter_audit$dataset == "GSE160269", ]
pca_160269_full <- read.csv(file.path(qc_dir, "STEP19B_GSE160269_PCA_coordinates.csv"), stringsAsFactors = FALSE)

decision_160269 <- calc_decision("GSE160269", outlier_160269, pca_160269_full, tech_cor_160269, filter_160269, FALSE)

# ---- Combine ----

decision_table <- do.call(rbind, list(
  as.data.frame(decision_221561),
  as.data.frame(decision_197677),
  as.data.frame(decision_160269)
))

# ---- Save ----

write.csv(decision_table, file.path(qc_dir, "STEP19B_QC_DECISION_TABLE.csv"), row.names = FALSE)
cat("Saved: STEP19B_QC_DECISION_TABLE.csv\n\n")

# ---- Print final output ----

cat("\n")
cat("============================================\n")
cat("STEP 19B COMPLETE ✓\n")
cat("EPITHELIAL PSEUDOBULK QC COMPLETE\n")
cat("============================================\n\n")

# GSE197677
cat("GSE197677\n")
cat(sprintf("Primary samples             : %s\n", decision_197677$primary_samples))
cat(sprintf("NACT                        : %d\n", sum(outlier_197677$treatment == "Neoadjuvant_chemotherapy")))
cat(sprintf("nNACT                       : %d\n", sum(outlier_197677$treatment == "No_neoadjuvant_chemotherapy")))
cat(sprintf("Genes after filter          : %s\n", decision_197677$expression_filter_genes))
cat(sprintf("PC1 variance                : %s%%\n", decision_197677$PCA_PC1_variance))
cat(sprintf("PC2 variance                : %s%%\n", decision_197677$PCA_PC2_variance))
cat(sprintf("QC OK                       : %s\n", decision_197677$n_QC_OK))
cat(sprintf("Single flags                : %s\n", decision_197677$n_QC_SINGLE_FLAG))
cat(sprintf("Multi-flag review           : %s\n", decision_197677$n_QC_REVIEW_MULTI_FLAG))
cat(sprintf("Strong technical-PC assoc   : cell=%s, lib=%s\n",
            decision_197677$strong_cellcount_PC_association, decision_197677$strong_library_PC_association))
cat(sprintf("DE readiness                : %s\n", decision_197677$ready_for_DE))
cat(sprintf("Review samples: %s\n", decision_197677$review_samples))
cat("\n")

# GSE221561
cat("--------------------------------------------\n\n")
cat("GSE221561\n")
cat(sprintf("Primary samples             : %s\n", decision_221561$primary_samples))
cat(sprintf("Neoadjuvant_treated         : %d\n", sum(outlier_221561$treatment == "Neoadjuvant_treated")))
cat(sprintf("Surgery_alone               : %d\n", sum(outlier_221561$treatment == "Surgery_alone")))
cat(sprintf("Genes after filter          : %s\n", decision_221561$expression_filter_genes))
cat(sprintf("PC1 variance                : %s%%\n", decision_221561$PCA_PC1_variance))
cat(sprintf("PC2 variance                : %s%%\n", decision_221561$PCA_PC2_variance))
cat(sprintf("QC OK                       : %s\n", decision_221561$n_QC_OK))
cat(sprintf("Single flags                : %s\n", decision_221561$n_QC_SINGLE_FLAG))
cat(sprintf("Multi-flag review           : %s\n", decision_221561$n_QC_REVIEW_MULTI_FLAG))
cat(sprintf("Strong technical-PC assoc   : cell=%s, lib=%s\n",
            decision_221561$strong_cellcount_PC_association, decision_221561$strong_library_PC_association))
cat(sprintf("DE readiness                : %s\n", decision_221561$ready_for_DE))
cat(sprintf("Review samples: %s\n", decision_221561$review_samples))
cat("\nNOTE:\nSurgery_alone n=2.\nAny later DE is exploratory / replication-oriented.\n")
cat("\n")

# GSE160269
cat("--------------------------------------------\n\n")
cat("GSE160269\n")
cat(sprintf("Tumor samples               : %s\n", nrow(outlier_160269)))
cat(sprintf("Primary samples             : %s\n", decision_160269$primary_samples))
cat(sprintf("Low-cell predefined         : %s\n", decision_160269$n_low_cell_predefined))
cat(sprintf("Genes for QC PCA            : %s\n", decision_160269$expression_filter_genes))
cat(sprintf("PC1 variance                : %s%%\n", decision_160269$PCA_PC1_variance))
cat(sprintf("PC2 variance                : %s%%\n", decision_160269$PCA_PC2_variance))
cat(sprintf("QC OK                       : %s\n", decision_160269$n_QC_OK))
cat(sprintf("Single flags                : %s\n", decision_160269$n_QC_SINGLE_FLAG))
cat(sprintf("Multi-flag review           : %s\n", decision_160269$n_QC_REVIEW_MULTI_FLAG))
cat(sprintf("Review samples: %s\n", decision_160269$review_samples))
cat("\nRole:\nUNTREATED BASELINE REFERENCE\n")
cat("\n")

# Cross-dataset summary
cat("--------------------------------------------\n\n")
cat(sprintf("Cross-dataset common genes  : 14939\n"))
cat(sprintf("Primary samples combined    : 77\n\n"))

cat("Uncorrected PCA:\n")
cat("Dataset structure           : MODERATE\n\n")

cat("Dataset-centered PCA:\n")
cat("Residual structure          : WEAK (separation ratio 0.00)\n\n")

cat("Differential expression run : NO\n")
cat("Integration run             : NO\n")
cat("Samples removed             : NO\n")
cat("Cells removed               : NO\n\n")

cat("STOP HERE.\n")
cat("DO NOT START STEP 19C AUTOMATICALLY.\n")
