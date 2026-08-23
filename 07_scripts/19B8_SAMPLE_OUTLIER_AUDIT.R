#!/usr/bin/env Rscript
# =============================================================================
# STEP 19B8: SAMPLE OUTLIER AUDIT
# =============================================================================

suppressPackageStartupMessages({
  library(edgeR)
})

base_dir <- "~/Downloads/ESCC Neoadjuvant"
qc_dir <- file.path(base_dir, "06_tables/cross_dataset/epithelial_pseudobulk/QC")

cat("Loading data...\n")

dge_221561 <- readRDS(file.path(qc_dir, "GSE221561_DGEList_TMM.rds"))
dge_197677 <- readRDS(file.path(qc_dir, "GSE197677_DGEList_TMM.rds"))
dge_160269 <- readRDS(file.path(qc_dir, "GSE160269_DGEList_TMM.rds"))

logcpm_221561 <- readRDS(file.path(qc_dir, "GSE221561_logCPM_QC.rds"))
logcpm_197677 <- readRDS(file.path(qc_dir, "GSE197677_logCPM_QC.rds"))
logcpm_160269 <- readRDS(file.path(qc_dir, "GSE160269_logCPM_QC.rds"))

cor_221561 <- read.csv(file.path(qc_dir, "STEP19B_GSE221561_sample_correlation.csv"), row.names = 1, check.names = FALSE)
cor_197677 <- read.csv(file.path(qc_dir, "STEP19B_GSE197677_sample_correlation.csv"), row.names = 1, check.names = FALSE)
cor_160269 <- read.csv(file.path(qc_dir, "STEP19B_GSE160269_sample_correlation.csv"), row.names = 1, check.names = FALSE)

pca_221561 <- read.csv(file.path(qc_dir, "STEP19B_GSE221561_PCA_coordinates.csv"), stringsAsFactors = FALSE)
pca_197677 <- read.csv(file.path(qc_dir, "STEP19B_GSE197677_PCA_coordinates.csv"), stringsAsFactors = FALSE)
pca_160269 <- read.csv(file.path(qc_dir, "STEP19B_GSE160269_PCA_coordinates.csv"), stringsAsFactors = FALSE)

cat("Loaded.\n\n")

# ---- Outlier detection function ----

detect_outliers <- function(dge, logcpm, cor_mat, pca_coords, dataset, predefined_low = character(0)) {
  cat(sprintf("=== %s OUTLIER AUDIT ===\n\n", dataset))
  
  samples <- dge$samples$sample
  n <- length(samples)
  
  # Build metrics matrix
  metrics <- data.frame(
    sample = samples,
    log10_epithelial_cells = log10(dge$samples$epithelial_cells),
    log10_library_size = log10(dge$samples$lib.size),
    genes_detected = colSums(logcpm > 0)[samples],
    PC1 = pca_coords$PC1[match(samples, pca_coords$sample)],
    PC2 = pca_coords$PC2[match(samples, pca_coords$sample)],
    stringsAsFactors = FALSE
  )
  
  # Mean sample correlation
  mean_cor <- numeric(n)
  for (i in 1:n) {
    sample_name <- samples[i]
    # Get correlations with all other samples
    cors <- as.numeric(cor_mat[sample_name, ])
    cors <- cors[-which(samples == sample_name)]  # exclude self
    mean_cor[i] <- mean(cors)
  }
  metrics$mean_sample_correlation <- mean_cor
  
  # Calculate robust z-scores
  robust_z <- data.frame(matrix(NA, nrow = n, ncol = 6,
                                 dimnames = list(samples,
                                                 c("log10_cells", "log10_libsize", "genes_detected",
                                                   "PC1", "PC2", "mean_cor"))))
  
  metric_cols <- c("log10_epithelial_cells", "log10_library_size", "genes_detected", "PC1", "PC2", "mean_sample_correlation")
  colnames(robust_z) <- c("log10_cells", "log10_libsize", "genes_detected", "PC1", "PC2", "mean_cor")
  
  for (j in 1:6) {
    x <- metrics[[metric_cols[j]]]
    med <- median(x)
    mad_val <- mad(x, constant = 1.4826)
    if (mad_val > 0) {
      robust_z[, j] <- (x - med) / mad_val
    }
  }
  
  # Flag outliers (abs(robust_z) > 3.5)
  # For mean_cor: only flag LOW correlation
  flags <- data.frame(matrix(FALSE, nrow = n, ncol = 6,
                              dimnames = list(samples, colnames(robust_z))))
  
  flags[, 1:5] <- abs(robust_z[, 1:5]) > 3.5
  flags[, "mean_cor"] <- robust_z[, "mean_cor"] < -3.5  # only low correlation
  
  # Count flags per sample
  n_flags <- rowSums(flags)
  flagged_metrics <- apply(flags, 1, function(x) paste(names(x)[x], collapse = ","))
  
  # Assign QC status
  qc_status <- character(n)
  for (i in 1:n) {
    if (samples[i] %in% predefined_low) {
      qc_status[i] <- "QC_LOW_CELL_PREDEFINED"
    } else if (n_flags[i] == 0) {
      qc_status[i] <- "QC_OK"
    } else if (n_flags[i] == 1) {
      qc_status[i] <- "QC_SINGLE_FLAG"
    } else {
      qc_status[i] <- "QC_REVIEW_MULTI_FLAG"
    }
  }
  
  # Build result table
  result <- data.frame(
    dataset = dataset,
    sample = samples,
    treatment = dge$samples$treatment_standardized,
    epithelial_cells = dge$samples$epithelial_cells,
    library_size = dge$samples$lib.size,
    genes_detected = metrics$genes_detected,
    PC1 = metrics$PC1,
    PC2 = metrics$PC2,
    mean_sample_correlation = round(metrics$mean_sample_correlation, 4),
    n_outlier_metrics = n_flags,
    flagged_metrics = flagged_metrics,
    QC_status = qc_status,
    stringsAsFactors = FALSE
  )
  
  # Print summary
  cat(sprintf("  QC_OK:                 %d\n", sum(qc_status == "QC_OK")))
  cat(sprintf("  QC_SINGLE_FLAG:        %d\n", sum(qc_status == "QC_SINGLE_FLAG")))
  cat(sprintf("  QC_REVIEW_MULTI_FLAG:  %d\n", sum(qc_status == "QC_REVIEW_MULTI_FLAG")))
  cat(sprintf("  QC_LOW_CELL_PREDEFINED:%d\n", sum(qc_status == "QC_LOW_CELL_PREDEFINED")))
  
  multi_flag <- samples[qc_status == "QC_REVIEW_MULTI_FLAG"]
  if (length(multi_flag) > 0) {
    cat("\n  Samples requiring review:\n")
    for (s in multi_flag) {
      cat(sprintf("    %s: %s\n", s, flagged_metrics[samples == s]))
    }
  }
  
  cat("\n")
  return(result)
}

# ---- Run outlier audit ----

outlier_221561 <- detect_outliers(dge_221561, logcpm_221561, cor_221561, pca_221561, "GSE221561")
outlier_197677 <- detect_outliers(dge_197677, logcpm_197677, cor_197677, pca_197677, "GSE197677")
outlier_160269 <- detect_outliers(dge_160269, logcpm_160269, cor_160269, pca_160269, "GSE160269",
                                   predefined_low = c("P63T", "P87T"))

# ---- Combine and save ----

outlier_all <- rbind(outlier_221561, outlier_197677, outlier_160269)

write.csv(outlier_all, file.path(qc_dir, "STEP19B_sample_outlier_audit.csv"), row.names = FALSE)
cat("Saved: STEP19B_sample_outlier_audit.csv\n")

cat("\nDone.\n")
