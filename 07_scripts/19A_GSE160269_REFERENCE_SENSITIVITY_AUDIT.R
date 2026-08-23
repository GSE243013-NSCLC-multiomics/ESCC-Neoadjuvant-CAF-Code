#!/usr/bin/env Rscript
# ==============================================================================
# 19A: GSE160269 NORMAL REFERENCE SENSITIVITY AUDIT
# Reruns P16T and P39T InferCNA with 5 reference configurations + balanced analysis
# ==============================================================================

library(Seurat)
library(SeuratObject)
library(infercna)
library(data.table)
library(ggplot2)
library(patchwork)

# ==============================================================================
# CONFIGURATION
# ==============================================================================

TUMORS <- c("P16T", "P39T")
CONFIGS <- c("FULL_REFERENCE", "LEAVE_P126N_OUT", "LEAVE_P127N_OUT",
             "LEAVE_P128N_OUT", "LEAVE_P130N_OUT")
ALL_NORMALS <- c("P126N", "P127N", "P128N", "P130N")

OUTDIR <- "06_tables/cross_dataset/malignant_epithelial/GSE160269_reference_sensitivity"
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# InferCNA parameters (identical to all previous runs)
CNA_PARAMS <- list(
  genome = "hg38",
  window = 100,
  range = c(-3, 3),
  n = 5000,
  noise = 0.1,
  center.method = "median",
  isLog = FALSE
)

# ==============================================================================
# FUNCTIONS
# ==============================================================================

compute_reference_config <- function(all_normals, leave_out) {
  setdiff(all_normals, leave_out)
}

create_balanced_refCells <- function(refCells_original, seed = 42) {
  sizes <- sapply(refCells_original, length)
  min_size <- min(sizes)
  set.seed(seed)
  balanced <- lapply(refCells_original, function(x) {
    if (length(x) > min_size) {
      sample(x, min_size)
    } else {
      x
    }
  })
  return(balanced)
}

compute_cna_metrics <- function(cna_mat, tumor_cells, ref_cells_list, gene_table) {
  # cna_mat is the output matrix from infercna (genes x cells)

  tumor_cna <- cna_mat[, tumor_cells, drop = FALSE]
  all_ref_cells <- unlist(ref_cells_list)
  ref_cna <- cna_mat[, all_ref_cells, drop = FALSE]

  # Use infercna's cnaSignal function
  signal <- infercna::cnaSignal(cna_mat)
  if (is.null(names(signal))) names(signal) <- colnames(cna_mat)

  tumor_signal <- as.numeric(signal[tumor_cells])
  ref_signal <- as.numeric(signal[all_ref_cells])

  tumor_median <- median(tumor_signal, na.rm = TRUE)
  ref_median <- median(ref_signal, na.rm = TRUE)
  ratio <- tumor_median / ref_median

  # cnaCor using infercna's function
  tumor_cor <- infercna::cnaCor(tumor_cna)
  if (is.null(names(tumor_cor))) names(tumor_cor) <- colnames(tumor_cna)
  tumor_cnaCor <- median(as.numeric(tumor_cor), na.rm = TRUE)

  # Normal-to-tumor correlation (QC)
  tumor_mean_profile <- rowMeans(tumor_cna)
  normal_cor_vals <- sapply(all_ref_cells, function(cell) {
    suppressWarnings(cor(cna_mat[, cell], tumor_mean_profile, method = "pearson", use = "pairwise.complete.obs"))
  })
  normal_cor_median <- median(normal_cor_vals, na.rm = TRUE)

  # Separation
  separation <- tumor_cnaCor - normal_cor_median

  # Percentiles
  q95 <- quantile(ref_signal, 0.95, na.rm = TRUE)
  q99 <- quantile(ref_signal, 0.99, na.rm = TRUE)
  pct_above_95 <- mean(tumor_signal > q95, na.rm = TRUE) * 100
  pct_above_99 <- mean(tumor_signal > q99, na.rm = TRUE) * 100

  # Coherent arms
  coherent_arms <- c()
  chroms <- unique(gene_table$chromosome_name)
  for (chrom in chroms) {
    genes_in_arm <- gene_table[gene_table$chromosome_name == chrom, ]
    arms_in_chrom <- unique(genes_in_arm$arm)
    for (arm in arms_in_chrom) {
      arm_genes <- genes_in_arm[genes_in_arm$arm == arm, "symbol"]
      arm_genes_in_mat <- intersect(arm_genes, rownames(cna_mat))
      if (length(arm_genes_in_mat) < 5) next

      arm_tumor <- tumor_cna[arm_genes_in_mat, , drop = FALSE]
      arm_ref <- ref_cna[arm_genes_in_mat, , drop = FALSE]

      tumor_mean <- rowMeans(arm_tumor, na.rm = TRUE)
      ref_mean <- rowMeans(arm_ref, na.rm = TRUE)
      delta <- tumor_mean - ref_mean

      gain_fraction <- mean(delta > 0.05, na.rm = TRUE)
      loss_fraction <- mean(delta < -0.05, na.rm = TRUE)
      abs_delta <- mean(abs(delta), na.rm = TRUE)

      if (abs_delta >= 0.05 && (gain_fraction >= 0.50 || loss_fraction >= 0.50)) {
        direction <- if (gain_fraction >= 0.50) "gain" else "loss"
        coherent_arms <- c(coherent_arms, paste0(chrom, arm, "_", direction))
      }
    }
  }

  list(
    tumor_cells = length(tumor_cells),
    total_normal = length(all_ref_cells),
    n_ref_groups = length(ref_cells_list),
    tumor_median_signal = tumor_median,
    ref_median_signal = ref_median,
    cna_signal_ratio = ratio,
    tumor_cnaCor = tumor_cnaCor,
    normal_cor_median = normal_cor_median,
    separation = separation,
    pct_above_95 = pct_above_95,
    pct_above_99 = pct_above_99,
    coherent_arm_count = length(coherent_arms),
    coherent_arm_names = paste(coherent_arms, collapse = "; ")
  )
}

# ==============================================================================
# MAIN
# ==============================================================================

results_list <- list()

for (tumor in TUMORS) {
  cat("\n============================================================\n")
  cat(sprintf("TUMOR: %s\n", tumor))
  cat("============================================================\n")

  # Load prepared data
  base_dir <- sprintf("06_tables/cross_dataset/malignant_epithelial/GSE160269_%s_pilot", tumor)
  cpm_file <- file.path(base_dir, sprintf("GSE160269_%s_infercna_CPM_hg38.rds", tumor))
  annot_file <- file.path(base_dir, sprintf("GSE160269_%s_infercna_cell_annotation.csv", tumor))
  gene_file <- file.path(base_dir, sprintf("GSE160269_%s_infercna_gene_table.csv", tumor))
  manifest_file <- file.path(base_dir, sprintf("GSE160269_%s_infercna_manifest.csv", tumor))

  cpm_combined <- readRDS(cpm_file)
  annot <- fread(annot_file)
  gene_table <- fread(gene_file)
  manifest <- fread(manifest_file)

  tumor_cells <- annot[infercna_role == "TUMOR_TEST", cell_barcode]
  all_ref_cells <- annot[infercna_role == "NORMAL_REFERENCE", cell_barcode]

  # Get normal cells by specimen
  refCells_original <- list()
  for (n in ALL_NORMALS) {
    refCells_original[[n]] <- annot[infercna_role == "NORMAL_REFERENCE" & sample_id == n, cell_barcode]
  }

  cat(sprintf("Tumor cells: %d\n", length(tumor_cells)))
  cat("Normal cells by group:\n")
  for (n in ALL_NORMALS) {
    cat(sprintf("  %s: %d\n", n, length(refCells_original[[n]])))
  }

  # Standard configurations (A-E)
  for (config in CONFIGS) {
    cat(sprintf("\n--- Config: %s ---\n", config))

    if (config == "FULL_REFERENCE") {
      ref_cells <- refCells_original
      label <- "FULL"
    } else {
      leave_out <- sub("LEAVE_(.*)_OUT", "\\1", config)
      ref_cells <- refCells_original[setdiff(ALL_NORMALS, leave_out)]
      label <- paste0("NO_", leave_out)
    }

    cat(sprintf("  Reference groups: %d\n", length(ref_cells)))
    cat(sprintf("  Total reference cells: %d\n", sum(sapply(ref_cells, length))))

    # Run InferCNA
    cna_result <- infercna::infercna(
      m = as.matrix(cpm_combined),
      refCells = ref_cells,
      window = CNA_PARAMS$window,
      range = CNA_PARAMS$range,
      n = CNA_PARAMS$n,
      noise = CNA_PARAMS$noise,
      center.method = CNA_PARAMS$center.method,
      isLog = CNA_PARAMS$isLog,
      verbose = TRUE
    )

    metrics <- compute_cna_metrics(cna_result, tumor_cells, ref_cells, gene_table)
    metrics$tumor_id <- tumor
    metrics$config <- config
    metrics$analysis_type <- "UNBALANCED"

    results_list[[length(results_list) + 1]] <- metrics

    cat(sprintf("  CNA signal ratio: %.3f\n", metrics$cna_signal_ratio))
    cat(sprintf("  cnaCor: %.3f\n", metrics$tumor_cnaCor))
    cat(sprintf("  Coherent arms: %d\n", metrics$coherent_arm_count))
  }

  # Balanced analysis
  cat("\n--- BALANCED_WITHIN_DATASET ---\n")
  balanced_refCells <- create_balanced_refCells(refCells_original, seed = 42)

  cat("  Balanced group sizes:\n")
  for (n in names(balanced_refCells)) {
    cat(sprintf("    %s: %d\n", n, length(balanced_refCells[[n]])))
  }
  cat(sprintf("  Total balanced reference cells: %d\n", sum(sapply(balanced_refCells, length))))

  cna_result_balanced <- infercna::infercna(
    m = as.matrix(cpm_combined),
    refCells = balanced_refCells,
    window = CNA_PARAMS$window,
    range = CNA_PARAMS$range,
    n = CNA_PARAMS$n,
    noise = CNA_PARAMS$noise,
    center.method = CNA_PARAMS$center.method,
    isLog = CNA_PARAMS$isLog,
    verbose = TRUE
  )

  metrics_balanced <- compute_cna_metrics(cna_result_balanced, tumor_cells, balanced_refCells, gene_table)
  metrics_balanced$tumor_id <- tumor
  metrics_balanced$config <- "BALANCED_WITHIN_DATASET"
  metrics_balanced$analysis_type <- "BALANCED"

  results_list[[length(results_list) + 1]] <- metrics_balanced

  cat(sprintf("  CNA signal ratio: %.3f\n", metrics_balanced$cna_signal_ratio))
  cat(sprintf("  cnaCor: %.3f\n", metrics_balanced$tumor_cnaCor))
  cat(sprintf("  Coherent arms: %d\n", metrics_balanced$coherent_arm_count))

  # Save CNA results for each config
  config_outdir <- file.path(OUTDIR, tumor)
  dir.create(config_outdir, recursive = TRUE, showWarnings = FALSE)

  # Save full CNA for balanced (most conservative)
  saveRDS(cna_result_balanced, file.path(config_outdir, "cna_balanced.rds"))

  cat(sprintf("\n  Saved balanced CNA results to: %s\n", config_outdir))
  rm(cna_result_balanced)
  gc(verbose = FALSE)
}

# ==============================================================================
# COMPILE RESULTS TABLE
# ==============================================================================

cat("\n\n============================================================\n")
cat("SENSITIVITY AUDIT RESULTS\n")
cat("============================================================\n")

results_dt <- rbindlist(results_list)
fwrite(results_dt, file.path(OUTDIR, "reference_sensitivity_metrics.csv"))

cat("\nFull metrics saved to: reference_sensitivity_metrics.csv\n")

# Print summary table
cat("\n--- P16T RESULTS ---\n")
p16t <- results_dt[tumor_id == "P16T"]
print(p16t[, .(config, tumor_cells, total_normal, n_ref_groups,
               tumor_median_signal, ref_median_signal, cna_signal_ratio,
               tumor_cnaCor, separation, pct_above_95, pct_above_99,
               coherent_arm_count)], nrows = 20)

cat("\n--- P39T RESULTS ---\n")
p39t <- results_dt[tumor_id == "P39T"]
print(p39t[, .(config, tumor_cells, total_normal, n_ref_groups,
               tumor_median_signal, ref_median_signal, cna_signal_ratio,
               tumor_cnaCor, separation, pct_above_95, pct_above_99,
               coherent_arm_count)], nrows = 20)

# ==============================================================================
# REFERENCE STABILITY CLASSIFICATION
# ==============================================================================

cat("\n\n============================================================\n")
cat("REFERENCE STABILITY CLASSIFICATION\n")
cat("============================================================\n")

classify_stability <- function(dt) {
  # Extract FULL and BALANCED results
  full <- dt[config == "FULL_REFERENCE"]
  balanced <- dt[config == "BALANCED_WITHIN_DATASET"]
  leave_out <- dt[!config %in% c("FULL_REFERENCE", "BALANCED_WITHIN_DATASET")]

  # Check consistency across leave-out configs
  ratios <- leave_out$cna_signal_ratio
  cnaCors <- leave_out$tumor_cnaCor
  arms <- leave_out$coherent_arm_count

  ratio_range <- max(ratios) - min(ratios)
  cnaCor_range <- max(cnaCors) - min(cnaCors)

  cat(sprintf("\n  Leave-out configs:\n"))
  cat(sprintf("    CNA signal ratio range: %.3f - %.3f (delta=%.3f)\n",
              min(ratios), max(ratios), ratio_range))
  cat(sprintf("    cnaCor range: %.3f - %.3f (delta=%.3f)\n",
              min(cnaCors), max(cnaCors), cnaCor_range))
  cat(sprintf("    Coherent arms range: %d - %d\n", min(arms), max(arms)))

  # Classify
  # REFERENCE_STABLE_STRONG: full ratio >= 3.0 AND cnaCor >= 0.6 AND all leave-out ratios >= 2.5
  # REFERENCE_STABLE_WEAK: full ratio < 3.0 OR cnaCor < 0.6, but leave-out configs consistent
  # REFERENCE_SENSITIVE: large variation across configs
  # INCONCLUSIVE: inconsistent

  if (full$cna_signal_ratio >= 3.0 && full$tumor_cnaCor >= 0.6 &&
      min(ratios) >= 2.5 && min(cnaCors) >= 0.5) {
    return("REFERENCE_STABLE_STRONG")
  } else if (ratio_range <= 1.0 && cnaCor_range <= 0.15 &&
             abs(full$cna_signal_ratio - balanced$cna_signal_ratio) <= 0.5) {
    return("REFERENCE_STABLE_WEAK")
  } else if (ratio_range > 2.0 || cnaCor_range > 0.2) {
    return("REFERENCE_SENSITIVE")
  } else {
    return("INCONCLUSIVE")
  }
}

for (tumor in TUMORS) {
  tumor_dt <- results_dt[tumor_id == tumor]
  stability <- classify_stability(tumor_dt)
  cat(sprintf("\n%s: %s\n", tumor, stability))
}

# ==============================================================================
# SENSITIVITY COMPARISON PLOTS
# ==============================================================================

cat("\n\nGenerating sensitivity comparison plots...\n")

for (tumor in TUMORS) {
  tumor_dt <- results_dt[tumor_id == tumor]

  # CNA signal ratio comparison
  p1 <- ggplot(tumor_dt, aes(x = config, y = cna_signal_ratio, fill = analysis_type)) +
    geom_col() +
    geom_hline(yintercept = 3.0, linetype = "dashed", color = "red") +
    coord_flip() +
    labs(title = sprintf("%s: CNA Signal Ratio by Reference Configuration", tumor),
         x = "Configuration", y = "CNA Signal Ratio (T/N)") +
    theme_minimal() +
    theme(legend.position = "bottom")

  # cnaCor comparison
  p2 <- ggplot(tumor_dt, aes(x = config, y = tumor_cnaCor, fill = analysis_type)) +
    geom_col() +
    geom_hline(yintercept = 0.6, linetype = "dashed", color = "red") +
    coord_flip() +
    labs(title = sprintf("%s: cnaCor by Reference Configuration", tumor),
         x = "Configuration", y = "Tumor Median cnaCor") +
    theme_minimal() +
    theme(legend.position = "bottom")

  # Coherent arms comparison
  p3 <- ggplot(tumor_dt, aes(x = config, y = coherent_arm_count, fill = analysis_type)) +
    geom_col() +
    coord_flip() +
    labs(title = sprintf("%s: Coherent CNA Arms by Reference Configuration", tumor),
         x = "Configuration", y = "Coherent Arm Count") +
    theme_minimal() +
    theme(legend.position = "bottom")

  # Combined plot
  combined <- (p1 / p2 / p3) + plot_layout(guides = "collect")

  ggsave(
    filename = file.path(OUTDIR, sprintf("%s_sensitivity_comparison.png", tumor)),
    plot = combined,
    width = 12, height = 12, dpi = 300
  )

  cat(sprintf("  Saved: %s_sensitivity_comparison.png\n", tumor))
}

cat("\n\n============================================================\n")
cat("REFERENCE SENSITIVITY AUDIT COMPLETE\n")
cat("============================================================\n")
cat(sprintf("Results saved to: %s\n", OUTDIR))
