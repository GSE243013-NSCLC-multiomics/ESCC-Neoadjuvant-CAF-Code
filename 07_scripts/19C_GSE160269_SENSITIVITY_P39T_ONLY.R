#!/usr/bin/env Rscript
# ==============================================================================
# 19C: GSE160269 REFERENCE SENSITIVITY - P39T ONLY
# ==============================================================================

library(Seurat)
library(SeuratObject)
library(infercna)
library(data.table)

# Configuration
ALL_NORMALS <- c("P126N", "P127N", "P128N", "P130N")
OUTDIR <- "06_tables/cross_dataset/malignant_epithelial/GSE160269_reference_sensitivity"

# Load prepared data
base_dir <- "06_tables/cross_dataset/malignant_epithelial/GSE160269_P39T_pilot"
cpm_combined <- readRDS(file.path(base_dir, "GSE160269_P39T_infercna_CPM_hg38.rds"))
annot <- fread(file.path(base_dir, "GSE160269_P39T_infercna_cell_annotation.csv"))
gene_table <- fread(file.path(base_dir, "GSE160269_P39T_infercna_gene_table.csv"))

tumor_cells <- annot[infercna_role == "TUMOR_TEST", cell_barcode]

# Get normal cells by specimen
refCells_original <- list()
for (n in ALL_NORMALS) {
  refCells_original[[n]] <- annot[infercna_role == "NORMAL_REFERENCE" & sample_id == n, cell_barcode]
}

cat("P39T tumor cells:", length(tumor_cells), "\n")
cat("Normal cells:\n")
for (n in ALL_NORMALS) cat(sprintf("  %s: %d\n", n, length(refCells_original[[n]])))

# Function to compute metrics
compute_metrics <- function(cna_mat, tumor_cells, ref_cells_list, gene_table) {
  tumor_cna <- cna_mat[, tumor_cells, drop = FALSE]
  all_ref_cells <- unlist(ref_cells_list)
  
  # Signal
  signal <- infercna::cnaSignal(cna_mat)
  if (is.null(names(signal))) names(signal) <- colnames(cna_mat)
  tumor_signal <- as.numeric(signal[tumor_cells])
  ref_signal <- as.numeric(signal[all_ref_cells])
  
  tumor_median <- median(tumor_signal, na.rm = TRUE)
  ref_median <- median(ref_signal, na.rm = TRUE)
  
  # cnaCor
  tumor_cor <- infercna::cnaCor(tumor_cna)
  tumor_cnaCor <- median(as.numeric(tumor_cor), na.rm = TRUE)
  
  # Normal-to-tumor correlation
  tumor_mean_profile <- rowMeans(tumor_cna)
  normal_cor_vals <- sapply(all_ref_cells, function(cell) {
    suppressWarnings(cor(cna_mat[, cell], tumor_mean_profile, method = "pearson", use = "pairwise.complete.obs"))
  })
  normal_cor_median <- median(normal_cor_vals, na.rm = TRUE)
  
  # Percentiles
  q95 <- quantile(ref_signal, 0.95, na.rm = TRUE)
  q99 <- quantile(ref_signal, 0.99, na.rm = TRUE)
  
  # Coherent arms
  coherent_arms <- c()
  chroms <- unique(gene_table$chromosome_name)
  for (chrom in chroms) {
    genes_in_arm <- gene_table[gene_table$chromosome_name == chrom, ]
    for (arm in unique(genes_in_arm$arm)) {
      arm_genes <- genes_in_arm[genes_in_arm$arm == arm, "symbol"]
      arm_genes_in_mat <- intersect(arm_genes, rownames(cna_mat))
      if (length(arm_genes_in_mat) < 5) next
      
      arm_tumor <- tumor_cna[arm_genes_in_mat, , drop = FALSE]
      arm_ref <- cna_mat[arm_genes_in_mat, all_ref_cells, drop = FALSE]
      
      delta <- rowMeans(arm_tumor, na.rm = TRUE) - rowMeans(arm_ref, na.rm = TRUE)
      gain_frac <- mean(delta > 0.05, na.rm = TRUE)
      loss_frac <- mean(delta < -0.05, na.rm = TRUE)
      abs_delta <- mean(abs(delta), na.rm = TRUE)
      
      if (abs_delta >= 0.05 && (gain_frac >= 0.50 || loss_frac >= 0.50)) {
        direction <- if (gain_frac >= 0.50) "gain" else "loss"
        coherent_arms <- c(coherent_arms, paste0(chrom, arm, "_", direction))
      }
    }
  }
  
  data.table(
    tumor_cells = length(tumor_cells),
    total_normal = length(all_ref_cells),
    n_ref_groups = length(ref_cells_list),
    tumor_median_signal = tumor_median,
    ref_median_signal = ref_median,
    cna_signal_ratio = tumor_median / ref_median,
    tumor_cnaCor = tumor_cnaCor,
    normal_cor_median = normal_cor_median,
    separation = tumor_cnaCor - normal_cor_median,
    pct_above_95 = mean(tumor_signal > q95, na.rm = TRUE) * 100,
    pct_above_99 = mean(tumor_signal > q99, na.rm = TRUE) * 100,
    coherent_arm_count = length(coherent_arms),
    coherent_arm_names = paste(coherent_arms, collapse = "; ")
  )
}

# Run all configurations
configs <- list(
  FULL_REFERENCE = ALL_NORMALS,
  LEAVE_P126N_OUT = c("P127N", "P128N", "P130N"),
  LEAVE_P127N_OUT = c("P126N", "P128N", "P130N"),
  LEAVE_P128N_OUT = c("P126N", "P127N", "P130N"),
  LEAVE_P130N_OUT = c("P126N", "P127N", "P128N")
)

results <- list()

for (config_name in names(configs)) {
  cat(sprintf("\n--- %s ---\n", config_name))
  ref_cells <- refCells_original[configs[[config_name]]]
  cat(sprintf("  Reference cells: %d\n", sum(sapply(ref_cells, length))))
  
  cna_result <- infercna::infercna(
    m = as.matrix(cpm_combined),
    refCells = ref_cells,
    window = 100, range = c(-3, 3), n = 5000,
    noise = 0.1, center.method = "median", isLog = FALSE,
    verbose = TRUE
  )
  
  metrics <- compute_metrics(cna_result, tumor_cells, ref_cells, gene_table)
  metrics$config <- config_name
  metrics$analysis_type <- "UNBALANCED"
  results[[config_name]] <- metrics
  
  cat(sprintf("  Ratio: %.3f, cnaCor: %.3f, Arms: %d\n",
              metrics$cna_signal_ratio, metrics$tumor_cnaCor, metrics$coherent_arm_count))
  
  rm(cna_result)
  gc(verbose = FALSE)
}

# Balanced analysis
cat("\n--- BALANCED_WITHIN_DATASET ---\n")
min_size <- min(sapply(refCells_original, length))
set.seed(42)
balanced <- lapply(refCells_original, function(x) sample(x, min_size))
cat(sprintf("  Balanced to %d cells per group\n", min_size))

cna_balanced <- infercna::infercna(
  m = as.matrix(cpm_combined),
  refCells = balanced,
  window = 100, range = c(-3, 3), n = 5000,
  noise = 0.1, center.method = "median", isLog = FALSE,
  verbose = TRUE
)

metrics_balanced <- compute_metrics(cna_balanced, tumor_cells, balanced, gene_table)
metrics_balanced$config <- "BALANCED_WITHIN_DATASET"
metrics_balanced$analysis_type <- "BALANCED"
results[["BALANCED"]] <- metrics_balanced

# Compile results
results_dt <- rbindlist(results)
results_dt$tumor_id <- "P39T"

fwrite(results_dt, file.path(OUTDIR, "P39T_sensitivity_metrics.csv"))

cat("\n\n============================================================\n")
cat("P39T SENSITIVITY RESULTS\n")
cat("============================================================\n\n")

print(results_dt[, .(config, tumor_cells, total_normal, n_ref_groups,
                      cna_signal_ratio, tumor_cnaCor, separation,
                      pct_above_95, pct_above_99, coherent_arm_count)])

cat("\nSaved: P39T_sensitivity_metrics.csv\n")
