#!/usr/bin/env Rscript
# ==============================================================================
# 18E2: COPYKAT REFERENCE SENSITIVITY - RUN B AND C
# ==============================================================================

library(copykat)
library(data.table)

cat("============================================================\n")
cat("STEP 18E2: COPYKAT REFERENCE SENSITIVITY\n")
cat("============================================================\n\n")

OUTDIR <- "06_tables/cross_dataset/malignant_epithelial/GSE160269_copykat_validation/reference_audit"
BASEDIR <- "06_tables/cross_dataset/malignant_epithelial/GSE160269_copykat_validation"

# ============================================================
# 1. LOAD RAW COUNTS AND ANNOTATIONS
# ============================================================

cat("1. Loading raw counts...\n\n")

# Load full raw counts from 18D preflight
raw_counts_full <- readRDS(file.path(BASEDIR, "P16T_raw_counts.rds"))
# Note: P16T and P39T have the same normal cells, so we can reuse the matrix structure

# Load tumor annotations
p16t_annot <- fread(file.path(BASEDIR, "P16T_cell_annotations.csv"))
p39t_annot <- fread(file.path(BASEDIR, "P39T_cell_annotations.csv"))

cat("P16T:", sum(p16t_annot$role == "TUMOR"), "tumor cells\n")
cat("P39T:", sum(p39t_annot$role == "TUMOR"), "tumor cells\n\n")

# ============================================================
# 2. LOAD REFERENCE STRATEGIES
# ============================================================

cat("2. Loading reference strategies...\n\n")

ref_same_sample_p16t <- fread(file.path(OUTDIR, "reference_P16T_SAME_SAMPLE_NON_EPITHELIAL.csv"))
ref_same_sample_p39t <- fread(file.path(OUTDIR, "reference_P39T_SAME_SAMPLE_NON_EPITHELIAL.csv"))
ref_adjacent <- fread(file.path(OUTDIR, "reference_ADJACENT_NORMAL_EPITHELIAL.csv"))

cat("P16T same-sample:", nrow(ref_same_sample_p16t), "cells\n")
cat("P39T same-sample:", nrow(ref_same_sample_p39t), "cells\n")
cat("Adjacent-normal epithelial:", nrow(ref_adjacent), "cells\n\n")

# ============================================================
# 3. HELPER FUNCTION TO RUN COPYKAT
# ============================================================

run_copykat_sensitivity <- function(tumor_name, strategy_name, norm_cell_barcodes) {
  cat(sprintf("\n--- %s x %s ---\n", tumor_name, strategy_name))
  
  # Create output directory
  out_subdir <- file.path(OUTDIR, paste0(tumor_name, "_", strategy_name))
  dir.create(out_subdir, recursive = TRUE, showWarnings = FALSE)
  
  # Load tumor raw counts
  tumor_annot <- if (tumor_name == "P16T") p16t_annot else p39t_annot
  tumor_cells <- tumor_annot[role == "TUMOR", cell_barcode]
  
  # Get all cells we need
  all_cells <- unique(c(tumor_cells, norm_cell_barcodes))
  
  # Load raw counts and subset
  raw_counts <- readRDS(file.path(BASEDIR, paste0(tumor_name, "_raw_counts.rds")))
  cells_in_data <- intersect(all_cells, colnames(raw_counts))
  
  cat(sprintf("   Tumor cells: %d\n", length(tumor_cells)))
  cat(sprintf("   Normal cells needed: %d\n", length(norm_cell_barcodes)))
  cat(sprintf("   Cells in data: %d\n", length(cells_in_data)))
  
  # Subset to needed cells
  raw_sub <- raw_counts[, cells_in_data, drop = FALSE]
  
  # Convert to data.frame
  raw_df <- as.data.frame(as.matrix(raw_sub))
  
  # Identify which normal cells are in the data
  norm_in_data <- intersect(norm_cell_barcodes, cells_in_data)
  cat(sprintf("   Normal cells found: %d\n", length(norm_in_data)))
  
  if (length(norm_in_data) < 50) {
    cat("   WARNING: Too few normal cells for reliable baseline\n")
  }
  
  # Run CopyKAT
  cat("   Running CopyKAT...\n")
  start_time <- Sys.time()
  
  result <- tryCatch({
    copykat(
      rawmat = raw_df,
      id.type = "S",
      cell.line = "no",
      ngene.chr = 5,
      LOW.DR = 0.05,
      UP.DR = 0.1,
      win.size = 25,
      norm.cell.names = norm_in_data,
      sam.name = paste0("GSE160269_", tumor_name, "_", strategy_name),
      n.cores = 1,
      output.seg = FALSE,
      plot.genes = FALSE,
      genome = "hg20"
    )
  }, error = function(e) {
    cat(sprintf("   ERROR: %s\n", conditionMessage(e)))
    return(NULL)
  })
  
  end_time <- Sys.time()
  runtime <- as.numeric(difftime(end_time, start_time, units = "mins"))
  cat(sprintf("   Runtime: %.2f min\n", runtime))
  
  if (is.null(result)) {
    cat("   FAILED\n")
    return(NULL)
  }
  
  # Save results
  saveRDS(result, file.path(out_subdir, paste0(tumor_name, "_", strategy_name, "_copykat_object.rds")))
  fwrite(as.data.table(result$prediction, keep.rownames = "cell_barcode"),
         file.path(out_subdir, paste0(tumor_name, "_", strategy_name, "_predictions.csv")))
  
  if (!is.null(result$CNV.matrix)) {
    saveRDS(result$CNV.matrix, file.path(out_subdir, paste0(tumor_name, "_", strategy_name, "_CNV_matrix.rds")))
  }
  
  cat(sprintf("   Results saved to: %s\n", out_subdir))
  
  return(result)
}

# ============================================================
# 4. RUN STRATEGY B: SAME_SAMPLE_NON_EPITHELIAL
# ============================================================

cat("\n============================================================\n")
cat("RUNNING STRATEGY B: SAME_SAMPLE_NON_EPITHELIAL\n")
cat("============================================================\n")

# P16T
result_p16t_b <- run_copykat_sensitivity("P16T", "SAME_SAMPLE_NON_EPITHELIAL", 
                                          ref_same_sample_p16t$cell_barcode)

# P39T
result_p39t_b <- run_copykat_sensitivity("P39T", "SAME_SAMPLE_NON_EPITHELIAL",
                                          ref_same_sample_p39t$cell_barcode)

# ============================================================
# 5. RUN STRATEGY C: ADJACENT_NORMAL_EPITHELIAL
# ============================================================

cat("\n============================================================\n")
cat("RUNNING STRATEGY C: ADJACENT_NORMAL_EPITHELIAL\n")
cat("============================================================\n")

# P16T
result_p16t_c <- run_copykat_sensitivity("P16T", "ADJACENT_NORMAL_EPITHELIAL",
                                          ref_adjacent$cell_barcode)

# P39T
result_p39t_c <- run_copykat_sensitivity("P39T", "ADJACENT_NORMAL_EPITHELIAL",
                                          ref_adjacent$cell_barcode)

# ============================================================
# 6. SUMMARY
# ============================================================

cat("\n============================================================\n")
cat("STEP 18E2 COPYKAT REFERENCE SENSITIVITY COMPLETE\n")
cat("============================================================\n\n")

cat("Strategy B (SAME_SAMPLE_NON_EPITHELIAL):\n")
cat(sprintf("  P16T: %s\n", ifelse(is.null(result_p16t_b), "FAILED", "SUCCESS")))
cat(sprintf("  P39T: %s\n", ifelse(is.null(result_p39t_b), "FAILED", "SUCCESS")))

cat("\nStrategy C (ADJACENT_NORMAL_EPITHELIAL):\n")
cat(sprintf("  P16T: %s\n", ifelse(is.null(result_p16t_c), "FAILED", "SUCCESS")))
cat(sprintf("  P39T: %s\n", ifelse(is.null(result_p39t_c), "FAILED", "SUCCESS")))

cat("\nFiles saved to:", OUTDIR, "\n")
