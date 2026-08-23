#!/usr/bin/env Rscript
# Step 13: Build Seurat objects for GSE197677 and GSE160269
# GSE197677: Load pre-processed Seurat RDS objects
# GSE160269: Read UMI matrices one at a time, create per-type Seurat, merge

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(readr)
})

PROJECT <- "~/Downloads/ESCC Neoadjuvant"
OBJ_DIR <- file.path(PROJECT, "03_objects")
RAW_DIR <- file.path(PROJECT, "01_raw/GEO")

dir.create(file.path(OBJ_DIR, "GSE197677"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OBJ_DIR, "GSE160269"), recursive = TRUE, showWarnings = FALSE)

cat("=== Step 13: Build Seurat Objects ===\n")
cat("Start:", format(Sys.time()), "\n\n")

# ============================================================
# PART A: GSE197677 — Load pre-processed Seurat objects
# ============================================================
cat("--- GSE197677: Loading pre-processed Seurat objects ---\n")

gse197677_dir <- file.path(RAW_DIR, "GSE197677")

cat("Loading integrated object...\n")
integrated <- readRDS(file.path(gse197677_dir, "GSE197677_ESCC.integrated.rds"))
cat("  Integrated class:", class(integrated)[1], "\n")

saveRDS(integrated, file.path(OBJ_DIR, "GSE197677", "GSE197677_integrated.rds"))
cat("  Saved: GSE197677_integrated.rds\n")
cat("  Note: CD45 object skipped (8.5 GB exceeds memory)\n")

rm(integrated)
gc()

# ============================================================
# PART B: GSE160269 — Build Seurat from UMI matrices
# ============================================================
cat("\n--- GSE160269: Building Seurat from UMI matrices ---\n")

gse160269_dir <- file.path(RAW_DIR, "GSE160269")

# Read cell metadata
cat("  Reading cell metadata...\n")
cells_neg <- read_tsv(file.path(gse160269_dir, "GSE160269_CD45neg_cells.txt.gz"),
                       col_types = cols(), progress = FALSE)
cells_pos <- read_tsv(file.path(gse160269_dir, "GSE160269_CD45pos_cells.txt.gz"),
                       col_types = cols(), progress = FALSE)

cells_all <- bind_rows(cells_neg, cells_pos)
cat("    Total cells:", nrow(cells_all), "\n")

# Process each cell-type matrix separately
umi_files <- list.files(gse160269_dir, pattern = "UMI_matrix_.*\\.txt\\.gz$", full.names = TRUE)
cat("  Found", length(umi_files), "cell-type UMI matrix files\n")

seu_list <- list()
for (f in umi_files) {
  fname <- basename(f)
  cell_type <- gsub("UMI_matrix_(.*)\\.txt\\.gz", "\\1", fname)
  cat("\n  Processing:", cell_type, "\n")
  
  # Read matrix (genes x cells, with header)
  mat <- read.table(f, header = TRUE, row.names = 1, 
                    check.names = FALSE, stringsAsFactors = FALSE)
  cat("    Raw dims:", nrow(mat), "genes x", ncol(mat), "cells\n")
  
  # Convert to sparse
  sparse_mat <- as(as.matrix(mat), "dgCMatrix")
  rm(mat)
  gc()
  
  # Create Seurat object
  seu <- CreateSeuratObject(
    counts = sparse_mat,
    project = paste0("GSE160269_", cell_type),
    min.cells = 0,
    min.features = 0
  )
  rm(sparse_mat)
  gc()
  
  # Add metadata
  meta <- data.frame(row.names = colnames(seu))
  meta$dataset <- "GSE160269"
  meta$technology <- "10x_5prime"
  meta$cell_type <- cell_type
  
  # Match cell metadata - Seurat may add prefix, so match by original barcodes
  # Extract original barcodes (remove prefix if present)
  orig_barcodes <- gsub("^.*_", "", colnames(seu))
  
  cell_info <- cells_all %>% filter(cell %in% orig_barcodes)
  if (nrow(cell_info) > 0) {
    meta$sample <- NA
    meta$annotated_type <- NA
    # Match by position
    match_idx <- match(orig_barcodes, cell_info$cell)
    valid_idx <- !is.na(match_idx)
    if (any(valid_idx)) {
      meta$sample[valid_idx] <- cell_info$sample[match_idx[valid_idx]]
      meta$annotated_type[valid_idx] <- cell_info$annotated_type[match_idx[valid_idx]]
    }
  }
  
  meta$patient_id <- gsub("^(P[0-9]+).*", "\\1", meta$sample)
  meta$tissue <- ifelse(grepl("T$", meta$sample), "tumor",
                        ifelse(grepl("N$", meta$sample), "normal", "unknown"))
  meta$treatment <- "Surgery_alone"
  
  seu <- AddMetaData(seu, meta)
  
  cat("    Seurat:", ncol(seu), "cells,", nrow(seu), "genes\n")
  
  seu_list[[cell_type]] <- seu
  rm(seu)
  gc()
}

# Merge all Seurat objects
cat("\n  Merging", length(seu_list), "Seurat objects...\n")
seu_160269 <- merge(seu_list[[1]], y = seu_list[-1], 
                     add.cell.ids = names(seu_list),
                     project = "GSE160269")

cat("  Merged:", ncol(seu_160269), "cells,", nrow(seu_160269), "genes\n")

# Save
saveRDS(seu_160269, file.path(OBJ_DIR, "GSE160269", "GSE160269_raw.rds"))
cat("\n  Saved: GSE160269_raw.rds\n")

# Print summary
cat("\n  Sample distribution:\n")
print(table(seu_160269$sample))

cat("\n  Cell type distribution:\n")
print(table(seu_160269$cell_type))

cat("\n=== Step 13 Complete ===\n")
cat("End:", format(Sys.time()), "\n")
