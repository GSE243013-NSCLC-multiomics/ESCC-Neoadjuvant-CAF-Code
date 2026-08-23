#!/usr/bin/env Rscript
# ==============================================================================
# 18D1: GSE160269 COPYKAT PREFLIGHT
# Extract raw UMI counts for P16T and P39T + known normal cells
# ==============================================================================

library(Seurat)
library(data.table)

cat("============================================================\n")
cat("STEP 18D1: COPYKAT PREFLIGHT\n")
cat("============================================================\n\n")

OUTDIR <- "06_tables/cross_dataset/malignant_epithelial/GSE160269_copykat_validation"
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 1. LOAD GSE160269 SEURAT OBJECT
# ============================================================

cat("1. Loading GSE160269 Seurat object...\n")

obj <- readRDS("03_objects/GSE160269/GSE160269_annotated_merged_raw.rds")
cat("   Loaded:", dim(obj)[1], "genes x", dim(obj)[2], "cells\n\n")

# ============================================================
# 2. EXTRACT METADATA
# ============================================================

cat("2. Extracting metadata...\n")

meta <- as.data.table(obj@meta.data, keep.rownames = "cell_barcode")
cat("   Total cells:", nrow(meta), "\n")
cat("   Unique specimens:", length(unique(meta$specimen_id)), "\n\n")

# ============================================================
# 3. IDENTIFY TUMOR EPITHELIAL CELLS (P16T and P39T)
# ============================================================

cat("3. Identifying tumor epithelial cells...\n")

tumor_specimens <- c("P16T", "P39T")
tumor_cells <- list()

for (spec in tumor_specimens) {
  cells <- meta[
    specimen_id == spec &
    tissue == "Tumor" &
    cell_type_author_source == "Epithelia",
    cell_barcode
  ]
  tumor_cells[[spec]] <- cells
  cat(sprintf("   %s: %d tumor epithelial cells\n", spec, length(cells)))
}

cat("\n")

# ============================================================
# 4. IDENTIFY KNOWN NORMAL CELLS
# ============================================================

cat("4. Identifying known normal cells (non-epithelial from GSE160269)...\n")

# Get normal cell types
normal_types <- c("Tcell", "Bcell", "Myeloid", "Endothelial", "Fibroblast", "FRC", "Pericytes")
normal_cells <- list()
normal_counts <- data.table()

for (cell_type in normal_types) {
  cells <- meta[cell_type_author_source == cell_type, cell_barcode]
  if (length(cells) > 0) {
    normal_cells[[cell_type]] <- cells
    normal_counts <- rbind(normal_counts, data.table(
      cell_type = cell_type,
      total_available = length(cells)
    ))
    cat(sprintf("   %s: %d cells available\n", cell_type, length(cells)))
  }
}

cat("\n")

# ============================================================
# 5. DOWNSAMPLE NORMAL CELLS (max 500 per lineage)
# ============================================================

cat("5. Downsampling normal cells (max 500 per lineage, seed=16026918)...\n")

set.seed(16026918)
MAX_CELLS <- 500
selected_normal <- list()

for (cell_type in names(normal_cells)) {
  n_available <- length(normal_cells[[cell_type]])
  n_sample <- min(n_available, MAX_CELLS)
  
  if (n_sample < n_available) {
    selected <- sample(normal_cells[[cell_type]], n_sample)
  } else {
    selected <- normal_cells[[cell_type]]
  }
  
  selected_normal[[cell_type]] <- selected
  cat(sprintf("   %s: selected %d / %d cells\n", cell_type, n_sample, n_available))
}

# Combine all selected normal cells
all_selected_normal <- unlist(selected_normal)
cat(sprintf("\n   Total selected normal cells: %d\n", length(all_selected_normal)))

# Save known normal cell IDs
normal_dt <- rbindlist(lapply(names(selected_normal), function(ct) {
  data.table(cell_barcode = selected_normal[[ct]], cell_type = ct)
}))
fwrite(normal_dt, file.path(OUTDIR, "GSE160269_COPYKAT_known_normal_cells.csv"))
cat("   Saved: GSE160269_COPYKAT_known_normal_cells.csv\n\n")

# ============================================================
# 6. EXTRACT RAW COUNTS
# ============================================================

cat("6. Extracting raw UMI counts...\n")

# Get cells we need (tumor + selected normal)
needed_cells <- c(
  tumor_cells[["P16T"]],
  tumor_cells[["P39T"]],
  all_selected_normal
)
needed_cells <- unique(needed_cells)
cat("   Cells needed:", length(needed_cells), "\n")

# Find cells in data
cells_in_data <- intersect(needed_cells, colnames(obj))
cat("   Cells found in data:", length(cells_in_data), "\n")

# Subset the object to only needed cells
obj_sub <- subset(obj, cells = cells_in_data)
cat("   Subset object created:", dim(obj_sub)[1], "genes x", dim(obj_sub)[2], "cells\n")

# Now extract counts from subset
layers <- Layers(obj_sub, assay = "RNA")
cat("   Layers in subset:", length(layers), "\n")

# Find common genes across all layers
common_genes <- rownames(obj_sub)
for (layer in layers) {
  layer_genes <- rownames(LayerData(obj_sub, layer = layer, assay = "RNA"))
  common_genes <- intersect(common_genes, layer_genes)
}
cat("   Common genes:", length(common_genes), "\n")

# Extract and combine counts from subset
cat("   Combining layers from subset...\n")
all_counts <- list()
for (layer in layers) {
  layer_data <- LayerData(obj_sub, layer = layer, assay = "RNA")
  layer_data <- layer_data[common_genes, , drop = FALSE]
  all_counts[[layer]] <- layer_data
}

raw_counts <- do.call(cbind, all_counts)
rm(all_counts, obj_sub)
gc(verbose = FALSE)

cat("   Raw count matrix:", nrow(raw_counts), "genes x", ncol(raw_counts), "cells\n\n")

# ============================================================
# 7. PREPARE P16T DATA
# ============================================================

cat("7. Preparing P16T data...\n")

p16t_cells <- c(tumor_cells[["P16T"]], all_selected_normal)
p16t_cells_in_data <- intersect(p16t_cells, colnames(raw_counts))
cat(sprintf("   Cells in data: %d / %d\n", length(p16t_cells_in_data), length(p16t_cells)))

p16t_raw <- raw_counts[, p16t_cells_in_data, drop = FALSE]
cat(sprintf("   P16T raw count matrix: %d genes x %d cells\n", nrow(p16t_raw), ncol(p16t_raw)))

# Check for empty genes
gene_sums <- rowSums(p16t_raw)
valid_genes <- gene_sums > 0
cat(sprintf("   Genes with counts > 0: %d / %d\n", sum(valid_genes), nrow(p16t_raw)))

# Save P16T data
saveRDS(p16t_raw, file.path(OUTDIR, "P16T_raw_counts.rds"))

# Save cell annotations
p16t_annot <- data.table(
  cell_barcode = colnames(p16t_raw),
  role = ifelse(colnames(p16t_raw) %in% tumor_cells[["P16T"]], "TUMOR", "NORMAL"),
  specimen_id = sapply(colnames(p16t_raw), function(cb) {
    meta[cell_barcode == cb, specimen_id]
  }),
  cell_type = sapply(colnames(p16t_raw), function(cb) {
    meta[cell_barcode == cb, cell_type_author_source]
  })
)
fwrite(p16t_annot, file.path(OUTDIR, "P16T_cell_annotations.csv"))

cat(sprintf("   P16T tumor cells: %d\n", sum(p16t_annot$role == "TUMOR")))
cat(sprintf("   P16T normal cells: %d\n", sum(p16t_annot$role == "NORMAL")))
cat("   Saved: P16T_raw_counts.rds\n")
cat("   Saved: P16T_cell_annotations.csv\n\n")

# ============================================================
# 8. PREPARE P39T DATA
# ============================================================

cat("8. Preparing P39T data...\n")

p39t_cells <- c(tumor_cells[["P39T"]], all_selected_normal)
p39t_cells_in_data <- intersect(p39t_cells, colnames(raw_counts))
cat(sprintf("   Cells in data: %d / %d\n", length(p39t_cells_in_data), length(p39t_cells)))

p39t_raw <- raw_counts[, p39t_cells_in_data, drop = FALSE]
cat(sprintf("   P39T raw count matrix: %d genes x %d cells\n", nrow(p39t_raw), ncol(p39t_raw)))

# Check for empty genes
gene_sums <- rowSums(p39t_raw)
valid_genes <- gene_sums > 0
cat(sprintf("   Genes with counts > 0: %d / %d\n", sum(valid_genes), nrow(p39t_raw)))

# Save P39T data
saveRDS(p39t_raw, file.path(OUTDIR, "P39T_raw_counts.rds"))

# Save cell annotations
p39t_annot <- data.table(
  cell_barcode = colnames(p39t_raw),
  role = ifelse(colnames(p39t_raw) %in% tumor_cells[["P39T"]], "TUMOR", "NORMAL"),
  specimen_id = sapply(colnames(p39t_raw), function(cb) {
    meta[cell_barcode == cb, specimen_id]
  }),
  cell_type = sapply(colnames(p39t_raw), function(cb) {
    meta[cell_barcode == cb, cell_type_author_source]
  })
)
fwrite(p39t_annot, file.path(OUTDIR, "P39T_cell_annotations.csv"))

cat(sprintf("   P39T tumor cells: %d\n", sum(p39t_annot$role == "TUMOR")))
cat(sprintf("   P39T normal cells: %d\n", sum(p39t_annot$role == "NORMAL")))
cat("   Saved: P39T_raw_counts.rds\n")
cat("   Saved: P39T_cell_annotations.csv\n\n")

# ============================================================
# 9. SUMMARY
# ============================================================

cat("============================================================\n")
cat("STEP 18D1 PREFLIGHT COMPLETE\n")
cat("============================================================\n\n")

cat("P16T:\n")
cat(sprintf("  Tumor cells: %d\n", sum(p16t_annot$role == "TUMOR")))
cat(sprintf("  Normal cells: %d\n", sum(p16t_annot$role == "NORMAL")))
cat(sprintf("  Total genes: %d\n", nrow(p16t_raw)))
cat(sprintf("  Matrix dimensions: %d x %d\n\n", nrow(p16t_raw), ncol(p16t_raw)))

cat("P39T:\n")
cat(sprintf("  Tumor cells: %d\n", sum(p39t_annot$role == "TUMOR")))
cat(sprintf("  Normal cells: %d\n", sum(p39t_annot$role == "NORMAL")))
cat(sprintf("  Total genes: %d\n", nrow(p39t_raw)))
cat(sprintf("  Matrix dimensions: %d x %d\n\n", nrow(p39t_raw), ncol(p39t_raw)))

cat("Known normal cells:\n")
for (ct in names(selected_normal)) {
  cat(sprintf("  %s: %d cells\n", ct, length(selected_normal[[ct]])))
}

cat("\nFiles saved to:", OUTDIR, "\n")

rm(obj, raw_counts)
gc(verbose = FALSE)
