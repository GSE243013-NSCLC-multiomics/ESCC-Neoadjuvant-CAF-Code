#!/usr/bin/env Rscript
# ============================================================
# STEP 19A2: RAW COUNTS ACCESS AUDIT
# ============================================================

library(Seurat)
library(SeuratObject)
library(data.table)
library(Matrix)

cat("============================================================\n")
cat("STEP 19A2: RAW COUNTS ACCESS AUDIT\n")
cat("============================================================\n\n")

OUTDIR <- "06_tables/cross_dataset/epithelial_pseudobulk"
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

inventory <- fread(file.path(OUTDIR, "STEP19A_tumor_epithelial_sample_inventory.csv"))

results <- list()

# ============================================================
# 1. GSE221561
# ============================================================

cat("1. GSE221561: Loading...\n\n")

obj221 <- readRDS("03_objects/GSE221561/GSE221561_counts_merged.rds")
rna221 <- obj221[["RNA"]]
md221 <- obj221@meta.data

# Get tumor samples
tumor_samples_221 <- unique(md221$sample_id[md221$sample_id %in% inventory$sample[inventory$dataset == "GSE221561"]])
cat("Tumor samples in object:", paste(sort(tumor_samples_221), collapse=", "), "\n\n")

for (smp in sort(tumor_samples_221)) {
  cat(sprintf("  %s: ", smp))
  
  # Find which layer has this sample
  layer_name <- NULL
  for (lyr in Layers(rna221)) {
    layer_cells <- colnames(LayerData(rna221, layer = lyr))
    if (smp %in% md221[layer_cells, "sample_id"]) {
      layer_name <- lyr
      break
    }
  }
  
  if (is.null(layer_name)) {
    cat("LAYER NOT FOUND\n")
    next
  }
  
  # Get raw counts for this layer
  layer_counts <- LayerData(rna221, layer = layer_name)
  
  # Get epithelial cells for this sample
  sample_cells <- rownames(md221[md221$sample_id == smp, ])
  epi_cells <- sample_cells[md221[sample_cells, "author_Cell_type"] == "Epithelial"]
  
  if (length(epi_cells) == 0) {
    # Try broad_celltype or other annotation
    epi_cells <- sample_cells[md221[sample_cells, "author_Cell_type_sub"] %in% c("Epithelial", "Malignant epithelial")]
  }
  
  # Extract epithelial counts from layer
  epi_in_layer <- intersect(epi_cells, colnames(layer_counts))
  epi_raw <- layer_counts[, epi_in_layer, drop = FALSE]
  
  n_expected <- inventory[sample == smp, epithelial_cells]
  n_found <- ncol(epi_raw)
  raw_genes <- nrow(epi_raw)
  raw_count_sum <- sum(epi_raw)
  zero_lib <- sum(colSums(epi_raw) == 0)
  
  # Validate
  if (n_found == n_expected) {
    validation <- "PASS"
  } else if (abs(n_found - n_expected) <= 5) {
    validation <- "PASS"  # minor discrepancy OK
  } else {
    validation <- "MISMATCH"
  }
  
  cat(sprintf("found=%d, expected=%d, genes=%d, sum=%d, zeros=%d, %s\n",
              n_found, n_expected, raw_genes, raw_count_sum, zero_lib, validation))
  
  results[[paste0("GSE221561_", smp)]] <- data.table(
    sample = smp,
    epithelial_cells_expected = n_expected,
    epithelial_cells_found = n_found,
    raw_genes = raw_genes,
    raw_count_sum = raw_count_sum,
    zero_library_cells = zero_lib,
    counts_source = "LayerData",
    counts_layer_or_slot = layer_name,
    counts_validation = validation
  )
  
  if (validation == "MISMATCH") {
    cat("  *** MISMATCH — STOPPING ***\n")
    stop(sprintf("GSE221569 %s: expected %d but found %d epithelial cells", smp, n_expected, n_found))
  }
}

# ============================================================
# 2. GSE197677
# ============================================================

cat("\n2. GSE197677: Loading...\n\n")

obj197 <- readRDS("03_objects/GSE197677/GSE197677_integrated_broadcelltype_annotated.rds")
rna197 <- obj197@assays[["RNA"]]
md197 <- obj197@meta.data

# Get tumor samples
tumor_samples_197 <- inventory$sample[inventory$dataset == "GSE197677"]
cat("Tumor samples:", paste(sort(tumor_samples_197), collapse=", "), "\n\n")

# Get raw counts matrix (Seurat v4: counts in rna@counts)
raw_counts_197 <- rna197@counts

for (smp in sort(tumor_samples_197)) {
  cat(sprintf("  %s: ", smp))
  
  # Get cells for this sample
  sample_cells <- rownames(md197[md197$sample == smp, ])
  
  # Get epithelial cells
  epi_cells <- sample_cells[md197[sample_cells, "broad_celltype_final"] == "Epithelial"]
  
  # Extract raw counts
  epi_in_data <- intersect(epi_cells, colnames(raw_counts_197))
  epi_raw <- raw_counts_197[, epi_in_data, drop = FALSE]
  
  n_expected <- inventory[sample == smp, epithelial_cells]
  n_found <- ncol(epi_raw)
  raw_genes <- nrow(epi_raw)
  raw_count_sum <- sum(epi_raw)
  zero_lib <- sum(colSums(epi_raw) == 0)
  
  if (n_found == n_expected) {
    validation <- "PASS"
  } else if (abs(n_found - n_expected) <= 5) {
    validation <- "PASS"
  } else {
    validation <- "MISMATCH"
  }
  
  cat(sprintf("found=%d, expected=%d, genes=%d, sum=%d, zeros=%d, %s\n",
              n_found, n_expected, raw_genes, raw_count_sum, zero_lib, validation))
  
  results[[paste0("GSE197677_", smp)]] <- data.table(
    sample = smp,
    epithelial_cells_expected = n_expected,
    epithelial_cells_found = n_found,
    raw_genes = raw_genes,
    raw_count_sum = raw_count_sum,
    zero_library_cells = zero_lib,
    counts_source = "rna@counts",
    counts_layer_or_slot = "counts",
    counts_validation = validation
  )
  
  if (validation == "MISMATCH") {
    cat("  *** MISMATCH — STOPPING ***\n")
    stop(sprintf("GSE197677 %s: expected %d but found %d epithelial cells", smp, n_expected, n_found))
  }
}

# ============================================================
# 3. GSE160269 (memory-efficient: load only epithelial layer)
# ============================================================

cat("\n3. GSE160269: Loading epithelial layer only...\n\n")

# Load pre-saved metadata
md160 <- readRDS("09_tmp/gse160269_metadata.rds")

# Load full object and extract only the epithelial layer
obj160 <- readRDS("03_objects/GSE160269/GSE160269_annotated_merged_raw.rds")
rna160 <- obj160[["RNA"]]
epi_layer_name <- "counts.GSE160269.SeuratProject.SeuratProject.SeuratProject.SeuratProject.SeuratProject"
epi_counts <- LayerData(rna160, layer = epi_layer_name)
rm(obj160, rna160); gc(verbose = FALSE)

cat("Epithelial count matrix:", nrow(epi_counts), "genes x", ncol(epi_counts), "cells\n")

# Get tumor specimens
specimen_to_sample <- inventory$sample[inventory$dataset == "GSE160269"]
cat("Tumor specimens:", length(specimen_to_sample), "\n\n")

for (smp in sort(specimen_to_sample)) {
  cat(sprintf("  %s: ", smp))
  
  # Get cells for this specimen
  sample_cells <- rownames(md160[md160$specimen_id == smp, ])
  
  if (length(sample_cells) == 0) {
    cat("NO CELLS FOUND\n")
    next
  }
  
  # Get epithelial cells (GSE160269 uses "Epithelia")
  epi_cells <- sample_cells[md160[sample_cells, "source_celltype"] == "Epithelia"]
  
  if (length(epi_cells) == 0) {
    cat("0 epithelial cells found\n")
    next
  }
  
  # Extract from the epithelial count matrix
  epi_in_data <- intersect(epi_cells, colnames(epi_counts))
  epi_raw <- epi_counts[, epi_in_data, drop = FALSE]
  
  n_expected <- inventory[sample == smp, epithelial_cells]
  n_found <- ncol(epi_raw)
  raw_genes <- nrow(epi_raw)
  raw_count_sum <- sum(epi_raw)
  zero_lib <- sum(colSums(epi_raw) == 0)
  
  if (n_found == n_expected) {
    validation <- "PASS"
  } else if (abs(n_found - n_expected) <= 5) {
    validation <- "PASS"
  } else {
    validation <- "MISMATCH"
  }
  
  cat(sprintf("found=%d, expected=%d, genes=%d, sum=%d, zeros=%d, %s\n",
              n_found, n_expected, raw_genes, raw_count_sum, zero_lib, validation))
  
  results[[paste0("GSE160269_", smp)]] <- data.table(
    sample = smp,
    epithelial_cells_expected = n_expected,
    epithelial_cells_found = n_found,
    raw_genes = raw_genes,
    raw_count_sum = raw_count_sum,
    zero_library_cells = zero_lib,
    counts_source = "LayerData",
    counts_layer_or_slot = epi_layer_name,
    counts_validation = validation
  )
  
  rm(epi_raw); gc(verbose = FALSE)
  
  if (validation == "MISMATCH") {
    cat("  *** MISMATCH — STOPPING ***\n")
    stop(sprintf("GSE160269 %s: expected %d but found %d epithelial cells", smp, n_expected, n_found))
  }
}

rm(epi_counts, md160); gc(verbose = FALSE)

# ============================================================
# 4. COMBINE AND SAVE
# ============================================================

cat("\n4. Combining results...\n\n")

audit <- rbindlist(results)

# Merge with inventory
inventory[, raw_counts_available := FALSE]
inventory[, raw_gene_features := NA_integer_]
inventory[, gene_id_format := "TO_BE_CHECKED"]

for (i in 1:nrow(audit)) {
  smp <- audit[i, sample]
  inv_idx <- inventory$sample == smp
  if (any(inv_idx)) {
    inventory[inv_idx, raw_counts_available := TRUE]
    inventory[inv_idx, raw_gene_features := audit[i, raw_genes]]
    inventory[inv_idx, gene_id_format := "TO_BE_CHECKED"]
  }
}

fwrite(inventory, file.path(OUTDIR, "STEP19A_tumor_epithelial_sample_inventory.csv"))
fwrite(audit, file.path(OUTDIR, "STEP19A_raw_counts_access_audit.csv"))

cat("Audit results:\n")
cat("  Total samples audited:", nrow(audit), "\n")
cat("  PASS:", sum(audit$counts_validation == "PASS"), "\n")
cat("  MISMATCH:", sum(audit$counts_validation == "MISMATCH"), "\n")
cat("  FAIL:", sum(audit$counts_validation == "FAIL"), "\n")

# Summary by dataset
cat("\nBy dataset:\n")
for (ds in c("GSE221561", "GSE197677", "GSE160269")) {
  ds_audit <- audit[grep(ds, audit$sample)]
  cat(sprintf("  %s: %d samples, all PASS=%s\n",
              ds, nrow(ds_audit), all(ds_audit$counts_validation == "PASS")))
}

cat("\nSaved: STEP19A_raw_counts_access_audit.csv\n")
cat("Saved: STEP19A_tumor_epithelial_sample_inventory.csv (updated)\n")
