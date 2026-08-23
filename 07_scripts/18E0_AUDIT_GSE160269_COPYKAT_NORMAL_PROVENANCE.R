#!/usr/bin/env Rscript
# ==============================================================================
# 18E0: AUDIT CURRENT COPYKAT NORMAL PROVENANCE
# ==============================================================================

library(data.table)
library(Seurat)

cat("============================================================\n")
cat("STEP 18E0: COPYKAT NORMAL PROVENANCE AUDIT\n")
cat("============================================================\n\n")

OUTDIR <- "06_tables/cross_dataset/malignant_epithelial/GSE160269_copykat_validation/reference_audit"
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

BASEDIR <- "06_tables/cross_dataset/malignant_epithelial/GSE160269_copykat_validation"

# ============================================================
# 1. LOAD KNOWN NORMAL CELLS FROM STEP 18D
# ============================================================

cat("1. Loading Step 18D known normal cells...\n")

known_normal <- fread(file.path(BASEDIR, "GSE160269_COPYKAT_known_normal_cells.csv"))
cat("   Total known normal cells:", nrow(known_normal), "\n")
cat("   Cell types:", paste(unique(known_normal$cell_type), collapse = ", "), "\n\n")

# ============================================================
# 2. LOAD GSE160269 METADATA
# ============================================================

cat("2. Loading GSE160269 metadata...\n")

obj <- readRDS("03_objects/GSE160269/GSE160269_annotated_merged_raw.rds")
meta <- as.data.table(obj@meta.data, keep.rownames = "cell_barcode")

cat("   Total cells:", nrow(meta), "\n\n")

# ============================================================
# 3. TRACE PROVENANCE FOR EACH TUMOR
# ============================================================

cat("3. Tracing normal cell provenance...\n\n")

results <- list()

for (tumor in c("P16T", "P39T")) {
  cat(sprintf("--- %s ---\n", tumor))
  
  # Get annotations for this tumor's CopyKAT run
  annot_file <- file.path(BASEDIR, paste0(tumor, "_cell_annotations.csv"))
  annot <- fread(annot_file)
  
  # Get normal cells used in this tumor's CopyKAT run
  tumor_normal_cells <- annot[role == "NORMAL", cell_barcode]
  
  # Match with metadata to get provenance
  provenance <- meta[cell_barcode %in% tumor_normal_cells, 
                     .(cell_barcode, specimen_id, tissue, cell_type_author_source)]
  
  # Classify provenance
  provenance[, provenance := fifelse(
    specimen_id == tumor, "SAME_TUMOR_SAMPLE",
    fifelse(tissue == "Adjacent_normal", "ADJACENT_NORMAL_SAMPLE", "OTHER_TUMOR_SAMPLE")
  )]
  
  # Summary by specimen
  by_specimen <- provenance[, .N, by = .(specimen_id, tissue, provenance)]
  setorder(by_specimen, specimen_id)
  
  cat(sprintf("   Total normal cells: %d\n", nrow(provenance)))
  cat("\n   By specimen:\n")
  for (i in 1:nrow(by_specimen)) {
    cat(sprintf("     %s (%s): %d cells [%s]\n",
                by_specimen$specimen_id[i],
                by_specimen$tissue[i],
                by_specimen$N[i],
                by_specimen$provenance[i]))
  }
  
  # Summary by lineage
  by_lineage <- provenance[, .N, by = cell_type_author_source]
  setorder(by_lineage, -N)
  
  cat("\n   By lineage:\n")
  for (i in 1:nrow(by_lineage)) {
    cat(sprintf("     %s: %d cells\n",
                by_lineage$cell_type_author_source[i],
                by_lineage$N[i]))
  }
  
  # Check for same-tumor-sample cells
  same_sample <- provenance[provenance == "SAME_TUMOR_SAMPLE"]
  cat(sprintf("\n   Same-tumor-sample normal cells: %d\n", nrow(same_sample)))
  
  # Add tumor label
  provenance[, tumor_sample := tumor]
  results[[tumor]] <- provenance
  
  cat("\n")
}

# ============================================================
# 4. SAVE PROVENANCE TABLE
# ============================================================

cat("4. Saving provenance table...\n")

all_provenance <- rbindlist(results)
fwrite(all_provenance, file.path(OUTDIR, "STEP18E_current_known_normal_provenance.csv"))
cat("   Saved: STEP18E_current_known_normal_provenance.csv\n\n")

# ============================================================
# 5. SUMMARY
# ============================================================

cat("============================================================\n")
cat("STEP 18E0 PROVENANCE AUDIT COMPLETE\n")
cat("============================================================\n\n")

for (tumor in c("P16T", "P39T")) {
  prov <- results[[tumor]]
  
  cat(sprintf("%s:\n", tumor))
  cat(sprintf("  Total known-normal cells: %d\n", nrow(prov)))
  
  # Provenance counts
  prov_counts <- prov[, .N, by = provenance]
  for (i in 1:nrow(prov_counts)) {
    cat(sprintf("  %s: %d\n", prov_counts$provenance[i], prov_counts$N[i]))
  }
  
  # Check same-sample
  same_count <- prov[provenance == "SAME_TUMOR_SAMPLE", .N]
  cat(sprintf("  Same-sample non-epithelial available: %s\n",
              ifelse(same_count > 0, "YES", "NO")))
  
  cat("\n")
}

cat("Files saved to:", OUTDIR, "\n")
