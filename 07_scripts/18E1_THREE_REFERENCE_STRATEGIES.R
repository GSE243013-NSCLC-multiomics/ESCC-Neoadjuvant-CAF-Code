#!/usr/bin/env Rscript
# ==============================================================================
# 18E1: DEFINE THREE COPYKAT REFERENCE STRATEGIES
# ==============================================================================

library(data.table)
library(Seurat)

cat("============================================================\n")
cat("STEP 18E1: DEFINE THREE REFERENCE STRATEGIES\n")
cat("============================================================\n\n")

OUTDIR <- "06_tables/cross_dataset/malignant_epithelial/GSE160269_copykat_validation/reference_audit"
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

BASEDIR <- "06_tables/cross_dataset/malignant_epithelial/GSE160269_copykat_validation"

# ============================================================
# 1. LOAD METADATA
# ============================================================

cat("1. Loading GSE160269 metadata...\n")

obj <- readRDS("03_objects/GSE160269/GSE160269_annotated_merged_raw.rds")
meta <- as.data.table(obj@meta.data, keep.rownames = "cell_barcode")
rm(obj)
gc(verbose = FALSE)

cat("   Loaded\n\n")

# ============================================================
# 2. DEFINE STRATEGIES
# ============================================================

cat("2. Defining reference strategies...\n\n")

strategies <- list()

# -----------------------------------------------------------
# STRATEGY A: CURRENT_18D (frozen baseline)
# -----------------------------------------------------------
cat("--- STRATEGY A: CURRENT_18D ---\n")
cat("Using Step 18D known-normal cells (frozen baseline)\n")

known_normal_18d <- fread(file.path(BASEDIR, "GSE160269_COPYKAT_known_normal_cells.csv"))
cat(sprintf("   Total cells: %d\n", nrow(known_normal_18d)))
cat(sprintf("   Cell types: %s\n\n", paste(unique(known_normal_18d$cell_type), collapse = ", ")))

strategies[["CURRENT_18D"]] <- known_normal_18d

# -----------------------------------------------------------
# STRATEGY B: SAME_SAMPLE_NON_EPITHELIAL
# -----------------------------------------------------------
cat("--- STRATEGY B: SAME_SAMPLE_NON_EPITHELIAL ---\n")

for (tumor in c("P16T", "P39T")) {
  cat(sprintf("\n%s:\n", tumor))
  
  # Get non-epithelial cells from same specimen
  same_sample_normals <- meta[
    specimen_id == tumor &
    cell_type_author_source != "Epithelia" &
    cell_type_author_source != "",
    .(cell_barcode, cell_type = cell_type_author_source)
  ]
  
  cat(sprintf("   Available same-sample non-epithelial: %d\n", nrow(same_sample_normals)))
  
  if (nrow(same_sample_normals) < 100) {
    cat(sprintf("   SAME_SAMPLE_REFERENCE_INSUFFICIENT (%d < 100)\n", nrow(same_sample_normals)))
    strategies[[paste0(tumor, "_SAME_SAMPLE_NON_EPITHELIAL")]] <- "INSUFFICIENT"
  } else {
    # Downsample to max 300 per lineage
    set.seed(16026918)
    by_lineage <- same_sample_normals[, .N, by = cell_type]
    cat("   Lineage counts:\n")
    for (i in 1:nrow(by_lineage)) {
      cat(sprintf("     %s: %d\n", by_lineage$cell_type[i], by_lineage$N[i]))
    }
    
    # Sample max 300 per lineage
    sampled <- same_sample_normals[, .SD[sample(.N, min(.N, 300))], by = cell_type]
    cat(sprintf("   Sampled: %d cells\n", nrow(sampled)))
    
    strategies[[paste0(tumor, "_SAME_SAMPLE_NON_EPITHELIAL")]] <- sampled
  }
}

# -----------------------------------------------------------
# STRATEGY C: ADJACENT_NORMAL_EPITHELIAL
# -----------------------------------------------------------
cat("\n--- STRATEGY C: ADJACENT_NORMAL_EPITHELIAL ---\n")

adj_normal_epithelial <- meta[
  tissue == "Adjacent_normal" &
  cell_type_author_source == "Epithelia",
  .(cell_barcode, specimen_id)
]

cat(sprintf("   Total adjacent-normal epithelial cells: %d\n", nrow(adj_normal_epithelial)))
cat("   By specimen:\n")
by_spec <- adj_normal_epithelial[, .N, by = specimen_id]
for (i in 1:nrow(by_spec)) {
  cat(sprintf("     %s: %d\n", by_spec$specimen_id[i], by_spec$N[i]))
}

strategies[["ADJACENT_NORMAL_EPITHELIAL"]] <- adj_normal_epithelial

# ============================================================
# 3. SAVE STRATEGIES
# ============================================================

cat("\n3. Saving reference strategies...\n")

for (name in names(strategies)) {
  if (is.data.table(strategies[[name]])) {
    fwrite(strategies[[name]], file.path(OUTDIR, paste0("reference_", name, ".csv")))
    cat(sprintf("   Saved: reference_%s.csv\n", name))
  } else {
    cat(sprintf("   %s: %s (not saved)\n", name, strategies[[name]]))
  }
}

# ============================================================
# 4. SUMMARY
# ============================================================

cat("\n============================================================\n")
cat("STEP 18E1 REFERENCE STRATEGIES DEFINED\n")
cat("============================================================\n\n")

cat("STRATEGY A: CURRENT_18D\n")
cat("  3500 cells from across all GSE160269 samples\n")
cat("  500 per lineage (7 lineages)\n\n")

cat("STRATEGY B: SAME_SAMPLE_NON_EPITHELIAL\n")
for (tumor in c("P16T", "P39T")) {
  key <- paste0(tumor, "_SAME_SAMPLE_NON_EPITHELIAL")
  val <- strategies[[key]]
  if (is.character(val)) {
    cat(sprintf("  %s: %s\n", tumor, val))
  } else {
    cat(sprintf("  %s: %d cells\n", tumor, nrow(val)))
  }
}

cat("\nSTRATEGY C: ADJACENT_NORMAL_EPITHELIAL\n")
cat(sprintf("  %d cells from P126N/P127N/P128N/P130N\n", nrow(adj_normal_epithelial)))

cat("\nFiles saved to:", OUTDIR, "\n")
