#!/usr/bin/env Rscript
# ============================================================
# STEP 19A3: GENE IDENTIFIER HARMONIZATION AUDIT
# ============================================================

library(data.table)

cat("============================================================\n")
cat("STEP 9A3: GENE IDENTIFIER HARMONIZATION AUDIT\n")
cat("============================================================\n\n")

OUTDIR <- "06_tables/cross_dataset/epithelial_pseudobulk"
MAPDIR <- file.path(OUTDIR, "gene_mapping")
dir.create(MAPDIR, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# HELPER: classify gene identifiers
# ============================================================

classify_genes <- function(features) {
  dt <- data.table(original_feature = features)
  
  # ENSG pattern: starts with ENSG followed by digits
  is_ensg <- grepl("^ENSG\\d+", dt$original_feature)
  
  # HGNC-like symbol: starts with a letter, contains letters/digits/dots/hyphens
  # Not an ENSG ID
  is_symbol <- !is_ensg & grepl("^[A-Z][A-Za-z0-9._-]*$", dt$original_feature)
  
  # Other: anything else (e.g., features with underscores, lowercase starts, etc.)
  is_other <- !is_ensg & !is_symbol
  
  dt[, mapping_status := fifelse(is_ensg, "ENSG_TO_SYMBOL",
                        fifelse(is_symbol, "DIRECT_SYMBOL",
                        "UNMAPPED"))]
  
  # For ENSG IDs, try to extract symbol from the feature name
  # Some ENSG features have format "ENSG00000123456.5" or "ENSG00000123456_gene_symbol"
  dt[, mapped_symbol := fifelse(is_ensg, NA_character_,
                        fifelse(is_symbol, original_feature, NA_character_))]
  
  # Check for duplicates
  dt[, duplicate_target := duplicated(dt$mapped_symbol) | duplicated(dt$mapped_symbol, fromLast = TRUE)]
  dt[duplicated(dt$mapped_symbol) | duplicated(dt$mapped_symbol, fromLast = TRUE), duplicate_target := TRUE]
  
  return(dt)
}

# ============================================================
# 1. GSE221561
# ============================================================

cat("1. GSE221561: Loading feature names...\n\n")

# Load the counts merged object to get feature names
obj221 <- readRDS("03_objects/GSE221561/GSE221561_counts_merged.rds")
features221 <- rownames(obj221[["RNA"]])
rm(obj221); gc(verbose = FALSE)

cat("  Total features:", length(features221), "\n")

mapping221 <- classify_genes(features221)

cat("  Classification:\n")
print(table(mapping221$mapping_status))
cat("  Duplicated symbols:", sum(mapping221$duplicate_target), "\n")

# For ENSG IDs, check if they have gene symbols appended
ensg_with_symbol <- mapping221[mapping_status == "ENSG_TO_SYMBOL" & grepl("_", original_feature)]
if (nrow(ensg_with_symbol) > 0) {
  cat("  ENSG features with appended symbols:", nrow(ensg_with_symbol), "\n")
  # Extract symbol from ENSG_SYMBOL format
  ensg_with_symbol[, mapped_symbol := sub("^ENSG[0-9]+[._]", "", original_feature)]
  cat("  Example extracted symbols:", paste(head(ensg_with_symbol$mapped_symbol, 5), collapse=", "), "\n")
}

# Mark keep_for_cross_dataset
mapping221[, keep_for_cross_dataset := mapping_status == "DIRECT_SYMBOL" & !duplicate_target]

fwrite(mapping221, file.path(MAPDIR, "GSE221561_gene_mapping.csv"))
cat("\n  Saved: GSE221561_gene_mapping.csv\n\n")

# ============================================================
# 2. GSE197677
# ============================================================

cat("2. GSE197677: Loading feature names...\n\n")

obj197 <- readRDS("03_objects/GSE197677/GSE197677_integrated_broadcelltype_annotated.rds")
features197 <- rownames(obj197@assays[["RNA"]])
rm(obj197); gc(verbose = FALSE)

cat("  Total features:", length(features197), "\n")

mapping197 <- classify_genes(features197)

cat("  Classification:\n")
print(table(mapping197$mapping_status))
cat("  Duplicated symbols:", sum(mapping197$duplicate_target), "\n")

mapping197[, keep_for_cross_dataset := mapping_status == "DIRECT_SYMBOL" & !duplicate_target]

fwrite(mapping197, file.path(MAPDIR, "GSE197677_gene_mapping.csv"))
cat("\n  Saved: GSE197677_gene_mapping.csv\n\n")

# ============================================================
# 3. GSE160269
# ============================================================

cat("3. GSE160269: Loading feature names...\n\n")

# Load from the saved epithelial counts
epi_counts <- readRDS("09_tmp/gse160269_epithelial_counts.rds")
features160 <- rownames(epi_counts)
rm(epi_counts); gc(verbose = FALSE)

cat("  Total features:", length(features160), "\n")

mapping160 <- classify_genes(features160)

cat("  Classification:\n")
print(table(mapping160$mapping_status))
cat("  Duplicated symbols:", sum(mapping160$duplicate_target), "\n")

mapping160[, keep_for_cross_dataset := mapping_status == "DIRECT_SYMBOL" & !duplicate_target]

fwrite(mapping160, file.path(MAPDIR, "GSE160269_gene_mapping.csv"))
cat("\n  Saved: GSE160269_gene_mapping.csv\n\n")

# ============================================================
# 4. SUMMARY
# ============================================================

cat("============================================================\n")
cat("GENE HARMONIZATION SUMMARY\n")
cat("============================================================\n\n")

for (ds in c("GSE221561", "GSE197677", "GSE160269")) {
  m <- get(paste0("mapping", substr(ds, 4, 8)))
  cat(sprintf("%s:\n", ds))
  cat(sprintf("  Total features: %d\n", nrow(m)))
  cat(sprintf("  DIRECT_SYMBOL: %d\n", sum(m$mapping_status == "DIRECT_SYMBOL")))
  cat(sprintf("  ENSG_TO_SYMBOL: %d\n", sum(m$mapping_status == "ENSG_TO_SYMBOL")))
  cat(sprintf("  UNMAPPED: %d\n", sum(m$mapping_status == "UNMAPPED")))
  cat(sprintf("  Duplicated: %d\n", sum(m$duplicate_target)))
  cat(sprintf("  Keep for cross-dataset: %d\n\n", sum(m$keep_for_cross_dataset)))
}
