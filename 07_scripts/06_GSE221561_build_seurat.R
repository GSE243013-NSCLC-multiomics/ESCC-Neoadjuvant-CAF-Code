# ============================================================
# ESCC Neoadjuvant Project - STEP 9
# Build GSE221561 Seurat object and match author metadata
# ============================================================

options(stringsAsFactors = FALSE, timeout = 1200)

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(Matrix)
})

cat("\n============================================================\n")
cat("STEP 9: BUILD GSE221561 SEURAT OBJECT\n")
cat("============================================================\n\n")

# Paths
data_dir <- "02_processed/GSE221561/01_10x_raw"
metadata_file <- "01_raw/GEO/GSE221561/GSE221561_metadata_celltype.csv.gz"
object_dir <- "03_objects/GSE221561"

dir.create(object_dir, recursive = TRUE, showWarnings = FALSE)
dir.create("00_metadata", recursive = TRUE, showWarnings = FALSE)
dir.create("08_logs", recursive = TRUE, showWarnings = FALSE)

# Known 11 samples
samples <- c("S0520T", "S1265T", "S1315TS", "S1535T", "S2423T", "S2487T", "S6417N", "S6417T", "S6829T", "S9478N", "S9478T")
cat("Samples:", length(samples), "\n", paste(samples, collapse = ", "), "\n\n")

# Helper
find_one_file <- function(sample_id, pattern) {
  files <- list.files(data_dir, recursive = TRUE, full.names = TRUE)
  hits <- files[grepl(sample_id, basename(files), fixed = TRUE) & grepl(pattern, basename(files), ignore.case = TRUE)]
  if (length(hits) != 1) stop("Expected exactly one file for ", sample_id, " / ", pattern, "; found ", length(hits))
  hits
}

# 1. Read barcodes
cat("===== 1. READING BARCODE INVENTORIES =====\n\n")
barcode_list <- list()
for (s in samples) {
  barcode_file <- find_one_file(s, "barcodes\\.tsv\\.gz$")
  bc <- fread(barcode_file, header = FALSE)[[1]]
  bc <- as.character(bc)
  barcode_list[[s]] <- bc
  cat(sprintf("%-10s %8d barcodes\n", s, length(bc)))
}
all_barcodes <- unique(unlist(barcode_list, use.names = FALSE))
cat("\nUnique barcodes:", length(all_barcodes), "\n")

# 2. Read author metadata
cat("\n===== 2. READING AUTHOR METADATA =====\n\n")
meta <- fread(metadata_file, encoding = "UTF-8")
cat("Rows:", nrow(meta), "| Columns:", ncol(meta), "\n")

# 3. Detect sample column
cat("\n===== 3. DETECTING SAMPLE COLUMN =====\n\n")
score_sample_column <- function(x) sum(as.character(x) %in% samples, na.rm = TRUE)
sample_scores <- data.table(column = names(meta), exact_matches = vapply(meta, score_sample_column, numeric(1)))
setorder(sample_scores, -exact_matches)
print(head(sample_scores, 15))

best_sample_col <- sample_scores$column[1]
best_sample_score <- sample_scores$exact_matches[1]
if (best_sample_score == 0) { cat("\nNo column contained exact sample IDs.\n"); best_sample_col <- NA_character_
} else { cat("\nSelected sample column:", best_sample_col, "\n") }

# 4. Detect barcode column
strip_sample_from_cell <- function(x) {
  x <- as.character(x)
  for (s in samples) {
    x <- sub(paste0("^", s, "[_:.\\-]?"), "", x)
    x <- sub(paste0("[_:.\\-]?", s, "$"), "", x)
  }
  x
}
score_barcode_column <- function(x) sum(strip_sample_from_cell(x) %in% all_barcodes, na.rm = TRUE)
barcode_scores <- data.table(column = names(meta), barcode_matches = vapply(meta, score_barcode_column, numeric(1)))
setorder(barcode_scores, -barcode_matches)

cat("\n===== 4. DETECTING BARCODE COLUMN =====\n\n")
print(head(barcode_scores, 15))

best_barcode_col <- barcode_scores$column[1]
best_barcode_score <- barcode_scores$barcode_matches[1]
if (best_barcode_score == 0) stop("Could not identify barcode column.")
cat("\nSelected barcode column:", best_barcode_col, "\n")

# 5. Create standardized metadata
meta[, barcode_audit := strip_sample_from_cell(get(best_barcode_col))]
if (!is.na(best_sample_col)) {
  meta[, sample_id_audit := as.character(get(best_sample_col))]
} else {
  cell_raw <- as.character(meta[[best_barcode_col]])
  inferred <- rep(NA_character_, length(cell_raw))
  for (s in samples) { hit <- grepl(s, cell_raw, fixed = TRUE); inferred[hit & is.na(inferred)] <- s }
  meta[, sample_id_audit := inferred]
}
meta[, cell_key := paste(sample_id_audit, barcode_audit, sep = "::")]

cat("\n===== 5. METADATA KEY AUDIT =====\n\n")
cat("Rows:", nrow(meta), "\n")
cat("Missing sample IDs:", sum(is.na(meta$sample_id_audit) | !meta$sample_id_audit %in% samples), "\n")
cat("Duplicated keys:", sum(duplicated(meta$cell_key)), "\n")

# 6. Build Seurat objects
cat("\n===== 6. BUILDING 11 SEURAT OBJECTS =====\n\n")
object_list <- list()
matrix_audit <- list()

for (s in samples) {
  cat("\n--------------------------------------------\nSample:", s, "\n")
  
  counts <- ReadMtx(mtx = find_one_file(s, "matrix\\.mtx\\.gz$"), cells = find_one_file(s, "barcodes\\.tsv\\.gz$"), features = find_one_file(s, "features\\.tsv\\.gz$"), feature.column = 2, cell.column = 1, unique.features = TRUE)
  
  cat("Genes:", nrow(counts), "| Cells:", ncol(counts), "\n")
  
  original_bc <- colnames(counts)
  colnames(counts) <- paste(s, original_bc, sep = "_")
  
  obj <- CreateSeuratObject(counts = counts, project = "GSE221561", min.cells = 0, min.features = 0)
  obj$sample_id <- s
  obj$barcode_10x <- original_bc
  obj$dataset <- "GSE221561"
  
  object_keys <- paste(s, original_bc, sep = "::")
  idx <- match(object_keys, meta$cell_key)
  matched <- !is.na(idx)
  obj$metadata_matched <- matched
  
  for (nm in names(meta)) {
    if (nm %in% c("sample_id_audit", "barcode_audit", "cell_key")) next
    new_nm <- paste0("author_", make.names(nm))
    values <- rep(NA, ncol(obj))
    values[matched] <- meta[[nm]][idx[matched]]
    obj[[new_nm]] <- values
  }
  
  object_list[[s]] <- obj
  matrix_audit[[s]] <- data.table(sample_id = s, n_genes = nrow(obj), matrix_cells = ncol(obj), metadata_matched = sum(matched), metadata_unmatched = sum(!matched), match_percent = round(100 * mean(matched), 3))
  
  cat("Metadata matched:", sum(matched), "/", ncol(obj), sprintf("(%.2f%%)\n", 100 * mean(matched)))
  rm(counts); gc(verbose = FALSE)
}

# 7. Audit
audit <- rbindlist(matrix_audit)
cat("\n===== 7. MATCH AUDIT =====\n\n")
print(audit)
fwrite(audit, "00_metadata/05_GSE221561_seurat_metadata_match_audit.csv", bom = TRUE)

overall_match <- sum(audit$metadata_matched) / sum(audit$matrix_cells) * 100
cat(sprintf("\nOVERALL MATCH RATE: %.3f%%\n", overall_match))

# 8. Metadata not in matrix
matrix_keys <- unlist(lapply(samples, function(s) paste(s, barcode_list[[s]], sep = "::")), use.names = FALSE)
meta_missing <- meta[!cell_key %in% matrix_keys]
cat("Metadata rows not in matrices:", nrow(meta_missing), "\n")

# 9. Merge
cat("\n===== 8. MERGING 11 OBJECTS =====\n\n")
merged <- merge(x = object_list[[1]], y = object_list[-1], project = "GSE221561", merge.data = FALSE)
cat("Merged genes:", nrow(merged), "| Cells:", ncol(merged), "\n")

# 10. Save raw
raw_object_file <- file.path(object_dir, "GSE221561_counts_merged.rds")
saveRDS(merged, raw_object_file, compress = FALSE)
cat("\nSaved:", raw_object_file, "\n")

# 11. Save annotated
matched_cells <- colnames(merged)[merged$metadata_matched %in% TRUE]
annotated <- subset(merged, cells = matched_cells)
annotated_file <- file.path(object_dir, "GSE221561_author_annotated.rds")
saveRDS(annotated, annotated_file, compress = FALSE)
cat("Saved:", annotated_file, "\nAnnotated cells:", ncol(annotated), "\n")

# 12. Final counts
final_sample_counts <- as.data.table(table(annotated$sample_id))
setnames(final_sample_counts, c("sample_id", "n_cells"))
fwrite(final_sample_counts, "00_metadata/05_GSE221561_final_sample_cell_counts.csv", bom = TRUE)
cat("\n===== 9. FINAL SAMPLE CELL COUNTS =====\n\n")
print(final_sample_counts)

# 13. Author field inventory
author_fields <- grep("^author_", colnames(annotated[[]]), value = TRUE)
field_table <- data.table(metadata_field = author_fields, non_missing = vapply(author_fields, function(x) sum(!is.na(annotated[[x, drop = TRUE]])), numeric(1)), unique_values = vapply(author_fields, function(x) uniqueN(annotated[[x, drop = TRUE]], na.rm = TRUE), numeric(1)))
fwrite(field_table, "00_metadata/05_GSE221561_author_metadata_fields.csv", bom = TRUE)

# 14. Object size
cat("\n===== 10. OBJECT FILES =====\n\n")
for (f in c(raw_object_file, annotated_file)) {
  cat(basename(f), ":", round(file.info(f)$size / 1024^2, 2), " MB\n")
}

# 15. Session info
sink("08_logs/06_GSE221561_build_seurat_sessionInfo.txt")
print(sessionInfo())
sink()

cat("\n============================================================\n")
if (overall_match >= 99 && ncol(annotated) == nrow(meta)) {
  cat("STEP 9 FINISHED ✓\nBARCODE/METADATA MATCH IS EXCELLENT ✓\n")
} else if (overall_match >= 95) {
  cat("STEP 9 FINISHED ✓\nMATCH RATE >=95%; REVIEW UNMATCHED CELLS NEXT.\n")
} else {
  cat("STEP 9 FINISHED WITH WARNING\nMATCH RATE <95%; DO NOT START QC YET.\n")
}
cat("============================================================\n")
cat("\nIMPORTANT:\n- No cells were filtered.\n- No normalization was performed.\n- No clustering was performed.\n- Original counts were preserved.\n")
cat("\nMain object: 03_objects/GSE221561/GSE221561_author_annotated.rds\n")
