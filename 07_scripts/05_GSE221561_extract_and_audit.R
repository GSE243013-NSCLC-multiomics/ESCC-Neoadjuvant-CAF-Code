# ============================================================
# ESCC Neoadjuvant Project
# Step 8: Extract GSE221561 and audit samples / metadata
# ============================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
})

cat("\n============================================================\n")
cat("STEP 8: GSE221561 EXTRACTION + METADATA AUDIT\n")
cat("============================================================\n\n")

# Paths
tar_file <- "01_raw/GEO/GSE221561/GSE221561_RAW.tar"
cell_meta_file <- "01_raw/GEO/GSE221561/GSE221561_metadata_celltype.csv.gz"
extract_dir <- "02_processed/GSE221561/01_10x_raw"

dir.create(extract_dir, recursive = TRUE, showWarnings = FALSE)
dir.create("00_metadata", recursive = TRUE, showWarnings = FALSE)
dir.create("08_logs", recursive = TRUE, showWarnings = FALSE)

# 1. Extract archive
cat("===== 1. Extracting RAW.tar =====\n")
if (!file.exists(tar_file)) stop("Missing archive: ", tar_file)

existing <- list.files(extract_dir, recursive = TRUE, full.names = TRUE)
if (length(existing) == 0) {
  untar(tar_file, exdir = extract_dir)
  cat("Extraction finished.\n")
} else {
  cat("Extraction directory already contains files. Skipping.\n")
}

# 2. Find all extracted files
all_files <- list.files(extract_dir, recursive = TRUE, full.names = TRUE)
all_files <- all_files[!file.info(all_files)$isdir]
cat("\nExtracted files:", length(all_files), "\n")
if (length(all_files) != 33) cat("WARNING: expected 33 files, found", length(all_files), "\n")

# 3. Classify 10x files
basename_vec <- basename(all_files)
file_type <- fifelse(grepl("barcodes\\.tsv\\.gz$", basename_vec, ignore.case = TRUE), "barcodes",
  fifelse(grepl("features\\.tsv\\.gz$", basename_vec, ignore.case = TRUE), "features",
    fifelse(grepl("matrix\\.mtx\\.gz$", basename_vec, ignore.case = TRUE), "matrix", "other")))

extract_sample_id <- function(x) {
  x <- basename(x)
  x <- sub("_barcodes\\.tsv\\.gz$", "", x, ignore.case = TRUE)
  x <- sub("_features\\.tsv\\.gz$", "", x, ignore.case = TRUE)
  x <- sub("_matrix\\.mtx\\.gz$", "", x, ignore.case = TRUE)
  sub("^.*_", "", x)
}

sample_id <- vapply(basename_vec, extract_sample_id, character(1))
file_inventory <- data.table(sample_id = sample_id, file_type = file_type, filename = basename_vec, full_path = all_files, size_MB = round(file.info(all_files)$size / 1024^2, 3))
setorder(file_inventory, sample_id, file_type)
fwrite(file_inventory, "00_metadata/04_GSE221561_extracted_file_inventory.csv", bom = TRUE)

# 4. Verify 3-file structure
triplet_check <- file_inventory[file_type != "other", .(n_files = .N, has_barcodes = any(file_type == "barcodes"), has_features = any(file_type == "features"), has_matrix = any(file_type == "matrix")), by = sample_id]
triplet_check[, complete_10x := n_files == 3 & has_barcodes & has_features & has_matrix]
setorder(triplet_check, sample_id)
fwrite(triplet_check, "00_metadata/04_GSE221561_10x_triplet_check.csv", bom = TRUE)

cat("\n===== 2. 10x SAMPLE CHECK =====\n\n")
print(triplet_check)
if (!all(triplet_check$complete_10x)) stop("At least one sample missing complete triplet.")
cat("\nAll", nrow(triplet_check), "samples contain complete 10x triplets.\n")

# 5. Read cell metadata
cat("\n===== 3. READING CELL METADATA =====\n")
meta <- fread(cell_meta_file, encoding = "UTF-8")
cat("Cells:", nrow(meta), "| Columns:", ncol(meta), "\n\n")
cat("Column names:\n")
for (i in seq_along(names(meta))) cat(sprintf("%02d. %s\n", i, names(meta)[i]))

# 6. Identify key columns
find_columns <- function(pattern) grep(pattern, names(meta), ignore.case = TRUE, value = TRUE)
candidate_columns <- list(
  sample = find_columns("sample|orig\\.ident|patient|subject|case"),
  celltype = find_columns("cell.?type|celltype|annotation|cluster"),
  treatment = find_columns("neoadjuvant|treat|therapy|nac"),
  tissue = find_columns("tissue|source|tumou?r|normal"),
  stage = find_columns("stage|T_stage|N_stage"),
  qc = find_columns("nCount|nFeature|percent\\.mt|mito")
)

cat("\n===== 4. CANDIDATE METADATA COLUMNS =====\n")
for (nm in names(candidate_columns)) {
  cat("\n", toupper(nm), ":\n", sep = "")
  if (length(candidate_columns[[nm]]) == 0) cat("  <none>\n")
  else for (x in candidate_columns[[nm]]) cat("  -", x, "\n")
}

# 7. Unique value audit
interesting <- unique(unlist(candidate_columns[c("sample", "celltype", "treatment", "tissue", "stage")]))
cat("\n===== 5. UNIQUE VALUE AUDIT =====\n")
for (col in interesting) {
  vals <- meta[[col]]
  nunique <- uniqueN(vals)
  cat("\n---", col, "| Unique:", nunique, "---\n")
  if (nunique <= 40) print(sort(table(vals, useNA = "ifany"), decreasing = TRUE))
  else print(head(sort(table(vals, useNA = "ifany"), decreasing = TRUE), 20))
}

# 8. Sample column matching
known_samples <- triplet_check$sample_id
score_sample_column <- function(col) sum(known_samples %in% unique(as.character(meta[[col]])))
sample_scores <- data.table(column = unique(candidate_columns$sample), exact_sample_matches = vapply(unique(candidate_columns$sample), score_sample_column, numeric(1)))
setorder(sample_scores, -exact_sample_matches)

cat("\n===== 6. SAMPLE COLUMN MATCHING =====\n\n")
print(sample_scores)

best_sample_col <- NA_character_
if (nrow(sample_scores) > 0 && max(sample_scores$exact_sample_matches) > 0) {
  best_sample_col <- sample_scores$column[1]
  cat("\nBest sample column:", best_sample_col, "\n")
}

# 9. Sample-level audit
if (!is.na(best_sample_col)) {
  meta[, sample_id_audit := as.character(get(best_sample_col))]
  cell_counts <- meta[, .N, by = sample_id_audit]
  setnames(cell_counts, "N", "n_cells")
  sample_audit <- merge(data.table(sample_id = known_samples), cell_counts, by.x = "sample_id", by.y = "sample_id_audit", all = TRUE)
  sample_audit[, in_10x_archive := sample_id %in% known_samples]
  sample_audit[, in_cell_metadata := !is.na(n_cells)]
  setorder(sample_audit, sample_id)
  fwrite(sample_audit, "00_metadata/04_GSE221561_sample_audit.csv", bom = TRUE)

  cat("\n===== 7. SAMPLE-LEVEL CELL COUNTS =====\n\n")
  print(sample_audit)

  # Cell-type counts
  celltype_cols <- candidate_columns$celltype
  if (length(celltype_cols) > 0) {
    preferred <- grep("^cell.?type$|celltype|annotation", celltype_cols, ignore.case = TRUE, value = TRUE)
    best_celltype_col <- if (length(preferred) > 0) preferred[1] else celltype_cols[1]
    cat("\nSelected cell-type column:", best_celltype_col, "\n")
    celltype_counts <- meta[, .N, by = c("sample_id_audit", best_celltype_col)]
    setnames(celltype_counts, "N", "n_cells")
    fwrite(celltype_counts, "00_metadata/04_GSE221561_celltype_by_sample.csv", bom = TRUE)
    cat("\n===== 8. TOP CELL TYPES =====\n\n")
    print(head(meta[, .N, by = best_celltype_col][order(-N)], 30))
  }
}

# 10. Treatment audit
treatment_cols <- candidate_columns$treatment
if (length(treatment_cols) > 0) {
  cat("\n===== 9. TREATMENT FIELD AUDIT =====\n")
  for (col in treatment_cols) { cat("\n", col, "\n"); print(sort(table(meta[[col]], useNA = "ifany"), decreasing = TRUE)) }
}

# 11. Column dictionary
column_dictionary <- data.table(column_number = seq_along(names(meta)), column_name = names(meta), class = vapply(meta, function(x) paste(class(x), collapse = ";"), character(1)), unique_values = vapply(meta, uniqueN, integer(1)), missing_values = vapply(meta, function(x) sum(is.na(x)), integer(1)))
fwrite(column_dictionary, "00_metadata/04_GSE221561_metadata_dictionary.csv", bom = TRUE)

# 12. Session info
sink("08_logs/05_GSE221561_extract_and_audit_sessionInfo.txt")
print(sessionInfo())
sink()

cat("\n============================================================\n")
cat("STEP 8 FINISHED\n")
cat("============================================================\n")
cat("\nFiles created:\n")
cat("  00_metadata/04_GSE221561_extracted_file_inventory.csv\n")
cat("  00_metadata/04_GSE221561_10x_triplet_check.csv\n")
cat("  00_metadata/04_GSE221561_sample_audit.csv\n")
cat("  00_metadata/04_GSE221561_celltype_by_sample.csv\n")
cat("  00_metadata/04_GSE221561_metadata_dictionary.csv\n")
cat("\nExtracted 10x data: 02_processed/GSE221561/01_10x_raw/\n")
