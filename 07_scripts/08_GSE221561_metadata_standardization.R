# ============================================================
# ESCC Neoadjuvant Project - STEP 11
# Standardize patient / sample / tissue / treatment / cell type
# ============================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
})

cat("\n============================================================\n")
cat("STEP 11: GSE221561 METADATA STANDARDIZATION\n")
cat("============================================================\n\n")

# Paths
object_file <- "03_objects/GSE221561/GSE221561_author_annotated.rds"
table_dir <- "06_tables/GSE221561/metadata"
log_dir <- "08_logs"

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

# Load
cat("===== 1. Loading annotated object =====\n\n")
obj <- readRDS(object_file)
md <- as.data.table(obj[[]], keep.rownames = "cell_id")
cat("Cells:", nrow(md), "| Samples:", uniqueN(md$sample_id), "| Fields:", ncol(md), "\n")

# Author fields
author_fields <- grep("^author_", names(md), value = TRUE)
cat("\n===== 2. AUTHOR METADATA FIELDS =====\n\n")
for (x in author_fields) cat(sprintf("%-45s unique=%d  missing=%d\n", x, uniqueN(md[[x]], na.rm = TRUE), sum(is.na(md[[x]]))))

# Detect fields by values
score_values <- function(column, expected) {
  vals <- tolower(trimws(as.character(md[[column]])))
  sum(tolower(expected) %in% vals)
}

celltype_expected <- c("Epithelial", "T cell", "Fibroblast", "Myeloid", "Endothelial", "B cell", "Mast cell")
treatment_expected <- c("Chemoradiotherapy", "Chemoradiotherapy and Immunotherapy", "Surgery alone", "Chemotherapy", "Adjacent normal", "Antiangiogenesis")
tissue_expected <- c("Tumor", "Adjacent normal")

detect_best_field <- function(expected) {
  scores <- data.table(field = author_fields, score = vapply(author_fields, score_values, numeric(1), expected = expected))
  setorder(scores, -score)
  scores
}

celltype_scores <- detect_best_field(celltype_expected)
treatment_scores <- detect_best_field(treatment_expected)
tissue_scores <- detect_best_field(tissue_expected)

cat("\n===== 3. FIELD DETECTION =====\n\n")
cat("Cell-type candidates:\n"); print(head(celltype_scores, 5))
cat("\nTreatment candidates:\n"); print(head(treatment_scores, 5))
cat("\nTissue candidates:\n"); print(head(tissue_scores, 5))

celltype_col <- celltype_scores$field[1]
treatment_col <- treatment_scores$field[1]
tissue_col <- tissue_scores$field[1]

cat("\nSelected cell-type:", celltype_col, "\n")
cat("Selected treatment:", treatment_col, "\n")
cat("Selected tissue:", tissue_col, "\n")

# Detect staging fields
find_field_name <- function(pattern) {
  hits <- author_fields[grepl(pattern, author_fields, ignore.case = TRUE)]
  if (length(hits) == 0) return(NA_character_)
  hits[1]
}

T_stage_col <- find_field_name("T[._]?stage")
N_stage_col <- find_field_name("N[._]?stage")
M_stage_col <- find_field_name("M[._]?stage")
response_col <- find_field_name("response|pCR|TRG|regression")

cat("\n===== 4. CLINICAL FIELD DETECTION =====\n\n")
cat("T stage:", T_stage_col, "\nN stage:", N_stage_col, "\nM stage:", M_stage_col, "\nResponse:", response_col, "\n")

# Create standardized metadata
std <- data.table(
  cell_id = md$cell_id,
  sample_id = as.character(md$sample_id),
  barcode_10x = as.character(md$barcode_10x),
  author_cell_type = as.character(md[[celltype_col]]),
  treatment_raw = as.character(md[[treatment_col]]),
  tissue_raw = as.character(md[[tissue_col]])
)

# Patient ID from sample
std[, patient_id := sub("^S([0-9]+).*$", "P\\1", sample_id)]

# Standardize tissue
std[, tissue_std := fifelse(grepl("adjacent|normal", tissue_raw, ignore.case = TRUE), "Adjacent_normal",
  fifelse(grepl("tumou?r", tissue_raw, ignore.case = TRUE), "Tumor", NA_character_))]

# Standardize treatment
std[, treatment_group := fifelse(tissue_std == "Adjacent_normal", "Adjacent_normal",
  fifelse(grepl("surgery alone", treatment_raw, ignore.case = TRUE), "Surgery_alone",
    fifelse(!is.na(treatment_raw) & trimws(treatment_raw) != "", "Neoadjuvant_treated", NA_character_)))]

# Standardize cell class
x <- tolower(trimws(std$author_cell_type))
std[, cell_class_std := fcase(
  grepl("epithelial", x), "Epithelial_unclassified",
  grepl("t cell|t-cell|nk", x), "T_NK",
  grepl("myeloid", x), "Myeloid",
  grepl("fibroblast|caf", x), "Fibroblast_CAF_candidate",
  grepl("endothelial", x), "Endothelial",
  grepl("b cell|b-cell|plasma", x), "B_plasma",
  grepl("mast", x), "Mast",
  default = "Other"
)]

# Add staging
std[, T_stage := if (!is.na(T_stage_col)) as.character(md[[T_stage_col]]) else NA_character_]
std[, N_stage := if (!is.na(N_stage_col)) as.character(md[[N_stage_col]]) else NA_character_]
std[, M_stage := if (!is.na(M_stage_col)) as.character(md[[M_stage_col]]) else NA_character_]
std[, response_raw := if (!is.na(response_col)) as.character(md[[response_col]]) else NA_character_]

# Eligibility
std[, metadata_complete := !is.na(patient_id) & !is.na(sample_id) & !is.na(tissue_std) & !is.na(treatment_group) & !is.na(cell_class_std)]

cat("\n===== 5. STANDARDIZED METADATA CHECK =====\n\n")
cat("Cells:", nrow(std), "\nComplete:", sum(std$metadata_complete), "\nIncomplete:", sum(!std$metadata_complete), "\n")

# Sample mapping
collapse_unique <- function(x) { vals <- unique(as.character(x)); vals <- vals[!is.na(vals) & trimws(vals) != ""]; if (length(vals) == 0) return(NA_character_); paste(sort(vals), collapse = "; ") }

sample_map <- std[, .(patient_id = collapse_unique(patient_id), n_cells = .N, tissue = collapse_unique(tissue_std), treatment_group = collapse_unique(treatment_group), treatment_raw = collapse_unique(treatment_raw), T_stage = collapse_unique(T_stage), N_stage = collapse_unique(N_stage), M_stage = collapse_unique(M_stage), response = collapse_unique(response_raw)), by = sample_id]
setorder(sample_map, patient_id, sample_id)

cat("\n===== 6. SAMPLE -> PATIENT MAP =====\n\n")
print(sample_map)

# Consistency check
consistency <- std[, .(n_patient_ids = uniqueN(patient_id, na.rm = TRUE), n_tissues = uniqueN(tissue_std, na.rm = TRUE), n_treatment_groups = uniqueN(treatment_group, na.rm = TRUE), n_T_stage = uniqueN(T_stage, na.rm = TRUE), n_N_stage = uniqueN(N_stage, na.rm = TRUE)), by = sample_id]
consistency[, basic_consistency := n_patient_ids == 1 & n_tissues == 1 & n_treatment_groups == 1]

cat("\n===== 7. SAMPLE CONSISTENCY =====\n\n")
print(consistency)

# Patient audit
patient_map <- std[, .(n_cells = .N, n_samples = uniqueN(sample_id), samples = collapse_unique(sample_id), tissues = collapse_unique(tissue_std), treatment_groups = collapse_unique(treatment_group), treatment_raw = collapse_unique(treatment_raw), T_stage = collapse_unique(T_stage), N_stage = collapse_unique(N_stage), M_stage = collapse_unique(M_stage)), by = patient_id]
setorder(patient_map, patient_id)

cat("\n===== 8. PATIENT-LEVEL AUDIT =====\n\n")
print(patient_map)

# Multi-sample patients
paired_patients <- patient_map[n_samples > 1]
cat("\n===== 9. MULTI-SAMPLE PATIENTS =====\n\n")
if (nrow(paired_patients) > 0) print(paired_patients) else cat("None detected.\n")

# Cell counts
cell_counts <- std[, .N, by = .(patient_id, sample_id, tissue_std, treatment_group, treatment_raw, cell_class_std)]
setnames(cell_counts, "N", "n_cells")
setorder(cell_counts, patient_id, sample_id, cell_class_std)

cat("\n===== 10. CELL COUNTS =====\n\n")
print(cell_counts)

# Pseudobulk eligibility
pb_eligibility <- cell_counts[, .(patient_id, sample_id, tissue_std, treatment_group, treatment_raw, cell_class_std, n_cells, eligible_20 = n_cells >= 20, eligible_30 = n_cells >= 30)]

cat("\n===== 11. PSEUDOBULK ELIGIBILITY =====\n\n")
print(pb_eligibility)

# Save outputs
fwrite(std, file.path(table_dir, "GSE221561_cell_metadata_standardized.csv.gz"))
fwrite(sample_map, file.path(table_dir, "GSE221561_sample_patient_mapping.csv"), bom = TRUE)
fwrite(consistency, file.path(table_dir, "GSE221561_sample_metadata_consistency.csv"), bom = TRUE)
fwrite(patient_map, file.path(table_dir, "GSE221561_patient_audit.csv"), bom = TRUE)
fwrite(cell_counts, file.path(table_dir, "GSE221561_cell_counts_by_sample_class.csv"), bom = TRUE)
fwrite(pb_eligibility, file.path(table_dir, "GSE221561_pseudobulk_cell_number_eligibility.csv"), bom = TRUE)

field_selection <- data.table(standardized_field = c("author_cell_type", "treatment_raw", "tissue_raw", "T_stage", "N_stage", "M_stage", "response_raw"), source_field = c(celltype_col, treatment_col, tissue_col, T_stage_col, N_stage_col, M_stage_col, response_col))
fwrite(field_selection, file.path(table_dir, "GSE221561_metadata_field_selection.csv"), bom = TRUE)

sink(file.path(log_dir, "08_GSE221561_metadata_standardization_sessionInfo.txt"))
print(sessionInfo())
sink()

cat("\n============================================================\n")
cat("STEP 11 FINISHED ✓\n")
cat("============================================================\n\n")
cat("No cells removed ✓\nNo normalization ✓\nNo clustering ✓\nNo malignant calls ✓\n")
cat("\nMain mapping: 06_tables/GSE221561/metadata/GSE221561_sample_patient_mapping.csv\n")
cat("Patient audit: 06_tables/GSE221561/metadata/GSE221561_patient_audit.csv\n")
cat("Pseudobulk: 06_tables/GSE221561/metadata/GSE221561_pseudobulk_cell_number_eligibility.csv\n")
