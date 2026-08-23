# ============================================================
# ESCC Neoadjuvant Project - STEP 12A
# Validate sample-patient-treatment mapping
# ============================================================

suppressPackageStartupMessages({ library(data.table) })

cat("\n============================================================\n")
cat("STEP 12A: SAMPLE / PATIENT / TREATMENT VALIDATION\n")
cat("============================================================\n\n")

meta_file <- "06_tables/GSE221561/metadata/GSE221561_cell_metadata_standardized.csv.gz"
if (!file.exists(meta_file)) stop("Cannot find: ", meta_file)
md <- fread(meta_file)
cat("Cells:", nrow(md), "| Samples:", uniqueN(md$sample_id), "| Patients:", uniqueN(md$patient_id), "\n\n")

collapse_unique <- function(x) { x <- unique(as.character(x)); x <- x[!is.na(x) & trimws(x) != ""]; if (length(x) == 0) return(NA_character_); paste(sort(x), collapse = "; ") }

# 1. Sample mapping
sample_map <- md[, .(patient_id = collapse_unique(patient_id), n_cells = .N, tissue_std = collapse_unique(tissue_std), treatment_group = collapse_unique(treatment_group), treatment_raw = collapse_unique(treatment_raw), T_stage = collapse_unique(T_stage), N_stage = collapse_unique(N_stage), M_stage = collapse_unique(M_stage)), by = sample_id]
setorder(sample_map, sample_id)
cat("===== 1. EXACT SAMPLE MAPPING =====\n\n")
print(sample_map)

# 2. Sample-patient combinations
sample_patient <- unique(md[, .(sample_id, patient_id, tissue_std, treatment_group, treatment_raw)])
setorder(sample_patient, patient_id, sample_id)
cat("\n===== 2. UNIQUE SAMPLE-PATIENT COMBINATIONS =====\n\n")
print(sample_patient)

# 3. Treatment counts
sample_treatment_summary <- sample_patient[, .(n_samples = uniqueN(sample_id), n_patients = uniqueN(patient_id), samples = paste(sort(unique(sample_id)), collapse = ", "), patients = paste(sort(unique(patient_id)), collapse = ", ")), by = treatment_group]
cat("\n===== 3. TREATMENT COUNTS =====\n\n")
print(sample_treatment_summary)

# 4. Exact treatment categories
raw_treatment_summary <- sample_patient[, .(n_samples = uniqueN(sample_id), n_patients = uniqueN(patient_id), samples = paste(sort(unique(sample_id)), collapse = ", "), patients = paste(sort(unique(patient_id)), collapse = ", ")), by = .(treatment_group, treatment_raw)]
setorder(raw_treatment_summary, treatment_group, treatment_raw)
cat("\n===== 4. EXACT TREATMENT CATEGORIES =====\n\n")
print(raw_treatment_summary)

# 5. Expected structure check
expected_neoadjuvant <- c("S0520T", "S1265T", "S1315TS", "S1535T", "S2423T", "S2487T", "S6829T")
expected_surgery_tumor <- c("S6417T", "S9478T")
expected_adjacent <- c("S6417N", "S9478N")

check_group <- function(samples, label) {
  x <- sample_patient[sample_id %in% samples]
  cat("\n", label, "\n--------------------------------------------\n")
  print(x)
  cat("Samples:", uniqueN(x$sample_id), "/", length(samples), "| Patients:", uniqueN(x$patient_id), "\n")
}

cat("\n===== 5. EXPECTED COHORT STRUCTURE =====\n")
check_group(expected_neoadjuvant, "NEOADJUVANT SAMPLES")
check_group(expected_surgery_tumor, "SURGERY-ALONE TUMOR")
check_group(expected_adjacent, "ADJACENT NORMAL")

# 6. Consistency
consistency <- md[, .(n_patient = uniqueN(patient_id), n_tissue = uniqueN(tissue_std), n_treatment_group = uniqueN(treatment_group), n_treatment_raw = uniqueN(treatment_raw)), by = sample_id]
consistency[, consistent := n_patient == 1 & n_tissue == 1 & n_treatment_group == 1 & n_treatment_raw == 1]
cat("\n===== 6. CONSISTENCY CHECK =====\n\n")
print(consistency)

# Save
outdir <- "06_tables/GSE221561/metadata"
fwrite(sample_patient, file.path(outdir, "GSE221561_sample_patient_treatment_validated.csv"), bom = TRUE)
fwrite(sample_treatment_summary, file.path(outdir, "GSE221561_treatment_group_summary_validated.csv"), bom = TRUE)
fwrite(raw_treatment_summary, file.path(outdir, "GSE221561_treatment_raw_summary_validated.csv"), bom = TRUE)

cat("\n============================================================\n")
cat("STEP 12A FINISHED\n")
cat("============================================================\n")
