# ============================================================
# ESCC Neoadjuvant Project
# STEP 18A
#
# GSE160269 CNV Gene-Layer Preflight
#
# PURPOSE:
#   Audit epithelial cells by sample/tissue
#   Confirm normal reference availability (4 adjacent-normal patients)
#   Check gene overlap with hg38 InferCNA gene table
#   Identify eligible tumor samples for CNV runs
#
# IMPORTANT:
#   NO CNV analysis yet
#   NO malignant calling yet
#   NO cells removed
# ============================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(data.table)
  library(infercna)
})

cat("\n============================================================\n")
cat("STEP 18A: GSE160269 CNV GENE-LAYER PREFLIGHT\n")
cat("============================================================\n\n")

outdir <- "06_tables/cross_dataset/malignant_epithelial/GSE160269_preflight"

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ============================================================
# 1. LOAD GSE160269 OBJECT
# ============================================================

cat("===== 1. LOAD GSE160269 =====\n\n")

obj_path <- "03_objects/GSE160269/GSE160269_annotated_merged_raw.rds"

if (!file.exists(obj_path)) {
  stop("Missing: ", obj_path)
}

obj <- readRDS(obj_path)

cat(
  "Object:",
  obj_path, "\n"
)
cat("Class:", class(obj), "\n")
cat("Cells:", ncol(obj), "\n")
cat("Features:", nrow(obj), "\n")
cat("Assays:", paste(names(obj@assays), collapse = ", "), "\n")

# ============================================================
# 2. AUDIT METADATA COLUMNS
# ============================================================

cat("\n===== 2. METADATA AUDIT =====\n\n")

required_cols <- c(
  "patient_id",
  "specimen_id",
  "tissue",
  "cell_type_author_source",
  "treatment_group"
)

missing_cols <- setdiff(
  required_cols,
  colnames(obj@meta.data)
)

if (length(missing_cols) > 0) {
  stop(
    "GSE160269 missing metadata columns: ",
    paste(missing_cols, collapse = ", ")
  )
}

cat(
  "Required columns present:",
  paste(required_cols, collapse = ", "),
  "\n"
)

# ============================================================
# 3. EPITHELIAL CELL AUDIT
# ============================================================

cat("\n===== 3. EPITHELIAL CELL AUDIT =====\n\n")

meta <- as.data.table(obj@meta.data, keep.rownames = "cell_barcode")

epi <- meta[
  cell_type_author_source == "Epithelia"
]

cat("Total epithelial cells:", nrow(epi), "\n")
cat("Tumor epithelial:", sum(epi$tissue == "Tumor"), "\n")
cat("Normal epithelial:", sum(epi$tissue == "Adjacent_normal"), "\n")

# By sample
epi_by_sample <- epi[
  ,
  .(
    cells = .N
  ),
  by = .(
    patient_id = as.character(patient_id),
    specimen_id = as.character(specimen_id),
    tissue = as.character(tissue),
    treatment_group = as.character(treatment_group)
  )
]

setorder(epi_by_sample, tissue, -cells)

cat("\n--- Epithelial cells by sample ---\n\n")
print(epi_by_sample, nrows = 100)

fwrite(
  epi_by_sample,
  file.path(outdir, "STEP18A_GSE160269_epithelial_sample_inventory.csv"),
  bom = TRUE
)

# ============================================================
# 4. NORMAL REFERENCE GROUP AUDIT
# ============================================================

cat("\n===== 4. NORMAL REFERENCE GROUP AUDIT =====\n\n")

normals <- epi_by_sample[
  tissue == "Adjacent_normal"
]

cat("Adjacent-normal samples:", nrow(normals), "\n")
cat("Total normal epithelial cells:", sum(normals$cells), "\n\n")

cat("Per-sample normal reference cells:\n")
for (i in seq_len(nrow(normals))) {
  cat(
    "  ",
    normals$patient_id[i],
    " (", normals$specimen_id[i], "): ",
    normals$cells[i],
    " cells\n",
    sep = ""
  )
}

# Classify reference adequacy
normals[
  ,
  ref_status := fifelse(
    cells >= 30,
    "ADEQUATE",
    fifelse(
      cells >= 20,
      "MARGINAL",
      "VERY_LOW"
    )
  )
]

cat("\nReference adequacy:\n")
print(normals[, .(patient_id, specimen_id, cells, ref_status)])

fwrite(
  normals,
  file.path(outdir, "STEP18A_GSE160269_normal_reference_audit.csv"),
  bom = TRUE
)

# ============================================================
# 5. HG38 GENE OVERLAP CHECK
# ============================================================

cat("\n===== 5. HG38 GENE OVERLAP =====\n\n")

hg38_dt <- as.data.table(infercna:::hg38)
hg38_genes <- hg38_dt$symbol

cat("InferCNA hg38 gene table:", length(hg38_genes), "genes\n")

# Get genes from GSE160269
obj_genes <- rownames(obj@assays$RNA)

cat("GSE160269 genes:", length(obj_genes), "\n")

# Check overlap
overlap <- intersect(obj_genes, hg38_genes)

cat("Overlap with hg38:", length(overlap), "\n")
cat(
  "Overlap rate:",
  round(100 * length(overlap) / length(hg38_genes), 1),
  "% of hg38\n"
)

# Duplicate check
dup_genes <- obj_genes[duplicated(obj_genes)]
cat("Duplicate gene symbols in GSE160269:", length(dup_genes), "\n")

if (length(dup_genes) > 0) {
  cat("First 20 duplicates:", paste(head(dup_genes, 20), collapse = ", "), "\n")
}

# ============================================================
# 6. TUMOR SAMPLE ELIGIBILITY
# ============================================================

cat("\n===== 6. TUMOR SAMPLE ELIGIBILITY =====\n\n")

tumors <- epi_by_sample[
  tissue == "Tumor"
]

tumors[
  ,
  eligible_ge30 := cells >= 30
]

tumors[
  ,
  eligible_ge100 := cells >= 100
]

cat("Tumor samples:", nrow(tumors), "\n")
cat("Tumor samples with >=30 cells:", sum(tumors$eligible_ge30), "\n")
cat("Tumor samples with >=100 cells:", sum(tumors$eligible_ge100), "\n")

# Top 10 by cell count
setorder(tumors, -cells)
cat("\nTop 10 tumor samples by epithelial cell count:\n")
print(tumors[1:min(10, nrow(tumors)), .(patient_id, specimen_id, cells, eligible_ge30, eligible_ge100)])

fwrite(
  tumors,
  file.path(outdir, "STEP18A_GSE160269_tumor_sample_eligibility.csv"),
  bom = TRUE
)

# ============================================================
# 7. SPECIMEN_ID FORMAT AUDIT
# ============================================================

cat("\n===== 7. SPECIMEN_ID FORMAT =====\n\n")

cat("All unique specimen_ids:\n")
print(sort(unique(epi_by_sample$specimen_id)))

# Check if T/N suffix is consistent
has_t_suffix <- grepl("T$", epi_by_sample$specimen_id)
has_n_suffix <- grepl("N$", epi_by_sample$specimen_id)

cat("\nSpecimens ending with T:", sum(has_t_suffix), "\n")
cat("Specimens ending with N:", sum(has_n_suffix), "\n")

# ============================================================
# 8. REFERENCE STRATEGY SUMMARY
# ============================================================

cat("\n===== 8. REFERENCE STRATEGY =====\n\n")

cat("GSE160269 CNV will use ALL adjacent-normal epithelial cells as reference.\n\n")

cat("Reference groups (4 separate, NOT pooled):\n")
for (i in seq_len(nrow(normals))) {
  cat(
    "  Group",
    i,
    ":",
    normals$specimen_id[i],
    "-",
    normals$cells[i],
    "cells\n"
  )
}

cat(
  "\nTotal reference cells:",
  sum(normals$cells),
  "\n"
)

cat(
  "\nNOTE: This is a LIMITED reference (183 cells total vs standard 1000).\n",
  "GSE160269 CNV results should be treated as LOWER CONFIDENCE.\n",
  "Sensitivity analysis will be performed after initial results.\n",
  sep = ""
)

# ============================================================
# 9. SAVE SUMMARY
# ============================================================

cat("\n===== 9. PREFLIGHT SUMMARY =====\n\n")

summary_dt <- data.table(
  metric = c(
    "total_cells",
    "total_epithelial",
    "tumor_epithelial",
    "normal_epithelial",
    "normal_samples",
    "normal_reference_cells",
    "tumor_samples",
    "tumor_samples_ge30",
    "tumor_samples_ge100",
    "hg38_genes",
    "gse160269_genes",
    "overlap_genes",
    "overlap_rate_pct",
    "reference_status"
  ),
  value = as.character(c(
    ncol(obj),
    nrow(epi),
    sum(epi$tissue == "Tumor"),
    sum(epi$tissue == "Adjacent_normal"),
    nrow(normals),
    sum(normals$cells),
    nrow(tumors),
    sum(tumors$eligible_ge30),
    sum(tumors$eligible_ge100),
    length(hg38_genes),
    length(obj_genes),
    length(overlap),
    round(100 * length(overlap) / length(hg38_genes), 1),
    "LIMITED_REFERENCE_LOWER_CONFIDENCE"
  ))
)

print(summary_dt)

fwrite(
  summary_dt,
  file.path(outdir, "STEP18A_GSE160269_preflight_summary.csv"),
  bom = TRUE
)

# ============================================================
# FINAL
# ============================================================

cat("\n============================================================\n")
cat("STEP 18A FINISHED\n")
cat("============================================================\n\n")

cat("Outputs:\n")
cat(
  "  ",
  file.path(outdir, "STEP18A_GSE160269_epithelial_sample_inventory.csv"),
  "\n",
  sep = ""
)
cat(
  "  ",
  file.path(outdir, "STEP18A_GSE160269_normal_reference_audit.csv"),
  "\n",
  sep = ""
)
cat(
  "  ",
  file.path(outdir, "STEP18A_GSE160269_tumor_sample_eligibility.csv"),
  "\n",
  sep = ""
)
cat(
  "  ",
  file.path(outdir, "STEP18A_GSE160269_preflight_summary.csv"),
  "\n",
  sep = ""
)

cat("\nNo CNV analysis run.\n")
cat("No malignant labels assigned.\n")
cat("No cells filtered.\n")
