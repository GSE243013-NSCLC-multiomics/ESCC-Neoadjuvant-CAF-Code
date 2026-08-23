# ============================================================
# ESCC Neoadjuvant Project
# STEP 13E
# Add validated patient/specimen/tissue metadata
# back to the 8 GSE160269 Seurat part objects
#
# NO MERGE
# NO QC
# NO NORMALIZATION
# NO INTEGRATION
# ============================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
})

cat("\n============================================================\n")
cat("STEP 13E: ADD VALIDATED METADATA TO GSE160269 PARTS\n")
cat("============================================================\n\n")

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------

input_dir <- "03_objects/GSE160269/parts"

output_dir <- "03_objects/GSE160269/annotated_parts"

mapping_file <- paste0(
  "06_tables/GSE160269/metadata/",
  "GSE160269_cell_index_patient_mapping.csv.gz"
)

audit_dir <- "06_tables/GSE160269/metadata"

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  audit_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# 1. Load validated mapping
# ------------------------------------------------------------

cat("===== 1. LOADING VALIDATED CELL MAPPING =====\n\n")

if (!file.exists(mapping_file)) {
  stop("Cannot find mapping file: ", mapping_file)
}

mapping <- fread(
  mapping_file
)

required_fields <- c(
  "cell_id",
  "original_barcode",
  "cell_type_source",
  "patient_id",
  "specimen_id",
  "tissue_std",
  "compartment_std",
  "barcode_core",
  "technical_cell_key"
)

missing_fields <- setdiff(
  required_fields,
  names(mapping)
)

if (length(missing_fields) > 0) {
  stop(
    "Mapping file is missing fields: ",
    paste(missing_fields, collapse = ", ")
  )
}

cat("Mapping rows:", nrow(mapping), "\n")

if (nrow(mapping) != 208659) {
  stop(
    "Expected 208659 mapping rows, found ",
    nrow(mapping)
  )
}

cat(
  "Unique cells:",
  uniqueN(mapping$cell_id),
  "\n"
)

if (uniqueN(mapping$cell_id) != 208659) {
  stop("cell_id is not globally unique.")
}

if (anyDuplicated(mapping$technical_cell_key)) {
  stop("Duplicate technical_cell_key detected.")
}

cat("✓ Mapping table validated\n")

# ------------------------------------------------------------
# 2. Find 8 original part objects
# ------------------------------------------------------------

part_files <- sort(
  list.files(
    input_dir,
    pattern = "\\.rds$",
    full.names = TRUE
  )
)

if (length(part_files) != 8) {
  stop(
    "Expected 8 part objects, found ",
    length(part_files)
  )
}

cat("\n===== 2. PART OBJECTS =====\n\n")

for (f in part_files) {
  cat(" - ", basename(f), "\n", sep = "")
}

# ------------------------------------------------------------
# 3. Process one object at a time
# ------------------------------------------------------------

audit_list <- list()

for (f in part_files) {

  cat("\n============================================================\n")
  cat("Processing: ", basename(f), "\n", sep = "")
  cat("============================================================\n")

  outfile <- file.path(
    output_dir,
    basename(f)
  )

  # Resume support
  if (file.exists(outfile)) {

    cat("✓ Annotated object already exists — skipping\n")

    existing <- readRDS(outfile)

    audit_list[[basename(f)]] <- data.table(
      object = basename(f),
      cells = ncol(existing),
      matched = sum(
        !is.na(existing$patient_id)
      ),
      unmatched = sum(
        is.na(existing$patient_id)
      ),
      patients = uniqueN(
        existing$patient_id,
        na.rm = TRUE
      ),
      specimens = uniqueN(
        existing$specimen_id,
        na.rm = TRUE
      ),
      status = "existing"
    )

    rm(existing)
    gc(verbose = FALSE)

    next
  }

  # ----------------------------------------------------------
  # Load part
  # ----------------------------------------------------------

  obj <- readRDS(f)

  cell_ids <- colnames(obj)

  cat("Cells:", length(cell_ids), "\n")

  # ----------------------------------------------------------
  # Match by globally unique cell_id
  # ----------------------------------------------------------

  idx <- match(
    cell_ids,
    mapping$cell_id
  )

  matched <- !is.na(idx)

  cat(
    "Matched:",
    sum(matched),
    "/",
    length(cell_ids),
    sprintf(
      " (%.3f%%)\n",
      100 * mean(matched)
    )
  )

  if (!all(matched)) {

    unmatched_ids <- cell_ids[
      !matched
    ]

    fwrite(
      data.table(
        object = basename(f),
        cell_id = unmatched_ids
      ),
      file.path(
        audit_dir,
        paste0(
          "GSE160269_unmatched_",
          sub("\\.rds$", "", basename(f)),
          ".csv"
        )
      ),
      bom = TRUE
    )

    stop(
      "Not all cells matched for ",
      basename(f),
      ". STOPPING before saving."
    )
  }

  # ----------------------------------------------------------
  # Safety check:
  # original_barcode must agree with mapping
  # ----------------------------------------------------------

  object_barcode <- as.character(
    obj$original_barcode
  )

  mapping_barcode <- mapping$original_barcode[
    idx
  ]

  barcode_agreement <- (
    object_barcode ==
    mapping_barcode
  )

  if (!all(barcode_agreement)) {

    stop(
      "original_barcode mismatch in ",
      basename(f)
    )
  }

  cat("✓ original_barcode agrees 100%\n")

  # ----------------------------------------------------------
  # Add validated metadata
  # ----------------------------------------------------------

  obj$patient_id <- mapping$patient_id[
    idx
  ]

  obj$specimen_id <- mapping$specimen_id[
    idx
  ]

  obj$tissue <- mapping$tissue_std[
    idx
  ]

  obj$compartment_code <- mapping$compartment_std[
    idx
  ]

  obj$barcode_core <- mapping$barcode_core[
    idx
  ]

  obj$technical_cell_key <-
    mapping$technical_cell_key[
      idx
    ]

  obj$cell_type_author_source <-
    mapping$cell_type_source[
      idx
    ]

  obj$dataset <- "GSE160269"

  # Untreated baseline cohort
  obj$treatment_group <- "Untreated_baseline"

  # ----------------------------------------------------------
  # Internal consistency checks
  # ----------------------------------------------------------

  if (any(is.na(obj$patient_id))) {
    stop("Missing patient_id after mapping.")
  }

  if (any(is.na(obj$specimen_id))) {
    stop("Missing specimen_id after mapping.")
  }

  if (any(is.na(obj$tissue))) {
    stop("Missing tissue after mapping.")
  }

  if (anyDuplicated(obj$technical_cell_key)) {
    stop(
      "Duplicate technical_cell_key within ",
      basename(f)
    )
  }

  # ----------------------------------------------------------
  # Summary
  # ----------------------------------------------------------

  n_patients <- uniqueN(
    obj$patient_id
  )

  n_specimens <- uniqueN(
    obj$specimen_id
  )

  cat("Patients :", n_patients, "\n")
  cat("Specimens:", n_specimens, "\n")

  cat("Tissue counts:\n")
  print(
    table(
      obj$tissue
    )
  )

  # ----------------------------------------------------------
  # Save immediately
  # ----------------------------------------------------------

  saveRDS(
    obj,
    outfile,
    compress = FALSE
  )

  size_mb <- file.info(
    outfile
  )$size / 1024^2

  cat(
    "✓ Saved: ",
    outfile,
    "\n",
    sep = ""
  )

  cat(
    sprintf(
      "Object size: %.2f MB\n",
      size_mb
    )
  )

  audit_list[[basename(f)]] <- data.table(
    object = basename(f),
    cells = ncol(obj),
    matched = sum(matched),
    unmatched = sum(!matched),
    patients = n_patients,
    specimens = n_specimens,
    tumor_cells = sum(
      obj$tissue == "Tumor"
    ),
    adjacent_normal_cells = sum(
      obj$tissue == "Adjacent_normal"
    ),
    size_MB = round(
      size_mb,
      2
    ),
    status = "built"
  )

  rm(obj)
  gc(verbose = FALSE)
}

# ------------------------------------------------------------
# 4. Combined audit
# ------------------------------------------------------------

audit <- rbindlist(
  audit_list,
  fill = TRUE
)

cat("\n============================================================\n")
cat("3. ANNOTATED PART AUDIT\n")
cat("============================================================\n\n")

print(audit)

fwrite(
  audit,
  file.path(
    audit_dir,
    "GSE160269_annotated_part_audit.csv"
  ),
  bom = TRUE
)

# ------------------------------------------------------------
# 5. Verify all annotated objects together
# ------------------------------------------------------------

annotated_files <- sort(
  list.files(
    output_dir,
    pattern = "\\.rds$",
    full.names = TRUE
  )
)

if (length(annotated_files) != 8) {
  stop(
    "Expected 8 annotated part objects, found ",
    length(annotated_files)
  )
}

total_cells <- 0

patient_pool <- character(0)
specimen_pool <- character(0)

tumor_specimens <- character(0)
normal_specimens <- character(0)

for (f in annotated_files) {

  obj <- readRDS(f)

  total_cells <- total_cells + ncol(obj)

  patient_pool <- c(
    patient_pool,
    obj$patient_id
  )

  specimen_pool <- c(
    specimen_pool,
    obj$specimen_id
  )

  tumor_specimens <- c(
    tumor_specimens,
    obj$specimen_id[
      obj$tissue == "Tumor"
    ]
  )

  normal_specimens <- c(
    normal_specimens,
    obj$specimen_id[
      obj$tissue == "Adjacent_normal"
    ]
  )

  rm(obj)
  gc(verbose = FALSE)
}

n_patients <- uniqueN(
  patient_pool
)

n_specimens <- uniqueN(
  specimen_pool
)

n_tumors <- uniqueN(
  tumor_specimens
)

n_normals <- uniqueN(
  normal_specimens
)

cat("\n============================================================\n")
cat("FINAL VALIDATION\n")
cat("============================================================\n\n")

cat(
  "Annotated objects :",
  length(annotated_files),
  "\n"
)

cat(
  "Total cells       :",
  total_cells,
  "\n"
)

cat(
  "Patients          :",
  n_patients,
  "\n"
)

cat(
  "Specimens         :",
  n_specimens,
  "\n"
)

cat(
  "Tumor specimens   :",
  n_tumors,
  "\n"
)

cat(
  "Adjacent normals  :",
  n_normals,
  "\n"
)

cat("\nExpected:\n")
cat("  8 annotated objects\n")
cat("  208659 cells\n")
cat("  60 patients\n")
cat("  64 specimens\n")
cat("  60 tumor specimens\n")
cat("  4 adjacent normals\n")

# ------------------------------------------------------------
# 6. Final decision
# ------------------------------------------------------------

cat("\n============================================================\n")

if (
  length(annotated_files) == 8 &&
  total_cells == 208659 &&
  n_patients == 60 &&
  n_specimens == 64 &&
  n_tumors == 60 &&
  n_normals == 4
) {

  cat("STEP 13E FINISHED ✓\n")
  cat("ALL METADATA WRITTEN BACK SUCCESSFULLY ✓\n")

} else {

  cat("STEP 13E FINISHED WITH REVIEW FLAG ⚠\n")
  cat("DO NOT MERGE YET.\n")
}

cat("============================================================\n\n")

cat("Annotated objects saved in:\n")
cat(
  "03_objects/GSE160269/annotated_parts/\n"
)

cat("\nIMPORTANT:\n")
cat("No merge performed ✓\n")
cat("No cells filtered ✓\n")
cat("No QC performed ✓\n")
cat("No normalization performed ✓\n")
cat("No integration performed ✓\n")
cat("STEP 14 NOT RUN ✓\n")

cat("\nSTOP HERE.\n")

