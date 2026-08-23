# ============================================================
# ESCC Neoadjuvant Project
# STEP 13D
# Parse GSE160269 patient / tissue / compartment from barcodes
#
# AUDIT ONLY
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
cat("STEP 13D: GSE160269 PATIENT / TISSUE MAPPING AUDIT\n")
cat("============================================================\n\n")

obj_dir <- "03_objects/GSE160269/parts"
out_dir <- "06_tables/GSE160269/metadata"

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# 1. Read only metadata from 8 objects
# ------------------------------------------------------------

files <- sort(
  list.files(
    obj_dir,
    pattern = "\\.rds$",
    full.names = TRUE
  )
)

if (length(files) != 8) {
  stop(
    "Expected 8 GSE160269 part objects, found ",
    length(files)
  )
}

index_list <- list()

cat("===== 1. READING 8 OBJECTS =====\n\n")

for (f in files) {

  obj <- readRDS(f)

  if (!"original_barcode" %in% colnames(obj[[]])) {
    stop(
      "Missing original_barcode in ",
      basename(f)
    )
  }

  tmp <- data.table(
    cell_id = colnames(obj),

    original_barcode =
      as.character(
        obj$original_barcode
      ),

    cell_type_source =
      as.character(
        obj$source_celltype
      )
  )

  index_list[[basename(f)]] <- tmp

  cat(
    sprintf(
      "%-45s %8d cells\n",
      basename(f),
      nrow(tmp)
    )
  )

  rm(obj)
  gc(verbose = FALSE)
}

cells <- rbindlist(
  index_list,
  fill = TRUE
)

cat(
  "\nTotal cells:",
  nrow(cells),
  "\n"
)

if (nrow(cells) != 208659) {
  stop(
    "Expected 208659 cells, found ",
    nrow(cells)
  )
}

# ------------------------------------------------------------
# 2. Parse barcode
#
# Expected:
# P1T-I-AAAA...
# P1T-E-AAAA...
#
# Components:
# P1 = patient
# T  = tissue
# I/E = compartment/library code
# ------------------------------------------------------------

pattern <- "^([Pp][0-9]+)([TtNn])-([IiEe])-(.+)$"

m <- regexec(
  pattern,
  cells$original_barcode
)

parts <- regmatches(
  cells$original_barcode,
  m
)

parsed <- rbindlist(
  lapply(
    parts,
    function(x) {

      if (length(x) == 5) {

        data.table(
          patient_code = toupper(x[2]),
          tissue_code = toupper(x[3]),
          compartment_code = toupper(x[4]),
          barcode_core = x[5],
          parse_success = TRUE
        )

      } else {

        data.table(
          patient_code = NA_character_,
          tissue_code = NA_character_,
          compartment_code = NA_character_,
          barcode_core = NA_character_,
          parse_success = FALSE
        )
      }
    }
  )
)

cells <- cbind(
  cells,
  parsed
)

# ------------------------------------------------------------
# 3. Standardize fields
# ------------------------------------------------------------

cells[
  ,
  patient_id :=
    patient_code
]

cells[
  ,
  tissue_std :=
    fcase(
      tissue_code == "T",
      "Tumor",

      tissue_code == "N",
      "Adjacent_normal",

      default = NA_character_
    )
]

cells[
  ,
  specimen_id :=
    paste0(
      patient_code,
      tissue_code
    )
]

cells[
  ,
  compartment_std :=
    fcase(
      compartment_code == "I",
      "I",

      compartment_code == "E",
      "E",

      default = NA_character_
    )
]

# ------------------------------------------------------------
# 4. Parsing success
# ------------------------------------------------------------

cat("\n===== 2. BARCODE PARSING =====\n\n")

cat(
  "Successfully parsed:",
  sum(cells$parse_success),
  "/",
  nrow(cells),
  "\n"
)

cat(
  sprintf(
    "Parse rate: %.3f%%\n",
    100 * mean(cells$parse_success)
  )
)

if (!all(cells$parse_success)) {

  bad <- cells[
    parse_success == FALSE
  ]

  fwrite(
    bad,
    file.path(
      out_dir,
      "GSE160269_unparsed_barcodes.csv"
    ),
    bom = TRUE
  )

  cat(
    "\nWARNING: unparsed barcodes detected.\n"
  )

  print(
    head(
      bad,
      30
    )
  )
}

# ------------------------------------------------------------
# 5. Overall cohort audit
# ------------------------------------------------------------

cat("\n===== 3. COHORT STRUCTURE =====\n\n")

cat(
  "Unique patients:",
  uniqueN(cells$patient_id),
  "\n"
)

cat(
  "Unique specimens:",
  uniqueN(cells$specimen_id),
  "\n"
)

cat(
  "Tumor specimens:",
  uniqueN(
    cells[
      tissue_std == "Tumor",
      specimen_id
    ]
  ),
  "\n"
)

cat(
  "Adjacent-normal specimens:",
  uniqueN(
    cells[
      tissue_std == "Adjacent_normal",
      specimen_id
    ]
  ),
  "\n"
)

# ------------------------------------------------------------
# 6. Patient/specimen table
# ------------------------------------------------------------

specimen_table <- cells[
  ,
  .(
    n_cells = .N,

    cell_types =
      paste(
        sort(
          unique(
            cell_type_source
          )
        ),
        collapse = "; "
      ),

    compartments =
      paste(
        sort(
          unique(
            compartment_std
          )
        ),
        collapse = "; "
      )
  ),
  by = .(
    patient_id,
    specimen_id,
    tissue_std
  )
]

# Natural patient number for sorting
specimen_table[
  ,
  patient_number :=
    as.integer(
      sub(
        "^P",
        "",
        patient_id
      )
    )
]

setorder(
  specimen_table,
  patient_number,
  tissue_std
)

specimen_table[
  ,
  patient_number := NULL
]

cat("\n===== 4. PATIENT / SPECIMEN TABLE =====\n\n")

print(specimen_table)

# ------------------------------------------------------------
# 7. Patient-level audit
# ------------------------------------------------------------

patient_table <- cells[
  ,
  .(
    n_cells = .N,

    n_specimens =
      uniqueN(
        specimen_id
      ),

    specimens =
      paste(
        sort(
          unique(
            specimen_id
          )
        ),
        collapse = "; "
      ),

    tissues =
      paste(
        sort(
          unique(
            tissue_std
          )
        ),
        collapse = "; "
      )
  ),
  by = patient_id
]

patient_table[
  ,
  patient_number :=
    as.integer(
      sub(
        "^P",
        "",
        patient_id
      )
    )
]

setorder(
  patient_table,
  patient_number
)

patient_table[
  ,
  patient_number := NULL
]

cat("\n===== 5. PATIENT-LEVEL AUDIT =====\n\n")

print(patient_table)

# ------------------------------------------------------------
# 8. Patients with paired tumor + normal
# ------------------------------------------------------------

paired <- patient_table[
  grepl(
    "Tumor",
    tissues,
    fixed = TRUE
  ) &
  grepl(
    "Adjacent_normal",
    tissues,
    fixed = TRUE
  )
]

cat("\n===== 6. PAIRED TUMOR / NORMAL PATIENTS =====\n\n")

if (nrow(paired) == 0) {

  cat(
    "No paired tumor/normal patients detected.\n"
  )

} else {

  print(paired)
}

cat(
  "\nNumber of paired patients:",
  nrow(paired),
  "\n"
)

# ------------------------------------------------------------
# 9. Cell type x tissue
# ------------------------------------------------------------

celltype_tissue <- cells[
  ,
  .N,
  by = .(
    cell_type_source,
    tissue_std
  )
]

setnames(
  celltype_tissue,
  "N",
  "n_cells"
)

setorder(
  celltype_tissue,
  cell_type_source,
  tissue_std
)

cat("\n===== 7. CELL TYPE x TISSUE =====\n\n")

print(celltype_tissue)

# ------------------------------------------------------------
# 10. Cell type x compartment
# ------------------------------------------------------------

celltype_compartment <- cells[
  ,
  .N,
  by = .(
    cell_type_source,
    compartment_std
  )
]

setnames(
  celltype_compartment,
  "N",
  "n_cells"
)

setorder(
  celltype_compartment,
  cell_type_source,
  compartment_std
)

cat("\n===== 8. CELL TYPE x I/E COMPARTMENT =====\n\n")

print(celltype_compartment)

# ------------------------------------------------------------
# 11. Per specimen / cell type cell counts
# ------------------------------------------------------------

pb_table <- cells[
  ,
  .N,
  by = .(
    patient_id,
    specimen_id,
    tissue_std,
    cell_type_source
  )
]

setnames(
  pb_table,
  "N",
  "n_cells"
)

pb_table[
  ,
  eligible_20 :=
    n_cells >= 20
]

pb_table[
  ,
  eligible_30 :=
    n_cells >= 30
]

setorder(
  pb_table,
  patient_id,
  specimen_id,
  cell_type_source
)

cat("\n===== 9. PSEUDOBULK CELL-NUMBER AUDIT =====\n\n")

print(pb_table)

# ------------------------------------------------------------
# 12. Duplicate core barcode audit
#
# Same barcode sequence may legitimately appear in different
# specimens. specimen + compartment + barcode must be unique.
# ------------------------------------------------------------

cells[
  ,
  technical_cell_key :=
    paste(
      specimen_id,
      compartment_std,
      barcode_core,
      sep = "::"
    )
]

dup_keys <- cells[
  ,
  .N,
  by = technical_cell_key
][
  N > 1
]

cat("\n===== 10. TECHNICAL CELL-KEY DUPLICATES =====\n\n")

cat(
  "Duplicated specimen+compartment+barcode keys:",
  nrow(dup_keys),
  "\n"
)

if (nrow(dup_keys) > 0) {

  print(
    head(
      dup_keys,
      30
    )
  )

  fwrite(
    dup_keys,
    file.path(
      out_dir,
      "GSE160269_duplicate_technical_cell_keys.csv"
    ),
    bom = TRUE
  )
}

# ------------------------------------------------------------
# 13. Save tables
# ------------------------------------------------------------

fwrite(
  cells,
  file.path(
    out_dir,
    "GSE160269_cell_index_patient_mapping.csv.gz"
  )
)

fwrite(
  specimen_table,
  file.path(
    out_dir,
    "GSE160269_specimen_audit.csv"
  ),
  bom = TRUE
)

fwrite(
  patient_table,
  file.path(
    out_dir,
    "GSE160269_patient_audit.csv"
  ),
  bom = TRUE
)

fwrite(
  paired,
  file.path(
    out_dir,
    "GSE160269_paired_tumor_normal_patients.csv"
  ),
  bom = TRUE
)

fwrite(
  celltype_tissue,
  file.path(
    out_dir,
    "GSE160269_celltype_by_tissue.csv"
  ),
  bom = TRUE
)

fwrite(
  celltype_compartment,
  file.path(
    out_dir,
    "GSE160269_celltype_by_compartment.csv"
  ),
  bom = TRUE
)

fwrite(
  pb_table,
  file.path(
    out_dir,
    "GSE160269_pseudobulk_cell_number_audit.csv"
  ),
  bom = TRUE
)

# ------------------------------------------------------------
# 14. Final validation
# ------------------------------------------------------------

n_patients <- uniqueN(
  cells$patient_id
)

n_tumor <- uniqueN(
  cells[
    tissue_std == "Tumor",
    specimen_id
  ]
)

n_normal <- uniqueN(
  cells[
    tissue_std == "Adjacent_normal",
    specimen_id
  ]
)

cat("\n============================================================\n")
cat("FINAL VALIDATION\n")
cat("============================================================\n\n")

cat(
  "Cells             :",
  nrow(cells),
  "\n"
)

cat(
  "Patients          :",
  n_patients,
  "\n"
)

cat(
  "Tumor specimens   :",
  n_tumor,
  "\n"
)

cat(
  "Adjacent normals  :",
  n_normal,
  "\n"
)

cat(
  "Paired patients   :",
  nrow(paired),
  "\n"
)

cat("\nExpected from study description:\n")
cat("  208,659 cells\n")
cat("  60 patients\n")
cat("  60 tumor specimens\n")
cat("  4 adjacent-normal specimens\n")

cat("\n============================================================\n")

if (
  nrow(cells) == 208659 &&
  n_patients == 60 &&
  n_tumor == 60 &&
  n_normal == 4 &&
  all(cells$parse_success)
) {

  cat("STEP 13D FINISHED ✓\n")
  cat("COHORT STRUCTURE MATCHES EXPECTATION ✓\n")

} else {

  cat("STEP 13D FINISHED WITH REVIEW FLAG ⚠\n")
  cat("DO NOT MERGE YET.\n")
}

cat("============================================================\n\n")

cat("AUDIT ONLY.\n")
cat("No Seurat objects modified ✓\n")
cat("No merge ✓\n")
cat("No QC ✓\n")
cat("No normalization ✓\n")
cat("No integration ✓\n")
cat("STEP 14 NOT RUN ✓\n")

cat("\nSTOP HERE.\n")

