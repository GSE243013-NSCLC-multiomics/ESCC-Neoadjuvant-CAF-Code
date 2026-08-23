# ============================================================
# ESCC Neoadjuvant Project
# STEP 13F
# Merge 8 annotated GSE160269 parts
#
# ONLY:
#   merge raw-count objects
#   validate metadata
#   save merged object
#
# NO QC
# NO NORMALIZATION
# NO PCA / UMAP
# NO CROSS-DATASET INTEGRATION
# ============================================================

options(
  stringsAsFactors = FALSE,
  future.globals.maxSize = 50 * 1024^3
)

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
})

cat("\n============================================================\n")
cat("STEP 13F: MERGE GSE160269 ANNOTATED PARTS\n")
cat("============================================================\n\n")

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------

input_dir <- "03_objects/GSE160269/annotated_parts"

output_dir <- "03_objects/GSE160269"

output_file <- file.path(
  output_dir,
  "GSE160269_annotated_merged_raw.rds"
)

checkpoint_dir <- file.path(
  output_dir,
  "merge_checkpoints"
)

table_dir <- "06_tables/GSE160269/metadata"

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  checkpoint_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# 1. Find the 8 annotated parts
# ------------------------------------------------------------

files <- sort(
  list.files(
    input_dir,
    pattern = "\\.rds$",
    full.names = TRUE
  )
)

cat("===== 1. INPUT OBJECTS =====\n\n")

if (length(files) != 8) {
  stop(
    "Expected 8 annotated objects, found ",
    length(files)
  )
}

for (f in files) {

  info <- file.info(f)

  cat(
    sprintf(
      "%-55s %8.2f MB\n",
      basename(f),
      info$size / 1024^2
    )
  )
}

# ------------------------------------------------------------
# 2. Fast pre-merge metadata validation
# ------------------------------------------------------------

cat("\n===== 2. PRE-MERGE VALIDATION =====\n\n")

audit <- list()
all_cell_ids <- character(0)
all_technical_keys <- character(0)

for (f in files) {

  obj <- readRDS(f)

  required <- c(
    "patient_id",
    "specimen_id",
    "tissue",
    "compartment_code",
    "original_barcode",
    "technical_cell_key",
    "cell_type_author_source",
    "dataset",
    "treatment_group"
  )

  missing <- setdiff(
    required,
    colnames(obj[[]])
  )

  if (length(missing) > 0) {
    stop(
      basename(f),
      " missing metadata fields: ",
      paste(missing, collapse = ", ")
    )
  }

  if (any(is.na(obj$patient_id))) {
    stop(
      "Missing patient_id in ",
      basename(f)
    )
  }

  ids <- colnames(obj)
  keys <- as.character(
    obj$technical_cell_key
  )

  audit[[basename(f)]] <- data.table(
    object = basename(f),
    genes = nrow(obj),
    cells = ncol(obj),
    patients = uniqueN(obj$patient_id),
    specimens = uniqueN(obj$specimen_id),
    tumor_cells = sum(
      obj$tissue == "Tumor"
    ),
    normal_cells = sum(
      obj$tissue == "Adjacent_normal"
    )
  )

  all_cell_ids <- c(
    all_cell_ids,
    ids
  )

  all_technical_keys <- c(
    all_technical_keys,
    keys
  )

  rm(obj)
  gc(verbose = FALSE)
}

audit_dt <- rbindlist(
  audit,
  fill = TRUE
)

print(audit_dt)

total_pre <- sum(
  audit_dt$cells
)

cat(
  "\nTotal cells before merge:",
  total_pre,
  "\n"
)

if (total_pre != 208659) {
  stop(
    "Expected 208659 cells before merge."
  )
}

if (anyDuplicated(all_cell_ids)) {
  stop(
    "Duplicated Seurat cell IDs detected before merge."
  )
}

if (anyDuplicated(all_technical_keys)) {
  stop(
    "Duplicated technical_cell_key detected before merge."
  )
}

cat("✓ Cell IDs globally unique\n")
cat("✓ Technical cell keys globally unique\n")

rm(
  all_cell_ids,
  all_technical_keys
)

gc(verbose = FALSE)

# ------------------------------------------------------------
# 3. Resume support
# ------------------------------------------------------------

if (file.exists(output_file)) {

  cat(
    "\nFinal merged object already exists:\n",
    output_file,
    "\n",
    sep = ""
  )

  cat(
    "Will validate existing object instead of rebuilding.\n"
  )

  merged <- readRDS(
    output_file
  )

} else {

  # ----------------------------------------------------------
  # 4. Sequential merge
  #
  # Only one new part is loaded at a time.
  # Checkpoint after every successful merge.
  # ----------------------------------------------------------

  cat("\n===== 3. SEQUENTIAL MERGE =====\n\n")

  # First object
  merged <- readRDS(
    files[1]
  )

  cat(
    "START:",
    basename(files[1]),
    "| cells:",
    ncol(merged),
    "\n"
  )

  if (length(files) > 1) {

    for (i in 2:length(files)) {

      cat(
        "\n------------------------------------------------------------\n"
      )

      cat(
        "Adding ",
        i,
        "/",
        length(files),
        ": ",
        basename(files[i]),
        "\n",
        sep = ""
      )

      next_obj <- readRDS(
        files[i]
      )

      old_n <- ncol(
        merged
      )

      add_n <- ncol(
        next_obj
      )

      # All cell IDs are already globally unique.
      merged_new <- merge(
        x = merged,
        y = next_obj,
        merge.data = FALSE
      )

      expected_n <- old_n + add_n

      if (ncol(merged_new) != expected_n) {

        stop(
          "Cell-count mismatch after merging ",
          basename(files[i]),
          ": expected ",
          expected_n,
          ", found ",
          ncol(merged_new)
        )
      }

      rm(
        merged,
        next_obj
      )

      merged <- merged_new

      rm(
        merged_new
      )

      gc(verbose = FALSE)

      cat(
        "Merged cells:",
        ncol(merged),
        "\n"
      )

      # ----------------------------------------------
      # Save checkpoint.
      # Keep only newest checkpoint to save disk.
      # ----------------------------------------------

      checkpoint <- file.path(
        checkpoint_dir,
        sprintf(
          "GSE160269_merge_%02d_of_%02d.rds",
          i,
          length(files)
        )
      )

      saveRDS(
        merged,
        checkpoint,
        compress = FALSE
      )

      cat(
        "✓ Checkpoint saved: ",
        checkpoint,
        "\n",
        sep = ""
      )

      old_checkpoints <- list.files(
        checkpoint_dir,
        pattern = "\\.rds$",
        full.names = TRUE
      )

      old_checkpoints <- setdiff(
        old_checkpoints,
        checkpoint
      )

      if (length(old_checkpoints) > 0) {

        file.remove(
          old_checkpoints
        )
      }
    }
  }

  # ----------------------------------------------------------
  # 5. Save final object
  # ----------------------------------------------------------

  cat("\n===== 4. SAVING FINAL MERGED OBJECT =====\n\n")

  saveRDS(
    merged,
    output_file,
    compress = FALSE
  )

  cat(
    "✓ Saved:\n",
    output_file,
    "\n",
    sep = ""
  )
}

# ------------------------------------------------------------
# 6. Final validation
# ------------------------------------------------------------

cat("\n===== 5. FINAL VALIDATION =====\n\n")

final_cells <- ncol(
  merged
)

final_genes <- nrow(
  merged
)

final_patients <- uniqueN(
  merged$patient_id
)

final_specimens <- uniqueN(
  merged$specimen_id
)

final_tumor_specimens <- uniqueN(
  merged$specimen_id[
    merged$tissue == "Tumor"
  ]
)

final_normal_specimens <- uniqueN(
  merged$specimen_id[
    merged$tissue == "Adjacent_normal"
  ]
)

final_celltypes <- sort(
  unique(
    as.character(
      merged$cell_type_author_source
    )
  )
)

duplicate_cell_ids <- anyDuplicated(
  colnames(merged)
)

duplicate_keys <- anyDuplicated(
  as.character(
    merged$technical_cell_key
  )
)

cat(
  "Genes              : ",
  final_genes,
  "\n",
  sep = ""
)

cat(
  "Cells              : ",
  final_cells,
  "\n",
  sep = ""
)

cat(
  "Patients           : ",
  final_patients,
  "\n",
  sep = ""
)

cat(
  "Specimens          : ",
  final_specimens,
  "\n",
  sep = ""
)

cat(
  "Tumor specimens    : ",
  final_tumor_specimens,
  "\n",
  sep = ""
)

cat(
  "Adjacent normals   : ",
  final_normal_specimens,
  "\n",
  sep = ""
)

cat(
  "Cell types         : ",
  length(final_celltypes),
  "\n",
  sep = ""
)

cat(
  "Cell-type labels   : ",
  paste(
    final_celltypes,
    collapse = ", "
  ),
  "\n",
  sep = ""
)

cat(
  "Duplicate cell IDs : ",
  duplicate_cell_ids,
  "\n",
  sep = ""
)

cat(
  "Duplicate keys     : ",
  duplicate_keys,
  "\n",
  sep = ""
)

# ------------------------------------------------------------
# 7. Cell type counts
# ------------------------------------------------------------

celltype_counts <- as.data.table(
  table(
    merged$cell_type_author_source
  )
)

setnames(
  celltype_counts,
  c(
    "cell_type",
    "n_cells"
  )
)

setorder(
  celltype_counts,
  -n_cells
)

cat("\n===== 6. CELL TYPE COUNTS =====\n\n")

print(
  celltype_counts
)

fwrite(
  celltype_counts,
  file.path(
    table_dir,
    "GSE160269_merged_celltype_counts.csv"
  ),
  bom = TRUE
)

# ------------------------------------------------------------
# 8. Patient / tissue summary
# ------------------------------------------------------------

metadata_dt <- as.data.table(
  merged[[]],
  keep.rownames = "cell_id"
)

patient_tissue <- metadata_dt[
  ,
  .(
    n_cells = .N
  ),
  by = .(
    patient_id,
    specimen_id,
    tissue,
    cell_type_author_source
  )
]

fwrite(
  patient_tissue,
  file.path(
    table_dir,
    "GSE160269_merged_patient_tissue_celltype_counts.csv"
  ),
  bom = TRUE
)

# ------------------------------------------------------------
# 9. Object size
# ------------------------------------------------------------

size_gb <- file.info(
  output_file
)$size / 1024^3

cat(
  "\nMerged RDS size: ",
  round(size_gb, 2),
  " GB\n",
  sep = ""
)

# ------------------------------------------------------------
# 10. Save merge audit
# ------------------------------------------------------------

merge_audit <- data.table(
  dataset = "GSE160269",

  genes = final_genes,

  cells = final_cells,

  patients = final_patients,

  specimens = final_specimens,

  tumor_specimens =
    final_tumor_specimens,

  adjacent_normal_specimens =
    final_normal_specimens,

  cell_types =
    length(final_celltypes),

  duplicate_cell_ids =
    duplicate_cell_ids,

  duplicate_technical_keys =
    duplicate_keys,

  object_size_GB =
    round(size_gb, 3)
)

fwrite(
  merge_audit,
  file.path(
    table_dir,
    "GSE160269_merge_audit.csv"
  ),
  bom = TRUE
)

# ------------------------------------------------------------
# 11. Final decision
# ------------------------------------------------------------

cat("\n============================================================\n")

if (
  final_cells == 208659 &&
  final_patients == 60 &&
  final_specimens == 64 &&
  final_tumor_specimens == 60 &&
  final_normal_specimens == 4 &&
  length(final_celltypes) == 8 &&
  duplicate_cell_ids == 0 &&
  duplicate_keys == 0
) {

  cat("STEP 13F FINISHED ✓\n")
  cat("GSE160269 MERGE VALIDATED ✓\n")

} else {

  cat("STEP 13F FINISHED WITH REVIEW FLAG ⚠\n")
  cat("DO NOT CONTINUE TO CROSS-DATASET ANALYSIS.\n")
}

cat("============================================================\n\n")

cat(
  "Merged object:\n",
  output_file,
  "\n\n",
  sep = ""
)

cat("IMPORTANT:\n")
cat("No cells filtered ✓\n")
cat("No normalization ✓\n")
cat("No PCA / UMAP ✓\n")
cat("No cross-dataset integration ✓\n")
cat("STEP 14 NOT RUN ✓\n")

cat("\nSTOP HERE.\n")

