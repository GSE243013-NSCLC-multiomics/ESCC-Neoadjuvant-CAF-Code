
options(
  stringsAsFactors = FALSE
)

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
})

proc_dir <- "02_processed/GSE160269/10x_parts"
obj_dir  <- "03_objects/GSE160269/parts"

dir.create(
  obj_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

manifest <- fread(
  "00_metadata/06_GSE160269_UMI_matrix_manifest.tsv"
)

audit <- list()


for (i in seq_len(nrow(manifest))) {

  label <- manifest$label[i]

  cat(
    "\n============================================================\n"
  )

  cat(
    "Building:",
    label,
    "\n"
  )

  cat(
    "============================================================\n"
  )

  part_dir <- file.path(
    proc_dir,
    label
  )

  out_file <- file.path(
    obj_dir,
    paste0(
      "GSE160269_",
      label,
      ".rds"
    )
  )


  # ----------------------------------------------------------
  # Resume support
  # ----------------------------------------------------------

  if (file.exists(out_file)) {

    cat(
      "✓ Object already exists — skipping.\n"
    )

    info <- file.info(
      out_file
    )

    audit[[label]] <- data.table(
      label = label,
      object_file = out_file,
      status = "existing",
      object_size_MB =
        round(
          info$size / 1024^2,
          2
        )
    )

    next
  }


  # ----------------------------------------------------------
  # Read sparse matrix
  # ----------------------------------------------------------

  counts <- ReadMtx(
    mtx = file.path(
      part_dir,
      "matrix.mtx.gz"
    ),

    cells = file.path(
      part_dir,
      "barcodes.tsv.gz"
    ),

    features = file.path(
      part_dir,
      "features.tsv.gz"
    ),

    cell.column = 1,
    feature.column = 1,
    unique.features = TRUE
  )


  original_barcode <- colnames(
    counts
  )


  # ----------------------------------------------------------
  # IMPORTANT:
  # Make globally unique cell names BEFORE any future merge.
  #
  # This prevents the exact problem that broke the old script.
  # ----------------------------------------------------------

  globally_unique_name <- paste(
    label,
    original_barcode,
    sep = "__"
  )

  colnames(
    counts
  ) <- globally_unique_name


  # ----------------------------------------------------------
  # Build Seurat object
  # ----------------------------------------------------------

  obj <- CreateSeuratObject(
    counts = counts,
    project = "GSE160269",
    min.cells = 0,
    min.features = 0
  )


  obj$dataset <- "GSE160269"

  obj$source_matrix <- label

  obj$source_celltype <- label

  obj$original_barcode <-
    original_barcode


  cat(
    "Genes:",
    nrow(obj),
    "\n"
  )

  cat(
    "Cells:",
    ncol(obj),
    "\n"
  )


  # ----------------------------------------------------------
  # Save immediately
  # ----------------------------------------------------------

  saveRDS(
    obj,
    out_file,
    compress = FALSE
  )


  info <- file.info(
    out_file
  )


  audit[[label]] <- data.table(
    label = label,
    genes = nrow(obj),
    cells = ncol(obj),
    object_file = out_file,
    status = "built",
    object_size_MB =
      round(
        info$size / 1024^2,
        2
      )
  )


  cat(
    "✓ Saved:",
    out_file,
    "\n"
  )


  # ----------------------------------------------------------
  # Explicitly release RAM before next matrix
  # ----------------------------------------------------------

  rm(
    obj,
    counts
  )

  invisible(
    gc()
  )
}


# ============================================================
# Audit
# ============================================================

audit_dt <- rbindlist(
  audit,
  fill = TRUE
)

fwrite(
  audit_dt,
  "00_metadata/06_GSE160269_part_object_audit.csv",
  bom = TRUE
)

cat(
  "\n============================================================\n"
)

cat(
  "GSE160269 PART OBJECT AUDIT\n"
)

cat(
  "============================================================\n\n"
)

print(
  audit_dt
)


cat(
  "\nObjects present:\n"
)

objects <- list.files(
  obj_dir,
  pattern = "\\.rds$",
  full.names = TRUE
)

for (f in objects) {

  info <- file.info(f)

  cat(
    sprintf(
      "  %-50s %8.2f MB\n",
      basename(f),
      info$size / 1024^2
    )
  )
}


cat(
  "\n============================================================\n"
)

cat(
  "STEP 13B FINISHED ✓\n"
)

cat(
  "============================================================\n\n"
)

cat(
  "IMPORTANT:\n"
)

cat(
  "  GSE160269 parts were built separately.\n"
)

cat(
  "  They were NOT merged.\n"
)

cat(
  "  Metadata matching was NOT attempted yet.\n"
)

cat(
  "  QC was NOT run.\n"
)

cat(
  "  Normalization was NOT run.\n"
)

cat(
  "  Integration was NOT run.\n"
)

cat(
  "  STEP 14 was NOT run.\n\n"
)

cat(
  "STOP HERE and review the audit before continuing.\n"
)

