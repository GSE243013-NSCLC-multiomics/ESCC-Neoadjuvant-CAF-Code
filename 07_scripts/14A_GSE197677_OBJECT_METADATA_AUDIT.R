# ============================================================
# ESCC Neoadjuvant Project
# STEP 14A
# Audit GSE197677 integrated Seurat object
#
# ALSO:
# - sanity-check GSE221561 established cell counts
# - sanity-check GSE160269 merge audit
#
# NO QC FILTERING
# NO NORMALIZATION
# NO PCA / UMAP
# NO CROSS-DATASET INTEGRATION
# ============================================================

options(
  stringsAsFactors = FALSE,
  timeout = 1200
)

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
})

cat("\n============================================================\n")
cat("STEP 14A: GSE197677 OBJECT + METADATA AUDIT\n")
cat("============================================================\n\n")

out_dir <- "06_tables/GSE197677/audit"
log_dir <- "08_logs"

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ============================================================
# 1. SANITY CHECK GSE221561
# ============================================================

cat("===== 1. GSE221561 SANITY CHECK =====\n\n")

g221_audit <- paste0(
  "00_metadata/",
  "05_GSE221561_seurat_metadata_match_audit.csv"
)

if (file.exists(g221_audit)) {

  x <- fread(g221_audit)

  raw_cells <- sum(
    x$matrix_cells
  )

  author_cells <- sum(
    x$metadata_matched
  )

  cat(
    "Raw 10x matrix cells     : ",
    raw_cells,
    "\n",
    sep = ""
  )

  cat(
    "Author-retained cells    : ",
    author_cells,
    "\n",
    sep = ""
  )

  cat(
    "Raw not in author metadata: ",
    raw_cells - author_cells,
    "\n",
    sep = ""
  )

  cat(
    sprintf(
      "Author retention         : %.2f%%\n",
      100 * author_cells / raw_cells
    )
  )

  if (
    raw_cells == 90378 &&
    author_cells == 55150
  ) {
    cat(
      "✓ GSE221561 established counts are unchanged.\n"
    )
  } else {
    cat(
      "⚠ GSE221561 counts differ from the previously ",
      "validated state.\n"
    )
  }

} else {

  cat(
    "GSE221561 match-audit CSV not found.\n"
  )
}

cat(
  "\nNOTE: The value '206,991 cells' shown in the recent\n",
  "Free Code summary should NOT be accepted unless independently\n",
  "verified. Our established audit was 90,378 raw and 55,150\n",
  "author-retained cells.\n",
  sep = ""
)

# ============================================================
# 2. SANITY CHECK GSE160269
# ============================================================

cat("\n===== 2. GSE160269 SANITY CHECK =====\n\n")

g160_audit <- paste0(
  "06_tables/GSE160269/metadata/",
  "GSE160269_merge_audit.csv"
)

if (file.exists(g160_audit)) {

  x160 <- fread(
    g160_audit
  )

  print(x160)

} else {

  cat(
    "GSE160269 merge audit CSV not found.\n"
  )
}

# ============================================================
# 3. FIND GSE197677 OBJECT
# ============================================================

cat("\n===== 3. LOCATING GSE197677 OBJECT =====\n\n")

candidates <- list.files(
  "03_objects/GSE197677",
  pattern = "\\.rds$",
  full.names = TRUE,
  ignore.case = TRUE
)

if (length(candidates) == 0) {
  stop(
    "No GSE197677 RDS object found in 03_objects/GSE197677/"
  )
}

for (f in candidates) {

  cat(
    sprintf(
      "%-60s %.2f GB\n",
      basename(f),
      file.info(f)$size / 1024^3
    )
  )
}

preferred <- candidates[
  grepl(
    "integrated",
    basename(candidates),
    ignore.case = TRUE
  )
]

if (length(preferred) >= 1) {

  object_file <- preferred[1]

} else {

  object_file <- candidates[1]
}

cat(
  "\nSelected object:\n",
  object_file,
  "\n",
  sep = ""
)

# ============================================================
# 4. LOAD OBJECT
# ============================================================

cat("\n===== 4. LOADING GSE197677 OBJECT =====\n\n")
cat(
  "This is a large RDS. Loading may take several minutes.\n\n"
)

obj <- readRDS(
  object_file
)

cat("Object loaded successfully.\n")
cat("Class: ", paste(class(obj), collapse = ", "), "\n", sep = "")

# Check what slots are available
available_slots <- slotNames(obj)
cat("Available slots: ", paste(head(available_slots, 10), collapse = ", "), "...\n", sep = "")

# Try to get basic info
obj_genes <- tryCatch(nrow(obj), error = function(e) {
  cat("nrow error: ", conditionMessage(e), "\n")
  # Try alternative: use the count matrix
  tryCatch({
    counts <- GetAssayData(obj, layer = "counts")
    nrow(counts)
  }, error = function(e2) {
    tryCatch({
      counts <- GetAssayData(obj, slot = "counts")
      nrow(counts)
    }, error = function(e3) NA_integer_)
  })
})

obj_cells <- tryCatch(ncol(obj), error = function(e) {
  cat("ncol error: ", conditionMessage(e), "\n")
  tryCatch({
    counts <- GetAssayData(obj, layer = "counts")
    ncol(counts)
  }, error = function(e2) {
    tryCatch({
      counts <- GetAssayData(obj, slot = "counts")
      ncol(counts)
    }, error = function(e3) NA_integer_)
  })
})

cat("Genes: ", obj_genes, "\n", sep = "")
cat("Cells: ", obj_cells, "\n", sep = "")

obj_default_assay <- tryCatch(DefaultAssay(obj), error = function(e) "unknown")
cat("Default assay: ", obj_default_assay, "\n", sep = "")

obj_assays <- tryCatch(Assays(obj), error = function(e) character(0))
cat("Assays: ", paste(obj_assays, collapse = ", "), "\n", sep = "")

obj_reductions <- tryCatch(Reductions(obj), error = function(e) character(0))
cat("Reductions: ", paste(obj_reductions, collapse = ", "), "\n", sep = "")

# ============================================================
# 5. ASSAY / LAYER AUDIT
# ============================================================

cat("\n===== 5. ASSAY / LAYER AUDIT =====\n\n")

assay_audit <- list()

for (assay_name in Assays(obj)) {

  assay_obj <- tryCatch(
    obj[[assay_name]],
    error = function(e) NULL
  )

  if (is.null(assay_obj)) next

  layers <- tryCatch(
    Layers(assay_obj),
    error = function(e) {
      character(0)
    }
  )

  cat(
    "\nAssay: ",
    assay_name,
    "\n",
    sep = ""
  )

  cat(
    "Layers: ",
    if (length(layers) == 0) {
      "<not available>"
    } else {
      paste(layers, collapse = ", ")
    },
    "\n",
    sep = ""
  )

  has_counts <- any(
    grepl(
      "^counts($|\\.)",
      layers,
      ignore.case = TRUE
    )
  )

  has_data <- any(
    grepl(
      "^data($|\\.)",
      layers,
      ignore.case = TRUE
    )
  )

  assay_audit[[assay_name]] <- data.table(
    assay = assay_name,
    layers = paste(
      layers,
      collapse = "; "
    ),
    has_counts_layer = has_counts,
    has_data_layer = has_data
  )
}

assay_audit <- rbindlist(
  assay_audit,
  fill = TRUE
)

fwrite(
  assay_audit,
  file.path(
    out_dir,
    "GSE197677_assay_layer_audit.csv"
  ),
  bom = TRUE
)

cat("\nAssay summary:\n")
print(assay_audit)

# ============================================================
# 6. METADATA FIELD INVENTORY
# ============================================================

cat("\n===== 6. METADATA FIELD INVENTORY =====\n\n")

md <- tryCatch(
  as.data.table(
    obj[[]],
    keep.rownames = "cell_id"
  ),
  error = function(e) {
    cat("Warning: obj[[]] failed, trying alternative method\n")
    # Try to get metadata directly from the slot
    tryCatch({
      md_raw <- slot(obj, "meta.data")
      md_dt <- as.data.table(md_raw, keep.rownames = "cell_id")
      md_dt
    }, error = function(e2) {
      stop("Cannot extract metadata from object: ", conditionMessage(e2))
    })
  }
)

field_table <- data.table(
  column_number = seq_along(names(md)),
  column_name = names(md),
  class = vapply(
    md,
    function(x) {
      paste(
        class(x),
        collapse = ";"
      )
    },
    character(1)
  ),
  unique_values = vapply(
    md,
    uniqueN,
    integer(1),
    na.rm = TRUE
  ),
  missing_values = vapply(
    md,
    function(x) {
      sum(is.na(x))
    },
    integer(1)
  )
)

print(field_table)

fwrite(
  field_table,
  file.path(
    out_dir,
    "GSE197677_metadata_field_inventory.csv"
  ),
  bom = TRUE
)

# ============================================================
# 7. DETECT IMPORTANT COLUMNS
# ============================================================

find_fields <- function(pattern) {

  grep(
    pattern,
    names(md),
    ignore.case = TRUE,
    value = TRUE
  )
}

candidate_fields <- list(

  sample = find_fields(
    "sample|orig\\.ident|library|subject"
  ),

  patient = find_fields(
    "patient|subject|donor|case"
  ),

  tissue = find_fields(
    "tissue|tumou?r|normal|source"
  ),

  treatment = find_fields(
    "NAC|treat|therapy|neoadjuvant"
  ),

  response = find_fields(
    "response|pCR|TRG|regression"
  ),

  celltype = find_fields(
    "cell.?type|celltype|annotation|cluster|major|minor"
  ),

  qc = find_fields(
    "nCount|nFeature|percent.*mt|mito|UMI"
  )
)

cat("\n===== 7. IMPORTANT METADATA CANDIDATES =====\n")

for (nm in names(candidate_fields)) {

  cat(
    "\n",
    toupper(nm),
    ":\n",
    sep = ""
  )

  hits <- candidate_fields[[nm]]

  if (length(hits) == 0) {

    cat("  <none detected>\n")

  } else {

    for (x in hits) {
      cat(
        "  - ",
        x,
        "\n",
        sep = ""
      )
    }
  }
}

# ============================================================
# 8. PRINT CATEGORICAL VALUES
# ============================================================

cat("\n===== 8. CATEGORICAL VALUE AUDIT =====\n")

interesting <- unique(
  unlist(
    candidate_fields[
      c(
        "sample",
        "patient",
        "tissue",
        "treatment",
        "response",
        "celltype"
      )
    ]
  )
)

categorical_summary <- list()

for (col in interesting) {

  n_unique <- uniqueN(
    md[[col]],
    na.rm = FALSE
  )

  cat(
    "\n------------------------------------------------------------\n"
  )

  cat(
    col,
    "\n"
  )

  cat(
    "Unique values: ",
    n_unique,
    "\n",
    sep = ""
  )

  if (n_unique <= 60) {

    tab <- sort(
      table(
        md[[col]],
        useNA = "ifany"
      ),
      decreasing = TRUE
    )

    print(tab)

    tmp <- as.data.table(
      tab
    )

    setnames(
      tmp,
      c(
        "value",
        "n_cells"
      )
    )

    tmp[
      ,
      field := col
    ]

    categorical_summary[[col]] <- tmp

  } else {

    cat(
      "Too many values; top 20 shown:\n"
    )

    print(
      head(
        sort(
          table(
            md[[col]],
            useNA = "ifany"
          ),
          decreasing = TRUE
        ),
        20
      )
    )
  }
}

if (length(categorical_summary) > 0) {

  fwrite(
    rbindlist(
      categorical_summary,
      fill = TRUE
    ),
    file.path(
      out_dir,
      "GSE197677_categorical_metadata_values.csv"
    ),
    bom = TRUE
  )
}

# ============================================================
# 9. QC FIELD SUMMARY — AUDIT ONLY
# ============================================================

cat("\n===== 9. QC FIELD SUMMARY =====\n\n")

qc_fields <- unique(
  candidate_fields$qc
)

if (length(qc_fields) == 0) {

  cat(
    "No obvious QC fields detected.\n"
  )

} else {

  for (col in qc_fields) {

    x <- suppressWarnings(
      as.numeric(
        md[[col]]
      )
    )

    keep <- is.finite(
      x
    )

    cat(
      "\n",
      col,
      "\n",
      sep = ""
    )

    if (sum(keep) == 0) {

      cat(
        "  Not numeric / no finite values\n"
      )

    } else {

      q <- quantile(
        x[keep],
        probs = c(
          0,
          0.01,
          0.05,
          0.25,
          0.5,
          0.75,
          0.95,
          0.99,
          1
        ),
        na.rm = TRUE
      )

      print(q)
    }
  }
}

# ============================================================
# 10. SEARCH FOR EXPECTED GSE197677 STRUCTURE
# ============================================================

cat("\n===== 10. EXPECTED COHORT STRUCTURE CHECK =====\n\n")

cat(
  "Study expectation from project plan:\n"
)

cat(
  "  18 ESCC tumor samples\n"
)

cat(
  "  12 normal esophageal samples\n"
)

cat(
  "  NAC(+) and NAC(-) represented\n\n"
)

# We intentionally do NOT guess which field encodes these.
# Only report obvious values detected from metadata.

for (col in names(md)) {

  vals <- unique(
    trimws(
      as.character(
        md[[col]]
      )
    )
  )

  vals_low <- tolower(
    vals
  )

  interesting_hit <- any(
    grepl(
      "nac|tumou?r|normal",
      vals_low
    ),
    na.rm = TRUE
  )

  if (interesting_hit) {

    cat(
      "Potential cohort field: ",
      col,
      "\n",
      sep = ""
    )

    tab <- sort(
      table(
        md[[col]],
        useNA = "ifany"
      ),
      decreasing = TRUE
    )

    print(tab)

    cat("\n")
  }
}

# ============================================================
# 11. SAVE OBJECT SUMMARY
# ============================================================

object_summary <- data.table(
  object_file = object_file,
  object_size_GB = round(
    file.info(object_file)$size / 1024^3,
    3
  ),
  genes = obj_genes,
  cells = obj_cells,
  default_assay = obj_default_assay,
  assays = paste(
    obj_assays,
    collapse = "; "
  ),
  reductions = paste(
    obj_reductions,
    collapse = "; "
  ),
  metadata_fields = ncol(md)
)

fwrite(
  object_summary,
  file.path(
    out_dir,
    "GSE197677_object_summary.csv"
  ),
  bom = TRUE
)

# ============================================================
# 12. SESSION INFO
# ============================================================

sink(
  file.path(
    log_dir,
    "14A_GSE197677_OBJECT_METADATA_AUDIT_sessionInfo.txt"
  )
)

print(
  sessionInfo()
)

sink()

# ============================================================
# FINAL
# ============================================================

cat("\n============================================================\n")
cat("STEP 14A FINISHED ✓\n")
cat("============================================================\n\n")

cat(
  "GSE197677 object was inspected, not modified.\n\n"
)

cat("No cells filtered ✓\n")
cat("No normalization ✓\n")
cat("No PCA / UMAP ✓\n")
cat("No cross-dataset integration ✓\n")

cat("\nSTOP HERE.\n")

