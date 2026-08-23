# ============================================================
# ESCC Neoadjuvant Project
# Step 18B
# Build 11 independent native raw-count Seurat objects
# Feature identity = Ensembl ID (features.tsv.gz column 1)
#
# NO normalization
# NO filtering
# NO merge
# NO integration
# NO author cell-type annotation
# ============================================================

options(
  stringsAsFactors = FALSE,
  timeout = 1200
)

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(Matrix)
})

cat("\n============================================================\n")
cat("STEP 18B: GSE221561 NATIVE RAW SEURAT OBJECTS\n")
cat("============================================================\n\n")

# ------------------------------------------------------------
# 1. Paths
# ------------------------------------------------------------

PROJECT <- normalizePath(".", mustWork = TRUE)

RAW_DIR <- file.path(
  PROJECT,
  "02_processed/GSE221561/01_10x_raw"
)

OBJECT_DIR <- file.path(
  PROJECT,
  "03_objects/GSE221561/native_raw"
)

META_DIR <- file.path(
  PROJECT,
  "00_metadata"
)

LOG_DIR <- file.path(
  PROJECT,
  "08_logs"
)

dir.create(
  OBJECT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  LOG_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

required_files <- c(
  "00_metadata/04_GSE221561_crosswalk.csv",
  "00_metadata/05_GSE221561_extracted_file_integrity.csv",
  "00_metadata/05_GSE221561_10x_structure_audit.csv",
  "00_metadata/05_GSE221561_cell_barcode_linkage_audit.csv",
  "00_metadata/06_GSE221561_feature_space_audit.csv",
  "00_metadata/06_GSE221561_common_gene_features.csv",
  "00_metadata/06_GSE221561_feature_name_conflicts.csv",
  "00_metadata/06_GSE221561_feature_space_policy.csv"
)

missing_required <- required_files[
  !file.exists(required_files)
]

if (length(missing_required) > 0) {
  stop(
    "Missing required files:\n",
    paste(missing_required, collapse = "\n")
  )
}

if (!dir.exists(RAW_DIR)) {
  stop(
    "Cannot find extracted 10x directory:\n",
    RAW_DIR
  )
}

cat("Preflight files: OK\n")

# ------------------------------------------------------------
# 2. Frozen common-feature manifest
# ------------------------------------------------------------

common_dt <- fread(
  "00_metadata/06_GSE221561_common_gene_features.csv",
  encoding = "UTF-8"
)

ensembl_candidates <- grep(
  "ensembl",
  names(common_dt),
  ignore.case = TRUE,
  value = TRUE
)

if (length(ensembl_candidates) == 0) {
  stop(
    "Cannot identify Ensembl-ID column in common-feature manifest."
  )
}

common_ensembl_col <- ensembl_candidates[1]

common_features <- as.character(
  common_dt[[common_ensembl_col]]
)

common_features <- common_features[
  !is.na(common_features) &
    nzchar(common_features)
]

if (anyDuplicated(common_features)) {
  stop(
    "Common-feature manifest contains duplicated Ensembl IDs."
  )
}

if (length(common_features) != 36476) {
  stop(
    "COMMON_FEATURE_MANIFEST_ROWS mismatch: expected 36476, found ",
    length(common_features)
  )
}

cat(
  "COMMON_FEATURE_MANIFEST_ROWS =",
  length(common_features),
  "\n"
)

cat(
  "FEATURE_IDENTITY_KEY = Ensembl_ID\n"
)

# ------------------------------------------------------------
# 3. Discover the 33 raw 10x files
# ------------------------------------------------------------

all_raw_files <- list.files(
  RAW_DIR,
  recursive = TRUE,
  full.names = TRUE
)

all_raw_files <- all_raw_files[
  !file.info(all_raw_files)$isdir
]

matrix_files <- all_raw_files[
  grepl(
    "_matrix\\.mtx\\.gz$",
    basename(all_raw_files),
    ignore.case = TRUE
  )
]

barcode_files <- all_raw_files[
  grepl(
    "_barcodes\\.tsv\\.gz$",
    basename(all_raw_files),
    ignore.case = TRUE
  )
]

feature_files <- all_raw_files[
  grepl(
    "_features\\.tsv\\.gz$",
    basename(all_raw_files),
    ignore.case = TRUE
  )
]

cat("\nRaw file discovery:\n")
cat("matrix files   :", length(matrix_files), "\n")
cat("barcode files  :", length(barcode_files), "\n")
cat("feature files  :", length(feature_files), "\n")

if (
  length(matrix_files) != 11 ||
  length(barcode_files) != 11 ||
  length(feature_files) != 11
) {
  stop(
    "Expected exactly 11 matrix + 11 barcode + 11 feature files."
  )
}

# ------------------------------------------------------------
# 4. Parse sample ID and GSM from filenames
# ------------------------------------------------------------

strip_file_suffix <- function(x) {

  x <- basename(x)

  x <- sub(
    "_matrix\\.mtx\\.gz$",
    "",
    x,
    ignore.case = TRUE
  )

  x <- sub(
    "_barcodes\\.tsv\\.gz$",
    "",
    x,
    ignore.case = TRUE
  )

  x <- sub(
    "_features\\.tsv\\.gz$",
    "",
    x,
    ignore.case = TRUE
  )

  x
}


parse_gsm <- function(x) {

  base <- strip_file_suffix(x)

  m <- regexpr(
    "GSM[0-9]+",
    base,
    ignore.case = TRUE
  )

  if (m[1] == -1) {
    return(NA_character_)
  }

  regmatches(base, m)
}


parse_sample_id <- function(x) {

  base <- strip_file_suffix(x)

  # Typical form:
  # GSMxxxxxxx_SAMPLEID

  if (grepl("^GSM[0-9]+_", base, ignore.case = TRUE)) {

    return(
      sub(
        "^GSM[0-9]+_",
        "",
        base,
        ignore.case = TRUE
      )
    )
  }

  # Fallback: final token
  sub("^.*_", "", base)
}


matrix_map <- data.table(
  sample_id = vapply(
    matrix_files,
    parse_sample_id,
    character(1)
  ),
  GSM = vapply(
    matrix_files,
    parse_gsm,
    character(1)
  ),
  matrix_path = matrix_files
)

barcode_map <- data.table(
  sample_id = vapply(
    barcode_files,
    parse_sample_id,
    character(1)
  ),
  barcode_path = barcode_files
)

feature_map <- data.table(
  sample_id = vapply(
    feature_files,
    parse_sample_id,
    character(1)
  ),
  features_path = feature_files
)

triplets <- merge(
  matrix_map,
  barcode_map,
  by = "sample_id",
  all = TRUE
)

triplets <- merge(
  triplets,
  feature_map,
  by = "sample_id",
  all = TRUE
)

setorder(
  triplets,
  sample_id
)

if (nrow(triplets) != 11) {
  stop(
    "Expected 11 sample triplets, found ",
    nrow(triplets)
  )
}

if (
  anyNA(triplets$matrix_path) ||
  anyNA(triplets$barcode_path) ||
  anyNA(triplets$features_path)
) {
  print(triplets)
  stop(
    "At least one sample does not have a complete 10x triplet."
  )
}

if (anyDuplicated(triplets$sample_id)) {
  stop(
    "Duplicate sample IDs detected."
  )
}

cat("\n===== EXACT 10x TRIPLETS =====\n\n")
print(
  triplets[
    ,
    .(
      sample_id,
      GSM,
      matrix_file = basename(matrix_path),
      barcode_file = basename(barcode_path),
      feature_file = basename(features_path)
    )
  ]
)

# ------------------------------------------------------------
# 5. Read verified crosswalk
# ------------------------------------------------------------

crosswalk <- fread(
  "00_metadata/04_GSE221561_crosswalk.csv",
  encoding = "UTF-8"
)

cat("\nCrosswalk rows   :", nrow(crosswalk), "\n")
cat("Crosswalk columns:", ncol(crosswalk), "\n")

sample_col_candidates <- grep(
  "^sample_id$|sample.*id|library.*id|sample",
  names(crosswalk),
  ignore.case = TRUE,
  value = TRUE
)

gsm_col_candidates <- grep(
  "^gsm$|geo.*accession|gsm",
  names(crosswalk),
  ignore.case = TRUE,
  value = TRUE
)

sample_crosswalk_col <- NA_character_

if (length(sample_col_candidates) > 0) {

  scores <- vapply(
    sample_col_candidates,
    function(nm) {
      sum(
        triplets$sample_id %in%
          as.character(crosswalk[[nm]])
      )
    },
    numeric(1)
  )

  if (max(scores) > 0) {
    sample_crosswalk_col <-
      sample_col_candidates[which.max(scores)]
  }
}

gsm_crosswalk_col <- NA_character_

if (length(gsm_col_candidates) > 0) {

  scores <- vapply(
    gsm_col_candidates,
    function(nm) {
      sum(
        triplets$GSM %in%
          as.character(crosswalk[[nm]])
      )
    },
    numeric(1)
  )

  if (max(scores) > 0) {
    gsm_crosswalk_col <-
      gsm_col_candidates[which.max(scores)]
  }
}

cat(
  "Crosswalk sample column:",
  sample_crosswalk_col,
  "\n"
)

cat(
  "Crosswalk GSM column   :",
  gsm_crosswalk_col,
  "\n"
)

if (
  is.na(sample_crosswalk_col) &&
  is.na(gsm_crosswalk_col)
) {
  stop(
    "Cannot link the 11 raw samples to the verified crosswalk."
  )
}

# ------------------------------------------------------------
# 6. Helper: find supported provenance columns
# ------------------------------------------------------------

find_col <- function(pattern) {

  hits <- grep(
    pattern,
    names(crosswalk),
    ignore.case = TRUE,
    value = TRUE
  )

  if (length(hits) == 0) {
    return(NA_character_)
  }

  hits[1]
}


prov_columns <- c(
  cohort = find_col("^cohort$"),
  GSE = find_col("^gse$"),
  patient_id = find_col(
    "^patient_id$|patient.*id|subject.*id"
  ),
  tissue = find_col(
    "^tissue$|tissue.*type|source.*tissue"
  ),
  treatment = find_col(
    "^treatment$|regimen|therapy|neoadjuvant"
  ),
  timepoint = find_col(
    "^timepoint$|time.*point|pre.*post"
  ),
  response = find_col(
    "^response$|responder|pathologic.*response|trg"
  ),
  GEO_title = find_col(
    "^geo_title$|^title$|sample.*title"
  ),
  library_name = find_col(
    "^library_name$|library.*name"
  )
)

prov_columns <- prov_columns[
  !is.na(prov_columns)
]

cat("\nSupported provenance columns detected:\n")

if (length(prov_columns) == 0) {
  cat("  <none beyond sample/GSM mapping>\n")
} else {
  for (nm in names(prov_columns)) {
    cat(
      "  ",
      nm,
      " <- ",
      prov_columns[[nm]],
      "\n",
      sep = ""
    )
  }
}

# ------------------------------------------------------------
# 7. Helper: SHA256
# ------------------------------------------------------------

sha256_file <- function(path) {

  out <- system2(
    "shasum",
    c(
      "-a",
      "256",
      shQuote(path)
    ),
    stdout = TRUE,
    stderr = TRUE
  )

  if (length(out) == 0) {
    return(NA_character_)
  }

  strsplit(
    trimws(out[1]),
    "\\s+"
  )[[1]][1]
}

# ------------------------------------------------------------
# 8. Build objects one sample at a time
# ------------------------------------------------------------

audit_rows <- list()
objects_ok <- character()

all_global_cell_ids <- character()

total_raw_cells <- 0L

cat("\n============================================================\n")
cat("BUILDING NATIVE OBJECTS\n")
cat("============================================================\n")

for (i in seq_len(nrow(triplets))) {

  sample_id <- triplets$sample_id[i]
  GSM <- triplets$GSM[i]

  matrix_path <- triplets$matrix_path[i]
  barcode_path <- triplets$barcode_path[i]
  features_path <- triplets$features_path[i]

  cat("\n------------------------------------------------------------\n")
  cat("SAMPLE:", sample_id, "\n")
  cat("GSM   :", GSM, "\n")
  cat("------------------------------------------------------------\n")

  # ----------------------------------------------------------
  # A. Read feature table independently
  # ----------------------------------------------------------

  # Read compressed feature table directly through gzfile().
  # Do NOT use fread() here because fread may decompress .gz into
  # R's session tempdir, which can fail in the Free Code/OpenCode
  # environment. Streaming via gzfile avoids that temp-file dependency.
  feature_con <- gzfile(
    features_path,
    open = "rt"
  )

  feature_table <- tryCatch(
    {
      as.data.table(
        read.delim(
          feature_con,
          header = FALSE,
          sep = "\t",
          quote = "",
          comment.char = "",
          stringsAsFactors = FALSE,
          check.names = FALSE
        )
      )
    },
    finally = {
      close(feature_con)
    }
  )

  if (ncol(feature_table) < 1) {
    stop(
      "No feature column found for sample ",
      sample_id
    )
  }

  ensembl_ids <- as.character(
    feature_table[[1]]
  )

  if (anyNA(ensembl_ids)) {
    stop(
      "NA Ensembl IDs found in ",
      sample_id
    )
  }

  if (anyDuplicated(ensembl_ids)) {
    stop(
      "Duplicated Ensembl IDs found in ",
      sample_id
    )
  }

  # ----------------------------------------------------------
  # B. Read raw 10x matrix
  # IMPORTANT: feature.column = 1
  # ----------------------------------------------------------

  counts <- ReadMtx(
    mtx = matrix_path,
    cells = barcode_path,
    features = features_path,
    cell.column = 1,
    feature.column = 1,
    unique.features = FALSE,
    strip.suffix = FALSE
  )

  if (!inherits(counts, "sparseMatrix")) {
    stop(
      "Matrix is not sparse for sample ",
      sample_id
    )
  }

  if (nrow(counts) != length(ensembl_ids)) {
    stop(
      "Matrix/features row-count mismatch for ",
      sample_id
    )
  }

  if (!identical(
    rownames(counts),
    ensembl_ids
  )) {
    stop(
      "Matrix rownames are not exact features.tsv.gz column-1 ",
      "Ensembl IDs for ",
      sample_id
    )
  }

  if (anyDuplicated(rownames(counts))) {
    stop(
      "Duplicated matrix rownames in ",
      sample_id
    )
  }

  raw_barcodes <- colnames(counts)

  if (anyDuplicated(raw_barcodes)) {
    stop(
      "Duplicated raw barcodes within ",
      sample_id
    )
  }

  native_features <- nrow(counts)
  raw_cells <- ncol(counts)

  expected_native_features <- if (
    sample_id == "S2487T"
  ) {
    61541L
  } else {
    36601L
  }

  if (
    native_features !=
    expected_native_features
  ) {

    cat("\n")
    cat(
      "STEP18B_NATIVE_SEURAT_GATE = FAIL\n"
    )
    cat(
      "REASON = NATIVE_FEATURE_COUNT_MISMATCH\n"
    )

    stop(
      sample_id,
      ": expected ",
      expected_native_features,
      " features, found ",
      native_features
    )
  }

  # ----------------------------------------------------------
  # C. Verify frozen common feature set
  # ----------------------------------------------------------

  common_present <- sum(
    common_features %in%
      rownames(counts)
  )

  common_missing <- length(
    setdiff(
      common_features,
      rownames(counts)
    )
  )

  outside_common <- length(
    setdiff(
      rownames(counts),
      common_features
    )
  )

  if (
    common_present != 36476 ||
    common_missing != 0
  ) {
    stop(
      sample_id,
      ": frozen common feature set incomplete."
    )
  }

  # ----------------------------------------------------------
  # D. Globally unique cell IDs
  # ----------------------------------------------------------

  global_cell_ids <- paste0(
    sample_id,
    "_",
    raw_barcodes
  )

  if (anyDuplicated(global_cell_ids)) {
    stop(
      "Duplicate global cell IDs within sample ",
      sample_id
    )
  }

  overlap_global <- intersect(
    global_cell_ids,
    all_global_cell_ids
  )

  if (length(overlap_global) > 0) {
    stop(
      "Global cell-ID collision detected."
    )
  }

  all_global_cell_ids <- c(
    all_global_cell_ids,
    global_cell_ids
  )

  colnames(counts) <- global_cell_ids

  # ----------------------------------------------------------
  # E. Find exactly one crosswalk record
  # ----------------------------------------------------------

  cw <- NULL

  if (!is.na(sample_crosswalk_col)) {

    m_idx <- which(
      as.character(
        crosswalk[[sample_crosswalk_col]]
      ) == sample_id
    )

    if (length(m_idx) == 1) {
      cw <- crosswalk[m_idx, ]
    }
  }

  if (
    is.null(cw) &&
    !is.na(gsm_crosswalk_col) &&
    !is.na(GSM)
  ) {

    m_idx <- which(
      as.character(
        crosswalk[[gsm_crosswalk_col]]
      ) == GSM
    )

    if (length(m_idx) == 1) {
      cw <- crosswalk[m_idx, ]
    }
  }

  if (is.null(cw)) {
    stop(
      "Crosswalk did not give exactly one record for ",
      sample_id
    )
  }

  # ----------------------------------------------------------
  # F. Cell metadata
  # ----------------------------------------------------------

  md <- data.frame(
    row.names = global_cell_ids,
    sample_id = rep(
      sample_id,
      raw_cells
    ),
    GSM = rep(
      GSM,
      raw_cells
    ),
    raw_barcode = raw_barcodes,
    stringsAsFactors = FALSE
  )

  if (length(prov_columns) > 0) {

    for (new_name in names(prov_columns)) {

      source_name <- prov_columns[[new_name]]

      md[[new_name]] <- rep(
        as.character(cw[[source_name]][1]),
        raw_cells
      )
    }
  }

  # ----------------------------------------------------------
  # G. Create raw Seurat object
  # ----------------------------------------------------------

  obj <- CreateSeuratObject(
    counts = counts,
    project = "GSE221561",
    min.cells = 0,
    min.features = 0,
    meta.data = md
  )

  if (ncol(obj) != raw_cells) {
    stop(
      "CreateSeuratObject changed cell count for ",
      sample_id
    )
  }

  if (nrow(obj) != native_features) {
    stop(
      "CreateSeuratObject changed feature count for ",
      sample_id
    )
  }

  # ----------------------------------------------------------
  # H. Frozen feature-space provenance
  # ----------------------------------------------------------

  obj@misc$GSE221561_feature_policy <- list(
    feature_identity_key = "Ensembl_ID",
    native_feature_count = native_features,
    cross_sample_feature_strategy = "INTERSECTION",
    cross_sample_common_features = 36476L,
    native_feature_space_preserved = "YES",
    normalized = "NO",
    qc_filtered = "NO",
    author_celltype_attached = "NO"
  )

  # ----------------------------------------------------------
  # I. Verify NO downstream analysis
  # ----------------------------------------------------------

  if (length(Reductions(obj)) != 0) {
    stop(
      "Unexpected dimensional reduction in ",
      sample_id
    )
  }

  if (length(Graphs(obj)) != 0) {
    stop(
      "Unexpected graph in ",
      sample_id
    )
  }

  vf <- VariableFeatures(obj)

  if (length(vf) != 0) {
    stop(
      "Unexpected variable features in ",
      sample_id
    )
  }

  # ----------------------------------------------------------
  # J. Save exactly one native object
  # ----------------------------------------------------------

  object_file <- file.path(
    OBJECT_DIR,
    paste0(
      sample_id,
      "_native_raw.rds"
    )
  )

  object_status <- "CREATED"

  if (file.exists(object_file)) {

    old <- readRDS(object_file)

    same_sample <- all(
      unique(old$sample_id) == sample_id
    )

    same_cells <- ncol(old) == raw_cells
    same_features <- nrow(old) == native_features

    old_feature_ids_ok <- identical(
      rownames(old),
      rownames(obj)
    )

    if (
      same_sample &&
      same_cells &&
      same_features &&
      old_feature_ids_ok
    ) {

      cat(
        "Existing object is already structurally valid.\n"
      )

      object_status <- "ALREADY_VALID"

      rm(old)

    } else {

      stop(
        "REVIEW_REQUIRED: conflicting existing object:\n",
        object_file
      )
    }

  } else {

    saveRDS(
      obj,
      object_file,
      compress = FALSE
    )
  }

  object_size <- file.info(
    object_file
  )$size

  object_sha256 <- sha256_file(
    object_file
  )

  # ----------------------------------------------------------
  # K. Pull sample provenance into audit
  # ----------------------------------------------------------

  patient_val <- NA_character_
  tissue_val <- NA_character_
  treatment_val <- NA_character_

  if ("patient_id" %in% names(md)) {
    patient_val <- unique(md$patient_id)[1]
  }

  if ("tissue" %in% names(md)) {
    tissue_val <- unique(md$tissue)[1]
  }

  if ("treatment" %in% names(md)) {
    treatment_val <- unique(md$treatment)[1]
  }

  audit_rows[[sample_id]] <- data.table(
    sample_id = sample_id,
    GSM = GSM,
    patient_id = patient_val,
    tissue = tissue_val,
    treatment = treatment_val,
    native_features = native_features,
    common_features_present = common_present,
    common_features_missing = common_missing,
    features_outside_common = outside_common,
    raw_cells = raw_cells,
    object_cells = ncol(obj),
    feature_identity_key = "Ensembl_ID",
    cell_id_policy = "<sample_id>_<raw_barcode>",
    counts_preserved = "YES",
    normalized = "NO",
    qc_filtered = "NO",
    author_celltype_attached = "NO",
    object_path = file.path(
      "03_objects/GSE221561/native_raw",
      basename(object_file)
    ),
    object_size_bytes = object_size,
    object_sha256 = object_sha256,
    status = object_status
  )

  total_raw_cells <- total_raw_cells +
    raw_cells

  objects_ok <- c(
    objects_ok,
    sample_id
  )

  cat(
    "Native features       :",
    native_features,
    "\n"
  )

  cat(
    "Raw cells             :",
    raw_cells,
    "\n"
  )

  cat(
    "Common present        :",
    common_present,
    "\n"
  )

  cat(
    "Outside common        :",
    outside_common,
    "\n"
  )

  cat(
    "Object status         :",
    object_status,
    "\n"
  )

  rm(
    counts,
    obj
  )

  gc()
}

# ------------------------------------------------------------
# 9. Combined audit
# ------------------------------------------------------------

audit <- rbindlist(
  audit_rows,
  fill = TRUE
)

setorder(
  audit,
  sample_id
)

audit_file <- paste0(
  "00_metadata/",
  "07_GSE221561_native_seurat_object_audit.csv"
)

fwrite(
  audit,
  audit_file,
  bom = TRUE
)

cat("\n============================================================\n")
cat("NATIVE OBJECT AUDIT\n")
cat("============================================================\n\n")

print(
  audit[
    ,
    .(
      sample_id,
      native_features,
      common_features_present,
      common_features_missing,
      features_outside_common,
      raw_cells,
      object_cells,
      status
    )
  ]
)

# ------------------------------------------------------------
# 10. Cross-object validation
# ------------------------------------------------------------

global_ids_unique <- (
  length(all_global_cell_ids) ==
    uniqueN(all_global_cell_ids)
)

native_objects <- list.files(
  OBJECT_DIR,
  pattern = "_native_raw\\.rds$",
  full.names = TRUE
)

all_common_ok <- all(
  audit$common_features_present == 36476 &
    audit$common_features_missing == 0
)

native_spaces_ok <- all(
  audit$native_features ==
    ifelse(
      audit$sample_id == "S2487T",
      61541,
      36601
    )
)

raw_counts_ok <- all(
  audit$counts_preserved == "YES"
)

# ------------------------------------------------------------
# 11. Provenance completeness checks
# ------------------------------------------------------------

sample_provenance_complete <-
  all(!is.na(audit$sample_id)) &&
  all(nzchar(audit$sample_id))

patient_mapping_preserved <-
  if ("patient_id" %in% names(audit)) {
    all(
      !is.na(audit$patient_id) &
        nzchar(audit$patient_id)
    )
  } else {
    FALSE
  }

tissue_mapping_preserved <-
  if ("tissue" %in% names(audit)) {
    all(
      !is.na(audit$tissue) &
        nzchar(audit$tissue)
    )
  } else {
    FALSE
  }

treatment_mapping_preserved <-
  if ("treatment" %in% names(audit)) {
    all(
      !is.na(audit$treatment) &
        nzchar(audit$treatment)
    )
  } else {
    FALSE
  }

# These must remain unresolved if that was the
# previously frozen crosswalk state.
timepoint_unresolved_preserved <- TRUE
response_unresolved_preserved <- TRUE

if ("timepoint" %in% names(crosswalk)) {

  vals <- unique(
    toupper(
      trimws(
        as.character(
          crosswalk$timepoint
        )
      )
    )
  )

  vals <- vals[
    !is.na(vals) &
      nzchar(vals)
  ]

  if (
    length(vals) > 0 &&
    !all(vals %in% c(
      "UNRESOLVED",
      "UNKNOWN",
      "NA",
      "N/A"
    ))
  ) {
    timepoint_unresolved_preserved <- FALSE
  }
}

if ("response" %in% names(crosswalk)) {

  vals <- unique(
    toupper(
      trimws(
        as.character(
          crosswalk$response
        )
      )
    )
  )

  vals <- vals[
    !is.na(vals) &
      nzchar(vals)
  ]

  if (
    length(vals) > 0 &&
    !all(vals %in% c(
      "UNRESOLVED",
      "UNKNOWN",
      "NA",
      "N/A"
    ))
  ) {
    response_unresolved_preserved <- FALSE
  }
}

# ------------------------------------------------------------
# 12. Save sessionInfo
# ------------------------------------------------------------

sink(
  "08_logs/07_GSE221561_native_seurat_sessionInfo.txt"
)

print(
  sessionInfo()
)

sink()

# ------------------------------------------------------------
# 13. Final gate
# ------------------------------------------------------------

gate <- "PASS"

if (
  length(native_objects) != 11 ||
  length(objects_ok) != 11 ||
  !all_common_ok ||
  !native_spaces_ok ||
  !raw_counts_ok ||
  !global_ids_unique ||
  !sample_provenance_complete ||
  !patient_mapping_preserved ||
  !tissue_mapping_preserved ||
  !treatment_mapping_preserved
) {
  gate <- "REVIEW_REQUIRED"
}

cat("\n")
cat("============================================================\n")
cat("STEP 18B FINAL GATE\n")
cat("============================================================\n\n")

cat(
  "STEP18B_NATIVE_SEURAT_GATE = ",
  gate,
  "\n",
  sep = ""
)

cat(
  "NATIVE_SEURAT_OBJECTS = ",
  length(native_objects),
  "/11\n",
  sep = ""
)

cat(
  "TOTAL_RAW_MATRIX_CELLS = ",
  total_raw_cells,
  "\n",
  sep = ""
)

cat(
  "AUTHOR_METADATA_CELLS = 55150\n"
)

cat("\n")

cat(
  "FEATURE_IDENTITY_KEY = Ensembl_ID\n"
)

cat(
  "S2487T_NATIVE_FEATURES = ",
  audit[
    sample_id == "S2487T",
    native_features
  ][1],
  "\n",
  sep = ""
)

other_native <- unique(
  audit[
    sample_id != "S2487T",
    native_features
  ]
)

cat(
  "OTHER_SAMPLE_NATIVE_FEATURES = ",
  paste(
    other_native,
    collapse = ","
  ),
  "\n",
  sep = ""
)

cat(
  "COMMON_FEATURES_EXPECTED = 36476\n"
)

cat(
  "ALL_OBJECTS_CONTAIN_ALL_COMMON_FEATURES = ",
  ifelse(all_common_ok, "YES", "NO"),
  "\n",
  sep = ""
)

cat("\n")

cat(
  "FULL_NATIVE_FEATURE_SPACES_PRESERVED = ",
  ifelse(native_spaces_ok, "YES", "NO"),
  "\n",
  sep = ""
)

cat(
  "RAW_COUNTS_PRESERVED = ",
  ifelse(raw_counts_ok, "YES", "NO"),
  "\n",
  sep = ""
)

cat("\n")

cat(
  "GLOBAL_CELL_IDS_UNIQUE = ",
  ifelse(global_ids_unique, "YES", "NO"),
  "\n",
  sep = ""
)

cat(
  "SAMPLE_PROVENANCE_COMPLETE = ",
  ifelse(
    sample_provenance_complete,
    "YES",
    "NO"
  ),
  "\n",
  sep = ""
)

cat(
  "PATIENT_MAPPING_PRESERVED = ",
  ifelse(
    patient_mapping_preserved,
    "YES",
    "NO"
  ),
  "\n",
  sep = ""
)

cat(
  "TISSUE_MAPPING_PRESERVED = ",
  ifelse(
    tissue_mapping_preserved,
    "YES",
    "NO"
  ),
  "\n",
  sep = ""
)

cat(
  "TREATMENT_REGIMEN_PRESERVED = ",
  ifelse(
    treatment_mapping_preserved,
    "YES",
    "NO"
  ),
  "\n",
  sep = ""
)

cat(
  "TIMEPOINT_UNRESOLVED_PRESERVED = ",
  ifelse(
    timepoint_unresolved_preserved,
    "YES",
    "NO"
  ),
  "\n",
  sep = ""
)

cat(
  "RESPONSE_UNRESOLVED_PRESERVED = ",
  ifelse(
    response_unresolved_preserved,
    "YES",
    "NO"
  ),
  "\n",
  sep = ""
)

cat("\n")

cat(
  "AUTHOR_CELLTYPE_ATTACHED = NO\n"
)

cat(
  "NORMALIZATION_PERFORMED = NO\n"
)

cat(
  "QC_FILTERING_PERFORMED = NO\n"
)

cat(
  "MERGED_OBJECT_CREATED = NO\n"
)

cat(
  "HARMONIZED_ASSAY_CREATED = NO\n"
)

cat("\n")

cat(
  "Audit saved: ",
  audit_file,
  "\n",
  sep = ""
)

cat(
  "Objects saved: ",
  "03_objects/GSE221561/native_raw/\n",
  sep = ""
)

cat(
  "SessionInfo saved: ",
  "08_logs/07_GSE221561_native_seurat_sessionInfo.txt\n",
  sep = ""
)

cat("\n============================================================\n")
cat("STEP 18B FINISHED\n")
cat("============================================================\n")
