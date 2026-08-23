# ============================================================
# ESCC Neoadjuvant Project
# STEP22A-6D
#
# BUILD ALL-27 FIBROBLAST RAW-COUNT PSEUDOBULKS
#
# Source of fibroblast identity:
#   STEP22A-6C frozen Fibroblast_CAF cell-ID list
#
# Raw expression source:
#   original 27 sparse count matrices
#
# Feature universe:
#   true all-27 common ENSG universe = 16,790
#
# IMPORTANT:
#   CORE genes are restored here as RAW COUNT ROWS only.
#
# This step DOES NOT:
#   - use CORE genes for cell selection
#   - inspect/print gene-specific CORE values
#   - normalize counts
#   - calculate CPM
#   - calculate log2(CPM+1)
#   - calculate z-scores
#   - calculate CORE2 / CORE4
#   - perform Before/After inference
#
# Frozen fibroblast technical threshold:
#   >=20 cells per patient-timepoint
#
# Histology remains unresolved.
# ============================================================

options(
  stringsAsFactors = FALSE,
  timeout = 1800
)

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
})

cat("\n")
cat("====================================================================\n")
cat("STEP22A-6D: ALL27 FIBROBLAST RAW PSEUDOBULK\n")
cat("====================================================================\n\n")

# ============================================================
# 0. Paths
# ============================================================

PROJECT <- getwd()

META_DIR <- file.path(
  PROJECT,
  "00_metadata",
  "OMIX005710"
)

OBJECT_DIR <- file.path(
  PROJECT,
  "03_objects",
  "OMIX005710"
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

FIBRO_IDS_FILE <- file.path(
  META_DIR,
  "STEP22A_06C_ALL27_FIBROBLAST_CELL_IDS.csv"
)

FIBRO_QC_FILE <- file.path(
  META_DIR,
  "STEP22A_06C_ALL27_FIBROBLAST_SAMPLE_QC.csv"
)

COMMON_FILE <- file.path(
  META_DIR,
  "STEP22A_06B_ALL27_COMMON_ENSG_FEATURES.csv"
)

COMBINED_INV <- file.path(
  META_DIR,
  "STEP22A_06A_COMBINED_27_TUMOR_TECHNICAL_INVENTORY.csv"
)

PAIRED_AUDIT <- file.path(
  META_DIR,
  "STEP22A_03B_MATRIX_AUDIT.csv"
)

REMAINING_AUDIT <- file.path(
  META_DIR,
  "STEP22A_06B_REMAINING_15_MATRIX_AUDIT.csv"
)

GATE_FILE <- file.path(
  META_DIR,
  "STEP22A_ANALYSIS_GATE_STATUS.csv"
)

FROZEN_INV <- file.path(
  META_DIR,
  "STEP22A_PREANALYSIS_FROZEN_FILE_INVENTORY.csv"
)

SAMPLE_META_OUT <- file.path(
  META_DIR,
  "STEP22A_06D_ALL27_FIBROBLAST_SAMPLE_PSEUDOBULK_METADATA.csv"
)

PT_META_OUT <- file.path(
  META_DIR,
  "STEP22A_06D_ALL27_FIBROBLAST_PATIENT_TIMEPOINT_METADATA.csv"
)

PAIR_STATUS_OUT <- file.path(
  META_DIR,
  "STEP22A_06D_PAIRED_FIBROBLAST_TECHNICAL_STATUS.csv"
)

PROVENANCE_OUT <- file.path(
  META_DIR,
  "STEP22A_06D_PSEUDOBULK_PROVENANCE.csv"
)

SUMMARY_OUT <- file.path(
  META_DIR,
  "STEP22A_06D_ALL27_FIBROBLAST_PSEUDOBULK_SUMMARY.txt"
)

BUNDLE_OUT <- file.path(
  OBJECT_DIR,
  paste0(
    "OMIX005710_ALL27_FIBROBLAST_RAW_PSEUDOBULK_",
    "HISTOLOGY_UNRESOLVED.rds"
  )
)

# ============================================================
# Helpers
# ============================================================

read_gene_table <- function(path) {

  x <- fread(
    path,
    header = FALSE,
    encoding = "UTF-8"
  )

  if (ncol(x) < 1) {

    stop(
      "Empty gene table: ",
      path
    )
  }

  gene_id <- as.character(
    x[[1]]
  )

  if (ncol(x) >= 2) {

    gene_symbol <- as.character(
      x[[2]]
    )

  } else {

    gene_symbol <- gene_id
  }

  bad <- (
    is.na(gene_symbol) |
    gene_symbol == ""
  )

  gene_symbol[bad] <- gene_id[bad]

  if (
    anyDuplicated(
      gene_id
    ) > 0
  ) {

    stop(
      "Duplicate ENSG IDs in ",
      path
    )
  }

  data.table(
    gene_id = gene_id,
    gene_symbol = gene_symbol
  )
}


patient_number <- function(x) {

  suppressWarnings(
    as.integer(
      sub(
        "^P",
        "",
        x
      )
    )
  )
}

# ============================================================
# 1. Gate check
# ============================================================

cat("===== 1. ANALYSIS GATE =====\n\n")

gate_dt <- fread(
  GATE_FILE,
  encoding = "UTF-8"
)

get_gate <- function(x) {

  z <- gate_dt[
    gate == x,
    allowed
  ]

  if (length(z) != 1) {

    stop(
      "Cannot resolve gate: ",
      x
    )
  }

  as.character(z)
}

technical_allowed <- get_gate(
  "technical_preprocessing"
)

core_allowed <- get_gate(
  "paired_CORE2_test"
)

manuscript_allowed <- get_gate(
  "manuscript_update"
)

cat(
  "Technical preprocessing :",
  technical_allowed,
  "\n"
)

cat(
  "Paired CORE2 inference  :",
  core_allowed,
  "\n"
)

cat(
  "Manuscript update       :",
  manuscript_allowed,
  "\n\n"
)

if (
  technical_allowed != "YES"
) {

  stop(
    "TECHNICAL_PREPROCESSING_BLOCKED"
  )
}

if (
  core_allowed == "YES"
) {

  stop(
    "CORE2 inference gate unexpectedly open"
  )
}

if (
  manuscript_allowed == "YES"
) {

  stop(
    "Manuscript gate unexpectedly open"
  )
}

cat(
  "✓ Raw-count pseudobulk preparation allowed.\n"
)

cat(
  "✓ CORE2 inference remains blocked.\n\n"
)

# ============================================================
# 2. Load frozen Step22A-6C fibroblast identities
# ============================================================

cat("===== 2. LOAD FROZEN FIBROBLAST CELL IDs =====\n\n")

fibro_ids <- fread(
  FIBRO_IDS_FILE,
  encoding = "UTF-8"
)

fibro_qc_6c <- fread(
  FIBRO_QC_FILE,
  encoding = "UTF-8"
)

if (
  nrow(
    fibro_ids
  ) != 33385
) {

  cat(
    "WARNING: Step22A-6C reported 33,385 fibroblasts; ",
    "current file contains ",
    nrow(fibro_ids),
    ".\n",
    sep = ""
  )
}

if (
  anyDuplicated(
    fibro_ids$cell_id
  ) > 0
) {

  stop(
    "DUPLICATE_FROZEN_FIBROBLAST_CELL_ID"
  )
}

if (
  nrow(
    fibro_qc_6c
  ) != 27
) {

  stop(
    "Expected 27 Step22A-6C fibroblast-QC rows"
  )
}

cat(
  "Frozen Fibroblast_CAF cells:",
  nrow(
    fibro_ids
  ),
  "\n"
)

cat(
  "Tumor observations:",
  nrow(
    fibro_qc_6c
  ),
  "\n\n"
)

cat(
  "Cell-selection source: STEP22A-6C ONLY\n"
)

cat(
  "CORE genes used for cell selection: NO\n\n"
)

# ============================================================
# 3. Load full 16,790 ENSG universe
# ============================================================

cat("===== 3. LOAD TRUE ALL27 COMMON ENSG UNIVERSE =====\n\n")

common <- fread(
  COMMON_FILE,
  encoding = "UTF-8"
)

setorder(
  common,
  all27_common_index
)

if (
  nrow(
    common
  ) != 16790
) {

  cat(
    "WARNING: expected 16,790 features; observed ",
    nrow(common),
    ".\n",
    sep = ""
  )
}

common_ids <- common$gene_id

if (
  anyDuplicated(
    common_ids
  ) > 0
) {

  stop(
    "DUPLICATED_ALL27_COMMON_ENSG"
  )
}

core_genes <- c(
  "CDKN1B",
  "GSN",
  "BGN",
  "TIMP1"
)

core_rows <- common[
  toupper(
    gene_symbol
  ) %in%
    core_genes
]

if (
  nrow(core_rows) != 4 ||
  uniqueN(
    toupper(
      core_rows$gene_symbol
    )
  ) != 4
) {

  stop(
    "CORE_GENE_MAPPING_NOT_4_OF_4"
  )
}

cat(
  "Full pseudobulk features:",
  nrow(common),
  "\n"
)

cat(
  "CORE rows restored in raw-count structure:",
  nrow(core_rows),
  "/ 4\n"
)

cat(
  "CORE gene-specific values inspected:",
  "NO\n\n"
)

# ============================================================
# 4. Build 27 matrix path table
# ============================================================

cat("===== 4. LOAD ALL27 MATRIX PATHS =====\n\n")

combined <- fread(
  COMBINED_INV,
  encoding = "UTF-8"
)

paired_audit <- fread(
  PAIRED_AUDIT,
  encoding = "UTF-8"
)[
  ,
  .(
    sample_name,
    matrix_file,
    barcodes_file,
    genes_file
  )
]

remaining_audit <- fread(
  REMAINING_AUDIT,
  encoding = "UTF-8"
)[
  ,
  .(
    sample_name,
    matrix_file,
    barcodes_file,
    genes_file
  )
]

audit <- rbindlist(
  list(
    paired_audit,
    remaining_audit
  ),
  use.names = TRUE
)

if (
  nrow(combined) != 27 ||
  nrow(audit) != 27
) {

  stop(
    "ALL27_INPUT_COUNT_INVALID"
  )
}

audit <- audit[
  match(
    combined$sample_name,
    audit$sample_name
  )
]

if (
  anyNA(
    audit$sample_name
  ) ||
  !identical(
    audit$sample_name,
    combined$sample_name
  )
) {

  stop(
    "ALL27_MATRIX_PATH_ALIGNMENT_FAILED"
  )
}

sample_table <- cbind(
  combined[
    ,
    .(
      technical_set,
      patient_id,
      timepoint,
      sample_name,
      hrs_accession
    )
  ],
  audit[
    ,
    .(
      matrix_file,
      barcodes_file,
      genes_file
    )
  ]
)

cat(
  "Tumor observations:",
  nrow(sample_table),
  "\n"
)

cat(
  "Patients:",
  uniqueN(
    sample_table$patient_id
  ),
  "\n\n"
)

# ============================================================
# 5. Rebuild raw-count pseudobulk one sample at a time
#
# Memory-safe:
#   only one raw sparse matrix is resident at a time.
# ============================================================

cat("===== 5. REBUILD ALL27 RAW-COUNT PSEUDOBULKS =====\n\n")

pb_vectors <- list()
sample_meta_list <- list()

for (
  i in seq_len(
    nrow(sample_table)
  )
) {

  row <- sample_table[i]

  sample_name <- row$sample_name

  matrix_path <- file.path(
    PROJECT,
    row$matrix_file
  )

  barcode_path <- file.path(
    PROJECT,
    row$barcodes_file
  )

  gene_path <- file.path(
    PROJECT,
    row$genes_file
  )

  for (
    f in c(
      matrix_path,
      barcode_path,
      gene_path
    )
  ) {

    if (!file.exists(f)) {

      stop(
        "Missing input file: ",
        f
      )
    }
  }

  gt <- read_gene_table(
    gene_path
  )

  barcodes <- as.character(
    fread(
      barcode_path,
      header = FALSE,
      encoding = "UTF-8"
    )[[1]]
  )

  if (
    anyDuplicated(
      barcodes
    ) > 0
  ) {

    stop(
      "Duplicate raw barcodes in ",
      sample_name
    )
  }

  sample_fib <- fibro_ids[
    sample_name ==
      row$sample_name
  ]

  n_fib <- nrow(
    sample_fib
  )

  frozen_qc_count <- fibro_qc_6c[
    sample_name ==
      row$sample_name,
    n_fibroblast_cells
  ]

  if (
    length(
      frozen_qc_count
    ) != 1
  ) {

    stop(
      "Cannot resolve frozen fibroblast count for ",
      sample_name
    )
  }

  if (
    n_fib !=
      frozen_qc_count
  ) {

    stop(
      "FROZEN_FIBROBLAST_ID_COUNT_MISMATCH: ",
      sample_name,
      " IDs=",
      n_fib,
      " QC=",
      frozen_qc_count
    )
  }

  # Verify cell IDs are internally consistent.
  expected_cell_ids <- paste0(
    sample_name,
    "__",
    sample_fib$barcode_raw
  )

  if (
    !identical(
      expected_cell_ids,
      sample_fib$cell_id
    )
  ) {

    stop(
      "FROZEN_CELL_ID_BARCODE_MISMATCH: ",
      sample_name
    )
  }

  barcode_idx <- match(
    sample_fib$barcode_raw,
    barcodes
  )

  if (
    anyNA(
      barcode_idx
    )
  ) {

    stop(
      "FROZEN_FIBROBLAST_BARCODE_NOT_FOUND: ",
      sample_name
    )
  }

  if (
    anyDuplicated(
      barcode_idx
    ) > 0
  ) {

    stop(
      "DUPLICATE_SELECTED_BARCODE_INDEX: ",
      sample_name
    )
  }

  common_idx <- match(
    common_ids,
    gt$gene_id
  )

  if (
    anyNA(
      common_idx
    )
  ) {

    stop(
      "ALL27_COMMON_ENSG_MISSING: ",
      sample_name
    )
  }

  counts_raw <- Matrix::readMM(
    matrix_path
  )

  counts_raw <- as(
    counts_raw,
    "dgCMatrix"
  )

  if (
    nrow(counts_raw) !=
      nrow(gt) ||
    ncol(counts_raw) !=
      length(barcodes)
  ) {

    stop(
      "RAW_MATRIX_DIMENSION_MISMATCH: ",
      sample_name
    )
  }

  # Reorder to the exact 16,790 ENSG common universe.
  counts_common <- counts_raw[
    common_idx,
    ,
    drop = FALSE
  ]

  rm(
    counts_raw
  )

  # Select ONLY frozen Step22A-6C fibroblast cells.
  fib_counts <- counts_common[
    ,
    barcode_idx,
    drop = FALSE
  ]

  rm(
    counts_common
  )

  # Raw integer pseudobulk sum.
  pb <- as.numeric(
    Matrix::rowSums(
      fib_counts
    )
  )

  rm(
    fib_counts
  )

  if (
    any(
      !is.finite(pb)
    ) ||
    any(
      pb < 0
    )
  ) {

    stop(
      "INVALID_PSEUDOBULK_COUNTS: ",
      sample_name
    )
  }

  if (
    any(
      abs(
        pb -
          round(pb)
      ) >
        1e-8
    )
  ) {

    stop(
      "NON_INTEGER_RAW_PSEUDOBULK_COUNTS: ",
      sample_name
    )
  }

  library_size <- sum(
    pb
  )

  detected_genes <- sum(
    pb > 0
  )

  technical_status <- if (
    n_fib >= 20
  ) {
    "PASS_GE20"
  } else {
    "FAIL_LT20"
  }

  pb_vectors[[
    sample_name
  ]] <- pb

  sample_meta_list[[
    sample_name
  ]] <- data.table(
    technical_set =
      row$technical_set,
    patient_id =
      row$patient_id,
    timepoint =
      row$timepoint,
    sample_name =
      sample_name,
    hrs_accession =
      row$hrs_accession,
    n_fibroblast_cells =
      n_fib,
    raw_library_size =
      library_size,
    detected_genes =
      detected_genes,
    frozen_technical_status =
      technical_status,
    histology_status =
      "UNRESOLVED",
    primary_scoring_allowed =
      "NO"
  )

  cat(
    sprintf(
      "[%02d/27] %-10s %-4s %-6s fibro=%5d  library=%10.0f  genes=%5d  %s\n",
      i,
      sample_name,
      row$patient_id,
      row$timepoint,
      n_fib,
      library_size,
      detected_genes,
      technical_status
    )
  )

  gc(
    verbose = FALSE
  )
}

cat("\n")

# ============================================================
# 6. Assemble sample-level raw-count matrix
# ============================================================

cat("===== 6. ASSEMBLE SAMPLE-LEVEL PSEUDOBULK MATRIX =====\n\n")

sample_names <- sample_table$sample_name

sample_pb <- do.call(
  cbind,
  pb_vectors[
    sample_names
  ]
)

rownames(
  sample_pb
) <- common_ids

colnames(
  sample_pb
) <- sample_names

rm(
  pb_vectors
)

sample_meta <- rbindlist(
  sample_meta_list[
    sample_names
  ],
  use.names = TRUE
)

rm(
  sample_meta_list
)

if (
  nrow(sample_pb) !=
    nrow(common) ||
  ncol(sample_pb) != 27
) {

  stop(
    "SAMPLE_PSEUDOBULK_DIMENSION_INVALID"
  )
}

if (
  !identical(
    colnames(sample_pb),
    sample_meta$sample_name
  )
) {

  stop(
    "SAMPLE_PSEUDOBULK_METADATA_ORDER_MISMATCH"
  )
}

cat(
  "Sample pseudobulk matrix:",
  nrow(sample_pb),
  "genes x",
  ncol(sample_pb),
  "tumor observations\n\n"
)

# ============================================================
# 7. Collapse to patient × timepoint
#
# If multiple source samples ever exist for the same
# patient-timepoint, raw counts are summed.
# ============================================================

cat("===== 7. PATIENT × TIMEPOINT COLLAPSE =====\n\n")

sample_meta[
  ,
  patient_timepoint :=
    paste(
      patient_id,
      timepoint,
      sep = "__"
    )
]

sample_meta[
  ,
  patient_order :=
    patient_number(
      patient_id
    )
]

sample_meta[
  ,
  timepoint_order :=
    fifelse(
      timepoint ==
        "Before",
      1L,
      2L
    )
]

pt_levels <- sample_meta[
  order(
    patient_order,
    timepoint_order,
    sample_name
  ),
  unique(
    patient_timepoint
  )
]

collapse_design <- Matrix::sparseMatrix(
  i = seq_len(
    nrow(sample_meta)
  ),
  j = match(
    sample_meta$patient_timepoint,
    pt_levels
  ),
  x = 1,
  dims = c(
    nrow(sample_meta),
    length(pt_levels)
  ),
  dimnames = list(
    sample_meta$sample_name,
    pt_levels
  )
)

if (
  !identical(
    rownames(collapse_design),
    colnames(sample_pb)
  )
) {

  stop(
    "PATIENT_TIMEPOINT_COLLAPSE_ORDER_MISMATCH"
  )
}

pt_pb <- sample_pb %*%
  collapse_design

pt_pb <- as.matrix(
  pt_pb
)

pt_meta <- sample_meta[
  ,
  .(
    technical_set =
      paste(
        unique(
          technical_set
        ),
        collapse = ";"
      ),
    number_of_source_samples =
      .N,
    source_samples =
      paste(
        sample_name,
        collapse = ";"
      ),
    source_HRS =
      paste(
        hrs_accession,
        collapse = ";"
      ),
    n_fibroblast_cells =
      sum(
        n_fibroblast_cells
      )
  ),
  by = .(
    patient_id,
    timepoint,
    patient_timepoint
  )
]

pt_meta[
  ,
  patient_order :=
    patient_number(
      patient_id
    )
]

pt_meta[
  ,
  timepoint_order :=
    fifelse(
      timepoint ==
        "Before",
      1L,
      2L
    )
]

setorder(
  pt_meta,
  patient_order,
  timepoint_order,
  patient_timepoint
)

pt_meta[
  ,
  c(
    "patient_order",
    "timepoint_order"
  ) := NULL
]

pt_meta <- pt_meta[
  match(
    colnames(pt_pb),
    patient_timepoint
  )
]

if (
  anyNA(
    pt_meta$patient_timepoint
  )
) {

  stop(
    "PATIENT_TIMEPOINT_METADATA_ALIGNMENT_FAILED"
  )
}

pt_meta[
  ,
  raw_library_size :=
    as.numeric(
      colSums(
        pt_pb
      )
    )
]

pt_meta[
  ,
  detected_genes :=
    as.integer(
      colSums(
        pt_pb > 0
      )
    )
]

pt_meta[
  ,
  frozen_technical_status :=
    fifelse(
      n_fibroblast_cells >= 20,
      "PASS_GE20",
      "FAIL_LT20"
    )
]

pt_meta[
  ,
  technical_reference_candidate :=
    fifelse(
      frozen_technical_status ==
        "PASS_GE20",
      "YES",
      "NO"
    )
]

pt_meta[
  ,
  histology_status :=
    "UNRESOLVED"
]

pt_meta[
  ,
  final_ESCC_reference_eligible :=
    "NO_HISTOLOGY_GATE_BLOCKED"
]

pt_meta[
  ,
  CORE_scoring_allowed :=
    "NO"
]

cat(
  "Patient-timepoint pseudobulks:",
  ncol(pt_pb),
  "\n"
)

cat(
  "Patients represented:",
  uniqueN(
    pt_meta$patient_id
  ),
  "\n\n"
)

if (
  ncol(pt_pb) != 27
) {

  cat(
    "NOTE: patient-timepoint collapse produced ",
    ncol(pt_pb),
    " columns rather than 27 because one or more ",
    "patient-timepoints had multiple source samples.\n\n",
    sep = ""
  )
}

# ============================================================
# 8. Frozen >=20 technical eligibility
# ============================================================

cat("===== 8. FROZEN >=20 FIBROBLAST TECHNICAL ELIGIBILITY =====\n\n")

technical_pass_n <- sum(
  pt_meta$
    frozen_technical_status ==
    "PASS_GE20"
)

technical_fail <- pt_meta[
  frozen_technical_status ==
    "FAIL_LT20"
]

cat(
  "Technical PASS patient-timepoints:",
  technical_pass_n,
  "/",
  nrow(pt_meta),
  "\n"
)

cat(
  "Technical FAIL patient-timepoints:",
  nrow(
    technical_fail
  ),
  "\n\n"
)

if (
  nrow(
    technical_fail
  ) > 0
) {

  print(
    technical_fail[
      ,
      .(
        patient_id,
        timepoint,
        patient_timepoint,
        n_fibroblast_cells,
        frozen_technical_status
      )
    ]
  )

  cat("\n")
}

# Frozen threshold must NOT be relaxed.
if (
  any(
    technical_fail$
      n_fibroblast_cells >= 20
  )
) {

  stop(
    "TECHNICAL_THRESHOLD_LOGIC_ERROR"
  )
}

# ============================================================
# 9. Paired-patient technical status
#
# Paired means a patient has both Before + After tumor
# observations in the all-27 source structure.
#
# A technically evaluable pair requires BOTH to PASS_GE20.
# ============================================================

cat("===== 9. PAIRED TECHNICAL STATUS =====\n\n")

pair_status <- pt_meta[
  ,
  .(
    has_before =
      any(
        timepoint ==
          "Before"
      ),
    has_after =
      any(
        timepoint ==
          "After"
      ),
    before_fibroblasts =
      sum(
        n_fibroblast_cells[
          timepoint ==
            "Before"
        ]
      ),
    after_fibroblasts =
      sum(
        n_fibroblast_cells[
          timepoint ==
            "After"
        ]
      ),
    before_technical_pass =
      any(
        timepoint ==
          "Before" &
        frozen_technical_status ==
          "PASS_GE20"
      ),
    after_technical_pass =
      any(
        timepoint ==
          "After" &
        frozen_technical_status ==
          "PASS_GE20"
      )
  ),
  by = patient_id
]

pair_status[
  ,
  complete_tumor_pair :=
    has_before &
    has_after
]

pair_status[
  ,
  technically_evaluable_pair :=
    complete_tumor_pair &
    before_technical_pass &
    after_technical_pass
]

pair_status[
  ,
  patient_order :=
    patient_number(
      patient_id
    )
]

setorder(
  pair_status,
  patient_order
)

pair_status[
  ,
  patient_order := NULL
]

paired_only <- pair_status[
  complete_tumor_pair ==
    TRUE
]

n_paired <- nrow(
  paired_only
)

n_paired_technical <- sum(
  paired_only$
    technically_evaluable_pair
)

print(
  paired_only
)

cat(
  "\nComplete tumor pairs:",
  n_paired,
  "\n"
)

cat(
  "Technically evaluable pairs:",
  n_paired_technical,
  "/",
  n_paired,
  "\n\n"
)

# ============================================================
# 10. Structural raw-count integrity
# ============================================================

cat("===== 10. RAW-COUNT STRUCTURAL INTEGRITY =====\n\n")

if (
  any(
    !is.finite(
      pt_pb
    )
  ) ||
  any(
    pt_pb < 0
  )
) {

  stop(
    "INVALID_PATIENT_TIMEPOINT_RAW_COUNTS"
  )
}

if (
  any(
    abs(
      pt_pb -
        round(
          pt_pb
        )
    ) >
      1e-8
  )
) {

  stop(
    "PATIENT_TIMEPOINT_COUNTS_NOT_INTEGER_SUMS"
  )
}

if (
  any(
    colSums(
      pt_pb
    ) <= 0
  )
) {

  stop(
    "ZERO_LIBRARY_PATIENT_TIMEPOINT"
  )
}

cat(
  "All pseudobulk libraries non-zero: YES\n"
)

cat(
  "Raw counts remain integer sums: YES\n"
)

cat(
  "Normalization performed: NO\n\n"
)

# ============================================================
# 11. CORE firewall
#
# We only verify structural row presence.
# We DO NOT print/read individual pseudobulk values.
# ============================================================

cat("===== 11. CORE FIREWALL =====\n\n")

core_symbols_present <- sort(
  unique(
    toupper(
      core_rows$gene_symbol
    )
  )
)

cat(
  "CORE rows retained structurally:",
  paste(
    core_symbols_present,
    collapse = " / "
  ),
  "\n"
)

cat(
  "CORE gene-specific pseudobulk values printed:",
  "NO\n"
)

cat(
  "CORE gene-specific treatment effects inspected:",
  "NO\n"
)

cat(
  "CORE2 scoring:",
  "NO\n"
)

cat(
  "CORE4 scoring:",
  "NO\n\n"
)

# ============================================================
# 12. Save metadata
# ============================================================

sample_meta[
  ,
  c(
    "patient_order",
    "timepoint_order"
  ) := NULL
]

fwrite(
  sample_meta,
  SAMPLE_META_OUT,
  bom = TRUE
)

fwrite(
  pt_meta,
  PT_META_OUT,
  bom = TRUE
)

fwrite(
  pair_status,
  PAIR_STATUS_OUT,
  bom = TRUE
)

# ============================================================
# 13. Save compact RDS bundle
#
# 16,790 x 27 raw counts is small.
#
# Both all observations and the technical-pass mask are stored.
# Histology eligibility remains FALSE.
# ============================================================

technical_pass_mask <- (
  pt_meta$
    frozen_technical_status ==
    "PASS_GE20"
)

bundle <- list(
  raw_counts_all_patient_timepoints =
    pt_pb,
  patient_timepoint_metadata =
    as.data.frame(
      pt_meta
    ),
  sample_level_raw_counts =
    sample_pb,
  sample_metadata =
    as.data.frame(
      sample_meta
    ),
  feature_metadata =
    as.data.frame(
      common[
        ,
        .(
          all27_common_index,
          gene_id,
          gene_symbol
        )
      ]
    ),
  paired_technical_status =
    as.data.frame(
      pair_status
    ),
  technical_pass_mask =
    technical_pass_mask,
  technical_pass_patient_timepoints =
    pt_meta$patient_timepoint[
      technical_pass_mask
    ],
  frozen_rules = list(
    fibroblast_identity_source =
      "STEP22A_06C_FROZEN_CELL_IDS",
    fibroblast_minimum_cells =
      20L,
    feature_universe =
      "TRUE_ALL27_COMMON_16790_ENSG",
    count_type =
      "RAW_INTEGER_SUM",
    normalization_performed =
      FALSE,
    CORE_genes_present_in_raw_pseudobulk =
      TRUE,
    CORE_gene_specific_values_inspected =
      FALSE,
    CORE2_scored =
      FALSE,
    CORE4_scored =
      FALSE,
    histology_mapping =
      "INCOMPLETE",
    EASC_identity =
      "UNRESOLVED",
    final_ESCC_reference_allowed =
      FALSE,
    primary_inference_allowed =
      FALSE
  )
)

saveRDS(
  bundle,
  BUNDLE_OUT,
  compress = "gzip"
)

cat(
  "✓ Raw-count pseudobulk bundle saved.\n\n"
)

# ============================================================
# 14. Provenance
# ============================================================

provenance <- data.table(
  item = c(
    "fibroblast_identity_source",
    "all27_common_features",
    "sample_level_pseudobulks",
    "patient_timepoint_pseudobulks",
    "technical_PASS_GE20",
    "technical_FAIL_LT20",
    "complete_paired_patients",
    "technically_evaluable_paired_patients",
    "CORE_rows_retained_in_raw_structure",
    "CORE_gene_specific_values_printed",
    "normalization_performed",
    "CPM_calculated",
    "log2CPM_calculated",
    "zscore_calculated",
    "CORE2_calculated",
    "CORE4_calculated",
    "Before_After_inference",
    "histology_mapping",
    "EASC_identity",
    "primary_inference_allowed"
  ),
  value = c(
    "STEP22A_06C_FROZEN_FIBROBLAST_CELL_IDS",
    as.character(
      nrow(common)
    ),
    as.character(
      ncol(sample_pb)
    ),
    as.character(
      ncol(pt_pb)
    ),
    as.character(
      technical_pass_n
    ),
    as.character(
      nrow(
        technical_fail
      )
    ),
    as.character(
      n_paired
    ),
    as.character(
      n_paired_technical
    ),
    "YES_4_OF_4",
    "NO",
    "NO",
    "NO",
    "NO",
    "NO",
    "NO",
    "NO",
    "NO",
    "INCOMPLETE",
    "UNRESOLVED",
    "NO"
  )
)

fwrite(
  provenance,
  PROVENANCE_OUT,
  bom = TRUE
)

# ============================================================
# 15. Frozen project integrity
# ============================================================

cat("===== 12. FROZEN PROJECT CHECK =====\n\n")

frozen_changes <- list()

if (
  file.exists(
    FROZEN_INV
  )
) {

  frozen <- fread(
    FROZEN_INV,
    encoding = "UTF-8"
  )

  for (
    i in seq_len(
      nrow(
        frozen
      )
    )
  ) {

    path <- frozen$path[i]

    if (!file.exists(path)) {

      frozen_changes[[
        length(
          frozen_changes
        ) + 1
      ]] <- c(
        path,
        "MISSING"
      )

      next
    }

    expected <- suppressWarnings(
      as.numeric(
        frozen$size_bytes[i]
      )
    )

    actual <- file.info(
      path
    )$size

    if (
      is.finite(
        expected
      ) &&
      expected != actual
    ) {

      frozen_changes[[
        length(
          frozen_changes
        ) + 1
      ]] <- c(
        path,
        "SIZE_CHANGED"
      )
    }
  }
}

cat(
  "Frozen original files changed:",
  length(
    frozen_changes
  ),
  "\n"
)

if (
  length(
    frozen_changes
  ) > 0
) {

  print(
    frozen_changes
  )

  stop(
    "STEP22A_FROZEN_PROJECT_CHANGED"
  )
}

cat(
  "✓ Frozen Step19-21 inventory unchanged.\n\n"
)

# ============================================================
# 16. Session info
# ============================================================

sink(
  file.path(
    LOG_DIR,
    "STEP22A_06D_sessionInfo.txt"
  )
)

print(
  sessionInfo()
)

sink()

# ============================================================
# 17. Summary
# ============================================================

summary_lines <- c(
  "STEP22A-6D ALL27 FIBROBLAST RAW-COUNT PSEUDOBULK",
  paste0(
    "Frozen fibroblast cells: ",
    nrow(fibro_ids)
  ),
  paste0(
    "Full common ENSG features: ",
    nrow(common)
  ),
  paste0(
    "Sample-level pseudobulks: ",
    ncol(sample_pb)
  ),
  paste0(
    "Patient-timepoint pseudobulks: ",
    ncol(pt_pb)
  ),
  paste0(
    "Technical PASS >=20: ",
    technical_pass_n,
    "/",
    nrow(pt_meta)
  ),
  paste0(
    "Technical FAIL <20: ",
    nrow(
      technical_fail
    )
  ),
  paste0(
    "Complete paired patients: ",
    n_paired
  ),
  paste0(
    "Technically evaluable paired patients: ",
    n_paired_technical,
    "/",
    n_paired
  ),
  "CORE rows restored structurally: 4/4",
  "CORE gene-specific values printed: NO",
  "Normalization performed: NO",
  "CPM calculated: NO",
  "log2(CPM+1) calculated: NO",
  "z-score calculated: NO",
  "CORE2 calculated: NO",
  "CORE4 calculated: NO",
  "Before/After inference: NO",
  "Histology mapping: INCOMPLETE",
  "EASC identity: UNRESOLVED",
  "Final ESCC reference: BLOCKED",
  "Primary CORE2 inference: BLOCKED",
  "Frozen Step19-21 changed: NO",
  "Manuscript changed: NO"
)

writeLines(
  summary_lines,
  SUMMARY_OUT
)

# ============================================================
# FINAL
# ============================================================

cat("\n")
cat("====================================================================\n")
cat("STEP22A-6D COMPLETE\n")
cat("====================================================================\n\n")

cat(
  "FROZEN FIBROBLAST CELLS       :",
  nrow(fibro_ids),
  "\n"
)

cat(
  "FULL COMMON ENSG FEATURES     :",
  nrow(common),
  "\n"
)

cat(
  "CORE ROWS RESTORED            : 4 / 4\n\n"
)

cat(
  "SAMPLE PSEUDOBULKS            :",
  ncol(sample_pb),
  "\n"
)

cat(
  "PATIENT-TIMEPOINT PSEUDOBULKS :",
  ncol(pt_pb),
  "\n\n"
)

cat(
  "TECHNICAL PASS >=20           :",
  technical_pass_n,
  "/",
  nrow(pt_meta),
  "\n"
)

cat(
  "TECHNICAL FAIL <20            :",
  nrow(
    technical_fail
  ),
  "\n"
)

if (
  nrow(
    technical_fail
  ) > 0
) {

  cat(
    "FAILED PATIENT-TIMEPOINT(S)  :"
  )

  cat(
    paste0(
      " ",
      technical_fail$patient_timepoint,
      " (",
      technical_fail$n_fibroblast_cells,
      " cells)"
    ),
    sep = ""
  )

  cat("\n")
}

cat("\n")

cat(
  "COMPLETE PAIRED PATIENTS      :",
  n_paired,
  "\n"
)

cat(
  "TECHNICALLY EVALUABLE PAIRS   :",
  n_paired_technical,
  "/",
  n_paired,
  "\n"
)

paired_fail <- paired_only[
  technically_evaluable_pair ==
    FALSE
]

if (
  nrow(
    paired_fail
  ) > 0
) {

  cat(
    "PAIRED TECHNICAL FAILURES    :",
    paste(
      paired_fail$patient_id,
      collapse = ", "
    ),
    "\n"
  )
}

cat("\n")

cat(
  "RAW INTEGER COUNTS            : YES\n"
)

cat(
  "NORMALIZATION                 : NOT PERFORMED\n"
)

cat(
  "CPM                           : NOT CALCULATED\n"
)

cat(
  "log2(CPM+1)                   : NOT CALCULATED\n"
)

cat(
  "Z-SCORE                       : NOT CALCULATED\n\n"
)

cat(
  "CORE GENE-SPECIFIC VALUES     : NOT PRINTED / NOT INSPECTED\n"
)

cat(
  "CORE2 SCORING                 : NOT PERFORMED\n"
)

cat(
  "CORE4 SCORING                 : NOT PERFORMED\n"
)

cat(
  "BEFORE/AFTER INFERENCE        : NOT PERFORMED\n\n"
)

cat(
  "HISTOLOGY MAPPING             : INCOMPLETE\n"
)

cat(
  "EASC IDENTITY                 : UNRESOLVED\n"
)

cat(
  "FINAL ESCC REFERENCE          : BLOCKED\n"
)

cat(
  "PRIMARY CORE2 INFERENCE       : BLOCKED\n\n"
)

cat(
  "PSEUDOBULK BUNDLE:\n"
)

cat(
  "  03_objects/OMIX005710/\n"
)

cat(
  "  OMIX005710_ALL27_FIBROBLAST_RAW_PSEUDOBULK_",
  "HISTOLOGY_UNRESOLVED.rds\n\n",
  sep = ""
)

cat(
  "FROZEN STEP19-21 MODIFIED     : NO\n"
)

cat(
  "MANUSCRIPT MODIFIED           : NO\n\n"
)

cat(
  "NEXT REQUIRED STEP:\n"
)

cat(
  "STEP22A-6E = AUDIT THE FINAL TECHNICAL REFERENCE\n"
)

cat(
  "AND FREEZE THE CONSEQUENCE OF ANY <20-CELL SAMPLE\n"
)

cat(
  "BEFORE ANY NORMALIZATION OR CORE2 CALCULATION.\n\n"
)

cat(
  "STOP HERE.\n"
)

cat(
  "DO NOT CALCULATE CORE2.\n\n"
)
