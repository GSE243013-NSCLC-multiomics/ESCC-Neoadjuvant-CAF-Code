# ============================================================
# ESCC Neoadjuvant Project
# STEP22A-5C
#
# FIBROBLAST PSEUDOBULK TECHNICAL PREPARATION
#
# Input:
#   - Step22A-5A v2 Seurat object
#   - Step22A-5B frozen broad-cell annotation
#
# This step:
#   1. selects cells ALREADY classified as Fibroblast_CAF
#   2. aggregates RAW counts to sample level
#   3. collapses to patient × timepoint level
#   4. performs technical integrity/QC only
#
# This step DOES NOT:
#   - normalize pseudobulk counts
#   - compute CORE2
#   - compute CORE4
#   - inspect treatment-associated expression changes
#   - perform DE
#   - perform paired statistics
#   - exclude any patient based on histology
#   - modify manuscript
#
# Histology mapping remains INCOMPLETE.
# ============================================================

options(
  stringsAsFactors = FALSE,
  timeout = 1200,
  future.globals.maxSize = 8 * 1024^3
)

set.seed(20260814)

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(Seurat)
})

cat("\n")
cat("====================================================================\n")
cat("STEP22A-5C: FIBROBLAST PSEUDOBULK TECHNICAL PREPARATION\n")
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

INPUT_OBJECT <- file.path(
  OBJECT_DIR,
  "OMIX005710_PAIRED_QC_CLUSTERED_PREANNOTATION_v2_COMMON_ENSG.rds"
)

ANNOTATION_FILE <- file.path(
  META_DIR,
  "STEP22A_05B_ANNOTATED_CELL_METADATA.csv"
)

FEATURE_MAP_FILE <- file.path(
  META_DIR,
  "STEP22A_05A_COMMON_ENSG_FEATURES.csv"
)

GATE_FILE <- file.path(
  META_DIR,
  "STEP22A_ANALYSIS_GATE_STATUS.csv"
)

FROZEN_INVENTORY <- file.path(
  META_DIR,
  "STEP22A_PREANALYSIS_FROZEN_FILE_INVENTORY.csv"
)

SAMPLE_OBJECT_OUT <- file.path(
  OBJECT_DIR,
  "OMIX005710_FIBROBLAST_SAMPLE_PSEUDOBULK_HISTOLOGY_UNRESOLVED.rds"
)

PATIENT_OBJECT_OUT <- file.path(
  OBJECT_DIR,
  "OMIX005710_FIBROBLAST_PATIENT_TIMEPOINT_PSEUDOBULK_HISTOLOGY_UNRESOLVED.rds"
)

# ============================================================
# 1. Gate check
# ============================================================

cat("===== 1. ANALYSIS GATE CHECK =====\n\n")

if (!file.exists(GATE_FILE)) {
  stop(
    "Missing gate file: ",
    GATE_FILE
  )
}

gate_dt <- fread(
  GATE_FILE,
  encoding = "UTF-8"
)

get_gate <- function(x) {

  z <- gate_dt[
    gate == x
  ]

  if (nrow(z) != 1) {
    stop(
      "Cannot uniquely resolve gate: ",
      x
    )
  }

  as.character(
    z$allowed[1]
  )
}

technical_allowed <- get_gate(
  "technical_preprocessing"
)

fibro_qc_allowed <- get_gate(
  "fibroblast_QC"
)

core2_allowed <- get_gate(
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
  "Fibroblast QC           :",
  fibro_qc_allowed,
  "\n"
)

cat(
  "Paired CORE2 inference  :",
  core2_allowed,
  "\n"
)

cat(
  "Manuscript update       :",
  manuscript_allowed,
  "\n\n"
)

if (technical_allowed != "YES") {
  stop(
    "TECHNICAL_PREPROCESSING_BLOCKED"
  )
}

if (fibro_qc_allowed != "YES") {
  stop(
    "FIBROBLAST_QC_BLOCKED"
  )
}

if (core2_allowed == "YES") {
  stop(
    "ERROR: CORE2 inference must remain blocked."
  )
}

if (manuscript_allowed == "YES") {
  stop(
    "ERROR: manuscript modification must remain blocked."
  )
}

cat(
  "✓ Technical pseudobulk preparation allowed.\n"
)

cat(
  "✓ CORE2 inference remains blocked.\n\n"
)

# ============================================================
# 2. Load object + frozen annotation metadata
# ============================================================

cat("===== 2. LOAD INPUT OBJECT + ANNOTATION =====\n\n")

for (f in c(
  INPUT_OBJECT,
  ANNOTATION_FILE,
  FEATURE_MAP_FILE
)) {

  if (!file.exists(f)) {
    stop(
      "Missing required file:\n",
      f
    )
  }
}

obj <- readRDS(
  INPUT_OBJECT
)

DefaultAssay(
  obj
) <- "RNA"

ann <- fread(
  ANNOTATION_FILE,
  encoding = "UTF-8"
)

cat(
  "Seurat cells     :",
  ncol(obj),
  "\n"
)

cat(
  "Annotation rows  :",
  nrow(ann),
  "\n"
)

if (anyDuplicated(
  ann$cell_id
) > 0) {

  stop(
    "DUPLICATE_CELL_IDS_IN_ANNOTATION"
  )
}

missing_from_annotation <- setdiff(
  colnames(obj),
  ann$cell_id
)

extra_annotation <- setdiff(
  ann$cell_id,
  colnames(obj)
)

cat(
  "Object cells missing annotation :",
  length(
    missing_from_annotation
  ),
  "\n"
)

cat(
  "Annotation cells absent object  :",
  length(
    extra_annotation
  ),
  "\n\n"
)

if (
  length(
    missing_from_annotation
  ) > 0 ||
  length(
    extra_annotation
  ) > 0
) {

  stop(
    "OBJECT_ANNOTATION_CELL_SET_MISMATCH"
  )
}

ann <- ann[
  match(
    colnames(obj),
    cell_id
  )
]

if (!identical(
  ann$cell_id,
  colnames(obj)
)) {

  stop(
    "OBJECT_ANNOTATION_ORDER_MISMATCH"
  )
}

# ============================================================
# 3. Select previously frozen Fibroblast_CAF cells
#
# IMPORTANT:
# No gene expression is used here for selection.
# Selection comes ONLY from Step22A-5B annotation.
# ============================================================

cat("===== 3. SELECT FROZEN Fibroblast_CAF COMPARTMENT =====\n\n")

fib_mask <- (
  ann$broad_celltype ==
    "Fibroblast_CAF"
)

fib_cells <- ann$cell_id[
  fib_mask
]

fib_meta <- copy(
  ann[
    fib_mask
  ]
)

cat(
  "Fibroblast_CAF cells:",
  length(
    fib_cells
  ),
  "\n"
)

if (
  length(
    fib_cells
  ) < 100
) {

  stop(
    "TOO_FEW_FIBROBLAST_CAF_CELLS"
  )
}

if (
  length(
    fib_cells
  ) != 13270
) {

  cat(
    "WARNING: Step22A-5B screenshot reported 13,270 fibroblast cells; ",
    "current frozen annotation contains ",
    length(fib_cells),
    ".\n",
    sep = ""
  )
}

cat(
  "Selection source: Step22A-5B broad annotation only.\n"
)

cat(
  "CORE genes used for selection: NO\n\n"
)

# ============================================================
# 4. Retrieve RAW count matrix
# ============================================================

cat("===== 4. RETRIEVE RAW COUNTS =====\n\n")

raw_counts <- tryCatch(

  {
    LayerData(
      obj,
      assay = "RNA",
      layer = "counts"
    )
  },

  error = function(e) {

    GetAssayData(
      obj,
      assay = "RNA",
      slot = "counts"
    )
  }
)

if (
  nrow(raw_counts) == 0 ||
  ncol(raw_counts) == 0
) {

  stop(
    "RAW_COUNT_LAYER_EMPTY"
  )
}

if (!identical(
  colnames(raw_counts),
  colnames(obj)
)) {

  stop(
    "RAW_COUNT_CELL_ORDER_MISMATCH"
  )
}

fib_counts <- raw_counts[
  ,
  fib_cells,
  drop = FALSE
]

cat(
  "Fibroblast raw-count matrix:",
  nrow(fib_counts),
  "features x",
  ncol(fib_counts),
  "cells\n\n"
)

# Release the full count matrix reference if possible.
rm(
  raw_counts
)

gc()

# ============================================================
# 5. Verify feature identity
# ============================================================

cat("===== 5. FEATURE IDENTITY CHECK =====\n\n")

feature_map <- fread(
  FEATURE_MAP_FILE,
  encoding = "UTF-8"
)

# Seurat replaces '_' with '-' in feature names at object creation.
# Reconcile the frozen feature map to the object-side names.
object_feature_names <- rownames(
  fib_counts
)

reconciled_feature_names <- gsub(
  "_",
  "-",
  feature_map$seurat_feature_name
)

if (!identical(
  reconciled_feature_names,
  object_feature_names
)) {

  stop(
    "PSEUDOBULK_FEATURE_NAME_RECONCILIATION_FAILED"
  )
}

feature_meta <- data.table(
  common_feature_index =
    feature_map$common_feature_index,

  gene_id =
    feature_map$gene_id,

  gene_symbol =
    feature_map$gene_symbol,

  seurat_feature_name =
    object_feature_names
)

if (!identical(
  feature_meta$seurat_feature_name,
  rownames(
    fib_counts
  )
)) {

  stop(
    "PSEUDOBULK_FEATURE_ORDER_MISMATCH"
  )
}

cat(
  "Common aligned features:",
  nrow(
    feature_meta
  ),
  "\n"
)

cat(
  "Alignment key: ENSG ID\n"
)

cat(
  "Union + zero filling: NO\n\n"
)

# ============================================================
# 6. Build SAMPLE-level fibroblast pseudobulk
#
# Raw integer counts are summed.
# ============================================================

cat("===== 6. SAMPLE-LEVEL FIBROBLAST PSEUDOBULK =====\n\n")

required_meta <- c(
  "patient_id",
  "timepoint",
  "sample_name",
  "hrs_accession"
)

missing_meta <- setdiff(
  required_meta,
  names(
    fib_meta
  )
)

if (
  length(
    missing_meta
  ) > 0
) {

  stop(
    "Missing fibroblast metadata columns: ",
    paste(
      missing_meta,
      collapse = ", "
    )
  )
}

sample_check <- fib_meta[
  ,
  .(
    patient_n =
      uniqueN(
        patient_id
      ),

    timepoint_n =
      uniqueN(
        timepoint
      ),

    hrs_n =
      uniqueN(
        hrs_accession
      ),

    n_fibroblasts =
      .N
  ),
  by = sample_name
]

if (
  any(
    sample_check$patient_n != 1
  ) ||
  any(
    sample_check$timepoint_n != 1
  ) ||
  any(
    sample_check$hrs_n != 1
  )
) {

  print(
    sample_check
  )

  stop(
    "SAMPLE_METADATA_NOT_UNIQUE"
  )
}

sample_levels <- unique(
  fib_meta$sample_name
)

sample_design <- sparseMatrix(

  i = seq_len(
    nrow(
      fib_meta
    )
  ),

  j = match(
    fib_meta$sample_name,
    sample_levels
  ),

  x = 1,

  dims = c(
    nrow(
      fib_meta
    ),
    length(
      sample_levels
    )
  ),

  dimnames = list(
    fib_meta$cell_id,
    sample_levels
  )
)

if (!identical(
  rownames(
    sample_design
  ),
  colnames(
    fib_counts
  )
)) {

  stop(
    "SAMPLE_AGGREGATION_DESIGN_ORDER_MISMATCH"
  )
}

sample_pb <- fib_counts %*%
  sample_design

sample_pb <- as(
  sample_pb,
  "dgCMatrix"
)

sample_meta <- fib_meta[
  ,
  .(
    patient_id =
      unique(
        patient_id
      ),

    timepoint =
      unique(
        timepoint
      ),

    hrs_accession =
      unique(
        hrs_accession
      ),

    n_fibroblast_cells =
      .N
  ),
  by = sample_name
]

sample_meta <- sample_meta[
  match(
    colnames(
      sample_pb
    ),
    sample_name
  )
]

sample_meta[
  ,
  library_size :=
    as.numeric(
      Matrix::colSums(
        sample_pb
      )
    )
]

sample_meta[
  ,
  detected_genes :=
    as.integer(
      Matrix::colSums(
        sample_pb > 0
      )
    )
]

sample_meta[
  ,
  technical_status :=
    fifelse(
      n_fibroblast_cells >= 20,
      "PASS_GE20",
      "INSUFFICIENT_FIBROBLAST_CELLS"
    )
]

cat(
  "Sample-level pseudobulks:",
  ncol(
    sample_pb
  ),
  "\n\n"
)

print(
  sample_meta
)

# ============================================================
# 7. Collapse sample-level pseudobulk to PATIENT × TIMEPOINT
#
# This follows the frozen rule:
# sample-level first, then collapse technical/specimen-level
# units belonging to the same patient-timepoint.
# ============================================================

cat("\n===== 7. PATIENT × TIMEPOINT COLLAPSE =====\n\n")

sample_meta[
  ,
  patient_timepoint :=
    paste(
      patient_id,
      timepoint,
      sep = "__"
    )
]

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

pt_order <- unique(
  sample_meta[
    order(
      patient_number(
        patient_id
      ),
      factor(
        timepoint,
        levels = c(
          "Before",
          "After"
        )
      )
    ),
    patient_timepoint
  ]
)

collapse_design <- sparseMatrix(

  i = seq_len(
    nrow(
      sample_meta
    )
  ),

  j = match(
    sample_meta$patient_timepoint,
    pt_order
  ),

  x = 1,

  dims = c(
    nrow(
      sample_meta
    ),
    length(
      pt_order
    )
  ),

  dimnames = list(
    sample_meta$sample_name,
    pt_order
  )
)

if (!identical(
  rownames(
    collapse_design
  ),
  colnames(
    sample_pb
  )
)) {

  stop(
    "PATIENT_TIMEPOINT_COLLAPSE_ORDER_MISMATCH"
  )
}

patient_pb <- sample_pb %*%
  collapse_design

patient_pb <- as(
  patient_pb,
  "dgCMatrix"
)

pt_meta <- sample_meta[
  ,
  .(
    n_fibroblast_cells =
      sum(
        n_fibroblast_cells
      ),

    number_of_source_specimens =
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
      )
  ),
  by = .(
    patient_id,
    timepoint,
    patient_timepoint
  )
]

pt_meta <- pt_meta[
  match(
    colnames(
      patient_pb
    ),
    patient_timepoint
  )
]

pt_meta[
  ,
  library_size :=
    as.numeric(
      Matrix::colSums(
        patient_pb
      )
    )
]

pt_meta[
  ,
  detected_genes :=
    as.integer(
      Matrix::colSums(
        patient_pb > 0
      )
    )
]

pt_meta[
  ,
  technical_status :=
    fifelse(
      n_fibroblast_cells >= 20,
      "PASS_GE20",
      "INSUFFICIENT_FIBROBLAST_CELLS"
    )
]

pt_meta[
  ,
  histology_mapping :=
    "INCOMPLETE"
]

pt_meta[
  ,
  primary_CORE2_inference_allowed :=
    "NO"
]

cat(
  "Patient-timepoint pseudobulks:",
  ncol(
    patient_pb
  ),
  "\n\n"
)

print(
  pt_meta
)

# ============================================================
# 8. Structural validation
# ============================================================

cat("\n===== 8. PSEUDOBULK STRUCTURAL VALIDATION =====\n\n")

n_patients <- uniqueN(
  pt_meta$patient_id
)

n_before <- sum(
  pt_meta$timepoint ==
    "Before"
)

n_after <- sum(
  pt_meta$timepoint ==
    "After"
)

paired_status <- pt_meta[
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
      )
  ),
  by = patient_id
]

paired_status[
  ,
  complete_pair :=
    has_before &
    has_after
]

paired_status[
  ,
  technical_pair_pass :=
    complete_pair &
    before_fibroblasts >= 20 &
    after_fibroblasts >= 20
]

paired_status[
  ,
  patient_order :=
    patient_number(
      patient_id
    )
]

setorder(
  paired_status,
  patient_order
)

paired_status[
  ,
  patient_order := NULL
]

cat(
  "Patients:",
  n_patients,
  "\n"
)

cat(
  "Patient-timepoints:",
  nrow(
    pt_meta
  ),
  "\n"
)

cat(
  "Before pseudobulks:",
  n_before,
  "\n"
)

cat(
  "After pseudobulks :",
  n_after,
  "\n"
)

cat(
  "Complete paired patients:",
  sum(
    paired_status$
      complete_pair
  ),
  "\n"
)

cat(
  "Technically passing pairs:",
  sum(
    paired_status$
      technical_pair_pass
  ),
  "\n\n"
)

print(
  paired_status
)

if (
  n_patients != 6
) {

  stop(
    "EXPECTED_6_PAIRED_PATIENTS"
  )
}

if (
  nrow(
    pt_meta
  ) != 12
) {

  stop(
    "EXPECTED_12_PATIENT_TIMEPOINT_PSEUDOBULKS"
  )
}

if (
  n_before != 6 ||
  n_after != 6
) {

  stop(
    "EXPECTED_6_BEFORE_AND_6_AFTER_PSEUDOBULKS"
  )
}

if (
  !all(
    paired_status$
      technical_pair_pass
  )
) {

  stop(
    "AT_LEAST_ONE_PAIRED_PATIENT_FAILS_FIBROBLAST_TECHNICAL_QC"
  )
}

if (
  any(
    Matrix::colSums(
      patient_pb
    ) <= 0
  )
) {

  stop(
    "ZERO_LIBRARY_SIZE_PSEUDOBULK"
  )
}

cat(
  "✓ 6/6 paired patients technically evaluable.\n"
)

cat(
  "✓ 12/12 patient-timepoint pseudobulks have non-zero libraries.\n\n"
)

# ============================================================
# 9. Explicitly DO NOT normalize / score
# ============================================================

cat("===== 9. INFERENCE FIREWALL =====\n\n")

cat(
  "Raw counts aggregated       : YES\n"
)

cat(
  "CPM normalization           : NO\n"
)

cat(
  "TMM normalization           : NO\n"
)

cat(
  "DESeq2 normalization        : NO\n"
)

cat(
  "logCPM calculation          : NO\n"
)

cat(
  "CORE2 scoring               : NO\n"
)

cat(
  "CORE4 scoring               : NO\n"
)

cat(
  "Before/After test           : NO\n"
)

cat(
  "Before/After gene inspection: NO\n"
)

cat(
  "Histology-based exclusion   : NO\n\n"
)

# ============================================================
# 10. Save sample-level metadata
# ============================================================

SAMPLE_META_OUT <- file.path(
  META_DIR,
  "STEP22A_05C_FIBROBLAST_SAMPLE_PSEUDOBULK_METADATA.csv"
)

fwrite(
  sample_meta,
  SAMPLE_META_OUT,
  bom = TRUE
)

# ============================================================
# 11. Save patient-timepoint metadata
# ============================================================

PT_META_OUT <- file.path(
  META_DIR,
  "STEP22A_05C_FIBROBLAST_PATIENT_TIMEPOINT_METADATA.csv"
)

fwrite(
  pt_meta,
  PT_META_OUT,
  bom = TRUE
)

PAIR_OUT <- file.path(
  META_DIR,
  "STEP22A_05C_FIBROBLAST_PAIRED_TECHNICAL_STATUS.csv"
)

fwrite(
  paired_status,
  PAIR_OUT,
  bom = TRUE
)

# ============================================================
# 12. Save feature metadata
# ============================================================

FEATURE_OUT <- file.path(
  META_DIR,
  "STEP22A_05C_PSEUDOBULK_FEATURE_METADATA.csv"
)

fwrite(
  feature_meta,
  FEATURE_OUT,
  bom = TRUE
)

# ============================================================
# 13. Save sample-level pseudobulk object
# ============================================================

sample_bundle <- list(

  counts =
    sample_pb,

  sample_metadata =
    as.data.frame(
      sample_meta
    ),

  feature_metadata =
    as.data.frame(
      feature_meta
    ),

  aggregation_level =
    "sample",

  count_type =
    "raw_integer_sum",

  cell_selection_source =
    "STEP22A-5B_FROZEN_Fibroblast_CAF",

  histology_mapping =
    "INCOMPLETE",

  CORE2_scored =
    FALSE,

  CORE4_scored =
    FALSE,

  inference_allowed =
    FALSE
)

saveRDS(
  sample_bundle,
  SAMPLE_OBJECT_OUT,
  compress = "gzip"
)

# ============================================================
# 14. Save patient-timepoint pseudobulk object
# ============================================================

patient_bundle <- list(

  counts =
    patient_pb,

  patient_timepoint_metadata =
    as.data.frame(
      pt_meta
    ),

  paired_technical_status =
    as.data.frame(
      paired_status
    ),

  feature_metadata =
    as.data.frame(
      feature_meta
    ),

  aggregation_level =
    "patient_timepoint",

  aggregation_path =
    "cells -> sample pseudobulk -> patient_timepoint",

  count_type =
    "raw_integer_sum",

  cell_selection_source =
    "STEP22A-5B_FROZEN_Fibroblast_CAF",

  minimum_fibroblast_cells =
    20L,

  histology_mapping =
    "INCOMPLETE",

  EASC_patient =
    "UNRESOLVED",

  CORE2_scored =
    FALSE,

  CORE4_scored =
    FALSE,

  treatment_inference_performed =
    FALSE,

  inference_allowed =
    FALSE
)

saveRDS(
  patient_bundle,
  PATIENT_OBJECT_OUT,
  compress = "gzip"
)

cat(
  "✓ Pseudobulk objects saved.\n\n"
)

# ============================================================
# 15. Small provenance manifest
# ============================================================

PROVENANCE_OUT <- file.path(
  META_DIR,
  "STEP22A_05C_PSEUDOBULK_PROVENANCE.csv"
)

provenance <- data.table(

  item = c(
    "source_Seurat_object",
    "source_annotation",
    "aggregation_path",
    "cell_type_selected",
    "cell_selection_uses_CORE_genes",
    "raw_counts_used",
    "sample_level_pseudobulk_created",
    "patient_timepoint_pseudobulk_created",
    "normalization_performed",
    "CORE2_scoring_performed",
    "CORE4_scoring_performed",
    "treatment_test_performed",
    "histology_mapping",
    "EASC_patient",
    "primary_inference_allowed"
  ),

  value = c(
    basename(
      INPUT_OBJECT
    ),
    basename(
      ANNOTATION_FILE
    ),
    "cells_to_sample_to_patient_timepoint",
    "Fibroblast_CAF",
    "NO",
    "YES",
    "YES",
    "YES",
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
# 16. Frozen Step19-21 check
# ============================================================

cat("===== 10. FROZEN PROJECT CHECK =====\n\n")

frozen_changes <- list()

if (
  file.exists(
    FROZEN_INVENTORY
  )
) {

  frozen <- fread(
    FROZEN_INVENTORY,
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
# 17. Session info
# ============================================================

sink(
  file.path(
    LOG_DIR,
    "STEP22A_05C_sessionInfo.txt"
  )
)

print(
  sessionInfo()
)

sink()

# ============================================================
# 18. Summary
# ============================================================

SUMMARY_OUT <- file.path(
  META_DIR,
  "STEP22A_05C_FIBROBLAST_PSEUDOBULK_TECHNICAL_SUMMARY.txt"
)

summary_lines <- c(

  "STEP22A-5C FIBROBLAST PSEUDOBULK TECHNICAL PREPARATION",

  paste0(
    "Fibroblast_CAF cells: ",
    length(
      fib_cells
    )
  ),

  paste0(
    "Features: ",
    nrow(
      patient_pb
    )
  ),

  paste0(
    "Sample-level pseudobulks: ",
    ncol(
      sample_pb
    )
  ),

  paste0(
    "Patient-timepoint pseudobulks: ",
    ncol(
      patient_pb
    )
  ),

  paste0(
    "Paired patients: ",
    n_patients
  ),

  paste0(
    "Technically passing paired patients: ",
    sum(
      paired_status$
        technical_pair_pass
    )
  ),

  "Raw counts only: YES",

  "Normalization performed: NO",

  "CORE2 scoring performed: NO",

  "CORE4 scoring performed: NO",

  "Treatment inference performed: NO",

  "Histology mapping: INCOMPLETE",

  "EASC patient: UNRESOLVED",

  "Primary CORE2 inference allowed: NO",

  "Frozen Step19-21 changed: NO",

  "Manuscript changed: NO"
)

writeLines(
  summary_lines,
  SUMMARY_OUT
)

# ============================================================
# FINAL PRINT
# ============================================================

cat("\n")
cat("====================================================================\n")
cat("STEP22A-5C COMPLETE\n")
cat("====================================================================\n\n")

cat(
  "FIBROBLAST_CAF CELLS         :",
  length(
    fib_cells
  ),
  "\n"
)

cat(
  "COMMON ENSG FEATURES         :",
  nrow(
    patient_pb
  ),
  "\n\n"
)

cat(
  "SAMPLE PSEUDOBULKS           :",
  ncol(
    sample_pb
  ),
  "\n"
)

cat(
  "PATIENT-TIMEPOINT PSEUDOBULKS:",
  ncol(
    patient_pb
  ),
  "\n"
)

cat(
  "PAIRED PATIENTS              :",
  n_patients,
  "\n"
)

cat(
  "TECHNICALLY PASSING PAIRS    :",
  sum(
    paired_status$
      technical_pair_pass
  ),
  "/",
  n_patients,
  "\n\n"
)

cat(
  "RAW COUNTS AGGREGATED        : YES\n"
)

cat(
  "NORMALIZATION PERFORMED      : NO\n"
)

cat(
  "CORE2 SCORING                : NOT PERFORMED\n"
)

cat(
  "CORE4 SCORING                : NOT PERFORMED\n"
)

cat(
  "BEFORE/AFTER INFERENCE       : NOT PERFORMED\n\n"
)

cat(
  "HISTOLOGY MAPPING            : INCOMPLETE\n"
)

cat(
  "EASC PATIENT                 : UNRESOLVED\n"
)

cat(
  "PRIMARY CORE2 INFERENCE      : BLOCKED\n\n"
)

cat(
  "FROZEN STEP19-21 MODIFIED    : NO\n"
)

cat(
  "MANUSCRIPT MODIFIED          : NO\n\n"
)

cat(
  "PATIENT-TIMEPOINT OBJECT:\n"
)

cat(
  "  03_objects/OMIX005710/\n"
)

cat(
  "  OMIX005710_FIBROBLAST_PATIENT_TIMEPOINT_PSEUDOBULK_",
  "HISTOLOGY_UNRESOLVED.rds\n\n",
  sep = ""
)

cat(
  "NEXT REQUIRED STEP:\n"
)

cat(
  "STEP22A-5D = PSEUDOBULK TECHNICAL QC / ",
  "FROZEN SCORING IMPLEMENTATION AUDIT\n"
)

cat(
  "CORE2 VALUES MUST STILL NOT BE CALCULATED.\n\n"
)

cat(
  "STOP HERE.\n"
)

cat(
  "DO NOT RUN BEFORE/AFTER STATISTICS.\n\n"
)
