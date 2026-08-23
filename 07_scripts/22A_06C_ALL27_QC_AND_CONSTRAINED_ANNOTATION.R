# ============================================================
# ESCC Neoadjuvant Project
# STEP22A-6C
#
# ALL 27 TUMORS:
#   - dynamic ENSG alignment to true all-27 common universe
#   - same frozen sample-specific technical QC
#   - joint PCA / graph clustering
#   - constrained broad cell-type annotation
#
# IMPORTANT ANTI-CIRCULARITY RULE:
#
# CDKN1B / GSN / BGN / TIMP1 are present in raw/common data,
# but are REMOVED BEFORE:
#   - QC calculation
#   - normalization
#   - HVG selection
#   - PCA
#   - clustering
#   - annotation
#
# Thus CORE genes cannot influence fibroblast-cell selection.
#
# Frozen biological marker dictionary remains 33.
# Joint measurable annotation markers = 31:
#   KRT5 / KRT14 structurally unavailable in all-27 common set.
#
# NO CORE scoring.
# NO Before/After inference.
# Histology remains unresolved.
# ============================================================

options(
  stringsAsFactors = FALSE,
  timeout = 1800,
  future.globals.maxSize = 16 * 1024^3
)

# 16 GiB physical RAM machine: reduce R vector heap ceiling to leave
# headroom for the OS and prevent silent OOM during ScaleData/PCA.
mem.maxVSize(20 * 1024^3)

set.seed(20260814)

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(Seurat)
  library(ggplot2)
})

cat("\n")
cat("====================================================================\n")
cat("STEP22A-6C: ALL27 QC + CONSTRAINED BROAD ANNOTATION\n")
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

FIGURE_DIR <- file.path(
  PROJECT,
  "05_figures",
  "external_validation",
  "Step22A"
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
  FIGURE_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  LOG_DIR,
  recursive = TRUE,
  showWarnings = FALSE
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

COMMON_FILE <- file.path(
  META_DIR,
  "STEP22A_06B_ALL27_COMMON_ENSG_FEATURES.csv"
)

MARKER_LOCK_FILE <- file.path(
  META_DIR,
  "STEP22A_06B1_JOINT_ANNOTATION_MARKER_LOCK.csv"
)

MISSINGNESS_STATUS <- file.path(
  META_DIR,
  "STEP22A_06B1_ALIGNMENT_FINAL_STATUS.csv"
)

GATE_FILE <- file.path(
  META_DIR,
  "STEP22A_ANALYSIS_GATE_STATUS.csv"
)

FROZEN_INV <- file.path(
  META_DIR,
  "STEP22A_PREANALYSIS_FROZEN_FILE_INVENTORY.csv"
)

QC_THRESHOLDS_OUT <- file.path(
  META_DIR,
  "STEP22A_06C_ALL27_SAMPLE_QC_THRESHOLDS.csv"
)

RETAINED_OUT <- file.path(
  META_DIR,
  "STEP22A_06C_ALL27_RETAINED_CELL_COUNTS.csv"
)

CLUSTER_EVIDENCE_OUT <- file.path(
  META_DIR,
  "STEP22A_06C_ALL27_CLUSTER_MODULE_EVIDENCE.csv"
)

CLUSTER_ANNOTATION_OUT <- file.path(
  META_DIR,
  "STEP22A_06C_ALL27_CLUSTER_BROAD_ANNOTATION.csv"
)

BROAD_COUNTS_OUT <- file.path(
  META_DIR,
  "STEP22A_06C_ALL27_BROAD_CELLTYPE_COUNTS.csv"
)

CELL_METADATA_OUT <- file.path(
  META_DIR,
  "STEP22A_06C_ALL27_ANNOTATED_CELL_METADATA.csv"
)

FIBRO_IDS_OUT <- file.path(
  META_DIR,
  "STEP22A_06C_ALL27_FIBROBLAST_CELL_IDS.csv"
)

FIBRO_QC_OUT <- file.path(
  META_DIR,
  "STEP22A_06C_ALL27_FIBROBLAST_SAMPLE_QC.csv"
)

MARKER_MAP_OUT <- file.path(
  META_DIR,
  "STEP22A_06C_ALL27_MARKER_TO_ENSG_MAPPING.csv"
)

ANTI_CIRC_OUT <- file.path(
  META_DIR,
  "STEP22A_06C_ANTI_CIRCULARITY_AUDIT.csv"
)

SUMMARY_OUT <- file.path(
  META_DIR,
  "STEP22A_06C_ALL27_QC_ANNOTATION_SUMMARY.txt"
)

BUNDLE_OUT <- file.path(
  OBJECT_DIR,
  "OMIX005710_STEP22A_06C_ALL27_ANNOTATION_BUNDLE.rds"
)

# ============================================================
# Helpers
# ============================================================

robust_scale <- function(x) {

  z <- mad(
    x,
    constant = 1.4826,
    na.rm = TRUE
  )

  if (
    !is.finite(z) ||
    z <= 0
  ) {

    z <- IQR(
      x,
      na.rm = TRUE
    ) / 1.349
  }

  if (
    !is.finite(z) ||
    z <= 0
  ) {

    z <- sd(
      x,
      na.rm = TRUE
    )
  }

  if (
    !is.finite(z) ||
    z <= 0
  ) {

    z <- 1
  }

  z
}


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

  if (anyDuplicated(gene_id) > 0) {

    stop(
      "Duplicate ENSG IDs in: ",
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
      "Cannot uniquely resolve gate: ",
      x
    )
  }

  as.character(z)
}

technical_allowed <- get_gate(
  "technical_preprocessing"
)

annotation_allowed <- get_gate(
  "broad_cell_annotation"
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
  "Broad annotation        :",
  annotation_allowed,
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

if (technical_allowed != "YES") {
  stop("TECHNICAL_PREPROCESSING_BLOCKED")
}

if (annotation_allowed != "YES") {
  stop("BROAD_ANNOTATION_BLOCKED")
}

if (core_allowed == "YES") {
  stop("CORE2 inference gate unexpectedly open")
}

if (manuscript_allowed == "YES") {
  stop("Manuscript gate unexpectedly open")
}

cat(
  "✓ QC / broad annotation allowed.\n"
)

cat(
  "✓ CORE2 inference remains blocked.\n\n"
)

# ============================================================
# 2. Verify Step22A-6B1 decision
# ============================================================

cat("===== 2. VERIFY STRUCTURAL-MISSINGNESS LOCK =====\n\n")

status_6b1 <- fread(
  MISSINGNESS_STATUS,
  encoding = "UTF-8"
)

final_status <- status_6b1[
  item ==
    "STEP22A_06B_final_status",
  value
]

if (
  length(final_status) != 1 ||
  final_status !=
    "PASS_WITH_FROZEN_STRUCTURAL_MARKER_MISSINGNESS"
) {

  stop(
    "STEP22A_06B1_STATUS_NOT_APPROVED"
  )
}

cat(
  "Step22A-6B status:",
  final_status,
  "\n\n"
)

# ============================================================
# 3. Load 27-sample inventory + matrix paths
# ============================================================

cat("===== 3. LOAD ALL 27 MATRIX PATHS =====\n\n")

combined <- fread(
  COMBINED_INV,
  encoding = "UTF-8"
)

if (nrow(combined) != 27) {
  stop("Expected 27 tumor observations")
}

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
  nrow(audit) != 27 ||
  anyDuplicated(
    audit$sample_name
  ) > 0
) {

  stop(
    "ALL27_MATRIX_AUDIT_INVALID"
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
  "Unique patients:",
  uniqueN(
    sample_table$patient_id
  ),
  "\n"
)

cat(
  "Before:",
  sum(
    sample_table$timepoint ==
      "Before"
  ),
  "\n"
)

cat(
  "After :",
  sum(
    sample_table$timepoint ==
      "After"
  ),
  "\n\n"
)

if (
  uniqueN(
    sample_table$patient_id
  ) != 21 ||
  sum(
    sample_table$timepoint ==
      "Before"
  ) != 15 ||
  sum(
    sample_table$timepoint ==
      "After"
  ) != 12
) {

  stop(
    "ALL27_COHORT_STRUCTURE_MISMATCH"
  )
}

# ============================================================
# 4. Load true all-27 common ENSG universe
# ============================================================

cat("===== 4. TRUE ALL-27 COMMON FEATURE UNIVERSE =====\n\n")

common <- fread(
  COMMON_FILE,
  encoding = "UTF-8"
)

required_common <- c(
  "all27_common_index",
  "gene_id",
  "gene_symbol"
)

if (
  !all(
    required_common %in%
      names(common)
  )
) {

  stop(
    "ALL27_COMMON_FILE_COLUMNS_INVALID"
  )
}

setorder(
  common,
  all27_common_index
)

if (
  nrow(common) != 16790
) {

  cat(
    "WARNING: screenshot reported 16,790 all-27 genes; file contains ",
    nrow(common),
    ".\n",
    sep = ""
  )
}

if (
  anyDuplicated(
    common$gene_id
  ) > 0
) {

  stop(
    "DUPLICATED_ALL27_COMMON_ENSG"
  )
}

common[
  ,
  symbol_upper :=
    toupper(
      gene_symbol
    )
]

core_genes <- c(
  "CDKN1B",
  "GSN",
  "BGN",
  "TIMP1"
)

core_map <- common[
  symbol_upper %in%
    core_genes
]

if (
  nrow(core_map) != 4 ||
  uniqueN(
    core_map$symbol_upper
  ) != 4
) {

  stop(
    "CORE_GENE_TO_ENSG_MAPPING_NOT_UNIQUE"
  )
}

core_ids <- core_map$gene_id

analysis_common <- common[
  !gene_id %in%
    core_ids
]

analysis_ids <- analysis_common$gene_id

cat(
  "All-27 common ENSG     :",
  nrow(common),
  "\n"
)

cat(
  "CORE genes withheld    :",
  length(core_ids),
  "\n"
)

cat(
  "Cell-analysis features :",
  length(analysis_ids),
  "\n\n"
)

cat(
  "CORE genes removed before QC/PCA/clustering/annotation:\n"
)

for (i in seq_len(nrow(core_map))) {

  cat(
    "  ",
    core_map$symbol_upper[i],
    " -> ",
    core_map$gene_id[i],
    "\n",
    sep = ""
  )
}

cat("\n")

# Mito gene IDs in non-CORE analysis universe.
mito_ids <- analysis_common[
  grepl(
    "^MT-",
    gene_symbol,
    ignore.case = TRUE
  ),
  gene_id
]

cat(
  "Mitochondrial genes in analysis universe:",
  length(mito_ids),
  "\n\n"
)

# ============================================================
# 5. Load 31-marker frozen joint annotation panel
# ============================================================

cat("===== 5. LOAD FROZEN JOINT MARKER PANEL =====\n\n")

marker_lock <- fread(
  MARKER_LOCK_FILE,
  encoding = "UTF-8"
)

marker_use <- marker_lock[
  marker_class ==
    "BROAD_ANNOTATION" &
  joint_annotation_use ==
    "YES"
]

if (
  nrow(marker_use) != 31
) {

  stop(
    "Expected exactly 31 measurable joint broad markers"
  )
}

if (
  any(
    marker_use$gene %in%
      core_genes
  )
) {

  stop(
    "CORE_GENE_FOUND_IN_ANNOTATION_PANEL"
  )
}

if (
  any(
    marker_use$gene %in%
      c(
        "KRT5",
        "KRT14"
      )
  )
) {

  stop(
    "STRUCTURALLY_UNAVAILABLE_KERATIN_IN_JOINT_PANEL"
  )
}

marker_map_list <- list()

for (
  i in seq_len(
    nrow(marker_use)
  )
) {

  gene <- marker_use$gene[i]

  module <- marker_use$marker_module[i]

  hits <- common[
    symbol_upper ==
      toupper(gene)
  ]

  if (nrow(hits) == 0) {

    stop(
      "Frozen measurable marker absent from all-27 universe: ",
      gene
    )
  }

  # Deterministic first common-order ENSG if a symbol maps more than once.
  selected <- hits[
    order(
      all27_common_index
    )
  ][1]

  if (
    selected$gene_id %in%
      core_ids
  ) {

    stop(
      "CORE feature unexpectedly selected as broad marker"
    )
  }

  marker_map_list[[
    length(marker_map_list) + 1
  ]] <- data.table(
    marker_module = module,
    marker = gene,
    gene_id = selected$gene_id,
    candidate_ENSG_count = nrow(hits),
    mapping_rule = "FIRST_ALL27_COMMON_ORDER"
  )
}

marker_map <- rbindlist(
  marker_map_list
)

if (nrow(marker_map) != 31) {
  stop("MARKER_MAPPING_INCOMPLETE")
}

fwrite(
  marker_map,
  MARKER_MAP_OUT,
  bom = TRUE
)

module_sizes <- marker_map[
  ,
  .N,
  by = marker_module
]

print(
  module_sizes
)

expected_module_sizes <- c(
  Fibroblast_CAF = 12,
  Epithelial = 4,
  Immune = 7,
  Endothelial = 4,
  Mural_Pericyte = 4
)

for (module in names(expected_module_sizes)) {

  n <- module_sizes[
    marker_module == module,
    N
  ]

  if (
    length(n) != 1 ||
    n != expected_module_sizes[
      module
    ]
  ) {

    stop(
      "Unexpected marker count for ",
      module
    )
  }
}

cat(
  "\n✓ Joint measurable panel = 31 markers.\n"
)

cat(
  "✓ KRT5/KRT14 remain structural NA, not zero.\n"
)

cat(
  "✓ CORE genes absent from annotation panel.\n\n"
)

# ============================================================
# 6. Per-sample QC BEFORE merging
#
# Same frozen algorithm as Step22A-5A:
#
#   nFeature >= max(200, median - 5 MAD)
#   nFeature <= median + 5 MAD
#
#   nCount >= max(500, median - 5 MAD)
#   nCount <= median + 5 MAD
#
#   percent.mt <= max(20, median + 5 MAD), capped at 40
#
# Treatment/timepoint is NOT used to determine thresholds.
#
# Difference from old exploratory 12-sample run:
# CORE genes are withheld uniformly from ALL 27 samples to
# eliminate any possible circular cell-selection pathway.
# ============================================================

cat("===== 6. ALL27 SAMPLE-SPECIFIC TECHNICAL QC =====\n\n")

matrix_list <- list()
meta_list <- list()
qc_threshold_list <- list()

total_raw_cells <- 0L
total_retained_cells <- 0L

common_ids <- common$gene_id

analysis_index_in_common <- match(
  analysis_ids,
  common_ids
)

mito_index_in_analysis <- match(
  mito_ids,
  analysis_ids
)

mito_index_in_analysis <- mito_index_in_analysis[
  !is.na(
    mito_index_in_analysis
  )
]

for (
  i in seq_len(
    nrow(sample_table)
  )
) {

  row <- sample_table[i]

  sample_name <- row$sample_name
  patient_id <- row$patient_id
  timepoint <- row$timepoint
  hrs <- row$hrs_accession

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
      "Duplicate barcodes: ",
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
      nrow(gt)
  ) {

    stop(
      "Gene dimension mismatch: ",
      sample_name
    )
  }

  if (
    ncol(counts_raw) !=
      length(barcodes)
  ) {

    stop(
      "Barcode dimension mismatch: ",
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
      "ALL27_COMMON_ENSG_MISSING_IN_SAMPLE: ",
      sample_name
    )
  }

  # Exact true all-27 common order.
  counts_common <- counts_raw[
    common_idx,
    ,
    drop = FALSE
  ]

  rm(
    counts_raw
  )

  # Explicitly withhold four CORE rows from cell analysis.
  counts <- counts_common[
    analysis_index_in_common,
    ,
    drop = FALSE
  ]

  rm(
    counts_common
  )

  rownames(
    counts
  ) <- analysis_ids

  cell_ids <- paste0(
    sample_name,
    "__",
    barcodes
  )

  colnames(
    counts
  ) <- cell_ids

  nCount <- as.numeric(
    Matrix::colSums(
      counts
    )
  )

  nFeature <- as.integer(
    Matrix::colSums(
      counts > 0
    )
  )

  if (
    length(
      mito_index_in_analysis
    ) > 0
  ) {

    mt_count <- as.numeric(
      Matrix::colSums(
        counts[
          mito_index_in_analysis,
          ,
          drop = FALSE
        ]
      )
    )

    percent_mt <- ifelse(
      nCount > 0,
      100 * mt_count /
        nCount,
      100
    )

  } else {

    percent_mt <- rep(
      0,
      length(
        nCount
      )
    )
  }

  feature_med <- median(
    nFeature
  )

  feature_mad <- robust_scale(
    nFeature
  )

  count_med <- median(
    nCount
  )

  count_mad <- robust_scale(
    nCount
  )

  mt_med <- median(
    percent_mt
  )

  mt_mad <- robust_scale(
    percent_mt
  )

  feature_low <- max(
    200,
    floor(
      feature_med -
        5 * feature_mad
    )
  )

  feature_high <- ceiling(
    feature_med +
      5 * feature_mad
  )

  count_low <- max(
    500,
    floor(
      count_med -
        5 * count_mad
    )
  )

  count_high <- ceiling(
    count_med +
      5 * count_mad
  )

  mt_high <- min(
    40,
    max(
      20,
      mt_med +
        5 * mt_mad
    )
  )

  keep <- (
    nFeature >=
      feature_low &
    nFeature <=
      feature_high &
    nCount >=
      count_low &
    nCount <=
      count_high &
    percent_mt <=
      mt_high
  )

  raw_n <- length(
    keep
  )

  retained_n <- sum(
    keep
  )

  retained_fraction <- (
    retained_n /
      raw_n
  )

  total_raw_cells <- (
    total_raw_cells +
      raw_n
  )

  total_retained_cells <- (
    total_retained_cells +
      retained_n
  )

  qc_threshold_list[[
    sample_name
  ]] <- data.table(
    technical_set = row$technical_set,
    patient_id = patient_id,
    timepoint = timepoint,
    sample_name = sample_name,
    hrs_accession = hrs,
    raw_cells = raw_n,
    retained_cells = retained_n,
    retained_fraction = retained_fraction,
    nFeature_low = feature_low,
    nFeature_high = feature_high,
    nCount_low = count_low,
    nCount_high = count_high,
    percent_mt_high = mt_high,
    nFeature_median = feature_med,
    nCount_median = count_med,
    percent_mt_median = mt_med
  )

  cat(
    sprintf(
      "[%02d/27] %-10s %-4s %-6s raw=%6d keep=%6d (%5.1f%%)  genes=%d-%d counts=%d-%d mt<=%.1f%%\n",
      i,
      sample_name,
      patient_id,
      timepoint,
      raw_n,
      retained_n,
      100 * retained_fraction,
      feature_low,
      feature_high,
      count_low,
      count_high,
      mt_high
    )
  )

  # Same frozen QC guards.
  if (
    retained_n < 500
  ) {

    stop(
      "\nQC_GATE_STOP: fewer than 500 retained cells in ",
      sample_name
    )
  }

  if (
    retained_fraction < 0.40
  ) {

    stop(
      "\nQC_GATE_STOP: >60% cells removed in ",
      sample_name
    )
  }

  counts <- counts[
    ,
    keep,
    drop = FALSE
  ]

  retained_ids <- cell_ids[
    keep
  ]

  matrix_list[[
    sample_name
  ]] <- counts

  meta_list[[
    sample_name
  ]] <- data.table(
    cell_id = retained_ids,
    barcode_raw = barcodes[
      keep
    ],
    technical_set = row$technical_set,
    patient_id = patient_id,
    timepoint = timepoint,
    sample_name = sample_name,
    hrs_accession = hrs,
    nCount_RNA = nCount[
      keep
    ],
    nFeature_RNA = nFeature[
      keep
    ],
    percent.mt = percent_mt[
      keep
    ]
  )

  rm(
    counts
  )

  gc(
    verbose = FALSE
  )
}

cat("\n")

qc_thresholds <- rbindlist(
  qc_threshold_list,
  use.names = TRUE
)

fwrite(
  qc_thresholds,
  QC_THRESHOLDS_OUT,
  bom = TRUE
)

cat(
  "Raw cells across 27 tumors      :",
  total_raw_cells,
  "\n"
)

cat(
  "QC-retained cells across 27     :",
  total_retained_cells,
  "\n"
)

cat(
  "Overall retained fraction       :",
  sprintf(
    "%.1f%%",
    100 *
      total_retained_cells /
      total_raw_cells
  ),
  "\n\n"
)

# ============================================================
# 7. Merge only QC-passed cells
# ============================================================

cat("===== 7. MERGE QC-PASSED ALL27 MATRICES =====\n\n")

counts_all <- Reduce(
  function(a, b) {

    cbind(
      a,
      b
    )
  },
  matrix_list
)

meta_all <- rbindlist(
  meta_list,
  use.names = TRUE
)

rm(
  matrix_list,
  meta_list
)

gc()

if (
  nrow(counts_all) !=
    length(
      analysis_ids
    )
) {

  stop(
    "MERGED_ANALYSIS_FEATURE_DIMENSION_MISMATCH"
  )
}

if (
  ncol(counts_all) !=
    nrow(meta_all)
) {

  stop(
    "MERGED_CELL_METADATA_DIMENSION_MISMATCH"
  )
}

if (
  !identical(
    colnames(counts_all),
    meta_all$cell_id
  )
) {

  stop(
    "MERGED_CELL_METADATA_ORDER_MISMATCH"
  )
}

cat(
  "Merged non-CORE matrix:",
  nrow(counts_all),
  "genes x",
  ncol(counts_all),
  "QC-passed cells\n\n"
)

# ============================================================
# 8. Create Seurat object
# ============================================================

cat("===== 8. CREATE ALL27 SEURAT OBJECT =====\n\n")

metadata_df <- as.data.frame(
  meta_all
)

rownames(
  metadata_df
) <- metadata_df$cell_id

obj <- CreateSeuratObject(
  counts = counts_all,
  assay = "RNA",
  project = "OMIX005710_ALL27",
  min.cells = 0,
  min.features = 0,
  meta.data = metadata_df
)

rm(
  counts_all,
  metadata_df
)

gc()

# Explicit overwrite with precomputed non-CORE QC values.
obj$nCount_RNA <- meta_all$nCount_RNA
obj$nFeature_RNA <- meta_all$nFeature_RNA
obj$percent.mt <- meta_all$percent.mt

cat(
  "Seurat cells:",
  ncol(obj),
  "\n"
)

cat(
  "Seurat features:",
  nrow(obj),
  "\n\n"
)

if (
  any(
    core_ids %in%
      rownames(obj)
  )
) {

  stop(
    "CORE_GENE_PRESENT_IN_CELL_ANALYSIS_OBJECT"
  )
}

cat(
  "✓ CORE gene rows absent from cell-analysis object.\n\n"
)

# ============================================================
# 9. Normalize + variable features
# ============================================================

cat("===== 9. NORMALIZATION + HVGs =====\n\n")

DefaultAssay(
  obj
) <- "RNA"

obj <- NormalizeData(
  obj,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

obj <- FindVariableFeatures(
  obj,
  selection.method = "vst",
  nfeatures = 2000,
  verbose = FALSE
)

hvg <- VariableFeatures(
  obj
)

if (
  any(
    core_ids %in%
      hvg
  )
) {

  stop(
    "CORE_GENE_ENTERED_HVG_SET"
  )
}

cat(
  "Variable features:",
  length(hvg),
  "\n"
)

cat(
  "CORE genes in HVGs:",
  0,
  "\n\n"
)

# ============================================================
# 10. PCA + clustering
# ============================================================

cat("===== 10. PCA + GRAPH CLUSTERING =====\n\n")

obj <- ScaleData(
  obj,
  features = hvg,
  block.size = 500,
  verbose = FALSE
)

obj <- RunPCA(
  obj,
  features = hvg,
  npcs = 30,
  seed.use = 20260814,
  verbose = FALSE
)

obj <- FindNeighbors(
  obj,
  reduction = "pca",
  dims = 1:30,
  verbose = FALSE
)

obj <- FindClusters(
  obj,
  resolution = 0.4,
  random.seed = 20260814,
  verbose = FALSE
)

cluster_id <- as.character(
  obj$seurat_clusters
)

cluster_levels <- sort(
  unique(
    as.integer(
      cluster_id
    )
  )
)

cluster_levels <- as.character(
  cluster_levels
)

cat(
  "Unsupervised clusters:",
  length(
    cluster_levels
  ),
  "\n\n"
)

# ============================================================
# 11. UMAP
# ============================================================

cat("===== 11. UMAP =====\n\n")

obj <- RunUMAP(
  obj,
  reduction = "pca",
  dims = 1:30,
  seed.use = 20260814,
  verbose = FALSE
)

cat(
  "✓ UMAP complete.\n\n"
)

# ============================================================
# 12. Extract 31-marker normalized expression
# ============================================================

cat("===== 12. FROZEN 31-MARKER EXPRESSION =====\n\n")

rna_data <- tryCatch(
  {
    LayerData(
      obj,
      assay = "RNA",
      layer = "data"
    )
  },
  error = function(e) {
    GetAssayData(
      obj,
      assay = "RNA",
      slot = "data"
    )
  }
)

marker_features <- marker_map$gene_id

if (
  !all(
    marker_features %in%
      rownames(
        rna_data
      )
  )
) {

  missing <- setdiff(
    marker_features,
    rownames(
      rna_data
    )
  )

  stop(
    "Mapped annotation feature missing from Seurat object: ",
    paste(
      missing,
      collapse = ", "
    )
  )
}

marker_expr <- rna_data[
  marker_features,
  ,
  drop = FALSE
]

rownames(
  marker_expr
) <- marker_map$marker

cat(
  "Frozen measurable marker matrix:",
  nrow(marker_expr),
  "markers x",
  ncol(marker_expr),
  "cells\n\n"
)

# ============================================================
# 13. Same constrained cluster-level annotation rule
#
# A cluster gets a broad label only when:
#
#   1. same module ranks #1 for:
#        - mean normalized expression
#        - mean marker detection fraction
#
#   2. supportive marker count >=
#        max(2, ceiling(25% of MEASURABLE module markers))
#
#   3. support marker detected in >=10% cluster cells
#
# Otherwise Unresolved.
#
# No treatment/timepoint information enters the decision.
# ============================================================

cat("===== 13. CONSTRAINED CLUSTER ANNOTATION =====\n\n")

marker_modules <- split(
  marker_map$marker,
  marker_map$marker_module
)

module_sizes_vector <- vapply(
  marker_modules,
  length,
  integer(1)
)

minimum_support <- pmax(
  2L,
  ceiling(
    module_sizes_vector *
      0.25
  )
)

minimum_support <- setNames(
  as.integer(
    minimum_support
  ),
  names(
    module_sizes_vector
  )
)

support_detection_fraction <- 0.10

cluster_evidence_list <- list()
cluster_decision_list <- list()

for (
  cl in cluster_levels
) {

  idx <- which(
    cluster_id == cl
  )

  module_rows <- list()

  for (
    module_name in names(
      marker_modules
    )
  ) {

    markers <- marker_modules[[
      module_name
    ]]

    m <- marker_expr[
      markers,
      idx,
      drop = FALSE
    ]

    marker_detection <- Matrix::rowMeans(
      m > 0
    )

    module_mean_expression <- (
      as.numeric(
        sum(m)
      ) /
        (
          nrow(m) *
            ncol(m)
        )
    )

    module_mean_detection <- mean(
      marker_detection
    )

    support_count <- sum(
      marker_detection >=
        support_detection_fraction
    )

    module_rows[[
      module_name
    ]] <- data.table(
      cluster = cl,
      module = module_name,
      n_cells = length(idx),
      mean_normalized_expression =
        module_mean_expression,
      mean_marker_detection_fraction =
        module_mean_detection,
      supportive_markers =
        as.integer(
          support_count
        ),
      required_supportive_markers =
        minimum_support[
          module_name
        ]
    )
  }

  module_table <- rbindlist(
    module_rows
  )

  expression_order <- module_table[
    order(
      -mean_normalized_expression,
      module
    )
  ]

  detection_order <- module_table[
    order(
      -mean_marker_detection_fraction,
      module
    )
  ]

  top_expression <- (
    expression_order$module[1]
  )

  top_detection <- (
    detection_order$module[1]
  )

  consensus <- (
    top_expression ==
      top_detection
  )

  top_row <- module_table[
    module ==
      top_expression
  ]

  support_pass <- (
    consensus &&
    top_row$supportive_markers >=
      top_row$required_supportive_markers &&
    top_row$mean_normalized_expression > 0
  )

  if (support_pass) {

    assigned <- top_expression

    status <- (
      "CONSENSUS_EXPRESSION_DETECTION"
    )

  } else {

    assigned <- "Unresolved"

    if (!consensus) {

      status <- (
        "EXPRESSION_DETECTION_DISAGREE"
      )

    } else {

      status <- (
        "INSUFFICIENT_MARKER_SUPPORT"
      )
    }
  }

  cluster_evidence_list[[
    cl
  ]] <- module_table

  cluster_decision_list[[
    cl
  ]] <- data.table(
    cluster = cl,
    n_cells = length(idx),
    n_samples = uniqueN(
      obj$sample_name[idx]
    ),
    n_patients = uniqueN(
      obj$patient_id[idx]
    ),
    before_cells = sum(
      obj$timepoint[idx] ==
        "Before"
    ),
    after_cells = sum(
      obj$timepoint[idx] ==
        "After"
    ),
    top_expression_module =
      top_expression,
    top_detection_module =
      top_detection,
    top_supportive_markers =
      top_row$supportive_markers,
    required_supportive_markers =
      top_row$required_supportive_markers,
    assigned_broad_type =
      assigned,
    decision_status =
      status
  )
}

cluster_evidence <- rbindlist(
  cluster_evidence_list
)

cluster_decisions <- rbindlist(
  cluster_decision_list
)

cluster_decisions[
  ,
  cluster_order :=
    as.integer(
      cluster
    )
]

setorder(
  cluster_decisions,
  cluster_order
)

cluster_decisions[
  ,
  cluster_order := NULL
]

fwrite(
  cluster_evidence,
  CLUSTER_EVIDENCE_OUT,
  bom = TRUE
)

fwrite(
  cluster_decisions,
  CLUSTER_ANNOTATION_OUT,
  bom = TRUE
)

print(
  cluster_decisions
)

# ============================================================
# 14. Apply labels to cells
# ============================================================

cat("\n===== 14. APPLY BROAD CELL LABELS =====\n\n")

lookup <- setNames(
  cluster_decisions$
    assigned_broad_type,
  cluster_decisions$
    cluster
)

obj$STEP22A_broad_celltype <- unname(
  lookup[
    as.character(
      obj$seurat_clusters
    )
  ]
)

if (
  anyNA(
    obj$STEP22A_broad_celltype
  )
) {

  stop(
    "BROAD_CELLTYPE_CELL_MAPPING_FAILED"
  )
}

broad_counts <- as.data.table(
  table(
    obj$STEP22A_broad_celltype
  )
)

setnames(
  broad_counts,
  c(
    "broad_celltype",
    "cells"
  )
)

broad_counts[
  ,
  fraction :=
    cells /
    sum(
      cells
    )
]

setorder(
  broad_counts,
  -cells
)

fwrite(
  broad_counts,
  BROAD_COUNTS_OUT,
  bom = TRUE
)

print(
  broad_counts
)

# ============================================================
# 15. Per-sample retained-cell summary
# ============================================================

cat("\n===== 15. ALL27 SAMPLE SUMMARY =====\n\n")

retained_summary <- as.data.table(
  obj@meta.data
)[
  ,
  .(
    retained_cells = .N,
    median_nFeature =
      median(
        nFeature_RNA
      ),
    median_nCount =
      median(
        nCount_RNA
      ),
    median_percent_mt =
      median(
        percent.mt
      )
  ),
  by = .(
    technical_set,
    patient_id,
    timepoint,
    sample_name,
    hrs_accession
  )
]

retained_summary[
  ,
  patient_order :=
    patient_number(
      patient_id
    )
]

retained_summary[
  ,
  timepoint_order :=
    factor(
      timepoint,
      levels = c(
        "Before",
        "After"
      )
    )
]

setorder(
  retained_summary,
  patient_order,
  timepoint_order
)

retained_summary[
  ,
  `:=`(
    patient_order = NULL,
    timepoint_order = NULL
  )
]

fwrite(
  retained_summary,
  RETAINED_OUT,
  bom = TRUE
)

# ============================================================
# 16. Fibroblast technical QC
# ============================================================

cat("===== 16. FIBROBLAST SAMPLE QC =====\n\n")

cell_dt <- data.table(
  cell_id = colnames(obj),
  barcode_raw =
    as.character(
      obj$barcode_raw
    ),
  technical_set =
    as.character(
      obj$technical_set
    ),
  patient_id =
    as.character(
      obj$patient_id
    ),
  timepoint =
    as.character(
      obj$timepoint
    ),
  sample_name =
    as.character(
      obj$sample_name
    ),
  hrs_accession =
    as.character(
      obj$hrs_accession
    ),
  cluster =
    as.character(
      obj$seurat_clusters
    ),
  broad_celltype =
    as.character(
      obj$STEP22A_broad_celltype
    ),
  nCount_RNA =
    as.numeric(
      obj$nCount_RNA
    ),
  nFeature_RNA =
    as.numeric(
      obj$nFeature_RNA
    ),
  percent_mt =
    as.numeric(
      obj$percent.mt
    )
)

fibro_qc <- cell_dt[
  ,
  .(
    n_total_QC_cells = .N,
    n_fibroblast_cells =
      sum(
        broad_celltype ==
          "Fibroblast_CAF"
      )
  ),
  by = .(
    technical_set,
    patient_id,
    timepoint,
    sample_name,
    hrs_accession
  )
]

fibro_qc[
  ,
  fibroblast_fraction :=
    n_fibroblast_cells /
      n_total_QC_cells
]

fibro_qc[
  ,
  frozen_fibroblast_QC :=
    fifelse(
      n_fibroblast_cells >= 20,
      "PASS_GE20",
      "FAIL_LT20"
    )
]

fibro_qc[
  ,
  patient_order :=
    patient_number(
      patient_id
    )
]

fibro_qc[
  ,
  timepoint_order :=
    factor(
      timepoint,
      levels = c(
        "Before",
        "After"
      )
    )
]

setorder(
  fibro_qc,
  patient_order,
  timepoint_order
)

fibro_qc[
  ,
  `:=`(
    patient_order = NULL,
    timepoint_order = NULL
  )
]

fwrite(
  fibro_qc,
  FIBRO_QC_OUT,
  bom = TRUE
)

print(
  fibro_qc
)

eligible_tumors <- sum(
  fibro_qc$
    frozen_fibroblast_QC ==
    "PASS_GE20"
)

cat(
  "\nTumor observations with >=20 fibroblasts:",
  eligible_tumors,
  "/ 27\n\n"
)

# ============================================================
# 17. Save cell annotation metadata
# ============================================================

cat("===== 17. SAVE CELL METADATA / FIBROBLAST IDs =====\n\n")

fwrite(
  cell_dt,
  CELL_METADATA_OUT,
  bom = TRUE
)

fibro_ids <- cell_dt[
  broad_celltype ==
    "Fibroblast_CAF",
  .(
    cell_id,
    barcode_raw,
    technical_set,
    patient_id,
    timepoint,
    sample_name,
    hrs_accession,
    cluster
  )
]

fwrite(
  fibro_ids,
  FIBRO_IDS_OUT,
  bom = TRUE
)

cat(
  "Annotated QC-passed cells:",
  nrow(
    cell_dt
  ),
  "\n"
)

cat(
  "Fibroblast_CAF cell IDs:",
  nrow(
    fibro_ids
  ),
  "\n\n"
)

# ============================================================
# 18. UMAP coordinates + figures
# ============================================================

cat("===== 18. SAVE FIGURES =====\n\n")

p_cluster <- DimPlot(
  obj,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE
) +
  labs(
    title = "All 27 Tumors: Unsupervised Clusters"
  )

ggsave(
  filename = file.path(
    FIGURE_DIR,
    "STEP22A_06C_ALL27_UMAP_CLUSTERS.png"
  ),
  plot = p_cluster,
  width = 10,
  height = 8,
  dpi = 180
)

p_type <- DimPlot(
  obj,
  reduction = "umap",
  group.by = "STEP22A_broad_celltype",
  label = TRUE,
  repel = TRUE
) +
  labs(
    title = "All 27 Tumors: Constrained Broad Cell Types"
  )

ggsave(
  filename = file.path(
    FIGURE_DIR,
    "STEP22A_06C_ALL27_UMAP_BROAD_CELLTYPE.png"
  ),
  plot = p_type,
  width = 10,
  height = 8,
  dpi = 180
)

umap <- as.data.table(
  Embeddings(
    obj,
    "umap"
  ),
  keep.rownames = "cell_id"
)

umap[
  ,
  broad_celltype :=
    cell_dt[
      match(
        cell_id,
        cell_dt$cell_id
      ),
      broad_celltype
    ]
]

umap[
  ,
  sample_name :=
    cell_dt[
      match(
        cell_id,
        cell_dt$cell_id
      ),
      sample_name
    ]
]

fwrite(
  umap,
  file.path(
    META_DIR,
    "STEP22A_06C_ALL27_UMAP_COORDINATES.csv"
  ),
  bom = TRUE
)

cat(
  "✓ UMAP figures/coordinates saved.\n\n"
)

# ============================================================
# 19. Explicit anti-circularity audit
# ============================================================

cat("===== 19. ANTI-CIRCULARITY AUDIT =====\n\n")

anti <- data.table(
  item = c(
    "true_all27_common_ENSG",
    "CORE_genes_present_in_raw_common_universe",
    "CORE_genes_removed_before_QC",
    "CORE_genes_removed_before_normalization",
    "CORE_genes_removed_before_HVG",
    "CORE_genes_removed_before_PCA",
    "CORE_genes_removed_before_clustering",
    "CORE_genes_removed_before_annotation",
    "CORE2_scoring_performed",
    "CORE4_scoring_performed",
    "treatment_used_for_QC_thresholds",
    "treatment_used_for_clustering",
    "treatment_used_for_annotation",
    "FindAllMarkers_used",
    "KRT5_KRT14_treated_as_zero",
    "histology_used_for_cell_annotation"
  ),
  status = c(
    as.character(
      nrow(common)
    ),
    "YES_4_OF_4",
    "YES",
    "YES",
    "YES",
    "YES",
    "YES",
    "YES",
    "NO",
    "NO",
    "NO",
    "NO",
    "NO",
    "NO",
    "NO_STRUCTURAL_NA",
    "NO"
  )
)

fwrite(
  anti,
  ANTI_CIRC_OUT,
  bom = TRUE
)

print(
  anti
)

# ============================================================
# 20. Small annotation bundle only
#
# Deliberately DO NOT save a second huge Seurat object.
# 6D will rebuild raw fibroblast pseudobulk directly from
# source matrices using the frozen fibroblast cell-ID list.
# ============================================================

cat("\n===== 20. SAVE SMALL ANNOTATION BUNDLE =====\n\n")

bundle <- list(
  all27_common_features =
    common[
      ,
      .(
        all27_common_index,
        gene_id,
        gene_symbol
      )
    ],
  core_withheld =
    core_map[
      ,
      .(
        gene_id,
        gene_symbol
      )
    ],
  analysis_feature_count =
    length(
      analysis_ids
    ),
  marker_mapping =
    marker_map,
  sample_QC =
    qc_thresholds,
  cluster_module_evidence =
    cluster_evidence,
  cluster_annotation =
    cluster_decisions,
  broad_celltype_counts =
    broad_counts,
  fibroblast_sample_QC =
    fibro_qc,
  rules = list(
    frozen_biological_markers = 33L,
    joint_measurable_markers = 31L,
    KRT5_KRT14 =
      "STRUCTURAL_NA_NOT_ZERO",
    CORE_genes_used_for_cell_selection =
      FALSE,
    treatment_used_for_annotation =
      FALSE,
    FindAllMarkers_used =
      FALSE,
    minimum_fibroblast_cells =
      20L,
    histology_mapping =
      "INCOMPLETE"
  )
)

saveRDS(
  bundle,
  BUNDLE_OUT,
  compress = "gzip"
)

cat(
  "✓ Small bundle saved.\n"
)

cat(
  "✓ Full Seurat object intentionally NOT saved.\n\n"
)

# ============================================================
# 21. Frozen Step19-21 integrity check
# ============================================================

cat("===== 21. FROZEN PROJECT CHECK =====\n\n")

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
      nrow(frozen)
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
      is.finite(expected) &&
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
# 22. Session info
# ============================================================

sink(
  file.path(
    LOG_DIR,
    "STEP22A_06C_sessionInfo.txt"
  )
)

print(
  sessionInfo()
)

sink()

# ============================================================
# 23. Summary
# ============================================================

fibro_total <- sum(
  broad_counts[
    broad_celltype ==
      "Fibroblast_CAF",
    cells
  ]
)

unresolved_total <- sum(
  broad_counts[
    broad_celltype ==
      "Unresolved",
    cells
  ]
)

summary_lines <- c(
  "STEP22A-6C ALL27 QC + CONSTRAINED BROAD ANNOTATION",
  paste0(
    "Tumor observations: ",
    nrow(sample_table)
  ),
  paste0(
    "Patients represented: ",
    uniqueN(sample_table$patient_id)
  ),
  paste0(
    "True all-27 common ENSG: ",
    nrow(common)
  ),
  paste0(
    "CORE genes withheld from cell analysis: ",
    length(core_ids)
  ),
  paste0(
    "Non-CORE cell-analysis features: ",
    length(analysis_ids)
  ),
  paste0(
    "Raw cells: ",
    total_raw_cells
  ),
  paste0(
    "QC-retained cells: ",
    total_retained_cells
  ),
  paste0(
    "Clusters: ",
    nrow(cluster_decisions)
  ),
  paste0(
    "Fibroblast_CAF cells: ",
    fibro_total
  ),
  paste0(
    "Unresolved cells: ",
    unresolved_total
  ),
  paste0(
    "Tumor observations >=20 fibroblasts: ",
    eligible_tumors,
    "/27"
  ),
  "Frozen biological marker dictionary: 33 retained",
  "Joint measurable annotation markers: 31",
  "KRT5/KRT14 interpreted as structural NA: YES",
  "CORE genes used for QC/HVG/PCA/clustering/annotation: NO",
  "CORE2 scoring performed: NO",
  "CORE4 scoring performed: NO",
  "Treatment inference performed: NO",
  "Histology mapping: INCOMPLETE",
  "Primary CORE2 inference: BLOCKED",
  "Full Seurat object saved: NO",
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
cat("STEP22A-6C COMPLETE\n")
cat("====================================================================\n\n")

cat(
  "TUMOR OBSERVATIONS            : 27\n"
)

cat(
  "PATIENTS REPRESENTED          : 21\n"
)

cat(
  "TRUE ALL27 COMMON ENSG        :",
  nrow(common),
  "\n"
)

cat(
  "CORE GENES WITHHELD           :",
  length(core_ids),
  "/ 4\n"
)

cat(
  "NON-CORE ANALYSIS FEATURES    :",
  length(analysis_ids),
  "\n\n"
)

cat(
  "RAW CELLS                     :",
  total_raw_cells,
  "\n"
)

cat(
  "QC-RETAINED CELLS             :",
  total_retained_cells,
  "\n"
)

cat(
  "RETAINED FRACTION             :",
  sprintf(
    "%.1f%%",
    100 *
      total_retained_cells /
      total_raw_cells
  ),
  "\n\n"
)

cat(
  "UNSUPERVISED CLUSTERS         :",
  nrow(cluster_decisions),
  "\n\n"
)

cat(
  "BROAD CELL-TYPE COUNTS:\n"
)

print(
  broad_counts
)

cat("\n")

cat(
  "FIBROBLAST_CAF CELLS          :",
  fibro_total,
  "\n"
)

cat(
  "TUMORS >=20 FIBROBLASTS      :",
  eligible_tumors,
  "/ 27\n\n"
)

cat(
  "FROZEN BIOLOGICAL MARKERS     : 33 / 33 RETAINED\n"
)

cat(
  "JOINT MEASURABLE MARKERS      : 31 / 33\n"
)

cat(
  "KRT5/KRT14                    : STRUCTURAL NA, NOT ZERO\n\n"
)

cat(
  "CORE USED FOR QC              : NO\n"
)

cat(
  "CORE USED FOR HVG/PCA         : NO\n"
)

cat(
  "CORE USED FOR CLUSTERING      : NO\n"
)

cat(
  "CORE USED FOR ANNOTATION      : NO\n"
)

cat(
  "FindAllMarkers USED           : NO\n"
)

cat(
  "TREATMENT USED FOR ANNOTATION : NO\n\n"
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
  "PRIMARY CORE2 INFERENCE       : BLOCKED\n\n"
)

cat(
  "FULL SEURAT OBJECT SAVED      : NO\n"
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
  "STEP22A-6D = BUILD FULL ALL27 FIBROBLAST\n"
)

cat(
  "RAW-COUNT PSEUDOBULKS FROM THE FROZEN\n"
)

cat(
  "STEP22A-6C FIBROBLAST CELL-ID LIST.\n"
)

cat(
  "STILL NO CORE2 SCORING.\n\n"
)

cat(
  "STOP HERE.\n"
)
