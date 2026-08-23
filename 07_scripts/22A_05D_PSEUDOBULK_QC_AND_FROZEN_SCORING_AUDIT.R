# ============================================================
# ESCC Neoadjuvant Project
# STEP22A-5D
#
# PSEUDOBULK TECHNICAL QC
# + FROZEN STEP19L SCORING IMPLEMENTATION AUDIT
#
# IMPORTANT:
#
# ALLOWED:
#   - inspect pseudobulk library sizes
#   - inspect fibroblast cell numbers
#   - inspect detected-gene counts
#   - flag technical depth outliers descriptively
#   - verify CORE genes exist in aligned feature universe
#   - read / hash / audit frozen Step19L script
#   - extract frozen scoring/normalization code snippets
#
# NOT ALLOWED:
#   - calculate CORE2
#   - calculate CORE4
#   - calculate CDKN1B/GSN treatment deltas
#   - inspect gene-level Before/After effects
#   - perform normalization on OMIX pseudobulk
#   - run paired tests
#   - exclude samples based on new thresholds
#   - modify Step19-21
#   - modify manuscript
#
# Histology remains unresolved.
# ============================================================

options(
  stringsAsFactors = FALSE,
  timeout = 1200
)

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
})

cat("\n")
cat("====================================================================\n")
cat("STEP22A-5D: PSEUDOBULK QC + FROZEN SCORING AUDIT\n")
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

PB_OBJECT <- file.path(
  OBJECT_DIR,
  paste0(
    "OMIX005710_FIBROBLAST_PATIENT_TIMEPOINT_",
    "PSEUDOBULK_HISTOLOGY_UNRESOLVED.rds"
  )
)

STEP19L_SCRIPT <- file.path(
  PROJECT,
  "07_scripts",
  "19L_CONSERVED_FIBROBLAST_CORE_VALIDATION.R"
)

GATE_FILE <- file.path(
  META_DIR,
  "STEP22A_ANALYSIS_GATE_STATUS.csv"
)

FROZEN_INV <- file.path(
  META_DIR,
  "STEP22A_PREANALYSIS_FROZEN_FILE_INVENTORY.csv"
)

# ============================================================
# Helpers
# ============================================================

sha256_file <- function(path) {

  cmd <- sprintf(
    "shasum -a 256 %s",
    shQuote(path)
  )

  out <- system(
    cmd,
    intern = TRUE
  )

  sub(
    "\\s+.*$",
    "",
    out[1]
  )
}

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
# 1. Analysis gate
# ============================================================

cat("===== 1. ANALYSIS GATE =====\n\n")

if (!file.exists(GATE_FILE)) {
  stop(
    "Missing analysis gate file."
  )
}

gate <- fread(
  GATE_FILE,
  encoding = "UTF-8"
)

get_allowed <- function(x) {

  z <- gate[
    gate == x
  ]

  if (nrow(z) != 1) {
    stop(
      "Cannot resolve gate: ",
      x
    )
  }

  as.character(
    z$allowed[1]
  )
}

technical_allowed <- get_allowed(
  "technical_preprocessing"
)

core2_allowed <- get_allowed(
  "paired_CORE2_test"
)

manuscript_allowed <- get_allowed(
  "manuscript_update"
)

cat(
  "Technical preprocessing :",
  technical_allowed,
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

if (core2_allowed == "YES") {
  stop(
    "CORE2 gate unexpectedly open."
  )
}

if (manuscript_allowed == "YES") {
  stop(
    "Manuscript gate unexpectedly open."
  )
}

cat(
  "✓ Technical audit allowed.\n"
)

cat(
  "✓ CORE2 inference remains blocked.\n\n"
)

# ============================================================
# 2. Load frozen pseudobulk technical object
# ============================================================

cat("===== 2. LOAD PSEUDOBULK OBJECT =====\n\n")

if (!file.exists(PB_OBJECT)) {
  stop(
    "Missing pseudobulk object:\n",
    PB_OBJECT
  )
}

pb <- readRDS(
  PB_OBJECT
)

required_elements <- c(
  "counts",
  "patient_timepoint_metadata",
  "feature_metadata"
)

missing_elements <- setdiff(
  required_elements,
  names(pb)
)

if (
  length(
    missing_elements
  ) > 0
) {

  stop(
    "Pseudobulk object missing: ",
    paste(
      missing_elements,
      collapse = ", "
    )
  )
}

counts <- pb$counts

meta <- as.data.table(
  pb$patient_timepoint_metadata
)

features <- as.data.table(
  pb$feature_metadata
)

cat(
  "Features             :",
  nrow(counts),
  "\n"
)

cat(
  "Patient-timepoints   :",
  ncol(counts),
  "\n"
)

cat(
  "Metadata rows        :",
  nrow(meta),
  "\n\n"
)

if (
  nrow(counts) != 18354
) {

  cat(
    "WARNING: expected 18,354 common ENSG features; observed ",
    nrow(counts),
    ".\n",
    sep = ""
  )
}

if (
  ncol(counts) != 12 ||
  nrow(meta) != 12
) {

  stop(
    "EXPECTED_12_PATIENT_TIMEPOINTS"
  )
}

if (
  !identical(
    colnames(counts),
    meta$patient_timepoint
  )
) {

  stop(
    "PSEUDOBULK_METADATA_ORDER_MISMATCH"
  )
}

cat(
  "✓ 12 patient-timepoint pseudobulks aligned.\n\n"
)

# ============================================================
# 3. Basic technical QC
# ============================================================

cat("===== 3. PSEUDOBULK TECHNICAL QC =====\n\n")

meta[
  ,
  observed_library_size :=
    as.numeric(
      Matrix::colSums(
        counts
      )
    )
]

meta[
  ,
  observed_detected_genes :=
    as.integer(
      Matrix::colSums(
        counts > 0
      )
    )
]

meta[
  ,
  counts_per_fibroblast :=
    observed_library_size /
    n_fibroblast_cells
]

meta[
  ,
  log10_library_size :=
    log10(
      observed_library_size + 1
    )
]

meta[
  ,
  log10_counts_per_fibroblast :=
    log10(
      counts_per_fibroblast + 1
    )
]

# ============================================================
# 4. Verify metadata library sizes
# ============================================================

cat("===== 4. LIBRARY SIZE CONSISTENCY =====\n\n")

meta[
  ,
  library_size_match :=
    observed_library_size ==
    as.numeric(
      library_size
    )
]

cat(
  "Library-size metadata matches raw matrix:",
  sum(
    meta$library_size_match
  ),
  "/",
  nrow(meta),
  "\n"
)

if (
  !all(
    meta$library_size_match
  )
) {

  print(
    meta[
      library_size_match == FALSE
    ]
  )

  stop(
    "LIBRARY_SIZE_METADATA_MISMATCH"
  )
}

cat(
  "✓ All library sizes verified directly from raw pseudobulk counts.\n\n"
)

# ============================================================
# 5. Frozen 20-cell QC only
#
# No new exclusion threshold is introduced.
# ============================================================

cat("===== 5. FROZEN FIBROBLAST CELL THRESHOLD =====\n\n")

meta[
  ,
  frozen_cell_threshold_status :=
    fifelse(
      n_fibroblast_cells >= 20,
      "PASS_GE20",
      "FAIL_LT20"
    )
]

print(
  meta[
    ,
    .(
      patient_id,
      timepoint,
      n_fibroblast_cells,
      observed_library_size,
      observed_detected_genes,
      counts_per_fibroblast,
      frozen_cell_threshold_status
    )
  ]
)

if (
  any(
    meta$frozen_cell_threshold_status !=
      "PASS_GE20"
  )
) {

  stop(
    "FROZEN_FIBROBLAST_CELL_THRESHOLD_FAILED"
  )
}

cat(
  "\n✓ All 12 patient-timepoints pass the frozen >=20 fibroblast rule.\n\n"
)

# ============================================================
# 6. Descriptive depth-outlier audit
#
# IMPORTANT:
# These are REVIEW FLAGS only.
#
# They are NOT exclusion criteria.
# They are NOT allowed to remove a sample.
#
# Robust comparison uses all 12 observations together and
# ignores treatment labels when defining thresholds.
# ============================================================

cat("===== 6. DESCRIPTIVE DEPTH OUTLIER AUDIT =====\n\n")

lib_med <- median(
  meta$log10_library_size
)

lib_mad <- robust_scale(
  meta$log10_library_size
)

cpf_med <- median(
  meta$log10_counts_per_fibroblast
)

cpf_mad <- robust_scale(
  meta$log10_counts_per_fibroblast
)

det_med <- median(
  meta$observed_detected_genes
)

det_mad <- robust_scale(
  meta$observed_detected_genes
)

meta[
  ,
  library_depth_review :=
    abs(
      log10_library_size -
        lib_med
    ) >
    3 * lib_mad
]

meta[
  ,
  per_cell_depth_review :=
    abs(
      log10_counts_per_fibroblast -
        cpf_med
    ) >
    3 * cpf_mad
]

meta[
  ,
  detected_gene_review :=
    abs(
      observed_detected_genes -
        det_med
    ) >
    3 * det_mad
]

meta[
  ,
  any_technical_review :=
    library_depth_review |
    per_cell_depth_review |
    detected_gene_review
]

meta[
  ,
  technical_review_status :=
    fifelse(
      any_technical_review,
      "REVIEW_ONLY_NO_EXCLUSION",
      "NO_ROBUST_OUTLIER_FLAG"
    )
]

meta[
  ,
  patient_order :=
    patient_number(
      patient_id
    )
]

meta[
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
  meta,
  patient_order,
  timepoint_order
)

meta[
  ,
  `:=`(
    patient_order = NULL,
    timepoint_order = NULL
  )
]

print(
  meta[
    ,
    .(
      patient_id,
      timepoint,
      n_fibroblast_cells,
      observed_library_size,
      counts_per_fibroblast,
      observed_detected_genes,
      technical_review_status
    )
  ]
)

cat(
  "\nTechnical review flags:",
  sum(
    meta$any_technical_review
  ),
  "/ 12\n"
)

cat(
  "Samples excluded because of these flags: 0\n\n"
)

# ============================================================
# 7. Library-depth range
# ============================================================

cat("===== 7. DEPTH RANGE SUMMARY =====\n\n")

min_lib_idx <- which.min(
  meta$observed_library_size
)

max_lib_idx <- which.max(
  meta$observed_library_size
)

min_lib <- meta[
  min_lib_idx
]

max_lib <- meta[
  max_lib_idx
]

library_ratio <- (
  max_lib$observed_library_size /
  min_lib$observed_library_size
)

cat(
  "Lowest library :",
  min_lib$patient_id,
  min_lib$timepoint,
  "=",
  min_lib$observed_library_size,
  "\n"
)

cat(
  "Highest library:",
  max_lib$patient_id,
  max_lib$timepoint,
  "=",
  max_lib$observed_library_size,
  "\n"
)

cat(
  "Max/min ratio  :",
  sprintf(
    "%.2f x",
    library_ratio
  ),
  "\n\n"
)

# ============================================================
# 8. CORE gene PRESENCE audit only
#
# NO COUNTS ARE READ OR PRINTED FOR THESE GENES.
# NO SCORE IS CALCULATED.
# ============================================================

cat("===== 8. FROZEN CORE GENE PRESENCE AUDIT =====\n\n")

core_gene_names <- c(
  "CDKN1B",
  "GSN",
  "BGN",
  "TIMP1"
)

core_presence <- rbindlist(
  lapply(
    core_gene_names,
    function(g) {

      hits <- features[
        toupper(
          gene_symbol
        ) ==
          toupper(
            g
          )
      ]

      data.table(
        gene = g,
        present = (
          nrow(hits) > 0
        ),
        ENSG_count = nrow(hits),
        ENSG_ids = paste(
          hits$gene_id,
          collapse = ";"
        )
      )
    }
  )
)

print(
  core_presence
)

if (
  !all(
    core_presence$present
  )
) {

  stop(
    "FROZEN_CORE_GENE_MISSING_FROM_EXTERNAL_FEATURE_UNIVERSE"
  )
}

cat(
  "\n✓ CDKN1B / GSN / BGN / TIMP1 are technically available.\n"
)

cat(
  "✓ Their pseudobulk expression values were NOT inspected.\n\n"
)

# ============================================================
# 9. Verify frozen Step19L script exists
# ============================================================

cat("===== 9. LOCATE FROZEN STEP19L IMPLEMENTATION =====\n\n")

if (!file.exists(STEP19L_SCRIPT)) {

  stop(
    "Missing frozen Step19L script:\n",
    STEP19L_SCRIPT
  )
}

step19_lines <- readLines(
  STEP19L_SCRIPT,
  warn = FALSE,
  encoding = "UTF-8"
)

cat(
  "Step19L script:",
  STEP19L_SCRIPT,
  "\n"
)

cat(
  "Lines:",
  length(
    step19_lines
  ),
  "\n"
)

step19_sha <- sha256_file(
  STEP19L_SCRIPT
)

cat(
  "SHA256:",
  step19_sha,
  "\n\n"
)

# ============================================================
# 10. Confirm key frozen definitions exist
# ============================================================

cat("===== 10. FROZEN DEFINITION KEYWORD CHECK =====\n\n")

required_terms <- c(
  "CORE2",
  "CORE4",
  "CDKN1B",
  "GSN",
  "BGN",
  "TIMP1"
)

term_check <- rbindlist(
  lapply(
    required_terms,
    function(term) {

      hits <- grep(
        term,
        step19_lines,
        ignore.case = TRUE,
        fixed = TRUE
      )

      data.table(
        term = term,
        line_hits = length(
          hits
        ),
        present = length(
          hits
        ) > 0
      )
    }
  )
)

print(
  term_check
)

if (
  !all(
    term_check$present
  )
) {

  stop(
    "FROZEN_STEP19L_CORE_DEFINITION_INCOMPLETE"
  )
}

cat(
  "\n✓ Frozen CORE2/CORE4/gene definitions are visible in Step19L.\n\n"
)

# ============================================================
# 11. Scoring / normalization implementation text audit
#
# We extract SOURCE CODE ONLY.
#
# We do not execute any extracted scoring expression.
# ============================================================

cat("===== 11. SCORING IMPLEMENTATION SOURCE AUDIT =====\n\n")

audit_patterns <- c(
  "CORE2",
  "CORE4",
  "module",
  "score",
  "scale\\(",
  "cpm\\(",
  "logCPM",
  "DGEList",
  "calcNormFactors",
  "voom",
  "normalize",
  "log2",
  "rowMeans",
  "colMeans",
  "sweep\\(",
  "standard"
)

matched_lines <- integer(0)

for (
  pattern in audit_patterns
) {

  matched_lines <- union(
    matched_lines,
    grep(
      pattern,
      step19_lines,
      ignore.case = TRUE,
      perl = TRUE
    )
  )
}

matched_lines <- sort(
  matched_lines
)

# Add +/-4 lines of context around each hit.
context_lines <- sort(
  unique(
    unlist(
      lapply(
        matched_lines,
        function(i) {

          seq(
            max(
              1,
              i - 4
            ),
            min(
              length(
                step19_lines
              ),
              i + 4
            )
          )
        }
      )
    )
  )
)

source_audit <- data.table(
  line_number = context_lines,
  code = step19_lines[
    context_lines
  ]
)

SCORING_SOURCE_OUT <- file.path(
  META_DIR,
  "STEP22A_05D_FROZEN_STEP19L_SCORING_SOURCE_AUDIT.csv"
)

fwrite(
  source_audit,
  SCORING_SOURCE_OUT,
  bom = TRUE
)

cat(
  "Matched scoring/normalization source lines:",
  length(
    matched_lines
  ),
  "\n"
)

cat(
  "Context lines saved:",
  nrow(
    source_audit
  ),
  "\n\n"
)

# ============================================================
# 12. Detect implementation vocabulary
#
# This is descriptive source-code detection ONLY.
# It does not choose or run a method.
# ============================================================

implementation_terms <- list(

  edgeR_cpm =
    "cpm\\(",

  DGEList =
    "DGEList",

  calcNormFactors =
    "calcNormFactors",

  voom =
    "\\bvoom\\b",

  scale =
    "scale\\(",

  log2 =
    "log2\\(",

  rowMeans =
    "rowMeans",

  colMeans =
    "colMeans",

  weighted_sum =
    "weight|weights",

  zscore_text =
    "z.?score|standardiz"
)

implementation_audit <- rbindlist(
  lapply(
    names(
      implementation_terms
    ),
    function(name) {

      pattern <- implementation_terms[[
        name
      ]]

      hits <- grep(
        pattern,
        step19_lines,
        ignore.case = TRUE,
        perl = TRUE
      )

      data.table(
        implementation_term = name,
        detected = (
          length(
            hits
          ) > 0
        ),
        hit_count = length(
          hits
        ),
        line_numbers = paste(
          hits,
          collapse = ";"
        )
      )
    }
  )
)

IMPLEMENTATION_OUT <- file.path(
  META_DIR,
  "STEP22A_05D_FROZEN_SCORING_IMPLEMENTATION_TERMS.csv"
)

fwrite(
  implementation_audit,
  IMPLEMENTATION_OUT,
  bom = TRUE
)

print(
  implementation_audit
)

cat("\n")

# ============================================================
# 13. Save exact relevant text excerpt
# ============================================================

SCORING_TEXT_OUT <- file.path(
  META_DIR,
  "STEP22A_05D_FROZEN_STEP19L_SCORING_SOURCE.txt"
)

con <- file(
  SCORING_TEXT_OUT,
  open = "wt",
  encoding = "UTF-8"
)

writeLines(
  c(
    "STEP22A-5D FROZEN STEP19L SOURCE AUDIT",
    paste0(
      "Source: ",
      STEP19L_SCRIPT
    ),
    paste0(
      "SHA256: ",
      step19_sha
    ),
    "",
    "IMPORTANT:",
    "This is a source-code audit only.",
    "No scoring expression was executed.",
    "",
    "========================================",
    ""
  ),
  con
)

for (
  i in context_lines
) {

  writeLines(
    sprintf(
      "%05d | %s",
      i,
      step19_lines[i]
    ),
    con
  )
}

close(
  con
)

cat(
  "✓ Frozen Step19L implementation excerpt saved.\n"
)

cat(
  "✓ No extracted scoring code was executed.\n\n"
)

# ============================================================
# 14. Audit frozen Step19L result files
#
# File names / sizes / hashes only.
# Do NOT read score values.
# ============================================================

cat("===== 12. FROZEN STEP19L OUTPUT PROVENANCE =====\n\n")

candidate_roots <- c(
  file.path(
    PROJECT,
    "04_results",
    "cross_dataset",
    "pathway_analysis",
    "Step19L"
  ),
  file.path(
    PROJECT,
    "06_tables",
    "cross_dataset",
    "pathway_analysis",
    "Step19L"
  ),
  file.path(
    PROJECT,
    "05_figures",
    "cross_dataset",
    "pathway_analysis",
    "Step19L"
  )
)

step19_files <- character(0)

for (
  root in candidate_roots
) {

  if (
    dir.exists(
      root
    )
  ) {

    step19_files <- c(
      step19_files,
      list.files(
        root,
        recursive = TRUE,
        full.names = TRUE
      )
    )
  }
}

step19_files <- unique(
  step19_files[
    file.exists(
      step19_files
    ) &
    !file.info(
      step19_files
    )$isdir
  ]
)

if (
  length(
    step19_files
  ) == 0
) {

  cat(
    "WARNING: no Step19L result files found in standard directories.\n"
  )

  frozen_output_manifest <- data.table(
    file = character(0),
    size_bytes = numeric(0),
    sha256 = character(0)
  )

} else {

  frozen_output_manifest <- rbindlist(
    lapply(
      step19_files,
      function(path) {

        data.table(
          file = file.path(
            ".",
            sub(
              paste0(
                "^",
                gsub(
                  "([.()+*?^$|{}\\[\\]\\\\])",
                  "\\\\\\1",
                  PROJECT
                ),
                "/?"
              ),
              "",
              path
            )
          ),
          size_bytes = file.info(
            path
          )$size,
          sha256 = sha256_file(
            path
          )
        )
      }
    )
  )
}

OUTPUT_MANIFEST <- file.path(
  META_DIR,
  "STEP22A_05D_FROZEN_STEP19L_OUTPUT_MANIFEST.csv"
)

fwrite(
  frozen_output_manifest,
  OUTPUT_MANIFEST,
  bom = TRUE
)

cat(
  "Frozen Step19L result files inventoried:",
  nrow(
    frozen_output_manifest
  ),
  "\n"
)

cat(
  "Result values read:",
  "NO\n\n"
)

# ============================================================
# 15. Save technical pseudobulk QC
# ============================================================

QC_OUT <- file.path(
  META_DIR,
  "STEP22A_05D_PSEUDOBULK_TECHNICAL_QC.csv"
)

fwrite(
  meta,
  QC_OUT,
  bom = TRUE
)

# ============================================================
# 16. Paired technical symmetry report
#
# Descriptive counts/depth only.
# NO gene expression or score.
# ============================================================

paired_depth <- dcast(
  meta[
    ,
    .(
      patient_id,
      timepoint,
      n_fibroblast_cells,
      observed_library_size,
      counts_per_fibroblast,
      observed_detected_genes
    )
  ],
  patient_id ~ timepoint,
  value.var = c(
    "n_fibroblast_cells",
    "observed_library_size",
    "counts_per_fibroblast",
    "observed_detected_genes"
  )
)

PAIRED_DEPTH_OUT <- file.path(
  META_DIR,
  "STEP22A_05D_PAIRED_TECHNICAL_DEPTH_AUDIT.csv"
)

fwrite(
  paired_depth,
  PAIRED_DEPTH_OUT,
  bom = TRUE
)

# ============================================================
# 17. Explicit inference firewall
# ============================================================

firewall <- data.table(

  action = c(
    "raw_pseudobulk_loaded",
    "library_size_inspected",
    "fibroblast_cell_counts_inspected",
    "detected_gene_counts_inspected",
    "CORE_gene_presence_checked",
    "CORE_gene_expression_inspected",
    "normalization_applied",
    "CORE2_score_calculated",
    "CORE4_score_calculated",
    "Before_After_gene_effect_inspected",
    "paired_statistical_test_run",
    "sample_excluded_by_new_QC_rule",
    "histology_mapping_complete",
    "primary_inference_allowed"
  ),

  status = c(
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
    "NO",
    "NO",
    "NO"
  )
)

FIREWALL_OUT <- file.path(
  META_DIR,
  "STEP22A_05D_INFERENCE_FIREWALL.csv"
)

fwrite(
  firewall,
  FIREWALL_OUT,
  bom = TRUE
)

# ============================================================
# 18. Frozen Step19-21 integrity check
# ============================================================

cat("===== 13. FROZEN PROJECT INTEGRITY =====\n\n")

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

    if (
      !file.exists(
        path
      )
    ) {

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
# 19. Summary
# ============================================================

SUMMARY_OUT <- file.path(
  META_DIR,
  "STEP22A_05D_PSEUDOBULK_QC_SCORING_AUDIT_SUMMARY.txt"
)

summary_lines <- c(

  "STEP22A-5D PSEUDOBULK TECHNICAL QC + FROZEN SCORING AUDIT",

  paste0(
    "Patient-timepoint pseudobulks: ",
    nrow(
      meta
    )
  ),

  paste0(
    "Common ENSG features: ",
    nrow(
      counts
    )
  ),

  paste0(
    "Technical review flags: ",
    sum(
      meta$any_technical_review
    )
  ),

  "Samples excluded by new technical rules: 0",

  paste0(
    "Library-size max/min ratio: ",
    sprintf(
      "%.3f",
      library_ratio
    )
  ),

  "All 12 pass frozen >=20 fibroblast rule: YES",

  "CDKN1B present: YES",
  "GSN present: YES",
  "BGN present: YES",
  "TIMP1 present: YES",

  paste0(
    "Frozen Step19L SHA256: ",
    step19_sha
  ),

  paste0(
    "Frozen Step19L scoring-source context lines saved: ",
    nrow(
      source_audit
    )
  ),

  "Frozen Step19L output VALUES read: NO",

  "External pseudobulk normalization performed: NO",

  "CORE2 scoring performed: NO",

  "CORE4 scoring performed: NO",

  "Before/After inference performed: NO",

  "Histology mapping: INCOMPLETE",

  "EASC patient: UNRESOLVED",

  "Primary CORE2 inference: BLOCKED",

  "Frozen Step19-21 modified: NO",

  "Manuscript modified: NO"
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
cat("STEP22A-5D COMPLETE\n")
cat("====================================================================\n\n")

cat(
  "PATIENT-TIMEPOINT PSEUDOBULKS :",
  nrow(
    meta
  ),
  "\n"
)

cat(
  "COMMON ENSG FEATURES          :",
  nrow(
    counts
  ),
  "\n"
)

cat(
  "ALL >=20 FIBROBLASTS          : YES\n\n"
)

cat(
  "LOWEST LIBRARY                :",
  min_lib$patient_id,
  min_lib$timepoint,
  "=",
  min_lib$observed_library_size,
  "\n"
)

cat(
  "HIGHEST LIBRARY               :",
  max_lib$patient_id,
  max_lib$timepoint,
  "=",
  max_lib$observed_library_size,
  "\n"
)

cat(
  "MAX/MIN LIBRARY RATIO         :",
  sprintf(
    "%.2f x",
    library_ratio
  ),
  "\n"
)

cat(
  "TECHNICAL REVIEW FLAGS        :",
  sum(
    meta$any_technical_review
  ),
  "/ 12\n"
)

cat(
  "SAMPLES EXCLUDED              : 0\n\n"
)

cat(
  "CORE GENE PRESENCE:\n"
)

for (
  i in seq_len(
    nrow(
      core_presence
    )
  )
) {

  cat(
    "  ",
    core_presence$gene[i],
    ": ",
    ifelse(
      core_presence$present[i],
      "PRESENT",
      "MISSING"
    ),
    "\n",
    sep = ""
  )
}

cat("\n")

cat(
  "CORE GENE EXPRESSION INSPECTED: NO\n\n"
)

cat(
  "FROZEN STEP19L SCRIPT         : FOUND\n"
)

cat(
  "FROZEN STEP19L SHA256         :",
  step19_sha,
  "\n"
)

cat(
  "SCORING SOURCE LINES SAVED    :",
  nrow(
    source_audit
  ),
  "\n\n"
)

cat(
  "SCORING IMPLEMENTATION TERMS:\n"
)

print(
  implementation_audit[
    detected == TRUE
  ]
)

cat("\n")

cat(
  "NORMALIZATION APPLIED TO OMIX : NO\n"
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
  "EASC PATIENT                  : UNRESOLVED\n"
)

cat(
  "PRIMARY CORE2 INFERENCE       : BLOCKED\n\n"
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
  "STEP22A-5E = LOCK THE EXACT FROZEN ",
  "NORMALIZATION / SCORE IMPLEMENTATION\n"
)

cat(
  "FROM STEP19L SOURCE WITHOUT CALCULATING ",
  "EXTERNAL CORE2 VALUES.\n\n"
)

cat(
  "STOP HERE.\n"
)

cat(
  "DO NOT SCORE CORE2 YET.\n\n"
)
