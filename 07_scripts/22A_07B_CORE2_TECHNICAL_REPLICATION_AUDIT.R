# ============================================================
# ESCC Neoadjuvant Project
# STEP22A-7B
#
# INDEPENDENT TECHNICAL REPLICATION AUDIT OF STEP22A-7A CORE2
#
# PURPOSE:
#   Determine whether the uniformly negative secondary CORE2
#   result could be explained by a computational / mapping error.
#
# AUDITS:
#   1. primary gate still CLOSED
#   2. Before/After sample-name consistency
#   3. P15 Before never enters technical reference
#   4. scenario-specific reference membership
#   5. library-size CPM
#   6. log2(CPM+1)
#   7. per-gene z-score
#   8. frozen CORE2_STRONG formula
#   9. After - Before delta orientation
#  10. independent results == Step22A-7A results
#
# ALSO decomposes CORE2 delta into:
#   -CDKN1B_z / 2
#   -GSN_z / 2
#
# This decomposition is TECHNICAL ONLY.
# No biological conclusion is made here.
#
# CORE4 is NOT interpreted in this step.
# ============================================================

options(
  stringsAsFactors = FALSE
)

suppressPackageStartupMessages({
  library(data.table)
})

cat("\n")
cat("====================================================================\n")
cat("STEP22A-7B: CORE2 TECHNICAL REPLICATION AUDIT\n")
cat("====================================================================\n\n")

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

TABLE_DIR <- file.path(
  PROJECT,
  "06_tables",
  "external_validation",
  "Step22A"
)

RESULT_DIR <- file.path(
  PROJECT,
  "04_results",
  "external_validation",
  "Step22A"
)

LOG_DIR <- file.path(
  PROJECT,
  "08_logs"
)

dir.create(
  RESULT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

RAW_BUNDLE <- file.path(
  OBJECT_DIR,
  "OMIX005710_ALL27_FIBROBLAST_RAW_PSEUDOBULK_HISTOLOGY_UNRESOLVED.rds"
)

GATE_FILE <- file.path(
  META_DIR,
  "STEP22A_ANALYSIS_GATE_STATUS.csv"
)

SCENARIO_FILE <- file.path(
  META_DIR,
  "STEP22A_06E_HISTOLOGY_SENSITIVITY_TECHNICAL_INTERSECTION.csv"
)

PAIR_LOCK_FILE <- file.path(
  META_DIR,
  "STEP22A_06E_FINAL_PAIRED_TECHNICAL_LOCK.csv"
)

SAMPLE_META_FILE <- file.path(
  META_DIR,
  "STEP22A_06D_ALL27_FIBROBLAST_SAMPLE_PSEUDOBULK_METADATA.csv"
)

PT_META_FILE <- file.path(
  META_DIR,
  "STEP22A_06D_ALL27_FIBROBLAST_PATIENT_TIMEPOINT_METADATA.csv"
)

GOVERNANCE_FILE <- file.path(
  META_DIR,
  "STEP22A_06F_FINAL_GOVERNANCE_AUDIT.csv"
)

FROZEN_INV <- file.path(
  META_DIR,
  "STEP22A_PREANALYSIS_FROZEN_FILE_INVENTORY.csv"
)

SCORES_7A_FILE <- file.path(
  TABLE_DIR,
  "STEP22A_07A_SECONDARY_SCENARIO_CORE_SCORES.csv"
)

DELTAS_7A_FILE <- file.path(
  TABLE_DIR,
  "STEP22A_07A_SECONDARY_PAIRED_DELTAS.csv"
)

RESULTS_7A_FILE <- file.path(
  TABLE_DIR,
  "STEP22A_07A_SECONDARY_PAIRED_SCENARIO_RESULTS.csv"
)

AUDIT_OUT <- file.path(
  META_DIR,
  "STEP22A_07B_CORE2_TECHNICAL_REPLICATION_AUDIT.csv"
)

TIMEPOINT_OUT <- file.path(
  META_DIR,
  "STEP22A_07B_TIMEPOINT_LABEL_AUDIT.csv"
)

COMPONENT_OUT <- file.path(
  TABLE_DIR,
  "STEP22A_07B_CORE2_PAIRWISE_COMPONENT_AUDIT.csv"
)

SCENARIO_RECALC_OUT <- file.path(
  TABLE_DIR,
  "STEP22A_07B_CORE2_INDEPENDENT_SCENARIO_RECALC.csv"
)

SUMMARY_OUT <- file.path(
  RESULT_DIR,
  "STEP22A_07B_CORE2_TECHNICAL_AUDIT_SUMMARY.txt"
)


# ============================================================
# Helpers
# ============================================================

safe_sd <- function(x) {

  z <- sd(x)

  if (
    !is.finite(z) ||
    z < 1e-10
  ) {

    z <- 1e-10
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


audit_rows <- list()


add_audit <- function(
  check,
  status,
  observed,
  expected,
  severity = "HARD"
) {

  audit_rows[[
    length(audit_rows) + 1
  ]] <<- data.table(
    check = check,
    status = status,
    observed = as.character(
      observed
    ),
    expected = as.character(
      expected
    ),
    severity = severity
  )
}


# ============================================================
# 1. Governance
# ============================================================

cat("===== 1. GOVERNANCE =====\n\n")

gate <- fread(
  GATE_FILE,
  encoding = "UTF-8"
)

primary_gate <- gate[
  gate ==
    "paired_CORE2_test",
  allowed
]

if (length(primary_gate) != 1) {

  stop(
    "Cannot uniquely resolve paired_CORE2_test gate"
  )
}

cat(
  "Primary CORE2 gate:",
  primary_gate,
  "\n"
)

if (primary_gate == "YES") {

  stop(
    "PRIMARY_CORE2_GATE_UNEXPECTEDLY_OPEN"
  )
}

gov <- fread(
  GOVERNANCE_FILE,
  encoding = "UTF-8"
)

gov_fail <- sum(
  gov$status != "PASS"
)

cat(
  "Step22A-6F governance failures:",
  gov_fail,
  "\n\n"
)

if (gov_fail != 0) {

  stop(
    "PREVIOUS_GOVERNANCE_AUDIT_FAILED"
  )
}

add_audit(
  "primary_gate_still_closed",
  "PASS",
  primary_gate,
  "NOT YES"
)

add_audit(
  "governance_6F_failures",
  "PASS",
  gov_fail,
  0
)


# ============================================================
# 2. Independent timepoint-label audit
# ============================================================

cat("===== 2. TIMEPOINT LABEL AUDIT =====\n\n")

sample_meta <- fread(
  SAMPLE_META_FILE,
  encoding = "UTF-8"
)

sample_meta[
  ,
  suffix_timepoint :=
    fifelse(
      grepl(
        "_B$",
        sample_name
      ),
      "Before",
      fifelse(
        grepl(
          "_A$",
          sample_name
        ),
        "After",
        "UNCLASSIFIED"
      )
    )
]

sample_meta[
  ,
  suffix_contradiction :=
    suffix_timepoint !=
      "UNCLASSIFIED" &
    suffix_timepoint !=
      timepoint
]

timepoint_contradictions <- sum(
  sample_meta$
    suffix_contradiction
)

timepoint_unclassified <- sum(
  sample_meta$
    suffix_timepoint ==
    "UNCLASSIFIED"
)

fwrite(
  sample_meta[
    ,
    .(
      patient_id,
      timepoint,
      sample_name,
      suffix_timepoint,
      suffix_contradiction
    )
  ],
  TIMEPOINT_OUT,
  bom = TRUE
)

cat(
  "Sample-name/timepoint contradictions:",
  timepoint_contradictions,
  "\n"
)

cat(
  "Unclassified sample-name suffixes:",
  timepoint_unclassified,
  "\n\n"
)

if (
  timepoint_contradictions > 0
) {

  print(
    sample_meta[
      suffix_contradiction ==
        TRUE
    ]
  )

  stop(
    "TIMEPOINT_LABEL_CONTRADICTION"
  )
}

add_audit(
  "sample_name_timepoint_contradictions",
  "PASS",
  timepoint_contradictions,
  0
)


# ============================================================
# 3. Load frozen raw pseudobulk structure
# ============================================================

cat("===== 3. LOAD RAW PSEUDOBULK =====\n\n")

bundle <- readRDS(
  RAW_BUNDLE
)

counts <- bundle$
  raw_counts_all_patient_timepoints

pt_meta_bundle <- as.data.table(
  bundle$
    patient_timepoint_metadata
)

features <- as.data.table(
  bundle$
    feature_metadata
)

pt_meta_file <- fread(
  PT_META_FILE,
  encoding = "UTF-8"
)

if (
  nrow(counts) !=
    nrow(features)
) {

  stop(
    "FEATURE_MATRIX_DIMENSION_MISMATCH"
  )
}

if (
  ncol(counts) !=
    nrow(pt_meta_bundle)
) {

  stop(
    "PT_MATRIX_DIMENSION_MISMATCH"
  )
}

if (
  !identical(
    colnames(counts),
    pt_meta_bundle$
      patient_timepoint
  )
) {

  stop(
    "PSEUDOBULK_COLUMN_ORDER_MISMATCH"
  )
}

if (
  !identical(
    pt_meta_bundle$
      patient_timepoint,
    pt_meta_file$
      patient_timepoint
  )
) {

  stop(
    "BUNDLE_VS_CSV_PATIENT_TIMEPOINT_ORDER_MISMATCH"
  )
}

cat(
  "Raw pseudobulk:",
  nrow(counts),
  "genes x",
  ncol(counts),
  "patient-timepoints\n\n"
)

add_audit(
  "raw_pseudobulk_dimensions",
  "PASS",
  paste(
    nrow(counts),
    ncol(counts),
    sep = "x"
  ),
  "16790x27"
)


# ============================================================
# 4. CORE2 gene mapping
# ============================================================

cat("===== 4. CORE2 GENE MAPPING =====\n\n")

features[
  ,
  symbol_upper :=
    toupper(
      gene_symbol
    )
]

core2_genes <- c(
  "CDKN1B",
  "GSN"
)

core_idx <- integer(
  length(
    core2_genes
  )
)

names(
  core_idx
) <- core2_genes

for (
  g in core2_genes
) {

  hit <- which(
    features$
      symbol_upper ==
      g
  )

  if (
    length(hit) != 1
  ) {

    stop(
      "CORE2_GENE_MAPPING_NOT_UNIQUE: ",
      g
    )
  }

  core_idx[
    g
  ] <- hit
}

cat(
  "CDKN1B ENSG:",
  features$gene_id[
    core_idx[
      "CDKN1B"
    ]
  ],
  "\n"
)

cat(
  "GSN ENSG:",
  features$gene_id[
    core_idx[
      "GSN"
    ]
  ],
  "\n\n"
)

add_audit(
  "CORE2_gene_mapping",
  "PASS",
  "CDKN1B=1;GSN=1",
  "CDKN1B=1;GSN=1"
)


# ============================================================
# 5. P15 firewall
# ============================================================

cat("===== 5. P15 BEFORE FIREWALL =====\n\n")

technical_pass <- (
  pt_meta_bundle$
    frozen_technical_status ==
    "PASS_GE20"
)

p15_before_idx <- which(
  pt_meta_bundle$patient_id ==
    "P15" &
  pt_meta_bundle$timepoint ==
    "Before"
)

if (
  length(
    p15_before_idx
  ) != 1
) {

  stop(
    "P15_BEFORE_NOT_UNIQUELY_IDENTIFIED"
  )
}

p15_before_pass <- technical_pass[
  p15_before_idx
]

p15_before_cells <- (
  pt_meta_bundle$
    n_fibroblast_cells[
      p15_before_idx
    ]
)

cat(
  "P15 Before fibroblasts:",
  p15_before_cells,
  "\n"
)

cat(
  "P15 Before technical pass:",
  p15_before_pass,
  "\n\n"
)

if (
  p15_before_pass
) {

  stop(
    "P15_BEFORE_UNEXPECTEDLY_IN_TECHNICAL_REFERENCE"
  )
}

add_audit(
  "P15_Before_reference_exclusion",
  "PASS",
  paste0(
    "cells=",
    p15_before_cells,
    ";pass=",
    p15_before_pass
  ),
  "cells=18;pass=FALSE"
)


# ============================================================
# 6. Frozen paired technical set
# ============================================================

cat("===== 6. PAIRED TECHNICAL SET =====\n\n")

pair_lock <- fread(
  PAIR_LOCK_FILE,
  encoding = "UTF-8"
)

baseline_pairs <- pair_lock[
  paired_candidate_status ==
    "PAIRED_TECHNICAL_CANDIDATE",
  patient_id
]

baseline_pairs <- baseline_pairs[
  order(
    patient_number(
      baseline_pairs
    )
  )
]

expected_pairs <- c(
  "P6",
  "P9",
  "P16",
  "P18",
  "P19"
)

cat(
  "Frozen paired IDs:",
  paste(
    baseline_pairs,
    collapse = ", "
  ),
  "\n\n"
)

if (
  !identical(
    baseline_pairs,
    expected_pairs
  )
) {

  stop(
    "PAIRED_TECHNICAL_SET_CHANGED"
  )
}

add_audit(
  "paired_technical_membership",
  "PASS",
  paste(
    baseline_pairs,
    collapse = ";"
  ),
  "P6;P9;P16;P18;P19"
)


# ============================================================
# 7. Load Step22A-7A output
# ============================================================

cat("===== 7. LOAD STEP22A-7A OUTPUTS =====\n\n")

scores_7a <- fread(
  SCORES_7A_FILE,
  encoding = "UTF-8"
)

deltas_7a <- fread(
  DELTAS_7A_FILE,
  encoding = "UTF-8"
)

results_7a <- fread(
  RESULTS_7A_FILE,
  encoding = "UTF-8"
)

scenarios <- fread(
  SCENARIO_FILE,
  encoding = "UTF-8"
)

if (
  nrow(
    scenarios
  ) != 23
) {

  stop(
    "SCENARIO_COUNT_NOT_23"
  )
}

cat(
  "Frozen scenarios:",
  nrow(
    scenarios
  ),
  "\n"
)

cat(
  "Step22A-7A score rows:",
  nrow(
    scores_7a
  ),
  "\n"
)

cat(
  "Step22A-7A paired delta rows:",
  nrow(
    deltas_7a
  ),
  "\n\n"
)


# ============================================================
# 8. Independent recomputation across all 23 scenarios
# ============================================================

cat("===== 8. INDEPENDENT CORE2 RECOMPUTATION =====\n\n")

scenario_audit <- list()
component_rows <- list()
independent_results <- list()

max_score_diff_global <- 0
max_z_diff_global <- 0
max_log2cpm_diff_global <- 0
max_delta_diff_global <- 0

for (
  i in seq_len(
    nrow(scenarios)
  )
) {

  s <- scenarios[i]

  sid <- as.character(
    s$scenario_id
  )

  hyp <- as.character(
    s$hypothetical_EASC_patient
  )

  if (
    is.na(hyp) ||
    hyp == ""
  ) {

    hyp <- "UNKNOWN"
  }


  # ----------------------------------------------------------
  # Independent reference construction
  # ----------------------------------------------------------

  ref_keep <- technical_pass

  if (
    hyp != "UNKNOWN"
  ) {

    ref_keep <- (
      ref_keep &
      pt_meta_bundle$
        patient_id != hyp
    )
  }

  ref_idx <- which(
    ref_keep
  )

  expected_ref_n <- as.integer(
    s$technical_reference_n
  )

  if (
    length(ref_idx) !=
      expected_ref_n
  ) {

    stop(
      "REFERENCE_SIZE_MISMATCH: ",
      sid
    )
  }


  # P15 Before must never enter any scenario.
  if (
    p15_before_idx %in%
      ref_idx
  ) {

    stop(
      "P15_BEFORE_ENTERED_REFERENCE: ",
      sid
    )
  }


  # ----------------------------------------------------------
  # Library size uses all 16,790 genes
  # ----------------------------------------------------------

  lib <- as.numeric(
    colSums(
      counts[
        ,
        ref_idx,
        drop = FALSE
      ]
    )
  )

  if (
    any(
      !is.finite(lib)
    ) ||
    any(
      lib <= 0
    )
  ) {

    stop(
      "INVALID_LIBRARY_SIZE: ",
      sid
    )
  }


  # ----------------------------------------------------------
  # Independent CORE2 transformation
  # ----------------------------------------------------------

  raw2 <- counts[
    core_idx,
    ref_idx,
    drop = FALSE
  ]

  l2c <- log2(
    sweep(
      raw2,
      2,
      lib,
      "/"
    ) *
      1e6 +
      1
  )

  rownames(
    l2c
  ) <- core2_genes


  mu <- rowMeans(
    l2c
  )

  sig <- apply(
    l2c,
    1,
    safe_sd
  )

  z <- sweep(
    l2c,
    1,
    mu,
    "-"
  )

  z <- sweep(
    z,
    1,
    sig,
    "/"
  )


  # Internal z-score sanity.
  z_mean_error <- max(
    abs(
      rowMeans(z)
    )
  )

  z_sd_error <- max(
    abs(
      apply(
        z,
        1,
        sd
      ) -
        1
    )
  )


  core2 <- (
    -z[
      "CDKN1B",
    ] -
    z[
      "GSN",
    ]
  ) / 2


  recalc <- data.table(
    scenario_id =
      sid,

    hypothetical_EASC_patient =
      hyp,

    patient_id =
      pt_meta_bundle$
        patient_id[
          ref_idx
        ],

    timepoint =
      pt_meta_bundle$
        timepoint[
          ref_idx
        ],

    patient_timepoint =
      pt_meta_bundle$
        patient_timepoint[
          ref_idx
        ],

    reference_n =
      length(
        ref_idx
      ),

    library_size_recalc =
      lib,

    CDKN1B_log2CPM_recalc =
      as.numeric(
        l2c[
          "CDKN1B",
        ]
      ),

    GSN_log2CPM_recalc =
      as.numeric(
        l2c[
          "GSN",
        ]
      ),

    CDKN1B_z_recalc =
      as.numeric(
        z[
          "CDKN1B",
        ]
      ),

    GSN_z_recalc =
      as.numeric(
        z[
          "GSN",
        ]
      ),

    CORE2_STRONG_recalc =
      as.numeric(
        core2
      )
  )


  # ----------------------------------------------------------
  # Compare against frozen Step22A-7A scores
  # ----------------------------------------------------------

  old <- scores_7a[
    scenario_id ==
      sid,
    .(
      patient_timepoint,
      library_size,
      CDKN1B_log2CPM,
      GSN_log2CPM,
      CDKN1B_z,
      GSN_z,
      CORE2_STRONG
    )
  ]


  cmp <- merge(
    recalc,
    old,
    by = "patient_timepoint",
    all = TRUE
  )


  if (
    nrow(cmp) !=
      length(ref_idx) ||
    anyNA(
      cmp$CORE2_STRONG
    ) ||
    anyNA(
      cmp$CORE2_STRONG_recalc
    )
  ) {

    stop(
      "7A_SCORE_MEMBERSHIP_MISMATCH: ",
      sid
    )
  }


  lib_diff <- max(
    abs(
      cmp$
        library_size_recalc -
      cmp$
        library_size
    )
  )

  l2c_diff <- max(
    abs(
      c(
        cmp$
          CDKN1B_log2CPM_recalc -
        cmp$
          CDKN1B_log2CPM,

        cmp$
          GSN_log2CPM_recalc -
        cmp$
          GSN_log2CPM
      )
    )
  )

  z_diff <- max(
    abs(
      c(
        cmp$
          CDKN1B_z_recalc -
        cmp$
          CDKN1B_z,

        cmp$
          GSN_z_recalc -
        cmp$
          GSN_z
      )
    )
  )

  score_diff <- max(
    abs(
      cmp$
        CORE2_STRONG_recalc -
      cmp$
        CORE2_STRONG
    )
  )


  max_log2cpm_diff_global <- max(
    max_log2cpm_diff_global,
    l2c_diff
  )

  max_z_diff_global <- max(
    max_z_diff_global,
    z_diff
  )

  max_score_diff_global <- max(
    max_score_diff_global,
    score_diff
  )


  # ----------------------------------------------------------
  # Independent paired delta
  # ----------------------------------------------------------

  paired_ids <- baseline_pairs

  if (
    hyp %in%
      paired_ids
  ) {

    paired_ids <- paired_ids[
      paired_ids != hyp
    ]
  }

  expected_pair_n <- as.integer(
    s$paired_candidate_n
  )

  if (
    length(
      paired_ids
    ) != expected_pair_n
  ) {

    stop(
      "PAIR_COUNT_MISMATCH: ",
      sid
    )
  }


  scenario_delta <- numeric()

  for (
    patient in paired_ids
  ) {

    b <- recalc[
      patient_id ==
        patient &
      timepoint ==
        "Before"
    ]

    a <- recalc[
      patient_id ==
        patient &
      timepoint ==
        "After"
    ]


    if (
      nrow(b) != 1 ||
      nrow(a) != 1
    ) {

      stop(
        "PAIR_MAPPING_FAILURE: ",
        sid,
        " / ",
        patient
      )
    }


    delta_cdkn1b_z <- (
      a$CDKN1B_z_recalc -
      b$CDKN1B_z_recalc
    )

    delta_gsn_z <- (
      a$GSN_z_recalc -
      b$GSN_z_recalc
    )

    contribution_cdkn1b <- (
      -delta_cdkn1b_z / 2
    )

    contribution_gsn <- (
      -delta_gsn_z / 2
    )

    delta_recalc <- (
      a$CORE2_STRONG_recalc -
      b$CORE2_STRONG_recalc
    )


    algebra_error <- abs(
      delta_recalc -
        (
          contribution_cdkn1b +
          contribution_gsn
        )
    )


    old_delta <- deltas_7a[
      scenario_id ==
        sid &
      patient_id ==
        patient,
      delta_CORE2
    ]


    if (
      length(
        old_delta
      ) != 1
    ) {

      stop(
        "7A_DELTA_NOT_UNIQUE: ",
        sid,
        " / ",
        patient
      )
    }


    delta_diff <- abs(
      delta_recalc -
        old_delta
    )


    max_delta_diff_global <- max(
      max_delta_diff_global,
      delta_diff
    )


    component_rows[[
      length(
        component_rows
      ) + 1
    ]] <- data.table(
      scenario_id =
        sid,

      hypothetical_EASC_patient =
        hyp,

      patient_id =
        patient,

      n_paired =
        length(
          paired_ids
        ),

      before_patient_timepoint =
        b$patient_timepoint,

      after_patient_timepoint =
        a$patient_timepoint,

      delta_definition =
        "AFTER_MINUS_BEFORE",

      delta_CDKN1B_z =
        delta_cdkn1b_z,

      contribution_minus_CDKN1B_over2 =
        contribution_cdkn1b,

      delta_GSN_z =
        delta_gsn_z,

      contribution_minus_GSN_over2 =
        contribution_gsn,

      CORE2_delta_recalc =
        delta_recalc,

      CORE2_delta_7A =
        old_delta,

      absolute_delta_difference =
        delta_diff,

      algebra_error =
        algebra_error
    )


    scenario_delta <- c(
      scenario_delta,
      delta_recalc
    )
  }


  independent_results[[
    length(
      independent_results
    ) + 1
  ]] <- data.table(
    scenario_id =
      sid,

    hypothetical_EASC_patient =
      hyp,

    reference_n =
      length(
        ref_idx
      ),

    paired_n =
      length(
        paired_ids
      ),

    paired_ids =
      paste(
        paired_ids,
        collapse = ";"
      ),

    mean_delta =
      mean(
        scenario_delta
      ),

    median_delta =
      median(
        scenario_delta
      ),

    positive_n =
      sum(
        scenario_delta > 0
      ),

    negative_n =
      sum(
        scenario_delta < 0
      ),

    positive_fraction =
      mean(
        scenario_delta > 0
      ),

    direction =
      ifelse(
        mean(
          scenario_delta
        ) > 0,
        "POSITIVE",
        ifelse(
          mean(
            scenario_delta
          ) < 0,
          "NEGATIVE",
          "ZERO"
        )
      )
  )


  scenario_audit[[
    length(
      scenario_audit
    ) + 1
  ]] <- data.table(
    scenario_id =
      sid,

    hypothetical_EASC_patient =
      hyp,

    reference_n =
      length(
        ref_idx
      ),

    paired_n =
      length(
        paired_ids
      ),

    library_size_max_abs_diff =
      lib_diff,

    log2CPM_max_abs_diff =
      l2c_diff,

    zscore_max_abs_diff =
      z_diff,

    CORE2_score_max_abs_diff =
      score_diff,

    z_mean_max_abs_error =
      z_mean_error,

    z_sd_max_abs_error =
      z_sd_error,

    P15_Before_in_reference =
      p15_before_idx %in%
        ref_idx
  )


  cat(
    sprintf(
      "[%02d/23] %-28s ref=%2d pair=%d score.diff=%.3g delta.mean=% .4f\n",
      i,
      sid,
      length(
        ref_idx
      ),
      length(
        paired_ids
      ),
      score_diff,
      mean(
        scenario_delta
      )
    )
  )
}

cat("\n")


scenario_audit_dt <- rbindlist(
  scenario_audit
)

component_dt <- rbindlist(
  component_rows
)

independent_dt <- rbindlist(
  independent_results
)


fwrite(
  scenario_audit_dt,
  AUDIT_OUT,
  bom = TRUE
)

fwrite(
  component_dt,
  COMPONENT_OUT,
  bom = TRUE
)

fwrite(
  independent_dt,
  SCENARIO_RECALC_OUT,
  bom = TRUE
)


# ============================================================
# 9. Global technical equivalence checks
# ============================================================

cat("===== 9. GLOBAL NUMERICAL EQUIVALENCE =====\n\n")

tolerance <- 1e-10

cat(
  "Max abs log2CPM difference :",
  format(
    max_log2cpm_diff_global,
    scientific = TRUE
  ),
  "\n"
)

cat(
  "Max abs z-score difference :",
  format(
    max_z_diff_global,
    scientific = TRUE
  ),
  "\n"
)

cat(
  "Max abs CORE2 difference   :",
  format(
    max_score_diff_global,
    scientific = TRUE
  ),
  "\n"
)

cat(
  "Max abs paired-delta diff  :",
  format(
    max_delta_diff_global,
    scientific = TRUE
  ),
  "\n\n"
)


numeric_pass <- (
  max_log2cpm_diff_global <=
    tolerance &&
  max_z_diff_global <=
    tolerance &&
  max_score_diff_global <=
    tolerance &&
  max_delta_diff_global <=
    tolerance
)


add_audit(
  "independent_recalculation_matches_7A",
  ifelse(
    numeric_pass,
    "PASS",
    "FAIL"
  ),
  paste0(
    "l2c=",
    signif(
      max_log2cpm_diff_global,
      5
    ),
    ";z=",
    signif(
      max_z_diff_global,
      5
    ),
    ";score=",
    signif(
      max_score_diff_global,
      5
    ),
    ";delta=",
    signif(
      max_delta_diff_global,
      5
    )
  ),
  "<=1e-10"
)


# ============================================================
# 10. Delta-orientation audit
# ============================================================

cat("===== 10. AFTER - BEFORE ORIENTATION =====\n\n")

algebra_max <- max(
  component_dt$
    algebra_error
)

delta_diff_max <- max(
  component_dt$
    absolute_delta_difference
)

cat(
  "CORE2 contribution algebra max error:",
  format(
    algebra_max,
    scientific = TRUE
  ),
  "\n"
)

cat(
  "Step22A-7A delta max difference:",
  format(
    delta_diff_max,
    scientific = TRUE
  ),
  "\n"
)

cat(
  "Delta definition:",
  unique(
    component_dt$
      delta_definition
  ),
  "\n\n"
)


if (
  algebra_max >
    tolerance
) {

  stop(
    "CORE2_COMPONENT_ALGEBRA_FAILURE"
  )
}


if (
  delta_diff_max >
    tolerance
) {

  stop(
    "AFTER_BEFORE_DELTA_MISMATCH"
  )
}


add_audit(
  "delta_orientation",
  "PASS",
  "AFTER_MINUS_BEFORE",
  "AFTER_MINUS_BEFORE"
)


# ============================================================
# 11. Independent direction result
# ============================================================

cat("===== 11. INDEPENDENT DIRECTION CHECK =====\n\n")

hyp_ind <- independent_dt[
  scenario_id !=
    "S00_UNRESOLVED"
]

negative_scenarios <- sum(
  hyp_ind$
    mean_delta < 0
)

median_negative_scenarios <- sum(
  hyp_ind$
    median_delta < 0
)

majority_negative_scenarios <- sum(
  hyp_ind$
    negative_n >
    hyp_ind$
      paired_n /
      2
)

all_patient_negative_scenarios <- sum(
  hyp_ind$
    negative_n ==
    hyp_ind$
      paired_n
)


cat(
  "Mean delta negative:",
  negative_scenarios,
  "/ 22\n"
)

cat(
  "Median delta negative:",
  median_negative_scenarios,
  "/ 22\n"
)

cat(
  "Majority patients negative:",
  majority_negative_scenarios,
  "/ 22\n"
)

cat(
  "All paired patients negative:",
  all_patient_negative_scenarios,
  "/ 22\n\n"
)


add_audit(
  "independent_22_scenario_direction",
  ifelse(
    negative_scenarios == 22,
    "PASS",
    "FAIL"
  ),
  paste0(
    negative_scenarios,
    "/22 negative"
  ),
  "22/22 negative"
)


# ============================================================
# 12. Compare scenario summary against Step22A-7A
# ============================================================

cat("===== 12. SCENARIO-SUMMARY CROSSCHECK =====\n\n")

old_core2 <- results_7a[
  endpoint ==
    "CORE2_STRONG",
  .(
    scenario_id,
    old_mean_delta =
      mean_delta,
    old_median_delta =
      median_delta,
    old_positive_n =
      positive_n,
    old_negative_n =
      negative_n
  )
]


summary_cmp <- merge(
  independent_dt,
  old_core2,
  by = "scenario_id"
)


summary_mean_diff <- max(
  abs(
    summary_cmp$
      mean_delta -
    summary_cmp$
      old_mean_delta
  )
)

summary_median_diff <- max(
  abs(
    summary_cmp$
      median_delta -
    summary_cmp$
      old_median_delta
  )
)

summary_count_mismatch <- sum(
  summary_cmp$
    positive_n !=
    summary_cmp$
      old_positive_n |
  summary_cmp$
    negative_n !=
    summary_cmp$
      old_negative_n
)


cat(
  "Mean summary max diff:",
  format(
    summary_mean_diff,
    scientific = TRUE
  ),
  "\n"
)

cat(
  "Median summary max diff:",
  format(
    summary_median_diff,
    scientific = TRUE
  ),
  "\n"
)

cat(
  "Direction-count mismatches:",
  summary_count_mismatch,
  "\n\n"
)


summary_pass <- (
  summary_mean_diff <=
    tolerance &&
  summary_median_diff <=
    tolerance &&
  summary_count_mismatch == 0
)


add_audit(
  "7A_scenario_summary_reproduced",
  ifelse(
    summary_pass,
    "PASS",
    "FAIL"
  ),
  paste0(
    "mean=",
    signif(
      summary_mean_diff,
      5
    ),
    ";median=",
    signif(
      summary_median_diff,
      5
    ),
    ";countMismatch=",
    summary_count_mismatch
  ),
  "all <=1e-10; countMismatch=0"
)


# ============================================================
# 13. Frozen Step19-21 integrity
# ============================================================

cat("===== 13. FROZEN PROJECT CHECK =====\n\n")

frozen <- fread(
  FROZEN_INV,
  encoding = "UTF-8"
)

changed <- character()


for (
  i in seq_len(
    nrow(frozen)
  )
) {

  path <- as.character(
    frozen$path[i]
  )

  if (!file.exists(path)) {

    changed <- c(
      changed,
      paste0(
        path,
        ":MISSING"
      )
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

    changed <- c(
      changed,
      paste0(
        path,
        ":SIZE_CHANGED"
      )
    )
  }
}


cat(
  "Frozen original files changed:",
  length(
    changed
  ),
  "\n\n"
)


if (
  length(
    changed
  ) > 0
) {

  print(
    changed
  )

  stop(
    "FROZEN_PROJECT_CHANGED"
  )
}


add_audit(
  "frozen_Step19_21_integrity",
  "PASS",
  length(
    changed
  ),
  0
)


# ============================================================
# 14. Final audit result
# ============================================================

final_audit <- rbindlist(
  audit_rows,
  use.names = TRUE
)

hard_failures <- final_audit[
  severity ==
    "HARD" &
  status ==
    "FAIL"
]


cat("===== 14. FINAL TECHNICAL AUDIT =====\n\n")

cat(
  "Governance/technical checks:",
  nrow(
    final_audit
  ),
  "\n"
)

cat(
  "Hard failures:",
  nrow(
    hard_failures
  ),
  "\n\n"
)


if (
  nrow(
    hard_failures
  ) > 0
) {

  print(
    hard_failures
  )

  stop(
    "STEP22A_07B_TECHNICAL_AUDIT_FAILED"
  )
}


technical_status <- (
  "TECHNICAL_REPLICATION_AUDIT_PASS_NO_SIGN_PAIR_REFERENCE_ERROR_DETECTED"
)


# ============================================================
# 15. Summary file
# ============================================================

summary_lines <- c(
  "STEP22A-7B CORE2 TECHNICAL REPLICATION AUDIT",
  "",
  paste0(
    "Technical audit failures: ",
    nrow(
      hard_failures
    )
  ),
  paste0(
    "Sample-name/timepoint contradictions: ",
    timepoint_contradictions
  ),
  paste0(
    "P15 Before fibroblast cells: ",
    p15_before_cells
  ),
  paste0(
    "P15 Before in technical reference: ",
    p15_before_pass
  ),
  paste0(
    "Max abs log2CPM difference vs 7A: ",
    format(
      max_log2cpm_diff_global,
      scientific = TRUE
    )
  ),
  paste0(
    "Max abs z-score difference vs 7A: ",
    format(
      max_z_diff_global,
      scientific = TRUE
    )
  ),
  paste0(
    "Max abs CORE2 score difference vs 7A: ",
    format(
      max_score_diff_global,
      scientific = TRUE
    )
  ),
  paste0(
    "Max abs paired delta difference vs 7A: ",
    format(
      max_delta_diff_global,
      scientific = TRUE
    )
  ),
  "Delta definition: After - Before",
  "",
  paste0(
    "Independent mean-negative scenarios: ",
    negative_scenarios,
    "/22"
  ),
  paste0(
    "Independent median-negative scenarios: ",
    median_negative_scenarios,
    "/22"
  ),
  paste0(
    "Independent majority-negative scenarios: ",
    majority_negative_scenarios,
    "/22"
  ),
  paste0(
    "Independent all-patient-negative scenarios: ",
    all_patient_negative_scenarios,
    "/22"
  ),
  "",
  paste0(
    "Technical status: ",
    technical_status
  ),
  "",
  "INTERPRETATION BOUNDARY:",
  "This audit establishes computational reproducibility only.",
  "It does NOT establish a biological explanation.",
  "Primary ESCC-only analysis remains blocked.",
  "CORE4 has not been used to rescue or redefine CORE2."
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
cat("STEP22A-7B COMPLETE\n")
cat("====================================================================\n\n")

cat(
  "TECHNICAL AUDIT FAILURES      : 0\n"
)

cat(
  "TIMEPOINT CONTRADICTIONS      :",
  timepoint_contradictions,
  "\n"
)

cat(
  "P15 BEFORE IN REFERENCE       : NO\n"
)

cat(
  "DELTA DEFINITION              : AFTER - BEFORE\n\n"
)

cat(
  "7A log2CPM MAX ABS DIFF       :",
  format(
    max_log2cpm_diff_global,
    scientific = TRUE
  ),
  "\n"
)

cat(
  "7A Z-SCORE MAX ABS DIFF       :",
  format(
    max_z_diff_global,
    scientific = TRUE
  ),
  "\n"
)

cat(
  "7A CORE2 MAX ABS DIFF         :",
  format(
    max_score_diff_global,
    scientific = TRUE
  ),
  "\n"
)

cat(
  "7A DELTA MAX ABS DIFF         :",
  format(
    max_delta_diff_global,
    scientific = TRUE
  ),
  "\n\n"
)

cat(
  "INDEPENDENT MEAN NEGATIVE     :",
  negative_scenarios,
  "/ 22\n"
)

cat(
  "INDEPENDENT MEDIAN NEGATIVE   :",
  median_negative_scenarios,
  "/ 22\n"
)

cat(
  "MAJORITY PATIENTS NEGATIVE    :",
  majority_negative_scenarios,
  "/ 22\n"
)

cat(
  "ALL PATIENTS NEGATIVE         :",
  all_patient_negative_scenarios,
  "/ 22\n\n"
)

cat(
  "TECHNICAL STATUS:\n"
)

cat(
  technical_status,
  "\n\n"
)

cat(
  "INTERPRETATION:\n"
)

cat(
  "  No computational sign / pair / reference-membership\n"
)

cat(
  "  error was detected in the CORE2 secondary result.\n"
)

cat(
  "  This does NOT by itself prove a biological explanation.\n\n"
)

cat(
  "PRIMARY ESCC-ONLY             : STILL BLOCKED\n"
)

cat(
  "CORE4 INTERPRETATION          : NOT YET PERFORMED\n"
)

cat(
  "MANUSCRIPT MODIFIED           : NO\n\n"
)

cat(
  "NEXT IF THIS PASSES:\n"
)

cat(
  "STEP22A-7C = PRE-SPECIFIED COMPARATOR AUDIT:\n"
)

cat(
  "CORE4 + CDKN1B/GSN CONTRIBUTIONS + PAIRED LOPO.\n\n"
)
