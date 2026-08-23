# ============================================================
# ESCC Neoadjuvant Project
# STEP22A-7C
#
# PRE-SPECIFIED COMPARATOR / COMPONENT / LOPO AUDIT
#
# Questions:
#
# 1. Is the negative CORE2_STRONG result specific to the
#    CDKN1B/GSN arm, or is the four-gene CORE module also
#    directionally discordant?
#
# 2. Which arithmetic components contribute to CORE2 delta?
#
# 3. Is the negative CORE2 result driven by one paired patient?
#
# Analyses:
#
# A. Independent CORE4 re-calculation in all 22 frozen
#    hypothetical-EASC scenarios.
#
# B. Decomposition:
#
#    CORE2_STRONG =
#      (-CDKN1B_z - GSN_z) / 2
#
#    therefore paired delta contribution:
#
#      -delta(CDKN1B_z)/2
#      -delta(GSN_z)/2
#
#    Positive-arm comparator:
#
#      BGN_TIMP1_component =
#        (BGN_z + TIMP1_z) / 2
#
#    and:
#
#      CORE4 =
#        (CORE2_STRONG + BGN_TIMP1_component) / 2
#
# C. Full leave-one-paired-patient-out (LOPO):
#    the left-out patient is removed BOTH from:
#      - z-score reference
#      - paired effect calculation
#
# Strict boundaries:
#
# - PRIMARY ESCC-only analysis remains BLOCKED.
# - No scenario selection.
# - No threshold modification.
# - No refitting of weights.
# - No new marker selection.
# - LOPO is influence analysis, not a new endpoint.
# - CORE4 cannot "rescue" CORE2.
# ============================================================

options(
  stringsAsFactors = FALSE
)

suppressPackageStartupMessages({
  library(data.table)
})

set.seed(20260818)

cat("\n")
cat("====================================================================\n")
cat("STEP22A-7C: CORE4 + COMPONENT + PAIRED LOPO AUDIT\n")
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
  TABLE_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  RESULT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

RAW_BUNDLE <- file.path(
  OBJECT_DIR,
  paste0(
    "OMIX005710_ALL27_FIBROBLAST_RAW_PSEUDOBULK_",
    "HISTOLOGY_UNRESOLVED.rds"
  )
)

GATE_FILE <- file.path(
  META_DIR,
  "STEP22A_ANALYSIS_GATE_STATUS.csv"
)

SECONDARY_AUTH <- file.path(
  META_DIR,
  "STEP22A_07A_SECONDARY_SENSITIVITY_AUTHORIZATION.csv"
)

SCENARIO_FILE <- file.path(
  META_DIR,
  "STEP22A_06E_HISTOLOGY_SENSITIVITY_TECHNICAL_INTERSECTION.csv"
)

PAIR_LOCK_FILE <- file.path(
  META_DIR,
  "STEP22A_06E_FINAL_PAIRED_TECHNICAL_LOCK.csv"
)

FROZEN_INV <- file.path(
  META_DIR,
  "STEP22A_PREANALYSIS_FROZEN_FILE_INVENTORY.csv"
)

# Step22A-7A
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

# Step22A-7B
AUDIT_7B_FILE <- file.path(
  META_DIR,
  "STEP22A_07B_CORE2_TECHNICAL_REPLICATION_AUDIT.csv"
)

TIMEPOINT_7B_FILE <- file.path(
  META_DIR,
  "STEP22A_07B_TIMEPOINT_LABEL_AUDIT.csv"
)

COMPONENT_7B_FILE <- file.path(
  TABLE_DIR,
  "STEP22A_07B_CORE2_PAIRWISE_COMPONENT_AUDIT.csv"
)

INDEPENDENT_7B_FILE <- file.path(
  TABLE_DIR,
  "STEP22A_07B_CORE2_INDEPENDENT_SCENARIO_RECALC.csv"
)

# Outputs
CORE4_AUDIT_OUT <- file.path(
  META_DIR,
  "STEP22A_07C_CORE4_INDEPENDENT_REPLICATION_AUDIT.csv"
)

SCENARIO_OUT <- file.path(
  TABLE_DIR,
  "STEP22A_07C_SCENARIO_COMPARATOR_SUMMARY.csv"
)

COMPONENT_OUT <- file.path(
  TABLE_DIR,
  "STEP22A_07C_COMPONENT_SCENARIO_SUMMARY.csv"
)

LOPO_OUT <- file.path(
  TABLE_DIR,
  "STEP22A_07C_PAIRED_LOPO_RESULTS.csv"
)

ROBUSTNESS_OUT <- file.path(
  TABLE_DIR,
  "STEP22A_07C_GLOBAL_ROBUSTNESS_SUMMARY.csv"
)

SUMMARY_OUT <- file.path(
  RESULT_DIR,
  "STEP22A_07C_CORE4_COMPONENT_LOPO_SUMMARY.txt"
)

# ============================================================
# Helpers
# ============================================================

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


direction_label <- function(x) {

  if (
    !is.finite(x)
  ) {

    return("NA")
  }

  if (x > 0) {

    return("POSITIVE")

  } else if (x < 0) {

    return("NEGATIVE")

  } else {

    return("ZERO")
  }
}


resolve_project_path <- function(x) {

  if (
    grepl(
      "^/",
      x
    )
  ) {

    return(x)
  }

  file.path(
    PROJECT,
    x
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

manuscript_gate <- gate[
  gate ==
    "manuscript_update",
  allowed
]

if (
  length(primary_gate) != 1 ||
  length(manuscript_gate) != 1
) {

  stop(
    "Cannot uniquely resolve governance gates"
  )
}

cat(
  "Primary CORE2 gate :",
  primary_gate,
  "\n"
)

cat(
  "Manuscript gate    :",
  manuscript_gate,
  "\n"
)

if (
  primary_gate == "YES"
) {

  stop(
    "PRIMARY_CORE2_GATE_UNEXPECTEDLY_OPEN"
  )
}

if (
  manuscript_gate == "YES"
) {

  stop(
    "MANUSCRIPT_GATE_UNEXPECTEDLY_OPEN"
  )
}

auth <- fread(
  SECONDARY_AUTH,
  encoding = "UTF-8"
)

secondary_status <- auth[
  item ==
    "SECONDARY_HISTOLOGY_SENSITIVITY",
  status
]

scenario_selection <- auth[
  item ==
    "SCENARIO_SELECTION_ALLOWED",
  status
]

if (
  length(secondary_status) != 1 ||
  secondary_status != "AUTHORIZED"
) {

  stop(
    "SECONDARY_SENSITIVITY_NOT_AUTHORIZED"
  )
}

if (
  length(scenario_selection) != 1 ||
  scenario_selection != "NO"
) {

  stop(
    "SCENARIO_SELECTION_FIREWALL_NOT_LOCKED"
  )
}

cat(
  "Secondary sensitivity:",
  secondary_status,
  "\n"
)

cat(
  "Scenario selection    :",
  scenario_selection,
  "\n\n"
)

cat(
  "✓ Primary remains blocked.\n"
)

cat(
  "✓ Secondary comparator audit allowed.\n\n"
)

# ============================================================
# 2. Verify Step22A-7B technical audit artifacts
# ============================================================

cat("===== 2. VERIFY STEP22A-7B TECHNICAL AUDIT =====\n\n")

audit7b <- fread(
  AUDIT_7B_FILE,
  encoding = "UTF-8"
)

time7b <- fread(
  TIMEPOINT_7B_FILE,
  encoding = "UTF-8"
)

component7b <- fread(
  COMPONENT_7B_FILE,
  encoding = "UTF-8"
)

independent7b <- fread(
  INDEPENDENT_7B_FILE,
  encoding = "UTF-8"
)

tol <- 1e-10

audit_numeric_pass <- (
  max(
    audit7b$log2CPM_max_abs_diff,
    na.rm = TRUE
  ) <= tol &&
  max(
    audit7b$zscore_max_abs_diff,
    na.rm = TRUE
  ) <= tol &&
  max(
    audit7b$CORE2_score_max_abs_diff,
    na.rm = TRUE
  ) <= tol
)

p15_reference_count <- sum(
  audit7b$P15_Before_in_reference %in%
    c(
      TRUE,
      "TRUE",
      "T",
      1,
      "1"
    )
)

timepoint_contradictions <- sum(
  time7b$suffix_contradiction %in%
    c(
      TRUE,
      "TRUE",
      "T",
      1,
      "1"
    )
)

component_algebra_max <- max(
  component7b$algebra_error,
  na.rm = TRUE
)

component_delta_diff_max <- max(
  component7b$absolute_delta_difference,
  na.rm = TRUE
)

hyp_independent7b <- independent7b[
  grepl(
    "^S_HYPOTHETICAL_EASC_",
    scenario_id
  )
]

negative7b <- sum(
  hyp_independent7b$mean_delta < 0
)

cat(
  "7B numerical audit pass     :",
  audit_numeric_pass,
  "\n"
)

cat(
  "P15 reference violations    :",
  p15_reference_count,
  "\n"
)

cat(
  "Timepoint contradictions    :",
  timepoint_contradictions,
  "\n"
)

cat(
  "Component algebra max error :",
  format(
    component_algebra_max,
    scientific = TRUE
  ),
  "\n"
)

cat(
  "Paired-delta max difference :",
  format(
    component_delta_diff_max,
    scientific = TRUE
  ),
  "\n"
)

cat(
  "Independent negative CORE2  :",
  negative7b,
  "/22\n\n"
)

if (
  !audit_numeric_pass ||
  p15_reference_count != 0 ||
  timepoint_contradictions != 0 ||
  component_algebra_max > tol ||
  component_delta_diff_max > tol ||
  negative7b != 22
) {

  stop(
    "STEP22A_07B_TECHNICAL_AUDIT_NOT_CLEAN"
  )
}

cat(
  "✓ Step22A-7B technical replication verified.\n\n"
)

# ============================================================
# 3. Load frozen scenarios + paired set
# ============================================================

cat("===== 3. LOAD FROZEN SCENARIOS / PAIRED SET =====\n\n")

scenario_dt <- fread(
  SCENARIO_FILE,
  encoding = "UTF-8"
)

scenario_dt <- scenario_dt[
  grepl(
    "^S_HYPOTHETICAL_EASC_",
    scenario_id
  )
]

if (
  nrow(
    scenario_dt
  ) != 22
) {

  stop(
    "EXPECTED_22_HYPOTHETICAL_EASC_SCENARIOS"
  )
}

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

if (
  !identical(
    baseline_pairs,
    expected_pairs
  )
) {

  stop(
    "FROZEN_PAIRED_SET_CHANGED"
  )
}

cat(
  "Frozen hypothetical EASC scenarios:",
  nrow(
    scenario_dt
  ),
  "\n"
)

cat(
  "Baseline paired technical IDs:",
  paste(
    baseline_pairs,
    collapse = ", "
  ),
  "\n\n"
)

# ============================================================
# 4. Load raw pseudobulk
# ============================================================

cat("===== 4. LOAD RAW PSEUDOBULK =====\n\n")

bundle <- readRDS(
  RAW_BUNDLE
)

counts <- bundle$
  raw_counts_all_patient_timepoints

pt_meta <- as.data.table(
  bundle$
    patient_timepoint_metadata
)

features <- as.data.table(
  bundle$
    feature_metadata
)

if (
  nrow(counts) !=
    nrow(features) ||
  ncol(counts) !=
    nrow(pt_meta)
) {

  stop(
    "RAW_PSEUDOBULK_DIMENSION_MISMATCH"
  )
}

if (
  !identical(
    colnames(counts),
    pt_meta$patient_timepoint
  )
) {

  stop(
    "RAW_PSEUDOBULK_METADATA_ORDER_MISMATCH"
  )
}

technical_pass <- (
  pt_meta$frozen_technical_status ==
    "PASS_GE20"
)

if (
  sum(
    technical_pass
  ) != 26
) {

  stop(
    "EXPECTED_26_TECHNICAL_PASS_PATIENT_TIMEPOINTS"
  )
}

cat(
  "Raw pseudobulk:",
  nrow(counts),
  "genes x",
  ncol(counts),
  "patient-timepoints\n"
)

cat(
  "Technical PASS:",
  sum(
    technical_pass
  ),
  "/27\n\n"
)

# ============================================================
# 5. Resolve frozen four-gene CORE module
# ============================================================

cat("===== 5. FOUR-GENE CORE MAPPING =====\n\n")

features[
  ,
  symbol_upper :=
    toupper(
      gene_symbol
    )
]

core_genes <- c(
  "BGN",
  "CDKN1B",
  "GSN",
  "TIMP1"
)

core_index <- integer(
  length(
    core_genes
  )
)

names(
  core_index
) <- core_genes

for (
  g in core_genes
) {

  hit <- which(
    features$symbol_upper ==
      g
  )

  if (
    length(hit) != 1
  ) {

    stop(
      "CORE_GENE_MAPPING_NOT_UNIQUE: ",
      g
    )
  }

  core_index[g] <- hit

  cat(
    sprintf(
      "%-7s -> %s\n",
      g,
      features$gene_id[
        hit
      ]
    )
  )
}

cat("\n")

# ============================================================
# 6. Frozen scoring function
# ============================================================

score_reference <- function(
  ref_idx
) {

  reference_counts <- counts[
    ,
    ref_idx,
    drop = FALSE
  ]

  library_size <- as.numeric(
    colSums(
      reference_counts
    )
  )

  if (
    any(
      !is.finite(
        library_size
      )
    ) ||
    any(
      library_size <= 0
    )
  ) {

    stop(
      "INVALID_LIBRARY_SIZE"
    )
  }

  core_raw <- counts[
    core_index,
    ref_idx,
    drop = FALSE
  ]

  l2c <- log2(
    sweep(
      core_raw,
      2,
      library_size,
      "/"
    ) *
      1e6 +
      1
  )

  rownames(
    l2c
  ) <- core_genes

  gene_mean <- rowMeans(
    l2c
  )

  gene_sd <- apply(
    l2c,
    1,
    safe_sd
  )

  z <- sweep(
    l2c,
    1,
    gene_mean,
    "-"
  )

  z <- sweep(
    z,
    1,
    gene_sd,
    "/"
  )

  CORE2_STRONG <- (
    -z[
      "CDKN1B",
    ] -
    z[
      "GSN",
    ]
  ) / 2

  BGN_TIMP1_COMPONENT <- (
    z[
      "BGN",
    ] +
    z[
      "TIMP1",
    ]
  ) / 2

  CORE4 <- (
    z[
      "BGN",
    ] -
    z[
      "CDKN1B",
    ] -
    z[
      "GSN",
    ] +
    z[
      "TIMP1",
    ]
  ) / 4

  # Algebraic identity:
  # CORE4 = (CORE2_STRONG + BGN_TIMP1_COMPONENT)/2
  identity_error <- max(
    abs(
      CORE4 -
        (
          CORE2_STRONG +
          BGN_TIMP1_COMPONENT
        ) /
        2
    )
  )

  if (
    identity_error > 1e-12
  ) {

    stop(
      "CORE4_ALGEBRA_IDENTITY_FAILURE"
    )
  }

  data.table(
    patient_id =
      pt_meta$patient_id[
        ref_idx
      ],

    timepoint =
      pt_meta$timepoint[
        ref_idx
      ],

    patient_timepoint =
      pt_meta$patient_timepoint[
        ref_idx
      ],

    library_size =
      library_size,

    BGN_z =
      as.numeric(
        z[
          "BGN",
        ]
      ),

    CDKN1B_z =
      as.numeric(
        z[
          "CDKN1B",
        ]
      ),

    GSN_z =
      as.numeric(
        z[
          "GSN",
        ]
      ),

    TIMP1_z =
      as.numeric(
        z[
          "TIMP1",
        ]
      ),

    CORE2_STRONG =
      as.numeric(
        CORE2_STRONG
      ),

    BGN_TIMP1_COMPONENT =
      as.numeric(
        BGN_TIMP1_COMPONENT
      ),

    CORE4 =
      as.numeric(
        CORE4
      )
  )
}


paired_delta_table <- function(
  score_dt,
  paired_ids
) {

  out <- list()

  for (
    patient in paired_ids
  ) {

    before <- score_dt[
      patient_id ==
        patient &
      timepoint ==
        "Before"
    ]

    after <- score_dt[
      patient_id ==
        patient &
      timepoint ==
        "After"
    ]

    if (
      nrow(before) != 1 ||
      nrow(after) != 1
    ) {

      stop(
        "PAIR_MAPPING_FAILURE: ",
        patient
      )
    }

    delta_cdkn1b <- (
      after$CDKN1B_z -
      before$CDKN1B_z
    )

    delta_gsn <- (
      after$GSN_z -
      before$GSN_z
    )

    contribution_cdkn1b <- (
      -delta_cdkn1b / 2
    )

    contribution_gsn <- (
      -delta_gsn / 2
    )

    out[[
      length(out) + 1
    ]] <- data.table(
      patient_id =
        patient,

      delta_CORE2 =
        after$CORE2_STRONG -
        before$CORE2_STRONG,

      delta_CORE4 =
        after$CORE4 -
        before$CORE4,

      delta_BGN_TIMP1_component =
        after$BGN_TIMP1_COMPONENT -
        before$BGN_TIMP1_COMPONENT,

      delta_CDKN1B_z =
        delta_cdkn1b,

      contribution_minus_CDKN1B_over2 =
        contribution_cdkn1b,

      delta_GSN_z =
        delta_gsn,

      contribution_minus_GSN_over2 =
        contribution_gsn
    )
  }

  x <- rbindlist(
    out
  )

  algebra_error <- max(
    abs(
      x$delta_CORE2 -
        (
          x$contribution_minus_CDKN1B_over2 +
          x$contribution_minus_GSN_over2
        )
    )
  )

  core4_error <- max(
    abs(
      x$delta_CORE4 -
        (
          x$delta_CORE2 +
          x$delta_BGN_TIMP1_component
        ) /
        2
    )
  )

  if (
    algebra_error > 1e-12 ||
    core4_error > 1e-12
  ) {

    stop(
      "PAIRED_COMPONENT_ALGEBRA_FAILURE"
    )
  }

  x
}

# ============================================================
# 7. Independent CORE4 audit across 22 scenarios
# ============================================================

cat("===== 6. INDEPENDENT CORE4 / COMPONENT SCENARIO AUDIT =====\n\n")

scores7a <- fread(
  SCORES_7A_FILE,
  encoding = "UTF-8"
)

deltas7a <- fread(
  DELTAS_7A_FILE,
  encoding = "UTF-8"
)

results7a <- fread(
  RESULTS_7A_FILE,
  encoding = "UTF-8"
)

scenario_summary_list <- list()
core4_audit_list <- list()
pairwise_recalc_list <- list()

max_core4_score_diff <- 0
max_core4_delta_diff <- 0
max_core2_score_diff <- 0
max_core2_delta_diff <- 0

for (
  i in seq_len(
    nrow(
      scenario_dt
    )
  )
) {

  s <- scenario_dt[i]

  sid <- as.character(
    s$scenario_id
  )

  hyp <- as.character(
    s$hypothetical_EASC_patient
  )

  ref_keep <- technical_pass

  ref_keep <- (
    ref_keep &
    pt_meta$patient_id != hyp
  )

  ref_idx <- which(
    ref_keep
  )

  if (
    length(
      ref_idx
    ) !=
      as.integer(
        s$technical_reference_n
      )
  ) {

    stop(
      "REFERENCE_N_MISMATCH: ",
      sid
    )
  }

  paired_ids <- baseline_pairs

  if (
    hyp %in%
      paired_ids
  ) {

    paired_ids <- paired_ids[
      paired_ids != hyp
    ]
  }

  if (
    length(
      paired_ids
    ) !=
      as.integer(
        s$paired_candidate_n
      )
  ) {

    stop(
      "PAIRED_N_MISMATCH: ",
      sid
    )
  }

  sc <- score_reference(
    ref_idx
  )

  pair_delta <- paired_delta_table(
    sc,
    paired_ids
  )

  pair_delta[
    ,
    `:=`(
      scenario_id =
        sid,
      hypothetical_EASC_patient =
        hyp,
      n_paired =
        length(
          paired_ids
        )
    )
  ]

  pairwise_recalc_list[[
    length(
      pairwise_recalc_list
    ) + 1
  ]] <- pair_delta

  # ----------------------------------------------------------
  # Compare independent scores with 7A.
  # ----------------------------------------------------------

  old_scores <- scores7a[
    scenario_id ==
      sid,
    .(
      patient_timepoint,
      CORE2_STRONG,
      CORE4
    )
  ]

  cmp_scores <- merge(
    sc[
      ,
      .(
        patient_timepoint,
        CORE2_recalc =
          CORE2_STRONG,
        CORE4_recalc =
          CORE4
      )
    ],
    old_scores,
    by = "patient_timepoint"
  )

  if (
    nrow(
      cmp_scores
    ) !=
      nrow(sc)
  ) {

    stop(
      "7A_SCORE_MEMBERSHIP_MISMATCH: ",
      sid
    )
  }

  core2_score_diff <- max(
    abs(
      cmp_scores$CORE2_recalc -
        cmp_scores$CORE2_STRONG
    )
  )

  core4_score_diff <- max(
    abs(
      cmp_scores$CORE4_recalc -
        cmp_scores$CORE4
    )
  )

  old_delta <- deltas7a[
    scenario_id ==
      sid,
    .(
      patient_id,
      delta_CORE2,
      delta_CORE4
    )
  ]

  cmp_delta <- merge(
    pair_delta[
      ,
      .(
        patient_id,
        CORE2_delta_recalc =
          delta_CORE2,
        CORE4_delta_recalc =
          delta_CORE4
      )
    ],
    old_delta,
    by = "patient_id"
  )

  if (
    nrow(
      cmp_delta
    ) !=
      length(
        paired_ids
      )
  ) {

    stop(
      "7A_DELTA_MEMBERSHIP_MISMATCH: ",
      sid
    )
  }

  core2_delta_diff <- max(
    abs(
      cmp_delta$CORE2_delta_recalc -
        cmp_delta$delta_CORE2
    )
  )

  core4_delta_diff <- max(
    abs(
      cmp_delta$CORE4_delta_recalc -
        cmp_delta$delta_CORE4
    )
  )

  max_core2_score_diff <- max(
    max_core2_score_diff,
    core2_score_diff
  )

  max_core4_score_diff <- max(
    max_core4_score_diff,
    core4_score_diff
  )

  max_core2_delta_diff <- max(
    max_core2_delta_diff,
    core2_delta_diff
  )

  max_core4_delta_diff <- max(
    max_core4_delta_diff,
    core4_delta_diff
  )

  core4_audit_list[[
    length(
      core4_audit_list
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

    n_paired =
      length(
        paired_ids
      ),

    CORE2_score_max_abs_diff =
      core2_score_diff,

    CORE4_score_max_abs_diff =
      core4_score_diff,

    CORE2_delta_max_abs_diff =
      core2_delta_diff,

    CORE4_delta_max_abs_diff =
      core4_delta_diff
  )

  # ----------------------------------------------------------
  # Scenario summary.
  # ----------------------------------------------------------

  scenario_summary_list[[
    length(
      scenario_summary_list
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

    n_paired =
      length(
        paired_ids
      ),

    paired_ids =
      paste(
        paired_ids,
        collapse = ";"
      ),

    CORE2_mean_delta =
      mean(
        pair_delta$delta_CORE2
      ),

    CORE2_median_delta =
      median(
        pair_delta$delta_CORE2
      ),

    CORE2_negative_n =
      sum(
        pair_delta$delta_CORE2 < 0
      ),

    CORE2_positive_n =
      sum(
        pair_delta$delta_CORE2 > 0
      ),

    CORE4_mean_delta =
      mean(
        pair_delta$delta_CORE4
      ),

    CORE4_median_delta =
      median(
        pair_delta$delta_CORE4
      ),

    CORE4_negative_n =
      sum(
        pair_delta$delta_CORE4 < 0
      ),

    CORE4_positive_n =
      sum(
        pair_delta$delta_CORE4 > 0
      ),

    BGN_TIMP1_mean_delta =
      mean(
        pair_delta$
          delta_BGN_TIMP1_component
      ),

    BGN_TIMP1_median_delta =
      median(
        pair_delta$
          delta_BGN_TIMP1_component
      ),

    BGN_TIMP1_negative_n =
      sum(
        pair_delta$
          delta_BGN_TIMP1_component < 0
      ),

    BGN_TIMP1_positive_n =
      sum(
        pair_delta$
          delta_BGN_TIMP1_component > 0
      ),

    CDKN1B_mean_CORE2_contribution =
      mean(
        pair_delta$
          contribution_minus_CDKN1B_over2
      ),

    GSN_mean_CORE2_contribution =
      mean(
        pair_delta$
          contribution_minus_GSN_over2
      ),

    CORE2_mean_direction =
      direction_label(
        mean(
          pair_delta$delta_CORE2
        )
      ),

    CORE4_mean_direction =
      direction_label(
        mean(
          pair_delta$delta_CORE4
        )
      ),

    BGN_TIMP1_mean_direction =
      direction_label(
        mean(
          pair_delta$
            delta_BGN_TIMP1_component
        )
      )
  )

  cat(
    sprintf(
      "[%02d/22] EASC=%-4s pairs=%d  CORE2=% .4f  CORE4=% .4f  BGN/TIMP1=% .4f\n",
      i,
      hyp,
      length(
        paired_ids
      ),
      mean(
        pair_delta$delta_CORE2
      ),
      mean(
        pair_delta$delta_CORE4
      ),
      mean(
        pair_delta$
          delta_BGN_TIMP1_component
      )
    )
  )
}

cat("\n")

core4_audit <- rbindlist(
  core4_audit_list
)

scenario_summary <- rbindlist(
  scenario_summary_list
)

pairwise_recalc <- rbindlist(
  pairwise_recalc_list
)

fwrite(
  core4_audit,
  CORE4_AUDIT_OUT,
  bom = TRUE
)

fwrite(
  scenario_summary,
  SCENARIO_OUT,
  bom = TRUE
)

cat(
  "Max CORE2 score difference vs 7A :",
  format(
    max_core2_score_diff,
    scientific = TRUE
  ),
  "\n"
)

cat(
  "Max CORE4 score difference vs 7A :",
  format(
    max_core4_score_diff,
    scientific = TRUE
  ),
  "\n"
)

cat(
  "Max CORE2 delta difference vs 7A :",
  format(
    max_core2_delta_diff,
    scientific = TRUE
  ),
  "\n"
)

cat(
  "Max CORE4 delta difference vs 7A :",
  format(
    max_core4_delta_diff,
    scientific = TRUE
  ),
  "\n\n"
)

if (
  max_core2_score_diff > tol ||
  max_core4_score_diff > tol ||
  max_core2_delta_diff > tol ||
  max_core4_delta_diff > tol
) {

  stop(
    "INDEPENDENT_CORE4_REPLICATION_FAILED"
  )
}

cat(
  "✓ CORE4 independently reproduced within tolerance.\n\n"
)

# ============================================================
# 8. CDKN1B / GSN component summary
# ============================================================

cat("===== 7. CDKN1B / GSN CONTRIBUTION AUDIT =====\n\n")

component_hyp <- component7b[
  grepl(
    "^S_HYPOTHETICAL_EASC_",
    scenario_id
  )
]

component_summary <- component_hyp[
  ,
  .(
    n_paired =
      .N,

    mean_CDKN1B_contribution =
      mean(
        contribution_minus_CDKN1B_over2
      ),

    median_CDKN1B_contribution =
      median(
        contribution_minus_CDKN1B_over2
      ),

    CDKN1B_negative_contribution_n =
      sum(
        contribution_minus_CDKN1B_over2 < 0
      ),

    CDKN1B_positive_contribution_n =
      sum(
        contribution_minus_CDKN1B_over2 > 0
      ),

    mean_GSN_contribution =
      mean(
        contribution_minus_GSN_over2
      ),

    median_GSN_contribution =
      median(
        contribution_minus_GSN_over2
      ),

    GSN_negative_contribution_n =
      sum(
        contribution_minus_GSN_over2 < 0
      ),

    GSN_positive_contribution_n =
      sum(
        contribution_minus_GSN_over2 > 0
      )
  ),
  by = .(
    scenario_id,
    hypothetical_EASC_patient
  )
]

component_summary[
  ,
  CDKN1B_mean_direction :=
    vapply(
      mean_CDKN1B_contribution,
      direction_label,
      character(1)
    )
]

component_summary[
  ,
  GSN_mean_direction :=
    vapply(
      mean_GSN_contribution,
      direction_label,
      character(1)
    )
]

component_check <- merge(
  component_summary,
  scenario_summary[
    ,
    .(
      scenario_id,
      CORE2_mean_delta
    )
  ],
  by = "scenario_id"
)

component_check[
  ,
  contribution_sum_error :=
    abs(
      CORE2_mean_delta -
        (
          mean_CDKN1B_contribution +
          mean_GSN_contribution
        )
    )
]

if (
  max(
    component_check$
      contribution_sum_error
  ) > tol
) {

  stop(
    "CORE2_COMPONENT_SUMMARY_ALGEBRA_FAILURE"
  )
}

fwrite(
  component_check,
  COMPONENT_OUT,
  bom = TRUE
)

cdkn1b_mean_negative_n <- sum(
  component_check$
    mean_CDKN1B_contribution < 0
)

gsn_mean_negative_n <- sum(
  component_check$
    mean_GSN_contribution < 0
)

cdkn1b_mean_positive_n <- sum(
  component_check$
    mean_CDKN1B_contribution > 0
)

gsn_mean_positive_n <- sum(
  component_check$
    mean_GSN_contribution > 0
)

both_component_negative_n <- sum(
  component_check$
    mean_CDKN1B_contribution < 0 &
  component_check$
    mean_GSN_contribution < 0
)

cat(
  "CDKN1B mean contribution negative:",
  cdkn1b_mean_negative_n,
  "/22\n"
)

cat(
  "GSN mean contribution negative    :",
  gsn_mean_negative_n,
  "/22\n"
)

cat(
  "Both mean contributions negative  :",
  both_component_negative_n,
  "/22\n\n"
)

# ============================================================
# 9. Scenario-level comparator envelope
# ============================================================

cat("===== 8. SCENARIO-LEVEL COMPARATOR ENVELOPE =====\n\n")

core2_mean_negative_n <- sum(
  scenario_summary$
    CORE2_mean_delta < 0
)

core2_median_negative_n <- sum(
  scenario_summary$
    CORE2_median_delta < 0
)

core4_mean_negative_n <- sum(
  scenario_summary$
    CORE4_mean_delta < 0
)

core4_mean_positive_n <- sum(
  scenario_summary$
    CORE4_mean_delta > 0
)

core4_median_negative_n <- sum(
  scenario_summary$
    CORE4_median_delta < 0
)

core4_median_positive_n <- sum(
  scenario_summary$
    CORE4_median_delta > 0
)

weak_mean_negative_n <- sum(
  scenario_summary$
    BGN_TIMP1_mean_delta < 0
)

weak_mean_positive_n <- sum(
  scenario_summary$
    BGN_TIMP1_mean_delta > 0
)

weak_median_negative_n <- sum(
  scenario_summary$
    BGN_TIMP1_median_delta < 0
)

weak_median_positive_n <- sum(
  scenario_summary$
    BGN_TIMP1_median_delta > 0
)

cat(
  "CORE2 mean negative        :",
  core2_mean_negative_n,
  "/22\n"
)

cat(
  "CORE4 mean negative        :",
  core4_mean_negative_n,
  "/22\n"
)

cat(
  "CORE4 mean positive        :",
  core4_mean_positive_n,
  "/22\n"
)

cat(
  "BGN/TIMP1 mean negative    :",
  weak_mean_negative_n,
  "/22\n"
)

cat(
  "BGN/TIMP1 mean positive    :",
  weak_mean_positive_n,
  "/22\n\n"
)

# ============================================================
# 10. Full paired LOPO
#
# Important:
# each left-out patient is removed from:
#   (a) z-score reference
#   (b) paired effect set
# ============================================================

cat("===== 9. FULL PAIRED LOPO RECOMPUTATION =====\n\n")

lopo_list <- list()

run_counter <- 0L

for (
  i in seq_len(
    nrow(
      scenario_dt
    )
  )
) {

  s <- scenario_dt[i]

  sid <- as.character(
    s$scenario_id
  )

  hyp <- as.character(
    s$hypothetical_EASC_patient
  )

  paired_ids <- baseline_pairs

  if (
    hyp %in%
      paired_ids
  ) {

    paired_ids <- paired_ids[
      paired_ids != hyp
    ]
  }

  base_ref_keep <- (
    technical_pass &
    pt_meta$patient_id != hyp
  )

  for (
    leave_out in paired_ids
  ) {

    run_counter <- (
      run_counter +
        1L
    )

    lopo_ref_keep <- (
      base_ref_keep &
      pt_meta$patient_id !=
        leave_out
    )

    ref_idx <- which(
      lopo_ref_keep
    )

    remaining_pairs <- paired_ids[
      paired_ids !=
        leave_out
    ]

    sc <- score_reference(
      ref_idx
    )

    d <- paired_delta_table(
      sc,
      remaining_pairs
    )

    lopo_list[[
      length(
        lopo_list
      ) + 1
    ]] <- data.table(
      scenario_id =
        sid,

      hypothetical_EASC_patient =
        hyp,

      left_out_patient =
        leave_out,

      original_pair_n =
        length(
          paired_ids
        ),

      remaining_pair_n =
        length(
          remaining_pairs
        ),

      reference_n_after_LOPO =
        length(
          ref_idx
        ),

      remaining_pair_ids =
        paste(
          remaining_pairs,
          collapse = ";"
        ),

      CORE2_mean_delta =
        mean(
          d$delta_CORE2
        ),

      CORE2_median_delta =
        median(
          d$delta_CORE2
        ),

      CORE2_negative_n =
        sum(
          d$delta_CORE2 < 0
        ),

      CORE2_positive_n =
        sum(
          d$delta_CORE2 > 0
        ),

      CORE2_all_remaining_negative =
        all(
          d$delta_CORE2 < 0
        ),

      CORE4_mean_delta =
        mean(
          d$delta_CORE4
        ),

      CORE4_median_delta =
        median(
          d$delta_CORE4
        ),

      CORE4_negative_n =
        sum(
          d$delta_CORE4 < 0
        ),

      CORE4_positive_n =
        sum(
          d$delta_CORE4 > 0
        ),

      BGN_TIMP1_mean_delta =
        mean(
          d$delta_BGN_TIMP1_component
        ),

      BGN_TIMP1_median_delta =
        median(
          d$delta_BGN_TIMP1_component
        )
    )
  }
}

lopo <- rbindlist(
  lopo_list
)

fwrite(
  lopo,
  LOPO_OUT,
  bom = TRUE
)

expected_lopo_runs <- (
  sum(
    scenario_dt$
      paired_candidate_n
  )
)

cat(
  "LOPO runs:",
  nrow(
    lopo
  ),
  "\n"
)

cat(
  "Expected LOPO runs:",
  expected_lopo_runs,
  "\n\n"
)

if (
  nrow(
    lopo
  ) !=
    expected_lopo_runs
) {

  stop(
    "LOPO_RUN_COUNT_MISMATCH"
  )
}

# 17 scenarios n=5 -> 85 LOPO runs with remaining n=4
# 5 scenarios n=4 -> 20 LOPO runs with remaining n=3
lopo_remaining_n3 <- sum(
  lopo$
    remaining_pair_n == 3
)

lopo_remaining_n4 <- sum(
  lopo$
    remaining_pair_n == 4
)

core2_lopo_mean_negative_n <- sum(
  lopo$
    CORE2_mean_delta < 0
)

core2_lopo_median_negative_n <- sum(
  lopo$
    CORE2_median_delta < 0
)

core2_lopo_majority_negative_n <- sum(
  lopo$
    CORE2_negative_n >
    lopo$
      remaining_pair_n /
      2
)

core2_lopo_all_negative_n <- sum(
  lopo$
    CORE2_all_remaining_negative
)

core4_lopo_mean_negative_n <- sum(
  lopo$
    CORE4_mean_delta < 0
)

core4_lopo_mean_positive_n <- sum(
  lopo$
    CORE4_mean_delta > 0
)

weak_lopo_mean_negative_n <- sum(
  lopo$
    BGN_TIMP1_mean_delta < 0
)

weak_lopo_mean_positive_n <- sum(
  lopo$
    BGN_TIMP1_mean_delta > 0
)

cat(
  "LOPO remaining n=3:",
  lopo_remaining_n3,
  "\n"
)

cat(
  "LOPO remaining n=4:",
  lopo_remaining_n4,
  "\n\n"
)

cat(
  "CORE2 LOPO mean negative       :",
  core2_lopo_mean_negative_n,
  "/",
  nrow(
    lopo
  ),
  "\n"
)

cat(
  "CORE2 LOPO median negative     :",
  core2_lopo_median_negative_n,
  "/",
  nrow(
    lopo
  ),
  "\n"
)

cat(
  "CORE2 LOPO majority negative   :",
  core2_lopo_majority_negative_n,
  "/",
  nrow(
    lopo
  ),
  "\n"
)

cat(
  "CORE2 LOPO all patients neg.   :",
  core2_lopo_all_negative_n,
  "/",
  nrow(
    lopo
  ),
  "\n\n"
)

# ============================================================
# 11. Deterministic interpretation pattern
# ============================================================

cat("===== 10. DETERMINISTIC PATTERN CLASSIFICATION =====\n\n")

if (
  core2_lopo_mean_negative_n ==
    nrow(lopo) &&
  core2_lopo_median_negative_n ==
    nrow(lopo) &&
  core2_lopo_majority_negative_n ==
    nrow(lopo)
) {

  core2_lopo_status <- (
    "CORE2_NEGATIVE_ROBUST_TO_ALL_PAIRED_LOPO_RUNS"
  )

} else if (
  core2_lopo_mean_negative_n >=
    0.90 *
    nrow(lopo)
) {

  core2_lopo_status <- (
    "CORE2_NEGATIVE_MOSTLY_ROBUST_TO_PAIRED_LOPO"
  )

} else {

  core2_lopo_status <- (
    "CORE2_NEGATIVE_SENSITIVE_TO_PAIRED_PATIENT_INFLUENCE"
  )
}


if (
  both_component_negative_n == 22
) {

  component_pattern <- (
    "CDKN1B_AND_GSN_BOTH_CONTRIBUTE_NEGATIVELY_IN_ALL_22_SCENARIOS"
  )

} else if (
  cdkn1b_mean_negative_n == 22 &&
  gsn_mean_negative_n < 22
) {

  component_pattern <- (
    "CDKN1B_CONTRIBUTION_CONSISTENTLY_NEGATIVE_GSN_MIXED"
  )

} else if (
  gsn_mean_negative_n == 22 &&
  cdkn1b_mean_negative_n < 22
) {

  component_pattern <- (
    "GSN_CONTRIBUTION_CONSISTENTLY_NEGATIVE_CDKN1B_MIXED"
  )

} else {

  component_pattern <- (
    "CDKN1B_GSN_CONTRIBUTIONS_MIXED_ACROSS_SCENARIOS"
  )
}


if (
  weak_mean_negative_n == 22 &&
  core4_mean_negative_n == 22
) {

  module_pattern <- (
    "BOTH_CORE_ARMS_AND_CORE4_NEGATIVE_ACROSS_ALL_22_SCENARIOS"
  )

} else if (
  weak_mean_positive_n == 22
) {

  module_pattern <- (
    "CORE2_STRONG_NEGATIVE_BUT_BGN_TIMP1_ARM_POSITIVE_ACROSS_ALL_22"
  )

} else if (
  core4_mean_negative_n == 22
) {

  module_pattern <- (
    "CORE4_NEGATIVE_ALL_22_WITH_MIXED_BGN_TIMP1_ARM"
  )

} else if (
  core4_mean_positive_n == 22
) {

  module_pattern <- (
    "CORE4_POSITIVE_ALL_22_DESPITE_NEGATIVE_CORE2_STRONG"
  )

} else {

  module_pattern <- (
    "CORE4_AND_BGN_TIMP1_PATTERN_MIXED_ACROSS_HISTOLOGY_SCENARIOS"
  )
}


cat(
  "CORE2 LOPO status:\n",
  core2_lopo_status,
  "\n\n",
  sep = ""
)

cat(
  "CORE2 component pattern:\n",
  component_pattern,
  "\n\n",
  sep = ""
)

cat(
  "Four-gene module pattern:\n",
  module_pattern,
  "\n\n",
  sep = ""
)

# ============================================================
# 12. Global robustness table
# ============================================================

robustness <- data.table(
  metric = c(
    "hypothetical_EASC_scenarios",
    "CORE2_scenario_mean_negative",
    "CORE2_scenario_median_negative",
    "CORE4_scenario_mean_negative",
    "CORE4_scenario_mean_positive",
    "CORE4_scenario_median_negative",
    "CORE4_scenario_median_positive",
    "BGN_TIMP1_scenario_mean_negative",
    "BGN_TIMP1_scenario_mean_positive",
    "CDKN1B_mean_contribution_negative",
    "CDKN1B_mean_contribution_positive",
    "GSN_mean_contribution_negative",
    "GSN_mean_contribution_positive",
    "both_CDKN1B_GSN_mean_contributions_negative",
    "LOPO_total_runs",
    "LOPO_remaining_n3",
    "LOPO_remaining_n4",
    "CORE2_LOPO_mean_negative",
    "CORE2_LOPO_median_negative",
    "CORE2_LOPO_majority_negative",
    "CORE2_LOPO_all_remaining_patients_negative",
    "CORE4_LOPO_mean_negative",
    "CORE4_LOPO_mean_positive",
    "BGN_TIMP1_LOPO_mean_negative",
    "BGN_TIMP1_LOPO_mean_positive",
    "CORE2_LOPO_status",
    "CORE2_component_pattern",
    "four_gene_module_pattern"
  ),

  value = as.character(
    c(
      22,
      core2_mean_negative_n,
      core2_median_negative_n,
      core4_mean_negative_n,
      core4_mean_positive_n,
      core4_median_negative_n,
      core4_median_positive_n,
      weak_mean_negative_n,
      weak_mean_positive_n,
      cdkn1b_mean_negative_n,
      cdkn1b_mean_positive_n,
      gsn_mean_negative_n,
      gsn_mean_positive_n,
      both_component_negative_n,
      nrow(
        lopo
      ),
      lopo_remaining_n3,
      lopo_remaining_n4,
      core2_lopo_mean_negative_n,
      core2_lopo_median_negative_n,
      core2_lopo_majority_negative_n,
      core2_lopo_all_negative_n,
      core4_lopo_mean_negative_n,
      core4_lopo_mean_positive_n,
      weak_lopo_mean_negative_n,
      weak_lopo_mean_positive_n,
      core2_lopo_status,
      component_pattern,
      module_pattern
    )
  )
)

fwrite(
  robustness,
  ROBUSTNESS_OUT,
  bom = TRUE
)

# ============================================================
# 13. Frozen Step19-21 integrity
# ============================================================

cat("===== 11. FROZEN STEP19-21 INTEGRITY =====\n\n")

frozen <- fread(
  FROZEN_INV,
  encoding = "UTF-8"
)

changed <- character()

for (
  i in seq_len(
    nrow(
      frozen
    )
  )
) {

  raw_path <- as.character(
    frozen$path[i]
  )

  path <- resolve_project_path(
    raw_path
  )

  if (
    !file.exists(
      path
    )
  ) {

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
    is.finite(
      expected
    ) &&
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
    "STEP22A_FROZEN_PROJECT_CHANGED"
  )
}

# ============================================================
# 14. Summary
# ============================================================

summary_lines <- c(
  "STEP22A-7C CORE4 + COMPONENT + PAIRED LOPO AUDIT",
  "",
  "GOVERNANCE",
  "Primary ESCC-only analysis: BLOCKED",
  "Secondary histology sensitivity: AUTHORIZED",
  "Scenario selection: FORBIDDEN",
  "",
  "INDEPENDENT COMPARATOR REPLICATION",
  paste0(
    "CORE2 score max abs diff vs 7A: ",
    format(
      max_core2_score_diff,
      scientific = TRUE
    )
  ),
  paste0(
    "CORE4 score max abs diff vs 7A: ",
    format(
      max_core4_score_diff,
      scientific = TRUE
    )
  ),
  paste0(
    "CORE2 delta max abs diff vs 7A: ",
    format(
      max_core2_delta_diff,
      scientific = TRUE
    )
  ),
  paste0(
    "CORE4 delta max abs diff vs 7A: ",
    format(
      max_core4_delta_diff,
      scientific = TRUE
    )
  ),
  "",
  "22-SCENARIO DIRECTIONS",
  paste0(
    "CORE2 mean negative: ",
    core2_mean_negative_n,
    "/22"
  ),
  paste0(
    "CORE4 mean negative: ",
    core4_mean_negative_n,
    "/22"
  ),
  paste0(
    "CORE4 mean positive: ",
    core4_mean_positive_n,
    "/22"
  ),
  paste0(
    "BGN/TIMP1 component mean negative: ",
    weak_mean_negative_n,
    "/22"
  ),
  paste0(
    "BGN/TIMP1 component mean positive: ",
    weak_mean_positive_n,
    "/22"
  ),
  "",
  "CORE2 COMPONENTS",
  paste0(
    "CDKN1B mean contribution negative: ",
    cdkn1b_mean_negative_n,
    "/22"
  ),
  paste0(
    "GSN mean contribution negative: ",
    gsn_mean_negative_n,
    "/22"
  ),
  paste0(
    "Both mean contributions negative: ",
    both_component_negative_n,
    "/22"
  ),
  "",
  "PAIRED LOPO",
  paste0(
    "LOPO runs: ",
    nrow(
      lopo
    )
  ),
  paste0(
    "Remaining n=3 runs: ",
    lopo_remaining_n3
  ),
  paste0(
    "Remaining n=4 runs: ",
    lopo_remaining_n4
  ),
  paste0(
    "CORE2 LOPO mean negative: ",
    core2_lopo_mean_negative_n,
    "/",
    nrow(
      lopo
    )
  ),
  paste0(
    "CORE2 LOPO median negative: ",
    core2_lopo_median_negative_n,
    "/",
    nrow(
      lopo
    )
  ),
  paste0(
    "CORE2 LOPO majority negative: ",
    core2_lopo_majority_negative_n,
    "/",
    nrow(
      lopo
    )
  ),
  paste0(
    "CORE2 LOPO all remaining patients negative: ",
    core2_lopo_all_negative_n,
    "/",
    nrow(
      lopo
    )
  ),
  "",
  paste0(
    "CORE2 LOPO status: ",
    core2_lopo_status
  ),
  paste0(
    "CORE2 component pattern: ",
    component_pattern
  ),
  paste0(
    "Four-gene module pattern: ",
    module_pattern
  ),
  "",
  "INTERPRETATION BOUNDARY",
  "CORE4 is a pre-specified comparator, not a rescue endpoint.",
  "LOPO tests influence, not significance.",
  "No scenario was selected.",
  "Primary ESCC-only validation remains unavailable.",
  "Frozen Step19-21 modified: NO",
  "Manuscript modified: NO"
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
cat("STEP22A-7C COMPLETE\n")
cat("====================================================================\n\n")

cat(
  "PRIMARY ESCC-ONLY             : BLOCKED\n"
)

cat(
  "SECONDARY COMPARATOR AUDIT    : COMPLETE\n\n"
)

cat(
  "CORE4 INDEPENDENT AUDIT       : PASS\n"
)

cat(
  "CORE4 SCORE MAX ABS DIFF      :",
  format(
    max_core4_score_diff,
    scientific = TRUE
  ),
  "\n"
)

cat(
  "CORE4 DELTA MAX ABS DIFF      :",
  format(
    max_core4_delta_diff,
    scientific = TRUE
  ),
  "\n\n"
)

cat(
  "22-SCENARIO DIRECTION:\n"
)

cat(
  "  CORE2 MEAN NEGATIVE         :",
  core2_mean_negative_n,
  "/ 22\n"
)

cat(
  "  CORE4 MEAN NEGATIVE         :",
  core4_mean_negative_n,
  "/ 22\n"
)

cat(
  "  CORE4 MEAN POSITIVE         :",
  core4_mean_positive_n,
  "/ 22\n"
)

cat(
  "  BGN/TIMP1 MEAN NEGATIVE     :",
  weak_mean_negative_n,
  "/ 22\n"
)

cat(
  "  BGN/TIMP1 MEAN POSITIVE     :",
  weak_mean_positive_n,
  "/ 22\n\n"
)

cat(
  "CORE2 ARITHMETIC COMPONENTS:\n"
)

cat(
  "  CDKN1B CONTRIB. NEGATIVE    :",
  cdkn1b_mean_negative_n,
  "/ 22\n"
)

cat(
  "  GSN CONTRIB. NEGATIVE       :",
  gsn_mean_negative_n,
  "/ 22\n"
)

cat(
  "  BOTH NEGATIVE               :",
  both_component_negative_n,
  "/ 22\n\n"
)

cat(
  "FULL PAIRED LOPO:\n"
)

cat(
  "  TOTAL RUNS                  :",
  nrow(
    lopo
  ),
  "\n"
)

cat(
  "  REMAINING n=3               :",
  lopo_remaining_n3,
  "\n"
)

cat(
  "  REMAINING n=4               :",
  lopo_remaining_n4,
  "\n"
)

cat(
  "  CORE2 MEAN NEGATIVE         :",
  core2_lopo_mean_negative_n,
  "/",
  nrow(
    lopo
  ),
  "\n"
)

cat(
  "  CORE2 MEDIAN NEGATIVE       :",
  core2_lopo_median_negative_n,
  "/",
  nrow(
    lopo
  ),
  "\n"
)

cat(
  "  MAJORITY NEGATIVE           :",
  core2_lopo_majority_negative_n,
  "/",
  nrow(
    lopo
  ),
  "\n"
)

cat(
  "  ALL REMAINING PAIRS NEG.    :",
  core2_lopo_all_negative_n,
  "/",
  nrow(
    lopo
  ),
  "\n\n"
)

cat(
  "CORE2 LOPO STATUS:\n",
  core2_lopo_status,
  "\n\n",
  sep = ""
)

cat(
  "CORE2 COMPONENT PATTERN:\n",
  component_pattern,
  "\n\n",
  sep = ""
)

cat(
  "FOUR-GENE MODULE PATTERN:\n",
  module_pattern,
  "\n\n",
  sep = ""
)

cat(
  "IMPORTANT:\n"
)

cat(
  "  CORE4 is a comparator, not a rescue endpoint.\n"
)

cat(
  "  LOPO is an influence audit, not a significance test.\n"
)

cat(
  "  Primary ESCC-only validation remains BLOCKED.\n\n"
)

cat(
  "FROZEN STEP19-21 MODIFIED     : NO\n"
)

cat(
  "MANUSCRIPT MODIFIED           : NO\n\n"
)

cat(
  "STOP HERE.\n"
)

cat(
  "REVIEW MODULE PATTERN BEFORE ANY FURTHER ANALYSIS.\n\n"
)
