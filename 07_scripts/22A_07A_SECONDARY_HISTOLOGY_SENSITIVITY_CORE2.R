# ============================================================
# ESCC Neoadjuvant Project
# STEP22A-7A
#
# SECONDARY HISTOLOGY SENSITIVITY CORE2 ANALYSIS
#
# IMPORTANT GOVERNANCE STATUS
#
# PRIMARY ESCC-ONLY ANALYSIS:
#   REMAINS BLOCKED
#
# THIS STEP:
#   - opens ONLY the pre-frozen secondary histology branch
#   - evaluates ALL 22 hypothetical EASC identities
#   - also calculates one unresolved technical baseline
#     for descriptive reference only
#
# NO scenario may be promoted to primary.
# NO scenario selection based on results.
#
# Frozen scoring:
#
# raw pseudobulk counts
# -> library-size CPM
# -> log2(CPM + 1)
# -> per-gene z-score across CURRENT scenario-specific
#    eligible tumor reference
#
# CORE4 =
#   (BGN_z - CDKN1B_z - GSN_z + TIMP1_z) / 4
#
# CORE2_STRONG =
#   (-CDKN1B_z - GSN_z) / 2
#
# CORE2_WEAK =
#   (BGN_z + TIMP1_z) / 2
#
# Paired secondary endpoint:
#   delta_CORE2 =
#     After - Before
#
# Positive = frozen replication direction.
#
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
cat("STEP22A-7A: SECONDARY HISTOLOGY SENSITIVITY CORE2\n")
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

RESULT_DIR <- file.path(
  PROJECT,
  "04_results",
  "external_validation",
  "Step22A"
)

TABLE_DIR <- file.path(
  PROJECT,
  "06_tables",
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

dir.create(
  TABLE_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  LOG_DIR,
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

GOVERNANCE_AUDIT <- file.path(
  META_DIR,
  "STEP22A_06F_FINAL_GOVERNANCE_AUDIT.csv"
)

PRIMARY_GATE <- file.path(
  META_DIR,
  "STEP22A_ANALYSIS_GATE_STATUS.csv"
)

HIST_TECH_FILE <- file.path(
  META_DIR,
  "STEP22A_06E_HISTOLOGY_SENSITIVITY_TECHNICAL_INTERSECTION.csv"
)

PAIRED_LOCK_FILE <- file.path(
  META_DIR,
  "STEP22A_06E_FINAL_PAIRED_TECHNICAL_LOCK.csv"
)

FROZEN_INV <- file.path(
  META_DIR,
  "STEP22A_PREANALYSIS_FROZEN_FILE_INVENTORY.csv"
)

AUTH_OUT <- file.path(
  META_DIR,
  "STEP22A_07A_SECONDARY_SENSITIVITY_AUTHORIZATION.csv"
)

SCORES_OUT <- file.path(
  TABLE_DIR,
  "STEP22A_07A_SECONDARY_SCENARIO_CORE_SCORES.csv"
)

DELTAS_OUT <- file.path(
  TABLE_DIR,
  "STEP22A_07A_SECONDARY_PAIRED_DELTAS.csv"
)

RESULTS_OUT <- file.path(
  TABLE_DIR,
  "STEP22A_07A_SECONDARY_PAIRED_SCENARIO_RESULTS.csv"
)

ROBUSTNESS_OUT <- file.path(
  TABLE_DIR,
  "STEP22A_07A_SECONDARY_HISTOLOGY_ROBUSTNESS_SUMMARY.csv"
)

SUMMARY_OUT <- file.path(
  RESULT_DIR,
  "STEP22A_07A_SECONDARY_HISTOLOGY_SENSITIVITY_SUMMARY.txt"
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


exact_signflip <- function(delta) {

  delta <- as.numeric(delta)

  n <- length(delta)

  if (n < 1) {

    return(
      list(
        observed = NA_real_,
        p_two_sided = NA_real_,
        p_directional = NA_real_,
        permutations = 0L
      )
    )
  }

  sign_grid <- expand.grid(
    rep(
      list(
        c(
          -1,
          1
        )
      ),
      n
    )
  )

  signs <- as.matrix(
    sign_grid
  )

  perm_effect <- as.numeric(
    signs %*% delta
  ) / n

  observed <- mean(delta)

  tol <- 1e-12

  p2 <- mean(
    abs(
      perm_effect
    ) >=
      abs(
        observed
      ) -
      tol
  )

  p1 <- mean(
    perm_effect >=
      observed -
      tol
  )

  list(
    observed = observed,
    p_two_sided = p2,
    p_directional = p1,
    permutations = nrow(
      signs
    )
  )
}


wilcox_signed_rank <- function(delta) {

  warning_text <- character()

  fit <- withCallingHandlers(
    tryCatch(
      wilcox.test(
        delta,
        mu = 0,
        alternative = "two.sided",
        paired = FALSE,
        exact = TRUE,
        correct = FALSE
      ),
      error = function(e) e
    ),
    warning = function(w) {

      warning_text <<- c(
        warning_text,
        conditionMessage(w)
      )

      invokeRestart(
        "muffleWarning"
      )
    }
  )

  if (
    inherits(
      fit,
      "error"
    )
  ) {

    return(
      list(
        statistic = NA_real_,
        p = NA_real_,
        warning = conditionMessage(
          fit
        )
      )
    )
  }

  list(
    statistic = unname(
      fit$statistic
    ),
    p = fit$p.value,
    warning = paste(
      unique(
        warning_text
      ),
      collapse = " | "
    )
  )
}


frozen_integrity <- function() {

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

  changed
}

# ============================================================
# 1. Primary gate MUST remain closed
# ============================================================

cat("===== 1. PRIMARY GOVERNANCE GATE =====\n\n")

gate <- fread(
  PRIMARY_GATE,
  encoding = "UTF-8"
)

get_gate <- function(x) {

  y <- gate[
    gate == x,
    allowed
  ]

  if (length(y) != 1) {

    stop(
      "Cannot uniquely resolve gate: ",
      x
    )
  }

  as.character(y)
}


primary_core_gate <- get_gate(
  "paired_CORE2_test"
)

manuscript_gate <- get_gate(
  "manuscript_update"
)


cat(
  "Primary paired CORE2 gate :",
  primary_core_gate,
  "\n"
)

cat(
  "Manuscript update gate    :",
  manuscript_gate,
  "\n\n"
)


if (
  primary_core_gate == "YES"
) {

  stop(
    "PRIMARY_CORE2_GATE_MUST_REMAIN_CLOSED"
  )
}


if (
  manuscript_gate == "YES"
) {

  stop(
    "MANUSCRIPT_GATE_MUST_REMAIN_CLOSED"
  )
}


cat(
  "✓ Primary ESCC-only analysis remains BLOCKED.\n"
)

cat(
  "✓ This run is SECONDARY SENSITIVITY ONLY.\n\n"
)

# ============================================================
# 2. Step22A-6F governance audit
# ============================================================

cat("===== 2. VERIFY STEP22A-6F GOVERNANCE =====\n\n")

gov <- fread(
  GOVERNANCE_AUDIT,
  encoding = "UTF-8"
)


if (
  any(
    gov$status != "PASS"
  )
) {

  print(
    gov[
      status != "PASS"
    ]
  )

  stop(
    "STEP22A_06F_GOVERNANCE_NOT_ALL_PASS"
  )
}


cat(
  "Governance checks:",
  nrow(gov),
  "\n"
)

cat(
  "Failures:",
  0,
  "\n\n"
)

# ============================================================
# 3. Create explicit SECONDARY authorization
#
# This does NOT change the primary gate file.
# ============================================================

cat("===== 3. SECONDARY-BRANCH AUTHORIZATION =====\n\n")

authorization <- data.table(
  item = c(
    "PRIMARY_ESCC_ONLY",
    "SECONDARY_HISTOLOGY_SENSITIVITY",
    "PRIMARY_GATE_FILE_MODIFIED",
    "SCENARIO_SELECTION_ALLOWED",
    "ALL_PRE_FROZEN_EASC_SCENARIOS_REQUIRED",
    "MANUSCRIPT_UPDATE_ALLOWED",
    "AUTHOR_NONRESPONSE_INTERPRETED_AS_HISTOLOGY_EVIDENCE"
  ),
  status = c(
    "BLOCKED",
    "AUTHORIZED",
    "NO",
    "NO",
    "YES_22_OF_22",
    "NO",
    "NO"
  )
)

fwrite(
  authorization,
  AUTH_OUT,
  bom = TRUE
)

print(
  authorization
)

cat(
  "\n✓ Secondary branch authorized without opening primary.\n\n"
)

# ============================================================
# 4. Load all 23 frozen scenarios
# ============================================================

cat("===== 4. LOAD PRE-FROZEN HISTOLOGY SCENARIOS =====\n\n")

scenarios <- fread(
  HIST_TECH_FILE,
  encoding = "UTF-8"
)


if (
  nrow(
    scenarios
  ) != 23
) {

  stop(
    "Expected 23 frozen histology scenarios"
  )
}


hypothetical <- scenarios[
  grepl(
    "^S_HYPOTHETICAL_EASC_",
    scenario_id
  )
]


baseline_scenario <- scenarios[
  scenario_id ==
    "S00_UNRESOLVED"
]


if (
  nrow(
    hypothetical
  ) != 22 ||
  nrow(
    baseline_scenario
  ) != 1
) {

  stop(
    "FROZEN_HISTOLOGY_SCENARIO_STRUCTURE_INVALID"
  )
}


cat(
  "Unresolved descriptive baseline:",
  nrow(
    baseline_scenario
  ),
  "\n"
)

cat(
  "Hypothetical EASC scenarios:",
  nrow(
    hypothetical
  ),
  "\n\n"
)

# ============================================================
# 5. Load frozen paired technical candidate IDs
# ============================================================

cat("===== 5. LOAD FROZEN PAIRED TECHNICAL SET =====\n\n")

paired_lock <- fread(
  PAIRED_LOCK_FILE,
  encoding = "UTF-8"
)

baseline_pairs <- paired_lock[
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
    "FROZEN_PAIRED_TECHNICAL_SET_CHANGED"
  )
}


cat(
  "Baseline paired technical candidates:",
  paste(
    baseline_pairs,
    collapse = ", "
  ),
  "\n\n"
)

# ============================================================
# 6. First opening of raw pseudobulk bundle for secondary work
# ============================================================

cat("===== 6. OPEN RAW PSEUDOBULK BUNDLE =====\n\n")

if (
  !file.exists(
    RAW_BUNDLE
  )
) {

  stop(
    "Missing raw pseudobulk bundle"
  )
}


bundle <- readRDS(
  RAW_BUNDLE
)


required_bundle <- c(
  "raw_counts_all_patient_timepoints",
  "patient_timepoint_metadata",
  "feature_metadata"
)


if (
  !all(
    required_bundle %in%
      names(bundle)
  )
) {

  stop(
    "RAW_PSEUDOBULK_BUNDLE_STRUCTURE_INVALID"
  )
}


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
    nrow(features)
) {

  stop(
    "FEATURE_COUNT_MATRIX_DIMENSION_MISMATCH"
  )
}


if (
  ncol(counts) !=
    nrow(pt_meta)
) {

  stop(
    "PATIENT_TIMEPOINT_MATRIX_DIMENSION_MISMATCH"
  )
}


if (
  !identical(
    colnames(counts),
    pt_meta$patient_timepoint
  )
) {

  stop(
    "PSEUDOBULK_METADATA_ORDER_MISMATCH"
  )
}


cat(
  "Raw pseudobulk matrix:",
  nrow(counts),
  "genes x",
  ncol(counts),
  "patient-timepoints\n"
)

cat(
  "This is SECONDARY sensitivity analysis only.\n\n"
)

# ============================================================
# 7. Resolve exact four CORE genes
# ============================================================

cat("===== 7. RESOLVE FROZEN CORE GENES =====\n\n")

features[
  ,
  symbol_upper :=
    toupper(
      gene_symbol
    )
]


core_order <- c(
  "BGN",
  "CDKN1B",
  "GSN",
  "TIMP1"
)


core_index <- integer(
  length(
    core_order
  )
)


for (
  i in seq_along(
    core_order
  )
) {

  hits <- which(
    features$symbol_upper ==
      core_order[i]
  )

  if (
    length(hits) != 1
  ) {

    stop(
      "CORE_GENE_MAPPING_NOT_UNIQUE: ",
      core_order[i]
    )
  }

  core_index[i] <- hits
}


names(
  core_index
) <- core_order


cat(
  "Frozen CORE genes resolved:",
  paste(
    core_order,
    collapse = " / "
  ),
  "\n"
)

cat(
  "Gene weights refit:",
  "NO\n\n"
)

# ============================================================
# 8. Technical-pass reference
# ============================================================

cat("===== 8. BASELINE TECHNICAL REFERENCE =====\n\n")

technical_pass <- (
  pt_meta$
    frozen_technical_status ==
    "PASS_GE20"
)


if (
  sum(
    technical_pass
  ) != 26
) {

  stop(
    "Expected 26 technical-pass patient-timepoints"
  )
}


cat(
  "Technical-pass patient-timepoints:",
  sum(
    technical_pass
  ),
  "/ 27\n\n"
)

# ============================================================
# 9. Scenario scoring function
# ============================================================

score_scenario <- function(
  scenario_id,
  hypothetical_patient,
  expected_reference_n,
  expected_pair_n,
  scenario_role
) {

  # ----------------------------------------------------------
  # Reference membership
  # ----------------------------------------------------------

  ref_keep <- technical_pass

  if (
    hypothetical_patient != "UNKNOWN"
  ) {

    ref_keep <- (
      ref_keep &
      pt_meta$patient_id !=
        hypothetical_patient
    )
  }


  ref_idx <- which(
    ref_keep
  )


  if (
    length(
      ref_idx
    ) !=
      expected_reference_n
  ) {

    stop(
      "SCENARIO_REFERENCE_N_MISMATCH: ",
      scenario_id,
      " observed=",
      length(ref_idx),
      " expected=",
      expected_reference_n
    )
  }

  # ----------------------------------------------------------
  # Frozen transform:
  # raw count -> CPM -> log2(CPM+1)
  #
  # IMPORTANT:
  # library size uses ALL 16,790 raw-count genes.
  # ----------------------------------------------------------

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
      "INVALID_LIBRARY_SIZE: ",
      scenario_id
    )
  }


  core_raw <- counts[
    core_index,
    ref_idx,
    drop = FALSE
  ]


  core_cpm <- sweep(
    core_raw,
    2,
    library_size,
    "/"
  ) * 1e6


  core_l2c <- log2(
    core_cpm + 1
  )


  rownames(
    core_l2c
  ) <- core_order


  # ----------------------------------------------------------
  # Frozen Step19L z-score:
  # each gene across current scenario reference.
  # ----------------------------------------------------------

  gene_mean <- rowMeans(
    core_l2c
  )

  gene_sd <- apply(
    core_l2c,
    1,
    safe_sd
  )


  z <- sweep(
    core_l2c,
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


  CORE2_STRONG <- (
    -z[
      "CDKN1B",
    ] -
    z[
      "GSN",
    ]
  ) / 2


  CORE2_WEAK <- (
    z[
      "BGN",
    ] +
    z[
      "TIMP1",
    ]
  ) / 2


  score_dt <- data.table(
    scenario_id =
      scenario_id,

    scenario_role =
      scenario_role,

    hypothetical_EASC_patient =
      hypothetical_patient,

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

    reference_n =
      length(
        ref_idx
      ),

    library_size =
      library_size,

    BGN_log2CPM =
      as.numeric(
        core_l2c[
          "BGN",
        ]
      ),

    CDKN1B_log2CPM =
      as.numeric(
        core_l2c[
          "CDKN1B",
        ]
      ),

    GSN_log2CPM =
      as.numeric(
        core_l2c[
          "GSN",
        ]
      ),

    TIMP1_log2CPM =
      as.numeric(
        core_l2c[
          "TIMP1",
        ]
      ),

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

    CORE4 =
      as.numeric(
        CORE4
      ),

    CORE2_STRONG =
      as.numeric(
        CORE2_STRONG
      ),

    CORE2_WEAK =
      as.numeric(
        CORE2_WEAK
      )
  )


  # ----------------------------------------------------------
  # Paired set:
  # start from frozen technical five;
  # remove hypothetical EASC if it is one of those five.
  # ----------------------------------------------------------

  paired_ids <- baseline_pairs

  if (
    hypothetical_patient %in%
      paired_ids
  ) {

    paired_ids <- paired_ids[
      paired_ids !=
        hypothetical_patient
    ]
  }


  if (
    length(
      paired_ids
    ) !=
      expected_pair_n
  ) {

    stop(
      "SCENARIO_PAIRED_N_MISMATCH: ",
      scenario_id,
      " observed=",
      length(
        paired_ids
      ),
      " expected=",
      expected_pair_n
    )
  }


  delta_list <- list()


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
        "INCOMPLETE_PAIRED_SCORE: ",
        scenario_id,
        " / ",
        patient
      )
    }


    delta_list[[
      patient
    ]] <- data.table(
      scenario_id =
        scenario_id,

      scenario_role =
        scenario_role,

      hypothetical_EASC_patient =
        hypothetical_patient,

      patient_id =
        patient,

      n_paired =
        length(
          paired_ids
        ),

      CORE2_before =
        before$CORE2_STRONG,

      CORE2_after =
        after$CORE2_STRONG,

      delta_CORE2 =
        after$CORE2_STRONG -
        before$CORE2_STRONG,

      CORE4_before =
        before$CORE4,

      CORE4_after =
        after$CORE4,

      delta_CORE4 =
        after$CORE4 -
        before$CORE4
    )
  }


  delta_dt <- rbindlist(
    delta_list
  )


  # ----------------------------------------------------------
  # Summary function for one endpoint
  # ----------------------------------------------------------

  summarize_endpoint <- function(
    endpoint,
    delta_vector
  ) {

    delta_vector <- as.numeric(
      delta_vector
    )

    n <- length(
      delta_vector
    )

    signfit <- exact_signflip(
      delta_vector
    )

    wilcoxfit <- wilcox_signed_rank(
      delta_vector
    )

    delta_sd <- if (
      n > 1
    ) {
      sd(
        delta_vector
      )
    } else {
      NA_real_
    }


    standardized <- if (
      is.finite(
        delta_sd
      ) &&
      delta_sd > 0
    ) {

      mean(
        delta_vector
      ) /
        delta_sd

    } else {

      NA_real_
    }


    data.table(
      scenario_id =
        scenario_id,

      scenario_role =
        scenario_role,

      hypothetical_EASC_patient =
        hypothetical_patient,

      endpoint =
        endpoint,

      reference_n =
        length(
          ref_idx
        ),

      n_paired =
        n,

      paired_ids =
        paste(
          paired_ids,
          collapse = ";"
        ),

      mean_delta =
        mean(
          delta_vector
        ),

      median_delta =
        median(
          delta_vector
        ),

      sd_delta =
        delta_sd,

      IQR_delta =
        IQR(
          delta_vector
        ),

      positive_n =
        sum(
          delta_vector > 0
        ),

      zero_n =
        sum(
          delta_vector == 0
        ),

      negative_n =
        sum(
          delta_vector < 0
        ),

      positive_fraction =
        mean(
          delta_vector > 0
        ),

      paired_standardized_effect =
        standardized,

      wilcoxon_W =
        wilcoxfit$statistic,

      wilcoxon_two_sided_P =
        wilcoxfit$p,

      wilcoxon_warning =
        wilcoxfit$warning,

      signflip_permutations =
        signfit$permutations,

      signflip_two_sided_P =
        signfit$p_two_sided,

      signflip_directional_P =
        signfit$p_directional,

      mean_direction =
        ifelse(
          mean(
            delta_vector
          ) > 0,
          "POSITIVE",
          ifelse(
            mean(
              delta_vector
            ) < 0,
            "NEGATIVE",
            "ZERO"
          )
        ),

      median_direction =
        ifelse(
          median(
            delta_vector
          ) > 0,
          "POSITIVE",
          ifelse(
            median(
              delta_vector
            ) < 0,
            "NEGATIVE",
            "ZERO"
          )
        ),

      majority_positive =
        (
          sum(
            delta_vector > 0
          ) >
            n / 2
        )
    )
  }


  result_core2 <- summarize_endpoint(
    "CORE2_STRONG",
    delta_dt$delta_CORE2
  )

  result_core4 <- summarize_endpoint(
    "CORE4",
    delta_dt$delta_CORE4
  )


  list(
    scores =
      score_dt,

    deltas =
      delta_dt,

    results =
      rbindlist(
        list(
          result_core2,
          result_core4
        )
      )
  )
}

# ============================================================
# 10. Run unresolved descriptive baseline
# ============================================================

cat("===== 9. UNRESOLVED TECHNICAL BASELINE =====\n\n")

base_row <- baseline_scenario[1]


baseline_run <- score_scenario(
  scenario_id =
    base_row$scenario_id,

  hypothetical_patient =
    "UNKNOWN",

  expected_reference_n =
    as.integer(
      base_row$technical_reference_n
    ),

  expected_pair_n =
    as.integer(
      base_row$paired_candidate_n
    ),

  scenario_role =
    "DESCRIPTIVE_UNRESOLVED_NOT_ESCC"
)


cat(
  "Descriptive unresolved baseline computed.\n"
)

cat(
  "Primary interpretation allowed: NO\n\n"
)

# ============================================================
# 11. Run ALL 22 hypothetical EASC scenarios
# ============================================================

cat("===== 10. RUN ALL 22 PRE-FROZEN EASC SCENARIOS =====\n\n")

all_runs <- list(
  baseline_run
)

for (
  i in seq_len(
    nrow(
      hypothetical
    )
  )
) {

  row <- hypothetical[i]

  patient <- as.character(
    row$hypothetical_EASC_patient
  )


  cat(
    sprintf(
      "[%02d/22] hypothetical EASC = %-4s  ref=%2d  pairs=%d\n",
      i,
      patient,
      as.integer(
        row$technical_reference_n
      ),
      as.integer(
        row$paired_candidate_n
      )
    )
  )


  all_runs[[
    length(
      all_runs
    ) + 1
  ]] <- score_scenario(
    scenario_id =
      row$scenario_id,

    hypothetical_patient =
      patient,

    expected_reference_n =
      as.integer(
        row$technical_reference_n
      ),

    expected_pair_n =
      as.integer(
        row$paired_candidate_n
      ),

    scenario_role =
      "SECONDARY_HISTOLOGY_SENSITIVITY"
  )
}

cat(
  "\n✓ All 22 hypothetical EASC scenarios completed.\n\n"
)

# ============================================================
# 12. Assemble outputs
# ============================================================

cat("===== 11. ASSEMBLE SECONDARY OUTPUTS =====\n\n")

scores_all <- rbindlist(
  lapply(
    all_runs,
    `[[`,
    "scores"
  ),
  use.names = TRUE
)

deltas_all <- rbindlist(
  lapply(
    all_runs,
    `[[`,
    "deltas"
  ),
  use.names = TRUE
)

results_all <- rbindlist(
  lapply(
    all_runs,
    `[[`,
    "results"
  ),
  use.names = TRUE
)


fwrite(
  scores_all,
  SCORES_OUT,
  bom = TRUE
)

fwrite(
  deltas_all,
  DELTAS_OUT,
  bom = TRUE
)

fwrite(
  results_all,
  RESULTS_OUT,
  bom = TRUE
)


cat(
  "Scenario score rows:",
  nrow(
    scores_all
  ),
  "\n"
)

cat(
  "Paired delta rows:",
  nrow(
    deltas_all
  ),
  "\n"
)

cat(
  "Scenario result rows:",
  nrow(
    results_all
  ),
  "\n\n"
)

# ============================================================
# 13. CORE2 robustness across the 22 hypothetical scenarios
#
# These 22 scenarios are NOT independent datasets.
# This is a robustness envelope, NOT meta-analysis.
# ============================================================

cat("===== 12. CORE2 HISTOLOGY-ROBUSTNESS ENVELOPE =====\n\n")

core2_hyp <- results_all[
  scenario_role ==
    "SECONDARY_HISTOLOGY_SENSITIVITY" &
  endpoint ==
    "CORE2_STRONG"
]


if (
  nrow(
    core2_hyp
  ) != 22
) {

  stop(
    "Expected exactly 22 CORE2 sensitivity results"
  )
}


mean_positive_n <- sum(
  core2_hyp$mean_delta > 0
)

median_positive_n <- sum(
  core2_hyp$median_delta > 0
)

majority_positive_n <- sum(
  core2_hyp$majority_positive
)

all_patient_positive_n <- sum(
  core2_hyp$positive_n ==
    core2_hyp$n_paired
)

signflip_p05_n <- sum(
  core2_hyp$signflip_two_sided_P <
    0.05,
  na.rm = TRUE
)

wilcox_p05_n <- sum(
  core2_hyp$wilcoxon_two_sided_P <
    0.05,
  na.rm = TRUE
)


pair_n4 <- sum(
  core2_hyp$n_paired ==
    4
)

pair_n5 <- sum(
  core2_hyp$n_paired ==
    5
)


robustness <- data.table(
  metric = c(
    "hypothetical_EASC_scenarios",
    "paired_n_equals_4",
    "paired_n_equals_5",
    "mean_CORE2_delta_positive",
    "median_CORE2_delta_positive",
    "majority_patients_positive",
    "all_paired_patients_positive",
    "exact_signflip_two_sided_P_lt_0.05",
    "wilcoxon_two_sided_P_lt_0.05",
    "minimum_mean_delta",
    "maximum_mean_delta",
    "minimum_median_delta",
    "maximum_median_delta",
    "minimum_positive_fraction",
    "maximum_positive_fraction"
  ),

  value = c(
    22,
    pair_n4,
    pair_n5,
    mean_positive_n,
    median_positive_n,
    majority_positive_n,
    all_patient_positive_n,
    signflip_p05_n,
    wilcox_p05_n,
    min(
      core2_hyp$mean_delta
    ),
    max(
      core2_hyp$mean_delta
    ),
    min(
      core2_hyp$median_delta
    ),
    max(
      core2_hyp$median_delta
    ),
    min(
      core2_hyp$positive_fraction
    ),
    max(
      core2_hyp$positive_fraction
    )
  )
)


fwrite(
  robustness,
  ROBUSTNESS_OUT,
  bom = TRUE
)

print(
  robustness
)

# ============================================================
# 14. Deterministic descriptive robustness classification
#
# This classification is explicitly SECONDARY.
#
# It does NOT replace the primary ESCC-only endpoint.
# ============================================================

all_mean_positive <- (
  mean_positive_n == 22
)

all_median_positive <- (
  median_positive_n == 22
)

all_majority_positive <- (
  majority_positive_n == 22
)


if (
  all_mean_positive &&
  all_median_positive &&
  all_majority_positive
) {

  secondary_status <- (
    "SECONDARY_DIRECTION_ROBUST_ACROSS_ALL_22_HISTOLOGY_SCENARIOS"
  )

} else if (
  mean_positive_n >= 18 &&
  median_positive_n >= 18 &&
  majority_positive_n >= 18
) {

  secondary_status <- (
    "SECONDARY_DIRECTION_MOSTLY_ROBUST_ACROSS_HISTOLOGY_SCENARIOS"
  )

} else if (
  mean_positive_n > 11
) {

  secondary_status <- (
    "SECONDARY_DIRECTION_PARTIALLY_ROBUST"
  )

} else {

  secondary_status <- (
    "SECONDARY_DIRECTION_NOT_ROBUST_TO_HISTOLOGY_UNCERTAINTY"
  )
}


cat("\n")

cat(
  "Secondary robustness status:\n"
)

cat(
  secondary_status,
  "\n\n"
)

cat(
  "IMPORTANT:\n"
)

cat(
  "This status is SECONDARY only.\n"
)

cat(
  "Primary ESCC-only validation remains unavailable.\n\n"
)

# ============================================================
# 15. Statistical-resolution warning
# ============================================================

cat("===== 13. SMALL-N EXACT-TEST RESOLUTION =====\n\n")

cat(
  "With n=5 paired patients:\n"
)

cat(
  "  minimum possible two-sided exact sign-flip P = 2/32 = 0.0625\n"
)

cat(
  "With n=4 paired patients:\n"
)

cat(
  "  minimum possible two-sided exact sign-flip P = 2/16 = 0.125\n\n"
)

cat(
  "Therefore P<0.05 is not an attainable requirement\n"
)

cat(
  "for the frozen exact sign-flip test in these scenarios.\n\n"
)

# ============================================================
# 16. Frozen project integrity
# ============================================================

cat("===== 14. FROZEN STEP19-21 INTEGRITY =====\n\n")

changes <- frozen_integrity()

cat(
  "Frozen original files changed:",
  length(
    changes
  ),
  "\n"
)


if (
  length(
    changes
  ) > 0
) {

  print(
    changes
  )

  stop(
    "STEP22A_FROZEN_PROJECT_CHANGED"
  )
}


cat(
  "✓ Frozen Step19-21 inventory unchanged.\n\n"
)

# ============================================================
# 17. Summary file
# ============================================================

base_core2 <- results_all[
  scenario_id ==
    "S00_UNRESOLVED" &
  endpoint ==
    "CORE2_STRONG"
]


summary_lines <- c(
  "STEP22A-7A SECONDARY HISTOLOGY SENSITIVITY CORE2",
  "",
  "GOVERNANCE",
  "Primary ESCC-only analysis: BLOCKED",
  "Secondary histology sensitivity: AUTHORIZED",
  "Scenario selection: FORBIDDEN",
  "All pre-frozen hypothetical EASC scenarios run: 22/22",
  "Author nonresponse treated as histology evidence: NO",
  "",
  "UNRESOLVED DESCRIPTIVE BASELINE",
  paste0(
    "paired n: ",
    base_core2$n_paired
  ),
  paste0(
    "mean delta CORE2: ",
    format(
      base_core2$mean_delta,
      digits = 6
    )
  ),
  paste0(
    "median delta CORE2: ",
    format(
      base_core2$median_delta,
      digits = 6
    )
  ),
  paste0(
    "positive patients: ",
    base_core2$positive_n,
    "/",
    base_core2$n_paired
  ),
  "Interpretation: descriptive only; unknown EASC retained",
  "",
  "22-SCENARIO HISTOLOGY ROBUSTNESS",
  paste0(
    "Mean delta positive: ",
    mean_positive_n,
    "/22"
  ),
  paste0(
    "Median delta positive: ",
    median_positive_n,
    "/22"
  ),
  paste0(
    "Majority paired patients positive: ",
    majority_positive_n,
    "/22"
  ),
  paste0(
    "All paired patients positive: ",
    all_patient_positive_n,
    "/22"
  ),
  paste0(
    "Paired n=4 scenarios: ",
    pair_n4
  ),
  paste0(
    "Paired n=5 scenarios: ",
    pair_n5
  ),
  paste0(
    "Minimum mean delta: ",
    format(
      min(
        core2_hyp$mean_delta
      ),
      digits = 6
    )
  ),
  paste0(
    "Maximum mean delta: ",
    format(
      max(
        core2_hyp$mean_delta
      ),
      digits = 6
    )
  ),
  paste0(
    "Minimum median delta: ",
    format(
      min(
        core2_hyp$median_delta
      ),
      digits = 6
    )
  ),
  paste0(
    "Maximum median delta: ",
    format(
      max(
        core2_hyp$median_delta
      ),
      digits = 6
    )
  ),
  "",
  paste0(
    "Secondary robustness classification: ",
    secondary_status
  ),
  "",
  "STATISTICAL RESOLUTION",
  "n=5 exact two-sided sign-flip minimum P: 0.0625",
  "n=4 exact two-sided sign-flip minimum P: 0.125",
  "",
  "INTERPRETATION FIREWALL",
  "These 22 scenarios are not 22 independent replications.",
  "Do not meta-analyze them.",
  "Do not choose the most favorable scenario.",
  "Do not promote a hypothetical EASC scenario to primary.",
  "Primary ESCC-only validation remains BLOCKED.",
  "",
  "Frozen Step19-21 modified: NO",
  "Manuscript modified: NO"
)


writeLines(
  summary_lines,
  SUMMARY_OUT
)

# ============================================================
# 18. Session info
# ============================================================

sink(
  file.path(
    LOG_DIR,
    "STEP22A_07A_sessionInfo.txt"
  )
)

print(
  sessionInfo()
)

sink()

# ============================================================
# FINAL
# ============================================================

cat("\n")
cat("====================================================================\n")
cat("STEP22A-7A COMPLETE\n")
cat("====================================================================\n\n")

cat(
  "PRIMARY ESCC-ONLY             : BLOCKED\n"
)

cat(
  "SECONDARY SENSITIVITY         : AUTHORIZED / COMPLETED\n"
)

cat(
  "HYPOTHETICAL EASC SCENARIOS   : 22 / 22\n"
)

cat(
  "UNRESOLVED BASELINE           : 1 DESCRIPTIVE ONLY\n\n"
)

cat(
  "PAIRED SAMPLE SIZE RANGE      :",
  min(
    core2_hyp$n_paired
  ),
  "-",
  max(
    core2_hyp$n_paired
  ),
  "\n\n"
)

cat(
  "CORE2 MEAN DELTA POSITIVE     :",
  mean_positive_n,
  "/ 22\n"
)

cat(
  "CORE2 MEDIAN DELTA POSITIVE   :",
  median_positive_n,
  "/ 22\n"
)

cat(
  "MAJORITY PATIENTS POSITIVE    :",
  majority_positive_n,
  "/ 22\n"
)

cat(
  "ALL PAIRED PATIENTS POSITIVE  :",
  all_patient_positive_n,
  "/ 22\n\n"
)

cat(
  "MEAN DELTA RANGE              :",
  sprintf(
    "%.4f to %.4f",
    min(
      core2_hyp$mean_delta
    ),
    max(
      core2_hyp$mean_delta
    )
  ),
  "\n"
)

cat(
  "MEDIAN DELTA RANGE            :",
  sprintf(
    "%.4f to %.4f",
    min(
      core2_hyp$median_delta
    ),
    max(
      core2_hyp$median_delta
    )
  ),
  "\n"
)

cat(
  "POSITIVE FRACTION RANGE       :",
  sprintf(
    "%.3f to %.3f",
    min(
      core2_hyp$positive_fraction
    ),
    max(
      core2_hyp$positive_fraction
    )
  ),
  "\n\n"
)

cat(
  "EXACT SIGN-FLIP P<0.05        :",
  signflip_p05_n,
  "/ 22\n"
)

cat(
  "WILCOXON P<0.05               :",
  wilcox_p05_n,
  "/ 22\n\n"
)

cat(
  "SECONDARY ROBUSTNESS STATUS:\n"
)

cat(
  secondary_status,
  "\n\n"
)

cat(
  "IMPORTANT:\n"
)

cat(
  "  These 22 scenarios are a sensitivity envelope,\n"
)

cat(
  "  NOT 22 independent validation cohorts.\n"
)

cat(
  "  Do NOT select the best scenario.\n"
)

cat(
  "  Do NOT meta-analyze the 22 scenarios.\n"
)

cat(
  "  Do NOT promote any scenario to primary.\n\n"
)

cat(
  "PRIMARY ESCC VALIDATION       : STILL UNAVAILABLE\n"
)

cat(
  "MANUSCRIPT MODIFIED           : NO\n"
)

cat(
  "FROZEN STEP19-21 MODIFIED     : NO\n\n"
)

cat(
  "STOP HERE.\n"
)

cat(
  "REVIEW THE 22-SCENARIO DIRECTIONAL ROBUSTNESS BEFORE STEP22A-7B.\n\n"
)
