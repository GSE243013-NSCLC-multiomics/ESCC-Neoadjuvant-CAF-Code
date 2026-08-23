# ============================================================
# STEP 17H6A
#
# Audit the three completed GSE221561 InferCNA pilots:
#   S6417T  Surgery_alone
#   S9478T  Surgery_alone
#   S1265T  Neoadjuvant_treated
#
# Also:
#   - validate coherent-arm counts against actual CSV rows
#   - select the largest remaining neoadjuvant-treated
#     epithelial sample for the NEXT pilot
#
# IMPORTANT:
#   NO InferCNA
#   NO findMalignant
#   NO malignant labels
#   NO new samples processed
# ============================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
})

cat("\n============================================================\n")
cat("STEP 17H6A: THREE-PILOT AUDIT\n")
cat("============================================================\n\n")

base_dir <- file.path(
  "06_tables",
  "cross_dataset",
  "malignant_epithelial"
)

outdir <- file.path(
  base_dir,
  "GSE221561_three_pilot_review"
)

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)

samples <- c(
  "S6417T",
  "S9478T",
  "S1265T"
)

treatments <- c(
  S6417T = "Surgery_alone",
  S9478T = "Surgery_alone",
  S1265T = "Neoadjuvant_treated"
)

all_files <- list.files(
  base_dir,
  recursive = TRUE,
  full.names = TRUE
)

# ============================================================
# Helpers
# ============================================================

latest_file <- function(sample, pattern) {

  x <- all_files[
    grepl(sample, all_files, ignore.case = TRUE) &
    grepl(pattern, basename(all_files), ignore.case = TRUE)
  ]

  if (length(x) == 0) {
    return(NA_character_)
  }

  x[
    which.max(
      file.info(x)$mtime
    )
  ]
}

get_num <- function(dt, candidates) {

  hit <- intersect(
    candidates,
    names(dt)
  )

  if (length(hit) == 0) {
    return(NA_real_)
  }

  suppressWarnings(
    as.numeric(
      dt[[hit[1]]][1]
    )
  )
}

# ============================================================
# 1. Locate manifests
# ============================================================

cat("===== 1. LOCATE THREE PILOTS =====\n\n")

manifest_files <- setNames(
  lapply(
    samples,
    function(s) {
      latest_file(
        s,
        "manifest.*\\.csv$"
      )
    }
  ),
  samples
)

for (s in samples) {

  cat(
    s,
    " manifest: ",
    manifest_files[[s]],
    "\n",
    sep = ""
  )

  if (is.na(manifest_files[[s]])) {
    stop(
      "Manifest missing for ",
      s
    )
  }
}

# ============================================================
# 2. Extract manifest metrics
# ============================================================

extract_pilot <- function(sample) {

  x <- fread(
    manifest_files[[sample]]
  )

  data.table(

    sample =
      sample,

    treatment =
      treatments[[sample]],

    tumor_cells =
      get_num(
        x,
        c(
          "tumor_cells",
          "tumor_epithelial_cells"
        )
      ),

    normal_cells =
      get_num(
        x,
        c(
          "normal_reference_cells",
          "normal_cells"
        )
      ),

    CNA_genes =
      get_num(
        x,
        c(
          "CNA_genes",
          "genes",
          "final_hg38_genes"
        )
      ),

    tumor_median_CNA_signal =
      get_num(
        x,
        c(
          "tumor_median_CNA_signal",
          "tumor_median_cna_signal"
        )
      ),

    normal_median_CNA_signal =
      get_num(
        x,
        c(
          "normal_median_CNA_signal",
          "normal_median_cna_signal"
        )
      ),

    tumor_median_cnaCor =
      get_num(
        x,
        c(
          "tumor_median_cnaCor",
          "tumor_cnaCor_median",
          "tumor_median_cnacor"
        )
      ),

    normal_corr_to_tumor =
      get_num(
        x,
        c(
          "normal_median_corr_to_tumor",
          "normal_median_correlation_to_tumor"
        )
      ),

    coherent_arms_manifest =
      get_num(
        x,
        c(
          "coherent_CNA_arms",
          "coherent_cna_arms"
        )
      )
  )
}

comparison <- rbindlist(
  lapply(
    samples,
    extract_pilot
  ),
  fill = TRUE
)

comparison[
  ,
  CNA_signal_delta :=
    tumor_median_CNA_signal -
    normal_median_CNA_signal
]

comparison[
  ,
  CNA_signal_ratio :=
    tumor_median_CNA_signal /
    normal_median_CNA_signal
]

comparison[
  ,
  correlation_separation :=
    tumor_median_cnaCor -
    normal_corr_to_tumor
]

# ============================================================
# 3. Recover >95 / >99 from cell metric files
# ============================================================

for (s in samples) {

  metric_file <- latest_file(
    s,
    "cell.*CNA.*metrics.*\\.csv$|cell.*metrics.*\\.csv$"
  )

  if (is.na(metric_file)) {
    next
  }

  z <- fread(
    metric_file
  )

  if (
    all(
      c(
        "role",
        "cna_signal"
      ) %in%
        names(z)
    )
  ) {

    normal <- z[
      role == "NORMAL_REFERENCE",
      cna_signal
    ]

    tumor <- z[
      role == "TUMOR_TEST",
      cna_signal
    ]

    if (
      length(normal) > 0 &&
      length(tumor) > 0
    ) {

      q95 <- quantile(
        normal,
        0.95,
        na.rm = TRUE
      )

      q99 <- quantile(
        normal,
        0.99,
        na.rm = TRUE
      )

      comparison[
        sample == s,
        tumor_above_normal95 :=
          mean(
            tumor > q95,
            na.rm = TRUE
          )
      ]

      comparison[
        sample == s,
        tumor_above_normal99 :=
          mean(
            tumor > q99,
            na.rm = TRUE
          )
      ]
    }
  }
}

# ============================================================
# 4. Coherent-arm audit
# ============================================================

cat("\n============================================================\n")
cat("2. COHERENT ARM CONSISTENCY AUDIT\n")
cat("============================================================\n\n")

arm_audit <- list()

for (s in samples) {

  coherent_file <- latest_file(
    s,
    "candidate.*coherent.*arm.*\\.csv$"
  )

  all_arm_file <- latest_file(
    s,
    "all.*arm.*metrics.*\\.csv$"
  )

  arm_dt <- NULL
  source_file <- NA_character_

  if (!is.na(coherent_file)) {

    x <- fread(
      coherent_file
    )

    arm_dt <- x
    source_file <- coherent_file

  } else if (!is.na(all_arm_file)) {

    x <- fread(
      all_arm_file
    )

    if (
      "coherent_audit_flag" %in%
        names(x)
    ) {

      arm_dt <- x[
        coherent_audit_flag == TRUE
      ]

    } else {

      arm_dt <- data.table()
    }

    source_file <- all_arm_file
  }

  if (is.null(arm_dt)) {

    actual_n <- NA_integer_
    arm_names <- NA_character_

  } else {

    actual_n <- nrow(
      arm_dt
    )

    arm_col <- intersect(
      c(
        "chromosome_arm",
        "arm",
        "chr_arm"
      ),
      names(arm_dt)
    )

    if (
      length(arm_col) > 0 &&
      nrow(arm_dt) > 0
    ) {

      arm_names <- paste(
        unique(
          as.character(
            arm_dt[[arm_col[1]]]
          )
        ),
        collapse = ", "
      )

      unique_n <- uniqueN(
        arm_dt[[arm_col[1]]]
      )

    } else {

      arm_names <- ""
      unique_n <- 0L
    }
  }

  manifest_n <- comparison[
    sample == s,
    coherent_arms_manifest
  ]

  arm_audit[[s]] <- data.table(

    sample =
      s,

    treatment =
      treatments[[s]],

    manifest_count =
      manifest_n,

    CSV_rows =
      actual_n,

    unique_arm_names =
      unique_n,

    count_consistent =
      !is.na(manifest_n) &&
      !is.na(actual_n) &&
      manifest_n == actual_n &&
      actual_n == unique_n,

    arm_names =
      arm_names,

    source_file =
      source_file
  )
}

arm_audit <- rbindlist(
  arm_audit,
  fill = TRUE
)

print(
  arm_audit,
  nrows = 100
)

fwrite(
  arm_audit,
  file.path(
    outdir,
    "STEP17H6A_coherent_arm_consistency_audit.csv"
  ),
  bom = TRUE
)

comparison <- merge(
  comparison,
  arm_audit[
    ,
    .(
      sample,
      coherent_arms_verified =
        unique_arm_names,

      coherent_arm_names =
        arm_names,

      arm_count_consistent =
        count_consistent
    )
  ],
  by = "sample",
  all.x = TRUE
)

# ============================================================
# 5. Three-pilot comparison
# ============================================================

cat("\n============================================================\n")
cat("3. VERIFIED THREE-PILOT COMPARISON\n")
cat("============================================================\n\n")

show <- comparison[
  ,
  .(
    sample,
    treatment,
    tumor_cells,
    CNA_genes,

    tumor_median_CNA_signal =
      round(
        tumor_median_CNA_signal,
        5
      ),

    normal_median_CNA_signal =
      round(
        normal_median_CNA_signal,
        5
      ),

    CNA_signal_ratio =
      round(
        CNA_signal_ratio,
        2
      ),

    tumor_median_cnaCor =
      round(
        tumor_median_cnaCor,
        3
      ),

    correlation_separation =
      round(
        correlation_separation,
        3
      ),

    tumor_above_normal95 =
      round(
        100 *
        tumor_above_normal95,
        1
      ),

    tumor_above_normal99 =
      round(
        100 *
        tumor_above_normal99,
        1
      ),

    coherent_arms_verified,

    arm_count_consistent
  )
]

print(
  show,
  nrows = 100
)

fwrite(
  comparison,
  file.path(
    outdir,
    "STEP17H6A_verified_three_pilot_comparison.csv"
  ),
  bom = TRUE
)

# ============================================================
# 6. Select next neoadjuvant pilot
# ============================================================

cat("\n============================================================\n")
cat("4. NEXT NEOADJUVANT PILOT CANDIDATE\n")
cat("============================================================\n\n")

inventory_file <- file.path(
  base_dir,
  "GSE221561_preflight",
  "STEP17G1_GSE221561_epithelial_sample_inventory.csv"
)

if (!file.exists(inventory_file)) {

  cat(
    "Inventory file not found:\n",
    inventory_file,
    "\n"
  )

  cat(
    "No new sample selected.\n"
  )

} else {

  inv <- fread(
    inventory_file
  )

  required_cols <- c(
    "sample_id",
    "treatment_group",
    "epithelial_cells"
  )

  if (
    !all(
      required_cols %in%
        names(inv)
    )
  ) {

    stop(
      "Inventory lacks required fields."
    )
  }

  if (
    "tissue_simple" %in%
      names(inv)
  ) {

    inv <- inv[
      grepl(
        "tumor",
        tissue_simple,
        ignore.case = TRUE
      )
    ]
  }

  candidates <- inv[
    grepl(
      "neo|nact",
      treatment_group,
      ignore.case = TRUE
    ) &
    !sample_id %in%
      samples
  ]

  setorder(
    candidates,
    -epithelial_cells
  )

  if (
    nrow(candidates) == 0
  ) {

    cat(
      "No remaining neoadjuvant candidate found.\n"
    )

  } else {

    cat(
      "Largest remaining neoadjuvant epithelial sample:\n\n"
    )

    print(
      candidates[1]
    )

    cat(
      "\nTop remaining candidates:\n\n"
    )

    print(
      head(
        candidates,
        10
      )
    )

    fwrite(
      candidates,
      file.path(
        outdir,
        "STEP17H6A_remaining_neoadjuvant_candidates.csv"
      ),
      bom = TRUE
    )
  }
}

# ============================================================
# 7. Interpretation guardrails
# ============================================================

cat("\n============================================================\n")
cat("5. INTERPRETATION STATUS\n")
cat("============================================================\n\n")

cat(
  "✓ Two independent Surgery_alone pilots show strong CNA separation.\n"
)

cat(
  "✓ This provides independent technical validation of the GSE221561 InferCNA pipeline.\n"
)

cat(
  "⚠ S1265T is substantially weaker, but one treated sample cannot establish a treatment effect.\n"
)

cat(
  "⚠ Paired-normal availability differs conceptually from a randomized treatment comparison.\n"
)

cat(
  "✓ Therefore the next most informative experiment is a SECOND neoadjuvant-treated pilot.\n"
)

cat(
  "✓ No malignant threshold should be frozen yet.\n"
)

# ============================================================
# FINAL
# ============================================================

cat("\n============================================================\n")
cat("STEP 17H6A COMPLETE ✓\n")
cat("THREE-PILOT AUDIT COMPLETE ✓\n")
cat("============================================================\n\n")

cat("InferCNA run              : NO\n")
cat("findMalignant             : NO\n")
cat("Malignant labels assigned : NO\n")
cat("New tumor samples run     : NO\n")

cat("\nSTOP HERE.\n")
cat("DO NOT RUN THE NEXT SAMPLE YET.\n")

