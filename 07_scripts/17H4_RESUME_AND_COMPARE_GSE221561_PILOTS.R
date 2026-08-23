# ============================================================
# ESCC Neoadjuvant Project
# STEP 17H4
#
# Resume after restart.
#
# Compare EXISTING GSE221561 InferCNA pilots:
#   S6417T  = Surgery_alone
#   S1265T  = Neoadjuvant_treated
#
# IMPORTANT:
#   - DO NOT rerun InferCNA
#   - DO NOT run findMalignant()
#   - DO NOT assign malignant labels
#   - DO NOT process additional samples
# ============================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
})

cat("\n============================================================\n")
cat("STEP 17H4: RESUME + COMPARE GSE221561 CNV PILOTS\n")
cat("============================================================\n\n")

base_dir <- file.path(
  "06_tables",
  "cross_dataset",
  "malignant_epithelial"
)

outdir <- file.path(
  base_dir,
  "GSE221561_two_pilot_review"
)

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ============================================================
# 1. Find existing outputs robustly
# ============================================================

all_files <- list.files(
  base_dir,
  recursive = TRUE,
  full.names = TRUE
)

cat("===== 1. SEARCH EXISTING PILOT OUTPUTS =====\n\n")

find_sample_files <- function(sample) {

  hits <- all_files[
    grepl(
      sample,
      all_files,
      ignore.case = TRUE
    )
  ]

  if (length(hits) == 0) {
    stop(
      "No files found for ",
      sample
    )
  }

  list(
    all = hits,

    manifest = hits[
      grepl(
        "manifest.*\\.csv$",
        basename(hits),
        ignore.case = TRUE
      )
    ],

    summary = hits[
      grepl(
        "tumor.*normal.*summary.*\\.csv$|summary.*tumor.*normal.*\\.csv$",
        basename(hits),
        ignore.case = TRUE
      )
    ],

    cell_metrics = hits[
      grepl(
        "cell.*CNA.*metrics.*\\.csv$|cell.*metrics.*\\.csv$",
        basename(hits),
        ignore.case = TRUE
      )
    ],

    arms = hits[
      grepl(
        "candidate.*coherent.*arm.*\\.csv$|coherent.*arm.*\\.csv$",
        basename(hits),
        ignore.case = TRUE
      )
    ],

    cna = hits[
      grepl(
        "CNA_all_cells\\.rds$",
        basename(hits),
        ignore.case = TRUE
      )
    ]
  )
}

pick_latest <- function(x, description, sample) {

  if (length(x) == 0) {
    return(NA_character_)
  }

  x[
    which.max(
      file.info(x)$mtime
    )
  ]
}

s6417 <- find_sample_files("S6417T")
s1265 <- find_sample_files("S1265T")

s6417_manifest <- pick_latest(
  s6417$manifest,
  "manifest",
  "S6417T"
)

s1265_manifest <- pick_latest(
  s1265$manifest,
  "manifest",
  "S1265T"
)

s6417_summary <- pick_latest(
  s6417$summary,
  "summary",
  "S6417T"
)

s1265_summary <- pick_latest(
  s1265$summary,
  "summary",
  "S1265T"
)

s6417_cna <- pick_latest(
  s6417$cna,
  "CNA RDS",
  "S6417T"
)

s1265_cna <- pick_latest(
  s1265$cna,
  "CNA RDS",
  "S1265T"
)

cat("S6417T manifest : ", s6417_manifest, "\n", sep = "")
cat("S6417T summary  : ", s6417_summary, "\n", sep = "")
cat("S6417T CNA      : ", s6417_cna, "\n\n", sep = "")

cat("S1265T manifest : ", s1265_manifest, "\n", sep = "")
cat("S1265T summary  : ", s1265_summary, "\n", sep = "")
cat("S1265T CNA      : ", s1265_cna, "\n", sep = "")

# ============================================================
# 2. Critical resume check
# ============================================================

cat("\n============================================================\n")
cat("2. RESUME CHECK\n")
cat("============================================================\n\n")

critical <- c(
  s6417_manifest,
  s6417_cna,
  s1265_manifest,
  s1265_cna
)

if (any(is.na(critical))) {

  cat("One or more critical completed-pilot files are missing.\n\n")

  cat("STOP HERE.\n")
  cat("DO NOT automatically rerun InferCNA.\n")

  stop(
    "Pilot completion cannot yet be confirmed."
  )
}

cat("✓ S6417T CNA output exists\n")
cat("✓ S1265T CNA output exists\n")
cat("✓ Both manifests exist\n")
cat("✓ No InferCNA rerun is required\n")

# ============================================================
# 3. Load manifests
# ============================================================

cat("\n============================================================\n")
cat("3. LOAD PILOT MANIFESTS\n")
cat("============================================================\n\n")

m6417 <- fread(
  s6417_manifest
)

m1265 <- fread(
  s1265_manifest
)

cat("S6417T manifest:\n")
print(m6417)

cat("\nS1265T manifest:\n")
print(m1265)

# ============================================================
# 4. Flexible metric extractor
# ============================================================

get_num <- function(dt, candidates) {

  nm <- intersect(
    candidates,
    names(dt)
  )

  if (length(nm) == 0) {
    return(NA_real_)
  }

  suppressWarnings(
    as.numeric(
      dt[[nm[1]]][1]
    )
  )
}

get_chr <- function(dt, candidates) {

  nm <- intersect(
    candidates,
    names(dt)
  )

  if (length(nm) == 0) {
    return(NA_character_)
  }

  as.character(
    dt[[nm[1]]][1]
  )
}

extract_manifest <- function(
  dt,
  sample,
  treatment
) {

  data.table(

    sample = sample,

    treatment = treatment,

    tumor_cells =
      get_num(
        dt,
        c(
          "tumor_cells",
          "tumor_epithelial_cells"
        )
      ),

    normal_reference_cells =
      get_num(
        dt,
        c(
          "normal_reference_cells",
          "normal_cells"
        )
      ),

    CNA_genes =
      get_num(
        dt,
        c(
          "CNA_genes",
          "genes",
          "final_hg38_genes"
        )
      ),

    runtime_minutes =
      get_num(
        dt,
        c(
          "runtime_minutes",
          "runtime_min"
        )
      ),

    tumor_median_CNA_signal =
      get_num(
        dt,
        c(
          "tumor_median_CNA_signal",
          "tumor_median_cna_signal"
        )
      ),

    normal_median_CNA_signal =
      get_num(
        dt,
        c(
          "normal_median_CNA_signal",
          "normal_median_cna_signal"
        )
      ),

    tumor_median_cnaCor =
      get_num(
        dt,
        c(
          "tumor_median_cnaCor",
          "tumor_median_cnacor",
          "tumor_cnaCor_median"
        )
      ),

    normal_median_corr_to_tumor =
      get_num(
        dt,
        c(
          "normal_median_corr_to_tumor",
          "normal_median_correlation_to_tumor"
        )
      ),

    coherent_CNA_arms =
      get_num(
        dt,
        c(
          "coherent_CNA_arms",
          "coherent_cna_arms"
        )
      )
  )
}

comparison <- rbindlist(
  list(

    extract_manifest(
      m6417,
      "S6417T",
      "Surgery_alone"
    ),

    extract_manifest(
      m1265,
      "S1265T",
      "Neoadjuvant_treated"
    )
  ),
  fill = TRUE
)

# ============================================================
# 5. Recompute derived metrics consistently
# ============================================================

comparison[
  ,
  CNA_signal_delta :=
    tumor_median_CNA_signal -
    normal_median_CNA_signal
]

comparison[
  ,
  CNA_signal_ratio :=
    fifelse(
      !is.na(normal_median_CNA_signal) &
      normal_median_CNA_signal != 0,
      tumor_median_CNA_signal /
        normal_median_CNA_signal,
      NA_real_
    )
]

comparison[
  ,
  correlation_separation :=
    tumor_median_cnaCor -
    normal_median_corr_to_tumor
]

cat("\n============================================================\n")
cat("4. GSE221561 TWO-PILOT COMPARISON\n")
cat("============================================================\n\n")

print(
  comparison,
  nrows = 100
)

fwrite(
  comparison,
  file.path(
    outdir,
    "STEP17H4_GSE221561_S6417T_vs_S1265T_manifest_comparison.csv"
  ),
  bom = TRUE
)

# ============================================================
# 6. Read summary tables if available
# ============================================================

cat("\n============================================================\n")
cat("5. ORIGINAL TUMOR/NORMAL SUMMARY TABLES\n")
cat("============================================================\n\n")

if (!is.na(s6417_summary)) {

  cat("--- S6417T ---\n")

  x <- fread(
    s6417_summary
  )

  print(
    x,
    nrows = 100
  )

} else {

  cat("S6417T summary CSV not found; manifest retained.\n")
}

cat("\n")

if (!is.na(s1265_summary)) {

  cat("--- S1265T ---\n")

  y <- fread(
    s1265_summary
  )

  print(
    y,
    nrows = 100
  )

} else {

  cat("S1265T summary CSV not found; manifest retained.\n")
}

# ============================================================
# 7. Inspect coherent-arm files
# ============================================================

cat("\n============================================================\n")
cat("6. COHERENT CNA ARM REVIEW\n")
cat("============================================================\n\n")

read_arm_file <- function(files, sample) {

  if (length(files) == 0) {

    cat(
      sample,
      ": no coherent-arm CSV located\n",
      sep = ""
    )

    return(
      data.table()
    )
  }

  f <- files[
    which.max(
      file.info(files)$mtime
    )
  ]

  cat(
    sample,
    " arm file: ",
    f,
    "\n",
    sep = ""
  )

  x <- fread(f)

  if (nrow(x) == 0) {

    cat(
      sample,
      ": 0 candidate coherent arms\n"
    )

  } else {

    print(
      x,
      nrows = 100
    )
  }

  x
}

a6417 <- read_arm_file(
  s6417$arms,
  "S6417T"
)

cat("\n")

a1265 <- read_arm_file(
  s1265$arms,
  "S1265T"
)

# ============================================================
# 8. CNA RDS dimension validation
# ============================================================

cat("\n============================================================\n")
cat("7. CNA OBJECT VALIDATION\n")
cat("============================================================\n\n")

cna6417 <- readRDS(
  s6417_cna
)

cna1265 <- readRDS(
  s1265_cna
)

cat(
  "S6417T CNA: ",
  nrow(cna6417),
  " genes x ",
  ncol(cna6417),
  " cells\n",
  sep = ""
)

cat(
  "S1265T CNA: ",
  nrow(cna1265),
  " genes x ",
  ncol(cna1265),
  " cells\n",
  sep = ""
)

cat(
  sprintf(
    "S6417T range: %.4f to %.4f\n",
    min(cna6417, na.rm = TRUE),
    max(cna6417, na.rm = TRUE)
  )
)

cat(
  sprintf(
    "S1265T range: %.4f to %.4f\n",
    min(cna1265, na.rm = TRUE),
    max(cna1265, na.rm = TRUE)
  )
)

rm(
  cna6417,
  cna1265
)

gc(verbose = FALSE)

# ============================================================
# 9. Interpretation guardrails
# ============================================================

cat("\n============================================================\n")
cat("8. INTERPRETATION GUARDRAILS\n")
cat("============================================================\n\n")

cat(
  "1. This step compares existing pilots only.\n"
)

cat(
  "2. Treatment group differences are descriptive, NOT causal.\n"
)

cat(
  "3. A larger CNA signal does NOT automatically mean malignant.\n"
)

cat(
  "4. cnaCor must be interpreted together with CNA signal and chromosome-scale coherence.\n"
)

cat(
  "5. No threshold is being used here for malignant calling.\n"
)

cat(
  "6. Do not run findMalignant() yet.\n"
)

cat(
  "7. Do not process the remaining GSE221561 samples yet.\n"
)

# ============================================================
# 10. Final concise comparison
# ============================================================

cat("\n============================================================\n")
cat("9. REVIEW GATE\n")
cat("============================================================\n\n")

for (i in seq_len(nrow(comparison))) {

  z <- comparison[i]

  cat(
    "--- ",
    z$sample,
    " [",
    z$treatment,
    "] ---\n",
    sep = ""
  )

  cat(
    "Tumor cells                 : ",
    z$tumor_cells,
    "\n",
    sep = ""
  )

  cat(
    "Normal reference cells      : ",
    z$normal_reference_cells,
    "\n",
    sep = ""
  )

  cat(
    "CNA genes                   : ",
    z$CNA_genes,
    "\n",
    sep = ""
  )

  cat(
    sprintf(
      "Tumor median CNA signal     : %.5f\n",
      z$tumor_median_CNA_signal
    )
  )

  cat(
    sprintf(
      "Normal median CNA signal    : %.5f\n",
      z$normal_median_CNA_signal
    )
  )

  cat(
    sprintf(
      "CNA signal delta            : %.5f\n",
      z$CNA_signal_delta
    )
  )

  cat(
    sprintf(
      "CNA signal ratio            : %.3f\n",
      z$CNA_signal_ratio
    )
  )

  cat(
    sprintf(
      "Tumor median cnaCor         : %.3f\n",
      z$tumor_median_cnaCor
    )
  )

  cat(
    sprintf(
      "Normal correlation to tumor : %.3f\n",
      z$normal_median_corr_to_tumor
    )
  )

  cat(
    sprintf(
      "Correlation separation      : %.3f\n",
      z$correlation_separation
    )
  )

  cat(
    "Coherent CNA arms           : ",
    z$coherent_CNA_arms,
    "\n\n",
    sep = ""
  )
}

# ============================================================
# FINAL
# ============================================================

cat("============================================================\n")
cat("STEP 17H4 COMPLETE ✓\n")
cat("RESTART RECOVERY COMPLETE ✓\n")
cat("============================================================\n\n")

cat("S6417T InferCNA rerun : NO ✓\n")
cat("S1265T InferCNA rerun : NO ✓\n")
cat("findMalignant()       : NO ✓\n")
cat("Malignant labels      : NO ✓\n")
cat("Other samples run     : NO ✓\n")

cat("\nSTOP HERE.\n")
cat("REVIEW THE TWO PILOTS BEFORE ANY EXPANSION.\n")

