# ============================================================
# ESCC Neoadjuvant Project
# STEP 17A
#
# Malignant epithelial strategy audit across:
#   GSE221561
#   GSE197677
#   GSE160269
#
# PURPOSE:
#   quantify epithelial cells by sample/patient/tissue
#   assess normal epithelial reference availability
#   identify tumor-normal paired subjects
#
# IMPORTANT:
#   NO malignant calling yet
#   NO CNV analysis yet
#   NO cells removed
#   NO DE
#   NO normalization
#   NO integration
# ============================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
})

cat("\n============================================================\n")
cat("STEP 17A: MALIGNANT EPITHELIAL STRATEGY AUDIT\n")
cat("============================================================\n\n")

outdir <- "06_tables/cross_dataset/malignant_epithelial"

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ============================================================
# INPUTS
# ============================================================

f221 <- paste0(
  "06_tables/GSE221561/metadata/",
  "GSE221561_cell_metadata_standardized.csv.gz"
)

f197 <- paste0(
  "06_tables/cross_dataset/",
  "STEP16A_GSE197677_sample_by_broad_celltype_REBUILT.csv"
)

f160 <- paste0(
  "06_tables/GSE160269/metadata/",
  "GSE160269_merged_patient_tissue_celltype_counts.csv"
)

inputs <- c(
  f221,
  f197,
  f160
)

missing <- inputs[
  !file.exists(inputs)
]

if (length(missing) > 0) {

  stop(
    "Missing input files:\n",
    paste(
      missing,
      collapse = "\n"
    )
  )
}

cat("===== 1. INPUT FILES =====\n\n")

for (f in inputs) {
  cat("✓ ", f, "\n", sep = "")
}

# ============================================================
# HELPER
# ============================================================

standardize_tissue <- function(x) {

  y <- tolower(
    trimws(
      as.character(x)
    )
  )

  fcase(

    grepl(
      "tumor|tumour|scc",
      y
    ),
    "Tumor",

    grepl(
      "normal|adjacent",
      y
    ),
    "Normal",

    default =
      NA_character_
  )
}

# ============================================================
# 2. GSE221561
# ============================================================

cat("\n===== 2. GSE221561 EPITHELIAL INVENTORY =====\n\n")

x221 <- fread(f221)

req221 <- c(
  "patient_id",
  "sample_id",
  "tissue_std",
  "treatment_group",
  "cell_class_std"
)

miss221 <- setdiff(
  req221,
  names(x221)
)

if (length(miss221) > 0) {
  stop(
    "GSE221561 missing: ",
    paste(miss221, collapse = ", ")
  )
}

# Complete sample inventory first
inv221 <- unique(
  x221[
    ,
    .(
      dataset = "GSE221561",
      patient_id =
        as.character(patient_id),
      sample_id =
        as.character(sample_id),
      tissue =
        standardize_tissue(tissue_std),
      treatment_group =
        as.character(treatment_group)
    )
  ]
)

epi221 <- x221[
  cell_class_std ==
    "Epithelial_unclassified",
  .(
    epithelial_cells = .N
  ),
  by = .(
    patient_id =
      as.character(patient_id),
    sample_id =
      as.character(sample_id)
  )
]

g221 <- merge(
  inv221,
  epi221,
  by = c(
    "patient_id",
    "sample_id"
  ),
  all.x = TRUE
)

g221[
  is.na(epithelial_cells),
  epithelial_cells := 0L
]

cat(
  "Samples:",
  nrow(g221),
  "\n"
)

cat(
  "Epithelial cells:",
  sum(g221$epithelial_cells),
  "\n"
)

# ============================================================
# 3. GSE197677
# ============================================================

cat("\n===== 3. GSE197677 EPITHELIAL INVENTORY =====\n\n")

x197 <- fread(f197)

req197 <- c(
  "patient_id",
  "sample_id",
  "tissue",
  "treatment_group",
  "broad_celltype",
  "n_cells"
)

miss197 <- setdiff(
  req197,
  names(x197)
)

if (length(miss197) > 0) {
  stop(
    "GSE197677 missing: ",
    paste(miss197, collapse = ", ")
  )
}

inv197 <- unique(
  x197[
    ,
    .(
      dataset = "GSE197677",
      patient_id =
        as.character(patient_id),
      sample_id =
        as.character(sample_id),
      tissue =
        standardize_tissue(tissue),
      treatment_group =
        as.character(treatment_group)
    )
  ]
)

epi197 <- x197[
  broad_celltype == "Epithelial",
  .(
    epithelial_cells =
      sum(
        n_cells,
        na.rm = TRUE
      )
  ),
  by = .(
    patient_id =
      as.character(patient_id),
    sample_id =
      as.character(sample_id)
  )
]

g197 <- merge(
  inv197,
  epi197,
  by = c(
    "patient_id",
    "sample_id"
  ),
  all.x = TRUE
)

g197[
  is.na(epithelial_cells),
  epithelial_cells := 0L
]

cat(
  "Samples:",
  nrow(g197),
  "\n"
)

cat(
  "Epithelial cells:",
  sum(g197$epithelial_cells),
  "\n"
)

# ============================================================
# 4. GSE160269
# ============================================================

cat("\n===== 4. GSE160269 EPITHELIAL INVENTORY =====\n\n")

x160 <- fread(f160)

req160 <- c(
  "patient_id",
  "specimen_id",
  "tissue",
  "cell_type_author_source",
  "n_cells"
)

miss160 <- setdiff(
  req160,
  names(x160)
)

if (length(miss160) > 0) {
  stop(
    "GSE160269 missing: ",
    paste(miss160, collapse = ", ")
  )
}

inv160 <- unique(
  x160[
    ,
    .(
      dataset = "GSE160269",
      patient_id =
        as.character(patient_id),
      sample_id =
        as.character(specimen_id),
      tissue =
        standardize_tissue(tissue),

      treatment_group =
        fifelse(
          standardize_tissue(tissue) ==
            "Tumor",
          "Untreated_baseline",
          "Adjacent_normal"
        )
    )
  ]
)

epi160 <- x160[
  cell_type_author_source ==
    "Epithelia",
  .(
    epithelial_cells =
      sum(
        n_cells,
        na.rm = TRUE
      )
  ),
  by = .(
    patient_id =
      as.character(patient_id),
    sample_id =
      as.character(specimen_id)
  )
]

g160 <- merge(
  inv160,
  epi160,
  by = c(
    "patient_id",
    "sample_id"
  ),
  all.x = TRUE
)

g160[
  is.na(epithelial_cells),
  epithelial_cells := 0L
]

cat(
  "Samples:",
  nrow(g160),
  "\n"
)

cat(
  "Epithelial cells:",
  sum(g160$epithelial_cells),
  "\n"
)

# ============================================================
# 5. COMBINE
# ============================================================

epi <- rbindlist(
  list(
    g221,
    g197,
    g160
  ),
  use.names = TRUE,
  fill = TRUE
)

epi[
  ,
  epithelial_eligible_20 :=
    epithelial_cells >= 20
]

epi[
  ,
  epithelial_eligible_30 :=
    epithelial_cells >= 30
]

epi[
  ,
  epithelial_eligible_100 :=
    epithelial_cells >= 100
]

setorder(
  epi,
  dataset,
  tissue,
  sample_id
)

cat("\n============================================================\n")
cat("5. ALL SAMPLE EPITHELIAL INVENTORY\n")
cat("============================================================\n\n")

print(
  epi,
  nrows = 200
)

fwrite(
  epi,
  file.path(
    outdir,
    "STEP17A_all_sample_epithelial_inventory.csv"
  ),
  bom = TRUE
)

# ============================================================
# 6. DATASET-LEVEL EPITHELIAL CAPACITY
# ============================================================

capacity <- epi[
  ,
  .(
    samples_total = .N,

    samples_epi_ge20 =
      sum(
        as.integer(epithelial_eligible_20)
      ),

    samples_epi_ge30 =
      sum(
        as.integer(epithelial_eligible_30)
      ),

    samples_epi_ge100 =
      sum(
        as.integer(epithelial_eligible_100)
      ),

    total_epithelial_cells =
      sum(
        epithelial_cells
      ),

    median_epithelial_cells =
      as.numeric(
        median(
          epithelial_cells
        )
      ),

    min_epithelial_cells =
      as.numeric(
        min(
          epithelial_cells
        )
      ),

    max_epithelial_cells =
      as.numeric(
        max(
          epithelial_cells
        )
      )
  ),
  by = .(
    dataset,
    tissue,
    treatment_group
  )
]

setorder(
  capacity,
  dataset,
  tissue,
  treatment_group
)

cat("\n============================================================\n")
cat("6. EPITHELIAL CAPACITY BY DATASET / GROUP\n")
cat("============================================================\n\n")

print(
  capacity,
  nrows = 100
)

fwrite(
  capacity,
  file.path(
    outdir,
    "STEP17A_epithelial_capacity_by_group.csv"
  ),
  bom = TRUE
)

# ============================================================
# 7. NORMAL EPITHELIAL REFERENCE CAPACITY
# ============================================================

normal_reference <- epi[
  tissue == "Normal",
  .(
    normal_samples_total = .N,

    normal_samples_ge20 =
      sum(
        as.integer(
          epithelial_cells >= 20
        )
      ),

    normal_samples_ge30 =
      sum(
        as.integer(
          epithelial_cells >= 30
        )
      ),

    normal_samples_ge100 =
      sum(
        as.integer(
          epithelial_cells >= 100
        )
      ),

    total_normal_epithelial_cells =
      sum(
        epithelial_cells
      ),

    median_normal_epithelial_cells =
      as.numeric(
        median(
          epithelial_cells
        )
      ),

    min_normal_epithelial_cells =
      as.numeric(
        min(
          epithelial_cells
        )
      ),

    max_normal_epithelial_cells =
      as.numeric(
        max(
          epithelial_cells
        )
      )
  ),
  by = dataset
]

# Technical feasibility flag only.
# This is NOT a biological validation threshold.
normal_reference[
  ,
  technical_reference_status :=
    fcase(

      normal_samples_ge30 >= 2 &
      total_normal_epithelial_cells >= 100,
      "REFERENCE_TECHNICALLY_FEASIBLE",

      normal_samples_ge20 >= 1,
      "LIMITED_REFERENCE_REVIEW_REQUIRED",

      default =
        "REFERENCE_INSUFFICIENT"
    )
]

cat("\n============================================================\n")
cat("7. NORMAL EPITHELIAL REFERENCE CAPACITY\n")
cat("============================================================\n\n")

print(normal_reference)

cat(
  "\nNOTE:\n",
  "The status above is a technical workflow flag only.\n",
  "It is NOT proof that a normal reference is biologically perfect.\n",
  sep = ""
)

fwrite(
  normal_reference,
  file.path(
    outdir,
    "STEP17A_normal_epithelial_reference_capacity.csv"
  ),
  bom = TRUE
)

# ============================================================
# 8. PAIRED TUMOR-NORMAL SUBJECTS
# ============================================================

pair_audit <- epi[
  !is.na(patient_id) &
    patient_id != "",
  .(
    has_tumor =
      any(tissue == "Tumor"),

    has_normal =
      any(tissue == "Normal"),

    tumor_samples =
      paste(
        sort(
          unique(
            sample_id[
              tissue == "Tumor"
            ]
          )
        ),
        collapse = "; "
      ),

    normal_samples =
      paste(
        sort(
          unique(
            sample_id[
              tissue == "Normal"
            ]
          )
        ),
        collapse = "; "
      ),

    tumor_epi_cells =
      sum(
        epithelial_cells[
          tissue == "Tumor"
        ]
      ),

    normal_epi_cells =
      sum(
        epithelial_cells[
          tissue == "Normal"
        ]
      )
  ),
  by = .(
    dataset,
    patient_id
  )
]

pair_audit[
  ,
  paired_tumor_normal :=
    has_tumor &
    has_normal
]

pair_audit[
  ,
  paired_reference_ge30 :=
    paired_tumor_normal &
    tumor_epi_cells >= 30 &
    normal_epi_cells >= 30
]

paired_only <- pair_audit[
  paired_tumor_normal == TRUE
]

cat("\n============================================================\n")
cat("8. PAIRED TUMOR / NORMAL EPITHELIAL SUBJECTS\n")
cat("============================================================\n\n")

if (nrow(paired_only) == 0) {

  cat(
    "No paired subjects detected.\n"
  )

} else {

  print(
    paired_only,
    nrows = 100
  )
}

fwrite(
  pair_audit,
  file.path(
    outdir,
    "STEP17A_tumor_normal_pairing_epithelial_audit.csv"
  ),
  bom = TRUE
)

# ============================================================
# 9. TUMOR EPITHELIAL ELIGIBILITY
# ============================================================

tumor_epi <- epi[
  tissue == "Tumor"
]

tumor_summary <- tumor_epi[
  ,
  .(
    tumor_samples_total =
      .N,

    tumor_samples_ge20 =
      sum(
        as.integer(
          epithelial_cells >= 20
        )
      ),

    tumor_samples_ge30 =
      sum(
        as.integer(
          epithelial_cells >= 30
        )
      ),

    tumor_samples_ge100 =
      sum(
        as.integer(
          epithelial_cells >= 100
        )
      ),

    tumor_epithelial_cells =
      sum(
        epithelial_cells
      ),

    median_tumor_epithelial =
      as.numeric(
        median(
          epithelial_cells
        )
      )
  ),
  by = .(
    dataset,
    treatment_group
  )
]

cat("\n============================================================\n")
cat("9. TUMOR EPITHELIAL ANALYSIS CAPACITY\n")
cat("============================================================\n\n")

print(tumor_summary)

fwrite(
  tumor_summary,
  file.path(
    outdir,
    "STEP17A_tumor_epithelial_analysis_capacity.csv"
  ),
  bom = TRUE
)

# ============================================================
# 10. MALIGNANCY-INFERENCE STRATEGY PLAN
# ============================================================

strategy <- data.table(

  dataset = c(
    "GSE221561",
    "GSE197677",
    "GSE160269"
  ),

  epithelial_source = c(
    "Author-retained broad epithelial cells",
    "Reviewed broad epithelial clusters",
    "Author epithelial UMI matrix"
  ),

  normal_reference_type = c(
    "Adjacent-normal epithelial cells",
    "Normal esophageal epithelial cells",
    "Adjacent-normal epithelial cells"
  ),

  proposed_next_method = c(
    "CNV-based malignant inference within epithelial compartment",
    "CNV-based malignant inference within epithelial compartment",
    "CNV-based malignant inference within epithelial compartment"
  ),

  important_rule = c(
    paste0(
      "Do not equate tumor-tissue epithelial cells ",
      "with malignant cells automatically"
    ),

    paste0(
      "Do not equate tumor-tissue epithelial cells ",
      "with malignant cells automatically"
    ),

    paste0(
      "Do not equate tumor-tissue epithelial cells ",
      "with malignant cells automatically"
    )
  )
)

strategy <- merge(
  strategy,
  normal_reference[
    ,
    .(
      dataset,
      normal_samples_ge30,
      total_normal_epithelial_cells,
      technical_reference_status
    )
  ],
  by = "dataset",
  all.x = TRUE
)

cat("\n============================================================\n")
cat("10. MALIGNANT EPITHELIAL STRATEGY PLAN\n")
cat("============================================================\n\n")

print(strategy)

fwrite(
  strategy,
  file.path(
    outdir,
    "STEP17A_malignant_epithelial_strategy_plan.csv"
  ),
  bom = TRUE
)

# ============================================================
# 11. PAIRING SUMMARY
# ============================================================

pair_summary <- pair_audit[
  ,
  .(
    subjects =
      uniqueN(patient_id),

    paired_subjects =
      sum(
        as.integer(paired_tumor_normal)
      ),

    paired_with_epi_ge30 =
      sum(
        as.integer(paired_reference_ge30)
      )
  ),
  by = dataset
]

cat("\n============================================================\n")
cat("11. TUMOR-NORMAL PAIRING SUMMARY\n")
cat("============================================================\n\n")

print(pair_summary)

fwrite(
  pair_summary,
  file.path(
    outdir,
    "STEP17A_pairing_summary.csv"
  ),
  bom = TRUE
)

# ============================================================
# FINAL
# ============================================================

cat("\n============================================================\n")
cat("FINAL INTERPRETATION GATE\n")
cat("============================================================\n\n")

cat(
  "This step did NOT call any epithelial cell malignant.\n\n"
)

cat(
  "Broad epithelial = epithelial compartment only.\n"
)

cat(
  "Tumor tissue origin = NOT sufficient evidence of malignancy.\n"
)

cat(
  "Normal epithelial references will be evaluated for CNV-based ",
  "malignancy inference in the next step.\n"
)

cat("\n============================================================\n")
cat("STEP 17A FINISHED ✓\n")
cat("============================================================\n\n")

cat("No large Seurat objects loaded ✓\n")
cat("No malignant labels assigned ✓\n")
cat("No CNV analysis run ✓\n")
cat("No cells filtered ✓\n")
cat("No pseudobulk DE run ✓\n")
cat("No cross-dataset integration ✓\n")

cat("\nMain outputs:\n")

cat(
  "  06_tables/cross_dataset/malignant_epithelial/",
  "STEP17A_normal_epithelial_reference_capacity.csv\n",
  sep = ""
)

cat(
  "  06_tables/cross_dataset/malignant_epithelial/",
  "STEP17A_tumor_epithelial_analysis_capacity.csv\n",
  sep = ""
)

cat(
  "  06_tables/cross_dataset/malignant_epithelial/",
  "STEP17A_tumor_normal_pairing_epithelial_audit.csv\n",
  sep = ""
)

cat(
  "  06_tables/cross_dataset/malignant_epithelial/",
  "STEP17A_malignant_epithelial_strategy_plan.csv\n",
  sep = ""
)

cat("\nSTOP HERE.\n")
