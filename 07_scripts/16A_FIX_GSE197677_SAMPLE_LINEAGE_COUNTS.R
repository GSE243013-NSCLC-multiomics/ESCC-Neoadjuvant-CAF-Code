# ============================================================
# ESCC Neoadjuvant Project
# STEP 16A
#
# Rebuild GSE197677 sample x broad-cell-type counts
# and correct cross-dataset lineage eligibility.
#
# IMPORTANT:
# - uses lightweight CSV files only
# - does NOT load large Seurat RDS
# - NO pseudobulk DE
# - NO filtering
# - NO normalization
# - NO integration
# ============================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
})

cat("\n============================================================\n")
cat("STEP 16A: FIX GSE197677 SAMPLE-LINEAGE COUNTS\n")
cat("============================================================\n\n")

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------

meta_file <- paste0(
  "06_tables/GSE197677/metadata/",
  "GSE197677_standardized_cell_metadata.csv.gz"
)

mapping_file <- paste0(
  "06_tables/GSE197677/celltype/",
  "GSE197677_FINAL_cluster_to_broad_celltype_mapping.csv"
)

old_eligibility_file <- paste0(
  "06_tables/cross_dataset/",
  "STEP16_primary_lineage_sample_eligibility.csv"
)

outdir <- "06_tables/cross_dataset"

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)

for (f in c(
  meta_file,
  mapping_file,
  old_eligibility_file
)) {

  if (!file.exists(f)) {
    stop("Missing required file: ", f)
  }
}

# ============================================================
# 1. LOAD LIGHTWEIGHT GSE197677 METADATA
# ============================================================

cat("===== 1. LOAD GSE197677 CELL METADATA =====\n\n")

md <- fread(meta_file)

cat(
  "Cell metadata rows: ",
  nrow(md),
  "\n",
  sep = ""
)

if (nrow(md) != 71116) {
  stop(
    "Expected 71,116 GSE197677 cells; found ",
    nrow(md)
  )
}

required_md <- c(
  "sample",
  "subject_id_internal",
  "tissue_std",
  "treatment_group",
  "seurat_clusters"
)

missing_md <- setdiff(
  required_md,
  names(md)
)

if (length(missing_md) > 0) {

  cat("Available metadata fields:\n")
  print(names(md))

  stop(
    "Missing metadata fields: ",
    paste(missing_md, collapse = ", ")
  )
}

# ============================================================
# 2. LOAD FINAL CLUSTER -> BROAD CELL TYPE MAPPING
# ============================================================

cat("\n===== 2. LOAD FINAL CLUSTER MAPPING =====\n\n")

mapping <- fread(mapping_file)

required_map <- c(
  "cluster",
  "broad_celltype_final"
)

missing_map <- setdiff(
  required_map,
  names(mapping)
)

if (length(missing_map) > 0) {
  stop(
    "Mapping fields missing: ",
    paste(missing_map, collapse = ", ")
  )
}

mapping[
  ,
  cluster_key :=
    as.character(cluster)
]

md[
  ,
  cluster_key :=
    as.character(seurat_clusters)
]

cat(
  "Clusters in cell metadata: ",
  uniqueN(md$cluster_key),
  "\n",
  sep = ""
)

cat(
  "Clusters in final mapping: ",
  uniqueN(mapping$cluster_key),
  "\n",
  sep = ""
)

# ============================================================
# 3. JOIN FINAL BROAD CELL TYPE TO EVERY CELL
# ============================================================

cat("\n===== 3. MAP BROAD CELL TYPE TO 71,116 CELLS =====\n\n")

md2 <- merge(
  md,
  mapping[
    ,
    .(
      cluster_key,
      broad_celltype_final
    )
  ],
  by = "cluster_key",
  all.x = TRUE,
  sort = FALSE
)

cat(
  "Rows after mapping: ",
  nrow(md2),
  "\n",
  sep = ""
)

unmapped <- sum(
  is.na(md2$broad_celltype_final)
)

cat(
  "Unmapped cells: ",
  unmapped,
  "\n",
  sep = ""
)

if (
  nrow(md2) != 71116 ||
  unmapped != 0
) {
  stop(
    "Broad-cell-type mapping is incomplete."
  )
}

# ============================================================
# 4. VERIFY FINAL BROAD CELL TOTALS
# ============================================================

cat("\n===== 4. GSE197677 BROAD CELL TOTALS =====\n\n")

broad_totals <- md2[
  ,
  .(
    n_cells = .N
  ),
  by = broad_celltype_final
]

setorder(
  broad_totals,
  -n_cells
)

print(broad_totals)

cat(
  "\nTotal: ",
  sum(broad_totals$n_cells),
  "\n",
  sep = ""
)

if (
  sum(broad_totals$n_cells) != 71116
) {
  stop(
    "Broad-cell total does not equal 71,116."
  )
}

# ============================================================
# 5. REBUILD SAMPLE x BROAD CELL TYPE COUNTS
# ============================================================

cat("\n===== 5. REBUILD SAMPLE x BROAD CELL TYPE COUNTS =====\n\n")

sample_counts_all <- md2[
  ,
  .(
    n_cells = .N
  ),
  by = .(
    patient_id =
      subject_id_internal,

    sample_id =
      sample,

    tissue =
      tissue_std,

    treatment_group,

    broad_celltype =
      broad_celltype_final
  )
]

setorder(
  sample_counts_all,
  sample_id,
  broad_celltype
)

print(
  sample_counts_all,
  nrows = 300
)

fwrite(
  sample_counts_all,
  paste0(
    outdir,
    "/STEP16A_GSE197677_sample_by_broad_celltype_REBUILT.csv"
  ),
  bom = TRUE
)

# ============================================================
# 6. CHECK SAMPLE COHORT STRUCTURE
# ============================================================

cat("\n===== 6. GSE197677 SAMPLE COHORT =====\n\n")

sample_inventory_197 <- unique(
  sample_counts_all[
    ,
    .(
      patient_id,
      sample_id,
      tissue,
      treatment_group
    )
  ]
)

sample_group_check <- sample_inventory_197[
  ,
  .(
    n_samples =
      uniqueN(sample_id),

    n_subjects =
      uniqueN(
        patient_id,
        na.rm = TRUE
      ),

    samples =
      paste(
        sort(
          unique(sample_id)
        ),
        collapse = ", "
      )
  ),
  by = .(
    tissue,
    treatment_group
  )
]

print(sample_group_check)

# ============================================================
# 7. BUILD COMPLETE GRID FOR THREE PRIMARY LINEAGES
# ============================================================

primary_lineages <- c(
  "Epithelial",
  "Myeloid",
  "Fibroblast_CAF"
)

grid197 <- sample_inventory_197[
  ,
  .(
    broad_celltype =
      primary_lineages
  ),
  by = .(
    patient_id,
    sample_id,
    tissue,
    treatment_group
  )
]

observed197 <- sample_counts_all[
  broad_celltype %in%
    primary_lineages
]

elig197 <- merge(
  grid197,
  observed197,
  by = c(
    "patient_id",
    "sample_id",
    "tissue",
    "treatment_group",
    "broad_celltype"
  ),
  all.x = TRUE
)

elig197[
  is.na(n_cells),
  n_cells := 0L
]

elig197[
  ,
  dataset := "GSE197677"
]

elig197[
  ,
  eligible_20 :=
    n_cells >= 20
]

elig197[
  ,
  eligible_30 :=
    n_cells >= 30
]

elig197[
  ,
  analysis_status :=
    fcase(

      n_cells >= 30,
      "ELIGIBLE_30",

      n_cells >= 20,
      "ELIGIBLE_20_ONLY",

      default =
        "MISSING_NOT_ZERO"
    )
]

elig197[
  ,
  important_note :=
    fifelse(
      analysis_status ==
        "MISSING_NOT_ZERO",

      paste0(
        "Do not encode as expression zero; ",
        "exclude sample-lineage from pseudobulk"
      ),

      "Eligible cell inventory"
    )
]

setcolorder(
  elig197,
  c(
    "dataset",
    "patient_id",
    "sample_id",
    "tissue",
    "treatment_group",
    "broad_celltype",
    "n_cells",
    "eligible_20",
    "eligible_30",
    "analysis_status",
    "important_note"
  )
)

# ============================================================
# 8. GSE197677 CORRECTED ELIGIBILITY
# ============================================================

cat("\n============================================================\n")
cat("7. GSE197677 CORRECTED TUMOR ELIGIBILITY\n")
cat("============================================================\n\n")

summary197 <- elig197[
  tissue == "Tumor",
  .(
    samples_total =
      uniqueN(sample_id),

    samples_eligible_20 =
      uniqueN(
        sample_id[
          eligible_20
        ]
      ),

    samples_eligible_30 =
      uniqueN(
        sample_id[
          eligible_30
        ]
      ),

    median_cells =
      median(n_cells),

    min_cells =
      min(n_cells),

    max_cells =
      max(n_cells)
  ),
  by = .(
    treatment_group,
    broad_celltype
  )
]

setorder(
  summary197,
  broad_celltype,
  treatment_group
)

print(summary197)

fwrite(
  summary197,
  paste0(
    outdir,
    "/STEP16A_GSE197677_corrected_tumor_eligibility.csv"
  ),
  bom = TRUE
)

# ============================================================
# 9. PATCH ORIGINAL STEP16 ELIGIBILITY
# ============================================================

cat("\n===== 8. PATCH CROSS-DATASET ELIGIBILITY =====\n\n")

old <- fread(
  old_eligibility_file
)

old_without_197 <- old[
  dataset != "GSE197677"
]

corrected <- rbindlist(
  list(
    old_without_197,
    elig197
  ),
  use.names = TRUE,
  fill = TRUE
)

setorder(
  corrected,
  dataset,
  sample_id,
  broad_celltype
)

fwrite(
  corrected,
  paste0(
    outdir,
    "/STEP16A_primary_lineage_sample_eligibility_CORRECTED.csv"
  ),
  bom = TRUE
)

# ============================================================
# 10. TREATMENT-GROUP-SPECIFIC CAPACITY
# ============================================================

cat("\n============================================================\n")
cat("9. CORRECTED TUMOR CAPACITY BY TREATMENT GROUP\n")
cat("============================================================\n\n")

tumor_capacity <- corrected[
  tissue == "Tumor",
  .(
    samples_total =
      uniqueN(sample_id),

    patients_total =
      uniqueN(
        patient_id,
        na.rm = TRUE
      ),

    samples_eligible_20 =
      sum(
        as.integer(eligible_20)
      ),

    samples_eligible_30 =
      sum(
        as.integer(eligible_30)
      ),

    patients_eligible_30 =
      uniqueN(
        patient_id[
          eligible_30
        ],
        na.rm = TRUE
      ),

    median_cells =
      as.numeric(
        median(n_cells)
      ),

    min_cells =
      as.numeric(
        min(n_cells)
      ),

    max_cells =
      as.numeric(
        max(n_cells)
      )
  ),
  by = .(
    dataset,
    treatment_group,
    broad_celltype
  )
]

setorder(
  tumor_capacity,
  broad_celltype,
  dataset,
  treatment_group
)

print(
  tumor_capacity,
  nrows = 200
)

fwrite(
  tumor_capacity,
  paste0(
    outdir,
    "/STEP16A_tumor_capacity_by_treatment_CORRECTED.csv"
  ),
  bom = TRUE
)

# ============================================================
# 11-13. DATASET-SPECIFIC CHECKS
# ============================================================

cat("\n============================================================\n")
cat("10. GSE221561 TUMOR GROUP CHECK\n")
cat("============================================================\n\n")

print(
  tumor_capacity[
    dataset == "GSE221561"
  ]
)

cat("\n============================================================\n")
cat("11. GSE197677 TUMOR GROUP CHECK\n")
cat("============================================================\n\n")

print(
  tumor_capacity[
    dataset == "GSE197677"
  ]
)

cat("\n============================================================\n")
cat("12. GSE160269 TUMOR GROUP CHECK\n")
cat("============================================================\n\n")

print(
  tumor_capacity[
    dataset == "GSE160269"
  ]
)

# ============================================================
# 14. WIDE HUMAN-READABLE GSE197677 TABLE
# ============================================================

wide197 <- dcast(
  elig197,
  patient_id +
    sample_id +
    tissue +
    treatment_group ~
    broad_celltype,
  value.var = "n_cells",
  fill = 0
)

cat("\n============================================================\n")
cat("13. GSE197677 SAMPLE x PRIMARY LINEAGE COUNTS\n")
cat("============================================================\n\n")

print(
  wide197,
  nrows = 100
)

fwrite(
  wide197,
  paste0(
    outdir,
    "/STEP16A_GSE197677_primary_lineage_counts_WIDE.csv"
  ),
  bom = TRUE
)

# ============================================================
# 15. FINAL VALIDATION
# ============================================================

cat("\n============================================================\n")
cat("FINAL VALIDATION\n")
cat("============================================================\n\n")

cat(
  "GSE197677 cells reconstructed : ",
  sum(broad_totals$n_cells),
  "\n",
  sep = ""
)

cat(
  "GSE197677 samples             : ",
  uniqueN(sample_inventory_197$sample_id),
  "\n",
  sep = ""
)

cat(
  "Tumor NACT samples            : ",
  uniqueN(
    sample_inventory_197[
      tissue == "Tumor" &
      treatment_group ==
        "Neoadjuvant_chemotherapy",
      sample_id
    ]
  ),
  "\n",
  sep = ""
)

cat(
  "Tumor nNACT samples           : ",
  uniqueN(
    sample_inventory_197[
      tissue == "Tumor" &
      treatment_group ==
        "No_neoadjuvant_chemotherapy",
      sample_id
    ]
  ),
  "\n",
  sep = ""
)

cat(
  "Normal samples                : ",
  uniqueN(
    sample_inventory_197[
      tissue == "Normal",
      sample_id
    ]
  ),
  "\n",
  sep = ""
)

ok <- (
  sum(broad_totals$n_cells) == 71116 &&
  uniqueN(sample_inventory_197$sample_id) == 14 &&
  uniqueN(
    sample_inventory_197[
      tissue == "Tumor" &
      treatment_group ==
        "Neoadjuvant_chemotherapy",
      sample_id
    ]
  ) == 6 &&
  uniqueN(
    sample_inventory_197[
      tissue == "Tumor" &
      treatment_group ==
        "No_neoadjuvant_chemotherapy",
      sample_id
    ]
  ) == 4 &&
  uniqueN(
    sample_inventory_197[
      tissue == "Normal",
      sample_id
    ]
  ) == 4
)

cat("\n============================================================\n")

if (ok) {

  cat("STEP 16A FINISHED ✓\n")
  cat("GSE197677 SAMPLE-LINEAGE COUNTS REBUILT ✓\n")
  cat("CROSS-DATASET ELIGIBILITY CORRECTED ✓\n")

} else {

  cat("STEP 16A FINISHED WITH REVIEW FLAG ⚠\n")
  cat("DO NOT START PSEUDOBULK YET.\n")
}

cat("============================================================\n\n")

cat("No large Seurat object loaded ✓\n")
cat("No cells filtered ✓\n")
cat("No pseudobulk expression generated ✓\n")
cat("No differential expression run ✓\n")
cat("No cross-dataset integration ✓\n")

cat("\nSTOP HERE.\n")
