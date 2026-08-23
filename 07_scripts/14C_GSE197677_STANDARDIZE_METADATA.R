# ============================================================
# ESCC Neoadjuvant Project
# STEP 14C
# Standardize GSE197677 sample / tissue / treatment metadata
# and verify RNA counts availability
#
# NO FILTERING
# NO NORMALIZATION
# NO PCA / UMAP
# NO CROSS-DATASET INTEGRATION
# ============================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
})

cat("\n============================================================\n")
cat("STEP 14C: GSE197677 METADATA STANDARDIZATION\n")
cat("============================================================\n\n")

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------

object_dir <- "03_objects/GSE197677"

out_dir <- "06_tables/GSE197677/metadata"

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# 1. Locate integrated RDS
# ------------------------------------------------------------

files <- list.files(
  object_dir,
  pattern = "\\.rds$",
  full.names = TRUE,
  ignore.case = TRUE
)

hits <- files[
  grepl(
    "integrated",
    basename(files),
    ignore.case = TRUE
  )
]

if (length(hits) == 0) {
  stop("Integrated GSE197677 RDS not found.")
}

object_file <- hits[1]

cat("===== 1. LOADING OBJECT =====\n\n")

cat(
  "Object: ",
  object_file,
  "\n",
  sep = ""
)

obj <- readRDS(
  object_file
)

md <- as.data.table(
  obj[[]],
  keep.rownames = "cell_id"
)

cat(
  "Cells: ",
  nrow(md),
  "\n",
  sep = ""
)

cat(
  "Metadata fields: ",
  ncol(md),
  "\n",
  sep = ""
)

# ------------------------------------------------------------
# 2. Verify required fields
# ------------------------------------------------------------

required <- c(
  "sample",
  "CellType",
  "NAC",
  "subtype"
)

missing <- setdiff(
  required,
  names(md)
)

if (length(missing) > 0) {
  stop(
    "Missing required fields: ",
    paste(missing, collapse = ", ")
  )
}

# ------------------------------------------------------------
# 3. RNA counts audit
# ------------------------------------------------------------

cat("\n===== 2. RNA COUNTS AUDIT =====\n\n")

assay_names <- tryCatch(
  Assays(obj),
  error = function(e) character(0)
)

cat("Available assays: ", paste(assay_names, collapse = ", "), "\n", sep = "")

if (!"RNA" %in% assay_names) {

  stop(
    "RNA assay is absent. ",
    "Patient-level raw-count pseudobulk cannot proceed."
  )
}

counts_available <- FALSE

# Try to get counts using GetAssayData with v4 syntax
counts_available <- tryCatch(
  {
    x <- GetAssayData(
      obj,
      assay = "RNA",
      slot = "counts"
    )
    nrow(x) > 0 && ncol(x) > 0
  },
  error = function(e) {
    cat("GetAssayData with slot='counts' failed: ", conditionMessage(e), "\n")
    FALSE
  }
)

# If that failed, try v5 syntax
if (!counts_available) {
  counts_available <- tryCatch(
    {
      x <- GetAssayData(
        obj,
        assay = "RNA",
        layer = "counts"
      )
      nrow(x) > 0 && ncol(x) > 0
    },
    error = function(e) {
      cat("GetAssayData with layer='counts' failed: ", conditionMessage(e), "\n")
      FALSE
    }
  )
}

cat(
  "RNA raw counts available: ",
  counts_available,
  "\n",
  sep = ""
)

if (!counts_available) {

  cat(
    "\nWARNING:\n",
    "RNA counts could not be verified.\n",
    "Do NOT use SCT/integrated values for pseudobulk DE.\n",
    sep = ""
  )
}

# ------------------------------------------------------------
# 4. Standardize tissue
# ------------------------------------------------------------

md[
  ,
  tissue_std :=
    fcase(

      tolower(
        trimws(
          as.character(CellType)
        )
      ) == "scc",
      "Tumor",

      tolower(
        trimws(
          as.character(CellType)
        )
      ) == "normal",
      "Normal",

      default = NA_character_
    )
]

# ------------------------------------------------------------
# 5. Standardize treatment
# ------------------------------------------------------------

md[
  ,
  treatment_raw :=
    as.character(NAC)
]

md[
  ,
  treatment_group :=
    fcase(

      tissue_std == "Normal",
      "Normal",

      toupper(
        trimws(treatment_raw)
      ) == "NACT",
      "Neoadjuvant_chemotherapy",

      toupper(
        trimws(treatment_raw)
      ) == "NNACT",
      "No_neoadjuvant_chemotherapy",

      default = NA_character_
    )
]

# ------------------------------------------------------------
# 6. Internal subject key candidate
#
# ESC15 / ESN15 -> subject code 15
#
# This is an INTERNAL pairing key only.
# We do NOT claim it maps to a GEO accession.
# ------------------------------------------------------------

md[
  ,
  subject_code_internal :=
    sub(
      "^ES[CN]",
      "",
      as.character(sample),
      ignore.case = TRUE
    )
]

valid_subject_code <- grepl(
  "^[0-9]+$",
  md$subject_code_internal
)

md[
  !valid_subject_code,
  subject_code_internal := NA_character_
]

md[
  ,
  subject_id_internal :=
    fifelse(
      !is.na(subject_code_internal),
      paste0(
        "SUBJ",
        subject_code_internal
      ),
      NA_character_
    )
]

# ------------------------------------------------------------
# 7. Cell-level completeness
# ------------------------------------------------------------

md[
  ,
  metadata_complete :=
    !is.na(sample) &
    !is.na(tissue_std) &
    !is.na(treatment_group)
]

cat("\n===== 3. STANDARDIZATION COMPLETENESS =====\n\n")

cat(
  "Complete cells : ",
  sum(md$metadata_complete),
  "\n",
  sep = ""
)

cat(
  "Incomplete     : ",
  sum(!md$metadata_complete),
  "\n",
  sep = ""
)

# ------------------------------------------------------------
# 8. Exact sample-level table
# ------------------------------------------------------------

collapse_unique <- function(x) {

  z <- unique(
    as.character(x)
  )

  z <- z[
    !is.na(z) &
    trimws(z) != ""
  ]

  if (length(z) == 0) {
    return(NA_character_)
  }

  paste(
    sort(z),
    collapse = "; "
  )
}

sample_table <- md[
  ,
  .(
    n_cells = .N,

    tissue =
      collapse_unique(
        tissue_std
      ),

    treatment =
      collapse_unique(
        treatment_group
      ),

    treatment_raw =
      collapse_unique(
        treatment_raw
      ),

    subtype =
      collapse_unique(
        subtype
      ),

    subject_id_internal =
      collapse_unique(
        subject_id_internal
      )
  ),
  by = sample
]

setorder(
  sample_table,
  sample
)

cat("\n===== 4. STANDARDIZED SAMPLE TABLE =====\n\n")

print(
  sample_table
)

# ------------------------------------------------------------
# 9. Tumor-only treatment comparison
# ------------------------------------------------------------

tumor_samples <- sample_table[
  tissue == "Tumor"
]

cat("\n===== 5. TUMOR TREATMENT COHORT =====\n\n")

print(
  tumor_samples
)

treatment_summary <- tumor_samples[
  ,
  .(
    n_samples =
      uniqueN(sample),

    samples =
      paste(
        sort(
          unique(sample)
        ),
        collapse = ", "
      )
  ),
  by = treatment
]

cat("\nTumor comparison groups:\n\n")

print(
  treatment_summary
)

# ------------------------------------------------------------
# 10. Internal tumor-normal pairing audit
# ------------------------------------------------------------

pairing <- sample_table[
  !is.na(subject_id_internal),
  .(
    n_samples =
      uniqueN(sample),

    samples =
      paste(
        sort(
          unique(sample)
        ),
        collapse = "; "
      ),

    tissues =
      paste(
        sort(
          unique(tissue)
        ),
        collapse = "; "
      ),

    treatments =
      paste(
        sort(
          unique(treatment)
        ),
        collapse = "; "
      )
  ),
  by = subject_id_internal
]

pairing[
  ,
  paired_tumor_normal :=
    grepl(
      "Tumor",
      tissues,
      fixed = TRUE
    ) &
    grepl(
      "Normal",
      tissues,
      fixed = TRUE
    )
]

setorder(
  pairing,
  subject_id_internal
)

cat("\n===== 6. INTERNAL TUMOR / NORMAL PAIRING AUDIT =====\n\n")

print(
  pairing
)

cat(
  "\nPotential paired tumor/normal internal IDs: ",
  sum(pairing$paired_tumor_normal),
  "\n",
  sep = ""
)

cat(
  "\nIMPORTANT:\n",
  "These are INTERNAL ESC/ESN suffix-based pairing keys.\n",
  "They are not claimed to be GEO accession mappings.\n",
  sep = ""
)

# ------------------------------------------------------------
# 11. Author-QC boundary audit
# ------------------------------------------------------------

cat("\n===== 7. AUTHOR QC BOUNDARY AUDIT =====\n\n")

qc_rows <- list()

for (
  field in intersect(
    c(
      "nCount_RNA",
      "nFeature_RNA",
      "percent.mt"
    ),
    names(md)
  )
) {

  x <- suppressWarnings(
    as.numeric(
      md[[field]]
    )
  )

  q <- quantile(
    x,
    probs = c(
      0,
      0.01,
      0.05,
      0.5,
      0.95,
      0.99,
      1
    ),
    na.rm = TRUE
  )

  tmp <- data.table(
    metric = field,
    min = q[1],
    p01 = q[2],
    p05 = q[3],
    median = q[4],
    p95 = q[5],
    p99 = q[6],
    max = q[7]
  )

  qc_rows[[field]] <- tmp
}

qc_audit <- rbindlist(
  qc_rows,
  fill = TRUE
)

print(
  qc_audit
)

# ------------------------------------------------------------
# 12. Explicit analysis decision
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("ANALYSIS DECISION FOR GSE197677\n")
cat("============================================================\n\n")

cat(
  "Integrated object cells : ",
  nrow(md),
  "\n",
  sep = ""
)

cat(
  "Integrated samples      : ",
  uniqueN(md$sample),
  "\n",
  sep = ""
)

cat(
  "Tumor samples           : ",
  uniqueN(
    md[
      tissue_std == "Tumor",
      sample
    ]
  ),
  "\n",
  sep = ""
)

cat(
  "Normal samples          : ",
  uniqueN(
    md[
      tissue_std == "Normal",
      sample
    ]
  ),
  "\n",
  sep = ""
)

cat(
  "NACT tumor samples      : ",
  uniqueN(
    md[
      tissue_std == "Tumor" &
      treatment_group ==
        "Neoadjuvant_chemotherapy",
      sample
    ]
  ),
  "\n",
  sep = ""
)

cat(
  "nNACT tumor samples     : ",
  uniqueN(
    md[
      tissue_std == "Tumor" &
      treatment_group ==
        "No_neoadjuvant_chemotherapy",
      sample
    ]
  ),
  "\n",
  sep = ""
)

cat(
  "RNA counts available   : ",
  counts_available,
  "\n",
  sep = ""
)

cat(
  "\nPRIMARY TREATMENT REPLICATION:\n",
  "  NACT tumor vs nNACT tumor\n",
  "  using the 10 tumor samples in THIS integrated object.\n",
  sep = ""
)

cat(
  "\nNormal cells are contextual/reference samples,\n",
  "not part of the primary NACT-vs-nNACT treatment contrast.\n"
)

cat(
  "\nNo additional hard QC filtering will be applied now.\n",
  "The published object already shows author-QC boundaries.\n",
  sep = ""
)

# ------------------------------------------------------------
# 13. Save tables
# ------------------------------------------------------------

fwrite(
  sample_table,
  file.path(
    out_dir,
    "GSE197677_standardized_sample_table.csv"
  ),
  bom = TRUE
)

fwrite(
  treatment_summary,
  file.path(
    out_dir,
    "GSE197677_tumor_treatment_group_summary.csv"
  ),
  bom = TRUE
)

fwrite(
  pairing,
  file.path(
    out_dir,
    "GSE197677_internal_pairing_audit.csv"
  ),
  bom = TRUE
)

fwrite(
  qc_audit,
  file.path(
    out_dir,
    "GSE197677_author_QC_boundary_audit.csv"
  ),
  bom = TRUE
)

# Lightweight cell metadata only
key_fields <- intersect(
  c(
    "cell_id",
    "sample",
    "subject_id_internal",
    "tissue_std",
    "treatment_group",
    "treatment_raw",
    "subtype",
    "nCount_RNA",
    "nFeature_RNA",
    "percent.mt",
    "seurat_clusters"
  ),
  names(md)
)

fwrite(
  md[
    ,
    ..key_fields
  ],
  file.path(
    out_dir,
    "GSE197677_standardized_cell_metadata.csv.gz"
  )
)

# ------------------------------------------------------------
# FINAL
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP 14C FINISHED ✓\n")
cat("============================================================\n\n")

cat("No cells filtered ✓\n")
cat("No second QC cutoff applied ✓\n")
cat("No normalization ✓\n")
cat("No PCA / UMAP ✓\n")
cat("No cross-dataset integration ✓\n")

cat("\nSTOP HERE.\n")

