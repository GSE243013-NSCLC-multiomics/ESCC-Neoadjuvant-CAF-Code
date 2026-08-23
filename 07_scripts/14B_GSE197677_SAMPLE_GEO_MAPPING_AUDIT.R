# ============================================================
# ESCC Neoadjuvant Project
# STEP 14B
# GSE197677 sample / GEO / treatment mapping audit
#
# AUDIT ONLY
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
cat("STEP 14B: GSE197677 SAMPLE / GEO MAPPING AUDIT\n")
cat("============================================================\n\n")

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------

object_dir <- "03_objects/GSE197677"

geo_meta_file <- paste0(
  "00_metadata/",
  "GSE197677_sample_metadata_raw.csv"
)

out_dir <- "06_tables/GSE197677/metadata"

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# 1. Locate integrated object
# ------------------------------------------------------------

rds_files <- list.files(
  object_dir,
  pattern = "\\.rds$",
  full.names = TRUE,
  ignore.case = TRUE
)

integrated_file <- rds_files[
  grepl(
    "integrated",
    basename(rds_files),
    ignore.case = TRUE
  )
]

if (length(integrated_file) == 0) {
  stop("GSE197677 integrated RDS not found.")
}

integrated_file <- integrated_file[1]

cat("===== 1. OBJECT =====\n\n")
cat("File:", integrated_file, "\n")

obj <- readRDS(
  integrated_file
)

md <- as.data.table(
  obj[[]],
  keep.rownames = "cell_id"
)

cat("Cells from metadata:", nrow(md), "\n")
cat("Metadata fields    :", ncol(md), "\n\n")

# ------------------------------------------------------------
# 2. Required known fields
# ------------------------------------------------------------

known_fields <- c(
  "sample",
  "CellType",
  "NAC",
  "subtype"
)

missing_known <- setdiff(
  known_fields,
  names(md)
)

if (length(missing_known) > 0) {

  cat(
    "WARNING: expected fields missing: ",
    paste(missing_known, collapse = ", "),
    "\n",
    sep = ""
  )
}

# ------------------------------------------------------------
# 3. Print ALL metadata columns
# ------------------------------------------------------------

cat("===== 2. ALL METADATA FIELDS =====\n\n")

for (i in seq_along(names(md))) {

  x <- md[[i]]

  cat(
    sprintf(
      "%02d. %-30s unique=%d missing=%d\n",
      i,
      names(md)[i],
      uniqueN(x, na.rm = TRUE),
      sum(is.na(x))
    )
  )
}

# ------------------------------------------------------------
# 4. Object-level sample audit
# ------------------------------------------------------------

cat("\n===== 3. OBJECT SAMPLE AUDIT =====\n\n")

if (!"sample" %in% names(md)) {
  stop("Field 'sample' not present.")
}

sample_ids <- sort(
  unique(
    as.character(md$sample)
  )
)

cat(
  "Unique object sample IDs:",
  length(sample_ids),
  "\n\n"
)

print(sample_ids)

# ------------------------------------------------------------
# 5. Sample x phenotype table
# ------------------------------------------------------------

collapse_unique <- function(x) {

  x <- unique(
    as.character(x)
  )

  x <- x[
    !is.na(x) &
    trimws(x) != ""
  ]

  if (length(x) == 0) {
    return(NA_character_)
  }

  paste(
    sort(x),
    collapse = "; "
  )
}

sample_summary <- md[
  ,
  .(
    n_cells = .N,

    CellType =
      if ("CellType" %in% names(md))
        collapse_unique(CellType)
      else NA_character_,

    NAC =
      if ("NAC" %in% names(md))
        collapse_unique(NAC)
      else NA_character_,

    subtype =
      if ("subtype" %in% names(md))
        collapse_unique(subtype)
      else NA_character_,

    orig_ident =
      if ("orig.ident" %in% names(md))
        collapse_unique(orig.ident)
      else NA_character_
  ),
  by = sample
]

setorder(
  sample_summary,
  sample
)

cat("\n===== 4. SAMPLE-LEVEL PHENOTYPES =====\n\n")

print(sample_summary)

fwrite(
  sample_summary,
  file.path(
    out_dir,
    "GSE197677_object_sample_summary.csv"
  ),
  bom = TRUE
)

# ------------------------------------------------------------
# 6. Cell counts by sample / CellType / NAC
# ------------------------------------------------------------

group_fields <- intersect(
  c(
    "sample",
    "CellType",
    "NAC",
    "subtype"
  ),
  names(md)
)

cell_counts <- md[
  ,
  .N,
  by = group_fields
]

setnames(
  cell_counts,
  "N",
  "n_cells"
)

setorderv(
  cell_counts,
  group_fields
)

cat("\n===== 5. SAMPLE x TISSUE x NAC COUNTS =====\n\n")

print(cell_counts)

fwrite(
  cell_counts,
  file.path(
    out_dir,
    "GSE197677_sample_tissue_NAC_cell_counts.csv"
  ),
  bom = TRUE
)

# ------------------------------------------------------------
# 7. Cohort counts based on object
# ------------------------------------------------------------

cat("\n===== 6. OBJECT COHORT STRUCTURE =====\n\n")

if ("CellType" %in% names(md)) {

  object_tissue <- unique(
    md[
      ,
      .(
        sample,
        CellType
      )
    ]
  )

  tissue_summary <- object_tissue[
    ,
    .(
      n_samples = uniqueN(sample)
    ),
    by = CellType
  ]

  print(tissue_summary)

} else {

  tissue_summary <- data.table()
}

if ("NAC" %in% names(md)) {

  object_nac <- unique(
    md[
      ,
      .(
        sample,
        NAC
      )
    ]
  )

  nac_summary <- object_nac[
    ,
    .(
      n_samples = uniqueN(sample),
      samples = paste(
        sort(
          unique(sample)
        ),
        collapse = ", "
      )
    ),
    by = NAC
  ]

  cat("\nNAC groups:\n")

  print(nac_summary)

} else {

  nac_summary <- data.table()
}

# ------------------------------------------------------------
# 8. Load the 30-record GEO metadata audit
# ------------------------------------------------------------

cat("\n===== 7. GEO SERIES METADATA =====\n\n")

if (!file.exists(geo_meta_file)) {

  stop(
    "Cannot find GEO metadata file: ",
    geo_meta_file
  )
}

geo <- fread(
  geo_meta_file,
  encoding = "UTF-8"
)

cat(
  "GEO metadata rows:",
  nrow(geo),
  "\n"
)

cat(
  "GEO metadata columns:",
  ncol(geo),
  "\n"
)

if (nrow(geo) != 30) {

  cat(
    "WARNING: expected 30 GEO sample records, found ",
    nrow(geo),
    "\n",
    sep = ""
  )
}

# ------------------------------------------------------------
# 9. GEO accession field
# ------------------------------------------------------------

geo_id_col <- NULL

for (
  candidate in c(
    "geo_accession",
    "geo_accession_source"
  )
) {

  if (candidate %in% names(geo)) {
    geo_id_col <- candidate
    break
  }
}

if (is.null(geo_id_col)) {

  geo[, GEO_record :=
        paste0(
          "GEO_ROW_",
          seq_len(.N)
        )]

  geo_id_col <- "GEO_record"
}

# ------------------------------------------------------------
# 10. Search ALL GEO metadata fields for object sample IDs
# ------------------------------------------------------------

cat("\n===== 8. OBJECT SAMPLE IDs vs GEO METADATA =====\n\n")

search_sample_in_geo <- function(sample_id) {

  hits <- rep(
    FALSE,
    nrow(geo)
  )

  hit_fields <- vector(
    "list",
    nrow(geo)
  )

  for (col in names(geo)) {

    x <- as.character(
      geo[[col]]
    )

    h <- grepl(
      sample_id,
      x,
      fixed = TRUE
    )

    hits <- hits | h

    if (any(h, na.rm = TRUE)) {

      for (i in which(h)) {

        hit_fields[[i]] <- c(
          hit_fields[[i]],
          col
        )
      }
    }
  }

  idx <- which(
    hits
  )

  if (length(idx) == 0) {

    return(
      data.table(
        sample = sample_id,
        GEO_accession = NA_character_,
        matched_fields = NA_character_
      )
    )
  }

  rbindlist(
    lapply(
      idx,
      function(i) {

        data.table(
          sample = sample_id,

          GEO_accession =
            as.character(
              geo[[geo_id_col]][i]
            ),

          matched_fields =
            paste(
              unique(
                hit_fields[[i]]
              ),
              collapse = "; "
            )
        )
      }
    )
  )
}

geo_matches <- rbindlist(
  lapply(
    sample_ids,
    search_sample_in_geo
  ),
  fill = TRUE
)

setorder(
  geo_matches,
  sample,
  GEO_accession
)

print(
  geo_matches
)

fwrite(
  geo_matches,
  file.path(
    out_dir,
    "GSE197677_object_sample_to_GEO_matches.csv"
  ),
  bom = TRUE
)

# ------------------------------------------------------------
# 11. Match count per object sample
# ------------------------------------------------------------

match_count <- geo_matches[
  ,
  .(
    GEO_records_matched =
      sum(
        !is.na(GEO_accession)
      ),

    GEO_accessions =
      paste(
        sort(
          unique(
            GEO_accession[
              !is.na(GEO_accession)
            ]
          )
        ),
        collapse = "; "
      )
  ),
  by = sample
]

object_geo_audit <- merge(
  sample_summary,
  match_count,
  by = "sample",
  all.x = TRUE
)

object_geo_audit[
  is.na(GEO_records_matched),
  GEO_records_matched := 0
]

cat("\n===== 9. OBJECT ↔ GEO MAPPING SUMMARY =====\n\n")

print(
  object_geo_audit
)

fwrite(
  object_geo_audit,
  file.path(
    out_dir,
    "GSE197677_object_GEO_mapping_audit.csv"
  ),
  bom = TRUE
)

# ------------------------------------------------------------
# 12. GEO records NOT represented in object
# ------------------------------------------------------------

matched_geo_ids <- unique(
  geo_matches$GEO_accession[
    !is.na(
      geo_matches$GEO_accession
    )
  ]
)

all_geo_ids <- unique(
  as.character(
    geo[[geo_id_col]]
  )
)

geo_not_in_object <- setdiff(
  all_geo_ids,
  matched_geo_ids
)

cat("\n===== 10. GEO RECORDS NOT MAPPED TO OBJECT SAMPLE IDs =====\n\n")

cat(
  "Total GEO records       :",
  length(all_geo_ids),
  "\n"
)

cat(
  "GEO records mapped      :",
  length(matched_geo_ids),
  "\n"
)

cat(
  "GEO records not mapped  :",
  length(geo_not_in_object),
  "\n\n"
)

print(
  geo_not_in_object
)

if (length(geo_not_in_object) > 0) {

  geo_missing_rows <- geo[
    get(geo_id_col) %in%
      geo_not_in_object
  ]

  fwrite(
    geo_missing_rows,
    file.path(
      out_dir,
      "GSE197677_GEO_records_not_in_integrated_object.csv"
    ),
    bom = TRUE
  )
}

# ------------------------------------------------------------
# 13. Search GEO columns for tumor/normal/NAC-related text
# ------------------------------------------------------------

cat("\n===== 11. GEO COHORT-RELATED FIELDS =====\n\n")

interesting_geo_fields <- list()

for (col in names(geo)) {

  vals <- unique(
    trimws(
      as.character(
        geo[[col]]
      )
    )
  )

  vals_low <- tolower(
    vals
  )

  if (
    any(
      grepl(
        "nac|neoadjuvant|tumou?r|normal|cancer",
        vals_low
      ),
      na.rm = TRUE
    )
  ) {

    interesting_geo_fields[[col]] <- vals

    cat(
      "\nFIELD: ",
      col,
      "\n",
      sep = ""
    )

    print(
      sort(
        table(
          geo[[col]],
          useNA = "ifany"
        ),
        decreasing = TRUE
      )
    )
  }
}

# ------------------------------------------------------------
# 14. Determine whether integrated RDS appears to be subset
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("INTEGRATED OBJECT COHORT SUMMARY\n")
cat("============================================================\n\n")

cat(
  "Cells in integrated object :",
  nrow(md),
  "\n"
)

cat(
  "Object sample IDs          :",
  length(sample_ids),
  "\n"
)

cat(
  "GEO Series records         :",
  nrow(geo),
  "\n"
)

if ("CellType" %in% names(md)) {

  scc_samples <- uniqueN(
    md[
      tolower(
        as.character(CellType)
      ) == "scc",
      sample
    ]
  )

  normal_samples <- uniqueN(
    md[
      tolower(
        as.character(CellType)
      ) == "normal",
      sample
    ]
  )

  cat(
    "Object SCC samples        :",
    scc_samples,
    "\n"
  )

  cat(
    "Object normal samples     :",
    normal_samples,
    "\n"
  )
}

if ("NAC" %in% names(md)) {

  cat("\nObject NAC sample counts:\n")

  print(
    unique(
      md[
        ,
        .(
          sample,
          NAC
        )
      ]
    )[
      ,
      .(
        n_samples =
          uniqueN(sample)
      ),
      by = NAC
    ]
  )
}

# ------------------------------------------------------------
# 15. Author-QC appearance
# ------------------------------------------------------------

cat("\n===== 12. AUTHOR-QC APPEARANCE =====\n\n")

for (
  col in intersect(
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
      md[[col]]
    )
  )

  cat(
    col,
    ":\n"
  )

  print(
    quantile(
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
  )

  cat("\n")
}

cat(
  "NOTE: If minima/maxima align with conventional QC cutoffs\n",
  "(e.g. nFeature >=200, percent.mt <=10), this strongly suggests\n",
  "the published integrated object is already QC-filtered.\n"
)

# ------------------------------------------------------------
# 16. Save lightweight cell-level key metadata
# ------------------------------------------------------------

keep_fields <- intersect(
  c(
    "cell_id",
    "sample",
    "CellType",
    "NAC",
    "subtype",
    "orig.ident",
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
    ..keep_fields
  ],
  file.path(
    out_dir,
    "GSE197677_cell_metadata_key_fields.csv.gz"
  )
)

# ------------------------------------------------------------
# FINAL
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP 14B FINISHED ✓\n")
cat("============================================================\n\n")

cat("AUDIT ONLY.\n\n")

cat("No cells filtered ✓\n")
cat("No QC thresholds applied ✓\n")
cat("No normalization ✓\n")
cat("No PCA / UMAP ✓\n")
cat("No cross-dataset integration ✓\n")

cat("\nSTOP HERE.\n")

