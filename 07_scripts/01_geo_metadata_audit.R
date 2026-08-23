# ============================================================
# ESCC Neoadjuvant Project
# Step 4: GEO metadata and supplementary-file audit
# ============================================================

options(
  timeout = 1200,
  stringsAsFactors = FALSE
)

suppressPackageStartupMessages({
  library(GEOquery)
  library(data.table)
})

cat("\n============================================\n")
cat("STEP 4: GEO METADATA AUDIT\n")
cat("============================================\n\n")

# ------------------------------------------------------------
# 1. GEO datasets
# ------------------------------------------------------------

gse_ids <- c(
  "GSE221561",
  "GSE197677",
  "GSE160269"
)

dir.create("00_metadata", showWarnings = FALSE, recursive = TRUE)
dir.create("08_logs", showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------
# 2. Helper: safely combine GEO phenotype data
# ------------------------------------------------------------

extract_metadata <- function(gse_object, gse_id) {

  if (!is.list(gse_object)) {
    gse_object <- list(gse_object)
  }

  metadata_list <- list()

  for (i in seq_along(gse_object)) {

    eset <- gse_object[[i]]

    pd <- Biobase::pData(eset)

    pd <- as.data.frame(
      pd,
      stringsAsFactors = FALSE
    )

    pd$geo_accession_source <- rownames(pd)
    pd$GSE <- gse_id
    pd$platform_index <- i

    # Put identifiers first
    first_cols <- c(
      "GSE",
      "geo_accession_source",
      "platform_index"
    )

    pd <- pd[
      ,
      c(
        first_cols,
        setdiff(colnames(pd), first_cols)
      ),
      drop = FALSE
    ]

    metadata_list[[i]] <- pd
  }

  data.table::rbindlist(
    metadata_list,
    fill = TRUE,
    use.names = TRUE
  )
}


# ------------------------------------------------------------
# 3. Audit each GEO dataset
# ------------------------------------------------------------

audit_list <- list()

for (gse_id in gse_ids) {

  cat("\n--------------------------------------------\n")
  cat("Processing:", gse_id, "\n")
  cat("--------------------------------------------\n")

  # ----------------------------------------------------------
  # A. GEO sample metadata
  # ----------------------------------------------------------

  cat("[1/3] Reading GEO Series Matrix metadata...\n")

  geo_result <- tryCatch(

    {
      getGEO(
        gse_id,
        GSEMatrix = TRUE,
        getGPL = FALSE
      )
    },

    error = function(e) {
      cat("ERROR reading metadata:", conditionMessage(e), "\n")
      return(NULL)
    }
  )

  if (!is.null(geo_result)) {

    meta <- extract_metadata(
      geo_result,
      gse_id
    )

    metadata_file <- file.path(
      "00_metadata",
      paste0(gse_id, "_sample_metadata_raw.csv")
    )

    fwrite(
      meta,
      metadata_file,
      bom = TRUE
    )

    cat(
      "Saved metadata:",
      metadata_file,
      "\n"
    )

    cat(
      "Number of GEO sample records:",
      nrow(meta),
      "\n"
    )

    n_platforms <- length(geo_result)

  } else {

    meta <- NULL
    n_platforms <- NA_integer_
  }


  # ----------------------------------------------------------
  # B. Supplementary file LIST ONLY
  # ----------------------------------------------------------

  cat("[2/3] Reading supplementary-file list...\n")
  cat("      Large supplementary files will NOT be downloaded.\n")

  supp <- tryCatch(

    {
      getGEOSuppFiles(
        gse_id,
        fetch_files = FALSE
      )
    },

    error = function(e) {
      cat(
        "ERROR reading supplementary file list:",
        conditionMessage(e),
        "\n"
      )
      return(NULL)
    }
  )

  if (!is.null(supp)) {

    supp <- as.data.frame(
      supp,
      stringsAsFactors = FALSE
    )

    supp_file <- file.path(
      "00_metadata",
      paste0(gse_id, "_supplementary_files.csv")
    )

    fwrite(
      supp,
      supp_file,
      bom = TRUE
    )

    cat(
      "Saved supplementary-file list:",
      supp_file,
      "\n"
    )

    cat(
      "Number of supplementary files:",
      nrow(supp),
      "\n"
    )

    n_supp <- nrow(supp)

  } else {

    n_supp <- NA_integer_
  }


  # ----------------------------------------------------------
  # C. Basic series summary
  # ----------------------------------------------------------

  cat("[3/3] Creating audit summary...\n")

  if (!is.null(meta)) {

    sample_column <- if (
      "geo_accession" %in% colnames(meta)
    ) {
      "geo_accession"
    } else {
      "geo_accession_source"
    }

    n_samples <- length(
      unique(meta[[sample_column]])
    )

  } else {

    n_samples <- NA_integer_
  }

  audit_list[[gse_id]] <- data.frame(
    GSE = gse_id,
    GEO_sample_records = n_samples,
    GEO_platform_objects = n_platforms,
    supplementary_files = n_supp,
    metadata_success = !is.null(meta),
    supplementary_list_success = !is.null(supp),
    audit_date = as.character(Sys.Date()),
    stringsAsFactors = FALSE
  )

  cat("Finished:", gse_id, "\n")
}


# ------------------------------------------------------------
# 4. Save combined audit table
# ------------------------------------------------------------

audit_table <- data.table::rbindlist(
  audit_list,
  fill = TRUE
)

fwrite(
  audit_table,
  "00_metadata/00_GEO_series_audit.csv",
  bom = TRUE
)

cat("\n============================================\n")
cat("GEO SERIES AUDIT SUMMARY\n")
cat("============================================\n\n")

print(audit_table)


# ------------------------------------------------------------
# 5. Show files created
# ------------------------------------------------------------

cat("\n============================================\n")
cat("FILES CREATED\n")
cat("============================================\n\n")

created_files <- list.files(
  "00_metadata",
  pattern = "GSE|GEO_series",
  full.names = TRUE
)

print(created_files)


# ------------------------------------------------------------
# 6. Session info
# ------------------------------------------------------------

sink("08_logs/01_geo_metadata_audit_sessionInfo.txt")
print(sessionInfo())
sink()

cat("\nSession information saved to:\n")
cat("08_logs/01_geo_metadata_audit_sessionInfo.txt\n")

cat("\n============================================\n")
cat("STEP 4 FINISHED\n")
cat("============================================\n")
