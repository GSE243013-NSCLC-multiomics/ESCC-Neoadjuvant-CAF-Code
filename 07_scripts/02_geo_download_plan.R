# ============================================================
# ESCC Neoadjuvant Project
# Step 5: Build GEO supplementary-file download plan
# NO large data will be downloaded in this step
# ============================================================

options(
  stringsAsFactors = FALSE,
  timeout = 1200
)

suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
})

cat("\n============================================\n")
cat("STEP 5: GEO DOWNLOAD PLAN\n")
cat("============================================\n\n")

gse_ids <- c(
  "GSE221561",
  "GSE197677",
  "GSE160269"
)

# ------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------

extract_url <- function(dt) {

  for (nm in names(dt)) {
    x <- as.character(dt[[nm]])
    hit <- grepl("^(ftp|https?)://", x, ignore.case = TRUE)
    if (any(hit, na.rm = TRUE)) {
      return(x)
    }
  }

  rn <- rownames(dt)
  if (!is.null(rn) && any(grepl("^(ftp|https?)://", rn, ignore.case = TRUE))) {
    return(rn)
  }

  rep(NA_character_, nrow(dt))
}


extract_filename <- function(dt, url) {
  candidate_names <- c("fname", "filename", "file_name", "name")

  for (nm in candidate_names) {
    if (nm %in% names(dt)) {
      x <- as.character(dt[[nm]])
      if (any(!is.na(x) & nzchar(x))) {
        return(basename(x))
      }
    }
  }

  out <- basename(url)
  bad <- is.na(out) | out == "." | out == ""
  if (any(bad)) {
    rn <- rownames(dt)
    if (!is.null(rn)) {
      out[bad] <- basename(rn[bad])
    }
  }
  out
}


extract_size <- function(dt) {
  possible <- c("size", "filesize", "file_size", "bytes")

  for (nm in possible) {
    if (nm %in% names(dt)) {
      raw <- dt[[nm]]
      if (is.numeric(raw)) return(as.numeric(raw))
      y <- suppressWarnings(as.numeric(gsub("[^0-9.]", "", as.character(raw))))
      if (any(!is.na(y))) return(y)
    }
  }
  rep(NA_real_, nrow(dt))
}


human_size <- function(bytes) {
  vapply(bytes, function(x) {
    if (is.na(x)) return(NA_character_)
    if (x >= 1024^3) sprintf("%.2f GB", x / 1024^3)
    else if (x >= 1024^2) sprintf("%.2f MB", x / 1024^2)
    else if (x >= 1024) sprintf("%.2f KB", x / 1024)
    else sprintf("%.0f B", x)
  }, character(1))
}


classify_file <- function(filename) {
  x <- tolower(filename)
  fifelse(grepl("\\.rds(\\.gz)?$", x), "RDS object",
  fifelse(grepl("\\.h5ad(\\.gz)?$", x), "AnnData/H5AD",
  fifelse(grepl("\\.(h5|hdf5)(\\.gz)?$", x), "HDF5 matrix",
  fifelse(grepl("matrix.*\\.mtx|\\.mtx(\\.gz)?$", x), "MTX matrix",
  fifelse(grepl("barcode", x), "Barcodes",
  fifelse(grepl("feature|gene", x) & grepl("\\.(tsv|txt|csv)(\\.gz)?$", x), "Features/genes",
  fifelse(grepl("annot|metadata|meta|celltype|cell_type|cluster", x), "Annotation/metadata",
  fifelse(grepl("\\.(csv|tsv|txt)(\\.gz)?$", x), "Tabular data",
  fifelse(grepl("\\.(tar|tar\\.gz|tgz|zip)$", x), "Archive",
  "Other")))))))))
}


assign_priority <- function(gse, filename, file_type) {
  x <- tolower(filename)

  if (gse == "GSE221561") {
    if (file_type %in% c("RDS object", "HDF5 matrix", "MTX matrix", "Barcodes", "Features/genes", "Annotation/metadata", "Tabular data", "Archive")) {
      return("1_FIRST")
    }
    return("2_REVIEW")
  }

  if (gse == "GSE197677") {
    if (grepl("rds", x) || file_type %in% c("RDS object", "Annotation/metadata", "Tabular data")) {
      return("2_SECOND")
    }
    return("3_REVIEW")
  }

  if (gse == "GSE160269") {
    if (file_type %in% c("RDS object", "HDF5 matrix", "MTX matrix", "Barcodes", "Features/genes", "Annotation/metadata", "Tabular data", "Archive")) {
      return("3_THIRD")
    }
  }

  "4_REVIEW"
}


# ------------------------------------------------------------
# Read supplementary-file inventories
# ------------------------------------------------------------

all_files <- list()

for (gse in gse_ids) {
  f <- file.path("00_metadata", paste0(gse, "_supplementary_files.csv"))
  if (!file.exists(f)) stop("Missing file: ", f)

  cat("Reading:", f, "\n")
  x <- fread(f, encoding = "UTF-8")

  url <- extract_url(x)
  filename <- extract_filename(x, url)
  size_bytes <- extract_size(x)

  tmp <- data.table(
    GSE = gse,
    filename = filename,
    file_type = classify_file(filename),
    size_bytes = size_bytes,
    size_human = human_size(size_bytes),
    source_url = url
  )

  tmp[, priority := mapply(assign_priority, GSE, filename, file_type, USE.NAMES = FALSE)]
  all_files[[gse]] <- tmp
}

plan <- rbindlist(all_files, fill = TRUE)

# Add download decision column
plan[, download_now := fifelse(priority %in% c("1_FIRST", "2_SECOND"), "CANDIDATE", "LATER")]
plan[file_type == "Other", download_now := "REVIEW"]

# Sort
priority_order <- c("1_FIRST", "2_SECOND", "3_THIRD", "4_REVIEW")
plan[, priority := factor(priority, levels = priority_order)]
setorder(plan, priority, GSE, filename)
plan[, priority := as.character(priority)]

# Save machine-readable plan
outfile <- "00_metadata/01_GEO_download_plan.csv"
fwrite(plan, outfile, bom = TRUE)
cat("\nSaved:\n", outfile, "\n")

# Create simplified review table
review <- plan[, .(GSE, filename, file_type, size_human, priority, download_now)]
review_file <- "00_metadata/01_GEO_download_plan_review.csv"
fwrite(review, review_file, bom = TRUE)

# Print exact inventory
cat("\n============================================\n")
cat("EXACT SUPPLEMENTARY FILE INVENTORY\n")
cat("============================================\n\n")

for (gse in gse_ids) {
  cat("\n--------------------------------------------\n")
  cat(gse, "\n")
  cat("--------------------------------------------\n")
  print(review[GSE == gse], row.names = FALSE)
}

# Summary
cat("\n============================================\n")
cat("DOWNLOAD PLAN SUMMARY\n")
cat("============================================\n\n")

summary_table <- plan[, .(files = .N, known_size_GB = round(sum(size_bytes, na.rm = TRUE) / 1024^3, 3)), by = .(GSE, priority)]
print(summary_table)

cat("\n============================================\n")
cat("IMPORTANT\n")
cat("============================================\n")
cat("No supplementary data files were downloaded.\n")
cat("This step only created the download inventory.\n")

# Session info
sink("08_logs/02_geo_download_plan_sessionInfo.txt")
print(sessionInfo())
sink()

cat("\n============================================\n")
cat("STEP 5 FINISHED\n")
cat("============================================\n")
