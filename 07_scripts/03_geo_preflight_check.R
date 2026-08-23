# ============================================================
# ESCC Neoadjuvant Project
# Step 6: GEO download preflight check
# Show exact filenames / URLs before downloading
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
})

cat("\n============================================\n")
cat("STEP 6: GEO DOWNLOAD PREFLIGHT CHECK\n")
cat("============================================\n\n")

plan_file <- "00_metadata/01_GEO_download_plan.csv"

if (!file.exists(plan_file)) {
  stop("Cannot find: ", plan_file)
}

plan <- fread(plan_file, encoding = "UTF-8")

needed <- c("GSE", "filename", "file_type", "size_human", "priority", "download_now", "source_url")
missing_cols <- setdiff(needed, names(plan))
if (length(missing_cols) > 0) {
  stop("Missing columns: ", paste(missing_cols, collapse = ", "))
}

# Print exact inventory
cat("\n============================================\n")
cat("ALL GEO SUPPLEMENTARY FILES\n")
cat("============================================\n")

for (gse in unique(plan$GSE)) {
  cat("\n--------------------------------------------\n")
  cat(gse, "\n")
  cat("--------------------------------------------\n")
  x <- plan[GSE == gse, .(filename, file_type, priority, source_url)]
  print(x, row.names = FALSE)
}

# First cohort: GSE221561
first <- plan[GSE == "GSE221561"]
cat("\n============================================\n")
cat("GSE221561 — FIRST DOWNLOAD CANDIDATES\n")
cat("============================================\n\n")

for (i in seq_len(nrow(first))) {
  cat("FILE", i, "\n")
  cat("Filename :", first$filename[i], "\n")
  cat("Type     :", first$file_type[i], "\n")
  cat("URL      :", first$source_url[i], "\n")
  cat("\n")
}

# Save compact preflight table
fwrite(plan[, .(GSE, filename, file_type, priority, source_url)], "00_metadata/02_GEO_preflight_inventory.csv", bom = TRUE)
cat("\nSaved:\n00_metadata/02_GEO_preflight_inventory.csv\n")

cat("\n============================================\n")
cat("STEP 6 R CHECK FINISHED\n")
cat("============================================\n")
