# ============================================================
# ESCC Neoadjuvant Project
# Step 3: Install core R packages
# ============================================================

options(
  repos = c(CRAN = "https://cloud.r-project.org"),
  timeout = 600
)

cat("\n===== R version =====\n")
print(R.version.string)

cat("\n===== Install CRAN packages =====\n")

cran_pkgs <- c(
  "BiocManager",
  "Seurat",
  "data.table",
  "dplyr",
  "readr",
  "stringr",
  "tibble",
  "ggplot2",
  "patchwork"
)

for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("\nInstalling:", pkg, "\n")
    install.packages(pkg)
  } else {
    cat("Already installed:", pkg, "\n")
  }
}

cat("\n===== Install Bioconductor packages =====\n")

bioc_pkgs <- c(
  "GEOquery",
  "SingleCellExperiment",
  "edgeR",
  "limma"
)

for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("\nInstalling:", pkg, "\n")
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  } else {
    cat("Already installed:", pkg, "\n")
  }
}

cat("\n===== Package verification =====\n")

check_pkgs <- c(
  "Seurat",
  "GEOquery",
  "SingleCellExperiment",
  "edgeR",
  "limma",
  "data.table",
  "dplyr",
  "ggplot2"
)

for (pkg in check_pkgs) {
  ok <- requireNamespace(pkg, quietly = TRUE)

  if (ok) {
    cat(
      sprintf(
        "[OK] %-22s version %s\n",
        pkg,
        as.character(packageVersion(pkg))
      )
    )
  } else {
    cat(sprintf("[FAILED] %s\n", pkg))
  }
}

cat("\n===== Session information =====\n")
print(sessionInfo())

cat("\n===== STEP 3 R PACKAGE INSTALLATION FINISHED =====\n")
