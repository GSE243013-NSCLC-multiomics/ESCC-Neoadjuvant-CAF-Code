# ============================================================
# STEP 17G2B
# Balance GSE221561 S6417T InferCNA pilot references
#
# Uses EXISTING 17G2 raw-count output.
#
# Keep:
#   S6417T tumor = ALL epithelial cells
#   S6417N normal = exactly 500 cells
#   S9478N normal = exactly 500 cells
#
# Fixed random seed for reproducibility.
#
# NO InferCNA run
# NO re-extraction from Seurat
# NO JoinLayers
# NO malignant calling
# ============================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
})

set.seed(22156117)

cat("\n============================================================\n")
cat("STEP 17G2B: BALANCE NORMAL REFERENCES\n")
cat("============================================================\n\n")

indir <- paste0(
  "06_tables/cross_dataset/",
  "malignant_epithelial/",
  "GSE221561_S6417T_pilot"
)

if (!dir.exists(indir)) {
  stop("Pilot directory not found: ", indir)
}

# ------------------------------------------------------------
# 1. Locate existing 17G2 outputs
# ------------------------------------------------------------

all_files <- list.files(
  indir,
  full.names = TRUE
)

raw_hits <- all_files[
  grepl(
    "raw_counts.*\\.rds$",
    basename(all_files),
    ignore.case = TRUE
  )
]

ann_hits <- all_files[
  grepl(
    "cell_annotation.*\\.csv$",
    basename(all_files),
    ignore.case = TRUE
  )
]

gene_hits <- all_files[
  grepl(
    "gene_table.*\\.csv$",
    basename(all_files),
    ignore.case = TRUE
  )
]

if (length(raw_hits) == 0) {
  stop("17G2 raw-count RDS not found.")
}

if (length(ann_hits) == 0) {
  stop("17G2 cell annotation CSV not found.")
}

if (length(gene_hits) == 0) {
  stop("17G2 gene table CSV not found.")
}

raw_file <- raw_hits[1]
ann_file <- ann_hits[1]
gene_file <- gene_hits[1]

cat("Raw counts : ", raw_file, "\n", sep = "")
cat("Annotation : ", ann_file, "\n", sep = "")
cat("Gene table : ", gene_file, "\n\n", sep = "")

# ------------------------------------------------------------
# 2. Load
# ------------------------------------------------------------

counts <- readRDS(raw_file)
ann <- fread(ann_file)
gene_table <- fread(gene_file)

cat(
  "Existing matrix: ",
  nrow(counts),
  " genes x ",
  ncol(counts),
  " cells\n",
  sep = ""
)

cat(
  "Annotation rows: ",
  nrow(ann),
  "\n",
  sep = ""
)

if (anyDuplicated(rownames(counts)) > 0) {
  stop("Existing raw matrix contains duplicate genes.")
}

# ------------------------------------------------------------
# 3. Resolve annotation column names robustly
# ------------------------------------------------------------

sample_col <- intersect(
  c("sample_id", "sample"),
  names(ann)
)

if (length(sample_col) == 0) {
  stop("Cannot find sample/sample_id in annotation.")
}

sample_col <- sample_col[1]

if (!"cell_id" %in% names(ann)) {
  stop("cell_id column missing.")
}

if (!"infercna_role" %in% names(ann)) {
  stop("infercna_role column missing.")
}

ann[
  ,
  sample_resolved :=
    as.character(
      get(sample_col)
    )
]

# ------------------------------------------------------------
# 4. Identify cells
# ------------------------------------------------------------

tumor_cells <- ann[
  infercna_role == "TUMOR_TEST" &
  sample_resolved == "S6417T",
  cell_id
]

s6417n_all <- ann[
  infercna_role == "NORMAL_REFERENCE" &
  sample_resolved == "S6417N",
  cell_id
]

s9478n_all <- ann[
  infercna_role == "NORMAL_REFERENCE" &
  sample_resolved == "S9478N",
  cell_id
]

cat("\n===== EXISTING CELL INVENTORY =====\n\n")

cat(
  "S6417T tumor : ",
  length(tumor_cells),
  "\n",
  sep = ""
)

cat(
  "S6417N normal: ",
  length(s6417n_all),
  "\n",
  sep = ""
)

cat(
  "S9478N normal: ",
  length(s9478n_all),
  "\n",
  sep = ""
)

if (length(tumor_cells) == 0) {
  stop("No S6417T tumor cells found.")
}

if (length(s6417n_all) < 500) {
  stop("S6417N has fewer than 500 available cells.")
}

if (length(s9478n_all) < 500) {
  stop("S9478N has fewer than 500 available cells.")
}

# ------------------------------------------------------------
# 5. Fixed-seed balanced subsampling
# ------------------------------------------------------------

s6417n_ref <- sample(
  s6417n_all,
  500,
  replace = FALSE
)

s9478n_ref <- sample(
  s9478n_all,
  500,
  replace = FALSE
)

refCells <- list(
  S6417N = s6417n_ref,
  S9478N = s9478n_ref
)

selected_cells <- c(
  tumor_cells,
  s6417n_ref,
  s9478n_ref
)

if (anyDuplicated(selected_cells) > 0) {
  stop("Duplicate selected cell IDs detected.")
}

missing_cells <- setdiff(
  selected_cells,
  colnames(counts)
)

if (length(missing_cells) > 0) {
  stop(
    "Selected cells absent from raw matrix: ",
    length(missing_cells)
  )
}

# ------------------------------------------------------------
# 6. Subset RAW counts
# ------------------------------------------------------------

balanced_counts <- counts[
  ,
  selected_cells,
  drop = FALSE
]

cat("\n===== BALANCED RAW MATRIX =====\n\n")

cat(
  nrow(balanced_counts),
  " genes x ",
  ncol(balanced_counts),
  " cells\n",
  sep = ""
)

cat(
  "Tumor  : ",
  length(tumor_cells),
  "\n",
  sep = ""
)

cat(
  "S6417N : ",
  length(s6417n_ref),
  "\n",
  sep = ""
)

cat(
  "S9478N : ",
  length(s9478n_ref),
  "\n",
  sep = ""
)

# ------------------------------------------------------------
# 7. Remove genes zero after balanced subset
# ------------------------------------------------------------

keep_gene <- Matrix::rowSums(
  balanced_counts
) > 0

removed_zero <- sum(!keep_gene)

balanced_counts <- balanced_counts[
  keep_gene,
  ,
  drop = FALSE
]

cat(
  "Zero genes removed after balancing: ",
  removed_zero,
  "\n",
  sep = ""
)

cat(
  "Final genes: ",
  nrow(balanced_counts),
  "\n",
  sep = ""
)

# ------------------------------------------------------------
# 8. Recompute CPM FROM RAW COUNTS
# ------------------------------------------------------------

cat("\n===== CPM RECOMPUTATION =====\n\n")

lib_size <- Matrix::colSums(
  balanced_counts
)

if (any(lib_size <= 0)) {
  stop("Zero-library cell detected.")
}

balanced_cpm <- balanced_counts %*%
  Matrix::Diagonal(
    x = 1e6 / lib_size
  )

rownames(balanced_cpm) <- rownames(balanced_counts)
colnames(balanced_cpm) <- colnames(balanced_counts)

cpm_sums <- Matrix::colSums(
  balanced_cpm
)

cat(
  sprintf(
    "CPM column sum: median %.2f; range %.2f - %.2f\n",
    median(cpm_sums),
    min(cpm_sums),
    max(cpm_sums)
  )
)

if (
  max(
    abs(
      cpm_sums - 1e6
    )
  ) > 1
) {
  stop("CPM validation failed.")
}

cat("CPM validated\n")

# ------------------------------------------------------------
# 9. Balanced annotation
# ------------------------------------------------------------

balanced_ann <- ann[
  cell_id %in% selected_cells
]

balanced_ann[
  ,
  infercna_role :=
    fifelse(
      cell_id %in% tumor_cells,
      "TUMOR_TEST",
      "NORMAL_REFERENCE"
    )
]

balanced_ann[
  ,
  reference_group :=
    fifelse(
      cell_id %in% s6417n_ref,
      "S6417N",
      fifelse(
        cell_id %in% s9478n_ref,
        "S9478N",
        NA_character_
      )
    )
]

# restore exact matrix cell order
balanced_ann[
  ,
  matrix_order :=
    match(
      cell_id,
      colnames(balanced_cpm)
    )
]

setorder(
  balanced_ann,
  matrix_order
)

# ------------------------------------------------------------
# 10. Gene table subset
# ------------------------------------------------------------

if (!"gene" %in% names(gene_table)) {
  stop("gene column missing from gene table.")
}

balanced_gene_table <- gene_table[
  gene %in%
    rownames(balanced_cpm)
]

balanced_gene_table[
  ,
  matrix_order :=
    match(
      gene,
      rownames(balanced_cpm)
    )
]

balanced_gene_table <- balanced_gene_table[
  !is.na(matrix_order)
]

setorder(
  balanced_gene_table,
  matrix_order
)

if (
  nrow(balanced_gene_table) !=
    nrow(balanced_cpm)
) {
  stop(
    "Gene table / matrix row count mismatch."
  )
}

# ------------------------------------------------------------
# 11. Memory estimate
# ------------------------------------------------------------

dense_gb <- (
  nrow(balanced_cpm) *
  ncol(balanced_cpm) *
  8
) / 1024^3

cat("\n===== MEMORY =====\n\n")

cat(
  sprintf(
    "One dense matrix : %.3f GB\n",
    dense_gb
  )
)

cat(
  sprintf(
    "4x working estimate: %.3f GB\n",
    dense_gb * 4
  )
)

# ------------------------------------------------------------
# 12. Save as BALANCED canonical pilot
# ------------------------------------------------------------

cat("\n===== SAVE BALANCED PILOT =====\n\n")

raw_out <- file.path(
  indir,
  "GSE221561_S6417T_BALANCED_infercna_raw_counts_hg38.rds"
)

cpm_out <- file.path(
  indir,
  "GSE221561_S6417T_BALANCED_infercna_CPM_hg38.rds"
)

refs_out <- file.path(
  indir,
  "GSE221561_S6417T_BALANCED_infercna_refCells.rds"
)

ann_out <- file.path(
  indir,
  "GSE221561_S6417T_BALANCED_infercna_cell_annotation.csv"
)

gene_out <- file.path(
  indir,
  "GSE221561_S6417T_BALANCED_infercna_gene_table.csv"
)

manifest_out <- file.path(
  indir,
  "GSE221561_S6417T_BALANCED_infercna_manifest.csv"
)

saveRDS(
  balanced_counts,
  raw_out,
  compress = FALSE
)

saveRDS(
  balanced_cpm,
  cpm_out,
  compress = FALSE
)

saveRDS(
  refCells,
  refs_out
)

fwrite(
  balanced_ann,
  ann_out,
  bom = TRUE
)

fwrite(
  balanced_gene_table,
  gene_out,
  bom = TRUE
)

manifest <- data.table(

  dataset =
    "GSE221561",

  pilot =
    "S6417T",

  tumor_cells =
    length(tumor_cells),

  S6417N_reference_cells =
    length(s6417n_ref),

  S9478N_reference_cells =
    length(s9478n_ref),

  total_reference_cells =
    length(s6417n_ref) +
    length(s9478n_ref),

  total_cells =
    ncol(balanced_cpm),

  genes =
    nrow(balanced_cpm),

  normalization =
    "CPM",

  random_seed =
    22156117L,

  reference_sampling =
    "500 cells per normal epithelial sample",

  infercna_run =
    FALSE
)

fwrite(
  manifest,
  manifest_out,
  bom = TRUE
)

# ------------------------------------------------------------
# 13. Reload validation
# ------------------------------------------------------------

check <- readRDS(cpm_out)
check_refs <- readRDS(refs_out)

if (
  dim(check)[2] !=
    length(tumor_cells) + 1000
) {
  stop("Unexpected final cell number.")
}

if (
  lengths(check_refs)[["S6417N"]] != 500 ||
  lengths(check_refs)[["S9478N"]] != 500
) {
  stop("Reference balancing failed.")
}

if (
  !all(
    unlist(
      check_refs,
      use.names = FALSE
    ) %in%
      colnames(check)
  )
) {
  stop("Reference IDs missing after reload.")
}

cat("\n============================================================\n")
cat("STEP 17G2B COMPLETE\n")
cat("BALANCED GSE221561 PILOT READY\n")
cat("============================================================\n\n")

cat(
  "Tumor S6417T           : ",
  length(tumor_cells),
  "\n",
  sep = ""
)

cat(
  "Normal S6417N          : ",
  length(s6417n_ref),
  "\n",
  sep = ""
)

cat(
  "Normal S9478N          : ",
  length(s9478n_ref),
  "\n",
  sep = ""
)

cat(
  "Total cells            : ",
  ncol(balanced_cpm),
  "\n",
  sep = ""
)

cat(
  "Final hg38 genes       : ",
  nrow(balanced_cpm),
  "\n",
  sep = ""
)

cat(
  sprintf(
    "Dense matrix estimate : %.3f GB\n",
    dense_gb
  )
)

cat("\nInferCNA run             : NO\n")
cat("JoinLayers               : NO\n")
cat("Malignant labels         : NO\n")

cat("\nSTOP HERE.\n")
cat("DO NOT RUN 17G3 YET.\n")
