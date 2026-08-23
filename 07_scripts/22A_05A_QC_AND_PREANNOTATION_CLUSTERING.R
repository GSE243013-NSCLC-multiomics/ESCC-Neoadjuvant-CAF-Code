# ============================================================
# ESCC Neoadjuvant Project
# STEP22A-5A
#
# Paired tumor:
#   - construct sparse Seurat object
#   - sample-specific technical QC
#   - normalization
#   - PCA
#   - graph clustering
#   - UMAP
#
# IMPORTANT:
#   NO cell-type annotation in this step
#   NO fibroblast selection
#   NO CORE2 / CORE4 scoring
#   NO CDKN1B / GSN / BGN / TIMP1 selection
#   NO treatment-effect inference
#   NO manuscript modification
# ============================================================

options(
  stringsAsFactors = FALSE,
  timeout = 1200,
  future.globals.maxSize = 8 * 1024^3
)

set.seed(20260814)

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(Seurat)
  library(ggplot2)
})

cat("\n")
cat("====================================================================\n")
cat("STEP22A-5A: QC + PRE-ANNOTATION CLUSTERING\n")
cat("====================================================================\n\n")

# ============================================================
# 0. Paths
# ============================================================

PROJECT <- getwd()

META_DIR <- file.path(
  PROJECT,
  "00_metadata",
  "OMIX005710"
)

OBJECT_DIR <- file.path(
  PROJECT,
  "03_objects",
  "OMIX005710"
)

RESULT_DIR <- file.path(
  PROJECT,
  "04_results",
  "external_validation",
  "Step22A"
)

FIGURE_DIR <- file.path(
  PROJECT,
  "05_figures",
  "external_validation",
  "Step22A"
)

LOG_DIR <- file.path(
  PROJECT,
  "08_logs"
)

dir.create(
  OBJECT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  RESULT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  FIGURE_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  LOG_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

MATRIX_AUDIT_FILE <- file.path(
  META_DIR,
  "STEP22A_03B_MATRIX_AUDIT.csv"
)

GATE_FILE <- file.path(
  META_DIR,
  "STEP22A_ANALYSIS_GATE_STATUS.csv"
)

FROZEN_INVENTORY <- file.path(
  META_DIR,
  "STEP22A_PREANALYSIS_FROZEN_FILE_INVENTORY.csv"
)

OUT_OBJECT <- file.path(
  OBJECT_DIR,
  "OMIX005710_PAIRED_QC_CLUSTERED_PREANNOTATION.rds"
)

# ============================================================
# 1. Histology / analysis gate
# ============================================================

cat("===== 1. ANALYSIS GATE CHECK =====\n\n")

if (!file.exists(GATE_FILE)) {
  stop(
    "Missing gate file: ",
    GATE_FILE
  )
}

gate <- fread(
  GATE_FILE,
  encoding = "UTF-8"
)

required_gate_rows <- c(
  "technical_preprocessing",
  "broad_cell_annotation",
  "paired_CORE2_test",
  "manuscript_update"
)

missing_gate_rows <- setdiff(
  required_gate_rows,
  gate$gate
)

if (length(missing_gate_rows) > 0) {
  stop(
    "Missing gate rows: ",
    paste(
      missing_gate_rows,
      collapse = ", "
    )
  )
}

get_allowed <- function(x) {
  gate[
    gate == x,
    allowed
  ][1]
}

technical_allowed <- get_allowed(
  "technical_preprocessing"
)

annotation_allowed <- get_allowed(
  "broad_cell_annotation"
)

core2_allowed <- get_allowed(
  "paired_CORE2_test"
)

manuscript_allowed <- get_allowed(
  "manuscript_update"
)

cat(
  "Technical preprocessing :",
  technical_allowed,
  "\n"
)

cat(
  "Broad annotation        :",
  annotation_allowed,
  "\n"
)

cat(
  "Paired CORE2 test       :",
  core2_allowed,
  "\n"
)

cat(
  "Manuscript update       :",
  manuscript_allowed,
  "\n\n"
)

if (technical_allowed != "YES") {
  stop(
    "TECHNICAL_PREPROCESSING_NOT_ALLOWED"
  )
}

if (annotation_allowed != "YES") {
  stop(
    "BROAD_ANNOTATION_GATE_NOT_OPEN"
  )
}

if (core2_allowed == "YES") {
  stop(
    "ERROR: CORE2 inference should still be blocked."
  )
}

if (manuscript_allowed == "YES") {
  stop(
    "ERROR: manuscript gate should still be blocked."
  )
}

cat(
  "✓ Histology gate respected.\n"
)

cat(
  "✓ Technical preprocessing may continue.\n"
)

cat(
  "✓ CORE2 inference remains blocked.\n\n"
)

# ============================================================
# 2. Load 12-sample matrix audit
# ============================================================

cat("===== 2. LOAD MATRIX AUDIT =====\n\n")

if (!file.exists(MATRIX_AUDIT_FILE)) {
  stop(
    "Missing matrix audit: ",
    MATRIX_AUDIT_FILE
  )
}

audit <- fread(
  MATRIX_AUDIT_FILE,
  encoding = "UTF-8"
)

if (nrow(audit) != 12) {
  stop(
    "Expected 12 paired tumor samples, found ",
    nrow(audit)
  )
}

expected_patients <- c(
  "P6",
  "P9",
  "P15",
  "P16",
  "P18",
  "P19"
)

observed_patients <- sort(
  unique(
    audit$patient_id
  )
)

if (!setequal(
  expected_patients,
  observed_patients
)) {
  stop(
    "PAIRED_PATIENT_SET_MISMATCH"
  )
}

cat(
  "Samples  :",
  nrow(audit),
  "\n"
)

cat(
  "Patients :",
  length(observed_patients),
  "\n"
)

cat(
  "Patients :",
  paste(
    expected_patients,
    collapse = ", "
  ),
  "\n\n"
)

# ============================================================
# 3. Read 12 sparse matrices
# ============================================================

cat("===== 3. READ 12 SPARSE COUNT MATRICES =====\n\n")

matrix_list <- list()
cell_meta_list <- list()

reference_gene_id <- NULL
reference_gene_symbol <- NULL
reference_feature_name <- NULL

read_gene_table <- function(path) {

  x <- fread(
    path,
    header = FALSE,
    encoding = "UTF-8"
  )

  if (ncol(x) < 1) {
    stop(
      "Empty gene table: ",
      path
    )
  }

  gene_id <- as.character(
    x[[1]]
  )

  if (ncol(x) >= 2) {

    gene_symbol <- as.character(
      x[[2]]
    )

  } else {

    gene_symbol <- gene_id
  }

  bad_symbol <- (
    is.na(gene_symbol) |
    gene_symbol == ""
  )

  gene_symbol[
    bad_symbol
  ] <- gene_id[
    bad_symbol
  ]

  feature_name <- make.unique(
    gene_symbol,
    sep = "__dup"
  )

  list(
    gene_id = gene_id,
    gene_symbol = gene_symbol,
    feature_name = feature_name
  )
}

for (i in seq_len(nrow(audit))) {

  row <- audit[i]

  sample_name <- row$sample_name
  patient_id <- row$patient_id
  timepoint <- row$timepoint
  hrs <- row$hrs_accession

  cat(
    sprintf(
      "[%02d/12] %-4s %-6s %-10s ",
      i,
      patient_id,
      timepoint,
      sample_name
    )
  )

  matrix_path <- file.path(
    PROJECT,
    row$matrix_file
  )

  barcode_path <- file.path(
    PROJECT,
    row$barcodes_file
  )

  gene_path <- file.path(
    PROJECT,
    row$genes_file
  )

  for (f in c(
    matrix_path,
    barcode_path,
    gene_path
  )) {

    if (!file.exists(f)) {
      stop(
        "\nMissing input file: ",
        f
      )
    }
  }

  gene_info <- read_gene_table(
    gene_path
  )

  if (is.null(reference_gene_id)) {

    reference_gene_id <- gene_info$gene_id
    reference_gene_symbol <- gene_info$gene_symbol
    reference_feature_name <- gene_info$feature_name

  } else {

    if (!identical(
      gene_info$gene_id,
      reference_gene_id
    )) {
      stop(
        "\nFEATURE_TABLE_MISMATCH: gene IDs differ across samples."
      )
    }

    if (!identical(
      gene_info$gene_symbol,
      reference_gene_symbol
    )) {
      stop(
        "\nFEATURE_TABLE_MISMATCH: gene symbols differ across samples."
      )
    }
  }

  barcodes <- fread(
    barcode_path,
    header = FALSE,
    encoding = "UTF-8"
  )[[1]]

  barcodes <- as.character(
    barcodes
  )

  counts <- Matrix::readMM(
    matrix_path
  )

  counts <- as(
    counts,
    "CsparseMatrix"
  )

  if (
    nrow(counts) != length(
      reference_feature_name
    )
  ) {
    stop(
      "\nGene dimension mismatch: ",
      sample_name
    )
  }

  if (
    ncol(counts) != length(
      barcodes
    )
  ) {
    stop(
      "\nBarcode dimension mismatch: ",
      sample_name
    )
  }

  rownames(counts) <- reference_feature_name

  cell_ids <- paste0(
    sample_name,
    "__",
    barcodes
  )

  colnames(counts) <- cell_ids

  matrix_list[[
    sample_name
  ]] <- counts

  cell_meta_list[[
    sample_name
  ]] <- data.table(
    cell_id = cell_ids,
    barcode_raw = barcodes,
    patient_id = patient_id,
    timepoint = timepoint,
    sample_name = sample_name,
    hrs_accession = hrs
  )

  cat(
    nrow(counts),
    "genes x",
    ncol(counts),
    "cells\n"
  )
}

cat("\n")
cat(
  "✓ All 12 matrices use an identical feature table.\n\n"
)

# ============================================================
# 4. Combine sparse matrices
# ============================================================

cat("===== 4. COMBINE PAIRED MATRICES =====\n\n")

counts_all <- Reduce(
  function(a, b) {
    cbind(
      a,
      b
    )
  },
  matrix_list
)

meta_all <- rbindlist(
  cell_meta_list,
  use.names = TRUE,
  fill = TRUE
)

setkey(
  meta_all,
  cell_id
)

meta_all <- meta_all[
  colnames(
    counts_all
  )
]

if (!identical(
  meta_all$cell_id,
  colnames(counts_all)
)) {
  stop(
    "CELL_METADATA_ORDER_MISMATCH"
  )
}

cat(
  "Combined matrix:",
  nrow(counts_all),
  "genes x",
  ncol(counts_all),
  "cells\n"
)

if (
  ncol(counts_all) != 90699
) {

  cat(
    "\nWARNING: expected 90,699 raw barcodes from Step22A-3B, observed ",
    ncol(counts_all),
    ".\n",
    sep = ""
  )
}

cat("\n")

# ============================================================
# 5. Create Seurat object
# ============================================================

cat("===== 5. CREATE SEURAT OBJECT =====\n\n")

metadata_df <- as.data.frame(
  meta_all
)

rownames(
  metadata_df
) <- metadata_df$cell_id

obj <- CreateSeuratObject(
  counts = counts_all,
  assay = "RNA",
  project = "OMIX005710",
  min.cells = 0,
  min.features = 0,
  meta.data = metadata_df
)

rm(
  counts_all,
  matrix_list,
  cell_meta_list
)

gc()

mito_features <- grep(
  "^MT-",
  rownames(obj),
  ignore.case = TRUE,
  value = TRUE
)

cat(
  "Mitochondrial genes detected:",
  length(mito_features),
  "\n"
)

if (length(mito_features) > 0) {

  obj[[
    "percent.mt"
  ]] <- PercentageFeatureSet(
    obj,
    features = mito_features
  )

} else {

  obj$percent.mt <- 0

  cat(
    "WARNING: no MT-* genes detected; percent.mt set to 0.\n"
  )
}

cat(
  "Raw cells in Seurat object:",
  ncol(obj),
  "\n\n"
)

# ============================================================
# 6. Frozen sample-specific technical QC
#
# This is intentionally LENIENT.
#
# No treatment group is used to choose thresholds.
#
# Lower technical limits:
#   nFeature >= 200
#   nCount   >= 500
#
# Sample-specific upper/outlier limits:
#   robust median + 5 MAD
#
# Mito threshold:
#   max(20%, median + 5 MAD), capped at 40%
#
# This performs only technical QC / obvious upper-outlier
# screening. It is NOT a claim of complete doublet or ambient
# RNA removal.
# ============================================================

cat("===== 6. SAMPLE-SPECIFIC TECHNICAL QC =====\n\n")

robust_scale <- function(x) {

  x <- x[
    is.finite(x)
  ]

  m <- mad(
    x,
    constant = 1.4826,
    na.rm = TRUE
  )

  if (
    !is.finite(m) ||
    m <= 0
  ) {

    m <- IQR(
      x,
      na.rm = TRUE
    ) / 1.349
  }

  if (
    !is.finite(m) ||
    m <= 0
  ) {

    m <- sd(
      x,
      na.rm = TRUE
    )
  }

  if (
    !is.finite(m) ||
    m <= 0
  ) {

    m <- 1
  }

  m
}

qc_thresholds <- list()
qc_keep <- rep(
  FALSE,
  ncol(obj)
)

names(
  qc_keep
) <- colnames(obj)

sample_order <- audit$sample_name

for (sample_name in sample_order) {

  idx <- which(
    obj$sample_name ==
    sample_name
  )

  md <- obj@meta.data[
    idx,
    ,
    drop = FALSE
  ]

  feature_med <- median(
    md$nFeature_RNA
  )

  feature_mad <- robust_scale(
    md$nFeature_RNA
  )

  count_med <- median(
    md$nCount_RNA
  )

  count_mad <- robust_scale(
    md$nCount_RNA
  )

  mt_med <- median(
    md$percent.mt
  )

  mt_mad <- robust_scale(
    md$percent.mt
  )

  feature_low <- max(
    200,
    floor(
      feature_med -
      5 * feature_mad
    )
  )

  feature_high <- ceiling(
    feature_med +
    5 * feature_mad
  )

  count_low <- max(
    500,
    floor(
      count_med -
      5 * count_mad
    )
  )

  count_high <- ceiling(
    count_med +
    5 * count_mad
  )

  mt_high <- min(
    40,
    max(
      20,
      mt_med +
      5 * mt_mad
    )
  )

  keep <- (
    md$nFeature_RNA >=
      feature_low &
    md$nFeature_RNA <=
      feature_high &
    md$nCount_RNA >=
      count_low &
    md$nCount_RNA <=
      count_high &
    md$percent.mt <=
      mt_high
  )

  qc_keep[
    rownames(md)
  ] <- keep

  retained <- sum(
    keep
  )

  raw_n <- nrow(
    md
  )

  retained_fraction <- (
    retained /
    raw_n
  )

  qc_thresholds[[
    sample_name
  ]] <- data.table(
    sample_name = sample_name,
    patient_id = unique(
      md$patient_id
    ),
    timepoint = unique(
      md$timepoint
    ),
    raw_cells = raw_n,
    retained_cells = retained,
    retained_fraction = retained_fraction,
    nFeature_low = feature_low,
    nFeature_high = feature_high,
    nCount_low = count_low,
    nCount_high = count_high,
    percent_mt_high = mt_high,
    nFeature_median = feature_med,
    nCount_median = count_med,
    percent_mt_median = mt_med
  )

  cat(
    sprintf(
      "%-10s %-4s %-6s raw=%5d keep=%5d (%5.1f%%)  features=%d-%d counts=%d-%d mt<=%.1f%%\n",
      sample_name,
      unique(md$patient_id),
      unique(md$timepoint),
      raw_n,
      retained,
      100 * retained_fraction,
      feature_low,
      feature_high,
      count_low,
      count_high,
      mt_high
    )
  )

  if (
    retained < 500
  ) {

    stop(
      "\nQC_GATE_STOP: fewer than 500 cells retained in ",
      sample_name
    )
  }

  if (
    retained_fraction < 0.40
  ) {

    stop(
      "\nQC_GATE_STOP: more than 60% of cells removed in ",
      sample_name,
      ". Review before proceeding."
    )
  }
}

qc_threshold_table <- rbindlist(
  qc_thresholds,
  use.names = TRUE
)

QC_THRESHOLD_OUT <- file.path(
  META_DIR,
  "STEP22A_05A_SAMPLE_QC_THRESHOLDS.csv"
)

fwrite(
  qc_threshold_table,
  QC_THRESHOLD_OUT,
  bom = TRUE
)

obj$STEP22A_QC_PASS <- qc_keep[
  colnames(obj)
]

# ============================================================
# 7. Save raw cell QC metadata before filtering
# ============================================================

raw_qc <- as.data.table(
  obj@meta.data,
  keep.rownames = "cell_id"
)

RAW_QC_OUT <- file.path(
  META_DIR,
  "STEP22A_05A_CELL_QC_METADATA_RAW.csv"
)

fwrite(
  raw_qc[
    ,
    .(
      cell_id,
      patient_id,
      timepoint,
      sample_name,
      hrs_accession,
      nCount_RNA,
      nFeature_RNA,
      percent.mt,
      STEP22A_QC_PASS
    )
  ],
  RAW_QC_OUT,
  bom = TRUE
)

# ============================================================
# 8. Apply technical QC
# ============================================================

cat("\n===== 7. APPLY TECHNICAL QC =====\n\n")

raw_cells <- ncol(
  obj
)

keep_cells <- colnames(
  obj
)[
  obj$STEP22A_QC_PASS
]

obj <- subset(
  obj,
  cells = keep_cells
)

retained_cells <- ncol(
  obj
)

cat(
  "Raw cells      :",
  raw_cells,
  "\n"
)

cat(
  "Retained cells :",
  retained_cells,
  "\n"
)

cat(
  "Removed cells  :",
  raw_cells -
  retained_cells,
  "\n"
)

cat(
  "Retained       :",
  sprintf(
    "%.1f%%",
    100 *
      retained_cells /
      raw_cells
  ),
  "\n\n"
)

before_cells <- sum(
  obj$timepoint ==
  "Before"
)

after_cells <- sum(
  obj$timepoint ==
  "After"
)

cat(
  "Retained Before cells :",
  before_cells,
  "\n"
)

cat(
  "Retained After cells  :",
  after_cells,
  "\n\n"
)

# ============================================================
# 9. Per-sample retained-cell summary
# ============================================================

retained_summary <- as.data.table(
  obj@meta.data
)[
  ,
  .(
    retained_cells = .N,
    median_nFeature = median(
      nFeature_RNA
    ),
    median_nCount = median(
      nCount_RNA
    ),
    median_percent_mt = median(
      percent.mt
    )
  ),
  by = .(
    patient_id,
    timepoint,
    sample_name,
    hrs_accession
  )
]

setorder(
  retained_summary,
  patient_id,
  timepoint
)

RETAINED_OUT <- file.path(
  META_DIR,
  "STEP22A_05A_RETAINED_CELL_COUNTS.csv"
)

fwrite(
  retained_summary,
  RETAINED_OUT,
  bom = TRUE
)

print(
  retained_summary
)

# ============================================================
# 10. Normalization + HVGs
# ============================================================

cat("\n===== 8. NORMALIZE DATA =====\n\n")

DefaultAssay(
  obj
) <- "RNA"

obj <- NormalizeData(
  obj,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

obj <- FindVariableFeatures(
  obj,
  selection.method = "vst",
  nfeatures = 3000,
  verbose = FALSE
)

cat(
  "Variable features:",
  length(
    VariableFeatures(
      obj
    )
  ),
  "\n\n"
)

# ============================================================
# 11. PCA
# ============================================================

cat("===== 9. PCA =====\n\n")

obj <- ScaleData(
  obj,
  features = VariableFeatures(
    obj
  ),
  verbose = FALSE
)

obj <- RunPCA(
  obj,
  features = VariableFeatures(
    obj
  ),
  npcs = 30,
  seed.use = 20260814,
  verbose = FALSE
)

cat(
  "PCA components:",
  ncol(
    Embeddings(
      obj,
      "pca"
    )
  ),
  "\n\n"
)

# ============================================================
# 12. Neighbors + clustering
#
# Broad-compartment discovery only.
# No treatment information enters clustering.
# ============================================================

cat("===== 10. GRAPH CLUSTERING =====\n\n")

obj <- FindNeighbors(
  obj,
  reduction = "pca",
  dims = 1:30,
  verbose = FALSE
)

obj <- FindClusters(
  obj,
  resolution = 0.4,
  random.seed = 20260814,
  verbose = FALSE
)

cluster_counts <- as.data.table(
  table(
    obj$seurat_clusters
  )
)

setnames(
  cluster_counts,
  c(
    "cluster",
    "cells"
  )
)

CLUSTER_COUNT_OUT <- file.path(
  META_DIR,
  "STEP22A_05A_PREANNOTATION_CLUSTER_COUNTS.csv"
)

fwrite(
  cluster_counts,
  CLUSTER_COUNT_OUT,
  bom = TRUE
)

cat(
  "Clusters detected:",
  nrow(
    cluster_counts
  ),
  "\n\n"
)

print(
  cluster_counts
)

# ============================================================
# 13. UMAP
# ============================================================

cat("\n===== 11. UMAP =====\n\n")

obj <- RunUMAP(
  obj,
  reduction = "pca",
  dims = 1:30,
  seed.use = 20260814,
  verbose = FALSE
)

# ============================================================
# 14. Save preannotation cell metadata
# ============================================================

cell_meta_out <- as.data.table(
  obj@meta.data,
  keep.rownames = "cell_id"
)

CELL_META_OUT <- file.path(
  META_DIR,
  "STEP22A_05A_PREANNOTATION_CELL_METADATA.csv"
)

fwrite(
  cell_meta_out[
    ,
    .(
      cell_id,
      patient_id,
      timepoint,
      sample_name,
      hrs_accession,
      nCount_RNA,
      nFeature_RNA,
      percent.mt,
      STEP22A_QC_PASS,
      seurat_clusters
    )
  ],
  CELL_META_OUT,
  bom = TRUE
)

# ============================================================
# 15. QC plots
# ============================================================

cat("===== 12. SAVE QC / PREANNOTATION FIGURES =====\n\n")

qc_plot_df <- as.data.frame(
  raw_qc
)

p_feature <- ggplot(
  qc_plot_df,
  aes(
    x = sample_name,
    y = nFeature_RNA,
    fill = STEP22A_QC_PASS
  )
) +
  geom_boxplot(
    outlier.size = 0.2
  ) +
  scale_y_log10() +
  labs(
    title = "Genes Detected per Cell Before QC",
    x = "Sample",
    y = "Detected Genes",
    fill = "QC Pass"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(
      angle = 60,
      hjust = 1
    )
  )

ggsave(
  filename = file.path(
    FIGURE_DIR,
    "STEP22A_05A_QC_NFEATURE_BY_SAMPLE.png"
  ),
  plot = p_feature,
  width = 12,
  height = 6,
  dpi = 180
)

p_mt <- ggplot(
  qc_plot_df,
  aes(
    x = sample_name,
    y = percent.mt,
    fill = STEP22A_QC_PASS
  )
) +
  geom_boxplot(
    outlier.size = 0.2
  ) +
  labs(
    title = "Mitochondrial Fraction Before QC",
    x = "Sample",
    y = "Mitochondrial Reads (%)",
    fill = "QC Pass"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(
      angle = 60,
      hjust = 1
    )
  )

ggsave(
  filename = file.path(
    FIGURE_DIR,
    "STEP22A_05A_QC_MITO_BY_SAMPLE.png"
  ),
  plot = p_mt,
  width = 12,
  height = 6,
  dpi = 180
)

p_umap_cluster <- DimPlot(
  obj,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE
) +
  labs(
    title = "Pre-annotation UMAP: Unsupervised Clusters"
  )

ggsave(
  filename = file.path(
    FIGURE_DIR,
    "STEP22A_05A_PREANNOTATION_UMAP_CLUSTERS.png"
  ),
  plot = p_umap_cluster,
  width = 10,
  height = 8,
  dpi = 180
)

p_umap_sample <- DimPlot(
  obj,
  reduction = "umap",
  group.by = "sample_name"
) +
  labs(
    title = "Pre-annotation UMAP by Sample"
  )

ggsave(
  filename = file.path(
    FIGURE_DIR,
    "STEP22A_05A_PREANNOTATION_UMAP_SAMPLE.png"
  ),
  plot = p_umap_sample,
  width = 10,
  height = 8,
  dpi = 180
)

cat(
  "✓ QC and preannotation plots saved.\n\n"
)

# ============================================================
# 16. Record explicit anti-circularity status
# ============================================================

guard_out <- file.path(
  META_DIR,
  "STEP22A_05A_ANTI_CIRCULARITY_STATUS.csv"
)

guard_table <- data.table(
  item = c(
    "cell_type_annotation_performed",
    "fibroblast_selection_performed",
    "CORE2_scoring_performed",
    "CORE4_scoring_performed",
    "CDKN1B_used_for_selection",
    "GSN_used_for_selection",
    "BGN_used_for_selection",
    "TIMP1_used_for_selection",
    "treatment_used_for_QC_thresholds",
    "treatment_used_for_clustering",
    "histology_inference_performed"
  ),
  status = c(
    "NO",
    "NO",
    "NO",
    "NO",
    "NO",
    "NO",
    "NO",
    "NO",
    "NO",
    "NO",
    "NO"
  )
)

fwrite(
  guard_table,
  guard_out,
  bom = TRUE
)

# ============================================================
# 17. Save Seurat object
# ============================================================

cat("===== 13. SAVE PREANNOTATION SEURAT OBJECT =====\n\n")

saveRDS(
  obj,
  OUT_OBJECT,
  compress = FALSE
)

object_size_gib <- file.info(
  OUT_OBJECT
)$size / 1024^3

cat(
  "Object saved:\n",
  OUT_OBJECT,
  "\n"
)

cat(
  "Object file size:",
  sprintf(
    "%.3f GiB",
    object_size_gib
  ),
  "\n\n"
)

# ============================================================
# 18. Frozen Step19-21 integrity check
# ============================================================

cat("===== 14. FROZEN PROJECT CHECK =====\n\n")

frozen_changes <- list()

if (
  file.exists(
    FROZEN_INVENTORY
  )
) {

  frozen <- fread(
    FROZEN_INVENTORY,
    encoding = "UTF-8"
  )

  for (
    i in seq_len(
      nrow(
        frozen
      )
    )
  ) {

    path <- frozen$path[i]

    if (!file.exists(path)) {

      frozen_changes[[
        length(
          frozen_changes
        ) + 1
      ]] <- c(
        path,
        "MISSING"
      )

      next
    }

    expected <- as.numeric(
      frozen$size_bytes[i]
    )

    actual <- file.info(
      path
    )$size

    if (
      is.finite(expected) &&
      expected != actual
    ) {

      frozen_changes[[
        length(
          frozen_changes
        ) + 1
      ]] <- c(
        path,
        "SIZE_CHANGED"
      )
    }
  }
}

cat(
  "Frozen original files changed:",
  length(
    frozen_changes
  ),
  "\n"
)

if (
  length(
    frozen_changes
  ) > 0
) {

  print(
    frozen_changes
  )

  stop(
    "STEP22A_FROZEN_PROJECT_CHANGED"
  )
}

cat(
  "✓ Frozen Step19-21 inventory unchanged.\n\n"
)

# ============================================================
# 19. Session info
# ============================================================

sink(
  file.path(
    LOG_DIR,
    "STEP22A_05A_sessionInfo.txt"
  )
)

print(
  sessionInfo()
)

sink()

# ============================================================
# 20. Summary
# ============================================================

summary_out <- file.path(
  META_DIR,
  "STEP22A_05A_QC_CLUSTERING_SUMMARY.txt"
)

summary_lines <- c(
  "STEP22A-5A QC + PREANNOTATION CLUSTERING",
  paste0(
    "Raw cells: ",
    raw_cells
  ),
  paste0(
    "Retained cells: ",
    retained_cells
  ),
  paste0(
    "Retained fraction: ",
    sprintf(
      "%.4f",
      retained_cells /
        raw_cells
    )
  ),
  paste0(
    "Before retained cells: ",
    before_cells
  ),
  paste0(
    "After retained cells: ",
    after_cells
  ),
  paste0(
    "Clusters: ",
    nrow(
      cluster_counts
    )
  ),
  "Cell type annotation performed: NO",
  "Fibroblast selection performed: NO",
  "CORE2 scoring performed: NO",
  "CORE4 scoring performed: NO",
  "Histology mapping: INCOMPLETE",
  "Primary CORE2 inference: BLOCKED",
  "Frozen Step19-21 changed: NO",
  "Manuscript changed: NO"
)

writeLines(
  summary_lines,
  summary_out
)

# ============================================================
# FINAL PRINT
# ============================================================

cat("\n")
cat("====================================================================\n")
cat("STEP22A-5A COMPLETE\n")
cat("====================================================================\n\n")

cat(
  "PAIRED PATIENTS              : 6\n"
)

cat(
  "PAIRED TUMOR SAMPLES         : 12\n\n"
)

cat(
  "RAW CELLS                    :",
  raw_cells,
  "\n"
)

cat(
  "QC-RETAINED CELLS            :",
  retained_cells,
  "\n"
)

cat(
  "RETAINED FRACTION            :",
  sprintf(
    "%.1f%%",
    100 *
      retained_cells /
      raw_cells
  ),
  "\n"
)

cat(
  "  BEFORE                     :",
  before_cells,
  "\n"
)

cat(
  "  AFTER                      :",
  after_cells,
  "\n\n"
)

cat(
  "PREANNOTATION CLUSTERS       :",
  nrow(
    cluster_counts
  ),
  "\n\n"
)

cat(
  "CELL TYPE ANNOTATION         : NOT PERFORMED\n"
)

cat(
  "FIBROBLAST SELECTION         : NOT PERFORMED\n"
)

cat(
  "CORE2 SCORING                : NOT PERFORMED\n"
)

cat(
  "CORE4 SCORING                : NOT PERFORMED\n\n"
)

cat(
  "CORE GENES USED FOR SELECTION: NO\n"
)

cat(
  "TREATMENT USED FOR QC        : NO\n"
)

cat(
  "TREATMENT USED FOR CLUSTERING: NO\n\n"
)

cat(
  "HISTOLOGY MAPPING            : INCOMPLETE\n"
)

cat(
  "PRIMARY CORE2 INFERENCE      : BLOCKED\n\n"
)

cat(
  "FROZEN STEP19-21 MODIFIED    : NO\n"
)

cat(
  "MANUSCRIPT MODIFIED          : NO\n\n"
)

cat(
  "SEURAT OBJECT:\n"
)

cat(
  "  03_objects/OMIX005710/\n"
)

cat(
  "  OMIX005710_PAIRED_QC_CLUSTERED_PREANNOTATION.rds\n\n"
)

cat(
  "NEXT REQUIRED STEP:\n"
)

cat(
  "STEP22A-5B = CONSTRAINED BROAD CELL-TYPE ANNOTATION\n"
)

cat(
  "USING ONLY THE FROZEN NON-CORE MARKER SET\n\n"
)

cat(
  "STOP HERE.\n"
)

cat(
  "DO NOT SCORE CORE2 YET.\n\n"
)
