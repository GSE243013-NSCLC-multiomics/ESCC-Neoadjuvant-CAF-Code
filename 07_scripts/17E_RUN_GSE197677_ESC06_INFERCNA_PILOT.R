# ============================================================
# ESCC Neoadjuvant Project
# STEP 17E
#
# FIRST REAL InferCNA pilot
#
# Dataset:
#   GSE197677
#
# Tumor:
#   ESC06 epithelial cells
#
# Normal reference:
#   ESN11 / ESN12 / ESN13 / ESN15 epithelial cells
#
# IMPORTANT:
#   THIS STEP DOES RUN CNA INFERENCE
#
# BUT:
#   NO malignant labels
#   NO findMalignant()
#   NO cell filtering
#   NO batch processing
#   NO other tumor samples
# ============================================================

options(
  stringsAsFactors = FALSE,
  timeout = 3600
)

suppressPackageStartupMessages({
  library(infercna)
  library(Matrix)
  library(data.table)
  library(ggplot2)
})

cat("\n============================================================\n")
cat("STEP 17E: GSE197677 ESC06 INFERCNA PILOT\n")
cat("============================================================\n\n")

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------

indir <- paste0(
  "06_tables/cross_dataset/",
  "malignant_epithelial/",
  "GSE197677_ESC06_pilot"
)

outdir <- paste0(
  "06_tables/cross_dataset/",
  "malignant_epithelial/",
  "GSE197677_ESC06_pilot/",
  "17E_results"
)

figdir <- paste0(
  "05_figures/GSE197677/",
  "malignant_epithelial/",
  "ESC06_pilot"
)

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  figdir,
  recursive = TRUE,
  showWarnings = FALSE
)

cpm_file <- file.path(
  indir,
  "GSE197677_ESC06_infercna_CPM_hg38.rds"
)

refs_file <- file.path(
  indir,
  "GSE197677_ESC06_infercna_refCells.rds"
)

annotation_file <- file.path(
  indir,
  "GSE197677_ESC06_infercna_cell_annotation.csv"
)

gene_file <- file.path(
  indir,
  "GSE197677_ESC06_infercna_gene_table.csv"
)

required <- c(
  cpm_file,
  refs_file,
  annotation_file,
  gene_file
)

missing <- required[
  !file.exists(required)
]

if (length(missing) > 0) {

  stop(
    "Missing Step 17D files:\n",
    paste(
      missing,
      collapse = "\n"
    )
  )
}

# ============================================================
# 1. LOAD PILOT INPUTS
# ============================================================

cat("===== 1. LOAD STEP 17D INPUTS =====\n\n")

cpm <- readRDS(
  cpm_file
)

refCells <- readRDS(
  refs_file
)

annotation <- fread(
  annotation_file
)

gene_table <- fread(
  gene_file
)

cat(
  "CPM matrix: ",
  nrow(cpm),
  " genes x ",
  ncol(cpm),
  " cells\n",
  sep = ""
)

cat(
  "Reference groups: ",
  length(refCells),
  "\n",
  sep = ""
)

print(
  lengths(refCells)
)

cat(
  "Cell annotation rows: ",
  nrow(annotation),
  "\n",
  sep = ""
)

# ============================================================
# 2. VALIDATE PILOT
# ============================================================

cat("\n===== 2. PILOT VALIDATION =====\n\n")

if (ncol(cpm) != 392) {
  stop(
    "Expected 392 pilot cells; found ",
    ncol(cpm)
  )
}

tumor_cells <- annotation[
  infercna_role == "TUMOR_TEST",
  cell_id
]

normal_cells <- unlist(
  refCells,
  use.names = FALSE
)

cat(
  "Tumor cells : ",
  length(tumor_cells),
  "\n",
  sep = ""
)

cat(
  "Normal cells: ",
  length(normal_cells),
  "\n",
  sep = ""
)

if (length(tumor_cells) != 81) {
  stop(
    "Expected 81 ESC06 tumor epithelial cells."
  )
}

if (length(normal_cells) != 311) {
  stop(
    "Expected 311 normal epithelial cells."
  )
}

if (
  !all(
    tumor_cells %in%
      colnames(cpm)
  )
) {
  stop(
    "Tumor cells missing from CPM matrix."
  )
}

if (
  !all(
    normal_cells %in%
      colnames(cpm)
  )
) {
  stop(
    "Normal reference cells missing from CPM matrix."
  )
}

if (
  length(
    intersect(
      tumor_cells,
      normal_cells
    )
  ) > 0
) {
  stop(
    "Tumor/reference cell overlap detected."
  )
}

cat("✓ Cell identities validated\n")

# ============================================================
# 3. CONFIGURE hg38
# ============================================================

cat("\n===== 3. CONFIGURE INFERCNA hg38 =====\n\n")

infercna::useGenome(
  "hg38"
)

genome <- infercna::retrieveGenome(
  name = "hg38"
)

cat(
  "Genome reference genes: ",
  nrow(genome),
  "\n",
  sep = ""
)

cat(
  "Input genes           : ",
  nrow(cpm),
  "\n",
  sep = ""
)

# ============================================================
# 4. CONVERT PILOT TO DENSE MATRIX
# ============================================================

cat("\n===== 4. PREPARE INFERCNA MATRIX =====\n\n")

dense_estimate_gb <- (
  nrow(cpm) *
  ncol(cpm) *
  8
) / 1024^3

cat(
  sprintf(
    "Dense matrix estimate: %.3f GB\n",
    dense_estimate_gb
  )
)

m <- as.matrix(
  cpm
)

storage.mode(m) <- "double"

cat(
  "Dense matrix: ",
  nrow(m),
  " x ",
  ncol(m),
  "\n",
  sep = ""
)

cat(
  "Range: ",
  sprintf(
    "%.3f",
    min(m)
  ),
  " to ",
  sprintf(
    "%.3f",
    max(m)
  ),
  "\n",
  sep = ""
)

if (
  any(
    !is.finite(m)
  )
) {
  stop(
    "CPM matrix contains NA / Inf."
  )
}

# ============================================================
# 5. PARAMETER FREEZE
# ============================================================

params <- data.table(

  parameter = c(
    "genome",
    "window",
    "range_lower",
    "range_upper",
    "n",
    "noise",
    "center.method",
    "isLog"
  ),

  value = c(
    "hg38",
    "100",
    "-3",
    "3",
    "5000",
    "0.1",
    "median",
    "FALSE"
  )
)

cat("\n===== 5. PILOT PARAMETERS =====\n\n")

print(params)

fwrite(
  params,
  file.path(
    outdir,
    "GSE197677_ESC06_17E_parameters.csv"
  ),
  bom = TRUE
)

# ============================================================
# 6. RUN INFERCNA
# ============================================================

cat("\n============================================================\n")
cat("6. RUNNING INFERCNA — FIRST REAL CNV PILOT\n")
cat("============================================================\n\n")

start_time <- Sys.time()

cna <- tryCatch(

  infercna::infercna(
    m = m,
    refCells = refCells,
    window = 100,
    range = c(-3, 3),
    n = 5000,
    noise = 0.1,
    center.method = "median",
    isLog = FALSE,
    verbose = TRUE
  ),

  error = function(e) {

    cat(
      "\nINFERCNA ERROR:\n",
      conditionMessage(e),
      "\n"
    )

    stop(
      "InferCNA pilot failed."
    )
  }
)

end_time <- Sys.time()

runtime_minutes <- as.numeric(
  difftime(
    end_time,
    start_time,
    units = "mins"
  )
)

cat(
  "\nRuntime: ",
  sprintf(
    "%.2f minutes",
    runtime_minutes
  ),
  "\n",
  sep = ""
)

# ============================================================
# 7. CNA MATRIX VALIDATION
# ============================================================

cat("\n============================================================\n")
cat("7. CNA MATRIX VALIDATION\n")
cat("============================================================\n\n")

cat(
  "CNA matrix: ",
  nrow(cna),
  " genes x ",
  ncol(cna),
  " cells\n",
  sep = ""
)

cat(
  "CNA range: ",
  sprintf(
    "%.4f",
    min(cna, na.rm = TRUE)
  ),
  " to ",
  sprintf(
    "%.4f",
    max(cna, na.rm = TRUE)
  ),
  "\n",
  sep = ""
)

cat(
  "NA values : ",
  sum(is.na(cna)),
  "\n",
  sep = ""
)

cat(
  "Inf values: ",
  sum(is.infinite(cna)),
  "\n",
  sep = ""
)

if (
  ncol(cna) != 392
) {
  stop(
    "Unexpected number of CNA columns."
  )
}

if (
  any(
    !is.finite(cna)
  )
) {
  stop(
    "CNA matrix contains non-finite values."
  )
}

if (
  !all(
    tumor_cells %in%
      colnames(cna)
  )
) {
  stop(
    "Tumor cells missing from CNA matrix."
  )
}

if (
  !all(
    normal_cells %in%
      colnames(cna)
  )
) {
  stop(
    "Reference cells missing from CNA matrix."
  )
}

cat("✓ CNA matrix validated\n")

# ============================================================
# 8. SAVE CNA MATRIX IMMEDIATELY
# ============================================================

cat("\n===== 8. SAVE CNA MATRICES =====\n\n")

cna_all_file <- file.path(
  outdir,
  "GSE197677_ESC06_infercna_CNA_all_cells.rds"
)

cna_tumor_file <- file.path(
  outdir,
  "GSE197677_ESC06_infercna_CNA_tumor_only.rds"
)

saveRDS(
  cna,
  cna_all_file,
  compress = FALSE
)

cna_tumor <- cna[
  ,
  tumor_cells,
  drop = FALSE
]

saveRDS(
  cna_tumor,
  cna_tumor_file,
  compress = FALSE
)

cat(
  "✓ ",
  cna_all_file,
  "\n",
  sep = ""
)

cat(
  "✓ ",
  cna_tumor_file,
  "\n",
  sep = ""
)

# ============================================================
# 9. CNA SIGNAL
# ============================================================

cat("\n============================================================\n")
cat("9. CNA SIGNAL DIAGNOSTICS\n")
cat("============================================================\n\n")

cna_signal <- infercna::cnaSignal(
  cna
)

if (
  is.null(
    names(cna_signal)
  )
) {
  names(cna_signal) <-
    colnames(cna)
}

mean_abs_cna <- colMeans(
  abs(cna)
)

fraction_nonzero <- colMeans(
  abs(cna) > 1e-8
)

cell_metrics <- data.table(

  cell_id =
    colnames(cna),

  cna_signal =
    as.numeric(
      cna_signal[
        colnames(cna)
      ]
    ),

  mean_abs_cna =
    as.numeric(
      mean_abs_cna
    ),

  fraction_nonzero =
    as.numeric(
      fraction_nonzero
    )
)

cell_metrics <- merge(
  cell_metrics,
  annotation[
    ,
    .(
      cell_id,
      sample,
      tissue,
      treatment_group,
      infercna_role,
      reference_group
    )
  ],
  by = "cell_id",
  all.x = TRUE
)

if (
  any(
    is.na(
      cell_metrics$infercna_role
    )
  )
) {
  stop(
    "Some CNA cells lack annotation."
  )
}

setorder(
  cell_metrics,
  infercna_role,
  sample,
  -cna_signal
)

fwrite(
  cell_metrics,
  file.path(
    outdir,
    "GSE197677_ESC06_17E_cell_CNA_metrics.csv"
  ),
  bom = TRUE
)

# ============================================================
# 10. QC SUMMARY BY ROLE / SAMPLE
# ============================================================

cat("\n============================================================\n")
cat("10. CNA QC SUMMARY BY ROLE / SAMPLE\n")
cat("============================================================\n\n")

qc_by_sample <- cell_metrics[
  ,
  .(
    cells = .N,

    median_CNA_signal =
      median(
        cna_signal,
        na.rm = TRUE
      ),

    mean_CNA_signal =
      mean(
        cna_signal,
        na.rm = TRUE
      ),

    median_mean_abs_CNA =
      median(
        mean_abs_cna,
        na.rm = TRUE
      ),

    median_fraction_nonzero =
      median(
        fraction_nonzero,
        na.rm = TRUE
      ),

    p25_CNA_signal =
      as.numeric(
        quantile(
          cna_signal,
          0.25,
          na.rm = TRUE
        )
      ),

    p75_CNA_signal =
      as.numeric(
        quantile(
          cna_signal,
          0.75,
          na.rm = TRUE
        )
      )
  ),
  by = .(
    infercna_role,
    sample
  )
]

setorder(
  qc_by_sample,
  infercna_role,
  sample
)

print(qc_by_sample)

fwrite(
  qc_by_sample,
  file.path(
    outdir,
    "GSE197677_ESC06_17E_QC_by_sample.csv"
  ),
  bom = TRUE
)

# ============================================================
# 11. GLOBAL TUMOR VS NORMAL QC
# ============================================================

cat("\n============================================================\n")
cat("11. TUMOR VS NORMAL CNA SIGNAL SUMMARY\n")
cat("============================================================\n\n")

qc_by_role <- cell_metrics[
  ,
  .(
    cells = .N,

    median_CNA_signal =
      median(
        cna_signal,
        na.rm = TRUE
      ),

    mean_CNA_signal =
      mean(
        cna_signal,
        na.rm = TRUE
      ),

    median_mean_abs_CNA =
      median(
        mean_abs_cna,
        na.rm = TRUE
      ),

    median_fraction_nonzero =
      median(
        fraction_nonzero,
        na.rm = TRUE
      ),

    min_CNA_signal =
      min(
        cna_signal,
        na.rm = TRUE
      ),

    max_CNA_signal =
      max(
        cna_signal,
        na.rm = TRUE
      )
  ),
  by = infercna_role
]

print(qc_by_role)

fwrite(
  qc_by_role,
  file.path(
    outdir,
    "GSE197677_ESC06_17E_QC_by_role.csv"
  ),
  bom = TRUE
)

# ============================================================
# 12. ESC06 TUMOR CNA SIGNAL DISTRIBUTION
# ============================================================

cat("\n============================================================\n")
cat("12. ESC06 TUMOR CNA SIGNAL QUANTILES\n")
cat("============================================================\n\n")

tumor_metric <- cell_metrics[
  infercna_role ==
    "TUMOR_TEST"
]

tumor_quantiles <- data.table(

  quantile = c(
    "min",
    "p05",
    "p25",
    "median",
    "p75",
    "p95",
    "max"
  ),

  CNA_signal = as.numeric(
    quantile(
      tumor_metric$cna_signal,
      probs = c(
        0,
        0.05,
        0.25,
        0.50,
        0.75,
        0.95,
        1
      ),
      na.rm = TRUE
    )
  )
)

print(tumor_quantiles)

fwrite(
  tumor_quantiles,
  file.path(
    outdir,
    "GSE197677_ESC06_17E_tumor_CNA_signal_quantiles.csv"
  ),
  bom = TRUE
)

# ============================================================
# 13. ORDER CELLS FOR HEATMAP
# ============================================================

cat("\n===== 13. PREPARE HEATMAP ORDER =====\n\n")

ordered_normal_cells <- unlist(
  lapply(
    names(refCells),
    function(x) {
      refCells[[x]]
    }
  ),
  use.names = FALSE
)

ordered_cells <- c(
  ordered_normal_cells,
  tumor_cells
)

ordered_cells <- ordered_cells[
  ordered_cells %in%
    colnames(cna)
]

cna_ordered <- cna[
  ,
  ordered_cells,
  drop = FALSE
]

cat(
  "Heatmap cells: ",
  ncol(cna_ordered),
  "\n",
  sep = ""
)

# ============================================================
# 14. FULL PILOT CNA HEATMAP
# ============================================================

cat("\n===== 14. FULL PILOT CNA HEATMAP =====\n\n")

p_all <- tryCatch(

  infercna::cnaPlot(
    cna = cna_ordered,
    limits = c(-1, 1),
    x.hide = NULL,
    orderCells = FALSE,
    title =
      "GSE197677 ESC06 InferCNA Pilot",
    subtitle =
      "4 normal epithelial reference groups followed by ESC06 tumor epithelial cells"
  ),

  error = function(e) {

    cat(
      "Full heatmap plotting failed:\n",
      conditionMessage(e),
      "\n"
    )

    NULL
  }
)

if (!is.null(p_all)) {

  ggsave(
    file.path(
      figdir,
      "GSE197677_ESC06_InferCNA_all_cells_heatmap.png"
    ),
    p_all,
    width = 13,
    height = 9,
    dpi = 300
  )

  ggsave(
    file.path(
      figdir,
      "GSE197677_ESC06_InferCNA_all_cells_heatmap.pdf"
    ),
    p_all,
    width = 13,
    height = 9
  )

  cat("✓ Full CNA heatmap saved\n")
}

# ============================================================
# 15. ESC06 TUMOR-ONLY HEATMAP
# ============================================================

cat("\n===== 15. ESC06 TUMOR-ONLY CNA HEATMAP =====\n\n")

p_tumor <- tryCatch(

  infercna::cnaPlot(
    cna = cna_tumor,
    limits = c(-1, 1),
    x.hide = NULL,
    orderCells = TRUE,
    title =
      "GSE197677 ESC06 Tumor Epithelial CNA",
    subtitle =
      "Tumor epithelial cells only; hierarchical cell ordering"
  ),

  error = function(e) {

    cat(
      "Tumor heatmap plotting failed:\n",
      conditionMessage(e),
      "\n"
    )

    NULL
  }
)

if (!is.null(p_tumor)) {

  ggsave(
    file.path(
      figdir,
      "GSE197677_ESC06_InferCNA_tumor_only_heatmap.png"
    ),
    p_tumor,
    width = 13,
    height = 7,
    dpi = 300
  )

  ggsave(
    file.path(
      figdir,
      "GSE197677_ESC06_InferCNA_tumor_only_heatmap.pdf"
    ),
    p_tumor,
    width = 13,
    height = 7
  )

  cat("✓ Tumor-only CNA heatmap saved\n")
}

# ============================================================
# 16. CNA SIGNAL DISTRIBUTION PLOT
# ============================================================

cat("\n===== 16. CNA SIGNAL DISTRIBUTION PLOT =====\n\n")

plot_dt <- copy(
  cell_metrics
)

plot_dt[
  ,
  display_group :=
    fifelse(
      infercna_role ==
        "TUMOR_TEST",
      "ESC06 tumor",
      sample
    )
]

plot_dt[
  ,
  display_group :=
    factor(
      display_group,
      levels = c(
        "ESN11",
        "ESN12",
        "ESN13",
        "ESN15",
        "ESC06 tumor"
      )
    )
]

p_signal <- ggplot(
  plot_dt,
  aes(
    x = display_group,
    y = cna_signal
  )
) +
  geom_boxplot(
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.18,
    height = 0,
    size = 0.8,
    alpha = 0.45
  ) +
  labs(
    title =
      "GSE197677 ESC06 InferCNA Pilot",
    subtitle =
      "CNA signal by normal reference sample and tumor epithelial cells",
    x = "Cell Group",
    y = "CNA Signal"
  ) +
  theme_classic() +
  theme(
    axis.text.x =
      element_text(
        angle = 45,
        hjust = 1
      )
  )

ggsave(
  file.path(
    figdir,
    "GSE197677_ESC06_CNA_signal_by_group.png"
  ),
  p_signal,
  width = 8,
  height = 6,
  dpi = 300
)

ggsave(
  file.path(
    figdir,
    "GSE197677_ESC06_CNA_signal_by_group.pdf"
  ),
  p_signal,
  width = 8,
  height = 6
)

# ============================================================
# 17. PARAMETER / OUTPUT MANIFEST
# ============================================================

manifest <- data.table(

  dataset =
    "GSE197677",

  tumor_sample =
    "ESC06",

  tumor_cells =
    length(tumor_cells),

  normal_reference_groups =
    length(refCells),

  normal_reference_cells =
    length(normal_cells),

  input_genes =
    nrow(m),

  CNA_genes =
    nrow(cna),

  total_cells =
    ncol(cna),

  window =
    100,

  n_top_genes =
    5000,

  noise =
    0.1,

  center_method =
    "median",

  input_normalization =
    "CPM",

  isLog =
    FALSE,

  runtime_minutes =
    runtime_minutes,

  malignant_calling =
    FALSE,

  pilot_status =
    "CNA_GENERATED_REVIEW_REQUIRED"
)

fwrite(
  manifest,
  file.path(
    outdir,
    "GSE197677_ESC06_17E_manifest.csv"
  ),
  bom = TRUE
)

# ============================================================
# 18. RELOAD CNA VALIDATION
# ============================================================

cat("\n===== 18. RELOAD CNA VALIDATION =====\n\n")

check <- readRDS(
  cna_all_file
)

cat(
  "Reloaded CNA: ",
  nrow(check),
  " x ",
  ncol(check),
  "\n",
  sep = ""
)

if (
  !identical(
    dim(check),
    dim(cna)
  )
) {
  stop(
    "Reloaded CNA dimensions changed."
  )
}

cat("✓ CNA file reload validated\n")

# ============================================================
# FINAL
# ============================================================

cat("\n============================================================\n")
cat("STEP 17E FINISHED ✓\n")
cat("FIRST GSE197677 INFERCNA PILOT COMPLETE ✓\n")
cat("============================================================\n\n")

cat(
  "Pilot tumor        : ESC06\n"
)

cat(
  "Tumor cells        : ",
  length(tumor_cells),
  "\n",
  sep = ""
)

cat(
  "Normal refs        : ",
  length(refCells),
  " groups / ",
  length(normal_cells),
  " cells\n",
  sep = ""
)

cat(
  "CNA matrix         : ",
  nrow(cna),
  " genes x ",
  ncol(cna),
  " cells\n",
  sep = ""
)

cat(
  "Runtime            : ",
  sprintf(
    "%.2f min",
    runtime_minutes
  ),
  "\n",
  sep = ""
)

cat("\nIMPORTANT:\n")
cat("CNA profiles generated ✓\n")
cat("Malignant labels NOT assigned ✓\n")
cat("findMalignant() NOT run ✓\n")
cat("No tumor cells removed ✓\n")
cat("No other samples processed ✓\n")

cat("\nMain figures:\n")

cat(
  "  ",
  figdir,
  "/GSE197677_ESC06_InferCNA_all_cells_heatmap.png\n",
  sep = ""
)

cat(
  "  ",
  figdir,
  "/GSE197677_ESC06_InferCNA_tumor_only_heatmap.png\n",
  sep = ""
)

cat(
  "  ",
  figdir,
  "/GSE197677_ESC06_CNA_signal_by_group.png\n",
  sep = ""
)

cat("\nSTOP HERE.\n")
cat("DO NOT RUN findMalignant().\n")
cat("DO NOT PROCESS THE OTHER 9 TUMOR SAMPLES YET.\n")
