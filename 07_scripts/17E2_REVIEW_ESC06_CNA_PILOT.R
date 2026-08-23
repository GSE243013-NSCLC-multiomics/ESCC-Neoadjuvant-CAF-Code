# ============================================================
# ESCC Neoadjuvant Project
# STEP 17E2
#
# Review existing ESC06 InferCNA pilot.
#
# DOES NOT rerun infercna().
# DOES NOT change parameters.
# DOES NOT call malignant cells.
#
# Outputs:
#   - manual chromosome-ordered CNA heatmap
#   - tumor-only CNA heatmap
#   - chromosome-arm tumor vs normal summary
#   - tumor CNA correlation metrics
# ============================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(infercna)
})

cat("\n============================================================\n")
cat("STEP 17E2: REVIEW ESC06 CNA PILOT\n")
cat("============================================================\n\n")

indir <- paste0(
  "06_tables/cross_dataset/",
  "malignant_epithelial/",
  "GSE197677_ESC06_pilot"
)

resultdir <- file.path(
  indir,
  "17E_results"
)

figdir <- paste0(
  "05_figures/GSE197677/",
  "malignant_epithelial/",
  "ESC06_pilot/"
)

dir.create(
  figdir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# 1. Input files
# ------------------------------------------------------------

cna_file <- file.path(
  resultdir,
  "GSE197677_ESC06_infercna_CNA_all_cells.rds"
)

cna_tumor_file <- file.path(
  resultdir,
  "GSE197677_ESC06_infercna_CNA_tumor_only.rds"
)

annotation_file <- file.path(
  indir,
  "GSE197677_ESC06_infercna_cell_annotation.csv"
)

gene_file <- file.path(
  indir,
  "GSE197677_ESC06_infercna_gene_table.csv"
)

files <- c(
  cna_file,
  cna_tumor_file,
  annotation_file,
  gene_file
)

missing <- files[
  !file.exists(files)
]

if (length(missing) > 0) {
  stop(
    "Missing files:\n",
    paste(missing, collapse = "\n")
  )
}

# ------------------------------------------------------------
# 2. Load
# ------------------------------------------------------------

cat("===== 1. LOAD EXISTING CNA OUTPUT =====\n\n")

cna <- readRDS(cna_file)

cna_tumor <- readRDS(
  cna_tumor_file
)

ann <- fread(
  annotation_file
)

genes <- fread(
  gene_file
)

cat(
  "All-cell CNA matrix: ",
  nrow(cna),
  " genes x ",
  ncol(cna),
  " cells\n",
  sep = ""
)

cat(
  "Tumor CNA matrix   : ",
  nrow(cna_tumor),
  " genes x ",
  ncol(cna_tumor),
  " cells\n",
  sep = ""
)

# ------------------------------------------------------------
# 3. Match gene coordinates
# ------------------------------------------------------------

cat("\n===== 2. MATCH CNA GENES TO GENOME =====\n\n")

gene_order <- genes[
  gene %in% rownames(cna)
]

gene_order <- gene_order[
  !is.na(chromosome) &
  chromosome %in%
    c(
      as.character(1:22),
      "X"
    )
]

gene_order[
  ,
  chr_num :=
    fifelse(
      chromosome == "X",
      23L,
      as.integer(chromosome)
    )
]

setorder(
  gene_order,
  chr_num,
  start_position
)

gene_order <- gene_order[
  !duplicated(gene)
]

ordered_genes <- gene_order$gene

cna <- cna[
  ordered_genes,
  ,
  drop = FALSE
]

cna_tumor <- cna_tumor[
  ordered_genes[
    ordered_genes %in%
      rownames(cna_tumor)
  ],
  ,
  drop = FALSE
]

cat(
  "Genome-ordered genes: ",
  nrow(cna),
  "\n",
  sep = ""
)

# ------------------------------------------------------------
# 4. Cell order
# ------------------------------------------------------------

normal_order <- c(
  "ESN11",
  "ESN12",
  "ESN13",
  "ESN15"
)

normal_cells <- unlist(
  lapply(
    normal_order,
    function(s) {
      ann[
        sample == s &
        infercna_role ==
          "NORMAL_REFERENCE",
        cell_id
      ]
    }
  ),
  use.names = FALSE
)

tumor_cells <- ann[
  infercna_role ==
    "TUMOR_TEST",
  cell_id
]

cell_order <- c(
  normal_cells,
  tumor_cells
)

cell_order <- cell_order[
  cell_order %in%
    colnames(cna)
]

cna <- cna[
  ,
  cell_order,
  drop = FALSE
]

cat(
  "Ordered cells: ",
  length(cell_order),
  "\n",
  sep = ""
)

# ------------------------------------------------------------
# 5. Build chromosome positions
# ------------------------------------------------------------

chr_info <- gene_order[
  gene %in% rownames(cna)
]

chr_info[
  ,
  matrix_row :=
    match(
      gene,
      rownames(cna)
    )
]

chr_bounds <- chr_info[
  ,
  .(
    start = min(matrix_row),
    end = max(matrix_row),
    midpoint =
      (
        min(matrix_row) +
        max(matrix_row)
      ) / 2
  ),
  by = chromosome
]

chr_bounds[
  ,
  chr_num :=
    fifelse(
      chromosome == "X",
      23L,
      as.integer(chromosome)
    )
]

setorder(
  chr_bounds,
  chr_num
)

# ------------------------------------------------------------
# 6. Manual all-cell heatmap
# ------------------------------------------------------------

cat("\n===== 3. DRAW MANUAL ALL-CELL CNA HEATMAP =====\n\n")

png_file <- file.path(
  figdir,
  "GSE197677_ESC06_CNA_manual_all_cells_heatmap.png"
)

png(
  png_file,
  width = 4200,
  height = 2600,
  res = 300
)

par(
  mar = c(5, 6, 4, 2)
)

z <- t(
  cna
)

plot_limit <- 0.6

z_plot <- pmax(
  pmin(
    z,
    plot_limit
  ),
  -plot_limit
)

image(
  x = seq_len(ncol(z_plot)),
  y = seq_len(nrow(z_plot)),
  z = t(z_plot),
  axes = FALSE,
  xlab = "Chromosome",
  ylab = "Cells",
  main = "GSE197677 ESC06 InferCNA Pilot",
  col = hcl.colors(
    101,
    palette = "Blue-Red 3"
  )
)

for (b in chr_bounds$end) {

  abline(
    v = b + 0.5,
    col = "grey70",
    lwd = 0.5
  )
}

axis(
  side = 1,
  at = chr_bounds$midpoint,
  labels = chr_bounds$chromosome,
  las = 2,
  cex.axis = 0.7
)

cell_group <- ann[
  match(
    cell_order,
    cell_id
  ),
  fifelse(
    infercna_role ==
      "TUMOR_TEST",
    "ESC06",
    sample
  )
]

group_run <- rle(
  cell_group
)

group_end <- cumsum(
  group_run$lengths
)

group_start <- c(
  1,
  head(
    group_end,
    -1
  ) + 1
)

group_mid <- (
  group_start +
  group_end
) / 2

for (b in head(group_end, -1)) {

  abline(
    h = b + 0.5,
    col = "black",
    lwd = 1
  )
}

axis(
  side = 2,
  at = group_mid,
  labels = group_run$values,
  las = 2,
  cex.axis = 0.8
)

dev.off()

cat(
  "✓ ",
  png_file,
  "\n",
  sep = ""
)

# ------------------------------------------------------------
# 7. Tumor-only heatmap
# ------------------------------------------------------------

cat("\n===== 4. DRAW ESC06 TUMOR-ONLY HEATMAP =====\n\n")

tumor_avg <- rowMeans(
  cna_tumor
)

tumor_cor <- apply(
  cna_tumor,
  2,
  function(x) {
    suppressWarnings(
      cor(
        x,
        tumor_avg,
        method = "pearson",
        use = "pairwise.complete.obs"
      )
    )
  }
)

tumor_order <- names(
  sort(
    tumor_cor,
    decreasing = TRUE
  )
)

tumor_plot <- cna_tumor[
  ,
  tumor_order,
  drop = FALSE
]

png_file_tumor <- file.path(
  figdir,
  "GSE197677_ESC06_CNA_manual_tumor_only_heatmap.png"
)

png(
  png_file_tumor,
  width = 4200,
  height = 1800,
  res = 300
)

par(
  mar = c(5, 5, 4, 2)
)

zt <- t(
  tumor_plot
)

zt <- pmax(
  pmin(
    zt,
    plot_limit
  ),
  -plot_limit
)

image(
  x = seq_len(ncol(zt)),
  y = seq_len(nrow(zt)),
  z = t(zt),
  axes = FALSE,
  xlab = "Chromosome",
  ylab = "ESC06 tumor epithelial cells",
  main =
    "ESC06 Tumor Epithelial CNA Profiles",
  col = hcl.colors(
    101,
    palette = "Blue-Red 3"
  )
)

for (b in chr_bounds$end) {

  abline(
    v = b + 0.5,
    col = "grey70",
    lwd = 0.5
  )
}

axis(
  side = 1,
  at = chr_bounds$midpoint,
  labels = chr_bounds$chromosome,
  las = 2,
  cex.axis = 0.7
)

dev.off()

cat(
  "✓ ",
  png_file_tumor,
  "\n",
  sep = ""
)

# ------------------------------------------------------------
# 8. CNA signal + CNA correlation
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("5. CNA SIGNAL + CNA CORRELATION\n")
cat("============================================================\n\n")

signal <- infercna::cnaSignal(
  cna
)

if (
  is.null(
    names(signal)
  )
) {
  names(signal) <- colnames(cna)
}

tumor_cor_all <- rep(
  NA_real_,
  ncol(cna)
)

names(tumor_cor_all) <-
  colnames(cna)

tumor_mean_profile <- rowMeans(
  cna[
    ,
    tumor_cells,
    drop = FALSE
  ]
)

for (cell in colnames(cna)) {

  tumor_cor_all[cell] <- suppressWarnings(
    cor(
      cna[, cell],
      tumor_mean_profile,
      method = "pearson",
      use = "pairwise.complete.obs"
    )
  )
}

metrics <- data.table(
  cell_id =
    colnames(cna),

  cna_signal =
    as.numeric(
      signal[
        colnames(cna)
      ]
    ),

  tumor_profile_correlation =
    as.numeric(
      tumor_cor_all[
        colnames(cna)
      ]
    ),

  mean_abs_CNA =
    colMeans(
      abs(cna)
    )
)

metrics <- merge(
  metrics,
  ann[
    ,
    .(
      cell_id,
      sample,
      infercna_role
    )
  ],
  by = "cell_id",
  all.x = TRUE
)

metrics[
  ,
  display_group :=
    fifelse(
      infercna_role ==
        "TUMOR_TEST",
      "ESC06",
      sample
    )
]

metric_summary <- metrics[
  ,
  .(
    cells = .N,

    median_CNA_signal =
      median(
        cna_signal,
        na.rm = TRUE
      ),

    median_tumor_correlation =
      median(
        tumor_profile_correlation,
        na.rm = TRUE
      ),

    mean_tumor_correlation =
      mean(
        tumor_profile_correlation,
        na.rm = TRUE
      ),

    median_mean_abs_CNA =
      median(
        mean_abs_CNA,
        na.rm = TRUE
      )
  ),
  by = display_group
]

setorder(
  metric_summary,
  display_group
)

print(metric_summary)

fwrite(
  metrics,
  file.path(
    resultdir,
    "GSE197677_ESC06_17E2_cell_CNA_review_metrics.csv"
  ),
  bom = TRUE
)

fwrite(
  metric_summary,
  file.path(
    resultdir,
    "GSE197677_ESC06_17E2_group_CNA_review_summary.csv"
  ),
  bom = TRUE
)

# ------------------------------------------------------------
# 9. Chromosome-arm summary
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("6. CHROMOSOME-ARM TUMOR VS NORMAL SUMMARY\n")
cat("============================================================\n\n")

gene_arm <- gene_order[
  gene %in% rownames(cna),
  .(
    gene,
    chromosome,
    arm
  )
]

gene_arm[
  ,
  row_index :=
    match(
      gene,
      rownames(cna)
    )
]

arm_list <- split(
  gene_arm$row_index,
  paste0(
    gene_arm$chromosome,
    gene_arm$arm
  )
)

arm_summary_list <- lapply(
  names(arm_list),
  function(a) {

    idx <- arm_list[[a]]

    if (length(idx) < 10) {
      return(NULL)
    }

    tumor_arm <- colMeans(
      cna[
        idx,
        tumor_cells,
        drop = FALSE
      ]
    )

    normal_arm <- colMeans(
      cna[
        idx,
        normal_cells,
        drop = FALSE
      ]
    )

    data.table(
      chromosome_arm = a,

      genes =
        length(idx),

      tumor_median =
        median(
          tumor_arm
        ),

      normal_median =
        median(
          normal_arm
        ),

      tumor_minus_normal =
        median(tumor_arm) -
        median(normal_arm),

      tumor_positive_fraction =
        mean(
          tumor_arm > 0.05
        ),

      tumor_negative_fraction =
        mean(
          tumor_arm < -0.05
        ),

      normal_positive_fraction =
        mean(
          normal_arm > 0.05
        ),

      normal_negative_fraction =
        mean(
          normal_arm < -0.05
        )
    )
  }
)

arm_summary <- rbindlist(
  arm_summary_list,
  fill = TRUE
)

arm_summary[
  ,
  abs_delta :=
    abs(
      tumor_minus_normal
    )
]

setorder(
  arm_summary,
  -abs_delta
)

print(
  arm_summary,
  nrows = 100
)

fwrite(
  arm_summary,
  file.path(
    resultdir,
    "GSE197677_ESC06_17E2_chromosome_arm_summary.csv"
  ),
  bom = TRUE
)

# ------------------------------------------------------------
# 10. Candidate coherent arms
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("7. CANDIDATE COHERENT CNA ARMS — AUDIT ONLY\n")
cat("============================================================\n\n")

candidate_arms <- arm_summary[
  abs_delta >= 0.05 &
  (
    tumor_positive_fraction >= 0.50 |
    tumor_negative_fraction >= 0.50
  )
]

if (nrow(candidate_arms) == 0) {

  cat(
    "No chromosome arms pass the exploratory coherence flag.\n"
  )

} else {

  print(
    candidate_arms,
    nrows = 100
  )
}

fwrite(
  candidate_arms,
  file.path(
    resultdir,
    "GSE197677_ESC06_17E2_candidate_coherent_CNA_arms.csv"
  ),
  bom = TRUE
)

# ------------------------------------------------------------
# 11. Signal/correlation scatter
# ------------------------------------------------------------

cat("\n===== 8. CNA SIGNAL / CORRELATION PLOT =====\n\n")

plot_metrics <- copy(
  metrics
)

plot_metrics[
  ,
  display_group :=
    factor(
      display_group,
      levels = c(
        "ESN11",
        "ESN12",
        "ESN13",
        "ESN15",
        "ESC06"
      )
    )
]

p <- ggplot(
  plot_metrics,
  aes(
    x = tumor_profile_correlation,
    y = cna_signal,
    shape = display_group
  )
) +
  geom_point(
    alpha = 0.65,
    size = 1.5
  ) +
  labs(
    title =
      "ESC06 CNA Pilot Review",
    subtitle =
      "CNA signal versus correlation with ESC06 mean CNA profile",
    x =
      "Correlation with ESC06 Mean CNA Profile",
    y =
      "CNA Signal",
    shape =
      "Sample"
  ) +
  theme_classic()

ggsave(
  file.path(
    figdir,
    "GSE197677_ESC06_CNA_signal_vs_tumor_correlation.png"
  ),
  p,
  width = 8,
  height = 6,
  dpi = 300
)

# ------------------------------------------------------------
# 12. Final decision gate
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("9. REVIEW GATE\n")
cat("============================================================\n\n")

tumor_metrics <- metric_summary[
  display_group == "ESC06"
]

normal_metrics <- metric_summary[
  display_group != "ESC06"
]

normal_signal_median <- median(
  normal_metrics$median_CNA_signal,
  na.rm = TRUE
)

normal_cor_median <- median(
  normal_metrics$median_tumor_correlation,
  na.rm = TRUE
)

cat(
  sprintf(
    "ESC06 median CNA signal            : %.4f\n",
    tumor_metrics$median_CNA_signal
  )
)

cat(
  sprintf(
    "Median normal-reference CNA signal : %.4f\n",
    normal_signal_median
  )
)

cat(
  sprintf(
    "ESC06 median tumor-profile corr.   : %.4f\n",
    tumor_metrics$median_tumor_correlation
  )
)

cat(
  sprintf(
    "Median normal tumor-profile corr.  : %.4f\n",
    normal_cor_median
  )
)

cat(
  "Candidate coherent chromosome arms : ",
  nrow(candidate_arms),
  "\n",
  sep = ""
)

cat("\nInterpretation rule for manual review:\n")

cat(
  "1. Look for chromosome-scale contiguous bands in ESC06.\n"
)

cat(
  "2. Check that the same bands are not equally prominent in normal references.\n"
)

cat(
  "3. CNA signal alone is NOT sufficient for malignant calling.\n"
)

cat(
  "4. Tumor-profile correlation provides complementary evidence.\n"
)

cat(
  "5. Do NOT change InferCNA parameters based only on the desired biological answer.\n"
)

# ------------------------------------------------------------
# FINAL
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP 17E2 FINISHED ✓\n")
cat("============================================================\n\n")

cat("InferCNA was NOT rerun ✓\n")
cat("No parameters changed ✓\n")
cat("No malignant labels assigned ✓\n")
cat("No cells removed ✓\n")
cat("No other samples processed ✓\n")

cat("\nMain figure:\n")
cat(
  png_file,
  "\n"
)

cat("\nSTOP HERE.\n")
