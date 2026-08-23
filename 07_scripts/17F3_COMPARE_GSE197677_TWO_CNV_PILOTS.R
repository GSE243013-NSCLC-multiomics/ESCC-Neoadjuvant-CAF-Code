# ============================================================
# ESCC Neoadjuvant Project
# STEP 17F3
#
# Compare two EXISTING GSE197677 InferCNA pilots:
#   ESC06
#   ESC15
#
# PURPOSE:
#   - recompute metrics with identical definitions
#   - compare tumor CNA signal vs normal references
#   - compute official-style tumor cnaCor
#   - compare chromosome-arm coherence
#   - decide whether GSE197677 is ready for expansion
#
# IMPORTANT:
#   NO infercna() rerun
#   NO parameter changes
#   NO findMalignant()
#   NO malignant labels
#   NO remaining samples processed
# ============================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(infercna)
  library(ggplot2)
})

cat("\n============================================================\n")
cat("STEP 17F3: COMPARE ESC06 + ESC15 CNV PILOTS\n")
cat("============================================================\n\n")

base_dir <- paste0(
  "06_tables/cross_dataset/",
  "malignant_epithelial"
)

outdir <- file.path(
  base_dir,
  "GSE197677_two_pilot_review"
)

figdir <- paste0(
  "05_figures/GSE197677/",
  "malignant_epithelial/",
  "two_pilot_review"
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

samples <- c(
  "ESC06",
  "ESC15"
)

# ============================================================
# 1. ROBUST FILE FINDER
# ============================================================

find_one <- function(
  sample,
  pattern,
  description
) {

  files <- list.files(
    base_dir,
    recursive = TRUE,
    full.names = TRUE
  )

  hit <- files[
    grepl(
      sample,
      files,
      ignore.case = TRUE
    ) &
    grepl(
      pattern,
      basename(files),
      ignore.case = TRUE
    )
  ]

  if (length(hit) == 0) {

    stop(
      "Cannot find ",
      description,
      " for ",
      sample
    )
  }

  # Prefer result files inside the corresponding pilot directory
  hit <- hit[
    order(
      nchar(hit)
    )
  ]

  hit[1]
}

# ============================================================
# 2. AUDIT FUNCTION
# ============================================================

audit_pilot <- function(sample) {

  cat("\n")
  cat("============================================================\n")
  cat("AUDIT ", sample, "\n", sep = "")
  cat("============================================================\n\n")

  # ----------------------------------------------------------
  # Locate files
  # ----------------------------------------------------------

  cna_file <- find_one(
    sample,
    "CNA_(all_cells|matrix)\\.rds$",
    "CNA matrix"
  )

  annotation_file <- find_one(
    sample,
    "cell_annotation\\.csv$",
    "cell annotation"
  )

  gene_file <- find_one(
    sample,
    "gene_table\\.csv$",
    "gene table"
  )

  cat(
    "CNA        : ",
    cna_file,
    "\n",
    sep = ""
  )

  cat(
    "Annotation : ",
    annotation_file,
    "\n",
    sep = ""
  )

  cat(
    "Gene table : ",
    gene_file,
    "\n",
    sep = ""
  )

  # ----------------------------------------------------------
  # Load
  # ----------------------------------------------------------

  cna <- readRDS(
    cna_file
  )

  ann <- fread(
    annotation_file
  )

  genes <- fread(
    gene_file
  )

  required_ann <- c(
    "cell_id",
    "sample",
    "infercna_role"
  )

  miss <- setdiff(
    required_ann,
    names(ann)
  )

  if (length(miss) > 0) {

    stop(
      sample,
      " annotation missing: ",
      paste(
        miss,
        collapse = ", "
      )
    )
  }

  tumor_cells <- ann[
    infercna_role ==
      "TUMOR_TEST",
    cell_id
  ]

  normal_cells <- ann[
    infercna_role ==
      "NORMAL_REFERENCE",
    cell_id
  ]

  tumor_cells <- intersect(
    tumor_cells,
    colnames(cna)
  )

  normal_cells <- intersect(
    normal_cells,
    colnames(cna)
  )

  if (length(tumor_cells) == 0) {
    stop(
      sample,
      ": no tumor cells found."
    )
  }

  if (length(normal_cells) == 0) {
    stop(
      sample,
      ": no normal reference cells found."
    )
  }

  cat(
    "\nCNA matrix   : ",
    nrow(cna),
    " genes x ",
    ncol(cna),
    " cells\n",
    sep = ""
  )

  cat(
    "Tumor cells  : ",
    length(tumor_cells),
    "\n",
    sep = ""
  )

  cat(
    "Normal cells : ",
    length(normal_cells),
    "\n",
    sep = ""
  )

  # ==========================================================
  # 3. CNA SIGNAL
  # ==========================================================

  signal <- infercna::cnaSignal(
    cna
  )

  if (is.null(names(signal))) {
    names(signal) <-
      colnames(cna)
  }

  # ==========================================================
  # 4. OFFICIAL-STYLE TUMOR CNA CORRELATION
  #
  # IMPORTANT:
  # calculate correlation WITHIN tumor cells only.
  # Normal references must not contribute to the tumor mean.
  # ==========================================================

  cna_tumor <- cna[
    ,
    tumor_cells,
    drop = FALSE
  ]

  tumor_cor <- infercna::cnaCor(
    cna_tumor
  )

  if (is.null(names(tumor_cor))) {
    names(tumor_cor) <-
      tumor_cells
  }

  # ----------------------------------------------------------
  # Normal-reference correlation to tumor mean
  #
  # QC only. This is NOT infercna cnaCor classification.
  # ----------------------------------------------------------

  tumor_mean_profile <- rowMeans(
    cna_tumor
  )

  normal_cor_to_tumor <- sapply(
    normal_cells,
    function(cell) {

      suppressWarnings(
        cor(
          cna[, cell],
          tumor_mean_profile,
          method = "pearson",
          use = "pairwise.complete.obs"
        )
      )
    }
  )

  # ==========================================================
  # 5. PER-CELL TABLE
  # ==========================================================

  tumor_dt <- data.table(

    dataset =
      "GSE197677",

    pilot =
      sample,

    cell_id =
      tumor_cells,

    role =
      "TUMOR_TEST",

    cna_signal =
      as.numeric(
        signal[
          tumor_cells
        ]
      ),

    cna_cor_to_tumor =
      as.numeric(
        tumor_cor[
          tumor_cells
        ]
      ),

    mean_abs_cna =
      colMeans(
        abs(
          cna_tumor
        )
      )
  )

  normal_dt <- data.table(

    dataset =
      "GSE197677",

    pilot =
      sample,

    cell_id =
      normal_cells,

    role =
      "NORMAL_REFERENCE",

    cna_signal =
      as.numeric(
        signal[
          normal_cells
        ]
      ),

    cna_cor_to_tumor =
      as.numeric(
        normal_cor_to_tumor[
          normal_cells
        ]
      ),

    mean_abs_cna =
      colMeans(
        abs(
          cna[
            ,
            normal_cells,
            drop = FALSE
          ]
        )
      )
  )

  cell_dt <- rbindlist(
    list(
      tumor_dt,
      normal_dt
    ),
    use.names = TRUE
  )

  # ==========================================================
  # 6. NORMAL-CALIBRATED DIAGNOSTICS
  #
  # These are QC summaries only.
  # NOT malignant thresholds.
  # ==========================================================

  normal_signal_95 <- as.numeric(
    quantile(
      normal_dt$cna_signal,
      0.95,
      na.rm = TRUE
    )
  )

  normal_signal_99 <- as.numeric(
    quantile(
      normal_dt$cna_signal,
      0.99,
      na.rm = TRUE
    )
  )

  normal_cor_95 <- as.numeric(
    quantile(
      normal_dt$cna_cor_to_tumor,
      0.95,
      na.rm = TRUE
    )
  )

  tumor_signal_median <- median(
    tumor_dt$cna_signal,
    na.rm = TRUE
  )

  normal_signal_median <- median(
    normal_dt$cna_signal,
    na.rm = TRUE
  )

  tumor_cor_median <- median(
    tumor_dt$cna_cor_to_tumor,
    na.rm = TRUE
  )

  normal_cor_median <- median(
    normal_dt$cna_cor_to_tumor,
    na.rm = TRUE
  )

  # ==========================================================
  # 7. CHROMOSOME-ARM AUDIT
  # ==========================================================

  required_gene <- c(
    "gene",
    "chromosome",
    "arm",
    "start_position"
  )

  miss_gene <- setdiff(
    required_gene,
    names(genes)
  )

  if (length(miss_gene) > 0) {

    stop(
      sample,
      " gene table missing: ",
      paste(
        miss_gene,
        collapse = ", "
      )
    )
  }

  gene_dt <- genes[
    gene %in%
      rownames(cna)
  ]

  gene_dt[
    ,
    row_index :=
      match(
        gene,
        rownames(cna)
      )
  ]

  gene_dt <- gene_dt[
    !is.na(row_index) &
    !is.na(chromosome) &
    !is.na(arm)
  ]

  gene_dt[
    ,
    chromosome_arm :=
      paste0(
        chromosome,
        arm
      )
  ]

  arm_indices <- split(
    gene_dt$row_index,
    gene_dt$chromosome_arm
  )

  arm_out <- lapply(
    names(arm_indices),
    function(a) {

      idx <- unique(
        arm_indices[[a]]
      )

      # avoid unstable very small arms
      if (length(idx) < 20) {
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

        pilot =
          sample,

        chromosome_arm =
          a,

        genes =
          length(idx),

        tumor_median =
          median(
            tumor_arm,
            na.rm = TRUE
          ),

        normal_median =
          median(
            normal_arm,
            na.rm = TRUE
          ),

        delta =
          median(
            tumor_arm,
            na.rm = TRUE
          ) -
          median(
            normal_arm,
            na.rm = TRUE
          ),

        tumor_gain_fraction =
          mean(
            tumor_arm > 0.05,
            na.rm = TRUE
          ),

        tumor_loss_fraction =
          mean(
            tumor_arm < -0.05,
            na.rm = TRUE
          ),

        normal_gain_fraction =
          mean(
            normal_arm > 0.05,
            na.rm = TRUE
          ),

        normal_loss_fraction =
          mean(
            normal_arm < -0.05,
            na.rm = TRUE
          )
      )
    }
  )

  arm_dt <- rbindlist(
    arm_out,
    fill = TRUE
  )

  arm_dt[
    ,
    abs_delta :=
      abs(delta)
  ]

  # Exploratory audit flag only
  arm_dt[
    ,
    coherent_audit_flag :=
      abs_delta >= 0.05 &
      (
        tumor_gain_fraction >= 0.50 |
        tumor_loss_fraction >= 0.50
      )
  ]

  setorder(
    arm_dt,
    -abs_delta
  )

  coherent <- arm_dt[
    coherent_audit_flag == TRUE
  ]

  # ==========================================================
  # 8. PILOT SUMMARY
  # ==========================================================

  summary <- data.table(

    dataset =
      "GSE197677",

    pilot =
      sample,

    tumor_cells =
      length(tumor_cells),

    normal_cells =
      length(normal_cells),

    genes =
      nrow(cna),

    tumor_median_CNA_signal =
      tumor_signal_median,

    normal_median_CNA_signal =
      normal_signal_median,

    CNA_signal_delta =
      tumor_signal_median -
      normal_signal_median,

    CNA_signal_ratio =
      ifelse(
        normal_signal_median > 0,
        tumor_signal_median /
          normal_signal_median,
        NA_real_
      ),

    tumor_fraction_above_normal95_signal =
      mean(
        tumor_dt$cna_signal >
          normal_signal_95,
        na.rm = TRUE
      ),

    tumor_fraction_above_normal99_signal =
      mean(
        tumor_dt$cna_signal >
          normal_signal_99,
        na.rm = TRUE
      ),

    tumor_median_cnaCor =
      tumor_cor_median,

    normal_median_correlation_to_tumor =
      normal_cor_median,

    correlation_separation =
      tumor_cor_median -
      normal_cor_median,

    tumor_fraction_above_normal95_correlation =
      mean(
        tumor_dt$cna_cor_to_tumor >
          normal_cor_95,
        na.rm = TRUE
      ),

    tumor_median_mean_abs_CNA =
      median(
        tumor_dt$mean_abs_cna,
        na.rm = TRUE
      ),

    normal_median_mean_abs_CNA =
      median(
        normal_dt$mean_abs_cna,
        na.rm = TRUE
      ),

    coherent_CNA_arms =
      nrow(coherent),

    coherent_arm_names =
      ifelse(
        nrow(coherent) == 0,
        "None",
        paste(
          coherent$chromosome_arm,
          collapse = "; "
        )
      )
  )

  list(
    summary = summary,
    cells = cell_dt,
    arms = arm_dt,
    coherent = coherent
  )
}

# ============================================================
# 9. RUN BOTH PILOTS
# ============================================================

results <- lapply(
  samples,
  audit_pilot
)

names(results) <-
  samples

summary_dt <- rbindlist(
  lapply(
    results,
    `[[`,
    "summary"
  ),
  fill = TRUE
)

cell_dt <- rbindlist(
  lapply(
    results,
    `[[`,
    "cells"
  ),
  fill = TRUE
)

arm_dt <- rbindlist(
  lapply(
    results,
    `[[`,
    "arms"
  ),
  fill = TRUE
)

coherent_dt <- rbindlist(
  lapply(
    results,
    `[[`,
    "coherent"
  ),
  fill = TRUE
)

# ============================================================
# 10. SIDE-BY-SIDE SUMMARY
# ============================================================

cat("\n")
cat("============================================================\n")
cat("10. TWO-PILOT COMPARISON — IDENTICAL DEFINITIONS\n")
cat("============================================================\n\n")

print(
  summary_dt,
  nrows = 100
)

fwrite(
  summary_dt,
  file.path(
    outdir,
    "STEP17F3_ESC06_ESC15_comparison.csv"
  ),
  bom = TRUE
)

fwrite(
  cell_dt,
  file.path(
    outdir,
    "STEP17F3_ESC06_ESC15_cell_metrics.csv"
  ),
  bom = TRUE
)

fwrite(
  arm_dt,
  file.path(
    outdir,
    "STEP17F3_ESC06_ESC15_all_arm_metrics.csv"
  ),
  bom = TRUE
)

fwrite(
  coherent_dt,
  file.path(
    outdir,
    "STEP17F3_ESC06_ESC15_coherent_arm_flags.csv"
  ),
  bom = TRUE
)

# ============================================================
# 11. COHERENT ARM COMPARISON
# ============================================================

cat("\n")
cat("============================================================\n")
cat("11. COHERENT CHROMOSOME-ARM FLAGS\n")
cat("============================================================\n\n")

if (nrow(coherent_dt) == 0) {

  cat(
    "No exploratory coherent arms detected in either pilot.\n"
  )

} else {

  print(
    coherent_dt[
      ,
      .(
        pilot,
        chromosome_arm,
        genes,
        tumor_median,
        normal_median,
        delta,
        tumor_gain_fraction,
        tumor_loss_fraction,
        normal_gain_fraction,
        normal_loss_fraction
      )
    ],
    nrows = 200
  )
}

# ============================================================
# 12. CNA SIGNAL PLOT
# ============================================================

plot_dt <- copy(
  cell_dt
)

plot_dt[
  ,
  group :=
    paste(
      pilot,
      role,
      sep = "_"
    )
]

p_signal <- ggplot(
  plot_dt,
  aes(
    x = group,
    y = cna_signal
  )
) +
  geom_boxplot(
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.15,
    size = 0.6,
    alpha = 0.25
  ) +
  labs(
    title =
      "GSE197677 Two-Pilot CNA Signal Review",
    subtitle =
      "ESC06 and ESC15 calculated with identical definitions",
    x =
      "Pilot / Cell Role",
    y =
      "CNA Signal"
  ) +
  theme_classic() +
  theme(
    axis.text.x =
      element_text(
        angle = 35,
        hjust = 1
      )
  )

ggsave(
  file.path(
    figdir,
    "STEP17F3_ESC06_ESC15_CNA_signal_comparison.png"
  ),
  p_signal,
  width = 8,
  height = 6,
  dpi = 300
)

# ============================================================
# 13. TUMOR-ONLY SIGNAL vs CNA COR
# ============================================================

tumor_plot_dt <- cell_dt[
  role ==
    "TUMOR_TEST"
]

p_cor <- ggplot(
  tumor_plot_dt,
  aes(
    x = cna_cor_to_tumor,
    y = cna_signal,
    shape = pilot
  )
) +
  geom_point(
    alpha = 0.55,
    size = 1.3
  ) +
  labs(
    title =
      "GSE197677 Tumor Epithelial CNA Review",
    subtitle =
      "Official-style cell-to-tumor CNA correlation",
    x =
      "Cell-to-Tumor CNA Correlation",
    y =
      "CNA Signal",
    shape =
      "Pilot"
  ) +
  theme_classic()

ggsave(
  file.path(
    figdir,
    "STEP17F3_ESC06_ESC15_tumor_signal_vs_cnaCor.png"
  ),
  p_cor,
  width = 7,
  height = 6,
  dpi = 300
)

# ============================================================
# 14. INTERPRETATION GATE
# ============================================================

cat("\n")
cat("============================================================\n")
cat("12. INTERPRETATION GATE\n")
cat("============================================================\n\n")

for (s in samples) {

  z <- summary_dt[
    pilot == s
  ]

  cat(
    "\n--- ",
    s,
    " ---\n",
    sep = ""
  )

  cat(
    sprintf(
      "Tumor median CNA signal              : %.5f\n",
      z$tumor_median_CNA_signal
    )
  )

  cat(
    sprintf(
      "Normal median CNA signal             : %.5f\n",
      z$normal_median_CNA_signal
    )
  )

  cat(
    sprintf(
      "Signal delta                         : %.5f\n",
      z$CNA_signal_delta
    )
  )

  cat(
    sprintf(
      "Tumor cells > normal 95th signal     : %.1f%%\n",
      100 *
        z$tumor_fraction_above_normal95_signal
    )
  )

  cat(
    sprintf(
      "Tumor median official-style cnaCor   : %.3f\n",
      z$tumor_median_cnaCor
    )
  )

  cat(
    sprintf(
      "Normal corr. to tumor median         : %.3f\n",
      z$normal_median_correlation_to_tumor
    )
  )

  cat(
    sprintf(
      "Correlation separation              : %.3f\n",
      z$correlation_separation
    )
  )

  cat(
    "Coherent CNA arms                    : ",
    z$coherent_CNA_arms,
    " [",
    z$coherent_arm_names,
    "]\n",
    sep = ""
  )
}

cat("\n")
cat("IMPORTANT INTERPRETATION RULES:\n")
cat("\n")

cat(
  "1. A negative normal-to-tumor correlation is NOT evidence of a treatment effect.\n"
)

cat(
  "2. CNA signal alone is NOT a malignant-cell classifier.\n"
)

cat(
  "3. Cell-to-tumor cnaCor should be interpreted together with CNA signal.\n"
)

cat(
  "4. Chromosome-arm coherence should be visible across contiguous genes and multiple tumor cells.\n"
)

cat(
  "5. These thresholds are audit flags only; no malignant labels are assigned.\n"
)

cat(
  "6. Do not expand to the remaining 8 samples until both pilot heatmaps and these metrics are reviewed.\n"
)

# ============================================================
# FINAL
# ============================================================

cat("\n")
cat("============================================================\n")
cat("STEP 17F3 FINISHED\n")
cat("TWO-PILOT CNA COMPARISON COMPLETE\n")
cat("============================================================\n\n")

cat("InferCNA rerun              : NO\n")
cat("Parameters changed          : NO\n")
cat("findMalignant()             : NO\n")
cat("Malignant labels assigned   : NO\n")
cat("Remaining samples processed : NO\n")

cat("\nMain outputs:\n")

cat(
  "  ",
  file.path(
    outdir,
    "STEP17F3_ESC06_ESC15_comparison.csv"
  ),
  "\n",
  sep = ""
)

cat(
  "  ",
  file.path(
    figdir,
    "STEP17F3_ESC06_ESC15_tumor_signal_vs_cnaCor.png"
  ),
  "\n",
  sep = ""
)

cat("\nSTOP HERE.\n")
