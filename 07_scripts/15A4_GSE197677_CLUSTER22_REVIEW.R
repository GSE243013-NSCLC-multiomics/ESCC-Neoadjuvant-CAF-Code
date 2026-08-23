# ============================================================
# ESCC Neoadjuvant Project
# STEP 15A4
#
# Targeted review of ambiguous GSE197677 Cluster 22
#
# Questions:
# 1. Epithelial or Fibroblast/CAF?
# 2. EMT/hybrid phenotype?
# 3. Potential doublet-like QC profile?
# 4. Restricted to one/few samples?
#
# NO annotation write-back
# NO filtering
# NO normalization
# NO reclustering
# NO integration
# ============================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(ggplot2)
})

cat("\n============================================================\n")
cat("STEP 15A4: CLUSTER 22 TARGETED REVIEW\n")
cat("============================================================\n\n")

object_dir <- "03_objects/GSE197677"
table_dir  <- "06_tables/GSE197677/celltype"
figure_dir <- "05_figures/GSE197677/celltype"

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 1. Load previous validated candidate table
# ------------------------------------------------------------

candidate_file <- file.path(
  table_dir,
  "GSE197677_cluster_broad_celltype_candidates_FIXED.csv"
)

if (!file.exists(candidate_file)) {
  stop(
    "Missing FIXED candidate file. ",
    "Step 15A2 must be completed first."
  )
}

cand <- fread(candidate_file)

cat("===== 1. EXISTING CLUSTER 22 RESULT =====\n\n")

print(
  cand[
    cluster == "22" | cluster == 22
  ]
)

# ------------------------------------------------------------
# 2. Identify clear epithelial/fibroblast reference clusters
# ------------------------------------------------------------

epi_refs <- as.character(
  cand[
    suggested_broad_type == "Epithelial" &
    as.character(cluster) != "22",
    cluster
  ]
)

fib_refs <- as.character(
  cand[
    suggested_broad_type == "Fibroblast_CAF" &
    as.character(cluster) != "22",
    cluster
  ]
)

cat("\nClear Epithelial reference clusters:\n")
print(epi_refs)

cat("\nClear Fibroblast_CAF reference clusters:\n")
print(fib_refs)

# ------------------------------------------------------------
# 3. Locate and load large Seurat object
# ------------------------------------------------------------

files <- list.files(
  object_dir,
  pattern = "\\.rds$",
  full.names = TRUE,
  ignore.case = TRUE
)

hits <- files[
  grepl(
    "integrated",
    basename(files),
    ignore.case = TRUE
  )
]

if (length(hits) == 0) {
  stop("GSE197677 integrated RDS not found.")
}

object_file <- hits[1]

cat("\n===== 2. LOADING OBJECT =====\n\n")
cat("File: ", object_file, "\n", sep = "")
cat("Large RDS — loading may take several minutes.\n\n")

obj <- readRDS(object_file)

cat("Extracting metadata from slot...\n")

meta <- as.data.table(
  obj@meta.data,
  keep.rownames = "cell_id"
)

if (!"seurat_clusters" %in% colnames(meta)) {
  stop("seurat_clusters missing.")
}

clusters <- as.character(
  meta$seurat_clusters
)

idx22 <- which(
  clusters == "22"
)

cat(
  "Cluster 22 cells: ",
  length(idx22),
  "\n",
  sep = ""
)

if (length(idx22) == 0) {
  stop("Cluster 22 contains no cells.")
}

# ------------------------------------------------------------
# 4. Sample / treatment distribution
# ------------------------------------------------------------

cat("\n===== 3. CLUSTER 22 SAMPLE DISTRIBUTION =====\n\n")

meta[
  ,
  cluster_audit :=
    as.character(
      seurat_clusters
    )
]

sample_fields <- intersect(
  c(
    "sample",
    "CellType",
    "NAC",
    "subtype"
  ),
  names(meta)
)

cluster22_samples <- meta[
  cluster_audit == "22",
  .N,
  by = sample_fields
]

setnames(
  cluster22_samples,
  "N",
  "n_cells"
)

if ("sample" %in% names(cluster22_samples)) {
  setorder(
    cluster22_samples,
    -n_cells
  )
}

print(cluster22_samples)

fwrite(
  cluster22_samples,
  file.path(
    table_dir,
    "GSE197677_cluster22_sample_distribution.csv"
  ),
  bom = TRUE
)

# ------------------------------------------------------------
# 5. Cluster 22 fraction per sample
# ------------------------------------------------------------

if ("sample" %in% names(meta)) {

  total_by_sample <- meta[
    ,
    .(
      total_cells = .N
    ),
    by = sample
  ]

  c22_by_sample <- meta[
    cluster_audit == "22",
    .(
      cluster22_cells = .N
    ),
    by = sample
  ]

  sample_fraction <- merge(
    total_by_sample,
    c22_by_sample,
    by = "sample",
    all.x = TRUE
  )

  sample_fraction[
    is.na(cluster22_cells),
    cluster22_cells := 0
  ]

  sample_fraction[
    ,
    cluster22_fraction :=
      cluster22_cells /
      total_cells
  ]

  setorder(
    sample_fraction,
    -cluster22_fraction
  )

  cat("\nCluster 22 fraction by sample:\n\n")
  print(sample_fraction)

  fwrite(
    sample_fraction,
    file.path(
      table_dir,
      "GSE197677_cluster22_fraction_by_sample.csv"
    ),
    bom = TRUE
  )
}

# ------------------------------------------------------------
# 6. QC profile: Cluster 22 vs all other cells
# ------------------------------------------------------------

cat("\n===== 4. CLUSTER 22 QC PROFILE =====\n\n")

qc_fields <- intersect(
  c(
    "nCount_RNA",
    "nFeature_RNA",
    "percent.mt"
  ),
  names(meta)
)

qc_comparison <- list()

for (field in qc_fields) {

  x22 <- suppressWarnings(
    as.numeric(
      meta[
        cluster_audit == "22",
        get(field)
      ]
    )
  )

  xother <- suppressWarnings(
    as.numeric(
      meta[
        cluster_audit != "22",
        get(field)
      ]
    )
  )

  qc_comparison[[field]] <- data.table(
    metric = field,

    cluster22_median =
      median(
        x22,
        na.rm = TRUE
      ),

    other_cells_median =
      median(
        xother,
        na.rm = TRUE
      ),

    cluster22_p95 =
      as.numeric(
        quantile(
          x22,
          0.95,
          na.rm = TRUE
        )
      ),

    other_cells_p95 =
      as.numeric(
        quantile(
          xother,
          0.95,
          na.rm = TRUE
        )
      )
  )
}

qc_comparison <- rbindlist(
  qc_comparison,
  fill = TRUE
)

print(qc_comparison)

fwrite(
  qc_comparison,
  file.path(
    table_dir,
    "GSE197677_cluster22_QC_comparison.csv"
  ),
  bom = TRUE
)

# ------------------------------------------------------------
# 7. Extended marker panels
#
# Broad identity + squamous epithelial + fibroblast/CAF + EMT
# ------------------------------------------------------------

panels <- list(

  Epithelial_core = c(
    "EPCAM",
    "KRT8",
    "KRT18",
    "KRT19",
    "KRT7",
    "CLDN4",
    "CLDN7",
    "KRT5",
    "KRT14",
    "TP63",
    "SFN"
  ),

  Squamous_epithelial = c(
    "KRT4",
    "KRT13",
    "KRT6A",
    "KRT6B",
    "KRT16",
    "KRT17",
    "DSG3",
    "DSC3",
    "SERPINB3",
    "SERPINB4"
  ),

  Fibroblast_core = c(
    "COL1A1",
    "COL1A2",
    "COL3A1",
    "DCN",
    "LUM",
    "COL6A1",
    "COL6A2",
    "COL6A3",
    "PDGFRA",
    "COL14A1",
    "C7"
  ),

  CAF_activation = c(
    "FAP",
    "ACTA2",
    "TAGLN",
    "POSTN",
    "THY1",
    "CXCL12",
    "COL11A1"
  ),

  EMT_hybrid = c(
    "VIM",
    "FN1",
    "ZEB1",
    "ZEB2",
    "SNAI1",
    "SNAI2",
    "TWIST1",
    "TGFBI",
    "ITGA5"
  ),

  Pericyte = c(
    "RGS5",
    "CSPG4",
    "MCAM",
    "RBP1",
    "NOTCH3",
    "DES",
    "PDGFRB"
  )
)

all_markers <- unique(
  unlist(
    panels,
    use.names = FALSE
  )
)

# ------------------------------------------------------------
# 8. Retrieve raw RNA counts
# ------------------------------------------------------------

cat("\n===== 5. RAW RNA MARKER EXTRACTION =====\n\n")

counts <- GetAssayData(
  obj,
  assay = "RNA",
  layer = "counts"
)

markers_present <- intersect(
  all_markers,
  rownames(counts)
)

markers_missing <- setdiff(
  all_markers,
  rownames(counts)
)

cat(
  "Extended markers present: ",
  length(markers_present),
  " / ",
  length(all_markers),
  "\n",
  sep = ""
)

if (length(markers_missing) > 0) {

  cat("\nMissing markers:\n")
  print(markers_missing)
}

marker_counts <- counts[
  markers_present,
  ,
  drop = FALSE
]

# ------------------------------------------------------------
# 9. Compare cluster 22 + epithelial/fibroblast references
# ------------------------------------------------------------

comparison_clusters <- unique(
  c(
    epi_refs,
    fib_refs,
    "22"
  )
)

comparison_clusters <- comparison_clusters[
  comparison_clusters %in%
    unique(clusters)
]

cat("\nClusters in targeted comparison:\n")
print(comparison_clusters)

result_list <- list()
detect_list <- list()

for (cl in comparison_clusters) {

  idx <- which(
    clusters == cl
  )

  sub <- marker_counts[
    ,
    idx,
    drop = FALSE
  ]

  marker_sum <- Matrix::rowSums(
    sub
  )

  detection <- Matrix::rowSums(
    sub > 0
  ) / length(idx)

  total_umi <- sum(
    as.numeric(
      meta$nCount_RNA[idx]
    ),
    na.rm = TRUE
  )

  tmp <- data.table(
    gene = names(marker_sum),
    cluster = cl,
    n_cells = length(idx),
    counts = as.numeric(marker_sum),
    total_UMI = total_umi
  )

  tmp[
    ,
    CPM :=
      counts /
      total_UMI *
      1e6
  ]

  tmp[
    ,
    log1p_CPM :=
      log1p(CPM)
  ]

  result_list[[cl]] <- tmp

  detect_list[[cl]] <- data.table(
    gene = names(detection),
    cluster = cl,
    detection_fraction =
      as.numeric(detection)
  )

  rm(
    sub,
    marker_sum,
    detection
  )

  gc(verbose = FALSE)
}

expr <- rbindlist(
  result_list
)

detect <- rbindlist(
  detect_list
)

# ------------------------------------------------------------
# 10. Gene-wise z-score within targeted clusters
# ------------------------------------------------------------

expr[
  ,
  zscore := {

    x <- log1p_CPM

    if (
      length(unique(x)) <= 1 ||
      sd(x, na.rm = TRUE) == 0
    ) {

      rep(0, .N)

    } else {

      as.numeric(
        scale(x)
      )
    }

  },
  by = gene
]

lookup <- rbindlist(
  lapply(
    names(panels),
    function(p) {

      data.table(
        panel = p,
        gene = panels[[p]]
      )
    }
  )
)

expr <- merge(
  expr,
  lookup,
  by = "gene",
  all.x = TRUE
)

detect <- merge(
  detect,
  lookup,
  by = "gene",
  all.x = TRUE
)

# ------------------------------------------------------------
# 11. Cluster 22 panel-level evidence
# ------------------------------------------------------------

panel_expr <- expr[
  !is.na(panel),
  .(
    mean_zscore =
      mean(
        zscore,
        na.rm = TRUE
      ),

    mean_logCPM =
      mean(
        log1p_CPM,
        na.rm = TRUE
      )
  ),
  by = .(
    cluster,
    panel
  )
]

panel_detect <- detect[
  !is.na(panel),
  .(
    mean_detection =
      mean(
        detection_fraction,
        na.rm = TRUE
      )
  ),
  by = .(
    cluster,
    panel
  )
]

panel_summary <- merge(
  panel_expr,
  panel_detect,
  by = c(
    "cluster",
    "panel"
  ),
  all = TRUE
)

cat("\n============================================================\n")
cat("6. CLUSTER 22 EXTENDED PANEL SUMMARY\n")
cat("============================================================\n\n")

print(
  panel_summary[
    cluster == "22"
  ][
    order(-mean_zscore)
  ]
)

# ------------------------------------------------------------
# 12. Individual Cluster 22 markers
# ------------------------------------------------------------

cluster22_markers <- merge(
  expr[
    cluster == "22"
  ],
  detect[
    cluster == "22",
    .(
      gene,
      detection_fraction
    )
  ],
  by = "gene",
  all.x = TRUE
)

setorder(
  cluster22_markers,
  panel,
  -zscore
)

cat("\n============================================================\n")
cat("7. CLUSTER 22 INDIVIDUAL MARKERS\n")
cat("============================================================\n\n")

print(
  cluster22_markers[
    ,
    .(
      panel,
      gene,
      log1p_CPM = round(log1p_CPM, 3),
      zscore = round(zscore, 3),
      detection_fraction =
        round(
          detection_fraction,
          3
        )
    )
  ],
  nrows = 200
)

# ------------------------------------------------------------
# 13. Reference similarity
#
# Correlate Cluster 22 marker profile against:
# clear epithelial references
# clear fibroblast references
# ------------------------------------------------------------

wide <- dcast(
  expr,
  gene ~ cluster,
  value.var = "log1p_CPM"
)

genes <- wide$gene

mat <- as.matrix(
  wide[
    ,
    -1
  ]
)

rownames(mat) <- genes

similarity_rows <- list()

if (
  "22" %in%
    colnames(mat)
) {

  target <- mat[
    ,
    "22"
  ]

  for (
    ref in setdiff(
      colnames(mat),
      "22"
    )
  ) {

    r <- suppressWarnings(
      cor(
        target,
        mat[, ref],
        method = "spearman",
        use = "pairwise.complete.obs"
      )
    )

    ref_type <- fifelse(
      ref %in% epi_refs,
      "Epithelial_reference",

      fifelse(
        ref %in% fib_refs,
        "Fibroblast_reference",
        "Other"
      )
    )

    similarity_rows[[ref]] <- data.table(
      cluster22 = "22",
      reference_cluster = ref,
      reference_type = ref_type,
      spearman_r = r
    )
  }
}

similarity <- rbindlist(
  similarity_rows,
  fill = TRUE
)

setorder(
  similarity,
  -spearman_r
)

cat("\n============================================================\n")
cat("8. CLUSTER 22 REFERENCE SIMILARITY\n")
cat("============================================================\n\n")

print(similarity)

# ------------------------------------------------------------
# 14. Mean similarity by lineage
# ------------------------------------------------------------

similarity_summary <- similarity[
  ,
  .(
    mean_spearman =
      mean(
        spearman_r,
        na.rm = TRUE
      ),

    max_spearman =
      max(
        spearman_r,
        na.rm = TRUE
      )
  ),
  by = reference_type
]

cat("\nMean reference similarity:\n\n")
print(similarity_summary)

# ------------------------------------------------------------
# 15. Save tables
# ------------------------------------------------------------

fwrite(
  panel_summary,
  file.path(
    table_dir,
    "GSE197677_cluster22_extended_panel_summary.csv"
  ),
  bom = TRUE
)

fwrite(
  cluster22_markers,
  file.path(
    table_dir,
    "GSE197677_cluster22_extended_marker_evidence.csv"
  ),
  bom = TRUE
)

fwrite(
  similarity,
  file.path(
    table_dir,
    "GSE197677_cluster22_reference_similarity.csv"
  ),
  bom = TRUE
)

fwrite(
  similarity_summary,
  file.path(
    table_dir,
    "GSE197677_cluster22_reference_similarity_summary.csv"
  ),
  bom = TRUE
)

# ------------------------------------------------------------
# 16. Targeted marker heatmap
# ------------------------------------------------------------

plot_dt <- copy(
  expr
)

plot_dt[
  ,
  cluster :=
    factor(
      cluster,
      levels = comparison_clusters
    )
]

gene_order <- unlist(
  panels,
  use.names = FALSE
)

gene_order <- gene_order[
  gene_order %in%
    unique(plot_dt$gene)
]

plot_dt[
  ,
  gene :=
    factor(
      gene,
      levels = rev(
        gene_order
      )
    )
]

p <- ggplot(
  plot_dt,
  aes(
    x = cluster,
    y = gene,
    fill = zscore
  )
) +
  geom_tile() +
  facet_grid(
    panel ~ .,
    scales = "free_y",
    space = "free_y"
  ) +
  labs(
    title =
      "GSE197677 Cluster 22: Epithelial vs Fibroblast Review",
    subtitle =
      "Raw RNA count pseudobulk; comparison with clear reference clusters",
    x = "Cluster",
    y = "Marker",
    fill = "Z-score"
  ) +
  theme_classic() +
  theme(
    axis.text.x =
      element_text(
        angle = 45,
        hjust = 1
      ),

    strip.text.y =
      element_text(
        angle = 0
      )
  )

ggsave(
  file.path(
    figure_dir,
    "GSE197677_cluster22_targeted_review_heatmap.png"
  ),
  p,
  width = 10,
  height = 12,
  dpi = 300
)

ggsave(
  file.path(
    figure_dir,
    "GSE197677_cluster22_targeted_review_heatmap.pdf"
  ),
  p,
  width = 10,
  height = 12
)

# ------------------------------------------------------------
# 17. Detection heatmap
# ------------------------------------------------------------

detect_plot <- copy(
  detect
)

detect_plot[
  ,
  cluster :=
    factor(
      cluster,
      levels = comparison_clusters
    )
]

detect_plot[
  ,
  gene :=
    factor(
      gene,
      levels = rev(
        gene_order
      )
    )
]

p2 <- ggplot(
  detect_plot,
  aes(
    x = cluster,
    y = gene,
    fill = detection_fraction
  )
) +
  geom_tile() +
  facet_grid(
    panel ~ .,
    scales = "free_y",
    space = "free_y"
  ) +
  labs(
    title =
      "GSE197677 Cluster 22 Marker Detection",
    subtitle =
      "Fraction of cells expressing each marker",
    x = "Cluster",
    y = "Marker",
    fill = "Detection\nFraction"
  ) +
  theme_classic() +
  theme(
    axis.text.x =
      element_text(
        angle = 45,
        hjust = 1
      ),

    strip.text.y =
      element_text(
        angle = 0
      )
  )

ggsave(
  file.path(
    figure_dir,
    "GSE197677_cluster22_targeted_detection_heatmap.png"
  ),
  p2,
  width = 10,
  height = 12,
  dpi = 300
)

# ------------------------------------------------------------
# FINAL
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP 15A4 FINISHED ✓\n")
cat("============================================================\n\n")

cat(
  "Cluster 22 remains UNASSIGNED until manual review.\n\n"
)

cat("No metadata modified ✓\n")
cat("No broad labels written back ✓\n")
cat("No malignant epithelial calls ✓\n")
cat("No cells filtered ✓\n")
cat("No normalization ✓\n")
cat("No integration ✓\n")

cat("\nMain figure:\n")

cat(
  "05_figures/GSE197677/celltype/",
  "GSE197677_cluster22_targeted_review_heatmap.png\n",
  sep = ""
)

cat("\nSTOP HERE.\n")
