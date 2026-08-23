# ============================================================
# ESCC Neoadjuvant Project
# STEP 15A2
#
# Fix GSE197677 broad cell-type scoring using RAW RNA counts.
#
# Strategy:
#   raw RNA counts
#      -> cluster pseudobulk marker counts
#      -> CPM
#      -> log1p(CPM)
#      -> gene-wise Z score across clusters
#
# ALSO:
#   marker detection fraction per cluster
#
# NO cell filtering
# NO normalization of individual cells
# NO reclustering
# NO metadata modification
# NO malignant calling
# NO cross-dataset integration
# ============================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(ggplot2)
})

cat("\n============================================================\n")
cat("STEP 15A2: GSE197677 RAW-COUNT LINEAGE SCORING\n")
cat("============================================================\n\n")

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------

object_dir <- "03_objects/GSE197677"

table_dir <- "06_tables/GSE197677/celltype"
figure_dir <- "05_figures/GSE197677/celltype"

dir.create(
  table_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# 1. Locate integrated object
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

cat("===== 1. LOADING OBJECT =====\n\n")
cat("File: ", object_file, "\n", sep = "")

obj <- readRDS(object_file)

# Extract metadata directly from slot (faster, no save/reload)
cat("Extracting metadata from slot...\n")

md <- tryCatch(
  as.data.table(
    obj@meta.data,
    keep.rownames = "cell_id"
  ),
  error = function(e) {
    stop("Cannot extract metadata: ", conditionMessage(e))
  }
)

cat(
  "Metadata cells: ",
  nrow(md),
  "\n",
  sep = ""
)

if (!"seurat_clusters" %in% colnames(md)) {
  stop("seurat_clusters not found.")
}

if (!"nCount_RNA" %in% colnames(md)) {
  stop("nCount_RNA not found.")
}

# ------------------------------------------------------------
# 2. Cluster definitions
# ------------------------------------------------------------

cluster_vector <- as.character(
  md$seurat_clusters
)

valid <- (
  !is.na(cluster_vector) &
  cluster_vector != ""
)

active_clusters <- sort(
  unique(
    as.integer(
      cluster_vector[valid]
    )
  )
)

active_clusters <- as.character(
  active_clusters
)

cat("\n===== 2. ACTIVE CLUSTERS =====\n\n")

cluster_counts <- data.table(
  cluster = cluster_vector
)[
  cluster %in% active_clusters,
  .(n_cells = .N),
  by = cluster
]

cluster_counts[
  ,
  cluster_numeric :=
    as.integer(cluster)
]

setorder(
  cluster_counts,
  cluster_numeric
)

cluster_counts[
  ,
  cluster_numeric := NULL
]

print(cluster_counts)

cat(
  "\nOccupied clusters: ",
  nrow(cluster_counts),
  "\n",
  sep = ""
)

if (nrow(cluster_counts) != 25) {
  cat(
    "WARNING: expected 25 occupied clusters (0-24).\n"
  )
}

# ------------------------------------------------------------
# 3. Broad canonical marker panels
# ------------------------------------------------------------

marker_panels <- list(

  Epithelial = c(
    "EPCAM",
    "KRT8",
    "KRT18",
    "KRT19",
    "KRT5",
    "KRT14",
    "TP63",
    "SFN"
  ),

  T_NK = c(
    "CD3D",
    "CD3E",
    "TRAC",
    "TRBC1",
    "TRBC2",
    "NKG7",
    "KLRD1",
    "GNLY",
    "CCL5"
  ),

  Myeloid = c(
    "LST1",
    "TYROBP",
    "FCER1G",
    "CTSS",
    "AIF1",
    "LILRB1",
    "C1QC",
    "C1QA",
    "S100A8",
    "S100A9"
  ),

  B_plasma = c(
    "CD79A",
    "CD79B",
    "MS4A1",
    "CD37",
    "CD74",
    "HLA-DRA",
    "MZB1",
    "JCHAIN",
    "SDC1"
  ),

  Fibroblast_CAF = c(
    "COL1A1",
    "COL1A2",
    "COL3A1",
    "DCN",
    "LUM",
    "COL6A1",
    "COL6A2",
    "PDGFRA",
    "FAP",
    "ACTA2",
    "TAGLN"
  ),

  Endothelial = c(
    "PECAM1",
    "VWF",
    "EMCN",
    "KDR",
    "RAMP2",
    "PLVAP",
    "ENG"
  ),

  Mast = c(
    "TPSAB1",
    "TPSB2",
    "KIT",
    "CPA3",
    "HDC"
  ),

  Pericyte_other_stromal = c(
    "RGS5",
    "CSPG4",
    "MCAM",
    "RBP1",
    "NOTCH3",
    "DES"
  )
)

all_markers <- unique(
  unlist(
    marker_panels,
    use.names = FALSE
  )
)

# ------------------------------------------------------------
# 4. Get RAW RNA counts
# ------------------------------------------------------------

cat("\n===== 3. RNA RAW COUNTS =====\n\n")

if (!"RNA" %in% Assays(obj)) {
  stop("RNA assay missing.")
}

counts <- tryCatch(

  GetAssayData(
    obj,
    assay = "RNA",
    layer = "counts"
  ),

  error = function(e) {
    stop(
      "Could not retrieve RNA counts: ",
      conditionMessage(e)
    )
  }
)

cat(
  "RNA count matrix: ",
  nrow(counts),
  " genes x ",
  ncol(counts),
  " cells\n",
  sep = ""
)

if (ncol(counts) != nrow(md)) {
  stop(
    "Counts columns and metadata rows differ."
  )
}

# ------------------------------------------------------------
# 5. Marker availability
# ------------------------------------------------------------

markers_present <- intersect(
  all_markers,
  rownames(counts)
)

markers_missing <- setdiff(
  all_markers,
  rownames(counts)
)

cat("\n===== 4. MARKER AVAILABILITY =====\n\n")

cat(
  "Requested: ",
  length(all_markers),
  "\n",
  sep = ""
)

cat(
  "Present  : ",
  length(markers_present),
  "\n",
  sep = ""
)

cat(
  "Missing  : ",
  length(markers_missing),
  "\n",
  sep = ""
)

if (length(markers_missing) > 0) {

  cat(
    "\nMissing markers:\n"
  )

  print(markers_missing)
}

marker_availability <- rbindlist(
  lapply(
    names(marker_panels),
    function(lineage) {

      genes <- marker_panels[[lineage]]

      data.table(
        lineage = lineage,
        gene = genes,
        present = genes %in% rownames(counts)
      )
    }
  )
)

fwrite(
  marker_availability,
  file.path(
    table_dir,
    "GSE197677_marker_availability_raw_counts.csv"
  ),
  bom = TRUE
)

# ------------------------------------------------------------
# 6. Extract ONLY marker rows
# ------------------------------------------------------------

marker_counts <- counts[
  markers_present,
  ,
  drop = FALSE
]

cat(
  "\nMarker count matrix: ",
  nrow(marker_counts),
  " x ",
  ncol(marker_counts),
  "\n",
  sep = ""
)

# ------------------------------------------------------------
# 7. Compute cluster pseudobulk marker counts
#
# Important:
# library size uses total RNA UMI from metadata,
# not only marker genes.
# ------------------------------------------------------------

cat("\n===== 5. CLUSTER PSEUDOBULK =====\n\n")

pb_count_list <- list()
detect_list <- list()
cluster_library <- list()

for (cl in active_clusters) {

  idx <- which(
    cluster_vector == cl
  )

  if (length(idx) == 0) {
    next
  }

  cat(
    "Cluster ",
    cl,
    ": ",
    length(idx),
    " cells\n",
    sep = ""
  )

  sub_counts <- marker_counts[
    ,
    idx,
    drop = FALSE
  ]

  marker_sum <- Matrix::rowSums(
    sub_counts
  )

  detection <- Matrix::rowSums(
    sub_counts > 0
  ) / length(idx)

  # Total RNA UMI library size for this cluster
  library_size <- sum(
    as.numeric(
      md$nCount_RNA[idx]
    ),
    na.rm = TRUE
  )

  pb_count_list[[cl]] <- data.table(
    gene = names(marker_sum),
    cluster = cl,
    marker_counts = as.numeric(marker_sum),
    cluster_library_size = library_size
  )

  detect_list[[cl]] <- data.table(
    gene = names(detection),
    cluster = cl,
    detection_fraction = as.numeric(detection)
  )

  cluster_library[[cl]] <- data.table(
    cluster = cl,
    n_cells = length(idx),
    total_UMI = library_size
  )

  rm(
    sub_counts,
    marker_sum,
    detection
  )

  gc(verbose = FALSE)
}

pb_counts <- rbindlist(
  pb_count_list,
  fill = TRUE
)

detect_dt <- rbindlist(
  detect_list,
  fill = TRUE
)

library_dt <- rbindlist(
  cluster_library,
  fill = TRUE
)

# ------------------------------------------------------------
# 8. CPM + log1p CPM
# ------------------------------------------------------------

pb_counts[
  ,
  CPM :=
    (marker_counts /
       cluster_library_size) *
    1e6
]

pb_counts[
  ,
  log1p_CPM :=
    log1p(CPM)
]

cat("\n===== 6. RAW PSEUDOBULK SANITY CHECK =====\n\n")

summary_by_gene <- pb_counts[
  ,
  .(
    min_logCPM = min(log1p_CPM),
    max_logCPM = max(log1p_CPM),
    range_logCPM =
      max(log1p_CPM) -
      min(log1p_CPM)
  ),
  by = gene
]

nonzero_variance_genes <- summary_by_gene[
  range_logCPM > 1e-8,
  .N
]

cat(
  "Markers with non-zero cluster variation: ",
  nonzero_variance_genes,
  " / ",
  nrow(summary_by_gene),
  "\n",
  sep = ""
)

print(
  summary_by_gene[
    order(-range_logCPM)
  ][1:min(20, .N)]
)

if (nonzero_variance_genes == 0) {

  stop(
    "All markers still have zero cluster variation. ",
    "Do NOT proceed."
  )
}

# ------------------------------------------------------------
# 9. Gene-wise z-score across clusters
# ------------------------------------------------------------

pb_counts[
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

cat(
  "\nZ-score range: ",
  round(min(pb_counts$zscore), 3),
  " to ",
  round(max(pb_counts$zscore), 3),
  "\n",
  sep = ""
)

# ------------------------------------------------------------
# 10. Attach lineage
# ------------------------------------------------------------

lookup <- rbindlist(
  lapply(
    names(marker_panels),
    function(lineage) {

      data.table(
        lineage = lineage,
        gene = marker_panels[[lineage]]
      )
    }
  )
)

pb_scored <- merge(
  pb_counts,
  lookup,
  by = "gene",
  all.x = TRUE
)

detect_scored <- merge(
  detect_dt,
  lookup,
  by = "gene",
  all.x = TRUE
)

# ------------------------------------------------------------
# 11. Lineage score = mean marker gene z-score
# ------------------------------------------------------------

lineage_scores <- pb_scored[
  !is.na(lineage),
  .(
    lineage_score =
      mean(
        zscore,
        na.rm = TRUE
      ),

    mean_logCPM =
      mean(
        log1p_CPM,
        na.rm = TRUE
      ),

    n_markers_used =
      uniqueN(gene)
  ),
  by = .(
    cluster,
    lineage
  )
]

# Detection evidence
detection_scores <- detect_scored[
  !is.na(lineage),
  .(
    mean_detection_fraction =
      mean(
        detection_fraction,
        na.rm = TRUE
      )
  ),
  by = .(
    cluster,
    lineage
  )
]

lineage_scores <- merge(
  lineage_scores,
  detection_scores,
  by = c(
    "cluster",
    "lineage"
  ),
  all.x = TRUE
)

# ------------------------------------------------------------
# 12. Rank broad lineage candidates
# ------------------------------------------------------------

setorder(
  lineage_scores,
  cluster,
  -lineage_score
)

ranked <- lineage_scores[
  ,
  {

    ord <- order(
      -lineage_score
    )

    L <- lineage[ord]
    S <- lineage_score[ord]
    D <- mean_detection_fraction[ord]

    .(
      suggested_broad_type =
        L[1],

      top_score =
        S[1],

      top_detection =
        D[1],

      second_type =
        if (length(L) >= 2)
          L[2]
        else
          NA_character_,

      second_score =
        if (length(S) >= 2)
          S[2]
        else
          NA_real_,

      score_margin =
        if (length(S) >= 2)
          S[1] - S[2]
        else
          NA_real_
    )
  },
  by = cluster
]

ranked <- merge(
  cluster_counts,
  ranked,
  by = "cluster",
  all.x = TRUE
)

ranked[
  ,
  review_status :=
    fcase(

      is.na(score_margin),
      "MANUAL_REVIEW_REQUIRED",

      score_margin >= 0.50 &
        top_detection >= 0.10,
      "STRONG_CANDIDATE",

      score_margin >= 0.25,
      "MODERATE_CANDIDATE",

      default =
        "MANUAL_REVIEW_REQUIRED"
    )
]

ranked[
  ,
  cluster_number :=
    as.integer(cluster)
]

setorder(
  ranked,
  cluster_number
)

ranked[
  ,
  cluster_number := NULL
]

cat("\n============================================================\n")
cat("7. BROAD CELL-TYPE CANDIDATES — FIXED\n")
cat("============================================================\n\n")

print(ranked)

# ------------------------------------------------------------
# 13. Top 3 lineage scores for every cluster
# ------------------------------------------------------------

top3 <- lineage_scores[
  order(
    as.integer(cluster),
    -lineage_score
  ),
  head(.SD, 3),
  by = cluster
]

cat("\n===== 8. TOP 3 LINEAGES PER CLUSTER =====\n\n")

print(top3)

# ------------------------------------------------------------
# 14. Save tables
# ------------------------------------------------------------

fwrite(
  pb_counts,
  file.path(
    table_dir,
    "GSE197677_rawcount_marker_pseudobulk.csv"
  ),
  bom = TRUE
)

fwrite(
  detect_dt,
  file.path(
    table_dir,
    "GSE197677_marker_detection_fraction.csv"
  ),
  bom = TRUE
)

fwrite(
  lineage_scores,
  file.path(
    table_dir,
    "GSE197677_cluster_lineage_scores_RAWCOUNTS.csv"
  ),
  bom = TRUE
)

fwrite(
  ranked,
  file.path(
    table_dir,
    "GSE197677_cluster_broad_celltype_candidates_FIXED.csv"
  ),
  bom = TRUE
)

fwrite(
  top3,
  file.path(
    table_dir,
    "GSE197677_cluster_top3_lineages.csv"
  ),
  bom = TRUE
)

fwrite(
  summary_by_gene,
  file.path(
    table_dir,
    "GSE197677_marker_cluster_variation_audit.csv"
  ),
  bom = TRUE
)

# ------------------------------------------------------------
# 15. Marker Z-score heatmap
# ------------------------------------------------------------

heat <- copy(
  pb_scored
)

heat[
  ,
  cluster_number :=
    as.integer(cluster)
]

heat[
  ,
  cluster :=
    factor(
      cluster,
      levels =
        as.character(
          sort(
            unique(
              cluster_number
            )
          )
        )
    )
]

gene_order <- unlist(
  marker_panels,
  use.names = FALSE
)

gene_order <- gene_order[
  gene_order %in%
    unique(heat$gene)
]

heat[
  ,
  gene :=
    factor(
      gene,
      levels = rev(
        gene_order
      )
    )
]

p_marker <- ggplot(
  heat,
  aes(
    x = cluster,
    y = gene,
    fill = zscore
  )
) +
  geom_tile() +
  facet_grid(
    lineage ~ .,
    scales = "free_y",
    space = "free_y"
  ) +
  labs(
    title =
      "GSE197677 Broad Cell-Type Marker Audit",
    subtitle =
      "RNA raw-count cluster pseudobulk; gene-wise Z-score",
    x = "Seurat Cluster",
    y = "Canonical Marker",
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
    "GSE197677_broad_celltype_marker_heatmap_FIXED.png"
  ),
  p_marker,
  width = 12,
  height = 13,
  dpi = 300
)

ggsave(
  file.path(
    figure_dir,
    "GSE197677_broad_celltype_marker_heatmap_FIXED.pdf"
  ),
  p_marker,
  width = 12,
  height = 13
)

# ------------------------------------------------------------
# 16. Detection-fraction heatmap
# ------------------------------------------------------------

detect_plot <- copy(
  detect_scored
)

detect_plot[
  ,
  cluster_number :=
    as.integer(cluster)
]

detect_plot[
  ,
  cluster :=
    factor(
      cluster,
      levels =
        as.character(
          sort(
            unique(
              cluster_number
            )
          )
        )
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

p_detect <- ggplot(
  detect_plot,
  aes(
    x = cluster,
    y = gene,
    fill = detection_fraction
  )
) +
  geom_tile() +
  facet_grid(
    lineage ~ .,
    scales = "free_y",
    space = "free_y"
  ) +
  labs(
    title =
      "GSE197677 Canonical Marker Detection Fraction",
    subtitle =
      "Fraction of cells expressing each marker within each cluster",
    x = "Seurat Cluster",
    y = "Canonical Marker",
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
    "GSE197677_marker_detection_heatmap_FIXED.png"
  ),
  p_detect,
  width = 12,
  height = 13,
  dpi = 300
)

ggsave(
  file.path(
    figure_dir,
    "GSE197677_marker_detection_heatmap_FIXED.pdf"
  ),
  p_detect,
  width = 12,
  height = 13
)

# ------------------------------------------------------------
# 17. Lineage score heatmap
# ------------------------------------------------------------

score_plot <- copy(
  lineage_scores
)

score_plot[
  ,
  cluster_number :=
    as.integer(cluster)
]

score_plot[
  ,
  cluster :=
    factor(
      cluster,
      levels =
        as.character(
          sort(
            unique(
              cluster_number
            )
          )
        )
    )
]

p_score <- ggplot(
  score_plot,
  aes(
    x = cluster,
    y = lineage,
    fill = lineage_score
  )
) +
  geom_tile() +
  labs(
    title =
      "GSE197677 Broad Lineage Scores",
    subtitle =
      "Raw RNA count pseudobulk marker Z-score",
    x = "Seurat Cluster",
    y = "Broad Lineage",
    fill = "Lineage\nScore"
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
    figure_dir,
    "GSE197677_cluster_lineage_score_heatmap_FIXED.png"
  ),
  p_score,
  width = 11,
  height = 6,
  dpi = 300
)

ggsave(
  file.path(
    figure_dir,
    "GSE197677_cluster_lineage_score_heatmap_FIXED.pdf"
  ),
  p_score,
  width = 11,
  height = 6
)

# ------------------------------------------------------------
# 18. Final sanity check
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("FINAL SANITY CHECK\n")
cat("============================================================\n\n")

score_range <- range(
  lineage_scores$lineage_score,
  na.rm = TRUE
)

cat(
  "Lineage score range: ",
  round(score_range[1], 3),
  " to ",
  round(score_range[2], 3),
  "\n",
  sep = ""
)

all_zero <- all(
  abs(
    lineage_scores$lineage_score
  ) < 1e-10
)

cat(
  "All lineage scores zero: ",
  all_zero,
  "\n",
  sep = ""
)

if (all_zero) {

  cat(
    "\nERROR: scores remain zero.\n"
  )

  cat(
    "DO NOT USE CANDIDATE LABELS.\n"
  )

} else {

  cat(
    "\n✓ Raw-count scoring produced non-zero lineage separation.\n"
  )
}

# ------------------------------------------------------------
# FINAL
# ------------------------------------------------------------

cat("\n============================================================\n")

if (!all_zero) {

  cat("STEP 15A2 FINISHED ✓\n")
  cat("ZERO-SCORE PROBLEM RESOLVED ✓\n")

} else {

  cat("STEP 15A2 FINISHED WITH ERROR FLAG ⚠\n")
}

cat("============================================================\n\n")

cat("No cells filtered ✓\n")
cat("No individual-cell normalization performed ✓\n")
cat("No reclustering ✓\n")
cat("No Seurat metadata modified ✓\n")
cat("No malignant epithelial calls ✓\n")
cat("No cross-dataset integration ✓\n")

cat("\nImportant figures:\n")

cat(
  "  05_figures/GSE197677/celltype/",
  "GSE197677_broad_celltype_marker_heatmap_FIXED.png\n",
  sep = ""
)

cat(
  "  05_figures/GSE197677/celltype/",
  "GSE197677_marker_detection_heatmap_FIXED.png\n",
  sep = ""
)

cat(
  "  05_figures/GSE197677/celltype/",
  "GSE197677_cluster_lineage_score_heatmap_FIXED.png\n",
  sep = ""
)

cat("\nSTOP HERE.\n")

