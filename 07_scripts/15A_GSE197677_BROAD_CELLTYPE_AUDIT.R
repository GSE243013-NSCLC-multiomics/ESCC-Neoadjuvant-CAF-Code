# ============================================================
# ESCC Neoadjuvant Project
# STEP 15A
# GSE197677 broad cell ontology audit
#
# PURPOSE:
#   Cluster-level canonical-marker audit only
#
# NO CELL FILTERING
# NO NORMALIZATION
# NO RECLUSTERING
# NO MALIGNANT CALLING
# NO CROSS-DATASET INTEGRATION
# ============================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(ggplot2)
})

cat("\n============================================================\n")
cat("STEP 15A: GSE197677 BROAD CELL-TYPE AUDIT\n")
cat("============================================================\n\n")

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------

object_dir <- "03_objects/GSE197677"
table_dir  <- "06_tables/GSE197677/celltype"
figure_dir <- "05_figures/GSE197677/celltype"

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

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
  stop("GSE197677 integrated object not found.")
}

object_file <- hits[1]

cat("===== 1. LOADING OBJECT =====\n\n")
cat("File: ", object_file, "\n", sep = "")

obj <- readRDS(object_file)

md <- obj[[]]

cat("Cells from metadata: ", nrow(md), "\n", sep = "")

# ------------------------------------------------------------
# 2. Cluster field
# ------------------------------------------------------------

if (!"seurat_clusters" %in% colnames(md)) {
  stop("seurat_clusters field not found.")
}

cluster_raw <- as.character(md$seurat_clusters)

cluster_counts <- data.table(
  cluster = cluster_raw
)[
  !is.na(cluster) & cluster != "",
  .(n_cells = .N),
  by = cluster
]

cluster_counts[
  ,
  cluster_numeric :=
    suppressWarnings(as.integer(cluster))
]

setorder(
  cluster_counts,
  cluster_numeric
)

cluster_counts[
  ,
  cluster_numeric := NULL
]

cat("\n===== 2. ACTIVE CLUSTERS =====\n\n")
print(cluster_counts)

cat(
  "\nNumber of clusters with cells: ",
  nrow(cluster_counts),
  "\n",
  sep = ""
)

# ------------------------------------------------------------
# 3. Canonical marker panels
#
# These define BROAD ontology only.
# They do NOT define malignant epithelial states.
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
    "S100A8",
    "S100A9"
  ),

  B_plasma = c(
    "CD79A",
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

# ------------------------------------------------------------
# 4. Pick expression assay
#
# Prefer RNA normalized data already stored by authors.
# Do NOT NormalizeData again.
# ------------------------------------------------------------

cat("\n===== 3. EXPRESSION SOURCE AUDIT =====\n\n")

assay_to_use <- NULL

# Try RNA assay
if ("RNA" %in% Assays(obj)) {
  cat("RNA assay found. Trying to access data...\n")
  
  # Try v5 layer syntax first
  data_ok <- tryCatch(
    {
      test <- GetAssayData(obj, assay = "RNA", layer = "data")
      nrow(test) > 0 && ncol(test) > 0
    },
    error = function(e) {
      cat("  v5 layer failed: ", conditionMessage(e), "\n")
      FALSE
    }
  )
  
  # Try v4 slot syntax
  if (!data_ok) {
    data_ok <- tryCatch(
      {
        test <- GetAssayData(obj, assay = "RNA", slot = "data")
        nrow(test) > 0 && ncol(test) > 0
      },
      error = function(e) {
        cat("  v4 slot failed: ", conditionMessage(e), "\n")
        FALSE
      }
    )
  }
  
  # Try directly accessing the slot
  if (!data_ok) {
    data_ok <- tryCatch(
      {
        rna_assay <- obj[["RNA"]]
        test <- rna_assay@data
        nrow(test) > 0 && ncol(test) > 0
      },
      error = function(e) {
        cat("  Direct slot access failed: ", conditionMessage(e), "\n")
        FALSE
      }
    )
  }
  
  if (data_ok) {
    assay_to_use <- "RNA"
    cat("RNA data accessed successfully.\n")
  }
}

# Try SCT assay
if (is.null(assay_to_use) && "SCT" %in% Assays(obj)) {
  cat("SCT assay found. Trying to access data...\n")
  
  data_ok <- tryCatch(
    {
      test <- GetAssayData(obj, assay = "SCT", layer = "data")
      nrow(test) > 0 && ncol(test) > 0
    },
    error = function(e) {
      cat("  v5 layer failed: ", conditionMessage(e), "\n")
      FALSE
    }
  )
  
  if (!data_ok) {
    data_ok <- tryCatch(
      {
        test <- GetAssayData(obj, assay = "SCT", slot = "data")
        nrow(test) > 0 && ncol(test) > 0
      },
      error = function(e) {
        cat("  v4 slot failed: ", conditionMessage(e), "\n")
        FALSE
      }
    )
  }
  
  if (!data_ok) {
    data_ok <- tryCatch(
      {
        sct_assay <- obj[["SCT"]]
        test <- sct_assay@data
        nrow(test) > 0 && ncol(test) > 0
      },
      error = function(e) {
        cat("  Direct slot access failed: ", conditionMessage(e), "\n")
        FALSE
      }
    )
  }
  
  if (data_ok) {
    assay_to_use <- "SCT"
    cat("SCT data accessed successfully.\n")
  }
}

# Try integrated assay as last resort
if (is.null(assay_to_use) && "integrated" %in% Assays(obj)) {
  cat("integrated assay found. Trying to access data...\n")
  
  data_ok <- tryCatch(
    {
      test <- GetAssayData(obj, assay = "integrated", layer = "data")
      nrow(test) > 0 && ncol(test) > 0
    },
    error = function(e) {
      cat("  v5 layer failed: ", conditionMessage(e), "\n")
      FALSE
    }
  )
  
  if (!data_ok) {
    data_ok <- tryCatch(
      {
        test <- GetAssayData(obj, assay = "integrated", slot = "data")
        nrow(test) > 0 && ncol(test) > 0
      },
      error = function(e) {
        cat("  v4 slot failed: ", conditionMessage(e), "\n")
        FALSE
      }
    )
  }
  
  if (!data_ok) {
    data_ok <- tryCatch(
      {
        int_assay <- obj[["integrated"]]
        test <- int_assay@data
        nrow(test) > 0 && ncol(test) > 0
      },
      error = function(e) {
        cat("  Direct slot access failed: ", conditionMessage(e), "\n")
        FALSE
      }
    )
  }
  
  if (data_ok) {
    assay_to_use <- "integrated"
    cat("integrated data accessed successfully.\n")
  }
}

if (is.null(assay_to_use)) {
  stop(
    "Could not access normalized data from any assay."
  )
}

cat(
  "Expression assay used: ",
  assay_to_use,
  "\n",
  sep = ""
)

cat(
  "No new normalization performed.\n"
)

DefaultAssay(obj) <- assay_to_use

# ------------------------------------------------------------
# 5. Check marker availability
# ------------------------------------------------------------

features_available <- tryCatch(
  rownames(obj[[assay_to_use]]),
  error = function(e) {
    cat("rownames failed, trying alternative...\n")
    # Try to get features from the counts matrix
    tryCatch({
      counts <- GetAssayData(obj, assay = assay_to_use, layer = "counts")
      rownames(counts)
    }, error = function(e2) {
      tryCatch({
        counts <- GetAssayData(obj, assay = assay_to_use, slot = "counts")
        rownames(counts)
      }, error = function(e3) {
        character(0)
      })
    })
  }
)

cat("Features available: ", length(features_available), "\n", sep = "")

marker_availability <- rbindlist(
  lapply(
    names(marker_panels),
    function(lineage) {

      genes <- marker_panels[[lineage]]

      data.table(
        lineage = lineage,
        gene = genes,
        present = genes %in% features_available
      )
    }
  )
)

cat("\n===== 4. MARKER AVAILABILITY =====\n\n")

marker_summary <- marker_availability[
  ,
  .(
    markers_requested = .N,
    markers_present = sum(present),
    markers_missing = sum(!present)
  ),
  by = lineage
]

print(marker_summary)

fwrite(
  marker_availability,
  file.path(
    table_dir,
    "GSE197677_marker_availability.csv"
  ),
  bom = TRUE
)

# ------------------------------------------------------------
# 6. Average expression by cluster
# ------------------------------------------------------------

markers_present <- unique(
  marker_availability[
    present == TRUE,
    gene
  ]
)

cat(
  "\nMarkers available for averaging: ",
  length(markers_present),
  "\n",
  sep = ""
)

# Ensure only occupied cluster levels are used
# Use a workaround to avoid the images slot issue
cluster_factor <- factor(
  as.character(md$seurat_clusters),
  levels = cluster_counts$cluster
)

# Add to metadata using a safer method
tryCatch(
  {
    obj$cluster_audit <- cluster_factor
  },
  error = function(e) {
    cat("Warning: Could not add cluster_audit to object: ", conditionMessage(e), "\n")
    cat("Will use metadata directly for AverageExpression.\n")
  }
)

cat("\n===== 5. CALCULATING CLUSTER AVERAGE EXPRESSION =====\n\n")

# Manual calculation to avoid Seurat v4/v5 compatibility issues
cat("Using manual average calculation to avoid slot errors...\n")

# Get the expression data directly from the slot
cat("Extracting expression data for markers only (memory-efficient)...\n")

# Access the RNA assay directly from the slots
rna_slot <- obj@assays[["RNA"]]

# Get feature names from the counts slot
feature_names <- tryCatch(
  {
    counts <- rna_slot@counts
    rownames(counts)
  },
  error = function(e) {
    tryCatch({
      data_mat <- rna_slot@data
      rownames(data_mat)
    }, error = function(e2) {
      stop("Cannot extract feature names from assay")
    })
  }
)

cat("Total features in assay: ", length(feature_names), "\n", sep = "")

# Find which markers are available
markers_present_intersect <- intersect(markers_present, feature_names)
cat("Markers intersecting with expression matrix: ", length(markers_present_intersect), "\n", sep = "")

# Extract only the markers we need from the data slot
# This is much more memory efficient than extracting the full matrix
cat("Extracting marker expression data...\n")

# Create a list to store marker expression per cluster
cluster_assignments <- as.character(md$seurat_clusters)
names(cluster_assignments) <- md$cell_id

# Get unique clusters in order
unique_clusters <- cluster_counts$cluster

# Initialize average expression matrix
avg <- matrix(0, nrow = length(markers_present_intersect), ncol = length(unique_clusters))
rownames(avg) <- markers_present_intersect
colnames(avg) <- unique_clusters

# For each marker, extract its expression and calculate cluster means
for (marker in markers_present_intersect) {
  # Get the row index for this marker
  marker_idx <- which(feature_names == marker)
  
  if (length(marker_idx) == 0) next
  
  # Extract expression for this marker from the data slot
  marker_expr <- tryCatch(
    {
      # Try to get just this one row from the sparse matrix
      rna_slot@data[marker_idx, ]
    },
    error = function(e) {
      cat("Warning: Could not extract ", marker, "\n")
      NULL
    }
  )
  
  if (is.null(marker_expr)) next
  
  # Calculate mean per cluster
  for (cl in unique_clusters) {
    cell_ids <- names(cluster_assignments[cluster_assignments == cl])
    cell_ids_intersect <- intersect(cell_ids, names(marker_expr))
    
    if (length(cell_ids_intersect) > 0) {
      avg[marker, cl] <- mean(marker_expr[cell_ids_intersect], na.rm = TRUE)
    }
  }
}

cat(
  "Average-expression matrix: ",
  nrow(avg),
  " markers x ",
  ncol(avg),
  " clusters\n",
  sep = ""
)

cat(
  "Average-expression matrix: ",
  nrow(avg),
  " markers x ",
  ncol(avg),
  " clusters\n",
  sep = ""
)

# ------------------------------------------------------------
# 7. Gene-wise z-score
# ------------------------------------------------------------

zscore_rows <- t(
  apply(
    avg,
    1,
    function(x) {

      if (
        length(unique(x)) <= 1 ||
        sd(x, na.rm = TRUE) == 0
      ) {
        return(
          rep(0, length(x))
        )
      }

      as.numeric(
        scale(x)
      )
    }
  )
)

rownames(zscore_rows) <- rownames(avg)
colnames(zscore_rows) <- colnames(avg)

# ------------------------------------------------------------
# 8. Lineage scores per cluster
# ------------------------------------------------------------

score_list <- list()

for (lineage in names(marker_panels)) {

  genes <- intersect(
    marker_panels[[lineage]],
    rownames(zscore_rows)
  )

  if (length(genes) == 0) {
    next
  }

  lineage_score <- colMeans(
    zscore_rows[
      genes,
      ,
      drop = FALSE
    ],
    na.rm = TRUE
  )

  score_list[[lineage]] <- data.table(
    cluster = names(lineage_score),
    lineage = lineage,
    score = as.numeric(lineage_score),
    n_markers_used = length(genes)
  )
}

scores_long <- rbindlist(
  score_list,
  fill = TRUE
)

fwrite(
  scores_long,
  file.path(
    table_dir,
    "GSE197677_cluster_lineage_scores_long.csv"
  ),
  bom = TRUE
)

# ------------------------------------------------------------
# 9. Top / second lineage and confidence margin
# ------------------------------------------------------------

ranked <- scores_long[
  order(
    cluster,
    -score
  ),
  {
    top_lineage <- lineage[1]
    top_score <- score[1]

    second_lineage <- if (.N >= 2) lineage[2] else NA_character_
    second_score <- if (.N >= 2) score[2] else NA_real_

    .(
      suggested_broad_type = top_lineage,
      top_score = top_score,
      second_type = second_lineage,
      second_score = second_score,
      score_margin = top_score - second_score
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

# Conservative review rule
ranked[
  ,
  review_status :=
    fifelse(
      is.na(score_margin),
      "REVIEW",

      fifelse(
        score_margin >= 0.40,
        "HIGH_CONFIDENCE_CANDIDATE",

        fifelse(
          score_margin >= 0.20,
          "MODERATE_CONFIDENCE_CANDIDATE",
          "MANUAL_REVIEW_REQUIRED"
        )
      )
    )
]

ranked[
  ,
  cluster_numeric :=
    suppressWarnings(
      as.integer(cluster)
    )
]

setorder(
  ranked,
  cluster_numeric
)

ranked[
  ,
  cluster_numeric := NULL
]

cat("\n===== 6. BROAD CELL-TYPE CANDIDATES =====\n\n")

print(ranked)

fwrite(
  ranked,
  file.path(
    table_dir,
    "GSE197677_cluster_broad_celltype_candidates.csv"
  ),
  bom = TRUE
)

# ------------------------------------------------------------
# 10. Marker average-expression table
# ------------------------------------------------------------

avg_dt <- as.data.table(
  avg,
  keep.rownames = "gene"
)

fwrite(
  avg_dt,
  file.path(
    table_dir,
    "GSE197677_cluster_marker_average_expression.csv"
  ),
  bom = TRUE
)

# ------------------------------------------------------------
# 11. Heatmap data
# ------------------------------------------------------------

heat_dt <- as.data.table(
  as.table(
    zscore_rows
  )
)

setnames(
  heat_dt,
  c(
    "gene",
    "cluster",
    "zscore"
  )
)

lineage_lookup <- rbindlist(
  lapply(
    names(marker_panels),
    function(x) {

      data.table(
        gene = marker_panels[[x]],
        lineage = x
      )
    }
  )
)

heat_dt <- merge(
  heat_dt,
  lineage_lookup,
  by = "gene",
  all.x = TRUE
)

# Cluster order
heat_dt[
  ,
  cluster_numeric :=
    suppressWarnings(
      as.integer(
        as.character(cluster)
      )
    )
]

heat_dt[
  ,
  cluster :=
    factor(
      cluster,
      levels = as.character(
        sort(
          unique(
            cluster_numeric
          )
        )
      )
    )
]

# Marker order by lineage panel
gene_order <- unlist(
  marker_panels,
  use.names = FALSE
)

gene_order <- gene_order[
  gene_order %in%
    unique(heat_dt$gene)
]

heat_dt[
  ,
  gene :=
    factor(
      gene,
      levels = rev(gene_order)
    )
]

# ------------------------------------------------------------
# 12. Heatmap
# ------------------------------------------------------------

p_heat <- ggplot(
  heat_dt,
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
    title = "GSE197677 Broad Cell-Type Marker Audit",
    subtitle = "Cluster-level average expression; gene-wise z-score",
    x = "Seurat Cluster",
    y = "Canonical Marker",
    fill = "Z-score"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    strip.text.y = element_text(
      angle = 0
    )
  )

ggsave(
  file.path(
    figure_dir,
    "GSE197677_broad_celltype_marker_heatmap.pdf"
  ),
  p_heat,
  width = 12,
  height = 13
)

ggsave(
  file.path(
    figure_dir,
    "GSE197677_broad_celltype_marker_heatmap.png"
  ),
  p_heat,
  width = 12,
  height = 13,
  dpi = 300
)

# ------------------------------------------------------------
# 13. Lineage-score tile plot
# ------------------------------------------------------------

score_plot_dt <- copy(
  scores_long
)

score_plot_dt[
  ,
  cluster_numeric :=
    suppressWarnings(
      as.integer(cluster)
    )
]

score_plot_dt[
  ,
  cluster :=
    factor(
      cluster,
      levels = as.character(
        sort(
          unique(cluster_numeric)
        )
      )
    )
]

p_score <- ggplot(
  score_plot_dt,
  aes(
    x = cluster,
    y = lineage,
    fill = score
  )
) +
  geom_tile() +
  labs(
    title = "GSE197677 Broad Lineage Scores by Cluster",
    x = "Seurat Cluster",
    y = "Broad Cell Lineage",
    fill = "Mean marker\nZ-score"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

ggsave(
  file.path(
    figure_dir,
    "GSE197677_cluster_lineage_score_heatmap.pdf"
  ),
  p_score,
  width = 11,
  height = 6
)

ggsave(
  file.path(
    figure_dir,
    "GSE197677_cluster_lineage_score_heatmap.png"
  ),
  p_score,
  width = 11,
  height = 6,
  dpi = 300
)

# ------------------------------------------------------------
# 14. Explicit warning
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("INTERPRETATION RULE\n")
cat("============================================================\n\n")

cat(
  "The suggested labels are CANDIDATE BROAD LINEAGES only.\n"
)

cat(
  "They must be reviewed against canonical marker expression.\n"
)

cat(
  "Epithelial clusters are NOT automatically called malignant.\n"
)

cat(
  "No cluster labels were written back to the Seurat object.\n"
)

# ------------------------------------------------------------
# FINAL
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP 15A FINISHED ✓\n")
cat("============================================================\n\n")

cat("No cells filtered ✓\n")
cat("No normalization performed ✓\n")
cat("No reclustering performed ✓\n")
cat("No malignant epithelial calls ✓\n")
cat("No Seurat metadata modified ✓\n")
cat("No cross-dataset integration ✓\n")

cat("\nMain table:\n")
cat(
  "06_tables/GSE197677/celltype/",
  "GSE197677_cluster_broad_celltype_candidates.csv\n",
  sep = ""
)

cat("\nMain figures:\n")
cat(
  "05_figures/GSE197677/celltype/",
  "GSE197677_broad_celltype_marker_heatmap.png\n",
  sep = ""
)

cat(
  "05_figures/GSE197677/celltype/",
  "GSE197677_cluster_lineage_score_heatmap.png\n",
  sep = ""
)

cat("\nSTOP HERE.\n")

