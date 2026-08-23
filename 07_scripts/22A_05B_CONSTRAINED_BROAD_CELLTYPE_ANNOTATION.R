# ============================================================
# ESCC Neoadjuvant Project
# STEP22A-5B
#
# CONSTRAINED BROAD CELL-TYPE ANNOTATION
#
# Input:
#   STEP22A-5A v2 preannotation Seurat object
#
# Allowed marker sets ONLY:
#   Fibroblast_CAF
#   Epithelial
#   Immune
#   Endothelial
#   Mural_Pericyte
#
# Explicitly FORBIDDEN for annotation:
#   CDKN1B
#   GSN
#   BGN
#   TIMP1
#   CORE2
#   CORE4
#
# NO treatment-guided annotation.
# NO differential expression marker discovery.
# NO CORE scoring.
# NO before/after inference.
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
cat("STEP22A-5B: CONSTRAINED BROAD CELL-TYPE ANNOTATION\n")
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

INPUT_OBJECT <- file.path(
  OBJECT_DIR,
  "OMIX005710_PAIRED_QC_CLUSTERED_PREANNOTATION_v2_COMMON_ENSG.rds"
)

FEATURE_MAP_FILE <- file.path(
  META_DIR,
  "STEP22A_05A_COMMON_ENSG_FEATURES.csv"
)

GATE_FILE <- file.path(
  META_DIR,
  "STEP22A_ANALYSIS_GATE_STATUS.csv"
)

FROZEN_INVENTORY <- file.path(
  META_DIR,
  "STEP22A_PREANALYSIS_FROZEN_FILE_INVENTORY.csv"
)

# ============================================================
# 1. Analysis-gate protection
# ============================================================

cat("===== 1. ANALYSIS GATE CHECK =====\n\n")

if (!file.exists(GATE_FILE)) {
  stop(
    "Missing analysis gate file: ",
    GATE_FILE
  )
}

gate <- fread(
  GATE_FILE,
  encoding = "UTF-8"
)

gate_allowed <- function(x) {

  out <- gate[
    gate == x,
    allowed
  ]

  if (length(out) != 1) {
    stop(
      "Cannot resolve gate: ",
      x
    )
  }

  out
}

annotation_allowed <- gate_allowed(
  "broad_cell_annotation"
)

core2_allowed <- gate_allowed(
  "paired_CORE2_test"
)

manuscript_allowed <- gate_allowed(
  "manuscript_update"
)

cat(
  "Broad annotation gate :",
  annotation_allowed,
  "\n"
)

cat(
  "Paired CORE2 gate     :",
  core2_allowed,
  "\n"
)

cat(
  "Manuscript gate       :",
  manuscript_allowed,
  "\n\n"
)

if (annotation_allowed != "YES") {
  stop(
    "BROAD_CELL_ANNOTATION_NOT_ALLOWED"
  )
}

if (core2_allowed == "YES") {
  stop(
    "ERROR: CORE2 inference must remain blocked."
  )
}

if (manuscript_allowed == "YES") {
  stop(
    "ERROR: manuscript modification must remain blocked."
  )
}

cat(
  "✓ Broad annotation allowed.\n"
)

cat(
  "✓ CORE2 inference remains blocked.\n\n"
)

# ============================================================
# 2. Load STEP22A-5A v2 object
# ============================================================

cat("===== 2. LOAD PREANNOTATION SEURAT OBJECT =====\n\n")

if (!file.exists(INPUT_OBJECT)) {
  stop(
    "Missing Seurat object:\n",
    INPUT_OBJECT
  )
}

obj <- readRDS(
  INPUT_OBJECT
)

DefaultAssay(
  obj
) <- "RNA"

cat(
  "Cells    :",
  ncol(obj),
  "\n"
)

cat(
  "Features :",
  nrow(obj),
  "\n"
)

cat(
  "Clusters :",
  length(
    unique(
      obj$seurat_clusters
    )
  ),
  "\n\n"
)

if (
  ncol(obj) != 73284
) {

  cat(
    "WARNING: screenshot reported 73,284 retained cells; ",
    "object contains ",
    ncol(obj),
    ".\n\n",
    sep = ""
  )
}

# ============================================================
# 3. Frozen non-CORE marker definition
# ============================================================

cat("===== 3. FREEZE BROAD MARKER SET =====\n\n")

marker_modules <- list(

  Fibroblast_CAF = c(
    "COL1A1",
    "COL1A2",
    "DCN",
    "LUM",
    "COL3A1",
    "COL6A1",
    "COL6A2",
    "COL6A3",
    "DPT",
    "CFD",
    "PDGFRA",
    "FAP"
  ),

  Epithelial = c(
    "EPCAM",
    "KRT8",
    "KRT18",
    "KRT19",
    "KRT5",
    "KRT14"
  ),

  Immune = c(
    "PTPRC",
    "CD3D",
    "CD3E",
    "NKG7",
    "LST1",
    "MS4A1",
    "CD79A"
  ),

  Endothelial = c(
    "PECAM1",
    "VWF",
    "EMCN",
    "KDR"
  ),

  Mural_Pericyte = c(
    "RGS5",
    "CSPG4",
    "MCAM",
    "NOTCH3"
  )
)

forbidden <- c(
  "CDKN1B",
  "GSN",
  "BGN",
  "TIMP1",
  "CORE2",
  "CORE4"
)

all_allowed_markers <- unique(
  unlist(
    marker_modules
  )
)

forbidden_overlap <- intersect(
  toupper(
    all_allowed_markers
  ),
  toupper(
    forbidden
  )
)

if (
  length(
    forbidden_overlap
  ) > 0
) {

  stop(
    "FORBIDDEN_CORE_MARKER_OVERLAP: ",
    paste(
      forbidden_overlap,
      collapse = ", "
    )
  )
}

cat(
  "Frozen broad markers:",
  length(
    all_allowed_markers
  ),
  "\n"
)

cat(
  "Forbidden CORE-related markers used:",
  length(
    forbidden_overlap
  ),
  "\n\n"
)

if (
  length(
    all_allowed_markers
  ) != 33
) {

  stop(
    "Expected exactly 33 frozen broad markers."
  )
}

# ============================================================
# 4. Map symbols -> actual Seurat feature names
#
# Alignment remains ENSG-based from Step22A-5A.
# Here gene symbols are used ONLY to retrieve prespecified
# canonical markers.
#
# If multiple ENSGs share the same symbol:
# use the FIRST entry in the frozen common-feature order.
# This is deterministic and outcome-independent.
# ============================================================

cat("===== 4. MAP FROZEN MARKERS TO COMMON ENSG FEATURES =====\n\n")

if (!file.exists(FEATURE_MAP_FILE)) {

  stop(
    "Missing common feature map:\n",
    FEATURE_MAP_FILE
  )
}

feature_map <- fread(
  FEATURE_MAP_FILE,
  encoding = "UTF-8"
)

needed_cols <- c(
  "common_feature_index",
  "gene_id",
  "gene_symbol",
  "seurat_feature_name"
)

missing_cols <- setdiff(
  needed_cols,
  names(
    feature_map
  )
)

if (
  length(
    missing_cols
  ) > 0
) {

  stop(
    "Feature map missing columns: ",
    paste(
      missing_cols,
      collapse = ", "
    )
  )
}

setorder(
  feature_map,
  common_feature_index
)

marker_map <- list()

for (
  module_name in names(
    marker_modules
  )
) {

  for (
    marker in marker_modules[[
      module_name
    ]]
  ) {

    hits <- feature_map[
      toupper(
        gene_symbol
      ) ==
        toupper(
          marker
        )
    ]

    if (
      nrow(
        hits
      ) == 0
    ) {

      stop(
        "Frozen marker missing from common feature universe: ",
        marker
      )
    }

    selected <- hits[1]

    if (
      !selected$seurat_feature_name
      %in%
      rownames(
        obj
      )
    ) {

      stop(
        "Mapped feature not present in Seurat object: ",
        selected$seurat_feature_name
      )
    }

    marker_map[[
      length(
        marker_map
      ) + 1
    ]] <- data.table(

      module =
        module_name,

      marker =
        marker,

      selected_gene_id =
        selected$gene_id,

      selected_feature =
        selected$seurat_feature_name,

      candidate_ENSG_count =
        nrow(
          hits
        ),

      mapping_rule =
        "FIRST_COMMON_ENSG_ORDER"
    )
  }
}

marker_map <- rbindlist(
  marker_map
)

MARKER_MAP_OUT <- file.path(
  META_DIR,
  "STEP22A_05B_FROZEN_MARKER_MAPPING.csv"
)

fwrite(
  marker_map,
  MARKER_MAP_OUT,
  bom = TRUE
)

cat(
  "Mapped markers:",
  nrow(
    marker_map
  ),
  "/ 33\n"
)

if (
  nrow(
    marker_map
  ) != 33
) {

  stop(
    "MARKER_MAPPING_INCOMPLETE"
  )
}

cat(
  "✓ All 33 frozen markers mapped.\n\n"
)

# ============================================================
# 5. Retrieve normalized expression layer
# ============================================================

cat("===== 5. EXTRACT NORMALIZED MARKER EXPRESSION =====\n\n")

rna_data <- tryCatch(

  {
    LayerData(
      obj,
      assay = "RNA",
      layer = "data"
    )
  },

  error = function(e) {

    GetAssayData(
      obj,
      assay = "RNA",
      slot = "data"
    )
  }
)

if (
  nrow(
    rna_data
  ) == 0
) {

  stop(
    "RNA normalized data layer is empty."
  )
}

selected_features <- marker_map$selected_feature

marker_expr <- rna_data[
  selected_features,
  ,
  drop = FALSE
]

rownames(
  marker_expr
) <- marker_map$marker

cat(
  "Marker matrix:",
  nrow(marker_expr),
  "markers x",
  ncol(marker_expr),
  "cells\n\n"
)

# ============================================================
# 6. Cluster structure
# ============================================================

cat("===== 6. BUILD CLUSTER STRUCTURE =====\n\n")

cluster_id <- as.character(
  obj$seurat_clusters
)

cluster_numeric <- suppressWarnings(
  as.integer(
    unique(
      cluster_id
    )
  )
)

if (
  all(
    !is.na(
      cluster_numeric
    )
  )
) {

  cluster_levels <- as.character(
    sort(
      unique(
        as.integer(
          cluster_id
        )
      )
    )
  )

} else {

  cluster_levels <- sort(
    unique(
      cluster_id
    )
  )
}

cat(
  "Clusters:",
  length(
    cluster_levels
  ),
  "\n\n"
)

# ============================================================
# 7. Constrained cluster-level evidence calculation
#
# Two independent marker summaries:
#
# A. Expression evidence
#    Mean normalized expression across module markers + cells.
#
# B. Detection evidence
#    For each marker:
#       fraction of cells with normalized expression > 0.
#    Module score:
#       mean marker detection fraction.
#
# A cluster receives a broad label ONLY when the SAME module
# ranks #1 for BOTH expression and detection.
#
# Support floor:
#   at least max(2, ceiling(25% of module markers))
#   markers must be detected in >=10% of cluster cells.
#
# Otherwise:
#   Unresolved
#
# No treatment/timepoint variable is used.
# ============================================================

cat("===== 7. CONSTRAINED CLUSTER EVIDENCE =====\n\n")

module_sizes <- vapply(
  marker_modules,
  length,
  integer(1)
)

minimum_support <- pmax(
  2L,
  ceiling(
    module_sizes * 0.25
  )
)

minimum_support <- setNames(
  minimum_support,
  names(
    module_sizes
  )
)

support_detection_fraction <- 0.10

cluster_evidence <- list()
cluster_decisions <- list()

for (
  cl in cluster_levels
) {

  idx <- which(
    cluster_id ==
      cl
  )

  if (
    length(
      idx
    ) == 0
  ) {
    next
  }

  sample_count <- uniqueN(
    obj$sample_name[
      idx
    ]
  )

  patient_count <- uniqueN(
    obj$patient_id[
      idx
    ]
  )

  n_before <- sum(
    obj$timepoint[
      idx
    ] ==
      "Before"
  )

  n_after <- sum(
    obj$timepoint[
      idx
    ] ==
      "After"
  )

  module_rows <- list()

  for (
    module_name in names(
      marker_modules
    )
  ) {

    markers <- marker_modules[[
      module_name
    ]]

    m <- marker_expr[
      markers,
      idx,
      drop = FALSE
    ]

    marker_detection <- Matrix::rowMeans(
      m > 0
    )

    module_mean_expression <- mean(
      m
    )

    module_mean_detection <- mean(
      marker_detection
    )

    support_count <- sum(
      marker_detection >=
        support_detection_fraction
    )

    module_rows[[
      module_name
    ]] <- data.table(

      cluster =
        cl,

      module =
        module_name,

      n_cells =
        length(
          idx
        ),

      mean_normalized_expression =
        as.numeric(
          module_mean_expression
        ),

      mean_marker_detection_fraction =
        as.numeric(
          module_mean_detection
        ),

      supportive_markers =
        as.integer(
          support_count
        ),

      required_supportive_markers =
        as.integer(
          minimum_support[
            module_name
          ]
        )
    )
  }

  module_table <- rbindlist(
    module_rows
  )

  expression_order <- module_table[
    order(
      -mean_normalized_expression,
      module
    )
  ]

  detection_order <- module_table[
    order(
      -mean_marker_detection_fraction,
      module
    )
  ]

  top_expression <- expression_order$module[1]
  top_detection <- detection_order$module[1]

  consensus <- (
    top_expression ==
      top_detection
  )

  top_row <- module_table[
    module ==
      top_expression
  ]

  support_pass <- (
    consensus &&
    top_row$supportive_markers >=
      top_row$required_supportive_markers &&
    top_row$mean_normalized_expression > 0
  )

  if (
    support_pass
  ) {

    assigned_type <- top_expression

    decision_status <- (
      "CONSENSUS_EXPRESSION_DETECTION"
    )

  } else {

    assigned_type <- "Unresolved"

    if (
      !consensus
    ) {

      decision_status <- (
        "EXPRESSION_DETECTION_DISAGREE"
      )

    } else {

      decision_status <- (
        "INSUFFICIENT_MARKER_SUPPORT"
      )
    }
  }

  cluster_decisions[[
    cl
  ]] <- data.table(

    cluster =
      cl,

    n_cells =
      length(
        idx
      ),

    n_samples =
      sample_count,

    n_patients =
      patient_count,

    before_cells =
      n_before,

    after_cells =
      n_after,

    top_expression_module =
      top_expression,

    top_detection_module =
      top_detection,

    top_supportive_markers =
      top_row$supportive_markers,

    required_supportive_markers =
      top_row$required_supportive_markers,

    assigned_broad_type =
      assigned_type,

    decision_status =
      decision_status
  )

  cluster_evidence[[
    cl
  ]] <- module_table
}

cluster_evidence <- rbindlist(
  cluster_evidence
)

cluster_decisions <- rbindlist(
  cluster_decisions
)

cluster_decisions[
  ,
  cluster_order := as.integer(
    cluster
  )
]

setorder(
  cluster_decisions,
  cluster_order
)

cluster_decisions[
  ,
  cluster_order := NULL
]

EVIDENCE_OUT <- file.path(
  META_DIR,
  "STEP22A_05B_CLUSTER_MODULE_EVIDENCE.csv"
)

DECISION_OUT <- file.path(
  META_DIR,
  "STEP22A_05B_CLUSTER_BROAD_ANNOTATION.csv"
)

fwrite(
  cluster_evidence,
  EVIDENCE_OUT,
  bom = TRUE
)

fwrite(
  cluster_decisions,
  DECISION_OUT,
  bom = TRUE
)

print(
  cluster_decisions
)

# ============================================================
# 8. Apply cluster labels to cells
# ============================================================

cat("\n===== 8. APPLY CLUSTER LABELS TO CELLS =====\n\n")

cluster_label_lookup <- setNames(
  cluster_decisions$
    assigned_broad_type,
  cluster_decisions$
    cluster
)

cell_broad_type <- unname(
  cluster_label_lookup[
    cluster_id
  ]
)

if (
  anyNA(
    cell_broad_type
  )
) {

  stop(
    "CELL_BROAD_TYPE_MAPPING_FAILED"
  )
}

obj$STEP22A_broad_celltype <- (
  cell_broad_type
)

broad_counts <- as.data.table(
  table(
    obj$STEP22A_broad_celltype
  )
)

setnames(
  broad_counts,
  c(
    "broad_celltype",
    "cells"
  )
)

broad_counts[
  ,
  fraction :=
    cells /
    sum(
      cells
    )
]

setorder(
  broad_counts,
  -cells
)

cat(
  "Broad cell counts:\n\n"
)

print(
  broad_counts
)

BROAD_COUNT_OUT <- file.path(
  META_DIR,
  "STEP22A_05B_BROAD_CELLTYPE_COUNTS.csv"
)

fwrite(
  broad_counts,
  BROAD_COUNT_OUT,
  bom = TRUE
)

# ============================================================
# 9. Fibroblast sample-level QC
#
# Frozen rule:
# fibroblast-specific pseudobulk requires >=20 fibroblast cells
# per patient-timepoint unless a stricter frozen threshold exists.
#
# This step ONLY evaluates technical eligibility.
# No pseudobulk and NO CORE scoring are performed.
# ============================================================

cat("\n===== 9. FIBROBLAST SAMPLE QC =====\n\n")

cell_dt <- data.table(

  cell_id =
    colnames(
      obj
    ),

  patient_id =
    as.character(
      obj$patient_id
    ),

  timepoint =
    as.character(
      obj$timepoint
    ),

  sample_name =
    as.character(
      obj$sample_name
    ),

  hrs_accession =
    as.character(
      obj$hrs_accession
    ),

  broad_celltype =
    as.character(
      obj$STEP22A_broad_celltype
    )
)

fibro_qc <- cell_dt[
  ,
  .(
    n_total_cells =
      .N,

    n_fibroblasts =
      sum(
        broad_celltype ==
          "Fibroblast_CAF"
      )
  ),
  by = .(
    patient_id,
    timepoint,
    sample_name,
    hrs_accession
  )
]

fibro_qc[
  ,
  fibroblast_fraction :=
    n_fibroblasts /
    n_total_cells
]

fibro_qc[
  ,
  fibroblast_pseudobulk_technical_status :=
    fifelse(
      n_fibroblasts >= 20,
      "PASS_GE20",
      "INSUFFICIENT_FIBROBLAST_CELLS"
    )
]

patient_number <- function(x) {

  as.integer(
    sub(
      "^P",
      "",
      x
    )
  )
}

fibro_qc[
  ,
  patient_order :=
    patient_number(
      patient_id
    )
]

setorder(
  fibro_qc,
  patient_order,
  timepoint
)

fibro_qc[
  ,
  patient_order := NULL
]

FIBRO_QC_OUT <- file.path(
  META_DIR,
  "STEP22A_05B_FIBROBLAST_SAMPLE_QC.csv"
)

fwrite(
  fibro_qc,
  FIBRO_QC_OUT,
  bom = TRUE
)

print(
  fibro_qc
)

# ============================================================
# 10. Paired technical evaluability
# ============================================================

cat("\n===== 10. PAIRED FIBROBLAST TECHNICAL EVALUABILITY =====\n\n")

paired_technical <- fibro_qc[
  ,
  .(
    before_fibroblasts =
      sum(
        n_fibroblasts[
          timepoint ==
            "Before"
        ]
      ),

    after_fibroblasts =
      sum(
        n_fibroblasts[
          timepoint ==
            "After"
        ]
      ),

    before_pass =
      any(
        timepoint ==
          "Before" &
        fibroblast_pseudobulk_technical_status ==
          "PASS_GE20"
      ),

    after_pass =
      any(
        timepoint ==
          "After" &
        fibroblast_pseudobulk_technical_status ==
          "PASS_GE20"
      )
  ),
  by = patient_id
]

paired_technical[
  ,
  technically_evaluable_pair :=
    before_pass &
    after_pass
]

paired_technical[
  ,
  patient_order :=
    patient_number(
      patient_id
    )
]

setorder(
  paired_technical,
  patient_order
)

paired_technical[
  ,
  patient_order := NULL
]

PAIR_QC_OUT <- file.path(
  META_DIR,
  "STEP22A_05B_PAIRED_FIBROBLAST_TECHNICAL_QC.csv"
)

fwrite(
  paired_technical,
  PAIR_QC_OUT,
  bom = TRUE
)

print(
  paired_technical
)

n_technically_evaluable <- sum(
  paired_technical$
    technically_evaluable_pair
)

cat(
  "\nTechnically evaluable paired patients:",
  n_technically_evaluable,
  "/",
  nrow(
    paired_technical
  ),
  "\n"
)

cat(
  "IMPORTANT: histology mapping remains incomplete, ",
  "so this does NOT authorize CORE2 inference.\n\n",
  sep = ""
)

# ============================================================
# 11. Save annotated CELL METADATA only
#
# Do NOT save another 4+ GiB Seurat object.
# ============================================================

cat("===== 11. SAVE ANNOTATION METADATA =====\n\n")

annotated_cell_metadata <- data.table(

  cell_id =
    colnames(
      obj
    ),

  patient_id =
    as.character(
      obj$patient_id
    ),

  timepoint =
    as.character(
      obj$timepoint
    ),

  sample_name =
    as.character(
      obj$sample_name
    ),

  hrs_accession =
    as.character(
      obj$hrs_accession
    ),

  seurat_cluster =
    as.character(
      obj$seurat_clusters
    ),

  broad_celltype =
    as.character(
      obj$STEP22A_broad_celltype
    ),

  nCount_RNA =
    as.numeric(
      obj$nCount_RNA
    ),

  nFeature_RNA =
    as.numeric(
      obj$nFeature_RNA
    ),

  percent_mt =
    as.numeric(
      obj$percent.mt
    )
)

CELL_META_OUT <- file.path(
  META_DIR,
  "STEP22A_05B_ANNOTATED_CELL_METADATA.csv"
)

fwrite(
  annotated_cell_metadata,
  CELL_META_OUT,
  bom = TRUE
)

annotation_bundle <- list(

  marker_modules =
    marker_modules,

  marker_mapping =
    marker_map,

  cluster_evidence =
    cluster_evidence,

  cluster_annotation =
    cluster_decisions,

  fibroblast_sample_qc =
    fibro_qc,

  paired_technical_qc =
    paired_technical,

  rules =
    list(

      treatment_used_for_annotation =
        FALSE,

      FindAllMarkers_used =
        FALSE,

      CORE2_used =
        FALSE,

      CORE4_used =
        FALSE,

      CDKN1B_used =
        FALSE,

      GSN_used =
        FALSE,

      BGN_used =
        FALSE,

      TIMP1_used =
        FALSE,

      cluster_assignment_rule =
        paste0(
          "Same module must rank first for ",
          "both mean normalized expression ",
          "and mean marker detection; ",
          "support >= max(2, ceiling(25% of markers)); ",
          "support marker detection >=10%."
        )
    )
)

BUNDLE_OUT <- file.path(
  OBJECT_DIR,
  "OMIX005710_STEP22A_05B_BROAD_ANNOTATION_BUNDLE.rds"
)

saveRDS(
  annotation_bundle,
  BUNDLE_OUT,
  compress = "gzip"
)

cat(
  "✓ Annotated cell metadata saved.\n"
)

cat(
  "✓ Small annotation bundle saved.\n"
)

cat(
  "✓ No second full Seurat object was created.\n\n"
)

# ============================================================
# 12. Visualization
# ============================================================

cat("===== 12. SAVE ANNOTATION FIGURES =====\n\n")

p_umap <- DimPlot(
  obj,
  reduction = "umap",
  group.by = "STEP22A_broad_celltype",
  label = TRUE,
  repel = TRUE
) +
  labs(
    title = "Constrained Broad Cell-Type Annotation"
  )

ggsave(
  filename = file.path(
    FIGURE_DIR,
    "STEP22A_05B_UMAP_BROAD_CELLTYPE.png"
  ),
  plot = p_umap,
  width = 10,
  height = 8,
  dpi = 180
)

fibro_plot_df <- as.data.frame(
  fibro_qc
)

fibro_plot_df$sample_label <- paste0(
  fibro_plot_df$patient_id,
  "_",
  fibro_plot_df$timepoint
)

p_fibro <- ggplot(
  fibro_plot_df,
  aes(
    x = sample_label,
    y = n_fibroblasts
  )
) +
  geom_col() +
  geom_hline(
    yintercept = 20,
    linetype = "dashed"
  ) +
  labs(
    title = "Fibroblast/CAF Cells per Patient-Timepoint",
    subtitle = "Dashed line indicates the frozen 20-cell technical threshold",
    x = "Patient-Timepoint",
    y = "Fibroblast/CAF Cells"
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
    "STEP22A_05B_FIBROBLAST_CELL_COUNTS.png"
  ),
  plot = p_fibro,
  width = 10,
  height = 6,
  dpi = 180
)

cat(
  "✓ Figures saved.\n\n"
)

# ============================================================
# 13. Anti-circularity audit
# ============================================================

cat("===== 13. ANTI-CIRCULARITY AUDIT =====\n\n")

anti_circularity <- data.table(

  item = c(
    "author_cell_annotation_available",
    "annotation_method",
    "annotation_unit",
    "treatment_used_for_annotation",
    "response_used_for_annotation",
    "FindAllMarkers_used_for_annotation",
    "CORE2_used_for_annotation",
    "CORE4_used_for_annotation",
    "CDKN1B_used_for_annotation",
    "GSN_used_for_annotation",
    "BGN_used_for_annotation",
    "TIMP1_used_for_annotation",
    "CORE2_scoring_performed",
    "CORE2_inference_performed",
    "histology_mapping_complete"
  ),

  status = c(
    "NO",
    "CONSTRAINED_FROZEN_MARKER_MODULES",
    "UNSUPERVISED_CLUSTER",
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
    "NO",
    "NO"
  )
)

ANTI_OUT <- file.path(
  META_DIR,
  "STEP22A_05B_ANTI_CIRCULARITY_AUDIT.csv"
)

fwrite(
  anti_circularity,
  ANTI_OUT,
  bom = TRUE
)

print(
  anti_circularity
)

# ============================================================
# 14. Frozen Step19-21 integrity check
# ============================================================

cat("\n===== 14. FROZEN PROJECT CHECK =====\n\n")

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

    if (
      !file.exists(
        path
      )
    ) {

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

    expected <- suppressWarnings(
      as.numeric(
        frozen$size_bytes[i]
      )
    )

    actual <- file.info(
      path
    )$size

    if (
      is.finite(
        expected
      ) &&
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
# 15. Session info
# ============================================================

sink(
  file.path(
    LOG_DIR,
    "STEP22A_05B_sessionInfo.txt"
  )
)

print(
  sessionInfo()
)

sink()

# ============================================================
# 16. Final summary
# ============================================================

fibro_cells <- broad_counts[
  broad_celltype ==
    "Fibroblast_CAF",
  cells
]

if (
  length(
    fibro_cells
  ) == 0
) {
  fibro_cells <- 0
}

fibro_clusters <- cluster_decisions[
  assigned_broad_type ==
    "Fibroblast_CAF",
  cluster
]

unresolved_cells <- broad_counts[
  broad_celltype ==
    "Unresolved",
  cells
]

if (
  length(
    unresolved_cells
  ) == 0
) {
  unresolved_cells <- 0
}

summary_out <- file.path(
  META_DIR,
  "STEP22A_05B_BROAD_ANNOTATION_SUMMARY.txt"
)

summary_lines <- c(

  "STEP22A-5B CONSTRAINED BROAD CELL-TYPE ANNOTATION",

  paste0(
    "Total QC-retained cells: ",
    ncol(obj)
  ),

  paste0(
    "Preannotation clusters: ",
    nrow(
      cluster_decisions
    )
  ),

  paste0(
    "Fibroblast_CAF clusters: ",
    length(
      fibro_clusters
    )
  ),

  paste0(
    "Fibroblast_CAF cells: ",
    fibro_cells
  ),

  paste0(
    "Unresolved cells: ",
    unresolved_cells
  ),

  paste0(
    "Technically evaluable paired patients (fibroblasts >=20 both timepoints): ",
    n_technically_evaluable,
    "/",
    nrow(
      paired_technical
    )
  ),

  "Author cell-level annotation available: NO",

  "Annotation method: constrained frozen marker modules",

  "Treatment used for annotation: NO",

  "FindAllMarkers used for annotation: NO",

  "CORE2 used for annotation: NO",

  "CORE2 scoring performed: NO",

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
cat("STEP22A-5B COMPLETE\n")
cat("====================================================================\n\n")

cat(
  "QC-RETAINED CELLS            :",
  ncol(obj),
  "\n"
)

cat(
  "PREANNOTATION CLUSTERS       :",
  nrow(
    cluster_decisions
  ),
  "\n\n"
)

cat(
  "BROAD CELL-TYPE COUNTS:\n"
)

print(
  broad_counts
)

cat("\n")

cat(
  "FIBROBLAST_CAF CLUSTERS      :",
  length(
    fibro_clusters
  ),
  "\n"
)

cat(
  "FIBROBLAST_CAF CLUSTER IDs   :",
  ifelse(
    length(
      fibro_clusters
    ) > 0,
    paste(
      fibro_clusters,
      collapse = ", "
    ),
    "<none>"
  ),
  "\n"
)

cat(
  "FIBROBLAST_CAF CELLS         :",
  fibro_cells,
  "\n\n"
)

cat(
  "PAIRED PATIENTS >=20 FIBROBLASTS BOTH TIMEPOINTS :",
  n_technically_evaluable,
  "/",
  nrow(
    paired_technical
  ),
  "\n\n"
)

cat(
  "ANNOTATION METHOD            : CONSTRAINED FROZEN MARKERS\n"
)

cat(
  "TREATMENT USED               : NO\n"
)

cat(
  "FindAllMarkers USED          : NO\n"
)

cat(
  "CORE2/CORE4 USED             : NO\n"
)

cat(
  "CDKN1B/GSN/BGN/TIMP1 USED    : NO\n\n"
)

cat(
  "CORE2 SCORING                : NOT PERFORMED\n"
)

cat(
  "CORE2 PRIMARY INFERENCE      : BLOCKED\n"
)

cat(
  "HISTOLOGY MAPPING            : INCOMPLETE\n\n"
)

cat(
  "FULL SEURAT COPY CREATED     : NO\n"
)

cat(
  "FROZEN STEP19-21 MODIFIED    : NO\n"
)

cat(
  "MANUSCRIPT MODIFIED          : NO\n\n"
)

cat(
  "NEXT REQUIRED STEP:\n"
)

if (
  fibro_cells == 0
) {

  cat(
    "STEP22A TECHNICAL VALIDATION NOT EVALUABLE:\n"
  )

  cat(
    "NO FIBROBLAST_CAF COMPARTMENT IDENTIFIED\n"
  )

} else {

  cat(
    "STEP22A-5C = FIBROBLAST PSEUDOBULK TECHNICAL PREPARATION\n"
  )

  cat(
    "NO CORE2 SCORING YET\n"
  )
}

cat("\nSTOP HERE.\n\n")
