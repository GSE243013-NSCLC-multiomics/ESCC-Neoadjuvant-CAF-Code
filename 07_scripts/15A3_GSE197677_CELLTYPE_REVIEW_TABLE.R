# ============================================================
# ESCC Neoadjuvant Project
# STEP 15A3
# Review GSE197677 broad cell-type candidates
#
# Uses outputs from STEP 15A2 only.
# Does NOT load the large Seurat RDS.
#
# NO metadata write-back
# NO filtering
# NO normalization
# NO integration
# ============================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
})

cat("\n============================================================\n")
cat("STEP 15A3: GSE197677 CELL-TYPE REVIEW TABLE\n")
cat("============================================================\n\n")

indir <- "06_tables/GSE197677/celltype"

candidate_file <- file.path(
  indir,
  "GSE197677_cluster_broad_celltype_candidates_FIXED.csv"
)

top3_file <- file.path(
  indir,
  "GSE197677_cluster_top3_lineages.csv"
)

pb_file <- file.path(
  indir,
  "GSE197677_rawcount_marker_pseudobulk.csv"
)

detect_file <- file.path(
  indir,
  "GSE197677_marker_detection_fraction.csv"
)

for (f in c(
  candidate_file,
  top3_file,
  pb_file,
  detect_file
)) {

  if (!file.exists(f)) {
    stop("Missing file: ", f)
  }
}

# ------------------------------------------------------------
# 1. Load STEP 15A2 results
# ------------------------------------------------------------

cand <- fread(candidate_file)
top3 <- fread(top3_file)
pb <- fread(pb_file)
det <- fread(detect_file)

cat("Candidate clusters:", nrow(cand), "\n")

if (nrow(cand) != 25) {
  cat("WARNING: expected 25 occupied clusters.\n")
}

# ------------------------------------------------------------
# 2. Canonical panels
# ------------------------------------------------------------

panels <- list(

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
    "C1QA",
    "C1QC",
    "S100A8",
    "S100A9"
  ),

  B_plasma = c(
    "CD79A",
    "MS4A1",
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
    "PLVAP"
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

gene_to_lineage <- rbindlist(
  lapply(
    names(panels),
    function(x) {
      data.table(
        gene = panels[[x]],
        lineage = x
      )
    }
  )
)

# ------------------------------------------------------------
# 3. Add lineage to marker pseudobulk
# ------------------------------------------------------------

pb2 <- merge(
  pb,
  gene_to_lineage,
  by = "gene",
  all.x = TRUE
)

det2 <- merge(
  det,
  gene_to_lineage,
  by = "gene",
  all.x = TRUE
)

# ------------------------------------------------------------
# 4. Helper: strongest markers within cluster
# ------------------------------------------------------------

get_top_markers <- function(
  cluster_id,
  lineage_id,
  n = 5
) {

  x <- pb2[
    cluster == cluster_id &
    lineage == lineage_id
  ]

  if (nrow(x) == 0) {
    return(NA_character_)
  }

  setorder(
    x,
    -zscore,
    -log1p_CPM
  )

  genes <- head(
    x$gene,
    n
  )

  paste(
    genes,
    collapse = ", "
  )
}

get_detection_markers <- function(
  cluster_id,
  lineage_id,
  n = 5
) {

  x <- det2[
    cluster == cluster_id &
    lineage == lineage_id
  ]

  if (nrow(x) == 0) {
    return(NA_character_)
  }

  setorder(
    x,
    -detection_fraction
  )

  x <- head(x, n)

  paste0(
    x$gene,
    "=",
    sprintf(
      "%.0f%%",
      100 * x$detection_fraction
    ),
    collapse = ", "
  )
}

# ------------------------------------------------------------
# 5. Top-3 lineage columns
# ------------------------------------------------------------

top3[
  ,
  rank :=
    seq_len(.N),
  by = cluster
]

top1 <- top3[
  rank == 1,
  .(
    cluster,
    lineage1 = lineage,
    score1 = lineage_score,
    detection1 = mean_detection_fraction
  )
]

top2 <- top3[
  rank == 2,
  .(
    cluster,
    lineage2 = lineage,
    score2 = lineage_score,
    detection2 = mean_detection_fraction
  )
]

top3b <- top3[
  rank == 3,
  .(
    cluster,
    lineage3 = lineage,
    score3 = lineage_score,
    detection3 = mean_detection_fraction
  )
]

review <- Reduce(
  function(x, y) {
    merge(
      x,
      y,
      by = "cluster",
      all.x = TRUE
    )
  },
  list(
    cand,
    top1,
    top2,
    top3b
  )
)

# ------------------------------------------------------------
# 6. Marker evidence for top lineage
# ------------------------------------------------------------

review[
  ,
  top_marker_zscore :=
    mapply(
      get_top_markers,
      cluster,
      lineage1,
      SIMPLIFY = TRUE
    )
]

review[
  ,
  top_marker_detection :=
    mapply(
      get_detection_markers,
      cluster,
      lineage1,
      SIMPLIFY = TRUE
    )
]

# ------------------------------------------------------------
# 7. Conservative recommendation
# ------------------------------------------------------------

review[
  ,
  audit_decision :=
    fcase(

      !is.na(score_margin) &
      score_margin >= 0.50 &
      top_detection >= 0.10,
      "LIKELY_ACCEPT",

      !is.na(score_margin) &
      score_margin >= 0.25,
      "REVIEW_HEATMAP",

      default =
        "MANUAL_REVIEW"
    )
]

# Epithelial is still epithelial only.
# DO NOT call malignant here.
review[
  suggested_broad_type == "Epithelial",
  audit_note :=
    "Broad epithelial only; malignant status NOT assigned"
]

review[
  suggested_broad_type != "Epithelial",
  audit_note :=
    "Broad lineage candidate only"
]

review[
  ,
  cluster_number :=
    as.integer(cluster)
]

setorder(
  review,
  cluster_number
)

review[
  ,
  cluster_number := NULL
]

# ------------------------------------------------------------
# 8. Print compact review table
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("CLUSTER REVIEW TABLE\n")
cat("============================================================\n\n")

compact <- review[
  ,
  .(
    cluster,
    n_cells,
    lineage1,
    score1 = round(score1, 3),
    detection1 = round(detection1, 3),

    lineage2,
    score2 = round(score2, 3),

    score_margin = round(score_margin, 3),

    audit_decision,

    top_marker_zscore,

    top_marker_detection
  )
]

print(
  compact,
  nrows = 100
)

# ------------------------------------------------------------
# 9. Flag ambiguous clusters
# ------------------------------------------------------------

ambiguous <- review[
  audit_decision != "LIKELY_ACCEPT"
]

cat("\n============================================================\n")
cat("AMBIGUOUS / REVIEW-REQUIRED CLUSTERS\n")
cat("============================================================\n\n")

if (nrow(ambiguous) == 0) {

  cat(
    "No ambiguous clusters under current rule.\n"
  )

} else {

  print(
    ambiguous[
      ,
      .(
        cluster,
        n_cells,
        lineage1,
        score1,
        lineage2,
        score2,
        score_margin,
        audit_decision,
        top_marker_zscore
      )
    ]
  )
}

# ------------------------------------------------------------
# 10. Lineage cell totals using candidate labels
#
# AUDIT ONLY; not written to Seurat.
# ------------------------------------------------------------

lineage_totals <- review[
  ,
  .(
    clusters = .N,
    cells = sum(n_cells)
  ),
  by = suggested_broad_type
]

setorder(
  lineage_totals,
  -cells
)

cat("\n============================================================\n")
cat("CANDIDATE BROAD-LINEAGE CELL TOTALS\n")
cat("============================================================\n\n")

print(
  lineage_totals
)

# ------------------------------------------------------------
# 11. Save
# ------------------------------------------------------------

fwrite(
  review,
  file.path(
    indir,
    "GSE197677_cluster_manual_review_table.csv"
  ),
  bom = TRUE
)

fwrite(
  compact,
  file.path(
    indir,
    "GSE197677_cluster_manual_review_compact.csv"
  ),
  bom = TRUE
)

fwrite(
  ambiguous,
  file.path(
    indir,
    "GSE197677_clusters_requiring_manual_review.csv"
  ),
  bom = TRUE
)

fwrite(
  lineage_totals,
  file.path(
    indir,
    "GSE197677_candidate_broad_lineage_totals.csv"
  ),
  bom = TRUE
)

cat("\n============================================================\n")
cat("STEP 15A3 FINISHED ✓\n")
cat("============================================================\n\n")

cat(
  "This step only organized evidence for review.\n\n"
)

cat("No Seurat object loaded ✓\n")
cat("No Seurat metadata modified ✓\n")
cat("No cell-type labels written back ✓\n")
cat("No malignant epithelial calls ✓\n")
cat("No integration ✓\n")

cat("\nMain review file:\n")
cat(
  "06_tables/GSE197677/celltype/",
  "GSE197677_cluster_manual_review_compact.csv\n",
  sep = ""
)

cat("\nSTOP HERE.\n")
