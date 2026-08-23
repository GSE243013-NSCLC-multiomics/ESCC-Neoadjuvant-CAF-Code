# ============================================================
# ESCC Neoadjuvant Project
# STEP 17G1
#
# GSE221561 CNV preflight
#
# PURPOSE:
#   1. audit epithelial sample capacity
#   2. select independent untreated technical pilot
#   3. audit RNA count layers
#   4. prove sample-wise raw counts are accessible
#   5. audit mixed SYMBOL / ENSG feature IDs against hg38
#   6. quantify duplicate hg38 target symbols
#
# IMPORTANT:
#   NO InferCNA run
#   NO JoinLayers
#   NO feature aggregation
#   NO gene renaming
#   NO malignant labels
#   NO cell filtering
# ============================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(data.table)
  library(infercna)
})

cat("\n============================================================\n")
cat("STEP 17G1: GSE221561 CNV GENE/LAYER PREFLIGHT\n")
cat("============================================================\n\n")

# ============================================================
# 1. PATHS
# ============================================================

object_dir <- "03_objects/GSE221561"

meta_file <- paste0(
  "06_tables/GSE221561/metadata/",
  "GSE221561_cell_metadata_standardized.csv.gz"
)

outdir <- paste0(
  "06_tables/cross_dataset/",
  "malignant_epithelial/",
  "GSE221561_preflight"
)

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)

if (!file.exists(meta_file)) {
  stop(
    "Missing metadata file: ",
    meta_file
  )
}

# ============================================================
# 2. LOCATE GSE221561 OBJECT
# ============================================================

cat("===== 1. LOCATE GSE221561 OBJECT =====\n\n")

rds_files <- list.files(
  object_dir,
  pattern = "\\.rds$",
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)

if (length(rds_files) == 0) {
  stop(
    "No RDS objects found under ",
    object_dir
  )
}

priority <- rep(
  99L,
  length(rds_files)
)

priority[
  grepl(
    "author.*annotated|annotated.*author",
    basename(rds_files),
    ignore.case = TRUE
  )
] <- 1L

priority[
  grepl(
    "counts.*merged|merged.*counts",
    basename(rds_files),
    ignore.case = TRUE
  )
] <- pmin(
  priority[
    grepl(
      "counts.*merged|merged.*counts",
      basename(rds_files),
      ignore.case = TRUE
    )
  ],
  2L
)

rds_files <- rds_files[
  order(
    priority,
    file.info(rds_files)$size
  )
]

cat("Candidate RDS files:\n")

for (f in rds_files) {
  cat(
    sprintf(
      "  %.2f GB  %s\n",
      file.info(f)$size / 1024^3,
      f
    )
  )
}

object_file <- rds_files[1]

cat(
  "\nSelected object: ",
  object_file,
  "\n",
  sep = ""
)

# ============================================================
# 3. LOAD STANDARDIZED CELL METADATA
# ============================================================

cat("\n===== 2. LOAD STANDARDIZED METADATA =====\n\n")

md <- fread(
  meta_file
)

required_md <- c(
  "cell_id",
  "patient_id",
  "sample_id",
  "tissue_std",
  "treatment_group",
  "cell_class_std"
)

missing_md <- setdiff(
  required_md,
  names(md)
)

if (length(missing_md) > 0) {

  cat("Available columns:\n")
  print(names(md))

  stop(
    "Missing metadata columns: ",
    paste(
      missing_md,
      collapse = ", "
    )
  )
}

cat(
  "Metadata cells: ",
  nrow(md),
  "\n",
  sep = ""
)

md[
  ,
  tissue_simple :=
    fcase(

      grepl(
        "tumor|tumour|scc",
        tissue_std,
        ignore.case = TRUE
      ),
      "Tumor",

      grepl(
        "normal|adjacent",
        tissue_std,
        ignore.case = TRUE
      ),
      "Normal",

      default =
        as.character(tissue_std)
    )
]

# ============================================================
# 4. EPITHELIAL INVENTORY
# ============================================================

cat("\n============================================================\n")
cat("3. GSE221561 EPITHELIAL SAMPLE INVENTORY\n")
cat("============================================================\n\n")

epi_md <- md[
  cell_class_std ==
    "Epithelial_unclassified"
]

cat(
  "Total broad epithelial cells: ",
  nrow(epi_md),
  "\n",
  sep = ""
)

epi_inventory <- epi_md[
  ,
  .(
    epithelial_cells = .N
  ),
  by = .(
    patient_id,
    sample_id,
    tissue_simple,
    treatment_group
  )
]

setorder(
  epi_inventory,
  tissue_simple,
  treatment_group,
  -epithelial_cells
)

print(
  epi_inventory,
  nrows = 100
)

fwrite(
  epi_inventory,
  file.path(
    outdir,
    "STEP17G1_GSE221561_epithelial_sample_inventory.csv"
  ),
  bom = TRUE
)

# ============================================================
# 5. NORMAL EPITHELIAL REFERENCES
# ============================================================

cat("\n============================================================\n")
cat("4. NORMAL EPITHELIAL REFERENCES\n")
cat("============================================================\n\n")

normal_inventory <- epi_inventory[
  tissue_simple == "Normal"
]

print(
  normal_inventory,
  nrows = 100
)

cat(
  "\nNormal epithelial samples: ",
  nrow(normal_inventory),
  "\n",
  sep = ""
)

cat(
  "Total normal epithelial cells: ",
  sum(
    normal_inventory$epithelial_cells
  ),
  "\n",
  sep = ""
)

if (nrow(normal_inventory) < 2) {
  stop(
    "Fewer than 2 independent normal epithelial samples."
  )
}

# ============================================================
# 6. SELECT TECHNICAL PILOT
# ============================================================

cat("\n============================================================\n")
cat("5. TUMOR PILOT CANDIDATES\n")
cat("============================================================\n\n")

tumor_inventory <- epi_inventory[
  tissue_simple == "Tumor"
]

setorder(
  tumor_inventory,
  -epithelial_cells
)

print(
  tumor_inventory,
  nrows = 100
)

surgery_candidates <- tumor_inventory[
  grepl(
    "surgery",
    treatment_group,
    ignore.case = TRUE
  )
]

if (nrow(surgery_candidates) > 0) {

  setorder(
    surgery_candidates,
    -epithelial_cells
  )

  pilot_row <-
    surgery_candidates[1]

  selection_rule <-
    "largest Surgery_alone epithelial sample"

} else {

  pilot_row <-
    tumor_inventory[1]

  selection_rule <-
    "largest tumor epithelial sample; no Surgery_alone candidate detected"
}

pilot_sample <-
  pilot_row$sample_id

pilot_patient <-
  pilot_row$patient_id

pilot_epi_cells <-
  pilot_row$epithelial_cells

pilot_treatment <-
  pilot_row$treatment_group

cat(
  "\nSelected pilot sample   : ",
  pilot_sample,
  "\n",
  sep = ""
)

cat(
  "Patient                 : ",
  pilot_patient,
  "\n",
  sep = ""
)

cat(
  "Treatment               : ",
  pilot_treatment,
  "\n",
  sep = ""
)

cat(
  "Tumor epithelial cells  : ",
  pilot_epi_cells,
  "\n",
  sep = ""
)

cat(
  "Selection rule          : ",
  selection_rule,
  "\n",
  sep = ""
)

# ============================================================
# 7. LOAD SEURAT OBJECT
# ============================================================

cat("\n===== 6. LOAD SEURAT OBJECT =====\n\n")

cat(
  "Loading object; this may take several minutes...\n"
)

obj <- readRDS(
  object_file
)

cat(
  "Object class : ",
  paste(
    class(obj),
    collapse = ", "
  ),
  "\n",
  sep = ""
)

cat(
  "Object cells : ",
  ncol(obj),
  "\n",
  sep = ""
)

cat(
  "Assays       : ",
  paste(
    Seurat::Assays(obj),
    collapse = ", "
  ),
  "\n",
  sep = ""
)

if (
  !"RNA" %in%
    Seurat::Assays(obj)
) {
  stop(
    "RNA assay not found."
  )
}

rna <- obj@assays[["RNA"]]

cat(
  "RNA class    : ",
  paste(
    class(rna),
    collapse = ", "
  ),
  "\n",
  sep = ""
)

cat(
  "RNA features : ",
  nrow(rna),
  "\n",
  sep = ""
)

# ============================================================
# 8. METADATA / OBJECT CELL MATCH
# ============================================================

cat("\n===== 7. CELL-ID MATCH =====\n\n")

object_cells <- colnames(obj)

metadata_in_object <- sum(
  md$cell_id %in%
    object_cells
)

object_in_metadata <- sum(
  object_cells %in%
    md$cell_id
)

cat(
  "Metadata cells found in object: ",
  metadata_in_object,
  " / ",
  nrow(md),
  "\n",
  sep = ""
)

cat(
  "Object cells found in metadata: ",
  object_in_metadata,
  " / ",
  length(object_cells),
  "\n",
  sep = ""
)

if (
  metadata_in_object <
    nrow(md)
) {

  missing_meta_cells <- md[
    !cell_id %in%
      object_cells,
    cell_id
  ]

  cat(
    "\nMetadata cells absent from object: ",
    length(missing_meta_cells),
    "\n",
    sep = ""
  )

  cat("Examples:\n")
  print(
    head(
      missing_meta_cells,
      20
    )
  )
}

pilot_cells <- epi_md[
  sample_id == pilot_sample &
  tissue_simple == "Tumor",
  cell_id
]

normal_cells <- epi_md[
  tissue_simple == "Normal",
  cell_id
]

missing_pilot <- setdiff(
  pilot_cells,
  object_cells
)

missing_normal <- setdiff(
  normal_cells,
  object_cells
)

cat(
  "\nPilot epithelial cells absent from object : ",
  length(missing_pilot),
  "\n",
  sep = ""
)

cat(
  "Normal epithelial cells absent from object: ",
  length(missing_normal),
  "\n",
  sep = ""
)

if (
  length(missing_pilot) > 0 ||
  length(missing_normal) > 0
) {
  stop(
    "Selected pilot/reference cells are not fully represented in object."
  )
}

# ============================================================
# 9. COUNT LAYERS
# ============================================================

cat("\n============================================================\n")
cat("8. RNA COUNT LAYER AUDIT\n")
cat("============================================================\n\n")

layers <- tryCatch(
  SeuratObject::Layers(
    rna,
    search = NA
  ),
  error = function(e) {
    character(0)
  }
)

cat("All RNA layers:\n")
print(layers)

count_layers <- layers[
  grepl(
    "^counts($|\\.)",
    layers,
    ignore.case = TRUE
  )
]

cat(
  "\nCount layers found: ",
  length(count_layers),
  "\n",
  sep = ""
)

print(count_layers)

if (length(count_layers) == 0) {
  stop(
    "No RNA counts layer detected."
  )
}

# ============================================================
# 10. MAP CELLS TO COUNT LAYERS
# ============================================================

cat("\n===== 9. MAP CELLS TO COUNT LAYERS =====\n\n")

get_layer_cells <- function(
  assay_object,
  layer_name
) {

  x <- tryCatch(
    SeuratObject::Cells(
      assay_object,
      layer = layer_name
    ),
    error = function(e) NULL
  )

  if (!is.null(x)) {
    return(x)
  }

  m <- SeuratObject::LayerData(
    assay_object,
    layer = layer_name,
    fast = FALSE
  )

  colnames(m)
}

layer_membership <- rbindlist(
  lapply(
    count_layers,
    function(lyr) {

      lc <- get_layer_cells(
        rna,
        lyr
      )

      cat(
        sprintf(
          "%-40s %7d cells\n",
          lyr,
          length(lc)
        )
      )

      data.table(
        cell_id = lc,
        count_layer = lyr
      )
    }
  ),
  use.names = TRUE
)

dup_layer_membership <- layer_membership[
  ,
  .N,
  by = cell_id
][
  N > 1
]

cat(
  "\nCells assigned to >1 count layer: ",
  nrow(dup_layer_membership),
  "\n",
  sep = ""
)

epi_layer <- merge(
  epi_md[
    ,
    .(
      cell_id,
      patient_id,
      sample_id,
      tissue_simple,
      treatment_group
    )
  ],
  layer_membership,
  by = "cell_id",
  all.x = TRUE
)

unmapped_epi_layers <- sum(
  is.na(
    epi_layer$count_layer
  )
)

cat(
  "Epithelial cells without count layer: ",
  unmapped_epi_layers,
  "\n",
  sep = ""
)

if (unmapped_epi_layers > 0) {
  stop(
    "Some epithelial cells do not map to a counts layer."
  )
}

layer_sample <- epi_layer[
  ,
  .(
    epithelial_cells = .N
  ),
  by = .(
    sample_id,
    tissue_simple,
    treatment_group,
    count_layer
  )
]

setorder(
  layer_sample,
  tissue_simple,
  sample_id,
  count_layer
)

cat("\nCount layer x epithelial sample:\n\n")

print(
  layer_sample,
  nrows = 200
)

fwrite(
  layer_sample,
  file.path(
    outdir,
    "STEP17G1_GSE221561_epithelial_count_layer_map.csv"
  ),
  bom = TRUE
)

# ============================================================
# 11. SAMPLE-WISE RAW COUNT EXTRACTION TEST
# ============================================================

cat("\n============================================================\n")
cat("10. SAMPLE-WISE RAW COUNT EXTRACTION TEST\n")
cat("============================================================\n\n")

set.seed(170221561)

normal_samples <- sort(
  unique(
    normal_inventory$sample_id
  )
)

test_cells <- head(
  pilot_cells,
  50
)

for (s in normal_samples) {

  cc <- epi_md[
    sample_id == s &
    tissue_simple == "Normal",
    cell_id
  ]

  test_cells <- c(
    test_cells,
    head(
      cc,
      25
    )
  )
}

test_cells <- unique(
  test_cells
)

test_membership <- layer_membership[
  cell_id %in%
    test_cells
]

test_parts <- list()

for (lyr in unique(test_membership$count_layer)) {

  lyr_test_cells <- test_membership[
    count_layer == lyr,
    cell_id
  ]

  cat(
    "Extracting ",
    length(lyr_test_cells),
    " test cells from ",
    lyr,
    "...\n",
    sep = ""
  )

  mm <- SeuratObject::LayerData(
    rna,
    layer = lyr,
    cells = lyr_test_cells,
    fast = FALSE
  )

  cat(
    "  ",
    nrow(mm),
    " features x ",
    ncol(mm),
    " cells; nnzero=",
    Matrix::nnzero(mm),
    "\n",
    sep = ""
  )

  if (
    ncol(mm) == 0 ||
    Matrix::nnzero(mm) == 0
  ) {
    stop(
      "Raw-count extraction returned an empty matrix."
    )
  }

  test_parts[[lyr]] <- mm
}

cat(
  "\nSAMPLE-WISE RAW COUNTS ACCESSIBLE WITHOUT JoinLayers\n"
)

rm(test_parts)
gc(verbose = FALSE)

# ============================================================
# 12. hg38 GENE-ID AUDIT
# ============================================================

cat("\n============================================================\n")
cat("11. hg38 FEATURE-ID AUDIT\n")
cat("============================================================\n\n")

infercna::useGenome(
  "hg38"
)

genome <- infercna::retrieveGenome(
  name = "hg38"
)

input_gene <- rownames(
  rna
)

if (is.null(input_gene)) {
  stop(
    "RNA feature names are unavailable."
  )
}

cat(
  "Input RNA features: ",
  length(input_gene),
  "\n",
  sep = ""
)

ensembl_clean <- sub(
  "\\.[0-9]+$",
  "",
  input_gene
)

is_ensg <- grepl(
  "^ENSG[0-9]+$",
  ensembl_clean,
  ignore.case = TRUE
)

symbol_match_index <- match(
  input_gene,
  as.character(
    genome$symbol
  )
)

ensembl_match_index <- match(
  ensembl_clean,
  as.character(
    genome$ensembl_gene_id
  )
)

gene_map <- data.table(
  input_feature =
    input_gene,

  ENSG_like =
    is_ensg,

  symbol_match =
    !is.na(
      symbol_match_index
    ),

  ensembl_match =
    !is.na(
      ensembl_match_index
    )
)

gene_map[
  ,
  mapping_type :=
    fcase(

      symbol_match,
      "EXACT_SYMBOL",

      !symbol_match &
      ensembl_match,
      "ENSEMBL_TO_SYMBOL",

      default =
        "UNMAPPED"
    )
]

gene_map[, hg38_symbol := NA_character_]
gene_map[, hg38_ensembl := NA_character_]

sym_idx <- which(gene_map$mapping_type == "EXACT_SYMBOL")
ens_idx <- which(gene_map$mapping_type == "ENSEMBL_TO_SYMBOL")

if (length(sym_idx) > 0) {
  gene_map[sym_idx, hg38_symbol := as.character(genome$symbol[symbol_match_index[sym_idx]])]
  gene_map[sym_idx, hg38_ensembl := as.character(genome$ensembl_gene_id[symbol_match_index[sym_idx]])]
}

if (length(ens_idx) > 0) {
  gene_map[ens_idx, hg38_symbol := as.character(genome$symbol[ensembl_match_index[ens_idx]])]
  gene_map[ens_idx, hg38_ensembl := as.character(genome$ensembl_gene_id[ensembl_match_index[ens_idx]])]
}

mapping_summary <- gene_map[
  ,
  .(
    features = .N
  ),
  by = mapping_type
]

mapping_summary[
  ,
  percent :=
    100 *
    features /
    sum(features)
]

cat("Mapping summary:\n\n")
print(mapping_summary)

mapped <- gene_map[
  mapping_type !=
    "UNMAPPED" &
  !is.na(hg38_symbol)
]

unique_target_symbols <- uniqueN(
  mapped$hg38_symbol
)

target_dup <- mapped[
  ,
  .(
    input_rows = .N,

    input_features =
      paste(
        head(
          input_feature,
          10
        ),
        collapse = "; "
      ),

    mapping_types =
      paste(
        sort(
          unique(
            mapping_type
          )
        ),
        collapse = "; "
      )
  ),
  by = hg38_symbol
][
  input_rows > 1
]

setorder(
  target_dup,
  -input_rows,
  hg38_symbol
)

cat(
  "\nMapped input rows       : ",
  nrow(mapped),
  "\n",
  sep = ""
)

cat(
  "Unique hg38 symbols    : ",
  unique_target_symbols,
  "\n",
  sep = ""
)

cat(
  "Target symbols with >1 input row: ",
  nrow(target_dup),
  "\n",
  sep = ""
)

cat(
  "Input rows involved in target duplicates: ",
  sum(target_dup$input_rows),
  "\n",
  sep = ""
)

if (nrow(target_dup) > 0) {

  cat(
    "\nExamples of duplicate target symbols:\n\n"
  )

  print(
    head(
      target_dup,
      30
    )
  )
}

fwrite(
  gene_map,
  file.path(
    outdir,
    "STEP17G1_GSE221561_feature_hg38_mapping.csv"
  ),
  bom = TRUE
)

fwrite(
  mapping_summary,
  file.path(
    outdir,
    "STEP17G1_GSE221561_feature_mapping_summary.csv"
  ),
  bom = TRUE
)

fwrite(
  target_dup,
  file.path(
    outdir,
    "STEP17G1_GSE221561_duplicate_target_symbols.csv"
  ),
  bom = TRUE
)

# ============================================================
# 13. INPUT-FEATURE CATEGORY DETAIL
# ============================================================

cat("\n============================================================\n")
cat("12. FEATURE FORMAT SUMMARY\n")
cat("============================================================\n\n")

format_summary <- data.table(

  metric = c(
    "Total input features",
    "ENSG-like input features",
    "Exact hg38 symbol matches",
    "Ensembl-to-symbol matches",
    "Unmapped input features",
    "Mapped input rows",
    "Unique mapped hg38 symbols",
    "Duplicated hg38 target symbols"
  ),

  value = c(
    length(input_gene),

    sum(
      is_ensg
    ),

    sum(
      gene_map$mapping_type ==
        "EXACT_SYMBOL"
    ),

    sum(
      gene_map$mapping_type ==
        "ENSEMBL_TO_SYMBOL"
    ),

    sum(
      gene_map$mapping_type ==
        "UNMAPPED"
    ),

    nrow(mapped),

    unique_target_symbols,

    nrow(target_dup)
  )
)

print(format_summary)

fwrite(
  format_summary,
  file.path(
    outdir,
    "STEP17G1_GSE221561_feature_format_summary.csv"
  ),
  bom = TRUE
)

# ============================================================
# 14. PILOT / REFERENCE LAYER PLAN
# ============================================================

cat("\n============================================================\n")
cat("13. PILOT / REFERENCE LAYER PLAN\n")
cat("============================================================\n\n")

pilot_layer_plan <- epi_layer[
  sample_id %in%
    c(
      pilot_sample,
      normal_samples
    ),
  .(
    epithelial_cells = .N
  ),
  by = .(
    sample_id,
    tissue_simple,
    treatment_group,
    count_layer
  )
]

setorder(
  pilot_layer_plan,
  tissue_simple,
  sample_id,
  count_layer
)

print(
  pilot_layer_plan,
  nrows = 100
)

fwrite(
  pilot_layer_plan,
  file.path(
    outdir,
    "STEP17G1_GSE221561_pilot_reference_layer_plan.csv"
  ),
  bom = TRUE
)

# ============================================================
# 15. PREFLIGHT DECISION
# ============================================================

cat("\n============================================================\n")
cat("14. PREFLIGHT DECISION\n")
cat("============================================================\n\n")

raw_counts_ok <-
  TRUE

gene_mapping_ok <-
  unique_target_symbols >= 10000

duplicate_review_needed <-
  nrow(target_dup) > 0

if (!raw_counts_ok) {

  status <-
    "BLOCKED_RAW_COUNTS"

} else if (!gene_mapping_ok) {

  status <-
    "BLOCKED_INSUFFICIENT_HG38_MAPPING"

} else if (duplicate_review_needed) {

  status <-
    "READY_AFTER_DUPLICATE_FEATURE_REVIEW"

} else {

  status <-
    "READY_FOR_GSE221561_PILOT_PREPARATION"
}

decision <- data.table(

  dataset =
    "GSE221561",

  object_file =
    object_file,

  object_cells =
    ncol(obj),

  metadata_cells =
    nrow(md),

  RNA_features =
    nrow(rna),

  count_layers =
    length(count_layers),

  raw_counts_samplewise_accessible =
    raw_counts_ok,

  mapped_input_rows =
    nrow(mapped),

  unique_hg38_symbols =
    unique_target_symbols,

  duplicated_target_symbols =
    nrow(target_dup),

  normal_epithelial_samples =
    nrow(normal_inventory),

  normal_epithelial_cells =
    sum(
      normal_inventory$epithelial_cells
    ),

  selected_pilot_sample =
    pilot_sample,

  selected_pilot_patient =
    pilot_patient,

  selected_pilot_treatment =
    pilot_treatment,

  selected_pilot_epithelial_cells =
    pilot_epi_cells,

  selection_rule =
    selection_rule,

  preflight_status =
    status
)

print(decision)

fwrite(
  decision,
  file.path(
    outdir,
    "STEP17G1_GSE221561_CNV_preflight_decision.csv"
  ),
  bom = TRUE
)

# ============================================================
# 16. RECORD GSE197677 FREEZE
# ============================================================

freeze <- data.table(

  dataset =
    "GSE197677",

  status =
    "CNV_INCONCLUSIVE_AFTER_TWO_PILOTS",

  pilot_1 =
    "ESC06",

  pilot_2 =
    "ESC15",

  malignant_labels_assigned =
    FALSE,

  remaining_samples_processed =
    FALSE,

  infercna_parameters_changed =
    FALSE,

  next_action =
    "Independent technical validation in GSE221561 before expansion"
)

fwrite(
  freeze,
  file.path(
    outdir,
    "STEP17G1_GSE197677_CNV_FREEZE_STATUS.csv"
  ),
  bom = TRUE
)

# ============================================================
# CLEANUP
# ============================================================

rm(
  obj,
  rna
)

gc(verbose = FALSE)

# ============================================================
# FINAL
# ============================================================

cat("\n============================================================\n")
cat("STEP 17G1 FINISHED\n")
cat("GSE221561 CNV PREFLIGHT COMPLETE\n")
cat("============================================================\n\n")

cat(
  "Pilot sample              : ",
  pilot_sample,
  "\n",
  sep = ""
)

cat(
  "Pilot treatment           : ",
  pilot_treatment,
  "\n",
  sep = ""
)

cat(
  "Pilot epithelial cells    : ",
  pilot_epi_cells,
  "\n",
  sep = ""
)

cat(
  "Normal epithelial samples : ",
  nrow(normal_inventory),
  "\n",
  sep = ""
)

cat(
  "Normal epithelial cells   : ",
  sum(
    normal_inventory$epithelial_cells
  ),
  "\n",
  sep = ""
)

cat(
  "RNA input features        : ",
  length(input_gene),
  "\n",
  sep = ""
)

cat(
  "Unique mapped hg38 symbols: ",
  unique_target_symbols,
  "\n",
  sep = ""
)

cat(
  "Duplicate target symbols  : ",
  nrow(target_dup),
  "\n",
  sep = ""
)

cat(
  "Counts layer number       : ",
  length(count_layers),
  "\n",
  sep = ""
)

cat(
  "Preflight status          : ",
  status,
  "\n",
  sep = ""
)

cat("\nNO InferCNA run\n")
cat("NO JoinLayers\n")
cat("NO gene aggregation\n")
cat("NO gene IDs changed\n")
cat("NO malignant labels\n")
cat("NO cells filtered\n")
cat("GSE197677 remaining 8 samples NOT processed\n")

cat("\nSTOP HERE.\n")
