# ============================================================
# ESCC Neoadjuvant Project
# STEP 17F1
#
# Prepare SECOND independent GSE197677 InferCNA pilot.
#
# Strategy:
#   - exclude ESC06
#   - identify remaining tumor sample with MOST epithelial cells
#   - use EXACTLY the same four normal epithelial references
#   - use EXACTLY the same hg38 / CPM preparation as ESC06
#
# IMPORTANT:
#   NO InferCNA run
#   NO parameter changes
#   NO malignant calling
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
cat("STEP 17F1: PREPARE SECOND GSE197677 CNV PILOT\n")
cat("============================================================\n\n")

# ============================================================
# 1. PATHS
# ============================================================

object_dir <- "03_objects/GSE197677"

meta_file <- paste0(
  "06_tables/GSE197677/metadata/",
  "GSE197677_standardized_cell_metadata.csv.gz"
)

mapping_file <- paste0(
  "06_tables/GSE197677/celltype/",
  "GSE197677_FINAL_cluster_to_broad_celltype_mapping.csv"
)

base_outdir <- paste0(
  "06_tables/cross_dataset/",
  "malignant_epithelial"
)

dir.create(
  base_outdir,
  recursive = TRUE,
  showWarnings = FALSE
)

if (!file.exists(meta_file)) {
  stop("Missing: ", meta_file)
}

if (!file.exists(mapping_file)) {
  stop("Missing: ", mapping_file)
}

# ============================================================
# 2. LOAD METADATA + FINAL CELL-TYPE MAPPING
# ============================================================

cat("===== 1. LOAD METADATA =====\n\n")

md <- fread(meta_file)
mapping <- fread(mapping_file)

required_md <- c(
  "cell_id",
  "sample",
  "tissue_std",
  "treatment_group",
  "seurat_clusters"
)

missing_md <- setdiff(
  required_md,
  names(md)
)

if (length(missing_md) > 0) {
  stop(
    "Missing metadata fields: ",
    paste(missing_md, collapse = ", ")
  )
}

mapping[
  ,
  cluster_key :=
    as.character(cluster)
]

md[
  ,
  cluster_key :=
    as.character(seurat_clusters)
]

md <- merge(
  md,
  mapping[
    ,
    .(
      cluster_key,
      broad_celltype_final
    )
  ],
  by = "cluster_key",
  all.x = TRUE,
  sort = FALSE
)

if (
  any(
    is.na(
      md$broad_celltype_final
    )
  )
) {
  stop(
    "Some cells lack final broad-cell annotation."
  )
}

# ============================================================
# 3. AUDIT ALL TUMOR EPITHELIAL SAMPLES
# ============================================================

cat("\n============================================================\n")
cat("2. ALL GSE197677 TUMOR EPITHELIAL SAMPLES\n")
cat("============================================================\n\n")

tumor_candidates <- md[
  tissue_std == "Tumor" &
  broad_celltype_final == "Epithelial",
  .(
    epithelial_cells = .N
  ),
  by = .(
    sample,
    treatment_group
  )
]

setorder(
  tumor_candidates,
  -epithelial_cells
)

tumor_candidates[
  ,
  previous_ESC06_pilot :=
    sample == "ESC06"
]

print(
  tumor_candidates,
  nrows = 100
)

fwrite(
  tumor_candidates,
  file.path(
    base_outdir,
    "STEP17F1_GSE197677_tumor_epithelial_candidates.csv"
  ),
  bom = TRUE
)

# ============================================================
# 4. SELECT SECOND PILOT
# ============================================================

second_candidates <- tumor_candidates[
  sample != "ESC06" &
  epithelial_cells >= 30
]

if (nrow(second_candidates) == 0) {
  stop(
    "No second tumor epithelial pilot available."
  )
}

pilot_sample <-
  second_candidates$sample[1]

pilot_n <-
  second_candidates$epithelial_cells[1]

pilot_treatment <-
  second_candidates$treatment_group[1]

cat("\n============================================================\n")
cat("3. SECOND PILOT SELECTED\n")
cat("============================================================\n\n")

cat(
  "Sample            : ",
  pilot_sample,
  "\n",
  sep = ""
)

cat(
  "Treatment group   : ",
  pilot_treatment,
  "\n",
  sep = ""
)

cat(
  "Epithelial cells  : ",
  pilot_n,
  "\n",
  sep = ""
)

cat(
  "Selection rule    : largest tumor epithelial compartment excluding ESC06\n"
)

# ============================================================
# 5. DEFINE NORMAL REFERENCES
# ============================================================

normal_samples <- c(
  "ESN11",
  "ESN12",
  "ESN13",
  "ESN15"
)

tumor_cells <- md[
  sample == pilot_sample &
  tissue_std == "Tumor" &
  broad_celltype_final == "Epithelial",
  cell_id
]

normal_meta <- md[
  sample %in% normal_samples &
  tissue_std == "Normal" &
  broad_celltype_final == "Epithelial"
]

normal_cells <- normal_meta$cell_id

normal_summary <- normal_meta[
  ,
  .(
    epithelial_cells = .N
  ),
  by = sample
]

setorder(
  normal_summary,
  sample
)

cat("\n===== 4. NORMAL REFERENCE GROUPS =====\n\n")

print(normal_summary)

cat(
  "\nTotal normal reference cells: ",
  length(normal_cells),
  "\n",
  sep = ""
)

if (length(normal_cells) != 311) {
  warning(
    "Normal epithelial total differs from ESC06 pilot's 311 cells."
  )
}

if (uniqueN(normal_meta$sample) != 4) {
  stop(
    "Expected four normal reference groups."
  )
}

# ============================================================
# 6. LOCATE ORIGINAL OBJECT
# ============================================================

cat("\n===== 5. LOCATE ORIGINAL GSE197677 OBJECT =====\n\n")

files <- list.files(
  object_dir,
  pattern = "\\.rds$",
  full.names = TRUE,
  ignore.case = TRUE
)

hits <- files[
  basename(files) ==
    "GSE197677_integrated.rds"
]

if (length(hits) == 0) {

  hits <- files[
    grepl(
      "integrated",
      basename(files),
      ignore.case = TRUE
    ) &
    !grepl(
      "broadcelltype",
      basename(files),
      ignore.case = TRUE
    )
  ]
}

if (length(hits) == 0) {
  stop(
    "Original GSE197677 integrated RDS not found."
  )
}

object_file <- hits[1]

cat(
  "Object: ",
  object_file,
  "\n",
  sep = ""
)

# ============================================================
# 7. CONFIGURE hg38
# ============================================================

cat("\n===== 6. CONFIGURE hg38 =====\n\n")

infercna::useGenome(
  "hg38"
)

genome <- infercna::retrieveGenome(
  name = "hg38"
)

hg38_symbols <- unique(
  as.character(
    genome$symbol
  )
)

cat(
  "hg38 reference genes: ",
  length(hg38_symbols),
  "\n",
  sep = ""
)

# ============================================================
# 8. LOAD OBJECT + RAW COUNTS
# ============================================================

cat("\n===== 7. LOAD OBJECT =====\n\n")

cat(
  "Loading original GSE197677 object...\n"
)

obj <- readRDS(
  object_file
)

if (!"RNA" %in% Assays(obj)) {
  stop(
    "RNA assay missing."
  )
}

counts <- tryCatch(
  Seurat::GetAssayData(
    object = obj,
    assay = "RNA",
    layer = "counts"
  ),
  error = function(e) {
    Seurat::GetAssayData(
      object = obj,
      assay = "RNA",
      slot = "counts"
    )
  }
)

cat(
  "Raw counts: ",
  nrow(counts),
  " genes x ",
  ncol(counts),
  " cells\n",
  sep = ""
)

# ============================================================
# 9. VALIDATE CELL IDS
# ============================================================

cat("\n===== 8. CELL-ID VALIDATION =====\n\n")

selected_cells <- c(
  tumor_cells,
  normal_cells
)

missing_cells <- setdiff(
  selected_cells,
  colnames(counts)
)

cat(
  "Selected cells missing from counts: ",
  length(missing_cells),
  "\n",
  sep = ""
)

if (length(missing_cells) > 0) {

  print(
    head(
      missing_cells,
      20
    )
  )

  stop(
    "Selected cell IDs are incomplete."
  )
}

# ============================================================
# 10. EXTRACT PILOT
# ============================================================

cat("\n===== 9. EXTRACT SECOND PILOT RAW COUNTS =====\n\n")

pilot_counts <- counts[
  ,
  selected_cells,
  drop = FALSE
]

cat(
  "Pilot matrix: ",
  nrow(pilot_counts),
  " genes x ",
  ncol(pilot_counts),
  " cells\n",
  sep = ""
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

# ============================================================
# 11. hg38 SYMBOL FILTER
# ============================================================

cat("\n===== 10. hg38 GENE OVERLAP =====\n\n")

input_genes <- rownames(
  pilot_counts
)

overlap <- intersect(
  input_genes,
  hg38_symbols
)

cat(
  "Input genes : ",
  length(input_genes),
  "\n",
  sep = ""
)

cat(
  "hg38 genes  : ",
  length(overlap),
  "\n",
  sep = ""
)

cat(
  sprintf(
    "Overlap     : %.2f%%\n",
    100 *
      length(overlap) /
      length(input_genes)
  )
)

if (length(overlap) < 5000) {
  stop(
    "Too few hg38 genes."
  )
}

overlap <- unique(
  overlap
)

pilot_counts <- pilot_counts[
  overlap,
  ,
  drop = FALSE
]

keep <- Matrix::rowSums(
  pilot_counts
) > 0

pilot_counts <- pilot_counts[
  keep,
  ,
  drop = FALSE
]

cat(
  "Non-zero hg38 genes: ",
  nrow(pilot_counts),
  "\n",
  sep = ""
)

# ============================================================
# 12. CPM
# ============================================================

cat("\n===== 11. CPM PREPARATION =====\n\n")

library_sizes <- Matrix::colSums(
  pilot_counts
)

if (any(library_sizes <= 0)) {
  stop(
    "Zero-library cell detected."
  )
}

pilot_cpm <- pilot_counts %*%
  Matrix::Diagonal(
    x = 1e6 / library_sizes
  )

rownames(pilot_cpm) <-
  rownames(pilot_counts)

colnames(pilot_cpm) <-
  colnames(pilot_counts)

cpm_sums <- Matrix::colSums(
  pilot_cpm
)

cat(
  sprintf(
    "Median CPM column sum: %.2f\n",
    median(cpm_sums)
  )
)

if (
  max(
    abs(
      cpm_sums - 1e6
    )
  ) > 1
) {
  stop(
    "CPM validation failed."
  )
}

cat("✓ CPM validated\n")

# ============================================================
# 13. NORMAL REFERENCE LIST
# ============================================================

cat("\n===== 12. BUILD REFERENCE LIST =====\n\n")

refCells <- lapply(
  normal_samples,
  function(s) {

    normal_meta[
      sample == s,
      cell_id
    ]
  }
)

names(refCells) <- normal_samples

print(
  lengths(refCells)
)

if (
  !all(
    unlist(
      refCells,
      use.names = FALSE
    ) %in%
      colnames(pilot_cpm)
  )
) {
  stop(
    "Reference IDs invalid."
  )
}

# ============================================================
# 14. ANNOTATION
# ============================================================

annotation <- md[
  cell_id %in%
    selected_cells,
  .(
    cell_id,
    sample,
    tissue =
      tissue_std,
    treatment_group,
    broad_celltype =
      broad_celltype_final
  )
]

annotation[
  ,
  infercna_role :=
    fifelse(
      cell_id %in%
        tumor_cells,
      "TUMOR_TEST",
      "NORMAL_REFERENCE"
    )
]

annotation[
  ,
  reference_group :=
    fifelse(
      infercna_role ==
        "NORMAL_REFERENCE",
      sample,
      NA_character_
    )
]

role_summary <- annotation[
  ,
  .(
    cells = .N
  ),
  by = .(
    infercna_role,
    sample,
    treatment_group
  )
]

cat("\n===== 13. PILOT CELL COMPOSITION =====\n\n")

print(role_summary)

# ============================================================
# 15. GENE TABLE
# ============================================================

gene_table <- data.table(
  gene =
    rownames(
      pilot_cpm
    )
)

genome_subset <- genome[
  match(
    gene_table$gene,
    genome$symbol
  ),
]

gene_table[
  ,
  chromosome :=
    as.character(
      genome_subset$chromosome_name
    )
]

gene_table[
  ,
  arm :=
    as.character(
      genome_subset$arm
    )
]

gene_table[
  ,
  start_position :=
    genome_subset$start_position
]

gene_table[
  ,
  end_position :=
    genome_subset$end_position
]

if (
  any(
    is.na(
      gene_table$chromosome
    )
  )
) {
  stop(
    "Genome mapping incomplete."
  )
}

# ============================================================
# 16. MEMORY
# ============================================================

cat("\n===== 14. MEMORY ESTIMATE =====\n\n")

dense_gb <- (
  nrow(pilot_cpm) *
  ncol(pilot_cpm) *
  8
) / 1024^3

cat(
  sprintf(
    "Dense matrix estimate       : %.3f GB\n",
    dense_gb
  )
)

cat(
  sprintf(
    "Conservative 4x working set : %.3f GB\n",
    dense_gb * 4
  )
)

# ============================================================
# 17. OUTPUT DIRECTORY
# ============================================================

pilot_dir <- file.path(
  base_outdir,
  paste0(
    "GSE197677_",
    pilot_sample,
    "_pilot"
  )
)

dir.create(
  pilot_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ============================================================
# 18. SAVE INPUTS
# ============================================================

cat("\n===== 15. SAVE SECOND PILOT INPUT =====\n\n")

raw_file <- file.path(
  pilot_dir,
  paste0(
    "GSE197677_",
    pilot_sample,
    "_infercna_raw_counts_hg38.rds"
  )
)

cpm_file <- file.path(
  pilot_dir,
  paste0(
    "GSE197677_",
    pilot_sample,
    "_infercna_CPM_hg38.rds"
  )
)

refs_file <- file.path(
  pilot_dir,
  paste0(
    "GSE197677_",
    pilot_sample,
    "_infercna_refCells.rds"
  )
)

annotation_file <- file.path(
  pilot_dir,
  paste0(
    "GSE197677_",
    pilot_sample,
    "_infercna_cell_annotation.csv"
  )
)

genes_file <- file.path(
  pilot_dir,
  paste0(
    "GSE197677_",
    pilot_sample,
    "_infercna_gene_table.csv"
  )
)

manifest_file <- file.path(
  pilot_dir,
  paste0(
    "GSE197677_",
    pilot_sample,
    "_infercna_manifest.csv"
  )
)

saveRDS(
  pilot_counts,
  raw_file,
  compress = FALSE
)

saveRDS(
  pilot_cpm,
  cpm_file,
  compress = FALSE
)

saveRDS(
  refCells,
  refs_file
)

fwrite(
  annotation,
  annotation_file,
  bom = TRUE
)

fwrite(
  gene_table,
  genes_file,
  bom = TRUE
)

manifest <- data.table(

  dataset =
    "GSE197677",

  pilot_number =
    2L,

  pilot_sample =
    pilot_sample,

  selection_rule =
    "largest tumor epithelial compartment excluding ESC06",

  treatment_group =
    pilot_treatment,

  tumor_epithelial_cells =
    length(tumor_cells),

  normal_reference_groups =
    length(refCells),

  normal_reference_cells =
    length(normal_cells),

  total_cells =
    ncol(pilot_cpm),

  hg38_genes =
    nrow(pilot_cpm),

  normalization =
    "CPM",

  isLog =
    FALSE,

  planned_window =
    100,

  planned_n =
    5000,

  planned_noise =
    0.1,

  planned_center_method =
    "median",

  parameter_change_from_ESC06 =
    FALSE,

  infercna_run =
    FALSE
)

fwrite(
  manifest,
  manifest_file,
  bom = TRUE
)

# ============================================================
# 19. RELOAD VALIDATION
# ============================================================

cat("\n===== 16. RELOAD VALIDATION =====\n\n")

check <- readRDS(
  cpm_file
)

check_refs <- readRDS(
  refs_file
)

cat(
  "Reloaded CPM: ",
  nrow(check),
  " genes x ",
  ncol(check),
  " cells\n",
  sep = ""
)

cat(
  "Reference groups: ",
  length(check_refs),
  "\n",
  sep = ""
)

if (
  !all(
    unlist(
      check_refs,
      use.names = FALSE
    ) %in%
      colnames(check)
  )
) {
  stop(
    "Reference validation failed."
  )
}

if (
  !all(
    tumor_cells %in%
      colnames(check)
  )
) {
  stop(
    "Tumor cell validation failed."
  )
}

# ============================================================
# CLEANUP
# ============================================================

rm(
  obj,
  counts,
  pilot_counts,
  pilot_cpm,
  check
)

gc(verbose = FALSE)

# ============================================================
# FINAL
# ============================================================

cat("\n============================================================\n")
cat("STEP 17F1 FINISHED ✓\n")
cat("SECOND GSE197677 CNV PILOT INPUT READY ✓\n")
cat("============================================================\n\n")

cat(
  "Second pilot sample      : ",
  pilot_sample,
  "\n",
  sep = ""
)

cat(
  "Treatment group          : ",
  pilot_treatment,
  "\n",
  sep = ""
)

cat(
  "Tumor epithelial cells   : ",
  length(tumor_cells),
  "\n",
  sep = ""
)

cat(
  "Normal reference cells   : ",
  length(normal_cells),
  "\n",
  sep = ""
)

cat(
  "Normal reference groups  : ",
  length(refCells),
  "\n",
  sep = ""
)

cat(
  "Total pilot cells        : ",
  length(selected_cells),
  "\n",
  sep = ""
)

cat(
  "Final hg38 genes         : ",
  nrow(gene_table),
  "\n",
  sep = ""
)

cat("\nSAME PARAMETER PLAN AS ESC06:\n")
cat("  window = 100\n")
cat("  n = 5000\n")
cat("  noise = 0.1\n")
cat("  center.method = median\n")
cat("  isLog = FALSE\n")

cat("\nNO InferCNA run ✓\n")
cat("NO parameters changed ✓\n")
cat("NO malignant labels assigned ✓\n")
cat("NO cells filtered ✓\n")

cat("\nSTOP HERE.\n")
cat("DO NOT RUN THE REMAINING SAMPLES.\n")
