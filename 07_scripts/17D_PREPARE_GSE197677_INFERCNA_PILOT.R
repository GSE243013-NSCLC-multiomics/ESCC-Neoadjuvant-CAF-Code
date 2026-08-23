# ============================================================
# ESCC Neoadjuvant Project
# STEP 17D
#
# Prepare GSE197677 InferCNA pilot input
#
# Pilot tumor:
#   ESC06 epithelial cells
#
# Normal epithelial references:
#   ESN11
#   ESN12
#   ESN13
#   ESN15
#
# OUTPUT:
#   raw sparse counts
#   CPM sparse matrix
#   reference cell list
#   cell metadata
#   gene audit
#
# IMPORTANT:
#   NO InferCNA run
#   NO CNA inference
#   NO malignant calling
#   NO filtering of cells
#   NO JoinLayers
# ============================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(infercna)
})

cat("\n============================================================\n")
cat("STEP 17D: PREPARE GSE197677 INFERCNA PILOT INPUT\n")
cat("============================================================\n\n")

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------

object_dir <- "03_objects/GSE197677"

outdir <- paste0(
  "06_tables/cross_dataset/",
  "malignant_epithelial/",
  "GSE197677_ESC06_pilot"
)

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# 1. Genome
# ------------------------------------------------------------

cat("===== 1. CONFIGURE hg38 =====\n\n")

infercna::useGenome("hg38")

genome <- infercna::retrieveGenome(
  name = "hg38"
)

cat(
  "hg38 genes: ",
  nrow(genome),
  "\n",
  sep = ""
)

hg38_symbols <- unique(
  as.character(
    genome$symbol
  )
)

# ------------------------------------------------------------
# 2. Locate original GSE197677 object
# ------------------------------------------------------------

cat("\n===== 2. LOCATE ORIGINAL OBJECT =====\n\n")

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
    "Original GSE197677 integrated object not found."
  )
}

object_file <- hits[1]

cat(
  "Object: ",
  object_file,
  "\n",
  sep = ""
)

# ------------------------------------------------------------
# 3. Load standardized metadata + final mapping
# ------------------------------------------------------------

cat("\n===== 3. LOAD METADATA =====\n\n")

meta_file <- paste0(
  "06_tables/GSE197677/metadata/",
  "GSE197677_standardized_cell_metadata.csv.gz"
)

mapping_file <- paste0(
  "06_tables/GSE197677/celltype/",
  "GSE197677_FINAL_cluster_to_broad_celltype_mapping.csv"
)

if (!file.exists(meta_file)) {
  stop("Missing: ", meta_file)
}

if (!file.exists(mapping_file)) {
  stop("Missing: ", mapping_file)
}

md <- fread(meta_file)
mapping <- fread(mapping_file)

needed <- c(
  "cell_id",
  "sample",
  "tissue_std",
  "treatment_group",
  "seurat_clusters"
)

miss <- setdiff(
  needed,
  names(md)
)

if (length(miss) > 0) {
  stop(
    "Metadata missing: ",
    paste(miss, collapse = ", ")
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
    "Some cells lack broad-cell annotation."
  )
}

# ------------------------------------------------------------
# 4. Define pilot cells
# ------------------------------------------------------------

cat("\n===== 4. DEFINE PILOT CELLS =====\n\n")

pilot_tumor <- "ESC06"

normal_samples <- c(
  "ESN11",
  "ESN12",
  "ESN13",
  "ESN15"
)

tumor_cells <- md[
  sample == pilot_tumor &
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

cat(
  "Pilot tumor: ",
  pilot_tumor,
  "\n",
  sep = ""
)

cat(
  "Tumor epithelial cells: ",
  length(tumor_cells),
  "\n",
  sep = ""
)

cat(
  "Normal epithelial cells: ",
  length(normal_cells),
  "\n\n",
  sep = ""
)

normal_counts <- normal_meta[
  ,
  .N,
  by = sample
]

setorder(
  normal_counts,
  sample
)

print(normal_counts)

if (length(tumor_cells) < 30) {
  stop(
    "Pilot tumor has <30 epithelial cells."
  )
}

if (uniqueN(normal_meta$sample) < 2) {
  stop(
    "Need at least two independent normal reference groups."
  )
}

all_selected_cells <- c(
  tumor_cells,
  normal_cells
)

if (
  anyDuplicated(
    all_selected_cells
  ) > 0
) {
  stop(
    "Duplicate selected cell IDs detected."
  )
}

cat(
  "\nTotal pilot cells: ",
  length(all_selected_cells),
  "\n",
  sep = ""
)

# ------------------------------------------------------------
# 5. Load object
# ------------------------------------------------------------

cat("\n===== 5. LOAD GSE197677 OBJECT =====\n\n")

cat(
  "Loading object; this may take several minutes...\n"
)

obj <- readRDS(
  object_file
)

if (!"RNA" %in% Assays(obj)) {
  stop(
    "RNA assay missing."
  )
}

rna_counts <- tryCatch(
  GetAssayData(
    obj,
    assay = "RNA",
    layer = "counts"
  ),
  error = function(e) {
    GetAssayData(
      obj,
      assay = "RNA",
      slot = "counts"
    )
  }
)

cat(
  "Full RNA counts: ",
  nrow(rna_counts),
  " genes x ",
  ncol(rna_counts),
  " cells\n",
  sep = ""
)

# ------------------------------------------------------------
# 6. Check selected cell IDs
# ------------------------------------------------------------

cat("\n===== 6. CELL-ID VALIDATION =====\n\n")

missing_cells <- setdiff(
  all_selected_cells,
  colnames(rna_counts)
)

cat(
  "Selected cells absent from counts: ",
  length(missing_cells),
  "\n",
  sep = ""
)

if (length(missing_cells) > 0) {

  cat(
    "Examples:\n"
  )

  print(
    head(
      missing_cells,
      20
    )
  )

  stop(
    "Selected cell IDs do not all exist in RNA counts."
  )
}

# ------------------------------------------------------------
# 7. Extract raw counts
# ------------------------------------------------------------

cat("\n===== 7. EXTRACT RAW COUNTS =====\n\n")

pilot_counts <- rna_counts[
  ,
  all_selected_cells,
  drop = FALSE
]

cat(
  "Pilot raw matrix: ",
  nrow(pilot_counts),
  " genes x ",
  ncol(pilot_counts),
  " cells\n",
  sep = ""
)

cat(
  "Matrix class: ",
  paste(
    class(pilot_counts),
    collapse = ", "
  ),
  "\n",
  sep = ""
)

cat(
  "Non-zero entries: ",
  Matrix::nnzero(
    pilot_counts
  ),
  "\n",
  sep = ""
)

# ------------------------------------------------------------
# 8. Gene symbol / hg38 filtering
# ------------------------------------------------------------

cat("\n===== 8. hg38 GENE OVERLAP =====\n\n")

genes_raw <- rownames(
  pilot_counts
)

overlap <- genes_raw[
  genes_raw %in%
    hg38_symbols
]

cat(
  "Input genes     : ",
  length(genes_raw),
  "\n",
  sep = ""
)

cat(
  "hg38 overlap    : ",
  length(overlap),
  "\n",
  sep = ""
)

cat(
  sprintf(
    "Overlap rate   : %.2f%%\n",
    100 *
      length(overlap) /
      length(genes_raw)
  )
)

if (length(overlap) < 5000) {
  stop(
    "Too few genes overlap InferCNA hg38 reference."
  )
}

overlap <- unique(
  overlap
)

pilot_counts_hg38 <- pilot_counts[
  overlap,
  ,
  drop = FALSE
]

gene_sums <- Matrix::rowSums(
  pilot_counts_hg38
)

pilot_counts_hg38 <- pilot_counts_hg38[
  gene_sums > 0,
  ,
  drop = FALSE
]

cat(
  "hg38 genes with non-zero pilot expression: ",
  nrow(pilot_counts_hg38),
  "\n",
  sep = ""
)

# ------------------------------------------------------------
# 9. Library sizes
# ------------------------------------------------------------

cat("\n===== 9. LIBRARY SIZE AUDIT =====\n\n")

library_sizes <- Matrix::colSums(
  pilot_counts_hg38
)

lib_summary <- data.table(
  cell_id =
    names(library_sizes),

  library_size_hg38 =
    as.numeric(library_sizes)
)

lib_summary <- merge(
  lib_summary,
  md[
    ,
    .(
      cell_id,
      sample,
      tissue_std,
      treatment_group,
      broad_celltype_final
    )
  ],
  by = "cell_id",
  all.x = TRUE
)

cat(
  "Minimum library size: ",
  min(
    lib_summary$library_size_hg38
  ),
  "\n",
  sep = ""
)

cat(
  "Median library size : ",
  median(
    lib_summary$library_size_hg38
  ),
  "\n",
  sep = ""
)

cat(
  "Maximum library size: ",
  max(
    lib_summary$library_size_hg38
  ),
  "\n",
  sep = ""
)

if (
  any(
    lib_summary$library_size_hg38 <= 0
  )
) {
  stop(
    "At least one selected cell has zero hg38-mapped library size."
  )
}

# ------------------------------------------------------------
# 10. CPM normalization
# ------------------------------------------------------------

cat("\n===== 10. PREPARE CPM MATRIX =====\n\n")

scale_factor <- 1e6 /
  library_sizes

pilot_cpm <- pilot_counts_hg38 %*%
  Matrix::Diagonal(
    x = scale_factor
  )

rownames(pilot_cpm) <-
  rownames(
    pilot_counts_hg38
  )

colnames(pilot_cpm) <-
  colnames(
    pilot_counts_hg38
  )

cat(
  "CPM matrix: ",
  nrow(pilot_cpm),
  " genes x ",
  ncol(pilot_cpm),
  " cells\n",
  sep = ""
)

cat(
  "CPM class : ",
  paste(
    class(pilot_cpm),
    collapse = ", "
  ),
  "\n",
  sep = ""
)

cpm_totals <- Matrix::colSums(
  pilot_cpm
)

cat(
  sprintf(
    "Median CPM column sum: %.2f\n",
    median(cpm_totals)
  )
)

if (
  max(
    abs(
      cpm_totals -
        1e6
    )
  ) > 1
) {

  stop(
    "CPM column sums failed validation."
  )
}

cat("✓ CPM column sums validated\n")

# ------------------------------------------------------------
# 11. Build InferCNA refCells list
# ------------------------------------------------------------

cat("\n===== 11. BUILD NORMAL REFERENCE LIST =====\n\n")

refCells <- lapply(
  normal_samples,
  function(s) {

    normal_meta[
      sample == s,
      cell_id
    ]
  }
)

names(refCells) <-
  normal_samples

ref_sizes <- lengths(
  refCells
)

print(ref_sizes)

if (length(refCells) < 2) {
  stop(
    "InferCNA reference list has <2 groups."
  )
}

if (
  any(
    ref_sizes == 0
  )
) {
  stop(
    "At least one reference group is empty."
  )
}

all_ref_cells <- unlist(
  refCells,
  use.names = FALSE
)

missing_ref <- setdiff(
  all_ref_cells,
  colnames(pilot_cpm)
)

if (length(missing_ref) > 0) {
  stop(
    "Some reference cells are absent from CPM matrix."
  )
}

tumor_ref_overlap <- intersect(
  tumor_cells,
  all_ref_cells
)

if (length(tumor_ref_overlap) > 0) {
  stop(
    "Tumor cells accidentally present in reference list."
  )
}

cat(
  "✓ ",
  length(refCells),
  " independent normal reference groups validated\n",
  sep = ""
)

# ------------------------------------------------------------
# 12. Pilot cell annotation
# ------------------------------------------------------------

cat("\n===== 12. PILOT CELL ANNOTATION =====\n\n")

pilot_annotation <- md[
  cell_id %in%
    all_selected_cells,
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

pilot_annotation[
  ,
  infercna_role :=
    fifelse(
      cell_id %in%
        tumor_cells,
      "TUMOR_TEST",
      "NORMAL_REFERENCE"
    )
]

pilot_annotation[
  ,
  reference_group :=
    fifelse(
      infercna_role ==
        "NORMAL_REFERENCE",
      sample,
      NA_character_
    )
]

setorder(
  pilot_annotation,
  infercna_role,
  sample,
  cell_id
)

role_summary <- pilot_annotation[
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

print(role_summary)

# ------------------------------------------------------------
# 13. Gene table
# ------------------------------------------------------------

cat("\n===== 13. PILOT GENE TABLE =====\n\n")

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
    "Genome annotation missing after hg38 filtering."
  )
}

cat(
  "Final pilot genes: ",
  nrow(gene_table),
  "\n",
  sep = ""
)

# ------------------------------------------------------------
# 14. Dense-memory estimate
# ------------------------------------------------------------

cat("\n===== 14. MEMORY ESTIMATE =====\n\n")

dense_gb <- (
  nrow(pilot_cpm) *
  ncol(pilot_cpm) *
  8
) / 1024^3

cat(
  sprintf(
    "One dense genes x cells matrix: %.3f GB\n",
    dense_gb
  )
)

cat(
  sprintf(
    "Conservative 4x working estimate: %.3f GB\n",
    dense_gb * 4
  )
)

# ------------------------------------------------------------
# 15. Save pilot inputs
# ------------------------------------------------------------

cat("\n===== 15. SAVE PILOT INPUTS =====\n\n")

counts_file <- file.path(
  outdir,
  "GSE197677_ESC06_infercna_raw_counts_hg38.rds"
)

cpm_file <- file.path(
  outdir,
  "GSE197677_ESC06_infercna_CPM_hg38.rds"
)

refs_file <- file.path(
  outdir,
  "GSE197677_ESC06_infercna_refCells.rds"
)

annotation_file <- file.path(
  outdir,
  "GSE197677_ESC06_infercna_cell_annotation.csv"
)

genes_file <- file.path(
  outdir,
  "GSE197677_ESC06_infercna_gene_table.csv"
)

manifest_file <- file.path(
  outdir,
  "GSE197677_ESC06_infercna_manifest.csv"
)

saveRDS(
  pilot_counts_hg38,
  counts_file,
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
  pilot_annotation,
  annotation_file,
  bom = TRUE
)

fwrite(
  gene_table,
  genes_file,
  bom = TRUE
)

manifest <- data.table(
  dataset = "GSE197677",
  pilot_sample = pilot_tumor,
  tumor_epithelial_cells =
    length(tumor_cells),
  normal_reference_groups =
    length(refCells),
  normal_reference_cells =
    length(all_ref_cells),
  total_cells =
    ncol(pilot_cpm),
  hg38_genes =
    nrow(pilot_cpm),
  normalization =
    "CPM",
  log_transformed =
    FALSE,
  infercna_run =
    FALSE
)

fwrite(
  manifest,
  manifest_file,
  bom = TRUE
)

cat(
  "✓ ",
  counts_file,
  "\n",
  sep = ""
)

cat(
  "✓ ",
  cpm_file,
  "\n",
  sep = ""
)

cat(
  "✓ ",
  refs_file,
  "\n",
  sep = ""
)

cat(
  "✓ ",
  annotation_file,
  "\n",
  sep = ""
)

cat(
  "✓ ",
  genes_file,
  "\n",
  sep = ""
)

# ------------------------------------------------------------
# 16. Reload validation
# ------------------------------------------------------------

cat("\n===== 16. RELOAD VALIDATION =====\n\n")

check_cpm <- readRDS(
  cpm_file
)

check_refs <- readRDS(
  refs_file
)

cat(
  "Reloaded CPM: ",
  nrow(check_cpm),
  " x ",
  ncol(check_cpm),
  "\n",
  sep = ""
)

cat(
  "Reloaded reference groups: ",
  length(check_refs),
  "\n",
  sep = ""
)

all_refs_valid <- all(
  unlist(
    check_refs,
    use.names = FALSE
  ) %in%
    colnames(
      check_cpm
    )
)

cat(
  "All reference IDs present: ",
  all_refs_valid,
  "\n",
  sep = ""
)

if (!all_refs_valid) {
  stop(
    "Reloaded reference validation failed."
  )
}

if (
  !all(
    tumor_cells %in%
      colnames(check_cpm)
  )
) {
  stop(
    "Pilot tumor cells missing after reload."
  )
}

# ------------------------------------------------------------
# Cleanup
# ------------------------------------------------------------

rm(
  obj,
  rna_counts,
  pilot_counts,
  pilot_counts_hg38,
  pilot_cpm,
  check_cpm
)

gc(verbose = FALSE)

# ------------------------------------------------------------
# FINAL
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP 17D FINISHED ✓\n")
cat("GSE197677 ESC06 INFERCNA INPUT READY ✓\n")
cat("============================================================\n\n")

cat(
  "Pilot tumor sample       : ",
  pilot_tumor,
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
  "Normal reference groups  : ",
  length(refCells),
  "\n",
  sep = ""
)

cat(
  "Normal epithelial cells  : ",
  length(all_ref_cells),
  "\n",
  sep = ""
)

cat(
  "Total pilot cells        : ",
  length(all_selected_cells),
  "\n",
  sep = ""
)

cat(
  "Final hg38 genes         : ",
  nrow(gene_table),
  "\n",
  sep = ""
)

cat("\nCPM validated ✓\n")
cat("Reference IDs validated ✓\n")
cat("Tumor/reference separation validated ✓\n")

cat("\nNO InferCNA run ✓\n")
cat("NO CNA inference ✓\n")
cat("NO malignant labels assigned ✓\n")
cat("NO cells filtered ✓\n")
cat("NO JoinLayers ✓\n")

cat("\nSTOP HERE.\n")
cat("DO NOT RUN CNV YET.\n")
