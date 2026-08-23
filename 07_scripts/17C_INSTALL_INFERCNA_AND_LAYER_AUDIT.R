# ============================================================
# ESCC Neoadjuvant Project
# STEP 17C
#
# 1. Install / verify infercna in current renv project
# 2. Verify hg38 genome reference
# 3. Audit GSE197677 RNA count layers
# 4. Map sample / epithelial cells to count layers
# 5. Prove raw counts can be extracted sample-wise
#
# NO CNV inference
# NO JoinLayers
# NO normalization
# NO malignant calling
# NO cell filtering
# ============================================================

options(
  stringsAsFactors = FALSE,
  timeout = 1200
)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(data.table)
})

cat("\n============================================================\n")
cat("STEP 17C: INFERCNA INSTALL + GSE197677 LAYER AUDIT\n")
cat("============================================================\n\n")

outdir <- "06_tables/cross_dataset/malignant_epithelial"

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ============================================================
# 1. INSTALL / VERIFY INFERCNA
# ============================================================

cat("===== 1. INFERCNA INSTALLATION =====\n\n")

if (!requireNamespace("renv", quietly = TRUE)) {
  stop(
    "renv is not available in this project."
  )
}

if (!requireNamespace("infercna", quietly = TRUE)) {

  cat(
    "infercna not installed.\n",
    "Installing into the current renv project...\n\n",
    sep = ""
  )

  tryCatch(
    {
      renv::install(
        "jlaffy/infercna"
      )
    },
    error = function(e) {

      cat(
        "\nINSTALLATION FAILED:\n",
        conditionMessage(e),
        "\n"
      )

      stop(
        "infercna installation failed. ",
        "STOP before any CNV analysis."
      )
    }
  )
}

if (!requireNamespace("infercna", quietly = TRUE)) {
  stop(
    "infercna still unavailable after installation."
  )
}

cat(
  "✓ infercna installed\n"
)

cat(
  "Version: ",
  as.character(
    packageVersion("infercna")
  ),
  "\n",
  sep = ""
)

# ============================================================
# 2. VERIFY hg38
# ============================================================

cat("\n===== 2. INFERCNA hg38 GENOME =====\n\n")

infercna::useGenome(
  "hg38"
)

genome <- infercna::retrieveGenome(
  name = "hg38"
)

cat(
  "Genome rows: ",
  nrow(genome),
  "\n",
  sep = ""
)

cat(
  "Genome columns: ",
  paste(
    names(genome),
    collapse = ", "
  ),
  "\n",
  sep = ""
)

required_genome_cols <- c(
  "symbol",
  "chromosome_name",
  "start_position",
  "arm"
)

missing_genome_cols <- setdiff(
  required_genome_cols,
  names(genome)
)

if (length(missing_genome_cols) > 0) {
  stop(
    "hg38 genome missing expected columns: ",
    paste(
      missing_genome_cols,
      collapse = ", "
    )
  )
}

cat("✓ hg38 reference available\n")

# ============================================================
# 3. LOCATE ORIGINAL GSE197677 OBJECT
# ============================================================

cat("\n===== 3. LOCATE GSE197677 OBJECT =====\n\n")

object_dir <- "03_objects/GSE197677"

all_rds <- list.files(
  object_dir,
  pattern = "\\.rds$",
  full.names = TRUE,
  ignore.case = TRUE
)

preferred <- all_rds[
  basename(all_rds) ==
    "GSE197677_integrated.rds"
]

if (length(preferred) == 0) {

  preferred <- all_rds[
    grepl(
      "integrated",
      basename(all_rds),
      ignore.case = TRUE
    ) &
    !grepl(
      "broadcelltype",
      basename(all_rds),
      ignore.case = TRUE
    )
  ]
}

if (length(preferred) == 0) {
  stop(
    "Original GSE197677 integrated RDS not found."
  )
}

object_file <- preferred[1]

cat(
  "Object: ",
  object_file,
  "\n",
  sep = ""
)

cat(
  sprintf(
    "Size: %.2f GB\n",
    file.info(object_file)$size /
      1024^3
  )
)

# ============================================================
# 4. LOAD LIGHTWEIGHT METADATA / FINAL MAPPING
# ============================================================

cat("\n===== 4. LOAD CELL METADATA + FINAL CLUSTER MAPPING =====\n\n")

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

cell_meta <- fread(
  meta_file
)

mapping <- fread(
  mapping_file
)

required_meta <- c(
  "cell_id",
  "sample",
  "subject_id_internal",
  "tissue_std",
  "treatment_group",
  "seurat_clusters"
)

miss <- setdiff(
  required_meta,
  names(cell_meta)
)

if (length(miss) > 0) {
  stop(
    "Cell metadata missing fields: ",
    paste(miss, collapse = ", ")
  )
}

mapping[
  ,
  cluster_key :=
    as.character(cluster)
]

cell_meta[
  ,
  cluster_key :=
    as.character(seurat_clusters)
]

cell_meta <- merge(
  cell_meta,
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
      cell_meta$broad_celltype_final
    )
  )
) {
  stop(
    "Some GSE197677 cells lack final broad-cell labels."
  )
}

cat(
  "Cell metadata rows: ",
  nrow(cell_meta),
  "\n",
  sep = ""
)

cat(
  "Epithelial cells : ",
  sum(
    cell_meta$broad_celltype_final ==
      "Epithelial"
  ),
  "\n",
  sep = ""
)

# ============================================================
# 5. LOAD SEURAT OBJECT
# ============================================================

cat("\n===== 5. LOAD ORIGINAL OBJECT =====\n\n")

cat(
  "Loading 5.8 GB object; this may take several minutes...\n"
)

obj <- readRDS(
  object_file
)

if (!"RNA" %in% Assays(obj)) {
  stop(
    "RNA assay not found."
  )
}

rna <- obj@assays[["RNA"]]

cat(
  "RNA assay class: ",
  paste(
    class(rna),
    collapse = ", "
  ),
  "\n",
  sep = ""
)

# ============================================================
# 6. LIST RNA LAYERS
# ============================================================

cat("\n===== 6. RNA LAYERS =====\n\n")

layers <- tryCatch(
  SeuratObject::Layers(
    rna,
    search = NA
  ),
  error = function(e) {
    cat(
      "SeuratObject::Layers failed, trying alternative...\n"
    )
    # Alternative: check slot names
    tryCatch({
      sn <- slotNames(rna)
      sn[grepl("counts", sn, ignore.case = TRUE)]
    }, error = function(e2) character(0))
  }
)

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
    "No RNA count layers found."
  )
}

# ============================================================
# 7. MAP CELL -> COUNT LAYER WITHOUT JOINING
# ============================================================

cat("\n===== 7. CELL MEMBERSHIP BY COUNT LAYER =====\n\n")

layer_membership_list <- list()

for (lyr in count_layers) {

  lyr_cells <- tryCatch(
    SeuratObject::Cells(
      rna,
      layer = lyr
    ),
    error = function(e) {
      cat(
        "  Cells() failed for",
        lyr,
        ":",
        conditionMessage(e),
        "\n"
      )
      character(0)
    }
  )

  cat(
    sprintf(
      "%-40s %7d cells\n",
      lyr,
      length(lyr_cells)
    )
  )

  if (length(lyr_cells) > 0) {
    layer_membership_list[[lyr]] <- data.table(
      cell_id = lyr_cells,
      count_layer = lyr
    )
  }
}

layer_membership <- rbindlist(
  layer_membership_list,
  use.names = TRUE
)

layer_dup <- layer_membership[
  ,
  .N,
  by = cell_id
][
  N > 1
]

cat(
  "\nCells appearing in >1 count layer: ",
  nrow(layer_dup),
  "\n",
  sep = ""
)

# ============================================================
# 8. JOIN LAYER MEMBERSHIP TO CELL METADATA
# ============================================================

layer_meta <- merge(
  cell_meta,
  layer_membership,
  by = "cell_id",
  all.x = TRUE
)

unassigned_layer_cells <- sum(
  is.na(
    layer_meta$count_layer
  )
)

cat(
  "Cells without a count layer: ",
  unassigned_layer_cells,
  "\n",
  sep = ""
)

if (unassigned_layer_cells > 0) {
  stop(
    "Some cells do not map to an RNA count layer."
  )
}

# ============================================================
# 9. LAYER x SAMPLE AUDIT
# ============================================================

cat("\n============================================================\n")
cat("8. COUNT LAYER x SAMPLE AUDIT\n")
cat("============================================================\n\n")

layer_sample <- layer_meta[
  ,
  .(
    all_cells = .N,

    epithelial_cells =
      sum(
        broad_celltype_final ==
          "Epithelial"
      ),

    tumor_epithelial =
      sum(
        broad_celltype_final ==
          "Epithelial" &
        tissue_std ==
          "Tumor"
      ),

    normal_epithelial =
      sum(
        broad_celltype_final ==
          "Epithelial" &
        tissue_std ==
          "Normal"
      )
  ),
  by = .(
    count_layer,
    sample,
    tissue_std,
    treatment_group
  )
]

setorder(
  layer_sample,
  sample,
  count_layer
)

print(
  layer_sample,
  nrows = 200
)

fwrite(
  layer_sample,
  file.path(
    outdir,
    "STEP17C_GSE197677_count_layer_sample_map.csv"
  ),
  bom = TRUE
)

# ============================================================
# 10. SAMPLE -> NUMBER OF COUNT LAYERS
# ============================================================

sample_layer_check <- layer_meta[
  ,
  .(
    cells = .N,

    epithelial_cells =
      sum(
        broad_celltype_final ==
          "Epithelial"
      ),

    count_layers =
      uniqueN(
        count_layer
      ),

    count_layer_names =
      paste(
        sort(
          unique(
            count_layer
          )
        ),
        collapse = "; "
      )
  ),
  by = .(
    sample,
    tissue_std,
    treatment_group
  )
]

setorder(
  sample_layer_check,
  sample
)

cat("\n============================================================\n")
cat("9. SAMPLE -> COUNT LAYER MAPPING\n")
cat("============================================================\n\n")

print(
  sample_layer_check,
  nrows = 100
)

fwrite(
  sample_layer_check,
  file.path(
    outdir,
    "STEP17C_GSE197677_sample_to_count_layers.csv"
  ),
  bom = TRUE
)

# ============================================================
# 11. RAW COUNT EXTRACTION TEST
# ============================================================

cat("\n============================================================\n")
cat("10. SAMPLE-WISE RAW COUNT EXTRACTION TEST\n")
cat("============================================================\n\n")

epi_meta <- layer_meta[
  broad_celltype_final ==
    "Epithelial"
]

if (nrow(epi_meta) == 0) {
  stop(
    "No epithelial cells found."
  )
}

test_row <- epi_meta[
  tissue_std == "Tumor"
][1]

if (nrow(test_row) == 0) {
  test_row <- epi_meta[1]
}

test_layer <- test_row$count_layer
test_sample <- test_row$sample

test_cells <- epi_meta[
  sample == test_sample &
  count_layer == test_layer,
  cell_id
]

test_cells <- head(
  test_cells,
  100
)

cat(
  "Test sample : ",
  test_sample,
  "\n",
  sep = ""
)

cat(
  "Test layer  : ",
  test_layer,
  "\n",
  sep = ""
)

cat(
  "Test cells  : ",
  length(test_cells),
  "\n",
  sep = ""
)

test_counts <- tryCatch(
  SeuratObject::LayerData(
    rna,
    layer = test_layer,
    cells = test_cells,
    fast = FALSE
  ),
  error = function(e) {

    cat(
      "LayerData extraction failed:\n",
      conditionMessage(e),
      "\n"
    )

    NULL
  }
)

if (is.null(test_counts)) {
  stop(
    "Sample-wise count extraction test failed."
  )
}

cat(
  "Extracted matrix: ",
  nrow(test_counts),
  " genes x ",
  ncol(test_counts),
  " cells\n",
  sep = ""
)

cat(
  "Matrix class: ",
  paste(
    class(test_counts),
    collapse = ", "
  ),
  "\n",
  sep = ""
)

cat(
  "Non-zero values: ",
  Matrix::nnzero(test_counts),
  "\n",
  sep = ""
)

if (
  nrow(test_counts) == 0 ||
  ncol(test_counts) == 0 ||
  Matrix::nnzero(test_counts) == 0
) {
  stop(
    "Count extraction returned an empty/zero matrix."
  )
}

cat(
  "\n✓ RAW COUNTS ARE ACCESSIBLE SAMPLE-WISE\n"
)

# ============================================================
# 12. GENE FORMAT / hg38 OVERLAP
# ============================================================

cat("\n============================================================\n")
cat("11. GENE FORMAT / hg38 OVERLAP\n")
cat("============================================================\n\n")

genes <- rownames(
  test_counts
)

ensg_pct <- 100 *
  mean(
    grepl(
      "^ENSG[0-9]+",
      genes,
      ignore.case = TRUE
    )
  )

symbol_overlap <- sum(
  unique(genes) %in%
    unique(genome$symbol)
)

symbol_overlap_pct <- 100 *
  symbol_overlap /
  length(
    unique(genes)
  )

ensembl_overlap <- NA_integer_
ensembl_overlap_pct <- NA_real_

if (
  "ensembl_gene_id" %in%
    names(genome)
) {

  stripped_genes <- sub(
    "\\.[0-9]+$",
    "",
    genes
  )

  ensembl_overlap <- sum(
    unique(stripped_genes) %in%
      unique(
        genome$ensembl_gene_id
      )
  )

  ensembl_overlap_pct <- 100 *
    ensembl_overlap /
    length(
      unique(stripped_genes)
    )
}

cat(
  sprintf(
    "ENSG-like input genes : %.2f%%\n",
    ensg_pct
  )
)

cat(
  sprintf(
    "hg38 symbol overlap   : %d / %d (%.2f%%)\n",
    symbol_overlap,
    length(unique(genes)),
    symbol_overlap_pct
  )
)

if (!is.na(ensembl_overlap_pct)) {

  cat(
    sprintf(
      "hg38 Ensembl overlap  : %d / %d (%.2f%%)\n",
      ensembl_overlap,
      length(unique(genes)),
      ensembl_overlap_pct
    )
  )
}

gene_audit <- data.table(
  dataset =
    "GSE197677",

  genes =
    length(
      unique(genes)
    ),

  ENSG_like_pct =
    round(
      ensg_pct,
      3
    ),

  hg38_symbol_overlap =
    symbol_overlap,

  hg38_symbol_overlap_pct =
    round(
      symbol_overlap_pct,
      3
    ),

  hg38_ensembl_overlap =
    ensembl_overlap,

  hg38_ensembl_overlap_pct =
    round(
      ensembl_overlap_pct,
      3
    )
)

fwrite(
  gene_audit,
  file.path(
    outdir,
    "STEP17C_GSE197677_gene_hg38_overlap.csv"
  ),
  bom = TRUE
)

# ============================================================
# 13. NORMAL REFERENCE GROUPS
# ============================================================

cat("\n============================================================\n")
cat("12. NORMAL EPITHELIAL REFERENCE GROUPS\n")
cat("============================================================\n\n")

normal_refs <- layer_meta[
  broad_celltype_final ==
    "Epithelial" &
  tissue_std ==
    "Normal",
  .(
    normal_epithelial_cells =
      .N,

    count_layers =
      uniqueN(
        count_layer
      ),

    count_layer_names =
      paste(
        sort(
          unique(count_layer)
        ),
        collapse = "; "
      )
  ),
  by = sample
]

setorder(
  normal_refs,
  -normal_epithelial_cells
)

print(normal_refs)

if (nrow(normal_refs) < 2) {
  stop(
    "Fewer than two normal epithelial reference samples."
  )
}

cat(
  "\n✓ ",
  nrow(normal_refs),
  " separate normal epithelial reference groups available\n",
  sep = ""
)

fwrite(
  normal_refs,
  file.path(
    outdir,
    "STEP17C_GSE197677_normal_reference_groups.csv"
  ),
  bom = TRUE
)

# ============================================================
# 14. CHOOSE FIRST PILOT TUMOR SAMPLE
# ============================================================

tumor_candidates <- layer_meta[
  broad_celltype_final ==
    "Epithelial" &
  tissue_std ==
    "Tumor",
  .(
    epithelial_cells = .N,
    count_layers =
      uniqueN(
        count_layer
      )
  ),
  by = .(
    sample,
    treatment_group
  )
][
  epithelial_cells >= 30
]

setorder(
  tumor_candidates,
  epithelial_cells
)

cat("\n============================================================\n")
cat("13. PILOT TUMOR SAMPLE CANDIDATES\n")
cat("============================================================\n\n")

print(tumor_candidates)

if (nrow(tumor_candidates) == 0) {
  stop(
    "No tumor epithelial sample has >=30 cells."
  )
}

pilot_sample <- tumor_candidates$sample[1]

cat(
  "\nSelected first pilot sample: ",
  pilot_sample,
  "\n",
  sep = ""
)

cat(
  "Reason: smallest tumor epithelial compartment with >=30 cells.\n"
)

pilot_plan <- data.table(
  pilot_dataset =
    "GSE197677",

  pilot_tumor_sample =
    pilot_sample,

  pilot_tumor_epithelial_cells =
    tumor_candidates$epithelial_cells[1],

  normal_reference_samples =
    paste(
      normal_refs$sample,
      collapse = "; "
    ),

  normal_reference_groups =
    nrow(normal_refs),

  extraction_mode =
    "sample-wise LayerData; NO global JoinLayers",

  cnv_run_status =
    "NOT_RUN"
)

fwrite(
  pilot_plan,
  file.path(
    outdir,
    "STEP17C_GSE197677_pilot_plan.csv"
  ),
  bom = TRUE
)

# ============================================================
# 15. SAVE SESSION INFO
# ============================================================

sink(
  "08_logs/17C_sessionInfo.txt"
)

print(
  sessionInfo()
)

sink()

# ============================================================
# CLEANUP
# ============================================================

rm(
  test_counts,
  obj,
  rna
)

gc(verbose = FALSE)

# ============================================================
# FINAL
# ============================================================

cat("\n============================================================\n")
cat("STEP 17C FINISHED ✓\n")
cat("============================================================\n\n")

cat("infercna installed / verified ✓\n")
cat("hg38 reference verified ✓\n")
cat("GSE197677 count layers mapped ✓\n")
cat("Sample-wise raw-count extraction tested ✓\n")
cat("Normal epithelial reference groups identified ✓\n")
cat("First pilot tumor sample selected ✓\n")

cat("\nNO JoinLayers performed ✓\n")
cat("NO expression normalization performed ✓\n")
cat("NO CNV inference performed ✓\n")
cat("NO malignant labels assigned ✓\n")
cat("NO cells filtered ✓\n")

cat("\nSTOP HERE.\n")
cat("DO NOT RUN CNV YET.\n")
