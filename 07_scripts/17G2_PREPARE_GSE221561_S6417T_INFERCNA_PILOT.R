# ============================================================
# ESCC Neoadjuvant Project
# STEP 17G2
#
# Prepare GSE221561 S6417T InferCNA pilot inputs.
#
# Strategy:
#   - S6417T tumor epithelial (Surgery_alone, 454 cells)
#   - S6417N + S9478N normal epithelial references
#   - Extract from split layers WITHOUT JoinLayers
#   - Resolve 9 duplicate hg38 target symbols
#   - Same InferCNA parameters as GSE197677 pilots
#
# IMPORTANT:
#   NO InferCNA run
#   NO JoinLayers
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
cat("STEP 17G2: PREPARE GSE221561 S6417T INFERCNA PILOT\n")
cat("============================================================\n\n")

# ============================================================
# 1. PATHS
# ============================================================

object_dir <- "03_objects/GSE221561"

meta_file <- paste0(
  "06_tables/GSE221561/metadata/",
  "GSE221561_cell_metadata_standardized.csv.gz"
)

preflight_dir <- paste0(
  "06_tables/cross_dataset/",
  "malignant_epithelial/",
  "GSE221561_preflight"
)

base_outdir <- paste0(
  "06_tables/cross_dataset/",
  "malignant_epithelial"
)

pilot_dir <- file.path(
  base_outdir,
  "GSE221561_S6417T_pilot"
)

dir.create(
  pilot_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

if (!file.exists(meta_file)) {
  stop("Missing: ", meta_file)
}

# ============================================================
# 2. LOAD METADATA
# ============================================================

cat("===== 1. LOAD METADATA =====\n\n")

md <- fread(meta_file)

md[
  ,
  tissue_simple :=
    fcase(
      grepl("tumor|tumour|scc", tissue_std, ignore.case = TRUE),
      "Tumor",
      grepl("normal|adjacent", tissue_std, ignore.case = TRUE),
      "Normal",
      default = as.character(tissue_std)
    )
]

cat("Metadata cells: ", nrow(md), "\n", sep = "")

# ============================================================
# 3. DEFINE PILOT AND REFERENCES
# ============================================================

cat("\n===== 2. DEFINE PILOT AND REFERENCES =====\n\n")

pilot_sample <- "S6417T"
pilot_patient <- "P6417"
normal_samples <- c("S6417N", "S9478N")

pilot_cells <- md[
  sample_id == pilot_sample &
  tissue_simple == "Tumor" &
  cell_class_std == "Epithelial_unclassified",
  cell_id
]

normal_meta <- md[
  sample_id %in% normal_samples &
  tissue_simple == "Normal" &
  cell_class_std == "Epithelial_unclassified"
]

normal_cells <- normal_meta$cell_id

cat("Pilot sample       : ", pilot_sample, "\n", sep = "")
cat("Pilot patient      : ", pilot_patient, "\n", sep = "")
cat("Pilot tumor cells  : ", length(pilot_cells), "\n", sep = "")
cat("Normal samples     : ", paste(normal_samples, collapse = ", "), "\n", sep = "")
cat("Normal ref cells   : ", length(normal_cells), "\n", sep = "")

normal_summary <- normal_meta[
  ,
  .(epithelial_cells = .N),
  by = sample_id
]
setorder(normal_summary, sample_id)
print(normal_summary)

# ============================================================
# 4. LOCATE AND LOAD SEURAT OBJECT
# ============================================================

cat("\n===== 3. LOAD SEURAT OBJECT =====\n\n")

rds_files <- list.files(
  object_dir,
  pattern = "\\.rds$",
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)

object_file <- rds_files[
  grepl("author_annotated", basename(rds_files), ignore.case = TRUE)
][1]

if (is.na(object_file)) {
  object_file <- rds_files[1]
}

cat("Object: ", object_file, "\n", sep = "")

obj <- readRDS(object_file)

cat("Object cells: ", ncol(obj), "\n", sep = "")
cat("Assays: ", paste(Seurat::Assays(obj), collapse = ", "), "\n", sep = "")

rna <- obj@assays[["RNA"]]

cat("RNA features: ", nrow(rna), "\n", sep = "")

# ============================================================
# 5. MAP CELLS TO COUNT LAYERS
# ============================================================

cat("\n===== 4. MAP CELLS TO COUNT LAYERS =====\n\n")

layers <- tryCatch(
  SeuratObject::Layers(rna, search = NA),
  error = function(e) character(0)
)

count_layers <- layers[
  grepl("^counts($|\\.)", layers, ignore.case = TRUE)
]

cat("Count layers: ", paste(count_layers, collapse = ", "), "\n", sep = "")

get_layer_cells <- function(assay_object, layer_name) {
  x <- tryCatch(
    SeuratObject::Cells(assay_object, layer = layer_name),
    error = function(e) NULL
  )
  if (!is.null(x)) return(x)
  m <- SeuratObject::LayerData(assay_object, layer = layer_name, fast = FALSE)
  colnames(m)
}

layer_membership <- rbindlist(lapply(count_layers, function(lyr) {
  lc <- get_layer_cells(rna, lyr)
  cat(sprintf("  %-20s %6d cells\n", lyr, length(lc)))
  data.table(cell_id = lc, count_layer = lyr)
}), use.names = TRUE)

selected_cells <- c(pilot_cells, normal_cells)
missing_cells <- setdiff(selected_cells, layer_membership$cell_id)

cat("\nSelected cells missing from layers: ", length(missing_cells), "\n", sep = "")

if (length(missing_cells) > 0) {
  stop("Selected cells not fully represented in count layers.")
}

# ============================================================
# 6. LOAD hg38 AND FEATURE MAPPING
# ============================================================

cat("\n===== 5. LOAD hg38 AND FEATURE MAPPING =====\n\n")

infercna::useGenome("hg38")
genome <- infercna::retrieveGenome(name = "hg38")

mapping_file <- file.path(
  preflight_dir,
  "STEP17G1_GSE221561_feature_hg38_mapping.csv"
)

if (!file.exists(mapping_file)) {
  stop("Missing: ", mapping_file)
}

gene_map <- fread(mapping_file)

mapped <- gene_map[
  mapping_type != "UNMAPPED" &
  !is.na(hg38_symbol)
]

# Resolve duplicates: sum counts for duplicate target symbols
dup_targets <- mapped[, .(input_rows = .N), by = hg38_symbol][input_rows > 1]

cat("Total mapped input rows: ", nrow(mapped), "\n", sep = "")
cat("Unique hg38 symbols: ", uniqueN(mapped$hg38_symbol), "\n", sep = "")
cat("Duplicate targets: ", nrow(dup_targets), "\n", sep = "")

# ============================================================
# 7. EXTRACT RAW COUNTS SAMPLE-WISE
# ============================================================

cat("\n===== 6. EXTRACT RAW COUNTS SAMPLE-WISE =====\n\n")

extract_counts <- function(cells, layer_membership, rna) {

  cell_layers <- layer_membership[cell_id %in% cells]

  parts <- lapply(unique(cell_layers$count_layer), function(lyr) {
    lyr_cells <- cell_layers[count_layer == lyr, cell_id]
    cat(sprintf("  Extracting %d cells from %s\n", length(lyr_cells), lyr))
    SeuratObject::LayerData(rna, layer = lyr, cells = lyr_cells, fast = FALSE)
  })

  if (length(parts) == 1) {
    return(parts[[1]])
  }

  common_genes <- Reduce(intersect, lapply(parts, rownames))
  do.call(cbind, lapply(parts, function(p) p[common_genes, , drop = FALSE]))
}

cat("Pilot tumor cells (", length(pilot_cells), "):\n", sep = "")
pilot_raw <- extract_counts(pilot_cells, layer_membership, rna)

cat("\nNormal cells (", length(normal_cells), "):\n", sep = "")
normal_raw <- extract_counts(normal_cells, layer_membership, rna)

cat("\nPilot raw: ", nrow(pilot_raw), " genes x ", ncol(pilot_raw), " cells\n", sep = "")
cat("Normal raw: ", nrow(normal_raw), " genes x ", ncol(normal_raw), " cells\n", sep = "")

# ============================================================
# 8. MAP TO hg38 (separate for pilot and normal)
# ============================================================

cat("\n===== 7. MAP TO hg38 =====\n\n")

mapped_lookup <- mapped[, .(input_feature, hg38_symbol)]
mapped_lookup <- mapped_lookup[!duplicated(mapped_lookup)]

map_and_deduplicate <- function(raw_mat, label) {

  cat("\nMapping ", label, "...\n", sep = "")

  genes <- rownames(raw_mat)
  mapped_idx <- match(genes, mapped_lookup$input_feature)
  mapped_syms <- mapped_lookup$hg38_symbol[mapped_idx]
  has_symbol <- !is.na(mapped_syms)

  cat("  Genes with hg38 symbol: ", sum(has_symbol), "/", length(genes), "\n", sep = "")

  mapped_mat <- raw_mat[has_symbol, , drop = FALSE]
  rownames(mapped_mat) <- mapped_syms[has_symbol]

  dup_rows <- which(duplicated(rownames(mapped_mat)))
  if (length(dup_rows) > 0) {
    cat("  Summing ", length(dup_rows), " duplicate rows\n", sep = "")

    dup_names <- rownames(mapped_mat)[dup_rows]
    for (dn in unique(dup_names)) {
      idx <- which(rownames(mapped_mat) == dn)
      mapped_mat[idx[1], ] <- Matrix::colSums(mapped_mat[idx, , drop = FALSE])
      mapped_mat <- mapped_mat[-idx[-1], , drop = FALSE]
    }
  }

  cat("  Final: ", nrow(mapped_mat), " genes x ", ncol(mapped_mat), " cells\n", sep = "")
  mapped_mat
}

pilot_mapped <- map_and_deduplicate(pilot_raw, "pilot")
normal_mapped <- map_and_deduplicate(normal_raw, "normal")

# ============================================================
# 9. CPM NORMALIZATION (separate, then combine)
# ============================================================

cat("\n===== 8. CPM NORMALIZATION =====\n\n")

cpm_pilot <- t(t(pilot_mapped) / Matrix::colSums(pilot_mapped) * 1e6)
cpm_normal <- t(t(normal_mapped) / Matrix::colSums(normal_mapped) * 1e6)

cat("Pilot CPM median colsum: ", median(Matrix::colSums(cpm_pilot)), "\n", sep = "")
cat("Normal CPM median colsum: ", median(Matrix::colSums(cpm_normal)), "\n", sep = "")

# Combine CPM matrices (cells are columns, genes are rows)
common_genes_cpm <- intersect(rownames(cpm_pilot), rownames(cpm_normal))
cat("Common CPM genes: ", length(common_genes_cpm), "\n", sep = "")

cpm <- cbind(
  cpm_pilot[common_genes_cpm, , drop = FALSE],
  cpm_normal[common_genes_cpm, , drop = FALSE]
)

cat("Combined CPM: ", nrow(cpm), " genes x ", ncol(cpm), " cells\n", sep = "")

cpm_sums <- Matrix::colSums(cpm)
cat("Final CPM median colsum: ", median(cpm_sums), "\n", sep = "")

if (max(abs(cpm_sums - 1e6)) > 1) {
  stop("CPM validation failed.")
}

cat("CPM validated\n")

# Save combined raw counts (hg38-mapped, before CPM)
raw_file <- file.path(pilot_dir, "GSE221561_S6417T_infercna_raw_counts_hg38.rds")
common_genes_raw <- intersect(rownames(pilot_mapped), rownames(normal_mapped))
combined_raw <- cbind(
  pilot_mapped[common_genes_raw, , drop = FALSE],
  normal_mapped[common_genes_raw, , drop = FALSE]
)
saveRDS(combined_raw, raw_file, compress = FALSE)
cat("Saved raw counts: ", nrow(combined_raw), " genes x ", ncol(combined_raw), " cells\n", sep = "")
rm(combined_raw)
gc(verbose = FALSE)

# ============================================================
# 10. BUILD REFERENCE LIST
# ============================================================

cat("\n===== 9. BUILD REFERENCE LIST =====\n\n")

refCells <- lapply(normal_samples, function(s) {
  normal_meta[sample_id == s, cell_id]
})
names(refCells) <- normal_samples

cat("Reference groups:\n")
print(lengths(refCells))

# ============================================================
# 11. ANNOTATION
# ============================================================

cat("\n===== 10. ANNOTATION =====\n\n")

annotation <- md[
  cell_id %in% colnames(cpm),
  .(cell_id, sample_id, patient_id, tissue = tissue_simple, treatment_group, cell_class_std)
]

annotation[, infercna_role := fifelse(
  cell_id %in% pilot_cells,
  "TUMOR_TEST",
  "NORMAL_REFERENCE"
)]

annotation[, reference_group := fifelse(
  infercna_role == "NORMAL_REFERENCE",
  sample_id,
  NA_character_
)]

role_summary <- annotation[
  ,
  .(cells = .N),
  by = .(infercna_role, sample_id, treatment_group)
]

print(role_summary)

# ============================================================
# 12. GENE TABLE
# ============================================================

cat("\n===== 11. GENE TABLE =====\n\n")

gene_table <- data.table(gene = rownames(cpm))

genome_sub <- genome[match(gene_table$gene, genome$symbol), ]

gene_table[, chromosome := as.character(genome_sub$chromosome_name)]
gene_table[, arm := as.character(genome_sub$arm)]
gene_table[, start_position := genome_sub$start_position]
gene_table[, end_position := genome_sub$end_position]

if (any(is.na(gene_table$chromosome))) {
  stop("Genome mapping incomplete.")
}

cat("Gene table: ", nrow(gene_table), " genes\n", sep = "")

# ============================================================
# 13. SAVE INPUTS
# ============================================================

cat("\n===== 12. SAVE INPUTS =====\n\n")

cpm_file <- file.path(pilot_dir, "GSE221561_S6417T_infercna_CPM_hg38.rds")
raw_file <- file.path(pilot_dir, "GSE221561_S6417T_infercna_raw_counts_hg38.rds")
refs_file <- file.path(pilot_dir, "GSE221561_S6417T_infercna_refCells.rds")
annotation_file <- file.path(pilot_dir, "GSE221561_S6417T_infercna_cell_annotation.csv")
genes_file <- file.path(pilot_dir, "GSE221561_S6417T_infercna_gene_table.csv")
manifest_file <- file.path(pilot_dir, "GSE221561_S6417T_infercna_manifest.csv")

saveRDS(cpm, cpm_file, compress = FALSE)
saveRDS(refCells, refs_file)
fwrite(annotation, annotation_file, bom = TRUE)
fwrite(gene_table, genes_file, bom = TRUE)

manifest <- data.table(
  dataset = "GSE221561",
  pilot_number = 1L,
  pilot_sample = pilot_sample,
  pilot_patient = pilot_patient,
  selection_rule = "largest Surgery_alone epithelial sample",
  treatment_group = "Surgery_alone",
  tumor_epithelial_cells = length(pilot_cells),
  normal_reference_groups = length(refCells),
  normal_reference_cells = length(normal_cells),
  total_cells = ncol(cpm),
  hg38_genes = nrow(cpm),
  normalization = "CPM",
  isLog = FALSE,
  planned_window = 100,
  planned_n = 5000,
  planned_noise = 0.1,
  planned_center_method = "median",
  parameter_change_from_GSE197677 = FALSE,
  infercna_run = FALSE
)

fwrite(manifest, manifest_file, bom = TRUE)

# ============================================================
# 14. RELOAD VALIDATION
# ============================================================

cat("\n===== 13. RELOAD VALIDATION =====\n\n")

check <- readRDS(cpm_file)
check_refs <- readRDS(refs_file)

cat("Reloaded CPM: ", nrow(check), " genes x ", ncol(check), " cells\n", sep = "")
cat("Reference groups: ", length(check_refs), "\n", sep = "")

if (!all(unlist(check_refs, use.names = FALSE) %in% colnames(check))) {
  stop("Reference validation failed.")
}

if (!all(pilot_cells %in% colnames(check))) {
  stop("Tumor cell validation failed.")
}

# ============================================================
# CLEANUP
# ============================================================

rm(obj, rna, pilot_raw, normal_raw, pilot_mapped, normal_mapped, cpm_pilot, cpm_normal, check)
gc(verbose = FALSE)

# ============================================================
# FINAL
# ============================================================

cat("\n============================================================\n")
cat("STEP 17G2 FINISHED\n")
cat("GSE221561 S6417T PILOT INPUT READY\n")
cat("============================================================\n\n")

cat("Pilot sample      : ", pilot_sample, "\n", sep = "")
cat("Pilot treatment   : Surgery_alone\n")
cat("Tumor cells       : ", length(pilot_cells), "\n", sep = "")
cat("Normal ref cells  : ", length(normal_cells), "\n", sep = "")
cat("Normal ref groups : ", length(refCells), "\n", sep = "")
cat("Total pilot cells : ", ncol(cpm), "\n", sep = "")
cat("Final hg38 genes  : ", nrow(cpm), "\n", sep = "")

cat("\nSAME PARAMETERS AS GSE197677:\n")
cat("  window = 100\n")
cat("  n = 5000\n")
cat("  noise = 0.1\n")
cat("  center.method = median\n")
cat("  isLog = FALSE\n")

cat("\nNO InferCNA run\n")
cat("NO parameters changed\n")
cat("NO malignant labels assigned\n")
cat("NO cells filtered\n")

cat("\nSTOP HERE.\n")
