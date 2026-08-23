# ============================================================
# ESCC Neoadjuvant Project
# STEP 18D
#
# Prepare GSE160269 P39T InferCNA pilot inputs.
#
# Strategy:
#   - P39T tumor epithelial (2884 cells, highest cell count)
#   - P126N + P127N + P128N + P130N as 4 separate reference groups
#   - All 183 adjacent-normal epithelial cells as reference
#   - Same InferCNA parameters as all prior pilots
#
# IMPORTANT:
#   NO InferCNA run
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
cat("STEP 18D: PREPARE GSE160269 P39T INFERCNA PILOT\n")
cat("============================================================\n\n")

# ============================================================
# 1. PATHS
# ============================================================

object_file <- "03_objects/GSE160269/GSE160269_annotated_merged_raw.rds"

base_outdir <- "06_tables/cross_dataset/malignant_epithelial"

pilot_dir <- file.path(
  base_outdir,
  "GSE160269_P39T_pilot"
)

dir.create(
  pilot_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

if (!file.exists(object_file)) {
  stop("Missing: ", object_file)
}

# ============================================================
# 2. DEFINE PILOT AND REFERENCES
# ============================================================

cat("===== 1. DEFINE PILOT AND REFERENCES =====\n\n")

pilot_specimen <- "P39T"
pilot_patient <- "P39"
normal_specimens <- c("P126N", "P127N", "P128N", "P130N")

cat("Pilot specimen     : ", pilot_specimen, "\n", sep = "")
cat("Pilot patient      : ", pilot_patient, "\n", sep = "")
cat("Normal specimens   : ", paste(normal_specimens, collapse = ", "), "\n", sep = "")

# ============================================================
# 3. LOAD SEURAT OBJECT
# ============================================================

cat("\n===== 2. LOAD SEURAT OBJECT =====\n\n")

obj <- readRDS(object_file)

cat("Object cells: ", ncol(obj), "\n", sep = "")

# ============================================================
# 4. DEFINE PILOT AND NORMAL CELLS
# ============================================================

cat("\n===== 3. DEFINE PILOT AND NORMAL CELLS =====\n\n")

meta <- as.data.table(obj@meta.data, keep.rownames = "cell_barcode")

pilot_cells <- meta[
  specimen_id == pilot_specimen &
  tissue == "Tumor" &
  cell_type_author_source == "Epithelia",
  cell_barcode
]

cat("Pilot tumor cells found: ", length(pilot_cells), "\n", sep = "")

if (length(pilot_cells) < 30) {
  stop("Pilot epithelial tumor cells < 30. STOPPING.")
}

normal_cells_by_specimen <- lapply(normal_specimens, function(s) {
  cells <- meta[
    specimen_id == s &
    tissue == "Adjacent_normal" &
    cell_type_author_source == "Epithelia",
    cell_barcode
  ]
  cat(sprintf("  %s: %d cells\n", s, length(cells)))
  cells
})
names(normal_cells_by_specimen) <- normal_specimens

all_normal_cells <- unlist(normal_cells_by_specimen, use.names = FALSE)
cat("Total normal reference cells: ", length(all_normal_cells), "\n", sep = "")

# ============================================================
# 5. EXTRACT RAW COUNTS
# ============================================================

cat("\n===== 4. EXTRACT RAW COUNTS =====\n\n")

rna <- obj@assays[["RNA"]]

layers <- tryCatch(
  SeuratObject::Layers(rna, search = NA),
  error = function(e) character(0)
)

count_layers <- layers[
  grepl("^counts($|\\.)", layers, ignore.case = TRUE)
]

cat("Count layers: ", length(count_layers), "\n", sep = "")

# Use the layer containing epithelial cells
# Find which layer has our cells
target_layer <- NULL
for (lyr in count_layers) {
  lyr_cells <- tryCatch(
    SeuratObject::Cells(rna, layer = lyr),
    error = function(e) NULL
  )
  if (is.null(lyr_cells)) {
    m <- SeuratObject::LayerData(rna, layer = lyr, fast = FALSE)
    lyr_cells <- colnames(m)
  }
  if (pilot_cells[1] %in% lyr_cells) {
    target_layer <- lyr
    cat("Found pilot cells in layer: ", lyr, "\n", sep = "")
    break
  }
}

if (is.null(target_layer)) {
  stop("Could not find pilot cells in any layer.")
}

# Extract from the identified layer
cat("Extracting pilot cells...\n")
pilot_raw <- SeuratObject::LayerData(rna, layer = target_layer, cells = pilot_cells, fast = FALSE)

cat("Extracting normal cells...\n")
normal_raw <- SeuratObject::LayerData(rna, layer = target_layer, cells = all_normal_cells, fast = FALSE)

cat("\nPilot raw: ", nrow(pilot_raw), " genes x ", ncol(pilot_raw), " cells\n", sep = "")
cat("Normal raw: ", nrow(normal_raw), " genes x ", ncol(normal_raw), " cells\n", sep = "")

# ============================================================
# 6. MAP TO hg38
# ============================================================

cat("\n===== 5. MAP TO hg38 =====\n\n")

hg38_dt <- as.data.table(infercna:::hg38)
hg38_genes <- hg38_dt$symbol

cat("InferCNA hg38 genes: ", length(hg38_genes), "\n", sep = "")

map_and_deduplicate <- function(raw_mat, label) {
  cat("\nMapping ", label, "...\n", sep = "")
  genes <- rownames(raw_mat)
  
  in_hg38 <- genes %in% hg38_genes
  cat("  Genes in hg38: ", sum(in_hg38), "/", length(genes), "\n", sep = "")
  
  mapped_mat <- raw_mat[in_hg38, , drop = FALSE]
  
  dup_genes <- rownames(mapped_mat)[duplicated(rownames(mapped_mat))]
  if (length(dup_genes) > 0) {
    cat("  Summing ", length(unique(dup_genes)), " duplicate gene groups\n", sep = "")
    for (dg in unique(dup_genes)) {
      idx <- which(rownames(mapped_mat) == dg)
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
# 7. CPM NORMALIZATION
# ============================================================

cat("\n===== 6. CPM NORMALIZATION =====\n\n")

cpm_pilot <- t(t(pilot_mapped) / Matrix::colSums(pilot_mapped) * 1e6)
cpm_normal <- t(t(normal_mapped) / Matrix::colSums(normal_mapped) * 1e6)

cat("Pilot CPM median colsum: ", median(Matrix::colSums(cpm_pilot)), "\n", sep = "")
cat("Normal CPM median colsum: ", median(Matrix::colSums(cpm_normal)), "\n", sep = "")

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

# ============================================================
# 8. BUILD REFERENCE LIST (4 SEPARATE GROUPS)
# ============================================================

cat("\n===== 7. BUILD REFERENCE LIST =====\n\n")

refCells <- normal_cells_by_specimen

cat("Reference groups:\n")
print(lengths(refCells))

# ============================================================
# 9. ANNOTATION
# ============================================================

cat("\n===== 8. ANNOTATION =====\n\n")

annotation <- meta[
  cell_barcode %in% colnames(cpm),
  .(cell_barcode, sample_id = specimen_id, patient_id, tissue, cell_type_author_source)
]

annotation[, infercna_role := fifelse(
  cell_barcode %in% pilot_cells,
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
  by = .(infercna_role, sample_id)
]

print(role_summary)

# ============================================================
# 10. GENE TABLE
# ============================================================

cat("\n===== 9. GENE TABLE =====\n\n")

gene_table <- data.table(gene = rownames(cpm))

genome_sub <- hg38_dt[match(gene_table$gene, hg38_dt$symbol), ]

gene_table[, chromosome := as.character(genome_sub$chromosome_name)]
gene_table[, arm := as.character(genome_sub$arm)]
gene_table[, start_position := genome_sub$start_position]
gene_table[, end_position := genome_sub$end_position]

if (any(is.na(gene_table$chromosome))) {
  n_na <- sum(is.na(gene_table$chromosome))
  cat("WARNING: ", n_na, " genes without chromosome info\n", sep = "")
  gene_table <- gene_table[!is.na(chromosome)]
  cpm <- cpm[gene_table$gene, , drop = FALSE]
}

cat("Gene table: ", nrow(gene_table), " genes\n", sep = "")

# ============================================================
# 11. SAVE INPUTS
# ============================================================

cat("\n===== 10. SAVE INPUTS =====\n\n")

cpm_file <- file.path(pilot_dir, "GSE160269_P39T_infercna_CPM_hg38.rds")
refs_file <- file.path(pilot_dir, "GSE160269_P39T_infercna_refCells.rds")
annotation_file <- file.path(pilot_dir, "GSE160269_P39T_infercna_cell_annotation.csv")
genes_file <- file.path(pilot_dir, "GSE160269_P39T_infercna_gene_table.csv")
manifest_file <- file.path(pilot_dir, "GSE160269_P39T_infercna_manifest.csv")

saveRDS(cpm, cpm_file, compress = FALSE)
saveRDS(refCells, refs_file)
fwrite(annotation, annotation_file, bom = TRUE)
fwrite(gene_table, genes_file, bom = TRUE)

manifest <- data.table(
  dataset = "GSE160269",
  pilot_number = 2L,
  pilot_sample = pilot_specimen,
  pilot_patient = pilot_patient,
  selection_rule = "Tumor epithelial, second pilot (2884 cells)",
  treatment_group = "Untreated_baseline",
  tumor_epithelial_cells = length(pilot_cells),
  normal_reference_groups = length(refCells),
  normal_reference_cells = length(all_normal_cells),
  total_cells = ncol(cpm),
  hg38_genes = nrow(cpm),
  normalization = "CPM",
  isLog = FALSE,
  planned_window = 100,
  planned_n = 5000,
  planned_noise = 0.1,
  planned_center_method = "median",
  reference_quality = "LIMITED_183_CELLS_LOWER_CONFIDENCE",
  infercna_run = FALSE
)

fwrite(manifest, manifest_file, bom = TRUE)

# ============================================================
# 12. RELOAD VALIDATION
# ============================================================

cat("\n===== 11. RELOAD VALIDATION =====\n\n")

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

cat("Validation PASSED\n")

# ============================================================
# CLEANUP
# ============================================================

rm(obj, rna, pilot_raw, normal_raw, pilot_mapped, normal_mapped, cpm_pilot, cpm_normal, check)
gc(verbose = FALSE)

# ============================================================
# FINAL
# ============================================================

cat("\n============================================================\n")
cat("STEP 18D FINISHED\n")
cat("GSE160269 P39T PILOT INPUT READY\n")
cat("============================================================\n\n")

cat("Pilot sample      : ", pilot_specimen, "\n", sep = "")
cat("Pilot treatment   : Untreated_baseline\n")
cat("Tumor cells       : ", length(pilot_cells), "\n", sep = "")
cat("Normal ref cells  : ", length(all_normal_cells), "\n", sep = "")
cat("Normal ref groups : ", length(refCells), "\n", sep = "")
cat("Total cells       : ", ncol(cpm), "\n", sep = "")
cat("Final hg38 genes  : ", nrow(cpm), "\n", sep = "")

cat("\nSAME PARAMETERS AS ALL PRIOR PILOTS:\n")
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
