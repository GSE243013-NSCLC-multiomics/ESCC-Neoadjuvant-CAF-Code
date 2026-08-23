# ============================================================
# ESCC Neoadjuvant Project
# STEP 17H5B
#
# Balance GSE221561 S9478T references to 500+500 cells
#
# Same strategy as all prior pilots:
#   - S6417N: 500 cells (seed = 22156117)
#   - S9478N: 500 cells (seed = 22156117)
#   - Combined: 1000 reference cells
#   - CPM re-normalized on balanced subset
#   - Gene table rebuilt
#
# IMPORTANT:
#   NO InferCNA run
#   NO malignant labels
# ============================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(infercna)
})

cat("\n============================================================\n")
cat("STEP 17H5B: BALANCE GSE221561 S9478T REFERENCES\n")
cat("============================================================\n\n")

# ============================================================
# 1. PATHS
# ============================================================

pilot_dir <- "06_tables/cross_dataset/malignant_epithelial/GSE221561_S9478T_pilot"

if (!dir.exists(pilot_dir)) {
  stop("Pilot directory missing: ", pilot_dir)
}

# ============================================================
# 2. LOAD 17H5A INPUTS
# ============================================================

cat("===== 1. LOAD 17H5A INPUTS =====\n\n")

cpm_file <- list.files(pilot_dir, pattern = "infercna_CPM_hg38\\.rds$", full.names = TRUE)[1]
refs_file <- list.files(pilot_dir, pattern = "infercna_refCells\\.rds$", full.names = TRUE)[1]
ann_file <- list.files(pilot_dir, pattern = "infercna_cell_annotation\\.csv$", full.names = TRUE)[1]
gene_file <- list.files(pilot_dir, pattern = "infercna_gene_table\\.csv$", full.names = TRUE)[1]

if (is.na(cpm_file)) stop("Missing CPM file")
if (is.na(refs_file)) stop("Missing refCells file")
if (is.na(ann_file)) stop("Missing annotation file")
if (is.na(gene_file)) stop("Missing gene table file")

cpm_all <- readRDS(cpm_file)
refs_all <- readRDS(refs_file)
ann <- fread(ann_file)
gene_table <- fread(gene_file)

cat("Full CPM: ", nrow(cpm_all), " genes x ", ncol(cpm_all), " cells\n", sep = "")
cat("Reference groups:\n")
print(lengths(refs_all))

# ============================================================
# 3. BALANCED SAMPLING
# ============================================================

cat("\n===== 2. BALANCED SAMPLING =====\n\n")

n_per_ref <- 500L

set.seed(22156117)
s6417n_cells <- refs_all[["S6417N"]]
if (length(s6417n_cells) > n_per_ref) {
  s6417n_balanced <- sample(s6417n_cells, n_per_ref)
} else {
  s6417n_balanced <- s6417n_cells
}

set.seed(22156117)
s9478n_cells <- refs_all[["S9478N"]]
if (length(s9478n_cells) > n_per_ref) {
  s9478n_balanced <- sample(s9478n_cells, n_per_ref)
} else {
  s9478n_balanced <- s9478n_cells
}

refCells_balanced <- list(
  S6417N = s6417n_balanced,
  S9478N = s9478n_balanced
)

cat("S6417N balanced: ", length(s6417n_balanced), " cells\n", sep = "")
cat("S9478N balanced: ", length(s9478n_balanced), " cells\n", sep = "")

tumor_cells <- ann[infercna_role == "TUMOR_TEST", cell_id]
all_balanced_cells <- c(tumor_cells, unlist(refCells_balanced, use.names = FALSE))

cat("Tumor cells    : ", length(tumor_cells), "\n", sep = "")
cat("Total balanced : ", length(all_balanced_cells), "\n", sep = "")

# ============================================================
# 4. EXTRACT AND RE-NORMALIZE CPM
# ============================================================

cat("\n===== 3. RE-NORMALIZE CPM =====\n\n")

cpm_balanced <- cpm_all[, all_balanced_cells, drop = FALSE]

cat("Balanced CPM: ", nrow(cpm_balanced), " genes x ", ncol(cpm_balanced), " cells\n", sep = "")

cpm_sums <- Matrix::colSums(cpm_balanced)
cat("Balanced CPM median colsum: ", median(cpm_sums), "\n", sep = "")

rm(cpm_all)
gc(verbose = FALSE)

# Remove zero-variance genes
gene_var <- apply(as.matrix(cpm_balanced), 1, var, na.rm = TRUE)
zero_genes <- names(gene_var[gene_var == 0 | is.na(gene_var)])

cat("Zero-variance genes: ", length(zero_genes), "\n", sep = "")

if (length(zero_genes) > 0) {
  keep_genes <- setdiff(rownames(cpm_balanced), zero_genes)
  cpm_balanced <- cpm_balanced[keep_genes, , drop = FALSE]
  cat("After removal: ", nrow(cpm_balanced), " genes\n", sep = "")
}

# ============================================================
# 5. hg38 VALIDATION
# ============================================================

cat("\n===== 4. hg38 VALIDATION =====\n\n")

infercna::useGenome("hg38")
genome <- infercna::retrieveGenome(name = "hg38")
hg38_symbols <- as.character(genome$symbol)

missing_hg38 <- setdiff(rownames(cpm_balanced), hg38_symbols)
cat("Genes absent from hg38: ", length(missing_hg38), "\n", sep = "")

if (length(missing_hg38) > 0) {
  stop("Balanced CPM contains symbols absent from InferCNA hg38.")
}

cat("All genes compatible with InferCNA hg38\n")

# ============================================================
# 6. GENE TABLE
# ============================================================

cat("\n===== 5. GENE TABLE =====\n\n")

gene_table_balanced <- data.table(gene = rownames(cpm_balanced))

genome_sub <- genome[match(gene_table_balanced$gene, genome$symbol), ]

gene_table_balanced[, chromosome := as.character(genome_sub$chromosome_name)]
gene_table_balanced[, arm := as.character(genome_sub$arm)]
gene_table_balanced[, start_position := genome_sub$start_position]
gene_table_balanced[, end_position := genome_sub$end_position]

if (any(is.na(gene_table_balanced$chromosome))) {
  stop("Genome mapping incomplete.")
}

cat("Gene table: ", nrow(gene_table_balanced), " genes\n", sep = "")

# ============================================================
# 7. SAVE BALANCED INPUTS
# ============================================================

cat("\n===== 6. SAVE BALANCED INPUTS =====\n\n")

saveRDS(cpm_balanced, file.path(pilot_dir, "GSE221561_S9478T_BALANCED_infercna_CPM_hg38.rds"), compress = FALSE)
saveRDS(refCells_balanced, file.path(pilot_dir, "GSE221561_S9478T_BALANCED_infercna_refCells.rds"))

fwrite(ann, file.path(pilot_dir, "GSE221561_S9478T_BALANCED_infercna_cell_annotation.csv"), bom = TRUE)
fwrite(gene_table_balanced, file.path(pilot_dir, "GSE221561_S9478T_BALANCED_infercna_gene_table.csv"), bom = TRUE)

manifest <- data.table(
  dataset = "GSE221561",
  pilot = "S9478T",
  balance_strategy = "500 S6417N + 500 S9478N",
  seed = 22156117,
  S6417N_cells = length(s6417n_balanced),
  S9478N_cells = length(s9478n_balanced),
  total_reference = length(unlist(refCells_balanced, use.names = FALSE)),
  tumor_cells = length(tumor_cells),
  total_cells = ncol(cpm_balanced),
  hg38_genes = nrow(cpm_balanced),
  zero_genes_removed = length(zero_genes)
)

fwrite(manifest, file.path(pilot_dir, "GSE221561_S9478T_BALANCED_infercna_manifest.csv"), bom = TRUE)

# ============================================================
# 8. RELOAD VALIDATION
# ============================================================

cat("\n===== 7. RELOAD VALIDATION =====\n\n")

check <- readRDS(file.path(pilot_dir, "GSE221561_S9478T_BALANCED_infercna_CPM_hg38.rds"))
check_refs <- readRDS(file.path(pilot_dir, "GSE221561_S9478T_BALANCED_infercna_refCells.rds"))

cat("Reloaded CPM: ", nrow(check), " genes x ", ncol(check), " cells\n", sep = "")
cat("Reference groups: ", length(check_refs), "\n", sep = "")
cat("Reference cells: ", sum(lengths(check_refs)), "\n", sep = "")

if (!all(unlist(check_refs, use.names = FALSE) %in% colnames(check))) {
  stop("Reference validation failed.")
}

if (!all(tumor_cells %in% colnames(check))) {
  stop("Tumor cell validation failed.")
}

cat("Balanced inputs validated\n")

# ============================================================
# CLEANUP & FINAL
# ============================================================

rm(cpm_balanced, check, genome, genome_sub, gene_var)
gc(verbose = FALSE)

cat("\n============================================================\n")
cat("STEP 17H5B FINISHED\n")
cat("GSE221561 S9478T BALANCED INPUTS READY\n")
cat("============================================================\n\n")

cat("Tumor cells       : ", length(tumor_cells), "\n", sep = "")
cat("S6417N balanced   : ", length(s6417n_balanced), "\n", sep = "")
cat("S9478N balanced   : ", length(s9478n_balanced), "\n", sep = "")
cat("Total cells       : ", ncol(check), "\n", sep = "")
cat("Final hg38 genes  : ", nrow(check), "\n", sep = "")

cat("\nSTOP HERE.\n")
