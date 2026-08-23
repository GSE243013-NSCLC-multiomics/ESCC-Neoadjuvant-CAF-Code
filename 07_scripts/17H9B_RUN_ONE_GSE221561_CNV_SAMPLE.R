# ============================================================
# ESCC Neoadjuvant Project
# STEP 17H9B
#
# Reusable worker: run full GSE221561 CNV pipeline for ONE sample.
#
# Usage:
#   Rscript 17H9B_RUN_ONE_GSE221561_CNV_SAMPLE.R SAMPLE_ID
#
# Example:
#   Rscript 17H9B_RUN_ONE_GSE221561_CNV_SAMPLE.R S1535T
#
# Combines: 17H7A (prep) + 17H7B (balance) + 17H7C (InferCNA)
#
# Frozen parameters (identical to all prior pilots):
#   genome        = hg38
#   window        = 100
#   range         = c(-3, 3)
#   n             = 5000
#   noise         = 0.1
#   center.method = median
#   isLog         = FALSE
#   seed          = 22156117
#   normal refs   = 500 S6417N + 500 S9478N
#
# IMPORTANT:
#   DOES run InferCNA
#   DOES NOT run findMalignant()
#   DOES NOT assign malignant labels
#   DOES NOT filter tumor cells
# ============================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Usage: Rscript 17H9B_RUN_ONE_GSE221561_CNV_SAMPLE.R SAMPLE_ID")
}

sample_id <- args[1]
sample_var <- sample_id  # avoid column name shadowing in data.table

options(stringsAsFactors = FALSE, timeout = 3600)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(infercna)
  library(Matrix)
  library(data.table)
})

cat("\n============================================================\n")
cat("STEP 17H9B: GSE221561 SINGLE-SAMPLE CNV PIPELINE\n")
cat("Sample: ", sample_id, "\n", sep = "")
cat("============================================================\n\n")

# ============================================================
# A. RAW COUNTS PREPARATION (from 17H7A template)
# ============================================================

cat("===== PHASE A: RAW COUNTS PREPARATION =====\n\n")

# --- Paths ---
object_dir <- "03_objects/GSE221561"
meta_file <- "06_tables/GSE221561/metadata/GSE221561_cell_metadata_standardized.csv.gz"
preflight_dir <- "06_tables/cross_dataset/malignant_epithelial/GSE221561_preflight"
base_outdir <- "06_tables/cross_dataset/malignant_epithelial"
pilot_dir <- file.path(base_outdir, paste0("GSE221561_", sample_id, "_pilot"))
results_dir <- file.path(pilot_dir, "17H9_results")
figdir <- file.path("05_figures/GSE221561/malignant_epithelial", paste0(sample_id, "_pilot"))

dir.create(pilot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(meta_file)) stop("Missing: ", meta_file)

# --- Load metadata ---
md <- fread(meta_file)
md[, tissue_simple := fcase(
  grepl("tumor|tumour|scc", tissue_std, ignore.case = TRUE), "Tumor",
  grepl("normal|adjacent", tissue_std, ignore.case = TRUE), "Normal",
  default = as.character(tissue_std)
)]

# --- Define pilot and references ---
patient_id <- md[sample_id == sample_var, patient_id][1]
if (is.na(patient_id)) stop("Sample not found in metadata: ", sample_var)

treatment <- md[sample_id == sample_var, treatment_group][1]
cat("Sample     : ", sample_var, "\n", sep = "")
cat("Patient    : ", patient_id, "\n", sep = "")
cat("Treatment  : ", treatment, "\n", sep = "")

normal_samples <- c("S6417N", "S9478N")

pilot_cells <- md[
  sample_id == sample_var &
  tissue_simple == "Tumor" &
  cell_class_std == "Epithelial_unclassified",
  cell_id
]

cat("Tumor epithelial cells: ", length(pilot_cells), "\n", sep = "")

if (length(pilot_cells) < 30) {
  stop("Epithelial tumor cells < 30. Cannot run InferCNA.")
}

normal_meta <- md[
  sample_id %in% normal_samples &
  tissue_simple == "Normal" &
  cell_class_std == "Epithelial_unclassified"
]
normal_cells <- normal_meta$cell_id

cat("Normal ref cells: ", length(normal_cells), "\n", sep = "")

# --- Load Seurat object ---
rds_files <- list.files(object_dir, pattern = "\\.rds$", recursive = TRUE,
                        full.names = TRUE, ignore.case = TRUE)
object_file <- rds_files[grepl("author_annotated", basename(rds_files), ignore.case = TRUE)][1]
if (is.na(object_file)) object_file <- rds_files[1]

cat("Loading: ", object_file, "\n", sep = "")
obj <- readRDS(object_file)
rna <- obj@assays[["RNA"]]

# --- Map cells to count layers ---
layers <- tryCatch(SeuratObject::Layers(rna, search = NA), error = function(e) character(0))
count_layers <- layers[grepl("^counts($|\\.)", layers, ignore.case = TRUE)]

get_layer_cells <- function(assay_object, layer_name) {
  x <- tryCatch(SeuratObject::Cells(assay_object, layer = layer_name), error = function(e) NULL)
  if (!is.null(x)) return(x)
  m <- SeuratObject::LayerData(assay_object, layer = layer_name, fast = FALSE)
  colnames(m)
}

layer_membership <- rbindlist(lapply(count_layers, function(lyr) {
  lc <- get_layer_cells(rna, lyr)
  data.table(cell_id = lc, count_layer = lyr)
}), use.names = TRUE)

selected_cells <- c(pilot_cells, normal_cells)
missing_cells <- setdiff(selected_cells, layer_membership$cell_id)
if (length(missing_cells) > 0) stop("Selected cells not in count layers.")

# --- Load hg38 feature mapping ---
infercna::useGenome("hg38")
genome <- infercna::retrieveGenome(name = "hg38")

mapping_file <- file.path(preflight_dir, "STEP17G1_GSE221561_feature_hg38_mapping.csv")
if (!file.exists(mapping_file)) stop("Missing: ", mapping_file)

gene_map <- fread(mapping_file)
mapped <- gene_map[mapping_type != "UNMAPPED" & !is.na(hg38_symbol)]

# --- Extract raw counts sample-wise ---
extract_counts <- function(cells, layer_membership, rna) {
  cell_layers <- layer_membership[cell_id %in% cells]
  parts <- lapply(unique(cell_layers$count_layer), function(lyr) {
    lyr_cells <- cell_layers[count_layer == lyr, cell_id]
    SeuratObject::LayerData(rna, layer = lyr, cells = lyr_cells, fast = FALSE)
  })
  if (length(parts) == 1) return(parts[[1]])
  common_genes <- Reduce(intersect, lapply(parts, rownames))
  do.call(cbind, lapply(parts, function(p) p[common_genes, , drop = FALSE]))
}

cat("Extracting pilot cells...\n")
pilot_raw <- extract_counts(pilot_cells, layer_membership, rna)
cat("Extracting normal cells...\n")
normal_raw <- extract_counts(normal_cells, layer_membership, rna)

# --- Map to hg38 ---
mapped_lookup <- mapped[, .(input_feature, hg38_symbol)]
mapped_lookup <- mapped_lookup[!duplicated(mapped_lookup)]

map_and_deduplicate <- function(raw_mat, label) {
  genes <- rownames(raw_mat)
  mapped_idx <- match(genes, mapped_lookup$input_feature)
  mapped_syms <- mapped_lookup$hg38_symbol[mapped_idx]
  has_symbol <- !is.na(mapped_syms)
  mapped_mat <- raw_mat[has_symbol, , drop = FALSE]
  rownames(mapped_mat) <- mapped_syms[has_symbol]
  dup_rows <- which(duplicated(rownames(mapped_mat)))
  if (length(dup_rows) > 0) {
    dup_names <- rownames(mapped_mat)[dup_rows]
    for (dn in unique(dup_names)) {
      idx <- which(rownames(mapped_mat) == dn)
      mapped_mat[idx[1], ] <- Matrix::colSums(mapped_mat[idx, , drop = FALSE])
      mapped_mat <- mapped_mat[-idx[-1], , drop = FALSE]
    }
  }
  cat(sprintf("  %s: %d genes x %d cells\n", label, nrow(mapped_mat), ncol(mapped_mat)))
  mapped_mat
}

pilot_mapped <- map_and_deduplicate(pilot_raw, "pilot")
normal_mapped <- map_and_deduplicate(normal_raw, "normal")

# --- CPM normalization on combined matrix ---
common_genes_raw <- intersect(rownames(pilot_mapped), rownames(normal_mapped))
combined_raw <- cbind(
  pilot_mapped[common_genes_raw, , drop = FALSE],
  normal_mapped[common_genes_raw, , drop = FALSE]
)
cpm <- t(t(combined_raw) / Matrix::colSums(combined_raw) * 1e6)

cpm_sums <- Matrix::colSums(cpm)
if (max(abs(cpm_sums - 1e6)) > 1) stop("CPM validation failed.")

cat("Combined CPM: ", nrow(cpm), " genes x ", ncol(cpm), " cells\n", sep = "")

# --- Build reference list ---
refCells <- lapply(normal_samples, function(s) normal_meta[sample_id == s, cell_id])
names(refCells) <- normal_samples

# --- Annotation ---
annotation <- md[
  cell_id %in% colnames(cpm),
  .(cell_id, sample_id, patient_id, tissue = tissue_simple, treatment_group, cell_class_std)
]
annotation[, infercna_role := fifelse(cell_id %in% pilot_cells, "TUMOR_TEST", "NORMAL_REFERENCE")]
annotation[, reference_group := fifelse(infercna_role == "NORMAL_REFERENCE", sample_id, NA_character_)]

# --- Gene table ---
gene_table <- data.table(gene = rownames(cpm))
genome_sub <- genome[match(gene_table$gene, genome$symbol), ]
gene_table[, chromosome := as.character(genome_sub$chromosome_name)]
gene_table[, arm := as.character(genome_sub$arm)]
gene_table[, start_position := genome_sub$start_position]
gene_table[, end_position := genome_sub$end_position]

# --- Save prep outputs ---
saveRDS(cpm, file.path(pilot_dir, paste0("GSE221561_", sample_id, "_infercna_CPM_hg38.rds")), compress = FALSE)
saveRDS(refCells, file.path(pilot_dir, paste0("GSE221561_", sample_id, "_infercna_refCells.rds")))
fwrite(annotation, file.path(pilot_dir, paste0("GSE221561_", sample_id, "_infercna_cell_annotation.csv")), bom = TRUE)
fwrite(gene_table, file.path(pilot_dir, paste0("GSE221561_", sample_id, "_infercna_gene_table.csv")), bom = TRUE)

cat("\nPhase A complete.\n\n")

# Cleanup
rm(obj, rna, pilot_raw, normal_raw, pilot_mapped, normal_mapped, cpm_pilot, cpm_normal)
gc(verbose = FALSE)

# ============================================================
# B. BALANCE REFERENCES (from 17H7B template)
# ============================================================

cat("===== PHASE B: BALANCE REFERENCES =====\n\n")

n_per_ref <- 500L

set.seed(22156117)
s6417n_cells <- refCells[["S6417N"]]
s6417n_balanced <- if (length(s6417n_cells) > n_per_ref) sample(s6417n_cells, n_per_ref) else s6417n_cells

set.seed(22156117)
s9478n_cells <- refCells[["S9478N"]]
s9478n_balanced <- if (length(s9478n_cells) > n_per_ref) sample(s9478n_cells, n_per_ref) else s9478n_cells

refCells_balanced <- list(S6417N = s6417n_balanced, S9478N = s9478n_balanced)

cat("S6417N balanced: ", length(s6417n_balanced), "\n", sep = "")
cat("S9478N balanced: ", length(s9478n_balanced), "\n", sep = "")

all_balanced_cells <- c(pilot_cells, unlist(refCells_balanced, use.names = FALSE))
cpm_balanced <- cpm[, all_balanced_cells, drop = FALSE]

# Remove zero-variance genes
gene_var <- apply(as.matrix(cpm_balanced), 1, var, na.rm = TRUE)
zero_genes <- names(gene_var[gene_var == 0 | is.na(gene_var)])
if (length(zero_genes) > 0) {
  cpm_balanced <- cpm_balanced[setdiff(rownames(cpm_balanced), zero_genes), , drop = FALSE]
  cat("Removed ", length(zero_genes), " zero-variance genes\n", sep = "")
}

# Validate hg38
hg38_symbols <- as.character(genome$symbol)
missing_hg38 <- setdiff(rownames(cpm_balanced), hg38_symbols)
if (length(missing_hg38) > 0) stop("Symbols absent from hg38.")

# Rebuild gene table for balanced
gene_table_balanced <- data.table(gene = rownames(cpm_balanced))
genome_sub2 <- genome[match(gene_table_balanced$gene, genome$symbol), ]
gene_table_balanced[, chromosome := as.character(genome_sub2$chromosome_name)]
gene_table_balanced[, arm := as.character(genome_sub2$arm)]
gene_table_balanced[, start_position := genome_sub2$start_position]
gene_table_balanced[, end_position := genome_sub2$end_position]

# Save balanced inputs
saveRDS(cpm_balanced, file.path(pilot_dir, paste0("GSE221561_", sample_id, "_BALANCED_infercna_CPM_hg38.rds")), compress = FALSE)
saveRDS(refCells_balanced, file.path(pilot_dir, paste0("GSE221561_", sample_id, "_BALANCED_infercna_refCells.rds")))
fwrite(annotation, file.path(pilot_dir, paste0("GSE221561_", sample_id, "_BALANCED_infercna_cell_annotation.csv")), bom = TRUE)
fwrite(gene_table_balanced, file.path(pilot_dir, paste0("GSE221561_", sample_id, "_BALANCED_infercna_gene_table.csv")), bom = TRUE)

cat("Balanced CPM: ", nrow(cpm_balanced), " genes x ", ncol(cpm_balanced), " cells\n", sep = "")

rm(cpm, gene_var)
gc(verbose = FALSE)

cat("\nPhase B complete.\n\n")

# ============================================================
# C. INFERCNA (from 17H7C template)
# ============================================================

cat("===== PHASE C: RUN INFERCNA =====\n\n")

# Locate balanced files
cpm_file <- file.path(pilot_dir, paste0("GSE221561_", sample_id, "_BALANCED_infercna_CPM_hg38.rds"))
refs_file <- file.path(pilot_dir, paste0("GSE221561_", sample_id, "_BALANCED_infercna_refCells.rds"))
annotation_file <- file.path(pilot_dir, paste0("GSE221561_", sample_id, "_BALANCED_infercna_cell_annotation.csv"))
gene_file <- file.path(pilot_dir, paste0("GSE221561_", sample_id, "_BALANCED_infercna_gene_table.csv"))

cpm <- readRDS(cpm_file)
refCells <- readRDS(refs_file)
ann <- fread(annotation_file)
gene_table <- fread(gene_file)

tumor_cells <- ann[infercna_role == "TUMOR_TEST", cell_id]
normal_cells <- unlist(refCells, use.names = FALSE)

cat("Tumor cells  : ", length(tumor_cells), "\n", sep = "")
cat("Normal cells : ", length(normal_cells), "\n", sep = "")

# Validation
if (ncol(cpm) != length(tumor_cells) + length(normal_cells)) stop("Cell count mismatch.")
if (!all(tumor_cells %in% colnames(cpm))) stop("Tumor cells missing.")
if (!all(normal_cells %in% colnames(cpm))) stop("Normal cells missing.")
if (length(intersect(tumor_cells, normal_cells)) > 0) stop("Overlap detected.")

# Dense matrix
m <- as.matrix(cpm)
storage.mode(m) <- "double"
if (any(!is.finite(m))) stop("Non-finite CPM values.")
rm(cpm); gc(verbose = FALSE)

# Frozen parameters
params <- data.table(
  parameter = c("genome", "window", "range", "n", "noise", "center.method", "isLog"),
  value = c("hg38", "100", "-3,3", "5000", "0.1", "median", "FALSE")
)
fwrite(params, file.path(results_dir, paste0("GSE221561_", sample_id, "_17H9_parameters.csv")), bom = TRUE)

# Run InferCNA
cat("\nRunning InferCNA...\n")
start_time <- Sys.time()

cna <- tryCatch(
  infercna::infercna(
    m = m, refCells = refCells,
    window = 100, range = c(-3, 3), n = 5000,
    noise = 0.1, center.method = "median", isLog = FALSE,
    verbose = TRUE
  ),
  error = function(e) {
    cat("INFERCNA ERROR: ", conditionMessage(e), "\n")
    stop("InferCNA failed for ", sample_id)
  }
)

runtime_min <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
cat("Runtime: ", sprintf("%.2f min", runtime_min), "\n", sep = "")

# Validate CNA
if (any(!is.finite(cna))) stop("Non-finite CNA values.")
if (!all(tumor_cells %in% colnames(cna))) stop("Tumor cells absent from CNA.")
if (!all(normal_cells %in% colnames(cna))) stop("Normal cells absent from CNA.")

# Save CNA
saveRDS(cna, file.path(results_dir, paste0("GSE221561_", sample_id, "_infercna_CNA_all_cells.rds")), compress = FALSE)
cna_tumor <- cna[, tumor_cells, drop = FALSE]
saveRDS(cna_tumor, file.path(results_dir, paste0("GSE221561_", sample_id, "_infercna_CNA_tumor_only.rds")), compress = FALSE)

# ============================================================
# D. METRICS
# ============================================================

cat("\n===== PHASE D: METRICS =====\n\n")

signal <- infercna::cnaSignal(cna)
if (is.null(names(signal))) names(signal) <- colnames(cna)

tumor_cor <- infercna::cnaCor(cna_tumor)
if (is.null(names(tumor_cor))) names(tumor_cor) <- colnames(cna_tumor)

tumor_mean_profile <- rowMeans(cna_tumor)
normal_cor <- sapply(normal_cells, function(cell) {
  suppressWarnings(cor(cna[, cell], tumor_mean_profile, method = "pearson", use = "pairwise.complete.obs"))
})

tumor_metrics <- data.table(
  cell_id = tumor_cells, role = "TUMOR_TEST",
  cna_signal = as.numeric(signal[tumor_cells]),
  cna_cor_to_tumor = as.numeric(tumor_cor[tumor_cells]),
  mean_abs_CNA = colMeans(abs(cna_tumor))
)
normal_metrics <- data.table(
  cell_id = normal_cells, role = "NORMAL_REFERENCE",
  cna_signal = as.numeric(signal[normal_cells]),
  cna_cor_to_tumor = as.numeric(normal_cor[normal_cells]),
  mean_abs_CNA = colMeans(abs(cna[, normal_cells, drop = FALSE]))
)
cell_metrics <- rbindlist(list(tumor_metrics, normal_metrics))
cell_metrics <- merge(cell_metrics, ann[, .(cell_id, sample = sample_var)], by = "cell_id", all.x = TRUE)

# Summary
normal_signal_95 <- as.numeric(quantile(normal_metrics$cna_signal, 0.95, na.rm = TRUE))
normal_signal_99 <- as.numeric(quantile(normal_metrics$cna_signal, 0.99, na.rm = TRUE))

tumor_signal_med <- median(tumor_metrics$cna_signal, na.rm = TRUE)
normal_signal_med <- median(normal_metrics$cna_signal, na.rm = TRUE)
signal_delta <- tumor_signal_med - normal_signal_med
signal_ratio <- tumor_signal_med / normal_signal_med
tumor_cor_med <- median(tumor_metrics$cna_cor_to_tumor, na.rm = TRUE)
normal_cor_med <- median(normal_metrics$cna_cor_to_tumor, na.rm = TRUE)
corr_separation <- tumor_cor_med - normal_cor_med
tumor_above_95 <- mean(tumor_metrics$cna_signal > normal_signal_95, na.rm = TRUE)
tumor_above_99 <- mean(tumor_metrics$cna_signal > normal_signal_99, na.rm = TRUE)

summary_dt <- data.table(
  metric = c("Tumor_cells", "Normal_cells", "Tumor_median_CNA_signal", "Normal_median_CNA_signal",
             "CNA_signal_delta", "CNA_signal_ratio", "Tumor_median_cnaCor",
             "Normal_median_corr_to_tumor", "Correlation_separation",
             "Tumor_above_normal_95th", "Tumor_above_normal_99th"),
  value = c(length(tumor_cells), length(normal_cells), tumor_signal_med, normal_signal_med,
            signal_delta, signal_ratio, tumor_cor_med, normal_cor_med, corr_separation,
            tumor_above_95, tumor_above_99)
)
fwrite(summary_dt, file.path(results_dir, paste0("GSE221561_", sample_id, "_17H9_tumor_normal_summary.csv")), bom = TRUE)
fwrite(cell_metrics, file.path(results_dir, paste0("GSE221561_", sample_id, "_17H9_cell_CNA_metrics.csv")), bom = TRUE)

# ============================================================
# E. ARM AUDIT
# ============================================================

cat("\n===== PHASE E: CHROMOSOME-ARM AUDIT =====\n\n")

g <- gene_table[gene %in% rownames(cna)]
g[, row_index := match(gene, rownames(cna))]
g <- g[!is.na(row_index) & !is.na(chromosome) & !is.na(arm)]
g[, chromosome_arm := paste0(chromosome, arm)]

arm_index <- split(g$row_index, g$chromosome_arm)
arm_list <- lapply(names(arm_index), function(a) {
  idx <- unique(arm_index[[a]])
  if (length(idx) < 20) return(NULL)
  ta <- colMeans(cna[idx, tumor_cells, drop = FALSE])
  na_vals <- colMeans(cna[idx, normal_cells, drop = FALSE])
  data.table(
    chromosome_arm = a, genes = length(idx),
    tumor_median = median(ta, na.rm = TRUE),
    normal_median = median(na_vals, na.rm = TRUE),
    delta = median(ta, na.rm = TRUE) - median(na_vals, na.rm = TRUE),
    tumor_gain_fraction = mean(ta > 0.05, na.rm = TRUE),
    tumor_loss_fraction = mean(ta < -0.05, na.rm = TRUE),
    normal_gain_fraction = mean(na_vals > 0.05, na.rm = TRUE),
    normal_loss_fraction = mean(na_vals < -0.05, na.rm = TRUE)
  )
})
arm_dt <- rbindlist(arm_list, fill = TRUE)
arm_dt[, abs_delta := abs(delta)]
arm_dt[, coherent_audit_flag := abs_delta >= 0.05 & (tumor_gain_fraction >= 0.50 | tumor_loss_fraction >= 0.50)]
setorder(arm_dt, -abs_delta)
coherent <- arm_dt[coherent_audit_flag == TRUE]

fwrite(arm_dt, file.path(results_dir, paste0("GSE221561_", sample_id, "_17H9_all_arm_metrics.csv")), bom = TRUE)
fwrite(coherent, file.path(results_dir, paste0("GSE221561_", sample_id, "_17H9_candidate_coherent_arms.csv")), bom = TRUE)

cat("Coherent arms: ", nrow(coherent), "\n", sep = "")
if (nrow(coherent) > 0) {
  cat("Arm names: ", paste(sort(coherent$chromosome_arm), collapse = ", "), "\n", sep = "")
}

# ============================================================
# F. FIGURES
# ============================================================

cat("\n===== PHASE F: FIGURES =====\n\n")

genome_dt <- as.data.table(genome)
heat_genes <- data.table(gene = rownames(cna))
heat_genes[, genome_idx := match(gene, as.character(genome_dt$symbol))]
heat_genes[, chromosome := as.character(genome_dt$chromosome_name[genome_idx])]
heat_genes[, start := genome_dt$start_position[genome_idx]]
heat_genes[, chr_order := fifelse(chromosome == "X", 23L, fifelse(chromosome == "Y", 24L, suppressWarnings(as.integer(chromosome))))]
heat_genes <- heat_genes[!is.na(chr_order) & !is.na(start)]
setorder(heat_genes, chr_order, start)
ordered_genes <- heat_genes$gene

cna_heat <- cna[ordered_genes, , drop = FALSE]
normal_order <- unlist(refCells, use.names = FALSE)
cell_order <- c(normal_order, tumor_cells)
cell_order <- cell_order[cell_order %in% colnames(cna_heat)]
cna_heat <- cna_heat[, cell_order, drop = FALSE]
chr_bounds <- heat_genes[, .(start_index = .I[1], end_index = .I[.N], midpoint = mean(c(.I[1], .I[.N]))), by = chromosome]

plot_limit <- 0.6

# Heatmap 1: all cells
heat_png <- file.path(figdir, paste0("GSE221561_", sample_id, "_InferCNA_manual_all_cells_heatmap.png"))
png(heat_png, width = 4200, height = 2800, res = 300)
par(mar = c(5, 6, 4, 2))
z <- t(cna_heat)
z <- pmax(pmin(z, plot_limit), -plot_limit)
image(x = seq_len(ncol(z)), y = seq_len(nrow(z)), z = t(z), axes = FALSE,
      xlab = "Chromosome", ylab = "Cells",
      main = paste0("GSE221561 ", sample_id, " InferCNA"),
      col = hcl.colors(101, "Blue-Red 3"))
for (b in chr_bounds$end_index) abline(v = b + 0.5, lwd = 0.4)
axis(1, at = chr_bounds$midpoint, labels = chr_bounds$chromosome, las = 2, cex.axis = 0.65)
cell_group <- rep(NA_character_, length(cell_order))
cell_group[cell_order %in% refCells[["S6417N"]]] <- "S6417N"
cell_group[cell_order %in% refCells[["S9478N"]]] <- "S9478N"
cell_group[cell_order %in% tumor_cells] <- sample_id
rr <- rle(cell_group)
group_end <- cumsum(rr$lengths)
group_start <- c(1, head(group_end, -1) + 1)
group_mid <- (group_start + group_end) / 2
for (b in head(group_end, -1)) abline(h = b + 0.5, lwd = 1.2)
axis(2, at = group_mid, labels = rr$values, las = 2, cex.axis = 0.8)
dev.off()

# Heatmap 2: tumor only
tumor_cor_order <- names(sort(tumor_cor, decreasing = TRUE))
tumor_cor_order <- tumor_cor_order[tumor_cor_order %in% colnames(cna_heat)]
tumor_heat <- cna_heat[, tumor_cor_order, drop = FALSE]

tumor_png <- file.path(figdir, paste0("GSE221561_", sample_id, "_InferCNA_manual_tumor_only_heatmap.png"))
png(tumor_png, width = 4200, height = 2000, res = 300)
par(mar = c(5, 5, 4, 2))
zt <- t(tumor_heat)
zt <- pmax(pmin(zt, plot_limit), -plot_limit)
image(x = seq_len(ncol(zt)), y = seq_len(nrow(zt)), z = t(zt), axes = FALSE,
      xlab = "Chromosome", ylab = paste0(sample_id, " Tumor Epithelial Cells"),
      main = paste0(sample_id, " Tumor Epithelial CNA Profiles"),
      col = hcl.colors(101, "Blue-Red 3"))
for (b in chr_bounds$end_index) abline(v = b + 0.5, lwd = 0.4)
axis(1, at = chr_bounds$midpoint, labels = chr_bounds$chromosome, las = 2, cex.axis = 0.65)
dev.off()

# Scatter plot
scatter_png <- file.path(figdir, paste0("GSE221561_", sample_id, "_CNA_signal_vs_cnaCor.png"))
png(scatter_png, width = 3200, height = 2800, res = 300)
par(mar = c(5, 5, 4, 2))
plot(normal_metrics$cna_signal, normal_metrics$cna_cor_to_tumor,
     col = adjustcolor("grey50", alpha.f = 0.5), pch = 16, cex = 0.8,
     xlab = "CNA Signal", ylab = "cnaCor",
     main = paste0(sample_id, ": CNA Signal vs cnaCor"),
     xlim = range(c(normal_metrics$cna_signal, tumor_metrics$cna_signal), na.rm = TRUE),
     ylim = range(c(normal_metrics$cna_cor_to_tumor, tumor_metrics$cna_cor_to_tumor), na.rm = TRUE))
points(tumor_metrics$cna_signal, tumor_metrics$cna_cor_to_tumor,
       col = adjustcolor("red", alpha.f = 0.6), pch = 16, cex = 0.8)
abline(h = 0, lty = 2, col = "grey60")
abline(v = normal_signal_95, lty = 2, col = "blue")
legend("topleft", legend = c("Normal references", paste0("Tumor (", sample_id, ")")),
       col = c("grey50", "red"), pch = 16, bty = "n", cex = 0.9)
dev.off()

cat("Figures saved.\n")

# ============================================================
# G. CANONICAL MANIFEST
# ============================================================

cat("\n===== PHASE G: MANIFEST =====\n\n")

# --- Per-sample validation ---
# 1. tumor cell count
cm_tumor_count <- sum(ann$infercna_role == "TUMOR_TEST")
cna_tumor_ncol <- ncol(cna_tumor)
valid_tumor_count <- (cm_tumor_count == cna_tumor_ncol) && (cna_tumor_ncol == length(tumor_cells))

# 2. gene count
valid_gene_count <- (nrow(cna) == nrow(cna_tumor))

# 3. arm count
valid_arm_count <- (nrow(coherent) == length(unique(coherent$chromosome_arm)))

# 4. balanced references
valid_balanced <- (length(refCells[["S6417N"]]) == 500) && (length(refCells[["S9478N"]]) == 500)

# 5. CPM column sums
valid_cpm <- max(abs(Matrix::colSums(readRDS(cpm_file)) - 1e6)) <= 1

# 6. no overlap
valid_no_overlap <- length(intersect(tumor_cells, normal_cells)) == 0

# 7. all CNA finite
valid_finite <- all(is.finite(cna))

all_valid <- all(valid_tumor_count, valid_gene_count, valid_arm_count,
                 valid_balanced, valid_cpm, valid_no_overlap, valid_finite)

validation_status <- if (all_valid) "PASS" else "FAIL"

cat("Validation: ", validation_status, "\n", sep = "")
if (!all_valid) {
  cat("  tumor_count: ", valid_tumor_count, "\n", sep = "")
  cat("  gene_count: ", valid_gene_count, "\n", sep = "")
  cat("  arm_count: ", valid_arm_count, "\n", sep = "")
  cat("  balanced: ", valid_balanced, "\n", sep = "")
  cat("  cpm: ", valid_cpm, "\n", sep = "")
  cat("  no_overlap: ", valid_no_overlap, "\n", sep = "")
  cat("  finite: ", valid_finite, "\n", sep = "")
  stop("VALIDATION FAILED for ", sample_id)
}

arm_names_str <- if (nrow(coherent) > 0) paste(sort(coherent$chromosome_arm), collapse = ", ") else "none"

manifest <- data.table(
  dataset = "GSE221561",
  sample = sample_id,
  treatment = treatment,
  tumor_cells = length(tumor_cells),
  normal_reference_cells = length(normal_cells),
  CNA_genes = nrow(cna),
  tumor_median_CNA_signal = tumor_signal_med,
  normal_median_CNA_signal = normal_signal_med,
  CNA_signal_delta = signal_delta,
  CNA_signal_ratio = signal_ratio,
  tumor_median_cnaCor = tumor_cor_med,
  normal_median_corr_to_tumor = normal_cor_med,
  correlation_separation = corr_separation,
  tumor_above_normal95 = tumor_above_95,
  tumor_above_normal99 = tumor_above_99,
  coherent_CNA_arms = nrow(coherent),
  coherent_CNA_arm_names = arm_names_str,
  runtime_minutes = runtime_min,
  findMalignant_run = FALSE,
  malignant_labels_assigned = FALSE,
  validation_status = validation_status,
  canonical_manifest_path = file.path(results_dir, paste0("GSE221561_", sample_id, "_17H9_manifest.csv"))
)

fwrite(manifest, file.path(results_dir, paste0("GSE221561_", sample_id, "_17H9_manifest.csv")), bom = TRUE)

# ============================================================
# H. FINAL SUMMARY
# ============================================================

cat("\n===== PHASE H: FINAL SUMMARY =====\n\n")

cat("============================================\n")
cat("17H9 SAMPLE COMPLETE\n")
cat("Sample: ", sample_id, "\n", sep = "")
cat("Treatment: ", treatment, "\n", sep = "")
cat("Tumor cells: ", length(tumor_cells), "\n", sep = "")
cat("CNA genes: ", nrow(cna), "\n", sep = "")
cat("Signal ratio: ", sprintf("%.2f", signal_ratio), "\n", sep = "")
cat("Tumor cnaCor: ", sprintf("%.3f", tumor_cor_med), "\n", sep = "")
cat("Correlation separation: ", sprintf("%.3f", corr_separation), "\n", sep = "")
cat("Tumor > normal 95th: ", sprintf("%.1f%%", 100 * tumor_above_95), "\n", sep = "")
cat("Tumor > normal 99th: ", sprintf("%.1f%%", 100 * tumor_above_99), "\n", sep = "")
cat("Coherent arms: ", nrow(coherent), "\n", sep = "")
if (nrow(coherent) > 0) cat("Coherent arm names: ", arm_names_str, "\n", sep = "")
cat("Validation: ", validation_status, "\n", sep = "")
cat("findMalignant: NO\n")
cat("Malignant labels: NO\n")
cat("============================================\n")
