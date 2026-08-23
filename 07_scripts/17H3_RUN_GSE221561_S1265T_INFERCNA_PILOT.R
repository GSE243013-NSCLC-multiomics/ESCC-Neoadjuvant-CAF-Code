# ============================================================
# ESCC Neoadjuvant Project
# STEP 17H3
#
# Run InferCNA on BALANCED GSE221561 S1265T pilot
#
# Tumor:  S1265T epithelial cells = 4972
# Normal: S6417N = 500, S9478N = 500
#
# Frozen parameters — identical to GSE197677 and S6417T:
#   genome        = hg38
#   window        = 100
#   range         = c(-3, 3)
#   n             = 5000
#   noise         = 0.1
#   center.method = median
#   isLog         = FALSE
#
# IMPORTANT:
#   DOES run InferCNA
#   DOES NOT run findMalignant()
#   DOES NOT assign malignant labels
#   DOES NOT filter tumor cells
#   DOES NOT process other samples
# ============================================================

options(
  stringsAsFactors = FALSE,
  timeout = 3600
)

suppressPackageStartupMessages({
  library(infercna)
  library(Matrix)
  library(data.table)
})

cat("\n============================================================\n")
cat("STEP 17H3: GSE221561 S1265T INFERCNA PILOT\n")
cat("============================================================\n\n")

# ============================================================
# 1. PATHS
# ============================================================

indir <- "06_tables/cross_dataset/malignant_epithelial/GSE221561_S1265T_pilot"

outdir <- file.path(indir, "17H3_results")

figdir <- "05_figures/GSE221561/malignant_epithelial/S1265T_pilot"

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

if (!dir.exists(indir)) {
  stop("Pilot directory missing: ", indir)
}

# ============================================================
# 2. LOCATE BALANCED INPUT FILES
# ============================================================

cat("===== 1. LOCATE BALANCED 17H2 INPUTS =====\n\n")

all_files <- list.files(indir, full.names = TRUE)

pick_balanced_file <- function(include, extension, description) {
  x <- all_files[
    grepl(include, basename(all_files), ignore.case = TRUE) &
    grepl(extension, basename(all_files), ignore.case = TRUE)
  ]
  xb <- x[grepl("BALANCED", basename(x), ignore.case = TRUE)]
  if (length(xb) > 0) x <- xb
  if (length(x) == 0) stop("Cannot find ", description, " in ", indir)
  if (length(x) > 1) {
    cat("\nMultiple candidates for ", description, ":\n", sep = "")
    print(x)
    x <- x[which.max(file.info(x)$mtime)]
  }
  x[1]
}

cpm_file <- pick_balanced_file("CPM", "\\.rds$", "balanced CPM matrix")
refs_file <- pick_balanced_file("refCells", "\\.rds$", "balanced reference list")
annotation_file <- pick_balanced_file("cell_annotation", "\\.csv$", "balanced cell annotation")
gene_file <- pick_balanced_file("gene_table", "\\.csv$", "balanced gene table")

cat("CPM        : ", cpm_file, "\n", sep = "")
cat("refCells   : ", refs_file, "\n", sep = "")
cat("Annotation : ", annotation_file, "\n", sep = "")
cat("Gene table : ", gene_file, "\n", sep = "")

# ============================================================
# 3. LOAD
# ============================================================

cat("\n===== 2. LOAD BALANCED PILOT =====\n\n")

cpm <- readRDS(cpm_file)
refCells <- readRDS(refs_file)
ann <- fread(annotation_file)
gene_table <- fread(gene_file)

cat("CPM matrix: ", nrow(cpm), " genes x ", ncol(cpm), " cells\n", sep = "")
cat("Matrix class: ", paste(class(cpm), collapse = ", "), "\n", sep = "")
cat("Reference groups:\n")
print(lengths(refCells))

# ============================================================
# 4. VALIDATION
# ============================================================

cat("\n===== 3. CELL / REFERENCE VALIDATION =====\n\n")

tumor_cells <- ann[infercna_role == "TUMOR_TEST", cell_id]
normal_cells <- unlist(refCells, use.names = FALSE)

cat("Tumor cells  : ", length(tumor_cells), "\n", sep = "")
cat("Normal cells : ", length(normal_cells), "\n", sep = "")
cat("Total        : ", length(tumor_cells) + length(normal_cells), "\n", sep = "")

if (ncol(cpm) != length(tumor_cells) + length(normal_cells)) {
  stop("CPM cell count mismatch.")
}
if (!all(tumor_cells %in% colnames(cpm))) stop("Tumor cells missing from CPM.")
if (!all(normal_cells %in% colnames(cpm))) stop("Normal cells missing from CPM.")
if (length(intersect(tumor_cells, normal_cells)) > 0) stop("Overlap detected.")
if (anyDuplicated(rownames(cpm)) > 0) stop("Duplicate gene symbols.")

cat("Validated\n")

# ============================================================
# 5. CPM VALIDATION
# ============================================================

cat("\n===== 4. CPM VALIDATION =====\n\n")

cpm_totals <- Matrix::colSums(cpm)
cat(sprintf("CPM column sums: median %.2f; range %.2f - %.2f\n", median(cpm_totals), min(cpm_totals), max(cpm_totals)))

if (max(abs(cpm_totals - 1e6)) > 1) stop("CPM validation failed.")
cat("CPM validated\n")

# ============================================================
# 6. hg38
# ============================================================

cat("\n===== 5. CONFIGURE hg38 =====\n\n")

infercna::useGenome("hg38")
genome <- infercna::retrieveGenome(name = "hg38")
hg38_symbols <- as.character(genome$symbol)

missing_hg38 <- setdiff(rownames(cpm), hg38_symbols)
cat("Input genes            : ", nrow(cpm), "\n", sep = "")
cat("Genes absent from hg38 : ", length(missing_hg38), "\n", sep = "")

if (length(missing_hg38) > 0) {
  stop("Balanced CPM contains symbols absent from InferCNA hg38.")
}

cat("All genes compatible with InferCNA hg38\n")

# ============================================================
# 7. DENSE MATRIX
# ============================================================

cat("\n===== 6. PREPARE DENSE MATRIX =====\n\n")

dense_gb <- (nrow(cpm) * ncol(cpm) * 8) / 1024^3
cat(sprintf("Expected dense size: %.3f GB\n", dense_gb))

m <- as.matrix(cpm)
storage.mode(m) <- "double"

cat("Dense matrix: ", nrow(m), " genes x ", ncol(m), " cells\n", sep = "")

if (any(!is.finite(m))) stop("Dense CPM contains NA/Inf.")

rm(cpm)
gc(verbose = FALSE)

# ============================================================
# 8. FROZEN PARAMETERS
# ============================================================

cat("\n============================================================\n")
cat("7. FROZEN INFERCNA PARAMETERS\n")
cat("============================================================\n\n")

params <- data.table(
  parameter = c("genome", "window", "range", "n", "noise", "center.method", "isLog"),
  value = c("hg38", "100", "-3,3", "5000", "0.1", "median", "FALSE")
)
print(params)
fwrite(params, file.path(outdir, "GSE221561_S1265T_17H3_parameters.csv"), bom = TRUE)

# ============================================================
# 9. RUN INFERCNA
# ============================================================

cat("\n============================================================\n")
cat("8. RUNNING INFERCNA — S1265T PILOT\n")
cat("============================================================\n\n")

cat("This is the second GSE221561 CNV run.\n")
cat("No malignant calls will be made in this step.\n\n")

start_time <- Sys.time()

cna <- tryCatch(
  infercna::infercna(
    m = m,
    refCells = refCells,
    window = 100,
    range = c(-3, 3),
    n = 5000,
    noise = 0.1,
    center.method = "median",
    isLog = FALSE,
    verbose = TRUE
  ),
  error = function(e) {
    cat("\n============================================\n")
    cat("INFERCNA ERROR\n")
    cat("============================================\n")
    cat(conditionMessage(e), "\n")
    stop("Step 17H3 stopped during InferCNA.")
  }
)

end_time <- Sys.time()
runtime_min <- as.numeric(difftime(end_time, start_time, units = "mins"))

cat("\nInferCNA runtime: ", sprintf("%.2f min", runtime_min), "\n", sep = "")

# ============================================================
# 10. CNA MATRIX VALIDATION
# ============================================================

cat("\n============================================================\n")
cat("9. CNA MATRIX VALIDATION\n")
cat("============================================================\n\n")

cat("CNA matrix: ", nrow(cna), " genes x ", ncol(cna), " cells\n", sep = "")
cat(sprintf("CNA range : %.4f to %.4f\n", min(cna, na.rm = TRUE), max(cna, na.rm = TRUE)))
cat("NA values  : ", sum(is.na(cna)), "\n", sep = "")
cat("Inf values : ", sum(is.infinite(cna)), "\n", sep = "")

if (any(!is.finite(cna))) stop("CNA matrix contains non-finite values.")
if (!all(tumor_cells %in% colnames(cna))) stop("Tumor cells absent from CNA output.")
if (!all(normal_cells %in% colnames(cna))) stop("Normal cells absent from CNA output.")

cat("CNA matrix validated\n")

# ============================================================
# 11. SAVE CNA
# ============================================================

cat("\n===== 10. SAVE CNA OUTPUT =====\n\n")

all_cna_file <- file.path(outdir, "GSE221561_S1265T_infercna_CNA_all_cells.rds")
tumor_cna_file <- file.path(outdir, "GSE221561_S1265T_infercna_CNA_tumor_only.rds")

saveRDS(cna, all_cna_file, compress = FALSE)
cna_tumor <- cna[, tumor_cells, drop = FALSE]
saveRDS(cna_tumor, tumor_cna_file, compress = FALSE)

cat("Saved: ", all_cna_file, "\n", sep = "")
cat("Saved: ", tumor_cna_file, "\n", sep = "")

# ============================================================
# 12. CNA SIGNAL
# ============================================================

cat("\n============================================================\n")
cat("11. CNA SIGNAL\n")
cat("============================================================\n\n")

signal <- infercna::cnaSignal(cna)
if (is.null(names(signal))) names(signal) <- colnames(cna)

# ============================================================
# 13. TUMOR cnaCor
# ============================================================

cat("\n============================================================\n")
cat("12. TUMOR CNA CORRELATION\n")
cat("============================================================\n\n")

tumor_cor <- infercna::cnaCor(cna_tumor)
if (is.null(names(tumor_cor))) names(tumor_cor) <- colnames(cna_tumor)

tumor_mean_profile <- rowMeans(cna_tumor)

normal_cor <- sapply(normal_cells, function(cell) {
  suppressWarnings(cor(cna[, cell], tumor_mean_profile, method = "pearson", use = "pairwise.complete.obs"))
})

# ============================================================
# 14. CELL METRICS
# ============================================================

tumor_metrics <- data.table(
  cell_id = tumor_cells,
  role = "TUMOR_TEST",
  cna_signal = as.numeric(signal[tumor_cells]),
  cna_cor_to_tumor = as.numeric(tumor_cor[tumor_cells]),
  mean_abs_CNA = colMeans(abs(cna_tumor))
)

normal_metrics <- data.table(
  cell_id = normal_cells,
  role = "NORMAL_REFERENCE",
  cna_signal = as.numeric(signal[normal_cells]),
  cna_cor_to_tumor = as.numeric(normal_cor[normal_cells]),
  mean_abs_CNA = colMeans(abs(cna[, normal_cells, drop = FALSE]))
)

cell_metrics <- rbindlist(list(tumor_metrics, normal_metrics))
cell_metrics <- merge(cell_metrics, ann[, .(cell_id, sample = sample_id)], by = "cell_id", all.x = TRUE)

# ============================================================
# 15. TUMOR VS NORMAL SUMMARY
# ============================================================

cat("\n============================================================\n")
cat("13. S1265T TUMOR VS NORMAL SUMMARY\n")
cat("============================================================\n\n")

normal_signal_95 <- as.numeric(quantile(normal_metrics$cna_signal, 0.95, na.rm = TRUE))
normal_signal_99 <- as.numeric(quantile(normal_metrics$cna_signal, 0.99, na.rm = TRUE))
normal_cor_95 <- as.numeric(quantile(normal_metrics$cna_cor_to_tumor, 0.95, na.rm = TRUE))

summary_dt <- data.table(
  metric = c(
    "Cells",
    "Median CNA signal",
    "Mean CNA signal",
    "Median mean abs CNA",
    "Median CNA correlation to tumor",
    "Tumor cells above normal 95th CNA signal",
    "Tumor cells above normal 99th CNA signal"
  ),
  Tumor_S1265T = c(
    length(tumor_cells),
    median(tumor_metrics$cna_signal, na.rm = TRUE),
    mean(tumor_metrics$cna_signal, na.rm = TRUE),
    median(tumor_metrics$mean_abs_CNA, na.rm = TRUE),
    median(tumor_metrics$cna_cor_to_tumor, na.rm = TRUE),
    mean(tumor_metrics$cna_signal > normal_signal_95, na.rm = TRUE),
    mean(tumor_metrics$cna_signal > normal_signal_99, na.rm = TRUE)
  ),
  Normal_references = c(
    length(normal_cells),
    median(normal_metrics$cna_signal, na.rm = TRUE),
    mean(normal_metrics$cna_signal, na.rm = TRUE),
    median(normal_metrics$mean_abs_CNA, na.rm = TRUE),
    median(normal_metrics$cna_cor_to_tumor, na.rm = TRUE),
    0.05, 0.01
  )
)

print(summary_dt, nrows = 100)
fwrite(summary_dt, file.path(outdir, "GSE221561_S1265T_17H3_tumor_normal_summary.csv"), bom = TRUE)
fwrite(cell_metrics, file.path(outdir, "GSE221561_S1265T_17H3_cell_CNA_metrics.csv"), bom = TRUE)

# ============================================================
# 16. REFERENCE SAMPLE SUMMARY
# ============================================================

cat("\n============================================================\n")
cat("14. REFERENCE-SAMPLE QC\n")
cat("============================================================\n\n")

sample_summary <- cell_metrics[
  ,
  .(
    cells = .N,
    median_CNA_signal = median(cna_signal, na.rm = TRUE),
    median_correlation_to_tumor = median(cna_cor_to_tumor, na.rm = TRUE),
    median_mean_abs_CNA = median(mean_abs_CNA, na.rm = TRUE)
  ),
  by = .(role, sample)
]
setorder(sample_summary, role, sample)
print(sample_summary)
fwrite(sample_summary, file.path(outdir, "GSE221561_S1265T_17H3_sample_QC.csv"), bom = TRUE)

# ============================================================
# 17. CHROMOSOME-ARM COHERENCE
# ============================================================

cat("\n============================================================\n")
cat("15. CHROMOSOME-ARM TUMOR VS NORMAL AUDIT\n")
cat("============================================================\n\n")

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
    chromosome_arm = a,
    genes = length(idx),
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

cat("Top chromosome-arm differences:\n\n")
print(head(arm_dt, 20))
cat("\nCandidate coherent CNA arms: ", nrow(coherent), "\n", sep = "")

if (nrow(coherent) > 0) print(coherent, nrows = 100)

fwrite(arm_dt, file.path(outdir, "GSE221561_S1265T_17H3_all_arm_metrics.csv"), bom = TRUE)
fwrite(coherent, file.path(outdir, "GSE221561_S1265T_17H3_candidate_coherent_arms.csv"), bom = TRUE)

# ============================================================
# 18. MANUAL CHROMOSOME-ORDERED HEATMAP
# ============================================================

cat("\n============================================================\n")
cat("16. MANUAL CNA HEATMAP\n")
cat("============================================================\n\n")

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

heat_png <- file.path(figdir, "GSE221561_S1265T_InferCNA_manual_all_cells_heatmap.png")

png(heat_png, width = 4200, height = 2800, res = 300)
par(mar = c(5, 6, 4, 2))

z <- t(cna_heat)
plot_limit <- 0.6
z <- pmax(pmin(z, plot_limit), -plot_limit)

image(x = seq_len(ncol(z)), y = seq_len(nrow(z)), z = t(z), axes = FALSE, xlab = "Chromosome", ylab = "Cells", main = "GSE221561 S1265T InferCNA Pilot", col = hcl.colors(101, "Blue-Red 3"))

for (b in chr_bounds$end_index) abline(v = b + 0.5, lwd = 0.4)
axis(1, at = chr_bounds$midpoint, labels = chr_bounds$chromosome, las = 2, cex.axis = 0.65)

cell_group <- rep(NA_character_, length(cell_order))
cell_group[cell_order %in% refCells[["S6417N"]]] <- "S6417N"
cell_group[cell_order %in% refCells[["S9478N"]]] <- "S9478N"
cell_group[cell_order %in% tumor_cells] <- "S1265T"

rr <- rle(cell_group)
group_end <- cumsum(rr$lengths)
group_start <- c(1, head(group_end, -1) + 1)
group_mid <- (group_start + group_end) / 2

for (b in head(group_end, -1)) abline(h = b + 0.5, lwd = 1.2)
axis(2, at = group_mid, labels = rr$values, las = 2, cex.axis = 0.8)

dev.off()
cat("Saved: ", heat_png, "\n", sep = "")

# ============================================================
# 19. TUMOR-ONLY HEATMAP
# ============================================================

tumor_cor_order <- names(sort(tumor_cor, decreasing = TRUE))
tumor_cor_order <- tumor_cor_order[tumor_cor_order %in% colnames(cna_heat)]

tumor_heat <- cna_heat[, tumor_cor_order, drop = FALSE]

tumor_png <- file.path(figdir, "GSE221561_S1265T_InferCNA_manual_tumor_only_heatmap.png")

png(tumor_png, width = 4200, height = 2000, res = 300)
par(mar = c(5, 5, 4, 2))

zt <- t(tumor_heat)
zt <- pmax(pmin(zt, plot_limit), -plot_limit)

image(x = seq_len(ncol(zt)), y = seq_len(nrow(zt)), z = t(zt), axes = FALSE, xlab = "Chromosome", ylab = "S1265T Tumor Epithelial Cells", main = "S1265T Tumor Epithelial CNA Profiles", col = hcl.colors(101, "Blue-Red 3"))

for (b in chr_bounds$end_index) abline(v = b + 0.5, lwd = 0.4)
axis(1, at = chr_bounds$midpoint, labels = chr_bounds$chromosome, las = 2, cex.axis = 0.65)

dev.off()
cat("Saved: ", tumor_png, "\n", sep = "")

# ============================================================
# 20. FINAL REVIEW GATE
# ============================================================

cat("\n============================================================\n")
cat("17. GSE221561 S1265T PILOT REVIEW GATE\n")
cat("============================================================\n\n")

tumor_signal_med <- median(tumor_metrics$cna_signal, na.rm = TRUE)
normal_signal_med <- median(normal_metrics$cna_signal, na.rm = TRUE)
tumor_cor_med <- median(tumor_metrics$cna_cor_to_tumor, na.rm = TRUE)
normal_cor_med <- median(normal_metrics$cna_cor_to_tumor, na.rm = TRUE)

cat(sprintf("S1265T tumor median CNA signal        : %.5f\n", tumor_signal_med))
cat(sprintf("Normal median CNA signal              : %.5f\n", normal_signal_med))
cat(sprintf("CNA signal delta (Tumor - Normal)      : %.5f\n", tumor_signal_med - normal_signal_med))
cat(sprintf("CNA signal ratio (Tumor / Normal)      : %.3f\n", tumor_signal_med / normal_signal_med))
cat(sprintf("Tumor median cnaCor                    : %.3f\n", tumor_cor_med))
cat(sprintf("Normal median correlation to tumor     : %.3f\n", normal_cor_med))
cat(sprintf("Correlation separation                 : %.3f\n", tumor_cor_med - normal_cor_med))
cat(sprintf("Tumor > normal 95th CNA signal         : %.1f%%\n", 100 * mean(tumor_metrics$cna_signal > normal_signal_95, na.rm = TRUE)))
cat("Candidate coherent chromosome arms      : ", nrow(coherent), "\n", sep = "")

if (nrow(coherent) > 0) {
  cat("Candidate arm names                    : ", paste(coherent$chromosome_arm, collapse = ", "), "\n", sep = "")
}

# ============================================================
# 21. MANIFEST
# ============================================================

manifest <- data.table(
  dataset = "GSE221561",
  pilot = "S1265T",
  treatment_group = "Neoadjuvant_treated",
  tumor_cells = length(tumor_cells),
  normal_reference_cells = length(normal_cells),
  S6417N_cells = length(refCells[["S6417N"]]),
  S9478N_cells = length(refCells[["S9478N"]]),
  CNA_genes = nrow(cna),
  runtime_minutes = runtime_min,
  genome = "hg38", window = 100, n = 5000, noise = 0.1,
  center_method = "median", isLog = FALSE,
  tumor_median_CNA_signal = tumor_signal_med,
  normal_median_CNA_signal = normal_signal_med,
  tumor_median_cnaCor = tumor_cor_med,
  normal_median_corr_to_tumor = normal_cor_med,
  coherent_CNA_arms = nrow(coherent),
  malignant_calling = FALSE,
  status = "CNA_GENERATED_REVIEW_REQUIRED"
)

fwrite(manifest, file.path(outdir, "GSE221561_S1265T_17H3_manifest.csv"), bom = TRUE)

# ============================================================
# 22. RELOAD VALIDATION
# ============================================================

cat("\n===== 18. RELOAD CNA VALIDATION =====\n\n")

check <- readRDS(all_cna_file)
cat("Reloaded CNA: ", nrow(check), " genes x ", ncol(check), " cells\n", sep = "")
if (!identical(dim(check), dim(cna))) stop("Reloaded CNA dimensions changed.")
cat("CNA reload validated\n")

# ============================================================
# FINAL
# ============================================================

cat("\n============================================================\n")
cat("STEP 17H3 FINISHED\n")
cat("GSE221561 S1265T INFERCNA PILOT COMPLETE\n")
cat("============================================================\n\n")

cat("Tumor cells        : ", length(tumor_cells), "\n", sep = "")
cat("Normal refs        : ", length(normal_cells), " cells in 2 balanced groups\n", sep = "")
cat("CNA matrix         : ", nrow(cna), " genes x ", ncol(cna), " cells\n", sep = "")
cat("Runtime            : ", sprintf("%.2f min", runtime_min), "\n", sep = "")

cat("\nCNA profiles generated          \n")
cat("Malignant labels NOT assigned  \n")
cat("findMalignant() NOT run        \n")
cat("No tumor cells removed         \n")
cat("No other GSE221561 samples run \n")
cat("GSE197677 remains frozen       \n")

cat("\nMain heatmap:\n")
cat(heat_png, "\n")
cat("\nTumor-only heatmap:\n")
cat(tumor_png, "\n")
