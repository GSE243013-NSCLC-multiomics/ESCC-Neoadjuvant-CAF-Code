# ============================================================
# ESCC Neoadjuvant Project
# STEP 17F2
#
# Run InferCNA on SECOND GSE197677 pilot (ESC15).
#
# Identical parameters to ESC06 pilot:
#   window = 100, n = 5000, noise = 0.1
#   center.method = median, isLog = FALSE
#
# IMPORTANT:
#   NO parameter changes after seeing results
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
cat("STEP 17F2: RUN INFERCNA ON GSE197677 ESC15 PILOT\n")
cat("============================================================\n\n")

# ============================================================
# 1. PATHS
# ============================================================

base_outdir <- paste0(
  "06_tables/cross_dataset/",
  "malignant_epithelial"
)

pilot_dir <- file.path(
  base_outdir,
  "GSE197677_ESC15_pilot"
)

if (!dir.exists(pilot_dir)) {
  stop("Pilot directory not found: ", pilot_dir)
}

# ============================================================
# 2. LOAD PREPARED INPUTS
# ============================================================

cat("===== 1. LOAD PREPARED INPUTS =====\n\n")

cpm_file <- file.path(
  pilot_dir,
  "GSE197677_ESC15_infercna_CPM_hg38.rds"
)

refs_file <- file.path(
  pilot_dir,
  "GSE197677_ESC15_infercna_refCells.rds"
)

annotation_file <- file.path(
  pilot_dir,
  "GSE197677_ESC15_infercna_cell_annotation.csv"
)

manifest_file <- file.path(
  pilot_dir,
  "GSE197677_ESC15_infercna_manifest.csv"
)

cat("CPM file   : ", cpm_file, "\n", sep = "")
cat("Ref file   : ", refs_file, "\n", sep = "")

if (!file.exists(cpm_file)) {
  stop("Missing CPM file.")
}

if (!file.exists(refs_file)) {
  stop("Missing reference file.")
}

cpm <- readRDS(cpm_file)
refCells <- readRDS(refs_file)
annotation <- fread(annotation_file)
manifest <- fread(manifest_file)

cat("CPM matrix : ", nrow(cpm), " genes x ", ncol(cpm), " cells\n", sep = "")
cat("Ref groups : ", length(refCells), "\n", sep = "")

if (inherits(cpm, "dgCMatrix")) {
  cat("Converting sparse CPM to dense matrix...\n")
  cpm <- as.matrix(cpm)
  cat("Conversion complete.\n")
}

# ============================================================
# 3. CONFIGURE hg38
# ============================================================

cat("\n===== 2. CONFIGURE hg38 =====\n\n")

infercna::useGenome("hg38")

# ============================================================
# 4. RUN INFERCNA
# ============================================================

cat("\n===== 3. RUN INFERCNA =====\n\n")

cat("Parameters (identical to ESC06):\n")
cat("  window          = 100\n")
cat("  n               = 5000\n")
cat("  noise           = 0.1\n")
cat("  center.method   = median\n")
cat("  isLog           = FALSE\n")
cat("  range           = c(-3, 3)\n\n")

t_start <- proc.time()

cna <- tryCatch(
  infercna::infercna(
    m = cpm,
    refCells = refCells,
    isLog = FALSE,
    window = 100,
    n = 5000,
    noise = 0.1,
    center.method = "median",
    range = c(-3, 3)
  ),
  error = function(e) {
    cat("ERROR: ", conditionMessage(e), "\n", sep = "")
    return(NULL)
  }
)

t_end <- proc.time()
elapsed <- (t_end - t_start)[[3]]

cat(
  sprintf("\nInferCNA completed in %.2f min\n", elapsed / 60)
)

if (is.null(cna)) {
  stop("InferCNA failed.")
}

# ============================================================
# 5. INSPECT RESULT
# ============================================================

cat("\n===== 4. INSPECT CNA RESULT =====\n\n")

cat("Class      : ", class(cna), "\n", sep = "")
cat("Dimensions : ", paste(dim(cna), collapse = " x "), "\n", sep = "")
cat("Range      : ", sprintf("%.3f to %.3f", min(cna), max(cna)), "\n", sep = "")

cat("\nFirst 5 genes (rows):\n")
print(head(rownames(cna), 5))

cat("\nFirst 10 cells (cols):\n")
print(head(colnames(cna), 10))

# ============================================================
# 6. TUMOR vs NORMAL CNA SIGNAL
# ============================================================

cat("\n===== 5. TUMOR vs NORMAL CNA SIGNAL =====\n\n")

tumor_cells <- annotation[infercna_role == "TUMOR_TEST", cell_id]
normal_cells <- annotation[infercna_role == "NORMAL_REFERENCE", cell_id]

cat("Tumor cells in CNA  : ", sum(tumor_cells %in% colnames(cna)), "/", length(tumor_cells), "\n", sep = "")
cat("Normal cells in CNA: ", sum(normal_cells %in% colnames(cna)), "/", length(normal_cells), "\n", sep = "")

tumor_cna <- cna[, tumor_cells[tumor_cells %in% colnames(cna)], drop = FALSE]
normal_cna <- cna[, normal_cells[normal_cells %in% colnames(cna)], drop = FALSE]

tumor_mean <- rowMeans(tumor_cna, na.rm = TRUE)
normal_mean <- rowMeans(normal_cna, na.rm = TRUE)

cat("\nTumor mean CNA (per gene):\n")
cat("  Mean   : ", sprintf("%.4f", mean(tumor_mean, na.rm = TRUE)), "\n", sep = "")
cat("  Median : ", sprintf("%.4f", median(tumor_mean, na.rm = TRUE)), "\n", sep = "")
cat("  Range  : ", sprintf("%.4f to %.4f", min(tumor_mean, na.rm = TRUE), max(tumor_mean, na.rm = TRUE)), "\n", sep = "")

cat("\nNormal mean CNA (per gene):\n")
cat("  Mean   : ", sprintf("%.4f", mean(normal_mean, na.rm = TRUE)), "\n", sep = "")
cat("  Median : ", sprintf("%.4f", median(normal_mean, na.rm = TRUE)), "\n", sep = "")
cat("  Range  : ", sprintf("%.4f to %.4f", min(normal_mean, na.rm = TRUE), max(normal_mean, na.rm = TRUE)), "\n", sep = "")

# ============================================================
# 7. CORRELATION TEST (TUMOR PROFILE)
# ============================================================

cat("\n===== 6. TUMOR-PROFILE CORRELATION =====\n\n")

tumor_profile <- rowMeans(tumor_cna, na.rm = TRUE)
normal_profile <- rowMeans(normal_cna, na.rm = TRUE)

cor_tumor <- cor(tumor_profile, normal_profile, use = "complete.obs")
cat("Tumor vs Normal profile correlation: ", sprintf("%.4f", cor_tumor), "\n", sep = "")

# Per-cell correlation with tumor mean
cell_cor_tumor <- sapply(
  seq_len(ncol(tumor_cna)),
  function(i) {
    cor(tumor_cna[, i], tumor_profile, use = "complete.obs")
  }
)

cat("\nPer-cell correlation with tumor mean:\n")
cat("  Mean   : ", sprintf("%.4f", mean(cell_cor_tumor, na.rm = TRUE)), "\n", sep = "")
cat("  Median : ", sprintf("%.4f", median(cell_cor_tumor, na.rm = TRUE)), "\n", sep = "")
cat("  Min    : ", sprintf("%.4f", min(cell_cor_tumor, na.rm = TRUE)), "\n", sep = "")
cat("  Max    : ", sprintf("%.4f", max(cell_cor_tumor, na.rm = TRUE)), "\n", sep = "")

# ============================================================
# 8. CNA HEATMAP (FULL, CHROMOSOME-ORDERED)
# ============================================================

cat("\n===== 7. CNA HEATMAP =====\n\n")

fig_dir <- paste0(
  "05_figures/GSE197677/malignant_epithelial/ESC15_pilot"
)

dir.create(
  fig_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

genome <- infercna::retrieveGenome(name = "hg38")

gene_info <- data.table(
  gene = rownames(cna)
)

genome_sub <- genome[
  match(gene_info$gene, genome$symbol),
]

gene_info[, chr := as.character(genome_sub$chromosome_name)]
gene_info[, arm := as.character(genome_sub$arm)]

chr_order <- paste0("chr", c(as.character(1:22), "X", "Y"))
gene_info[, chr_factor := factor(paste0("chr", chr), levels = chr_order)]

setorder(gene_info, chr_factor, arm, gene)

cna_ordered <- cna[gene_info$gene, , drop = FALSE]

# Save heatmap
png(
  file.path(
    fig_dir,
    "ESC15_CNA_heatmap_full.png"
  ),
  width = 1600,
  height = 1200,
  res = 150
)

image(
  t(cna_ordered),
  col = colorRampPalette(
    c("blue", "white", "red")
  )(100),
  zlim = c(-1.5, 1.5),
  xaxt = "n",
  yaxt = "n",
  main = paste0(
    "ESC15 InferCNA Heatmap\n",
    nrow(cna_ordered), " genes x ",
    ncol(cna_ordered), " cells"
  ),
  xlab = "Cells",
  ylab = "Genes (chromosome order)"
)

# Add chromosome labels
chr_starts <- which(
  !duplicated(gene_info$chr_factor)
)

chr_midpoints <- tapply(
  seq_len(nrow(gene_info)),
  gene_info$chr_factor,
  function(x) {
    (min(x) + max(x)) / 2 / nrow(gene_info)
  }
)

axis(
  2,
  at = as.numeric(chr_midpoints),
  labels = names(chr_midpoints),
  cex.axis = 0.5,
  las = 2
)

# Add tumor/normal color bar
annotation_ordered <- annotation[
  match(colnames(cna_ordered), cell_id),
]

cell_colors <- ifelse(
  annotation_ordered$infercna_role == "TUMOR_TEST",
  "red",
  "steelblue"
)

par(xpd = TRUE)
rect(
  0,
  -0.02,
  1,
  0
)

for (i in seq_along(cell_colors)) {
  rect(
    (i - 1) / length(cell_colors),
    -0.02,
    i / length(cell_colors),
    0,
    col = cell_colors[i],
    border = NA
  )
}

legend(
  "bottomright",
  legend = c("Tumor", "Normal ref"),
  fill = c("red", "steelblue"),
  cex = 0.6,
  inset = c(0, -0.05)
)

dev.off()

cat("Heatmap saved: ", file.path(fig_dir, "ESC15_CNA_heatmap_full.png"), "\n", sep = "")

# ============================================================
# 9. PER-CHROMOSOME CNA SUMMARY
# ============================================================

cat("\n===== 8. PER-CHROMOSOME CNA SUMMARY =====\n\n")

chr_summary <- data.table(
  chromosome = character(),
  mean_tumor = numeric(),
  mean_normal = numeric(),
  n_genes = integer()
)

for (ch in unique(gene_info$chr_factor)) {

  if (is.na(ch)) next

  genes_in_chr <- gene_info[chr_factor == ch, gene]

  if (length(genes_in_chr) < 10) next

  tumor_vals <- colMeans(
    cna[genes_in_chr, tumor_cells[tumor_cells %in% colnames(cna)], drop = FALSE],
    na.rm = TRUE
  )

  normal_vals <- colMeans(
    cna[genes_in_chr, normal_cells[normal_cells %in% colnames(cna)], drop = FALSE],
    na.rm = TRUE
  )

  chr_summary <- rbind(
    chr_summary,
    data.table(
      chromosome = as.character(ch),
      mean_tumor = mean(tumor_vals),
      mean_normal = mean(normal_vals),
      n_genes = length(genes_in_chr)
    )
  )
}

chr_summary[, delta := mean_tumor - mean_normal]

chr_summary[, abs_delta := abs(delta)]
setorder(chr_summary, -abs_delta)

cat("Top 10 chromosomes by |delta|:\n")
print(head(chr_summary, 10), digits = 4)

# ============================================================
# 10. SAVE RESULTS
# ============================================================

cat("\n===== 9. SAVE RESULTS =====\n\n")

cna_rds <- file.path(
  pilot_dir,
  "GSE197677_ESC15_infercna_CNA_matrix.rds"
)

chr_rds <- file.path(
  pilot_dir,
  "GSE197677_ESC15_CNA_chr_summary.csv"
)

cor_file <- file.path(
  pilot_dir,
  "GSE197677_ESC15_CNA_tumor_correlations.csv"
)

saveRDS(cna, cna_rds)
fwrite(chr_summary, chr_rds, bom = TRUE)

cor_dt <- data.table(
  cell_id = names(cell_cor_tumor),
  correlation_with_tumor_mean = cell_cor_tumor
)

fwrite(cor_dt, cor_file, bom = TRUE)

# ============================================================
# 11. UPDATE MANIFEST
# ============================================================

cat("\n===== 10. UPDATE MANIFEST =====\n\n")

manifest[, infercna_run := TRUE]
manifest[, infercna_runtime_min := round(elapsed / 60, 2)]
manifest[, cna_genes := nrow(cna)]
manifest[, tumor_profile_correlation := round(cor_tumor, 4)]
manifest[, output_cna_matrix := cna_rds]

fwrite(
  manifest,
  manifest_file,
  bom = TRUE
)

# ============================================================
# CLEANUP
# ============================================================

rm(cpm, cna, tumor_cna, normal_cna)
gc(verbose = FALSE)

# ============================================================
# FINAL
# ============================================================

cat("\n============================================================\n")
cat("STEP 17F2 FINISHED\n")
cat("SECOND GSE197677 CNV PILOT COMPLETE\n")
cat("============================================================\n\n")

cat("Pilot sample       : ESC15\n")
cat("Treatment group    : Neoadjuvant_chemotherapy\n")
cat("Tumor cells        : ", length(tumor_cells), "\n", sep = "")
cat("Normal ref cells   : ", length(normal_cells), "\n", sep = "")
cat("CNA genes          : ", nrow(cna), "\n", sep = "")
cat("Runtime            : ", sprintf("%.2f min", elapsed / 60), "\n", sep = "")
cat("Tumor profile corr : ", sprintf("%.4f", cor_tumor), "\n", sep = "")

cat("\nCAUTION:\n")
cat("  This is a SECOND independent pilot.\n")
cat("  Do NOT combine results yet.\n")
cat("  Do NOT call malignant.\n")
cat("  STOP HERE for manual review.\n")
