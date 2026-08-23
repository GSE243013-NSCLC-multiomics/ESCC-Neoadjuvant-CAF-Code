# ============================================================
# ESCC Neoadjuvant Project - STEP 10
# QC audit of author-retained GSE221561 cells
# ============================================================

options(stringsAsFactors = FALSE, timeout = 1200)

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(ggplot2)
})

cat("\n============================================================\n")
cat("STEP 10: GSE221561 QC AUDIT\n")
cat("============================================================\n\n")

# Paths
object_file <- "03_objects/GSE221561/GSE221561_author_annotated.rds"
result_dir <- "04_results/GSE221561/QC"
figure_dir <- "05_figures/GSE221561/QC"
table_dir  <- "06_tables/GSE221561/QC"
log_dir    <- "08_logs"

for (d in c(result_dir, figure_dir, table_dir, log_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# 1. Load
cat("===== 1. Loading annotated object =====\n\n")
if (!file.exists(object_file)) stop("Cannot find: ", object_file)
obj <- readRDS(object_file)
cat("Genes:", nrow(obj), "| Cells:", ncol(obj), "| Samples:", length(unique(obj$sample_id)), "\n")

# 2. Raw to author retention
cat("\n===== 2. RAW -> AUTHOR RETENTION =====\n\n")
match_file <- "00_metadata/05_GSE221561_seurat_metadata_match_audit.csv"
if (file.exists(match_file)) {
  retention <- fread(match_file)
  retention[, author_retention_percent := round(100 * metadata_matched / matrix_cells, 2)]
  retention_out <- retention[, .(sample_id, raw_matrix_cells = matrix_cells, author_retained_cells = metadata_matched, raw_not_in_author_metadata = metadata_unmatched, author_retention_percent)]
  fwrite(retention_out, file.path(table_dir, "GSE221561_raw_to_author_retention.csv"), bom = TRUE)
  print(retention_out)
  cat(sprintf("\nOverall author retention: %.2f%%\n", 100 * sum(retention_out$author_retained_cells) / sum(retention_out$raw_matrix_cells)))
}

# 3. Check mitochondrial genes
cat("\n===== 3. MITOCHONDRIAL GENE CHECK =====\n\n")
gene_names <- rownames(obj)
mt_upper <- grep("^MT-", gene_names, value = TRUE)
mt_lower <- grep("^mt-", gene_names, value = TRUE)
if (length(mt_upper) > 0) { mt_pattern <- "^MT-"; mt_genes <- mt_upper
} else if (length(mt_lower) > 0) { mt_pattern <- "^mt-"; mt_genes <- mt_lower
} else { mt_pattern <- NA_character_; mt_genes <- character(0) }
cat("Mitochondrial genes:", length(mt_genes), "\n")
if (length(mt_genes) > 0) { cat("Pattern:", mt_pattern, "\n"); obj$percent.mt.recomputed <- PercentageFeatureSet(obj, pattern = mt_pattern)
} else { cat("WARNING: no mitochondrial genes detected.\n"); obj$percent.mt.recomputed <- NA_real_ }

# 4. Build QC table
cat("\n===== 4. BUILDING QC TABLE =====\n\n")
qc <- as.data.table(obj[[]], keep.rownames = "cell_id")
cat("QC rows:", nrow(qc), "\n")

# 5. Detect author QC fields
cat("\n===== 5. AUTHOR QC FIELD DETECTION =====\n\n")
author_fields <- grep("^author_", names(qc), value = TRUE)
find_author_field <- function(patterns) author_fields[Reduce(`|`, lapply(patterns, function(p) grepl(p, author_fields, ignore.case = TRUE)))]

author_ncount <- find_author_field(c("nCount", "UMI"))
author_nfeature <- find_author_field(c("nFeature", "gene"))
author_mt <- author_fields[grepl("percent.*mt|mito", author_fields, ignore.case = TRUE)]

cat("Author nCount:", paste(author_ncount, collapse = ", "), "\n")
cat("Author nFeature:", paste(author_nfeature, collapse = ", "), "\n")
cat("Author mt:", paste(author_mt, collapse = ", "), "\n")

# 6. Compare metrics
cat("\n===== 6. QC RECOMPUTATION CHECK =====\n\n")
comparison_results <- list()
compare_metric <- function(current_col, author_col, metric_name) {
  x <- suppressWarnings(as.numeric(qc[[current_col]]))
  y <- suppressWarnings(as.numeric(qc[[author_col]]))
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) == 0) return(NULL)
  data.table(metric = metric_name, current_column = current_col, author_column = author_col, compared_cells = sum(keep), pearson_r = cor(x[keep], y[keep], method = "pearson"), spearman_r = cor(x[keep], y[keep], method = "spearman"), median_abs_diff = median(abs(x[keep] - y[keep])), max_abs_diff = max(abs(x[keep] - y[keep])))
}

if (length(author_ncount) > 0) { tmp <- compare_metric("nCount_RNA", author_ncount[1], "nCount_RNA"); if (!is.null(tmp)) comparison_results[["nCount"]] <- tmp }
if (length(author_nfeature) > 0) { tmp <- compare_metric("nFeature_RNA", author_nfeature[1], "nFeature_RNA"); if (!is.null(tmp)) comparison_results[["nFeature"]] <- tmp }
if (length(author_mt) > 0 && any(is.finite(qc$percent.mt.recomputed))) { tmp <- compare_metric("percent.mt.recomputed", author_mt[1], "percent_mt"); if (!is.null(tmp)) comparison_results[["percent_mt"]] <- tmp }

if (length(comparison_results) > 0) { comparison_table <- rbindlist(comparison_results, fill = TRUE); fwrite(comparison_table, file.path(table_dir, "GSE221561_author_vs_recomputed_QC.csv"), bom = TRUE); print(comparison_table) }

# 7. Per-sample QC summary
cat("\n===== 7. QC SUMMARY BY SAMPLE =====\n\n")
qc_summary <- qc[, .(n_cells = .N, median_nCount_RNA = median(nCount_RNA, na.rm = TRUE), median_nFeature_RNA = median(nFeature_RNA, na.rm = TRUE), median_percent_mt = median(percent.mt.recomputed, na.rm = TRUE), mean_nCount_RNA = mean(nCount_RNA, na.rm = TRUE), mean_nFeature_RNA = mean(nFeature_RNA, na.rm = TRUE), mean_percent_mt = mean(percent.mt.recomputed, na.rm = TRUE)), by = sample_id]
setorder(qc_summary, sample_id)
print(qc_summary)
fwrite(qc_summary, file.path(table_dir, "GSE221561_QC_summary_by_sample.csv"), bom = TRUE)

# 8. Quantiles
cat("\n===== 8. QC QUANTILES =====\n\n")
quantile_levels <- c(0, 0.01, 0.05, 0.25, 0.50, 0.75, 0.95, 0.99, 1)
quantile_names <- c("min", "p01", "p05", "p25", "p50", "p75", "p95", "p99", "max")
make_quantiles <- function(x) { q <- quantile(x, probs = quantile_levels, na.rm = TRUE, names = FALSE); as.list(setNames(as.numeric(q), quantile_names)) }

ncount_q <- qc[, make_quantiles(nCount_RNA), by = sample_id]; ncount_q[, metric := "nCount_RNA"]
nfeature_q <- qc[, make_quantiles(nFeature_RNA), by = sample_id]; nfeature_q[, metric := "nFeature_RNA"]
if (any(is.finite(qc$percent.mt.recomputed))) { mt_q <- qc[, make_quantiles(percent.mt.recomputed), by = sample_id]; mt_q[, metric := "percent.mt"]; quantile_table <- rbindlist(list(ncount_q, nfeature_q, mt_q), fill = TRUE)
} else { quantile_table <- rbindlist(list(ncount_q, nfeature_q), fill = TRUE) }
setcolorder(quantile_table, c("sample_id", "metric", quantile_names))
fwrite(quantile_table, file.path(table_dir, "GSE221561_QC_quantiles_by_sample.csv"), bom = TRUE)
print(quantile_table)

# 9. Overall quantiles
overall <- rbindlist(list(data.table(metric = "nCount_RNA", t(as.data.frame(make_quantiles(qc$nCount_RNA)))), data.table(metric = "nFeature_RNA", t(as.data.frame(make_quantiles(qc$nFeature_RNA))))), fill = TRUE)
if (any(is.finite(qc$percent.mt.recomputed))) { overall_mt <- data.table(metric = "percent.mt", t(as.data.frame(make_quantiles(qc$percent.mt.recomputed)))); overall <- rbindlist(list(overall, overall_mt), fill = TRUE) }
fwrite(overall, file.path(table_dir, "GSE221561_QC_quantiles_overall.csv"), bom = TRUE)
cat("\n===== 9. OVERALL QC QUANTILES =====\n\n")
print(overall)

# 10. Save cell-level QC
fwrite(qc[, .(cell_id, sample_id, nCount_RNA, nFeature_RNA, percent.mt.recomputed)], file.path(table_dir, "GSE221561_cell_level_QC.csv.gz"))

# 11. QC plots
cat("\n===== 10. GENERATING QC FIGURES =====\n\n")
sample_order <- sort(unique(qc$sample_id))
qc[, sample_id_plot := factor(sample_id, levels = sample_order)]

p_feature <- ggplot(qc, aes(x = sample_id_plot, y = nFeature_RNA)) + geom_violin(scale = "width", trim = TRUE) + geom_boxplot(width = 0.12, outlier.shape = NA) + labs(title = "GSE221561: Detected Genes per Cell", x = "Sample", y = "nFeature_RNA") + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(figure_dir, "GSE221561_QC_nFeature_by_sample.pdf"), p_feature, width = 10, height = 6)
ggsave(file.path(figure_dir, "GSE221561_QC_nFeature_by_sample.png"), p_feature, width = 10, height = 6, dpi = 300)

p_count <- ggplot(qc, aes(x = sample_id_plot, y = nCount_RNA)) + geom_violin(scale = "width", trim = TRUE) + geom_boxplot(width = 0.12, outlier.shape = NA) + labs(title = "GSE221561: UMI Counts per Cell", x = "Sample", y = "nCount_RNA") + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(figure_dir, "GSE221561_QC_nCount_by_sample.pdf"), p_count, width = 10, height = 6)
ggsave(file.path(figure_dir, "GSE221561_QC_nCount_by_sample.png"), p_count, width = 10, height = 6, dpi = 300)

if (any(is.finite(qc$percent.mt.recomputed))) {
  p_mt <- ggplot(qc, aes(x = sample_id_plot, y = percent.mt.recomputed)) + geom_violin(scale = "width", trim = TRUE) + geom_boxplot(width = 0.12, outlier.shape = NA) + labs(title = "GSE221561: Mitochondrial RNA Fraction", x = "Sample", y = "Mitochondrial RNA (%)") + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(figure_dir, "GSE221561_QC_percentMT_by_sample.pdf"), p_mt, width = 10, height = 6)
  ggsave(file.path(figure_dir, "GSE221561_QC_percentMT_by_sample.png"), p_mt, width = 10, height = 6, dpi = 300)
}

p_scatter <- ggplot(qc, aes(x = nCount_RNA, y = nFeature_RNA)) + geom_point(alpha = 0.15, size = 0.25) + labs(title = "GSE221561: UMI Counts vs Detected Genes", x = "nCount_RNA", y = "nFeature_RNA") + theme_classic()
ggsave(file.path(figure_dir, "GSE221561_QC_nCount_vs_nFeature.pdf"), p_scatter, width = 7, height = 6)
ggsave(file.path(figure_dir, "GSE221561_QC_nCount_vs_nFeature.png"), p_scatter, width = 7, height = 6, dpi = 300)

gc(verbose = FALSE)

# 12. Session info
sink(file.path(log_dir, "07_GSE221561_qc_audit_sessionInfo.txt"))
print(sessionInfo())
sink()

cat("\n============================================================\n")
cat("STEP 10 FINISHED ✓\n")
cat("============================================================\n\n")
cat("Cells audited:", nrow(qc), "\n")
cat("No cells removed ✓\nNo normalization performed ✓\nNo clustering performed ✓\n")
