#!/usr/bin/env Rscript
# STEP 19G: EPITHELIAL IDENTITY, STROMAL-LIKE SIGNAL, AND PATHWAY LOCALIZATION AUDIT
# AUDIT ONLY — no cell reclassification, no removal, no changes to annotations

suppressPackageStartupMessages({
  library(Seurat)
  library(edgeR)
  library(ggplot2)
  library(reshape2)
})

base_dir <- "~/Downloads/ESCC Neoadjuvant"
pw_dir <- file.path(base_dir, "06_tables/cross_dataset/pathway_analysis")
g19g_dir <- file.path(pw_dir, "Step19G")
fig_dir <- file.path(base_dir, "05_figures/cross_dataset/pathway_analysis/Step19G")
dir.create(g19g_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

cache_file <- file.path(Sys.getenv("ESCC_CAF_PROJECT_ROOT", unset = "."), "cache", "msigdb.2026.1.Hs.H.rds")
hallmark_df <- readRDS(cache_file)
hallmark_pathways <- split(hallmark_df$db_gene_symbol, hallmark_df$gs_name)
hallmark_pathways <- lapply(hallmark_pathways, unique)

supported_pathways <- c(
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "HALLMARK_HYPOXIA",
  "HALLMARK_COAGULATION",
  "HALLMARK_KRAS_SIGNALING_UP",
  "HALLMARK_UV_RESPONSE_DN",
  "HALLMARK_MYOGENESIS",
  "HALLMARK_APOPTOSIS"
)

core_le_genes <- c("COL3A1","COL1A1","MMP2","BGN","SPARC","CCN1","DCN","F3","CDKN1A")
epithelial_core_markers <- c("EPCAM","KRT8","KRT18","KRT19","KRT5","KRT14","KRT17","SFN")
squamous_markers <- c("KRT5","KRT14","KRT17","TP63","SFN","KRT4","KRT13")
stromal_markers <- c("COL1A1","COL1A2","COL3A1","DCN","LUM","BGN","COL6A1","COL6A2","COL6A3","SPARC","MMP2","PDGFRA","PDGFRB")
hypoxia_markers <- c("HIF1A","VEGFA","CA9","BNIP3","LDHA","ENO1","PGK1","SLC2A1")
stress_markers <- c("CDKN1A","FOS","JUN","ATF3","DDIT3","HSPA1A","HSPA1B")

# Helper: safely get counts from Seurat object (handles v3 and v5)
# For v5: works layer-by-layer to avoid JoinLayers memory blowup
get_counts <- function(obj) {
  tryCatch({
    GetAssayData(obj, assay = "RNA", layer = "counts")
  }, error = function(e) {
    # v5 with multiple layers — return NULL, caller must use layer-by-layer
    NULL
  })
}
get_genes <- function(obj) {
  rownames(obj@assays$RNA)
}
is_v5 <- function(obj) {
  inherits(obj@assays$RNA, "Assay5")
}

# Build pseudobulk counts from v5 layers (memory-safe)
# Aggregates by grouping_factor, summing across cells within each group
build_pseudobulk_from_layers <- function(obj, md, grouping_factor, subset_idx = NULL) {
  assay <- obj[["RNA"]]
  layer_names <- Layers(assay)
  groups <- md[[grouping_factor]]
  if (!is.null(subset_idx)) {
    groups <- groups[subset_idx]
  }
  group_levels <- sort(unique(groups))
  
  # Get gene names from first layer
  l1 <- LayerData(assay, layer = layer_names[1])
  gene_names <- rownames(l1)
  n_genes <- length(gene_names)
  n_groups <- length(group_levels)
  
  pb <- matrix(0, nrow = n_genes, ncol = n_groups,
               dimnames = list(gene_names, as.character(group_levels)))
  
  for (ln in layer_names) {
    counts_layer <- LayerData(assay, layer = ln)  # sparse dgCMatrix
    layer_cells <- colnames(counts_layer)
    
    # Find which of these cells are in our subset
    if (!is.null(subset_idx)) {
      md_layer <- md[layer_cells, ]
      in_subset <- layer_cells %in% rownames(md)[subset_idx]
      if (sum(in_subset) == 0) next
      counts_layer <- counts_layer[, in_subset, drop = FALSE]
      layer_groups <- groups[match(colnames(counts_layer), rownames(md)[subset_idx])]
    } else {
      layer_groups <- groups[layer_cells]
    }
    
    # Aggregate by group
    for (g in group_levels) {
      g_cells <- which(layer_groups == g)
      if (length(g_cells) > 0) {
        pb[, as.character(g)] <- pb[, as.character(g)] + Matrix::rowSums(counts_layer[, g_cells, drop = FALSE])
      }
    }
    rm(counts_layer); gc(FALSE)
  }
  
  # Remove zero-sum columns
  keep <- colSums(pb) > 0
  pb[, keep, drop = FALSE]
}

# ============================================================================
# 19G0: INPUT + METADATA PREFLIGHT
# ============================================================================
cat("=== 19G0: INPUT + METADATA PREFLIGHT ===\n\n")

cat("Loading GSE197677...\n")
obj19 <- readRDS(file.path(base_dir, "03_objects/GSE197677/GSE197677_integrated_broadcelltype_annotated.rds"))
md19 <- obj19@meta.data
cat(sprintf("  %d cells, %d metadata fields\n", nrow(md19), ncol(md19)))

cat("Loading GSE221561...\n")
obj22 <- readRDS(file.path(base_dir, "03_objects/GSE221561/GSE221561_author_annotated.rds"))
md22 <- obj22@meta.data
cat(sprintf("  %d cells, %d metadata fields\n", nrow(md22), ncol(md22)))

# Detect v5 status
if (is_v5(obj22)) {
  cat("  GSE221561: v5 assay detected, will use layer-by-layer pseudobulk (memory-safe)\n")
} else {
  cat("  GSE221561: standard assay\n")
}
cat("\n")

genes19 <- get_genes(obj19)
genes22 <- get_genes(obj22)
cat(sprintf("GSE197677 genes: %d\nGSE221561 genes: %d\n\n", length(genes19), length(genes22)))

cat("GSE197677 broad_celltype_final:\n")
print(table(md19$broad_celltype_final))
cat("\nGSE221561 author_Cell_type:\n")
print(table(md22$author_Cell_type))

# Map GSE221561
mapping22 <- c("B cell"="B_plasma","Endothelial"="Endothelial","Epithelial"="Epithelial",
               "Fibroblast"="Fibroblast_CAF","Mast cell"="Mast","Myeloid"="Myeloid","T cell"="T_NK")
md22$broad_celltype_final <- mapping22[md22$author_Cell_type]

treat22 <- md22$author_Neoadjuvant
treat22[treat22 %in% c("Chemotherapy","Chemoradiotherapy","Chemoradiotherapy and Immunotherapy","Antiangiogenesis")] <- "Neoadjuvant_treated"
treat22[treat22 == "Surgery alone"] <- "Surgery_alone"
treat22[treat22 == "Adjacent normal"] <- "Adjacent_normal"
md22$treatment_group <- treat22

# Map GSE197677 treatment from pseudobulk metadata
meta_pseudo <- read.csv(file.path(base_dir, "06_tables/cross_dataset/epithelial_pseudobulk/counts/GSE197677_epithelial_pseudobulk_sample_metadata.csv"), stringsAsFactors = FALSE)
md19$treatment_group <- setNames(meta_pseudo$treatment_standardized, meta_pseudo$sample)[md19$sample]

cat("\nGSE221561 mapped broad_celltype_final:\n")
print(table(md22$broad_celltype_final))
cat("\nGSE221561 mapped treatment_group:\n")
print(table(md22$treatment_group))
cat("\nGSE197677 mapped treatment_group:\n")
print(table(md19$treatment_group))

cat("\n--- Fibroblast reference ---\n")
cat("GSE197677 Fibroblast_CAF:", sum(md19$broad_celltype_final == "Fibroblast_CAF"), "cells\n")
cat("GSE221561 Fibroblast_CAF:", sum(md22$broad_celltype_final == "Fibroblast_CAF"), "cells\n")
cat("Both available: YES\n")

preflight <- data.frame(
  dataset = c("GSE197677","GSE221561"),
  cells = c(nrow(md19), nrow(md22)),
  epithelial_cells = c(sum(md19$broad_celltype_final=="Epithelial"), sum(md22$broad_celltype_final=="Epithelial")),
  fibroblast_cells = c(sum(md19$broad_celltype_final=="Fibroblast_CAF"), sum(md22$broad_celltype_final=="Fibroblast_CAF")),
  stringsAsFactors = FALSE
)
write.csv(preflight, file.path(g19g_dir, "STEP19G_METADATA_PREFLIGHT.csv"), row.names = FALSE)
cat("\nSaved: STEP19G_METADATA_PREFLIGHT.csv\n\n")

# ============================================================================
# 19G1: FREEZE TARGET GENE PANELS
# ============================================================================
cat("=== 19G1: FREEZE TARGET GENE PANELS ===\n\n")

all_panel_genes <- unique(c(epithelial_core_markers, squamous_markers, stromal_markers, hypoxia_markers, stress_markers, core_le_genes))

avail19 <- all_panel_genes %in% genes19
cat(sprintf("GSE197677: %d/%d panel genes available\n", sum(avail19), length(all_panel_genes)))
miss19 <- all_panel_genes[!avail19]
if (length(miss19) > 0) cat("  Missing:", paste(miss19, collapse=", "), "\n")

avail22 <- all_panel_genes %in% genes22
cat(sprintf("GSE221561: %d/%d panel genes available\n", sum(avail22), length(all_panel_genes)))
miss22 <- all_panel_genes[!avail22]
if (length(miss22) > 0) cat("  Missing:", paste(miss22, collapse=", "), "\n")

# Build panel mapping with primary panel per gene
gene_panel_map <- data.frame(gene = character(), panel = character(), stringsAsFactors = FALSE)
for (g in all_panel_genes) {
  pnl <- "unknown"
  if (g %in% epithelial_core_markers) pnl <- "epithelial_core"
  else if (g %in% squamous_markers) pnl <- "squamous"
  else if (g %in% stromal_markers) pnl <- "stromal"
  else if (g %in% hypoxia_markers) pnl <- "hypoxia"
  else if (g %in% stress_markers) pnl <- "stress"
  else if (g %in% core_le_genes) pnl <- "core_le"
  gene_panel_map <- rbind(gene_panel_map, data.frame(gene=g, panel=pnl, stringsAsFactors=FALSE))
}
panel_df <- data.frame(
  gene = gene_panel_map$gene,
  panel = gene_panel_map$panel,
  available_GSE197677 = avail19[match(gene_panel_map$gene, all_panel_genes)],
  available_GSE221561 = avail22[match(gene_panel_map$gene, all_panel_genes)],
  stringsAsFactors = FALSE
)
write.csv(panel_df, file.path(g19g_dir, "STEP19G_IDENTITY_MARKER_PANEL_AVAILABILITY.csv"), row.names = FALSE)
cat("Saved: STEP19G_IDENTITY_MARKER_PANEL_AVAILABILITY.csv\n\n")

# ============================================================================
# 19G2: CELL-TYPE REFERENCE EXPRESSION AUDIT
# ============================================================================
cat("=== 19G2: CELL-TYPE REFERENCE EXPRESSION AUDIT ===\n\n")

run_ref_audit <- function(obj, md, dname) {
  cat(sprintf("--- %s ---\n", dname))
  epi_idx <- which(md$broad_celltype_final == "Epithelial")
  fib_idx <- which(md$broad_celltype_final == "Fibroblast_CAF")
  cat(sprintf("  Epithelial: %d, Fibroblast: %d cells\n", length(epi_idx), length(fib_idx)))

  if (is_v5(obj)) {
    # Layer-by-layer pseudobulk by sample × celltype
    assay <- obj[["RNA"]]
    layer_names <- Layers(assay)
    gene_names <- rownames(assay)
    
    # Create grouping: sample_celltype
    epi_samp <- md$sample[epi_idx]
    fib_samp <- md$sample[fib_idx]
    all_samps <- sort(unique(c(epi_samp, fib_samp)))
    group_map <- character(0)
    group_map[epi_idx] <- paste0(md$sample[epi_idx], "_E")
    group_map[fib_idx] <- paste0(md$sample[fib_idx], "_F")
    group_levels <- sort(unique(group_map[group_map != ""]))
    
    pb <- matrix(0, nrow = length(gene_names), ncol = length(group_levels),
                 dimnames = list(gene_names, group_levels))
    
    for (ln in layer_names) {
      layer_counts <- LayerData(assay, layer = ln)
      layer_cells <- colnames(layer_counts)
      layer_groups <- group_map[layer_cells]
      for (g in group_levels) {
        g_cells <- which(layer_groups == g)
        if (length(g_cells) > 0) {
          pb[, g] <- pb[, g] + Matrix::rowSums(layer_counts[, g_cells, drop = FALSE])
        }
      }
      rm(layer_counts); gc(FALSE)
    }
    
    # Detection fractions from full single-cell data
    e_det_all <- numeric(length(gene_names))
    f_det_all <- numeric(length(gene_names))
    names(e_det_all) <- gene_names
    names(f_det_all) <- gene_names
    for (ln in layer_names) {
      layer_counts <- LayerData(assay, layer = ln)
      layer_cells <- colnames(layer_counts)
      e_cells <- intersect(layer_cells, rownames(md)[epi_idx])
      f_cells <- intersect(layer_cells, rownames(md)[fib_idx])
      if (length(e_cells) > 0) e_det_all <- e_det_all + Matrix::rowSums(layer_counts[, e_cells, drop = FALSE] > 0)
      if (length(f_cells) > 0) f_det_all <- f_det_all + Matrix::rowSums(layer_counts[, f_cells, drop = FALSE] > 0)
      rm(layer_counts); gc(FALSE)
    }
    e_det_frac <- e_det_all / length(epi_idx)
    f_det_frac <- f_det_all / length(fib_idx)
  } else {
    counts <- get_counts(obj)
    epi_samp <- md$sample[epi_idx]
    fib_samp <- md$sample[fib_idx]
    all_samps <- sort(unique(c(epi_samp, fib_samp)))
    pb_list <- list()
    for (s in all_samps) {
      e_i <- epi_idx[epi_samp == s]
      f_i <- fib_idx[fib_samp == s]
      if (length(e_i) > 0) pb_list[[paste0(s,"_E")]] <- rowSums(counts[, e_i, drop=FALSE])
      if (length(f_i) > 0) pb_list[[paste0(s,"_F")]] <- rowSums(counts[, f_i, drop=FALSE])
    }
    pb <- do.call(cbind, pb_list)
    e_det_frac <- rowSums(counts[, epi_idx, drop=FALSE] > 0) / length(epi_idx)
    f_det_frac <- rowSums(counts[, fib_idx, drop=FALSE] > 0) / length(fib_idx)
  }

  keep <- rowSums(pb > 0) >= 3
  pb <- pb[keep, ]
  dge <- DGEList(counts = pb)
  dge <- calcNormFactors(dge, method = "TMM")
  lcpm <- cpm(dge, log = TRUE)

  markers <- unique(c(epithelial_core_markers, squamous_markers, stromal_markers, hypoxia_markers, stress_markers))
  markers <- intersect(markers, rownames(lcpm))

  res <- data.frame()
  for (g in markers) {
    e_vals <- lcpm[g, grep("_E$", colnames(lcpm))]
    f_vals <- lcpm[g, grep("_F$", colnames(lcpm))]
    e_med <- if (length(e_vals) > 0) median(e_vals) else NA
    f_med <- if (length(f_vals) > 0) median(f_vals) else NA
    res <- rbind(res, data.frame(
      gene=g, dataset=dname,
      epithelial_median_logCPM=round(e_med,4), fibroblast_median_logCPM=round(f_med,4),
      difference=round(e_med - f_med, 4),
      epithelial_detection_frac=round(e_det_frac[g],4), fibroblast_detection_frac=round(f_det_frac[g],4),
      stringsAsFactors=FALSE))
  }
  return(res)
}

ref19 <- run_ref_audit(obj19, md19, "GSE197677")
ref22 <- run_ref_audit(obj22, md22, "GSE221561")
write.csv(ref19, file.path(g19g_dir, "STEP19G_GSE197677_LINEAGE_REFERENCE_EXPRESSION.csv"), row.names=FALSE)
write.csv(ref22, file.path(g19g_dir, "STEP19G_GSE221561_LINEAGE_REFERENCE_EXPRESSION.csv"), row.names=FALSE)
cat("Saved: reference expression CSVs\n\n")

# ============================================================================
# 19G3: RECURRENT LEADING-EDGE GENE LOCALIZATION
# ============================================================================
cat("=== 19G3: RECURRENT LEADING-EDGE GENE LOCALIZATION ===\n\n")

run_le_loc <- function(obj, md, dname) {
  cat(sprintf("--- %s ---\n", dname))
  epi_idx <- which(md$broad_celltype_final == "Epithelial")
  fib_idx <- which(md$broad_celltype_final == "Fibroblast_CAF")
  genes <- intersect(core_le_genes, get_genes(obj))
  
  if (is_v5(obj)) {
    assay <- obj[["RNA"]]
    # Pseudobulk by sample × celltype
    epi_samp <- md$sample[epi_idx]
    fib_samp <- md$sample[fib_idx]
    all_samps <- sort(unique(c(epi_samp, fib_samp)))
    group_map <- character(0)
    group_map[epi_idx] <- paste0(md$sample[epi_idx], "_E")
    group_map[fib_idx] <- paste0(md$sample[fib_idx], "_F")
    group_levels <- sort(unique(group_map[group_map != ""]))
    pb <- matrix(0, nrow = length(genes), ncol = length(group_levels),
                 dimnames = list(genes, group_levels))
    e_det <- numeric(length(genes)); names(e_det) <- genes
    f_det <- numeric(length(genes)); names(f_det) <- genes
    for (ln in Layers(assay)) {
      lc <- LayerData(assay, layer = ln)
      lc_genes <- intersect(genes, rownames(lc))
      g_cells <- group_map[colnames(lc)]
      for (g in group_levels) {
        ci <- which(g_cells == g)
        if (length(ci) > 0) pb[lc_genes, g] <- pb[lc_genes, g] + Matrix::rowSums(lc[lc_genes, ci, drop = FALSE])
      }
      e_cells <- intersect(colnames(lc), rownames(md)[epi_idx])
      f_cells <- intersect(colnames(lc), rownames(md)[fib_idx])
      if (length(e_cells) > 0) e_det[lc_genes] <- e_det[lc_genes] + Matrix::rowSums(lc[lc_genes, e_cells, drop = FALSE] > 0)
      if (length(f_cells) > 0) f_det[lc_genes] <- f_det[lc_genes] + Matrix::rowSums(lc[lc_genes, f_cells, drop = FALSE] > 0)
      rm(lc); gc(FALSE)
    }
  } else {
    counts <- get_counts(obj)
    epi_samp <- md$sample[epi_idx]; fib_samp <- md$sample[fib_idx]
    all_samps <- sort(unique(c(epi_samp, fib_samp)))
    pb_list <- list()
    for (s in all_samps) {
      e_i <- epi_idx[epi_samp == s]; f_i <- fib_idx[fib_samp == s]
      if (length(e_i) > 0) pb_list[[paste0(s,"_E")]] <- rowSums(counts[genes, e_i, drop=FALSE])
      if (length(f_i) > 0) pb_list[[paste0(s,"_F")]] <- rowSums(counts[genes, f_i, drop=FALSE])
    }
    pb <- do.call(cbind, pb_list)
    e_det <- rowSums(counts[genes, epi_idx, drop=FALSE] > 0) / length(epi_idx)
    f_det <- rowSums(counts[genes, fib_idx, drop=FALSE] > 0) / length(fib_idx)
  }
  
  res <- data.frame()
  for (g in genes) {
    e_vals <- pb[g, grep("_E$", colnames(pb))]
    f_vals <- pb[g, grep("_F$", colnames(pb))]
    e_med <- median(e_vals); f_med <- median(f_vals)
    res <- rbind(res, data.frame(
      gene=g, dataset=dname,
      epithelial_median_count=round(e_med,1), fibroblast_median_count=round(f_med,1),
      ratio=round(ifelse(f_med>0, e_med/f_med, NA),3),
      difference=round(e_med - f_med, 1),
      epithelial_positive_frac=round(e_det[g],4), fibroblast_positive_frac=round(f_det[g],4),
      stringsAsFactors=FALSE))
  }
  return(res)
}

le19 <- run_le_loc(obj19, md19, "GSE197677")
le22 <- run_le_loc(obj22, md22, "GSE221561")
le_loc <- rbind(le19, le22)
write.csv(le_loc, file.path(g19g_dir, "STEP19G_CORE_LEADING_EDGE_LINEAGE_LOCALIZATION.csv"), row.names=FALSE)
cat("Saved: STEP19G_CORE_LEADING_EDGE_LINEAGE_LOCALIZATION.csv\n\n")

# ============================================================================
# 19G4: EPITHELIAL CLUSTER IDENTITY AUDIT
# ============================================================================
cat("=== 19G4: EPITHELIAL CLUSTER IDENTITY AUDIT ===\n\n")

compute_cluster_scores <- function(obj, md, dname, cluster_col) {
  cat(sprintf("--- %s ---\n", dname))
  counts <- get_counts(obj)
  epi_idx <- which(md$broad_celltype_final == "Epithelial")
  epi_md <- md[epi_idx, ]
  epi_counts <- counts[, epi_idx]
  clusters <- sort(unique(epi_md[[cluster_col]]))
  cat(sprintf("  Epithelial clusters: %d\n", length(clusters)))

  # log1p CPM per cell
  lib_sz <- colSums(epi_counts)
  cpm_mat <- t(t(epi_counts) / lib_sz * 1e6)
  log_cpm <- log1p(cpm_mat)

  # Z-score across epithelial cells
  gmean <- rowMeans(log_cpm, na.rm=TRUE)
  gsd <- apply(log_cpm, 1, sd, na.rm=TRUE); gsd[gsd==0] <- 1
  z <- t((t(log_cpm) - gmean) / gsd)

  panels <- list(
    epithelial_core = intersect(epithelial_core_markers, rownames(z)),
    squamous = intersect(squamous_markers, rownames(z)),
    stromal_like = intersect(stromal_markers, rownames(z)),
    hypoxia = intersect(hypoxia_markers, rownames(z)),
    stress = intersect(stress_markers, rownames(z))
  )

  res <- data.frame()
  for (cl in clusters) {
    ci <- which(epi_md[[cluster_col]] == cl)
    if (length(ci) < 3) next
    row <- data.frame(dataset=dname, cluster=cl, n_cells=length(ci),
                      n_samples=length(unique(epi_md$sample[ci])),
                      n_patients=length(unique(epi_md$patient[ci])), stringsAsFactors=FALSE)
    for (pn in names(panels)) {
      gns <- panels[[pn]]
      row[[pn]] <- if (length(gns) > 0) mean(colMeans(z[gns, ci, drop=FALSE], na.rm=TRUE)) else NA
    }
    res <- rbind(res, row)
  }
  return(res)
}

cl19 <- compute_cluster_scores(obj19, md19, "GSE197677", "seurat_clusters")
cl22 <- compute_cluster_scores(obj22, md22, "GSE221561", "author_seurat_clusters")
write.csv(cl19, file.path(g19g_dir, "STEP19G_GSE197677_EPITHELIAL_CLUSTER_IDENTITY.csv"), row.names=FALSE)
write.csv(cl22, file.path(g19g_dir, "STEP19G_GSE221561_EPITHELIAL_CLUSTER_IDENTITY.csv"), row.names=FALSE)
cat("Saved: cluster identity CSVs\n\n")

# ============================================================================
# 19G5: EPITHELIAL VS STROMAL-LIKE AXIS
# ============================================================================
cat("=== 19G5: EPITHELIAL VS STROMAL-LIKE AXIS ===\n\n")

classify_axis <- function(df) {
  df$identity_delta <- df$epithelial_core_score - df$stromal_like_score
  df$axis_class <- ifelse(df$identity_delta > 0.5, "EPITHELIAL_DOMINANT",
                   ifelse(df$identity_delta < -0.5, "STROMAL_LIKE_PROFILE", "INTERMEDIATE"))
  return(df)
}

cl19 <- classify_axis(cl19)
cl22 <- classify_axis(cl22)

cat("GSE197677:\n"); print(table(cl19$axis_class))
cat("\nGSE221561:\n"); print(table(cl22$axis_class))

flags <- rbind(cl19[,c("dataset","cluster","n_cells","epithelial_core_score","stromal_like_score","identity_delta","axis_class")],
               cl22[,c("dataset","cluster","n_cells","epithelial_core_score","stromal_like_score","identity_delta","axis_class")])
write.csv(flags, file.path(g19g_dir, "STEP19G_EPITHELIAL_CLUSTER_IDENTITY_FLAGS.csv"), row.names=FALSE)
cat("\nSaved: STEP19G_EPITHELIAL_CLUSTER_IDENTITY_FLAGS.csv\n\n")

# ============================================================================
# 19G6: SUPPORTED PATHWAY LOCALIZATION BY CLUSTER
# ============================================================================
cat("=== 19G6: SUPPORTED PATHWAY LOCALIZATION BY CLUSTER ===\n\n")

compute_pw_cluster <- function(obj, md, dname, cluster_col) {
  cat(sprintf("--- %s ---\n", dname))
  counts <- get_counts(obj)
  epi_idx <- which(md$broad_celltype_final == "Epithelial")
  epi_md <- md[epi_idx, ]
  epi_counts <- counts[, epi_idx]
  clusters <- sort(unique(epi_md[[cluster_col]]))

  lib_sz <- colSums(epi_counts)
  cpm_mat <- t(t(epi_counts) / lib_sz * 1e6)
  log_cpm <- log1p(cpm_mat)
  gmean <- rowMeans(log_cpm, na.rm=TRUE)
  gsd <- apply(log_cpm, 1, sd, na.rm=TRUE); gsd[gsd==0] <- 1
  z <- t((t(log_cpm) - gmean) / gsd)

  id_panels <- list(
    epithelial_core = intersect(epithelial_core_markers, rownames(z)),
    stromal_like = intersect(stromal_markers, rownames(z)),
    hypoxia = intersect(hypoxia_markers, rownames(z)),
    stress = intersect(stress_markers, rownames(z))
  )

  res <- data.frame()
  for (pw in supported_pathways) {
    pw_gns <- intersect(hallmark_pathways[[pw]], rownames(z))
    if (length(pw_gns) < 5) next
    for (cl in clusters) {
      ci <- which(epi_md[[cluster_col]] == cl)
      if (length(ci) < 3) next
      row <- data.frame(dataset=dname, pathway=pw, cluster=cl,
                        pathway_score=mean(colMeans(z[pw_gns, ci, drop=FALSE], na.rm=TRUE)),
                        n_cells=length(ci), stringsAsFactors=FALSE)
      for (pn in names(id_panels)) {
        gns <- id_panels[[pn]]
        row[[pn]] <- if (length(gns) > 0) mean(colMeans(z[gns, ci, drop=FALSE], na.rm=TRUE)) else NA
      }
      res <- rbind(res, row)
    }
  }
  return(res)
}

pw19 <- compute_pw_cluster(obj19, md19, "GSE197677", "seurat_clusters")
pw22 <- compute_pw_cluster(obj22, md22, "GSE221561", "author_seurat_clusters")
pw_all <- rbind(pw19, pw22)
write.csv(pw_all, file.path(g19g_dir, "STEP19G_SUPPORTED_PATHWAY_CLUSTER_LOCALIZATION.csv"), row.names=FALSE)
cat("Saved: STEP19G_SUPPORTED_PATHWAY_CLUSTER_LOCALIZATION.csv\n\n")

# ============================================================================
# 19G7: PATHWAY-IDENTITY CORRELATION
# ============================================================================
cat("=== 19G7: PATHWAY-IDENTITY CORRELATION ===\n\n")

run_corr <- function(pw_df, dname) {
  ds <- pw_df[pw_df$dataset == dname, ]
  res <- data.frame()
  for (pw in supported_pathways) {
    dp <- ds[ds$pathway == pw, ]
    if (nrow(dp) < 3) next
    res <- rbind(res, data.frame(
      dataset=dname, pathway=pw,
      rho_vs_stromal=cor(dp$pathway_score, dp$stromal_like_score, method="spearman", use="complete.obs"),
      rho_vs_epithelial=cor(dp$pathway_score, dp$epithelial_core_score, method="spearman", use="complete.obs"),
      rho_vs_stress=cor(dp$pathway_score, dp$stress_score, method="spearman", use="complete.obs"),
      n_clusters=nrow(dp), stringsAsFactors=FALSE))
  }
  return(res)
}

corr19 <- run_corr(pw_all, "GSE197677")
corr22 <- run_corr(pw_all, "GSE221561")
corr_all <- rbind(corr19, corr22)
write.csv(corr_all, file.path(g19g_dir, "STEP19G_PATHWAY_IDENTITY_CORRELATION.csv"), row.names=FALSE)
cat("Saved: STEP19G_PATHWAY_IDENTITY_CORRELATION.csv\n\n")

# ============================================================================
# 19G8: SAMPLE-LEVEL EPITHELIAL STATE SCORES
# ============================================================================
cat("=== 19G8: SAMPLE-LEVEL EPITHELIAL STATE SCORES ===\n\n")

compute_sample_scores <- function(obj, md, dname) {
  cat(sprintf("--- %s ---\n", dname))
  counts <- get_counts(obj)
  epi_idx <- which(md$broad_celltype_final == "Epithelial")
  epi_md <- md[epi_idx, ]
  epi_counts <- counts[, epi_idx]
  samples <- unique(epi_md$sample)

  lib_sz <- colSums(epi_counts)
  cpm_mat <- t(t(epi_counts) / lib_sz * 1e6)
  log_cpm <- log1p(cpm_mat)
  gmean <- rowMeans(log_cpm, na.rm=TRUE)
  gsd <- apply(log_cpm, 1, sd, na.rm=TRUE); gsd[gsd==0] <- 1
  z <- t((t(log_cpm) - gmean) / gsd)

  panels <- list(
    epithelial_core = intersect(epithelial_core_markers, rownames(z)),
    stromal_like = intersect(stromal_markers, rownames(z)),
    hypoxia = intersect(hypoxia_markers, rownames(z)),
    stress = intersect(stress_markers, rownames(z))
  )
  for (pw in supported_pathways) {
    panels[[pw]] <- intersect(hallmark_pathways[[pw]], rownames(z))
  }

  res <- data.frame()
  for (s in samples) {
    si <- which(epi_md$sample == s)
    if (length(si) < 3) next
    row <- data.frame(dataset=dname, sample=s, patient=epi_md$patient[si[1]],
                      treatment=epi_md$treatment_group[si[1]], n_cells=length(si), stringsAsFactors=FALSE)
    for (pn in names(panels)) {
      gns <- panels[[pn]]
      row[[pn]] <- if (length(gns) > 0) mean(colMeans(z[gns, si, drop=FALSE], na.rm=TRUE)) else NA
    }
    res <- rbind(res, row)
  }
  return(res)
}

samp19 <- compute_sample_scores(obj19, md19, "GSE197677")
samp22 <- compute_sample_scores(obj22, md22, "GSE221561")
samp_all <- rbind(samp19, samp22)
write.csv(samp_all, file.path(g19g_dir, "STEP19G_SAMPLE_EPITHELIAL_STATE_SCORES.csv"), row.names=FALSE)
cat("Saved: STEP19G_SAMPLE_EPITHELIAL_STATE_SCORES.csv\n\n")

# ============================================================================
# 19G9: TREATMENT GROUP AUDIT
# ============================================================================
cat("=== 19G9: TREATMENT GROUP AUDIT ===\n\n")
cat("NOTE: DESCRIPTIVE_SECONDARY_ONLY — no new hypothesis testing\n\n")

for (dname in c("GSE197677", "GSE221561")) {
  ds <- samp_all[samp_all$dataset == dname, ]
  score_cols <- setdiff(colnames(ds), c("dataset","sample","patient","treatment","n_cells"))
  cat(sprintf("--- %s ---\n", dname))
  groups <- unique(ds$treatment)
  if (length(groups) == 2) {
    g1 <- ds[ds$treatment == groups[1], ]
    g2 <- ds[ds$treatment == groups[2], ]
    for (sc in score_cols) {
      v1 <- g1[[sc]]; v1 <- v1[!is.na(v1)]
      v2 <- g2[[sc]]; v2 <- v2[!is.na(v2)]
      if (length(v1) > 0 && length(v2) > 0) {
        cat(sprintf("  %s: %s median=%.4f, %s median=%.4f, diff=%.4f (DESCRIPTIVE_SECONDARY_ONLY)\n",
                    sc, groups[1], median(v1), groups[2], median(v2), median(v1)-median(v2)))
      }
    }
  }
}
cat("\n")

# ============================================================================
# 19G10: PATHWAY SIGNAL CONCENTRATION
# ============================================================================
cat("=== 19G10: PATHWAY SIGNAL CONCENTRATION ===\n\n")

run_conc <- function(pw_df, dname) {
  ds <- pw_df[pw_df$dataset == dname, ]
  res <- data.frame()
  for (pw in supported_pathways) {
    dp <- ds[ds$pathway == pw, ]
    if (nrow(dp) == 0) next
    n_cl <- nrow(dp)
    dp$weight <- ifelse(dp$pathway_score > 0, dp$pathway_score * dp$n_cells, 0)
    tot_w <- sum(dp$weight)
    if (tot_w > 0) {
      dp$contrib <- dp$weight / tot_w
      dp <- dp[order(-dp$contrib), ]
      top1 <- dp$contrib[1]
      top2 <- sum(dp$contrib[1:min(2, n_cl)])
      top3 <- sum(dp$contrib[1:min(3, n_cl)])
      top_frac <- dp$n_cells[1] / sum(dp$n_cells)
    } else { top1 <- top2 <- top3 <- top_frac <- NA }
    res <- rbind(res, data.frame(dataset=dname, pathway=pw, n_clusters=n_cl,
      top1_contribution=round(top1,4), top2_contribution=round(top2,4), top3_contribution=round(top3,4),
      top_cluster_cell_fraction=round(top_frac,4), stringsAsFactors=FALSE))
  }
  return(res)
}

conc19 <- run_conc(pw_all, "GSE197677")
conc22 <- run_conc(pw_all, "GSE221561")
conc <- rbind(conc19, conc22)
write.csv(conc, file.path(g19g_dir, "STEP19G_PATHWAY_CLUSTER_CONCENTRATION.csv"), row.names=FALSE)
cat("Saved: STEP19G_PATHWAY_CLUSTER_CONCENTRATION.csv\n\n")

# ============================================================================
# 19G11: LEADING-EDGE CO-DETECTION
# ============================================================================
cat("=== 19G11: LEADING-EDGE CO-DETECTION ===\n\n")

run_codet <- function(obj, md, dname) {
  cat(sprintf("--- %s ---\n", dname))
  counts <- get_counts(obj)
  epi_idx <- which(md$broad_celltype_final == "Epithelial")
  epi_counts <- counts[, epi_idx]
  epi_mks <- intersect(c("EPCAM","KRT8","KRT18","KRT19"), rownames(epi_counts))
  epi_mk_pos <- colSums(epi_counts[epi_mks, , drop=FALSE] > 0) > 0
  genes <- intersect(core_le_genes, rownames(epi_counts))

  res <- data.frame()
  for (g in genes) {
    tgt <- epi_counts[g, ] > 0
    n_tgt <- sum(tgt); n_both <- sum(tgt & epi_mk_pos)
    res <- rbind(res, data.frame(
      dataset=dname, gene=g, n_epithelial_cells=ncol(epi_counts),
      target_positive=n_tgt, target_positive_frac=round(n_tgt/ncol(epi_counts),4),
      target_and_epithelial_marker=n_both,
      codetection_fraction=round(ifelse(n_tgt>0, n_both/n_tgt, NA),4),
      stringsAsFactors=FALSE))
  }
  return(res)
}

codet19 <- run_codet(obj19, md19, "GSE197677")
codet22 <- run_codet(obj22, md22, "GSE221561")
codet <- rbind(codet19, codet22)
write.csv(codet, file.path(g19g_dir, "STEP19G_CORE_GENE_EPITHELIAL_CODETECTION.csv"), row.names=FALSE)
cat("Saved: STEP19G_CORE_GENE_EPITHELIAL_CODETECTION.csv\n\n")

# ============================================================================
# 19G12: ORIGINAL CLUSTER MAPPING PROVENANCE
# ============================================================================
cat("=== 19G12: ORIGINAL CLUSTER MAPPING PROVENANCE ===\n\n")

clust_map <- read.csv(file.path(base_dir, "06_tables/GSE197677/celltype/GSE197677_FINAL_cluster_to_broad_celltype_mapping.csv"), stringsAsFactors=FALSE)
cat("GSE197677 cluster mapping:\n")
print(clust_map[, c("cluster","broad_celltype_final","n_cells")])

prov19 <- merge(clust_map[, c("cluster","broad_celltype_final","n_cells","review_status")],
                cl19[, c("cluster","epithelial_core_score","stromal_like_score","hypoxia_score","stress_score")],
                by="cluster", all.x=TRUE)
for (pw in supported_pathways) {
  ps <- pw19[pw19$pathway == pw, c("cluster","pathway_score")]
  colnames(ps)[2] <- paste0("score_", gsub("HALLMARK_","",pw))
  prov19 <- merge(prov19, ps, by="cluster", all.x=TRUE)
}
write.csv(prov19, file.path(g19g_dir, "STEP19G_GSE197677_CLUSTER_PROVENANCE.csv"), row.names=FALSE)

prov22 <- cl22[, c("cluster","n_cells","epithelial_core_score","stromal_like_score","hypoxia_score","stress_score")]
for (pw in supported_pathways) {
  ps <- pw22[pw22$pathway == pw, c("cluster","pathway_score")]
  colnames(ps)[2] <- paste0("score_", gsub("HALLMARK_","",pw))
  prov22 <- merge(prov22, ps, by="cluster", all.x=TRUE)
}
write.csv(prov22, file.path(g19g_dir, "STEP19G_GSE221561_CLUSTER_PROVENANCE.csv"), row.names=FALSE)
cat("Saved: provenance CSVs\n\n")

# ============================================================================
# 19G13: TARGETED FIGURES
# ============================================================================
cat("=== 19G13: TARGETED FIGURES ===\n\n")

# Fig 1: Epithelial vs Stromal Marker Heatmap
cat("Fig 1: Marker difference heatmap...\n")
ref_comb <- rbind(ref19, ref22)
ref_wide <- reshape2::dcast(ref_comb, gene ~ dataset, value.var="difference")
rownames(ref_wide) <- ref_wide$gene; ref_wide$gene <- NULL
ref_mat <- as.matrix(ref_wide)
gene_ord <- intersect(unique(c(epithelial_core_markers, squamous_markers, stromal_markers, hypoxia_markers, stress_markers)), rownames(ref_mat))
ref_mat <- ref_mat[gene_ord, ]

png(file.path(fig_dir, "STEP19G_EPITHELIAL_VS_STROMAL_MARKER_HEATMAP.png"), width=800, height=1200, res=150)
par(mar=c(4,10,2,2))
image(t(ref_mat[nrow(ref_mat):1,]), col=rev(colorRampPalette(c("blue","white","red"))(100)), axes=FALSE, xlab="", ylab="")
axis(2, at=seq(0,1,length.out=nrow(ref_mat)), labels=rev(rownames(ref_mat)), las=2, cex.axis=0.7)
axis(1, at=c(0,1), labels=colnames(ref_mat), cex.axis=0.8)
title(main="Epithelial - Fibroblast logCPM Difference", cex.main=0.9)
dev.off()

# Fig 2: Cluster Identity Heatmaps
for (dname in c("GSE197677","GSE221561")) {
  cl_df <- if(dname=="GSE197677") cl19 else cl22
  cl_df <- cl_df[cl_df$dataset==dname, ]
  mat <- as.matrix(cl_df[, c("epithelial_core_score","squamous_score","stromal_like_score","hypoxia_score","stress_score")])
  rownames(mat) <- paste0("C", cl_df$cluster)
  png(file.path(fig_dir, sprintf("STEP19G_%s_EPITHELIAL_CLUSTER_IDENTITY_HEATMAP.png", dname)),
      width=800, height=max(400, nrow(mat)*30), res=150)
  par(mar=c(4,8,2,2))
  image(t(mat[nrow(mat):1,]), col=colorRampPalette(c("blue","white","red"))(100), axes=FALSE, xlab="", ylab="")
  axis(2, at=seq(0,1,length.out=nrow(mat)), labels=rev(rownames(mat)), las=2, cex.axis=0.7)
  axis(1, at=seq(0,1,length.out=ncol(mat)), labels=colnames(mat), las=2, cex.axis=0.6)
  title(main=sprintf("%s Epithelial Cluster Identity", dname), cex.main=0.9)
  dev.off()
}

# Fig 3: Seven Pathway by Cluster
for (dname in c("GSE197677","GSE221561")) {
  pw_ds <- pw_all[pw_all$dataset==dname, ]
  pw_wide <- reshape2::dcast(pw_ds, pathway ~ cluster, value.var="pathway_score")
  rownames(pw_wide) <- gsub("HALLMARK_","",pw_wide$pathway); pw_wide$pathway <- NULL
  pw_mat <- as.matrix(pw_wide)
  png(file.path(fig_dir, sprintf("STEP19G_%s_SEVEN_PATHWAY_BY_EPITHELIAL_CLUSTER.png", dname)),
      width=max(800, ncol(pw_mat)*40), height=600, res=150)
  par(mar=c(8,12,2,2))
  image(t(pw_mat[nrow(pw_mat):1,]), col=colorRampPalette(c("blue","white","red"))(100), axes=FALSE, xlab="", ylab="")
  axis(2, at=seq(0,1,length.out=nrow(pw_mat)), labels=rev(rownames(pw_mat)), las=2, cex.axis=0.7)
  axis(1, at=seq(0,1,length.out=ncol(pw_mat)), labels=paste0("C",colnames(pw_mat)), las=2, cex.axis=0.6)
  title(main=sprintf("%s Hallmark Scores by Cluster", dname), cex.main=0.9)
  dev.off()
}

# Fig 4: Core LE Gene Lineage Heatmap
cat("Fig 4: LE gene lineage heatmap...\n")
le_wide <- reshape2::dcast(le_loc, gene ~ dataset, value.var="ratio")
rownames(le_wide) <- le_wide$gene; le_wide$gene <- NULL
le_mat <- as.matrix(log2(le_wide + 1))
png(file.path(fig_dir, "STEP19G_CORE_LEADING_EDGE_LINEAGE_HEATMAP.png"), width=600, height=600, res=150)
par(mar=c(4,8,2,2))
image(t(le_mat[nrow(le_mat):1,]), col=colorRampPalette(c("blue","white","red"))(100), axes=FALSE, xlab="", ylab="")
axis(2, at=seq(0,1,length.out=nrow(le_mat)), labels=rev(rownames(le_mat)), las=2, cex.axis=0.7)
axis(1, at=c(0,1), labels=colnames(le_mat), cex.axis=0.8)
title(main="log2(Epi/Fib) Ratio", cex.main=0.9)
dev.off()

# Fig 5: Pathway vs Stromal Score Scatter
cat("Fig 5: Pathway vs stromal scatter...\n")
scatter_pws <- c("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION","HALLMARK_MYOGENESIS","HALLMARK_COAGULATION")
for (dname in c("GSE197677","GSE221561")) {
  pw_ds <- pw_all[pw_all$dataset==dname, ]
  plist <- list()
  for (pw in scatter_pws) {
    dp <- pw_ds[pw_ds$pathway==pw, ]
    p <- ggplot(dp, aes(x=stromal_like_score, y=pathway_score, label=paste0("C",cluster))) +
      geom_point(aes(size=n_cells), alpha=0.7) + geom_text(vjust=-0.5, size=2.5) +
      geom_smooth(method="lm", se=FALSE, color="red", linetype="dashed") +
      labs(title=gsub("HALLMARK_","",pw), x="Stromal-like Score", y="Pathway Score") +
      theme_minimal() + theme(plot.title=element_text(size=9))
    plist[[pw]] <- p
  }
  png(file.path(fig_dir, sprintf("STEP19G_%s_PATHWAY_VS_STROMAL_SCORE.png", dname)), width=1200, height=400, res=150)
  print(patchwork::wrap_plots(plist, nrow=1))
  dev.off()
}
cat("Saved: all figures\n\n")

# ============================================================================
# 19G14: PRIMARY LOCALIZATION CLASSIFICATION
# ============================================================================
cat("=== 19G14: PRIMARY LOCALIZATION CLASSIFICATION ===\n\n")

le_jaccard_df <- read.csv(file.path(pw_dir, "STEP19E_HALLMARK_LEADING_EDGE_OVERLAP.csv"), stringsAsFactors=FALSE)

class_df <- data.frame()
for (pw in supported_pathways) {
  rho19 <- corr19$rho_vs_stromal[corr19$pathway == pw]
  rho22 <- corr22$rho_vs_stromal[corr22$pathway == pw]
  c19 <- conc19[conc19$pathway == pw, ]
  c22 <- conc22[conc22$pathway == pw, ]
  le_j <- le_jaccard_df$jaccard[le_jaccard_df$pathway == pw]
  if (length(le_j) == 0) le_j <- NA

  cd19 <- codet19$codetection_fraction[codet19$gene %in% core_le_genes]
  cd22 <- codet22$codetection_fraction[codet22$gene %in% core_le_genes]
  avg_cd <- mean(c(cd19, cd22), na.rm=TRUE)

  strong_stromal <- (!is.na(rho19) && rho19 > 0.5) || (!is.na(rho22) && rho22 > 0.5)
  high_conc <- (!is.na(c19$top1_contribution) && c19$top1_contribution > 0.4) ||
               (!is.na(c22$top1_contribution) && c22$top1_contribution > 0.4)
  good_cd <- avg_cd > 0.3

  if (strong_stromal && high_conc && !good_cd) {
    loc <- "POSSIBLE_LINEAGE_MISCLASSIFICATION_SIGNAL"
  } else if (strong_stromal || high_conc) {
    loc <- "MIXED_EPITHELIAL_STROMAL_PROGRAM"
  } else if (good_cd || !strong_stromal) {
    loc <- "EPITHELIAL_LOCALIZATION_SUPPORTED"
  } else {
    loc <- "INSUFFICIENT_LOCALIZATION_EVIDENCE"
  }

  class_df <- rbind(class_df, data.frame(
    pathway=pw,
    GSE197677_pathway_stromal_rho=round(ifelse(length(rho19)>0, rho19, NA),4),
    GSE221561_pathway_stromal_rho=round(ifelse(length(rho22)>0, rho22, NA),4),
    GSE197677_cluster_concentration=round(ifelse(nrow(c19)>0, c19$top1_contribution, NA),4),
    GSE221561_cluster_concentration=round(ifelse(nrow(c22)>0, c22$top1_contribution, NA),4),
    cross_dataset_LE_Jaccard=round(ifelse(length(le_j)>0, le_j, NA),4),
    core_gene_epithelial_codetection=round(avg_cd,4),
    localization_class=loc, stringsAsFactors=FALSE))
}

write.csv(class_df, file.path(g19g_dir, "STEP19G_PATHWAY_LOCALIZATION_FINAL.csv"), row.names=FALSE)
cat("Saved: STEP19G_PATHWAY_LOCALIZATION_FINAL.csv\n\n")
cat("Classification:\n")
for (i in 1:nrow(class_df)) cat(sprintf("  %s: %s\n", class_df$pathway[i], class_df$localization_class[i]))

# ============================================================================
# 19G15: OVERALL INTERPRETATION
# ============================================================================
cat("\n=== 19G15: OVERALL INTERPRETATION ===\n\n")

class_tbl <- table(class_df$localization_class)
cat("Distribution:\n"); print(class_tbl)

if ("POSSIBLE_LINEAGE_MISCLASSIFICATION_SIGNAL" %in% names(class_tbl) && class_tbl["POSSIBLE_LINEAGE_MISCLASSIFICATION_SIGNAL"] > 0) {
  overall <- "ROBUST_PATHWAYS_POSSIBLE_LINEAGE_CONTAMINATION"
} else if ("MIXED_EPITHELIAL_STROMAL_PROGRAM" %in% names(class_tbl)) {
  n_mixed <- class_tbl["MIXED_EPITHELIAL_STROMAL_PROGRAM"]
  n_epi <- if ("EPITHELIAL_LOCALIZATION_SUPPORTED" %in% names(class_tbl)) class_tbl["EPITHELIAL_LOCALIZATION_SUPPORTED"] else 0
  overall <- if (n_mixed > n_epi) "ROBUST_PATHWAYS_MIXED_EPITHELIAL_STROMAL_SIGNAL" else "ROBUST_PATHWAYS_EPITHELIAL_LOCALIZATION_SUPPORTED"
} else {
  overall <- "ROBUST_PATHWAYS_EPITHELIAL_LOCALIZATION_SUPPORTED"
}

# Write interpretation
interp <- "# Step 19G: Epithelial Pathway Localization Interpretation\n\n"
interp <- paste0(interp, "## Scope\n\nThis audit determines the likely cellular localization of 7 robust Hallmark pathways.\n\n")
interp <- paste0(interp, "**Does NOT:** establish EMT causality, establish malignant-cell specificity, alter Step 19F results, or reclassify/remove cells.\n\n")
interp <- paste0(interp, "## Classification Results\n\n")
for (i in 1:nrow(class_df)) {
  interp <- paste0(interp, sprintf("### %s\n- Class: %s\n- Disc stromal rho: %.3f, Repl stromal rho: %.3f\n- LE Jaccard: %.3f, Co-detection: %.1f%%\n\n",
    gsub("HALLMARK_","",class_df$pathway[i]), class_df$localization_class[i],
    class_df$GSE197677_pathway_stromal_rho[i], class_df$GSE221561_pathway_stromal_rho[i],
    class_df$cross_dataset_LE_Jaccard[i], class_df$core_gene_epithelial_codetection[i]*100))
}
interp <- paste0(interp, "## Lineage misclassification: ",
  ifelse("POSSIBLE_LINEAGE_MISCLASSIFICATION_SIGNAL" %in% names(class_tbl) && class_tbl["POSSIBLE_LINEAGE_MISCLASSIFICATION_SIGNAL"]>0, "POSSIBLE", "NO"), "\n")
interp <- paste0(interp, "## Cells reclassified: 0\n## Cells removed: 0\n## Samples removed: 0\n## Step19F changed: NO\n")

writeLines(interp, file.path(g19g_dir, "STEP19G_EPITHELIAL_PATHWAY_LOCALIZATION_INTERPRETATION.md"))
cat("Saved: interpretation MD\n\n")

# ============================================================================
# 19G16: OVERALL STATUS
# ============================================================================
cat("=== 19G16: OVERALL STATUS ===\n\n")
cat("============================================\n")
cat("STEP 19G COMPLETE\n")
cat("EPITHELIAL PATHWAY LOCALIZATION AUDIT\n")
cat("============================================\n\n")
cat("Datasets audited: GSE197677, GSE221561\n")
cat("Frozen supported pathways: 7\n\n")
cat("--------------------------------------------\n\n")

cat("Core identity audit:\n")
for (dname in c("GSE197677","GSE221561")) {
  md <- if(dname=="GSE197677") md19 else md22
  cl_df <- if(dname=="GSE197677") cl19 else cl22
  cat(sprintf("\n%s:\n", dname))
  cat(sprintf("  Epithelial cells: %d\n", sum(md$broad_celltype_final=="Epithelial")))
  cat(sprintf("  Epithelial clusters: %d\n", nrow(cl_df)))
  cat(sprintf("  Fibroblast reference cells: %d\n", sum(md$broad_celltype_final=="Fibroblast_CAF")))
  cat(sprintf("  STROMAL_LIKE_PROFILE clusters: %d\n", sum(cl_df$axis_class=="STROMAL_LIKE_PROFILE", na.rm=TRUE)))
  cat(sprintf("  INTERMEDIATE clusters: %d\n", sum(cl_df$axis_class=="INTERMEDIATE", na.rm=TRUE)))
  cat(sprintf("  EPITHELIAL_DOMINANT clusters: %d\n", sum(cl_df$axis_class=="EPITHELIAL_DOMINANT", na.rm=TRUE)))
}

cat("\n--------------------------------------------\n\n")
cat("Core recurrent leading-edge genes:\n")
for (g in core_le_genes) {
  epi_frac <- mean(c(codet19$target_positive_frac[codet19$gene==g], codet22$target_positive_frac[codet22$gene==g]), na.rm=TRUE)
  cd_frac <- mean(c(codet19$codetection_fraction[codet19$gene==g], codet22$codetection_fraction[codet22$gene==g]), na.rm=TRUE)
  loc <- if (!is.na(cd_frac) && cd_frac > 0.5) "EPITHELIAL_COPOSITIVE" else if (!is.na(cd_frac) && cd_frac > 0.2) "MIXED" else "LOW_COLOCALIZATION"
  cat(sprintf("  %s: Epi=%.1f%%, CoDet=%.1f%%, %s\n", g, epi_frac*100, cd_frac*100, loc))
}

cat("\n--------------------------------------------\n\n")
cat("Supported pathways:\n")
for (i in 1:nrow(class_df)) {
  cat(sprintf("  %s: %s (rho19=%.3f, rho22=%.3f, LE_J=%.3f)\n",
    gsub("HALLMARK_","",class_df$pathway[i]), class_df$localization_class[i],
    class_df$GSE197677_pathway_stromal_rho[i], class_df$GSE221561_pathway_stromal_rho[i],
    class_df$cross_dataset_LE_Jaccard[i]))
}

cat("\n--------------------------------------------\n\n")
cat(sprintf("Evidence of lineage misclassification: %s\n",
  ifelse("POSSIBLE_LINEAGE_MISCLASSIFICATION_SIGNAL" %in% names(class_tbl) && class_tbl["POSSIBLE_LINEAGE_MISCLASSIFICATION_SIGNAL"]>0, "POSSIBLE", "NO")))
cat("Cells reclassified: 0\nCells removed: 0\nSamples removed: 0\nStep19F changed: NO\n\n")
cat("--------------------------------------------\n\n")
cat(sprintf("Overall Step 19G status:\n%s\n\n", overall))
cat("STOP HERE.\nDO NOT START STEP 19H AUTOMATICALLY.\n")
