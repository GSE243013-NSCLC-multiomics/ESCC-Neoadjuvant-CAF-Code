#!/usr/bin/env Rscript
# STEP 19G: MEMORY-SAFE EPITHELIAL LOCALIZATION AUDIT
# NO JoinLayers. NO full dense matrices. Layer-wise pseudobulk only.

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(edgeR)
  library(Matrix)
  library(ggplot2)
  library(reshape2)
})

base_dir <- "~/Downloads/ESCC Neoadjuvant"
pw_dir <- file.path(base_dir, "06_tables/cross_dataset/pathway_analysis")
g19g_dir <- file.path(pw_dir, "Step19G")
ckpt_dir <- file.path(g19g_dir, "checkpoints")
fig_dir <- file.path(base_dir, "05_figures/cross_dataset/pathway_analysis/Step19G")
dir.create(g19g_dir, recursive=TRUE, showWarnings=FALSE)
dir.create(ckpt_dir, recursive=TRUE, showWarnings=FALSE)
dir.create(fig_dir, recursive=TRUE, showWarnings=FALSE)

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

# ============================================================================
# MEMORY-SAFE HELPERS
# ============================================================================

aggregate_pseudobulk <- function(counts_mat, group_ids) {
  groups <- sort(unique(group_ids))
  if (length(groups) == 0) return(NULL)
  if (length(groups) == 1) {
    # Single group: just sum columns
    pb <- Matrix(rowSums(counts_mat), ncol=1, sparse=TRUE)
    colnames(pb) <- groups
    rownames(pb) <- rownames(counts_mat)
    return(pb)
  }
  group_factor <- factor(group_ids, levels=groups)
  design <- Matrix::sparse.model.matrix(~0 + group_factor)
  colnames(design) <- groups
  pb <- counts_mat %*% design
  pb <- as(pb, "dgCMatrix")
  return(pb)
}

get_count_layers <- function(obj) {
  rna <- obj@assays[["RNA"]]
  layers <- Layers(rna)
  count_layers <- layers[grepl("^counts($|\\.)", layers)]
  return(count_layers)
}

pseudobulk_layerwise <- function(obj, md, group_field, target_genes=NULL) {
  is_v5 <- inherits(obj@assays[["RNA"]], "Assay5")
  if (!is_v5) {
    counts <- GetAssayData(obj, assay="RNA", layer="counts")
    if (!is.null(target_genes)) {
      shared <- intersect(target_genes, rownames(counts))
      counts <- counts[shared, , drop=FALSE]
    }
    groups <- as.character(md[[group_field]])
    names(groups) <- rownames(md)
    groups <- groups[colnames(counts)]
    pb <- aggregate_pseudobulk(counts, groups)
    return(list(pb=pb, genes=rownames(pb)))
  }
  layer_names <- get_count_layers(obj)
  cat(sprintf("  Processing %d count layers...\n", length(layer_names)))

  # Collect per-group sums across layers
  all_group_sums <- list()  # group -> named numeric vector of gene sums
  all_genes_set <- character()

  for (ln in layer_names) {
    mat <- LayerData(obj@assays[["RNA"]], layer=ln)
    cell_ids <- colnames(mat)
    cell_md <- md[cell_ids, , drop=FALSE]
    groups <- as.character(cell_md[[group_field]])
    names(groups) <- cell_ids

    # Restrict to target genes if provided
    if (!is.null(target_genes)) {
      shared <- intersect(target_genes, rownames(mat))
      if (length(shared) == 0) { rm(mat); gc(); next }
      mat <- mat[shared, , drop=FALSE]
    }

    # Aggregate: for each group, sum columns
    layer_groups <- sort(unique(groups))
    for (g in layer_groups) {
      cells_g <- names(groups)[groups == g]
      if (length(cells_g) == 0) next
      sums_g <- as.numeric(rowSums(mat[, cells_g, drop=FALSE]))
      names(sums_g) <- rownames(mat)
      if (is.null(all_group_sums[[g]])) {
        all_group_sums[[g]] <- sums_g
      } else {
        # Align and sum
        shared_genes <- intersect(names(all_group_sums[[g]]), names(sums_g))
        if (length(shared_genes) > 0) {
          all_group_sums[[g]][shared_genes] <- all_group_sums[[g]][shared_genes] + sums_g[shared_genes]
        }
        new_genes <- setdiff(names(sums_g), names(all_group_sums[[g]]))
        if (length(new_genes) > 0) {
          all_group_sums[[g]] <- c(all_group_sums[[g]], sums_g[new_genes])
        }
      }
      all_genes_set <- union(all_genes_set, names(sums_g))
    }
    rm(mat); gc()
  }

  # Build final pseudobulk matrix
  groups <- sort(names(all_group_sums))
  pb <- Matrix(0, nrow=length(all_genes_set), ncol=length(groups), sparse=TRUE)
  rownames(pb) <- all_genes_set
  colnames(pb) <- groups
  for (g in groups) {
    genes_g <- names(all_group_sums[[g]])
    pb[genes_g, g] <- all_group_sums[[g]]
  }
  # Remove zero rows
  keep <- rowSums(pb) > 0
  pb <- pb[keep, , drop=FALSE]
  return(list(pb=pb, genes=rownames(pb)))
}

normalize_pseudobulk <- function(pb_mat) {
  # Convert to dense for edgeR (pseudobulk is small)
  pb_dense <- as.matrix(pb_mat)
  keep <- rowSums(pb_dense > 0) >= 2
  pb_dense <- pb_dense[keep, , drop=FALSE]
  dge <- DGEList(counts=pb_dense)
  dge <- calcNormFactors(dge, method="TMM")
  lcpm <- cpm(dge, log=TRUE, prior.count=2)
  return(lcpm)
}

# ============================================================================
# 19G0: INPUT + METADATA PREFLIGHT
# ============================================================================
cat("=== 19G0: INPUT + METADATA PREFLIGHT ===\n\n")

cat("Loading GSE197677...\n")
obj19 <- readRDS(file.path(base_dir, "03_objects/GSE197677/GSE197677_integrated_broadcelltype_annotated.rds"))
md19 <- obj19@meta.data
cat(sprintf("  %d cells, %d metadata fields\n", nrow(md19), ncol(md19)))
genes19 <- rownames(obj19@assays$RNA)

cat("Loading GSE221561...\n")
obj22 <- readRDS(file.path(base_dir, "03_objects/GSE221561/GSE221561_author_annotated.rds"))
md22 <- obj22@meta.data
cat(sprintf("  %d cells, %d metadata fields\n", nrow(md22), ncol(md22)))
cat("  NO JoinLayers called.\n")
layers22 <- get_count_layers(obj22)
cat(sprintf("  Count layers: %d\n", length(layers22)))
cat("  Layer names:", paste(layers22, collapse=", "), "\n")
for (ln in layers22) {
  mat <- LayerData(obj22[["RNA"]], layer=ln)
  cat(sprintf("    %s: %d genes x %d cells (%s)\n", ln, nrow(mat), ncol(mat), class(mat)[1]))
  rm(mat); gc()
}
genes22 <- rownames(obj22@assays$RNA)
cat(sprintf("  Total features: %d\n\n", length(genes22)))

cat("GSE197677 broad_celltype_final:\n")
print(table(md19$broad_celltype_final))
cat("\nGSE221561 author_Cell_type:\n")
print(table(md22$author_Cell_type))

mapping22 <- c("B cell"="B_plasma","Endothelial"="Endothelial","Epithelial"="Epithelial",
               "Fibroblast"="Fibroblast_CAF","Mast cell"="Mast","Myeloid"="Myeloid","T cell"="T_NK")
md22$broad_celltype_final <- mapping22[md22$author_Cell_type]
treat22 <- md22$author_Neoadjuvant
treat22[treat22 %in% c("Chemotherapy","Chemoradiotherapy","Chemoradiotherapy and Immunotherapy","Antiangiogenesis")] <- "Neoadjuvant_treated"
treat22[treat22 == "Surgery alone"] <- "Surgery_alone"
treat22[treat22 == "Adjacent normal"] <- "Adjacent_normal"
md22$treatment_group <- treat22

meta_pseudo <- read.csv(file.path(base_dir, "06_tables/cross_dataset/epithelial_pseudobulk/counts/GSE197677_epithelial_pseudobulk_sample_metadata.csv"), stringsAsFactors=FALSE)
md19$treatment_group <- setNames(meta_pseudo$treatment_standardized, meta_pseudo$sample)[md19$sample]

# Unify sample column name
md22$sample <- md22$sample_id
md22$patient <- md22$sample_id  # No separate patient column in GSE221561

cat("\nGSE221561 treatment:\n"); print(table(md22$treatment_group))
cat("GSE197677 treatment:\n"); print(table(md19$treatment_group))

cat("\nCell inventories:\n")
cat(sprintf("  GSE197677: %d total, %d epithelial, %d fibroblast\n",
    nrow(md19), sum(md19$broad_celltype_final=="Epithelial"), sum(md19$broad_celltype_final=="Fibroblast_CAF")))
cat(sprintf("  GSE221561: %d total, %d epithelial, %d fibroblast\n",
    nrow(md22), sum(md22$broad_celltype_final=="Epithelial"), sum(md22$broad_celltype_final=="Fibroblast_CAF")))

preflight <- data.frame(
  dataset=c("GSE197677","GSE221561"), cells=c(nrow(md19),nrow(md22)),
  epithelial_cells=c(sum(md19$broad_celltype_final=="Epithelial"),sum(md22$broad_celltype_final=="Epithelial")),
  fibroblast_cells=c(sum(md19$broad_celltype_final=="Fibroblast_CAF"),sum(md22$broad_celltype_final=="Fibroblast_CAF")),
  stringsAsFactors=FALSE)
write.csv(preflight, file.path(g19g_dir, "STEP19G_METADATA_PREFLIGHT.csv"), row.names=FALSE)
cat("\nSaved: STEP19G_METADATA_PREFLIGHT.csv\n\n")

# ============================================================================
# 19G1: FREEZE GENE PANELS
# ============================================================================
cat("=== 19G1: FREEZE TARGET GENE PANELS ===\n\n")
all_panel_genes <- unique(c(epithelial_core_markers, squamous_markers, stromal_markers, hypoxia_markers, stress_markers, core_le_genes))
avail19 <- all_panel_genes %in% genes19
avail22 <- all_panel_genes %in% genes22
cat(sprintf("GSE197677: %d/%d panel genes\n", sum(avail19), length(all_panel_genes)))
cat(sprintf("GSE221561: %d/%d panel genes\n", sum(avail22), length(all_panel_genes)))

gene_panel_map <- data.frame(gene=all_panel_genes, panel="unknown", stringsAsFactors=FALSE)
for (i in seq_along(all_panel_genes)) {
  g <- all_panel_genes[i]
  if (g %in% epithelial_core_markers) gene_panel_map$panel[i] <- "epithelial_core"
  else if (g %in% squamous_markers) gene_panel_map$panel[i] <- "squamous"
  else if (g %in% stromal_markers) gene_panel_map$panel[i] <- "stromal"
  else if (g %in% hypoxia_markers) gene_panel_map$panel[i] <- "hypoxia"
  else if (g %in% stress_markers) gene_panel_map$panel[i] <- "stress"
  else if (g %in% core_le_genes) gene_panel_map$panel[i] <- "core_le"
}
panel_df <- data.frame(gene=gene_panel_map$gene, panel=gene_panel_map$panel,
  available_GSE197677=avail19, available_GSE221561=avail22, stringsAsFactors=FALSE)
write.csv(panel_df, file.path(g19g_dir, "STEP19G_IDENTITY_MARKER_PANEL_AVAILABILITY.csv"), row.names=FALSE)
cat("Saved: STEP19G_IDENTITY_MARKER_PANEL_AVAILABILITY.csv\n\n")

# ============================================================================
# 19G2-3: LINEAGE REFERENCE + LE GENE LOCALIZATION (MEMORY-SAFE)
# ============================================================================
cat("=== 19G2-3: LINEAGE REFERENCE + LE LOCALIZATION ===\n\n")

run_lineage_audit <- function(obj, md, dname) {
  cat(sprintf("--- %s ---\n", dname))
  epi_idx <- which(md$broad_celltype_final == "Epithelial")
  fib_idx <- which(md$broad_celltype_final == "Fibroblast_CAF")
  md$group_id <- paste0(md$sample, "_", md$broad_celltype_final)

  # Only need marker genes + LE genes for this section
  needed_genes <- unique(c(epithelial_core_markers, squamous_markers, stromal_markers,
                           hypoxia_markers, stress_markers, core_le_genes))
  needed_genes <- intersect(needed_genes, rownames(obj@assays$RNA))
  cat(sprintf("  Target genes available: %d\n", length(needed_genes)))

  # Pseudobulk by sample x celltype
  pb <- pseudobulk_layerwise(obj, md, "group_id", target_genes=needed_genes)
  cat(sprintf("  Pseudobulk: %d genes x %d groups\n", nrow(pb$pb), ncol(pb$pb)))

  # Normalize
  lcpm <- normalize_pseudobulk(pb$pb)

  # Reference expression (19G2)
  ref_res <- data.frame()
  markers <- unique(c(epithelial_core_markers, squamous_markers, stromal_markers, hypoxia_markers, stress_markers))
  markers <- intersect(markers, rownames(lcpm))
  for (g in markers) {
    e_vals <- lcpm[g, grep("_Epithelial$", colnames(lcpm))]
    f_vals <- lcpm[g, grep("_Fibroblast_CAF$", colnames(lcpm))]
    e_med <- if (length(e_vals) > 0) median(e_vals) else NA
    f_med <- if (length(f_vals) > 0) median(f_vals) else NA
    # Detection from raw counts (cell-level, but sparse and selected genes only)
    e_det_sum <- 0; e_total <- 0; f_det_sum <- 0; f_total <- 0
    if (!inherits(obj@assays[["RNA"]], "Assay5")) {
      # v3/v4: get counts directly from the assay slot
      counts <- obj@assays[["RNA"]]@counts
      if (g %in% rownames(counts)) {
        e_det_sum <- sum(counts[g, epi_idx] > 0)
        e_total <- length(epi_idx)
        f_det_sum <- sum(counts[g, fib_idx] > 0)
        f_total <- length(fib_idx)
      }
    } else {
      for (ln in get_count_layers(obj)) {
        mat <- LayerData(obj@assays[["RNA"]], layer=ln)
        if (!(g %in% rownames(mat))) { rm(mat); next }
        cells <- colnames(mat)
        e_cells <- intersect(cells, rownames(md)[md$broad_celltype_final=="Epithelial"])
        f_cells <- intersect(cells, rownames(md)[md$broad_celltype_final=="Fibroblast_CAF"])
        if (length(e_cells) > 0) { e_det_sum <- e_det_sum + sum(mat[g, e_cells] > 0); e_total <- e_total + length(e_cells) }
        if (length(f_cells) > 0) { f_det_sum <- f_det_sum + sum(mat[g, f_cells] > 0); f_total <- f_total + length(f_cells) }
        rm(mat)
      }
    }
    e_det <- if (e_total > 0) e_det_sum / e_total else NA
    f_det <- if (f_total > 0) f_det_sum / f_total else NA
    ref_res <- rbind(ref_res, data.frame(gene=g, dataset=dname,
      epithelial_median_logCPM=round(e_med,4), fibroblast_median_logCPM=round(f_med,4),
      difference=round(e_med - f_med, 4),
      epithelial_detection_frac=round(e_det,4), fibroblast_detection_frac=round(f_det,4),
      stringsAsFactors=FALSE))
  }

  # LE gene localization (19G3)
  le_res <- data.frame()
  le_genes <- intersect(core_le_genes, rownames(lcpm))
  for (g in le_genes) {
    e_vals <- lcpm[g, grep("_Epithelial$", colnames(lcpm))]
    f_vals <- lcpm[g, grep("_Fibroblast_CAF$", colnames(lcpm))]
    e_med <- if (length(e_vals) > 0) median(e_vals) else 0
    f_med <- if (length(f_vals) > 0) median(f_vals) else 0
    e_det_sum <- 0; e_total <- 0; f_det_sum <- 0; f_total <- 0
    if (!inherits(obj@assays[["RNA"]], "Assay5")) {
      counts <- GetAssayData(obj, assay="RNA", layer="counts")
      if (g %in% rownames(counts)) {
        e_det_sum <- sum(counts[g, epi_idx] > 0); e_total <- length(epi_idx)
        f_det_sum <- sum(counts[g, fib_idx] > 0); f_total <- length(fib_idx)
      }
    } else {
      for (ln in get_count_layers(obj)) {
        mat <- LayerData(obj@assays[["RNA"]], layer=ln)
        if (!(g %in% rownames(mat))) { rm(mat); next }
        cells <- colnames(mat)
        e_cells <- intersect(cells, rownames(md)[md$broad_celltype_final=="Epithelial"])
        f_cells <- intersect(cells, rownames(md)[md$broad_celltype_final=="Fibroblast_CAF"])
        if (length(e_cells) > 0) { e_det_sum <- e_det_sum + sum(mat[g, e_cells] > 0); e_total <- e_total + length(e_cells) }
        if (length(f_cells) > 0) { f_det_sum <- f_det_sum + sum(mat[g, f_cells] > 0); f_total <- f_total + length(f_cells) }
        rm(mat)
      }
    }
    e_frac <- if (e_total > 0) e_det_sum / e_total else NA
    f_frac <- if (f_total > 0) f_det_sum / f_total else NA
    # Convert logCPM to approximate median count (reverse logCPM)
    le_res <- rbind(le_res, data.frame(gene=g, dataset=dname,
      epithelial_median_logCPM=round(e_med,4), fibroblast_median_logCPM=round(f_med,4),
      ratio=round(ifelse(f_med != 0, e_med - f_med, NA), 3),
      epithelial_positive_frac=round(e_frac,4), fibroblast_positive_frac=round(f_frac,4),
      stringsAsFactors=FALSE))
  }

  return(list(ref=ref_res, le=le_res))
}

# ============================================================================
# 19G4-5: CLUSTER IDENTITY AUDIT (PSEUDOBULK-BASED)
# ============================================================================
cat("\n=== 19G4-5: CLUSTER IDENTITY + AXIS CLASSIFICATION ===\n\n")

compute_cluster_scores_pseudobulk <- function(obj, md, dname, cluster_col) {
  cat(sprintf("--- %s ---\n", dname))
  md$cluster_group <- paste0(md$sample, "_C", md[[cluster_col]])

  all_scores_genes <- unique(c(epithelial_core_markers, squamous_markers, stromal_markers, hypoxia_markers, stress_markers))
  all_scores_genes <- intersect(all_scores_genes, rownames(obj@assays$RNA))

  pb <- pseudobulk_layerwise(obj, md, "cluster_group", target_genes=all_scores_genes)
  cat(sprintf("  Pseudobulk: %d genes x %d sample-cluster groups\n", nrow(pb$pb), ncol(pb$pb)))

  lcpm <- normalize_pseudobulk(pb$pb)

  # Identify epithelial clusters (those with broad_celltype_final=Epithelial in their metadata)
  epi_clusters <- sort(unique(md[[cluster_col]][md$broad_celltype_final == "Epithelial"]))
  cat(sprintf("  Epithelial clusters: %d\n", length(epi_clusters)))

  panels <- list(
    epithelial_core = intersect(epithelial_core_markers, rownames(lcpm)),
    squamous = intersect(squamous_markers, rownames(lcpm)),
    stromal_like = intersect(stromal_markers, rownames(lcpm)),
    hypoxia = intersect(hypoxia_markers, rownames(lcpm)),
    stress = intersect(stress_markers, rownames(lcpm))
  )

  res <- data.frame()
  for (cl in epi_clusters) {
    cl_cols <- grep(paste0("_C", cl, "$"), colnames(lcpm), value=TRUE)
    if (length(cl_cols) < 2) next
    cl_mat <- lcpm[, cl_cols, drop=FALSE]
    row_info <- data.frame(dataset=dname, cluster=cl, n_samples=length(cl_cols), stringsAsFactors=FALSE)
    for (pn in names(panels)) {
      gns <- panels[[pn]]
      if (length(gns) > 0) {
        # Mean across genes, then mean across samples
        gene_means <- rowMeans(cl_mat[gns, , drop=FALSE], na.rm=TRUE)
        row_info[[pn]] <- mean(gene_means, na.rm=TRUE)
      } else {
        row_info[[pn]] <- NA
      }
    }
    # Also count total cells in this cluster
    cl_cells <- sum(md[[cluster_col]] == cl & md$broad_celltype_final == "Epithelial")
    row_info$n_cells <- cl_cells
    res <- rbind(res, row_info)
  }
  # Rename score columns to include _score suffix for downstream consistency
  score_cols <- intersect(c("epithelial_core","squamous","stromal_like","hypoxia","stress"), colnames(res))
  colnames(res)[match(score_cols, colnames(res))] <- paste0(score_cols, "_score")
  return(res)
}

# ============================================================================
# 19G6-8: PATHWAY CLUSTER + CORRELATION + SAMPLE SCORES
# ============================================================================
cat("\n=== 19G6-8: PATHWAY CLUSTER, CORRELATION, SAMPLE SCORES ===\n\n")

compute_pw_cluster_pseudobulk <- function(obj, md, dname, cluster_col) {
  cat(sprintf("--- %s ---\n", dname))
  md$cluster_group <- paste0(md$sample, "_C", md[[cluster_col]])

  pw_genes <- character()
  for (pw in supported_pathways) {
    pw_genes <- c(pw_genes, intersect(hallmark_pathways[[pw]], rownames(obj@assays$RNA)))
  }
  pw_genes <- unique(pw_genes)
  cat(sprintf("  Pathway genes available: %d\n", length(pw_genes)))

  pb <- pseudobulk_layerwise(obj, md, "cluster_group", target_genes=pw_genes)
  lcpm <- normalize_pseudobulk(pb$pb)

  epi_clusters <- sort(unique(md[[cluster_col]][md$broad_celltype_final == "Epithelial"]))

  id_panels <- list(
    epithelial_core = intersect(epithelial_core_markers, rownames(lcpm)),
    stromal_like = intersect(stromal_markers, rownames(lcpm)),
    hypoxia = intersect(hypoxia_markers, rownames(lcpm)),
    stress = intersect(stress_markers, rownames(lcpm))
  )

  res <- data.frame()
  for (pw in supported_pathways) {
    pw_gns <- intersect(hallmark_pathways[[pw]], rownames(lcpm))
    if (length(pw_gns) < 5) next
    for (cl in epi_clusters) {
      cl_cols <- grep(paste0("_C", cl, "$"), colnames(lcpm), value=TRUE)
      if (length(cl_cols) < 2) next
      cl_mat <- lcpm[, cl_cols, drop=FALSE]
      row_info <- data.frame(dataset=dname, pathway=pw, cluster=cl, stringsAsFactors=FALSE)
      row_info$pathway_score <- mean(rowMeans(cl_mat[pw_gns, , drop=FALSE], na.rm=TRUE), na.rm=TRUE)
      for (pn in names(id_panels)) {
        gns <- id_panels[[pn]]
        row_info[[pn]] <- if (length(gns) > 0) mean(rowMeans(cl_mat[gns, , drop=FALSE], na.rm=TRUE), na.rm=TRUE) else NA
      }
      row_info$n_cells <- sum(md[[cluster_col]] == cl & md$broad_celltype_final == "Epithelial")
      res <- rbind(res, row_info)
    }
  }
  # Rename score columns to include _score suffix for downstream consistency
  score_cols <- intersect(c("epithelial_core","stromal_like","hypoxia","stress"), colnames(res))
  colnames(res)[match(score_cols, colnames(res))] <- paste0(score_cols, "_score")
  return(res)
}

compute_sample_scores_pseudobulk <- function(obj, md, dname) {
  cat(sprintf("--- %s sample scores ---\n", dname))
  md$sample_group <- md$sample

  all_genes <- unique(c(epithelial_core_markers, stromal_markers, hypoxia_markers, stress_markers))
  for (pw in supported_pathways) all_genes <- c(all_genes, hallmark_pathways[[pw]])
  all_genes <- intersect(unique(all_genes), rownames(obj@assays$RNA))

  pb <- pseudobulk_layerwise(obj, md, "sample_group", target_genes=all_genes)
  lcpm <- normalize_pseudobulk(pb$pb)

  # Restrict to epithelial samples (filter by epithelial cell count)
  epi_counts <- table(md$sample[md$broad_celltype_final == "Epithelial"])
  epi_samples <- names(epi_counts[epi_counts >= 3])

  panels <- list(
    epithelial_core = intersect(epithelial_core_markers, rownames(lcpm)),
    stromal_like = intersect(stromal_markers, rownames(lcpm)),
    hypoxia = intersect(hypoxia_markers, rownames(lcpm)),
    stress = intersect(stress_markers, rownames(lcpm))
  )
  for (pw in supported_pathways) panels[[pw]] <- intersect(hallmark_pathways[[pw]], rownames(lcpm))

  res <- data.frame()
  for (s in epi_samples) {
    if (!(s %in% colnames(lcpm))) next
    pat <- if ("patient" %in% colnames(md)) md$patient[md$sample==s][1] else s
    treat <- if ("treatment_group" %in% colnames(md)) md$treatment_group[md$sample==s][1] else NA
    row <- data.frame(dataset=dname, sample=s,
      patient=pat,
      treatment=treat,
      n_epi_cells=epi_counts[s], stringsAsFactors=FALSE)
    for (pn in names(panels)) {
      gns <- panels[[pn]]
      row[[pn]] <- if (length(gns) > 0 && s %in% colnames(lcpm)) mean(lcpm[gns, s], na.rm=TRUE) else NA
    }
    res <- rbind(res, row)
  }
  return(res)
}

# ============================================================================
# 19G11: EPITHELIAL CO-DETECTION (LAYER-WISE, SELECTED GENES ONLY)
# ============================================================================
cat("\n=== 19G11: EPITHELIAL CO-DETECTION ===\n\n")

run_codet_layerwise <- function(obj, md, dname) {
  cat(sprintf("--- %s ---\n", dname))
  epi_markers <- intersect(c("EPCAM","KRT8","KRT18","KRT19"), rownames(obj@assays$RNA))
  target_genes <- intersect(core_le_genes, rownames(obj@assays$RNA))
  needed <- unique(c(epi_markers, target_genes))

  # Accumulate counts across layers
  n_cells_total <- 0
  gene_pos <- setNames(rep(0, length(target_genes)), target_genes)
  both_pos <- setNames(rep(0, length(target_genes)), target_genes)
  epi_pos_total <- 0

  if (!inherits(obj@assays[["RNA"]], "Assay5")) {
    counts <- obj@assays[["RNA"]]@counts
    epi_cells <- colnames(counts)[md[colnames(counts), "broad_celltype_final"] == "Epithelial"]
    epi_mat <- counts[needed, epi_cells, drop=FALSE]
    epi_mk_pos <- colSums(epi_mat[intersect(epi_markers, rownames(epi_mat)), , drop=FALSE] > 0) > 0
    n_cells_total <- length(epi_cells)
    epi_pos_total <- sum(epi_mk_pos)
    for (g in target_genes) {
      if (g %in% rownames(epi_mat)) {
        tgt <- epi_mat[g, ] > 0
        gene_pos[g] <- sum(tgt)
        both_pos[g] <- sum(tgt & epi_mk_pos)
      }
    }
  } else {
    for (ln in get_count_layers(obj)) {
      mat <- LayerData(obj@assays[["RNA"]], layer=ln)
      cells <- colnames(mat)
      cell_types <- md[cells, "broad_celltype_final"]
      epi_cells <- cells[cell_types == "Epithelial"]
      if (length(epi_cells) == 0) { rm(mat); next }
      shared <- intersect(needed, rownames(mat))
      if (length(shared) == 0) { rm(mat); next }
      epi_mat <- mat[shared, epi_cells, drop=FALSE]
      epi_mk <- intersect(epi_markers, rownames(epi_mat))
      epi_mk_pos <- if (length(epi_mk) > 0) colSums(epi_mat[epi_mk, , drop=FALSE] > 0) > 0 else rep(FALSE, ncol(epi_mat))
      n_cells_total <- n_cells_total + length(epi_cells)
      epi_pos_total <- epi_pos_total + sum(epi_mk_pos)
      for (g in target_genes) {
        if (g %in% rownames(epi_mat)) {
          tgt <- epi_mat[g, ] > 0
          gene_pos[g] <- gene_pos[g] + sum(tgt)
          both_pos[g] <- both_pos[g] + sum(tgt & epi_mk_pos)
        }
      }
      rm(mat, epi_mat); gc()
    }
  }

  res <- data.frame()
  for (g in target_genes) {
    n <- gene_pos[g]
    b <- both_pos[g]
    res <- rbind(res, data.frame(dataset=dname, gene=g, n_epithelial_cells=n_cells_total,
      target_positive=n, target_positive_frac=round(ifelse(n_cells_total>0, n/n_cells_total, NA),4),
      target_and_epithelial_marker=b,
      codetection_fraction=round(ifelse(n>0, b/n, NA),4), stringsAsFactors=FALSE))
  }
  return(res)
}

# ============================================================================
# MAIN EXECUTION
# ============================================================================
cat("\n========================================\n")
cat("EXECUTING MEMORY-SAFE STEP 19G\n")
cat("========================================\n\n")

# --- Run 19G2-3 ---
audit19 <- run_lineage_audit(obj19, md19, "GSE197677")
audit22 <- run_lineage_audit(obj22, md22, "GSE221561")
write.csv(audit19$ref, file.path(g19g_dir, "STEP19G_GSE197677_LINEAGE_REFERENCE_EXPRESSION.csv"), row.names=FALSE)
write.csv(audit22$ref, file.path(g19g_dir, "STEP19G_GSE221561_LINEAGE_REFERENCE_EXPRESSION.csv"), row.names=FALSE)
le_loc <- rbind(audit19$le, audit22$le)
write.csv(le_loc, file.path(g19g_dir, "STEP19G_CORE_LEADING_EDGE_LINEAGE_LOCALIZATION.csv"), row.names=FALSE)
cat("Saved: lineage reference + LE localization CSVs\n\n")
saveRDS(list(audit19=audit19, audit22=audit22), file.path(ckpt_dir, "19G_CHECKPOINT_lineage_reference.rds"))

# --- Run 19G4-5 ---
cl19 <- compute_cluster_scores_pseudobulk(obj19, md19, "GSE197677", "seurat_clusters")
cl22 <- compute_cluster_scores_pseudobulk(obj22, md22, "GSE221561", "author_seurat_clusters")
write.csv(cl19, file.path(g19g_dir, "STEP19G_GSE197677_EPITHELIAL_CLUSTER_IDENTITY.csv"), row.names=FALSE)
write.csv(cl22, file.path(g19g_dir, "STEP19G_GSE221561_EPITHELIAL_CLUSTER_IDENTITY.csv"), row.names=FALSE)

# Axis classification
for (df_name in c("cl19","cl22")) {
  df <- get(df_name)
  df$identity_delta <- df$epithelial_core_score - df$stromal_like_score
  df$axis_class <- ifelse(df$identity_delta > 0.5, "EPITHELIAL_DOMINANT",
                   ifelse(df$identity_delta < -0.5, "STROMAL_LIKE_PROFILE", "INTERMEDIATE"))
  assign(df_name, df)
}
flags <- rbind(cl19[,c("dataset","cluster","n_cells","epithelial_core_score","stromal_like_score","identity_delta","axis_class")],
               cl22[,c("dataset","cluster","n_cells","epithelial_core_score","stromal_like_score","identity_delta","axis_class")])
write.csv(flags, file.path(g19g_dir, "STEP19G_EPITHELIAL_CLUSTER_IDENTITY_FLAGS.csv"), row.names=FALSE)
cat("Saved: cluster identity + flags CSVs\n")
cat("GSE197677 axis:\n"); print(table(cl19$axis_class))
cat("GSE221561 axis:\n"); print(table(cl22$axis_class))
saveRDS(list(cl19=cl19, cl22=cl22), file.path(ckpt_dir, "19G_CHECKPOINT_cluster_identity.rds"))
cat("\n")

# --- Run 19G6-8 ---
pw19 <- compute_pw_cluster_pseudobulk(obj19, md19, "GSE197677", "seurat_clusters")
pw22 <- compute_pw_cluster_pseudobulk(obj22, md22, "GSE221561", "author_seurat_clusters")
pw_all <- rbind(pw19, pw22)
write.csv(pw_all, file.path(g19g_dir, "STEP19G_SUPPORTED_PATHWAY_CLUSTER_LOCALIZATION.csv"), row.names=FALSE)

# Correlation
corr_all <- data.frame()
for (dname in c("GSE197677","GSE221561")) {
  ds <- pw_all[pw_all$dataset==dname, ]
  for (pw in supported_pathways) {
    dp <- ds[ds$pathway==pw, ]
    if (nrow(dp) < 2) next
    corr_all <- rbind(corr_all, data.frame(dataset=dname, pathway=pw,
      rho_vs_stromal=round(cor(dp$pathway_score, dp$stromal_like_score, method="spearman", use="complete.obs"),4),
      rho_vs_epithelial=round(cor(dp$pathway_score, dp$epithelial_core_score, method="spearman", use="complete.obs"),4),
      rho_vs_stress=round(cor(dp$pathway_score, dp$stress_score, method="spearman", use="complete.obs"),4),
      n_clusters=nrow(dp), stringsAsFactors=FALSE))
  }
}
write.csv(corr_all, file.path(g19g_dir, "STEP19G_PATHWAY_IDENTITY_CORRELATION.csv"), row.names=FALSE)

# Sample scores
samp19 <- compute_sample_scores_pseudobulk(obj19, md19, "GSE197677")
samp22 <- compute_sample_scores_pseudobulk(obj22, md22, "GSE221561")
samp_all <- rbind(samp19, samp22)
write.csv(samp_all, file.path(g19g_dir, "STEP19G_SAMPLE_EPITHELIAL_STATE_SCORES.csv"), row.names=FALSE)

# Treatment group summary
cat("\nTreatment group summary (DESCRIPTIVE_SECONDARY_ONLY):\n")
for (dname in c("GSE197677","GSE221561")) {
  ds <- samp_all[samp_all$dataset==dname, ]
  score_cols <- setdiff(colnames(ds), c("dataset","sample","patient","treatment","n_epi_cells"))
  groups <- unique(ds$treatment)
  if (length(groups) == 2) {
    g1 <- ds[ds$treatment==groups[1], ]; g2 <- ds[ds$treatment==groups[2], ]
    cat(sprintf("  %s:\n", dname))
    for (sc in score_cols) {
      v1 <- g1[[sc]]; v1 <- v1[!is.na(v1)]
      v2 <- g2[[sc]]; v2 <- v2[!is.na(v2)]
      if (length(v1)>0 && length(v2)>0) {
        cat(sprintf("    %s: %s median=%.4f, %s median=%.4f, diff=%.4f\n",
            sc, groups[1], median(v1), groups[2], median(v2), median(v1)-median(v2)))
      }
    }
  }
}
saveRDS(list(pw_all=pw_all, corr_all=corr_all, samp_all=samp_all), file.path(ckpt_dir, "19G_CHECKPOINT_pathway_localization.rds"))
cat("\nSaved: pathway + sample score CSVs\n\n")

# --- 19G10: Concentration ---
cat("=== 19G10: PATHWAY SIGNAL CONCENTRATION ===\n\n")
conc <- data.frame()
for (dname in c("GSE197677","GSE221561")) {
  ds <- pw_all[pw_all$dataset==dname, ]
  for (pw in supported_pathways) {
    dp <- ds[ds$pathway==pw, ]
    if (nrow(dp)==0) next
    dp$weight <- ifelse(dp$pathway_score>0, dp$pathway_score * dp$n_cells, 0)
    tot_w <- sum(dp$weight)
    if (tot_w > 0) {
      dp$contrib <- dp$weight / tot_w
      dp <- dp[order(-dp$contrib), ]
      top1 <- dp$contrib[1]; top2 <- sum(dp$contrib[1:min(2,nrow(dp))]); top3 <- sum(dp$contrib[1:min(3,nrow(dp))])
      top_frac <- dp$n_cells[1] / sum(dp$n_cells)
    } else { top1 <- top2 <- top3 <- top_frac <- NA }
    conc <- rbind(conc, data.frame(dataset=dname, pathway=pw, n_clusters=nrow(dp),
      top1_contribution=round(top1,4), top2_contribution=round(top2,4), top3_contribution=round(top3,4),
      top_cluster_cell_fraction=round(top_frac,4), stringsAsFactors=FALSE))
  }
}
write.csv(conc, file.path(g19g_dir, "STEP19G_PATHWAY_CLUSTER_CONCENTRATION.csv"), row.names=FALSE)
cat("Saved: STEP19G_PATHWAY_CLUSTER_CONCENTRATION.csv\n\n")

# --- 19G11: Co-detection ---
codet19 <- run_codet_layerwise(obj19, md19, "GSE197677")
codet22 <- run_codet_layerwise(obj22, md22, "GSE221561")
codet <- rbind(codet19, codet22)
write.csv(codet, file.path(g19g_dir, "STEP19G_CORE_GENE_EPITHELIAL_CODETECTION.csv"), row.names=FALSE)
cat("Saved: STEP19G_CORE_GENE_EPITHELIAL_CODETECTION.csv\n\n")
saveRDS(list(codet19=codet19, codet22=codet22), file.path(ckpt_dir, "19G_CHECKPOINT_core_gene_codetection.rds"))

# ============================================================================
# 19G12-16: PROVENANCE, FIGURES, CLASSIFICATION, FINAL OUTPUT
# ============================================================================
cat("\n=== 19G12: CLUSTER PROVENANCE ===\n\n")
clust_map <- read.csv(file.path(base_dir, "06_tables/GSE197677/celltype/GSE197677_FINAL_cluster_to_broad_celltype_mapping.csv"), stringsAsFactors=FALSE)
prov19 <- merge(clust_map[, c("cluster","broad_celltype_final","n_cells","review_status")],
                cl19[, c("cluster","epithelial_core_score","stromal_like_score","hypoxia_score","stress_score")],
                by="cluster", all.x=TRUE)
for (pw in supported_pathways) {
  ps <- pw19[pw19$pathway==pw, c("cluster","pathway_score")]
  colnames(ps)[2] <- paste0("score_", gsub("HALLMARK_","",pw))
  prov19 <- merge(prov19, ps, by="cluster", all.x=TRUE)
}
write.csv(prov19, file.path(g19g_dir, "STEP19G_GSE197677_CLUSTER_PROVENANCE.csv"), row.names=FALSE)
prov22 <- cl22[, c("cluster","n_cells","epithelial_core_score","stromal_like_score","hypoxia_score","stress_score")]
for (pw in supported_pathways) {
  ps <- pw22[pw22$pathway==pw, c("cluster","pathway_score")]
  colnames(ps)[2] <- paste0("score_", gsub("HALLMARK_","",pw))
  prov22 <- merge(prov22, ps, by="cluster", all.x=TRUE)
}
write.csv(prov22, file.path(g19g_dir, "STEP19G_GSE221561_CLUSTER_PROVENANCE.csv"), row.names=FALSE)
cat("Saved: cluster provenance CSVs\n\n")

# --- 19G13: Figures ---
cat("=== 19G13: FIGURES ===\n\n")

# Fig 1: Marker heatmap
ref19 <- audit19$ref; ref22 <- audit22$ref
ref_comb <- rbind(ref19, ref22)
ref_wide <- reshape2::dcast(ref_comb, gene ~ dataset, value.var="difference")
rownames(ref_wide) <- ref_wide$gene; ref_wide$gene <- NULL
ref_mat <- as.matrix(ref_wide)
gene_ord <- intersect(unique(c(epithelial_core_markers, squamous_markers, stromal_markers, hypoxia_markers, stress_markers)), rownames(ref_mat))
ref_mat <- ref_mat[gene_ord, ]
png(file.path(fig_dir, "STEP19G_EPITHELIAL_VS_STROMAL_MARKER_HEATMAP.png"), width=800, height=1200, res=150)
par(mar=c(4,10,2,2))
image(t(ref_mat[nrow(ref_mat):1,]), col=rev(colorRampPalette(c("blue","white","red"))(100)), axes=FALSE)
axis(2, at=seq(0,1,length.out=nrow(ref_mat)), labels=rev(rownames(ref_mat)), las=2, cex.axis=0.7)
axis(1, at=c(0,1), labels=colnames(ref_mat), cex.axis=0.8)
title(main="Epithelial - Fibroblast logCPM Difference", cex.main=0.9)
dev.off()

# Fig 2: Cluster identity heatmaps
for (dname in c("GSE197677","GSE221561")) {
  cl_df <- if(dname=="GSE197677") cl19 else cl22
  mat <- as.matrix(cl_df[, c("epithelial_core_score","squamous_score","stromal_like_score","hypoxia_score","stress_score")])
  rownames(mat) <- paste0("C", cl_df$cluster)
  png(file.path(fig_dir, sprintf("STEP19G_%s_EPITHELIAL_CLUSTER_IDENTITY_HEATMAP.png", dname)),
      width=800, height=max(400, nrow(mat)*30), res=150)
  par(mar=c(4,8,2,2))
  image(t(mat[nrow(mat):1,]), col=colorRampPalette(c("blue","white","red"))(100), axes=FALSE)
  axis(2, at=seq(0,1,length.out=nrow(mat)), labels=rev(rownames(mat)), las=2, cex.axis=0.7)
  axis(1, at=seq(0,1,length.out=ncol(mat)), labels=colnames(mat), las=2, cex.axis=0.6)
  title(main=sprintf("%s Cluster Identity", dname), cex.main=0.9)
  dev.off()
}

# Fig 3: Seven pathways by cluster
for (dname in c("GSE197677","GSE221561")) {
  pw_ds <- pw_all[pw_all$dataset==dname, ]
  pw_wide <- reshape2::dcast(pw_ds, pathway ~ cluster, value.var="pathway_score")
  rownames(pw_wide) <- gsub("HALLMARK_","",pw_wide$pathway); pw_wide$pathway <- NULL
  pw_mat <- as.matrix(pw_wide)
  png(file.path(fig_dir, sprintf("STEP19G_%s_SEVEN_PATHWAY_BY_EPITHELIAL_CLUSTER.png", dname)),
      width=max(800, ncol(pw_mat)*40), height=600, res=150)
  par(mar=c(8,12,2,2))
  image(t(pw_mat[nrow(pw_mat):1,]), col=colorRampPalette(c("blue","white","red"))(100), axes=FALSE)
  axis(2, at=seq(0,1,length.out=nrow(pw_mat)), labels=rev(rownames(pw_mat)), las=2, cex.axis=0.7)
  axis(1, at=seq(0,1,length.out=ncol(pw_mat)), labels=paste0("C",colnames(pw_mat)), las=2, cex.axis=0.6)
  title(main=sprintf("%s Hallmark by Cluster", dname), cex.main=0.9)
  dev.off()
}

# Fig 4: LE gene lineage heatmap
le_wide <- reshape2::dcast(le_loc, gene ~ dataset, value.var="ratio")
rownames(le_wide) <- le_wide$gene; le_wide$gene <- NULL
le_mat <- as.matrix(le_wide)
png(file.path(fig_dir, "STEP19G_CORE_LEADING_EDGE_LINEAGE_HEATMAP.png"), width=600, height=600, res=150)
par(mar=c(4,8,2,2))
image(t(le_mat[nrow(le_mat):1,]), col=colorRampPalette(c("blue","white","red"))(100), axes=FALSE)
axis(2, at=seq(0,1,length.out=nrow(le_mat)), labels=rev(rownames(le_mat)), las=2, cex.axis=0.7)
axis(1, at=c(0,1), labels=colnames(le_mat), cex.axis=0.8)
title(main="Epi logCPM - Fib logCPM", cex.main=0.9)
dev.off()
cat("Saved: all figures\n\n")

# ============================================================================
# 19G14: PRIMARY LOCALIZATION CLASSIFICATION
# ============================================================================
cat("=== 19G14: PRIMARY LOCALIZATION CLASSIFICATION ===\n\n")

le_jaccard_df <- read.csv(file.path(pw_dir, "STEP19E_HALLMARK_LEADING_EDGE_OVERLAP.csv"), stringsAsFactors=FALSE)

corr19 <- corr_all[corr_all$dataset=="GSE197677", ]
corr22 <- corr_all[corr_all$dataset=="GSE221561", ]
conc19 <- conc[conc$dataset=="GSE197677", ]
conc22 <- conc[conc$dataset=="GSE221561", ]

class_df <- data.frame()
for (pw in supported_pathways) {
  rho19_val <- corr19$rho_vs_stromal[corr19$pathway==pw]
  rho22_val <- corr22$rho_vs_stromal[corr22$pathway==pw]
  c19_row <- conc19[conc19$pathway==pw, ]
  c22_row <- conc22[conc22$pathway==pw, ]
  le_j <- le_jaccard_df$jaccard[le_jaccard_df$pathway==pw]
  if (length(le_j)==0) le_j <- NA

  cd19 <- codet19$codetection_fraction[codet19$gene %in% core_le_genes]
  cd22 <- codet22$codetection_fraction[codet22$gene %in% core_le_genes]
  avg_cd <- mean(c(cd19, cd22), na.rm=TRUE)

  # Safe extraction with length checks
  rho19_safe <- if (length(rho19_val) > 0) rho19_val else NA
  rho22_safe <- if (length(rho22_val) > 0) rho22_val else NA
  c19_top1 <- if (nrow(c19_row) > 0) c19_row$top1_contribution else NA
  c22_top1 <- if (nrow(c22_row) > 0) c22_row$top1_contribution else NA
  le_j_safe <- if (length(le_j) > 0) le_j else NA

  strong_stromal <- (!is.na(rho19_val) && rho19_val > 0.5) || (!is.na(rho22_val) && rho22_val > 0.5)
  high_conc <- (!is.na(c19_top1) && c19_top1 > 0.4) ||
               (!is.na(c22_top1) && c22_top1 > 0.4)
  good_cd <- avg_cd > 0.3

  if (strong_stromal && high_conc && !good_cd) {
    loc <- "POSSIBLE_LINEAGE_MISCLASSIFICATION_SIGNAL"
  } else if (strong_stromal || high_conc) {
    loc <- "MIXED_EPITHELIAL_STROMAL_PROGRAM"
  } else {
    loc <- "EPITHELIAL_LOCALIZATION_SUPPORTED"
  }

  class_df <- rbind(class_df, data.frame(
    pathway=pw,
    GSE197677_pathway_stromal_rho=round(rho19_safe,4),
    GSE221561_pathway_stromal_rho=round(rho22_safe,4),
    GSE197677_cluster_concentration=round(c19_top1,4),
    GSE221561_cluster_concentration=round(c22_top1,4),
    cross_dataset_LE_Jaccard=round(le_j_safe,4),
    core_gene_epithelial_codetection=round(avg_cd,4),
    localization_class=loc, stringsAsFactors=FALSE))
}
write.csv(class_df, file.path(g19g_dir, "STEP19G_PATHWAY_LOCALIZATION_FINAL.csv"), row.names=FALSE)
cat("Saved: STEP19G_PATHWAY_LOCALIZATION_FINAL.csv\n\n")
cat("Classification:\n")
for (i in 1:nrow(class_df)) cat(sprintf("  %s: %s\n", class_df$pathway[i], class_df$localization_class[i]))

# ============================================================================
# 19G15: INTERPRETATION
# ============================================================================
cat("\n=== 19G15: INTERPRETATION ===\n\n")
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
cat(sprintf("\nOverall: %s\n", overall))

interp <- "# Step 19G: Epithelial Pathway Localization Interpretation\n\n"
interp <- paste0(interp, "## Scope\nThis audit determines the likely cellular localization of 7 robust Hallmark pathways from Step 19F.\n\n")
interp <- paste0(interp, "**Does NOT:** establish EMT causality, malignant-cell specificity, alter Step 19F, or reclassify/remove cells.\n\n")
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
# 19G16: FINAL OUTPUT
# ============================================================================
cat("\n============================================\n")
cat("STEP 19G COMPLETE\n")
cat("EPITHELIAL PATHWAY LOCALIZATION AUDIT\n")
cat("============================================\n\n")
cat("JoinLayers used: NO\n")
cat("Full dense matrix created: NO\n")
cat(sprintf("GSE221561 count layers processed: %d\n", length(layers22)))
cat("GSE197677: 71116 cells loaded, v3 assay\n")
cat("GSE221561: 55150 cells loaded, v5 assay (NOT joined)\n\n")

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
  cat(sprintf("  Fibroblast reference: %d\n", sum(md$broad_celltype_final=="Fibroblast_CAF")))
  cat(sprintf("  EPITHELIAL_DOMINANT: %d\n", sum(cl_df$axis_class=="EPITHELIAL_DOMINANT", na.rm=TRUE)))
  cat(sprintf("  INTERMEDIATE: %d\n", sum(cl_df$axis_class=="INTERMEDIATE", na.rm=TRUE)))
  cat(sprintf("  STROMAL_LIKE_PROFILE: %d\n", sum(cl_df$axis_class=="STROMAL_LIKE_PROFILE", na.rm=TRUE)))
}

cat("\n--------------------------------------------\n\n")
cat("Core recurrent leading-edge genes:\n")
for (g in core_le_genes) {
  e19 <- codet19$codetection_fraction[codet19$gene==g]
  e22 <- codet22$codetection_fraction[codet22$gene==g]
  cd <- mean(c(e19, e22), na.rm=TRUE)
  # Reference localization from le_loc
  le19 <- le_loc[le_loc$gene==g & le_loc$dataset=="GSE197677", ]
  le22 <- le_loc[le_loc$gene==g & le_loc$dataset=="GSE221561", ]
  e19_diff <- if(nrow(le19)>0) le19$epithelial_median_logCPM - le19$fibroblast_median_logCPM else NA
  e22_diff <- if(nrow(le22)>0) le22$epithelial_median_logCPM - le22$fibroblast_median_logCPM else NA
  avg_diff <- mean(c(e19_diff, e22_diff), na.rm=TRUE)
  loc <- if (!is.na(cd) && cd > 0.5) "EPITHelial_COPOSITIVE" else if (!is.na(cd) && cd > 0.2) "MIXED" else "LOW_COLOCALIZATION"
  cat(sprintf("  %s: Epi-Fib diff=%.2f, CoDet=%.1f%%, %s\n", g, avg_diff, cd*100, loc))
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
cat("Cells reclassified: 0\n")
cat("Cells removed: 0\n")
cat("Samples removed: 0\n")
cat("Step19F results changed: NO\n\n")
cat("--------------------------------------------\n\n")
cat(sprintf("Overall Step 19G status:\n%s\n\n", overall))
cat("STOP HERE.\nDO NOT START STEP 19H AUTOMATICALLY.\n")
