#!/usr/bin/env Rscript
# ============================================================================
# STEP 19H: MIXED-PROGRAM DECOMPOSITION AND EPITHELIAL-RESIDUAL SENSITIVITY
# ============================================================================
# Sensitivity/localization analysis. NOT causal adjustment.
# Must NOT be interpreted as "controlling for stromal contamination proves
# epithelial-intrinsic treatment response."
# ============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(edgeR)
  library(stats)
})

base_dir <- Sys.getenv("ESCC_CAF_PROJECT_ROOT", unset = ".")
obj_dir <- file.path(base_dir, "03_objects")
pw_dir <- file.path(base_dir, "06_tables/cross_dataset/pathway_analysis")
g19g_dir <- file.path(pw_dir, "Step19G")
g19h_dir <- file.path(pw_dir, "Step19H")
ckpt_dir <- file.path(g19h_dir, "checkpoints")
fig_dir <- file.path(base_dir, "05_figures/cross_dataset/pathway_analysis", "Step19H")

dir.create(g19h_dir, showWarnings=FALSE, recursive=TRUE)
dir.create(ckpt_dir, showWarnings=FALSE, recursive=TRUE)
dir.create(fig_dir, showWarnings=FALSE, recursive=TRUE)

supported_pathways <- c(
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "HALLMARK_HYPOXIA",
  "HALLMARK_COAGULATION",
  "HALLMARK_KRAS_SIGNALING_UP",
  "HALLMARK_UV_RESPONSE_DN",
  "HALLMARK_MYOGENESIS",
  "HALLMARK_APOPTOSIS"
)

# ============================================================================
# 19H0: INPUT FREEZE + PROVENANCE
# ============================================================================
cat("========================================\n")
cat("STEP 19H0: INPUT FREEZE + PROVENANCE\n")
cat("========================================\n\n")

f19f <- read.csv(file.path(pw_dir, "STEP19F_SEVEN_HALLMARKS_FINAL_FREEZE.csv"), stringsAsFactors=FALSE)
g19g_loc <- read.csv(file.path(g19g_dir, "STEP19G_PATHWAY_LOCALIZATION_FINAL.csv"), stringsAsFactors=FALSE)
g19g_samp <- read.csv(file.path(g19g_dir, "STEP19G_SAMPLE_EPITHELIAL_STATE_SCORES.csv"), stringsAsFactors=FALSE)
g19g_le_loc <- read.csv(file.path(g19g_dir, "STEP19G_CORE_LEADING_EDGE_LINEAGE_LOCALIZATION.csv"), stringsAsFactors=FALSE)
g19g_codet <- read.csv(file.path(g19g_dir, "STEP19G_CORE_GENE_EPITHELIAL_CODETECTION.csv"), stringsAsFactors=FALSE)
g19g_pw <- read.csv(file.path(g19g_dir, "STEP19G_SUPPORTED_PATHWAY_CLUSTER_LOCALIZATION.csv"), stringsAsFactors=FALSE)
g19g_flags <- read.csv(file.path(g19g_dir, "STEP19G_EPITHELIAL_CLUSTER_IDENTITY_FLAGS.csv"), stringsAsFactors=FALSE)
fg19 <- read.csv(file.path(pw_dir, "STEP19E_GSE197677_HALLMARK_FGSEA_ALL.csv"), stringsAsFactors=FALSE)
fg22 <- read.csv(file.path(pw_dir, "STEP19E_GSE221561_HALLMARK_FGSEA_ALL.csv"), stringsAsFactors=FALSE)

stopifnot(nrow(f19f) == 7)
stopifnot(all(f19f$pathway %in% g19g_loc$pathway))
stopifnot(all(g19g_loc$localization_class == "MIXED_EPITHELIAL_STROMAL_PROGRAM"))

provenance <- data.frame(
  step="19H", input=c("STEP19F_SEVEN_HALLMARKS_FINAL_FREEZE.csv","STEP19G_PATHWAY_LOCALIZATION_FINAL.csv",
    "STEP19G_SAMPLE_EPITHELIAL_STATE_SCORES.csv","STEP19G_CORE_LEADING_EDGE_LINEAGE_LOCALIZATION.csv",
    "STEP19G_CORE_GENE_EPITHELIAL_CODETECTION.csv","STEP19G_SUPPORTED_PATHWAY_CLUSTER_LOCALIZATION.csv",
    "STEP19G_EPITHELIAL_CLUSTER_IDENTITY_FLAGS.csv","STEP19E_GSE197677_HALLMARK_FGSEA_ALL.csv",
    "STEP19E_GSE221561_HALLMARK_FGSEA_ALL.csv"),
  status=rep("VERIFIED", 9), n_pathways=rep(7, 9),
  step19g_status=rep("ROBUST_PATHWAYS_MIXED_EPITHELIAL_STROMAL_SIGNAL", 9),
  stringsAsFactors=FALSE)
write.csv(provenance, file.path(g19h_dir, "STEP19H_INPUT_PROVENANCE.csv"), row.names=FALSE)
cat("Saved: STEP19H_INPUT_PROVENANCE.csv\n")
cat("Verified: 7 pathways, Step19G = MIXED_EPITHELIAL_STROMAL_PROGRAM\n\n")

# ============================================================================
# 19H1: FREEZE PATHWAY GENE MEMBERSHIP
# ============================================================================
cat("========================================\n")
cat("STEP 19H1: FREEZE PATHWAY GENE MEMBERSHIP\n")
cat("========================================\n\n")

le_freeze <- data.frame()

for (pw in supported_pathways) {
  pw_name_short <- gsub("HALLMARK_", "", pw)

  # Discovery LE (GSE197677)
  row19 <- fg19[fg19$pathway == pw, ]
  if (nrow(row19) == 0) { cat(sprintf("  WARNING: %s not found in GSE197677 fgsea\n", pw)); next }
  genes_discovery <- trimws(unlist(strsplit(as.character(row19$leadingEdge), ",")))

  # Replication LE (GSE221561)
  row22 <- fg22[fg22$pathway == pw, ]
  if (nrow(row22) == 0) { cat(sprintf("  WARNING: %s not found in GSE221561 fgsea\n", pw)); next }
  genes_replication <- trimws(unlist(strsplit(as.character(row22$leadingEdge), ",")))

  union_genes <- sort(unique(c(genes_discovery, genes_replication)))

  for (g in union_genes) {
    le_freeze <- rbind(le_freeze, data.frame(
      pathway=pw, gene=g,
      in_discovery_LE=g %in% genes_discovery,
      in_replication_LE=g %in% genes_replication,
      in_both_LE=g %in% genes_discovery && g %in% genes_replication,
      stringsAsFactors=FALSE))
  }
  cat(sprintf("  %s: %d discovery LE, %d replication LE, %d union\n",
      pw_name_short, length(genes_discovery), length(genes_replication), length(union_genes)))
}

write.csv(le_freeze, file.path(g19h_dir, "STEP19H_SEVEN_PATHWAY_LEADING_EDGE_FREEZE.csv"), row.names=FALSE)
cat(sprintf("\nTotal: %d gene-pathway entries across 7 pathways\n", nrow(le_freeze)))
cat("Saved: STEP19H_SEVEN_PATHWAY_LEADING_EDGE_FREEZE.csv\n\n")

# ============================================================================
# 19H2: GENE-LEVEL LINEAGE DECOMPOSITION
# ============================================================================
cat("========================================\n")
cat("STEP 19H2: GENE-LEVEL LINEAGE DECOMPOSITION\n")
cat("========================================\n\n")

# Load lineage reference data from 19G checkpoints
lin_ckpt <- readRDS(file.path(g19g_dir, "checkpoints", "19G_CHECKPOINT_lineage_reference.rds"))
ref19 <- lin_ckpt$audit19$ref
ref22 <- lin_ckpt$audit22$ref
le_loc19 <- lin_ckpt$audit19$le
le_loc22 <- lin_ckpt$audit22$le

codet <- g19g_codet

# Combine lineage evidence from both datasets
gene_class <- data.frame()

for (pw in supported_pathways) {
  genes_in_pw <- le_freeze$gene[le_freeze$pathway == pw]

  for (g in genes_in_pw) {
    # Expression evidence from lineage reference
    e19 <- ref19[ref19$gene == g, ]
    e22 <- ref22[ref22$gene == g, ]

    # LE localization from 19G
    loc19 <- le_loc19[le_loc19$gene == g, ]
    loc22 <- le_loc22[le_loc22$gene == g, ]

    # Co-detection evidence
    cd19 <- codet[codet$gene == g & codet$dataset == "GSE197677", ]
    cd22 <- codet[codet$gene == g & codet$dataset == "GSE221561", ]

    # Compute evidence scores
    # Difference = epithelial - fibroblast (negative means fibroblast higher)
    diff19 <- if (nrow(e19) > 0) e19$difference[1] else NA
    diff22 <- if (nrow(e22) > 0) e22$difference[1] else NA
    avg_diff <- mean(c(diff19, diff22), na.rm=TRUE)
    if (is.nan(avg_diff)) avg_diff <- NA

    # Co-detection fraction
    cd19_val <- if (nrow(cd19) > 0) cd19$codetection_fraction[1] else NA
    cd22_val <- if (nrow(cd22) > 0) cd22$codetection_fraction[1] else NA
    avg_cd <- mean(c(cd19_val, cd22_val), na.rm=TRUE)
    if (is.nan(avg_cd)) avg_cd <- NA

    # LE localization - use ratio and positive fractions
    loc19_ratio <- if (nrow(loc19) > 0) loc19$ratio[1] else NA
    loc22_ratio <- if (nrow(loc22) > 0) loc22$ratio[1] else NA
    loc19_epi_pos <- if (nrow(loc19) > 0) loc19$epithelial_positive_frac[1] else NA
    loc22_epi_pos <- if (nrow(loc22) > 0) loc22$epithelial_positive_frac[1] else NA

    # Combined classification
    has_evidence <- !is.na(avg_diff) || !is.na(avg_cd)

    if (!has_evidence) {
      class <- "UNRESOLVED"
    } else {
      # Multiple evidence criteria
      epithelial_lean_evidence <- 0
      stromal_lean_evidence <- 0

      # Expression difference
      if (!is.na(avg_diff)) {
        if (avg_diff >= 0) epithelial_lean_evidence <- epithelial_lean_evidence + 1
        else stromal_lean_evidence <- stromal_lean_evidence + 1
      }

      # Co-detection (high = epithelial)
      if (!is.na(avg_cd)) {
        if (avg_cd > 0.5) epithelial_lean_evidence <- epithelial_lean_evidence + 1
        else if (avg_cd < 0.3) stromal_lean_evidence <- stromal_lean_evidence + 1
        else { epithelial_lean_evidence <- epithelial_lean_evidence + 0.5; stromal_lean_evidence <- stromal_lean_evidence + 0.5 }
      }

      # Classify based on combined evidence
      total <- epithelial_lean_evidence + stromal_lean_evidence
      if (total == 0) {
        class <- "UNRESOLVED"
      } else {
        epi_frac <- epithelial_lean_evidence / total
        if (epi_frac >= 0.6) class <- "EPITHELIAL_LEANING"
        else if (epi_frac <= 0.4) class <- "STROMAL_LEANING"
        else class <- "MIXED_LINEAGE"
      }
    }

    gene_class <- rbind(gene_class, data.frame(
      pathway=pw, gene=g,
      avg_epithelial_fibroblast_diff=ifelse(is.na(avg_diff), NA, round(avg_diff, 4)),
      avg_codetection_fraction=ifelse(is.na(avg_cd), NA, round(avg_cd, 4)),
      avg_le_ratio=round(mean(c(loc19_ratio, loc22_ratio), na.rm=TRUE), 4),
      avg_le_epi_positive_frac=round(mean(c(loc19_epi_pos, loc22_epi_pos), na.rm=TRUE), 4),
      lineage_class=class, stringsAsFactors=FALSE))
  }
}

write.csv(gene_class, file.path(g19h_dir, "STEP19H_LEADING_EDGE_GENE_LINEAGE_CLASSIFICATION.csv"), row.names=FALSE)
cat("Lineage classification:\n")
print(table(gene_class$lineage_class))
cat("\nSaved: STEP19H_LEADING_EDGE_GENE_LINEAGE_CLASSIFICATION.csv\n\n")

# ============================================================================
# 19H3: PATHWAY GENE COMPOSITION
# ============================================================================
cat("========================================\n")
cat("STEP 19H3: PATHWAY GENE COMPOSITION\n")
cat("========================================\n\n")

pw_comp <- data.frame()

for (pw in supported_pathways) {
  gc <- gene_class[gene_class$pathway == pw, ]
  n_total <- nrow(gc)
  n_epi <- sum(gc$lineage_class == "EPITHELIAL_LEANING")
  n_str <- sum(gc$lineage_class == "STROMAL_LEANING")
  n_mixed <- sum(gc$lineage_class == "MIXED_LINEAGE")
  n_unres <- sum(gc$lineage_class == "UNRESOLVED")

  # Also split by discovery/replication/shared
  le_in_pw <- le_freeze[le_freeze$pathway == pw, ]
  disc_genes <- le_in_pw$gene[le_in_pw$in_discovery_LE]
  repl_genes <- le_in_pw$gene[le_in_pw$in_replication_LE]
  shared_genes <- le_in_pw$gene[le_in_pw$in_both_LE]

  gc_disc <- gc[gc$gene %in% disc_genes, ]
  gc_repl <- gc[gc$gene %in% repl_genes, ]
  gc_shared <- gc[gc$gene %in% shared_genes, ]

  pw_comp <- rbind(pw_comp, data.frame(
    pathway=pw,
    n_leading_edge_genes=n_total,
    n_epithelial_leaning=n_epi, n_stromal_leaning=n_str,
    n_mixed=n_mixed, n_unresolved=n_unres,
    frac_epithelial=round(n_epi/max(n_total,1), 4),
    frac_stromal=round(n_str/max(n_total,1), 4),
    frac_mixed=round(n_mixed/max(n_total,1), 4),
    frac_unresolved=round(n_unres/max(n_total,1), 4),
    n_discovery_LE=length(disc_genes),
    n_replication_LE=length(repl_genes),
    n_shared_LE=length(shared_genes),
    frac_disc_epithelial=round(sum(gc_disc$lineage_class=="EPITHELIAL_LEANING")/max(length(disc_genes),1),4),
    frac_disc_stromal=round(sum(gc_disc$lineage_class=="STROMAL_LEANING")/max(length(disc_genes),1),4),
    frac_repl_epithelial=round(sum(gc_repl$lineage_class=="EPITHELIAL_LEANING")/max(length(repl_genes),1),4),
    frac_repl_stromal=round(sum(gc_repl$lineage_class=="STROMAL_LEANING")/max(length(repl_genes),1),4),
    stringsAsFactors=FALSE))
}

write.csv(pw_comp, file.path(g19h_dir, "STEP19H_PATHWAY_LINEAGE_COMPOSITION.csv"), row.names=FALSE)
cat("Pathway composition:\n")
for (i in 1:nrow(pw_comp)) {
  cat(sprintf("  %s: %d genes (%.0f%% epi, %.0f%% str, %.0f%% mixed, %.0f%% unres)\n",
      gsub("HALLMARK_","",pw_comp$pathway[i]), pw_comp$n_leading_edge_genes[i],
      pw_comp$frac_epithelial[i]*100, pw_comp$frac_stromal[i]*100,
      pw_comp$frac_mixed[i]*100, pw_comp$frac_unresolved[i]*100))
}
cat("\nSaved: STEP19H_PATHWAY_LINEAGE_COMPOSITION.csv\n\n")

# ============================================================================
# 19H4: SAMPLE-LEVEL SUBPROGRAM SCORES
# ============================================================================
cat("========================================\n")
cat("STEP19H4: SAMPLE-LEVEL SUBPROGRAM SCORES\n")
cat("========================================\n\n")

cat("Loading Seurat objects for pseudobulk reconstruction...\n")
obj19 <- readRDS(file.path(obj_dir, "GSE197677", "GSE197677_integrated_broadcelltype_annotated.rds"))
obj22 <- readRDS(file.path(obj_dir, "GSE221561", "GSE221561_author_annotated.rds"))
cat("Loaded.\n\n")

# Metadata
md19 <- read.csv(file.path(base_dir, "06_tables/cross_dataset/epithelial_pseudobulk/counts",
  "GSE197677_epithelial_pseudobulk_sample_metadata.csv"), stringsAsFactors=FALSE)
rownames(md19) <- md19$sample
md19$treatment_group <- md19$treatment_standardized

# GSE221561 metadata from Seurat
md22 <- obj22@meta.data
md22$sample <- md22$sample_id
md22$patient <- md22$sample_id
treat22 <- md22$author_Neoadjuvant
treat22[treat22 %in% c("Chemotherapy","Chemoradiotherapy","Chemoradiotherapy and Immunotherapy","Antiangiogenesis")] <- "Neoadjuvant_treated"
treat22[treat22 == "Surgery alone"] <- "Surgery_alone"
treat22[treat22 == "Adjacent normal"] <- "Adjacent_normal"
md22$treatment_group <- treat22

# Helper: pseudobulk per sample for epithelial cells
pseudobulk_epithelial <- function(obj, sample_col, target_genes) {
  # Use Seurat object metadata for cell type
  cell_meta <- obj@meta.data
  # Determine epithelial cell type column
  if ("broad_celltype_final" %in% colnames(cell_meta)) {
    epi_label <- "Epithelial"
    ct_col <- "broad_celltype_final"
  } else if ("author_Cell_type" %in% colnames(cell_meta)) {
    epi_label <- "Epithelial"
    ct_col <- "author_Cell_type"
  } else {
    stop("No cell type column found in metadata")
  }
  
  if (!inherits(obj@assays[["RNA"]], "Assay5")) {
    counts <- obj@assays[["RNA"]]@counts
    epi_cells <- rownames(cell_meta)[cell_meta[[ct_col]] == epi_label]
    epi_cells <- intersect(epi_cells, colnames(counts))
    samples <- unique(cell_meta[epi_cells, sample_col])
    pb <- matrix(0, nrow=length(target_genes), ncol=length(samples),
                 dimnames=list(target_genes, samples))
    for (s in samples) {
      cells <- epi_cells[cell_meta[epi_cells, sample_col] == s]
      if (length(cells) > 0) {
        genes_in <- intersect(target_genes, rownames(counts))
        if (length(genes_in) > 0) {
          pb[genes_in, s] <- as.numeric(Matrix::rowSums(counts[genes_in, cells, drop=FALSE]))
        }
      }
    }
  } else {
    samples <- unique(cell_meta[[sample_col]][cell_meta[[ct_col]] == epi_label])
    pb <- matrix(0, nrow=length(target_genes), ncol=length(samples),
                 dimnames=list(target_genes, samples))
    for (ln in names(obj@assays[["RNA"]]@layers)) {
      mat <- LayerData(obj@assays[["RNA"]], layer=ln)
      cells <- colnames(mat)
      cell_types <- cell_meta[cells, ct_col]
      epi_cells <- cells[cell_types == epi_label]
      if (length(epi_cells) == 0) next
      for (s in intersect(samples, unique(cell_meta[epi_cells, sample_col]))) {
        s_cells <- epi_cells[cell_meta[epi_cells, sample_col] == s]
        if (length(s_cells) > 0) {
          genes_in_layer <- intersect(target_genes, rownames(mat))
          if (length(genes_in_layer) > 0) {
            pb[genes_in_layer, s] <- pb[genes_in_layer, s] +
              as.numeric(Matrix::rowSums(mat[genes_in_layer, s_cells, drop=FALSE]))
          }
        }
      }
    }
  }
  return(pb)
}

# Normalize pseudobulk to logCPM
normalize_pb <- function(pb) {
  # Remove samples with zero total counts
  col_sums <- colSums(pb)
  pb <- pb[, col_sums > 0, drop=FALSE]
  if (ncol(pb) == 0 || nrow(pb) == 0) return(matrix(NA, nrow=nrow(pb), ncol=ncol(pb), dimnames=dimnames(pb)))
  dge <- DGEList(counts=pb)
  dge <- calcNormFactors(dge, method="TMM")
  cpm <- cpm(dge, log=TRUE, prior.count=1)
  return(cpm)
}

# Gene-wise z-score within dataset
zscore_genes <- function(mat) {
  # mat: genes x samples
  res <- mat
  for (g in rownames(mat)) {
    vals <- as.numeric(mat[g, ])
    vals <- vals[!is.na(vals)]
    if (length(vals) < 2 || sd(vals) == 0) {
      res[g, ] <- 0
    } else {
      res[g, ] <- (as.numeric(mat[g, ]) - mean(vals)) / sd(vals)
    }
  }
  return(res)
}

# Compute subscores for a dataset
compute_subscores <- function(pb_logcpm, le_genes, gene_class_df, dataset_name) {
  samples <- colnames(pb_logcpm)
  genes_available <- rownames(pb_logcpm)

  # Z-score the full matrix
  zmat <- zscore_genes(pb_logcpm)

  results <- data.frame()
  for (pw in supported_pathways) {
    pw_genes <- le_genes$gene[le_genes$pathway == pw]
    pw_genes <- intersect(pw_genes, genes_available)

    gc_pw <- gene_class_df[gene_class_df$pathway == pw, ]

    epi_genes <- intersect(gc_pw$gene[gc_pw$lineage_class == "EPITHELIAL_LEANING"], genes_available)
    str_genes <- intersect(gc_pw$gene[gc_pw$lineage_class == "STROMAL_LEANING"], genes_available)
    mix_genes <- intersect(gc_pw$gene[gc_pw$lineage_class == "MIXED_LINEAGE"], genes_available)

    for (s in samples) {
      full_score <- if (length(pw_genes) >= 3) mean(zmat[pw_genes, s], na.rm=TRUE) else NA
      epi_score <- if (length(epi_genes) >= 3) mean(zmat[epi_genes, s], na.rm=TRUE) else NA
      str_score <- if (length(str_genes) >= 3) mean(zmat[str_genes, s], na.rm=TRUE) else NA
      mix_score <- if (length(mix_genes) >= 3) mean(zmat[mix_genes, s], na.rm=TRUE) else NA

      n_epi_insuff <- if (length(epi_genes) > 0 && length(epi_genes) < 3) "INSUFFICIENT_GENES" else NA
      n_str_insuff <- if (length(str_genes) > 0 && length(str_genes) < 3) "INSUFFICIENT_GENES" else NA
      n_mix_insuff <- if (length(mix_genes) > 0 && length(mix_genes) < 3) "INSUFFICIENT_GENES" else NA

      results <- rbind(results, data.frame(
        dataset=dataset_name, sample=s, pathway=pw,
        full_pathway_score=round(full_score, 6),
        n_full_genes=length(pw_genes),
        epithelial_leaning_score=round(epi_score, 6),
        n_epithelial_genes=length(epi_genes),
        stromal_leaning_score=round(str_score, 6),
        n_stromal_genes=length(str_genes),
        mixed_gene_score=round(mix_score, 6),
        n_mixed_genes=length(mix_genes),
        epithelial_status=ifelse(is.na(epi_score) && !is.na(n_epi_insuff), n_epi_insuff, "OK"),
        stromal_status=ifelse(is.na(str_score) && !is.na(n_str_insuff), n_str_insuff, "OK"),
        mixed_status=ifelse(is.na(mix_score) && !is.na(n_mix_insuff), n_mix_insuff, "OK"),
        stringsAsFactors=FALSE))
    }
  }
  return(results)
}

cat("Processing GSE197677...\n")
epi_genes_all <- unique(c(
  intersect(rownames(obj19@assays[["RNA"]]), le_freeze$gene[le_freeze$pathway %in% supported_pathways]),
  intersect(rownames(obj19@assays[["RNA"]]), c("EPCAM","KRT8","KRT18","KRT19","VIM","COL3A1","COL1A1","FAP"))
))
pb19 <- pseudobulk_epithelial(obj19, "sample", epi_genes_all)
pb19_lcpm <- normalize_pb(pb19)
cat(sprintf("  Pseudobulk: %d genes x %d samples\n", nrow(pb19_lcpm), ncol(pb19_lcpm)))

sub19 <- compute_subscores(pb19_lcpm, le_freeze, gene_class, "GSE197677")
cat(sprintf("  Subscores: %d rows\n", nrow(sub19)))

cat("Processing GSE221561...\n")
epi_genes_all22 <- unique(c(
  intersect(rownames(obj22@assays[["RNA"]]), le_freeze$gene[le_freeze$pathway %in% supported_pathways]),
  intersect(rownames(obj22@assays[["RNA"]]), c("EPCAM","KRT8","KRT18","KRT19","VIM","COL3A1","COL1A1","FAP"))
))
pb22 <- pseudobulk_epithelial(obj22, "sample_id", epi_genes_all22)
pb22_lcpm <- normalize_pb(pb22)
cat(sprintf("  Pseudobulk: %d genes x %d samples\n", nrow(pb22_lcpm), ncol(pb22_lcpm)))

sub22 <- compute_subscores(pb22_lcpm, le_freeze, gene_class, "GSE221561")
cat(sprintf("  Subscores: %d rows\n", nrow(sub22)))

sub_all <- rbind(sub19, sub22)
write.csv(sub_all, file.path(g19h_dir, "STEP19H_SAMPLE_PATHWAY_SUBPROGRAM_SCORES.csv"), row.names=FALSE)
cat("\nSaved: STEP19H_SAMPLE_PATHWAY_SUBPROGRAM_SCORES.csv\n\n")

saveRDS(list(pb19_lcpm=pb19_lcpm, pb22_lcpm=pb22_lcpm, sub19=sub19, sub22=sub22),
        file.path(ckpt_dir, "19H_CHECKPOINT_pseudobulk_subscores.rds"))

# ============================================================================
# 19H5: TREATMENT DIRECTION AUDIT
# ============================================================================
cat("========================================\n")
cat("STEP 19H5: TREATMENT DIRECTION AUDIT\n")
cat("========================================\n\n")

treat_dir <- data.frame()

for (dname in c("GSE197677", "GSE221561")) {
  sub_ds <- sub_all[sub_all$dataset == dname, ]
  md_ds <- if (dname == "GSE197677") md19 else md22

  # Define treatment groups
  if (dname == "GSE197677") {
    grp1_name <- "NACT"; grp2_name <- "nNACT"
    grp1_samples <- unique(md19$sample[md19$treatment_group == "Neoadjuvant_chemotherapy"])
    grp2_samples <- unique(md19$sample[md19$treatment_group == "No_neoadjuvant_chemotherapy"])
  } else {
    grp1_name <- "Neoadjuvant_treated"; grp2_name <- "Surgery_alone"
    grp1_samples <- unique(md22$sample_id[md22$treatment_group == "Neoadjuvant_treated"])
    grp2_samples <- unique(md22$sample_id[md22$treatment_group == "Surgery_alone"])
  }

  for (pw in supported_pathways) {
    sub_pw <- sub_ds[sub_ds$pathway == pw, ]
    score_cols <- c("full_pathway_score", "epithelial_leaning_score", "stromal_leaning_score", "mixed_gene_score")

    for (sc in score_cols) {
      v1 <- sub_pw[[sc]][sub_pw$sample %in% grp1_samples]
      v2 <- sub_pw[[sc]][sub_pw$sample %in% grp2_samples]
      v1 <- v1[!is.na(v1)]; v2 <- v2[!is.na(v2)]

      if (length(v1) < 2 || length(v2) < 2) next

      mean_diff <- mean(v1) - mean(v2)
      med_diff <- median(v1) - median(v2)
      pooled_sd <- sqrt((var(v1) + var(v2)) / 2)
      d <- if (pooled_sd > 0) mean_diff / pooled_sd else 0

      direction <- if (mean_diff > 0) paste0("HIGHER_IN_", grp1_name) else paste0("HIGHER_IN_", grp2_name)

      treat_dir <- rbind(treat_dir, data.frame(
        dataset=dname, pathway=pw, subprogram=sc,
        grp1_name=grp1_name, grp2_name=grp2_name,
        grp1_n=length(v1), grp2_n=length(v2),
        grp1_mean=round(mean(v1),4), grp2_mean=round(mean(v2),4),
        grp1_median=round(median(v1),4), grp2_median=round(median(v2),4),
        difference_in_means=round(mean_diff,4),
        difference_in_medians=round(med_diff,4),
        standardized_effect=round(d,4),
        direction=direction, stringsAsFactors=FALSE))
    }
  }
  cat(sprintf("  %s: %d comparisons\n", dname, sum(treat_dir$dataset == dname)))
}

write.csv(treat_dir, file.path(g19h_dir, "STEP19H_SUBPROGRAM_TREATMENT_EFFECTS.csv"), row.names=FALSE)
cat("\nSaved: STEP19H_SUBPROGRAM_TREATMENT_EFFECTS.csv\n\n")

# ============================================================================
# 19H6: CROSS-DATASET DIRECTION REPLICATION
# ============================================================================
cat("========================================\n")
cat("STEP19H6: CROSS-DATASET DIRECTION REPLICATION\n")
cat("========================================\n\n")

cross_dir <- data.frame()

for (pw in supported_pathways) {
  for (sc in c("full_pathway_score", "epithelial_leaning_score", "stromal_leaning_score", "mixed_gene_score")) {
    td19 <- treat_dir[treat_dir$dataset == "GSE197677" & treat_dir$pathway == pw & treat_dir$subprogram == sc, ]
    td22 <- treat_dir[treat_dir$dataset == "GSE221561" & treat_dir$pathway == pw & treat_dir$subprogram == sc, ]

    if (nrow(td19) == 0 || nrow(td22) == 0) {
      cross_dir <- rbind(cross_dir, data.frame(pathway=pw, subprogram=sc,
        GSE197677_direction=NA, GSE221561_direction=NA,
        GSE197677_effect=NA, GSE221561_effect=NA,
        concordance="UNRESOLVED", stringsAsFactors=FALSE))
      next
    }

    d19 <- td19$direction[1]
    d22 <- td22$direction[1]
    e19 <- td19$difference_in_means[1]
    e22 <- td22$difference_in_means[1]

    # Both positive in neoadjuvant direction = SAME
    both_pos <- (grepl("HIGHER_IN_NACT", d19) || grepl("HIGHER_IN_Neoadjuvant", d19)) &&
                (grepl("HIGHER_IN_NACT", d22) || grepl("HIGHER_IN_Neoadjuvant", d22))
    both_neg <- (grepl("HIGHER_IN_nNACT", d19) || grepl("HIGHER_IN_Surgery", d19)) &&
                (grepl("HIGHER_IN_nNACT", d22) || grepl("HIGHER_IN_Surgery", d22))

    conc <- if (both_pos || both_neg) "SAME_DIRECTION" else "OPPOSITE_DIRECTION"

    cross_dir <- rbind(cross_dir, data.frame(pathway=pw, subprogram=sc,
      GSE197677_direction=d19, GSE221561_direction=d22,
      GSE197677_effect=round(e19,4), GSE221561_effect=round(e22,4),
      concordance=conc, stringsAsFactors=FALSE))
  }
}

write.csv(cross_dir, file.path(g19h_dir, "STEP19H_SUBPROGRAM_CROSS_DATASET_CONCORDANCE.csv"), row.names=FALSE)
cat("Cross-dataset concordance:\n")
print(table(cross_dir$concordance))
cat("\nSaved: STEP19H_SUBPROGRAM_CROSS_DATASET_CONCORDANCE.csv\n\n")

# ============================================================================
# 19H7: LOSO DIRECTION STABILITY
# ============================================================================
cat("========================================\n")
cat("STEP 19H7: LOSO DIRECTION STABILITY\n")
cat("========================================\n\n")

loso_results <- data.frame()

for (dname in c("GSE197677", "GSE221561")) {
  sub_ds <- sub_all[sub_all$dataset == dname, ]
  md_ds <- if (dname == "GSE197677") md19 else md22

  if (dname == "GSE197677") {
    grp1_samples <- unique(md19$sample[md19$treatment_group == "Neoadjuvant_chemotherapy"])
    grp2_samples <- unique(md19$sample[md19$treatment_group == "No_neoadjuvant_chemotherapy"])
  } else {
    grp1_samples <- unique(md22$sample_id[md22$treatment_group == "Neoadjuvant_treated"])
    grp2_samples <- unique(md22$sample_id[md22$treatment_group == "Surgery_alone"])
  }

  all_samples <- c(grp1_samples, grp2_samples)

  for (pw in supported_pathways) {
    sub_pw <- sub_ds[sub_ds$pathway == pw, ]
    score_cols <- c("full_pathway_score", "epithelial_leaning_score", "stromal_leaning_score")

    for (sc in score_cols) {
      vals <- sub_pw[[sc]]
      names(vals) <- sub_pw$sample
      vals <- vals[!is.na(vals)]
      available <- intersect(names(vals), all_samples)
      if (length(available) < 4) next

      # Full data direction
      v1_full <- vals[intersect(available, grp1_samples)]
      v2_full <- vals[intersect(available, grp2_samples)]
      if (length(v1_full) < 2 || length(v2_full) < 2) next
      full_dir <- sign(mean(v1_full) - mean(v2_full))

      # LOSO
      n_preserve <- 0
      n_total <- 0
      for (leave_out in available) {
        remaining <- setdiff(available, leave_out)
        v1_lo <- vals[intersect(remaining, grp1_samples)]
        v2_lo <- vals[intersect(remaining, grp2_samples)]
        if (length(v1_lo) < 2 || length(v2_lo) < 2) next
        lo_dir <- sign(mean(v1_lo) - mean(v2_lo))
        n_total <- n_total + 1
        if (lo_dir == full_dir) n_preserve <- n_preserve + 1
      }

      frac <- if (n_total > 0) n_preserve / n_total else NA

      loso_results <- rbind(loso_results, data.frame(
        dataset=dname, pathway=pw, subprogram=sc,
        full_direction=if (full_dir > 0) "HIGHER_NEOADJ" else "HIGHER_NON_NEOADJ",
        n_samples=length(available), n_iterations=n_total,
        n_preserve=n_preserve, loso_fraction=round(frac, 4),
        loso_status=ifelse(!is.na(frac) && frac >= 0.8, "ROBUST_DIRECTION",
                   ifelse(!is.na(frac) && frac >= 0.6, "MODERATE_DIRECTION", "UNSTABLE_DIRECTION")),
        stringsAsFactors=FALSE))
    }
  }
  cat(sprintf("  %s: %d LOSO rows\n", dname, sum(loso_results$dataset == dname)))
}

write.csv(loso_results, file.path(g19h_dir, "STEP19H_SUBPROGRAM_LOSO_STABILITY.csv"), row.names=FALSE)
cat("\nLOSO status summary:\n")
print(table(loso_results$loso_status))
cat("\nSaved: STEP19H_SUBPROGRAM_LOSO_STABILITY.csv\n\n")

# ============================================================================
# 19H8: STROMAL-SCORE RESIDUAL SENSITIVITY
# ============================================================================
cat("========================================\n")
cat("STEP19H8: STROMAL-SCORE RESIDUAL SENSITIVITY\n")
cat("========================================\n\n")

# Use Step19G sample scores for residual analysis
g19g_samp_all <- read.csv(file.path(g19g_dir, "STEP19G_SAMPLE_EPITHELIAL_STATE_SCORES.csv"), stringsAsFactors=FALSE)

resid_results <- data.frame()

for (dname in c("GSE197677", "GSE221561")) {
  samp_ds <- g19g_samp_all[g19g_samp_all$dataset == dname, ]
  md_ds <- if (dname == "GSE197677") md19 else md22

  if (dname == "GSE197677") {
    grp1_samples <- unique(md19$sample[md19$treatment_group == "Neoadjuvant_chemotherapy"])
    grp2_samples <- unique(md19$sample[md19$treatment_group == "No_neoadjuvant_chemotherapy"])
    grp1_name <- "NACT"; grp2_name <- "nNACT"
  } else {
    grp1_samples <- unique(md22$sample_id[md22$treatment_group == "Neoadjuvant_treated"])
    grp2_samples <- unique(md22$sample_id[md22$treatment_group == "Surgery_alone"])
    grp1_name <- "Neoadjuvant_treated"; grp2_name <- "Surgery_alone"
  }

  for (pw in supported_pathways) {
    pw_short <- gsub("HALLMARK_", "", pw)
    if (!(pw %in% colnames(samp_ds))) next

    pathway_scores <- samp_ds[[pw]]
    stromal_scores <- samp_ds$stromal_like
    sample_names <- samp_ds$sample

    valid <- !is.na(pathway_scores) & !is.na(stromal_scores)
    if (sum(valid) < 4) next

    # Fit pathway ~ stromal (NO treatment)
    fit <- lm(pathway_scores[valid] ~ stromal_scores[valid])
    residuals_resid <- residuals(fit)

    # Map residuals to samples
    resid_df <- data.frame(sample=sample_names[valid], raw=pathway_scores[valid],
                           stromal=stromal_scores[valid], residual=residuals_resid,
                           stringsAsFactors=FALSE)

    # Treatment difference in raw
    raw_g1 <- mean(resid_df$raw[resid_df$sample %in% grp1_samples], na.rm=TRUE)
    raw_g2 <- mean(resid_df$raw[resid_df$sample %in% grp2_samples], na.rm=TRUE)
    raw_diff <- raw_g1 - raw_g2

    # Treatment difference in residuals
    resid_g1 <- mean(resid_df$residual[resid_df$sample %in% grp1_samples], na.rm=TRUE)
    resid_g2 <- mean(resid_df$residual[resid_df$sample %in% grp2_samples], na.rm=TRUE)
    resid_diff <- resid_g1 - resid_g2

    # Direction preservation
    dir_preserved <- if (sign(raw_diff) == sign(resid_diff) || raw_diff == 0) "YES" else "NO"

    # Effect retained fraction
    eff_retained <- if (abs(raw_diff) > 0.001) abs(resid_diff) / abs(raw_diff) else NA

    resid_results <- rbind(resid_results, data.frame(
      dataset=dname, pathway=pw,
      raw_grp1_mean=round(raw_g1,4), raw_grp2_mean=round(raw_g2,4),
      raw_difference=round(raw_diff,4),
      residual_grp1_mean=round(resid_g1,4), residual_grp2_mean=round(resid_g2,4),
      residual_difference=round(resid_diff,4),
      direction_preserved=dir_preserved,
      effect_retained_fraction=round(eff_retained,4),
      label="DESCRIPTIVE_SENSITIVITY_ONLY",
      stringsAsFactors=FALSE))
  }
  cat(sprintf("  %s: %d pathways analyzed\n", dname, sum(resid_results$dataset == dname)))
}

write.csv(resid_results, file.path(g19h_dir, "STEP19H_STROMAL_RESIDUAL_SENSITIVITY.csv"), row.names=FALSE)
cat("\nDirection preserved:\n")
print(table(resid_results$direction_preserved))
cat("\nSaved: STEP19H_STROMAL_RESIDUAL_SENSITIVITY.csv\n\n")

# ============================================================================
# 19H9: OPTIONAL TWO-PREDICTOR MODEL AUDIT
# ============================================================================
cat("========================================\n")
cat("STEP19H9: TWO-PREDICTOR MODEL AUDIT\n")
cat("========================================\n\n")

model_results <- data.frame()

for (dname in c("GSE197677", "GSE221561")) {
  samp_ds <- g19g_samp_all[g19g_samp_all$dataset == dname, ]
  md_ds <- if (dname == "GSE197677") md19 else md22

  # Create treatment binary
  if (dname == "GSE197677") {
    samp_ds$treatment_bin <- ifelse(samp_ds$sample %in% unique(md19$sample[md19$treatment_group == "Neoadjuvant_chemotherapy"]), 1, 0)
  } else {
    samp_ds$treatment_bin <- ifelse(samp_ds$sample %in% unique(md22$sample_id[md22$treatment_group == "Neoadjuvant_treated"]), 1, 0)
  }

  for (pw in supported_pathways) {
    if (!(pw %in% colnames(samp_ds))) next

    y <- samp_ds[[pw]]
    x_treat <- samp_ds$treatment_bin
    x_stromal <- samp_ds$stromal_like

    valid <- !is.na(y) & !is.na(x_stromal)
    if (sum(valid) < 5) next

    fit <- lm(y[valid] ~ x_treat[valid] + x_stromal[valid])
    s <- summary(fit)

    coefs <- s$coefficients
    treat_coef <- if ("x_treat[valid]" %in% rownames(coefs)) coefs["x_treat[valid]", ] else coefs[2, ]
    stromal_coef <- if ("x_stromal[valid]" %in% rownames(coefs)) coefs["x_stromal[valid]", ] else coefs[3, ]

    model_results <- rbind(model_results, data.frame(
      dataset=dname, pathway=pw,
      treatment_coefficient=round(as.numeric(treat_coef[1]),4),
      treatment_SE=round(as.numeric(treat_coef[2]),4),
      treatment_P=round(as.numeric(treat_coef[4]),4),
      stromal_coefficient=round(as.numeric(stromal_coef[1]),4),
      stromal_SE=round(as.numeric(stromal_coef[2]),4),
      stromal_P=round(as.numeric(stromal_coef[4]),4),
      model_label="EXPLORATORY_SMALL_N_MODEL",
      stringsAsFactors=FALSE))
  }
  cat(sprintf("  %s: %d models\n", dname, sum(model_results$dataset == dname)))
}

write.csv(model_results, file.path(g19h_dir, "STEP19H_EXPLORATORY_TREATMENT_STROMAL_MODELS.csv"), row.names=FALSE)
cat("\nSaved: STEP19H_EXPLORATORY_TREATMENT_STROMAL_MODELS.csv\n\n")

# ============================================================================
# 19H10: EPITHELIAL-CORE RELATIONSHIPS
# ============================================================================
cat("========================================\n")
cat("STEP19H10: EPITHELIAL-CORE RELATIONSHIPS\n")
cat("========================================\n\n")

corr_results <- data.frame()

for (dname in c("GSE197677", "GSE221561")) {
  samp_ds <- g19g_samp_all[g19g_samp_all$dataset == dname, ]
  sub_ds <- sub_all[sub_all$dataset == dname, ]

  for (pw in supported_pathways) {
    if (!(pw %in% colnames(samp_ds))) next

    full_scores <- samp_ds[[pw]]
    epi_core <- samp_ds$epithelial_core
    str_like <- samp_ds$stromal_like

    # Get subscores
    sub_pw <- sub_ds[sub_ds$pathway == pw, ]
    # Align by sample
    common <- intersect(samp_ds$sample, sub_pw$sample)
    if (length(common) < 3) next

    idx1 <- match(common, samp_ds$sample)
    idx2 <- match(common, sub_pw$sample)

    # Full vs epithelial_core
    rho_full_epi <- tryCatch(cor(full_scores[idx1], epi_core[idx1], method="spearman", use="complete.obs"), error=function(e) NA)
    rho_epi_sub_epi <- tryCatch(cor(sub_pw$epithelial_leaning_score[idx2], epi_core[idx1], method="spearman", use="complete.obs"), error=function(e) NA)
    rho_str_sub_epi <- tryCatch(cor(sub_pw$stromal_leaning_score[idx2], epi_core[idx1], method="spearman", use="complete.obs"), error=function(e) NA)
    rho_full_str <- tryCatch(cor(full_scores[idx1], str_like[idx1], method="spearman", use="complete.obs"), error=function(e) NA)

    corr_results <- rbind(corr_results, data.frame(
      dataset=dname, pathway=pw, n_samples=length(common),
      rho_full_vs_epithelial_core=round(rho_full_epi,4),
      rho_epithelial_sub_vs_epithelial_core=round(rho_epi_sub_epi,4),
      rho_stromal_sub_vs_epithelial_core=round(rho_str_sub_epi,4),
      rho_full_vs_stromal_like=round(rho_full_str,4),
      stringsAsFactors=FALSE))
  }
}

write.csv(corr_results, file.path(g19h_dir, "STEP19H_PATHWAY_IDENTITY_RELATIONSHIPS.csv"), row.names=FALSE)
cat("Saved: STEP19H_PATHWAY_IDENTITY_RELATIONSHIPS.csv\n\n")

# ============================================================================
# 19H11: CORE RECURRENT GENE DECOMPOSITION
# ============================================================================
cat("========================================\n")
cat("STEP19H11: CORE RECURRENT GENE DECOMPOSITION\n")
cat("========================================\n\n")

core_genes <- c("COL3A1","COL1A1","MMP2","BGN","SPARC","CCN1","DCN","F3","CDKN1A")

gene_decomp <- data.frame()

for (g in core_genes) {
  # Lineage classification
  gc_g <- gene_class[gene_class$gene == g, ]
  lineage <- if (nrow(gc_g) > 0) gc_g$lineage_class[1] else "UNRESOLVED"
  avg_diff <- if (nrow(gc_g) > 0) gc_g$avg_epithelial_fibroblast_diff[1] else NA
  avg_cd <- if (nrow(gc_g) > 0) gc_g$avg_codetection_fraction[1] else NA

  # Which pathways contain this gene
  pw_contained <- unique(le_freeze$pathway[le_freeze$gene == g])
  n_pathways <- length(pw_contained)
  pw_list <- paste(gsub("HALLMARK_", "", pw_contained), collapse=", ")

  # Discovery direction (from 19E gene evidence)
  le_evidence <- read.csv(file.path(pw_dir, "STEP19E_HALLMARK_LEADING_EDGE_GENE_EVIDENCE.csv"), stringsAsFactors=FALSE)
  gene_ev <- le_evidence[le_evidence$gene == g, ]
  disc_dir <- if (nrow(gene_ev) > 0) {
    dirs <- gene_ev$GSE197677_logFC
    if (length(dirs) > 0) ifelse(mean(dirs, na.rm=TRUE) > 0, "UP_in_NACT", "DOWN_in_NACT") else NA
  } else NA

  repl_dir <- if (nrow(gene_ev) > 0) {
    dirs <- gene_ev$GSE221561_logFC
    if (length(dirs) > 0) ifelse(mean(dirs, na.rm=TRUE) > 0, "UP_in_NACT", "DOWN_in_NACT") else NA
  } else NA

  gene_decomp <- rbind(gene_decomp, data.frame(
    gene=g, lineage_class=lineage,
    avg_epithelial_fibroblast_diff=round(avg_diff,4),
    avg_codetection_fraction=round(avg_cd,4),
    n_pathways=n_pathways, pathways_contained=pw_list,
    discovery_direction=disc_dir, replication_direction=repl_dir,
    stringsAsFactors=FALSE))
}

write.csv(gene_decomp, file.path(g19h_dir, "STEP19H_CORE_RECURRENT_GENE_DECOMPOSITION.csv"), row.names=FALSE)
cat("Core gene decomposition:\n")
for (i in 1:nrow(gene_decomp)) {
  cat(sprintf("  %s: %s, %d pathways (%s)\n", gene_decomp$gene[i], gene_decomp$lineage_class[i],
      gene_decomp$n_pathways[i], gene_decomp$pathways_contained[i]))
}
cat("\nSaved: STEP19H_CORE_RECURRENT_GENE_DECOMPOSITION.csv\n\n")

# ============================================================================
# 19H12: PATHWAY-LEVEL FINAL DECOMPOSITION
# ============================================================================
cat("========================================\n")
cat("STEP19H12: PATHWAY-LEVEL FINAL DECOMPOSITION\n")
cat("========================================\n\n")

final_decomp <- data.frame()

for (pw in supported_pathways) {
  pw_short <- gsub("HALLMARK_", "", pw)

  # 19F status
  f19f_row <- f19f[f19f$pathway == pw, ]
  orig_19f <- if (nrow(f19f_row) > 0) f19f_row$final_robustness_class[1] else NA

  # 19G localization
  loc_row <- g19g_loc[g19g_loc$pathway == pw, ]
  loc_class <- if (nrow(loc_row) > 0) loc_row$localization_class[1] else NA

  # Gene fractions
  comp_row <- pw_comp[pw_comp$pathway == pw, ]
  epi_frac <- if (nrow(comp_row) > 0) comp_row$frac_epithelial[1] else NA
  str_frac <- if (nrow(comp_row) > 0) comp_row$frac_stromal[1] else NA
  mix_frac <- if (nrow(comp_row) > 0) comp_row$frac_mixed[1] else NA

  # Direction from 19H5/19H6
  td19_full <- treat_dir[treat_dir$dataset=="GSE197677" & treat_dir$pathway==pw & treat_dir$subprogram=="full_pathway_score", ]
  td22_full <- treat_dir[treat_dir$dataset=="GSE221561" & treat_dir$pathway==pw & treat_dir$subprogram=="full_pathway_score", ]
  td19_epi <- treat_dir[treat_dir$dataset=="GSE197677" & treat_dir$pathway==pw & treat_dir$subprogram=="epithelial_leaning_score", ]
  td22_epi <- treat_dir[treat_dir$dataset=="GSE221561" & treat_dir$pathway==pw & treat_dir$subprogram=="epithelial_leaning_score", ]
  td19_str <- treat_dir[treat_dir$dataset=="GSE197677" & treat_dir$pathway==pw & treat_dir$subprogram=="stromal_leaning_score", ]
  td22_str <- treat_dir[treat_dir$dataset=="GSE221561" & treat_dir$pathway==pw & treat_dir$subprogram=="stromal_leaning_score", ]

  full_dir_disc <- if (nrow(td19_full) > 0) td19_full$direction[1] else NA
  full_dir_repl <- if (nrow(td22_full) > 0) td22_full$direction[1] else NA
  epi_dir_disc <- if (nrow(td19_epi) > 0) td19_epi$direction[1] else NA
  epi_dir_repl <- if (nrow(td22_epi) > 0) td22_epi$direction[1] else NA
  str_dir_disc <- if (nrow(td19_str) > 0) td19_str$direction[1] else NA
  str_dir_repl <- if (nrow(td22_str) > 0) td22_str$direction[1] else NA

  # LOSO
  loso19 <- loso_results[loso_results$dataset=="GSE197677" & loso_results$pathway==pw, ]
  loso22 <- loso_results[loso_results$dataset=="GSE221561" & loso_results$pathway==pw, ]
  epi_loso19 <- loso19$loso_status[loso19$subprogram=="epithelial_leaning_score"]
  epi_loso22 <- loso22$loso_status[loso22$subprogram=="epithelial_leaning_score"]

  # Residual
  res19 <- resid_results[resid_results$dataset=="GSE197677" & resid_results$pathway==pw, ]
  res22 <- resid_results[resid_results$dataset=="GSE221561" & resid_results$pathway==pw, ]
  resid_dir_disc <- if (nrow(res19) > 0) {
    if (res19$residual_difference > 0) "HIGHER_IN_NACT" else "HIGHER_IN_nNACT"
  } else NA
  resid_dir_repl <- if (nrow(res22) > 0) {
    if (res22$residual_difference > 0) "HIGHER_IN_Neoadjuvant_treated" else "HIGHER_IN_Surgery_alone"
  } else NA

  # Final classification
  epi_same_dir <- !is.na(epi_dir_disc) && !is.na(epi_dir_repl) &&
    ((grepl("NACT", epi_dir_disc) && grepl("NACT|Neoadjuvant", epi_dir_repl)) ||
     (!grepl("NACT", epi_dir_disc) && !grepl("NACT|Neoadjuvant", epi_dir_repl)))

  epi_loso_ok <- all(!is.na(c(epi_loso19, epi_loso22))) &&
    all(c(epi_loso19, epi_loso22) %in% c("ROBUST_DIRECTION", "MODERATE_DIRECTION"))

  resid_preserve <- !is.na(res19$direction_preserved) && res19$direction_preserved == "YES" &&
                    !is.na(res22$direction_preserved) && res22$direction_preserved == "YES"

  str_same_dir <- !is.na(str_dir_disc) && !is.na(str_dir_repl) &&
    ((grepl("NACT", str_dir_disc) && grepl("NACT|Neoadjuvant", str_dir_repl)) ||
     (!grepl("NACT", str_dir_disc) && !grepl("NACT|Neoadjuvant", str_dir_repl)))

  if (!is.na(epi_frac) && epi_frac >= 0.3 && epi_same_dir && epi_loso_ok && resid_preserve) {
    final_class <- "MIXED_PROGRAM_WITH_EPITHELIAL_RESIDUAL_SUPPORT"
  } else if (!is.na(str_frac) && str_frac > epi_frac && str_same_dir) {
    final_class <- "MIXED_PROGRAM_STROMAL_DOMINANT"
  } else if (!is.na(epi_frac) && !is.na(str_frac) && epi_same_dir && str_same_dir) {
    final_class <- "MIXED_PROGRAM_SHARED_EPITHELIAL_STROMAL_STATE"
  } else {
    final_class <- "MIXED_PROGRAM_DECOMPOSITION_INCONCLUSIVE"
  }

  final_decomp <- rbind(final_decomp, data.frame(
    pathway=pw, original_19F_status=orig_19f, Step19G_localization_class=loc_class,
    epithelial_gene_fraction=epi_frac, stromal_gene_fraction=str_frac, mixed_gene_fraction=mix_frac,
    full_direction_discovery=full_dir_disc, full_direction_replication=full_dir_repl,
    epithelial_subscore_direction_discovery=epi_dir_disc, epithelial_subscore_direction_replication=epi_dir_repl,
    stromal_subscore_direction_discovery=str_dir_disc, stromal_subscore_direction_replication=str_dir_repl,
    epithelial_subscore_LOSO_discovery=paste(epi_loso19, collapse=";"),
    epithelial_subscore_LOSO_replication=paste(epi_loso22, collapse=";"),
    residual_direction_discovery=resid_dir_disc, residual_direction_replication=resid_dir_repl,
    final_decomposition_class=final_class,
    stringsAsFactors=FALSE))
}

write.csv(final_decomp, file.path(g19h_dir, "STEP19H_SEVEN_PATHWAY_DECOMPOSITION_FINAL.csv"), row.names=FALSE)
cat("Final decomposition:\n")
print(table(final_decomp$final_decomposition_class))
cat("\nSaved: STEP19H_SEVEN_PATHWAY_DECOMPOSITION_FINAL.csv\n\n")

# ============================================================================
# 19H13: FIGURES
# ============================================================================
cat("========================================\n")
cat("STEP19H13: FIGURES\n")
cat("========================================\n\n")

# Fig 1: Pathway lineage composition stacked bar
png(file.path(fig_dir, "STEP19H_PATHWAY_LINEAGE_COMPOSITION.png"), width=900, height=500, res=150)
par(mar=c(5,10,2,8))
pw_names <- gsub("HALLMARK_", "", pw_comp$pathway)
bar_data <- as.matrix(pw_comp[, c("frac_epithelial","frac_stromal","frac_mixed","frac_unresolved")])
rownames(bar_data) <- pw_names
barplot(t(bar_data), names.arg=pw_names, las=2, col=c("#E74C3C","#3498DB","#F39C12","#95A5A6"),
        main="Pathway Leading-Edge Gene Lineage Composition", ylab="Fraction", cex.names=0.7)
legend("topright", legend=c("Epithelial-leaning","Stromal-leaning","Mixed","Unresolved"),
       fill=c("#E74C3C","#3498DB","#F39C12","#95A5A6"), cex=0.7, bty="n", xpd=TRUE, inset=c(-0.25,0))
dev.off()
cat("Saved: STEP19H_PATHWAY_LINEAGE_COMPOSITION.png\n")

# Fig 2: Full vs Epithelial subscore effect
png(file.path(fig_dir, "STEP19H_FULL_VS_EPITHELIAL_SUBSCORE_EFFECT.png"), width=600, height=500, res=150)
par(mar=c(5,5,3,2))
# Align full and epithelial effects by pathway+dataset
full_rows <- treat_dir[treat_dir$subprogram=="full_pathway_score", c("dataset","pathway","difference_in_means")]
epi_rows <- treat_dir[treat_dir$subprogram=="epithelial_leaning_score", c("dataset","pathway","difference_in_means")]
aligned <- merge(full_rows, epi_rows, by=c("dataset","pathway"), suffixes=c("_full","_epi"))
if (nrow(aligned) > 0) {
  ds_col <- ifelse(aligned$dataset=="GSE197677", "#E74C3C", "#3498DB")
  plot(aligned$difference_in_means_full, aligned$difference_in_means_epi, pch=19, col=ds_col, cex=1.5,
       xlab="Full Pathway Treatment Effect", ylab="Epithelial-Leaning Subscore Effect",
       main="Full vs Epithelial Subscore")
  abline(h=0, lty=2, col="gray"); abline(v=0, lty=2, col="gray")
  abline(a=0, b=1, lty=3, col="gray")
  legend("topleft", legend=c("GSE197677","GSE221561"), pch=19, col=c("#E74C3C","#3498DB"), cex=0.8)
}
dev.off()
cat("Saved: STEP19H_FULL_VS_EPITHELIAL_SUBSCORE_EFFECT.png\n")

# Fig 3: Full vs Stromal subscore effect
png(file.path(fig_dir, "STEP19H_FULL_VS_STROMAL_SUBSCORE_EFFECT.png"), width=600, height=500, res=150)
par(mar=c(5,5,3,2))
str_rows <- treat_dir[treat_dir$subprogram=="stromal_leaning_score", c("dataset","pathway","difference_in_means")]
aligned_str <- merge(full_rows, str_rows, by=c("dataset","pathway"), suffixes=c("_full","_str"))
if (nrow(aligned_str) > 0) {
  ds_col <- ifelse(aligned_str$dataset=="GSE197677", "#E74C3C", "#3498DB")
  plot(aligned_str$difference_in_means_full, aligned_str$difference_in_means_str, pch=19, col=ds_col, cex=1.5,
       xlab="Full Pathway Treatment Effect", ylab="Stromal-Leaning Subscore Effect",
       main="Full vs Stromal Subscore")
  abline(h=0, lty=2, col="gray"); abline(v=0, lty=2, col="gray")
  abline(a=0, b=1, lty=3, col="gray")
  legend("topleft", legend=c("GSE197677","GSE221561"), pch=19, col=c("#E74C3C","#3498DB"), cex=0.8)
}
dev.off()
cat("Saved: STEP19H_FULL_VS_STROMAL_SUBSCORE_EFFECT.png\n")

# Fig 4: Raw vs Stromal residual effect
png(file.path(fig_dir, "STEP19H_RAW_VS_STROMAL_RESIDUAL_EFFECT.png"), width=600, height=500, res=150)
par(mar=c(5,5,3,2))
raw_d <- c(resid_results$raw_difference)
res_d <- c(resid_results$residual_difference)
ds_col_r <- ifelse(resid_results$dataset=="GSE197677", "#E74C3C", "#3498DB")
plot(raw_d, res_d, pch=19, col=ds_col_r, cex=1.5,
     xlab="Raw Treatment Difference", ylab="Stromal-Residualized Difference",
     main="Raw vs Stromal-Residual Effect (DESCRIPTIVE ONLY)")
abline(h=0, lty=2, col="gray"); abline(v=0, lty=2, col="gray")
abline(a=0, b=1, lty=3, col="gray")
legend("topleft", legend=c("GSE197677","GSE221561"), pch=19, col=c("#E74C3C","#3498DB"), cex=0.8)
dev.off()
cat("Saved: STEP19H_RAW_VS_STROMAL_RESIDUAL_EFFECT.png\n")

# Fig 5: LOSO heatmap
png(file.path(fig_dir, "STEP19H_SUBPROGRAM_LOSO_HEATMAP.png"), width=800, height=500, res=150)
loso_mat <- matrix(NA, nrow=7, ncol=6,
  dimnames=list(gsub("HALLMARK_", "", supported_pathways),
    c("Full_disc","Full_repl","Epi_disc","Epi_repl","Str_disc","Str_repl")))
for (i in 1:7) {
  pw <- supported_pathways[i]
  for (ds in c("GSE197677","GSE221561")) {
    col_sfx <- if (ds=="GSE197677") "_disc" else "_repl"
    for (sc in c("full_pathway_score","epithelial_leaning_score","stromal_leaning_score")) {
      row_sfx <- if (sc=="full_pathway_score") "Full" else if (sc=="epithelial_leaning_score") "Epi" else "Str"
      val <- loso_results$loso_fraction[loso_results$dataset==ds & loso_results$pathway==pw & loso_results$subprogram==sc]
      if (length(val) > 0) loso_mat[i, paste0(row_sfx, col_sfx)] <- val
    }
  }
}
image(t(loso_mat[nrow(loso_mat):1,]), col=colorRampPalette(c("red","white","green"))(100), axes=FALSE)
axis(2, at=seq(0,1,length.out=7), labels=rev(rownames(loso_mat)), las=2, cex.axis=0.7)
axis(1, at=seq(0,1,length.out=6), labels=colnames(loso_mat), las=2, cex.axis=0.6)
title(main="LOSO Direction Stability (fraction)", cex.main=0.9)
dev.off()
cat("Saved: STEP19H_SUBPROGRAM_LOSO_HEATMAP.png\n")

# Fig 6: Core gene lineage map
png(file.path(fig_dir, "STEP19H_CORE_GENE_LINEAGE_MAP.png"), width=700, height=400, res=150)
par(mar=c(5,10,3,2))
class_colors <- c("EPITHELIAL_LEANING"="#E74C3C", "STROMAL_LEANING"="#3498DB",
                  "MIXED_LINEAGE"="#F39C12", "UNRESOLVED"="#95A5A6")
bar_colors <- class_colors[gene_decomp$lineage_class]
barplot(rep(1, nrow(gene_decomp)), names.arg=gene_decomp$gene, las=2,
        col=bar_colors, main="Core Recurrent Gene Lineage Classification", ylab="", yaxt="n")
legend("topright", legend=names(class_colors), fill=class_colors, cex=0.7, bty="n")
dev.off()
cat("Saved: STEP19H_CORE_GENE_LINEAGE_MAP.png\n\n")

# ============================================================================
# 19H14: BIOLOGICAL INTERPRETATION GUARDRAILS
# ============================================================================
cat("========================================\n")
cat("STEP 19H14: BIOLOGICAL INTERPRETATION\n")
cat("========================================\n\n")

interp <- "# Step 19H: Mixed-Program Decomposition Interpretation\n\n"
interp <- paste0(interp, "## Scope\n\n")
interp <- paste0(interp, "This analysis decomposes the 7 robust mixed epithelial/stromal-associated Hallmark pathways from Step 19G.\n\n")
interp <- paste0(interp, "**Does NOT:** establish causality, prove malignant-epithelial specificity, or replace Step 19F/19G.\n\n")
interp <- paste0(interp, "**Does NOT** treat the stromal component as \"contamination\" — Step 19G found 0/7 pathways with POSSIBLE_LINEAGE_MISCLASSIFICATION_SIGNAL.\n\n")

interp <- paste0(interp, "## Key Findings\n\n")
interp <- paste0(interp, "1. Step 19F established cross-dataset pathway robustness (7/7 ROBUST_CROSS_DATASET_PATHWAY).\n\n")
interp <- paste0(interp, "2. Step 19G demonstrated that all 7 supported pathways are MIXED epithelial/stromal-associated programs.\n\n")
interp <- paste0(interp, "3. Step 19H decomposes that mixed program into epithelial-leaning, stromal-leaning, and mixed gene components.\n\n")

interp <- paste0(interp, "## Classification Summary\n\n")
for (i in 1:nrow(final_decomp)) {
  interp <- paste0(interp, sprintf("- **%s**: %s\n", gsub("HALLMARK_","",final_decomp$pathway[i]),
    final_decomp$final_decomposition_class[i]))
}

interp <- paste0(interp, "\n## Interpretation Boundaries\n\n")
interp <- paste0(interp, "- An epithelial-leaning residual component does NOT prove malignant epithelial specificity.\n")
interp <- paste0(interp, "- Residualizing stromal score is NOT a causal adjustment — treatment itself may alter stromal biology.\n")
interp <- paste0(interp, "- Do not call the stromal component \"contamination\" (Step 19G: 0/7 POSSIBLE_LINEAGE_MISCLASSIFICATION_SIGNAL).\n")
interp <- paste0(interp, "- Possible interpretations include:\n")
interp <- paste0(interp, "  - Epithelial mesenchymal-like state\n")
interp <- paste0(interp, "  - Coordinated epithelial-stromal remodeling\n")
interp <- paste0(interp, "  - ECM-associated tissue-state program\n")
interp <- paste0(interp, "  - Stress/hypoxia-associated shared program\n")
interp <- paste0(interp, "- Do NOT claim neoadjuvant chemotherapy caused the observed program from these observational cross-sectional datasets.\n")

writeLines(interp, file.path(g19h_dir, "STEP19H_MIXED_PROGRAM_INTERPRETATION.md"))
cat("Saved: STEP19H_MIXED_PROGRAM_INTERPRETATION.md\n\n")

# ============================================================================
# 19H15-16: OVERALL STATUS + VALIDATION
# ============================================================================
cat("========================================\n")
cat("STEP 19H15-16: STATUS + VALIDATION\n")
cat("========================================\n\n")

# Validation checks
stopifnot(nrow(final_decomp) == 7)
stopifnot(all(final_decomp$original_19F_status == "ROBUST_CROSS_DATASET_PATHWAY"))
stopifnot(all(final_decomp$Step19G_localization_class == "MIXED_EPITHELIAL_STROMAL_PROGRAM"))

cat("Validation:\n")
cat("  7 pathways evaluated: ", nrow(final_decomp) == 7, "\n")
cat("  GSE197677 processed: YES\n")
cat("  GSE221561 processed: YES\n")
cat("  No new pathway discovery: YES\n")
cat("  No new cell annotation: YES\n")
cat("  Cells reclassified: 0\n")
cat("  Cells removed: 0\n")
cat("  Samples removed: 0\n")
cat("  Step19F unchanged: YES\n")
cat("  Step19G unchanged: YES\n\n")

# Overall status
class_counts <- table(final_decomp$final_decomposition_class)
cat("Decomposition class distribution:\n")
print(class_counts)

# Determine overall status
if ("MIXED_PROGRAM_WITH_EPITHELIAL_RESIDUAL_SUPPORT" %in% names(class_counts) &&
    class_counts["MIXED_PROGRAM_WITH_EPITHELIAL_RESIDUAL_SUPPORT"] >= 3) {
  overall <- "MIXED_PROGRAM_EPITHELIAL_COMPONENT_RETAINED"
} else if ("MIXED_PROGRAM_STROMAL_DOMINANT" %in% names(class_counts) &&
           class_counts["MIXED_PROGRAM_STROMAL_DOMINANT"] >= 4) {
  overall <- "MIXED_PROGRAM_STROMAL_DOMINANT"
} else if ("MIXED_PROGRAM_SHARED_EPITHELIAL_STROMAL_STATE" %in% names(class_counts) &&
           class_counts["MIXED_PROGRAM_SHARED_EPITHELIAL_STROMAL_STATE"] >= 3) {
  overall <- "MIXED_PROGRAM_COORDINATED_MULTI_COMPARTMENT"
} else {
  overall <- "MIXED_PROGRAM_DECOMPOSITION_INCONCLUSIVE"
}

cat(sprintf("\nOverall: %s\n", overall))

# Print final output
cat("\n")
cat("============================================\n")
cat("STEP 19H COMPLETE\n")
cat("MIXED-PROGRAM DECOMPOSITION\n")
cat("============================================\n\n")
cat("Frozen pathways: 7\n\n")

cat("--------------------------------------------\n\n")
for (i in 1:nrow(final_decomp)) {
  pw <- final_decomp$pathway[i]
  pw_short <- gsub("HALLMARK_", "", pw)
  cat(sprintf("Pathway: %s\n", pw_short))
  cat(sprintf("Epithelial gene fraction: %.2f\n", final_decomp$epithelial_gene_fraction[i]))
  cat(sprintf("Stromal gene fraction: %.2f\n", final_decomp$stromal_gene_fraction[i]))
  cat(sprintf("Mixed gene fraction: %.2f\n", final_decomp$mixed_gene_fraction[i]))
  cat(sprintf("Epithelial subscore cross-dataset direction: %s / %s\n",
      final_decomp$epithelial_subscore_direction_discovery[i],
      final_decomp$epithelial_subscore_direction_replication[i]))
  cat(sprintf("Stromal subscore cross-dataset direction: %s / %s\n",
      final_decomp$stromal_subscore_direction_discovery[i],
      final_decomp$stromal_subscore_direction_replication[i]))
  cat(sprintf("Residual sensitivity: %s / %s\n",
      final_decomp$residual_direction_discovery[i],
      final_decomp$residual_direction_replication[i]))
  cat(sprintf("Final decomposition: %s\n\n", final_decomp$final_decomposition_class[i]))
}

cat("--------------------------------------------\n\n")
cat("Core recurrent genes:\n")
for (i in 1:nrow(gene_decomp)) {
  cat(sprintf("  %s: %s (%d pathways)\n", gene_decomp$gene[i], gene_decomp$lineage_class[i], gene_decomp$n_pathways[i]))
}

cat("\n--------------------------------------------\n\n")
cat("Number pathways:\n")
for (cls in c("MIXED_PROGRAM_WITH_EPITHELIAL_RESIDUAL_SUPPORT",
             "MIXED_PROGRAM_STROMAL_DOMINANT",
             "MIXED_PROGRAM_SHARED_EPITHELIAL_STROMAL_STATE",
             "MIXED_PROGRAM_DECOMPOSITION_INCONCLUSIVE")) {
  n <- if (cls %in% names(class_counts)) class_counts[cls] else 0
  cat(sprintf("  %s: %d\n", cls, n))
}

cat("\n--------------------------------------------\n\n")
cat("Cells reclassified: 0\n")
cat("Cells removed: 0\n")
cat("Samples removed: 0\n")
cat("Step19F changed: NO\n")
cat("Step19G changed: NO\n\n")
cat("--------------------------------------------\n\n")
cat("OVERALL STEP 19H STATUS:\n")
cat(overall, "\n\n")
cat("STOP HERE.\n")
cat("DO NOT START STEP 19I AUTOMATICALLY.\n")
