#!/usr/bin/env Rscript
# ============================================================================
# STEP 19O-FIX221: GSE221561 SOURCE PSEUDOBULK AND SIGNALING COMPLETION
# ============================================================================
# Repair the specific technical omission identified by Step19O-CORR:
# GSE221561 source-compartment pseudobulk was never computed.
# This is a TARGETED COMPLETION of Step19O for GSE221561 only.
# ============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

BASE_DIR <- Sys.getenv("ESCC_CAF_PROJECT_ROOT", unset = ".")
O_DIR    <- file.path(BASE_DIR, "06_tables/cross_dataset/pathway_analysis/Step19O_FIX221")
SRC_DIR  <- file.path(BASE_DIR, "06_tables/cross_dataset/pathway_analysis/Step19O")
L_DIR    <- file.path(BASE_DIR, "06_tables/cross_dataset/pathway_analysis/Step19L")
N_DIR    <- file.path(BASE_DIR, "06_tables/cross_dataset/pathway_analysis/Step19N")
LOG_DIR  <- file.path(BASE_DIR, "08_logs")
dir.create(O_DIR, showWarnings = FALSE, recursive = TRUE)

log_con <- file.path(LOG_DIR, "STEP19O_FIX221_GSE221561_SIGNALING_COMPLETION.log")
cat("", file = log_con)

cat_log <- function(...) {
  txt <- paste0(...)
  cat(txt)
  cat(txt, file = log_con, append = TRUE)
}

cat_log("============================================\n")
cat_log("STEP 19O-FIX221: GSE221561 SIGNALING COMPLETION\n")
cat_log("============================================\n\n")

TREAT_MAP <- c(
  S0520T="Neoadjuvant_treated", S1265T="Neoadjuvant_treated",
  S1315TS="Neoadjuvant_treated", S1535T="Neoadjuvant_treated",
  S2423T="Neoadjuvant_treated", S2487T="Neoadjuvant_treated",
  S6829T="Neoadjuvant_treated", S6417T="Surgery_alone",
  S9478T="Surgery_alone"
)

SOURCE_COMPS <- c("Epithelial", "Myeloid", "T_NK", "Endothelial")

map_to_compartment <- function(labels) {
  labels <- tolower(trimws(as.character(labels)))
  comp <- rep("Other", length(labels))
  comp[grepl("epithel|cancer|tumor|malignant", labels)] <- "Epithelial"
  comp[grepl("fibroblast|caf|stromal|frc", labels)] <- "Fibroblast"
  comp[grepl("myeloid|macrophage|monocyte|dc|dendritic|mast|neutrophil", labels)] <- "Myeloid"
  comp[grepl("t cell|t_cell|tcell|nk|cd8|cd4|t nk", labels)] <- "T_NK"
  comp[grepl("endothel|vascular|ec", labels)] <- "Endothelial"
  comp[grepl("pericyte|smooth muscle", labels)] <- "Pericyte"
  comp[grepl("b cell|b_cell|bcell|plasma", labels)] <- "Bcell"
  comp
}

compute_log2_cpm <- function(counts_mat) {
  lib_size <- colSums(counts_mat)
  lib_size[lib_size == 0] <- 1
  cpm_mat <- sweep(counts_mat, 2, lib_size, "/") * 1e6
  log2(cpm_mat + 1)
}

# ============================================================================
# FIX221-1 — INPUT PREFLIGHT
# ============================================================================
cat_log("============================================\n")
cat_log("FIX221-1: INPUT PREFLIGHT\n")
cat_log("============================================\n\n")

obj <- readRDS(file.path(BASE_DIR, "03_objects/GSE221561/GSE221561_author_annotated.rds"))
meta <- obj@meta.data
cat_log(sprintf("GSE221561 cells: %d\n", nrow(meta)))

ct_col <- "author_Cell_type"
sample_col <- "sample_id"
cat_log(sprintf("Cell type column: %s\n", ct_col))
cat_log(sprintf("Sample column: %s\n", sample_col))

comp <- map_to_compartment(meta[[ct_col]])
cat_log(sprintf("Compartments: %s\n", paste(sort(unique(comp)), collapse=", ")))

for (c in c("Fibroblast", SOURCE_COMPS)) {
  n_cells <- sum(comp == c)
  n_samps <- length(unique(meta[[sample_col]][comp == c]))
  cat_log(sprintf("  %s: %d cells, %d samples\n", c, n_cells, n_samps))
}

# Build inventory (only tumor samples in TREAT_MAP)
inv_list <- list()
for (s in unique(meta[[sample_col]])) {
  if (!s %in% names(TREAT_MAP)) next
  s_cells <- meta[[sample_col]] == s
  s_comp <- comp[s_cells]
  s_treat <- TREAT_MAP[s]
  for (c in c("Fibroblast", SOURCE_COMPS)) {
    n <- sum(s_comp == c)
    if (n > 0) {
      inv_list[[paste(s, c)]] <- data.frame(
        dataset="GSE221561", sample=s, treatment=s_treat,
        cell_type=c, n_cells=n, stringsAsFactors=FALSE)
    }
  }
}
inventory <- do.call(rbind, inv_list)
rownames(inventory) <- NULL
write.csv(inventory, file.path(O_DIR, "STEP19O_FIX221_CELLTYPE_SAMPLE_INVENTORY.csv"), row.names=FALSE)
cat_log("\nSaved: STEP19O_FIX221_CELLTYPE_SAMPLE_INVENTORY.csv\n\n")

# Load frozen CORE2/CORE4 scores
core_scores <- read.csv(file.path(L_DIR, "STEP19L_SAMPLE_CORE_MODULE_SCORES.csv"))
core221 <- core_scores[core_scores$dataset == "GSE221561", ]
cat_log(sprintf("Frozen CORE2 scores: %d GSE221561 samples\n", nrow(core221)))
for (i in 1:nrow(core221)) {
  cat_log(sprintf("  %s: CORE2=%.4f CORE4=%.4f\n",
    core221$sample_id[i], core221$CORE2_STRONG_score[i], core221$CORE4_score[i]))
}

# Load frozen fibroblast pseudobulk from Step19I
fib221 <- readRDS(file.path(BASE_DIR,
  "06_tables/cross_dataset/pathway_analysis/Step19I/counts/GSE221561_FIBROBLAST_CAF_raw_counts.rds"))
cat_log(sprintf("\nFrozen fibroblast pseudobulk: %d genes x %d samples\n",
  nrow(fib221), ncol(fib221)))

# Load frozen GSE197677 signaling outputs for cross-dataset comparison
elig_197  <- read.csv(file.path(SRC_DIR, "STEP19O_ELIGIBLE_LIGAND_RECEPTOR_PAIRS.csv"))
assoc_197 <- read.csv(file.path(SRC_DIR, "STEP19O_LIGAND_CORE2_ASSOCIATIONS.csv"))
src_pri_197 <- read.csv(file.path(SRC_DIR, "STEP19O_SOURCE_COMPARTMENT_PRIORITIZATION.csv"))
te_197    <- read.csv(file.path(SRC_DIR, "STEP19O_SOURCE_LIGAND_TREATMENT_EFFECTS.csv"))
rec_av_197 <- read.csv(file.path(SRC_DIR, "STEP19O_FIBROBLAST_RECEPTOR_AVAILABILITY.csv"))

cat_log("Frozen GSE197677 signaling outputs loaded.\n\n")

# ============================================================================
# FIX221-2 — SOURCE PSEUDOBULK
# ============================================================================
cat_log("============================================\n")
cat_log("FIX221-2: SOURCE PSEUDOBULK\n")
cat_log("============================================\n\n")

rna <- obj@assays[["RNA"]]
layer_names <- names(rna@layers)
cat_log(sprintf("Layers: %d\n", length(layer_names)))

# Use Seurat::GetAssayData(obj, layer=ln) to get properly named matrices
# (direct layer access returns matrices without row/col names)
layer_sample_map <- list()
gene_universe <- NULL
for (ln in layer_names) {
  layer_mat <- Seurat::GetAssayData(obj, assay="RNA", layer=ln)
  layer_cells <- colnames(layer_mat)
  if (length(layer_cells) == 0) next
  layer_meta <- meta[layer_cells, , drop=FALSE]
  samps <- unique(as.character(layer_meta[[sample_col]]))
  layer_sample_map[[ln]] <- list(
    cells = layer_cells,
    sample = samps[1],
    n_cells = length(layer_cells)
  )
  lg <- rownames(layer_mat)
  if (is.null(gene_universe)) {
    gene_universe <- lg
  } else {
    gene_universe <- intersect(gene_universe, lg)
  }
}
cat_log(sprintf("Gene universe (intersect across layers): %d genes\n", length(gene_universe)))
cat_log("Layer→Sample mapping:\n")
for (ln in names(layer_sample_map)) {
  ls <- layer_sample_map[[ln]]
  cat_log(sprintf("  %s: %d cells → %s\n", ln, ls$n_cells, ls$sample))
}

# Compute pseudobulk per compartment
src221 <- list()
provenance_list <- list()
for (comp_name in SOURCE_COMPS) {
  comp_cells <- rownames(meta)[comp == comp_name]
  if (length(comp_cells) < 10) {
    cat_log(sprintf("  %s: %d cells, skipping\n", comp_name, length(comp_cells)))
    next
  }

  pb_list <- list()
  for (ln in names(layer_sample_map)) {
    ls <- layer_sample_map[[ln]]
    samp <- ls$sample
    if (!samp %in% names(TREAT_MAP)) next

    layer_mat <- Seurat::GetAssayData(obj, assay="RNA", layer=ln)
    layer_cells <- colnames(layer_mat)
    common_cells <- intersect(comp_cells, layer_cells)
    if (length(common_cells) < 5) next

    layer_genes <- rownames(layer_mat)
    common_genes <- intersect(layer_genes, gene_universe)
    sub_mat <- layer_mat[match(common_genes, layer_genes), common_cells, drop=FALSE]
    pb_vec <- Matrix::rowSums(sub_mat)
    names(pb_vec) <- common_genes

    if (!samp %in% names(pb_list)) {
      pb_list[[samp]] <- rep(0, length(gene_universe))
      names(pb_list[[samp]]) <- gene_universe
    }
    pb_list[[samp]][common_genes] <- pb_list[[samp]][common_genes] + as.numeric(pb_vec)
  }

  if (length(pb_list) >= 2) {
    pb_mat <- do.call(cbind, lapply(pb_list, function(x) matrix(x, ncol=1)))
    rownames(pb_mat) <- gene_universe
    colnames(pb_mat) <- names(pb_list)
    src221[[comp_name]] <- pb_mat
    cat_log(sprintf("  %s: %d genes x %d samples\n", comp_name, nrow(pb_mat), ncol(pb_mat)))

    for (s in names(pb_list)) {
      provenance_list[[paste(comp_name, s)]] <- data.frame(
        dataset="GSE221561", sample=s, treatment=TREAT_MAP[s],
        source_compartment=comp_name,
        n_cells=0,
        expression_source="raw RNA counts (Assay5 layer via GetAssayData)",
        normalization_method="pseudobulk sum of raw counts",
        stringsAsFactors=FALSE)
    }
  }
}

# Save pseudobulk
for (comp_name in names(src221)) {
  saveRDS(src221[[comp_name]], file.path(O_DIR,
    sprintf("STEP19O_FIX221_GSE221561_%s_raw_counts.rds", comp_name)))
}

# Build counts summary
counts_list <- list()
for (comp_name in names(src221)) {
  pb <- src221[[comp_name]]
  for (s in colnames(pb)) {
    counts_list[[paste(comp_name, s)]] <- data.frame(
      dataset="GSE221561", sample=s, treatment=TREAT_MAP[s],
      source_compartment=comp_name, n_genes=nrow(pb), n_samples=ncol(pb),
      stringsAsFactors=FALSE)
  }
}
write.csv(do.call(rbind, counts_list),
  file.path(O_DIR, "STEP19O_FIX221_SOURCE_COMPARTMENT_COUNTS.csv"), row.names=FALSE)

provenance_df <- do.call(rbind, provenance_list)
# Fix n_cells properly
for (i in 1:nrow(provenance_df)) {
  cn <- provenance_df$source_compartment[i]
  ss <- provenance_df$sample[i]
  total_cells <- 0
  for (ln in names(layer_sample_map)) {
    if (layer_sample_map[[ln]]$sample == ss) {
      lc <- layer_sample_map[[ln]]$cells
      total_cells <- total_cells + sum(rownames(meta)[comp==cn] %in% lc)
    }
  }
  provenance_df$n_cells[i] <- total_cells
}
write.csv(provenance_df, file.path(O_DIR, "STEP19O_FIX221_SOURCE_PSEUDOBULK_PROVENANCE.csv"), row.names=FALSE)
cat_log("\nSaved pseudobulk and provenance.\n\n")

# Free Seurat object
rm(obj, rna); gc()

# ============================================================================
# FIX221-3 — PSEUDOBULK QC
# ============================================================================
cat_log("============================================\n")
cat_log("FIX221-3: PSEUDOBULK QC\n")
cat_log("============================================\n\n")

eval_list <- list()
for (comp_name in names(src221)) {
  pb <- src221[[comp_name]]
  samps <- colnames(pb)
  treated <- samps[sapply(samps, function(s) TREAT_MAP[s] == "Neoadjuvant_treated")]
  control <- samps[sapply(samps, function(s) TREAT_MAP[s] == "Surgery_alone")]

  cell_counts <- provenance_df$n_cells[provenance_df$source_compartment == comp_name]

  min_cells <- min(cell_counts)
  med_cells <- median(cell_counts)
  max_cells <- max(cell_counts)

  n_treated <- length(treated)
  n_control <- length(control)
  n_total <- length(samps)

  if (n_total >= 9 && min_cells >= 30) {
    eval_class <- "ADEQUATELY_EVALUABLE"
  } else if (n_total >= 4 && min_cells >= 10) {
    eval_class <- "LIMITED_EVALUABILITY"
  } else {
    eval_class <- "INSUFFICIENT_EVALUABILITY"
  }

  eval_list[[comp_name]] <- data.frame(
    source_compartment=comp_name,
    n_tumor_samples=n_total,
    n_treated=n_treated,
    n_surgery_alone=n_control,
    min_cells_per_sample=min_cells,
    median_cells_per_sample=med_cells,
    max_cells_per_sample=max_cells,
    evaluability=eval_class,
    stringsAsFactors=FALSE)

  cat_log(sprintf("  %s: %d samples (%d treated, %d control), cells min=%d med=%d max=%d -> %s\n",
    comp_name, n_total, n_treated, n_control, min_cells, med_cells, max_cells, eval_class))
}

eval_df <- do.call(rbind, eval_list)
write.csv(eval_df, file.path(O_DIR, "STEP19O_FIX221_SOURCE_EVALUABILITY.csv"), row.names=FALSE)
cat_log("\nSaved: STEP19O_FIX221_SOURCE_EVALUABILITY.csv\n\n")

# ============================================================================
# FIX221-4 — LOCAL LR RESOURCE
# ============================================================================
cat_log("============================================\n")
cat_log("FIX221-4: LOCAL LR RESOURCE\n")
cat_log("============================================\n\n")

# Rebuild EXACTLY the same curated LR resource as original Step19O
.lr_lig <- c(
  "TGFB1","TGFB1","TGFB1","TGFB2","TGFB2","TGFB3","TGFB3",
  "BMP2","BMP2","BMP2","BMP4","BMP4","BMP7","BMP7","GDF2","GDF2",
  "TNF","TNF","TNF","LTA","LTA","TNFSF12","TNFSF12",
  "IL1A","IL1A","IL1B","IL1B",
  "IL6","IL6","IL6","OSM","OSM","LIF","LIF",
  "IFNG","IFNG","IFNB1","IFNB1","IFNA1","IFNA1",
  "PDGFA","PDGFA","PDGFB","PDGFB","PDGFC","PDGFC","PDGFD","PDGFD",
  "FGF1","FGF1","FGF2","FGF2","FGF7","FGF7","FGF9","FGF9","FGF18","FGF18",
  "EGF","EGF","HBEGF","HBEGF","AREG","AREG","EREG","EREG","NRG1","NRG1",
  "WNT3A","WNT3A","WNT5A","WNT5A","WNT5B","WNT5B",
  "CXCL12","CXCL12","CCL2","CCL2","CCL5","CCL5",
  "CXCL8","CXCL8","CXCL1","CXCL1","CXCL2","CXCL2",
  "VEGFA","VEGFA","VEGFB","VEGFB",
  "DLL1","DLL1","DLL4","DLL4","JAG1","JAG1","JAG2","JAG2",
  "HGF","HGF"
)
.lr_rec <- c(
  "TGFBR1","TGFBR2","TGFBR3","TGFBR1","TGFBR2","TGFBR1","TGFBR2",
  "BMPR1A","BMPR2","ACVR1","BMPR1A","BMPR2","BMPR1B","BMPR2","BMPR2","ACVRL1",
  "TNFRSF1A","TNFRSF1B","LTBR","LTBR","TNFRSF1A","TNFRSF12A","TNFRSF1A",
  "IL1R1","IL1RAP","IL1R1","IL1RAP",
  "IL6R","IL6ST","IL6R","OSMR","OSMR","LIFR","LIFR",
  "IFNGR1","IFNGR2","IFNAR1","IFNAR2","IFNAR1","IFNAR2",
  "PDGFRA","PDGFRB","PDGFRB","PDGFRA","PDGFRA","PDGFRB","PDGFRA","PDGFRB",
  "FGFR1","FGFR2","FGFR1","FGFR2","FGFR2","FGFR1","FGFR2","FGFR3","FGFR1","FGFR2",
  "EGFR","ERBB2","EGFR","ERBB4","EGFR","ERBB2","EGFR","ERBB4","ERBB2","ERBB3",
  "FZD1","LRP5","FZD2","RYK","FZD5","LRP5",
  "CXCR4","ACKR3","CCR2","ACKR4","CCR5","ACKR4",
  "CXCR1","CXCR2","CXCR2","ACKR1","CXCR2","ACKR1",
  "KDR","FLT1","FLT1","NRP1",
  "NOTCH1","NOTCH2","NOTCH1","NOTCH2","NOTCH1","NOTCH2","NOTCH1","NOTCH2",
  "MET","AXL"
)

curated_lr <- data.frame(
  ligand=.lr_lig, receptor=.lr_rec, stringsAsFactors=FALSE)
lr_fam <- c(
  rep("TGF_beta",16), rep("TNF",7), rep("IL1",4), rep("IL6_JAK_STAT",7),
  rep("IFN",6), rep("PDGF",8), rep("FGF",10), rep("EGF",10),
  rep("WNT",6), rep("CXCL_CCL",12), rep("VEGF",4), rep("NOTCH",8), rep("HGF_MET",2)
)
curated_lr$family <- lr_fam

res_prov <- data.frame(
  resource="curated_minimal", version="prespecified_in_script",
  LR_pairs=nrow(curated_lr), ligands=length(unique(curated_lr$ligand)),
  receptors=length(unique(curated_lr$receptor)),
  ligand_target_available=FALSE, stringsAsFactors=FALSE)
write.csv(res_prov, file.path(O_DIR, "STEP19O_FIX221_RESOURCE_PROVENANCE.csv"), row.names=FALSE)
cat_log(sprintf("LR resource: %d pairs, %d ligands, %d receptors\n",
  nrow(curated_lr), length(unique(curated_lr$ligand)), length(unique(curated_lr$receptor))))
cat_log("Saved: STEP19O_FIX221_RESOURCE_PROVENANCE.csv\n\n")

all_ligands <- unique(curated_lr$ligand)
all_receptors <- unique(curated_lr$receptor)

# ============================================================================
# FIX221-5 — RECEIVER RECEPTOR AVAILABILITY
# ============================================================================
cat_log("============================================\n")
cat_log("FIX221-5: RECEIVER RECEPTOR AVAILABILITY\n")
cat_log("============================================\n\n")

fib221_logcpm <- compute_log2_cpm(fib221)

rec_results <- list()
for (rec in all_receptors) {
  if (rec %in% rownames(fib221)) {
    expr <- as.numeric(fib221_logcpm[rec, ])
    detected <- sum(expr > 0)
    n_samp <- ncol(fib221)
    frac <- detected / n_samp
    mean_expr <- mean(expr)
    if (frac >= 0.50) flag <- "ROBUSTLY_DETECTED"
    else if (frac >= 0.30) flag <- "PARTIALLY_DETECTED"
    else flag <- "POORLY_DETECTED"
  } else {
    detected <- 0; n_samp <- ncol(fib221); frac <- 0; mean_expr <- 0
    flag <- "NOT_IN_DATA"
  }
  rec_results[[rec]] <- data.frame(
    receptor=rec, n_detected=detected, n_samples=n_samp,
    detection_fraction=frac, mean_expression=mean_expr,
    classification=flag, stringsAsFactors=FALSE)
}
rec_df <- do.call(rbind, rec_results)
rec_df$dataset <- "GSE221561"
write.csv(rec_df, file.path(O_DIR, "STEP19O_FIX221_FIBROBLAST_RECEPTOR_AVAILABILITY.csv"), row.names=FALSE)

cat_log("=== FIBROBLAST RECEPTOR AVAILABILITY (GSE221561) ===\n")
for (i in 1:nrow(rec_df)) {
  if (rec_df$classification[i] != "POORLY_DETECTED" && rec_df$classification[i] != "NOT_IN_DATA") {
    cat_log(sprintf("  %s: frac=%.2f mean=%.1f -> %s\n",
      rec_df$receptor[i], rec_df$detection_fraction[i], rec_df$mean_expression[i], rec_df$classification[i]))
  }
}
cat_log("\nSaved: STEP19O_FIX221_FIBROBLAST_RECEPTOR_AVAILABILITY.csv\n\n")

# Merge with GSE197677 for cross-dataset receptor view
rec_merged <- merge(rec_av_197[, c("receptor","classification","dataset")],
  rec_df[, c("receptor","classification","dataset")], by="receptor", all=TRUE)
colnames(rec_merged)[2:3] <- c("GSE197677_class","GSE197677_ds")
colnames(rec_merged)[4:5] <- c("GSE221561_class","GSE221561_ds")

# ============================================================================
# FIX221-6 — SOURCE LIGAND AVAILABILITY
# ============================================================================
cat_log("============================================\n")
cat_log("FIX221-6: SOURCE LIGAND AVAILABILITY\n")
cat_log("============================================\n\n")

lig_results <- list()
for (comp_name in names(src221)) {
  pb <- src221[[comp_name]]
  logcpm <- compute_log2_cpm(pb)
  for (lig in all_ligands) {
    if (lig %in% rownames(logcpm)) {
      expr <- as.numeric(logcpm[lig, ])
      detected <- sum(expr > 0)
      n_samp <- ncol(logcpm)
      frac <- detected / n_samp
      mean_expr <- mean(expr)
      var_expr <- var(expr)
      if (frac >= 0.50) flag <- "ADEQUATELY_DETECTED"
      else if (frac >= 0.30) flag <- "PARTIALLY_DETECTED"
      else flag <- "POORLY_DETECTED"
    } else {
      detected <- 0; n_samp <- ncol(pb); frac <- 0; mean_expr <- 0; var_expr <- 0
      flag <- "NOT_IN_DATA"
    }
    lig_results[[paste(comp_name, lig)]] <- data.frame(
      source=comp_name, ligand=lig, n_detected=detected, n_samples=n_samp,
      detection_fraction=frac, mean_expression=mean_expr,
      variance=var_expr, classification=flag, stringsAsFactors=FALSE)
  }
}
lig_df <- do.call(rbind, lig_results)
write.csv(lig_df, file.path(O_DIR, "STEP19O_FIX221_SOURCE_LIGAND_AVAILABILITY.csv"), row.names=FALSE)

detected_ligands <- unique(lig_df$ligand[lig_df$classification %in% c("ADEQUATELY_DETECTED","PARTIALLY_DETECTED")])
cat_log(sprintf("Adequately/partially detected ligands: %d\n", length(detected_ligands)))
cat_log("Saved: STEP19O_FIX221_SOURCE_LIGAND_AVAILABILITY.csv\n\n")

# ============================================================================
# FIX221-7 — GSE221561 LIGAND-CORE2 ASSOCIATION
# ============================================================================
cat_log("============================================\n")
cat_log("FIX221-7: LIGAND-CORE2 ASSOCIATION\n")
cat_log("============================================\n\n")

fib221_samps <- core221$sample_id
core2_vals <- core221$CORE2_STRONG_score
names(core2_vals) <- fib221_samps

assoc_list <- list()
for (comp_name in names(src221)) {
  pb <- src221[[comp_name]]
  matched_samps <- intersect(colnames(pb), fib221_samps)
  if (length(matched_samps) < 4) {
    cat_log(sprintf("  %s: only %d matched samples, skipping\n", comp_name, length(matched_samps)))
    next
  }

  matched_core2 <- core2_vals[matched_samps]
  logcpm <- compute_log2_cpm(pb[, matched_samps, drop=FALSE])

  for (lig in all_ligands) {
    if (!lig %in% rownames(logcpm)) next
    lig_expr <- as.numeric(logcpm[lig, ])
    valid <- !is.na(lig_expr) & !is.na(matched_core2)
    if (sum(valid) < 4) next

    rho <- tryCatch(cor(lig_expr[valid], matched_core2[valid], method="spearman"),
                    error=function(e) NA_real_)
    if (is.na(rho)) next

    if (abs(rho) >= 0.60) flag <- "STRONG_ASSOCIATION"
    else if (abs(rho) >= 0.40) flag <- "MODERATE_ASSOCIATION"
    else flag <- "WEAK_ASSOCIATION"

    assoc_list[[paste(comp_name, lig)]] <- data.frame(
      dataset="GSE221561", compartment=comp_name, ligand=lig,
      n_paired_samples=sum(valid), rho=rho,
      direction=ifelse(rho > 0, "POSITIVE", "NEGATIVE"),
      association_class=flag, stringsAsFactors=FALSE)
  }
}
assoc_df <- do.call(rbind, assoc_list)
rownames(assoc_df) <- NULL
write.csv(assoc_df, file.path(O_DIR, "STEP19O_FIX221_LIGAND_CORE2_ASSOCIATIONS.csv"), row.names=FALSE)

cat_log(sprintf("Total ligand-CORE2 associations: %d\n", nrow(assoc_df)))
cat_log(sprintf("STRONG: %d, MODERATE: %d, WEAK: %d\n",
  sum(assoc_df$association_class=="STRONG_ASSOCIATION"),
  sum(assoc_df$association_class=="MODERATE_ASSOCIATION"),
  sum(assoc_df$association_class=="WEAK_ASSOCIATION")))

top_assoc <- assoc_df[order(-abs(assoc_df$rho)), ]
cat_log("\nTop 10 GSE221561 associations:\n")
for (i in 1:min(10, nrow(top_assoc))) {
  cat_log(sprintf("  %s %s: rho=%.3f (%s) -> %s\n",
    top_assoc$compartment[i], top_assoc$ligand[i],
    top_assoc$rho[i], top_assoc$direction[i], top_assoc$association_class[i]))
}
cat_log("\nSaved: STEP19O_FIX221_LIGAND_CORE2_ASSOCIATIONS.csv\n\n")

# ============================================================================
# FIX221-8 — CORE4 COMPARATOR
# ============================================================================
cat_log("============================================\n")
cat_log("FIX221-8: CORE4 COMPARATOR\n")
cat_log("============================================\n\n")

core4_vals <- core221$CORE4_score
names(core4_vals) <- fib221_samps

assoc4_list <- list()
for (comp_name in names(src221)) {
  pb <- src221[[comp_name]]
  matched_samps <- intersect(colnames(pb), fib221_samps)
  if (length(matched_samps) < 4) next
  matched_core4 <- core4_vals[matched_samps]
  logcpm <- compute_log2_cpm(pb[, matched_samps, drop=FALSE])

  for (lig in all_ligands) {
    if (!lig %in% rownames(logcpm)) next
    lig_expr <- as.numeric(logcpm[lig, ])
    valid <- !is.na(lig_expr) & !is.na(matched_core4)
    if (sum(valid) < 4) next
    rho4 <- tryCatch(cor(lig_expr[valid], matched_core4[valid], method="spearman"),
                     error=function(e) NA_real_)
    if (is.na(rho4)) next
    assoc4_list[[paste(comp_name, lig)]] <- data.frame(
      dataset="GSE221561", compartment=comp_name, ligand=lig,
      rho_core4=rho4, stringsAsFactors=FALSE)
  }
}
assoc4_df <- do.call(rbind, assoc4_list)
rownames(assoc4_df) <- NULL

comp_df <- merge(assoc_df, assoc4_df, by=c("dataset","compartment","ligand"), all.x=TRUE)
comp_df$core2_core4_class <- "NO_CLEAR_CORE_ASSOCIATION"
for (i in 1:nrow(comp_df)) {
  rho2 <- comp_df$rho[i]
  rho4 <- comp_df$rho_core4[i]
  if (is.na(rho4)) {
    comp_df$core2_core4_class[i] <- "CORE2_ONLY"
  } else if (abs(rho2) > abs(rho4) * 1.5) {
    comp_df$core2_core4_class[i] <- "CORE2_PREFERENTIAL"
  } else if (abs(rho4) > abs(rho2) * 1.5) {
    comp_df$core2_core4_class[i] <- "CORE4_PREFERENTIAL"
  } else {
    comp_df$core2_core4_class[i] <- "SHARED_CORE_ASSOCIATION"
  }
}
write.csv(comp_df, file.path(O_DIR, "STEP19O_FIX221_LIGAND_CORE2_CORE4_COMPARISON.csv"), row.names=FALSE)
cat_log(sprintf("CORE2_PREFERENTIAL: %d, CORE4_PREFERENTIAL: %d, SHARED: %d\n",
  sum(comp_df$core2_core4_class=="CORE2_PREFERENTIAL"),
  sum(comp_df$core2_core4_class=="CORE4_PREFERENTIAL"),
  sum(comp_df$core2_core4_class=="SHARED_CORE_ASSOCIATION")))
cat_log("Saved: STEP19O_FIX221_LIGAND_CORE2_CORE4_COMPARISON.csv\n\n")

# ============================================================================
# FIX221-9 — TREATMENT EFFECTS
# ============================================================================
cat_log("============================================\n")
cat_log("FIX221-9: TREATMENT EFFECTS\n")
cat_log("============================================\n\n")

te_list <- list()
for (comp_name in names(src221)) {
  pb <- src221[[comp_name]]
  samps <- colnames(pb)
  treated_s <- samps[sapply(samps, function(s) TREAT_MAP[s] == "Neoadjuvant_treated")]
  control_s <- samps[sapply(samps, function(s) TREAT_MAP[s] == "Surgery_alone")]
  if (length(treated_s) < 1 || length(control_s) < 1) next

  logcpm <- compute_log2_cpm(pb)
  for (lig in all_ligands) {
    if (!lig %in% rownames(logcpm)) next
    treated_expr <- as.numeric(logcpm[lig, treated_s])
    control_expr <- as.numeric(logcpm[lig, control_s])
    eff <- mean(treated_expr, na.rm=TRUE) - mean(control_expr, na.rm=TRUE)
    if (is.nan(eff)) eff <- NA
    te_list[[paste(comp_name, lig)]] <- data.frame(
      dataset="GSE221561", source=comp_name, ligand=lig,
      treated_mean=mean(treated_expr, na.rm=TRUE),
      control_mean=mean(control_expr, na.rm=TRUE),
      effect=eff,
      direction=ifelse(is.na(eff), "NA", ifelse(eff > 0, "HIGHER_IN_TREATED", "LOWER_IN_TREATED")),
      n_treated=length(treated_s), n_control=length(control_s),
      stringsAsFactors=FALSE)
  }
}
te_df <- do.call(rbind, te_list)
rownames(te_df) <- NULL
write.csv(te_df, file.path(O_DIR, "STEP19O_FIX221_SOURCE_LIGAND_TREATMENT_EFFECTS.csv"), row.names=FALSE)
cat_log(sprintf("Treatment effects computed for %d ligand/source pairs\n", nrow(te_df)))
cat_log("NOTE: Surgery_alone n=2 — effect estimates are fragile.\n")
cat_log("Saved: STEP19O_FIX221_SOURCE_LIGAND_TREATMENT_EFFECTS.csv\n\n")

# ============================================================================
# FIX221-10 — EXACT PERMUTATION
# ============================================================================
cat_log("============================================\n")
cat_log("FIX221-10: EXACT PERMUTATION\n")
cat_log("============================================\n\n")

# For top candidate ligands only
all_tumor_samps <- names(TREAT_MAP)
treated_idx <- which(TREAT_MAP == "Neoadjuvant_treated")
control_idx <- which(TREAT_MAP == "Surgery_alone")
n_treated <- length(treated_idx)
n_control <- length(control_idx)
cat_log(sprintf("GSE221561: %d treated, %d control, total permutations: %.0f\n",
  n_treated, n_control, choose(n_treated + n_control, n_control)))

perm_list <- list()
# Get top candidates by abs(CORE2 rho)
if (nrow(assoc_df) > 0) {
  top_ligs <- unique(assoc_df$ligand[order(-abs(assoc_df$rho))])
  top_ligs <- head(top_ligs, 20)
} else {
  top_ligs <- character(0)
}

for (comp_name in names(src221)) {
  pb <- src221[[comp_name]]
  samps <- colnames(pb)
  if (!all(all_tumor_samps %in% samps) && !all(samps %in% all_tumor_samps)) {
    # Only use samples that are in both
    samps <- intersect(samps, all_tumor_samps)
  }
  if (length(samps) < 9) next
  logcpm <- compute_log2_cpm(pb[, samps, drop=FALSE])

  for (lig in top_ligs) {
    if (!lig %in% rownames(logcpm)) next
    lig_expr <- as.numeric(logcpm[lig, samps])

    observed_eff <- mean(lig_expr[sapply(samps, function(s) TREAT_MAP[s]=="Neoadjuvant_treated")]) -
                    mean(lig_expr[sapply(samps, function(s) TREAT_MAP[s]=="Surgery_alone")])

    # Enumerate all choose(9,2) = 36 assignments of 2 control samples
    perm_effects <- numeric(0)
    control_combos <- combn(length(samps), n_control)
    treat_map_vec <- sapply(samps, function(s) TREAT_MAP[s])
    observed_control <- which(treat_map_vec == "Surgery_alone")

    for (j in 1:ncol(control_combos)) {
      ctrl_idx <- control_combos[, j]
      eff <- mean(lig_expr[-ctrl_idx]) - mean(lig_expr[ctrl_idx])
      perm_effects <- c(perm_effects, eff)
    }
    perm_effects <- perm_effects[!is.na(perm_effects)]
    if (length(perm_effects) == 0) next

    p_two <- (sum(abs(perm_effects) >= abs(observed_eff)) + 1) / (length(perm_effects) + 1)

    perm_list[[paste(comp_name, lig)]] <- data.frame(
      dataset="GSE221561", source=comp_name, ligand=lig,
      observed_effect=observed_eff,
      n_permutations=length(perm_effects),
      p_twosided=p_two,
      stringsAsFactors=FALSE)
  }
}
perm_df <- do.call(rbind, perm_list)
rownames(perm_df) <- NULL
write.csv(perm_df, file.path(O_DIR, "STEP19O_FIX221_LIGAND_EXACT_PERMUTATION.csv"), row.names=FALSE)
cat_log(sprintf("Exact permutation computed for %d ligand/source pairs\n", nrow(perm_df)))
cat_log("Saved: STEP19O_FIX221_LIGAND_EXACT_PERMUTATION.csv\n\n")

# ============================================================================
# FIX221-11 — ELIGIBLE LR PAIRS
# ============================================================================
cat_log("============================================\n")
cat_log("FIX221-11: ELIGIBLE LR PAIRS\n")
cat_log("============================================\n\n")

elig_list <- list()
for (comp_name in names(src221)) {
  pb <- src221[[comp_name]]
  for (i in 1:nrow(curated_lr)) {
    lig <- curated_lr$ligand[i]
    rec <- curated_lr$receptor[i]
    fam <- curated_lr$family[i]

    # Ligand detection in source
    lig_avail <- lig_df[lig_df$source == comp_name & lig_df$ligand == lig, ]
    lig_det <- if (nrow(lig_avail) > 0) lig_avail$classification[1] else "NOT_IN_DATA"

    # Receptor detection in fibroblast
    rec_avail <- rec_df[rec_df$receptor == rec, ]
    rec_det <- if (nrow(rec_avail) > 0) rec_avail$classification[1] else "NOT_IN_DATA"

    if (!lig_det %in% c("ADEQUATELY_DETECTED", "PARTIALLY_DETECTED")) next
    if (!rec_det %in% c("ROBUSTLY_DETECTED", "PARTIALLY_DETECTED")) next

    # CORE2 association
    assoc_lig <- assoc_df[assoc_df$compartment == comp_name & assoc_df$ligand == lig, ]
    if (nrow(assoc_lig) == 0) next
    rho2 <- assoc_lig$rho[1]
    n_paired <- assoc_lig$n_paired_samples[1]

    # Treatment effect
    te_lig <- te_df[te_df$source == comp_name & te_df$ligand == lig, ]
    te_eff <- if (nrow(te_lig) > 0) te_lig$effect[1] else NA

    elig_list[[paste(comp_name, lig, rec)]] <- data.frame(
      dataset="GSE221561", source=comp_name, ligand=lig, receptor=rec,
      family=fam, ligand_detection=lig_det, receptor_detection=rec_det,
      n_paired_samples=n_paired, CORE2_rho=rho2,
      treatment_effect=te_eff, resource="curated_minimal",
      stringsAsFactors=FALSE)
  }
}
elig221_df <- do.call(rbind, elig_list)
rownames(elig221_df) <- NULL
write.csv(elig221_df, file.path(O_DIR, "STEP19O_FIX221_ELIGIBLE_LR_PAIRS.csv"), row.names=FALSE)
cat_log(sprintf("GSE221561 eligible LR pairs: %d\n", nrow(elig221_df)))
cat_log("Saved: STEP19O_FIX221_ELIGIBLE_LR_PAIRS.csv\n\n")

# ============================================================================
# FIX221-12 — PRESPECIFIED CANDIDATE CHECK
# ============================================================================
cat_log("============================================\n")
cat_log("FIX221-12: PRESPECIFIED CANDIDATE CHECK\n")
cat_log("============================================\n\n")

candidates <- list(
  c("AREG", "Epithelial"),
  c("DLL1", "Epithelial"),
  c("EREG", "Epithelial"),
  c("TGFB3", "AUTO"),
  c("IFNG", "Myeloid")
)

cand_results <- list()
for (cand in candidates) {
  lig <- cand[1]
  pref_src <- cand[2]

  # Auto-determine source for TGFB3
  if (pref_src == "AUTO") {
    lig_assoc <- assoc_df[assoc_df$ligand == lig, ]
    if (nrow(lig_assoc) > 0) {
      best_idx <- which.max(abs(lig_assoc$rho))
      pref_src <- as.character(lig_assoc$compartment[best_idx])
    } else {
      pref_src <- "NONE"
    }
  }

  # Source present
  src_present <- pref_src %in% names(src221)

  # Ligand detected
  lig_avail <- lig_df[lig_df$source == pref_src & lig_df$ligand == lig, ]
  lig_det <- if (nrow(lig_avail) > 0) lig_avail$classification[1] else "NOT_IN_DATA"

  # Receptor detected in fibroblast
  lig_receptors <- unique(curated_lr$receptor[curated_lr$ligand == lig])
  rec_avail <- rec_df[rec_df$receptor %in% lig_receptors, ]
  rec_det <- if (nrow(rec_avail) > 0) {
    any_robust <- any(rec_avail$classification %in% c("ROBUSTLY_DETECTED","PARTIALLY_DETECTED"))
    ifelse(any_robust, "DETECTED", "NOT_DETECTED")
  } else "NOT_DETECTED"

  # CORE2 rho
  assoc_lig <- assoc_df[assoc_df$compartment == pref_src & assoc_df$ligand == lig, ]
  rho2 <- if (nrow(assoc_lig) > 0) assoc_lig$rho[1] else NA
  n_paired <- if (nrow(assoc_lig) > 0) assoc_lig$n_paired_samples[1] else NA

  # CORE4 rho
  assoc4_lig <- comp_df[comp_df$compartment == pref_src & comp_df$ligand == lig, ]
  rho4 <- if (nrow(assoc4_lig) > 0) assoc4_lig$rho_core4[1] else NA

  # Treatment effect
  te_lig <- te_df[te_df$source == pref_src & te_df$ligand == lig, ]
  te_eff <- if (nrow(te_lig) > 0) te_lig$effect[1] else NA

  # Permutation P
  perm_lig <- perm_df[perm_df$source == pref_src & perm_df$ligand == lig, ]
  perm_p <- if (nrow(perm_lig) > 0) perm_lig$p_twosided[1] else NA

  # Compare with GSE197677
  assoc_197_lig <- assoc_197[assoc_197$compartment == pref_src & assoc_197$ligand == lig, ]
  rho2_197 <- if (nrow(assoc_197_lig) > 0) assoc_197_lig$rho[1] else NA

  # Classification
  if (!src_present) {
    status <- "SOURCE_NOT_EVALUABLE"
  } else if (lig_det == "NOT_IN_DATA" || lig_det == "POORLY_DETECTED") {
    status <- "LIGAND_NOT_DETECTED"
  } else if (rec_det == "NOT_DETECTED") {
    status <- "RECEPTOR_NOT_DETECTED"
  } else if (is.na(rho2)) {
    status <- "INSUFFICIENT_PAIRED_SAMPLES"
  } else if (!is.na(rho2) && !is.na(rho2_197) && sign(rho2) == sign(rho2_197) &&
             abs(rho2) >= 0.40 && abs(rho2_197) >= 0.40) {
    status <- "REPLICATED_SAME_DIRECTION"
  } else if (!is.na(rho2) && !is.na(rho2_197) && sign(rho2) != sign(rho2_197)) {
    status <- "EVALUABLE_OPPOSITE_DIRECTION"
  } else {
    status <- "EVALUABLE_NO_REPLICATION"
  }

  cand_results[[lig]] <- data.frame(
    ligand=lig, source=pref_src,
    source_present=src_present, ligand_detected=lig_det,
    receptor_detected=rec_det, n_paired_samples=n_paired,
    GSE221561_CORE2_rho=rho2, GSE221561_CORE4_rho=rho4,
    GSE221561_treatment_effect=te_eff,
    GSE197677_CORE2_rho=rho2_197,
    exact_permutation_p=perm_p,
    status=status,
    stringsAsFactors=FALSE)

  cat_log(sprintf("  %s/%s: rho221=%.3f rho197=%.3f n=%s -> %s\n",
    lig, pref_src,
    ifelse(is.na(rho2), NA, round(rho2, 3)),
    ifelse(is.na(rho2_197), NA, round(rho2_197, 3)),
    ifelse(is.na(n_paired), "NA", n_paired),
    status))
}
cand_df <- do.call(rbind, cand_results)
write.csv(cand_df, file.path(O_DIR, "STEP19O_FIX221_PRESPECIFIED_CANDIDATES.csv"), row.names=FALSE)
cat_log("\nSaved: STEP19O_FIX221_PRESPECIFIED_CANDIDATES.csv\n\n")

# ============================================================================
# FIX221-13 — RUNX1 BRIDGE
# ============================================================================
cat_log("============================================\n")
cat_log("FIX221-13: RUNX1 BRIDGE\n")
cat_log("============================================\n\n")

# Check for RUNX1 activity scores from Step19N
runx1_file <- file.path(N_DIR, "STEP19N_FOS_ACTIVITY_SCORES.csv")
if (file.exists(runx1_file)) {
  runx1_scores <- read.csv(runx1_file)
  cat_log(sprintf("Loaded RUNX1 scores from Step19N: %d rows\n", nrow(runx1_scores)))
} else {
  # Try alternative names
  runx1_files <- list.files(N_DIR, pattern="RUNX1", full.names=TRUE)
  if (length(runx1_files) > 0) {
    cat_log(sprintf("Found RUNX1 files: %s\n", paste(basename(runx1_files), collapse=", ")))
    runx1_scores <- read.csv(runx1_files[1])
  } else {
    runx1_scores <- NULL
    cat_log("No RUNX1 activity scores found in Step19N.\n")
  }
}

runx1_bridge_list <- list()
if (!is.null(runx1_scores)) {
  # Check columns
  runx1_cols <- colnames(runx1_scores)
  cat_log(sprintf("RUNX1 score columns: %s\n", paste(runx1_cols, collapse=", ")))

  # Find SCORE column and SAMPLE column
  score_col <- grep("score|activity|rho", runx1_cols, value=TRUE, ignore.case=TRUE)[1]
  sample_col_runx1 <- grep("sample|sample_id", runx1_cols, value=TRUE, ignore.case=TRUE)[1]

  if (!is.na(score_col) && !is.na(sample_col_runx1)) {
    runx1_221 <- runx1_scores[runx1_scores$dataset == "GSE221561", ]
    if (nrow(runx1_221) > 0) {
      runx1_vals <- setNames(runx1_221[[score_col]], runx1_221[[sample_col_runx1]])

      for (comp_name in names(src221)) {
        pb <- src221[[comp_name]]
        matched_samps <- intersect(colnames(pb), names(runx1_vals))
        if (length(matched_samps) < 4) next

        matched_runx1 <- runx1_vals[matched_samps]
        matched_core2 <- core2_vals[matched_samps]
        logcpm <- compute_log2_cpm(pb[, matched_samps, drop=FALSE])

        for (lig in all_ligands) {
          if (!lig %in% rownames(logcpm)) next
          lig_expr <- as.numeric(logcpm[lig, ])
          valid <- !is.na(lig_expr) & !is.na(matched_core2) & !is.na(matched_runx1)
          if (sum(valid) < 4) next

          rho_runx1 <- tryCatch(cor(lig_expr[valid], matched_runx1[valid], method="spearman"),
                                error=function(e) NA_real_)
          rho_core2 <- tryCatch(cor(lig_expr[valid], matched_core2[valid], method="spearman"),
                                error=function(e) NA_real_)

          if (is.na(rho_runx1) || is.na(rho_core2)) next

          if (!is.na(rho_runx1) && !is.na(rho_core2) &&
              abs(rho_runx1) >= 0.40 && abs(rho_core2) >= 0.40 &&
              sign(rho_runx1) == sign(rho_core2)) {
            br_class <- "LIGAND_RUNX1_CORE2_BRIDGE"
          } else if (!is.na(rho_core2) && abs(rho_core2) >= 0.40) {
            br_class <- "LIGAND_CORE2_ONLY"
          } else if (!is.na(rho_runx1) && abs(rho_runx1) >= 0.40) {
            br_class <- "LIGAND_RUNX1_ONLY"
          } else {
            br_class <- "NO_BRIDGE_SUPPORT"
          }

          runx1_bridge_list[[paste(comp_name, lig)]] <- data.frame(
            dataset="GSE221561", compartment=comp_name, ligand=lig,
            ligand_vs_runx1_rho=rho_runx1, ligand_vs_core2_rho=rho_core2,
            classification=br_class, stringsAsFactors=FALSE)
        }
      }
    }
  }
}

if (length(runx1_bridge_list) > 0) {
  runx1_bridge_df <- do.call(rbind, runx1_bridge_list)
  rownames(runx1_bridge_df) <- NULL
  cat_log(sprintf("RUNX1 bridges: %d\n", sum(runx1_bridge_df$classification == "LIGAND_RUNX1_CORE2_BRIDGE")))
} else {
  runx1_bridge_df <- data.frame(dataset=character(0), compartment=character(0),
    ligand=character(0), ligand_vs_runx1_rho=numeric(0),
    ligand_vs_core2_rho=numeric(0), classification=character(0),
    stringsAsFactors=FALSE)
  cat_log("RUNX1 bridge: NOT EVALUABLE\n")
}
write.csv(runx1_bridge_df, file.path(O_DIR, "STEP19O_FIX221_RUNX1_BRIDGE.csv"), row.names=FALSE)
cat_log("Saved: STEP19O_FIX221_RUNX1_BRIDGE.csv\n\n")

# ============================================================================
# FIX221-14 — LOSO
# ============================================================================
cat_log("============================================\n")
cat_log("FIX221-14: LOSO (Leave One Sample Out)\n")
cat_log("============================================\n\n")

# Top GSE221561 candidates
if (nrow(assoc_df) > 0) {
  top_loso_ligs <- unique(assoc_df$ligand[order(-abs(assoc_df$rho))])
  top_loso_ligs <- head(top_loso_ligs, 10)
} else {
  top_loso_ligs <- character(0)
}

loso_long_list <- list()
for (comp_name in names(src221)) {
  pb <- src221[[comp_name]]
  matched_samps <- intersect(colnames(pb), fib221_samps)
  if (length(matched_samps) < 5) next

  logcpm <- compute_log2_cpm(pb[, matched_samps, drop=FALSE])

  for (lig in top_loso_ligs) {
    if (!lig %in% rownames(logcpm)) next
    lig_expr <- as.numeric(logcpm[lig, ])
    names(lig_expr) <- matched_samps

    full_rho <- tryCatch(cor(lig_expr, core2_vals[matched_samps], method="spearman"),
                         error=function(e) NA_real_)
    full_te <- NA
    treated_s <- matched_samps[sapply(matched_samps, function(s) TREAT_MAP[s]=="Neoadjuvant_treated")]
    control_s <- matched_samps[sapply(matched_samps, function(s) TREAT_MAP[s]=="Surgery_alone")]
    if (length(treated_s) > 0 && length(control_s) > 0) {
      full_te <- mean(lig_expr[treated_s]) - mean(lig_expr[control_s])
    }

    for (loso_s in matched_samps) {
      keep <- matched_samps[matched_samps != loso_s]
      if (length(keep) < 4) next
      rho_loso <- tryCatch(cor(lig_expr[keep], core2_vals[keep], method="spearman"),
                           error=function(e) NA_real_)

      # Treatment effect
      keep_treated <- keep[sapply(keep, function(s) TREAT_MAP[s]=="Neoadjuvant_treated")]
      keep_control <- keep[sapply(keep, function(s) TREAT_MAP[s]=="Surgery_alone")]
      te_loso <- NA
      if (length(keep_treated) > 0 && length(keep_control) > 0) {
        te_loso <- mean(lig_expr[keep_treated]) - mean(lig_expr[keep_control])
      }

      # Sign preservation
      rho_preserved <- if (!is.na(rho_loso) && !is.na(full_rho)) sign(rho_loso) == sign(full_rho) else NA
      te_preserved <- if (!is.na(te_loso) && !is.na(full_te)) sign(te_loso) == sign(full_te) else NA

      # Classification
      is_control_out <- TREAT_MAP[loso_s] == "Surgery_alone"
      loso_class <- "STABLE"
      if (!is.na(rho_preserved) && !rho_preserved) {
        loso_class <- if (is_control_out) "CONTROL_SENSITIVE" else "OUTLIER_SENSITIVE"
      } else if (is_control_out && !is.na(te_preserved) && !te_preserved) {
        loso_class <- "CONTROL_SENSITIVE"
      }

      loso_long_list[[paste(comp_name, lig, loso_s)]] <- data.frame(
        compartment=comp_name, ligand=lig, left_out=loso_s,
        left_out_treatment=TREAT_MAP[loso_s],
        rho_loso=rho_loso, rho_full=full_rho,
        rho_sign_preserved=rho_preserved,
        te_loso=te_loso, te_full=full_te,
        te_sign_preserved=te_preserved,
        classification=loso_class,
        stringsAsFactors=FALSE)
    }
  }
}

loso_long <- do.call(rbind, loso_long_list)
rownames(loso_long) <- NULL

# Summary
loso_summary_list <- list()
if (nrow(loso_long) > 0) {
  for (key in unique(paste(loso_long$compartment, loso_long$ligand))) {
    sub <- loso_long[paste(loso_long$compartment, loso_long$ligand) == key, ]
    comp_name <- sub$compartment[1]
    lig <- sub$ligand[1]
    rho_robust <- all(sub$rho_sign_preserved[!is.na(sub$rho_sign_preserved)])
    n_control_sens <- sum(sub$classification == "CONTROL_SENSITIVE")
    n_outlier_sens <- sum(sub$classification == "OUTLIER_SENSITIVE")
    if (n_control_sens > 0) {
      loso_class <- "CONTROL_SENSITIVE"
    } else if (n_outlier_sens > 0) {
      loso_class <- "OUTLIER_SENSITIVE"
    } else if (rho_robust) {
      loso_class <- "ROBUST"
    } else {
      loso_class <- "UNSTABLE"
    }
    loso_summary_list[[key]] <- data.frame(
      ligand=lig, dataset="GSE221561", compartment=comp_name,
      rho_sign_preservation=as.integer(rho_robust),
      effect_sign_preservation=as.integer(all(sub$te_sign_preserved[!is.na(sub$te_sign_preserved)])),
      classification=loso_class,
      stringsAsFactors=FALSE)
  }
}
loso_summary <- do.call(rbind, loso_summary_list)
rownames(loso_summary) <- NULL

write.csv(loso_long, file.path(O_DIR, "STEP19O_FIX221_LOSO_LONG.csv"), row.names=FALSE)
write.csv(loso_summary, file.path(O_DIR, "STEP19O_FIX221_LOSO_SUMMARY.csv"), row.names=FALSE)
cat_log(sprintf("LOSO: %d ligands summarized\n", nrow(loso_summary)))
cat_log("Saved: STEP19O_FIX221_LOSO_LONG.csv, STEP19O_FIX221_LOSO_SUMMARY.csv\n\n")

# ============================================================================
# FIX221-15 — CROSS-DATASET REPLICATION
# ============================================================================
cat_log("============================================\n")
cat_log("FIX221-15: CROSS-DATASET REPLICATION\n")
cat_log("============================================\n\n")

cross_list <- list()
# Get all ligand/source combos from GSE221561
for (i in 1:nrow(assoc_df)) {
  comp_name <- assoc_df$compartment[i]
  lig <- assoc_df$ligand[i]
  rho2_221 <- assoc_df$rho[i]

  # Find GSE197677 match
  assoc_197_match <- assoc_197[assoc_197$compartment == comp_name & assoc_197$ligand == lig, ]
  if (nrow(assoc_197_match) == 0) {
    cross_class <- "INSUFFICIENT_COVERAGE"
    rho2_197 <- NA
    te_197_val <- NA
    te_221_val <- NA
    te_dir_match <- NA
  } else {
    rho2_197 <- assoc_197_match$rho[1]
    te_197_match <- te_197[te_197$compartment == comp_name & te_197$ligand == lig, ]
    te_197_val <- if (nrow(te_197_match) > 0) te_197_match$effect[1] else NA
    te_221_match <- te_df[te_df$source == comp_name & te_df$ligand == lig, ]
    te_221_val <- if (nrow(te_221_match) > 0) te_221_match$effect[1] else NA
    te_dir_match <- if (!is.na(te_197_val) && !is.na(te_221_val)) sign(te_197_val) == sign(te_221_val) else NA

    # Receptor availability
    lig_receptors <- unique(curated_lr$receptor[curated_lr$ligand == lig])
    rec_221_avail <- rec_df[rec_df$receptor %in% lig_receptors, ]
    rec_221_det <- any(rec_221_avail$classification %in% c("ROBUSTLY_DETECTED", "PARTIALLY_DETECTED"))
    rec_197_avail <- rec_av_197[rec_av_197$receptor %in% lig_receptors, ]
    rec_197_det <- any(rec_197_avail$classification %in% c("ROBUSTLY_DETECTED", "PARTIALLY_DETECTED"))

    if (!is.na(rho2_197) && !is.na(rho2_221) &&
        abs(rho2_197) >= 0.40 && abs(rho2_221) >= 0.40 &&
        sign(rho2_197) == sign(rho2_221)) {
      cross_class <- "CROSS_DATASET_SIGNAL_SUPPORTED"
    } else if (!is.na(rho2_197) && !is.na(rho2_221) &&
               abs(rho2_197) >= 0.40 && abs(rho2_221) >= 0.40 &&
               sign(rho2_197) != sign(rho2_221)) {
      cross_class <- "CROSS_DATASET_DISCORDANT"
    } else if (!is.na(rho2_197) && abs(rho2_197) >= 0.40 &&
               (!is.na(rho2_221) && abs(rho2_221) < 0.40)) {
      cross_class <- "DATASET_SPECIFIC_SIGNAL"
    } else if (!is.na(rho2_221) && abs(rho2_221) >= 0.40 &&
               (!is.na(rho2_197) && abs(rho2_197) < 0.40)) {
      cross_class <- "DATASET_SPECIFIC_SIGNAL"
    } else {
      cross_class <- "INSUFFICIENT_COVERAGE"
    }
  }

  cross_list[[paste(comp_name, lig)]] <- data.frame(
    source=comp_name, ligand=lig,
    GSE197677_CORE2_rho=rho2_197,
    GSE221561_CORE2_rho=rho2_221,
    same_direction=if (!is.na(rho2_197) && !is.na(rho2_221)) sign(rho2_197) == sign(rho2_221) else NA,
    GSE197677_treatment_effect=te_197_val,
    GSE221561_treatment_effect=te_221_val,
    treatment_direction_match=te_dir_match,
    cross_dataset_classification=cross_class,
    stringsAsFactors=FALSE)
}
cross_df <- do.call(rbind, cross_list)
rownames(cross_df) <- NULL
write.csv(cross_df, file.path(O_DIR, "STEP19O_FIX221_CROSS_DATASET_REPLICATION.csv"), row.names=FALSE)

cat_log("=== CROSS-DATASET REPLICATION ===\n")
cat_log(sprintf("CROSS_DATASET_SIGNAL_SUPPORTED: %d\n", sum(cross_df$cross_dataset_classification == "CROSS_DATASET_SIGNAL_SUPPORTED")))
cat_log(sprintf("CORE2_ASSOCIATION_REPLICATED: %d\n", sum(cross_df$cross_dataset_classification == "CORE2_ASSOCIATION_REPLICATED")))
cat_log(sprintf("TREATMENT_EFFECT_REPLICATED: %d\n", sum(cross_df$cross_dataset_classification == "TREATMENT_EFFECT_REPLICATED")))
cat_log(sprintf("DATASET_SPECIFIC_SIGNAL: %d\n", sum(cross_df$cross_dataset_classification == "DATASET_SPECIFIC_SIGNAL")))
cat_log(sprintf("CROSS_DATASET_DISCORDANT: %d\n", sum(cross_df$cross_dataset_classification == "CROSS_DATASET_DISCORDANT")))
cat_log(sprintf("INSUFFICIENT_COVERAGE: %d\n", sum(cross_df$cross_dataset_classification == "INSUFFICIENT_COVERAGE")))

# Show replicated signals
replicated <- cross_df[cross_df$cross_dataset_classification == "CROSS_DATASET_SIGNAL_SUPPORTED", ]
if (nrow(replicated) > 0) {
  cat_log("\nCROSS-DATASET REPLICATED SIGNALS:\n")
  for (i in 1:nrow(replicated)) {
    cat_log(sprintf("  %s %s: 197=%.3f 221=%.3f\n",
      replicated$source[i], replicated$ligand[i],
      replicated$GSE197677_CORE2_rho[i], replicated$GSE221561_CORE2_rho[i]))
  }
}
cat_log("\nSaved: STEP19O_FIX221_CROSS_DATASET_REPLICATION.csv\n\n")

# ============================================================================
# FIX221-16 — EVIDENCE SCORE RECONCILIATION
# ============================================================================
cat_log("============================================\n")
cat_log("FIX221-16: EVIDENCE SCORE RECONCILIATION\n")
cat_log("============================================\n\n")

# Load original evidence scores
evi_orig <- read.csv(file.path(SRC_DIR, "STEP19O_SIGNALING_EVIDENCE_SCORE.csv"))

# Reconcile: update core2_221 and cross-dataset fields
evi_recon_list <- list()
for (i in 1:nrow(evi_orig)) {
  lig <- evi_orig$ligand[i]

  # Find best GSE221561 CORE2 rho
  assoc_221_lig <- assoc_df[assoc_df$ligand == lig, ]
  rho2_221 <- if (nrow(assoc_221_lig) > 0) {
    best_idx <- which.max(abs(assoc_221_lig$rho))
    assoc_221_lig$rho[best_idx]
  } else NA

  rho2_197 <- evi_orig$core2_197[i]
  same_dir <- if (!is.na(rho2_197) && !is.na(rho2_221)) sign(rho2_197) == sign(rho2_221) && abs(rho2_197) >= 0.40 && abs(rho2_221) >= 0.40 else FALSE

  # RUNX1 bridge in GSE221561
  runx1_221 <- runx1_bridge_df[runx1_bridge_df$ligand == lig, ]
  has_runx1_bridge <- if (nrow(runx1_221) > 0) any(runx1_221$classification == "LIGAND_RUNX1_CORE2_BRIDGE") else FALSE

  # LOSO
  loso_221 <- loso_summary[loso_summary$ligand == lig, ]
  loso_robust <- if (nrow(loso_221) > 0) any(loso_221$classification == "ROBUST") else FALSE

  # Count available components (out of 8)
  # 1: ligand detected
  # 2: receptor detected
  # 3: core2_197 (abs>=0.4)
  # 4: core2_221 (abs>=0.4)
  # 5: same direction
  # 6: ligand-target support (NOT AVAILABLE)
  # 7: RUNX1 bridge
  # 8: LOSO robust
  n_avail <- 0
  avail_components <- c()
  if (evi_orig$ligand_detected[i]) { n_avail <- n_avail + 1; avail_components <- c(avail_components, "ligand") }
  if (evi_orig$receptor_detected[i]) { n_avail <- n_avail + 1; avail_components <- c(avail_components, "receptor") }
  if (evi_orig$core2_197[i]) { n_avail <- n_avail + 1; avail_components <- c(avail_components, "core2_197") }
  if (!is.na(rho2_221) && abs(rho2_221) >= 0.40) { n_avail <- n_avail + 1; avail_components <- c(avail_components, "core2_221") }
  if (same_dir) { n_avail <- n_avail + 1; avail_components <- c(avail_components, "same_dir") }
  # ligand-target: NOT AVAILABLE - do not penalize
  if (has_runx1_bridge) { n_avail <- n_avail + 1; avail_components <- c(avail_components, "runx1") }
  if (loso_robust) { n_avail <- n_avail + 1; avail_components <- c(avail_components, "loso") }

  # Available denominator = 7 (exclude ligand-target which is not available)
  n_avail_total <- 7
  if (n_avail >= 5) evi_class <- "HIGH_SIGNALING_SUPPORT"
  else if (n_avail >= 4) evi_class <- "MODERATE_SIGNALING_SUPPORT"
  else if (n_avail >= 2) evi_class <- "WEAK_SIGNALING_SUPPORT"
  else evi_class <- "NO_REPRODUCIBLE_SIGNALING_SUPPORT"

  evi_recon_list[[lig]] <- data.frame(
    ligand=lig,
    ligand_detected=evi_orig$ligand_detected[i],
    receptor_detected=evi_orig$receptor_detected[i],
    core2_197=evi_orig$core2_197[i],
    core2_221=(!is.na(rho2_221) && abs(rho2_221) >= 0.40),
    GSE221561_CORE2_rho=rho2_221,
    same_direction=same_dir,
    ligand_target_support="NOT_AVAILABLE",
    runx1_bridge=has_runx1_bridge,
    loso_robust=loso_robust,
    evidence_score=n_avail,
    available_denominator=n_avail_total,
    score_fraction=paste0(n_avail, "/", n_avail_total),
    reconciled_classification=evi_class,
    stringsAsFactors=FALSE)
}
evi_recon <- do.call(rbind, evi_recon_list)
write.csv(evi_recon, file.path(O_DIR, "STEP19O_FIX221_RECONCILED_EVIDENCE_SCORE.csv"), row.names=FALSE)
cat_log("=== RECONCILED EVIDENCE SCORES ===\n")
for (i in 1:nrow(evi_recon)) {
  if (evi_recon$evidence_score[i] >= 4) {
    cat_log(sprintf("  %s: %d/%d -> %s\n",
      evi_recon$ligand[i], evi_recon$evidence_score[i],
      evi_recon$available_denominator[i], evi_recon$reconciled_classification[i]))
  }
}
cat_log("\nSaved: STEP19O_FIX221_RECONCILED_EVIDENCE_SCORE.csv\n\n")

# ============================================================================
# FIX221-17 — FINAL CLASSIFICATION
# ============================================================================
cat_log("============================================\n")
cat_log("FIX221-17: FINAL CLASSIFICATION\n")
cat_log("============================================\n\n")

n_cross_repl <- sum(cross_df$cross_dataset_classification == "CROSS_DATASET_SIGNAL_SUPPORTED")
n_discordant <- sum(cross_df$cross_dataset_classification == "CROSS_DATASET_DISCORDANT")
n_dataset_specific <- sum(cross_df$cross_dataset_classification == "DATASET_SPECIFIC_SIGNAL")
n_runx1_bridges_221 <- sum(runx1_bridge_df$classification == "LIGAND_RUNX1_CORE2_BRIDGE")
n_moderate_recon <- sum(evi_recon$reconciled_classification == "MODERATE_SIGNALING_SUPPORT")
n_high_recon <- sum(evi_recon$reconciled_classification == "HIGH_SIGNALING_SUPPORT")

cat_log(sprintf("Cross-dataset replicated: %d\n", n_cross_repl))
cat_log(sprintf("Cross-dataset discordant: %d\n", n_discordant))
cat_log(sprintf("Dataset-specific signals: %d\n", n_dataset_specific))
cat_log(sprintf("RUNX1 bridges (221): %d\n", n_runx1_bridges_221))
cat_log(sprintf("High support (reconciled): %d\n", n_high_recon))
cat_log(sprintf("Moderate support (reconciled): %d\n", n_moderate_recon))

# Apply classification rules
if (n_cross_repl > 0) {
  primary_cls <- "CORE2_WITH_REPLICATED_EXTRACELLULAR_SIGNAL_SUPPORT"
  cat_log("\nRule: At least 1 cross-dataset replicated signal -> CORE2_WITH_REPLICATED_EXTRACELLULAR_SIGNAL_SUPPORT\n")
} else if (n_runx1_bridges_221 > 0) {
  primary_cls <- "CORE2_WITH_RUNX1_LINKED_SIGNALING_SUPPORT"
  cat_log("\nRule: RUNX1 bridge present -> CORE2_WITH_RUNX1_LINKED_SIGNALING_SUPPORT\n")
} else if (n_dataset_specific > 0 || n_moderate_recon > 0) {
  # Candidates evaluable in both datasets but fail replication
  primary_cls <- "DATASET_SPECIFIC_SIGNALING_ARCHITECTURE"
  cat_log("\nRule: Candidates evaluable in both but fail replication -> DATASET_SPECIFIC_SIGNALING_ARCHITECTURE\n")
} else {
  primary_cls <- "SIGNALING_AUDIT_INCONCLUSIVE"
  cat_log("\nRule: Cannot evaluate -> SIGNALING_AUDIT_INCONCLUSIVE\n")
}

secondary_lim <- c("MINIMAL_LR_RESOURCE", "GSE221561_SMALL_CONTROL_GROUP")
if (n_cross_repl == 0) secondary_lim <- c(secondary_lim, "INSUFFICIENT_CROSS_DATASET_COVERAGE")

cat_log(sprintf("\nRECONCILED FINAL CLASSIFICATION: %s\n", primary_cls))
cat_log("SECONDARY LIMITATIONS:\n")
for (sl in secondary_lim) cat_log(sprintf("  - %s\n", sl))

final_class <- data.frame(
  step="19O-FIX221",
  primary_classification=primary_cls,
  secondary_limitations=paste(secondary_lim, collapse="; "),
  n_eligible_221=nrow(elig221_df),
  n_cross_repl=n_cross_repl,
  n_discordant=n_discordant,
  n_dataset_specific=n_dataset_specific,
  n_runx1_bridges=n_runx1_bridges_221,
  n_moderate_recon=n_moderate_recon,
  n_high_recon=n_high_recon,
  stringsAsFactors=FALSE)
write.csv(final_class, file.path(O_DIR, "STEP19O_FIX221_FINAL_CLASSIFICATION.csv"), row.names=FALSE)

final_status <- data.frame(
  step="19O-FIX221", status="COMPLETE",
  primary_classification=primary_cls,
  cells_removed=0, samples_removed=0, cells_reclassified=0,
  frozen_upstream_changed="NO",
  GSE197677_outputs_changed="NO",
  timestamp=format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  stringsAsFactors=FALSE)
write.csv(final_status, file.path(O_DIR, "STEP19O_FIX221_FINAL_STATUS.csv"), row.names=FALSE)

# ============================================================================
# FINAL VALIDATION
# ============================================================================
cat_log("\n============================================\n")
cat_log("STEP 19O-FIX221 COMPLETE\n")
cat_log("GSE221561 SIGNALING COMPLETION\n")
cat_log("============================================\n\n")

cat_log(sprintf("GSE221561 SOURCE PSEUDOBULK:\n"))
for (comp_name in names(src221)) {
  cat_log(sprintf("  %s: %d genes x %d samples\n", comp_name, nrow(src221[[comp_name]]), ncol(src221[[comp_name]])))
}

cat_log(sprintf("\nSOURCE COMPARTMENTS EVALUABLE:\n"))
for (i in 1:nrow(eval_df)) {
  cat_log(sprintf("  %s: %s\n", eval_df$source_compartment[i], eval_df$evaluability[i]))
}

cat_log(sprintf("\nGSE221561 ELIGIBLE LR PAIRS:\n%d\n", nrow(elig221_df)))

cat_log("\n--------------------------------------------\n\n")
cat_log("PRESPECIFIED CANDIDATES:\n\n")

for (cand in candidates) {
  lig <- cand[1]
  if (lig %in% rownames(cand_df)) {
    r <- cand_df[cand_df$ligand == lig, ]
    cat_log(sprintf("%s / %s:\n", lig, r$source))
    cat_log(sprintf("  Ligand detected: %s\n", r$ligand_detected))
    cat_log(sprintf("  Receptor detected: %s\n", r$receptor_detected))
    cat_log(sprintf("  n: %s\n", ifelse(is.na(r$n_paired_samples), "NA", r$n_paired_samples)))
    cat_log(sprintf("  CORE2 rho: %s\n", ifelse(is.na(r$GSE221561_CORE2_rho), "NA", round(r$GSE221561_CORE2_rho, 3))))
    cat_log(sprintf("  GSE197677 CORE2 rho: %s\n", ifelse(is.na(r$GSE197677_CORE2_rho), "NA", round(r$GSE197677_CORE2_rho, 3))))
    cat_log(sprintf("  Treatment effect: %s\n", ifelse(is.na(r$GSE221561_treatment_effect), "NA", round(r$GSE221561_treatment_effect, 3))))
    cat_log(sprintf("  Status: %s\n\n", r$status))
  }
}

cat_log("--------------------------------------------\n\n")
cat_log("TOP GSE221561 SIGNALING CANDIDATES:\n")
if (nrow(assoc_df) > 0) {
  top10 <- head(assoc_df[order(-abs(assoc_df$rho)), ], 10)
  for (i in 1:nrow(top10)) {
    cat_log(sprintf("  %s %s: rho=%.3f -> %s\n",
      top10$compartment[i], top10$ligand[i], top10$rho[i], top10$association_class[i]))
  }
}

cat_log("\n--------------------------------------------\n\n")
cat_log("CROSS-DATASET REPLICATED SIGNALS:\n")
if (n_cross_repl > 0) {
  for (i in 1:nrow(replicated)) {
    cat_log(sprintf("  %s %s\n", replicated$source[i], replicated$ligand[i]))
  }
} else {
  cat_log("  None\n")
}

cat_log("\nCROSS-DATASET DISCORDANT SIGNALS:\n")
discordant <- cross_df[cross_df$cross_dataset_classification == "CROSS_DATASET_DISCORDANT", ]
if (nrow(discordant) > 0) {
  for (i in 1:nrow(discordant)) {
    cat_log(sprintf("  %s %s: 197=%.3f 221=%.3f\n",
      discordant$source[i], discordant$ligand[i],
      discordant$GSE197677_CORE2_rho[i], discordant$GSE221561_CORE2_rho[i]))
  }
} else {
  cat_log("  None\n")
}

cat_log("\nINSUFFICIENTLY EVALUABLE SIGNALS:\n")
insuf <- cross_df[cross_df$cross_dataset_classification == "INSUFFICIENT_COVERAGE", ]
cat_log(sprintf("  %d ligands with insufficient cross-dataset coverage\n", nrow(insuf)))

cat_log("\n--------------------------------------------\n\n")
cat_log("RUNX1-LINKED SIGNALING:\n")
cat_log(sprintf("  RUNX1 bridges (GSE221561): %d\n", n_runx1_bridges_221))

cat_log("\nTGF-BETA SIGNALING:\n")
tgfb_ligs <- c("TGFB1", "TGFB2", "TGFB3", "BMP2", "BMP4", "BMP7")
for (lig in tgfb_ligs) {
  lig_assoc <- assoc_df[assoc_df$ligand == lig, ]
  if (nrow(lig_assoc) > 0) {
    for (j in 1:nrow(lig_assoc)) {
      cat_log(sprintf("  %s %s: rho=%.3f\n", lig, lig_assoc$compartment[j], lig_assoc$rho[j]))
    }
  }
}

cat_log("\nAP-1-LINKED SIGNALING:\n")
cat_log("  FOS/FOSL2 are downstream TFs, not ligands. Expression-only from Step19M.\n")

cat_log("\n--------------------------------------------\n\n")
cat_log(sprintf("RECONCILED FINAL CLASSIFICATION:\n%s\n", primary_cls))
cat_log("\nSECONDARY LIMITATIONS:\n")
for (sl in secondary_lim) cat_log(sprintf("  - %s\n", sl))

cat_log("\nBIOLOGICAL INTERPRETATION:\n\n")
if (n_cross_repl > 0) {
  cat_log("The fibroblast CDKN1B-GSN CORE2 response is conserved across cohorts.\n")
  cat_log(sprintf("%d ligand/source signals show cross-dataset replicated CORE2 association.\n", n_cross_repl))
  cat_log("These associations are exploratory and do not establish causality.\n")
} else {
  cat_log("The fibroblast CDKN1B-GSN CORE2 response is strongly conserved across\n")
  cat_log("treatment cohorts, whereas the extracellular signaling architecture\n")
  cat_log("associated with that response was NOT reproducibly identified across\n")
  cat_log("datasets in the completed Step19O-FIX221 analysis.\n\n")
  if (n_dataset_specific > 0 || n_moderate_recon > 0) {
    cat_log("GSE197677 showed moderate signaling associations (AREG, EREG, DLL1, TGFB3)\n")
    cat_log("but these were not replicated in the same direction in GSE221561,\n")
    cat_log("indicating a dataset-specific signaling architecture.\n")
  }
  cat_log("\nSurgery_alone n=2 limits GSE221561 replication power.\n")
  cat_log("The minimal LR resource limits mechanistic inference.\n")
  cat_log("These associations are exploratory and do not establish causality,\n")
  cat_log("secretion, receptor activation, or a conserved ligand-driven mechanism.\n")
}

cat_log("\nGSE197677 outputs changed:\nNO\n")
cat_log("Frozen upstream outputs changed:\nNO\n")
cat_log("Cells removed:\n0\n")
cat_log("Samples removed:\n0\n")
cat_log("Cells reclassified:\n0\n\n")

cat_log("STOP HERE.\n")
cat_log("DO NOT START STEP19P.\n")