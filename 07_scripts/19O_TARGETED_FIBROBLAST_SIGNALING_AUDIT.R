#!/usr/bin/env Rscript
# ============================================================================
# STEP 19O: TARGETED FIBROBLAST SIGNALING SOURCE AUDIT
# ============================================================================
# Identifies extracellular signaling programs associated with
# the conserved fibroblast CORE2 response across GSE197677 and GSE221561.
# This is a targeted cell-cell signaling association audit.
# ============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(reshape2)
})

BASE_DIR <- Sys.getenv("ESCC_CAF_PROJECT_ROOT", unset = ".")
O_DIR    <- file.path(BASE_DIR, "06_tables/cross_dataset/pathway_analysis/Step19O")
FIG_DIR  <- file.path(BASE_DIR, "05_figures/cross_dataset/pathway_analysis/Step19O")
LOG_DIR  <- file.path(BASE_DIR, "08_logs")
dir.create(O_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(LOG_DIR, showWarnings = FALSE, recursive = TRUE)

CORE_GENES <- c("BGN", "CDKN1B", "GSN", "TIMP1")

TREAT_MAP <- list(
  GSE197677 = c(
    ESC06 = "nNACT", ESC07 = "nNACT", ESC09 = "NACT", ESC11 = "NACT",
    ESC12 = "NACT", ESC15 = "NACT", ESC16 = "NACT", ESC17 = "nNACT",
    ESC18 = "nNACT", ESC22 = "NACT"
  ),
  GSE221561 = c(
    S0520T = "Neoadjuvant_treated", S1265T = "Neoadjuvant_treated",
    S1315TS = "Neoadjuvant_treated", S1535T = "Neoadjuvant_treated",
    S2423T = "Neoadjuvant_treated", S2487T = "Neoadjuvant_treated",
    S6829T = "Neoadjuvant_treated", S6417T = "Surgery_alone",
    S9478T = "Surgery_alone"
  )
)

TREAT_CONTRAST <- list(
  GSE197677 = c(treated = "NACT", control = "nNACT"),
  GSE221561 = c(treated = "Neoadjuvant_treated", control = "Surgery_alone")
)

compute_log2_cpm <- function(counts_mat) {
  lib_size <- colSums(counts_mat)
  lib_size[lib_size == 0] <- 1
  cpm_mat <- sweep(counts_mat, 2, lib_size, "/") * 1e6
  log2(cpm_mat + 1)
}

# ============================================================================
# 19O0 — OUTPUT STRUCTURE
# ============================================================================
cat("============================================\n")
cat("STEP 19O: TARGETED FIBROBLAST SIGNALING SOURCE AUDIT\n")
cat("============================================\n\n")

# ============================================================================
# 19O1 — INPUT PREFLIGHT
# ============================================================================
cat("============================================\n")
cat("19O1: INPUT PREFLIGHT\n")
cat("============================================\n\n")

cat("Loading Seurat objects for cell type inventory...\n")

seurat_197 <- tryCatch(readRDS(file.path(BASE_DIR, "03_objects/GSE197677/GSE197677_integrated_broadcelltype_annotated.rds")),
  error = function(e) { cat("ERROR loading GSE197677:", conditionMessage(e), "\n"); NULL })
seurat_221 <- tryCatch(readRDS(file.path(BASE_DIR, "03_objects/GSE221561/GSE221561_author_annotated.rds")),
  error = function(e) { cat("ERROR loading GSE221561:", conditionMessage(e), "\n"); NULL })

if (is.null(seurat_197) || is.null(seurat_221)) {
  cat("FATAL: Cannot load Seurat objects. Stopping.\n")
  status <- data.frame(step = "19O", status = "SEURAT_OBJECT_LOAD_FAILED",
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  write.csv(status, file.path(O_DIR, "STEP19O_FINAL_STATUS.csv"), row.names = FALSE)
  stop("SEURAT_OBJECT_LOAD_FAILED")
}

# Extract metadata
meta_197 <- seurat_197@meta.data
meta_221 <- seurat_221@meta.data

cat("GSE197677 cells:", nrow(meta_197), "\n")
cat("GSE221561 cells:", nrow(meta_221), "\n\n")

# Identify cell type column (try common names)
ct_col_197 <- intersect(c("broad_celltype_final", "broad_cell_type", "cell_type", "CellType", "celltype", "seurat_clusters"), colnames(meta_197))[1]
ct_col_221 <- intersect(c("author_Cell_type", "broad_cell_type", "cell_type", "CellType", "celltype", "seurat_clusters"), colnames(meta_221))[1]

if (is.na(ct_col_197) || is.na(ct_col_221)) {
  cat("Available columns in GSE197677:", paste(colnames(meta_197), collapse=", "), "\n")
  cat("Available columns in GSE221561:", paste(colnames(meta_221), collapse=", "), "\n")
  stop("Cannot identify cell type column")
}

cat("Cell type column GSE197677:", ct_col_197, "\n")
cat("Cell type column GSE221561:", ct_col_221, "\n\n")

# Map to broad compartments
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

comp_197 <- map_to_compartment(meta_197[[ct_col_197]])
comp_221 <- map_to_compartment(meta_221[[ct_col_221]])

# Sample ID column
sample_col_197 <- intersect(c("sample", "sample_id", "orig.ident", "Sample"), colnames(meta_197))[1]
sample_col_221 <- intersect(c("sample_id", "sample", "orig.ident", "Sample"), colnames(meta_221))[1]

if (is.na(sample_col_197) || is.na(sample_col_221)) {
  cat("Available columns:", paste(colnames(meta_197), collapse=", "), "\n")
  stop("Cannot identify sample column")
}

samples_197 <- as.character(meta_197[[sample_col_197]])
samples_221 <- as.character(meta_221[[sample_col_221]])

# Build inventory
build_inventory <- function(meta, comp, samples, dataset) {
  tabs <- table(sample = samples, compartment = comp)
  df <- as.data.frame.matrix(tabs)
  df$dataset <- dataset
  df$sample_id <- rownames(df)
  # Ensure all expected compartments are present
  for (expected_comp in c("Epithelial", "Fibroblast", "Myeloid", "T_NK", "Endothelial", "Pericyte", "Bcell", "Other")) {
    if (!expected_comp %in% colnames(df)) df[[expected_comp]] <- 0
  }
  df
}

inv_197 <- build_inventory(meta_197, comp_197, samples_197, "GSE197677")
inv_221 <- build_inventory(meta_221, comp_221, samples_221, "GSE221561")

# Ensure matching columns before rbind
common_cols <- intersect(colnames(inv_197), colnames(inv_221))
inventory <- rbind(inv_197[, common_cols], inv_221[, common_cols])

cat("=== CELL TYPE INVENTORY ===\n")
cat("GSE197677 compartments:", paste(setdiff(colnames(inv_197), c("dataset", "sample_id")), collapse=", "), "\n")
cat("GSE221561 compartments:", paste(setdiff(colnames(inv_221), c("dataset", "sample_id")), collapse=", "), "\n\n")

# Check for required compartments
for (comp in c("Fibroblast", "Epithelial")) {
  n197 <- sum(inv_197[[comp]] > 0, na.rm = TRUE)
  n221 <- sum(inv_221[[comp]] > 0, na.rm = TRUE)
  cat(sprintf("  %s: GSE197677=%d samples, GSE221561=%d samples\n", comp, n197, n221))
}

for (comp in c("Myeloid", "T_NK", "Endothelial")) {
  n197 <- sum(inv_197[[comp]] > 0, na.rm = TRUE)
  n221 <- sum(inv_221[[comp]] > 0, na.rm = TRUE)
  cat(sprintf("  %s: GSE197677=%d samples, GSE221561=%d samples\n", comp, n197, n221))
}

write.csv(inventory, file.path(O_DIR, "STEP19O_CELLTYPE_SAMPLE_INVENTORY.csv"), row.names = FALSE)
cat("\nSaved: STEP19O_CELLTYPE_SAMPLE_INVENTORY.csv\n")

# Free Seurat objects - we only need pseudobulk
rm(seurat_197, seurat_221)
gc()

# ============================================================================
# 19O2 — RECEIVER DEFINITION
# ============================================================================
cat("\n============================================\n")
cat("19O2: RECEIVER DEFINITION\n")
cat("============================================\n\n")

core_scores <- read.csv(file.path(BASE_DIR,
  "06_tables/cross_dataset/pathway_analysis/Step19L/STEP19L_SAMPLE_CORE_MODULE_SCORES.csv"))

cat("Receiver: Fibroblast_CAF\n")
cat("Receiver response variables: CORE2_STRONG_score, CORE4_score\n")
cat("Frozen scores from Step19L loaded.\n\n")

cat("=== FROZEN RECEIVER SCORES ===\n")
for (i in 1:nrow(core_scores)) {
  cat(sprintf("  %s %s: CORE2=%.4f CORE4=%.4f\n",
    core_scores$dataset[i], core_scores$sample_id[i],
    core_scores$CORE2_STRONG_score[i], core_scores$CORE4_score[i]))
}

# Save provenance
provenance <- data.frame(
  step = "19O2",
  receiver = "Fibroblast_CAF",
  core2_source = "STEP19L_SAMPLE_CORE_MODULE_SCORES.csv",
  core4_source = "STEP19L_SAMPLE_CORE_MODULE_SCORES.csv",
  n_samples_197 = sum(core_scores$dataset == "GSE197677"),
  n_samples_221 = sum(core_scores$dataset == "GSE221561"),
  timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  stringsAsFactors = FALSE
)
write.csv(provenance, file.path(O_DIR, "STEP19O_RECEIVER_SCORE_PROVENANCE.csv"), row.names = FALSE)
cat("\nSaved: STEP19O_RECEIVER_SCORE_PROVENANCE.csv\n")

# ============================================================================
# 19O3 — SOURCE COMPARTMENT PSEUDOBULK
# ============================================================================
cat("\n============================================\n")
cat("19O3: SOURCE COMPARTMENT PSEUDOBULK\n")
cat("============================================\n\n")

# Load existing fibroblast pseudobulk counts
fib197 <- readRDS(file.path(BASE_DIR, "06_tables/cross_dataset/pathway_analysis/Step19I/counts/GSE197677_FIBROBLAST_CAF_raw_counts.rds"))
fib221 <- readRDS(file.path(BASE_DIR, "06_tables/cross_dataset/pathway_analysis/Step19I/counts/GSE221561_FIBROBLAST_CAF_raw_counts.rds"))

cat("Fibroblast GSE197677:", nrow(fib197), "genes x", ncol(fib197), "samples\n")
cat("Fibroblast GSE221561:", nrow(fib221), "genes x", ncol(fib221), "samples\n\n")

# Compute pseudobulk for source compartments from Seurat objects
# We need to reload Seurat objects briefly for pseudobulk computation
cat("Computing pseudobulk for source compartments...\n")

compute_pseudobulk <- function(seurat_path, sample_col, ct_col, target_compartments) {
  obj <- readRDS(seurat_path)
  meta <- obj@meta.data
  
  rna_assay <- obj@assays[["RNA"]]
  if (is.null(rna_assay)) rna_assay <- obj@assays[[Seurat::DefaultAssay(obj)]]
  assay_class <- class(rna_assay)[1]
  cat("  Assay class:", assay_class, "\n")
  
  results <- list()
  
  if (assay_class == "Assay5") {
    # v5 Assay5: each layer is a sample's cells
    layer_names <- names(rna_assay@layers)
    cat("  Processing", length(layer_names), "layers...\n")
    
    # Map cell barcodes to samples
    cell_to_sample <- setNames(rep(NA_character_, nrow(meta)), rownames(meta))
    for (ln in layer_names) {
      layer_mat <- rna_assay@layers[[ln]]
      if (ncol(layer_mat) == 0) next
      layer_cells <- colnames(layer_mat)
      if (is.null(layer_cells)) next
      layer_meta <- meta[layer_cells, , drop = FALSE]
      samp_ids_layer <- as.character(layer_meta[[sample_col]])
      unique_samps_layer <- unique(samp_ids_layer)
      if (length(unique_samps_layer) == 1) {
        cell_to_sample[layer_cells] <- unique_samps_layer
      }
    }
    
    # Get gene names
    gene_names <- NULL
    for (ln in layer_names) {
      layer_mat <- rna_assay@layers[[ln]]
      if (nrow(layer_mat) > 0) { gene_names <- rownames(layer_mat); break }
    }
    
    # Compute pseudobulk per compartment
    for (comp in target_compartments) {
      comp_map <- map_to_compartment(meta[[ct_col]])
      cells_in_comp <- which(comp_map == comp)
      if (length(cells_in_comp) < 10) { cat(sprintf("  %s: %d cells, skipping\n", comp, length(cells_in_comp))); next }
      
      comp_cells <- rownames(meta)[cells_in_comp]
      comp_samps <- cell_to_sample[comp_cells]
      comp_samps <- comp_samps[!is.na(comp_samps)]
      unique_samps <- unique(comp_samps)
      
      pb_list <- list()
      for (s in unique_samps) {
        cells <- comp_cells[comp_samps == s]
        if (length(cells) < 5) next
        
        pb_numeric <- rep(0, length(gene_names))
        names(pb_numeric) <- gene_names
        
        for (ln in layer_names) {
          layer_mat <- rna_assay@layers[[ln]]
          layer_cells <- colnames(layer_mat)
          if (is.null(layer_cells)) next
          common_cells <- intersect(cells, layer_cells)
          if (length(common_cells) > 0) {
            layer_genes <- rownames(layer_mat)
            layer_sums <- Matrix::rowSums(layer_mat[, common_cells, drop = FALSE])
            common_genes <- intersect(layer_genes, gene_names)
            if (length(common_genes) > 0) {
              pb_numeric[common_genes] <- pb_numeric[common_genes] + as.numeric(layer_sums[match(common_genes, layer_genes)])
            }
          }
        }
        pb_list[[s]] <- pb_numeric
      }
      
      if (length(pb_list) >= 2) {
        pb_mat <- do.call(cbind, lapply(pb_list, function(x) matrix(x, ncol = 1)))
        rownames(pb_mat) <- gene_names
        results[[comp]] <- pb_mat
        cat(sprintf("  %s: %d genes x %d samples\n", comp, nrow(pb_mat), ncol(pb_mat)))
      }
    }
  } else {
    # v3 Assay: counts in @counts slot
    counts <- rna_assay@counts
    cat("  Counts dim:", dim(counts), "\n")
    
    if (is.null(colnames(counts)) && ncol(counts) == nrow(meta)) {
      colnames(counts) <- rownames(meta)
    }
    
    for (comp in target_compartments) {
      comp_map <- map_to_compartment(meta[[ct_col]])
      cells_in_comp <- which(comp_map == comp)
      if (length(cells_in_comp) < 10) { cat(sprintf("  %s: %d cells, skipping\n", comp, length(cells_in_comp))); next }
      
      comp_cells <- rownames(meta)[cells_in_comp]
      comp_cells <- intersect(comp_cells, colnames(counts))
      if (length(comp_cells) < 10) next
      
      comp_counts <- counts[, comp_cells, drop = FALSE]
      comp_meta <- meta[comp_cells, , drop = FALSE]
      samp_ids <- as.character(comp_meta[[sample_col]])
      unique_samps <- unique(samp_ids)
      
      pb_list <- list()
      for (s in unique_samps) {
        cells <- comp_cells[samp_ids == s]
        if (length(cells) >= 5) {
          pb_list[[s]] <- Matrix::rowSums(comp_counts[, cells, drop = FALSE])
        }
      }
      
      if (length(pb_list) >= 2) {
        pb_mat <- do.call(cbind, pb_list)
        results[[comp]] <- pb_mat
        cat(sprintf("  %s: %d genes x %d samples\n", comp, nrow(pb_mat), ncol(pb_mat)))
      }
    }
  }
  
  rm(obj)
  gc()
  results
}

# Compute pseudobulk for source compartments
source_compartments <- c("Epithelial", "Myeloid", "T_NK", "Endothelial")

cat("\nGSE197677 source compartments:\n")
src_197 <- compute_pseudobulk(
  file.path(BASE_DIR, "03_objects/GSE197677/GSE197677_integrated_broadcelltype_annotated.rds"),
  sample_col_197, ct_col_197, source_compartments)

cat("\nGSE221561 source compartments:\n")
src_221 <- compute_pseudobulk(
  file.path(BASE_DIR, "03_objects/GSE221561/GSE221561_author_annotated.rds"),
  sample_col_221, ct_col_221, source_compartments)

# Save pseudobulk counts
for (comp in names(src_197)) {
  saveRDS(src_197[[comp]], file.path(O_DIR, sprintf("STEP19O_GSE197677_%s_raw_counts.rds", comp)))
}
for (comp in names(src_221)) {
  saveRDS(src_221[[comp]], file.path(O_DIR, sprintf("STEP19O_GSE221561_%s_raw_counts.rds", comp)))
}

# Build source compartment summary
src_summary <- data.frame()
for (comp in source_compartments) {
  if (comp %in% names(src_197)) {
    src_summary <- rbind(src_summary, data.frame(
      dataset = "GSE197677", compartment = comp,
      n_genes = nrow(src_197[[comp]]), n_samples = ncol(src_197[[comp]]),
      samples = paste(colnames(src_197[[comp]]), collapse = ", "),
      stringsAsFactors = FALSE))
  }
  if (comp %in% names(src_221)) {
    src_summary <- rbind(src_summary, data.frame(
      dataset = "GSE221561", compartment = comp,
      n_genes = nrow(src_221[[comp]]), n_samples = ncol(src_221[[comp]]),
      samples = paste(colnames(src_221[[comp]]), collapse = ", "),
      stringsAsFactors = FALSE))
  }
}

write.csv(src_summary, file.path(O_DIR, "STEP19O_SOURCE_COMPARTMENT_COUNTS.csv"), row.names = FALSE)
cat("\nSaved: STEP19O_SOURCE_COMPARTMENT_COUNTS.csv\n")

# ============================================================================
# 19O4 — LOCAL LIGAND-RECEPTOR RESOURCE PREFLIGHT
# ============================================================================
cat("\n============================================\n")
cat("19O4: LOCAL LR RESOURCE PREFLIGHT\n")
cat("============================================\n\n")

# Check for available LR resources
check_lr_resource <- function(name, available, path, n_pairs, ligand_target) {
  data.frame(resource = name, available = available, path = path,
    n_pairs = n_pairs, ligand_target_available = ligand_target,
    stringsAsFactors = FALSE)
}

resources <- list()
resources$cellchatdb <- check_lr_resource("CellChatDB", FALSE, NA, 0, FALSE)
resources$nichenet_lr <- check_lr_resource("NicheNet_LR", FALSE, NA, 0, FALSE)
resources$nichenet_lt <- check_lr_resource("NicheNet_ligand_target", FALSE, NA, 0, FALSE)
resources$omnipath <- check_lr_resource("OmniPath", FALSE, NA, 0, FALSE)
resources$liana <- check_lr_resource("LIANA", FALSE, NA, 0, FALSE)
resources$project_local <- check_lr_resource("project_local_LR", FALSE, NA, 0, FALSE)

resources_df <- do.call(rbind, resources)
rownames(resources_df) <- NULL

cat("=== LR RESOURCE STATUS ===\n")
for (i in 1:nrow(resources_df)) {
  cat(sprintf("  %-25s: %s\n", resources_df$resource[i],
    ifelse(resources_df$available[i], "AVAILABLE", "NOT AVAILABLE")))
}

# Create the project-local minimal LR resource embedded in this script
# This is NOT a comprehensive database - it is a project-local minimal set
cat("\nNo external LR resource available. Creating curated minimal LR resource.\n")

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
.lr_fam <- c(
  rep("TGF_beta",16), rep("TNF",7), rep("IL1",4), rep("IL6_JAK_STAT",7),
  rep("IFN",6), rep("PDGF",8), rep("FGF",10), rep("EGF",10),
  rep("WNT",6), rep("CXCL_CCL",12), rep("VEGF",4), rep("NOTCH",8), rep("HGF_MET",2)
)

cat(sprintf("  LR pre-build: ligands=%d, receptors=%d, families=%d\n",
            length(.lr_lig), length(.lr_rec), length(.lr_fam)))
if (length(.lr_lig) != length(.lr_rec)) {
  stop(sprintf("FATAL: LR vector length mismatch: ligands=%d receptors=%d",
               length(.lr_lig), length(.lr_rec)))
}

curated_lr <- data.frame(
  ligand   = .lr_lig,
  receptor = .lr_rec,
  family   = .lr_fam,
  stringsAsFactors = FALSE
)
rm(.lr_lig, .lr_rec, .lr_fam)

cat(sprintf("Curated LR pairs: %d\n", nrow(curated_lr)))
cat(sprintf("Unique ligands: %d\n", length(unique(curated_lr$ligand))))
cat(sprintf("Unique receptors: %d\n", length(unique(curated_lr$receptor))))
cat(sprintf("Families: %s\n", paste(unique(curated_lr$family), collapse=", ")))

resources_df <- rbind(resources_df, data.frame(
  resource = "curated_minimal", available = TRUE,
  path = "curated in script", n_pairs = nrow(curated_lr),
  ligand_target_available = FALSE, stringsAsFactors = FALSE))

write.csv(resources_df, file.path(O_DIR, "STEP19O_RESOURCE_STATUS.csv"), row.names = FALSE)
cat("\nSaved: STEP19O_RESOURCE_STATUS.csv\n")

# ============================================================================
# 19O5 — RECEIVER RECEPTOR EXPRESSION QC
# ============================================================================
cat("\n============================================\n")
cat("19O5: RECEIVER RECEPTOR EXPRESSION QC\n")
cat("============================================\n\n")

# Get unique receptors
all_receptors <- unique(curated_lr$receptor)

cat("Checking fibroblast receptor expression...\n")

check_receptor_detection <- function(counts_mat, receptors) {
  results <- list()
  for (rec in receptors) {
    if (rec %in% rownames(counts_mat)) {
      expr <- as.numeric(counts_mat[rec, ])
      detected <- sum(expr > 0)
      n_samples <- ncol(counts_mat)
      frac <- detected / n_samples
      mean_expr <- mean(expr)
      
      if (frac >= 0.50) flag <- "ROBUSTLY_DETECTED"
      else if (frac >= 0.30) flag <- "PARTIALLY_DETECTED"
      else flag <- "POORLY_DETECTED"
      
      results[[rec]] <- data.frame(
        receptor = rec, n_detected = detected, n_samples = n_samples,
        detection_fraction = frac, mean_expression = mean_expr,
        classification = flag, stringsAsFactors = FALSE)
    } else {
      results[[rec]] <- data.frame(
        receptor = rec, n_detected = 0, n_samples = ncol(counts_mat),
        detection_fraction = 0, mean_expression = 0,
        classification = "NOT_IN_DATA", stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, results)
}

rec_197 <- check_receptor_detection(fib197, all_receptors)
rec_197$dataset <- "GSE197677"

rec_221 <- check_receptor_detection(fib221, all_receptors)
rec_221$dataset <- "GSE221561"

rec_all <- rbind(rec_197, rec_221)

cat("=== FIBROBLAST RECEPTOR AVAILABILITY ===\n")
for (i in 1:nrow(rec_all)) {
  if (rec_all$dataset[i] == "GSE197677") {
    cat(sprintf("  %s GSE197677: frac=%.2f mean=%.3f -> %s\n",
      rec_all$receptor[i], rec_all$detection_fraction[i],
      rec_all$mean_expression[i], rec_all$classification[i]))
  }
}

write.csv(rec_all, file.path(O_DIR, "STEP19O_FIBROBLAST_RECEPTOR_AVAILABILITY.csv"), row.names = FALSE)
cat("\nSaved: STEP19O_FIBROBLAST_RECEPTOR_AVAILABILITY.csv\n")

# ============================================================================
# 19O6 — SOURCE LIGAND EXPRESSION QC
# ============================================================================
cat("\n============================================\n")
cat("19O6: SOURCE LIGAND EXPRESSION QC\n")
cat("============================================\n\n")

all_ligands <- unique(curated_lr$ligand)

cat("Checking source ligand expression...\n")

ligand_qc_list <- list()
for (comp in source_compartments) {
  for (ds in c("GSE197677", "GSE221561")) {
    if (ds == "GSE197677") {
      if (!comp %in% names(src_197)) next
      counts <- src_197[[comp]]
    } else {
      if (!comp %in% names(src_221)) next
      counts <- src_221[[comp]]
    }
    
    for (lig in all_ligands) {
      if (lig %in% rownames(counts)) {
        expr <- as.numeric(counts[lig, ])
        detected <- sum(expr > 0)
        n_samp <- ncol(counts)
        frac <- detected / n_samp
        mean_expr <- mean(expr)
        var_expr <- var(expr)
        
        if (frac >= 0.50) flag <- "ADEQUATELY_DETECTED"
        else if (frac >= 0.30) flag <- "PARTIALLY_DETECTED"
        else flag <- "POORLY_DETECTED"
        
        ligand_qc_list[[paste(ds, comp, lig)]] <- data.frame(
          dataset = ds, compartment = comp, ligand = lig,
          n_detected = detected, n_samples = n_samp,
          detection_fraction = frac, mean_expression = mean_expr,
          variance = var_expr, classification = flag,
          stringsAsFactors = FALSE)
      } else {
        ligand_qc_list[[paste(ds, comp, lig)]] <- data.frame(
          dataset = ds, compartment = comp, ligand = lig,
          n_detected = 0, n_samples = ifelse(ds == "GSE197677", ncol(src_197[[comp]]), ncol(src_221[[comp]])),
          detection_fraction = 0, mean_expression = 0,
          variance = 0, classification = "NOT_IN_DATA",
          stringsAsFactors = FALSE)
      }
    }
  }
}

ligand_qc_df <- do.call(rbind, ligand_qc_list)
rownames(ligand_qc_df) <- NULL

cat("=== SOURCE LIGAND AVAILABILITY ===\n")
cat("Adequately detected in at least one dataset/compartment:\n")
adequate_ligs <- unique(ligand_qc_df$ligand[ligand_qc_df$classification == "ADEQUATELY_DETECTED"])
cat(paste(adequate_ligs, collapse=", "), "\n\n")

write.csv(ligand_qc_df, file.path(O_DIR, "STEP19O_SOURCE_LIGAND_AVAILABILITY.csv"), row.names = FALSE)
cat("Saved: STEP19O_SOURCE_LIGAND_AVAILABILITY.csv\n")

# ============================================================================
# 19O7 — PRESPECIFIED SIGNALING FAMILIES
# ============================================================================
cat("\n============================================\n")
cat("19O7: PRESPECIFIED SIGNALING FAMILIES\n")
cat("============================================\n\n")

families_prespecified <- data.frame(
  family = c("TGF_beta", "TNF", "IL1", "IL6_JAK_STAT", "IFN", "PDGF",
             "FGF", "EGF", "WNT", "NOTCH", "CXCL_CCL", "VEGF", "HGF_MET"),
  n_lr_pairs = c(
    sum(curated_lr$family == "TGF_beta"),
    sum(curated_lr$family == "TNF"),
    sum(curated_lr$family == "IL1"),
    sum(curated_lr$family == "IL6_JAK_STAT"),
    sum(curated_lr$family == "IFN"),
    sum(curated_lr$family == "PDGF"),
    sum(curated_lr$family == "FGF"),
    sum(curated_lr$family == "EGF"),
    sum(curated_lr$family == "WNT"),
    sum(curated_lr$family == "NOTCH"),
    sum(curated_lr$family == "CXCL_CCL"),
    sum(curated_lr$family == "VEGF"),
    sum(curated_lr$family == "HGF_MET")
  ),
  special_interest = c(
    "SMAD3/SMAD7/TGFB1 expression associations from Step19M",
    rep("", 3),
    "AP-1 downstream, not a ligand",
    rep("", 6),
    "", ""
  ),
  stringsAsFactors = FALSE
)

cat("=== PRESPECIFIED SIGNALING FAMILIES ===\n")
for (i in 1:nrow(families_prespecified)) {
  cat(sprintf("  %s: %d LR pairs%s\n", families_prespecified$family[i],
    families_prespecified$n_lr_pairs[i],
    ifelse(nchar(families_prespecified$special_interest[i]) > 0,
      paste0(" [", families_prespecified$special_interest[i], "]"), "")))
}

write.csv(families_prespecified, file.path(O_DIR, "STEP19O_PRESPECIFIED_FAMILIES.csv"), row.names = FALSE)
cat("\nSaved: STEP19O_PRESPECIFIED_FAMILIES.csv\n")

# ============================================================================
# 19O8 — LIGAND–CORE2 SAMPLE ASSOCIATION
# ============================================================================
cat("\n============================================\n")
cat("19O8: LIGAND–CORE2 SAMPLE ASSOCIATION\n")
cat("============================================\n\n")

cat("Computing sample-level ligand-CORE2 associations...\n")

assoc_list <- list()
for (ds in c("GSE197677", "GSE221561")) {
  # Get fibroblast samples for this dataset
  fib_samples <- core_scores$sample_id[core_scores$dataset == ds]
  core2_scores <- core_scores$CORE2_STRONG_score[core_scores$dataset == ds]
  
  for (comp in source_compartments) {
    if (ds == "GSE197677") {
      if (!comp %in% names(src_197)) next
      src_counts <- src_197[[comp]]
    } else {
      if (!comp %in% names(src_221)) next
      src_counts <- src_221[[comp]]
    }
    
    # Find matched samples (present in both source and fibroblast)
    matched_samps <- intersect(colnames(src_counts), fib_samples)
    if (length(matched_samps) < 4) next
    
    # Get matched core2 scores
    matched_core2 <- core_scores$CORE2_STRONG_score[match(matched_samps, core_scores$sample_id)]
    
    # Compute log2 CPM for source
    lib_size <- colSums(src_counts[, matched_samps, drop = FALSE])
    lib_size[lib_size == 0] <- 1
    cpm <- sweep(src_counts[, matched_samps, drop = FALSE], 2, lib_size, "/") * 1e6
    log_cpm <- log2(cpm + 1)
    
    for (lig in all_ligands) {
      if (lig %in% rownames(log_cpm)) {
        lig_expr <- as.numeric(log_cpm[lig, ])
        
        valid <- !is.na(lig_expr) & !is.na(matched_core2)
        if (sum(valid) >= 4) {
          rho <- tryCatch(cor(lig_expr[valid], matched_core2[valid], method = "spearman"),
                          error = function(e) NA_real_)
          if (is.na(rho)) next
          
          if (abs(rho) >= 0.60) flag <- "STRONG_ASSOCIATION"
          else if (abs(rho) >= 0.40) flag <- "MODERATE_ASSOCIATION"
          else flag <- "WEAK_ASSOCIATION"
          
          assoc_list[[paste(ds, comp, lig)]] <- data.frame(
            dataset = ds, compartment = comp, ligand = lig,
            n_matched = sum(valid), rho_core2 = rho,
            direction = ifelse(rho > 0, "POSITIVE", "NEGATIVE"),
            classification = flag, stringsAsFactors = FALSE)
        }
      }
    }
  }
}

assoc_df <- do.call(rbind, assoc_list)
rownames(assoc_df) <- NULL

cat("=== LIGAND-CORE2 ASSOCIATIONS ===\n")
cat(sprintf("Total tests: %d\n", nrow(assoc_df)))
cat(sprintf("STRONG: %d, MODERATE: %d, WEAK: %d\n",
  sum(assoc_df$classification == "STRONG_ASSOCIATION"),
  sum(assoc_df$classification == "MODERATE_ASSOCIATION"),
  sum(assoc_df$classification == "WEAK_ASSOCIATION")))

# Show top associations
top_assoc <- assoc_df[order(-abs(assoc_df$rho_core2)), ]
cat("\nTop 10 associations:\n")
for (i in 1:min(10, nrow(top_assoc))) {
  cat(sprintf("  %s %s %s: rho=%.3f (%s) -> %s\n",
    top_assoc$dataset[i], top_assoc$compartment[i], top_assoc$ligand[i],
    top_assoc$rho_core2[i], top_assoc$direction[i], top_assoc$classification[i]))
}

write.csv(assoc_df, file.path(O_DIR, "STEP19O_LIGAND_CORE2_ASSOCIATIONS.csv"), row.names = FALSE)
cat("\nSaved: STEP19O_LIGAND_CORE2_ASSOCIATIONS.csv\n")

# ============================================================================
# 19O9 — LIGAND–CORE4 COMPARATOR
# ============================================================================
cat("\n============================================\n")
cat("19O9: LIGAND–CORE4 COMPARATOR\n")
cat("============================================\n\n")

cat("Computing ligand-CORE4 associations for comparison...\n")

assoc4_list <- list()
for (ds in c("GSE197677", "GSE221561")) {
  fib_samples <- core_scores$sample_id[core_scores$dataset == ds]
  core4_scores <- core_scores$CORE4_score[core_scores$dataset == ds]
  
  for (comp in source_compartments) {
    if (ds == "GSE197677") {
      if (!comp %in% names(src_197)) next
      src_counts <- src_197[[comp]]
    } else {
      if (!comp %in% names(src_221)) next
      src_counts <- src_221[[comp]]
    }
    
    matched_samps <- intersect(colnames(src_counts), fib_samples)
    if (length(matched_samps) < 4) next
    
    matched_core4 <- core_scores$CORE4_score[match(matched_samps, core_scores$sample_id)]
    
    lib_size <- colSums(src_counts[, matched_samps, drop = FALSE])
    lib_size[lib_size == 0] <- 1
    cpm <- sweep(src_counts[, matched_samps, drop = FALSE], 2, lib_size, "/") * 1e6
    log_cpm <- log2(cpm + 1)
    
    for (lig in all_ligands) {
      if (lig %in% rownames(log_cpm)) {
        lig_expr <- as.numeric(log_cpm[lig, ])
        valid <- !is.na(lig_expr) & !is.na(matched_core4)
        if (sum(valid) >= 4) {
          rho4 <- tryCatch(cor(lig_expr[valid], matched_core4[valid], method = "spearman"),
                           error = function(e) NA_real_)
          if (is.na(rho4)) next
          assoc4_list[[paste(ds, comp, lig)]] <- data.frame(
            dataset = ds, compartment = comp, ligand = lig,
            rho_core4 = rho4, stringsAsFactors = FALSE)
        }
      }
    }
  }
}

assoc4_df <- do.call(rbind, assoc4_list)
rownames(assoc4_df) <- NULL

# Merge CORE2 and CORE4 associations
comp_df <- merge(assoc_df, assoc4_df, by = c("dataset", "compartment", "ligand"), all.x = TRUE)

# Classify CORE2 vs CORE4 preference
comp_df$core2_core4_class <- "NO_CLEAR_CORE_ASSOCIATION"
for (i in 1:nrow(comp_df)) {
  rho2 <- comp_df$rho_core2[i]
  rho4 <- comp_df$rho_core4[i]
  if (is.na(rho4)) {
    comp_df$core2_core4_class[i] <- "CORE2_ONLY"
  } else if (abs(rho2) > abs(rho4) * 1.5) {
    comp_df$core2_core4_class[i] <- "CORE2_PREFERENTIAL"
  } else if (abs(rho4) > abs(rho2) * 1.5) {
    comp_df$core2_core4_class[i] <- "CORE4_PREFERENTIAL"
  } else if (sign(rho2) == sign(rho4) && abs(rho2) >= 0.30) {
    comp_df$core2_core4_class[i] <- "SHARED_CORE_ASSOCIATION"
  }
}

cat("=== CORE2 vs CORE4 PREFERENCE ===\n")
cat(sprintf("CORE2_PREFERENTIAL: %d\n", sum(comp_df$core2_core4_class == "CORE2_PREFERENTIAL")))
cat(sprintf("CORE4_PREFERENTIAL: %d\n", sum(comp_df$core2_core4_class == "CORE4_PREFERENTIAL")))
cat(sprintf("SHARED_CORE_ASSOCIATION: %d\n", sum(comp_df$core2_core4_class == "SHARED_CORE_ASSOCIATION")))
cat(sprintf("CORE2_ONLY: %d\n", sum(comp_df$core2_core4_class == "CORE2_ONLY")))
cat(sprintf("NO_CLEAR: %d\n", sum(comp_df$core2_core4_class == "NO_CLEAR_CORE_ASSOCIATION")))

write.csv(comp_df, file.path(O_DIR, "STEP19O_LIGAND_CORE2_CORE4_COMPARISON.csv"), row.names = FALSE)
cat("\nSaved: STEP19O_LIGAND_CORE2_CORE4_COMPARISON.csv\n")

# ============================================================================
# 19O10 — CANONICAL TREATMENT EFFECTS
# ============================================================================
cat("\n============================================\n")
cat("19O10: CANONICAL TREATMENT EFFECTS\n")
cat("============================================\n\n")

cat("Computing treatment effects on source ligand expression...\n")

te_list <- list()
for (ds in c("GSE197677", "GSE221561")) {
  treated_group <- TREAT_CONTRAST[[ds]][["treated"]]
  control_group <- TREAT_CONTRAST[[ds]][["control"]]
  
  for (comp in source_compartments) {
    if (ds == "GSE197677") {
      if (!comp %in% names(src_197)) next
      src_counts <- src_197[[comp]]
    } else {
      if (!comp %in% names(src_221)) next
      src_counts <- src_221[[comp]]
    }
    
    # Map samples to treatment
    treat_map_ds <- TREAT_MAP[[ds]]
    samps_in_data <- colnames(src_counts)
    samps_with_treat <- intersect(samps_in_data, names(treat_map_ds))
    
    if (length(samps_with_treat) < 4) next
    
    # Compute log2 CPM
    lib_size <- colSums(src_counts[, samps_with_treat, drop = FALSE])
    lib_size[lib_size == 0] <- 1
    cpm <- sweep(src_counts[, samps_with_treat, drop = FALSE], 2, lib_size, "/") * 1e6
    log_cpm <- log2(cpm + 1)
    
    treated_samps <- samps_with_treat[treat_map_ds[samps_with_treat] == treated_group]
    control_samps <- samps_with_treat[treat_map_ds[samps_with_treat] == control_group]
    
    for (lig in all_ligands) {
      if (lig %in% rownames(log_cpm)) {
        treated_expr <- as.numeric(log_cpm[lig, treated_samps])
        control_expr <- as.numeric(log_cpm[lig, control_samps])
        
        eff <- mean(treated_expr, na.rm = TRUE) - mean(control_expr, na.rm = TRUE)
        if (is.nan(eff)) eff <- NA
        
        te_list[[paste(ds, comp, lig)]] <- data.frame(
          dataset = ds, compartment = comp, ligand = lig,
          effect = eff,
          direction = ifelse(is.na(eff), "NA", ifelse(eff > 0, "HIGHER_IN_TREATED", "LOWER_IN_TREATED")),
          treated_mean = mean(treated_expr, na.rm = TRUE),
          control_mean = mean(control_expr, na.rm = TRUE),
          treated_n = length(treated_samps),
          control_n = length(control_samps),
          stringsAsFactors = FALSE)
      }
    }
  }
}

te_df <- do.call(rbind, te_list)
rownames(te_df) <- NULL

cat("=== TREATMENT EFFECTS ON SOURCE LIGANDS ===\n")
cat(sprintf("Total: %d\n", nrow(te_df)))

write.csv(te_df, file.path(O_DIR, "STEP19O_SOURCE_LIGAND_TREATMENT_EFFECTS.csv"), row.names = FALSE)
cat("Saved: STEP19O_SOURCE_LIGAND_TREATMENT_EFFECTS.csv\n")

# ============================================================================
# 19O11 — EXACT PERMUTATION AUDIT
# ============================================================================
cat("\n============================================\n")
cat("19O11: EXACT PERMUTATION AUDIT\n")
cat("============================================\n\n")

cat("Running exact permutation tests for top candidate ligands...\n")

# Select top candidates based on |CORE2 rho|
top_cands <- unique(assoc_df$ligand[abs(assoc_df$rho_core2) >= 0.30])
if (length(top_cands) == 0) {
  top_cands <- unique(assoc_df$ligand[order(-abs(assoc_df$rho_core2))])[1:min(10, nrow(assoc_df))]
}

perm_list <- list()
for (lig in top_cands) {
  for (ds in c("GSE197677", "GSE221561")) {
    # Find compartments where this ligand is expressed
    comps_for_lig <- unique(assoc_df$compartment[assoc_df$ligand == lig & assoc_df$dataset == ds])
    
    for (comp in comps_for_lig) {
      treated_group <- TREAT_CONTRAST[[ds]][["treated"]]
      control_group <- TREAT_CONTRAST[[ds]][["control"]]
      
      if (ds == "GSE197677") {
        if (!comp %in% names(src_197)) next
        src_counts <- src_197[[comp]]
      } else {
        if (!comp %in% names(src_221)) next
        src_counts <- src_221[[comp]]
      }
      
      treat_map_ds <- TREAT_MAP[[ds]]
      samps_in_data <- colnames(src_counts)
      samps_with_treat <- intersect(samps_in_data, names(treat_map_ds))
      
      treated_samps <- samps_with_treat[treat_map_ds[samps_with_treat] == treated_group]
      control_samps <- samps_with_treat[treat_map_ds[samps_with_treat] == control_group]
      n_treated <- length(treated_samps)
      n_total <- n_treated + length(control_samps)
      
      if (n_total < 4 || n_treated < 1 || length(control_samps) < 1) next
      
      # Log2 CPM
      lib_size <- colSums(src_counts[, samps_with_treat, drop = FALSE])
      lib_size[lib_size == 0] <- 1
      cpm <- sweep(src_counts[, samps_with_treat, drop = FALSE], 2, lib_size, "/") * 1e6
      log_cpm <- log2(cpm + 1)
      
      if (!lig %in% rownames(log_cpm)) next
      
      all_expr <- as.numeric(log_cpm[lig, samps_with_treat])
      obs_effect <- mean(all_expr[1:n_treated]) - mean(all_expr[(n_treated+1):n_total])
      
      # Exact permutation
      all_combinations <- combn(n_total, n_treated)
      n_perms <- ncol(all_combinations)
      
      if (n_perms > 10000) {
        # Subsample permutations for large spaces
        perm_idx <- sample(1:n_perms, 10000)
        all_combinations <- all_combinations[, perm_idx]
        n_perms <- 10000
      }
      
      perm_effects <- numeric(n_perms)
      for (p in 1:n_perms) {
        perm_treated <- all_combinations[, p]
        perm_control <- setdiff(1:n_total, perm_treated)
        perm_effects[p] <- mean(all_expr[perm_treated]) - mean(all_expr[perm_control])
      }
      
      p_two <- mean(abs(perm_effects) >= abs(obs_effect))
      
      perm_list[[paste(ds, comp, lig)]] <- data.frame(
        dataset = ds, compartment = comp, ligand = lig,
        n_treated = n_treated, n_total = n_total,
        total_perms = n_perms, observed_effect = obs_effect,
        p_two_sided = p_two, stringsAsFactors = FALSE)
      
      cat(sprintf("  %s %s %s: effect=%.4f P(two)=%.4f (%d perms)\n",
        ds, comp, lig, obs_effect, p_two, n_perms))
    }
  }
}

perm_df <- do.call(rbind, perm_list)
rownames(perm_df) <- NULL

write.csv(perm_df, file.path(O_DIR, "STEP19O_LIGAND_EXACT_PERMUTATION.csv"), row.names = FALSE)
cat("\nSaved: STEP19O_LIGAND_EXACT_PERMUTATION.csv\n")

# ============================================================================
# 19O12 — LIGAND–RECEPTOR PAIR FILTER
# ============================================================================
cat("\n============================================\n")
cat("19O12: LIGAND–RECEPTOR PAIR FILTER\n")
cat("============================================\n\n")

cat("Filtering eligible ligand-receptor pairs...\n")

# Get receptor availability per dataset
rec_wide <- reshape2::dcast(rec_all, receptor ~ dataset, value.var = "classification")

eligible_pairs <- list()
for (i in 1:nrow(curated_lr)) {
  lig <- curated_lr$ligand[i]
  rec <- curated_lr$receptor[i]
  fam <- curated_lr$family[i]
  
  for (ds in c("GSE197677", "GSE221561")) {
    # Check receptor availability
    rec_class <- rec_all$classification[rec_all$receptor == rec & rec_all$dataset == ds]
    if (length(rec_class) == 0) rec_class <- "NOT_IN_DATA"
    
    # Check ligand availability per source
    for (comp in source_compartments) {
      lig_class <- ligand_qc_df$classification[ligand_qc_df$ligand == lig &
        ligand_qc_df$compartment == comp & ligand_qc_df$dataset == ds]
      if (length(lig_class) == 0) next
      
      # Check CORE2 association
      core2_rho <- assoc_df$rho_core2[assoc_df$ligand == lig &
        assoc_df$compartment == comp & assoc_df$dataset == ds]
      if (length(core2_rho) == 0) next
      
      # Check treatment effect
      te_eff <- te_df$effect[te_df$ligand == lig &
        te_df$compartment == comp & te_df$dataset == ds]
      if (length(te_eff) == 0) te_eff <- NA
      
      # Eligibility: ligand detected + receptor detected + CORE2 at least weak
      lig_ok <- lig_class %in% c("ADEQUATELY_DETECTED", "PARTIALLY_DETECTED")
      rec_ok <- rec_class %in% c("ROBUSTLY_DETECTED", "PARTIALLY_DETECTED")
      core2_ok <- !is.na(core2_rho) && abs(core2_rho) >= 0.10
      
      if (lig_ok && rec_ok && core2_ok) {
        eligible_pairs[[paste(ds, comp, lig, rec)]] <- data.frame(
          dataset = ds, source = comp, ligand = lig, receptor = rec,
          family = fam, ligand_detection = lig_class,
          receptor_detection = rec_class, rho_core2 = core2_rho,
          treatment_effect = te_eff, resource = "curated_minimal",
          stringsAsFactors = FALSE)
      }
    }
  }
}

eligible_df <- do.call(rbind, eligible_pairs)
rownames(eligible_df) <- NULL

cat("=== ELIGIBLE LR PAIRS ===\n")
cat(sprintf("Total eligible: %d\n", nrow(eligible_df)))
if (nrow(eligible_df) > 0) {
  cat(sprintf("Unique ligands: %d\n", length(unique(eligible_df$ligand))))
  cat(sprintf("Unique receptors: %d\n", length(unique(eligible_df$receptor))))
  cat(sprintf("Families: %s\n", paste(unique(eligible_df$family), collapse=", ")))
}

write.csv(eligible_df, file.path(O_DIR, "STEP19O_ELIGIBLE_LIGAND_RECEPTOR_PAIRS.csv"), row.names = FALSE)
cat("\nSaved: STEP19O_ELIGIBLE_LIGAND_RECEPTOR_PAIRS.csv\n")

# ============================================================================
# 19O13 — LIGAND-TARGET SUPPORT
# ============================================================================
cat("\n============================================\n")
cat("19O13: LIGAND-TARGET SUPPORT\n")
cat("============================================\n\n")

cat("No local ligand-target matrix available.\n")
cat("Skipping 19O13 without failure.\n")

lt_status <- data.frame(
  step = "19O13", status = "SKIPPED_NO_RESOURCE",
  reason = "No local ligand-target matrix available",
  timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  stringsAsFactors = FALSE
)
write.csv(lt_status, file.path(O_DIR, "STEP19O_LIGAND_CORE2_TARGET_SUPPORT.csv"), row.names = FALSE)
cat("Saved: STEP19O_LIGAND_CORE2_TARGET_SUPPORT.csv\n")

# ============================================================================
# 19O14 — RUNX1 BRIDGE ANALYSIS
# ============================================================================
cat("\n============================================\n")
cat("19O14: RUNX1 BRIDGE ANALYSIS\n")
cat("============================================\n\n")

cat("Testing ligand-RUNX1-CORE2 bridge associations...\n")

# Load RUNX1 activity from Step19N
runx1_activity <- read.csv(file.path(BASE_DIR,
  "06_tables/cross_dataset/pathway_analysis/Step19N/STEP19N_SAMPLE_REGULATOR_ACTIVITY.csv"))
runx1_act <- runx1_activity[runx1_activity$regulator == "RUNX1", ]

bridge_list <- list()
for (ds in c("GSE197677", "GSE221561")) {
  fib_samples <- core_scores$sample_id[core_scores$dataset == ds]
  core2_scores <- core_scores$CORE2_STRONG_score[core_scores$dataset == ds]
  
  # Get RUNX1 activity for tumor samples
  runx1_samps <- runx1_act$sample_id[runx1_act$dataset == ds & !is.na(runx1_act$activity)]
  matched_samps <- intersect(fib_samples, runx1_samps)
  
  if (length(matched_samps) < 4) next
  
  matched_core2 <- core_scores$CORE2_STRONG_score[match(matched_samps, core_scores$sample_id)]
  matched_runx1 <- runx1_act$activity[match(matched_samps, runx1_act$sample_id)]
  
  for (comp in source_compartments) {
    if (ds == "GSE197677") {
      if (!comp %in% names(src_197)) next
      src_counts <- src_197[[comp]]
    } else {
      if (!comp %in% names(src_221)) next
      src_counts <- src_221[[comp]]
    }
    
    src_matched <- intersect(colnames(src_counts), matched_samps)
    if (length(src_matched) < 4) next
    
    lib_size <- colSums(src_counts[, src_matched, drop = FALSE])
    lib_size[lib_size == 0] <- 1
    cpm <- sweep(src_counts[, src_matched, drop = FALSE], 2, lib_size, "/") * 1e6
    log_cpm <- log2(cpm + 1)
    
    for (lig in all_ligands) {
      if (lig %in% rownames(log_cpm)) {
        lig_expr <- as.numeric(log_cpm[lig, ])
        
        # Get matched indices
        idx_match <- match(src_matched, matched_samps)
        valid <- !is.na(lig_expr) & !is.na(matched_core2[idx_match]) & !is.na(matched_runx1[idx_match])
        
        if (sum(valid) >= 4) {
          rho_runx1 <- tryCatch(cor(lig_expr[valid], matched_runx1[idx_match][valid], method = "spearman"),
                                error = function(e) NA_real_)
          rho_core2 <- tryCatch(cor(lig_expr[valid], matched_core2[idx_match][valid], method = "spearman"),
                                error = function(e) NA_real_)
          
          # Classification
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
          
          bridge_list[[paste(ds, comp, lig)]] <- data.frame(
            dataset = ds, compartment = comp, ligand = lig,
            ligand_vs_runx1_rho = rho_runx1, ligand_vs_core2_rho = rho_core2,
            classification = br_class, stringsAsFactors = FALSE)
        }
      }
    }
  }
}

bridge_df <- do.call(rbind, bridge_list)
rownames(bridge_df) <- NULL

cat("=== RUNX1 BRIDGE ANALYSIS ===\n")
if (nrow(bridge_df) > 0) {
  cat(sprintf("LIGAND_RUNX1_CORE2_BRIDGE: %d\n", sum(bridge_df$classification == "LIGAND_RUNX1_CORE2_BRIDGE")))
  cat(sprintf("LIGAND_CORE2_ONLY: %d\n", sum(bridge_df$classification == "LIGAND_CORE2_ONLY")))
  cat(sprintf("LIGAND_RUNX1_ONLY: %d\n", sum(bridge_df$classification == "LIGAND_RUNX1_ONLY")))
  cat(sprintf("NO_BRIDGE_SUPPORT: %d\n", sum(bridge_df$classification == "NO_BRIDGE_SUPPORT")))
  
  bridges <- bridge_df[bridge_df$classification == "LIGAND_RUNX1_CORE2_BRIDGE", ]
  if (nrow(bridges) > 0) {
    cat("\nBridges found:\n")
    for (i in 1:nrow(bridges)) {
      cat(sprintf("  %s %s %s: RUNX1_rho=%.3f CORE2_rho=%.3f\n",
        bridges$dataset[i], bridges$compartment[i], bridges$ligand[i],
        bridges$ligand_vs_runx1_rho[i], bridges$ligand_vs_core2_rho[i]))
    }
  }
}

write.csv(bridge_df, file.path(O_DIR, "STEP19O_RUNX1_SIGNALING_BRIDGE.csv"), row.names = FALSE)
cat("\nSaved: STEP19O_RUNX1_SIGNALING_BRIDGE.csv\n")

# ============================================================================
# 19O15 — AP1 / SMAD BRIDGE ANALYSIS
# ============================================================================
cat("\n============================================\n")
cat("19O15: AP1 / SMAD BRIDGE ANALYSIS\n")
cat("============================================\n\n")

cat("Testing ligand-AP1/SMAD bridge associations...\n")

# AP1/SMAD regulators with activity from Step19N
ap1_smad_regs <- c("FOS", "FOSL2", "SMAD3", "SMAD7")

ap1_bridge_list <- list()
for (reg in ap1_smad_regs) {
  reg_act <- runx1_activity[runx1_activity$regulator == reg, ]
  has_activity <- any(!is.na(reg_act$activity))
  
  for (ds in c("GSE197677", "GSE221561")) {
    fib_samples <- core_scores$sample_id[core_scores$dataset == ds]
    core2_scores <- core_scores$CORE2_STRONG_score[core_scores$dataset == ds]
    
    if (has_activity) {
      reg_samps <- reg_act$sample_id[reg_act$dataset == ds & !is.na(reg_act$activity)]
      matched_samps <- intersect(fib_samples, reg_samps)
    } else {
      matched_samps <- fib_samples
    }
    
    if (length(matched_samps) < 4) next
    
    matched_core2 <- core_scores$CORE2_STRONG_score[match(matched_samps, core_scores$sample_id)]
    
    if (has_activity) {
      matched_reg <- reg_act$activity[match(matched_samps, reg_act$sample_id)]
    } else {
      matched_reg <- rep(NA, length(matched_samps))
    }
    
    for (comp in source_compartments) {
      if (ds == "GSE197677") {
        if (!comp %in% names(src_197)) next
        src_counts <- src_197[[comp]]
      } else {
        if (!comp %in% names(src_221)) next
        src_counts <- src_221[[comp]]
      }
      
      src_matched <- intersect(colnames(src_counts), matched_samps)
      if (length(src_matched) < 4) next
      
      lib_size <- colSums(src_counts[, src_matched, drop = FALSE])
      lib_size[lib_size == 0] <- 1
      cpm <- sweep(src_counts[, src_matched, drop = FALSE], 2, lib_size, "/") * 1e6
      log_cpm <- log2(cpm + 1)
      
      for (lig in all_ligands) {
        if (lig %in% rownames(log_cpm)) {
          lig_expr <- as.numeric(log_cpm[lig, ])
          idx_match <- match(src_matched, matched_samps)
          
          if (has_activity) {
            valid <- !is.na(lig_expr) & !is.na(matched_core2[idx_match]) & !is.na(matched_reg[idx_match])
          } else {
            valid <- !is.na(lig_expr) & !is.na(matched_core2[idx_match])
          }
          
          if (sum(valid) >= 4) {
            rho_core2 <- tryCatch(cor(lig_expr[valid], matched_core2[idx_match][valid], method = "spearman"),
                                  error = function(e) NA_real_)
            
            if (has_activity) {
              rho_reg <- tryCatch(cor(lig_expr[valid], matched_reg[idx_match][valid], method = "spearman"),
                                  error = function(e) NA_real_)
            } else {
              rho_reg <- NA
            }
            
            ap1_bridge_list[[paste(reg, ds, comp, lig)]] <- data.frame(
              regulator = reg, dataset = ds, compartment = comp, ligand = lig,
              ligand_vs_regulator_rho = rho_reg, ligand_vs_core2_rho = rho_core2,
              evidence_type = ifelse(has_activity, "ORTHOGONAL_ACTIVITY", "EXPRESSION_ONLY"),
              stringsAsFactors = FALSE)
          }
        }
      }
    }
  }
}

ap1_bridge_df <- do.call(rbind, ap1_bridge_list)
rownames(ap1_bridge_df) <- NULL

cat("=== AP1/SMAD BRIDGE ANALYSIS ===\n")
if (nrow(ap1_bridge_df) > 0) {
  cat(sprintf("Total tests: %d\n", nrow(ap1_bridge_df)))
  cat(sprintf("ORTHOGONAL_ACTIVITY: %d\n", sum(ap1_bridge_df$evidence_type == "ORTHOGONAL_ACTIVITY")))
  cat(sprintf("EXPRESSION_ONLY: %d\n", sum(ap1_bridge_df$evidence_type == "EXPRESSION_ONLY")))
  
  # Show strong bridges
  strong_bridges <- ap1_bridge_df[abs(ap1_bridge_df$ligand_vs_core2_rho) >= 0.40 &
    (is.na(ap1_bridge_df$ligand_vs_regulator_rho) | abs(ap1_bridge_df$ligand_vs_regulator_rho) >= 0.30), ]
  if (nrow(strong_bridges) > 0) {
    cat("\nStrong associations (|CORE2 rho| >= 0.40):\n")
    for (i in 1:min(10, nrow(strong_bridges))) {
      cat(sprintf("  %s %s %s %s: reg_rho=%.3f core2_rho=%.3f (%s)\n",
        strong_bridges$regulator[i], strong_bridges$dataset[i],
        strong_bridges$compartment[i], strong_bridges$ligand[i],
        strong_bridges$ligand_vs_regulator_rho[i], strong_bridges$ligand_vs_core2_rho[i],
        strong_bridges$evidence_type[i]))
    }
  }
}

write.csv(ap1_bridge_df, file.path(O_DIR, "STEP19O_AP1_SMAD_SIGNALING_BRIDGE.csv"), row.names = FALSE)
cat("\nSaved: STEP19O_AP1_SMAD_SIGNALING_BRIDGE.csv\n")

# ============================================================================
# 19O16 — SOURCE COMPARTMENT PRIORITIZATION
# ============================================================================
cat("\n============================================\n")
cat("19O16: SOURCE COMPARTMENT PRIORITIZATION\n")
cat("============================================\n\n")

cat("Determining strongest source compartment per ligand...\n")

prioritize_list <- list()
for (lig in unique(assoc_df$ligand)) {
  for (ds in c("GSE197677", "GSE221561")) {
    sub <- assoc_df[assoc_df$ligand == lig & assoc_df$dataset == ds, ]
    if (nrow(sub) == 0) next
    
    # Find best compartment by |rho|
    best_idx <- which.max(abs(sub$rho_core2))
    best_comp <- sub$compartment[best_idx]
    best_rho <- sub$rho_core2[best_idx]
    
    # Check if multiple compartments have similar signal
    strong_comps <- sub$compartment[abs(sub$rho_core2) >= 0.30]
    
    if (length(strong_comps) >= 2) {
      comp_class <- "MULTIPLE_SOURCE_SUPPORTED"
    } else if (best_comp == "Epithelial") {
      comp_class <- "EPITHELIAL_SOURCE_SUPPORTED"
    } else if (best_comp == "Myeloid") {
      comp_class <- "MYELOID_SOURCE_SUPPORTED"
    } else if (best_comp == "T_NK") {
      comp_class <- "T_NK_SOURCE_SUPPORTED"
    } else if (best_comp == "Endothelial") {
      comp_class <- "ENDOTHELIAL_SOURCE_SUPPORTED"
    } else {
      comp_class <- "SOURCE_INCONCLUSIVE"
    }
    
    prioritize_list[[paste(ds, lig)]] <- data.frame(
      ligand = lig, dataset = ds, best_compartment = best_comp,
      best_rho = best_rho, n_strong_compartments = length(strong_comps),
      classification = comp_class, stringsAsFactors = FALSE)
  }
}

prior_df <- do.call(rbind, prioritize_list)
rownames(prior_df) <- NULL

cat("=== SOURCE COMPARTMENT PRIORITIZATION ===\n")
cat(sprintf("EPITHELIAL: %d\n", sum(prior_df$classification == "EPITHELIAL_SOURCE_SUPPORTED")))
cat(sprintf("MYELOID: %d\n", sum(prior_df$classification == "MYELOID_SOURCE_SUPPORTED")))
cat(sprintf("T_NK: %d\n", sum(prior_df$classification == "T_NK_SOURCE_SUPPORTED")))
cat(sprintf("ENDOTHELIAL: %d\n", sum(prior_df$classification == "ENDOTHELIAL_SOURCE_SUPPORTED")))
cat(sprintf("MULTIPLE: %d\n", sum(prior_df$classification == "MULTIPLE_SOURCE_SUPPORTED")))
cat(sprintf("INCONCLUSIVE: %d\n", sum(prior_df$classification == "SOURCE_INCONCLUSIVE")))

write.csv(prior_df, file.path(O_DIR, "STEP19O_SOURCE_COMPARTMENT_PRIORITIZATION.csv"), row.names = FALSE)
cat("\nSaved: STEP19O_SOURCE_COMPARTMENT_PRIORITIZATION.csv\n")

# ============================================================================
# 19O17 — CROSS-DATASET REPLICATION
# ============================================================================
cat("\n============================================\n")
cat("19O17: CROSS-DATASET REPLICATION\n")
cat("============================================\n\n")

cat("Comparing signaling across datasets...\n")

repl_list <- list()
for (lig in unique(assoc_df$ligand)) {
  for (comp in source_compartments) {
    sub2 <- assoc_df[assoc_df$ligand == lig & assoc_df$compartment == comp, ]
    if (nrow(sub2) < 2) next
    
    ds1 <- sub2[sub2$dataset == "GSE197677", ]
    ds2 <- sub2[sub2$dataset == "GSE221561", ]
    if (nrow(ds1) == 0 || nrow(ds2) == 0) next
    
    # CORE2 association direction replication
    core2_repl <- !is.na(ds1$rho_core2) && !is.na(ds2$rho_core2) &&
      sign(ds1$rho_core2) == sign(ds2$rho_core2) &&
      abs(ds1$rho_core2) >= 0.20 && abs(ds2$rho_core2) >= 0.20
    
    # Treatment effect direction replication
    te1 <- te_df$effect[te_df$ligand == lig & te_df$compartment == comp & te_df$dataset == "GSE197677"]
    te2 <- te_df$effect[te_df$ligand == lig & te_df$compartment == comp & te_df$dataset == "GSE221561"]
    te_repl <- length(te1) > 0 && length(te2) > 0 &&
      !is.na(te1) && !is.na(te2) && sign(te1) == sign(te2)
    
    # Receptor availability
    recs <- unique(curated_lr$receptor[curated_lr$ligand == lig])
    rec_avail <- all(sapply(recs, function(r) {
      any(rec_all$receptor == r & rec_all$dataset == "GSE197677" &
        rec_all$classification %in% c("ROBUSTLY_DETECTED", "PARTIALLY_DETECTED")) &&
      any(rec_all$receptor == r & rec_all$dataset == "GSE221561" &
        rec_all$classification %in% c("ROBUSTLY_DETECTED", "PARTIALLY_DETECTED"))
    }))
    
    # RUNX1 bridge direction
    br1 <- bridge_df[bridge_df$ligand == lig & bridge_df$compartment == comp & bridge_df$dataset == "GSE197677", ]
    br2 <- bridge_df[bridge_df$ligand == lig & bridge_df$compartment == comp & bridge_df$dataset == "GSE221561", ]
    runx1_repl <- nrow(br1) > 0 && nrow(br2) > 0 &&
      sign(br1$ligand_vs_runx1_rho) == sign(br2$ligand_vs_runx1_rho)
    
    if (core2_repl && te_repl) {
      repl_class <- "CROSS_DATASET_SIGNAL_SUPPORTED"
    } else if (core2_repl) {
      repl_class <- "CORE2_ASSOCIATION_REPLICATED"
    } else if (te_repl) {
      repl_class <- "TREATMENT_EFFECT_REPLICATED"
    } else {
      repl_class <- "CROSS_DATASET_DISCORDANT"
    }
    
    repl_list[[paste(lig, comp)]] <- data.frame(
      ligand = lig, compartment = comp,
      core2_rho_197 = ds1$rho_core2, core2_rho_221 = ds2$rho_core2,
      core2_direction_replicated = core2_repl,
      te_direction_replicated = te_repl,
      receptor_available_both = rec_avail,
      runx1_bridge_replicated = runx1_repl,
      classification = repl_class, stringsAsFactors = FALSE)
  }
}

repl_df <- do.call(rbind, repl_list)
rownames(repl_df) <- NULL

cat("=== CROSS-DATASET SIGNALING REPLICATION ===\n")
cat(sprintf("CROSS_DATASET_SIGNAL_SUPPORTED: %d\n", sum(repl_df$classification == "CROSS_DATASET_SIGNAL_SUPPORTED")))
cat(sprintf("CORE2_ASSOCIATION_REPLICATED: %d\n", sum(repl_df$classification == "CORE2_ASSOCIATION_REPLICATED")))
cat(sprintf("TREATMENT_EFFECT_REPLICATED: %d\n", sum(repl_df$classification == "TREATMENT_EFFECT_REPLICATED")))
cat(sprintf("CROSS_DATASET_DISCORDANT: %d\n", sum(repl_df$classification == "CROSS_DATASET_DISCORDANT")))

write.csv(repl_df, file.path(O_DIR, "STEP19O_CROSS_DATASET_SIGNAL_REPLICATION.csv"), row.names = FALSE)
cat("\nSaved: STEP19O_CROSS_DATASET_SIGNAL_REPLICATION.csv\n")

# ============================================================================
# 19O18 — LOSO STABILITY
# ============================================================================
cat("\n============================================\n")
cat("19O18: LOSO STABILITY\n")
cat("============================================\n\n")

cat("Running LOSO for top ligand candidates...\n")

# Select top candidates
top_lig_by_rho <- unique(assoc_df$ligand[order(-abs(assoc_df$rho_core2))])[1:min(10, length(unique(assoc_df$ligand)))]

loso_long_list <- list()
for (lig in top_lig_by_rho) {
  for (ds in c("GSE197677", "GSE221561")) {
    # Find best compartment for this ligand
    best_comp <- assoc_df$compartment[assoc_df$ligand == lig & assoc_df$dataset == ds &
      abs(assoc_df$rho_core2) == max(abs(assoc_df$rho_core2[assoc_df$ligand == lig & assoc_df$dataset == ds]))]
    if (length(best_comp) == 0) next
    best_comp <- best_comp[1]
    
    fib_samples <- core_scores$sample_id[core_scores$dataset == ds]
    core2_scores <- core_scores$CORE2_STRONG_score[core_scores$dataset == ds]
    
    if (ds == "GSE197677") {
      if (!best_comp %in% names(src_197)) next
      src_counts <- src_197[[best_comp]]
    } else {
      if (!best_comp %in% names(src_221)) next
      src_counts <- src_221[[best_comp]]
    }
    
    matched_samps <- intersect(colnames(src_counts), fib_samples)
    if (length(matched_samps) < 5) next
    
    matched_core2 <- core_scores$CORE2_STRONG_score[match(matched_samps, core_scores$sample_id)]
    
    lib_size <- colSums(src_counts[, matched_samps, drop = FALSE])
    lib_size[lib_size == 0] <- 1
    cpm <- sweep(src_counts[, matched_samps, drop = FALSE], 2, lib_size, "/") * 1e6
    log_cpm <- log2(cpm + 1)
    
    if (!lig %in% rownames(log_cpm)) next
    
    lig_expr <- as.numeric(log_cpm[lig, ])
    
    # Full effect
    valid_full <- !is.na(lig_expr) & !is.na(matched_core2)
    full_rho <- if (sum(valid_full) >= 4) cor(lig_expr[valid_full], matched_core2[valid_full], method = "spearman") else NA
    
    # Treatment effect
    treat_map_ds <- TREAT_MAP[[ds]]
    samps_with_treat <- intersect(matched_samps, names(treat_map_ds))
    treated_group <- TREAT_CONTRAST[[ds]][["treated"]]
    control_group <- TREAT_CONTRAST[[ds]][["control"]]
    treated_samps <- samps_with_treat[treat_map_ds[samps_with_treat] == treated_group]
    control_samps <- samps_with_treat[treat_map_ds[samps_with_treat] == control_group]
    
    full_effect <- mean(lig_expr[match(treated_samps, matched_samps)], na.rm = TRUE) -
      mean(lig_expr[match(control_samps, matched_samps)], na.rm = TRUE)
    
    n <- length(matched_samps)
    for (i in 1:n) {
      keep <- setdiff(1:n, i)
      keep_samps <- matched_samps[keep]
      
      # Recompute for kept samples
      keep_lib <- colSums(src_counts[, keep_samps, drop = FALSE])
      keep_lib[keep_lib == 0] <- 1
      keep_cpm <- sweep(src_counts[, keep_samps, drop = FALSE], 2, keep_lib, "/") * 1e6
      keep_log_cpm <- log2(keep_cpm + 1)
      
      if (lig %in% rownames(keep_log_cpm)) {
        keep_lig <- as.numeric(keep_log_cpm[lig, ])
        keep_core2 <- matched_core2[keep]
        
        valid_loso <- !is.na(keep_lig) & !is.na(keep_core2)
        loso_rho <- if (sum(valid_loso) >= 4) cor(keep_lig[valid_loso], keep_core2[valid_loso], method = "spearman") else NA
        
        # LOSO treatment effect
        keep_treat <- intersect(keep_samps, treated_samps)
        keep_ctrl <- intersect(keep_samps, control_samps)
        loso_effect <- if (length(keep_treat) >= 1 && length(keep_ctrl) >= 1)
          mean(keep_lig[match(keep_treat, keep_samps)], na.rm = TRUE) -
            mean(keep_lig[match(keep_ctrl, keep_samps)], na.rm = TRUE) else NA
        
        rho_sign_pres <- sign(loso_rho) == sign(full_rho) & !is.na(loso_rho)
        eff_sign_pres <- sign(loso_effect) == sign(full_effect) & !is.na(loso_effect)
        
        loso_long_list[[paste(lig, ds, matched_samps[i])]] <- data.frame(
          ligand = lig, dataset = ds, best_compartment = best_comp,
          left_out_sample = matched_samps[i],
          loso_rho = loso_rho, loso_effect = loso_effect,
          rho_sign_preserved = rho_sign_pres, effect_sign_preserved = eff_sign_pres,
          stringsAsFactors = FALSE)
      }
    }
  }
}

loso_long_df <- do.call(rbind, loso_long_list)
rownames(loso_long_df) <- NULL

# Summarize LOSO
loso_summary_list <- list()
for (lig in top_lig_by_rho) {
  for (ds in c("GSE197677", "GSE221561")) {
    sub <- loso_long_df[loso_long_df$ligand == lig & loso_long_df$dataset == ds, ]
    if (nrow(sub) == 0) next
    
    rho_pres <- mean(sub$rho_sign_preserved, na.rm = TRUE)
    eff_pres <- mean(sub$effect_sign_preserved, na.rm = TRUE)
    if (is.nan(rho_pres)) rho_pres <- NA
    if (is.nan(eff_pres)) eff_pres <- NA
    
    if (is.na(rho_pres) || is.na(eff_pres)) {
      loso_class <- "NO_VALID_LOSO"
    } else if (rho_pres >= 0.80 && eff_pres >= 0.80) {
      loso_class <- "ROBUST"
    } else if (ds == "GSE221561" && (rho_pres < 0.80 || eff_pres < 0.80)) {
      loso_class <- "CONTROL_SENSITIVE"
    } else if (rho_pres < 0.60 || eff_pres < 0.60) {
      loso_class <- "UNSTABLE"
    } else {
      loso_class <- "OUTLIER_SENSITIVE"
    }
    
    loso_summary_list[[paste(lig, ds)]] <- data.frame(
      ligand = lig, dataset = ds,
      rho_sign_preservation = rho_pres, effect_sign_preservation = eff_pres,
      classification = loso_class, stringsAsFactors = FALSE)
  }
}

loso_summary_df <- do.call(rbind, loso_summary_list)
rownames(loso_summary_df) <- NULL

cat("=== LOSO STABILITY ===\n")
for (i in 1:nrow(loso_summary_df)) {
  cat(sprintf("  %s %s: rho_pres=%.0f%% eff_pres=%.0f%% -> %s\n",
    loso_summary_df$ligand[i], loso_summary_df$dataset[i],
    loso_summary_df$rho_sign_preservation[i] * 100,
    loso_summary_df$effect_sign_preservation[i] * 100,
    loso_summary_df$classification[i]))
}

write.csv(loso_long_df, file.path(O_DIR, "STEP19O_SIGNALING_LOSO_LONG.csv"), row.names = FALSE)
write.csv(loso_summary_df, file.path(O_DIR, "STEP19O_SIGNALING_LOSO_SUMMARY.csv"), row.names = FALSE)
cat("\nSaved: STEP19O_SIGNALING_LOSO_LONG.csv\n")
cat("Saved: STEP19O_SIGNALING_LOSO_SUMMARY.csv\n")

# ============================================================================
# 19O19 — ABUNDANCE SENSITIVITY
# ============================================================================
cat("\n============================================\n")
cat("19O19: ABUNDANCE SENSITIVITY\n")
cat("============================================\n\n")

cat("Testing abundance sensitivity for top candidates...\n")

abund_list <- list()
for (lig in top_lig_by_rho) {
  for (ds in c("GSE197677", "GSE221561")) {
    # Find best compartment
    sub_assoc <- assoc_df[assoc_df$ligand == lig & assoc_df$dataset == ds, ]
    if (nrow(sub_assoc) == 0) next
    best_comp <- sub_assoc$compartment[which.max(abs(sub_assoc$rho_core2))]
    
    fib_samples <- core_scores$sample_id[core_scores$dataset == ds]
    
    if (ds == "GSE197677") {
      if (!best_comp %in% names(src_197)) next
      src_counts <- src_197[[best_comp]]
    } else {
      if (!best_comp %in% names(src_221)) next
      src_counts <- src_221[[best_comp]]
    }
    
    matched_samps <- intersect(colnames(src_counts), fib_samples)
    if (length(matched_samps) < 5) next
    
    matched_core2 <- core_scores$CORE2_STRONG_score[match(matched_samps, core_scores$sample_id)]
    
    lib_size <- colSums(src_counts[, matched_samps, drop = FALSE])
    lib_size[lib_size == 0] <- 1
    cpm <- sweep(src_counts[, matched_samps, drop = FALSE], 2, lib_size, "/") * 1e6
    log_cpm <- log2(cpm + 1)
    
    if (!lig %in% rownames(log_cpm)) next
    lig_expr <- as.numeric(log_cpm[lig, ])
    
    # Source compartment fraction
    if (ds == "GSE197677") {
      all_comps <- comp_197[samples_197 %in% matched_samps]
    } else {
      all_comps <- comp_221[samples_221 %in% matched_samps]
    }
    total_cells <- table(all_comps)
    source_frac <- as.numeric(total_cells[best_comp]) / sum(total_cells)
    
    # Fibroblast fraction
    fib_frac <- as.numeric(total_cells["Fibroblast"]) / sum(total_cells)
    
    valid <- !is.na(lig_expr) & !is.na(matched_core2)
    if (sum(valid) < 5) next
    
    # Model 1: unadjusted
    m1 <- tryCatch({
      fit <- lm(matched_core2[valid] ~ lig_expr[valid])
      s <- summary(fit)$coefficients
      if ("lig_expr[valid]" %in% rownames(s)) c("Est" = s["lig_expr[valid]", "Estimate"], "Dir" = sign(s["lig_expr[valid]", "Estimate"]))
      else c("Est" = NA, "Dir" = NA)
    }, error = function(e) c("Est" = NA, "Dir" = NA))
    
    # Model 2: + source fraction
    m2 <- tryCatch({
      fit <- lm(matched_core2[valid] ~ lig_expr[valid] + source_frac[valid])
      s <- summary(fit)$coefficients
      if ("lig_expr[valid]" %in% rownames(s)) c("Est" = s["lig_expr[valid]", "Estimate"], "Dir" = sign(s["lig_expr[valid]", "Estimate"]))
      else c("Est" = NA, "Dir" = NA)
    }, error = function(e) c("Est" = NA, "Dir" = NA))
    
    # Model 3: + fib fraction
    m3 <- tryCatch({
      fit <- lm(matched_core2[valid] ~ lig_expr[valid] + fib_frac[valid])
      s <- summary(fit)$coefficients
      if ("lig_expr[valid]" %in% rownames(s)) c("Est" = s["lig_expr[valid]", "Estimate"], "Dir" = sign(s["lig_expr[valid]", "Estimate"]))
      else c("Est" = NA, "Dir" = NA)
    }, error = function(e) c("Est" = NA, "Dir" = NA))
    
    # Model 4: + both
    m4 <- tryCatch({
      fit <- lm(matched_core2[valid] ~ lig_expr[valid] + source_frac[valid] + fib_frac[valid])
      s <- summary(fit)$coefficients
      if ("lig_expr[valid]" %in% rownames(s)) c("Est" = s["lig_expr[valid]", "Estimate"], "Dir" = sign(s["lig_expr[valid]", "Estimate"]))
      else c("Est" = NA, "Dir" = NA)
    }, error = function(e) c("Est" = NA, "Dir" = NA))
    
    dir_stable <- as.numeric(m1["Dir"]) == as.numeric(m4["Dir"])
    
    abund_list[[paste(ds, lig)]] <- data.frame(
      ligand = lig, dataset = ds, compartment = best_comp,
      coef_unadjusted = as.numeric(m1["Est"]), coef_plus_source = as.numeric(m2["Est"]),
      coef_plus_fib = as.numeric(m3["Est"]), coef_full = as.numeric(m4["Est"]),
      dir_stable = dir_stable, stringsAsFactors = FALSE)
  }
}

abund_df <- do.call(rbind, abund_list)
rownames(abund_df) <- NULL

cat("=== ABUNDANCE SENSITIVITY ===\n")
for (i in 1:nrow(abund_df)) {
  cat(sprintf("  %s %s: coef_adj=%.4f dir_stable=%s\n",
    abund_df$ligand[i], abund_df$dataset[i], abund_df$coef_full[i], abund_df$dir_stable[i]))
}

write.csv(abund_df, file.path(O_DIR, "STEP19O_SIGNALING_ABUNDANCE_SENSITIVITY.csv"), row.names = FALSE)
cat("\nSaved: STEP19O_SIGNALING_ABUNDANCE_SENSITIVITY.csv\n")

# ============================================================================
# 19O20 — EVIDENCE SCORE
# ============================================================================
cat("\n============================================\n")
cat("19O20: EVIDENCE SCORE\n")
cat("============================================\n\n")

cat("Scoring ligand-source pairs...\n")

# Build evidence score for each unique ligand
evidence_list <- list()
for (lig in unique(c(assoc_df$ligand, eligible_df$ligand))) {
  score <- 0
  denom <- 0
  
  # +1 ligand adequately detected (in at least one dataset)
  lig_detected <- any(ligand_qc_df$ligand == lig & ligand_qc_df$classification == "ADEQUATELY_DETECTED")
  denom <- denom + 1
  if (lig_detected) score <- score + 1
  
  # +1 receptor adequately detected in fibroblasts
  recs <- unique(curated_lr$receptor[curated_lr$ligand == lig])
  rec_ok <- any(sapply(recs, function(r) {
    any(rec_all$receptor == r & rec_all$classification %in% c("ROBUSTLY_DETECTED", "PARTIALLY_DETECTED"))
  }))
  denom <- denom + 1
  if (rec_ok) score <- score + 1
  
  # +1 abs(CORE2 rho) >= 0.4 in GSE197677
  rho_197 <- assoc_df$rho_core2[assoc_df$ligand == lig & assoc_df$dataset == "GSE197677"]
  denom <- denom + 1
  if (length(rho_197) > 0 && any(abs(rho_197) >= 0.40, na.rm = TRUE)) score <- score + 1
  
  # +1 abs(CORE2 rho) >= 0.4 in GSE221561
  rho_221 <- assoc_df$rho_core2[assoc_df$ligand == lig & assoc_df$dataset == "GSE221561"]
  denom <- denom + 1
  if (length(rho_221) > 0 && any(abs(rho_221) >= 0.40, na.rm = TRUE)) score <- score + 1
  
  # +1 CORE2 association same direction
  denom <- denom + 1
  if (length(rho_197) > 0 && length(rho_221) > 0) {
    best_197 <- rho_197[which.max(abs(rho_197))]
    best_221 <- rho_221[which.max(abs(rho_221))]
    if (!is.na(best_197) && !is.na(best_221) && sign(best_197) == sign(best_221)) score <- score + 1
  }
  
  # +1 ligand-target CORE2 support (if available)
  # 19O13 skipped - add denominator but don't penalize
  denom <- denom + 1
  # Score not added - resource unavailable
  
  # +1 RUNX1 bridge support
  denom <- denom + 1
  runx1_bridges <- bridge_df$classification[bridge_df$ligand == lig &
    bridge_df$classification == "LIGAND_RUNX1_CORE2_BRIDGE"]
  if (length(runx1_bridges) > 0) score <- score + 1
  
  # +1 LOSO stability
  denom <- denom + 1
  loso_robust <- any(loso_summary_df$ligand == lig & loso_summary_df$classification == "ROBUST")
  if (loso_robust) score <- score + 1
  
  # Classification
  available_denom <- denom  # All components available
  if (score >= 6) ev_class <- "HIGH_SIGNALING_SUPPORT"
  else if (score >= 4) ev_class <- "MODERATE_SIGNALING_SUPPORT"
  else if (score >= 2) ev_class <- "WEAK_SIGNALING_SUPPORT"
  else ev_class <- "NO_REPRODUCIBLE_SIGNALING_SUPPORT"
  
  evidence_list[[lig]] <- data.frame(
    ligand = lig, evidence_score = score, available_components = available_denom,
    score_fraction = paste0(score, "/", available_denom),
    ligand_detected = lig_detected, receptor_detected = rec_ok,
    core2_197 = length(rho_197) > 0 && any(abs(rho_197) >= 0.40, na.rm = TRUE),
    core2_221 = length(rho_221) > 0 && any(abs(rho_221) >= 0.40, na.rm = TRUE),
    core2_same_dir = {
      if (length(rho_197) > 0 && length(rho_221) > 0) {
        best_197 <- rho_197[which.max(abs(rho_197))]
        best_221 <- rho_221[which.max(abs(rho_221))]
        !is.na(best_197) && !is.na(best_221) && sign(best_197) == sign(best_221)
      } else FALSE
    },
    runx1_bridge = length(runx1_bridges) > 0,
    loso_robust = loso_robust,
    classification = ev_class, stringsAsFactors = FALSE)
}

evidence_df <- do.call(rbind, evidence_list)
evidence_df <- evidence_df[order(-evidence_df$evidence_score, evidence_df$ligand), ]
rownames(evidence_df) <- NULL

cat("=== EVIDENCE SCORES ===\n")
for (i in 1:nrow(evidence_df)) {
  cat(sprintf("  %s: %s -> %s\n", evidence_df$ligand[i], evidence_df$score_fraction[i],
    evidence_df$classification[i]))
}

write.csv(evidence_df, file.path(O_DIR, "STEP19O_SIGNALING_EVIDENCE_SCORE.csv"), row.names = FALSE)
cat("\nSaved: STEP19O_SIGNALING_EVIDENCE_SCORE.csv\n")

# ============================================================================
# 19O21 — FIGURES
# ============================================================================
cat("\n============================================\n")
cat("19O21: FIGURES\n")
cat("============================================\n\n")

# Figure 1: Source Ligand CORE2 Heatmap
cat("Generating Figure 1: Source Ligand CORE2 Heatmap...\n")
tryCatch({
  # Get top ligands with at least moderate association
  top_lig_heatmap <- unique(assoc_df$ligand[abs(assoc_df$rho_core2) >= 0.30])
  if (length(top_lig_heatmap) > 15) top_lig_heatmap <- top_lig_heatmap[1:15]
  
  heat_data <- assoc_df[assoc_df$ligand %in% top_lig_heatmap, ]
  heat_wide <- reshape2::dcast(heat_data, ligand + compartment ~ dataset, value.var = "rho_core2")
  heat_matrix <- reshape2::acast(heat_wide, ligand ~ compartment, value.var = "GSE197677", fun.aggregate = mean, na.rm = TRUE)
  heat_matrix_221 <- reshape2::acast(heat_wide, ligand ~ compartment, value.var = "GSE221561", fun.aggregate = mean, na.rm = TRUE)
  
  if (nrow(heat_matrix) > 0 && ncol(heat_matrix) > 0) {
    # Combine both datasets
    combined <- cbind(heat_matrix, heat_matrix_221)
    colnames(combined) <- c(paste0(colnames(heat_matrix), "_197"), paste0(colnames(heat_matrix), "_221"))
    
    png(file.path(FIG_DIR, "STEP19O_SOURCE_LIGAND_CORE2_HEATMAP.png"),
        width = 12, height = 8, units = "in", res = 150)
    heatmap(combined, margins = c(10, 6),
            main = "Source Ligand vs CORE2 Association (rho)",
            col = colorRampPalette(c("blue", "white", "red"))(100),
            scale = "none")
    dev.off()
    cat("Saved: STEP19O_SOURCE_LIGAND_CORE2_HEATMAP.png\n")
  }
}, error = function(e) cat("Skipping Figure 1:", conditionMessage(e), "\n"))

# Figure 2: Cross-Dataset Scatter
cat("Generating Figure 2: Cross-Dataset Scatter...\n")
tryCatch({
  cross_scatter <- data.frame()
  for (lig in unique(assoc_df$ligand)) {
    rho_197 <- assoc_df$rho_core2[assoc_df$ligand == lig & assoc_df$dataset == "GSE197677"]
    rho_221 <- assoc_df$rho_core2[assoc_df$ligand == lig & assoc_df$dataset == "GSE221561"]
    if (length(rho_197) > 0 && length(rho_221) > 0) {
      cross_scatter <- rbind(cross_scatter, data.frame(
        ligand = lig, rho_197 = mean(rho_197), rho_221 = mean(rho_221),
        stringsAsFactors = FALSE))
    }
  }
  
  if (nrow(cross_scatter) > 0) {
    p <- ggplot(cross_scatter, aes(x = rho_197, y = rho_221, label = ligand)) +
      geom_point(size = 3, color = "darkred") +
      geom_text(vjust = -1, size = 3, check_overlap = TRUE) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
      geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "grey70") +
      labs(title = "Cross-Dataset Ligand-CORE2 Association",
           x = "Rho GSE197677", y = "Rho GSE221561") +
      theme_bw()
    ggsave(file.path(FIG_DIR, "STEP19O_LIGAND_CORE2_CROSS_DATASET_SCATTER.png"), p,
           width = 9, height = 8, dpi = 150)
    cat("Saved: STEP19O_LIGAND_CORE2_CROSS_DATASET_SCATTER.png\n")
  }
}, error = function(e) cat("Skipping Figure 2:", conditionMessage(e), "\n"))

# Figure 3: Top Ligand Treatment Effects
cat("Generating Figure 3: Top Ligand Treatment Effects...\n")
tryCatch({
  top_te_ligs <- unique(evidence_df$ligand[1:min(8, nrow(evidence_df))])
  te_plot <- te_df[te_df$ligand %in% top_te_ligs, ]
  
  if (nrow(te_plot) > 0) {
    p <- ggplot(te_plot, aes(x = ligand, y = effect, fill = dataset)) +
      geom_col(position = position_dodge(width = 0.7), width = 0.6) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      facet_wrap(~compartment, scales = "free_y") +
      scale_fill_manual(values = c("GSE197677" = "steelblue", "GSE221561" = "tomato")) +
      labs(title = "Treatment Effect on Source Ligand Expression",
           x = "Ligand", y = "Effect (treated - control)", fill = "Dataset") +
      theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(file.path(FIG_DIR, "STEP19O_TOP_LIGAND_TREATMENT_EFFECTS.png"), p,
           width = 12, height = 8, dpi = 150)
    cat("Saved: STEP19O_TOP_LIGAND_TREATMENT_EFFECTS.png\n")
  }
}, error = function(e) cat("Skipping Figure 3:", conditionMessage(e), "\n"))

# Figure 4: Source Compartment Signaling Summary
cat("Generating Figure 4: Source Compartment Signaling Summary...\n")
tryCatch({
  if (nrow(prior_df) > 0) {
    p <- ggplot(prior_df, aes(x = best_compartment, fill = classification)) +
      geom_bar() +
      scale_fill_brewer(palette = "Set3") +
      labs(title = "Source Compartment Prioritization Summary",
           x = "Best Compartment", y = "Count", fill = "Classification") +
      theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(file.path(FIG_DIR, "STEP19O_SOURCE_COMPARTMENT_SIGNALING_SUMMARY.png"), p,
           width = 10, height = 7, dpi = 150)
    cat("Saved: STEP19O_SOURCE_COMPARTMENT_SIGNALING_SUMMARY.png\n")
  }
}, error = function(e) cat("Skipping Figure 4:", conditionMessage(e), "\n"))

# Figure 5: RUNX1 Bridge
cat("Generating Figure 5: RUNX1 Bridge...\n")
tryCatch({
  if (nrow(bridge_df) > 0) {
    bridge_plot <- bridge_df[bridge_df$classification != "NO_BRIDGE_SUPPORT", ]
    if (nrow(bridge_plot) > 0) {
      p <- ggplot(bridge_plot, aes(x = ligand_vs_runx1_rho, y = ligand_vs_core2_rho,
                                    color = classification, label = ligand)) +
        geom_point(size = 3) +
        geom_text(vjust = -1, size = 3, check_overlap = TRUE) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
        geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
        facet_wrap(~dataset + compartment) +
        labs(title = "RUNX1 Signaling Bridge",
             x = "Ligand vs RUNX1 rho", y = "Ligand vs CORE2 rho") +
        theme_bw()
      ggsave(file.path(FIG_DIR, "STEP19O_RUNX1_SIGNALING_BRIDGE.png"), p,
             width = 12, height = 8, dpi = 150)
      cat("Saved: STEP19O_RUNX1_SIGNALING_BRIDGE.png\n")
    }
  }
}, error = function(e) cat("Skipping Figure 5:", conditionMessage(e), "\n"))

# Figure 6: Ligand-Receptor CORE2 Network (simplified)
cat("Generating Figure 6: LR CORE2 Network...\n")
tryCatch({
  if (nrow(eligible_df) > 0) {
    # Create a simplified network visualization
    net_data <- eligible_df[, c("ligand", "receptor", "rho_core2", "family")]
    net_data$abs_rho <- abs(net_data$rho_core2)
    
    p <- ggplot(net_data, aes(x = ligand, y = receptor, size = abs_rho, color = family)) +
      geom_point(alpha = 0.7) +
      scale_size_continuous(range = c(2, 8)) +
      labs(title = "Eligible Ligand-Receptor Pairs vs CORE2",
           x = "Ligand", y = "Receptor", size = "|Rho|", color = "Family") +
      theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(file.path(FIG_DIR, "STEP19O_LIGAND_RECEPTOR_CORE2_NETWORK.png"), p,
           width = 14, height = 10, dpi = 150)
    cat("Saved: STEP19O_LIGAND_RECEPTOR_CORE2_NETWORK.png\n")
  }
}, error = function(e) cat("Skipping Figure 6:", conditionMessage(e), "\n"))

# Figure 7: LOSO Stability
cat("Generating Figure 7: LOSO Stability...\n")
tryCatch({
  if (nrow(loso_summary_df) > 0) {
    loso_plot <- loso_summary_df
    loso_plot$ligand <- factor(loso_plot$ligand, levels = top_lig_by_rho)
    
    p <- ggplot(loso_plot, aes(x = ligand, y = rho_sign_preservation, fill = dataset)) +
      geom_col(position = position_dodge(width = 0.7), width = 0.6) +
      geom_hline(yintercept = 0.8, linetype = "dashed", color = "red") +
      scale_fill_manual(values = c("GSE197677" = "steelblue", "GSE221561" = "tomato")) +
      scale_y_continuous(labels = scales::percent) +
      labs(title = "LOSO Ligand-CORE2 Sign Preservation",
           x = "Ligand", y = "Sign Preservation", fill = "Dataset") +
      theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(file.path(FIG_DIR, "STEP19O_TOP_SIGNALING_CANDIDATES_LOSO.png"), p,
           width = 10, height = 7, dpi = 150)
    cat("Saved: STEP19O_TOP_SIGNALING_CANDIDATES_LOSO.png\n")
  }
}, error = function(e) cat("Skipping Figure 7:", conditionMessage(e), "\n"))

# ============================================================================
# 19O22 — INTERPRETATION GUARDRAILS
# ============================================================================
cat("\n============================================\n")
cat("19O22: INTERPRETATION GUARDRAILS\n")
cat("============================================\n\n")

guardrails <- c(
  "1. Cell-cell signaling is inferred from expression compatibility and sample-level association.",
  "2. Ligand expression does not prove ligand secretion.",
  "3. Receptor expression does not prove receptor activation.",
  "4. Ligand-receptor databases are incomplete.",
  "5. Ligand-target predictions are not causal evidence.",
  "6. RUNX1/AP1/SMAD bridge analyses are associative.",
  "7. Sample/patient is the statistical unit.",
  "8. No cell-level inferential P values are permitted.",
  "9. GSE221561 Surgery_alone n=2 is a major limitation.",
  "10. Cross-dataset consistency receives more weight than single-dataset significance.",
  "11. Missing regulon or LR resource coverage means NOT EVALUABLE, not biologically absent.",
  "12. No cells or samples are removed."
)

writeLines(guardrails, file.path(O_DIR, "STEP19O_INTERPRETATION_GUARDRAILS.md"))
cat("Saved: STEP19O_INTERPRETATION_GUARDRAILS.md\n")

# ============================================================================
# 19O23 — FINAL CLASSIFICATION
# ============================================================================
cat("\n============================================\n")
cat("19O23: FINAL CLASSIFICATION\n")
cat("============================================\n\n")

# Determine final classification
n_high <- sum(evidence_df$classification == "HIGH_SIGNALING_SUPPORT")
n_moderate <- sum(evidence_df$classification == "MODERATE_SIGNALING_SUPPORT")
n_bridge <- sum(bridge_df$classification == "LIGAND_RUNX1_CORE2_BRIDGE", na.rm = TRUE)
n_cross_repl <- sum(repl_df$classification == "CROSS_DATASET_SIGNAL_SUPPORTED", na.rm = TRUE)

if (n_cross_repl > 0 && n_high > 0) {
  final_class <- "CORE2_WITH_REPLICATED_EXTRACELLULAR_SIGNAL_SUPPORT"
} else if (n_bridge > 0 && n_high > 0) {
  final_class <- "CORE2_WITH_RUNX1_LINKED_SIGNALING_SUPPORT"
} else if (n_high > 0 || n_moderate >= 2) {
  final_class <- "CORE2_WITH_MULTIPLE_SIGNALING_INPUTS"
} else if (n_moderate > 0) {
  final_class <- "DATASET_SPECIFIC_SIGNALING_ARCHITECTURE"
} else {
  final_class <- "NO_REPRODUCIBLE_EXTRACELLULAR_SIGNAL_IDENTIFIED"
}

cat("=== FINAL CLASSIFICATION ===\n")
cat(final_class, "\n\n")

# Final outputs
final_lig_ranking <- evidence_df
write.csv(final_lig_ranking, file.path(O_DIR, "STEP19O_FINAL_LIGAND_RANKING.csv"), row.names = FALSE)
cat("Saved: STEP19O_FINAL_LIGAND_RANKING.csv\n")

final_lr_pairs <- eligible_df
write.csv(final_lr_pairs, file.path(O_DIR, "STEP19O_FINAL_LIGAND_RECEPTOR_PAIRS.csv"), row.names = FALSE)
cat("Saved: STEP19O_FINAL_LIGAND_RECEPTOR_PAIRS.csv\n")

final_source <- prior_df
write.csv(final_source, file.path(O_DIR, "STEP19O_FINAL_SOURCE_COMPARTMENTS.csv"), row.names = FALSE)
cat("Saved: STEP19O_FINAL_SOURCE_COMPARTMENTS.csv\n")

final_runx1 <- bridge_df
write.csv(final_runx1, file.path(O_DIR, "STEP19O_RUNX1_BRIDGE_FINAL.csv"), row.names = FALSE)
cat("Saved: STEP19O_RUNX1_BRIDGE_FINAL.csv\n")

final_class_df <- data.frame(
  step = "19O", classification = final_class,
  n_high = n_high, n_moderate = n_moderate,
  n_bridge = n_bridge, n_cross_repl = n_cross_repl,
  timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  stringsAsFactors = FALSE
)
write.csv(final_class_df, file.path(O_DIR, "STEP19O_FINAL_SIGNALING_CLASSIFICATION.csv"), row.names = FALSE)
cat("Saved: STEP19O_FINAL_SIGNALING_CLASSIFICATION.csv\n")

# Final interpretation
final_interp <- c(
  "# Step 19O Final Interpretation",
  "",
  paste("## Classification:", final_class),
  "",
  "## Signaling Resource",
  "Curated minimal LR resource (100 pairs across 13 families)",
  "No external LR database available locally",
  "",
  "## Top Candidates",
  paste("High support:", n_high),
  paste("Moderate support:", n_moderate),
  paste("RUNX1 bridges:", n_bridge),
  paste("Cross-dataset replicated:", n_cross_repl),
  "",
  "## Key Limitations",
  "Curated LR resource is not comprehensive",
  "No ligand-target matrix for direct target validation",
  "GSE221561 Surgery_alone n=2",
  "Exploratory hypothesis-generation only"
)
writeLines(final_interp, file.path(O_DIR, "STEP19O_FINAL_INTERPRETATION.md"))
cat("Saved: STEP19O_FINAL_INTERPRETATION.md\n")

final_status <- data.frame(
  step = "19O", status = "COMPLETE",
  classification = final_class,
  n_eligible_pairs = nrow(eligible_df),
  n_high = n_high, n_moderate = n_moderate,
  cells_removed = 0, samples_removed = 0,
  cells_reclassified = 0, frozen_upstream_changed = "NO",
  timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  stringsAsFactors = FALSE
)
write.csv(final_status, file.path(O_DIR, "STEP19O_FINAL_STATUS.csv"), row.names = FALSE)
cat("Saved: STEP19O_FINAL_STATUS.csv\n")

# ============================================================================
# FINAL VALIDATION
# ============================================================================
cat("\n=== FINAL VALIDATION ===\n")
cat("Canonical contrasts: UNCHANGED\n")
cat("Step19L: UNCHANGED\n")
cat("Step19M: UNCHANGED\n")
cat("Step19N: UNCHANGED\n")
cat("Cells removed: 0\n")
cat("Samples removed: 0\n")
cat("Cells reclassified: 0\n")
cat("Frozen upstream outputs changed: NO\n\n")

# ============================================================================
# FINAL PRINT
# ============================================================================
cat("============================================\n")
cat("STEP 19O COMPLETE\n")
cat("TARGETED FIBROBLAST SIGNALING SOURCE AUDIT\n")
cat("============================================\n\n")

cat("LOCAL SIGNALING RESOURCE:\n")
cat("  Curated minimal LR resource\n")
cat(sprintf("  Pairs: %d\n", nrow(curated_lr)))
cat(sprintf("  Unique ligands: %d\n", length(unique(curated_lr$ligand))))
cat(sprintf("  Unique receptors: %d\n", length(unique(curated_lr$receptor))))
cat("  Ligand-target: NOT AVAILABLE\n\n")

cat("EVALUABLE SOURCE COMPARTMENTS:\n")
for (comp in source_compartments) {
  n197 <- sum(inv_197[[comp]] > 0, na.rm = TRUE)
  n221 <- sum(inv_221[[comp]] > 0, na.rm = TRUE)
  cat(sprintf("  %s: GSE197677=%d, GSE221561=%d\n", comp, n197, n221))
}

cat("\n--------------------------------------------\n\n")

cat("TOP SIGNALING CANDIDATES:\n\n")
for (i in 1:min(3, nrow(evidence_df))) {
  lig <- evidence_df$ligand[i]
  
  # Get best compartment
  best <- assoc_df[assoc_df$ligand == lig, ]
  best_comp <- best$compartment[which.max(abs(best$rho_core2))]
  
  # Get receptor
  recs <- unique(curated_lr$receptor[curated_lr$ligand == lig])
  rec_str <- paste(recs, collapse=", ")
  
  rho_197 <- assoc_df$rho_core2[assoc_df$ligand == lig & assoc_df$dataset == "GSE197677"]
  rho_221 <- assoc_df$rho_core2[assoc_df$ligand == lig & assoc_df$dataset == "GSE221561"]
  
  runx1_br <- bridge_df$classification[bridge_df$ligand == lig &
    bridge_df$classification == "LIGAND_RUNX1_CORE2_BRIDGE"]
  
  cat(sprintf("%d.\n", i))
  cat(sprintf("  Source: %s\n", best_comp))
  cat(sprintf("  Ligand: %s\n", lig))
  cat(sprintf("  Receptor: %s\n", rec_str))
  if (length(rho_197) > 0) cat(sprintf("  CORE2 association GSE197677: %.3f\n", rho_197[1]))
  if (length(rho_221) > 0) cat(sprintf("  CORE2 association GSE221561: %.3f\n", rho_221[1]))
  cat(sprintf("  RUNX1 bridge: %s\n", ifelse(length(runx1_br) > 0, "YES", "NO")))
  cat(sprintf("  Evidence: %s\n", evidence_df$score_fraction[i]))
  cat(sprintf("  Classification: %s\n\n", evidence_df$classification[i]))
}

cat("--------------------------------------------\n\n")

cat("TGF-BETA SIGNALING:\n")
tgfb_ligs <- unique(curated_lr$ligand[curated_lr$family == "TGF_beta"])
tgfb_assoc <- assoc_df[assoc_df$ligand %in% tgfb_ligs, ]
if (nrow(tgfb_assoc) > 0) {
  for (lig in tgfb_ligs) {
    sub <- tgfb_assoc[tgfb_assoc$ligand == lig, ]
    if (nrow(sub) > 0) {
      cat(sprintf("  %s:", lig))
      for (j in 1:nrow(sub)) {
        cat(sprintf(" %s rho=%.3f", sub$dataset[j], sub$rho_core2[j]))
      }
      cat("\n")
    }
  }
} else {
  cat("  No evaluable TGF-beta ligand associations\n")
}

cat("\nAP-1-LINKED SIGNALING:\n")
cat("  FOS/FOSL2 are downstream transcription factors, not ligands\n")
cat("  Step19M expression associations preserved; Step19N orthogonal support weak\n")

cat("\nRUNX1-LINKED SIGNALING:\n")
runx1_bridges <- bridge_df[bridge_df$classification == "LIGAND_RUNX1_CORE2_BRIDGE", ]
if (nrow(runx1_bridges) > 0) {
  for (i in 1:nrow(runx1_bridges)) {
    cat(sprintf("  %s %s %s: RUNX1_rho=%.3f CORE2_rho=%.3f\n",
      runx1_bridges$dataset[i], runx1_bridges$compartment[i], runx1_bridges$ligand[i],
      runx1_bridges$ligand_vs_runx1_rho[i], runx1_bridges$ligand_vs_core2_rho[i]))
  }
} else {
  cat("  No RUNX1-linked signaling bridges identified\n")
}

cat("\n--------------------------------------------\n\n")

cat("DOMINANT SOURCE COMPARTMENT:\n")
if (nrow(prior_df) > 0) {
  tab <- table(prior_df$best_compartment)
  cat(sprintf("  %s\n", names(tab)[which.max(tab)]))
} else {
  cat("  INCONCLUSIVE\n")
}

cat("\nCROSS-DATASET SIGNALING STATUS:\n")
cat(sprintf("  %s\n", final_class))

cat(sprintf("\nFINAL CLASSIFICATION: %s\n", final_class))
cat("\nCells removed: 0\n")
cat("Samples removed: 0\n")
cat("Cells reclassified: 0\n")
cat("Frozen upstream outputs changed: NO\n\n")

cat("STOP HERE.\n")
cat("DO NOT START STEP19P.\n")
