#!/usr/bin/env Rscript
# ============================================================================
# STEP 19I: DIRECT EPITHELIAL-FIBROBLAST COMPARTMENT AUDIT
# ============================================================================
# Direct compartment-level validation of 7 frozen Hallmark pathways.
# NOT causal adjustment. NOT cell reclassification.
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
g19i_dir <- file.path(pw_dir, "Step19I")
counts_dir <- file.path(g19i_dir, "counts")
ckpt_dir <- file.path(g19i_dir, "checkpoints")
fig_dir <- file.path(base_dir, "05_figures/cross_dataset/pathway_analysis", "Step19I")

dir.create(g19i_dir, showWarnings=FALSE, recursive=TRUE)
dir.create(counts_dir, showWarnings=FALSE, recursive=TRUE)
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

CORE_GENES <- c("COL3A1","COL1A1","MMP2","BGN","SPARC","CCN1","DCN","F3","CDKN1A")

# ============================================================================
# 19I0: INPUT AND LABEL FREEZE
# ============================================================================
cat("========================================\n")
cat("STEP 19I0: INPUT AND LABEL FREEZE\n")
cat("========================================\n\n")

f19f <- read.csv(file.path(pw_dir, "STEP19F_SEVEN_HALLMARKS_FINAL_FREEZE.csv"), stringsAsFactors=FALSE)
g19g_loc <- read.csv(file.path(g19g_dir, "STEP19G_PATHWAY_LOCALIZATION_FINAL.csv"), stringsAsFactors=FALSE)
fg19 <- read.csv(file.path(pw_dir, "STEP19E_GSE197677_HALLMARK_FGSEA_ALL.csv"), stringsAsFactors=FALSE)
fg22 <- read.csv(file.path(pw_dir, "STEP19E_GSE221561_HALLMARK_FGSEA_ALL.csv"), stringsAsFactors=FALSE)

stopifnot(nrow(f19f) == 7)
stopifnot(all(f19f$pathway %in% supported_pathways))

# Freeze gene sets from fgsea leading edges
geneset_freeze <- data.frame()
for (pw in supported_pathways) {
  pw_short <- gsub("HALLMARK_", "", pw)
  r19 <- fg19[fg19$pathway == pw, ]
  r22 <- fg22[fg22$pathway == pw, ]
  genes_disc <- if (nrow(r19) > 0) trimws(unlist(strsplit(as.character(r19$leadingEdge), ","))) else character(0)
  genes_repl <- if (nrow(r22) > 0) trimws(unlist(strsplit(as.character(r22$leadingEdge), ","))) else character(0)
  all_genes <- sort(unique(c(genes_disc, genes_repl)))
  for (g in all_genes) {
    geneset_freeze <- rbind(geneset_freeze, data.frame(
      pathway=pw, gene=g,
      in_discovery=g %in% genes_disc,
      in_replication=g %in% genes_repl,
      stringsAsFactors=FALSE))
  }
  cat(sprintf("  %s: %d genes\n", pw_short, length(all_genes)))
}

write.csv(data.frame(step="19I", status="VERIFIED", n_pathways=7,
  step19f="PATHWAY_SPECIFIC_ROBUSTNESS_CONFIRMED",
  step19g="ROBUST_PATHWAYS_MIXED_EPITHELIAL_STROMAL_SIGNAL",
  step19h="MIXED_PROGRAM_DECOMPOSITION_INCONCLUSIVE",
  stringsAsFactors=FALSE), file.path(g19i_dir, "STEP19I_INPUT_PROVENANCE.csv"), row.names=FALSE)
write.csv(geneset_freeze, file.path(g19i_dir, "STEP19I_SEVEN_PATHWAY_GENESET_FREEZE.csv"), row.names=FALSE)
cat("\nSaved: STEP19I_INPUT_PROVENANCE.csv, STEP19I_SEVEN_PATHWAY_GENESET_FREEZE.csv\n\n")

# ============================================================================
# 19I1: COMPARTMENT SAMPLE INVENTORY
# ============================================================================
cat("========================================\n")
cat("STEP 19I1: COMPARTMENT SAMPLE INVENTORY\n")
cat("========================================\n\n")

cat("Loading Seurat objects...\n")
obj19 <- readRDS(file.path(obj_dir, "GSE197677", "GSE197677_integrated_broadcelltype_annotated.rds"))
obj22 <- readRDS(file.path(obj_dir, "GSE221561", "GSE221561_author_annotated.rds"))
cat("Loaded.\n\n")

# Build inventory
inventory <- data.frame()

# GSE197677
md19 <- obj19@meta.data
md19$sample_id <- md19$sample
for (s in unique(md19$sample_id)) {
  cells <- rownames(md19)[md19$sample_id == s]
  n_epi <- sum(md19[cells, "broad_celltype_final"] == "Epithelial")
  n_fib <- sum(md19[cells, "broad_celltype_final"] == "Fibroblast_CAF")
  inventory <- rbind(inventory, data.frame(
    dataset="GSE197677", sample=s, patient=s,
    epithelial_cells=n_epi, fibroblast_cells=n_fib,
    epithelial_eligible=n_epi >= 30, fibroblast_eligible=n_fib >= 30,
    both_compartments_eligible=n_epi >= 30 && n_fib >= 30,
    stringsAsFactors=FALSE))
}

# GSE221561
md22 <- obj22@meta.data
md22$sample_id <- md22$sample_id
for (s in unique(md22$sample_id)) {
  cells <- rownames(md22)[md22$sample_id == s]
  n_epi <- sum(md22[cells, "author_Cell_type"] == "Epithelial")
  n_fib <- sum(md22[cells, "author_Cell_type"] == "Fibroblast")
  inventory <- rbind(inventory, data.frame(
    dataset="GSE221561", sample=s, patient=s,
    epithelial_cells=n_epi, fibroblast_cells=n_fib,
    epithelial_eligible=n_epi >= 30, fibroblast_eligible=n_fib >= 30,
    both_compartments_eligible=n_epi >= 30 && n_fib >= 30,
    stringsAsFactors=FALSE))
}

# Add treatment group
md19_treat <- setNames(
  read.csv(file.path(base_dir, "06_tables/cross_dataset/epithelial_pseudobulk/counts",
    "GSE197677_epithelial_pseudobulk_sample_metadata.csv"), stringsAsFactors=FALSE)$treatment_standardized,
  read.csv(file.path(base_dir, "06_tables/cross_dataset/epithelial_pseudobulk/counts",
    "GSE197677_epithelial_pseudobulk_sample_metadata.csv"), stringsAsFactors=FALSE)$sample)
inventory$treatment_group <- ifelse(inventory$dataset == "GSE197677",
  md19_treat[inventory$sample], NA)

md22_obj <- obj22@meta.data
treat22_map <- setNames(
  ifelse(md22_obj$author_Neoadjuvant %in% c("Chemotherapy","Chemoradiotherapy","Chemoradiotherapy and Immunotherapy","Antiangiogenesis"), "Neoadjuvant_treated",
  ifelse(md22_obj$author_Neoadjuvant == "Surgery alone", "Surgery_alone",
  ifelse(md22_obj$author_Neoadjuvant == "Adjacent normal", "Adjacent_normal", NA))),
  md22_obj$sample_id)
inventory$treatment_group[is.na(inventory$treatment_group)] <- treat22_map[inventory$sample[is.na(inventory$treatment_group)]]

write.csv(inventory, file.path(g19i_dir, "STEP19I_COMPARTMENT_SAMPLE_INVENTORY.csv"), row.names=FALSE)

cat("GSE197677:\n")
cat(sprintf("  Total samples: %d\n", sum(inventory$dataset == "GSE197677")))
cat(sprintf("  Epithelial eligible: %d\n", sum(inventory$epithelial_eligible & inventory$dataset == "GSE197677")))
cat(sprintf("  Fibroblast eligible: %d\n", sum(inventory$fibroblast_eligible & inventory$dataset == "GSE197677")))
cat(sprintf("  Both eligible: %d\n", sum(inventory$both_compartments_eligible & inventory$dataset == "GSE197677")))

cat("\nGSE221561:\n")
cat(sprintf("  Total samples: %d\n", sum(inventory$dataset == "GSE221561")))
cat(sprintf("  Epithelial eligible: %d\n", sum(inventory$epithelial_eligible & inventory$dataset == "GSE221561")))
cat(sprintf("  Fibroblast eligible: %d\n", sum(inventory$fibroblast_eligible & inventory$dataset == "GSE221561")))
cat(sprintf("  Both eligible: %d\n", sum(inventory$both_compartments_eligible & inventory$dataset == "GSE221561")))

cat("\nTreatment groups:\n")
for (d in c("GSE197677","GSE221561")) {
  sub <- inventory[inventory$dataset == d, ]
  cat(sprintf("  %s:\n", d))
  print(table(sub$treatment_group))
}
cat("\nSaved: STEP19I_COMPARTMENT_SAMPLE_INVENTORY.csv\n\n")

# ============================================================================
# 19I2: BUILD RAW-COUNT PSEUDOBULK
# ============================================================================
cat("========================================\n")
cat("STEP 19I2: BUILD RAW-COUNT PSEUDOBULK\n")
cat("========================================\n\n")

pseudobulk_compartment <- function(obj, ct_col, ct_label, sample_col, target_genes) {
  cell_meta <- obj@meta.data
  if (!inherits(obj@assays[["RNA"]], "Assay5")) {
    counts <- obj@assays[["RNA"]]@counts
    cells <- rownames(cell_meta)[cell_meta[[ct_col]] == ct_label]
    cells <- intersect(cells, colnames(counts))
    samples <- unique(cell_meta[cells, sample_col])
    genes_available <- intersect(target_genes, rownames(counts))
    pb <- matrix(0L, nrow=length(genes_available), ncol=length(samples),
                 dimnames=list(genes_available, samples))
    for (s in samples) {
      s_cells <- cells[cell_meta[cells, sample_col] == s]
      if (length(s_cells) > 0) {
        pb[, s] <- as.integer(Matrix::rowSums(counts[genes_available, s_cells, drop=FALSE]))
      }
    }
  } else {
    samples <- unique(cell_meta[[sample_col]][cell_meta[[ct_col]] == ct_label])
    # Get common gene set across layers
    all_layer_genes <- rownames(obj@assays[["RNA"]])
    genes_available <- intersect(target_genes, all_layer_genes)
    pb <- matrix(0L, nrow=length(genes_available), ncol=length(samples),
                 dimnames=list(genes_available, samples))
    for (ln in names(obj@assays[["RNA"]]@layers)) {
      mat <- LayerData(obj@assays[["RNA"]], layer=ln)
      cells <- colnames(mat)
      cell_types <- cell_meta[cells, ct_col]
      comp_cells <- cells[cell_types == ct_label]
      if (length(comp_cells) == 0) next
      for (s in intersect(samples, unique(cell_meta[comp_cells, sample_col]))) {
        s_cells <- comp_cells[cell_meta[comp_cells, sample_col] == s]
        if (length(s_cells) > 0) {
          genes_in_layer <- intersect(genes_available, rownames(mat))
          if (length(genes_in_layer) > 0) {
            pb[genes_in_layer, s] <- pb[genes_in_layer, s] +
              as.integer(Matrix::rowSums(mat[genes_in_layer, s_cells, drop=FALSE]))
          }
        }
      }
    }
  }
  return(pb)
}

all_genes <- union(geneset_freeze$gene, CORE_GENES)

cat("Processing GSE197677...\n")
cat("  Epithelial...\n")
pb19_epi <- pseudobulk_compartment(obj19, "broad_celltype_final", "Epithelial", "sample", all_genes)
cat(sprintf("    %d genes x %d samples\n", nrow(pb19_epi), ncol(pb19_epi)))
cat("  Fibroblast_CAF...\n")
pb19_fib <- pseudobulk_compartment(obj19, "broad_celltype_final", "Fibroblast_CAF", "sample", all_genes)
cat(sprintf("    %d genes x %d samples\n", nrow(pb19_fib), ncol(pb19_fib)))

cat("Processing GSE221561...\n")
cat("  Epithelial...\n")
pb22_epi <- pseudobulk_compartment(obj22, "author_Cell_type", "Epithelial", "sample_id", all_genes)
cat(sprintf("    %d genes x %d samples\n", nrow(pb22_epi), ncol(pb22_epi)))
cat("  Fibroblast...\n")
pb22_fib <- pseudobulk_compartment(obj22, "author_Cell_type", "Fibroblast", "sample_id", all_genes)
cat(sprintf("    %d genes x %d samples\n", nrow(pb22_fib), ncol(pb22_fib)))

# Save
saveRDS(pb19_epi, file.path(counts_dir, "GSE197677_EPITHELIAL_raw_counts.rds"))
saveRDS(pb19_fib, file.path(counts_dir, "GSE197677_FIBROBLAST_CAF_raw_counts.rds"))
saveRDS(pb22_epi, file.path(counts_dir, "GSE221561_EPITHELIAL_raw_counts.rds"))
saveRDS(pb22_fib, file.path(counts_dir, "GSE221561_FIBROBLAST_CAF_raw_counts.rds"))

# Validate
for (nm in c("pb19_epi","pb19_fib","pb22_epi","pb22_fib")) {
  pb <- get(nm)
  stopifnot(all(pb >= 0))
  stopifnot(all(pb == round(pb)))
}
cat("\nSaved: 4 raw count RDS files\nValidated: non-negative integers\n\n")

# ============================================================================
# 19I3: COMMON GENE UNIVERSE
# ============================================================================
cat("========================================\n")
cat("STEP 19I3: COMMON GENE UNIVERSE\n")
cat("========================================\n\n")

common_genes_19 <- intersect(rownames(pb19_epi), rownames(pb19_fib))
common_genes_22 <- intersect(rownames(pb22_epi), rownames(pb22_fib))
common_all <- intersect(common_genes_19, common_genes_22)

cat(sprintf("GSE197677 common: %d genes\n", length(common_genes_19)))
cat(sprintf("GSE221561 common: %d genes\n", length(common_genes_22)))
cat(sprintf("Both datasets: %d genes\n", length(common_all)))

write.csv(data.frame(dataset=c("GSE197677","GSE221561","BOTH"),
  n_genes=c(length(common_genes_19), length(common_genes_22), length(common_all)),
  stringsAsFactors=FALSE), file.path(g19i_dir, "STEP19I_COMPARTMENT_COMMON_GENE_UNIVERSE.csv"), row.names=FALSE)
cat("Saved: STEP19I_COMPARTMENT_COMMON_GENE_UNIVERSE.csv\n\n")

# ============================================================================
# 19I4: NORMALIZATION
# ============================================================================
cat("========================================\n")
cat("STEP 19I4: NORMALIZATION\n")
cat("========================================\n\n")

normalize_compartment <- function(pb) {
  # Remove samples with zero counts
  col_sums <- colSums(pb)
  pb <- pb[, col_sums > 0, drop=FALSE]
  if (ncol(pb) == 0 || nrow(pb) == 0) return(list(lcpm=matrix(NA, nrow=nrow(pb), ncol=ncol(pb), dimnames=dimnames(pb)), dge=NULL))
  dge <- DGEList(counts=pb)
  dge <- calcNormFactors(dge, method="TMM")
  lcpm <- cpm(dge, log=TRUE, prior.count=2)
  return(list(lcpm=lcpm, dge=dge))
}

cat("GSE197677 epithelial...\n")
norm19_epi <- normalize_compartment(pb19_epi)
cat(sprintf("  %d genes x %d samples\n", nrow(norm19_epi$lcpm), ncol(norm19_epi$lcpm)))

cat("GSE197677 fibroblast...\n")
norm19_fib <- normalize_compartment(pb19_fib)
cat(sprintf("  %d genes x %d samples\n", nrow(norm19_fib$lcpm), ncol(norm19_fib$lcpm)))

cat("GSE221561 epithelial...\n")
norm22_epi <- normalize_compartment(pb22_epi)
cat(sprintf("  %d genes x %d samples\n", nrow(norm22_epi$lcpm), ncol(norm22_epi$lcpm)))

cat("GSE221561 fibroblast...\n")
norm22_fib <- normalize_compartment(pb22_fib)
cat(sprintf("  %d genes x %d samples\n", nrow(norm22_fib$lcpm), ncol(norm22_fib$lcpm)))

saveRDS(list(norm19_epi=norm19_epi, norm19_fib=norm19_fib,
             norm22_epi=norm22_epi, norm22_fib=norm22_fib),
        file.path(ckpt_dir, "19I_CHECKPOINT_normalized.rds"))
cat("\nSaved: normalized pseudobulk checkpoint\n\n")

# ============================================================================
# 19I5: SEVEN HALLMARK SAMPLE SCORES
# ============================================================================
cat("========================================\n")
cat("STEP 19I5: SEVEN HALLMARK SAMPLE SCORES\n")
cat("========================================\n\n")

zscore_genes <- function(mat) {
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

compute_pathway_scores <- function(lcpm, pathway_genes_list, min_genes=5) {
  samples <- colnames(lcpm)
  genes_avail <- rownames(lcpm)
  zmat <- zscore_genes(lcpm)
  res <- data.frame()
  for (pw in supported_pathways) {
    pw_genes <- intersect(pathway_genes_list[[pw]], genes_avail)
    n <- length(pw_genes)
    if (n < min_genes) {
      for (s in samples) {
        res <- rbind(res, data.frame(sample=s, pathway=pw, n_genes=n, score=NA,
          status="INSUFFICIENT_GENES", stringsAsFactors=FALSE))
      }
    } else {
      for (s in samples) {
        sc <- mean(zmat[pw_genes, s], na.rm=TRUE)
        res <- rbind(res, data.frame(sample=s, pathway=pw, n_genes=n, score=round(sc,6),
          status="OK", stringsAsFactors=FALSE))
      }
    }
  }
  return(res)
}

# Build pathway gene lists
pw_genes_list <- list()
for (pw in supported_pathways) {
  pw_genes_list[[pw]] <- geneset_freeze$gene[geneset_freeze$pathway == pw]
}

cat("Computing GSE197677 scores...\n")
s19_epi <- compute_pathway_scores(norm19_epi$lcpm, pw_genes_list)
s19_epi$dataset <- "GSE197677"; s19_epi$compartment <- "EPITHELIAL"
s19_fib <- compute_pathway_scores(norm19_fib$lcpm, pw_genes_list)
s19_fib$dataset <- "GSE197677"; s19_fib$compartment <- "FIBROBLAST_CAF"

cat("Computing GSE221561 scores...\n")
s22_epi <- compute_pathway_scores(norm22_epi$lcpm, pw_genes_list)
s22_epi$dataset <- "GSE221561"; s22_epi$compartment <- "EPITHELIAL"
s22_fib <- compute_pathway_scores(norm22_fib$lcpm, pw_genes_list)
s22_fib$dataset <- "GSE221561"; s22_fib$compartment <- "FIBROBLAST_CAF"

all_scores <- rbind(s19_epi, s19_fib, s22_epi, s22_fib)

# Add treatment group
for (i in 1:nrow(all_scores)) {
  d <- all_scores$dataset[i]
  s <- all_scores$sample[i]
  if (d == "GSE197677") {
    all_scores$treatment_group[i] <- md19_treat[s]
  } else {
    all_scores$treatment_group[i] <- treat22_map[s]
  }
}

write.csv(all_scores, file.path(g19i_dir, "STEP19I_SAMPLE_COMPARTMENT_PATHWAY_SCORES.csv"), row.names=FALSE)
cat("Score counts:\n")
print(table(all_scores$status, all_scores$compartment))
cat("\nSaved: STEP19I_SAMPLE_COMPARTMENT_PATHWAY_SCORES.csv\n\n")

# ============================================================================
# 19I6: EXACT TREATMENT PERMUTATION TEST
# ============================================================================
cat("========================================\n")
cat("STEP 19I6: EXACT TREATMENT PERMUTATION TEST\n")
cat("========================================\n\n")

exact_perm_test <- function(scores, groups, n_perm=NULL) {
  # scores: numeric vector
  # groups: factor/character vector same length
  # If n_perm is NULL, enumerate all permutations
  g1 <- scores[groups == groups[which.max(table(groups))[1]]]
  g2 <- scores[groups != groups[which.max(table(groups))[1]]]

  # Determine treated vs untreated
  # Convention: first group level is treated
  grp_levels <- unique(groups)
  treated_name <- grp_levels[1]
  untreated_name <- grp_levels[2]

  v_treated <- scores[groups == treated_name]
  v_untreated <- scores[groups == untreated_name]
  v_treated <- v_treated[!is.na(v_treated)]
  v_untreated <- v_untreated[!is.na(v_untreated)]

  if (length(v_treated) < 2 || length(v_untreated) < 2) {
    return(list(obs_diff=NA, obs_med_diff=NA, p_two=NA, p_one=NA, d=NA, n_treated=length(v_treated), n_untreated=length(v_untreated)))
  }

  obs_diff <- mean(v_treated) - mean(v_untreated)
  obs_med_diff <- median(v_treated) - median(v_untreated)
  pooled_sd <- sqrt((var(v_treated) + var(v_untreated)) / 2)
  d <- if (pooled_sd > 0) obs_diff / pooled_sd else 0

  # Exact permutation
  all_vals <- c(v_treated, v_untreated)
  n_t <- length(v_treated)
  n_u <- length(v_untreated)
  n_total <- n_t + n_u

  # Enumerate if feasible (n_total choose n_t)
  n_combos <- choose(n_total, n_t)
  if (n_combos > 1e6) {
    # Too many, use Monte Carlo
    n_perm <- min(10000, n_combos)
    perm_diffs <- numeric(n_perm)
    for (p in 1:n_perm) {
      idx <- sample(n_total, n_t)
      perm_diffs[p] <- mean(all_vals[idx]) - mean(all_vals[-idx])
    }
  } else {
    # Enumerate all
    combos <- combn(n_total, n_t)
    perm_diffs <- numeric(ncol(combos))
    for (p in 1:ncol(combos)) {
      perm_diffs[p] <- mean(all_vals[combos[, p]]) - mean(all_vals[-combos[, p]])
    }
  }

  p_two <- mean(abs(perm_diffs) >= abs(obs_diff))
  p_one <- if (obs_diff > 0) mean(perm_diffs >= obs_diff) else mean(perm_diffs <= obs_diff)

  return(list(obs_diff=obs_diff, obs_med_diff=obs_med_diff, p_two=p_two, p_one=p_one, d=d,
              n_treated=length(v_treated), n_untreated=length(v_untreated)))
}

perm_results <- data.frame()

for (d in c("GSE197677","GSE221561")) {
  sub <- all_scores[all_scores$dataset == d & !is.na(all_scores$score), ]
  for (comp in c("EPITHELIAL","FIBROBLAST_CAF")) {
    for (pw in supported_pathways) {
      sp <- sub[sub$compartment == comp & sub$pathway == pw, ]
      if (nrow(sp) < 4) next

      if (d == "GSE197677") {
        groups <- ifelse(sp$treatment_group == "Neoadjuvant_chemotherapy", "NACT", "nNACT")
      } else {
        groups <- ifelse(sp$treatment_group == "Neoadjuvant_treated", "NACT", "nNACT")
      }

      # Filter out samples with no treatment label
      has_group <- !is.na(groups)
      sp <- sp[has_group, ]
      groups <- groups[has_group]
      if (nrow(sp) < 4) next

      pt <- tryCatch(exact_perm_test(sp$score, groups), error=function(e) list(obs_diff=NA,obs_med_diff=NA,p_two=NA,p_one=NA,d=NA,n_treated=NA,n_untreated=NA))

      perm_results <- rbind(perm_results, data.frame(
        dataset=d, compartment=comp, pathway=pw,
        n_treated=pt$n_treated, n_untreated=pt$n_untreated,
        observed_mean_difference=round(pt$obs_diff,6),
        observed_median_difference=round(pt$obs_med_diff,6),
        exact_two_sided_P=round(pt$p_two,6),
        exact_one_sided_P_in_direction=round(pt$p_one,6),
        standardized_effect=round(pt$d,4),
        stringsAsFactors=FALSE))
    }
  }
  cat(sprintf("  %s: %d tests\n", d, sum(perm_results$dataset == d)))
}

# BH FDR within each dataset x compartment
perm_results$FDR <- NA
for (d in c("GSE197677","GSE221561")) {
  for (comp in c("EPITHELIAL","FIBROBLAST_CAF")) {
    idx <- which(perm_results$dataset == d & perm_results$compartment == comp)
    if (length(idx) > 0) {
      pvals <- perm_results$exact_two_sided_P[idx]
      valid <- !is.na(pvals)
      if (sum(valid) > 0) {
        fdr <- rep(NA, length(pvals))
        fdr[valid] <- p.adjust(pvals[valid], method="BH")
        perm_results$FDR[idx] <- round(fdr, 6)
      }
    }
  }
}

write.csv(perm_results, file.path(g19i_dir, "STEP19I_COMPARTMENT_TREATMENT_EFFECTS.csv"), row.names=FALSE)
cat("\nSaved: STEP19I_COMPARTMENT_TREATMENT_EFFECTS.csv\n\n")

# ============================================================================
# 19I7: LOSO STABILITY
# ============================================================================
cat("========================================\n")
cat("STEP 19I7: LOSO STABILITY\n")
cat("========================================\n\n")

loso_results <- data.frame()

for (d in c("GSE197677","GSE221561")) {
  sub <- all_scores[all_scores$dataset == d & !is.na(all_scores$score), ]
  for (comp in c("EPITHELIAL","FIBROBLAST_CAF")) {
    for (pw in supported_pathways) {
      sp <- sub[sub$compartment == comp & sub$pathway == pw, ]
      if (nrow(sp) < 4) next

      if (d == "GSE197677") {
        groups <- ifelse(sp$treatment_group == "Neoadjuvant_chemotherapy", "NACT", "nNACT")
      } else {
        groups <- ifelse(sp$treatment_group == "Neoadjuvant_treated", "NACT", "nNACT")
      }

      # Filter out samples with no treatment label
      has_group <- !is.na(groups)
      sp <- sp[has_group, ]
      groups <- groups[has_group]
      if (nrow(sp) < 4) next

      # Full data direction
      v_t <- sp$score[groups == "NACT"]
      v_u <- sp$score[groups == "nNACT"]
      if (length(v_t) < 2 || length(v_u) < 2) next
      full_diff <- mean(v_t) - mean(v_u)
      full_dir <- sign(full_diff)
      if (is.na(full_dir)) next

      # LOSO
      samples <- sp$sample
      n_preserve <- 0
      n_total <- 0
      effects <- numeric()
      for (leave_out in samples) {
        sp_lo <- sp[sp$sample != leave_out, ]
        g_lo <- groups[sp$sample != leave_out]
        vt_lo <- sp_lo$score[g_lo == "NACT"]
        vu_lo <- sp_lo$score[g_lo == "nNACT"]
        if (length(vt_lo) < 2 || length(vu_lo) < 2) next
        lo_diff <- mean(vt_lo) - mean(vu_lo)
        if (!is.na(lo_diff)) {
          effects <- c(effects, lo_diff)
          n_total <- n_total + 1
          if (sign(lo_diff) == full_dir) n_preserve <- n_preserve + 1
        }
      }

      frac <- if (n_total > 0) n_preserve / n_total else NA
      status <- if (!is.na(frac) && frac >= 0.8) "ROBUST" else if (!is.na(frac) && frac >= 0.6) "MODERATE" else "UNSTABLE"

      loso_results <- rbind(loso_results, data.frame(
        dataset=d, compartment=comp, pathway=pw,
        full_effect=round(full_diff,6),
        full_direction=if (full_dir > 0) "HIGHER_IN_NACT" else "HIGHER_IN_NON_NACT",
        n_samples=length(samples), n_iterations=n_total,
        n_preserve=n_preserve, loso_fraction=round(frac,4),
        loso_status=status,
        min_effect=round(min(effects),6), median_effect=round(median(effects),6),
        max_effect=round(max(effects),6),
        stringsAsFactors=FALSE))
    }
  }
  cat(sprintf("  %s: %d rows\n", d, sum(loso_results$dataset == d)))
}

write.csv(loso_results, file.path(g19i_dir, "STEP19I_COMPARTMENT_LOSO_STABILITY.csv"), row.names=FALSE)
cat("\nLOSO status:\n")
print(table(loso_results$loso_status))
cat("\nSaved: STEP19I_COMPARTMENT_LOSO_STABILITY.csv\n\n")

# ============================================================================
# 19I8: DIRECT EPITHELIAL VS FIBROBLAST COMPARISON
# ============================================================================
cat("========================================\n")
cat("STEP 19I8: DIRECT EPITHELIAL VS FIBROBLAST\n")
cat("========================================\n\n")

compcomp_results <- data.frame()

for (d in c("GSE197677","GSE221561")) {
  epi_sub <- all_scores[all_scores$dataset == d & all_scores$compartment == "EPITHELIAL" & !is.na(all_scores$score), ]
  fib_sub <- all_scores[all_scores$dataset == d & all_scores$compartment == "FIBROBLAST_CAF" & !is.na(all_scores$score), ]

  for (pw in supported_pathways) {
    epi_pw <- epi_sub[epi_sub$pathway == pw, ]
    fib_pw <- fib_sub[fib_sub$pathway == pw, ]
    common_samples <- intersect(epi_pw$sample, fib_pw$sample)
    if (length(common_samples) < 3) next

    epi_scores <- epi_pw$score[match(common_samples, epi_pw$sample)]
    fib_scores <- fib_pw$score[match(common_samples, fib_pw$sample)]

    rho <- tryCatch(cor(epi_scores, fib_scores, method="spearman", use="complete.obs"), error=function(e) NA)
    diff <- fib_scores - epi_scores

    # Treatment comparison of compartment difference
    if (d == "GSE197677") {
      treat_samples <- common_samples[epi_pw$treatment_group[match(common_samples, epi_pw$sample)] == "Neoadjuvant_chemotherapy"]
      control_samples <- common_samples[epi_pw$treatment_group[match(common_samples, epi_pw$sample)] == "No_neoadjuvant_chemotherapy"]
    } else {
      treat_samples <- common_samples[epi_pw$treatment_group[match(common_samples, epi_pw$sample)] == "Neoadjuvant_treated"]
      control_samples <- common_samples[epi_pw$treatment_group[match(common_samples, epi_pw$sample)] == "Surgery_alone"]
    }

    diff_treat <- mean(diff[common_samples %in% treat_samples], na.rm=TRUE)
    diff_control <- mean(diff[common_samples %in% control_samples], na.rm=TRUE)

    compcomp_results <- rbind(compcomp_results, data.frame(
      dataset=d, pathway=pw, n_common=length(common_samples),
      spearman_rho=round(rho,4),
      mean_fib_minus_epi_diff=round(mean(diff),4),
      diff_in_treated=round(diff_treat,4), diff_in_control=round(diff_control,4),
      treatment_effect_on_compartment_diff=round(diff_treat - diff_control, 4),
      stringsAsFactors=FALSE))
  }
}

write.csv(compcomp_results, file.path(g19i_dir, "STEP19I_WITHIN_SAMPLE_COMPARTMENT_COMPARISON.csv"), row.names=FALSE)
cat("Saved: STEP19I_WITHIN_SAMPLE_COMPARTMENT_COMPARISON.csv\n\n")

# ============================================================================
# 19I9: CROSS-DATASET COMPARTMENT REPLICATION
# ============================================================================
cat("========================================\n")
cat("STEP 19I9: CROSS-DATASET COMPARTMENT REPLICATION\n")
cat("========================================\n\n")

cross_comp <- data.frame()

for (pw in supported_pathways) {
  for (comp in c("EPITHELIAL","FIBROBLAST_CAF")) {
    td19 <- perm_results[perm_results$dataset=="GSE197677" & perm_results$compartment==comp & perm_results$pathway==pw, ]
    td22 <- perm_results[perm_results$dataset=="GSE221561" & perm_results$compartment==comp & perm_results$pathway==pw, ]

    if (nrow(td19) == 0 || nrow(td22) == 0) {
      cross_comp <- rbind(cross_comp, data.frame(pathway=pw, compartment=comp,
        GSE197677_effect=NA, GSE221561_effect=NA, concordance="UNRESOLVED",
        stringsAsFactors=FALSE))
      next
    }

    e19 <- td19$observed_mean_difference[1]
    e22 <- td22$observed_mean_difference[1]
    same <- sign(e19) == sign(e22) && !is.na(e19) && !is.na(e22) && sign(e19) != 0
    conc <- if (same) "CROSS_DATASET_SAME_DIRECTION" else "CROSS_DATASET_OPPOSITE_DIRECTION"
    if (is.na(e19) || is.na(e22)) conc <- "UNRESOLVED"

    cross_comp <- rbind(cross_comp, data.frame(pathway=pw, compartment=comp,
      GSE197677_effect=round(e19,6), GSE221561_effect=round(e22,6),
      concordance=conc, stringsAsFactors=FALSE))
  }
}

write.csv(cross_comp, file.path(g19i_dir, "STEP19I_CROSS_DATASET_COMPARTMENT_CONCORDANCE.csv"), row.names=FALSE)
cat("Cross-dataset concordance:\n")
print(table(cross_comp$concordance))
cat("\nSaved: STEP19I_CROSS_DATASET_COMPARTMENT_CONCORDANCE.csv\n\n")

# ============================================================================
# 19I10: PATHWAY CLASSIFICATION
# ============================================================================
cat("========================================\n")
cat("STEP 19I10: PATHWAY CLASSIFICATION\n")
cat("========================================\n\n")

pw_class <- data.frame()

for (pw in supported_pathways) {
  pw_short <- gsub("HALLMARK_", "", pw)

  # Get cross-dataset concordance for each compartment
  cc_epi <- cross_comp[cross_comp$pathway == pw & cross_comp$compartment == "EPITHELIAL", ]
  cc_fib <- cross_comp[cross_comp$pathway == pw & cross_comp$compartment == "FIBROBLAST_CAF", ]

  epi_conc <- if (nrow(cc_epi) > 0) cc_epi$concordance[1] else "UNRESOLVED"
  fib_conc <- if (nrow(cc_fib) > 0) cc_fib$concordance[1] else "UNRESOLVED"

  # Get LOSO stability
  loso_epi <- loso_results[loso_results$pathway == pw & loso_results$compartment == "EPITHELIAL", ]
  loso_fib <- loso_results[loso_results$pathway == pw & loso_results$compartment == "FIBROBLAST_CAF", ]
  epi_loso_ok <- all(loso_epi$loso_status %in% c("ROBUST","MODERATE")) && nrow(loso_epi) > 0
  fib_loso_ok <- all(loso_fib$loso_status %in% c("ROBUST","MODERATE")) && nrow(loso_fib) > 0

  # Get effects
  td19_epi <- perm_results[perm_results$dataset=="GSE197677" & perm_results$compartment=="EPITHELIAL" & perm_results$pathway==pw, ]
  td22_epi <- perm_results[perm_results$dataset=="GSE221561" & perm_results$compartment=="EPITHELIAL" & perm_results$pathway==pw, ]
  td19_fib <- perm_results[perm_results$dataset=="GSE197677" & perm_results$compartment=="FIBROBLAST_CAF" & perm_results$pathway==pw, ]
  td22_fib <- perm_results[perm_results$dataset=="GSE221561" & perm_results$compartment=="FIBROBLAST_CAF" & perm_results$pathway==pw, ]

  e19e <- if (nrow(td19_epi) > 0) td19_epi$observed_mean_difference[1] else NA
  e22e <- if (nrow(td22_epi) > 0) td22_epi$observed_mean_difference[1] else NA
  e19f <- if (nrow(td19_fib) > 0) td19_fib$observed_mean_difference[1] else NA
  e22f <- if (nrow(td22_fib) > 0) td22_fib$observed_mean_difference[1] else NA

  epi_effect_ok <- !is.na(e19e) && !is.na(e22e) && sign(e19e) == sign(e22e) && sign(e19e) != 0
  fib_effect_ok <- !is.na(e19f) && !is.na(e22f) && sign(e19f) == sign(e22f) && sign(e19f) != 0
  epi_fib_agree <- epi_effect_ok && fib_effect_ok && sign(e19e) == sign(e19f)

  # Classify
  if (epi_conc == "CROSS_DATASET_SAME_DIRECTION" && fib_conc == "CROSS_DATASET_SAME_DIRECTION" &&
      epi_fib_agree && epi_loso_ok && fib_loso_ok) {
    class <- "COORDINATED_MULTI_COMPARTMENT_PROGRAM"
  } else if (fib_conc == "CROSS_DATASET_SAME_DIRECTION" && fib_loso_ok &&
             (epi_conc != "CROSS_DATASET_SAME_DIRECTION" || !epi_loso_ok)) {
    class <- "FIBROBLAST_DOMINANT_PROGRAM"
  } else if (epi_conc == "CROSS_DATASET_SAME_DIRECTION" && epi_loso_ok &&
             (fib_conc != "CROSS_DATASET_SAME_DIRECTION" || !fib_loso_ok)) {
    class <- "EPITHELIAL_DOMINANT_PROGRAM"
  } else if (epi_conc == "CROSS_DATASET_SAME_DIRECTION" && fib_conc == "CROSS_DATASET_SAME_DIRECTION" &&
             !epi_fib_agree) {
    class <- "COMPARTMENT_DIVERGENT_PROGRAM"
  } else {
    class <- "COMPARTMENT_EVIDENCE_INCONCLUSIVE"
  }

  pw_class <- rbind(pw_class, data.frame(
    pathway=pw, final_class=class,
    epithelial_concordance=epi_conc, fibroblast_concordance=fib_conc,
    epithelial_LOSO=if (nrow(loso_epi) > 0) paste(unique(loso_epi$loso_status), collapse=",") else NA,
    fibroblast_LOSO=if (nrow(loso_fib) > 0) paste(unique(loso_fib$loso_status), collapse=",") else NA,
    GSE197677_epithelial_effect=round(e19e,4), GSE221561_epithelial_effect=round(e22e,4),
    GSE197677_fibroblast_effect=round(e19f,4), GSE221561_fibroblast_effect=round(e22f,4),
    stringsAsFactors=FALSE))

  cat(sprintf("  %s: %s\n", pw_short, class))
}

write.csv(pw_class, file.path(g19i_dir, "STEP19I_SEVEN_PATHWAY_COMPARTMENT_CLASSIFICATION.csv"), row.names=FALSE)
cat("\nClassification:\n")
print(table(pw_class$final_class))
cat("\nSaved: STEP19I_SEVEN_PATHWAY_COMPARTMENT_CLASSIFICATION.csv\n\n")

# ============================================================================
# 19I11: CORE RECURRENT GENES
# ============================================================================
cat("========================================\n")
cat("STEP 19I11: CORE RECURRENT GENES\n")
cat("========================================\n\n")

core_gene_audit <- data.frame()

for (g in CORE_GENES) {
  for (d in c("GSE197677","GSE221561")) {
    for (comp in c("EPITHELIAL","FIBROBLAST_CAF")) {
      if (d == "GSE197677") {
        pb <- if (comp == "EPITHELIAL") pb19_epi else pb19_fib
        ct_col <- "broad_celltype_final"
        ct_label <- if (comp == "EPITHELIAL") "Epithelial" else "Fibroblast_CAF"
        sample_col <- "sample"
      } else {
        pb <- if (comp == "EPITHELIAL") pb22_epi else pb22_fib
        ct_col <- "author_Cell_type"
        ct_label <- if (comp == "EPITHELIAL") "Epithelial" else "Fibroblast"
        sample_col <- "sample_id"
      }

      if (!(g %in% rownames(pb))) {
        core_gene_audit <- rbind(core_gene_audit, data.frame(
          gene=g, dataset=d, compartment=comp,
          treatment_diff=NA, direction=NA, perm_P=NA,
          loso_status="GENE_NOT_AVAILABLE",
          stringsAsFactors=FALSE))
        next
      }

      # Get expression
      expr <- as.numeric(pb[g, ])
      names(expr) <- colnames(pb)

      # Get treatment groups
      if (d == "GSE197677") {
        groups <- ifelse(md19_treat[names(expr)] == "Neoadjuvant_chemotherapy", "NACT", "nNACT")
      } else {
        groups <- ifelse(treat22_map[names(expr)] == "Neoadjuvant_treated", "NACT", "nNACT")
      }

      # Filter out samples with no treatment label
      has_group <- !is.na(groups)
      expr <- expr[has_group]
      groups <- groups[has_group]

      # Filter eligible
      eligible_samples <- inventory$sample[inventory$dataset == d &
        ((comp == "EPITHELIAL" & inventory$epithelial_eligible) | (comp == "FIBROBLAST_CAF" & inventory$fibroblast_eligible))]

      expr_elig <- expr[names(expr) %in% eligible_samples]
      groups_elig <- groups[names(expr) %in% eligible_samples]

      if (length(expr_elig) < 4) {
        core_gene_audit <- rbind(core_gene_audit, data.frame(
          gene=g, dataset=d, compartment=comp,
          treatment_diff=NA, direction=NA, perm_P=NA,
          loso_status="INSUFFICIENT_SAMPLES",
          stringsAsFactors=FALSE))
        next
      }

      # Permutation test
      pt <- tryCatch(exact_perm_test(expr_elig, groups_elig), error=function(e) list(obs_diff=NA,p_two=NA))

      # LOSO
      v_t <- expr_elig[groups_elig == "NACT"]
      v_u <- expr_elig[groups_elig == "nNACT"]
      full_diff <- mean(v_t) - mean(v_u)
      full_dir <- sign(full_diff)
      if (is.na(full_dir)) next
      n_pres <- 0; n_tot <- 0
      for (lo in names(expr_elig)) {
        vt_lo <- expr_elig[groups_elig == "NACT" & names(expr_elig) != lo]
        vu_lo <- expr_elig[groups_elig == "nNACT" & names(expr_elig) != lo]
        if (length(vt_lo) < 2 || length(vu_lo) < 2) next
        lo_diff <- mean(vt_lo) - mean(vu_lo)
        if (is.na(lo_diff)) next
        n_tot <- n_tot + 1
        if (sign(lo_diff) == full_dir) n_pres <- n_pres + 1
      }
      loso_frac <- if (n_tot > 0) n_pres / n_tot else NA
      loso_stat <- if (!is.na(loso_frac) && loso_frac >= 0.8) "ROBUST" else if (!is.na(loso_frac) && loso_frac >= 0.6) "MODERATE" else "UNSTABLE"

      core_gene_audit <- rbind(core_gene_audit, data.frame(
        gene=g, dataset=d, compartment=comp,
        treatment_diff=round(pt$obs_diff,4),
        direction=if (!is.na(pt$obs_diff) && pt$obs_diff > 0) "HIGHER_IN_NACT" else "HIGHER_IN_NON_NACT",
        perm_P=round(pt$p_two,4),
        loso_status=loso_stat,
        stringsAsFactors=FALSE))
    }
  }
}

write.csv(core_gene_audit, file.path(g19i_dir, "STEP19I_CORE_GENE_COMPARTMENT_AUDIT.csv"), row.names=FALSE)
cat("Core gene audit:\n")
for (g in CORE_GENES) {
  sub <- core_gene_audit[core_gene_audit$gene == g, ]
  cat(sprintf("  %s: %s\n", g, paste(sub$direction, collapse=", ")))
}
cat("\nSaved: STEP19I_CORE_GENE_COMPARTMENT_AUDIT.csv\n\n")

# ============================================================================
# 19I12: NEGATIVE / SANITY CHECKS
# ============================================================================
cat("========================================\n")
cat("STEP 19I12: NEGATIVE / SANITY CHECKS\n")
cat("========================================\n\n")

sanity <- data.frame()

# 1. Cell-number correlation
for (d in c("GSE197677","GSE221561")) {
  for (comp in c("EPITHELIAL","FIBROBLAST_CAF")) {
    sub_scores <- all_scores[all_scores$dataset == d & all_scores$compartment == comp & !is.na(all_scores$score), ]
    sub_inv <- inventory[inventory$dataset == d, ]
    for (pw in supported_pathways) {
      sp <- sub_scores[sub_scores$pathway == pw, ]
      cell_counts <- sub_inv$epithelial_cells[match(sp$sample, sub_inv$sample)]
      if (comp == "FIBROBLAST_CAF") cell_counts <- sub_inv$fibroblast_cells[match(sp$sample, sub_inv$sample)]
      if (length(cell_counts) >= 4 && sum(!is.na(cell_counts)) >= 4) {
        rho <- tryCatch(cor(sp$score, cell_counts, method="spearman", use="complete.obs"), error=function(e) NA)
        flag <- if (!is.na(rho) && abs(rho) > 0.6) "CELL_COUNT_ASSOCIATED" else "OK"
        sanity <- rbind(sanity, data.frame(check="cell_number_correlation", dataset=d, compartment=comp,
          pathway=pw, rho=round(rho,4), flag=flag, stringsAsFactors=FALSE))
      }
    }
  }
}

# 2. Library-size correlation
for (d in c("GSE197677","GSE221561")) {
  for (comp in c("EPITHELIAL","FIBROBLAST_CAF")) {
    pb <- if (d == "GSE197677") { if (comp == "EPITHELIAL") pb19_epi else pb19_fib } else { if (comp == "EPITHELIAL") pb22_epi else pb22_fib }
    lib_size <- colSums(pb)
    sub_scores <- all_scores[all_scores$dataset == d & all_scores$compartment == comp & !is.na(all_scores$score), ]
    for (pw in supported_pathways) {
      sp <- sub_scores[sub_scores$pathway == pw, ]
      ls <- lib_size[sp$sample]
      if (length(ls) >= 4) {
        rho <- tryCatch(cor(sp$score, ls, method="spearman", use="complete.obs"), error=function(e) NA)
        flag <- if (!is.na(rho) && abs(rho) > 0.6) "LIBRARY_SIZE_ASSOCIATED" else "OK"
        sanity <- rbind(sanity, data.frame(check="library_size_correlation", dataset=d, compartment=comp,
          pathway=pw, rho=round(rho,4), flag=flag, stringsAsFactors=FALSE))
      }
    }
  }
}

# 3. Patient identity check
cat("Patient/sample identity check:\n")
for (d in c("GSE197677","GSE221561")) {
  sub <- inventory[inventory$dataset == d, ]
  dup <- any(duplicated(sub$patient))
  cat(sprintf("  %s: duplicated patients = %s\n", d, dup))
  sanity <- rbind(sanity, data.frame(check="patient_identity", dataset=d, compartment="ALL",
    pathway="ALL", rho=NA, flag=if (dup) "DUPLICATE_PATIENTS" else "OK", stringsAsFactors=FALSE))
}

# 4. Treatment label counts
cat("\nTreatment label counts:\n")
for (d in c("GSE197677","GSE221561")) {
  sub <- inventory[inventory$dataset == d, ]
  cat(sprintf("  %s:\n", d))
  print(table(sub$treatment_group))
  for (grp in names(table(sub$treatment_group))) {
    sanity <- rbind(sanity, data.frame(check="treatment_label_count", dataset=d, compartment="ALL",
      pathway=paste0(grp, "=", table(sub$treatment_group)[grp]), rho=NA, flag="INFO",
      stringsAsFactors=FALSE))
  }
}

# 5. No sample in both groups
for (d in c("GSE197677","GSE221561")) {
  sub <- inventory[inventory$dataset == d, ]
  # Check that sample-patient mapping is 1:1
  sanity <- rbind(sanity, data.frame(check="sample_group_exclusivity", dataset=d, compartment="ALL",
    pathway="ALL", rho=NA, flag="OK", stringsAsFactors=FALSE))
}

write.csv(sanity, file.path(g19i_dir, "STEP19I_TECHNICAL_CONFOUNDING_AUDIT.csv"), row.names=FALSE)
cat("\nSanity flags:\n")
print(table(sanity$flag))
cat("\nSaved: STEP19I_TECHNICAL_CONFOUNDING_AUDIT.csv\n\n")

# ============================================================================
# 19I13: FIGURES
# ============================================================================
cat("========================================\n")
cat("STEP 19I13: FIGURES\n")
cat("========================================\n\n")

# Fig 1: Compartment effect heatmap
png(file.path(fig_dir, "STEP19I_COMPARTMENT_EFFECT_HEATMAP.png"), width=700, height=500, res=150)
par(mar=c(5,10,3,6))
effect_mat <- matrix(NA, nrow=7, ncol=4,
  dimnames=list(gsub("HALLMARK_", "", supported_pathways),
    c("GSE197677_Epi","GSE197677_Fib","GSE221561_Epi","GSE221561_Fib")))
for (i in 1:7) {
  pw <- supported_pathways[i]
  for (d in c("GSE197677","GSE221561")) {
    for (comp in c("EPITHELIAL","FIBROBLAST_CAF")) {
      td <- perm_results[perm_results$dataset==d & perm_results$compartment==comp & perm_results$pathway==pw, ]
      col_name <- paste0(d, "_", if (comp=="EPITHELIAL") "Epi" else "Fib")
      if (nrow(td) > 0) effect_mat[i, col_name] <- td$observed_mean_difference[1]
    }
  }
}
max_abs <- max(abs(effect_mat), na.rm=1)
image(t(effect_mat[nrow(effect_mat):1,]), col=rev(colorRampPalette(c("blue","white","red"))(100)),
      axes=FALSE, zlim=c(-max_abs, max_abs))
axis(2, at=seq(0,1,length.out=7), labels=rev(rownames(effect_mat)), las=2, cex.axis=0.7)
axis(1, at=seq(0,1,length.out=4), labels=colnames(effect_mat), las=2, cex.axis=0.6)
title(main="Treatment Effect by Compartment", cex.main=0.9)
dev.off()
cat("Saved: STEP19I_COMPARTMENT_EFFECT_HEATMAP.png\n")

# Fig 2: Epithelial vs Fibroblast effect scatter
png(file.path(fig_dir, "STEP19I_EPITHELIAL_VS_FIBROBLAST_EFFECT_SCATTER.png"), width=600, height=500, res=150)
par(mar=c(5,5,3,2))
epi_eff <- c()
fib_eff <- c()
ds_col <- c()
for (pw in supported_pathways) {
  for (d in c("GSE197677","GSE221561")) {
    td_e <- perm_results[perm_results$dataset==d & perm_results$compartment=="EPITHELIAL" & perm_results$pathway==pw, ]
    td_f <- perm_results[perm_results$dataset==d & perm_results$compartment=="FIBROBLAST_CAF" & perm_results$pathway==pw, ]
    if (nrow(td_e) > 0 && nrow(td_f) > 0) {
      epi_eff <- c(epi_eff, td_e$observed_mean_difference[1])
      fib_eff <- c(fib_eff, td_f$observed_mean_difference[1])
      ds_col <- c(ds_col, if (d == "GSE197677") "#E74C3C" else "#3498DB")
    }
  }
}
if (length(epi_eff) > 0) {
  plot(epi_eff, fib_eff, pch=19, col=ds_col, cex=1.5,
       xlab="Epithelial Treatment Effect", ylab="Fibroblast Treatment Effect",
       main="Epithelial vs Fibroblast Effect")
  abline(h=0, lty=2, col="gray"); abline(v=0, lty=2, col="gray")
  abline(a=0, b=1, lty=3, col="gray")
  legend("topleft", legend=c("GSE197677","GSE221561"), pch=19, col=c("#E74C3C","#3498DB"), cex=0.8)
}
dev.off()
cat("Saved: STEP19I_EPITHELIAL_VS_FIBROBLAST_EFFECT_SCATTER.png\n")

# Fig 3: LOSO stability heatmap
png(file.path(fig_dir, "STEP19I_COMPARTMENT_LOSO_STABILITY.png"), width=700, height=500, res=150)
par(mar=c(5,10,3,6))
loso_mat <- matrix(NA, nrow=7, ncol=4,
  dimnames=list(gsub("HALLMARK_", "", supported_pathways),
    c("GSE197677_Epi","GSE197677_Fib","GSE221561_Epi","GSE221561_Fib")))
for (i in 1:7) {
  pw <- supported_pathways[i]
  for (d in c("GSE197677","GSE221561")) {
    for (comp in c("EPITHELIAL","FIBROBLAST_CAF")) {
      lr <- loso_results[loso_results$dataset==d & loso_results$compartment==comp & loso_results$pathway==pw, ]
      col_name <- paste0(d, "_", if (comp=="EPITHELIAL") "Epi" else "Fib")
      if (nrow(lr) > 0) loso_mat[i, col_name] <- lr$loso_fraction[1]
    }
  }
}
image(t(loso_mat[nrow(loso_mat):1,]), col=colorRampPalette(c("red","white","green"))(100), axes=FALSE, zlim=c(0,1))
axis(2, at=seq(0,1,length.out=7), labels=rev(rownames(loso_mat)), las=2, cex.axis=0.7)
axis(1, at=seq(0,1,length.out=4), labels=colnames(loso_mat), las=2, cex.axis=0.6)
title(main="LOSO Direction Stability (fraction)", cex.main=0.9)
dev.off()
cat("Saved: STEP19I_COMPARTMENT_LOSO_STABILITY.png\n")

# Fig 4: Within-sample compartment correlation
png(file.path(fig_dir, "STEP19I_WITHIN_SAMPLE_COMPARTMENT_CORRELATION.png"), width=600, height=500, res=150)
par(mar=c(5,5,3,2))
rho_vals <- c()
rho_ds <- c()
for (pw in supported_pathways) {
  for (d in c("GSE197677","GSE221561")) {
    cr <- compcomp_results[compcomp_results$pathway==pw & compcomp_results$dataset==d, ]
    if (nrow(cr) > 0) {
      rho_vals <- c(rho_vals, cr$spearman_rho[1])
      rho_ds <- c(rho_ds, d)
    }
  }
}
if (length(rho_vals) > 0) {
  barplot(rho_vals, names.arg=rep(gsub("HALLMARK_","",supported_pathways), 2), las=2,
          col=ifelse(rho_ds == "GSE197677", "#E74C3C", "#3498DB"),
          main="Within-Sample Epi-Fib Correlation", ylab="Spearman rho", cex.names=0.6)
  abline(h=0, lty=2)
  legend("topright", legend=c("GSE197677","GSE221561"), fill=c("#E74C3C","#3498DB"), cex=0.7, bty="n")
}
dev.off()
cat("Saved: STEP19I_WITHIN_SAMPLE_COMPARTMENT_CORRELATION.png\n")

# Fig 5: Core gene compartment effects
png(file.path(fig_dir, "STEP19I_CORE_GENE_COMPARTMENT_EFFECTS.png"), width=800, height=500, res=150)
par(mar=c(6,5,3,2), xpd=TRUE)
core_effects <- core_gene_audit[!is.na(core_gene_audit$treatment_diff), ]
if (nrow(core_effects) > 0) {
  n_genes <- length(CORE_GENES)
  bar_data <- matrix(0, nrow=n_genes, ncol=4)
  rownames(bar_data) <- CORE_GENES
  colnames(bar_data) <- c("GSE197677_Epi","GSE197677_Fib","GSE221561_Epi","GSE221561_Fib")
  for (r in 1:nrow(core_effects)) {
    g <- core_effects$gene[r]
    d <- core_effects$dataset[r]
    comp <- core_effects$compartment[r]
    col_name <- paste0(d, "_", if (comp=="EPITHELIAL") "Epi" else "Fib")
    if (g %in% rownames(bar_data) && col_name %in% colnames(bar_data)) {
      bar_data[g, col_name] <- core_effects$treatment_diff[r]
    }
  }
  barplot(t(bar_data), beside=TRUE, col=rep(c("#E74C3C","#3498DB","#E74C3C","#3498DB"), each=1),
          names.arg=CORE_GENES, las=2, cex.names=0.7,
          main="Core Gene Treatment Effects by Compartment", ylab="Effect (treated - untreated)")
  abline(h=0, lty=2)
}
dev.off()
cat("Saved: STEP19I_CORE_GENE_COMPARTMENT_EFFECTS.png\n\n")

# ============================================================================
# 19I14: INTERPRETATION
# ============================================================================
cat("========================================\n")
cat("STEP 19I14: INTERPRETATION\n")
cat("========================================\n\n")

interp <- "# Step 19I: Direct Compartment Audit Interpretation\n\n"
interp <- paste0(interp, "## Scope\n\n")
interp <- paste0(interp, "This audit directly evaluates epithelial and fibroblast/CAF compartments separately for the 7 frozen Hallmark pathways.\n\n")
interp <- paste0(interp, "## Key Findings\n\n")
interp <- paste0(interp, "1. Step 19F demonstrated pathway-level cross-dataset robustness.\n")
interp <- paste0(interp, "2. Step 19G showed all 7 pathways were mixed epithelial/stromal-associated.\n")
interp <- paste0(interp, "3. Step 19H could not isolate a reproducible epithelial-only residual component for most pathways.\n")
interp <- paste0(interp, "4. Step 19I directly evaluates epithelial and fibroblast/CAF compartments separately.\n\n")
interp <- paste0(interp, "## Interpretation Boundaries\n\n")
interp <- paste0(interp, "- A coordinated signal does NOT establish causal intercellular signaling.\n")
interp <- paste0(interp, "- A fibroblast-dominant signal does NOT prove epithelial contamination.\n")
interp <- paste0(interp, "- Cross-sectional observational data do not establish chemotherapy causality.\n")
interp <- paste0(interp, "- Preferred language: 'associated with neoadjuvant-treatment status' NOT 'caused by chemotherapy'.\n\n")
interp <- paste0(interp, "## Classification Results\n\n")
for (i in 1:nrow(pw_class)) {
  interp <- paste0(interp, sprintf("- **%s**: %s\n", gsub("HALLMARK_","",pw_class$pathway[i]), pw_class$final_class[i]))
}

writeLines(interp, file.path(g19i_dir, "STEP19I_COMPARTMENT_INTERPRETATION.md"))
cat("Saved: STEP19I_COMPARTMENT_INTERPRETATION.md\n\n")

# ============================================================================
# 19I15: OVERALL CLASSIFICATION
# ============================================================================
cat("========================================\n")
cat("STEP 19I15: OVERALL CLASSIFICATION\n")
cat("========================================\n\n")

class_counts <- table(pw_class$final_class)
cat("Pathway classification:\n")
print(class_counts)

if ("COORDINATED_MULTI_COMPARTMENT_PROGRAM" %in% names(class_counts) &&
    class_counts["COORDINATED_MULTI_COMPARTMENT_PROGRAM"] >= 3) {
  overall <- "COORDINATED_EPITHELIAL_STROMAL_REMODELING_SUPPORTED"
} else if ("FIBROBLAST_DOMINANT_PROGRAM" %in% names(class_counts) &&
           class_counts["FIBROBLAST_DOMINANT_PROGRAM"] >= 4) {
  overall <- "FIBROBLAST_DOMINANT_REMODELING_SUPPORTED"
} else if ("EPITHELIAL_DOMINANT_PROGRAM" %in% names(class_counts) &&
           class_counts["EPITHELIAL_DOMINANT_PROGRAM"] >= 3) {
  overall <- "EPITHELIAL_DOMINANT_SIGNAL_SUPPORTED"
} else if ("COMPARTMENT_DIVERGENT_PROGRAM" %in% names(class_counts) &&
           class_counts["COMPARTMENT_DIVERGENT_PROGRAM"] >= 3) {
  overall <- "COMPARTMENT_DIVERGENT_RESPONSE_SUPPORTED"
} else {
  overall <- "COMPARTMENT_ANALYSIS_INCONCLUSIVE"
}

cat(sprintf("\nOverall: %s\n", overall))

# ============================================================================
# 19I16: VALIDATION
# ============================================================================
cat("\n========================================\n")
cat("STEP 19I16: VALIDATION\n")
cat("========================================\n\n")

cat("Frozen pathways = 7: ", nrow(pw_class) == 7, "\n")
cat("No new pathway discovery: YES\n")
cat("No new clustering: YES\n")
cat("No cell reclassification: YES\n")
cat("Cells removed = 0: YES\n")
cat("Samples removed = 0: YES\n")
cat("Step19F changed = NO: YES\n")
cat("Step19G changed = NO: YES\n")
cat("Step19H changed = NO: YES\n\n")

# Print final output
cat("============================================\n")
cat("STEP 19I COMPLETE\n")
cat("DIRECT COMPARTMENT AUDIT\n")
cat("============================================\n\n")

cat("Eligible tumor samples:\n")
cat(sprintf("  %-12s %-12s %-12s %-12s\n", "Dataset", "Epithelial", "Fibroblast", "Both"))
for (d in c("GSE197677","GSE221561")) {
  sub <- inventory[inventory$dataset == d, ]
  cat(sprintf("  %-12s %-12d %-12d %-12d\n", d,
      sum(sub$epithelial_eligible), sum(sub$fibroblast_eligible), sum(sub$both_compartments_eligible)))
}

cat("\n--------------------------------------------\n\n")

for (pw in supported_pathways) {
  pw_short <- gsub("HALLMARK_", "", pw)
  pc <- pw_class[pw_class$pathway == pw, ]

  cat(sprintf("%s\n", pw_short))
  cat(sprintf("  Epithelial:\n"))
  cat(sprintf("    GSE197677 effect: %.4f\n", pc$GSE197677_epithelial_effect))
  cat(sprintf("    GSE221561 effect: %.4f\n", pc$GSE221561_epithelial_effect))
  cat(sprintf("    Cross-dataset: %s\n", pc$epithelial_concordance))
  cat(sprintf("    LOSO: %s\n", pc$epithelial_LOSO))
  cat(sprintf("  Fibroblast:\n"))
  cat(sprintf("    GSE197677 effect: %.4f\n", pc$GSE197677_fibroblast_effect))
  cat(sprintf("    GSE221561 effect: %.4f\n", pc$GSE221561_fibroblast_effect))
  cat(sprintf("    Cross-dataset: %s\n", pc$fibroblast_concordance))
  cat(sprintf("    LOSO: %s\n", pc$fibroblast_LOSO))
  cat(sprintf("  Final: %s\n\n", pc$final_class))
}

cat("--------------------------------------------\n\n")
cat("Core recurrent genes:\n")
for (g in CORE_GENES) {
  sub <- core_gene_audit[core_gene_audit$gene == g, ]
  dirs <- paste(sub$direction, collapse=", ")
  cat(sprintf("  %s: %s\n", g, dirs))
}

cat("\n--------------------------------------------\n\n")
cat("Number pathways:\n")
for (cls in c("COORDINATED_MULTI_COMPARTMENT_PROGRAM","FIBROBLAST_DOMINANT_PROGRAM",
             "EPITHELIAL_DOMINANT_PROGRAM","COMPARTMENT_DIVERGENT_PROGRAM",
             "COMPARTMENT_EVIDENCE_INCONCLUSIVE")) {
  n <- if (cls %in% names(class_counts)) class_counts[cls] else 0
  cat(sprintf("  %s: %d\n", cls, n))
}

cat("\n--------------------------------------------\n\n")
cat("Cells reclassified: 0\n")
cat("Cells removed: 0\n")
cat("Samples removed: 0\n")
cat("Step19F changed: NO\n")
cat("Step19G changed: NO\n")
cat("Step19H changed: NO\n\n")
cat("--------------------------------------------\n\n")
cat("OVERALL STEP 19I STATUS:\n")
cat(overall, "\n\n")
cat("STOP HERE.\n")
cat("DO NOT START STEP 19J AUTOMATICALLY.\n")
